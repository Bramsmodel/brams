!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module ModGasPart

  use ModRcio, only: &
       cio

  use ModGaspartFields, only: &
       GaspartFields
  
  use node_mod, only: &
       mchnum,        & ! INTENT(IN)
       master_num, &    ! INTENT(IN)
       nodeibcon,     & !INTENT(IN)
       mynum            !INTENT(IN)

  use grid_dims, only: &
       maxgrds ! INTENT(IN)

  use mem_grid, only: &    !INTENT(IN)
       ngrids,        &    !INTENT(IN)
       grid_g,        &    !INTENT(IN)
       dzt,  &             !INTENT(IN)
       ngridsh          ! INTENT(OUT)
  
  use mem_emiss, only: &
       EINDNO,         & !INTENT(IN)
       EINDNO2,        & !INTENT(IN)
       EINDPM,         & !INTENT(IN)
       EINDCO,         & !INTENT(IN)
       EINDSO2,        & !INTENT(IN)
       EINDVOC,        & !INTENT(IN)
       EVEINO,         & !INTENT(IN)
       EVEINO2,        & !INTENT(IN)
       EVEIPM,         & !INTENT(IN)
       EVEICO,         & !INTENT(IN)
       EVEISO2,        & !INTENT(IN)
       EVEIVOC,        & !INTENT(IN)
       EFSAT,          & !INTENT(IN)
       EFSUN,          & !INTENT(IN)
       WEEKDAYIN, &         !INTENT(IN)
       chemdata_in      ! INTENT(IN)

  use ParLib, only: &
       parf_bcast       ! Subroutine

  use an_header, only: &
       head_table,     &  ! Type
       nvbtab             ! INTENT(OUT)
  
  use ModVarTable, only: &
       VarTable

  use teb_vars_const, only : &
       RUSHH1,               & !INTENT(IN)
       RUSHH2,               & !INTENT(IN)
       DAYLIGHT                !INTENT(IN)

  use mem_leaf, only: &
       leaf_g              !INTENT(IN)
  
  use ModBasicFields, only: &
       BasicFields            !INTENT(IN)

  implicit none

  private
  public :: init_conc1
  public :: init_conc2
  public :: init_conc_prev
  public :: le_fontes
  public :: sources_teb
  public :: emfactor
  
contains



  subroutine le_fontes(ng, n1, n2, n3, np, ia, iz, ja, jz, time, &
       oneBasicFields, oneGaspartFields)
    ! Arguments:
    integer, intent(in) :: ng, n1, n2, n3, np, ia, iz, ja, jz
    real, intent(in)    :: time
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(GaspartFields), pointer, intent(in) :: oneGaspartFields
    
    ! Local Variables:
    integer :: ig, kl, i, j
    integer, parameter :: ngases=6
    integer :: kgas(ngases)
    integer :: len1
    character(len=10) :: gas(ngases), tracer

    ! Init arrays: gas, kgas
    gas(:)  = (/'NO  ', 'NO2 ', 'PM25', 'CO  ', 'SO2 ', 'VOCS'/)
    kgas(:) = (/  1,    4,     7,     10,   13,     16/)

    !------------ Allocation table of qsc = dum1------------------------
    ! Distribute sources in order to avoid highly localized sources:
    ! (somente definidas em 1 ponto de grade)
    ! kgas(1) = 1,2,3    => NO
    ! kgas(2) = 4,5,6    => NO2
    ! kgas(3) = 7,8,9    => PM25
    ! kgas(4) = 10,11,12 => CO
    ! kgas(5) = 13,14,15 => SO2 
    ! kgas(6) = 16,17,18 => VOCS
    !----------------------------------------------------------------------


    !sources only in the inner grid
    if (ng/=ngrids) return

    !Allocate memory in DUM1 vector for sources:
    ! first kl 3D vector oneGaspartFields%gasr(1,1,1) levels are left for
    ! emission sources calculation
    kl = 3*ngases 

    oneGaspartFields%gasr(:,:,:) = 0.

    do ig=1,ngases

       !Call dum1_zero(n1,n2,n3,ia,iz,ja,jz,kl,scalar_g(ig,ng)%dum1(1,1,1))  

       tracer = gas(ig)

       len1   = LEN_trim(tracer) + 1

       call read_sources_teb(ng, n1, n2, n3, np, ia, iz, ja, jz, &
            oneGaspartFields%gasr, oneGaspartFields%fusog, &
            tracer(1:len1), kgas(ig), leaf_g(ng)%G_URBAN, &
            grid_g(ng)%dxt, grid_g(ng)%dyt, time       )

       call reorganize_sources_teb(ng, n1, n2, n3, ia, iz, ja, jz, &
            oneGaspartFields%gasr, kgas(ig)                )

       call convert_to_misture_ratio_teb(ng, n1, n2, n3, ia, iz, ja, jz,     &
            kgas(ig), oneGaspartFields%gasr,                             &
            oneBasicFields%dn0, grid_g(ng)%rtgt,                    &
            grid_g(ng)%dxt, grid_g(ng)%dyt, dzt                    )

    enddo !end of gases' looping
  end subroutine le_fontes

  !--------------------------------------------------

  subroutine read_sources_teb(ng, n1, n2, n3, np, ia, iz, ja, jz, &
       qsc, fuso, gas, kgas, schar, dxt, dyt, time)
    ! Arguments:
    integer, intent(IN)          :: ng, n1, n2, n3, np, ia, iz, ja, jz
    real, intent(INOUT)          :: qsc(n1,n2,n3)
    real, intent(IN)             :: fuso(n2,n3)
    character(len=*), intent(IN) :: gas
    integer, intent(IN)          :: kgas
    real, intent(IN)             :: schar(n2,n3,np)
    real, intent(IN)             :: dxt(n2,n3)
    real, intent(IN)             :: dyt(n2,n3)
    real, intent(IN)             :: time
    ! Local Variables:
    integer :: i, j, idays
    character(len=3) :: cday
    real :: pfat, pfat2, emiss, emiind, area, pft, r_q, tign, ax1, ax2, bx1, bx2
    real :: timeq1, timeq2

    !****************************************************************************
    !defining the emission rate for each gas/particle
    !****************************************************************************

    !The following values are provided by CETESB's report for 2001
    !these emission rates are relative to an area of 8 million m2, wich 
    !represent the total area of RMSP. However, the urbanized part of this
    !area is about 1.5 and it will be considered here.

    !    veicular emissions rate in ug/min/m2 (for an area of 8 million of m2)

    !             CO   = 390.6
    !             NOx  =  89.3
    !             NO   = NOx*0.9
    !             NO2  = NOx*0.1
    !             VOC's=  91.0
    !             SO2  =   5.3

    ! In order to consider the diurnal cycle of veicular emission, the emission
    ! rates will have units of kg/m2/day

    !    veicular emissions rate in kg/day/m2 (for an area of 1.5 thousand of m2)

    !             CO   = 3.0173515E-03
    !             NOx  = 6.8566209E-04
    !             NO   = NOx*0.9 =  6.1709585E-04
    !             NO2  = NOx*0.1 =  6.8566209E-05
    !             VOC's= 6.9954334E-04
    !             SO2  = 4.0730592E-05
    !             PM15 = 6.2648396E-06
    !    industrial emissions rate in kg/s/m2 (for an area of 1.5 thousand of m2)
    !
    !  emiind
    !             CO   = 8.1599860E-10
    !             NOx  = 6.8566209E-04
    !             NO   = 2.6636227E-10
    !             NO2  = 2.9595805E-11
    !             VOC's= 2.5367833E-10
    !             SO2  = 3.6149164E-10
    !             PM25 = 4.3421278E-10

    if     (gas=='CO')   then
       emiss  = EVEICO 
       emiind = EINDCO
    elseif (gas=='NO')   then
       emiss  = EVEINO
       emiind = EINDNO
    elseif (gas=='NO2')  then
       emiss  = EVEINO2
       emiind = EINDNO2
    elseif (gas=='VOCS') then
       emiss  = EVEIVOC
       emiind = EINDVOC
    elseif (gas=='SO2')  then
       emiss  = EVEISO2
       emiind = EINDSO2
    elseif (gas=='PM25') then
       emiss  = EVEIPM
       emiind = EINDPM
    endif

    pft   = emiss/273234.9
    !the value of 273234.9 was obtained in order to have the integral for
    !one day = emiss
    !it is used to distribute emissions following the diurnal cycle 

    tign  = 0.0*3600.                            !UTC time of ignition 

    !print*,'valor de time=',time

    idays = int((time/3600.)/24.0 + 0.00001)  !number of days of simulation
    tign  = tign + real(idays)*24.0*3600.0

    pfat  = 1.0

    call EMFACTOR(WEEKDAYIN,idays,cday)
    if     (cday=='SAT') then
       pfat=EFSAT
    elseif (cday=='SUN') then
       pfat=EFSUN
    endif

    !******************************************************************
    !gaussian distribution                                            *
    !                                     timeqi                      *
    !                                      ____                       *
    !                               _     |    |      _               *
    !              1               |      (x-mi)**2    |              *
    !f(x)= ----------------- * EXP | -  -------------- |  i=1,2       *
    !      sigmai* sqrt(2*PI)      |    2*(sigmai**2)  |              *
    !     |__________________|      -  |____________| -               *
    !             axi                    8.5 or 10.6                  *
    !******************************************************************

    ax1 = 4.5
    ax2 = 9.2

    do i = 1,n2
       do j= 1,n3

          if (nint(SCHAR(I,J,2))/=0) then

             !bx1=10.81
             !bx2=19.0
             bx1 = RUSHH1 - fuso(i,j) + DAYLIGHT
             bx2 = RUSHH2 - fuso(i,j) + DAYLIGHT

             timeq1 = (time/3600.0 - tign/3600.0) - bx1
             timeq2 = (time/3600.0 - tign/3600.0) - bx2

             r_q = pft*(ax1*exp(-(timeq1)**2/8.5) + &
                  ax2*exp(-(timeq2)**2/10.6))*pfat 

             area = 1./(dxt(i,j)*dyt(i,j))

             !    Identificando pontos de grade com a classificacao urbano 1 
             if     (nint(SCHAR(I,J,2))==1) then
                pfat2 = 1.
                ! Identificando pontos de grade com a classificacao urbano 2 
             elseif (nint(SCHAR(I,J,2))==2) then
                pfat2 = 0.33333333
                ! Identificando pontos de grade com a classificacao urbano 3 
             elseif (nint(SCHAR(I,J,2))==3) then
                pfat2 = 0.2
             endif

             qsc(kgas,i,j) = (r_q*area*pfat*pfat2) + (emiind*area)

          endif

       enddo
    enddo
  end subroutine read_sources_teb

  !-------------------------------------------------------

  subroutine reorganize_sources_teb(ngrid, n1, n2, n3, ia, iz, ja, jz, qsc, kgas)
    ! Arguments:
    integer, intent(IN) :: ngrid, n1, n2, n3, ia, iz, ja, jz
    real, intent(INOUT) :: qsc(n1,n2,n3)
    integer, intent(IN) :: kgas
    ! Local Variables:
    integer :: i,j, ii, jj
    real    :: f
    integer         :: jInit, jEnd, iInit, iEnd

    !fator de distribuicao 20% para cada um do 9 primeiros vizinhos
    !( incluindo o proprio site i,j)

    f = 0.2

    ! Defining bounds
    if (btest(nodeibcon(mynum,ngrid),0)) then !OESTE
       iInit = 3
    else
       iInit = 1
    endif
    if (btest(nodeibcon(mynum,ngrid),1)) then !LESTE
       iEnd = n2-2
    else
       iEnd = n2
    endif
    if (btest(nodeibcon(mynum,ngrid),2)) then !NORTE
       jInit = 3
    else
       jInit = 1
    endif
    if (btest(nodeibcon(mynum,ngrid),3)) then !SUL
       jEnd = n3-2
    else
       jEnd = n3
    endif

    do j=jInit,jEnd
       do i=iInit,iEnd

          ! ponto grade i,j,k(=2)
          qsc(kgas+1,i,j) = qsc(kgas+1,i,j) + (1.-f)*qsc(kgas,i,j)

          !distribuicao nos 9 sites em torno de i,j    !  j+1  .   .   .
          do jj = max(j-1,1),min(j+1,n3)
             do ii = max(i-1,1),min(i+1,n2)

                qsc(kgas+1,ii,jj) = qsc(kgas+1,ii,jj) + &
                     (1./9.)*f*qsc(kgas,i,j)  ! ponto grade i,j,k= 2 e 3

             enddo
          enddo

       enddo
    enddo
  end subroutine reorganize_sources_teb

  !-------------------------------------------------------------

  subroutine convert_to_misture_ratio_teb(ng, n1, n2, n3, ia, iz, ja, jz, &
       kgas, qsc, dn0, rtgt, dxt, dyt, dzt) !zt
    ! Arguments:
    integer, intent(IN) :: ng, n1, n2, n3, ia,iz,ja,jz, kgas
    real, intent(INOUT) :: qsc(n1,n2,n3)
    real, intent(IN)    :: dn0(n1,n2,n3)
    real, intent(IN)    :: rtgt(n2,n3)
    real, intent(IN)    :: dxt(n2,n3)
    real, intent(IN)    :: dyt(n2,n3)
    real, intent(IN)    :: dzt(n1)
    ! Local variables:
    real    :: fcu, vol
    integer :: i, j, k
    !

    !Fator de conversao de unidades
    fcu = 1.        !=> kg [gas/part] /kg [ar]
    !!fcu =1.e+12   !=> ng [gas/part] /kg [ar]
    !fcu =1.e+6      !=> mg [gas/part] /kg [ar]  

    do j = 1,n3
       do i = 1,n2
          do k = 2,2

             vol             = 1./(dxt(i,j)*dyt(i,j)*dzt(k)*rtgt(i,j))

             qsc(kgas+1,i,j) = fcu*qsc(kgas+1,i,j)/(vol*dn0(k,i,j))

          enddo
       enddo
    enddo

  end subroutine convert_to_misture_ratio_teb

  !------------------------------------------------------------

  subroutine sources_teb(n1, n2, n3, ia, iz, ja, jz, ig, ngrids,&
       oneGaspartFields)
    ! Arguments:
    integer, intent(in) :: n1, n2, n3, ia, iz, ja, jz, ig, ngrids
    type(GaspartFields), pointer, intent(in) :: oneGaspartFields

    if (ig==ngrids) then
       call EMISSAO(n1, n2, n3, ia, iz, ja, jz, oneGaspartFields)
    endif

  end subroutine sources_teb

  !------------------------------------------------------------

  subroutine EMISSAO(n1, n2, n3, ia, iz, ja, jz, &
       oneGaspartFields)
    ! Arguments:
    integer, intent(in)               :: n1, n2, n3, ia, iz, ja, jz
    type(GaspartFields), pointer, intent(in) :: oneGaspartFields

    call tendgas(n1,n2,n3,ia,iz,ja,jz,oneGaspartFields%pnot  ,oneGaspartFields%gasr,2)
    call tendgas(n1,n2,n3,ia,iz,ja,jz,oneGaspartFields%pno2t ,oneGaspartFields%gasr,5)
    call tendgas(n1,n2,n3,ia,iz,ja,jz,oneGaspartFields%ppm25t,oneGaspartFields%gasr,8)
    call tendgas(n1,n2,n3,ia,iz,ja,jz,oneGaspartFields%pcot  ,oneGaspartFields%gasr,11)
    call tendgas(n1,n2,n3,ia,iz,ja,jz,oneGaspartFields%pso2t ,oneGaspartFields%gasr,14)
    call tendgas(n1,n2,n3,ia,iz,ja,jz,oneGaspartFields%pvoct ,oneGaspartFields%gasr,17)

  end subroutine EMISSAO

  !---------------------------------------------------------------

  subroutine TENDGAS(n1, n2, n3, ia, iz, ja, jz, tende, qsc, k)
    ! Arguments:
    integer, intent(in) :: n1, n2, n3, ia, iz, ja, jz, k
    real, intent(inout) :: tende(n1,n2,n3)
    real, intent(in)    :: qsc(n1,n2,n3)
    ! Local Variables:
    integer :: i, j

    do i=ia,iz
       do j=ja,jz
          tende(2,i,j) = tende(2,i,j) + qsc(k,i,j)
       enddo
    enddo

  end subroutine TENDGAS

  !----------------------------------------------------------------

  subroutine init_conc1(ictrl, ngrid, n1, n2, n3, np, &
       G_URBAN, no,                                   &
       no2, pm25,                                     &
       co, voc,                                       &
       so2, so4, aer, zt                              )
    ! Arguments:
    integer, intent(in) :: ictrl, ngrid, n1, n2, n3, np
    real, intent(inout) :: G_URBAN(n2,n3,np)
    real, intent(inout) :: no(n1,n2,n3)
    real, intent(inout) :: no2(n1,n2,n3)
    real, intent(inout) :: pm25(n1,n2,n3)
    real, intent(inout) :: co(n1,n2,n3)
    real, intent(inout) :: voc(n1,n2,n3)
    real, intent(inout) :: so2(n1,n2,n3)
    real, intent(inout) :: so4(n1,n2,n3)
    real, intent(inout) :: aer(n1,n2,n3)
    real, intent(inout) :: zt(n1)
    ! Local Variables:
    real, parameter :: PMAR=28.97
    integer         :: i, j, ii, jj, k
    real            :: pfat, f, expo
    real            :: nox(n1,n2,n3)
    real            :: no2x(n1,n2,n3)
    real            :: pm25x(n1,n2,n3)
    real            :: cox(n1,n2,n3)
    real            :: vocx(n1,n2,n3)
    real            :: so2x(n1,n2,n3)
    real            :: no0, no20, pm250, co0, voc0, so20
    integer         :: jInit, jEnd, iInit, iEnd

    no(:,:,:)   = 0.       !NO
    no2(:,:,:)  = 0.       !NO2
    pm25(:,:,:) = 0.       !PM25
    co(:,:,:)   = 0.       !CO
    so2(:,:,:)  = 0.       !SO2
    so4(:,:,:)  = 0.       !SO4
    voc(:,:,:)  = 0.       !VOCS
    aer(:,:,:)  = 0.       !VOCS

    if (ictrl==1) then
       !--------
       !Flag which defines initial condition (IF RUNTYPE = INITIAL)
       do j=1,n3
          do i=1,n2

             pfat = 0.0

             !identify surface type by using G_URBAN 
             if (nint(G_URBAN(i,j,2))/=0.) then

                if (nint(G_URBAN(i,j,2))==1) then
                   pfat = 1.
                endif

                if (nint(G_URBAN(i,j,2))==2) then
                   pfat = 0.3     !using only 30% of values for urban regions
                endif
                if (nint(G_URBAN(i,j,2))==3) then
                   pfat = 0.2     !using only 20% of values for urban regions
                endif
             else
                pfat = 0.1
             endif

             !defining surface backgroud concentrations (typical values for 00 Z in urban regions)
             !can be better difined if using real values (eg. read a concentration file)

             no0   = 248.04 * pfat  !NO  (estacao congonhas dia 31/07/99 ug/m3)
             no20  =  62.30 * pfat  !NO2 (estacao congonhas dia 31/07/99 ug/m3)
             pm250 =  17.28 * pfat  !PM25(estacao santana   dia 31/07/99 ug/m3)
             co0   =   1.20 * pfat  !CO  (estacao congonhas dia 31/07/99 ppm)
             so20  =  17.50 * pfat  !SO2 (estacao congonhas dia 31/07/99 ug/m3)
             voc0  =   0.31 * pfat  !VOCS (Leila)

             do k=2,n1-1
                expo         = exp(-zt(k)/zt(n1-1))

                nox  (k,i,j) = (no0*expo)  *(1.e-9)/1.275     !NO
                no2x (k,i,j) = (no20*expo) *(1.e-9)/1.275     !NO2
                pm25x(k,i,j) = (pm250*expo)*(1.e-9)/1.275     !PM25
                cox  (k,i,j) = (co0*expo)  *(28.00/PMAR)*1.e-6 !CO
                so2x (k,i,j) = (so20*expo) *(1.e-9)/1.275     !SO2
                vocx (k,i,j) = (voc0*expo) *(42.08/PMAR)*1.e-6 !VOCS
             enddo  !end of looping in k(vertical levels)

          enddo  !end of looping in i (longitude)
       enddo   !end of looping in j (latitute)

    endif  !end of if in ictrl

    f = 0.2

    !re-distributing sources
    ! Defining bounds
    if (btest(nodeibcon(mynum,ngrid),0)) then !OESTE
       iInit = 3
    else
       iInit = 1
    endif
    if (btest(nodeibcon(mynum,ngrid),1)) then !LESTE
       iEnd = n2-2
    else
       iEnd = n2
    endif
    if (btest(nodeibcon(mynum,ngrid),2)) then !NORTE
       jInit = 3
    else
       jInit = 1
    endif
    if (btest(nodeibcon(mynum,ngrid),3)) then !SUL
       jEnd = n3-2
    else
       jEnd = n3
    endif

    do j=jInit,jEnd
       do i=iInit,iEnd

          do k=2,n1-1

             no  (k,i,j) = no  (k,i,j) + (1. - f)*nox  (k,i,j)
             no2 (k,i,j) = no2 (k,i,j) + (1. - f)*no2x (k,i,j)
             pm25(k,i,j) = pm25(k,i,j) + (1. - f)*pm25x(k,i,j)
             co  (k,i,j) = co  (k,i,j) + (1. - f)*cox  (k,i,j)
             so2 (k,i,j) = so2 (k,i,j) + (1. - f)*so2x (k,i,j)
             voc (k,i,j) = voc (k,i,j) + (1. - f)*vocx (k,i,j)

             !distribution into the 9 sites around i,j    !  j+1  .   .   .
             !do jj = j-1,j+1                          !   j   .   .   .       
             !do ii = i-1,i+1                        !  j-1  .   .   .
             !      i-1  i  i+1
             do jj = max(j-1,1),min(j+1,n3)
                do ii = max(i-1,1),min(i+1,n2)

                   no  (k,ii,jj) = no  (k,ii,jj) + (1./9.)*f*nox  (k,i,j)
                   no2 (k,ii,jj) = no2 (k,ii,jj) + (1./9.)*f*no2x (k,i,j)
                   pm25(k,ii,jj) = pm25(k,ii,jj) + (1./9.)*f*pm25x(k,i,j)
                   co  (k,ii,jj) = co  (k,ii,jj) + (1./9.)*f*cox  (k,i,j)
                   so2 (k,ii,jj) = so2 (k,ii,jj) + (1./9.)*f*so2x (k,i,j)
                   voc (k,ii,jj) = voc (k,ii,jj) + (1./9.)*f*vocx (k,i,j)
                   so4 (k,ii,jj) = so2 (k,ii,jj)*0.2  
                   aer (k,ii,jj) = so2 (k,ii,jj)*0.0 

                enddo
             enddo

          enddo !vertical
       enddo  !longitude
    enddo   !latitude
  end subroutine init_conc1

  !----------------------------------------------------------------

  subroutine init_conc2(ictrl, ngrid, n1, n2, n3, np,    &
       G_URBAN, s7p, s8p, s9p, s10p, s11p, s12p, s13p, zt)
    ! Arguments:
    integer, intent(in) :: ictrl, ngrid
    integer, intent(in) :: n1, n2, n3, np
    real, intent(in)    :: G_URBAN(n2,n3,np)
    real, intent(out)   :: s7p(n1,n2,n3)
    real, intent(out)   :: s8p(n1,n2,n3)
    real, intent(out)   :: s9p(n1,n2,n3)
    real, intent(out)   :: s10p(n1,n2,n3)
    real, intent(out)   :: s11p(n1,n2,n3)
    real, intent(out)   :: s12p(n1,n2,n3)
    real, intent(out)   :: s13p(n1,n2,n3)
    real, intent(in)    :: zt(n1)
    ! Local variables:
    integer         :: i, j, ii, jj, k
    real, parameter :: PMAR=28.97
    real            :: pfat, f, expo
    real            :: s7p0, s8p0, s9p0, s10p0, s11p0, s12p0, s13p0
    real            :: s7p2(n1,n2,n3)
    real            :: s8p2(n1,n2,n3)
    real            :: s9p2(n1,n2,n3)
    real            :: s10p2(n1,n2,n3)
    real            :: s11p2(n1,n2,n3)
    real            :: s12p2(n1,n2,n3)
    real            :: s13p2(n1,n2,n3)
    integer         :: jInit, jEnd, iInit, iEnd


    if (ictrl==1) then
       !--------
       !Flag which defines initial condition (IF RUNTYPE = INITIAL)

       s7p(:,:,:)  = 0.       !O3
       s8p(:,:,:)  = 0.       !RHCO
       s9p(:,:,:)  = 0.                          !HO2
       s10p(:,:,:) = 0.                          !O3P
       s11p(:,:,:) = 0.                          !O1D
       s12p(:,:,:) = 0.                          !HO
       s13p(:,:,:) = 0.                          !RO2

       do j=1,n3
          do i=1,n2

             !identify surface type by using G_URBAN 
             if (nint(G_URBAN(i,j,2))/=0.) then

                if (nint(G_URBAN(i,j,2))==1) then
                   pfat=1.
                endif

                if (nint(G_URBAN(i,j,2))==2) then
                   pfat=0.3     !using only 30% of values for urban regions
                endif

                if (nint(G_URBAN(i,j,2))==3) then
                   pfat=0.2     !using only 20% of values for urban regions
                endif
             else
                pfat=0.1     
             endif

             ! defining surface backgroud concentrations
             ! (typical values for 00 Z in urban regions)
             ! can be better difined if using real values
             ! (eg. read a concentration file)

             s7p0  = 8.00  *pfat        !O3 (estacao santana dia 31/07/99 ug/m3)
             s8p0  = 0.0187*pfat        !RHCO (valor aproximado segundo Leila,2003, comunicacao pessoal)
             s9p0  = (4. *1.0e-7) *pfat !HO2 ppm
             s10p0 = (4.8*1.0e-10)*pfat !o3p ppm
             s11p0 = (4.8*1.0e-11)*pfat !o1d ppm
             s12p0 = (8. *1.0e-7) *pfat !HO  ppm
             s13p0 = (4. *1.0e-7) *pfat !RO2 ppm

             do k=2,n1-1
                expo         = exp(-zt(k)/zt(n1-1))

                s7p2(k,i,j)  = (s7p0 *expo)*(1.0e-9)/1.275      !O3
                s8p2(k,i,j)  = (s8p0 *expo)*(44.05/PMAR)*1.e-6  !RHCO
                s9p2(k,i,j)  = (s9p0 *expo)*(33.0/PMAR) *1.e-6  !HO2
                s10p2(k,i,j) = (s10p0*expo)*(16.0/PMAR) *1.e-6  !O3P
                s11p2(k,i,j) = (s11p0*expo)*(16.0/PMAR) *1.e-6  !O1D
                s12p2(k,i,j) = (s12p0*expo)*(17.0/PMAR) *1.e-6  !HO
                s13p2(k,i,j) = (s13p0*expo)*(47.0/PMAR) *1.e-6  !RO2

             enddo

          enddo
       enddo

    endif

    f=0.2

    ! Defining bounds
    if (btest(nodeibcon(mynum,ngrid),0)) then !OESTE
       iInit = 3
    else
       iInit = 1
    endif
    if (btest(nodeibcon(mynum,ngrid),1)) then !LESTE
       iEnd = n2-2
    else
       iEnd = n2
    endif
    if (btest(nodeibcon(mynum,ngrid),2)) then !NORTE
       jInit = 3
    else
       jInit = 1
    endif
    if (btest(nodeibcon(mynum,ngrid),3)) then !SUL
       jEnd = n3-2
    else
       jEnd = n3
    endif

    do j=jInit,jEnd
       do i=iInit,iEnd

          do k=2,n1-1

             s7p(k,i,j)  =  s7p(k,i,j) + (1.-f)*s7p2(k,i,j)
             s8p(k,i,j)  =  s8p(k,i,j) + (1.-f)*s8p2(k,i,j)
             s9p(k,i,j)  =  s9p(k,i,j) + (1.-f)*s9p2(k,i,j)
             s10p(k,i,j) = s10p(k,i,j) + (1.-f)*s10p2(k,i,j)
             s11p(k,i,j) = s11p(k,i,j) + (1.-f)*s11p2(k,i,j)
             s12p(k,i,j) = s12p(k,i,j) + (1.-f)*s12p2(k,i,j)
             s13p(k,i,j) = s13p(k,i,j) + (1.-f)*s13p2(k,i,j)

             !distribution into the 9 sites around i,j    !  j+1  .   .   .
             !do jj = j-1,j+1                          !   j   .   .   .       
             !do ii = i-1,i+1                        !  j-1  .   .   .
             !      i-1  i  i+1
             do jj = max(j-1,1),min(j+1,n3)
                do ii = max(i-1,1),min(i+1,n2)

                   s7p(k,ii,jj)  =  s7p(k,ii,jj)+(1./9.)*f* s7p2(k,i,j)
                   s8p(k,ii,jj)  =  s8p(k,ii,jj)+(1./9.)*f* s8p2(k,i,j)
                   s9p(k,ii,jj)  =  s9p(k,ii,jj)+(1./9.)*f* s9p2(k,i,j)
                   s10p(k,ii,jj) = s10p(k,ii,jj)+(1./9.)*f*s10p2(k,i,j)
                   s11p(k,ii,jj) = s11p(k,ii,jj)+(1./9.)*f*s11p2(k,i,j)
                   s12p(k,ii,jj) = s12p(k,ii,jj)+(1./9.)*f*s12p2(k,i,j)
                   s13p(k,ii,jj) = s13p(k,ii,jj)+(1./9.)*f*s13p2(k,i,j)

                enddo
             enddo

          enddo

       enddo
    enddo
  end subroutine init_conc2

  !----------------------------------------------------------------

  subroutine EMFACTOR(dayin, idays, dayout)
    ! Arguments:
    character(len=3), intent(IN)  :: dayin
    integer, intent(IN)           :: idays
    character(len=3), intent(OUT) :: dayout
    ! Local Variables:
    integer ::i, j, k, id, ndays, irest
    character(len=3) :: cday(7)

    ! Init array: cday
    cday(:) = (/'SUN', 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT'/)

    ndays = (idays/7)*7    ! (21)
    irest = idays - ndays  !  (3)

    do i=1,7
       if (dayin==cday(i)) j = i
    enddo

    k = j + irest

    if (k>7) k = k-7

    dayout = cday(k)

  end subroutine EMFACTOR

  !###########################################################################

  subroutine init_conc_prev(oneVarTable, oneVarTableSize)

    ! This routine initializes gas variables from a previous (day-1) history file
    include "constants.h"
    include "files.h"
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(in) :: oneVarTableSize

    ! Local variables:
    integer :: ngrids1, ioutput1,  &
         nnxp1(maxgrds), nnyp1(maxgrds), nnzp1(maxgrds), nzg1, nzs1, npatch1
    integer :: ie, maxarr, ngr, nc
    character (len=f_name_length) :: hnameinh
    integer, parameter :: iunhd=11


    ! Open the input history header file and read some of the info.

    if (mchnum==master_num) then !io-process only

       write(*,*)'chemdata_in=',chemdata_in
       !pause

       nc       = LEN_trim(chemdata_in)
       hnameinh = chemdata_in(1:nc-9)//'.vfm'

       call rams_f_open(iunhd, chemdata_in(1:len_trim(chemdata_in)), 'FORMATTED', 'OLD', 'READ', 0)

       ie      = cio(iunhd, 1, 'ngrids',  ngrids1)
       ngridsh = ngrids1
       ie      = cio(iunhd, 1, 'nnxp',    nnxp1(1:ngrids1))
       ie      = cio(iunhd, 1, 'nnyp',    nnyp1(1:ngrids1))
       ie      = cio(iunhd, 1, 'nnzp',    nnzp1(1:ngrids1))
       ie      = cio(iunhd, 1, 'npatch',  npatch1)
       ie      = cio(iunhd, 1, 'nzg',     nzg1)
       ie      = cio(iunhd, 1, 'nzs',     nzs1)
       ie      = cio(iunhd, 1, 'ioutput', ioutput1)

       ! Find maximum size of any array on history file.
       ! Allocate scratch array of this size.

       maxarr = 0
       do ngr=1,ngridsh
          maxarr = max(                            &
               maxarr,                             &
               nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr),   &
               nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1, &
               nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1  )
       enddo

    endif

    ! Broadcasting Data
    call parf_bcast(ngrids1, master_num)
    call parf_bcast(nnxp1, int(size(nnxp1,1),i8), master_num)
    call parf_bcast(nnyp1, int(size(nnyp1,1),i8), master_num)
    call parf_bcast(nnzp1, int(size(nnzp1,1),i8), master_num)
    call parf_bcast(npatch1, master_num)
    call parf_bcast(nzg1, master_num)
    call parf_bcast(nzs1, master_num)
    call parf_bcast(ioutput1, master_num)
    call parf_bcast(maxarr, master_num)

    ! read stuff here

    call hist_pol_read(maxarr, hnameinh(1:len_trim(hnameinh)), iunhd, &
         oneVarTable, oneVarTableSize)

    if (mchnum==master_num) then !io-process only
       print*, 'back from read'
       close(iunhd)
    endif


  end subroutine init_conc_prev


  !******************************************************************************

  subroutine hist_pol_read(maxarr, hnamein, iunhd, oneVarTable, oneVarTableSize)
    include "constants.h"
    include "files.h"

    ! Arguments:
    integer, intent(in)           :: maxarr
    character(len=f_name_length), intent(in) :: hnamein
    integer, intent(in)           :: iunhd
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(in) :: oneVarTableSize

    ! Local variables:
    integer            :: ngr, npts, nptsh, nv, nvh, i
    real, allocatable  :: scr(:)
    integer, parameter :: inhunt=10
    type(head_table), allocatable :: hr_table(:)
    integer            :: iCopyFlg


    !  Read variable header info
    if (mchnum==master_num) then !io-process only
       allocate (scr(maxarr))
       rewind(iunhd)
       read(iunhd,*) nvbtab
       allocate (hr_table(nvbtab))
       do nv=1,nvbtab
          read(iunhd,*)                &
               hr_table(nv)%string,    &
               hr_table(nv)%npointer,  &
               hr_table(nv)%idim_type, &
               hr_table(nv)%ngrid,     &
               hr_table(nv)%nvalues
       enddo

    endif

    ! Broadcasting data
    call parf_bcast(nvbtab, master_num)

    if (mchnum==master_num) then !io-process only
       ! Open history data file
       call rams_f_open(inhunt, hnamein(1:len_trim(hnamein)), 'UNFORMATTED', 'OLD', 'READ', 0)
    endif

    ! Loop through all variables
    do nvh=1,nvbtab

       if (mchnum==master_num) then !io-process only
          ! Read a variable
          nptsh = hr_table(nvh)%nvalues
          read(inhunt) (scr(i), i=1,nptsh)

          !  See if this variable is active in the current run
          ngr = hr_table(nvh)%ngrid
          if (ngr>ngrids) cycle

          do nv=1,oneVarTableSize
             npts = oneVarTable(nv)%npts
             if (hr_table(nvh)%string==oneVarTable(nv)%name) then
                if (nptsh/=npts) then
                   print*, 'Grid point number mismatch on history field:',  &
                        oneVarTable(nv)%name,npts,nptsh
                   call fatal_error('History read number points error')
                endif

                iCopyFlg = 0

                if(  oneVarTable(nv)%name=='PNO'   .or. &
                     oneVarTable(nv)%name=='PNO2'  .or. &
                     oneVarTable(nv)%name=='PPM25' .or. &
                     oneVarTable(nv)%name=='PCO'   .or. &
                     oneVarTable(nv)%name=='PVOC'  .or. &
                     oneVarTable(nv)%name=='PSO2'  .or. &
                     oneVarTable(nv)%name=='PSO4'  .or. &
                     oneVarTable(nv)%name=='PAER'  .or. &
                     oneVarTable(nv)%name=='PVOC'  .or. &
                     oneVarTable(nv)%name=='PSO2'  .or. &
                     oneVarTable(nv)%name=='PO3'   .or. &
                     oneVarTable(nv)%name=='PRHCO' .or. &
                     oneVarTable(nv)%name=='PHO2'  .or. &
                     oneVarTable(nv)%name=='PO3P'  .or. &
                     oneVarTable(nv)%name=='PO1D'  .or. &
                     oneVarTable(nv)%name=='PHO'   .or. &
                     oneVarTable(nv)%name=='PROO'       ) then

                   write (UNIT=6, FMT='(a25,2i5,3x,a18,i10)') &
                        'Polutants History filling grid: ',   &
                        ngr, nv, oneVarTable(nv)%name, npts
                   iCopyFlg = 1
                   exit
                endif

             endif
          enddo

       endif

       ! Broadcasting last array read
       if (iCopyFlg==1) then
          call parf_bcast(nv, master_num)
          call parf_bcast(ngr, master_num)
          call parf_bcast(npts, master_num)
          call parf_bcast(scr, int(npts, i8), master_num)

          if (oneVarTable(nv)%idim_type == 2) then
             call atob(npts, scr, oneVarTable(nv)%var_p_2D)
          else if (oneVarTable(nv)%idim_type == 3) then
             call atob(npts, scr, oneVarTable(nv)%var_p_3D)
          else if (oneVarTable(nv)%idim_type == 4) then
             call atob(npts, scr, oneVarTable(nv)%var_p_4D)
          else if (oneVarTable(nv)%idim_type == 5) then
             call atob(npts, scr, oneVarTable(nv)%var_p_4D)
          else if (oneVarTable(nv)%idim_type == 6) then
             call atob(npts, scr, oneVarTable(nv)%var_p_3D)
          else if (oneVarTable(nv)%idim_type == 7) then
             call atob(npts, scr, oneVarTable(nv)%var_p_3D)
          end if
       endif

    enddo

    if (mchnum==master_num) then !io-process only
       ! Close the input history file
       close(inhunt)
       deallocate(scr, hr_table)
    endif

  end subroutine hist_pol_read
end module ModGasPart
