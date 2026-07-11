!==========================================================================================
!- GFDL cloud microphysics 
!- Adapted to BRAMS 5.0+ by Saulo Freitas 2018/2021
!==========================================================================================
module ModMicGfdlDriver

  use ModNamelistFile, only: &
       NamelistFile
  
  use ModBasicFields, only: &
       BasicFields

  use ModMicroFields, only: &
       MicroFields
  
  use mem_grid, only:   &
       grid_vars, &
       ngrids,          & ! INTENT(IN)
       ngrid,           & ! INTENT(IN)
       zm,              & ! INTENT(IN)
       dzt,             & ! INTENT(IN)
       dtlt,            & ! INTENT(IN)
       jdim,            & ! INTENT(IN)
       maxnzp,          & ! INTENT(IN)
       deltaxn,         & ! INTENT(IN)
       deltayn,         & ! INTENT(IN)
       time,            & ! INTENT(IN)
       zt,              & ! INTENT(IN)
       itime1,          & ! INTENT(IN)
       if_adap,         & ! INTENT(IN)
       grid_g,          & ! INTENT(IN)
       nnzp, &
       npatch, &
       imonth1, &
       dtlongn, &
       timmax

  use node_mod, only:  &
       mzp,             & ! INTENT(IN)
       mxp,             & ! INTENT(IN)
       myp,             & ! INTENT(IN)
       ja,              & ! INTENT(IN)
       jz,              & ! INTENT(IN)
       ia,              & ! INTENT(IN)
       iz,              & ! INTENT(IN)
       mynum, &
       i0, &
       j0        ! INTENT(IN)

  use gfdl_cloud_microphys_mod, only: &
       gfdl_cloud_microphys_driver,  &
       gfdl_cloud_microphys_init

  use rconstants, only: &
       p00, &
       cp, &
       cpor, &
       alvl, &
       alvi, &
       cpi, &
       cpi4, &
       cp253i

  use mem_leaf, only: &
       leaf_g, &
       isfcl

  implicit none

  private
  public :: micro_gfdl

contains

  subroutine micro_gfdl(oneNamelistFile, oneBasicFields, oneMicroFields)
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    
    integer :: ims, ime, jms, jme, kms, kme  

    !---bounds for the GFDL microphysics
    ims = 1; ime = mxp
    jms = 1; jme = myp 
    kms = 1; kme = mzp-2

    call brams_to_mic_gfdl(ia,ja,iz,jz, &
         oneNamelistFile%ilwrtyp,     &
         oneNamelistFile%iswrtyp     &
         ,ims, ime, jms, jme, kms, kme&
         ,mzp    &
         ,ngrid    &
         ,mynum    &
                                !
         ,deltaxn(ngrid) &   ! INTENT(IN)
         ,deltayn(ngrid) &   ! INTENT(IN)
                                !
         ,dtlt  &
         ,time  &
         ,zm  &
         ,dzt  &       
         ,zt  &
         ,oneBasicFields &
         ,grid_g (ngrid) &
         ,oneMicroFields &
                                !
         )

    !- for consistency with surface and radiation schemes, the total
    !- precip will be also stored in the pcpg array
    oneMicroFields%pcpg(:,:)=oneMicroFields%pcprr(:,:)

  end subroutine micro_gfdl


  
  !=======================================================================================
  !
  subroutine brams_to_mic_gfdl(ia,ja,iz,jz  &
       ,ilwrtyp     &
       ,iswrtyp     &
       ,ims, ime, jms, jme, kms, kme   &
       ,m1    &
       ,ngrid &
       ,mynum &
       ,deltaxn&   ! INTENT(IN)
       ,deltayn&   ! INTENT(IN)
       ,dtlt  &
       ,time  &
       ,zm    &
       ,dzt   &
       ,zt    &
       ,oneBasicFields &
       ,grd &
       ,oneMicroFields &
       )

    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(grid_vars)  ::grd
    type(MicroFields), pointer, intent(in) :: oneMicroFields

    integer, intent(IN) ::  &          
         ilwrtyp       &
         ,iswrtyp       &
         ,ims, ime, jms, jme, kms, kme   &
         ,m1    &
         ,ngrid &
         ,mynum &
         ,ia,ja,iz,jz

    real, intent(IN) ::   &
         dtlt        &
         ,time        &
         ,deltaxn     &  
         ,deltayn   

    real,  intent(IN)   ,dimension(m1) :: &
         zm    &
         ,dzt   &
         ,zt    

    real, dimension( m1, ims:ime , jms:jme )  ::   &
         exner


    real, dimension( ims:ime , jms:jme, kms:kme )  ::   &
         temp,dz,w,u,v, dp

    real, dimension( ims:ime, jms:jme , kms:kme) ::                  &
         qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr , &
         re_cloud, re_ice ,re_snow, rad_cf , &
         DQVDT_micro , & 
         DQLDT_micro , & 
         DQRDT_micro , & 
         DQIDT_micro , & 
         DQSDT_micro , & 
         DQGDT_micro , & 
         DQADT_micro , & 
         DUDT_micro , & 
         DVDT_micro , & 
         DTDT_micro , &
         NACTLI     , &
         REV_MC_X, RSU_MC_X, EVAPC_X

    real, dimension( ims:ime , jms:jme )  ::  &
         RAINNC     &
         ,RAINNCV    &
         ,SNOWNC     &
         ,SNOWNCV    &
         ,GRAUPELNC  &
         ,GRAUPELNCV 

    real, dimension( ims:ime , jms:jme )  ::  &
         area    &
         ,frland  &
         ,cnv_fraction

    real, dimension(ims:ime , jms:jme ,0:kme) :: &
         PFI_LS_X &
         ,PFL_LS_X 

    real, dimension(ims:ime , jms:jme ) :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
    !----------------------------------------------------------------------
    ! qv              water vapor    mixing ratio (kg/kg)
    ! qc              cloud water    mixing ratio (kg/kg)
    ! qr              rain water     mixing ratio (kg/kg)
    ! qi              cloud ice      mixing ratio (kg/kg)
    ! qs              snow            mixing ratio (kg/kg)
    ! qg              graupel            mixing ratio (kg/kg)
    !
    ! qnc             cloud water number concentration (#/kg)
    ! qni             cloud ice   number concentration (#/kg)
    ! qnr             rain        number concentration (#/kg)
    ! qnwfa      water friendly aerosol number concentration (#/kg) - CCN
    ! qnifa      ice   friendly aerosol number concentration (#/kg) - IFN
    !
    !-- th            potential temperature    (K)
    !-- w             vertical velocity (cartesian) (m/s)
    !-- rho           density of air           (kg/m^3)
    !-- press         pressure                 (Pa)
    !-- RAINNC        grid scale precipitation (mm)
    !-- RAINNCV       one time step grid scale precipitation (mm/step)
    !-- SNOWNC        grid scale snow and ice (mm)
    !-- SNOWNCV       one time step grid scale snow and ice (mm/step)
    !-- GRAUPELNC     grid scale graupel (mm)
    !-- GRAUPELNCV    one time step grid scale graupel (mm/step)
    !-- HAILNC        grid scale hail (mm)
    !-- HAILNCV       one time step grid scale hail (mm/step)
    !-- z             Height above sea level   (m)
    !-- dt            Time step              (s)


    real :: tempk,rliq,rice,til,qhydm,tairstr,exner_tmp,rcond
    real :: dx, dy    ! grid spacing (m)
    real :: dt_moist  ! model timestep (s)
    logical :: start_of_simulation =.true. 
    integer, allocatable, dimension(:),save:: flip
    real   , dimension(m1):: dp_tmp

    logical  :: lhydrostatic   =.false.
    logical  :: lphys_hydrostatic =.false.
    integer :: i,j,k,m
    real, parameter :: ANV_ICEFALL=1.0, LS_ICEFALL=1.0


    if(start_of_simulation) then
       !-initialization
       call gfdl_cloud_microphys_init()
       !- define the vector "flip" to invert the z-axis orientation
       allocate(flip(m1)); flip(:)=-9999
       call flipz_minus1(flip,kme+1)
       !-
       start_of_simulation =.false.
    endif

    !---- zero-out microphysics tendencies
    DQVDT_micro = 0.
    DQLDT_micro = 0.
    DQRDT_micro = 0.
    DQIDT_micro = 0.
    DQSDT_micro = 0.
    DQGDT_micro = 0.
    DQADT_micro = 0.
    DUDT_micro = 0.
    DVDT_micro = 0.
    DTDT_micro = 0.
    !--- zero-out 3D Precipitation Fluxes 
    PFI_LS_X = 0. ! Ice
    PFL_LS_X = 0. ! Liquid
    !--- output rain re-evaporation and sublimation
    REV_MC_X= 0. ;RSU_MC_X= 0. ;EVAPC_X= 0. 

    dt_moist= dtlt  ! time step   (s)

    do j = jms,jme
       do i = ims,ime 

          area      (i,j) = deltaxn*deltayn
          frland   (i,j) = 1.-leaf_g(ngrid)%patch_area(i,j,1)
          cnv_fraction(i,j) = 0.1

          do k=1,kme
             qc_curr (i,j,k)= max(0.0,oneMicroFields%rcp(flip(k),i,j))      ! QC     
             qr_curr (i,j,k)= max(0.0,oneMicroFields%rrp(flip(k),i,j))      ! QR   
             qi_curr (i,j,k)= max(0.0,oneMicroFields%rpp(flip(k),i,j))      ! QI   
             qs_curr (i,j,k)= max(0.0,oneMicroFields%rsp(flip(k),i,j))      ! QS   
             qg_curr (i,j,k)= max(0.0,oneMicroFields%rgp(flip(k),i,j))      ! QG   

             rad_cf  (i,j,k)= max(0.0,oneMicroFields%cldfr(flip(k),i,j))      ! cloud fraction

             rcond = qc_curr (i,j,k) + qr_curr (i,j,k) + &
                  qi_curr (i,j,k) + qs_curr (i,j,k) + qg_curr (i,j,k)

             qv_curr (i,j,k)= max(1.e-12, oneBasicFields%rtp(flip(k),i,j) - rcond)      

             W       (i,j,k)= oneBasicFields%wp(flip(k),i,j)    ! vertical velocity (m/s) 
             U       (i,j,k)= oneBasicFields%up(flip(k),i,j)    ! U velocity (m/s) 
             V       (i,j,k)= oneBasicFields%vp(flip(k),i,j)    ! V velocity (m/s) 

             ! Delta-Z layer thickness (gfdl expects this to be negative)
             DZ      (i,j,k)= - grd%rtgt(i,j)/dzt(flip(k)) ! layer thickness (m) 

          enddo
          !
          !--special treatment for pressure layers
          do k=1,kme+2     ! note that kme = mzp-2
             ! Exner function (brams data structure)
             exner  (k,i,j)= ((oneBasicFields%pp(k,i,j)+oneBasicFields%pi0(k,i,j))*cpi) 
          enddo
          do k=1,kme      
             ! pressure thickness (Pa)
             dp_tmp(k) = 0.5*(exner(k+2,i,j)** cpor - exner(k,i,j)** cpor) * p00
             !if(mynum==1 .and. i==ims .and. j==jms) print*,"x=",dp_tmp(k) !, press(k+2,i,j) ,press(k,i,j),k
          enddo
          m=kme
          do k=1,kme      
             DP(i,j,k) = - dp_tmp(m) ; m=m-1 !invert press
          enddo
          !
          !--special treatment for CCN/ICN numbers conce
          do k=1,kme             
             NACTLI (i,j,k)= 3.e+2 !cm-2 !???????????????
             !qni_curr(i,j,k)= max(0.0,oneMicroFields%cpp(flip(k),i,j))      ! NI  !??????????????? 
             !qnr_curr(i,j,k)= max(0.0,oneMicroFields%crp(flip(k),i,j))      ! NR  !???????????????
          enddo
          !
          !- get  temperature (temp) from theta_il (thp) and condensates
          do k=1,kme
             exner_tmp = exner (flip(k),i,j)
             tempK   = oneBasicFields%theta(flip(k),i,j)*exner_tmp
             til     = oneBasicFields%thp  (flip(k),i,j)*exner_tmp

             rliq=  qc_curr(i,j,k) + qr_curr(i,j,k)  
             rice=  qi_curr(i,j,k) + qs_curr(i,j,k) + qg_curr(i,j,k)
             qhydm=  alvl * rliq + alvi * rice

             if (tempK .gt. 253.) then
                tairstr = 0.5 * (til + sqrt(til * (til + cpi4 * qhydm)))
             else
                tairstr = til * (1. + qhydm * cp253i)
             endif
             !- updated  temperature TEMP in Kelvin (adv+dif+rad+conv+)
             TEMP (i,j,k) = tairstr
          enddo
       enddo
    enddo


    !--- Execute GFDL microphysics
    call gfdl_cloud_microphys_driver( &
                                ! input water/cloud species and liquid+ice CCN [NACTL+NACTI]
                                ! RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, (NACTL+NACTI)/1.e6, &
         qv_curr,   &! QV=qv_curr,     
         qc_curr,   &! QC=qc_curr,     
         qr_curr,   &! QR=qr_curr,     
         qi_curr,   &! QI=qi_curr,     
         qs_curr,   &! QS=qs_curr,     
         qg_curr,   &! QG=qg_curr,     
         rad_cf,    &! cloud fraction  
         NACTLI,    &!@NR=qnr_curr, 
                                ! Output tendencies
         DQVDT_micro, DQLDT_micro, DQRDT_micro, DQIDT_micro, &
         DQSDT_micro, DQGDT_micro, DQADT_micro, DTDT_micro, &
                                ! Input fields
         TEMP, W, U, V, DUDT_micro, DVDT_micro, DZ, DP, &
                                ! constant inputs
         AREA, dt_moist, FRLAND, CNV_FRACTION, &
         
         ANV_ICEFALL, LS_ICEFALL, &
                                ! Output rain re-evaporation and sublimation
         REV_MC_X, RSU_MC_X, &   ! EVAPC_X, & 
         
                                ! Output precipitates
         PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
                                ! Output mass flux during sedimentation (Pa kg/kg)
         PFL_LS_X(:,:,1:kme), PFI_LS_X(:,:,1:kme), &
                                ! constant grid/time information
         LHYDROSTATIC, LPHYS_HYDROSTATIC, &
         ims,ime, jms,jme, 1,kme, 1, kme, &
         re_cloud, re_ice, re_snow)

    !-------------------- print ---------->
    ! if(mynum==10000 ) then
    !         print*,"-----------------------------------------------------------------------"
    !         print*,"WX=",maxval(V), maxval(W),  maxval(U),mynum
    !         print*,"WN=",minval(V), minval(W),  minval(U),mynum
    !         print*,"TX=",maxval(TEMP), maxval(qv_curr),  maxval(qc_curr),mynum
    !         print*,"TN=",minval(TEMP), minval(qv_curr),  minval(qc_curr),mynum
    !         print*,"PX=",maxval(dp), maxval(dz),  maxval(qc_curr),mynum
    !         print*,"PN=",minval(dp), minval(dz),  minval(qc_curr),mynum
    !         call flush(6)
    ! endif
    !-------------------- print ----------<


    !--- update variables after cloud microphysics processes (only in the halo)
    do j = ja,jz
       do i = ia,iz

          !- surface precipitation (units are kg/m^2 = mm and mm/s)
          !- accpr,pcprr   kg/m2 - rain+ice+snow+graupel ! p = for each dt  (or per time step)
          oneMicroFields%pcprr(i,j) = dt_moist*(PRCP_RAIN(i,j) + PRCP_SNOW   (i,j) + &
               PRCP_ICE (i,j) + PRCP_GRAUPEL(i,j)) / 86400. 

          oneMicroFields%accpr(i,j) = oneMicroFields%accpr(i,j) + oneMicroFields%pcprr(i,j)  

          !- accps,pcprs  kg/m2 - ice+snow
          oneMicroFields%pcprs(i,j) = dt_moist*(PRCP_SNOW(i,j) + PRCP_ICE(i,j))/ 86400.    
          oneMicroFields%accps(i,j) = oneMicroFields%accps(i,j) + oneMicroFields%pcprs(i,j)

          !- accpg, pcprg  &! kg/m2 - graupel
          oneMicroFields%pcprg(i,j) = dt_moist* PRCP_GRAUPEL(i,j)/ 86400.
          oneMicroFields%accpg(i,j) = oneMicroFields%accpg(i,j) + oneMicroFields%pcprg(i,j)


          do k=kme,1,-1
             !-srf    
             !if(i==89 .and. j==22 .and. minval(oneMicroFields%rcp(:,i,j))< 0.)  then
             !      print*,'rrp',k, oneMicroFields%rcp(flip(k),i,j), qi_curr(i,j,k) , qi_curr(i,j,k) + DQIDT_micro(i,j,k) * DT_MOIST
             !      print*,'rv ',k, oneBasicFields%rv(flip(k),i,j), qv_curr(i,j,k), qv_curr(i,j,k) + DQVDT_micro(i,j,k) * DT_MOIST
             !call flush(6)
             !endif
             !-srf 
             oneMicroFields%rcp(flip(k),i,j)= max(0., qc_curr(i,j,k) + DQLDT_micro(i,j,k) * DT_MOIST)
             oneMicroFields%rrp(flip(k),i,j)= max(0., qr_curr(i,j,k) + DQRDT_micro(i,j,k) * DT_MOIST)
             oneMicroFields%rpp(flip(k),i,j)= max(0., qi_curr(i,j,k) + DQIDT_micro(i,j,k) * DT_MOIST)
             oneMicroFields%rsp(flip(k),i,j)= max(0., qs_curr(i,j,k) + DQSDT_micro(i,j,k) * DT_MOIST)
             oneMicroFields%rgp(flip(k),i,j)= max(0., qg_curr(i,j,k) + DQGDT_micro(i,j,k) * DT_MOIST)

             oneMicroFields%cldfr(flip(k),i,j) = max(0., rad_cf(i,j,k) + DQADT_micro(i,j,k) * DT_MOIST)

             !- cloud ice and liq effective radius,  RRTMG requires in micrometer
             oneMicroFields%rei (flip(k),i,j)= re_ice  (i,j,k)  
             oneMicroFields%rel (flip(k),i,j)= re_cloud(i,j,k)
             !oneMicroFields%snow(flip(k),i,j)= re_snow (i,j,k)  !-- srf check this latter.


             oneBasicFields%rv(flip(k),i,j)= max(1.e-12, qv_curr(i,j,k) + DQVDT_micro(i,j,k) * DT_MOIST)

             tempK                =   temp(i,j,k) + DTDT_micro (i,j,k) * DT_MOIST

             oneBasicFields%theta(flip(k),i,j) = tempK/ exner(flip(k),i,j)

          enddo

          do k=2,kme
             !- update liq-ice potential temperature THP in Kelvin including microphysics processes
             rliq     = oneMicroFields%rcp(k,i,j) + oneMicroFields%rrp(k,i,j)   
             rice     = oneMicroFields%rpp(k,i,j) + oneMicroFields%rsp(k,i,j) + oneMicroFields%rgp(k,i,j)  

             tempK    = oneBasicFields%theta(k,i,j) * exner(k,i,j)

             oneBasicFields%thp(k,i,j) = oneBasicFields%theta(k,i,j)*(1. + alvl * rliq/(cp * max(tempK,253.))  &
                  + alvi * rice/(cp * max(tempK,253.)) ) **(-1.0)
             !- update total water 
             oneBasicFields%rtp(k,i,j) = oneBasicFields%rv(k,i,j) + rliq + rice
          enddo
          !reserved   
          !  cpp(k)= qni_curr(1,k,1)
          !  crp(k)= qnr_curr(1,k,1)
          !         U1   = U1   + DUDT_micro  * DT_MOIST
          !         V1   = V1   + DVDT_micro  * DT_MOIST   
          !         W1 was updated in gfdl_microphsycis, fill WI tendency export
          !         WI =  (W1 - W)/DT_MOIST
          !         RAD_CF = RAD_CF + DQADT_micro * DT_MOIST
          !reserved

          !- definition for k=1
          oneBasicFields%rtp(1,i,j)  = oneBasicFields%rtp(2,i,j)  
          oneMicroFields%rcp  (1,i,j)  = oneMicroFields%rcp  (2,i,j)  
          oneMicroFields%rrp  (1,i,j)  = oneMicroFields%rrp  (2,i,j)  
          oneMicroFields%rpp  (1,i,j)  = oneMicroFields%rpp  (2,i,j)  
          oneMicroFields%rsp  (1,i,j)  = oneMicroFields%rsp  (2,i,j)  
          oneMicroFields%rgp  (1,i,j)  = oneMicroFields%rgp  (2,i,j)  
          oneMicroFields%cldfr(1,i,j)  = oneMicroFields%cldfr(2,i,j)

          oneBasicFields%rv   (1,i,j) = oneBasicFields%rv   (2,i,j)  
          oneBasicFields%theta(1,i,j) = oneBasicFields%theta(2,i,j)
          oneBasicFields%thp  (1,i,j) = oneBasicFields%thp  (2,i,j)

          !- cloud ice and liq effective radius 
          oneMicroFields%rei(1,i,j)= oneMicroFields%rei(2,i,j) 
          oneMicroFields%rel(1,i,j)= oneMicroFields%rel(2,i,j)
          !oneMicroFields%snow(1,i,j)=oneMicroFields%snow(2,i,j)

          !temporary 
          !oneMicroFields%rei   (1:m1,i,j) =5.01E-6 * 1.e+6 ! RRTM requires in micrometer
          !oneMicroFields%rel   (1:m1,i,j) =2.51E-6 * 1.e+6 ! RRTM requires in micromete


          !reserved
          !  oneMicroFields%cpp (1,i,j)  = oneMicroFields%cpp (2,i,j)  
          !  oneMicroFields%crp (1,i,j)  = oneMicroFields%crp (2,i,j)  
          !  oneMicroFields%ccp (1,i,j)  = oneMicroFields%ccp (2,i,j) 

       enddo
    enddo

  end  subroutine brams_to_mic_gfdl



  
  !------------------------------------------------------------------------------------

  subroutine flipz_minus1(flip,mzp)
    integer, intent(In) :: mzp
    integer, dimension(mzp), intent(inout) :: flip
    integer :: m,k
    m=mzp
    do k=1,mzp
       flip(k)=m
       m=m-1
    enddo
    flip(mzp)=1
  end subroutine flipz_minus1
end module ModMicGfdlDriver
