!##############################################################
!
! Programmed by  Enio Pereira de Souza
! Adapted    by  Alvaro Luiz Fazenda    (for V.5.04)
! Revised    by  Saulo Freitas (for BRAMS version 5.2/2015)
!##############################################################
module ModRShCuPar

  use ModBasicFields, only : &
       BasicFields

  use ModMicroFields, only : &
       MicroFields

  use ModShcuFields, only: &
       ShcuFields
  
  use mem_grid, only : &
       time,    &   ! INTENT(IN)
       zt,      &   !INTENT(IN)
       zm,      &   !INTENT(IN/OUT)
       ngrid,   &   ! INTENT(IN)
       TIME,    &   ! INTENT(IN)
       INITIAL, &   ! INTENT(IN)
       dtlt,    &   ! INTENT(IN)
       grid_g       ! INTENT(IN)

  use ModTurbFields, only: &
       TurbFields

  use mem_tend, only: &
       tend

  use node_mod, only : &
       mxp,   &   ! intent(in)
       myp,                   &   ! intent(in)
       mzp,                   &   ! intent(in)
       ia,                    &   ! intent(in)
       iz,                    &   ! intent(in)
       ja,                    &   ! intent(in)
       jz,                    &   ! intent(in)
       i0,                    &   ! intent(in)
       j0                         ! intent(in)

  use shcu_vars_const, only : &
       shcufrq, &       ! intent(in)
       entf,   &   ! intent(out)
       alhf,                    &   ! intent(out)
       qvcon,                   &   ! intent(out)
       akvd,                    &   ! intent(out)
       cl_con,                  &   ! intent(out)
       dtdt,                    &   ! intent(in)
       wc,                      &   ! intent(in)
       drdt                         ! intent(IN)

  use ModConvComs, only : &
       nkp,    &   ! intent(in)  ! parameter
       zc,                      &   ! intent(out)
       ze,                      &   ! intent(out)
       pke,                     &   ! intent(out)
       the,                     &   ! intent(in)
       thve,                    &   ! intent(out)
       te,                      &   ! intent(out)
       pe,                      &   ! intent(out)
       rhoe,                    &   ! intent(out)
       dzlow,                   &   ! intent(out) ! talvez var.local
       dzhigh,                  &   ! intent(out) ! talvez var.local
       zmid,                    &   ! intent(out) ! talvez var.local
       cdzmin,                  &   ! intent(out) ! rever funcao.
       zcon,                    &   ! intent(in)
       kmt,                     &   ! intent(out)
       wcon,                    &   ! intent(in/out)
       zzcon,                   &   ! intent(in)
       wpe,                     &   ! intent(in/out)
       thtcon,                  &   ! intent(in)
       picon,                   &   ! intent(in)
       kcon,                    &   ! intent(out)  ! talvez var.local
       tlcl,                    &   ! intent(in/out)
       plcl,                    &   ! intent(in/out)
       dzlcl,                   &   ! intent(in/out)
       klcl,                    &   ! intent(out)
       igo, &                          ! intent(out)
       qvct1,                 &   ! intent(out)
       qvct2,                 &   ! intent(out)
       qvct3,                 &   ! intent(out)
       qvct4,                 &   ! intent(out)
       dncon,                   &   ! intent(out)
       icprtfl,                 &   ! intent(out)   ! maybe local variable
       icpltfl

  use ModRConv, only: &
       lcl, &
       vertmap2

  use shcu_vars_const, only : &
       dtdt, &   ! intent(in/out)
       alvl,                        &   ! intent(in) ! parameter
       drdt, &                             ! intent(in/out)
       dscv, &   ! intent(out)
       dsc,                         &   ! intent(in)
       cp,                          &   ! intent(in) ! parameter
       qvc,                         &   ! intent(in)
       wlc,                         &   ! intent(in)
       dsev,                        &   ! intent(out)
       dse,                         &   ! intent(in)
       qve,                         &   ! intent(in)
       dsc0v,                       &   ! intent(out) ! maybe local var.?
       dsc0,                        &   ! intent(in)
       dsc0vm,                      &   ! intent(out) ! maybe local var.?
       entf,                        &   ! intent(in)
       g,                           &   ! intent(in)  ! parameter
       cape,                        &   ! intent(out)
       delz,                        &   ! intent(in)
       ktop, &                             ! intent(out)
       efic,                        &   ! intent(out) ! maybe local var.?
       alhf,                        &   ! intent(in)
       dcape,                       &   ! intent(out) ! maybe local var.?
       tcape,                       &   ! intent(out) ! maybe local var.?
       wc, &                               ! intent(out)
       dqdt,                        &   ! intent(out) ! maybe local var.?
       uhe,                     &   ! intent(out)
       evaps,                   &   ! intent(out) ! maybe local var.?
       qvse,                    &   ! intent(out) ! maybe local var.?
       uhes,                    &   ! intent(out)
       rhe,                     &   ! intent(out) ! **
       gamma,                   &   ! intent(out)
       uhc,                     &   ! intent(out)
       dldzby2,                 &   ! intent(out) ! maybe local var.?
       qvcon,                   &   ! intent(in/out)
       akvd,                    &   ! intent(in/out)
       cl_con,                  &   ! intent(in/out)
       cl_pe,                   &   ! intent(in/out)
       cpr,                     &   ! intent(in)  ! parameter
       p00,                     &   ! intent(in)  ! parameter
       r,                       &   ! intent(in)  ! parameter
       kzi,                     &   ! intent(out) ! maybe local var.?
       akvde                        ! intent(in/out)

  use mem_scratch, only : &
       vctr5,  & ! intent(in/out)
       vctr6                        ! intent(in/out)

  implicit none

  include "constants.h"

  private
  public :: shcupa
  
contains

  
  subroutine shcupa(oneBasicFields, oneTurbFields, oneMicroFields, oneShcuFields)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    type(ShcuFields), pointer, intent(in) :: oneShcuFields


    real    :: cptime = 0. !7200.

    integer :: i, j

    !  if(initial.eq.2.and.time.lt.cptime-dtlt)return

    if(mod(time+dtlt+.001,shcufrq).le.DTLT.or.time.lt..01)then

       oneShcuFields%THSRCSH = 0.
       oneShcuFields%RTSRCSH = 0.
       oneShcuFields%SHMF    = 0.

       call SHCUPAR(mzp,mxp,myp,ia,iz,ja,jz,i0,j0,                   &
            oneBasicFields%wp, oneBasicFields%theta,   &
            oneBasicFields%pp, oneBasicFields%pi0,     &
            oneBasicFields%dn0, oneBasicFields%rv,     &
            oneShcuFields%THSRCSH,                            &
            oneShcuFields%RTSRCSH,                            & 
            oneShcuFields%SHMF, grid_g(ngrid)%rtgt,        &
            oneTurbFields%sflux_t,                              & 
            oneTurbFields%sflux_r,                              &
            oneTurbFields%vkh,                                &
            oneMicroFields%rcp)

    endif

    call ACCUM(int(mxp*myp*mzp,i8), tend%tht, oneShcuFields%THSRCSH)
    call ACCUM(int(mxp*myp*mzp,i8), tend%rtt, oneShcuFields%RTSRCSH)

    return
  end subroutine SHCUPA

  !##############################################################
  subroutine SHCUPAR(m1,m2,m3,ia,iz,ja,jz,i0,j0,                  &
       WP,THETA,PP,PI0,DN0,RV,                      &
       THSRCSH,RTSRCSH,SHMF,RTGT,TFZ,QFZ,KHV,RCLOUD)

    integer, intent(IN) :: m1, m2, m3, ia, iz, ja, jz, i0, j0
    real, intent(IN)    :: WP(m1,m2,m3), THETA(m1,m2,m3),        &
         PP(m1,m2,m3), PI0(m1,m2,m3), DN0(m1,m2,m3),             &
         RV(m1,m2,m3)
    real, intent(OUT)   :: THSRCSH(m1,m2,m3), RTSRCSH(m1,m2,m3), &
         SHMF(m2,m3)
    real, intent(IN)    :: RTGT(m2,m3), TFZ(m2,m3), QFZ(m2,m3),  &
         KHV(m1,m2,m3), RCLOUD(m1,m2,m3)

    integer :: ICPCNT = 0
    integer :: IPRTFRQ, I, J, K

    logical,parameter :: make_conv_grid_interpol=.false.

    !  TFZMIN = TFZ(1,1)
    !  TFZMAX = TFZ(1,1)


    ICPRTFL=0
    IPRTFRQ=8
    ICPLTFL=0
    ICPCNT=ICPCNT+1
    if(mod(ICPCNT-IPRTFRQ+1,IPRTFRQ).eq.0) then
       ICPRTFL=1
    endif

    do J=ja,jz
       do I=ia,iz
          !        
          !do initialization
          SHMF(i,j)=0.0
          THSRCSH(:,I,J)=0.
          RTSRCSH(:,I,J)=0.
          !
          ! Definition of the enthalpy flux ENTF (K*m/s) and
          ! latent energy flux (m/s)*(g/g)
          !
          ENTF=min(.50,TFZ(I,J))
          ALHF=min(.00028,QFZ(I,J))
          !
          if(ENTF.le.0.0) cycle
          !
          do K=1,m1
             WCON(K)=WP(K,I,J)
             THTCON(K)=THETA(K,I,J)
             PICON(K)=(PP(K,I,J)+PI0(K,I,J))
             DNCON(K)=DN0(K,I,J) 
             !EPS
             !At this point, we change rv by qv as required by the
             !basic equations.
             !EPS            
             QVCON(K) =RV(K,I,J)/(1.0+RV(K,I,J))            
             ZCON(K)  =ZT(K) *RTGT(I,J)
             ZZCON(K) =ZM(K) *RTGT(I,J)
             AKVD(K)  =KHV(K,I,J)
             CL_CON(K)=RCLOUD(K,I,J)
          enddo
          !
          IGO=1
          !
          call SHCU_ENV(m1-1,make_conv_grid_interpol) 
          !
          if(IGO.ne.0) call CL_TOP
          ! 
          if(IGO.ne.0) call W_SHALLOW(I,J,TIME)
          !
          if(IGO.ne.0) call SH_RATES
          ! 
          if(IGO.ne.0) then
             !
             if (make_conv_grid_interpol) call SH2MOD(m1)

             do K=2,m1-1
                THSRCSH(K,I,J)=DTDT(K)
                RTSRCSH(K,I,J)=DRDT(K) 
             enddo
             SHMF(i,j) = DNCON(KLCL)*WC(KLCL)
             !
          endif
          !
          !
       enddo
    enddo
    !
    return
  end subroutine SHCUPAR
  !
  !     ******************************************************************
  !
  subroutine SHCU_ENV(NZ,make_conv_grid_interpol)  

    integer, intent(IN) :: NZ
    logical, intent(IN) :: make_conv_grid_interpol
    real :: DLAMB(NKP)  ! , HZ(NKP)  ! not used
    !       Basic constants
    real :: CONST1, CONST2, ES00, EPSLON, UMMEPS, TA0, CONST3,  &
         C0, DLAMB0, ZREF, ZNZ, TLLL, PLLL, ZLLL, RLLL, DZLLL, DZDD, DTHV
    integer :: NKMID, K


    !  R=287.
    !  CP=1004.
    !  RCP=.286
    !  CPR=3.4965
    !  ALVL=2.5E6
    !  ALIV=2.834E6
    !  P00=1E5
    !  G=9.8
    !  AKLV=2340.6
    !  AKIV=2825.7
    !  ------------All variables above have been defined in USE shcu_vars_const
    DZLOW=200.
    DZHIGH=500.
    !srf      ZMID=3000.
    ZMID=4000.
    !  GAMMAD=.00976
    CONST1=17223003.15
    CONST2=29.65
    ES00=611.2
    EPSLON=.622
    UMMEPS=.378
    TA0=273.15
    CONST3=17.67
    !
    CDZMIN=3000.
    !EPS
    C0=0.0
    !_modified in 01/22/99      DLAMB=0.003
    DLAMB0=1.E-06
    ZREF=700.
    !
    !
    !-srf avoiding interpolation to the finer/convective grid, old RAMS approach
    if(make_conv_grid_interpol) then
       ZC(1)=0.
       do K=2,nz
          ZC(K)=ZCON(K)
       enddo
       ZE(1)=0.
       do K=2,nz
          ZE(K)=ZZCON(K)
       enddo
       !FIND MODEL TOP ON NORMAL GRID
       KMT=NZ 

    else
       !           INTERPOLATE MODEL SOUNDING (ENVIRONMENT) TO HIGHER
       !             RESOLUTION GRID

       !- original version
       NKMID=ZMID/DZLOW+1
       ZC(1)=0.
       do K=2,NKMID
          ZC(K)=ZC(K-1)+DZLOW
       enddo
       do K=NKMID+1,NKP
          ZC(K)=ZC(K-1)+DZHIGH
       enddo
       ZE(1)=0.
       do K=2,NKP
          ZE(K)=(ZC(K)+ZC(K-1))*.5
       enddo
       !                   FIND MODEL TOP ON CONVECTIVE GRID
       ZNZ=ZCON(NZ)
       do K=NKP,1,-1
          if(ZE(K).lt.ZNZ)GO TO 13
       enddo
       stop ' ENVIR STOP 12'
13     continue
       KMT=K

    endif
    !-srf 

    !-   DO ACTUAL INTERPOLATION to convective or normal grid
    call HTINT(NZ,WCON,ZZCON,KMT,WPE,ZE)
    call HTINT(NZ,THTCON,ZCON,KMT,THE,ZE)
    call HTINT(NZ,QVCON,ZCON,KMT,QVE,ZE)
    call HTINT(NZ,AKVD,ZCON,KMT,AKVDE,ZE)
    call HTINT(NZ,CL_CON,ZCON,KMT,CL_PE,ZE)


    !
    do K=1,KMT
       QVE(K)=max(QVE(K),1E-8)
    enddo
    !
    !         COMPUTE THETA V, THETA E, AND GET PRESSURE PROFILE
    !
    PKE(1)=PICON(1)
    do K=1,KMT
       THVE(K)=THE(K)*(1.+.61*QVE(K))
    enddo
    !
    do K=2,KMT
       PKE(K)=PKE(K-1)-G*2.*(ZE(K)-ZE(K-1))/(THVE(K)+THVE(K-1))
    enddo
    do K=1,KMT
       TE(K)=THE(K)*PKE(K)/CP
       PE(K)=(PKE(K)/CP)**CPR*P00
       RHOE(K)=PE(K)/(R*TE(K)*(1.+.61*QVE(K)))
    enddo
    !
    !EPS    FIND THE MAIN SOURCE LEVEL OF THE UPDRAFT. WE WILL ASSUME
    !EPSTHAT PARCELS START FROM THE FIRST LAYER DEFINED BY K=2
    !
    KCON=2
    !
    !         FIND THE LCL OF A LAYER AVERAGE AROUND THE SOURCE LEVEL
    !
    !77 CONTINUE !Never used
    !      TLLL=(TE(KCON)+TE(KCON+1)+TE(KCON-1))/3.
    TLLL=(TE(KCON)+TE(KCON-1))/2.      
    PLLL=PE(KCON)
    !      RLLL=(QVE(KCON)+QVE(KCON+1)+QVE(KCON-1))/3.
    RLLL=(QVE(KCON)+QVE(KCON-1))/2.      
    ZLLL=ZE(KCON)
    !
    call LCL(TLLL,PLLL,RLLL,TLCL,PLCL,DZLCL)
    !
    !         FIND THE CLOSEST LEVEL ON THE CONVECTIVE GRID TO THE LCL
    !
    DZLLL=1E20
    do K=1,KMT
       DZDD=abs(ZE(K)-(ZLLL+DZLCL))
       if(DZDD.lt.DZLLL)then
          DZLLL=DZDD
          KLCL=K
       endif
    enddo
    !
    !     Determination of Zi, the PBL top as the level where
    !     the turbulent diffusion coefficient (AKVDE) vanishes. The loop
    !     starts at k=3 because sometimes AKVDE is zero at k=1 or 2.
    !
    !      AKVMIN=1.
    !      DO K=3,KMT
    !        IF(AKVDE(K).LT.AKVMIN.or.abs(cl_pe(k)).gt.1.e-5) THEN
    !          KZI=K-1
    !          GO TO 78
    !        ENDIF
    !      ENDDO
    !
    !      Determination of Zi, the PBL top as the level where
    !      the virtual potential temperature THVE increases by
    !      more than DTHV as compared to the previous level
    !
    DTHV=0.5
    do K=3,KMT
       if((THVE(K)-(THVE(K-1))).gt.DTHV) then
          !eps      KZI=K-1
          KZI=K
          exit
       endif
    enddo
    !
    !78 CONTINUE !Not used, commented line 509
    !
    !     Following Wilde et. al. (1985) shallow cumulus onset occurs
    !     when the top of the entrainment zone overlaps the base of 
    !     the LCL zone. 
    !
    !_Comentado em 10/03/99 para teste  de sensibilidade
    if(KZI.lt.KLCL) then
       IGO=0
       return
    endif
    !
    !EPS
    !determination of environment variables
    !
    do K=1,KMT
       DSE(K)=CP*TE(K)+G*ZE(K)
       UHE(K)=DSE(K)+ALVL*QVE(K)
       EVAPS(K)=ES00*exp(CONST3*(TE(K)-TA0)/(TE(K)-CONST2))
       QVSE(K)=EPSLON*EVAPS(K)/(PE(K)-UMMEPS*EVAPS(K))
       UHES(K)=DSE(K)+ALVL*QVSE(K)
       RHE(K)=QVE(K)/QVSE(K)
       !
       gamma(K)=CONST1*PE(K)*QVSE(K)**2/EVAPS(K)
       gamma(K)=gamma(K)/((TE(K)-CONST2)*(TE(K)-CONST2))
    enddo
    !
    !     calculating the cloud moist static energy profile. We assume
    !     that entrainment occurs only above the LCL.
    !      
    UHC(1)=UHE(1)
    UHC(2)=UHE(2)
    DELZ(2)=ZE(2)-ZE(1)
    !
    !     Option 1: the parcel is entrained since its source level
    !
    !C      DO K=3,KMT
    !C          DLAMB(K)=EXP(LOG(DLAMB0)+2.3*(ZE(K)-ZE(1))/ZREF)            
    !C          DLDZBY2(K)=MIN(1.,DLAMB(K))
    !C          DLDZBY2(K)=DLDZBY2(K)*DELZ(K)/2
    !     
    !C        UHC(K)=(UHC(K-1)-DLDZBY2(K)*(UHC(K-1)-UHE(K)-UHE(K-1)))/  &
    !C       (1+DLDZBY2(K))
    !C      ENDDO
    !CCC
    !
    !     Option 2: the parcel rises without mixing till its lcl and 
    !     starts being entrained from the lcl up to the cloud top.
    !
    if(KLCL.ge.3) then
       do K=3,KLCL
          !      IF(KLCL.GE.2) THEN
          !         DO K=2,KLCL         
          UHC(K)=UHE(2)
          DELZ(K)=ZE(K)-ZE(K-1)
       enddo
       !
       do K=KLCL+1,KMT
          DELZ(K)=ZE(K)-ZE(K-1)
          DLAMB(K)=exp(log(DLAMB0)+2.3*(ZE(K)-ZE(1))/ZREF)            
          DLDZBY2(K)=min(1.,DLAMB(K))
          DLDZBY2(K)=DLDZBY2(K)*DELZ(K)/2
          !     
          UHC(K)=(UHC(K-1)-DLDZBY2(K)*(UHC(K-1)-UHE(K)-UHE(K-1)))   &
               /(1+DLDZBY2(K))
       enddo
    else
       IGO=0
       !        WRITE(11,*)'C',KLCL
       return
    endif
    !
    !calculating the in-cloud variables
    !
    do K=1,KMT
       DSC(K) =DSE(K)+(UHC(K)-UHES(K))/(1+gamma(K))
       DSC0(K)=DSE(K)+(UHE(2)-UHES(K))/(1+gamma(K))
       QVC(K) =QVSE(K)+gamma(K)*(UHC(K)-UHES(K))/(ALVL*(1+gamma(K)))
       WLC(K)=0.0
    enddo
    !
    do K=KLCL+1,KMT    
       WLC(K)=WLC(K-1)-(QVC(K)-QVC(K-1))-DLAMB(K)*              &
            (QVC(K)-QVE(K))*DELZ(K)+(C0-DLAMB(K))*WLC(K-1)*DELZ(K)
       ! FOR DEBUG:
       !* 250 Floating-point data overflow PROG=shcu_env ELN=599(4003dc44c)
       !* 253 Invalid operation PROG=shcu_env ELN=599(4003dc44c)
       !WLC(K)=WLC(K-1)
       !WLC(K)=WLC(K) - (QVC(K)-QVC(K-1))
       !TEMPOR1 = DLAMB(K)*(QVC(K)-QVE(K))
       !TEMPOR1 = TEMPOR1*DELZ(K)
       !WLC(K)=WLC(K) - TEMPOR1
       !WLC(K)=WLC(K) + (C0-DLAMB(K))*WLC(K-1)*DELZ(K)

       !IF (WLC(K)<0) EXIT

       WLC(K)=min(QVC(KLCL),WLC(K))
       WLC(K)=max(.1E-12,WLC(K))

    enddo
    !
    return
  end subroutine SHCU_ENV
  !
  !     ******************************************************************
  !
  subroutine SH2MOD(m1)

    integer, intent(IN) :: m1

    ! Local variables
    integer :: K
    real :: TFTC, TFRC, TFTM, TFRM, FTRES, FRRES

!!$  ! External Function to Sum a array
!!$  real, external :: ssum

    !        Compute integrated heating and moistening tendencies
    !
    do K=2,KMT
       QVCT1(K) = RHOE(K)*DTDT(K)*PKE(K)
       QVCT2(K) = RHOE(K)*ALVL*DRDT(K)
       QVCT3(K) = (ZC(K)-ZC(K-1))*QVCT1(K)
       QVCT4(K) = (ZC(K)-ZC(K-1))*QVCT2(K)
       !        print*,'0',QVCT1(K),QVCT2(k),k,DRDT(K),zc(k)
    enddo
!!$  TFTC=SSUM(KMT-1,QVCT3(2),1)
!!$  TFRC=SSUM(KMT-1,QVCT4(2),1)
    TFTC = sum(QVCT3(2:KMT))
    TFRC = sum(QVCT4(2:KMT))
    !
    !         Transfer tendencies to model grid
    !
    !new--------
    call vertmap2(qvct1,zc,kmt,vctr5,zzcon,m1-1)
    call vertmap2(qvct2,zc,kmt,vctr6,zzcon,m1-1)

    do K=2,m1-1
       VCTR5(K) = VCTR5(K)*(ZZCON(K)-ZZCON(K-1))
       VCTR6(K) = VCTR6(K)*(ZZCON(K)-ZZCON(K-1))
    enddo
    !new--------

    !old--------
    !      DO K=1,m1
    !        VCTR5(K)=0.
    !        VCTR6(K)=0.
    !      ENDDO
    !
    !      DZLFT=0.
    !      L=2
    !      DO K=2,m1-1
    !        IF(DZLFT.NE.0.) THEN
    !          VCTR5(K)=VCTR5(K)+QVCT1(L)*DZLFT
    !          VCTR6(K)=VCTR6(K)+QVCT2(L)*DZLFT
    !          L=L+1
    !        ENDIF
    !   60   CONTINUE
    !        IF(ZC(L).LE.ZZCON(K)) THEN
    !          VCTR5(K)=VCTR5(K)+QVCT1(L)*(ZC(L)-ZC(L-1))
    !          VCTR6(K)=VCTR6(K)+QVCT2(L)*(ZC(L)-ZC(L-1))
    !          L=L+1
    !          DZLFT=0.
    !          GO TO 60
    !        ELSE
    !          VCTR5(K)=VCTR5(K)+QVCT1(L)*(ZZCON(K)-ZC(L-1))
    !          VCTR6(K)=VCTR6(K)+QVCT2(L)*(ZZCON(K)-ZC(L-1))
    !          DZLFT=ZC(L)-ZZCON(K)
    !        ENDIF
    !      ENDDO
    !
    !old--------
    !
    !         Make sure the transfer from the convective grid to the model
    !           grid happened correctly.
    !
!!$  TFTM=SSUM(m1-2,VCTR5(2),1)
!!$  TFRM=SSUM(m1-2,VCTR6(2),1)
    TFTM = sum(VCTR5(2:m1))
    TFRM = sum(VCTR6(2:m1))
    !
    FTRES=TFTM-TFTC
    FRRES=TFRM-TFRC
    if(abs(FTRES).gt.(0.01*abs(TFTC))) then
       print*,' Energy error in grid tranfser in convective param.'
       print*,' TFTM,TFTC ',TFTM,TFTC,100.*ftres/TFTC
    endif
    !
    !         Change energy tendencies to temperature and mixing ratio
    !           tendencies.
    !
    do K=2,m1-1
       DTDT(K)=VCTR5(K)/((ZZCON(K)-ZZCON(K-1))*DNCON(K)*PICON(K))
       DRDT(K)=VCTR6(K)/((ZZCON(K)-ZZCON(K-1))*DNCON(K)*ALVL)
       !print*,'1',k,DTDT(k)*86400.,DRDT(k)*86400.,zzcon(k)
    enddo

    !
    return
  end subroutine SH2MOD

  !
  !     ******************************************************************
  !      
  subroutine CL_TOP

    !
    !EPS      
    !
    !

    ! Local variables
    real :: DSCVM(NKP), DSEVM(NKP), TEM(NKP), EMP(NKP)
    real :: GAMMAD, BUOY1, BUOY2
    integer :: K, L

    !       Basic constants

    !  R=287.
    !  CP=1004.5
    !  RCP=.286
    !  CPR=3.4965
    !  ALVL=2.50E6
    !  ALIV=2.837E6
    !  P00=1E5
    !  G=9.806
    !  AKLV=2340.6
    !  AKIV=2825.7
    !  All variables above have been defined in USE shcu_vars_const

    GAMMAD=.00976
    !EPS
    !     Determination of the cloud top based on integrated cloud buoyancy. 
    !     First, determination of virtual static energy profiles 
    !     Sv=S+0.608Cp<T>q-Cp<T>l
    !
    do K=1,KMT
       DSCV(K)=DSC(K)+CP*TE(K)*(.608*QVC(K)-WLC(K))
       DSEV(K)=DSE(K)+.608*CP*TE(K)*QVE(K)
       DSC0V(K)=DSC0(K)+.608*CP*TE(K)*QVE(K)
       !       DSC0V(K)=DSC0(K)+CP*TE(K)*(.608*QVC(K)-WLC(K))
    enddo
    !
    do K=2,KMT
       DSCVM(K)=(DSCV(K)+DSCV(K-1))/2
       DSEVM(K)=(DSEV(K)+DSEV(K-1))/2
       DSC0VM(K)=(DSC0V(K)+DSC0V(K-1))/2
       TEM(K)=(TE(K)+TE(K-1))/2
    enddo
    !
    !     determination of the integrated bouyancy between the surface and 
    !      the LCL (see Albrecht et al. 1986, Eq. A5)
    !     ENTF is the enthalpy flux at surface in (K*m/s)
    !
    if(ENTF.gt.0.0) then
       BUOY1=ENTF*(1+.608*QVE(1))
       BUOY1=G*ZE(KLCL)*BUOY1/TE(1)
       BUOY1=1.3333*BUOY1**.6667
    else
       BUOY1=.0
       IGO=0
       !        WRITE(11,*)'D'
       return
    endif
    !
    !     checking if the parcel is able to sustain positive buoyancy one
    !     level above the LCL
    !
    CAPE=.0
    BUOY2=.0
    L=KLCL+1
    BUOY2=BUOY2+GAMMAD*(DSCVM(L)-DSEVM(L))*DELZ(L)/TEM(L)
    !                           
    if((BUOY1+BUOY2).le.0.0) then
       IGO=0 
       !       WRITE(11,*)'BUOY'      
       return
    endif
    !           
    !     Integration continues till the level of zero buoyancy is found.
    !     The cloud top is assumed to be one level below.
    !
    !      
    do K=KLCL+2,KMT                 
       BUOY2=BUOY2+GAMMAD*(DSCVM(K)-DSEVM(K))*DELZ(K)/TEM(K)
       if((BUOY1+BUOY2).le.0.0) then
          KTOP=K-1
          GO TO 88
       endif
    enddo
88  continue      
    !
    do K=KLCL,KTOP
       EMP(K)=DSC0VM(K)-DSEVM(K)
       EMP(K)=max(0.,EMP(K))
       CAPE=CAPE+GAMMAD*EMP(K)*DELZ(K)/TEM(K)
    enddo
    !
    return
  end subroutine CL_TOP
  !
  !     ******************************************************************
  ! 
  subroutine W_SHALLOW(IP,JP,TIME)

    integer, intent(IN) :: IP, JP
    real, intent(IN)    :: TIME

    !      
    !     The vertical velocity at cloud base is calculated according
    !     to the heat engine framework as defined by Renno'  and
    !     Ingersoll, 1996 Eq(42)
    !

    ! Local variables
    real :: CMI, TCOLD, THOT, FIN, SIGWSHB
    integer :: ICOUNT, K

    !  CP=1004.   ! defined in USE shcu_vars_const
    CMI=16.
    !  EMIS=.9         ! not used
    !  STEFAN=5.67E-8  ! not used
    !
    TCOLD=TE(KLCL)
    !*      TCOLD=TE(2)
    ICOUNT=1
    do K=KLCL+1,KTOP
       !*      DO K=3,KTOP
       TCOLD=TCOLD+TE(K)
       ICOUNT=ICOUNT+1
    enddo
    !
    !     TCOLD is the temperature of the environment, averaged
    !     between the first level and the cloud top.
    !      
    TCOLD=TCOLD/float(ICOUNT)
    THOT=TE(2)
    !
    EFIC=(THOT-TCOLD)/THOT
    !
    if(EFIC.le.0.0.or.CAPE.lt.20.0)then
       IGO=0 
       !        WRITE(11,*)'Eta= ',EFIC,'CAP= ',CAPE       
       return
    endif
    !
    !     the effective vertical velocity at cloud base is calculated 
    !     according to the heat engine framework as deffined by Renno'
    !     and Ingersoll, 1996 Eq.(34)
    !
    FIN=RHOE(2)*(CP*ENTF+ALVL*ALHF)
    !
    if(FIN.le.50.0) then
       IGO=0
       !        WRITE(11,*)'F'
       return
    endif
    !      
    !     TCAPE=2*CAPE
    !
    !     Calculating BCAPE
    !
    !      WSTAR=ENTF*(1+0.608*QVE(1))
    !      WSTAR=(G*ZE(KLCL)*WSTAR/TE(1))**.333
    !      BCAPE=CMI*WSTAR**2
    !
    !     Calculating DCAPE. Here we simply assume DCAPE=0.5CAPE
    !
    !eps DCAPE=1.0*CAPE
    DCAPE=0.
    !
    !      TCAPE=CAPE+DCAPE+BCAPE
    TCAPE=CAPE+DCAPE     
    !
    if(TCAPE.lt.50.0) then
       IGO=0
       return
    endif
    !
    SIGWSHB=EFIC*FIN/(RHOE(KLCL)*TCAPE)
    !      SIGWSHB=SIGWSHB+WPE(KLCL)
    !
    if(SIGWSHB.le.0.0) then
       IGO=0
       !WRITE(11,*)'G'
       return
    endif
    !
    !
    !     calculating the effetive vertical velocity inside the cloud
    !     interpolating from SIGWSHB in the cloud base to zero
    !     on the cloud top.
    ! 
    !   
    do K=KLCL,KTOP
       WC(K)=SIGWSHB*((ZE(KTOP)-ZE(K))/(ZE(KTOP)-ZE(KLCL)))
    enddo
    !
    return
  end subroutine W_SHALLOW
  !
  !     ******************************************************************
  !
  subroutine SH_RATES

    ! Local variables                                                
    real :: WSSC(NKP), WQSC(NKP), DSDT(NKP)
    real :: WRSUP, WRBASE
    integer :: K, ISUBCL

    do K=1,KMT
       WSSC(K)=0.
       WQSC(K)=0.
       DSDT(K)=0.
       DRDT(K)=0.
       DTDT(K)=0.
    enddo
    !
    !     calculating the transports w's' and w'r'
    !      
    do K=KLCL+1,KTOP
       WSSC(K)=WC(K)*(DSC(K)-ALVL*WLC(K)-DSE(K))
       WQSC(K)=WC(K)*(QVC(K)+WLC(K)-QVE(K))
    enddo
    !     calculating heating and moistening rates due to
    !     shallow non-precipitating cumulus
    !
    do K=KLCL+1,KTOP-1
       DQDT(K)=-(WQSC(K+1)-WQSC(K-1))/(ZE(K+1)-ZE(K-1))
       DRDT(K)=DQDT(K)/(1-QVE(K))**2
       DSDT(K)=-(WSSC(K+1)-WSSC(K-1))/(ZE(K+1)-ZE(K-1))
       DTDT(K)=DSDT(K)/PKE(K)
    enddo
    !
    !     ISUBCL is a flag corresponding to removal (=1) or not (=0)
    !     of moisture from the boudary layer by shallow clouds
    !
    ISUBCL=0
    !                                                           
    if(ISUBCL.eq.0)return
    !
    !     for calculating the removal of moisture from the PBL, we 
    !     consider that the latent heat flux is linearly interpolated between
    !     the surface and the cloud base
    !                            
    WRSUP=ALHF
    WRBASE=WC(KLCL)*(QVC(KLCL)-QVE(KLCL))
    !
    do K=1,KLCL
       WC(K)  =(WRBASE*ZE(K)+(WRSUP*(ZE(KLCL)-ZE(K))))/(QVE(K)*ZE(KLCL))
       WQSC(K)= WC(K)*(QVC(K)-QVE(K))
    enddo
    !
    do K=2,KLCL
       DRDT(K)=-(WQSC(K)-WQSC(K-1))/DELZ(K)
    enddo
    !

    !  PRINT *, "Saindo de SH_RATES"

    return
  end subroutine SH_RATES
  !                           
  !*********************************************************************
end module ModRShCuPar
