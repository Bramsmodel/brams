!############################# Change Log ##################################
! 5.0.2
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003, 2006 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################
module ModLeaf3Teb

  use ModUrban, only: &
       urban

  use ModGasPart, only: &
       emfactor

  use teb_vars_const, only: &
       tminbld,             & ! intent(in)
       d_road,              & ! intent(in)
       tc_road,             & ! intent(in)
       d_wall,              & ! intent(in)
       tc_wall,             & ! intent(in)
       d_roof,              & ! intent(in)
       tc_roof,             & ! intent(in)
       nurbtype,            & ! intent(in)
       ileafcod, &               ! intent(in)
       daylight,    &
       rushh1,      &
       rushh2,      &
       xcpd,        &
       xrd,         &
       xrv,         &
       z0_town,     &
       bld,         &
       bld_height,  &
       bld_hl_ratio,&
       aroof,       &
       eroof,       &
       aroad,       &
       eroad,       &
       awall,       &
       htraf,       &
       hindu,       &
       pletraf,     &
       daylight,    &
       ewall,       &
       hc_roof,     &
       hc_road,     &
       hc_wall,     &
       pleindu

  use mem_emiss, only: &
       efsat, &
       efsun, &
       weekdayin

  implicit none
  
  private

  public :: tebc_init
  public :: teb_init
  public :: leaf3_teb_interface

contains

  !Subroutine to link leaf3 to TEB for landclass = urban type 1 and 
  !urban type 2
  !Adapted by Edmilson Freitas 
  !DATE : Jul 7th 2006
  !Last modification:

  ! TEB
  subroutine leaf3_teb_interface(             &
       istp, ztstepfrc, ztstep, cosz, zzref,  &
       zrat, solar,                           &
       zpa, zta, zua, zva, rv,                &
       zpcp_in, fuso,                         &
       zt_canyon, zr_canyon,                  &
       zts_roof, zts_road, zts_wall,          &
       zti_road, zti_bld,                     &
       zws_roof, zws_road,                    &
       zt_roof, zt_road, zt_wall,             &
       zh_town, zle_town, zemis_town,         &
       zsfu_town, zsfv_town, zts_town,        &
       zalb_town, ig_urban,                   &
       zh_traffic, zh_industry,               &
       zle_traffic, zle_industry,             &
       t2m, r2m,                              &
       time, itime1, dpdz, dens               )
    ! constants:
    integer, parameter  :: inteb = 3    ! teb vertical dimension
    ! arguments:
    integer, intent(in) :: istp ! not used
    real, intent(in)    :: ztstepfrc !time step of the input forcing ! not used
    real, intent(in)    :: ZTSTEP    ! Time step of integration (sec.)
    real, intent(in)    :: COSZ      ! Cos(Zenith Angle)
    real, intent(in)    :: ZZREF     ! reference height of the first atm. level
    real, intent(in)    :: ZRAT      ! atmospheric infrared radiation
    real, intent(in)    :: SOLAR
    real, intent(in)    :: ZPA       ! atmospheric level pressure
    real, intent(in)    :: ZTA       ! atmospheric temperature at level za
    real, intent(in)    :: ZUA, ZVA  ! wind speeds in the x and y dirs at level za
    real, intent(in)    :: RV        ! atmospheric mixing ratio at level za
    real, intent(in)    :: ZPCP_IN   ! Precip. Rate (kg/m^2/s)
    real, intent(in)    :: fuso      ! local time correction
    real, intent(inout) :: ZT_CANYON ! canyon air temperature
    real, intent(inout) :: ZR_CANYON ! canyon air vapor mixing ratio
    real, intent(inout) :: ZTS_ROOF  ! roof surface temperature
    real, intent(inout) :: ZTS_ROAD  ! road surface temperature
    real, intent(inout) :: ZTS_WALL  ! wall surface temperature
    real, intent(inout) :: ZTI_ROAD  ! road deep temperature
    real, intent(inout) :: ZTI_BLD   ! INTERNAL BLDING TEMP (K)
    real, intent(inout) :: ZWS_ROOF  ! roof water reservoir
    real, intent(inout) :: ZWS_ROAD  ! road water reservoir
    real, intent(inout) :: ZT_ROOF(INTEB) ! roof layers temperatures
    real, intent(inout) :: ZT_ROAD(INTEB) ! road layers temperatures
    real, intent(inout) :: ZT_WALL(INTEB) ! wall layers temperatures
    real, intent(out)   :: ZH_TOWN        ! TOWN AVE. SENS. HEAT FLUX
    real, intent(out)   :: ZLE_TOWN       !TOWN AVE. LAT. HEAT FLUX
    real, intent(out)   :: ZEMIS_TOWN     ! town equivalent emissivity
    real, intent(out)   :: ZSFU_TOWN      !TOWN SCALE EDDY U-MOM FLUX
    real, intent(out)   :: ZSFV_TOWN      !TOWN SCALE EDDY V-MOM FLUX
    real, intent(out)   :: ZTS_TOWN       !TOWN SFC TEMP
    real, intent(out)   :: ZALB_TOWN      !TOWN EQV. ALBEDO
    integer, intent(in) :: IG_URBAN    ! Flag for urban(1) and suburban(2)
    real, intent(out)   :: ZH_TRAFFIC
    real, intent(out)   :: ZH_INDUSTRY
    real, intent(out)   :: ZLE_TRAFFIC
    real, intent(out)   :: ZLE_INDUSTRY
    real, intent(out)   :: T2M            ! Extrapolated 2 m temperature
    real, intent(out)   :: R2M            ! Extrapolated 2 m specif humidity
    integer, intent(in) :: itime1
    real, intent(in)    :: time
    real, intent(in)    :: dpdz
    real, intent(in)    :: dens
    ! Local Variables:
    ! Declarations of variables
    !INPUT / OUTPUT VARIABLES
    real :: pfat

    real :: ZRG        ! incoming solar radiation
    real :: ZPS        ! surface pressure
    real :: ZQA        ! atmospheric specific humidity at level za
    real :: ZRR        ! rain rate
!!$  real :: ZPRECIP    ! precipitation rate ! not used
    !    0.9.D   TEB Diagnostics:
    !           ----------------
    real :: ZRN_ROOF     ! net radiation over roof
    real :: ZH_ROOF      ! sensible heat flux over roof
    real :: ZLE_ROOF     ! latent heat flux over roof
    real :: ZGFLUX_ROOF  ! flux through the roof
    real :: ZRUNOFF_ROOF ! runoff over the ground
    real :: ZRN_ROAD     ! net radiation over road
    real :: ZH_ROAD      ! sensible heat flux over road
    real :: ZLE_ROAD     ! latent heat flux over road
    real :: ZGFLUX_ROAD  ! flux through the road
    real :: ZRUNOFF_ROAD ! runoff over the ground
    real :: ZRN_WALL     ! net radiation over wall
    real :: ZH_WALL      ! sensible heat flux over wall
    real :: ZLE_WALL     ! latent heat flux over wall
    real :: ZGFLUX_WALL  ! flux through the wall
    !
    !only OUTPUT VARIABLES
    real :: ZU_CANYON    ! CANYON HOR. WIND
    real :: ZRN_TOWN     ! TOWN SCALE NET RAD
    real :: ZCH_TOWN     ! TOWN AVERAGED HEAT TRANSFER
    real :: ZGFLUX_TOWN  ! TOWN SCALE GROUND HEAT STORAGE
    real :: ZRUNOFF_TOWN ! TOWN SCALE RUNOFF
    !  MODEL FORCING VARIABLES
    real :: ZSVF_ROAD     ! ROAD SKY VIEW FACTOR
    real :: ZSVF_WALL     ! WALL SKY VIEW FACTOR
    real :: ZCAN_HW_RATIO ! CANYON HEIGHT TO WIDTH RATIO
    real :: ZDIR_SW_RAD   ! incoming direct solar radiation on an horiz.surface
    real :: ZSCA_SW_RAD   ! scattered incoming solar rad.

    integer :: I  !!, JDT, JINDT
    real ::       ZZ0_TOWN    &
         ,ZBLD    &
         ,ZBLD_HEIGHT    &
         ,ZBLD_HL_RATIO    &
         ,ZALB_ROOF    &
         ,ZEMIS_ROOF    &
         ,ZALB_ROAD    &
         ,ZEMIS_ROAD    &
         ,ZALB_WALL    &
         ,ZEMIS_WALL    &
         ,        ax1, ax2, tmp_teb, bx1, bx2 &
         ,        timeq1, timeq2, tign
    integer :: idays
    character(len=3) :: cday
    real :: ZDIRCOSZW, ZTANZEN
    !                  ZDIRCOSZW = Cosinus of the angle between the 
    !                  normal to the surface and the vertical
    !                  ZTANZEN   = tangent of solar zenith angle
    !
    real   :: ZRHOA,  ZEXNA, ZEXNS, ZVMOD, ZTVI, ZAZENIT, ZEXN2, ZP2
    !                  ZRHOA   = air density
    !                  ZEXNA   = Exner function
    !                  ZEXNS   = Exner function
    !                  ZVMOD   = modulus of the wind parallel to the orography
    !                  ZTVI    = virtual temperature 
    !                  ZAZENIT = solar zenith angle
    real, external ::  rslf ! Function in "therm_lib.f90"


    ZRG = SOLAR

    ZDIRCOSZW = 1.0

    ZQA = (1./((1000./rv) + 1.))*1000.


    ZZ0_TOWN        =  Z0_TOWN(IG_URBAN)    
    ZBLD            =  BLD (IG_URBAN) 
    ZBLD_HEIGHT     =  BLD_HEIGHT (IG_URBAN)    
    ZBLD_HL_RATIO   =  BLD_HL_RATIO(IG_URBAN)
    ZALB_ROOF       =  AROOF(IG_URBAN)
    ZEMIS_ROOF      =  EROOF(IG_URBAN)
    ZALB_ROAD       =  AROAD(IG_URBAN)
    ZEMIS_ROAD      =  EROAD(IG_URBAN)
    ZALB_WALL       =  AWALL(IG_URBAN)
    ZEMIS_WALL      =  EWALL(IG_URBAN)

    ZH_TRAFFIC      =  HTRAF(IG_URBAN)
    ZH_INDUSTRY     =  HINDU(IG_URBAN)
    ZLE_TRAFFIC     =  PLETRAF(IG_URBAN) 
    ZLE_INDUSTRY    =  PLEINDU(IG_URBAN)

    !
    !*      1.     Calculate the canyon shape
    !              --------------------------
    !
    ZCAN_HW_RATIO = ZBLD_HL_RATIO*ZBLD/(1.0 - ZBLD)
    !write(*,*)'ZCAN_HW_RATIO=',ZCAN_HW_RATIO
    !write(*,*)'ZBLD_HL_RATIO=',ZBLD_HL_RATIO
    !pause
    !!*      2.     Calculate sky view factors
    !              --------------------------
    !
    ZSVF_ROAD  = sqrt(ZCAN_HW_RATIO*ZCAN_HW_RATIO + 1.0) - ZCAN_HW_RATIO
    ZSVF_WALL  = 0.5*(1.0 - ZSVF_ROAD)/ZCAN_HW_RATIO
    !*      3.     Calculate town albedo and emissivity
    !              ------------------------------------
    !
    ZALB_TOWN  = (1.0 - ZBLD)*ZALB_ROAD + ZBLD*ZALB_ROOF
    ZEMIS_TOWN = (1.0 - ZBLD)*ZEMIS_ROAD*ZSVF_ROAD + &
         (1.0 - ZBLD)*ZEMIS_WALL*(1.0 - ZSVF_ROAD) + ZBLD*ZEMIS_ROOF
    ZRR       = ZPCP_IN

    if (ZRG < 0.0) ZRG = 0.0

    ! Initialize forcing values:
    ZAZENIT = acos(COSZ)
    ZTANZEN = tan(ZAZENIT)
    if (ZRG>0.0 .and. ZTANZEN<0.0) then
       ZTANZEN = sqrt(9999.0)
    endif

    ! Perform some conversions and final calculations using forcing variables:
    !
    ZVMOD   = sqrt(ZUA*ZUA + ZVA*ZVA)
    !
    ! Make sure wind magnitude is above a minimum threshold (original =1.0, 
    ! for test purpose = 0.5) :
    !
    ZVMOD   = max(0.5, ZVMOD)
    ZTVI    = ZTA*(1. + ((XRV/XRD) - 1.)*ZQA )
    ZRHOA   = ZPA/XRD/ZTVI
    !            ZRHOA   = dens
    ZEXNA   = (ZPA/100000.)**(XRD/XCPD)

    ZPS     = DPDZ*(ZBLD_HEIGHT-ZZREF) + ZPA
    ZEXNS   = (ZPS/100000.)**(XRD/XCPD)

    ! Calculate forcing needed by TEB from existing forcing
    ! for now assume all incoming radiation is direct
    !
    ZDIR_SW_RAD = ZRG
    ZSCA_SW_RAD = 0.0
    !

    tmp_teb = time + (itime1/100 + mod(itime1,100)/60.)*3600
    idays   = int((tmp_teb/3600.)/24.)  !number of days of simulation
    tign    = real(idays)*24.*3600.

    pfat=1.

    call EMFACTOR(WEEKDAYIN, idays, cday)
    if (cday=='SAT') pfat = EFSAT
    if (cday=='SUN') pfat = EFSUN

    !Fonte urbana (castanho, 1999 - tese de mestrado)
    bx1    = RUSHH1 - fuso + DAYLIGHT
    bx2    = RUSHH2 - fuso + DAYLIGHT
    timeq1 = (tmp_teb - tign)/3600. - bx1
    timeq2 = (tmp_teb - tign)/3600. - bx2

    !fonte de calor sensivel veicular
    ax2        = ZH_TRAFFIC - 5.
    ax1        = 0.63*ax2
    ZH_TRAFFIC = ((ax1*exp(-(timeq1)**2/8.5) + &
         ax2*exp(-(timeq2)**2/10.6))*pfat + 5.)

    !fonte de calor latente veicular
    ax2         = ZLE_TRAFFIC - 5.
    ax1         = 0.63*ax2
    ZLE_TRAFFIC = ((ax1*exp(-(timeq1)**2/8.5) + &
         ax2*exp(-(timeq2)**2/10.6))*pfat + 5.)

    !                              -------------------------------
    ZRN_ROOF     = 0.
    ZH_ROOF      = 0.
    ZLE_ROOF     = 0.
    ZGFLUX_ROOF  = 0.
    ZRUNOFF_ROOF = 0.
    ZRN_ROAD     = 0.
    ZH_ROAD      = 0.
    ZLE_ROAD     = 0.
    ZGFLUX_ROAD  = 0.
    ZRUNOFF_ROAD = 0.
    ZRN_WALL     = 0.
    ZH_WALL      = 0
    ZLE_WALL     = 0.
    ZGFLUX_WALL  = 0.
    ZRN_TOWN     = 0.
    ZH_TOWN      = 0.
    ZLE_TOWN     = 0.
    ZGFLUX_TOWN  = 0.
    ZRUNOFF_TOWN = 0.
    ZSFU_TOWN    = 0.
    ZSFV_TOWN    = 0.
    ZCH_TOWN     = 0.

    !*                        14.B   CALL THE MAIN SUBROUTINE OF TEB
    call URBAN(ZTS_TOWN, ZEMIS_TOWN, ZALB_TOWN,                   &
         ZT_CANYON, ZR_CANYON, ZU_CANYON,                         &
         ZTS_ROOF,ZTS_ROAD,ZTS_WALL,ZTI_ROAD,ZTI_BLD,             &
         ZT_ROOF, ZT_ROAD, ZT_WALL, ZWS_ROOF,ZWS_ROAD,            &
         ZPS, ZPA, ZEXNS, ZEXNA, ZTA, ZQA, ZRHOA,                 &
         ZRAT, ZDIR_SW_RAD, ZSCA_SW_RAD, ZTANZEN,                 &
         ZRR, ZZREF, ZDIRCOSZW, ZUA, ZVA, ZVMOD,                  &
         ZH_TRAFFIC, ZLE_TRAFFIC, ZH_INDUSTRY, ZLE_INDUSTRY,      &
         ZTSTEP,                                                  &
         ZZ0_TOWN,                                                &
         ZBLD, ZBLD_HEIGHT, (2.*ZBLD_HL_RATIO*ZBLD),              &
         ZCAN_HW_RATIO, ZALB_ROOF, ZEMIS_ROOF,                    &
         HC_ROOF(1:3), TC_ROOF(1:3), D_ROOF(1:3),                 &
         ZALB_ROAD, ZEMIS_ROAD, ZSVF_ROAD,                        &
         HC_ROAD(1:3), TC_ROAD(1:3), D_ROAD(1:3),                 &
         ZALB_WALL, ZEMIS_WALL, ZSVF_WALL,                        &
         HC_WALL(1:3), TC_WALL(1:3), D_WALL(1:3),                 &
         ZRN_ROOF, ZH_ROOF, ZLE_ROOF, ZGFLUX_ROOF,                &
         ZRUNOFF_ROOF,                                            &
         ZRN_ROAD, ZH_ROAD, ZLE_ROAD, ZGFLUX_ROAD,                &
         ZRUNOFF_ROAD,                                            &
         ZRN_WALL, ZH_WALL, ZLE_WALL, ZGFLUX_WALL,                &
         ZRN_TOWN, ZH_TOWN, ZLE_TOWN, ZGFLUX_TOWN,                &
         ZRUNOFF_TOWN, ZSFU_TOWN, ZSFV_TOWN, ZCH_TOWN)

    !Extrapolating temperature and specific humidity to 2 m height

    ZP2   = DPDZ*(2. - ZBLD_HEIGHT) + ZPS

    ZEXN2 = (ZP2/100000.)**(XRD/XCPD)

    T2M   = ZT_WALL(1)*ZEXN2/ZEXNS

    R2M   = ZR_CANYON*(rslf(ZP2,T2M)/rslf(ZPS, ZT_WALL(1)))

  end subroutine LEAF3_TEB_INTERFACE

  ! TEB
  !     ############################################
  subroutine INI_TG_PROFILE(PTS, PTI, PTC, PD, PT)
    !   ############################################
    !
    !!****  *INI_TG_PROFILE*  
    !!
    !!    PURPOSE
    !!    -------
    !
    !     Computes the equilibrium of a temperature profile through a
    !     conductive material, given the two extreme surface temperatures.
    !         
    !!      
    !!    AUTHOR
    !!    ------
    !!
    !!V. Masson           * Meteo-France *
    !!
    !!    MODIFICATIONS
    !!    -------------
    !!      Original 02/11/98 
    !!               17/05/2002 Edmilson Freitas (IAG-USP). Elimination of the
    !!                          the declarations that are not needed. Change from
    !!                          Module to normal subroutine.
    !----------------------------------------------------------------------------
    !
    !*       0.     DECLARATIONS
    !               ------------
    !
    !
    !
    !*      0.1    declarations of arguments
    !
    !INPUT VARIABLES:
    real, intent(in)  :: PTS      ! surface temperature
    real, intent(in)  :: PTI      ! internal temperature
    real, intent(in)  :: PTC(3)   ! thermal conductivity for roof layers
    real, intent(in)  :: PD(3)    ! depth of roof layers
    !
    !OUTPUT VARIABLES
    real, intent(out) :: PT(3)    ! layers temperatures
    !
    !*      0.2    declarations of local variables
    !
    real, dimension(size(PT))   :: ZA ! lower diag.
    real, dimension(size(PT))   :: ZB ! main  diag.
    real, dimension(size(PT))   :: ZC ! upper diag.
    real, dimension(size(PT))   :: ZY ! r.h.s.
    real, dimension(size(PT))   :: ZX ! solution
    !
    real, dimension(0:size(PT)) :: ZMTC_O_D
    ! mean thermal conductivity over distance between 2 layers
    !
    integer                     :: ILAYER ! number of roof,road or wall layers
    integer                     :: JLAYER ! loop counter
    !----------------------------------------------------------------------------
    !
    !*      1.     Layer thermal properties
    !              ------------------------
    !
    ILAYER      = size(PT)
    ZA(:)       = 0.
    ZB(:)       = 0.
    ZC(:)       = 0.
    ZX(:)       = 0.
    ZY(:)       = 0.
    ZMTC_O_D(:) = 0.
    !
    ZMTC_O_D(0) = 2.*PTC(1)/PD(1)
    !
    do JLAYER=1,ILAYER-1
       ZMTC_O_D(JLAYER) = 2./( PD(JLAYER)/PTC(JLAYER) + &
            PD(JLAYER+1)/PTC(JLAYER+1) )
    end do
    !
    ZMTC_O_D(ILAYER) = 2.*PTC(ILAYER)/PD(ILAYER)
    !
    !----------------------------------------------------------------------------
    !
    !*      2.     Surface layer coefficients
    !              --++++++------------------
    !
    !
    ZA(1) =   0.

    ZB(1) =   ZMTC_O_D(0) + ZMTC_O_D(1)

    ZC(1) = - ZMTC_O_D(1)
    !
    ZY(1) =   ZMTC_O_D(0)*PTS
    !
    !
    !----------------------------------------------------------------------------
    !
    !*      3.     Other layers coefficients
    !              -------------------------
    !
    do JLAYER=2,ILAYER-1
       ZA(JLAYER) = - ZMTC_O_D(JLAYER - 1)

       ZB(JLAYER) =   ZMTC_O_D(JLAYER - 1) + ZMTC_O_D(JLAYER)

       ZC(JLAYER) = - ZMTC_O_D(JLAYER)

       ZY(JLAYER) =   0.
    end do
    !
    !----------------------------------------------------------------------------
    !
    !*      4.     Inside layer coefficients
    !              -------------------------
    !
    ZA(ILAYER) = - ZMTC_O_D(ILAYER - 1)

    ZB(ILAYER) =   ZMTC_O_D(ILAYER - 1) + ZMTC_O_D(ILAYER)

    ZC(ILAYER) =   0.

    ZY(ILAYER) =   ZMTC_O_D(ILAYER)*PTI
    !
    !
    !----------------------------------------------------------------------------
    !
    !*      5.     Tri-diagonal system resolution
    !              ------------------------------
    !
    call TRID(ZX, ZA, ZB, ZC, ZY, ILAYER)
    !
    PT(:) = ZX(:)
    !
    !----------------------------------------------------------------------------

  end subroutine INI_TG_PROFILE

  !-------------------------------------------------------------------------------
  ! TEB
  subroutine TEB_INIT(n1, n2, n3, np, vegt, theta, rv, pi, pp,       &
       TROOF, TROAD, TWALL, TIBLD, TIROAD, TCANYON, RCANYON, TSROOF, &
       TSROAD, TSWALL, HT, LET, HIN, LEIN, WSROOF, WSROAD, EMISTOWN, &
       ALBTOWN, TSTOWN, G_URBAN)

    ! Arguments:
    integer, intent(IN) :: n1, n2, n3, np
    real, intent(IN)    :: vegt(n2,n3,np)
    real, intent(IN)    :: theta(n1,n2,n3)
    real, intent(IN)    :: rv(n1,n2,n3)
    real, intent(IN)    :: pi(n1,n2,n3)
    real, intent(IN)    :: pp(n1,n2,n3)
    real, intent(OUT)   :: TROOF(n1,n2,n3)
    real, intent(OUT)   :: TROAD(n1,n2,n3)
    real, intent(OUT)   :: TWALL(n1,n2,n3)
    real, intent(OUT)   :: TIBLD(n2,n3)
    real, intent(OUT)   :: TIROAD(n2,n3)
    real, intent(OUT)   :: TCANYON(n2,n3)
    real, intent(OUT)   :: RCANYON(n2,n3)
    real, intent(OUT)   :: TSROOF(n2,n3)
    real, intent(OUT)   :: TSROAD(n2,n3)
    real, intent(OUT)   :: TSWALL(n2,n3)
    real, intent(OUT)   :: HT(n2,n3)
    real, intent(OUT)   :: LET(n2,n3)
    real, intent(OUT)   :: HIN(n2,n3)
    real, intent(OUT)   :: LEIN(n2,n3)
    real, intent(OUT)   :: WSROOF(n2,n3)
    real, intent(OUT)   :: WSROAD(n2,n3)
    real, intent(INOUT) :: EMISTOWN(n2,n3)
    real, intent(INOUT) :: ALBTOWN(n2,n3)
    real, intent(INOUT) :: TSTOWN(n2,n3)
    real, intent(OUT)   :: G_URBAN(n2,n3,np)

    ! Local Variables:
    integer :: i, j, k, ilf, inp 
    real    :: cpi, hcpi, pis, pl2

    ! Initiating data
    HT      = 0.
    LET     = 0.
    HIN     = 0.
    LEIN    = 0.
    WSROOF  = 0.
    WSROAD  = 0.
    TIBLD   = 0.
    TIROAD  = 0.
    TCANYON = 0.
    RCANYON = 0.
    TSROOF  = 0.
    TSROAD  = 0.
    TSWALL  = 0.
    TROOF   = 0.
    TROAD   = 0.
    TWALL   = 0.
    G_URBAN = 0.

    cpi = 1./1004.
    hcpi = 0.5*cpi

    do i=1,n2

       do j=1,n3

          pis = (pp(1,i,j) + pi(1,i,j) + pp(2,i,j) + pi(2,i,j)) * hcpi

          pl2 = (pp(2,i,j) + pi(2,i,j)) * cpi

          do inp=2,np !patch looping

             do ilf=1,NURBTYPE
                if (nint(vegt(i,j,inp))==ILEAFCOD(ILF)) &
                     G_URBAN(i,j,inp)=float(ilf)
             enddo

             if (nint(G_URBAN(i,j,inp))/=0.) then

                !internal temperature defined in RAMSIN
                TIBLD(i,j)    = 273.16 + TMINBLD
                !surface
                !TIROAD(i,j)   = (theta(2,i,j)+theta(1,i,j))*0.5*pis
                TIROAD(i,j)   = 281.16
                !model's first level
                TSROOF(i,j)   = (theta(2,i,j))*pl2
                !surface
                TSROAD(i,j)   = (theta(2,i,j)+theta(1,i,j))*0.5*pis
                !average surface and first level
                TSWALL(i,j)   = (TSROOF(i,j)+TSROAD(i,j))*0.5
                !average surface and first level
                TCANYON(i,j)  = TSWALL(i,j)
                !average surface and first level
                TSTOWN(i,j)   = TSWALL(i,j)
                EMISTOWN(i,j) =   0.8878604
                ALBTOWN(i,j)  =   0.129
                RCANYON(i,j)  = (1./((1000./rv(2,i,j))+1.))*1000.
                HT(i,j)       =   5.
                LET(i,j)      =   5.
                HIN(i,j)      =  10.
                LEIN(i,j)     =  40.

                call INI_TG_PROFILE(TSROOF(i,j), TIBLD(i,j), &
                     TC_ROOF(1:3),D_ROOF(1:3), TROOF(2:4,i,j))

                call INI_TG_PROFILE(TSROAD(i,j), TIROAD(i,j), &
                     TC_ROAD(1:3),D_ROAD(1:3), TROAD(2:4,i,j))

                call INI_TG_PROFILE(TSWALL(i,j), TIBLD(i,j), &
                     TC_WALL(1:3),D_WALL(1:3), TWALL(2:4,i,j))

             endif

          enddo

       enddo
    enddo

    return
  end subroutine TEB_INIT


  ! TEB
  subroutine TEBC_INIT(n2, n3, np, G_URBAN, EMIS, ALB, TS)
    ! Arguments:
    integer, intent(in) :: n2, n3, np
    real, intent(out)   :: G_URBAN(n2,n3,np)
    real, intent(out)   :: EMIS(n2,n3)
    real, intent(out)   :: ALB(n2,n3)
    real, intent(out)   :: TS(n2,n3)

    G_URBAN = 0.
    EMIS    = 0.
    ALB     = 0.
    TS      = 0.

  end subroutine TEBC_INIT
end module ModLeaf3Teb
