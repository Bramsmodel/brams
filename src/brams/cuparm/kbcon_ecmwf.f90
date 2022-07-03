module kbcon_ecmwf

contains 

  subroutine trigg_ecmwf(mentr_rate2,PO,DNO,KIDIA, KFDIA,mgmxp,mgmzp, &
       KLEV, PTENH,  PRMISH, PRMISSH, PGEOH, PAPH,&
       PQHFL,PAHFS, PTEN, PRMIS,PRMISS, PGEO,ierr,k22,kbcon,ktop,&
       m1,m2,m3,j,wu)  
    implicit none

    !  ROTINA CUBASEN - MECANISMO DE DISPARO DE CONVECCAO UTILIZADO NO MODELO DO ECMWF
    !  Referencia: Jakob, K. and Siebesma, P. A new subcloud Model for Mass-Flux Convection 
    !  Schemes: Influence on Triggering, updrafts properties, and Climate Model.
    !  Monthly Weather Review, v. 131, p. 2765-2778, 2003.
    !   
    !  Adaptada para o modelo CATT-BRAMS por Claudio Silva, Saulo Freitas, Nelson Pereira Filho e Rafael Mello
    !
    !  claudio.moises@cptec.inpe.br
    !
    integer           ::  KLON, KLEV, KIDIA, KFDIA, KTDIA,mgmxp,mgmzp
    integer           ::  JK, JL, JKE
    integer 	   ::  ierr(mgmxp), kbcon(mgmxp), ktop(mgmxp), aux(mgmxp), k22(mgmxp)
    integer           ::  m1, m2, m3, j!intent in
    real              ::  wu(mgmxp,mgmzp) !intent out
    real              ::  PO(mgmxp), DNO(mgmxp) 
    real		   ::  PTEN(mgmxp,mgmzp)   ! TEMPERATURA DO MCGAS (oC)
    real		   ::  PQEN(mgmxp,mgmzp)   ! UMIDADE ESPECIFICA (KG/KG)
    real		   ::  PTENL(mgmxp,mgmzp)  ! TEMPERATURA DO MCGAS (oC)
    real		   ::  PQENL(mgmxp,mgmzp)  ! UMIDADE ESPECIFICA (KG/KG)
    real		   ::  PQSEN(mgmxp,mgmzp)  ! UMIDADE ESPECIFICA DE SATURACAO (KG/KG)
    real		   ::  PQSENL(mgmxp,mgmzp) ! UMIDADE ESPECIFICA DE SATURACAO (KG/KG)
    real		   ::  PGEO(mgmxp,mgmzp)   ! ALTURA GEOPOTENCIAL (M) 
    real		   ::  PGEOL(mgmxp,mgmzp)  ! ALTURA GEOPOTENCIAL (M) 

    !srf REAL	   ::  PGEOH(mgmxp,mgmzp+1)  ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)
    real		   ::  PGEOH(mgmxp,mgmzp  )  ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)

    real		   ::  PGEOHL(mgmxp,mgmzp+1) ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)

    ! REAL		   ::  PGEOH(mgmxp,mgmzp)  ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)
    ! REAL		   ::  PGEOHL(mgmxp,mgmzp) ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)

    !srf REAL	   ::  PAPH(mgmxp,mgmzp+1)  ! PRESSAO (hPa)
    real		   ::  PAPH(mgmxp,mgmzp  )  ! PRESSAO (hPa)

    real		   ::  PAPHL(mgmxp,mgmzp+1) ! PRESSAO (hPa)
    real		   ::  PTENH(mgmxp,mgmzp)   ! TEMPERATURA (oC) NA GRADE DO MODELO DE NUVEM
    real		   ::  PTENHL(mgmxp,mgmzp)  ! TEMPERATURA (oC) NA GRADE DO MODELO DE NUVEM
    real		   ::  PQENH(mgmxp,mgmzp)   ! UMIDADE ESPECIFICA (KG/KG) NA GRADE DO MODELO DE NUVEM
    real		   ::  PQENHL(mgmxp,mgmzp)  ! UMIDADE ESPECIFICA (KG/KG) NA GRADE DO MODELO DE NUVEM
    real		   ::  PQSENH(mgmxp,mgmzp)  ! UMIDADE ESPECIFICA DE SATURACAO (KG/KG) NA GRADE DO MODELO DE NUVEM
    real		   ::  PQSENHL(mgmxp,mgmzp) ! UMIDADE ESPECIFICA DE SATURACAO (KG/KG) NA GRADE DO MODELO DE NUVEM
    real		   ::  PRMIS(mgmxp,mgmzp)   ! RAZAO DE MISTURA DE VAPOR (KG/KG)
    real		   ::  PRMISL(mgmxp,mgmzp)  ! RAZAO DE MISTURA DE VAPOR (KG/KG)
    real	           ::  PRMISS(mgmxp,mgmzp)  ! RAZAO DE MISTURA DE SATURACAO 
    real	           ::  PRMISSL(mgmxp,mgmzp) ! RAZAO DE MISTURA DE SATURACAO 
    real              ::  PSAT(mgmxp,mgmzp)    ! PRESSAO DE SATURACAO DO VAPOR (hPa)
    real              ::  PSATL(mgmxp,mgmzp)   ! PRESSAO DE SATURACAO DO VAPOR (hPa)
    logical	   ::  LLFLAG(mgmxp)       ! FLAG QUE INDICA A OCORRENCIA DE CONVECCAO CUMULUS
    real		   ::  LATI(mgmxp,mgmzp)    ! LATITUDE INICIAL EM GRAUS (PARA LER OS ARQUIVOS DE RADIOSSONDAGENS DA
    real		   ::  PUREL(mgmxp,mgmzp)   ! UMIDADE RELATIVA DO AR (%)
    real		   ::  PURELL(mgmxp,mgmzp)  ! UMIDADE RELATIVA DO AR (%)
    real  		   ::  PQHFL(mgmxp)        ! FLUXO DE UMIDADE (KG/M2S)   
    real		   ::  PAHFS(mgmxp)        ! FLUXO DE CALOR SENSIVEL (W/M2)
    logical	   ::  LDCUM(mgmxp)        ! FLAG PARA OCORRENCIA DE CONVECCAO DEEP-CUMULUS 
    logical	   ::  LDSC(mgmxp)         ! FLAG PARA OCORRENCIA DE CONVECAO SC
    integer	   ::  KCBOT(mgmxp)        ! NIVEL DA BASE DE NUVENS CUMULUS
    integer	   ::  KBOTSC(mgmxp)       ! NIVEL DA BASE DE NUVENS STRATO CUMULUS
    integer	   ::  KCTOP(mgmxp)        ! NIVEL DO TOPO DE NUVEM CUMULUS
    integer	   ::  KDPL(mgmxp)         ! EQUIVALENTE AO k22
    real		   ::  PCAPE(mgmxp)        ! CAPE (J/KG)
    real, parameter   ::  RCPD = 1004.64
    real              ::  PRMISSH(mgmxp,mgmzp)
    real              ::  PRMISH(mgmxp,mgmzp)
    real              ::  PRMISSHL(mgmxp,mgmzp)
    real              ::  PRMISHL(mgmxp,mgmzp)
    real              ::  mentr_rate2L(mgmxp,mgmzp)
    real              ::  mentr_rate2(mgmxp,mgmzp)
    real              ::  ZWU2H(mgmxp,mgmzp)

    do JL=KIDIA,KFDIA
       do JK=KLEV, 1, -1  ! INVERTIDO PARA FICAR COMPATIVEL COM O ESQUEMA DE GRELL !
          JKE = KLEV+1-JK
          PTENHL   (JL,JK)     =  PTENH     (JL,JKE)
          PRMISHL  (JL,JK)     =  PRMISH    (JL,JKE)
          PRMISSHL (JL,JK)     =  PRMISSH   (JL,JKE)
          PQENHL   (JL,JK)     =  PRMISHL   (JL,JK)/(1.+PRMISHL(JL,JK))   !CONVERTENDO DE "r" PARA "q" 
          PQSENHL  (JL,JK)     =  PRMISSHL  (JL,JK)/(1.+PRMISSHL(JL,JK)) 
          PGEOHL   (JL,JK)     =  9.8*PGEOH (JL,JKE)  !convertendo z para geopotencial *g
          PAPHL    (JL,JK)     =  100*PAPH  (JL,JKE)   ! convertende de hPa para Pa
          PTENL    (JL,JK)     =  PTEN      (JL,JKE)
          PRMISL   (JL,JK)     =  PRMIS     (JL,JKE)
          PRMISSL  (JL,JK)     =  PRMISS    (JL,JKE)
          PQENL    (JL,JK)     =  PRMISL    (JL,JK)/(1.+PRMISL (JL,JK))  
          PQSENL   (JL,JK)     =  PRMISSL   (JL,JK)/(1.+PRMISSL(JL,JK))
          PGEOL    (JL,JK)     =  9.8*PGEO  (JL,JKE)
          mentr_rate2L(JL,JK)  =  mentr_rate2(JL,JKE)
          !print*,'1',JL,JK, PTENHL   (JL,JK), PRMISHL  (JL,JK),PRMISSHL (JL,JK),PQENHL   (JL,JK)  


       enddo
       PAPHL(JL,KLEV+1)    =  PO(JL)
       PGEOHL(JL,KLEV+1)   =  0.9*PGEOHL(JL,klev)
       PQHFL(JL)           =  PQHFL(JL) 
       PAHFS(JL)           =  RCPD*DNO(JL)*PAHFS(JL)
       !print*,'11',JL,JK, PAPHL(JL,KLEV+1),PGEOHL(JL,KLEV+1) , PQHFL(JL),  PAHFS(JL) 

    enddo


    call CUBASEN(mentr_rate2L, PO, DNO, KIDIA, KFDIA, mgmxp,mgmzp,  KLEV, PTENHL,PQENHL,PQSENHL,PGEOHL, PAPHL,&
         & PQHFL, PAHFS, PTENL, PQENL,PQSENL, PGEOL,ierr,k22,kbcon,ktop,m1,m2,m3,j,wu)  



  end subroutine trigg_ecmwf

  !------------------------------------ subrotina cubasen (Jakob e Siebesma, 2003)


  subroutine CUBASEN(mentr_rate2,PO, DNO, KIDIA, KFDIA, mgmxp,mgmzp, &
       KLEV, PTENH,  PQENH, PQSENH, PGEOH, PAPH,&
       PQHFL, PAHFS, PTEN, PQEN,PQSEN, PGEO,ierr,k22,kcbot,kctop,m1,m2,m3,j,wu)  
    implicit none


    !          THIS ROUTINE CALCULATES CLOUD BASE FIELDS
    !          CLOUD BASE HEIGHT AND CLOUD TOP HEIGHT

    !          A. Pier Siebesma   KNMI ********      
    !          modified C Jakob (ECMWF) (01/2001) 
    !          modified P Bechtold (ECMWF) (08/2002) 
    !          (include cycling over levels to find unstable departure/base level+
    !           mixed layer properties +w Trigger)

    !          PURPOSE.
    !          --------
    !          TO PRODUCE CLOUD BASE AND CLOUD TOP VALUES FOR CU-PARAMETRIZATION

    !          INTERFACE
    !          ---------
    !          THIS ROUTINE IS CALLED FROM *CUMASTR*.
    !          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
    !          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
    !                 KLAB=0 FOR STABLE LAYERS
    !                 KLAB=1 FOR SUBCLOUD LEVELS
    !                 KLAB=2 FOR CLOUD LEVELS LEVEL

    !          METHOD.
    !          --------
    !          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP
    !          (ENTRAINING PLUME, WITH ENTRAINMENT PROPORTIONAL TO (1/Z))

    !     PARAMETER     DESCRIPTION                                   UNITS
    !     ---------     -----------                                   -----
    !     INPUT PARAMETERS (INTEGER):

    !    *KIDIA*        START POINT
    !    *KFDIA*        END POINT
    !    *KLON*         NUMBER OF GRID POINTS PER PACKET
    !    *KTDIA*        START OF THE VERTICAL LOOP
    !    *KLEV*         NUMBER OF LEVELS

    !    INPUT PARAMETERS (REAL):

    !    not used at the moment because we want to use linear intepolation
    !    for fields on the half levels.

    !    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
    !    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
    !    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
    !    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
    !    *PSSTRU*       KINEMATIC surface U-MOMENTUM FLUX             (M/S)^2
    !    *PSSTRV*       KINEMATIC surface V-MOMENTUM FLUX             (M/S)^2
    !    *PWN*          NORMALIZED LARGE-SCALE VERTICAL VELOCITY      (M/S)
    !    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
    !    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
    !    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
    !    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
    !    *PGEO*         GEOPOTENTIAL                                  M2/S2
    !    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
    !    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
    !    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
    !    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2

    !    UPDATED PARAMETERS (REAL):

    !    *PTU*          TEMPERATURE IN UPDRAFTS                         K
    !    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
    !    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
    !    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
    !    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S

    !    UPDATED PARAMETERS (INTEGER):

    !    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
    !                        KLAB=2 FOR CLOUD LEVELS

    !    OUTPUT PARAMETERS (LOGICAL):

    !    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
    !    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST

    !    OUTPUT PARAMETERS (INTEGER):

    !    *KCBOT*       CLOUD BASE LEVEL !    
    !    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL 
    !                  WITH A NON-ZERO CLOUD UPDRAFT.
    !    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS
    !    *KDPL*        DEPARTURE LEVEL
    !    *PCAPE*       PSEUDOADIABATIQUE max CAPE (J/KG)

    !          EXTERNALS
    !          ---------
    !          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT

    !          MODIFICATIONS
    !          -------------
    !             92-09-21 : Update to Cy44      J.-J. MORCRETTE
    !             02-11-02 : Use fixed last possible departure level and 
    !                        last updraft computation level for bit-reproducibility
    !                                            D.Salmond &  J. Hague
    !             03-07-03 : Tuning for p690     J. Hague
    !        M.Hamrud      01-Oct-2003 CY28 Cleaning
    !----------------------------------------------------------------------

    real, parameter  :: RCPD   = 1004.64
    real, parameter  :: RETV   = 0.608 
    real, parameter  :: RD     = 287.04 
    real, parameter  :: RV     = 461.50  
    real, parameter  :: RG     = 9.8200
    real, parameter  :: RKAP   = 0.41
    integer          :: NJKT1  
    integer          :: NJKT2 
    real, parameter  :: RDEPTHS = 20000.
    real, parameter  :: ENTRPEN = 1.0E-4
    real, parameter  :: ENTRSCV = 1.E-3
    logical          :: LMFDUDV = .false. 
    real, parameter  :: RLMIN = 1E-5 
    real, parameter  :: RESTT = 611.14
    real, parameter  :: RMD = 28.9644
    real, parameter  :: RMV = 18.0153
    real, parameter  :: RTT = 273.16
    real, parameter  :: R2ES = RESTT*(RD/RV)
    real, parameter  :: R3LES = 17.269 
    real, parameter  :: R3IES = 21.875
    real, parameter  :: R4LES = 35.86 
    real, parameter  :: R4IES = 7.66
    real, parameter  :: R5LES = R3LES*(273.16-R4LES)
    real, parameter  :: R5IES = R3IES*(273.16-R4IES)
    real, parameter  :: RKBOL = 1.380658E-23
    real, parameter  :: RNAVO = 6.0221367E+23
    real, parameter  :: R = RNAVO*RKBOL
    real, parameter  :: RLVTT = 2.5008E+6
    real, parameter  :: RLSTT = 2.8345E+6
    real, parameter  :: R5ALVCP = R5LES*RLVTT/RCPD
    real, parameter  :: R5ALSCP = R5IES*RLSTT/RCPD 
    real, parameter  :: RALVDCP = RLVTT/RCPD
    real, parameter  :: RLMLT=RLSTT-RLVTT
    real, parameter  :: RALSDCP = RLSTT/RCPD
    real, parameter  :: RALFDCP = RLMLT/RCPD
    real, parameter  :: RTWAT = 273.16 
    real, parameter  :: RTBER = RTWAT-5.
    real, parameter  :: RTICE = RTT-0.1
    real, parameter  :: RTICECU = RTT-23.0
    real, parameter  :: RTWAT_RTICECU_R  = 1./(RTWAT-RTICECU)
    real, parameter  :: RTWAT_RTICE_R = 1./(RTWAT-RTICE)


    integer :: KLON, mgmxp,mgmzp
    integer :: KLEV 
    integer :: KIDIA 
    integer :: KFDIA 
    integer :: KTDIA ! Argument NOT used
    real    :: PTENH(mgmxp,mgmzp) 

    real    :: PQENH(mgmxp,mgmzp)

    real    :: PGEOH(mgmxp,mgmzp+1) 
    real    :: PAPH(mgmxp,mgmzp+1) 
    real    :: PQHFL(mgmxp) 
    real    :: PAHFS(mgmxp) 
    real    :: PSSTRU(mgmxp) ! Argument NOT used
    real    :: PSSTRV(mgmxp) ! Argument NOT used
    real    :: PWN(mgmxp,mgmzp) ! Argument NOT used 
    real    :: PTEN(mgmxp,mgmzp)  
    real    :: PQEN(mgmxp,mgmzp) 
    real    :: PGEO(mgmxp,mgmzp) 
    real    :: PUEN(mgmxp,mgmzp) 
    real    :: PVEN(mgmxp,mgmzp)  
    real    :: PTU(mgmxp,mgmzp) 
    real    :: PQU(mgmxp,mgmzp)  
    real    :: PLU(mgmxp,mgmzp) 
    real    :: PUU(mgmxp,mgmzp) 
    real    :: PVU(mgmxp,mgmzp) 
    real    :: PWUBASE(mgmxp) 
    integer :: KLAB(mgmxp,mgmzp) 
    logical :: LDCUM(mgmxp) 
    logical :: LDSC(mgmxp) 
    integer :: KCBOT(mgmxp) 
    integer :: KBOTSC(mgmxp) 
    integer :: KCTOP(mgmxp) 
    integer :: KDPL(mgmxp) 
    real    :: PCAPE(mgmxp) 
    integer :: ICTOP(mgmxp)
    integer :: ICBOT(mgmxp)
    integer :: IBOTSC(mgmxp)
    integer :: ILAB(mgmxp,mgmzp)
    integer :: IDPL(mgmxp)  
    integer :: ierr(mgmxp)   ! flag que passa ldcum (logica) para o inteiro correspondente 
    ! correspondente do esquema de grell, if ierr = 0, LDCUM = .T. --
    ! CLAUDIO SILVA
    integer  :: AUX(mgmxp) ! para inverter os niveis kctop e kcbot
    integer  :: kbcon(mgmxp), ktop(mgmxp), k22(mgmxp)

    real  :: PQSEN(mgmxp,mgmzp)
    real  :: PQSENH(mgmxp,mgmzp)
    real  :: PRMISS(mgmxp,mgmzp)
    real  :: PRMIS(mgmxp,mgmzp)
    real  :: PSAT(mgmxp,mgmzp)
    real  :: PVAP(mgmxp,mgmzp)
    real  :: PO(mgmxp), DNO(mgmxp) 
    real  :: PRMISSH(mgmxp,mgmzp)
    real  :: PRMISH(mgmxp,mgmzp)

    real  :: mentr_rate2(mgmxp,mgmzp) ! TAXA DE ENTRANHAMENTO DE GRELL

    integer           ::  m1, m2, m3, j, jke, jkes !intent in
    real              ::  wu(mgmxp,mgmzp)
    !real              ::  wu_C(m1,m2,m3) ! matriz auxiliar para calcular a 
    !                                             ! velocidade vertical - e a energia cinetica

    !intent out
    !             LOCAL STORAGE
    !             ----- -------

    logical ::         LL_LDBASE(mgmxp),&
         & LLGO_ON(mgmxp),&
         & LLDEEP(mgmxp),    LLDCUM(mgmxp), &!*UPG PB
         & LLDSC(mgmxp),     LLFIRST(mgmxp)  

    logical ::     LLRESET,        LLRESETJL(mgmxp)
    integer :: ICALL, IK, IKB, IS, JK, JL, JKK, JKT, JKB !*UPG PB
    integer :: JKT1
    integer :: JKT2   
    real    :: ZS(mgmxp,mgmzp)
    real    :: ZSENH(mgmxp,mgmzp+1)
    real    :: ZQENH(mgmxp,mgmzp+1)
    real    :: ZSUH (mgmxp,mgmzp)
    real    :: ZWU2H(mgmxp,mgmzp)
    real    :: ZBUOH(mgmxp,mgmzp)  
    real    :: ZQOLD(mgmxp),ZPH(mgmxp)
    real    :: ZMIX(mgmxp)
    real    :: ZDZ(mgmxp)
    real    :: ZCBASE(mgmxp)
    real    :: ZLU(mgmxp,mgmzp)
    real    :: ZQU(mgmxp,mgmzp)
    real    :: ZTU(mgmxp,mgmzp)
    real    :: ZUU(mgmxp,mgmzp)
    real    :: ZVU(mgmxp,mgmzp)  
    real    :: WVEL(mgmxp,mgmzp)
    real    :: ZCAPE(mgmxp,mgmzp) ! local for CAPE at every departure level
    real    :: ZBUOF, ZZ, ZC2, ZEPSADD
    real    :: ZRHO      ! DENSITY AT SURFACE (KG/M^3) 
    real    :: ZKHVFL    ! SURFACE BUOYANCY FLUX (K M/S)
    real    :: ZWS       ! SIGMA_W AT LOWEST MODEL HALFLEVEL (M/S)
    real    :: ZQEXC,ZQEXCS(mgmxp)! HUMIDITY EXCESS AT LOWEST MODEL HALFLEVEL (KG/KG)
    real    :: ZTEXC,ZTEXCS(mgmxp)! TEMPERATURE EXCESS AT LOWEST MODEL HALFLEVEL (K)
    real    :: ZEPS      ! FRACTIONAL ENTRAINMENT RATE   [M^-1]
    real    :: ZTVENH    ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)  
    real    :: ZTVUH     ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
    real    :: ZLGLAC    ! UPDRAFT LIQUID WATER FROZEN IN ONE LAYER
    real    :: zqsu
    real    :: ZCOR
    real    :: zdq
    real    :: zalfaw
    real    :: zfacw
    real    :: zfaci
    real    :: zfac
    real    :: zesdp
    real    :: zdqsdt
    real    :: zdtdp
    real    :: zdp
    real    :: zpdifftop
    real    :: zpdiffbot
    real    :: ZSF
    real    :: ZQF
    real    :: zaw
    real    :: zbw  
    real    :: ZTVEN1
    real    :: ZTVEN2
    real    :: ZTVU1
    real    :: ZTVU2 ! pseudoadiabatique T_v
    real    :: ZDTVTRIG(mgmxp) ! virtual temperatures
    real    :: ZWORK1
    real    :: ZWORK2! work arrays for T and w perturbations
    real    :: ZRCPD
    real    :: ZRG
    real    :: ZTMP
    real    :: ZTINY
    real    :: TINY
    real    :: ZLUMIN
    integer :: INDB, INDT ! PARA LIMITAR O NJKT1 E O NJKT2

    real :: ZWS1D(mgmxp),ZRHO1D(mgmxp)!srf 

    !----------------------------------------------------------------------
    !     0.           INITIALIZE CONSTANTS AND FIELDS
    !                  -------------------------------
    !----------------------------------------------------------------------
    !      CONTEUDO MINIMO DE AGUA PARA ENCONTRAR A BASE DA NUVEM! CLAUDIO SILVA


    ZLUMIN = 5.e-4
    ZC2    = 0.55
    ZAW    = 1.0
    ZBW    = 2.0
    ZEPSADD= 1.E-4

    do JL=KIDIA,KFDIA
       PWUBASE(JL)=0.0
       LLGO_ON(JL)=.true.
       LLFIRST(JL)=.true.
       KDPL(JL)=KLEV
       ZTEXCS(JL)=0.2
       ZQEXCS(JL)=1.E-4
    enddo

!!!!!!!!!!!!!!! Para determinar os limites de busca da base e topo ! 

    INDB = 0
    INDT = 0
    do JL = KIDIA, KFDIA 
       do JK = 1, KLEV
	  if(PAPH(JL,JK).gt.60.and.INDB.eq.0) then
             JKT2 = JK   
             INDB = 1
	  endif

          if(PAPH(JL,JK).gt.350E2.and.INDT.eq.0) then  ! o limite inferior de busca da base
             JKT1 = JK                                 ! da nuvem esta muito profundo,
             ! testar 600
             INDT = 1
          endif

       enddo
    enddo

!!!!!!!!!!!!!!!!!!!

    ZRG=1.0/RG
    ZRCPD=1.0/RCPD
    ZTINY=tiny(ZRG)

    do JK=1,KLEV
       do JL=KIDIA,KFDIA
          ZTU(JL,JK) = PTENH(JL,JK) !PTU(JL,JK)
          ZQU(JL,JK) = PQENH(JL,JK) !PQU(JL,JK)
          ZLU(JL,JK) = PLU(JL,JK)
          ZUU(JL,JK) = PUU(JL,JK)
          ZVU(JL,JK) = PVU(JL,JK)
          ILAB(JL,JK)= KLAB(JL,JK)
          ZCAPE(JL,JK)= 0.0
       enddo
    enddo
    !----------------------------------------------------------------------
    !       -----------------------------------------------------------
    !       1.1  PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
    !             OF SPECIFIC HUMIDITY AND STATIC ENERGY
    !       -----------------------------------------------------------

    do JK=1,KLEV
       do JL=KIDIA,KFDIA
          ZWU2H(JL,JK) = 0.0
          ZS   (JL,JK) = RCPD*PTEN(JL,JK) + PGEO(JL,JK)
          ZQENH(JL,JK) = PQENH(JL,JK)
          ZSENH(JL,JK) = RCPD*PTENH(JL,JK)+PGEOH(JL,JK)
       enddo
    enddo


    do JKK = KLEV,JKT1,-1  ! Big external loop for level testing:
       ! find first departure level that produces deepest cloud top
       ! or take surface level for shallow convection and Sc
       ! 
       !        ---------------------------------------------------------
       !        1.2    INITIALISE FIELDS AT DEPARTURE HALF MODEL LEVEL
       !        ---------------------------------------------------------
       !
       IS=0
       do JL=KIDIA,KFDIA
          if (LLGO_ON(JL)) then
             IS=IS+1
             IDPL(JL)      = JKK      ! departure level
             ICBOT  (JL)   = JKK      ! cloud base level for convection, (-1 if not found)
             IBOTSC (JL)   = KLEV-1   ! sc    base level for sc-clouds , (-1 if not found)
             ICTOP(JL)     = KLEV-1   ! cloud top for convection (-1 if not found) --->>> Mas o topo fica embaixo ??
             LLDCUM(JL)    = .false.  ! on exit: true if cloudbase=found
             LLDSC (JL)    = .false.  ! on exit: true if cloudbase=found
             LL_LDBASE(JL) =.false.   ! on exit: true if cloudbase=found
             ZDTVTRIG(JL)  = 0.0
             ZUU(JL,JKK)   = PUEN(JL,JKK)*(PAPH(JL,JKK+1)-PAPH(JL,JKK))
             ZVU(JL,JKK)   = PVEN(JL,JKK)*(PAPH(JL,JKK+1)-PAPH(JL,JKK))
          endif
       enddo

       if(IS /= 0) then

          if(JKK == KLEV) then

             do JL=KIDIA,KFDIA
                if (LLGO_ON(JL)) then
                   ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1.+RETV*PQEN(JL,JKK))))
                   !- buoyancy flux (H+LE)
                   ZKHVFL= (PAHFS(JL)/RCPD+RETV*PTEN(JL,JKK)*PQHFL(JL))/ZRHO
                   !-convective-scale velocity w*
                   ZWS=0.001-1.5*RKAP*ZKHVFL*PGEOH(JL,KLEV)/PTEN(JL,KLEV)
                   !srf
                   ZWS1D(JL)=1.2*ZWS**.3333
                   ZRHO1D(JL)=ZRHO
                   !srf

                   if(ZWS >= ZTINY) then
                      !-convective-scale velocity w*
                      ZWS=1.2*ZWS**.3333
                      ILAB(JL,JKK)= 1
                      !- temperature excess 
                      ZTEXC     = max(-1.5*PAHFS(JL)/(ZRHO*ZWS*RCPD),0.0)
                      !- moisture  excess
                      ZQEXC     = max(-1.5*PQHFL(JL)/(ZRHO*ZWS),0.0)

                      ZTEXCS(JL)=ZTEXC
                      ZQEXCS(JL)=ZQEXC

                      !- initial values for updrafts
                      !- humidty
                      ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEXC
                      !- static energy (?)
                      ZSUH (JL,JKK) = ZSENH(JL,JKK) + RCPD*ZTEXC
                      !- temperature
                      ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD + ZTEXC
                      !- kinetic energy
                      ZWU2H(JL,JKK) = ZWS**2				   
                      !-  determine buoyancy at lowest half level
                      ZTVENH            = (1.0+RETV*ZQENH(JL,JKK)) &
                           & *(ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD  
                      ZTVUH             = (1.0+RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
                      ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH

                   else
                      LLGO_ON(JL)=.false.      ! non-convective point

                   endif
                endif
             enddo

          else ! se jkk /= klev

             do JL=KIDIA,KFDIA
                if (LLGO_ON(JL)) then
                   !- air density
                   ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1.+RETV*PQEN(JL,JKK))))
                   ILAB(JL,JKK)= 1
                   !- temperature excess
                   ZTEXC=.2
                   !- moisture  excess
                   ZQEXC=1.E-4

                   !srf--------------------------------------------------- 
                   !- temperature excess 			     
                   if(jkk == klev -1) then
                      ZTEXC     = min(-1.5*PAHFS(JL)/(ZRHO1D(JL)*ZWS1d(JL)*RCPD),0.2)
                      !- moisture  excess
                      ZQEXC     = min(-1.5*PQHFL(JL)/(ZRHO1D(JL)*ZWS1D(JL)),1.E-4)
                   endif
                   !srf--------------------------------------------------- 

                   !- initial values for updrafts
                   !- humidty
                   ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEXC
                   !- static energy (?)
                   ZSUH (JL,JKK) = ZSENH(JL,JKK) + RCPD*ZTEXC
                   !- temperature
                   ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))*ZRCPD + ZTEXC
                   ! construct mixed layer for parcels emanating in lowest 60 hPa
                   ! Esabilizando ainda mais as parcelas, fazendo a Camada de mistura ateh 80 hPa
                   if (PAPH(JL,KLEV+1)-PAPH(JL,JKK-1)<80.E2) then
                      ZQU(JL,JKK) =0.0
                      ZSUH(JL,JKK)=0.0
                      ZWORK1      =0.0
                      do JK=JKK+1,JKK-1,-1
                         if( ZWORK1 < 50.E2 ) then
                            ZWORK2=PAPH(JL,JK)-PAPH(JL,JK-1)
                            ZWORK1      =ZWORK1+ZWORK2
                            ZQU(JL,JKK) =ZQU(JL,JKK) +ZQENH(JL,JK)*ZWORK2
                            ZSUH(JL,JKK)=ZSUH(JL,JKK)+ZSENH(JL,JK)*ZWORK2
                         endif
                      enddo

                      ZQU(JL,JKK) =ZQU(JL,JKK) /ZWORK1+ZQEXC
                      ZSUH(JL,JKK)=ZSUH(JL,JKK)/ZWORK1+RCPD*ZTEXC
                      ZTU(JL,JKK) =(ZSUH(JL,JKK)-PGEOH(JL,JKK))/RCPD+ZTEXC
                   endif
                   !- kinetic energy
                   ZWU2H(JL,JKK) = 1.0 
                   !
                   !  determine buoyancy at lowest half level
                   !
                   ZTVENH            = (1.0 +RETV*ZQENH(JL,JKK)) &
                        & *(ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD  
                   ZTVUH             = (1.0 +RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
                   ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH
                endif
             enddo
          endif
       endif

       !----------------------------------------------------------------------
       !     2.0          DO ASCENT IN SUBCLOUD AND LAYER,
       !                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
       !                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
       !                  CHECK FOR BUOYANCY AND SET FLAGS
       !                  -------------------------------------
       !       ------------------------------------------------------------
       !        1.2  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
       !       ------------------------------------------------------------
       do JK=JKK-1,JKT2,-1
          IS=0
          if(JKK==KLEV) then ! 1/z mixing for shallow
             do JL=KIDIA,KFDIA
                if (LLGO_ON(JL)) then
                   IS         = IS+1
                   ZDZ(JL)        = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG

                   zeps = mentr_rate2(jl,jk)

                   !ZEPS = ZC2/(PGEO(JL,JK)*ZRG + ZDZ(jl)) + ZEPSADD
                   !mentr_rate2(jl,jk) = ZC2/(PGEO(JL,JK)*ZRG + ZDZ(jl)) + ZEPSADD
                   !mentr_rate2(jl,jk)=zeps

                   ZMIX(JL)   = 0.5*ZDZ(JL)*ZEPS
                   ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5
                   ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5
                   ZTMP = 1.0/(1.0+ZMIX(JL))
                   ZQU(JL,JK)= (ZQU(JL,JK+1)*(1.0-ZMIX(JL))&
                        & +2.0*ZMIX(jl)*ZQF) * ZTMP  
                   ZSUH (JL,JK)= (ZSUH(JL,JK+1)*(1.0-ZMIX(JL))&
                        & +2.0*ZMIX(jl)*ZSF) * ZTMP  
                   ZQOLD(JL)  = ZQU(JL,JK)
                   ZTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
                   ZPH  (JL)    = PAPH(JL,JK)
                endif
             enddo
          else

             do JL=KIDIA,KFDIA
                if (LLGO_ON(JL)) then
                   IS         = IS+1
                   ZDZ(JL)        = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
                   ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5
                   ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5

                   !srf			      ZMIX(JL)= ENTRPEN           *(PGEOH(JL,JK)-PGEOH(JL,JK+1))/RG
                   ZMIX(JL)= mentr_rate2(jl,jk)*(PGEOH(JL,JK)-PGEOH(JL,JK+1))/RG

                   ZQU(JL,JK)= ZQU(JL,JK+1)*(1.0 -ZMIX(JL))+ ZQF*ZMIX(JL)
                   ZSUH(JL,JK)= ZSUH(JL,JK+1)*(1.0 -ZMIX(JL))+ ZSF*ZMIX(JL)
                   ZQOLD(JL)  = ZQU(JL,JK)
                   ZTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
                   ZPH  (JL)    = PAPH(JL,JK)
                endif
             enddo
          endif


          if (IS == 0) exit
          IK=JK
          ICALL=1
          call CUADJTQ(KIDIA,KFDIA,mgmxp,mgmzp,KTDIA,KLEV,IK,ZPH,ZTU,ZQU,LLGO_ON,ICALL)  
          do JL=KIDIA,KFDIA
             if(LLGO_ON(JL)) then

                ! add condensation to water
                ZDQ=max(ZQOLD(JL)-ZQU(JL,JK),0.0)

                ! O CALCULO DA BASE DA NUVEM DEPENDE EXCLUSIAMENTE DE ZLU !
                ZLU(JL,JK)=ZLU(JL,JK+1)+ZDQ

                ! freezing
                ZLGLAC=ZDQ*((1.0-FOEALFCU(ZTU(JL,JK)))-&
                     & (1.0-FOEALFCU(ZTU(JL,JK+1))))  

                ! pseudo-microphysics
                if(JKK==KLEV) then  ! no precip for shallow
                   ZLU(JL,JK)=min(ZLU(JL,JK),5.E-3)

                   !* chose a more pseudo-adiabatic formulation as original overestimates
                   !* water loading efect and therefore strongly underestimates cloud thickness

                else
                   ZLU(JL,JK)=0.5*ZLU(JL,JK) 
                endif

                ! update dry static energy after condensation + freezing
                ZSUH(JL,JK) = RCPD*(ZTU(JL,JK)+RALFDCP*ZLGLAC)+PGEOH(JL,JK)

                ! Buoyancy on half and full levels
                ZTVUH           = (1.0+RETV*ZQU(JL,JK)-ZLU(JL,JK))*ZTU(JL,JK)&
                     & +RALFDCP*ZLGLAC  

                ZTVENH          = (1.0+RETV*ZQENH(JL,JK)) &
                     & *(ZSENH(JL,JK)-PGEOH(JL,JK))*ZRCPD  

                ZBUOH(JL,JK)   = (ZTVUH-ZTVENH)*RG/ZTVENH
                ZBUOF          = (ZBUOH(JL,JK) + ZBUOH(JL,JK+1))*0.5

                ! solve kinetic energy equation
                ZTMP=1.0/(1.0+2.0*ZBW*ZMIX(jl))
                ZWU2H(JL,JK) = (ZWU2H(JL,JK+1)*(1.0-2.0*ZBW*ZMIX(jl))&
                     & +2.0*ZAW*ZBUOF*ZDZ(jl))*ZTMP 

                ! compute pseudoadiabatique CAPE for diagnostics
                ZTVU2 = ZTU(JL,JK)  *(1.0 +RETV*ZQU(JL,JK))
                ZTVEN2= PTENH(JL,JK)*(1.0 +RETV*PQENH(JL,JK))

                if (JK == JKK-1) then
                   ZTVU1  = ZTVU2
                   ZTVEN1 = ZTVEN2
                endif
                ZBUOF = (ZTVU2+ZTVU1-ZTVEN1-ZTVEN2)/ZTVEN2
                ZBUOF = ZBUOF*ZDZ(JL)*RG
                ZCAPE(JL,JKK)  = ZCAPE(JL,JKK) + max(0.0,ZBUOF)
                ZTVU1=ZTVU2
                ZTVEN1=ZTVEN2


                if(ZLU(JL,JK) >ZLUMIN.and.ILAB(JL,JK+1)==1) then
                   IK=JK+1
                   ZQSU=FOEEWM(ZTU(JL,IK))/PAPH(JL,IK)
                   ZQSU=min(0.5,ZQSU)
                   ZCOR=1.0/(1.0-RETV*ZQSU)
                   ZQSU=ZQSU*ZCOR
                   ZDQ=min(0.,ZQU(JL,IK)-ZQSU)
                   ZALFAW=FOEALFA(ZTU(JL,IK))
                   ZFACW=R5LES/((ZTU(JL,IK)-R4LES)**2)
                   ZFACI=R5IES/((ZTU(JL,IK)-R4IES)**2)
                   ZFAC=ZALFAW*ZFACW+(1.-ZALFAW)*ZFACI
                   ZESDP=FOEEWM(ZTU(JL,IK))/PAPH(JL,IK)
                   ZCOR=1.0/(1.0-RETV*ZESDP)
                   ZDQSDT=ZFAC*ZCOR*ZQSU
                   ZDTDP=RD*ZTU(JL,IK)/(RCPD*PAPH(JL,IK))
                   ZDP=ZDQ/(ZDQSDT*ZDTDP)
                   ZCBASE(JL)=PAPH(JL,IK)+ZDP

                   ! chose nearest half level as cloud base
                   ZPDIFFTOP=ZCBASE(JL)-PAPH(JL,JK)
                   ZPDIFFBOT=PAPH(JL,JK+1)-ZCBASE(JL)

                   if(ZPDIFFTOP > ZPDIFFBOT.and.ZWU2H(JL,JK+1)>0.0) then
                      JKB=min(KLEV-1,JK+1)
                      ILAB(JL,JKB)=2 !*UPG
                      ILAB(JL,JK)=2
                      LL_LDBASE(JL) =.true.
                      LLDSC(JL)   =.true.
                      IBOTSC(JL) =JKB
                      ICBOT(JL)  =JKB
                      ZLU(JL,JK+1) = RLMIN

                   elseif(ZPDIFFTOP <= ZPDIFFBOT.and.ZWU2H(JL,JK)>0.0) then
                      ILAB(JL,JK)=2
                      LL_LDBASE(JL) =.true.
                      LLDSC(JL)   =.true.
                      IBOTSC(JL) =JK
                      ICBOT(JL)  =JK

                   endif
                   JKB=ICBOT(JL)
        	endif

         ! decide on presence of convection, cloud base and cloud top based on
         ! kinetic energy
        	if (ZWU2H(JL,JK) < 0.0) then
                   LLGO_ON(JL) = .false. 
                   if (ZLU(JL,JK+1)>0.0) then  ! AQUI ZLU PRECISA SER MAIOR QUE ZERO
                      ICTOP(JL)   = JK
                      LLDCUM(JL)   = .true.
                   else
                      LLDCUM(JL)   = .false.
                   endif
        	else
                   if (ZLU(JL,JK)>ZLUMIN) then
                      ILAB(JL,JK) = 2
                   else
                      ILAB(JL,JK) = 1
                   endif
        	endif
             endif
          enddo

          if(LMFDUDV.and.JKK==KLEV) then
             do JL=KIDIA,KFDIA
                if(.not.LL_LDBASE(JL).and.LLGO_ON(JL)) then
                   ZUU(JL,JKK)=ZUU(JL,JKK)+PUEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
                   ZVU(JL,JKK)=ZVU(JL,JKK)+PVEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
                endif
             enddo
          endif

       enddo


       if( JKK==KLEV) then
          ! set values for departure level for PBL clouds = first model level
          do JL=KIDIA,KFDIA
             LDSC(JL)  = LLDSC(JL)
             if(LDSC(JL)) then
                KBOTSC(JL)= IBOTSC(JL)
             else
                KBOTSC(JL)=-1
             endif

             LLGO_ON(JL) = .false.
             JKT=ICTOP(JL)
             JKB=ICBOT(JL)
             LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>RDEPTHS

             if(LLDEEP(JL)) LLDCUM(JL)=.false. ! no deep allowed for KLEV
             lldeep(jl)=.false.                ! for deep convection start only at level KLEV-1
             ! and form mixed layer, so go on
             ! test further for deep convective columns as not yet found

             if ( LLDEEP(JL) ) LLFIRST(JL)=.false.
             LLGO_ON(JL) = .not.LLDEEP(JL)

             if(LLDCUM(JL)) then
                KCBOT(JL)= ICBOT(JL)
                KCTOP(JL)= ICTOP(JL)
                KDPL(JL)  = IDPL(JL)
                LDCUM(JL) = LLDCUM(JL)
                PWUBASE(JL)=sqrt(max(ZWU2H(JL,JKB),0.0))
             else
                KCTOP(JL)=-1
                KCBOT(JL)=-1
                KDPL(JL) =KLEV-1
                LDCUM(JL)=.false.
                PWUBASE(JL)=0.0
             endif
          enddo
          do JK=KLEV,1,-1
             do JL=KIDIA,KFDIA
                JKT=ICTOP(JL)
                if ( JK>=JKT ) then
                   KLAB(JL,JK)=ILAB(JL,JK)
                   PTU(JL,JK)=ZTU(JL,JK)
                   PQU(JL,JK)=ZQU(JL,JK)
                   PLU(JL,JK)=ZLU(JL,JK)
                endif
             enddo
          enddo
          if(LMFDUDV) then
             do JL=KIDIA,KFDIA
                if(LDCUM(JL)) then
                   IKB=KCBOT(JL)
                   ZZ=1.0 /(PAPH(JL,JKK+1)-PAPH(JL,IKB))
                   PUU(JL,JKK)=ZUU(JL,JKK)*ZZ
                   PVU(JL,JKK)=ZVU(JL,JKK)*ZZ
                endif
             enddo
          endif
       endif

       if( JKK < KLEV ) then
          LLRESET=.false.
          do JL=KIDIA,KFDIA
             if ( .not.LLDEEP(JL) ) then
                JKT=ICTOP(JL)
                JKB=ICBOT(JL)
                LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>=RDEPTHS 
             endif
             LLRESETJL(JL)=LLDEEP(JL).and.LLFIRST(JL)
             LLRESET=LLRESET.or.LLRESETJL(JL)
          enddo

          if(LLRESET) then
             do JK=KLEV,1,-1
                do JL=KIDIA,KFDIA
                   ! keep first departure level that produces deep cloud
                   if ( LLRESETJL(JL) ) then 
                      JKT=ICTOP(JL)
                      JKB=IDPL(JL)
                      if ( JK<=JKB .and. JK>=JKT ) then
                         KLAB(JL,JK)=ILAB(JL,JK)
                         PTU(JL,JK)=ZTU(JL,JK)
                         PQU(JL,JK)=ZQU(JL,JK)
                         PLU(JL,JK)=ZLU(JL,JK)
                      else 
                         KLAB(JL,JK)=1
                         PTU(JL,JK)=PTENH(JL,JK)
                         PQU(JL,JK)=PQENH(JL,JK)
                         PLU(JL,JK)=0.0
                      endif
                      if ( JK<JKT ) KLAB(JL,JK)=0
                   endif
                enddo
             enddo
          endif

          do JL=KIDIA,KFDIA
             if (LLDEEP(JL) .and. LLFIRST(JL)) then
                KDPL(JL)  = IDPL(JL)
                KCTOP(JL) = ICTOP(JL)
                KCBOT(JL) = ICBOT(JL)
                LDCUM(JL) = LLDCUM(JL)
                LDSC(JL)  = .false.
                KBOTSC(JL)= -1
                JKB=KCBOT(JL)
                PWUBASE(JL)=sqrt(max(ZWU2H(JL,JKB),0.0))

                !  no initialization of wind for deep here, this is done in
                !  CUINI and CUASCN
                LLFIRST(JL)=.false.
             endif
             LLGO_ON(JL) = .not.LLDEEP(JL)
          enddo
       endif
    enddo ! end of big loop for search of departure level     



    ! chose maximum CAPE value
    do JL=KIDIA,KFDIA
       PCAPE(JL) = maxval(ZCAPE(JL,:))
    enddo

!!!! PASSANDO O LDCUM DO ESQUEMA "CUBASEN" PARA IERR DO GRELL

    do JL = KIDIA, KFDIA
       if (LDCUM(JL)) then
          ierr(jl) = 0
       else
          ierr(jl) = 2
       endif

    enddo


!!!!!! INVERTENDO A ORDEM DE KCBOT E KCTOP

    do JL = KIDIA, KFDIA
       if(ierr(jl).eq.0) then 
          KCTOP(JL) = KLEV+1-KCTOP(JL)
          KCBOT(JL) = KLEV+1-KCBOT(JL)
          K22(JL)   = KLEV+1-KDPL(JL)
       else
          KCTOP(JL) = -1
          KCBOT(JL) = -1  
          K22 (JL)  = -1    
       endif
    enddo


    ! A variavel ZWU2H(JL,JK) corresponde a energia cinetica, 
    ! como pode ser visto pelo calculo da velocidade na base das nuvens
    ! PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0))
    ! portanto, faz-se necessario calcular esta velocidade para toda a coluna.
    ! mas e importante notar que nao podemos ter valores negativos da energia cinetica!

    do JL=KIDIA,KFDIA 
       if(ierr(jl).eq.0) then 
          do JK=KLEV,1,-1
             JKE = KLEV+1-JK  

             !wu_C(JKE,JL,J) = ZWU2H(JL,JK)

             !if(jke.ge.kctop(jl)) then ! zero a matriz do topo para cima, com isso 
             !wu_C(JKE,JL,J) = 0. ! excluimos a energia cinetica negativa (que o 
             !endif                     ! algoritmo permite, embora nao seja possivel 
             ! na realidade.
             !wu(JKE,JL,J) = sqrt(max(2*wu_C(JKE,JL,J),0.0))
             wu(JL,JKE) = sqrt(max(2*ZWU2H(JL,JK),0.0))

          enddo
       endif
       if(ierr(jl).eq.2) then 
          do JK=KLEV,1,-1
             JKE = KLEV+1-JK
             !wu_C(JKE,JL,J) = 0.0
             !wu(JKE,JL,J) = sqrt(max(2*wu_C(JKE,JL,J),0.0))
             wu(JL,JKE) = 0.
          enddo
       endif
    enddo





  end subroutine CUBASEN

  !!----------------------------------------------------------------
  !!-----------------SUBROTINA QUE AJUSTA OS CAMPOS DE UMIDADE E ETEMPERATURA
  !!----------------- DEVIDO A CONDENSACAO
  !!--------------------------------------------------------------------------------
  subroutine CUADJTQ(KIDIA,KFDIA,mgmxp,mgmzp,KTDIA,KLEV,KK,PSP,PT,PQ,LDFLAG,KCALL)
    implicit none

    !          M.TIEDTKE         E.C.M.W.F.     12/89 
    !          MODIFICATIONS
    !          -------------
    !          D.SALMOND         CRAY(UK))      12/8/91
    !          J.J. MORCRETTE    ECMWF          92-09-18   Update to Cy44
    !          J.F. MAHFOUF      ECMWF          96-06-11   Smoothing option
    !          PURPOSE.
    !          --------
    !          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT

    !          INTERFACE
    !          ---------
    !          THIS ROUTINE IS CALLED FROM SUBROUTINES:
    !              *COND*     (T AND Q AT CONDENSATION LEVEL)
    !              *CUBASE*   (T AND Q AT CONDENSATION LEVEL)
    !              *CUASC*    (T AND Q AT CLOUD LEVELS)
    !              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
    !              *CUSTRAT*  (T AND Q AT CONDENSATION LEVEL)
    !          INPUT ARE UNADJUSTED T AND Q VALUES,
    !          IT RETURNS ADJUSTED VALUES OF T AND Q

    !     PARAMETER     DESCRIPTION                                   UNITS
    !     ---------     -----------                                   -----
    !     INPUT PARAMETERS (INTEGER):

    !    *KIDIA*        START POINT
    !    *KFDIA*        END POINT
    !    *KLON*         NUMBER OF GRID POINTS PER PACKET
    !    *KTDIA*        START OF THE VERTICAL LOOP
    !    *KLEV*         NUMBER OF LEVELS
    !    *KK*           LEVEL
    !    *KCALL*        DEFINES CALCULATION AS
    !                      KCALL=0  ENV. T AND QS IN*CUINI*
    !                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
    !                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)

    !     INPUT PARAMETERS (LOGICAL):

    !    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)

    !     INPUT PARAMETERS (REAL):

    !    *PSP*          PRESSURE                                        PA

    !     UPDATED PARAMETERS (REAL):

    !    *PT*           TEMPERATURE                                     K
    !    *PQ*           SPECIFIC HUMIDITY                             KG/KG


    !          EXTERNALS   
    !          ---------
    !          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
    !          FOR CONDENSATION CALCULATIONS.
    !          THE TABLES ARE INITIALISED IN *SUPHEC*.

    !----------------------------------------------------------------------
    integer :: kidia,kfdia,mgmxp,mgmzp,ktdia,klev,kk,kcall
!    implicit logical (L)
    real, parameter :: RETV    = 0.608 
    real, parameter :: RTT     = 273.16 
    real, parameter :: RESTT   = 611.14
    real, parameter :: RMD     = 28.9644
    real, parameter :: RMV     = 18.0153
    real, parameter :: RV      = 461.50
    real, parameter :: RD      = 287.04
    real, parameter :: R2ES    = RESTT*(RD/RV)
    real, parameter :: R3LES   = 17.269 
    real, parameter :: R3IES   = 21.875
    real, parameter :: R4LES   = 35.86 
    real, parameter :: R4IES   = 7.66
    real, parameter :: R5LES   = R3LES*(273.16-R4LES)
    real, parameter :: R5IES   = R3IES*(273.16-R4IES)
    real, parameter :: RCPD    = 3.5*RD
    real, parameter :: RLVTT   = 2.5008E+6
    real, parameter :: RLSTT   = 2.8345E+6
    real, parameter :: R5ALVCP = R5LES*RLVTT/RCPD
    real, parameter :: R5ALSCP = R5IES*RLSTT/RCPD 
    real, parameter :: RALVDCP = RLVTT/RCPD
    real, parameter :: RLMLT   = RLSTT-RLVTT
    real, parameter :: RALSDCP = RLSTT/RCPD
    real, parameter :: RALFDCP = RLMLT/RCPD
    real, parameter :: RTWAT   = 273.16 
    real, parameter :: RTBER   = RTWAT-5.
    real, parameter :: RTICE   = RTT-0.1
    real, parameter :: RTICECU = RTT-23.0
    logical         :: LPHYLIN = .true. ! if not linearized physics is used
    real, parameter :: RLPTRC  = 266.425 
    real, parameter :: RLPAL1  = 0.15 
    real, parameter :: RLPAL2  = 20. 
    real, parameter :: ZQMAX   = 0.5
    logical    LDFLAG(mgmxp)
    real ::    PT(mgmxp,mgmzp)
    real ::    PQ(mgmxp,mgmzp)
    real ::    PSP(mgmxp)
    real ::    ZCOND(mgmxp)
    real ::    ZQP(mgmxp)
    real ::    ZQSAT1

    integer :: isum, jl
    real :: z1s, z2s, zcond1, zcor, zcoeewi, z0alfa
    real :: zfoeewi, zfoeewl, zoealfa, zqsat, ztarg

    !*********************************************
    if (LPHYLIN) then ! para usar a fIsica nAo linearizada (ClAudio Silva)
       !*********************************************                 
       !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY, 
       !                  -----------------------------------------------------
       if (KCALL.eq.1 ) then
          ISUM=0
          do JL=KIDIA,KFDIA
             if(LDFLAG(JL)) then
                ZQP(JL)=1./PSP(JL)
                ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                ZQSAT=min(0.5,ZQSAT)
                ZCOR=1./(1.-RETV *ZQSAT)
                ZQSAT=ZQSAT*ZCOR
                ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                ZCOND(JL)=max(ZCOND(JL),0.)
                PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND(JL)
                PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
                if(ZCOND(JL).ne.0.0) ISUM=ISUM+1
             else
                ZCOND(JL)=0.0
             endif
          enddo

          if(ISUM.eq.0) GO TO 230
          do JL=KIDIA,KFDIA
             if(LDFLAG(JL).and.ZCOND(JL).ne.0.) then
                ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                ZQSAT=min(0.5,ZQSAT)
                ZCOR=1./(1.-RETV*ZQSAT)
                ZQSAT=ZQSAT*ZCOR
                ZCOND1=(PQ(JL,KK)-ZQSAT) /(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
                PQ(JL,KK)=PQ(JL,KK)-ZCOND1
             endif
          enddo

230       continue

       endif

       if(KCALL.eq.2) then
          ISUM=0
          do JL=KIDIA,KFDIA
             if(LDFLAG(JL)) then
                ZQP(JL)=1./PSP(JL)
                ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                ZQSAT=min(0.5,ZQSAT)
                ZCOR=1./(1.-RETV  *ZQSAT)
                ZQSAT=ZQSAT*ZCOR
                ZCOND(JL)=(PQ(JL,KK)-ZQSAT) /(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                ZCOND(JL)=min(ZCOND(JL),0.)
                PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND(JL)
                PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
                if(ZCOND(JL).ne.0.0) ISUM=ISUM+1
             else
                ZCOND(JL)=0.0
             endif
          enddo

          if(ISUM.eq.0) GO TO 330

          do JL=KIDIA,KFDIA
             if(LDFLAG(JL).and.ZCOND(JL).ne.0.) then
                ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                ZQSAT=min(0.5,ZQSAT)
                ZCOR=1./(1.-RETV  *ZQSAT)
                ZQSAT=ZQSAT*ZCOR
                ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
                PQ(JL,KK)=PQ(JL,KK)-ZCOND1
             endif
          enddo

330       continue

       endif

       if(KCALL.eq.0) then
          do JL=KIDIA,KFDIA
             ZQP(JL)=1./PSP(JL)
             ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
             ZQSAT=min(0.5,ZQSAT)
             ZCOR=1./(1.-RETV  *ZQSAT)
             ZQSAT=ZQSAT*ZCOR
             ZCOND1=(PQ(JL,KK)-ZQSAT) /(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
             ZCOND1=max(0.,ZCOND1)
             PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
             PQ(JL,KK)=PQ(JL,KK)-ZCOND1
          enddo

          do JL=KIDIA,KFDIA
             ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
             ZQSAT=min(0.5,ZQSAT)
             ZCOR=1./(1.-RETV  *ZQSAT)
             ZQSAT=ZQSAT*ZCOR
             ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
             PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
             PQ(JL,KK)=PQ(JL,KK)-ZCOND1
          enddo

       endif
       if(KCALL.eq.4) then
          do JL=KIDIA,KFDIA
             ZQP(JL)=1./PSP(JL)
             ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
             ZQSAT=min(0.5,ZQSAT)
             ZCOR=1./(1.-RETV  *ZQSAT)
             ZQSAT=ZQSAT*ZCOR
             ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
             PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND(JL)
             PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
          enddo

          do JL=KIDIA,KFDIA
             ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
             ZQSAT=min(0.5,ZQSAT)
             ZCOR=1./(1.-RETV  *ZQSAT)
             ZQSAT=ZQSAT*ZCOR
             ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
             PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
             PQ(JL,KK)=PQ(JL,KK)-ZCOND1
          enddo

       endif

       !*********************************************
    else ! line 103
       !*********************************************                 
       !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
       !                  -----------------------------------------------------

       if (KCALL.eq.1 ) then

          ISUM=0
          do JL=KIDIA,KFDIA
             if(LDFLAG(JL)) then
                ZQP(JL)=1./PSP(JL)
                ZTARG=PT(JL,KK)
                ZOEALFA=0.5*(tanh(RLPAL1*(ZTARG-RLPTRC))+1.)
                ZFOEEWL=R2ES*exp(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                ZFOEEWI=R2ES*exp(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                Z1S=tanh(RLPAL2*(ZQSAT-ZQMAX))
                ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                ZCOR=1./(1.-RETV  *ZQSAT)
                ZQSAT=ZQSAT*ZCOR
                Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                ZCOND(JL)=max(ZCOND(JL),0.)
                PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND(JL)
                PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
                if(ZCOND(JL).ne.0.0) ISUM=ISUM+1
             else
                ZCOND(JL)=0.0
             endif
          enddo

          if(ISUM.eq.0) GO TO 231


          do JL=KIDIA,KFDIA
             if(LDFLAG(JL).and.ZCOND(JL).ne.0.) then

                ZTARG=PT(JL,KK)
                ZOEALFA=0.5*(tanh(RLPAL1*(ZTARG-RLPTRC))+1.)
                ZFOEEWL=R2ES*exp(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                ZFOEEWI=R2ES*exp(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                Z1S=tanh(RLPAL2*(ZQSAT-ZQMAX))
                ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                ZCOR=1./(1.-RETV  *ZQSAT)
                ZQSAT=ZQSAT*ZCOR
                Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+ (1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
                PQ(JL,KK)=PQ(JL,KK)-ZCOND1
             endif
          enddo

231       continue

       endif

       if(KCALL.eq.2) then

          ISUM=0
          do JL=KIDIA,KFDIA
             if(LDFLAG(JL)) then
                ZQP(JL)=1./PSP(JL)
                ZTARG=PT(JL,KK)
                ZOEALFA=0.5*(tanh(RLPAL1*(ZTARG-RLPTRC))+1.)
                ZFOEEWL=R2ES*exp(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                ZFOEEWI=R2ES*exp(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                Z1S=tanh(RLPAL2*(ZQSAT-ZQMAX))
                ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                ZCOR=1./(1.-RETV  *ZQSAT)
                ZQSAT=ZQSAT*ZCOR
                Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                ZCOND(JL)=min(ZCOND(JL),0.)
                PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND(JL)
                PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
                if(ZCOND(JL).ne.0.0) ISUM=ISUM+1
             else
                ZCOND(JL)=0.0
             endif
          enddo

          if(ISUM.eq.0) GO TO 331

          do JL=KIDIA,KFDIA
             if(LDFLAG(JL).and.ZCOND(JL).ne.0.) then
                ZTARG=PT(JL,KK)
                ZOEALFA=0.5*(tanh(RLPAL1*(ZTARG-RLPTRC))+1.)
                ZFOEEWL=R2ES*exp(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                ZFOEEWI=R2ES*exp(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                Z1S=tanh(RLPAL2*(ZQSAT-ZQMAX))
                ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                ZCOR=1./(1.-RETV  *ZQSAT)
                ZQSAT=ZQSAT*ZCOR
                Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
                PQ(JL,KK)=PQ(JL,KK)-ZCOND1
             endif
          enddo

331       continue

       endif

       if(KCALL.eq.0) then

          do JL=KIDIA,KFDIA
             ZQP(JL)=1./PSP(JL)
             ZTARG=PT(JL,KK)
             ZOEALFA=0.5*(tanh(RLPAL1*(ZTARG-RLPTRC))+1.)
             ZFOEEWL=R2ES*exp(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
             ZFOEEWI=R2ES*exp(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
             ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
             Z1S=tanh(RLPAL2*(ZQSAT-ZQMAX))
             ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
             ZCOR=1./(1.-RETV  *ZQSAT)
             ZQSAT=ZQSAT*ZCOR
             Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
             ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
             PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
             PQ(JL,KK)=PQ(JL,KK)-ZCOND1
             PT(JL,KK)=PT(JL,KK)
          enddo


          do JL=KIDIA,KFDIA
             ZTARG=PT(JL,KK)
             ZOEALFA=0.5*(tanh(RLPAL1*(ZTARG-RLPTRC))+1.)
             ZFOEEWL=R2ES*exp(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
             ZFOEEWI=R2ES*exp(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
             ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
             Z1S=tanh(RLPAL2*(ZQSAT-ZQMAX))
             ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
             ZCOR=1./(1.-RETV  *ZQSAT)
             ZQSAT=ZQSAT*ZCOR
             Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
             ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
             PT(JL,KK)=PT(JL,KK)+ (ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
             PQ(JL,KK)=PQ(JL,KK)-ZCOND1
          enddo


       endif

       if(KCALL.eq.4) then

          do JL=KIDIA,KFDIA
             ZQP(JL)=1./PSP(JL)
             ZTARG=PT(JL,KK)
             ZOEALFA=0.5*(tanh(RLPAL1*(ZTARG-RLPTRC))+1.)
             ZFOEEWL=R2ES*exp(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
             ZFOEEWI=R2ES*exp(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
             ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
             Z1S=tanh(RLPAL2*(ZQSAT-ZQMAX))
             ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
             ZCOR=1./(1.-RETV  *ZQSAT)
             ZQSAT=ZQSAT*ZCOR
             Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
             ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
             PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
          enddo

          do JL=KIDIA,KFDIA
             ZTARG=PT(JL,KK)
             ZOEALFA=0.5*(tanh(RLPAL1*(ZTARG-RLPTRC))+1.)
             ZFOEEWL=R2ES*exp(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
             ZFOEEWI=R2ES*exp(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
             ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
             Z1S=tanh(RLPAL2*(ZQSAT-ZQMAX))
             ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
             ZQSAT=min(ZQMAX,ZQSAT)
             ZCOR=1./(1.-RETV  *ZQSAT)
             ZQSAT=ZQSAT*ZCOR
             Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
             ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
             PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
             PQ(JL,KK)=PQ(JL,KK)-ZCOND1
          enddo
       endif

       !*********************************************
    endif
    !*********************************************          


  end subroutine CUADJTQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function FOEALFA (PTARG) 
    implicit none
    real  PTARG
    real, parameter :: RTICE = 273.16-0.1
    real, parameter :: RTWAT = 273.16 
    FOEALFA = min(1.,((max(RTICE,min(RTWAT,PTARG))-RTICE)/(RTWAT-RTICE))**2) 
    return
  end function FOEALFA


  real function FOEDELTA (PTARG) 
    implicit none
    real  PTARG
    real, parameter :: RTT = 273.16
    FOEDELTA = max (0.,sign(1.,PTARG-RTT))
    return
  end function FOEDELTA


  real function FOEEWM (PTARG)
    implicit none
    real  PTARG
    real, parameter :: RTT     = 273.16
    real, parameter :: RESTT   = 611.14
    real, parameter :: RD      = 287.04
    real, parameter :: RV      = 461.50
    real, parameter :: R2ES    = RESTT*(RD/RV)
    real, parameter :: R3LES   = 17.269 
    real, parameter :: R3IES   = 21.875
    real, parameter :: R4LES   = 35.86 
    real, parameter :: R4IES   = 7.66
    FOEEWM =(R2ES*(FOEALFA(PTARG)*exp(R3LES*(PTARG-RTT)/(PTARG-R4LES))+(1.-FOEALFA(PTARG))*exp(R3IES*(PTARG-RTT)/(PTARG-R4IES))))
    return
  end function FOEEWM


  real function FOEALFCU (PTARG) 
    implicit none
    real  PTARG
    real, parameter :: RTWAT = 273.16 
    real, parameter :: RTT = 273.16
    real, parameter :: RTICECU = RTT-23.0
    real, parameter :: RTWAT_RTICECU_R  = 1./(RTWAT-RTICECU)
    FOEALFCU = min(1.,((max(RTICECU,min(RTWAT,PTARG))-RTICECU)*RTWAT_RTICECU_R)**2)
    return
  end function FOEALFCU


  real function FOEDEM (PTARG)
    implicit none
    real  PTARG
    real, parameter :: R3LES   = 17.269 
    real, parameter :: R3IES   = 21.875
    real, parameter :: R4LES   = 35.86 
    real, parameter :: R4IES   = 7.66
    real, parameter :: R5LES   = R3LES*(273.16-R4LES)
    real, parameter :: R5IES   = R3IES*(273.16-R4IES)
    real, parameter :: RD      = 287.04
    real, parameter :: RV      = 461.50
    real, parameter :: RCPD    = 3.5*RD
    real, parameter :: RLVTT   = 2.5008E+6
    real, parameter :: RLSTT   = 2.8345E+6
    real, parameter :: R5ALVCP = R5LES*RLVTT/RCPD
    real, parameter :: R5ALSCP = R5IES*RLSTT/RCPD 
    real, parameter :: RALVDCP = RLVTT/RCPD
    FOEDEM = FOEALFA(PTARG)*R5ALVCP*(1./(PTARG-R4LES)**2)+(1.-FOEALFA(PTARG))*R5ALSCP*(1./(PTARG-R4IES)**2)
    return
  end function FOEDEM


  real function FOELDCPM(PTARG)
    implicit none
    real  PTARG
    real, parameter :: RD      = 287.04
    real, parameter :: RCPD    = 3.5*RD
    real, parameter :: RLSTT   = 2.8345E+6
    real, parameter :: RLVTT   = 2.5008E+6
    real, parameter :: RALSDCP = RLSTT/RCPD
    real, parameter :: RALVDCP = RLVTT/RCPD
    FOELDCPM = FOEALFA(PTARG)*RALVDCP+(1.-FOEALFA(PTARG))*RALSDCP
    return
  end function FOELDCPM

end module kbcon_ecmwf



