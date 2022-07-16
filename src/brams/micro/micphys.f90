!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module micphys

  use ModMicControl, only: MicControl

  use ModNamelistFile, only: namelistFile

  use grid_dims, only: &
       NZPMAX, maxgrds

  implicit none

  include "files.h"

  character(len=32) :: lastCopyTo, lastCopyFrom

  !for rams 2M microphysics
  integer :: irime   ,&
       iplaws  ,&
       idust   ,&
       isalt   ,&
       imbudget,&
       imbudtot,&
       iccnlev ,&
       imd1flg ,&
       imd2flg 
  !--------------------------------------------------------------------------
  !     The product [(nthz-1)  * dthz ] must equal 25.0.
  !     The product [(nrhhz-1) * drhhz] must equal 0.18.
  !     The product [(ntc-1)   * dtc  ] must equal 20.0.
  !     The product [(ndnc-1)  * ddnc ] must equal 20.e-6.

  integer, parameter :: nthz=26
  integer, parameter :: nrhhz=10
  integer, parameter :: ngam=5000
  integer, parameter :: ninc=201
  integer, parameter :: ndns=15
  integer, parameter :: ntc=21
  integer, parameter :: ndnc=11
  integer, parameter :: nd1cc=30
  integer, parameter :: nd1cr=15
  integer, parameter :: nr2cr=10
  integer, parameter :: nd2cr=30
  integer, parameter :: nr2rr=20
  integer, parameter :: nd2rr=20
  integer, parameter :: ncat=7
  integer, parameter :: nhcat=15
  integer, parameter :: npairc=93
  integer, parameter :: npairr=131
  integer, parameter :: nembc=20

  real, parameter    :: dtc=1.
  real, parameter    :: ddnc=2.e-6
  real, parameter    :: dthz=1.
  real, parameter    :: drhhz=.02
  !--------------------------------------------------------------------------
  integer            :: mcphys_type ! from RAMSIN
  integer            :: level ! from RAMSIN
  integer            :: icloud ! from RAMSIN
  integer            :: idriz ! from RAMSIN
  integer            :: irain ! from RAMSIN
  integer            :: ipris ! from RAMSIN
  integer            :: isnow ! from RAMSIN
  integer            :: iaggr ! from RAMSIN
  integer            :: igraup ! from RAMSIN
  integer            :: ihail ! from RAMSIN
  integer            :: mkcoltab ! from RAMSIN
  !---for micro-2m
  !integer           :: jnmb(ncat)
  integer            :: jnmb(8)

  integer            :: ipairc(nhcat,nhcat)
  integer            :: ipairr(nhcat,nhcat)
  integer            :: jhabtab(31,100,2)
  integer            :: jhcat(nzpmax,ncat)
  integer            :: ict1(nzpmax,ncat)
  integer            :: ict2(nzpmax,ncat)

  real               :: cparm ! from RAMSIN
  real               :: rparm ! from RAMSIN
  real               :: pparm ! from RAMSIN
  real               :: sparm ! from RAMSIN
  real               :: aparm ! from RAMSIN
  real               :: gparm ! from RAMSIN
  real               :: hparm ! from RAMSIN
  real               :: dparm ! from RAMSIN
  real               :: cnparm! from RAMSIN
  real               :: gnparm! from RAMSIN

  real               :: rictmin
  real               :: rictmax
  real               :: dps
  real               :: dps2
  real               :: d1min
  real               :: r2min
  real               :: d2min
  real               :: d1max
  real               :: r2max
  real               :: d2max
  real               :: d1ecc
  real               :: d1ecr
  real               :: r2ecr
  real               :: r2err
  real               :: colf
  real               :: pi4dt
  real               :: sedtime0
  real               :: sedtime1
  real               :: epsil

  real               :: emb0(ncat)
  real               :: emb1(ncat)
  !---for micro-2m
  !real              :: gnu(ncat) ! from RAMSIN
  real               :: gnu(8) ! from RAMSIN
  !---
  real               :: parm(ncat)
  real               :: emb0log(ncat)
  real               :: emb1log(ncat)
  real               :: dict(ncat)


  ! Changing name of variabel SHAPE to VAR_SHAPE
  ! To avoid confusing with SHAPE FUNCTION (intrinsic)
  ! ALF

  real               :: var_shape(nhcat)
  real               :: cfmas(nhcat)
  real               :: pwmas(nhcat)
  real               :: cfvt(nhcat)
  real               :: pwvt(nhcat)
  real               :: dpsmi(nhcat)
  real               :: cfden(nhcat)
  real               :: pwden(nhcat)
  real               :: cfemb0(nhcat)
  real               :: cfen0(nhcat)
  real               :: pwemb0(nhcat)
  real               :: pwen0(nhcat)
  real               :: vtfac(nhcat)
  real               :: frefac1(nhcat)
  real               :: frefac2(nhcat)
  real               :: cfmasft(nhcat)
  real               :: dnfac(nhcat)
  real               :: sipfac(nhcat)
  real               :: pwmasi(nhcat)
  real               :: ch1(nhcat)
  real               :: ch3(nhcat)
  real               :: cdp1(nhcat)
  real               :: pwvtmasi(nhcat)

  !--(DMK-CARRIO-INI)-----------------------------------------------------
  !change_MP for chem
  real        :: xcoll(nzpmax)
  real        :: rsedim(nzpmax)
  !end change_MP
  !--(DMK-CARRIO-FIM)-----------------------------------------------------

  real         :: tair(nzpmax)
  real         :: tairc(nzpmax)
  real         :: tairstrc(nzpmax)
  real         :: til(nzpmax)
  real         :: rvstr(nzpmax)
  real         :: press(nzpmax)
  real         :: pitot(nzpmax)
  real         :: rliq(nzpmax)
  real         :: rice(nzpmax)
  real         :: qhydm(nzpmax)
  real         :: rvlsair(nzpmax)
  real         :: rvisair(nzpmax)
  real         :: rvs0(nzpmax)
  real         :: thrmcon(nzpmax)
  real         :: vapdif(nzpmax)
  real         :: dynvisc(nzpmax)
  real         :: rdynvsci(nzpmax)
  real         :: denfac(nzpmax)
  real         :: dn0i(nzpmax)
  real         :: colfacr(nzpmax)
  real         :: colfacr2(nzpmax)
  real         :: colfacc(nzpmax)
  real         :: colfacc2(nzpmax)
  real         :: sumuy(nzpmax)
  real         :: sumuz(nzpmax)
  real         :: sumvr(nzpmax)
  real         :: scrmic1(nzpmax)
  real         :: scrmic2(nzpmax)
  real         :: scrmic3(nzpmax)
  real         :: cccnx(nzpmax)
  real         :: cifnx(nzpmax)

  real         :: rx(nzpmax,ncat)
  real         :: cx(nzpmax,ncat)
  real         :: qr(nzpmax,ncat)
  real         :: qx(nzpmax,ncat)
  real         :: tx(nzpmax,ncat)
  real         :: emb(nzpmax,ncat)
  real         :: vterm(nzpmax,ncat)
  real         :: vap(nzpmax,ncat)
  real         :: ttest(nzpmax,ncat)
  real         :: wct1(nzpmax,ncat)
  real         :: wct2(nzpmax,ncat)
  real         :: sb(nzpmax,ncat)
  real         :: sd(nzpmax,ncat)
  real         :: se(nzpmax,ncat)
  real         :: sf(nzpmax,ncat)
  real         :: sg(nzpmax,ncat)
  real         :: sh(nzpmax,ncat)
  real         :: sm(nzpmax,ncat)
  real         :: ss(nzpmax,ncat)
  real         :: su(nzpmax,ncat)
  real         :: sw(nzpmax,ncat)
  real         :: sy(nzpmax,ncat)
  real         :: sz(nzpmax,ncat)
  real         :: tref(nzpmax,2)
  real         :: rvsref(nzpmax,2)
  real         :: rvsrefp(nzpmax,2)
  real         :: sa(nzpmax,9)
  real         :: eff(nzpmax,10)

  real         :: rxfer(nzpmax,ncat,ncat)
  real         :: qrxfer(nzpmax,ncat,ncat)
  real         :: enxfer(nzpmax,ncat,ncat)
  real         :: dispemb0(nhcat,maxgrds)
  real         :: dispemb1(nhcat,maxgrds)
  real         :: ch2(nhcat,maxgrds)

  real         :: coltabc(nembc,nembc,npairc)
  real         :: coltabr(nembc,nembc,npairr)

  real         :: frachz(nrhhz,nthz)
  real         :: fracc(ndnc,ntc,maxgrds)
  real         :: gamm(4)
  real         :: gamn1(4)
  real         :: gam(ngam,3)
  real         :: gaminc(ngam,2)
  real         :: gamsip13(ngam)
  real         :: gamsip24(ngam)
  real         :: rmlttab(ninc)
  real         :: enmlttab(ninc,nhcat)
  real         :: shedtab(ninc,ndns)
  real         :: sc(2)
  real         :: sk(2)
  real         :: sl(2)
  real         :: sj(7)
  real         :: pcprx(7)
  real         :: accpx(7)

  real         :: r1tabcc(nd1cc)
  real         :: c1tabcc(nd1cc)
  real         :: c2tabcc(nd1cc)
  real         :: r1tabcr(nd1cr,nr2cr,nd2cr)
  real         :: c1tabcr(nd1cr,nr2cr,nd2cr)
  real         :: c2tabrr(nr2rr,nd2rr)

  character(len=f_name_length) :: coltabfn ! from RAMSIN
  ! Modif. by ALF

  !---------------------------------------------------------------------------

contains
  subroutine StoreNamelistFileAtMicphys(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    mcphys_type = oneNamelistFile%mcphys_type
    aparm = oneNamelistFile%aparm
    coltabfn = oneNamelistFile%coltabfn
    cparm = oneNamelistFile%cparm
    if    (mcphys_type==0) then
       gnu(1:7) = oneNamelistFile%gnu(1:7)
    elseif(mcphys_type==1) then
       gnu(1:8) = oneNamelistFile%gnu(1:8)
    endif
    gparm = oneNamelistFile%gparm
    hparm = oneNamelistFile%hparm
    iaggr = oneNamelistFile%iaggr
    icloud = oneNamelistFile%icloud
    idriz = oneNamelistFile%idriz
    igraup = oneNamelistFile%igraup
    ihail = oneNamelistFile%ihail
    ipris = oneNamelistFile%ipris
    irain = oneNamelistFile%irain
    isnow = oneNamelistFile%isnow
    level = oneNamelistFile%level
    mkcoltab = oneNamelistFile%mkcoltab
    pparm = oneNamelistFile%pparm
    rparm = oneNamelistFile%rparm
    sparm = oneNamelistFile%sparm

    !-for rams 2M microphysics
    !-special settings - needed further checking/moving to RAMSIN file
    dparm = oneNamelistFile%dparm 
    cnparm = oneNamelistFile%cnparm
    epsil   =oneNamelistFile%epsil
    gnparm = oneNamelistFile%gnparm
    irime   =oneNamelistFile%irime
    iplaws  =oneNamelistFile%iplaws
    iccnlev =oneNamelistFile%iccnlev
    !- local definition
    idust   =0! Dust: 0=dustmodeloff, 1=dustmodelon
    isalt   =0! Sea Salt model: 0 = off, 1 = on
    imbudget=0! Micro budget instant: 0=Off,1=partial,2=all
    imbudtot=0! Micro budget totals:  0=Off,1=partial,2=all
    imd1flg =0! Mineral Dust (small mode)
    imd2flg =0! Mineral Dust (large mode)
    !-
  end subroutine StoreNamelistFileAtMicphys

  subroutine DeepCopyFromMicControl(oneMicControl, proc)
    type(MicControl), pointer, intent(in) :: oneMicControl
    character(len=*), intent(in) :: proc

    character(len=*), parameter :: h="**(DeepCopyFromMicControl)**"

    if (lastCopyTo == "") then
       call fatal_error(h//" wrong order; last CopyTo by "//lastCopyTo)
    end if

    lastCopyTo=""
    lastCopyFrom=proc


    irime=oneMicControl%irime
    iplaws=oneMicControl%iplaws
    idust=oneMicControl%idust
    isalt=oneMicControl%isalt
    imbudget=oneMicControl%imbudget
    imbudtot=oneMicControl%imbudtot
    iccnlev=oneMicControl%iccnlev
    imd1flg=oneMicControl%imd1flg
    imd2flg =oneMicControl%imd2flg 
    mcphys_type =oneMicControl%mcphys_type 
    level =oneMicControl%level 
    icloud =oneMicControl%icloud 
    idriz =oneMicControl%idriz 
    irain =oneMicControl%irain 
    ipris =oneMicControl%ipris 
    isnow =oneMicControl%isnow 
    iaggr =oneMicControl%iaggr 
    igraup =oneMicControl%igraup 
    ihail =oneMicControl%ihail 
    mkcoltab =oneMicControl%mkcoltab 
    jnmb=oneMicControl%jnmb
    jnmb=oneMicControl%jnmb
    ipairc=oneMicControl%ipairc
    ipairr=oneMicControl%ipairr
    jhabtab=oneMicControl%jhabtab
    jhcat=oneMicControl%jhcat
    ict1=oneMicControl%ict1
    ict2=oneMicControl%ict2
    cparm =oneMicControl%cparm 
    rparm =oneMicControl%rparm 
    pparm =oneMicControl%pparm 
    sparm =oneMicControl%sparm 
    aparm =oneMicControl%aparm 
    gparm =oneMicControl%gparm 
    hparm =oneMicControl%hparm 
    dparm =oneMicControl%dparm 
    cnparm=oneMicControl%cnparm
    gnparm=oneMicControl%gnparm
    rictmin=oneMicControl%rictmin
    rictmax=oneMicControl%rictmax
    dps=oneMicControl%dps
    dps2=oneMicControl%dps2
    d1min=oneMicControl%d1min
    r2min=oneMicControl%r2min
    d2min=oneMicControl%d2min
    d1max=oneMicControl%d1max
    r2max=oneMicControl%r2max
    d2max=oneMicControl%d2max
    d1ecc=oneMicControl%d1ecc
    d1ecr=oneMicControl%d1ecr
    r2ecr=oneMicControl%r2ecr
    r2err=oneMicControl%r2err
    colf=oneMicControl%colf
    pi4dt=oneMicControl%pi4dt
    sedtime0=oneMicControl%sedtime0
    sedtime1=oneMicControl%sedtime1
    epsil=oneMicControl%epsil
    emb0=oneMicControl%emb0
    emb1=oneMicControl%emb1
    gnu=oneMicControl%gnu
    parm=oneMicControl%parm
    emb0log=oneMicControl%emb0log
    emb1log=oneMicControl%emb1log
    dict=oneMicControl%dict
    var_shape=oneMicControl%var_shape
    cfmas=oneMicControl%cfmas
    pwmas=oneMicControl%pwmas
    cfvt=oneMicControl%cfvt
    pwvt=oneMicControl%pwvt
    dpsmi=oneMicControl%dpsmi
    cfden=oneMicControl%cfden
    pwden=oneMicControl%pwden
    cfemb0=oneMicControl%cfemb0
    cfen0=oneMicControl%cfen0
    pwemb0=oneMicControl%pwemb0
    pwen0=oneMicControl%pwen0
    vtfac=oneMicControl%vtfac
    frefac1=oneMicControl%frefac1
    frefac2=oneMicControl%frefac2
    cfmasft=oneMicControl%cfmasft
    dnfac=oneMicControl%dnfac
    sipfac=oneMicControl%sipfac
    pwmasi=oneMicControl%pwmasi
    ch1=oneMicControl%ch1
    ch3=oneMicControl%ch3
    cdp1=oneMicControl%cdp1
    pwvtmasi=oneMicControl%pwvtmasi
    xcoll=oneMicControl%xcoll
    rsedim=oneMicControl%rsedim
    tair=oneMicControl%tair
    tairc=oneMicControl%tairc
    tairstrc=oneMicControl%tairstrc
    til=oneMicControl%til
    rvstr=oneMicControl%rvstr
    press=oneMicControl%press
    pitot=oneMicControl%pitot
    rliq=oneMicControl%rliq
    rice=oneMicControl%rice
    qhydm=oneMicControl%qhydm
    rvlsair=oneMicControl%rvlsair
    rvisair=oneMicControl%rvisair
    rvs0=oneMicControl%rvs0
    thrmcon=oneMicControl%thrmcon
    vapdif=oneMicControl%vapdif
    dynvisc=oneMicControl%dynvisc
    rdynvsci=oneMicControl%rdynvsci
    denfac=oneMicControl%denfac
    dn0i=oneMicControl%dn0i
    colfacr=oneMicControl%colfacr
    colfacr2=oneMicControl%colfacr2
    colfacc=oneMicControl%colfacc
    colfacc2=oneMicControl%colfacc2
    sumuy=oneMicControl%sumuy
    sumuz=oneMicControl%sumuz
    sumvr=oneMicControl%sumvr
    scrmic1=oneMicControl%scrmic1
    scrmic2=oneMicControl%scrmic2
    scrmic3=oneMicControl%scrmic3
    cccnx=oneMicControl%cccnx
    cifnx=oneMicControl%cifnx
    rx=oneMicControl%rx
    cx=oneMicControl%cx
    qr=oneMicControl%qr
    qx=oneMicControl%qx
    tx=oneMicControl%tx
    emb=oneMicControl%emb
    vterm=oneMicControl%vterm
    vap=oneMicControl%vap
    ttest=oneMicControl%ttest
    wct1=oneMicControl%wct1
    wct2=oneMicControl%wct2
    sb=oneMicControl%sb
    sd=oneMicControl%sd
    se=oneMicControl%se
    sf=oneMicControl%sf
    sg=oneMicControl%sg
    sh=oneMicControl%sh
    sm=oneMicControl%sm
    ss=oneMicControl%ss
    su=oneMicControl%su
    sw=oneMicControl%sw
    sy=oneMicControl%sy
    sz=oneMicControl%sz
    tref=oneMicControl%tref
    rvsref=oneMicControl%rvsref
    rvsrefp=oneMicControl%rvsrefp
    sa=oneMicControl%sa
    eff=oneMicControl%eff
    rxfer=oneMicControl%rxfer
    qrxfer=oneMicControl%qrxfer
    enxfer=oneMicControl%enxfer
    dispemb0=oneMicControl%dispemb0
    dispemb1=oneMicControl%dispemb1
    ch2=oneMicControl%ch2
    coltabc=oneMicControl%coltabc
    coltabr=oneMicControl%coltabr
    frachz=oneMicControl%frachz
    fracc=oneMicControl%fracc
    gamm=oneMicControl%gamm
    gamn1=oneMicControl%gamn1
    gam=oneMicControl%gam
    gaminc=oneMicControl%gaminc
    gamsip13=oneMicControl%gamsip13
    gamsip24=oneMicControl%gamsip24
    rmlttab=oneMicControl%rmlttab
    enmlttab=oneMicControl%enmlttab
    shedtab=oneMicControl%shedtab
    sc=oneMicControl%sc
    sk=oneMicControl%sk
    sl=oneMicControl%sl
    sj=oneMicControl%sj
    pcprx=oneMicControl%pcprx
    accpx=oneMicControl%accpx
    r1tabcc=oneMicControl%r1tabcc
    c1tabcc=oneMicControl%c1tabcc
    c2tabcc=oneMicControl%c2tabcc
    r1tabcr=oneMicControl%r1tabcr
    c1tabcr=oneMicControl%c1tabcr
    c2tabrr=oneMicControl%c2tabrr
    coltabfn=oneMicControl%coltabfn
  end subroutine DeepCopyFromMicControl



  subroutine DeepCopyToMicControl(oneMicControl, proc)
    type(MicControl), pointer, intent(in) :: oneMicControl
    character(len=*), intent(in) :: proc

    character(len=*), parameter :: h="**(DeepCopyToMicControl)**"

    if (lastCopyFrom == "") then
       call fatal_error(h//" wrong order; last CopyFrom by "//lastCopyFrom)
    end if

    lastCopyTo=proc
    lastCopyFrom=""

    oneMicControl%irime=irime
    oneMicControl%iplaws=iplaws
    oneMicControl%idust=idust
    oneMicControl%isalt=isalt
    oneMicControl%imbudget=imbudget
    oneMicControl%imbudtot=imbudtot
    oneMicControl%iccnlev=iccnlev
    oneMicControl%imd1flg=imd1flg
    oneMicControl%imd2flg =imd2flg 
    oneMicControl%mcphys_type =mcphys_type 
    oneMicControl%level =level 
    oneMicControl%icloud =icloud 
    oneMicControl%idriz =idriz 
    oneMicControl%irain =irain 
    oneMicControl%ipris =ipris 
    oneMicControl%isnow =isnow 
    oneMicControl%iaggr =iaggr 
    oneMicControl%igraup =igraup 
    oneMicControl%ihail =ihail 
    oneMicControl%mkcoltab =mkcoltab 
    oneMicControl%jnmb=jnmb
    oneMicControl%jnmb=jnmb
    oneMicControl%ipairc=ipairc
    oneMicControl%ipairr=ipairr
    oneMicControl%jhabtab=jhabtab
    oneMicControl%jhcat=jhcat
    oneMicControl%ict1=ict1
    oneMicControl%ict2=ict2
    oneMicControl%cparm =cparm 
    oneMicControl%rparm =rparm 
    oneMicControl%pparm =pparm 
    oneMicControl%sparm =sparm 
    oneMicControl%aparm =aparm 
    oneMicControl%gparm =gparm 
    oneMicControl%hparm =hparm 
    oneMicControl%dparm =dparm 
    oneMicControl%cnparm=cnparm
    oneMicControl%gnparm=gnparm
    oneMicControl%rictmin=rictmin
    oneMicControl%rictmax=rictmax
    oneMicControl%dps=dps
    oneMicControl%dps2=dps2
    oneMicControl%d1min=d1min
    oneMicControl%r2min=r2min
    oneMicControl%d2min=d2min
    oneMicControl%d1max=d1max
    oneMicControl%r2max=r2max
    oneMicControl%d2max=d2max
    oneMicControl%d1ecc=d1ecc
    oneMicControl%d1ecr=d1ecr
    oneMicControl%r2ecr=r2ecr
    oneMicControl%r2err=r2err
    oneMicControl%colf=colf
    oneMicControl%pi4dt=pi4dt
    oneMicControl%sedtime0=sedtime0
    oneMicControl%sedtime1=sedtime1
    oneMicControl%epsil=epsil
    oneMicControl%emb0=emb0
    oneMicControl%emb1=emb1
    oneMicControl%gnu=gnu
    oneMicControl%parm=parm
    oneMicControl%emb0log=emb0log
    oneMicControl%emb1log=emb1log
    oneMicControl%dict=dict
    oneMicControl%var_shape=var_shape
    oneMicControl%cfmas=cfmas
    oneMicControl%pwmas=pwmas
    oneMicControl%cfvt=cfvt
    oneMicControl%pwvt=pwvt
    oneMicControl%dpsmi=dpsmi
    oneMicControl%cfden=cfden
    oneMicControl%pwden=pwden
    oneMicControl%cfemb0=cfemb0
    oneMicControl%cfen0=cfen0
    oneMicControl%pwemb0=pwemb0
    oneMicControl%pwen0=pwen0
    oneMicControl%vtfac=vtfac
    oneMicControl%frefac1=frefac1
    oneMicControl%frefac2=frefac2
    oneMicControl%cfmasft=cfmasft
    oneMicControl%dnfac=dnfac
    oneMicControl%sipfac=sipfac
    oneMicControl%pwmasi=pwmasi
    oneMicControl%ch1=ch1
    oneMicControl%ch3=ch3
    oneMicControl%cdp1=cdp1
    oneMicControl%pwvtmasi=pwvtmasi
    oneMicControl%xcoll=xcoll
    oneMicControl%rsedim=rsedim
    oneMicControl%tair=tair
    oneMicControl%tairc=tairc
    oneMicControl%tairstrc=tairstrc
    oneMicControl%til=til
    oneMicControl%rvstr=rvstr
    oneMicControl%press=press
    oneMicControl%pitot=pitot
    oneMicControl%rliq=rliq
    oneMicControl%rice=rice
    oneMicControl%qhydm=qhydm
    oneMicControl%rvlsair=rvlsair
    oneMicControl%rvisair=rvisair
    oneMicControl%rvs0=rvs0
    oneMicControl%thrmcon=thrmcon
    oneMicControl%vapdif=vapdif
    oneMicControl%dynvisc=dynvisc
    oneMicControl%rdynvsci=rdynvsci
    oneMicControl%denfac=denfac
    oneMicControl%dn0i=dn0i
    oneMicControl%colfacr=colfacr
    oneMicControl%colfacr2=colfacr2
    oneMicControl%colfacc=colfacc
    oneMicControl%colfacc2=colfacc2
    oneMicControl%sumuy=sumuy
    oneMicControl%sumuz=sumuz
    oneMicControl%sumvr=sumvr
    oneMicControl%scrmic1=scrmic1
    oneMicControl%scrmic2=scrmic2
    oneMicControl%scrmic3=scrmic3
    oneMicControl%cccnx=cccnx
    oneMicControl%cifnx=cifnx
    oneMicControl%rx=rx
    oneMicControl%cx=cx
    oneMicControl%qr=qr
    oneMicControl%qx=qx
    oneMicControl%tx=tx
    oneMicControl%emb=emb
    oneMicControl%vterm=vterm
    oneMicControl%vap=vap
    oneMicControl%ttest=ttest
    oneMicControl%wct1=wct1
    oneMicControl%wct2=wct2
    oneMicControl%sb=sb
    oneMicControl%sd=sd
    oneMicControl%se=se
    oneMicControl%sf=sf
    oneMicControl%sg=sg
    oneMicControl%sh=sh
    oneMicControl%sm=sm
    oneMicControl%ss=ss
    oneMicControl%su=su
    oneMicControl%sw=sw
    oneMicControl%sy=sy
    oneMicControl%sz=sz
    oneMicControl%tref=tref
    oneMicControl%rvsref=rvsref
    oneMicControl%rvsrefp=rvsrefp
    oneMicControl%sa=sa
    oneMicControl%eff=eff
    oneMicControl%rxfer=rxfer
    oneMicControl%qrxfer=qrxfer
    oneMicControl%enxfer=enxfer
    oneMicControl%dispemb0=dispemb0
    oneMicControl%dispemb1=dispemb1
    oneMicControl%ch2=ch2
    oneMicControl%coltabc=coltabc
    oneMicControl%coltabr=coltabr
    oneMicControl%frachz=frachz
    oneMicControl%fracc=fracc
    oneMicControl%gamm=gamm
    oneMicControl%gamn1=gamn1
    oneMicControl%gam=gam
    oneMicControl%gaminc=gaminc
    oneMicControl%gamsip13=gamsip13
    oneMicControl%gamsip24=gamsip24
    oneMicControl%rmlttab=rmlttab
    oneMicControl%enmlttab=enmlttab
    oneMicControl%shedtab=shedtab
    oneMicControl%sc=sc
    oneMicControl%sk=sk
    oneMicControl%sl=sl
    oneMicControl%sj=sj
    oneMicControl%pcprx=pcprx
    oneMicControl%accpx=accpx
    oneMicControl%r1tabcc=r1tabcc
    oneMicControl%c1tabcc=c1tabcc
    oneMicControl%c2tabcc=c2tabcc
    oneMicControl%r1tabcr=r1tabcr
    oneMicControl%c1tabcr=c1tabcr
    oneMicControl%c2tabrr=c2tabrr
    oneMicControl%coltabfn=coltabfn
  end subroutine DeepCopyToMicControl
end Module micphys
