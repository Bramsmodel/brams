!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModMicControl

  use ModNamelistFile, only: &
       namelistFile

  use ModParallelEnvironment, only: &
       MsgDump

  use grid_dims, only: &
       NZPMAX, maxgrds

  implicit none

  include "files.h"
  include "MicConstants.h"

  private

  public :: MicControl
  public :: CreateMicControl
  public :: DestroyMicControl
  public :: DumpMicControl
  
  type MicControl
     !for rams 2M microphysics
     integer :: irime
     integer :: iplaws
     integer :: idust
     integer :: isalt
     integer :: imbudget
     integer :: imbudtot
     integer :: iccnlev
     integer :: imd1flg
     integer :: imd2flg 
     !--------------------------------------------------------------------------
     integer :: mcphys_type ! from RAMSIN
     integer :: level ! from RAMSIN
     integer :: icloud ! from RAMSIN
     integer :: idriz ! from RAMSIN
     integer :: irain ! from RAMSIN
     integer :: ipris ! from RAMSIN
     integer :: isnow ! from RAMSIN
     integer :: iaggr ! from RAMSIN
     integer :: igraup ! from RAMSIN
     integer :: ihail ! from RAMSIN
     integer :: mkcoltab ! from RAMSIN
     !---for micro-2m
     !integer   :: jnmb(ncat)
     integer :: jnmb(8)
     integer :: ipairc(nhcat,nhcat)
     integer :: ipairr(nhcat,nhcat)
     integer :: jhabtab(31,100,2)
     integer :: jhcat(nzpmax,ncat)
     integer :: ict1(nzpmax,ncat)
     integer :: ict2(nzpmax,ncat)
     real :: cparm ! from RAMSIN
     real :: rparm ! from RAMSIN
     real :: pparm ! from RAMSIN
     real :: sparm ! from RAMSIN
     real :: aparm ! from RAMSIN
     real :: gparm ! from RAMSIN
     real :: hparm ! from RAMSIN
     real :: dparm ! from RAMSIN
     real :: cnparm! from RAMSIN
     real :: gnparm! from RAMSIN
     real :: rictmin
     real :: rictmax
     real :: dps
     real :: dps2
     real :: d1min
     real :: r2min
     real :: d2min
     real :: d1max
     real :: r2max
     real :: d2max
     real :: d1ecc
     real :: d1ecr
     real :: r2ecr
     real :: r2err
     real :: colf
     real :: pi4dt
     real :: sedtime0
     real :: sedtime1
     real :: epsil
     real :: emb0(ncat)
     real :: emb1(ncat)
     !---for micro-2m
     real :: gnu(8) ! from RAMSIN
     real :: parm(ncat)
     real :: emb0log(ncat)
     real :: emb1log(ncat)
     real :: dict(ncat)
     real :: var_shape(nhcat)
     real :: cfmas(nhcat)
     real :: pwmas(nhcat)
     real :: cfvt(nhcat)
     real :: pwvt(nhcat)
     real :: dpsmi(nhcat)
     real :: cfden(nhcat)
     real :: pwden(nhcat)
     real :: cfemb0(nhcat)
     real :: cfen0(nhcat)
     real :: pwemb0(nhcat)
     real :: pwen0(nhcat)
     real :: vtfac(nhcat)
     real :: frefac1(nhcat)
     real :: frefac2(nhcat)
     real :: cfmasft(nhcat)
     real :: dnfac(nhcat)
     real :: sipfac(nhcat)
     real :: pwmasi(nhcat)
     real :: ch1(nhcat)
     real :: ch3(nhcat)
     real :: cdp1(nhcat)
     real :: pwvtmasi(nhcat)
     real :: xcoll(nzpmax)
     real :: rsedim(nzpmax)
     real :: tair(nzpmax)
     real :: tairc(nzpmax)
     real :: tairstrc(nzpmax)
     real :: til(nzpmax)
     real :: rvstr(nzpmax)
     real :: press(nzpmax)
     real :: pitot(nzpmax)
     real :: rliq(nzpmax)
     real :: rice(nzpmax)
     real :: qhydm(nzpmax)
     real :: rvlsair(nzpmax)
     real :: rvisair(nzpmax)
     real :: rvs0(nzpmax)
     real :: thrmcon(nzpmax)
     real :: vapdif(nzpmax)
     real :: dynvisc(nzpmax)
     real :: rdynvsci(nzpmax)
     real :: denfac(nzpmax)
     real :: dn0i(nzpmax)
     real :: colfacr(nzpmax)
     real :: colfacr2(nzpmax)
     real :: colfacc(nzpmax)
     real :: colfacc2(nzpmax)
     real :: sumuy(nzpmax)
     real :: sumuz(nzpmax)
     real :: sumvr(nzpmax)
     real :: scrmic1(nzpmax)
     real :: scrmic2(nzpmax)
     real :: scrmic3(nzpmax)
     real :: cccnx(nzpmax)
     real :: cifnx(nzpmax)
     real :: rx(nzpmax,ncat)
     real :: cx(nzpmax,ncat)
     real :: qr(nzpmax,ncat)
     real :: qx(nzpmax,ncat)
     real :: tx(nzpmax,ncat)
     real :: emb(nzpmax,ncat)
     real :: vterm(nzpmax,ncat)
     real :: vap(nzpmax,ncat)
     real :: ttest(nzpmax,ncat)
     real :: wct1(nzpmax,ncat)
     real :: wct2(nzpmax,ncat)
     real :: sb(nzpmax,ncat)
     real :: sd(nzpmax,ncat)
     real :: se(nzpmax,ncat)
     real :: sf(nzpmax,ncat)
     real :: sg(nzpmax,ncat)
     real :: sh(nzpmax,ncat)
     real :: sm(nzpmax,ncat)
     real :: ss(nzpmax,ncat)
     real :: su(nzpmax,ncat)
     real :: sw(nzpmax,ncat)
     real :: sy(nzpmax,ncat)
     real :: sz(nzpmax,ncat)
     real :: tref(nzpmax,2)
     real :: rvsref(nzpmax,2)
     real :: rvsrefp(nzpmax,2)
     real :: sa(nzpmax,9)
     real :: eff(nzpmax,10)
     real :: rxfer(nzpmax,ncat,ncat)
     real :: qrxfer(nzpmax,ncat,ncat)
     real :: enxfer(nzpmax,ncat,ncat)
     real :: dispemb0(nhcat,maxgrds)
     real :: dispemb1(nhcat,maxgrds)
     real :: ch2(nhcat,maxgrds)
     real :: coltabc(nembc,nembc,npairc)
     real :: coltabr(nembc,nembc,npairr)
     real :: frachz(nrhhz,nthz)
     real :: fracc(ndnc,ntc,maxgrds)
     real :: gamm(4)
     real :: gamn1(4)
     real :: gam(ngam,3)
     real :: gaminc(ngam,2)
     real :: gamsip13(ngam)
     real :: gamsip24(ngam)
     real :: rmlttab(ninc)
     real :: enmlttab(ninc,nhcat)
     real :: shedtab(ninc,ndns)
     real :: sc(2)
     real :: sk(2)
     real :: sl(2)
     real :: sj(7)
     real :: pcprx(7)
     real :: accpx(7)
     real :: r1tabcc(nd1cc)
     real :: c1tabcc(nd1cc)
     real :: c2tabcc(nd1cc)
     real :: r1tabcr(nd1cr,nr2cr,nd2cr)
     real :: c1tabcr(nd1cr,nr2cr,nd2cr)
     real :: c2tabrr(nr2rr,nd2rr)
     character(len=f_name_length) :: coltabfn ! from RAMSIN
  end type MicControl


  !---------------------------------------------------------------------------

contains



  function CreateMicControl(oneNamelistFile) result(res)
    ! Allocate a MicControl variable
    ! Initialize components originated by namelist file
    ! Remaining components should be initialized elsewhere
    type(namelistFile), pointer, intent(in) :: oneNamelistFile
    type(MicControl), pointer :: res

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateMicControl)**"

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if
    
    res%mcphys_type = oneNamelistFile%mcphys_type
    res%aparm = oneNamelistFile%aparm
    res%coltabfn = oneNamelistFile%coltabfn
    res%cparm = oneNamelistFile%cparm
    if (res%mcphys_type==0) then
       res%gnu(1:7) = oneNamelistFile%gnu(1:7)
    else if (res%mcphys_type==1) then
       res%gnu(1:8) = oneNamelistFile%gnu(1:8)
    end if
    res%gparm = oneNamelistFile%gparm
    res%hparm = oneNamelistFile%hparm
    res%iaggr = oneNamelistFile%iaggr
    res%icloud = oneNamelistFile%icloud
    res%idriz = oneNamelistFile%idriz
    res%igraup = oneNamelistFile%igraup
    res%ihail = oneNamelistFile%ihail
    res%ipris = oneNamelistFile%ipris
    res%irain = oneNamelistFile%irain
    res%isnow = oneNamelistFile%isnow
    res%level = oneNamelistFile%level
    res%mkcoltab = oneNamelistFile%mkcoltab
    res%pparm = oneNamelistFile%pparm
    res%rparm = oneNamelistFile%rparm
    res%sparm = oneNamelistFile%sparm

    !-for rams 2M microphysics
    !-special settings - needed further checking/moving to RAMSIN file

    res%dparm = oneNamelistFile%dparm 
    res%cnparm = oneNamelistFile%cnparm
    res%epsil   =oneNamelistFile%epsil
    res%gnparm = oneNamelistFile%gnparm
    res%irime   =oneNamelistFile%irime
    res%iplaws  =oneNamelistFile%iplaws
    res%iccnlev =oneNamelistFile%iccnlev

    !- local definition

    res%idust   =0! Dust: 0=dustmodeloff, 1=dustmodelon
    res%isalt   =0! Sea Salt model: 0 = off, 1 = on
    res%imbudget=0! Micro budget instant: 0=Off,1=partial,2=all
    res%imbudtot=0! Micro budget totals:  0=Off,1=partial,2=all
    res%imd1flg =0! Mineral Dust (small mode)
    res%imd2flg =0! Mineral Dust (large mode)
    !-
  end function CreateMicControl




  subroutine DestroyMicControl(oneMicControl)
    type(MicControl), pointer, intent(inout) :: oneMicControl

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyMicControl)**"

    if (associated(oneMicControl)) then
       deallocate(oneMicControl, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneMicControl fails with stat="//&
               trim(adjustl(str(1))))
       end if
       nullify(oneMicControl)
    end if
  end subroutine DestroyMicControl




  subroutine DumpMicControl(oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpMicControl)**"

    if (associated(oneMicControl)) then
       call MsgDump(h//" oneMicControl allocated")
    else
       call MsgDump(h//" oneMicControl deallocated")
    end if
  end subroutine DumpMicControl
end module ModMicControl
