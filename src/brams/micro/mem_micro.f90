!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################


module mem_micro

  use mem_radiate, only: &
       ilwrtyp, &
       iswrtyp       ! INTENT(IN)

  use mem_cuparm , only: &
       nnqparm                ! INTENT(IN)

  use ModVarTables, only: &
       InsertVTab

  use ModMicControl, only: &
       MicControl

  use ModMicroFields, only: &
       MicroFields

  implicit none

  include "constants.h"

  character(len=32) :: copyTo=""
  character(len=32) :: copyFrom=""
  private

  public :: micro_vars
  public :: micro_g
  public :: microm_g
  public :: alloc_micro
  public :: nullify_micro
  public :: dealloc_micro
  public :: filltab_micro
  public :: DeepCopyToMicroFields
  public :: DeepCopyFromMicroFields

  type micro_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, contiguous, pointer :: rcp(:,:,:)
     real, contiguous, pointer :: rdp(:,:,:)
     real, contiguous, pointer :: rrp(:,:,:)
     real, contiguous, pointer :: rpp(:,:,:)
     real, contiguous, pointer :: rsp(:,:,:)
     real, contiguous, pointer :: rap(:,:,:)
     real, contiguous, pointer :: rgp(:,:,:)
     real, contiguous, pointer :: rhp(:,:,:)
     real, contiguous, pointer :: ccp(:,:,:)
     real, contiguous, pointer :: cdp(:,:,:)
     real, contiguous, pointer :: crp(:,:,:)
     real, contiguous, pointer :: cpp(:,:,:)
     real, contiguous, pointer :: csp(:,:,:)
     real, contiguous, pointer :: cap(:,:,:)
     real, contiguous, pointer :: cgp(:,:,:)
     real, contiguous, pointer :: chp(:,:,:)
     real, contiguous, pointer :: cccnp(:,:,:)
     real, contiguous, pointer :: gccnp(:,:,:)
     real, contiguous, pointer :: cifnp(:,:,:)
     real, contiguous, pointer :: q2(:,:,:)
     real, contiguous, pointer :: q6(:,:,:)
     real, contiguous, pointer :: q7(:,:,:)
     real, contiguous, pointer :: rei(:,:,:)
     real, contiguous, pointer :: rel(:,:,:)
     real, contiguous, pointer :: cldfr(:,:,:)
     real, contiguous, pointer :: cccmp(:,:,:)
     real, contiguous, pointer :: gccmp(:,:,:)
     real, contiguous, pointer :: cnm1p(:,:,:)
     real, contiguous, pointer :: cnm2p(:,:,:)
     real, contiguous, pointer :: cnm3p(:,:,:)
     real, contiguous, pointer :: cnm8p(:,:,:)
     real, contiguous, pointer :: md1np(:,:,:)
     real, contiguous, pointer :: md2np(:,:,:)
     real, contiguous, pointer :: salt_filmp(:,:,:)
     real, contiguous, pointer :: salt_jetp(:,:,:)
     real, contiguous, pointer :: salt_spmp(:,:,:)
     real, contiguous, pointer :: pcpvr(:,:,:)
     real, contiguous, pointer :: pcpvp(:,:,:)
     real, contiguous, pointer :: pcpvs(:,:,:)
     real, contiguous, pointer :: pcpva(:,:,:)
     real, contiguous, pointer :: pcpvg(:,:,:)
     real, contiguous, pointer :: pcpvh(:,:,:)
     real, contiguous, pointer :: pcpvd(:,:,:)
     real, contiguous, pointer :: nuccldr(:,:,:)
     real, contiguous, pointer :: nuccldc(:,:,:)
     real, contiguous, pointer :: nucicer(:,:,:)
     real, contiguous, pointer :: nucicec(:,:,:)
     real, contiguous, pointer :: inuchomr(:,:,:)
     real, contiguous, pointer :: inuchomc(:,:,:)
     real, contiguous, pointer :: inuccontr(:,:,:)
     real, contiguous, pointer :: inuccontc(:,:,:)
     real, contiguous, pointer :: inucifnr(:,:,:)
     real, contiguous, pointer :: inucifnc(:,:,:)
     real, contiguous, pointer :: inuchazr(:,:,:)
     real, contiguous, pointer :: inuchazc(:,:,:)
     real, contiguous, pointer :: vapliq(:,:,:)
     real, contiguous, pointer :: vapice(:,:,:)
     real, contiguous, pointer :: vapcld(:,:,:)
     real, contiguous, pointer :: vaprain(:,:,:)
     real, contiguous, pointer :: vappris(:,:,:)
     real, contiguous, pointer :: vapsnow(:,:,:)
     real, contiguous, pointer :: vapaggr(:,:,:)
     real, contiguous, pointer :: vapgrau(:,:,:)
     real, contiguous, pointer :: vaphail(:,:,:)
     real, contiguous, pointer :: vapdriz(:,:,:)
     real, contiguous, pointer :: meltice(:,:,:)
     real, contiguous, pointer :: meltpris(:,:,:)
     real, contiguous, pointer :: meltsnow(:,:,:)
     real, contiguous, pointer :: meltaggr(:,:,:)
     real, contiguous, pointer :: meltgrau(:,:,:)
     real, contiguous, pointer :: melthail(:,:,:)
     real, contiguous, pointer :: cld2rain(:,:,:)
     real, contiguous, pointer :: rimecld(:,:,:)
     real, contiguous, pointer :: rimecldsnow(:,:,:)
     real, contiguous, pointer :: rimecldaggr(:,:,:)
     real, contiguous, pointer :: rimecldgrau(:,:,:)
     real, contiguous, pointer :: rimecldhail(:,:,:)
     real, contiguous, pointer :: rain2ice(:,:,:)
     real, contiguous, pointer :: rain2pr(:,:,:)
     real, contiguous, pointer :: rain2sn(:,:,:)
     real, contiguous, pointer :: rain2ag(:,:,:)
     real, contiguous, pointer :: rain2gr(:,:,:)
     real, contiguous, pointer :: rain2ha(:,:,:)
     real, contiguous, pointer :: rain2ha_xtra(:,:,:)
     real, contiguous, pointer :: ice2rain(:,:,:)
     real, contiguous, pointer :: aggregate(:,:,:)
     real, contiguous, pointer :: aggrselfpris(:,:,:)
     real, contiguous, pointer :: aggrselfsnow(:,:,:)
     real, contiguous, pointer :: aggrprissnow(:,:,:)
     real, contiguous, pointer :: latheatvap(:,:,:)
     real, contiguous, pointer :: latheatfrz(:,:,:)
     real, contiguous, pointer :: nuccldrt(:,:,:)
     real, contiguous, pointer :: nuccldct(:,:,:)
     real, contiguous, pointer :: nucicert(:,:,:)
     real, contiguous, pointer :: nucicect(:,:,:)
     real, contiguous, pointer :: inuchomrt(:,:,:)
     real, contiguous, pointer :: inuchomct(:,:,:)
     real, contiguous, pointer :: inuccontrt(:,:,:)
     real, contiguous, pointer :: inuccontct(:,:,:)
     real, contiguous, pointer :: inucifnrt(:,:,:)
     real, contiguous, pointer :: inucifnct(:,:,:)
     real, contiguous, pointer :: inuchazrt(:,:,:)
     real, contiguous, pointer :: inuchazct(:,:,:)
     real, contiguous, pointer :: vapliqt(:,:,:)
     real, contiguous, pointer :: vapicet(:,:,:)
     real, contiguous, pointer :: vapcldt(:,:,:)
     real, contiguous, pointer :: vapraint(:,:,:)
     real, contiguous, pointer :: vapprist(:,:,:)
     real, contiguous, pointer :: vapsnowt(:,:,:)
     real, contiguous, pointer :: vapaggrt(:,:,:)
     real, contiguous, pointer :: vapgraut(:,:,:)
     real, contiguous, pointer :: vaphailt(:,:,:)
     real, contiguous, pointer :: vapdrizt(:,:,:)
     real, contiguous, pointer :: melticet(:,:,:)
     real, contiguous, pointer :: meltprist(:,:,:)
     real, contiguous, pointer :: meltsnowt(:,:,:)
     real, contiguous, pointer :: meltaggrt(:,:,:)
     real, contiguous, pointer :: meltgraut(:,:,:)
     real, contiguous, pointer :: melthailt(:,:,:)
     real, contiguous, pointer :: cld2raint(:,:,:)
     real, contiguous, pointer :: rimecldt(:,:,:)
     real, contiguous, pointer :: rimecldsnowt(:,:,:)
     real, contiguous, pointer :: rimecldaggrt(:,:,:)
     real, contiguous, pointer :: rimecldgraut(:,:,:)
     real, contiguous, pointer :: rimecldhailt(:,:,:)
     real, contiguous, pointer :: rain2icet(:,:,:)
     real, contiguous, pointer :: rain2prt(:,:,:)
     real, contiguous, pointer :: rain2snt(:,:,:)
     real, contiguous, pointer :: rain2agt(:,:,:)
     real, contiguous, pointer :: rain2grt(:,:,:)
     real, contiguous, pointer :: rain2hat(:,:,:)
     real, contiguous, pointer :: rain2ha_xtrat(:,:,:)
     real, contiguous, pointer :: ice2raint(:,:,:)
     real, contiguous, pointer :: aggregatet(:,:,:)
     real, contiguous, pointer :: aggrselfprist(:,:,:)
     real, contiguous, pointer :: aggrselfsnowt(:,:,:)
     real, contiguous, pointer :: aggrprissnowt(:,:,:)
     real, contiguous, pointer :: latheatvapt(:,:,:)
     real, contiguous, pointer :: latheatfrzt(:,:,:)

     ! Variables to be dimensioned by (nnxp,nyp)
     real, contiguous, pointer :: accpr(:,:)
     real, contiguous, pointer :: accpp(:,:)
     real, contiguous, pointer :: accps(:,:)
     real, contiguous, pointer :: accpa(:,:)
     real, contiguous, pointer :: accpg(:,:)
     real, contiguous, pointer :: accph(:,:)
     real, contiguous, pointer :: accpd(:,:)
     real, contiguous, pointer :: pcprr(:,:)
     real, contiguous, pointer :: pcprp(:,:)
     real, contiguous, pointer :: pcprs(:,:)
     real, contiguous, pointer :: pcpra(:,:)
     real, contiguous, pointer :: pcprg(:,:)
     real, contiguous, pointer :: pcprh(:,:)
     real, contiguous, pointer :: pcprd(:,:)
     real, contiguous, pointer :: pcpg(:,:)
     real, contiguous, pointer :: qpcpg(:,:)
     real, contiguous, pointer :: dpcpg(:,:)

  end type micro_vars

  type (micro_vars), allocatable :: micro_g(:), microm_g(:)


  interface ToCopy
     module procedure ToCopy_2D
     module procedure ToCopy_3D
  end interface ToCopy

contains                  

  subroutine ToCopy_2D(src, dest)
    real, pointer, intent(in) :: src(:,:)
    real, pointer, intent(in) :: dest(:,:)
    dest=src
  end subroutine ToCopy_2D


  subroutine ToCopy_3D(src, dest)
    real, pointer, intent(in) :: src(:,:,:)
    real, pointer, intent(in) :: dest(:,:,:)
    dest=src
  end subroutine ToCopy_3D


  subroutine alloc_micro(micro,n1,n2,n3,ng,oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    type (micro_vars) :: micro
    integer, intent(in) :: n1,n2,n3,ng

    ! Allocate arrays based on options (if necessary)


    if(oneMicControl%mcphys_type == 2 .or. &
         oneMicControl%mcphys_type == 3 .or. &
         oneMicControl%mcphys_type == 4) then ! gthompson microphysics

       oneMicControl%level = 3
       oneMicControl%idriz = 0
       oneMicControl%icloud = 1
       oneMicControl%irain = 1
       oneMicControl%ipris  = 1 
       oneMicControl%isnow = 1
       oneMicControl%igraup = 1
       oneMicControl%ihail = 0
       oneMicControl%iaggr  = 0
       oneMicControl%jnmb(1) = 1 !cloud
       oneMicControl%jnmb(2) = 1 !rain
       oneMicControl%jnmb(3) = 1 !pristine
       oneMicControl%jnmb(4) = 1 !snow
       oneMicControl%jnmb(5) = 0 !agg
       oneMicControl%jnmb(6) = 1 !graupel
       oneMicControl%jnmb(7) = 0 !ihail
       oneMicControl%jnmb(8) = 0 !idriz

       !- cloud liq water
       allocate (micro%rcp(n1,n2,n3))
       micro%rcp  =0.0

       !- rain
       allocate (micro%rrp  (n1,n2,n3))
       micro%rrp  =0.0
       !- for this scheme, the rain rate below will
       !- account for rain+ice+snow+graupel
       allocate (micro%accpr(n2,n3))
       micro%accpr=0.0
       allocate (micro%pcprr(n2,n3))
       micro%pcprr=0.0

       !- ice
       allocate (micro%rpp  (n1,n2,n3))
       micro%rpp  =0.0
       !- don t need to be allocated, see coments above
       !allocate (micro%accpp(n2,n3)) ;   micro%accpp=0.0
       !allocate (micro%pcprp(n2,n3)) ;   micro%pcprp=0.0

       !- snow
       allocate (micro%rsp  (n1,n2,n3))
       micro%rsp=0.0
       !- the rates bellow will account for snow and ice
       allocate (micro%accps(n2,n3))
       micro%accps =0.0
       allocate (micro%pcprs(n2,n3))
       micro%pcprs =0.0

       !- graupel
       allocate (micro%rgp  (n1,n2,n3))
       micro%rgp  =0.0
       !- the rates bellow will account for only graupel
       allocate (micro%accpg(n2,n3))
       micro%accpg=0.0
       allocate (micro%pcprg(n2,n3))
       micro%pcprg=0.0

       if(oneMicControl%mcphys_type  == 2 .or. &
            oneMicControl%mcphys_type  == 3) then ! only for double-moment and 
          !- number concentration for cloud/rain/ice
          !- obs : ccp don t need to be allocated for the single-moment
          !- cloud water scheme (the same for CCN and IFN).
          allocate(micro%crp  (n1,n2,n3))
          micro%crp  =0.0 
          allocate(micro%cpp  (n1,n2,n3))
          micro%cpp  =0.0 
          !---these should not be allocated for mcphys_type  == 2 because
          !---they are not used for this option
          !            !ST
          !          allocate(micro%ccp  (n1,n2,n3)) ;micro%ccp  =0.0 
          !            allocate(micro%cccnp(n1,n2,n3)) ;micro%cccnp=0.0 !;endif 
          !            allocate(micro%cifnp(n1,n2,n3)) ;micro%cifnp=0.0 !;endif 
          !            !ST
       endif
       !- only for cloud water double-moment and aerosol aware microphysics         
       if(oneMicControl%mcphys_type  == 3) then ! only for double-moment and 
          allocate(micro%ccp  (n1,n2,n3))
          micro%ccp  =0.0 
          allocate(micro%cccnp(n1,n2,n3))
          micro%cccnp=0.0 !;endif 
          allocate(micro%cifnp(n1,n2,n3))
          micro%cifnp=0.0 !;endif 
       endif

       !- 3D cloud fraction from GFDL cloud microphysics and GF convection
       if(oneMicControl%mcphys_type  == 4 .or. nnqparm(ng) == 8) then 
          allocate(micro%cldfr  (n1,n2,n3))
          micro%cldfr  =0.0 
       endif

       !- for consistency with the other parts of BRAMS
       !- pgcp will be the total precipitation rate
       allocate (micro%pcpg (n2,n3))
       micro%pcpg =0.0
       !-the allocations below are tmp for leaf-3
       allocate (micro%qpcpg(n2,n3))
       micro%qpcpg=0.0         
       allocate (micro%dpcpg(n2,n3))
       micro%dpcpg=0.0

       !- allocation of memory for effective radius for RRTMG
       if(ilwrtyp==6 .or. iswrtyp==6 ) then
          allocate (micro%rei  (n1,n2,n3))
          micro%rei  =0.0  
          allocate (micro%rel  (n1,n2,n3))
          micro%rel  =0.0
       endif

    else  ! for the traditional RAMS microphysics

       if (oneMicControl%level >= 2 ) then
          allocate (micro%rcp(n1,n2,n3))
          micro%rcp    =0.0
       endif
       if (oneMicControl%level >= 3) then         
          if(oneMicControl%irain >= 1)  then
             allocate (micro%rrp  (n1,n2,n3))
             micro%rrp  =0.0
             allocate (micro%accpr(n2,n3))
             micro%accpr=0.0
             allocate (micro%pcprr(n2,n3))
             micro%pcprr=0.0
             allocate (micro%pcpvr(n1,n2,n3))
             micro%pcpvr=0.0
             allocate (micro%q2   (n1,n2,n3))
             micro%q2   =0.0
          endif
          if(oneMicControl%ipris >= 1)  then
             allocate (micro%rpp  (n1,n2,n3))
             micro%rpp  =0.0
             allocate (micro%accpp(n2,n3))
             micro%accpp=0.0
             allocate (micro%pcprp(n2,n3))
             micro%pcprp=0.0
             allocate (micro%pcpvp(n1,n2,n3))
             micro%pcpvp=0.0
          endif
          if(oneMicControl%isnow >= 1)  then
             allocate (micro%rsp  (n1,n2,n3))
             micro%rsp   =0.0
             allocate (micro%accps(n2,n3))
             micro%accps =0.0
             allocate (micro%pcprs(n2,n3))
             micro%pcprs =0.0
             allocate (micro%pcpvs(n1,n2,n3))
             micro%pcpvs =0.0
          endif
          if(oneMicControl%iaggr >= 1)  then
             allocate (micro%rap  (n1,n2,n3))
             micro%rap  =0.0
             allocate (micro%accpa(n2,n3))
             micro%accpa=0.0
             allocate (micro%pcpra(n2,n3))
             micro%pcpra=0.0
             allocate (micro%pcpva(n1,n2,n3))
             micro%pcpva=0.0
          endif
          if(oneMicControl%igraup >= 1) then
             allocate (micro%rgp  (n1,n2,n3))
             micro%rgp  =0.0
             allocate (micro%accpg(n2,n3))
             micro%accpg=0.0
             allocate (micro%pcprg(n2,n3))
             micro%pcprg=0.0
             allocate (micro%pcpvg(n1,n2,n3))
             micro%pcpvg=0.0
             allocate (micro%q6   (n1,n2,n3))
             micro%q6   =0.0
          endif
          if(oneMicControl%ihail >= 1)  then
             allocate (micro%rhp  (n1,n2,n3))
             micro%rhp  =0.0
             allocate (micro%accph(n2,n3))
             micro%accph=0.0
             allocate (micro%pcprh(n2,n3))
             micro%pcprh=0.0
             allocate (micro%pcpvh(n1,n2,n3))
             micro%pcpvh=0.0
             allocate (micro%q7   (n1,n2,n3))
             micro%q7   =0.0
          endif
          if(oneMicControl%jnmb(1) >= 5)  then
             allocate(micro%ccp  (n1,n2,n3))
             micro%ccp  =0.0
          endif
          if(oneMicControl%jnmb(2) == 5)  then
             allocate(micro%crp  (n1,n2,n3))
             micro%crp  =0.0
          endif
          if(oneMicControl%jnmb(3) >= 5)  then
             allocate(micro%cpp  (n1,n2,n3))
             micro%cpp  =0.0
          endif
          if(oneMicControl%jnmb(4) == 5)  then
             allocate(micro%csp  (n1,n2,n3))
             micro%csp  =0.0
          endif
          if(oneMicControl%jnmb(5) == 5)  then
             allocate(micro%cap  (n1,n2,n3))
             micro%cap  =0.0
          endif
          if(oneMicControl%jnmb(6) == 5)  then
             allocate(micro%cgp  (n1,n2,n3))
             micro%cgp  =0.0
          endif
          if(oneMicControl%jnmb(7) == 5)  then
             allocate(micro%chp  (n1,n2,n3))
             micro%chp  =0.0
          endif
          if(oneMicControl%icloud  >= 5)  then
             allocate(micro%cccnp(n1,n2,n3))
             micro%cccnp=0.0
          endif
          if(oneMicControl%ipris   >= 5)  then
             allocate(micro%cifnp(n1,n2,n3))
             micro%cifnp=0.0
          endif

          if(oneMicControl%icloud >= 5)   then
             allocate(micro%cccmp(n1,n2,n3))
             micro%cccmp=0.0
          endif

          allocate (micro%pcpg (n2,n3))
          micro%pcpg =0.0
          allocate (micro%qpcpg(n2,n3))
          micro%qpcpg=0.0
          allocate (micro%dpcpg(n2,n3))
          micro%dpcpg=0.0

          !- only for 2M microphysics
          if(oneMicControl%mcphys_type == 1)  then
             if(oneMicControl%idriz >= 1 )  then
                allocate (micro%rdp  (n1,n2,n3))
                micro%rdp  =0.0
                allocate (micro%accpd(n2,n3))
                micro%accpd=0.0
                allocate (micro%pcprd(n2,n3))
                micro%pcprd=0.0
                allocate (micro%pcpvd(n1,n2,n3))
                micro%pcpvd=0.0
             endif

             if(oneMicControl%jnmb(8) >= 5)  then
                allocate(micro%cdp  (n1,n2,n3))
                micro%cdp  =0.0
             endif
             if(oneMicControl%idriz>= 5)  then
                allocate(micro%gccnp(n1,n2,n3))
                micro%gccnp=0.0
             endif
             if(oneMicControl%idriz   >= 5)  then
                allocate(micro%gccmp(n1,n2,n3))
                micro%gccmp=0.0
             endif

             if(oneMicControl%iccnlev >= 2 .and. oneMicControl%jnmb(1) >= 5)  then
                allocate(micro%cnm1p(n1,n2,n3))
                micro%cnm1p=0.0
             endif
             if(oneMicControl%iccnlev >= 2 .and. oneMicControl%jnmb(2) >= 1)  then
                allocate(micro%cnm2p(n1,n2,n3))
                micro%cnm2p=0.0
             endif
             if(oneMicControl%iccnlev >= 2 .and. oneMicControl%jnmb(3) >= 1)  then
                allocate(micro%cnm3p(n1,n2,n3))
                micro%cnm3p=0.0
             endif
             if(oneMicControl%iccnlev >= 2 .and. oneMicControl%jnmb(8) >= 1)  then
                allocate(micro%cnm8p(n1,n2,n3))
                micro%cnm8p=0.0
             endif

             if(oneMicControl%idust == 1 .or. oneMicControl%imd1flg == 1)  then
                allocate(micro%md1np(n1,n2,n3))
                micro%md1np=0.0
             endif
             if(oneMicControl%idust == 1 .or. oneMicControl%imd2flg == 1)  then
                allocate(micro%md2np(n1,n2,n3))
                micro%md2np=0.0
             endif
             if(oneMicControl%isalt == 1) then
                allocate(micro%salt_filmp(n1,n2,n3))
                micro%salt_filmp =0.0
             endif
             if(oneMicControl%isalt == 1) then
                allocate(micro%salt_jetp (n1,n2,n3))
                micro%salt_jetp  =0.0
             endif
             if(oneMicControl%isalt == 1) then
                allocate(micro%salt_spmp (n1,n2,n3))
                micro%salt_spmp  =0.0
             endif


             !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
             if(oneMicControl%imbudget>=1 .or. oneMicControl%imbudtot>=1) then
                allocate (micro%latheatvap(n1,n2,n3))
                micro%latheatvap=0.0
                allocate (micro%latheatfrz(n1,n2,n3))
                micro%latheatfrz=0.0
             endif
             if(oneMicControl%imbudget>=1) then
                allocate (micro%nuccldr  (n1,n2,n3))
                micro%nuccldr  =0.0
                allocate (micro%nuccldc  (n1,n2,n3))
                micro%nuccldc  =0.0
                allocate (micro%cld2rain (n1,n2,n3))
                micro%cld2rain =0.0
                allocate (micro%ice2rain (n1,n2,n3))
                micro%ice2rain =0.0
                allocate (micro%nucicer  (n1,n2,n3))
                micro%nucicer  =0.0
                allocate (micro%nucicec  (n1,n2,n3))
                micro%nucicec  =0.0
                allocate (micro%vapliq(n1,n2,n3))
                micro%vapliq   =0.0   
                allocate (micro%vapice(n1,n2,n3))
                micro%vapice   =0.0
                allocate (micro%meltice  (n1,n2,n3))
                micro%meltice  =0.0
                allocate (micro%rimecld  (n1,n2,n3))
                micro%rimecld  =0.0
                allocate (micro%rain2ice (n1,n2,n3))
                micro%rain2ice =0.0
                allocate (micro%aggregate(n1,n2,n3))
                micro%aggregate=0.0
             endif
             if(oneMicControl%imbudget==2) then
                allocate (micro%inuchomr    (n1,n2,n3))
                micro%inuchomr    =0.0 
                allocate (micro%inuchomc    (n1,n2,n3))
                micro%inuchomc    =0.0
                allocate (micro%inuccontr   (n1,n2,n3))
                micro%inuccontr   =0.0
                allocate (micro%inuccontc   (n1,n2,n3))
                micro%inuccontc   =0.0
                allocate (micro%inucifnr    (n1,n2,n3))
                micro%inucifnr    =0.0
                allocate (micro%inucifnc    (n1,n2,n3))
                micro%inucifnc    =0.0
                allocate (micro%inuchazr    (n1,n2,n3))
                micro%inuchazr    =0.0
                allocate (micro%inuchazc    (n1,n2,n3))
                micro%inuchazc    =0.0
                allocate (micro%vapcld   (n1,n2,n3))
                micro%vapcld  =0.0
                allocate (micro%vaprain     (n1,n2,n3))
                micro%vaprain  =0.0
                allocate (micro%vappris     (n1,n2,n3))
                micro%vappris  =0.0
                allocate (micro%vapsnow     (n1,n2,n3))
                micro%vapsnow  =0.0   
                allocate (micro%vapaggr     (n1,n2,n3))
                micro%vapaggr  =0.0
                allocate (micro%vapgrau     (n1,n2,n3))
                micro%vapgrau  =0.0
                allocate (micro%vaphail     (n1,n2,n3))
                micro%vaphail  =0.0
                allocate (micro%vapdriz     (n1,n2,n3))
                micro%vapdriz  =0.0
                allocate (micro%meltpris    (n1,n2,n3))
                micro%meltpris    =0.0
                allocate (micro%meltsnow    (n1,n2,n3))
                micro%meltsnow    =0.0
                allocate (micro%meltaggr    (n1,n2,n3))
                micro%meltaggr    =0.0
                allocate (micro%meltgrau    (n1,n2,n3))
                micro%meltgrau    =0.0
                allocate (micro%melthail    (n1,n2,n3))
                micro%melthail    =0.0
                allocate (micro%rimecldsnow (n1,n2,n3))
                micro%rimecldsnow =0.0
                allocate (micro%rimecldaggr (n1,n2,n3))
                micro%rimecldaggr =0.0
                allocate (micro%rimecldgrau (n1,n2,n3))
                micro%rimecldgrau =0.0
                allocate (micro%rimecldhail (n1,n2,n3))
                micro%rimecldhail =0.0
                allocate (micro%rain2pr     (n1,n2,n3))
                micro%rain2pr  =0.0
                allocate (micro%rain2sn     (n1,n2,n3))
                micro%rain2sn  =0.0
                allocate (micro%rain2ag     (n1,n2,n3))
                micro%rain2ag  =0.0
                allocate (micro%rain2gr     (n1,n2,n3))
                micro%rain2gr  =0.0
                allocate (micro%rain2ha     (n1,n2,n3))
                micro%rain2ha  =0.0
                allocate (micro%rain2ha_xtra(n1,n2,n3))
                micro%rain2ha_xtra=0.0
                allocate (micro%aggrselfpris(n1,n2,n3))
                micro%aggrselfpris=0.0
                allocate (micro%aggrselfsnow(n1,n2,n3))
                micro%aggrselfsnow=0.0
                allocate (micro%aggrprissnow(n1,n2,n3))
                micro%aggrprissnow=0.0
             endif
             !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
             if(oneMicControl%imbudtot>=1) then
                allocate (micro%nuccldrt   (n1,n2,n3))
                micro%nuccldrt=0.0
                allocate (micro%nuccldct   (n1,n2,n3))
                micro%nuccldct=0.0
                allocate (micro%cld2raint  (n1,n2,n3))
                micro%cld2raint  =0.0
                allocate (micro%ice2raint  (n1,n2,n3))
                micro%ice2raint  =0.0
                allocate (micro%nucicert   (n1,n2,n3))
                micro%nucicert=0.0
                allocate (micro%nucicect   (n1,n2,n3))
                micro%nucicect=0.0
                allocate (micro%vapliqt    (n1,n2,n3))
                micro%vapliqt=0.0
                allocate (micro%vapicet    (n1,n2,n3))
                micro%vapicet=0.0
                allocate (micro%melticet   (n1,n2,n3))
                micro%melticet=0.0
                allocate (micro%rimecldt   (n1,n2,n3))
                micro%rimecldt=0.0
                allocate (micro%rain2icet  (n1,n2,n3))
                micro%rain2icet  =0.0
                allocate (micro%aggregatet (n1,n2,n3))
                micro%aggregatet =0.0
                allocate (micro%latheatvapt(n1,n2,n3))
                micro%latheatvapt=0.0
                allocate (micro%latheatfrzt(n1,n2,n3))
                micro%latheatfrzt=0.0
             endif
             if(oneMicControl%imbudtot==2) then
                allocate (micro%inuchomrt(n1,n2,n3))
                micro%inuchomrt    =0.0
                allocate (micro%inuchomct(n1,n2,n3))
                micro%inuchomct    =0.0
                allocate (micro%inuccontrt(n1,n2,n3))
                micro%inuccontrt   =0.0
                allocate (micro%inuccontct(n1,n2,n3))
                micro%inuccontct   =0.0
                allocate (micro%inucifnrt(n1,n2,n3))
                micro%inucifnrt    =0.0
                allocate (micro%inucifnct(n1,n2,n3))
                micro%inucifnct    =0.0
                allocate (micro%inuchazrt(n1,n2,n3))
                micro%inuchazrt    =0.0
                allocate (micro%inuchazct(n1,n2,n3))
                micro%inuchazct    =0.0
                allocate (micro%vapcldt(n1,n2,n3))
                micro%vapcldt      =0.0
                allocate (micro%vapraint(n1,n2,n3))
                micro%vapraint     =0.0
                allocate (micro%vapprist(n1,n2,n3))
                micro%vapprist     =0.0
                allocate (micro%vapsnowt(n1,n2,n3))
                micro%vapsnowt     =0.0
                allocate (micro%vapaggrt(n1,n2,n3))
                micro%vapaggrt     =0.0
                allocate (micro%vapgraut(n1,n2,n3))
                micro%vapgraut     =0.0
                allocate (micro%vaphailt(n1,n2,n3))
                micro%vaphailt     =0.0
                allocate (micro%vapdrizt(n1,n2,n3))
                micro%vapdrizt     =0.0
                allocate (micro%meltprist(n1,n2,n3))
                micro%meltprist    =0.0
                allocate (micro%meltsnowt(n1,n2,n3))
                micro%meltsnowt    =0.0
                allocate (micro%meltaggrt(n1,n2,n3))
                micro%meltaggrt    =0.0
                allocate (micro%meltgraut(n1,n2,n3))
                micro%meltgraut    =0.0
                allocate (micro%melthailt(n1,n2,n3))
                micro%melthailt    =0.0
                allocate (micro%rimecldsnowt(n1,n2,n3))
                micro%rimecldsnowt =0.0
                allocate (micro%rimecldaggrt(n1,n2,n3))
                micro%rimecldaggrt =0.0
                allocate (micro%rimecldgraut(n1,n2,n3))
                micro%rimecldgraut =0.0
                allocate (micro%rimecldhailt(n1,n2,n3))
                micro%rimecldhailt =0.0
                allocate (micro%rain2prt(n1,n2,n3))
                micro%rain2prt     =0.0
                allocate (micro%rain2snt(n1,n2,n3))
                micro%rain2snt     =0.0
                allocate (micro%rain2agt(n1,n2,n3))
                micro%rain2agt     =0.0
                allocate (micro%rain2grt(n1,n2,n3))
                micro%rain2grt     =0.0
                allocate (micro%rain2hat(n1,n2,n3))
                micro%rain2hat     =0.0
                allocate (micro%rain2ha_xtrat(n1,n2,n3))
                micro%rain2ha_xtrat=0.0
                allocate (micro%aggrselfprist(n1,n2,n3))
                micro%aggrselfprist=0.0
                allocate (micro%aggrselfsnowt(n1,n2,n3))
                micro%aggrselfsnowt=0.0
                allocate (micro%aggrprissnowt(n1,n2,n3))
                micro%aggrprissnowt=0.0
             endif
          endif! mcphys_type=1     
          !- allocation of memory for effective radius for RRTMG
          if(ilwrtyp==6 .or. iswrtyp==6 ) then
             allocate (micro%rei  (n1,n2,n3))
             micro%rei  =0.0  
             allocate (micro%rel  (n1,n2,n3))
             micro%rel  =0.0
          endif
       endif ! oneMicControl%level >=3 
    endif
    return 
  end subroutine alloc_micro




  subroutine nullify_micro(micro) 
    type (micro_vars) :: micro

    if (associated(micro%rcp))     nullify (micro%rcp)
    if (associated(micro%rdp))     nullify (micro%rdp)
    if (associated(micro%rrp))     nullify (micro%rrp)
    if (associated(micro%rpp))     nullify (micro%rpp)
    if (associated(micro%rsp))     nullify (micro%rsp)
    if (associated(micro%rap))     nullify (micro%rap)
    if (associated(micro%rgp))     nullify (micro%rgp)
    if (associated(micro%rhp))     nullify (micro%rhp)
    if (associated(micro%ccp))     nullify (micro%ccp)
    if (associated(micro%cdp))     nullify (micro%cdp)
    if (associated(micro%crp))     nullify (micro%crp)
    if (associated(micro%cpp))     nullify (micro%cpp)
    if (associated(micro%csp))     nullify (micro%csp)
    if (associated(micro%cap))     nullify (micro%cap)
    if (associated(micro%cgp))     nullify (micro%cgp)
    if (associated(micro%chp))     nullify (micro%chp)
    if (associated(micro%cccnp))   nullify (micro%cccnp)
    if (associated(micro%gccnp))   nullify (micro%gccnp)
    if (associated(micro%cifnp))   nullify (micro%cifnp)
    if (associated(micro%q2))      nullify (micro%q2)
    if (associated(micro%q6))      nullify (micro%q6)
    if (associated(micro%q7))      nullify (micro%q7)

    if (associated(micro%cccmp))   nullify (micro%cccmp)
    if (associated(micro%gccmp))   nullify (micro%gccmp)
    if (associated(micro%cnm1p))   nullify (micro%cnm1p)
    if (associated(micro%cnm2p))   nullify (micro%cnm2p)
    if (associated(micro%cnm3p))   nullify (micro%cnm3p)
    if (associated(micro%cnm8p))   nullify (micro%cnm8p)
    if (associated(micro%md1np))   nullify (micro%md1np)
    if (associated(micro%md2np))   nullify (micro%md2np)
    if (associated(micro%salt_filmp)) nullify (micro%salt_filmp)
    if (associated(micro%salt_jetp))  nullify (micro%salt_jetp)
    if (associated(micro%salt_spmp))  nullify (micro%salt_spmp)
    if (associated(micro%pcpvr))   nullify (micro%pcpvr)
    if (associated(micro%pcpvp))   nullify (micro%pcpvp)
    if (associated(micro%pcpvs))   nullify (micro%pcpvs)
    if (associated(micro%pcpva))   nullify (micro%pcpva)
    if (associated(micro%pcpvg))   nullify (micro%pcpvg)
    if (associated(micro%pcpvh))   nullify (micro%pcpvh)
    if (associated(micro%pcpvd))   nullify (micro%pcpvd)

    if (associated(micro%accpr))   nullify (micro%accpr)
    if (associated(micro%accpp))   nullify (micro%accpp)
    if (associated(micro%accps))   nullify (micro%accps)
    if (associated(micro%accpa))   nullify (micro%accpa)
    if (associated(micro%accpg))   nullify (micro%accpg)
    if (associated(micro%accph))   nullify (micro%accph)
    if (associated(micro%accpd))   nullify (micro%accpd)
    if (associated(micro%pcprr))   nullify (micro%pcprr)
    if (associated(micro%pcprp))   nullify (micro%pcprp)
    if (associated(micro%pcprs))   nullify (micro%pcprs)
    if (associated(micro%pcpra))   nullify (micro%pcpra)
    if (associated(micro%pcprg))   nullify (micro%pcprg)
    if (associated(micro%pcprh))   nullify (micro%pcprh)
    if (associated(micro%pcprd))   nullify (micro%pcprd)
    if (associated(micro%pcpg))    nullify (micro%pcpg)
    if (associated(micro%qpcpg))   nullify (micro%qpcpg)
    if (associated(micro%dpcpg))   nullify (micro%dpcpg)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
    if (associated(micro%nuccldr))      nullify (micro%nuccldr)
    if (associated(micro%nuccldc))      nullify (micro%nuccldc)
    if (associated(micro%nucicer))      nullify (micro%nucicer)
    if (associated(micro%nucicec))      nullify (micro%nucicec)
    if (associated(micro%inuchomr))     nullify (micro%inuchomr)
    if (associated(micro%inuchomc))     nullify (micro%inuchomc)
    if (associated(micro%inuccontr))    nullify (micro%inuccontr)
    if (associated(micro%inuccontc))    nullify (micro%inuccontc)
    if (associated(micro%inucifnr))     nullify (micro%inucifnr)
    if (associated(micro%inucifnc))     nullify (micro%inucifnc)
    if (associated(micro%inuchazr))     nullify (micro%inuchazr)
    if (associated(micro%inuchazc))     nullify (micro%inuchazc)
    if (associated(micro%vapliq))       nullify (micro%vapliq)
    if (associated(micro%vapice))       nullify (micro%vapice)
    if (associated(micro%vapcld))       nullify (micro%vapcld)
    if (associated(micro%vaprain))      nullify (micro%vaprain)
    if (associated(micro%vappris))      nullify (micro%vappris)
    if (associated(micro%vapsnow))      nullify (micro%vapsnow)
    if (associated(micro%vapaggr))      nullify (micro%vapaggr)
    if (associated(micro%vapgrau))      nullify (micro%vapgrau)
    if (associated(micro%vaphail))      nullify (micro%vaphail)
    if (associated(micro%vapdriz))      nullify (micro%vapdriz)
    if (associated(micro%meltice))      nullify (micro%meltice)
    if (associated(micro%meltpris))     nullify (micro%meltpris)
    if (associated(micro%meltsnow))     nullify (micro%meltsnow)
    if (associated(micro%meltaggr))     nullify (micro%meltaggr)
    if (associated(micro%meltgrau))     nullify (micro%meltgrau)
    if (associated(micro%melthail))     nullify (micro%melthail)
    if (associated(micro%cld2rain))     nullify (micro%cld2rain)
    if (associated(micro%rimecld))      nullify (micro%rimecld)
    if (associated(micro%rimecldsnow))  nullify (micro%rimecldsnow)
    if (associated(micro%rimecldaggr))  nullify (micro%rimecldaggr)
    if (associated(micro%rimecldgrau))  nullify (micro%rimecldgrau)
    if (associated(micro%rimecldhail))  nullify (micro%rimecldhail)
    if (associated(micro%rain2ice))     nullify (micro%rain2ice)
    if (associated(micro%rain2pr))      nullify (micro%rain2pr)
    if (associated(micro%rain2sn))      nullify (micro%rain2sn)
    if (associated(micro%rain2ag))      nullify (micro%rain2ag)
    if (associated(micro%rain2gr))      nullify (micro%rain2gr)
    if (associated(micro%rain2ha))      nullify (micro%rain2ha)
    if (associated(micro%rain2ha_xtra)) nullify (micro%rain2ha_xtra)
    if (associated(micro%ice2rain))     nullify (micro%ice2rain)
    if (associated(micro%aggregate))    nullify (micro%aggregate)
    if (associated(micro%aggrselfpris)) nullify (micro%aggrselfpris)
    if (associated(micro%aggrselfsnow)) nullify (micro%aggrselfsnow)
    if (associated(micro%aggrprissnow)) nullify (micro%aggrprissnow)
    if (associated(micro%latheatvap))   nullify (micro%latheatvap)
    if (associated(micro%latheatfrz))   nullify (micro%latheatfrz)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
    if (associated(micro%nuccldrt))      nullify (micro%nuccldrt)
    if (associated(micro%nuccldct))      nullify (micro%nuccldct)
    if (associated(micro%nucicert))      nullify (micro%nucicert)
    if (associated(micro%nucicect))      nullify (micro%nucicect)
    if (associated(micro%inuchomrt))     nullify (micro%inuchomrt)
    if (associated(micro%inuchomct))     nullify (micro%inuchomct)
    if (associated(micro%inuccontrt))    nullify (micro%inuccontrt)
    if (associated(micro%inuccontct))    nullify (micro%inuccontct)
    if (associated(micro%inucifnrt))     nullify (micro%inucifnrt)
    if (associated(micro%inucifnct))     nullify (micro%inucifnct)
    if (associated(micro%inuchazrt))     nullify (micro%inuchazrt)
    if (associated(micro%inuchazct))     nullify (micro%inuchazct)
    if (associated(micro%vapliqt))       nullify (micro%vapliqt)
    if (associated(micro%vapicet))       nullify (micro%vapicet)
    if (associated(micro%vapcldt))       nullify (micro%vapcldt)
    if (associated(micro%vapraint))      nullify (micro%vapraint)
    if (associated(micro%vapprist))      nullify (micro%vapprist)
    if (associated(micro%vapsnowt))      nullify (micro%vapsnowt)
    if (associated(micro%vapaggrt))      nullify (micro%vapaggrt)
    if (associated(micro%vapgraut))      nullify (micro%vapgraut)
    if (associated(micro%vaphailt))      nullify (micro%vaphailt)
    if (associated(micro%vapdrizt))      nullify (micro%vapdrizt)
    if (associated(micro%melticet))      nullify (micro%melticet)
    if (associated(micro%meltprist))     nullify (micro%meltprist)
    if (associated(micro%meltsnowt))     nullify (micro%meltsnowt)
    if (associated(micro%meltaggrt))     nullify (micro%meltaggrt)
    if (associated(micro%meltgraut))     nullify (micro%meltgraut)
    if (associated(micro%melthailt))     nullify (micro%melthailt)
    if (associated(micro%cld2raint))     nullify (micro%cld2raint)
    if (associated(micro%rimecldt))      nullify (micro%rimecldt)
    if (associated(micro%rimecldsnowt))  nullify (micro%rimecldsnowt)
    if (associated(micro%rimecldaggrt))  nullify (micro%rimecldaggrt)
    if (associated(micro%rimecldgraut))  nullify (micro%rimecldgraut)
    if (associated(micro%rimecldhailt))  nullify (micro%rimecldhailt)
    if (associated(micro%rain2icet))     nullify (micro%rain2icet)
    if (associated(micro%rain2prt))      nullify (micro%rain2prt)
    if (associated(micro%rain2snt))      nullify (micro%rain2snt)
    if (associated(micro%rain2agt))      nullify (micro%rain2agt)
    if (associated(micro%rain2grt))      nullify (micro%rain2grt)
    if (associated(micro%rain2hat))      nullify (micro%rain2hat)
    if (associated(micro%rain2ha_xtrat)) nullify (micro%rain2ha_xtrat)
    if (associated(micro%ice2raint))     nullify (micro%ice2raint)
    if (associated(micro%aggregatet))    nullify (micro%aggregatet)
    if (associated(micro%aggrselfprist)) nullify (micro%aggrselfprist)
    if (associated(micro%aggrselfsnowt)) nullify (micro%aggrselfsnowt)
    if (associated(micro%aggrprissnowt)) nullify (micro%aggrprissnowt)
    if (associated(micro%latheatvapt))   nullify (micro%latheatvapt)
    if (associated(micro%latheatfrzt))   nullify (micro%latheatfrzt)

    if (associated(micro%rei))     nullify (micro%rei)
    if (associated(micro%rel))     nullify (micro%rel)
    if (associated(micro%cldfr))   nullify (micro%cldfr)
    return
  end subroutine nullify_micro

  subroutine dealloc_micro(micro)
    type (micro_vars) :: micro

    if (associated(micro%rcp))     deallocate (micro%rcp)
    if (associated(micro%rdp))     deallocate (micro%rdp)
    if (associated(micro%rrp))     deallocate (micro%rrp)
    if (associated(micro%rpp))     deallocate (micro%rpp)
    if (associated(micro%rsp))     deallocate (micro%rsp)
    if (associated(micro%rap))     deallocate (micro%rap)
    if (associated(micro%rgp))     deallocate (micro%rgp)
    if (associated(micro%rhp))     deallocate (micro%rhp)
    if (associated(micro%ccp))     deallocate (micro%ccp)
    if (associated(micro%cdp))     deallocate (micro%cdp)
    if (associated(micro%crp))     deallocate (micro%crp)
    if (associated(micro%cpp))     deallocate (micro%cpp)
    if (associated(micro%csp))     deallocate (micro%csp)
    if (associated(micro%cap))     deallocate (micro%cap)
    if (associated(micro%cgp))     deallocate (micro%cgp)
    if (associated(micro%chp))     deallocate (micro%chp)
    if (associated(micro%cccnp))   deallocate (micro%cccnp)
    if (associated(micro%gccnp))   deallocate (micro%gccnp)
    if (associated(micro%cifnp))   deallocate (micro%cifnp)
    if (associated(micro%q2))      deallocate (micro%q2)
    if (associated(micro%q6))      deallocate (micro%q6)
    if (associated(micro%q7))      deallocate (micro%q7)

    if (associated(micro%cccmp))   deallocate (micro%cccmp)
    if (associated(micro%gccmp))   deallocate (micro%gccmp)
    if (associated(micro%cnm1p))   deallocate (micro%cnm1p)
    if (associated(micro%cnm2p))   deallocate (micro%cnm2p)
    if (associated(micro%cnm3p))   deallocate (micro%cnm3p)
    if (associated(micro%cnm8p))   deallocate (micro%cnm8p)
    if (associated(micro%md1np))   deallocate (micro%md1np)
    if (associated(micro%md2np))   deallocate (micro%md2np)
    if (associated(micro%salt_filmp)) deallocate (micro%salt_filmp)
    if (associated(micro%salt_jetp))  deallocate (micro%salt_jetp)
    if (associated(micro%salt_spmp))  deallocate (micro%salt_spmp)
    if (associated(micro%pcpvr))   deallocate (micro%pcpvr)
    if (associated(micro%pcpvp))   deallocate (micro%pcpvp)
    if (associated(micro%pcpvs))   deallocate (micro%pcpvs)
    if (associated(micro%pcpva))   deallocate (micro%pcpva)
    if (associated(micro%pcpvg))   deallocate (micro%pcpvg)
    if (associated(micro%pcpvh))   deallocate (micro%pcpvh)
    if (associated(micro%pcpvd))   deallocate (micro%pcpvd)

    if (associated(micro%accpr))   deallocate (micro%accpr)
    if (associated(micro%accpp))   deallocate (micro%accpp)
    if (associated(micro%accps))   deallocate (micro%accps)
    if (associated(micro%accpa))   deallocate (micro%accpa)
    if (associated(micro%accpg))   deallocate (micro%accpg)
    if (associated(micro%accph))   deallocate (micro%accph)
    if (associated(micro%accpd))   deallocate (micro%accpd)
    if (associated(micro%pcprr))   deallocate (micro%pcprr)
    if (associated(micro%pcprp))   deallocate (micro%pcprp)
    if (associated(micro%pcprs))   deallocate (micro%pcprs)
    if (associated(micro%pcpra))   deallocate (micro%pcpra)
    if (associated(micro%pcprg))   deallocate (micro%pcprg)
    if (associated(micro%pcprh))   deallocate (micro%pcprh)
    if (associated(micro%pcprd))   deallocate (micro%pcprd)
    if (associated(micro%pcpg))    deallocate (micro%pcpg)
    if (associated(micro%qpcpg))   deallocate (micro%qpcpg)
    if (associated(micro%dpcpg))   deallocate (micro%dpcpg)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
    if (associated(micro%nuccldr))      deallocate (micro%nuccldr)
    if (associated(micro%nuccldc))      deallocate (micro%nuccldc)
    if (associated(micro%nucicer))      deallocate (micro%nucicer)
    if (associated(micro%nucicec))      deallocate (micro%nucicec)
    if (associated(micro%inuchomr))     deallocate (micro%inuchomr)     
    if (associated(micro%inuchomc))     deallocate (micro%inuchomc) 
    if (associated(micro%inuccontr))    deallocate (micro%inuccontr)     
    if (associated(micro%inuccontc))    deallocate (micro%inuccontc)
    if (associated(micro%inucifnr))     deallocate (micro%inucifnr)     
    if (associated(micro%inucifnc))     deallocate (micro%inucifnc) 
    if (associated(micro%inuchazr))     deallocate (micro%inuchazr)     
    if (associated(micro%inuchazc))     deallocate (micro%inuchazc)
    if (associated(micro%vapliq))       deallocate (micro%vapliq)
    if (associated(micro%vapice))       deallocate (micro%vapice)
    if (associated(micro%vapcld))       deallocate (micro%vapcld)
    if (associated(micro%vaprain))      deallocate (micro%vaprain)
    if (associated(micro%vappris))      deallocate (micro%vappris)
    if (associated(micro%vapsnow))      deallocate (micro%vapsnow)
    if (associated(micro%vapaggr))      deallocate (micro%vapaggr)
    if (associated(micro%vapgrau))      deallocate (micro%vapgrau)
    if (associated(micro%vaphail))      deallocate (micro%vaphail)
    if (associated(micro%vapdriz))      deallocate (micro%vapdriz)      
    if (associated(micro%meltice))      deallocate (micro%meltice) 
    if (associated(micro%meltpris))     deallocate (micro%meltpris) 
    if (associated(micro%meltsnow))     deallocate (micro%meltsnow) 
    if (associated(micro%meltaggr))     deallocate (micro%meltaggr)
    if (associated(micro%meltgrau))     deallocate (micro%meltgrau) 
    if (associated(micro%melthail))     deallocate (micro%melthail) 
    if (associated(micro%cld2rain))     deallocate (micro%cld2rain)
    if (associated(micro%rimecld))      deallocate (micro%rimecld)
    if (associated(micro%rimecldsnow))  deallocate (micro%rimecldsnow)
    if (associated(micro%rimecldaggr))  deallocate (micro%rimecldaggr)
    if (associated(micro%rimecldgrau))  deallocate (micro%rimecldgrau)
    if (associated(micro%rimecldhail))  deallocate (micro%rimecldhail)
    if (associated(micro%rain2ice))     deallocate (micro%rain2ice)
    if (associated(micro%rain2pr))      deallocate (micro%rain2pr) 
    if (associated(micro%rain2sn))      deallocate (micro%rain2sn) 
    if (associated(micro%rain2ag))      deallocate (micro%rain2ag) 
    if (associated(micro%rain2gr))      deallocate (micro%rain2gr) 
    if (associated(micro%rain2ha))      deallocate (micro%rain2ha)
    if (associated(micro%rain2ha_xtra)) deallocate (micro%rain2ha_xtra)
    if (associated(micro%ice2rain))     deallocate (micro%ice2rain)
    if (associated(micro%aggregate))    deallocate (micro%aggregate) 
    if (associated(micro%aggrselfpris)) deallocate (micro%aggrselfpris)
    if (associated(micro%aggrselfsnow)) deallocate (micro%aggrselfsnow)
    if (associated(micro%aggrprissnow)) deallocate (micro%aggrprissnow)
    if (associated(micro%latheatvap))   deallocate (micro%latheatvap)
    if (associated(micro%latheatfrz))   deallocate (micro%latheatfrz)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
    if (associated(micro%nuccldrt))      deallocate (micro%nuccldrt)
    if (associated(micro%nuccldct))      deallocate (micro%nuccldct)
    if (associated(micro%nucicert))      deallocate (micro%nucicert)
    if (associated(micro%nucicect))      deallocate (micro%nucicect)
    if (associated(micro%inuchomrt))     deallocate (micro%inuchomrt)     
    if (associated(micro%inuchomct))     deallocate (micro%inuchomct) 
    if (associated(micro%inuccontrt))    deallocate (micro%inuccontrt)     
    if (associated(micro%inuccontct))    deallocate (micro%inuccontct)
    if (associated(micro%inucifnrt))     deallocate (micro%inucifnrt)     
    if (associated(micro%inucifnct))     deallocate (micro%inucifnct) 
    if (associated(micro%inuchazrt))     deallocate (micro%inuchazrt)     
    if (associated(micro%inuchazct))     deallocate (micro%inuchazct)
    if (associated(micro%vapliqt))       deallocate (micro%vapliqt)
    if (associated(micro%vapicet))       deallocate (micro%vapicet)
    if (associated(micro%vapcldt))       deallocate (micro%vapcldt)
    if (associated(micro%vapraint))      deallocate (micro%vapraint)
    if (associated(micro%vapprist))      deallocate (micro%vapprist)
    if (associated(micro%vapsnowt))      deallocate (micro%vapsnowt)
    if (associated(micro%vapaggrt))      deallocate (micro%vapaggrt)
    if (associated(micro%vapgraut))      deallocate (micro%vapgraut)
    if (associated(micro%vaphailt))      deallocate (micro%vaphailt)
    if (associated(micro%vapdrizt))      deallocate (micro%vapdrizt)      
    if (associated(micro%melticet))      deallocate (micro%melticet) 
    if (associated(micro%meltprist))     deallocate (micro%meltprist) 
    if (associated(micro%meltsnowt))     deallocate (micro%meltsnowt) 
    if (associated(micro%meltaggrt))     deallocate (micro%meltaggrt)
    if (associated(micro%meltgraut))     deallocate (micro%meltgraut) 
    if (associated(micro%melthailt))     deallocate (micro%melthailt) 
    if (associated(micro%cld2raint))     deallocate (micro%cld2raint)
    if (associated(micro%rimecldt))      deallocate (micro%rimecldt)
    if (associated(micro%rimecldsnowt))  deallocate (micro%rimecldsnowt)
    if (associated(micro%rimecldaggrt))  deallocate (micro%rimecldaggrt)
    if (associated(micro%rimecldgraut))  deallocate (micro%rimecldgraut)
    if (associated(micro%rimecldhailt))  deallocate (micro%rimecldhailt)
    if (associated(micro%rain2icet))     deallocate (micro%rain2icet)
    if (associated(micro%rain2prt))      deallocate (micro%rain2prt) 
    if (associated(micro%rain2snt))      deallocate (micro%rain2snt) 
    if (associated(micro%rain2agt))      deallocate (micro%rain2agt) 
    if (associated(micro%rain2grt))      deallocate (micro%rain2grt) 
    if (associated(micro%rain2hat))      deallocate (micro%rain2hat)
    if (associated(micro%rain2ha_xtrat)) deallocate (micro%rain2ha_xtrat)
    if (associated(micro%ice2raint))     deallocate (micro%ice2raint)
    if (associated(micro%aggregatet))    deallocate (micro%aggregatet) 
    if (associated(micro%aggrselfprist)) deallocate (micro%aggrselfprist)
    if (associated(micro%aggrselfsnowt)) deallocate (micro%aggrselfsnowt)
    if (associated(micro%aggrprissnowt)) deallocate (micro%aggrprissnowt)
    if (associated(micro%latheatvapt))   deallocate (micro%latheatvapt)
    if (associated(micro%latheatfrzt))   deallocate (micro%latheatfrzt)

    if (associated(micro%rei))     deallocate (micro%rei)
    if (associated(micro%rel))     deallocate (micro%rel)
    if (associated(micro%cldfr))   deallocate (micro%cldfr)

    return
  end subroutine dealloc_micro


  subroutine filltab_micro(micro,microm,imean,n1,n2,n3,ng)
    type (micro_vars) :: micro,microm
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer(kind=i8) :: npts
    real, pointer :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n1*n2*n3
    if (associated(micro%rcp))   &
         call InsertVTab (micro%rcp,microm%rcp  &
         ,ng, npts, imean,  &
         'RCP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rdp))   &
         call InsertVTab (micro%rdp,microm%rdp  &
         ,ng, npts, imean,  &
         'RDP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rrp))   &
         call InsertVTab (micro%rrp,microm%rrp  &
         ,ng, npts, imean,  &
         'RRP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rpp))   &
         call InsertVTab (micro%rpp,microm%rpp  &
         ,ng, npts, imean,  &
         'RPP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rsp))   &
         call InsertVTab (micro%rsp,microm%rsp  &
         ,ng, npts, imean,  &
         'RSP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rap))   &
         call InsertVTab (micro%rap,microm%rap  &
         ,ng, npts, imean,  &
         'RAP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rgp))   &
         call InsertVTab (micro%rgp,microm%rgp  &
         ,ng, npts, imean,  &
         'RGP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rhp))   &
         call InsertVTab (micro%rhp,microm%rhp  &
         ,ng, npts, imean,  &
         'RHP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%ccp))   &
         call InsertVTab (micro%ccp,microm%ccp  &
         ,ng, npts, imean,  &
         'CCP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cdp))   &
         call InsertVTab (micro%cdp,microm%cdp  &
         ,ng, npts, imean,  &
         'CDP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%crp))   &
         call InsertVTab (micro%crp,microm%crp  &
         ,ng, npts, imean,  &
         'CRP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cpp))   &
         call InsertVTab (micro%cpp,microm%cpp  &
         ,ng, npts, imean,  &
         'CPP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%csp))   &
         call InsertVTab (micro%csp,microm%csp  &
         ,ng, npts, imean,  &
         'CSP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cap))   &
         call InsertVTab (micro%cap,microm%cap  &
         ,ng, npts, imean,  &
         'CAP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cgp))   &
         call InsertVTab (micro%cgp,microm%cgp  &
         ,ng, npts, imean,  &
         'CGP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%chp))   &
         call InsertVTab (micro%chp,microm%chp  &
         ,ng, npts, imean,  &
         'CHP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cccnp)) &
         call InsertVTab (micro%cccnp,microm%cccnp  &
         ,ng, npts, imean,  &
         'CCCNP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%gccnp)) &
         call InsertVTab (micro%gccnp,microm%gccnp  &
         ,ng, npts, imean,  &
         'GCCNP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cifnp)) &
         call InsertVTab (micro%cifnp,microm%cifnp  &
         ,ng, npts, imean,  &
         'CIFNP :3:hist:anal:mpti:mpt3:mpt1')

    if (associated(micro%q2))   &
         call InsertVTab (micro%q2,microm%q2  &
         ,ng, npts, imean,  &
         'Q2 :3:hist:anal:mpti:mpt3')
    if (associated(micro%q6)) &
         call InsertVTab (micro%q6,microm%q6  &
         ,ng, npts, imean,  &
         'Q6 :3:hist:anal:mpti:mpt3')
    if (associated(micro%q7)) &
         call InsertVTab (micro%q7,microm%q7  &
         ,ng, npts, imean,  &
         'Q7 :3:hist:anal:mpti:mpt3')

    if (associated(micro%cccmp)) &
         call InsertVTab (micro%cccmp,microm%cccmp  &
         ,ng, npts, imean,  &
         'CCCMP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%gccmp)) &
         call InsertVTab (micro%gccmp,microm%gccmp  &
         ,ng, npts, imean,  &
         'GCCMP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cnm1p)) &
         call InsertVTab (micro%cnm1p,microm%cnm1p  &
         ,ng, npts, imean,  &
         'CNM1P :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cnm2p)) &
         call InsertVTab (micro%cnm2p,microm%cnm2p  &
         ,ng, npts, imean,  &
         'CNM2P :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cnm3p)) &
         call InsertVTab (micro%cnm3p,microm%cnm3p  &
         ,ng, npts, imean,  &
         'CNM3P :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cnm8p)) &
         call InsertVTab (micro%cnm8p,microm%cnm8p  &
         ,ng, npts, imean,  &
         'CNM8P :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%md1np)) &
         call InsertVTab (micro%md1np,microm%md1np  &
         ,ng, npts, imean,  &
         'MD1NP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%md2np)) &
         call InsertVTab (micro%md2np,microm%md2np  &
         ,ng, npts, imean,  &
         'MD2NP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%salt_filmp)) &
         call InsertVTab (micro%salt_filmp,microm%salt_filmp  &
         ,ng, npts, imean,  &
         'SALT_FILMP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%salt_jetp)) &
         call InsertVTab (micro%salt_jetp,microm%salt_jetp  &
         ,ng, npts, imean,  &
         'SALT_JETP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%salt_spmp)) &
         call InsertVTab (micro%salt_spmp,microm%salt_spmp  &
         ,ng, npts, imean,  &
         'SALT_SPMP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rei))   &
         call InsertVTab (micro%rei,microm%rei  &
         ,ng, npts, imean,  &
         'REI :3:hist:anal:mpti:mpt3')
    if (associated(micro%rel))   &
         call InsertVTab (micro%rel,microm%rel  &
         ,ng, npts, imean,  &
         'REL :3:hist:anal:mpti:mpt3')
    if (associated(micro%cldfr))   &
         call InsertVTab (micro%cldfr,microm%cldfr  &
         ,ng, npts, imean,  &
         'CLDFR :3:hist:anal:mpti:mpt3')


    !VERTICAL PRECIPITATION RATES
    if (associated(micro%pcpvr)) &
         call InsertVTab (micro%pcpvr,microm%pcpvr  &
         ,ng, npts, imean,  &
         'PCPVR :3:hist:anal:mpti:mpt3')
    if (associated(micro%pcpvp)) &
         call InsertVTab (micro%pcpvp,microm%pcpvp  &
         ,ng, npts, imean,  &
         'PCPVP :3:hist:anal:mpti:mpt3')
    if (associated(micro%pcpvs)) &
         call InsertVTab (micro%pcpvs,microm%pcpvs  &
         ,ng, npts, imean,  &
         'PCPVS :3:hist:anal:mpti:mpt3')
    if (associated(micro%pcpva)) &
         call InsertVTab (micro%pcpva,microm%pcpva  &
         ,ng, npts, imean,  &
         'PCPVA :3:hist:anal:mpti:mpt3')
    if (associated(micro%pcpvg)) &
         call InsertVTab (micro%pcpvg,microm%pcpvg  &
         ,ng, npts, imean,  &
         'PCPVG :3:hist:anal:mpti:mpt3')
    if (associated(micro%pcpvh)) &
         call InsertVTab (micro%pcpvh,microm%pcpvh  &
         ,ng, npts, imean,  &
         'PCPVH :3:hist:anal:mpti:mpt3')
    if (associated(micro%pcpvd)) &
         call InsertVTab (micro%pcpvd,microm%pcpvd  &
         ,ng, npts, imean,  &
         'PCPVD :3:hist:anal:mpti:mpt3')

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (instantaneous)
    if (associated(micro%nuccldr)) &
         call InsertVTab (micro%nuccldr,microm%nuccldr  &
         ,ng, npts, imean,  &
         'NUCCLDR :3:hist:anal:mpt3')
    if (associated(micro%nuccldc)) &
         call InsertVTab (micro%nuccldc,microm%nuccldc  &
         ,ng, npts, imean,  &
         'NUCCLDC :3:hist:anal:mpt3')

    if (associated(micro%nucicer)) &
         call InsertVTab (micro%nucicer,microm%nucicer  &
         ,ng, npts, imean,  &
         'NUCICER :3:hist:anal:mpt3')
    if (associated(micro%nucicec)) &
         call InsertVTab (micro%nucicec,microm%nucicec  &
         ,ng, npts, imean,  &
         'NUCICEC :3:hist:anal:mpt3')
    if (associated(micro%inuchomr)) &
         call InsertVTab (micro%inuchomr,microm%inuchomr  &
         ,ng, npts, imean,  &
         'INUCHOMR :3:hist:anal:mpt3')
    if (associated(micro%inuchomc)) &
         call InsertVTab (micro%inuchomc,microm%inuchomc  &
         ,ng, npts, imean,  &
         'INUCHOMC :3:hist:anal:mpt3')
    if (associated(micro%inuccontr)) &
         call InsertVTab (micro%inuccontr,microm%inuccontr  &
         ,ng, npts, imean,  &
         'INUCCONTR :3:hist:anal:mpt3')
    if (associated(micro%inuccontc)) &
         call InsertVTab (micro%inuccontc,microm%inuccontc  &
         ,ng, npts, imean,  &
         'INUCCONTC :3:hist:anal:mpt3')
    if (associated(micro%inucifnr)) &
         call InsertVTab (micro%inucifnr,microm%inucifnr  &
         ,ng, npts, imean,  &
         'INUCIFNR :3:hist:anal:mpt3')
    if (associated(micro%inucifnc)) &
         call InsertVTab (micro%inucifnc,microm%inucifnc  &
         ,ng, npts, imean,  &
         'INUCIFNC :3:hist:anal:mpt3')
    if (associated(micro%inuchazr)) &
         call InsertVTab (micro%inuchazr,microm%inuchazr  &
         ,ng, npts, imean,  &
         'INUCHAZR :3:hist:anal:mpt3')
    if (associated(micro%inuchazc)) &
         call InsertVTab (micro%inuchazc,microm%inuchazc  &
         ,ng, npts, imean,  &
         'INUCHAZC :3:hist:anal:mpt3')

    if (associated(micro%vapliq)) &
         call InsertVTab (micro%vapliq,microm%vapliq  &
         ,ng, npts, imean,  &
         'VAPLIQ :3:hist:anal:mpt3')
    if (associated(micro%vapice)) &
         call InsertVTab (micro%vapice,microm%vapice  &
         ,ng, npts, imean,  &
         'VAPICE :3:hist:anal:mpt3')
    if (associated(micro%vapcld)) &
         call InsertVTab (micro%vapcld,microm%vapcld  &
         ,ng, npts, imean,  &
         'VAPCLD :3:hist:anal:mpt3')
    if (associated(micro%vaprain)) &
         call InsertVTab (micro%vaprain,microm%vaprain  &
         ,ng, npts, imean,  &
         'VAPRAIN :3:hist:anal:mpt3')
    if (associated(micro%vappris)) &
         call InsertVTab (micro%vappris,microm%vappris  &
         ,ng, npts, imean,  &
         'VAPPRIS :3:hist:anal:mpt3')
    if (associated(micro%vapsnow)) &
         call InsertVTab (micro%vapsnow,microm%vapsnow  &
         ,ng, npts, imean,  &
         'VAPSNOW :3:hist:anal:mpt3')
    if (associated(micro%vapaggr)) &
         call InsertVTab (micro%vapaggr,microm%vapaggr  &
         ,ng, npts, imean,  &
         'VAPAGGR :3:hist:anal:mpt3')
    if (associated(micro%vapgrau)) &
         call InsertVTab (micro%vapgrau,microm%vapgrau  &
         ,ng, npts, imean,  &
         'VAPGRAU :3:hist:anal:mpt3')
    if (associated(micro%vaphail)) &
         call InsertVTab (micro%vaphail,microm%vaphail  &
         ,ng, npts, imean,  &
         'VAPHAIL :3:hist:anal:mpt3')
    if (associated(micro%vapdriz)) &
         call InsertVTab (micro%vapdriz,microm%vapdriz  &
         ,ng, npts, imean,  &
         'VAPDRIZ :3:hist:anal:mpt3')

    if (associated(micro%meltice)) &
         call InsertVTab (micro%meltice,microm%meltice  &
         ,ng, npts, imean,  &
         'MELTICE :3:hist:anal:mpt3')
    if (associated(micro%meltpris)) &
         call InsertVTab (micro%meltpris,microm%meltpris  &
         ,ng, npts, imean,  &
         'MELTPRIS :3:hist:anal:mpt3')
    if (associated(micro%meltsnow)) &
         call InsertVTab (micro%meltsnow,microm%meltsnow  &
         ,ng, npts, imean,  &
         'MELTSNOW :3:hist:anal:mpt3')
    if (associated(micro%meltaggr)) &
         call InsertVTab (micro%meltaggr,microm%meltaggr  &
         ,ng, npts, imean,  &
         'MELTAGGR :3:hist:anal:mpt3')
    if (associated(micro%meltgrau)) &
         call InsertVTab (micro%meltgrau,microm%meltgrau  &
         ,ng, npts, imean,  &
         'MELTGRAU :3:hist:anal:mpt3')
    if (associated(micro%melthail)) &
         call InsertVTab (micro%melthail,microm%melthail  &
         ,ng, npts, imean,  &
         'MELTHAIL :3:hist:anal:mpt3')

    if (associated(micro%cld2rain)) &
         call InsertVTab (micro%cld2rain,microm%cld2rain  &
         ,ng, npts, imean,  &
         'CLD2RAIN :3:hist:anal:mpt3')
    if (associated(micro%rimecld)) &
         call InsertVTab (micro%rimecld,microm%rimecld  &
         ,ng, npts, imean,  &
         'RIMECLD :3:hist:anal:mpt3')
    if (associated(micro%rimecldsnow)) &
         call InsertVTab (micro%rimecldsnow,microm%rimecldsnow  &
         ,ng, npts, imean,  &
         'RIMECLDSNOW :3:hist:anal:mpt3')
    if (associated(micro%rimecldaggr)) &
         call InsertVTab (micro%rimecldaggr,microm%rimecldaggr  &
         ,ng, npts, imean,  &
         'RIMECLDAGGR :3:hist:anal:mpt3')
    if (associated(micro%rimecldgrau)) &
         call InsertVTab (micro%rimecldgrau,microm%rimecldgrau  &
         ,ng, npts, imean,  &
         'RIMECLDGRAU :3:hist:anal:mpt3')
    if (associated(micro%rimecldhail)) &
         call InsertVTab (micro%rimecldhail,microm%rimecldhail  &
         ,ng, npts, imean,  &
         'RIMECLDHAIL :3:hist:anal:mpt3')

    if (associated(micro%rain2ice)) &
         call InsertVTab (micro%rain2ice,microm%rain2ice  &
         ,ng, npts, imean,  &
         'RAIN2ICE :3:hist:anal:mpt3')
    if (associated(micro%rain2pr)) &
         call InsertVTab (micro%rain2pr,microm%rain2pr  &
         ,ng, npts, imean,  &
         'RAIN2PR :3:hist:anal:mpt3')
    if (associated(micro%rain2sn)) &
         call InsertVTab (micro%rain2sn,microm%rain2sn  &
         ,ng, npts, imean,  &
         'RAIN2SN :3:hist:anal:mpt3')
    if (associated(micro%rain2ag)) &
         call InsertVTab (micro%rain2ag,microm%rain2ag  &
         ,ng, npts, imean,  &
         'RAIN2AG :3:hist:anal:mpt3')
    if (associated(micro%rain2gr)) &
         call InsertVTab (micro%rain2gr,microm%rain2gr  &
         ,ng, npts, imean,  &
         'RAIN2GR :3:hist:anal:mpt3')
    if (associated(micro%rain2ha)) &
         call InsertVTab (micro%rain2ha,microm%rain2ha  &
         ,ng, npts, imean,  &
         'RAIN2HA :3:hist:anal:mpt3')
    if (associated(micro%rain2ha_xtra)) &
         call InsertVTab (micro%rain2ha_xtra,microm%rain2ha_xtra  &
         ,ng, npts, imean,  &
         'RAIN2HA_XTRA :3:hist:anal:mpt3')
    if (associated(micro%ice2rain)) &
         call InsertVTab (micro%ice2rain,microm%ice2rain  &
         ,ng, npts, imean,  &
         'ICE2RAIN :3:hist:anal:mpt3')

    if (associated(micro%aggregate)) &
         call InsertVTab (micro%aggregate,microm%aggregate  &
         ,ng, npts, imean,  &
         'AGGREGATE :3:hist:anal:mpt3')
    if (associated(micro%aggrselfpris)) &
         call InsertVTab (micro%aggrselfpris,microm%aggrselfpris  &
         ,ng, npts, imean,  &
         'AGGRSELFPRIS :3:hist:anal:mpt3')
    if (associated(micro%aggrselfsnow)) &
         call InsertVTab (micro%aggrselfsnow,microm%aggrselfsnow  &
         ,ng, npts, imean,  &
         'AGGRSELFSNOW :3:hist:anal:mpt3')
    if (associated(micro%aggrprissnow)) &
         call InsertVTab (micro%aggrprissnow,microm%aggrprissnow  &
         ,ng, npts, imean,  &
         'AGGRPRISSNOW :3:hist:anal:mpt3')

    if (associated(micro%latheatvap)) &
         call InsertVTab (micro%latheatvap,microm%latheatvap  &
         ,ng, npts, imean,  &
         'LATHEATVAP :3:hist:anal:mpt3')
    if (associated(micro%latheatfrz)) &
         call InsertVTab (micro%latheatfrz,microm%latheatfrz  &
         ,ng, npts, imean,  &
         'LATHEATFRZ :3:hist:anal:mpt3')
    !END MICRO BUDGET PROCESSES (instantaneous)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
    if (associated(micro%nuccldrt)) &
         call InsertVTab (micro%nuccldrt,microm%nuccldrt  &
         ,ng, npts, imean,  &
         'NUCCLDRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%nuccldct)) &
         call InsertVTab (micro%nuccldct,microm%nuccldct  &
         ,ng, npts, imean,  &
         'NUCCLDCT :3:hist:anal:mpti:mpt3')

    if (associated(micro%nucicert)) &
         call InsertVTab (micro%nucicert,microm%nucicert  &
         ,ng, npts, imean,  &
         'NUCICERT :3:hist:anal:mpti:mpt3')
    if (associated(micro%nucicect)) &
         call InsertVTab (micro%nucicect,microm%nucicect  &
         ,ng, npts, imean,  &
         'NUCICECT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuchomrt)) &
         call InsertVTab (micro%inuchomrt,microm%inuchomrt  &
         ,ng, npts, imean,  &
         'INUCHOMRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuchomct)) &
         call InsertVTab (micro%inuchomct,microm%inuchomct  &
         ,ng, npts, imean,  &
         'INUCHOMCT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuccontrt)) &
         call InsertVTab (micro%inuccontrt,microm%inuccontrt  &
         ,ng, npts, imean,  &
         'INUCCONTRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuccontct)) &
         call InsertVTab (micro%inuccontct,microm%inuccontct  &
         ,ng, npts, imean,  &
         'INUCCONTCT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inucifnrt)) &
         call InsertVTab (micro%inucifnrt,microm%inucifnrt  &
         ,ng, npts, imean,  &
         'INUCIFNRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inucifnct)) &
         call InsertVTab (micro%inucifnct,microm%inucifnct  &
         ,ng, npts, imean,  &
         'INUCIFNCT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuchazrt)) &
         call InsertVTab (micro%inuchazrt,microm%inuchazrt  &
         ,ng, npts, imean,  &
         'INUCHAZRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuchazct)) &
         call InsertVTab (micro%inuchazct,microm%inuchazct  &
         ,ng, npts, imean,  &
         'INUCHAZCT :3:hist:anal:mpti:mpt3')

    if (associated(micro%vapliqt)) &
         call InsertVTab (micro%vapliqt,microm%vapliqt  &
         ,ng, npts, imean,  &
         'VAPLIQT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapicet)) &
         call InsertVTab (micro%vapicet,microm%vapicet  &
         ,ng, npts, imean,  &
         'VAPICET :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapcldt)) &
         call InsertVTab (micro%vapcldt,microm%vapcldt  &
         ,ng, npts, imean,  &
         'VAPCLDT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapraint)) &
         call InsertVTab (micro%vapraint,microm%vapraint  &
         ,ng, npts, imean,  &
         'VAPRAINT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapprist)) &
         call InsertVTab (micro%vapprist,microm%vapprist  &
         ,ng, npts, imean,  &
         'VAPPRIST :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapsnowt)) &
         call InsertVTab (micro%vapsnowt,microm%vapsnowt  &
         ,ng, npts, imean,  &
         'VAPSNOWT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapaggrt)) &
         call InsertVTab (micro%vapaggrt,microm%vapaggrt  &
         ,ng, npts, imean,  &
         'VAPAGGRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapgraut)) &
         call InsertVTab (micro%vapgraut,microm%vapgraut  &
         ,ng, npts, imean,  &
         'VAPGRAUT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vaphailt)) &
         call InsertVTab (micro%vaphailt,microm%vaphailt  &
         ,ng, npts, imean,  &
         'VAPHAILT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapdrizt)) &
         call InsertVTab (micro%vapdrizt,microm%vapdrizt  &
         ,ng, npts, imean,  &
         'VAPDRIZT :3:hist:anal:mpti:mpt3')

    if (associated(micro%melticet)) &
         call InsertVTab (micro%melticet,microm%melticet  &
         ,ng, npts, imean,  &
         'MELTICET :3:hist:anal:mpti:mpt3')
    if (associated(micro%meltprist)) &
         call InsertVTab (micro%meltprist,microm%meltprist  &
         ,ng, npts, imean,  &
         'MELTPRIST :3:hist:anal:mpti:mpt3')
    if (associated(micro%meltsnowt)) &
         call InsertVTab (micro%meltsnowt,microm%meltsnowt  &
         ,ng, npts, imean,  &
         'MELTSNOWT :3:hist:anal:mpti:mpt3')
    if (associated(micro%meltaggrt)) &
         call InsertVTab (micro%meltaggrt,microm%meltaggrt  &
         ,ng, npts, imean,  &
         'MELTAGGRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%meltgraut)) &
         call InsertVTab (micro%meltgraut,microm%meltgraut  &
         ,ng, npts, imean,  &
         'MELTGRAUT :3:hist:anal:mpti:mpt3')
    if (associated(micro%melthailt)) &
         call InsertVTab (micro%melthailt,microm%melthailt  &
         ,ng, npts, imean,  &
         'MELTHAILT :3:hist:anal:mpti:mpt3')

    if (associated(micro%cld2raint)) &
         call InsertVTab (micro%cld2raint,microm%cld2raint  &
         ,ng, npts, imean,  &
         'CLD2RAINT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldt)) &
         call InsertVTab (micro%rimecldt,microm%rimecldt  &
         ,ng, npts, imean,  &
         'RIMECLDT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldsnowt)) &
         call InsertVTab (micro%rimecldsnowt,microm%rimecldsnowt  &
         ,ng, npts, imean,  &
         'RIMECLDSNOWT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldaggrt)) &
         call InsertVTab (micro%rimecldaggrt,microm%rimecldaggrt  &
         ,ng, npts, imean,  &
         'RIMECLDAGGRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldgraut)) &
         call InsertVTab (micro%rimecldgraut,microm%rimecldgraut  &
         ,ng, npts, imean,  &
         'RIMECLDGRAUT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldhailt)) &
         call InsertVTab (micro%rimecldhailt,microm%rimecldhailt  &
         ,ng, npts, imean,  &
         'RIMECLDHAILT :3:hist:anal:mpti:mpt3')

    if (associated(micro%rain2icet)) &
         call InsertVTab (micro%rain2icet,microm%rain2icet  &
         ,ng, npts, imean,  &
         'RAIN2ICET :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2prt)) &
         call InsertVTab (micro%rain2prt,microm%rain2prt  &
         ,ng, npts, imean,  &
         'RAIN2PRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2snt)) &
         call InsertVTab (micro%rain2snt,microm%rain2snt  &
         ,ng, npts, imean,  &
         'RAIN2SNT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2agt)) &
         call InsertVTab (micro%rain2agt,microm%rain2agt  &
         ,ng, npts, imean,  &
         'RAIN2AGT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2grt)) &
         call InsertVTab (micro%rain2grt,microm%rain2grt  &
         ,ng, npts, imean,  &
         'RAIN2GRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2hat)) &
         call InsertVTab (micro%rain2hat,microm%rain2hat  &
         ,ng, npts, imean,  &
         'RAIN2HAT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2ha_xtrat)) &
         call InsertVTab (micro%rain2ha_xtrat,microm%rain2ha_xtrat  &
         ,ng, npts, imean,  &
         'RAIN2HA_XTRAT :3:hist:anal:mpti:mpt3')
    if (associated(micro%ice2raint)) &
         call InsertVTab (micro%ice2raint,microm%ice2raint  &
         ,ng, npts, imean,  &
         'ICE2RAINT :3:hist:anal:mpti:mpt3')

    if (associated(micro%aggregatet)) &
         call InsertVTab (micro%aggregatet,microm%aggregatet  &
         ,ng, npts, imean,  &
         'AGGREGATET :3:hist:anal:mpti:mpt3')
    if (associated(micro%aggrselfprist)) &
         call InsertVTab (micro%aggrselfprist,microm%aggrselfprist  &
         ,ng, npts, imean,  &
         'AGGRSELFPRIST :3:hist:anal:mpti:mpt3')
    if (associated(micro%aggrselfsnowt)) &
         call InsertVTab (micro%aggrselfsnowt,microm%aggrselfsnowt  &
         ,ng, npts, imean,  &
         'AGGRSELFSNOWT :3:hist:anal:mpti:mpt3')
    if (associated(micro%aggrprissnowt)) &
         call InsertVTab (micro%aggrprissnowt,microm%aggrprissnowt  &
         ,ng, npts, imean,  &
         'AGGRPRISSNOWT :3:hist:anal:mpti:mpt3')

    if (associated(micro%latheatvapt)) &
         call InsertVTab (micro%latheatvapt,microm%latheatvapt  &
         ,ng, npts, imean,  &
         'LATHEATVAPT :3:hist:anal:mpti:mpt3')
    if (associated(micro%latheatfrzt)) &
         call InsertVTab (micro%latheatfrzt,microm%latheatfrzt  &
         ,ng, npts, imean,  &
         'LATHEATFRZT :3:hist:anal:mpti:mpt3')
    !END MICRO BUDGET PROCECCES (totals)

    npts=n2*n3
    if (associated(micro%accpr)) &
         call InsertVTab (micro%accpr,microm%accpr  &
         ,ng, npts, imean,  &
         'ACCPR :2:hist:anal:mpti:mpt3')
    if (associated(micro%accpp)) &
         call InsertVTab (micro%accpp,microm%accpp  &
         ,ng, npts, imean,  &
         'ACCPP :2:hist:anal:mpti:mpt3')
    if (associated(micro%accps)) &
         call InsertVTab (micro%accps,microm%accps  &
         ,ng, npts, imean,  &
         'ACCPS :2:hist:anal:mpti:mpt3')
    if (associated(micro%accpa)) &
         call InsertVTab (micro%accpa,microm%accpa  &
         ,ng, npts, imean,  &
         'ACCPA :2:hist:anal:mpti:mpt3')
    if (associated(micro%accpg)) &
         call InsertVTab (micro%accpg,microm%accpg  &
         ,ng, npts, imean,  &
         'ACCPG :2:hist:anal:mpti:mpt3')
    if (associated(micro%accph)) &
         call InsertVTab (micro%accph,microm%accph  &
         ,ng, npts, imean,  &
         'ACCPH :2:hist:anal:mpti:mpt3')
    if (associated(micro%accpd)) &
         call InsertVTab (micro%accpd,microm%accpd  &
         ,ng, npts, imean,  &
         'ACCPD :2:hist:anal:mpti:mpt3')
    if (associated(micro%pcprr)) &
         call InsertVTab (micro%pcprr,microm%pcprr  &
         ,ng, npts, imean,  &
         'PCPRR :2:hist:anal:mpt3')
    if (associated(micro%pcprp)) &
         call InsertVTab (micro%pcprp,microm%pcprp  &
         ,ng, npts, imean,  &
         'PCPRP :2:hist:anal:mpt3')
    if (associated(micro%pcprs)) &
         call InsertVTab (micro%pcprs,microm%pcprs  &
         ,ng, npts, imean,  &
         'PCPRS :2:hist:anal:mpt3')
    if (associated(micro%pcpra)) &
         call InsertVTab (micro%pcpra,microm%pcpra  &
         ,ng, npts, imean,  &
         'PCPRA :2:hist:anal:mpt3')
    if (associated(micro%pcprg)) &
         call InsertVTab (micro%pcprg,microm%pcprg  &
         ,ng, npts, imean,  &
         'PCPRG :2:hist:anal:mpt3')
    if (associated(micro%pcprh)) &
         call InsertVTab (micro%pcprh,microm%pcprh  &
         ,ng, npts, imean,  &
         'PCPRH :2:hist:anal:mpt3')
    if (associated(micro%pcprd)) &
         call InsertVTab (micro%pcprd,microm%pcprd  &
         ,ng, npts, imean,  &
         'PCPRD :2:hist:anal:mpt3')
    if (associated(micro%pcpg)) &
         call InsertVTab (micro%pcpg,microm%pcpg  &
         ,ng, npts, imean,  &
         'PCPG :2:hist:mpti:mpt3')
    if (associated(micro%qpcpg)) &
         call InsertVTab (micro%qpcpg,microm%qpcpg  &
         ,ng, npts, imean,  &
         'QPCPG :2:hist:mpti:mpt3')
    if (associated(micro%dpcpg)) &
         call InsertVTab (micro%dpcpg,microm%dpcpg  &
         ,ng, npts, imean,  &
         'DPCPG :2:hist:mpti:mpt3')

    return
  end subroutine filltab_micro




  subroutine DeepCopyToMicroFields(oneMicroFields, from)
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    character(len=*), intent(in) :: from

    logical :: assOld, assNew
    character(len=*), parameter :: h="**(DeepCopyToMicroFields)**"

    if (copyTo /= "") then
       call fatal_error(h//" invoked from "//trim(adjustl(from))//&
            " just after invoked from "//trim(adjustl(copyTo)))
    end if

    copyTo=from
    copyFrom=""

    assOld=associated(micro_g(1)%rcp)
    assNew=associated(oneMicroFields%rcp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rcp, oneMicroFields%rcp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rcp associated but oneMicroFields%rcp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rcp not associated and oneMicroFields%rcp associated")
    end if

    assOld=associated(micro_g(1)%rdp)
    assNew=associated(oneMicroFields%rdp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rdp, oneMicroFields%rdp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rdp associated but oneMicroFields%rdp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rdp not associated and oneMicroFields%rdp associated")
    end if

    assOld=associated(micro_g(1)%rrp)
    assNew=associated(oneMicroFields%rrp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rrp, oneMicroFields%rrp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rrp associated but oneMicroFields%rrp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rrp not associated and oneMicroFields%rrp associated")
    end if

    assOld=associated(micro_g(1)%rpp)
    assNew=associated(oneMicroFields%rpp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rpp, oneMicroFields%rpp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rpp associated but oneMicroFields%rpp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rpp not associated and oneMicroFields%rpp associated")
    end if

    assOld=associated(micro_g(1)%rsp)
    assNew=associated(oneMicroFields%rsp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rsp, oneMicroFields%rsp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rsp associated but oneMicroFields%rsp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rsp not associated and oneMicroFields%rsp associated")
    end if
    assOld=associated(micro_g(1)%rap)
    assNew=associated(oneMicroFields%rap)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rap, oneMicroFields%rap)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rap associated but oneMicroFields%rap not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rap not associated and oneMicroFields%rap associated")
    end if

    assOld=associated(micro_g(1)%rgp)
    assNew=associated(oneMicroFields%rgp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rgp, oneMicroFields%rgp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rgp associated but oneMicroFields%rgp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rgp not associated and oneMicroFields%rgp associated")
    end if

    assOld=associated(micro_g(1)%rhp)
    assNew=associated(oneMicroFields%rhp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rhp, oneMicroFields%rhp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rhp associated but oneMicroFields%rhp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rhp not associated and oneMicroFields%rhp associated")
    end if

    assOld=associated(micro_g(1)%ccp)
    assNew=associated(oneMicroFields%ccp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%ccp, oneMicroFields%ccp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%ccp associated but oneMicroFields%ccp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%ccp not associated and oneMicroFields%ccp associated")
    end if

    assOld=associated(micro_g(1)%cdp)
    assNew=associated(oneMicroFields%cdp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cdp, oneMicroFields%cdp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cdp associated but oneMicroFields%cdp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cdp not associated and oneMicroFields%cdp associated")
    end if

    assOld=associated(micro_g(1)%crp)
    assNew=associated(oneMicroFields%crp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%crp, oneMicroFields%crp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%crp associated but oneMicroFields%crp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%crp not associated and oneMicroFields%crp associated")
    end if

    assOld=associated(micro_g(1)%cpp)
    assNew=associated(oneMicroFields%cpp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cpp, oneMicroFields%cpp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cpp associated but oneMicroFields%cpp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cpp not associated and oneMicroFields%cpp associated")
    end if

    assOld=associated(micro_g(1)%csp)
    assNew=associated(oneMicroFields%csp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%csp, oneMicroFields%csp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%csp associated but oneMicroFields%csp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%csp not associated and oneMicroFields%csp associated")
    end if

    assOld=associated(micro_g(1)%cap)
    assNew=associated(oneMicroFields%cap)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cap, oneMicroFields%cap)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cap associated but oneMicroFields%cap not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cap not associated and oneMicroFields%cap associated")
    end if

    assOld=associated(micro_g(1)%cgp)
    assNew=associated(oneMicroFields%cgp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cgp, oneMicroFields%cgp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cgp associated but oneMicroFields%cgp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cgp not associated and oneMicroFields%cgp associated")
    end if

    assOld=associated(micro_g(1)%chp)
    assNew=associated(oneMicroFields%chp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%chp, oneMicroFields%chp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%chp associated but oneMicroFields%chp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%chp not associated and oneMicroFields%chp associated")
    end if

    assOld=associated(micro_g(1)%cccnp)
    assNew=associated(oneMicroFields%cccnp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cccnp, oneMicroFields%cccnp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cccnp associated but oneMicroFields%cccnp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cccnp not associated and oneMicroFields%cccnp associated")
    end if

    assOld=associated(micro_g(1)%gccnp)
    assNew=associated(oneMicroFields%gccnp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%gccnp, oneMicroFields%gccnp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%gccnp associated but oneMicroFields%gccnp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%gccnp not associated and oneMicroFields%gccnp associated")
    end if

    assOld=associated(micro_g(1)%cifnp)
    assNew=associated(oneMicroFields%cifnp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cifnp, oneMicroFields%cifnp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cifnp associated but oneMicroFields%cifnp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cifnp not associated and oneMicroFields%cifnp associated")
    end if

    assOld=associated(micro_g(1)%q2)
    assNew=associated(oneMicroFields%q2)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%q2, oneMicroFields%q2)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%q2 associated but oneMicroFields%q2 not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%q2 not associated and oneMicroFields%q2 associated")
    end if

    assOld=associated(micro_g(1)%q6)
    assNew=associated(oneMicroFields%q6)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%q6, oneMicroFields%q6)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%q6 associated but oneMicroFields%q6 not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%q6 not associated and oneMicroFields%q6 associated")
    end if

    assOld=associated(micro_g(1)%q7)
    assNew=associated(oneMicroFields%q7)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%q7, oneMicroFields%q7)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%q7 associated but oneMicroFields%q7 not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%q7 not associated and oneMicroFields%q7 associated")
    end if

    assOld=associated(micro_g(1)%rei)
    assNew=associated(oneMicroFields%rei)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rei, oneMicroFields%rei)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rei associated but oneMicroFields%rei not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rei not associated and oneMicroFields%rei associated")
    end if

    assOld=associated(micro_g(1)%rel)
    assNew=associated(oneMicroFields%rel)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rel, oneMicroFields%rel)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rel associated but oneMicroFields%rel not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rel not associated and oneMicroFields%rel associated")
    end if

    assOld=associated(micro_g(1)%cldfr)
    assNew=associated(oneMicroFields%cldfr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cldfr, oneMicroFields%cldfr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cldfr associated but oneMicroFields%cldfr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cldfr not associated and oneMicroFields%cldfr associated")
    end if

    assOld=associated(micro_g(1)%cccmp)
    assNew=associated(oneMicroFields%cccmp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cccmp, oneMicroFields%cccmp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cccmp associated but oneMicroFields%cccmp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cccmp not associated and oneMicroFields%cccmp associated")
    end if

    assOld=associated(micro_g(1)%gccmp)
    assNew=associated(oneMicroFields%gccmp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%gccmp, oneMicroFields%gccmp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%gccmp associated but oneMicroFields%gccmp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%gccmp not associated and oneMicroFields%gccmp associated")
    end if

    assOld=associated(micro_g(1)%cnm1p)
    assNew=associated(oneMicroFields%cnm1p)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cnm1p, oneMicroFields%cnm1p)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cnm1p associated but oneMicroFields%cnm1p not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cnm1p not associated and oneMicroFields%cnm1p associated")
    end if

    assOld=associated(micro_g(1)%cnm2p)
    assNew=associated(oneMicroFields%cnm2p)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cnm2p, oneMicroFields%cnm2p)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cnm2p associated but oneMicroFields%cnm2p not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cnm2p not associated and oneMicroFields%cnm2p associated")
    end if

    assOld=associated(micro_g(1)%cnm3p)
    assNew=associated(oneMicroFields%cnm3p)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cnm3p, oneMicroFields%cnm3p)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cnm3p associated but oneMicroFields%cnm3p not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cnm3p not associated and oneMicroFields%cnm3p associated")
    end if

    assOld=associated(micro_g(1)%cnm8p)
    assNew=associated(oneMicroFields%cnm8p)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cnm8p, oneMicroFields%cnm8p)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cnm8p associated but oneMicroFields%cnm8p not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cnm8p not associated and oneMicroFields%cnm8p associated")
    end if

    assOld=associated(micro_g(1)%md1np)
    assNew=associated(oneMicroFields%md1np)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%md1np, oneMicroFields%md1np)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%md1np associated but oneMicroFields%md1np not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%md1np not associated and oneMicroFields%md1np associated")
    end if

    assOld=associated(micro_g(1)%md2np)
    assNew=associated(oneMicroFields%md2np)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%md2np, oneMicroFields%md2np)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%md2np associated but oneMicroFields%md2np not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%md2np not associated and oneMicroFields%md2np associated")
    end if

    assOld=associated(micro_g(1)%salt_filmp)
    assNew=associated(oneMicroFields%salt_filmp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%salt_filmp, oneMicroFields%salt_filmp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%salt_filmp associated but oneMicroFields%salt_filmp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%salt_filmp not associated and oneMicroFields%salt_filmp associated")
    end if

    assOld=associated(micro_g(1)%salt_jetp)
    assNew=associated(oneMicroFields%salt_jetp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%salt_jetp, oneMicroFields%salt_jetp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%salt_jetp associated but oneMicroFields%salt_jetp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%salt_jetp not associated and oneMicroFields%salt_jetp associated")
    end if

    assOld=associated(micro_g(1)%salt_spmp)
    assNew=associated(oneMicroFields%salt_spmp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%salt_spmp, oneMicroFields%salt_spmp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%salt_spmp associated but oneMicroFields%salt_spmp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%salt_spmp not associated and oneMicroFields%salt_spmp associated")
    end if

    assOld=associated(micro_g(1)%pcpvr)
    assNew=associated(oneMicroFields%pcpvr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcpvr, oneMicroFields%pcpvr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvr associated but oneMicroFields%pcpvr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvr not associated and oneMicroFields%pcpvr associated")
    end if

    assOld=associated(micro_g(1)%pcpvp)
    assNew=associated(oneMicroFields%pcpvp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcpvp, oneMicroFields%pcpvp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvp associated but oneMicroFields%pcpvp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvp not associated and oneMicroFields%pcpvp associated")
    end if

    assOld=associated(micro_g(1)%pcpvs)
    assNew=associated(oneMicroFields%pcpvs)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcpvs, oneMicroFields%pcpvs)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvs associated but oneMicroFields%pcpvs not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvs not associated and oneMicroFields%pcpvs associated")
    end if

    assOld=associated(micro_g(1)%pcpva)
    assNew=associated(oneMicroFields%pcpva)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcpva, oneMicroFields%pcpva)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpva associated but oneMicroFields%pcpva not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpva not associated and oneMicroFields%pcpva associated")
    end if

    assOld=associated(micro_g(1)%pcpvg)
    assNew=associated(oneMicroFields%pcpvg)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcpvg, oneMicroFields%pcpvg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvg associated but oneMicroFields%pcpvg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvg not associated and oneMicroFields%pcpvg associated")
    end if

    assOld=associated(micro_g(1)%pcpvh)
    assNew=associated(oneMicroFields%pcpvh)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcpvh, oneMicroFields%pcpvh)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvh associated but oneMicroFields%pcpvh not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvh not associated and oneMicroFields%pcpvh associated")
    end if

    assOld=associated(micro_g(1)%pcpvd)
    assNew=associated(oneMicroFields%pcpvd)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcpvd, oneMicroFields%pcpvd)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvd associated but oneMicroFields%pcpvd not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvd not associated and oneMicroFields%pcpvd associated")
    end if

    assOld=associated(micro_g(1)%nuccldr)
    assNew=associated(oneMicroFields%nuccldr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%nuccldr, oneMicroFields%nuccldr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nuccldr associated but oneMicroFields%nuccldr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nuccldr not associated and oneMicroFields%nuccldr associated")
    end if

    assOld=associated(micro_g(1)%nuccldc)
    assNew=associated(oneMicroFields%nuccldc)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%nuccldc, oneMicroFields%nuccldc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nuccldc associated but oneMicroFields%nuccldc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nuccldc not associated and oneMicroFields%nuccldc associated")
    end if

    assOld=associated(micro_g(1)%nucicer)
    assNew=associated(oneMicroFields%nucicer)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%nucicer, oneMicroFields%nucicer)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nucicer associated but oneMicroFields%nucicer not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nucicer not associated and oneMicroFields%nucicer associated")
    end if

    assOld=associated(micro_g(1)%nucicec)
    assNew=associated(oneMicroFields%nucicec)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%nucicec, oneMicroFields%nucicec)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nucicec associated but oneMicroFields%nucicec not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nucicec not associated and oneMicroFields%nucicec associated")
    end if

    assOld=associated(micro_g(1)%inuchomr)
    assNew=associated(oneMicroFields%inuchomr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuchomr, oneMicroFields%inuchomr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchomr associated but oneMicroFields%inuchomr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchomr not associated and oneMicroFields%inuchomr associated")
    end if

    assOld=associated(micro_g(1)%inuchomc)
    assNew=associated(oneMicroFields%inuchomc)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuchomc, oneMicroFields%inuchomc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchomc associated but oneMicroFields%inuchomc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchomc not associated and oneMicroFields%inuchomc associated")
    end if

    assOld=associated(micro_g(1)%inuccontr)
    assNew=associated(oneMicroFields%inuccontr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuccontr, oneMicroFields%inuccontr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuccontr associated but oneMicroFields%inuccontr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuccontr not associated and oneMicroFields%inuccontr associated")
    end if

    assOld=associated(micro_g(1)%inuccontc)
    assNew=associated(oneMicroFields%inuccontc)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuccontc, oneMicroFields%inuccontc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuccontc associated but oneMicroFields%inuccontc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuccontc not associated and oneMicroFields%inuccontc associated")
    end if

    assOld=associated(micro_g(1)%inucifnr)
    assNew=associated(oneMicroFields%inucifnr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inucifnr, oneMicroFields%inucifnr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inucifnr associated but oneMicroFields%inucifnr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inucifnr not associated and oneMicroFields%inucifnr associated")
    end if

    assOld=associated(micro_g(1)%inucifnc)
    assNew=associated(oneMicroFields%inucifnc)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inucifnc, oneMicroFields%inucifnc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inucifnc associated but oneMicroFields%inucifnc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inucifnc not associated and oneMicroFields%inucifnc associated")
    end if

    assOld=associated(micro_g(1)%inuchazr)
    assNew=associated(oneMicroFields%inuchazr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuchazr, oneMicroFields%inuchazr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchazr associated but oneMicroFields%inuchazr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchazr not associated and oneMicroFields%inuchazr associated")
    end if

    assOld=associated(micro_g(1)%inuchazc)
    assNew=associated(oneMicroFields%inuchazc)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuchazc, oneMicroFields%inuchazc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchazc associated but oneMicroFields%inuchazc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchazc not associated and oneMicroFields%inuchazc associated")
    end if

    assOld=associated(micro_g(1)%vapliq)
    assNew=associated(oneMicroFields%vapliq)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapliq, oneMicroFields%vapliq)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapliq associated but oneMicroFields%vapliq not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapliq not associated and oneMicroFields%vapliq associated")
    end if

    assOld=associated(micro_g(1)%vapice)
    assNew=associated(oneMicroFields%vapice)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapice, oneMicroFields%vapice)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapice associated but oneMicroFields%vapice not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapice not associated and oneMicroFields%vapice associated")
    end if

    assOld=associated(micro_g(1)%vapcld)
    assNew=associated(oneMicroFields%vapcld)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapcld, oneMicroFields%vapcld)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapcld associated but oneMicroFields%vapcld not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapcld not associated and oneMicroFields%vapcld associated")
    end if

    assOld=associated(micro_g(1)%vaprain)
    assNew=associated(oneMicroFields%vaprain)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vaprain, oneMicroFields%vaprain)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vaprain associated but oneMicroFields%vaprain not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vaprain not associated and oneMicroFields%vaprain associated")
    end if

    assOld=associated(micro_g(1)%vappris)
    assNew=associated(oneMicroFields%vappris)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vappris, oneMicroFields%vappris)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vappris associated but oneMicroFields%vappris not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vappris not associated and oneMicroFields%vappris associated")
    end if

    assOld=associated(micro_g(1)%vapsnow)
    assNew=associated(oneMicroFields%vapsnow)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapsnow, oneMicroFields%vapsnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapsnow associated but oneMicroFields%vapsnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapsnow not associated and oneMicroFields%vapsnow associated")
    end if

    assOld=associated(micro_g(1)%vapaggr)
    assNew=associated(oneMicroFields%vapaggr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapaggr, oneMicroFields%vapaggr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapaggr associated but oneMicroFields%vapaggr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapaggr not associated and oneMicroFields%vapaggr associated")
    end if

    assOld=associated(micro_g(1)%vapgrau)
    assNew=associated(oneMicroFields%vapgrau)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapgrau, oneMicroFields%vapgrau)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapgrau associated but oneMicroFields%vapgrau not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapgrau not associated and oneMicroFields%vapgrau associated")
    end if

    assOld=associated(micro_g(1)%vaphail)
    assNew=associated(oneMicroFields%vaphail)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vaphail, oneMicroFields%vaphail)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vaphail associated but oneMicroFields%vaphail not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vaphail not associated and oneMicroFields%vaphail associated")
    end if

    assOld=associated(micro_g(1)%vapdriz)
    assNew=associated(oneMicroFields%vapdriz)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapdriz, oneMicroFields%vapdriz)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapdriz associated but oneMicroFields%vapdriz not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapdriz not associated and oneMicroFields%vapdriz associated")
    end if

    assOld=associated(micro_g(1)%meltice)
    assNew=associated(oneMicroFields%meltice)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%meltice, oneMicroFields%meltice)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltice associated but oneMicroFields%meltice not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltice not associated and oneMicroFields%meltice associated")
    end if

    assOld=associated(micro_g(1)%meltpris)
    assNew=associated(oneMicroFields%meltpris)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%meltpris, oneMicroFields%meltpris)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltpris associated but oneMicroFields%meltpris not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltpris not associated and oneMicroFields%meltpris associated")
    end if

    assOld=associated(micro_g(1)%meltsnow)
    assNew=associated(oneMicroFields%meltsnow)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%meltsnow, oneMicroFields%meltsnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltsnow associated but oneMicroFields%meltsnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltsnow not associated and oneMicroFields%meltsnow associated")
    end if

    assOld=associated(micro_g(1)%meltaggr)
    assNew=associated(oneMicroFields%meltaggr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%meltaggr, oneMicroFields%meltaggr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltaggr associated but oneMicroFields%meltaggr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltaggr not associated and oneMicroFields%meltaggr associated")
    end if

    assOld=associated(micro_g(1)%meltgrau)
    assNew=associated(oneMicroFields%meltgrau)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%meltgrau, oneMicroFields%meltgrau)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltgrau associated but oneMicroFields%meltgrau not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltgrau not associated and oneMicroFields%meltgrau associated")
    end if

    assOld=associated(micro_g(1)%melthail)
    assNew=associated(oneMicroFields%melthail)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%melthail, oneMicroFields%melthail)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%melthail associated but oneMicroFields%melthail not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%melthail not associated and oneMicroFields%melthail associated")
    end if

    assOld=associated(micro_g(1)%cld2rain)
    assNew=associated(oneMicroFields%cld2rain)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cld2rain, oneMicroFields%cld2rain)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cld2rain associated but oneMicroFields%cld2rain not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cld2rain not associated and oneMicroFields%cld2rain associated")
    end if

    assOld=associated(micro_g(1)%rimecld)
    assNew=associated(oneMicroFields%rimecld)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecld, oneMicroFields%rimecld)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecld associated but oneMicroFields%rimecld not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecld not associated and oneMicroFields%rimecld associated")
    end if

    assOld=associated(micro_g(1)%rimecldsnow)
    assNew=associated(oneMicroFields%rimecldsnow)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecldsnow, oneMicroFields%rimecldsnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldsnow associated but oneMicroFields%rimecldsnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldsnow not associated and oneMicroFields%rimecldsnow associated")
    end if

    assOld=associated(micro_g(1)%rimecldaggr)
    assNew=associated(oneMicroFields%rimecldaggr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecldaggr, oneMicroFields%rimecldaggr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldaggr associated but oneMicroFields%rimecldaggr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldaggr not associated and oneMicroFields%rimecldaggr associated")
    end if

    assOld=associated(micro_g(1)%rimecldgrau)
    assNew=associated(oneMicroFields%rimecldgrau)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecldgrau, oneMicroFields%rimecldgrau)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldgrau associated but oneMicroFields%rimecldgrau not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldgrau not associated and oneMicroFields%rimecldgrau associated")
    end if

    assOld=associated(micro_g(1)%rimecldhail)
    assNew=associated(oneMicroFields%rimecldhail)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecldhail, oneMicroFields%rimecldhail)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldhail associated but oneMicroFields%rimecldhail not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldhail not associated and oneMicroFields%rimecldhail associated")
    end if

    assOld=associated(micro_g(1)%rain2ice)
    assNew=associated(oneMicroFields%rain2ice)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2ice, oneMicroFields%rain2ice)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ice associated but oneMicroFields%rain2ice not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ice not associated and oneMicroFields%rain2ice associated")
    end if

    assOld=associated(micro_g(1)%rain2pr)
    assNew=associated(oneMicroFields%rain2pr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2pr, oneMicroFields%rain2pr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2pr associated but oneMicroFields%rain2pr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2pr not associated and oneMicroFields%rain2pr associated")
    end if

    assOld=associated(micro_g(1)%rain2sn)
    assNew=associated(oneMicroFields%rain2sn)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2sn, oneMicroFields%rain2sn)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2sn associated but oneMicroFields%rain2sn not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2sn not associated and oneMicroFields%rain2sn associated")
    end if

    assOld=associated(micro_g(1)%rain2ag)
    assNew=associated(oneMicroFields%rain2ag)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2ag, oneMicroFields%rain2ag)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ag associated but oneMicroFields%rain2ag not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ag not associated and oneMicroFields%rain2ag associated")
    end if

    assOld=associated(micro_g(1)%rain2gr)
    assNew=associated(oneMicroFields%rain2gr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2gr, oneMicroFields%rain2gr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2gr associated but oneMicroFields%rain2gr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2gr not associated and oneMicroFields%rain2gr associated")
    end if

    assOld=associated(micro_g(1)%rain2ha)
    assNew=associated(oneMicroFields%rain2ha)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2ha, oneMicroFields%rain2ha)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ha associated but oneMicroFields%rain2ha not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ha not associated and oneMicroFields%rain2ha associated")
    end if

    assOld=associated(micro_g(1)%rain2ha_xtra)
    assNew=associated(oneMicroFields%rain2ha_xtra)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2ha_xtra, oneMicroFields%rain2ha_xtra)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ha_xtra associated but oneMicroFields%rain2ha_xtra not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ha_xtra not associated and oneMicroFields%rain2ha_xtra associated")
    end if

    assOld=associated(micro_g(1)%ice2rain)
    assNew=associated(oneMicroFields%ice2rain)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%ice2rain, oneMicroFields%ice2rain)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%ice2rain associated but oneMicroFields%ice2rain not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%ice2rain not associated and oneMicroFields%ice2rain associated")
    end if

    assOld=associated(micro_g(1)%aggregate)
    assNew=associated(oneMicroFields%aggregate)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%aggregate, oneMicroFields%aggregate)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggregate associated but oneMicroFields%aggregate not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggregate not associated and oneMicroFields%aggregate associated")
    end if

    assOld=associated(micro_g(1)%aggrselfpris)
    assNew=associated(oneMicroFields%aggrselfpris)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%aggrselfpris, oneMicroFields%aggrselfpris)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrselfpris associated but oneMicroFields%aggrselfpris not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrselfpris not associated and oneMicroFields%aggrselfpris associated")
    end if

    assOld=associated(micro_g(1)%aggrselfsnow)
    assNew=associated(oneMicroFields%aggrselfsnow)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%aggrselfsnow, oneMicroFields%aggrselfsnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrselfsnow associated but oneMicroFields%aggrselfsnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrselfsnow not associated and oneMicroFields%aggrselfsnow associated")
    end if

    assOld=associated(micro_g(1)%aggrprissnow)
    assNew=associated(oneMicroFields%aggrprissnow)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%aggrprissnow, oneMicroFields%aggrprissnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrprissnow associated but oneMicroFields%aggrprissnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrprissnow not associated and oneMicroFields%aggrprissnow associated")
    end if

    assOld=associated(micro_g(1)%latheatvap)
    assNew=associated(oneMicroFields%latheatvap)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%latheatvap, oneMicroFields%latheatvap)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%latheatvap associated but oneMicroFields%latheatvap not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%latheatvap not associated and oneMicroFields%latheatvap associated")
    end if

    assOld=associated(micro_g(1)%latheatfrz)
    assNew=associated(oneMicroFields%latheatfrz)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%latheatfrz, oneMicroFields%latheatfrz)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%latheatfrz associated but oneMicroFields%latheatfrz not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%latheatfrz not associated and oneMicroFields%latheatfrz associated")
    end if

    assOld=associated(micro_g(1)%nuccldrt)
    assNew=associated(oneMicroFields%nuccldrt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%nuccldrt, oneMicroFields%nuccldrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nuccldrt associated but oneMicroFields%nuccldrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nuccldrt not associated and oneMicroFields%nuccldrt associated")
    end if

    assOld=associated(micro_g(1)%nuccldct)
    assNew=associated(oneMicroFields%nuccldct)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%nuccldct, oneMicroFields%nuccldct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nuccldct associated but oneMicroFields%nuccldct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nuccldct not associated and oneMicroFields%nuccldct associated")
    end if

    assOld=associated(micro_g(1)%nucicert)
    assNew=associated(oneMicroFields%nucicert)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%nucicert, oneMicroFields%nucicert)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nucicert associated but oneMicroFields%nucicert not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nucicert not associated and oneMicroFields%nucicert associated")
    end if

    assOld=associated(micro_g(1)%nucicect)
    assNew=associated(oneMicroFields%nucicect)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%nucicect, oneMicroFields%nucicect)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nucicect associated but oneMicroFields%nucicect not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nucicect not associated and oneMicroFields%nucicect associated")
    end if

    assOld=associated(micro_g(1)%inuchomrt)
    assNew=associated(oneMicroFields%inuchomrt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuchomrt, oneMicroFields%inuchomrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchomrt associated but oneMicroFields%inuchomrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchomrt not associated and oneMicroFields%inuchomrt associated")
    end if

    assOld=associated(micro_g(1)%inuchomct)
    assNew=associated(oneMicroFields%inuchomct)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuchomct, oneMicroFields%inuchomct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchomct associated but oneMicroFields%inuchomct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchomct not associated and oneMicroFields%inuchomct associated")
    end if

    assOld=associated(micro_g(1)%inuccontrt)
    assNew=associated(oneMicroFields%inuccontrt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuccontrt, oneMicroFields%inuccontrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuccontrt associated but oneMicroFields%inuccontrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuccontrt not associated and oneMicroFields%inuccontrt associated")
    end if

    assOld=associated(micro_g(1)%inuccontct)
    assNew=associated(oneMicroFields%inuccontct)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuccontct, oneMicroFields%inuccontct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuccontct associated but oneMicroFields%inuccontct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuccontct not associated and oneMicroFields%inuccontct associated")
    end if

    assOld=associated(micro_g(1)%inucifnrt)
    assNew=associated(oneMicroFields%inucifnrt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inucifnrt, oneMicroFields%inucifnrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inucifnrt associated but oneMicroFields%inucifnrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inucifnrt not associated and oneMicroFields%inucifnrt associated")
    end if

    assOld=associated(micro_g(1)%inucifnct)
    assNew=associated(oneMicroFields%inucifnct)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inucifnct, oneMicroFields%inucifnct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inucifnct associated but oneMicroFields%inucifnct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inucifnct not associated and oneMicroFields%inucifnct associated")
    end if

    assOld=associated(micro_g(1)%inuchazrt)
    assNew=associated(oneMicroFields%inuchazrt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuchazrt, oneMicroFields%inuchazrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchazrt associated but oneMicroFields%inuchazrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchazrt not associated and oneMicroFields%inuchazrt associated")
    end if

    assOld=associated(micro_g(1)%inuchazct)
    assNew=associated(oneMicroFields%inuchazct)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%inuchazct, oneMicroFields%inuchazct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchazct associated but oneMicroFields%inuchazct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchazct not associated and oneMicroFields%inuchazct associated")
    end if

    assOld=associated(micro_g(1)%vapliqt)
    assNew=associated(oneMicroFields%vapliqt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapliqt, oneMicroFields%vapliqt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapliqt associated but oneMicroFields%vapliqt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapliqt not associated and oneMicroFields%vapliqt associated")
    end if

    assOld=associated(micro_g(1)%vapicet)
    assNew=associated(oneMicroFields%vapicet)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapicet, oneMicroFields%vapicet)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapicet associated but oneMicroFields%vapicet not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapicet not associated and oneMicroFields%vapicet associated")
    end if

    assOld=associated(micro_g(1)%vapcldt)
    assNew=associated(oneMicroFields%vapcldt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapcldt, oneMicroFields%vapcldt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapcldt associated but oneMicroFields%vapcldt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapcldt not associated and oneMicroFields%vapcldt associated")
    end if

    assOld=associated(micro_g(1)%vapraint)
    assNew=associated(oneMicroFields%vapraint)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapraint, oneMicroFields%vapraint)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapraint associated but oneMicroFields%vapraint not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapraint not associated and oneMicroFields%vapraint associated")
    end if

    assOld=associated(micro_g(1)%vapprist)
    assNew=associated(oneMicroFields%vapprist)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapprist, oneMicroFields%vapprist)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapprist associated but oneMicroFields%vapprist not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapprist not associated and oneMicroFields%vapprist associated")
    end if

    assOld=associated(micro_g(1)%vapsnowt)
    assNew=associated(oneMicroFields%vapsnowt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapsnowt, oneMicroFields%vapsnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapsnowt associated but oneMicroFields%vapsnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapsnowt not associated and oneMicroFields%vapsnowt associated")
    end if

    assOld=associated(micro_g(1)%vapaggrt)
    assNew=associated(oneMicroFields%vapaggrt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapaggrt, oneMicroFields%vapaggrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapaggrt associated but oneMicroFields%vapaggrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapaggrt not associated and oneMicroFields%vapaggrt associated")
    end if

    assOld=associated(micro_g(1)%vapgraut)
    assNew=associated(oneMicroFields%vapgraut)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapgraut, oneMicroFields%vapgraut)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapgraut associated but oneMicroFields%vapgraut not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapgraut not associated and oneMicroFields%vapgraut associated")
    end if

    assOld=associated(micro_g(1)%vaphailt)
    assNew=associated(oneMicroFields%vaphailt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vaphailt, oneMicroFields%vaphailt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vaphailt associated but oneMicroFields%vaphailt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vaphailt not associated and oneMicroFields%vaphailt associated")
    end if

    assOld=associated(micro_g(1)%vapdrizt)
    assNew=associated(oneMicroFields%vapdrizt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%vapdrizt, oneMicroFields%vapdrizt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapdrizt associated but oneMicroFields%vapdrizt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapdrizt not associated and oneMicroFields%vapdrizt associated")
    end if

    assOld=associated(micro_g(1)%melticet)
    assNew=associated(oneMicroFields%melticet)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%melticet, oneMicroFields%melticet)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%melticet associated but oneMicroFields%melticet not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%melticet not associated and oneMicroFields%melticet associated")
    end if

    assOld=associated(micro_g(1)%meltprist)
    assNew=associated(oneMicroFields%meltprist)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%meltprist, oneMicroFields%meltprist)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltprist associated but oneMicroFields%meltprist not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltprist not associated and oneMicroFields%meltprist associated")
    end if

    assOld=associated(micro_g(1)%meltsnowt)
    assNew=associated(oneMicroFields%meltsnowt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%meltsnowt, oneMicroFields%meltsnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltsnowt associated but oneMicroFields%meltsnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltsnowt not associated and oneMicroFields%meltsnowt associated")
    end if

    assOld=associated(micro_g(1)%meltaggrt)
    assNew=associated(oneMicroFields%meltaggrt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%meltaggrt, oneMicroFields%meltaggrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltaggrt associated but oneMicroFields%meltaggrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltaggrt not associated and oneMicroFields%meltaggrt associated")
    end if

    assOld=associated(micro_g(1)%meltgraut)
    assNew=associated(oneMicroFields%meltgraut)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%meltgraut, oneMicroFields%meltgraut)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltgraut associated but oneMicroFields%meltgraut not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltgraut not associated and oneMicroFields%meltgraut associated")
    end if

    assOld=associated(micro_g(1)%melthailt)
    assNew=associated(oneMicroFields%melthailt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%melthailt, oneMicroFields%melthailt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%melthailt associated but oneMicroFields%melthailt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%melthailt not associated and oneMicroFields%melthailt associated")
    end if

    assOld=associated(micro_g(1)%cld2raint)
    assNew=associated(oneMicroFields%cld2raint)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%cld2raint, oneMicroFields%cld2raint)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cld2raint associated but oneMicroFields%cld2raint not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cld2raint not associated and oneMicroFields%cld2raint associated")
    end if

    assOld=associated(micro_g(1)%rimecldt)
    assNew=associated(oneMicroFields%rimecldt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecldt, oneMicroFields%rimecldt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldt associated but oneMicroFields%rimecldt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldt not associated and oneMicroFields%rimecldt associated")
    end if

    assOld=associated(micro_g(1)%rimecldsnowt)
    assNew=associated(oneMicroFields%rimecldsnowt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecldsnowt, oneMicroFields%rimecldsnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldsnowt associated but oneMicroFields%rimecldsnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldsnowt not associated and oneMicroFields%rimecldsnowt associated")
    end if

    assOld=associated(micro_g(1)%rimecldaggrt)
    assNew=associated(oneMicroFields%rimecldaggrt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecldaggrt, oneMicroFields%rimecldaggrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldaggrt associated but oneMicroFields%rimecldaggrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldaggrt not associated and oneMicroFields%rimecldaggrt associated")
    end if

    assOld=associated(micro_g(1)%rimecldgraut)
    assNew=associated(oneMicroFields%rimecldgraut)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecldgraut, oneMicroFields%rimecldgraut)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldgraut associated but oneMicroFields%rimecldgraut not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldgraut not associated and oneMicroFields%rimecldgraut associated")
    end if

    assOld=associated(micro_g(1)%rimecldhailt)
    assNew=associated(oneMicroFields%rimecldhailt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rimecldhailt, oneMicroFields%rimecldhailt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldhailt associated but oneMicroFields%rimecldhailt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldhailt not associated and oneMicroFields%rimecldhailt associated")
    end if

    assOld=associated(micro_g(1)%rain2icet)
    assNew=associated(oneMicroFields%rain2icet)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2icet, oneMicroFields%rain2icet)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2icet associated but oneMicroFields%rain2icet not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2icet not associated and oneMicroFields%rain2icet associated")
    end if

    assOld=associated(micro_g(1)%rain2prt)
    assNew=associated(oneMicroFields%rain2prt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2prt, oneMicroFields%rain2prt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2prt associated but oneMicroFields%rain2prt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2prt not associated and oneMicroFields%rain2prt associated")
    end if

    assOld=associated(micro_g(1)%rain2snt)
    assNew=associated(oneMicroFields%rain2snt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2snt, oneMicroFields%rain2snt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2snt associated but oneMicroFields%rain2snt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2snt not associated and oneMicroFields%rain2snt associated")
    end if

    assOld=associated(micro_g(1)%rain2agt)
    assNew=associated(oneMicroFields%rain2agt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2agt, oneMicroFields%rain2agt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2agt associated but oneMicroFields%rain2agt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2agt not associated and oneMicroFields%rain2agt associated")
    end if

    assOld=associated(micro_g(1)%rain2grt)
    assNew=associated(oneMicroFields%rain2grt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2grt, oneMicroFields%rain2grt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2grt associated but oneMicroFields%rain2grt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2grt not associated and oneMicroFields%rain2grt associated")
    end if

    assOld=associated(micro_g(1)%rain2hat)
    assNew=associated(oneMicroFields%rain2hat)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2hat, oneMicroFields%rain2hat)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2hat associated but oneMicroFields%rain2hat not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2hat not associated and oneMicroFields%rain2hat associated")
    end if

    assOld=associated(micro_g(1)%rain2ha_xtrat)
    assNew=associated(oneMicroFields%rain2ha_xtrat)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%rain2ha_xtrat, oneMicroFields%rain2ha_xtrat)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ha_xtrat associated but oneMicroFields%rain2ha_xtrat not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ha_xtrat not associated and oneMicroFields%rain2ha_xtrat associated")
    end if

    assOld=associated(micro_g(1)%ice2raint)
    assNew=associated(oneMicroFields%ice2raint)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%ice2raint, oneMicroFields%ice2raint)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%ice2raint associated but oneMicroFields%ice2raint not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%ice2raint not associated and oneMicroFields%ice2raint associated")
    end if

    assOld=associated(micro_g(1)%aggregatet)
    assNew=associated(oneMicroFields%aggregatet)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%aggregatet, oneMicroFields%aggregatet)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggregatet associated but oneMicroFields%aggregatet not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggregatet not associated and oneMicroFields%aggregatet associated")
    end if

    assOld=associated(micro_g(1)%aggrselfprist)
    assNew=associated(oneMicroFields%aggrselfprist)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%aggrselfprist, oneMicroFields%aggrselfprist)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrselfprist associated but oneMicroFields%aggrselfprist not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrselfprist not associated and oneMicroFields%aggrselfprist associated")
    end if

    assOld=associated(micro_g(1)%aggrselfsnowt)
    assNew=associated(oneMicroFields%aggrselfsnowt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%aggrselfsnowt, oneMicroFields%aggrselfsnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrselfsnowt associated but oneMicroFields%aggrselfsnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrselfsnowt not associated and oneMicroFields%aggrselfsnowt associated")
    end if

    assOld=associated(micro_g(1)%aggrprissnowt)
    assNew=associated(oneMicroFields%aggrprissnowt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%aggrprissnowt, oneMicroFields%aggrprissnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrprissnowt associated but oneMicroFields%aggrprissnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrprissnowt not associated and oneMicroFields%aggrprissnowt associated")
    end if

    assOld=associated(micro_g(1)%latheatvapt)
    assNew=associated(oneMicroFields%latheatvapt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%latheatvapt, oneMicroFields%latheatvapt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%latheatvapt associated but oneMicroFields%latheatvapt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%latheatvapt not associated and oneMicroFields%latheatvapt associated")
    end if

    assOld=associated(micro_g(1)%latheatfrzt)
    assNew=associated(oneMicroFields%latheatfrzt)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%latheatfrzt, oneMicroFields%latheatfrzt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%latheatfrzt associated but oneMicroFields%latheatfrzt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%latheatfrzt not associated and oneMicroFields%latheatfrzt associated")
    end if

    assOld=associated(micro_g(1)%accpr)
    assNew=associated(oneMicroFields%accpr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%accpr, oneMicroFields%accpr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpr associated but oneMicroFields%accpr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpr not associated and oneMicroFields%accpr associated")
    end if

    assOld=associated(micro_g(1)%accpp)
    assNew=associated(oneMicroFields%accpp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%accpp, oneMicroFields%accpp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpp associated but oneMicroFields%accpp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpp not associated and oneMicroFields%accpp associated")
    end if

    assOld=associated(micro_g(1)%accps)
    assNew=associated(oneMicroFields%accps)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%accps, oneMicroFields%accps)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accps associated but oneMicroFields%accps not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accps not associated and oneMicroFields%accps associated")
    end if

    assOld=associated(micro_g(1)%accpa)
    assNew=associated(oneMicroFields%accpa)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%accpa, oneMicroFields%accpa)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpa associated but oneMicroFields%accpa not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpa not associated and oneMicroFields%accpa associated")
    end if

    assOld=associated(micro_g(1)%accpg)
    assNew=associated(oneMicroFields%accpg)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%accpg, oneMicroFields%accpg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpg associated but oneMicroFields%accpg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpg not associated and oneMicroFields%accpg associated")
    end if

    assOld=associated(micro_g(1)%accph)
    assNew=associated(oneMicroFields%accph)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%accph, oneMicroFields%accph)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accph associated but oneMicroFields%accph not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accph not associated and oneMicroFields%accph associated")
    end if

    assOld=associated(micro_g(1)%accpd)
    assNew=associated(oneMicroFields%accpd)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%accpd, oneMicroFields%accpd)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpd associated but oneMicroFields%accpd not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpd not associated and oneMicroFields%accpd associated")
    end if

    assOld=associated(micro_g(1)%pcprr)
    assNew=associated(oneMicroFields%pcprr)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcprr, oneMicroFields%pcprr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprr associated but oneMicroFields%pcprr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprr not associated and oneMicroFields%pcprr associated")
    end if

    assOld=associated(micro_g(1)%pcprp)
    assNew=associated(oneMicroFields%pcprp)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcprp, oneMicroFields%pcprp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprp associated but oneMicroFields%pcprp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprp not associated and oneMicroFields%pcprp associated")
    end if

    assOld=associated(micro_g(1)%pcprs)
    assNew=associated(oneMicroFields%pcprs)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcprs, oneMicroFields%pcprs)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprs associated but oneMicroFields%pcprs not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprs not associated and oneMicroFields%pcprs associated")
    end if

    assOld=associated(micro_g(1)%pcpra)
    assNew=associated(oneMicroFields%pcpra)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcpra, oneMicroFields%pcpra)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpra associated but oneMicroFields%pcpra not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpra not associated and oneMicroFields%pcpra associated")
    end if

    assOld=associated(micro_g(1)%pcprg)
    assNew=associated(oneMicroFields%pcprg)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcprg, oneMicroFields%pcprg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprg associated but oneMicroFields%pcprg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprg not associated and oneMicroFields%pcprg associated")
    end if

    assOld=associated(micro_g(1)%pcprh)
    assNew=associated(oneMicroFields%pcprh)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcprh, oneMicroFields%pcprh)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprh associated but oneMicroFields%pcprh not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprh not associated and oneMicroFields%pcprh associated")
    end if

    assOld=associated(micro_g(1)%pcprd)
    assNew=associated(oneMicroFields%pcprd)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcprd, oneMicroFields%pcprd)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprd associated but oneMicroFields%pcprd not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprd not associated and oneMicroFields%pcprd associated")
    end if

    assOld=associated(micro_g(1)%pcpg)
    assNew=associated(oneMicroFields%pcpg)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%pcpg, oneMicroFields%pcpg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpg associated but oneMicroFields%pcpg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpg not associated and oneMicroFields%pcpg associated")
    end if

    assOld=associated(micro_g(1)%qpcpg)
    assNew=associated(oneMicroFields%qpcpg)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%qpcpg, oneMicroFields%qpcpg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%qpcpg associated but oneMicroFields%qpcpg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%qpcpg not associated and oneMicroFields%qpcpg associated")
    end if

    assOld=associated(micro_g(1)%dpcpg)
    assNew=associated(oneMicroFields%dpcpg)
    if (assOld .and. assNew) then
       call ToCopy(micro_g(1)%dpcpg, oneMicroFields%dpcpg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%dpcpg associated but oneMicroFields%dpcpg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%dpcpg not associated and oneMicroFields%dpcpg associated")
    end if

  end subroutine DeepCopyToMicroFields




  subroutine DeepCopyFromMicroFields(oneMicroFields, from)
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    character(len=*), intent(in) :: from

    logical :: assOld, assNew
    character(len=*), parameter :: h="**(DeepCopyFromMicroFields)**"

    if (copyFrom /= "") then
       call fatal_error(h//" invoked from "//trim(adjustl(from))//&
            " just after invoked from "//trim(adjustl(copyFrom)))
    end if

    copyTo=""
    copyFrom=from

    assOld=associated(micro_g(1)%rcp)
    assNew=associated(oneMicroFields%rcp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rcp, micro_g(1)%rcp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rcp associated but oneMicroFields%rcp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rcp not associated and oneMicroFields%rcp associated")
    end if

    assOld=associated(micro_g(1)%rdp)
    assNew=associated(oneMicroFields%rdp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rdp, micro_g(1)%rdp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rdp associated but oneMicroFields%rdp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rdp not associated and oneMicroFields%rdp associated")
    end if

    assOld=associated(micro_g(1)%rrp)
    assNew=associated(oneMicroFields%rrp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rrp, micro_g(1)%rrp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rrp associated but oneMicroFields%rrp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rrp not associated and oneMicroFields%rrp associated")
    end if

    assOld=associated(micro_g(1)%rpp)
    assNew=associated(oneMicroFields%rpp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rpp, micro_g(1)%rpp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rpp associated but oneMicroFields%rpp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rpp not associated and oneMicroFields%rpp associated")
    end if

    assOld=associated(micro_g(1)%rsp)
    assNew=associated(oneMicroFields%rsp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rsp, micro_g(1)%rsp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rsp associated but oneMicroFields%rsp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rsp not associated and oneMicroFields%rsp associated")
    end if
    assOld=associated(micro_g(1)%rap)
    assNew=associated(oneMicroFields%rap)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rap, micro_g(1)%rap)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rap associated but oneMicroFields%rap not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rap not associated and oneMicroFields%rap associated")
    end if

    assOld=associated(micro_g(1)%rgp)
    assNew=associated(oneMicroFields%rgp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rgp, micro_g(1)%rgp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rgp associated but oneMicroFields%rgp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rgp not associated and oneMicroFields%rgp associated")
    end if

    assOld=associated(micro_g(1)%rhp)
    assNew=associated(oneMicroFields%rhp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rhp, micro_g(1)%rhp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rhp associated but oneMicroFields%rhp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rhp not associated and oneMicroFields%rhp associated")
    end if

    assOld=associated(micro_g(1)%ccp)
    assNew=associated(oneMicroFields%ccp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%ccp, micro_g(1)%ccp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%ccp associated but oneMicroFields%ccp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%ccp not associated and oneMicroFields%ccp associated")
    end if

    assOld=associated(micro_g(1)%cdp)
    assNew=associated(oneMicroFields%cdp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cdp, micro_g(1)%cdp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cdp associated but oneMicroFields%cdp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cdp not associated and oneMicroFields%cdp associated")
    end if

    assOld=associated(micro_g(1)%crp)
    assNew=associated(oneMicroFields%crp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%crp, micro_g(1)%crp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%crp associated but oneMicroFields%crp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%crp not associated and oneMicroFields%crp associated")
    end if

    assOld=associated(micro_g(1)%cpp)
    assNew=associated(oneMicroFields%cpp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cpp, micro_g(1)%cpp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cpp associated but oneMicroFields%cpp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cpp not associated and oneMicroFields%cpp associated")
    end if

    assOld=associated(micro_g(1)%csp)
    assNew=associated(oneMicroFields%csp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%csp, micro_g(1)%csp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%csp associated but oneMicroFields%csp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%csp not associated and oneMicroFields%csp associated")
    end if

    assOld=associated(micro_g(1)%cap)
    assNew=associated(oneMicroFields%cap)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cap, micro_g(1)%cap)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cap associated but oneMicroFields%cap not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cap not associated and oneMicroFields%cap associated")
    end if

    assOld=associated(micro_g(1)%cgp)
    assNew=associated(oneMicroFields%cgp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cgp, micro_g(1)%cgp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cgp associated but oneMicroFields%cgp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cgp not associated and oneMicroFields%cgp associated")
    end if

    assOld=associated(micro_g(1)%chp)
    assNew=associated(oneMicroFields%chp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%chp, micro_g(1)%chp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%chp associated but oneMicroFields%chp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%chp not associated and oneMicroFields%chp associated")
    end if

    assOld=associated(micro_g(1)%cccnp)
    assNew=associated(oneMicroFields%cccnp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cccnp, micro_g(1)%cccnp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cccnp associated but oneMicroFields%cccnp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cccnp not associated and oneMicroFields%cccnp associated")
    end if

    assOld=associated(micro_g(1)%gccnp)
    assNew=associated(oneMicroFields%gccnp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%gccnp, micro_g(1)%gccnp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%gccnp associated but oneMicroFields%gccnp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%gccnp not associated and oneMicroFields%gccnp associated")
    end if

    assOld=associated(micro_g(1)%cifnp)
    assNew=associated(oneMicroFields%cifnp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cifnp, micro_g(1)%cifnp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cifnp associated but oneMicroFields%cifnp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cifnp not associated and oneMicroFields%cifnp associated")
    end if

    assOld=associated(micro_g(1)%q2)
    assNew=associated(oneMicroFields%q2)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%q2, micro_g(1)%q2)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%q2 associated but oneMicroFields%q2 not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%q2 not associated and oneMicroFields%q2 associated")
    end if

    assOld=associated(micro_g(1)%q6)
    assNew=associated(oneMicroFields%q6)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%q6, micro_g(1)%q6)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%q6 associated but oneMicroFields%q6 not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%q6 not associated and oneMicroFields%q6 associated")
    end if

    assOld=associated(micro_g(1)%q7)
    assNew=associated(oneMicroFields%q7)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%q7, micro_g(1)%q7)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%q7 associated but oneMicroFields%q7 not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%q7 not associated and oneMicroFields%q7 associated")
    end if

    assOld=associated(micro_g(1)%rei)
    assNew=associated(oneMicroFields%rei)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rei, micro_g(1)%rei)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rei associated but oneMicroFields%rei not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rei not associated and oneMicroFields%rei associated")
    end if

    assOld=associated(micro_g(1)%rel)
    assNew=associated(oneMicroFields%rel)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rel, micro_g(1)%rel)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rel associated but oneMicroFields%rel not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rel not associated and oneMicroFields%rel associated")
    end if

    assOld=associated(micro_g(1)%cldfr)
    assNew=associated(oneMicroFields%cldfr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cldfr, micro_g(1)%cldfr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cldfr associated but oneMicroFields%cldfr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cldfr not associated and oneMicroFields%cldfr associated")
    end if

    assOld=associated(micro_g(1)%cccmp)
    assNew=associated(oneMicroFields%cccmp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cccmp, micro_g(1)%cccmp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cccmp associated but oneMicroFields%cccmp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cccmp not associated and oneMicroFields%cccmp associated")
    end if

    assOld=associated(micro_g(1)%gccmp)
    assNew=associated(oneMicroFields%gccmp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%gccmp, micro_g(1)%gccmp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%gccmp associated but oneMicroFields%gccmp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%gccmp not associated and oneMicroFields%gccmp associated")
    end if

    assOld=associated(micro_g(1)%cnm1p)
    assNew=associated(oneMicroFields%cnm1p)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cnm1p, micro_g(1)%cnm1p)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cnm1p associated but oneMicroFields%cnm1p not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cnm1p not associated and oneMicroFields%cnm1p associated")
    end if

    assOld=associated(micro_g(1)%cnm2p)
    assNew=associated(oneMicroFields%cnm2p)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cnm2p, micro_g(1)%cnm2p)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cnm2p associated but oneMicroFields%cnm2p not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cnm2p not associated and oneMicroFields%cnm2p associated")
    end if

    assOld=associated(micro_g(1)%cnm3p)
    assNew=associated(oneMicroFields%cnm3p)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cnm3p, micro_g(1)%cnm3p)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cnm3p associated but oneMicroFields%cnm3p not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cnm3p not associated and oneMicroFields%cnm3p associated")
    end if

    assOld=associated(micro_g(1)%cnm8p)
    assNew=associated(oneMicroFields%cnm8p)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cnm8p, micro_g(1)%cnm8p)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cnm8p associated but oneMicroFields%cnm8p not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cnm8p not associated and oneMicroFields%cnm8p associated")
    end if

    assOld=associated(micro_g(1)%md1np)
    assNew=associated(oneMicroFields%md1np)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%md1np, micro_g(1)%md1np)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%md1np associated but oneMicroFields%md1np not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%md1np not associated and oneMicroFields%md1np associated")
    end if

    assOld=associated(micro_g(1)%md2np)
    assNew=associated(oneMicroFields%md2np)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%md2np, micro_g(1)%md2np)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%md2np associated but oneMicroFields%md2np not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%md2np not associated and oneMicroFields%md2np associated")
    end if

    assOld=associated(micro_g(1)%salt_filmp)
    assNew=associated(oneMicroFields%salt_filmp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%salt_filmp, micro_g(1)%salt_filmp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%salt_filmp associated but oneMicroFields%salt_filmp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%salt_filmp not associated and oneMicroFields%salt_filmp associated")
    end if

    assOld=associated(micro_g(1)%salt_jetp)
    assNew=associated(oneMicroFields%salt_jetp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%salt_jetp, micro_g(1)%salt_jetp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%salt_jetp associated but oneMicroFields%salt_jetp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%salt_jetp not associated and oneMicroFields%salt_jetp associated")
    end if

    assOld=associated(micro_g(1)%salt_spmp)
    assNew=associated(oneMicroFields%salt_spmp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%salt_spmp, micro_g(1)%salt_spmp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%salt_spmp associated but oneMicroFields%salt_spmp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%salt_spmp not associated and oneMicroFields%salt_spmp associated")
    end if

    assOld=associated(micro_g(1)%pcpvr)
    assNew=associated(oneMicroFields%pcpvr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcpvr, micro_g(1)%pcpvr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvr associated but oneMicroFields%pcpvr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvr not associated and oneMicroFields%pcpvr associated")
    end if

    assOld=associated(micro_g(1)%pcpvp)
    assNew=associated(oneMicroFields%pcpvp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcpvp, micro_g(1)%pcpvp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvp associated but oneMicroFields%pcpvp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvp not associated and oneMicroFields%pcpvp associated")
    end if

    assOld=associated(micro_g(1)%pcpvs)
    assNew=associated(oneMicroFields%pcpvs)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcpvs, micro_g(1)%pcpvs)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvs associated but oneMicroFields%pcpvs not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvs not associated and oneMicroFields%pcpvs associated")
    end if

    assOld=associated(micro_g(1)%pcpva)
    assNew=associated(oneMicroFields%pcpva)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcpva, micro_g(1)%pcpva)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpva associated but oneMicroFields%pcpva not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpva not associated and oneMicroFields%pcpva associated")
    end if

    assOld=associated(micro_g(1)%pcpvg)
    assNew=associated(oneMicroFields%pcpvg)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcpvg, micro_g(1)%pcpvg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvg associated but oneMicroFields%pcpvg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvg not associated and oneMicroFields%pcpvg associated")
    end if

    assOld=associated(micro_g(1)%pcpvh)
    assNew=associated(oneMicroFields%pcpvh)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcpvh, micro_g(1)%pcpvh)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvh associated but oneMicroFields%pcpvh not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvh not associated and oneMicroFields%pcpvh associated")
    end if

    assOld=associated(micro_g(1)%pcpvd)
    assNew=associated(oneMicroFields%pcpvd)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcpvd, micro_g(1)%pcpvd)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpvd associated but oneMicroFields%pcpvd not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpvd not associated and oneMicroFields%pcpvd associated")
    end if

    assOld=associated(micro_g(1)%nuccldr)
    assNew=associated(oneMicroFields%nuccldr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%nuccldr, micro_g(1)%nuccldr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nuccldr associated but oneMicroFields%nuccldr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nuccldr not associated and oneMicroFields%nuccldr associated")
    end if

    assOld=associated(micro_g(1)%nuccldc)
    assNew=associated(oneMicroFields%nuccldc)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%nuccldc, micro_g(1)%nuccldc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nuccldc associated but oneMicroFields%nuccldc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nuccldc not associated and oneMicroFields%nuccldc associated")
    end if

    assOld=associated(micro_g(1)%nucicer)
    assNew=associated(oneMicroFields%nucicer)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%nucicer, micro_g(1)%nucicer)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nucicer associated but oneMicroFields%nucicer not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nucicer not associated and oneMicroFields%nucicer associated")
    end if

    assOld=associated(micro_g(1)%nucicec)
    assNew=associated(oneMicroFields%nucicec)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%nucicec, micro_g(1)%nucicec)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nucicec associated but oneMicroFields%nucicec not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nucicec not associated and oneMicroFields%nucicec associated")
    end if

    assOld=associated(micro_g(1)%inuchomr)
    assNew=associated(oneMicroFields%inuchomr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuchomr, micro_g(1)%inuchomr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchomr associated but oneMicroFields%inuchomr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchomr not associated and oneMicroFields%inuchomr associated")
    end if

    assOld=associated(micro_g(1)%inuchomc)
    assNew=associated(oneMicroFields%inuchomc)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuchomc, micro_g(1)%inuchomc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchomc associated but oneMicroFields%inuchomc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchomc not associated and oneMicroFields%inuchomc associated")
    end if

    assOld=associated(micro_g(1)%inuccontr)
    assNew=associated(oneMicroFields%inuccontr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuccontr, micro_g(1)%inuccontr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuccontr associated but oneMicroFields%inuccontr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuccontr not associated and oneMicroFields%inuccontr associated")
    end if

    assOld=associated(micro_g(1)%inuccontc)
    assNew=associated(oneMicroFields%inuccontc)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuccontc, micro_g(1)%inuccontc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuccontc associated but oneMicroFields%inuccontc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuccontc not associated and oneMicroFields%inuccontc associated")
    end if

    assOld=associated(micro_g(1)%inucifnr)
    assNew=associated(oneMicroFields%inucifnr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inucifnr, micro_g(1)%inucifnr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inucifnr associated but oneMicroFields%inucifnr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inucifnr not associated and oneMicroFields%inucifnr associated")
    end if

    assOld=associated(micro_g(1)%inucifnc)
    assNew=associated(oneMicroFields%inucifnc)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inucifnc, micro_g(1)%inucifnc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inucifnc associated but oneMicroFields%inucifnc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inucifnc not associated and oneMicroFields%inucifnc associated")
    end if

    assOld=associated(micro_g(1)%inuchazr)
    assNew=associated(oneMicroFields%inuchazr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuchazr, micro_g(1)%inuchazr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchazr associated but oneMicroFields%inuchazr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchazr not associated and oneMicroFields%inuchazr associated")
    end if

    assOld=associated(micro_g(1)%inuchazc)
    assNew=associated(oneMicroFields%inuchazc)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuchazc, micro_g(1)%inuchazc)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchazc associated but oneMicroFields%inuchazc not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchazc not associated and oneMicroFields%inuchazc associated")
    end if

    assOld=associated(micro_g(1)%vapliq)
    assNew=associated(oneMicroFields%vapliq)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapliq, micro_g(1)%vapliq)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapliq associated but oneMicroFields%vapliq not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapliq not associated and oneMicroFields%vapliq associated")
    end if

    assOld=associated(micro_g(1)%vapice)
    assNew=associated(oneMicroFields%vapice)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapice, micro_g(1)%vapice)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapice associated but oneMicroFields%vapice not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapice not associated and oneMicroFields%vapice associated")
    end if

    assOld=associated(micro_g(1)%vapcld)
    assNew=associated(oneMicroFields%vapcld)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapcld, micro_g(1)%vapcld)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapcld associated but oneMicroFields%vapcld not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapcld not associated and oneMicroFields%vapcld associated")
    end if

    assOld=associated(micro_g(1)%vaprain)
    assNew=associated(oneMicroFields%vaprain)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vaprain, micro_g(1)%vaprain)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vaprain associated but oneMicroFields%vaprain not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vaprain not associated and oneMicroFields%vaprain associated")
    end if

    assOld=associated(micro_g(1)%vappris)
    assNew=associated(oneMicroFields%vappris)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vappris, micro_g(1)%vappris)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vappris associated but oneMicroFields%vappris not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vappris not associated and oneMicroFields%vappris associated")
    end if

    assOld=associated(micro_g(1)%vapsnow)
    assNew=associated(oneMicroFields%vapsnow)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapsnow, micro_g(1)%vapsnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapsnow associated but oneMicroFields%vapsnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapsnow not associated and oneMicroFields%vapsnow associated")
    end if

    assOld=associated(micro_g(1)%vapaggr)
    assNew=associated(oneMicroFields%vapaggr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapaggr, micro_g(1)%vapaggr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapaggr associated but oneMicroFields%vapaggr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapaggr not associated and oneMicroFields%vapaggr associated")
    end if

    assOld=associated(micro_g(1)%vapgrau)
    assNew=associated(oneMicroFields%vapgrau)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapgrau, micro_g(1)%vapgrau)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapgrau associated but oneMicroFields%vapgrau not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapgrau not associated and oneMicroFields%vapgrau associated")
    end if

    assOld=associated(micro_g(1)%vaphail)
    assNew=associated(oneMicroFields%vaphail)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vaphail, micro_g(1)%vaphail)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vaphail associated but oneMicroFields%vaphail not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vaphail not associated and oneMicroFields%vaphail associated")
    end if

    assOld=associated(micro_g(1)%vapdriz)
    assNew=associated(oneMicroFields%vapdriz)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapdriz, micro_g(1)%vapdriz)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapdriz associated but oneMicroFields%vapdriz not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapdriz not associated and oneMicroFields%vapdriz associated")
    end if

    assOld=associated(micro_g(1)%meltice)
    assNew=associated(oneMicroFields%meltice)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%meltice, micro_g(1)%meltice)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltice associated but oneMicroFields%meltice not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltice not associated and oneMicroFields%meltice associated")
    end if

    assOld=associated(micro_g(1)%meltpris)
    assNew=associated(oneMicroFields%meltpris)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%meltpris, micro_g(1)%meltpris)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltpris associated but oneMicroFields%meltpris not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltpris not associated and oneMicroFields%meltpris associated")
    end if

    assOld=associated(micro_g(1)%meltsnow)
    assNew=associated(oneMicroFields%meltsnow)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%meltsnow, micro_g(1)%meltsnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltsnow associated but oneMicroFields%meltsnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltsnow not associated and oneMicroFields%meltsnow associated")
    end if

    assOld=associated(micro_g(1)%meltaggr)
    assNew=associated(oneMicroFields%meltaggr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%meltaggr, micro_g(1)%meltaggr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltaggr associated but oneMicroFields%meltaggr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltaggr not associated and oneMicroFields%meltaggr associated")
    end if

    assOld=associated(micro_g(1)%meltgrau)
    assNew=associated(oneMicroFields%meltgrau)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%meltgrau, micro_g(1)%meltgrau)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltgrau associated but oneMicroFields%meltgrau not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltgrau not associated and oneMicroFields%meltgrau associated")
    end if

    assOld=associated(micro_g(1)%melthail)
    assNew=associated(oneMicroFields%melthail)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%melthail, micro_g(1)%melthail)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%melthail associated but oneMicroFields%melthail not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%melthail not associated and oneMicroFields%melthail associated")
    end if

    assOld=associated(micro_g(1)%cld2rain)
    assNew=associated(oneMicroFields%cld2rain)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cld2rain, micro_g(1)%cld2rain)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cld2rain associated but oneMicroFields%cld2rain not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cld2rain not associated and oneMicroFields%cld2rain associated")
    end if

    assOld=associated(micro_g(1)%rimecld)
    assNew=associated(oneMicroFields%rimecld)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecld, micro_g(1)%rimecld)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecld associated but oneMicroFields%rimecld not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecld not associated and oneMicroFields%rimecld associated")
    end if

    assOld=associated(micro_g(1)%rimecldsnow)
    assNew=associated(oneMicroFields%rimecldsnow)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecldsnow, micro_g(1)%rimecldsnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldsnow associated but oneMicroFields%rimecldsnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldsnow not associated and oneMicroFields%rimecldsnow associated")
    end if

    assOld=associated(micro_g(1)%rimecldaggr)
    assNew=associated(oneMicroFields%rimecldaggr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecldaggr, micro_g(1)%rimecldaggr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldaggr associated but oneMicroFields%rimecldaggr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldaggr not associated and oneMicroFields%rimecldaggr associated")
    end if

    assOld=associated(micro_g(1)%rimecldgrau)
    assNew=associated(oneMicroFields%rimecldgrau)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecldgrau, micro_g(1)%rimecldgrau)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldgrau associated but oneMicroFields%rimecldgrau not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldgrau not associated and oneMicroFields%rimecldgrau associated")
    end if

    assOld=associated(micro_g(1)%rimecldhail)
    assNew=associated(oneMicroFields%rimecldhail)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecldhail, micro_g(1)%rimecldhail)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldhail associated but oneMicroFields%rimecldhail not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldhail not associated and oneMicroFields%rimecldhail associated")
    end if

    assOld=associated(micro_g(1)%rain2ice)
    assNew=associated(oneMicroFields%rain2ice)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2ice, micro_g(1)%rain2ice)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ice associated but oneMicroFields%rain2ice not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ice not associated and oneMicroFields%rain2ice associated")
    end if

    assOld=associated(micro_g(1)%rain2pr)
    assNew=associated(oneMicroFields%rain2pr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2pr, micro_g(1)%rain2pr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2pr associated but oneMicroFields%rain2pr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2pr not associated and oneMicroFields%rain2pr associated")
    end if

    assOld=associated(micro_g(1)%rain2sn)
    assNew=associated(oneMicroFields%rain2sn)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2sn, micro_g(1)%rain2sn)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2sn associated but oneMicroFields%rain2sn not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2sn not associated and oneMicroFields%rain2sn associated")
    end if

    assOld=associated(micro_g(1)%rain2ag)
    assNew=associated(oneMicroFields%rain2ag)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2ag, micro_g(1)%rain2ag)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ag associated but oneMicroFields%rain2ag not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ag not associated and oneMicroFields%rain2ag associated")
    end if

    assOld=associated(micro_g(1)%rain2gr)
    assNew=associated(oneMicroFields%rain2gr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2gr, micro_g(1)%rain2gr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2gr associated but oneMicroFields%rain2gr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2gr not associated and oneMicroFields%rain2gr associated")
    end if

    assOld=associated(micro_g(1)%rain2ha)
    assNew=associated(oneMicroFields%rain2ha)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2ha, micro_g(1)%rain2ha)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ha associated but oneMicroFields%rain2ha not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ha not associated and oneMicroFields%rain2ha associated")
    end if

    assOld=associated(micro_g(1)%rain2ha_xtra)
    assNew=associated(oneMicroFields%rain2ha_xtra)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2ha_xtra, micro_g(1)%rain2ha_xtra)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ha_xtra associated but oneMicroFields%rain2ha_xtra not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ha_xtra not associated and oneMicroFields%rain2ha_xtra associated")
    end if

    assOld=associated(micro_g(1)%ice2rain)
    assNew=associated(oneMicroFields%ice2rain)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%ice2rain, micro_g(1)%ice2rain)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%ice2rain associated but oneMicroFields%ice2rain not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%ice2rain not associated and oneMicroFields%ice2rain associated")
    end if

    assOld=associated(micro_g(1)%aggregate)
    assNew=associated(oneMicroFields%aggregate)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%aggregate, micro_g(1)%aggregate)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggregate associated but oneMicroFields%aggregate not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggregate not associated and oneMicroFields%aggregate associated")
    end if

    assOld=associated(micro_g(1)%aggrselfpris)
    assNew=associated(oneMicroFields%aggrselfpris)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%aggrselfpris, micro_g(1)%aggrselfpris)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrselfpris associated but oneMicroFields%aggrselfpris not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrselfpris not associated and oneMicroFields%aggrselfpris associated")
    end if

    assOld=associated(micro_g(1)%aggrselfsnow)
    assNew=associated(oneMicroFields%aggrselfsnow)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%aggrselfsnow, micro_g(1)%aggrselfsnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrselfsnow associated but oneMicroFields%aggrselfsnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrselfsnow not associated and oneMicroFields%aggrselfsnow associated")
    end if

    assOld=associated(micro_g(1)%aggrprissnow)
    assNew=associated(oneMicroFields%aggrprissnow)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%aggrprissnow, micro_g(1)%aggrprissnow)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrprissnow associated but oneMicroFields%aggrprissnow not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrprissnow not associated and oneMicroFields%aggrprissnow associated")
    end if

    assOld=associated(micro_g(1)%latheatvap)
    assNew=associated(oneMicroFields%latheatvap)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%latheatvap, micro_g(1)%latheatvap)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%latheatvap associated but oneMicroFields%latheatvap not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%latheatvap not associated and oneMicroFields%latheatvap associated")
    end if

    assOld=associated(micro_g(1)%latheatfrz)
    assNew=associated(oneMicroFields%latheatfrz)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%latheatfrz, micro_g(1)%latheatfrz)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%latheatfrz associated but oneMicroFields%latheatfrz not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%latheatfrz not associated and oneMicroFields%latheatfrz associated")
    end if

    assOld=associated(micro_g(1)%nuccldrt)
    assNew=associated(oneMicroFields%nuccldrt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%nuccldrt, micro_g(1)%nuccldrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nuccldrt associated but oneMicroFields%nuccldrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nuccldrt not associated and oneMicroFields%nuccldrt associated")
    end if

    assOld=associated(micro_g(1)%nuccldct)
    assNew=associated(oneMicroFields%nuccldct)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%nuccldct, micro_g(1)%nuccldct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nuccldct associated but oneMicroFields%nuccldct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nuccldct not associated and oneMicroFields%nuccldct associated")
    end if

    assOld=associated(micro_g(1)%nucicert)
    assNew=associated(oneMicroFields%nucicert)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%nucicert, micro_g(1)%nucicert)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nucicert associated but oneMicroFields%nucicert not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nucicert not associated and oneMicroFields%nucicert associated")
    end if

    assOld=associated(micro_g(1)%nucicect)
    assNew=associated(oneMicroFields%nucicect)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%nucicect, micro_g(1)%nucicect)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%nucicect associated but oneMicroFields%nucicect not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%nucicect not associated and oneMicroFields%nucicect associated")
    end if

    assOld=associated(micro_g(1)%inuchomrt)
    assNew=associated(oneMicroFields%inuchomrt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuchomrt, micro_g(1)%inuchomrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchomrt associated but oneMicroFields%inuchomrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchomrt not associated and oneMicroFields%inuchomrt associated")
    end if

    assOld=associated(micro_g(1)%inuchomct)
    assNew=associated(oneMicroFields%inuchomct)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuchomct, micro_g(1)%inuchomct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchomct associated but oneMicroFields%inuchomct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchomct not associated and oneMicroFields%inuchomct associated")
    end if

    assOld=associated(micro_g(1)%inuccontrt)
    assNew=associated(oneMicroFields%inuccontrt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuccontrt, micro_g(1)%inuccontrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuccontrt associated but oneMicroFields%inuccontrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuccontrt not associated and oneMicroFields%inuccontrt associated")
    end if

    assOld=associated(micro_g(1)%inuccontct)
    assNew=associated(oneMicroFields%inuccontct)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuccontct, micro_g(1)%inuccontct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuccontct associated but oneMicroFields%inuccontct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuccontct not associated and oneMicroFields%inuccontct associated")
    end if

    assOld=associated(micro_g(1)%inucifnrt)
    assNew=associated(oneMicroFields%inucifnrt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inucifnrt, micro_g(1)%inucifnrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inucifnrt associated but oneMicroFields%inucifnrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inucifnrt not associated and oneMicroFields%inucifnrt associated")
    end if

    assOld=associated(micro_g(1)%inucifnct)
    assNew=associated(oneMicroFields%inucifnct)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inucifnct, micro_g(1)%inucifnct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inucifnct associated but oneMicroFields%inucifnct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inucifnct not associated and oneMicroFields%inucifnct associated")
    end if

    assOld=associated(micro_g(1)%inuchazrt)
    assNew=associated(oneMicroFields%inuchazrt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuchazrt, micro_g(1)%inuchazrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchazrt associated but oneMicroFields%inuchazrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchazrt not associated and oneMicroFields%inuchazrt associated")
    end if

    assOld=associated(micro_g(1)%inuchazct)
    assNew=associated(oneMicroFields%inuchazct)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%inuchazct, micro_g(1)%inuchazct)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%inuchazct associated but oneMicroFields%inuchazct not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%inuchazct not associated and oneMicroFields%inuchazct associated")
    end if

    assOld=associated(micro_g(1)%vapliqt)
    assNew=associated(oneMicroFields%vapliqt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapliqt, micro_g(1)%vapliqt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapliqt associated but oneMicroFields%vapliqt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapliqt not associated and oneMicroFields%vapliqt associated")
    end if

    assOld=associated(micro_g(1)%vapicet)
    assNew=associated(oneMicroFields%vapicet)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapicet, micro_g(1)%vapicet)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapicet associated but oneMicroFields%vapicet not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapicet not associated and oneMicroFields%vapicet associated")
    end if

    assOld=associated(micro_g(1)%vapcldt)
    assNew=associated(oneMicroFields%vapcldt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapcldt, micro_g(1)%vapcldt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapcldt associated but oneMicroFields%vapcldt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapcldt not associated and oneMicroFields%vapcldt associated")
    end if

    assOld=associated(micro_g(1)%vapraint)
    assNew=associated(oneMicroFields%vapraint)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapraint, micro_g(1)%vapraint)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapraint associated but oneMicroFields%vapraint not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapraint not associated and oneMicroFields%vapraint associated")
    end if

    assOld=associated(micro_g(1)%vapprist)
    assNew=associated(oneMicroFields%vapprist)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapprist, micro_g(1)%vapprist)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapprist associated but oneMicroFields%vapprist not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapprist not associated and oneMicroFields%vapprist associated")
    end if

    assOld=associated(micro_g(1)%vapsnowt)
    assNew=associated(oneMicroFields%vapsnowt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapsnowt, micro_g(1)%vapsnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapsnowt associated but oneMicroFields%vapsnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapsnowt not associated and oneMicroFields%vapsnowt associated")
    end if

    assOld=associated(micro_g(1)%vapaggrt)
    assNew=associated(oneMicroFields%vapaggrt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapaggrt, micro_g(1)%vapaggrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapaggrt associated but oneMicroFields%vapaggrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapaggrt not associated and oneMicroFields%vapaggrt associated")
    end if

    assOld=associated(micro_g(1)%vapgraut)
    assNew=associated(oneMicroFields%vapgraut)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapgraut, micro_g(1)%vapgraut)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapgraut associated but oneMicroFields%vapgraut not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapgraut not associated and oneMicroFields%vapgraut associated")
    end if

    assOld=associated(micro_g(1)%vaphailt)
    assNew=associated(oneMicroFields%vaphailt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vaphailt, micro_g(1)%vaphailt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vaphailt associated but oneMicroFields%vaphailt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vaphailt not associated and oneMicroFields%vaphailt associated")
    end if

    assOld=associated(micro_g(1)%vapdrizt)
    assNew=associated(oneMicroFields%vapdrizt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%vapdrizt, micro_g(1)%vapdrizt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%vapdrizt associated but oneMicroFields%vapdrizt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%vapdrizt not associated and oneMicroFields%vapdrizt associated")
    end if

    assOld=associated(micro_g(1)%melticet)
    assNew=associated(oneMicroFields%melticet)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%melticet, micro_g(1)%melticet)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%melticet associated but oneMicroFields%melticet not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%melticet not associated and oneMicroFields%melticet associated")
    end if

    assOld=associated(micro_g(1)%meltprist)
    assNew=associated(oneMicroFields%meltprist)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%meltprist, micro_g(1)%meltprist)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltprist associated but oneMicroFields%meltprist not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltprist not associated and oneMicroFields%meltprist associated")
    end if

    assOld=associated(micro_g(1)%meltsnowt)
    assNew=associated(oneMicroFields%meltsnowt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%meltsnowt, micro_g(1)%meltsnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltsnowt associated but oneMicroFields%meltsnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltsnowt not associated and oneMicroFields%meltsnowt associated")
    end if

    assOld=associated(micro_g(1)%meltaggrt)
    assNew=associated(oneMicroFields%meltaggrt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%meltaggrt, micro_g(1)%meltaggrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltaggrt associated but oneMicroFields%meltaggrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltaggrt not associated and oneMicroFields%meltaggrt associated")
    end if

    assOld=associated(micro_g(1)%meltgraut)
    assNew=associated(oneMicroFields%meltgraut)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%meltgraut, micro_g(1)%meltgraut)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%meltgraut associated but oneMicroFields%meltgraut not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%meltgraut not associated and oneMicroFields%meltgraut associated")
    end if

    assOld=associated(micro_g(1)%melthailt)
    assNew=associated(oneMicroFields%melthailt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%melthailt, micro_g(1)%melthailt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%melthailt associated but oneMicroFields%melthailt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%melthailt not associated and oneMicroFields%melthailt associated")
    end if

    assOld=associated(micro_g(1)%cld2raint)
    assNew=associated(oneMicroFields%cld2raint)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%cld2raint, micro_g(1)%cld2raint)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%cld2raint associated but oneMicroFields%cld2raint not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%cld2raint not associated and oneMicroFields%cld2raint associated")
    end if

    assOld=associated(micro_g(1)%rimecldt)
    assNew=associated(oneMicroFields%rimecldt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecldt, micro_g(1)%rimecldt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldt associated but oneMicroFields%rimecldt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldt not associated and oneMicroFields%rimecldt associated")
    end if

    assOld=associated(micro_g(1)%rimecldsnowt)
    assNew=associated(oneMicroFields%rimecldsnowt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecldsnowt, micro_g(1)%rimecldsnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldsnowt associated but oneMicroFields%rimecldsnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldsnowt not associated and oneMicroFields%rimecldsnowt associated")
    end if

    assOld=associated(micro_g(1)%rimecldaggrt)
    assNew=associated(oneMicroFields%rimecldaggrt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecldaggrt, micro_g(1)%rimecldaggrt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldaggrt associated but oneMicroFields%rimecldaggrt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldaggrt not associated and oneMicroFields%rimecldaggrt associated")
    end if

    assOld=associated(micro_g(1)%rimecldgraut)
    assNew=associated(oneMicroFields%rimecldgraut)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecldgraut, micro_g(1)%rimecldgraut)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldgraut associated but oneMicroFields%rimecldgraut not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldgraut not associated and oneMicroFields%rimecldgraut associated")
    end if

    assOld=associated(micro_g(1)%rimecldhailt)
    assNew=associated(oneMicroFields%rimecldhailt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rimecldhailt, micro_g(1)%rimecldhailt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rimecldhailt associated but oneMicroFields%rimecldhailt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rimecldhailt not associated and oneMicroFields%rimecldhailt associated")
    end if

    assOld=associated(micro_g(1)%rain2icet)
    assNew=associated(oneMicroFields%rain2icet)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2icet, micro_g(1)%rain2icet)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2icet associated but oneMicroFields%rain2icet not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2icet not associated and oneMicroFields%rain2icet associated")
    end if

    assOld=associated(micro_g(1)%rain2prt)
    assNew=associated(oneMicroFields%rain2prt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2prt, micro_g(1)%rain2prt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2prt associated but oneMicroFields%rain2prt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2prt not associated and oneMicroFields%rain2prt associated")
    end if

    assOld=associated(micro_g(1)%rain2snt)
    assNew=associated(oneMicroFields%rain2snt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2snt, micro_g(1)%rain2snt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2snt associated but oneMicroFields%rain2snt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2snt not associated and oneMicroFields%rain2snt associated")
    end if

    assOld=associated(micro_g(1)%rain2agt)
    assNew=associated(oneMicroFields%rain2agt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2agt, micro_g(1)%rain2agt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2agt associated but oneMicroFields%rain2agt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2agt not associated and oneMicroFields%rain2agt associated")
    end if

    assOld=associated(micro_g(1)%rain2grt)
    assNew=associated(oneMicroFields%rain2grt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2grt, micro_g(1)%rain2grt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2grt associated but oneMicroFields%rain2grt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2grt not associated and oneMicroFields%rain2grt associated")
    end if

    assOld=associated(micro_g(1)%rain2hat)
    assNew=associated(oneMicroFields%rain2hat)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2hat, micro_g(1)%rain2hat)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2hat associated but oneMicroFields%rain2hat not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2hat not associated and oneMicroFields%rain2hat associated")
    end if

    assOld=associated(micro_g(1)%rain2ha_xtrat)
    assNew=associated(oneMicroFields%rain2ha_xtrat)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%rain2ha_xtrat, micro_g(1)%rain2ha_xtrat)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%rain2ha_xtrat associated but oneMicroFields%rain2ha_xtrat not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%rain2ha_xtrat not associated and oneMicroFields%rain2ha_xtrat associated")
    end if

    assOld=associated(micro_g(1)%ice2raint)
    assNew=associated(oneMicroFields%ice2raint)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%ice2raint, micro_g(1)%ice2raint)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%ice2raint associated but oneMicroFields%ice2raint not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%ice2raint not associated and oneMicroFields%ice2raint associated")
    end if

    assOld=associated(micro_g(1)%aggregatet)
    assNew=associated(oneMicroFields%aggregatet)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%aggregatet, micro_g(1)%aggregatet)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggregatet associated but oneMicroFields%aggregatet not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggregatet not associated and oneMicroFields%aggregatet associated")
    end if

    assOld=associated(micro_g(1)%aggrselfprist)
    assNew=associated(oneMicroFields%aggrselfprist)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%aggrselfprist, micro_g(1)%aggrselfprist)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrselfprist associated but oneMicroFields%aggrselfprist not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrselfprist not associated and oneMicroFields%aggrselfprist associated")
    end if

    assOld=associated(micro_g(1)%aggrselfsnowt)
    assNew=associated(oneMicroFields%aggrselfsnowt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%aggrselfsnowt, micro_g(1)%aggrselfsnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrselfsnowt associated but oneMicroFields%aggrselfsnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrselfsnowt not associated and oneMicroFields%aggrselfsnowt associated")
    end if

    assOld=associated(micro_g(1)%aggrprissnowt)
    assNew=associated(oneMicroFields%aggrprissnowt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%aggrprissnowt, micro_g(1)%aggrprissnowt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%aggrprissnowt associated but oneMicroFields%aggrprissnowt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%aggrprissnowt not associated and oneMicroFields%aggrprissnowt associated")
    end if

    assOld=associated(micro_g(1)%latheatvapt)
    assNew=associated(oneMicroFields%latheatvapt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%latheatvapt, micro_g(1)%latheatvapt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%latheatvapt associated but oneMicroFields%latheatvapt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%latheatvapt not associated and oneMicroFields%latheatvapt associated")
    end if

    assOld=associated(micro_g(1)%latheatfrzt)
    assNew=associated(oneMicroFields%latheatfrzt)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%latheatfrzt, micro_g(1)%latheatfrzt)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%latheatfrzt associated but oneMicroFields%latheatfrzt not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%latheatfrzt not associated and oneMicroFields%latheatfrzt associated")
    end if

    assOld=associated(micro_g(1)%accpr)
    assNew=associated(oneMicroFields%accpr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%accpr, micro_g(1)%accpr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpr associated but oneMicroFields%accpr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpr not associated and oneMicroFields%accpr associated")
    end if

    assOld=associated(micro_g(1)%accpp)
    assNew=associated(oneMicroFields%accpp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%accpp, micro_g(1)%accpp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpp associated but oneMicroFields%accpp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpp not associated and oneMicroFields%accpp associated")
    end if

    assOld=associated(micro_g(1)%accps)
    assNew=associated(oneMicroFields%accps)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%accps, micro_g(1)%accps)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accps associated but oneMicroFields%accps not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accps not associated and oneMicroFields%accps associated")
    end if

    assOld=associated(micro_g(1)%accpa)
    assNew=associated(oneMicroFields%accpa)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%accpa, micro_g(1)%accpa)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpa associated but oneMicroFields%accpa not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpa not associated and oneMicroFields%accpa associated")
    end if

    assOld=associated(micro_g(1)%accpg)
    assNew=associated(oneMicroFields%accpg)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%accpg, micro_g(1)%accpg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpg associated but oneMicroFields%accpg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpg not associated and oneMicroFields%accpg associated")
    end if

    assOld=associated(micro_g(1)%accph)
    assNew=associated(oneMicroFields%accph)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%accph, micro_g(1)%accph)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accph associated but oneMicroFields%accph not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accph not associated and oneMicroFields%accph associated")
    end if

    assOld=associated(micro_g(1)%accpd)
    assNew=associated(oneMicroFields%accpd)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%accpd, micro_g(1)%accpd)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%accpd associated but oneMicroFields%accpd not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%accpd not associated and oneMicroFields%accpd associated")
    end if

    assOld=associated(micro_g(1)%pcprr)
    assNew=associated(oneMicroFields%pcprr)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcprr, micro_g(1)%pcprr)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprr associated but oneMicroFields%pcprr not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprr not associated and oneMicroFields%pcprr associated")
    end if

    assOld=associated(micro_g(1)%pcprp)
    assNew=associated(oneMicroFields%pcprp)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcprp, micro_g(1)%pcprp)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprp associated but oneMicroFields%pcprp not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprp not associated and oneMicroFields%pcprp associated")
    end if

    assOld=associated(micro_g(1)%pcprs)
    assNew=associated(oneMicroFields%pcprs)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcprs, micro_g(1)%pcprs)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprs associated but oneMicroFields%pcprs not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprs not associated and oneMicroFields%pcprs associated")
    end if

    assOld=associated(micro_g(1)%pcpra)
    assNew=associated(oneMicroFields%pcpra)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcpra, micro_g(1)%pcpra)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpra associated but oneMicroFields%pcpra not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpra not associated and oneMicroFields%pcpra associated")
    end if

    assOld=associated(micro_g(1)%pcprg)
    assNew=associated(oneMicroFields%pcprg)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcprg, micro_g(1)%pcprg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprg associated but oneMicroFields%pcprg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprg not associated and oneMicroFields%pcprg associated")
    end if

    assOld=associated(micro_g(1)%pcprh)
    assNew=associated(oneMicroFields%pcprh)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcprh, micro_g(1)%pcprh)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprh associated but oneMicroFields%pcprh not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprh not associated and oneMicroFields%pcprh associated")
    end if

    assOld=associated(micro_g(1)%pcprd)
    assNew=associated(oneMicroFields%pcprd)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcprd, micro_g(1)%pcprd)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcprd associated but oneMicroFields%pcprd not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcprd not associated and oneMicroFields%pcprd associated")
    end if

    assOld=associated(micro_g(1)%pcpg)
    assNew=associated(oneMicroFields%pcpg)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%pcpg, micro_g(1)%pcpg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%pcpg associated but oneMicroFields%pcpg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%pcpg not associated and oneMicroFields%pcpg associated")
    end if

    assOld=associated(micro_g(1)%qpcpg)
    assNew=associated(oneMicroFields%qpcpg)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%qpcpg, micro_g(1)%qpcpg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%qpcpg associated but oneMicroFields%qpcpg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%qpcpg not associated and oneMicroFields%qpcpg associated")
    end if

    assOld=associated(micro_g(1)%dpcpg)
    assNew=associated(oneMicroFields%dpcpg)
    if (assOld .and. assNew) then
       call ToCopy(oneMicroFields%dpcpg, micro_g(1)%dpcpg)
    else if (assOld .and. (.not. assNew)) then
       call fatal_error(h//" micro_g(1)%dpcpg associated but oneMicroFields%dpcpg not associated")
    else if ((.not. assOld) .and. assNew) then
       call fatal_error(h//" micro_g(1)%dpcpg not associated and oneMicroFields%dpcpg associated")
    end if
  end subroutine DeepCopyFromMicroFields

end module mem_micro
