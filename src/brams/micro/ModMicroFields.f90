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


module ModMicroFields

  use iso_fortran_env, only: &
       int64
  
  use ModNamelistFile, only: &
       NamelistFile
  
  use ModNodeDimensions, only: &
       NodeDimensions

  use ModMicControl, only: &
       MicControl

  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModVarTable, only: &
       VarTable, &
       InsertVarTable
  
  implicit none

  private

  public :: MicroFields
  public :: CreateMicroFields
  public :: DestroyMicroFields
  public :: DumpMicroFields
  public :: InsertMicroFieldsAtVarTable

  type MicroFields

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

  end type MicroFields


contains                  



  function CreateMicroFields(gridId, oneNamelistFile, oneNodeDimensions, oneMicControl) result(res)
    integer, intent(in) :: gridId
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(MicControl), pointer, intent(in) :: oneMicControl
    type(MicroFields), pointer :: res

    integer :: ierr
    integer :: mzp
    integer :: mxp
    integer :: myp
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateMicroFields)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneNodeDimensions)) then
       call fatal_error(h//" oneNodeDimensions not associated")
    else if (.not. associated(oneMicControl)) then
       call fatal_error(h//" oneMicControl not associated")
    end if
    
    mzp=oneNodeDimensions%mzp
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if

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
       allocate (res%rcp(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate rcp fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%rcp  =0.0

       !- rain
       allocate (res%rrp  (mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate rrp fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%rrp  =0.0
       !- for this scheme, the rain rate below will
       !- account for rain+ice+snow+graupel
       allocate (res%accpr(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate accpr fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%accpr=0.0
       allocate (res%pcprr(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate pcprr fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%pcprr=0.0

       !- ice
       allocate (res%rpp  (mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate rpp fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%rpp  =0.0
       !- don t need to be allocated, see coments above

       !- snow
       allocate (res%rsp  (mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate rsp fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%rsp=0.0
       !- the rates bellow will account for snow and ice
       allocate (res%accps(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate accps fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%accps =0.0
       allocate (res%pcprs(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate pcprs fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%pcprs =0.0

       !- graupel
       allocate (res%rgp  (mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate rgp fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%rgp  =0.0
       !- the rates bellow will account for only graupel
       allocate (res%accpg(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate accpg fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%accpg=0.0
       allocate (res%pcprg(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate pcprg fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%pcprg=0.0

       if(oneMicControl%mcphys_type  == 2 .or. oneMicControl%mcphys_type  == 3) then ! only for double-moment and 
          !- number concentration for cloud/rain/ice
          !- obs : ccp don t need to be allocated for the single-moment
          !- cloud water scheme (the same for CCN and IFN).
          allocate (res%crp  (mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate crp fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%crp  =0.0 
          allocate (res%cpp  (mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate cpp fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%cpp  =0.0 
          !---these should not be allocated for oneMicControl%mcphys_type  == 2 because
          !---they are not used for this option
          !            !ST
          !          allocate(res%ccp  (mzp,mxp,myp)) ;res%ccp  =0.0 
          !            allocate(res%cccnp(mzp,mxp,myp)) ;res%cccnp=0.0 !;endif 
          !            allocate(res%cifnp(mzp,mxp,myp)) ;res%cifnp=0.0 !;endif 
          !            !ST
       endif
       !- only for cloud water double-moment and aerosol aware microphysics         
       if(oneMicControl%mcphys_type  == 3) then ! only for double-moment and 
          allocate (res%ccp  (mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate ccp fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%ccp  =0.0 
          allocate (res%cccnp(mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate cccnp fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%cccnp=0.0 !;endif 
          allocate (res%cifnp(mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate cifnp fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%cifnp=0.0 !;endif 
       endif

       !- 3D cloud fraction from GFDL cloud microphysics and GF convection
       if(oneMicControl%mcphys_type  == 4 .or. oneNamelistFile%nnqparm(gridId) == 8) then 
          allocate (res%cldfr  (mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate cldfr fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%cldfr  =0.0 
       endif

       !- for consistency with the other parts of BRAMS
       !- pgcp will be the total precipitation rate
       allocate (res%pcpg (mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate pcpg fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%pcpg =0.0
       !-the allocations below are tmp for leaf-3
       allocate (res%qpcpg(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate qpcpg fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%qpcpg=0.0         
       allocate (res%dpcpg(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate dpcpg fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%dpcpg=0.0

       !- allocation of memory for effective radius for RRTMG
       if(oneNamelistFile%ilwrtyp==6 .or. oneNamelistFile%iswrtyp==6 ) then
          allocate (res%rei  (mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate rei fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%rei  =0.0  
          allocate (res%rel  (mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate rel fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%rel  =0.0
       endif

    else  ! for the traditional RAMS microphysics

       if (oneMicControl%level >= 2 ) then
          allocate (res%rcp(mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate rcp fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%rcp    =0.0
       endif
       if (oneMicControl%level >= 3) then         
          if(oneMicControl%irain >= 1)  then
             allocate (res%rrp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate rrp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%rrp  =0.0
             allocate (res%accpr(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate accpr fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%accpr=0.0
             allocate (res%pcprr(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcprr fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcprr=0.0
             allocate (res%pcpvr(mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcpvr fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcpvr=0.0
             allocate (res%q2   (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate q2 fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%q2   =0.0
          endif
          if(oneMicControl%ipris >= 1)  then
             allocate (res%rpp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate rpp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%rpp  =0.0
             allocate (res%accpp(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate accpp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%accpp=0.0
             allocate (res%pcprp(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcprp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcprp=0.0
             allocate (res%pcpvp(mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcpvp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcpvp=0.0
          endif
          if(oneMicControl%isnow >= 1)  then
             allocate (res%rsp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate rsp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%rsp   =0.0
             allocate (res%accps(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate accps fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%accps =0.0
             allocate (res%pcprs(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcprs fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcprs =0.0
             allocate (res%pcpvs(mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcpvs fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcpvs =0.0
          endif
          if(oneMicControl%iaggr >= 1)  then
             allocate (res%rap  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate rap fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%rap  =0.0
             allocate (res%accpa(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate accpa fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%accpa=0.0
             allocate (res%pcpra(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcpra fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcpra=0.0
             allocate (res%pcpva(mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcpva fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcpva=0.0
          endif
          if(oneMicControl%igraup >= 1) then
             allocate (res%rgp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate rgp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%rgp  =0.0
             allocate (res%accpg(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate accpg fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%accpg=0.0
             allocate (res%pcprg(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcprg fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcprg=0.0
             allocate (res%pcpvg(mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcpvg fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcpvg=0.0
             allocate (res%q6   (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate q6 fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%q6   =0.0
          endif
          if(oneMicControl%ihail >= 1)  then
             allocate (res%rhp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate rhp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%rhp  =0.0
             allocate (res%accph(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate accph fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%accph=0.0
             allocate (res%pcprh(mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcprh fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcprh=0.0
             allocate (res%pcpvh(mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate pcpvh fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%pcpvh=0.0
             allocate (res%q7   (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate q7 fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%q7   =0.0
          endif
          if(oneMicControl%jnmb(1) >= 5)  then
             allocate (res%ccp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate ccp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%ccp  =0.0
          endif
          if(oneMicControl%jnmb(2) == 5)  then
             allocate (res%crp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate crp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%crp  =0.0
          endif
          if(oneMicControl%jnmb(3) >= 5)  then
             allocate (res%cpp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate cpp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%cpp  =0.0
          endif
          if(oneMicControl%jnmb(4) == 5)  then
             allocate (res%csp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate csp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%csp  =0.0
          endif
          if(oneMicControl%jnmb(5) == 5)  then
             allocate (res%cap  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate cap fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%cap  =0.0
          endif
          if(oneMicControl%jnmb(6) == 5)  then
             allocate (res%cgp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate cgp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%cgp  =0.0
          endif
          if(oneMicControl%jnmb(7) == 5)  then
             allocate (res%chp  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate chp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%chp  =0.0
          endif
          if(oneMicControl%icloud  >= 5)  then
             allocate (res%cccnp(mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate cccnp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%cccnp=0.0
          endif
          if(oneMicControl%ipris   >= 5)  then
             allocate (res%cifnp(mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate cifnp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%cifnp=0.0
          endif

          if(oneMicControl%icloud >= 5)   then
             allocate (res%cccmp(mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate cccmp fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%cccmp=0.0
          endif

          allocate (res%pcpg (mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate pcpg fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%pcpg =0.0
          allocate (res%qpcpg(mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate qpcpg fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%qpcpg=0.0
          allocate (res%dpcpg(mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate dpcpg fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%dpcpg=0.0

          !- only for 2M microphysics
          if(oneMicControl%mcphys_type == 1)  then
             if(oneMicControl%idriz >= 1 )  then
                allocate (res%rdp  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rdp fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rdp  =0.0
                allocate (res%accpd(mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate accpd fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%accpd=0.0
                allocate (res%pcprd(mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate pcprd fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%pcprd=0.0
                allocate (res%pcpvd(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate pcpvd fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%pcpvd=0.0
             endif

             if(oneMicControl%jnmb(8) >= 5)  then
                allocate (res%cdp  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate cdp fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%cdp  =0.0
             endif
             if(oneMicControl%idriz>= 5)  then
                allocate (res%gccnp(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate gccnp fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%gccnp=0.0
             endif
             if(oneMicControl%idriz   >= 5)  then
                allocate (res%gccmp(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate gccmp fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%gccmp=0.0
             endif

             if(oneMicControl%iccnlev >= 2 .and. oneMicControl%jnmb(1) >= 5)  then
                allocate (res%cnm1p(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate cnm1p fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%cnm1p=0.0
             endif
             if(oneMicControl%iccnlev >= 2 .and. oneMicControl%jnmb(2) >= 1)  then
                allocate (res%cnm2p(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate cnm2p fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%cnm2p=0.0
             endif
             if(oneMicControl%iccnlev >= 2 .and. oneMicControl%jnmb(3) >= 1)  then
                allocate (res%cnm3p(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate cnm3p fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%cnm3p=0.0
             endif
             if(oneMicControl%iccnlev >= 2 .and. oneMicControl%jnmb(8) >= 1)  then
                allocate (res%cnm8p(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate cnm8p fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%cnm8p=0.0
             endif

             if(oneMicControl%idust == 1 .or. oneMicControl%imd1flg == 1)  then
                allocate (res%md1np(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate md1np fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%md1np=0.0
             endif
             if(oneMicControl%idust == 1 .or. oneMicControl%imd2flg == 1)  then
                allocate (res%md2np(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate md2np fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%md2np=0.0
             endif
             if(oneMicControl%isalt == 1) then
                allocate (res%salt_filmp(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate salt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%salt_filmp =0.0
             endif
             if(oneMicControl%isalt == 1) then
                allocate (res%salt_jetp (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate salt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%salt_jetp  =0.0
             endif
             if(oneMicControl%isalt == 1) then
                allocate (res%salt_spmp (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate salt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%salt_spmp  =0.0
             endif


             !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
             if(oneMicControl%imbudget>=1 .or. oneMicControl%imbudtot>=1) then
                allocate (res%latheatvap(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate latheatvap fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%latheatvap=0.0
                allocate (res%latheatfrz(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate latheatfrz fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%latheatfrz=0.0
             endif
             if(oneMicControl%imbudget>=1) then
                allocate (res%nuccldr  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate nuccldr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%nuccldr  =0.0
                allocate (res%nuccldc  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate nuccldc fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%nuccldc  =0.0
                allocate (res%cld2rain (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate cld2rain fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%cld2rain =0.0
                allocate (res%ice2rain (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate ice2rain fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%ice2rain =0.0
                allocate (res%nucicer  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate nucicer fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%nucicer  =0.0
                allocate (res%nucicec  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate nucicec fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%nucicec  =0.0
                allocate (res%vapliq(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapliq fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapliq   =0.0   
                allocate (res%vapice(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapice fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapice   =0.0
                allocate (res%meltice  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate meltice fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%meltice  =0.0
                allocate (res%rimecld  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecld fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecld  =0.0
                allocate (res%rain2ice (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2ice fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2ice =0.0
                allocate (res%aggregate(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate aggregate fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%aggregate=0.0
             endif
             if(oneMicControl%imbudget==2) then
                allocate (res%inuchomr    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuchomr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuchomr    =0.0 
                allocate (res%inuchomc    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuchomc fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuchomc    =0.0
                allocate (res%inuccontr   (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuccontr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuccontr   =0.0
                allocate (res%inuccontc   (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuccontc fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuccontc   =0.0
                allocate (res%inucifnr    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inucifnr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inucifnr    =0.0
                allocate (res%inucifnc    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inucifnc fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inucifnc    =0.0
                allocate (res%inuchazr    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuchazr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuchazr    =0.0
                allocate (res%inuchazc    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuchazc fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuchazc    =0.0
                allocate (res%vapcld   (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapcld fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapcld  =0.0
                allocate (res%vaprain     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vaprain fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vaprain  =0.0
                allocate (res%vappris     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vappris fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vappris  =0.0
                allocate (res%vapsnow     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapsnow fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapsnow  =0.0   
                allocate (res%vapaggr     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapaggr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapaggr  =0.0
                allocate (res%vapgrau     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapgrau fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapgrau  =0.0
                allocate (res%vaphail     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vaphail fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vaphail  =0.0
                allocate (res%vapdriz     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapdriz fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapdriz  =0.0
                allocate (res%meltpris    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate meltpris fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%meltpris    =0.0
                allocate (res%meltsnow    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate meltsnow fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%meltsnow    =0.0
                allocate (res%meltaggr    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate meltaggr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%meltaggr    =0.0
                allocate (res%meltgrau    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate meltgrau fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%meltgrau    =0.0
                allocate (res%melthail    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate melthail fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%melthail    =0.0
                allocate (res%rimecldsnow (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecldsnow fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecldsnow =0.0
                allocate (res%rimecldaggr (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecldaggr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecldaggr =0.0
                allocate (res%rimecldgrau (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecldgrau fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecldgrau =0.0
                allocate (res%rimecldhail (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecldhail fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecldhail =0.0
                allocate (res%rain2pr     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2pr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2pr  =0.0
                allocate (res%rain2sn     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2sn fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2sn  =0.0
                allocate (res%rain2ag     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2ag fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2ag  =0.0
                allocate (res%rain2gr     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2gr fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2gr  =0.0
                allocate (res%rain2ha     (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2ha fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2ha  =0.0
                allocate (res%rain2ha_xtra(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2ha fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2ha_xtra=0.0
                allocate (res%aggrselfpris(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate aggrselfpris fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%aggrselfpris=0.0
                allocate (res%aggrselfsnow(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate aggrselfsnow fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%aggrselfsnow=0.0
                allocate (res%aggrprissnow(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate aggrprissnow fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%aggrprissnow=0.0
             endif
             !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
             if(oneMicControl%imbudtot>=1) then
                allocate (res%nuccldrt   (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate nuccldrt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%nuccldrt=0.0
                allocate (res%nuccldct   (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate nuccldct fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%nuccldct=0.0
                allocate (res%cld2raint  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate cld2raint fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%cld2raint  =0.0
                allocate (res%ice2raint  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate ice2raint fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%ice2raint  =0.0
                allocate (res%nucicert   (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate nucicert fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%nucicert=0.0
                allocate (res%nucicect   (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate nucicect fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%nucicect=0.0
                allocate (res%vapliqt    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapliqt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapliqt=0.0
                allocate (res%vapicet    (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapicet fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapicet=0.0
                allocate (res%melticet   (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate melticet fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%melticet=0.0
                allocate (res%rimecldt   (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecldt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecldt=0.0
                allocate (res%rain2icet  (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2icet fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2icet  =0.0
                allocate (res%aggregatet (mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate aggregatet fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%aggregatet =0.0
                allocate (res%latheatvapt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate latheatvapt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%latheatvapt=0.0
                allocate (res%latheatfrzt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate latheatfrzt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%latheatfrzt=0.0
             endif
             if(oneMicControl%imbudtot==2) then
                allocate (res%inuchomrt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuchomrt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuchomrt    =0.0
                allocate (res%inuchomct(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuchomct fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuchomct    =0.0
                allocate (res%inuccontrt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuccontrt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuccontrt   =0.0
                allocate (res%inuccontct(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuccontct fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuccontct   =0.0
                allocate (res%inucifnrt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inucifnrt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inucifnrt    =0.0
                allocate (res%inucifnct(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inucifnct fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inucifnct    =0.0
                allocate (res%inuchazrt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuchazrt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuchazrt    =0.0
                allocate (res%inuchazct(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate inuchazct fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%inuchazct    =0.0
                allocate (res%vapcldt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapcldt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapcldt      =0.0
                allocate (res%vapraint(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapraint fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapraint     =0.0
                allocate (res%vapprist(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapprist fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapprist     =0.0
                allocate (res%vapsnowt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapsnowt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapsnowt     =0.0
                allocate (res%vapaggrt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapaggrt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapaggrt     =0.0
                allocate (res%vapgraut(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapgraut fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapgraut     =0.0
                allocate (res%vaphailt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vaphailt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vaphailt     =0.0
                allocate (res%vapdrizt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate vapdrizt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%vapdrizt     =0.0
                allocate (res%meltprist(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate meltprist fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%meltprist    =0.0
                allocate (res%meltsnowt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate meltsnowt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%meltsnowt    =0.0
                allocate (res%meltaggrt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate meltaggrt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%meltaggrt    =0.0
                allocate (res%meltgraut(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate meltgraut fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%meltgraut    =0.0
                allocate (res%melthailt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate melthailt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%melthailt    =0.0
                allocate (res%rimecldsnowt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecldsnowt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecldsnowt =0.0
                allocate (res%rimecldaggrt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecldaggrt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecldaggrt =0.0
                allocate (res%rimecldgraut(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecldgraut fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecldgraut =0.0
                allocate (res%rimecldhailt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rimecldhailt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rimecldhailt =0.0
                allocate (res%rain2prt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2prt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2prt     =0.0
                allocate (res%rain2snt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2snt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2snt     =0.0
                allocate (res%rain2agt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2agt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2agt     =0.0
                allocate (res%rain2grt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2grt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2grt     =0.0
                allocate (res%rain2hat(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2hat fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2hat     =0.0
                allocate (res%rain2ha_xtrat(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate rain2ha fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%rain2ha_xtrat=0.0
                allocate (res%aggrselfprist(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate aggrselfprist fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%aggrselfprist=0.0
                allocate (res%aggrselfsnowt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate aggrselfsnowt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%aggrselfsnowt=0.0
                allocate (res%aggrprissnowt(mzp,mxp,myp), stat=ierr)
                if (ierr /= 0) then
                   write(str(1),"(i8)") ierr
                   call fatal_error(h//" allocate aggrprissnowt fails with stat="//&
                        trim(adjustl(str(1))))
                end if
                res%aggrprissnowt=0.0
             endif
          endif! oneMicControl%mcphys_type=1     
          !- allocation of memory for effective radius for RRTMG
          if(oneNamelistFile%ilwrtyp==6 .or. oneNamelistFile%iswrtyp==6 ) then
             allocate (res%rei  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate rei fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%rei  =0.0  
             allocate (res%rel  (mzp,mxp,myp), stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                call fatal_error(h//" allocate rel fails with stat="//&
                     trim(adjustl(str(1))))
             end if
             res%rel  =0.0
          endif
       endif ! oneMicControl%level >=3 
    endif
  end function CreateMicroFields



!!$  subroutine filltab_micro(micro,microm,imean,n1,n2,n3,ng)
!!$
!!$    implicit none
!!$    include "constants.h"
!!$    type (MicroFields) :: micro,microm
!!$    integer, intent(in) :: imean,n1,n2,n3,ng
!!$    integer(kind=i8) :: npts
!!$    real, pointer :: var,varm
!!$
!!$    ! Fill pointers to arrays into variable tables
!!$
!!$    npts=n1*n2*n3
!!$    if (associated(micro%rcp))   &
!!$         call InsertVTab (micro%rcp,microm%rcp  &
!!$         ,ng, npts, imean,  &
!!$         'RCP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%rdp))   &
!!$         call InsertVTab (micro%rdp,microm%rdp  &
!!$         ,ng, npts, imean,  &
!!$         'RDP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%rrp))   &
!!$         call InsertVTab (micro%rrp,microm%rrp  &
!!$         ,ng, npts, imean,  &
!!$         'RRP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%rpp))   &
!!$         call InsertVTab (micro%rpp,microm%rpp  &
!!$         ,ng, npts, imean,  &
!!$         'RPP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%rsp))   &
!!$         call InsertVTab (micro%rsp,microm%rsp  &
!!$         ,ng, npts, imean,  &
!!$         'RSP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%rap))   &
!!$         call InsertVTab (micro%rap,microm%rap  &
!!$         ,ng, npts, imean,  &
!!$         'RAP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%rgp))   &
!!$         call InsertVTab (micro%rgp,microm%rgp  &
!!$         ,ng, npts, imean,  &
!!$         'RGP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%rhp))   &
!!$         call InsertVTab (micro%rhp,microm%rhp  &
!!$         ,ng, npts, imean,  &
!!$         'RHP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%ccp))   &
!!$         call InsertVTab (micro%ccp,microm%ccp  &
!!$         ,ng, npts, imean,  &
!!$         'CCP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cdp))   &
!!$         call InsertVTab (micro%cdp,microm%cdp  &
!!$         ,ng, npts, imean,  &
!!$         'CDP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%crp))   &
!!$         call InsertVTab (micro%crp,microm%crp  &
!!$         ,ng, npts, imean,  &
!!$         'CRP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cpp))   &
!!$         call InsertVTab (micro%cpp,microm%cpp  &
!!$         ,ng, npts, imean,  &
!!$         'CPP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%csp))   &
!!$         call InsertVTab (micro%csp,microm%csp  &
!!$         ,ng, npts, imean,  &
!!$         'CSP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cap))   &
!!$         call InsertVTab (micro%cap,microm%cap  &
!!$         ,ng, npts, imean,  &
!!$         'CAP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cgp))   &
!!$         call InsertVTab (micro%cgp,microm%cgp  &
!!$         ,ng, npts, imean,  &
!!$         'CGP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%chp))   &
!!$         call InsertVTab (micro%chp,microm%chp  &
!!$         ,ng, npts, imean,  &
!!$         'CHP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cccnp)) &
!!$         call InsertVTab (micro%cccnp,microm%cccnp  &
!!$         ,ng, npts, imean,  &
!!$         'CCCNP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%gccnp)) &
!!$         call InsertVTab (micro%gccnp,microm%gccnp  &
!!$         ,ng, npts, imean,  &
!!$         'GCCNP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cifnp)) &
!!$         call InsertVTab (micro%cifnp,microm%cifnp  &
!!$         ,ng, npts, imean,  &
!!$         'CIFNP :3:hist:anal:mpti:mpt3:mpt1')
!!$
!!$    if (associated(micro%q2))   &
!!$         call InsertVTab (micro%q2,microm%q2  &
!!$         ,ng, npts, imean,  &
!!$         'Q2 :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%q6)) &
!!$         call InsertVTab (micro%q6,microm%q6  &
!!$         ,ng, npts, imean,  &
!!$         'Q6 :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%q7)) &
!!$         call InsertVTab (micro%q7,microm%q7  &
!!$         ,ng, npts, imean,  &
!!$         'Q7 :3:hist:anal:mpti:mpt3')
!!$
!!$    if (associated(micro%cccmp)) &
!!$         call InsertVTab (micro%cccmp,microm%cccmp  &
!!$         ,ng, npts, imean,  &
!!$         'CCCMP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%gccmp)) &
!!$         call InsertVTab (micro%gccmp,microm%gccmp  &
!!$         ,ng, npts, imean,  &
!!$         'GCCMP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cnm1p)) &
!!$         call InsertVTab (micro%cnm1p,microm%cnm1p  &
!!$         ,ng, npts, imean,  &
!!$         'CNM1P :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cnm2p)) &
!!$         call InsertVTab (micro%cnm2p,microm%cnm2p  &
!!$         ,ng, npts, imean,  &
!!$         'CNM2P :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cnm3p)) &
!!$         call InsertVTab (micro%cnm3p,microm%cnm3p  &
!!$         ,ng, npts, imean,  &
!!$         'CNM3P :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%cnm8p)) &
!!$         call InsertVTab (micro%cnm8p,microm%cnm8p  &
!!$         ,ng, npts, imean,  &
!!$         'CNM8P :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%md1np)) &
!!$         call InsertVTab (micro%md1np,microm%md1np  &
!!$         ,ng, npts, imean,  &
!!$         'MD1NP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%md2np)) &
!!$         call InsertVTab (micro%md2np,microm%md2np  &
!!$         ,ng, npts, imean,  &
!!$         'MD2NP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%salt_filmp)) &
!!$         call InsertVTab (micro%salt_filmp,microm%salt_filmp  &
!!$         ,ng, npts, imean,  &
!!$         'SALT_FILMP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%salt_jetp)) &
!!$         call InsertVTab (micro%salt_jetp,microm%salt_jetp  &
!!$         ,ng, npts, imean,  &
!!$         'SALT_JETP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%salt_spmp)) &
!!$         call InsertVTab (micro%salt_spmp,microm%salt_spmp  &
!!$         ,ng, npts, imean,  &
!!$         'SALT_SPMP :3:hist:anal:mpti:mpt3:mpt1')
!!$    if (associated(micro%rei))   &
!!$         call InsertVTab (micro%rei,microm%rei  &
!!$         ,ng, npts, imean,  &
!!$         'REI :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rel))   &
!!$         call InsertVTab (micro%rel,microm%rel  &
!!$         ,ng, npts, imean,  &
!!$         'REL :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%cldfr))   &
!!$         call InsertVTab (micro%cldfr,microm%cldfr  &
!!$         ,ng, npts, imean,  &
!!$         'CLDFR :3:hist:anal:mpti:mpt3')
!!$
!!$
!!$    !VERTICAL PRECIPITATION RATES
!!$    if (associated(micro%pcpvr)) &
!!$         call InsertVTab (micro%pcpvr,microm%pcpvr  &
!!$         ,ng, npts, imean,  &
!!$         'PCPVR :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%pcpvp)) &
!!$         call InsertVTab (micro%pcpvp,microm%pcpvp  &
!!$         ,ng, npts, imean,  &
!!$         'PCPVP :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%pcpvs)) &
!!$         call InsertVTab (micro%pcpvs,microm%pcpvs  &
!!$         ,ng, npts, imean,  &
!!$         'PCPVS :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%pcpva)) &
!!$         call InsertVTab (micro%pcpva,microm%pcpva  &
!!$         ,ng, npts, imean,  &
!!$         'PCPVA :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%pcpvg)) &
!!$         call InsertVTab (micro%pcpvg,microm%pcpvg  &
!!$         ,ng, npts, imean,  &
!!$         'PCPVG :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%pcpvh)) &
!!$         call InsertVTab (micro%pcpvh,microm%pcpvh  &
!!$         ,ng, npts, imean,  &
!!$         'PCPVH :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%pcpvd)) &
!!$         call InsertVTab (micro%pcpvd,microm%pcpvd  &
!!$         ,ng, npts, imean,  &
!!$         'PCPVD :3:hist:anal:mpti:mpt3')
!!$
!!$    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (instantaneous)
!!$    if (associated(micro%nuccldr)) &
!!$         call InsertVTab (micro%nuccldr,microm%nuccldr  &
!!$         ,ng, npts, imean,  &
!!$         'NUCCLDR :3:hist:anal:mpt3')
!!$    if (associated(micro%nuccldc)) &
!!$         call InsertVTab (micro%nuccldc,microm%nuccldc  &
!!$         ,ng, npts, imean,  &
!!$         'NUCCLDC :3:hist:anal:mpt3')
!!$
!!$    if (associated(micro%nucicer)) &
!!$         call InsertVTab (micro%nucicer,microm%nucicer  &
!!$         ,ng, npts, imean,  &
!!$         'NUCICER :3:hist:anal:mpt3')
!!$    if (associated(micro%nucicec)) &
!!$         call InsertVTab (micro%nucicec,microm%nucicec  &
!!$         ,ng, npts, imean,  &
!!$         'NUCICEC :3:hist:anal:mpt3')
!!$    if (associated(micro%inuchomr)) &
!!$         call InsertVTab (micro%inuchomr,microm%inuchomr  &
!!$         ,ng, npts, imean,  &
!!$         'INUCHOMR :3:hist:anal:mpt3')
!!$    if (associated(micro%inuchomc)) &
!!$         call InsertVTab (micro%inuchomc,microm%inuchomc  &
!!$         ,ng, npts, imean,  &
!!$         'INUCHOMC :3:hist:anal:mpt3')
!!$    if (associated(micro%inuccontr)) &
!!$         call InsertVTab (micro%inuccontr,microm%inuccontr  &
!!$         ,ng, npts, imean,  &
!!$         'INUCCONTR :3:hist:anal:mpt3')
!!$    if (associated(micro%inuccontc)) &
!!$         call InsertVTab (micro%inuccontc,microm%inuccontc  &
!!$         ,ng, npts, imean,  &
!!$         'INUCCONTC :3:hist:anal:mpt3')
!!$    if (associated(micro%inucifnr)) &
!!$         call InsertVTab (micro%inucifnr,microm%inucifnr  &
!!$         ,ng, npts, imean,  &
!!$         'INUCIFNR :3:hist:anal:mpt3')
!!$    if (associated(micro%inucifnc)) &
!!$         call InsertVTab (micro%inucifnc,microm%inucifnc  &
!!$         ,ng, npts, imean,  &
!!$         'INUCIFNC :3:hist:anal:mpt3')
!!$    if (associated(micro%inuchazr)) &
!!$         call InsertVTab (micro%inuchazr,microm%inuchazr  &
!!$         ,ng, npts, imean,  &
!!$         'INUCHAZR :3:hist:anal:mpt3')
!!$    if (associated(micro%inuchazc)) &
!!$         call InsertVTab (micro%inuchazc,microm%inuchazc  &
!!$         ,ng, npts, imean,  &
!!$         'INUCHAZC :3:hist:anal:mpt3')
!!$
!!$    if (associated(micro%vapliq)) &
!!$         call InsertVTab (micro%vapliq,microm%vapliq  &
!!$         ,ng, npts, imean,  &
!!$         'VAPLIQ :3:hist:anal:mpt3')
!!$    if (associated(micro%vapice)) &
!!$         call InsertVTab (micro%vapice,microm%vapice  &
!!$         ,ng, npts, imean,  &
!!$         'VAPICE :3:hist:anal:mpt3')
!!$    if (associated(micro%vapcld)) &
!!$         call InsertVTab (micro%vapcld,microm%vapcld  &
!!$         ,ng, npts, imean,  &
!!$         'VAPCLD :3:hist:anal:mpt3')
!!$    if (associated(micro%vaprain)) &
!!$         call InsertVTab (micro%vaprain,microm%vaprain  &
!!$         ,ng, npts, imean,  &
!!$         'VAPRAIN :3:hist:anal:mpt3')
!!$    if (associated(micro%vappris)) &
!!$         call InsertVTab (micro%vappris,microm%vappris  &
!!$         ,ng, npts, imean,  &
!!$         'VAPPRIS :3:hist:anal:mpt3')
!!$    if (associated(micro%vapsnow)) &
!!$         call InsertVTab (micro%vapsnow,microm%vapsnow  &
!!$         ,ng, npts, imean,  &
!!$         'VAPSNOW :3:hist:anal:mpt3')
!!$    if (associated(micro%vapaggr)) &
!!$         call InsertVTab (micro%vapaggr,microm%vapaggr  &
!!$         ,ng, npts, imean,  &
!!$         'VAPAGGR :3:hist:anal:mpt3')
!!$    if (associated(micro%vapgrau)) &
!!$         call InsertVTab (micro%vapgrau,microm%vapgrau  &
!!$         ,ng, npts, imean,  &
!!$         'VAPGRAU :3:hist:anal:mpt3')
!!$    if (associated(micro%vaphail)) &
!!$         call InsertVTab (micro%vaphail,microm%vaphail  &
!!$         ,ng, npts, imean,  &
!!$         'VAPHAIL :3:hist:anal:mpt3')
!!$    if (associated(micro%vapdriz)) &
!!$         call InsertVTab (micro%vapdriz,microm%vapdriz  &
!!$         ,ng, npts, imean,  &
!!$         'VAPDRIZ :3:hist:anal:mpt3')
!!$
!!$    if (associated(micro%meltice)) &
!!$         call InsertVTab (micro%meltice,microm%meltice  &
!!$         ,ng, npts, imean,  &
!!$         'MELTICE :3:hist:anal:mpt3')
!!$    if (associated(micro%meltpris)) &
!!$         call InsertVTab (micro%meltpris,microm%meltpris  &
!!$         ,ng, npts, imean,  &
!!$         'MELTPRIS :3:hist:anal:mpt3')
!!$    if (associated(micro%meltsnow)) &
!!$         call InsertVTab (micro%meltsnow,microm%meltsnow  &
!!$         ,ng, npts, imean,  &
!!$         'MELTSNOW :3:hist:anal:mpt3')
!!$    if (associated(micro%meltaggr)) &
!!$         call InsertVTab (micro%meltaggr,microm%meltaggr  &
!!$         ,ng, npts, imean,  &
!!$         'MELTAGGR :3:hist:anal:mpt3')
!!$    if (associated(micro%meltgrau)) &
!!$         call InsertVTab (micro%meltgrau,microm%meltgrau  &
!!$         ,ng, npts, imean,  &
!!$         'MELTGRAU :3:hist:anal:mpt3')
!!$    if (associated(micro%melthail)) &
!!$         call InsertVTab (micro%melthail,microm%melthail  &
!!$         ,ng, npts, imean,  &
!!$         'MELTHAIL :3:hist:anal:mpt3')
!!$
!!$    if (associated(micro%cld2rain)) &
!!$         call InsertVTab (micro%cld2rain,microm%cld2rain  &
!!$         ,ng, npts, imean,  &
!!$         'CLD2RAIN :3:hist:anal:mpt3')
!!$    if (associated(micro%rimecld)) &
!!$         call InsertVTab (micro%rimecld,microm%rimecld  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLD :3:hist:anal:mpt3')
!!$    if (associated(micro%rimecldsnow)) &
!!$         call InsertVTab (micro%rimecldsnow,microm%rimecldsnow  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLDSNOW :3:hist:anal:mpt3')
!!$    if (associated(micro%rimecldaggr)) &
!!$         call InsertVTab (micro%rimecldaggr,microm%rimecldaggr  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLDAGGR :3:hist:anal:mpt3')
!!$    if (associated(micro%rimecldgrau)) &
!!$         call InsertVTab (micro%rimecldgrau,microm%rimecldgrau  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLDGRAU :3:hist:anal:mpt3')
!!$    if (associated(micro%rimecldhail)) &
!!$         call InsertVTab (micro%rimecldhail,microm%rimecldhail  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLDHAIL :3:hist:anal:mpt3')
!!$
!!$    if (associated(micro%rain2ice)) &
!!$         call InsertVTab (micro%rain2ice,microm%rain2ice  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2ICE :3:hist:anal:mpt3')
!!$    if (associated(micro%rain2pr)) &
!!$         call InsertVTab (micro%rain2pr,microm%rain2pr  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2PR :3:hist:anal:mpt3')
!!$    if (associated(micro%rain2sn)) &
!!$         call InsertVTab (micro%rain2sn,microm%rain2sn  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2SN :3:hist:anal:mpt3')
!!$    if (associated(micro%rain2ag)) &
!!$         call InsertVTab (micro%rain2ag,microm%rain2ag  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2AG :3:hist:anal:mpt3')
!!$    if (associated(micro%rain2gr)) &
!!$         call InsertVTab (micro%rain2gr,microm%rain2gr  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2GR :3:hist:anal:mpt3')
!!$    if (associated(micro%rain2ha)) &
!!$         call InsertVTab (micro%rain2ha,microm%rain2ha  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2HA :3:hist:anal:mpt3')
!!$    if (associated(micro%rain2ha_xtra)) &
!!$         call InsertVTab (micro%rain2ha_xtra,microm%rain2ha_xtra  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2HA_XTRA :3:hist:anal:mpt3')
!!$    if (associated(micro%ice2rain)) &
!!$         call InsertVTab (micro%ice2rain,microm%ice2rain  &
!!$         ,ng, npts, imean,  &
!!$         'ICE2RAIN :3:hist:anal:mpt3')
!!$
!!$    if (associated(micro%aggregate)) &
!!$         call InsertVTab (micro%aggregate,microm%aggregate  &
!!$         ,ng, npts, imean,  &
!!$         'AGGREGATE :3:hist:anal:mpt3')
!!$    if (associated(micro%aggrselfpris)) &
!!$         call InsertVTab (micro%aggrselfpris,microm%aggrselfpris  &
!!$         ,ng, npts, imean,  &
!!$         'AGGRSELFPRIS :3:hist:anal:mpt3')
!!$    if (associated(micro%aggrselfsnow)) &
!!$         call InsertVTab (micro%aggrselfsnow,microm%aggrselfsnow  &
!!$         ,ng, npts, imean,  &
!!$         'AGGRSELFSNOW :3:hist:anal:mpt3')
!!$    if (associated(micro%aggrprissnow)) &
!!$         call InsertVTab (micro%aggrprissnow,microm%aggrprissnow  &
!!$         ,ng, npts, imean,  &
!!$         'AGGRPRISSNOW :3:hist:anal:mpt3')
!!$
!!$    if (associated(micro%latheatvap)) &
!!$         call InsertVTab (micro%latheatvap,microm%latheatvap  &
!!$         ,ng, npts, imean,  &
!!$         'LATHEATVAP :3:hist:anal:mpt3')
!!$    if (associated(micro%latheatfrz)) &
!!$         call InsertVTab (micro%latheatfrz,microm%latheatfrz  &
!!$         ,ng, npts, imean,  &
!!$         'LATHEATFRZ :3:hist:anal:mpt3')
!!$    !END MICRO BUDGET PROCESSES (instantaneous)
!!$
!!$    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
!!$    if (associated(micro%nuccldrt)) &
!!$         call InsertVTab (micro%nuccldrt,microm%nuccldrt  &
!!$         ,ng, npts, imean,  &
!!$         'NUCCLDRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%nuccldct)) &
!!$         call InsertVTab (micro%nuccldct,microm%nuccldct  &
!!$         ,ng, npts, imean,  &
!!$         'NUCCLDCT :3:hist:anal:mpti:mpt3')
!!$
!!$    if (associated(micro%nucicert)) &
!!$         call InsertVTab (micro%nucicert,microm%nucicert  &
!!$         ,ng, npts, imean,  &
!!$         'NUCICERT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%nucicect)) &
!!$         call InsertVTab (micro%nucicect,microm%nucicect  &
!!$         ,ng, npts, imean,  &
!!$         'NUCICECT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%inuchomrt)) &
!!$         call InsertVTab (micro%inuchomrt,microm%inuchomrt  &
!!$         ,ng, npts, imean,  &
!!$         'INUCHOMRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%inuchomct)) &
!!$         call InsertVTab (micro%inuchomct,microm%inuchomct  &
!!$         ,ng, npts, imean,  &
!!$         'INUCHOMCT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%inuccontrt)) &
!!$         call InsertVTab (micro%inuccontrt,microm%inuccontrt  &
!!$         ,ng, npts, imean,  &
!!$         'INUCCONTRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%inuccontct)) &
!!$         call InsertVTab (micro%inuccontct,microm%inuccontct  &
!!$         ,ng, npts, imean,  &
!!$         'INUCCONTCT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%inucifnrt)) &
!!$         call InsertVTab (micro%inucifnrt,microm%inucifnrt  &
!!$         ,ng, npts, imean,  &
!!$         'INUCIFNRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%inucifnct)) &
!!$         call InsertVTab (micro%inucifnct,microm%inucifnct  &
!!$         ,ng, npts, imean,  &
!!$         'INUCIFNCT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%inuchazrt)) &
!!$         call InsertVTab (micro%inuchazrt,microm%inuchazrt  &
!!$         ,ng, npts, imean,  &
!!$         'INUCHAZRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%inuchazct)) &
!!$         call InsertVTab (micro%inuchazct,microm%inuchazct  &
!!$         ,ng, npts, imean,  &
!!$         'INUCHAZCT :3:hist:anal:mpti:mpt3')
!!$
!!$    if (associated(micro%vapliqt)) &
!!$         call InsertVTab (micro%vapliqt,microm%vapliqt  &
!!$         ,ng, npts, imean,  &
!!$         'VAPLIQT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%vapicet)) &
!!$         call InsertVTab (micro%vapicet,microm%vapicet  &
!!$         ,ng, npts, imean,  &
!!$         'VAPICET :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%vapcldt)) &
!!$         call InsertVTab (micro%vapcldt,microm%vapcldt  &
!!$         ,ng, npts, imean,  &
!!$         'VAPCLDT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%vapraint)) &
!!$         call InsertVTab (micro%vapraint,microm%vapraint  &
!!$         ,ng, npts, imean,  &
!!$         'VAPRAINT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%vapprist)) &
!!$         call InsertVTab (micro%vapprist,microm%vapprist  &
!!$         ,ng, npts, imean,  &
!!$         'VAPPRIST :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%vapsnowt)) &
!!$         call InsertVTab (micro%vapsnowt,microm%vapsnowt  &
!!$         ,ng, npts, imean,  &
!!$         'VAPSNOWT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%vapaggrt)) &
!!$         call InsertVTab (micro%vapaggrt,microm%vapaggrt  &
!!$         ,ng, npts, imean,  &
!!$         'VAPAGGRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%vapgraut)) &
!!$         call InsertVTab (micro%vapgraut,microm%vapgraut  &
!!$         ,ng, npts, imean,  &
!!$         'VAPGRAUT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%vaphailt)) &
!!$         call InsertVTab (micro%vaphailt,microm%vaphailt  &
!!$         ,ng, npts, imean,  &
!!$         'VAPHAILT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%vapdrizt)) &
!!$         call InsertVTab (micro%vapdrizt,microm%vapdrizt  &
!!$         ,ng, npts, imean,  &
!!$         'VAPDRIZT :3:hist:anal:mpti:mpt3')
!!$
!!$    if (associated(micro%melticet)) &
!!$         call InsertVTab (micro%melticet,microm%melticet  &
!!$         ,ng, npts, imean,  &
!!$         'MELTICET :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%meltprist)) &
!!$         call InsertVTab (micro%meltprist,microm%meltprist  &
!!$         ,ng, npts, imean,  &
!!$         'MELTPRIST :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%meltsnowt)) &
!!$         call InsertVTab (micro%meltsnowt,microm%meltsnowt  &
!!$         ,ng, npts, imean,  &
!!$         'MELTSNOWT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%meltaggrt)) &
!!$         call InsertVTab (micro%meltaggrt,microm%meltaggrt  &
!!$         ,ng, npts, imean,  &
!!$         'MELTAGGRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%meltgraut)) &
!!$         call InsertVTab (micro%meltgraut,microm%meltgraut  &
!!$         ,ng, npts, imean,  &
!!$         'MELTGRAUT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%melthailt)) &
!!$         call InsertVTab (micro%melthailt,microm%melthailt  &
!!$         ,ng, npts, imean,  &
!!$         'MELTHAILT :3:hist:anal:mpti:mpt3')
!!$
!!$    if (associated(micro%cld2raint)) &
!!$         call InsertVTab (micro%cld2raint,microm%cld2raint  &
!!$         ,ng, npts, imean,  &
!!$         'CLD2RAINT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rimecldt)) &
!!$         call InsertVTab (micro%rimecldt,microm%rimecldt  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLDT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rimecldsnowt)) &
!!$         call InsertVTab (micro%rimecldsnowt,microm%rimecldsnowt  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLDSNOWT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rimecldaggrt)) &
!!$         call InsertVTab (micro%rimecldaggrt,microm%rimecldaggrt  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLDAGGRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rimecldgraut)) &
!!$         call InsertVTab (micro%rimecldgraut,microm%rimecldgraut  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLDGRAUT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rimecldhailt)) &
!!$         call InsertVTab (micro%rimecldhailt,microm%rimecldhailt  &
!!$         ,ng, npts, imean,  &
!!$         'RIMECLDHAILT :3:hist:anal:mpti:mpt3')
!!$
!!$    if (associated(micro%rain2icet)) &
!!$         call InsertVTab (micro%rain2icet,microm%rain2icet  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2ICET :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rain2prt)) &
!!$         call InsertVTab (micro%rain2prt,microm%rain2prt  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2PRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rain2snt)) &
!!$         call InsertVTab (micro%rain2snt,microm%rain2snt  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2SNT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rain2agt)) &
!!$         call InsertVTab (micro%rain2agt,microm%rain2agt  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2AGT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rain2grt)) &
!!$         call InsertVTab (micro%rain2grt,microm%rain2grt  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2GRT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rain2hat)) &
!!$         call InsertVTab (micro%rain2hat,microm%rain2hat  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2HAT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%rain2ha_xtrat)) &
!!$         call InsertVTab (micro%rain2ha_xtrat,microm%rain2ha_xtrat  &
!!$         ,ng, npts, imean,  &
!!$         'RAIN2HA_XTRAT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%ice2raint)) &
!!$         call InsertVTab (micro%ice2raint,microm%ice2raint  &
!!$         ,ng, npts, imean,  &
!!$         'ICE2RAINT :3:hist:anal:mpti:mpt3')
!!$
!!$    if (associated(micro%aggregatet)) &
!!$         call InsertVTab (micro%aggregatet,microm%aggregatet  &
!!$         ,ng, npts, imean,  &
!!$         'AGGREGATET :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%aggrselfprist)) &
!!$         call InsertVTab (micro%aggrselfprist,microm%aggrselfprist  &
!!$         ,ng, npts, imean,  &
!!$         'AGGRSELFPRIST :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%aggrselfsnowt)) &
!!$         call InsertVTab (micro%aggrselfsnowt,microm%aggrselfsnowt  &
!!$         ,ng, npts, imean,  &
!!$         'AGGRSELFSNOWT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%aggrprissnowt)) &
!!$         call InsertVTab (micro%aggrprissnowt,microm%aggrprissnowt  &
!!$         ,ng, npts, imean,  &
!!$         'AGGRPRISSNOWT :3:hist:anal:mpti:mpt3')
!!$
!!$    if (associated(micro%latheatvapt)) &
!!$         call InsertVTab (micro%latheatvapt,microm%latheatvapt  &
!!$         ,ng, npts, imean,  &
!!$         'LATHEATVAPT :3:hist:anal:mpti:mpt3')
!!$    if (associated(micro%latheatfrzt)) &
!!$         call InsertVTab (micro%latheatfrzt,microm%latheatfrzt  &
!!$         ,ng, npts, imean,  &
!!$         'LATHEATFRZT :3:hist:anal:mpti:mpt3')
!!$    !END MICRO BUDGET PROCECCES (totals)
!!$
!!$    npts=n2*n3
!!$    if (associated(micro%accpr)) &
!!$         call InsertVTab (micro%accpr,microm%accpr  &
!!$         ,ng, npts, imean,  &
!!$         'ACCPR :2:hist:anal:mpti:mpt3')
!!$    if (associated(micro%accpp)) &
!!$         call InsertVTab (micro%accpp,microm%accpp  &
!!$         ,ng, npts, imean,  &
!!$         'ACCPP :2:hist:anal:mpti:mpt3')
!!$    if (associated(micro%accps)) &
!!$         call InsertVTab (micro%accps,microm%accps  &
!!$         ,ng, npts, imean,  &
!!$         'ACCPS :2:hist:anal:mpti:mpt3')
!!$    if (associated(micro%accpa)) &
!!$         call InsertVTab (micro%accpa,microm%accpa  &
!!$         ,ng, npts, imean,  &
!!$         'ACCPA :2:hist:anal:mpti:mpt3')
!!$    if (associated(micro%accpg)) &
!!$         call InsertVTab (micro%accpg,microm%accpg  &
!!$         ,ng, npts, imean,  &
!!$         'ACCPG :2:hist:anal:mpti:mpt3')
!!$    if (associated(micro%accph)) &
!!$         call InsertVTab (micro%accph,microm%accph  &
!!$         ,ng, npts, imean,  &
!!$         'ACCPH :2:hist:anal:mpti:mpt3')
!!$    if (associated(micro%accpd)) &
!!$         call InsertVTab (micro%accpd,microm%accpd  &
!!$         ,ng, npts, imean,  &
!!$         'ACCPD :2:hist:anal:mpti:mpt3')
!!$    if (associated(micro%pcprr)) &
!!$         call InsertVTab (micro%pcprr,microm%pcprr  &
!!$         ,ng, npts, imean,  &
!!$         'PCPRR :2:hist:anal:mpt3')
!!$    if (associated(micro%pcprp)) &
!!$         call InsertVTab (micro%pcprp,microm%pcprp  &
!!$         ,ng, npts, imean,  &
!!$         'PCPRP :2:hist:anal:mpt3')
!!$    if (associated(micro%pcprs)) &
!!$         call InsertVTab (micro%pcprs,microm%pcprs  &
!!$         ,ng, npts, imean,  &
!!$         'PCPRS :2:hist:anal:mpt3')
!!$    if (associated(micro%pcpra)) &
!!$         call InsertVTab (micro%pcpra,microm%pcpra  &
!!$         ,ng, npts, imean,  &
!!$         'PCPRA :2:hist:anal:mpt3')
!!$    if (associated(micro%pcprg)) &
!!$         call InsertVTab (micro%pcprg,microm%pcprg  &
!!$         ,ng, npts, imean,  &
!!$         'PCPRG :2:hist:anal:mpt3')
!!$    if (associated(micro%pcprh)) &
!!$         call InsertVTab (micro%pcprh,microm%pcprh  &
!!$         ,ng, npts, imean,  &
!!$         'PCPRH :2:hist:anal:mpt3')
!!$    if (associated(micro%pcprd)) &
!!$         call InsertVTab (micro%pcprd,microm%pcprd  &
!!$         ,ng, npts, imean,  &
!!$         'PCPRD :2:hist:anal:mpt3')
!!$    if (associated(micro%pcpg)) &
!!$         call InsertVTab (micro%pcpg,microm%pcpg  &
!!$         ,ng, npts, imean,  &
!!$         'PCPG :2:hist:mpti:mpt3')
!!$    if (associated(micro%qpcpg)) &
!!$         call InsertVTab (micro%qpcpg,microm%qpcpg  &
!!$         ,ng, npts, imean,  &
!!$         'QPCPG :2:hist:mpti:mpt3')
!!$    if (associated(micro%dpcpg)) &
!!$         call InsertVTab (micro%dpcpg,microm%dpcpg  &
!!$         ,ng, npts, imean,  &
!!$         'DPCPG :2:hist:mpti:mpt3')
!!$
!!$    return
!!$  end subroutine filltab_micro

  subroutine DestroyMicroFields(oneMicroFields)
    type(MicroFields), pointer, intent(inout) :: oneMicroFields

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyMicroFields)**"

    if (.not. associated(oneMicroFields)) then
       return
    end if
    
    if (associated(oneMicroFields%rcp)) then
       deallocate(oneMicroFields%rcp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rcp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rdp)) then
       deallocate(oneMicroFields%rdp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rdp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rrp)) then
       deallocate(oneMicroFields%rrp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rrp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rpp)) then
       deallocate(oneMicroFields%rpp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rpp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rsp)) then
       deallocate(oneMicroFields%rsp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rsp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rap)) then
       deallocate(oneMicroFields%rap, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rap fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rgp)) then
       deallocate(oneMicroFields%rgp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rgp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rhp)) then
       deallocate(oneMicroFields%rhp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rhp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%ccp)) then
       deallocate(oneMicroFields%ccp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate ccp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cdp)) then
       deallocate(oneMicroFields%cdp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cdp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%crp)) then
       deallocate(oneMicroFields%crp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate crp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cpp)) then
       deallocate(oneMicroFields%cpp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cpp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%csp)) then
       deallocate(oneMicroFields%csp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate csp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cap)) then
       deallocate(oneMicroFields%cap, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cap fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cgp)) then
       deallocate(oneMicroFields%cgp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cgp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%chp)) then
       deallocate(oneMicroFields%chp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate chp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cccnp)) then
       deallocate(oneMicroFields%cccnp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cccnp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%gccnp)) then
       deallocate(oneMicroFields%gccnp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate gccnp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cifnp)) then
       deallocate(oneMicroFields%cifnp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cifnp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%q2)) then
       deallocate(oneMicroFields%q2, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate q2 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%q6)) then
       deallocate(oneMicroFields%q6, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate q6 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%q7)) then
       deallocate(oneMicroFields%q7, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate q7 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rei)) then
       deallocate(oneMicroFields%rei, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rei fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rel)) then
       deallocate(oneMicroFields%rel, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rel fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cldfr)) then
       deallocate(oneMicroFields%cldfr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cldfr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cccmp)) then
       deallocate(oneMicroFields%cccmp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cccmp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%gccmp)) then
       deallocate(oneMicroFields%gccmp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate gccmp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cnm1p)) then
       deallocate(oneMicroFields%cnm1p, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cnm1p fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cnm2p)) then
       deallocate(oneMicroFields%cnm2p, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cnm2p fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cnm3p)) then
       deallocate(oneMicroFields%cnm3p, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cnm3p fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cnm8p)) then
       deallocate(oneMicroFields%cnm8p, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cnm8p fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%md1np)) then
       deallocate(oneMicroFields%md1np, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate md1np fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%md2np)) then
       deallocate(oneMicroFields%md2np, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate md2np fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%salt_filmp)) then
       deallocate(oneMicroFields%salt_filmp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate salt_filmp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%salt_jetp)) then
       deallocate(oneMicroFields%salt_jetp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate salt_jetp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%salt_spmp)) then
       deallocate(oneMicroFields%salt_spmp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate salt_spmp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcpvr)) then
       deallocate(oneMicroFields%pcpvr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcpvr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcpvp)) then
       deallocate(oneMicroFields%pcpvp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcpvp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcpvs)) then
       deallocate(oneMicroFields%pcpvs, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcpvs fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcpva)) then
       deallocate(oneMicroFields%pcpva, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcpva fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcpvg)) then
       deallocate(oneMicroFields%pcpvg, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcpvg fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcpvh)) then
       deallocate(oneMicroFields%pcpvh, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcpvh fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcpvd)) then
       deallocate(oneMicroFields%pcpvd, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcpvd fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%nuccldr)) then
       deallocate(oneMicroFields%nuccldr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate nuccldr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%nuccldc)) then
       deallocate(oneMicroFields%nuccldc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate nuccldc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%nucicer)) then
       deallocate(oneMicroFields%nucicer, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate nucicer fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%nucicec)) then
       deallocate(oneMicroFields%nucicec, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate nucicec fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuchomr)) then
       deallocate(oneMicroFields%inuchomr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuchomr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuchomc)) then
       deallocate(oneMicroFields%inuchomc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuchomc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuccontr)) then
       deallocate(oneMicroFields%inuccontr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuccontr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuccontc)) then
       deallocate(oneMicroFields%inuccontc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuccontc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inucifnr)) then
       deallocate(oneMicroFields%inucifnr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inucifnr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inucifnc)) then
       deallocate(oneMicroFields%inucifnc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inucifnc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuchazr)) then
       deallocate(oneMicroFields%inuchazr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuchazr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuchazc)) then
       deallocate(oneMicroFields%inuchazc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuchazc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapliq)) then
       deallocate(oneMicroFields%vapliq, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapliq fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapice)) then
       deallocate(oneMicroFields%vapice, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapice fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapcld)) then
       deallocate(oneMicroFields%vapcld, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapcld fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vaprain)) then
       deallocate(oneMicroFields%vaprain, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vaprain fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vappris)) then
       deallocate(oneMicroFields%vappris, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vappris fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapsnow)) then
       deallocate(oneMicroFields%vapsnow, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapsnow fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapaggr)) then
       deallocate(oneMicroFields%vapaggr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapaggr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapgrau)) then
       deallocate(oneMicroFields%vapgrau, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapgrau fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vaphail)) then
       deallocate(oneMicroFields%vaphail, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vaphail fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapdriz)) then
       deallocate(oneMicroFields%vapdriz, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapdriz fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%meltice)) then
       deallocate(oneMicroFields%meltice, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate meltice fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%meltpris)) then
       deallocate(oneMicroFields%meltpris, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate meltpris fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%meltsnow)) then
       deallocate(oneMicroFields%meltsnow, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate meltsnow fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%meltaggr)) then
       deallocate(oneMicroFields%meltaggr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate meltaggr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%meltgrau)) then
       deallocate(oneMicroFields%meltgrau, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate meltgrau fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%melthail)) then
       deallocate(oneMicroFields%melthail, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate melthail fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cld2rain)) then
       deallocate(oneMicroFields%cld2rain, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cld2rain fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecld)) then
       deallocate(oneMicroFields%rimecld, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecld fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecldsnow)) then
       deallocate(oneMicroFields%rimecldsnow, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecldsnow fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecldaggr)) then
       deallocate(oneMicroFields%rimecldaggr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecldaggr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecldgrau)) then
       deallocate(oneMicroFields%rimecldgrau, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecldgrau fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecldhail)) then
       deallocate(oneMicroFields%rimecldhail, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecldhail fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2ice)) then
       deallocate(oneMicroFields%rain2ice, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2ice fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2pr)) then
       deallocate(oneMicroFields%rain2pr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2pr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2sn)) then
       deallocate(oneMicroFields%rain2sn, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2sn fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2ag)) then
       deallocate(oneMicroFields%rain2ag, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2ag fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2gr)) then
       deallocate(oneMicroFields%rain2gr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2gr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2ha)) then
       deallocate(oneMicroFields%rain2ha, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2ha fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2ha_xtra)) then
       deallocate(oneMicroFields%rain2ha_xtra, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2ha_xtra fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%ice2rain)) then
       deallocate(oneMicroFields%ice2rain, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate ice2rain fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%aggregate)) then
       deallocate(oneMicroFields%aggregate, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate aggregate fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%aggrselfpris)) then
       deallocate(oneMicroFields%aggrselfpris, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate aggrselfpris fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%aggrselfsnow)) then
       deallocate(oneMicroFields%aggrselfsnow, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate aggrselfsnow fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%aggrprissnow)) then
       deallocate(oneMicroFields%aggrprissnow, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate aggrprissnow fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%latheatvap)) then
       deallocate(oneMicroFields%latheatvap, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate latheatvap fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%latheatfrz)) then
       deallocate(oneMicroFields%latheatfrz, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate latheatfrz fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%nuccldrt)) then
       deallocate(oneMicroFields%nuccldrt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate nuccldrt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%nuccldct)) then
       deallocate(oneMicroFields%nuccldct, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate nuccldct fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%nucicert)) then
       deallocate(oneMicroFields%nucicert, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate nucicert fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%nucicect)) then
       deallocate(oneMicroFields%nucicect, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate nucicect fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuchomrt)) then
       deallocate(oneMicroFields%inuchomrt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuchomrt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuchomct)) then
       deallocate(oneMicroFields%inuchomct, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuchomct fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuccontrt)) then
       deallocate(oneMicroFields%inuccontrt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuccontrt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuccontct)) then
       deallocate(oneMicroFields%inuccontct, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuccontct fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inucifnrt)) then
       deallocate(oneMicroFields%inucifnrt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inucifnrt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inucifnct)) then
       deallocate(oneMicroFields%inucifnct, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inucifnct fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuchazrt)) then
       deallocate(oneMicroFields%inuchazrt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuchazrt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%inuchazct)) then
       deallocate(oneMicroFields%inuchazct, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate inuchazct fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapliqt)) then
       deallocate(oneMicroFields%vapliqt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapliqt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapicet)) then
       deallocate(oneMicroFields%vapicet, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapicet fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapcldt)) then
       deallocate(oneMicroFields%vapcldt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapcldt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapraint)) then
       deallocate(oneMicroFields%vapraint, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapraint fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapprist)) then
       deallocate(oneMicroFields%vapprist, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapprist fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapsnowt)) then
       deallocate(oneMicroFields%vapsnowt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapsnowt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapaggrt)) then
       deallocate(oneMicroFields%vapaggrt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapaggrt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapgraut)) then
       deallocate(oneMicroFields%vapgraut, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapgraut fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vaphailt)) then
       deallocate(oneMicroFields%vaphailt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vaphailt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%vapdrizt)) then
       deallocate(oneMicroFields%vapdrizt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vapdrizt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%melticet)) then
       deallocate(oneMicroFields%melticet, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate melticet fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%meltprist)) then
       deallocate(oneMicroFields%meltprist, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate meltprist fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%meltsnowt)) then
       deallocate(oneMicroFields%meltsnowt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate meltsnowt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%meltaggrt)) then
       deallocate(oneMicroFields%meltaggrt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate meltaggrt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%meltgraut)) then
       deallocate(oneMicroFields%meltgraut, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate meltgraut fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%melthailt)) then
       deallocate(oneMicroFields%melthailt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate melthailt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%cld2raint)) then
       deallocate(oneMicroFields%cld2raint, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cld2raint fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecldt)) then
       deallocate(oneMicroFields%rimecldt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecldt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecldsnowt)) then
       deallocate(oneMicroFields%rimecldsnowt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecldsnowt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecldaggrt)) then
       deallocate(oneMicroFields%rimecldaggrt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecldaggrt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecldgraut)) then
       deallocate(oneMicroFields%rimecldgraut, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecldgraut fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rimecldhailt)) then
       deallocate(oneMicroFields%rimecldhailt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rimecldhailt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2icet)) then
       deallocate(oneMicroFields%rain2icet, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2icet fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2prt)) then
       deallocate(oneMicroFields%rain2prt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2prt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2snt)) then
       deallocate(oneMicroFields%rain2snt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2snt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2agt)) then
       deallocate(oneMicroFields%rain2agt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2agt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2grt)) then
       deallocate(oneMicroFields%rain2grt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2grt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2hat)) then
       deallocate(oneMicroFields%rain2hat, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2hat fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%rain2ha_xtrat)) then
       deallocate(oneMicroFields%rain2ha_xtrat, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2ha_xtrat fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%ice2raint)) then
       deallocate(oneMicroFields%ice2raint, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate ice2raint fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%aggregatet)) then
       deallocate(oneMicroFields%aggregatet, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate aggregatet fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%aggrselfprist)) then
       deallocate(oneMicroFields%aggrselfprist, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate aggrselfprist fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%aggrselfsnowt)) then
       deallocate(oneMicroFields%aggrselfsnowt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate aggrselfsnowt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%aggrprissnowt)) then
       deallocate(oneMicroFields%aggrprissnowt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate aggrprissnowt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%latheatvapt)) then
       deallocate(oneMicroFields%latheatvapt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate latheatvapt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%latheatfrzt)) then
       deallocate(oneMicroFields%latheatfrzt, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate latheatfrzt fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%accpr)) then
       deallocate(oneMicroFields%accpr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate accpr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%accpp)) then
       deallocate(oneMicroFields%accpp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate accpp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%accps)) then
       deallocate(oneMicroFields%accps, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate accps fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%accpa)) then
       deallocate(oneMicroFields%accpa, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate accpa fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%accpg)) then
       deallocate(oneMicroFields%accpg, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate accpg fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%accph)) then
       deallocate(oneMicroFields%accph, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate accph fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%accpd)) then
       deallocate(oneMicroFields%accpd, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate accpd fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcprr)) then
       deallocate(oneMicroFields%pcprr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcprr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcprp)) then
       deallocate(oneMicroFields%pcprp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcprp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcprs)) then
       deallocate(oneMicroFields%pcprs, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcprs fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcpra)) then
       deallocate(oneMicroFields%pcpra, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcpra fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcprg)) then
       deallocate(oneMicroFields%pcprg, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcprg fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcprh)) then
       deallocate(oneMicroFields%pcprh, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcprh fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcprd)) then
       deallocate(oneMicroFields%pcprd, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcprd fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%pcpg)) then
       deallocate(oneMicroFields%pcpg, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pcpg fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%qpcpg)) then
       deallocate(oneMicroFields%qpcpg, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate qpcpg fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneMicroFields%dpcpg)) then
       deallocate(oneMicroFields%dpcpg, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate dpcpg fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    deallocate(oneMicroFields, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneMicroFields fails with stat="//&
            trim(adjustl(str(1))))
    end if
    nullify(oneMicroFields)
  end subroutine DestroyMicroFields


  subroutine DumpMicroFields(oneMicroFields, name)
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpMicroFields)**"

    if (associated(oneMicroFields)) then
       call MsgDump(h//" "//name//" associated")
    else
       call MsgDump(h//" "//name//" not associated")
    end if
  end subroutine DumpMicroFields


 
  subroutine InsertMicroFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneMicroFields, oneAveMicroFields, imean)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    type(MicroFields), pointer, intent(in) :: oneAveMicroFields
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(InsertMicroFieldsAtVarTable)**"

    if (.not. associated(oneMicroFields)) then
       call fatal_error(h//" oneMicroFields not associated")
    else if (.not. associated(oneAveMicroFields)) then
       call fatal_error(h//" oneAveMicroFields not associated")
    else if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    end if

    ! Fill pointers to arrays into variable tables

    if (associated(oneMicroFields%rcp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rcp, &
            'RCP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%rcp, imean)
    end if

    if (associated(oneMicroFields%rdp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rdp, &
            'RDP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%rdp, imean)
    end if

    if (associated(oneMicroFields%rrp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rrp, &
            'RRP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%rrp, imean)
    end if

    if (associated(oneMicroFields%rpp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rpp, &
            'RPP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%rpp, imean)
    end if

    if (associated(oneMicroFields%rsp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rsp, &
            'RSP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%rsp, imean)
    end if

    if (associated(oneMicroFields%rap)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rap, &
            'RAP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%rap, imean)
    end if

    if (associated(oneMicroFields%rgp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rgp, &
            'RGP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%rgp, imean)
    end if

    if (associated(oneMicroFields%rhp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rhp, &
            'RHP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%rhp, imean)
    end if

    if (associated(oneMicroFields%ccp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%ccp, &
            'CCP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%ccp, imean)
    end if

    if (associated(oneMicroFields%cdp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cdp, &
            'CDP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cdp, imean)
    end if

    if (associated(oneMicroFields%crp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%crp, &
            'CRP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%crp, imean)
    end if

    if (associated(oneMicroFields%cpp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cpp, &
            'CPP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cpp, imean)
    end if

    if (associated(oneMicroFields%csp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%csp, &
            'CSP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%csp, imean)
    end if

    if (associated(oneMicroFields%cap)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cap, &
            'CAP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cap, imean)
    end if

    if (associated(oneMicroFields%cgp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cgp, &
            'CGP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cgp, imean)
    end if

    if (associated(oneMicroFields%chp)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%chp, &
            'CHP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%chp, imean)
    end if

    if (associated(oneMicroFields%cccnp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cccnp, &
            'CCCNP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cccnp, imean)
    end if

    if (associated(oneMicroFields%gccnp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%gccnp, &
            'GCCNP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%gccnp, imean)
    end if

    if (associated(oneMicroFields%cifnp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cifnp, &
            'CIFNP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cifnp, imean)
    end if


    if (associated(oneMicroFields%q2)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%q2, &
            'Q2 :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%q2, imean)
    end if

    if (associated(oneMicroFields%q6)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%q6, &
            'Q6 :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%q6, imean)
    end if

    if (associated(oneMicroFields%q7)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%q7, &
            'Q7 :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%q7, imean)
    end if

    if (associated(oneMicroFields%cccmp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cccmp, &
            'CCCMP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cccmp, imean)
    end if

    if (associated(oneMicroFields%gccmp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%gccmp, &
            'GCCMP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%gccmp, imean)
    end if

    if (associated(oneMicroFields%cnm1p)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cnm1p, &
            'CNM1P :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cnm1p, imean)
    end if

    if (associated(oneMicroFields%cnm2p)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cnm2p, &
            'CNM2P :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cnm2p, imean)
    end if

    if (associated(oneMicroFields%cnm3p)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cnm3p, &
            'CNM3P :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cnm3p, imean)
    end if

    if (associated(oneMicroFields%cnm8p)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cnm8p, &
            'CNM8P :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%cnm8p, imean)
    end if

    if (associated(oneMicroFields%md1np)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%md1np, &
            'MD1NP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%md1np, imean)
    end if

    if (associated(oneMicroFields%md2np)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%md2np, &
            'MD2NP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%md2np, imean)
    end if

    if (associated(oneMicroFields%salt_filmp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%salt_filmp, &
            'SALT_FILMP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%salt_filmp, imean)
    end if

    if (associated(oneMicroFields%salt_jetp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%salt_jetp, &
            'SALT_JETP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%salt_jetp, imean)
    end if

    if (associated(oneMicroFields%salt_spmp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%salt_spmp, &
            'SALT_SPMP :3:hist:anal:mpti:mpt3:mpt1', &
            oneAveMicroFields%salt_spmp, imean)
    end if

    if (associated(oneMicroFields%rei)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rei, &
            'REI :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rei, imean)
    end if

    if (associated(oneMicroFields%rel)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rel, &
            'REL :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rel, imean)
    end if

    if (associated(oneMicroFields%cldfr)) then   
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cldfr, &
            'CLDFR :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%cldfr, imean)
    end if

    !VERTICAL PRECIPITATION RATES
    if (associated(oneMicroFields%pcpvr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcpvr, &
            'PCPVR :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%pcpvr, imean)
    end if

    if (associated(oneMicroFields%pcpvp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcpvp, &
            'PCPVP :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%pcpvp, imean)
    end if

    if (associated(oneMicroFields%pcpvs)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcpvs, &
            'PCPVS :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%pcpvs, imean)
    end if

    if (associated(oneMicroFields%pcpva)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcpva, &
            'PCPVA :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%pcpva, imean)
    end if

    if (associated(oneMicroFields%pcpvg)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcpvg, &
            'PCPVG :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%pcpvg, imean)
    end if

    if (associated(oneMicroFields%pcpvh)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcpvh, &
            'PCPVH :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%pcpvh, imean)
    end if

    if (associated(oneMicroFields%pcpvd)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcpvd, &
            'PCPVD :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%pcpvd, imean)
    end if


    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (instantaneous)
    if (associated(oneMicroFields%nuccldr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%nuccldr, &
            'NUCCLDR :3:hist:anal:mpt3', &
            oneAveMicroFields%nuccldr, imean)
    end if

    if (associated(oneMicroFields%nuccldc)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%nuccldc, &
            'NUCCLDC :3:hist:anal:mpt3', &
            oneAveMicroFields%nuccldc, imean)
    end if

    if (associated(oneMicroFields%nucicer)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%nucicer, &
            'NUCICER :3:hist:anal:mpt3', &
            oneAveMicroFields%nucicer, imean)
    end if

    if (associated(oneMicroFields%nucicec)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%nucicec, &
            'NUCICEC :3:hist:anal:mpt3', &
            oneAveMicroFields%nucicec, imean)
    end if

    if (associated(oneMicroFields%inuchomr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuchomr, &
            'INUCHOMR :3:hist:anal:mpt3', &
            oneAveMicroFields%inuchomr, imean)
    end if

    if (associated(oneMicroFields%inuchomc)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuchomc, &
            'INUCHOMC :3:hist:anal:mpt3', &
            oneAveMicroFields%inuchomc, imean)
    end if

    if (associated(oneMicroFields%inuccontr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuccontr, &
            'INUCCONTR :3:hist:anal:mpt3', &
            oneAveMicroFields%inuccontr, imean)
    end if

    if (associated(oneMicroFields%inuccontc)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuccontc, &
            'INUCCONTC :3:hist:anal:mpt3', &
            oneAveMicroFields%inuccontc, imean)
    end if

    if (associated(oneMicroFields%inucifnr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inucifnr, &
            'INUCIFNR :3:hist:anal:mpt3', &
            oneAveMicroFields%inucifnr, imean)
    end if

    if (associated(oneMicroFields%inucifnc)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inucifnc, &
            'INUCIFNC :3:hist:anal:mpt3', &
            oneAveMicroFields%inucifnc, imean)
    end if

    if (associated(oneMicroFields%inuchazr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuchazr, &
            'INUCHAZR :3:hist:anal:mpt3', &
            oneAveMicroFields%inuchazr, imean)
    end if

    if (associated(oneMicroFields%inuchazc)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuchazc, &
            'INUCHAZC :3:hist:anal:mpt3', &
            oneAveMicroFields%inuchazc, imean)
    end if

    if (associated(oneMicroFields%vapliq)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapliq, &
            'VAPLIQ :3:hist:anal:mpt3', &
            oneAveMicroFields%vapliq, imean)
    end if

    if (associated(oneMicroFields%vapice)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapice, &
            'VAPICE :3:hist:anal:mpt3', &
            oneAveMicroFields%vapice, imean)
    end if

    if (associated(oneMicroFields%vapcld)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapcld, &
            'VAPCLD :3:hist:anal:mpt3', &
            oneAveMicroFields%vapcld, imean)
    end if

    if (associated(oneMicroFields%vaprain)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vaprain, &
            'VAPRAIN :3:hist:anal:mpt3', &
            oneAveMicroFields%vaprain, imean)
    end if

    if (associated(oneMicroFields%vappris)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vappris, &
            'VAPPRIS :3:hist:anal:mpt3', &
            oneAveMicroFields%vappris, imean)
    end if

    if (associated(oneMicroFields%vapsnow)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapsnow, &
            'VAPSNOW :3:hist:anal:mpt3', &
            oneAveMicroFields%vapsnow, imean)
    end if

    if (associated(oneMicroFields%vapaggr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapaggr, &
            'VAPAGGR :3:hist:anal:mpt3', &
            oneAveMicroFields%vapaggr, imean)
    end if

    if (associated(oneMicroFields%vapgrau)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapgrau, &
            'VAPGRAU :3:hist:anal:mpt3', &
            oneAveMicroFields%vapgrau, imean)
    end if

    if (associated(oneMicroFields%vaphail)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vaphail, &
            'VAPHAIL :3:hist:anal:mpt3', &
            oneAveMicroFields%vaphail, imean)
    end if

    if (associated(oneMicroFields%vapdriz)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapdriz, &
            'VAPDRIZ :3:hist:anal:mpt3', &
            oneAveMicroFields%vapdriz, imean)
    end if


    if (associated(oneMicroFields%meltice)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%meltice, &
            'MELTICE :3:hist:anal:mpt3', &
            oneAveMicroFields%meltice, imean)
    end if

    if (associated(oneMicroFields%meltpris)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%meltpris, &
            'MELTPRIS :3:hist:anal:mpt3', &
            oneAveMicroFields%meltpris, imean)
    end if

    if (associated(oneMicroFields%meltsnow)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%meltsnow, &
            'MELTSNOW :3:hist:anal:mpt3', &
            oneAveMicroFields%meltsnow, imean)
    end if

    if (associated(oneMicroFields%meltaggr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%meltaggr, &
            'MELTAGGR :3:hist:anal:mpt3', &
            oneAveMicroFields%meltaggr, imean)
    end if

    if (associated(oneMicroFields%meltgrau)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%meltgrau, &
            'MELTGRAU :3:hist:anal:mpt3', &
            oneAveMicroFields%meltgrau, imean)
    end if

    if (associated(oneMicroFields%melthail)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%melthail, &
            'MELTHAIL :3:hist:anal:mpt3', &
            oneAveMicroFields%melthail, imean)
    end if


    if (associated(oneMicroFields%cld2rain)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cld2rain, &
            'CLD2RAIN :3:hist:anal:mpt3', &
            oneAveMicroFields%cld2rain, imean)
    end if

    if (associated(oneMicroFields%rimecld)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecld, &
            'RIMECLD :3:hist:anal:mpt3', &
            oneAveMicroFields%rimecld, imean)
    end if

    if (associated(oneMicroFields%rimecldsnow)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecldsnow, &
            'RIMECLDSNOW :3:hist:anal:mpt3', &
            oneAveMicroFields%rimecldsnow, imean)
    end if

    if (associated(oneMicroFields%rimecldaggr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecldaggr, &
            'RIMECLDAGGR :3:hist:anal:mpt3', &
            oneAveMicroFields%rimecldaggr, imean)
    end if

    if (associated(oneMicroFields%rimecldgrau)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecldgrau, &
            'RIMECLDGRAU :3:hist:anal:mpt3', &
            oneAveMicroFields%rimecldgrau, imean)
    end if

    if (associated(oneMicroFields%rimecldhail)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecldhail, &
            'RIMECLDHAIL :3:hist:anal:mpt3', &
            oneAveMicroFields%rimecldhail, imean)
    end if


    if (associated(oneMicroFields%rain2ice)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2ice, &
            'RAIN2ICE :3:hist:anal:mpt3', &
            oneAveMicroFields%rain2ice, imean)
    end if

    if (associated(oneMicroFields%rain2pr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2pr, &
            'RAIN2PR :3:hist:anal:mpt3', &
            oneAveMicroFields%rain2pr, imean)
    end if

    if (associated(oneMicroFields%rain2sn)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2sn, &
            'RAIN2SN :3:hist:anal:mpt3', &
            oneAveMicroFields%rain2sn, imean)
    end if

    if (associated(oneMicroFields%rain2ag)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2ag, &
            'RAIN2AG :3:hist:anal:mpt3', &
            oneAveMicroFields%rain2ag, imean)
    end if

    if (associated(oneMicroFields%rain2gr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2gr, &
            'RAIN2GR :3:hist:anal:mpt3', &
            oneAveMicroFields%rain2gr, imean)
    end if

    if (associated(oneMicroFields%rain2ha)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2ha, &
            'RAIN2HA :3:hist:anal:mpt3', &
            oneAveMicroFields%rain2ha, imean)
    end if

    if (associated(oneMicroFields%rain2ha_xtra)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2ha_xtra, &
            'RAIN2HA_XTRA :3:hist:anal:mpt3', &
            oneAveMicroFields%rain2ha_xtra, imean)
    end if

    if (associated(oneMicroFields%ice2rain)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%ice2rain, &
            'ICE2RAIN :3:hist:anal:mpt3', &
            oneAveMicroFields%ice2rain, imean)
    end if

    if (associated(oneMicroFields%aggregate)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%aggregate, &
            'AGGREGATE :3:hist:anal:mpt3', &
            oneAveMicroFields%aggregate, imean)
    end if

    if (associated(oneMicroFields%aggrselfpris)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%aggrselfpris, &
            'AGGRSELFPRIS :3:hist:anal:mpt3', &
            oneAveMicroFields%aggrselfpris, imean)
    end if

    if (associated(oneMicroFields%aggrselfsnow)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%aggrselfsnow, &
            'AGGRSELFSNOW :3:hist:anal:mpt3', &
            oneAveMicroFields%aggrselfsnow, imean)
    end if

    if (associated(oneMicroFields%aggrprissnow)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%aggrprissnow, &
            'AGGRPRISSNOW :3:hist:anal:mpt3', &
            oneAveMicroFields%aggrprissnow, imean)
    end if


    if (associated(oneMicroFields%latheatvap)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%latheatvap, &
            'LATHEATVAP :3:hist:anal:mpt3', &
            oneAveMicroFields%latheatvap, imean)
    end if

    if (associated(oneMicroFields%latheatfrz)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%latheatfrz, &
            'LATHEATFRZ :3:hist:anal:mpt3', &
            oneAveMicroFields%latheatfrz, imean)
    end if

    !END MICRO BUDGET PROCESSES (instantaneous)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
    if (associated(oneMicroFields%nuccldrt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%nuccldrt, &
            'NUCCLDRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%nuccldrt, imean)
    end if

    if (associated(oneMicroFields%nuccldct)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%nuccldct, &
            'NUCCLDCT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%nuccldct, imean)
    end if

    if (associated(oneMicroFields%nucicert)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%nucicert, &
            'NUCICERT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%nucicert, imean)
    end if

    if (associated(oneMicroFields%nucicect)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%nucicect, &
            'NUCICECT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%nucicect, imean)
    end if

    if (associated(oneMicroFields%inuchomrt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuchomrt, &
            'INUCHOMRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%inuchomrt, imean)
    end if

    if (associated(oneMicroFields%inuchomct)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuchomct, &
            'INUCHOMCT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%inuchomct, imean)
    end if

    if (associated(oneMicroFields%inuccontrt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuccontrt, &
            'INUCCONTRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%inuccontrt, imean)
    end if

    if (associated(oneMicroFields%inuccontct)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuccontct, &
            'INUCCONTCT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%inuccontct, imean)
    end if

    if (associated(oneMicroFields%inucifnrt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inucifnrt, &
            'INUCIFNRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%inucifnrt, imean)
    end if

    if (associated(oneMicroFields%inucifnct)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inucifnct, &
            'INUCIFNCT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%inucifnct, imean)
    end if

    if (associated(oneMicroFields%inuchazrt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuchazrt, &
            'INUCHAZRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%inuchazrt, imean)
    end if

    if (associated(oneMicroFields%inuchazct)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%inuchazct, &
            'INUCHAZCT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%inuchazct, imean)
    end if


    if (associated(oneMicroFields%vapliqt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapliqt, &
            'VAPLIQT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vapliqt, imean)
    end if

    if (associated(oneMicroFields%vapicet)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapicet, &
            'VAPICET :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vapicet, imean)
    end if

    if (associated(oneMicroFields%vapcldt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapcldt, &
            'VAPCLDT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vapcldt, imean)
    end if

    if (associated(oneMicroFields%vapraint)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapraint, &
            'VAPRAINT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vapraint, imean)
    end if

    if (associated(oneMicroFields%vapprist)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapprist, &
            'VAPPRIST :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vapprist, imean)
    end if

    if (associated(oneMicroFields%vapsnowt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapsnowt, &
            'VAPSNOWT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vapsnowt, imean)
    end if

    if (associated(oneMicroFields%vapaggrt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapaggrt, &
            'VAPAGGRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vapaggrt, imean)
    end if

    if (associated(oneMicroFields%vapgraut)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapgraut, &
            'VAPGRAUT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vapgraut, imean)
    end if

    if (associated(oneMicroFields%vaphailt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vaphailt, &
            'VAPHAILT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vaphailt, imean)
    end if

    if (associated(oneMicroFields%vapdrizt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%vapdrizt, &
            'VAPDRIZT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%vapdrizt, imean)
    end if


    if (associated(oneMicroFields%melticet)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%melticet, &
            'MELTICET :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%melticet, imean)
    end if

    if (associated(oneMicroFields%meltprist)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%meltprist, &
            'MELTPRIST :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%meltprist, imean)
    end if

    if (associated(oneMicroFields%meltsnowt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%meltsnowt, &
            'MELTSNOWT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%meltsnowt, imean)
    end if

    if (associated(oneMicroFields%meltaggrt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%meltaggrt, &
            'MELTAGGRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%meltaggrt, imean)
    end if

    if (associated(oneMicroFields%meltgraut)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%meltgraut, &
            'MELTGRAUT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%meltgraut, imean)
    end if

    if (associated(oneMicroFields%melthailt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%melthailt, &
            'MELTHAILT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%melthailt, imean)
    end if


    if (associated(oneMicroFields%cld2raint)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%cld2raint, &
            'CLD2RAINT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%cld2raint, imean)
    end if

    if (associated(oneMicroFields%rimecldt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecldt, &
            'RIMECLDT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rimecldt, imean)
    end if

    if (associated(oneMicroFields%rimecldsnowt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecldsnowt, &
            'RIMECLDSNOWT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rimecldsnowt, imean)
    end if

    if (associated(oneMicroFields%rimecldaggrt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecldaggrt, &
            'RIMECLDAGGRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rimecldaggrt, imean)
    end if

    if (associated(oneMicroFields%rimecldgraut)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecldgraut, &
            'RIMECLDGRAUT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rimecldgraut, imean)
    end if

    if (associated(oneMicroFields%rimecldhailt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rimecldhailt, &
            'RIMECLDHAILT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rimecldhailt, imean)
    end if


    if (associated(oneMicroFields%rain2icet)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2icet, &
            'RAIN2ICET :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rain2icet, imean)
    end if

    if (associated(oneMicroFields%rain2prt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2prt, &
            'RAIN2PRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rain2prt, imean)
    end if

    if (associated(oneMicroFields%rain2snt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2snt, &
            'RAIN2SNT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rain2snt, imean)
    end if

    if (associated(oneMicroFields%rain2agt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2agt, &
            'RAIN2AGT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rain2agt, imean)
    end if

    if (associated(oneMicroFields%rain2grt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2grt, &
            'RAIN2GRT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rain2grt, imean)
    end if

    if (associated(oneMicroFields%rain2hat)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2hat, &
            'RAIN2HAT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rain2hat, imean)
    end if

    if (associated(oneMicroFields%rain2ha_xtrat)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2ha_xtrat, &
            'RAIN2HA_XTRAT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rain2ha_xtrat, imean)
    end if

    if (associated(oneMicroFields%ice2raint)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%ice2raint, &
            'ICE2RAINT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%ice2raint, imean)
    end if


    if (associated(oneMicroFields%aggregatet)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%aggregatet, &
            'AGGREGATET :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%aggregatet, imean)
    end if

    if (associated(oneMicroFields%aggrselfprist)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%aggrselfprist, &
            'AGGRSELFPRIST :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%aggrselfprist, imean)
    end if

    if (associated(oneMicroFields%aggrselfsnowt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%aggrselfsnowt, &
            'AGGRSELFSNOWT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%aggrselfsnowt, imean)
    end if

    if (associated(oneMicroFields%aggrprissnowt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%aggrprissnowt, &
            'AGGRPRISSNOWT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%aggrprissnowt, imean)
    end if


    if (associated(oneMicroFields%latheatvapt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%latheatvapt, &
            'LATHEATVAPT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%latheatvapt, imean)
    end if

    if (associated(oneMicroFields%latheatfrzt)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%latheatfrzt, &
            'LATHEATFRZT :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%latheatfrzt, imean)
    end if

    !END MICRO BUDGET PROCECCES (totals)

    if (associated(oneMicroFields%accpr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%accpr, &
            'ACCPR :2:hist:anal:mpti:mpt3', &
            oneAveMicroFields%accpr, imean)
    end if

    if (associated(oneMicroFields%accpp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%accpp, &
            'ACCPP :2:hist:anal:mpti:mpt3', &
            oneAveMicroFields%accpp, imean)
    end if

    if (associated(oneMicroFields%accps)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%accps, &
            'ACCPS :2:hist:anal:mpti:mpt3', &
            oneAveMicroFields%accps, imean)
    end if

    if (associated(oneMicroFields%accpa)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%accpa, &
            'ACCPA :2:hist:anal:mpti:mpt3', &
            oneAveMicroFields%accpa, imean)
    end if

    if (associated(oneMicroFields%accpg)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%accpg, &
            'ACCPG :2:hist:anal:mpti:mpt3', &
            oneAveMicroFields%accpg, imean)
    end if

    if (associated(oneMicroFields%accph)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%accph, &
            'ACCPH :2:hist:anal:mpti:mpt3', &
            oneAveMicroFields%accph, imean)
    end if

    if (associated(oneMicroFields%accpd)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%accpd, &
            'ACCPD :2:hist:anal:mpti:mpt3', &
            oneAveMicroFields%accpd, imean)
    end if

    if (associated(oneMicroFields%pcprr)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcprr, &
            'PCPRR :2:hist:anal:mpt3', &
            oneAveMicroFields%pcprr, imean)
    end if

    if (associated(oneMicroFields%pcprp)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcprp, &
            'PCPRP :2:hist:anal:mpt3', &
            oneAveMicroFields%pcprp, imean)
    end if

    if (associated(oneMicroFields%pcprs)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcprs, &
            'PCPRS :2:hist:anal:mpt3', &
            oneAveMicroFields%pcprs, imean)
    end if

    if (associated(oneMicroFields%pcpra)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcpra, &
            'PCPRA :2:hist:anal:mpt3', &
            oneAveMicroFields%pcpra, imean)
    end if

    if (associated(oneMicroFields%pcprg)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcprg, &
            'PCPRG :2:hist:anal:mpt3', &
            oneAveMicroFields%pcprg, imean)
    end if

    if (associated(oneMicroFields%pcprh)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcprh, &
            'PCPRH :2:hist:anal:mpt3', &
            oneAveMicroFields%pcprh, imean)
    end if

    if (associated(oneMicroFields%pcprd)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcprd, &
            'PCPRD :2:hist:anal:mpt3', &
            oneAveMicroFields%pcprd, imean)
    end if

    if (associated(oneMicroFields%pcpg)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%pcpg, &
            'PCPG :2:hist:mpti:mpt3', &
            oneAveMicroFields%pcpg, imean)
    end if

    if (associated(oneMicroFields%qpcpg)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%qpcpg, &
            'QPCPG :2:hist:mpti:mpt3', &
            oneAveMicroFields%qpcpg, imean)
    end if

    if (associated(oneMicroFields%dpcpg)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%dpcpg, &
            'DPCPG :2:hist:mpti:mpt3', &
            oneAveMicroFields%dpcpg, imean)
    end if

  end subroutine InsertMicroFieldsAtVarTable

end module ModMicroFields
