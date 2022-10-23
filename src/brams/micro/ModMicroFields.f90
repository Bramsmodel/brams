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
     real, contiguous, pointer :: vapliqt(:,:,:)
     real, contiguous, pointer :: vapicet(:,:,:)
     real, contiguous, pointer :: melticet(:,:,:)
     real, contiguous, pointer :: cld2raint(:,:,:)
     real, contiguous, pointer :: rimecldt(:,:,:)
     real, contiguous, pointer :: rain2icet(:,:,:)
     real, contiguous, pointer :: ice2raint(:,:,:)
     real, contiguous, pointer :: aggregatet(:,:,:)
     real, contiguous, pointer :: latheatvapt(:,:,:)
     real, contiguous, pointer :: latheatfrzt(:,:,:)

     ! Variables to be dimensioned by (nnxp,nyp)
     real, contiguous, pointer :: accpr(:,:)
     real, contiguous, pointer :: accpp(:,:)
     real, contiguous, pointer :: accps(:,:)
     real, contiguous, pointer :: accpa(:,:)
     real, contiguous, pointer :: accpg(:,:)
     real, contiguous, pointer :: accph(:,:)
     real, contiguous, pointer :: pcprr(:,:)
     real, contiguous, pointer :: pcprp(:,:)
     real, contiguous, pointer :: pcprs(:,:)
     real, contiguous, pointer :: pcpra(:,:)
     real, contiguous, pointer :: pcprg(:,:)
     real, contiguous, pointer :: pcprh(:,:)
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
    if (associated(oneMicroFields%melticet)) then
       deallocate(oneMicroFields%melticet, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate melticet fails with stat="//&
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
    if (associated(oneMicroFields%rain2icet)) then
       deallocate(oneMicroFields%rain2icet, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rain2icet fails with stat="//&
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

    if (associated(oneMicroFields%melticet)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%melticet, &
            'MELTICET :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%melticet, imean)
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

    if (associated(oneMicroFields%rain2icet)) then 
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneMicroFields%rain2icet, &
            'RAIN2ICET :3:hist:anal:mpti:mpt3', &
            oneAveMicroFields%rain2icet, imean)
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
