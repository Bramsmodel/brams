!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModTurbFields

  use iso_fortran_env, only: &
       int64
  
  use ModNamelistFile, only: &
       NamelistFile

  use ModNodeDimensions, only: &
       NodeDimensions
  
  use ModParallelEnvironment, only: &
       MsgDump

  use var_tables, only: &
       InsertVTab

  implicit none

  private
  public :: TurbFields
  public :: CreateTurbFields
  public :: DestroyTurbFields
  public :: DumpTurbFields
  public :: InsertTurbFieldsAtVarTable

  character(len=32) :: previousCall=""
  
  type TurbFields

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, contiguous, pointer :: tkep(:,:,:) => null()
     real, contiguous, pointer :: epsp(:,:,:) => null()
     real, contiguous, pointer :: hkm(:,:,:) => null()
     real, contiguous, pointer :: vkm(:,:,:) => null()
     real, contiguous, pointer :: vkh(:,:,:) => null()
     real, contiguous, pointer :: cdrag(:,:,:) => null()

     ! Variables to be dimensioned by (nxp,nyp)
     real, contiguous, pointer :: sflux_u(:,:) => null()
     real, contiguous, pointer :: sflux_v(:,:) => null()
     real, contiguous, pointer :: sflux_w(:,:) => null()
     real, contiguous, pointer :: sflux_t(:,:) => null()
     real, contiguous, pointer :: sflux_r(:,:) => null()

     ![MLO - For Nakanishi/Niino
     real, contiguous, pointer :: kpbl(:,:) => null()
     !MLO]

  end type TurbFields

contains




  function CreateTurbFields(oneNodeDims, oneNamelistFile, gridId, createAverage) result(res)
    ! createAverage iff creating average fields
    type(NodeDimensions), pointer, intent(in) :: oneNodeDims
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    integer, intent(in) :: gridId
    logical, intent(in) :: createAverage
    type(TurbFields), pointer :: res

    integer :: ierr
    integer :: mzp
    integer :: mxp
    integer :: myp
    integer :: idiffk
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateTurbFields)**"
    logical, parameter :: dumpLocal=.true.

    mzp=oneNodeDims%mzp
    mxp=oneNodeDims%mxp
    myp=oneNodeDims%myp

    ! whenever creating average fields, allocate full size fields
    ! only when average fields will be output
    if (createAverage .and. oneNamelistFile%avgtim == 0.0) then
       mzp=1
       mxp=1
       myp=1
    end if
    
    if (dumpLocal) then
       write(str(1),"(i8)") mzp
       write(str(2),"(i8)") mxp
       write(str(3),"(i8)") myp
       call MsgDump(h//" starts with"//&
            " mzp="//trim(adjustl(str(1)))//&
            ", mxp="//trim(adjustl(str(2)))//&
            ", myp="//trim(adjustl(str(3))))
    end if

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if


    idiffk=oneNamelistFile%idiffk(gridId)

    if (idiffk == 1 .or. idiffk == 4 .or.  &
         idiffk == 5 .or. idiffk == 6 .or.  &
         idiffk == 7) then
       allocate (res%tkep(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate tkep fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%tkep=0.0
    end if

    if (idiffk == 5 .or. idiffk == 6) then
       allocate (res%epsp(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate epsp fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%epsp=0.0
    end if

    allocate (res%hkm(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate hkm fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%hkm=0.0


    allocate (res%vkm(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate vkm fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%vkm=0.0

    allocate (res%vkh(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate vkh fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%vkh=0.0

    if (oneNamelistFile%if_urban_canopy > 0) then
       allocate (res%cdrag(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate cdrag fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%cdrag=0.0
    end if

    allocate (res%sflux_u(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate sflux_u fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%sflux_u=0.0

    allocate (res%sflux_v(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate sflux_v fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%sflux_v=0.0

    allocate (res%sflux_w(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate sflux_w fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%sflux_w=0.0

    allocate (res%sflux_t(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate sflux_t fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%sflux_t=0.0

    allocate (res%sflux_r(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate sflux_r fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%sflux_r=0.0

    if (idiffk == 7) then
       allocate (res%kpbl(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate kpbl fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%kpbl=0.0
    end if

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end function CreateTurbFields





  subroutine DestroyTurbFields(oneTurbFields)
    type(TurbFields), pointer, intent(inout) :: oneTurbFields

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyTurbFields)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(oneTurbFields)) then
       return
    end if

    if (associated(oneTurbFields%tkep)) then
       deallocate(oneTurbFields%tkep, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate tkep fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%epsp)) then
       deallocate(oneTurbFields%epsp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate epsp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%hkm)) then
       deallocate(oneTurbFields%hkm, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate hkm fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%vkm)) then
       deallocate(oneTurbFields%vkm, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vkm fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%vkh)) then
       deallocate(oneTurbFields%vkh, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vkh fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%cdrag)) then
       deallocate(oneTurbFields%cdrag, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cdrag fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%sflux_u)) then
       deallocate(oneTurbFields%sflux_u, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate sflux_u fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%sflux_v)) then
       deallocate(oneTurbFields%sflux_v, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate sflux_v fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%sflux_w)) then
       deallocate(oneTurbFields%sflux_w, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate sflux_w fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%sflux_t)) then
       deallocate(oneTurbFields%sflux_t, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate sflux_t fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%sflux_r)) then
       deallocate(oneTurbFields%sflux_r, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate sflux_r fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    if (associated(oneTurbFields%kpbl)) then
       deallocate(oneTurbFields%kpbl, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate kpbl fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    deallocate(oneTurbFields, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneTurbFields fails with stat="//&
            trim(adjustl(str(1))))
    end if
    nullify(oneTurbFields)

  end subroutine DestroyTurbFields





  subroutine DumpTurbFields(oneTurbFields, name)
    type(TurbFields), pointer, intent(inout) :: oneTurbFields
    character(len=*), intent(in) :: name

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpTurbFields)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(oneTurbFields)) then
       call MsgDump(h//" "//name//" not associated")
       return
    end if

    if (associated(oneTurbFields%tkep)) then
       write(str(1),"(i8)") size(oneTurbFields%tkep,1)
       write(str(2),"(i8)") size(oneTurbFields%tkep,2)
       write(str(3),"(i8)") size(oneTurbFields%tkep,3)
       call MsgDump(h//" "//name//"%tkep associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    else
       call MsgDump(h//" "//name//"%tkep not associated")
    end if

    if (associated(oneTurbFields%epsp)) then
       write(str(1),"(i8)") size(oneTurbFields%epsp,1)
       write(str(2),"(i8)") size(oneTurbFields%epsp,2)
       write(str(3),"(i8)") size(oneTurbFields%epsp,3)
       call MsgDump(h//" "//name//"%epsp associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    else
       call MsgDump(h//" "//name//"%epsp not associated")
    end if

    if (associated(oneTurbFields%hkm)) then
       write(str(1),"(i8)") size(oneTurbFields%hkm,1)
       write(str(2),"(i8)") size(oneTurbFields%hkm,2)
       write(str(3),"(i8)") size(oneTurbFields%hkm,3)
       call MsgDump(h//" "//name//"%hkm associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    else
       call MsgDump(h//" "//name//"%hkm not associated")
    end if

    if (associated(oneTurbFields%vkm)) then
       write(str(1),"(i8)") size(oneTurbFields%vkm,1)
       write(str(2),"(i8)") size(oneTurbFields%vkm,2)
       write(str(3),"(i8)") size(oneTurbFields%vkm,3)
       call MsgDump(h//" "//name//"%vkm associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    else
       call MsgDump(h//" "//name//"%vkm not associated")
    end if

    if (associated(oneTurbFields%vkh)) then
       write(str(1),"(i8)") size(oneTurbFields%vkh,1)
       write(str(2),"(i8)") size(oneTurbFields%vkh,2)
       write(str(3),"(i8)") size(oneTurbFields%vkh,3)
       call MsgDump(h//" "//name//"%vkh associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    else
       call MsgDump(h//" "//name//"%vkh not associated")
    end if

    if (associated(oneTurbFields%cdrag)) then
       write(str(1),"(i8)") size(oneTurbFields%cdrag,1)
       write(str(2),"(i8)") size(oneTurbFields%cdrag,2)
       write(str(3),"(i8)") size(oneTurbFields%cdrag,3)
       call MsgDump(h//" "//name//"%cdrag associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    else
       call MsgDump(h//" "//name//"%cdrag not associated")
    end if

    if (associated(oneTurbFields%sflux_u)) then
       write(str(1),"(i8)") size(oneTurbFields%sflux_u,1)
       write(str(2),"(i8)") size(oneTurbFields%sflux_u,2)
       call MsgDump(h//" "//name//"%sflux_u associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//")")
    else
       call MsgDump(h//" "//name//"%sflux_u not associated")
    end if

    if (associated(oneTurbFields%sflux_v)) then
       write(str(1),"(i8)") size(oneTurbFields%sflux_v,1)
       write(str(2),"(i8)") size(oneTurbFields%sflux_v,2)
       call MsgDump(h//" "//name//"%sflux_v associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//")")
    else
       call MsgDump(h//" "//name//"%sflux_v not associated")
    end if

    if (associated(oneTurbFields%sflux_w)) then
       write(str(1),"(i8)") size(oneTurbFields%sflux_w,1)
       write(str(2),"(i8)") size(oneTurbFields%sflux_w,2)
       call MsgDump(h//" "//name//"%sflux_w associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//")")
    else
       call MsgDump(h//" "//name//"%sflux_w not associated")
    end if

    if (associated(oneTurbFields%sflux_t)) then
       write(str(1),"(i8)") size(oneTurbFields%sflux_t,1)
       write(str(2),"(i8)") size(oneTurbFields%sflux_t,2)
       call MsgDump(h//" "//name//"%sflux_t associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//")")
    else
       call MsgDump(h//" "//name//"%sflux_t not associated")
    end if

    if (associated(oneTurbFields%sflux_r)) then
       write(str(1),"(i8)") size(oneTurbFields%sflux_r,1)
       write(str(2),"(i8)") size(oneTurbFields%sflux_r,2)
       call MsgDump(h//" "//name//"%sflux_r associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//")")
    else
       call MsgDump(h//" "//name//"%sflux_r not associated")
    end if

    if (associated(oneTurbFields%kpbl)) then
       write(str(1),"(i8)") size(oneTurbFields%kpbl,1)
       write(str(2),"(i8)") size(oneTurbFields%kpbl,2)
       call MsgDump(h//" "//name//"%kpbl associated with size ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//")")
    else
       call MsgDump(h//" "//name//"%kpbl not associated")
    end if

  end subroutine DumpTurbFields




  subroutine InsertTurbFieldsAtVarTable(oneTurbFields, oneAveTurbFields, &
       oneNamelistFile, gridId)
    type(TurbFields), pointer, intent(in) :: oneTurbFields
    type(TurbFields), pointer, intent(in) :: oneAveTurbFields
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    integer, intent(in) :: gridId

    ! Local Variables:
    integer :: imean
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertTurbFieldsAtVarTable)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(oneTurbFields)) then
       call fatal_error(h//" oneTurbFields not associated")
    else if (.not. associated(oneAveTurbFields)) then
       call fatal_error(h//" oneAveTurbFields not associated")
    else if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if

    ! Should average fields be stored at variable tables?

    if (oneNamelistFile%avgtim == 0) then
       imean=0 ! do not store
    else
       imean=1 ! store
    end if

    ! Fill pointers to arrays into variable tables

    if (associated(oneTurbFields%tkep)) then
       call InsertVTab(oneTurbFields%tkep, oneAveTurbFields%tkep,  &
            gridId, int(size(oneTurbFields%tkep),int64), imean, 'TKEP :3:hist:anal:mpti:mpt3:mpt1')
    end if

    if (associated(oneTurbFields%epsp)) then
       call InsertVTab(oneTurbFields%epsp, oneAveTurbFields%epsp,  &
            gridId, int(size(oneTurbFields%epsp),int64), imean, 'EPSP :3:hist:anal:mpti:mpt3:mpt1')
    end if

    if (associated(oneTurbFields%hkm)) then
       call InsertVTab(oneTurbFields%hkm, oneAveTurbFields%hkm,  &
            gridId, int(size(oneTurbFields%hkm),int64), imean, 'HKM :3:hist:anal:mpti:mpt3:mpt1')
    end if

    if (associated(oneTurbFields%vkm)) then
       call InsertVTab(oneTurbFields%vkm, oneAveTurbFields%vkm,  &
            gridId, int(size(oneTurbFields%vkm),int64), imean, 'VKM :3:hist:mpti:mpt3:mpt1')
    end if

    if (associated(oneTurbFields%vkh)) then
       call InsertVTab(oneTurbFields%vkh, oneAveTurbFields%vkh,  &
            gridId, int(size(oneTurbFields%vkh),int64), imean, 'VKH :3:hist:anal:mpti:mpt3:mpt1')
    end if

    if (associated(oneTurbFields%cdrag)) then
       call InsertVTab(oneTurbFields%cdrag, oneAveTurbFields%cdrag,  &
            gridId, int(size(oneTurbFields%cdrag),int64), imean, 'CDRAG :3:hist:anal:mpti')
    end if

    if (associated(oneTurbFields%sflux_u)) then
       call InsertVTab(oneTurbFields%sflux_u, oneAveTurbFields%sflux_u,  &
            gridId, int(size(oneTurbFields%sflux_u),int64), imean, 'SFLUX_U :2:anal:mpt3:mpt1')
    end if

    if (associated(oneTurbFields%sflux_v)) then
       call InsertVTab(oneTurbFields%sflux_v, oneAveTurbFields%sflux_v,  &
            gridId, int(size(oneTurbFields%sflux_v),int64), imean, 'SFLUX_V :2:anal:mpt3:mpt1')
    end if

    if (associated(oneTurbFields%sflux_w)) then
       call InsertVTab (oneTurbFields%sflux_w, oneAveTurbFields%sflux_w,  &
            gridId, int(size(oneTurbFields%sflux_w),int64), imean, 'SFLUX_W :2:anal:mpt3')
    end if

    if (associated(oneTurbFields%sflux_t)) then
       call InsertVTab (oneTurbFields%sflux_t, oneAveTurbFields%sflux_t,  &
            gridId, int(size(oneTurbFields%sflux_t),int64), imean, 'SFLUX_T :2:anal:mpt3')
    end if

    if (associated(oneTurbFields%sflux_r)) then
       call InsertVTab (oneTurbFields%sflux_r, oneAveTurbFields%sflux_r,  &
            gridId, int(size(oneTurbFields%sflux_r),int64), imean, 'SFLUX_R :2:anal:mpt3')
    end if

    if (associated(oneTurbFields%kpbl)) then
       call InsertVTab (oneTurbFields%kpbl, oneAveTurbFields%kpbl,  &
            gridId, int(size(oneTurbFields%kpbl),int64), imean, 'KPBL :2:anal:mpt3')
    end if

  end subroutine InsertTurbFieldsAtVarTable
end module ModTurbFields
