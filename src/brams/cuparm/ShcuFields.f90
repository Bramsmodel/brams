! Module necessary to Shallow Cumulus param.

module ModShcuFields

  use iso_fortran_env, only: &
       int64
  
  use ModNamelistFile, only: &
       NamelistFile

  use ModControlVars, only: &
       ControlVars

  use ModNodeDimensions, only: &
       NodeDimensions

  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModVarTables, only: &
       InsertVTab
  
  implicit none

  private
  public :: ShcuFields
  public :: CreateShcuFields
  public :: DestroyShcuFields
  public :: DumpShcuFields
  public :: InsertShcuFieldsAtVarTable 

  type ShcuFields
     real, pointer, contiguous :: thsrcsh(:,:,:) ! dimensioned (nzp,nxp,nyp)
     real, pointer, contiguous :: rtsrcsh(:,:,:) ! dimensioned (nzp,nxp,nyp)
     real, pointer, contiguous :: shmf(:,:) ! dimensioned (mxp,myp)
  end type ShcuFields

  
contains


  function CreateShcuFields(oneNodeDimensions, oneControlVars) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(ControlVars), pointer, intent(in) :: oneControlVars
    type(ShcuFields), pointer :: res

    integer :: mxp
    integer :: myp
    integer :: mzp
    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateShcuFields)**"

    if (.not. associated(oneNodeDimensions)) then
       call fatal_error(h//" oneNodeDimensions not associated")
    else if (.not. associated(oneControlVars)) then
       call fatal_error(h//" oneControlVars not associated")
    end if
       

    if (oneControlVars%nnshcu /= 1) then
       nullify(res)
       return
    end if
    
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp
    mzp=oneNodeDimensions%mzp

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if

    allocate(res%thsrcsh(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate thsrcsh fails with stat="// &
            trim(adjustl(str(1))))
    end if
    
    allocate(res%rtsrcsh(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate rtsrcsh fails with stat="// &
            trim(adjustl(str(1))))
    end if
    
    allocate(res%shmf(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate shmf fails with stat="// &
            trim(adjustl(str(1))))
    end if
  end function CreateShcuFields
    



  subroutine DestroyShcuFields(oneShcuFields)
    type(ShcuFields), pointer, intent(inout) :: oneShcuFields

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyShcuFields)**"
    
    if (.not. associated(oneShcuFields)) then
       return
    end if

    deallocate(oneShcuFields%thsrcsh, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate thsrcsh fails with stat="// &
            trim(adjustl(str(1))))
    end if
    
    deallocate(oneShcuFields%rtsrcsh, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate rtsrcsh fails with stat="// &
            trim(adjustl(str(1))))
    end if
    
    deallocate(oneShcuFields%shmf, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate shmf fails with stat="// &
            trim(adjustl(str(1))))
    end if

    deallocate(oneShcuFields, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneShcuFields fails with stat="//&
            trim(adjustl(str(1))))
    end if

    nullify(oneShcuFields)
  end subroutine DestroyShcuFields
    

    



  subroutine DumpShcuFields(oneShcuFields, name)
    type(ShcuFields), pointer, intent(in) :: oneShcuFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpShcuFields)**"
    
    if (associated(oneShcuFields)) then
       call MsgDump(h//" "//name//" associated")
    else
       call MsgDump(h//" "//name//" not associated")
    end if
  end subroutine DumpShcuFields




  subroutine InsertShcuFieldsAtVarTable(oneShcuFields, oneAveShcuFields, &
       oneControlVars, oneNamelistFile, gridId)
    type(ShcuFields), pointer, intent(in) :: oneShcuFields
    type(ShcuFields), pointer, intent(in) :: oneAveShcuFields
    type(ControlVars), pointer, intent(in) :: oneControlVars
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    integer, intent(in) :: gridId

    integer :: imean
    integer(kind=int64) :: npts
    character(len=*), parameter :: h="**(InsertShcuFieldsAtVarTable)**" 

    if (.not. associated(oneControlVars)) then
       call fatal_error(h//" oneControlVars not associated")
    end if
    
    if (oneControlVars%nnshcu /= 1) then
       return
    end if
    
    if (.not. associated(oneShcuFields)) then
       call fatal_error(h//" oneShcuFields not associated")
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

    
    if (associated(oneShcuFields%thsrcsh)) then
       npts=int(size(oneShcuFields%thsrcsh,1),int64) * &
            int(size(oneShcuFields%thsrcsh,2),int64) * &
            int(size(oneShcuFields%thsrcsh,3),int64)
       call InsertVTab (oneShcuFields%thsrcsh,oneAveShcuFields%thsrcsh &
            ,gridId, npts, imean,  &
            'THSRCSH :3:hist:anal:mpti:mpt3')
    end if
    
    if (associated(oneShcuFields%rtsrcsh)) then
       npts=int(size(oneShcuFields%rtsrcsh,1),int64) * &
            int(size(oneShcuFields%rtsrcsh,2),int64) * &
            int(size(oneShcuFields%rtsrcsh,3),int64)
       call InsertVTab (oneShcuFields%rtsrcsh,oneAveShcuFields%rtsrcsh &
            ,gridId, npts, imean,  &
            'RTSRCSH :3:hist:anal:mpti:mpt3')
    end if

    if (associated(oneShcuFields%shmf)) then
       npts=int(size(oneShcuFields%shmf,1),int64) * &
            int(size(oneShcuFields%shmf,2),int64) 
       call InsertVTab (oneShcuFields%shmf,oneAveShcuFields%shmf  &
            ,gridId, npts, imean,  &
            'SHMF :2:hist:anal:mpti:mpt3')
    end if
  end subroutine InsertShcuFieldsAtVarTable

end module ModShcuFields
