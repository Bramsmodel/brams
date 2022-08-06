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
  
  use ModVarTable, only: &
       VarTable, &
       InsertAtVarTable
  
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
    res%thsrcsh = 0.0
    
    allocate(res%rtsrcsh(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate rtsrcsh fails with stat="// &
            trim(adjustl(str(1))))
    end if
    res%rtsrcsh = 0.0
    
    allocate(res%shmf(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate shmf fails with stat="// &
            trim(adjustl(str(1))))
    end if
    res%shmf = 0.0
    
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




  subroutine InsertShcuFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneShcuFields, oneAveShcuFields)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(ShcuFields), pointer, intent(in) :: oneShcuFields
    type(ShcuFields), pointer, intent(in) :: oneAveShcuFields

    logical :: assAve
    logical :: assThis
    character(len=*), parameter :: h="**(InsertShcuFieldsAtVarTable)**" 

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(oneShcuFields)) then
       call fatal_error(h//" oneShcuFields not associated")
    end if
    

    ! Fill pointers to arrays into variable tables

    assAve=associated(oneAveShcuFields)
    
    if (associated(oneShcuFields%thsrcsh)) then
       if (assAve) then
          assThis=associated(oneAveShcuFields%thsrcsh)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneShcuFields%thsrcsh, &
               'THSRCSH :3:hist:anal:mpti:mpt3', &
               oneAveShcuFields%thsrcsh)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneShcuFields%thsrcsh, &
               'THSRCSH :3:hist:anal:mpti:mpt3')
       end if
    end if
    
    if (associated(oneShcuFields%rtsrcsh)) then
       if (assAve) then
          assThis=associated(oneAveShcuFields%rtsrcsh)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneShcuFields%rtsrcsh, &
               'RTSRCSH :3:hist:anal:mpti:mpt3', &
               oneAveShcuFields%rtsrcsh)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneShcuFields%rtsrcsh, &
               'RTSRCSH :3:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(oneShcuFields%shmf)) then
       if (assAve) then
          assThis=associated(oneAveShcuFields%shmf)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneShcuFields%shmf, &
               'SHMF :2:hist:anal:mpti:mpt3', &
               oneAveShcuFields%shmf)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneShcuFields%shmf, &
               'SHMF :2:hist:anal:mpti:mpt3')
       end if
    end if
  end subroutine InsertShcuFieldsAtVarTable

end module ModShcuFields
