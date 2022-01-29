module ModGridDims

  ! ModGridDims: stores grid dimensions of a single grid
  !          at a variable of type "GridDims"

  use ModNamelistFile, only: NamelistFile
  use ModParallelEnvironment, only: ParallelEnvironment, MsgDump
  implicit none
  private
  public :: GridDims
  public :: CreateGridDims
  public :: DumpGridDims
  public :: DestroyGridDims

  type GridDims
     integer :: nnxp  ! # x points
     integer :: nnyp  ! # y points
     integer :: nnzp  ! # z points
  end type GridDims

contains



  ! CreateGridDims: create and fill variable of this type,
  !                 extracting info from the Namelist File.



  function CreateGridDims(gridId, oneNamelistFile)
    integer, intent(in) :: gridId
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(GridDims), pointer :: CreateGridDims
    
    character(len=16) :: c0, c1
    character(len=*), parameter :: h="**(CreateGridDims)**"
    logical, parameter :: dumpLocal=.false.

    ! correctness of input arguments

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneNamelistFile")
    else if (oneNamelistFile%ngrids <= 0) then
       write(c0,"(i16)") oneNamelistFile%ngrids 
       call fatal_error(h//" ngrids on namelist file "//&
            "should be >= 1 but is"//trim(adjustl(c0)))
    else if (gridId > oneNamelistFile%ngrids) then
       write(c0,"(i16)") oneNamelistFile%ngrids 
       write(c1,"(i16)") gridId 
       call fatal_error(h//" invoked with gridId ("//&
            trim(adjustl(c1))//") > number of grids on Namelist File ("//&
            trim(adjustl(c0))//")")
    else if (gridId < 1) then
       write(c1,"(i16)") gridId 
       call fatal_error(h//" invoked with gridId ("//&
            trim(adjustl(c1))//") < 1 ")
    end if

    ! create a variable of type GridDims and fill entries

    allocate(CreateGridDims)
    
    CreateGridDims%nnxp = oneNamelistFile%nnxp(gridId)
    CreateGridDims%nnyp = oneNamelistFile%nnyp(gridId)
    CreateGridDims%nnzp = oneNamelistFile%nnzp(gridId)

    if (dumpLocal) then
       call DumpGridDims(CreateGridDims, h)
    end if
  end function CreateGridDims



  ! DumpGridDims: prints info of a variable of type grid on Dump file



  subroutine DumpGridDims(oneGridDims, name)
    type(GridDims), pointer, intent(in) :: oneGridDims
    character(len=*), intent(in) :: name

    character(len=16) :: c0, c1, c2
    character(len=512) :: str
    character(len=*), parameter :: h="**(DumpGridDims)**"

    if (.not. associated(oneGridDims)) then
       call MsgDump(h//" null GridDims")
    else
       write(c0,"(i8)") oneGridDims%nnxp
       write(c1,"(i8)") oneGridDims%nnyp
       write(c2,"(i8)") oneGridDims%nnzp
       str = " "//name//" produces grid with dimensions "//&
            "["//trim(adjustl(c0))//","//trim(adjustl(c1))//&
            ","//trim(adjustl(c2))//"]"
       call MsgDump(h//trim(str))
    end if
  end subroutine DumpGridDims



  ! DestroyGridDims: deallocate area of a variable of type grid



  subroutine DestroyGridDims(oneGridDims)
    type(GridDims), pointer, intent(inout) :: oneGridDims
    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DestroyGridDims)**"
    
    if (associated(oneGridDims)) then
       deallocate(oneGridDims, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate fails with stat="//&
               trim(adjustl(c0)))
       end if
    end if
    nullify(oneGridDims)
  end subroutine DestroyGridDims
end module ModGridDims
