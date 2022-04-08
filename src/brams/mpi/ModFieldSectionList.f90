module ModFieldSectionList

  use ModParallelEnvironment, only: &
       MsgDump

  use ModFieldSection, only: &
       FieldSection, &
       DumpFieldSection

  implicit none
  private

  public :: FieldSectionNode
  public :: FieldSectionList
  public :: CreateFieldSectionNode
  public :: CreateFieldSectionList
  public :: DumpFieldSectionList
  public :: AppendNodeToFieldSectionList

  type FieldSectionNode
     type(FieldSection), pointer :: entry => null()
     ! field variable name
     type(FieldSectionNode), pointer :: next => null()
     type(FieldSectionNode), pointer :: previous => null()
     ! double linked list of FieldSection, since
     ! one communication may have multiple fields
     ! to communicate
  end type FieldSectionNode

  type FieldSectionList
     type(FieldSectionNode), pointer :: head => null()
     type(FieldSectionNode), pointer :: tail => null()
  end type FieldSectionList
contains



  function CreateFieldSectionList() result(new)
    type(FieldSectionList), pointer :: new

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSectionList)**"

    allocate(new, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate new fails with stat="//&
            trim(adjustl(c0)))
    end if
  end function CreateFieldSectionList





  function CreateFieldSectionNode(oneFieldSection) result(new)
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    type(FieldSectionNode), pointer :: new

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSectionNode)**"

    allocate(new, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate new fails with stat="//&
            trim(adjustl(c0)))
    end if

    new%entry => oneFieldSection
  end function CreateFieldSectionNode



  subroutine DumpFieldSectionList(oneFieldSectionList)
    type(FieldSectionList), pointer, intent(in) :: oneFieldSectionList

    type(FieldSectionNode), pointer :: this
    character(len=*), parameter :: h="**(DumpFieldSectionList)**"

    if (.not. associated(oneFieldSectionList)) then
       call MsgDump(h//" oneFieldSectionList not associated")
    else
       this => oneFieldSectionList%head
       if (.not. associated(this)) then
          call MsgDump(h//" empty oneFieldSectionList")
       else
          do while (associated(this))
             call DumpFieldSection(this%entry, h)
             this => this%next
          end do
       end if
    end if
  end subroutine DumpFieldSectionList




  subroutine AppendNodeToFieldSectionList(oneFieldSectionNode, oneFieldSectionList)

    type(FieldSectionNode), pointer, intent(in) :: oneFieldSectionNode
    type(FieldSectionList), pointer, intent(in) :: oneFieldSectionList

    integer :: ierr
    character(len=8) :: c0
    type(FieldSectionNode), pointer :: this
    character(len=*), parameter :: h="**(AppendNodeToFieldSectionList)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSectionNode)) then
       call fatal_error(h//" null oneFieldSectionNode")
    end if
    if (.not. associated(oneFieldSectionList)) then
       call fatal_error(h//" null oneFieldSectionList")
    end if

    if (associated(oneFieldSectionList%tail)) then
       oneFieldSectionNode%previous => oneFieldSectionList%tail
       oneFieldSectionList%tail%next => oneFieldSectionNode
       oneFieldSectionList%tail => oneFieldSectionNode
    else
       oneFieldSectionList%head => oneFieldSectionNode
       oneFieldSectionList%tail => oneFieldSectionNode
    end if
  end subroutine AppendNodeToFieldSectionList
end module ModFieldSectionList
