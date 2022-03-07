module ModFieldSectionList

  use ModParallelEnvironment, only: &
       MsgDump

  use ModFieldSection, only: &
       FieldSection, &
       FieldSectionName, &
       DumpFieldSection
  
  implicit none
  private

  public :: FieldSectionNode
  public :: FieldSectionList
  public :: CreateFieldSectionNode
  public :: CreateFieldSectionList
  public :: DumpFieldSectionList
  public :: AppendNodeToFieldSectionList
  public :: FindFieldNamed
  public :: FieldSectionListHeadNode
  public :: NextFieldSectionNodeAtList
  public :: FieldSectionListTailNode
  public :: FieldSectionAtNode

  interface FindFieldNamed
     module procedure FindFieldSectionAtList
  end interface FindFieldNamed

  type FieldSectionNode
     private
     type(FieldSection), pointer :: entry => null()
     ! field variable name
     type(FieldSectionNode), pointer :: next => null()
     type(FieldSectionNode), pointer :: previous => null()
     ! double linked list of FieldSection, since
     ! one communication may have multiple fields
     ! to communicate
  end type FieldSectionNode

  type FieldSectionList
     private
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
  
  

  
  function FindFieldSectionAtList(oneFieldSectionList, name) result(node)

    ! Finds the Field Section with name in the list;
    ! returns null Field Section if not in the list
    
    type(FieldSectionList), pointer, intent(in) :: oneFieldSectionList
    character(len=*), intent(in) :: name
    type(FieldSection), pointer :: node

    type(FieldSectionNode), pointer :: this
    type(FieldSection), pointer :: thisSection
    character(len=*), parameter :: h="**(FindFieldSectionAtList)**"
    logical, parameter :: dumpLocal=.false.

    node => null()
    if (associated(oneFieldSectionList)) then
       this => oneFieldSectionList%head
       do while (associated(this))
          thisSection => this%entry
          if (trim(adjustl(FieldSectionName(thisSection))) == &
               trim(adjustl(name))) then
             node => thisSection
             exit
          else
             this => this%next
          end if
       end do
    end if
  end function FindFieldSectionAtList


  
  function FieldSectionListHeadNode(oneFieldSectionList) result(this)
    type(FieldSectionList), pointer, intent(in) :: oneFieldSectionList
    type(FieldSectionNode), pointer :: this

    this => null()
    if (associated(oneFieldSectionList)) then
       this => oneFieldSectionList%head
    end if
  end function FieldSectionListHeadNode


  function NextFieldSectionNodeAtList(oneFieldSectionNode) result(this)
    type(FieldSectionNode), pointer, intent(in) :: oneFieldSectionNode
    type(FieldSectionNode), pointer :: this

    this => null()
    if (associated(oneFieldSectionNode)) then
       this => oneFieldSectionNode%next
    end if
  end function NextFieldSectionNodeAtList


  
  function FieldSectionListTailNode(oneFieldSectionList) result(this)
    type(FieldSectionList), pointer, intent(in) :: oneFieldSectionList
    type(FieldSectionNode), pointer :: this

    this => null()
    if (associated(oneFieldSectionList)) then
       this => oneFieldSectionList%tail
    end if
  end function FieldSectionListTailNode



  function FieldSectionAtNode(oneFieldSectionNode) result(this)
    type(FieldSectionNode), pointer, intent(in) :: oneFieldSectionNode
    type(FieldSection), pointer :: this

    this => null()
    if (associated(oneFieldSectionNode)) then
       this => oneFieldSectionNode%entry
    end if
  end function FieldSectionAtNode
end module ModFieldSectionList
