module ModFieldSectionList


  ! field sections to be communicated
  ! among processes in a single message


  use ModParallelEnvironment, only: MsgDump


  implicit none
  include "ranks.h" ! for kind=i8

  private
  public :: FieldSection
  public :: CreateFieldSection
  public :: DumpFieldSection
  public :: FieldSectionList
  public :: CreateFieldSectionList
  public :: InsertAtFieldSectionList
  public :: DestroyFieldSectionList
  public :: DumpFieldSectionList
  public :: NextFieldSection


  ! FieldSection: one entry of a list of fields
  !               to be communicated to a single process
  !               in a single message passing operation.
  !               Data to communicate is the horizontal
  !               section [xStart:xEnd,yStart:yEnd] (local indices)
  !               of field_XXX (XXX=2D, 3D, 4D or I2D).
  !               If the field has more than 2 dimensions, then
  !               the remaining dimensions of each pair (x,y) of
  !               the section should be fully communicated.
  !               Component idim_type informs which are the remaining
  !               dimensions to be communicated, in a coded scheme.
  !               Component name has the field name to be communicated.
  !               Component fieldSectionSize is the size of the field
  !               to be communicated.


  type FieldSection
     real, pointer :: field_2D(:,:) => null()
     real, pointer :: field_3D(:,:,:) => null()
     real, pointer :: field_4D(:,:,:,:) => null()
     integer, pointer :: field_I2D(:,:) => null()
     ! field_XXX points to the array to extract
     ! the section to be communicated
     integer :: xStart = -1
     integer :: xEnd = -1
     integer :: yStart = -1
     integer :: yEnd = -1
     ! the 2D section to be communicated is, in local indices,
     ! [xStart:xEnd,yStart:yEnd]
     integer :: idim_type = -1
     ! field dimensioning code, to know which other dimensions
     ! should be communicated:
     ! idim_type == 2 means (nmxp, nmyp)
     !   no other dimensions to communicate
     ! idim_type == 3 means (nmzp, nmxp, nmyp)
     !   communicate first dimension for each (x,y)
     ! idim_type == 4 means (nzg, nmxp, nmyp, npatch)
     !   communicate first and last dimension for each (x,y)
     ! idim_type == 5 means (nzs, nmxp, nmyp, npatch)
     !   communicate first and last dimension for each (x,y)
     ! idim_type == 6 means (nmxp, nmyp, npatch)
     !   communicate  last dimension for each (x,y)
     ! idim_type == 7 means (nmxp, nmyp, nwave)
     !   communicate  last dimension for each (x,y)
     integer(kind=i8) :: fieldSectionSize = -1_i8
     ! number of data elements to communicate
     character(len=16) :: name = ""
     ! field variable name
     type(FieldSection), pointer :: next => null()
     type(FieldSection), pointer :: previous => null()
     ! double linked list of FieldSection, since
     ! one communication may have multiple fields
     ! to communicate
  end type FieldSection


  ! FieldSectionList: list of FieldSection


  type FieldSectionList
     type(FieldSection), pointer :: head=>null()
     type(FieldSection), pointer :: tail=>null()
  end type FieldSectionList


  interface CreateFieldSection
     module procedure CreateFieldSection_I2D
     module procedure CreateFieldSection_2D
     module procedure CreateFieldSection_3D
     module procedure CreateFieldSection_4D
  end interface CreateFieldSection
contains





  function CreateFieldSection_I2D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the integer 2D field should be communicated

    integer, pointer, intent(in) :: field(:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_I2D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(c0)))
    end if
    oneFieldSection%field_I2d => field
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type = idim_type
    if (idim_type == 2) then
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1)
    else
       write(c0,"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(c0)))
    end if
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneFieldSection, h)
    end if
  end function CreateFieldSection_I2D





  function CreateFieldSection_2D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 2D field should be communicated

    real, pointer, intent(in) :: field(:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_2D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(c0)))
    end if
    oneFieldSection%field_2d => field
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type = idim_type
    if (idim_type == 2) then
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1)
    else
       write(c0,"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(c0)))
    end if
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneFieldSection, h)
    end if
  end function CreateFieldSection_2D





  function CreateFieldSection_3D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 3D field should be communicated

    real, pointer, intent(in) :: field(:,:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_3D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(c0)))
    end if
    oneFieldSection%field_3d => field
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type=idim_type
    select case (idim_type)
    case (3)
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,1)
    case (6,7)
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,3)
    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(c0)))
    end select
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneFieldSection, h)
    end if
  end function CreateFieldSection_3D





  function CreateFieldSection_4D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 4D field should be communicated

    real, pointer, intent(in) :: field(:,:,:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_4D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(c0)))
    end if
    oneFieldSection%field_4d => field
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type=idim_type
    select case (idim_type)
    case (4, 5)
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,1) * &
            size(field,4)
    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(c0)))
    end select
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneFieldSection, h)
    end if
  end function CreateFieldSection_4D





  function StringFieldSection(oneFieldSection) result(res)

    ! String with the fields of a type FieldSection variable

    type(FieldSection), pointer :: oneFieldSection
    character(len=256) :: res

    character(len=128) :: string
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(StringFieldSection)**"

    if (.not. associated(oneFieldSection)) then
       res = " null FieldSection"
    else
       write(c0,"(i8)") oneFieldSection%xStart
       write(c1,"(i8)") oneFieldSection%xEnd
       write(c2,"(i8)") oneFieldSection%yStart
       write(c3,"(i8)") oneFieldSection%yEnd
       select case (oneFieldSection%idim_type)
       case(2)
          string="("//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
       case(3)
          write(c4,"(i8)") size(oneFieldSection%field_3D,1)
          string="(1:"//trim(adjustl(c4))//","//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
       case(4:5)
          write(c4,"(i8)") size(oneFieldSection%field_4D,1)
          write(c5,"(i8)") size(oneFieldSection%field_4D,4)
          string="(1:"//trim(adjustl(c4))//","//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               "1:"//trim(adjustl(c5))//")"
       case(6:7)
          write(c4,"(i8)") size(oneFieldSection%field_3D,3)
          string="("//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               "1:"//trim(adjustl(c4))//")"
       case default
          write(c0,"(i8)") oneFieldSection%idim_type
          call fatal_error(h//" field section "//trim(oneFieldSection%name)//&
               " with unknown idim_type="//trim(adjustl(c0)))
       end select
       write(c0,"(i8)") oneFieldSection%FieldSectionSize
       res = "field section "//trim(oneFieldSection%name)//&
            trim(string)//" of size "//trim(adjustl(c0))
    end if
  end function StringFieldSection





  subroutine DumpFieldSection(oneFieldSection, strMsg)

    ! Dumps a variable of type FieldSection with own header
    ! or with header "strMsg", if present

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    character(len=*), intent(in), optional :: strMsg

    character(len=*), parameter :: h="**(DumpFieldSection)**"

    if (present(strMsg)) then
       call MsgDump(trim(strMsg)//" "//&
            trim(adjustl(StringFieldSection(oneFieldSection))))
    else
       call MsgDump(h//" "//&
            trim(adjustl(StringFieldSection(oneFieldSection))))
    end if
  end subroutine DumpFieldSection





  subroutine DestroyFieldSection(oneEntry)

    ! reclaims memory area and returns null pointer

    type(FieldSection), pointer, intent(inout) :: oneEntry

    integer :: ierr
    character(len=128) :: name
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DestroyFieldSection)**"
    logical, parameter :: dumpLocal=.false.

    if (associated(oneEntry)) then
       name = oneEntry%name
       deallocate(oneEntry, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate fails with stat="//&
               trim(adjustl(c0)))
       end if
       if (dumpLocal) then
          call MsgDump(h//" named "//trim(adjustl(name)))
       end if
    end if
    nullify(oneEntry)
  end subroutine DestroyFieldSection


  ! CreateFieldSectionList: Creates empty list from a null pointer


  function CreateFieldSectionList() result (OneList)
    type (FieldSectionList), pointer :: OneList

    character(len=*), parameter :: h="**(CreateFieldSectionList)**"

    allocate(OneList)
    OneList%head=>null()
    OneList%tail=>null()
  end function CreateFieldSectionList


  ! InsertAtFieldSectionList: Given a FieldSectionList variable "list",
  !                           append FieldSection variable "node" to the
  !                           list


  subroutine InsertAtFieldSectionList(node, list)
    type(FieldSection), pointer :: node
    type(FieldSectionList), pointer :: list

    character(len=*), parameter :: h="**(InsertAtFieldSectionList)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(node)) then
       call fatal_error(h//" node not associated")
    else if (.not. associated(list)) then
       call fatal_error(h//" list not associated")
    end if

    if (.not. associated(list%head)) then

       ! if empty list, just point head and tail to the node

       list%head => node
       list%tail => node
    else if (.not. associated(list%tail)) then

       ! badly formed list if head is associated and tail is not associated

       call fatal_error(h//" bad list: head associated but null tail ")
    else

       ! append node to the end of the list

       list%tail%next => node
       node%previous => list%tail
       node%next => null()
       list%tail => node
    end if
    if (dumpLocal) then
       call MsgDump(h//" inserted the node "//&
            trim(adjustl(StringFieldSection(node))))
    end if
  end subroutine InsertAtFieldSectionList


  ! RemoveFromFieldSectionList: removes a node from the list,


  subroutine RemoveFromFieldSectionList(node, list)
    type(FieldSection), pointer :: node
    type(FieldSectionList), pointer :: list

    logical :: found
    type(FieldSection), pointer :: isThis
    character(len=*), parameter :: h="**(RemoveFromFieldSectionList)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(node)) then
       call fatal_error(h//" node not associated")
    else if (.not. associated(list)) then
       call fatal_error(h//" list not associated")
    end if

    found = .false.
    isThis => list%head
    do while (associated(isThis))
       found = associated(isThis, node)
       if (found) exit
       isThis => isThis%next
    end do
    if (.not. found) then
       call fatal_error(h//" node "//&
            trim(adjustl(StringFieldSection(node)))//&
            " is not in list")
    end if

    if (associated(list%head, node)) then
       list%head => node%next
       if (associated(list%head)) then
          list%head%previous => null()
       end if
    end if
    if (associated(list%tail, node)) then
       list%tail => node%previous
       if (associated(list%tail)) then
          list%tail%next => null()
       end if
    end if
    if (associated(node%previous)) then
       node%previous%next => node%next
    end if
    if (associated(node%next)) then
       node%next%previous => node%previous
    end if
    node%next => null()
    node%previous => null()
    if (dumpLocal) then
       call MsgDump(h//" removed the node "//&
            trim(adjustl(StringFieldSection(node))))
    end if
    call DestroyFieldSection(node)
  end subroutine RemoveFromFieldSectionList


  ! DestroyFieldSectionList: deallocates all nodes and the list


  subroutine DestroyFieldSectionList(list)
    type(FieldSectionList), pointer :: list

    type(FieldSection), pointer :: node
    character(len=*), parameter :: h="**(DestroyFieldSectionList)**"

    if (.not. associated(list)) then
       call MsgDump(h//" list not associated")
       return
    end if

    do
       if (.not. associated(list%head) .and. &
            .not. associated(list%tail)) then
          deallocate(list)
          list => null()
          exit
       else if (associated(list%head) .and. &
            associated(list%tail)) then
          node => list%head
          call RemoveFromFieldSectionList(node, list)
       else if (associated(list%head)) then
          call fatal_error(h//&
               " head associated but tail not associated")
       else if (associated(list%tail)) then
          call fatal_error(h//&
               " head not associated but tail associated")
       end if
    end do
  end subroutine DestroyFieldSectionList


  ! DumpFieldSectionList: Dumps a variable of type FieldSectionList at 
  !                       this processor dump file


  subroutine DumpFieldSectionList(list)
    type(FieldSectionList), pointer :: list

    integer :: cnt
    character(len=8) :: c0
    type(FieldSection), pointer :: node
    character(len=*), parameter :: h="**(DumpFieldSectionList)**"

    if (.not. associated(list)) then
       call MsgDump(h//" list is not associated")
    else

       ! list length

       cnt = 0
       node=>list%head
       do
          if (.not. associated(node)) exit
          cnt = cnt+1
          node => node%next
       end do
       if (cnt == 0) then
          call MsgDump(h//" empty list")
       else
          node=>list%head
          do
             if (.not. associated(node)) then
                exit
             end if
             call MsgDump(h//" "//trim(adjustl(StringFieldSection(node))))
             node => node%next
          end do
       end if
    end if
  end subroutine DumpFieldSectionList


  ! NextFieldSection: returns node following "node" at the list;
  !                   if input "node" is empty, returns list head;
  !                   if no more nodes in the list, returns null


  function NextFieldSection(node, list) result(next)
    type(FieldSection), pointer :: node
    type(FieldSectionList), pointer :: list
    type(FieldSection), pointer :: next

    character(len=*), parameter :: h="**(NextFieldSection)**"
    logical, parameter :: dumpLocal=.false.

    next => null()
    if (.not. associated(list)) then
       call fatal_error(h//" invoked with not associated list")
    else
       if (associated(list%head)) then
          if (.not. associated(node)) then
             next => list%head
          else
             next => node%next
          end if
       end if
    end if

    if (dumpLocal) then
       call MsgDump(h//" returns "//&
            trim(adjustl(StringFieldSection(next))))
    end if
  end function NextFieldSection
end module ModFieldSectionList
