module ModFieldSectionList


  ! List of field sections to be communicated
  ! among processes in a single message


  use var_tables, only: &
       var_tables_r, &
       GetVTabSectionSize, &
       VerifyVTabEntry

  use ModParallelEnvironment, only: &
       MsgDump


  use ModDomainDecomp, only: &
       DomainDecomp

  implicit none

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


  ! FieldSection: one entry of a list of field sections to be sent/received 
  !               in a single message.
  !               It represents the section [xStart:xEnd,yStart:yEnd]
  !               of the field pointed by field_XXX (XXX=2D, 3D, 4D or I2D).
  !               Number of points (reals) to be communicated is fieldSectionSize.

  integer, parameter :: UNDEF=-1
  type FieldSection
     real, pointer :: field_2D(:,:) => null()
     real, pointer :: field_3D(:,:,:) => null()
     real, pointer :: field_4D(:,:,:,:) => null()
     integer, pointer :: field_I2D(:,:) => null()
     integer :: xStart=UNDEF   ! local index
     integer :: xEnd=UNDEF     ! local index
     integer :: yStart=UNDEF   ! local index
     integer :: yEnd=UNDEF     ! local index
     integer :: fieldSectionSize=UNDEF            ! # reals to communicate
     integer :: idim_type=UNDEF
     character(len=16) :: name=""
     type(FieldSection), pointer :: next=>null()
     type(FieldSection), pointer :: previous=>null()
  end type FieldSection


  ! FieldSectionList: list of FieldSection


  type FieldSectionList
     type(FieldSection), pointer :: head=>null()
     type(FieldSection), pointer :: tail=>null()
  end type FieldSectionList


  interface CreateFieldSection
     module procedure CreateFieldSectionVarTable
     module procedure CreateFieldSection_I2D
     module procedure CreateFieldSection_2D
     module procedure CreateFieldSection_3D
     module procedure CreateFieldSection_4D
  end interface CreateFieldSection
contains


  ! CreateFieldSectionVarTable: returns pointer to a newly created 
  !                     FieldSection variable


  function CreateFieldSectionVarTable(vTabPtr, &
       xStart, xEnd, yStart, yEnd, GlobalWithGhost) result(oneEntry)
    type(var_tables_r), pointer :: vTabPtr
    integer, intent(in) :: xStart
    integer, intent(in) :: xEnd
    integer, intent(in) :: yStart
    integer, intent(in) :: yEnd
    type(DomainDecomp), pointer :: GlobalWithGhost
    type(FieldSection), pointer :: oneEntry

    character(len=*), parameter :: h="**(CreateFieldSectionVarTable)**"
    logical, parameter :: dumpLocal=.false.
    character(len=8) :: c0

    ! verify input arguments

    call VerifyVTabEntry(vTabPtr)

    ! allocate and fill entry

    allocate(oneEntry)
    select case (vTabPtr%idim_type)
    case (2)
       oneEntry%field_2D => vTabPtr%var_p_2D
    case (3)
       oneEntry%field_3D => vTabPtr%var_p_3D
    case (4)
       oneEntry%field_4D => vTabPtr%var_p_4D
    case (5)
       oneEntry%field_4D => vTabPtr%var_p_4D
    case (6)
       oneEntry%field_3D => vTabPtr%var_p_3D
    case (7)
       oneEntry%field_3D => vTabPtr%var_p_3D
    case default
       write(c0,"(i8)") vTabPtr%idim_type
       call fatal_error(h//" unknown idim_type="//trim(adjustl(c0)))
    end select
!!$    oneEntry%vTabPtr => vTabPtr
    oneEntry%xStart = xStart
    oneEntry%xEnd = xEnd
    oneEntry%yStart = yStart
    oneEntry%yEnd = yEnd
    oneEntry%fieldSectionSize = GetVTabSectionSize(vTabPtr, &
         xStart, xEnd, yStart, yEnd)
    oneEntry%idim_type=vTabPtr%idim_type
    oneEntry%name=vTabPtr%name
    oneEntry%next=>null()
    oneEntry%previous=>null()
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneEntry)
    end if
  end function CreateFieldSectionVarTable



  function CreateFieldSection_I2D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd, &
       GlobalWithGhost) result(oneEntry)
    integer, pointer, intent(in) :: field(:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart
    integer, intent(in) :: xEnd
    integer, intent(in) :: yStart
    integer, intent(in) :: yEnd
    type(DomainDecomp), pointer :: GlobalWithGhost
    type(FieldSection), pointer :: oneEntry

    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_I2D)**"
    logical, parameter :: dumpLocal=.false.

    ! allocate and fill entry

    allocate(oneEntry)
    oneEntry%field_I2d => field
    oneEntry%xStart = xStart
    oneEntry%xEnd = xEnd
    oneEntry%yStart = yStart
    oneEntry%yEnd = yEnd
    if (idim_type == 2) then
       oneEntry%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1)
    else
       write(c0,"(i8)") idim_type
       call fatal_error(h//" unknown idim_type="//&
            trim(adjustl(c0)))
    end if
    oneEntry%fieldSectionSize = &
         (yEnd - yStart +1) * &
         (xEnd - xStart +1) 
    oneEntry%idim_type=idim_type
    oneEntry%name=name
    oneEntry%next=>null()
    oneEntry%previous=>null()
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneEntry)
    end if
  end function CreateFieldSection_I2D



  function CreateFieldSection_2D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd, &
       GlobalWithGhost) result(oneEntry)
    real, pointer, intent(in) :: field(:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart
    integer, intent(in) :: xEnd
    integer, intent(in) :: yStart
    integer, intent(in) :: yEnd
    type(DomainDecomp), pointer :: GlobalWithGhost
    type(FieldSection), pointer :: oneEntry

    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_2D)**"
    logical, parameter :: dumpLocal=.false.

    ! allocate and fill entry

    allocate(oneEntry)
    oneEntry%field_2d => field
    oneEntry%xStart = xStart
    oneEntry%xEnd = xEnd
    oneEntry%yStart = yStart
    oneEntry%yEnd = yEnd
    if (idim_type == 2) then
       oneEntry%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1)
    else
       write(c0,"(i8)") idim_type
       call fatal_error(h//" unknown idim_type="//&
            trim(adjustl(c0)))
    end if
    oneEntry%fieldSectionSize = &
         (yEnd - yStart +1) * &
         (xEnd - xStart +1) 
    oneEntry%idim_type=idim_type
    oneEntry%name=name
    oneEntry%next=>null()
    oneEntry%previous=>null()
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneEntry)
    end if
  end function CreateFieldSection_2D
  



  function CreateFieldSection_3D(field, idim_type, name, &
       xStart, xEnd, yStart, yEnd, &
       GlobalWithGhost) result(oneEntry)
    real, pointer, intent(in) :: field(:,:,:)
    integer, intent(in) :: idim_type
    character(len=*), intent(in) :: name
    integer, intent(in) :: xStart
    integer, intent(in) :: xEnd
    integer, intent(in) :: yStart
    integer, intent(in) :: yEnd
    type(DomainDecomp), pointer :: GlobalWithGhost
    type(FieldSection), pointer :: oneEntry

    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_3D)**"
    logical, parameter :: dumpLocal=.false.

    ! allocate and fill entry

    allocate(oneEntry)
    oneEntry%field_3d => field
    oneEntry%xStart = xStart
    oneEntry%xEnd = xEnd
    oneEntry%yStart = yStart
    oneEntry%yEnd = yEnd
    if (idim_type == 3) then
       oneEntry%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,1)
    else if (idim_type == 6) then
       oneEntry%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,1)
    else if (idim_type == 7) then
       oneEntry%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,3)
    else
       write(c0,"(i8)") idim_type
       call fatal_error(h//" unknown idim_type="//&
            trim(adjustl(c0)))
    end if
    oneEntry%idim_type=idim_type
    oneEntry%name=name
    oneEntry%next=>null()
    oneEntry%previous=>null()
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneEntry)
    end if
  end function CreateFieldSection_3D
  



  function CreateFieldSection_4D(field, idim_type, name, &
       xStart, xEnd, yStart, yEnd, &
       GlobalWithGhost) result(oneEntry)
    real, pointer, intent(in) :: field(:,:,:,:)
    integer, intent(in) :: idim_type
    character(len=*), intent(in) :: name
    integer, intent(in) :: xStart
    integer, intent(in) :: xEnd
    integer, intent(in) :: yStart
    integer, intent(in) :: yEnd
    type(DomainDecomp), pointer :: GlobalWithGhost
    type(FieldSection), pointer :: oneEntry

    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_4D)**"
    logical, parameter :: dumpLocal=.false.

    ! allocate and fill entry

    allocate(oneEntry)
    oneEntry%field_4d => field
    oneEntry%xStart = xStart
    oneEntry%xEnd = xEnd
    oneEntry%yStart = yStart
    oneEntry%yEnd = yEnd
    if (idim_type == 4) then
       oneEntry%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,1) *&
            size(field,4)
    else if (idim_type == 5) then
       oneEntry%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,1) *&
            size(field,4)
    else
       write(c0,"(i8)") idim_type
       call fatal_error(h//" unknown idim_type="//&
            trim(adjustl(c0)))
    end if
    oneEntry%idim_type=idim_type
    oneEntry%name=name
    oneEntry%next=>null()
    oneEntry%previous=>null()
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneEntry)
    end if
  end function CreateFieldSection_4D
  


  ! StringFieldSection: Returns a string with the fields of 
  !                     a variable of type FieldSection


  function StringFieldSection(oneEntry) result(res)
    type(FieldSection), pointer :: oneEntry
    character(len=256) :: res

    character(len=128) :: string
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(StringFieldSection)**"

    if (.not. associated(oneEntry)) then
       res = " null FieldSection"
    else if (oneEntry%idim_type == UNDEF) then
       call fatal_error(h//" undefined idim_type: FieldSection not created")
    else
       write(c0,"(i8)") oneEntry%xStart
       write(c1,"(i8)") oneEntry%xEnd
       write(c2,"(i8)") oneEntry%yStart
       write(c3,"(i8)") oneEntry%yEnd
       select case (oneEntry%idim_type)
       case(2)
          string="("//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
       case(3)
          write(c4,"(i8)") size(oneEntry%field_3D,1)
          string="(1:"//trim(adjustl(c4))//","//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
       case(4:5)
          write(c4,"(i8)") size(oneEntry%field_4D,1)
          write(c5,"(i8)") size(oneEntry%field_4D,4)
          string="(1:"//trim(adjustl(c4))//","//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               "1:"//trim(adjustl(c5))//")"
       case(6:7)
          write(c4,"(i8)") size(oneEntry%field_3D,3)
          string="("//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               "1:"//trim(adjustl(c4))//")"
       case default
          write(c0,"(i8)") oneEntry%idim_type
          call fatal_error(h//" field section "//trim(oneEntry%name)//&
               " with unknown idim_type="//trim(adjustl(c0)))
       end select
       write(c0,"(i8)") oneEntry%FieldSectionSize
       res = "field section (local indices) "//trim(oneEntry%name)//&
            trim(string)//" of size "//trim(adjustl(c0))
    end if
  end function StringFieldSection

  
  ! DumpFieldSection: Dumps a variable of type FieldSection at 
  !                   this processor dump file


  subroutine DumpFieldSection(oneEntry)
    type(FieldSection), pointer :: oneEntry

    character(len=*), parameter :: h="**(DumpFieldSection)**"

    call MsgDump(h//trim(adjustl(StringFieldSection(oneEntry))))
  end subroutine DumpFieldSection
            

  ! DestroyFieldSection: reclaims memory area and returns
  !                      null pointer


  subroutine DestroyFieldSection(oneEntry)
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
