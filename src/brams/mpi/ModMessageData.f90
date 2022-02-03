module ModMessageData
  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModNeighbourNodes, only: &
       NeighbourNodes

  use ModDomainDecomp, only: &
       DomainDecomp

  use ModFieldSectionList, only: &
       FieldSection, &
       CreateFieldSection, &
       FieldSectionList, &
       CreateFieldSectionList, &
       DestroyFieldSectionList, &
       InsertAtFieldSectionList, &
       DumpFieldSectionList

  use var_tables, only: &
       var_tables_r
  

  implicit none
  private
  public :: MessageData
  public :: InitializeMessageData
  public :: TransferMessageData
  public :: InsertFieldSectionAtSendRecvMessageData
  public :: CleanMessageData
  public :: DestroyMessageData

  ! data to send/receive to/from one node in one message

  type MessageData
     real, allocatable :: buf(:)  ! message buffer
     integer :: bufSize=0         ! message buffer size
     type (FieldSectionList), pointer :: fieldList=> null() ! field sections to communicate
  end type MessageData

contains






  subroutine InitializeMessageData (MsgData)

    ! Create MessageData components of each entry
    ! of a previously allocated MesageData array 

    type(MessageData), intent(inout) :: MsgData(:)

    integer :: i
    character(len=*), parameter :: h="**(InitializeMessageData)**"
    logical, parameter :: dumpLocal=.false.
    
    do i = 1, size(MsgData)
       MsgData(i)%bufSize = 0
       MsgData(i)%fieldList => CreateFieldSectionList()
    end do
    if (dumpLocal) then
       call MsgDump(h//" done")
    end if
  end subroutine InitializeMessageData






  subroutine TransferMessageData(left, right)

    ! assignment of MessageData variables

    type(MessageData) :: left
    type(MessageData) :: right
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(TransferMessageData)**"

    left%bufSize = right%bufSize
    left%fieldList => right%fieldList
    if (dumpLocal) then
       call MsgDump(h//" left => right")
    end if
  end subroutine TransferMessageData






  subroutine InsertFieldSectionAtSendRecvMessageData(&
       vTabPtr, ParEnv, Neigh, GlobalWithGhost, &
       xbSend, xeSend, ybSend, yeSend, willSend, SendMsgData, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMsgData)

    ! Insert, at previosly existing message data variable,
    ! field sections to communicate. 
    ! The field sections are arrays indexed by number of neighbours.
    ! Logical arrays indicate if each neighbour has field
    ! sections to communicate.
    ! Insertion is performed only for neighbours that have
    ! field sections to communicate.

    type(var_tables_r), pointer, intent(in) :: vTabPtr
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! all remaining arguments are dimensioned by number of neighbours 
    ! and indexed by neighbour number

    ! region to be sent to each neighbour (global indices)

    integer, intent(in) :: xbSend(:)
    integer, intent(in) :: xeSend(:)
    integer, intent(in) :: ybSend(:)
    integer, intent(in) :: yeSend(:)

    ! which neighbours will receive msgs from this node

    logical, intent(in) :: willSend(:)

    ! potential sending message data

    type(MessageData), intent(inout) :: SendMsgData(:)

    ! region to be received from each neighbour (global indices)

    integer, intent(in) :: xbRecv(:)
    integer, intent(in) :: xeRecv(:)
    integer, intent(in) :: ybRecv(:)
    integer, intent(in) :: yeRecv(:)

    ! which neighbours will send msgs to this node

    logical, intent(in) :: willRecv(:)

    ! potential receiving message data

    type(MessageData), intent(inout) :: RecvMsgData(:)

    integer :: i
    integer :: thisNode
    integer :: x0, y0
    type(FieldSection), pointer :: oneFieldSection
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(InsertFieldSectionAtSendRecvMessageData)**"
    logical, parameter :: dumpLocal=.false.

    ! check arguments

    if (.not. associated(vTabPtr)) then
       call fatal_error(h//" vTabPtr not associated")
    else if (.not. associated(ParEnv)) then
       call fatal_error(h//" ParEnv not associated")
    else if (.not. associated(Neigh)) then
       call fatal_error(h//" Neigh not associated")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" GlobalWithGhost not associated")
    end if

    ! offsets to convert global indices to local indices at this proc

    thisNode = ParEnv%myNum
    x0 = GlobalWithGhost%xb(thisNode) - 1
    y0 = GlobalWithGhost%yb(thisNode) - 1

    ! create list of Field Sections to send and to receive

    do i = 1, Neigh%nNeigh
       if (willSend(i)) then
          select case (vTabPtr%idim_type)
          case (2)
             oneFieldSection =>  CreateFieldSection(&
                  vTabPtr%var_p_2D, &
                  vTabPtr%name, &
                  vTabPtr%idim_type, &
                  xbSend(i)-x0, xeSend(i)-x0, &
                  ybSend(i)-y0, yeSend(i)-y0)
          case (3,6,7)
             oneFieldSection =>  CreateFieldSection(&
                  vTabPtr%var_p_3D, &
                  vTabPtr%name, &
                  vTabPtr%idim_type, &
                  xbSend(i)-x0, xeSend(i)-x0, &
                  ybSend(i)-y0, yeSend(i)-y0)
          case (4,5)
             oneFieldSection =>  CreateFieldSection(&
                  vTabPtr%var_p_4D, &
                  vTabPtr%name, &
                  vTabPtr%idim_type, &
                  xbSend(i)-x0, xeSend(i)-x0, &
                  ybSend(i)-y0, yeSend(i)-y0)
          case default
             write(c0,"(i8)") vTabPtr%idim_type
             call fatal_error(h//" unknown idim_type="//trim(adjustl(c0)))
          end select
          call InsertAtFieldSectionList(oneFieldSection, &
               SendMsgData(i)%fieldList)
          SendMsgData(i)%bufSize = SendMsgData(i)%bufSize+oneFieldSection%fieldSectionSize
          if (dumpLocal) then
             write(c0,"(i8)") i
             write(c1,"(i8)") SendMsgData(i)%bufSize
             call MsgDump(h//" increased buffer size of SendMsgData("//&
                  trim(adjustl(c0))//") to "//trim(adjustl(c1)))
          end if
       end if
       if (willRecv(i)) then
          select case (vTabPtr%idim_type)
          case (2)
             oneFieldSection =>  CreateFieldSection(&
                  vTabPtr%var_p_2D, &
                  vTabPtr%name, &
                  vTabPtr%idim_type, &
                  xbRecv(i)-x0, xeRecv(i)-x0, &
                  ybRecv(i)-y0, yeRecv(i)-y0)
          case (3,6,7)
             oneFieldSection =>  CreateFieldSection(&
                  vTabPtr%var_p_3D, &
                  vTabPtr%name, &
                  vTabPtr%idim_type, &
                  xbRecv(i)-x0, xeRecv(i)-x0, &
                  ybRecv(i)-y0, yeRecv(i)-y0)
          case (4,5)
             oneFieldSection =>  CreateFieldSection(&
                  vTabPtr%var_p_4D, &
                  vTabPtr%name, &
                  vTabPtr%idim_type, &
                  xbRecv(i)-x0, xeRecv(i)-x0, &
                  ybRecv(i)-y0, yeRecv(i)-y0)
          case default
             write(c0,"(i8)") vTabPtr%idim_type
             call fatal_error(h//" unknown idim_type="//trim(adjustl(c0)))
          end select
          call InsertAtFieldSectionList(oneFieldSection, &
               RecvMsgData(i)%fieldList)
          RecvMsgData(i)%bufSize = RecvMsgData(i)%bufSize+oneFieldSection%fieldSectionSize
          if (dumpLocal) then
             write(c0,"(i8)") i
             write(c1,"(i8)") RecvMsgData(i)%bufSize
             call MsgDump(h//" increased buffer size of RecvMsgData("//&
                  trim(adjustl(c0))//") to "//trim(adjustl(c1)))
          end if
       end if
    end do
  end subroutine InsertFieldSectionAtSendRecvMessageData





  subroutine CleanMessageData (MsgData)
    type(MessageData), intent(inout) :: MsgData(:)

    ! To clean a message data array,
    ! destroy field section lists of null entries,
    ! create fresh areas for non-empty field section lists,
    ! since these non-empty entries are pointed by someone else

    integer :: i
    character(len=*), parameter :: h="**(CleanMessageData)**"
    character(len=8) :: c0, c1
    logical, parameter :: dumpLocal=.false.

    do i = 1, size(MsgData)
       if (MsgData(i)%bufSize == 0) then
          if (dumpLocal) then
             write(c0,"(i8)") i
             write(c1,"(i8)") MsgData(i)%bufSize
             call MsgDump(h//" will destroy field section list of MsgData("//&
                  trim(adjustl(c0))//") of size "//trim(adjustl(c1)))
          end if
          call DestroyFieldSectionList(MsgData(i)%fieldList)
       else
          if (dumpLocal) then
             write(c0,"(i8)") i
             write(c1,"(i8)") MsgData(i)%bufSize
             call MsgDump(h//" sets buffer size of MsgData("//&
                  trim(adjustl(c0))//") of size "//trim(adjustl(c1))//&
                  " to zero and nullify fieldList pointer")
          end if
          MsgData(i)%bufSize = 0
          MsgData(i)%fieldList => null()
       end if
    end do
  end subroutine CleanMessageData






  subroutine DestroyMessageData(msgData)

    ! destroy a variable of this type

    type(MessageData) :: msgData

    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(DestroyMessageData)**"
    logical, parameter :: dumpLocal=.false.

    msgData%bufSize = 0
    if (allocated(msgData%buf)) then
       deallocate(msgData%buf, stat=ierr)
          if (ierr /= 0) then
             write(c0,"(i8)") ierr
             write(c1,"(i8)") msgData%bufSize
             call fatal_error(h//" deallocate msgData%buf("//&
                  trim(adjustl(c1))//") fails with stat="//&
                  trim(adjustl(c0)))
          end if
    end if
    if (dumpLocal) then
       call MsgDump(h//" nullify bufSize and deallocate buf of MsgData;"//&
            " will destroy field section list")
    end if
    call DestroyFieldSectionList(msgData%fieldList)
  end subroutine DestroyMessageData
end module ModMessageData
