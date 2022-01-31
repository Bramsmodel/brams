module ModMessageSetSendRecv

  ! At each timestep, the ghost zone of some fixed fields 
  ! needs updating at fixed points of the computation.
  ! Ghost zone update requires message passing.
  ! The set of fields to update varies with the computation point.
  !
  ! For each point of the computation the set of fields
  ! and the ranks to access are stored at a pair of variables
  ! of type MessageSet, one for send and one for receive operations.
  !
  ! Type MessageSet and procedures to create, destroy, dump
  ! and perform general message passing operations are implemented
  ! by this module. Procedures to create and destroy specific
  ! MessageSet operations required by timestep are implemented
  ! by ModMessageSet, that uses the general operations defined
  ! at this module.
  !
  ! MessageSet variables contain data for a set of message
  ! passing operations. These variables are created during
  ! initialization and used throughout integration.
  !
  ! All message passing operations are nonblocking.
  ! Procedure PostRecvSendMsgs initiate all send/recv
  ! operations of one set. Procedure WaitSendRecvMsgs finalize
  ! all send/recv operations of the same set.
  !
  ! MessageSet type components are described bellow.
  
  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       Brams2MpiProcNbr, &
       Mpi2BramsProcNbr, &
       MsgDump

  use ModMessageData, only: &
       MessageData, &
       TransferMessageData, &
       DestroyMessageData

  use ModNeighbourNodes, only: &
       NeighbourNodes

  use ModFieldSectionList, only: &
       FieldSection, &
       CreateFieldSection, &
       NextFieldSection, &
       DumpFieldSection, &
       CreateFieldSectionList, &
       InsertAtFieldSectionList, &
       DumpFieldSectionList

  use ParLib, only: &
       parf_get_noblock_real, &
       parf_send_noblock_real, &
       parf_wait_any_nostatus, &
       parf_wait_all_nostatus

  use ModBuffering, only: &
       FieldSection2Buffer, &
       Buffer2FieldSection

  use ModDomainDecomp, only: &
       DomainDecomp

  use var_tables, only: &
       var_tables_r

!  use CUPARM_GRELL3, only: g3d_g

  implicit none
  private
  public :: MessageSet
  public :: CreateMessageSet
  public :: InsertFieldSectionAtMessageSet
  public :: DumpMessageSet
  public :: DestroyMessageSet
  public :: PostRecvSendMsgs
  public :: WaitSendRecvMsgs

  include "mpif.h"

  integer, parameter :: UNDEFINED=-1

  ! all messages to send/receive from one process
  ! to update a field.
  ! arrays are indexed 1:nMsgs

  type MessageSet

     ! A set of message passing operations in
     ! one direction (send or recv):

     character(len=64) :: name
     ! message passing set name
     integer :: nMsgs=UNDEFINED
     ! number of messages on the set,
     ! one message for each rank
     integer :: tag
     ! same tag for all messages in the set
     type(MessageData), allocatable :: oneMsg(:)
     ! data to communicate at each message;
     ! array of MessageData is indexed by message
     ! number on the set. See file ModMessageSet.f90
     ! for the declaration of MessageSet type
     integer, allocatable :: request(:)
     ! communication request; array of requests is
     ! indexed by message number on the set
     integer, allocatable :: otherProc(:)
     ! MPI rank to communicate; array of ranks
     ! is indexed by message number on the set
  end type MessageSet

contains





  function CreateMessageSet (name, tag, hasMsg, Neigh) result(Msgs)

    ! Generates a variable of type MessageSet containing
    ! all message envelopes and no message data
    ! to be sent by this node to neighbour nodes or
    ! to be received by this node from neighbour nodes
    ! during the communication denoted by "name".
    ! Input includes "hasMsg", a logical array indexed by
    ! neighbour number that stores if a neighbour will or
    ! will not communicate with this node in this communication.
    ! Same tag should be used on send and recieve.
    ! Returning variable has as many messages as the number
    ! of true values in hasMsg. Arrays in returning variable
    ! are indexed by true value count, at range 1:nMsgs.

    character(len=*), intent(in) :: name
    integer, intent(in) :: tag
    logical, intent(in) :: hasMsg(:)
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(MessageSet), pointer :: Msgs

    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(CreateMessageSet)**"

    integer :: nNeigh
    integer :: nMsgs
    integer :: iNeigh
    integer :: ierr
    integer :: cntMsg

    ! no message if no neighbours or no node has messages

    nNeigh = Neigh%nNeigh
    nMsgs = count(hasMsg)
    if (nNeigh == 0 .or. nMsgs == 0) then
       Msgs => null()
       return
    end if

    ! there are messages: allocate area

    allocate(Msgs, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate Msgs fails with stat="//&
            trim(adjustl(c0)))
    end if

    allocate(Msgs%oneMsg(nMsgs), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") nMsgs
       call fatal_error(h//" allocate Msgs%oneMsg("//trim(adjustl(c1))//&
            ") fails with stat="//trim(adjustl(c0)))
    end if

    allocate(Msgs%request(nMsgs), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") nMsgs
       call fatal_error(h//" allocate Msgs%request("//trim(adjustl(c1))//&
            ") fails with stat="//trim(adjustl(c0)))
    end if

    allocate(Msgs%otherProc(nMsgs), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") nMsgs
       call fatal_error(h//" allocate Msgs%otherProc("//trim(adjustl(c1))//&
            ") fails with stat="//trim(adjustl(c0)))
    end if

    ! store basic info

    Msgs%name = name
    Msgs%nMsgs = nMsgs
    Msgs%tag = tag

    ! for each neighbour node that will communicate,
    ! build message envelop to/from neighbour and
    ! create empty message data.

    cntMsg=0
    do iNeigh = 1, nNeigh
       if (hasMsg(iNeigh)) then
          cntMsg = cntMsg + 1
          Msgs%oneMsg(cntMsg)%bufSize=0
          Msgs%oneMsg(cntMsg)%fieldList => CreateFieldSectionList()
          Msgs%request(cntMsg) = MPI_REQUEST_NULL
          Msgs%otherProc(cntMsg)= Brams2MpiProcNbr(Neigh%neigh(iNeigh))
       end if
    end do

    if (cntMsg /= nMsgs) then
       write(c0,"(i8)") cntMsg
       write(c1,"(i8)") nMsgs
       call fatal_error(h//" inconsistency: cntMsg="//trim(adjustl(c0))//&
            " while nMsgs="//trim(adjustl(c1)))
    end if
  end function CreateMessageSet





  subroutine InsertFieldSectionAtMessageSet(&
       myNum, vTabPtr, Neigh, GlobalWithGhost, &
       xbComm, xeComm, ybComm, yeComm, willComm, &
       Msgs)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable
    
    integer, intent(in) :: myNum
    
    type(var_tables_r), pointer, intent(in) :: vTabPtr

    ! Neigh is an array of neighbours (to this process)
    ! BRAMS process numbers to communicate;
    ! Each neighbour process may (or may not) exchange
    ! messages with this process on a MessageSet variable;
    ! The number of processes for potential communication
    ! is the size of the Neigh array

    type(NeighbourNodes), pointer, intent(in) :: Neigh

    ! Global indices of domain partition with Ghost Zones
    
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! all remaining arguments are dimensioned by number of neighbours
    ! and indexed up to the size of Neigh array

    ! global indices of the region of this process field
    ! to be sent on a send MessageSet variable or to be received
    ! on a receive MessageSet variable;
    ! The arrays of global indices have the same size of the Neigh
    ! array (number of processes to communicate) and are indexed
    ! accordingly

    integer, intent(in) :: xbComm(:)
    integer, intent(in) :: xeComm(:)
    integer, intent(in) :: ybComm(:)
    integer, intent(in) :: yeComm(:)

    ! which neighbours (BRAMS process number) will
    ! receive msgs from this node on a send MessageSet variable
    ! or will send msgs to this node on a receive Message set
    ! variable

    logical, intent(in) :: willComm(:)

    ! MessageSet variable to be updated by field section inclusion

    type(MessageSet), pointer, intent(inout) :: Msgs

    integer :: nMsgs
    integer :: x0, y0
    integer :: cntMsg
    integer :: iNeigh
    type(FieldSection), pointer :: oneFieldSection
    character(len=8) :: c0
    character(len=*), parameter :: h="**(InsertFieldSectionAtMessageSet)**"

    ! check arguments

    if (.not. associated(vTabPtr)) then
       call fatal_error(h//" vTabPtr not associated")
    else if (.not. associated(Neigh)) then
       call fatal_error(h//" Neigh not associated")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" GlobalWithGhost not associated")
    end if

    ! return if no messages to send

    if (.not. associated(Msgs)) then
       return
    end if
    nMsgs = Msgs%nMsgs

    ! offsets to convert global indices to local indices at this proc

    x0 = GlobalWithGhost%xb(myNum) - 1
    y0 = GlobalWithGhost%yb(myNum) - 1

    ! create list of Field Sections to communicate, one for
    ! each process to communicate and insert at this MessageSet
    ! field section list

    cntMsg = 0
    do iNeigh = 1, Neigh%nNeigh
       if (willComm(iNeigh)) then
          cntMsg = cntMsg + 1
          if (cntMsg > nMsgs) then
             write(c0,"(i8)") nMsgs
             call fatal_error(h//" nMsgs ("//&
                  trim(adjustl(c0))//") exceeded while inserting field "//&
                  trim(adjustl(vTabPtr%name))//&
                  " at message "//trim(adjustl(Msgs%name)))
          end if
          oneFieldSection =>  CreateFieldSection(&
               vTabPtr, &
               xbComm(iNeigh)-x0, xeComm(iNeigh)-x0, &
               ybComm(iNeigh)-y0, yeComm(iNeigh)-y0, &
               GlobalWithGhost)
          call InsertAtFieldSectionList(&
               oneFieldSection, &
               Msgs%oneMsg(cntMsg)%fieldList)
          Msgs%oneMsg(cntMsg)%bufSize = Msgs%oneMsg(cntMsg)%bufSize + &
               oneFieldSection%fieldSectionSize
       end if
    end do
  end subroutine InsertFieldSectionAtMessageSet




  subroutine DumpMessageSet(Msgs)

    ! dumps a message set variable
    
    type(MessageSet), pointer, intent(in) :: Msgs

    character(len=8) :: c0, c1, c2, c3, c4
    character(len=*), parameter :: h="**(DumpMessageSet)**"
    integer :: i

    if (.not. associated(Msgs)) then
       call MsgDump(h//" empty messages")
    else
       write(c0,"(i8)") Msgs%nMsgs
       if (Msgs%tag == UNDEFINED) then
          c2="UNDEF"
       else
          write(c2,"(i8)") Msgs%tag
       end if
       call MsgDump(h//" "//trim(adjustl(Msgs%name))//&
            " has "//trim(adjustl(c0))//" messages and tag "//&
            trim(adjustl(c2)))
       do i = 1, Msgs%nMsgs
          write(c0,"(i8)") i
          if (Msgs%otherProc(i) == UNDEFINED) then
             c1="UNDEF"
          else
             write(c1,"(i8)") Msgs%otherProc(i)
          end if
          if (Msgs%oneMsg(i)%bufSize == UNDEFINED) then
             c3="UNDEF"
          else
             write(c3,"(i8)") Msgs%oneMsg(i)%bufSize
          end if
          if (Msgs%request(i) == MPI_REQUEST_NULL) then
             c4="NULL"
          else
             write(c4,"(Z8)") Msgs%request(i)
          end if
          call MsgDump(h//" message number "//trim(adjustl(c0))//&
               " to/from MPI node "//trim(adjustl(c1))//&
               ", request "//trim(adjustl(c4))//&
               ", size "//trim(adjustl(c3))//&
               " and field sections:")
          call DumpFieldSectionList(Msgs%oneMsg(i)%fieldList)
       end do
    end if
  end subroutine DumpMessageSet





  subroutine DestroyMessageSet (Msgs)

    ! deallocate memory of a message set variable
    
    type(MessageSet), pointer, intent(inout) :: Msgs

    integer :: msg
    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DestroyMessageSet)**"

    if (associated(Msgs)) then
       do msg = 1, Msgs%nMsgs
          call DestroyMessageData(Msgs%oneMsg(msg))
       end do
       deallocate(Msgs%oneMsg, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate Msgs%oneMsg fails with stat="//&
               trim(adjustl(c0)))
       end if

       deallocate(Msgs%request, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate Msgs%request fails with stat="//&
               trim(adjustl(c0)))
       end if

       deallocate(Msgs%otherProc, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate Msgs%otherProc fails with stat="//&
               trim(adjustl(c0)))
       end if

       deallocate(Msgs, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate Msgs fails with stat="//&
               trim(adjustl(c0)))
       end if

    end if
    Msgs => null()
  end subroutine DestroyMessageSet




  subroutine PostRecvSendMsgs(SendMsg, RecvMsg)

    ! posts all nonblocking send and recv operations of
    ! a message set pair of variables
    
    type(MessageSet), pointer, intent(in) :: SendMsg
    type(MessageSet), pointer, intent(in) :: RecvMsg

    integer :: iSend
    integer :: iRecv
    integer :: firstBuffer
    integer :: lastBuffer
    integer :: ierr
    type(MessageData), pointer :: msgData => null()
    type(FieldSection), pointer :: node => null()
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(PostRecvSendMsgs)**"
    logical, parameter :: dumpLocal=.false.

    ! post nonblocking receive for each receiving message;
    ! a single receive msg from each process

    if (associated(RecvMsg)) then
       if (dumpLocal) then
          if (RecvMsg%nMsgs > 0) then
             write(c0,"(i8)") RecvMsg%nMsgs
             call MsgDump(h//" for "//trim(adjustl(RecvMsg%name))//&
                  " will post "//trim(adjustl(c0))//&
                  " nonblocking receives")
          end if
       end if
       do iRecv= 1,RecvMsg%nMsgs

          ! data to receive from one process
          
          msgData => RecvMsg%oneMsg(iRecv)

          ! allocate receive buffer for this receive

          allocate(msgData%buf(msgData%bufSize), stat=ierr)
          if (ierr /= 0) then
             write(c0,"(i8)") ierr
             write(c1,"(i8)") msgData%bufSize
             call fatal_error(h//" allocate msgData%buf("//&
                  trim(adjustl(c1))//") to recv fails with stat="//&
                  trim(adjustl(c0)))
          end if

          ! post receive

          call parf_get_noblock_real(msgData%buf, msgData%bufSize, &
               RecvMsg%otherProc(iRecv), &
               RecvMsg%tag, RecvMsg%request(iRecv))

          if (dumpLocal) then
             write(c0,"(i8)") iRecv
             write(c1,"(i8)") RecvMsg%otherProc(iRecv)
             write(c2,"(i8)") size(msgData%buf)
             write(c3,"(i8)") RecvMsg%tag
             if (RecvMsg%request(iRecv) == MPI_REQUEST_NULL) then
                c4="NULL"
             else
                write(c4,"(Z8)") RecvMsg%request(iRecv)
             end if
             call MsgDump(h//" for "//trim(adjustl(RecvMsg%name))//&
                  " post recv number "//trim(adjustl(c0))//&
                  " from MPI node "//trim(adjustl(c1))//&
                  " with buffer size "//trim(adjustl(c2))//&
                  " tag "//trim(adjustl(c3))//" and request "//trim(adjustl(c4)))
          end if
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" no message to receive")
       end if
    end if

    ! for each sending message,
    ! build send buffer and copy field sections to the buffer;
    ! post nonblocking send;
    ! A single send message to each process

    if (associated(SendMsg)) then
       if (dumpLocal) then
          if (SendMsg%nMsgs > 0) then
             write(c0,"(i8)") SendMsg%nMsgs
             call MsgDump(h//" for "//trim(adjustl(SendMsg%name))//&
                  " will post "//trim(adjustl(c0))//&
                  " nonblocking sends")
          end if
       end if
       do iSend = 1,SendMsg%nMsgs

          ! data to send
          
          msgData => SendMsg%oneMsg(iSend)

          ! allocate send buffer
          
          allocate(msgData%buf(msgData%bufSize), stat=ierr)
          if (ierr /= 0) then
             write(c0,"(i8)") ierr
             write(c1,"(i8)") msgData%bufSize
             call fatal_error(h//" allocate msgData%buf("//&
                  trim(adjustl(c1))//") to send fails with stat="//&
                  trim(adjustl(c0)))
          end if

          ! fill send buffer with field sections to send
          
          node => null()
          lastBuffer=0
          if (dumpLocal) then
             write(c0,"(i8)") iSend
             call MsgDump(h//" starts building send buffer for send number "//&
                  trim(adjustl(c0)))
          end if
          do
             node => NextFieldSection(node, msgData%fieldList)
             if (associated(node)) then
                if (dumpLocal) then
                   call MsgDump(h//" insert FieldSection at send buffer")
                   call DumpFieldSection(node)
                end if
                if (associated(node%vTabPtr%var_p_2D)) then
                   firstBuffer = lastBuffer+1
                   call FieldSection2Buffer(node%vTabPtr%var_p_2D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      call MsgDump(h//" filled buf["//trim(adjustl(c0))//&
                           ":"//trim(adjustl(c1))//&
                           "] with 2D field "//trim(adjustl(node%vTabPtr%name))//"["//&
                           trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                           trim(adjustl(c4))//":"//trim(adjustl(c5))//"]")
                   end if

                else if (associated(node%vTabPtr%var_p_3D)) then
                   firstBuffer = lastBuffer+1
                   call FieldSection2Buffer(node%vTabPtr%var_p_3D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      select case (node%vTabPtr%idim_type)
                      case(3)
                         call MsgDump(h//" filled buf["//trim(adjustl(c0))//&
                              ":"//trim(adjustl(c1))//&
                              "] with 3D field "//trim(adjustl(node%vTabPtr%name))//"[:,"//&
                              trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                              trim(adjustl(c4))//":"//trim(adjustl(c5))//"]")
                      case(6:7)
                         call MsgDump(h//" filled buf["//trim(adjustl(c0))//&
                              ":"//trim(adjustl(c1))//&
                              "] with 3D field "//trim(adjustl(node%vTabPtr%name))//"["//&
                              trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                              trim(adjustl(c4))//":"//trim(adjustl(c5))//",:]")
                      end select
                   end if

                else if (associated(node%vTabPtr%var_p_4D)) then
                   firstBuffer = lastBuffer+1
                   call FieldSection2Buffer(node%vTabPtr%var_p_4D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      call MsgDump(h//" filled buf["//trim(adjustl(c0))//&
                           ":"//trim(adjustl(c1))//&
                           "] with 4D field "//trim(adjustl(node%vTabPtr%name))//"[:,"//&
                           trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                           trim(adjustl(c4))//":"//trim(adjustl(c5))//",:]")
                   end if
                else
                   call fatal_error(h//" inconsistent var_tables entry named "//&
                        trim(adjustl(node%vTabPtr%name)))
                end if
             else
                if (dumpLocal) then
                   call MsgDump(h//" finishes building sending buffer")
                end if
                exit
             end if
          end do
          if (lastBuffer /= size(msgData%buf)) then
             write(c0,"(i8)") lastBuffer
             write(c1,"(i8)") size(msgData%buf)
             call fatal_error(h//" the send buffer of size "//trim(adjustl(c1))//&
                  " was filled with "//trim(adjustl(c0))//" entries")
          end if

          ! post send message
          
          call parf_send_noblock_real(msgData%buf, msgData%bufSize, &
               SendMsg%otherProc(iSend), &
               SendMsg%tag, SendMsg%request(iSend))
          if (dumpLocal) then
             write(c1,"(i8)") SendMsg%otherProc(iSend)
             write(c2,"(i8)") size(msgData%buf)
             write(c3,"(i8)") SendMsg%tag
             if (SendMsg%request(iSend) == MPI_REQUEST_NULL) then
                c4="NULL"
             else
                write(c4,"(Z8)") SendMsg%request(iSend)
             end if
             call MsgDump(h//" sends to MPI node "//trim(adjustl(c1))//&
                  " buffer of size "//trim(adjustl(c2))//&
                  " tag "//trim(adjustl(c3))//" and request "//trim(adjustl(c4)))
          end if
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" no message to send")
       end if
    end if
  end subroutine PostRecvSendMsgs




  subroutine WaitSendRecvMsgs(SendMsg, RecvMsg)
    type(MessageSet), pointer, intent(in) :: SendMsg
    type(MessageSet), pointer, intent(in) :: RecvMsg

    ! waits for all nonblocking send and recv operations of
    ! a message set pair of variables
    
    integer :: i
    integer :: iSend
    integer :: iRecv
    integer :: firstBuffer
    integer :: lastBuffer
    integer :: recvNbr
    integer :: sendNbr
    integer :: ierr
    type(MessageData), pointer :: msgData => null()
    type(FieldSection), pointer :: node => null()
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(WaitSendRecvMsgs)**"
    logical, parameter :: dumpLocal=.false.

    ! for each receive message:
    ! build send buffer and copy field sections to the buffer;
    ! post nonblocking send;
    ! A single send message to each process

    if (associated(RecvMsg)) then
       if (dumpLocal) then
          write(c0,"(i8)") RecvMsg%nMsgs
          call MsgDump(h//" for "//trim(adjustl(RecvMsg%name))//&
               " waits on "//trim(adjustl(c0))//" receives")
       end if

       do iRecv= 1,RecvMsg%nMsgs

          ! wait on any arrived message

          call parf_wait_any_nostatus(RecvMsg%nMsgs, &
               RecvMsg%request, recvNbr)
          msgData => RecvMsg%oneMsg(recvNbr)
          if (dumpLocal) then
             write(c0,"(i8)") recvNbr
             write(c1,"(i8)") RecvMsg%otherProc(recvNbr)
             call MsgDump(h//" received message #"//trim(adjustl(c0))//&
                  " from MPI proc "//trim(adjustl(c1)))
          end if

          ! extract field sections from incoming buffer
          ! and store at destination fields

          node => null()
          lastBuffer=0
          do
             node => NextFieldSection(node, msgData%fieldList)
             if (associated(node)) then
                if (associated(node%vTabPtr%var_p_2D)) then
                   firstBuffer = lastBuffer+1
                   call Buffer2FieldSection(node%vTabPtr%var_p_2D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)
                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      call MsgDump(h//&
                           "filled 2D field "//trim(adjustl(node%vTabPtr%name))//"["//&
                           trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                           trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                           " with buf["//trim(adjustl(c0))//&
                           ":"//trim(adjustl(c1))//"]")
                   end if
                else if (associated(node%vTabPtr%var_p_3D)) then
                   firstBuffer = lastBuffer+1
                   call Buffer2FieldSection(node%vTabPtr%var_p_3D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      select case (node%vTabPtr%idim_type)
                      case(3)
                         call MsgDump(h//&
                              " filled 3D field "//trim(adjustl(node%vTabPtr%name))//"[:,"//&
                              trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                              trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                              " with buf["//trim(adjustl(c0))//&
                              ":"//trim(adjustl(c1))//"]")
                      case(6:7)
                         call MsgDump(h//&
                              " filled 3D field "//trim(adjustl(node%vTabPtr%name))//"["//&
                              trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                              trim(adjustl(c4))//":"//trim(adjustl(c5))//",:]"//&
                              " with buf["//trim(adjustl(c0))//&
                              ":"//trim(adjustl(c1))//"]")
                      end select
                   end if
                else if (associated(node%vTabPtr%var_p_4D)) then
                   firstBuffer = lastBuffer+1
                   call Buffer2FieldSection(node%vTabPtr%var_p_4D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      call MsgDump(h//&
                           " filled 4D field "//trim(adjustl(node%vTabPtr%name))//"[:,"//&
                           trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                           trim(adjustl(c4))//":"//trim(adjustl(c5))//",:]"//&
                           " with buf["//trim(adjustl(c0))//&
                           ":"//trim(adjustl(c1))//"]")
                   end if
                else
                   call fatal_error(h//" inconsistent var_tables entry")
                end if
             else
                exit
             end if
          end do

          ! done with this message; deallocate buffer

          deallocate(msgData%buf, stat=ierr)
          if (ierr /= 0) then
             write(c0,"(i8)") ierr
             write(c1,"(i8)") msgData%bufSize
             call fatal_error(h//" deallocate msgData%buf("//&
                  trim(adjustl(c1))//") to recv fails with stat="//&
                  trim(adjustl(c0)))
          end if

       end do
    end if

    ! for all posted send messages, wait on pending request,
    ! deallocate buffer and empty request

    if (associated(SendMsg)) then
!CDIR$ NOVECTOR
       do iSend = 1,SendMsg%nMsgs
          call parf_wait_any_nostatus(SendMsg%nMsgs, &
               SendMsg%request, sendNbr)
          msgData => SendMsg%oneMsg(sendNbr)
          deallocate(msgData%buf, stat=ierr)
          if (ierr /= 0) then
             write(c0,"(i8)") ierr
             write(c1,"(i8)") msgData%bufSize
             call fatal_error(h//" deallocate msgData%buf("//&
                  trim(adjustl(c1))//") to send fails with stat="//&
                  trim(adjustl(c0)))
          end if

       end do
    end if
  end subroutine WaitSendRecvMsgs
end module ModMessageSetSendRecv
