module ModMessageSet

  ! At each timestep, the ghost zone of some fields has to be
  ! updated. Ghost zone update requires message passing.
  ! There are a few code locations that need ghost zone
  ! update. The set of fields to update is fixed at one
  ! location, but varies with the location.
  !
  ! For each location that requires message passing, message
  ! data and message envelop are stored at a pair of variables
  ! of type MessageSet, one for send and one for receive operations.
  !
  ! Each MessageSet variable contains message data and message envelop
  ! for an array of message passing operations. Each array element
  ! contains all data to pass to a single rank, on a single
  ! message passing operation. Consequently, field sections of
  ! all fields that require ghost zone update involving two fixed
  ! ranks are stored at a single array entry of a MessageSet variable
  ! and used on a single message passing operation. 
  !
  ! Message passing operates on contiguous buffer. Field sections
  ! to update are not contiguous in memory and have to be copied
  ! to/from the buffer. Mapping of a field section to buffer locations
  ! and inverse mapping are computed once, at MessageSet creation
  ! and stored at the corresponding MessageSet array entry.
  !
  ! Type MessageSet and procedures to create, destroy, dump
  ! and perform general message passing operations are implemented
  ! by this module.
  !
  ! Procedures to create and destroy specific MessageSet operations
  ! (those required by current code version) are implemented by this module.
  !
  ! Message passing operations are also implemented by this module.
  ! All message passing operations are nonblocking.
  ! Procedure PostSendRecvMsgs initiate all send/recv
  ! operations of one set. Procedure WaitSendRecvMsgs finalize
  ! all send/recv operations of the same set.
  ! Copy field sections to the contiguous buffer and the reverse
  ! opetation are included on both message passing procedures.
  !
  ! To overlap computation with communication, procedure
  ! PostSendRecvMsgs should be invoked as early as possible,
  ! mainly when all fields to send are fully computed,
  ! and procedure WaitSendRecvMsgs should be invoked as
  ! late as possible, mainly when the ghost zone of the
  ! updated fields will be used in the computation.


  use ModGridDims, only: &
       GridDims

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       Brams2MpiProcNbr, &
       MsgDump

  use ModNeighbourNodes, only: &
       NeighbourNodes, &
       NodesToSendRecvMessages, &
       IncludeDomainBoundaries, &
       DumpNeighbourNodes

  use ModDomainDecomp, only: &
       DomainDecomp

  use ModFieldSection, only: &
       FieldSection, &
       CreateFieldSection

  use ModMessageData, only: &
       MessageData, &
       CreateMessageData, &
       DumpMessageData, &
       AppendFieldSectionToMessageData, &
       PostRecvMessageData, &
       PostSendMessageData, &
       FillMessageDataBuffer, &
       FillMessageDataBufferVariableAdressArr, &
       FillMessageDataBufferVariableAdressScalar, &
       DecomposeMessageDataBuffer, &
       AllocateMessageDataBuffer, &
       DeallocateMessageDataBuffer, &
       DestroyMessageData, &
       UpdateFieldAdress

  use var_tables, only: &
       var_tables_r, &
       GetVTabEntry

  use ModNamelistFile, only: &
       NamelistFile

  use mem_grid, only : &
       dyncore_flag

  use ParLib, only: &
       parf_wait_any_nostatus

  implicit none
  include "mpif.h"

  private

  public :: MessageSet
  public :: DumpMessageSet

  public :: CreateAcousticMessageSet
  public :: DestroyAcousticMessageSet

  public :: CreateDn0MessageSet
  public :: DestroyDn0MessageSet

  public :: CreateG3DMessageSet
  public :: DestroyG3DMessageSet

  public :: CreateSelectedGhostZoneMessageSet
  public :: DestroySelectedGhostZoneMessageSet

  public :: CreateAllGhostZoneMessageSet
  public :: DestroyAllGhostZoneMessageSet

  public :: CreateAcouDampOneMessageSet
  public :: DestroyAcouDampOneMessageSet

  public :: CreateWideGhostZoneMessageSet
  public :: DestroyWideGhostZoneMessageSet

  public :: UpdateFieldAdress

  public :: PostSendRecvMsgs
  public :: PostSendRecvMsgsVariableAdress
  public :: WaitSendRecvMsgs
  
  character(len=*), parameter :: sendDirection="send"
  character(len=*), parameter :: recvDirection="recv"

  ! all messages to send/receive from one process
  ! to update a field.
  ! arrays are indexed 1:nMsgs

  type MessageSet
     private
     ! A set of message passing operations in
     ! one direction (send or recv):

     character(len=64) :: name
     ! message passing set name
     integer :: nMsgs=-1
     ! number of messages on the set,
     ! one message for each other rank involved
     ! in message passing with this rank
     integer :: tag
     ! same tag for all messages in the set; part
     ! of message passing envelop
     integer, allocatable :: request(:)
     ! communication request; array of requests is
     ! indexed by message number on the set; each
     ! entry is part of message passing envelop
     integer, allocatable :: otherProc(:)
     ! MPI rank to communicate; array of ranks
     ! is indexed by message number on the set; each
     ! entry is part of message passing envelop
     type(MessageData), allocatable :: msgData(:)
     ! data to communicate at each message;
     ! array of MessageData is indexed by message
     ! number on the set. Each entry contains all message
     ! data to exchange on a single message passing
     ! operation
  end type MessageSet


  interface InsertFieldSectionAtSendRecvMessageSet
     module procedure InsertFieldSectionAtSendRecvMessageSetFromVTab
     module procedure InsertFieldSectionAtSendRecvMessageSet_2D
     module procedure InsertFieldSectionAtSendRecvMessageSet_3D
     module procedure InsertFieldSectionAtSendRecvMessageSet_4D
  end interface InsertFieldSectionAtSendRecvMessageSet

  interface InsertFieldSectionAtMessageSet
     module procedure InsertFieldSectionAtMessageSetFromVTab
     module procedure InsertFieldSectionAtMessageSet_2D  
     module procedure InsertFieldSectionAtMessageSet_3D  
     module procedure InsertFieldSectionAtMessageSet_4D
  end interface InsertFieldSectionAtMessageSet

  interface UpdateFieldAdress
     module procedure UpdateFieldAdressAtMessageSet_2D
     module procedure UpdateFieldAdressAtMessageSet_3D
     module procedure UpdateFieldAdressAtMessageSet_4D
  end interface UpdateFieldAdress

  interface PostSendRecvMsgs
     module procedure PostSendRecvMsgsFixedAdress
  end interface PostSendRecvMsgs

  interface PostSendRecvMsgsVariableAdress
     module procedure PostSendRecvMsgsVariableAdressArr
     module procedure PostSendRecvMsgsVariableAdressScalar
  end interface PostSendRecvMsgsVariableAdress


  interface WaitSendRecvMsgs
     module procedure WaitSendRecvMsgsFixedAdress
     module procedure WaitSendRecvMsgsVariableAdress
  end interface WaitSendRecvMsgs

contains





  function CreateMessageSet (name, direction, tag, hasMsg, Neigh) result(Msgs)

    ! Generates a variable of type MessageSet containing
    ! all message envelopes and no message data
    ! to be sent by this node to all neighbour nodes or
    ! to be received by this node from all neighbour nodes
    ! during the communication denoted by "name".
    !
    ! Input variable Neigh contains all neighbour ranks
    ! that are neighbours to this node, and as so, may
    ! be used on ghost zone update message exchange.
    ! See ModNeighbourNodes for further information.
    !
    ! Input variable hasMsg stores if a neighbour node
    ! will or will not communicate with this node in
    ! this communication. Use of a neighbour node in a
    ! communication depends on the field that requires
    ! ghost zone update. 
    !
    ! Same tag should be used on send and recieve.
    !
    ! Input variable direction is just for dumping.
    !
    ! Returning variable has as many messages as the number
    ! of true values in hasMsg. Arrays in returning variable
    ! are indexed by true value count, at range 1:nMsgs.

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: direction
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

    ! no message if no neighbours or no message to exchange

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

    allocate(Msgs%msgData(nMsgs), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") nMsgs
       call fatal_error(h//" allocate Msgs%msgData("//trim(adjustl(c1))//&
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
          Msgs%request(cntMsg) = MPI_REQUEST_NULL
          Msgs%otherProc(cntMsg)= Brams2MpiProcNbr(Neigh%neigh(iNeigh))
          call CreateMessageData(Msgs%msgData(cntMsg),name,direction)
       end if
    end do

    if (cntMsg /= nMsgs) then
       write(c0,"(i8)") cntMsg
       write(c1,"(i8)") nMsgs
       call fatal_error(h//" inconsistency: cntMsg="//trim(adjustl(c0))//&
            " while nMsgs="//trim(adjustl(c1)))
    end if
  end function CreateMessageSet





  subroutine DestroyMessageSet (Msgs)

    ! deallocate memory of a message set variable

    type(MessageSet), pointer, intent(inout) :: Msgs

    integer :: msg
    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DestroyMessageSet)**"

    if (associated(Msgs)) then
       do msg = 1, Msgs%nMsgs
          call DestroyMessageData(Msgs%msgData(msg))
       end do
       deallocate(Msgs%msgData, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate Msgs%msgData fails with stat="//&
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





  subroutine DumpMessageSet(Msgs)

    ! dumps a message set variable

    type(MessageSet), pointer, intent(in) :: Msgs

    character(len=8) :: c0, c1, c2, c3, c4
    character(len=*), parameter :: h="**(DumpMessageSet)**"
    integer :: i

    if (.not. associated(Msgs)) then
       call MsgDump(h//" empty Message Set")
    else
       write(c0,"(i8)") Msgs%nMsgs
       if (Msgs%tag == -1) then
          c2="UNDEF"
       else
          write(c2,"(i8)") Msgs%tag
       end if
       call MsgDump(h//" "//trim(adjustl(Msgs%name))//&
            " has "//trim(adjustl(c0))//" messages and tag "//&
            trim(adjustl(c2)))
       do i = 1, Msgs%nMsgs
          write(c0,"(i8)") i
          if (Msgs%otherProc(i) == -1) then
             c1="UNDEF"
          else
             write(c1,"(i8)") Msgs%otherProc(i)
          end if
          if (Msgs%request(i) == MPI_REQUEST_NULL) then
             c2="NULL"
          else
             write(c2,"(Z8)") Msgs%request(i)
          end if
          call MsgDump(h//" message number "//trim(adjustl(c0))//&
               " to/from MPI node "//trim(adjustl(c1))//&
               ", has request "//trim(adjustl(c2))//&
               " and message data:")
          call DumpMessageData(Msgs%msgData(i),h)
       end do
    end if
  end subroutine DumpMessageSet



  subroutine NodesRegionsSendRecv(&
       nMachs, nNeigh, myNum, tag, &
       Neigh, GlobalOwn, NameSend, NameRecv, &
       xbToUpdate, xeToUpdate, ybToUpdate, yeToUpdate, &
       xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet)

    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum
    integer, intent(in) :: tag
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    character(len=*), intent(in) :: NameSend
    character(len=*), intent(in) :: NameRecv
    ! global indices to update at each rank
    integer, intent(in) :: xbToUpdate(nMachs)
    integer, intent(in) :: xeToUpdate(nMachs)
    integer, intent(in) :: ybToUpdate(nMachs)
    integer, intent(in) :: yeToUpdate(nMachs)
    ! this rank will send messages to which neighbours
    logical, intent(out) :: willSend(nNeigh)
    ! this rank will send this rectangular region to each neighbour
    integer, intent(out) :: xbSend(nNeigh)
    integer, intent(out) :: xeSend(nNeigh)
    integer, intent(out) :: ybSend(nNeigh)
    integer, intent(out) :: yeSend(nNeigh)
    ! send message set
    type(MessageSet), pointer, intent(inout) :: SendMessageSet
    ! this rank will recv messsages from which neighbours
    logical, intent(out) :: willRecv(nNeigh)
    ! this rank will recv this rectangular region from each neighbour
    integer, intent(out) :: xbRecv(nNeigh)
    integer, intent(out) :: xeRecv(nNeigh)
    integer, intent(out) :: ybRecv(nNeigh)
    integer, intent(out) :: yeRecv(nNeigh)
    ! recv message set
    type(MessageSet), pointer, intent(inout) :: RecvMessageSet

    character(len=*), parameter :: h="**(NodesRegionsSendRecv)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h// "starts")
    end if

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalOwn, &
         xbToUpdate, xeToUpdate, ybToUpdate, yeToUpdate, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
         trim(adjustl(NameSend))//trim(adjustl(NameRecv)))

    ! build message set

    SendMessageSet => CreateMessageSet(&
         NameSend, &
         sendDirection, &
         tag, &
         willSend, &
         Neigh)
    RecvMessageSet => CreateMessageSet(&
         NameRecv, &
         recvDirection, &
         tag, &
         willRecv, &
         Neigh)

    if (dumpLocal) then
       call MsgDump(h//" finishes creating Message Sets:")
       call DumpMessageSet(SendMessageSet)
       call DumpMessageSet(RecvMessageSet)
    end if
  end subroutine NodesRegionsSendRecv



  subroutine InsertFieldSectionAtMessageSetFromVTab(&
       myNum, vTabPtr, nNeigh, GlobalWithGhost, &
       xbComm, xeComm, ybComm, yeComm, willComm, &
       Msgs)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable

    integer, intent(in) :: myNum

    type(var_tables_r), pointer, intent(in) :: vTabPtr

    ! nNeigh is number of processes for potential communication

    integer, intent(in) :: nNeigh

    ! Global indices of domain partition with Ghost Zones

    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! all remaining arguments are dimensioned by nNeigh

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
    character(len=*), parameter :: h="**(InsertFieldSectionAtMessageSetFromVTab)**"
    real, pointer :: PNull(:,:) => null()
    ! check arguments

    if (.not. associated(vTabPtr)) then
       call fatal_error(h//" vTabPtr not associated")
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
    do iNeigh = 1, nNeigh
       if (willComm(iNeigh)) then
          cntMsg = cntMsg + 1
          if (cntMsg > nMsgs) then
             write(c0,"(i8)") nMsgs
             call fatal_error(h//" nMsgs ("//&
                  trim(adjustl(c0))//") exceeded while inserting field "//&
                  trim(adjustl(vTabPtr%name))//&
                  " at message "//trim(adjustl(Msgs%name)))
          end if

          select case (vTabPtr%idim_type)
          case (2)
             oneFieldSection =>  CreateFieldSection(&
                  vTabPtr%var_p_2D, &
                  vTabPtr%name, &
                  vTabPtr%idim_type, &
                  xbComm(iNeigh)-x0, xeComm(iNeigh)-x0, &
                  ybComm(iNeigh)-y0, yeComm(iNeigh)-y0)
          case (3,6,7)
             oneFieldSection =>  CreateFieldSection(&
                  vTabPtr%var_p_3D, &
                  vTabPtr%name, &
                  vTabPtr%idim_type, &
                  xbComm(iNeigh)-x0, xeComm(iNeigh)-x0, &
                  ybComm(iNeigh)-y0, yeComm(iNeigh)-y0)
          case (4,5)
             oneFieldSection =>  CreateFieldSection(&
                  vTabPtr%var_p_4D, &
                  vTabPtr%name, &
                  vTabPtr%idim_type, &
                  xbComm(iNeigh)-x0, xeComm(iNeigh)-x0, &
                  ybComm(iNeigh)-y0, yeComm(iNeigh)-y0)
          case default
             write(c0,"(i8)") vTabPtr%idim_type
             call fatal_error(h//" unknown idim_type="//trim(adjustl(c0)))
          end select
          call AppendFieldSectionToMessageData(oneFieldSection, Msgs%msgData(cntMsg))
       end if
    end do
  end subroutine InsertFieldSectionAtMessageSetFromVTab





  subroutine InsertFieldSectionAtMessageSet_2D(&
       myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
       xbComm, xeComm, ybComm, yeComm, willComm, &
       Msgs)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable

    integer, intent(in) :: myNum

    ! field memory address

    real, pointer, intent(in) :: field(:,:)

    ! field name

    character(len=*), intent(in) :: fieldName

    ! idim_type codes the remained dimensions to communicate 

    integer, intent(in) :: idim_type

    ! nNeigh is number of processes for potential communication

    integer, intent(in) :: nNeigh

    ! Global indices of domain partition with Ghost Zones

    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! all remaining arguments are dimensioned by nNeigh

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
    character(len=*), parameter :: h="**(InsertFieldSectionAtMessageSet_2D)**"
    real, pointer :: PNull(:,:) => null()
    ! check arguments

    if (.not. associated(GlobalWithGhost)) then
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
    do iNeigh = 1, nNeigh
       if (willComm(iNeigh)) then
          cntMsg = cntMsg + 1
          if (cntMsg > nMsgs) then
             write(c0,"(i8)") nMsgs
             call fatal_error(h//" nMsgs ("//&
                  trim(adjustl(c0))//") exceeded while inserting field "//&
                  trim(adjustl(fieldName))//&
                  " at message "//trim(adjustl(Msgs%name)))
          end if

          select case (idim_type)
          case (2)
             oneFieldSection =>  CreateFieldSection(&
                  field, &
                  fieldName, &
                  idim_type, &
                  xbComm(iNeigh)-x0, xeComm(iNeigh)-x0, &
                  ybComm(iNeigh)-y0, yeComm(iNeigh)-y0)
          case default
             write(c0,"(i8)") idim_type
             call fatal_error(h//" idim_type ("//trim(adjustl(c0))//&
                  ") incompatible with a 2D field")
          end select
          call AppendFieldSectionToMessageData(oneFieldSection, Msgs%msgData(cntMsg))
       end if
    end do
  end subroutine InsertFieldSectionAtMessageSet_2D





  subroutine InsertFieldSectionAtMessageSet_3D(&
       myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
       xbComm, xeComm, ybComm, yeComm, willComm, &
       Msgs)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable

    integer, intent(in) :: myNum

    ! field memory address

    real, pointer, intent(in) :: field(:,:,:)

    ! field name

    character(len=*), intent(in) :: fieldName

    ! idim_type codes the remained dimensions to communicate 

    integer, intent(in) :: idim_type

    ! nNeigh is number of processes for potential communication

    integer, intent(in) :: nNeigh

    ! Global indices of domain partition with Ghost Zones

    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! all remaining arguments are dimensioned by nNeigh

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
    character(len=*), parameter :: h="**(InsertFieldSectionAtMessageSet_3D)**"
    real, pointer :: PNull(:,:) => null()
    ! check arguments

    if (.not. associated(GlobalWithGhost)) then
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
    do iNeigh = 1, nNeigh
       if (willComm(iNeigh)) then
          cntMsg = cntMsg + 1
          if (cntMsg > nMsgs) then
             write(c0,"(i8)") nMsgs
             call fatal_error(h//" nMsgs ("//&
                  trim(adjustl(c0))//") exceeded while inserting field "//&
                  trim(adjustl(fieldName))//&
                  " at message "//trim(adjustl(Msgs%name)))
          end if

          select case (idim_type)
          case (3,6,7)
             oneFieldSection =>  CreateFieldSection(&
                  field, &
                  fieldName, &
                  idim_type, &
                  xbComm(iNeigh)-x0, xeComm(iNeigh)-x0, &
                  ybComm(iNeigh)-y0, yeComm(iNeigh)-y0)
          case default
             write(c0,"(i8)") idim_type
             call fatal_error(h//" idim_type ("//trim(adjustl(c0))//&
                  ") incompatible with a 2D field")
          end select
          call AppendFieldSectionToMessageData(oneFieldSection, Msgs%msgData(cntMsg))
       end if
    end do
  end subroutine InsertFieldSectionAtMessageSet_3D




  subroutine InsertFieldSectionAtMessageSet_4D(&
       myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
       xbComm, xeComm, ybComm, yeComm, willComm, &
       Msgs)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable

    integer, intent(in) :: myNum

    ! field memory address

    real, pointer, intent(in) :: field(:,:,:,:)

    ! field name

    character(len=*), intent(in) :: fieldName

    ! idim_type codes the remained dimensions to communicate 

    integer, intent(in) :: idim_type

    ! nNeigh is number of processes for potential communication

    integer, intent(in) :: nNeigh

    ! Global indices of domain partition with Ghost Zones

    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! all remaining arguments are dimensioned by nNeigh

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
    character(len=*), parameter :: h="**(InsertFieldSectionAtMessageSet_4D)**"
    real, pointer :: PNull(:,:) => null()
    ! check arguments

    if (.not. associated(GlobalWithGhost)) then
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
    do iNeigh = 1, nNeigh
       if (willComm(iNeigh)) then
          cntMsg = cntMsg + 1
          if (cntMsg > nMsgs) then
             write(c0,"(i8)") nMsgs
             call fatal_error(h//" nMsgs ("//&
                  trim(adjustl(c0))//") exceeded while inserting field "//&
                  trim(adjustl(fieldName))//&
                  " at message "//trim(adjustl(Msgs%name)))
          end if

          select case (idim_type)
          case (4,5)
             oneFieldSection =>  CreateFieldSection(&
                  field, &
                  fieldName, &
                  idim_type, &
                  xbComm(iNeigh)-x0, xeComm(iNeigh)-x0, &
                  ybComm(iNeigh)-y0, yeComm(iNeigh)-y0)
          case default
             write(c0,"(i8)") idim_type
             call fatal_error(h//" idim_type ("//trim(adjustl(c0))//&
                  ") incompatible with a 2D field")
          end select
          call AppendFieldSectionToMessageData(oneFieldSection, Msgs%msgData(cntMsg))
       end if
    end do
  end subroutine InsertFieldSectionAtMessageSet_4D



  subroutine InsertFieldSectionAtSendRecvMessageSetFromVTab(&
       varName, myNum, nNeigh, gridId, GlobalWithGhost, &
       xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    character(len=*), intent(in) :: varName

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable

    integer, intent(in) :: myNum

    ! nNeigh is the number of processes for potential communication

    integer, intent(in) :: nNeigh

    ! gridId is this grid number (required while vTabPtr is not included in type(Grid)

    integer, intent(in) :: gridId

    ! Global indices of domain partition with Ghost Zones

    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! this rank will send this rectangular region to each neighbour
    integer, intent(in) :: xbSend(nNeigh)
    integer, intent(in) :: xeSend(nNeigh)
    integer, intent(in) :: ybSend(nNeigh)
    integer, intent(in) :: yeSend(nNeigh)

    ! this rank will send messages to which neighbours
    logical, intent(in) :: willSend(nNeigh)

    ! send message set
    type(MessageSet), pointer, intent(inout) :: SendMessageSet

    ! this rank will recv messsages from which neighbours
    logical, intent(in) :: willRecv(nNeigh)

    ! this rank will recv this rectangular region from each neighbour
    integer, intent(in) :: xbRecv(nNeigh)
    integer, intent(in) :: xeRecv(nNeigh)
    integer, intent(in) :: ybRecv(nNeigh)
    integer, intent(in) :: yeRecv(nNeigh)

    ! recv message set
    type(MessageSet), pointer, intent(inout) :: RecvMessageSet


    type(var_tables_r), pointer :: vTabPtr => null()

    character(len=*), parameter :: h="**(InsertFieldSectionAtSendRecvMessageSetFromVTab)**"

    ! check arguments

    if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" GlobalWithGhost not associated")
    end if

    call GetVTabEntry(trim(varName), gridId, vTabPtr)

    ! include vTab field on field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, nNeigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, nNeigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet)
  end subroutine InsertFieldSectionAtSendRecvMessageSetFromVTab




  subroutine InsertFieldSectionAtSendRecvMessageSet_2D(&
       field, fieldName, idim_type, myNum, nNeigh, GlobalWithGhost, &
       xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! field memory address

    real, pointer, intent(in) :: field(:,:)

    ! field name

    character(len=*), intent(in) :: fieldName

    ! idim_type codes the remained dimensions to communicate 

    integer, intent(in) :: idim_type

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable

    integer, intent(in) :: myNum

    ! nNeigh is the number of processes for potential communication

    integer, intent(in) :: nNeigh

    ! Global indices of domain partition with Ghost Zones

    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! this rank will send this rectangular region to each neighbour
    integer, intent(in) :: xbSend(nNeigh)
    integer, intent(in) :: xeSend(nNeigh)
    integer, intent(in) :: ybSend(nNeigh)
    integer, intent(in) :: yeSend(nNeigh)

    ! this rank will send messages to which neighbours
    logical, intent(in) :: willSend(nNeigh)

    ! send message set
    type(MessageSet), pointer, intent(inout) :: SendMessageSet

    ! this rank will recv messsages from which neighbours
    logical, intent(in) :: willRecv(nNeigh)

    ! this rank will recv this rectangular region from each neighbour
    integer, intent(in) :: xbRecv(nNeigh)
    integer, intent(in) :: xeRecv(nNeigh)
    integer, intent(in) :: ybRecv(nNeigh)
    integer, intent(in) :: yeRecv(nNeigh)

    ! recv message set
    type(MessageSet), pointer, intent(inout) :: RecvMessageSet

    character(len=*), parameter :: h="**(InsertFieldSectionAtSendRecvMessageSet_2D)**"

    ! check arguments

    if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" GlobalWithGhost not associated")
    end if

    ! include field on field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet)
    call InsertFieldSectionAtMessageSet(&
         myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet)
  end subroutine InsertFieldSectionAtSendRecvMessageSet_2D





  subroutine InsertFieldSectionAtSendRecvMessageSet_3D(&
       field, fieldName, idim_type, myNum, nNeigh, GlobalWithGhost, &
       xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! field memory address

    real, pointer, intent(in) :: field(:,:,:)

    ! field name

    character(len=*), intent(in) :: fieldName

    ! idim_type codes the remained dimensions to communicate 

    integer, intent(in) :: idim_type

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable

    integer, intent(in) :: myNum

    ! nNeigh is the number of processes for potential communication

    integer, intent(in) :: nNeigh

    ! Global indices of domain partition with Ghost Zones

    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! this rank will send this rectangular region to each neighbour
    integer, intent(in) :: xbSend(nNeigh)
    integer, intent(in) :: xeSend(nNeigh)
    integer, intent(in) :: ybSend(nNeigh)
    integer, intent(in) :: yeSend(nNeigh)

    ! this rank will send messages to which neighbours
    logical, intent(in) :: willSend(nNeigh)

    ! send message set
    type(MessageSet), pointer, intent(inout) :: SendMessageSet

    ! this rank will recv messsages from which neighbours
    logical, intent(in) :: willRecv(nNeigh)

    ! this rank will recv this rectangular region from each neighbour
    integer, intent(in) :: xbRecv(nNeigh)
    integer, intent(in) :: xeRecv(nNeigh)
    integer, intent(in) :: ybRecv(nNeigh)
    integer, intent(in) :: yeRecv(nNeigh)

    ! recv message set
    type(MessageSet), pointer, intent(inout) :: RecvMessageSet

    character(len=*), parameter :: h="**(InsertFieldSectionAtSendRecvMessageSet_3D)**"

    ! check arguments

    if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" GlobalWithGhost not associated")
    end if

    ! include field on field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet)
    call InsertFieldSectionAtMessageSet(&
         myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet)
  end subroutine InsertFieldSectionAtSendRecvMessageSet_3D





  subroutine InsertFieldSectionAtSendRecvMessageSet_4D(&
       field, fieldName, idim_type, myNum, nNeigh, GlobalWithGhost, &
       xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! field memory address

    real, pointer, intent(in) :: field(:,:,:,:)

    ! field name

    character(len=*), intent(in) :: fieldName

    ! idim_type codes the remained dimensions to communicate 

    integer, intent(in) :: idim_type

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable

    integer, intent(in) :: myNum

    ! nNeigh is the number of processes for potential communication

    integer, intent(in) :: nNeigh

    ! Global indices of domain partition with Ghost Zones

    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost

    ! this rank will send this rectangular region to each neighbour
    integer, intent(in) :: xbSend(nNeigh)
    integer, intent(in) :: xeSend(nNeigh)
    integer, intent(in) :: ybSend(nNeigh)
    integer, intent(in) :: yeSend(nNeigh)

    ! this rank will send messages to which neighbours
    logical, intent(in) :: willSend(nNeigh)

    ! send message set
    type(MessageSet), pointer, intent(inout) :: SendMessageSet

    ! this rank will recv messsages from which neighbours
    logical, intent(in) :: willRecv(nNeigh)

    ! this rank will recv this rectangular region from each neighbour
    integer, intent(in) :: xbRecv(nNeigh)
    integer, intent(in) :: xeRecv(nNeigh)
    integer, intent(in) :: ybRecv(nNeigh)
    integer, intent(in) :: yeRecv(nNeigh)

    ! recv message set
    type(MessageSet), pointer, intent(inout) :: RecvMessageSet

    character(len=*), parameter :: h="**(InsertFieldSectionAtSendRecvMessageSet_4D)**"

    ! check arguments

    if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" GlobalWithGhost not associated")
    end if

    ! include field on field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet)
    call InsertFieldSectionAtMessageSet(&
         myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet)
  end subroutine InsertFieldSectionAtSendRecvMessageSet_4D





  subroutine CreateAcousticMessageSet(&
       gridId, GridSize, ParEnv, Neigh, &
       GlobalOwn, &
       GlobalOwnWithBC, &
       GlobalWithGhost, &
       AcouSendU, AcouRecvU, TagU, &
       AcouSendV, AcouRecvV, TagV, &
       AcouSendPNorth, AcouRecvPNorth, TagPNorth, &
       AcouSendPEast, AcouRecvPEast, TagPEast, &
       AcouSendUV, AcouRecvUV, TagUV, &
       AcouSendWP, AcouRecvWP, TagWP)

    integer, intent(in) :: gridId
    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: AcouSendU
    type(MessageSet), pointer, intent(inout) :: AcouRecvU
    type(MessageSet), pointer, intent(inout) :: AcouSendV
    type(MessageSet), pointer, intent(inout) :: AcouRecvV
    type(MessageSet), pointer, intent(inout) :: AcouSendPNorth
    type(MessageSet), pointer, intent(inout) :: AcouRecvPNorth
    type(MessageSet), pointer, intent(inout) :: AcouSendPEast
    type(MessageSet), pointer, intent(inout) :: AcouRecvPEast
    type(MessageSet), pointer, intent(inout) :: AcouSendUV
    type(MessageSet), pointer, intent(inout) :: AcouRecvUV
    type(MessageSet), pointer, intent(inout) :: AcouSendWP
    type(MessageSet), pointer, intent(inout) :: AcouRecvWP
    integer, intent(in) :: TagU
    integer, intent(in) :: TagV
    integer, intent(in) :: TagPNorth
    integer, intent(in) :: TagPEast
    integer, intent(in) :: TagUV
    integer, intent(in) :: TagWP

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    character(len=30) :: tmp_name

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(parEnv%nMachs)
    integer :: xeSend(parEnv%nMachs)
    integer :: ybSend(parEnv%nMachs)
    integer :: yeSend(parEnv%nMachs)
    integer :: xbRecv(parEnv%nMachs)
    integer :: xeRecv(parEnv%nMachs)
    integer :: ybRecv(parEnv%nMachs)
    integer :: yeRecv(parEnv%nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(parEnv%nMachs)
    logical :: willRecv(parEnv%nMachs)

    character(len=*), parameter :: NameSendU="AcouSendU"
    character(len=*), parameter :: NameRecvU="AcouRecvU"

    character(len=*), parameter :: NameSendV="AcouSendV"
    character(len=*), parameter :: NameRecvV="AcouRecvV"

    character(len=*), parameter :: NameSendPNorth="AcouSendPNorth"
    character(len=*), parameter :: NameRecvPNorth="AcouRecvPNorth"

    character(len=*), parameter :: NameSendPEast="AcouSendPEast"
    character(len=*), parameter :: NameRecvPEast="AcouRecvPEast"

    character(len=*), parameter :: NameSendUV="AcouSendUV"
    character(len=*), parameter :: NameRecvUV="AcouRecvUV"

    character(len=*), parameter :: NameSendWP="AcouSendWP"
    character(len=*), parameter :: NameRecvWP="AcouRecvWP"

    character(len=*), parameter :: h="**(CreateAcousticMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" starts with null GlobalOwn")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    if (.not. associated(Neigh)) then
       AcouSendU => null()
       AcouRecvU => null()
       AcouSendV => null()
       AcouRecvV => null()
       AcouSendPNorth => null()
       AcouRecvPNorth => null()
       AcouSendPEast => null()
       AcouRecvPEast => null()
       AcouSendUV => null()
       AcouRecvUV => null()
       AcouSendWP => null()
       AcouRecvWP => null()
       return
    end if

    if (dumpLocal) then
       call MsgDump(h//" starts creating AcousSend/RecvU")
    end if

    myNum  = ParEnv%myNum
    nMachs = ParEnv%nMachs
    nNeigh = Neigh%nNeigh

    ! AcouSendU, AcouRecvU:
    ! messages update GlobalOwn [xb-1:xb-1,yb:ye]

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagU, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwn, &
         NameSend=NameSendU, &
         NameRecv=NameRecvU, &
         xbToUpdate=GlobalOwn%xb-1, &
         xeToUpdate=GlobalOwn%xb-1, &
         ybToUpdate=GlobalOwn%yb, &
         yeToUpdate=GlobalOwn%ye, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=AcouSendU, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AcouRecvU)

    ! get field

    if(dyncore_flag==2) then
       tmp_name='UC'
    else
       tmp_name='UP'
    endif

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendU, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvU)

    if (dumpLocal) then
       call MsgDump(h//" finishes AcouSendU MessageSet:")
       call DumpMessageSet(AcouSendU)
       call MsgDump(h//" finishes AcouRecvU MessageSet:")
       call DumpMessageSet(AcouRecvU)
       call MsgDump(h//" starts creating AcousSend/RecvV")
    end if

    ! AcouSendV, AcouRecvV:
    ! messages update GlobalOwn [xb:xe,yb-1:yb-1]

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagV, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwn, &
         NameSend=NameSendV, &
         NameRecv=NameRecvV, &
         xbToUpdate=GlobalOwn%xb, &
         xeToUpdate=GlobalOwn%xe, &
         ybToUpdate=GlobalOwn%yb-1, &
         yeToUpdate=GlobalOwn%yb-1, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=AcouSendV, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AcouRecvV)

    ! get field

    if(dyncore_flag==2) then
       tmp_name='VC'
    else
       tmp_name='VP'
    endif

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendV, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvV)

    if (dumpLocal) then
       call MsgDump(h//" finishes AcouSendV MessageSet:")
       call DumpMessageSet(AcouSendV)
       call MsgDump(h//" finishes AcouRecvV MessageSet:")
       call DumpMessageSet(AcouRecvV)
       call MsgDump(h//" starts creating AcousSend/RecvPNorth")
    end if

    ! AcouSendPNorth, AcouRecvPNorth: union of
    !               GlobalOwn [xe+1:xe+1,yb:ye] with
    !               GlobalOwn [xb:xe,ye+1:ye+1]

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagPNorth, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwn, &
         NameSend=NameSendPNorth, &
         NameRecv=NameRecvPNorth, &
         xbToUpdate=GlobalOwn%xe+1, &
         xeToUpdate=GlobalOwn%xe+1, &
         ybToUpdate=GlobalOwn%yb, &
         yeToUpdate=GlobalOwn%ye, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=AcouSendPNorth, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AcouRecvPNorth)

    ! get field

    if(dyncore_flag==2) then
       tmp_name='PC'
    else
       tmp_name='PP'
    endif

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendPNorth, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvPNorth)

    if (dumpLocal) then
       call MsgDump(h//" finishes AcouSendPNorth MessageSet:")
       call DumpMessageSet(AcouSendPNorth)
       call MsgDump(h//" finishes AcouRecvPNorth MessageSet:")
       call DumpMessageSet(AcouRecvPNorth)
       call MsgDump(h//" starts creating AcousSend/RecvPEast")
    end if

    ! AcouSendPEast, AcouRecvPEast: 
    !               GlobalOwn [xb:xe,ye+1:ye+1]

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagPEast, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwn, &
         NameSend=NameSendPEast, &
         NameRecv=NameRecvPEast, &
         xbToUpdate=GlobalOwn%xb, &
         xeToUpdate=GlobalOwn%xe, &
         ybToUpdate=GlobalOwn%ye+1, &
         yeToUpdate=GlobalOwn%ye+1, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=AcouSendPEast, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AcouRecvPEast)

    ! get field

    if(dyncore_flag==2) then
       tmp_name='PC'
    else
       tmp_name='PP'
    endif

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendPEast, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvPEast)

    if (dumpLocal) then
       call MsgDump(h//" finishes AcouSendPEast MessageSet:")
       call DumpMessageSet(AcouSendPEast)
       call MsgDump(h//" finishes AcouRecvPEast MessageSet:")
       call DumpMessageSet(AcouRecvPEast)
       call MsgDump(h//" starts creating AcousSend/RecvUV")
    end if

    ! AcouSendUV, AcouRecvUV:
    ! messages update entire GhostZone

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagUV, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwnWithBC, &
         NameSend=NameSendUV, &
         NameRecv=NameRecvUV, &
         xbToUpdate=GlobalWithGhost%xb, &
         xeToUpdate=GlobalWithGhost%xe, &
         ybToUpdate=GlobalWithGhost%yb, &
         yeToUpdate=GlobalWithGhost%ye, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=AcouSendUV, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AcouRecvUV)

    ! get field UP

    if(dyncore_flag==2) then
       tmp_name='UC'
    else
       tmp_name='UP'
    endif

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendUV, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvUV)

    ! get field VP

    if(dyncore_flag==2) then
       tmp_name='VC'
    else
       tmp_name='VP'
    endif

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendUV, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvUV)

    if (dumpLocal) then
       call MsgDump(h//" finishes with AcouSendUV MessageSet:")
       call DumpMessageSet(AcouSendUV)
       call MsgDump(h//" finishes with AcouRecvUV MessageSet:")
       call DumpMessageSet(AcouRecvUV)
       call MsgDump(h//" starts creating AcousSend/RecvWP")
    end if

    ! AcouSendWP, AcouRecvWP:
    ! messages update entire GhostZone

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagWP, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwnWithBC, &
         NameSend=NameSendWP, &
         NameRecv=NameRecvWP, &
         xbToUpdate=GlobalWithGhost%xb, &
         xeToUpdate=GlobalWithGhost%xe, &
         ybToUpdate=GlobalWithGhost%yb, &
         yeToUpdate=GlobalWithGhost%ye, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=AcouSendWP, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AcouRecvWP)

    ! get field UP

    if(dyncore_flag==2) then
       tmp_name='WC'
    else
       tmp_name='WP'
    endif

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendWP, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvWP)

    ! get field VP

    if(dyncore_flag==2) then
       tmp_name='PC'
    else
       tmp_name='PP'
    endif

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendWP, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvWP)

    if (dumpLocal) then
       call MsgDump(h//" finishes AcouSendWP MessageSet:")
       call DumpMessageSet(AcouSendWP)
       call MsgDump(h//" finishes AcouRecvWP MessageSet:")
       call DumpMessageSet(AcouRecvWP)
    end if

  end subroutine CreateAcousticMessageSet





  subroutine DestroyAcousticMessageSet( &
       AcouSendU, AcouRecvU, &
       AcouSendV, AcouRecvV, &
       AcouSendPNorth, AcouRecvPNorth, &
       AcouSendPEast, AcouRecvPEast, &
       AcouSendUV, AcouRecvUV, &
       AcouSendWP, AcouRecvWP)

    type(MessageSet), pointer, intent(inout) :: AcouSendU
    type(MessageSet), pointer, intent(inout) :: AcouRecvU
    type(MessageSet), pointer, intent(inout) :: AcouSendV
    type(MessageSet), pointer, intent(inout) :: AcouRecvV
    type(MessageSet), pointer, intent(inout) :: AcouSendPNorth
    type(MessageSet), pointer, intent(inout) :: AcouRecvPNorth
    type(MessageSet), pointer, intent(inout) :: AcouSendPEast
    type(MessageSet), pointer, intent(inout) :: AcouRecvPEast
    type(MessageSet), pointer, intent(inout) :: AcouSendUV
    type(MessageSet), pointer, intent(inout) :: AcouRecvUV
    type(MessageSet), pointer, intent(inout) :: AcouSendWP
    type(MessageSet), pointer, intent(inout) :: AcouRecvWP

    character(len=*), parameter :: h="**(DestroyAcousticMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" will destroy "//&
            " AcouSend/RecvU, AcouSend/RecvV, AcouSend/RecvP,"//&
            " AcouSend/RecvUV and AcouSend/RecvWP")
    end if

    call DestroyMessageSet(AcouSendU)
    call DestroyMessageSet(AcouRecvU)

    call DestroyMessageSet(AcouSendV)
    call DestroyMessageSet(AcouRecvV)

    call DestroyMessageSet(AcouSendPNorth)
    call DestroyMessageSet(AcouRecvPNorth)

    call DestroyMessageSet(AcouSendPEast)
    call DestroyMessageSet(AcouRecvPEast)

    call DestroyMessageSet(AcouSendUV)
    call DestroyMessageSet(AcouRecvUV)

    call DestroyMessageSet(AcouSendWP)
    call DestroyMessageSet(AcouRecvWP)
  end subroutine DestroyAcousticMessageSet





  subroutine CreateDn0MessageSet(&
       gridId, GridSize, ParEnv, Neigh, &
       GlobalOwn, GlobalWithGhost, &
       SendDn0u, RecvDn0u, TagDn0u, &
       SendDn0v, RecvDn0v, TagDn0v)

    integer, intent(in) :: gridId
    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: SendDn0u
    type(MessageSet), pointer, intent(inout) :: RecvDn0u
    type(MessageSet), pointer, intent(inout) :: SendDn0v
    type(MessageSet), pointer, intent(inout) :: RecvDn0v
    integer, intent(in) :: TagDn0u
    integer, intent(in) :: TagDn0v

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(parEnv%nMachs)
    integer :: xeSend(parEnv%nMachs)
    integer :: ybSend(parEnv%nMachs)
    integer :: yeSend(parEnv%nMachs)
    integer :: xbRecv(parEnv%nMachs)
    integer :: xeRecv(parEnv%nMachs)
    integer :: ybRecv(parEnv%nMachs)
    integer :: yeRecv(parEnv%nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(parEnv%nMachs)
    logical :: willRecv(parEnv%nMachs)

    character(len=*), parameter :: NameSendDn0u="SendDn0u"
    character(len=*), parameter :: NameRecvDn0u="RecvDn0u"

    character(len=*), parameter :: NameSendDn0v="SendDn0v"
    character(len=*), parameter :: NameRecvDn0v="RecvDn0v"

    character(len=30) :: tmp_name
    character(len=*), parameter :: h="**(CreateDn0MessageSet)**"
    logical, parameter :: dumpLocal=.false.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" starts with null GlobalOwn")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    if (.not. associated(Neigh)) then

       ! default output (case no neighbours)

       SendDn0u => null()
       RecvDn0u => null()
       SendDn0v => null()
       RecvDn0v => null()
       return
    end if

    if (dumpLocal) then
       call MsgDump(h//" will create Send/RecvDn0u and Send/RecvDn0v")
    end if

    myNum  = ParEnv%myNum
    nNeigh = Neigh%nNeigh

    ! SendDn0u, RecvDn0u:
    ! messages update GlobalOwn [xe+1:xe+1,yb:ye]

    call NodesRegionsSendRecv(&
         nMachs=ParEnv%nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagDn0u, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwn, &
         NameSend=NameSendDn0u, &
         NameRecv=NameRecvDn0u, &
         xbToUpdate=GlobalOwn%xe + 1, &
         xeToUpdate=GlobalOwn%xe + 1, &
         ybToUpdate=GlobalOwn%yb, &
         yeToUpdate=GlobalOwn%ye, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=SendDn0u, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=RecvDn0u)

    ! get field

    tmp_name='DN0U'

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendDn0u, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvDn0u)

    if (dumpLocal) then
       call MsgDump(h//" done building SendDn0u MessageSet:")
       call DumpMessageSet(SendDn0u)
       call MsgDump(h//" done building RecvDn0u MessageSet:")
       call DumpMessageSet(RecvDn0u)
    end if

    ! SendDn0v, RecvDn0v:
    ! messages update GlobalOwn [xb:xe,ye+1:ye+1]

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagDn0v, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwn, &
         NameSend=NameSendDn0v, &
         NameRecv=NameRecvDn0v, &
         xbToUpdate=GlobalOwn%xb, &
         xeToUpdate=GlobalOwn%xe, &
         ybToUpdate=GlobalOwn%ye+1, &
         yeToUpdate=GlobalOwn%ye+1, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=SendDn0v, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=RecvDn0v)

    ! get field

    tmp_name='DN0V'

    call InsertFieldSectionAtSendRecvMessageSet(&
         tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendDn0v, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvDn0v)

    if (dumpLocal) then
       call MsgDump(h//" done building SendDn0v MessageSet:")
       call DumpMessageSet(SendDn0v)
       call MsgDump(h//" done building RecvDn0v MessageSet:")
       call DumpMessageSet(RecvDn0v)
    end if
  end subroutine CreateDn0MessageSet





  subroutine DestroyDn0MessageSet( &
       SendDn0u, RecvDn0u, SendDn0v, RecvDn0v)

    type(MessageSet), pointer, intent(inout) :: SendDn0u
    type(MessageSet), pointer, intent(inout) :: RecvDn0u
    type(MessageSet), pointer, intent(inout) :: SendDn0v
    type(MessageSet), pointer, intent(inout) :: RecvDn0v
    character(len=*), parameter :: h="**(DestroyDn0MessageSet)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" will destroy "//&
            " Send/RecvDn0u and Send/RecvDn0v")
    end if

    call DestroyMessageSet(SendDn0u)
    call DestroyMessageSet(RecvDn0u)
    call DestroyMessageSet(SendDn0v)
    call DestroyMessageSet(RecvDn0v)

  end subroutine DestroyDn0MessageSet





  subroutine CreateG3DMessageSet(&
       gridId, GridSize, ParEnv, Neigh, &
       GlobalOwnWithBC, GlobalWithGhost, &
       Ramsin, &
       SendG3D, RecvG3D, TagG3D)

    integer, intent(in) :: gridId
    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(NamelistFile), pointer, intent(in) :: Ramsin
    type(MessageSet), pointer, intent(inout) :: SendG3D
    type(MessageSet), pointer, intent(inout) :: RecvG3D
    integer, intent(in) :: TagG3D

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    integer :: g3d_spread
    integer :: g3d_smoothh

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(parEnv%nMachs)
    integer :: xeSend(parEnv%nMachs)
    integer :: ybSend(parEnv%nMachs)
    integer :: yeSend(parEnv%nMachs)
    integer :: xbRecv(parEnv%nMachs)
    integer :: xeRecv(parEnv%nMachs)
    integer :: ybRecv(parEnv%nMachs)
    integer :: yeRecv(parEnv%nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(parEnv%nMachs)
    logical :: willRecv(parEnv%nMachs)

    character(len=*), parameter :: NameSendG3D="SendG3D"
    character(len=*), parameter :: NameRecvG3D="RecvG3D"

    integer :: vTabNbr
    character(len=30) :: tmp_name

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateG3DMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwnWithBC)) then
       call fatal_error(h//" starts with null GlobalOwnWithBC")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours or
    ! selected namelist variables were not set)

    SendG3D => null()
    RecvG3D => null()

    ! there will be messages if there are neighbour nodes
    ! and any of the namelist variables
    ! g3d_spread or g3d_smoothh were set

    g3d_spread = Ramsin%g3d_spread
    g3d_smoothh = Ramsin%g3d_smoothh

    if (associated(Neigh) .and. &
         (g3d_spread == 1 .or. g3d_smoothh == 1)) then

       if (dumpLocal) then
          call MsgDump(h//" will create Send/RecvG3D")
       end if

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       ! SendG3D, RecvG3D:
       ! messages update entire GhostZone

       call NodesRegionsSendRecv(&
            nMachs=nMachs, &
            nNeigh=nNeigh, &
            myNum=myNum, &
            tag=TagG3D, &
            Neigh=Neigh, &
            GlobalOwn=GlobalOwnWithBC, &
            NameSend=NameSendG3D, &
            NameRecv=NameRecvG3d, &
            xbToUpdate=GlobalOwnWithBC%xb, &
            xeToUpdate=GlobalOwnWithBC%xe, &
            ybToUpdate=GlobalOwnWithBC%yb, &
            yeToUpdate=GlobalOwnWithBC%ye, &
            xbSend=xbSend, &
            xeSend=xeSend, &
            ybSend=ybSend, &
            yeSend=yeSend, &
            willSend=willSend, &
            SendMessageSet=SendG3D, &
            xbRecv=xbRecv, &
            xeRecv=xeRecv, &
            ybRecv=ybRecv, &
            yeRecv=yeRecv, &
            willRecv=willRecv, &
            RecvMessageSet=RecvG3D)

       ! when g3d_spread is selected, send and receive fields TTENS and QVTTENS

       if (g3d_spread == 1) then

          tmp_name='TTENS'

          call InsertFieldSectionAtSendRecvMessageSet(&
               tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, SendG3D, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvG3D)

          tmp_name='QVTTENS'

          call InsertFieldSectionAtSendRecvMessageSet(&
               tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, SendG3D, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvG3D)
       end if

       ! when g3d_smoothh is selected, send and receive fields THSRC and RTSRC

       if (g3d_smoothh == 1) then

          tmp_name='THSRC'

          call InsertFieldSectionAtSendRecvMessageSet(&
               tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, SendG3D, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvG3D)

          tmp_name='RTSRC'

          call InsertFieldSectionAtSendRecvMessageSet(&
               tmp_name, myNum, nNeigh, gridId, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, SendG3D, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvG3D)

       end if
       if (dumpLocal) then
          call MsgDump(h//" done building SendG3D MessageSet:")
          call DumpMessageSet(SendG3D)
          call MsgDump(h//" done building RecvG3D MessageSet:")
          call DumpMessageSet(RecvG3D)
       end if
    end if
  end subroutine CreateG3DMessageSet





  subroutine DestroyG3DMessageSet( &
       SendG3D, RecvG3D)

    type(MessageSet), pointer, intent(inout) :: SendG3D
    type(MessageSet), pointer, intent(inout) :: RecvG3D
    character(len=*), parameter :: h="**(DestroyG3DMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" will destroy Send/RecvG3D")
    end if

    call DestroyMessageSet(SendG3D)
    call DestroyMessageSet(RecvG3D)

  end subroutine DestroyG3DMessageSet





  subroutine CreateSelectedGhostZoneMessageSet(&
       gridId, num_var, vtab_r, &
       GridSize, ParEnv, Neigh, &
       GlobalOwnWithBC, GlobalWithGhost, &
       SelectedGhostZoneSend, SelectedGhostZoneRecv, TagSelectedGhostZone)

    integer, intent(in) :: gridId
    integer, intent(in) :: num_var(:)
    type(var_tables_r), target, intent(in) ::  vtab_r(:,:)
    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: SelectedGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: SelectedGhostZoneRecv
    integer, intent(in) :: TagSelectedGhostZone


    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    integer :: vTabNbr

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(parEnv%nMachs)
    integer :: xeSend(parEnv%nMachs)
    integer :: ybSend(parEnv%nMachs)
    integer :: yeSend(parEnv%nMachs)
    integer :: xbRecv(parEnv%nMachs)
    integer :: xeRecv(parEnv%nMachs)
    integer :: ybRecv(parEnv%nMachs)
    integer :: yeRecv(parEnv%nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(parEnv%nMachs)
    logical :: willRecv(parEnv%nMachs)

    character(len=*), parameter :: NameSendSelectedGhostZone="SelectedGhostZoneSend"
    character(len=*), parameter :: NameRecvSelectedGhostZone="SelectedGhostZoneRecv"

    character(len=*), parameter :: h="**(CreateSelectedGhostZoneMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwnWithBC)) then
       call fatal_error(h//" starts with null GlobalOwnWithBC")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    if (.not. associated(Neigh)) then
       SelectedGhostZoneSend => null()
       SelectedGhostZoneRecv => null()
       return
    end if

    if (dumpLocal) then
       call MsgDump(h//" will create SelectedGhostZoneSend/Recv")
    end if

    myNum  = ParEnv%myNum
    nMachs = ParEnv%nMachs
    nNeigh = Neigh%nNeigh

    ! SelectedGhostZoneSend, SelectedGhostZoneRecv:
    ! messages update entire GhostZone

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagSelectedGhostZone, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwnWithBC, &
         NameSend=NameSendSelectedGhostZone, &
         NameRecv=NameRecvSelectedGhostZone, &
         xbToUpdate=GlobalWithGhost%xb, &
         xeToUpdate=GlobalWithGhost%xe, &
         ybToUpdate=GlobalWithGhost%yb, &
         yeToUpdate=GlobalWithGhost%ye, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=SelectedGhostZoneSend, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=SelectedGhostZoneRecv)

    ! take all var_tables field that should be communicated

    do vTabNbr = 1, num_var(gridId)

       if (vtab_r(vTabNbr,gridId)%impt1 == 1) then

          call InsertFieldSectionAtSendRecvMessageSet(&
               vtab_r(vTabNbr,gridId)%name, myNum, nNeigh, gridId, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, SelectedGhostZoneSend, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, SelectedGhostZoneRecv)
       end if
    end do

    if (dumpLocal) then
       call MsgDump(h//" finishes with SelectedGhostZoneSend MessageSet:")
       call DumpMessageSet(SelectedGhostZoneSend)
       call MsgDump(h//" finishes with SelectedGhostZoneRecv MessageSet:")
       call DumpMessageSet(SelectedGhostZoneRecv)
    end if
  end subroutine CreateSelectedGhostZoneMessageSet





  subroutine DestroySelectedGhostZoneMessageSet( &
       SelectedGhostZoneSend, SelectedGhostZoneRecv)

    type(MessageSet), pointer, intent(inout) :: SelectedGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: SelectedGhostZoneRecv
    character(len=*), parameter :: h="**(DestroySelectedGhostZoneMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" will destroy SelectedGhostZoneSend/Recv")
    end if
    call DestroyMessageSet(SelectedGhostZoneSend)
    call DestroyMessageSet(SelectedGhostZoneRecv)

  end subroutine DestroySelectedGhostZoneMessageSet





  subroutine CreateAllGhostZoneMessageSet(&
       gridId, num_var, vtab_r, &
       GridSize, ParEnv, Neigh, &
       GlobalOwnWithBC, GlobalWithGhost, &
       AllGhostZoneSend, AllGhostZoneRecv, TagAllGhostZone)

    integer, intent(in) :: gridId
    integer, intent(in) :: num_var(:)
    type(var_tables_r), target, intent(in) ::  vtab_r(:,:)
    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: AllGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: AllGhostZoneRecv
    integer, intent(in) :: TagAllGhostZone


    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    integer :: vTabNbr
    character(len=32) :: vTabName

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(parEnv%nMachs)
    integer :: xeSend(parEnv%nMachs)
    integer :: ybSend(parEnv%nMachs)
    integer :: yeSend(parEnv%nMachs)
    integer :: xbRecv(parEnv%nMachs)
    integer :: xeRecv(parEnv%nMachs)
    integer :: ybRecv(parEnv%nMachs)
    integer :: yeRecv(parEnv%nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(parEnv%nMachs)
    logical :: willRecv(parEnv%nMachs)

    character(len=*), parameter :: NameSendAllGhostZone="AllGhostZoneSend"
    character(len=*), parameter :: NameRecvAllGhostZone="AllGhostZoneRecv"

    character(len=*), parameter :: h="**(CreateAllGhostZoneMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwnWithBC)) then
       call fatal_error(h//" starts with null GlobalOwnWithBC")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    if (.not. associated(Neigh)) then
       AllGhostZoneSend => null()
       AllGhostZoneRecv => null()
       return
    end if


    if (dumpLocal) then
       call MsgDump(h//" will create AllGhostZoneSend/Recv")
    end if

    myNum  = ParEnv%myNum
    nMachs = ParEnv%nMachs
    nNeigh = Neigh%nNeigh

    ! AllGhostZoneSend, AllGhostZoneRecv:
    ! messages update entire GhostZone

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=TagAllGhostZone, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwnWithBC, &
         NameSend=NameSendAllGhostZone, &
         NameRecv=NameRecvAllGhostZone, &
         xbToUpdate=GlobalWithGhost%xb, &
         xeToUpdate=GlobalWithGhost%xe, &
         ybToUpdate=GlobalWithGhost%yb, &
         yeToUpdate=GlobalWithGhost%ye, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=AllGhostZoneSend, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AllGhostZoneRecv)

    ! take all var_tables field that should be communicated

    do vTabNbr = 1, num_var(gridId)

       vTabName = vtab_r(vTabNbr,gridId)%name

       if (&
            trim(adjustl(vTabName)) /= "LPU" .and. &
            trim(adjustl(vTabName)) /= "LPV" .and. &
            trim(adjustl(vTabName)) /= "LPW" ) then

          call InsertFieldSectionAtSendRecvMessageSet(&
               vTabName, myNum, nNeigh, gridId, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, AllGhostZoneSend, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AllGhostZoneRecv)
       end if
    end do

    if (dumpLocal) then
       call MsgDump(h//" finishes with AllGhostZoneSend MessageSet:")
       call DumpMessageSet(AllGhostZoneSend)
       call MsgDump(h//" finishes with AllGhostZoneRecv MessageSet:")
       call DumpMessageSet(AllGhostZoneRecv)
    end if
  end subroutine CreateAllGhostZoneMessageSet





  subroutine DestroyAllGhostZoneMessageSet( &
       AllGhostZoneSend, AllGhostZoneRecv)

    type(MessageSet), pointer, intent(inout) :: AllGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: AllGhostZoneRecv
    character(len=*), parameter :: h="**(DestroyAllGhostZoneMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" will destroy AllGhostZoneSend/Recv")
    end if

    call DestroyMessageSet(AllGhostZoneSend)
    call DestroyMessageSet(AllGhostZoneRecv)

  end subroutine DestroyAllGhostZoneMessageSet




  subroutine CreateAcouDampOneMessageSet(&
       field, fieldName, idim_type,  &
       ParEnv, Neigh, GlobalOwn, GlobalWithGhost, &
       Tag, NameSend, NameRecv, &
       AcouDampOneSend, AcouDampOneRecv)

    real, pointer, intent(in) :: field(:,:,:)
    character(len=*), intent(in) :: fieldName
    integer, intent(in) :: idim_type
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    integer, intent(in) :: Tag
    character(len=*), intent(in) :: NameSend
    character(len=*), intent(in) :: NameRecv
    type(MessageSet), pointer, intent(inout) :: AcouDampOneSend
    type(MessageSet), pointer, intent(inout) :: AcouDampOneRecv


    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    integer :: vTabNbr
    character(len=32) :: vTabName

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(parEnv%nMachs)
    integer :: xeSend(parEnv%nMachs)
    integer :: ybSend(parEnv%nMachs)
    integer :: yeSend(parEnv%nMachs)
    integer :: xbRecv(parEnv%nMachs)
    integer :: xeRecv(parEnv%nMachs)
    integer :: ybRecv(parEnv%nMachs)
    integer :: yeRecv(parEnv%nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(parEnv%nMachs)
    logical :: willRecv(parEnv%nMachs)

    character(len=*), parameter :: h="**(CreateAcouDampOneMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" starts with null GlobalOwn")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    if (.not. associated(Neigh)) then
       AcouDampOneSend => null()
       AcouDampOneRecv => null()
       return
    end if

    if (dumpLocal) then
       call MsgDump(h//" will create AcouDampOneSend/Recv")
    end if

    myNum  = ParEnv%myNum
    nMachs = ParEnv%nMachs
    nNeigh = Neigh%nNeigh

    ! AcouDampOneSend, AcouDampOneRecv:
    ! messages update entire GhostZone

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
         nNeigh=nNeigh, &
         myNum=myNum, &
         tag=Tag, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwn, &
         NameSend=NameSend, &
         NameRecv=NameRecv, &
         xbToUpdate=GlobalWithGhost%xb, &
         xeToUpdate=GlobalWithGhost%xe, &
         ybToUpdate=GlobalWithGhost%yb, &
         yeToUpdate=GlobalWithGhost%ye, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         SendMessageSet=AcouDampOneSend, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AcouDampOneRecv)

    !

    call InsertFieldSectionAtSendRecvMessageSet(&
         field, fieldName, idim_type, myNum, nNeigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouDampOneSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouDampOneRecv)

    if (dumpLocal) then
       call MsgDump(h//" finishes with AcouDampOneSend MessageSet:")
       call DumpMessageSet(AcouDampOneSend)
       call MsgDump(h//" finishes with AcouDampOneRecv MessageSet:")
       call DumpMessageSet(AcouDampOneRecv)
    end if
  end subroutine CreateAcouDampOneMessageSet





  subroutine DestroyAcouDampOneMessageSet( &
       AcouDampOneSend, AcouDampOneRecv)

    type(MessageSet), pointer, intent(inout) :: AcouDampOneSend
    type(MessageSet), pointer, intent(inout) :: AcouDampOneRecv
    character(len=*), parameter :: h="**(DestroyAcouDampOneMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" will destroy AcouDampOneSend/Recv")
    end if

    call DestroyMessageSet(AcouDampOneSend)
    call DestroyMessageSet(AcouDampOneRecv)

  end subroutine DestroyAcouDampOneMessageSet




  
  subroutine CreateWideGhostZoneMessageSet(&
       GridSize, ParEnv, Neigh, &
       GlobalOwnWithBC, GlobalWithGhost, LocalOwn, &
       TagWideGhostZone, WideGhostZoneSend, WideGhostZoneRecv)

    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(DomainDecomp), pointer, intent(in) :: LocalOwn
    integer, intent(in) :: TagWideGhostZone
    type(MessageSet), pointer, intent(inout) :: WideGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: WideGhostZoneRecv

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions to send and receive

    integer :: xbSendNorth(parEnv%nMachs)
    integer :: xeSendNorth(parEnv%nMachs)
    integer :: ybSendNorth(parEnv%nMachs)
    integer :: yeSendNorth(parEnv%nMachs)
    integer :: xbRecvNorth(parEnv%nMachs)
    integer :: xeRecvNorth(parEnv%nMachs)
    integer :: ybRecvNorth(parEnv%nMachs)
    integer :: yeRecvNorth(parEnv%nMachs)

    integer :: xbSendSouth(parEnv%nMachs)
    integer :: xeSendSouth(parEnv%nMachs)
    integer :: ybSendSouth(parEnv%nMachs)
    integer :: yeSendSouth(parEnv%nMachs)
    integer :: xbRecvSouth(parEnv%nMachs)
    integer :: xeRecvSouth(parEnv%nMachs)
    integer :: ybRecvSouth(parEnv%nMachs)
    integer :: yeRecvSouth(parEnv%nMachs)

    integer :: xbSendEast(parEnv%nMachs)
    integer :: xeSendEast(parEnv%nMachs)
    integer :: ybSendEast(parEnv%nMachs)
    integer :: yeSendEast(parEnv%nMachs)
    integer :: xbRecvEast(parEnv%nMachs)
    integer :: xeRecvEast(parEnv%nMachs)
    integer :: ybRecvEast(parEnv%nMachs)
    integer :: yeRecvEast(parEnv%nMachs)

    integer :: xbSendWest(parEnv%nMachs)
    integer :: xeSendWest(parEnv%nMachs)
    integer :: ybSendWest(parEnv%nMachs)
    integer :: yeSendWest(parEnv%nMachs)
    integer :: xbRecvWest(parEnv%nMachs)
    integer :: xeRecvWest(parEnv%nMachs)
    integer :: ybRecvWest(parEnv%nMachs)
    integer :: yeRecvWest(parEnv%nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSendNorth(parEnv%nMachs)
    logical :: willRecvNorth(parEnv%nMachs)

    logical :: willSendSouth(parEnv%nMachs)
    logical :: willRecvSouth(parEnv%nMachs)

    logical :: willSendEast(parEnv%nMachs)
    logical :: willRecvEast(parEnv%nMachs)

    logical :: willSendWest(parEnv%nMachs)
    logical :: willRecvWest(parEnv%nMachs)

    logical :: willSend(parEnv%nMachs)
    logical :: willRecv(parEnv%nMachs)

    character(len=*), parameter :: NameSend="SendWideGhostZone"
    character(len=*), parameter :: NameRecv="RecvWideGhostZone"
    character(len=*), parameter :: NameSendRecvNorth="Send/RecvWideGhostZoneNorth"
    character(len=*), parameter :: NameSendRecvSouth="Send/RecvWideGhostZoneSouth"
    character(len=*), parameter :: NameSendRecvEast="Send/RecvWideGhostZoneEast"
    character(len=*), parameter :: NameSendRecvWest="Send/RecvWideGhostZoneWest"

    character(len=*), parameter :: scrName="SCR"
    character(len=*), parameter :: ufxName="UFX"
    character(len=*), parameter :: vfxName="VFX"
    character(len=*), parameter :: wfxName="WFX"

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    integer :: ierr
    integer :: iNeigh
    integer :: iNode
    integer :: x0
    integer :: y0
    integer :: nMsgs
    integer :: cntMsg
    integer :: fieldSectionSize
    integer, parameter :: ghostZoneWidth=3
    integer, parameter :: idim_type=3
    type(FieldSection), pointer :: oneFieldSection

    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateWideGhostZoneMessageSet)**"


    nMachs=ParEnv%nmachs
    myNum=ParEnv%mynum
    nNeigh=Neigh%nNeigh

    if (dumpLocal) then
       write(str(1),"(i8)") nMachs
       write(str(2),"(i8)") myNum
       call MsgDump(h//" enter with nMachs="//trim(adjustl(str(1)))//&
            "; myNum="//trim(adjustl(str(2))))
       call DumpNeighbourNodes(Neigh,"WideGhostZone")
    end if

    ! offsets to convert global indices to local indices at this proc

    x0 = GlobalWithGhost%xb(myNum) - 1
    y0 = GlobalWithGhost%yb(myNum) - 1

    ! neighbour communication

    ! which neighbour nodes will send and receive, given by
    ! willSend and willRecv

    ! which intervals will be send, given by
    ! (xbSend:xeSend,ybSend:yeSend)

    !  which intervals will be received, given by
    ! (xbRecv:xeRecv,ybRecv:yeRecv)

    ! intervals are computed as the intersection of
    ! (xbToUpdate:xeToUpdate, ybToUpdate:yeToUpdate)
    ! with GlobalOwnWithBC of neighbour ranks (given by Neigh)

    ! north neighbour communication

    if (dumpLocal) then
       call MsgDump(h//" compute NodesToSendRecvMessages to update North Ghost Zone")
    end if

    call NodesToSendRecvMessages( &
         thisNode=myNum, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwnWithBC, &
         xbToUpdate=GlobalOwnWithBC%xb, &
         xeToUpdate=GlobalOwnWithBC%xe, &
         ybToUpdate=GlobalOwnWithBC%ye+1, &
         yeToUpdate=GlobalOwnWithBC%ye+ghostZoneWidth, &
         xbSend=xbSendNorth, &
         xeSend=xeSendNorth, &
         ybSend=ybSendNorth, &
         yeSend=yeSendNorth, &
         willSend=willSendNorth, &
         xbRecv=xbRecvNorth, &
         xeRecv=xeRecvNorth, &
         ybRecv=ybRecvNorth, &
         yeRecv=yeRecvNorth, &
         willRecv=willRecvNorth, &
         varName=NameSendRecvNorth)

    ! south neighbour communication

    if (dumpLocal) then
       call MsgDump(h//" compute NodesToSendRecvMessages to update South Ghost Zone")
    end if

    call NodesToSendRecvMessages( &
         thisNode=myNum, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwnWithBC, &
         xbToUpdate=GlobalOwnWithBC%xb, &
         xeToUpdate=GlobalOwnWithBC%xe, &
         ybToUpdate=GlobalOwnWithBC%yb-ghostZoneWidth, &
         yeToUpdate=GlobalOwnWithBC%yb-1, &
         xbSend=xbSendSouth, &
         xeSend=xeSendSouth, &
         ybSend=ybSendSouth, &
         yeSend=yeSendSouth, &
         willSend=willSendSouth, &
         xbRecv=xbRecvSouth, &
         xeRecv=xeRecvSouth, &
         ybRecv=ybRecvSouth, &
         yeRecv=yeRecvSouth, &
         willRecv=willRecvSouth, &
         varName=NameSendRecvSouth)

    ! east neighbour communication

    if (dumpLocal) then
       call MsgDump(h//" compute NodesToSendRecvMessages to update East Ghost Zone")
    end if

    call NodesToSendRecvMessages( &
         thisNode=myNum, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwnWithBC, &
         xbToUpdate=GlobalOwnWithBC%xe+1, &
         xeToUpdate=GlobalOwnWithBC%xe+ghostZoneWidth, &
         ybToUpdate=GlobalOwnWithBC%yb, &
         yeToUpdate=GlobalOwnWithBC%ye, &
         xbSend=xbSendEast, &
         xeSend=xeSendEast, &
         ybSend=ybSendEast, &
         yeSend=yeSendEast, &
         willSend=willSendEast, &
         xbRecv=xbRecvEast, &
         xeRecv=xeRecvEast, &
         ybRecv=ybRecvEast, &
         yeRecv=yeRecvEast, &
         willRecv=willRecvEast, &
         varName=NameSendRecvEast)

    ! west neighbour communication

    if (dumpLocal) then
       call MsgDump(h//" compute NodesToSendRecvMessages to update West Ghost Zone")
    end if

    call NodesToSendRecvMessages( &
         thisNode=myNum, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwnWithBC, &
         xbToUpdate=GlobalOwnWithBC%xb-ghostZoneWidth, &
         xeToUpdate=GlobalOwnWithBC%xb-1, &
         ybToUpdate=GlobalOwnWithBC%yb, &
         yeToUpdate=GlobalOwnWithBC%ye, &
         xbSend=xbSendWest, &
         xeSend=xeSendWest, &
         ybSend=ybSendWest, &
         yeSend=yeSendWest, &
         willSend=willSendWest, &
         xbRecv=xbRecvWest, &
         xeRecv=xeRecvWest, &
         ybRecv=ybRecvWest, &
         yeRecv=yeRecvWest, &
         willRecv=willRecvWest, &
         varName=NameSendRecvWest)

    ! send message set will contain sends for all four directions

    willSend = willSendNorth .or. willSendSouth .or. willSendEast .or. willSendWest
    if (dumpLocal) then
       do iNeigh = 1, nNeigh
          write(str(1),"(l)") willSendNorth(iNeigh)
          write(str(2),"(l)") willSendSouth(iNeigh)
          write(str(3),"(l)") willSendEast(iNeigh)
          write(str(4),"(l)") willSendWest(iNeigh)
          write(str(5),"(i8)") Brams2MpiProcNbr(Neigh%neigh(iNeigh))
          call MsgDump(h//" send to MPI proc #"//trim(adjustl(str(5)))//&
               " willSendNorth="//trim(adjustl(str(1)))//&
               " willSendSouth="//trim(adjustl(str(2)))//&
               " willSendEast="//trim(adjustl(str(3)))//&
               " willSendWest="//trim(adjustl(str(4))))
       end do
    end if

    ! create message set for all sends

    WideGhostZoneSend => CreateMessageSet(&
         NameSend, &
         sendDirection, &
         TagWideGhostZone, &
         willSend, &
         Neigh)

    ! insert field sections named SCR, UFX, VFX, WFX at each direction
    ! with null field addresses, to be updated whenever
    ! real addresses are known

    if ( associated(WideGhostZoneSend)) then
       nMsgs = WideGhostZoneSend%nMsgs

       ! create list of Field Sections to send, one for
       ! each process to communicate and insert at the send MessageSet
       ! field section list

       ! since there is at most one neighbour node at each direction,
       ! there will be at most one MessageSet at each direction

       cntMsg = 0
       do iNeigh = 1, nNeigh

          if (willSendNorth(iNeigh)) then

             ! insert send communications to north

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting four fields "//&
                     " at message "//trim(adjustl(WideGhostZoneSend%name)))
             end if

             fieldSectionSize= &
                  (xeSendNorth(iNeigh)-xbSendNorth(iNeigh)+1) * &
                  (yeSendNorth(iNeigh)-ybSendNorth(iNeigh)+1) * &
                  (GridSize%nnzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                  ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                  ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                  ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                  ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))

          else if (willSendSouth(iNeigh)) then

             ! insert send communications to south

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting four fields "//&
                     " at message "//trim(adjustl(WideGhostZoneSend%name)))
             end if

             fieldSectionSize= &
                  (xeSendSouth(iNeigh)-xbSendSouth(iNeigh)+1) * &
                  (yeSendSouth(iNeigh)-ybSendSouth(iNeigh)+1) * &
                  (GridSize%nnzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                  ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                  ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                  ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                  ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))

          else if (willSendEast(iNeigh)) then

             ! insert send communications to east

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting four fields "//&
                     " at message "//trim(adjustl(WideGhostZoneSend%name)))
             end if

             fieldSectionSize= &
                  (xeSendEast(iNeigh)-xbSendEast(iNeigh)+1) * &
                  (yeSendEast(iNeigh)-ybSendEast(iNeigh)+1) * &
                  (GridSize%nnzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                  ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                  ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                  ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                  ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))

          else if (willSendWest(iNeigh)) then

             ! insert send communications to west

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting four fields "//&
                     " at message "//trim(adjustl(WideGhostZoneSend%name)))
             end if

             fieldSectionSize= &
                  (xeSendWest(iNeigh)-xbSendWest(iNeigh)+1) * &
                  (yeSendWest(iNeigh)-ybSendWest(iNeigh)+1) * &
                  (GridSize%nnzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                  ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                  ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                  ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                  ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
          end if

       end do
    end if

    ! recv message set will contain recvs for all four directions

    willRecv = willRecvNorth .or. willRecvSouth .or. willRecvEast .or. willRecvWest
    if (dumpLocal) then
       do iNeigh = 1, nNeigh
          write(str(1),"(l)") willRecvNorth(iNeigh)
          write(str(2),"(l)") willRecvSouth(iNeigh)
          write(str(3),"(l)") willRecvEast(iNeigh)
          write(str(4),"(l)") willRecvWest(iNeigh)
          write(str(5),"(i8)") Brams2MpiProcNbr(Neigh%neigh(iNeigh))
          call MsgDump(h//" recv from MPI proc #"//trim(adjustl(str(5)))//&
               " willRecvNorth="//trim(adjustl(str(1)))//&
               " willRecvSouth="//trim(adjustl(str(2)))//&
               " willRecvEast="//trim(adjustl(str(3)))//&
               " willRecvWest="//trim(adjustl(str(4))))
       end do
    end if

    ! create message set for all recvs

    WideGhostZoneRecv => CreateMessageSet(&
         NameRecv, &
         recvDirection, &
         TagWideGhostZone, &
         willRecv, &
         Neigh)

    if ( associated(WideGhostZoneRecv)) then
       nMsgs = WideGhostZoneRecv%nMsgs
       cntMsg = 0
       do iNeigh = 1, nNeigh
          if (willRecvNorth(iNeigh)) then

             ! insert recv communications from north

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting four fields "//&
                     " at message "//trim(adjustl(WideGhostZoneRecv%name)))
             end if

             fieldSectionSize= &
                  (xeRecvNorth(iNeigh)-xbRecvNorth(iNeigh)+1) * &
                  (yeRecvNorth(iNeigh)-ybRecvNorth(iNeigh)+1) * &
                  (GridSize%nnzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                  ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                  ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                  ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                  ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))

          else if (willRecvSouth(iNeigh)) then

             ! insert recv communications from south

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting four fields "//&
                     " at message "//trim(adjustl(WideGhostZoneRecv%name)))
             end if

             fieldSectionSize= &
                  (xeRecvSouth(iNeigh)-xbRecvSouth(iNeigh)+1) * &
                  (yeRecvSouth(iNeigh)-ybRecvSouth(iNeigh)+1) * &
                  (GridSize%nnzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                  ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                  ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                  ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                  ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))

          else if (willRecvEast(iNeigh)) then

             ! insert recv communications from east

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting four fields "//&
                     " at message "//trim(adjustl(WideGhostZoneRecv%name)))
             end if

             fieldSectionSize= &
                  (xeRecvEast(iNeigh)-xbRecvEast(iNeigh)+1) * &
                  (yeRecvEast(iNeigh)-ybRecvEast(iNeigh)+1) * &
                  (GridSize%nnzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                  ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                  ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                  ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                  ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))

          else if (willRecvWest(iNeigh)) then

             ! insert recv communications from west

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting four fields "//&
                     " at message "//trim(adjustl(WideGhostZoneRecv%name)))
             end if

             fieldSectionSize= &
                  (xeRecvWest(iNeigh)-xbRecvWest(iNeigh)+1) * &
                  (yeRecvWest(iNeigh)-ybRecvWest(iNeigh)+1) * &
                  (GridSize%nnzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                  ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                  ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                  ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                  ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
          end if
       end do
    end if

    if (dumpLocal) then
       call MsgDump(h//" finishes with WideGhostZoneSend MessageSet:")
       call DumpMessageSet(WideGhostZoneSend)
       call MsgDump(h//" finishes with WideGhostZoneRecv MessageSet:")
       call DumpMessageSet(WideGhostZoneRecv)
    end if

  end subroutine CreateWideGhostZoneMessageSet



  subroutine DestroyWideGhostZoneMessageSet(&
       WideGhostZoneSend, WideGhostZoneRecv)

    type(MessageSet), pointer, intent(inout) :: WideGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: WideGhostZoneRecv

    character(len=*), parameter :: h="**(DestroyWideGhostZoneMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" will destroy WideGhostZoneSend/Recv")
    end if
    call DestroyMessageSet(WideGhostZoneSend)
    call DestroyMessageSet(WideGhostZoneRecv)
  end subroutine DestroyWideGhostZoneMessageSet




  subroutine UpdateFieldAdressAtMessageSet_2D(oneMessageSet, field, name)
    type(MessageSet), pointer, intent(in) :: oneMessageSet
    real, pointer, intent(in) :: field(:,:)
    character(len=*), intent(in) :: name
    character(len=*), parameter :: h="**(UpdateFieldAdressAtMessageSet_2D)**"

    include "updateFieldAdressBody.f90"

  end subroutine UpdateFieldAdressAtMessageSet_2D





  subroutine UpdateFieldAdressAtMessageSet_3D(oneMessageSet, field, name)
    type(MessageSet), pointer, intent(in) :: oneMessageSet
    real, pointer, intent(in) :: field(:,:,:)
    character(len=*), intent(in) :: name
    character(len=*), parameter :: h="**(UpdateFieldAdressAtMessageSet_3D)**"

    include "updateFieldAdressBody.f90"

  end subroutine UpdateFieldAdressAtMessageSet_3D





  subroutine UpdateFieldAdressAtMessageSet_4D(oneMessageSet, field, name)
    type(MessageSet), pointer, intent(in) :: oneMessageSet
    real, pointer, intent(in) :: field(:,:,:,:)
    character(len=*), intent(in) :: name
    character(len=*), parameter :: h="**(UpdateFieldAdressAtMessageSet_4D)**"

    include "updateFieldAdressBody.f90"

  end subroutine UpdateFieldAdressAtMessageSet_4D




  subroutine PostSendRecvMsgsFixedAdress(SendMsg, RecvMsg)

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
    character(len=*), parameter :: h="**(PostSendRecvMsgsFixedAdress)**"
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

          call AllocateMessageDataBuffer(RecvMsg%msgData(iRecv))

          ! post receive

          write(c0,"(i8)") iRecv
          call PostRecvMessageData(RecvMsg%msgData(iRecv), &
               RecvMsg%otherProc(iRecv), RecvMsg%tag, &
               RecvMsg%request(iRecv), &
               h//" for iRecv="//trim(adjustl(c0)))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty receive message set")
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

          ! allocate and fill send buffer with field sections to send

          call AllocateMessageDataBuffer(SendMsg%msgData(iSend))
          call FillMessageDataBuffer(SendMsg%msgData(iSend))

          ! post send message

          write(c0,"(i8)") iSend
          call PostSendMessageData(SendMsg%msgData(iSend), &
               SendMsg%otherProc(iSend), SendMsg%tag, &
               SendMsg%request(iSend), &
               h//" for iSend="//trim(adjustl(c0)))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty send message set")
       end if
    end if
  end subroutine PostSendRecvMsgsFixedAdress




  
  subroutine PostSendRecvMsgsVariableAdressArr(SendMsg, RecvMsg, scp, ufx, vfx, wfx)

    ! posts all nonblocking send and recv operations of
    ! a message set pair of variables

    type(MessageSet), pointer, intent(in) :: SendMsg
    type(MessageSet), pointer, intent(in) :: RecvMsg
    real, pointer, intent(in):: scp(:,:,:)
    real, intent(in):: ufx(:,:,:)
    real, intent(in):: vfx(:,:,:)
    real, intent(in):: wfx(:,:,:)

    integer :: iSend
    integer :: iRecv
    integer :: firstBuffer
    integer :: lastBuffer
    integer :: ierr
    type(MessageData), pointer :: msgData => null()
    type(FieldSection), pointer :: node => null()
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(PostSendRecvMsgsVariableAdressArr)**"
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

          call AllocateMessageDataBuffer(RecvMsg%msgData(iRecv))

          ! post receive

          write(c0,"(i8)") iRecv
          call PostRecvMessageData(RecvMsg%msgData(iRecv), &
               RecvMsg%otherProc(iRecv), RecvMsg%tag, &
               RecvMsg%request(iRecv), &
               h//" for recv #"//trim(adjustl(c0)))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty receive message set")
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

          ! allocate and fill send buffer with field sections to send

          call AllocateMessageDataBuffer(SendMsg%msgData(iSend))
          call FillMessageDataBufferVariableAdressArr(SendMsg%msgData(iSend), &
               scp, ufx, vfx, wfx)

          ! post send message

          write(c0,"(i8)") iSend
          call PostSendMessageData(SendMsg%msgData(iSend), &
               SendMsg%otherProc(iSend), SendMsg%tag, &
               SendMsg%request(iSend), &
               h//" for send #"//trim(adjustl(c0)))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty send message set")
       end if
    end if
  end subroutine PostSendRecvMsgsVariableAdressArr



  subroutine PostSendRecvMsgsVariableAdressScalar(SendMsg, RecvMsg, scp, ufx, vfx, wfx)

    ! posts all nonblocking send and recv operations of
    ! a message set pair of variables

    type(MessageSet), pointer, intent(in) :: SendMsg
    type(MessageSet), pointer, intent(in) :: RecvMsg
    real, intent(in):: scp
    real, intent(in):: ufx(:,:,:)
    real, intent(in):: vfx(:,:,:)
    real, intent(in):: wfx(:,:,:)

    integer :: iSend
    integer :: iRecv
    integer :: firstBuffer
    integer :: lastBuffer
    integer :: ierr
    type(MessageData), pointer :: msgData => null()
    type(FieldSection), pointer :: node => null()
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(PostSendRecvMsgsVariableAdressScalar)**"
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

          call AllocateMessageDataBuffer(RecvMsg%msgData(iRecv))

          ! post receive

          write(c0,"(i8)") iRecv
          call PostRecvMessageData(RecvMsg%msgData(iRecv), &
               RecvMsg%otherProc(iRecv), RecvMsg%tag, &
               RecvMsg%request(iRecv), &
               h//" for recv #"//trim(adjustl(c0)))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty receive message set")
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

          ! allocate and fill send buffer with field sections to send

          call AllocateMessageDataBuffer(SendMsg%msgData(iSend))
          call FillMessageDataBufferVariableAdressScalar(SendMsg%msgData(iSend), &
               scp, ufx, vfx, wfx)

          ! post send message

          write(c0,"(i8)") iSend
          call PostSendMessageData(SendMsg%msgData(iSend), &
               SendMsg%otherProc(iSend), SendMsg%tag, &
               SendMsg%request(iSend), &
               h//" for send #"//trim(adjustl(c0)))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty send message set")
       end if
    end if
  end subroutine PostSendRecvMsgsVariableAdressScalar  




  subroutine WaitSendRecvMsgsFixedAdress(SendMsg, RecvMsg)
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
    character(len=*), parameter :: h="**(WaitSendRecvMsgsFixedAdress)**"
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
          msgData => RecvMsg%msgData(recvNbr)
          if (dumpLocal) then
             write(c0,"(i8)") recvNbr
             write(c1,"(i8)") RecvMsg%otherProc(recvNbr)
             call MsgDump(h//" received message #"//trim(adjustl(c0))//&
                  " from MPI proc "//trim(adjustl(c1)))
          end if

          ! extract field sections from incoming buffer
          ! and store at destination fields

          call DecomposeMessageDataBuffer(RecvMsg%msgData(recvNbr))
          call DeallocateMessageDataBuffer(RecvMsg%msgData(recvNbr))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty receive message set")
       end if
    end if

    ! for all posted send messages, wait on pending request,
    ! deallocate buffer and empty request

    if (associated(SendMsg)) then
       !CDIR$ NOVECTOR
       do iSend = 1,SendMsg%nMsgs
          call parf_wait_any_nostatus(SendMsg%nMsgs, &
               SendMsg%request, sendNbr)
          call DeallocateMessageDataBuffer(SendMsg%msgData(sendNbr))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty send message set")
       end if
    end if
  end subroutine WaitSendRecvMsgsFixedAdress




  subroutine WaitSendRecvMsgsVariableAdress(SendMsg, RecvMsg, &
       msgStartZ, msgEndZ, scr, ufx_local, wfx_local, vfx_local)

    type(MessageSet), pointer, intent(in) :: SendMsg
    type(MessageSet), pointer, intent(in) :: RecvMsg
    integer, intent(in) :: msgStartZ
    integer, intent(in) :: msgEndZ
    real, pointer, intent(in) :: scr(:,:,:)
    real, pointer, intent(in) :: ufx_local(:,:,:)
    real, pointer, intent(in) :: vfx_local(:,:,:)
    real, pointer, intent(in) :: wfx_local(:,:,:)

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
    character(len=*), parameter :: h="**(WaitSendRecvMsgsVariableAdress)**"
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
          msgData => RecvMsg%msgData(recvNbr)
          if (dumpLocal) then
             write(c0,"(i8)") recvNbr
             write(c1,"(i8)") RecvMsg%otherProc(recvNbr)
             call MsgDump(h//" received message #"//trim(adjustl(c0))//&
                  " from MPI proc "//trim(adjustl(c1)))
          end if

          ! extract field sections from incoming buffer
          ! and store at destination fields

          call DecomposeMessageDataBuffer(RecvMsg%msgData(recvNbr), &
               msgStartZ, msgEndZ, scr, ufx_local, vfx_local, wfx_local)
          call DeallocateMessageDataBuffer(RecvMsg%msgData(recvNbr))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty receive message set")
       end if
    end if

    ! for all posted send messages, wait on pending request,
    ! deallocate buffer and empty request

    if (associated(SendMsg)) then
       !CDIR$ NOVECTOR
       do iSend = 1,SendMsg%nMsgs
          call parf_wait_any_nostatus(SendMsg%nMsgs, &
               SendMsg%request, sendNbr)
          call DeallocateMessageDataBuffer(SendMsg%msgData(sendNbr))
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" empty send message set")
       end if
    end if
  end subroutine WaitSendRecvMsgsVariableAdress
end module ModMessageSet


