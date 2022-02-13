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
  ! Procedure PostRecvSendMsgs initiate all send/recv
  ! operations of one set. Procedure WaitSendRecvMsgs finalize
  ! all send/recv operations of the same set.
  ! Copy field sections to the contiguous buffer and the reverse
  ! opetation are included on both message passing procedures.
  !
  ! To overlap computation with communication, procedure
  ! PostRecvSendMsgs should be invoked as early as possible,
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
       IncludeDomainBoundaries

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
       FillMessageDataBufferWithFieldSectionData, &
       ExtractFieldSectionDataFromMessageDataBuffer, &
       AllocateMessageDataBuffer, &
       DeallocateMessageDataBuffer, &
       DestroyMessageData

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

  public :: PostRecvSendMsgs
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
    logical, parameter :: dumpLocal=.true.
 
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
    



  subroutine CreateAcousticMessageSet(&
       gridId, GridSize, ParEnv, Neigh, &
       GlobalOwn, &
       GlobalOwnWithBC, &
       GlobalWithGhost, &
       AcouSendU, AcouRecvU, &
       AcouSendV, AcouRecvV, &
       AcouSendP, AcouRecvP, &
       AcouSendUV, AcouRecvUV, &
       AcouSendWP, AcouRecvWP)

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
    type(MessageSet), pointer, intent(inout) :: AcouSendP
    type(MessageSet), pointer, intent(inout) :: AcouRecvP
    type(MessageSet), pointer, intent(inout) :: AcouSendUV
    type(MessageSet), pointer, intent(inout) :: AcouRecvUV
    type(MessageSet), pointer, intent(inout) :: AcouSendWP
    type(MessageSet), pointer, intent(inout) :: AcouRecvWP

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    character(len=*), parameter :: h="**(CreateAcousticMessageSet)**"
    logical, parameter :: dumpLocal=.true.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" starts with null GlobalOwn")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    AcouSendU => null()
    AcouRecvU => null()
    AcouSendV => null()
    AcouRecvV => null()
    AcouSendP => null()
    AcouRecvP => null()
    AcouSendUV => null()
    AcouRecvUV => null()
    AcouSendWP => null()
    AcouRecvWP => null()

    if (associated(Neigh)) then

       if (dumpLocal) then
          call MsgDump(h//" will create "//&
               " AcouSend/RecvU, AcouSend/RecvV, AcouSend/RecvP,"//&
               " AcouSend/RecvUV and AcouSend/RecvWP")
       end if

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateAcousticSendRecvU(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalOwn, GlobalWithGhost, &
            AcouSendU, AcouRecvU)

       call CreateAcousticSendRecvV(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalOwn, GlobalWithGhost, &
            AcouSendV, AcouRecvV)

       call CreateAcousticSendRecvP(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalOwn, GlobalWithGhost, &
            AcouSendP, AcouRecvP)

       call CreateAcousticSendRecvUV(&
            gridId, nMachs, nNeigh, myNum, &
            GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
            AcouSendUV, AcouRecvUV)

       call CreateAcousticSendRecvWP(&
            gridId, nMachs, nNeigh, myNum, &
            GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
            AcouSendWP, AcouRecvWP)
    end if
  end subroutine CreateAcousticMessageSet





  subroutine DestroyAcousticMessageSet( &
       AcouSendU, AcouRecvU, &
       AcouSendV, AcouRecvV, &
       AcouSendP, AcouRecvP, &
       AcouSendUV, AcouRecvUV, &
       AcouSendWP, AcouRecvWP)

    type(MessageSet), pointer, intent(inout) :: AcouSendU
    type(MessageSet), pointer, intent(inout) :: AcouRecvU
    type(MessageSet), pointer, intent(inout) :: AcouSendV
    type(MessageSet), pointer, intent(inout) :: AcouRecvV
    type(MessageSet), pointer, intent(inout) :: AcouSendP
    type(MessageSet), pointer, intent(inout) :: AcouRecvP
    type(MessageSet), pointer, intent(inout) :: AcouSendUV
    type(MessageSet), pointer, intent(inout) :: AcouRecvUV
    type(MessageSet), pointer, intent(inout) :: AcouSendWP
    type(MessageSet), pointer, intent(inout) :: AcouRecvWP

    character(len=*), parameter :: h="**(DestroyAcousticMessageSet)**"
    logical, parameter :: dumpLocal=.true.

    if (dumpLocal) then
       call MsgDump(h//" will destroy "//&
            " AcouSend/RecvU, AcouSend/RecvV, AcouSend/RecvP,"//&
            " AcouSend/RecvUV and AcouSend/RecvWP")
    end if

    call DestroyMessageSet(AcouSendU)
    call DestroyMessageSet(AcouRecvU)

    call DestroyMessageSet(AcouSendV)
    call DestroyMessageSet(AcouRecvV)

    call DestroyMessageSet(AcouSendP)
    call DestroyMessageSet(AcouRecvP)

    call DestroyMessageSet(AcouSendUV)
    call DestroyMessageSet(AcouRecvUV)

    call DestroyMessageSet(AcouSendWP)
    call DestroyMessageSet(AcouRecvWP)
  end subroutine DestroyAcousticMessageSet





  subroutine CreateDn0MessageSet(&
       gridId, GridSize, ParEnv, Neigh, &
       GlobalOwn, GlobalWithGhost, &
       SendDn0u, RecvDn0u, &
       SendDn0v, RecvDn0v)

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

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    character(len=*), parameter :: h="**(CreateDn0MessageSet)**"
    logical, parameter :: dumpLocal=.true.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" starts with null GlobalOwn")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    SendDn0u => null()
    RecvDn0u => null()
    SendDn0v => null()
    RecvDn0v => null()

    if (associated(Neigh)) then

       if (dumpLocal) then
          call MsgDump(h//" will create Send/RecvDn0u and Send/RecvDn0v")
       end if

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateSendRecvDn0u(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalOwn, GlobalWithGhost, &
            SendDn0u, RecvDn0u)

       call CreateSendRecvDn0v(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalOwn, GlobalWithGhost, &
            SendDn0v, RecvDn0v)
    end if
  end subroutine CreateDn0MessageSet





  subroutine DestroyDn0MessageSet( &
       SendDn0u, RecvDn0u, SendDn0v, RecvDn0v)

    type(MessageSet), pointer, intent(inout) :: SendDn0u
    type(MessageSet), pointer, intent(inout) :: RecvDn0u
    type(MessageSet), pointer, intent(inout) :: SendDn0v
    type(MessageSet), pointer, intent(inout) :: RecvDn0v
    character(len=*), parameter :: h="**(DestroyDn0MessageSet)**"
    logical, parameter :: dumpLocal=.true.

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
       SendG3D, RecvG3D)

    integer, intent(in) :: gridId
    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(NamelistFile), pointer, intent(in) :: Ramsin
    type(MessageSet), pointer, intent(inout) :: SendG3D
    type(MessageSet), pointer, intent(inout) :: RecvG3D

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    integer :: g3d_spread
    integer :: g3d_smoothh
    character(len=*), parameter :: h="**(CreateG3DMessageSet)**"
    logical, parameter :: dumpLocal=.true.

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
         (g3d_spread /= 0 .or. g3d_smoothh /= 0)) then

       if (dumpLocal) then
          call MsgDump(h//" will create Send/RecvG3D")
       end if

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateG3DSendRecv(&
            gridId, nMachs, nNeigh, myNum, &
            g3d_spread, g3d_smoothh, &
            GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
            SendG3D, RecvG3D)
    end if
  end subroutine CreateG3DMessageSet





  subroutine DestroyG3DMessageSet( &
       SendG3D, RecvG3D)

    type(MessageSet), pointer, intent(inout) :: SendG3D
    type(MessageSet), pointer, intent(inout) :: RecvG3D
    character(len=*), parameter :: h="**(DestroyG3DMessageSet)**"
    logical, parameter :: dumpLocal=.true.

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
       SelectedGhostZoneSend, SelectedGhostZoneRecv)

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


    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    character(len=*), parameter :: h="**(CreateSelectedGhostZoneMessageSet)**"
    logical, parameter :: dumpLocal=.true.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwnWithBC)) then
       call fatal_error(h//" starts with null GlobalOwnWithBC")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    SelectedGhostZoneSend => null()
    SelectedGhostZoneRecv => null()

    if (associated(Neigh)) then

       if (dumpLocal) then
          call MsgDump(h//" will create SelectedGhostZoneSend/Recv")
       end if

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateSelectedGhostZoneSendRecv(&
            gridId, nMachs, nNeigh, myNum, num_var, vtab_r, &
            GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
            SelectedGhostZoneSend, SelectedGhostZoneRecv)
    end if
  end subroutine CreateSelectedGhostZoneMessageSet





  subroutine DestroySelectedGhostZoneMessageSet( &
       SelectedGhostZoneSend, SelectedGhostZoneRecv)

    type(MessageSet), pointer, intent(inout) :: SelectedGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: SelectedGhostZoneRecv
    character(len=*), parameter :: h="**(DestroySelectedGhostZoneMessageSet)**"
    logical, parameter :: dumpLocal=.true.

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
       AllGhostZoneSend, AllGhostZoneRecv)

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


    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    character(len=*), parameter :: h="**(CreateAllGhostZoneMessageSet)**"
    logical, parameter :: dumpLocal=.true.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwnWithBC)) then
       call fatal_error(h//" starts with null GlobalOwnWithBC")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    AllGhostZoneSend => null()
    AllGhostZoneRecv => null()

    if (associated(Neigh)) then

       if (dumpLocal) then
          call MsgDump(h//" will create AllGhostZoneSend/Recv")
       end if

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateAllGhostZoneSendRecv(&
            gridId, nMachs, nNeigh, myNum, num_var, vtab_r, &
            GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
            AllGhostZoneSend, AllGhostZoneRecv)
    end if
  end subroutine CreateAllGhostZoneMessageSet





  subroutine DestroyAllGhostZoneMessageSet( &
       AllGhostZoneSend, AllGhostZoneRecv)

    type(MessageSet), pointer, intent(inout) :: AllGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: AllGhostZoneRecv
    character(len=*), parameter :: h="**(DestroyAllGhostZoneMessageSet)**"
    logical, parameter :: dumpLocal=.true.

    if (dumpLocal) then
       call MsgDump(h//" will destroy AllGhostZoneSend/Recv")
    end if

    call DestroyMessageSet(AllGhostZoneSend)
    call DestroyMessageSet(AllGhostZoneRecv)

  end subroutine DestroyAllGhostZoneMessageSet





  subroutine BuildUnion(nNeigh, &
       xbComm1, xeComm1, ybComm1, yeComm1, willComm1, &
       xbComm2, xeComm2, ybComm2, yeComm2, willComm2)
    integer, intent(in) :: nNeigh
    integer, intent(in) :: xbComm1(nNeigh)
    integer, intent(in) :: xeComm1(nNeigh)
    integer, intent(in) :: ybComm1(nNeigh)
    integer, intent(in) :: yeComm1(nNeigh)
    logical, intent(in) :: willComm1(nNeigh)

    integer, intent(inout) :: xbComm2(nNeigh)
    integer, intent(inout) :: xeComm2(nNeigh)
    integer, intent(inout) :: ybComm2(nNeigh)
    integer, intent(inout) :: yeComm2(nNeigh)
    logical, intent(inout) :: willComm2(nNeigh)

    integer :: iNeigh
    character(len=8) :: c0, c1, c2, c3, c4
    character(len=128) :: inter1, inter2
    character(len=*), parameter :: h="**(BuildUnion)**"
    logical, parameter :: dumpLocal=.true.

    do iNeigh = 1, nNeigh
       if (willComm1(iNeigh) .and. willComm2(iNeigh)) then
          write(c0,"(i8)") iNeigh
          inter1="inter1"
          write(c1,"(i8)") xbComm1(iNeigh)
          write(c2,"(i8)") xeComm1(iNeigh)
          inter1=trim(inter1)//"("//trim(adjustl(c1))//":"//trim(adjustl(c2))
          write(c1,"(i8)") ybComm1(iNeigh)
          write(c2,"(i8)") yeComm1(iNeigh)
          inter1=trim(inter1)//","//trim(adjustl(c1))//":"//trim(adjustl(c2))//")"
          inter2="inter2"
          write(c1,"(i8)") xbComm2(iNeigh)
          write(c2,"(i8)") xeComm2(iNeigh)
          inter2=trim(inter2)//"("//trim(adjustl(c1))//":"//trim(adjustl(c2))
          write(c1,"(i8)") ybComm2(iNeigh)
          write(c2,"(i8)") yeComm2(iNeigh)
          inter2=trim(inter2)//","//trim(adjustl(c1))//":"//trim(adjustl(c2))//")"
          call fatal_error(h//" both willComm("//trim(adjustl(c0))//&
               "); intervals="//trim(adjustl(inter1))//", "//trim(adjustl(inter2)))
       else if (willComm1(iNeigh)) then
          xbComm2(iNeigh) = xbComm1(iNeigh)
          xeComm2(iNeigh) = xeComm1(iNeigh)
          ybComm2(iNeigh) = ybComm1(iNeigh)
          yeComm2(iNeigh) = yeComm1(iNeigh)
          willComm2(iNeigh) = willComm1(iNeigh)
          if (dumpLocal) then
             write(c0,"(i8)") iNeigh
             write(c1,"(i8)") xbComm1(iNeigh)
             write(c2,"(i8)") xeComm1(iNeigh)
             write(c3,"(i8)") ybComm1(iNeigh)
             write(c4,"(i8)") yeComm1(iNeigh)
             call MsgDump(h//" updated Comm2("//trim(adjustl(c0))//") to ["//&
                  trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
                  trim(adjustl(c3))//":"//trim(adjustl(c4))//"]")
          end if
       end if
    end do
  end subroutine BuildUnion





  subroutine CreateAcousticSendRecvUV(&
       gridId, nMachs, nNeigh, myNum, &
       GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
       AcouSendUV, AcouRecvUV)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(GridDims), pointer, intent(in) :: GridSize
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: AcouSendUV
    type(MessageSet), pointer, intent(inout) :: AcouRecvUV

    integer, parameter :: TagUV=25
    character(len=*), parameter :: NameSendUV="AcouSendUV"
    character(len=*), parameter :: NameRecvUV="AcouRecvUV"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvUV)**"
    character(len=30) :: tmp_name
    character(len=8) :: str(10)
    logical, parameter :: dumpLocal=.true.

    if (dumpLocal) then
       call MsgDump(h// "starts")
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

    vTabPtr => null()
    if(dyncore_flag==2) then
       tmp_name='UC'
    else
       tmp_name='UP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! include UP on field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendUV)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvUV)

    ! get field VP

    vTabPtr => null()
    if(dyncore_flag==2) then
       tmp_name='VC'
    else
       tmp_name='VP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! include VP on field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendUV)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvUV)
    if (dumpLocal) then
       call MsgDump(h//" finishes with AcouSendUV MessageSet:")
       call DumpMessageSet(AcouSendUV)
       call MsgDump(h//" finishes with AcouRecvUV MessageSet:")
       call DumpMessageSet(AcouRecvUV)
    end if
  end subroutine CreateAcousticSendRecvUV





  subroutine CreateAcousticSendRecvWP(&
       gridId, nMachs, nNeigh, myNum, &
       GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
       AcouSendWP, AcouRecvWP)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(GridDims), pointer, intent(in) :: GridSize
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: AcouSendWP
    type(MessageSet), pointer, intent(inout) :: AcouRecvWP

    integer, parameter :: TagWP=26
    character(len=*), parameter :: NameSendWP="AcouSendWP"
    character(len=*), parameter :: NameRecvWP="AcouRecvWP"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvWP)**"
    character(len=30) :: tmp_name
    character(len=8) :: str(10)
    logical, parameter :: dumpLocal=.true.

 
    if (dumpLocal) then
       call MsgDump(h// "starts")
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

    vTabPtr => null()
    if(dyncore_flag==2) then
       tmp_name='WC'
    else
       tmp_name='WP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendWP)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvWP)

    ! get field VP

    vTabPtr => null()
    if(dyncore_flag==2) then
       tmp_name='PC'
    else
       tmp_name='PP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendWP)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvWP)
    if (dumpLocal) then
       call MsgDump(h//" finishes with AcouSendWP MessageSet:")
       call DumpMessageSet(AcouSendWP)
       call MsgDump(h//" finishes with AcouRecvWP MessageSet:")
       call DumpMessageSet(AcouRecvWP)
    end if
  end subroutine CreateAcousticSendRecvWP





  subroutine CreateSelectedGhostZoneSendRecv(&
       gridId, nMachs, nNeigh, myNum, num_var, vtab_r, &
       GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
       SelectedGhostZoneSend, SelectedGhostZoneRecv)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum
    integer, intent(in) :: num_var(:)
    type(var_tables_r), target, intent(in) ::  vtab_r(:,:)

    type(GridDims), pointer, intent(in) :: GridSize
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: SelectedGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: SelectedGhostZoneRecv

    integer :: vTabNbr
    integer, parameter :: TagSelectedGhostZone=27
    character(len=*), parameter :: NameSendSelectedGhostZone="SelectedGhostZoneSend"
    character(len=*), parameter :: NameRecvSelectedGhostZone="SelectedGhostZoneRecv"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateSelectedGhostZoneSendRecv)**"
    character(len=8) :: str(10)
    logical, parameter :: dumpLocal=.true.

 
    if (dumpLocal) then
       call MsgDump(h// "starts")
    end if

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
          vTabPtr => vtab_r(vTabNbr,gridId)

          ! build field sections to be sent and received

          call InsertFieldSectionAtMessageSet(&
               myNum, vTabPtr, Neigh, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, &
               SelectedGhostZoneSend)
          call InsertFieldSectionAtMessageSet(&
               myNum, vTabPtr, Neigh, GlobalWithGhost, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
               SelectedGhostZoneRecv)
       end if
    end do
    if (dumpLocal) then
       call MsgDump(h//" finishes with SelectedGhostZoneSend MessageSet:")
       call DumpMessageSet(SelectedGhostZoneSend)
       call MsgDump(h//" finishes with SelectedGhostZoneRecv MessageSet:")
       call DumpMessageSet(SelectedGhostZoneRecv)
    end if
  end subroutine CreateSelectedGhostZoneSendRecv





  subroutine CreateAllGhostZoneSendRecv(&
       gridId, nMachs, nNeigh, myNum, num_var, vtab_r, &
       GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
       AllGhostZoneSend, AllGhostZoneRecv)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum
    integer, intent(in) :: num_var(:)
    type(var_tables_r), target, intent(in) ::  vtab_r(:,:)

    type(GridDims), pointer, intent(in) :: GridSize
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: AllGhostZoneSend
    type(MessageSet), pointer, intent(inout) :: AllGhostZoneRecv

    integer :: vTabNbr
    integer, parameter :: TagAllGhostZone=28
    character(len=*), parameter :: NameSendAllGhostZone="AllGhostZoneSend"
    character(len=*), parameter :: NameRecvAllGhostZone="AllGhostZoneRecv"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAllGhostZoneSendRecv)**"
    character(len=8) :: str(10)
    logical, parameter :: dumpLocal=.true.

 
    if (dumpLocal) then
       call MsgDump(h// "starts")
    end if

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

       vTabPtr => vtab_r(vTabNbr,gridId)

       if (&
            trim(adjustl(vTabPtr%name)) /= "LPU" .and. &
            trim(adjustl(vTabPtr%name)) /= "LPV" .and. &
            trim(adjustl(vTabPtr%name)) /= "LPW" ) then

          ! build field sections to be sent and received

          call InsertFieldSectionAtMessageSet(&
               myNum, vTabPtr, Neigh, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, AllGhostZoneSend)
          call InsertFieldSectionAtMessageSet(&
               myNum, vTabPtr, Neigh, GlobalWithGhost, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AllGhostZoneRecv)
       end if
    end do
    if (dumpLocal) then
       call MsgDump(h//" finishes with AllGhostZoneSend MessageSet:")
       call DumpMessageSet(AllGhostZoneSend)
       call MsgDump(h//" finishes with AllGhostZoneRecv MessageSet:")
       call DumpMessageSet(AllGhostZoneRecv)
    end if
  end subroutine CreateAllGhostZoneSendRecv





  subroutine CreateSendRecvDn0u(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalOwn, GlobalWithGhost, &
       SendDn0u, RecvDn0u)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: SendDn0u
    type(MessageSet), pointer, intent(inout) :: RecvDn0u

    integer, parameter :: TagDn0u=29
    character(len=*), parameter :: NameSendDn0u="SendDn0u"
    character(len=*), parameter :: NameRecvDn0u="RecvDn0u"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateSendRecvDn0u)**"
    character(len=30) :: tmp_name
    character(len=8) :: str(10)
    logical, parameter :: dumpLocal=.true.

 
    if (dumpLocal) then
       call MsgDump(h// "starts")
    end if

    ! SendDn0u, RecvDn0u:
    ! messages update GlobalOwn [xe+1:xe+1,yb:ye]

    call NodesRegionsSendRecv(&
         nMachs=nMachs, &
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

    vTabPtr => null()
    tmp_name='DN0U'
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendDn0u)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvDn0u)
    if (dumpLocal) then
       call MsgDump(h//" finishes with SendDn0u MessageSet:")
       call DumpMessageSet(SendDn0u)
       call MsgDump(h//" finishes with RecvDn0u MessageSet:")
       call DumpMessageSet(RecvDn0u)
    end if
  end subroutine CreateSendRecvDn0u





  subroutine CreateSendRecvDn0v(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalOwn, GlobalWithGhost, &
       SendDn0v, RecvDn0v)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: SendDn0v
    type(MessageSet), pointer, intent(inout) :: RecvDn0v

    integer, parameter :: TagDn0v=30
    character(len=*), parameter :: NameSendDn0v="SendDn0v"
    character(len=*), parameter :: NameRecvDn0v="RecvDn0v"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateSendRecvDn0v)**"
    character(len=30) :: tmp_name
    character(len=8) :: str(10)
    logical, parameter :: dumpLocal=.true.

 
    if (dumpLocal) then
       call MsgDump(h// "starts")
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

    vTabPtr => null()
    tmp_name='DN0V'
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendDn0v)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvDn0v)
    if (dumpLocal) then
       call MsgDump(h//" finishes with SendDn0v MessageSet:")
       call DumpMessageSet(SendDn0v)
       call MsgDump(h//" finishes with RecvDn0v MessageSet:")
       call DumpMessageSet(RecvDn0v)
    end if
  end subroutine CreateSendRecvDn0v




  subroutine CreateG3DSendRecv(&
       gridId, nMachs, nNeigh, myNum, &
       g3d_spread, g3d_smoothh, &
       GridSize, Neigh, GlobalOwnWithBC, GlobalWithGhost, &
       SendG3D, RecvG3D)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum
    integer, intent(in) :: g3d_spread
    integer, intent(in) :: g3d_smoothh

    type(GridDims), pointer, intent(in) :: GridSize
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: SendG3D
    type(MessageSet), pointer, intent(inout) :: RecvG3D

    integer :: vTabNbr
    integer, parameter :: TagG3D=31
    character(len=*), parameter :: NameSendG3D="SendG3D"
    character(len=*), parameter :: NameRecvG3D="RecvG3D"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateG3DSendRecv)**"
    character(len=30) :: tmp_name
    character(len=8) :: str(10)
    logical, parameter :: dumpLocal=.true.

 
    if (dumpLocal) then
       call MsgDump(h// "starts")
    end if

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

       vTabPtr => null()
       tmp_name='TTENS'
       call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbSend, xeSend, ybSend, yeSend, willSend, &
            SendG3D)
       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
            RecvG3D)

       vTabPtr => null()
       tmp_name='QVTTENS'
       call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbSend, xeSend, ybSend, yeSend, willSend, &
            SendG3D)
       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
            RecvG3D)
    end if

    ! when g3d_smoothh is selected, send and receive fields THSRC and RTSRC

    if (g3d_smoothh == 1) then

       vTabPtr => null()
       tmp_name='THSRC'
       call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbSend, xeSend, ybSend, yeSend, willSend, &
            SendG3D)
       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
            RecvG3D)

       vTabPtr => null()
       tmp_name='RTSRC'
       call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbSend, xeSend, ybSend, yeSend, willSend, &
            SendG3D)
       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
            RecvG3D)
    end if
    if (dumpLocal) then
       call MsgDump(h//" finishes with SendG3D MessageSet:")
       call DumpMessageSet(SendG3D)
       call MsgDump(h//" finishes with RecvG3D MessageSet:")
       call DumpMessageSet(RecvG3D)
    end if
  end subroutine CreateG3DSendRecv





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
    real, pointer :: PNull(:,:) => null()
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
  end subroutine InsertFieldSectionAtMessageSet





  subroutine CreateAcousticSendRecvU(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalOwn, GlobalWithGhost, &
       AcouSendU, AcouRecvU)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: AcouSendU
    type(MessageSet), pointer, intent(inout) :: AcouRecvU

    integer, parameter :: TagU=22
    character(len=*), parameter :: NameSendU="AcouSendU"
    character(len=*), parameter :: NameRecvU="AcouRecvU"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvU)**"
    logical, parameter :: dumpLocal=.true.
    character(len=8) :: str(10)
    character(len=30) :: tmp_name

 
    if (dumpLocal) then
       call MsgDump(h// "starts")
    end if

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

    vTabPtr => null()
    if(dyncore_flag==2) then
       tmp_name='UC'
    else
       tmp_name='UP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendU)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvU)
    if (dumpLocal) then
       call MsgDump(h//" finishes with AcouSendU MessageSet:")
       call DumpMessageSet(AcouSendU)
       call MsgDump(h//" finishes with AcouRecvU MessageSet:")
       call DumpMessageSet(AcouRecvU)
    end if
  end subroutine CreateAcousticSendRecvU





  subroutine CreateAcousticSendRecvV(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalOwn, GlobalWithGhost, &
       AcouSendV, AcouRecvV)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: AcouSendV
    type(MessageSet), pointer, intent(inout) :: AcouRecvV

    integer, parameter :: TagV=23
    character(len=*), parameter :: NameSendV="AcouSendV"
    character(len=*), parameter :: NameRecvV="AcouRecvV"
    character(len=30) :: tmp_name
    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)


    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvV)**"
    logical, parameter :: dumpLocal=.true.
    character(len=8) :: str(10)

 
    if (dumpLocal) then
       call MsgDump(h// "starts")
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

    vTabPtr => null()
    if(dyncore_flag==2) then
       tmp_name='VC'
    else
       tmp_name='VP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendV)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvV)
    if (dumpLocal) then
       call MsgDump(h//" finishes with AcouSendV MessageSet:")
       call DumpMessageSet(AcouSendV)
       call MsgDump(h//" finishes with AcouRecvV MessageSet:")
       call DumpMessageSet(AcouRecvV)
    end if
  end subroutine CreateAcousticSendRecvV





  subroutine CreateAcousticSendRecvP(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalOwn, GlobalWithGhost, &
       AcouSendP, AcouRecvP)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(MessageSet), pointer, intent(inout) :: AcouSendP
    type(MessageSet), pointer, intent(inout) :: AcouRecvP

    integer, parameter :: TagP=24
    character(len=*), parameter :: NameSendP="AcouSendP"
    character(len=*), parameter :: NameRecvP="AcouRecvP"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! second region to build union of regions

    integer :: xbSend_2(nNeigh)
    integer :: xeSend_2(nNeigh)
    integer :: ybSend_2(nNeigh)
    integer :: yeSend_2(nNeigh)
    integer :: xbRecv_2(nNeigh)
    integer :: xeRecv_2(nNeigh)
    integer :: ybRecv_2(nNeigh)
    integer :: yeRecv_2(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    logical :: willSend_2(nNeigh)
    logical :: willRecv_2(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvP)**"
    character(len=30) :: tmp_name
    character(len=8) :: str(10)
    logical, parameter :: dumpLocal=.true.

 
    if (dumpLocal) then
       call MsgDump(h// "starts")
    end if

    ! AcouSendP, AcouRecvP: union of
    !               GlobalOwn [xe+1:xe+1,yb:ye] with
    !               GlobalOwn [xb:xe,ye+1:ye+1]

    ! first part of the union ([xe+1:xe+1,yb:ye])

    xbToBeUpdated = GlobalOwn%xe+1
    xeToBeUpdated = xbToBeUpdated
    ybToBeUpdated = GlobalOwn%yb
    yeToBeUpdated = GlobalOwn%ye

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalOwn, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
         "AcouSend/RecvP")

    ! second part of the union ([xb:xe,ye+1:ye+1])

    xbToBeUpdated = GlobalOwn%xb
    xeToBeUpdated = GlobalOwn%xe
    ybToBeUpdated = GlobalOwn%ye + 1
    yeToBeUpdated = ybToBeUpdated

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalOwn, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend_2, xeSend_2, ybSend_2, yeSend_2, willSend_2, &
         xbRecv_2, xeRecv_2, ybRecv_2, yeRecv_2, willRecv_2, &
         "AcouSend/RecvP")

    ! make union

    call BuildUnion(nNeigh, &
         xbSend_2, xeSend_2, ybSend_2, yeSend_2, willSend_2, &
         xbSend, xeSend, ybSend, yeSend, willSend)
    call BuildUnion(nNeigh, &
         xbRecv_2, xeRecv_2, ybRecv_2, yeRecv_2, willRecv_2, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    AcouSendP => CreateMessageSet(&
         NameSendP, &
         sendDirection, &
         TagP, &
         willSend, &
         Neigh)
    AcouRecvP => CreateMessageSet(&
         NameRecvP, &
         recvDirection, &
         TagP, &
         willRecv, &
         Neigh)

    ! get field

    vTabPtr => null()
    if(dyncore_flag==2) then
       tmp_name='PC'
    else
       tmp_name='PP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendP)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvP)
    if (dumpLocal) then
       call MsgDump(h//" finishes with AcouSendP MessageSet:")
       call DumpMessageSet(AcouSendP)
       call MsgDump(h//" finishes with AcouRecvP MessageSet:")
       call DumpMessageSet(AcouRecvP)
    end if
  end subroutine CreateAcousticSendRecvP





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
    logical, parameter :: dumpLocal=.true.

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

          ! allocate and fill send buffer with field sections to send

          call AllocateMessageDataBuffer(SendMsg%msgData(iSend))
          call FillMessageDataBufferWithFieldSectionData(SendMsg%msgData(iSend))

          ! post send message

          write(c0,"(i8)") iSend
          call PostSendMessageData(SendMsg%msgData(iSend), &
               SendMsg%otherProc(iSend), SendMsg%tag, &
               SendMsg%request(iSend), &
               h//" for iSend="//trim(adjustl(c0)))
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
    logical, parameter :: dumpLocal=.true.

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

          call ExtractFieldSectionDataFromMessageDataBuffer(RecvMsg%msgData(recvNbr))
          call DeallocateMessageDataBuffer(RecvMsg%msgData(recvNbr))
       end do
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
    end if
  end subroutine WaitSendRecvMsgs
end module ModMessageSet
