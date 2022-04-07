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

  use ModNodeDimensions, only: &
       NodeDimensions
  
  use ModNeighbourNodes, only: &
       NeighbourNodes, &
       NodesToSendRecvMessages, &
       IncludeDomainBoundaries, &
       DumpNeighbourNodes

  use ModDomainDecomp, only: &
       DomainDecomp

  use ModFieldSection, only: &
       FieldSection, &
       CreateFieldSection, &
       UpdateFieldAdress

  use ModFieldSectionList, only: &
       FieldSectionNode
  
  use ModMessageData, only: &
       MessageData, &
       CreateMessageData, &
       DumpMessageData, &
       AppendFieldSectionToMessageData, &
       PostRecvMessageData, &
       PostSendMessageData, &
       FillMessageDataBuffer, &
       DecomposeMessageDataBuffer, &
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

  public :: CreateAcoustNewMessageSet
  public :: DestroyAcoustNewMessageSet

  public :: CreateWideGhostZoneMessageSet
  public :: DestroyWideGhostZoneMessageSet

  public :: CreateAdvMntMessageSet
  public :: DestroyAdvMntMessageSet
  
  public :: PostSendRecvMsgs
  public :: WaitSendRecvMsgs

  public :: UpdateFieldAdressAtAdvMnt
  public :: UpdateFieldAdressAtAcoustNew  
  public :: UpdateSendFieldAdressAtAdvectcRk
  public :: UpdateRecvFieldAdressAtAdvectcRk
  
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
     module procedure InsertFieldSectionAtSendRecvMessageSet_1D
     module procedure InsertFieldSectionAtSendRecvMessageSet_2D
     module procedure InsertFieldSectionAtSendRecvMessageSet_3D
     module procedure InsertFieldSectionAtSendRecvMessageSet_4D
  end interface InsertFieldSectionAtSendRecvMessageSet

  interface InsertFieldSectionAtMessageSet
     module procedure InsertFieldSectionAtMessageSetFromVTab
     module procedure InsertFieldSectionAtMessageSet_1D  
     module procedure InsertFieldSectionAtMessageSet_2D  
     module procedure InsertFieldSectionAtMessageSet_3D  
     module procedure InsertFieldSectionAtMessageSet_4D
  end interface InsertFieldSectionAtMessageSet

  interface PostSendRecvMsgs
     module procedure PostSendRecvMsgsFixedAdress
     module procedure PostSendRecvMsgsFixedAdress1D
  end interface PostSendRecvMsgs

  interface WaitSendRecvMsgs
     module procedure WaitSendRecvMsgsFixedAdress
     module procedure WaitSendRecvMsgsFixedAdress1D
  end interface WaitSendRecvMsgs

  interface CreateAcoustNewOneMessageSet
     module procedure CreateAcoustNewOneMessageSet3D
     module procedure CreateAcoustNewOneMessageSet1D
  end interface CreateAcoustNewOneMessageSet

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
       call MsgDump(h// "starts computing field sections for "//&
            trim(adjustl(NameSend))//" and "//trim(adjustl(NameRecv)))
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
       call MsgDump(h//" finishes creating Message Sets with empty FieldSections:")
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
             call fatal_error(h//" unknown idim_type="//trim(adjustl(c0))//&
                  " for field "//trim(adjustl(vTabPtr%name)))
          end select
          call AppendFieldSectionToMessageData(oneFieldSection, Msgs%msgData(cntMsg))
       end if
    end do
  end subroutine InsertFieldSectionAtMessageSetFromVTab



  subroutine InsertFieldSectionAtMessageSet_1D(&
       myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
       zStart, zEnd, xbComm, xeComm, ybComm, yeComm, willComm, &
       Msgs, kMax)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! mynum is this BRAMS process number;
    ! It sends data on a send MessageSet variable or
    ! it receives data on a receive MessageSet variable

    integer, intent(in) :: myNum

    ! field memory address

    real, pointer, intent(in) :: field(:)

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

    integer, intent(in) :: zStart
    integer, intent(in) :: zEnd
    integer, intent(in) :: xbComm(:)
    integer, intent(in) :: xeComm(:)
    integer, intent(in) :: ybComm(:)
    integer, intent(in) :: yeComm(:)
    integer, intent(in) :: kMax

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
    character(len=*), parameter :: h="**(InsertFieldSectionAtMessageSet_1D)**"

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
          case (1)
             oneFieldSection => CreateFieldSection(&
                  field, &
                  fieldName, &
                  idim_type, &
                  zStart, zEnd, &
                  xbComm(iNeigh)-x0, xeComm(iNeigh)-x0, &
                  ybComm(iNeigh)-y0, yeComm(iNeigh)-y0)
          case default
             write(c0,"(i8)") idim_type
             call fatal_error(h//" idim_type ("//trim(adjustl(c0))//&
                  ") incompatible with a 1D field, for field "//trim(adjustl(fieldName)))
          end select
          call AppendFieldSectionToMessageData(oneFieldSection, Msgs%msgData(cntMsg))
       end if
    end do
  end subroutine InsertFieldSectionAtMessageSet_1D


  


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
                  ") incompatible with a 2D field, for field "//trim(adjustl(fieldName)))
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
                  ") incompatible with a 3D field, for field "//trim(adjustl(fieldName)))
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
                  ") incompatible with a 4D field, for field "//trim(adjustl(fieldName)))
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




  subroutine InsertFieldSectionAtSendRecvMessageSet_1D(&
       field, fieldName, idim_type, myNum, nNeigh, GlobalWithGhost, &
       zbSend, zeSend, zbRecv, zeRecv, &
       xbSend, xeSend, ybSend, yeSend, willSend, SendMessageSet, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMessageSet, &
       kMax)

    ! Inserts a section of a field to be communicated
    ! on a MessageSet variable

    ! field memory address

    real, pointer, intent(in) :: field(:)

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
    integer, intent(in) :: zbSend
    integer, intent(in) :: zeSend
    integer, intent(in) :: zbRecv
    integer, intent(in) :: zeRecv
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
    integer, intent(in) :: kMax

    ! recv message set
    type(MessageSet), pointer, intent(inout) :: RecvMessageSet

    character(len=*), parameter :: h="**(InsertFieldSectionAtSendRecvMessageSet_1D)**"

    ! check arguments

    if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" GlobalWithGhost not associated")
    end if

    ! include field on field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
         zbSend, zeSend, xbSend, xeSend, ybSend, yeSend, willSend, &
         SendMessageSet, kMax)
    call InsertFieldSectionAtMessageSet(&
         myNum, field, fieldName, idim_type, nNeigh, GlobalWithGhost, &
         zbRecv, zeRecv, xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
         RecvMessageSet, kMax)
  end subroutine InsertFieldSectionAtSendRecvMessageSet_1D




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




  subroutine CreateAcoustNewOneMessageSet3D(&
       field, fieldName, idim_type,  &
       ParEnv, Neigh, GlobalOwn, GlobalWithGhost, &
       Tag, NameSend, NameRecv, &
       AcoustNewOneSend, AcoustNewOneRecv)

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
    type(MessageSet), pointer, intent(inout) :: AcoustNewOneSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewOneRecv


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

    character(len=*), parameter :: h="**(CreateAcoustNewOneMessageSet3D)**"
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
       AcoustNewOneSend => null()
       AcoustNewOneRecv => null()
       return
    end if

    if (dumpLocal) then
       call MsgDump(h//" will create AcoustNewOneSend/Recv")
    end if

    myNum  = ParEnv%myNum
    nMachs = ParEnv%nMachs
    nNeigh = Neigh%nNeigh

    ! AcoustNewOneSend, AcoustNewOneRecv:
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
         SendMessageSet=AcoustNewOneSend, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AcoustNewOneRecv)

    !

    call InsertFieldSectionAtSendRecvMessageSet(&
         field, fieldName, idim_type, myNum, nNeigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcoustNewOneSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcoustNewOneRecv)

    if (dumpLocal) then
       call MsgDump(h//" finishes with AcoustNewOneSend MessageSet:")
       call DumpMessageSet(AcoustNewOneSend)
       call MsgDump(h//" finishes with AcoustNewOneRecv MessageSet:")
       call DumpMessageSet(AcoustNewOneRecv)
    end if
  end subroutine CreateAcoustNewOneMessageSet3D





  subroutine CreateAcoustNewOneMessageSet1D(&
       field, fieldName, idim_type,  &
       ParEnv, Neigh, GlobalOwn, GlobalWithGhost, &
       Tag, NameSend, NameRecv, &
       AcoustNewOneSend, AcoustNewOneRecv, kMax)

    real, pointer, intent(in) :: field(:)
    character(len=*), intent(in) :: fieldName
    integer, intent(in) :: idim_type
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    integer, intent(in) :: Tag
    character(len=*), intent(in) :: NameSend
    character(len=*), intent(in) :: NameRecv
    type(MessageSet), pointer, intent(inout) :: AcoustNewOneSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewOneRecv
    integer, intent(in) :: kMax


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

    character(len=*), parameter :: h="**(CreateAcoustNewOneMessageSet1D)**"
    logical, parameter :: dumpLocal=.false.

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" starts with null GlobalOwn")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    if (dumpLocal) then
       call MsgDump(h//" for field "//trim(adjustl(fieldName))//&
            " will create "//trim(adjustl(NameSend))//" and "//&
            trim(adjustl(NameRecv)))
    end if

    ! default output (case no neighbours)

    if (.not. associated(Neigh)) then
       if (dumpLocal) then
          call MsgDump(h//" no neighbours for this Message Set")
       end if
       AcoustNewOneSend => null()
       AcoustNewOneRecv => null()
       return
    end if

    myNum  = ParEnv%myNum
    nMachs = ParEnv%nMachs
    nNeigh = Neigh%nNeigh

    ! AcoustNewOneSend, AcoustNewOneRecv:
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
         SendMessageSet=AcoustNewOneSend, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         RecvMessageSet=AcoustNewOneRecv)

    !

    call InsertFieldSectionAtSendRecvMessageSet(&
         field, fieldName, idim_type, myNum, nNeigh, GlobalWithGhost, &
         1, kMax, 1, kMax, xbSend, xeSend, ybSend, yeSend, willSend, AcoustNewOneSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcoustNewOneRecv, kMax)

    if (dumpLocal) then
       call MsgDump(h//" finishes with AcoustNewOneSend MessageSet:")
       call DumpMessageSet(AcoustNewOneSend)
       call MsgDump(h//" finishes with AcoustNewOneRecv MessageSet:")
       call DumpMessageSet(AcoustNewOneRecv)
    end if
  end subroutine CreateAcoustNewOneMessageSet1D  





  
  subroutine CreateWideGhostZoneMessageSet(&
       ParEnv, Neigh, &
       GlobalOwnWithBC, GlobalWithGhost, mzp, zStartAdvMnt, zEndAdvMnt, &
       TagWideGhostZone, WideGhostZoneSend, WideGhostZoneRecv)

    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    integer, intent(in) :: mzp
    integer, intent(in) :: zStartAdvMnt
    integer, intent(in) :: zendAdvMnt
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

    character(len=*), parameter :: scpName="SCP"
    character(len=*), parameter :: ufxName="UFX"
    character(len=*), parameter :: vfxName="VFX"
    character(len=*), parameter :: wfxName="WFX"
    character(len=*), parameter :: scrName="SCR"
    character(len=*), parameter :: ufxLocalName="UFXLOC"
    character(len=*), parameter :: vfxLocalName="VFXLOC"
    character(len=*), parameter :: wfxLocalName="WFXLOC"

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
                  (mzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scpName, idim_type, &
                  1, mzp, &
                  xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                  ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  1, mzp, &
                  xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                  ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  1, mzp, &
                  xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                  ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  1, mzp, &
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
                  (mzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scpName, idim_type, &
                  1, mzp, &
                  xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                  ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  1, mzp, &
                  xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                  ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  1, mzp, &
                  xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                  ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  1, mzp, &
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
                  (mzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scpName, idim_type, &
                  1, mzp, &
                  xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                  ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  1, mzp, &
                  xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                  ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  1, mzp, &
                  xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                  ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  1, mzp, &
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
                  (mzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scpName, idim_type, &
                  1, mzp, &
                  xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                  ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxName, idim_type, &
                  1, mzp, &
                  xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                  ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxName, idim_type, &
                  1, mzp, &
                  xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                  ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneSend%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxName, idim_type, &
                  1, mzp, &
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
                  (mzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                  ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                  ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                  ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
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
                  (mzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                  ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                  ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                  ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
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
                  (mzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                  ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                  ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                  ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
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
                  (mzp + 2* ghostZoneWidth)

             oneFieldSection =>  CreateFieldSection(scrName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                  ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(ufxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                  ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(vfxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
                  xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                  ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                  fieldSectionSize)
             call AppendFieldSectionToMessageData(oneFieldSection, WideGhostZoneRecv%msgData(cntMsg))
             oneFieldSection =>  CreateFieldSection(wfxLocalName, idim_type, &
                  zStartAdvMnt, zEndAdvMnt, &
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


  subroutine UpdateSendFieldAdressAtAdvectcRk(&
       WideGhostZoneSend, &
       scp, ufx, vfx, wfx)
    type(MessageSet), pointer, intent(in) :: WideGhostZoneSend
    real, pointer, intent(in) :: scp(:,:,:)
    real, pointer, intent(in) :: ufx(:,:,:)
    real, pointer, intent(in) :: vfx(:,:,:)
    real, pointer, intent(in) :: wfx(:,:,:)

    integer :: iMsg
    type(FieldSectionNode), pointer :: fsnode
    character(len=*), parameter :: h="**(UpdateSendFieldAdressAtAdvectcRk)**"
    logical, parameter :: dumpLocal=.false.
    
    if (.not. associated(WideGhostZoneSend)) then
       call fatal_error(h//" WideGhostZoneSend not associated")
    end if

    if (dumpLocal) then
       call MsgDump(h//" of message set "//trim(adjustl(WideGhostZoneSend%name)))
       call MsgDump(h//" of field SCP")
       call MsgDump(h//" of field UFX")
       call MsgDump(h//" of field VFX")
       call MsgDump(h//" of field WFX")
    end if

    do iMsg = 1, WideGhostZoneSend%nMsgs
       fsnode => WideGhostZoneSend%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, scp, "SCP")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, ufx, "UFX")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, vfx, "VFX")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, wfx, "WFX")
    end do
  end subroutine UpdateSendFieldAdressAtAdvectcRk




  subroutine UpdateRecvFieldAdressAtAdvectcRk(&
       WideGhostZoneRecv, &
       scr, ufx_local, vfx_local, wfx_local)
    type(MessageSet), pointer, intent(in) :: WideGhostZoneRecv
    real, pointer, intent(in) :: scr(:,:,:)
    real, pointer, intent(in) :: ufx_local(:,:,:)
    real, pointer, intent(in) :: vfx_local(:,:,:)
    real, pointer, intent(in) :: wfx_local(:,:,:)

    integer :: iMsg
    type(FieldSectionNode), pointer :: fsnode
    character(len=*), parameter :: h="**(UpdateRecvFieldAdressAtAdvectcRk)**"
    logical, parameter :: dumpLocal=.false.
    
    if (.not. associated(WideGhostZoneRecv)) then
       call fatal_error(h//" WideGhostZoneRecv not associated")
    end if
    
    if (dumpLocal) then
       call MsgDump(h//" of message set "//trim(adjustl(WideGhostZoneRecv%name)))
       call MsgDump(h//" of field SCR")
       call MsgDump(h//" of field UFXLOC")
       call MsgDump(h//" of field VFXLOC")
       call MsgDump(h//" of field WFXLOC")
    end if

    do iMsg = 1, WideGhostZoneRecv%nMsgs
       fsnode => WideGhostZoneRecv%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, scr, "SCR")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, ufx_local, "UFXLOC")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, vfx_local, "VFXLOC")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, wfx_local, "WFXLOC")
    end do
  end subroutine UpdateRecvFieldAdressAtAdvectcRk
  


  subroutine OneAdvMntSendRecv(TwoD, Neigh, nNeigh, &
       willSendNorth, willSendSouth, willSendEast, willSendWest, &
       willRecvNorth, willRecvSouth, willRecvEast, willRecvWest, &
       NameSendX, NameSendY, sendDirection, TagX, &
       NameRecvX, NameRecvY, recvDirection, TagY, &
       zStartAdvMnt, zEndAdvMnt, &
       xbSendEast, xeSendEast, xbSendWest, xeSendWest, &
       xbSendNorth, xeSendNorth, xbSendSouth, xeSendSouth, &
       xbRecvEast, xeRecvEast, xbRecvWest, xeRecvWest, &
       xbRecvNorth, xeRecvNorth, xbRecvSouth, xeRecvSouth, &
       ybSendEast, yeSendEast, ybSendWest, yeSendWest, &
       ybSendNorth, yeSendNorth, ybSendSouth, yeSendSouth, &
       ybRecvEast, yeRecvEast, ybRecvWest, yeRecvWest, &
       ybRecvNorth, yeRecvNorth, ybRecvSouth, yeRecvSouth, &
       x0, y0, mzp, ghostZoneWidth, idim_type, &
       AdvMntSendX, AdvMntRecvX, AdvMntSendY, AdvMntRecvY, &
       fldName_1, fldName_2, fldName_3, fldName_4, &
       fldName_X, fldName_Y)

    logical, intent(in) :: TwoD
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    integer, intent(in) :: nNeigh
    logical, intent(in) :: willSendNorth(:)
    logical, intent(in) :: willSendSouth(:)
    logical, intent(in) :: willSendEast(:)
    logical, intent(in) :: willSendWest(:)
    logical, intent(in) :: willRecvNorth(:)
    logical, intent(in) :: willRecvSouth(:)
    logical, intent(in) :: willRecvEast(:)
    logical, intent(in) :: willRecvWest(:)
    character(len=*), intent(in) :: NameSendX
    character(len=*), intent(in) :: NameSendY
    character(len=*), intent(in) :: sendDirection
    integer, intent(in) :: TagX
    character(len=*), intent(in) :: NameRecvX
    character(len=*), intent(in) :: NameRecvY
    character(len=*), intent(in) :: recvDirection
    integer, intent(in) :: TagY
    integer, intent(in) :: zStartAdvMnt
    integer, intent(in) :: zEndAdvMnt
    integer, intent(in) :: xbSendEast(:)
    integer, intent(in) :: xeSendEast(:)
    integer, intent(in) :: xbSendWest(:)
    integer, intent(in) :: xeSendWest(:)
    integer, intent(in) :: xbSendNorth(:)
    integer, intent(in) :: xeSendNorth(:)
    integer, intent(in) :: xbSendSouth(:)
    integer, intent(in) :: xeSendSouth(:)
    integer, intent(in) :: xbRecvEast(:)
    integer, intent(in) :: xeRecvEast(:)
    integer, intent(in) :: xbRecvWest(:)
    integer, intent(in) :: xeRecvWest(:)
    integer, intent(in) :: xbRecvNorth(:)
    integer, intent(in) :: xeRecvNorth(:)
    integer, intent(in) :: xbRecvSouth(:)
    integer, intent(in) :: xeRecvSouth(:)

    integer, intent(in) :: ybSendEast(:)
    integer, intent(in) :: yeSendEast(:)
    integer, intent(in) :: ybSendWest(:)
    integer, intent(in) :: yeSendWest(:)
    integer, intent(in) :: ybSendNorth(:)
    integer, intent(in) :: yeSendNorth(:)
    integer, intent(in) :: ybSendSouth(:)
    integer, intent(in) :: yeSendSouth(:)
    integer, intent(in) :: ybRecvEast(:)
    integer, intent(in) :: yeRecvEast(:)
    integer, intent(in) :: ybRecvWest(:)
    integer, intent(in) :: yeRecvWest(:)
    integer, intent(in) :: ybRecvNorth(:)
    integer, intent(in) :: yeRecvNorth(:)
    integer, intent(in) :: ybRecvSouth(:)
    integer, intent(in) :: yeRecvSouth(:)

    integer, intent(in) :: x0
    integer, intent(in) :: y0
    integer, intent(in) :: mzp
    integer, intent(in) :: ghostZoneWidth
    integer, intent(in) :: idim_type
    type(MessageSet), pointer, intent(out) :: AdvMntSendX
    type(MessageSet), pointer, intent(out) :: AdvMntRecvX
    type(MessageSet), pointer, intent(out) :: AdvMntSendY
    type(MessageSet), pointer, intent(out) :: AdvMntRecvY
    character(len=*), intent(in), optional :: fldName_1
    character(len=*), intent(in), optional :: fldName_2
    character(len=*), intent(in), optional :: fldName_3
    character(len=*), intent(in), optional :: fldName_4
    character(len=*), intent(in), optional :: fldName_X
    character(len=*), intent(in), optional :: fldName_Y

    integer :: nMsgs
    integer :: cntMsg
    integer :: iNeigh
    integer :: fieldSectionSizeZ
    integer :: fieldSectionSize
    logical :: willSend(size(willSendWest))
    logical :: willRecv(size(willSendWest))
    type(FieldSection), pointer :: oneFieldSection

    character(len=8) :: str(10)
    character(len=128) :: strLong
    character(len=*), parameter :: h="**(OneAdvMntSendRecv)**"
    logical, parameter :: dumpLocal=.false.

    ! wrong call if present(fldName_X) and any of fldName_nbr

    if (present(fldName_X) .and. (&
         present(fldName_1) .or. &
         present(fldName_2) .or. &
         present(fldName_3) .or. &
         present(fldName_4))) then
       call fatal_error(h//" inconsistent call: fldName_X cannot be present "//&
            "with any of fldName_1, fldName_2, fldName_3, fldName_4")
    end if
    if (present(fldName_Y) .and. (&
         present(fldName_1) .or. &
         present(fldName_2) .or. &
         present(fldName_3) .or. &
         present(fldName_4))) then
       call fatal_error(h//" inconsistent call: fldName_Y cannot be present "//&
            "with any of fldName_1, fldName_2, fldName_3, fldName_4")
    end if


    if (dumpLocal) then
       call MsgDump(h//" building "//&
            trim(NameSendX)//", "//&
            trim(NameRecvX)//", "//&
            trim(NameSendY)//", "//&
            trim(NameRecvY))
       strLong=h//" invoked with"
       if (present(fldName_1)) then
          strLong=trim(strLong)//" fldName_1="//trim(fldName_1)
       end if
       if (present(fldName_2)) then
          strLong=trim(strLong)//" fldName_2="//trim(fldName_2)
       end if
       if (present(fldName_3)) then
          strLong=trim(strLong)//" fldName_3="//trim(fldName_3)
       end if
       if (present(fldName_4)) then
          strLong=trim(strLong)//" fldName_4="//trim(fldName_4)
       end if
       if (present(fldName_X)) then
          strLong=trim(strLong)//" fldName_X="//trim(fldName_X)
       end if
       if (present(fldName_Y)) then
          strLong=trim(strLong)//" fldName_Y="//trim(fldName_Y)
       end if
       call MsgDump(trim(strLong))
    end if
    
    ! z direction field section size (to accomodate 2D fields)

    if (TwoD) then
       fieldSectionSizeZ = 1
    else
       fieldSectionSizeZ = mzp + 2* ghostZoneWidth
    end if
    
    ! send message set will contain sends for east-west directions

    willSend = willSendEast .or. willSendWest

    ! create message set for east-west sends

    AdvMntSendX => CreateMessageSet(&
         NameSendX, &
         sendDirection, &
         TagX, &
         willSend, &
         Neigh)

    ! insert field sections named u3d, v3d at each direction
    ! with null field addresses, to be updated whenever
    ! real addresses are known

    if ( associated(AdvMntSendX)) then
       nMsgs = AdvMntSendX%nMsgs

       ! create list of Field Sections to send, one for
       ! each process to communicate and insert at the send MessageSet
       ! field section list

       ! since there is at most one neighbour node at each direction,
       ! there will be at most one MessageSet at each direction

       cntMsg = 0
       do iNeigh = 1, nNeigh
          if (willSendEast(iNeigh)) then

             ! insert send communications to east

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AdvMntSendX%name)))
             end if

             fieldSectionSize= &
                  (xeSendEast(iNeigh)-xbSendEast(iNeigh)+1) * &
                  (yeSendEast(iNeigh)-ybSendEast(iNeigh)+1) * &
                  fieldSectionSizeZ

             if (present(fldName_X)) then
                oneFieldSection =>  CreateFieldSection(fldName_X, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                     ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if
             if (present(fldName_1)) then
                oneFieldSection =>  CreateFieldSection(fldName_1, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                     ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if
             if (present(fldName_2)) then
                oneFieldSection =>  CreateFieldSection(fldName_2, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                     ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if
             if (present(fldName_3)) then
                oneFieldSection =>  CreateFieldSection(fldName_3, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                     ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if
             if (present(fldName_4)) then
                oneFieldSection =>  CreateFieldSection(fldName_4, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendEast(iNeigh)-x0, xeSendEast(iNeigh)-x0, &
                     ybSendEast(iNeigh)-y0, yeSendEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if

          else if (willSendWest(iNeigh)) then

             ! insert send communications to west

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AdvMntSendX%name)))
             end if

             fieldSectionSize= &
                  (xeSendWest(iNeigh)-xbSendWest(iNeigh)+1) * &
                  (yeSendWest(iNeigh)-ybSendWest(iNeigh)+1) * &
                  fieldSectionSizeZ

             if (present(fldName_X)) then
                oneFieldSection =>  CreateFieldSection(fldName_X, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                     ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if
             if (present(fldName_1)) then
                oneFieldSection =>  CreateFieldSection(fldName_1, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                     ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if
             if (present(fldName_2)) then
                oneFieldSection =>  CreateFieldSection(fldName_2, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                     ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if
             if (present(fldName_3)) then
                oneFieldSection =>  CreateFieldSection(fldName_3, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                     ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if
             if (present(fldName_4)) then
                oneFieldSection =>  CreateFieldSection(fldName_4, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendWest(iNeigh)-x0, xeSendWest(iNeigh)-x0, &
                     ybSendWest(iNeigh)-y0, yeSendWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendX%msgData(cntMsg))
             end if
          end if

       end do
    end if

    ! recv message set will contain recvs for east-west directions

    willRecv =  willRecvEast .or. willRecvWest

    ! create message set for  east-west recvs

    AdvMntRecvX => CreateMessageSet(&
         NameRecvX, &
         recvDirection, &
         TagX, &
         willRecv, &
         Neigh)

    if ( associated(AdvMntRecvX)) then
       nMsgs = AdvMntRecvX%nMsgs
       cntMsg = 0
       do iNeigh = 1, nNeigh
          if (willRecvEast(iNeigh)) then

             ! insert recv communications from east

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AdvMntRecvX%name)))
             end if

             fieldSectionSize= &
                  (xeRecvEast(iNeigh)-xbRecvEast(iNeigh)+1) * &
                  (yeRecvEast(iNeigh)-ybRecvEast(iNeigh)+1) * &
                  fieldSectionSizeZ
             if (present(fldName_X)) then
                oneFieldSection =>  CreateFieldSection(fldName_X, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                     ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if
             if (present(fldName_1)) then
                oneFieldSection =>  CreateFieldSection(fldName_1, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                     ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if
             if (present(fldName_2)) then
                oneFieldSection =>  CreateFieldSection(fldName_2, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                     ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if
             if (present(fldName_3)) then
                oneFieldSection =>  CreateFieldSection(fldName_3, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                     ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if
             if (present(fldName_4)) then
                oneFieldSection =>  CreateFieldSection(fldName_4, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvEast(iNeigh)-x0, xeRecvEast(iNeigh)-x0, &
                     ybRecvEast(iNeigh)-y0, yeRecvEast(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if

          else if (willRecvWest(iNeigh)) then

             ! insert recv communications from west

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AdvMntRecvX%name)))
             end if

             fieldSectionSize= &
                  (xeRecvWest(iNeigh)-xbRecvWest(iNeigh)+1) * &
                  (yeRecvWest(iNeigh)-ybRecvWest(iNeigh)+1) * &
                  fieldSectionSizeZ
             if (present(fldName_X)) then
                oneFieldSection =>  CreateFieldSection(fldName_X, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                     ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if
             if (present(fldName_1)) then
                oneFieldSection =>  CreateFieldSection(fldName_1, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                     ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if
             if (present(fldName_2)) then
                oneFieldSection =>  CreateFieldSection(fldName_2, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                     ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if
             if (present(fldName_3)) then
                oneFieldSection =>  CreateFieldSection(fldName_3, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                     ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if
             if (present(fldName_4)) then
                oneFieldSection =>  CreateFieldSection(fldName_4, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvWest(iNeigh)-x0, xeRecvWest(iNeigh)-x0, &
                     ybRecvWest(iNeigh)-y0, yeRecvWest(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvX%msgData(cntMsg))
             end if
          end if
       end do
    end if

    ! send message set will contain sends for north/south directions

    willSend = willSendNorth .or. willSendSouth

    ! create message set for  north-south sends

    AdvMntSendY => CreateMessageSet(&
         NameSendY, &
         sendDirection, &
         TagY, &
         willSend, &
         Neigh)

    ! insert field sections named u3d, v3d at each direction
    ! with null field addresses, to be updated whenever
    ! real addresses are known

    if ( associated(AdvMntSendY)) then
       nMsgs = AdvMntSendY%nMsgs

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
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AdvMntSendY%name)))
             end if

             fieldSectionSize= &
                  (xeSendNorth(iNeigh)-xbSendNorth(iNeigh)+1) * &
                  (yeSendNorth(iNeigh)-ybSendNorth(iNeigh)+1) * &
                  fieldSectionSizeZ
             if (present(fldName_Y)) then
                oneFieldSection =>  CreateFieldSection(fldName_Y, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                     ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if
             if (present(fldName_1)) then
                oneFieldSection =>  CreateFieldSection(fldName_1, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                     ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if
             if (present(fldName_2)) then
                oneFieldSection =>  CreateFieldSection(fldName_2, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                     ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if
             if (present(fldName_3)) then
                oneFieldSection =>  CreateFieldSection(fldName_3, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                     ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if
             if (present(fldName_4)) then
                oneFieldSection =>  CreateFieldSection(fldName_4, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendNorth(iNeigh)-x0, xeSendNorth(iNeigh)-x0, &
                     ybSendNorth(iNeigh)-y0, yeSendNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if

          else if (willSendSouth(iNeigh)) then

             ! insert send communications to south

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AdvMntSendY%name)))
             end if

             fieldSectionSize= &
                  (xeSendSouth(iNeigh)-xbSendSouth(iNeigh)+1) * &
                  (yeSendSouth(iNeigh)-ybSendSouth(iNeigh)+1) * &
                  fieldSectionSizeZ
             if (present(fldName_Y)) then
                oneFieldSection =>  CreateFieldSection(fldName_Y, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                     ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if
             if (present(fldName_1)) then
                oneFieldSection =>  CreateFieldSection(fldName_1, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                     ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if
             if (present(fldName_2)) then
                oneFieldSection =>  CreateFieldSection(fldName_2, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                     ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if
             if (present(fldName_3)) then
                oneFieldSection =>  CreateFieldSection(fldName_3, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                     ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if
             if (present(fldName_4)) then
                oneFieldSection =>  CreateFieldSection(fldName_4, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbSendSouth(iNeigh)-x0, xeSendSouth(iNeigh)-x0, &
                     ybSendSouth(iNeigh)-y0, yeSendSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntSendY%msgData(cntMsg))
             end if
          end if
       end do
    end if

    ! recv message set will contain recvs for north/south directions

    willRecv = willRecvNorth .or. willRecvSouth

    ! create message set for  north-south recvs

    AdvMntRecvY => CreateMessageSet(&
         NameRecvY, &
         recvDirection, &
         TagY, &
         willRecv, &
         Neigh)

    if ( associated(AdvMntRecvY)) then
       nMsgs = AdvMntRecvY%nMsgs
       cntMsg = 0
       do iNeigh = 1, nNeigh
          if (willRecvNorth(iNeigh)) then

             ! insert recv communications from north

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AdvMntRecvY%name)))
             end if

             fieldSectionSize= &
                  (xeRecvNorth(iNeigh)-xbRecvNorth(iNeigh)+1) * &
                  (yeRecvNorth(iNeigh)-ybRecvNorth(iNeigh)+1) * &
                  fieldSectionSizeZ
             if (present(fldName_Y)) then
                oneFieldSection =>  CreateFieldSection(fldName_Y, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                     ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if
             if (present(fldName_1)) then
                oneFieldSection =>  CreateFieldSection(fldName_1, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                     ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if
             if (present(fldName_2)) then
                oneFieldSection =>  CreateFieldSection(fldName_2, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                     ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if
             if (present(fldName_3)) then
                oneFieldSection =>  CreateFieldSection(fldName_3, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                     ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if
             if (present(fldName_4)) then
                oneFieldSection =>  CreateFieldSection(fldName_4, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvNorth(iNeigh)-x0, xeRecvNorth(iNeigh)-x0, &
                     ybRecvNorth(iNeigh)-y0, yeRecvNorth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if

          else if (willRecvSouth(iNeigh)) then

             ! insert recv communications from south

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AdvMntRecvY%name)))
             end if

             fieldSectionSize= &
                  (xeRecvSouth(iNeigh)-xbRecvSouth(iNeigh)+1) * &
                  (yeRecvSouth(iNeigh)-ybRecvSouth(iNeigh)+1) * &
                  fieldSectionSizeZ
             if (present(fldName_Y)) then
                oneFieldSection =>  CreateFieldSection(fldName_Y, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                     ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if
             if (present(fldName_1)) then
                oneFieldSection =>  CreateFieldSection(fldName_1, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                     ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if
             if (present(fldName_2)) then
                oneFieldSection =>  CreateFieldSection(fldName_2, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                     ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if
             if (present(fldName_3)) then
                oneFieldSection =>  CreateFieldSection(fldName_3, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                     ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if
             if (present(fldName_4)) then
                oneFieldSection =>  CreateFieldSection(fldName_4, idim_type, &
                     zStartAdvMnt, zEndAdvMnt, &
                     xbRecvSouth(iNeigh)-x0, xeRecvSouth(iNeigh)-x0, &
                     ybRecvSouth(iNeigh)-y0, yeRecvSouth(iNeigh)-y0, &
                     fieldSectionSize)
                call AppendFieldSectionToMessageData(oneFieldSection, AdvMntRecvY%msgData(cntMsg))
             end if
          end if
       end do
    end if

    if (dumpLocal) then
       call MsgDump(h//" built "//trim(NameSendX)//" MessageSet:")
       call DumpMessageSet(AdvMntSendX)
       call MsgDump(h//" built "//trim(NameRecvX)//" MessageSet:")
       call DumpMessageSet(AdvMntRecvX)
       call MsgDump(h//" built "//trim(NameSendY)//" MessageSet:")
       call DumpMessageSet(AdvMntSendY)
       call MsgDump(h//" built "//trim(NameRecvY)//" MessageSet:")
       call DumpMessageSet(AdvMntRecvY)
    end if
  end subroutine OneAdvMntSendRecv
  
  



  
  subroutine CreateAdvMntMessageSet(&
       ParEnv, Neigh, &
       GlobalOwnWithBC, GlobalWithGhostAdvMnt, &
       NodeDims, NodeDimsAdvMnt, &
       TagAdvMntUVX, AdvMntUVSendX, AdvMntUVRecvX, &
       TagAdvMntUVY, AdvMntUVSendY, AdvMntUVRecvY, &
       TagAdvMntDxDyX, AdvMntDxDySendX, AdvMntDxDyRecvX, &
       TagAdvMntDxDyY, AdvMntDxDySendY, AdvMntDxDyRecvY, &
       TagAdvMntDd0X, AdvMntDd0SendX, AdvMntDd0RecvX, &
       TagAdvMntDd0Y, AdvMntDd0SendY, AdvMntDd0RecvY, &
       TagAdvMntDenX, AdvMntDenSendX, AdvMntDenRecvX, &
       TagAdvMntDenY, AdvMntDenSendY, AdvMntDenRecvY, &
       TagAdvMntScaX, AdvMntScaSendX, AdvMntScaRecvX, &
       TagAdvMntScaY, AdvMntScaSendY, AdvMntScaRecvY)

    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwnWithBC
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhostAdvMnt
    type(NodeDimensions), pointer, intent(in) :: NodeDims
    type(NodeDimensions), pointer, intent(in) :: NodeDimsAdvMnt
    integer, intent(in) :: TagAdvMntUVX
    type(MessageSet), pointer, intent(inout) :: AdvMntUVSendX
    type(MessageSet), pointer, intent(inout) :: AdvMntUVRecvX
    integer, intent(in) :: TagAdvMntUVY
    type(MessageSet), pointer, intent(inout) :: AdvMntUVSendY
    type(MessageSet), pointer, intent(inout) :: AdvMntUVRecvY
    integer, intent(in) :: TagAdvMntDxDyX
    type(MessageSet), pointer, intent(inout) :: AdvMntDxDySendX
    type(MessageSet), pointer, intent(inout) :: AdvMntDxDyRecvX
    integer, intent(in) :: TagAdvMntDxDyY
    type(MessageSet), pointer, intent(inout) :: AdvMntDxDySendY
    type(MessageSet), pointer, intent(inout) :: AdvMntDxDyRecvY
    integer, intent(in) :: TagAdvMntDd0X
    type(MessageSet), pointer, intent(inout) :: AdvMntDd0SendX
    type(MessageSet), pointer, intent(inout) :: AdvMntDd0RecvX
    integer, intent(in) :: TagAdvMntDd0Y
    type(MessageSet), pointer, intent(inout) :: AdvMntDd0SendY
    type(MessageSet), pointer, intent(inout) :: AdvMntDd0RecvY
    integer, intent(in) :: TagAdvMntDenX
    type(MessageSet), pointer, intent(inout) :: AdvMntDenSendX
    type(MessageSet), pointer, intent(inout) :: AdvMntDenRecvX
    integer, intent(in) :: TagAdvMntDenY
    type(MessageSet), pointer, intent(inout) :: AdvMntDenSendY
    type(MessageSet), pointer, intent(inout) :: AdvMntDenRecvY
    integer, intent(in) :: TagAdvMntScaX
    type(MessageSet), pointer, intent(inout) :: AdvMntScaSendX
    type(MessageSet), pointer, intent(inout) :: AdvMntScaRecvX
    integer, intent(in) :: TagAdvMntScaY
    type(MessageSet), pointer, intent(inout) :: AdvMntScaSendY
    type(MessageSet), pointer, intent(inout) :: AdvMntScaRecvY

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

    character(len=*), parameter :: NameSendUVX="AdvMntUVSendX"
    character(len=*), parameter :: NameRecvUVX="AdvMntUVRecvX"
    character(len=*), parameter :: NameSendUVY="AdvMntUVSendY"
    character(len=*), parameter :: NameRecvUVY="AdvMntUVRecvY"
    character(len=*), parameter :: NameSendDxDyX="AdvMntDxDySendX"
    character(len=*), parameter :: NameRecvDxDyX="AdvMntDxDyRecvX"
    character(len=*), parameter :: NameSendDxDyY="AdvMntDxDySendY"
    character(len=*), parameter :: NameRecvDxDyY="AdvMntDxDyRecvY"
    character(len=*), parameter :: NameSendDd0X="AdvMntDd0SendX"
    character(len=*), parameter :: NameRecvDd0X="AdvMntDd0RecvX"
    character(len=*), parameter :: NameSendDd0Y="AdvMntDd0SendY"
    character(len=*), parameter :: NameRecvDd0Y="AdvMntDd0RecvY"
    character(len=*), parameter :: NameSendDenX="AdvMntDenSendX"
    character(len=*), parameter :: NameRecvDenX="AdvMntDenRecvX"
    character(len=*), parameter :: NameSendDenY="AdvMntDenSendY"
    character(len=*), parameter :: NameRecvDenY="AdvMntDenRecvY"
    character(len=*), parameter :: NameSendScaX="AdvMntScaSendX"
    character(len=*), parameter :: NameRecvScaX="AdvMntScaRecvX"
    character(len=*), parameter :: NameSendScaY="AdvMntScaSendY"
    character(len=*), parameter :: NameRecvScaY="AdvMntScaRecvY"
    character(len=*), parameter :: NameSendRecvNorth="AdvMntSend/RecvNorth"
    character(len=*), parameter :: NameSendRecvSouth="AdvMntSend/RecvSouth"
    character(len=*), parameter :: NameSendRecvEast="AdvMntSend/RecvEast"
    character(len=*), parameter :: NameSendRecvWest="AdvMntSend/RecvWest"

    character(len=*), parameter :: u3dName="U3D"
    character(len=*), parameter :: v3dName="V3D"
    character(len=*), parameter :: dxName="DXTW"
    character(len=*), parameter :: dyName="DYTW"
    character(len=*), parameter :: dd0_3dName="DDO_3D"
    character(len=*), parameter :: dd0_3duName="DDO_3DU"
    character(len=*), parameter :: dd0_3dvName="DDO_3DV"
    character(len=*), parameter :: dd0_3dwName="DDO_3DW"
    character(len=*), parameter :: den0_3dName="DEN0_3D"
    character(len=*), parameter :: den1_3dName="DEN1_3D"
    character(len=*), parameter :: den2_3dName="DEN2_3D"
    character(len=*), parameter :: den3_3dName="DEN3_3D"
    character(len=*), parameter :: vc3d_inName="VC3D_IN"
    character(len=*), parameter :: vc3d_outName="VC3D_OUT"
    
    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    integer :: mzp
    integer :: ierr
    integer :: iNeigh
    integer :: iNode
    integer :: x0
    integer :: y0
    integer :: nMsgs
    integer :: cntMsg
    integer :: fieldSectionSize
    integer :: bramsProcNbr
    integer, parameter :: ghostZoneWidth=3
    integer, parameter :: idim_type_3D=3
    integer, parameter :: idim_type_2D=2
    logical, parameter :: TwoD=.true.
    logical, parameter :: ThreeD=.false.
    type(FieldSection), pointer :: oneFieldSection

    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateAdvMntMessageSet)**"

    nMachs=ParEnv%nmachs
    myNum=ParEnv%mynum
    nNeigh=Neigh%nNeigh
    mzp=NodeDims%mzp

    if (dumpLocal) then
       write(str(1),"(i8)") nMachs
       write(str(2),"(i8)") myNum
       call MsgDump(h//" enter with nMachs="//trim(adjustl(str(1)))//&
            "; myNum="//trim(adjustl(str(2))))
       call DumpNeighbourNodes(Neigh,"AdvMnt")
    end if

    ! offsets to convert global indices to local indices at this proc

    x0 = GlobalWithGhostAdvMnt%xb(myNum) - 1
    y0 = GlobalWithGhostAdvMnt%yb(myNum) - 1

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

    ! extend send/recv region to include ghost zone
    ! this is technically wrong but required to replicate
    ! previous message passsing routines, to be excluded from the code

    do iNeigh = 1, nNeigh
       if (willSendNorth(iNeigh)) then
          bramsProcNbr = Neigh%neigh(iNeigh)
          if (dumpLocal) then
             write(str(1),"(i8)") Brams2MpiProcNbr(bramsProcNbr)
             write(str(2),"(i8)") xbSendNorth(iNeigh)
             write(str(3),"(i8)") xeSendNorth(iNeigh)
          end if
          ! west boundary
          if (.not. NodeDimsAdvMnt%borderWest) then
             xbSendNorth(iNeigh)=xbSendNorth(iNeigh)-ghostZoneWidth
          end if
          ! east boundary
          if (.not. NodeDimsAdvMnt%borderEast) then
             xeSendNorth(iNeigh)=xeSendNorth(iNeigh)+ghostZoneWidth
          end if
          if (dumpLocal) then
             write(str(4),"(i8)") xbSendNorth(iNeigh)
             write(str(5),"(i8)") xeSendNorth(iNeigh)
             call MsgDump(h//" send north to MPI #"//trim(adjustl(str(1)))//&
                  " expanded x interval from "//&
                  "("//trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")"//&
                  " to "//&
                  "("//trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
          end if
       end if
       if (willRecvNorth(iNeigh)) then
          bramsProcNbr = Neigh%neigh(iNeigh)
             write(str(1),"(i8)") Brams2MpiProcNbr(bramsProcNbr)
             write(str(2),"(i8)") xbRecvNorth(iNeigh)
             write(str(3),"(i8)") xeRecvNorth(iNeigh)
          ! west boundary
          if (.not. NodeDimsAdvMnt%borderWest) then
             xbRecvNorth(iNeigh)=xbRecvNorth(iNeigh)-ghostZoneWidth
          end if
          ! east boundary
          if (.not. NodeDimsAdvMnt%borderEast) then
             xeRecvNorth(iNeigh)=xeRecvNorth(iNeigh)+ghostZoneWidth
          end if
          if (dumpLocal) then
             write(str(4),"(i8)") xbRecvNorth(iNeigh)
             write(str(5),"(i8)") xeRecvNorth(iNeigh)
             call MsgDump(h//" recv north from MPI #"//trim(adjustl(str(1)))//&
                  " expanded x interval from "//&
                  "("//trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")"//&
                  " to "//&
                  "("//trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
          end if
       end if
    end do
    
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

    ! extend send/recv region to include ghost zone
    ! this is technically wrong but required to replicate
    ! previous message passsing routines, to be excluded from the code

    do iNeigh = 1, nNeigh
       if (willSendSouth(iNeigh)) then
          bramsProcNbr = Neigh%neigh(iNeigh)
          if (dumpLocal) then
             write(str(1),"(i8)") Brams2MpiProcNbr(bramsProcNbr)
             write(str(2),"(i8)") xbSendSouth(iNeigh)
             write(str(3),"(i8)") xeSendSouth(iNeigh)
          end if
          ! west boundary
          if (.not. NodeDimsAdvMnt%borderWest) then
             xbSendSouth(iNeigh)=xbSendSouth(iNeigh)-ghostZoneWidth
          end if
          ! east boundary
          if (.not. NodeDimsAdvMnt%borderEast) then
             xeSendSouth(iNeigh)=xeSendSouth(iNeigh)+ghostZoneWidth
          end if
          if (dumpLocal) then
             write(str(4),"(i8)") xbSendSouth(iNeigh)
             write(str(5),"(i8)") xeSendSouth(iNeigh)
             call MsgDump(h//" send south to MPI #"//trim(adjustl(str(1)))//&
                  " expanded x interval from "//&
                  "("//trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")"//&
                  " to "//&
                  "("//trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
          end if
       end if
       if (willRecvSouth(iNeigh)) then
          bramsProcNbr = Neigh%neigh(iNeigh)
             write(str(1),"(i8)") Brams2MpiProcNbr(bramsProcNbr)
             write(str(2),"(i8)") xbRecvSouth(iNeigh)
             write(str(3),"(i8)") xeRecvSouth(iNeigh)
          ! west boundary
          if (.not. NodeDimsAdvMnt%borderWest) then
             xbRecvSouth(iNeigh)=xbRecvSouth(iNeigh)-ghostZoneWidth
          end if
          ! east boundary
          if (.not. NodeDimsAdvMnt%borderEast) then
             xeRecvSouth(iNeigh)=xeRecvSouth(iNeigh)+ghostZoneWidth
          end if
          if (dumpLocal) then
             write(str(4),"(i8)") xbRecvSouth(iNeigh)
             write(str(5),"(i8)") xeRecvSouth(iNeigh)
             call MsgDump(h//" recv south from MPI #"//trim(adjustl(str(1)))//&
                  " expanded x interval from "//&
                  "("//trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")"//&
                  " to "//&
                  "("//trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
          end if
       end if
    end do

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

    ! extend send/recv region to include ghost zone
    ! this is technically wrong but required to replicate
    ! previous message passsing routines, to be excluded from the code

    do iNeigh = 1, nNeigh
       if (willSendEast(iNeigh)) then
          bramsProcNbr = Neigh%neigh(iNeigh)
          if (dumpLocal) then
             write(str(1),"(i8)") Brams2MpiProcNbr(bramsProcNbr)
             write(str(2),"(i8)") ybSendEast(iNeigh)
             write(str(3),"(i8)") yeSendEast(iNeigh)
          end if
          ! south boundary
          if (.not. NodeDimsAdvMnt%borderSouth) then
             ybSendEast(iNeigh)=ybSendEast(iNeigh)-ghostZoneWidth
          end if
          ! north boundary
          if (.not. NodeDimsAdvMnt%borderNorth) then
             yeSendEast(iNeigh)=yeSendEast(iNeigh)+ghostZoneWidth
          end if
          if (dumpLocal) then
             write(str(4),"(i8)") ybSendEast(iNeigh)
             write(str(5),"(i8)") yeSendEast(iNeigh)
             call MsgDump(h//" send east to MPI #"//trim(adjustl(str(1)))//&
                  " expanded y interval from "//&
                  "("//trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")"//&
                  " to "//&
                  "("//trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
          end if
       end if
       if (willRecvEast(iNeigh)) then
          bramsProcNbr = Neigh%neigh(iNeigh)
             write(str(1),"(i8)") Brams2MpiProcNbr(bramsProcNbr)
             write(str(2),"(i8)") ybRecvEast(iNeigh)
             write(str(3),"(i8)") yeRecvEast(iNeigh)
          ! south boundary
          if (.not. NodeDimsAdvMnt%borderSouth) then
             ybRecvEast(iNeigh)=ybRecvEast(iNeigh)-ghostZoneWidth
          end if
          ! north boundary
          if (.not. NodeDimsAdvMnt%borderNorth) then
             yeRecvEast(iNeigh)=yeRecvEast(iNeigh)+ghostZoneWidth
          end if
          if (dumpLocal) then
             write(str(4),"(i8)") ybRecvEast(iNeigh)
             write(str(5),"(i8)") yeRecvEast(iNeigh)
             call MsgDump(h//" recv east from MPI #"//trim(adjustl(str(1)))//&
                  " expanded y interval from "//&
                  "("//trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")"//&
                  " to "//&
                  "("//trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
          end if
       end if
    end do

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

    ! extend send/recv region to include ghost zone
    ! this is technically wrong but required to replicate
    ! previous message passsing routines, to be excluded from the code

    do iNeigh = 1, nNeigh
       if (willSendWest(iNeigh)) then
          bramsProcNbr = Neigh%neigh(iNeigh)
          if (dumpLocal) then
             write(str(1),"(i8)") Brams2MpiProcNbr(bramsProcNbr)
             write(str(2),"(i8)") ybSendWest(iNeigh)
             write(str(3),"(i8)") yeSendWest(iNeigh)
          end if
          ! south boundary
          if (.not. NodeDimsAdvMnt%borderSouth) then
             ybSendWest(iNeigh)=ybSendWest(iNeigh)-ghostZoneWidth
          end if
          ! north boundary
          if (.not. NodeDimsAdvMnt%borderNorth) then
             yeSendWest(iNeigh)=yeSendWest(iNeigh)+ghostZoneWidth
          end if
          if (dumpLocal) then
             write(str(4),"(i8)") ybSendWest(iNeigh)
             write(str(5),"(i8)") yeSendWest(iNeigh)
             call MsgDump(h//" send west to MPI #"//trim(adjustl(str(1)))//&
                  " expanded y interval from "//&
                  "("//trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")"//&
                  " to "//&
                  "("//trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
          end if
       end if
       if (willRecvWest(iNeigh)) then
          bramsProcNbr = Neigh%neigh(iNeigh)
             write(str(1),"(i8)") Brams2MpiProcNbr(bramsProcNbr)
             write(str(2),"(i8)") ybRecvWest(iNeigh)
             write(str(3),"(i8)") yeRecvWest(iNeigh)
          ! south boundary
          if (.not. NodeDimsAdvMnt%borderSouth) then
             ybRecvWest(iNeigh)=ybRecvWest(iNeigh)-ghostZoneWidth
          end if
          ! north boundary
          if (.not. NodeDimsAdvMnt%borderNorth) then
             yeRecvWest(iNeigh)=yeRecvWest(iNeigh)+ghostZoneWidth
          end if
          if (dumpLocal) then
             write(str(4),"(i8)") ybRecvWest(iNeigh)
             write(str(5),"(i8)") yeRecvWest(iNeigh)
             call MsgDump(h//" recv west from MPI #"//trim(adjustl(str(1)))//&
                  " expanded y interval from "//&
                  "("//trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")"//&
                  " to "//&
                  "("//trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
          end if
       end if
    end do

    ! create message set for UV 

    call OneAdvMntSendRecv(ThreeD, Neigh, nNeigh, &
       willSendNorth, willSendSouth, willSendEast, willSendWest, &
       willRecvNorth, willRecvSouth, willRecvEast, willRecvWest, &
       NameSendUVX, NameSendUVY, sendDirection, TagAdvMntUVX, &
       NameRecvUVX, NameRecvUVY, recvDirection, TagAdvMntUVY, &
       1, mzp, &
       xbSendEast, xeSendEast, xbSendWest, xeSendWest, &
       xbSendNorth, xeSendNorth, xbSendSouth, xeSendSouth, &
       xbRecvEast, xeRecvEast, xbRecvWest, xeRecvWest, &
       xbRecvNorth, xeRecvNorth, xbRecvSouth, xeRecvSouth, &
       ybSendEast, yeSendEast, ybSendWest, yeSendWest, &
       ybSendNorth, yeSendNorth, ybSendSouth, yeSendSouth, &
       ybRecvEast, yeRecvEast, ybRecvWest, yeRecvWest, &
       ybRecvNorth, yeRecvNorth, ybRecvSouth, yeRecvSouth, &
       x0, y0, NodeDimsAdvMnt%mzp, ghostZoneWidth, idim_type_3D, &
       AdvMntUVSendX, AdvMntUVRecvX, AdvMntUVSendY, AdvMntUVRecvY, &
       fldName_1=u3dName, &
       fldName_2=v3dName)

    ! create message set for Dd0

    call OneAdvMntSendRecv(ThreeD, Neigh, nNeigh, &
       willSendNorth, willSendSouth, willSendEast, willSendWest, &
       willRecvNorth, willRecvSouth, willRecvEast, willRecvWest, &
       NameSendDd0X, NameSendDd0Y, sendDirection, TagAdvMntDd0X, &
       NameRecvDd0X, NameRecvDd0Y, recvDirection, TagAdvMntDd0Y, &
       1, mzp, &
       xbSendEast, xeSendEast, xbSendWest, xeSendWest, &
       xbSendNorth, xeSendNorth, xbSendSouth, xeSendSouth, &
       xbRecvEast, xeRecvEast, xbRecvWest, xeRecvWest, &
       xbRecvNorth, xeRecvNorth, xbRecvSouth, xeRecvSouth, &
       ybSendEast, yeSendEast, ybSendWest, yeSendWest, &
       ybSendNorth, yeSendNorth, ybSendSouth, yeSendSouth, &
       ybRecvEast, yeRecvEast, ybRecvWest, yeRecvWest, &
       ybRecvNorth, yeRecvNorth, ybRecvSouth, yeRecvSouth, &
       x0, y0, NodeDimsAdvMnt%mzp, ghostZoneWidth, idim_type_3D, &
       AdvMntDd0SendX, AdvMntDd0RecvX, AdvMntDd0SendY, AdvMntDd0RecvY, &
       fldName_1=dd0_3dName, &
       fldName_2=dd0_3duName, &
       fldName_3=dd0_3dvName, &
       fldName_4=dd0_3dwName)

    ! create message set for Den

    call OneAdvMntSendRecv(ThreeD, Neigh, nNeigh, &
       willSendNorth, willSendSouth, willSendEast, willSendWest, &
       willRecvNorth, willRecvSouth, willRecvEast, willRecvWest, &
       NameSendDenX, NameSendDenY, sendDirection, TagAdvMntDenX, &
       NameRecvDenX, NameRecvDenY, recvDirection, TagAdvMntDenY, &
       1, mzp, &
       xbSendEast, xeSendEast, xbSendWest, xeSendWest, &
       xbSendNorth, xeSendNorth, xbSendSouth, xeSendSouth, &
       xbRecvEast, xeRecvEast, xbRecvWest, xeRecvWest, &
       xbRecvNorth, xeRecvNorth, xbRecvSouth, xeRecvSouth, &
       ybSendEast, yeSendEast, ybSendWest, yeSendWest, &
       ybSendNorth, yeSendNorth, ybSendSouth, yeSendSouth, &
       ybRecvEast, yeRecvEast, ybRecvWest, yeRecvWest, &
       ybRecvNorth, yeRecvNorth, ybRecvSouth, yeRecvSouth, &
       x0, y0, NodeDimsAdvMnt%mzp, ghostZoneWidth, idim_type_3D, &
       AdvMntDenSendX, AdvMntDenRecvX, AdvMntDenSendY, AdvMntDenRecvY, &
       fldName_1=den0_3dName, &
       fldName_2=den1_3dName, &
       fldName_3=den2_3dName, &
       fldName_4=den3_3dName)

    ! create message set for DxDy

    call OneAdvMntSendRecv(TwoD, Neigh, nNeigh, &
       willSendNorth, willSendSouth, willSendEast, willSendWest, &
       willRecvNorth, willRecvSouth, willRecvEast, willRecvWest, &
       NameSendDxDyX, NameSendDxDyY, sendDirection, TagAdvMntDxDyX, &
       NameRecvDxDyX, NameRecvDxDyY, recvDirection, TagAdvMntDxDyY, &
       1, mzp, &
       xbSendEast, xeSendEast, xbSendWest, xeSendWest, &
       xbSendNorth, xeSendNorth, xbSendSouth, xeSendSouth, &
       xbRecvEast, xeRecvEast, xbRecvWest, xeRecvWest, &
       xbRecvNorth, xeRecvNorth, xbRecvSouth, xeRecvSouth, &
       ybSendEast, yeSendEast, ybSendWest, yeSendWest, &
       ybSendNorth, yeSendNorth, ybSendSouth, yeSendSouth, &
       ybRecvEast, yeRecvEast, ybRecvWest, yeRecvWest, &
       ybRecvNorth, yeRecvNorth, ybRecvSouth, yeRecvSouth, &
       x0, y0, 1, ghostZoneWidth, idim_type_2D, &
       AdvMntDxDySendX, AdvMntDxDyRecvX, AdvMntDxDySendY, AdvMntDxDyRecvY, &
       fldName_1=dxName, &
       fldName_2=dyName)

    ! create message set for Sca

    call OneAdvMntSendRecv(ThreeD, Neigh, nNeigh, &
       willSendNorth, willSendSouth, willSendEast, willSendWest, &
       willRecvNorth, willRecvSouth, willRecvEast, willRecvWest, &
       NameSendScaX, NameSendScaY, sendDirection, TagAdvMntScaX, &
       NameRecvScaX, NameRecvScaY, recvDirection, TagAdvMntScaY, &
       1, mzp, &
       xbSendEast, xeSendEast, xbSendWest, xeSendWest, &
       xbSendNorth, xeSendNorth, xbSendSouth, xeSendSouth, &
       xbRecvEast, xeRecvEast, xbRecvWest, xeRecvWest, &
       xbRecvNorth, xeRecvNorth, xbRecvSouth, xeRecvSouth, &
       ybSendEast, yeSendEast, ybSendWest, yeSendWest, &
       ybSendNorth, yeSendNorth, ybSendSouth, yeSendSouth, &
       ybRecvEast, yeRecvEast, ybRecvWest, yeRecvWest, &
       ybRecvNorth, yeRecvNorth, ybRecvSouth, yeRecvSouth, &
       x0, y0, NodeDimsAdvMnt%mzp, ghostZoneWidth, idim_type_3D, &
       AdvMntScaSendX, AdvMntScaRecvX, AdvMntScaSendY, AdvMntScaRecvY, &
       fldName_X=vc3d_inName, &
       fldName_Y=vc3d_outName)

    if (dumpLocal) then
       call MsgDump(h//" finishes with AdvMntUVSendX MessageSet:")
       call DumpMessageSet(AdvMntUVSendX)
       call MsgDump(h//" finishes with AdvMntUVRecvX MessageSet:")
       call DumpMessageSet(AdvMntUVRecvX)
       call MsgDump(h//" finishes with AdvMntDxDySendX MessageSet:")
       call DumpMessageSet(AdvMntDxDySendX)
       call MsgDump(h//" finishes with AdvMntDxDyRecvX MessageSet:")
       call DumpMessageSet(AdvMntDxDyRecvX)
       call MsgDump(h//" finishes with AdvMntDd0SendX MessageSet:")
       call DumpMessageSet(AdvMntDd0SendX)
       call MsgDump(h//" finishes with AdvMntDd0RecvX MessageSet:")
       call DumpMessageSet(AdvMntDd0RecvX)
       call MsgDump(h//" finishes with AdvMntDenSendX MessageSet:")
       call DumpMessageSet(AdvMntDenSendX)
       call MsgDump(h//" finishes with AdvMntDenRecvX MessageSet:")
       call DumpMessageSet(AdvMntDenRecvX)
       call MsgDump(h//" finishes with AdvMntScaSendX MessageSet:")
       call DumpMessageSet(AdvMntScaSendX)
       call MsgDump(h//" finishes with AdvMntScaRecvX MessageSet:")
       call DumpMessageSet(AdvMntScaRecvX)
       call MsgDump(h//" finishes with AdvMntUVSendY MessageSet:")
       call DumpMessageSet(AdvMntUVSendY)
       call MsgDump(h//" finishes with AdvMntUVRecvY MessageSet:")
       call DumpMessageSet(AdvMntUVRecvY)
       call MsgDump(h//" finishes with AdvMntDxDySendY MessageSet:")
       call DumpMessageSet(AdvMntDxDySendY)
       call MsgDump(h//" finishes with AdvMntDxDyRecvY MessageSet:")
       call DumpMessageSet(AdvMntDxDyRecvY)
       call MsgDump(h//" finishes with AdvMntDd0SendY MessageSet:")
       call DumpMessageSet(AdvMntDd0SendY)
       call MsgDump(h//" finishes with AdvMntDd0RecvY MessageSet:")
       call DumpMessageSet(AdvMntDd0RecvY)
       call MsgDump(h//" finishes with AdvMntDenSendY MessageSet:")
       call DumpMessageSet(AdvMntDenSendY)
       call MsgDump(h//" finishes with AdvMntDenRecvY MessageSet:")
       call DumpMessageSet(AdvMntDenRecvY)
       call MsgDump(h//" finishes with AdvMntScaSendY MessageSet:")
       call DumpMessageSet(AdvMntScaSendY)
       call MsgDump(h//" finishes with AdvMntScaRecvY MessageSet:")
       call DumpMessageSet(AdvMntScaRecvY)
    end if
  end subroutine CreateAdvMntMessageSet



    subroutine DestroyAdvMntMessageSet(&
       AdvMntUVSendX, AdvMntUVRecvX, &
       AdvMntUVSendY, AdvMntUVRecvY, &
       AdvMntDxDySendX, AdvMntDxDyRecvX, &
       AdvMntDxDySendY, AdvMntDxDyRecvY, &
       AdvMntDd0SendX, AdvMntDd0RecvX, &
       AdvMntDd0SendY, AdvMntDd0RecvY, &
       AdvMntDenSendX, AdvMntDenRecvX, &
       AdvMntDenSendY, AdvMntDenRecvY, &
       AdvMntScaSendX, AdvMntScaRecvX, &
       AdvMntScaSendY, AdvMntScaRecvY)

    type(MessageSet), pointer, intent(inout) :: AdvMntUVSendX
    type(MessageSet), pointer, intent(inout) :: AdvMntUVRecvX
    type(MessageSet), pointer, intent(inout) :: AdvMntUVSendY
    type(MessageSet), pointer, intent(inout) :: AdvMntUVRecvY
    type(MessageSet), pointer, intent(inout) :: AdvMntDxDySendX
    type(MessageSet), pointer, intent(inout) :: AdvMntDxDyRecvX
    type(MessageSet), pointer, intent(inout) :: AdvMntDxDySendY
    type(MessageSet), pointer, intent(inout) :: AdvMntDxDyRecvY
    type(MessageSet), pointer, intent(inout) :: AdvMntDd0SendX
    type(MessageSet), pointer, intent(inout) :: AdvMntDd0RecvX
    type(MessageSet), pointer, intent(inout) :: AdvMntDd0SendY
    type(MessageSet), pointer, intent(inout) :: AdvMntDd0RecvY
    type(MessageSet), pointer, intent(inout) :: AdvMntDenSendX
    type(MessageSet), pointer, intent(inout) :: AdvMntDenRecvX
    type(MessageSet), pointer, intent(inout) :: AdvMntDenSendY
    type(MessageSet), pointer, intent(inout) :: AdvMntDenRecvY
    type(MessageSet), pointer, intent(inout) :: AdvMntScaSendX
    type(MessageSet), pointer, intent(inout) :: AdvMntScaRecvX
    type(MessageSet), pointer, intent(inout) :: AdvMntScaSendY
    type(MessageSet), pointer, intent(inout) :: AdvMntScaRecvY

    character(len=*), parameter :: h="**(DestroyAdvMntMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" will destroy all AdvMntSend/RecvX/Y")
    end if
    call DestroyMessageSet(AdvMntUVSendX)
    call DestroyMessageSet(AdvMntUVRecvX)
    call DestroyMessageSet(AdvMntUVSendY)
    call DestroyMessageSet(AdvMntUVRecvY)
    call DestroyMessageSet(AdvMntDxDySendX)
    call DestroyMessageSet(AdvMntDxDyRecvX)
    call DestroyMessageSet(AdvMntDxDySendY)
    call DestroyMessageSet(AdvMntDxDyRecvY)
    call DestroyMessageSet(AdvMntDd0SendX)
    call DestroyMessageSet(AdvMntDd0RecvX)
    call DestroyMessageSet(AdvMntDd0SendY)
    call DestroyMessageSet(AdvMntDd0RecvY)
    call DestroyMessageSet(AdvMntDenSendX)
    call DestroyMessageSet(AdvMntDenRecvX)
    call DestroyMessageSet(AdvMntDenSendY)
    call DestroyMessageSet(AdvMntDenRecvY)
    call DestroyMessageSet(AdvMntScaSendX)
    call DestroyMessageSet(AdvMntScaRecvX)
    call DestroyMessageSet(AdvMntScaSendY)
    call DestroyMessageSet(AdvMntScaRecvY)
  end subroutine DestroyAdvMntMessageSet




  subroutine UpdateFieldAdressAtAcoustNew(&
       AcoustNewSend, &
       AcoustNewRecv, &
       fieldName, &
       field)
       
    type(MessageSet), pointer, intent(in) :: AcoustNewSend
    type(MessageSet), pointer, intent(in) :: AcoustNewRecv
    character(len=*), intent(in) :: fieldName
    real, pointer, intent(in) :: field(:,:,:)

    integer :: iMsg
    type(FieldSectionNode), pointer :: fsnode
    character(len=*), parameter :: h="**(UpdateFieldAdressAtAcoustNew)**"
    logical, parameter :: dumpLocal=.false.
    
    if (.not. associated(AcoustNewSend)) then
       call fatal_error(h//" AcoustNewSend not associated")
    else if (.not. associated(AcoustNewRecv)) then
       call fatal_error(h//" AcoustNewRecv not associated")
    end if

    do iMsg = 1, AcoustNewSend%nMsgs
       fsnode => AcoustNewSend%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, field, fieldName)
    end do
    
    do iMsg = 1, AcoustNewRecv%nMsgs
       fsnode => AcoustNewRecv%msgData(iMsg)%list%head
    call UpdateFieldAdress(fsnode%entry, field, fieldName)
    end do
  end subroutine UpdateFieldAdressAtAcoustNew
  
  
  subroutine UpdateFieldAdressAtAdvMnt(&
       AdvMntUVSendX, AdvMntUVRecvX, &
       AdvMntUVSendY, AdvMntUVRecvY, &
       AdvMntDxDySendX, AdvMntDxDyRecvX, &
       AdvMntDxDySendY, AdvMntDxDyRecvY, &
       AdvMntDd0SendX, AdvMntDd0RecvX, &
       AdvMntDd0SendY, AdvMntDd0RecvY, &
       AdvMntDenSendX, AdvMntDenRecvX, &
       AdvMntDenSendY, AdvMntDenRecvY, &
       AdvMntScaSendX, AdvMntScaRecvX, &
       AdvMntScaSendY, AdvMntScaRecvY, &
       u3d, v3d, dxtW, dytW, &
       dd0_3d, dd0_3du, dd0_3dv, dd0_3dw, &
       den0_3d, den1_3d, den2_3d, den3_3d, &
       vc3d_in, vc3d_out)
    type(MessageSet), pointer, intent(in) :: AdvMntUVSendX
    type(MessageSet), pointer, intent(in) :: AdvMntUVRecvX
    type(MessageSet), pointer, intent(in) :: AdvMntUVSendY
    type(MessageSet), pointer, intent(in) :: AdvMntUVRecvY
    type(MessageSet), pointer, intent(in) :: AdvMntDxDySendX
    type(MessageSet), pointer, intent(in) :: AdvMntDxDyRecvX
    type(MessageSet), pointer, intent(in) :: AdvMntDxDySendY
    type(MessageSet), pointer, intent(in) :: AdvMntDxDyRecvY
    type(MessageSet), pointer, intent(in) :: AdvMntDd0SendX
    type(MessageSet), pointer, intent(in) :: AdvMntDd0RecvX
    type(MessageSet), pointer, intent(in) :: AdvMntDd0SendY
    type(MessageSet), pointer, intent(in) :: AdvMntDd0RecvY
    type(MessageSet), pointer, intent(in) :: AdvMntDenSendX
    type(MessageSet), pointer, intent(in) :: AdvMntDenRecvX
    type(MessageSet), pointer, intent(in) :: AdvMntDenSendY
    type(MessageSet), pointer, intent(in) :: AdvMntDenRecvY
    type(MessageSet), pointer, intent(in) :: AdvMntScaSendX
    type(MessageSet), pointer, intent(in) :: AdvMntScaRecvX
    type(MessageSet), pointer, intent(in) :: AdvMntScaSendY
    type(MessageSet), pointer, intent(in) :: AdvMntScaRecvY
    real, pointer, intent(in) :: u3d(:,:,:)
    real, pointer, intent(in) :: v3d(:,:,:)
    real, pointer, intent(in) :: dxtW(:,:)
    real, pointer, intent(in) :: dytW(:,:)
    real, pointer, intent(in) :: dd0_3d(:,:,:)
    real, pointer, intent(in) :: dd0_3du(:,:,:)
    real, pointer, intent(in) :: dd0_3dv(:,:,:)
    real, pointer, intent(in) :: dd0_3dw(:,:,:)
    real, pointer, intent(in) :: den0_3d(:,:,:)
    real, pointer, intent(in) :: den1_3d(:,:,:)
    real, pointer, intent(in) :: den2_3d(:,:,:)
    real, pointer, intent(in) :: den3_3d(:,:,:)
    real, pointer, intent(in) :: vc3d_in(:,:,:)
    real, pointer, intent(in) :: vc3d_out(:,:,:)

    integer :: iMsg
    type(FieldSectionNode), pointer :: fsnode
    character(len=*), parameter :: h="**(UpdateFieldAdressAtAdvMnt)**"
    
    if (.not. associated(AdvMntUVSendX)) then
       call fatal_error(h//" AdvMntUVSendX not associated")
    else if (.not. associated(AdvMntUVRecvX)) then
       call fatal_error(h//" AdvMntUVRecvX not associated")
    else if (.not. associated(AdvMntUVSendY)) then
       call fatal_error(h//" AdvMntUVSendY not associated")
    else if (.not. associated(AdvMntUVRecvY)) then
       call fatal_error(h//" AdvMntUVRecvY not associated")
    else if (.not. associated(AdvMntDxDySendX)) then
       call fatal_error(h//" AdvMntDxDySendX not associated")
    else if (.not. associated(AdvMntDxDyRecvX)) then
       call fatal_error(h//" AdvMntDxDyRecvX not associated")
    else if (.not. associated(AdvMntDxDySendY)) then
       call fatal_error(h//" AdvMntDxDySendY not associated")
    else if (.not. associated(AdvMntDxDyRecvY)) then
       call fatal_error(h//" AdvMntDxDyRecvY not associated")
    else if (.not. associated(AdvMntDd0SendX)) then
       call fatal_error(h//" AdvMntDd0SendX not associated")
    else if (.not. associated(AdvMntDd0RecvX)) then
       call fatal_error(h//" AdvMntDd0RecvX not associated")
    else if (.not. associated(AdvMntDd0SendY)) then
       call fatal_error(h//" AdvMntDd0SendY not associated")
    else if (.not. associated(AdvMntDd0RecvY)) then
       call fatal_error(h//" AdvMntDd0RecvY not associated")
    else if (.not. associated(AdvMntDenSendX)) then
       call fatal_error(h//" AdvMntDenSendX not associated")
    else if (.not. associated(AdvMntDenRecvX)) then
       call fatal_error(h//" AdvMntDenRecvX not associated")
    else if (.not. associated(AdvMntDenSendY)) then
       call fatal_error(h//" AdvMntDenSendY not associated")
    else if (.not. associated(AdvMntDenRecvY)) then
       call fatal_error(h//" AdvMntDenRecvY not associated")
    else if (.not. associated(AdvMntScaSendX)) then
       call fatal_error(h//" AdvMntScaSendX not associated")
    else if (.not. associated(AdvMntScaRecvX)) then
       call fatal_error(h//" AdvMntScaRecvX not associated")
    else if (.not. associated(AdvMntScaSendY)) then
       call fatal_error(h//" AdvMntScaSendY not associated")
    else if (.not. associated(AdvMntScaRecvY)) then
       call fatal_error(h//" AdvMntScaRecvY not associated")
    end if
    
    do iMsg = 1, AdvMntUVSendX%nMsgs
       fsnode => AdvMntUVSendX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, u3d, "U3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, v3d, "V3D")
    end do
    do iMsg = 1, AdvMntUVRecvX%nMsgs
       fsnode => AdvMntUVRecvX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, u3d, "U3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, v3d, "V3D")
    end do
    do iMsg = 1, AdvMntUVSendY%nMsgs
       fsnode => AdvMntUVSendY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, u3d, "U3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, v3d, "V3D")
    end do
    do iMsg = 1, AdvMntUVRecvY%nMsgs
       fsnode => AdvMntUVRecvY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, u3d, "U3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, v3d, "V3D")
    end do


    
    do iMsg = 1, AdvMntDxDySendX%nMsgs
       fsnode => AdvMntDxDySendX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, dxtW, "DXTW")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dytW, "DYTW")
    end do
    do iMsg = 1, AdvMntDxDyRecvX%nMsgs
       fsnode => AdvMntDxDyRecvX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, dxtW, "DXTW")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dytW, "DYTW")
    end do
    do iMsg = 1, AdvMntDxDySendY%nMsgs
       fsnode => AdvMntDxDySendY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, dxtW, "DXTW")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dytW, "DYTW")
    end do
    do iMsg = 1, AdvMntDxDyRecvY%nMsgs
       fsnode => AdvMntDxDyRecvY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, dxtW, "DXTW")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dytW, "DYTW")
    end do

    
    do iMsg = 1, AdvMntDd0SendX%nMsgs
       fsnode => AdvMntDd0SendX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, dd0_3d, "DD0_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3du, "DD0_3DU")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3dv, "DD0_3DV")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3dw, "DD0_3DW")
    end do
    do iMsg = 1, AdvMntDd0RecvX%nMsgs
       fsnode => AdvMntDd0RecvX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, dd0_3d, "DD0_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3du, "DD0_3DU")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3dv, "DD0_3DV")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3dw, "DD0_3DW")
    end do
    do iMsg = 1, AdvMntDd0SendY%nMsgs
       fsnode => AdvMntDd0SendY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, dd0_3d, "DD0_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3du, "DD0_3DU")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3dv, "DD0_3DV")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3dw, "DD0_3DW")
    end do
    do iMsg = 1, AdvMntDd0RecvY%nMsgs
       fsnode => AdvMntDd0RecvY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, dd0_3d, "DD0_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3du, "DD0_3DU")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3dv, "DD0_3DV")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, dd0_3dw, "DD0_3DW")
    end do

    do iMsg = 1, AdvMntDenSendX%nMsgs
       fsnode => AdvMntDenSendX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, den0_3d, "DEN0_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den1_3d, "DEN1_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den2_3d, "DEN2_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den3_3d, "DEN3_3D")
    end do
    do iMsg = 1, AdvMntDenRecvX%nMsgs
       fsnode => AdvMntDenRecvX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, den0_3d, "DEN0_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den1_3d, "DEN1_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den2_3d, "DEN2_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den3_3d, "DEN3_3D")
    end do
    do iMsg = 1, AdvMntDenSendY%nMsgs
       fsnode => AdvMntDenSendY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, den0_3d, "DEN0_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den1_3d, "DEN1_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den2_3d, "DEN2_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den3_3d, "DEN3_3D")
    end do
    do iMsg = 1, AdvMntDenRecvY%nMsgs
       fsnode => AdvMntDenRecvY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, den0_3d, "DEN0_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den1_3d, "DEN1_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den2_3d, "DEN2_3D")
       fsnode => fsnode%next
       call UpdateFieldAdress(fsnode%entry, den3_3d, "DEN3_3D")
    end do

    do iMsg = 1, AdvMntScaSendX%nMsgs
       fsnode => AdvMntScaSendX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, vc3d_in, "VC3D_IN")
    end do
    do iMsg = 1, AdvMntScaRecvX%nMsgs
       fsnode => AdvMntScaRecvX%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, vc3d_in, "VC3D_IN")
    end do
    do iMsg = 1, AdvMntScaSendY%nMsgs
       fsnode => AdvMntScaSendY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, vc3d_out, "VC3D_OUT")
    end do
    do iMsg = 1, AdvMntScaRecvY%nMsgs
       fsnode => AdvMntScaRecvY%msgData(iMsg)%list%head
       call UpdateFieldAdress(fsnode%entry, vc3d_out, "VC3D_OUT")
    end do
  end subroutine UpdateFieldAdressAtAdvMnt
  



  
  
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





  subroutine PostSendRecvMsgsFixedAdress1D(SendMsg, RecvMsg, nzp, nxp, nyp)

    ! posts all nonblocking send and recv operations of
    ! a message set pair of variables

    type(MessageSet), pointer, intent(in) :: SendMsg
    type(MessageSet), pointer, intent(in) :: RecvMsg
    integer, intent(in) :: nzp
    integer, intent(in) :: nxp
    integer, intent(in) :: nyp

    integer :: iSend
    integer :: iRecv
    integer :: firstBuffer
    integer :: lastBuffer
    integer :: ierr
    type(MessageData), pointer :: msgData => null()
    type(FieldSection), pointer :: node => null()
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(PostSendRecvMsgsFixedAdress1D)**"
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
          call FillMessageDataBuffer(SendMsg%msgData(iSend), nzp, nxp, nyp)

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
  end subroutine PostSendRecvMsgsFixedAdress1D  




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





  subroutine WaitSendRecvMsgsFixedAdress1D(SendMsg, RecvMsg, nzp, nxp, nyp)
    type(MessageSet), pointer, intent(in) :: SendMsg
    type(MessageSet), pointer, intent(in) :: RecvMsg
    integer, intent(in) :: nzp
    integer, intent(in) :: nxp
    integer, intent(in) :: nyp

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
    character(len=*), parameter :: h="**(WaitSendRecvMsgsFixedAdress1D)**"
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

          call DecomposeMessageDataBuffer(RecvMsg%msgData(recvNbr), nzp, nxp, nyp)
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
  end subroutine WaitSendRecvMsgsFixedAdress1D  



  subroutine OneAcoustNewSendRecv(OneD, Neigh, nNeigh, &
       willSend, willRecv, NameSend, NameRecv, &
       sendDirection, recvDirection, Tag, &
       zbSend, zeSend, zbRecv, zeRecv, &
       xbSend, xeSend, xbRecv, xeRecv, &
       ybSend, yeSend, ybRecv, yeRecv, &
       x0, y0, mzp, fldName, &
       AcoustNewSend, AcoustNewRecv, field)

    logical, intent(in) :: OneD
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    integer, intent(in) :: nNeigh
    logical, intent(in) :: willSend(:)
    logical, intent(in) :: willRecv(:)
    character(len=*), intent(in) :: NameSend
    character(len=*), intent(in) :: NameRecv
    character(len=*), intent(in) :: sendDirection
    character(len=*), intent(in) :: recvDirection
    integer, intent(in) :: Tag
    integer, intent(in) :: zbSend
    integer, intent(in) :: zeSend
    integer, intent(in) :: zbRecv
    integer, intent(in) :: zeRecv
    integer, intent(in) :: xbSend(:)
    integer, intent(in) :: xeSend(:)
    integer, intent(in) :: xbRecv(:)
    integer, intent(in) :: xeRecv(:)
    integer, intent(in) :: ybSend(:)
    integer, intent(in) :: yeSend(:)
    integer, intent(in) :: ybRecv(:)
    integer, intent(in) :: yeRecv(:)
    integer, intent(in) :: x0
    integer, intent(in) :: y0
    integer, intent(in) :: mzp
    character(len=*), intent(in) :: fldName
    type(MessageSet), pointer, intent(out) :: AcoustNewSend
    type(MessageSet), pointer, intent(out) :: AcoustNewRecv
    real, pointer, optional, intent(in) :: field(:)

    integer :: nMsgs
    integer :: cntMsg
    integer :: iNeigh
    integer :: idim_type
    integer :: fieldSectionSize
    type(FieldSection), pointer :: oneFieldSection

    character(len=8) :: str(10)
    character(len=128) :: strLong
    character(len=*), parameter :: h="**(OneAcoustNewSendRecv)**"
    logical, parameter :: dumpLocal=.false.

    if (OneD) then
       idim_type=1
    else
       idim_type=3
    end if
    
    if (dumpLocal) then
       write(str(1),"(i8)") idim_type
       call MsgDump(h//" building "//&
            trim(NameSend)//", "//&
            trim(NameRecv)//" with idim_type="//&
            trim(adjustl(str(1))))
    end if

    ! create message set sends

    AcoustNewSend => CreateMessageSet(&
         NameSend, &
         sendDirection, &
         Tag, &
         willSend, &
         Neigh)

    ! insert field sections at each direction
    ! with null field addresses, to be updated whenever
    ! real addresses are known

    if (associated(AcoustNewSend)) then
       nMsgs = AcoustNewSend%nMsgs

       ! create list of Field Sections to send, one for
       ! each process to communicate and insert at the send MessageSet
       ! field section list

       ! since there is at most one neighbour node at each direction,
       ! there will be at most one MessageSet at each direction

       cntMsg = 0
       do iNeigh = 1, nNeigh
          if (willSend(iNeigh)) then

             ! insert send communications to east

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AcoustNewSend%name)))
             end if

             if (present(field)) then
                oneFieldSection =>  CreateFieldSection(field, fldName, idim_type, &
                     zbSend, zeSend, &
                     xbSend(iNeigh)-x0, xeSend(iNeigh)-x0, &
                     ybSend(iNeigh)-y0, yeSend(iNeigh)-y0)
             else
                fieldSectionSize=&
                     (zeSend-zbSend+1)* &
                     (xeSend(iNeigh)-xbSend(iNeigh)+1) * &
                     (yeSend(iNeigh)-ybSend(iNeigh)+1)
                oneFieldSection =>  CreateFieldSection(fldName, idim_type, &
                     zbSend, zeSend, &
                     xbSend(iNeigh)-x0, xeSend(iNeigh)-x0, &
                     ybSend(iNeigh)-y0, yeSend(iNeigh)-y0, &
                     fieldSectionSize)
             end if
             call AppendFieldSectionToMessageData(oneFieldSection, AcoustNewSend%msgData(cntMsg))

          end if
       end do
    end if

    ! create message set for recvs

    AcoustNewRecv => CreateMessageSet(&
         NameRecv, &
         recvDirection, &
         Tag, &
         willRecv, &
         Neigh)

    if ( associated(AcoustNewRecv)) then
       nMsgs = AcoustNewRecv%nMsgs
       cntMsg = 0
       do iNeigh = 1, nNeigh
          if (willRecv(iNeigh)) then

             ! insert recv communications from east

             cntMsg = cntMsg + 1
             if (cntMsg > nMsgs) then
                write(str(1),"(i8)") nMsgs
                call fatal_error(h//" nMsgs ("//&
                     trim(adjustl(str(1)))//") exceeded while inserting fields "//&
                     " at message "//trim(adjustl(AcoustNewRecv%name)))
             end if

             if (present(field)) then
                oneFieldSection =>  CreateFieldSection(field, fldName, idim_type, &
                     zbRecv, zeRecv, &
                     xbRecv(iNeigh)-x0, xeRecv(iNeigh)-x0, &
                     ybRecv(iNeigh)-y0, yeRecv(iNeigh)-y0)
             else
                fieldSectionSize=&
                     (zeRecv-zbRecv+1)* &
                     (xeRecv(iNeigh)-xbRecv(iNeigh)+1) * &
                     (yeRecv(iNeigh)-ybRecv(iNeigh)+1)
                oneFieldSection =>  CreateFieldSection(fldName, idim_type, &
                     zbRecv, zeRecv, &
                     xbRecv(iNeigh)-x0, xeRecv(iNeigh)-x0, &
                     ybRecv(iNeigh)-y0, yeRecv(iNeigh)-y0, &
                     fieldSectionSize)
             end if
             call AppendFieldSectionToMessageData(oneFieldSection, AcoustNewRecv%msgData(cntMsg))
          end if
       end do
    end if

    if (dumpLocal) then
       call MsgDump(h//" built "//trim(NameSend)//" MessageSet:")
       call DumpMessageSet(AcoustNewSend)
       call MsgDump(h//" built "//trim(NameRecv)//" MessageSet:")
       call DumpMessageSet(AcoustNewRecv)
    end if
  end subroutine OneAcoustNewSendRecv


  subroutine CreateAcoustNewMessageSet(&
     GridSize, ParEnv, Neigh, GlobalOwn, GlobalWithGhost, &
     TagDiv, AcoustNewDivSend, AcoustNewDivRecv, &
     TagPP, AcoustNewPPSend, AcoustNewPPRecv, &
     TagAlpha, AcoustNewAlphaSend, AcoustNewAlphaRecv, &
     TagTht, AcoustNewThtSend, AcoustNewThtRecv, tht)

    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    integer, intent(in) :: TagDiv
    type(MessageSet), pointer, intent(inout) :: AcoustNewDivSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewDivRecv
    integer, intent(in) :: TagPP
    type(MessageSet), pointer, intent(inout) :: AcoustNewPPSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewPPRecv
    integer, intent(in) :: TagAlpha
    type(MessageSet), pointer, intent(inout) :: AcoustNewAlphaSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewAlphaRecv
    integer, intent(in) :: TagTht
    type(MessageSet), pointer, intent(inout) :: AcoustNewThtSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewThtRecv
    real, pointer, intent(in) :: tht(:)
    
    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    integer :: mzp
    integer :: x0
    integer :: y0
    
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

    character(len=*), parameter :: sendDirection="send"
    character(len=*), parameter :: recvDirection="recv"
    
    character(len=*), parameter :: NameAcoustNew="AcoustNew"
    character(len=*), parameter :: NameSendDiv="AcoustNewDivSend"
    character(len=*), parameter :: NameRecvDiv="AcoustNewDivRecv"
    character(len=*), parameter :: NameSendPP="AcoustNewPPSend"
    character(len=*), parameter :: NameRecvPP="AcoustNewPPRecv"
    character(len=*), parameter :: NameSendAlpha="AcoustNewAlphaSend"
    character(len=*), parameter :: NameRecvAlpha="AcoustNewAlphaRecv"
    character(len=*), parameter :: NameSendTht="AcoustNewThtSend"
    character(len=*), parameter :: NameRecvTht="AcoustNewThtRecv"
    
    character(len=*), parameter :: h="**(CreateAcoustNewMessageSet)**"
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
       AcoustNewDivSend => null()
       AcoustNewDivRecv => null()
       AcoustNewPPSend => null()
       AcoustNewPPRecv => null()
       AcoustNewAlphaSend => null()
       AcoustNewAlphaRecv => null()
       AcoustNewThtSend => null()
       AcoustNewThtRecv => null()
       return
    end if

    if (dumpLocal) then
       call MsgDump(h//" will create AcoustNewOneSend/Recv")
    end if

    myNum  = ParEnv%myNum
    nMachs = ParEnv%nMachs
    nNeigh = Neigh%nNeigh
    mzp = GridSize%nnzp
    x0 = GlobalWithGhost%xb(myNum) - 1
    y0 = GlobalWithGhost%yb(myNum) - 1
    
    call NodesToSendRecvMessages( &
         thisNode=myNum, &
         Neigh=Neigh, &
         GlobalOwn=GlobalOwn, &
         xbToUpdate=GlobalWithGhost%xb, &
         xeToUpdate=GlobalWithGhost%xe, &
         ybToUpdate=GlobalWithGhost%yb, &
         yeToUpdate=GlobalWithGhost%ye, &
         xbSend=xbSend, &
         xeSend=xeSend, &
         ybSend=ybSend, &
         yeSend=yeSend, &
         willSend=willSend, &
         xbRecv=xbRecv, &
         xeRecv=xeRecv, &
         ybRecv=ybRecv, &
         yeRecv=yeRecv, &
         willRecv=willRecv, &
         varName=NameAcoustNew)

    ! create message sets

    call OneAcoustNewSendRecv(.false., Neigh, nNeigh, &
         willSend, willRecv, NameSendDiv, NameRecvDiv, &
         sendDirection, recvDirection, TagDiv, &
         1, mzp, 1, mzp, &
         xbSend, xeSend, xbRecv, xeRecv, &
         ybSend, yeSend, ybRecv, yeRecv, &
         x0, y0, mzp, "DIV", &
         AcoustNewDivSend, AcoustNewDivRecv)

    call OneAcoustNewSendRecv(.false., Neigh, nNeigh, &
         willSend, willRecv, NameSendPP, NameRecvPP, &
         sendDirection, recvDirection, TagPP, &
         1, mzp, 1, mzp, &
         xbSend, xeSend, xbRecv, xeRecv, &
         ybSend, yeSend, ybRecv, yeRecv, &
         x0, y0, mzp, "PP", &
         AcoustNewPPSend, AcoustNewPPRecv)

    call OneAcoustNewSendRecv(.false., Neigh, nNeigh, &
         willSend, willRecv, NameSendAlpha, NameRecvAlpha, &
         sendDirection, recvDirection, TagAlpha, &
         1, mzp, 1, mzp, &
         xbSend, xeSend, xbRecv, xeRecv, &
         ybSend, yeSend, ybRecv, yeRecv, &
         x0, y0, mzp, "ALPHA", &
         AcoustNewAlphaSend, AcoustNewAlphaRecv)

    call OneAcoustNewSendRecv(.true., Neigh, nNeigh, &
         willSend, willRecv, NameSendTht, NameRecvTht, &
         sendDirection, recvDirection, TagTht, &
         1, mzp, 1, mzp, &
         xbSend, xeSend, xbRecv, xeRecv, &
         ybSend, yeSend, ybRecv, yeRecv, &
         x0, y0, mzp, "THT", &
         AcoustNewThtSend, AcoustNewThtRecv, tht)
    
    if (dumpLocal) then
       call MsgDump(h//" finishes with AcoustNewDivSend MessageSet:")
       call DumpMessageSet(AcoustNewDivSend)
       call MsgDump(h//" finishes with AcoustNewDivRecv MessageSet:")
       call DumpMessageSet(AcoustNewDivRecv)
       call MsgDump(h//" finishes with AcoustNewPPSend MessageSet:")
       call DumpMessageSet(AcoustNewPPSend)
       call MsgDump(h//" finishes with AcoustNewPPRecv MessageSet:")
       call DumpMessageSet(AcoustNewPPRecv)
       call MsgDump(h//" finishes with AcoustNewAlphaSend MessageSet:")
       call DumpMessageSet(AcoustNewAlphaSend)
       call MsgDump(h//" finishes with AcoustNewAlphaRecv MessageSet:")
       call DumpMessageSet(AcoustNewAlphaRecv)
       call MsgDump(h//" finishes with AcoustNewThtSend MessageSet:")
       call DumpMessageSet(AcoustNewThtSend)
       call MsgDump(h//" finishes with AcoustNewThtRecv MessageSet:")
       call DumpMessageSet(AcoustNewThtRecv)
    end if
  end subroutine CreateAcoustNewMessageSet
  




  subroutine DestroyAcoustNewMessageSet( &
       AcoustNewDivSend, AcoustNewDivRecv, &
       AcoustNewPPSend, AcoustNewPPRecv, &
       AcoustNewAlphaSend, AcoustNewAlphaRecv, &
       AcoustNewThtSend, AcoustNewThtRecv)

    type(MessageSet), pointer, intent(inout) :: AcoustNewDivSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewDivRecv
    type(MessageSet), pointer, intent(inout) :: AcoustNewPPSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewPPRecv
    type(MessageSet), pointer, intent(inout) :: AcoustNewAlphaSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewAlphaRecv
    type(MessageSet), pointer, intent(inout) :: AcoustNewThtSend
    type(MessageSet), pointer, intent(inout) :: AcoustNewThtRecv
    character(len=*), parameter :: h="**(DestroyAcoustNewMessageSet)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" will destroy AcoustNew...Send/Recv")
    end if

    call DestroyMessageSet(AcoustNewDivSend)
    call DestroyMessageSet(AcoustNewDivRecv)
    call DestroyMessageSet(AcoustNewPPSend)
    call DestroyMessageSet(AcoustNewPPRecv)
    call DestroyMessageSet(AcoustNewAlphaSend)
    call DestroyMessageSet(AcoustNewAlphaRecv)
    call DestroyMessageSet(AcoustNewThtSend)
    call DestroyMessageSet(AcoustNewThtRecv)

  end subroutine DestroyAcoustNewMessageSet
  
       
end module ModMessageSet
