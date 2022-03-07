module ModMessageData

  ! one message of a message passing operation is composed
  ! by a message data and a message envelop.
  !
  ! This module implements the message data of a message
  ! passing operation. Consequently, it contains all data
  ! to communicate to a single rank on a single operation.
  !
  ! Data to communicate is a list of field sections. A single
  ! message passing operation collapses all field sections
  ! into a single communication buffer.
  !
  ! This module provides operations to insert/remove field sections at
  ! the field section list, to copy in and copy out all entries 
  ! of the field section list to/from the communication buffer,
  ! procedures to create and destroy the communication buffer
  ! and procedures to send and receive the buffer.
  !
  ! Each field section list is constant throughout the computation.
  ! The list should be build at initialization and used whenever
  ! required. This organization avoids recomputing indices of field
  ! sections to copy in / copy out.

  use ModParallelEnvironment, only: &
       MsgDump

  use ModFieldSection, only: &
       FieldSection, &
       FieldSectionSize, &
       DestroyFieldSection, &
       StringFieldSection, &
       DumpFieldSection, &
       FieldSectionData2Buffer, &
       FieldSectionData2BufferVariableAdressArr, &
       FieldSectionData2BufferVariableAdressScalar, &
       Buffer2FieldSectionData


  use ParLib, only: &
       parf_get_noblock_real, &
       parf_send_noblock_real


  use ModFieldSectionList, only: &
       FieldSectionNode, &
       FieldSectionList, &
       CreateFieldSectionNode, &
       CreateFieldSectionList, &
       DumpFieldSectionList, &
       AppendNodeToFieldSectionList, &
       FieldSectionListHeadNode, &
       FieldSectionAtNode, &
       NextFieldSectionNodeAtList, &
       FindFieldNamed

  implicit none
  include "mpif.h"

  private
  public :: MessageData
  public :: CreateMessageData
  public :: DestroyMessageData
  public :: DumpMessageData
  public :: AppendFieldSectionToMessageData
  public :: AllocateMessageDataBuffer
  public :: FillMessageDataBuffer
  public :: PostSendMessageData
  public :: PostRecvMessageData
  public :: Buffer2FieldSectionData
  public :: DecomposeMessageDataBuffer
  public :: DeallocateMessageDataBuffer
  public :: FindFieldNamed
  public :: FillMessageDataBufferVariableAdressArr
  public :: FillMessageDataBufferVariableAdressScalar
  public :: FillMessageDataBufferVariableAdressOneArr

  type MessageData
     private
     ! data to communicate to one node in one message

     ! message communicates all field values stored
     ! at the field section list

     ! message send/receive the buffer "buf"

     ! prior to a send communication, procedure
     ! FillMessageDataBuffer
     ! copies field values stored at all field sections
     ! on the field section list to the buffer

     ! after a receive communication, procedure
     ! Buffer2FieldSectionData
     ! copies field values stored at the buffer
     ! to all field sections of the field section list

     ! a field section is appended to the field section list
     ! by procedure AppendFieldSectionToMessageData

     integer :: bufSize=0
     ! message buffer size
     real, allocatable :: buf(:)
     ! message buffer
     character(len=64) :: name=""
     ! message data name
     character(len=4) :: direction=""
     ! message data direction, used only for
     ! documentation; one of "send" or "recv"
     type (FieldSectionList), pointer :: list => null()
     ! list of Field Sections to communicate
  end type MessageData

  interface FindFieldNamed
     module procedure FindFieldSectionAtMessageData
  end interface FindFieldNamed

  interface FillMessageDataBuffer
     module procedure FillMessageDataBufferFixedAdress
     module procedure FillMessageDataBufferFixedAdress1D
  end interface FillMessageDataBuffer

  interface DecomposeMessageDataBuffer
     module procedure DecomposeMessageDataBufferFixedAdress
     module procedure DecomposeMessageDataBufferFixedAdress1D
     module procedure DecomposeMessageDataBufferVariableAdress
     module procedure DecomposeMessageDataBufferVariableAdressOneArr
  end interface DecomposeMessageDataBuffer
contains






  subroutine CreateMessageData (oneMessageData, name, direction)

    ! initialize components of a MessageData variable

    type(MessageData), intent(inout) :: oneMessageData
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: direction

    character(len=*), parameter :: h="**(CreateMessageData)**"
    logical, parameter :: dumpLocal=.false.

    oneMessageData%bufSize = 0
    oneMessageData%name = trim(adjustl(name))
    oneMessageData%direction = trim(adjustl(direction))
    oneMessageData%list => CreateFieldSectionList()
    if (dumpLocal) then
       call DumpMessageData(oneMessageData, h)
    end if
  end subroutine CreateMessageData






  subroutine DestroyMessageData(oneMessageData)

    ! deallocate components of a MessageData variable

    type(MessageData) :: oneMessageData

    integer :: ierr
    type(FieldSection), pointer:: this
    type(FieldSection), pointer:: next
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(DestroyMessageData)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" of "//trim(adjustl(oneMessageData%name)))
    end if
    if (allocated(oneMessageData%buf)) then
       call DeallocateMessageDataBuffer(oneMessageData)
    end if
    oneMessageData%bufSize = 0
    nullify(oneMessageData%list)
  end subroutine DestroyMessageData




  subroutine DumpMessageData(oneMessageData, strMsg)

    ! dumps components of a MessageData variable

    type(MessageData), intent(in) :: oneMessageData
    character(len=*), intent(in), optional :: strMsg

    type(FieldSection), pointer :: this
    character(len=128) :: msgHead
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DumpMessageData)**"

    if (present(strMsg)) then
       msgHead=strMsg
    else
       msgHead=h
    end if
    if (oneMessageData%name == "") then
       call MsgDump(trim(adjustl(msgHead))//" empty message data")
    else
       write(c0,"(i8)") oneMessageData%bufSize
       call MsgDump(trim(adjustl(msgHead))//" "//&
            trim(adjustl(oneMessageData%direction))//" message named "//&
            trim(adjustl(oneMessageData%name))//"; bufSize="//&
            trim(adjustl(c0))//"; field section list=")
       call DumpFieldSectionList(oneMessageData%list)
    end if
  end subroutine DumpMessageData




  subroutine AppendFieldSectionToMessageData (oneFieldSection, oneMessageData)

    ! appends a Field Section pointer to the tail of the Field Section list of
    ! a Message Data variable

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    type(MessageData), intent(inout) :: oneMessageData

    type(FieldSectionNode), pointer :: this
    character(len=8) :: c0
    character(len=*), parameter :: h="**(AppendFieldSectionToMessageData)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" oneFieldSection not associated")
    end if

    this => CreateFieldSectionNode(oneFieldSection)
    call AppendNodeToFieldSectionList(this, oneMessageData%list)
    oneMessageData%bufSize = oneMessageData%bufSize + FieldSectionSize(oneFieldSection)
    if (dumpLocal) then
       write(c0,"(i8)") oneMessageData%bufSize
       call MsgDump(h//" appended "//trim(adjustl(StringFieldSection(oneFieldSection)))//&
            " to Message Data "//trim(adjustl(oneMessageData%name))//"; current bufSize is "//&
            trim(adjustl(c0)))
    end if
  end subroutine AppendFieldSectionToMessageData






  subroutine AllocateMessageDataBuffer(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    ! allocates buf component of a Message Set variable with
    ! size of Message Set component bufSize

    integer :: bufSize
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(AllocateMessageDataBuffer)**"
    logical, parameter :: dumpLocal=.false.

    bufSize=oneMessageData%bufSize

    if (allocated(oneMessageData%buf)) then
       call fatal_error(h//" buf of "//trim(adjustl(oneMessageData%name))//&
            " already allocated")
    end if

    allocate(oneMessageData%buf(bufSize), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") bufSize
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocation of oneMessageData("//&
            trim(adjustl(c0))//") fails with stat="//&
            trim(adjustl(c1)))
    end if
    if (dumpLocal) then
       write(c0,"(i8)") oneMessageData%bufSize
       call MsgDump(h//" allocated buf of size "//trim(adjustl(c0))//&
            " to Message Data "//trim(adjustl(oneMessageData%name)))
    end if
  end subroutine AllocateMessageDataBuffer






  subroutine FillMessageDataBufferFixedAdress(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    ! copy field section values of the entire field section list of
    ! the Message Data variable to the buffer of the Message Data variable

    integer :: bufStart
    type(FieldSectionNode), pointer :: thisNode
    type(FieldSection), pointer :: this
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(FillMessageDataBufferFixedAdress)**"

    if (.not. allocated(oneMessageData%buf)) then
       call fatal_error(h//" not allocated buf for message data "//&
            trim(adjustl(oneMessageData%name)))
    end if

    if (dumpLocal) then
       call MsgDump(h//" to Message Data "//trim(adjustl(oneMessageData%name)))
    end if
    bufStart=1
    thisNode => FieldSectionListHeadNode(oneMessageData%list)
    do while (associated(thisNode))
       this => FieldSectionAtNode(thisNode)
       call FieldSectionData2Buffer(&
            this, &
            oneMessageData%buf, &
            bufStart, &
            oneMessageData%bufsize)
       thisNode => NextFieldSectionNodeAtList(thisNode)
    end do
  end subroutine FillMessageDataBufferFixedAdress





  subroutine FillMessageDataBufferFixedAdress1D(oneMessageData, nzp, nxp, nyp)
    type(MessageData), intent(inout) :: oneMessageData
    integer, intent(in) :: nzp
    integer, intent(in) :: nxp
    integer, intent(in) :: nyp

    ! copy field section values of the entire field section list of
    ! the Message Data variable to the buffer of the Message Data variable

    integer :: bufStart
    type(FieldSectionNode), pointer :: thisNode
    type(FieldSection), pointer :: this
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(FillMessageDataBufferFixedAdress1D)**"

    if (.not. allocated(oneMessageData%buf)) then
       call fatal_error(h//" not allocated buf for message data "//&
            trim(adjustl(oneMessageData%name)))
    end if

    if (dumpLocal) then
       call MsgDump(h//" to Message Data "//trim(adjustl(oneMessageData%name)))
    end if
    bufStart=1
    thisNode => FieldSectionListHeadNode(oneMessageData%list)
    do while (associated(thisNode))
       this => FieldSectionAtNode(thisNode)
       call FieldSectionData2Buffer(&
            this, &
            nzp, nxp, nyp, &
            oneMessageData%buf, &
            bufStart, &
            oneMessageData%bufsize)
       thisNode => NextFieldSectionNodeAtList(thisNode)
    end do
  end subroutine FillMessageDataBufferFixedAdress1D  





  subroutine FillMessageDataBufferVariableAdressArr(oneMessageData, scp, ufx, vfx, wfx)

    type(MessageData), intent(inout) :: oneMessageData
    real, pointer, intent(in):: scp(:,:,:)
    real, intent(in):: ufx(:,:,:)
    real, intent(in):: vfx(:,:,:)
    real, intent(in):: wfx(:,:,:)

    ! copy field section values of the entire field section list of
    ! the Message Data variable to the buffer of the Message Data variable

    integer :: bufStart
    type(FieldSectionNode), pointer :: oneNode
    type(FieldSection), pointer :: oneFieldSection
    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(FillMessageDataBufferVariableAdressArr)**"

    if (.not. associated(scp)) then
       call fatal_error(h//" invoked with null scp")
    end if

    if (.not. allocated(oneMessageData%buf)) then
       call fatal_error(h//" not allocated buf for message data "//&
            trim(adjustl(oneMessageData%name)))
    end if

    if (dumpLocal) then
       call MsgDump(h//" to Message Data "//trim(adjustl(oneMessageData%name)))
       write(str(1),"(i8)") lbound(scp,1)
       write(str(2),"(i8)") ubound(scp,1)
       write(str(3),"(i8)") lbound(scp,2)
       write(str(4),"(i8)") ubound(scp,2)
       write(str(5),"(i8)") lbound(scp,3)
       write(str(6),"(i8)") ubound(scp,3)
       call MsgDump(h//" scp dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") lbound(ufx,1)
       write(str(2),"(i8)") ubound(ufx,1)
       write(str(3),"(i8)") lbound(ufx,2)
       write(str(4),"(i8)") ubound(ufx,2)
       write(str(5),"(i8)") lbound(ufx,3)
       write(str(6),"(i8)") ubound(ufx,3)
       call MsgDump(h//" ufx dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") lbound(vfx,1)
       write(str(2),"(i8)") ubound(vfx,1)
       write(str(3),"(i8)") lbound(vfx,2)
       write(str(4),"(i8)") ubound(vfx,2)
       write(str(5),"(i8)") lbound(vfx,3)
       write(str(6),"(i8)") ubound(vfx,3)
       call MsgDump(h//" vfx dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") lbound(wfx,1)
       write(str(2),"(i8)") ubound(wfx,1)
       write(str(3),"(i8)") lbound(wfx,2)
       write(str(4),"(i8)") ubound(wfx,2)
       write(str(5),"(i8)") lbound(wfx,3)
       write(str(6),"(i8)") ubound(wfx,3)
       call MsgDump(h//" wfx dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
    end if

    bufStart=1
    oneNode => FieldSectionListHeadNode(oneMessageData%list)
    oneFieldSection => FieldSectionAtNode(oneNode)
    call FieldSectionData2BufferVariableAdressArr(scp, size(scp,1), oneFieldSection, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
    oneNode => NextFieldSectionNodeAtList(oneNode)
    oneFieldSection => FieldSectionAtNode(oneNode)
    call FieldSectionData2BufferVariableAdressArr(ufx, size(ufx,1), oneFieldSection, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
    oneNode => NextFieldSectionNodeAtList(oneNode)
    oneFieldSection => FieldSectionAtNode(oneNode)
    call FieldSectionData2BufferVariableAdressArr(vfx, size(vfx,1), oneFieldSection, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
    oneNode => NextFieldSectionNodeAtList(oneNode)
    oneFieldSection => FieldSectionAtNode(oneNode)
    call FieldSectionData2BufferVariableAdressArr(wfx, size(wfx,1), oneFieldSection, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
  end subroutine FillMessageDataBufferVariableAdressArr








  subroutine FillMessageDataBufferVariableAdressScalar(oneMessageData, scp, ufx, vfx, wfx)

    type(MessageData), intent(inout) :: oneMessageData
    real, intent(in):: scp
    real, intent(in):: ufx(:,:,:)
    real, intent(in):: vfx(:,:,:)
    real, intent(in):: wfx(:,:,:)

    ! copy field section values of the entire field section list of
    ! the Message Data variable to the buffer of the Message Data variable

    ! it is assumed that the "first dimension" of the field pointed by
    ! scalar scp is the same as the first dimension of ufx

    integer :: nxp
    integer :: nyp
    integer :: nzp
    integer :: bufStart
    type(FieldSectionNode), pointer :: oneNode
    type(FieldSection), pointer :: oneFieldSection
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(FillMessageDataBufferVariableAdressScalar)**"


    if (.not. allocated(oneMessageData%buf)) then
       call fatal_error(h//" not allocated buf for message data "//&
            trim(adjustl(oneMessageData%name)))
    end if

    if (dumpLocal) then
       call MsgDump(h//" to Message Data "//trim(adjustl(oneMessageData%name)))
    end if

    ! assume that field pointed by scp
    ! have the same shape as ufx
    
    nzp=size(ufx,1)
    nxp=size(ufx,2)
    nyp=size(ufx,3)
    
    bufStart=1
    oneNode => FieldSectionListHeadNode(oneMessageData%list)
    oneFieldSection => FieldSectionAtNode(oneNode)
    call FieldSectionData2BufferVariableAdressScalar(scp, &
         nzp, nxp, nyp, oneFieldSection, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
    oneNode => NextFieldSectionNodeAtList(oneNode)
    oneFieldSection => FieldSectionAtNode(oneNode)
    call FieldSectionData2BufferVariableAdressArr(ufx, nzp, oneFieldSection, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
    oneNode => NextFieldSectionNodeAtList(oneNode)
    oneFieldSection => FieldSectionAtNode(oneNode)
    call FieldSectionData2BufferVariableAdressArr(vfx, nzp, oneFieldSection, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
    oneNode => NextFieldSectionNodeAtList(oneNode)
    oneFieldSection => FieldSectionAtNode(oneNode)
    call FieldSectionData2BufferVariableAdressArr(wfx, nzp, oneFieldSection, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
  end subroutine FillMessageDataBufferVariableAdressScalar





  subroutine DecomposeMessageDataBufferFixedAdress(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    ! copy all field section values from the buffer of the Message Data variable
    ! to the field section pointed by the Message Data field section list 

    integer :: bufStart
    type(FieldSectionNode), pointer :: thisNode
    type(FieldSection), pointer :: this
    character(len=*), parameter :: h="**(DecomposeMessageDataBufferFixedAdress)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//"  of Message Data "//trim(adjustl(oneMessageData%name)))
    end if
    bufStart=1
    thisNode => FieldSectionListHeadNode(oneMessageData%list)
    do while (associated(thisNode))
       this => FieldSectionAtNode(thisNode)
       call Buffer2FieldSectionData(&
            this, &
            oneMessageData%buf, &
            bufStart, &
            oneMessageData%bufsize)
       thisNode => NextFieldSectionNodeAtList(thisNode)
    end do
  end subroutine DecomposeMessageDataBufferFixedAdress





  subroutine DecomposeMessageDataBufferFixedAdress1D(oneMessageData, nzp, nxp, nyp)
    type(MessageData), intent(inout) :: oneMessageData
    integer, intent(in) :: nzp
    integer, intent(in) :: nxp
    integer, intent(in) :: nyp

    ! copy all field section values from the buffer of the Message Data variable
    ! to the field section pointed by the Message Data field section list 

    integer :: bufStart
    type(FieldSectionNode), pointer :: thisNode
    type(FieldSection), pointer :: this
    character(len=*), parameter :: h="**(DecomposeMessageDataBufferFixedAdress1D)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//"  of Message Data "//trim(adjustl(oneMessageData%name)))
    end if
    bufStart=1
    thisNode => FieldSectionListHeadNode(oneMessageData%list)
    do while (associated(thisNode))
       this => FieldSectionAtNode(thisNode)
       call Buffer2FieldSectionData(&
            this, &
            nzp, nxp, nyp, &
            oneMessageData%buf, &
            bufStart, &
            oneMessageData%bufsize)
       thisNode => NextFieldSectionNodeAtList(thisNode)
    end do
  end subroutine DecomposeMessageDataBufferFixedAdress1D





  subroutine DecomposeMessageDataBufferVariableAdress(oneMessageData, &
       msgStartZ, msgEndZ, scr, ufx_local, wfx_local, vfx_local)
    type(MessageData), intent(inout) :: oneMessageData
    integer, intent(in) :: msgStartZ
    integer, intent(in) :: msgEndZ
    real, pointer, intent(in):: scr(:,:,:)
    real, pointer, intent(in):: ufx_local(:,:,:)
    real, pointer, intent(in):: vfx_local(:,:,:)
    real, pointer, intent(in):: wfx_local(:,:,:)

    ! copy all field section values from the buffer of the Message Data variable
    ! to the field section pointed by the Message Data field section list 

    integer :: bufStart
    type(FieldSectionNode), pointer :: oneNode
    type(FieldSection), pointer :: oneSection
    character(len=*), parameter :: h="**(DecomposeMessageDataBufferVariableAdress)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//"  of Message Data "//trim(adjustl(oneMessageData%name)))
    end if
    bufStart=1
    oneNode => FieldSectionListHeadNode(oneMessageData%list)
    oneSection => FieldSectionAtNode(oneNode)
    call Buffer2FieldSectionData(oneSection, &
         scr, msgStartZ, msgEndZ, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
    oneNode => NextFieldSectionNodeAtList(oneNode)
    oneSection => FieldSectionAtNode(oneNode)
    call Buffer2FieldSectionData(oneSection, &
         ufx_local, msgStartZ, msgEndZ, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
    oneNode => NextFieldSectionNodeAtList(oneNode)
    oneSection => FieldSectionAtNode(oneNode)
    call Buffer2FieldSectionData(oneSection, &
         vfx_local, msgStartZ, msgEndZ, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
    oneNode => NextFieldSectionNodeAtList(oneNode)
    oneSection => FieldSectionAtNode(oneNode)
    call Buffer2FieldSectionData(oneSection, &
         wfx_local, msgStartZ, msgEndZ, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
  end subroutine DecomposeMessageDataBufferVariableAdress






  subroutine DeallocateMessageDataBuffer(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    ! deallocates buf componenet of a Message Set variable

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DeallocateMessageDataBuffer)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" deallocate buf of "//&
            trim(adjustl(oneMessageData%name)))
    end if
    if (.not. allocated(oneMessageData%buf)) then
       call fatal_error(h//" "//&
            trim(adjustl(oneMessageData%name))//&
            " not allocated")
    else
       deallocate(oneMessageData%buf, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocation of oneMessageData "//&
               trim(adjustl(oneMessageData%name))//&
               " fails with stat="//trim(adjustl(c0)))
       end if
    end if
  end subroutine DeallocateMessageDataBuffer





  subroutine PostRecvMessageData(oneMessageData, &
       otherProc, tag, request, preMsgString)
    type(MessageData), intent(inout) :: oneMessageData
    integer, intent(in) :: otherProc
    integer, intent(in) :: tag
    integer, intent(inout) :: request
    character(len=*), optional :: preMsgString

    character(len=*), parameter :: h="**(PostRecvMessageData)**"
    character(len=8) :: c0, c1, c2, c3
    character(len=128) :: msgString
    logical, parameter :: dumpLocal=.false.

    if (.not. allocated(oneMessageData%buf)) then
       call fatal_error(h//" buf not allocated")
    end if

    call parf_get_noblock_real(&
         oneMessageData%buf, &
         oneMessageData%bufSize, &
         otherProc, &
         tag, &
         request)

    if (dumpLocal) then
       write(c0,"(i8)") oneMessageData%bufSize
       write(c1,"(i8)") otherProc
       write(c2,"(i8)") tag
       if (request == MPI_REQUEST_NULL) then
          c3="NULL"
       else
          write(c3,"(Z8)") request
       end if
       if (present(preMsgString)) then
          msgString=preMsgString
       else
          msgString=h
       end if
       call MsgDump(&
            trim(adjustl(msgString))//&
            " post recv from MPI rank "//trim(adjustl(c1))//&
            " with buffer size "//trim(adjustl(c0))//&
            " tag "//trim(adjustl(c2))//&
            " and request "//trim(adjustl(c3)))
    end if
  end subroutine PostRecvMessageData





  subroutine PostSendMessageData(oneMessageData, &
       otherProc, tag, request, preMsgString)
    type(MessageData), intent(in) :: oneMessageData
    integer, intent(in) :: otherProc
    integer, intent(in) :: tag
    integer, intent(inout) :: request
    character(len=*), optional :: preMsgString

    character(len=*), parameter :: h="**(PostSendMessageData)**"
    character(len=8) :: c0, c1, c2, c3
    character(len=128) :: msgString
    logical, parameter :: dumpLocal=.false.

    if (.not. allocated(oneMessageData%buf)) then
       call fatal_error(h//" buf not allocated")
    end if

    call parf_send_noblock_real(&
         oneMessageData%buf, &
         oneMessageData%bufSize, &
         otherProc, &
         tag, &
         request)

    if (dumpLocal) then
       write(c0,"(i8)") oneMessageData%bufSize
       write(c1,"(i8)") otherProc
       write(c2,"(i8)") tag
       if (request == MPI_REQUEST_NULL) then
          c3="NULL"
       else
          write(c3,"(Z8)") request
       end if
       if (present(preMsgString)) then
          msgString=preMsgString
       else
          msgString=h
       end if
       call MsgDump(&
            trim(adjustl(msgString))//&
            " post send to MPI rank "//trim(adjustl(c1))//&
            " with buffer size "//trim(adjustl(c0))//&
            " tag "//trim(adjustl(c2))//&
            " and request "//trim(adjustl(c3)))
    end if
  end subroutine PostSendMessageData





  function FindFieldSectionAtMessageData(oneMessageData, fieldName) result(node)
    type(MessageData), pointer, intent(in) :: oneMessageData
    character(len=*), intent(in) :: fieldName
    type(FieldSection), pointer :: node

    character(len=*), parameter :: h="**(FindFieldSectionAtMessageData)**"
    if (.not. associated(oneMessageData)) then
       call fatal_error(h//" oneMessageData not associated")
    end if

    node => FindFieldNamed(oneMessageData%list, fieldName)
  end function FindFieldSectionAtMessageData



  
  subroutine FillMessageDataBufferVariableAdressOneArr(oneMessageData, field)

    type(MessageData), intent(inout) :: oneMessageData
    real, target, intent(in):: field(:,:,:)

    ! copy field section values of the entire field section list of
    ! the Message Data variable to the buffer of the Message Data variable

    integer :: bufStart
    type(FieldSectionNode), pointer :: oneNode
    type(FieldSection), pointer :: oneFieldSection
    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(FillMessageDataBufferVariableAdressOneArr)**"

    if (.not. allocated(oneMessageData%buf)) then
       call fatal_error(h//" not allocated buf for message data "//&
            trim(adjustl(oneMessageData%name)))
    end if

    if (dumpLocal) then
       call MsgDump(h//" to Message Data "//trim(adjustl(oneMessageData%name)))
       write(str(1),"(i8)") lbound(field,1)
       write(str(2),"(i8)") ubound(field,1)
       write(str(3),"(i8)") lbound(field,2)
       write(str(4),"(i8)") ubound(field,2)
       write(str(5),"(i8)") lbound(field,3)
       write(str(6),"(i8)") ubound(field,3)
       call MsgDump(h//" field dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
    end if

    bufStart=1
    oneNode => FieldSectionListHeadNode(oneMessageData%list)
    oneFieldSection => FieldSectionAtNode(oneNode)
    call FieldSectionData2BufferVariableAdressArr(field, size(field,1), oneFieldSection, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
  end subroutine FillMessageDataBufferVariableAdressOneArr


  subroutine DecomposeMessageDataBufferVariableAdressOneArr(oneMessageData, &
       msgStartZ, msgEndZ, field)
    type(MessageData), intent(inout) :: oneMessageData
    integer, intent(in) :: msgStartZ
    integer, intent(in) :: msgEndZ
    real, target, intent(in):: field(:,:,:)

    ! copy all field section values from the buffer of the Message Data variable
    ! to the field section pointed by the Message Data field section list 

    integer :: bufStart
    type(FieldSectionNode), pointer :: oneNode
    type(FieldSection), pointer :: oneSection
    real, pointer :: fieldPtr(:,:,:)
    character(len=*), parameter :: h="**(DecomposeMessageDataBufferVariableAdressOneArr)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//"  of Message Data "//trim(adjustl(oneMessageData%name)))
    end if
    fieldPtr => field
    bufStart=1
    oneNode => FieldSectionListHeadNode(oneMessageData%list)
    oneSection => FieldSectionAtNode(oneNode)
    call Buffer2FieldSectionData(oneSection, &
         fieldPtr, msgStartZ, msgEndZ, &
         oneMessageData%buf, bufStart, oneMessageData%bufsize)
  end subroutine DecomposeMessageDataBufferVariableAdressOneArr
end module ModMessageData
