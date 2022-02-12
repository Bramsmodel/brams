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
       NextFieldSectionNodeAtList
  
  implicit none
  include "mpif.h"
  
  private
  public :: MessageData
  public :: CreateMessageData
  public :: DestroyMessageData
  public :: DumpMessageData
  public :: AppendFieldSectionToMessageData
  public :: AllocateMessageDataBuffer
  public :: FillMessageDataBufferWithFieldSectionData
  public :: PostSendMessageData
  public :: PostRecvMessageData
  public :: ExtractFieldSectionDataFromMessageDataBuffer
  public :: DeallocateMessageDataBuffer
  
  type MessageData
     private
     ! data to communicate to one node in one message

     ! message communicates all field values stored
     ! at the field section list

     ! message send/receive the buffer "buf"

     ! prior to a send communication, procedure
     ! FillMessageDataBufferWithFieldSectionData
     ! copies field values stored at all field sections
     ! on the field section list to the buffer

     ! after a receive communication, procedure
     ! ExtractFieldSectionDataFromMessageDataBuffer
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
    logical, parameter :: dumpLocal=.true.

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
    character(len=*), parameter :: h="**(DumpMessageData)**"

    if (present(strMsg)) then
       msgHead=strMsg
    else
       msgHead=h
    end if
    if (oneMessageData%name == "") then
       call MsgDump(trim(adjustl(msgHead))//" empty message data")
    else
       call MsgDump(trim(adjustl(msgHead))//" "//&
            trim(adjustl(oneMessageData%direction))//" message named "//&
            trim(adjustl(oneMessageData%name))//" with field section list")
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
    logical, parameter :: dumpLocal=.true.

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
    logical, parameter :: dumpLocal=.true.

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





  
  subroutine FillMessageDataBufferWithFieldSectionData(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    ! copy field section values of the entire field section list of
    ! the Message Data variable to the buffer of the Message Data variable

    integer :: bufStart
    type(FieldSectionNode), pointer :: thisNode
    type(FieldSection), pointer :: this
    logical, parameter :: dumpLocal=.true.
    character(len=*), parameter :: h="**(FillMessageDataBufferWithFieldSectionData)**"
    
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
  end subroutine FillMessageDataBufferWithFieldSectionData



  

  subroutine ExtractFieldSectionDataFromMessageDataBuffer(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    ! copy all field section values from the buffer of the Message Data variable
    ! to the field section pointed by the Message Data field section list 

    integer :: bufStart
    type(FieldSectionNode), pointer :: thisNode
    type(FieldSection), pointer :: this
    character(len=*), parameter :: h="**(ExtractFieldSectionDataFromMessageDataBuffer)**"
    logical, parameter :: dumpLocal=.true.

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
  end subroutine ExtractFieldSectionDataFromMessageDataBuffer



  


  subroutine DeallocateMessageDataBuffer(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    ! deallocates buf componenet of a Message Set variable

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DeallocateMessageDataBuffer)**"
    logical, parameter :: dumpLocal=.true.

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
    logical, parameter :: dumpLocal=.true.
    
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
    logical, parameter :: dumpLocal=.true.
    
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
end module ModMessageData
