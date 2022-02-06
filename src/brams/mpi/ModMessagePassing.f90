module ModMessagePassing

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
       MsgDump

  use ModMessageData, only: &
       MessageData, &
       AllocateMessageDataBuffer, &
       FillMessageDataBufferWithFieldSectionData, &
       ExtractFieldSectionDataFromMessageDataBuffer, &
       DeallocateMessageDataBuffer

  use ModMessageSet, only: &
       MessageSet
  
  use ModFieldSection, only: &
       FieldSection
       
  use ParLib, only: &
       parf_get_noblock_real, &
       parf_send_noblock_real, &
       parf_wait_any_nostatus

  implicit none
  include "mpif.h"
  include "ranks.h" ! for kind=i8

  private
  public :: PostRecvSendMsgs
  public :: WaitSendRecvMsgs


contains




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
    integer(kind=i8), parameter :: huge_i4=int(huge(1),kind=i8)
    integer :: bufSize_i4
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

          if (RecvMsg%msgData(iRecv)%bufSize > huge_i4) then
             write(c0,"(i8)") RecvMsg%msgData(iRecv)%bufSize
             write(c1,"(i8)") huge_i4
             call fatal_error(h//" receive buffer size ("//&
                  trim(adjustl(c0))//") cannot be represented as default kind,"//&
                  " since default is limited to "//trim(adjustl(c1)))
          else
             bufSize_i4=int(RecvMsg%msgData(iRecv)%bufSize)
          end if
          call parf_get_noblock_real(&
               RecvMsg%msgData(iRecv)%buf, &
               bufSize_i4, &
               RecvMsg%otherProc(iRecv), &
               RecvMsg%tag, &
               RecvMsg%request(iRecv))

          if (dumpLocal) then
             write(c0,"(i8)") iRecv
             write(c1,"(i8)") RecvMsg%otherProc(iRecv)
             write(c2,"(i8)") size(RecvMsg%msgData(iRecv)%buf)
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

          ! allocate and fill send buffer with field sections to send

          call AllocateMessageDataBuffer(SendMsg%msgData(iSend))
          call FillMessageDataBufferWithFieldSectionData(SendMsg%msgData(iSend))

          ! post send message

          if (SendMsg%msgData(iSend)%bufSize > huge_i4) then
             write(c0,"(i8)") SendMsg%msgData(iSend)%bufSize
             write(c1,"(i8)") huge_i4
             call fatal_error(h//" send buffer size ("//&
                  trim(adjustl(c0))//") cannot be represented as default kind,"//&
                  " since default is limited to "//trim(adjustl(c1)))
          else
             bufSize_i4=int(SendMsg%msgData(iSend)%bufSize)
          end if
          call parf_send_noblock_real(&
               SendMsg%msgData(iSend)%buf, &
               bufSize_i4, &
               SendMsg%otherProc(iSend), &
               SendMsg%tag, &
               SendMsg%request(iSend))
          if (dumpLocal) then
             write(c1,"(i8)") SendMsg%otherProc(iSend)
             write(c2,"(i8)") size(SendMsg%msgData(iSend)%buf)
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
end module ModMessagePassing
