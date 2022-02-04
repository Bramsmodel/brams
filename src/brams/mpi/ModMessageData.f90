module ModMessageData
  use ModParallelEnvironment, only: &
       MsgDump

  use ModFieldSection, only: &
       FieldSection, &
       AppendFieldSection, &
       FieldSectionSize, &
       DestroyFieldSection, &
       FieldSectionData2Buffer, &
       Buffer2FieldSectionData, &
       NextFieldSection
  

  implicit none
  include "ranks.h" ! for kind=i8
  
  private
  public :: MessageData
  public :: InitializeMessageData
  public :: AppendFieldSectionToMessageData
  public :: DestroyMessageData
  public :: AllocateMessageDataBuffer
  public :: FillMessageDataBufferWithFieldSectionData
  public :: ExtractFieldSectionDataFromMessageDataBuffer
  public :: DeallocateMessageDataBuffer
  
  ! data to send/receive to/from one node in one message

  type MessageData
     integer(kind=i8) :: bufSize=0_i8         ! message buffer size
     real, allocatable :: buf(:)  ! message buffer
     type (FieldSection), pointer :: head => null()
     type (FieldSection), pointer :: tail => null()
  end type MessageData

contains






  subroutine InitializeMessageData (MsgData)

    ! Create MessageData components of each entry
    ! of a MesageData 

    type(MessageData), intent(inout) :: MsgData

    character(len=*), parameter :: h="**(InitializeMessageData)**"
    logical, parameter :: dumpLocal=.true.
    
    MsgData%bufSize = 0
    MsgData%head => null()
    MsgData%tail => null()
    if (dumpLocal) then
       call MsgDump(h//" done")
    end if
  end subroutine InitializeMessageData






  subroutine DestroyMessageData(msgData)

    ! destroy a variable of this type

    type(MessageData) :: msgData

    integer :: ierr
    type(FieldSection), pointer:: this
    type(FieldSection), pointer:: next
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(DestroyMessageData)**"
    logical, parameter :: dumpLocal=.true.

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
    this => msgData%head
    do while (associated(this))
       next => this%next
       call DestroyFieldSection(this)
       this => next
    end do
    nullify(msgData%head)
    nullify(msgData%tail)
  end subroutine DestroyMessageData




  subroutine AppendFieldSectionToMessageData (oneFieldSection, oneMessageData)

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    type(MessageData), intent(inout) :: oneMessageData

    type(FieldSection), pointer :: previous
    character(len=*), parameter :: h="**(AppendFieldSectionToMessageData)**"

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" oneFieldSection not associated")
    end if

    if (.not. associated(oneMessageData%head)) then
       oneMessageData%head => oneFieldSection
    else
       call AppendFieldSection(oneMessageData%tail, oneFieldSection)
    end if
    oneMessageData%tail => oneFieldSection

    oneMessageData%bufSize = oneMessageData%bufSize + FieldSectionSize(oneFieldSection)
  end subroutine AppendFieldSectionToMessageData





  
  subroutine FillMessageDataBufferWithFieldSectionData(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    integer(kind=i8) :: bufStart
    type(FieldSection), pointer :: this

    bufStart=1_i8
    this => oneMessageData%head
    do while (associated(this))
       call FieldSectionData2Buffer(&
            this, &
            oneMessageData%buf, &
            bufStart, &
            oneMessageData%bufsize)
       this => NextFieldSection(this)
    end do
  end subroutine FillMessageDataBufferWithFieldSectionData



  

  subroutine ExtractFieldSectionDataFromMessageDataBuffer(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    integer(kind=i8) :: bufStart
    type(FieldSection), pointer :: this

    bufStart=1_i8
    this => oneMessageData%head
    do while (associated(this))
       call Buffer2FieldSectionData(&
            this, &
            oneMessageData%buf, &
            bufStart, &
            oneMessageData%bufsize)
       this => NextFieldSection(this)
    end do
  end subroutine ExtractFieldSectionDataFromMessageDataBuffer



  


  subroutine AllocateMessageDataBuffer(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    integer :: bufSize
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(AllocateMessageDataBuffer)**"

    bufSize=oneMessageData%bufSize

    allocate(oneMessageData%buf(bufSize), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") bufSize
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocation of oneMessageData("//&
            trim(adjustl(c0))//") fails with stat="//&
            trim(adjustl(c1)))
    end if
  end subroutine AllocateMessageDataBuffer



  


  subroutine DeallocateMessageDataBuffer(oneMessageData)
    type(MessageData), intent(inout) :: oneMessageData

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DeallocateMessageDataBuffer)**"


    deallocate(oneMessageData%buf, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocation of oneMessageData"//&
            " fails with stat="//trim(adjustl(c0)))
    end if
  end subroutine DeallocateMessageDataBuffer
end module ModMessageData
