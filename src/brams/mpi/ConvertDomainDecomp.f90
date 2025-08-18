module ModConvertDomainDecomp

  ! MPI send/receive operations to convert a field
  ! from one domain decomposition to other domain decomposition

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModGridDims, only: &
       GridDims
  
  use ModNeighbourNodes, only: &
       Inter
  
  use ModDomainDecomp, only: &
       DomainDecomp

  use ModNodeDimensions, only: &
       NodeDimensions
  
  implicit none

  private

  include "mpif.h"
  
  ! ConvertDomainDecomp: send and receive data for MPI
  !                      operations to convert field from
  !                      one domain decomposition (From)
  !                      to other domain decomposition (To)
  !
  ! arrays are indexed from 1 to nSend or nRecv
  ! destRankSend and srcRankRecv contain MPI ranks
  ! xb, xe, yb, ye contain local indices
  
  type ConvertDomainDecomp
     character(len=64) :: FromName
     character(len=64) :: ToName
     character(len=64) :: VarName
     integer :: nSend ! count send messages
     integer, pointer, contiguous :: destRankSend(:) ! MPI rank to receive each message
     integer, pointer, contiguous :: xbSend(:) ! x send region first local index
     integer, pointer, contiguous :: xeSend(:) ! x send region last local index
     integer, pointer, contiguous :: ybSend(:) ! y send region first local index
     integer, pointer, contiguous :: yeSend(:) ! y send region last local index
     integer :: nRecv ! count recv messages
     integer, pointer, contiguous :: srcRankRecv(:) ! MPI rank that send each message
     integer, pointer, contiguous :: xbRecv(:) ! x recv region first local index
     integer, pointer, contiguous :: xeRecv(:) ! x recv region last local index
     integer, pointer, contiguous :: ybRecv(:) ! y recv region first local index
     integer, pointer, contiguous :: yeRecv(:) ! y recv region last local index
  end type ConvertDomainDecomp

  integer, parameter :: tag = 98
  integer, parameter :: comm = MPI_COMM_WORLD

  interface SendRecvConvertDomainDecomp
     module procedure SendRecvConvertDomainDecomp_2D
     module procedure SendRecvConvertDomainDecomp_3D
  end interface SendRecvConvertDomainDecomp
  
  public :: ConvertDomainDecomp
  public :: CreateConvertDomainDecomp
  public :: DestroyConvertDomainDecomp
  public :: DumpConvertDomainDecomp
  public :: SendRecvConvertDomainDecomp
  public :: TestSendRecvConvertDomainDecomp
contains





  function CreateConvertDomainDecomp(oneParallelEnvironment, oneGridDims, &
       GlobalWithGhostFrom, FromName, GlobalWithGhostTo, ToName, &
       varName) result(oneConvertDomainDecomp)
    type(ParallelEnvironment), pointer, intent(in) :: oneParallelEnvironment
    type(GridDims), pointer, intent(in) :: oneGridDims
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhostFrom
    character(len=*), intent(in) :: FromName
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhostTo
    character(len=*), intent(in) :: ToName
    character(len=*), intent(in) :: varName
    type(ConvertDomainDecomp), pointer :: oneConvertDomainDecomp

    integer :: ierr
    integer :: rank
    integer :: nRanks ! how many processes at this computation
    integer :: myRank ! this rank at Brams enumeration
    integer :: xbOwns
    integer :: xeOwns
    integer :: ybOwns
    integer :: yeOwns
    integer :: xbInter ! x begin intersection
    integer :: xeInter ! x end intersection
    integer :: ybInter ! y begin intersection
    integer :: yeInter ! y end intersection
    logical :: hasInter ! has intersection or not
    character(len=512) :: message

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateConvertDomainDecomp)**"
    logical, parameter :: dumpLocal=.false.

    ! count processes and this process BRAMS rank
    
    nRanks=oneParallelEnvironment%nmachs
    myRank=oneParallelEnvironment%myNum

    if (dumpLocal) then
       call MsgDump(h//" starts for variable "//trim(varName)//&
            " that converts from "//trim(FromName)//&
            " to "//trim(ToName))
    end if
    
    ! otherwise, allocate result and their pointer components

    allocate(oneConvertDomainDecomp, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneConvertDomainDecomp fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%destRankSend(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate destRankSend fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%xbSend(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate xbSend fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%xeSend(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate xeSend fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%ybSend(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ybSend fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%yeSend(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate yeSend fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%srcRankRecv(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate destRankRecv fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%xbRecv(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate xbRecv fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%xeRecv(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate xeRecv fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%ybRecv(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ybRecv fails with message "//trim(message))
    end if

    allocate(oneConvertDomainDecomp%yeRecv(nRanks), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate yeRecv fails with message "//trim(message))
    end if

    ! initialize scalar components

    oneConvertDomainDecomp%FromName = trim(FromName)
    oneConvertDomainDecomp%ToName = trim(ToName)
    oneConvertDomainDecomp%VarName = trim(varName)
    oneConvertDomainDecomp%nSend=0
    oneConvertDomainDecomp%nRecv=0

    ! sections to send from this node to all nodes

    ! this node should send all points that it owns,
    ! which are GlobalWithGhost without ghost zone
    ! but including boundary conditions

    xbOwns = GlobalWithGhostFrom%xb(myRank)
    xeOwns = GlobalWithGhostFrom%xe(myRank)
    ybOwns = GlobalWithGhostFrom%yb(myRank)
    yeOwns = GlobalWithGhostFrom%ye(myRank)
    ! do not send ghost zone in Brams Domain Decomposition
    if (GlobalWithGhostFrom%FullDirection == "B") then
       if (xbOwns /= 1) then
          xbOwns = xbOwns + 1
       end if
       if (xeOwns /= oneGridDims%nnxp) then
          xeOwns = xeOwns - 1
       end if
       if (ybOwns /= 1) then
          ybOwns = ybOwns + 1
       end if
       if (yeOwns /= oneGridDims%nnyp) then
          yeOwns = yeOwns - 1
       end if
    end if
    do rank = 1, oneParallelEnvironment%nmachs

       ! intersection of points owned by this node "From" domain decomposition
       ! with points owned by node "rank" in the "To" domain decomposition;
       ! intersection is computed in global indices
       
       call Inter(&
            xbOwns, xeOwns, &
            ybOwns, yeOwns, &
            GlobalWithGhostTo%xb(rank), GlobalWithGhostTo%xe(rank), &
            GlobalWithGhostTo%yb(rank), GlobalWithGhostTo%ye(rank), &
            xbInter, xeInter, ybInter, yeInter, &
            hasInter)

       ! if there is intersection, store domain region,
       ! converting global indices into local indices and
       ! Brams rank enumeration into MPI rank enumeration
       
       if (hasInter) then
          oneConvertDomainDecomp%nSend = oneConvertDomainDecomp%nSend+1
          oneConvertDomainDecomp%xbSend(oneConvertDomainDecomp%nSend) = xbInter+1-GlobalWithGhostFrom%xb(myRank)
          oneConvertDomainDecomp%xeSend(oneConvertDomainDecomp%nSend) = xeInter+1-GlobalWithGhostFrom%xb(myRank)
          oneConvertDomainDecomp%ybSend(oneConvertDomainDecomp%nSend) = ybInter+1-GlobalWithGhostFrom%yb(myRank)
          oneConvertDomainDecomp%yeSend(oneConvertDomainDecomp%nSend) = yeInter+1-GlobalWithGhostFrom%yb(myRank)
          oneConvertDomainDecomp%destRankSend(oneConvertDomainDecomp%nSend) = rank-1

          if (dumpLocal) then
             write(str(1),"(i8)") xbInter
             write(str(2),"(i8)") xeInter
             write(str(3),"(i8)") ybInter
             write(str(4),"(i8)") yeInter
             write(str(5),"(i8)") myRank-1
             write(str(6),"(i8)") rank-1
             call MsgDump(h//" MPI rank "//trim(adjustl(str(5)))//" sends global ["//&
                  trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
                  trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//"]"//&
                  " to MPI rank "//trim(adjustl(str(6))))
             write(str(1),"(i8)") oneConvertDomainDecomp%xbSend(oneConvertDomainDecomp%nSend)
             write(str(2),"(i8)") oneConvertDomainDecomp%xeSend(oneConvertDomainDecomp%nSend)
             write(str(3),"(i8)") oneConvertDomainDecomp%ybSend(oneConvertDomainDecomp%nSend)
             write(str(4),"(i8)") oneConvertDomainDecomp%yeSend(oneConvertDomainDecomp%nSend)
             call MsgDump(h//" these are locals ["//&
                  trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
                  trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//"]")
          end if
       end if
    end do

    ! sections received by this node from all nodes

    ! this node will receive all points that it owns in "To" domain decomposition,
    ! including boundary conditions and ghost zone, from all ranks in the
    ! "From" domain decomposition excluding ghost zone

    do rank = 1, oneParallelEnvironment%nmachs

       ! the other node should send all points that it owns,
       ! which are GlobalWithGhost without ghost zone
       ! but including boundary conditions

       xbOwns = GlobalWithGhostFrom%xb(rank)
       xeOwns = GlobalWithGhostFrom%xe(rank)
       ybOwns = GlobalWithGhostFrom%yb(rank)
       yeOwns = GlobalWithGhostFrom%ye(rank)
       ! do not send ghost zone in Brams Domain Decomposition
       if (GlobalWithGhostFrom%FullDirection == "B") then
          if (xbOwns /= 1) then
             xbOwns = xbOwns + 1
          end if
          if (xeOwns /= oneGridDims%nnxp) then
             xeOwns = xeOwns - 1
          end if
          if (ybOwns /= 1) then
             ybOwns = ybOwns + 1
          end if
          if (yeOwns /= oneGridDims%nnyp) then
             yeOwns = yeOwns - 1
          end if
       end if
       
       ! intersection of points owned by this node "To" domain decomposition
       ! with points owned by node "rank" in the "From" domain decomposition
       
       call Inter(&
            xbOwns, xeOwns, &
            ybOwns, yeOwns, &
            GlobalWithGhostTo%xb(myRank), GlobalWithGhostTo%xe(myRank), &
            GlobalWithGhostTo%yb(myRank), GlobalWithGhostTo%ye(myRank), &
            xbInter, xeInter, ybInter, yeInter, &
            hasInter)

       ! if there is intersection, store domain region,
       ! converting global indices into local indices and
       ! Brams rank enumeration into MPI rank enumeration
       
       if (hasInter) then
          oneConvertDomainDecomp%nRecv = oneConvertDomainDecomp%nRecv+1
          oneConvertDomainDecomp%xbRecv(oneConvertDomainDecomp%nRecv) = xbInter+1-GlobalWIthGhostTo%xb(myRank)
          oneConvertDomainDecomp%xeRecv(oneConvertDomainDecomp%nRecv) = xeInter+1-GlobalWIthGhostTo%xb(myRank)
          oneConvertDomainDecomp%ybRecv(oneConvertDomainDecomp%nRecv) = ybInter+1-GlobalWIthGhostTo%yb(myRank)
          oneConvertDomainDecomp%yeRecv(oneConvertDomainDecomp%nRecv) = yeInter+1-GlobalWIthGhostTo%yb(myRank)
          oneConvertDomainDecomp%srcRankRecv(oneConvertDomainDecomp%nRecv) = rank-1
               
          if (dumpLocal) then
             write(str(1),"(i8)") xbInter
             write(str(2),"(i8)") xeInter
             write(str(3),"(i8)") ybInter
             write(str(4),"(i8)") yeInter
             write(str(5),"(i8)") myRank-1
             write(str(6),"(i8)") rank-1
             call MsgDump(h//" MPI rank "//trim(adjustl(str(5)))//" receives globals ["//&
                  trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
                  trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//"]"//&
                  " from MPI rank "//trim(adjustl(str(6))))
             write(str(1),"(i8)") oneConvertDomainDecomp%xbRecv(oneConvertDomainDecomp%nRecv)
             write(str(2),"(i8)") oneConvertDomainDecomp%xeRecv(oneConvertDomainDecomp%nRecv)
             write(str(3),"(i8)") oneConvertDomainDecomp%ybRecv(oneConvertDomainDecomp%nRecv)
             write(str(4),"(i8)") oneConvertDomainDecomp%yeRecv(oneConvertDomainDecomp%nRecv)
             call MsgDump(h//" these are locals ["//&
                  trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
                  trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//"]")
          end if
       end if
    end do

    if (dumpLocal) then
       call DumpConvertDomainDecomp(oneConvertDomainDecomp)
    end if
  end function CreateConvertDomainDecomp 



  

  subroutine DestroyConvertDomainDecomp(oneConvertDomainDecomp)
    type(ConvertDomainDecomp), pointer, intent(inout) :: oneConvertDomainDecomp

    integer :: ierr
    character(len=512) :: message
    
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyConvertDomainDecomp)**"
    logical, parameter :: dumpLocal=.false.

    ! nothing to do if null pointer

    if (.not. associated(oneConvertDomainDecomp)) then
       if (dumpLocal) then
          call MsgDump(h//" input pointer not associated")
       end if
       return
    end if
    
    ! otherwise, deallocate result and their pointer components

    deallocate(oneConvertDomainDecomp%destRankSend, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate destRankSend fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp%xbSend, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate xbSend fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp%xeSend, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate xeSend fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp%ybSend, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate ybSend fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp%yeSend, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate yeSend fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp%srcRankRecv, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate destRankRecv fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp%xbRecv, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate xbRecv fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp%xeRecv, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate xeRecv fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp%ybRecv, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate ybRecv fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp%yeRecv, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate yeRecv fails with message "//trim(message))
    end if

    deallocate(oneConvertDomainDecomp, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate oneConvertDomainDecomp fails with message "//trim(message))
    end if

    nullify(oneConvertDomainDecomp)
  end subroutine DestroyConvertDomainDecomp





  subroutine DumpConvertDomainDecomp(oneConvertDomainDecomp)
    type(ConvertDomainDecomp), pointer, intent(in) :: oneConvertDomainDecomp

    integer :: i
    
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpConvertDomainDecomp)**"

    if (.not. associated(oneConvertDomainDecomp)) then
       call MsgDump(h//" empty ")
       return
    end if

    call MsgDump(h//" variable "//trim(oneConvertDomainDecomp%VarName)//&
         " converts fields from domain "//&
         trim(oneConvertDomainDecomp%FromName)//" to domain "//&
         trim(oneConvertDomainDecomp%ToName))

    if (oneConvertDomainDecomp%nSend == 0) then
       call MsgDump(h//" no send operations")
    else
       write(str(1),"(i8)") oneConvertDomainDecomp%nSend
       call MsgDump(h//" conversion has "//trim(adjustl(str(1)))//" send operations:")
    end if

    do i = 1, oneConvertDomainDecomp%nSend
       write(str(1),"(i8)") oneConvertDomainDecomp%xbSend(i)
       write(str(2),"(i8)") oneConvertDomainDecomp%xeSend(i)
       write(str(3),"(i8)") oneConvertDomainDecomp%ybSend(i)
       write(str(4),"(i8)") oneConvertDomainDecomp%yeSend(i)
       write(str(5),"(i8)") oneConvertDomainDecomp%destRankSend(i)
       call MsgDump(h//" send local indices ["//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//"]"//&
            " to MPI rank "//trim(adjustl(str(5))))
    end do

    if (oneConvertDomainDecomp%nRecv == 0) then
       call MsgDump(h//" no receive operations")
    else
       write(str(1),"(i8)") oneConvertDomainDecomp%nRecv
       call MsgDump(h//" conversion has "//trim(adjustl(str(1)))//" receive operations:")
    end if

    do i = 1, oneConvertDomainDecomp%nRecv
       write(str(1),"(i8)") oneConvertDomainDecomp%xbRecv(i)
       write(str(2),"(i8)") oneConvertDomainDecomp%xeRecv(i)
       write(str(3),"(i8)") oneConvertDomainDecomp%ybRecv(i)
       write(str(4),"(i8)") oneConvertDomainDecomp%yeRecv(i)
       write(str(5),"(i8)") oneConvertDomainDecomp%srcRankRecv(i)
       call MsgDump(h//" receive local indices ["//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//"]"//&
            " from MPI rank "//trim(adjustl(str(5))))
    end do
  end subroutine DumpConvertDomainDecomp





  subroutine SendRecvConvertDomainDecomp_2D(fieldSend, fieldRecv, oneConvertDomainDecomp)

    ! posts all nonblocking send and recv operations of
    ! a message set pair of variables

    real, pointer, contiguous, intent(in) :: fieldSend(:,:)
    ! pointer intent(in), values intent(in)
    real, pointer, contiguous, intent(in) :: fieldRecv(:,:)
    ! pointer intent(in), values intent(out)
    type(ConvertDomainDecomp), pointer, intent(in) :: oneConvertDomainDecomp
    ! pointer intent(in), values intent(inout)

    type Buff
       real, allocatable :: buffer(:)
    end type Buff

    type(Buff), allocatable :: allBuffSend(:)
    type(Buff), allocatable :: allBuffRecv(:)
    
    integer :: iSend
    integer :: nSend
    integer :: iRecv
    integer :: nRecv
    integer :: dataSize
    integer :: ierr
    integer :: ierr2
    integer :: messageLen
    integer :: x
    integer :: y
    integer :: cnt
    integer, allocatable :: requestRecv(:)
    integer, allocatable :: statusRecv(:,:)
    character(len=512) :: message
    
    character(len=8) :: str(10)
    character(len=16) :: strL(10)
    character(len=*), parameter :: h="**(SendRecvConvertDomainDecomp_2D)**"
    logical, parameter :: dumpLocal=.false.

    ! return if empty conversion
    
    if (.not. associated(oneConvertDomainDecomp)) then
       if (dumpLocal) then
          call MsgDump(h//" oneConvertDomainDecomp not associated")
       end if
       return
    end if

    nSend = oneConvertDomainDecomp%nSend
    nRecv = oneConvertDomainDecomp%nRecv

    if (nSend > 0) then
       allocate(allBuffSend(nSend), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate allBuffSend fails with message "//trim(message))
       end if
    end if
    if (nRecv > 0) then
       allocate(allBuffRecv(nRecv), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate allBuffRecv fails with message "//trim(message))
       end if
       allocate(requestRecv(nRecv), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate requestRecv fails with message "//trim(message))
       end if
    end if
    
    ! post nonblocking receive for each receiving message

    if (dumpLocal) then
       write(str(1),"(i8)") nRecv
       call MsgDump(h//" for "//trim(oneConvertDomainDecomp%VarName)//&
            " will post "//trim(adjustl(str(1)))//&
            " nonblocking receives")
    end if
    do iRecv= 1, nRecv
       
       dataSize = &
            (oneConvertDomainDecomp%xeRecv(iRecv) - oneConvertDomainDecomp%xbRecv(iRecv) + 1) * &
            (oneConvertDomainDecomp%yeRecv(iRecv) - oneConvertDomainDecomp%ybRecv(iRecv) + 1)

       allocate(allBuffRecv(iRecv)%buffer(dataSize), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate allBuffRecv%buffer fails with message "//trim(message))
       end if

       call MPI_Irecv(&
            allBuffRecv(iRecv)%buffer, &
            dataSize, &
            MPI_REAL, &
            oneConvertDomainDecomp%srcRankRecv(iRecv), &
            tag, &
            comm, &
            requestRecv(iRecv), &
            ierr)
       if (ierr /= MPI_SUCCESS) then
          write(str(1),"(i8)") oneConvertDomainDecomp%xbRecv(iRecv)
          write(str(2),"(i8)") oneConvertDomainDecomp%xeRecv(iRecv)
          write(str(3),"(i8)") oneConvertDomainDecomp%ybRecv(iRecv)
          write(str(4),"(i8)") oneConvertDomainDecomp%yeRecv(iRecv)
          call MPI_Error_string(ierr, message, messageLen, ierr2)
          call fatal_error(h//" Irecv convertion "//trim(oneConvertDomainDecomp%VarName)//&
               " for variable "//trim(oneConvertDomainDecomp%ToName)//"("//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//") fails with message "//&
               message(1:messageLen))
       end if
       
       if (dumpLocal) then
          write(str(1),"(i8)") oneConvertDomainDecomp%srcRankRecv(iRecv)
          write(str(3),"(i8)") oneConvertDomainDecomp%xbRecv(iRecv)
          write(str(4),"(i8)") oneConvertDomainDecomp%xeRecv(iRecv)
          write(str(5),"(i8)") oneConvertDomainDecomp%ybRecv(iRecv)
          write(str(6),"(i8)") oneConvertDomainDecomp%yeRecv(iRecv)
          call MsgDump(h//" will recv from rank "//trim(adjustl(str(1)))//&
               " "//trim(oneConvertDomainDecomp%ToName)//"("//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       end if

    end do

    ! post blocking send for each send message

    if (dumpLocal) then
       write(str(1),"(i8)") nSend
       call MsgDump(h//" for "//trim(oneConvertDomainDecomp%VarName)//&
            " will post "//trim(adjustl(str(1)))//&
            " blocking sends")
    end if
    do iSend = 1, nSend

       dataSize = &
            (oneConvertDomainDecomp%xeSend(iSend) - oneConvertDomainDecomp%xbSend(iSend) + 1) * &
            (oneConvertDomainDecomp%yeSend(iSend) - oneConvertDomainDecomp%ybSend(iSend) + 1)

       allocate(allBuffSend(iSend)%buffer(dataSize), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate allBuffSend%buffer fails with message "//trim(message))
       end if

       cnt = 1
       do y = oneConvertDomainDecomp%ybSend(iSend), oneConvertDomainDecomp%yeSend(iSend)
          do x = oneConvertDomainDecomp%xbSend(iSend), oneConvertDomainDecomp%xeSend(iSend)
             allBuffSend(iSend)%buffer(cnt) = fieldSend(x,y)
             cnt = cnt + 1
          end do
       end do

       call MPI_Send(&
            allBuffSend(iSend)%buffer, &
            dataSize, &
            MPI_REAL, &
            oneConvertDomainDecomp%destRankSend(iSend), &
            tag, &
            comm, &
            ierr)
       if (ierr /= MPI_SUCCESS) then
          write(str(1),"(i8)") oneConvertDomainDecomp%xbSend(iSend)
          write(str(2),"(i8)") oneConvertDomainDecomp%xeSend(iSend)
          write(str(3),"(i8)") oneConvertDomainDecomp%ybSend(iSend)
          write(str(4),"(i8)") oneConvertDomainDecomp%yeSend(iSend)
          call MPI_Error_string(ierr, message, messageLen, ierr2)
          call fatal_error(h//" Send convertion "//trim(oneConvertDomainDecomp%VarName)//&
               " for variable "//trim(oneConvertDomainDecomp%ToName)//"("//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//") fails with message "//&
               message(1:messageLen))
       end if

       if (dumpLocal) then
          write(str(1),"(i8)") oneConvertDomainDecomp%destRankSend(iSend)
          write(str(3),"(i8)") oneConvertDomainDecomp%xbSend(iSend)
          write(str(4),"(i8)") oneConvertDomainDecomp%xeSend(iSend)
          write(str(5),"(i8)") oneConvertDomainDecomp%ybSend(iSend)
          write(str(6),"(i8)") oneConvertDomainDecomp%yeSend(iSend)
          call MsgDump(h//" sent to rank "//trim(adjustl(str(1)))//&
               " "//trim(oneConvertDomainDecomp%FromName)//"("//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
          write(strL(1),"(f16.8)") allBuffSend(iSend)%buffer(1)
          write(strL(2),"(f16.8)") allBuffSend(iSend)%buffer(dataSize)
          call MsgDump(h//" samples: "//&
               "first value sent="//trim(adjustl(strL(1)))//&
               " last value sent="//trim(adjustl(strL(2))))
       end if

    end do
    
    ! wait for all receive messages

    if (nRecv > 0) then

       allocate(statusRecv(MPI_STATUS_SIZE, nRecv), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate statusRecv fails with message "//trim(message))
       end if
       
       if (dumpLocal) then
          write(str(1),"(i8)") nRecv
          call MsgDump(h//" will wait for "//trim(adjustl(str(1)))//" receives")
       end if
    
       call MPI_Waitall(&
            nRecv, &
            requestRecv, &
            statusRecv, &
            ierr)

       do iRecv= 1, nRecv
          cnt = 1
          do y = oneConvertDomainDecomp%ybRecv(iRecv), oneConvertDomainDecomp%yeRecv(iRecv)
             do x = oneConvertDomainDecomp%xbRecv(iRecv), oneConvertDomainDecomp%xeRecv(iRecv)
                fieldRecv(x,y) = allBuffRecv(iRecv)%buffer(cnt)
                cnt = cnt + 1
             end do
          end do

          if (dumpLocal) then
             write(str(1),"(i8)") statusRecv(MPI_SOURCE,iRecv)
             write(str(3),"(i8)") oneConvertDomainDecomp%xbRecv(iRecv)
             write(str(4),"(i8)") oneConvertDomainDecomp%xeRecv(iRecv)
             write(str(5),"(i8)") oneConvertDomainDecomp%ybRecv(iRecv)
             write(str(6),"(i8)") oneConvertDomainDecomp%yeRecv(iRecv)
             call MsgDump(h//" recv from rank "//trim(adjustl(str(1)))//&
                  " "//trim(oneConvertDomainDecomp%ToName)//"("//&
                  trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
                  trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
             dataSize=size(allBuffRecv(iRecv)%buffer)
             write(strL(1),"(f16.8)") allBuffRecv(iRecv)%buffer(1)
             write(strL(2),"(f16.8)") allBuffRecv(iRecv)%buffer(dataSize)
             write(str(1),"(i8)") dataSize
             call MsgDump(h//" recv "//trim(adjustl(str(1)))//" values; "//&
                  "first value recv="//trim(adjustl(strL(1)))//&
                  " last value recv="//trim(adjustl(strL(2))))
          end if
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") nRecv
          call MsgDump(h//" done with "//trim(adjustl(str(1)))//" receives")
       end if

       do iRecv= 1, nRecv
          deallocate(allBuffRecv(iRecv)%buffer, stat=ierr, errmsg=message)
          if (ierr /= 0) then
             call fatal_error(h//" deallocate allBuffRecv%buffer fails with message "//trim(message))
          end if
       end do
    
       deallocate(statusRecv, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate statusRecv fails with message "//trim(message))
       end if

       deallocate(allBuffRecv, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate allBuffRecv fails with message "//trim(message))
       end if
       deallocate(requestRecv, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate requestRecv fails with message "//trim(message))
       end if
    end if

    do iSend = 1, nSend
       deallocate(allBuffSend(iSend)%buffer, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate allBuffSend%buffer fails with message "//trim(message))
       end if
    end do
    
    if (nSend > 0) then
       deallocate(allBuffSend, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate allBuffSend fails with message "//trim(message))
       end if
    end if
  end subroutine SendRecvConvertDomainDecomp_2D


  

  
  subroutine SendRecvConvertDomainDecomp_3D(fieldSend, fieldRecv, oneConvertDomainDecomp)

    ! posts all nonblocking send and recv operations of
    ! a message set pair of variables

    real, pointer, contiguous, intent(in) :: fieldSend(:,:,:)
    ! pointer intent(in), values intent(in)
    real, pointer, contiguous, intent(in) :: fieldRecv(:,:,:)
    ! pointer intent(in), values intent(out)
    type(ConvertDomainDecomp), pointer, intent(in) :: oneConvertDomainDecomp
    ! pointer intent(in), values intent(inout)

    type Buff
       real, allocatable :: buffer(:)
    end type Buff

    type(Buff), allocatable :: allBuffSend(:)
    type(Buff), allocatable :: allBuffRecv(:)
    
    integer :: iSend
    integer :: nSend
    integer :: iRecv
    integer :: nRecv
    integer :: dataSize
    integer :: ierr
    integer :: ierr2
    integer :: messageLen
    integer :: x
    integer :: y
    integer :: z
    integer :: nz
    integer :: cnt
    integer, allocatable :: requestRecv(:)
    integer, allocatable :: statusRecv(:,:)
    character(len=512) :: message
    
    character(len=8) :: str(10)
    character(len=16) :: strL(10)
    character(len=*), parameter :: h="**(SendRecvConvertDomainDecomp_3D)**"
    logical, parameter :: dumpLocal=.false.

    ! return if empty conversion
    
    if (.not. associated(oneConvertDomainDecomp)) then
       if (dumpLocal) then
          call MsgDump(h//" oneConvertDomainDecomp not associated")
       end if
       return
    end if

    nz = size(fieldRecv,1)
    if (nz /= size(fieldSend,1)) then
       call fatal_error(h//" first dimension of fieldSend and fieldRecv dissagree")
    end if

    nSend = oneConvertDomainDecomp%nSend
    if (nSend > 0) then
       allocate(allBuffSend(nSend), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate allBuffSend fails with message "//trim(message))
       end if
    end if

    nRecv = oneConvertDomainDecomp%nRecv
    if (nRecv > 0) then
       allocate(allBuffRecv(nRecv), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate allBuffRecv fails with message "//trim(message))
       end if
       allocate(requestRecv(nRecv), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate requestRecv fails with message "//trim(message))
       end if
    end if
    
    ! post nonblocking receive for each receiving message

    if (dumpLocal) then
       write(str(1),"(i8)") nRecv
       call MsgDump(h//" for "//trim(oneConvertDomainDecomp%VarName)//&
            " will post "//trim(adjustl(str(1)))//&
            " nonblocking receives")
    end if
    do iRecv= 1, nRecv

       dataSize = &
            nz * &
            (oneConvertDomainDecomp%xeRecv(iRecv) - oneConvertDomainDecomp%xbRecv(iRecv) + 1) * &
            (oneConvertDomainDecomp%yeRecv(iRecv) - oneConvertDomainDecomp%ybRecv(iRecv) + 1)

       allocate(allBuffRecv(iRecv)%buffer(dataSize), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          write(str(1),"(i8)") dataSize
          call fatal_error(h//" allocate allBuffRecv%buffer("//&
               trim(adjustl(str(1)))//") fails with message "//trim(message))
       end if

       call MPI_Irecv(&
            allBuffRecv(iRecv)%buffer, &
            dataSize, &
            MPI_REAL, &
            oneConvertDomainDecomp%srcRankRecv(iRecv), &
            tag, &
            comm, &
            requestRecv(iRecv), &
            ierr)
       if (ierr /= MPI_SUCCESS) then
          write(str(1),"(i8)") oneConvertDomainDecomp%xbRecv(iRecv)
          write(str(2),"(i8)") oneConvertDomainDecomp%xeRecv(iRecv)
          write(str(3),"(i8)") oneConvertDomainDecomp%ybRecv(iRecv)
          write(str(4),"(i8)") oneConvertDomainDecomp%yeRecv(iRecv)
          write(str(5),"(i8)") nz
          call MPI_Error_string(ierr, message, messageLen, ierr2)
          call fatal_error(h//" Irecv convertion "//trim(oneConvertDomainDecomp%VarName)//&
               " for variable "//trim(oneConvertDomainDecomp%ToName)//"("//&
               "1:"//trim(adjustl(str(5)))//","//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//") fails with message "//&
               message(1:messageLen))
       end if

       if (dumpLocal) then
          write(str(1),"(i8)") oneConvertDomainDecomp%srcRankRecv(iRecv)
          write(str(2),"(i8)") nz
          write(str(3),"(i8)") oneConvertDomainDecomp%xbRecv(iRecv)
          write(str(4),"(i8)") oneConvertDomainDecomp%xeRecv(iRecv)
          write(str(5),"(i8)") oneConvertDomainDecomp%ybRecv(iRecv)
          write(str(6),"(i8)") oneConvertDomainDecomp%yeRecv(iRecv)
          call MsgDump(h//" will recv from rank "//trim(adjustl(str(1)))//&
               " "//trim(oneConvertDomainDecomp%ToName)//"("//&
               "1:"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       end if
    end do

    ! post blocking send for each send message

    if (dumpLocal) then
       write(str(1),"(i8)") nSend
       call MsgDump(h//" for "//trim(oneConvertDomainDecomp%VarName)//&
            " will post "//trim(adjustl(str(1)))//&
            " blocking sends")
    end if
    do iSend = 1, nSend

       dataSize = &
            nz * &
            (oneConvertDomainDecomp%xeSend(iSend) - oneConvertDomainDecomp%xbSend(iSend) + 1) * &
            (oneConvertDomainDecomp%yeSend(iSend) - oneConvertDomainDecomp%ybSend(iSend) + 1)

       allocate(allBuffSend(iSend)%buffer(dataSize), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate allBuffSend%buffer fails with message "//trim(message))
       end if
       

       cnt = 0
       do y = oneConvertDomainDecomp%ybSend(iSend), oneConvertDomainDecomp%yeSend(iSend)
          do x = oneConvertDomainDecomp%xbSend(iSend), oneConvertDomainDecomp%xeSend(iSend)
             allBuffSend(iSend)%buffer(cnt+1:cnt+nz) = fieldSend(1:nz,x,y)
             cnt = cnt + nz
          end do
       end do

       call MPI_Send(&
            allBuffSend(iSend)%buffer, &
            dataSize, &
            MPI_REAL, &
            oneConvertDomainDecomp%destRankSend(iSend), &
            tag, &
            comm, &
            ierr)
       if (ierr /= MPI_SUCCESS) then
          write(str(1),"(i8)") oneConvertDomainDecomp%xbSend(iSend)
          write(str(2),"(i8)") oneConvertDomainDecomp%xeSend(iSend)
          write(str(3),"(i8)") oneConvertDomainDecomp%ybSend(iSend)
          write(str(4),"(i8)") oneConvertDomainDecomp%yeSend(iSend)
          call MPI_Error_string(ierr, message, messageLen, ierr2)
          call fatal_error(h//" Send convertion "//trim(oneConvertDomainDecomp%VarName)//&
               " for variable "//trim(oneConvertDomainDecomp%ToName)//"("//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//") fails with message "//&
               message(1:messageLen))
       end if

       if (dumpLocal) then
          write(str(1),"(i8)") oneConvertDomainDecomp%destRankSend(iSend)
          write(str(2),"(i8)") nz
          write(str(3),"(i8)") oneConvertDomainDecomp%xbSend(iSend)
          write(str(4),"(i8)") oneConvertDomainDecomp%xeSend(iSend)
          write(str(5),"(i8)") oneConvertDomainDecomp%ybSend(iSend)
          write(str(6),"(i8)") oneConvertDomainDecomp%yeSend(iSend)
          call MsgDump(h//" sent to rank "//trim(adjustl(str(1)))//&
               " "//trim(oneConvertDomainDecomp%FromName)//"("//&
               "1:"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
          write(strL(1),"(f16.8)") allBuffSend(iSend)%buffer(1)
          write(strL(2),"(f16.8)") allBuffSend(iSend)%buffer(dataSize)
          call MsgDump(h//" samples: "//&
               "first value sent="//trim(adjustl(strL(1)))//&
               " last value sent="//trim(adjustl(strL(2))))
       end if

    end do
    
    if (dumpLocal) then
       call MsgDump(h//" all sends finished")
    end if

    ! wait for all receive messages

    if (dumpLocal) then
       write(str(1),"(i8)") nRecv
       call MsgDump(h//" will wait for "//trim(adjustl(str(1)))//&
            " nonblocking receives")
    end if
    
    if (nRecv > 0) then

       allocate(statusRecv(MPI_STATUS_SIZE, nRecv), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate statusRecv fails with message "//trim(message))
       end if
       
       if (dumpLocal) then
          write(str(1),"(i8)") nRecv
          call MsgDump(h//" will wait for "//trim(adjustl(str(1)))//" receives")
       end if
       
       call MPI_Waitall(&
            nRecv, &
            requestRecv, &
            statusRecv, &
            ierr)

       do iRecv= 1, nRecv
          cnt = 0
          do y = oneConvertDomainDecomp%ybRecv(iRecv), oneConvertDomainDecomp%yeRecv(iRecv)
             do x = oneConvertDomainDecomp%xbRecv(iRecv), oneConvertDomainDecomp%xeRecv(iRecv)
                fieldRecv(1:nz,x,y) = allBuffRecv(iRecv)%buffer(cnt+1:cnt+nz)
                cnt = cnt + nz
             end do
          end do

          if (dumpLocal) then
             write(str(1),"(i8)") statusRecv(MPI_SOURCE,iRecv)
             write(str(2),"(i8)") nz
             write(str(3),"(i8)") oneConvertDomainDecomp%xbRecv(iRecv)
             write(str(4),"(i8)") oneConvertDomainDecomp%xeRecv(iRecv)
             write(str(5),"(i8)") oneConvertDomainDecomp%ybRecv(iRecv)
             write(str(6),"(i8)") oneConvertDomainDecomp%yeRecv(iRecv)
             call MsgDump(h//" recv from rank "//trim(adjustl(str(1)))//&
                  " "//trim(oneConvertDomainDecomp%ToName)//"("//&
                  "1:"//trim(adjustl(str(2)))//","//&
                  trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
                  trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
             dataSize=size(allBuffRecv(iRecv)%buffer)
             write(strL(1),"(f16.8)") allBuffRecv(iRecv)%buffer(1)
             write(strL(2),"(f16.8)") allBuffRecv(iRecv)%buffer(dataSize)
             write(str(1),"(i8)") dataSize
             call MsgDump(h//" recv "//trim(adjustl(str(1)))//" values; "//&
                  "first value recv="//trim(adjustl(strL(1)))//&
                  " last value recv="//trim(adjustl(strL(2))))
          end if
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") nRecv
          call MsgDump(h//" done with "//trim(adjustl(str(1)))//" receives")
       end if

       do iRecv= 1, nRecv
          deallocate(allBuffRecv(iRecv)%buffer, stat=ierr, errmsg=message)
          if (ierr /= 0) then
             call fatal_error(h//" deallocate allBuffRecv%buffer fails with message "//trim(message))
          end if
       end do
    
       deallocate(statusRecv, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate statusRecv fails with message "//trim(message))
       end if

       deallocate(allBuffRecv, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate allBuffRecv fails with message "//trim(message))
       end if
       deallocate(requestRecv, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate requestRecv fails with message "//trim(message))
       end if
    end if

    do iSend = 1, nSend
       deallocate(allBuffSend(iSend)%buffer, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate allBuffSend%buffer fails with message "//trim(message))
       end if
    end do
    
    if (nSend > 0) then
       deallocate(allBuffSend, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate allBuffSend fails with message "//trim(message))
       end if
    end if
  end subroutine SendRecvConvertDomainDecomp_3D

  


  
  subroutine TestSendRecvConvertDomainDecomp(testDim, &
       oneNodeDimensionsFrom, oneNodeDimensionsTo, oneConvertDomainDecomp)

    integer, intent(in) :: testDim ! 2 or 3 to test 2D or 3D array
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensionsFrom
    ! pointer intent(in), values intent(in)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensionsTo
    ! pointer intent(in), values intent(in)
    type(ConvertDomainDecomp), pointer, intent(in) :: oneConvertDomainDecomp
    ! pointer intent(in), values intent(in)

    integer :: ierr
    integer :: nz
    integer :: sizeXFrom
    integer :: sizeYFrom
    integer :: xoFrom
    integer :: yoFrom
    integer :: sizeXTo
    integer :: sizeYTo
    integer :: xoTo
    integer :: yoTo
    integer :: x
    integer :: y
    integer :: z
    real, pointer, contiguous :: fieldFrom_2D(:,:)
    real, pointer, contiguous :: fieldTo_2D(:,:)
    real, pointer, contiguous :: fieldFrom_3D(:,:,:)
    real, pointer, contiguous :: fieldTo_3D(:,:,:)
    character(len=512) :: message

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(TestSendRecvConvertDomainDecomp)**"
    logical, parameter :: dumpLocal=.false.

    if ((testDim /= 2) .and. (testDim /= 3)) then
       call fatal_error(h//" requires testDim==2 ou 3")
    end if

    if (.not. associated(oneConvertDomainDecomp)) then
       if (dumpLocal) then
          call MsgDump(h//" message passing correction is not tested since oneConvertDomainDecomp is not associated")
       end if
       return
    end if

    if (dumpLocal) then
       call MsgDump(h//" test conversion from "//&
            trim(oneConvertDomainDecomp%FromName)//" to "//&
            trim(oneConvertDomainDecomp%ToName))
    end if

    ! allocate and initialize field to be sent

    nz = oneNodeDimensionsFrom%mzp
    sizeXFrom = oneNodeDimensionsFrom%mxp
    sizeYFrom = oneNodeDimensionsFrom%myp
    xoFrom = oneNodeDimensionsFrom%i0 ! global x index offset
    yoFrom = oneNodeDimensionsFrom%j0 ! global y index offset

    if (dumpLocal) then
       write(str(1),"(i8)") xoFrom
       write(str(2),"(i8)") yoFrom
       call MsgDump(h//" From: (xo,yo)=("//&
            trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//")")
    end if

    if (testDim == 2) then
       allocate(fieldFrom_2D(sizeXFrom,sizeYFrom), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate fieldFrom_2D fails with message "//trim(message))
       end if
       
       do y = 1, sizeYFrom
          do x = 1, sizeXFrom
             fieldFrom_2D(x,y) = real(100*(x+xoFrom) + (y+yoFrom))
          end do
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") int(fieldFrom_2d(1,1))
          write(str(2),"(i8)") sizeXFrom
          write(str(3),"(i8)") sizeYFrom
          write(str(4),"(i8)") int(fieldFrom_2d(sizeXFrom,sizeYFrom))
          call MsgDump(h//&
               " fieldFrom(1,1)="//trim(adjustl(str(1)))//&
               " fieldFrom("//trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//")="//&
               trim(adjustl(str(4))))
       end if
               
    else
       allocate(fieldFrom_3D(nz,sizeXFrom,sizeYFrom), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate fieldFrom_3D fails with message "//trim(message))
       end if
       
       do y = 1, sizeYFrom
          do x = 1, sizeXFrom
             do z = 1, nz
                fieldFrom_3D(z,x,y) = real(10000*z + 100*(x+xoFrom) + (y+yoFrom))
             end do
          end do
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") int(fieldFrom_3d(1,1,1))
          write(str(2),"(i8)") sizeXFrom
          write(str(3),"(i8)") sizeYFrom
          write(str(4),"(i8)") int(fieldFrom_3d(nz,sizeXFrom,sizeYFrom))
          write(str(5),"(i8)") nz
          call MsgDump(h//&
               " fieldFrom(1,1,1)="//trim(adjustl(str(2)))//&
               " fieldFrom("//&
               trim(adjustl(str(5)))//","//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//")="//&
               trim(adjustl(str(4))))
       end if
               
    end if
    
    ! allocate and initialize field to be received

    sizeXTo = oneNodeDimensionsTo%mxp
    sizeYTo = oneNodeDimensionsTo%myp
    xoTo = oneNodeDimensionsTo%i0 ! global x index offset
    yoTo = oneNodeDimensionsTo%j0 ! global y index offset

    if (testDim == 2) then
       allocate(fieldTo_2D(sizeXTo,sizeYTo), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate fieldTo_2D fails with message "//trim(message))
       end if
       fieldTo_2D=0.0
    else
       allocate(fieldTo_3D(nz,sizeXTo,sizeYTo), stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" allocate fieldTo_3D fails with message "//trim(message))
       end if
       fieldTo_3D=0.0
    end if

    ! message passing to test
    
    if (testDim == 2) then
       call SendRecvConvertDomainDecomp(fieldFrom_2D, fieldTo_2D, oneConvertDomainDecomp)
    else
       call SendRecvConvertDomainDecomp(fieldFrom_3D, fieldTo_3D, oneConvertDomainDecomp)
    end if

    ! verify received field correctness

    if (dumpLocal) then
       write(str(1),"(i8)") xoTo
       write(str(2),"(i8)") yoTo
       call MsgDump(h//" To: (xo,yo)=("//&
            trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//")")
    end if

    if (testDim == 2) then
       do y = 1, sizeYTo
          do x = 1, sizeXTo
             if (int(fieldTo_2D(x,y)) /= 100*(x+xoTo) + (y+yoTo)) then
                write(str(1),"(i8)") x
                write(str(2),"(i8)") y
                write(str(3),"(i8)") 100*(x+xoTo) + (y+yoTo)
                write(str(4),"(i8)") int(fieldTo_2D(x,y))
                call fatal_error(h//" for local indices ("//&
                     trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//")"//&
                     " fieldTo_2D is "//trim(adjustl(str(4)))//&
                     " but should be "//trim(adjustl(str(3))))
             end if
          end do
       end do

       deallocate(fieldFrom_2D, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate fieldFrom_2D fails with message "//trim(message))
       end if
       
       deallocate(fieldTo_2D, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate fieldTo_2D fails with message "//trim(message))
       end if
       
    else
       do y = 1, sizeYTo
          do x = 1, sizeXTo
             do z = 1, nz
                if (int(fieldTo_3D(z,x,y)) /= 10000*z + 100*(x+xoTo) + (y+yoTo)) then
                   write(str(1),"(i8)") x
                   write(str(2),"(i8)") y
                   write(str(3),"(i8)") 10000*z + 100*(x+xoTo) + (y+yoTo)
                   write(str(4),"(i8)") int(fieldTo_3D(z,x,y))
                   write(str(5),"(i8)") z
                   call fatal_error(h//" for local indices ("//&
                        trim(adjustl(str(5)))//","//trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//")"//&
                        " fieldTo_3D is "//trim(adjustl(str(4)))//&
                        " but should be "//trim(adjustl(str(3))))
                end if
             end do
          end do
       end do

       deallocate(fieldFrom_3D, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate fieldFrom_3D fails with message "//trim(message))
       end if
       
       deallocate(fieldTo_3D, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate fieldTo_3D fails with message "//trim(message))
       end if
    end if

    if (dumpLocal) then
       call MsgDump(h//" finished correctly")
    end if
  end subroutine TestSendRecvConvertDomainDecomp
end module ModConvertDomainDecomp
