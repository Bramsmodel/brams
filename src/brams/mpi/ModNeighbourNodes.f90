module ModNeighbourNodes

  ! stores info of all BRAMS nodes that are neighbours to this node.
  !
  ! There are two kinds of neighbours to this node:
  !
  ! 1. another node which grid domain 
  ! (own grid points and ghost zone grid points) intersects this node
  ! own grid points (no ghost zone)
  !
  ! 2.another node which own grid points
  ! intersects this node grid domain (own grid points and ghost zone
  ! points)
  !
  ! For neighbours of the first kind, this node sould send its own
  ! points that reside on the ghost zone of the other node; the other
  ! node should receive the update.
  !
  ! For neighbours of the second kind, the other node sould send its own
  ! points that reside on the ghost zone of this node; this node
  ! should receive the update.


  use ModParallelEnvironment, only: ParallelEnvironment
  use ModParallelEnvironment, only: Brams2MpiProcNbr
  use ModParallelEnvironment, only: MsgDump

  use ModGridDims, only: &
       GridDims

  use ModDomainDecomp, only: &
       DomainDecomp

  implicit none
  private
  public :: NeighbourNodes
  public :: CreateNeighbourNodes
  public :: DumpNeighbourNodes
  public :: DestroyNeighbourNodes
  public :: NodesToSendRecvMessages
  public :: IncludeDomainBoundaries
  public :: GetNumberOfNeighbours

  type NeighbourNodes
     integer :: nNeigh
     ! how many neighbours
     integer, allocatable :: neigh(:)
     ! array of neighbours BRAMS process number
  end type NeighbourNodes

contains






  subroutine Inter(&
       xs1, xe1, ys1, ye1, &
       xs2, xe2, ys2, ye2, &
       xsInter, xeInter, ysInter, yeInter, &
       hasInter)

    ! intersection of areas [xs1:xe1,ys1:ye1] and [xs2:xe2,ys2:ye2].
    ! Returns in hasInter if there is (or not) intersection and the
    ! intersection itself at [xsInter:xeInter,ysInter:yeInter].
    
    integer, intent(in) :: xs1
    integer, intent(in) :: xe1
    integer, intent(in) :: ys1
    integer, intent(in) :: ye1
    integer, intent(in) :: xs2
    integer, intent(in) :: xe2
    integer, intent(in) :: ys2
    integer, intent(in) :: ye2
    integer, intent(out) :: xsInter
    integer, intent(out) :: xeInter
    integer, intent(out) :: ysInter
    integer, intent(out) :: yeInter
    logical, intent(out) :: hasInter

    character(len=8) :: str(10)
    character(len=256) :: strOut
    character(len=*), parameter :: h="**(Inter)**"
    logical, parameter :: dumpLocal=.false.
    
    xsInter = max(xs1,xs2)
    xeInter = min(xe1,xe2)
    ysInter = max(ys1,ys2)
    yeInter = min(ye1,ye2)
    hasInter = &
         xsInter <= xeInter .and. &
         ysInter <= yeInter
    if (dumpLocal) then
       write(str(1),"(i8)") xs1
       write(str(2),"(i8)") xe1
       write(str(3),"(i8)") ys1
       write(str(4),"(i8)") ye1
       write(str(5),"(i8)") xs2
       write(str(6),"(i8)") xe2
       write(str(7),"(i8)") ys2
       write(str(8),"(i8)") ye2
       strOut=&
            " ["//trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//"] inter"//&
            " ["//trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//","//&
            trim(adjustl(str(7)))//":"//trim(adjustl(str(8)))//"] ="
       if (hasInter) then
          write(str(1),"(i8)") xsInter
          write(str(2),"(i8)") xeInter
          write(str(3),"(i8)") ysInter
          write(str(4),"(i8)") yeInter
          strOut=trim(strOut)//&
               " ["//trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//"]"
       else
          strOut=trim(strOut)//" empty"
       end if
       call MsgDump(h//trim(strOut))
    end if
  end subroutine Inter






  function CreateNeighbourNodes(ParEnv, GlobalOwn, &
       GlobalWithGhost, varName)

    ! Finds all neighbours of this node
    ! that own grid cells to update this node ghost zone
    ! united with all neighbours of this node
    ! which ghost zone cells are owned by this node.
    ! That is, all neighbour nodes that exchange messages
    ! with this node for ghost zone updates.
    ! Returns a null pointer if no ghost zone updates.

    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    character(len=*), intent(in) :: varName
    type(NeighbourNodes), pointer :: CreateNeighbourNodes

    integer :: myNum
    integer :: nmachs
    integer :: node
    integer :: nNeigh
    integer :: cnt
    integer :: xsInter, xeInter, ysInter, yeInter
    integer :: ierr
    logical :: myNumSend, myNumRecv
    logical, allocatable :: isNeighbour(:)
    logical, parameter :: dumpLocal=.false.
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(CreateNeighbourNodes)**"

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" invoked with null GlobalWithGhost ")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" invoked with null GlobalOwn ")
    end if

    nmachs = ParEnv%nmachs
    myNum = ParEnv%myNum

    ! scratch grid partition storing which nodes are neighbour of this node:
    ! these are nodes where extended partition (with ghost zone) intersect
    ! this node partition, and as so, this node has to
    ! send owned cells and the other node has to receive these cells to
    ! update ghost zone,
    ! or 
    ! nodes where the partition intersects with this node
    ! extended partition (with ghost zone), ans as so, this node has to
    ! receive ghost zone cells owned by the other nodes, and the other nodes
    ! has to send owned dells

    allocate(isNeighbour(nmachs), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate isNeighbour fails with stat="//&
            trim(adjustl(c0)))
    end if
    
    do node = 1, nmachs

       ! this node send owned cells to update ghost zones of other nodes
       call Inter(&
            GlobalOwn%xb(myNum), GlobalOwn%xe(myNum), &
            GlobalOwn%yb(myNum), GlobalOwn%ye(myNum), &
            GlobalWithGhost%xb(node), GlobalWithGhost%xe(node), & 
            GlobalWithGhost%yb(node), GlobalWithGhost%ye(node), &
            xsInter, xeInter, ysInter, yeInter, myNumSend)

       ! this node receives cells owned by other nodes to update
       ! this node ghost zone
       call Inter(&
            GlobalOwn%xb(node), GlobalOwn%xe(node), &
            GlobalOwn%yb(node), GlobalOwn%ye(node), &
            GlobalWithGhost%xb(myNum), GlobalWithGhost%xe(myNum), & 
            GlobalWithGhost%yb(myNum), GlobalWithGhost%ye(myNum), &
            xsInter, xeInter, ysInter, yeInter, myNumRecv)
       isNeighbour(node) = myNumSend .or. myNumRecv
    end do

    ! exclude own node

    isNeighbour(myNum) = .false.

    ! how many neighbours; default is no neighbour

    nNeigh = count(isNeighbour)
    CreateNeighbourNodes => null()
    if (nNeigh /= 0) then

       ! build CreateNeighbourNodes, by storing
       ! the number and the list of
       ! neighbour nodes

       allocate(CreateNeighbourNodes, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" allocate CreateNeighbourNodes fails with stat="//&
               trim(adjustl(c0)))
       end if
       CreateNeighbourNodes%nNeigh = nNeigh
       allocate(CreateNeighbourNodes%neigh(nNeigh), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" allocate CreateNeighbourNodes%neigh fails with stat="//&
               trim(adjustl(c0)))
       end if

       cnt = 0
       do node =  1, nmachs
          if (isNeighbour(node)) then
             cnt = cnt + 1
             CreateNeighbourNodes%neigh(cnt) = node
          end if
       end do
       if (cnt /= CreateNeighbourNodes%nNeigh) then
          write(c0,"(i8)") cnt
          write(c1,"(i8)") CreateNeighbourNodes%nNeigh
          call fatal_error(h//" inconsistence: cnt ("//trim(adjustl(c0))//&
               ") differs from nNeigh ("//trim(adjustl(c1))//")")
       end if
    end if
    deallocate(isNeighbour, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate isNeighbour fails with stat="//&
            trim(adjustl(c0)))
    end if

    if (dumpLocal) then
       call DumpNeighbourNodes(CreateNeighbourNodes, varName)
    end if
  end function CreateNeighbourNodes



  ! DumpNeighbourNodes: at this process dump file



  subroutine DumpNeighbourNodes(OneNeighbourNodes, varName)
    type(NeighbourNodes), pointer :: OneNeighbourNodes
    character(len=*), intent(in) :: varName
    
    integer :: neigh
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DumpNeighbourNodes)**"

    if (.not. associated(OneNeighbourNodes)) then
       call MsgDump(h//" variable "//trim(varName)//" is not associated")
    else
       write(c0,"(i8)") OneNeighbourNodes%nNeigh
       call MsgDump(h//" for variable "//trim(varName)//&
            " this node has "//trim(adjustl(c0))//&
            " neighbour MPI ranks:", .true.)
       do neigh = 1, OneNeighbourNodes%nNeigh-1
          write(c0,"(i8)") Brams2MpiProcNbr(OneNeighbourNodes%neigh(neigh))
          call MsgDump(" "//trim(adjustl(c0))//",",.true.)
       end do
       write(c0,"(i8)") Brams2MpiProcNbr(OneNeighbourNodes%neigh(neigh))
       call MsgDump(" "//trim(adjustl(c0)))
    end if
  end subroutine DumpNeighbourNodes



  ! DestroyNeighbourNodes: returns allocated area, if any



  subroutine DestroyNeighbourNodes(OneNeighbourNodes)
    type(NeighbourNodes), pointer, intent(inout) :: OneNeighbourNodes

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DestroyNeighbourNodes)**"
    
    if (associated(OneNeighbourNodes)) then
       deallocate(OneNeighbourNodes%neigh, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate OneNeighbourNodes%neigh fails with stat="//&
               trim(adjustl(c0)))
       end if
       deallocate(OneNeighbourNodes, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate OneNeighbourNodes fails with stat="//&
               trim(adjustl(c0)))
       end if
    end if
    nullify(OneNeighbourNodes)
  end subroutine DestroyNeighbourNodes






  subroutine NodesToSendRecvMessages(thisNode, Neigh, GlobalOwn, &
       xbToUpdate, xeToUpdate, ybToUpdate, yeToUpdate, &
       xbSend, xeSend, ybSend, yeSend, willSend, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
       varName)

    ! given a rectangular region of grid points to be updated
    ! at each rank and the rectangular region of grid points
    ! owned by each rank,
    ! find to which rank this node will send messages
    ! and find from which other rank this rank will receive messages
    ! to update desired points;
    ! also finds the rectangular region of each rank to be
    ! sent or received.
    
    integer, intent(in) :: thisNode
    ! this rank BRAMS process number
    type(NeighbourNodes), pointer, intent(in) :: Neigh
    ! ranks that are neighbour to this rank
    type(DomainDecomp), pointer, intent(in)  :: GlobalOwn
    ! rectangular grid point regions owned by each rank
    integer, intent(in) :: xbToUpdate(:)
    integer, intent(in) :: xeToUpdate(:)
    integer, intent(in) :: ybToUpdate(:)
    integer, intent(in) :: yeToUpdate(:)
    ! rectangular grid point region to be updated at each rank;
    ! contains global indices; arrays are indexed by
    ! BRAMS process number
    integer, intent(out) :: xbSend(:)
    integer, intent(out) :: xeSend(:)
    integer, intent(out) :: ybSend(:)
    integer, intent(out) :: yeSend(:)
    ! rectangular grid point region to be sent to each neighbour;
    ! contains global indices; arrays indexed by
    ! neighbour number index
    logical, intent(out) :: willSend(:)
    ! to which BRAMS process number this rank
    ! will send msgs; array indexed by neighbour number index
    integer, intent(out) :: xbRecv(:)
    integer, intent(out) :: xeRecv(:)
    integer, intent(out) :: ybRecv(:)
    integer, intent(out) :: yeRecv(:)
    ! rectangular grid point region to be received from each neighbour;
    ! contains global indices; arrays indexed by
    ! neighbour number index
    logical, intent(out) :: willRecv(:)
    ! from which BRAMS process number this rank
    ! will receive msgs; array indexed by neighbour number index
    character(len=*), intent(in) :: varName
    ! message set variable name, only for dumpLocal

    integer :: otherNode
    integer :: nNeigh
    integer :: i, indMsg
    logical, parameter :: dumpLocal=.false.
    character(len=8) :: c0, c1, c2, c3, c4
    character(len=*), parameter :: h="**(NodesToSendRecvMessages)**"

    ! check arguments

    if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" starts with null GlobalOwn")
    end if

    ! default output is no message to send or receive

    xbSend=0
    xeSend=0
    ybSend=0
    yeSend=0
    willSend=.false.

    xbRecv=0
    xeRecv=0
    ybRecv=0
    yeRecv=0
    willRecv=.false.

    ! no messages if no neighbours

    if (.not. associated(Neigh)) then
       if (dumpLocal) then
          call MsgDump(h//" no neighbour for "//trim(adjustl(varName)))
       end if
       return
    end if

    ! auxiliar variables

    nNeigh = Neigh%nNeigh

    ! this node will send messages to other node
    ! if the region of the other node to be updated
    ! intersects with the region owned by this node.
    ! find which neighbours this node will send messages to
    ! and the global indices that this node will send

    do i = 1, nNeigh
       otherNode = Neigh%neigh(i)
       if (dumpLocal) then
          write(c0,"(i8)") Brams2MpiProcNbr(otherNode)
          write(c1,"(i8)") Brams2MpiProcNbr(thisNode)
          call MsgDump(h//" intersection region of GlobalOwn of thisNode "//&
               "(MPIrank="//trim(adjustl(c1))//") with region to update "//&
               " of otherNode (MPIrank="//trim(adjustl(c0))//")")
       end if
       call Inter(&
            xbToUpdate(otherNode), &
            xeToUpdate(otherNode), &
            ybToUpdate(otherNode), &
            yeToUpdate(otherNode), &
            GlobalOwn%xb(thisNode), &
            GlobalOwn%xe(thisNode), &
            GlobalOwn%yb(thisNode), &
            GlobalOwn%ye(thisNode), &
            xbSend(i), xeSend(i), &
            ybSend(i), yeSend(i), &
            willSend(i))

       if (dumpLocal) then
          write(c0,"(i8)") Brams2MpiProcNbr(otherNode)
          if (willSend(i)) then
             write(c1,"(i8)") xbSend(i)
             write(c2,"(i8)") xeSend(i)
             write(c3,"(i8)") ybSend(i)
             write(c4,"(i8)") yeSend(i)
             call MsgDump(h//" for "//trim(adjustl(varName))//&
                  " this node will send to MPI process "//trim(adjustl(c0))//&
                  " the region ["//&
                  trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
                  trim(adjustl(c3))//":"//trim(adjustl(c4))//"]")
          else
             call MsgDump(h//" no send message for MPI proc "//trim(adjustl(c0)))
          end if
       end if
    end do

    ! this node will receive messages from other nodes
    ! if the region of this node to be updated
    ! intersects with the region owned by the other node.
    ! find from which neighbours this node will receive messages
    ! and the global indices to be received from each neighbour

    do i = 1, nNeigh
       otherNode = Neigh%neigh(i)
       call Inter(&
            xbToUpdate(thisNode), &
            xeToUpdate(thisNode), &
            ybToUpdate(thisNode), &
            yeToUpdate(thisNode), &
            GlobalOwn%xb(otherNode), &
            GlobalOwn%xe(otherNode), &
            GlobalOwn%yb(otherNode), &
            GlobalOwn%ye(otherNode), &
            xbRecv(i), xeRecv(i), &
            ybRecv(i), yeRecv(i), &
            willRecv(i))

            
       if (dumpLocal) then
          write(c0,"(i8)") Brams2MpiProcNbr(otherNode)
          if (willRecv(i) .and. dumpLocal) then
             write(c1,"(i8)") xbRecv(i)
             write(c2,"(i8)") xeRecv(i)
             write(c3,"(i8)") ybRecv(i)
             write(c4,"(i8)") yeRecv(i)
             call MsgDump(h//" for "//trim(adjustl(varName))//&
                  " this node will recv from MPI process "//trim(adjustl(c0))//&
                  " the region ["//&
                  trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
                  trim(adjustl(c3))//":"//trim(adjustl(c4))//"]")
          else
             call MsgDump(h//" no recv message from MPI proc "//trim(adjustl(c0)))
          end if
       end if
    end do
  end subroutine NodesToSendRecvMessages






  subroutine IncludeDomainBoundaries(Neigh, GridSize, GlobalOwn, &
       xbComm, xeComm, ybComm, yeComm, willComm, &
       varName)

    type(NeighbourNodes), pointer, intent(in) :: Neigh
    type(GridDims), pointer, intent(in) :: GridSize
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn

    integer, intent(inout) :: xbComm(:)          ! global index, dimensioned Neigh%nNeigh
    integer, intent(inout) :: xeComm(:)          ! global index, dimensioned Neigh%nNeigh
    integer, intent(inout) :: ybComm(:)          ! global index, dimensioned Neigh%nNeigh
    integer, intent(inout) :: yeComm(:)          ! global index, dimensioned Neigh%nNeigh
    logical, intent(inout) :: willComm(:)        ! global index, dimensioned Neigh%nNeigh

    ! communication exchange name (variable of type message set)

    character(len=*), intent(in) :: varName

    integer :: indNeigh
    integer :: proc
    logical :: extend
    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(IncludeDomainBoundaries)**"

    ! check arguments

    if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" starts with null GlobalOwn")
    end if

    if (.not. associated(Neigh)) then
       return
    end if

    do indNeigh = 1, Neigh%nNeigh
       if (willComm(indNeigh)) then
          extend=.false.
          proc = Neigh%neigh(indNeigh)
          if (xbComm(indNeigh) == 2) then
             xbComm(indNeigh) = 1
             extend=.true.
          end if
          if (xeComm(indNeigh) == GridSize%nnxp-1) then
             xeComm(indNeigh) = GridSize%nnxp
             extend=.true.
          end if
          if (ybComm(indNeigh) == 2) then
             ybComm(indNeigh) = 1
             extend=.true.
          end if
          if (yeComm(indNeigh) == GridSize%nnyp-1) then
             yeComm(indNeigh) = GridSize%nnyp
             extend=.true.
          end if
          if (extend .and. dumpLocal) then
             write(str(1),"(i8)") proc
             write(str(2),"(i8)") xbComm(indNeigh)
             write(str(3),"(i8)") xeComm(indNeigh)
             write(str(4),"(i8)") ybComm(indNeigh)
             write(str(5),"(i8)") yeComm(indNeigh)
             call MsgDump(h//" communication region for proc "//trim(adjustl(str(1)))//&
                  " extended to ["//&
                  trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
                  trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//"] "//&
                  " for communication "//trim(adjustl(varName)))
             end if
       end if
    end do
  end subroutine IncludeDomainBoundaries






  integer function GetNumberOfNeighbours(OneNeighbourNodes)
    type(NeighbourNodes), pointer, intent(in) :: OneNeighbourNodes
    integer :: neigh
    character(len=8) :: c0
    character(len=*), parameter :: h="**(GetNumberOfNeighbours)**"

    if (.not. associated(OneNeighbourNodes)) then
       call MsgDump(h//" empty")
    else
       GetNumberOfNeighbours = OneNeighbourNodes%nNeigh
    end if
  end function GetNumberOfNeighbours
end module ModNeighbourNodes

