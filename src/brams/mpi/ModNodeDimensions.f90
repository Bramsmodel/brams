module ModNodeDimensions
  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModGridDims, only: &
       GridDims

  use ModDomainDecomp, only: &
       DomainDecomp
  
  implicit none
  private
  public :: NodeDimensions
  public :: CreateNodeDimensions
  public :: DestroyNodeDimensions
  public :: DumpNodeDimensions

  ! NodeDimensions: stores indices of this process domain decomposed sub-domain

  type NodeDimensions

     ! GhostZoneWidth: Width of Ghost Zone (in grid cells)

     integer :: GhostZoneWidth

     ! nxp, nyp, nzp are sizes (in grid points) of this process sub-domain,
     ! including ghost zone;
     ! fields with this GhostZoneWidth should be dimensioned (nzp,nxp,nyp)

     integer :: mxp
     !#x points (with ghost zone) on this process sub-domain
     integer :: myp
     !#y points (with ghost zone) on this process sub-domain
     integer :: mzp
     !#z points (with ghost zone) on this process sub-domain

     ! bounds of the own sub-domain part, that is, lower and upper
     ! index at x and y dimensions in this process sub-domain
     ! that are exclusively owned by this process, excluding ghost zones.
     ! These variables should be the lower and upper bounds of 
     ! loops intended to run through all sub-domain points

     integer :: ia
     ! first x index of the sub-domain exclusively owned by this process
     integer :: iz
     ! last x index of the sub-domain exclusively owned by this process
     integer :: izu
     ! last x index of the sub-domainexclusively owned by this process
     ! for wind x-component (u) - due to staggered grid

     integer :: ja
     ! first y index of the sub-domain exclusively owned by this process
     integer :: jz
     ! last y index of the sub-domain exclusively owned by this process
     integer :: jzv
     ! last y index of the sub-domainexclusively owned by this process
     ! for wind y-component (v) - due to staggered grid

     ! local index <--> global index conversion:
     ! global index = local index + offset

     integer :: i0
     ! x offset
     integer :: j0
     ! y offset

     ! wether this process sub-domain boundary is also a full grid boundary

     logical :: boundNorth
     ! on y axis upper limit
     logical :: boundSouth
     ! on y axis lower limit
     logical :: boundEast
     ! on x axis upper limit
     logical :: boundWest
     ! on x axis lower limit
     
  end type NodeDimensions

contains





  function CreateNodeDimensions(GridSize, ParEnv, LocalOwn, GlobalOwn, &
       verticalGhostZoneWidth, surfaceGhostZoneWidth, varName) result(res)

    ! Creates a pointer to a variable of type NodeDimensions named
    ! varName for given grid and parallel environment. Performs domain
    ! decomposition, filling all components of the created variable

    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(DomainDecomp), pointer, intent(in) :: LocalOwn
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    integer, intent(in) :: verticalGhostZoneWidth
    integer, intent(in) :: surfaceGhostZoneWidth
    character(len=*), intent(in) :: varName
    type(NodeDimensions), pointer :: res

    integer :: myNum
    integer :: ierr

    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateNodeDimensions)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(GridSize)) then
       call fatal_error(h//" invoked with null GridSize")
    else if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    else if (.not. associated(LocalOwn)) then
       call fatal_error(h//" invoked with null LocalOwn")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" invoked with null GlobalOwn")
    end if

    myNum = ParEnv%myNum

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate CreateNodeDimensions fails with stat="//&
            trim(adjustl(c0)))
    end if

    res%GhostZoneWidth = surfaceGhostZoneWidth

    res%boundNorth = btest(LocalOwn%ibcon(myNum),4)
    res%boundSouth = btest(LocalOwn%ibcon(myNum),3)
    res%boundEast = btest(LocalOwn%ibcon(myNum),2)
    res%boundWest = btest(LocalOwn%ibcon(myNum),1)
    
    res%mxp = LocalOwn%xe(myNum)-LocalOwn%xb(myNum)+1+2*surfaceGhostZoneWidth
    res%myp = LocalOwn%ye(myNum)-LocalOwn%yb(myNum)+1+2*surfaceGhostZoneWidth
    res%mzp = GridSize%nnzp + 2*verticalGhostZoneWidth
    
    res%ia=LocalOwn%xb(myNum)
    res%iz=LocalOwn%xe(myNum)
    if (res%boundEast) then
       res%izu=res%iz-1
    else
       res%izu=res%iz
    end if

    res%ja=LocalOwn%yb(myNum)
    res%jz=LocalOwn%ye(myNum)
    if (res%boundNorth) then
       res%jzv=res%jz-1
    else
       res%jzv=res%jz
    end if

    res%i0=GlobalOwn%xb(myNum) - surfaceGhostZoneWidth - 1
    res%j0=GlobalOwn%yb(myNum) - surfaceGhostZoneWidth - 1
     
    if (dumpLocal) then
       call DumpNodeDimensions(res, varName)
    end if
  end function CreateNodeDimensions





  subroutine DestroyNodeDimensions(oneNodeDimensions)

    ! Removes storage allocated to a variable of type NodeDimensions

    type(NodeDimensions), pointer, intent(inout) :: oneNodeDimensions

    integer :: ierr
    character(len=*), parameter :: h="**(DestroyNodeDimensions)**"
    character(len=8) :: c0

    if (associated(oneNodeDimensions)) then
       deallocate(oneNodeDimensions, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate oneNodeDimensions fails with stat="//&
               trim(adjustl(c0)))
       end if
    end if
    nullify(oneNodeDimensions)
  end subroutine DestroyNodeDimensions





  subroutine DumpNodeDimensions(oneNodeDimensions, varName)

    ! Dumps info stored at a variable of type NodeDimensions
    
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    character(len=*), intent(in) :: varName


    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpNodeDimensions)**"

    ! dumps at selected unit

    if (.not. associated(oneNodeDimensions)) then
       call MsgDump (h//" empty NodeDimensions named "//trim(varName))
       return
    else
       write(str(1),"(i8)") oneNodeDimensions%GhostZoneWidth
       call MsgDump (h//" of variable "//trim(adjustl(varName))//&
            " with ghost zone of width "//trim(adjustl(str(1))))
    end if

    write(str(1),"(i8)") oneNodeDimensions%mzp
    write(str(2),"(i8)") oneNodeDimensions%mxp
    write(str(3),"(i8)") oneNodeDimensions%myp
    call MsgDump(h//" Sub-domain dimensions"//&
         " mzp="//trim(adjustl(str(1)))//&
         ", mxp="//trim(adjustl(str(2)))//&
         ", myp="//trim(adjustl(str(3))))

    write(str(1),"(i8)") oneNodeDimensions%ia
    write(str(2),"(i8)") oneNodeDimensions%iz
    write(str(3),"(i8)") oneNodeDimensions%izu
    call MsgDump(h//" x axis indexing"//&
         " ia="//trim(adjustl(str(1)))//&
         ", iz="//trim(adjustl(str(2)))//&
         ", izu="//trim(adjustl(str(3))))

    write(str(1),"(i8)") oneNodeDimensions%ja
    write(str(2),"(i8)") oneNodeDimensions%jz
    write(str(3),"(i8)") oneNodeDimensions%jzv
    call MsgDump(h//" y axis indexing"//&
         " ja="//trim(adjustl(str(1)))//&
         ", jz="//trim(adjustl(str(2)))//&
         ", jzv="//trim(adjustl(str(3))))

    write(str(1),"(i8)") oneNodeDimensions%i0
    call MsgDump(h//" on x indexing, <global index> = <local index> + "//&
         trim(adjustl(str(1))))

    write(str(1),"(i8)") oneNodeDimensions%j0
    call MsgDump(h//" on y indexing, <global index> = <local index> + "//&
         trim(adjustl(str(1))))
    
    str(1)=""
    if (oneNodeDimensions%boundWest) str(1)=trim(str(1))//"X-"
    if (oneNodeDimensions%boundEast) str(1)=trim(str(1))//"X+"
    if (oneNodeDimensions%boundSouth) str(1)=trim(str(1))//"Y-"
    if (oneNodeDimensions%boundNorth) str(1)=trim(str(1))//"Y+"
    if (len(str(1)) == 0) then
       call MsgDump(h//" this sub-domain has none full grid boundaries")
    else
       call MsgDump(h//" this sub-domain at full grid boundaries "//&
            trim(adjustl(str(1))))
    end if
  end subroutine DumpNodeDimensions
end module ModNodeDimensions
