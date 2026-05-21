module ModDomainDecomp
  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump, &
       Brams2MpiProcNbr

  use ModGridDims, only: &
       GridDims

  implicit none


  ! DomainDecomp: stores indices of a domain decomposed grid for all MPI ranks.
  !               Arrays are indexed BRAMS rank (MPI rank+1), to avoid starting from 0.
  !               Refering to BRAMS rank i (MPI rank i-1), the indices of the sub-domain
  !               stored at this rank are [xb(i):xe(i), yb(i),ye(i)] and
  !               ibcon(i) stores if any sub-domain boundary is
  !               a full domain boundary.
  !               GhostZoneWidth is the width of the ghost zone
  !               (0 if no ghost zone).
  !               
  !               Starting at BRAMS 7.0, code uses 3 domain decompositions:
  !               (1) RAMS original, that partitions x-y grid surface into x-y rectangles
  !                   as evenly in the number of grid points as possible;
  !               (2) Monotonic Advection in the x direction, partitioning the grid keeping
  !                   all points of the x direction in a rank;
  !               (3) Monotonic Advection in the y direction, partitioning the grid keeping
  !                   all points of the y direction in a rank;
  !
  !               Domain decompositions (2) and (3) are required to obtain binary reproducibility
  !               and are used only when monotonic advection is selected at RAMSIN.
  !               In such a case, if the number of MPI ranks is greater than the number of points in
  !               x (or y) direction, the exceeding ranks will have empty partition.
  !               Ranks with empty partitions have component "Empty" set to .true.

  type DomainDecomp

     ! Empty is true iff the rank has no partition
     ! Empty is always false for BRAMS domain decomposition, and maybe true
     ! for some ranks in both X and Y Monotonic Advection decompositions

     logical, allocatable :: Empty(:)
     
     ! GhostZoneWidth: Width of Ghost Zone in surface points, 0 if no Ghost Zone

     integer :: GhostZoneWidth

     ! Full Direction: Which type of domain decomposition. It is one of
     ! "B" (for Brams domain decomposition),
     ! "X" (for Monotonic Advection in the X direction),
     ! "Y" (for Monotonic Advection in the Y direction)

     character :: FullDirection
     
     ! all arrays are indexed by BRAMS process number
     
     ! x axis indexed from xb to xe with nx indices;
     ! y axis indexed from yb to ye with ny indices;
     ! global indices at GlobalOwn and GlobalWithGhost;
     ! local indices at LocalOwn;

     ! To convert local indices of LocalOwn to global indices and vice-versa:
     ! global x index = local x index + (GlobalWithGhost%xb-1)
     ! global y index = local y index + (GlobalWithGhost%yb-1)

     ! Empty domains have nx=ny=0 and xb>nnxp, xe<1, yb>nnyp, ye<1
     ! at the ranks where domains are empty
     integer, allocatable :: xb(:)   ! first x
     integer, allocatable :: xe(:)   ! last x 
     integer, allocatable :: nx(:)   ! # points x
     integer, allocatable :: yb(:)   ! first y
     integer, allocatable :: ye(:)   ! last y 
     integer, allocatable :: ny(:)   ! # points y

     ! Array ibcon (boundary condition), also indexed by MPI rank + 1, stores
     ! if sub-domain boundary is a global boundary or not. Info is coded
     ! on a bit structure:
     ! Bit 1 is set iff west sub-domain boundary (low x value) is a global boundary
     ! Bit 2 is set iff east sub-domain boundary (high x value) is a global boundary
     ! Bit 3 is set iff south sub-domain boundary (low y value) is a global boundary
     ! Bit 4 is set iff north sub-domain boundary (high y value) is a global boundary

     ! To check if a sub-domain boundary is a global domain boundary, use:
     !  if btest(ibcon(mach),1) is true, sub-domain west boundary is a global boundary
     !  if btest(ibcon(mach),2) is true, sub-domain east boundary is a global boundary
     !  if btest(ibcon(mach),3) is true, sub-domain south boundary is a global boundary
     !  if btest(ibcon(mach),4) is true, sub-domain north boundary is a global boundary

     ! Empty domains have ibcon=0
     integer, allocatable :: ibcon(:) ! full domain boundary flag

  end type DomainDecomp

  
  ! Public procedures create specific domain decompositions and return pointers
  ! to instances of type DomainDecomposition. These procedures are invoked by
  ! ModGrid.f90, and stored as components of type Grid:
  !
  ! GlobalOwn: Brams X-Y partition of domain [2:nnxp-1,2:nnyp-1]
  !            which is the full domain except boundary conditions.
  !            Boundary conditions are marked as so.
  !            Stores global indices, not including ghost zone.
  !            Created by calling CreateGlobalOwn with FullDirection="B"
  !
  ! GlobalOwnMonAdvX and GlobalOwnMonAdvY:
  !            Partition of domain [2:nnxp-1,2:nnyp-1] for monotonic advection,
  !            which is the full domain except boundary conditions.
  !            For MonAdvX, only the Y direction is partitioned across ranks; each
  !            rank keeps the entire X domain. except boundary conditions extremes.
  !            For MonAdvY, only the X direction is partitioned across ranks; each
  !            rank keeps the entire Y domain, except boundary conditions extremes.
  !            Boundary conditions are marked as so.
  !            Stores global indices, not including ghost zone.
  !            If there are more ranks than points in X or Y direction, exceding
  !            ranks have Empty set to true. In such a case, remaining components
  !            are set as described above.
  !            Created by calling CreateGlobalOwn, passing "X" or "Y" as 
  !            argument FullDirection (the direction that will not be partitioned).
  !  
  ! GlobalOwnWithBC: It extends GlobalOwn by including boundary conditions. 
  !            Stores global indices. Created by calling CreateGlobalOwnWithBC, given GlobalOwn.
  !            Procedure CreateGlobalOwnWithBC applies to Brams or Monotonic Advection partitions.
  !  
  ! GlobalWithGhost: It extends GlobalOwn by including ghost zone and boundary conditions. 
  !            Not a real partition, since ghost zone of one sub-domain intersects with 
  !            interior points of other sub-domains. Ghost zone width is given. 
  !            Stores global indices. Created by calling CreateGlobalWithGhost,
  !            given GlobalOwn and GhostZoneWidth. Procedure CreateGlobalWithGhost applies
  !            only to Brams domain decompositions, since Monotonic Advection domain
  !            decomposition has no ghost zone.
  !
  ! LocalOwn: Local indices of interior columns of GlobalWithGhost for Brams domain decomposition
  !           or of GlobalOwnWithBC for Monotonic Advection domain decomposition.
  !           Stores indices used throughout BRAMS (both dynamics and physics);
  !           xb, xe, yb and ye are the local indices converted from GlobalWithGhost (Brams)
  !           or GlobalWithBC (Monotonic Advection) using GlobalOwn info;
  !           nx and ny are total sub-domain size, including boundary conditions and
  !           ghost zone (same as GlobalWithGhost)


  private
  public :: DomainDecomp
  public :: CreateGlobalOwn
  public :: CreateGlobalOwnWithBC
  public :: CreateGlobalWithGhost
  public :: CreateLocalOwn
  public :: DestroyDomainDecomp
  public :: DumpDomainDecomp

  ! scratch areas for dumping
  character(len=8) :: str(10)
  character(len=512) :: message
contains





  function CreateDomainDecomp(nmachs) result(ptr)

    ! Allocates all components of a pointer of 
    ! type DomainDecomp, returning the pointer
    ! Components are not initialized
    
    integer, intent(in) :: nmachs
    type(DomainDecomp), pointer :: ptr

    integer :: ierr
    character(len=*), parameter :: h="**(CreateDomainDecomp)**"

    allocate(ptr, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ptr fails with message="//trim(message))
    end if

    write(str(1),"(i8)") nmachs
    str(2)="("//trim(adjustl(str(1)))//")"

    allocate(ptr%Empty(nmachs), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ptr%Empty"//trim(adjustl(str(2)))//" fails with message="//trim(message))
    end if
    
    allocate(ptr%xb(nmachs), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ptr%xb"//trim(adjustl(str(2)))//" fails with message="//trim(message))
    end if

    allocate(ptr%xe(nmachs), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ptr%xe"//&
            trim(adjustl(str(2)))//" fails with message="//trim(message))
    end if

    allocate(ptr%nx(nmachs), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ptr%nx"//&
            trim(adjustl(str(2)))//" fails with message="//trim(message))
    end if

    allocate(ptr%yb(nmachs), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ptr%yb"//&
            trim(adjustl(str(2)))//" fails with message="//trim(message))
    end if

    allocate(ptr%ye(nmachs), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ptr%ye"//&
            trim(adjustl(str(2)))//" fails with message="//trim(message))
    end if

    allocate(ptr%ny(nmachs), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ptr%ny"//&
            trim(adjustl(str(2)))//" fails with message="//trim(message))
    end if

    allocate(ptr%ibcon(nmachs), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate ptr%ibcon"//&
            trim(adjustl(str(2)))//" fails with message="//trim(message))
    end if
  end function CreateDomainDecomp





  subroutine DestroyDomainDecomp(OneDomainDecomp)

    ! Removes storage allocated to a variable of type DomainDecomp

    type(DomainDecomp), pointer, intent(inout) :: OneDomainDecomp

    integer :: ierr
    character(len=*), parameter :: h="**(DestroyDomainDecomp)**"

    if (associated(oneDomainDecomp)) then

       deallocate(oneDomainDecomp%Empty, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate Empty fails with message="//trim(message))
       end if

       deallocate(oneDomainDecomp%xb, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate xb fails with message="//trim(message))
       end if

       deallocate(oneDomainDecomp%xe, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate xe fails with message="//trim(message))
       end if

       deallocate(oneDomainDecomp%nx, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate nx fails with message="//trim(message))
       end if

       deallocate(oneDomainDecomp%yb, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate yb fails with message="//trim(message))
       end if

       deallocate(oneDomainDecomp%ye, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate ye fails with message="//trim(message))
       end if

       deallocate(oneDomainDecomp%ny, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate ny fails with message="//trim(message))
       end if

       deallocate(oneDomainDecomp%ibcon, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate ibcon fails with message="//trim(message))
       end if

       deallocate(oneDomainDecomp, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneDomainDecomp fails with message="//trim(message))
       end if
    end if
    nullify(oneDomainDecomp)
  end subroutine DestroyDomainDecomp





  subroutine DumpDomainDecomp(oneDomainDecomp, varName)

    ! Dumps info stored at a variable of type DomainDecomp
    
    type(DomainDecomp), pointer, intent(in) :: oneDomainDecomp
    character(len=*), intent(in) :: varName

    integer :: nmachs
    integer :: mach
    character(len=*), parameter :: h="**(DumpDomainDecomp)**"
    character(len=256) :: strOut
    character(len=8) :: bnd

    ! dumps at selected unit

    if (.not. associated(oneDomainDecomp)) then
       call MsgDump (h//" empty DomainDecomp named "//trim(varName))
       return
    else if (.not. allocated(oneDomainDecomp%Empty)) then
       call fatal_error(h//" oneDomainDecomp%Empty not allocated")
    else if (.not. allocated(oneDomainDecomp%xb)) then
       call fatal_error(h//" oneDomainDecomp%xb not allocated")
    else if (.not. allocated(oneDomainDecomp%xe)) then
       call fatal_error(h//" oneDomainDecomp%xe not allocated")
    else if (.not. allocated(oneDomainDecomp%nx)) then
       call fatal_error(h//" oneDomainDecomp%nx not allocated")
    else if (.not. allocated(oneDomainDecomp%yb)) then
       call fatal_error(h//" oneDomainDecomp%yb not allocated")
    else if (.not. allocated(oneDomainDecomp%ye)) then
       call fatal_error(h//" oneDomainDecomp%ye not allocated")
    else if (.not. allocated(oneDomainDecomp%ny)) then
       call fatal_error(h//" oneDomainDecomp%ny not allocated")
    else if (.not. allocated(oneDomainDecomp%ibcon)) then
       call fatal_error(h//" oneDomainDecomp%ibcon not allocated")
    end if

    nmachs = size(oneDomainDecomp%xb)

    strOut=" named "//trim(adjustl(varName))//" is"

    select case (oneDomainDecomp%FullDirection)
    case ("B")
       strOut = trim(strOut)//" Brams"
    case ("X")
       strOut = trim(strOut)//" Monotonic Advection in X"
    case ("Y")
       strOut = trim(strOut)//" Monotonic Advection in Y"
    case default
       call fatal_error(h//" unknown FullDirection **"//trim(oneDomainDecomp%FullDirection)//"**")
    end select

    write(str(1),"(i8)") nmachs
    write(str(2),"(i8)") oneDomainDecomp%GhostZoneWidth
    strOut = trim(strOut)//" domain decomposition with "//&
         trim(adjustl(str(1)))//" sub-domains and ghost zone of width "//&
         trim(adjustl(str(2)))

    call MsgDump (h//trim(strOut))

    call MsgDump(' MPIrank   Empty   x-beg   x-end      nx   y-beg   y-end      ny    cols    ibcon')
    do mach = 1,nmachs
       bnd=""
       if (btest(oneDomainDecomp%ibcon(mach),1)) bnd=trim(bnd)//"X-"
       if (btest(oneDomainDecomp%ibcon(mach),2)) bnd=trim(bnd)//"X+"
       if (btest(oneDomainDecomp%ibcon(mach),3)) bnd=trim(bnd)//"Y-"
       if (btest(oneDomainDecomp%ibcon(mach),4)) bnd=trim(bnd)//"Y+"
       write(message,"(i8,l8,7i8,4x,a8)") &
            Brams2MpiProcNbr(mach), &
            oneDomainDecomp%Empty(mach), &
            oneDomainDecomp%xb(mach), &
            oneDomainDecomp%xe(mach), &
            oneDomainDecomp%nx(mach), &
            oneDomainDecomp%yb(mach), &
            oneDomainDecomp%ye(mach), &
            oneDomainDecomp%ny(mach), &
            oneDomainDecomp%nx(mach)*oneDomainDecomp%ny(mach), &
            adjustl(bnd)
       call MsgDump (trim(message))
    end do
  end subroutine DumpDomainDecomp




  
  function CreateGlobalOwn(GridSize, ParEnv, varName, FullDirection)

    ! Creates a pointer to a variable of type DomainDecomp named
    ! varName for given grid and parallel environment. Performs Brams 
    ! domain decomposition, filling all components of the created variable

    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    character(len=*), intent(in) :: varName
    character(len=*), intent(in) :: FullDirection
    type(DomainDecomp), pointer :: CreateGlobalOwn

    character(len=*), parameter :: h="**(CreateGlobalOwn)**"
    logical, parameter :: dumpLocal=.false.

    integer :: nxp
    integer :: nyp
    integer :: nmachs
    integer :: maxMachs
    integer :: mach
    integer :: quo
    integer :: rem

    if (.not. associated(GridSize)) then
       call fatal_error(h//" invoked with null GridSize")
    else if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    end if

    nxp = GridSize%nnxp
    nyp = GridSize%nnyp
    nmachs = ParEnv%nmachs
    

    if (dumpLocal) then
       write(str(3),"(i8)") nxp
       write(str(1),"(i8)") nyp
       write(str(2),"(i8)") nmachs
       select case (FullDirection)
       case ("X","Y")
          call MsgDump (h//" partitions domain ["//&
               trim(adjustl(str(3)))//" x "//trim(adjustl(str(1)))//&
               "] into "//trim(adjustl(str(2)))//" sub-domains"//&
               " for Monotonic Advection in the "//&
               trim(FullDirection)//" direction")
       case ("B")
          call MsgDump (h//" partitions domain ["//&
               trim(adjustl(str(3)))//" x "//trim(adjustl(str(1)))//&
               "] into "//trim(adjustl(str(2)))//" sub-domains"//&
               " for X-Y regular Brams domain decomposition")
       case default
          call fatal_error(h//" unknown FullDirection **"//trim(FullDirection)//"**")
       end select
    end if

    ! no ghost zone info

    CreateGlobalOwn => CreateDomainDecomp(nmachs)
    CreateGlobalOwn%GhostZoneWidth = 0
    CreateGlobalOwn%FullDirection = FullDirection

    ! select domain decomposition
       
    select case (trim(FullDirection))
    case ("Y")

       ! maximum number of non-empty ranks

       if (nmachs > nxp-2) then
          maxMachs = nxp-2
       else
          maxMachs = nmachs
       end if
       
       ! Monotonic Advection at Y only partitions the X direction
       
       quo = (nxp-2)/maxMachs
       rem = mod(nxp-2,maxMachs)

       ! ranks with non-empty partition
       
       do mach = 1, maxMachs
          CreateGlobalOwn%Empty(mach) = .false.
          CreateGlobalOwn%yb(mach) = 2
          CreateGlobalOwn%ye(mach) = nyp-1
          if (mach <= rem) then
             CreateGlobalOwn%xb(mach) = (mach-1)*(quo+1) + 2
             CreateGlobalOwn%xe(mach) = (mach  )*(quo+1) + 1
          else
             CreateGlobalOwn%xb(mach) = (mach-1)*(quo) + rem + 2
             CreateGlobalOwn%xe(mach) = (mach  )*(quo) + rem + 1
          end if
       end do

       ! ranks with empty partition
       
       do mach = maxMachs + 1, nmachs
          CreateGlobalOwn%Empty(mach) = .true.
          CreateGlobalOwn%yb(mach) = nyp+1
          CreateGlobalOwn%ye(mach) = 0
          CreateGlobalOwn%xb(mach) = nxp+1
          CreateGlobalOwn%xe(mach) = 0
       end do

    case ("X") 

        ! maximum number of non-empty ranks

       if (nmachs > nyp-2) then
          maxMachs = nyp-2
       else
          maxMachs = nmachs
       end if
       
       ! Monotonic Advection at X only partitions the Y direction
       
       quo = (nyp-2)/nmachs
       rem = mod(nyp-2,nmachs)

       ! ranks with non-empty partition
       
       do mach = 1, maxMachs
          CreateGlobalOwn%Empty(mach) = .false.
          CreateGlobalOwn%xb(mach) = 2
          CreateGlobalOwn%xe(mach) = nxp-1
          if (mach <= rem) then
             CreateGlobalOwn%yb(mach) = (mach-1)*(quo+1) + 2
             CreateGlobalOwn%ye(mach) = (mach  )*(quo+1) + 1
          else
             CreateGlobalOwn%yb(mach) = (mach-1)*(quo) + rem + 2
             CreateGlobalOwn%ye(mach) = (mach  )*(quo) + rem + 1
          end if
       end do

       ! ranks with empty partition
       
       do mach = maxMachs + 1, nmachs
          CreateGlobalOwn%Empty(mach) = .true.
          CreateGlobalOwn%yb(mach) = nyp+1
          CreateGlobalOwn%ye(mach) = 0
          CreateGlobalOwn%xb(mach) = nxp+1
          CreateGlobalOwn%xe(mach) = 0
       end do
       
    case ("B") 
       
       ! Brams domain decomposition in the X-Y directions

       call DomainDecompRams(nxp, nyp, nmachs, &
            CreateGlobalOwn%xb, CreateGlobalOwn%xe, &
            CreateGlobalOwn%yb, CreateGlobalOwn%ye)

       CreateGlobalOwn%Empty=.false.
       
    case default
       call fatal_error(h//" unknown FullDirection **"//trim(FullDirection)//"**")
    end select

    ! fill boundary condition ibcon
    
    call MarkBoundary(nxp, nyp, nmachs, CreateGlobalOwn%Empty, &
         CreateGlobalOwn%xb, CreateGlobalOwn%xe, &
         CreateGlobalOwn%yb, CreateGlobalOwn%ye, &
         CreateGlobalOwn%ibcon)
    
    ! fill number of points at each sub-domain 
    
    CreateGlobalOwn%nx = max(0, CreateGlobalOwn%xe - CreateGlobalOwn%xb + 1)
    CreateGlobalOwn%ny = max(0, CreateGlobalOwn%ye - CreateGlobalOwn%yb + 1)
    
    if (dumpLocal) then
       call DumpDomainDecomp(CreateGlobalOwn, varName)
    end if
    
    ! Verify partition correctness
    
    call CheckPartition(nxp, nyp, nmachs, CreateGlobalOwn%empty, &
         CreateGlobalOwn%xb, CreateGlobalOwn%xe, &
         CreateGlobalOwn%yb, CreateGlobalOwn%ye, &
         CreateGlobalOwn%ibcon)
  end function CreateGlobalOwn

  


  
  function CreateGlobalOwnWithBC(GridSize, ParEnv, GlobalOwn, varName)

    ! Creates a pointer to a variable of type DomainDecomp named
    ! varName for given grid and parallel environment. The variable
    ! extends the domain decomposition GlobalOwn by including
    ! boundary conditions. 

    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    character(len=*), intent(in) :: varName
    type(DomainDecomp), pointer :: CreateGlobalOwnWithBC

    character(len=*), parameter :: h="**(CreateGlobalOwnWithBC)**"
    logical, parameter :: dumpLocal=.false.

    integer :: nxp
    integer :: nyp
    integer :: nmachs
    integer :: i

    if (.not. associated(GridSize)) then
       call fatal_error(h//" invoked with null GridSize")
    else if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    end if

    nxp = GridSize%nnxp
    nyp = GridSize%nnyp
    nmachs = ParEnv%nmachs

    if (dumpLocal) then
       write(str(1),"(i8)") nxp
       write(str(2),"(i8)") nyp
       write(str(3),"(i8)") nmachs
       call MsgDump (h//" partitions domain ["//&
            trim(adjustl(str(1)))//" x "//trim(adjustl(str(2)))//&
            "] into "//trim(adjustl(str(3)))//" sub-domains including boundary conditions")
    end if

    CreateGlobalOwnWithBC => CreateDomainDecomp(nmachs)

    CreateGlobalOwnWithBC%GhostZoneWidth = GlobalOwn%GhostZoneWidth
    CreateGlobalOwnWithBC%FullDirection = GlobalOwn%FullDirection
    CreateGlobalOwnWithBC%Empty = GlobalOwn%Empty

    ! include boundary conditions if at border

    do i = 1, nmachs
       if (GlobalOwn%Empty(i)) then
          CreateGlobalOwnWithBC%xb(i) = nxp+1
          CreateGlobalOwnWithBC%xe(i) = 0
          CreateGlobalOwnWithBC%yb(i) = nyp+1
          CreateGlobalOwnWithBC%ye(i) = 0
       else
          if (btest(GlobalOwn%ibcon(i),1)) then
             CreateGlobalOwnWithBC%xb(i) = 1
          else
             CreateGlobalOwnWithBC%xb(i) = GlobalOwn%xb(i)
          end if
          if (btest(GlobalOwn%ibcon(i),2)) then
             CreateGlobalOwnWithBC%xe(i) = nxp
          else
             CreateGlobalOwnWithBC%xe(i) = GlobalOwn%xe(i)
          end if
          if (btest(GlobalOwn%ibcon(i),3)) then
             CreateGlobalOwnWithBC%yb(i) = 1
          else
             CreateGlobalOwnWithBC%yb(i) = GlobalOwn%yb(i)
          end if
          if (btest(GlobalOwn%ibcon(i),4)) then
             CreateGlobalOwnWithBC%ye(i) = nyp
          else
             CreateGlobalOwnWithBC%ye(i) = GlobalOwn%ye(i)
          end if
       end if
    end do

    CreateGlobalOwnWithBC%nx = max(0, &
         CreateGlobalOwnWithBC%xe - &
         CreateGlobalOwnWithBC%xb + 1)

    CreateGlobalOwnWithBC%ny = max(0,&
         CreateGlobalOwnWithBC%ye - &
         CreateGlobalOwnWithBC%yb + 1)

    CreateGlobalOwnWithBC%ibcon = GlobalOwn%ibcon

    ! Verify partition correctness

    call CheckPartition(nxp, nyp, nmachs, CreateGlobalOwnWithBC%empty, &
         CreateGlobalOwnWithBC%xb, CreateGlobalOwnWithBC%xe, &
         CreateGlobalOwnWithBC%yb, CreateGlobalOwnWithBC%ye, &
         CreateGlobalOwnWithBC%ibcon)

    if (dumpLocal) then
       call DumpDomainDecomp(CreateGlobalOwnWithBC, varName)
    end if
  end function CreateGlobalOwnWithBC



  

  function CreateGlobalWithGhost(GridSize, ParEnv, &
       GlobalOwn, GhostZoneWidth, varName)

    ! Creates a pointer to a variable of type DomainDecomp
    ! named varName that extends
    ! the sub-domains of GlobalOwn to include ghost zone
    ! of width GhostZoneWidth. Sub-domain [xb:xe] is limited
    ! by [1:nxp], same with [yb:ye]. Number of points nx, ny
    ! include ghost zone width. Boundary conditions ibcon are copied
    ! from GlobalOwn

    type(GridDims), pointer, intent(in) :: GridSize
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    integer, intent(in) :: GhostZoneWidth
    character(len=*), intent(in) :: varName
    type(DomainDecomp), pointer :: CreateGlobalWithGhost

    integer :: nxp
    integer :: nyp
    integer :: nmachs
    integer :: cell
    integer :: ierr
    character(len=*), parameter :: h="**(CreateGlobalWithGhost)**"
    logical, parameter :: dumpLocal = .false.

    if (.not. associated(GridSize)) then
       call fatal_error(h//" invoked with null GridSize")
    else if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" invoked with null GlobalOwn ")
    else if (GhostZoneWidth < 1) then
       write (str(1),"(i8)") GhostZoneWidth
       call fatal_error(h//" GhostZoneWidth sould be >= 1 but it is "//&
            trim(adjustl(str(1))))
    end if


    nmachs = ParEnv%nmachs
    nxp = GridSize%nnxp
    nyp = GridSize%nnyp

    if (dumpLocal) then
       write(str(1),"(i8)") nxp
       write(str(2),"(i8)") nyp
       write(str(3),"(i8)") GhostZoneWidth
       call MsgDump (h//" "//varName//&
            " inserts ghost zone of width "//trim(adjustl(str(3)))//&
            " into GlobalOwn ["//&
            trim(adjustl(str(1)))//" x "//trim(adjustl(str(2)))//"]")
    end if

    ! allocate return function value
    
    CreateGlobalWithGhost => CreateDomainDecomp(nmachs)

    ! set DomainDecomp components
    
    CreateGlobalWithGhost%GhostZoneWidth = GhostZoneWidth
    CreateGlobalWithGhost%FullDirection = GlobalOwn%FullDirection
    CreateGlobalWithGhost%Empty = GlobalOwn%Empty

    do cell = 1, nmachs
       CreateGlobalWithGhost%ibcon(cell) = GlobalOwn%ibcon(cell)
       if (.not. CreateGlobalWithGhost%Empty(cell)) then
          ! west boundary (lower x axis)
          if (btest(CreateGlobalWithGhost%ibcon(cell),1)) then
             CreateGlobalWithGhost%xb(cell) = 1
          else
             CreateGlobalWithGhost%xb(cell) = max(1,GlobalOwn%xb(cell) - GhostZoneWidth)
          end if
          ! east boundary (higher x axis)
          if (btest(CreateGlobalWithGhost%ibcon(cell),2)) then
             CreateGlobalWithGhost%xe(cell) = nxp
          else
             CreateGlobalWithGhost%xe(cell) = min(nxp,GlobalOwn%xe(cell) + GhostZoneWidth)
          end if
          ! south boundary (lower y axis)
          if (btest(CreateGlobalWithGhost%ibcon(cell),3)) then
             CreateGlobalWithGhost%yb(cell) = 1
          else
             CreateGlobalWithGhost%yb(cell) = max(1,GlobalOwn%yb(cell) - GhostZoneWidth)
          end if
          ! north boundary (higher y axis)
          if (btest(CreateGlobalWithGhost%ibcon(cell),4)) then
             CreateGlobalWithGhost%ye(cell) = nyp
          else
             CreateGlobalWithGhost%ye(cell) = min(nyp,GlobalOwn%ye(cell) + GhostZoneWidth)
          end if
       end if
       CreateGlobalWithGhost%nx(cell) = max(0,CreateGlobalWithGhost%xe(cell)-CreateGlobalWithGhost%xb(cell)+1)
       CreateGlobalWithGhost%ny(cell) = max(0,CreateGlobalWithGhost%ye(cell)-CreateGlobalWithGhost%yb(cell)+1)
    end do

    if (dumpLocal) then
       call DumpDomainDecomp(CreateGlobalWithGhost, varName)
    end if
  end function CreateGlobalWithGhost




  
  function CreateLocalOwn(ParEnv, GlobalWithGhost, &
       GlobalOwn, varName)

    ! Creates a pointer to a variable of type DomainDecomp that defines
    ! the **interior** points of GlobalWithGhost input variable,
    ! that is, the local indices of the sub-domain except ghost zone
    ! and boundary conditions.
    ! Returned pointer stores at [xb:xe,yb:ye] the local indices
    ! of interior points that should be computed at most physics
    ! and dynamics routines. Variables nx, ny should be used to
    ! dimension fields, since they include ghost zones and boundary
    ! conditions.
    
    
    type(ParallelEnvironment), pointer, intent(in) :: ParEnv
    type(DomainDecomp), pointer, intent(in) :: GlobalWithGhost
    type(DomainDecomp), pointer, intent(in) :: GlobalOwn
    character(len=*), intent(in) :: varName
    type(DomainDecomp), pointer :: CreateLocalOwn

    integer :: nmachs
    integer :: cell
    integer :: x0
    integer :: y0
    character(len=*), parameter :: h="**(CreateLocalOwn)**"
    logical, parameter :: dumpLocal = .false.

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" invoked with null GlobalWithGhost ")
    else if (.not. associated(GlobalOwn)) then
       call fatal_error(h//" invoked with null GlobalOwn ")
    end if

    if (dumpLocal) then
       call MsgDump (h//" "//trim(adjustl(varName))//&
            " stores local indices from global indices")
    end if

    nmachs=ParEnv%nmachs

    ! Local Interior info
    
    CreateLocalOwn => CreateDomainDecomp(nmachs)
    CreateLocalOwn%GhostZoneWidth = GlobalWithGhost%GhostZoneWidth 
    CreateLocalOwn%FullDirection = GlobalWithGhost%FullDirection 
    CreateLocalOwn%Empty = GlobalWithGhost%Empty

    do cell = 1, nmachs
       x0 = GlobalWithGhost%xb(cell)-1
       CreateLocalOwn%xb(cell) = GlobalOwn%xb(cell) - x0
       CreateLocalOwn%xe(cell) = GlobalOwn%xe(cell) - x0
       CreateLocalOwn%nx(cell) = GlobalWithGhost%nx(cell)
       y0 = GlobalWithGhost%yb(cell)-1
       CreateLocalOwn%yb(cell) = GlobalOwn%yb(cell) - y0
       CreateLocalOwn%ye(cell) = GlobalOwn%ye(cell) - y0
       CreateLocalOwn%ny(cell) = GlobalWithGhost%ny(cell)
       CreateLocalOwn%ibcon(cell) = GlobalOwn%ibcon(cell)
    end do

    if (dumpLocal) then
       call DumpDomainDecomp(CreateLocalOwn, varName)
    end if
  end function CreateLocalOwn




  
  subroutine DomainDecompRams(nxp,nyp,ranks,&
       xb,xe,yb,ye)

    ! original RAMS domain decomposition

    ! Arguments:
    integer, intent(in)  :: nxp
    integer, intent(in)  :: nyp
    integer, intent(in)  :: ranks
    integer, intent(out) :: xb(ranks)
    integer, intent(out) :: xe(ranks)
    integer, intent(out) :: yb(ranks)
    integer, intent(out) :: ye(ranks)

    ! Local variables:
    real  :: work(nxp,nyp)
    real :: workrow(nyp)
    real :: workcol(nxp)
    real :: workload(ranks)
    real :: workblock(ranks)
    integer :: jrows(ranks)
    integer :: jrow(ranks)
    integer :: nblocks(ranks)
    real :: relspeed(ranks)
    real :: bfact
    logical :: insuficient
    
    character(len=*), parameter :: h="**(DomainDecompRams)**"

    integer :: irank,i,j,islab,jranks,nslabs,min_blocks,nbigslabs,iblock &
         ,jrank,krank
    real :: aranks,aslabs,totspeed,workdom,workaccum,worksofar &
         ,slabspeed,workslab

    ! default relspeed = 1.0 for ranks of uniform speed.

    relspeed = 1.0

    ! work factor of 1.0 at interior points and .2 at boundary points

    work = 1.0
    bfact=.2
    do j = 1,nyp
       work(1,j) = bfact
       work(nxp,j) = bfact
    enddo
    do i = 1,nxp
       work(i,1) = bfact
       work(i,nyp) = bfact
    enddo

    ! This routine decomposes grid domains of size (nxp,nyp) into a number,
    ! specified by ranks, of rectangular subdomains.  The convention is followed
    ! that any internal boundaries (between subdomains) that are parallel to
    ! the x-axis run continuously across the full domain, while boundaries
    ! parallel to the y-axis may or may not run the full distance across the
    ! domain.  For convenience, regions of the domain bounded by adjacent
    ! east-west internal boundaries are termed "slabs", while smaller divisions
    ! within each slab are termed "blocks".  Each block is required to have
    ! a minimum dimension of 6 by 6 grid cells.  If this cannot be satisfied
    ! with the given input parameters, the subroutine stops.


    ! Estimate the number of slabs to be used (aslabs), and compute a final
    ! nearest integer value (nslabs) which is limited to allowable values.
    ! Zero out array for accumulating number of columns for each rank.

    aranks = float(ranks)
    aslabs = sqrt(aranks * float(nyp) / float(nxp))
    nslabs = min(ranks,max(1,nint(aslabs)))

    totspeed = 0.
    do irank = 1,ranks
       xe(irank) = 0
       totspeed = totspeed + relspeed(irank)
    enddo

    ! Compute total work load over each row and over entire domain.

    workdom = 0.
    do j = 1,nyp
       workrow(j) = 0.
       do i = 1,nxp
          workrow(j) = workrow(j) + work(i,j)
       enddo
       workdom = workdom + workrow(j)
    enddo
    workrow(2) = workrow(2) + workrow(1)
    workrow(nyp-1) = workrow(nyp-1) + workrow(nyp)

    ! Determine number of blocks and the average workload for each slab.

    min_blocks = ranks / nslabs
    nbigslabs = ranks - min_blocks * nslabs
    irank = 0
    do islab = 1,nslabs
       workload(islab) = 0.
       nblocks(islab) = min_blocks
       if (islab .le. nbigslabs) nblocks(islab) = min_blocks + 1
       do iblock = 1,nblocks(islab)
          irank = irank + 1
          workload(islab) = workload(islab)  &
               + workdom * relspeed(irank) / totspeed
       enddo
    enddo

    ! Assign all j-rows to their respective slabs in a way that balances the work
    ! load among slabs according to their respective numbers of ranks (blocks).
    ! The array jrows counts the number of rows in each slab, and the array
    ! jrow is the index of the southernmost row in each slab.

    do islab = 1,nslabs
       jrows(islab) = 0
    enddo

    workaccum = 0.
    worksofar = 0.
    islab = 0

    do j = 2,nyp-1
       workaccum = workaccum + workrow(j)
       if (workaccum - .5 * workrow(j) .gt. worksofar .and.  &
            islab .lt. nslabs) then
          islab = islab + 1
          jrow(islab) = j
          worksofar = worksofar + workload(islab)
       endif
       jrows(islab) = jrows(islab) + 1
    enddo

    irank = 0
    jrank = 0
    krank = 0
    do islab = 1,nslabs

       ! Compute the total work load for each slab and for each i-column in the
       ! slab.

       slabspeed = 0.
       workslab = 0.
       do i = 1,nxp
          workcol(i) = 0.
          do j = jrow(islab),jrow(islab)+jrows(islab)-1
             workcol(i) = workcol(i) + work(i,j)
          enddo
          workslab = workslab + workcol(i)
       enddo
       workcol(2) = workcol(2) + workcol(1)
       workcol(nxp-1) = workcol(nxp-1) + workcol(nxp)

       ! Determine average workload for each block.

       do iblock = 1,nblocks(islab)
          jrank = jrank + 1
          slabspeed = slabspeed + relspeed(jrank)
       enddo
       do iblock = 1,nblocks(islab)
          krank = krank + 1
          workblock(iblock) = workslab  &
               * relspeed(krank) / slabspeed
       enddo

       ! Assign the i-columns of each slab to their respective blocks in a way that
       ! balances the work load among the blocks.  The array ncols counts the number
       ! of i-columns on each rank, and the array ncol is the index of the
       ! westernmost i-column on each rank.

       workaccum = 0.
       worksofar = 0.

       iblock = 0
       do i = 2,nxp-1
          workaccum = workaccum + workcol(i)
          if (workaccum - .5 * workcol(i) .gt. worksofar .and.  &
               iblock .lt. nblocks(islab)) then
             iblock = iblock + 1
             irank = irank + 1
             yb(irank) = jrow(islab)
             xb(irank) = i
             ye(irank) = yb(irank) + jrows(islab) - 1
             worksofar = worksofar + workblock(iblock)
          endif
          xe(irank) = xe(irank) + 1
       enddo
    enddo

    do jrank = 1,ranks
       xe(jrank) = xb(jrank) + xe(jrank) - 1
    enddo

    ! Check to make sure that each subdomain has at least 2 interior 
    ! rows and columns.

    insuficient=.false.
    do jrank = 1,ranks
       if (ye(jrank) - yb(jrank) .lt. 1 .or.  &
            xe(jrank) - xb(jrank) .lt. 1) then
          insuficient=.true.
          write(str(1),"(i8)") nxp
          write(str(2),"(i8)") nyp
          write(str(3),"(i8)") ranks
          write(str(4),"(i8)") jrank
          write(str(5),"(i8)") xb(jrank)
          write(str(6),"(i8)") xe(jrank)
          write(str(7),"(i8)") yb(jrank)
          write(str(8),"(i8)") ye(jrank)
          call MsgDump(h//" process "//trim(adjustl(str(4)))//" received insuficient domain ["//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//" x "//&
               trim(adjustl(str(7)))//":"//trim(adjustl(str(8)))//"]")
       end if
    end do
    if (insuficient) then
       call fatal_error(h//" Partition domain ["//trim(adjustl(str(1)))//" x "//trim(adjustl(str(2)))//"]"//&
            " for "//trim(adjustl(str(3)))//&
            " processes produce sub-domains with less than 2 points in x or y; reduce number of processes")
    end if
  end subroutine DomainDecompRams





  subroutine CheckPartition(nxp, nyp, ranks, empty, &
       xb, xe, yb, ye, ibcon)

    ! verifies if arrays xb, xe, yb, ye
    ! partition domain [2:nxp-1,2:nyp-1] across ranks

    integer, intent(in)  :: nxp
    integer, intent(in)  :: nyp
    integer, intent(in)  :: ranks 
    logical, intent(in)  :: empty(:)
    integer, intent(in)  :: xb(:)
    integer, intent(in)  :: xe(:)
    integer, intent(in)  :: yb(:)
    integer, intent(in)  :: ye(:)
    integer, intent(in) :: ibcon(:)

    integer :: rank
    integer :: ix, xstart, xend
    integer :: iy, ystart, yend
    integer, parameter :: UNASSIGNED=-1
    integer :: owner(nxp,nyp)
    character(len=*), parameter :: h="**(CheckPartition)**"

    owner = UNASSIGNED
    do rank = 1, ranks
       if (.not. empty(rank)) then
          if (btest(ibcon(rank),1)) then
             xstart = 1
          else
             xstart = xb(rank)
          end if
          
          if (btest(ibcon(rank),2)) then
             xend = nxp
          else
             xend = xe(rank)
          end if
          
          if (btest(ibcon(rank),3)) then
             ystart = 1
          else
             ystart = yb(rank)
          end if
          
          if (btest(ibcon(rank),4)) then
             yend = nyp
          else
             yend = ye(rank)
          end if
          
          do iy = ystart, yend
             do ix = xstart, xend
                if (owner(ix,iy) == UNASSIGNED) then
                   owner(ix,iy) = rank
                else
                   write(str(1),"(i8)") ix
                   write(str(2),"(i8)") iy
                   write(str(3),"(i8)") rank
                   write(str(4),"(i8)") owner(ix,iy)
                   call fatal_error(h//"  point "//&
                        "("//trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//")"//&
                        " is assigned to rank "//trim(adjustl(str(3)))//&
                        " and to rank "//trim(adjustl(str(4))))
                end if
             end do
          end do
       end if
    end do
    if (any(owner == UNASSIGNED)) then
       do iy = ystart, yend
          do ix = xstart, xend
             if (owner(ix,iy) == UNASSIGNED) then
                write(str(1),"(i8)") ix
                write(str(2),"(i8)") iy
                call MsgDump(h//" point ("//trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//" is unasigned")
             end if
          end do
       end do
       call fatal_error(h//" there are unassigned domain points")
    end if
  end subroutine CheckPartition





  subroutine MarkBoundary(nxp, nyp, ranks, empty, &
       xb, xe, yb, ye, ibcon)

    ! verify if sub-domain boundaries are global domain 
    ! boundaries or not; fills array ibcon
    
    integer, intent(in)  :: nxp
    integer, intent(in)  :: nyp
    integer, intent(in)  :: ranks
    logical, intent(in)  :: empty(:)
    integer, intent(in)  :: xb(:)
    integer, intent(in)  :: xe(:)
    integer, intent(in)  :: yb(:)
    integer, intent(in)  :: ye(:)
    integer, intent(out) :: ibcon(:)

    integer :: rank
    character(len=*), parameter :: h="**(MarkBoundary)**"

    ! position ranks into domain points [1:nxp,1:nyp] making sure
    ! that ranks do not intersect

    do rank = 1, ranks
       ibcon(rank) = 0
       if (.not. empty(rank)) then
          if (xb(rank) == 2) then
             ibcon(rank)=ibset(ibcon(rank),1)
          end if
          if (xe(rank) == nxp-1) then
             ibcon(rank)=ibset(ibcon(rank),2)
          end if
          if (yb(rank) == 2) then
             ibcon(rank)=ibset(ibcon(rank),3)
          end if
          if (ye(rank) == nyp-1) then
             ibcon(rank)=ibset(ibcon(rank),4)
          end if
       end if
    end do
  end subroutine MarkBoundary
end module ModDomainDecomp
