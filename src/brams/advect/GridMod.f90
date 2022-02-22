module GridMod

  implicit none
  private
  type grid_type
     integer :: nx
     integer :: ny
     integer :: nz
     integer, pointer :: processors(:)
     integer :: lastProcessor   
  end type grid_type
  type(grid_type),dimension(:),allocatable :: grid
  integer :: maxGrids
  integer, parameter :: error=0

  public :: createGrid
  public :: destroyGrid
  public :: setGrid
  public :: getNx
  public :: getNy
  public :: getNz
  public :: addProcessor
  public :: getProcessor
  public :: getLastProcessor
  public :: removeProcessor
  public :: findProcessor
  public :: getMaxGrids


contains





  integer function createGrid(ngrids)
    integer,intent(in) :: ngrids

    integer :: i

    if(allocated(grid)) then
       createGrid=error
    else
       allocate(grid(ngrids))
       maxGrids=ngrids
       do i=1,ngrids
          grid(i)%nx=-1
          grid(i)%ny=-1
          grid(i)%nz=-1
          grid(i)%lastProcessor=0
          if(associated(grid(i)%processors)) then
             nullify(grid(i)%processors)
          end if
       end do
       createGrid=-1
    end if
  end function createGrid





  integer function addProcessor(np,ngrid)
    integer, intent(in) :: np
    integer, intent(in) ::ngrid

    integer,allocatable,dimension(:) :: auxProcessors
    integer :: i

    if(.not. allocated(grid) .or. ngrid>maxGrids) then
       addProcessor=error
    else
       if(.not. associated(grid(ngrid)%processors)) then
          allocate(grid(ngrid)%processors(1))
          grid(ngrid)%lastProcessor=grid(ngrid)%lastProcessor+1
          grid(ngrid)%processors(grid(ngrid)%lastProcessor)=np
          addProcessor=-1	    
       else
          allocate(auxProcessors(grid(ngrid)%lastProcessor))
          do i=1,grid(ngrid)%lastProcessor
             auxProcessors(i)=grid(ngrid)%processors(i)
          end do
          deallocate(grid(ngrid)%processors)
          grid(ngrid)%lastProcessor=grid(ngrid)%lastProcessor+1
          allocate(grid(ngrid)%processors(grid(ngrid)%lastProcessor))
          do i=1,grid(ngrid)%lastProcessor-1
             grid(ngrid)%processors(i)=auxProcessors(i)
          end do
          deallocate(auxProcessors)
          grid(ngrid)%processors(grid(ngrid)%lastProcessor)=np
          addProcessor=-1
       end if
    end if
  end function addProcessor





  integer function removeProcessor(position,ngrid)
    integer, intent(in) :: position
    integer, intent(in) :: ngrid

    integer,allocatable,dimension(:) :: auxProcessors
    integer :: i,k

    if(.not. allocated(grid) .or. &
         .not. associated(grid(ngrid)%processors) .or. &
         position<0 .or. &
         position>grid(ngrid)%lastProcessor .or. &
         ngrid>maxGrids) then
       removeProcessor=error
    else
       allocate(auxProcessors(grid(ngrid)%lastProcessor))
       do i=1,grid(ngrid)%lastProcessor
          auxProcessors(i)=grid(ngrid)%processors(i)
       end do
       deallocate(grid(ngrid)%processors)
       grid(ngrid)%lastProcessor=grid(ngrid)%lastProcessor-1
       allocate(grid(ngrid)%processors(grid(ngrid)%lastProcessor))
       k=0
       do i=1,grid(ngrid)%lastProcessor+1	 
          if(i/=position) then
             k=k+1
             print *, 'Copiando: ',auxProcessors(i),', posicao: ',position
             grid(ngrid)%processors(k)=auxProcessors(i)
          end if
       end do
       removeProcessor=-1
       deallocate(auxProcessors)
    end if
  end function removeProcessor





  integer function findProcessor(np,ngrid)   
    integer,intent(in) :: np
    integer,intent(in) :: ngrid

    integer :: i

    if(.not. allocated(grid) .or. &
         .not. associated(grid(ngrid)%processors) .or. &
         ngrid>maxGrids) then
       findProcessor=-1
    else
       findProcessor=0
       do i=1,grid(ngrid)%lastProcessor
          if(np==grid(ngrid)%processors(i)) then
             findProcessor=i
             exit
          end if
       end do
    end if
  end function findProcessor





  integer function getLastProcessor(ngrid)
    integer,intent(in) :: ngrid

    if(.not. allocated(grid)) then
       getLastProcessor=-1
    else
       getLastProcessor=grid(ngrid)%lastProcessor
    end if
  end function getLastProcessor

  integer function getProcessor(ngrid,position)
    integer, intent(in) :: ngrid
    integer, intent(in) :: position

    if(position<0 .or. &
         position>grid(ngrid)%lastProcessor .or. &
         .not. allocated(grid) &
         .or. ngrid>maxGrids) then
       getProcessor=-1
    else
       getProcessor=grid(ngrid)%processors(position)
    end if
  end function getProcessor





  integer function destroyGrid()
    if(.not. allocated(grid)) then
       destroyGrid=error
    else
       deallocate(grid)
       maxGrids=0
       destroyGrid=-1
    end if
  end function destroyGrid





  integer function setGrid(ngrid,n1,n2,n3)
    integer,intent(in) :: ngrid
    integer,intent(in) :: n1
    integer,intent(in) :: n2
    integer,intent(in) :: n3

    if(.not. allocated(grid) .or. ngrid>maxGrids) then
       setGrid=error
    else
       grid(ngrid)%nx=n1
       grid(ngrid)%ny=n2
       grid(ngrid)%nz=n3
       setGrid=-1
    end if
  end function setGrid





  integer function getNx(ngrid)
    integer,intent(in) :: ngrid

    if(.not. allocated(grid) .or. ngrid>maxGrids) then
       getNx=-1
    else
       getNx=grid(ngrid)%nx
    end if
  end function  getNx





  integer function getNy(ngrid)
    integer,intent(in) :: ngrid

    if(.not. allocated(grid) .or. ngrid>maxGrids) then
       getNy=-1
    else
       getNy=grid(ngrid)%ny
    end if
  end function  getNy





  integer function getNz(ngrid)
    integer,intent(in) :: ngrid

    if(.not. allocated(grid) .or. ngrid>maxGrids) then
       getNz=-1
    else
       getNz=grid(ngrid)%nz
    end if
  end function  getNz





  integer function getMaxGrids()
    getMaxGrids=maxGrids
  end function getMaxGrids
end module GridMod
