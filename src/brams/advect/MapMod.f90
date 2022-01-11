module mapMod

  private

  type map_type
     integer, pointer :: processor(:,:)     
     integer :: maxX
     integer :: maxY
  end type map_type
  
  type(map_type),allocatable :: map(:)
  
  integer,parameter :: error=0
  integer :: maxMaps

  public :: createmap
  public :: destroymap
  public :: setMap
  public :: setMapProcessor
  public :: getMapProcessor
  public :: getMaxX
  public :: getMaxY
  public :: dumpMap


contains





  integer function createmap(ngrids)
    integer,intent(IN) :: ngrids

    integer :: i

    if(allocated(map)) then
       createmap=error
    else
       allocate(map(ngrids))
       maxmaps=ngrids
       do i=1,ngrids 
          if(associated(map(i)%processor)) nullify(map(i)%processor)
       end do
       createmap=-1
    end if
  end function createmap





  integer function destroymap()
    if(.not. allocated(map)) then
       destroymap=error
    else
       deallocate(map)
       maxmaps=0
       destroymap=-1
    end if
  end function destroymap





  integer function setMap(ngrid,x,y)
    integer,intent(in) :: ngrid
    integer,intent(in) :: x
    integer,intent(in) :: y

    if(.not. allocated(map) .or. associated(map(ngrid)%processor)) then
       setMap=error
    else
       allocate(map(ngrid)%processor(x,y))
       map(ngrid)%maxX=x
       map(ngrid)%maxY=y
       setMap=-1
       map(ngrid)%processor=0
    end if
  end function setMap





  integer function getMaxX(ngrid)
    integer, intent(in) :: ngrid
    if(.not. allocated(map)) then
       getMaxX=-1
    else
       getMaxX= map(ngrid)%maxX
    end if
  end function getMaxX





  integer function getMaxY(ngrid)
    integer, intent(in) :: ngrid
    if(.not. allocated(map)) then
       getMaxY=-1
    else
       getMaxY= map(ngrid)%maxY
    end if
  end function getMaxY




  integer function setMapProcessor(ngrid,xpos,ypos,nprocessor)
    integer, intent(in) :: ngrid
    integer, intent(in) :: xpos
    integer, intent(in) :: ypos
    integer, intent(in) :: nprocessor

    if(.not. allocated(map) .or. .not. associated(map(ngrid)%processor) .or. xpos<0 .or. xpos>map(ngrid)%maxX &
         .or. ypos<0 .or. ypos>map(ngrid)%maxY) then
       setMapProcessor=error
    else
       map(ngrid)%processor(xpos,ypos)=nprocessor
       setMapProcessor=-1
    end if
  end function setMapProcessor





  integer function getMapProcessor(ngrid,xpos,ypos)
    integer, intent(in) :: ngrid
    integer, intent(in) :: xpos
    integer, intent(in) :: ypos

    if(.not. allocated(map) .or. &
         .not. associated(map(ngrid)%processor) .or. &
         xpos<0 .or. &
         xpos>map(ngrid)%maxX .or. &
         ypos<0 .or. &
         ypos>map(ngrid)%maxY) then
       getMapProcessor=-1
    else
       getMapProcessor=map(ngrid)%processor(xpos,ypos)
    end if
  end function getMapProcessor





  subroutine dumpMap()
    integer :: nm,i,j
    character(len=*),parameter :: str1="----------------------------"
    character(len=*),parameter :: str2="-------------"

    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       do nm=1,maxMaps
          write (*,"(a,i3.3,a)") str1//" Map No ",nm," "//str1
          write (*,*) ''
          write (*,'(a,"  ",$)') 'i=    '
          do i=1,map(nm)%maxX
             write(*,'(i3.3,a,$)') i,' '
          end do
          write (*,*) ''
          write (*,*) ''
          do j=1,map(nm)%maxY
             write(*,'(a,i3.3,"  ",$)') 'j= ',j
             do i=1,map(nm)%maxX
                write(*,'(i3.3,a,$)') map(nm)%processor(i,j),' '
             end do
             write(*,'(a)') '';call flush(6)
          end do
          write (*,'(a)') str1//str2//str1
          call flush(6)
       end do
    end if
  end subroutine dumpMap
end module mapMod
