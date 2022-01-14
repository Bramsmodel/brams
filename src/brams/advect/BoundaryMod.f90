module BoundaryMod
  implicit none
  private

  type bound_data
     integer :: neighbour
     integer :: firstLine
     integer :: lastLine
     integer :: firstColumn
     integer :: lastColumn
     integer :: firstLevel
     integer :: lastLevel
     integer :: itag
     integer :: operation !1= send, 2=recv
     integer :: axys !1 = x , 2 =y
  end type bound_data

  type bound_type
     type(bound_data),pointer :: boundLimits(:)
     integer :: lastBoundary
     integer :: m2
     integer :: m3
     integer :: ia
     integer :: ja
     integer :: iz
     integer :: jz
  end type bound_type
  type(bound_type),allocatable :: boundary(:,:) !dimension: ngrid,nprocessors

  integer :: maxGrids
  integer :: maxProcessors
  integer,parameter :: error=0
  integer,parameter :: welldone=-1

  public :: createBoundary
  public :: destroyBoundary
  public :: addBoundary
  public :: getfirstLine
  public :: getLastLine
  public :: getFirstColumn
  public :: getLastColumn
  public :: getFirstLevel
  public :: getLastLevel
  public :: dumpBounds
  public :: bound_data
  public :: getLastBoundary
  public :: getNeighbour
  public :: getSize
  public :: getTag
  public :: getOperation
  public :: getAxys
  public :: adjust
  public :: getlimits
  public :: setLimits

contains 





  integer function createBoundary(nGrids,nProcessors)
    integer,intent(in) :: nGrids
    integer,intent(in) :: nProcessors

    integer :: i,j

    if(allocated(boundary)) then
       createBoundary=error
    else
       allocate(boundary(nGrids,nProcessors))
       maxProcessors=nprocessors
       maxGrids=nGrids
       do i=1,nprocessors 
          do j=1,nGrids
             boundary(j,i)%lastBoundary=0
             if(associated(boundary(j,i)%boundLimits))    nullify(boundary(j,i)%boundLimits)
          end do
       end do
       createBoundary=-1
    end if
  end function createBoundary





  integer function destroyBoundary()  
    integer :: i,j

    if(.not. allocated(boundary)) then
       destroyBoundary=error
    else
       do i=1,maxprocessors 
          do j=1,maxGrids
             boundary(j,i)%lastBoundary=0
             if(associated(boundary(j,i)%boundLimits))    nullify(boundary(j,i)%boundLimits)
          end do
       end do
       deallocate(boundary)
       destroyBoundary=-1
    end if
  end function destroyBoundary





  subroutine setLimits(ngrid,np,ia,iz,ja,jz)
    integer, intent(in) :: ngrid
    integer, intent(in) :: np
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz


    boundary(ngrid,np)%ia=ia
    boundary(ngrid,np)%iz=iz
    boundary(ngrid,np)%ja=ja
    boundary(ngrid,np)%jz=jz
  end subroutine setLimits





  subroutine adjust(incx,incy,ngrid,nmachs,m2,m3)
    integer, intent(in) :: ngrid
    integer, intent(in) :: nmachs
    integer, intent(in) :: incx(maxProcessors)
    integer, intent(in) :: incy(maxProcessors)
    integer, intent(in) :: m2(maxProcessors)
    integer, intent(in) :: m3(maxProcessors)

    integer :: np,nb,ia,iz,ja,jz
    do np=1,maxProcessors
       if(nmachs>1) then
          boundary(ngrid,np)%m2=m2(np)
          boundary(ngrid,np)%m3=m3(np)
       else
          boundary(ngrid,np)%m2=boundary(ngrid,np)%iz-boundary(ngrid,np)%ia+3
          boundary(ngrid,np)%m3=boundary(ngrid,np)%jz-boundary(ngrid,np)%ja+3
       end if
    end do
    do np=1,maxProcessors
       boundary(ngrid,np)%ia=boundary(ngrid,np)%ia-incx(np)
       boundary(ngrid,np)%iz=boundary(ngrid,np)%iz-incx(np)
       boundary(ngrid,np)%ja=boundary(ngrid,np)%ja-incy(np)
       boundary(ngrid,np)%jz=boundary(ngrid,np)%jz-incy(np)
       do nb=1,boundary(ngrid,np)%lastBoundary
          boundary(ngrid,np)%boundLimits(nb)%firstLine=boundary(ngrid,np)%boundLimits(nb)%firstLine-incx(np)
          boundary(ngrid,np)%boundLimits(nb)%lastLine=boundary(ngrid,np)%boundLimits(nb)%lastLine-incx(np)
          boundary(ngrid,np)%boundLimits(nb)%firstColumn=boundary(ngrid,np)%boundLimits(nb)%firstColumn-incy(np)
          boundary(ngrid,np)%boundLimits(nb)%lastColumn=boundary(ngrid,np)%boundLimits(nb)%lastColumn-incy(np)
       end do
    end do
  end subroutine adjust





  integer function getLAstBoundary(ngrid,nprocessor)
    integer,intent(in) :: ngrid
    integer,intent(in) :: nProcessor
    if(.not. allocated(boundary) .or. nProcessor>maxProcessors .or. nGrid>maxGrids) then
       getLAstBoundary=0
    else
       getLAstBoundary=boundary(ngrid,nProcessor)%lastBoundary
    end if
  end function getLAstBoundary





  subroutine getlimits(ngrid,nprocessor,m2,m3,ia,iz,ja,jz)
    integer, intent(in) :: ngrid
    integer, intent(in) :: nProcessor
    integer, intent(out) :: m2
    integer, intent(out) :: m3
    integer, intent(out) :: ia
    integer, intent(out) :: iz
    integer, intent(out) :: ja
    integer, intent(out) :: jz

    if(.not. allocated(boundary) .or. nProcessor>maxProcessors .or. nGrid>maxGrids) then
       m2=0
       m3=0
       ia=0
       iz=0
       ja=0
       jz=0
    else
       m2=boundary(ngrid,nprocessor)%m2
       m3=boundary(ngrid,nprocessor)%m3
       ia=boundary(ngrid,nprocessor)%ia
       iz=boundary(ngrid,nprocessor)%iz
       ja=boundary(ngrid,nprocessor)%ja
       jz=boundary(ngrid,nprocessor)%jz
    end if
  end subroutine getLImits





  integer function addBoundary(nGrid,nProcessor,neighbour,localFirstLine,localLastLine, &
       localFirstColumn,localLastColumn,localFirstLevel,localLastLevel, &
       tag,oper,ax)
    integer, intent(in) :: localFirstLine
    integer, intent(in) :: localLastLine
    integer, intent(in) :: neighbour
    integer, intent(in) :: ngrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: localFirstColumn
    integer, intent(in) :: localLastColumn
    integer, intent(in) :: localFirstLevel
    integer, intent(in) :: localLastLevel
    integer, intent(in) :: tag
    integer, intent(in) :: oper
    integer, intent(in) :: ax

    integer :: i
    type(bound_data) :: auxBoundary(boundary(ngrid,nProcessor)%lastBoundary+1)


    if(.not. allocated(boundary) .or. nProcessor>maxProcessors .or. nGrid>maxGrids) then
       addBoundary=error
    else
       if(.not. associated(boundary(ngrid,nProcessor)%boundLimits)) then
          allocate(boundary(ngrid,nProcessor)%boundLimits(1))
          boundary(ngrid,nProcessor)%boundLimits(1)%neighbour=neighbour
          boundary(ngrid,nProcessor)%boundLimits(1)%firstLine=localFirstLine
          boundary(ngrid,nProcessor)%boundLimits(1)%lastLine=localLastLine
          boundary(ngrid,nProcessor)%boundLimits(1)%firstColumn=localFirstColumn
          boundary(ngrid,nProcessor)%boundLimits(1)%lastColumn=localLastColumn
          boundary(ngrid,nProcessor)%boundLimits(1)%firstLevel=localFirstLevel
          boundary(ngrid,nProcessor)%boundLimits(1)%lastLevel=localLastLevel
          boundary(ngrid,nProcessor)%boundLimits(1)%itag=tag
          boundary(ngrid,nProcessor)%boundLimits(1)%operation=oper
          boundary(ngrid,nProcessor)%boundLimits(1)%axys=ax
          boundary(ngrid,nProcessor)%lastBoundary=1

       else
          do i=1,boundary(ngrid,nProcessor)%lastBoundary
             auxBoundary(i)= boundary(ngrid,nProcessor)%boundLimits(i)
          end do
          deallocate(boundary(ngrid,nProcessor)%boundLimits)
          boundary(ngrid,nProcessor)%lastBoundary=boundary(ngrid,nProcessor)%lastBoundary+1
          allocate(boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary))
          do i=1,boundary(ngrid,nProcessor)%lastBoundary-1
             boundary(ngrid,nProcessor)%boundLimits(i)=auxBoundary(i)
          end do
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%neighbour=neighbour
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%firstLine=localFirstLine
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%lastLine=localLastLine
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%firstColumn=localFirstColumn
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%lastColumn=localLastColumn	    
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%firstLevel=localFirstLevel
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%lastLevel=localLastLevel	    
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%itag=tag	    
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%operation=oper	
          boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%axys=ax	
       end if
       addBoundary=-1
    end if
  end function addBoundary





  integer function getNeighbour(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getNeighbour=0
    else
       getNeighbour=boundary(ngrid,nProcessor)%boundLimits(nBound)%neighbour
    end if
  end function getNeighbour





  integer function getAxys(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getAxys=0
    else
       getAxys=boundary(ngrid,nProcessor)%boundLimits(nBound)%axys
    end if
  end function getAxys





  integer function getOperation(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getOperation=-1
    else
       getOperation=boundary(ngrid,nProcessor)%boundLimits(nBound)%Operation
    end if
  end function getOperation





  integer function getTag(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getTag=-1
    else
       getTag=boundary(ngrid,nProcessor)%boundLimits(nBound)%itag
    end if
  end function getTag





  integer function getfirstLine(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getFirstLine=-1
    else
       getFirstLine=boundary(ngrid,nProcessor)%boundLimits(nBound)%firstLine
    end if
  end function getFirstLine





  integer function getLastLine(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getLastLine=-1
    else
       getLastLine=boundary(ngrid,nProcessor)%boundLimits(nBound)%LastLine
    end if
  end function getLastLine





  integer function getFirstColumn(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getFirstColumn=-1
    else
       getFirstColumn=boundary(ngrid,nProcessor)%boundLimits(nBound)%FirstColumn
    end if
  end function getFirstColumn





  integer function getLastColumn(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getLastColumn=-1
    else
       getLastColumn=boundary(ngrid,nProcessor)%boundLimits(nBound)%LastColumn
    end if
  end function getLastColumn





  integer function getFirstLevel(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getFirstLevel=-1
    else
       getFirstLevel=boundary(ngrid,nProcessor)%boundLimits(nBound)%FirstLevel
    end if
  end function getFirstLevel





  integer function getLastLevel(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getLastLevel=-1
    else
       getLastLevel=boundary(ngrid,nProcessor)%boundLimits(nBound)%LastLevel
    end if
  end function getLastLevel





  subroutine dumpBounds()
    integer :: ng,np,nb

    do ng=1,maxGrids
       write (*,FMT='(A,I3.3,A)') '================== Grid ',ng,'=================='; call flush(6)
       write (*,FMT='(8(A3,1X),1X,A5,1X,A8,1X,A1,1X,A1)') 'Npr','Nei','ibg','ien','jbg','jen','zbg','zen','itag ','  Total ','O','D'
       do np=1,maxProcessors
          do nb=1,boundary(ng,np)%lastBoundary
             write (*,FMT='(8(I5,1X),1X,I5.5,1X,I8,1X,I1,1X,I1)') np, &
                  boundary(ng,np)%boundLimits(nb)%neighbour, &
                  boundary(ng,np)%boundLimits(nb)%firstLine, &
                  boundary(ng,np)%boundLimits(nb)%lastLine, &
                  boundary(ng,np)%boundLimits(nb)%firstColumn, &
                  boundary(ng,np)%boundLimits(nb)%lastColumn, &	    
                  boundary(ng,np)%boundLimits(nb)%firstLevel, &
                  boundary(ng,np)%boundLimits(nb)%lastLevel, &
                  boundary(ng,np)%boundLimits(nb)%itag, &
                  getsize(ng,np,nb),getOperation(ng,np,nb), &
                  getAxys(ng,np,nb) ; call flush(6)	    
          end do
       end do
    end do
  end subroutine dumpBounds





  integer function getSize(nGrid,nProcessor,nBound)
    integer, intent(in) :: nGrid
    integer, intent(in) :: nProcessor
    integer, intent(in) :: nBound

    if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
         nBound>boundary(ngrid,nProcessor)%lastBoundary) then
       getSize=0
    else
       getSize=((boundary(ngrid,nProcessor)% boundLimits(nBound)%lastLine - &
            boundary(ngrid,nProcessor)% boundLimits(nBound)%firstLine +1) * &
            (boundary(ngrid,nProcessor)% boundLimits(nBound)%lastColumn - &
            boundary(ngrid,nProcessor)% boundLimits(nBound)%firstColumn +1) * &
            (boundary(ngrid,nProcessor)% boundLimits(nBound)%lastlevel - &
            boundary(ngrid,nProcessor)% boundLimits(nBound)%firstlevel +1))
    end if
  end function getSize
end module BoundaryMod
