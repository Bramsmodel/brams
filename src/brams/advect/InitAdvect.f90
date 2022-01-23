module AdvectData
  implicit none
  private
  public InitAdvect

contains





  subroutine showNeighbour(ngrids)
    use ProcessorMod, only: dumpNeighbour
    integer, intent(in) :: ngrids
    call dumpNeighbour(ngrids)
  end subroutine showNeighbour





  subroutine showMap()
    use MapMod, only: dumpMap
    call dumpMap()
  end subroutine showMap





  subroutine showBounds()
    use BoundaryMod, only: dumpBounds
    call dumpBounds()
  end subroutine showBounds





  subroutine showMessages(fout,ngrids)
    use advMessageMod, only: dumpMessages
    integer, intent(in) :: fout
    integer, intent(in) :: ngrids
    call dumpMessages(fout,ngrids)
  end subroutine showMessages





  subroutine InitGrid(ngrids,nnxp,nnyp,nnzp)
    use GridMod, only: createGrid
    use GridMod, only: setGrid
    use errorMod, only: printError
    integer, intent(in) :: ngrids
    integer, intent(in) :: nnxp(ngrids)
    integer, intent(in) :: nnyp(ngrids)
    integer, intent(in) :: nnzp(ngrids)

    integer :: i,istat

    if(creategrid(ngrids)==0) istat=printError(1,'Creating grid on Initialization',6)
    do i=1,ngrids
       if(setGrid(i,nnxp(i),nnyp(i),nnzp(i))==0) istat=printError(1,'Setting grid on Initialization',6)
    end do
  end subroutine InitGrid





  subroutine InitAdvect(ngrids,nmachs,mynum,ghostzone,nnxp,nnyp,nnzp,ixb,ixe,iyb,iye)
    use advMessageMod, only: dumpMessages

    integer, intent(in) :: ngrids
    integer, intent(in) :: nmachs
    integer, intent(in) :: mynum
    integer, intent(in) :: ghostZone
    integer, intent(in), optional :: nnxp(ngrids)
    integer, intent(in), optional :: nnyp(ngrids)
    integer, intent(in), optional :: nnzp(ngrids)
    integer, intent(in), optional :: ixb(nmachs,ngrids)
    integer, intent(in), optional :: ixe(nmachs,ngrids)
    integer, intent(in), optional :: iyb(nmachs,ngrids)
    integer, intent(in), optional :: iye(nmachs,ngrids)

    if (mynum==1) then
       call InitGrid(ngrids,nnxp,nnyp,nnzp)
       call InitMap(nmachs,ngrids,ixb,ixe,iyb,iye,nnxp,nnyp,nnzp)
       call showMap()
       call InitBounds(ixb,ixe,iyb,iye,nmachs,ngrids,ghostZone)
    end if
    call sendInit(nmachs,mynum,ngrids)
  end subroutine InitAdvect





  subroutine InitMap(nmachs,ngrids,ixb,ixe,iyb,iye,nnxp,nnyp,nnzp)
    use MapMod, only: createmap,setMap,setMapProcessor
    use ProcessorMod, only: createProcessor
    use GridMod, only: addProcessor,findProcessor
    use errorMod, only: printError

    integer, intent(in) :: nmachs
    integer, intent(in) :: ngrids
    integer, intent(in) :: nnxp(ngrids)
    integer, intent(in) :: nnyp(ngrids)
    integer, intent(in) :: nnzp(ngrids)
    integer, intent(in) :: ixb(nmachs,ngrids)
    integer, intent(in) :: ixe(nmachs,ngrids)
    integer, intent(in) :: iyb(nmachs,ngrids)
    integer, intent(in) :: iye(nmachs,ngrids)

    integer :: i,j,ngrid,nproc,istat

    if(createProcessor(nmachs,ngrids)==0) then
       istat=printError(1,'Creating processor on Initialization',0)
    end if
    if(createmap(ngrids)==0) then
       istat=printError(1,'Creating map on Initialization',0)
    end if
    do i=1,ngrids
       if(setMap(i,nnxp(i),nnyp(i))==0) then
          istat=printError(1,'Setting map on Initialization',0)
       end if
    end do
    do ngrid=1,ngrids
       do nproc=1,nmachs
          do i=ixb(nproc,ngrid),ixe(nproc,ngrid)
             do j=iyb(nproc,ngrid),iye(nproc,ngrid)
                !Putting processor in especific point of grid
                if(setMapProcessor(ngrid,i,j,nproc)==0) then
                   istat=printError(1,'Attributing processor to map on Initialization',6)
                end if
             end do
          end do
       end do
    end do
  end subroutine InitMap





  subroutine InitBounds(ixb,ixe,iyb,iye,nmachs,ngrids,ghostZone)
    use GridMod, only: getMaxGrids,getNx,getNy,getNz
    use MapMod, only: getMapProcessor
    use ProcessorMod, only: pNorth,pSouth,pEast,pWest,setNeighbour,getPNeighbour,getMaxProcessors
    use BoundaryMod, only: createBoundary,destroyBoundary,addBoundary,adjust,setLimits,dumpBounds
    use errorMod, only: printError
    use dump

    include "constants.f90"

    integer, intent(in) :: nmachs
    integer, intent(in) :: ngrids
    integer, intent(in) :: ghostZone
    integer, intent(in) :: ixb(nmachs,ngrids)
    integer, intent(in) :: ixe(nmachs,ngrids)
    integer, intent(in) :: iyb(nmachs,ngrids)
    integer, intent(in) :: iye(nmachs,ngrids)

    integer :: ngrid,i,istep,j,ini,fin,istat,p1,p2,ii,nextJ,jj,jl,il,yini,yfin,xini,xfin,gz,gzi
    integer :: interval,err
    logical,allocatable,dimension(:) :: foundp
    integer :: tagCount
    integer :: incx(nmachs),incy(nmachs)
    integer :: m2(nmachs),m3(nmachs)
    logical :: ThereAreNeighbour(4,nmachs)

    integer,parameter :: pSend=1,pRecv=2,we=1,ew=2,ns=3,sn=4

    !   +------------------------------------------------------------------------+
    !   | 	   N					  ^		       |
    !   | 	   |					  |		       |
    !   |        W <-.-> E   ------> x direction = i	  | y direction = j    |
    !   | 	   |					  |		       |
    !   | 	   S					  |		       |
    !   +------------------------------------------------------------------------+


    if(createBoundary(getMaxGrids(),getMaxProcessors())==0) then
       istat=printError(1,'Creating processor on Initialization',6)
    end if
    allocate(foundp(getMaxProcessors()))
    foundp=.false.
    tagCount=0
    ThereAreNeighbour=.false.
    gz=ghostZone
    gzi=gz
    do ngrid=1,ngrids
       do p1=1,nmachs
          interval=ixe(p1,ngrid)-ixb(p1,ngrid)+2
          gz=min(interval,gz)
          interval=iye(p1,ngrid)-iyb(p1,ngrid)+2
          gz=min(interval,gz)	   	   
       end do
    end do
    if(gzi/=gz) then
       err=dumpMessage(c_tty,c_yes,'InitBounds','ADVC_MNT',c_notice,&
            'Ghostzone required is different of ghostzone possible:',(/gzi,gz/),'I2')
    end if
    do ngrid=1,ngrids
       do p1=1,nmachs
          do i=ixb(p1,ngrid),ixe(p1,ngrid)
	     p2=getMapProcessor(ngrid,i,iye(p1,ngrid)+1)
	     if(p2/=0) then 
	        ThereAreNeighbour(pSouth,p1)=.true.
	        ThereAreNeighbour(pNorth,p2)=.true.
	     end if
	  end do
          do j=iyb(p1,ngrid),iye(p1,ngrid)
	     p2=getMapProcessor(ngrid,ixe(p1,ngrid)+1,j)
	     if(p2/=0) then 
	        ThereAreNeighbour(pEast,p1)=.true.
	        ThereAreNeighbour(pWest,p2)=.true.
	     end if
	  end do
       end do

       do p1=1,nmachs
          incx(p1)=max(ixb(p1,ngrid)-gz-1,0)
          incy(p1)=max(iyb(p1,ngrid)-gz-1,0)

          p2=getMapProcessor(ngrid,ixe(p1,ngrid)+1,iyb(p1,ngrid))
          if(p2/=0) then
             istat=setNeighbour(p1,p2,pEast,ngrid)
             istat=setNeighbour(p2,p1,pWest,ngrid)	    
          end if

          p2=getMapProcessor(ngrid,ixb(p1,ngrid),iye(p1,ngrid)+1)
          if(p2/=0) then
             istat=setNeighbour(p1,p2,pSouth,ngrid)
             istat=setNeighbour(p2,p1,pNorth,ngrid)
          end if
       end do

       do p1=1,nmachs
          if(ThereAreNeighbour(1,p1)) then
             yini=iyb(p1,ngrid)-gz
          else
             yini=iyb(p1,ngrid)-1
          end if

          if(ThereAreNeighbour(2,p1)) then
             yfin=iye(p1,ngrid)+gz
          else
             yfin=iye(p1,ngrid)+1
          end if

          if(ThereAreNeighbour(4,p1)) then
             xini=ixb(p1,ngrid)-gz
          else
             xini=ixb(p1,ngrid)-1
          end if

          if(ThereAreNeighbour(3,p1)) then
             xfin=ixe(p1,ngrid)+gz
          else
             xfin=ixe(p1,ngrid)+1
          end if

          m2(p1)=xfin-xini+1
          m3(p1)=yfin-yini+1
       end do

       do p1=1,nmachs
          !right (east) and left (west) boundaries
          if(getPNeighbour(p1,1,ngrid)/=0) then
             yini=iyb(p1,ngrid)-gz
          else
             yini=iyb(p1,ngrid)-1
          end if
          p2=getMapProcessor(ngrid,ixe(p1,ngrid)+1,iyb(p1,ngrid))
          if(p2/=0) then
             if(iye(p2,ngrid)>iye(p1,ngrid)) then
                if(getPNeighbour(p2,2,ngrid)/=0) then
                   yfin=iye(p2,ngrid)+gz
                else
                   yfin=iye(p2,ngrid)+1
                end if
                !P1->P2
                tagCount=tagCount+1
                istat=addBoundary(nGrid,p1,p2,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid), &
                     yini,yfin,1,getNz(ngrid),tagCount,pSend,we)
                istat=addBoundary(nGrid,p2,p1,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid), &
                     yini,yfin,1,getNz(ngrid),tagCount,pRecv,we)		       
                !P2->P1
                tagCount=tagCount+1
                istat=addBoundary(nGrid,p2,p1,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz, &
                     yini,yfin,1,getNz(ngrid),tagCount,pSend,ew)
                istat=addBoundary(nGrid,p1,p2,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz, &
                     yini,yfin,1,getNz(ngrid),tagCount,pRecv,ew)		       
             else
                tagCount=tagCount+1
                if(getPNeighbour(p1,2,ngrid)/=0) then
                   yfin=iye(p1,ngrid)+gz
                else
                   yfin=iye(p1,ngrid)+1
                end if
                istat=addBoundary(nGrid,p1,p2,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid), &
                     yini,yfin,1,getNz(ngrid),tagCount,pSend,we)
                istat=addBoundary(nGrid,p2,p1,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid), &
                     yini,yfin,1,getNz(ngrid),tagCount,pRecv,we)
                tagCount=tagCount+1     
                istat=addBoundary(nGrid,p2,p1,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz, &
                     yini,yfin,1,getNz(ngrid),tagCount,pSend,ew)
                istat=addBoundary(nGrid,p1,p2,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz, &
                     yini,yfin,1,getNz(ngrid),tagCount,pRecv,ew)
             end if
          end if

          !upper (north) and lower (south) boundaries
          if(getPNeighbour(p1,4,ngrid)/=0) then
             xini=ixb(p1,ngrid)-gz
          else
             xini=ixb(p1,ngrid)-1
          end if
          p2=getMapProcessor(ngrid,ixb(p1,ngrid),iye(p1,ngrid)+1)
          if(p2/=0) then
             if(ixe(p2,ngrid)<ixe(p1,ngrid)) then
	        tagCount=tagCount+1
                if(getPNeighbour(p2,3,ngrid)/=0) then
	           xfin=ixe(p2,ngrid)+gz
		else
		   xfin=ixe(p2,ngrid)+1
		end if
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1-gz,iye(p1,ngrid), &
                     1,getNz(ngrid),tagCount,pSend,ns)
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1-gz, &
                     iye(p1,ngrid),1,getNz(ngrid),tagCount,pRecv,ns)
		tagCount=tagCount+1
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1, &
                     iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pSend,sn)
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1,iye(p1,ngrid)+gz, &
                     1,getNz(ngrid),tagCount,pRecv,sn)
	        p2=getMapProcessor(ngrid,ixe(p2,ngrid)+1,iye(p1,ngrid)+1)
                if(getPNeighbour(p2,4,ngrid)/=0) then
	           xini=ixb(p2,ngrid)-gz
		else
	           xini=ixb(p2,ngrid)+1
		end if
                if(getPNeighbour(p1,3,ngrid)/=0) then
	           xfin=ixe(p1,ngrid)+gz
		else
	           xfin=ixe(p1,ngrid)+1
		end if
	        tagCount=tagCount+1
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1-gz, &
                     iye(p1,ngrid),1,getNz(ngrid),tagCount,pSend,ns)
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1-gz, &
                     iye(p1,ngrid),1,getNz(ngrid),tagCount,pRecv,ns)
		tagCount=tagCount+1       
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1, &
                     iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pSend,sn)
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1, &
                     iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pRecv,sn)
             else
                if(getPNeighbour(p1,3,ngrid)/=0) then
	           xfin=ixe(p1,ngrid)+gz
		else
		   xfin=ixe(p1,ngrid)+1
		end if
	        tagCount=tagCount+1
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1-gz, &
                     iye(p1,ngrid),1,getNz(ngrid),tagCount,pSend,ns)
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1-gz, &
                     iye(p1,ngrid),1,getNz(ngrid),tagCount,pRecv,ns)
		tagCount=tagCount+1
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1, &
                     iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pSend,sn)
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1, &
                     iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pRecv,sn)
             end if
	  end if
       end do

       do p1=1,nmachs
          call SetLimits(ngrid,p1,ixb(p1,ngrid),ixe(p1,ngrid),iyb(p1,ngrid),iye(p1,ngrid))
       end do

       call adjust(incx,incy,ngrid,nmachs,m2,m3)
    end do
  end subroutine InitBounds





  subroutine sendInit(nmachs,mynum, ngrids)
    use ProcessorMod, only: getPNeighbour,setNeighbour
    use BoundaryMod, only: bound_data,getLAstBoundary,getNeighbour,getfirstLine,getLastLine, &
         getFirstColumn,getLastColumn,getFirstLevel,getLastLevel, &
         getTag,getSize,getOperation,getAxys,getLimits
    use advMessageMod, only: SendMessageI,RecvMessageI,createMessage,allocMessages, &
         newM2,newM3,newIa,newIz,newJa,newJz,SendMessageJ,RecvMessageJ, &
         totalSendI,totalSendJ,TotalRecvI,totalRecvJ			   
    include "mpif.h"

    integer, intent(in) :: nmachs
    integer, intent(in) :: mynum
    integer, intent(in) :: ngrids

    integer :: np,dir,istat,ierr,soc,i,ng,lb,nb,isi,iri,gz
    integer :: socl,iia,iiz,ija,ijz,isj,irj,lbl
    integer :: nSEndI,nsendJ,nRecvI,nRecvJ
    logical :: ing(4)
    integer :: handle(nmachs)
    integer :: status(MPI_STATUS_SIZE)
    integer,allocatable,dimension(:) :: datacomm,datacomm_local
    integer :: l_mynum

    integer :: info(9)
    integer :: newM2l(ngrids),newM3l(ngrids),newIal(ngrids)
    integer :: newIzl(ngrids),newJal(ngrids),newJzl(ngrids)



    l_mynum=mynum-1
    call createMessage(ngrids)
    if(l_mynum==0) then  !IF I am the master
       do np=0,nmachs-1
          do ng=1,ngrids
             lb=getLAstBoundary(ng,np+1)
             soc=lb*11
             if(np==0) socl=soc
             info(1)=np+1
             info(2)=ng
             info(3)=lb
             call getLimits(ng,np+1,info(4),info(5),info(6),info(7),info(8),info(9))
             if(np/=0) & !LFR-5.0  
                  call MPI_send(info,9, MPI_INTEGER,np,20000+np, &
                  MPI_COMM_WORLD, ierr)
             allocate(dataComm(soc))
             if(np==0) allocate(dataComm_local(socl))               
             i=1
             do nb=1,lb
                dataComm(i)=getNeighbour(ng,np+1,nb)
                dataComm(i+1)=getfirstLine(ng,np+1,nb)
                dataComm(i+2)=getLAstLine(ng,np+1,nb)
                dataComm(i+3)=getfirstColumn(ng,np+1,nb)
                dataComm(i+4)=getLastColumn(ng,np+1,nb)
                dataComm(i+5)=getfirstLevel(ng,np+1,nb)
                dataComm(i+6)=getLastLevel(ng,np+1,nb)
                dataComm(i+7)=getSize(ng,np+1,nb)
                dataComm(i+8)=getTag(ng,np+1,nb)
                dataComm(i+9)=getOperation(ng,np+1,nb)
                dataComm(i+10)=getAxys(ng,np+1,nb)
                if(np==0) then
                   dataComm_local(i)=getNeighbour(ng,np+1,nb)
                   dataComm_local(i+1)=getfirstLine(ng,np+1,nb)
                   dataComm_local(i+2)=getLAstLine(ng,np+1,nb)
                   dataComm_local(i+3)=getfirstColumn(ng,np+1,nb)
                   dataComm_local(i+4)=getLastColumn(ng,np+1,nb)
                   dataComm_local(i+5)=getfirstLevel(ng,np+1,nb)
                   dataComm_local(i+6)=getLastLevel(ng,np+1,nb)
                   dataComm_local(i+7)=getSize(ng,np+1,nb)
                   dataComm_local(i+8)=getTag(ng,np+1,nb)
                   dataComm_local(i+9)=getOperation(ng,np+1,nb)
                   dataComm_local(i+10)=getAxys(ng,np+1,nb)
                end if
                i=i+11
	     end do
             if(np/=0) & 
                  call MPI_send(dataComm,soc, MPI_INTEGER,np,30000+np, &
                  MPI_COMM_WORLD, ierr)
             if(np==0) then !copia para variaveis locais ao inves de enviar
                lbl=info(3)
                newM2l(ng)=info(4)
                newM3l(ng)=info(5)
                newIal(ng)=info(6)
                newIzl(ng)=info(7)
                newJal(ng)=info(8)
                newJzl(ng)=info(9)
             end if
	     deallocate(dataComm)
	  end do
       end do
    end if

    do ng=1,ngrids
       if(l_mynum/=0) then  
          call MPI_recv(info,9, MPI_INTEGER,0,20000+l_mynum, & 
               MPI_COMM_WORLD, status, ierr)
	  newM2(ng)=info(4)
	  newM3(ng)=info(5)
	  newIa(ng)=info(6)
	  newIz(ng)=info(7)
	  newJa(ng)=info(8)
	  newJz(ng)=info(9)
	  soc=info(3)*11 !Number of boundaries
          lb=info(3)
       else
          newM2(ng)=newM2l(ng)
          newM3(ng)=newM3l(ng)
          newIa(ng)=newIal(ng)
          newIz(ng)=newIzl(ng)
          newJa(ng)=newJal(ng)
          newJz(ng)=newJzl(ng)
          lb=lbl
       end if

       allocate(datacomm(soc))

       if(l_mynum/=0) then 
          call MPI_recv(datacomm,soc, MPI_INTEGER,0,30000+l_mynum, & 
               MPI_COMM_WORLD, status, ierr)
       else !Se eh o processador 1 apenas copia os dados pois eh local
          do i=1,socl
             datacomm(i)=datacomm_local(i)
          end do
       end if

       !Counting the messages by types
       nRecvI=0;nRecvJ=0;nSendI=0;nSendJ=0
       i=1
       do nb=1,lb
          if(dataComm(i+9)==1) then !Send
             if(dataComm(i+10)==1 .or. dataComm(i+10)==2) then !x axys
                nsendI=nsendI+1
             else ! y axys
                nsendJ=nsendJ+1
             end if
          else
             if(dataComm(i+10)==1 .or. dataComm(i+10)==2) then !x axys
                nRecvI=nRecvI+1
             else ! y axys
                nRecvJ=nRecvJ+1
             end if
          end if
          i=i+11
       end do

       call allocMessages(ng,nSendI,.true.,.true.) !IstoSend, IsI
       call allocMessages(ng,nRecvI,.false.,.true.) !IstoRecv, IsI
       call allocMessages(ng,nSEndJ,.true.,.false.) !IstoSend, IsJ
       call allocMessages(ng,nRecvJ,.false.,.false.) !IstoREcv, IsJ
       i=1
       isi=0;isj=0
       iri=0;irj=0
       do nb=1,lb
          if(dataComm(i+9)==1) then !SEnd
             if(dataComm(i+10)==1 .or. dataComm(i+10)==2) then !x axys
                isi=isi+1
                SendMessageI(ng)%Proc(isi)= dataComm(i+0)
                SendMessageI(ng)%ia(isi)=	  dataComm(i+1)
                SendMessageI(ng)%iz(isi)=	  dataComm(i+2)
                SendMessageI(ng)%ja(isi)=	  dataComm(i+3)
                SendMessageI(ng)%jz(isi)=	  dataComm(i+4)
                SendMessageI(ng)%za(isi)=	  dataComm(i+5)
                SendMessageI(ng)%zz(isi)=	  dataComm(i+6)
                SendMessageI(ng)%mSize(isi)=dataComm(i+7)
                SendMessageI(ng)%tag(isi)=  dataComm(i+8)	
                SendMessageI(ng)%mSize(isi)=(SendMessageI(ng)%iz(isi)- &
                     SendMessageI(ng)%ia(isi)+1)* &
                     (SendMessageI(ng)%jz(isi) - &
                     SendMessageI(ng)%ja(isi)+1) * &
                     ( SendMessageI(ng)%zz(isi) - &
                     SendMessageI(ng)%za(isi)+1)
                SendMessageI(ng)%start(isi)=totalsendI(ng)+1
                totalsendI(ng)=totalsendI(ng)+SendMessageI(ng)%mSize(isi)
             else

                isj=isj+1
                SendMessageJ(ng)%Proc(isj)= dataComm(i+0)
                SendMessageJ(ng)%ia(isj)=	  dataComm(i+1)
                SendMessageJ(ng)%iz(isj)=	  dataComm(i+2)
                SendMessageJ(ng)%ja(isj)=	  dataComm(i+3)
                SendMessageJ(ng)%jz(isj)=	  dataComm(i+4)
                SendMessageJ(ng)%za(isj)=	  dataComm(i+5)
                SendMessageJ(ng)%zz(isj)=	  dataComm(i+6)
                SendMessageJ(ng)%mSize(isj)=dataComm(i+7)
                SendMessageJ(ng)%tag(isj)=  dataComm(i+8)	
                SendMessageJ(ng)%mSize(isj)=(SendMessageJ(ng)%iz(isj)- &
                     SendMessageJ(ng)%ia(isj)+1)* &
                     (SendMessageJ(ng)%jz(isj) - &
                     SendMessageJ(ng)%ja(isj)+1) * &
                     ( SendMessageJ(ng)%zz(isj) - &
                     SendMessageJ(ng)%za(isj)+1)
                SendMessageJ(ng)%start(isj)=totalsendJ(ng)+1
                totalsendJ(ng)=totalsendJ(ng)+SendMessageJ(ng)%mSize(isj)		
             end if
          else if(dataComm(i+9)==2) then !Recv
             if(dataComm(i+10)==1 .or. dataComm(i+10)==2) then !x axys
                iri=iri+1
                RecvMessageI(ng)%Proc(iri)= dataComm(i+0)
                RecvMessageI(ng)%ia(iri)=	dataComm(i+1)
                RecvMessageI(ng)%iz(iri)=	dataComm(i+2)
                RecvMessageI(ng)%ja(iri)=	dataComm(i+3)
                RecvMessageI(ng)%jz(iri)=	dataComm(i+4)
                RecvMessageI(ng)%za(iri)=	     dataComm(i+5)
                RecvMessageI(ng)%zz(iri)=	     dataComm(i+6)
                RecvMessageI(ng)%mSize(iri)=dataComm(i+7)
                RecvMessageI(ng)%tag(iri)=  dataComm(i+8)			   
                RecvMessageI(ng)%mSize(iri)=(RecvMessageI(ng)%iz(iri)- &
                     RecvMessageI(ng)%ia(iri)+1)* &
                     (RecvMessageI(ng)%jz(iri) - &
                     RecvMessageI(ng)%ja(iri)+1) * &
                     ( RecvMessageI(ng)%zz(iri) - &
                     RecvMessageI(ng)%za(iri)+1)
                RecvMessageI(ng)%start(iri)=totalRecvI(ng)+1
                totalRecvI(ng)=totalRecvI(ng)+RecvMessageI(ng)%mSize(iri)
             else
                irj=irj+1
                RecvMessageJ(ng)%Proc(irj)= dataComm(i+0)
                RecvMessageJ(ng)%ia(irj)=	dataComm(i+1)
                RecvMessageJ(ng)%iz(irj)=	dataComm(i+2)
                RecvMessageJ(ng)%ja(irj)=	dataComm(i+3)
                RecvMessageJ(ng)%jz(irj)=	dataComm(i+4)
                RecvMessageJ(ng)%za(irj)=	     dataComm(i+5)
                RecvMessageJ(ng)%zz(irj)=	     dataComm(i+6)
                RecvMessageJ(ng)%mSize(irj)=dataComm(i+7)
                RecvMessageJ(ng)%tag(irj)=  dataComm(i+8)			   
                RecvMessageJ(ng)%mSize(irj)=(RecvMessageJ(ng)%iz(irj)- &
                     RecvMessageJ(ng)%ia(irj)+1)* &
                     (RecvMessageJ(ng)%jz(irj) - &
                     RecvMessageJ(ng)%ja(irj)+1) * &
                     ( RecvMessageJ(ng)%zz(irj) - &
                     RecvMessageJ(ng)%za(irj)+1)

                RecvMessageJ(ng)%start(irj)=totalRecvJ(ng)+1
                totalRecvJ(ng)=totalRecvJ(ng)+RecvMessageJ(ng)%mSize(irj)
             end if
          end if
          i=i+11
       end do
       deallocate(datacomm)
    end do
  end subroutine sendInit
end module AdvectData
