module advMessageMod

  ! repository of variables to be used by message passing procedures
  ! to exchange borders of enlarged ghost zone fields

  ! declares and exports required variables
  ! procedures to allocate, deallocate and dump such variables

  implicit none
  private

  type t_message

     ! info for each border exchange
     
     integer, pointer :: Proc(:)
     integer, pointer :: ia(:)
     integer, pointer :: iz(:)
     integer, pointer :: ja(:)
     integer, pointer :: jz(:)
     integer, pointer :: za(:)
     integer, pointer :: zz(:)
     integer, pointer :: mSize(:)
     integer, pointer :: tag(:)
     integer, pointer :: start(:)
  end type t_message

  public :: t_message

  type(t_message), allocatable, public :: SendMessageI(:) !ngrids,nrecv
  type(t_message), allocatable, public :: RecvMessageI(:) !ngrids,nsend
  type(t_message), allocatable, public :: SendMessageJ(:) !ngrids,nrecv
  type(t_message), allocatable, public :: RecvMessageJ(:) !ngrids,nsend

  integer, allocatable, public :: nsendI(:) !ngrids
  integer, allocatable, public :: nrecvI(:) !ngrids
  integer, allocatable, public :: nsendJ(:) !ngrids
  integer, allocatable, public :: nrecvJ(:) !ngrids
  integer, allocatable, public :: totalSendI(:) !ngrids
  integer, allocatable, public :: totalSendJ(:) !ngrids
  integer, allocatable, public :: totalRecvI(:) !ngrids
  integer, allocatable, public :: totalRecvJ(:) !ngrids
  integer, allocatable, public :: newIa(:) !ngrids
  integer, allocatable, public :: newJa(:) !ngrids
  integer, allocatable, public :: newIz(:) !ngrids
  integer, allocatable, public :: newJz(:) !ngrids
  integer, allocatable, public :: newM2(:) !ngrids
  integer, allocatable, public :: newM3(:) !ngrids
  
  logical :: isAllocated
  integer :: TotalMessages

  public :: createMessage
  public :: allocMessages
  public :: dumpMessages
  

contains

  


  
  subroutine createMessage(ngrids)
    integer,intent(in) :: ngrids

    integer :: i   

    if(isAllocated) then
       print *,'Error: Messages already allocated!'
       stop
    end if
    allocate(nsendI(ngrids))
    allocate(nrecvI(ngrids))
    allocate(nsendJ(ngrids))
    allocate(nrecvJ(ngrids))
    allocate(SendMessageI(ngrids))
    allocate(RecvMessageI(ngrids))
    allocate(SendMessageJ(ngrids))
    allocate(RecvMessageJ(ngrids))
    allocate(newM2(ngrids))
    allocate(newM3(ngrids))
    allocate(newIa(ngrids))
    allocate(newIz(ngrids))
    allocate(newJa(ngrids))
    allocate(newJz(ngrids))
    allocate(totalSendI(ngrids))
    allocate(totalSendJ(ngrids))
    allocate(totalRecvI(ngrids))
    allocate(totalRecvJ(ngrids))

    totalSendI=0
    totalSendJ=0
    TotalRecvI=0
    totalRecvJ=0

    do i=1,ngrids
       if(associated(SendMessageI(i)%Proc))  nullify(SendMessageI(i)%Proc)
       if(associated(SendMessageI(i)%ia))    nullify(SendMessageI(i)%ia)    
       if(associated(SendMessageI(i)%iz))	 nullify(SendMessageI(i)%iz)
       if(associated(SendMessageI(i)%ja))	 nullify(SendMessageI(i)%ja)
       if(associated(SendMessageI(i)%jz))	 nullify(SendMessageI(i)%jz)
       if(associated(SendMessageI(i)%za))	 nullify(SendMessageI(i)%za)
       if(associated(SendMessageI(i)%zz))	 nullify(SendMessageI(i)%zz)
       if(associated(SendMessageI(i)%mSize)) nullify(SendMessageI(i)%mSize)
       if(associated(SendMessageI(i)%tag))	 nullify(SendMessageI(i)%tag)
       if(associated(SendMessageI(i)%start))	 nullify(SendMessageI(i)%start)
       if(associated(RecvMessageI(i)%Proc))  nullify(RecvMessageI(i)%Proc)
       if(associated(RecvMessageI(i)%ia))    nullify(RecvMessageI(i)%ia)    
       if(associated(RecvMessageI(i)%iz))	 nullify(RecvMessageI(i)%iz)
       if(associated(RecvMessageI(i)%ja))	 nullify(RecvMessageI(i)%ja)
       if(associated(RecvMessageI(i)%jz))	 nullify(RecvMessageI(i)%jz)
       if(associated(RecvMessageI(i)%za))	 nullify(RecvMessageI(i)%za)
       if(associated(RecvMessageI(i)%zz))	 nullify(RecvMessageI(i)%zz)
       if(associated(RecvMessageI(i)%mSize)) nullify(RecvMessageI(i)%mSize)
       if(associated(RecvMessageI(i)%tag))	 nullify(RecvMessageI(i)%tag)
       if(associated(RecvMessageI(i)%start))	 nullify(RecvMessageI(i)%start)
       if(associated(SendMessageJ(i)%Proc))  nullify(SendMessageJ(i)%Proc)
       if(associated(SendMessageJ(i)%ia))    nullify(SendMessageJ(i)%ia)    
       if(associated(SendMessageJ(i)%iz))	 nullify(SendMessageJ(i)%iz)
       if(associated(SendMessageJ(i)%ja))	 nullify(SendMessageJ(i)%ja)
       if(associated(SendMessageJ(i)%jz))	 nullify(SendMessageJ(i)%jz)
       if(associated(SendMessageJ(i)%za))	 nullify(SendMessageJ(i)%za)
       if(associated(SendMessageJ(i)%zz))	 nullify(SendMessageJ(i)%zz)
       if(associated(SendMessageJ(i)%mSize)) nullify(SendMessageJ(i)%mSize)
       if(associated(SendMessageJ(i)%tag))	 nullify(SendMessageJ(i)%tag)
       if(associated(SendMessageJ(i)%start))	 nullify(SendMessageJ(i)%start)
       if(associated(RecvMessageJ(i)%Proc))  nullify(RecvMessageJ(i)%Proc)
       if(associated(RecvMessageJ(i)%ia))    nullify(RecvMessageJ(i)%ia)    
       if(associated(RecvMessageJ(i)%iz))	 nullify(RecvMessageJ(i)%iz)
       if(associated(RecvMessageJ(i)%ja))	 nullify(RecvMessageJ(i)%ja)
       if(associated(RecvMessageJ(i)%jz))	 nullify(RecvMessageJ(i)%jz)
       if(associated(RecvMessageJ(i)%za))	 nullify(RecvMessageJ(i)%za)
       if(associated(RecvMessageJ(i)%zz))	 nullify(RecvMessageJ(i)%zz)
       if(associated(RecvMessageJ(i)%mSize)) nullify(RecvMessageJ(i)%mSize)
       if(associated(RecvMessageJ(i)%tag))	 nullify(RecvMessageJ(i)%tag)
       if(associated(RecvMessageJ(i)%start))	 nullify(RecvMessageJ(i)%start)
    end do
  end subroutine createMessage




  
  subroutine allocMessages(ngrid,nmess,is2Send,isI)
    integer, intent(in) :: ngrid
    integer, intent(in) :: nmess
    logical, intent(in) :: is2Send
    logical, intent(in) :: isI

    if(is2Send) then
       if(isI) then
          nsendI(ngrid)=nmess
          allocate(SendMessageI(ngrid)%Proc(nmess)) 
          allocate(SendMessageI(ngrid)%ia(nmess))   
          allocate(SendMessageI(ngrid)%iz(nmess))   
          allocate(SendMessageI(ngrid)%ja(nmess))   
          allocate(SendMessageI(ngrid)%jz(nmess))   
          allocate(SendMessageI(ngrid)%za(nmess))   
          allocate(SendMessageI(ngrid)%zz(nmess))   
          allocate(SendMessageI(ngrid)%mSize(nmess))
          allocate(SendMessageI(ngrid)%tag(nmess))   
          allocate(SendMessageI(ngrid)%start(nmess))   
       else
          nsendJ(ngrid)=nmess
          allocate(SendMessageJ(ngrid)%Proc(nmess)) 
          allocate(SendMessageJ(ngrid)%ia(nmess))   
          allocate(SendMessageJ(ngrid)%iz(nmess))   
          allocate(SendMessageJ(ngrid)%ja(nmess))   
          allocate(SendMessageJ(ngrid)%jz(nmess))   
          allocate(SendMessageJ(ngrid)%za(nmess))   
          allocate(SendMessageJ(ngrid)%zz(nmess))   
          allocate(SendMessageJ(ngrid)%mSize(nmess))
          allocate(SendMessageJ(ngrid)%tag(nmess))   
          allocate(SendMessageJ(ngrid)%start(nmess))   
       end if
    else
       if(isI) then
          nrecvI(ngrid)=nmess
          allocate(RecvMessageI(ngrid)%Proc(nmess)) 
          allocate(RecvMessageI(ngrid)%ia(nmess))   
          allocate(RecvMessageI(ngrid)%iz(nmess))   
          allocate(RecvMessageI(ngrid)%ja(nmess))   
          allocate(RecvMessageI(ngrid)%jz(nmess))   
          allocate(RecvMessageI(ngrid)%za(nmess))   
          allocate(RecvMessageI(ngrid)%zz(nmess))   
          allocate(RecvMessageI(ngrid)%mSize(nmess))
          allocate(RecvMessageI(ngrid)%tag(nmess))   
          allocate(RecvMessageI(ngrid)%start(nmess))   
       else
          nrecvJ(ngrid)=nmess
          allocate(RecvMessageJ(ngrid)%Proc(nmess)) 
          allocate(RecvMessageJ(ngrid)%ia(nmess))   
          allocate(RecvMessageJ(ngrid)%iz(nmess))   
          allocate(RecvMessageJ(ngrid)%ja(nmess))   
          allocate(RecvMessageJ(ngrid)%jz(nmess))   
          allocate(RecvMessageJ(ngrid)%za(nmess))   
          allocate(RecvMessageJ(ngrid)%zz(nmess))   
          allocate(RecvMessageJ(ngrid)%mSize(nmess))
          allocate(RecvMessageJ(ngrid)%tag(nmess))   
          allocate(RecvMessageJ(ngrid)%start(nmess))   
       end if
    end if
  end subroutine allocMessages




  
  subroutine dumpMessages(fout,ngrids)
    integer, intent(in) :: fout
    integer, intent(in) :: ngrids

    integer :: ng,m
    do ng=1,ngrids
       write (fout,FMT='(A,I3.3,A)') '============================ Grid ',ng,'=============================='

       write (fout,FMT='("#Sends I: ",I5.5," #Sends J: ",I5.5," #Recs I:		   ",I5.5," #Recs J: ",I5.5 )') &
            nsendI(ng),nsendJ(ng),nrecvI(ng),nrecvJ(ng)	  

       write (fout,FMT='(A)') '=================================================================================='
       if(nSendI(ng)>0) then
	  write (fout,FMT='(A,I8.8,A)') '==== Send data ==== total ',totalSendI(ng),' ==== I direction ===='
	  do m=1,nSendI(ng)
	     write(fout,FMT='("#",I3.3," - Send     to  ",I5.5,1X,6(I3.3,1X),I8.8,1X,I5.5,1X,I5.5)') & 
                  m,SendMessageI(ng)%Proc(m), &
                  SendMessageI(ng)%ia(m) , &    
                  SendMessageI(ng)%iz(m) , &  
                  SendMessageI(ng)%ja(m) , &  
                  SendMessageI(ng)%jz(m) , &  
                  SendMessageI(ng)%za(m) , &  
                  SendMessageI(ng)%zz(m) , &  
                  SendMessageI(ng)%mSize(m), &
                  SendMessageI(ng)%tag(m), &
                  SendMessageI(ng)%start(m)
          end do
       end if
       if(nsendJ(ng)>0) then
          write (fout,FMT='(A,I8.8,A)') '==== Send data ==== total ',totalSendJ(ng),' ==== J direction ===='	  
          do m=1,nSendJ(ng)
	     write(fout,FMT='("#",I3.3," - Send     to  ",I5.5,1X,6(I3.3,1X),I8.8,1X,I5.5,1X,I5.5)') & 
                  m,SendMessageJ(ng)%Proc(m), &
                  SendMessageJ(ng)%ia(m) , &    
                  SendMessageJ(ng)%iz(m) , &  
                  SendMessageJ(ng)%ja(m) , &  
                  SendMessageJ(ng)%jz(m) , &  
                  SendMessageJ(ng)%za(m) , &  
                  SendMessageJ(ng)%zz(m) , &  
                  SendMessageJ(ng)%mSize(m), &
                  SendMessageJ(ng)%tag(m), &
                  SendMessageJ(ng)%start(m)
          end do
       end if
       if(nRecvI(ng)>0) then
          write (fout,FMT='(A,I8.8,A)') '==== Recv data ==== total ',totalRecvI(ng),' ==== I direction ===='	  
	  do m=1,nRecvI(ng)
	     write(fout,FMT='("#",I3.3," - Receive from ",I5.5,1X,6(I3.3,1X),I8.8,1X,I5.5,1X,I5.5)') & 
                  m,RecvMessageI(ng)%Proc(m), &
                  RecvMessageI(ng)%ia(m) , &    
                  RecvMessageI(ng)%iz(m) , &  
                  RecvMessageI(ng)%ja(m) , &  
                  RecvMessageI(ng)%jz(m) , &  
                  RecvMessageI(ng)%za(m) , &  
                  RecvMessageI(ng)%zz(m) , &  
                  RecvMessageI(ng)%mSize(m), &
                  RecvMessageI(ng)%tag(m), &
                  RecvMessageI(ng)%start(m)
          end do
       end if
       if(nRecvJ(ng)>0) then
          write (fout,FMT='(A,I8.8,A)') '==== Recv data ==== total ',totalRecvJ(ng),' ==== J direction ===='	  
	  do m=1,nRecvJ(ng)
	     write(fout,FMT='("#",I3.3," - Receive from ",I5.5,1X,6(I3.3,1X),I8.8,1X,I5.5,1X,I5.5)') & 
                  m,RecvMessageJ(ng)%Proc(m), &
                  RecvMessageJ(ng)%ia(m) , &    
                  RecvMessageJ(ng)%iz(m) , &  
                  RecvMessageJ(ng)%ja(m) , &  
                  RecvMessageJ(ng)%jz(m) , &  
                  RecvMessageJ(ng)%za(m) , &  
                  RecvMessageJ(ng)%zz(m) , &  
                  RecvMessageJ(ng)%mSize(m), &
                  RecvMessageJ(ng)%tag(m), &
                  RecvMessageJ(ng)%start(m)
          end do
       end if
    end do
    call flush(fout)
  end subroutine dumpMessages

end module advMessageMod
