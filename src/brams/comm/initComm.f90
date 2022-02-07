module ModComm
  !# Auxiliary Comunication Module for MPI extra halo comm
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**:Auxiliary Comunication Module for MPI extra halo comm:
  !# The original BRAM'S Models ghostzone (halo) is with 1
  !# position for each side (N,S,E,W). To use the models with high
  !# order advection schemes is need to extend the ghostzone for more
  !# positions. In this case the module extend in 3 more positions in
  !# order to BRAMS use the Runge-Kupta scheme.
  !# This module contains routines that initialize a new types for communication, copy the
  !# original vars inside the new borders extended vars, extend the
  !# first and last originals position in extend positions (Because
  !# the model has border without communication) and do the parallel
  !# MPI communication among processors
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
  !#
  !# **Date**: 2018Jan
  !# @endnote
  !#
  !# @changes
  !#
  !# +
  !# @endchanges
  !# @bug
  !# No active bugs reported now
  !# @endbug
  !#
  !# @todo
  !# 1. Convert this module to use the original BRAMS MPI communication &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !# @endwarning
  !#
  !#--- ----------------------------------------------------------------------------------------
  use ReadBcst, only: & !Just for broadcast comm infos
       Broadcast
  use ParLib  , only: & !Subroutines for parallel comunications
       parf_send_noblock_real, &
       parf_get_noblock_real , &
       parf_wait_any_nostatus, &
       parf_wait_all_nostatus
  use mem_grid, only: time
  use ModParallelEnvironment, only: MsgDump
  
  implicit none
  private
  public initExtraComm
  public copyMyPart
  public commHalo
  public commHaloAcou
  public expandBorder
  
  integer, public, parameter :: west=1
  !# number for west direction
  integer, public, parameter :: east=2
  !# number for east direction
  integer, public, parameter :: north=3
  !# number for north direction
  integer, public, parameter :: south=4
  !# number for south direction
  logical, public, allocatable  :: border(:,:)
  !# if exist border in that direction (E,W,N,S) is .true.
  integer, public :: i0LGZ
  ! offset on x from local to global index on partitions with GhostZoneLength
  integer, public :: j0LGZ
  ! offset on y from local to global index on partitions with GhostZoneLength
  
  ! module private variables
  
  type snd
     integer :: dest
     !# message destination processor to send
     integer :: tag
     !# tag that identify the message to send
     integer :: length
     !# total length of one var message to send
     integer :: i0
     !# local to global index offset at x direction for this proc
     integer :: xbeg
     !# local index of initial x position (lon) to send
     integer :: xend
     !# local index of final x position to send
     integer :: j0
     !# local to global index offset at y direction for this proc
     integer :: yBeg
     !# local index of initial y position (lat) to send
     integer :: yEnd
     !# local index of final y position to send
     real,pointer :: dataMess(:)
     !# Message to be sent
  end type snd

  type rcv
     integer :: from
     !# message origin processor to receive
     integer :: tag
     !# tag that identify the message to receive
     integer :: length
     !# total length of one var message to receive
     integer :: i0
     !# local to global index offset at x direction for this proc
     integer :: xbeg
     !# initial x position (Lon) to receive
     integer :: xend
     !# final x position to receive
     integer :: j0
     !# local to global index offset at y direction for this proc
     integer :: yBeg
     !# initial y position (lat) to receive
     integer :: yEnd
     !# final y position to receive
     real,pointer :: dataMess(:)
     !# Message to be received to receive
  end type rcv

  type(snd),allocatable :: send(:,:)
  !# The composit send info for runge-kupta (nproc,nmess/proc)
  type(snd),allocatable :: send_acou(:,:)
  !# The composit send info for acoustic (nproc,nmess/proc)
  type(rcv),allocatable :: receive(:,:)
  !# The composit receive info for runge-kupta (nproc,nmess/proc)
  type(rcv),allocatable :: receive_acou(:,:)
  !# The composit receive info for acoustic (nproc,nmess/proc)
  integer, allocatable  :: nMess(:)
  !# the total of messages per processor

  include "mpif.h"

contains



  
  subroutine initExtraComm(nmachs,mynum,GhostZoneLength,&
       nnxp,nnyp,nnzp,&
       ixb,ixe,iyb,iye,master_num, &
       nodei0,nodej0,nodemxp,nodemyp,nodemzp)
    !# Initialize auxiliary Comunication
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine perform the initialization of auxiliary communication. Count the
    !# total of message for each send/receive and prepare the types
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !# 1. Convert this module to use the original BRAMS MPI communication &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, parameter  :: ngrid=1
    !# total of grids - set to 1 by now
    integer, intent(in) :: nmachs
    !# total of processor
    integer, intent(in) :: mynum
    !# # of this processor
    integer, intent(in) :: GhostZoneLength
    !# size of ghostzone
    integer, intent(in) :: master_num
    !# # of master processor
    integer, intent(in) :: nnxp(ngrid)
    !# global total of points in x
    integer, intent(in) :: nnyp(ngrid)
    !# global total of points in y
    integer, intent(in) :: nnzp(ngrid)
    !# global total of points in z
    integer, intent(in) :: ixb(nmachs,ngrid)
    !# First position in x dir
    integer, intent(in) :: ixe(nmachs,ngrid)
    !# Last position in x dir
    integer, intent(in) :: iyb(nmachs,ngrid)
    !# First position in y dir
    integer, intent(in) :: iye(nmachs,ngrid)
    !# Last position in x dir
    integer, intent(in) :: nodei0(nmachs,ngrid)
    !# increment for i in this proc
    integer, intent(in) :: nodej0(nmachs,ngrid)
    !# increment for j in this proc
    integer, intent(in) :: nodemxp(nmachs,ngrid)
    !# total of points in x dir
    integer, intent(in) :: nodemyp(nmachs,ngrid)
    !# total of points in y dir
    integer, intent(in) :: nodemzp(nmachs,ngrid)
    !# total of points in z dir

    integer :: p1
    !# local for proc #1
    integer :: p2
    !# local for proc #2
    character :: cmyn

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(initExtraComm)**"
    character(len=8) :: str(10)
    
    if (nmachs<=7) write(cmyn,fmt='(I1)') mynum
    allocate (nmess(nmachs)) ! Allocate the array of count messages with total proc
    allocate (border(nmachs,4)) ! Allocate logical array to determine CC in directions
    border=.true. !In case of 1 proc there are all borders

    if(dumpLocal) then 
       write(str(1),"(i8)") nmachs
       call MsgDump(h//" allocates nmess("//trim(adjustl(str(1)))//")"//&
            " and border("//trim(adjustl(str(1)))//",4)")
    endif

    if(nmachs==1) return

    if(mynum==1) then !If I am the 'master'
       allocate(send(nmachs,nmachs)) !allocate the send
       allocate(receive(nmachs,nmachs)) !allocate the receive
       allocate(send_acou(nmachs,nmachs)) !allocate the send
       allocate(receive_acou(nmachs,nmachs)) !allocate the receive
       call fillSendReceive(nmachs,ngrid,ixb,ixe,iyb,iye,&
            nMess,nodei0,nodej0,nnxp,nnyp) !fill the arrays with neighbour
    endif
    call SendandSetAll(master_num,mynum,nmachs) !Send the info to all processors


    ! global index = local index + offset
    ! offset for partition with GhostZoneLength only for this MPI rank:
    if (border(mynum, west)) then
       i0LGZ = nodei0(mynum,1)
    else
       i0LGZ = nodei0(mynum,1)-GhostZoneLength+1
    end if
    if (border(mynum, north)) then
       j0LGZ = nodej0(mynum,1)
    else
       j0LGZ = nodej0(mynum,1)-GhostZoneLength+1
    end if
    
!!$    call dumpComm(nodemxp(mynum,1),nodemyp(mynum,1),nodemzp(mynum,1), &
!!$         myNum,nmachs,nodei0,nodej0) !dump for test (only for less than 7 processors)

    if(dumpLocal) then 
       write(str(1),"(i8)") nMess(mynum)
       call MsgDump(h//" this processor has "//trim(adjustl(str(1)))//&
            " MPI border exchanges with extended ghost zone for each variable")
       write(str(1),"(i8)") i0LGZ
       write(str(2),"(i8)") j0LGZ
       call MsgDump(h//" offsets local to global:"//&
            " i0LGZ="//trim(adjustl(str(1)))//&
            "; j0LGZ="//trim(adjustl(str(2))))
    endif
  end subroutine initExtraComm




  
  subroutine copyMyPart(scp,scr,ufx_local,vfx_local,wfx_local, &
       ufx,vfx,wfx,mxp,myp,mzp,isi,js,ks, &
       mzi,mzpp3,mxi,mxpp3,myi,mypp3,vname)
    !# Copy the internal part of arrays inside extended arrays
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine perform a copy of the internal parts of arrays (with no extended halo)
    !# in the extend arrays.
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(in) :: mxp
    !# points in x dir
    integer, intent(in) :: myp
    !# points in y dir
    integer, intent(in) :: mzp
    !# points in z dir
    integer, intent(in) :: isi
    !# is=1 -> i direction of advection
    integer, intent(in) :: js
    !# js=1 -> j direction of advection
    integer, intent(in) :: ks
    !# ks=1 -> k direction of advection
    integer, intent(in) :: mzi
    !# initial z position of extended array (-2)
    integer, intent(in) :: mzpp3
    !# final z position of extended array (mzp+3)
    integer, intent(in) :: mxi
    !# initial x position of extended array (-2)
    integer, intent(in) :: mxpp3
    !# final x position of extended array (mxp+3)
    integer, intent(in) :: myi
    !# initial y position of extended array (-2)
    integer, intent(in) :: mypp3
    !# final y position of extended array (myp+3)
    real,intent(in)     :: scp(mzp,mxp,myp)
    !# Original Brams' size var to be advected
    real,intent(in)     :: ufx(mzp,mxp,myp)
    !# Original Brams' size rhou*U
    real,intent(in)     :: vfx(mzp,mxp,myp)
    !# Original Brams' size rhou*V
    real,intent(in)     :: wfx(mzp,mxp,myp)
    !# Original Brams' size rhou*W
    character(len=*), intent(in) :: vname
    !# name of variable that will be advected
    real,intent(out)    :: scr(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var
    real,intent(out)    :: ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var rhou*U
    real,intent(out)    :: vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var rhou*V
    real,intent(out)    :: wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var rhou*W

    integer :: i,j,k !

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(copyMyPart)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") size(scp,1)
       write(str(2),"(i8)") size(scp,2)
       write(str(3),"(i8)") size(scp,3)
       write(str(4),"(i8)") size(scr,1)
       write(str(5),"(i8)") size(scr,2)
       write(str(6),"(i8)") size(scr,3)
       call MsgDump(h//" copy the internal part of "//vname//"("//&
            trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//") into scr("//&
            trim(adjustl(str(4)))//","//trim(adjustl(str(5)))//","//trim(adjustl(str(6)))//")")
    end if

    ! Copy the original (mzp,mxp,myp) vars to extended vars (-2:mzp+3,-2:mxp+3,-2:myp+3)
    ! only the internal values (without extra borders)
    do j=1,myp !-1+js
       do i=1,mxp !-1+isi
          do k=1,mzp !-1+ks
             scr(k,i,j)=scp(k,i,j)
          enddo
       enddo
    enddo
    do j=1,myp
       do i=1,mxp
          do k=1,mzp
             ufx_local(k,i,j)=ufx(k,i,j)
             vfx_local(k,i,j)=vfx(k,i,j)
             wfx_local(k,i,j)=wfx(k,i,j)
          enddo
          ufx_local(mzp+1,i,j)=0.
          vfx_local(mzp+1,i,j)=0.
          wfx_local(mzp+1,i,j)=0.
       enddo
    enddo
  end subroutine copyMyPart




  subroutine expandBorder(mxp,myp,mzp,is,js,ks, &
       mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
       scr,ufx_local,vfx_local,wfx_local)
    !# Copy the original halo (1) for extended halo if in model's border
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine perform a copy of the original halo of 1 (first and last position)
    !# if it is in model border and if variable=.false. (variable is in radvc_rk routine)
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(in) :: mxp
    !# points in x dir
    integer, intent(in) :: myp
    !# points in y dir
    integer, intent(in) :: mzp
    !# points in z dir
    integer, intent(in) :: is
    !# is=1 -> i direction of advection
    integer, intent(in) :: js
    !# js=1 -> j direction of advection
    integer, intent(in) :: ks
    !# ks=1 -> k direction of advection
    integer, intent(in) :: mzi
    !# initial z position of extended array (-2)
    integer, intent(in) :: mzpp3
    !# final z position of extended array (mzp+3)
    integer, intent(in) :: mxi
    !# initial x position of extended array (-2)
    integer, intent(in) :: mxpp3
    !# final x position of extended array (mxp+3)
    integer, intent(in) :: myi
    !# initial y position of extended array (-2)
    integer, intent(in) :: mypp3
    !# final y position of extended array (myp+3)
    real,intent(inout)    :: scr(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var
    real,intent(inout)    :: ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var rhou*U
    real,intent(inout)    :: vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var rhou*V
    real,intent(inout)    :: wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var rhou*W

    integer :: n

    ! direction X -
    do n=1,3
       scr(:,       1-n,:)=scr(:,       1,:)
       scr(:,mxp-1+n+is,:)=scr(:,mxp-1+is,:)
    enddo
    ufx_local(:,    0,:)=ufx_local(:,  1,:)
    ufx_local(:,mxp+1,:)=ufx_local(:,mxp,:)
    vfx_local(:,    0,:)=vfx_local(:,  1,:)
    wfx_local(:,    0,:)=wfx_local(:,  1,:)

    ! direction Y -
    do n = 1,3
       scr(:,:,       1-n)=scr(:,:,       1)
       scr(:,:,myp-1+n+js)=scr(:,:,myp-1+js)
    enddo
    ufx_local(:,:,    0)=ufx_local(:,:,  1)
    vfx_local(:,:,    0)=vfx_local(:,:,  1)
    vfx_local(:,:,myp+1)=vfx_local(:,:,myp)
    wfx_local(:,:,    0)=wfx_local(:,:,  1)

    ! direction z-
    do n = 1,3
       scr(1-n       ,:,:)=scr(1       ,:,:)
       scr(mzp-1+n+ks,:,:)=scr(mzp-1+ks,:,:)
    enddo
    ufx_local(0    ,:,:)=ufx_local(1,  :,:)
    vfx_local(0    ,:,:)=vfx_local(1,  :,:)
    wfx_local(0    ,:,:)=wfx_local(1  ,:,:)
    wfx_local(mzp+1,:,:)=wfx_local(mzp,:,:)
  end subroutine expandBorder




  
  subroutine commHalo(scr,ufx_local,vfx_local,wfx_local, &
       mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
       mxp,myp,mzp,myNum,nmachs,nodei0,nodej0,vname)
    !# Subroutine to communicated the extend halo
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine perform the communication among processors
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(in) :: mxp
    !# points in x dir
    integer, intent(in) :: myp
    !# points in y dir
    integer, intent(in) :: mzp
    !# points in z dir
    integer, intent(in) :: mzi
    !# initial z position of extended array (-2)
    integer, intent(in) :: mzpp3
    !# final z position of extended array (mzp+3)
    integer, intent(in) :: mxi
    !# initial x position of extended array (-2)
    integer, intent(in) :: mxpp3
    !# final x position of extended array (mxp+3)
    integer, intent(in) :: myi
    !# initial y position of extended array (-2)
    integer, intent(in) :: mypp3
    !# final y position of extended array (myp+3)
    integer, intent(in) :: nmachs
    !# total of processors
    integer, intent(in) :: myNum
    !# local processor
    integer, intent(in) :: nodei0(nmachs,1)
    !# increment for i in this proc
    integer, intent(in) :: nodej0(nmachs,1)
    !# increment for j in this proc
    character(len=*),intent(in) :: vname
    !# name of var to be communicated
    real,intent(inout)    :: scr(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var
    real,intent(inout)    :: ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var rhou*U
    real,intent(inout)    :: vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var rhou*V
    real,intent(inout)    :: wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    !# Oversided (extended) var rhou*W

    integer, parameter :: nOfVars=4
    !# total of vars to be communicated
    integer, parameter :: sizeOfHalo=6
    !# size of the halo on communication(both)
    integer :: iSend
    !# Count of sends
    integer :: iRecv
    !# Count of receives
    integer :: iRecS
    !# Count of waits
    integer :: iCnt
    !# Count of position for message array
    integer :: inext
    !# Nexto first position of array
    integer :: SizeMess
    !# Total size of message
    integer :: reqRecv(nMess(myNum))
    !# Request for MPI message receive
    integer :: reqSend(nMess(myNum))
    !# Request for MPI message send
    integer :: i
    !# loop count
    integer :: j
    !# loop count
    integer :: k
    !# loop count
    integer :: recNum
    !# returned position of receive wait messages
    character :: cn
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(commHalo)**"
    character(len=8) :: str(10)

    if(nmachs==1) return ! Only for parallel comm

    !Send request for get messages from neighbour
    do iRecv = 1, nMess(myNum)
       ! Size of complete message with nOfVars + halo
       sizeMess=nOfVars*size(scr,1)*&
            (receive(mynum,iRecv)%yEnd-receive(mynum,iRecv)%yBeg+1)*&
            (receive(mynum,iRecv)%xEnd-receive(mynum,iRecv)%xBeg+1)
       allocate(receive(mynum,iRecv)%dataMess(sizeMess))
       call parf_get_noblock_real(receive(mynum,iRecv)%dataMess,sizeMess, &
            receive(mynum,iRecv)%from-1,receive(mynum,iRecv)%tag,reqRecv(iRecv))

       if (dumpLocal) then
          write(str(1),"(i8)") iRecv
          write(str(2),"(i8)") sizeMess
          write(str(3),"(i8)") receive(mynum,iRecv)%from-1
          write(str(4),"(i8)") receive(mynum,iRecv)%tag
          call MsgDump(h//" dispatch recv #"//trim(adjustl(str(1)))//&
               " of section of variables "//trim(adjustl(vname))//&
               ", ufx_local, vfx_local, wfx_local"//&
               " with "//trim(adjustl(str(2)))//" reals"//&
               " from MPI proc "//trim(adjustl(str(3)))//&
               " with tag "//trim(adjustl(str(4))))
       end if
    end do

    ! Send messages for all neighbour
    do iSend = 1, nMess(myNum)

       ! Size of complete message with nOfVars + halo
       sizeMess=nOfVars*size(scr,1)*&
            (send(mynum,iSend)%yEnd-send(mynum,iSend)%yBeg+1)*&
            (send(mynum,iSend)%xEnd-send(mynum,iSend)%xBeg+1)
       allocate(send(mynum,iSend)%dataMess(sizeMess))

       !Before send fills the %dataMess with nOfVars in sequence)
       iCnt = 1
       !CDIR NOVECT
       do j = send(mynum,iSend)%yBeg,send(mynum,iSend)%yEnd
          do i = send(mynum,iSend)%xBeg,send(mynum,iSend)%xEnd
             do k=-2,mzp+3
                send(mynum,iSend)%dataMess(iCnt) = scr(k,i,j)      ; iCnt=iCnt+1
                send(mynum,iSend)%dataMess(iCnt) = ufx_local(k,i,j); iCnt=iCnt+1
                send(mynum,iSend)%dataMess(iCnt) = vfx_local(k,i,j); iCnt=iCnt+1
                send(mynum,iSend)%dataMess(iCnt) = wfx_local(k,i,j); iCnt=iCnt+1
             end do
          end do
       end do

       !Send all data filled above to neighbour
       call parf_send_noblock_real(send(mynum,iSend)%dataMess,sizeMess,send(mynum,iSend)%dest-1, &
            send(mynum,iSend)%tag,reqSend(iSend))
       if (dumpLocal) then
          write(str(1),"(i8)") iSend
          write(str(2),"(i8)") sizeMess
          write(str(3),"(i8)") send(mynum,iSend)%dest-1
          write(str(4),"(i8)") send(mynum,iSend)%tag
          write(str(5),"(i8)") send(mynum,iSend)%yBeg+send(mynum,iSend)%j0
          write(str(6),"(i8)") send(mynum,iSend)%yEnd+send(mynum,iSend)%j0
          write(str(7),"(i8)") send(mynum,iSend)%xBeg+send(mynum,iSend)%i0
          write(str(8),"(i8)") send(mynum,iSend)%xEnd+send(mynum,iSend)%i0
          write(str(9),"(i8)") -2
          write(str(10),"(i8)") mzp+3
          call MsgDump(h//" dispatch send #"//trim(adjustl(str(1)))//&
               " of "//trim(adjustl(vname))//&
               "("//trim(adjustl(str(9)))//":"//trim(adjustl(str(10)))//","//&
               trim(adjustl(str(7)))//":"//trim(adjustl(str(8)))//","//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")"//&
               ", ufx_local, vfx_local, wfx_local"//&
               " with "//trim(adjustl(str(2)))//" reals"//&
               " to MPI proc "//trim(adjustl(str(3)))//&
               " with tag "//trim(adjustl(str(4))))
       end if
    end do
    
    !Waiting (assync) and receive of neighbour messages
    do iRecV = 1, nMess(myNum)

       !Receive the message recNum from request
       call parf_wait_any_nostatus(nMess(myNum),reqRecv,recNum)

       ! Filling back the vars with message %dataMess from specific recnum
       iCnt = 1
       !CDIR NOVECT
       do j = receive(mynum,recNum)%yBeg,receive(mynum,recNum)%yEnd
          do i = receive(mynum,recNum)%xBeg,receive(mynum,recNum)%xEnd
             do k=-2,mzp+3
                scr(k,i,j)      =receive(mynum,recNum)%dataMess(iCnt)
                iCnt=iCnt+1
                ufx_local(k,i,j)=receive(mynum,recNum)%dataMess(iCnt)
                iCnt=iCnt+1
                vfx_local(k,i,j)=receive(mynum,recNum)%dataMess(iCnt)
                iCnt=iCnt+1 
                wfx_local(k,i,j)=receive(mynum,recNum)%dataMess(iCnt)
                iCnt=iCnt+1
             end do
          end do
       end do
       if (dumpLocal) then
          write(str(1),"(i8)") recNum
          write(str(2),"(i8)") sizeMess
          write(str(3),"(i8)") receive(mynum,recNum)%from-1
          write(str(4),"(i8)") receive(mynum,recNum)%tag
          write(str(5),"(i8)") receive(mynum,recNum)%yBeg+receive(mynum,recNum)%j0
          write(str(6),"(i8)") receive(mynum,recNum)%yEnd+receive(mynum,recNum)%j0
          write(str(7),"(i8)") receive(mynum,recNum)%xBeg+receive(mynum,recNum)%i0
          write(str(8),"(i8)") receive(mynum,recNum)%xEnd+receive(mynum,recNum)%i0
          write(str(9),"(i8)") -2
          write(str(10),"(i8)") mzp+3
          call MsgDump(h//" recv #"//trim(adjustl(str(1)))//&
               " of "//trim(adjustl(vname))//&
               "("//trim(adjustl(str(9)))//":"//trim(adjustl(str(10)))//","//&
               trim(adjustl(str(7)))//":"//trim(adjustl(str(8)))//","//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")"//&
               ", ufx_local, vfx_local, wfx_local"//&
               " with "//trim(adjustl(str(2)))//" reals"//&
               " from MPI proc "//trim(adjustl(str(3)))//&
               " with tag "//trim(adjustl(str(4))))
       end if
    end do

    !Waiting for conclusion of all sends
    do iRecS=1,nMess(myNum)
       call parf_wait_any_nostatus(nMess(myNum),reqSend,recNum)
    enddo

    !Deallocating the local buffers data
    do iRecV = 1, nMess(myNum)
       deallocate(receive(mynum,iRecv)%dataMess)
       deallocate(send(mynum,iRecv)%dataMess)
    enddo
  end subroutine commHalo




  
  subroutine countMess(nmachs,ngrid,ixb,ixe,iyb,iye,nMess)
    !# Count the total of message send/receive for each processor
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine perform a count of total of message send/receive for each processor
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(in) :: nMachs
    !# total of processors
    integer, intent(in) :: nGrid
    !# total of grids
    integer, intent(in) :: ixb(nmachs,ngrid)
    !# First position in x dir
    integer, intent(in) :: ixe(nmachs,ngrid)
    !# Last position in x dir
    integer, intent(in) :: iyb(nmachs,ngrid)
    !# First position in y dir
    integer, intent(in) :: iye(nmachs,ngrid)
    !# Last position in x dir
    integer, intent(out) :: nMess(nmachs)
    !# Total of messages per processor

    integer :: p1
    !# Local count processor 1
    integer :: p2
    !# Local count processor 2

    nMess=0
    ! Run over all processors P1. For each P1 run over P2 processors.
    ! if last position+1 of P1 is equal the first position of P2 then
    ! they are neighbourhood. In this case increment nMess for each one
    ! In x direction:
    do p1=1,nmachs
       do p2=p1+1,nmachs
          if(ixe(p1,ngrid)+1==ixb(p2,ngrid)) then
             if (iye(p1,ngrid)<=iye(p2,ngrid) .and. &
                  iyb(p1,ngrid)<=iye(p2,ngrid)) then
                nMess(p1)=nMess(p1)+1
                nMess(p2)=nMess(p2)+1
             endif
          endif
       enddo
    enddo
    !In Y direction:
    do p1=1,nmachs
       do p2=p1+1,nmachs
          if(iye(p1,ngrid)+1==iyb(p2,ngrid)) then
             if(ixe(p1,ngrid)>=ixb(p2,ngrid) .and. &
                  ixb(p1,ngrid)<=ixe(p2,ngrid)) then
                nMess(p1)=nMess(p1)+1
                nMess(p2)=nMess(p2)+1
             endif
          endif
       enddo
    enddo
  end subroutine countMess

  


  
  subroutine fillSendReceive(nmachs,ngrid,ixb,ixe,iyb,iye,nMess,nodei0,nodej0, &
       nnxp,nnyp)
    !# Fill the send receive types with communication information
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine put all information needed by MPI comunication in types send & receive
    !# in order to be used by communication routine.
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(in) :: nMachs
    !# total of processors
    integer, intent(in) :: nGrid
    !# total of grids
    integer, intent(in) :: ixb(nmachs,ngrid)
    !# First position in x dir
    integer, intent(in) :: ixe(nmachs,ngrid)
    !# Last position in x dir
    integer, intent(in) :: iyb(nmachs,ngrid)
    !# First position in y dir
    integer, intent(in) :: iye(nmachs,ngrid)
    !# Last position in x dir
    integer, intent(in) :: nodei0(nmachs,ngrid)
    !# increment for i in this proc
    integer, intent(in) :: nodej0(nmachs,ngrid)
    !# increment for j in this proc
    integer, intent(in) :: nnxp(ngrid)
    !# total of points in x dir
    integer, intent(in) :: nnyp(ngrid)
    !# total of points in y dir
    integer, intent(out) :: nMess(nmachs)
    !# Total of messages per processor

    integer :: p1
    !# Local count processor 1
    integer :: p2
    !# Local count processor 2
    integer :: itag
    !# tag count to control communication
    integer :: isize
    !# total size of the message
    integer :: ib
    !# first position of message in dir
    integer :: ie
    !# last position of message in dir

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(fillSendReceive)**"
    character(len=8) :: str

    if (dumpLocal) then
       call MsgDump(h//" computes msg passing send, send_acou, receive, receive_acou")
    end if
    
    itag=0
    nMess=0

    !Fill the border array with presence of neighbour in 4 directions
    border=.false.
    do p1=1,nmachs
       if(nodei0(p1,ngrid)==0)      border(p1, west)=.true.
       if(ixe(p1,ngrid)+1==nnxp(1)) border(p1, east)=.true.
       if(nodej0(p1,ngrid)==0)      border(p1,north)=.true.
       if(iye(p1,ngrid)+1==nnyp(1)) border(p1,south)=.true.
    enddo
    ! Run over all processors P1. For each P1 run over P2 processors.
    ! if last position+1 of P1 is equal the first position of P2 then
    ! they are neighbourhood. In this case increment nMess for each one
    ! In x direction:
    do p1=1,nmachs
       do p2=p1+1,nmachs
          ! X directions
          if(ixe(p1,ngrid)+1==ixb(p2,ngrid)) then
             if (iye(p1,ngrid)>=iye(p2,ngrid) .and. &
                  iyb(p1,ngrid)>=iyb(p2,ngrid)) then
                nMess(p1)=nMess(p1)+1

                !3 positions halo from p1 to p2
                itag=itag+1
                send(p1,nMess(p1))%dest     =p2 !Destination processor
                send_acou(p1,nmess(p1))%dest=p2
                send(p1,nMess(p1))%tag      =itag+1000 ! Tag used to send
                send_acou(p1,nMess(p1))%tag =itag+5000
                ! The lenght is  the size of common border 3 times (3 positions)
                send(p1,nMess(p1))%length   =(iye(p1,ngrid)-iyb(p1,ngrid)+3)*3

                ! The begin and end positions are adjust for local node position
                send(p1,nMess(p1))%i0       =nodei0(p1,ngrid)
                send(p1,nMess(p1))%xBeg     =ixe(p1,ngrid)-2-nodei0(p1,ngrid)
                send(p1,nMess(p1))%xEnd     =ixe(p1,ngrid)  -nodei0(p1,ngrid)

                send(p1,nMess(p1))%j0       =nodej0(p1,ngrid)
                if(border(p1, north)) then
                   send(p1,nMess(p1))%yBeg     =iyb(p1,ngrid)  -nodej0(p1,ngrid)-1
                else
                   send(p1,nMess(p1))%yBeg     =iyb(p1,ngrid)  -nodej0(p1,ngrid)
                endif
                if(border(p1, south)) then
                   send(p1,nMess(p1))%yEnd     =iye(p1,ngrid)  -nodej0(p1,ngrid)+1
                else
                   send(p1,nMess(p1))%yEnd     =iye(p1,ngrid)  -nodej0(p1,ngrid)
                endif
                !
                send_acou(p1,nMess(p1))%i0  =nodei0(p1,ngrid)
                send_acou(p1,nMess(p1))%xBeg=ixe(p1,ngrid)  -nodei0(p1,ngrid)
                send_acou(p1,nMess(p1))%xEnd=ixe(p1,ngrid)  -nodei0(p1,ngrid)

                send_acou(p1,nMess(p1))%j0  =nodej0(p1,ngrid)
                if(border(p1, north)) then
                   send_acou(p1,nMess(p1))%yBeg=iyb(p1,ngrid)  -nodej0(p1,ngrid)-1
                else
                   send_acou(p1,nMess(p1))%yBeg=iyb(p1,ngrid)  -nodej0(p1,ngrid)
                endif
                if(border(p1, south)) then
                   send_acou(p1,nMess(p1))%yEnd=iye(p1,ngrid)  -nodej0(p1,ngrid)+1
                else
                   send_acou(p1,nMess(p1))%yEnd=iye(p1,ngrid)  -nodej0(p1,ngrid)
                endif

                ! The algorithm used above is the same for the others halo below
                itag=itag+1
                receive(p1,nMess(p1))%from       =p2
                receive_acou(p1,nMess(p1))%from  =p2
                receive(p1,nMess(p1))%tag        =itag+2000
                receive_acou(p1,nMess(p1))%tag   =itag+6000
                receive(p1,nMess(p1))%length     =(iye(p1,ngrid)-iyb(p1,ngrid)+3)*3

                receive(p1,nMess(p1))%i0         =nodei0(p1,ngrid)
                receive(p1,nMess(p1))%xBeg       =ixe(p1,ngrid)+1-nodei0(p1,ngrid)
                receive(p1,nMess(p1))%xEnd       =ixe(p1,ngrid)+3-nodei0(p1,ngrid)

                receive(p1,nMess(p1))%j0         =nodej0(p1,ngrid)
                if(border(p1, north)) then
                   receive(p1,nMess(p1))%yBeg       =iyb(p1,ngrid)  -nodej0(p1,ngrid)-1
                else
                   receive(p1,nMess(p1))%yBeg       =iyb(p1,ngrid)  -nodej0(p1,ngrid)
                endif
                if(border(p1, south)) then
                   receive(p1,nMess(p1))%yEnd       =iye(p1,ngrid)  -nodej0(p1,ngrid)+1
                else
                   receive(p1,nMess(p1))%yEnd       =iye(p1,ngrid)  -nodej0(p1,ngrid)
                endif

                receive_acou(p1,nMess(p1))%i0    =nodei0(p1,ngrid)
                receive_acou(p1,nMess(p1))%xBeg  =ixe(p1,ngrid)+1-nodei0(p1,ngrid)
                receive_acou(p1,nMess(p1))%xEnd  =ixe(p1,ngrid)+1-nodei0(p1,ngrid)

                receive_acou(p1,nMess(p1))%j0    =nodej0(p1,ngrid)
                if(border(p1, north)) then
                   receive_acou(p1,nMess(p1))%yBeg  =iyb(p1,ngrid)  -nodej0(p1,ngrid)-1
                else
                   receive_acou(p1,nMess(p1))%yBeg  =iyb(p1,ngrid)  -nodej0(p1,ngrid)
                endif
                if(border(p1, south)) then
                   receive_acou(p1,nMess(p1))%yEnd  =iye(p1,ngrid)  -nodej0(p1,ngrid)+1
                else
                   receive_acou(p1,nMess(p1))%yEnd  =iye(p1,ngrid)  -nodej0(p1,ngrid)
                endif
                nMess(p2)=nMess(p2)+1

                !Fronteira de 3 indo de p2 para p1
                send(p2,nMess(p2))%dest       =p1
                send_acou(p2,nMess(p2))%dest  =p1
                send(p2,nMess(p2))%tag        =receive(p1,nMess(p1))%tag
                send_acou(p2,nMess(p2))%tag   =receive_acou(p1,nMess(p1))%tag
                send(p2,nMess(p2))%length   =(iye(p1,ngrid)-iyb(p1,ngrid)+3)*3

                send(p2,nMess(p2))%i0       =nodei0(p2,ngrid)
                send(p2,nMess(p2))%xBeg     =ixb(p2,ngrid)  -nodei0(p2,ngrid)
                send(p2,nMess(p2))%xEnd     =ixb(p2,ngrid)+2-nodei0(p2,ngrid)

                send(p2,nMess(p2))%j0       =nodej0(p2,ngrid)
                if(border(p1, north)) then
                   send(p2,nMess(p2))%yBeg     =iyb(p2,ngrid)  -nodej0(p2,ngrid)-1
                else
                   send(p2,nMess(p2))%yBeg     =iyb(p2,ngrid)  -nodej0(p2,ngrid)
                endif
                if(border(p1, south)) then
                   send(p2,nMess(p2))%yEnd     =iye(p2,ngrid)  -nodej0(p2,ngrid)+1
                else
                   send(p2,nMess(p2))%yEnd     =iye(p2,ngrid)  -nodej0(p2,ngrid)
                endif

                send_acou(p2,nMess(p2))%i0       =nodei0(p2,ngrid)
                send_acou(p2,nMess(p2))%xBeg     =ixb(p2,ngrid)  -nodei0(p2,ngrid)
                send_acou(p2,nMess(p2))%xEnd     =ixb(p2,ngrid)  -nodei0(p2,ngrid)

                send_acou(p2,nMess(p2))%j0       =nodej0(p2,ngrid)
                if(border(p1, north)) then
                   send_acou(p2,nMess(p2))%yBeg     =iyb(p2,ngrid)  -nodej0(p2,ngrid)-1
                else
                   send_acou(p2,nMess(p2))%yBeg     =iyb(p2,ngrid)  -nodej0(p2,ngrid)
                endif
                if(border(p1, south)) then
                   send_Acou(p2,nMess(p2))%yEnd     =iye(p2,ngrid)  -nodej0(p2,ngrid)+1
                else
                   send_Acou(p2,nMess(p2))%yEnd     =iye(p2,ngrid)  -nodej0(p2,ngrid)
                endif

                !Fronteira de 3 vindo de p1 para p2
                receive(p2,nMess(p2))%from       =p1
                receive_acou(p2,nMess(p2))%from  =p1
                receive(p2,nMess(p2))%tag   =send(p1,nMess(p1))%tag
                receive_acou(p2,nMess(p2))%tag   =send_acou(p1,nMess(p1))%tag
                receive(p2,nMess(p2))%length=(iye(p1,ngrid)-iyb(p1,ngrid)+3)*3

                receive(p2,nMess(p2))%i0    =nodei0(p2,ngrid)
                receive(p2,nMess(p2))%xBeg  =ixb(p2,ngrid)-3-nodei0(p2,ngrid)
                receive(p2,nMess(p2))%xEnd  =ixb(p2,ngrid)-1-nodei0(p2,ngrid)

                receive(p2,nMess(p2))%j0    =nodej0(p2,ngrid)
                if(border(p1, north)) then
                   receive(p2,nMess(p2))%yBeg  =iyb(p2,ngrid)  -nodej0(p2,ngrid)-1
                else
                   receive(p2,nMess(p2))%yBeg  =iyb(p2,ngrid)  -nodej0(p2,ngrid)
                endif
                if(border(p1, south)) then
                   receive(p2,nMess(p2))%yEnd  =iye(p2,ngrid)-  nodej0(p2,ngrid)+1
                else
                   receive(p2,nMess(p2))%yEnd  =iye(p2,ngrid)-  nodej0(p2,ngrid)
                endif

                receive_acou(p2,nMess(p2))%i0    =nodei0(p2,ngrid)
                receive_acou(p2,nMess(p2))%xBeg  =ixb(p2,ngrid)-1-nodei0(p2,ngrid)
                receive_acou(p2,nMess(p2))%xEnd  =ixb(p2,ngrid)-1-nodei0(p2,ngrid)

                receive_acou(p2,nMess(p2))%j0    =nodej0(p2,ngrid)
                if(border(p1, north)) then
                   receive_acou(p2,nMess(p2))%yBeg  =iyb(p2,ngrid)  -nodej0(p2,ngrid)-1
                else
                   receive_acou(p2,nMess(p2))%yBeg  =iyb(p2,ngrid)  -nodej0(p2,ngrid)
                endif
                if(border(p1, south)) then
                   receive_acou(p2,nMess(p2))%yEnd  =iye(p2,ngrid)  -nodej0(p2,ngrid)+1
                else
                   receive_acou(p2,nMess(p2))%yEnd  =iye(p2,ngrid)  -nodej0(p2,ngrid)
                endif
             endif
          endif
       enddo
    enddo
    do p1=1,nmachs
       do p2=p1+1,nmachs
          ! Y directions
          if(iye(p1,ngrid)+1==iyb(p2,ngrid)) then
             if(ixe(p1,ngrid)>=ixb(p2,ngrid) .and. ixb(p1,ngrid)<=ixe(p2,ngrid)) then
                nMess(p1)=nMess(p1)+1
                ib=max(ixb(p2,ngrid),ixb(p1,ngrid))
                ie=min(ixe(p2,ngrid),ixe(p1,ngrid))
                isize=ie-ib+3
                !Fronteira de 3 indo de p1 para p2
                itag=itag+1
                send(p1,nMess(p1))%dest       =p2
                send_acou(p1,nMess(p1))%dest       =p2
                send(p1,nMess(p1))%tag      =itag+3000
                send_acou(p1,nMess(p1))%tag =itag+7000
                send(p1,nMess(p1))%length   =isize*3

                send(p1,nMess(p1))%i0       =nodei0(p1,ngrid)
                if(border(p1, west)) then
                   send(p1,nMess(p1))%xBeg     =ib             -nodei0(p1,ngrid)-1
                else
                   send(p1,nMess(p1))%xBeg     =ib             -nodei0(p1,ngrid)
                endif
                if(border(p1, east)) then
                   send(p1,nMess(p1))%xEnd     =ie             -nodei0(p1,ngrid)+1
                else
                   send(p1,nMess(p1))%xEnd     =ie             -nodei0(p1,ngrid)
                endif

                send(p1,nMess(p1))%j0       =nodej0(p1,ngrid)
                send(p1,nMess(p1))%yBeg     =iye(p1,ngrid)-2-nodej0(p1,ngrid)
                send(p1,nMess(p1))%yEnd     =iye(p1,ngrid)  -nodej0(p1,ngrid)

                send_acou(p1,nMess(p1))%i0  =nodei0(p1,ngrid)
                if(border(p1, west)) then
                   send_acou(p1,nMess(p1))%xBeg     =ib-nodei0(p1,ngrid)-1
                else
                   send_acou(p1,nMess(p1))%xBeg     =ib-nodei0(p1,ngrid)
                endif
                if(border(p1, east)) then
                   send_acou(p1,nMess(p1))%xEnd     =ie  -nodei0(p1,ngrid)+1
                else
                   send_acou(p1,nMess(p1))%xEnd     =ie  -nodei0(p1,ngrid)
                endif

                send_acou(p1,nMess(p1))%j0       =nodej0(p1,ngrid)
                send_acou(p1,nMess(p1))%yBeg     =iye(p1,ngrid) -nodej0(p1,ngrid)
                send_acou(p1,nMess(p1))%yEnd     =iye(p1,ngrid) -nodej0(p1,ngrid)

                !Fronteira de 3 vindo de p2 para p1
                itag=itag+1
                receive(p1,nMess(p1))%from  =p2
                receive_acou(p1,nMess(p1))%from  =p2
                receive(p1,nMess(p1))%tag   =itag+4000
                receive_acou(p1,nMess(p1))%tag   =itag+8000
                receive(p1,nMess(p1))%length=isize*3

                receive(p1,nMess(p1))%i0 = nodei0(p1,ngrid)
                if(border(p1, west)) then
                   receive(p1,nMess(p1))%xBeg  =ib             -nodei0(p1,ngrid)-1
                else
                   receive(p1,nMess(p1))%xBeg  =ib             -nodei0(p1,ngrid)
                endif
                if(border(p1, east)) then
                   receive(p1,nMess(p1))%xEnd  =ie             -nodei0(p1,ngrid)+1
                else
                   receive(p1,nMess(p1))%xEnd  =ie             -nodei0(p1,ngrid)
                endif
                !

                receive(p1,nMess(p1))%j0 = nodej0(p1,ngrid)
                receive(p1,nMess(p1))%yBeg  =iye(p1,ngrid)+1-nodej0(p1,ngrid)
                receive(p1,nMess(p1))%yEnd  =iye(p1,ngrid)+3-nodej0(p1,ngrid)

                receive_acou(p1,nMess(p1))%i0 = nodei0(p1,ngrid)
                if(border(p1, west)) then
                   receive_acou(p1,nMess(p1))%xBeg  =ib             -nodei0(p1,ngrid)-1
                else
                   receive_acou(p1,nMess(p1))%xBeg  =ib             -nodei0(p1,ngrid)
                endif
                if(border(p1, east)) then
                   receive_acou(p1,nMess(p1))%xEnd  =ie             -nodei0(p1,ngrid)+1
                else
                   receive_acou(p1,nMess(p1))%xEnd  =ie             -nodei0(p1,ngrid)
                endif

                receive_acou(p1,nMess(p1))%j0 = nodej0(p1,ngrid)
                receive_acou(p1,nMess(p1))%yBeg  =iye(p1,ngrid)  -nodej0(p1,ngrid)+1
                receive_acou(p1,nMess(p1))%yEnd  =iye(p1,ngrid)  -nodej0(p1,ngrid)+1
                nMess(p2)=nMess(p2)+1

                !Fronteira de 3 indo de p2 para p1
                send(p2,nMess(p2))%dest       =p1
                send_acou(p2,nMess(p2))%dest  =p1
                send(p2,nMess(p2))%tag      =receive(p1,nMess(p1))%tag
                send_acou(p2,nMess(p2))%tag =receive_Acou(p1,nMess(p1))%tag
                send(p2,nMess(p2))%length   =isize*3

                send(p2,nMess(p2))%i0 = nodei0(p2,ngrid)
                if(border(p1, west)) then
                   send(p2,nMess(p2))%xBeg     =ib             -nodei0(p2,ngrid)-1
                else
                   send(p2,nMess(p2))%xBeg     =ib             -nodei0(p2,ngrid)
                endif
                if(border(p1, east)) then
                   send(p2,nMess(p2))%xEnd     =ie             -nodei0(p2,ngrid)+1
                else
                   send(p2,nMess(p2))%xEnd     =ie             -nodei0(p2,ngrid)
                endif

                send(p2,nMess(p2))%j0 = nodej0(p2,ngrid)
                send(p2,nMess(p2))%yBeg     =iyb(p2,ngrid)  -nodej0(p2,ngrid)
                send(p2,nMess(p2))%yEnd     =iyb(p2,ngrid)+2-nodej0(p2,ngrid)

                send_acou(p2,nMess(p2))%i0 = nodei0(p2,ngrid)
                if(border(p1, west)) then
                   send_acou(p2,nMess(p2))%xBeg     =ib             -nodei0(p2,ngrid)-1
                else
                   send_acou(p2,nMess(p2))%xBeg     =ib             -nodei0(p2,ngrid)
                endif
                if(border(p1, east)) then
                   send_acou(p2,nMess(p2))%xEnd     =ie             -nodei0(p2,ngrid)+1
                else
                   send_acou(p2,nMess(p2))%xEnd     =ie             -nodei0(p2,ngrid)
                endif

                send_acou(p2,nMess(p2))%j0 = nodej0(p2,ngrid)
                send_acou(p2,nMess(p2))%yBeg     =iyb(p2,ngrid)  -nodej0(p2,ngrid)
                send_acou(p2,nMess(p2))%yEnd     =iyb(p2,ngrid)  -nodej0(p2,ngrid)

                !Fronteira de 3 vindo de p1 para p2
                receive(p2,nMess(p2))%from  =p1
                receive_acou(p2,nMess(p2))%from  =p1
                receive(p2,nMess(p2))%tag   =send(p1,nMess(p1))%tag
                receive_acou(p2,nMess(p2))%tag   =send_acou(p1,nMess(p1))%tag
                receive(p2,nMess(p2))%length=isize*3

                receive(p2,nMess(p2))%i0 = nodei0(p2,ngrid)
                if(border(p1, west)) then
                   receive(p2,nMess(p2))%xBeg  =ib             -nodei0(p2,ngrid)-1
                else
                   receive(p2,nMess(p2))%xBeg  =ib             -nodei0(p2,ngrid)
                endif
                if(border(p1, east)) then
                   receive(p2,nMess(p2))%xEnd  =ie             -nodei0(p2,ngrid)+1
                else
                   receive(p2,nMess(p2))%xEnd  =ie             -nodei0(p2,ngrid)
                endif

                receive(p2,nMess(p2))%j0 = nodej0(p2,ngrid)
                receive(p2,nMess(p2))%yBeg  =iyb(p2,ngrid)-3-nodej0(p2,ngrid)
                receive(p2,nMess(p2))%yEnd  =iyb(p2,ngrid)-1-nodej0(p2,ngrid)

                receive_acou(p2,nMess(p2))%i0 = nodei0(p2,ngrid)
                if(border(p1, west)) then
                   receive_acou(p2,nMess(p2))%xBeg  =ib             -nodei0(p2,ngrid)-1
                else
                   receive_acou(p2,nMess(p2))%xBeg  =ib             -nodei0(p2,ngrid)
                endif
                if(border(p1, east)) then
                   receive_acou(p2,nMess(p2))%xEnd  =ie             -nodei0(p2,ngrid)+1
                else
                   receive_acou(p2,nMess(p2))%xEnd  =ie             -nodei0(p2,ngrid)
                endif

                receive_acou(p2,nMess(p2))%j0 = nodej0(p2,ngrid)
                receive_acou(p2,nMess(p2))%yBeg  =iyb(p2,ngrid)  -nodej0(p2,ngrid)-1
                receive_acou(p2,nMess(p2))%yEnd  =iyb(p2,ngrid)  -nodej0(p2,ngrid)-1
             endif
          endif
       enddo
    enddo
    receive_acou(:,:)%length=0
    send_acou(:,:)%length=0
  end subroutine fillSendReceive




  
  subroutine SendandSetAll(master_num,mynum,nmachs)
    !# Broadcast send receive types to all processors
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine broadcast all information needed by MPI comunication in types send & receive
    !# to all processors in order to be used by communication routine in each one.
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(in) :: master_num
    !# The master processor number
    integer, intent(in) :: mynum
    !# my local (this) processor
    integer, intent(in) :: nmachs
    !# total os processors

    integer :: p1
    !# local count of processors
    integer :: n
    !# count of messages in each processor
    integer :: l
    !# Auxiliary to communicate logical vars


    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(SendandSetAll)**"
    character(len=8) :: str

    if (dumpLocal) then
       call MsgDump(h//" broadcast msg passing send, send_acou, receive, receive_acou")
    end if

    ! In order to allocate the total of messages on each processor
    ! send the count of messages to all
    call Broadcast(nMess, master_num, "nMess")
    !If I am a 'slave' allocate send and receive local to number of messages
    if(mynum/=1) then
       allocate(send(nmachs,maxval(nMess)),receive(nmachs,maxval(nMess)))
       allocate(send_acou(nmachs,maxval(nMess)),receive_acou(nmachs,maxval(nMess)))
    endif
    ! Broadcast type message info to all
    do p1=1,nmachs
       do n=1,nMess(p1)
          call Broadcast(send(p1,n)%dest    ,master_num,'sto    ')
          call Broadcast(send(p1,n)%tag   ,master_num,'stag   ')
          call Broadcast(send(p1,n)%length,master_num,'slength')
          call Broadcast(send(p1,n)%i0  ,master_num,'i0  ')
          call Broadcast(send(p1,n)%xBeg  ,master_num,'sxBeg  ')
          call Broadcast(send(p1,n)%xEnd  ,master_num,'sxEnd  ')
          call Broadcast(send(p1,n)%j0  ,master_num,'j0  ')
          call Broadcast(send(p1,n)%yBeg  ,master_num,'syBeg  ')
          call Broadcast(send(p1,n)%yEnd  ,master_num,'syEnd  ')

          call Broadcast(receive(p1,n)%from  ,master_num,'rfrom    ')
          call Broadcast(receive(p1,n)%tag   ,master_num,'rtag   ')
          call Broadcast(receive(p1,n)%length,master_num,'rlength')
          call Broadcast(receive(p1,n)%i0    ,master_num,'i0  ')
          call Broadcast(receive(p1,n)%xBeg  ,master_num,'rxBeg  ')
          call Broadcast(receive(p1,n)%xEnd  ,master_num,'rxEnd  ')
          call Broadcast(receive(p1,n)%j0    ,master_num,'j0  ')
          call Broadcast(receive(p1,n)%yBeg  ,master_num,'ryBeg  ')
          call Broadcast(receive(p1,n)%yEnd  ,master_num,'ryEnd  ')

          call Broadcast(send_acou(p1,n)%dest  ,master_num,'Asto    ')
          call Broadcast(send_acou(p1,n)%tag   ,master_num,'Astag   ')
          call Broadcast(send_acou(p1,n)%length,master_num,'Aslength')
          call Broadcast(send_acou(p1,n)%i0    ,master_num,'i0  ')
          call Broadcast(send_acou(p1,n)%xBeg  ,master_num,'AsxBeg  ')
          call Broadcast(send_acou(p1,n)%xEnd  ,master_num,'AsxEnd  ')
          call Broadcast(send_acou(p1,n)%j0    ,master_num,'j0  ')
          call Broadcast(send_acou(p1,n)%yBeg  ,master_num,'AsyBeg  ')
          call Broadcast(send_acou(p1,n)%yEnd  ,master_num,'AsyEnd  ')

          call Broadcast(receive_acou(p1,n)%from  ,master_num,'Arfrom    ')
          call Broadcast(receive_acou(p1,n)%tag   ,master_num,'Artag   ')
          call Broadcast(receive_acou(p1,n)%length,master_num,'Arlength')
          call Broadcast(receive_acou(p1,n)%i0    ,master_num,'Ari0  ')
          call Broadcast(receive_acou(p1,n)%xBeg  ,master_num,'ArxBeg  ')
          call Broadcast(receive_acou(p1,n)%xEnd  ,master_num,'ArxEnd  ')
          call Broadcast(receive_acou(p1,n)%j0    ,master_num,'Arj0  ')
          call Broadcast(receive_acou(p1,n)%yBeg  ,master_num,'AryBeg  ')
          call Broadcast(receive_acou(p1,n)%yEnd  ,master_num,'AryEnd  ')

       enddo
       do n=1,4
          if(mynum==1 .and. border(p1,n)) then
             l=1
          else
             l=0
          end if
          call Broadcast(l,master_num,'bdir')
          if(l==1) then
             border(p1,n)=.true.
          else
             border(p1,n)=.false.
          end if
       enddo
    end do


  end subroutine SendandSetAll




  
  subroutine dumpComm(mxp,myp,mzp,myNum,nmachs,nodei0,nodej0)
    !# Dump communication initialization in warm up of model
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine dump the initialization of comm info if the number of processor
    !# is less than 8. Is usefull for debug of comm problems.
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    include "constants.h"
    integer, intent(in) :: mxp
    !# points in x dir
    integer, intent(in) :: myp
    !# points in y dir
    integer, intent(in) :: mzp
    !# points in z dir
    integer, intent(in) :: myNum
    !# Number of local processor
    integer, intent(in) :: nMachs
    !# Total of processors
    integer, intent(in) :: nodei0(nmachs,1)
    !# increment for i in this proc
    integer, intent(in) :: nodej0(nmachs,1)
    !# increment for j in this proc

    integer :: p1,p2,n

    if(nMachs>7) return !Because the size the dump works only for <7 procs
    open(unit=22,file='brams.log',position='append',action='write')
    if(mynum==1) then !Just one proc do the dump
       !Sends
       write (unit=22,fmt='(A2, A, A6, 1X, 6(A3,1X),A2,1X,4(A1,1X))') '','              ',  &
            'tag','to','lng','xbg','xnd','ybg','ynd','#m','W','E','N','S'
       do p1=1,nmachs
          do n=1,nMess(p1)
             write (unit=22,fmt='(I2, A, I6,1X,6(I3,1X),I2.2,1X,4(L1,1X))') p1,' ..send to... ', &
                  send(p1,n)%tag,send(p1,n)%dest,send(p1,n)%length, &
                  send(p1,n)%xBeg,send(p1,n)%xEnd,send(p1,n)%yBeg, &
                  send(p1,n)%yEnd,n,border(p1,1),border(p1,2),border(p1,3),border(p1,4)
          end do
       end do
       !Receives
       write (unit=22,fmt='(A2, A, A6, 1X, 6(A3,1X),A2)') '','              ',  &
            'tag','frm','lng','xbg','xnd','ybg','ynd','#m'
       do p1=1,nmachs
          !Receives
          do n=1,nMess(p1)
             write (unit=22,fmt='(I2, A, I6, 1X, 6(I3,1X),I2.2)') p1,' receive from ', &
                  receive(p1,n)%tag,receive(p1,n)%from,receive(p1,n)%length, &
                  receive(p1,n)%xBeg,receive(p1,n)%xEnd,receive(p1,n)%yBeg, &
                  receive(p1,n)%yEnd,n
          end do
       end do
       write(unit=22,fmt='(A)') '------ Comm for acoustic ----'
       !Sends Acoustic
       write (unit=22,fmt='(A2, A, A6, 1X, 6(A3,1X),A2,1X,4(A1,1X))') '','              ',  &
            'tag','to','lng','xbg','xnd','ybg','ynd','#m','W','E','N','S'
       do p1=1,nmachs
          do n=1,nMess(p1)
             write (unit=22,fmt='(I2, A, I6,1X,6(I3,1X),I2.2,1X,4(L1,1X))') p1,' ..send to... ', &
                  send_acou(p1,n)%tag,send_acou(p1,n)%dest,send_acou(p1,n)%length, &
                  send_acou(p1,n)%xBeg,send_acou(p1,n)%xEnd,send_acou(p1,n)%yBeg, &
                  send_acou(p1,n)%yEnd,n,border(p1,1),border(p1,2),border(p1,3),border(p1,4)
          end do
       end do
       !Receives
       write (unit=22,fmt='(A2, A, A6, 1X, 6(A3,1X),A2)') '','              ',  &
            'tag','frm','lng','xbg','xnd','ybg','ynd','#m'
       do p1=1,nmachs
          !Receives Acoustic
          do n=1,nMess(p1)
             write (unit=22,fmt='(I2, A, I6, 1X, 6(I3,1X),I2.2)') p1,' receive from ', &
                  receive_acou(p1,n)%tag,receive_acou(p1,n)%from,receive_acou(p1,n)%length, &
                  receive_acou(p1,n)%xBeg,receive_acou(p1,n)%xEnd,receive_acou(p1,n)%yBeg, &
                  receive_acou(p1,n)%yEnd,n
          end do
       end do
    endif
    close(unit=22)
  end subroutine dumpComm




  
  subroutine dumpMess(cpos,fdump,sizeMess,destination,tag,req,nm,vname)
    !# Dump a communication message - only for debug proposal
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine perform a dump in communication message - Used by develop's team
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(in) :: sizeMess,destination,tag,req,fdump,nm
    character(len=*),intent(in) :: cpos,vname
    character(len=16) :: csiz
    character(len=256)  :: wform

    wform='(A,A,A,I2.2,A,I4,A,I4,A,Z8,A,I8,A,F8.1)'
    if(trim(cpos)=='get') then
       write(fdump,fmt=wform) 'I asked for var ',trim(vname),', message #',nm, ', with tag ',tag,', from proc ',destination, &
            ', using request ',req,', size= ',sizemess,' at time= ',time
    elseif(trim(cpos)=='snd') then
       write(fdump,fmt=wform) 'I sent the  var ',trim(vname),', message #',nm, ', with tag ',tag,',  to  proc ',destination, &
            ', using request ',req,', size= ',sizemess,' at time= ',time
    elseif(trim(cpos)=='wit') then
       write(fdump,fmt=wform) 'I received  var ',trim(vname),', message #',nm, ', with tag ',tag,', from proc ',destination, &
            ', using request ',req,', size= ',sizemess,' at time= ',time
    else
       write(fdump,fmt='(A)') 'Option '//cpos//' invalid! ("get", "snd" or "wit" are valid)'
    endif
    call flush(fdump)
  end subroutine dumpMess





  subroutine commHaloAcou(scr,mxp,myp,mzp,mynum,vname)
    !# Subroutine to communicated the the normal halo for 1 var
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This routine perform the communication among processors like the original Brams
    !# communication. The communication is only 1 halo and used for cases that the original comms
    !# not works.
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2018Jan
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !# 1. Must be transported for original communication
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(in) :: mxp
    !# points in x dir
    integer, intent(in) :: myp
    !# points in y dir
    integer, intent(in) :: mzp
    !# points in z dir
    integer, intent(in) :: myNum
    !# local processor
    character(len=*),intent(in) :: vname
    !# name of var to be communicated
    real,intent(inout) :: scr(mzp,mxp,myp)
    !# var to be communicated

    integer, parameter :: nOfVars=4
    !# total of vars to be communicated
    integer, parameter :: sizeOfHalo=6
    !# size of the halo on communication(both)
    integer :: iSend
    !# Count of sends
    integer :: iRecv
    !# Count of receives
    integer :: iRecS
    !# Count of waits
    integer :: iCnt
    !# Count of position for message array
    integer :: inext
    !# Nexto first position of array
    integer :: SizeMess
    !# Total size of message
    integer :: reqRecv_acou(nMess(myNum))
    !# Request for MPI message receive
    integer :: reqsend_acou(nMess(myNum))
    !# Request for MPI message send
    integer :: i
    !# loop count
    integer :: j
    !# loop count
    integer :: k
    !# loop count
    integer :: recNum
    !# returned position of receive wait messages
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(commHaloAcou)**"
    character(len=8) :: str(10)

    !Send request for get messages from neighbour
    do iRecv = 1, nMess(myNum)
       !
       ! Size of complete message with nOfVars + halo
       sizeMess=size(scr,1)*&
            (receive_acou(mynum,iRecv)%yEnd-receive_acou(mynum,iRecv)%yBeg+1)*&
            (receive_acou(mynum,iRecv)%xEnd-receive_acou(mynum,iRecv)%xBeg+1)
       allocate(receive_acou(mynum,iRecv)%dataMess(sizeMess))
       !
       call parf_get_noblock_real(receive_acou(mynum,iRecv)%dataMess,sizeMess, &
            receive_acou(mynum,iRecv)%from-1,receive_acou(mynum,iRecv)%tag,reqRecv_acou(iRecv))

       if (dumpLocal) then
          write(str(1),"(i8)") iRecv
          write(str(2),"(i8)") sizeMess
          write(str(3),"(i8)") receive_acou(mynum,iRecv)%from-1
          write(str(4),"(i8)") receive_acou(mynum,iRecv)%tag
          call MsgDump(h//" dispatch recv #"//trim(adjustl(str(1)))//&
               " of section of variable "//trim(adjustl(vname))//&
               " with "//trim(adjustl(str(2)))//" reals"//&
               " from MPI proc "//trim(adjustl(str(3)))//&
               " with tag "//trim(adjustl(str(4))))
       end if
    end do

    ! Send messages for all neighbour
    do iSend = 1, nMess(myNum)

       ! Size of complete message with nOfVars + halo
       sizeMess=size(scr,1)*&
            (send_acou(mynum,iSend)%yEnd-send_acou(mynum,iSend)%yBeg+1)*&
            (send_acou(mynum,iSend)%xEnd-send_acou(mynum,iSend)%xBeg+1)
       allocate(send_acou(mynum,iSend)%dataMess(sizeMess))

       !Before send it fills the %dataMess with nOfVars in sequence)
       iCnt = 1
       !CDIR NOVECT
       do j = send_acou(mynum,iSend)%yBeg,send_acou(mynum,iSend)%yEnd
          do i = send_acou(mynum,iSend)%xBeg,send_acou(mynum,iSend)%xEnd
             do k=1,mzp
                send_acou(mynum,iSend)%dataMess(iCnt) = scr(k,i,j); iCnt=iCnt+1
             end do
          end do
       end do

       !Send all data filled above to neighbour
       call parf_send_noblock_real(send_acou(mynum,iSend)%dataMess,sizeMess,&
            send_acou(mynum,iSend)%dest-1, &
            send_acou(mynum,iSend)%tag,reqsend_acou(iSend))

       if (dumpLocal) then
          write(str(1),"(i8)") iSend
          write(str(2),"(i8)") sizeMess
          write(str(3),"(i8)") send_acou(mynum,iSend)%dest-1
          write(str(4),"(i8)") send_acou(mynum,iSend)%tag
          write(str(5),"(i8)") send_acou(mynum,iSend)%yBeg+send_acou(mynum,iSend)%j0
          write(str(6),"(i8)") send_acou(mynum,iSend)%yEnd+send_acou(mynum,iSend)%j0
          write(str(7),"(i8)") send_acou(mynum,iSend)%xBeg+send_acou(mynum,iSend)%i0
          write(str(8),"(i8)") send_acou(mynum,iSend)%xEnd+send_acou(mynum,iSend)%i0
          write(str(9),"(i8)") 1
          write(str(10),"(i8)") mzp
          call MsgDump(h//" dispatch send #"//trim(adjustl(str(1)))//&
               " of "//trim(adjustl(vname))//&
               "("//trim(adjustl(str(9)))//":"//trim(adjustl(str(10)))//","//&
               trim(adjustl(str(7)))//":"//trim(adjustl(str(8)))//","//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")"//&
               " with "//trim(adjustl(str(2)))//" reals"//&
               " to MPI proc "//trim(adjustl(str(3)))//&
               " with tag "//trim(adjustl(str(4))))
       end if
    end do


    !Waiting (assync) and receive_acou of neighbour messages
    do iRecV = 1, nMess(myNum)
       !receive_acou the message recNum from request
       call parf_wait_any_nostatus(nMess(myNum),reqRecv_acou,recNum)

       ! Filling back the vars with message %dataMess from specific recnum
       iCnt = 1
       !CDIR NOVECT
       do j = receive_acou(mynum,recNum)%yBeg,receive_acou(mynum,recNum)%yEnd
          do i = receive_acou(mynum,recNum)%xBeg,receive_acou(mynum,recNum)%xEnd
             do k=1,mzp
                scr(k,i,j)      =receive_acou(mynum,recNum)%dataMess(iCnt); iCnt=iCnt+1
             end do
          end do
       end do
       if (dumpLocal) then
          write(str(1),"(i8)") recNum
          write(str(2),"(i8)") sizeMess
          write(str(3),"(i8)") receive_acou(mynum,recNum)%from-1
          write(str(4),"(i8)") receive_acou(mynum,recNum)%tag
          write(str(5),"(i8)") receive_acou(mynum,recNum)%yBeg+receive_acou(mynum,recNum)%j0
          write(str(6),"(i8)") receive_acou(mynum,recNum)%yEnd+receive_acou(mynum,recNum)%j0
          write(str(7),"(i8)") receive_acou(mynum,recNum)%xBeg+receive_acou(mynum,recNum)%i0
          write(str(8),"(i8)") receive_acou(mynum,recNum)%xEnd+receive_acou(mynum,recNum)%i0
          write(str(9),"(i8)") 1
          write(str(10),"(i8)") mzp
          call MsgDump(h//" recv #"//trim(adjustl(str(1)))//&
               " of "//trim(adjustl(vname))//&
               "("//trim(adjustl(str(9)))//":"//trim(adjustl(str(10)))//","//&
               trim(adjustl(str(7)))//":"//trim(adjustl(str(8)))//","//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")"//&
               " with "//trim(adjustl(str(2)))//" reals"//&
               " from MPI proc "//trim(adjustl(str(3)))//&
               " with tag "//trim(adjustl(str(4))))
       end if
    end do

    !Waiting for conclusion of all sends
    do iRecS=1,nMess(myNum)
       call parf_wait_any_nostatus(nMess(myNum),reqSend_acou,recNum)
    enddo

    !Deallocating the local buffers data
    do iRecV = 1, nMess(myNum)
       deallocate(receive_acou(mynum,iRecv)%dataMess)
       deallocate(send_acou(mynum,iRecv)%dataMess)
    enddo
  end subroutine commHaloAcou
end module modComm
