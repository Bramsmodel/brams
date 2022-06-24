!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModMPassFull

  implicit none

  include "mpif.h"

  private
  public :: mk_2_buff
  public :: mk_3_buff
  public :: mk_4_buff
  public :: ex_2_buff
  public :: ex_2p_buff
  public :: ex_3_buff
  public :: ex_4_buff
  public :: ex_buff_carma

contains

  !---------------------------------------------------------------------------

  !For Carma - CATT
  subroutine mk_buff_carma(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
    real, intent(in) :: a(:,:,:)
    real, intent(inout) :: b(:,:,:)

    b(1:m1,1:m2,1:m3)=a(i1:i2,j1:j2,1:n3)

    return
  end subroutine mk_buff_carma

  !---------------------------------------------------------------------------

  subroutine mk_2_buff(a,b,n1,n2,m1,m2,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,m1,m2,i1,i2,j1,j2
    real, intent(in) :: a(:,:)
    real, intent(inout) :: b(:,:)

    b(1:m1,1:m2)=a(i1:i2,j1:j2)

    return
  end subroutine mk_2_buff

  !---------------------------------------------------------------------------

  subroutine mk_2p_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
    real, intent(in) :: a(:,:,:)
    real, intent(inout) :: b(:,:,:)

    b(1:m1,1:m2,1:m3)=a(i1:i2,j1:j2,1:n3)

    return
  end subroutine mk_2p_buff

  !---------------------------------------------------------------------------

  subroutine mk_3_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
    real, intent(in) :: a(:,:,:)
    real, intent(inout) :: b(:,:,:)

    b(1:m1,1:m2,1:m3)=a(1:n1,i1:i2,j1:j2)

    return
  end subroutine mk_3_buff

  !---------------------------------------------------------------------------

  subroutine mk_4_buff(a,b,n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2
    real, intent(in) :: a(:,:,:,:)
    real, intent(inout) :: b(:,:,:,:)

    b(1:m1,1:m2,1:m3,1:m4)=a(1:n1,i1:i2,j1:j2,1:n4)

    return
  end subroutine mk_4_buff

  !---------------------------------------------------------------------------

  subroutine ex_buff_carma(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
    real, intent(inout) :: a(n1,n2,n3)
    real, intent(in) :: b(m1,m2,m3)
!!$    real, intent(inout) :: a(:,:,:)
!!$    real, intent(in) :: b(:,:,:)

    a(i1+i0:i2+i0,j1+j0:j2+j0,1:n3) = b(i1:i2,j1:j2,1:m3)

    return
  end subroutine ex_buff_carma

  !---------------------------------------------------------------------------

  subroutine ex_2_buff(a,b,n1,n2,m1,m2,i0,j0,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,m1,m2,i0,j0,i1,i2,j1,j2
    real, intent(inout) :: a(n1,n2)
    real, intent(in) :: b(m1,m2)
!!$    real, intent(inout) :: a(:,:)
!!$    real, intent(in) :: b(:,:)

    a(i1+i0:i2+i0,j1+j0:j2+j0) = b(i1:i2,j1:j2)

    return
  end subroutine ex_2_buff

  !---------------------------------------------------------------------------

  subroutine ex_3_buff(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
    real, intent(inout) :: a(n1,n2,n3)
    real, intent(in) :: b(m1,m2,m3)
!!$    real, intent(inout) :: a(:,:,:)
!!$    real, intent(in) :: b(:,:,:)

    a(1:n1,i1+i0:i2+i0,j1+j0:j2+j0) = b(1:m1,i1:i2,j1:j2)

    return
  end subroutine ex_3_buff

  !---------------------------------------------------------------------------

  subroutine ex_4_buff(a,b,n1,n2,n3,n4,m1,m2,m3,m4,i0,j0,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,n3,n4,m1,m2,m3,m4,i0,j0,i1,i2,j1,j2
    real, intent(inout) :: a(n1,n2,n3,n4)
    real, intent(in) :: b(m1,m2,m3,m4)
!!$    real, intent(inout) :: a(:,:,:,:)
!!$    real, intent(in) :: b(:,:,:,:)

    a(1:n1,i1+i0:i2+i0,j1+j0:j2+j0,1:n4) = b(1:m1,i1:i2,j1:j2,1:m4)

    return
  end subroutine ex_4_buff

  !---------------------------------------------------------------------------

  subroutine ex_2p_buff(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
    implicit none
    integer, intent(in) :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
    real, intent(inout) :: a(n1,n2,n3)
    real, intent(in) :: b(m1,m2,m3)
!!$    real, intent(inout) :: a(:,:,:)
!!$    real, intent(in) :: b(:,:,:)

    a(i1+i0:i2+i0,j1+j0:j2+j0,1:n3) = b(i1:i2,j1:j2,1:m3)

    return
  end subroutine ex_2p_buff

!!$  !---------------------------------------------------------------------------
!!$
!!$  subroutine GatherAllChunks (nmachs, master_num, toGather, sizeToGather, &
!!$       gathered, sizeGathered, eachSizeGathered, displacement)
!!$    implicit none
!!$    integer, intent(in) :: nmachs
!!$    integer, intent(in) :: master_num
!!$    integer, intent(in) :: sizeToGather
!!$    integer, intent(in) :: sizeGathered
!!$    real, intent(in) :: toGather(sizeToGather)
!!$    real, intent(out) :: gathered(sizeGathered)
!!$    integer, intent(in) :: eachSizeGathered(nmachs)
!!$    integer, intent(in) :: displacement(nmachs)
!!$    integer :: ierr
!!$    character(len=8) :: c0, c1, c2, c3
!!$    character(len=*), parameter :: h="**(GatherAllChunks)**" 
!!$    logical, parameter :: dumpLocal=.false.
!!$
!!$    if (dumpLocal) then
!!$       write(c0,"(i8)") sizeToGather
!!$       write(c1,"(i8)") sizeGathered
!!$       write(c2,"(i8)") eachSizeGathered(1)
!!$       write(c3,"(i8)") displacement(1)
!!$       write(*,"(a)") h//" sizeToGather="//trim(adjustl(c0))//&
!!$            "; sizeGathered="//trim(adjustl(c1))//&
!!$            "; eachSizeGathered(1)="//trim(adjustl(c2))//&
!!$            "; displacement(1)="//trim(adjustl(c3))
!!$       call flush(6)
!!$       write(c0,"(e8.1)") sum(toGather)
!!$       write(*,"(a)") h//" sum(toGather)="//trim(adjustl(c0))
!!$       call flush(6)
!!$       write(c0,"(e8.1)") sum(gathered)
!!$       write(*,"(a)") h//" sum(gathered)="//trim(adjustl(c0))
!!$       call flush(6)
!!$    end if
!!$    ! gather a field
!!$
!!$    call MPI_GATHERV(toGather, sizeToGather, MPI_REAL, &
!!$         gathered, eachSizeGathered, displacement, MPI_REAL, &
!!$         master_num, MPI_COMM_WORLD, ierr)
!!$
!!$    if (ierr /= MPI_SUCCESS) THEN
!!$       write(c0,"(i8)") ierr
!!$       call fatal_error(h//" gatherv fails with ierr="//trim(adjustl(c0)))
!!$    else if (dumpLocal) then
!!$       write(*,"(a)") h//" done"
!!$       call flush(6)
!!$    end if
!!$  end subroutine GatherAllChunks
!!$
!!$  !---------------------------------------------------------------------------
!!$
!!$  subroutine StoreGathered (idim_type, nmachs, sizeGathered, gathered, &
!!$       nnxp, nnyp, nnzp, nzs, nzg, npatch, nwave, &
!!$       nodemxp, nodemyp, nodemzp, nodeia, nodeiz, nodeja, nodejz, &
!!$       nodei0, nodej0, il1, ir2, jb1, jt2, displacement, eachSizeGathered, &
!!$       sizeStored, stored)
!!$    implicit none
!!$    integer, intent(in) :: idim_type
!!$    integer, intent(in) :: nmachs
!!$    integer, intent(in) :: sizeGathered
!!$    real,    intent(in) :: gathered(sizeGathered)
!!$    integer, intent(in) :: nnxp
!!$    integer, intent(in) :: nnyp
!!$    integer, intent(in) :: nnzp
!!$    integer, intent(in) :: nzs
!!$    integer, intent(in) :: nzg
!!$    integer, intent(in) :: npatch
!!$    integer, intent(in) :: nwave
!!$    integer, intent(in) :: nodemxp(nmachs)
!!$    integer, intent(in) :: nodemyp(nmachs)
!!$    integer, intent(in) :: nodemzp(nmachs)
!!$    integer, intent(in) :: nodeia(nmachs)
!!$    integer, intent(in) :: nodeiz(nmachs)
!!$    integer, intent(in) :: nodeja(nmachs)
!!$    integer, intent(in) :: nodejz(nmachs)
!!$    integer, intent(in) :: nodei0(nmachs)
!!$    integer, intent(in) :: nodej0(nmachs)
!!$    integer, intent(in) :: il1(nmachs)
!!$    integer, intent(in) :: ir2(nmachs)
!!$    integer, intent(in) :: jb1(nmachs)
!!$    integer, intent(in) :: jt2(nmachs)
!!$    integer, intent(in) :: displacement(nmachs)
!!$    integer, intent(in) :: eachSizeGathered(nmachs)
!!$    integer, intent(in) :: sizeStored
!!$    real,    intent(out) :: stored(sizeStored)
!!$
!!$    integer :: proc
!!$    character(len=8) :: c0
!!$    character(len=*), parameter :: h="**(StoreGathered)**"
!!$
!!$    ! unpack gathered field, removing unnecessary ghost zones and
!!$    ! placing entries at correct field positions
!!$
!!$    select case (idim_type)
!!$    case (2)
!!$       do proc = 1, nmachs
!!$          if (eachSizeGathered(proc)/=0) then
!!$             call ex_1_buff(stored, gathered(displacement(proc)+1), &
!!$                  nnxp, nnyp, nodemxp(proc), nodemyp(proc), &
!!$                  nodei0(proc), nodej0(proc), &
!!$                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
!!$          end if
!!$       end do
!!$    case (3)
!!$       do proc = 1, nmachs
!!$          if (eachSizeGathered(proc)/=0) then
!!$             call ex_3_buff(stored, gathered(displacement(proc)+1), &
!!$                  nnzp, nnxp, nnyp, nnzp, nodemxp(proc), nodemyp(proc), &
!!$                  nodei0(proc), nodej0(proc), &
!!$                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
!!$          end if
!!$       end do
!!$    case (4)
!!$       do proc = 1, nmachs
!!$          if (eachSizeGathered(proc)/=0) then
!!$             call ex_4_buff(stored, gathered(displacement(proc)+1), &
!!$                  nzg, nnxp, nnyp, npatch, nzg, nodemxp(proc), nodemyp(proc), npatch, &
!!$                  nodei0(proc), nodej0(proc), &
!!$                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
!!$          end if
!!$       end do
!!$    case (5)
!!$       do proc = 1, nmachs
!!$          if (eachSizeGathered(proc)/=0) then
!!$             call ex_4_buff(stored, gathered(displacement(proc)+1), &
!!$                  nzs, nnxp, nnyp, npatch, nzs, nodemxp(proc), nodemyp(proc), npatch, &
!!$                  nodei0(proc), nodej0(proc), &
!!$                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
!!$          end if
!!$       end do
!!$    case (6)
!!$       do proc = 1, nmachs
!!$          if (eachSizeGathered(proc)/=0) then
!!$             call ex_2p_buff(stored, gathered(displacement(proc)+1), &
!!$                  nnxp, nnyp, npatch, nodemxp(proc), nodemyp(proc), npatch, &
!!$                  nodei0(proc), nodej0(proc), &
!!$                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
!!$          end if
!!$       end do
!!$    case (7)
!!$       do proc = 1, nmachs
!!$          if (eachSizeGathered(proc)/=0) then
!!$             call ex_buff_carma(stored, gathered(displacement(proc)+1), &
!!$                  nnxp, nnyp, nwave, nodemxp(proc), nodemyp(proc), nwave, &
!!$                  nodei0(proc), nodej0(proc), &
!!$                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
!!$          end if
!!$       end do
!!$    case default
!!$       write(c0,"(i8)") idim_type
!!$       call fatal_error(h//" invoked with unknown idim_type ("//&
!!$            trim(adjustl(c0))//")")
!!$    end select
!!$  end subroutine StoreGathered
!!$
!!$  !---------------------------------------------------------------------------
!!$
!!$  subroutine SizesDisplacements (idim_type, nmachs, &
!!$       nnxp, nnyp, nnzp, nzs, nzg, npatch, nwave, &
!!$       nodemxp, nodemyp, nodemzp, nodeia, nodeiz, nodeja, nodejz, nodeibcon, &
!!$       il1, ir2, jb1, jt2, eachSizeGathered, displacement)
!!$    implicit none
!!$    integer, intent(in) :: idim_type 
!!$    integer, intent(in) :: nmachs
!!$    integer, intent(in) :: nnxp
!!$    integer, intent(in) :: nnyp
!!$    integer, intent(in) :: nnzp
!!$    integer, intent(in) :: nzs
!!$    integer, intent(in) :: nzg
!!$    integer, intent(in) :: npatch
!!$    integer, intent(in) :: nwave
!!$    integer, intent(in) :: nodemxp(nmachs)
!!$    integer, intent(in) :: nodemyp(nmachs)
!!$    integer, intent(in) :: nodemzp(nmachs)
!!$    integer, intent(in) :: nodeia(nmachs)
!!$    integer, intent(in) :: nodeiz(nmachs)
!!$    integer, intent(in) :: nodeja(nmachs)
!!$    integer, intent(in) :: nodejz(nmachs)
!!$    integer, intent(in) :: nodeibcon(nmachs)
!!$    integer, intent(out) :: il1(nmachs)
!!$    integer, intent(out) :: ir2(nmachs)
!!$    integer, intent(out) :: jb1(nmachs)
!!$    integer, intent(out) :: jt2(nmachs)
!!$    integer, intent(out) :: eachSizeGathered(nmachs)
!!$    integer, intent(out) :: displacement(nmachs)
!!$
!!$    integer :: proc
!!$    character(len=8) :: c0
!!$    character(len=*), parameter :: h="**(SizesDisplacements)**" 
!!$
!!$    ! field portion to be gathered at each process
!!$
!!$    select case (idim_type)
!!$    case (2)
!!$       do proc = 1, nmachs
!!$          eachSizeGathered(proc) = nodemxp(proc)*nodemyp(proc)
!!$       end do
!!$    case (3)
!!$       do proc = 1, nmachs
!!$          eachSizeGathered(proc) = nodemzp(proc)*nodemxp(proc)*nodemyp(proc)
!!$       end do
!!$    case (4)
!!$       do proc = 1, nmachs
!!$          eachSizeGathered(proc) = nzg*nodemxp(proc)*nodemyp(proc)*npatch
!!$       end do
!!$    case (5)
!!$       do proc = 1, nmachs
!!$          eachSizeGathered(proc) = nzs*nodemxp(proc)*nodemyp(proc)*npatch
!!$       end do
!!$    case (6)
!!$       do proc = 1, nmachs
!!$          eachSizeGathered(proc) = nodemxp(proc)*nodemyp(proc)*npatch
!!$       end do
!!$    case (7)
!!$       do proc = 1, nmachs
!!$          eachSizeGathered(proc) = nodemxp(proc)*nodemyp(proc)*nwave
!!$       end do
!!$    case default
!!$       write(c0,"(i8)") idim_type
!!$       call fatal_error(h//" invoked with unknown idim_type ("//&
!!$            trim(adjustl(c0))//")")
!!$    end select
!!$
!!$    ! where to place field portion at gathering result
!!$
!!$    displacement(1)=0
!!$    do proc = 2, nmachs
!!$       displacement(proc)=displacement(proc-1)+eachSizeGathered(proc-1)
!!$    end do
!!$
!!$    ! eliminate unnecessary ghost zones while unpacking gathered result
!!$
!!$    do proc = 1, nmachs
!!$       if (btest(nodeibcon(proc),0)) then
!!$          il1(proc) = nodeia(proc) - 1
!!$       else
!!$          il1(proc) = nodeia(proc)
!!$       end if
!!$       if (btest(nodeibcon(proc),1)) then
!!$          ir2(proc) = nodeiz(proc) + 1
!!$       else
!!$          ir2(proc) = nodeiz(proc)
!!$       end if
!!$       if (btest(nodeibcon(proc),2)) then
!!$          jb1(proc) = nodeja(proc) - 1
!!$       else
!!$          jb1(proc) = nodeja(proc)
!!$       end if
!!$       if (btest(nodeibcon(proc),3)) then
!!$          jt2(proc) = nodejz(proc) + 1
!!$       else
!!$          jt2(proc) = nodejz(proc)
!!$       end if
!!$    end do
!!$  end subroutine SizesDisplacements
!!$
!!$  !---------------------------------------------------------------------------
!!$
!!$  subroutine OneGatherStoreAllChunks(mchnum, mynum, nmachs, master_num, &
!!$       toGather, idim_type, eachSizeGathered, displacement, &
!!$       nnxp, nnyp, nnzp, nzs, nzg, npatch, nwave, &
!!$       nodemxp, nodemyp, nodemzp, nodeia, nodeiz, nodeja, nodejz, &
!!$       nodei0, nodej0, il1, ir2, jb1, jt2, &
!!$       sizeStored, stored)
!!$    implicit none
!!$    integer, intent(in) :: mchnum
!!$    integer, intent(in) :: mynum
!!$    integer, intent(in) :: nmachs
!!$    integer, intent(in) :: master_num
!!$    real,    intent(in) :: toGather(*)
!!$    integer, intent(in) :: idim_type
!!$    integer, intent(in) :: eachSizeGathered(nmachs)
!!$    integer, intent(in) :: displacement(nmachs)
!!$    integer, intent(in) :: nnxp
!!$    integer, intent(in) :: nnyp
!!$    integer, intent(in) :: nnzp
!!$    integer, intent(in) :: nzs
!!$    integer, intent(in) :: nzg
!!$    integer, intent(in) :: npatch
!!$    integer, intent(in) :: nwave
!!$    integer, intent(in) :: nodemxp(nmachs)
!!$    integer, intent(in) :: nodemyp(nmachs)
!!$    integer, intent(in) :: nodemzp(nmachs)
!!$    integer, intent(in) :: nodeia(nmachs)
!!$    integer, intent(in) :: nodeiz(nmachs)
!!$    integer, intent(in) :: nodeja(nmachs)
!!$    integer, intent(in) :: nodejz(nmachs)
!!$    integer, intent(in) :: il1(nmachs)
!!$    integer, intent(in) :: ir2(nmachs)
!!$    integer, intent(in) :: jb1(nmachs)
!!$    integer, intent(in) :: jt2(nmachs)
!!$    integer, intent(in) :: nodei0(nmachs)
!!$    integer, intent(in) :: nodej0(nmachs)
!!$    integer, intent(in) :: sizeStored
!!$    real,    intent(out) :: stored(sizeStored)
!!$
!!$    integer :: ierr
!!$    integer :: sizeGathered
!!$    real, allocatable :: gathered(:)
!!$    character(len=8) :: c0, c1
!!$    character(len=*), parameter :: h="**(OneGatherStoreAllChunks)**"
!!$    logical, parameter :: dumpLocal=.false.
!!$
!!$    ! size of gathered field (contains all ghost zones of all processes)
!!$
!!$    sizeGathered = displacement(nmachs)+eachSizeGathered(nmachs)
!!$    allocate(gathered(sizeGathered), stat=ierr)
!!$    if (ierr /= 0) then
!!$       write(c0,"(i8)") sizeGathered
!!$       write(c1,"(i8)") ierr
!!$       call fatal_error(h//" allocate gathered("//trim(adjustl(c0))//&
!!$            ") failed with stat="//trim(adjustl(c1)))
!!$    else if (dumpLocal) then
!!$       write(c0,"(i8)") sizeGathered
!!$       write(*,"(a)") h//" allocated gathered("//trim(adjustl(c0))//")"
!!$       call flush(6)
!!$    end if
!!$
!!$    ! gathered field at master_num
!!$
!!$    if (dumpLocal) then
!!$       write(c0,"(i8)") eachSizeGathered(mynum)
!!$       write(*,"(a)") h//" will gather with local size "//trim(adjustl(c0))
!!$       call flush(6)
!!$    end if
!!$    call GatherAllChunks (nmachs, master_num, toGather, eachSizeGathered(mynum), &
!!$         gathered(1), sizeGathered, eachSizeGathered, displacement)
!!$    if (dumpLocal) then
!!$       write(*,"(a)") h//" done gathering"
!!$    end if
!!$
!!$    ! master_num unpacks fields (removes unnecessary ghost zones and positions entries)
!!$
!!$    if (mchnum == master_num) then
!!$       if (dumpLocal) then
!!$          write(*,"(a)") h//" master will StoreGathered"
!!$       end if
!!$       call StoreGathered (idim_type, nmachs, sizeGathered, gathered, &
!!$            nnxp, nnyp, nnzp, nzs, nzg, npatch, nwave, &
!!$            nodemxp, nodemyp, nodemzp, nodeia, nodeiz, nodeja, nodejz, &
!!$            nodei0, nodej0, il1, ir2, jb1, jt2, displacement, eachSizeGathered, &
!!$            sizeStored, stored)
!!$       if (dumpLocal) then
!!$          write(*,"(a)") h//" done StoreGathered"
!!$       end if
!!$    end if
!!$
!!$    deallocate(gathered, stat=ierr)
!!$    if (ierr /= 0) then
!!$       write(c1,"(i8)") ierr
!!$       call fatal_error(h//" deallocate gathered failed with stat="//trim(adjustl(c1)))
!!$    end if
!!$  end subroutine OneGatherStoreAllChunks
!!$
!!$  !---------------------------------------------------------------------------
end module ModMPassFull
