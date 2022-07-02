!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModRnode

  use ParLib, only: &
       parf_pack_max_size

  use grid_dims, only: &
       maxmach

  use mem_grid, only: &
       ngrids, &
       nnxp, &
       nnyp, &
       nnzp, &
       ibnd, &
       jbnd

  use node_mod, only: &
       f_ndmd_size,   &
       nbuff_feed,    &
       node_buffs_st, &
       nmachs,        &
       ipaths,        &
       node_buffs,    &
       mynum, &
       ia,   &
       ia1,  &
       ia_1,  &
       ibcon,  &
       iz,  &
       iz1,  &
       iz_1,  &
       izu,  &
       ja,  &
       ja1,  &
       ja_1,  &
       jz,  &
       jz1,  &
       jz_1,  &
       jzv,  &
       mxp,  &
       myp

  use ModVarTables, only: &
       num_var,         &
       vtab_r

  use mem_leaf, only : &
       ISFCL ! For SiB


  implicit none

  private

  public :: node_index
  public :: InitFields
contains

  !     ****************************************************************

  subroutine node_index()

    ia_1=max(ia-1,1)
    ia1=ia+1
    iz_1=iz-1
    iz1=min(iz+1,mxp)

    izu=iz
    if(iand(ibcon,2).ne.0) izu=iz-1

    if(myp.gt.1) then
       ja_1=max(ja-1,1)
       ja1=ja+1
       jz_1=jz-1
       jz1=min(jz+1,myp)

       jzv=jz
       if(iand(ibcon,8).ne.0) jzv=jz-1

    else
       print*,'Trying to do 2-dimensional run ??????'
       stop 'no parallel 2d'
       ja_1=1
       ja1=1
       jz_1=1
       jz1=1
       jzv=1
    endif
  end subroutine node_index

  ! ---------------------------------------------------------------------------

  subroutine InitFields(oneScalarTabSize, init)
    integer, intent(in) :: init
    integer, intent(in) :: oneScalarTabSize
    ! Local variables:

    include "constants.h"
    integer(i8) :: silly_i8
    integer :: ng, nm, itype, i1, j1, i2, j2, memf, npvar, nv


    !     Can we use existing memory for the nesting communication buffers?
    !       If not, allocate new buffers or compute buffer sizes.

    !       Check feedback buffer.

    ! itype=6 !Changed for reproducibility - Saulo Barros
    itype=7

    nbuff_feed=0
    do ng=1,ngrids
       do nm=1,nmachs
          i1=ipaths(1,itype,ng,nm)
          i2=ipaths(2,itype,ng,nm)
          j1=ipaths(3,itype,ng,nm)
          j2=ipaths(4,itype,ng,nm)
          memf=(i2-i1+1)*(j2-j1+1)*(nnzp(ng))  &
               *(4+oneScalarTabSize)
          nbuff_feed=max(nbuff_feed,memf)
       enddo
    enddo

    !____________________________________________
    !
    !    Allocate long time step send and receive buffers

    call parf_pack_max_size(1_i8, silly_i8)
    f_ndmd_size=silly_i8

    if(init /= 1) then
       do nm=1,nmachs
          if (associated(node_buffs(nm)%lbc_send_buff) ) &
               deallocate(node_buffs(nm)%lbc_send_buff)
          if (associated(node_buffs(nm)%lbc_recv_buff) ) &
               deallocate(node_buffs(nm)%lbc_recv_buff)
          if (associated(node_buffs_st(nm)%lbc_send_buff) ) &
               deallocate(node_buffs_st(nm)%lbc_send_buff)
          if (associated(node_buffs_st(nm)%lbc_recv_buff) ) &
               deallocate(node_buffs_st(nm)%lbc_recv_buff)
       enddo
    endif

    do nm=1,nmachs
       node_buffs(nm)%nsend = max(node_buffs(nm)%nsend ,nbuff_feed)
       node_buffs(nm)%nrecv = max(node_buffs(nm)%nrecv ,nbuff_feed)
       node_buffs_st(nm)%nsend = node_buffs(nm)%nsend
       node_buffs_st(nm)%nrecv = node_buffs(nm)%nrecv

       if (node_buffs(nm)%nsend > 0) &
            allocate(node_buffs(nm)%lbc_send_buff(node_buffs(nm)%nsend))
       if (node_buffs(nm)%nrecv > 0) &
            allocate(node_buffs(nm)%lbc_recv_buff(node_buffs(nm)%nrecv))
       if (node_buffs_st(nm)%nsend > 0) &
            allocate(node_buffs_st(nm)%lbc_send_buff(node_buffs_st(nm)%nsend*f_ndmd_size))
       if (node_buffs_st(nm)%nrecv > 0) &
            allocate(node_buffs_st(nm)%lbc_recv_buff(node_buffs_st(nm)%nrecv*f_ndmd_size))
    enddo

    !  Find number of lbc variables to be communicated.
    npvar=0
    do nv = 1,num_var(1)
       if(vtab_r(nv,1)%impt1 == 1 ) then
          npvar=npvar+1
       endif
    enddo
  end subroutine InitFields
end module ModRnode
