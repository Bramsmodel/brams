!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_cuparm

  use ModNamelistFile, only: &
       NamelistFile

  use grid_dims, only: &
       maxgrds,        & ! INTENT(IN)
       maxfiles          ! INTENT(IN)

  implicit none

  type cuparm_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, contiguous :: thsrc(:,:,:)
     real, pointer, contiguous :: rtsrc(:,:,:)
     real, pointer, contiguous :: thsrcf(:,:,:)
     real, pointer, contiguous :: rtsrcf(:,:,:)
     real, pointer, contiguous :: thsrcp(:,:,:)
     real, pointer, contiguous :: rtsrcp(:,:,:)
     real, pointer, contiguous :: clsrc(:,:,:) !srf-cloud/ice source term

     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, contiguous :: aconpr(:,:)
     real, pointer, contiguous :: conprr(:,:)
     real, pointer, contiguous :: conprrp(:,:)
     real, pointer, contiguous :: conprrf(:,:)

  end type cuparm_vars

  type (cuparm_vars), pointer :: cuparm_g(:) => null()
  type (cuparm_vars), pointer :: cuparmm_g(:) => null()
  type (cuparm_vars), pointer :: cuparm_g_sh(:) => null()
  type (cuparm_vars), pointer :: cuparmm_g_sh(:) => null()

  include "files.h"

  integer, parameter :: maxcufiles = maxfiles
  integer, parameter :: maxcugrids = 10

  character(len=f_name_length) :: fnames_cu(maxcufiles)
  character(len=14)  :: itotdate_cu(maxcufiles)
  real :: cu_times(maxcufiles)

  integer :: ncufiles
  integer :: ncufl
  real :: cutime1
  real :: cutime2

contains

  logical function hasAconpr(oneNamelistFile, ng)
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    integer, intent(IN) :: ng

    hasAconpr = oneNamelistFile%nnqparm(ng)>= 1 .or. oneNamelistFile%if_cuinv == 1

  end function hasAconpr

  subroutine alloc_cuparm(oneNamelistFile, cuparm, n1, n2, n3, ng)
    ! Arguments:
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type (cuparm_vars), intent(INOUT) :: cuparm
    integer, intent(IN)               :: n1, n2, n3, ng
    ! Local Variables:
    character(len=*), parameter :: h="**(alloc_cuparm)**"
    integer :: ierr

    ! Allocate arrays based on options (if necessary)

    if( hasAconpr(oneNamelistFile, ng) )  then
       if(oneNamelistFile%nnqparm(ng) < 3) then !srf: for conv 3 and 4, special arrays
          !srf:  are allocated instead of these ones

          allocate (cuparm%thsrc(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%thsrc")
          cuparm%thsrc=0.0
          allocate (cuparm%rtsrc(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%rtsrc")
          cuparm%rtsrc=0.0
          allocate (cuparm%clsrc(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%clsrc")
          cuparm%clsrc=0.0 !srf-cloud/ice source term

       endif
       allocate (cuparm%aconpr(n2,n3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating cuparm%aconpr")
       cuparm%aconpr = 0.
       allocate (cuparm%conprr(n2,n3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating cuparm%conprr")
       cuparm%conprr = 0.  

       if (oneNamelistFile%if_cuinv == 1) then
          allocate (cuparm%thsrcp(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%thsrcp")
          cuparm%thsrcp=0.0
          allocate (cuparm%rtsrcp(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%rtsrcp")
          cuparm%rtsrcp=0.0
          allocate (cuparm%thsrcf(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%thsrcf")
          cuparm%thsrcf=0.0
          allocate (cuparm%rtsrcf(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%rtsrcf")
          cuparm%rtsrcf=0.0
          allocate (cuparm%conprrp(n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%conprp")
          cuparm%conprrp=0.0
          allocate (cuparm%conprrf(n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%conprrf")
          cuparm%conprrf=0.0
       endif
    endif

  end subroutine alloc_cuparm

  ! ----------------------------------------------------------------------
  subroutine alloc_cuparm_sh(cuparm,n1,n2,n3,ng)
    type (cuparm_vars) :: cuparm
    integer, intent(in) :: n1,n2,n3,ng
    character(len=*), parameter :: h="**(alloc_cuparm)**"
    integer :: ierr
    ! Allocate arrays for shallow cum feedback

    allocate (cuparm%thsrc(n1,n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating cuparm%thsrc")
    cuparm%thsrc=0.0
    allocate (cuparm%rtsrc(n1,n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating cuparm%rtsrc")
    cuparm%rtsrc=0.0

  end subroutine alloc_cuparm_sh

  ! ----------------------------------------------------------------------

  subroutine nullify_cuparm(cuparm)
    ! Arguments:
    type (cuparm_vars), intent(inout) :: cuparm

    if (associated(cuparm%thsrc))    nullify (cuparm%thsrc)
    if (associated(cuparm%rtsrc))    nullify (cuparm%rtsrc)
    if (associated(cuparm%clsrc))    nullify (cuparm%clsrc)
    if (associated(cuparm%thsrcp))    nullify (cuparm%thsrcp)
    if (associated(cuparm%rtsrcp))    nullify (cuparm%rtsrcp)
    if (associated(cuparm%thsrcf))    nullify (cuparm%thsrcf)
    if (associated(cuparm%rtsrcf))    nullify (cuparm%rtsrcf)
    if (associated(cuparm%aconpr))   nullify (cuparm%aconpr)
    if (associated(cuparm%conprr))   nullify (cuparm%conprr)
    if (associated(cuparm%conprrp))   nullify (cuparm%conprrp)
    if (associated(cuparm%conprrf))   nullify (cuparm%conprrf)

  end subroutine nullify_cuparm

  ! ----------------------------------------------------------------------

  subroutine dealloc_cuparm(cuparm)
    ! Arguments:
    type (cuparm_vars), pointer, intent(inout) :: cuparm(:)

    integer :: igrid
    
    if (associated(cuparm)) then
       do igrid=1,size(cuparm)
          if (associated(cuparm(igrid)%thsrc))    deallocate (cuparm(igrid)%thsrc)
          if (associated(cuparm(igrid)%rtsrc))    deallocate (cuparm(igrid)%rtsrc)
          if (associated(cuparm(igrid)%clsrc))    deallocate (cuparm(igrid)%clsrc)
          if (associated(cuparm(igrid)%thsrcp))    deallocate (cuparm(igrid)%thsrcp)
          if (associated(cuparm(igrid)%rtsrcp))    deallocate (cuparm(igrid)%rtsrcp)
          if (associated(cuparm(igrid)%thsrcf))    deallocate (cuparm(igrid)%thsrcf)
          if (associated(cuparm(igrid)%rtsrcf))    deallocate (cuparm(igrid)%rtsrcf)
          if (associated(cuparm(igrid)%aconpr))   deallocate (cuparm(igrid)%aconpr)
          if (associated(cuparm(igrid)%conprr))   deallocate (cuparm(igrid)%conprr)
          if (associated(cuparm(igrid)%conprrp))   deallocate (cuparm(igrid)%conprrp)
          if (associated(cuparm(igrid)%conprrf))   deallocate (cuparm(igrid)%conprrf)
       end do
       nullify(cuparm)
    end if

  end subroutine dealloc_cuparm

  ! ----------------------------------------------------------------------

  subroutine filltab_cuparm_sh(oneVarTable, oneVarTableSize, cuparm, cuparmm, &
       ng, imean)

    use iso_fortran_env, only: &
         int64

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    ! Arguments:
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type (cuparm_vars), pointer, intent(in) :: cuparm(:)
    type (cuparm_vars), pointer, intent(in) :: cuparmm(:)
    integer, intent(in) :: ng
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(filltab_cuparm_sh)**"

    if (.not. associated(cuparm)) then
       call fatal_error(h//" invoked with unasociated cuparm")
    else if (.not. associated(cuparmm)) then
       call fatal_error(h//" invoked with unasociated cuparmm")
    else if (.not. associated(oneVarTable)) then
       call fatal_error(h//" invoked with unasociated oneVarTable")
    end if

    ! Fill pointers to arrays into variable tables

    if (associated(cuparm(ng)%thsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%thsrc, 'THSRC_SH :3:hist:anal:mpti:mpt3', &
            cuparmm(ng)%thsrc, imean)
    end if

    if (associated(cuparm(ng)%rtsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%rtsrc, 'RTSRC_SH :3:hist:anal:mpti:mpt3', &
            cuparmm(ng)%rtsrc, imean)
    end if

  end subroutine filltab_cuparm_sh

  ! ----------------------------------------------------------------------

  subroutine filltab_cuparm(oneVarTable, oneVarTableSize, cuparm, cuparmm, &
       ng, imean)

    use iso_fortran_env, only: &
         int64

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    ! Arguments:
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type (cuparm_vars), pointer, intent(in) :: cuparm(:)
    type (cuparm_vars), pointer, intent(in) :: cuparmm(:)
    integer, intent(in) :: ng
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(filltab_cuparm)**"

    if (.not. associated(cuparm)) then
       call fatal_error(h//" invoked with unasociated cuparm")
    else if (.not. associated(cuparmm)) then
       call fatal_error(h//" invoked with unasociated cuparmm")
    else if (.not. associated(oneVarTable)) then
       call fatal_error(h//" invoked with unasociated oneVarTable")
    end if

    ! Fill pointers to arrays into variable tables

    if (associated(cuparm(ng)%thsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%thsrc, 'THSRC :3:hist:anal:mpti:mpt3', &
            cuparmm(ng)%thsrc, imean)
    end if

    if (associated(cuparm(ng)%rtsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%rtsrc, 'RTSRC :3:hist:anal:mpti:mpt3', &
            cuparmm(ng)%rtsrc, imean)
    end if

    if (associated(cuparm(ng)%clsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%clsrc, 'CLSRC :3:hist:anal:mpti:mpt3', &
            cuparmm(ng)%clsrc, imean)
    end if

    if (associated(cuparm(ng)%thsrcp)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%thsrcp, 'THSRCP :3:mpti:', &
            cuparmm(ng)%thsrcp, imean)
    end if

    if (associated(cuparm(ng)%rtsrcp)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%rtsrcp, 'RTSRCP :3:mpti:', &
            cuparmm(ng)%rtsrcp, imean)
    end if

    if (associated(cuparm(ng)%thsrcf)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%thsrcf, 'THSRCF :3:mpti:', &
            cuparmm(ng)%thsrcf, imean)
    end if

    if (associated(cuparm(ng)%rtsrcf)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%rtsrcf, 'RTSRCF :3:mpti:', &
            cuparmm(ng)%rtsrcf, imean)
    end if

    if (associated(cuparm(ng)%aconpr)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%aconpr, 'ACONPR :2:hist:anal:mpti:mpt3', &
            cuparmm(ng)%aconpr, imean)
    end if

    if (associated(cuparm(ng)%conprr)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%conprr, 'CONPRR :2:hist:anal:mpt3', &
            cuparmm(ng)%conprr, imean)
    end if

    if (associated(cuparm(ng)%conprrp)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%conprrp, 'CONPRRP :2:mpti', &
            cuparmm(ng)%conprrp, imean)
    end if

    if (associated(cuparm(ng)%conprrf)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            cuparm(ng)%conprrf, 'CONPRRF :2:mpti', &
            cuparmm(ng)%conprrf, imean)
    end if

  end subroutine filltab_cuparm
end module mem_cuparm
