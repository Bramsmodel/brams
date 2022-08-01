!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_cuparm

  use ModNamelistFile, only: &
       namelistFile

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

!!$  type (cuparm_vars), allocatable :: cuparm_g(:)
!!$  type (cuparm_vars), allocatable :: cuparmm_g(:)
!!$  type (cuparm_vars), allocatable :: cuparm_g_sh(:)
!!$  type (cuparm_vars), allocatable :: cuparmm_g_sh(:)

  type (cuparm_vars), pointer :: cuparm_g(:) => null()
  type (cuparm_vars), pointer :: cuparmm_g(:) => null()
  type (cuparm_vars), pointer :: cuparm_g_sh(:) => null()
  type (cuparm_vars), pointer :: cuparmm_g_sh(:) => null()

  include "files.h"

  integer, parameter :: maxcufiles = maxfiles
  integer, parameter :: maxcugrids = 10

  integer :: if_cuinv             ! from RAMSIN
  real :: tcu_beg                 ! from RAMSIN
  real :: tcu_end                 ! from RAMSIN
  real :: cu_til                  ! from RAMSIN
  real :: cu_tel                  ! from RAMSIN
  real :: tnudcu                  ! from RAMSIN
  real :: wt_cu_grid(maxcugrids)  ! from RAMSIN
  character(len=128) :: cu_prefix ! from RAMSIN

  character(len=f_name_length) :: fnames_cu(maxcufiles)
  character(len=14)  :: itotdate_cu(maxcufiles)
  real :: cu_times(maxcufiles)

  integer :: ncufiles
  integer :: ncufl
  real :: cutime1
  real :: cutime2

  integer :: nnqparm(maxgrds) ! from RAMSIN
  real :: wcldbs              ! from RAMSIN
  real :: confrq              ! from RAMSIN

  public :: hasAconpr

contains

  logical function hasAconpr(ng)
    integer, intent(IN) :: ng

    hasAconpr = nnqparm(ng)>= 1 .or. if_cuinv == 1

  end function hasAconpr

  subroutine alloc_cuparm(cuparm, n1, n2, n3, ng)
    ! Arguments:
    type (cuparm_vars), intent(INOUT) :: cuparm
    integer, intent(IN)               :: n1, n2, n3, ng
    ! Local Variables:
    character(len=*), parameter :: h="**(alloc_cuparm)**"
    integer :: ierr

    ! Allocate arrays based on options (if necessary)

    if( hasAconpr(ng) )  then
       if(nnqparm(ng) < 3) then !srf: for conv 3 and 4, special arrays
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

       if (if_cuinv == 1) then
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
    type (cuparm_vars), intent(INOUT) :: cuparm

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
    type (cuparm_vars), intent(INOUT) :: cuparm

    if (associated(cuparm%thsrc))    deallocate (cuparm%thsrc)
    if (associated(cuparm%rtsrc))    deallocate (cuparm%rtsrc)
    if (associated(cuparm%clsrc))    deallocate (cuparm%clsrc)
    if (associated(cuparm%thsrcp))    deallocate (cuparm%thsrcp)
    if (associated(cuparm%rtsrcp))    deallocate (cuparm%rtsrcp)
    if (associated(cuparm%thsrcf))    deallocate (cuparm%thsrcf)
    if (associated(cuparm%rtsrcf))    deallocate (cuparm%rtsrcf)
    if (associated(cuparm%aconpr))   deallocate (cuparm%aconpr)
    if (associated(cuparm%conprr))   deallocate (cuparm%conprr)
    if (associated(cuparm%conprrp))   deallocate (cuparm%conprrp)
    if (associated(cuparm%conprrf))   deallocate (cuparm%conprrf)

  end subroutine dealloc_cuparm

  ! ----------------------------------------------------------------------

  subroutine filltab_cuparm_sh(oneVarTable, oneVarTableSize, cuparm, cuparmm)

    use iso_fortran_env, only: &
         int64
    
    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable

    ! Arguments:
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type (cuparm_vars), pointer, intent(in) :: cuparm, cuparmm


    ! Fill pointers to arrays into variable tables

    if (associated(cuparm%thsrc)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%thsrc, 'THSRC_SH :3:hist:anal:mpti:mpt3', &
               cuparmm%thsrc)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%thsrc, 'THSRC_SH :3:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(cuparm%rtsrc)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%rtsrc, 'RTSRC_SH :3:hist:anal:mpti:mpt3', &
               cuparmm%rtsrc)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%rtsrc, 'RTSRC_SH :3:hist:anal:mpti:mpt3')
       end if
    end if

  end subroutine filltab_cuparm_sh

  ! ----------------------------------------------------------------------

  subroutine filltab_cuparm(oneVarTable, oneVarTableSize, cuparm, cuparmm)

    use iso_fortran_env, only: &
         int64
    
    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable

    ! Arguments:
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type (cuparm_vars), pointer, intent(in) :: cuparm, cuparmm

    ! Fill pointers to arrays into variable tables

    if (associated(cuparm%thsrc)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%thsrc, 'THSRC :3:hist:anal:mpti:mpt3', &
               cuparmm%thsrc)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%thsrc, 'THSRC :3:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(cuparm%rtsrc)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%rtsrc, 'RTSRC :3:hist:anal:mpti:mpt3', &
               cuparmm%rtsrc)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%rtsrc, 'RTSRC :3:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(cuparm%clsrc)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%clsrc, 'CLSRC :3:hist:anal:mpti:mpt3', &
               cuparmm%clsrc)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%clsrc, 'CLSRC :3:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(cuparm%thsrcp)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%thsrcp, 'THSRCP :3:mpti:', &
               cuparmm%thsrcp)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%thsrcp, 'THSRCP :3:mpti:')
       end if
    end if

    if (associated(cuparm%rtsrcp)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%rtsrcp, 'RTSRCP :3:mpti:', &
               cuparmm%rtsrcp)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%rtsrcp, 'RTSRCP :3:mpti:')
       end if
    end if

    if (associated(cuparm%thsrcf)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%thsrcf, 'THSRCF :3:mpti:', &
               cuparmm%thsrcf)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%thsrcf, 'THSRCF :3:mpti:')
       end if
    end if

    if (associated(cuparm%rtsrcf)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%rtsrcf, 'RTSRCF :3:mpti:', &
               cuparmm%rtsrcf)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%rtsrcf, 'RTSRCF :3:mpti:')
       end if
    end if

    if (associated(cuparm%aconpr)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%aconpr, 'ACONPR :2:hist:anal:mpti:mpt3', &
               cuparmm%aconpr)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%aconpr, 'ACONPR :2:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(cuparm%conprr)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%conprr, 'CONPRR :2:hist:anal:mpt3', &
               cuparmm%conprr)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%conprr, 'CONPRR :2:hist:anal:mpt3')
       end if
    end if

    if (associated(cuparm%conprrp)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%conprrp, 'CONPRRP :2:mpti', &
               cuparmm%conprrp)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%conprrp, 'CONPRRP :2:mpti')
       end if
    end if

    if (associated(cuparm%conprrf)) then
       if (associated(cuparmm)) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%conprrf, 'CONPRRF :2:mpti', &
               cuparmm%conprrf)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuparm%conprrf, 'CONPRRF :2:mpti')
       end if
    end if

  end subroutine filltab_cuparm




  subroutine StoreNamelistFileAtMem_cuparm(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    confrq = oneNamelistFile%confrq
    cu_prefix = oneNamelistFile%cu_prefix
    cu_tel = oneNamelistFile%cu_tel
    cu_til = oneNamelistFile%cu_til
    if_cuinv = oneNamelistFile%if_cuinv
    nnqparm = oneNamelistFile%nnqparm
    tcu_beg = oneNamelistFile%tcu_beg
    tcu_end = oneNamelistFile%tcu_end
    tnudcu = oneNamelistFile%tnudcu
    wcldbs = oneNamelistFile%wcldbs
    wt_cu_grid = oneNamelistFile%wt_cu_grid
  end subroutine StoreNamelistFileAtMem_cuparm

end module mem_cuparm
