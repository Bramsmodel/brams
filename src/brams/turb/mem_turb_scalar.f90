!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

! For CATT

module mem_turb_scalar

  implicit none

  type turb_s_vars
     real, pointer, contiguous ::  hksc(:,:,:) => null()
  end type turb_s_vars

  type (turb_s_vars), allocatable, target :: turb_s(:)
  type (turb_s_vars), allocatable, target :: turbm_s(:)


contains

  subroutine alloc_turb_s(turb_s_local, n1, n2, n3)

    implicit none

    type (turb_s_vars)  :: turb_s_local
    integer, intent(in) :: n1,n2,n3

    allocate (turb_s_local%hksc(n1,n2,n3))
    turb_s_local%hksc=0.0

    return
  end subroutine alloc_turb_s

  !---------------------------------------------------------------

  subroutine nullify_turb_s(turb_s_local)

    implicit none

    type (turb_s_vars) :: turb_s_local

    integer :: nsc

    ! Deallocate all scratch arrays

    if (associated(turb_s_local%hksc ))  nullify (turb_s_local%hksc )

    return
  end subroutine nullify_turb_s
  !---------------------------------------------------------------

  subroutine dealloc_turb_s(turb_s_local)

    implicit none

    type (turb_s_vars) :: turb_s_local

    integer :: nsc

    ! Deallocate all scratch arrays

    if (associated(turb_s_local%hksc ))  deallocate (turb_s_local%hksc )

    return
  end subroutine dealloc_turb_s

  !---------------------------------------------------------------

  subroutine filltab_turb_s(oneVarTable, oneVarTableSize, &
       turb_s, turbm_s, imean)

    ! Build VarTable entry with turb_s_vars components

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable
    
    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(turb_s_vars), pointer, intent(in) :: turb_s
    type(turb_s_vars), pointer, intent(in) :: turbm_s
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(filltab_turb_s)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(turb_s)) then
       call fatal_error(h//" turb_s not associated")
    end if
    
    if (associated(turb_s%hksc)) then
       call InsertVarTable (OneVarTable, oneVarTableSize, &
            turb_s%hksc, &
            'HKSC :3:hist:anal:mpti:mpt3:mpt1', &
            turbm_s%hksc, imean)
    end if
  end subroutine filltab_turb_s

end module mem_turb_scalar
