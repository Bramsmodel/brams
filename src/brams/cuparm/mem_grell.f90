! Module necessary to GRELL param.

module mem_grell

  use ModNamelistFile, only: &
       NamelistFile
  
  type cuforc_vars
     real, pointer, contiguous :: lsfth(:,:,:) => null()
     real, pointer, contiguous :: lsfrt(:,:,:) => null()
  end type cuforc_vars

  type (cuforc_vars), pointer, contiguous :: cuforc_g(:) => null()
  type (cuforc_vars), pointer, contiguous :: cuforcm_g(:) => null()
  type (cuforc_vars), pointer, contiguous :: cuforc_sh_g(:) => null()
  type (cuforc_vars), pointer, contiguous :: cuforcm_sh_g(:) => null()



contains

  subroutine alloc_cu_forcings(cuforc, m1, m2, m3, ng)

    !    use shcu_vars_const, only : nnshcu   

    implicit none
    type (cuforc_vars) :: cuforc
    integer, intent(in) :: m1, m2, m3, ng

    allocate (cuforc%lsfth(m1, m2, m3));cuforc%lsfth=0.0
    allocate (cuforc%lsfrt(m1, m2, m3));cuforc%lsfrt=0.0

  end subroutine alloc_cu_forcings
  !---------------------------------------------------------------
  subroutine nullify_cuforc(cuforc)

    implicit none
    type (cuforc_vars) :: cuforc

    if (associated(cuforc%lsfth))   nullify (cuforc%lsfth)
    if (associated(cuforc%lsfrt))   nullify (cuforc%lsfrt)

  end subroutine nullify_cuforc
  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine filltab_cuforc_sh(oneVarTable, oneVarTableSize, &
       cuforc, cuforcm, imean)

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(cuforc_vars), pointer, intent(in) :: cuforc
    type(cuforc_vars), pointer, intent(in) :: cuforcm
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(filltab_cuforc_sh)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(cuforc)) then
       call fatal_error(h//" cuforc not associated")
    end if

!TO    if (associated(cuforc%lsfth)) then
!TO       call InsertVarTable(oneVarTable, oneVarTableSize, &
!TO            cuforc%lsfth, &
!TO            'LSFTH_SH :3:hist:anal:mpti:mpt3', &
!TO            cuforcm%lsfth, imean)
!TO    end if

!TO    if (associated(cuforc%lsfrt)) then
!TO       call InsertVarTable(oneVarTable, oneVarTableSize, &
!TO            cuforc%lsfrt, &
!TO            'LSFRT_SH :3:hist:anal:mpti:mpt3', &
!TO            cuforcm%lsfrt, imean)
!TO    end if
  end subroutine filltab_cuforc_sh

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine filltab_cuforc(oneVarTable, oneVarTableSize, &
       cuforc, cuforcm, imean)

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(cuforc_vars), pointer, intent(in) :: cuforc
    type(cuforc_vars), pointer, intent(in) :: cuforcm
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(filltab_cuforc)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(cuforc)) then
       call fatal_error(h//" cuforc not associated")
    end if

!TO    if (associated(cuforc%lsfth)) then
!TO       call InsertVarTable (oneVarTable, oneVarTableSize, &
!TO            cuforc%lsfth, &
!TO            'LSFTH :3:hist:anal:mpti:mpt3', &
!TO            cuforcm%lsfth, imean)
!TO    end if

!TO    if (associated(cuforc%lsfrt)) then
!TO       call InsertVarTable (oneVarTable, oneVarTableSize, &
!TO            cuforc%lsfrt, &
!TO            'LSFRT :3:hist:anal:mpti:mpt3', &
!TO            cuforcm%lsfrt, imean)
!TO    end if
  end subroutine filltab_cuforc

end module mem_grell
