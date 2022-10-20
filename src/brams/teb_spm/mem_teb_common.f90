module mem_teb_common

  type teb_common
     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, contiguous :: EMIS_TOWN(:,:) => null()
     real, pointer, contiguous :: ALB_TOWN(:,:) => null()
     real, pointer, contiguous :: TS_TOWN(:,:) => null()
  end type teb_common

  type(teb_common), allocatable, target :: tebc_g(:)
  type(teb_common), allocatable, target :: tebcm_g(:)

contains


  
  subroutine alloc_tebc(tebc,n1,n2,n3,ng)

    implicit none
    type(teb_common), intent(inout) :: tebc
    integer, intent(in) :: n1,n2,n3,ng

    allocate (tebc%EMIS_TOWN(n2,n3),tebc%ALB_TOWN(n2,n3), &
         tebc%TS_TOWN(n2,n3))

    return
  end subroutine alloc_tebc


  subroutine nullify_tebc(tebc)

    implicit none
    type(teb_common), intent(inout) :: tebc

    if (associated(tebc%EMIS_TOWN))  nullify (tebc%EMIS_TOWN)
    if (associated(tebc%ALB_TOWN))   nullify (tebc%ALB_TOWN)
    if (associated(tebc%TS_TOWN))    nullify (tebc%TS_TOWN)

    return
  end subroutine nullify_tebc

  subroutine dealloc_tebc(tebc)

    implicit none

    type(teb_common), intent(inout) :: tebc

    if (associated(tebc%EMIS_TOWN))  deallocate (tebc%EMIS_TOWN)
    if (associated(tebc%ALB_TOWN))   deallocate (tebc%ALB_TOWN)
    if (associated(tebc%TS_TOWN))    deallocate (tebc%TS_TOWN)

    return
  end subroutine dealloc_tebc


  subroutine filltab_tebc (oneVarTable, oneVarTableSize, &
       tebc, tebcm, imean)

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    implicit none

    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(teb_common), intent(in) :: tebc
    type(teb_common), intent(in) :: tebcm
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(filltab_tebc)**"


    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    end if
       
    ! Fill pointers to arrays into variable tables

    if (associated(tebc%EMIS_TOWN)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            tebc%EMIS_TOWN, &
            'EMIS_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebcm%EMIS_TOWN, &
            imean)
    end if
    
    if (associated(tebc%ALB_TOWN)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            tebc%ALB_TOWN, &
            'ALB_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebcm%ALB_TOWN, &
            imean)
    end if
    
    if (associated(tebc%TS_TOWN)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            tebc%TS_TOWN, &
            'TS_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebcm%TS_TOWN, &
            imean)
    end if
  end subroutine filltab_tebc

end module mem_teb_common
