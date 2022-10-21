! Module for urban canopy parameterization (land classes 19 and 21).

module mem_teb

  type teb_vars
     ! Variables to be dimensioned by (3,nxp,nyp)
     real, pointer, contiguous :: T_ROOF(:,:,:)
     real, pointer, contiguous :: T_ROAD(:,:,:)
     real, pointer, contiguous :: T_WALL(:,:,:)
     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, contiguous :: T_CANYON(:,:)
     real, pointer, contiguous :: R_CANYON(:,:)
     real, pointer, contiguous :: TS_ROOF(:,:)
     real, pointer, contiguous :: TS_ROAD(:,:)
     real, pointer, contiguous :: TS_WALL(:,:)
     real, pointer, contiguous :: TI_ROAD(:,:)
     real, pointer, contiguous :: WS_ROOF(:,:)
     real, pointer, contiguous :: WS_ROAD(:,:)
     real, pointer, contiguous :: TI_BLD(:,:)
     real, pointer, contiguous :: LE_TRAFFIC(:,:)
     real, pointer, contiguous :: H_TRAFFIC(:,:)
     real, pointer, contiguous :: LE_INDUSTRY(:,:)
     real, pointer, contiguous :: H_INDUSTRY(:,:)
     real, pointer, contiguous :: T2M_TOWN(:,:)
     real, pointer, contiguous :: R2M_TOWN(:,:)
     real, pointer, contiguous :: fuso(:,:)
  end type teb_vars

  type (teb_vars), allocatable, target :: teb_g(:), tebm_g(:)

contains

  subroutine alloc_teb(teb,n1,n2,n3,ng)

    implicit none
    type (teb_vars) :: teb
    integer, intent(in) :: n1,n2,n3,ng


    allocate (teb%T_ROOF(n1,n2,n3),teb%T_ROAD(n1,n2,n3), &
         teb%T_WALL(n1,n2,n3) )

    allocate (teb%T_CANYON(n2,n3),teb%R_CANYON(n2,n3), &
         teb%TS_ROOF(n2,n3),teb%TS_ROAD(n2,n3),teb%TS_WALL(n2,n3), &
         teb%TI_ROAD(n2,n3),teb%WS_ROOF(n2,n3),teb%WS_ROAD(n2,n3), &
         teb%TI_BLD(n2,n3),teb%LE_TRAFFIC(n2,n3),teb%H_TRAFFIC(n2,n3), &
         teb%LE_INDUSTRY(n2,n3),teb%H_INDUSTRY(n2,n3),              &
         teb%T2M_TOWN(n2,n3),teb%R2M_TOWN(n2,n3),teb%fuso(n2,n3))


    return
  end subroutine alloc_teb


  subroutine nullify_teb(teb)

    implicit none
    type (teb_vars) :: teb


    if (associated(teb%T_ROOF))  nullify (teb%T_ROOF)
    if (associated(teb%T_ROAD))  nullify (teb%T_ROAD)
    if (associated(teb%T_WALL))  nullify (teb%T_WALL)

    if (associated(teb%T_CANYON))  nullify  (teb%T_CANYON)
    if (associated(teb%R_CANYON))  nullify  (teb%R_CANYON)
    if (associated(teb%TS_ROOF))   nullify  (teb%TS_ROOF)
    if (associated(teb%TS_ROAD))   nullify  (teb%TS_ROAD)
    if (associated(teb%TS_WALL))   nullify  (teb%TS_WALL)
    if (associated(teb%TI_ROAD))   nullify  (teb%TI_ROAD)
    if (associated(teb%WS_ROOF))   nullify  (teb%WS_ROOF)
    if (associated(teb%WS_ROAD))   nullify  (teb%WS_ROAD)
    if (associated(teb%TI_BLD))    nullify  (teb%TI_BLD)
    if (associated(teb%LE_TRAFFIC)) nullify (teb%LE_TRAFFIC)
    if (associated(teb%H_TRAFFIC))  nullify (teb%H_TRAFFIC)
    if (associated(teb%LE_INDUSTRY)) nullify (teb%LE_INDUSTRY)
    if (associated(teb%H_INDUSTRY))  nullify (teb%H_INDUSTRY)
    if (associated(teb%T2M_TOWN))    nullify (teb%T2M_TOWN)
    if (associated(teb%R2M_TOWN))    nullify (teb%R2M_TOWN)
    if (associated(teb%fuso))    nullify (teb%fuso)

    return
  end subroutine nullify_teb

  subroutine dealloc_teb(teb)

    implicit none
    type (teb_vars) :: teb

    if (associated(teb%T_ROOF))  deallocate (teb%T_ROOF)
    if (associated(teb%T_ROAD))  deallocate (teb%T_ROAD)
    if (associated(teb%T_WALL))  deallocate (teb%T_WALL)

    if (associated(teb%T_CANYON))  deallocate  (teb%T_CANYON)
    if (associated(teb%R_CANYON))  deallocate  (teb%R_CANYON)
    if (associated(teb%TS_ROOF))   deallocate  (teb%TS_ROOF)
    if (associated(teb%TS_ROAD))   deallocate  (teb%TS_ROAD)
    if (associated(teb%TS_WALL))   deallocate  (teb%TS_WALL)
    if (associated(teb%TI_ROAD))   deallocate  (teb%TI_ROAD)
    if (associated(teb%WS_ROOF))   deallocate  (teb%WS_ROOF)
    if (associated(teb%WS_ROAD))   deallocate  (teb%WS_ROAD)
    if (associated(teb%TI_BLD))    deallocate  (teb%TI_BLD)
    if (associated(teb%LE_TRAFFIC)) deallocate (teb%LE_TRAFFIC)
    if (associated(teb%H_TRAFFIC))  deallocate (teb%H_TRAFFIC)
    if (associated(teb%LE_INDUSTRY)) deallocate (teb%LE_INDUSTRY)
    if (associated(teb%H_INDUSTRY))  deallocate (teb%H_INDUSTRY)

    if (associated(teb%T2M_TOWN))    deallocate (teb%T2M_TOWN)
    if (associated(teb%R2M_TOWN))    deallocate (teb%R2M_TOWN)
    if (associated(teb%fuso))    deallocate (teb%fuso)

    return
  end subroutine dealloc_teb


  subroutine filltab_teb(oneVarTable, oneVarTableSize, &
       teb, tebm, imean)

    ! Build VarTable entry with teb_vars components

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(teb_vars), pointer, intent(in) :: teb
    type(teb_vars), pointer, intent(in) :: tebm
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(filltab_ted)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    end if

    if (associated(teb%T_ROOF)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%T_ROOF, & 
            'T_ROOF :3:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%T_ROOF, imean)
    end if

    if (associated(teb%T_ROAD)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%T_ROAD, & 
            'T_ROAD :3:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%T_ROAD, imean)
    end if

    if (associated(teb%T_WALL)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%T_WALL, & 
            'T_WALL :3:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%T_WALL, imean)
    end if

    if (associated(teb%T_CANYON)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%T_CANYON, & 
            'T_CANYON :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%T_CANYON, imean)
    end if

    if (associated(teb%R_CANYON)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%R_CANYON, & 
            'R_CANYON :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%R_CANYON, imean)
    end if

    if (associated(teb%TS_ROOF)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%TS_ROOF, & 
            'TS_ROOF :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%TS_ROOF, imean)
    end if

    if (associated(teb%TS_ROAD)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%TS_ROAD, & 
            'TS_ROAD :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%TS_ROAD, imean)
    end if

    if (associated(teb%TS_WALL)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%TS_WALL, & 
            'TS_WALL :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%TS_WALL, imean)
    end if

    if (associated(teb%TI_ROAD)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%TI_ROAD, & 
            'TI_ROAD :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%TI_ROAD, imean)
    end if

    if (associated(teb%WS_ROOF)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%WS_ROOF, & 
            'WS_ROOF :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%WS_ROOF, imean)
    end if

    if (associated(teb%WS_ROAD)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%WS_ROAD, & 
            'WS_ROAD :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%WS_ROAD, imean)
    end if

    if (associated(teb%TI_BLD)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%TI_BLD, & 
            'TI_BLD :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%TI_BLD, imean)
    end if

    if (associated(teb%LE_TRAFFIC)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%LE_TRAFFIC, & 
            'LE_TRAFFIC :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%LE_TRAFFIC, imean)
    end if

    if (associated(teb%H_TRAFFIC)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%H_TRAFFIC, & 
            'H_TRAFFIC :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%H_TRAFFIC, imean)
    end if

    if (associated(teb%LE_INDUSTRY)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%LE_INDUSTRY, & 
            'LE_INDUSTRY :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%LE_INDUSTRY, imean)
    end if

    if (associated(teb%H_INDUSTRY)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%H_INDUSTRY, & 
            'H_INDUSTRY :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%H_INDUSTRY, imean)
    end if

    if (associated(teb%T2M_TOWN)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%T2M_TOWN, & 
            'T2M_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%T2M_TOWN, imean)
    end if

    if (associated(teb%R2M_TOWN)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%R2M_TOWN, & 
            'R2M_TOWN :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%R2M_TOWN, imean)
    end if

    if (associated(teb%fuso)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            teb%fuso, & 
            'FUSO :2:hist:anal:lite:mpti:mpt3:mpt1', &
            tebm%fuso, imean)
    end if
  end subroutine filltab_teb

end module mem_teb
