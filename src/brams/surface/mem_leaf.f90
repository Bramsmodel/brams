!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module mem_leaf

  use ModNamelistFile, only: &
       namelistFile

  use grid_dims, only: &
       nzgmax

  Type leaf_vars

     ! Variables to be dimensioned by (nxp,nyp,nzg,npatch)

     real, pointer, contiguous :: soil_water(:,:,:,:)
     real, pointer, contiguous :: soil_energy(:,:,:,:)
     real, pointer, contiguous :: soil_text(:,:,:,:)

     ! Variables to be dimensioned by (nxp,nyp,nzs,npatch)

     real, pointer, contiguous :: sfcwater_mass(:,:,:,:)
     real, pointer, contiguous :: sfcwater_energy(:,:,:,:)
     real, pointer, contiguous :: sfcwater_depth(:,:,:,:)

     ! Variables to be dimensioned by (nxp,nyp,npatch)

     real, pointer, contiguous :: ustar(:,:,:)
     real, pointer, contiguous :: tstar(:,:,:)
     real, pointer, contiguous :: rstar(:,:,:)
     real, pointer, contiguous :: veg_fracarea(:,:,:)
     real, pointer, contiguous :: veg_lai(:,:,:)
     real, pointer, contiguous :: veg_rough(:,:,:)
     real, pointer, contiguous :: veg_height(:,:,:)
     real, pointer, contiguous :: veg_albedo(:,:,:)
     real, pointer, contiguous :: veg_tai(:,:,:)
     real, pointer, contiguous :: patch_area(:,:,:)
     real, pointer, contiguous :: patch_rough(:,:,:)
     real, pointer, contiguous :: patch_wetind(:,:,:)
     real, pointer, contiguous :: leaf_class(:,:,:)
     real, pointer, contiguous :: soil_rough(:,:,:)
     real, pointer, contiguous :: sfcwater_nlev(:,:,:)
     real, pointer, contiguous :: stom_resist(:,:,:)
     real, pointer, contiguous :: ground_rsat(:,:,:)
     real, pointer, contiguous :: ground_rvap(:,:,:)
     real, pointer, contiguous :: veg_water(:,:,:)
     real, pointer, contiguous :: veg_temp(:,:,:)
     real, pointer, contiguous :: can_rvap(:,:,:)
     real, pointer, contiguous :: can_temp(:,:,:)
     real, pointer, contiguous :: veg_ndvip(:,:,:)
     real, pointer, contiguous :: veg_ndvic(:,:,:)
     real, pointer, contiguous :: veg_ndvif(:,:,:)

     real, pointer, contiguous :: R_aer(:,:,:)   !kml drydep

     ! Variables to be dimensioned by (nxp,nyp)

     real, pointer, contiguous :: snow_mass(:,:)
     real, pointer, contiguous :: snow_depth(:,:)
     real, pointer, contiguous :: seatp(:,:)
     real, pointer, contiguous :: seatf(:,:)
  End Type leaf_vars

  type (leaf_vars), pointer :: leaf_g(:) => null()
  type (leaf_vars), pointer :: leafm_g(:) => null()

  !----------------------------------------------------------------------------
  integer :: nslcon ! from RAMSIN
  integer :: nvgcon ! from RAMSIN
  integer :: nvegpat ! from RAMSIN
  integer :: isfcl ! from RAMSIN
  integer :: isfcl_ocean ! from RAMSIN
  real    :: zrough ! from RAMSIN
  real    :: pctlcon ! from RAMSIN
  real    :: ubmin
  real    :: albedo ! from RAMSIN
  real    :: drtcon ! from RAMSIN
  real    :: dthcon ! from RAMSIN
  real    :: seatmp ! from RAMSIN
  real    :: stgoff(nzgmax) ! from RAMSIN
  real    :: slmstr(nzgmax) ! from RAMSIN
  real    :: slz(nzgmax) ! from RAMSIN

Contains

  subroutine alloc_leaf(leaf,nz,nx,ny,nzg,nzs,np,ng)

    implicit none
    type (leaf_vars) :: leaf
    integer, intent(in) :: nz,nx,ny,nzg,nzs,np,ng

    ! Allocate arrays based on options (if necessary)

    allocate (leaf%soil_water     (nzg,nx,ny,np));leaf%soil_water =0.0
    allocate (leaf%soil_energy    (nzg,nx,ny,np));leaf%soil_energy=0.0
    allocate (leaf%soil_text      (nzg,nx,ny,np));leaf%soil_text  =0.0

    allocate (leaf%sfcwater_mass  (nzs,nx,ny,np));leaf%sfcwater_mass  =0.0
    allocate (leaf%sfcwater_energy(nzs,nx,ny,np));leaf%sfcwater_energy=0.0
    allocate (leaf%sfcwater_depth (nzs,nx,ny,np));leaf%sfcwater_depth =0.0

    allocate (leaf%ustar        (nx,ny,np));leaf%ustar=0.0
    allocate (leaf%tstar        (nx,ny,np));leaf%tstar=0.0
    allocate (leaf%rstar        (nx,ny,np));leaf%rstar=0.0

    allocate (leaf%veg_fracarea (nx,ny,np));leaf%veg_fracarea=0.0
    allocate (leaf%veg_lai      (nx,ny,np));leaf%veg_lai     =0.0
    allocate (leaf%veg_rough    (nx,ny,np));leaf%veg_rough   =0.0
    allocate (leaf%veg_height   (nx,ny,np));leaf%veg_height  =0.0
    allocate (leaf%veg_albedo   (nx,ny,np));leaf%veg_albedo  =0.0
    allocate (leaf%veg_tai      (nx,ny,np));leaf%veg_tai     =0.0

    allocate (leaf%patch_area   (nx,ny,np));leaf%patch_area  =0.0
    allocate (leaf%patch_rough  (nx,ny,np));leaf%patch_rough =0.0
    allocate (leaf%patch_wetind (nx,ny,np));leaf%patch_wetind=0.0
    allocate (leaf%leaf_class   (nx,ny,np));leaf%leaf_class  =0.0

    allocate (leaf%soil_rough   (nx,ny,np));leaf%soil_rough   =0.0
    allocate (leaf%sfcwater_nlev(nx,ny,np));leaf%sfcwater_nlev=0.0
    allocate (leaf%stom_resist  (nx,ny,np));leaf%stom_resist  =0.0

    allocate (leaf%ground_rsat  (nx,ny,np));leaf%ground_rsat=0.0
    allocate (leaf%ground_rvap  (nx,ny,np));leaf%ground_rvap=0.0
    ;
    allocate (leaf%veg_water    (nx,ny,np));leaf%veg_water  =0.0
    allocate (leaf%veg_temp     (nx,ny,np));leaf%veg_temp   =0.0
    ;                                          
    allocate (leaf%can_rvap     (nx,ny,np));leaf%can_rvap   =0.0
    allocate (leaf%can_temp     (nx,ny,np));leaf%can_temp   =0.0
    ;                                          
    allocate (leaf%veg_ndvip    (nx,ny,np));leaf%veg_ndvip  =0.0
    allocate (leaf%veg_ndvic    (nx,ny,np));leaf%veg_ndvic  =0.0
    allocate (leaf%veg_ndvif    (nx,ny,np));leaf%veg_ndvif  =0.0

    allocate (leaf%R_aer        (nx,ny,np));leaf%R_aer =0.0  !kml drydep

    allocate (leaf%snow_mass    (nx,ny));leaf%snow_mass =0.0
    allocate (leaf%snow_depth   (nx,ny));leaf%snow_depth=0.0
    allocate (leaf%seatp        (nx,ny));leaf%seatp     =0.0
    allocate (leaf%seatf        (nx,ny));leaf%seatf     =0.0

  end subroutine alloc_leaf

  !************************************************************************

  subroutine nullify_leaf(leaf)

    implicit none
    type (leaf_vars) :: leaf

    if(associated(leaf%soil_water))      nullify (leaf%soil_water)
    if(associated(leaf%soil_energy))     nullify (leaf%soil_energy)
    if(associated(leaf%soil_text))       nullify (leaf%soil_text)

    if(associated(leaf%sfcwater_mass))   nullify (leaf%sfcwater_mass)
    if(associated(leaf%sfcwater_energy)) nullify (leaf%sfcwater_energy)
    if(associated(leaf%sfcwater_depth))  nullify (leaf%sfcwater_depth)

    if(associated(leaf%ustar))           nullify (leaf%ustar)
    if(associated(leaf%tstar))           nullify (leaf%tstar)
    if(associated(leaf%rstar))           nullify (leaf%rstar)

    if(associated(leaf%veg_fracarea))    nullify (leaf%veg_fracarea)
    if(associated(leaf%veg_lai))         nullify (leaf%veg_lai)
    if(associated(leaf%veg_rough))       nullify (leaf%veg_rough)
    if(associated(leaf%veg_height))      nullify (leaf%veg_height)
    if(associated(leaf%veg_albedo))      nullify (leaf%veg_albedo)
    if(associated(leaf%veg_tai))         nullify (leaf%veg_tai)

    if(associated(leaf%patch_area))      nullify (leaf%patch_area)
    if(associated(leaf%patch_rough))     nullify (leaf%patch_rough)
    if(associated(leaf%patch_wetind))    nullify (leaf%patch_wetind)
    if(associated(leaf%leaf_class))      nullify (leaf%leaf_class)

    if(associated(leaf%soil_rough))      nullify (leaf%soil_rough)
    if(associated(leaf%sfcwater_nlev))   nullify (leaf%sfcwater_nlev)
    if(associated(leaf%stom_resist))     nullify (leaf%stom_resist)

    if(associated(leaf%ground_rsat))     nullify (leaf%ground_rsat)
    if(associated(leaf%ground_rvap))     nullify (leaf%ground_rvap)

    if(associated(leaf%veg_water))       nullify (leaf%veg_water)
    if(associated(leaf%veg_temp))        nullify (leaf%veg_temp)

    if(associated(leaf%can_rvap))        nullify (leaf%can_rvap)
    if(associated(leaf%can_temp))        nullify (leaf%can_temp)

    if(associated(leaf%veg_ndvip))       nullify (leaf%veg_ndvip)
    if(associated(leaf%veg_ndvic))       nullify (leaf%veg_ndvic)
    if(associated(leaf%veg_ndvif))       nullify (leaf%veg_ndvif)

    if(associated(leaf%R_aer))           nullify (leaf%R_aer)   !kml drydep

    if(associated(leaf%snow_mass))       nullify (leaf%snow_mass)
    if(associated(leaf%snow_depth))      nullify (leaf%snow_depth)
    if(associated(leaf%seatp))           nullify (leaf%seatp)
    if(associated(leaf%seatf))           nullify (leaf%seatf)

  end subroutine nullify_leaf

  ! ********************************************************************

  subroutine dealloc_leaf(leaf)

    implicit none
    type (leaf_vars) :: leaf

    if(associated(leaf%soil_water))      deallocate (leaf%soil_water)
    if(associated(leaf%soil_energy))     deallocate (leaf%soil_energy)
    if(associated(leaf%soil_text))       deallocate (leaf%soil_text)

    if(associated(leaf%sfcwater_mass))   deallocate (leaf%sfcwater_mass)
    if(associated(leaf%sfcwater_energy)) deallocate (leaf%sfcwater_energy)
    if(associated(leaf%sfcwater_depth))  deallocate (leaf%sfcwater_depth)

    if(associated(leaf%ustar))           deallocate (leaf%ustar)
    if(associated(leaf%tstar))           deallocate (leaf%tstar)
    if(associated(leaf%rstar))           deallocate (leaf%rstar)

    if(associated(leaf%veg_fracarea))    deallocate (leaf%veg_fracarea)
    if(associated(leaf%veg_lai))         deallocate (leaf%veg_lai)
    if(associated(leaf%veg_rough))       deallocate (leaf%veg_rough)
    if(associated(leaf%veg_height))      deallocate (leaf%veg_height)
    if(associated(leaf%veg_albedo))      deallocate (leaf%veg_albedo)
    if(associated(leaf%veg_tai))         deallocate (leaf%veg_tai)

    if(associated(leaf%patch_area))      deallocate (leaf%patch_area)
    if(associated(leaf%patch_rough))     deallocate (leaf%patch_rough)
    if(associated(leaf%patch_wetind))    deallocate (leaf%patch_wetind)
    if(associated(leaf%leaf_class))      deallocate (leaf%leaf_class)

    if(associated(leaf%soil_rough))      deallocate (leaf%soil_rough)
    if(associated(leaf%sfcwater_nlev))   deallocate (leaf%sfcwater_nlev)
    if(associated(leaf%stom_resist))     deallocate (leaf%stom_resist)

    if(associated(leaf%ground_rsat))     deallocate (leaf%ground_rsat)
    if(associated(leaf%ground_rvap))     deallocate (leaf%ground_rvap)

    if(associated(leaf%veg_water))       deallocate (leaf%veg_water)
    if(associated(leaf%veg_temp))        deallocate (leaf%veg_temp)

    if(associated(leaf%can_rvap))        deallocate (leaf%can_rvap)
    if(associated(leaf%can_temp))        deallocate (leaf%can_temp)

    if(associated(leaf%veg_ndvip))       deallocate (leaf%veg_ndvip)
    if(associated(leaf%veg_ndvic))       deallocate (leaf%veg_ndvic)
    if(associated(leaf%veg_ndvif))       deallocate (leaf%veg_ndvif)

    if(associated(leaf%R_aer))           deallocate (leaf%R_aer)   !kml drydep

    if(associated(leaf%snow_mass))       deallocate (leaf%snow_mass)
    if(associated(leaf%snow_depth))      deallocate (leaf%snow_depth)
    if(associated(leaf%seatp))           deallocate (leaf%seatp)
    if(associated(leaf%seatf))           deallocate (leaf%seatf)

  end subroutine dealloc_leaf

  ! ********************************************************************

  subroutine filltab_leaf(oneVarTable, oneVarTableSize, &
       leaf, leafm, imean)

    use io_params, only: &
         ipastin

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:) 
    integer, intent(inout) :: oneVarTableSize
    type (leaf_vars), pointer, intent(in) :: leaf
    type (leaf_vars), pointer, intent(in) :: leafm
    integer, intent(in) :: imean

    character(len=8) :: str_recycle

    str_recycle = ''
    if (ipastin == 1) then
       str_recycle = ':recycle'
    endif

    ! Fill pointers to arrays into variable tables

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%soil_water, &
         'SOIL_WATER :4:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%soil_water, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%soil_energy, &
         'SOIL_ENERGY :4:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%soil_energy, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%soil_text, &
         'SOIL_TEXT :4:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%soil_text, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%sfcwater_mass, &
         'SFCWATER_MASS :5:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%sfcwater_mass, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%sfcwater_energy,  &
         'SFCWATER_ENERGY :5:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%sfcwater_energy, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%sfcwater_depth, &
         'SFCWATER_DEPTH :5:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%sfcwater_depth, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%ustar, &
         'USTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%ustar, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%tstar, &
         'TSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%tstar, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%rstar, &
         'RSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%rstar, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_fracarea, &
         'VEG_FRACAREA :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%veg_fracarea, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_lai, &
         'VEG_LAI :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%veg_lai, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_rough, &
         'VEG_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%veg_rough, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_height, &
         'VEG_HEIGHT :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%veg_height, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_albedo, &
         'VEG_ALBEDO :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%veg_albedo, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_tai, &
         'VEG_TAI :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%veg_tai, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%patch_area, &
         'PATCH_AREA :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%patch_area, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%patch_rough, &
         'PATCH_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%patch_rough, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%patch_wetind, &
         'PATCH_WETIND :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%patch_wetind, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%leaf_class, &
         'LEAF_CLASS :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%leaf_class, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%soil_rough, &
         'SOIL_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%soil_rough, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%sfcwater_nlev, &
         'SFCWATER_NLEV :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%sfcwater_nlev, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%stom_resist, &
         'STOM_RESIST :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%stom_resist, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%ground_rsat, &
         'GROUND_RSAT :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%ground_rsat, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%ground_rvap, &
         'GROUND_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%ground_rvap, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_water, &
         'VEG_WATER :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%veg_water, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_temp, &
         'VEG_TEMP :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%veg_temp, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%can_rvap, &
         'CAN_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%can_rvap, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%can_temp, &
         'CAN_TEMP :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%can_temp, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_ndvip, &
         'VEG_NDVIP :6:hist:mpti', &
         leafm%veg_ndvip, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_ndvic, &
         'VEG_NDVIC :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         leafm%veg_ndvic, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%veg_ndvif,  &
         'VEG_NDVIF :6:hist:mpti', &
         leafm%veg_ndvif, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%R_aer, &
         'R_aer :6:hist:mpti', &
         leafm%R_aer, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%snow_mass, &
         'SNOW_MASS :2:mpti', &
         leafm%snow_mass, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%snow_depth, &
         'SNOW_DEPTH :2:mpti', &
         leafm%snow_depth, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%seatp, &
         'SEATP :2:mpti', &
         leafm%seatp, imean)

    call InsertVarTable (oneVarTable, oneVarTableSize, &
         leaf%seatf, &
         'SEATF :2:mpti', &
         leafm%seatf, imean)

  end subroutine filltab_leaf

  subroutine StoreNamelistFileAtMem_leaf(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    albedo = oneNamelistFile%albedo
    drtcon = oneNamelistFile%drtcon
    dthcon = oneNamelistFile%dthcon
    isfcl = oneNamelistFile%isfcl
    isfcl_ocean = oneNamelistFile%isfcl_ocean
    nslcon = oneNamelistFile%nslcon
    nvegpat = oneNamelistFile%nvegpat
    nvgcon = oneNamelistFile%nvgcon
    pctlcon = oneNamelistFile%pctlcon
    seatmp = oneNamelistFile%seatmp
    slmstr = oneNamelistFile%slmstr
    slz = oneNamelistFile%slz
    stgoff = oneNamelistFile%stgoff
    zrough = oneNamelistFile%zrough
    !-- if the ocean model is undefined, use the one 
    !-- in JULES or LEAF.
    if(isfcl_ocean == -999 ) then
       isfcl_ocean = isfcl
    endif
  end subroutine StoreNamelistFileAtMem_leaf
End Module mem_leaf
