!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module mem_leaf

  use ReadBcst, only: &
       ReadStoreOwnChunk

  use ModControlVars, only: &
       ControlVars

  use ModNamelistFile, only: &
       NamelistFile

  use ModParallelEnvironment, only: &
       MsgDump

  use grid_dims, only: &
       nzgmax

  use node_mod, only: &
       mynum

  Type leaf_vars

     ! Variables to be dimensioned by (nxp,nyp,nzg,npatch)

     real, pointer, contiguous :: soil_water(:,:,:,:) => null()
     real, pointer, contiguous :: soil_energy(:,:,:,:) => null()
     real, pointer, contiguous :: soil_text(:,:,:,:) => null()

     ! Variables to be dimensioned by (nxp,nyp,nzs,npatch)

     real, pointer, contiguous :: sfcwater_mass(:,:,:,:) => null()
     real, pointer, contiguous :: sfcwater_energy(:,:,:,:) => null()
     real, pointer, contiguous :: sfcwater_depth(:,:,:,:) => null()

     ! Variables to be dimensioned by (nxp,nyp,npatch)

     real, pointer, contiguous :: ustar(:,:,:) => null()
     real, pointer, contiguous :: tstar(:,:,:) => null()
     real, pointer, contiguous :: rstar(:,:,:) => null()
     real, pointer, contiguous :: veg_fracarea(:,:,:) => null()
     real, pointer, contiguous :: veg_lai(:,:,:) => null()
     real, pointer, contiguous :: veg_rough(:,:,:)  => null()
     real, pointer, contiguous :: veg_height(:,:,:)  => null()
     real, pointer, contiguous :: veg_albedo(:,:,:)  => null()
     real, pointer, contiguous :: veg_tai(:,:,:)  => null()
     real, pointer, contiguous :: patch_area(:,:,:)  => null()
     real, pointer, contiguous :: patch_rough(:,:,:)  => null()
     real, pointer, contiguous :: patch_wetind(:,:,:)  => null()
     real, pointer, contiguous :: leaf_class(:,:,:)  => null()
     real, pointer, contiguous :: soil_rough(:,:,:)  => null()
     real, pointer, contiguous :: sfcwater_nlev(:,:,:)  => null()
     real, pointer, contiguous :: stom_resist(:,:,:)  => null()
     real, pointer, contiguous :: ground_rsat(:,:,:)  => null()
     real, pointer, contiguous :: ground_rvap(:,:,:)  => null()
     real, pointer, contiguous :: veg_water(:,:,:)  => null()
     real, pointer, contiguous :: veg_temp(:,:,:)  => null()
     real, pointer, contiguous :: can_rvap(:,:,:)  => null()
     real, pointer, contiguous :: can_temp(:,:,:)  => null()
     real, pointer, contiguous :: veg_ndvip(:,:,:)  => null()
     real, pointer, contiguous :: veg_ndvic(:,:,:)  => null()
     real, pointer, contiguous :: veg_ndvif(:,:,:)  => null()

     real, pointer, contiguous :: R_aer(:,:,:) => null()   !kml drydep

     ! Variables to be dimensioned by (nxp,nyp)

     real, pointer, contiguous :: snow_mass(:,:)  => null()
     real, pointer, contiguous :: snow_depth(:,:)  => null()
     real, pointer, contiguous :: seatp(:,:)  => null()
     real, pointer, contiguous :: seatf(:,:)  => null()
  End Type leaf_vars

  integer, parameter :: UNDEF_INT = huge(1)
  real, parameter :: UNDEF_REAL = huge(1.0)


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

  include "UseVfm.h"

Contains

  subroutine alloc_leaf(leaf,nz,nx,ny,nzg,nzs,np,ng)

    implicit none
    type (leaf_vars) :: leaf
    integer, intent(in) :: nz,nx,ny,nzg,nzs,np,ng

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(alloc_leaf)**"
    logical, parameter :: dumpLocal=.false.

    ! Allocate arrays based on options (if necessary)

    if (dumpLocal) then
       write(str(1),"(i8)") nzg
       write(str(2),"(i8)") nx
       write(str(3),"(i8)") ny
       write(str(4),"(i8)") np
       call MsgDump(h//" aloca componentes 4D de leaf dimensionados ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//")")
    end if

    allocate (leaf%soil_water     (nzg,nx,ny,np));leaf%soil_water =0.0
    allocate (leaf%soil_energy    (nzg,nx,ny,np));leaf%soil_energy=0.0
    allocate (leaf%soil_text      (nzg,nx,ny,np));leaf%soil_text  =0.0


    if (dumpLocal) then
       write(str(1),"(i8)") nzs
       write(str(2),"(i8)") nx
       write(str(3),"(i8)") ny
       write(str(4),"(i8)") np
       call MsgDump(h//" aloca componentes 4D de leaf dimensionados ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//")")
    end if

    allocate (leaf%sfcwater_mass  (nzs,nx,ny,np));leaf%sfcwater_mass  =0.0
    allocate (leaf%sfcwater_energy(nzs,nx,ny,np));leaf%sfcwater_energy=0.0
    allocate (leaf%sfcwater_depth (nzs,nx,ny,np));leaf%sfcwater_depth =0.0


    if (dumpLocal) then
       write(str(1),"(i8)") nx
       write(str(2),"(i8)") ny
       write(str(3),"(i8)") np
       call MsgDump(h//" aloca componentes 3D de leaf dimensionados ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

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

    if (dumpLocal) then
       write(str(1),"(i8)") nx
       write(str(2),"(i8)") ny
       call MsgDump(h//" aloca componentes 2D de leaf dimensionados ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//")")
    end if

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
    type (leaf_vars), target, intent(in) :: leaf
    type (leaf_vars), target, intent(in) :: leafm
    integer, intent(in) :: imean

    real, pointer, contiguous :: pLeaf4D(:,:,:,:) => null()
    real, pointer, contiguous :: pLeafM4D(:,:,:,:) => null()
    real, pointer, contiguous :: pLeaf3D(:,:,:) => null()
    real, pointer, contiguous :: pLeafM3D(:,:,:) => null()
    real, pointer, contiguous :: pLeaf2D(:,:) => null()
    real, pointer, contiguous :: pLeafM2D(:,:) => null()
    character(len=8) :: str_recycle
    character(len=*), parameter :: h="**(filltab_leaf)**"
    logical, parameter :: dumpLocal=.false.

    str_recycle = ''
    if (ipastin == 1) then
       str_recycle = ':recycle'
    endif


    ! Fill pointers to arrays into variable tables

    pLeaf4D => leaf%soil_water
    if (imean == 1) then
       pLeafM4D => leafm%soil_water
    else
       pLeafM4D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf4D, &
         'SOIL_WATER :4:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM4D, imean)

    if (dumpLocal) then
       call MsgDump(h//" insert soil_energy at var tables")
    end if
    
    pLeaf4D => leaf%soil_energy
    if (imean == 1) then
       pLeafM4D => leafm%soil_energy
    else
       pLeafM4D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf4D, &
         'SOIL_ENERGY :4:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM4D, imean)

    pLeaf4D => leaf%soil_text
    if (imean == 1) then
       pLeafM4D => leafm%soil_text
    else
       pLeafM4D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf4D, &
         'SOIL_TEXT :4:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM4D, imean)

    pLeaf4D => leaf%sfcwater_mass
    if (imean == 1) then
       pLeafM4D => leafm%sfcwater_mass
    else
       pLeafM4D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf4D, &
         'SFCWATER_MASS :5:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM4D, imean)

    pLeaf4D => leaf%sfcwater_energy
    if (imean == 1) then
       pLeafM4D => leafm%sfcwater_energy
    else
       pLeafM4D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf4D, &
         'SFCWATER_ENERGY :5:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM4D, imean)

    pLeaf4D => leaf%sfcwater_depth
    if (imean == 1) then
       pLeafM4D => leafm%sfcwater_depth
    else
       pLeafM4D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf4D, &
         'SFCWATER_DEPTH :5:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM4D, imean)

    pLeaf3D => leaf%ustar
    if (imean == 1) then
       pLeafM3D => leafm%ustar
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'USTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%tstar
    if (imean == 1) then
       pLeafM3D => leafm%tstar
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'TSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%rstar
    if (imean == 1) then
       pLeafM3D => leafm%rstar
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'RSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_fracarea
    if (imean == 1) then
       pLeafM3D => leafm%veg_fracarea
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_FRACAREA :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_lai
    if (imean == 1) then
       pLeafM3D => leafm%veg_lai
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_LAI :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_rough
    if (imean == 1) then
       pLeafM3D => leafm%veg_rough
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_height
    if (imean == 1) then
       pLeafM3D => leafm%veg_height
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_HEIGHT :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_albedo
    if (imean == 1) then
       pLeafM3D => leafm%veg_albedo
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_ALBEDO :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_tai
    if (imean == 1) then
       pLeafM3D => leafm%veg_tai
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_TAI :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%patch_area
    if (imean == 1) then
       pLeafM3D => leafm%patch_area
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'PATCH_AREA :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%patch_rough
    if (imean == 1) then
       pLeafM3D => leafm%patch_rough
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'PATCH_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%patch_wetind
    if (imean == 1) then
       pLeafM3D => leafm%patch_wetind
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'PATCH_WETIND :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%leaf_class
    if (imean == 1) then
       pLeafM3D => leafm%leaf_class
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'LEAF_CLASS :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%soil_rough
    if (imean == 1) then
       pLeafM3D => leafm%soil_rough
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'SOIL_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%sfcwater_nlev
    if (imean == 1) then
       pLeafM3D => leafm%sfcwater_nlev
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'SFCWATER_NLEV :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%stom_resist
    if (imean == 1) then
       pLeafM3D => leafm%stom_resist
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'STOM_RESIST :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%ground_rsat
    if (imean == 1) then
       pLeafM3D => leafm%ground_rsat
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'GROUND_RSAT :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%ground_rvap
    if (imean == 1) then
       pLeafM3D => leafm%ground_rvap
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'GROUND_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_water
    if (imean == 1) then
       pLeafM3D => leafm%veg_water
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_WATER :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_temp
    if (imean == 1) then
       pLeafM3D => leafm%veg_temp
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_TEMP :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%can_rvap
    if (imean == 1) then
       pLeafM3D => leafm%can_rvap
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'CAN_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%can_temp
    if (imean == 1) then
       pLeafM3D => leafm%can_temp
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'CAN_TEMP :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_ndvip
    if (imean == 1) then
       pLeafM3D => leafm%veg_ndvip
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_NDVIP :6:hist:mpti', &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_ndvic
    if (imean == 1) then
       pLeafM3D => leafm%veg_ndvic
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_NDVIC :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
         pLeafM3D, imean)

    pLeaf3D => leaf%veg_ndvif
    if (imean == 1) then
       pLeafM3D => leafm%veg_ndvif
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'VEG_NDVIF :6:hist:mpti', &
         pLeafM3D, imean)

    pLeaf3D => leaf%R_aer
    if (imean == 1) then
       pLeafM3D => leafm%R_aer
    else
       pLeafM3D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf3D, &
         'R_aer :6:hist:mpti', &
         pLeafM3D, imean)

    pLeaf2D => leaf%snow_mass
    if (imean == 1) then
       pLeafM2D => leafm%snow_mass
    else
       pLeafM2D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf2D, &
         'SNOW_MASS :2:mpti', &
         pLeafM2D, imean)

    pLeaf2D => leaf%snow_depth
    if (imean == 1) then
       pLeafM2D => leafm%snow_depth
    else
       pLeafM2D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf2D, &
         'SNOW_DEPTH :2:mpti', &
         pLeafM2D, imean)

    pLeaf2D => leaf%seatp
    if (imean == 1) then
       pLeafM2D => leafm%seatp
    else
       pLeafM2D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf2D, &
         'SEATP :2:mpti', &
         pLeafM2D, imean)

    pLeaf2D => leaf%seatf
    if (imean == 1) then
       pLeafM2D => leafm%seatf
    else
       pLeafM2D => null()
    end if
    call InsertVarTable (oneVarTable, oneVarTableSize, &
         pLeaf2D, &
         'SEATF :2:mpti', &
         pLeafM2D, imean)

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
