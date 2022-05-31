module ChemDryDepDriver


  use rconstants, only: &
       cpi,             &
       cpor,            &
       p00,             &
       g,               &
       vonk

  use mem_grid, only: &
       grid_g,        &
       jdim,          &
       dzt,           &
       zt,            &
       nzpmax,        &
       npatch,        &
       dtlt,          &
       imonth1,       &
       idate1,        &
       iyear1,        &
       ngrid

  use micphys, only: &
       level

  use mem_cuparm, only: &
       cuparm_g,        &
       nnqparm

  use ModBasicFields, only: &
       BasicFields

  use ModTurbFields, only: &
       TurbFields

  use mem_leaf, only: &
       leaf_g

  use mem_micro, only: &
       micro_g

  use mem_radiate, only: &
       radiate_g

  use mem_chem1, only: &
       chem1_g,        &
       chemistry

  use mem_aer1, only: &
       aerosol, &
       aer1_g

  use module_dry_dep, only: &
       dd_sedim,            &
       dry_dep                 ! Subroutine



  implicit none

  private


  public :: drydep_driver



contains


  !========================================================================
  subroutine drydep_driver(m1,m2,m3,ia,iz,ja,jz, &
       oneBasicFields, oneTurbFields)

    integer,              intent(IN)    :: m1
    integer,              intent(IN)    :: m2
    integer,              intent(IN)    :: m3
    integer,              intent(IN)    :: ia
    integer,              intent(IN)    :: iz
    integer,              intent(IN)    :: ja
    integer,              intent(IN)    :: jz
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields

    real, pointer :: conprr_dummy(:,:)    


    if (associated(cuparm_g(ngrid)%conprr)) then
       conprr_dummy=>cuparm_g(ngrid)%conprr
    else
       allocate(conprr_dummy(m2,m3))
       nullify(conprr_dummy)
    endif

    call dry_dep(m1,m2,m3,ia,iz,ja,jz           & 
         ,cpi,cpor,p00,g,vonk            &
         ,jdim,dzt,zt,nzpmax,npatch,dtlt &
         ,level,nnqparm(ngrid)           &
         ,imonth1,idate1,iyear1          &
         ,chemistry,aerosol              &
         ,oneBasicFields%theta        &
         ,oneBasicFields%rv        &
         ,oneBasicFields%pp        &
         ,oneBasicFields%dn0        &
         ,oneBasicFields%pi0        &
         ,oneBasicFields%up              &
         ,oneBasicFields%vp           &
         ,oneTurbFields%tkep        &
         ,oneTurbFields%sflux_t          &
         ,oneTurbFields%sflux_r          &
         ,oneTurbFields%sflux_u          &
         ,oneTurbFields%sflux_v          &
         ,leaf_g(ngrid)%r_aer            &
         ,leaf_g(ngrid)%ustar        &
         ,leaf_g(ngrid)%tstar        &
         ,leaf_g(ngrid)%patch_area       &
         ,leaf_g(ngrid)%leaf_class       &
         ,leaf_g(ngrid)%patch_rough      & 
         ,micro_g(ngrid)%rcp        & 
         ,micro_g(ngrid)%pcpg            & 
         ,grid_g(ngrid)%rtgt        & 
         ,radiate_g(ngrid)%rshort        & 
         !                ,cuparm_g(ngrid)%conprr         & 
         ,conprr_dummy                   & 
         !-srf-27jan2015
         !-    changed (:,:,ngrid) to (:,:,:) to avoid
         !-    seg violation when AEROSOL=0
         !               ,aer1_g(:,:,ngrid)             & 
         !               ,chem1_g(:,ngrid)              & 
         !               ,dd_sedim(:,ngrid))    
         ,aer1_g  (:,:,:)              & 
         ,chem1_g (:,:)                & 
         ,dd_sedim(:,:)                )    

    return
  end subroutine drydep_driver
  !========================================================================


end module ChemDryDepDriver
