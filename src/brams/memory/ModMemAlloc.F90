!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModMemAlloc

  use ModLeafComs, only: &
       alloc_leafcol

  use ModCuParmFields, only: &
       InsertCuParmFieldsAtVarTable, &
       InsertCuParmShFieldsAtVarTable
  
  use mem_grid, only: &
       nxtnest, &
       oneGlobalGridData, &
       grid_g, &
       gridm_g, &
       naddsc, &
       maxnxp, &
       maxnyp, &
       maxnzp, &
       if_adap, &
       npatch, &
       nnxp, &
       nnyp, &
       nnzp, &
       nzg,  &
       nzs,  &
       ngrids, &
       nullify_grid, &
       alloc_grid, &
       nullify_GlobalGridData, &
       alloc_GlobalGridData, &
       dtlong, &
       nndtrat, &
       filltab_grid, &
       timmax, &
       dealloc_grid

  use mem_leaf, only: &
       leaf_g, &
       leafm_g, &
       nullify_leaf, &
       alloc_leaf, &
       filltab_leaf,&
       isfcl, &
       dealloc_leaf

  use ModMicroFields, only: &
       InsertMicroFieldsAtVarTable

  use mem_oda, only: &
       oda_g, &
       odam_g, &
       nullify_oda, &
       alloc_oda, &
       filltab_oda, &
       dealloc_oda

  use ModRadiateFields, only: &
       InsertRadiateFieldsAtVarTable
  
  use ModScalarFields, only: &
       InsertScalarFieldsAtVarTable 

  use mem_varinit, only: &
       varinit_g, &
       nullify_varinit, &
       alloc_varinit, &
       filltab_varinit, &
       dealloc_varinit

  use ModVarTable, only: &
       MarkLiteVarsAtVarTable

#ifdef JULES
  use ModJulesFields, only: &
       InsertJulesFieldsAtVarTable
#endif

  use ModShcuFields, only: &
       InsertShcuFieldsAtVarTable

  use mem_opt, only: &
       nullify_opt_scratch, &
       alloc_opt_scratch, &
       dealloc_opt_scratch

  use mem_aerad, only: &
       nwave,          &         !INTENT(IN)
       initial_definitions_aerad, & !Subroutine
       final_definitions_aerad   !Subroutine

  use mem_globaer, only: &
       initial_definitions_globaer, & !Subroutine
       final_definitions_globaer !Subroutine

  use mem_globrad, only : &
       ntotal, &
       nlayer, &
       initial_definitions_globrad, & !Subroutine
       final_definitions_globrad !Subroutine

  use mem_scratch, only: &
       alloc_scratch,    &
       nullify_scratch,  &
       createvctr, &
       dealloc_scratch

  use mem_tend, only: &
       nullify_tend, &
       alloc_tend, &
       filltab_tend, &
       dealloc_tend

  use ModParallelEnvironment, only: &
       MsgDump

  use ModGrid, only: &
       Grid

  use grid_dims, only: &
       maxmach, &
       maxsclr, &
       maxgrds

  use mem_scratch1, only: &
       alloc_scratch1, &
       nullify_scratch1

  use mem_nestb, only: &
       alloc_nestb

  use io_params, only: &
       avgtim, &
       frqanl

  use ModBasicFields, only: &
       InsertBasicFieldsAtVarTable

  use ModTurbFields, only: &
       InsertTurbFieldsAtVarTable

  use node_mod, only: &
       alloc_paths,   & !Subroutine
       nmachs,        &
       ipaths,        &
       iget_paths,    &
       nodemxp,       &
       nodemyp,       &
       nodebounds,    &
       mchnum,         &
       master_num,     &
       mynum

  use shcu_vars_const, only : &
       nnshcu                           ! INTENT(IN)

  use mem_grell_param, only: &
       ngrids_cp,            &
       flag_grell,           &
       closure_type,         &
       icoic,                &
       icoic_sh,             &
       define_memory

  use mem_grell, only: &
       grell_g,        &
       grellm_g,       &
       grell_g_sh,     &
       grellm_g_sh,    &
       alloc_grell,    &
       nullify_grell,  &
       filltab_grell,  &
       filltab_grell_sh, &
       cuforc_g, &
       cuforc_sh_g, &
       cuforcm_g, &
       cuforcm_sh_g, &
       nullify_cuforc, &
       alloc_cu_forcings, &
       filltab_cuforc_sh, &
       filltab_cuforc, &
       alloc_grell_sh

  use mem_scratch1_grell, only: &
       alloc_scratch1_grell

  use mem_carma, only: &
       carma,          &
       carma_m,        &
       nullify_carma,  &
       filltab_carma,  &
       alloc_carma,    &
       zero_carma,     &
       filltab_aotMap, &
       alloc_aotMap,   &
       carma_aotMap, &
       carma_aotMapm,  &
       nullify_aotMap

  use Extras, only: &
       extra2d,     &
       extra3d,     &
       extra2dm,    &
       extra3dm,    &
       na_extra2d,  &
       na_extra3d,  &
       nullify_extra2d, &
       alloc_extra2d,   &
       zero_extra2d,    &
       filltab_extra2d, &
       nullify_extra3d, &
       alloc_extra3d,   &
       filltab_extra3d, &
       zero_extra3d

  use mem_turb_scalar, only: &
       turb_s,         &
       turbm_s,        &
       nullify_turb_s, &
       alloc_turb_s,   &
       filltab_turb_s

  use machine_arq, only: &
       machine ! INTENT(IN)

  use mem_grid_dim_defs, only: &
       define_grid_dim_pointer ! subroutine

  use ModCuParGrell3, only: &
       alloc_grell3, &
       nullify_grell3, &
       filltab_grell3, &
       g3d_ens_g, &
       g3d_ensm_g, &
       g3d_g, &
       g3dm_g, &
       train_dim

  use chem1_list, only:  &
       PhotojMethod,    &       ! INTENT(IN)
       nspecies_chem=> nspecies ! INTENT(IN)

  use aer1_list, only: &        ! INTENT(IN)
       on,             &        ! INTENT(IN)
       mode_alloc,     &        ! INTENT(IN)
       nmodes,         &
       nspecies_aer=> nspecies , & ! INTENT(IN)
       spc_name, &
       numb_alloc, & ! INTENT(IN)
       aerosol_mechanism, NINORG ! INTENT(IN)

  use chem1aq_list, only: &
       nspeciesaq_chem=> nspeciesaq ! INTENT(IN)

  use ccatt_start, only: &
       ccatt                ! CHEM_RAMSIN

  use mem_chem1, only:     &
       nullify_chem1,      & ! Subroutine
       alloc_chem1,        & ! Subroutine
       filltab_chem1,      & ! Subroutine
       define_n_dyn_chem,  & ! Subroutine
       nullify_tend_chem1, & ! Subroutine
       alloc_tend_chem1,   & ! Subroutine
       filltab_tend_chem1, & ! Subroutine
       nsrc,               & ! INTENT(IN)
       max_ntimes_src,     & ! INTENT(IN)
       chemistry,          & ! CHEM_RAMSIN
       chem1_g,            &
       chem1m_g,           &
       chem1_src_g,        &
       chem1m_src_g,       &
       chem1_src_z_dim_g, &
       define_chem1_src_zdim

  use mem_chemic, only: &
       nullify_chemic,  & ! Subroutine
       alloc_chemic,    & ! Subroutine
       chemic_g

  use mem_chem1aq, only:     &
       nullify_chem1aq,      & ! Subroutine
       alloc_chem1aq,        & ! Subroutine
       filltab_chem1aq,      & ! Subroutine
       nullify_tend_chem1aq, & ! Subroutine
       alloc_tend_chem1aq,   & ! Subroutine
       filltab_tend_chem1aq, & ! Subroutine
       chemistry_aq,         & ! CHEM_RAMSIN
       chem1aq_g,            &
       chem1maq_g

  use mem_aer1, only:        &
       nullify_aer1,         & ! Subroutine
       nullify_aer2,         & !
       alloc_aer1,           & ! Subroutine
       alloc_aer2,           & ! Subroutine
       filltab_aer1,         & ! Subroutine
       filltab_aer2,         & ! Subroutine
       nullify_tend_aer1,    & ! Subroutine
       nullify_tend_aer2,    & ! Subroutine
       alloc_tend_aer1,      & ! Subroutine
       alloc_tend_aer2,      & ! Subroutine
       filltab_tend_aer1,    & ! Subroutine
       filltab_tend_aer2,    & ! Subroutine
       aer1_g,               &
       aer1m_g,              &
       aer1_src_z_dim_g,     &
       aer2_src_z_dim_g,     &
       aerosol         ,     &
       aer2_g          ,     &
       aer2m_g         ,     &
       nullify_aer1_inorg,         & ! Subroutine
       alloc_aer1_inorg,           & ! Subroutine
       filltab_aer1_inorg,         & ! Subroutine
       nullify_tend_aer1_inorg,    & ! Subroutine
       alloc_tend_aer1_inorg,      & ! Subroutine
       filltab_tend_aer1_inorg,    & ! Subroutine
       aer1_inorg_g,               &
       aer1m_inorg_g,              &
       define_aer1_src_zdim

  use mem_plume_chem1, only: &
       nullify_plume_chem1,  & ! Subroutine
       alloc_plume_chem1,    & ! Subroutine
       filltab_plume_chem1,  & ! Subroutine
       nveg_agreg,           & ! INTENT(IN)
       plumerise,            & ! CHEM_RAMSIN
       plume_g,              &
       plumem_g,             &
       plume_mean_g,         &
       plume_meanm_g,        &
       plume_fre_g,          &
       plumem_fre_g

  use mem_volc_chem1, only: &
       nullify_volc_chem1,  & ! Subroutine
       alloc_volc_chem1,    & ! Subroutine
       filltab_volc_chem1,  & ! Subroutine
       volcanoes,           & ! CHEM_RAMSIN
       volc_mean_g,         &
       volc_meanm_g

  use chem_sources, only:       &
       oneGlobalEmissData,      &
       alloc_GlobalEmissData,   &
       nullify_GlobalEmissData

  use module_dry_dep, only: &
       alloc_aer_sedim,     & ! Subroutine
       alloc_aer_sedim_numb

  use mem_stilt, only: &
       nullify_stilt,  & ! Subroutine
       alloc_stilt,    & ! Subroutine
       filltab_stilt,  & ! Subroutine
       stilt_g,        &
       stiltm_g

  use carma_fastjx, only: &
       alloc_carma_fastjx   ! Subroutine

  use mem_tuv, only : &
       alloc_carma_tuv, &     ! Subroutine
       nullify_carma_tuv, &
       tuv2carma, &
       carma_tuv, &
       carma_tuvm, &
       alloc_tuv_bio, &
       nullify_tuv_bio, &
       filltab_tuv_bio, &
       tuv_bio, &
       tuv_biom, &
       nbio

  use modtuv, only: &
       ks, &
       nw

  use digitalFilter, only:     &
       initDigitalFilter,     & ! subroutine
       dfVars,                & ! intent(out) - initializing
       applyDF

  use ModOptical, only: &
       setOptMemory

  use ModEvaluation, only: &
       allocStatistic, &
       timeCount, &
       nTimes

  use modIau, only : &
       applyIAU         &
       ,tend_iau_g       &
       ,nullifyIAU       &
       ,allocIau

  use ModAero2McphysFields, only: &
       InsertAero2McphysFieldsAtVarTable
  
  use parrrsw, only : &
       nbndsw

  implicit none

  private

  public :: MemAlloc

contains

  subroutine MemAlloc(oneGrid, proc_type)

    ! Argumenst:
    type(Grid), pointer, intent(in) :: oneGrid
    integer, intent(IN) :: proc_type

    ! Local Variables:
    integer, pointer :: nmzp(:), nmxp(:), nmyp(:)
    integer :: ng, nv, imean, na
    ! Local variable to control Shallow Cumulus memory allocation
    integer :: Alloc_ShCu_Flag
    ! Local variable to control Grell memory allocation
    integer :: Alloc_Grell_Flag
    integer :: ng_cp
    ! Flag to control new Grell MEmory allocation
    integer :: Alloc_Grell3_Flag

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(MemAlloc)**"
    integer :: ierr,n
    integer :: ne2d, ne3d, nsa
    integer :: idiffk
    real, pointer :: v_p

    integer, parameter :: maxvars=1000

    call alloc_paths(ngrids, nmachs)

    !-------------
    ! First, depending on type of process, define grid point pointers correctly..
    if (proc_type==0 .or. proc_type==1) then
       !  This is the call for either a single processor run or
       !    for the master process
       nmzp => nnzp
       nmxp => nnxp
       nmyp => nnyp
    elseif (proc_type==2) then
       !  This is the call for a initial compute node process
       nmzp => nnzp
       nmxp => nodemxp(mynum,:)
       nmyp => nodemyp(mynum,:)
    elseif (proc_type==3) then
       !  This is the call for a dynamic balance compute node process
       nmzp => nnzp
       nmxp => nodemxp(mynum,:)
       nmyp => nodemyp(mynum,:)
    endif

    ! Call global grid dimension definitions
    !**(JP)** To be eliminated
    call define_grid_dim_pointer(proc_type, ngrids, maxgrds, &
         nnzp, nnxp, nnyp, nnzp, nodemxp(mynum,:), nodemyp(mynum,:))

    !  If we are doing time-averaging for output, set flag ...
    imean = 0
    if (avgtim/=0.) imean = 1
    !-------------

    !-------------
    ! insert Basic Field variables at var_table
    do ng=1,ngrids
       call InsertBasicFieldsAtVarTable(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            oneGrid%oneBasicFields, oneGrid%oneAveBasicFields, imean)
    enddo
    !
    !Allocate and prepare optical properties memory
    if (oneGrid%oneNamelistFile%iswrtyp==6) then
       call setOptMemory(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            ngrids,imean,nmzp,nmxp,nmyp)
    end if

    !-------------
    ! Allocate Leaf surface scheme type
    allocate(leaf_g(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating leaf_g")
    allocate(leafm_g(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating leafm_g")
    do ng=1,ngrids
       call nullify_leaf(leaf_g(ng))
       call alloc_leaf(leaf_g(ng), nmzp(ng), nmxp(ng), nmyp(ng),  &
            nzg, nzs, npatch, ng)
       call nullify_leaf(leafm_g(ng))
       if (imean==1) then
          call alloc_leaf(leafm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng),  &
               nzg, nzs, npatch, ng)
       endif

       call filltab_leaf(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            leaf_g(ng), leafm_g(ng), imean)
    enddo
    ! Bob (1/10/2002) added the following line.  Is this the right place for
    ! the long term??
    call alloc_leafcol(nzg, nzs)
    !-------------
    !
#ifdef JULES
    ! insert Jules Fields variables at var_table
    call InsertJulesFieldsAtVarTable(&
         oneGrid%oneVarTable, &
         oneGrid%oneVarTableSize, &
         oneGrid%oneJulesFields, &
         oneGrid%oneAveJulesFields, &
         oneGrid%oneNamelistFile, &
         imean)
#endif
    !-------------

    !-------------
    ! insert Micro Fields variables at var_table
    call InsertMicroFieldsAtVarTable(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
         oneGrid%oneMicroFields, oneGrid%oneAveMicroFields, imean)
    !-------------

    !-------------
    ! insert Turb Field variables at var_table
    do ng=1,ngrids
       call InsertTurbFieldsAtVarTable(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            oneGrid%oneTurbFields, oneGrid%oneAveTurbFields, imean)
    enddo

    if (CCATT==1 .and. chemistry >= 0) then
       allocate(turb_s(ngrids), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating turb_s")
       allocate(turbm_s(ngrids), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating turbm_s")
       do ng=1,ngrids
          call nullify_turb_s(turb_s (ng))
          call alloc_turb_s(turb_s(ng), nmzp(ng), nmxp(ng), nmyp(ng))
          call nullify_turb_s(turbm_s (ng))
          if (imean==1) then
             call alloc_turb_s(turbm_s(ng), nmzp(ng), nmxp(ng), nmyp(ng))
          endif
          call filltab_turb_s(&
               oneGrid%oneVarTable, &
               oneGrid%oneVarTableSize, &
               turb_s(ng), &
               turbm_s(ng), imean)
       enddo
    endif
    !-------------

    !-------------
    ! Allocate varinit variables data type.
    ! These do not need "mean" type ever; that's why varinitm_g(ng) does not exist
    allocate(varinit_g(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating varinit_g")
    do ng=1,ngrids
       call nullify_varinit(varinit_g(ng))
       call alloc_varinit(varinit_g(ng), nmzp(ng), nmxp(ng), nmyp(ng))
       call filltab_varinit(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            varinit_g, imean, ng)
    enddo
    !-------------

    !-------------
    ! Allocate grid variables data type.
    allocate(grid_g(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating grid_g")
    allocate(gridm_g(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating gridm_g")
    do ng=1,ngrids
       call nullify_grid(grid_g(ng))
       call alloc_grid(grid_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng, if_adap)
       call nullify_grid(gridm_g(ng))
       if (imean == 1) then
          call alloc_grid(gridm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng, if_adap)
       end if

       call filltab_grid(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            grid_g(ng), gridm_g(ng), imean)
    enddo
    !-------------

    !-------------
    ! Allocate global grid variables data type.
    allocate(oneGlobalGridData(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating oneGlobalGridData")
    do ng=1,ngrids
       call nullify_GlobalGridData(oneGlobalGridData(ng))
       call alloc_GlobalGridData(oneGlobalGridData(ng), nnxp(ng), nnyp(ng))
    end do
    !-------------

    !-------------
    ! insert radiate fields at var table
    do ng=1,ngrids
       call InsertRadiateFieldsAtVarTable(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            oneGrid%oneRadiateFields, oneGrid%oneAveRadiateFields, imean)
    enddo

    !- only for CARMA/RRTM Radiations schems
    if (oneGrid%oneNamelistFile%ilwrtyp==4 .or. oneGrid%oneNamelistFile%iswrtyp==4 .or. oneGrid%oneNamelistFile%ilwrtyp==6 .or. oneGrid%oneNamelistFile%iswrtyp==6 ) then
       call initial_definitions_aerad()
       call initial_definitions_globrad()
       call initial_definitions_globaer()

       if (trim(PhotojMethod) == 'FAST-JX' .and. chemistry >= 0) then
          call alloc_carma_fastjx(maxval(nmzp(1:ngrids)),maxval(nmxp(1:ngrids)),maxval(nmyp(1:ngrids)))
       endif

       if ((trim(PhotojMethod) == 'FAST-TUV' .or. trim(PhotojMethod) == 'LUT').and. chemistry >= 0) then
          allocate(tuv2carma(nw),stat=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating tuv2carma fails")

          !--(BRAMS-5.0-INI)---------------------------------------------------
          tuv2carma = 0
          !--(BRAMS-5.0-FIM)---------------------------------------------------

          allocate(carma_tuv(ngrids), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating carma_tuv")
          allocate(carma_tuvm(ngrids), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating carma_tuvm")

          allocate(tuv_bio(ngrids,nbio),STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating tuv_bio")
          allocate(tuv_biom(ngrids,nbio),STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating tuv_biom")

          do ng=1,ngrids
             call nullify_carma_tuv(carma_tuv(ng))
             call nullify_carma_tuv(carma_tuvm(ng))
             !
             do n=1,nbio
                call nullify_tuv_bio(tuv_bio(ng,n))
                call nullify_tuv_bio(tuv_biom(ng,n))
             end do
             !
             do n=1,nbio
                call alloc_tuv_bio(tuv_bio(ng,n),nmxp(ng), nmyp(ng))
                if (imean == 1) then
                   call alloc_tuv_bio(tuv_biom(ng,n),nmxp(ng), nmyp(ng))
                end if
             end do
             !
             call alloc_carma_tuv(carma_tuv(ng),ntotal,nlayer,nmxp(ng), nmyp(ng), ng)
             call alloc_carma_tuv(carma_tuvm(ng),ntotal,nlayer,nmxp(ng), nmyp(ng), ng)
             !
             do n=1,nbio
                call filltab_tuv_bio(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
                     tuv_bio(ng,n), tuv_bio(ng,n), n, imean)
             end do
          end do
       endif

       allocate(carma_aotMap(ngrids), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating carma_aotMap")
       allocate(carma_aotMapm(ngrids), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating carma_aotMapm")
       if(.not. allocated(carma_aotMap)) allocate(carma_aotMap(ngrids))
       if(.not. allocated(carma_aotMapm)) allocate(carma_aotMapm(ngrids))
       do ng=1,ngrids
          call nullify_aotMap(carma_aotMap,ng)
          call nullify_aotMap(carma_aotMapm,ng)
          call alloc_aotMap(carma_aotMap(ng),nmxp(ng), nmyp(ng))
          if (imean == 1) then
             call alloc_aotMap(carma_aotMapm(ng),nmxp(ng), nmyp(ng))
          end if

          call filltab_aotMap(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
               carma_aotMap, carma_aotMapm, ng, imean)
       end do

       !-only CARMA
       if (oneGrid%oneNamelistFile%ilwrtyp==4 .or. oneGrid%oneNamelistFile%iswrtyp==4 .or. oneGrid%oneNamelistFile%iswrtyp==6) then !Colocando o RTM para alocar AOT

          allocate(carma(ngrids), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating carma/AOT rtm")
          allocate(carma_m(ngrids), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating carma_m")
          !if(oneGrid%oneNamelistFile%iswrtyp==4) then
          do ng=1,ngrids
             call nullify_carma(carma,ng)
             call alloc_carma(carma, ng, nmxp(ng), nmyp(ng), nwave)
             call zero_carma(carma, ng)
             call nullify_carma(carma_m, ng)

             if(imean==1) then
                call alloc_carma(carma_m, ng, nmxp(ng), nmyp(ng), nwave)
                call zero_carma(carma_m, ng)
             end if

             call filltab_carma(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
                  oneGrid%oneNamelistFile, &
                  carma, carma_m, ng, imean)
          end do
       end if
    endif
    !-------------

    !-------------
    ! Allocate cumulus parameterizations variables data type
    !
    ! Initiating values:
    Alloc_ShCu_Flag  = 0
    Alloc_Grell_Flag = 0
    Alloc_Grell3_Flag =0

    do ng=1,ngrids
       call InsertCuParmFieldsAtVarTable(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            oneGrid%oneCuParmFields, oneGrid%oneAveCuParmFields, imean)

       !-srf-feb2012: for shallow cumulus
       if (nnshcu(ng) == 2) then
          call InsertCuParmShFieldsAtVarTable(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
               oneGrid%oneCuParmFields, oneGrid%oneAveCuParmFields, imean)
       end if

    enddo
    !--------------------------------------------------------------------------
    ! insert Shallow Cumulus at var table
    do ng=1, ngrids
       if (nnshcu(ng)==1) then
          call InsertShcuFieldsAtVarTable(&
               oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
               oneGrid%oneShcuFields, oneGrid%oneAveShcuFields, imean)
       end if
    end do
    !--------------------------------------------------------------------------
    ! Allocate data for Shallow/Deep Cumulus options 2 and above
    !
    Alloc_SHCU_Flag   =0
    Alloc_Grell_Flag  =0
    Alloc_Grell3_Flag =0
    do ng=1, ngrids
       if (NNSHCU (ng)== 1 .or. NNSHCU (ng)== 2) Alloc_SHCU_Flag   = 1
       if (oneGrid%oneNamelistFile%nnqparm(ng) >  2) Alloc_Grell3_Flag = 1
    enddo

    if (Alloc_Grell_Flag == 1 .or. Alloc_SHCU_Flag == 1 .or. Alloc_Grell3_Flag == 1) then
       ! Calculating the necessary space for scratch data
       call define_memory(nmxp, nmyp, nmzp, ngrids, oneGrid%oneNamelistFile%nnqparm, nnshcu)
       ! Allocating data for scratch data
       call alloc_scratch1_grell()

       allocate(cuforc_g(ngrids_cp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(2),"(i8)") ngrids_cp
          call fatal_error(h//" allocate cuforc_g("//&
               trim(adjustl(str(2)))//") "//&
               "fails with stat="//trim(adjustl(str(1))))
       end if

       allocate(cuforcm_g(ngrids_cp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(2),"(i8)") ngrids_cp
          call fatal_error(h//" allocate cuforcm_g("//&
               trim(adjustl(str(2)))//") "//&
               "fails with stat="//trim(adjustl(str(1))))
       end if

       allocate(cuforc_sh_g(ngrids_cp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(2),"(i8)") ngrids_cp
          call fatal_error(h//" allocate cuforc_sh_g("//&
               trim(adjustl(str(2)))//") "//&
               "fails with stat="//trim(adjustl(str(1))))
       end if

       allocate(cuforcm_sh_g(ngrids_cp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(2),"(i8)") ngrids_cp
          call fatal_error(h//" allocate cuforcm_sh_g("//&
               trim(adjustl(str(2)))//") "//&
               "fails with stat="//trim(adjustl(str(1))))
       end if

       do ng=1,ngrids_cp
          call nullify_cuforc(cuforc_g(ng))
          call nullify_cuforc(cuforcm_g(ng))
          call nullify_cuforc(cuforc_sh_g(ng))
          call nullify_cuforc(cuforcm_sh_g(ng))
          call alloc_cu_forcings(cuforc_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
          call alloc_cu_forcings(cuforc_sh_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
          if (imean == 1) then
             call alloc_cu_forcings(cuforcm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
             call alloc_cu_forcings(cuforcm_sh_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
          endif

          call filltab_cuforc_sh(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
               cuforc_sh_g(ng), cuforcm_sh_g(ng), imean)

          call filltab_cuforc(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
               cuforc_g(ng), cuforcm_g(ng), imean)
       enddo
    endif

    !- Allocate data for Grell Cumulus version 3d,GD-FIM and GF schemes
    if (Alloc_Grell3_Flag == 1) then
       allocate(g3d_ens_g(train_dim,ngrids_cp))
       allocate(g3d_ensm_g(train_dim,ngrids_cp))
       allocate(g3d_g(ngrids_cp))
       allocate(g3dm_g(ngrids_cp))

       ng_cp = 1
       do ng=1,ngrids
          if (oneGrid%oneNamelistFile%nnqparm(ng) > 2) then
             !-- arrays needed for G3d , GD-FIM and GF schemes

             call nullify_grell3(g3d_ens_g (:,ng_cp) , g3d_g(ng_cp),train_dim)
             call nullify_grell3(g3d_ensm_g(:,ng_cp) ,g3dm_g(ng_cp),train_dim)

             call alloc_grell3(oneGrid%oneNamelistFile, &
                  g3d_ens_g(:,ng_cp),g3d_g(ng_cp),nmzp(ng),nmxp(ng),nmyp(ng),ng,train_dim)

             if (imean == 1) then
                call alloc_grell3(oneGrid%oneNamelistFile, &
                     g3d_ensm_g(:,ng_cp),g3dm_g(ng_cp),nmzp(ng),nmxp(ng),nmyp(ng),ng,train_dim)
             endif

             call filltab_grell3(&
                  oneGrid%oneVarTable, &
                  oneGrid%oneVarTableSize, &
                  g3d_ens_g, &
                  g3d_g,&
                  g3d_ensm_g, &
                  g3dm_g,&
                  train_dim, oneGrid%oneNamelistFile%nnqparm(ng), ng_cp, imean)

             ng_cp = ng_cp + 1
             if  (CLOSURE_TYPE == 'EN') then
                icoic = 0
             elseif (CLOSURE_TYPE == 'GR') then
                icoic = 2
             elseif (CLOSURE_TYPE == 'LO') then
                icoic = 5
             elseif (CLOSURE_TYPE == 'MC') then
                icoic = 8
             elseif (CLOSURE_TYPE == 'SC') then
                icoic = 11
                !elseif (CLOSURE_TYPE == 'AS') then
                !  stop 'CLOSURE_TYPE == AS is not allowed for GF scheme'
             elseif (CLOSURE_TYPE == 'PB' ) then
                icoic = 13
             else
                print *, "****Grell Closure type ERROR for GF scheme",CLOSURE_TYPE
                ! the subroutine opspec3 stop the program before this point.
             endif
          endif
          icoic_sh=0
       enddo
    endif
    !-------------

    !-------------
    !-srf tmp - not using ODA features for now
    ! Allocate oda variables data type.
    !    These do not need "mean" type ever.
    allocate(oda_g(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating oda_g")
    allocate(odam_g(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating odam_g")

    !-------------
    ! Allocate any added Scalar types
    ! NOT ALLOWING DIFFERENT NUMBERS OF SCALARS ON DIFFERENT NESTS
    !   Allocate length 1 of these datatypes by default
    if (naddsc>0) then
       call InsertScalarFieldsAtVarTable(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            oneGrid%oneScalarFields, oneGrid%oneAveScalarFields, &
            oneGrid%oneNamelistFile, imean)
    end if

    !-------------

    !-------------
    ! Allocate Tendency data type,  filltab_tendency is responsible
    !   for filling the main atmospheric model variables in the scalar table,
    !   so make sure to call any routines that define scalar variables first.
    ! Assuming same scalars on all grids!!!!!

    call nullify_tend(naddsc)

    call alloc_tend(nmzp, nmxp, nmyp, ngrids, naddsc, proc_type, &
         oneGrid%oneBasicFields, oneGrid%oneTurbFields, oneGrid%oneMicroFields)
    !-------------

    !-------------
    ! - CCATT chemistry/aerosol modules
    if (ccatt == 1  .and. chemistry >= 0) then

       ! Allocate global grid variables data type.
       allocate(oneGlobalEmissData(ngrids), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating oneGlobalEmissData")
       do ng=1,ngrids
          call nullify_GlobalEmissData(oneGlobalEmissData(ng))
          call alloc_GlobalEmissData(oneGlobalEmissData(ng), nnxp(ng), nnyp(ng))
       end do


       call define_n_dyn_chem(ngrids,dtlong,nndtrat,mynum)

       allocate(chem1_g(nspecies_chem,ngrids), chem1m_g(nspecies_chem,ngrids))
       allocate(chem1_src_g (max_ntimes_src,nsrc,nspecies_chem,ngrids),&
            chem1m_src_g(max_ntimes_src,nsrc,nspecies_chem,ngrids))

       do ng=1,ngrids

          call define_chem1_src_zdim(chem1_src_z_dim_g(:,ng),nmzp(ng))
          call nullify_chem1(chem1_g (:,ng),chem1_src_g (:,:,:,ng), nspecies_chem)
          call nullify_chem1(chem1m_g(:,ng),chem1m_src_g(:,:,:,ng), nspecies_chem)

          call alloc_chem1(chem1_g(:,ng),chem1_src_g(:,:,:,ng),chem1_src_z_dim_g(:,ng) &
               ,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem, ng,volcanoes)

          if (imean == 1) then
             call alloc_chem1(chem1m_g(:,ng),chem1m_src_g(:,:,:,ng),chem1_src_z_dim_g(:,ng) &
                  ,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem, ng,volcanoes)
          endif

          call filltab_chem1(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
               chem1_g(:,ng) ,chem1m_g(:,ng),  &
               chem1_src_g (:,:,:,ng),chem1m_src_g(:,:,:,ng),  &
               chem1_src_z_dim_g(:,ng), &
               nmzp(ng),nspecies_chem,volcanoes, imean)
       end do

       call nullify_tend_chem1(nspecies_chem)

       call alloc_tend_chem1(nmzp,nmxp,nmyp,ngrids,nspecies_chem,proc_type)

       if( plumerise > 0) then
          !-- plumerise  section
          allocate(plume_g     (nveg_agreg,nspecies_chem,ngrids),plumem_g     (nveg_agreg,nspecies_chem,ngrids))
          allocate(plume_mean_g(nveg_agreg,ngrids)      ,plume_meanm_g(nveg_agreg,ngrids))
          !- this is for FRP methodology
          allocate(plume_fre_g (5,ngrids)                       ,plumem_fre_g (5,ngrids))

          do ng=1,ngrids

             call nullify_plume_chem1(plume_g  (:,:,ng),plume_mean_g  (:,ng),plume_fre_g   (:,ng),nspecies_chem)
             call nullify_plume_chem1(plumem_g (:,:,ng),plume_meanm_g (:,ng),plumem_fre_g  (:,ng),nspecies_chem)

             call alloc_plume_chem1(plume_g(:,:,ng),plume_mean_g  (:,ng),plume_fre_g  (:,ng)&
                  ,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem)
             if (imean == 1) then
                call alloc_plume_chem1(plumem_g(:,:,ng),plume_meanm_g (:,ng),plumem_fre_g  (:,ng)&
                     ,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem)
             endif

             call filltab_plume_chem1(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
                  plume_g(:,:,ng), plumem_g(:,:,ng),  &
                  plume_mean_g(:,ng), plume_meanm_g(:,ng),    &
                  plume_fre_g(:,ng), plumem_fre_g(:,ng), nspecies_chem, ng, imean)
          enddo
       endif
       !-- aerosol section----------------------------------
       if( aerosol >= 1) then

          !--- allocation of mass concentration

          allocate(aer1_g(nmodes,nspecies_aer,ngrids),aer1m_g(nmodes,nspecies_aer,ngrids))

          call alloc_aer_sedim(npatch,ngrids, &
               nmodes,nspecies_aer,mode_alloc, &
               on,nmzp,nmxp,nmyp)

          do ng=1,ngrids
             call define_aer1_src_zdim(aer1_src_z_dim_g(:,ng),nmzp(ng))

             call nullify_aer1(aer1_g (:,:,ng),nmodes, nspecies_aer)
             call nullify_aer1(aer1m_g(:,:,ng),nmodes, nspecies_aer)

             call alloc_aer1(aer1_g(:,:,ng),aer1_src_z_dim_g(:,ng) &
                  ,nmzp(ng),nmxp(ng),nmyp(ng),nmodes,nspecies_aer)

             if (imean == 1) then
                call alloc_aer1(aer1m_g(:,:,ng),aer1_src_z_dim_g(:,ng) &
                     ,nmzp(ng),nmxp(ng),nmyp(ng),nmodes,nspecies_aer)
             endif

             call filltab_aer1(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
                  aer1_g(:,:,ng), aer1m_g(:,:,ng), aer1_src_z_dim_g(:,ng), imean)
          enddo

          call nullify_tend_aer1(nmodes,nspecies_aer)

          call alloc_tend_aer1(nmzp,nmxp,nmyp,ngrids,nmodes,nspecies_aer,proc_type)
          !
          !----- for Matrix only  --------------------
          if(aerosol==2 .and. aerosol_mechanism(1:6)=='MATRIX') then

             !--- allocation inorganics mass concentration
             allocate(aer1_inorg_g(ninorg,ngrids),aer1m_inorg_g(ninorg,ngrids))

             do ng=1,ngrids

                call nullify_aer1_inorg(aer1_inorg_g (:,ng),ninorg)
                call nullify_aer1_inorg(aer1m_inorg_g(:,ng),ninorg)

                call alloc_aer1_inorg(aer1_inorg_g(:,ng) &
                     ,nmzp(ng),nmxp(ng),nmyp(ng),ninorg)

                if (imean == 1) then
                   call alloc_aer1_inorg(aer1m_inorg_g(:,ng),nmzp(ng),nmxp(ng),nmyp(ng),ninorg)
                endif

                call filltab_aer1_inorg(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
                     aer1_inorg_g(:,ng),aer1m_inorg_g(:,ng), imean)

             enddo

             call nullify_tend_aer1_inorg(ninorg)

             call alloc_tend_aer1_inorg(nmzp,nmxp,nmyp,ngrids,ninorg,proc_type)

             !--- allocation for number concentration
             allocate(aer2_g(nmodes,ngrids),aer2m_g(nmodes,ngrids))
             !
             !- sedimentation is done inside matrix module
             !
             !CALL alloc_aer_sedim_numb(npatch,ngrids, &
             !                           nmodes,numb_alloc,on,nmzp,nmxp,nmyp)

             !--- allocation for aerosol to microphysics arrays
             !
             do ng=1,ngrids

                !-for now, the vertical dimentions of the source array will be nz
                aer2_src_z_dim_g(:,ng)=nmzp(ng)

                call nullify_aer2(aer2_g(:,ng),nmodes)

                call alloc_aer2(aer2_g(:,ng),nmzp(ng),nmxp(ng),nmyp(ng),nmodes)

                call nullify_aer2(aer2m_g(:,ng),nmodes)

                if (imean == 1) then
                   call alloc_aer2(aer2m_g(:,ng),nmzp(ng),nmxp(ng),nmyp(ng),nmodes)
                endif

                call filltab_aer2(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
                     aer2_g(:,ng), aer2m_g(:,ng), aer2_src_z_dim_g(:,ng), imean)

                call InsertAero2McphysFieldsAtVarTable(&
                     oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
                     oneGrid%oneAero2McphysFields, &
                     oneGrid%oneAveAero2McphysFields, imean)
             enddo

             call nullify_tend_aer2(nmodes)

             call alloc_tend_aer2(nmzp,nmxp,nmyp,ngrids,nmodes,proc_type)

          endif ! - for MATRIX only

       endif
       !-- end of aerosol section ----------------------------------

       !-- volcanoes section
       allocate(volc_mean_g(ngrids),volc_meanm_g(ngrids))
       if( volcanoes == 1) then
          !         allocate(volc_mean_g(ngrids),volc_meanm_g(ngrids))

          do ng=1,ngrids

             call nullify_volc_chem1(volc_mean_g(ng))
             call nullify_volc_chem1(volc_meanm_g(ng))

             call alloc_volc_chem1(volc_mean_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))

             if (imean == 1) then
                call alloc_volc_chem1(volc_meanm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))
             endif

             call filltab_volc_chem1(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
                  volc_mean_g(ng), volc_meanm_g(ng), imean)
          enddo
       endif
       !- end of volcanoes section ---
       !
       !-- CHEM1AQ section
       if (chemistry_aq > 0) then
          allocate(chem1aq_g  (nspeciesaq_chem,ngrids),chem1maq_g(nspeciesaq_chem,ngrids))
          allocate(chemic_g(ngrids))

          do ng=1,ngrids

             call nullify_chem1aq(chem1aq_g (:,ng), nspeciesaq_chem)
             call nullify_chem1aq(chem1maq_g(:,ng), nspeciesaq_chem)
             call nullify_chemic(chemic_g(ng))

             call alloc_chem1aq(chem1aq_g(:,ng) &
                  ,nmzp(ng),nmxp(ng),nmyp(ng),nspeciesaq_chem)

             call alloc_chemic(chemic_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),oneGrid%oneMicVars)

             if (imean == 1) then
                call alloc_chem1aq(chem1maq_g(:,ng) &
                     ,nmzp(ng),nmxp(ng),nmyp(ng),nspeciesaq_chem)
             endif

             call filltab_chem1aq(&
                  oneGrid%oneVarTable, &
                  oneGrid%oneVarTableSize, &
                  chem1aq_g(:,ng), &
                  chem1maq_g(:,ng), imean)
          enddo

          call nullify_tend_chem1aq(nspeciesaq_chem)

          call alloc_tend_chem1aq(nmzp,nmxp,nmyp,ngrids,nspeciesaq_chem,proc_type)
       endif
    endif
    !--------------

    !--------------
    ! filltab tendencies for TEB and CCATT submodels
    do ng=1,ngrids

       call filltab_tend(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, &
            oneGrid%oneBasicFields, oneGrid%oneMicroFields, oneGrid%oneTurbFields,  &
            oneGrid%oneScalarFields, naddsc, ng)

       if (ccatt == 1  .and. chemistry >= 0)  then

          call filltab_tend_chem1(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, &
               nspecies_chem, ng)

          !change MP ---chem1aq
          if(chemistry_aq >= 1) call filltab_tend_chem1aq(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, &
               nspeciesaq_chem, ng)
          !end change MP --chem1aq- END

          if (aerosol >= 1)  then
             call filltab_tend_aer1(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, &
                  nmodes,nspecies_aer,ng)
          end if
          if (aerosol == 2)  then
             call filltab_tend_aer1_inorg(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, &
                  ninorg,ng)
             call filltab_tend_aer2      (oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, &
                  nmodes,ng)
          endif

       endif

    enddo
    !--------------

    !--------------
    ! Allocate Scratch data type, This also fills the max's that are needed
    !    by nesting stuff.
    call createVctr(ngrids, nodemxp(mynum,:), nodemyp(mynum,:), nnzp)

    call nullify_scratch()

    call alloc_scratch(nmzp, nmxp, nmyp, nnzp, nnxp, nnyp, maxgrds, ngrids,  &
         nzg, nzs, npatch, proc_type, maxnxp, maxnyp, maxnzp, &
         oneGrid%oneNamelistFile)
    ! Reproducibility - Saulo Barros
    call nullify_scratch1()
    call alloc_scratch1(oneGrid%oneScalarTableSize, &
         nodebounds, maxgrds, ngrids, nnzp, mynum)
    ! For optmization - ALF
    call nullify_opt_scratch()
    if ((if_adap==0) .and. (OneGrid%oneNamelistFile%ihorgrad==2)) &
         call alloc_opt_scratch(proc_type, ngrids, nnzp, nnxp, nnyp, 1000, 1000)

    ! Allocate nested boundary interpolation arrays. All grids will be allocated.
    ! Changed by Alvaro L.Fazenda
    ! To correct a problem when running in a NEC SX-6
    ! Master process needs allocation for nesting in a parallel run
    if (proc_type==0 .or. proc_type==2 .or. proc_type==1) then
       do ng=1,ngrids
          if (nxtnest(ng)==0 ) then
             call alloc_nestb(oneGrid%oneScalarTableSize, ng,        1,        1,        1)
          else
             call alloc_nestb(oneGrid%oneScalarTableSize, ng, nnxp(ng), nnyp(ng), nnzp(ng))
          endif
       enddo
    endif
    !--------------

    !--------------
    !Allocate stilt variables data type
    allocate(stilt_g(ngrids), stiltm_g(ngrids))
    do ng=1, ngrids
       idiffk = OneGrid%oneNamelistFile%idiffk(ng)
       call nullify_stilt(stilt_g(ng))
       call alloc_stilt(idiffk,stilt_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)

       call nullify_stilt(stiltm_g(ng))       
       if (imean == 1) then
          call alloc_stilt(idiffk,stiltm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
       endif

       call filltab_stilt(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
            stilt_g(ng), stiltm_g(ng), imean)
    end do
    !--------------

    !--------------
    !
    allocate(extra2d(na_extra2d,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating extra2d")
    allocate(extra3d(na_extra3d,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating extra3d")
    allocate(extra2dm(na_extra2d,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating extra2dm")
    allocate(extra3dm(na_extra3d,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating extra3dm")
    call nullify_extra2d(extra2d,  na_extra2d, ngrids)
    call nullify_extra2d(extra2dm, na_extra2d, ngrids)
    call nullify_extra3d(extra3d,  na_extra3d, ngrids)
    call nullify_extra3d(extra3dm, na_extra3d, ngrids)
    ! Allocating itens of the structure
    do ng=1,ngrids
       call alloc_extra2d(extra2d, nmxp(ng), nmyp(ng), na_extra2d, ng)
       call zero_extra2d(extra2d, na_extra2d, ng)
       if (imean==1) then
          call alloc_extra2d(extra2dm, nmxp(ng), nmyp(ng), na_extra2d, ng)
          call zero_extra2d(extra2dm, na_extra2d, ng)
       end if
       call alloc_extra3d(extra3d, nmzp(ng), nmxp(ng), nmyp(ng), na_extra3d, ng)
       call zero_extra3d(extra3d, na_extra3d, ng)
       if (imean==1) then
          call alloc_extra3d(extra3dm,  &
               nmzp(ng), nmxp(ng), nmyp(ng), na_extra3d, ng)
          call zero_extra3d(extra3dm, na_extra3d, ng)
       end if
    end do

    do ng=1,ngrids
       do na=1,na_extra2d
          call filltab_extra2d(&
               oneGrid%oneVarTable, &
               oneGrid%oneVarTableSize, &
               extra2d(na,ng), &
               extra2dm(na,ng), &
               na, imean)
       end do
       do na=1,na_extra3d
          call filltab_extra3d(&
               oneGrid%oneVarTable, &
               oneGrid%oneVarTableSize, &
               extra3d(na,ng), &
               extra3dm(na,ng), &
               na, imean)
       end do
    end do
    !--------------

    !-------------
    ! allocate date for digital filter - rmf
    if(applyDF) call initDigitalFilter(dfVars, ngrids, nmzp, nmxp, nmyp)
    !-------------


    !-------------
    ! allocate data for IAU procedure 
    if(applyIAU > 0 ) then
       allocate(tend_iau_g(ngrids), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating tend_iau_g")
       do ng=1,ngrids
          call nullifyIau(tend_iau_g(ng))
          call   allocIau(tend_iau_g(ng), nmzp(ng), nmxp(ng), nmyp(ng))
       enddo
    endif
    !-------------
    ! Set "Lite" variable flags according to namelist input LITE_VARS.
    call MarkLiteVarsAtVarTable(oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
         oneGrid%oneParallelEnvironment)


    !Calling the statistic for allocate the amount of timesteps
    nTimes=int(timmax/dtlong)
    call allocStatistic(nTimes)
    !Initializa time Count
    timeCount=0

  end subroutine MemAlloc
end module ModMemAlloc
