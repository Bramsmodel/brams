
!############################# Change Log ##################################
! 5.2.0
!
!MB: new file for the Runge-Kutta dynamical core
!    (after Wicker, Skamarock, 2002, MWR)
!    Saulo Freitas (INPE), Michael Baldauf (DWD)
!
!############################# Change Log ##################################

module ModTimestepRK

#ifdef JULES
  use ModSfcLyrJules, only: &
       sfclyr_jules
#endif
  
  use ModNudAnalysis, only: &
       datassim
  
  use ModOdaNudge, only: &
       oda_nudge

  use ModSeaSalt, only: &
       SeaSaltDriver

  use ModRShCuPar, only: &
       shcupa

  use ModRConv, only: &
       cuparm

  use ModMicThompsonDriver, only: &
       micro_thompson

  use ModMicGfdlDriver, only: &
       micro_gfdl

  use ModMicrophysicsMisc, only: &
       negadj1

  use ModMicrophysicsDrive, only: &
       micro

  use ModRstilt, only: &
       prep_advflx_to_stilt

  use ModUrbanCanopy, only: &
       urban_canopy

  use ModLeaf3, only: &
       sfclyr

  use ModOzone, only: &
       ozone

  use ModGasPart, only: &
       le_fontes, &
       sources_teb

  use ModTimestep, only: &
       w_damping

  use ModCoriolis, only: &
       corlos

  use ModRThrm, only: &
       thermo, &
       thermo_boundary_driver, &
       theta_thp_rk

  use ModRexev, only: &
       exevolve, &
       get_true_air_density

  use ModCuParGrell3, only: &
       prepare_lsf, &
       cuparm_grell3_catt, &  ! subroutine
       g3d_g


  use ModDiffuse, only: &
       diffuse_brams31

  use ModTurbK, only: &
       diffuse

  use ModRtimi, only: &
       tend0, &
       predtr

  use ModRadvc, only: &
       advectc

  use ModRbnd, only: latbnd, &
       vpsets, &
       rayft,  &
       trsets

  use ModTurbFields, only: &
       DeepCopyToTurbFields,&
       DeepCopyFromTurbFields
  
  use grid_dims, only: &
       nzpmax

  use ModParallelEnvironment, only: &
       MsgDump

  use ModMessageSet, only: &
       PostSendRecvMsgs, &
       WaitSendRecvMsgs

  use ModAcoust, only:         &
       acoustic_new,            &
       init_div_damping_coeff,  &
       deallocate_alpha_div,    &
       apply_div_damping, &
       buoyancy

  use ModGrid, only: &
       Grid

  use node_mod, only: &
       mzp, mxp, myp,  & ! INTENT(IN)
       ia, iz, ja, jz, & ! INTENT(IN)
       i0, j0,         & ! INTENT(IN)
       izu, jzv,       & ! INTENT(IN)
       mynum,          & ! INTENT(IN)
       ibcon,          & ! INTENT(IN)
       nmachs, &
       mchnum,        &
       master_num, &
       nodemyp,  &  !intent(in)
       nodemxp,  &  !intent(in)
       nodemzp      !intent(in)

  use mem_cuparm, only: &
       NNQPARM, & ! INTENT(IN)
       IF_CUINV   ! INTENT(IN)

  use mem_varinit, only: &
       NUD_TYPE   &! INTENT(IN)
       ,varinit_g  ! INTENT(IN)

  use mem_oda,   only: &
       if_oda ! INTENT(IN)

  use micphys,   only: &
       mcphys_type,  &! INTENT(IN)
       level          ! INTENT(IN)

  use mem_grid, only: &
       hw4, &
       itopo, &
       ngrids,     & ! INTENT(IN)
       ngrid,      & ! INTENT(IN)
       npatch,     & ! INTENT(IN)
       time,       & ! INTENT(IN)
       dts,        &
       dtlong,     & ! INTENT(IN)
       dtlongn,    & ! INTENT(IN)
       ideltat,    &
       nnacoust,   &
       iyear1,     & ! INTENT(IN)
       imonth1,    & ! INTENT(IN)
       idate1,     & ! INTENT(IN)
       grid_g,     & ! INTENT(INOUT)
       nxtnest,    & ! INTENT(IN)
       if_adap,    & ! INTENT(IN)
       dtlt,       & ! INTENT(IN)
       istp,       & ! INTENT(IN)
       jdim,       & ! INTENT(IN)
       nzp,        & ! INTENT(IN)
       f_thermo_e, & ! INTENT(IN)
       f_thermo_w, & ! INTENT(IN)
       f_thermo_s, & ! INTENT(IN)
       f_thermo_n, &   ! INTENT(IN)
       zt,         &
       zm,         &
       dzt,        &
       itime1,     &
       vveldamp,   & ! INTENT(IN)
       ibnd,       &
       jbnd,       &
       nstbot,     &
       nnzp

  use shcu_vars_const, only: & ! For Shallow Cumulus Paramet.
       NNSHCU ! INTENT(IN)

  use mem_scalar, only: & ! For SiB
       scalar_g ! INTENT(IN)

  use mem_leaf, only: & 
       ISFCL          & ! INTENT(IN)
       ,ISFCL_OCEAN      ! INTENT(IN)

  ! TEB_SPM
  use teb_spm_start, only: &
       TEB_SPM ! INTENT(IN)
  use mem_emiss, only: &
       ichemi,         & ! INTENT(IN)
       isource           ! INTENT(IN)

  ! For specific optimization depending the type of machine
  use machine_arq, only: &
       machine ! INTENT(IN)

  use rconstants, only:  &
       g,                & ! (IN)
       cp,               & ! (IN)
       cpor,             & ! (IN)
       p00,              & ! (IN)
       rgas,             & ! (IN)
       pi180               ! (IN)

  use ccatt_start, only: &
       ccatt               ! (IN)

  use mem_tend, only: &
       tend

  use utilsMod, only: &
       Copy1DTo3D

  use mem_chem1, only: &
       nvert_src=>chem1_src_z_dim_g, & ! (IN)
       chem1_g,                      & ! (INOUT)
       nsrc,                         & ! (IN)
       chem1_src_g,                  & ! %sc_src(INOUT)
       chemistry,                    & ! (IN)
       split_method,                 & ! (IN)
       n_dyn_chem,                   &
       ntimes_src

  use mem_aer1, only:                  &
       aerosol,                        &! (IN)
       aer1_g,                         &! %sc_src(INOUT)
       aer2_g,                         &! %sc_src(INOUT)
       aer_nvert_src=>aer1_src_z_dim_g  ! (IN)

  use mem_plume_chem1, only:  plume_mean_g &   ! %flam_frac(IN), %fire_size(IN)
       ,plume_fre_g
  use mem_stilt, only: &
       iexev,          &  ! (IN)
       imassflx,       &  ! (IN)
       stilt_g            ! %dnp (IN)

  use mem_radiate, only: radiate_g

  use chem_sources, only :     &
       alloc_emiss_cycle,      &  ! Subroutine
       init_actual_time_index, &  ! Subroutine
       emiss_cycle,            &  ! (INOUT)
       emiss_cycle_alloc,      &
       srcmapfn                   ! (IN)

  use ChemSourcesDriver, only:  sources_driver            ! Subroutine

  use ChemDryDepDriver , only:  drydep_driver             ! Subroutine

  use ModChemistryDriver, only: &
       chemistry_driver, &
       aer_background

  use radiation, only: radiate ! Subroutine

  use ModTimeStamp, only: SynchronizedTimeStamp, TimeStamp

  use digitalfilter, only:         &
       applyDigitalFilter, & ! subroutine
       fileNameDF,& ! intent(inout) - file control
       dfVars,             &
       applyDF

  use ModMonotonicAdvection, only:                 &
       advmnt_driver,  &        ! subroutine
       advmnt

  use ModMatrixDriver, only: &
       MatrixDriver  !Matrix Aerosol Model


  use ModRamsMicrophysics2M, only: &
       micro_2M_rams60, &
       negadj1_2M_rams60

  use mem_radiate, only: &
       ilwrtyp, iswrtyp

  use MODCUPARGRELL3, only: g3d_g

  use ModWindFarm, only: &
       wind_farm_driver

  use ModOptical, only: &
       aodDriver

  use ModAerClim, only: &
       no_months, &
       no_src_types, &
       specieName, &
       aercam

  use modIau, only:    &
       CreateIauTendency &
       ,readIauTendency   &
       ,GetIauTendency    &
       ,timeWindowIAU     &
       ,applyIAU

  use ModLeaf3OceanOnly, only: &
       sfclyr_ocean_only

  use ModAdvectc_rk, only: &
       advectc_rk

  use mem_scratch, only: &
       scratch    

  use iso_fortran_env, only: &
       int64

  implicit none

  private
  public :: timestep_rk

contains



  subroutine timestep_rk(oneGrid)
    implicit none

    type(Grid), pointer :: oneGrid

    ! execution time instrumentation
    !    include "constants.h"
    include "tsNames.h"

    logical, parameter :: flag_Coriolis_in_every_RK_step = .false.

    integer :: l_rk
    real    :: rk_beta(3)
    integer :: rk_nmbr_small_timesteps(3)
    integer, parameter ::  &!rk_order = 2 ! for testing purposes (but works only with upwind1 and upwind3 advection)
         rk_order = 3 ! for Wicker, Skamarock (2002) MWR-scheme

    integer :: n
    character(len=2) :: crk
    logical :: singleProcRun
    character(len=2) :: cnzp
    integer :: nv,itime,i,j
    integer :: lenCopy
    character(len=256) :: julesFile

    logical,parameter :: stepDebug=.true.
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(timestep_rk)**"
    character(len=8) :: str(10)

    !MB: only for testing
    !integer :: nmbr_gpts
    !real    :: pm,tm

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    singleProcRun = nmachs == 1
    julesFile=oneGrid%Ramsin%julesin

!!$    call SynchronizedTimeStamp(TS_RESTO) ! Exper1.2, 2021_12

    !        +------------------------------------------------------------------+
    !        |   Timestep driver for the Runge-Kutta non-hydrostatic time-split |
    !        |      model.                                                      |
    !        +------------------------------------------------------------------+
    !
    !  Zero out all tendency arrays.
    !--------------------------------
    call tend0(oneGrid%ScalarTab, oneGrid%ScalarTabSize)  

    ! Implements the Incremental Analysis Update procedure -
    ! phase 2: add the IAU tendencies
    !-------------------------------
    if(applyIAU == 2 ) then
       call readIauTendency(ngrid, mzp*mxp*myp, mzp, mxp, myp,ia,iz,ja,jz ,time,dtlt)
       call GetIauTendency (ngrid, mzp*mxp*myp, mzp, mxp, myp,ia,iz,ja,jz ,time,dtlt)      
    endif

    !  Thermodynamic diagnosis
    !--------------------------------
    if (mcphys_type <= 1 .and. level/=3) then
       call thermo(mzp, mxp, myp, ia, iz, ja, jz, oneGrid%Basic, oneGrid%AveBasic)
    endif

    ! evolution of the Exner pressure: compression term
    if (iexev == 2) then
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'ADV', &
            oneGrid%Basic, oneGrid%AveBasic)
    end if

!!$    call SynchronizedTimeStamp(TS_DYNAMICS) ! Exper1.2, 2021_12

    if (CCATT==1 .and. chemistry >= 0) then
       call aodDriver(mzp,mxp,myp,ia,iz,ja,jz,ngrids,oneGrid%Basic)
    end if

    !  Radiation parameterization
    !--------------------------------
    call radiate(mzp,mxp,myp,ia,iz,ja,jz,mynum, oneGrid%Basic)

    !  Surface layer, soil and veggie model
    !----------------------------------------
    if (isfcl<=2) then
       call sfclyr(mzp,mxp,myp,ia,iz,ja,jz,ibcon, oneGrid%Basic)
#ifdef JULES
    elseif (isfcl == 5) then
       if (time==0.) then
          call sfclyr(mzp,mxp,myp,ia,iz,ja,jz,ibcon, oneGrid%Basic)
       end if
       call DeepCopyToTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
       call sfclyr_jules(mzp,mxp,myp,ia,iz,ja,jz,jdim,julesFile, &
            oneGrid%Basic, oneGrid%Turb)
       call DeepCopyFromTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
       !--- this combines the JULES land + LEAF ocean models.
       if (isfcl_ocean == 1) then
          call DeepCopyToTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
          call sfclyr_ocean_only  (mzp,mxp,myp,ia,iz,ja,jz,ibcon, &
               oneGrid%Basic, oneGrid%Turb)
          call DeepCopyFromTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
       end if
#endif
    endif

    !- Sea salt Aerossol inline source
    call SeaSaltDriver(ia,iz,ja,jz,ngrid,mxp,myp, oneGrid%Basic)

    !-- emission/deposition for CCATT chemistry models
    if (CCATT==1 .and. chemistry >= 0) then

       !- emission
       !- allocation for diurnal cycle of emission arrays   and
       !- actual_time_index on nodes
       if( (.not. emiss_cycle_alloc) .and. (chemistry >= 0) .and. &
            trim(srcmapfn) .ne.  'NONE' .and. trim(srcmapfn) .ne.  'none') then

          call alloc_emiss_cycle(mxp,myp,ngrids,nsrc)

          call init_actual_time_index(nsrc,ntimes_src)

       endif
       !

       !plume_mean_g(:,:) instead of plume_mean_g(:,ngrid) to avoid memory errors.
       !emiss_cycle(:,:)  instead of emiss_cycle(:,ngrid)  to avoid memory errors.
       !the same for the others var
       call sources_driver(ngrid, mzp,mxp,myp,ia,iz,ja,jz,                          &
            g,cp,cpor,p00,rgas,pi180,                                &
            radiate_g(ngrid)%cosz,oneGrid%Basic%theta,              &
            oneGrid%Basic%pp,oneGrid%Basic%pi0,oneGrid%Basic%rv,  &
            oneGrid%Basic%dn0,oneGrid%Basic%up,oneGrid%Basic%vp,  &
            time,iyear1,imonth1,idate1,itime1,dtlt,                  &
            grid_g(ngrid)%rtgt,grid_g(ngrid)%lpw,grid_g(ngrid)%glat, &
            grid_g(ngrid)%glon,zt,zm,dzt,nzpmax,                     &
            nvert_src    (:,:),                                  &
            chem1_g      (:,:),                                  &
            chem1_src_g  (:,:,:,:),                              &
            aer1_g       (:,:,:),                                &
            aer_nvert_src(:,:),                                  &
            plume_mean_g (:,:),                                  &
            stilt_g(ngrid)%dnp,                                  &
            emiss_cycle  (:,:),                                  &
            aer2_g       (:,:),                                  &
            plume_fre_g  (:,:)                                   )


       !- call dry deposition and sedimentation routines
       call drydep_driver(mzp,mxp,myp,ia,iz,ja,jz, oneGrid%Basic)

    endif

!!$    call SynchronizedTimeStamp(TS_PHYSICS) ! Exper1.2, 2021_12

    !  Send boundaries to adjoining nodes
    !-------------------------------------------
    call PostSendRecvMsgs(oneGrid%SelectedGhostZoneSend, oneGrid%SelectedGhostZoneRecv)

    !  Coriolis terms
    !  ----------------------------------------
    if ( .not. flag_Coriolis_in_every_RK_step ) then
       call corlos(mzp, mxp, myp, i0, j0, ia, iz, ja, jz, izu, jzv, &
            tend%ut, tend%vt, oneGrid%Basic, oneGrid%AveBasic)
    end if

!!$    call SynchronizedTimeStamp(TS_DYNAMICS) ! Exper1.2, 2021_12

    !  Cumulus parameterization version 1
    !----------------------------------------
    if (NNQPARM(ngrid)==1 .or. IF_CUINV==1) then
       call cuparm(oneGrid%Basic)
    end if

    !  Urban canopy parameterization
    !----------------------------------------
    if (OneGrid%Ramsin%if_urban_canopy==1) then
       call DeepCopyToTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
       call urban_canopy(oneGrid%Basic, oneGrid%Turb)
       call DeepCopyFromTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
    end if

    !  Analysis nudging and boundary condition
    !------------------------------------------
    if (NUD_TYPE>0) then
       call datassim(oneGrid%Basic)
    end if

    !  Observation data assimilation
    !----------------------------------------
    if (IF_ODA==1) then
       call oda_nudge(oneGrid%Basic)
    end if

    !  Nested grid boundaries
    !----------------------------------------
!!$    ! **(JP)** Comment out to avoid compilation errors 
!!$    ! at procedure nstbdriv, which is also commented out
!!$    if (nxtnest(ngrid)>=1) call nstbdriv()

    !  Rayleigh friction for theta
    !----------------------------------------
    call rayft(mxp,myp,mzp,mynum,ngrid,nnzp,if_adap,level,nodemyp,nodemxp,&
         scratch%vt3da,oneGrid%Basic)

    !  Get the overlap region between parallel nodes
    !---------------------------------------------------

!!$    call SynchronizedTimeStamp(TS_PHYSICS) ! Exper1.2, 2021_12

    call WaitSendRecvMsgs(oneGrid%SelectedGhostZoneSend, oneGrid%SelectedGhostZoneRecv)

    if (iexev == 2) then
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THA', &
            oneGrid%Basic, oneGrid%AveBasic)
    end if

    !- task 2:  NO production by "eclair"

!!$    call SynchronizedTimeStamp(TS_DYNAMICS) ! Exper1.2, 2021_12

    if (ccatt == 1) then
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,2,50,&
            oneGrid%Basic)
    end if

    !- CATT & Chemistry == CCATT
    !----------------------------------------
    if (ccatt==1 .and. split_method== 'PARALLEL' .and. n_dyn_chem==1) then
       ! task 3 : production/loss by chemical processes and inclusion of the
       ! chemistry tendency at the total tendency
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50,&
            oneGrid%Basic)
    endif
    if (ccatt==1 ) then
       ! task 4 : mass transfer between gas and liquid
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,4,50,&
            oneGrid%Basic)
    endif

    !---------------------------------------------------
    ! Shallow  cumulus parameterization by Souza
    if (NNSHCU(ngrid)==1) then
       call DeepCopyToTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
       call shcupa(oneGrid%Basic, oneGrid%Turb)
       call DeepCopyFromTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
    end if
    !---------------------------------------------------

    if (TEB_SPM==1) then
       ! Update urban emissions
       !----------------------------------------
       if (isource==1) then
          call sources_teb(mzp, mxp, myp, ia, iz, ja, jz, ngrid, ngrids)
       endif
       !  Update chemistry
       !----------------------------------------
       if (ichemi==1) then
          call ozone(mzp, mxp, myp, ia, iz, ja, jz, ngrid, dtlt, oneGrid%Basic)
       endif
    endif

!!$    call SynchronizedTimeStamp(TS_PHYSICS) ! Exper1.2, 2021_12

    if (iexev == 2) then
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THS', &
            oneGrid%Basic, oneGrid%AveBasic)
    end if

    !  Sub-grid diffusion terms
    !----------------------------------------
    if ((if_adap==0) .and. (OneGrid%Ramsin%ihorgrad==2)) then
       call diffuse_brams31(oneGrid%ScalarTab, oneGrid%ScalarTabSize, &
            oneGrid%Basic, oneGrid%Ramsin, oneGrid%Turb, oneGrid%AveTurb, oneGrid%Id)
    else
       call diffuse(oneGrid%ScalarTab, oneGrid%ScalarTabSize, oneGrid%Basic, &
            oneGrid%Turb, oneGrid%AveTurb, oneGrid%Ramsin, oneGrid%Id)
    endif

!!$    call SynchronizedTimeStamp(TS_DYNAMICS) ! Exper1.2, 2021_12

    !- STILT-BRAMS coupling (ML)
    if (imassflx == 1) then
       call prep_advflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz,ngrid, oneGrid%Basic)
    end if

    !- large and subgrid scale forcing for shallow and deep cumulus
    if( NNQPARM(ngrid) >=2  ) then
       call prepare_lsf(NNQPARM(ngrid), NNSHCU(ngrid),1, oneGrid%Basic, oneGrid%AveBasic)
    end if

    !- cumulus parameterizations options: G3d - GD-FIM and GF
    if (NNQPARM(ngrid)>=3) then
       call DeepCopyToTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
       call CUPARM_GRELL3_CATT(oneGrid,1,NNQPARM(ngrid),NNSHCU(ngrid))
       call DeepCopyFromTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
    end if

    !------------------------------------------------------------------------------
    ! init preparations for Runge-Kutta  -loop
    !------------------------------------------------------------------------------

!!$    call SynchronizedTimeStamp(TS_PHYSICS) ! Exper1.2, 2021_12

    if ( rk_order == 2 ) then
       ! Wicker, Skamarock (1998)-RK-scheme
       rk_beta(1) = 1.0 / 2.0    ! = beta(2,1) of Butcher tableau
       rk_beta(2) = 1.0          ! = beta(3,2) of Butcher tableau

       if ( mod( nnacoust(ngrid), 2) /= 0 ) then
          call fatal_error("ERROR in timestep_rk: nnacoust(ngrid) must be an integer multiple of 2")
       end if

       rk_nmbr_small_timesteps(1) = nnacoust(ngrid) / 2
       rk_nmbr_small_timesteps(2) = nnacoust(ngrid)

    else if ( rk_order == 3 ) then
       ! Wicker, Skamarock (2002)-RK-scheme
       rk_beta(1) = 1.0 / 3.0    ! = beta(2,1) of Butcher tableau
       rk_beta(2) = 1.0 / 2.0    ! = beta(3,2) of Butcher tableau
       rk_beta(3) = 1.0          ! = beta(4,3) of Butcher tableau

       if ( mod( nnacoust(ngrid), 6) /= 0 .and.  mod( nnacoust(ngrid), 4) /= 0 ) then
          call fatal_error("ERROR in timestep_rk: nnacoust(ngrid) must be an integer multiple of 6 or 4")
       end if
       rk_nmbr_small_timesteps(1) = nnacoust(ngrid) / 3
       rk_nmbr_small_timesteps(2) = nnacoust(ngrid) / 2
       rk_nmbr_small_timesteps(3) = nnacoust(ngrid)
       !------
       !- optimization : rk_nmbr_small_timesteps(1)=1 e dts=2* dtlt / nnacoust(ngrid) - uncoment the line below
       rk_nmbr_small_timesteps(1) = 1.
       !------

    else
       call fatal_error("ERROR in timestep_rk: false value for rk_order")
    end if

    dts = dtlt / nnacoust(ngrid)

    if ( apply_div_damping ) then
       ! calculation of divergence damping coeff:
       !MB: ATTENTION: dts different for different nests?!!
       if ( ideltat == 0 ) then
          ! it is sufficient to calculate alpha_div only once:
          if ( istp == 1 ) then
             if (dumpLocal) then
                call MsgDump(h//" invokes init_div_damping_coeff to calculate alpha_div on ideltat=0")
             end if
             call init_div_damping_coeff (oneGrid, dts)
          end if
       else
          if (dumpLocal) then
             call MsgDump(h//" invokes init_div_damping_coeff to calculate alpha_div on ideltat!=0")
          end if
          call init_div_damping_coeff (oneGrid, dts)
       end if
    end if


    ! begin of Runge-Kutta loop
    !---------------------------------------------------

    !  Lateral velocity boundaries - radiative
    !-------------------------------------------
    call latbnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon,nxtnest,ngrid,ibnd,jbnd, &
         grid_g(ngrid)%lpu,grid_g(ngrid)%lpv,oneGrid%Basic%up,&
         oneGrid%Basic%uc,tend%ut,oneGrid%Basic%vp,oneGrid%Basic%vc,&
         tend%vt,grid_g(ngrid)%dxt,grid_g(ngrid)%dyt)

!!$    call SynchronizedTimeStamp(TS_RK_RESTO) ! Exper1.2, 2021_12


    do l_rk = 1, rk_order

       !initialize the tendencies with the physics tendencies
       tend%ut_rk (:) =tend%ut (:)
       tend%vt_rk (:) =tend%vt (:)
       tend%wt_rk (:) =tend%wt (:)
       tend%pt_rk (:) =tend%pt (:)
       tend%tht_rk(:) =tend%tht(:)

!!$       call SynchronizedTimeStamp(TS_RK_RESTO) ! Exper1.2, 2021_12

       ! advection should give back tendencies
       ! ut_rk, vt_rk, wt_rk, pt_rk, tht_rk = physics tend + advection tendency

       !  Velocity advection
       !----------------------------------------
       if (dumpLocal) then
          call MsgDump(h//" invokes advectc_rk for V")
       end if
       call advectc_rk(oneGrid,'V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)

       !  advection of pi and theta_il
       if (dumpLocal) then
          call MsgDump(h//" invokes advectc_rk for THETAIL")
       end if
       call advectc_rk(oneGrid,'THETAIL',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)
       if (dumpLocal) then
          call MsgDump(h//" invokes advectc_rk for PI")
       end if
       call advectc_rk(oneGrid,'PI'     ,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)

!!$       call SynchronizedTimeStamp(TS_RK_ADV) ! Exper1.2, 2021_12

       if ( flag_Coriolis_in_every_RK_step ) then
          call corlos(mzp, mxp, myp, i0, j0, ia, iz, ja, jz, izu, jzv, &
               tend%ut_rk, tend%vt_rk, oneGrid%Basic, oneGrid%AveBasic)
       end if

       !  Buoyancy term for w equation
       !----------------------------------------
       call buoyancy(tend%wt_rk, oneGrid%Basic)

       if (dumpLocal) then
          call MsgDump(h//" starts exchanging borders of tend%tht_rk")
       end if
       call PostSendRecvMsgs(oneGrid%AcoustNewThtSend, oneGrid%AcoustNewThtRecv)

       if ( l_rk > 1 ) then
          ! (not necessary in the first RK substep)
          oneGrid%Basic%uc (:,:,:) = oneGrid%Basic%up (:,:,:)
          oneGrid%Basic%vc (:,:,:) = oneGrid%Basic%vp (:,:,:)
          oneGrid%Basic%wc (:,:,:) = oneGrid%Basic%wp (:,:,:)
          oneGrid%Basic%pc (:,:,:) = oneGrid%Basic%pp (:,:,:)
          oneGrid%Basic%thc(:,:,:) = oneGrid%Basic%thp(:,:,:)
       end if

       !-  Acoustic small timesteps
       call acoustic_new(oneGrid, rk_nmbr_small_timesteps(l_rk),l_rk )

       !- update thp -> thc (theta_il is not contained in acoustic_new)
       if (dumpLocal) then
          call MsgDump(h//" waits exchanging borders of tend%tht_rk")
       end if
       call WaitSendRecvMsgs(oneGrid%AcoustNewThtSend, oneGrid%AcoustNewThtRecv)

       call update_long_rk(int(mxp*myp*mzp,int64),dtlt,rk_beta(l_rk) &
            ,oneGrid%Basic%thc,oneGrid%Basic%thp  &
            ,tend%tht_rk)

       !- determine theta (dry potential temp.) for the buoyancy term:
       call theta_thp_rk(mzp,mxp,myp,ia,iz,ja,jz,"get_theta", &
            oneGrid%Basic, oneGrid%AveBasic)

       !-damping on vertical velocity to keep stability
       !MB: does this act on wc???
       if (vveldamp == 1) then
          call w_damping(mzp,mxp,myp,ia,iz,ja,jz,mynum, &
               oneGrid%Basic)
       end if
!!$       call SynchronizedTimeStamp(TS_RK_RESTO) ! Exper1.2, 2021_12

    end do
    ! end of Runge-Kutta loop
    ! i.e. the fields
    !    basic_g%uc, ..%vc, ..%wc, ..%pc, ..%thc
    ! contain the fields at the new time level n+1
    !---------------------------------------------------
    !
    if ( apply_div_damping ) then
       if ( ideltat /= 0 ) then
          call deallocate_alpha_div
       end if
    end if

!!$    call SynchronizedTimeStamp(TS_RK_RESTO) ! Exper1.2, 2021_12

    !
    !
    !  water species, tke and tracers advection
    !----------------------------------------
    if(advmnt == 1) then
       !- monotonic advection scheme
       call advmnt_driver(oneGrid, 'SCALAR',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,&
            i0,j0,nodemxp,nodemyp,nodemzp,mynum)
    elseif(advmnt == 0) then
       !- using the 2nd order forward upstream
       call advectc(oneGrid%ScalarTab, oneGrid%ScalarTabSize, oneGrid%Basic, oneGrid%AveBasic, &
            'SCALAR',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
    elseif(advmnt == 3) then
       !- using the WS advection
       call advectc_rk(oneGrid,'SCALAR',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)

    endif

!!$    call SynchronizedTimeStamp(TS_RK_ADVMON) ! Exper1.2, 2021_12

    !  Update scalars (water species, tke and tracers)
    !----------------------------------------
    call predtr(oneGrid%ScalarTab, oneGrid%ScalarTabSize)

    !-  copy current time into past time (u,v,w,pi,thetail)
    !---> thp   must be changed to THC for microphysics/bc/theta update
    !---> pp    must be changed to PC  for microphysics
    !---> wp    must be changed to WC  for microphysics
    !---> up,vp must be changed to UC,VC for output
    oneGrid%Basic%up (:,:,:) = oneGrid%Basic%uc (:,:,:)
    oneGrid%Basic%vp (:,:,:) = oneGrid%Basic%vc (:,:,:)
    oneGrid%Basic%wp (:,:,:) = oneGrid%Basic%wc (:,:,:)
    oneGrid%Basic%pp (:,:,:) = oneGrid%Basic%pc (:,:,:)
    oneGrid%Basic%thp(:,:,:) = oneGrid%Basic%thc(:,:,:)
    !---->
    !---->
    !
    !  Moisture variables positive definite
    !----------------------------------------
    if     (mcphys_type == 0) then
       call negadj1(mzp,mxp,myp, oneGrid%Basic)

    elseif(mcphys_type == 1) then
       call negadj1_2M_rams60(mzp,mxp,myp, oneGrid%Basic)
    endif

!!$    call SynchronizedTimeStamp(TS_RK_RESTO) ! Exper1.2, 2021_12

    !  Microphysics (applied on THP, just updated)
    !----------------------------------------
    if (mcphys_type == 0 .and. level==3) then
!!$       if (machine==1 .and. TEB_SPM==0) then
!!$          !- optimized version only for SX-6
!!$          call micro_opt()
!!$       else
       !- original Version used in a Generic IA32 machine
       call micro(oneGrid%Basic)
!!$       endif
    endif

    if (mcphys_type == 1 .and. level==3) then
       !- 2M rams microphysics
       call micro_2M_rams60(oneGrid%Basic)

    elseif (mcphys_type == 2 .or. mcphys_type == 3 ) then
       !- G. Thompson microphysics
       call micro_thompson(oneGrid%Basic)

    elseif(mcphys_type == 4 ) then
       call micro_gfdl(oneGrid%Basic)
    endif
    !----------------------------------------

!!$    call SynchronizedTimeStamp(TS_PHYSICS) ! Exper1.2, 2021_12

    !- Thermodynamic diagnosis
    if (mcphys_type <= 1 .and. level==3)  then
       call thermo(mzp, mxp, myp, 1, mxp, 1, myp, oneGrid%Basic, oneGrid%AveBasic)
    endif

    !  Apply scalar b.c.'s (THP is changed here)
    !----------------------------------------
    call trsets(oneGrid%ScalarTab, oneGrid%ScalarTabSize, oneGrid%Basic, &
         oneGrid%Turb, oneGrid%AveTurb)

    !---> THC must be changed to THP to include microphysics/trsets changes
    !---> for the next timestep
    oneGrid%Basic%thc(:,:,:) = oneGrid%Basic%thp(:,:,:)
    !--->

    !  Lateral velocity boundaries - radiative
    !-------------------------------------------
    !srf  call LATBND()

    !  Velocity/pressure boundary conditions
    !----------------------------------------
    call vpsets(mzp,mxp,myp,ia,iz,ja,jz,ibcon,nstbot, &
         oneGrid%Basic%up,oneGrid%Basic%vp,oneGrid%Basic%wp,&
         oneGrid%Basic%pp,oneGrid%Basic%uc,oneGrid%Basic%vc,&
         oneGrid%Basic%wc,oneGrid%Basic%pc,grid_g(ngrid)%dxu,&
         grid_g(ngrid)%dxm,grid_g(ngrid)%dyv,grid_g(ngrid)%dym,&
         grid_g(ngrid)%lpu,grid_g(ngrid)%lpv,grid_g(ngrid)%lpw, &
         oneGrid%Basic)

    !- call THERMO on the boundaries
    call thermo_boundary_driver((time+dtlongn(ngrid)), dtlong, &
         f_thermo_e(ngrid), f_thermo_w(ngrid), &
         f_thermo_s(ngrid), f_thermo_n(ngrid), &
         nzp, mxp, myp, jdim, oneGrid%Basic, oneGrid%AveBasic)

    if (iexev == 2) then
       call get_true_air_density(mzp,mxp,myp,ia,iz,ja,jz,oneGrid%Basic)
    end if

!!$    call SynchronizedTimeStamp(TS_DYNAMICS) ! Exper1.2, 2021_12

    !----------------------------------------
    !- chemistry - microphysics tranfers - sedimentation and tranfer from clouds to rain
    if (ccatt==1) then
       ! task 5 : sedimentation and mass transfer between clouds and rain
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,5,50,&
            oneGrid%Basic)
    endif

    !----------------------------------------
    !- chemistry/aerosol solvers
    if (ccatt==1) then
       if ( (split_method== 'PARALLEL' .and. N_DYN_CHEM > 1) .or.  &
            (split_method== 'SEQUENTIAL'                   ) .or.  &
            (split_method== 'SYMMETRIC'                    )       )then

          ! task 3 : production/loss by chemical processes and final updated
          !  of each specie
          call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50,&
               oneGrid%Basic)
       endif

       !- call Matrix Aerosol Model
       !- using symmetric/sequential spliting operator
       if(AEROSOL==2) then
          call MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp, oneGrid%Basic, oneGrid%Turb)
       endif
    endif
    if (ccatt==1 .and. aerosol == 1) then
       call aer_background(ngrid,mzp,mxp,myp,ia,iz,ja,jz)
    endif
    !----------------------------------------

    if (TEB_SPM==1) then
       !EDF  emission module
       if (isource==1) then
          ! Apply only for last finner grid
          if (ngrid==ngrids) then
             call le_fontes(ngrid, mzp, mxp, myp, &
                  npatch, ia, iz, ja, jz, (time+dtlongn(1)), oneGrid%Basic)
          endif
       endif
       !EDF
    endif

    !- windfarm
    call DeepCopyToTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)
    call wind_farm_driver(ngrid,mzp,mxp,myp,ia,iz,ja,jz, &
         oneGrid%Basic, oneGrid%Turb)
    call DeepCopyFromTurbFields(oneGrid%Turb, oneGrid%AveTurb,h)

    !- apply digital filter
    if (applyDF) then
       call applyDigitalFilter(fileNameDF, dfVars, oneGrid%Basic, oneGrid%AveBasic)
    end if


    ! Implements the Incremental Analysis Update procedure -
    !  phase 1: create and store the IAU tendencies
    !----------------------------------------
    !
    if(applyIAU == 1 .and. abs ( time + dtlt - timeWindowIAU*0.5) < 0.01) then

       if(mynum==1) print*,"timeIAU=",time,timeWindowIAU*0.5,&
            abs ( time - dtlt - timeWindowIAU*0.5),applyIAU
       if(mynum==1)call flush(6)

       call CreateIauTendency(ngrid, mzp*mxp*myp, mzp, mxp, myp,ia,iz,ja,jz&
            ,varinit_g(ngrid)%varup(:,:,:),varinit_g(ngrid)%varvp(:,:,:)  &
            ,varinit_g(ngrid)%varpp(:,:,:),varinit_g(ngrid)%vartp(:,:,:)  &
            ,varinit_g(ngrid)%varrp(:,:,:)                                &
            
            ,varinit_g(ngrid)%varuf(:,:,:),varinit_g(ngrid)%varvf(:,:,:)  &
            ,varinit_g(ngrid)%varpf(:,:,:),varinit_g(ngrid)%vartf(:,:,:)  &
            ,varinit_g(ngrid)%varrf(:,:,:)                                &
            
            ,oneGrid%Basic%up     (:,:,:)   ,oneGrid%Basic%vp  (:,:,:)  &
            ,oneGrid%Basic%theta  (:,:,:)   ,oneGrid%Basic%rtp (:,:,:)  &
            ,oneGrid%Basic%pp     (:,:,:)                                )
    endif

!!$    call SynchronizedTimeStamp(TS_PHYSICS) ! Exper1.2, 2021_12

  end subroutine timestep_rk
end module ModTimestepRK

