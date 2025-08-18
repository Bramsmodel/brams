!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module ModTimestep

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

  use ModLeaf3, only: &
       sfclyr

  use ModCoriolis, only: &
       corlos

  use ModRThrm, only: &
       thermo, &
       thermo_boundary_driver

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
       hadvance, &
       predtr

  use ModRadvc, only: &
       advectc

  use grid_dims, only: &
       nzpmax

  use ModMessageSet, only: &
       PostSendRecvMsgs, &
       WaitSendRecvMsgs

  use ModAcoust, only: &
       acoustic_new, &
       buoyancy

  use ModGrid, only: &
       Grid

  use ModBasicFields, only: &
       BasicFields

  use node_mod, only: &
       mzp, mxp, myp,  & ! INTENT(IN)
       ia, iz, ja, jz, & ! INTENT(IN)
       i0, j0,         & ! INTENT(IN)
       izu, jzv,       & ! INTENT(IN)
       mynum,          & ! INTENT(IN)
       ibcon,          & ! INTENT(IN)
       nmachs, &            ! INTENT(IN)
       nodemyp,  &  !intent(in)
       nodemxp,  &  !intent(in)
       nodemzp      !intent(in)

  use mem_varinit, only: &
       NUD_TYPE ! INTENT(IN)

  use mem_oda,   only: &
       if_oda ! INTENT(IN)

  use mem_grid, only: &
       ht,         & ! intent(in)
       ngrids,     & ! INTENT(IN)
       ngrid,      & ! INTENT(IN)
       npatch,     & ! INTENT(IN)
       time,       & ! INTENT(IN)
       dts,        &
       dtlong,     & ! INTENT(IN)
       dtlongn,    & ! INTENT(IN)
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
       vveldamp,   &
       ibnd,       &
       jbnd,       &
       nstbot,     &
       nnzp

  use shcu_vars_const, only: & ! For Shallow Cumulus Paramet.
       NNSHCU ! INTENT(IN)

  use mem_leaf, only: & ! For SiB
       ISFCL ! INTENT(IN)

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
       ,plume_fre_g      !

  use mem_stilt, only: &
       iexev,          &  ! (IN)
       imassflx,       &  ! (IN)
       stilt_g            ! %dnp (IN)

  use chem_sources, only :     &
       alloc_emiss_cycle,      &  ! Subroutine
       init_actual_time_index, &  ! Subroutine
       emiss_cycle,            &  ! (INOUT)
       emiss_cycle_alloc,      &
       srcmapfn                   ! (IN)

  use ChemSourcesDriver, only: &
       sources_driver            ! Subroutine

  use ChemDryDepDriver , only: &
       drydep_driver             ! Subroutine

  use ModChemistryDriver, only: &
       chemistry_driver, &
       aer_background

  use radiation, only: radiate ! Subroutine

  use ModTimeStamp, only: SynchronizedTimeStamp, TimeStamp


  use digitalFilter, only:         &
       applyDigitalFilter, & ! subroutine
       fileNameDF,& ! intent(inout) - file control
       dfVars,             &
       applyDF

  use ModMonotonicAdvection, only:                 &
       advmnt_driver,  &        ! subroutine
       advmnt

  use ModMatrixDriver, only: &
       MatrixDriver  !Matrix Aerosol Model

  use MODCUPARGRELL3, only: g3d_g

  use ModWindFarm, only: &
       wind_farm_driver

  use ModOptical, only: &
       aodDriver

  use ModRbnd, only: latbnd, &
       vpsets, &
       rayft,  &
       trsets

  implicit none

  private
  public :: timestep
  public :: w_damping

contains




  subroutine timestep(oneGrid)
    type(Grid), pointer :: oneGrid

    ! execution time instrumentation
    include "tsNames.h"

    integer, parameter :: acoshdp = 0
    character(len=256) :: julesFile

    character(len=*), parameter :: h="**(timestep)**"
    real, pointer :: vt3da(:)      
    julesFile=oneGrid%oneNamelistFile%julesin

    !        +-------------------------------------------------------------+
    !        |   Timestep driver for the hybrid non-hydrostatic time-split |
    !        |      model.                                                 |
    !        +-------------------------------------------------------------+

    !  Zero out all tendency arrays.
    !--------------------------------
    call tend0(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize)


    !  Thermodynamic diagnosis
    !--------------------------------
    if (oneGrid%oneMicVars%mcphys_type == 0 .and. oneGrid%oneMicVars%level/=3) then
       call thermo(mzp, mxp, myp, ia, iz, ja, jz, &
            oneGrid%oneBasicFields, oneGrid%oneMicVars, oneGrid%oneMicroFields)
    endif

    if (iexev == 2) then
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'ADV', &
            oneGrid%oneBasicFields)
    end if

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_DYNAMICS)
    if (CCATT==1 .and. chemistry >= 1) then
       call aodDriver(mzp,mxp,myp,ia,iz,ja,jz,ngrids,&
            oneGrid%oneNamelistFile, oneGrid%oneBasicFields, oneGrid%oneTurbFields, oneGrid%Id, oneGrid%oneControlVars)
    end if

    !  Radiation parameterization
    !--------------------------------
    call radiate(mzp,mxp,myp,ia,iz,ja,jz,mynum, &
         oneGrid%oneNamelistFile, oneGrid%oneBasicFields, &
         oneGrid%oneMicVars, oneGrid%oneMicroFields, &
         oneGrid%oneRadiateFields, oneGrid%oneCuParmFields)

    !  Surface layer, soil and veggie model
    !----------------------------------------
    if (isfcl<=2) then
       call sfclyr(mzp,mxp,myp,ia,iz,ja,jz,ibcon, &
            oneGrid%oneNamelistFile, oneGrid%oneBasicFields, oneGrid%oneTurbFields, &
            oneGrid%oneMicVars, oneGrid%oneMicroFields, oneGrid%oneRadiateFields, &
            oneGrid%oneCuParmFields)
#ifdef JULES
    elseif (isfcl == 5) then
       if (time==0.) then
          call sfclyr(mzp,mxp,myp,ia,iz,ja,jz,ibcon, &
               oneGrid%oneNamelistFile, oneGrid%oneBasicFields, oneGrid%oneTurbFields, &
               oneGrid%oneMicVars, oneGrid%oneMicroFields, oneGrid%oneRadiateFields, &
               oneGrid%oneCuParmFields)
       end if

       call sfclyr_jules(mzp,mxp,myp,ia,iz,ja,jz,jdim,julesFile, &
            oneGrid%oneNamelistFile, oneGrid%oneBasicFields, oneGrid%oneTurbFields, oneGrid%oneMicVars, &
            oneGrid%oneMicroFields, oneGrid%oneJulesFields, oneGrid%oneRadiateFields, oneGrid%oneCuParmFields)
#endif
    endif

    !-LFR Sea salt Aerossol inline source
    call SeaSaltDriver(ia,iz,ja,jz,ngrid,mxp,myp, oneGrid%oneBasicFields)

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
            oneGrid%oneRadiateFields%cosz,oneGrid%oneBasicFields%theta,              &
            oneGrid%oneBasicFields%pp,oneGrid%oneBasicFields%pi0,oneGrid%oneBasicFields%rv,  &
            oneGrid%oneBasicFields%dn0,oneGrid%oneBasicFields%up,oneGrid%oneBasicFields%vp,  &
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
       call drydep_driver(mzp,mxp,myp,ia,iz,ja,jz, &
            oneGrid%oneNamelistFile, &
            oneGrid%oneBasicFields, &
            oneGrid%oneTurbFields, &
            oneGrid%oneMicVars, &
            oneGrid%oneMicroFields, &
            oneGrid%oneRadiateFields, &
            oneGrid%oneCuParmFields)

       !- call Matrix Aerosol Model
       !----------------------------------------
       if(AEROSOL==2) then
          call MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp, oneGrid%oneBasicFields, oneGrid%oneTurbFields, &
               oneGrid%oneMicVars, oneGrid%oneMicroFields, oneGrid%oneAero2McPhysFields)
       endif

    endif

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_PHYSICS)

    !  Send boundaries to adjoining nodes
    !-------------------------------------------
    call PostSendRecvMsgs(oneGrid%SelectedGhostZoneSend, oneGrid%SelectedGhostZoneRecv)

    !  Coriolis terms
    !  ----------------------------------------
    call corlos(mzp, mxp, myp, i0, j0, ia, iz, ja, jz, izu, jzv,&
         tend%ut, tend%vt, oneGrid%oneBasicFields)

    !  Velocity advection
    !----------------------------------------
    call advectc(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, oneGrid%oneBasicFields, &
         'V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)


    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_DYNAMICS)

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_PHYSICS)

    !  Analysis nudging and boundary condition
    !------------------------------------------

    if (NUD_TYPE>0) then
       call datassim(oneGrid%oneBasicFields)
    end if


    !  Observation data assimilation
    !----------------------------------------
    if (IF_ODA==1) then
       call oda_nudge(oneGrid%oneBasicFields)
    end if

    !  Nested grid boundaries
    !----------------------------------------
!!$  ! **(JP)** Comment out to avoid compilation errors 
!!$  ! at procedure nstbdriv, which is also commented out
!!$  if (nxtnest(ngrid)>=1) call nstbdriv()

    !  Rayleigh friction for theta
    !----------------------------------------
    call rayft(mxp,myp,mzp,mynum,ngrid,nnzp,if_adap,oneGrid%oneMicVars%level,nodemyp,nodemxp,&
         vt3da,oneGrid%oneBasicFields)

    !  Get the overlap region between parallel nodes
    !---------------------------------------------------
    call WaitSendRecvMsgs(oneGrid%SelectedGhostZoneSend, oneGrid%SelectedGhostZoneRecv)

    if (iexev == 2) then
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THA', &
            oneGrid%oneBasicFields)
    end if

    !  Sub-grid diffusion terms
    !----------------------------------------
    if ((if_adap==0) .and. (OneGrid%oneNamelistFile%ihorgrad==2)) then
       call diffuse_brams31(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, &
            oneGrid%oneBasicFields, oneGrid%oneNamelistFile, oneGrid%oneTurbFields, oneGrid%Id, &
            oneGrid%oneMicVars, oneGrid%oneMicroFields)
    else
       call diffuse(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, oneGrid%oneBasicFields, &
            oneGrid%oneTurbFields, oneGrid%oneNamelistFile, oneGrid%Id, &
            oneGrid%oneMicVars, oneGrid%oneMicroFields)
    endif

    !  Velocity advection
    !----------------------------------------
    if(advmnt >= 1) then
       !-srf monotonic advection scheme
       call advmnt_driver(&
            oneNodeDimensions=oneGrid%oneNodeDimensions, &
            oneNodeDimensionsMonAdvX=oneGrid%oneNodeDimensionsMonAdvX, &
            oneNodeDimensionsMonAdvY=oneGrid%oneNodeDimensionsMonAdvY, &
            oneBasicFields=oneGrid%oneBasicFields, &
            oneScalarTableSize=oneGrid%oneScalarTableSize, &
            oneScalarTable=oneGrid%oneScalarTable, &
            Id=oneGrid%Id, &
            varn='T', &
            oneMicVars=oneGrid%oneMicVars, &
            MonAdvInputSend=oneGrid%MonAdvInputSend, &
            MonAdvInputRecv=oneGrid%MonAdvInputRecv, &
            ConvertBramsToMonAdvX=oneGrid%ConvertBramsToMonAdvX, &
            ConvertBramsToMonAdvY=oneGrid%ConvertBramsToMonAdvY, &
            ConvertMonAdvXToMonAdvY=oneGrid%ConvertMonAdvXToMonAdvY, &
            ConvertBramsToBrams=oneGrid%ConvertBramsToBrams, &
            ConvertMonAdvYToBrams=oneGrid%ConvertMonAdvYToBrams)
       if(advmnt >= 2) &
            call advectc(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, oneGrid%oneBasicFields, &
            'T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
    else

       call advectc(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, oneGrid%oneBasicFields, &
            'T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)

    endif     ! If Generic IA32 use old Advction Scheme

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_DYNAMICS)

    !- STILT-BRAMS coupling (ML)
    if (imassflx == 1) then
       call prep_advflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz,ngrid, oneGrid%oneBasicFields)
    end if

    !- large and subgrid scale forcing for shallow and deep cumulus
    if( oneGrid%oneNamelistFile%nnqparm(ngrid) >=3 .or. NNSHCU(ngrid)>=2 ) then
       call prepare_lsf(oneGrid%oneNamelistFile%nnqparm(ngrid), NNSHCU(ngrid),1, &
            oneGrid%oneNamelistFile, oneGrid%oneBasicFields, oneGrid%oneRadiateFields, &
            oneGrid%oneCuParmShFields)
    end if

    !-   Cumulus parameterization options 2->6:
    !
    !- G3d - GD-FIM and GF
    if (oneGrid%oneNamelistFile%nnqparm(ngrid)>=3) then
       call cuparm_grell3_catt(onegrid,1,&
            oneGrid%oneNamelistFile%nnqparm(ngrid),nnshcu(ngrid))
    end if

    !- task 2:  NO production by "eclair"
    if (ccatt == 1) then
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,2,50,&
            oneGrid%oneNamelistFile, &
            oneGrid%oneBasicFields, &
            oneGrid%oneMicroFields, &
            oneGrid%oneRadiateFields)
    end if

    !- CATT & Chemistry == CCATT
    !----------------------------------------
    if (ccatt==1 .and. split_method== 'PARALLEL' .and. n_dyn_chem==1) then
       ! task 3 : production/loss by chemical processes and inclusion of the
       ! chemistry tendency at the total tendency
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50,&
            oneGrid%oneNamelistFile, &
            oneGrid%oneBasicFields, &
            oneGrid%oneMicroFields, &
            oneGrid%oneRadiateFields)
    endif
    if (ccatt==1 ) then
       ! task 4 : mass transfer between gas and liquid
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,4,50,&
            oneGrid%oneNamelistFile, &
            oneGrid%oneBasicFields, &
            oneGrid%oneMicroFields, &
            oneGrid%oneRadiateFields)
    endif

    !---------------------------------------------------
    ! Shallow  cumulus parameterization by Souza
    if (NNSHCU(ngrid)==1) then
       call shcupa(oneGrid%oneBasicFields, oneGrid%oneTurbFields, &
            oneGrid%oneMicroFields, oneGrid%oneShcuFields)
    end if
    !---------------------------------------------------

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_PHYSICS)

    !  Update scalars
    !----------------------------------------
    !
    call predtr(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize)
    !
    !  Moisture variables positive definite
    !----------------------------------------
    if    (oneGrid%oneMicVars%mcphys_type == 0) then
       call negadj1(mzp,mxp,myp, oneGrid%oneBasicFields, oneGrid%oneMicVars, oneGrid%oneMicroFields)
    endif

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_DYNAMICS)

    !  Microphysics
    !----------------------------------------
    if (oneGrid%oneMicVars%mcphys_type == 0 .and. oneGrid%oneMicVars%level==3) then
       ! Original Version used in a Generic IA32 machine
       call micro(oneGrid%oneBasicFields, oneGrid%oneMicVars, oneGrid%oneMicroFields)
    else if (oneGrid%oneMicVars%mcphys_type == 2 .or. oneGrid%oneMicVars%mcphys_type == 3 ) then
       ! G. Thompson microphysics
       call micro_thompson(oneGrid%oneNamelistFile, oneGrid%oneBasicFields, &
            oneGrid%oneMicVars, oneGrid%oneMicroFields)
    else if (oneGrid%oneMicVars%mcphys_type == 4 ) then
       call micro_gfdl(oneGrid%oneNamelistFile, oneGrid%oneBasicFields, oneGrid%oneMicroFields)
    endif

    !----------------------------------------
    !- chemistry - microphysics tranfers - sedimentation and tranfer from clouds to rain
    if (ccatt==1) then
       ! task 5 : sedimentation and mass transfer between clouds and rain
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,5,50,&
            oneGrid%oneNamelistFile, &
            oneGrid%oneBasicFields, &
            oneGrid%oneMicroFields, &
            oneGrid%oneRadiateFields)
    endif

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_PHYSICS)

    !
    !  Thermodynamic diagnosis
    !----------------------------------------
    if (oneGrid%oneMicVars%mcphys_type == 0 .and. oneGrid%oneMicVars%level==3)  then
       call thermo(mzp,mxp,myp,1,mxp,1,myp, &
            oneGrid%oneBasicFields, oneGrid%oneMicVars, oneGrid%oneMicroFields)
    endif

    if (iexev == 2) then
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THS', &
            oneGrid%oneBasicFields)
    end if

    !-damping on vertical velocity to keep stability
    if (vveldamp == 1) then
       call w_damping(mzp,mxp,myp,ia,iz,ja,jz,mynum, &
            oneGrid%oneBasicFields)
    end if

    !  Apply scalar b.c.'s
    !----------------------------------------
    call trsets(oneGrid%oneScalarTable, oneGrid%oneScalarTableSize, oneGrid%oneBasicFields,&
         oneGrid%oneTurbFields,oneGrid%oneMicVars,oneGrid%oneMicroFields)

    !  Lateral velocity boundaries - radiative
    !-------------------------------------------
    call latbnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon,nxtnest,ngrid,ibnd,jbnd, &
         grid_g(ngrid)%lpu,grid_g(ngrid)%lpv,oneGrid%oneBasicFields%up,&
         oneGrid%oneBasicFields%uc,tend%ut,oneGrid%oneBasicFields%vp,oneGrid%oneBasicFields%vc,&
         tend%vt,grid_g(ngrid)%dxt,grid_g(ngrid)%dyt)

    !  Apply Asselin time filter
    !---------------------------------------------------
    !
    !  Af(t)=A(t)+gama*(Af(t-deltat)-2*A(t)+A(t+deltat))
    !
    !  Where A denotes quantities and,
    !  Af is the filtered quantities
    !---------------------------------------------------

    !  First stage Asselin filter
    !------------------------------------------
    !
    !  scratch=A(t)+gama*(Af(t-deltat)-2*A(t))
    !
    !       +--------------+--------------+
    !       | IN           | OUT          |
    !  +----+--------------+--------------|
    !  | AP | Af(t-deltat) | Af(t-deltat) |
    !  |----+--------------+--------------|
    !  | AC | A(t)         | scratch      |
    !  +----+--------------+--------------+
    !
    !------------------------------------------
    call hadvance(1, oneGrid%oneBasicFields)
    !  Buoyancy term for w equation
    !----------------------------------------
    call buoyancy(tend%wt, oneGrid%oneBasicFields, oneGrid%oneMicVars)

    !  Acoustic small timesteps
    !----------------------------------------
    ! AF: OUT: A(t+deltat)

    dts = 2. * dtlt / nnacoust(ngrid)

    call acoustic_new(oneGrid, nnacoust(ngrid),0 )   !MB:

    !  Last stage of Asselin filter
    !------------------------------------------
    !
    !  Af(t)=scratch+gama*A(t+deltat)
    !
    !       +--------------+--------------+
    !       | IN           | OUT          |
    !  +----+--------------+--------------|
    !  | AP ! A(t+deltat)  | Af(t)        |
    !  |----+--------------+--------------|
    !  | AC ! scratch      | A(t+deltat)  |
    !  +----+--------------+--------------+
    !
    !------------------------------------------
    call hadvance(2, oneGrid%oneBasicFields)

    !  Velocity/pressure boundary conditions
    !----------------------------------------
    call vpsets(mzp,mxp,myp,ia,iz,ja,jz,ibcon,nstbot, &
         !TO oneGrid%oneBasicFields%up,oneGrid%oneBasicFields%vp,oneGrid%oneBasicFields%wp,&
         !TO oneGrid%oneBasicFields%pp,oneGrid%oneBasicFields%uc,oneGrid%oneBasicFields%vc,&
         !TO oneGrid%oneBasicFields%wc,oneGrid%oneBasicFields%pc,grid_g(ngrid)%dxu,&
         grid_g(ngrid)%dxu,&
         grid_g(ngrid)%dxm,grid_g(ngrid)%dyv,grid_g(ngrid)%dym,&
         grid_g(ngrid)%lpu,grid_g(ngrid)%lpv,grid_g(ngrid)%lpw, &
         oneGrid%oneBasicFields)

    if (iexev == 2) then
       call get_true_air_density(mzp,mxp,myp,ia,iz,ja,jz,&
            oneGrid%oneBasicFields,oneGrid%oneMicVars)
    end if

    ! Call thermo on the boundaries
    call thermo_boundary_driver((time+dtlongn(ngrid)), dtlong, &
         f_thermo_e(ngrid), f_thermo_w(ngrid), &
         f_thermo_s(ngrid), f_thermo_n(ngrid), &
         nzp, mxp, myp, jdim, oneGrid%oneBasicFields, oneGrid%oneMicVars, oneGrid%oneMicroFields)

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_DYNAMICS)

    !
    !----------------------------------------
    !- chemistry/aerosol solvers
    if (ccatt==1) then
       if ( (split_method== 'PARALLEL' .and. N_DYN_CHEM > 1) .or.  &
            (split_method== 'SEQUENTIAL'                   ) .or.  &
            (split_method== 'SYMMETRIC'                    )       )then

          ! task 3 : production/loss by chemical processes and final updated
          !  of each specie
          call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50,&
            oneGrid%oneNamelistFile, &
            oneGrid%oneBasicFields, &
            oneGrid%oneMicroFields, &
            oneGrid%oneRadiateFields)
       endif

       !- call Matrix Aerosol Model
       !- using symmetric/sequential spliting operator
       !if(AEROSOL==2) then
       !
       !   CALL MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp)
       !endif
    endif
    if (ccatt==1 .and. aerosol == 1) then
       call aer_background(ngrid,mzp,mxp,myp,ia,iz,ja,jz)
    endif
    !----------------------------------------

    !windfarm
    call wind_farm_driver(ngrid,mzp,mxp,myp,ia,iz,ja,jz, &
         oneGrid%oneBasicFields, oneGrid%oneTurbFields)

    !- apply digital filter
    if (applyDF) then
       call applyDigitalFilter(fileNameDF, dfVars, &
            oneGrid%oneBasicFields, oneGrid%oneControlVars, &
            oneGrid%oneVarTable, oneGrid%oneVarTableSize)
    end if

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    ! call SynchronizedTimeStamp(TS_PHYSICS)

  end subroutine timestep

  !     *****************************************************************

  subroutine w_damping(mzp,mxp,myp,ia,iz,ja,jz,mynum,&
       oneBasicFields)
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    type(BasicFields), pointer, intent(in) :: oneBasicFields

    call apply_wdamp(mzp,mxp,myp,ia,iz,ja,jz,mynum   &
         ,oneBasicFields%up   ,oneBasicFields%vp     &
         ,oneBasicFields%wp   ,grid_g(ngrid)%rtgt    &
         ,grid_g(ngrid)%f13t  ,grid_g(ngrid)%f23t    &
         ,grid_g(ngrid)%dxt   ,grid_g(ngrid)%dyt     &
         ,tend%ut             ,tend%vt               &
         ,tend%wt                                            )

  end subroutine w_damping

  ! **********************************************************************

  subroutine apply_wdamp(m1,m2,m3,ia,iz,ja,jz,mynum,up,vp,wp,rtgt,f13t,f23t,dxt,dyt,ut,vt,wt)

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real,    intent(in) :: up(m1,m2,m3)
    real,    intent(in) :: vp(m1,m2,m3)
    real,    intent(in) :: wp(m1,m2,m3)
    real,    intent(in) :: rtgt(m2,m3)
    real,    intent(in) :: f13t(m2,m3)
    real,    intent(in) :: f23t(m2,m3)
    real,    intent(in) :: dxt(m2,m3)
    real,    intent(in) :: dyt(m2,m3)
    real,    intent(inout) :: ut(m1,m2,m3)
    real,    intent(inout) :: vt(m1,m2,m3)
    real,    intent(inout) :: wt(m1,m2,m3)

    integer :: i,j,k,ifm,icm,innest
    real :: c1x,c1y,c1z,cflnumh,cflnumv,cflz
    real :: vctr1(m1)
    real :: vctr2(m1)
    real :: vctr3(m1)
    real , parameter ::gama_w=0.3 !m/s²
    !     This routine damps the vertical velocity when CFLZ is exceed
    !     (Actually check on 80% of CFL)

    cflnumh = .80
    cflnumv = .50
    cflz=0.0
    !
    !c1x=0.0
    !vctr3(:)=0.0

    do j = ja,jz
       do i = ia,iz
          do k = 2,m1-1

             !cflx
             !vctr1(k) = .5*(up(k,i,j)+up(k,i-1,j))*dtlt*dxt(i,j)
             !cfly
             !vctr2(k) = .5*(vp(k,i,j)+vp(k,i,j-jdim))*dtlt*dyt(i,j)
             !cflz
             vctr3(k) = ((wp(k,i,j)+wp(k-1,i,j))  &
                  +(up(k,i,j)+up(k,i-1,j))*f13t(i,j)*ht(k)*rtgt(i,j)  &
                  +(vp(k,i,j)+vp(k,i,j-jdim))*f23t(i,j)*ht(k)*rtgt(i,j)  &
                  )*.5*dtlt*dzt(k)
             c1z=abs(vctr3(k))
             if(c1z > cflnumv) then
                wt(k,i,j) = wt(k,i,j) -gama_w*sign(1.,wp(k,i,j))*(c1z-cflnumv)
                print*,'wdamp applied at=',k,i,j,mynum,c1z,wp(k,i,j)
                call flush(6)
             endif
          enddo
          !do k = 2,m1-1
          !   c1x = abs(vctr1(k))
          !   c1y = abs(vctr2(k))
          !   c1z = abs(vctr3(k))
          !
          !   if (c1x .gt. cflxy(ngrid)) cflxy(ngrid) = c1x
          !   if (c1y .gt. cflxy(ngrid)) cflxy(ngrid) = c1y
          !    if (c1z .gt. cflz) cflz = c1z
          !enddo
       enddo
    enddo
    !print*,'at wdamp2-max cflz',mynum,cflz
    !call flush(6)

  end subroutine apply_wdamp

  ! ***************************************************************************
end module ModTimestep
