!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module ModTimestep

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
  
  use ModCoriolis, only: &
       corlos

  use ModRThrm, only: &
       thermo, &
       thermo_boundary_driver

  use ModRexev, only: &
       exevolve, &
       get_true_air_density

  use cuparm_grell3, only: &
       prepare_lsf

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
       BasicFields, &
       DeepCopyToBasicFields, &
       DeepCopyFromBasicFields

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


  use mem_cuparm, only: &
       NNQPARM, & ! INTENT(IN)
       IF_CUINV   ! INTENT(IN)

  use mem_varinit, only: &
       NUD_TYPE ! INTENT(IN)

  use mem_turb, only: &
       IF_URBAN_CANOPY, & ! INTENT(IN)
       ihorgrad           ! INTENT(IN)

  use mem_oda,   only: &
       if_oda ! INTENT(IN)

  use micphys,   only: &
       mcphys_type,  &! INTENT(IN)
       level          ! INTENT(IN)

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

  use mem_scalar, only: & ! For SiB
       scalar_g ! INTENT(IN)

  use mem_leaf, only: & ! For SiB
       ISFCL ! INTENT(IN)

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

  use mem_radiate, only: radiate_g

  use chem_sources, only :     &
       alloc_emiss_cycle,      &  ! Subroutine
       init_actual_time_index, &  ! Subroutine
       emiss_cycle,            &  ! (INOUT)
       emiss_cycle_alloc,      &
       srcmapfn                   ! (IN)

  use ChemSourcesDriver, only:  sources_driver            ! Subroutine

  use ChemDryDepDriver , only:  drydep_driver             ! Subroutine

  use module_chemistry_driver, only: chemistry_driver ! Subroutine

  use radiation, only: radiate ! Subroutine

  use ModTimeStamp, only: SynchronizedTimeStamp, TimeStamp

  use cuparm_grell3, only: cuparm_grell3_catt &  ! subroutine
       ,g3d_g

  use digitalFilter, only:         &
       applyDigitalFilter, & ! subroutine
       fileNameDF,& ! intent(inout) - file control
       dfVars,             &
       applyDF

  use ModMonotonicAdvection, only:                 &
       advmnt_driver,  &        ! subroutine
       advmnt

  use DriverMatrix, only: &
       MatrixDriver  !Matrix Aerosol Model


  use ModRamsMicrophysics2M, only: &
       micro_2M_rams60,&
       negadj1_2M_rams60

  use mem_radiate, only: &
       ilwrtyp, iswrtyp

  use CUPARM_GRELL3, only: g3d_g

  use ModWindFarm, only: &
       wind_farm_driver

  use ModOptical, only: &
       aodDriver

  use ModRbnd, only: latbnd, &
       vpsets, &
       rayft,  &
       trsets

  use mem_scratch, only: scratch

  implicit none

  private
  public :: timestep
  public :: w_damping

contains




  subroutine timestep(oneGrid)
    type(Grid), pointer, intent(in) :: oneGrid

    ! execution time instrumentation
    include "tsNames.h"

    integer, parameter :: acoshdp = 0
    character(len=256) :: julesFile

    character(len=*), parameter :: h="**(timestep)**"

    julesFile=oneGrid%Ramsin%julesin

    !        +-------------------------------------------------------------+
    !        |   Timestep driver for the hybrid non-hydrostatic time-split |
    !        |      model.                                                 |
    !        +-------------------------------------------------------------+

    !  Zero out all tendency arrays.
    !--------------------------------
    call tend0(oneGrid%ScalarTab, oneGrid%ScalarTabSize)

!!!!  IF( NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid) >=2 )CALL prepare_lsf_OLD(NNQPARM(ngrid), NNSHCU(ngrid),1)


    !  Thermodynamic diagnosis
    !--------------------------------
    if (mcphys_type <= 1 .and. level/=3) then
       call thermo(mzp, mxp, myp, ia, iz, ja, jz, oneGrid%Basic, oneGrid%AveBasic)
    endif

    if (iexev == 2) then
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'ADV', &
            oneGrid%Basic, oneGrid%AveBasic)
    end if

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_DYNAMICS)
    if (CCATT==1 .and. chemistry >= 1) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call aodDriver(mzp,mxp,myp,ia,iz,ja,jz,ngrids,oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    end if

    !  Radiation parameterization
    !--------------------------------
    call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
    call radiate(mzp,mxp,myp,ia,iz,ja,jz,mynum, oneGrid%Basic)
    call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)

    !  Surface layer, soil and veggie model
    !----------------------------------------
    if (isfcl<=2) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call sfclyr(mzp,mxp,myp,ia,iz,ja,jz,ibcon, oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
#ifdef JULES
    elseif (isfcl == 5) then
       if (time==0.) then
          call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
          call sfclyr(mzp,mxp,myp,ia,iz,ja,jz,ibcon, oneGrid%Basic)
          call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
       end if
       call SFCLYR_JULES(mzp,mxp,myp,ia,iz,ja,jz,jdim,julesFile)
       !DSM}
#endif
    endif

    !-LFR Sea salt Aerossol inline source
    call SeaSaltDriver(ia,iz,ja,jz,ngrid,mxp,myp)

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
       call drydep_driver(mzp,mxp,myp,ia,iz,ja,jz)

       !- call Matrix Aerosol Model
       !----------------------------------------
       if(AEROSOL==2) then
          call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
          call MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp, oneGrid%Basic)
          call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
       endif

    endif

!!!  !srf- large and subgrid scale forcing for shallow and deep cumulus
!!!  IF( NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid) >=2 )CALL prepare_lsf_OLD(NNQPARM(ngrid), NNSHCU(ngrid),2)

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_PHYSICS)

    !  Send boundaries to adjoining nodes
    !-------------------------------------------
    call PostSendRecvMsgs(oneGrid%SelectedGhostZoneSend, oneGrid%SelectedGhostZoneRecv)

    !  Coriolis terms
    !  ----------------------------------------
    call corlos(mzp, mxp, myp, i0, j0, ia, iz, ja, jz, izu, jzv,&
         tend%ut, tend%vt, oneGrid%Basic, oneGrid%AveBasic)

    !  Velocity advection
    !----------------------------------------
    call advectc(oneGrid%ScalarTab, oneGrid%ScalarTabSize, oneGrid%Basic, oneGrid%AveBasic, &
         'V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)


    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_DYNAMICS)

    !  Cumulus parameterization version 1
    !----------------------------------------
    if (NNQPARM(ngrid)==1 .or. IF_CUINV==1) then
       call cuparm()
    end if

    !  Urban canopy parameterization
    !----------------------------------------
    if (IF_URBAN_CANOPY==1) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call urban_canopy(oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    end if

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_PHYSICS)

    !  Analysis nudging and boundary condition
    !------------------------------------------

    if (NUD_TYPE>0) call DATASSIM()


    !  Observation data assimilation
    !----------------------------------------
    if (IF_ODA==1) call oda_nudge()

    !  Nested grid boundaries
    !----------------------------------------
!!$  ! **(JP)** Comment out to avoid compilation errors 
!!$  ! at procedure nstbdriv, which is also commented out
!!$  if (nxtnest(ngrid)>=1) call nstbdriv()

    !  Rayleigh friction for theta
    !----------------------------------------
    call rayft(mxp,myp,mzp,mynum,ngrid,nnzp,if_adap,level,nodemyp,nodemxp,&
         scratch%vt3da,oneGrid%Basic%theta,oneGrid%Basic%rv)

    !  Get the overlap region between parallel nodes
    !---------------------------------------------------
    call WaitSendRecvMsgs(oneGrid%SelectedGhostZoneSend, oneGrid%SelectedGhostZoneRecv)

    if (iexev == 2) then
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THA', &
            oneGrid%Basic, oneGrid%AveBasic)
    end if

    !  Sub-grid diffusion terms
    !----------------------------------------
    call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
    if ((if_adap==0) .and. (ihorgrad==2)) then
       call diffuse_brams31(oneGrid%ScalarTab, oneGrid%ScalarTabSize, oneGrid%Basic)
    else
       call diffuse(oneGrid%ScalarTab, oneGrid%ScalarTabSize, oneGrid%Basic)
    endif
    call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)

!!!!!  IF( NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid)>=2 ) CALL prepare_lsf_OLD(NNQPARM(ngrid), NNSHCU(ngrid),3)

    !  Velocity advection
    !----------------------------------------
    if(advmnt >= 1) then
       !-srf monotonic advection scheme
       call advmnt_driver(oneGrid, 'T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,&
            i0,j0,nodemxp,nodemyp,nodemzp,mynum)
       if(advmnt >= 2) &
            call advectc(oneGrid%ScalarTab, oneGrid%ScalarTabSize, oneGrid%Basic, oneGrid%AveBasic, &
            'T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
    else

       call advectc(oneGrid%ScalarTab, oneGrid%ScalarTabSize, oneGrid%Basic, oneGrid%AveBasic, &
            'T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)

    endif     ! If Generic IA32 use old Advction Scheme

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_DYNAMICS)

    !- STILT-BRAMS coupling (ML)
    if (imassflx == 1) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call prep_advflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz,ngrid, oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    end if

    !- large and subgrid scale forcing for shallow and deep cumulus
    !!1  IF(  NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid)>=2 ) CALL prepare_lsf_OLD(NNQPARM(ngrid), NNSHCU(ngrid),4)
    if( NNQPARM(ngrid) >=2 .or. NNSHCU(ngrid)>=2 ) then
       call prepare_lsf(NNQPARM(ngrid), NNSHCU(ngrid),1, oneGrid%Basic, oneGrid%AveBasic)
    end if

    !-   Cumulus parameterization options 2->6:
    !                    Deep Convection scheme
    !- call deep first, if there is deep convection , turn off shallow.
    if(NNQPARM(ngrid)==2) call CUPARM_GRELL_CATT(1)
    !
    !                    Shallow Convection scheme
    if(NNSHCU(ngrid)==2 ) call CUPARM_GRELL_CATT(2)
    !
    !- G3d - GD-FIM and GF
    if(NNQPARM(ngrid)>=3) call CUPARM_GRELL3_CATT(oneGrid,1,NNQPARM(ngrid),NNSHCU(ngrid))

    !- task 2:  NO production by "eclair"
    if (ccatt == 1) &
         call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,2,50)

    !- CATT & Chemistry == CCATT
    !----------------------------------------
    if (ccatt==1 .and. split_method== 'PARALLEL' .and. n_dyn_chem==1) then
       ! task 3 : production/loss by chemical processes and inclusion of the
       ! chemistry tendency at the total tendency
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50)
    endif
    if (ccatt==1 ) then
       ! task 4 : mass transfer between gas and liquid
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,4,50)
    endif

    !---------------------------------------------------
    ! Shallow  cumulus parameterization by Souza
    if (NNSHCU(ngrid)==1) call SHCUPA()
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
          call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
          call ozone(mzp, mxp, myp, ia, iz, ja, jz, ngrid, dtlt, oneGrid%Basic)
          call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
       endif
    endif

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_PHYSICS)

    !  Update scalars
    !----------------------------------------
    !
    call predtr(oneGrid%ScalarTab, oneGrid%ScalarTabSize)
    !
    !  Moisture variables positive definite
    !----------------------------------------
    if    (mcphys_type == 0) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call negadj1(mzp,mxp,myp, oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)


    elseif(mcphys_type == 1) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call negadj1_2M_rams60(mzp,mxp,myp, oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    endif

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_DYNAMICS)

    !  Microphysics
    !----------------------------------------
    if (mcphys_type == 0 .and. level==3) then
!!$       if (machine==1 .and. TEB_SPM==0) then
!!$          ! Optimized version only for SX-6
!!$          call micro_opt()
!!$       else
          ! Original Version used in a Generic IA32 machine
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call micro(oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
!!$       endif
    endif
    if (mcphys_type == 1 .and. level==3) then
       ! 2M rams microphysics
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call micro_2M_rams60(oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    endif
    if (mcphys_type == 2 .or. mcphys_type == 3 ) then
       ! G. Thompson microphysics
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call micro_thompson(oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    endif
    if (mcphys_type == 4 ) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call micro_gfdl(oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    endif

    !----------------------------------------
    !- chemistry - microphysics tranfers - sedimentation and tranfer from clouds to rain
    if (ccatt==1) then
       ! task 5 : sedimentation and mass transfer between clouds and rain
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,5,50)
    endif

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    !  call SynchronizedTimeStamp(TS_PHYSICS)

    !
    !  Thermodynamic diagnosis
    !----------------------------------------
    if (mcphys_type <= 1 .and. level==3)  then
       call thermo(mzp,mxp,myp,1,mxp,1,myp, oneGrid%Basic, oneGrid%AveBasic)
    endif

    if (iexev == 2) then
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THS', &
            oneGrid%Basic, oneGrid%AveBasic)
    end if

    !-damping on vertical velocity to keep stability
    if (vveldamp == 1) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call w_damping(mzp,mxp,myp,ia,iz,ja,jz,mynum, &
            oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    end if

    !  Apply scalar b.c.'s
    !----------------------------------------
    call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
    call trsets(oneGrid%ScalarTab, oneGrid%ScalarTabSize, oneGrid%Basic)
    call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)

    !  Lateral velocity boundaries - radiative
    !-------------------------------------------
    call latbnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon,nxtnest,ngrid,ibnd,jbnd, &
         grid_g(ngrid)%lpu,grid_g(ngrid)%lpv,oneGrid%Basic%up,&
         oneGrid%Basic%uc,tend%ut,oneGrid%Basic%vp,oneGrid%Basic%vc,&
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
    call hadvance(1, oneGrid%Basic, oneGrid%AveBasic)
    !  Buoyancy term for w equation
    !----------------------------------------
    call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
    call buoyancy(tend%wt, oneGrid%Basic)
    call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)

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
    call hadvance(2, oneGrid%Basic, oneGrid%AveBasic)

    !  Velocity/pressure boundary conditions
    !----------------------------------------
    call vpsets(mzp,mxp,myp,ia,iz,ja,jz,ibcon,nstbot, &
         oneGrid%Basic%up,oneGrid%Basic%vp,oneGrid%Basic%wp,&
         oneGrid%Basic%pp,oneGrid%Basic%uc,oneGrid%Basic%vc,&
         oneGrid%Basic%wc,oneGrid%Basic%pc,grid_g(ngrid)%dxu,&
         grid_g(ngrid)%dxm,grid_g(ngrid)%dyv,grid_g(ngrid)%dym,&
         grid_g(ngrid)%lpu,grid_g(ngrid)%lpv,grid_g(ngrid)%lpw)

    if (iexev == 2) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call get_true_air_density(mzp,mxp,myp,ia,iz,ja,jz,oneGrid%Basic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    end if

    ! Call thermo on the boundaries
    call thermo_boundary_driver((time+dtlongn(ngrid)), dtlong, &
         f_thermo_e(ngrid), f_thermo_w(ngrid), &
         f_thermo_s(ngrid), f_thermo_n(ngrid), &
         nzp, mxp, myp, jdim, oneGrid%Basic, oneGrid%AveBasic)

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
          call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50)
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

    if (TEB_SPM==1) then
       !EDF  emission module
       if (isource==1) then
          ! Apply only for last finner grid
          if (ngrid==ngrids) then
             call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
             call le_fontes(ngrid, mzp, mxp, myp, &
                  npatch, ia, iz, ja, jz, (time+dtlongn(1)), oneGrid%Basic)
             call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
          endif
       endif
       !EDF
    endif

    !windfarm
    call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
    call wind_farm_driver(ngrid,mzp,mxp,myp,ia,iz,ja,jz, &
         oneGrid%Basic)
    call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)

    !- apply digital filter
    if (applyDF) then
       call DeepCopyToBasicFields(oneGrid%Basic, oneGrid%AveBasic, h)
       call applyDigitalFilter(fileNameDF, dfVars, oneGrid%Basic, oneGrid%AveBasic)
       call DeepCopyFromBasicFields(oneGrid%Basic, oneGrid%AveBasic)
    end if

    !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
    ! call SynchronizedTimeStamp(TS_PHYSICS)

  end subroutine timestep

  !*************************************************************************

!!$  subroutine mass_flux(n1,n2,n3,m1,m2,m3,up,vp,wp  &
!!$       ,dn0,rtgu,rtgv,dyu,dxv,pp,pi0)
!!$
!!$    use mem_grid
!!$    use rconstants
!!$
!!$    implicit none
!!$    integer :: n1,n2,n3,m1,m2,m3
!!$    real :: up(m1,m2,m3),vp(m1,m2,m3),wp(m1,m2,m3)  &
!!$         ,dn0(n1,n2,n3),rtgu(n2,n3),dyu(n2,n3),dxv(n2,n3)  &
!!$         ,rtgv(n2,n3),pp(m1,m2,m3),pi0(n1,n2,n3)
!!$
!!$    real, save :: aintmass=0.
!!$
!!$    integer :: i,j,k
!!$    real :: wmass,emass,smass,nmass,prtot,tmass,ppp,area
!!$
!!$    !cc      if (mod(time,300.).gt..1) return
!!$
!!$    !  west/east bound
!!$    wmass=0.
!!$    emass=0.
!!$    do j=2,nyp-1
!!$       do k=2,nzp-1
!!$          i=1
!!$          wmass=wmass +  &
!!$               up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
!!$               *(dn0(k,i,j)+dn0(k,i+1,j))*.5
!!$          i=nxp-1
!!$          emass=emass -  &
!!$               up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
!!$               *(dn0(k,i,j)+dn0(k,i+1,j))*.5
!!$       enddo
!!$    enddo
!!$
!!$    !  north/south bound
!!$    smass=0.
!!$    nmass=0.
!!$    do i=2,nxp-1
!!$       do k=2,nzp-1
!!$          j=1
!!$          smass=smass +  &
!!$               vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
!!$               *(dn0(k,i,j)+dn0(k,i,j+1))*.5
!!$          j=nyp-1
!!$          nmass=nmass -  &
!!$               vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
!!$               *(dn0(k,i,j)+dn0(k,i,j+1))*.5
!!$       enddo
!!$    enddo
!!$
!!$    k=2
!!$    prtot=0.
!!$    do j=2,nyp-1
!!$       do i=2,nxp-1
!!$          ppp= ( (pp(k,i,j)+pi0(k,i,j))/cp )**cpor*p00
!!$          prtot=prtot+ppp/(dyu(i,j)*dxv(i,j))
!!$       enddo
!!$    enddo
!!$
!!$
!!$    tmass=wmass+emass+smass+nmass
!!$    aintmass=aintmass+tmass*dtlong
!!$    area=(nxp-2)*deltax*(nyp-2)*deltay
!!$
!!$
!!$    print*,'==============================='
!!$    print*,' Mass flux - W, E, S, N'
!!$    print*,  wmass,emass,smass,nmass
!!$    print*, 'total (kg/(m2 s):',tmass/area
!!$    print*, 'total (kg/m2):',aintmass/area
!!$    print*, 'total pr change (pa):',aintmass/area*9.8
!!$    print*, 'computed mean press:',prtot/area
!!$    print*,'==============================='
!!$
!!$    return
!!$  end subroutine mass_flux
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
