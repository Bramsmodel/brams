!###########################################################################
! CCATT- B - Regional Atmospheric Modeling System - RAMS
!###########################################################################
module chem_sources

  use ModControlVars, only: &
       ControlVars
  
  use mem_chem1, only:       &
       chem1_vars,           & ! Type
       chem1_src_vars          ! Type

  use ModChem1Constants, only: &
       maxsrcfiles
  
  use mem_aer1, only:        &
       aer1_vars

  use mem_plume_chem1, only: &
       plume_vars,           & ! Type
       plume_mean_vars,      & ! Type
       plume_fre_vars,       &
       iflam_frac ,          &
       imean_frp  ,     &
       istd_frp   ,     &
       imean_size ,     &
       istd_size  

  use mem_volc_chem1, only:  &
       volc_mean_vars          ! Type

  use ModNamelistFile, only: &
       namelistFile

  use ModDateUtils, only: &
       date_abs_secs2,    &    ! Subroutine
       date_add_to,       &    ! Subroutine
       date_make_big,     &    ! Subroutine
       julday                  ! Function

  use ReadBcst, only: &
       ReadStoreFullFieldAndOwnChunk, & ! Subroutine
       Broadcast                        ! Subroutine

  use ParLib, only: parf_bcast ! Subroutine

  use io_params, only:    &
       srctime1,          &
       srctime2,          &
       fnames_src,        &
       itotdate_src,      &
       src_times,         &
       actual_time_index, &
       nsrcfiles,         &
       next_srcfile


  implicit none

  include "constants.h"


  character(len=256)  :: srcmapfn

  character(len=32)   :: def_proc_src ! 'last_sources' 
  !    if does not exist => keep the current sources
  ! 'stop' 
  !    if does not exist => stop the execution

  logical             :: got_srcfiles_inv = .false.

  integer             :: src_swap=0

  !- arrays for the diurnal cycle of emission
  type, public ::  cycle_emission
     double precision, pointer, dimension(:,:) :: dcnorma_inv,   &
          emission_rate, &
          emission_rate_NOX
  end type cycle_emission
  real, parameter                                   :: emiss_cycle_time = 86400.
  logical                                           :: emiss_cycle_alloc = .false.
  type(cycle_emission), allocatable, dimension(:,:) :: emiss_cycle

  private

  type GlobalEmissData
     real, pointer :: sc_src(:,:)        ! chem1
     real, pointer :: src_dummy_2d(:,:)  ! chem1
     real, pointer :: flam_frac(:,:)     ! plume_chem1
     real, pointer :: fire_size(:,:)     ! plume_chem1
     real, pointer :: plum_heigth(:,:)   ! volc_chem1
     real, pointer :: vent_elev(:,:)     ! volc_chem1
     real, pointer :: duration (:,:)     ! volc_chem1
     real, pointer :: mean_frp (:,:)     ! plume_chem1
     real, pointer :: std_frp  (:,:)     ! plume_chem1
     real, pointer :: mean_size(:,:)     ! plume_chem1
     real, pointer :: std_size (:,:)     ! plume_chem1
  end type GlobalEmissData

  type(GlobalEmissData), allocatable, target :: oneGlobalEmissData(:)

  integer,parameter :: nucle = 1 ! nucleation mode
  integer,parameter :: accum = 2 ! accumulation mode
  integer,parameter :: coarse = 3 ! coarse mode
  integer,parameter :: aer_bburn=002
  integer,parameter :: urban=003
  integer,parameter :: bioge=004
  integer,parameter :: marin=005
  integer,parameter :: v_ash=006

  public :: read_sourcemaps,                & ! Subroutine
       sources,                        & ! Subroutine
       init_actual_time_index,         & ! Subroutine
       alloc_emiss_cycle,              & ! Subroutine
       get_diurnal_cycle_normalized,   & ! Subroutine
       oneGlobalEmissData,             & ! Type(GlobalEmissData)
       StoreNamelistFileAtChemSources, & ! Subroutine
       alloc_GlobalEmissData,          & ! Subroutine
       nullify_GlobalEmissData,        & ! Subroutine
       dealloc_GlobalEmissData,        & ! Subroutine
       def_proc_src,                   &
       srcmapfn,                       &
       emiss_cycle_alloc,              &
       emiss_cycle_time,               &
       emiss_cycle

contains



  !---------------------------------------------------------------------------- 
  subroutine alloc_GlobalEmissData(global, n2, n3)
    implicit none
    type (GlobalEmissData), intent (inout) :: global
    integer,                intent(in)     :: n2
    integer,                intent(in)     :: n3

    allocate (global%sc_src(n2,n3))
    allocate (global%src_dummy_2d(n2,n3))
    allocate (global%flam_frac(n2,n3))
    allocate (global%fire_size(n2,n3))
    allocate (global%plum_heigth(n2,n3))
    allocate (global%vent_elev(n2,n3))
    allocate (global%duration(n2,n3))
    allocate (global%mean_frp(n2,n3))
    allocate (global%std_frp(n2,n3))
    allocate (global%mean_size(n2,n3))
    allocate (global%std_size(n2,n3))

  end subroutine alloc_GlobalEmissData



  !---------------------------------------------------------------------------- 
  subroutine nullify_GlobalEmissData(global)
    implicit none
    type (GlobalEmissData), intent (inout) :: global

    nullify (global%sc_src)
    nullify (global%src_dummy_2d)
    nullify (global%flam_frac)
    nullify (global%fire_size)
    nullify (global%plum_heigth)
    nullify (global%vent_elev)
    nullify (global%duration)
    nullify (global%mean_frp )
    nullify (global%std_frp  )
    nullify (global%mean_size)
    nullify (global%std_size )

  end subroutine nullify_GlobalEmissData



  !---------------------------------------------------------------------------- 
  subroutine dealloc_GlobalEmissData(global)
    implicit none
    type (GlobalEmissData), intent (inout) :: global

    if (associated(global%sc_src))      deallocate (global%sc_src)
    if (associated(global%src_dummy_2d)) deallocate (global%src_dummy_2d)
    if (associated(global%flam_frac))   deallocate (global%flam_frac)
    if (associated(global%fire_size))   deallocate (global%fire_size)
    if (associated(global%plum_heigth)) deallocate (global%plum_heigth)
    if (associated(global%vent_elev))   deallocate (global%vent_elev)
    if (associated(global%duration))   deallocate (global%duration)
    if (associated(global%mean_frp))   deallocate (global%mean_frp)
    if (associated(global%std_frp))   deallocate (global%std_frp)
    if (associated(global%mean_size))   deallocate (global%mean_size)
    if (associated(global%std_size))   deallocate (global%std_size)

  end subroutine dealloc_GlobalEmissData



  !--(DMK-BRAMS 5.0)-----------------------------------------------------------
  ! Somente para BRAMS 5.0, 1 grade
  subroutine alloc_emiss_cycle(mxp,myp,ngrids,nsrc)

    ! node_mod
    integer , intent(IN) :: mxp
    integer , intent(IN) :: myp

    ! mem_grid
    integer , intent(IN) :: ngrids

    ! mem_chem1
    integer , intent(IN) :: nsrc

    integer ng,isrc

    allocate (emiss_cycle(nsrc,ngrids))
    do ng=1,ngrids
       do isrc=1,nsrc      
          allocate(emiss_cycle(isrc,ng)%dcnorma_inv(mxp,myp) )    
          emiss_cycle(isrc,ng)%dcnorma_inv = 0.
          allocate(emiss_cycle(isrc,ng)%emission_rate(mxp,myp) )    
          emiss_cycle(isrc,ng)%emission_rate = 0.
          allocate(emiss_cycle(isrc,ng)%emission_rate_NOX(mxp,myp) )    
          emiss_cycle(isrc,ng)%emission_rate_NOX = 0.
       enddo
    enddo

    emiss_cycle_alloc=.true.

  end subroutine alloc_emiss_cycle


  !---------------------------------------------------------------------------- 
  subroutine read_sourcemaps(ng,m1,m2,m3,ia,iz,ja,jz, &
       time,iyear1,imonth1,idate1,itime1,&
       ngrids,timmax,chem_nspecies,spc_chem_alloc,   &
       src,off,nsrc,nvert_src,chem1_src_g,bburn,     & 
       antro, bioge,geoge,                           &
       spc_chem_name,on,chemical_mechanism,          &
       emiss_ajust,co,aer_nspecies,spc_aer_alloc,    &
       spc_aer_name,src_name,      &
       chemistry,ntimes_src,aer1_g,nmodes,aerosol,   &
       plumerise,nveg_agreg,plume_mean_g,nzpmax,dzt, &
       rtgt,topt,transport,plume_g,        &
       tropical_forest,boreal_forest,savannah,       &
       grassland,diur_cycle,volcanoes,volc_mean_g,   &
       dn0,zt,zm,mchnum,master_num,mass_bin_dist,CO2,ISFCL,&
       aerosol_mechanism,plume_fre_g,emiss_ajust_aer, oneControlVars)

    ! original
    integer , intent(IN) :: ng
    integer , intent(IN) :: m1
    integer , intent(IN) :: m2
    integer , intent(IN) :: m3
    integer , intent(IN) :: ia
    integer , intent(IN) :: iz
    integer , intent(IN) :: ja
    integer , intent(IN) :: jz
    integer , intent(IN) :: isfcl,CO2  !DSM

    ! grid_dims
    integer, intent(IN) :: nzpmax

    ! mem_grid
    real    , intent(IN) :: time
    integer , intent(IN) :: iyear1
    integer , intent(IN) :: imonth1
    integer , intent(IN) :: idate1
    integer , intent(IN) :: itime1
    integer , intent(IN) :: ngrids
    real    , intent(IN) :: timmax
    real    , intent(IN) :: dzt(nzpmax)
    real    , intent(IN) :: zt(nzpmax)
    real    , intent(IN) :: zm(nzpmax)
    real    , intent(IN) :: rtgt(m2,m3)
    real    , intent(IN) :: topt(m2,m3)
    real    , intent(IN) :: dn0(m1,m2,m3)

    ! chem1_list
    integer          , intent(IN) :: chem_nspecies
    integer          , intent(IN) :: spc_chem_alloc(6,chem_nspecies)
    integer          , intent(IN) :: src
    character(LEN=8) , intent(IN) :: spc_chem_name(chem_nspecies)
    integer          , intent(IN) :: on
    integer          , intent(IN) :: off
    character(LEN=24), intent(IN) :: chemical_mechanism
    real             , intent(IN) :: emiss_ajust(chem_nspecies)
    integer          , intent(IN) :: CO
    integer          , intent(IN) :: transport

    ! mem_chem1
    integer              , intent(IN)    :: nsrc
    integer              , intent(IN)    :: nvert_src(nsrc,ngrids)
    type(chem1_src_vars) , intent(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies,ngrids)
    integer              , intent(IN)    :: bburn,antro, bioge,  geoge
    character(LEN=20)    , intent(IN)    :: src_name(nsrc)
    integer              , intent(IN)    :: chemistry
    integer              , intent(IN)    :: ntimes_src(nsrc)
    integer              , intent(IN)    :: diur_cycle(nsrc)

    ! aer1_list
    integer          , intent(IN) :: aer_nspecies
    integer          , intent(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    character(LEN=8) , intent(IN) :: spc_aer_name(nmodes,aer_nspecies)
    integer          , intent(IN) :: nmodes
    real             , intent(INOUT) :: mass_bin_dist(nmodes)
    character(LEN=24), intent(IN) :: aerosol_mechanism
    real             , intent(IN) :: emiss_ajust_aer(nmodes,aer_nspecies) 

    ! mem_aer1
    type (aer1_vars) , intent(INOUT) :: aer1_g(nmodes,aer_nspecies,ngrids)
    integer          , intent(IN)    :: aerosol

    ! mem_plume_chem1
    integer               , intent(IN)    :: nveg_agreg
    integer               , intent(IN)    :: plumerise
    type (plume_mean_vars), intent(INOUT) :: plume_mean_g(nveg_agreg,ngrids)
    type (plume_fre_vars) , intent(INOUT) :: plume_fre_g(5,ngrids)
    type (plume_vars)     , intent(INOUT) :: plume_g(nveg_agreg,chem_nspecies,ngrids)
    integer               , intent(IN)    :: tropical_forest
    integer               , intent(IN)    :: boreal_forest
    integer               , intent(IN)    :: savannah
    integer               , intent(IN)    :: grassland

    ! mem_volc_chem1
    type (volc_mean_vars), intent(INOUT) :: volc_mean_g(ngrids)
    integer              , intent(IN)    :: volcanoes    

    ! node_mod
    integer, intent(IN)   :: mchnum
    integer, intent(IN)   :: master_num

    type(ControlVars), pointer, intent(in) :: oneControlVars
    
    !- local var
    integer :: iunit,isrc,k1,k2,isrctime,iwtdo
    character(len=256)  :: fname
    character(len=2)   :: cgrid


    if(srcmapfn(1:LEN_trim(srcmapfn)) == 'NONE' .or. srcmapfn(1:LEN_trim(srcmapfn)) == 'none') return

    if( .not. got_srcfiles_inv ) then 

       ! Out: nsrcfiles(maxgrds), 
       !      fnames_src(maxsrcfiles,maxgrds), 
       !      itotdate_src(maxsrcfiles,maxgrds),
       !      src_times(maxsrcfiles,maxgrds), 
       !      next_srcfile(maxgrds), 
       !      srctime1, 
       !      srctime2, 
       !      got_srcfiles_inv

       call src_file_inv(srcmapfn,iyear1,imonth1,idate1,itime1,ngrids,time,timmax,nsrc,diur_cycle, &
            mchnum,master_num)

       ! Out: actual_time_index(max_ntimes_src,nsrc)

       call init_actual_time_index(nsrc,ntimes_src)

    endif

    ! In:  next_srcfile(maxgrds), 
    !      nsrcfiles(maxgrds)
    !
    ! Out: iwtdo, 
    !      srctime1, 
    !      srctime2

    !- check if files required does exist and decide what to do in case not
    if(mchnum==master_num) call check_src_files(next_srcfile(ng),nsrcfiles(ng),iwtdo)

    !PRINT *, 'Aqui 1: LFR',mchnum,iwtdo; call flush(6)
    call Broadcast(iwtdo, master_num, "iwtdo")
    !PRINT *, 'Aqui 2: LFR',mchnum,iwtdo;call flush(6)


    if(iwtdo == 0) return 

    call Broadcast(nsrcfiles, master_num, "nsrcfiles")

    !- loop at source files (isrctime will be > 1, only if will use linterp for any source)
    do isrctime=1,maxval(ntimes_src(1:nsrc))


       !- swap sources:  copy time level 2->1 and read the next data for the time level 2
       if(src_swap == 1 .and. maxval(ntimes_src(1:nsrc)) > 1 .and. isrctime==1) then 

          ! Out: chem1_src_g
          call swap_sources(m1,m2,m3,time,chem_nspecies,spc_chem_alloc, &
               src,off,nsrc,nvert_src(:,ng),chem1_src_g(:,:,:,ng),bburn,geoge)

          next_srcfile(ng) = next_srcfile(ng) + 1 
          cycle
       endif


       !- next file to open
       fname=fnames_src(next_srcfile(ng),ng)

       ! InOut: chem1_src_g, 
       !        aer1_g,
       !        plume_mean_g
       !        volc_mean_g

       !- read emission dataset using V-format
       call read_sources_vfm(ng,m1,m2,m3,iyear1,imonth1,idate1,isrctime,    &
            fname(1:LEN_trim(fname)),chem_nspecies,        &
            spc_chem_alloc,spc_chem_name,src,on,off,       &
            chemical_mechanism,emiss_ajust,co,aer_nspecies,&
            spc_aer_alloc,spc_aer_name,urban,nucle,accum,  &
            nsrc,chem1_src_g(:,:,:,ng),src_name,chemistry,aer1_g(:,:,ng), &
            nmodes,aerosol,plumerise,nveg_agreg,plume_mean_g(:,ng),   &
            volc_mean_g(ng),volcanoes,mchnum,master_num,mass_bin_dist,v_ash, &
            aerosol_mechanism,plume_fre_g(:,ng),emiss_ajust_aer, oneControlVars)

       !- split bburn emissions into flaming/smoldering parts
       if(plumerise /= 0) &
            call emis_flam_smold(m1,m2,m3,isrctime,                   &
            nsrc,chem1_src_g(:,:,:,ng),bburn,chem_nspecies, &
            spc_chem_alloc,src,on,off,transport,  &
            aer1_g(:,:,ng),aerosol,aer_nspecies,          &
            
            spc_aer_alloc,nmodes,aer_bburn,       &
            
            plume_mean_g(:,ng),plume_g(:,:,ng),tropical_forest, &
            boreal_forest,savannah,grassland,     &
            nveg_agreg,spc_aer_name,plume_fre_g(:,ng),plumerise)

       !-----
       !- to use the mass conservation fix
       !- convert from [kg/m^2]to mixing ratio expressed in [ppbm = 1.e9 kg/kg]
       !call convert_to_mixing_ratio(ng,m1,m2,m3,isrctime)

       !- convert from [kg/m^2]to tracer density expressed in [ 1.e9 kg/m3]
       call convert_to_tracer_density(m1,m2,m3,ia,iz,ja,jz,isrctime,        &
            nzpmax,dzt,rtgt,nsrc,                 &
            chem1_src_g(:,:,:,ng),bburn,chem_nspecies,  &
            spc_chem_alloc,src,off,transport,     &
            aer1_g(:,:,ng),aerosol,aer_nspecies,  &
            spc_aer_alloc,nmodes,aer_bburn,geoge, &
            volcanoes,bioge,CO2,ISFCL,spc_aer_name) 

       !- for now volcanic emissions only works with SIMPLE aerosol model
       if(AEROSOL == 1 .and. volcanoes == 1) then
          call vert_dist_of_volcanic_emission(m1,m2,m3,ia,iz,ja,jz,isrctime,&
               nzpmax,dzt,zt,zm,rtgt,topt,dn0,nsrc,&
               chem1_src_g(:,:,:,ng),geoge,chem_nspecies,  &
               spc_chem_alloc,src,off,transport, &
               aer1_g(:,:,ng),aerosol,aer_nspecies,&
               spc_aer_alloc,nmodes,volc_mean_g(ng), v_ash)
       endif
       !-----

       !- next time (don't change the position of this line)
       if(isrctime==1) next_srcfile(ng) = next_srcfile(ng) + 1 

    enddo

    if(mchnum==master_num) then
       !- update srctime for the next time of reading/update    
       srctime1 =   src_times(next_srcfile(ng)-1,1)
       srctime2 =   src_times(next_srcfile(ng)  ,1)
    end if

    !CALL Broadcast(srctime1, master_num, "srctime1") !LFR
    call Broadcast(srctime2, master_num, "srctime2") 
    !******* debug *****
    !write(*,fmt='(I3.3,1X,A,F12.6,1X,F12.6)') &
    !    mchnum,'LFR:srctime1 and 2',srctime1,srctime2
    !CALL flush(6)
    !*******************
    if(ng==ngrids) src_swap = 1


    !PRINT*,'next_srcfile,srctime1-2=',next_srcfile(ng),srctime1,srctime2
    !CALL flush(6)


  end subroutine read_sourcemaps




  !----------------------------------------------------------------------
  subroutine src_file_inv(srcpref,iyear1,imonth1,idate1,itime1,ngrids, &
       time,timmax,nsrc,diur_cycle,mchnum,master_num)

    ! original
    character(len=*) , intent(IN) :: srcpref
    integer          , intent(IN) :: iyear1
    integer          , intent(IN) :: imonth1
    integer          , intent(IN) :: idate1
    integer          , intent(IN) :: itime1
    integer          , intent(IN) :: ngrids
    real             , intent(IN) :: time
    real             , intent(IN) :: timmax

    ! mem_chem1
    integer , intent(IN) :: nsrc
    integer , intent(IN) :: diur_cycle(nsrc)

    ! node_mod
    integer, intent(IN) :: mchnum
    integer, intent(IN) :: master_num

    integer :: nc,nf,lnf,nvftot, ng,it,isrc,nfstart
    integer :: inyear,inmonth,indate,inhour
    integer :: index
    real    :: localInc
    logical :: there

    character(len=256), dimension(maxsrcfiles) :: fnames
    character(len=256) :: vpref
    character(len=14)  :: itotdate,itotdate_current
    character(len=2)   :: cgrid

    integer                :: lenFnames_src
    integer                :: lenItot_src
    integer                :: sizeCharVec
    integer                :: sizeIntVec
    integer                :: sizeRealVec
    integer                :: lastChar
    integer                :: ierr2
    character(len=8)       :: c0, c1
    integer,   allocatable :: intVec(:)
    character, allocatable :: charVec(:)
    real,      allocatable :: realVec(:)
    real(kind=r8)          :: secs_init,secs_src
    character(len=*), parameter :: h="**(src_file_inv)**"
    !    REAL(kind=8) :: secs_init,secs_src ! CCATT-BRAMS 4.3

    character(len=256) :: sVarName


    sizeIntVec     = 2*ngrids
    lenFnames_src  = len(fnames_src)
    lenItot_src    = len(itotdate_src)


    if (mchnum==master_num) then

       do ng=1,ngrids

          fnames(1:maxsrcfiles)= 'XXXXXXXXXXXXXXXX'

          !Get abs seconds of run start
          call date_abs_secs2(iyear1,imonth1,idate1,itime1*100,secs_init)

          ! get the current time  
          call date_add_to(iyear1,imonth1,idate1,itime1*100  &
               ,time,'s',inyear,inmonth,indate,inhour)

          call date_make_big (inyear,inmonth,indate,inhour,itotdate_current)

          ! Go through src files and make inventory

          nc=LEN_trim(srcpref)
          nvftot=-1
          vpref=srcpref

          write(cgrid,'(a1,i1)') 'g',ng

          !### TEMPORARY CALL SYSTEM ALTERNATIVE ##### RMF 
          !RMF: 
          !with this struct you can only select between daily or hourly emissions.

          index = 0
          there = .false.

          select case(minval(diur_cycle))
          case(0)
             nvftot = ceiling(((timmax/3600))) + 1 ; localInc = 3600.   
          case(1)
             nvftot = ceiling(((timmax/3600))/24.) ; localInc = 86400.   
          end select

          do nf = 1, nvftot

             call makefnam (sVarName,srcpref,0,inyear,inmonth,indate,inhour,'T',cgrid,'vfm')
             print*,'Looking for sources files -->: ',trim(sVarName)
             inquire(file=trim(sVarName),exist=there)

             if (there) then  
                index = index + 1 
                fnames(index) = trim(sVarName)           
             endif

             if(nvftot > maxsrcfiles .or. index > maxsrcfiles) then
                call fatal_error('Too many sources files')
             endif

             call date_add_to(inyear,inmonth,indate,inhour  &
                  ,localInc,'s',inyear,inmonth,indate,inhour)
          end do

          nvftot = index

          if (nvftot .eq. 0) then
             call fatal_error('Sources files not found!')
          endif

          !### TEMPORARY CALL SYSTEM ALTERNATIVE #####  RMF

          nsrcfiles(ng)=0
          do nf=1,nvftot
             lnf=LEN_trim(fnames(nf))
             !print*,lnf,fnames(nf)

             read(fnames(nf)(lnf-23:lnf-6),20) inyear,inmonth,indate,inhour
20           format(i4,1x,i2,1x,i2,1x,i6)

             call date_make_big(inyear,inmonth,indate,inhour,itotdate)

             nsrcfiles(ng)=nsrcfiles(ng)+1
             fnames_src(nsrcfiles(ng),ng)=fnames(nf)
             itotdate_src(nsrcfiles(ng),ng)=itotdate

             call date_abs_secs2(inyear,inmonth,indate,inhour,secs_src)
             src_times(nsrcfiles(ng),ng)=secs_src - secs_init

          enddo

          call RAMS_dintsort(nsrcfiles(ng),itotdate_src(:,ng),fnames_src(:,ng))

          !  start printing section
          !--------------------------------------------------------------

          print*,' '
          print*,' '
          print*,' '
          print*,'-------------------------------------------------------------'
          print*,'-----------  Sources File Input Inventory --for --- GRID=', ng
          print*,'-------------------------------------------------------------'
          do nf=1,nsrcfiles(ng)
             print 8,  nf, itotdate_src(nf,ng),src_times(nf,ng) ,trim(fnames_src(nf,ng))
          enddo
8         format(i4,1x,a16,1x,f10.0,2x,a)
          print*,'------------------------------------------------------'

       enddo ! ngrids

       if(ngrids > 1 .and. int((sum(nsrcfiles(1:ngrids)))/ nsrcfiles(1)) .ne. ngrids) then
          call fatal_error("The number of src files for each grid, MUST be the same")
       endif

       !- Are there enough src files available(for now only 1 grid is considered)

       ng=1
       if(src_times(nsrcfiles(ng),ng) < timmax) then 
          print*,'=============================================================================='
          print*,'Warning:'
          print*,'Not enough source files for the entire time integration were found, model will continue.'
          print*,'=============================================================================='
       endif

       !- perform some initializations (for now only 1 grid is considered)
       ng = 1 
       next_srcfile(1:ngrids) = 0
       do nf=1,nvftot
          if(itotdate_src(nf,ng) == itotdate_current) then
             next_srcfile(1:ngrids) = nf
             exit
          elseif(itotdate_src(nf,ng) > itotdate_current) then
             next_srcfile(1:ngrids) = nf-1
             exit
          endif
       enddo

       if (next_srcfile(ng) < 1) then
          call fatal_error(' next_srcfile < 1 ')
       endif

       if (next_srcfile(ng) > nvftot-1 .and. sum(diur_cycle(1:nsrc)) < 4 ) then
          call fatal_error('next_srcfile > nvftot-1')
       endif

       srctime1     =   src_times(next_srcfile(ng)  ,ng)
       srctime2     =   src_times(next_srcfile(ng)+1,ng)

       if (next_srcfile(ng) > nvftot-1 .and. sum(diur_cycle(1:nsrc)) == 4 ) &
            srctime2 = srctime1 + 86400. 

       !- fill src_times arrays above nvftot with valid numbers (to use in 
       !- case of forecast or not more available data)
       do ng=1,ngrids
          do nf=nsrcfiles(ng)+1,maxsrcfiles
             src_times(nf,ng) = src_times(nf-1,ng) + (srctime2-srctime1)
          enddo
       enddo

       got_srcfiles_inv = .true.

       ng=1
       print*,'next_srcfile=',next_srcfile(ng),srctime1,srctime2,itotdate_current
       print*,'next_srcfile= ',trim(fnames_src(next_srcfile(ng),ng))
       call flush(6)

    end if ! IF (mchnum == master_num)

    ! allocate 'int' broadcast area
    allocate(intVec(sizeIntVec), stat=ierr2)
    if (ierr2/=0) then
       write(c0,"(i8)") ierr2
       write(c1,"(i8)") sizeIntVec
       call fatal_error(h//" allocate intVec("//trim(adjustl(c1))// &
            ") fails with stat="//trim(adjustl(c0)))
    end if

    ! master process gathers data for broadcasting

    if (mchnum==master_num) then
       intVec(1:ngrids) = next_srcfile(1:ngrids)
       intVec(ngrids+1:sizeIntVec) = nsrcfiles(1:ngrids)
    end if

    ! broadcast integer data to remaining processes
    ! envia: nsrcfiles(1:ngrids), 
    !        next_srcfile(1:ngrids)
    call Broadcast(intVec, master_num, "intVec")

    ! scatter broadcasted data
    next_srcfile(1:ngrids) = intVec(1:ngrids)
    nsrcfiles(1:ngrids) = intVec(ngrids+1:sizeIntVec)

    ! deallocate broadcast area
    deallocate(intVec, stat=ierr2)
    if (ierr2/=0) then
       write(c0,"(i8)") ierr2
       call fatal_error(h//" deallocate intVec fails with stat="//&
            trim(adjustl(c0)))
    end if

    ! allocate broadcast area
    sizeCharVec = sum(nsrcfiles(1:ngrids)*(lenFnames_src + lenItot_src))
    allocate(charVec(sizeCharVec), stat=ierr2)
    if (ierr2/=0) then
       write(c0,"(i8)") ierr2
       write(c1,"(i8)") sizeCharVec
       call fatal_error(h//" allocate charVec("//trim(adjustl(c1))//&
            ") fails with stat="//trim(adjustl(c0)))
    end if

    ! master process prepares broadcast data

    if (mchnum==master_num) then
       lastChar = 0
       do ng=1,ngrids
          do nf=1,nsrcfiles(ng)
             do nc=1,lenFnames_src
                charVec(lastChar+nc) = fnames_src(nf,ng)(nc:nc)
             end do
             lastChar = lastChar + lenFnames_src
             do nc=1,lenItot_src
                charVec(lastChar+nc) = itotdate_src(nf,ng)(nc:nc)
             end do
             lastChar = lastChar + lenItot_src
          end do
       end do
    end if

    ! broadcast character data to remaining processes
    ! envia fnames_src(1:nsrcfiles(ng),1:ngrids),  
    !       itotdate_src(1:nsrcfiles(ng),1:ngrids)
    call Broadcast(charVec, master_num, "charVec")

    ! scatter broadcasted data
    lastChar=0
    do ng=1,ngrids
       do nf=1,nsrcfiles(ng)
          do nc=1,lenFnames_src
             fnames_src(nf,ng)(nc:nc) = charVec(lastChar+nc)
          end do
          lastChar = lastChar + lenFnames_src
          do nc=1,lenItot_src
             itotdate_src(nf,ng)(nc:nc) = charVec(lastChar+nc)
          end do
          lastChar = lastChar + lenItot_src
       end do
    end do

    ! deallocate broadcast area
    deallocate(charVec, stat=ierr2)
    if (ierr2/=0) then
       write(c0,"(i8)") ierr2
       call fatal_error(h//" deallocate charVec fails with stat="//&
            trim(adjustl(c0)))
    end if
    got_srcfiles_inv = .true.   !Modificacao do Massaru


  end subroutine src_file_inv





  !--------------------------------------------------------------

  subroutine init_actual_time_index(nsrc,ntimes_src)

    integer , intent(IN) :: nsrc
    integer , intent(IN) :: ntimes_src(nsrc)


    integer :: it

    !- index to control memory access of src arrays, because some
    !- arrays have 2 time levels and others only 1 time level
    !- for time 1, the memory position is always allocated for all sources types
    it = 1
    actual_time_index(it,1:nsrc) = 1
    !- for the second, will depend on if linterp is wanted or not.
    !- if linterp is not desired for any sources, actual_time_index = 1
    !- and will have actually one only memory position 
    it = 2
    actual_time_index(it,1:nsrc) = ntimes_src(1:nsrc) 

  end subroutine init_actual_time_index




  !--------------------------------------------------------------

  subroutine check_src_files(next_srcfile,nsrcfiles,iwtdo)

    ! original
    integer , intent(INOUT) :: next_srcfile
    integer , intent(IN)    :: nsrcfiles
    integer , intent(INOUT) :: iwtdo

    iwtdo = 1

    !-check if is not greater the max number defined
    if(next_srcfile+1 > maxsrcfiles) then
       call fatal_error("next_srcfile(ng)+1 > maxsrcfiles")
    endif
    !PRINT*,' next_srcfile(ng) , nsrcfiles(ng)',next_srcfile, nsrcfiles
    !CALL flush(6)

    if(next_srcfile > nsrcfiles) then 
       !    IF(next_srcfile(ng)+1 > nsrcfiles(ng)) THEN 

       !- situation 1 : stop model execution
       if(def_proc_src(1:LEN_trim(def_proc_src)) == 'STOP' ) then
          call fatal_error("Not src files available!")

          !- situation 2 :keep the current sources
       elseif(def_proc_src(1:LEN_trim(def_proc_src)) == 'LAST_SOURCES' )then

          iwtdo=0

          next_srcfile=next_srcfile+1

          !- update srctime for the next time of reading/update    
          srctime1 =   src_times(next_srcfile-1,1)
          srctime2 =   src_times(next_srcfile  ,1)

          !print*,'-----------------------------------------------------'
          print*,'Not src files available:'
          print*,'using previous day sources, model will continue ...'
          print*,'srctime1 and 2=',srctime1,srctime2
          call flush(6)

          return
       endif

    endif


  end subroutine check_src_files




  !--------------------------------------------------------------

  subroutine swap_sources(m1,m2,m3,time,chem_nspecies,spc_chem_alloc, &
       src,off,nsrc,nvert_src,chem1_src_g,bburn,geoge)

    ! original
    integer , intent(IN) :: m1
    integer , intent(IN) :: m2
    integer , intent(IN) :: m3
    real    , intent(IN) :: time

    ! chem1_list
    integer , intent(IN) :: chem_nspecies
    integer , intent(IN) :: spc_chem_alloc(6,chem_nspecies)
    integer , intent(IN) :: src
    integer , intent(IN) :: off

    ! mem_chem1
    integer             , intent(IN)    :: nsrc
    integer             , intent(IN)    :: nvert_src(nsrc)
    type(chem1_src_vars), intent(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    integer             , intent(IN)    :: bburn,geoge

    integer ispc,isrc
    integer, parameter :: it1=1, & ! time level 1
         it2=2    ! time level 2

    do isrc=1,nsrc  

       !- if the time level 2 uses the same memory allocation area
       !- of time level 1 => nothing to do
       if(actual_time_index(it2,isrc) == 1) cycle
       if(isrc == bburn) stop 444
       if(isrc == geoge) stop 444

       !- else: make the swap     
       do ispc=1,chem_nspecies
          if(spc_chem_alloc(src,ispc) == off) cycle

          chem1_src_g(it1,isrc,ispc)%sc_src(1:nvert_src(isrc),1:m2,1:m3) = &
               chem1_src_g(it2,isrc,ispc)%sc_src(1:nvert_src(isrc),1:m2,1:m3) 
       enddo
    enddo

    !- aerosol section still need to be done

    print*,'--> source swapped done at time (h)=',time/3600.


  end subroutine swap_sources





  !--------------------------------------------------------------

  subroutine read_sources_vfm(ng,m1,m2,m3,iyear,imon,iday,isrctime,fname, &
       chem_nspecies,spc_chem_alloc,spc_chem_name, &
       src,on,off,chemical_mechanism,emiss_ajust,  &
       co,aer_nspecies,spc_aer_alloc,spc_aer_name, &
       urban,nucle,accum,nsrc,chem1_src_g,src_name,&
       chemistry,aer1_g,nmodes,aerosol,plumerise,  &
       nveg_agreg,plume_mean_g,volc_mean_g,volcanoes, &
       mchnum,master_num,mass_bin_dist,v_ash,aerosol_mechanism,&
       plume_fre_g,emiss_ajust_aer, oneControlVars)


    use mem_grid, only: grid_g

    ! original
    integer       , intent(IN) :: ng
    integer       , intent(IN) :: m1
    integer       , intent(IN) :: m2
    integer       , intent(IN) :: m3
    integer       , intent(IN) :: iyear
    integer       , intent(IN) :: imon
    integer       , intent(IN) :: iday
    integer       , intent(IN) :: isrctime
    character*(*) , intent(IN) :: fname

    ! chem1_list
    integer          , intent(IN) :: chem_nspecies
    integer          , intent(IN) :: spc_chem_alloc(6,chem_nspecies)
    character(LEN=8) , intent(IN) :: spc_chem_name(chem_nspecies)
    integer          , intent(IN) :: src
    integer          , intent(IN) :: on
    integer          , intent(IN) :: off
    character(LEN=24), intent(IN) :: chemical_mechanism
    real             , intent(IN) :: emiss_ajust(chem_nspecies)
    integer          , intent(IN) :: CO

    ! aer1_list
    integer          , intent(IN) :: aer_nspecies
    integer          , intent(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    character(LEN=8) , intent(IN) :: spc_aer_name(nmodes,aer_nspecies)
    integer          , intent(IN) :: urban
    integer          , intent(IN) :: nucle
    integer          , intent(IN) :: accum
    integer          , intent(IN) :: nmodes
    integer          , intent(IN) :: v_ash
    real             , intent(INOUT) :: mass_bin_dist(nmodes)
    character(LEN=24), intent(IN) :: aerosol_mechanism
    real             , intent(IN) :: emiss_ajust_aer(nmodes,aer_nspecies)

    ! mem_chem1
    integer              , intent(IN)    :: nsrc
    type(chem1_src_vars) , intent(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    character(LEN=20)    , intent(IN)    :: src_name(nsrc)
    integer              , intent(IN)    :: chemistry

    ! mem_aer1
    type (aer1_vars) , intent(INOUT) :: aer1_g(nmodes,aer_nspecies)
    integer          , intent(IN)    :: aerosol
    integer          , intent(IN)    :: plumerise

    ! mem_plume_chem1
    integer               , intent(IN)    :: nveg_agreg
    type (plume_mean_vars), intent(INOUT) :: plume_mean_g(nveg_agreg)
    type (plume_fre_vars) , intent(INOUT) :: plume_fre_g(5)

    ! mem_volc_chem1
    type (volc_mean_vars), intent(INOUT) :: volc_mean_g
    integer              , intent(IN)    :: volcanoes    
    type(ControlVars), pointer, intent(in) :: oneControlVars
    
    ! node_mod
    integer, intent(in) :: mchnum
    integer, intent(in) :: master_num

    !- local var
    integer :: iunit,ispc,isrc,iveg_agreg,nvert,ihour,itim
    character(len=20) read_spc_name,read_src_name,section,read_units,date
    character(len=20) read_mean,read_veg_name
    integer read_ident_chem_mec,read_ident_src,dummy,read_ident_veg
    integer read_ident_aer,read_aer_mode, imode
    integer nspecies,nxp,nyp,i,j,ii,jj,iax,ipx,ivx
    real dep_glon(2), dep_glat(2)
    !    REAL, ALLOCATABLE, DIMENSION(:,:) :: src_dummy_2d
    real, pointer, dimension(:,:) :: src_dummy_2d
    character(len=32) :: chemical_mechanism_test,aerosol_mechanism_test
    logical :: there

    integer :: exDo


    integer                :: nc
    integer                :: ierr2
    character(len=8)       :: c0, c1, c2
    integer                :: sizeIntVec
    integer,   allocatable :: intVec(:)
    real, pointer          :: dummy_sc_src(:,:)
    integer                :: lastChar
    integer                :: sizeCharVec
    character, allocatable :: charVec(:)
    character(len=*), parameter :: h="**(read_sources_vfm)**"

    character(len=64) :: cdummy
    character(len=12) :: type_volc_process

    integer :: recn,recordLen

    !- initialization of local vars
    type_volc_process="XXXXXXXX"
    read_spc_name="XXXXXXXX"
    read_src_name="XXXXXXXX"
    section      ="XXXXXXXX"
    read_units   ="XXXXXXXX"
    date         ="XXXXXXXX"
    read_mean    ="XXXXXXXX"
    read_veg_name="XXXXXXXX"
    iax=0;ipx=0;ivx=0
    !- initial attributions/allocations
    nspecies = chem_nspecies + aer_nspecies*nmodes

    ihour = 0
    nvert = 1
    iunit = 2
    !- This routine does not allow negative fluxes.
    allocate (src_dummy_2d(m2,m3));src_dummy_2d=0.

    !################### LFR
    !  recordLen=4*m2*m3
    !  open(unit=33,file='aer.gra',&
    !      action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
    !      recl=recordLen)
    !  recn=1
    !
    !
    !################### LFR

    !-need to zerout aerosol sources, before reading the emission of the new day
    if(AEROSOL > 0) then
       mass_bin_dist = 0.0
       do ispc=1,aer_nspecies
          do imode=1,nmodes
             if(spc_aer_alloc(src,imode,ispc) == on) aer1_g(imode,ispc)%sc_src = 0.
          enddo
       enddo
    endif

    if (mchnum==master_num) then
       print*,'---------------------------------------------------------------------------'
       print*,'opening emission file= ',fname(1:len_trim(fname))
       inquire(file=fname,exist=there)
       if(.not.there) then
          call fatal_error("emission file not found!")
       endif
       !- open the source file
       open(unit=iunit,file=fname,form='formatted',status='old') 
       ispc = 0 ; isrc = 0; iveg_agreg = 0
       read(iunit,*) nxp,(dep_glon(i),i=1,2)
       read(iunit,*) nyp,(dep_glat(i),i=1,2)
       read(iunit,*) date 
       write (*,fmt='(a)') '=== sources header (nxpoints,nypoint,lon,lat,date) ==='
       write (*,fmt='(4x,2(i4.4,1x),4(f8.3,1x),a)') &
            nxp,nyp,(dep_glon(i),i=1,2),(dep_glat(i),i=1,2),date
       !- test if the source data is for the chemical mechanism that will be used:
       read(iunit,*)  chemical_mechanism_test,aerosol_mechanism_test
       if(trim( chemical_mechanism_test ) /=  trim(chemical_mechanism)) then
          call fatal_error("wrong chem mechanism at chem_sources. expected="// &
               trim(chemical_mechanism(1:len_trim(chemical_mechanism)))//" read="// &
               trim(chemical_mechanism_test(1:len_trim(chemical_mechanism_test))))
       else
          print*,'   chem mechanism= ',trim(chemical_mechanism(1:len_trim(chemical_mechanism)))
       endif
       if(trim( aerosol_mechanism_test ) /=  trim(aerosol_mechanism)) then
          call fatal_error("wrong aer mechanism at chem_sources. expected="// &
               trim(aerosol_mechanism(1:len_trim(aerosol_mechanism)))//" read="// &
               trim(aerosol_mechanism_test(1:len_trim(aerosol_mechanism_test))))
       else
          print*,'   aer  mechanism= ',trim(aerosol_mechanism(1:len_trim(aerosol_mechanism)))
          if(trim(aerosol_mechanism) == "matrix" ) then
             print*,"   level of matrix aer model is not checked in the read sources routine"
          endif
       endif
    end if ! (mchnum==master_num)

    do i=1,nspecies*nsrc*5
       if (mchnum==master_num) then
          read(iunit,*,iostat=exdo) section
          if(trim(section) == 'aerosol' .and. aerosol == off .and. iax==0) then
             print *,'warning: aerosol is present in emission file but turned off in ramsin.'
             print *,'the related emission fields will be ignored.'
             iax=1
          end if
          if( (trim(section) == 'plume' .or. trim(section) == 'plumefre') &
               .and. plumerise == 0 .and. ipx==0) then
             print *,'warning: plume information is present in emission file but plumerise' 
             print *,'is turned off in ramsin. the data will be ignored.'
             ipx=1
          end if
          if(trim(section) == 'volcanoes' .and. volcanoes /= on .and. ivx==0) then
             print *,'warning: volcanic emission data is present but turned off in ramsin.'
             print *, 'the data will be ignored.'
             ivx=0
          end if
       end if
       !
       call broadcast(exdo, master_num,"exdo")
       call broadcast(section, master_num, "section")
       !
       if (exdo<0) exit ! eof

       !emission section  --------------------------------

       !%%%%%%%%%%     chemistry section  %%%%%%%%%%%%%
       if(trim(section) == 'chemistry') then
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if (mchnum==master_num) then
             read(iunit,*)   read_spc_name &
                  , read_ident_chem_mec &
                  , read_src_name       &
                  , read_ident_src      &
                  , read_units
             ispc = read_ident_chem_mec
             isrc = read_ident_src
             itim = actual_time_index(isrctime,isrc)
             print*,"   chem spc/src: ", read_spc_name , read_src_name 
          endif
          sizeintvec=3
          !allocate 'int' broadcast area
          allocate(intvec(sizeintvec), stat=ierr2)
          if (ierr2/=0) then
             write(c0,"(i8)") ierr2
             write(c1,"(i8)") sizeintvec
             call fatal_error(h//" allocate intvec("//trim(adjustl(c1))// &
                  ") fails with stat="//trim(adjustl(c0)))
          endif
          ! master process gathers data for broadcasting   
          if (mchnum==master_num) then
             intvec(1) = ispc
             intvec(2) = isrc
             intvec(3) = itim
          endif
          ! broadcast integer data to remaining processes
          ! envia: ispc, isrc, itim
          call broadcast(intvec, master_num, "intvec")
          ! scatter broadcasted data
          ispc = intvec(1)
          isrc = intvec(2)
          itim = intvec(3)
          ! deallocate broadcast area
          deallocate(intvec, stat=ierr2)
          if (ierr2/=0) then
             write(c0,"(i8)") ierr2
             call fatal_error(h//" deallocate intvec fails with stat=" &
                  //trim(adjustl(c0)))
          endif
          if(.not. associated( chem1_src_g(itim,isrc,ispc)%sc_src)) then
             call fatal_error("chem source memory not allocated for specie "//trim(read_spc_name)//&
                  " and source "//trim(read_src_name))
          endif
          allocate(dummy_sc_src(m2,m3), stat=ierr2)
          if (ierr2/=0) then
             write(c0,"(i8)") ierr2
             call fatal_error(h//" allocate dummy_sc_src fails with stat="//trim(adjustl(c0)))
          endif

          call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%sc_src, &
               dummy_sc_src,                         &
               trim(spc_chem_name(ispc))//"_"//trim(src_name(isrc)), oneControlVars)

          where( dummy_sc_src(:,:) < 0.) dummy_sc_src(:,:)  = 0.0
          chem1_src_g(itim,isrc,ispc)%sc_src(1,1:m2,1:m3)=dummy_sc_src(1:m2,1:m3)
          if (mchnum==master_num) then 
             if (mchnum==master_num) write(*,fmt='("    Maxval,MinVal: ",2(E13.6,1X))') & 
                  maxval(chem1_src_g(itim,isrc,ispc)%sc_src(1,1:m2,1:m3)),&
                  minval(chem1_src_g(itim,isrc,ispc)%sc_src(1,1:m2,1:m3))
          endif
          deallocate(dummy_sc_src, stat=ierr2)
          if (ierr2/=0) then
             write(c0,"(i8)") ierr2
             call fatal_error(h//" deallocate dummy_sc_src fails with stat=" &
                  //trim(adjustl(c0)))
          endif
          !- biogenic co only for tracer runs         
          if(chemistry > 0 .and.  ispc == co .and. spc_chem_alloc(3  ,co) == on ) then           
             chem1_src_g(itim,3,co)%sc_src(1,:,:) = 0.
          endif
          !
          !- ajust emissions 
          chem1_src_g(itim,isrc,ispc)%sc_src(1,:,:) = emiss_ajust(ispc) &
               *chem1_src_g(itim,isrc,ispc)%sc_src(1,:,:)
          !srf- especial para o barca - voc x 0.6 somente para urbano:
          !- ajustando emissoes antropicas
          !if(emiss_ajust(ispc) < 1. .and. trim(src_name(isrc))=='antro') then
          !     chem1_src_g(itim,isrc,ispc,ng)%sc_src(1,:,:) = emiss_ajust(ispc)*&
          !     chem1_src_g(itim,isrc,ispc,ng)%sc_src(1,:,:)
          !     print*,'barca voc urbanos=',trim(src_name(isrc)),emiss_ajust(ispc),'x',trim(read_spc_name)
          !endif
          !srf- especial para o barca - end

          !%%%%%%%%%%     aerosol section  %%%%%%%%%%%%%
       elseif(trim(section) == 'aerosol' ) then     
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if (mchnum==master_num) then
             read(iunit,*)   read_spc_name       &
                  , read_ident_aer      &
                  , read_aer_mode       &
                  , read_src_name       &
                  , read_ident_src      &
                  , read_units
             ispc = read_ident_aer
             imode= read_aer_mode
             isrc = read_ident_src
             print*,"   aer  spc/src: ",read_spc_name,read_src_name,ispc,imode    

             if(trim(read_spc_name)=='v_ash1' .and. ispc==v_ash) then
                print*,'-> reading size distr for ', trim(spc_aer_name(imode,ispc)), ' source: ' &
                     , trim(src_name(isrc))
                read(iunit,*)   cdummy
                read(iunit,*)   mass_bin_dist(:)
                print*,'-> reading volcanic mass_bin_dist = ',mass_bin_dist(:)
             endif
          endif
          sizeintvec=3
          ! allocate 'int' broadcast area
          allocate(intvec(sizeintvec), stat=ierr2)
          if (ierr2/=0) then
             write(c0,"(i8)") ierr2
             write(c1,"(i8)") sizeintvec
             call fatal_error(h//" allocate intvec("//trim(adjustl(c1)) &
                  //") fails with stat="//trim(adjustl(c0)))
          endif
          ! master process gathers data for broadcasting   
          if (mchnum==master_num) then
             intvec(1) = ispc
             intvec(2) = imode
             intvec(3) = isrc
          endif
          ! broadcast integer data to remaining processes
          ! envia: ispc, imode, isrc
          call broadcast(intvec, master_num, "intvec")
          ! scatter broadcasted data
          ispc  = intvec(1)
          imode = intvec(2)
          isrc  = intvec(3)
          ! deallocate broadcast area
          deallocate(intvec, stat=ierr2)
          if (ierr2/=0) then
             write(c0,"(i8)") ierr2
             call fatal_error(h//" deallocate intvec fails with stat=" &
                  //trim(adjustl(c0)))
          endif
          !- broadcast of the mass distribution of volcanic ash
          call broadcast(mass_bin_dist, master_num, "mass_bin_dist")
          write(c2,"(i8)") mchnum      
          call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%src_dummy_2d &
               ,src_dummy_2d &
               ,trim(spc_aer_name(imode,ispc))//"_" &
               //trim(src_name(isrc)), oneControlVars)

          if( aerosol > 0) then
             if (spc_aer_alloc(src,imode,ispc) == on) then
                where( src_dummy_2d(:,:) < 0.) src_dummy_2d(:,:)  = 0.

                !-19OCT2020 srf - check this later for MATRIX aer mechanism 
                !           aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3) = src_dummy_2d(1:m2,1:m3)
                aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3) = emiss_ajust_aer(imode,ispc)*src_dummy_2d(1:m2,1:m3)
                !
                if (mchnum==master_num) write(*,fmt='("    Maxval,MinVal: ",2(E13.6,1X))') &
                     maxval(aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3)) &
                     ,minval(aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3))
                !################### LFR
                !            write(33,rec=recn) aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3)
                !            recn=recn+1
                !################### LFR


             endif

          else
             call fatal_error('aer memory not allocated')
          endif
          !%%%%%%%%%%     plume section for frp methodology %%%%%%%%%%%%%
       elseif(trim(section) == 'plumefre') then
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if(plumerise == 1.and. mchnum==master_num) then
             print*,"emission file prepared for plumerise 2 but plumerise flag in ramsin is not 2"
             call fatal_error('stop at routine read_sources 2')
          endif
          if (mchnum==master_num) then 
             print*,"   reading fre related products for plumerise model version 2"
             read(iunit,*)   read_mean      
          endif
          call broadcast(read_mean, master_num, "read_mean")
          if(trim(read_mean) == 'mean_fct' ) then     
             call readstorefullfieldandownchunk(ng,iunit,oneGlobalEmissData(ng)%flam_frac &
                  ,src_dummy_2d,trim(read_mean), oneControlVars)
             if(plumerise == 2) then
                where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
                plume_fre_g(iflam_frac)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
                if (mchnum==master_num) then 
                   print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(iflam_frac)%pvar) &
                        ,minval(plume_fre_g(iflam_frac)%pvar)
                endif
             endif
          elseif(trim(read_mean) == 'mean_frp' ) then     
             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%mean_frp &
                  ,src_dummy_2d,trim(read_mean), oneControlVars)
             if(plumerise == 2) then
                where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
                plume_fre_g(imean_frp)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)*1.e+6 ! convert from mw to w
                if (mchnum==master_num) then 
                   print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(imean_frp)%pvar) &
                        ,minval(plume_fre_g(imean_frp)%pvar)
                endif
             endif
          elseif(trim(read_mean) == 'std_frp' ) then
             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%std_frp &
                  ,src_dummy_2d,trim(read_mean), oneControlVars)
             if(plumerise == 2) then
                where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
                plume_fre_g(istd_frp)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)*1.e+6! convert from mw to w 
                if (mchnum==master_num) then 
                   print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(istd_frp)%pvar) &
                        ,minval(plume_fre_g(istd_frp)%pvar)
                endif
             endif
          elseif(trim(read_mean) == 'mean_size' ) then
             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%mean_size &
                  ,src_dummy_2d,trim(read_mean), oneControlVars)
             if(plumerise == 2) then
                where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
                plume_fre_g(imean_size)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)*1.e+6! convert from km2 to m2
                if (mchnum==master_num) then 
                   print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(imean_size)%pvar)& 
                        ,minval(plume_fre_g(imean_size)%pvar)
                endif
             endif
          elseif(trim(read_mean) == 'std_size' ) then
             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%std_size &
                  ,src_dummy_2d,trim(read_mean), oneControlVars)
             if(plumerise == 2) then
                where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
                plume_fre_g(istd_size)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)*1.e+6! convert from km2 to m2
                if (mchnum==master_num) then 
                   print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(istd_size)%pvar) &
                        ,minval(plume_fre_g(istd_size)%pvar)
                endif
             endif
          else
             call fatal_error('unknow error in frp methodology')
          endif
          !%%%%%%%%%%     plume section   %%%%%%%%%%%%%
       elseif(trim(section) == 'plume') then
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if(plumerise == 2 .and. mchnum==master_num) then
             print*,"emission file prepared for plumerise 1 but plumerise flag in ramsin is not 1"
             call fatal_error('stop at routine read_sources 1')
          endif
          if (mchnum==master_num) then
             read(iunit,*)   read_mean &
                  , dummy &
                  , read_veg_name       &
                  , read_ident_veg      &
                  , read_units
             iveg_agreg = read_ident_veg       
          endif
          call broadcast(iveg_agreg, master_num, "iveg_agreg")
          ! allocate broadcast area
          sizecharvec = len(read_mean)+len(read_veg_name)
          allocate(charvec(sizecharvec), stat=ierr2)
          if (ierr2/=0) then
             write(c0,"(i8)") ierr2
             write(c1,"(i8)") sizecharvec
             call fatal_error(h//" allocate charvec("//trim(adjustl(c1))//&
                  ") fails with stat="//trim(adjustl(c0)))
          endif
          ! master process prepares broadcast data
          if (mchnum==master_num) then
             lastchar = 0
             do nc=1,len(read_mean)
                charvec(lastchar+nc) = read_mean(nc:nc)
             enddo
             lastchar = lastchar+len(read_mean)
             do nc=1,len(read_veg_name)
                charvec(lastchar+nc) = read_veg_name(nc:nc)
             enddo
          endif
          ! broadcast character data to remaining processes
          ! envia read_mean
          call broadcast(charvec, master_num, "charvec")
          ! scatter broadcasted data
          lastchar=0
          do nc=1,len(read_mean)
             read_mean(nc:nc) = charvec(lastchar+nc)
          end do
          lastchar = lastchar+len(read_mean)
          do nc=1,len(read_veg_name)
             read_veg_name(nc:nc) = charvec(lastchar+nc)
          end do
          ! deallocate broadcast area
          deallocate(charvec, stat=ierr2)
          if (ierr2/=0) then
             write(c0,"(i8)") ierr2
             call fatal_error(h//" deallocate charvec fails with stat="//&
                  trim(adjustl(c0)))
          end if
          if(trim(read_mean) == 'mean_fct' ) then
             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%flam_frac, &
                  src_dummy_2d,                              &
                  trim(read_mean)//"_"//trim(read_veg_name), oneControlVars)

             if(plumerise == 1) then
                where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
                plume_mean_g(iveg_agreg)%flam_frac(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
             endif

          elseif(trim(read_mean) == 'firesize' ) then

             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%fire_size, &
                  src_dummy_2d,                              &
                  trim(read_mean)//"_"//trim(read_veg_name), oneControlVars)
             if(plumerise == 1) then
                where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
                plume_mean_g(iveg_agreg)%fire_size(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
             endif

          else
             call fatal_error('model not yet prepared for flam-frac of each one of species')
          endif


          !- volcanoes section  --------------------------------
       elseif(   trim(section) == 'volcanic-eruption' .or. &
            trim(section) == 'volcanic-degassing'     ) then

          if (mchnum==master_num) then
             read(iunit,*)  read_mean , &
                  read_units
             if(trim(section) == 'volcanic-eruption' ) type_volc_process="eruption"
             if(trim(section) == 'volcanic-degassing') type_volc_process="degassing"
          endif

          call broadcast(read_mean, master_num, "read_mean")
          call broadcast(read_units, master_num, "read_units")
          call broadcast(type_volc_process, master_num, "type_volc_process")

          !write(c2,"(i8)") mchnum
          !print *,'proc ',trim(adjustl(c2)),' volc: reading ',trim(read_mean),' units= ',trim(read_units)

          if(trim(read_mean) == 'inject_height' ) then

             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%plum_heigth, &
                  src_dummy_2d,                                &
                  trim(read_mean)//"_"//trim(read_units), oneControlVars)
             !             call vfirec(iunit,volc_mean_g%plum_heigth(:,:),m2*m3,'lin')
             if(volcanoes == on) then

                where( src_dummy_2d(:,:) < 0.) src_dummy_2d(:,:)  = 0.
                volc_mean_g%plum_heigth(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)

                if (mchnum==master_num)print*,'max inject_height (m) =',maxval(oneglobalemissdata(ng)%plum_heigth)
             endif
          elseif(trim(read_mean) == 'vent_elevation' ) then   

             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%vent_elev, &
                  src_dummy_2d,                              &
                  trim(read_mean)//"_"//trim(read_units), oneControlVars)

             if(volcanoes == on) then                
                where( src_dummy_2d(:,:) < 0.) src_dummy_2d(:,:)  = 0.    
                volc_mean_g%vent_elev(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)

                if (mchnum==master_num)print*,'max vent_elevation (m) =',maxval(oneglobalemissdata(ng)%vent_elev(:,:))

             endif
          elseif(trim(read_mean) == 'time_duration' ) then   

             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%duration, &
                  src_dummy_2d,                             &
                  trim(read_mean)//"_"//trim(read_units), oneControlVars)


             if(volcanoes == on) then                
                where( src_dummy_2d(:,:) < 0.) src_dummy_2d(:,:)  = 0.    

                volc_mean_g%duration(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
                if (mchnum==master_num)print*,'max time duration (m) =',maxval(oneglobalemissdata(ng)%duration(:,:))
             endif
          else
             call fatal_error('volc sources not prepared yet for: '//trim(read_mean))
          endif

       else 
          write (*,fmt='(a)') 'read_mean="'//trim(read_mean)//'"'
          call fatal_error('unknown error!')
       endif

    enddo

100 continue   

    if (mchnum==master_num) then 
       print*,'---------------------------------------------------------------------------'
       print*,'---------------------------------------------------------------------------'
       close (iunit) 
    endif

    !- volcanic ash mass distribution:
    if(AEROSOL == 1 .and. trim(type_volc_process)=="eruption") then
       !- inverting the array mass_bin_dist:
       !- bin 1 will be the smallest in size, bin 10 will be the bigger
       !
       !!do i=1,int(nbins/2)
       !!  rdummy                          =mass_bin_dist(nbins-(i-1))
       !!  mass_bin_dist(nbins-(i-1))=mass_bin_dist(i          )
       !!  mass_bin_dist(i          )=rdummy
       !!enddo 
       !
       ! 2nd way
       mass_bin_dist(:)=mass_bin_dist(nmodes:1:-1)
       print*, ' volc ash total mass = ', maxval(aer1_g(1,v_ash)%sc_src(1,:,:)),&
            mass_bin_dist(:); call flush(6)

       !
       !
       ! -- assuming the src array of ash bin = 1 (V_ASH1) contains the _total_ ash mass
       do imode=2,nmodes
          aer1_g(imode,v_ash)%sc_src(1,:,:)= mass_bin_dist(imode)*aer1_g(1,v_ash)%sc_src(1,:,:)

       enddo
       ! - corrected ash mass for bin 1
       aer1_g(1,v_ash)%sc_src(1,:,:)= mass_bin_dist(1)*aer1_g(1,v_ash)%sc_src(1,:,:)

       !print*," volc ash=",v_ash&
       !       ,1,mass_bin_dist(1),maxval(aer1_g(1,v_ash)%sc_src(1,:,:))

       ! -alterando a distribuicao de concentracao de cinzas 29/12/2013
       ! -aplicado ao caso do vulcano chileno PYHUE
       !       do imode=1,nmodes
       !         if(imode .le. 5) then
       !   aer1_g(imode,v_ash)%sc_src(1,:,:)=aer1_g(imode,v_ash)%sc_src(1,:,:)/60.
       !  Else
       !   aer1_g(imode,v_ash)%sc_src(1,:,:)=aer1_g(imode,v_ash)%sc_src(1,:,:)*1.65
       !         ENDIF         
       !       enddo
       !
    endif


    if(trim(type_volc_process)=="degassing" .and. volcanoes == on) then
       ! converting the plume_heigth of degass volcanic from meters above sea level
       ! to meters above the vent, to allow an unified treatment with the eruption
       ! volcanoes at the vertical mass distribution routine:
       volc_mean_g%plum_heigth(:,:) = volc_mean_g%plum_heigth(:,:) &
            - volc_mean_g%vent_elev  (:,:)
    endif

    !--(DMK-BRAMS-5.0-INI)-----------------------------------------------------------------------------------
    !    IF(AEROSOL > 0) THEN
    !       !- special section for aerosol sulfate 
    !       !- mode Aitken = 50% mode Accumulation
    !       !* ref: Stier et al., The aerosol-climate model ECHAM5-HAM. Atmos.
    !       !  Chem. Phys., 5,1125-1156,2005.
    !       IF(spc_aer_alloc(src,nucle,urban) == on) &
    !            aer1_g(nucle,urban)%sc_src(1,:,:)=0.5*aer1_g(accum,urban)%sc_src(1,:,:)
    !
    !       IF(spc_aer_alloc(src,accum,urban) == on) &
    !            aer1_g(accum,urban)%sc_src(1,:,:)=0.5*aer1_g(accum,urban)%sc_src(1,:,:)
    !
    !    ENDIF
    !--(DMK-BRAMS-5.0-INI)-----------------------------------------------------------------------------------

    deallocate (src_dummy_2d)


    !################### LFR
    !
    !  close(33)
    !
    !  open(unit=33,file='aer.ctl' &
    !       ,action='WRITE',status='replace',form='FORMATTED')
    !
    !  !writing the name of grads file
    !  write(33,*) 'dset ^aer.gra'
    !  !writing others infos to ctl
    !  write(33,*) 'undef -0.9990000E+34'
    !  write(33,*) 'title AerTest'
    !  write(33,*) 'xdef ',m2,' linear ',grid_g(1)%glon(1,1),grid_g(1)%glon(2,1)-grid_g(1)%glon(1,1)
    !  write(33,*) 'ydef ',m3,' linear ',grid_g(1)%glat(1,1),grid_g(1)%glat(1,2)-grid_g(1)%glat(1,1)
    !  write(33,*) 'zdef ',1,'levels',1000
    !  write(33,*) 'tdef 1 linear 00:00z01jan2018     1mo'
    !  write(33,*) 'vars ',4
    !  write(33,*) 'urban2 1 99 urban2'
    !  write(33,*) 'urban3 1 99 urban3' 
    !  write(33,*) 'bburn2 1 99 bburn2'
    !  write(33,*) 'bburn3 1 99 bburn3'
    !  write(33,*) 'endvars'!
    !
    !  close(33)
    !
    !################### LFR   


  end subroutine read_sources_vfm

  subroutine convert_to_tracer_density(m1,m2,m3,ia,iz,ja,jz,isrctime,      &
       nzpmax,dzt,rtgt,nsrc,chem1_src_g,   &
       bburn,chem_nspecies,spc_chem_alloc, &
       src,off,transport,aer1_g,aerosol,   &
       aer_nspecies,spc_aer_alloc,nmodes,  &
       aer_bburn,geoge,volcanoes,bioge,CO2,&
       ISFCL,spc_aer_name) 

    ! original
    integer              , intent(IN)    :: bioge
    integer , intent(IN) :: m1
    integer , intent(IN) :: m2
    integer , intent(IN) :: m3
    integer , intent(IN) :: ia
    integer , intent(IN) :: iz
    integer , intent(IN) :: ja
    integer , intent(IN) :: jz
    integer , intent(IN) :: isrctime
    integer , intent(IN) :: isfcl,CO2  !DSM

    ! grid_dims
    integer , intent(IN) :: nzpmax

    ! mem_grid
    real    , intent(IN) :: dzt(nzpmax)
    real    , intent(IN) :: rtgt(m2,m3)

    ! chem1_list
    integer , intent(IN) :: chem_nspecies
    integer , intent(IN) :: spc_chem_alloc(6,chem_nspecies)
    integer , intent(IN) :: src
    integer , intent(IN) :: off
    integer , intent(IN) :: transport

    ! aer1_list
    integer , intent(IN) :: aer_nspecies
    integer , intent(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    integer , intent(IN) :: nmodes
    integer , intent(IN) :: aer_bburn
    character(LEN=8) , intent(IN) :: spc_aer_name(nmodes,aer_nspecies)

    ! mem_chem1
    integer              , intent(IN)    :: nsrc
    type(chem1_src_vars) , intent(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    integer              , intent(IN)    :: bburn,geoge

    ! mem_aer1
    type (aer1_vars) , intent(INOUT) :: aer1_g(nmodes,aer_nspecies)
    integer          , intent(IN)    :: aerosol

    ! mem_volc_chem1
    integer          , intent(IN)    :: volcanoes    


    ! local var
    !Fator de conversao de unidades    
    real,parameter :: fcu =1.e+9 !=> ppbm 
    !LPCE
    real :: dz_inv,dz
    integer :: i,j,ksrc,isrc,ispc,imode,itim

    ksrc=2 !surface level of emission in the model

    !- chemistry section

    do j=ja,jz
       do i=ia,iz

          ! Todas as unidades estao em kg/m2/dia => use 'dz' em vez de 'vol'
          !    vol = 1./(dxt(i,j)*dyt(i,j)*dzt(k))*rtgt(i,j)

          !         dz        = grid_g(ng)%rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))
          dz        = rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))
          dz_inv    = 1./dz
          do ispc=1,chem_nspecies

             if(spc_chem_alloc(src,ispc) == off) cycle

             !- convert from kg/m^2  to  density (kg[gas]/m^3*1.e9)
             do isrc=1,nsrc

                !- only call geoge emissions if volcanoes is ON.
                if(isrc==geoge .and. volcanoes == off) cycle

                if(ISFCL == 5 .and. ispc == CO2 .and. isrc == bioge) cycle  !DSM e Saulo

                itim = actual_time_index(isrctime,isrc)
                chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) = &
                     chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) * fcu * dz_inv
             enddo
             ! copy smoldering emission from  bburn to ksrc
             itim = actual_time_index(isrctime,bburn)

             chem1_src_g(itim,bburn,ispc)%sc_src(ksrc,i,j) = &
                  chem1_src_g(itim,bburn,ispc)%sc_src(1   ,i,j)

          enddo
       enddo
    enddo

    !- aerosol section 
    if(AEROSOL > 0 ) then
       do j=ja,jz
          do i=ia,iz

             ! Todas as unidades estao em kg/m2/dia => use 'dz' em vez de 'vol'
             !    vol = 1./(dxt(i,j)*dyt(i,j)*dzt(k))*rtgt(i,j)

             !            dz        = grid_g(ng)%rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))
             dz        = rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))
             dz_inv    = 1./dz
             do ispc=1,aer_nspecies
                do imode=1,nmodes

                   if(spc_aer_alloc(src      ,imode,ispc) == off .or. &
                        spc_aer_alloc(transport,imode,ispc) == off) cycle

                   !- convert from kg/m^2  to  density (kg[aer]/m^3*1.e9)

                   aer1_g(imode,ispc)%sc_src(1,i,j) = &
                        aer1_g(imode,ispc)%sc_src(1,i,j) * fcu * dz_inv



                enddo
             enddo


             ! copy smoldering emission from aerosol bburn to ksrc
             if(AEROSOL==1) then

                do imode=1,nmodes

                   if(spc_aer_alloc(src      ,imode,aer_bburn) == off .or. &
                        spc_aer_alloc(transport,imode,aer_bburn) == off) cycle

                   aer1_g(imode,aer_bburn)%sc_src(ksrc,i,j) = aer1_g(imode,aer_bburn)%sc_src(1,i,j)
                enddo
             elseif(AEROSOL==2) then ! For MATRIX
                do ispc=1,aer_nspecies
                   do imode=1,nmodes

                      if(spc_aer_alloc(src      ,imode,ispc) == off .or. &
                           spc_aer_alloc(transport,imode,ispc) == off) cycle

                      if(spc_aer_name(imode,ispc)=="boc_bcar" .or. &
                           spc_aer_name(imode,ispc)=="boc_ocar" )    then
                         aer1_g(imode,ispc)%sc_src(ksrc,i,j) = aer1_g(imode,ispc)%sc_src(1,i,j)
                      endif
                   enddo
                enddo
             endif

          enddo
       enddo
    endif

  end subroutine convert_to_tracer_density

  !----------------------------------------------------------------------

  subroutine  emis_flam_smold(n1,n2,n3,isrctime,                 &
       nsrc,chem1_src_g,bburn,chem_nspecies, &
       spc_chem_alloc,src,on,off,transport,  &
       aer1_g,aerosol,aer_nspecies,          &
       spc_aer_alloc,nmodes,aer_bburn,       &
       plume_mean_g,plume_g,tropical_forest, &
       boreal_forest,savannah,grassland,     &
       nveg_agreg,spc_aer_name,plume_fre_g,plumerise)

    ! original
    integer , intent(IN) :: n1
    integer , intent(IN) :: n2
    integer , intent(IN) :: n3
    integer , intent(IN) :: isrctime,plumerise

    ! chem1_list
    integer          , intent(IN) :: chem_nspecies
    integer          , intent(IN) :: spc_chem_alloc(6,chem_nspecies)
    integer          , intent(IN) :: src
    integer          , intent(IN) :: on
    integer          , intent(IN) :: off
    integer          , intent(IN) :: transport

    ! aer1_list
    integer          , intent(IN) :: nmodes
    integer          , intent(IN) :: aer_nspecies
    integer          , intent(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    integer          , intent(IN) :: aer_bburn

    ! mem_chem1
    integer              , intent(IN)    :: nsrc
    type(chem1_src_vars) , intent(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    integer              , intent(IN)    :: bburn


    ! mem_aer1
    type (aer1_vars) , intent(INOUT) :: aer1_g(nmodes,aer_nspecies)
    integer          , intent(IN)    :: aerosol

    ! mem_plume_chem1
    integer               , intent(IN)    :: nveg_agreg
    type (plume_mean_vars), intent(INOUT) :: plume_mean_g(nveg_agreg)
    type (plume_fre_vars) , intent(INOUT) :: plume_fre_g(5)
    type (plume_vars)     , intent(INOUT) :: plume_g(nveg_agreg,chem_nspecies)
    integer               , intent(IN)    :: tropical_forest
    integer               , intent(IN)    :: boreal_forest
    integer               , intent(IN)    :: savannah
    integer               , intent(IN)    :: grassland
    character(LEN=8) , intent(IN) :: spc_aer_name(nmodes,aer_nspecies)

    real,dimension(n2,n3) :: smold_frac 
    integer iv,ispc,i,j,imode,itim
    integer:: imean_plume
    imean_plume = 1 !change this at alloc_plume_chem1 routine also


    !- time index of memory allocation position        
    itim = actual_time_index(isrctime,bburn)
    if(itim > 1) then
       call fatal_error('Time level 2 not allowed when plumerise is used!')
    endif

    if(imean_plume == on) then
       !-----  
       !- calcula a emissao smoldering e fatores para obtencao da fracao
       !- flaming em funcao da emissao smoldering
       if(plumerise == 1) then
          smold_frac(1:n2,1:n3) = 1.- ( plume_mean_g(tropical_forest)%flam_frac(1:n2,1:n3) + &
               plume_mean_g(boreal_forest)%flam_frac(1:n2,1:n3) + &
               plume_mean_g(savannah     )%flam_frac(1:n2,1:n3) + &
               plume_mean_g(grassland    )%flam_frac(1:n2,1:n3))    
       elseif(plumerise == 2) then

          smold_frac(1:n2,1:n3) = 1.- plume_fre_g(iflam_frac)%pvar(1:n2,1:n3)

       endif

       !- chemistry section (only for bburn source)
       do ispc = 1,chem_nspecies

          if(spc_chem_alloc(src,ispc) /= on) cycle 

          !- convert from 'total' emisson to 'smoldering' part
          chem1_src_g(itim,bburn,ispc)%sc_src(1,:,:) = smold_frac(:,:) * &  
               chem1_src_g(itim,bburn,ispc)%sc_src(1,:,:)
       enddo
       if(AEROSOL == 1 ) then
          !- aerosol section (only for bburn aerosols)
          do imode=1,nmodes

             if(spc_aer_alloc(src,imode,aer_bburn)       == off .or. &
                  spc_aer_alloc(transport,imode,aer_bburn) == off) cycle

             !- convert from 'total' emisson to 'smoldering' part
             aer1_g(imode,aer_bburn)%sc_src(1,:,:) = smold_frac(:,:) * & 
                  aer1_g(imode,aer_bburn)%sc_src(1,:,:)

          enddo
          !- this for MATRIX
       elseif(AEROSOL ==  2 ) then
          !- aerosol section (only for bburn aerosols)
          do ispc = 1,aer_nspecies

             do imode=1,nmodes

                if(spc_aer_alloc(src,imode,ispc)       == off .or. &
                     spc_aer_alloc(transport,imode,ispc) == off) cycle

                !-only for bburn aerosols) 
                if(spc_aer_name(imode,ispc)=="boc_bcar" .or. &
                     spc_aer_name(imode,ispc)=="boc_ocar" )    then
                   !print*,"bburn=",spc_aer_name(imode,ispc),maxval(aer1_g(imode,ispc)%sc_src(1,:,:))

                   !- convert from 'total' emisson to 'smoldering'
                   aer1_g(imode,ispc)%sc_src(1,:,:) = smold_frac(:,:) * & 
                        aer1_g(imode,ispc)%sc_src(1,:,:)
                endif
             enddo
          enddo
       endif

       !- convert from flaming fraction to relationship with phase smoldering emission
       if(plumerise == 1) then
          do iv = 1, nveg_agreg
             plume_mean_g(iv)%flam_frac(:,:) = plume_mean_g(iv)%flam_frac(:,:)/ &
                  (1.e-8+smold_frac(:,:))
          enddo
       elseif(plumerise == 2) then

          plume_fre_g(iflam_frac)%pvar(:,:) = plume_fre_g(iflam_frac)%pvar(:,:)/ &
               (1.e-8+smold_frac(:,:))       
       endif
       !-----

    else

       !-----      case where each specie has his own flaming fraction ----------------
       call fatal_error('Aerosol emission not ready for this option!')

       do ispc = 1,chem_nspecies
          if(spc_chem_alloc(src,ispc) /= on) cycle 
          smold_frac(:,:) = 1.- ( plume_g(tropical_forest,ispc)%fct(:,:) + &
               plume_g(boreal_forest,ispc)%fct(:,:) + &
               plume_g(savannah    ,ispc)%fct(:,:) + &
               plume_g(grassland    ,ispc)%fct(:,:)   )

          !- convert from 'total' emisson to 'smoldering'
          chem1_src_g(itim,bburn,ispc)%sc_src(1,:,:) = smold_frac(:,:) * &  
               chem1_src_g(itim,bburn,ispc)%sc_src(1,:,:)

          !- convert from flaming fraction to relationship with phase smoldering emission
          do iv = 1, nveg_agreg
             plume_g(iv,ispc)%fct(:,:) = plume_g(iv,ispc)%fct(:,:)/ &
                  (1.e-8+smold_frac(:,:))
             !- flamming emission =  plume_g(iv,iscp,ng)%fct(:,:) * &
             !             chem1_src_g(itim,bburn,ispc,ng)%sc_src(1,:,:)
          enddo

       enddo
       !-----
    endif
  end subroutine  emis_flam_smold
  !----------------------------------------------------------------------


  !----------------------------------------------------------------------

  subroutine sources(m1,m2,m3,ia,iz,ja,jz,itime1,time,imonth1,idate1,iyear1,glon,    &
       chem_nspecies,spc_chem_alloc,src,on,off,transport,nsrc,nvert_src,  &
       chem1_src_g,chem1_g,bburn,bioge,antro,geoge,diur_cycle,aer_nvert_src,    &
       aer1_g,aerosol,aer_nspecies,spc_aer_alloc,nmodes,aer_bburn,        &
       aer_sdust,aer_urban,aer_bioge,aer_marin,dnp,iexev,dn0,cosz,spc_chem_name, &
       emiss_cycle,volcanoes,aer_v_ash,CO2,isfcl,spc_aer_name,matrix_level,&
       aer2_g,SO2)

    use aer1_list, only :  akk,sulf,acc,bc1,bcar,occ,ocar,dd1,dust,dd2,boc, numb_alloc  

    implicit none

    ! original
    integer , intent(IN) :: m1
    integer , intent(IN) :: m2
    integer , intent(IN) :: m3
    integer , intent(IN) :: ia
    integer , intent(IN) :: iz
    integer , intent(IN) :: ja
    integer , intent(IN) :: jz
    real    , intent(IN) :: time
    integer , intent(IN) :: imonth1
    integer , intent(IN) :: idate1
    integer , intent(IN) :: iyear1
    integer , intent(IN) :: itime1
    integer , intent(IN) :: isfcl,CO2  !DSM
    integer , intent(IN) :: SO2  ! for matrix

    ! chem1_list
    integer , intent(IN) :: chem_nspecies
    integer , intent(IN) :: spc_chem_alloc(6,chem_nspecies)
    integer , intent(IN) :: src
    integer , intent(IN) :: on
    integer , intent(IN) :: off
    integer , intent(IN) :: transport
    character(LEN=8), intent(IN),dimension(chem_nspecies) :: spc_chem_name

    ! aer1_list
    integer , intent(IN) :: aer_nspecies
    integer , intent(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    integer , intent(IN) :: nmodes
    integer , intent(IN) :: aer_bburn
    integer , intent(IN) :: aer_sdust
    integer , intent(IN) :: aer_urban
    integer , intent(IN) :: aer_bioge
    integer , intent(IN) :: aer_marin
    integer , intent(IN) :: aer_v_ash
    character(LEN=8) , intent(IN) :: spc_aer_name(nmodes,aer_nspecies)
    character(LEN=1 ), intent(IN) :: matrix_level

    ! mem_chem1
    integer              , intent(IN)    :: nsrc
    integer              , intent(IN)    :: nvert_src(nsrc)
    type(chem1_src_vars) , intent(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    type(chem1_vars)     , intent(INOUT) :: chem1_g(chem_nspecies)
    integer              , intent(IN)    :: bburn
    integer              , intent(IN)    :: bioge
    integer              , intent(IN)    :: geoge
    integer              , intent(IN)    :: antro
    integer              , intent(IN)    :: diur_cycle(nsrc)

    ! mem_aer1
    integer         , intent(IN)    :: aer_nvert_src(aer_nspecies)
    type(aer1_vars) , intent(INOUT) :: aer1_g(nmodes,aer_nspecies)
    type(aer1_vars) , intent(INOUT) :: aer2_g(nmodes)
    integer         , intent(IN)    :: aerosol

    ! mem_stilt
    real    , pointer    :: dnp(:,:,:) ! in
    integer , intent(IN) :: iexev

    ! mem_basic
    real , pointer :: dn0(:,:,:) ! in

    ! mem_radiate
    real , intent(IN) :: cosz(m2,m3)

    type(cycle_emission), intent(INOUT) :: emiss_cycle(nsrc)

    ! mem_volc_chem1
    integer, intent(IN) :: volcanoes

    integer :: k_src,k_tend,it1,it2

    double precision :: tlinterp
    !-MFA
    real :: timeq2,timeq3,gglon,fuso,alfa(nsrc)
    !-for Lagrange
    real :: g,no,src1,src2,srcn1,srcn2,time1,time2,ztime1,ztime2,local_hour,htime1,htime2
    double precision :: tlinterp2
    !- ANTRO FROM CETESB-
    integer, parameter :: diur_cetesb_flag=0
    !-MFA

    integer :: iweek,idays,j,i,k,ispc,isrc,k2,imode
    real :: tign,strtim,timeq,r_q,r_antro,real_time,jd
    real, dimension(7) :: week_CYCLE
    !                     dia da semana:    SEG   TERQUA   QUI   SEX   SAB  DOM  
    !                            iweek=     1      2  34     5     6 7
    !- dados cetesb/campinas/2005
    data (week_CYCLE(iweek),iweek=1,7) /1.1, 1.1, 1.1, 1.1, 1.1, 0.83,0.67/ !total = 7

    real rt(nsrc),rt_aer(aer_nspecies)
    real, parameter :: bx_bburn  = 18.041288 * 3600., & !- pico em 18 UTC
         cx        = 2.5 * 3600.,       & ! 2.184936 * 3600., &
         rinti     = 0.8941851* 2.1813936e-8    , & ! 1/integral
         ax        = 2000.6038        , &
         bx_antro  = 9. *3600.  ,&  ! local time of peak 1 
         cx_antro  = 16.*3600.  ,&  ! local time of peak 2
         rsum      = 1.3848466E-05  !2.1311416E-05   ! 1/integral

    real, pointer, save,dimension(:,:,:) :: rho_air
    real, dimension(m2,m3) :: glon

    real local_cosz(m2,m3), local_emiss_bioge_diur_cycle(m2,m3)

    !- nocturnal/background/constant emission for biogenic/urban-industrial-transp processes
    real, parameter :: f_nct=0.15                  &! 15% per day
         , f_nct_dvd86400=f_nct/86400. &
         , um86400=1./86400.

    !- parameters for converting emission in mass to emission in number 
    !- only for matrix 
    integer, parameter :: nemis_spcs        = 10 
    real(4), parameter :: pi                = 3.141592654
    real(4), parameter :: pi6               = pi/6.0
    real(4), parameter :: emis_dens_sulf    = 1.770e+00 !g/cm^3]
    real(4), parameter :: emis_dens_bcar_bb = 1.390e+00 ![g/cm^3]
    real(4), parameter :: emis_dens_bcar_ur = 1.700e+00 ![g/cm^3]
    real(4), parameter :: emis_dens_ocar    = 1.000e+00 ![g/cm^3]
    real(4), parameter :: emis_dens_seas    = 2.165e+00 ![g/cm^3]
    real(4), parameter :: emis_dens_dust    = 2.600e+00 ![g/cm^3]
    real, dimension(nemis_spcs) :: emis_dens = (/ emis_dens_sulf, emis_dens_sulf, &
         emis_dens_bcar_ur, emis_dens_bcar_bb,emis_dens_bcar_bb, emis_dens_ocar,  &
         emis_dens_seas   , emis_dens_seas   , emis_dens_dust  , emis_dens_dust /)



    real(4), parameter :: dg_akk_sulf = 0.013   ![um]
    real(4), parameter :: dg_acc_sulf = 0.068   ![um]
    real(4), parameter :: dg_bc1_bcar = 0.030   ![um]
    real(4), parameter :: dg_boc_bcar = 0.021   ![um]
    real(4), parameter :: dg_boc_ocar = 0.021   ![um]
    real(4), parameter :: dg_occ_ocar = 0.030   ![um]
    real(4), parameter :: dg_ssa_seas = 0.370   ![um]
    real(4), parameter :: dg_ssc_seas = 3.930   ![um]
    real(4), parameter :: dg_dd1_dust = 0.580   ![um]
    real(4), parameter :: dg_dd2_dust = 5.400   ![um]


    real, dimension(nemis_spcs) :: dgn0_emis = (/ dg_akk_sulf, dg_acc_sulf, &
         dg_bc1_bcar, dg_boc_bcar, dg_boc_ocar, dg_occ_ocar, dg_ssa_seas, &
         dg_ssc_seas, dg_dd1_dust, dg_dd2_dust /)

    real(4), parameter :: sg_akk_sulf = 1.600d+00
    real(4), parameter :: sg_acc_sulf = 1.800d+00
    real(4), parameter :: sg_bc1_bcar = 1.800d+00
    real(4), parameter :: sg_boc_bcar = 1.550d+00
    real(4), parameter :: sg_boc_ocar = 1.555d+00
    real(4), parameter :: sg_occ_ocar = 1.800d+00
    real(4), parameter :: sg_ssa_seas = 1.800d+00
    real(4), parameter :: sg_ssc_seas = 2.000d+00
    real(4), parameter :: sg_dd1_dust = 1.800d+00
    real(4), parameter :: sg_dd2_dust = 1.800d+00 

    real, dimension(nemis_spcs) :: sig0_emis = (/ sg_akk_sulf, sg_acc_sulf, &
         sg_bc1_bcar, sg_boc_bcar, sg_boc_ocar, sg_occ_ocar, sg_ssa_seas, &
         sg_ssc_seas, sg_dd1_dust, sg_dd2_dust /)

    real, dimension(nemis_spcs) ::recip_part_mass,dp0_emis
    real, dimension(m1,m2,m3) :: dummy_src
    real :: factor
    !-end of parameters for matrix.

    dummy_src=0.0

    !- if using mass conservation fix : air dens changes with the time evolution
    !   IF( iexev == 2 )  rho_air => stilt_g(ng)%dnp(:,:,:) 
    if( iexev == 2 )  rho_air => dnp(:,:,:) 

    !- if not, air dens = air dens of basic state and  need to define for each when 
    !- have nested grids
    !   IF( iexev == 1  ) rho_air => basic_g(ng)%dn0(:,:,:) 
    if( iexev == 1  ) rho_air => dn0(:,:,:) 
    !
    !number of days of simulation
    idays = int(( float(itime1)/100. + time/3600.)/24.+.00001)  



    if(diur_cycle(bburn) == ON) then
       !-------------biomass burning diurnal cycle --------------------
       tign  = real(idays)*24.*3600.

       ! Modulacao da queimada media durante o ciclo diurno(unidade: 1/s)
       ! com a int( r_q dt) (0 - 24h)= 1.
       timeq= ( time + float(itime1)*0.01*3600. - tign )

       r_q  = rinti*( ax * exp( -(timeq-bx_bburn)**2/(2.*cx**2) ) + 100. -  &
            5.6712963e-4*( timeq ))

       emiss_cycle(bburn)%emission_rate(:,:)=r_q
       rt(bburn)= r_q
       alfa(bburn) = 0.
    else
       emiss_cycle(bburn)%emission_rate(:,:)= um86400 ! = 1./86400.
       rt(bburn)= um86400 ! = 1./86400.
       alfa(bburn) = 1.
    endif


    if(diur_cycle(antro) == ON) then

       !------------- anthropogenic diurnal cycle (industrial,urban, ...)
       real_time = time + float(itime1)*0.01*3600. !UTC

       ! weekly cycle
       ! week day
       !v1
       iweek= int(((float(julday(imonth1,idate1,iyear1))/7. - &
            int(julday(imonth1,idate1,iyear1)/7))*7.)) + 3
       if(iweek.gt.7) iweek = iweek-7
       !v2
       !call Greg2Jul(0, idate1, imonth1, iyear1, jd)
       !jd=jd+float(nint(time/86400.))
       !iweek = int(AMOD(jd, 7.)) 
       !if(iweek < 1 .or. iweek > 7) stop 315

       !- diurnal cycle
       if (diur_cetesb_flag == 0) then
          do j=ja,jz
             do i=ia,iz

                gglon = glon(i,j)
                fuso = int(gglon/15)
                !-to better keep continuity
                idays = int(( real_time +fuso*3600. )/86400.+.00001)  
                tign  = real(idays)*86400.

                timeq2= (real_time  -tign + fuso*3600.) - bx_antro
                timeq3= (real_time  -tign + fuso*3600.) - cx_antro
                !r_antro  =(exp(-((timeq2)**2)/((3.*3600) **2))+(exp(-((timeq3)**2)/((3.*3600) **2))+0.1))*rsum
                r_antro=(exp(-(timeq2**2)/(18400.**2))+exp(-(timeq3**2)/(18500.**2))+0.1)*rsum

                !- weekly + diurnal cycle
                r_antro = r_antro * week_CYCLE(iweek)

                emiss_cycle(antro)%emission_rate(i,j)=r_antro 
                alfa(antro) = 0.
             enddo
          enddo

       endif
       if (diur_cetesb_flag == 1) then

          do j=ja,jz
             do i=ia,iz

                gglon = glon(i,j)
                fuso  = int(gglon/15)

                !-to better keep continuity

                idays = int(( real_time )/86400.+.00001)  
                tign  = real(idays)*86400.


                local_hour = ((real_time - tign)/3600)+fuso


                if (local_hour.lt.0) then
                   local_hour = (24) + local_hour
                endif


                call interpolation_antro(local_hour,src1,src2,srcn1,srcn2,time1,time2)


                htime1 = time1 - fuso
                if (htime1.ge.24) then
                   htime1 = htime1 - 24
                endif

                htime2 = time2 - fuso
                if (htime2.ge.24) then
                   htime2 = htime2 - 24
                endif


                ztime1 = ((htime1*3600) + tign)
                ztime2 = ((htime2*3600) + tign)


                tlinterp=dble(src1 + (((time - ztime1)*(src2 - src1))/(ztime2 - ztime1)))
                tlinterp2=dble(srcn1 + (((time - ztime1)*(srcn2 - srcn1))/(ztime2 - ztime1)))

                !PRINT*,'TESTE',time,local_hour,time1,time2,ztime1,ztime2,tign,src2,src1,tlinterp,fuso

                emiss_cycle(antro)%emission_rate    (:,:)=tlinterp/3600
                emiss_cycle(antro)%emission_rate_NOX(:,:)=tlinterp2/3600

                alfa(antro)= 0.    
                r_antro=tlinterp/3600.
             enddo
          enddo
       endif
       if (diur_cetesb_flag.ne.1.and.diur_cetesb_flag.ne.0) then
          call fatal_error('No definition for diur_cetesb_flag!')
       endif

    else
       !----------- sources linearly time interpolated
       if(srctime2<=srctime1) then
          call fatal_error('srctime2<=srctime1! linterp')
       endif
       tlinterp=dble(time-srctime1)/dble(srctime2-srctime1)
       emiss_cycle(antro)%emission_rate(:,:)=tlinterp
       alfa(antro)= 1.
    endif

    if(diur_cycle(bioge) == ON) then

       !---------- sources with diurnal cycle and spatial dependence
       ! - using zenital angle from radiate routine/including constant background emission (f_cnt%)
       !       emiss_cycle(bioge,ng)%emission_rate(:,:)= f_nct_dvd86400 + emiss_cycle(bioge,ng)%dcnorma_inv(:,:)&    
       !            *MAX(0.,radiate_g(ng)%cosz(:,:)) * (1.-f_nct) 
       emiss_cycle(bioge)%emission_rate(:,:)= f_nct_dvd86400 + emiss_cycle(bioge)%dcnorma_inv(:,:)&    
            *max(0.,cosz(:,:)) * (1.-f_nct) 
       alfa(bioge) = 0.
    else
       !----------- sources linearly time interpolated    
       if(srctime2<=srctime1) then
          call fatal_error('srctime2==srctime1! linterp')
       endif
       tlinterp=dble(time-srctime1)/dble(srctime2-srctime1)
       emiss_cycle(bioge)%emission_rate(:,:) = tlinterp
       alfa(bioge) = 1.
    endif

    if(diur_cycle(geoge) == ON  ) then
       emiss_cycle(geoge)%emission_rate(:,:)=um86400 ! 1./86400 sec
       alfa(geoge) = 0.
    else
       !----------- sources linearly time interpolated    
       if(srctime2<=srctime1) then
          call fatal_error('srctime2==srctime1! linterp')
       endif
       tlinterp=dble(time-srctime1)/dble(srctime2-srctime1)
       emiss_cycle(geoge)%emission_rate(:,:) = tlinterp
       alfa(geoge) = 1.
    endif

    !-srf temporary array to store emission cycle for biogenic species.
    !-    (CO2 from JULES is already instantaneous while others biogenic emissions 
    !-    from MEGAN is per day)
    if(ISFCL == 5) local_emiss_bioge_diur_cycle(:,:)=emiss_cycle(bioge)%emission_rate(:,:)
    !

    ! print*,'emiss_cycle(antro,ng)%emission_rate=',emiss_cycle(antro,ng)%emission_rate(ia,ja)&
    ! ,srctime1,srctime2,tlinterp

    !-------------------------- perform emissions
    !- chemistry section 
    do ispc=1,chem_nspecies

       if(spc_chem_alloc(src,ispc) == off .or. spc_chem_alloc(transport,ispc) == off)  cycle

       do isrc=1,nsrc

          !- DSM/SRF - special treatment for biogenic/jules (CO2) emissions
          if(ISFCL == 5 .and. isrc == bioge) then
             if( ispc == CO2 ) then
                emiss_cycle(bioge)%emission_rate(:,:)=1.
             else
                emiss_cycle(bioge)%emission_rate(:,:)=local_emiss_bioge_diur_cycle(:,:)
             end if
          end if

          !- only call geoge emissions if volcanoes is ON.
          if(isrc==geoge .and. volcanoes == off) cycle

          !- memory position of source array for each time level
          it1 = 1 ! always 1
          it2 = actual_time_index(2,isrc)!might be 1 or 2, will be 1 if linterp = OFF
          !=>chem1_src_g(it1,:,:,:)% = chem1_src_g(it2,:,:,:)%

          k_src = 1 
          !- for anthropogenic and biogenic emissions
          k2 = 2                                          ! control for 2-dim src emission field

          !- for biomass burning and volcanic emissions          
          if(isrc == bburn .or. isrc == geoge)  k2 = m1-1  ! control for 3-dim src emission field

          do k=2,k2

             if(isrc == bburn .or. isrc == geoge) k_src=k  ! control for 3-dim source emission field


             if (diur_cetesb_flag == 1 .and. isrc == antro) then

                call source_to_tend_cycle_cetesb ( m1,m2,m3,ia,iz,ja,jz   &
                     ,chem1_src_g(it1,isrc,ispc)%sc_src &! source data at time level 1
                     ,chem1_src_g(it2,isrc,ispc)%sc_src &! source data at time level 2
                     ,chem1_g (   ispc)%sc_t&! tendency array
                     ,nvert_src(isrc)        &! vertical size of source array
                     ,emiss_cycle(isrc)%emission_rate   &! diurnal cycle of emission
                     ,alfa(isrc)&! alfa cte 
                     ,k_src&! vertical level of source array (where data is stored)
                     ,k &! vertical level of tendency
                     ,rho_air&! air density (to convert to mixing ratio tendency)
                     ,emiss_cycle(isrc)%emission_rate_NOX &! peso hora atual NOX     
                     ,ispc,spc_chem_name,chem_nspecies  )      

             else

                call source_to_tend_cycle ( m1,m2,m3,ia,iz,ja,jz   &
                     ,chem1_src_g(it1,isrc,ispc)%sc_src   &! source data at time level 1
                     ,chem1_src_g(it2,isrc,ispc)%sc_src   &! source data at time level 2
                     ,chem1_g (   ispc)%sc_t      &! tendency array
                     ,nvert_src(isrc)   &! vertical size of source array
                     ,emiss_cycle(isrc)%emission_rate   &! diurnal cycle of emission
                     ,alfa(isrc)   &! alfa cte 
                     ,k_src   &! vertical level of source array (where data is stored)
                     ,k    &! vertical level of tendency
                     ,rho_air)      ! air density (to convert to mixing ratio tendency)
             endif

          enddo
       enddo
    enddo

    !- aerosol section ----------------------------------------
    if(AEROSOL == 1)  then
       !- still need implementation of the emission cycle with space dependence
       rt_aer(aer_bburn) = rt(bburn)
       rt_aer(aer_sdust) = 1.0 ! on-line emission
       rt_aer(aer_urban) = 1.157407e-5 !rt(antro)  ! < must be fixed later with actual diurnal cycle
       rt_aer(aer_bioge) = 1.157407e-5 
       rt_aer(aer_marin) = 1.0 ! on-line emission 
       rt_aer(aer_v_ash) = 1.157407e-5 ! should be 1./duration

       do ispc=1,aer_nspecies

          do imode=1,nmodes

             if(spc_aer_alloc(src      ,imode,ispc) == off .or. &
                  spc_aer_alloc(transport,imode,ispc) == off) cycle


             k_src = 1;  k2 = 2     ! control for 2-dim aerosol (bioge, antro, marin, sdust)
             if(ispc == aer_bburn .or. ispc == aer_v_ash)  k2 = m1-1  ! control for 3-dim bburn aerosol

             do k=2,k2

                if(ispc == aer_bburn .or. ispc == aer_v_ash) k_src=k     ! control  for 3-dim bburn aerosol

                call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                     ,aer1_g(imode,ispc)%sc_src           &! source data
                     ,aer1_g(imode,ispc)%sc_t  &! tendency array
                     ,aer_nvert_src(ispc)  &! vertical size of source array
                     ,rt_aer(ispc)  &! diurnal cycle of emission
                     ,k_src   &! vertical level of source array (where data is stored)
                     ,k  &! vertical level of tendency
                     ,rho_air)   ! air density (to convert to mixing ratio tendency)

             enddo
          enddo
       enddo
    elseif(AEROSOL == 2)  then

       if(matrix_level .ne. "1") stop "Emissions are only configured for matrix_level=1"

       !- section for emission of mass --------------------
       do ispc=1,aer_nspecies
          do imode=1,nmodes
             if(spc_aer_alloc(src      ,imode,ispc) == off .or. &
                  spc_aer_alloc(transport,imode,ispc) == off) cycle

             !- case anthropogenic emssions
             if(   spc_aer_name(imode,ispc)=="occ_ocar" .or. &
                  spc_aer_name(imode,ispc)=="bc1_bcar" ) then
                k2   = 2 ! control for 2-dim aerosol (bioge, antro, marin, sdust)
                rt_aer(ispc) = r_antro
                k_src = 1
                !print*,"anthro",maxval(aer1_g(imode,ispc)%sc_src),spc_aer_name(imode,ispc), rt_aer(ispc)

                !DO k=2,k2
                call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                     ,aer1_g(imode,ispc)%sc_src           &! source data
                     ,aer1_g(imode,ispc)%sc_t  &! tendency array
                     ,aer_nvert_src(ispc)          &! vertical size of source array
                     ,rt_aer(ispc)  &! diurnal cycle of emission
                     ,k_src           &! vertical level of source array (where data is stored)
                     ,k2  &! vertical level of tendency
                     ,rho_air)   ! air density (to convert to mixing ratio tendency)
                ! ENDDO

                !- case sea salt emssions
             elseif( spc_aer_name(imode,ispc)=="ssa_seas" .or. &
                  spc_aer_name(imode,ispc)=="ssc_seas" )    then
                k2   = 2 ! control for 2-dim aerosol (bioge, antro, marin, sdust)
                rt_aer(ispc) = 1.0
                k_src = 1
                !print*,"seas",maxval(aer1_g(imode,ispc)%sc_src),spc_aer_name(imode,ispc)

                !DO k=2,k2
                call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                     ,aer1_g(imode,ispc)%sc_src           &! source data
                     ,aer1_g(imode,ispc)%sc_t  &! tendency array
                     ,aer_nvert_src(ispc)          &! vertical size of source array
                     ,rt_aer(ispc)  &! diurnal cycle of emission
                     ,k_src           &! vertical level of source array (where data is stored)
                     ,k2  &! vertical level of tendency
                     ,rho_air)   ! air density (to convert to mixing ratio tendency)
                !ENDDO

                !- case biomass burning emissions !the same must be for volcanic ASH
             elseif(spc_aer_name(imode,ispc)=="boc_bcar" .or. &
                  spc_aer_name(imode,ispc)=="boc_ocar" )    then
                !print*,"bb",maxval(aer1_g(imode,ispc)%sc_src),spc_aer_name(imode,ispc), rt(bburn)
                rt_aer(ispc) = rt(bburn)
                k2 = m1-1  !
                do k=2,k2
                   k_src=k

                   call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                        ,aer1_g(imode,ispc)%sc_src           &! source data
                        ,aer1_g(imode,ispc)%sc_t  &! tendency array
                        ,aer_nvert_src(ispc)          &! vertical size of source array
                        ,rt_aer(ispc)  &! diurnal cycle of emission
                        ,k_src           &! vertical level of source array (where data is stored)
                        ,k  &! vertical level of tendency
                        ,rho_air)   ! air density (to convert to mixing ratio tendency)
                enddo
             endif
          enddo
       enddo
       !- section for emission of mass 
       !- special section for on-line emissions (not from prep-chem-src)
       !- this section will provide emissions for aerosol species akk_sulf and acc_sulf
       !- for antro and bio-burn processes. These emissions are based on the SO2 gas emission
       !===> be sure emissions for the SO2 gas is provided (arrays chem1_src_g(*,isrc,SO2)%sc_src)
       if(spc_chem_alloc(src,SO2) == ON)  then

          do ispc=sulf,sulf
             do imode=1,nmodes
                if(spc_aer_alloc(transport,imode,ispc) == off) cycle


                if(   spc_aer_name(imode,ispc)=="akk_sulf" .or. &
                     spc_aer_name(imode,ispc)=="acc_sulf"      ) then
                   !print*,"src emissions for aer=",imode, ispc,spc_aer_name(imode,ispc)

                   if(   spc_aer_name(imode,ispc)=="akk_sulf" ) factor = 0.01 * 0.025* 96.07/62.66  

                   if(   spc_aer_name(imode,ispc)=="acc_sulf" ) factor = 0.99 * 0.025* 96.07/62.66  


                   ! 1) case anthropogenic emissions
                   isrc  = antro
                   dummy_src(1:m1,1:m2,1:m3) = 0.5*( chem1_src_g(it1,isrc,SO2)%sc_src(1:m1,1:m2,1:m3) + &
                        chem1_src_g(it2,isrc,SO2)%sc_src(1:m1,1:m2,1:m3)  )

                   k2   = 2 ! control for 2-dim aerosol (bioge, antro, marin, sdust)
                   rt_aer(ispc) = r_antro * factor
                   k_src = 1
                   call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                        ,dummy_src                &! source data
                        ,aer1_g(imode,ispc)%sc_t  &! tendency array
                        ,aer_nvert_src(ispc)&! vertical size of source array
                        ,rt_aer(ispc)&! diurnal cycle of emission
                        ,k_src&! vertical level of source array (where data is stored)
                        ,k2&! vertical level of tendency
                        ,rho_air)  ! air density (to convert to mixing ratio tendency)

                   ! 2) case biomass burning emissions
                   isrc  = bburn
                   dummy_src(1:m1,1:m2,1:m3) = 0.5*( chem1_src_g(it1,isrc,SO2)%sc_src(1:m1,1:m2,1:m3) + &
                        chem1_src_g(it2,isrc,SO2)%sc_src(1:m1,1:m2,1:m3)  )
                   rt_aer(ispc) = rt(bburn)* factor
                   k2 = m1-1  !
                   do k=2,k2
                      k_src=k

                      call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                           ,dummy_src                         &! source data
                           ,aer1_g(imode,ispc)%sc_t  &! tendency array
                           ,aer_nvert_src(ispc)          &! vertical size of source array
                           ,rt_aer(ispc)  &! diurnal cycle of emission
                           ,k_src           &! vertical level of source array (where data is stored)
                           ,k  &! vertical level of tendency
                           ,rho_air)   ! air density (to convert to mixing ratio tendency)
                   enddo

                endif

             enddo
          enddo
       else
          stop "Emissions for mass of aer akk/acc are required, but SO2 emissions are not set"

       endif
       !
       !print*,"done emissions for mass of akk acc"; call flush(6)
       !
       !--- section for emission of number concentration
       !- A) getting parameter recip_part_mass
       do i =1,nemis_spcs

          dp0_emis(i) = 1.0e-06 * dgn0_emis(i) * exp( 1.5e+00 * ( log(sig0_emis(i)) )**2 )  ! convert from [um] to [m]

          !recip_part_mass(i) = 1.0e+00 / ( 1.0e+12 * emis_dens(i) * pi6 * dp0_emis(i)**3 )
          !- unit 1/kg
          recip_part_mass(i) = 1.0e+00 / ( 1.0e+3 * emis_dens(i) * pi6 * dp0_emis(i)**3 )

       enddo
       !
       !- B) getting number emissions in terms of mass emissions and performing the emissions
       !
       do i =1,nemis_spcs
          if(i==1 ) then
             imode=akk
             ispc =sulf
          endif  !   
          if(i==2 ) then
             imode=acc
             ispc =sulf
          endif
          if(i==3 ) then
             imode=bc1
             ispc =bcar
          endif
          if(i==4 ) then
             imode=boc
             ispc =bcar
          endif  ! boc bcar + boc ocar
          if(i==5 ) then
             cycle
          endif  ! DUMMY  
          if(i==6 ) then
             imode=occ
             ispc =ocar
          endif
          if(i==7 ) then
             cycle
          endif  ! sea salt has not number emission  
          if(i==8 ) then
             cycle
          endif  ! sea salt has not number emission   
          if(i==9 ) then
             imode=dd1
             ispc =dust
          endif
          if(i==10) then
             imode=dd2
             ispc =dust
          endif
          !print*,"nemis=",i,imode,ispccall flush(6) 
          if(numb_alloc(transport,imode) == off) cycle

          !-special treatment for akk/acc - sulf 
          if( (imode == akk .or. imode==acc) .and. ispc==sulf) then

             if(  imode == akk ) factor = 0.01 * 0.025* 96.07/62.66 * 1.e-9*recip_part_mass(i)   !-unit  #/m3/s
             if(  imode == acc ) factor = 0.99 * 0.025* 96.07/62.66 * 1.e-9*recip_part_mass(i)   !-unit  #/m3/s

             ! 1) case anthropogenic emissions
             isrc  = antro
             dummy_src(1:m1,1:m2,1:m3) = 0.5*( chem1_src_g(it1,isrc,SO2)%sc_src(1:m1,1:m2,1:m3) + &
                  chem1_src_g(it2,isrc,SO2)%sc_src(1:m1,1:m2,1:m3)  )

             k2 = 2! control for 2-dim aerosol (bioge, antro, marin, sdust)
             rt_aer(ispc) = r_antro * factor 
             k_src = 1
             call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                  ,dummy_src        &! source data
                  ,aer2_g(imode)%sc_t       &! tendency array
                  ,aer_nvert_src(ispc)       &! vertical size of source array
                  ,rt_aer(ispc)       &! diurnal cycle of emission
                  ,k_src       &! vertical level of source array (where data is stored)
                  ,k2       &! vertical level of tendency
                  ,rho_air)  ! air density (to convert to mixing ratio tendency)


             ! 2) case biomass burning emissions
             isrc  = bburn
             dummy_src(1:m1,1:m2,1:m3) = 0.5*( chem1_src_g(it1,isrc,SO2)%sc_src(1:m1,1:m2,1:m3) + &
                  chem1_src_g(it2,isrc,SO2)%sc_src(1:m1,1:m2,1:m3)  )
             rt_aer(ispc) = rt(bburn)* factor
             k2 = m1-1  !
             do k=2,k2
                k_src=k

                call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                     ,dummy_src                         &! source data
                     ,aer2_g(imode)%sc_t          &! tendency array
                     ,aer_nvert_src(ispc)          &! vertical size of source array
                     ,rt_aer(ispc)  &! diurnal cycle of emission
                     ,k_src           &! vertical level of source array (where data is stored)
                     ,k  &! vertical level of tendency
                     ,rho_air)   ! air density (to convert to mixing ratio tendency)
             enddo
             !-special treatment for boc - bcar 
          elseif( imode == boc .and. ispc==bcar) then

             !-unit  #/m3/s
             dummy_src(:,:,:) =1.e-9*(recip_part_mass(4) * aer1_g(boc,bcar)%sc_src(:,:,:) +&
                  recip_part_mass(5) * aer1_g(boc,ocar)%sc_src(:,:,:) )
             ! only biomass burning emissions
             isrc  = bburn
             rt_aer(ispc) = rt(bburn)
             k2 = m1-1  !
             do k=2,k2
                k_src=k
                call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                     ,dummy_src                         &! source data
                     ,aer2_g(imode)%sc_t          &! tendency array
                     ,aer_nvert_src(ispc)          &! vertical size of source array
                     ,rt_aer(ispc)  &! diurnal cycle of emission
                     ,k_src           &! vertical level of source array (where data is stored)
                     ,k  &! vertical level of tendency
                     ,rho_air)   ! air density (to convert to mixing ratio tendency)
             enddo


          else
             if(spc_aer_alloc(src      ,imode,ispc) == off) then
                print*," not aloc for",imode,ispc
                stop 7777
             endif

             !-unit #/m3/s
             dummy_src(:,:,:) =1.e-9*recip_part_mass(i) * aer1_g(imode,ispc)%sc_src(:,:,:)
             !       
             ! only anthropogenic emissions
             isrc  = antro
             k2     = 2     ! control for 2-dim aerosol (antro, marin, sdust)
             rt_aer(ispc) = r_antro 
             k_src = 1
             call source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                  ,dummy_src    &! source data
                  ,aer2_g(imode)%sc_t    &! tendency array
                  ,aer_nvert_src(ispc)      &! vertical size of source array
                  ,rt_aer(ispc)     &! diurnal cycle of emission
                  ,k_src    &! vertical level of source array (where data is stored)
                  ,k2    &! vertical level of tendency
                  ,rho_air)     ! air density (to convert to mixing ratio tendency)

          endif
       enddo
    endif ! end of aerosol section (aerosol model) 
    !
  end subroutine sources
  !------------------------------------------------------------------

  subroutine source_to_tend(m1,m2,m3,ia,iz,ja,jz,sc_src,sc_t,nvert,rt,k_src,k_tend,rho_air)

    ! original
    integer , intent(IN)    :: m1
    integer,  intent(IN)    :: m2
    integer,  intent(IN)    :: m3
    integer,  intent(IN)    :: ia
    integer,  intent(IN)    :: iz
    integer,  intent(IN)    :: ja
    integer,  intent(IN)    :: jz
    real,     intent(IN)    :: sc_src(nvert,m2,m3)
    real,     intent(INOUT) :: sc_t(m1,m2,m3)
    integer,  intent(IN)    :: nvert
    real,     intent(IN)    :: rt
    integer,  intent(IN)    :: k_src
    integer,  intent(IN)    :: k_tend
    real,     intent(IN)    :: rho_air(m1,m2,m3)

    !- also this routine assumes that the source term is expressed in terms of  density (kg/m3)
    !- and, then, the conversion to mixing ratio (division by air density) is needed.
    !- k_src express the level where the src data is stored in the sc_src array
    sc_t(k_tend,ia:iz,ja:jz) = sc_t(k_tend,ia:iz,ja:jz) +  sc_src(k_src,ia:iz,ja:jz)*rt /  &
                                !- the level of air density
                                !- must be the same of tendency
         rho_air(k_tend,ia:iz,ja:jz) 

  end subroutine source_to_tend


  !------------------------------------------------------------------

  subroutine source_to_tend_cycle(m1,m2,m3,ia,iz,ja,jz,sc_src1,sc_src2,sc_t,nvert, &
       rt,alfa,k_src,k_tend,rho_air)

    ! original
    integer          , intent(IN)    :: m1
    integer          , intent(IN)    :: m2
    integer          , intent(IN)    :: m3
    integer          , intent(IN)    :: ia
    integer          , intent(IN)    :: iz
    integer          , intent(IN)    :: ja
    integer          , intent(IN)    :: jz
    real             , intent(IN)    :: sc_src1(nvert,m2,m3)
    real             , intent(IN)    :: sc_src2(nvert,m2,m3)
    real             , intent(INOUT) :: sc_t(m1,m2,m3)
    integer          , intent(IN)    :: nvert
    double precision , intent(IN)    :: rt(m2,m3)
    real             , intent(IN)    :: alfa
    integer          , intent(IN)    :: k_src
    integer          , intent(IN)    :: k_tend
    real             , intent(IN)    :: rho_air(m1,m2,m3)

    integer i,j

    !- also this routine assumes that the source term is expressed in terms of  density (kg/m3)
    !- and, then, the conversion to mixing ratio (division by air density) is needed.
    !- k_src express the level where the src data is stored in the sc_src array
    ! if(nint(alfa) == 0) then

    !   do j=ja,jz ; do i=ia,iz
    !  
    !
    !   sc_t(k_tend,i,j) = sc_t(k_tend,i,j) +  (                        &
    !                      sc_src1(k_src,i,j)* ( 1.0D0- rt(i,j) )*alfa  +  &
    !                      sc_src2(k_src,i,j)* rt(i,j)               )/ &
    !                      rho_air(k_tend,i,j) !- the level of air density
    !                          !- must be the same of tendency
    !                     
    !   enddo;enddo

    !original
    if(alfa > 0.) then 


       sc_t(k_tend,ia:iz,ja:jz) = sc_t   (k_tend ,ia:iz,ja:jz) + (  &
            sc_src1(k_src  ,ia:iz,ja:jz) * ( 1.0D0- rt(ia:iz,ja:jz) )*alfa +     &
            sc_src2(k_src  ,ia:iz,ja:jz) *          rt(ia:iz,ja:jz)        )  /  &
            rho_air(k_tend ,ia:iz,ja:jz) !- the level of air density
       !- must be the same of tendency
    else

       sc_t(k_tend,ia:iz,ja:jz) = sc_t   (k_tend ,ia:iz,ja:jz) + (  &
!!!    sc_src1(k_src  ,ia:iz,ja:jz) * ( 1.0D0- rt(ia:iz,ja:jz) )*alfa  +  &
            sc_src2(k_src  ,ia:iz,ja:jz) *          real(rt(ia:iz,ja:jz))         )  /  &
            rho_air(k_tend ,ia:iz,ja:jz) !- the level of air density
       !- must be the same of tendency

    endif

  end subroutine source_to_tend_cycle
  !------------------------------------------------------------------
  subroutine source_to_tend_cycle_cetesb(m1,m2,m3,ia,iz,ja,jz,sc_src1,sc_src2,sc_t,nvert, &
       rt,alfa,k_src,k_tend,rho_air,rt2,ispc,spc_chem_name,chem_nspecies)

    !USE chem1_list, ONLY: spc_name

    ! original
    integer          , intent(IN)    :: chem_nspecies
    integer          , intent(IN)    :: m1
    integer          , intent(IN)    :: m2
    integer          , intent(IN)    :: m3
    integer          , intent(IN)    :: ia
    integer          , intent(IN)    :: iz
    integer          , intent(IN)    :: ja
    integer          , intent(IN)    :: jz
    real             , intent(IN)    :: sc_src1(nvert,m2,m3)
    real             , intent(IN)    :: sc_src2(nvert,m2,m3)
    real             , intent(INOUT) :: sc_t(m1,m2,m3)
    integer          , intent(IN)    :: nvert
    double precision , intent(IN)    :: rt(m2,m3)
    real             , intent(IN)    :: alfa
    integer          , intent(IN)    :: k_src
    integer          , intent(IN)    :: k_tend
    real             , intent(IN)    :: rho_air(m1,m2,m3)
    double precision , intent(IN)    :: rt2(m2,m3) 
    integer          , intent(IN)    :: ispc
    character(LEN=8), intent(IN),dimension(chem_nspecies) :: spc_chem_name
    integer i,j

    if (spc_chem_name(ispc)=='NO'.or.spc_chem_name(ispc)=='NO2') then
       sc_t(k_tend,ia:iz,ja:jz) = sc_t   (k_tend ,ia:iz,ja:jz) + (  &
!!!                           sc_src1(k_src  ,ia:iz,ja:jz) * ( 1.0D0- rt(ia:iz,ja:jz) )*alfa  +  &
            sc_src2(k_src  ,ia:iz,ja:jz) *          real(rt2(ia:iz,ja:jz))         )  /  &
            rho_air(k_tend ,ia:iz,ja:jz) !- the level of air density
       !- must be the same of tendency

    else

       sc_t(k_tend,ia:iz,ja:jz) = sc_t   (k_tend ,ia:iz,ja:jz) + (  &
!!!                           sc_src1(k_src  ,ia:iz,ja:jz) * ( 1.0D0- rt(ia:iz,ja:jz) )*alfa  +  &
            sc_src2(k_src  ,ia:iz,ja:jz) *          real(rt(ia:iz,ja:jz))         )  /  &
            rho_air(k_tend ,ia:iz,ja:jz) !- the level of air density
       !- must be the same of tendency

    endif

  end subroutine source_to_tend_cycle_cetesb
  !------------------------------------------------------------------
  ! NOT USED
  !
  !  subroutine invert_air_dens(m1,rho_air)
  !   implicit none
  !   integer m1
  !   real rho_air(m1)
  !   rho_air(:)=1./rho_air(:)
  ! 
  !  end subroutine invert_air_dens
  !------------------------------------------------------------------


  !------------------------------------------------------------------
  subroutine get_diurnal_cycle_normalized(m2,m3,ia,iz,ja,jz,dtlt,glat, &
       glon,imonth1,idate1,iyear1,itime1, &
       antro,bioge,pi180,nsrc,emiss_cycle)

    ! original
    integer , intent(IN) :: m2
    integer , intent(IN) :: m3
    integer , intent(IN) :: ia
    integer , intent(IN) :: iz
    integer , intent(IN) :: ja
    integer , intent(IN) :: jz
    real    , intent(IN) :: dtlt
    real    , intent(IN) :: glat(m2,m3)
    real    , intent(IN) :: glon(m2,m3)

    ! mem_grid
    integer , intent(IN) :: imonth1
    integer , intent(IN) :: idate1
    integer , intent(IN) :: iyear1
    integer , intent(IN) :: itime1

    ! mem_chem1
    integer , intent(IN) :: antro
    integer , intent(IN) :: bioge

    ! rconstants
    real    , intent(IN) :: pi180


    integer , intent(IN) :: nsrc
    type(cycle_emission), intent(INOUT) :: emiss_cycle(nsrc)

    integer :: jday,i,j

    !--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    !    INTEGER :: julday
    !--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    real :: solfac,tdec,sdec,cdec,declin,d0,d02,dayhr,radlat,cslcsd,snlsnd, &
         gglon,dayhrr,hrangl

    !- local var
    integer it
    real time_x,dt_x,cosz
    real dcnorma(m2,m3)
    real xxx

    !print*,'-----------------------------------------------------' 
    !print*,'getting diurnal/space variable emission rates - time=',time/3600. ;call flush(6)
    !print*,'-----------------------------------------------------'


    dt_x=0.
    time_x=0.
    dcnorma(:,:)=0.


    do it=1,nint(86400./dtlt) 

       jday = julday(imonth1,idate1,iyear1) ! don't change the position of this line

       jday = jday + nint(time_x/86400.)
       !      sdec - sine of declination, cdec - cosine of declination
       declin = -23.5 * cos(6.283 / 365. * (jday + 9)) * pi180
       sdec = sin(declin)
       cdec = cos(declin)

       ! Find the factor, solfac, to multiply the solar constant to correct
       ! for Earth's varying distance to the sun.

       !d0 = 6.2831853 * float(jday-1) / 365.
       !d02 = d0 * 2.
       !solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0)  &
       !     + 0.000719 * cos(d02) + 0.000077 * sin(d02)

       ! Find the hour angle, THEN get cosine of zenith angle.

       dayhr = time_x / 3600. + float(itime1/100) + float(mod(itime1,100)) / 60.

       do j = ja,jz
          do i = ia,iz
             radlat = glat(i,j) * pi180
             !IF (lonrad .eq. 0) radlat = centlat(1) * pi180
             !IF (radlat .eq. declin) radlat = radlat + 1.e-5
             cslcsd = cos(radlat) * cdec
             snlsnd = sin(radlat) * sdec
             !gglon = glon(i,j)
             !IF (lonrad .eq. 0) gglon = centlon(1)
             dayhrr = mod(dayhr+glon(i,j)/15.+24.,24.)
             hrangl = 15. * (dayhrr - 12.) * pi180
             cosz = snlsnd + cslcsd * cos(hrangl)
             !cosz = min(cosz+1.0E-10, 1.0) 

             dcnorma(i,j)=dcnorma(i,j)+max(0.,cosz)

          end do
       end do
       time_x=time_x+dtlt
    end do
    do j = ja,jz
       do i = ia,iz
          !- invert dcnorma to save computation time
          !      dcnorma(:,:)=1./(dcnorma(:,:)*dtlt)
          dcnorma(i,j)=1./(dcnorma(i,j)*dtlt)
       enddo
    enddo
    !- transfer the emission cycle to bioge and antro arrays only, for now  
    emiss_cycle(bioge)%dcnorma_inv( :,:) = dcnorma(:,:)
    emiss_cycle(antro)%dcnorma_inv( :,:) = dcnorma(:,:)

    !print*,'dcnorma_INV=',dcnorma(3,4); call flush(6)

  end subroutine get_diurnal_cycle_normalized
  !------------------------------------------------------------------
  subroutine interpolation_antro(h_hour,src1,src2,srcn1,srcn2,time1,time2)

    integer :: i, j

    real               :: newX, g, gTemp, no
    real               :: src1,src2,srcn1,srcn2,time1,time2
    real, dimension(25) :: f, x, c
    integer  :: INDICE
    real :: h_hour,diff

    !--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    !- Marcelo - curva para SaoPaulo
    !DATA f /0.027,0.016,0.019,0.023,0.027,0.031,0.035,0.040,0.044,0.048,0.051,0.054,0.056,0.057,0.057,&
    !                0.057,0.056,0.054,0.051,0.048,0.044,0.040,0.035,0.031,0.027/
    !DATA c /0.018,0.015,0.020,0.025,0.031,0.036,0.041,0.045,0.048,0.050,0.050,0.049,0.047,0.047,0.048,&
    !                0.050,0.053,0.056,0.058,0.056,0.052,0.045,0.035,0.026,0.018/ 
    !DATA x /0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24/
    !
    !- Madeleine - curva para Lima
    data f /0.007, 0.007, 0.003, 0.003, 0.003, 0.010, 0.041, 0.051, 0.082, 0.071, 0.055, 0.049, 0.058, &
         0.068, 0.063, 0.070, 0.081, 0.082, 0.063, 0.063, 0.026, 0.022, 0.014, 0.010, 0.007/
    data c /0.007, 0.007, 0.000, 0.000, 0.000, 0.009, 0.021, 0.034, 0.093, 0.071, 0.059, 0.057, 0.062, &
         0.076, 0.071, 0.077, 0.077, 0.074, 0.065, 0.062, 0.028, 0.026, 0.017, 0.009, 0.007/ 
    data x /0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24/
    !--(DMK-CCATT-BRAMS-5.0-OLD)------------------------------------------------------------------
    !DATA f /0.027,0.016,0.019,0.023,0.027,0.031,0.035,0.040,0.044,0.048,0.051,0.054,0.056,0.057,0.057,0.057,0.056,0.054,0.051,0.048,0.044,0.040,0.035,0.031,0.027/
    !DATA c /0.018,0.015,0.020,0.025,0.031,0.036,0.041,0.045,0.048,0.050,0.050,0.049,0.047,0.047,0.048,0.050,0.053,0.056,0.058,0.056,0.052,0.045,0.035,0.026,0.018/ 
    !DATA x /0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24/
    !--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------     

    call sortZ(24,x,h_hour,1,INDICE)


    diff = x(INDICE) - h_hour



    if (diff.gt.0) then
       time2 = x(INDICE)
       time1 = x(INDICE-1)
       src2 = c(INDICE)
       src1 = c(INDICE-1)
       srcn2 = f(INDICE)
       srcn1 = f(INDICE-1)

    else

       time1 = x(INDICE)
       time2 = x(INDICE+1)
       src1 = c(INDICE)
       src2 = c(INDICE+1)
       srcn1 = f(INDICE)
       srcn2 = f(INDICE+1)

    endif



  end subroutine interpolation_antro
  !------------------------------------------------------------------
  subroutine sortZ(qtLevels, levels, pLevel, order, outIndex)

    integer, intent(IN)                    :: qtLevels
    real, dimension(qtLevels), intent(IN)  :: levels
    real, intent(IN)                       :: pLevel
    integer, intent(IN)                    :: order
    integer, intent(OUT) :: outIndex

    integer                   :: i,     &
         j,     &
         minIdx
    real                      :: currMin
    real, dimension(qtLevels) :: ztemp


    ztemp = levels

    !A SIMPLE SORT THAT RETURNS THE 'order' INDEXES NEAREST OF 'pLevel'

    do j = 1, order

       currMin = abs(ztemp(1) - pLevel)
       minIdx  = 1

       do i = 2, qtLevels

          if( abs(ztemp(i) - pLevel) .le. currMin)then

             currMin = abs(ztemp(i) - pLevel)
             minIdx  = i

          end if
       end do

       outIndex = minIdx
       ztemp(minIdx) = -9E5

    end do

  end subroutine sortZ
  !------------------------------------------------------------------
  subroutine vert_dist_of_volcanic_emission(m1,m2,m3,ia,iz,ja,jz,isrctime,   &
       nzpmax,dzt,zt,zm,rtgt,topt,dn0,nsrc,     &
       chem1_src_g,                             &
       geoge,chem_nspecies,spc_chem_alloc, &
       src,off,transport,aer1_g,aerosol,   &
       aer_nspecies,spc_aer_alloc,nmodes,  &
       volc_mean_g,v_ash) 

    integer , intent(IN) :: m1
    integer , intent(IN) :: m2
    integer , intent(IN) :: m3
    integer , intent(IN) :: ia
    integer , intent(IN) :: iz
    integer , intent(IN) :: ja
    integer , intent(IN) :: jz
    integer , intent(IN) :: isrctime

    ! grid_dims
    integer , intent(IN) :: nzpmax

    ! mem_grid
    real    , intent(IN) :: dzt(nzpmax),zt(nzpmax),zm(nzpmax)
    real    , intent(IN) :: rtgt(m2,m3)
    real    , intent(IN) :: topt(m2,m3)
    real    , intent(IN) :: dn0(m1,m2,m3)

    ! chem1_list
    integer , intent(IN) :: chem_nspecies
    integer , intent(IN) :: spc_chem_alloc(6,chem_nspecies)
    integer , intent(IN) :: src
    integer , intent(IN) :: off
    integer , intent(IN) :: transport

    ! aer1_list
    integer , intent(IN) :: aer_nspecies
    integer , intent(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    integer , intent(IN) :: nmodes
    integer , intent(IN) :: v_ash


    ! mem_chem1
    integer              , intent(IN)    :: nsrc
    type(chem1_src_vars) , intent(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    integer              , intent(IN)    :: geoge

    ! mem_aer1
    type (aer1_vars) , intent(INOUT) :: aer1_g(nmodes,aer_nspecies)
    integer          , intent(IN)    :: aerosol

    ! mem_volc
    type (volc_mean_vars), intent(INOUT) :: volc_mean_g

    ! local var
    real :: dz_inv,dz,x1,actual_vent_elev
    integer :: i,j,ksrc,isrc,ispc,imode,itim,k_initial,k_final,kk4,ko,kl,k
    integer :: k_init_umbr
    !
    real, parameter :: percen_mass_umbrel=.75 
    real, parameter :: base_umbrel=.25 
    real, allocatable :: vert_mass_dist(:)


    ksrc=2          !surface level of emission in the model
    isrc = geoge    !this routine is only for volcanoes
    itim = actual_time_index(isrctime,isrc)


    allocate(vert_mass_dist(m1))

    !- chemistry section
    do j=ja,jz
       do i=ia,iz

          !- if there is not mass to distribute = > cycle
          x1=0.
          do ispc=1,chem_nspecies
             if(spc_chem_alloc(src,ispc) == off) cycle
             x1=max(x1,chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j))
          enddo
          do imode=1,nmodes         
             if(spc_aer_alloc(src,imode,v_ash) == off ) cycle
             x1=max(x1,aer1_g(imode,v_ash)%sc_src(1,i,j))
          enddo
          !if(maxval( chem1_src_g(itim,isrc,1:chem_nspecies,ng)%sc_src(1,i,j)) <1.e-10) cycle
          if(x1<1.e-16) cycle
          print*,' ==============================================================================' 
          print*,' -> active volcano found at grid point=',i,j
          print*,' -> vent topo plum-heigth=',volc_mean_g%vent_elev(i,j),topt(i,j),volc_mean_g%plum_heigth(i,j)

          !- convert units back to 1.e+9 kg/m^2
          !- (this routine is called after 'convert_to_tracer_density' routine
          !-  where the emission field was converted to kg/m^3*1.e9)
          !- so, we need to convert back to create the vertical emission field
          dz = rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1)) 
          do ispc=1,chem_nspecies
             if(spc_chem_alloc(src,ispc) == off) cycle

             chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) = &
                  chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) * dz
          enddo
          !-do the same for aerosols
          !-aerosol section (only for volcanoes)
          do imode=1,nmodes
             if(spc_aer_alloc(src,imode,v_ash) == off ) cycle
             aer1_g(imode,v_ash)%sc_src(1,i,j) =  &
                  aer1_g(imode,v_ash)%sc_src(1,i,j) * dz
          enddo

          !-
          !- performs the vertical distribution

          ! - initial settings 
          k_final    =2     !- volc cloud top    
          k_init_umbr=2     !- volc cloud base
          k_initial  =2     !- volc vent 

          !- check vertical position with the actual model topography
          !
          volc_mean_g%vent_elev  (i,j)=max( volc_mean_g%vent_elev  (i,j),topt(i,j) )

          !- here plum_height is converted to height above sea level:
          !
          !volc_mean_g%plum_heigth(i,j)=MAX( volc_mean_g%plum_heigth(i,j)                             ,topt(i,j) )
          volc_mean_g%plum_heigth(i,j)=max( volc_mean_g%plum_heigth(i,j)+volc_mean_g%vent_elev  (i,j),topt(i,j) )

          !- 
          if( abs(volc_mean_g%plum_heigth(i,j)-volc_mean_g%vent_elev(i,j)) .le. dz) then
             !- case of pseudo 'surface emission', but possibly elevated due mismatch between
             !- model grid-scale topography and the actual elevation of the volcano vent
             do k=m1-1,2,-1
                !print*,'x',zm(k)*rtgt(i,j)+topt(i,j) , (1.-base_umbrel),volc_mean_g(ng)%plum_heigth(i,j)
                if(zm(k)*rtgt(i,j)+topt(i,j) < volc_mean_g%plum_heigth(i,j))then
                   k_initial=k
                   exit
                endif
             enddo
             k_initial=max(2,k_initial)
             k_final  =k_initial
             vert_mass_dist(1:m1) = 0.
             vert_mass_dist(k_final) = 1.
             !print*,'pseudo sfc emis=',k_initial,k_final

          else
             !- possibly more than 1 vert emission layer, we will use the 'umbrella' shape for this case:
             !- top volc cloud above sea level   (meters)     
             do k=m1,2,-1

                if(zm(k)*rtgt(i,j)+topt(i,j) < volc_mean_g%plum_heigth(i,j))then
                   k_final=k+1
                   exit
                endif
             enddo
             k_final=max(2,k_final)

             if( k_final == 2 ) then                     ! surface emission
                k_init_umbr = 2
                k_initial = 2
                vert_mass_dist(1:m1) = 0.
                vert_mass_dist(k_final)=1.

             else ! aerial emission

                !- bottom of volc cloud above sea level
                do k=m1-1,2,-1
                   !print*,'x',zm(k)*rtgt(i,j)+topt(i,j) , (1.-base_umbrel),volc_mean_g(ng)%plum_heigth(i,j)
                   if(zm(k)*rtgt(i,j)+topt(i,j) < (1.-base_umbrel)*volc_mean_g%plum_heigth(i,j))then
                      k_init_umbr=k
                      exit
                   endif
                enddo
                k_init_umbr=max(2,k_init_umbr)
                !print*,'kin-kfin=' ,k_init_umbr,k_final

                !- vertical mass distribution (initialization)
                vert_mass_dist(1:m1) = 0.

                !- part 1
                !- parabolic vertical distribution between k_init_umbr and k_final
                kk4 = k_final-k_init_umbr+2
                do ko=2,kk4-1
                   kl=ko+k_init_umbr-1
                   vert_mass_dist(kl) = 6.*percen_mass_umbrel* float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
                enddo

                if(sum(vert_mass_dist(1:m1)) .ne. percen_mass_umbrel) then
                   x1= ( percen_mass_umbrel- sum(vert_mass_dist(1:m1)) )/float(k_final-k_init_umbr+1)
                   do ko=k_init_umbr,k_final
                      vert_mass_dist(ko) = vert_mass_dist(ko)+ x1 !- values between 0 and 1.
                   enddo
                   !print*,'new mass=',sum(vmd)*100.,x1
                   !pause
                endif
                !
                !- part 2
                !- determine the actual level of the volcanoe vent
                do k=m1-1,2,-1
                   if(zm(k)*rtgt(i,j)+topt(i,j) < volc_mean_g%vent_elev(i,j))then
                      k_initial=k
                      exit
                   endif
                enddo
                k_initial=max(2,k_initial)

                !- linear detrainment from vent to base of umbrella
                do ko=k_initial,k_init_umbr
                   !vert_mass_dist(ko)=float(ko)/float(k_init_umbr)
                   vert_mass_dist(ko)=(log((float(ko)/float((k_initial-1)))))**2./float(k_init_umbr)
                enddo
                x1=sum(vert_mass_dist(2:k_init_umbr))
                do ko=2,k_init_umbr
                   vert_mass_dist(ko)=(1.-percen_mass_umbrel)*vert_mass_dist(ko)/x1
                enddo

                ! if(k_init_umbr > 2 ) vert_mass_dist(k_init_umbr)=max(vert_mass_dist(k_init_umbr),vert_mass_dist(k_init_umbr-1))
                !- check final mass conservation
                if(sum(vert_mass_dist(1:m1)) .ne. 1) then
                   x1= ( 1- sum(vert_mass_dist(1:m1)) )/float(k_final-k_init_umbr+1)
                   do ko=k_init_umbr,k_final
                      vert_mass_dist(ko) = vert_mass_dist(ko)+ x1 !- values between 0 and 1.
                   enddo
                   !print*,'new mass=',sum(vmd)*100.,x1
                   !pause
                endif

             endif
          endif

          print*,'vert mass distr=',sum(vert_mass_dist)*100.,i,j,volc_mean_g%plum_heigth(i,j),volc_mean_g%vent_elev(i,j)

          do ispc=1,chem_nspecies

             if(spc_chem_alloc(src,ispc) == off) cycle

             do k=k_initial,k_final

                !dz= rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1)) 
                dz_inv = dzt(k)/rtgt(i,j)

                !- convert from kg/m^2  to  density (kg[gas]/m^3*1.e9)

                chem1_src_g(itim,isrc,ispc)%sc_src(k,i,j) = &
                     chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) * dz_inv * vert_mass_dist(k)

                !IF(chem1_src_g(itim,isrc,ispc)%sc_src(k,i,j)>0.) &
                !  PRINT*,'ijk sc chem=',i,j,k,vert_mass_dist(k),chem1_src_g(itim,isrc,ispc)%sc_src(k,i,j)


             enddo
          enddo

          !-aerosol section (only for volcanoes)
          do imode=1,nmodes

             if(spc_aer_alloc(src,imode,v_ash) == off ) cycle

             !!    if(imode==1) print*, 'total ash mass kg/m2= ',sum( aer1_g(1:nmodes,v_ash)%sc_src(1,i,j) )*1.e-9
             print*,'----------- imode------ k ----- vert dist ------- mass of ash kg/m^3*1.e9 --------' 
             do k=k_initial,k_final

                !dz= rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1)) 
                dz_inv = dzt(k)/rtgt(i,j)

                !- convert from kg/m^2  to  density (kg[gas]/m^3*1.e9)
                aer1_g(imode,v_ash)%sc_src(k,i,j) =  &
                     aer1_g(imode,v_ash)%sc_src(1,i,j) * dz_inv * vert_mass_dist(k)

                !IF(aer1_g(imode,v_ash)%sc_src(k,i,j)>0. .and. imode==1) &
                print*,imode,k,vert_mass_dist(k),aer1_g(imode,v_ash)%sc_src(k,i,j)


             enddo
          enddo
          print*,' ==============================================================================' 

          !- big-loop X-Y
       enddo
    enddo
    deallocate(vert_mass_dist)

  end subroutine vert_dist_of_volcanic_emission



  subroutine StoreNamelistFileAtChemSources(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    srcmapfn = oneNamelistFile%srcmapfn
    def_proc_src = oneNamelistFile%def_proc_src
  end subroutine StoreNamelistFileAtChemSources
end module chem_sources
