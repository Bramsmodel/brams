!--------------------------------------------------------------------------------!
! Cumulus Parameterizations G3d and GF                                          !
! Implemented in BRAMS by Saulo Freitas @ Feb/2012 /  Feb/2021                               !
! Rafael Mello: included parallelization for spread/g3d_smoothh arrays               !
!--------------------------------------------------------------------------------!
module ModCuParGrell3

  use ModChemConvTransp, only: &
       trans_conv_mflx_GF
  
  use ModRadiateFields, only: &
       RadiateFields
  
  use ModGrid, only: &
       Grid

  use mem_stilt, only: &
       imassflx

  use ModRstilt, only: &
       prep_convflx_to_stilt

  use ModNamelistFile, only: &
       NamelistFile

  use ModMessageSet, only: &
       PostSendRecvMsgs,    &
       WaitSendRecvMsgs

  use mem_tend, only: &
       tend

  use mem_cuparm, only: &
       cuparm_g,&
       cuparm_g_sh

  use node_mod, only: &
       mynum,   &   ! INTENT(IN)
       mxp,     &   ! INTENT(IN)
       myp,     &   ! INTENT(IN)
       mzp,     &   ! INTENT(IN)
       ia,      &   ! INTENT(IN)
       iz,      &   ! INTENT(IN)
       ja,      &   ! INTENT(IN)
       jz,      &   ! INTENT(IN)
       i0,      &   ! INTENT(IN)
       j0,      &   ! INTENT(IN)
       ibcon        ! INTENT(IN)

  use mem_grid, only: &
       if_adap, &
       jdim, &
       zt, &
       zm, &
       dzm, &
       dzt, &
       hw4, &
       itopo, &
       dyncore_flag, &
       time,    &   ! INTENT(IN)
       initial, &   ! INTENT(IN)
       dtlt, &   ! INTENT(IN)
       itime1,  &   ! INTENT(IN)
       ngrid,   &   ! INTENT(IN)
       grid_g,  &   ! INTENT(IN)
       dtlongn, &   ! INTENT(IN)
       deltaxn, &   ! INTENT(IN)
       deltayn, &   ! INTENT(IN)
       npatch,  &   ! INTENT(IN)
       ztn,     &   ! INTENT(IN)
       zmn,     &   ! INTENT(IN)
       akminvar,&   ! INTENT(IN)
       nxtnest

  use mem_varinit, only: &
       nudlat

  use rconstants, only: &
       tkmin, &
       cpi

  use ModMicroFields, only: &
       MicroFields

  use io_params, only: &
       frqanl

  use mem_leaf, only: &
       leaf_g, &
       isfcl

  use ModMicControl, only: &
       MicControl
  
  use mem_grell_param, only: &
       mgmxp, &
       mgmyp, &
       mgmzp, &
       maxiens, &
       ngrids_cp,&
       maxens_g3d ,                        & !INTENT(IN)
       maxens2_g3d,                        & !INTENT(IN)
       maxens3_g3d,                        & !INTENT(IN)
       ensdim_g3d ,                        & !INTENT(IN)
       maxens ,                            & !INTENT(IN)
       maxens2,                            & !INTENT(IN)
       maxens3,                            & !INTENT(IN)
       ensdim ,                            & !INTENT(IN)
       icoic

  use mem_scratch1_grell, only: &
       ierr4d, &
       jmin4d, &
       kdet4d, &
       k224d, &
       kbcon4d, &
       ktop4d, &
       kpbl4d, &
       kstabi4d, &
       kstabm4d, &
       xmb4d, &
       edt4d, &
       pwav4d, &
       zup5d, &
       zdn5d, &
       iruncon, &
       pcup5d, &
       prup5d, &
       prdn5d, &
       clwup5d, &
       tup5d, &
       enup5d, &
       endn5d, &
       deup5d, &
       dedn5d, &
       zcup5d, &
       up_massdetr5d, &
       up_massentr5d, &
       dd_massdetr5d, &
       dd_massentr5d, &
       conv_cld_fr5d, &
       sigma4d, &
       klcl4d, &
       cprr4d


  use mem_grell, only: &
       cuforc_g, &
       cuforc_sh_g

  use mem_carma, only: &
       carma

  use module_cu_g3, only: &
       G3DRV

  use module_cu_gf, only: &
       GFDRV

  use module_cu_gf2, only: &
       GFDRV2

  use Phys_const, only: &
       cp, &
       p00, &
       tcrit, &
       g, &
       cpor , &
       XL, &
       rm, &
       rgas

  use ConvPar_GF_GEOS5, only: &
       make_IceNumber, &
       fract_liq_f, &
       make_DropletNumber, & 
       GF_GEOS5_DRV, &
       deep, &
       shal, &
       mid, &
       nmp, &
       lsmp, &
       cnmp,  &
       GF_convpar_init, &
       apply_sub_mp, &
       icumulus_gf, &
       liq_ice_number_conc

  use ccatt_start, only: &
       ccatt

  use mem_chem1, only: &
       chemistry

  use ModVarTable, only: &
       VarTable, &
       InsertVarTable

  use grid_dims, only: &
       nzpmax

  use mem_scratch, only: &
       vctr1, &
       vctr2

  use ModBasicFields, only: &
       BasicFields

  use ModRadvc, only: &
       advtndc, &
       fa_preptc, &
       fa_xc, &
       fa_yc, &
       fa_zc


  implicit none

  include "constants.h"

  private
  public :: cupar2mcphysics
  public :: prepare_lsf
  public :: init_weights
  public :: StoreNamelistFileAtCup_grell3
  public :: cuparm_grell3_catt
  public :: alloc_grell3
  public :: nullify_grell3
  public :: filltab_grell3
  public :: g3d_ens_g
  public :: g3d_ensm_g
  public :: g3d_g
  public :: g3dm_g
  public :: train_dim
  
  type g3d_ens_vars
     real, pointer, contiguous :: apr(:,:) => null()
     real, pointer, contiguous :: accapr(:,:) => null()
     real, pointer, contiguous :: weight(:,:) => null()
  end type g3d_ens_vars

  type(g3d_ens_vars), pointer, contiguous :: g3d_ens_g(:,:)
  type(g3d_ens_vars), pointer, contiguous :: g3d_ensm_g(:,:)

  type g3d_vars
     real, pointer, contiguous :: xmb_deep(:,:) => null()
     real, pointer, contiguous :: xmb_deep_dd(:,:) => null()
     real, pointer, contiguous :: err_deep(:,:) => null()
     real, pointer, contiguous :: xmb_shallow(:,:) => null()
     real, pointer, contiguous :: cugd_ttens(:,:,:) => null()
     real, pointer, contiguous :: cugd_qvtens(:,:,:) => null()
     real, pointer, contiguous :: thsrc(:,:,:) => null()
     real, pointer, contiguous :: rtsrc(:,:,:) => null()
     real, pointer, contiguous :: clsrc(:,:,:) => null()
     real, pointer, contiguous :: nlsrc(:,:,:) => null()
     real, pointer, contiguous :: nisrc(:,:,:) => null()
     real, pointer, contiguous :: usrc(:,:,:) => null()
     real, pointer, contiguous :: vsrc(:,:,:) => null()
     real, pointer, contiguous :: mup(:,:,:) => null()
     real, pointer, contiguous :: mdd(:,:,:) => null()
     real, pointer, contiguous :: mupsh(:,:,:) => null()
  end type g3d_vars

  type(g3d_vars), pointer, contiguous :: g3d_g(:)
  type(g3d_vars), pointer, contiguous :: g3dm_g(:)

  integer ::    ids,ide, jds,jde, kds,kde            &
       ,ims,ime, jms,jme, kms,kme            &
       ,ips,ipe, jps,jpe, kps,kpe            &
       ,its,ite, jts,jte, kts,kte

  integer           ::  imomentum=1
  integer           ::  ishallow_g3

  !- define if the training will be used or not
  integer,parameter:: training=0
  character(len=255) :: g3d_training_file

  !- define if the lateral subsidence spread will be done or not
  integer :: g3d_spread ! 1=ON, 0=OFF
  integer :: cugd_avedx

  !- define if the horizontal smoothing is to be done or not
  integer :: g3d_smoothh! 1=ON, 0=OFF

  !- define if the vertical smoothing is to be done or not
  integer :: g3d_smoothv! 1=ON, 0=OFF

  !- number of members of prec ensemble
  integer,parameter :: train_dim= 5

  character(LEN=6),parameter,dimension(train_dim) :: pre_name=(/ &
       'apr_gr' & !
       ,'apr_w ' & !
       ,'apr_mc' & !
       ,'apr_st' & !
       ,'apr_as' & !
       /)

  integer,parameter :: apr_gr=001
  integer,parameter :: apr_w =002
  integer,parameter :: apr_mc=003
  integer,parameter :: apr_st=004
  integer,parameter :: apr_as=005


  integer,parameter :: CPTIME = 0. !orig: CPTIME = 7200.

  integer,parameter :: i_forcing = 1

  integer,parameter :: autoconv = 1!2 ! =1, Kessler
  ! =2, Berry
  integer,parameter :: aerovap = 1!3  ! =1, orig
  ! =2, mix orig+new
  ! =3, new
  !- direct link cupar-microphysics
  logical,parameter :: do_cupar_mcphys_coupling = .true.
  !- read GF namelist
  logical :: read_GF_ConvPar_nml =.true. 

  character(len=*), parameter :: header="**(cup_grell3)**"

contains
  !-----------------------------------------
  subroutine nullify_grell3(g3d_ens,g3d,ndim_train)
    integer, intent(in) ::ndim_train
    type (g3d_ens_vars),dimension(ndim_train) :: g3d_ens
    type (g3d_vars) :: g3d
    integer i

    do i=1,ndim_train
       if (associated(g3d_ens(i)%apr))    nullify (g3d_ens(i)%apr)
       if (associated(g3d_ens(i)%accapr)) nullify (g3d_ens(i)%accapr)
       if (associated(g3d_ens(i)%weight)) nullify (g3d_ens(i)%weight)
    enddo

    if (associated(g3d%xmb_deep))    nullify (g3d%xmb_deep)
    if (associated(g3d%xmb_deep_dd)) nullify (g3d%xmb_deep_dd)
    if (associated(g3d%err_deep))    nullify (g3d%err_deep)
    if (associated(g3d%xmb_shallow)) nullify (g3d%xmb_shallow)

    if (associated(g3d%cugd_ttens))  nullify (g3d%cugd_ttens)
    if (associated(g3d%cugd_qvtens)) nullify (g3d%cugd_qvtens)

    if (associated(g3d%thsrc)) nullify (g3d%thsrc)
    if (associated(g3d%rtsrc)) nullify (g3d%rtsrc)
    if (associated(g3d%clsrc)) nullify (g3d%clsrc)
    if (associated(g3d%nlsrc)) nullify (g3d%nlsrc)
    if (associated(g3d%nisrc)) nullify (g3d%nisrc)
    if (associated(g3d%usrc )) nullify (g3d%usrc)
    if (associated(g3d%vsrc )) nullify (g3d%vsrc)
    if (associated(g3d%mup  )) nullify (g3d%mup)
    if (associated(g3d%mupsh)) nullify (g3d%mupsh)
    if (associated(g3d%mdd  )) nullify (g3d%mdd)

  end subroutine nullify_grell3
  !-----------------------------------------
  subroutine alloc_grell3(oneNamelistFile, &
       g3d_ens,g3d, m1, m2, m3, ng,ndim_train)
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type (g3d_ens_vars),dimension(ndim_train) :: g3d_ens
    type (g3d_vars) :: g3d
    integer, intent(in) :: m1, m2, m3, ng,ndim_train
    integer :: i

    if(oneNamelistFile%nnqparm(ng) == 3 .or. &
         oneNamelistFile%nnqparm(ng) == 6 .or. &
         oneNamelistFile%nnqparm(ng) == 5  .or. &
         oneNamelistFile%nnqparm(ng) == 8) then
       do i=1,ndim_train
          allocate(g3d_ens(i)%apr   (m2,m3))
          g3d_ens(i)%apr    =0.0
          allocate(g3d_ens(i)%accapr(m2,m3))
          g3d_ens(i)%accapr =0.0
          allocate(g3d_ens(i)%weight(m2,m3))
          g3d_ens(i)%weight =0.0
       enddo
       allocate (g3d%cugd_ttens (m1, m2, m3))
       g3d%cugd_ttens =0.0
       allocate (g3d%cugd_qvtens(m1, m2, m3))
       g3d%cugd_qvtens=0.0
    endif

    allocate (g3d%xmb_deep   (m2,m3))
    g3d%xmb_deep   =0.0
    allocate (g3d%xmb_deep_dd(m2,m3))
    g3d%xmb_deep_dd=0.0
    allocate (g3d%err_deep   (m2,m3))
    g3d%err_deep   =0.0
    allocate (g3d%xmb_shallow(m2,m3))
    g3d%xmb_shallow=0.0

    allocate (g3d%thsrc(m1, m2, m3))
    g3d%thsrc=0.0
    allocate (g3d%rtsrc(m1, m2, m3))
    g3d%rtsrc=0.0
    allocate (g3d%clsrc(m1, m2, m3))
    g3d%clsrc=0.0
    allocate (g3d%nlsrc(m1, m2, m3))
    g3d%nlsrc=0.0
    allocate (g3d%nisrc(m1, m2, m3))
    g3d%nisrc=0.0

    if( (imomentum==0 .or. imomentum==1 ) .and. oneNamelistFile%nnqparm(ng) >= 4) then
       allocate (g3d%usrc(m1, m2, m3))
       g3d%usrc=0.0
       allocate (g3d%vsrc(m1, m2, m3))
       g3d%vsrc=0.0
    endif
    allocate (g3d%mupsh(m1, m2, m3))
    g3d%mupsh=0.0
    allocate (g3d%mup  (m1, m2, m3))
    g3d%mup  =0.0
    allocate (g3d%mdd  (m1, m2, m3))
    g3d%mdd  =0.0

  end subroutine alloc_grell3

  !-----------------------------------------
  subroutine filltab_grell3(oneVarTable, oneVarTableSize, &
       g3d_ens, g3d, g3d_ensm, g3dm, ndim_train, nnqparm, &
       ng_cp, imean)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(g3d_ens_vars), pointer, intent(in) :: g3d_ens(:,:)
    type(g3d_ens_vars), pointer, intent(in) :: g3d_ensm(:,:)
    type(g3d_vars), pointer, intent(in) :: g3d(:)
    type(g3d_vars), pointer, intent(in) :: g3dm(:)
    integer, intent(in) :: ndim_train
    integer, intent(in) :: nnqparm
    integer, intent(in) :: ng_cp
    integer, intent(in) :: imean
    
    integer :: i
    character (len=4) :: arrprop
    character(len=*), parameter :: h="**(filltab_grell3)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(g3d)) then
       call fatal_error(h//" g3d not associated")
    else if (.not. associated(g3dm)) then
       call fatal_error(h//" g3dm not associated")
    else if (.not. associated(g3d_ens)) then
       call fatal_error(h//" g3d_ens not associated")
    else if (.not. associated(g3d_ensm)) then
       call fatal_error(h//" g3d_ensm not associated")
    end if
    
    ! Fill pointers to arrays into variable tables

    do i=1,ndim_train
       if (associated(g3d_ens(i,ng_cp)%apr)) then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  g3d_ens(i,ng_cp)%apr, &
                  trim(pre_name(i))//' :2:hist:mpti:mpt3', &
                  g3d_ensm(i,ng_cp)%apr, imean)
       end if

       if (associated(g3d_ens(i,ng_cp)%accapr)) then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  g3d_ens(i,ng_cp)%accapr, &
                  'acc'//trim(pre_name(i)(2:len_trim(pre_name(i))))//' :2:hist:anal:mpti:mpt3', &
                  g3d_ensm(i,ng_cp)%accapr, imean)
       end if
    end do

    do i=1,ndim_train
       if (associated(g3d_ens(i,ng_cp)%weight)) then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  g3d_ens(i,ng_cp)%weight, &
                  'weight'//trim(pre_name(i)(4:len_trim(pre_name(i))))//' :2:hist:anal:mpti:mpt3', &
                  g3d_ensm(i,ng_cp)%weight, imean)
       end if
    end do

    if (associated(g3d(ng_cp)%xmb_deep)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%xmb_deep, &
               'MFUP :2:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%xmb_deep, imean)
    end if

    if (associated(g3d(ng_cp)%xmb_deep_dd)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%xmb_deep_dd, &
               'MFDD :2:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%xmb_deep_dd, imean)
    end if

    if (associated(g3d(ng_cp)%err_deep)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%err_deep, &
               'XIERR :2:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%err_deep, imean)
    end if

    if (associated(g3d(ng_cp)%xmb_shallow)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%xmb_shallow, &
               'MFSH :2:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%xmb_shallow, imean)
    end if


    !- 3D Arrays

    !- define if the arrays will exchange 1 row x 1 line (not in use anymore)
    arrprop=''

    if (associated(g3d(ng_cp)%cugd_ttens)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%cugd_ttens, &
               'TTENS :3:hist:anal:mpti:mpt3'//trim(arrprop), &
               g3dm(ng_cp)%cugd_ttens, imean)
    end if

    if (associated(g3d(ng_cp)%cugd_qvtens)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%cugd_qvtens , &
               'QVTTENS :3:hist:anal:mpti:mpt3'//trim(arrprop), &
               g3dm(ng_cp)%cugd_qvtens, imean)
    end if

    if (associated(g3d(ng_cp)%thsrc)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%thsrc, &
               'THSRC :3:hist:anal:mpti:mpt3'//trim(arrprop), &
               g3dm(ng_cp)%thsrc, imean)
    end if

    if (associated(g3d(ng_cp)%rtsrc)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%rtsrc, &
               'RTSRC :3:hist:anal:mpti:mpt3'//trim(arrprop), &
               g3dm(ng_cp)%rtsrc, imean)
    end if

    !- this array does not need to be parallelized (only column)
    if (associated(g3d(ng_cp)%clsrc)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%clsrc, &
               'CLSRC :3:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%clsrc, imean)
    end if

    if (associated(g3d(ng_cp)%nlsrc)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%nlsrc, &
               'NLSRC :3:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%nlsrc, imean)
    end if

    if (associated(g3d(ng_cp)%nisrc)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%nisrc, &
               'NISRC :3:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%nisrc, imean)
    end if


    if(imomentum==1 .and. nnqparm >= 4) then
       !- these arrays does not need to be parallelized (only column)
       if (associated(g3d(ng_cp)%usrc)) then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  g3d(ng_cp)%usrc, &
                  'USRC :3:hist:anal:mpti:mpt3', &
                  g3dm(ng_cp)%usrc, imean)
       end if

       if (associated(g3d(ng_cp)%vsrc)) then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  g3d(ng_cp)%vsrc, &
                  'VSRC :3:hist:anal:mpti:mpt3', &
                  g3dm(ng_cp)%vsrc, imean)
       end if

    endif

    if (associated(g3d(ng_cp)%mupsh)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%mupsh, &
               'ZMFSH :3:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%mupsh, imean)
    end if

    if (associated(g3d(ng_cp)%mup)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%mup, &
               'ZMFUP :3:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%mup, imean)
    end if

    if (associated(g3d(ng_cp)%mdd)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               g3d(ng_cp)%mdd, &
               'ZMFDD :3:hist:anal:mpti:mpt3', &
               g3dm(ng_cp)%mdd, imean)
    end if

  end subroutine filltab_grell3
  !-------------------------------------------------------------

  subroutine CUPARM_GRELL3_CATT(oneGrid, iens,iinqparm,iinshcu)
    integer, intent(IN) :: iens,iinqparm,iinshcu
    type(Grid), pointer :: oneGrid ! intent(in)
    
    integer :: i,j,k
    real :: grid_length,theta2temp

    real, dimension( mxp , myp ) :: aot500,temp2m

    !----GF-GEOS-5 -------------------------------->
    integer, parameter :: mtp=1 !tmp
    integer, dimension(mzp)   :: flip
    integer, dimension(mxp,myp) :: kpbl,do_this_column
    real   , dimension(mtp)   :: FSCAV_INT
    real   , dimension(mxp,myp) :: CNV_FRC,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC
    real   , dimension(mzp,mxp,myp ) ::  zm3d&
         ,zt3d&
         ,dm3d&
         ,up&
         ,vp&
         ,wp&
         ,rvap&
         ,temp&
         ,press&
         , gsf_t& ! grid-scale forcing for temp
         , gsf_q& ! grid-scale forcing fo rv
         ,sgsf_t& ! sub-grid scale forcing for temp
         ,sgsf_q  ! sub-grid scale forcing for rv

    real,  dimension(mzp , mxp, myp ) ::          &
         buoy_exc    &
         ,advf_t   &
         ,SRC_BUOY    &
         ,REVSU_GF    &
         ,PRFIL_GF    & 
         ,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF

    real,  dimension(nmp, mzp , mxp, myp ) ::       &
         mp_ice   &
         ,mp_liq   &
         ,mp_cf    

    real,  dimension(nmp, mzp , mxp, myp ) ::       &
         SUB_MPQI & ! subsidence transport applied to ice mix ratio
         ,SUB_MPQL & ! subsidence transport applied to cloud mix ratio
         ,SUB_MPCF   ! subsidence transport applied to cloud fraction

    real   , dimension(mxp,myp)  :: var2d,col_sat ,LIGHTN_DENS,stochastic_sig

    real   ,dimension(mtp,mzp,mxp,myp)  :: SRC_CHEM
    !---
    real   ,dimension(mxp,myp,mzp,mtp)  :: TRACER  !geos-5 data structure
    !---
    real   , dimension(mxp,myp)  :: dx2d ,lons,lats,sfc_press,xland

    integer :: kr,n,i1,i2,j1,j2
    integer :: gridID
    integer :: idiffk
    real :: akmin
    integer, dimension(mxp,myp,maxiens) ::     &
         ierr4d_tmp      &
         ,jmin4d_tmp      &
         ,klcl4d_tmp      &
         ,k224d_tmp       &
         ,kbcon4d_tmp     &
         ,ktop4d_tmp      &
         ,kstabi4d_tmp     &
         ,kstabm4d_tmp

    real,dimension(mxp,myp,maxiens)     ::     &
         cprr4d_tmp      &
         ,xmb4d_tmp       &
         ,edt4d_tmp       &
         ,pwav4d_tmp      &
         ,sigma4d_tmp
    real,dimension(mxp,myp,mzp,maxiens) ::     &
         pcup5d_tmp      &
         ,up_massentr5d_tmp    &
         ,up_massdetr5d_tmp    &
         ,dd_massentr5d_tmp    &
         ,dd_massdetr5d_tmp    &
         ,zup5d_tmp       &
         ,zdn5d_tmp       &
         ,prup5d_tmp      &
         ,prdn5d_tmp      &
         ,clwup5d_tmp     &
         ,tup5d_tmp       &
         ,conv_cld_fr5d_tmp

    !if(initial.eq.2.and.time.lt.cptime) return
    !if(initial.eq.2.and.time.lt.dtlt) return

    gridId=oneGrid%Id
    idiffk=oneGrid%oneNamelistFile%idiffk(gridId)
    akmin=oneGrid%oneNamelistFile%akmin(gridId)

    if(mod(time,oneGrid%oneNamelistFile%confrq) < dtlt  .or. &
         time < 0.01 .or. &
         abs(time-cptime) < 0.01) then

       !-start convective transport of tracers
       iruncon=1
       g3d_g(ngrid)%thsrc       = 0.0
       g3d_g(ngrid)%rtsrc       = 0.0
       g3d_g(ngrid)%clsrc       = 0.0
       g3d_g(ngrid)%cugd_ttens  = 0.0
       g3d_g(ngrid)%cugd_qvtens = 0.0
       g3d_g(ngrid)%mupsh       = 0.0
       g3d_g(ngrid)%mup         = 0.0
       g3d_g(ngrid)%mdd         = 0.0
       cuparm_g(ngrid)%conprr   = 0.0

       if(liq_ice_number_conc > 0) then 
          g3d_g(ngrid)%nlsrc     = 0.0
          g3d_g(ngrid)%nisrc     = 0.0
       endif
       if(imomentum == 1 .and. oneGrid%oneNamelistFile%nnqparm(ngrid) >= 4) then
          g3d_g(ngrid)%usrc      = 0.0
          g3d_g(ngrid)%vsrc      = 0.0
       endif
       ishallow_g3=0

       if(i_forcing /= 1) then
          !call check(mzp * myp * mzp,tend%THT ,cuforc_g(ngrid)%lsfth ,tend%RTT  ,cuforc_g(ngrid)%lsfrt)
          call atob(mxp * myp * mzp,tend%THT  ,cuforc_g(ngrid)%lsfth     )
          call atob(mxp * myp * mzp,tend%RTT  ,cuforc_g(ngrid)%lsfrt     )
       endif

       !- converting WRF setting to BRAMS
       ids=1
       ide=mxp
       jds=1
       jde=myp
       kds=1
       kde=mzp
       ims=1
       ime=mxp
       jms=1
       jme=myp
       kms=1
       kme=mzp
       ips=ia+1
       ipe=iz-2
       jps=ja+1
       jpe=jz-2
       kps=1
       kpe=mzp
       its=ia
       ite=iz
       jts=ja
       jte=jz
       kts=1
       kte=mzp-1

       grid_length=sqrt(deltaxn(ngrid)*deltayn(ngrid))

       !-------------------------------------------------------------

       if(iinqparm==3) then  ! G3d scheme
          !
          !- lateral spreading
          if(g3d_spread == 0 )cugd_avedx=1
          if(g3d_spread == 1 )cugd_avedx=3

          if(oneGrid%oneNamelistFile%ilwrtyp==4 .or. oneGrid%oneNamelistFile%iswrtyp==4) then
             aot500(:,:)=carma(ngrid)%aot(:,:,11)
          else
             aot500(:,:)=0.0
          end if

          call G3DRV( mynum,i0,j0,time          &
               ,dtlt                         & !
               ,grid_length                  & !
               ,autoconv                     & !
               ,aerovap                      & !
               ,oneGrid%oneBasicFields%dn0           & !
               ,cuparm_g(ngrid)%CONPRR       & !
               ,oneGrid%oneBasicFields%up            & !
               ,oneGrid%oneBasicFields%vp            & !
               ,oneGrid%oneBasicFields%theta         & !
               ,oneGrid%oneBasicFields%thp           & !
               ,oneGrid%oneBasicFields%pp            & !
               ,oneGrid%oneBasicFields%pi0           & !
               ,oneGrid%oneBasicFields%wp            & !
               ,oneGrid%oneBasicFields%rv            & !
               ,grid_g(ngrid)%RTGT           & !
               ,tend%PT                      & !
               ,XL                    & !
               ,CP    & !
               ,G    & !
               ,rm                           &
               ,p00                          &
               ,cpor                         & !
               ,g3d_ens_g(apr_gr,ngrid)%apr  &
               ,g3d_ens_g(apr_w ,ngrid)%apr  &
               ,g3d_ens_g(apr_mc,ngrid)%apr  &
               ,g3d_ens_g(apr_st,ngrid)%apr  &
               ,g3d_ens_g(apr_as,ngrid)%apr  &
               ,g3d_g(ngrid)%xmb_deep        &
               ,g3d_g(ngrid)%xmb_shallow     &
                                !
               ,g3d_ens_g(apr_gr,ngrid)%weight &
               ,g3d_ens_g(apr_w ,ngrid)%weight &
               ,g3d_ens_g(apr_mc,ngrid)%weight &
               ,g3d_ens_g(apr_st,ngrid)%weight &
               ,g3d_ens_g(apr_as,ngrid)%weight &
                                !
               ,training &
                                !
               ,grid_g(ngrid)%topt              &
               ,leaf_g(ngrid)%patch_area        &
               ,npatch                          &
               ,oneGrid%oneRadiateFields%rshort         &
               
               ,cugd_avedx&
               ,imomentum          &
               ,ensdim,maxiens,maxens,maxens2,maxens3,icoic      &
               ,ishallow_g3                                      &
               ,ids,ide, jds,jde, kds,kde                        &
               ,ims,ime, jms,jme, kms,kme                        &
               ,ips,ipe, jps,jpe, kps,kpe                        &
               ,its,ite, jts,jte, kts,kte                        &
               ,g3d_g(ngrid)%THSRC      & ! temp tendency
               ,g3d_g(ngrid)%RTSRC      & ! rv tendency
               ,g3d_g(ngrid)%CLSRC      & ! cloud/ice tendency
               ,g3d_g(ngrid)%cugd_ttens    &
               ,g3d_g(ngrid)%cugd_qvtens   &
                                ! forcings -  for deep/shallow
               ,cuforc_g(ngrid)%lsfth    & ! forcing for theta deep
               ,cuforc_g(ngrid)%lsfrt    & ! forcing for rv deep
               ,cuforc_sh_g(ngrid)%lsfth   & ! forcing for theta shallow
               ,cuforc_sh_g(ngrid)%lsfrt   & ! forcing for rv shallow
               ,oneGrid%oneMicVars%level                      &
               ,oneGrid%oneMicroFields%rcp         & ! liquid water
               ,oneGrid%oneMicroFields%rrp         & ! pristine
               ,oneGrid%oneMicroFields%rpp         &
               ,oneGrid%oneMicroFields%rsp         &
               ,oneGrid%oneMicroFields%rap         &
               ,oneGrid%oneMicroFields%rgp         &
               ,oneGrid%oneMicroFields%rhp         &
               ,aot500                     & ! aot at 500nm
               )
          !
          !- exchange border information for parallel run
          if( g3d_spread == 1 .or. g3d_smoothh == 1) then
             call PostSendRecvMsgs(oneGrid%SendG3D, oneGrid%RecvG3D)
             call WaitSendRecvMsgs    (oneGrid%SendG3D, oneGrid%RecvG3D)
          endif
          !
          !
          !- call routine to do the lateral spread, smooths and limiters/fixers
          call conv_grell_spread3d_brams(mzp,mxp,myp,ia,iz,ja,jz,dtlt,&
               oneGrid%oneMicVars%level,cugd_avedx&
               ,XL    &
               ,CP    &
               ,G    &
               ,rm                           &
               ,p00                          &
               ,cpor                         &
               ,cuparm_g(ngrid)%CONPRR       &!preci rate
               ,oneGrid%oneBasicFields%theta         &
               ,oneGrid%oneBasicFields%thp           &
               ,oneGrid%oneBasicFields%pp            &
               ,oneGrid%oneBasicFields%pi0           &
               ,oneGrid%oneBasicFields%rv            &
               ,tend%PT                      &
               ,oneGrid%oneMicroFields%rcp           & ! liquid water
               ,oneGrid%oneMicroFields%rrp           & ! pristine
               ,oneGrid%oneMicroFields%rpp           &
               ,oneGrid%oneMicroFields%rsp           &
               ,oneGrid%oneMicroFields%rap           &
               ,oneGrid%oneMicroFields%rgp           &
               ,oneGrid%oneMicroFields%rhp           &
                                !      
               ,g3d_g(ngrid)%THSRC           & ! temp tendency
               ,g3d_g(ngrid)%RTSRC           & ! rv tendency
               ,g3d_g(ngrid)%CLSRC           & ! cloud/ice tendency
               ,g3d_g(ngrid)%cugd_ttens      &
               ,g3d_g(ngrid)%cugd_qvtens     &
               ,g3d_ens_g(apr_gr,ngrid)%apr  &
               ,g3d_ens_g(apr_w ,ngrid)%apr  &
               ,g3d_ens_g(apr_mc,ngrid)%apr  &
               ,g3d_ens_g(apr_st,ngrid)%apr  &
               ,g3d_ens_g(apr_as,ngrid)%apr  )

          !-------------------------------------------------------------
       elseif(iinqparm==5) then  ! GF 2014 scheme
          !
          !- no lateral spreading
          cugd_avedx=1

          if(oneGrid%oneNamelistFile%ilwrtyp==4 .or. oneGrid%oneNamelistFile%iswrtyp==4) then
             aot500(:,:)=carma(ngrid)%aot(:,:,11)
          else
             aot500(:,:)=0.0
          endif

          if(iinshcu == 3) ishallow_g3=1

          call GFDRV( CCATT                        &
               ,mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
               ,mynum,i0,j0,time,mzp,mxp,myp &
               ,dtlt                           & !
               ,grid_length                    & !
               ,autoconv                       & ! Const
               ,aerovap                        & ! Const
               ,oneGrid%oneBasicFields%dn0             & !3d ok
               ,cuparm_g(ngrid)%CONPRR         & !2d ok
               ,oneGrid%oneBasicFields%up              & !3d ok
               ,oneGrid%oneBasicFields%vp              & !3d ok
               ,oneGrid%oneBasicFields%theta           & !3d ok
               ,oneGrid%oneBasicFields%thp             & !3d ok
               ,oneGrid%oneBasicFields%pp              & !3d ok
               ,oneGrid%oneBasicFields%pi0             & !3d ok
               ,oneGrid%oneBasicFields%wp              & !3d ok
               ,oneGrid%oneBasicFields%rv              & !3d ok
               ,oneGrid%oneBasicFields%rtp             & !3d ok
               ,grid_g(ngrid)%rtgt             & !2d ok
               ,tend%pt                        & !3d !*** borda
               ,xl                             & ! Const
               ,cp                             & ! Const
               ,g                              & ! Const
               ,rm                             & ! Const
               ,p00                            & ! Const
               ,cpor                           & ! Const
               ,rgas                           & ! Const
               ,zmn(:,ngrid)                   & !
               ,ztn(:,ngrid)                   & !
               ,g3d_ens_g(apr_gr,ngrid)%apr    & !2d ok
               ,g3d_ens_g(apr_w ,ngrid)%apr    & !2d ok
               ,g3d_ens_g(apr_mc,ngrid)%apr    & !2d ok
               ,g3d_ens_g(apr_st,ngrid)%apr    & !2d ok
               ,g3d_ens_g(apr_as,ngrid)%apr    & !2d ok
               ,g3d_g(ngrid)%xmb_deep          & !2d ok
               ,g3d_g(ngrid)%err_deep          & !2d ok
               ,g3d_g(ngrid)%xmb_shallow       & !2d ok
               ,g3d_ens_g(apr_gr,ngrid)%weight & !2d ok
               ,g3d_ens_g(apr_w ,ngrid)%weight & !2d ok
               ,g3d_ens_g(apr_mc,ngrid)%weight & !2d ok
               ,g3d_ens_g(apr_st,ngrid)%weight & !2d ok
               ,g3d_ens_g(apr_as,ngrid)%weight & !2d ok
               ,training                       & !
               ,grid_g(ngrid)%topt             & !2d ok
               ,leaf_g(ngrid)%patch_area       & !3d *** Borda
               ,npatch                         & !
               ,oneGrid%oneRadiateFields%rshort        & !2d ok
               ,cugd_avedx                     & !
               ,imomentum                      & !
               ,ensdim_g3d                     & !
               ,maxiens                        & !
               ,maxens_g3d                     & !
               ,maxens2_g3d                    & !
               ,maxens3_g3d                    & !
               ,icoic                          & !
               ,ishallow_g3                    & !
               ,ids,ide, jds,jde, kds,kde      & !
               ,ims,ime, jms,jme, kms,kme      & !
               ,ips,ipe, jps,jpe, kps,kpe      & !
               ,its,ite, jts,jte, kts,kte      & !
               ,g3d_g(ngrid)%THSRC             & !3d ok ! temp tendency
               ,g3d_g(ngrid)%RTSRC             & !3d ok ! rv tendency
               ,g3d_g(ngrid)%CLSRC             & !3d ok ! cloud/ice tendency
               ,g3d_g(ngrid)%cugd_ttens        & !3d ok
               ,g3d_g(ngrid)%cugd_qvtens       & !3d ok
               ,cuforc_g(ngrid)%lsfth        & !3d *** borda forcing for theta deep
               ,cuforc_g(ngrid)%lsfrt        & !3d *** borda forcing for rv deep
               ,cuforc_sh_g(ngrid)%lsfth       & !3d *** borda forcing for theta shallow
               ,cuforc_sh_g(ngrid)%lsfrt       & !3d *** borda forcing for rv shallow
               ,oneGrid%oneMicVars%level   &
               ,oneGrid%oneMicroFields%rcp             & !3d ok ! liquid water
               ,oneGrid%oneMicroFields%rrp             & !3d ok ! pristine
               ,oneGrid%oneMicroFields%rpp             & !3d ok
               ,oneGrid%oneMicroFields%rsp             & !3d ok
               ,oneGrid%oneMicroFields%rgp             & !3d ok
               ,aot500                         &! aot at 500nm
               ,oneGrid%oneTurbFields%sflux_r          & !2d *** borda
               ,oneGrid%oneTurbFields%sflux_t          & !2d *** borda
               ,oneGrid%oneTurbFields%tkep             & !3d ok
               ,TKMIN                          &
               ,akmin                   &
                                !- for convective transport-start
               ,ierr4d                         & !4d *** borda
               ,jmin4d                         & !4d ok
               ,kdet4d                         & !4d ok
               ,k224d                          & !4d ok
               ,kbcon4d                        & !4d ok
               ,ktop4d                         & !4d ok
               ,kpbl4d                         & !4d ok
               ,kstabi4d                       & !4d ok
               ,kstabm4d                       & !4d ok
               ,xmb4d                          & !4d ok
               ,edt4d                          & !4d ok
               ,pwav4d                         & !4d ok
               ,pcup5d                         & !5d ok
               ,up_massentr5d                  & !5d ok
               ,up_massdetr5d                  & !5d ok
               ,dd_massentr5d                  & !5d ok
               ,dd_massdetr5d                  & !5d ok
               ,zup5d                          & !5d *** ERRO !!!!!
               ,zdn5d                          & !5d ok
               ,prup5d                         & !5d ok
               ,prdn5d                         & !5d ok
               ,clwup5d                        & !5d ok
               ,tup5d                          & !5d ok
                                !- for convective transport- end
               )

          !.. call dumpVarAllLatLonK(g3d_g(ngrid)%THSRC,'THSRC',938,0,0,1,mxp,1,myp,1,mzp,0.,660.,header)
          !.. call dumpVarAllLatLonK(g3d_g(ngrid)%RTSRC,'RTSRC',938,0,0,1,mxp,1,myp,1,mzp,0.,660.,header)
          !.. call dumpVarAllLatLonK(g3d_g(ngrid)%CLSRC,'CLSRC',938,0,0,1,mxp,1,myp,1,mzp,0.,660.,header)

          !-------------------------------------------------------------
       elseif(iinqparm==6) then  ! GF 2015 scheme
          !
          !- no lateral spreading
          cugd_avedx=1

          if(oneGrid%oneNamelistFile%ilwrtyp==4 .or. oneGrid%oneNamelistFile%iswrtyp==4) then
             aot500(:,:)=carma(ngrid)%aot(:,:,11)
          else
             aot500(:,:)=0.0
          endif

          if(isfcl == 5) then
             temp2m(:,:) = oneGrid%oneJulesFields%t2mj(:,:)
          else
             temp2m(:,:) =0.5*( oneGrid%oneBasicFields%theta(1,:,:)* &
                  (oneGrid%oneBasicFields%pp(1,:,:)+oneGrid%oneBasicFields%pi0(1,:,:))/cp + &
                  oneGrid%oneBasicFields%theta(2,:,:)*&
                  (oneGrid%oneBasicFields%pp(2,:,:)+oneGrid%oneBasicFields%pi0(2,:,:))/cp )
          endif

          if(iinshcu == 3) ishallow_g3=1

          call GFDRV2(mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
               ,mynum,i0,j0,time,mzp,mxp,myp &
               ,dtlt             & !
               ,grid_length                  & !
               ,autoconv                     & !
               ,aerovap                      & !
               ,oneGrid%oneBasicFields%dn0           & !
               ,cuparm_g(ngrid)%CONPRR       & !
               ,oneGrid%oneBasicFields%up            & !
               ,oneGrid%oneBasicFields%vp            & !
               ,oneGrid%oneBasicFields%theta         & !
               ,oneGrid%oneBasicFields%thp           & !
               ,oneGrid%oneBasicFields%pp            & !
               ,oneGrid%oneBasicFields%pi0           & !
               ,oneGrid%oneBasicFields%wp            & !
               ,oneGrid%oneBasicFields%rv            & !
               ,oneGrid%oneBasicFields%rtp           & !
               ,grid_g(ngrid)%RTGT           & !
               ,tend%PT                      & !
               ,XL                    & !
               ,CP                    & !
               ,G                    & !
               ,rm                           &
               ,p00                          &
               ,cpor                         & !
               ,rgas                         & !
               ,zmn(:,ngrid)                 & !
               ,ztn(:,ngrid)                 & !
               ,g3d_ens_g(apr_gr,ngrid)%apr  &
               ,g3d_ens_g(apr_w ,ngrid)%apr  &
               ,g3d_ens_g(apr_mc,ngrid)%apr  &
               ,g3d_ens_g(apr_st,ngrid)%apr  &
               ,g3d_ens_g(apr_as,ngrid)%apr  &
               ,g3d_g(ngrid)%xmb_deep        &
               ,g3d_g(ngrid)%err_deep        &
               ,g3d_g(ngrid)%xmb_shallow     &
               ,g3d_ens_g(apr_gr,ngrid)%weight &
               ,g3d_ens_g(apr_w ,ngrid)%weight &
               ,g3d_ens_g(apr_mc,ngrid)%weight &
               ,g3d_ens_g(apr_st,ngrid)%weight &
               ,g3d_ens_g(apr_as,ngrid)%weight &
               ,training &
               ,grid_g(ngrid)%topt              &
               ,leaf_g(ngrid)%patch_area        &
               ,npatch                          &
               ,oneGrid%oneRadiateFields%rshort         &
               ,cugd_avedx        &
               ,imomentum            &
               ,ensdim_g3d,maxiens,maxens_g3d,maxens2_g3d,maxens3_g3d,icoic      &
               ,ishallow_g3                                      &
               ,ids,ide, jds,jde, kds,kde                        &
               ,ims,ime, jms,jme, kms,kme                        &
               ,ips,ipe, jps,jpe, kps,kpe                        &
               ,its,ite, jts,jte, kts,kte                        &
               ,g3d_g(ngrid)%THSRC      & ! temp tendency
               ,g3d_g(ngrid)%RTSRC      & ! rv tendency
               ,g3d_g(ngrid)%CLSRC      & ! cloud/ice tendency
               ,g3d_g(ngrid)%USRC       & ! U tendency
               ,g3d_g(ngrid)%VSRC       & ! V tendency
               ,g3d_g(ngrid)%MUP        & ! updraft mass flux
               ,cuforc_g(ngrid)%lsfth   & ! forcing for theta deep
               ,cuforc_g(ngrid)%lsfrt   & ! forcing for rv deep
               ,cuforc_sh_g(ngrid)%lsfth   & ! forcing for theta shallow
               ,cuforc_sh_g(ngrid)%lsfrt  & ! forcing for rv shallow
               ,oneGrid%oneMicVars%level                     &
               ,oneGrid%oneMicroFields%rcp        & ! liquid water
               ,aot500                    &! aot at 500nm
               ,temp2m                    &! aot at 500nm
               ,oneGrid%oneTurbFields%sflux_r     &
               ,oneGrid%oneTurbFields%sflux_t     &
               ,oneGrid%oneTurbFields%tkep        &
               ,TKMIN                     &
               ,akmin              &
               ,do_cupar_mcphys_coupling  &
                                !- for convective transport-start
               ,ierr4d    &
               ,jmin4d       &
               ,kdet4d       &
               ,k224d           &
               ,kbcon4d      &
               ,ktop4d       &
               ,kpbl4d       &
               ,kstabi4d     &
               ,kstabm4d     &
               ,xmb4d        &
               ,edt4d        &
               ,pwav4d        &
               ,pcup5d       &
               ,up_massentr5d     &
               ,up_massdetr5d     &
               ,dd_massentr5d     &
               ,dd_massdetr5d     &
               ,zup5d        &
               ,zdn5d        &
               ,prup5d       &
               ,prdn5d       &
               ,clwup5d      &
               ,tup5d        &
                                !- for convective transport- end
               )

          !--------------------------------------------------------------------------------------------------------------------------

       elseif(iinqparm==8) then  ! GF 2021 scheme

          if(read_GF_ConvPar_nml) then
             !-- read the GF namelist
             call GF_convpar_init(mynum)
             read_GF_ConvPar_nml = .false.
          endif

          !--- these arrays must be set to zero every timestep.
          ierr4d_tmp        = 0.0
          jmin4d_tmp        = 0.0
          klcl4d_tmp        = 0.0
          k224d_tmp       = 0.0
          kbcon4d_tmp       = 0.0
          ktop4d_tmp        = 0.0
          kstabi4d_tmp      = 0.0
          kstabm4d_tmp      = 0.0
          cprr4d_tmp        = 0.0
          xmb4d_tmp       = 0.0
          edt4d_tmp       = 0.0
          pwav4d_tmp        = 0.0
          sigma4d_tmp       = 0.0
          pcup5d_tmp        = 0.0
          up_massentr5d_tmp = 0.0
          up_massdetr5d_tmp = 0.0
          dd_massentr5d_tmp = 0.0
          dd_massdetr5d_tmp = 0.0
          zup5d_tmp       = 0.0
          zdn5d_tmp       = 0.0
          prup5d_tmp        = 0.0
          prdn5d_tmp        = 0.0
          clwup5d_tmp       = 0.0
          tup5d_tmp       = 0.0
          conv_cld_fr5d_tmp = 0.0
          CNV_FRC (:,:)     = 0.0
          TRACER  (:,:,:,:) = 0.0 
          SRC_CHEM(:,:,:,:) = 0.0 

          call set_index_loops( ims,ime, jms,jme, kms,kme,    &
               its,ite, jts,jte, kts,kte,    &
               mxp,myp,mzp                   )

          FSCAV_INT(:)   = 0.1

          do k=1,mzp
             flip   (k)   = k
          enddo
          if( idiffk /= 2 .and. idiffk /= 3) then 
             if(idiffk == 7 ) then          
                kpbl (:,:) = nint(oneGrid%oneTurbFields%kpbl(:,:))
             else
                do j=1,myp
                   do i=1,mxp
                      call get_zi_gf2018(mzp,tkmin,oneGrid%oneTurbFields%tkep(:,i,j),zmn(:,ngrid) &
                           ,grid_g(ngrid)%rtgt(i,j)                          &
                           ,grid_g(ngrid)%topt(i,j),kpbl(i,j) )
                      kpbl (i,j) = max(1,min(kpbl (i,j),mzp-1))
                   enddo
                enddo
             end if
          else
             kpbl = 5  ! check later (introduce better formulation for Zi )
          endif
          !
          do j=1,myp
             do i=1,mxp
                do_this_column(i,j)=0

                do k=1,mzp-1
                   kr=k+1
                   zm3d   (k,i,j) = zmn(kr,ngrid)*grid_g(ngrid)%rtgt(i,j) !m - height above local terrain
                   zt3d   (k,i,j) = ztn(kr,ngrid)*grid_g(ngrid)%rtgt(i,j) !m
                   dm3d   (k,i,j) = oneGrid%oneBasicFields%dn0  (kr,i,j) !kg/m3
                   rvap   (k,i,j) = oneGrid%oneBasicFields%rv   (kr,i,j) !kg/kg

                   theta2temp     = (oneGrid%oneBasicFields%pp(kr,i,j)+oneGrid%oneBasicFields%pi0(kr,i,j))/cp   !K
                   temp   (k,i,j) = oneGrid%oneBasicFields%theta(kr,i,j)* theta2temp
                   press  (k,i,j) = ((oneGrid%oneBasicFields%pp(kr,i,j)+oneGrid%oneBasicFields%pi0(kr,i,j))/cp)**cpor*p00 !Pa
                   up     (k,i,j) = oneGrid%oneBasicFields%up(kr,i,j) !m/s
                   vp     (k,i,j) = oneGrid%oneBasicFields%vp(kr,i,j) !m/s
                   wp     (k,i,j) = oneGrid%oneBasicFields%wp(kr,i,j)*(-g*oneGrid%oneBasicFields%dn0(kr,i,j)) ! omega Pa/s

                   gsf_t (k,i,j) = (cuforc_g   (ngrid)%lsfth(kr,i,j) + oneGrid%oneRadiateFields%fthrd(kr,i,j))* theta2temp ! Adv+Rad, K/s
                   gsf_q (k,i,j) =  cuforc_g   (ngrid)%lsfrt(kr,i,j)              !kg/kg/s  Adv only
                   sgsf_t (k,i,j) =  cuforc_sh_g(ngrid)%lsfth(kr,i,j) * theta2temp !K/s     PBL only 
                   sgsf_q (k,i,j) =  cuforc_sh_g(ngrid)%lsfrt(kr,i,j)              !kg/kg/s PBL only 
                   advf_t (k,i,j) =  cuforc_g   (ngrid)%lsfth(kr,i,j) * theta2temp ! advection only, see 'prepare_lsl' routine
                   !falta---
                   buoy_exc(k,i,j)= 0.
                   !falta---    
                enddo
             enddo
          enddo

          if(APPLY_SUB_MP == 1) then
             do j=1,myp
                do i=1,mxp
                   do k=1,mzp-1
                      kr=k+1
                      mp_ice   (:,k,i,j) = 0. ! microg%ice(kr,i,j) in the future includes ice mix ratio 
                      mp_liq   (:,k,i,j) = 0. ! microg%liq(kr,i,j) in the future includes liq ratio 
                      mp_cf    (:,k,i,j) = 0. ! cloud fraction  
                   enddo
                enddo
             enddo
          endif

          do j=1,myp
             do i=1,mxp
                sfc_press(i,j) = 0.5*( ((oneGrid%oneBasicFields%pp(1,i,j)+oneGrid%oneBasicFields%pi0(1,i,j))/cp)**cpor*p00 +  &
                     ((oneGrid%oneBasicFields%pp(2,i,j)+oneGrid%oneBasicFields%pi0(2,i,j))/cp)**cpor*p00 ) !Pa

                xland(i,j) = leaf_g(ngrid)%patch_area(i,j,1) ! water = 1, land < 1
                lons (i,j) = grid_g(ngrid)%glon(i,j)
                lats (i,j) = grid_g(ngrid)%glat(i,j)
                !dx2d (i,j) = sqrt(deltaxn(ngrid)*deltayn(ngrid))
                dx2d (i,j) = sqrt(1./(grid_g(ngrid)%dxt(i,j)*grid_g(ngrid)%dyt(i,j)))
                col_sat       (i,j) = 0.
                stochastic_sig(i,j) = 1.
             enddo
          enddo

          if(oneGrid%oneNamelistFile%ilwrtyp==4 .or. oneGrid%oneNamelistFile%iswrtyp==4) then
             aot500(:,:)=carma(ngrid)%aot(:,:,11)
          else
             aot500(:,:)=0.0
          endif

          if(isfcl == 5 .and.  time > dtlt ) then
             temp2m(:,:) = oneGrid%oneJulesFields%t2mj(:,:) !K
          else
             temp2m(:,:) =0.5*(oneGrid%oneBasicFields%theta(1,:,:)* &
                  (oneGrid%oneBasicFields%pp(1,:,:)+oneGrid%oneBasicFields%pi0(1,:,:))/cp + &
                  oneGrid%oneBasicFields%theta(2,:,:)*&
                  (oneGrid%oneBasicFields%pp(2,:,:)+oneGrid%oneBasicFields%pi0(2,:,:))/cp ) !Kelvin
          endif

          !- call the driver routine to apply the parameterization
          call GF_GEOS5_DRV(mxp,myp,mzp,mtp ,nmp        &
               ,ims,ime, jms,jme, kms,kme   &
               ,its,ite, jts,jte, kts,kte   &
               ,flip        &
               ,fscav_int   &
               ,mynum       &
               ,dtlt        &
               ,dx2d        &
               ,stochastic_sig &
               ,zm3d        &
               ,zt3d        &
               ,dm3d        &
                                !--- sfc inputs
               ,lons        &
               ,lats        &
               ,aot500      &
               ,temp2m      &
               ,oneGrid%oneTurbFields%sflux_r &
               ,oneGrid%oneTurbFields%sflux_t &
               ,grid_g(ngrid)%topt    &
               ,xland                 &
               ,sfc_press   &
               ,kpbl        &
                                !--- atmos state
               ,col_sat     &
               ,up     &
               ,vp     &
               ,wp     &
               ,temp   &
               ,press  &
               ,rvap   &
               ,mp_ice &
               ,mp_liq &
               ,mp_cf  &
               ,rvap   &
                                !--- atmos composition state
               ,TRACER   & !- note: uses GEOS-5 data structure
                                !---- forcings---
               ,buoy_exc &
               , gsf_t   & ! forcing for theta adv+rad
               , gsf_q   & ! forcing for rv    adv
               ,advf_t   &
               ,sgsf_t   & ! forcing for theta pbl
               ,sgsf_q   & ! forcing for rv    pbl
                                !---- output ----
               ,cuparm_g(ngrid)%CONPRR  &
               ,LIGHTN_DENS             &
               ,g3d_g(ngrid)%THSRC      & ! temp tendency
               ,g3d_g(ngrid)%RTSRC      & ! rv tendency
               ,g3d_g(ngrid)%CLSRC      & ! cloud/ice  mass   mix ratio tendency
               ,g3d_g(ngrid)%NLSRC      & ! cloud drop number mix ratio tendency
               ,g3d_g(ngrid)%NISRC      & ! ice        number mix ratio tendency
               ,g3d_g(ngrid)%USRC       & ! U tendency
               ,g3d_g(ngrid)%VSRC       & ! V tendency
               ,SUB_MPQI    & 
               ,SUB_MPQL    & 
               ,SUB_MPCF    & 
               ,SRC_BUOY    &
               ,SRC_CHEM    & ! tracer tendency
               ,REVSU_GF    &
               ,PRFIL_GF    & 
                                !
               ,do_this_column    &
               ,ierr4d_tmp        & 
               ,jmin4d_tmp        &
               ,klcl4d_tmp        &
               ,k224d_tmp         &
               ,kbcon4d_tmp       &
               ,ktop4d_tmp        &
               ,kstabi4d_tmp      &
               ,kstabm4d_tmp      &
               ,cprr4d_tmp        &
               ,xmb4d_tmp         &
               ,edt4d_tmp         &
               ,pwav4d_tmp        &
               ,sigma4d_tmp       &
               ,pcup5d_tmp        &
               ,up_massentr5d_tmp &
               ,up_massdetr5d_tmp &
               ,dd_massentr5d_tmp &
               ,dd_massdetr5d_tmp &
               ,zup5d_tmp         &
               ,zdn5d_tmp         &
               ,prup5d_tmp        &
               ,prdn5d_tmp        &
               ,clwup5d_tmp       &
               ,tup5d_tmp         &
               ,conv_cld_fr5d_tmp &
                                !-- for debug/diagnostic
               ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC  &
               ,VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF &
               )
          !
          !-- outputs ....

          !if( icumulus_gf(deep) == 1) then 
          ! g3d_ens_g(1,ngrid)%accapr(:,:)=AA0(:,:) !Tpert_h(2, :,:)
          ! g3d_ens_g(2,ngrid)%accapr(:,:)=AA1(:,:) !Tpert_h(5,:,:)
          ! g3d_ens_g(3,ngrid)%accapr(:,:)=AA2(:,:) !Tpert_h(10,:,:)
          ! g3d_ens_g(4,ngrid)%accapr(:,:)=AA3(:,:) !Tpert_v(2,:,:)
          !endif

          !-- saving the precip of each mode in the array g3d_ens_g(1,ngrid)%accapr(
          if( icumulus_gf(deep) == 1)  g3d_ens_g(1,ngrid)%apr(:,:)= cprr4d_tmp(:,:,deep)  
          if( icumulus_gf(mid)  == 1)  g3d_ens_g(3,ngrid)%apr(:,:)= cprr4d_tmp(:,:,mid)  
          if( icumulus_gf(shal) == 1)  g3d_ens_g(2,ngrid)%apr(:,:)= cprr4d_tmp(:,:,shal)  


          if( icumulus_gf(deep) == 1) then 
             ! updraft mass flux at cloud base
             g3d_g(ngrid)%xmb_deep   (:,:) = xmb4d_tmp       (:,:,deep)  
             ! downdraft mass flux at cloud base
             g3d_g(ngrid)%xmb_deep_dd(:,:) = xmb4d_tmp       (:,:,deep) * edt4d_tmp(:,:,deep)
             g3d_g(ngrid)%err_deep   (:,:) = float(ierr4d_tmp(:,:,deep))
          endif
          if( icumulus_gf(shal) == 1) then 
             g3d_g(ngrid)%xmb_shallow(:,:) = xmb4d_tmp       (:,:,shal)
          endif

          g3d_g(ngrid)%THSRC=eoshift(g3d_g(ngrid)%THSRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%THSRC(1,:,:), DIM=1)
          g3d_g(ngrid)%RTSRC=eoshift(g3d_g(ngrid)%RTSRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%RTSRC(1,:,:), DIM=1)
          g3d_g(ngrid)%CLSRC=eoshift(g3d_g(ngrid)%CLSRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%CLSRC(1,:,:), DIM=1)
          g3d_g(ngrid)%USRC =eoshift(g3d_g(ngrid)%USRC,  SHIFT=-1, BOUNDARY=g3d_g(ngrid)%USRC (1,:,:), DIM=1)
          g3d_g(ngrid)%VSRC =eoshift(g3d_g(ngrid)%VSRC,  SHIFT=-1, BOUNDARY=g3d_g(ngrid)%VSRC (1,:,:), DIM=1)


          if(liq_ice_number_conc > 0) then
             g3d_g(ngrid)%NLSRC=eoshift(g3d_g(ngrid)%NLSRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%NLSRC(1,:,:), DIM=1)
             g3d_g(ngrid)%NISRC=eoshift(g3d_g(ngrid)%NISRC, SHIFT=-1, BOUNDARY=g3d_g(ngrid)%NISRC(1,:,:), DIM=1)
          endif

          !-- converting Dtemp/Dt to Dtheta/ Dt (temp = cp * theta/exner function) 
          g3d_g(ngrid)%THSRC = g3d_g(ngrid)%THSRC * cp / (oneGrid%oneBasicFields%pp + oneGrid%oneBasicFields%pi0)

          if( icumulus_gf(deep) == 1) then 
             do j=1,myp
                do i=1,mxp
                   temp(:,i,j)= zup5d_tmp(i,j,:,deep) ! zup already includes the XMB
                enddo
             enddo
             g3d_g(ngrid)%MUP=eoshift(temp,  SHIFT=-1, BOUNDARY=temp (1,:,:), DIM=1)
             do j=1,myp
                do i=1,mxp
                   temp(:,i,j)= zdn5d_tmp(i,j,:,deep)* edt4d_tmp(i,j,deep) ! zdo already includes the XMB
                enddo
             enddo
             g3d_g(ngrid)%MDD=eoshift(temp,  SHIFT=-1, BOUNDARY=temp (1,:,:), DIM=1)
          endif

          if( icumulus_gf(shal) == 1) then 
             do j=1,myp
                do i=1,mxp
                   temp(:,i,j)= zup5d_tmp(i,j,:,shal) ! zup already includes the XMB
                enddo
             enddo
             g3d_g(ngrid)%MUPSH=eoshift(temp,  SHIFT=-1, BOUNDARY=temp (1,:,:), DIM=1)
          endif

          !--- output for RRTM/CARMA and convective transport
          if((oneGrid%oneNamelistFile%ilwrtyp>0 .or. oneGrid%oneNamelistFile%iswrtyp>0) .or. chemistry >= 0) then
             do j=1,myp
                do i=1,mxp
                   if(do_this_column(i,j)==0) cycle
                   do n=1,maxiens
                      if( icumulus_gf(n) /= 1) cycle 
                      xmb4d   (i,j,n,ngrid) = xmb4d_tmp   (i,j,n)
                      ierr4d  (i,j,n,ngrid) = ierr4d_tmp  (i,j,n)
                   enddo

                   do k=1,mzp
                      do n=1,maxiens
                         if( icumulus_gf(n) /= 1) cycle 
                         zup5d        (k,i,j,n,ngrid) =  zup5d_tmp        (i,j,k,n)
                         clwup5d      (k,i,j,n,ngrid) =  clwup5d_tmp      (i,j,k,n)
                         up_massdetr5d(k,i,j,n,ngrid) =  up_massdetr5d_tmp(i,j,k,n)
                      enddo
                   enddo
                   if( chemistry >= 0) then ! - for convective transport only
                      do n=1,maxiens
                         if( icumulus_gf(n) /= 1) cycle 
                         jmin4d  (i,j,n,ngrid) = jmin4d_tmp  (i,j,n)
                         k224d   (i,j,n,ngrid) = k224d_tmp   (i,j,n)
                         kbcon4d (i,j,n,ngrid) = kbcon4d_tmp (i,j,n)
                         ktop4d  (i,j,n,ngrid) = ktop4d_tmp  (i,j,n)
                         kstabi4d(i,j,n,ngrid) = kstabi4d_tmp(i,j,n)
                         kstabm4d(i,j,n,ngrid) = kstabm4d_tmp(i,j,n)
                         edt4d   (i,j,n,ngrid) = edt4d_tmp   (i,j,n)
                         pwav4d  (i,j,n,ngrid) = pwav4d_tmp  (i,j,n)
                      enddo

                      do k=1,mzp
                         do n=1,maxiens
                            if( icumulus_gf(n) /= 1) cycle 
                            zdn5d        (k,i,j,n,ngrid) =  zdn5d_tmp        (i,j,k,n)
                            up_massentr5d(k,i,j,n,ngrid)=  up_massentr5d_tmp(i,j,k,n)
                            dd_massentr5d(k,i,j,n,ngrid)=  dd_massentr5d_tmp(i,j,k,n)
                            dd_massdetr5d(k,i,j,n,ngrid)=  dd_massdetr5d_tmp(i,j,k,n)
                            pcup5d       (k,i,j,n,ngrid)=  pcup5d_tmp       (i,j,k,n)
                            prup5d       (k,i,j,n,ngrid)=  prup5d_tmp       (i,j,k,n)
                            prdn5d       (k,i,j,n,ngrid)=  prdn5d_tmp       (i,j,k,n)
                            tup5d        (k,i,j,n,ngrid)=  tup5d_tmp        (i,j,k,n)
                         enddo
                      enddo
                   endif
                enddo
             enddo
          endif

          !
          !for checking
          !do j=1,myp
          !  do i=1,mxp
          ! if(do_this_column(i,j)==0) cycle
          !
          ! call moveup(mzp,g3d_g(ngrid)%THSRC (:,i,j))
          ! call moveup(mzp,g3d_g(ngrid)%RTSRC (:,i,j))
          ! call moveup(mzp,g3d_g(ngrid)%CLSRC (:,i,j))
          ! call moveup(mzp,g3d_g(ngrid)%USRC  (:,i,j))
          !  call moveup(mzp,g3d_g(ngrid)%VSRC  (:,i,j))
          !enddo;enddo
          !--------------------------------------------------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------------------------------------------
       endif

       !--- filling the output tendencies for the level k =1 and k=mzp
       g3d_g(ngrid)%THSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%THSRC(2    ,1:mxp,1:myp)
       g3d_g(ngrid)%RTSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%RTSRC(2    ,1:mxp,1:myp)
       g3d_g(ngrid)%CLSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%CLSRC(2    ,1:mxp,1:myp)
       g3d_g(ngrid)%THSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%THSRC(mzp-1,1:mxp,1:myp)
       g3d_g(ngrid)%RTSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%RTSRC(mzp-1,1:mxp,1:myp)
       g3d_g(ngrid)%CLSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%CLSRC(mzp-1,1:mxp,1:myp)
       if(imomentum==1 .and. oneGrid%oneNamelistFile%nnqparm(ngrid) >= 4) then
          g3d_g(ngrid)%USRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%USRC (2    ,1:mxp,1:myp)
          g3d_g(ngrid)%VSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%VSRC (2    ,1:mxp,1:myp)
          g3d_g(ngrid)%USRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%USRC (mzp-1,1:mxp,1:myp)
          g3d_g(ngrid)%VSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%VSRC (mzp-1,1:mxp,1:myp)
       endif
       if(liq_ice_number_conc>0) then
          g3d_g(ngrid)%NLSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%NLSRC(mzp-1,1:mxp,1:myp)
          g3d_g(ngrid)%NISRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%NISRC(mzp-1,1:mxp,1:myp)
       endif
       !-------------------------------------------------------------
    endif! 002
    !-------------------------------------------------------------
    ! stores precipitation rate for each closure, only for output/training
    if(oneGrid%oneNamelistFile%nnqparm(ngrid) == 3) then
       if (training > 0) then
          do i=1,train_dim
             call update(mxp*myp, g3d_ens_g(i,ngrid)%accapr,g3d_ens_g(i,ngrid)%apr,dtlt)
          enddo
       endif
    endif
    !--- for output only 
    if(oneGrid%oneNamelistFile%nnqparm(ngrid) == 8) then
       do i=1,train_dim
          call update(mxp*myp, g3d_ens_g(i,ngrid)%accapr,g3d_ens_g(i,ngrid)%apr,dtlt)
       enddo
    endif
    !----------------------------------------------------------

    call update(mxp*myp, cuparm_g(ngrid)%aconpr   ,cuparm_g(ngrid)%conprr   ,dtlt)

    call accum(int(mxp*myp*mzp,i8), tend%tht, g3d_g(ngrid)%thsrc)
    call accum(int(mxp*myp*mzp,i8), tend%rtt, g3d_g(ngrid)%rtsrc)

    if(imomentum == 1 .and. oneGrid%oneNamelistFile%nnqparm(ngrid) >= 4) then
       call accum(int(mxp*myp*mzp,i8), tend%ut, g3d_g(ngrid)%usrc)
       call accum(int(mxp*myp*mzp,i8), tend%vt, g3d_g(ngrid)%vsrc)
    endif

    if(do_cupar_mcphys_coupling) then
       call cupar2mcphysics(mzp,mxp,myp,ia,iz,ja,jz,ngrid,dtlt,&
            g3d_g  (ngrid)%clsrc   ,&
            oneGrid%oneBasicFields%theta   ,&
            oneGrid%oneBasicFields%pp      ,&
            oneGrid%oneBasicFields%pi0     ,&
            oneGrid%oneBasicFields%dn0, &
            oneGrid%oneMicVars,&
            oneGrid%oneMicroFields)
    else
       !if there is not direct coupling, send cloud/ice source to rtotal tendency
       call accum(int(mxp*myp*mzp,i8), tend%rtt, g3d_g(ngrid)%clsrc)

    endif
    !
    !--------- Convective Transport based on mass flux scheme -
    if (CCATT == 1 .and. iruncon == 1 .and. (iinqparm==5 .or. iinqparm==6) ) then

       if(iinqparm==5 .and. iens .ne. 1 ) &
            stop 'conv transp with GF scheme version 2014 only for deep convection'

       !- this call convective transport for deep convection
       call trans_conv_mflx_GF(1,iinqparm)

       !- if shallow convection was solved by GF version 2015, call again
       !- the convective transport routine to include the transport
       !- by shallow convection scheme.
       if(iinqparm==6 .and. iinshcu == 3)  call trans_conv_mflx_GF(2,iinqparm)

    endif

    !- this calls convective transport for deep/mid/shallow convection
    if (CCATT == 1 .and. iruncon == 1 .and. (iinqparm==8) ) then
       do n = 1, maxiens
          if(n == shal .and. iinshcu == 0) cycle
          call trans_conv_mflx_GF(n,iinqparm)
       enddo
    endif


    ! [ML------------- Stilt - BRAMS coupling  ------------------
    if (imassflx == 1 ) then

       !-srf -  mass fluxes from deep convection
       if( iinqparm==5 .or. iinqparm==6 )                                 &
            call prep_convflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz              &
            ,mgmxp,mgmyp,mgmzp,maxiens,ngrid,ngrids_cp  &
            ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d   &
            ,kstabi4d,kstabm4d,xmb4d,edt4d  &
            ,zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d  &
            ,1, oneGrid%oneNamelistFile%nnqparm)! = iens
       !-srf if shallow convection was solved by GF version 2015, call again
       !-    the convective transport routine to include the mass fluxes
       !-    from the shallow convection scheme.
       if(  iinqparm==6 .and. iinshcu == 3 )                              &
            call prep_convflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz              &
            ,mgmxp,mgmyp,mgmzp,maxiens,ngrid,ngrids_cp  &
            ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d   &
            ,kstabi4d,kstabm4d,xmb4d,edt4d  &
            ,zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d  &
            ,2, oneGrid%oneNamelistFile%nnqparm)! = iens

    endif
    ! ------------- Stilt - BRAMS coupling  ------------------ ML]
    !
  end subroutine CUPARM_GRELL3_CATT
  !
  !-------------------------------------------------------------------------------------------------
  !
  subroutine init_weights(ng,n2,n3,nnqparm)
    integer, intent(in)::ng,n2,n3,nnqparm
    integer :: it,i,j
    real sumx,hweight
    !- ordem dos pesos
    !apr_gr=001
    !apr_w =002
    !apr_mc=003
    !apr_st=004
    !apr_as=005

    if(training == 0) return

    !-- training on closures
    if(training == 1) then
       if(nnqparm==3) hweight = 0.2
       if(nnqparm==5) hweight = 0.25
       do j=1,n3
          do i=1,n2
             do it=1,train_dim

                g3d_ens_g(it,ng)%weight(i,j)=hweight
                !print*,'weights=', it,i,j, g3d_ens_g(it,ng)%weight(i,j)

                !if(it==apr_st) g3d_ens_g(it,ng)%weight(i,j)=0.175
                !if(it==apr_as) g3d_ens_g(it,ng)%weight(i,j)=0.25
                !if(it==apr_w ) g3d_ens_g(it,ng)%weight(i,j)=0.25
                !if(it==apr_mc) g3d_ens_g(it,ng)%weight(i,j)=0.25
                if(nnqparm==5) then

                   !-special treatment over the ocean
                   if(leaf_g(ng)%patch_area(i,j,1) .gt. 0.9 ) then ! water
                      if(it==apr_as) g3d_ens_g(it,ng)%weight(i,j)=0.0
                      if(it==apr_st) g3d_ens_g(it,ng)%weight(i,j)=0.425
                      if(it==apr_gr) g3d_ens_g(it,ng)%weight(i,j)=0.1667
                      if(it==apr_w ) g3d_ens_g(it,ng)%weight(i,j)=0.1667
                      if(it==apr_mc) g3d_ens_g(it,ng)%weight(i,j)=0.1667

                   else ! land

                      if(it==apr_as) g3d_ens_g(it,ng)%weight(i,j)=0.0
                      if(it==apr_st) g3d_ens_g(it,ng)%weight(i,j)=0.175
                      if(it==apr_gr) g3d_ens_g(it,ng)%weight(i,j)=0.25
                      if(it==apr_w ) g3d_ens_g(it,ng)%weight(i,j)=0.25
                      if(it==apr_mc) g3d_ens_g(it,ng)%weight(i,j)=0.25
                   endif
                endif
                !g3d_ens_g(it,ng)%weight(i,j)=float(i+j)*exp(-(float(it-2))**2)*float(i*j)

             enddo
          enddo
       enddo

       !-- training on CAPS
    elseif(training == 2) then

       do j=1,n3
          do i=1,n2

             g3d_ens_g(apr_gr,ng)%weight(i,j)=0.3333
             g3d_ens_g(apr_w ,ng)%weight(i,j)=0.3333
             g3d_ens_g(apr_mc,ng)%weight(i,j)=0.3333
             g3d_ens_g(apr_st,ng)%weight(i,j)=0.0
             g3d_ens_g(apr_as,ng)%weight(i,j)=0.0

          enddo
       enddo

    endif



    return! <<<<
    if(training == 1) then
       !- normalize a 1
       do j=1,n3
          do i=1,n2
             sumx=0.
             do it=1,train_dim
                sumx=sumx+g3d_ens_g(it,ng)%weight(i,j)
             enddo
             do it=1,train_dim
                g3d_ens_g(it,ng)%weight(i,j) = g3d_ens_g(it,ng)%weight(i,j)/sumx
             enddo
          enddo
       enddo
    endif

  end subroutine init_weights
  !-------------------------------------------------------------
  !-------------------------------------------------------------
  subroutine conv_grell_spread3d_brams(m1,m2,m3,ia,iz,ja,jz,dt, &
       level,cugd_avedx,                              &
       XLV,CP,G,r_v,p00,cpor,                         &
       conprr, theta,thetail,pp,pi0,                  &
       rv,pt,rcp,rrp,rpp,rsp,rap,rgp,rhp,             &
       RTHcuten,      &
       RQVcuten,      &
       RQCcuten,      &
       cugd_ttens,      &
       cugd_qvtens,      &
       apr_gr,      &
       apr_w,      &
       apr_mc,      &
       apr_st,      &
       apr_as      )

    integer,      intent(IN   ) :: m1,m2,m3,ia,iz,ja,jz,level,cugd_avedx
    real,         intent(IN   ) :: dt
    real,         intent(IN   ) :: XLV, R_v
    real,         intent(IN   ) :: CP,G, cpor, p00

    real, dimension(m1,m2,m3),intent(IN   ) ::     &
         theta   ,&
         thetail ,&
         pp ,&
         pi0 ,&
         pt ,&
         rv      ,rcp,rrp,rpp,rsp,rap,rgp,rhp


    real, dimension(m2,m3),intent(INOUT) ::   &
         conprr,                       &
         apr_gr,     &
         apr_w ,     &
         apr_mc,     &
         apr_st,     &
         apr_as



    real, dimension(m1,m2,m3),intent(INOUT) ::    &
         RTHcuten,    &
         RQVcuten,    &
         RQCcuten,    &
         cugd_ttens,  &
         cugd_qvtens


    ! local var
    real  ::   exner,r_sol,r_liq,fxc,tempk,dfxcdt,outt
    integer :: j,i,k,kk,jfs,jfe,ifs,ife,kts,kte,ii,jj
    integer :: cugd_spread

    real, dimension (m1,m2,m3) :: & ! orig (its-2:ite+2,kts:kte,jts-2:jte+2) ::     &
         RTHcuten_tmp, &  ! tmp RTHcuten
         RQVcuten_tmp     ! tmp RQVcuten

    real, dimension (m2,m3) :: & ! orig (its-2:ite+2,jts-2:jte+2) ::
         Qmem

    real   :: & ! orig (its-1:ite+1,jts-1:jte+1) ::
         smTT,smTQ

    real, dimension (m1) :: & ! orig (kts:kte) ::
         conv_TRASHT,conv_TRASHQ

    real :: Qmem1,Qmem2,Qmemf,Thresh

    !-initial settings
    ! g3d_smoothh=0  ! 0 or 1: do horizontal smoothing
    ! g3d_smoothv=0  ! 0 or 1: do vertical smoothing
    cugd_spread=cugd_avedx/2 ! = 0/1 => no/do spreading

    RTHcuten_tmp  = 0.0
    RQVcuten_tmp  = 0.0
    Qmem       = 1.0
    smTT       = 0.0
    smTQ       = 0.0
    conv_TRASHT= 0.0
    conv_TRASHQ= 0.0
    jfs=ja
    jfe=jz
    ifs=ia
    ife=iz
    kts=2
    kte=m1-1 !check if this correct or should be kte=m1

    !if(g3d_smoothh ==1 .or. cugd_spread > 0) then
    !  jfs=1
    !  jfe=m3
    !  ifs=1
    !  ife=m2
    !endif

    !- store input tendencies
    ! *** jm note -- for smoothing this goes out one row/column beyond tile in i and j
    do j=1,m3
       do i=1,m2
          RTHcuten_tmp(:,i,j)=RTHcuten (:,i,j)
          RQVcuten_tmp(:,i,j)=RQVcuten (:,i,j)
       enddo
    enddo



    ! ---------------- spreading   section --------------
    do j=ja,jz
       do i=ia,iz
          !
          ! for high res run, spread the subsidence
          ! this is tricky......only consider grid points where there was no rain,
          ! so cugd_tten and such are zero!
          !
          !      if do spreading
          if(cugd_spread > 0)then
             do k=kts,kte
                do jj=j-1,j+1 ! only 3 neighboors
                   do ii=i-1,i+1 ! only 3 neighboors

                      RTHcuten_tmp(k,i,j)=RTHcuten_tmp(k,i,j)     &
                           +Qmem(ii,jj)*cugd_ttens(k,ii,jj)

                      RQVcuten_tmp(k,i,j)=RQVcuten_tmp(k,i,j)     &
                           +Qmem(ii,jj)*cugd_qvtens(k,ii,jj)
                   enddo
                enddo
             enddo
             !      end spreading
             !
             !      if not spreading
          elseif(cugd_spread == 0)then
             do k=kts,kte
                RTHcuten_tmp(k,i,j)=RTHcuten_tmp(k,i,j)+cugd_ttens (k,i,j)
                RQVcuten_tmp(k,i,j)=RQVcuten_tmp(k,i,j)+cugd_qvtens(k,i,j)
             enddo
          endif
          !
       enddo  ! end i
    enddo  ! end j

    ! ----------------horizontal smoothing  section --------------

    !- if not doing horizontal smoothing, get the final tendencies
    if(g3d_smoothh == 0)then
       do j=ja,jz
          do i=ia,iz
             do k=kts,kte
                RTHcuten(k,i,j)=RTHcuten_tmp(k,i,j)
                RQVcuten(k,i,j)=RQVcuten_tmp(k,i,j)
             enddo ! enf k
          enddo  ! end j
       enddo  ! end j

       !- if doing horizontal smoothing ...
    else if(g3d_smoothh == 1)then
       do k=kts,kte
          do j=ja,jz
             do i=ia,iz

                smTT = 0.0
                smTQ = 0.0
                do jj=j-1,j+1 ! only 3 neighboors
                   do ii=i-1,i+1 ! only 3 neighboors
                      smTT = smTT +RTHcuten_tmp(k,ii,jj)
                      smTQ = smTQ +RQVcuten_tmp(k,ii,jj)

                   enddo  ! end ii
                enddo  ! end jj

                RTHcuten(k,i,j)=(3.*RTHcuten_tmp(k,i,j) + smTT)/12.
                RQVcuten(k,i,j)=(3.*RQVcuten_tmp(k,i,j) + smTQ)/12.

             enddo  ! end i
          enddo  ! end j
       enddo  ! end k

    endif  ! g3d_smoothh
    !
    ! - checking and limiting moistening/heating rates
    !
    do j=ja,jz
       do i=ia,iz
          !--- moistening section ------
          Qmemf  = 1.0
          Thresh = 1.e-20
          do k=kts,kte

             if(RQVcuten(k,i,j) < 0.0) then
                Qmem1 = rv(k,i,j)+RQVcuten(k,i,j)*dt
                if(Qmem1 < Thresh)then
                   Qmem1 = RQVcuten(k,i,j)
                   Qmem2 = (Thresh-rv(k,i,j))/dt
                   Qmemf = min(Qmemf,Qmem2/Qmem1)
                   Qmemf = max(0.,Qmemf)
                   Qmemf = min(1.,Qmemf)
                endif
             endif

          enddo  ! end k
          ! - limiting moistening
          do k=kts,kte
             RQVcuten   (k,i,j) = RQVcuten   (k,i,j)*Qmemf
             RQCcuten   (k,i,j) = RQCcuten   (k,i,j)*Qmemf
             RTHcuten   (k,i,j) = RTHcuten   (k,i,j)*Qmemf
             cugd_ttens (k,i,j) = cugd_ttens (k,i,j)*Qmemf
             cugd_qvtens(k,i,j) = cugd_qvtens(k,i,j)*Qmemf

          enddo ! end k

          !- limiting precip for consistency
          conprr (i,j) = conprr (i,j)*Qmemf
          apr_gr (i,j) = apr_gr (i,j)*Qmemf
          apr_w  (i,j) = apr_w  (i,j)*Qmemf
          apr_mc (i,j) = apr_mc (i,j)*Qmemf
          apr_st (i,j) = apr_st (i,j)*Qmemf
          apr_as (i,j) = apr_as (i,j)*Qmemf
          ! no futuro inclua tambem o limting para o fluxo de massa
          ! xmb(i,j)=xmb(i,j)*qmemf

          !--- heating section ------
          Thresh=200. ! max heating/cooling rate allowed  K/day
          !srf         Thresh=100. ! max heating/cooling rate allowed  K/day
          Qmemf=1.
          Qmem1=0.

          do k=kts,kte
             Qmem1=abs(RTHcuten(k,i,j))*86400.

             if(Qmem1 > Thresh)then
                Qmem2 = Thresh/Qmem1
                Qmemf = min(Qmemf,Qmem2)
                Qmemf = max(0.,Qmemf)
             endif
          enddo

          ! - limiting heating/cooling
          do k=kts,kte
             RTHcuten   (k,i,j) = RTHcuten   (k,i,j)*Qmemf
             RQVcuten   (k,i,j) = RQVcuten   (k,i,j)*Qmemf
             RQCcuten   (k,i,j) = RQCcuten   (k,i,j)*Qmemf
             cugd_ttens (k,i,j) = cugd_ttens (k,i,j)*Qmemf
             cugd_qvtens(k,i,j) = cugd_qvtens(k,i,j)*Qmemf
          enddo ! end k

          !- limiting precip for consistency
          conprr (i,j) = conprr (i,j)*Qmemf
          apr_gr (i,j) = apr_gr (i,j)*Qmemf
          apr_w  (i,j) = apr_w  (i,j)*Qmemf
          apr_mc (i,j) = apr_mc (i,j)*Qmemf
          apr_st (i,j) = apr_st (i,j)*Qmemf
          apr_as (i,j) = apr_as (i,j)*Qmemf
          ! no futuro inclua tambem o limting para o fluxo de massa
          ! xmb(i,j)=xmb(i,j)*qmemf

       enddo  ! end i
    enddo  ! end j
    !
    ! ---  vertical smooth ------------
    !
    if (g3d_smoothv == 1)then

       do j=ja,jz
          do i=ia,iz

             do k=kts+2,kte-2
                conv_TRASHT(k)= .25*(RTHcuten(k-1,i,j)+2.*RTHcuten(k,i,j)+RTHcuten(k+1,i,j))
                conv_TRASHQ(k)= .25*(RQVcuten(k-1,i,j)+2.*RQVcuten(k,i,j)+RQVcuten(k+1,i,j))
             enddo
             do k=kts+2,kte-2
                RTHcuten(k,i,j)=conv_TRASHT(k)
                RQVcuten(k,i,j)=conv_TRASHQ(k)
             enddo
          enddo  ! end i
       enddo  ! end j

    endif

    ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
    ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,

    !if(level <=2) then

    do j=ja,jz
       do i=ia,iz
          do k=kts,kte

             !if(RTHCUTEN (k,i,j) /= 0.0) then

             ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
             ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
             ! Exner's function = pp(k,i,j)+pi0(k,i,j)
             exner   = pp(k,i,j) + pi0(k,i,j)
             ! tendencia do theta devida a conv profunda
             RTHcuten(k,i,J) = cp/exner * RTHCUTEN(k,i,J) !- theta(k,i,j)*pt(k,i,j)/exner
             !endif
             !RQVcuten(k,i,J) = RQVCUTEN(k,i,J)+ RQCCUTEN(k,i,J)
             !RQCcuten(k,i,J) = 0.

          enddo
       enddo
    enddo


    !elseif(level > 2) then
    !
    ! do j=ja,jz; do i=ia,iz; do k=kts,kte
    !    !
    ! ! - tend na temperatura (para uso na convers�o do thetail
    ! outt=RTHCUTEN (k,i,j)
    ! ! Exner's function = pp(k,i,j)+pi0(k,i,j)
    ! exner= pp(k,i,j) + pi0(k,i,j)
    ! if(outt /= 0.0 ) then
    !   !
    !   ! converte tend da temperatura (outt) em tend de theta (outtem)
    !   ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
    !
    !   ! tendencia do theta  devida a conv profunda
    !   RTHCUTEN (k,i,j) = cp/exner * RTHCUTEN(k,i,j) - theta(k,i,j)*pt(k,i,j)/exner
    !
    ! endif
    !
    ! ! tendencia do theta_il devida a conv profunda
    ! r_liq= max(0.,rcp(k,i,j) + rrp(k,i,j))
    !
    ! r_sol= max(0.,rsp(k,i,j)+rpp(k,i,j)+ &
    !       rap(k,i,j)+rgp(k,i,j)+  &
    !       rhp(k,i,j))
    !
    ! tempk = theta(k,i,j)*(exner)/cp ! air temp (Kelvin)
    !
    ! if(tempk.le.253) then
    !   fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.))
    !
    !   dfxcdt = 2.83e6*RQCCUTEN(k,i,J)/(cp*amax1(tempk,253.))
    !
    !  RTHCUTEN (k,i,j) = (1./(1.+fxc))*( RTHCUTEN (k,i,j) - thetail(k,i,j)*dfxcdt )
    !
    !     else
    !
    !       fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.))
    !
    !!orig     dfxcdt = 2.5e6*OUTQC(I,K)*cuten(i)/(cp*amax1(tempk,253.)) - &
    !!orig         fxc/(cp*amax1(tempk,253.)) * cp * OUTT(I,K)
    !  !
    !  dfxcdt = 2.5e6*RQCCUTEN(k,i,J)/(cp*amax1(tempk,253.)) - &
    !             fxc/(cp*amax1(tempk,253.)) * cp * OUTT
    !
    !       RTHCUTEN (k,i,j) = (1./(1.+fxc))*( RTHCUTEN (k,i,j) - thetail(k,i,j)*dfxcdt )
    !
    !     endif
    !
    ! enddo; enddo; enddo
    !endif
    !- tendencies at boundaries
    RTHcuten(1,ia:iz,ja:jz)=RTHcuten(2,ia:iz,ja:jz)
    RQVcuten(1,ia:iz,ja:jz)=RQVcuten(2,ia:iz,ja:jz)
    RQCcuten(1,ia:iz,ja:jz)=RQCcuten(2,ia:iz,ja:jz)

    RTHcuten(m1,ia:iz,ja:jz)=RTHcuten(m1-1,ia:iz,ja:jz)
    RQVcuten(m1,ia:iz,ja:jz)=RQVcuten(m1-1,ia:iz,ja:jz)
    RQCcuten(m1,ia:iz,ja:jz)=RQCcuten(m1-1,ia:iz,ja:jz)

  end subroutine conv_grell_spread3d_brams

  !-------------------------------------------------------------

  subroutine StoreNamelistFileAtCup_grell3(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile

    g3d_spread = oneNamelistFile%g3d_spread
    g3d_smoothh = oneNamelistFile%g3d_smoothh
    g3d_smoothv = oneNamelistFile%g3d_smoothv

  end subroutine StoreNamelistFileAtCup_grell3

  !-----------------------------------------------------------------------------------

  subroutine moveup(m1,A)
    integer, intent(in) :: m1
    real, dimension(m1) :: A,B

    real :: dummy
    integer :: k,kr

    B=A
    do k=1,m1-1
       kr=k+1
       A(kr) = B(k)
    enddo
    A(1) = A(2)
  end subroutine moveup
  !-----------------------------------------------------------------------------------


  subroutine set_index_loops( ims,ime, jms,jme, kms,kme,    &
       its,ite, jts,jte, kts,kte,    &
       mxp,myp,mzp                   )

    integer, intent(IN)         :: mxp,myp,mzp
    integer, intent(INOUT)      :: ims,ime, jms,jme, kms,kme,&
         its,ite, jts,jte, kts,kte


    ims=1
    ime=mxp
    jms=1
    jme=myp
    kms=1
    kme=mzp
    its=1
    ite=mxp
    jts=1
    jte=myp
    kts=1
    kte=mzp

  end subroutine set_index_loops
  !*************************************************************************************
  subroutine check (m1,tht,ath,rtt,artt)
    integer, intent(in) :: m1
    real, dimension(m1), intent(in) :: tht,ath,rtt,artt

    integer k
    do k=1,m1
       print*,"check",k, tht(k),ath(k),rtt(k),artt(k)
    enddo
    print*,"mx1",maxval(tht),maxval(ath),minval(tht),minval(ath)
    print*,"mx2",maxval(rtt),maxval(artt),minval(rtt),minval(artt)

  end subroutine check
  !------------------------------------------------------------
  subroutine cupar2mcphysics(m1,m2,m3,ia,iz,ja,jz,ngrid,dtlt &
       ,clsrc  ,theta,pp,pi0,dn0,oneMicControl,oneMicroFields)
    type(MicControl), pointer, intent(in) :: oneMicControl
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    integer m1,m2,m3,ia,iz,ja,jz,k,i,j,ngrid
    real dtlt
    real, dimension(m1,m2,m3),intent(in) :: theta, pp, pi0,dn0
    real, dimension(m1,m2,m3),intent(in) :: clsrc! liquid/ice tendency from
    ! cumulus parameterization
    if(oneMicControl%level < 2  .and. oneMicControl%mcphys_type < 2 ) return

    if(oneMicControl%level == 2 .and. oneMicControl%mcphys_type < 2) then
       call mcphysics0(m1,m2,m3,ia,iz,ja,jz,dtlt &
            ,clsrc                  &
            ,tend%rct            &
            ,tend%rtt            )

    elseif(oneMicControl%level == 3 .and. (oneMicControl%mcphys_type >= 0))  then

       call mcphysics1(oneMicControl%mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi  &
            ,theta, pp, pi0,dn0 &
            ,clsrc          &! cumulus tendency
            ,tend%rct    &! cloud water mass mix ratio tendency 
            ,tend%rpt    &! pristine mass mix ratio tendency 
            ,tend%rtt    &! total water mass mix ratio tendency
            )

       if(oneMicControl%mcphys_type == 2) &
            call mcphysics2(oneMicControl%mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi  &
            ,theta, pp, pi0,dn0 &
            ,clsrc          &! cumulus tendency
            ,tend%rct    &! cloud water mass mix ratio tendency 
            ,tend%rpt    &! pristine mass mix ratio tendency 
            ,tend%rtt    &! total water mass mix ratio tendency
            ,tend%cpt    &! pristine number conc tendency 
            )

       if(oneMicControl%mcphys_type == 3) &
            call mcphysics3(oneMicControl%mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi  &
            ,theta, pp, pi0,dn0,oneMicroFields%ccp &
            ,clsrc          &! cumulus tendency
            ,tend%rct    &! cloud water mass mix ratio tendency 
            ,tend%rpt    &! pristine mass mix ratio tendency 
            ,tend%rtt    &! total water mass mix ratio tendency
            ,tend%cpt    &! pristine number conc tendency 
            ,tend%cct    &! cloud water  number conc tendency 
            )

    endif
    return
  end  subroutine cupar2mcphysics

  !------------------------------------------------------------------------
  subroutine mcphysics0(m1,m2,m3,ia,iz,ja,jz,dtlt,clsrc,rct,rtt)
    integer m1,m2,m3,ia,iz,ja,jz,k,i,j
    real dtlt
    real, dimension(m1,m2,m3),intent(in   ) :: clsrc
    real, dimension(m1,m2,m3),intent(inout) :: rct,rtt

    do j = ja,jz
       do i = ia,iz
          do k = 1,m1
             rct(k,i,j)=rct(k,i,j)+clsrc(k,i,j)
             rtt(k,i,j)=rtt(k,i,j)+clsrc(k,i,j)
          enddo
       enddo
    enddo
  end subroutine mcphysics0

  !------------------------------------------------------------------------
  subroutine mcphysics1(mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi,theta,pp,pi0,dn0 &
       ,clsrc,rct,rpt,rtt) 

    integer :: m1,m2,m3,ia,iz,ja,jz,k,i,j,mcphys_type
    real dtlt,cpi
    real, dimension(m1,m2,m3),intent(in)    :: clsrc
    real, dimension(m1,m2,m3),intent(in)    :: theta,pp,pi0,dn0
    real, dimension(m1,m2,m3),intent(inout) :: rct,rpt,rtt 
    real ::tempk,tem1,tqice,tqliq,add_ncp,add_npp
    real, parameter :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf)

    do j = ja,jz
       do i = ia,iz
          do k = 1,m1
             tempk = theta(k,i,j)*(pp(k,i,j)+pi0(k,i,j))*cpi ! air temp (Kelvin)

             tem1 = fract_liq_f(tempk)

             !- splitting cumulus tendency into water and ice tendencies
             rct(k,i,j) = rct(k,i,j)+clsrc(k,i,j)* tem1 ! cloud water

             rpt(k,i,j) = rpt(k,i,j)+clsrc(k,i,j)*(1.-tem1) ! pristine ice

             !- it must include also the ice/liq tendencies at rtt for
             !- consistency, since rtt includes ice and liq mixing ratios
             rtt(k,i,j) = rtt(k,i,j)+clsrc(k,i,j)

          enddo
       enddo
    enddo

  end subroutine mcphysics1
  !------------------------------------------------------------------------
  subroutine mcphysics2(mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi,theta,pp,pi0,dn0 &
       ,clsrc,rct,rpt,rtt,cpt)

    integer :: m1,m2,m3,ia,iz,ja,jz,k,i,j,mcphys_type
    real dtlt,cpi
    real, dimension(m1,m2,m3),intent(in)    :: clsrc
    real, dimension(m1,m2,m3),intent(in)    :: theta,pp,pi0,dn0
    real, dimension(m1,m2,m3),intent(inout) :: rct,rpt,rtt,cpt
    real ::tempk,tem1,tqice,tqliq,add_ncp,add_npp
    real, parameter :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf)

    do j = ja,jz
       do i = ia,iz
          do k = 1,m1
             tempk = theta(k,i,j)*(pp(k,i,j)+pi0(k,i,j))*cpi ! air temp (Kelvin)

             tem1 = fract_liq_f(tempk)

             !-- detrained pristine mass mixing ratio
             tqice = (1.-tem1) * clsrc(k,i,j) * dn0(k,i,j)* dtlt

             !-- detrained ICN ice number concenration in the time "dtlt"
             add_npp = max(0.0, make_IceNumber(tqice, tempk)/dn0(k,i,j))

             !- update tendency 
             cpt(k,i,j) = cpt(k,i,j)+ add_npp/dtlt

          enddo
       enddo
    enddo
  end subroutine mcphysics2

  !------------------------------------------------------------------------
  subroutine mcphysics3(mcphys_type,m1,m2,m3,ia,iz,ja,jz,dtlt,cpi,theta,pp,pi0,dn0 &
       ,ccp,clsrc,rct,rpt,rtt,cpt,cct)

    integer :: m1,m2,m3,ia,iz,ja,jz,k,i,j,mcphys_type
    real dtlt,cpi
    real, dimension(m1,m2,m3),intent(in)    :: clsrc
    real, dimension(m1,m2,m3),intent(in)    :: theta,pp,pi0,dn0,ccp
    real, dimension(m1,m2,m3),intent(inout) :: rct,rpt,rtt,cpt,cct
    real ::tempk,tem1,tqice,tqliq,add_ncp,add_npp
    real, parameter :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf)

    do j = ja,jz
       do i = ia,iz
          do k = 1,m1
             tempk = theta(k,i,j)*(pp(k,i,j)+pi0(k,i,j))*cpi ! air temp (Kelvin)

             tem1 = fract_liq_f(tempk)

             !-- detrained pristine mass mixing ratio
             tqice = (1.-tem1) * clsrc(k,i,j) * dn0(k,i,j)* dtlt

             !-- detrained ICN ice number concenration in the time "dtlt"
             add_npp = max(0.0, make_IceNumber(tqice, tempk)/dn0(k,i,j))

             !- update tendency 
             cpt(k,i,j) = cpt(k,i,j)+ add_npp/dtlt

             !--cloud number concentration
             tqliq =  tem1 * clsrc(k,i,j) * dn0(k,i,j) * dtlt
             !--check if

             add_ncp = make_DropletNumber(tqliq, ccp(k,i,j))/dn0(k,i,j)
             !or
             !    add_ncp = make_DropletNumber(tqliq, nwfa(k,i,j))/dn0(k,i,j)

             !- update tendency 
             cct(k,i,j) = cct(k,i,j)+ max(0.0, add_ncp/dtlt)

          enddo
       enddo
    enddo

  end subroutine mcphysics3

  !------------------------------------------------------------------------
  subroutine prepare_lsf(nnqparm,nnshcu,iwork, &
       oneNamelistFile, oneBasicFields, oneRadiateFields)
    character(len=3) :: forcing
    integer,intent(IN) :: nnqparm,nnshcu,iwork
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields

    !- scratchs (local arrays)
    real :: vt3da(mzp,mxp,myp)
    real :: vt3db(mzp,mxp,myp)
    real :: vt3dc(mzp,mxp,myp)
    real :: vt3dh(mzp,mxp,myp)
    real :: vt3dj(mzp,mxp,myp)
    real :: vt3dk(mzp,mxp,myp)
    real :: vt3di(mzp,mxp,myp)
    real :: vt3df(mzp,mxp,myp)
    real :: vt3dg(mzp,mxp,myp)
    real :: vt3de(mzp,mxp,myp)
    real :: vt3dd(mzp,mxp,myp)
    ! real :: vctr1(mzp)
    ! real :: vctr2(mzp)
    real :: scr1(mzp,mxp,myp)
    integer :: i,j,k
    !- parameter to define if include or not diffusion tendencies at forcing for deep convection
    logical,parameter :: forc_deep_pbl = .false.


    if(mod(time,oneNamelistFile%confrq).lt.dtlt .or. &
         time .lt. dtlt+.01) then

       !-
       !  the forcing for shallow is only due to diffusion in PBL only (which is calculated in turb routines)
       !  the forcing for deep is due to radiation + 3dim advection
       if(iwork.eq.1) then

          !----------- include radiation for theta
          if(oneNamelistFile%ilwrtyp + oneNamelistFile%iswrtyp > 0  .and. nnqparm /= 8 ) then
             cuforc_g(ngrid)%lsfth(1:mzp,1:mxp,1:myp)= oneRadiateFields%fthrd(1:mzp,1:mxp,1:myp)
          else
             cuforc_g(ngrid)%lsfth(1:mzp,1:mxp,1:myp)= 0.
          endif

          !-reset lsf for water vapor
          cuforc_g(ngrid)%lsfrt(1:mzp,1:mxp,1:myp)= 0.

          !----------- include advection for theta and rv (or should be rtp?)
          vt3dd=0.0
          vt3de=0.0
          vt3df=0.0
          vt3dg=0.0
          vt3dh=0.0
          vt3di=0.0
          vt3dj=0.0
          vt3dk=0.0
          vctr1=0.0
          vctr2=0.0
          if(dyncore_flag == 0) then
             do j = 1,myp
                do i = 1,mxp
                   do k = 1,mzp
                      vt3da(k,i,j) = (oneBasicFields%up(k,i,j)+ oneBasicFields%uc(k,i,j))*dtlt*0.5
                      vt3db(k,i,j) = (oneBasicFields%vp(k,i,j)+ oneBasicFields%vc(k,i,j))*dtlt*0.5
                      vt3dc(k,i,j) = (oneBasicFields%wp(k,i,j)+ oneBasicFields%wc(k,i,j))*dtlt*0.5
                   end do
                end do
             end do
          else
             do j = 1,myp
                do i = 1,mxp
                   do k = 1,mzp
                      vt3da(k,i,j) = oneBasicFields%uc(k,i,j)*dtlt
                      vt3db(k,i,j) = oneBasicFields%vc(k,i,j)*dtlt
                      vt3dc(k,i,j) = oneBasicFields%wc(k,i,j)*dtlt
                   end do
                end do
             end do
          endif
          call fa_preptc(mzp,mxp,myp            &
               ,vt3da     ,vt3db       &
               ,vt3dc     ,vt3dd       &
               ,vt3de     ,vt3df       &
               ,vt3dh     ,vt3di       &
               ,vt3dj     ,vt3dk       &
               ,mynum, &
               oneBasicFields)

          if(dyncore_flag == 0) then
             !---- thp
             scr1(1:mzp,1:mxp,1:myp) = oneBasicFields%thp(1:mzp,1:mxp,1:myp)

             ! output: scr1,vt3dg
             call fa_xc(mzp,mxp,myp,ia,iz,1,myp,oneBasicFields%thp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)

             ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
             ! output: scr1,vt3dg
             if (jdim .eq. 1)  &
                  call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,oneBasicFields%thp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)

             ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
             ! output: scr1,vt3dg
             call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,oneBasicFields%thp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)

             ! input:  thetap , lsfth,scr1, dtlt
             ! output: lsfth
             call advtndc(mzp,mxp,myp,ia,iz,ja,jz,oneBasicFields%thp,scr1,cuforc_g(ngrid)%lsfth,dtlt,mynum)
             !
          else
             !---- thc
             scr1(1:mzp,1:mxp,1:myp) = oneBasicFields%thc(1:mzp,1:mxp,1:myp)

             ! output: scr1,vt3dg
             call fa_xc(mzp,mxp,myp,ia,iz,1,myp,oneBasicFields%thc,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)

             ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
             ! output: scr1,vt3dg
             if (jdim .eq. 1)  &
                  call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,oneBasicFields%thc,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)

             ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
             ! output: scr1,vt3dg
             call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,oneBasicFields%thc,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)

             ! input:  thetac , lsfth,scr1, dtlt
             ! output: lsfth
             call advtndc(mzp,mxp,myp,ia,iz,ja,jz,oneBasicFields%thc,scr1,cuforc_g(ngrid)%lsfth,dtlt,mynum)
          endif

          !---- water vapor
          scr1(1:mzp,1:mxp,1:myp) = oneBasicFields%rv(1:mzp,1:mxp,1:myp)

          ! output: scr1,vt3dg
          call fa_xc(mzp,mxp,myp,ia,iz,1,myp,oneBasicFields%rv,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)

          ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
          ! output: scr1,vt3dg
          if (jdim .eq. 1)  &
               call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,oneBasicFields%rv,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)

          ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
          ! output: scr1,vt3dg
          call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,oneBasicFields%rv,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)

          ! input: basic(ngrid)%rv, scalart,scr1, dtlt
          ! output: lsfrt = rad + adv
          call advtndc(mzp,mxp,myp,ia,iz,ja,jz,oneBasicFields%rv,scr1,cuforc_g(ngrid)%lsfrt,dtlt,mynum)

          !- here the forcings contain rad+adv for temp and adv for water vapor
          !-end of inclusion of the advection forcings

          !- flag to include diffusion in pbl (only vertical) in deep convection forcing.
          !- the pbl forcing (cuforc_sh_g(ngrid)%lsfth and %lsfrt) has been calculated
          !- in the turbulence routines.
          if(forc_deep_pbl) then
             ! calcula o forcing para conveccao profunda = rad + pbl turb + adv
             ! for deep convection    LSF =  radiation + pbl_turb + advection
             cuforc_g(ngrid)%lsfth(:,:,:) = cuforc_g(ngrid)%lsfth(:,:,:)+cuforc_sh_g(ngrid)%lsfth(:,:,:)
             cuforc_g(ngrid)%lsfrt(:,:,:) = cuforc_g(ngrid)%lsfrt(:,:,:)+cuforc_sh_g(ngrid)%lsfrt(:,:,:)
          endif
          !if(mynum==1) print*,"2:lsfth",maxval(cuforc_g(ngrid)%lsfth),minval(cuforc_g(ngrid)%lsfth)
          !if(mynum==1) print*,"2:lsfrt",maxval(cuforc_g(ngrid)%lsfrt),minval(cuforc_g(ngrid)%lsfrt)

       endif
       if(iwork.eq.2) then

          call accum(int(mxp*myp*mzp,i8), cuforc_g(ngrid)%lsfth, cuparm_g_sh(ngrid)%thsrc)
          call accum(int(mxp*myp*mzp,i8), cuforc_g(ngrid)%lsfrt, cuparm_g_sh(ngrid)%rtsrc)

       endif

    endif

  end subroutine prepare_lsf



  !-------------------------------------------------------------------
  subroutine get_zi_gf2018(m1,tkmin,tkeg,z,rtgt,ztop,kzi)

    integer,intent(in):: m1
    integer :: kzimax,ktke_max,i,k
    real tkmin,tke_tmp
    real,intent(in), dimension(m1) :: tkeg,z
    real,intent(in)    :: ztop,rtgt
    integer,intent(out) :: kzi
    real, parameter :: rcpmin=1.e-5 , pblhmax=3000.

    kzi      = 1
    ktke_max = 1
    kzimax   = m1-1
    !---  max level for kzi
    do k=1,m1
       if(z(k).ge. pblhmax+ztop) then
          kzimax = min(k,m1-1)
          !if(j==8 .and. i==10) print*,"1",z(i,k), pblhmax,ztop(i),kzimax
          exit
       endif
    enddo
    !---
    do k=ktke_max,kzimax
       if(tkeg(k) .gt. 1.1*tkmin )  then
          kzi = k
          cycle
       else
          kzi = max(1,k-1)
          exit
       endif
    enddo
    kzi = max(1  ,kzi)
    kzi = min(kzimax,kzi)
    !print*,"2",kzi(i),i;call flush(6)
    !pbl(i) = max( z(i,kzi(i))-ztop(i), z(i,1)-ztop(i) )
  end subroutine get_zi_gf2018
  !-------------------------------------------------------------------
end module ModCuParGrell3
