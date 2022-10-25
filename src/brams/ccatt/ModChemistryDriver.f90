module ModChemistryDriver

  use rconstants , only: &
       cp,               & ! (IN)
       cpor,             & ! (IN)
       pi180,            & ! (IN)
       p00

  use mem_grid, only: &
       ibnd, &
       jdim, &
       jbnd, &
       nz, &
       time,          & ! (IN)
       iyear1,        & ! (IN)
       imonth1,       & ! (IN)
       idate1,        & ! (IN)
       itime1,        & ! (IN)
       dtlt,          & ! (IN)
       ngrid,         & ! (IN)
       ngrids,        & ! (IN)
       dtlongn,       & ! (IN)
       dzt,           & ! (IN)
       zm,            & ! (IN)
       initial,       & ! (IN)
       zt,            & ! (IN)
       centlat,       & ! (IN)
       centlon,       & ! (IN)
       grid_g,        & ! %glat(IN), %glon(IN), %rtgt(IN), %lpw(IN), %dxt(IN)
       platn,         & ! (IN)
       plonn,         & ! (IN)
       deltaxn,       & ! (IN)
       deltayn,       & ! (IN)
       xt,            & ! (IN)
       yt               ! (IN)

  use chem1_list, only :          &
       nx=>nspecies, &
       init_ajust, &
       CO2, &
       chem_nspecies =>nspecies,  & ! (IN)
       spc_chem_alloc=>spc_alloc, & ! (IN)
       spc_chem_name =>spc_name,  & ! (IN)
       src,                       & ! (IN)
       on,                        & ! (IN)
       off,                       & ! (IN)
       transport,                 & ! (IN)
       nspecies,                  & ! (IN)
       no,                        & ! (IN)
       weight,                    & ! (IN)
       co,                        & ! (IN)
       photojmethod,              & ! (IN)
       nr,                        & ! (IN)
       nr_photo                     ! (IN)

  use mem_chem1, only:               &
       CHEM_ASSIM, &
       CHEMISTRY,                    & ! (IN)
       maxnspecies,                  & ! (IN)
       nspecies_chem_transported,    & ! (IN)
       nspecies_chem_no_transported, & ! (IN)
       transp_chem_index,            & ! (IN)
       no_transp_chem_index,         & ! (IN)
       nsrc,                         & ! (IN)
       src_name,                     & ! (IN)
       antro,                        & ! (IN)
       bburn,                        & ! (IN)
       bioge,                        & ! (IN)
       geoge,                        & ! (IN)
       chem1_src_z_dim_g,            & ! (IN)
       chem1_g,                      & ! %sc_t(INOUT), %sc_p(INOUT)
       chem1_src_g,                  & ! %sc_src(INOUT)
       chem1_src_vars,               & ! Type
       chem1_vars,                   & ! Type
       N_DYN_CHEM_N,                 & ! (IN)
       N_DYN_CHEM,                   & ! (INOUT)
       split_method,                 & ! (IN)
       isplit                          ! (IN)


  use mem_chem1aq, only: &
       chem1aq_g, &
       CHEMISTRY_AQ        ! (IN)

  use ModBasicFields, only: &
       BasicFields

  use ModRadiateFields, only: &
       RadiateFields
  
  use FastJX, only :  &
       FastJX_driver, & ! Subroutine
       Initialize_Fast_JX, &
       fast_JX_initialized, &
       fast_JX_g        ! %jphoto(IN)

  use node_mod, only : &
       mynum,          & ! (IN)
       nodemxp,        &
       nodemyp,        &
       nodemzp,        &
       nodeia,         & ! (IN)
       nodeiz,         & ! (IN)
       nodeja,         & ! (IN)
       nodejz,         & ! (IN)
       i0,             & ! (IN)
       j0,             & ! (IN)
       ibcon             ! (IN)

  use ModNamelistFile, only: &
       NamelistFile
  
  use ModMicroFields, only: &
       MicroFields

  use grid_dims, only: &
       nzpmax,         &  ! (IN)
       nxpmax,         &  ! (IN)
       nypmax             ! (IN)

  use mem_scratch1_grell, only: &
       ierr4d                     ! (IN)

  use mem_grell_param, only: &
       mgmxp,                & ! (IN)
       mgmyp,                & ! (IN)
       maxiens,              & ! (IN)
       ngrids_cp               ! (IN)

  use mem_stilt, only: &
       iexev,          & ! (IN)
       stilt_g           ! %dnp(IN)

  use mem_globrad, only: &
       raddatfn            ! (IN)

  use carma_fastjx, only: &
       do3,               & ! (IN)
       daer,              & ! (IN)
       na                   ! (IN)

  use mem_aer1, only:                   &
       aer1_g, &
       AEROSOL, &
       aerosol, &
       aer_nvert_src=>aer1_src_z_dim_g, & ! (IN)
       aer1_vars                          ! Type

  use spack_utils, only: &
       allocindex,       & ! Subroutine
       index_g,          & ! IN
       nob,              & ! IN
       maxblock_size       ! IN

  use mem_spack, only: &
       alloc_spack,    & ! Subroutine
       spack_alloc,    & ! IN
       spack,          & ! OUT
       spack_2d,       & ! INOUT
       atol,           & ! IN
       rtol              ! IN

  use solve_sparse, only: &
       get_non_zeros

  use mod_chem_spack_ros_dyndt, only: &
       chem_ros_dyndt    ! Subroutine

  use mod_chem_spack_ros, only: &
       chem_ros          ! Subroutine

  use mod_chem_spack_qssa, only: &
       chem_qssa         ! Subroutine

  use mod_chem_trans_gasaq, only: &
       trans_gaq         ! Subroutine

  use mod_chem_spack_rodas3_dyndt, only: &
       chem_rodas3_dyndt ! Subroutine

  use uv_atten, only: &
       initialize_uv_atten,  & ! Subroutine
       uv_attenuation,       & ! Subroutine
       uv_atten_g,           & ! (IN)
       uv_atten_initialized    ! (INOUT)

  use mod_chem_trans_liq, only: &
       trans_liq         ! Subroutine

  use chem1aq_list, only: &
       nspeciesaq,        & ! (IN)
       ind_gas              ! (IN)

  use mem_chemic, only: &
       chemic_g            ! (IN)

  use mem_aerad, only: &
       nwave              ! (IN)

  use mem_carma, only: &
       carma              ! (IN)

  use Extras, only: &
       na_extra3d,  &     ! (IN)
       extra3d            ! (INOUT)

  use FastJX, only:  &
       fast_JX_g         ! %jphoto(IN)

  use ModtuvDriver, only: &
       tuvDriver! (IN)

  use mem_rrtm, only: &
       aot_rrtm_lw

  use parrrtm, only : &
       nbndlw

  use aer1_list, only: &
       coarse, &
       urban, &
       nspecies_aer   =>nspecies,  &
       spc_alloc_aer  =>spc_alloc

  implicit none

  private

  public :: chemistry_driver ! Subroutine
  public :: aer_background
  public :: initial_condition

contains

  !------------------------------------------------------------------
  subroutine chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,task,nt, &
       oneNamelistFile, oneBasicFields, oneMicroFields, oneRadiateFields)

    integer , intent(IN) :: mzp
    integer , intent(IN) :: mxp
    integer , intent(IN) :: myp
    integer , intent(IN) :: ia
    integer , intent(IN) :: iz
    integer , intent(IN) :: ja
    integer , intent(IN) :: jz
    integer , intent(IN) :: task
    integer , intent(IN) :: nt
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields
    
    integer, dimension(ngrids) :: itchim
    integer, dimension(ngrids) :: ncycle

    integer,parameter :: block_size_qssa=45
    integer,parameter :: block_size_ros=64 ! use 64 for pgi / 48 for intel compilers
    character(len=*), parameter :: h="**(chemistry_driver)**"

    !     INTEGER,PARAMETER :: blockSize=1 !tmp


    if( task == 3 .and. CHEMISTRY >= 0) then

       if (CHEMISTRY == 0) then

          call include_tracer_lifetime(mzp,mxp,myp,spc_chem_alloc,nspecies,transport,on,co, &
               chem1_g(CO,ngrid)%sc_p(:,:,:), chem1_g(CO,ngrid)%sc_t(:))

       else

          !-  call chemistry production/loss

          !- Determining number of blocks and masks, allocating spack type
          if(.not. spack_alloc) then

             !--(DMK-BRAMS-5.0-INI)----------------------------------------------------------------------
             call allocIndex(block_size_ros,nodemxp(mynum,:),nodemyp(mynum,:),nodemzp(mynum,:),&
                  nodeia(mynum,:),nodeiz(mynum,:),nodeja(mynum,:),nodejz(mynum,:),&
                  !--(DMK-BRAMS-5.0-ORI)----------------------------------------------------------------------
                  !             CALL allocIndex(block_size_ros,mmxp,mmyp,mmzp,mia,miz,mja,mjz,mibcon,&
                  !--(DMK-BRAMS-5.0-FIM)----------------------------------------------------------------------
                  ngrids,dtlongn,N_DYN_CHEM)
             call alloc_spack(CHEMISTRY)
          end if

          if(uv_atten_initialized == 0) then

             !--(DMK-BRAMS-5.0-INI)----------------------------------------------------------------------
             call initialize_uv_atten(ngrids,nodemxp(mynum,:),nodemyp(mynum,:))
             !--(DMK-BRAMS-5.0-OLD)----------------------------------------------------------------------
             !              CALL initialize_uv_atten(ngrids,mmxp,mmyp)
             !--(DMK-BRAMS-5.0-FIM)----------------------------------------------------------------------

             uv_atten_initialized=1
          end if

          if(Fast_JX_initialized == 0) then
             call initialize_fast_JX(nodemxp(mynum,:),nodemyp(mynum,:),nodemzp(mynum,:),ngrids)
             Fast_JX_initialized = 1
          end if

          !- photolysis or attenuation calculations
          if (mod(time + .001,oneNamelistFile%radfrq) .lt. dtlt .or. time .lt. 0.001) then
             if (trim(PhotojMethod) == 'FAST-JX') then

                call FastJX_driver(mzp,mxp,myp,ia,iz,ja,jz,nzpmax,ngrid,ngrids,grid_g(ngrid)%lpw, &
                     !--(DMK-BRAMS-5.0-INI)----------------------------------------------------------------------
                     grid_g(ngrid)%rtgt,zm,imonth1,idate1,iyear1,time,dzt, &
                     nodemzp(mynum,:),nodemxp(mynum,:),nodemyp(mynum,:),i0,j0,ngrids,&
                     oneBasicFields%pp,oneBasicFields%pi0, &
                     oneBasicFields%theta,oneBasicFields%dn0,oneRadiateFields%rlongup, &
                     oneRadiateFields%cosz,oneRadiateFields%albedt,raddatfn,do3,daer,na)

             elseif(trim(PhotojMethod) == 'LUT') then

                if(oneNamelistFile%ilwrtyp==4 .or. oneNamelistFile%iswrtyp==4) then
                   call uv_attenuation(mxp,myp,ia,iz,ja,jz,i0,j0, &
                        platn(ngrid),plonn(ngrid),deltaxn(ngrid),deltayn(ngrid), &
                        nxpmax,nypmax,xt,yt,nwave,carma(ngrid)%aot, &
                        oneNamelistFile%ilwrtyp,oneNamelistFile%iswrtyp,uv_atten_g(ngrid)%att,11)
                elseif(oneNamelistFile%ilwrtyp==6 .or. oneNamelistFile%iswrtyp==6) then
                   call uv_attenuation(mxp,myp,ia,iz,ja,jz,i0,j0, &
                        platn(ngrid),plonn(ngrid),deltaxn(ngrid),deltayn(ngrid), &
                        nxpmax,nypmax,xt,yt,nbndlw,real(aot_rrtm_lw), &
                        oneNamelistFile%ilwrtyp,oneNamelistFile%iswrtyp,uv_atten_g(ngrid)%att,2)
                end if

             elseif(trim(PhotojMethod) == 'FAST-TUV') then
                call tuvDriver(mzp,mxp,myp,ia,iz,ja,jz, &
                     oneBasicFields, oneNamelistFile, oneRadiateFields)

             else
                stop ' unknonw photolysis calculation method '

             endif
             if(trim(PhotojMethod) == 'LUT') then
                call tuvDriver(mzp,mxp,myp,ia,iz,ja,jz, &
                     oneBasicFields, oneNamelistFile, oneRadiateFields)
             endif
          endif

          !-srf: including options for operator splitting
          !- current dynamic splitting
          N_DYN_CHEM=N_DYN_CHEM_N(ngrid)

          if(split_method == 'PARALLEL' .and. N_DYN_CHEM > 1) &
               call chem_accum(mzp,mxp,myp,ia,iz,ja,jz,1,nspecies_chem_transported,chem1_g(:,ngrid), &
               transp_chem_index,N_DYN_CHEM)

          if (mod(time-(N_DYN_CHEM/ISPLIT-1)*dtlt , N_DYN_CHEM*dtlt) == 0. ) then
             !if(mynum==1)print*,'call chem: time=',time,ngrid; call flush(6)

             if(CHEMISTRY == 1) then
                itchim(ngrid)= 75  ! teste 60
                !STOP 'for now use only CHEMISTRY = 4'
                !Determining number of blocks and masks
                !CALL AllocIndex(block_size_qssa,mmxp,mmyp,mmzp,mia,miz,mja,mjz,mibcon, &
                !     ngrids,dtlongn,N_DYN_CHEM)

                call chem_qssa(mzp,mxp,myp,itchim(ngrid),nr,mynum,dtlt,oneBasicFields%pp, &
                     oneBasicFields%pi0, oneBasicFields%theta,oneBasicFields%rv, &
                     oneRadiateFields%cosz,fast_JX_g(ngrid)%jphoto, &
                     nspecies,nr_photo,weight,PhotoJMethod,chem1_g(:,ngrid), &
                     nob(ngrid),maxblock_size, &
                     index_g(ngrid)%block_end,index_g(ngrid)%indexk,index_g(ngrid)%indexi, &
                     index_g(ngrid)%indexj,index_g(ngrid)%kij_index, &
                     cp,cpor,p00, &
                     uv_atten_g(ngrid)%att)

                !Deallocating blocks data
                !CALL AllocIndex(block_size_qssa,mmxp,mmyp,mmzp,mia,miz,mja,mjz,mibcon, &
                !     ngrids,dtlongn,N_DYN_CHEM)


             elseif(CHEMISTRY == 2) then

                call chem_ros(mzp,mxp,myp,dtlt,oneBasicFields%pp,oneBasicFields%pi0, &
                     oneBasicFields%theta,oneBasicFields%rv,oneBasicFields%dn0, &
                     oneRadiateFields%cosz,oneMicroFields%rcp,nspecies,nr,nr_photo,weight, &
                     PhotojMethod,fast_JX_g(ngrid)%jphoto,maxnspecies,nspecies_chem_transported, &
                     nspecies_chem_no_transported,transp_chem_index,no_transp_chem_index, &
                     chem1_g(:,ngrid),nob(ngrid),maxblock_size, &
                     index_g(ngrid)%block_end,index_g(ngrid)%indexk,index_g(ngrid)%indexi, &
                     index_g(ngrid)%indexj,index_g(ngrid)%kij_index, &
                     spack,spack_2d,get_non_zeros,&
                     N_DYN_CHEM,split_method,cp,cpor,p00, &
                     uv_atten_g(ngrid)%att)

             elseif(CHEMISTRY == 3) then

                call chem_ros_dyndt(mzp,mxp,myp,dtlt,oneBasicFields%pp,oneBasicFields%pi0, &
                     oneBasicFields%theta,oneBasicFields%rv,oneBasicFields%dn0, &
                     oneRadiateFields%cosz,oneMicroFields%rcp,nspecies,nr,nr_photo,weight, &
                     PhotojMethod,maxnspecies,nspecies_chem_transported, &
                     nspecies_chem_no_transported,transp_chem_index,no_transp_chem_index, &
                     chem1_g(:,ngrid),nob(ngrid),maxblock_size, &
                     index_g(ngrid)%block_end,index_g(ngrid)%indexk,index_g(ngrid)%indexi, &
                     index_g(ngrid)%indexj,index_g(ngrid)%kij_index,index_g(ngrid)%last_accepted_dt, &
                     spack,spack_2d,atol,rtol,get_non_zeros,&
                     N_DYN_CHEM,split_method,cp,cpor,p00,fast_JX_g(ngrid)%jphoto, &
                     uv_atten_g(ngrid)%att)
             elseif(CHEMISTRY == 4) then

                call chem_rodas3_dyndt(mzp,mxp,myp,dtlt,oneBasicFields%pp,oneBasicFields%pi0, &
                     oneBasicFields%theta,oneBasicFields%rv,oneBasicFields%dn0, &
                     oneRadiateFields%cosz,oneMicroFields%rcp,nspecies,nr,nr_photo,weight, &
                     PhotojMethod,maxnspecies,nspecies_chem_transported, &
                     nspecies_chem_no_transported,transp_chem_index,no_transp_chem_index, &
                     chem1_g(:,ngrid),nob(ngrid),maxblock_size, &
                     index_g(ngrid)%block_end,index_g(ngrid)%indexk,index_g(ngrid)%indexi, &
                     index_g(ngrid)%indexj,index_g(ngrid)%kij_index,index_g(ngrid)%last_accepted_dt, &
                     spack,spack_2d,atol,rtol,get_non_zeros,&
                     N_DYN_CHEM,split_method,cp,cpor,p00,na_extra3d,extra3d(:,ngrid), &
                     fast_JX_g(ngrid)%jphoto,uv_atten_g(ngrid)%att)

             endif

             !- set dyn tend to zero, only for parallel splitting method
             if(split_method == 'PARALLEL' .and. N_DYN_CHEM > 1) &
                  call chem_accum(mzp,mxp,myp,ia,iz,ja,jz,2,nspecies_chem_transported,chem1_g(:,ngrid), &
                  transp_chem_index,N_DYN_CHEM)

          endif

       endif
       ! change MP 11/02/08
    elseif(task == 4  .and. CHEMISTRY  >= 1 .and. CHEMISTRY_AQ >= 1 ) then
       call trans_gaq(mzp,mxp,myp,ia,iz,ja,jz)

    elseif(task == 5  .and. CHEMISTRY  >= 1 .and. CHEMISTRY_AQ >= 1 ) then
       call trans_liq(mzp,mxp,myp,ia,iz,ja,jz,chem1_g(:,ngrid),chem1aq_g(:,ngrid), &
            nspeciesaq,ind_gas,chemic_g(ngrid)%coll,chemic_g(ngrid)%sedimr, &
            oneMicroFields%rcp,oneMicroFields%rrp)
       ! change MP 11/02/08 - end

    endif

  end subroutine chemistry_driver

  !------------------------------------------------------------------------------------------

  subroutine include_tracer_lifetime(m1,m2,m3,spc_alloc_chem,nspecies,&
       transport,on,co,sc_p,sc_t)


    integer , intent(IN)    :: m1
    integer , intent(IN)    :: m2
    integer , intent(IN)    :: m3
    integer , intent(IN)    :: spc_alloc_chem(6,nspecies)
    integer , intent(IN)    :: nspecies
    integer , intent(IN)    :: transport
    integer , intent(IN)    :: on
    integer , intent(IN)    :: co
    real    , intent(IN)    :: sc_p(m1,m2,m3)
    real    , intent(INOUT) :: sc_t(m1*m2*m3)

    ! local declarations
    real,parameter :: vm_CO_i   = 1./(30.*24.*3600.) ! 1/30 days

    ! - only for CO
    if(spc_alloc_chem(transport,CO) == ON) then

       call accum_lifetime(m1*m2*m3     &
            ,sc_t     & ! tend
            ,sc_p & ! mix ratio
            ,vm_CO_i   )   ! inv of lifetime (1/seconds)
    endif

  end subroutine include_tracer_lifetime

  !------------------------------------------------------------------------------------------

  subroutine accum_lifetime(m1,sc_t,sc_p,lf_inv)

    integer , intent(IN)    :: m1
    real    , intent(INOUT) :: sc_t(m1)
    real    , intent(IN)    :: sc_p(m1)
    real    , intent(IN)    :: lf_inv

    integer i

    do i=1,m1
       sc_t(i) = sc_t(i) - lf_inv*sc_p(i)
    enddo
  end subroutine accum_lifetime

  !------------------------------------------------------------------------------------------

  subroutine chem_accum(m1,m2,m3,ia,iz,ja,jz,ntask,nspecies_chem_transported,chem1_g, &
       transp_chem_index,N_DYN_CHEM)

    integer         , intent(IN) :: ntask
    integer         , intent(IN) :: m1
    integer         , intent(IN) :: m2
    integer         , intent(IN) :: m3
    integer         , intent(IN) :: ia
    integer         , intent(IN) :: iz
    integer         , intent(IN) :: ja
    integer         , intent(IN) :: jz

    ! mem_chem1
    integer          , intent(IN)    :: nspecies_chem_transported
    type (chem1_vars), intent(INOUT) :: chem1_g(nspecies_chem_transported)
    integer          , intent(INOUT) :: transp_chem_index(nspecies_chem_transported)
    integer          , intent(IN)    :: N_DYN_CHEM


    !- local var
    integer ispc,n,ixyz,ntps,i,j,k
    real N_DYN_CHEM_i

    ntps = m1 * m2 * m3

    if(ntask == 1) then

       N_DYN_CHEM_i= real(1./N_DYN_CHEM)

       do ispc=1,nspecies_chem_transported

          !- map the species to transported ones
          n=transp_chem_index(ispc)

          !- calculate the mean dynamic tendency for the entire chemistry timestep
          do ixyz= 1,ntps
             chem1_g(n)%sc_t_dyn(ixyz) = chem1_g(n)%sc_t_dyn(ixyz) + &
                  N_DYN_CHEM_i * chem1_g(n)%sc_t(ixyz)
             !call accum(ntps, chem1_g(n)%sc_t_dyn(1), chem1_g(n)%sc_t(1) )
          enddo

       enddo

    else
       do ispc=1,nspecies_chem_transported

          n=transp_chem_index(ispc) !- map the species to transported ones

          !- set to zero arrays for  accumulation over the next chem timestep
          chem1_g(n)%sc_t_dyn(1:ntps) = 0.

       enddo
       !do j=ja,jz ; do i=ia,iz ;do k=2,m1
       !extra3d(9 ,ngrid)%d3(k,i,j)=0.
       !extra3d(10,ngrid)%d3(k,i,j)=0.
       !extra3d(11,ngrid)%d3(k,i,j)=0.
       !enddo;enddo;enddo

    endif
  end  subroutine chem_accum
  !------------------------------------------------------------------


  subroutine latset_tracer(m1,m2,m3,ia,iz,ja,jz,ibcon,i0,j0,mynum,ap,uc,vc,dxu, &
       dxm,dyv,dym)
    integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k,lbw,lbe,lbs,lbn,i0,j0,mynum
    real :: thresh,dtlx,c1,dxr,dyr
    real, dimension(m1,m2,m3) :: ap,uc,vc
    real, dimension(m2,m3) :: dxu,dxm,dyv,dym
    !character(len=*) :: vnam
    real :: vctr17(m1)
    real :: vctr18(m1)

    if (iand(ibcon,1) .gt. 0) lbw = ia - 1
    if (iand(ibcon,2) .gt. 0) lbe = iz + 1
    if (iand(ibcon,4) .gt. 0) lbs = ja - 1
    if (iand(ibcon,8) .gt. 0) lbn = jz + 1

    thresh = 0.
    dtlx = dtlt


    if (ibnd .ne. 4) then

       ! Western boundary for lsflg = 2  == constant inflow, radiative b.c. outflow
       ! Veja que o campo ap(k,i,j) na borda so e' modificado se houver outflow, isto e'
       ! em regioes de inflow o valor permanece constante e igual ao valor inicial (t=0).
       ! Para introduzir uma especie de "nudging" na condicao de contorno, modifique o
       ! vetor ap(k,lbw,j) para a 'western boundary', por exemplo. O mesmo para as outras
       ! bordas.

       if (iand(ibcon,1) .gt. 0) then
          do j = 1,m3
             dxr = dxu(ia,j) / dxu(lbw,j)
             c1 = dtlx * dxu(lbw,j)
             do k = 1,m1
                vctr17(k) = -c1 * uc(k,lbw,j)
                vctr18(k) = ap(k,ia,j) + dxr * (ap(k,ia,j) - ap(k,ia+1,j))
             enddo
             do k = 1,m1
                if (vctr17(k) .ge. thresh) ap(k,lbw,j) = vctr18(k)
                !      !---srf-tmp
                !      if(i+i0==12 .and. j+j0==1 .and. k==2) then
                !       print*,'1 ap=',ap(k,lbw,j),lbw,vctr17(k)
                !       call flush(6)
                !      endif
                !      !---srf-tmp
             enddo
          enddo
       endif

       !     Eastern Boundary for LSFLG =  2

       if (iand(ibcon,2) .gt. 0) then
          do j = 1,m3
             dxr = dxu(iz-1,j) / dxu(iz,j)
             c1 = dtlx * dxu(iz,j)
             do k = 1,m1
                vctr17(k) = c1 * uc(k,iz,j)
                vctr18(k) = ap(k,iz,j) + dxr * (ap(k,iz,j) - ap(k,iz-1,j))
             enddo
             do k = 1,m1
                if (vctr17(k) .ge. thresh)  ap(k,lbe,j) = vctr18(k)
                !      !---srf-tmp
                !      if(i+i0==12 .and. j+j0==1 .and. k==2) then
                !       print*,'2 ap=',ap(k,lbe,j),lbe,vctr17(k)
                !       call flush(6)
                !      endif
                !      !---srf-tmp
             enddo
          enddo
       endif
    endif

    if(jdim.eq.1.and.jbnd.ne.4)then

       !     Southern boundary for LSFLG  2

       if (iand(ibcon,4) .gt. 0) then
          do i = 1,m2
             dyr = dyv(i,ja) / dyv(i,lbs)
             c1 = dtlx * dyv(i,lbs)
             !srf - fix from rams 60
             !             do k = 1,nz
             do k = 1,m1
                vctr17(k) = -c1 * vc(k,i,lbs)
                vctr18(k) = ap(k,i,ja) + dyr * (ap(k,i,ja) - ap(k,i,ja+1))
             enddo
             do k = 1,m1
                if (vctr17(k) .ge. thresh) ap(k,i,lbs) = vctr18(k)
                !      !---srf-tmp
                !      if(i+i0==12 .and. k==2) then
                !       print*,'4 ap=',ap(k,i,lbs),lbs,vctr17(k),mynum; call flush(6)
                !       write(mynum*10,100) lbs,i,i0,ja,mynum,ap(k,i,lbs),vctr17(k)*100000,ap(k,i,ja), ap(k,i,ja+1)
                !       100 format(1x,'4 ap=',5i4,5F10.2)
                !       call flush(mynum*10)
                !      endif
                !      !---srf-tmp
             enddo
          enddo
       endif

       !     Northern Boundary for LSFLG =  2

       if (iand(ibcon,8) .gt. 0) then
          do i = 1,m2
             dyr = dyv(i,jz-1) / dyv(i,jz)
             c1 = dtlx * dyv(i,jz)
             do k = 1,m1
                vctr17(k) = c1 * vc(k,i,jz)
                vctr18(k) = ap(k,i,jz) + dyr * (ap(k,i,jz) - ap(k,i,jz-1))
             enddo
             do k = 1,m1
                if (vctr17(k) .ge. thresh) ap(k,i,lbn) = vctr18(k)
                !---srf-tmp
                !if(i+i0==12  .and. k==2) then
                ! print*,'8 ap=',ap(k,i,lbn),lbn,vctr17(k)
                ! call flush(6)
                !endif
                !---srf-tmp
             enddo
          enddo
       endif
    endif

  end subroutine latset_tracer
  !-----------------------------------------------------------------------


  !------------------------------------------------------------------
  subroutine initial_condition(ng,m1,m2,m3)
    integer :: ng,m1,m2,m3,i,j,k,ispc
    real,parameter :: fcu =1.e+9 !=> ppbm
    real dummy,dum1,dum2,b
    integer ncount, icx1,jcx1,icx2,jcx2

    return
    chem1_g(CO2,ng)%sc_p(:,:,:)=390.* 44.0 / 28.97

    !--(DMK-CCATT-INI)-----------------------------------------------------
    return    ! Corrige distorcao no CO e NO quando CHEM_ASSIM==0
    !--(DMK-CCATT-OLD)-----------------------------------------------------
    !    IF(CHEM_ASSIM > 0 ) RETURN
    !--(DMK-CCATT-FIM)-----------------------------------------------------


    !-attention: values must be in ppbv (check if this correct, should
    !                                    not be in ppbm? SRF)
    !PRINT*,'============================= homogeneous initial. ====='
    !PRINT*,'CO= ',chem1_g(CO ,ng)%sc_p(1:m1,72,68)*fcu
    !PRINT*,'========================================================'

    !   aer1_g(accum,bburn,ng)%sc_p(10:15,1:m2,1:m3) = 300.
    !do k=11,m1
    !  chem1_g(CO,ng)%sc_p(k,:,:)=0.
    !enddo
    !chem1_g(O3,ng)%sc_p=chem1_g(CO,ng)%sc_p
    !----------------------------------  1 square   ----------------------------------
    !chem1_g(CO,ng)%sc_p= 0.
    !chem1_g(CO,ng)%sc_p(5:20,25:50,10:30) = 100.*28./28.96
    !chem1_g(NO,ng)%sc_p= 20.*28./28.96
    !chem1_g(NO,ng)%sc_p(5:20,25:50,10:30) = 100.*28./28.96

    !----------------------------------  1 square fim---------------------------
    !icx1=50
    !jcx1=50
    !icx2=50
    !jcx2=50
    !b=1./90.
    ! b=1./150.
    !b=1./450.
    !b=1./280.
    !do j=1,m3; do i=1,m2
    !  !gaussian
    !  !dum1 = (max( 0., exp(-b*( real(i-icx2)**2+ real(j-jcx2)**2 ) ) )*100)*28./28.96
    !  dum2 = (max( 0., exp(-b*( real(i-icx1)**2+ real(j-jcx1)**2 ) ) )*100)*28./28.96
    !  !
    !  !if(dum1>0.01) chem1_g(CO,ng)%sc_p(:,i,j)=dum1
    !  !if(dum2>0.01)
    !  chem1_g(CO,ng)%sc_p(:,i,j)=dum2
    !----------------------------
    !enddo;enddo
    !do j=1,m3; do i=1,m2;do k=1,m1
    !
    !  if(chem1_g(CO,ng)%sc_p(k,i,j)<20.*28./28.96)  chem1_g(CO,ng)%sc_p(k,i,j)=20.*28./28.96
    !enddo;enddo;enddo
    !return
    !----------------------------------  squares   ----------------------------------
    ! chem1_g(CO,ng)%sc_p= 20.*28./28.96
    !- dois quadrados
    ! chem1_g(CO,ng)%sc_p(:,40:60,75:95) = 100.*28./28.96
    ! chem1_g(CO,ng)%sc_p(:,40:60,05:25) = 100.*28./28.96
    !- for linear correlations
    ! chem1_g(Ch4,ng)%sc_p(:,:,:)= 12.4 + 1.53*chem1_g(CO,ng)%sc_p(:,:,:)
    !- for quadr correlation
    ! chem1_g(nO,ng)%sc_p(:,:,:)= 12.4 + 0.01* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**2.
    !- for 4th correlation
    ! chem1_g(O3,ng)%sc_p(:,:,:)= 12.4 + 0.00025* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**4.
    !----------------------------------  squares   fim---------------------------
    !return
    !
    icx1=50
    jcx1=25
    icx2=50
    jcx2=75
    dummy= 100.*28./28.96
    !go to 333!<<<<<<<<<<<<<<<
    ! go to 555!<<<<<<<<<<<<<<<
    !----------------------------------
    !----------------------------------  gaussianas   ----------------------------------
    !b=1./90.
    ! b=1./150.
    !b=1./450.
    !b=1./280.
    do j=1,m3/2; do i=1,m2
       !gaussian
       !dum1 = (max( 0., exp(-b*( real(i-icx2)**2+ real(j-jcx2)**2 ) ) )*100)*28./28.96
       !dum2 = (max( 0., exp(-b*( real(i-icx1)**2+ real(j-jcx1)**2 ) ) )*100)*28./28.96
       !
       !if(dum1>0.01) chem1_g(CO,ng)%sc_p(:,i,j)=dum1
       !if(dum2>0.01)
       !chem1_g(CO,ng)%sc_p(:,i,j)=dum2
       !--- for sum conservation
       !dum1 = (max( 0., exp(-0.5*b*( real(i-icx1-10)**2+ real(j-jcx1-3)**2 ) ) )*280)*28./28.96
       !chem1_g(form,ng)%sc_p(:,i,j)=dum1
       !----------------------------

    enddo;enddo
    do j=m3/2,m3; do i=1,m2
       !gaussian
       !dum1 = (max( 0., exp(-b*( real(i-icx2)**2+ real(j-jcx2)**2 ) ) )*100)*28./28.96
       !dum2 = (max( 0., exp(-b*( real(i-icx1)**2+ real(j-jcx1)**2 ) ) )*100)*28./28.96
       !
       !chem1_g(CO,ng)%sc_p(:,i,j)=chem1_g(CO,ng)%sc_p(:,i,j-m3/2)
       !--- for sum conservation
       !chem1_g(form,ng)%sc_p(:,i,j)=chem1_g(form,ng)%sc_p(:,i,j-m3/2)
       !----------------------------
    enddo;enddo

    ! ass backgrd
    do j=1,m3; do i=1,m2;do k=1,m1

       !if(chem1_g(CO,ng)%sc_p(k,i,j)<20.*28./28.96)  chem1_g(CO,ng)%sc_p(k,i,j)=20.*28./28.96

       !--- for sum conservation
       !if(chem1_g(form,ng)%sc_p(k,i,j)<30.*28./28.96) chem1_g(form,ng)%sc_p(k,i,j)=30.*28./28.96
       !----------------------------

    enddo;enddo;enddo

    ! - two tracers
    !  chem1_g(CH4,ng)%sc_p(:,:,:)=120.*28./28.96 - chem1_g(CO,ng)%sc_p(:,:,:)
    ! 3  tracers
    !  chem1_g(NO2,ng)%sc_p(:,:,:)=360.*28./28.96 - chem1_g(form,ng)%sc_p(:,:,:)- chem1_g(CO,ng)%sc_p(:,:,:)
    ! 4  tracers
    ! chem1_g(PAN,ng)%sc_p(:,:,:)=40.*28./28.96 + 0.75*chem1_g(form,ng)%sc_p(:,:,:)+ 1.4* chem1_g(CO,ng)%sc_p(:,:,:)
    ! chem1_g(NO3,ng)%sc_p(:,:,:)=chem1_g(PAN,ng)%sc_p(:,:,:) + chem1_g(FORM,ng)%sc_p(:,:,:)+  chem1_g(CO,ng)%sc_p(:,:,:)
    !
    !- for linear correlations
    ! chem1_g(CH4,ng)%sc_p(:,:,:)= 22.4 + 1.53*chem1_g(CO,ng)%sc_p(:,:,:)
    !- for quadr correlation
    ! chem1_g(NO,ng)%sc_p(:,:,:)= 32.4 + 0.01* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**2.
    !- for 4th correlation
    ! chem1_g(O3,ng)%sc_p(:,:,:)= 12.4 + 0.5*0.000025* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**4.
    !----------------------------------  gaussianas fim  ----------------------------------

    return 

    !----------------------------------  triangulos      ----------------------------------
333 continue
    ! dois triangulos
    do j=1,m3; do i=1,m2
       !gaussian
       dum1 = (max( 0., 1. - (1./15.)*sqrt( real(i-icx2)**2+ real(j-jcx2)**2 ))*100)*28./28.96
       dum2 = (max( 0., 1. - (1./15.)*sqrt( real(i-icx1)**2+ real(j-jcx1)**2 ))*100)*28./28.96
       !
       !if(dum1>0.01) chem1_g(CO,ng)%sc_p(:,i,j)=dum1
       !if(dum2>0.01) chem1_g(CO,ng)%sc_p(:,i,j)=dum2
    enddo;enddo
    ! ass backgrd
    do j=1,m3; do i=1,m2;do k=1,m1

       !if(chem1_g(CO,ng)%sc_p(k,i,j)<20.*28./28.96)  chem1_g(CO,ng)%sc_p(k,i,j)=20.*28./28.96
    enddo;enddo;enddo
    !----------------------------------  triangulos fim  ----------------------------------
    !- for linear correlations
    !chem1_g(Ch4,ng)%sc_p(:,:,:)= 12.4 + 1.53*chem1_g(CO,ng)%sc_p(:,:,:)
    !- for quadr correlation
    !chem1_g(nO,ng)%sc_p(:,:,:)= 12.4 + 0.01* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**2.
    !- for 4th correlation
    !chem1_g(O3,ng)%sc_p(:,:,:)= 12.4 + 0.00025* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**4.


    return
    !    333 continue
    !--------------------------------------------------------------------
    ! dois circulos para slot
    do j=1,m3; do i=1,m2
       !if( sqrt( (real (i-icx1))**2+(real(j-jcx1))**2 ) .le. 15.) chem1_g(CO,ng)%sc_p(:,i,j)=dummy
       !if( sqrt( (real (i-icx2))**2+(real(j-jcx2))**2 ) .le. 15.) chem1_g(CO,ng)%sc_p(:,i,j)=dummy
       !chem1_g(CO,ng)%sc_p(5:20,10:30,10:30) = 100.*28./28.96
    enddo;enddo
    ! dois retangulos para o slot

    dummy= 20.*28./28.96
    do j=jcx1-2,jcx1+2; do i=icx1-2,icx1+15
       ! 1st retangle
       ! chem1_g(CO,ng)%sc_p(:,i,j)=dummy
    enddo;enddo
    do j=jcx2-2,jcx2+2; do i=icx2-15,icx2+2
       ! 1st retangle
       ! chem1_g(CO,ng)%sc_p(:,i,j)=dummy
    enddo;enddo

    !- for linear correlations
    !chem1_g(Ch4,ng)%sc_p(:,:,:)= 12.4 + 1.53*chem1_g(CO,ng)%sc_p(:,:,:)
    !- for quadr correlation
    !chem1_g(nO,ng)%sc_p(:,:,:)= 12.4 + 0.01* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**2.
    !- for 4th correlation
    !chem1_g(O3,ng)%sc_p(:,:,:)= 12.4 + 0.00025* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**4.

    !----------------------------------    slot  fim ----------------------

555 continue






    return

    !DO ispc=O3,O3

    ! DO k=1,m1
    !  dummy=0.;ncount=0
    !
    ! DO j=1,m3; DO i=1,m2
    !   dummy= dummy+  chem1_g(ispc,ng)%sc_p(k,i,j)
    !ncount=ncount+1
    ! enddo;enddo
    ! chem1_g(ispc,ng)%sc_p(k,:,:)    = dummy/ncount-10.
    !chem1_g(ispc,ng)%sc_p(k,:,:)=minval(chem1_g(ispc,ng)%sc_p(k,:,:))
    ! DO j=1,m3
    !  DO i=1,m2
    !   chem1_g(ispc,ng)%sc_p(k,i,j)=5.+minval(chem1_g(ispc,ng)%sc_p(k,2:m2-1,2:m3-1))
    !   if(i==2.and.j==2) print*,'o3=',o3,k,chem1_g(ispc,ng)%sc_p(k,i,j)
    ! DO j=7,13
    !  DO i=7,12
    !   DO k=3,20
    !       chem1_g(ispc,ng)%sc_p(k,i,j)=100.
    !ENDDO;ENDDO!;ENDDO;ENDDO

    !return

    do ispc=1,nx
       !  chem1_g(ispc,ng)%sc_p(:,:,:)=init_ajust(ispc)*chem1_g(ispc ,ng)%sc_p(:,:,:)
    enddo

    !stop 33
    !RETURN

    !
    !do i=1,m2
    ! do j=1,m3
    ! chem1_g(O3 ,ng)%sc_p(1:13,i,j)=10. + i*j/100
    !  chem1_g(no ,ng)%sc_p(1:20,i,j)=20. + i*j**2/1000
    !  chem1_g(no2 ,ng)%sc_p(1:10,i,j)=1. + (i**2) * j*2/1000
    !enddo
    ! enddo
    !chem1_g(H2O,ng)%sc_p(:,:,:)=basic_g(ng)%RV(:,:,:)*fcu

    !RETURN

    ! chem1_g(CO,ng)%sc_p(:,:,:) = 60.
    !chem1_g(O3,ng)%sc_p(:,:,:) = 50.
    !chem1_g(O3,ng)%sc_p(25:,:,:) = 100.

    !chem1_g(H2O2,ng)%sc_p(1:13,:,:) = 9.

    !aer1_g(accum,bburn,ng)%sc_p(18:20,5:m2-5,5:m3-5) = 300.
    !aer1_g(accum,bburn,ng)%sc_p(1:3,3:m2-3,3:m3-3) = 300.

    return
    !aer1_g(coarse,bburn,ng)%sc_p(1:10,3:m2-3,3:m3-3) = 10.
    !aer1_g(nucle,urban,ng)%sc_p(1:10,3:m2-3,3:m3-3) = 10.
    !aer1_g(accum,urban,ng)%sc_p(1:10,3:m2-3,3:m3-3) = 10.

    !aer1_g(accum,bburn,ng)%sc_p(17:19,3:m2-3,3:m3-3) = 8.
    !aer1_g(coarse,bburn,ng)%sc_p(17:19,3:m2-3,3:m3-3) = 8.
    !aer1_g(nucle,urban,ng)%sc_p(17:19,3:m2-3,3:m3-3) = 8.
    !aer1_g(accum,urban,ng)%sc_p(17:19,3:m2-3,3:m3-3) = 8.


  end subroutine initial_condition


  subroutine aer_background(ngrid,m1,m2,m3,ia,iz,ja,jz)
    integer,intent(IN) :: m1,m2,m3,ia,iz,ja,jz,ngrid
    real,parameter :: fcui=1.e-9     !de billion g [gas/part] /kg [ar] para kg/kg
    real, parameter :: aer_pbl = 10. , aer_ft=5. ! microg/m^3
    integer  m1pbl,i,j,k,ispc,imode
    real aer_pbl_ppb,aer_ft_ppb

    !-- we are not longer using this approach
    return

    if(AEROSOL == 0) return

    !#ifdef SIMPLE
    !-
    if(AEROSOL == 1) then
       ispc =urban
       imode=coarse
       if(spc_alloc_aer(transport,imode,ispc) /= ON) return

       m1pbl = min(10,m1)

       do j=ja,jz; do i=ia,iz
          ! pbl
          do k=1,m1pbl
             aer_pbl_ppb=aer_pbl!*1.e-9 & ! from ug/m3 to kg/m3
             ! /basic_g(ngrid)%dn0(k,i,j) & ! to kg/kg
             ! *1.e9 ! to ppbv

             !aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=max(aer_pbl_ppb,&
             !   aer1_g(imode,ispc,ngrid)%sc_p(k,i,j))
             aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=aer_pbl_ppb

          enddo
          !free-troposphere
          do k=m1pbl+1,15 !m1
             !
             aer_ft_ppb=aer_ft!*1.e-9 & ! from ug/m3 to kg/m3
             !/basic_g(ngrid)%dn0(k,i,j) & ! to kg/kg
             !*1.e9 ! to ppbv
             !aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=max(aer_ft_ppb,&
             !   aer1_g(imode,ispc,ngrid)%sc_p(k,i,j))
             aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=aer_ft_ppb
             !
          enddo
          do k=16, m1!m1
             !
             aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=0.
             !
          enddo

       enddo;enddo
    elseif(AEROSOL == 2) then
       !#elif MATRIX
       !---srf incluir instrucoes para o MATRIX

       !#endif
    endif
  end subroutine aer_background
end module ModChemistryDriver
