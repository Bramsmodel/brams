module ModDiffuse

  use ModDiffSclr, only: &
       diffsclr_brams31
  
  use ModTurbKE, only: &
       tkescl, &
       tkeeps, &
       mxtked, &
       tkemy
  
  use ModTurbK, only: &
       strain, &
       bruvais, &
       mxdefm, &
       klbnd

  use ModTurbDiff, only: &
       diffvel
  
  use ModScalarTable, only:  &
       ScalarTable

  use ModBasicFields, only: &
       BasicFields
  
  use mem_tend,  only:    &
       tend                   ! %tket, %epst, %ut, %vt, %wt

  use ModTurbFields, only:     &
       TurbFields


  use mem_grid, only:     &
       jdim,              &   !        INTENT(IN)
       ngrid,             &   !        INTENT(IN)
       grid_g,            &   ! %rtgt  INTENT(IN)
       npatch,            &   !        INTENT(IN)
       zt,                &   !        INTENT(IN)
       nstbot,            &   !        INTENT(IN)
       nnzp

  use mem_leaf, only:     &
       leaf_g                 ! %ustar, %patch_area

  use ModMicroFields, only:    &
       MicroFields

  use node_mod, only:     &
       mxp,     &     !INTENT(IN)
       myp,     &     !INTENT(IN)
       mzp,     &     !INTENT(IN)
       ia,      &     !INTENT(IN)
       iz,      &     !INTENT(IN)
       ja,      &     !INTENT(IN)
       jz,      &     !INTENT(IN)
       ia_1,    &     !INTENT(IN)
       ja_1,    &     !INTENT(IN)
       iz1,     &     !INTENT(IN)
       jz1,     &     !INTENT(IN)
       ibcon,   &     !INTENT(IN)
       mynum,   &     !INTENT(IN)
       nodei0,  &     !INTENT(IN)
       nodej0,  &     !INTENT(IN)
       nodemyp, &     !INTENT(IN)
       nodemxp, &     !INTENT(IN)
       ia1,     &     !INTENT(IN)
       ja1,     &     !INTENT(IN)
       iz_1,    &     !INTENT(IN)
       jz_1,    &     !INTENT(IN)
       izu,     &     !INTENT(IN)
       jzv,     &     !INTENT(IN)
       mynum


  use ke_coms, only:      &
       alf_eps,   &
       alf_tke


  use mem_opt, only:      &      ! For optmization
       jstep,   &
       istep,   &
       opt                       ! %ind1_x_a,

  use ModNamelistFile, only: &
       NamelistFile

  use ModMicControl, only: &
       MicControl

  implicit none
  
  private
  public :: diffuse_brams31

contains



  subroutine diffuse_brams31(oneScalarTab, oneScalarTabSize, oneBasicFields, &
       oneNamelistFile, oneTurbFields, gridId, oneMicControl, oneMicroFields)

    ! +-----------------------------------------------------------------+
    ! \     this routine is the subdriver to compute tendencies due to  \
    ! \       subgrid-scale turbulence.                                 \
    ! +-----------------------------------------------------------------+

    !**(JP)** uses 1D arrays and scalar pointers representing 3D arrays;
    !**(JP)** requires full rewriting, affecting diffscalar_brams31, that is
    !**(JP)** called only here

    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(inout) :: oneScalarTabSize
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(TurbFields), pointer, intent(in) :: oneTurbFields
    integer, intent(in) :: gridId
    type(MicControl), pointer, intent(in) :: oneMicControl
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    
    include "constants.h"

    integer(kind=i8) :: mxyzp, ind
    integer :: n
    real :: s1,s2,s3


!!$    real, pointer :: vkh_p,hkh_p
    real, pointer :: vkh_p(:)
    real, pointer :: hkh_p(:)


    integer :: i,j,k,ksf
    integer :: idiffk
    real :: xkhkm

    !! For Optimization
    integer      :: htint_i, htint_j, iia, jja, iiz, jjz, iistep, jjstep
    !! ALF
    !
    ! (JP) removing scr2 from scratch
    real, target :: scr1(mxp*myp*mzp)
    real, target :: scr2(mxp*myp*mzp)
    real, target :: vt2da(mxp*myp)
    real, target :: vt3dp(mxp,myp,mzp)
    real, target :: vt3dn(mxp,myp,mzp)
    real, target :: vt3do(mxp,myp,mzp)
    real, target :: vt3dm(mxp,myp,mzp)
    real, target :: vt3dl(mxp,myp,mzp)
    real, target :: vt3dk(mxp,myp,mzp)
    real, target :: vt3dj(mxp,myp,mzp)
    real, target :: vt3di(mxp*myp*mzp)
    real, target :: vt3dh(mxp*myp*mzp)
    real, target :: vt3dg(mxp,myp,mzp)
    real, target :: vt3df(mxp,myp,mzp)
    real, target :: vt3de(mxp,myp,mzp)
    real, target :: vt3dd(mxp,myp,mzp)
    real, target :: vt3dc(mxp,myp,mzp)
    real, target :: vt3db(mxp,myp,mzp)
    real, target :: vt3da(mxp,myp,mzp)

    character(len=*), parameter :: h="**(diffuse_brams31)**" 
    real :: vctr1(mzp)
    real :: vctr2(mzp)
    real :: vctr3(mzp)
    real :: vctr34(mzp)
    
    idiffk=oneNamelistFile%idiffk(gridId)
    xkhkm=oneNamelistFile%xkhkm(gridId)
    
    mxyzp = mxp*myp*mzp

    ! Shaved-Eta Coordinate not available in this optimization method


    call strain(mzp,mxp,myp,ia,iz,ja,jz                       &
         ,ia_1,ja_1,iz1,jz1,jdim                                &
         ,oneBasicFields%up ,oneBasicFields%vp  &
         ,oneBasicFields%wp ,vt3da      &
         ,vt3db     ,vt3dc      &
         ,vt3dd     ,vt3de      &
         ,vt3df     ,vt3dg      &
         ,vt3dh     ,vt3di      &
         ,vt3dn     ,scr2           &
         ,idiffk)

    if (oneMicControl%level<=1) &
         vt3dp = 0.

    if (oneMicControl%level>=2) &
         call ae1_l(mxyzp, vt3dp, oneMicroFields%rcp)


    call bruvais(mzp,mxp,myp,ia,iz,ja,jz                          &
         ,oneBasicFields%theta ,oneBasicFields%rtp  &
         ,oneBasicFields%rv    ,vt3dp       &
         ,oneBasicFields%pp    ,oneBasicFields%pi0  &
         ,vt3dj        ,grid_g(ngrid)%rtgt  &
         ,grid_g(ngrid)%lpw, &
         oneMicControl)

    if (idiffk <= 3) then
       call mxdefm(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim            &
            ,vt3dh      ,vt3di      &
            ,vt3dj      ,vt3dk      &
            ,scr1       ,scr2           &
            ,oneBasicFields%dn0 ,grid_g(ngrid)%rtgt &
            ,grid_g(ngrid)%dxt  ,grid_g(ngrid)%dyt  &
            ,grid_g(ngrid)%lpw  ,mynum  &
            ,oneNamelistFile, gridId)
    endif

    if (idiffk == 1) then
       call tkemy(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim,nodei0(mynum,ngrid),nodej0(mynum,ngrid)  &
            ,oneTurbFields%tkep   ,tend%tket            &
            ,vt3dh        ,vt3di        &
            ,vt3dj        ,scr1         &
            ,grid_g(ngrid)%rtgt   ,oneBasicFields%theta &
            ,oneBasicFields%dn0   ,oneBasicFields%up    &
            ,oneBasicFields%vp    ,oneBasicFields%wp    &
            ,oneTurbFields%sflux_u   ,oneTurbFields%sflux_v    &
            ,oneTurbFields%sflux_w   ,oneTurbFields%sflux_t,vctr34 &
            ,grid_g(ngrid)%lpw       ,grid_g(ngrid)%lpu       &
            ,grid_g(ngrid)%lpv    )
    endif

    if (idiffk == 4) then
       call mxtked(mzp,mxp,myp,ia,iz,ja,jz  &
            ,ibcon,jdim  &
            ,oneTurbFields%tkep   ,tend%tket            &
            ,oneBasicFields%up    ,oneBasicFields%vp    &
            ,oneBasicFields%wp    ,oneBasicFields%rtp   &
            ,oneBasicFields%rv    ,oneBasicFields%theta &
            ,vt3da        ,vt3dc        &
            ,vt3dh        ,vt3dj        &
            ,scr1         ,scr2             &
            ,oneTurbFields%sflux_u   ,oneTurbFields%sflux_v    &
            ,oneTurbFields%sflux_w   ,oneTurbFields%sflux_t    &
            ,grid_g(ngrid)%dxt       ,grid_g(ngrid)%rtgt       &
            ,grid_g(ngrid)%lpw       )
    endif

    !_STC............................................................
    !_STC Call to subroutine tkescl for E-l closure
    !_STC (S. Trini Castelli)
    !_STC............................................................
    if (idiffk == 5) then
       call tkescl(mzp,mxp,myp,npatch,ia,iz,ja,jz  &
            ,oneTurbFields%tkep,tend%tket &
            ,oneTurbFields%epsp,tend%epst &
            ,vt3da,vt3dc  &
            ,vt3dh,vt3di  &
            ,vt3dj,scr1  &
            ,scr2 ,grid_g(ngrid)%rtgt  &
            ,vt3dd,vt3de,grid_g(ngrid)%dxt  &
            ,leaf_g(ngrid)%ustar,leaf_g(ngrid)%patch_area &
            ,grid_g(ngrid)%lpw,oneBasicFields%dn0  )
    endif
    !_STC............................................................
    !_STC Call to subroutine tkeeps for E-eps closure
    !_STC (S. Trini Castelli)
    !_STC............................................................
    if (idiffk == 6) then
       call tkeeps(mzp,mxp,myp,npatch,ia,iz,ja,jz  &
            ,oneTurbFields%tkep,tend%tket  &
            ,oneTurbFields%epsp,tend%epst  &
            ,vt3da,vt3dc  &
            ,vt3dh,vt3di  &
            ,vt3dj,scr1  &
            ,scr2 ,grid_g(ngrid)%rtgt  &
            ,leaf_g(ngrid)%ustar,leaf_g(ngrid)%patch_area &
            ,grid_g(ngrid)%lpw,oneBasicFields%dn0  )
    endif
    !_STC..................................................
    !_STC    Note: from subroutines TKESCL, TKEEPS :
    !_STC           VT3DI=Ke
    !_STC           SCR1=Km
    !_STC           VT3DH = Kh
    !_STC           SCR2 = SCR1 = Km
    !_STC..................................................
    !_STC............................................................

    call klbnd(mzp,mxp,myp,ibcon,jdim  &
         ,scr1 ,oneBasicFields%dn0,grid_g(ngrid)%lpw)
    call klbnd(mzp,mxp,myp,ibcon,jdim  &
         ,scr2 ,oneBasicFields%dn0,grid_g(ngrid)%lpw)
    call klbnd(mzp,mxp,myp,ibcon,jdim  &
         ,vt3dh,oneBasicFields%dn0,grid_g(ngrid)%lpw)
    !_STC ....... boundary conditions even on Ke diffusion coefficient
    if(idiffk ==  5 .or. idiffk == 6) &
         call klbnd(mzp,mxp,myp,ibcon,jdim  &
         ,vt3di,oneBasicFields%dn0,grid_g(ngrid)%lpw)

    !bob  swap new hkm, vkm, and vkh with past time level:  lagged K's have
    !bob  internal lateral boundary values from neighboring nodes

    ind = 0
    do j = 1,nodemyp(mynum,ngrid)
       do i = 1,nodemxp(mynum,ngrid)
          do k = 1,nnzp(ngrid)
             ind = ind + 1
             s1 = scr2(ind)
             s2 = scr1(ind)
             s3 = vt3dh(ind)
             scr2(ind) = oneTurbFields%hkm(k,i,j)
             scr1(ind) = oneTurbFields%vkm(k,i,j)
             vt3dh(ind) = oneTurbFields%vkh(k,i,j)
             !! also for vt3di = K(tke) ?????    22 March 02
             !!         vt3di(ind) = oneTurbFields%vke(k,i,j)
             oneTurbFields%hkm(k,i,j) = s1
             oneTurbFields%vkm(k,i,j) = s2
             oneTurbFields%vkh(k,i,j) = s3
          enddo
       enddo
    enddo

    call diffvel(mzp,mxp,myp,ia,iz,ja,jz,jdim,ia_1,ja_1             &
         ,ia1,ja1,iz_1,jz_1,iz1,jz1,izu,jzv,idiffk             &
         ,oneBasicFields%up     ,oneBasicFields%vp      &
         ,oneBasicFields%wp     ,tend%ut                &
         ,tend%vt               ,tend%wt                &
         ,vt3da         ,vt3db          &
         ,vt3dc         ,vt3dd          &
         ,vt3de         ,vt3df          &
         ,vt3dg         ,vt3dj          &
         ,vt3dk         ,vt3dl          &
         ,vt3dm         ,vt3dn          &
         ,vt3do         ,grid_g(ngrid)%rtgu     &
         ,grid_g(ngrid)%rtgv    ,grid_g(ngrid)%rtgt     &
         ,oneTurbFields%sflux_u ,oneTurbFields%sflux_v  &
         ,oneTurbFields%sflux_w ,oneBasicFields%dn0     &
         ,oneBasicFields%dn0u   ,oneBasicFields%dn0v    &
         ,scr1          ,scr2                   &
         ,ibcon,mynum,oneNamelistFile%ihorgrad)

    ! Convert momentum K's to scalar K's, if necessary

    if (idiffk <= 3) then
       do ind = 1,mxyzp
          scr2(ind) = scr2(ind) * xkhkm
       enddo
    elseif (idiffk == 4) then
       do ind = 1,mxyzp
          vt3di(ind) = 2. * scr1(ind)
       enddo
    endif


    ! Calculating Index and Weights for interpolation (htint)

    ! Setting indexes and Weights for HTINT with Optmization

    ! Loop for a block of spatial local region

    jjstep = min(jstep,(jz-ja))
    iistep = min(istep,(iz-ia))

    do jja = ja, jz, jjstep

       jjz = min(jja+jjstep-1, jz) ! Calculating last element

       do iia = ia, iz, iistep

          iiz = min(iia+iistep-1, iz) ! Calculating last element

          ! Compute the Indexes and Weights for x_dir - ALF
          do j=jja,jjz
             do i=iia,iiz
                do k=1,mzp !m1
                   vctr1(k)=grid_g(ngrid)%topt(i-1,j)+  &
                        zt(k)*grid_g(ngrid)%rtgt(i-1,j)

                   vctr2(k)=grid_g(ngrid)%topt(i  ,j)+  &
                        zt(k)*grid_g(ngrid)%rtgt(i ,j)

                   vctr3(k)=grid_g(ngrid)%topt(i+1,j)+  &
                        zt(k)*grid_g(ngrid)%rtgt(i+1,j)
                enddo

                htint_i = i-iia+1
                htint_j = j-jja+1

                call htint_index(mzp, vctr1, mzp, vctr2,  &
                     opt%ind1_x_a(1,htint_i,htint_j),   &
                     opt%ind2_x_a(1,htint_i,htint_j),   &
                     opt%weight_x_a(1,htint_i,htint_j))

                call htint_index(mzp, vctr3, mzp, vctr2,  &
                     opt%ind1_x_b(1,htint_i,htint_j),   &
                     opt%ind2_x_b(1,htint_i,htint_j),   &
                     opt%weight_x_b(1,htint_i,htint_j))


                ! Compute the Indexes and Weights for y_dir - ALF

                do k=1,mzp !m1
                   vctr1(k)=grid_g(ngrid)%topt(i,j-jdim)+  &
                        zt(k)*grid_g(ngrid)%rtgt(i,j-jdim)

                   vctr3(k)=grid_g(ngrid)%topt(i,j+jdim)+  &
                        zt(k)*grid_g(ngrid)%rtgt(i,j+jdim)
                enddo

                call htint_index(mzp, vctr1, mzp, vctr2,  &
                     opt%ind1_y_a(1,htint_i,htint_j),   &
                     opt%ind2_y_a(1,htint_i,htint_j),   &
                     opt%weight_y_a(1,htint_i,htint_j))

                call htint_index(mzp, vctr3, mzp, vctr2,  &
                     opt%ind1_y_b(1,htint_i,htint_j),   &
                     opt%ind2_y_b(1,htint_i,htint_j),   &
                     opt%weight_y_b(1,htint_i,htint_j))

             enddo
          enddo
          ! ALF

          !! End of calculating Index and Weights for interpolation (htint)


          do n = 1, oneScalarTabSize

             vt2da = 0.
             if (nstbot == 1) then
                if (oneScalarTab(n)%name == 'THP' .or. &
                     oneScalarTab(n)%name == 'THC' ) then
                   call atob(mxp*myp, oneTurbFields%sflux_t, vt2da)
                elseif (oneScalarTab(n)%name == 'RTP') then
                   call atob(mxp*myp, oneTurbFields%sflux_r, vt2da)
                endif
             endif

             ! 3/10/01 - Define ksf below, the "K scalar flag", to let subroutine diffsclr
             ! know which vertical K is being passed to it.  If diffsclr sees that it's
             ! a different K from the previous one, diffsclr will re-compute the tridiff
             ! matrix coefficients.  In order to use vertical scalar K's other than
             ! vt3dh and vt3di, use ksf = 3, ksf = 4, etc. for each different K.

             !_STC..................................................
             !_STC Corrections to account for the new idiffk options
             !_STC for E-l and E-eps closure. Isotropy hypothesis.
             !_STC (S. Trini Castelli)
             !_STC..................................................

!!$             if (oneScalarTab(n)%name == 'TKEP') then
!!$                vkh_p => vt3di(1)
!!$                hkh_p => scr2(1)
!!$                if (idiffk >= 4) hkh_p => vt3di(1)
!!$                ksf = 1   
!!$             elseif (oneScalarTab(n)%name == 'EPSP') then
!!$                vkh_p => vt3di(1)
!!$                hkh_p => scr2(1)
!!$                if (idiffk >= 4)  hkh_p => vt3di(1)
!!$                ksf = 3
!!$                ! Convert Ktke to Keps; it will be converted back after use below
!!$                call ae1t0_l(mxyzp, vkh_p, vkh_p, (ALF_EPS/ALF_TKE))
!!$                call ae1t0_l(mxyzp, hkh_p, hkh_p, (ALF_EPS/ALF_TKE))
!!$             else
!!$                vkh_p => vt3dh(1)
!!$                hkh_p => scr2(1)
!!$                if (idiffk >= 4) hkh_p => vt3dh(1)
!!$                ksf = 2
!!$             endif


             if (oneScalarTab(n)%name == 'TKEP') then
                vkh_p => vt3di
                hkh_p => scr2
                if (idiffk >= 4) hkh_p => vt3di
                ksf = 1   
             elseif (oneScalarTab(n)%name == 'EPSP') then
                vkh_p => vt3di
                hkh_p => scr2
                if (idiffk >= 4)  hkh_p => vt3di
                ksf = 3
                ! Convert Ktke to Keps; it will be converted back after use below
                call ae1t0_l(mxyzp, vkh_p, vkh_p, (ALF_EPS/ALF_TKE))
                call ae1t0_l(mxyzp, hkh_p, hkh_p, (ALF_EPS/ALF_TKE))
             else
                vkh_p => vt3dh
                hkh_p => scr2
                if (idiffk >= 4) hkh_p => vt3dh
                ksf = 2
             endif


             call diffsclr_brams31(mzp,mxp,myp,iia,iiz,jja,jjz,jdim        &
                  ,ia_1,ja_1,ia1,ja1,iz_1,jz_1,iz1,jz1,n,ksf               &
                  ,oneScalarTab(n)%var_p_3D,oneScalarTab(n)%var_t_1D,vt3da            &
                  ,vt3db           ,vt3df            &
                  ,vt3dg           ,vt3dj            &
                  ,vt3dk           ,vt3do            &
                  ,vt3dc           ,vt3dd            &
                  ,vt3dl           ,vt3dm            &
                  ,grid_g(ngrid)%rtgt     &
                  ,vt2da           ,oneBasicFields%dn0   &
                  ,vkh_p                   ,hkh_p                       )

             if (oneScalarTab(n)%name == 'EPSP') then
                call ae1t0_l(mxyzp, vkh_p, vkh_p, (ALF_TKE/ALF_EPS))
                call ae1t0_l(mxyzp, hkh_p, hkh_p, (ALF_TKE/ALF_EPS))
             endif

          enddo

       enddo

    enddo
  end subroutine diffuse_brams31
end module ModDiffuse
