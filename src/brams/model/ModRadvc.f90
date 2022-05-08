!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModRadvc
  use iso_fortran_env, only: &
       real64
  
  use ModBasicFields, only: &
       BasicFields, &
       DeepCopyToBasicFields, &
       DeepCopyFromBasicFields
  
  use ModMonotonicAdvection, only: &
       advmnt 

  use ModNamelistFile, only: &
       namelistfile

  use ModParallelEnvironment, only: &
       MsgDump

  use mem_grid, only: &
       ngrid, &
       grid_g, &
       hw4, &
       dzt, &
       zm, &
       zt, &
       itopo, &
       jdim, &
       dzm, &
       dtlt, &
       if_adap, &
       time, &
       dyncore_flag

  use mem_basic, only: &
       basic_g

  use mem_scratch, only: &
       scratch, &
       vctr1, &
       vctr2, &
       vctr3

  use ModScalarTable, only: &
       ScalarTable

  use grid_dims, only: &
       maxgrds, &
       nzpmax

  use mem_tend, only: &
       tend

  use ccatt_start, only: &
       ccatt      ! intent(in)

  use mem_aer1, only: &
       aerosol,    &
       num_scalar_aer_1st

  use mem_chem1, only: &
       nspecies_transported

  use module_dry_dep, only: &
       dd_sedim,         &
       fa_preptc_with_sedim

  implicit none
  
  private
  public :: advectc
  public :: advtndc
  public :: fa_preptc
  public :: fa_xc
  public :: fa_yc
  public :: fa_zc
contains



  
  subroutine advectc(oneScalarTab, oneScalarTabSize, oneBasic, oneAveBasic, &
       varn,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
    !> @brief: advectc
    !! @date:  18/Nov/2015
    !! @version:  5.2
    !! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
    !! 
    !! @copyright Under CC-GPL License by INPE/CPTEC
    !! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt

    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(inout) :: oneScalarTabSize
    type(BasicFields), pointer, intent(in) :: oneBasic
    type(BasicFields), pointer, intent(in) :: oneAveBasic
    integer :: mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,n
    integer(kind=real64) :: mxyzp
    character(len=*) :: varn

    real :: dtlto2
    real, dimension(maxgrds), save :: save_dtlt
    integer :: i,j,k,ind

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(advectc)**"

    integer :: i_scl

    integer, dimension(maxgrds), save :: ncall
    data ncall/maxgrds*0/

    if (ncall(ngrid) == 0 .or. dtlt .ne. save_dtlt(ngrid) ) then
       ncall(ngrid) = 1
       save_dtlt(ngrid) = dtlt
    endif
    mxyzp = mxp * myp * mzp

    if (trim(varn) .eq. 'V' .or. trim(varn) .eq. 'ALL') then
       ! Advect  U, V, and W

       if (if_adap == 0) then

          if ( (dyncore_flag == 0) .or. (dyncore_flag == 1) ) then 

             call vel_advectc(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv  &
                  ,tend%ut                    &
                  ,tend%vt                   ,tend%wt                    &
                  ,scratch%vt3da              &
                  ,scratch%vt3db             ,scratch%vt3dc, &
                  oneBasic)

          else if (dyncore_flag == 2) then 
             !MB: this branch must probably be transferred also to the other if-statements in this subroutine!
             print*, "velocity advection"
             call vel_advectc(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv  &
                  ,tend%ut_rk                 &
                  ,tend%vt_rk                ,tend%wt_rk                 &
                  ,scratch%vt3da              &
                  ,scratch%vt3db             ,scratch%vt3dc, &
                  oneBasic)
          end if

       else

          call DeepCopyToBasicFields(oneBasic, oneAveBasic)
          call vel_advectc_adap(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,jdim    &
               ,grid_g(ngrid)%lpu    ,grid_g(ngrid)%lpv       &
               ,grid_g(ngrid)%lpw    ,oneBasic%uc     &
               ,oneBasic%vc    ,oneBasic%wc     &
               ,tend%ut              ,tend%vt               &
               ,tend%wt              ,oneBasic%dn0    &
               ,oneBasic%dn0u  ,oneBasic%dn0v   &
               ,grid_g(ngrid)%aru    ,grid_g(ngrid)%arv     &
               ,grid_g(ngrid)%arw    ,grid_g(ngrid)%volu    &
               ,grid_g(ngrid)%volv   ,grid_g(ngrid)%volw    &
               ,scratch%vt3da        ,scratch%vt3db         &
               ,scratch%vt3dc        ,scratch%vt3dd         &
               ,scratch%vt3de        ,scratch%vt3df        ,time)
          call DeepCopyFromBasicFields(oneBasic, oneAveBasic)

       endif

    endif


    if (trim(varn) .eq. 'T' .or. trim(varn) .eq. 'SCALAR') then

       ! Advect  scalars

       call DeepCopyToBasicFields(oneBasic, oneAveBasic)
       dtlto2 = .5 * dtlt
       ind = 0
       do j = 1,myp
          do i = 1,mxp
             do k = 1,mzp
                ind = ind + 1
                scratch%vt3da(ind) = (oneBasic%up(k,i,j)  &
                     + oneBasic%uc(k,i,j)) * dtlto2
                scratch%vt3db(ind) = (oneBasic%vp(k,i,j)  &
                     + oneBasic%vc(k,i,j)) * dtlto2
                scratch%vt3dc(ind) = (oneBasic%wp(k,i,j)  &
                     + oneBasic%wc(k,i,j)) * dtlto2
             enddo
          enddo
       enddo
       call DeepCopyFromBasicFields(oneBasic, oneAveBasic)

       if (if_adap == 0) then

          call fa_preptc(mzp,mxp,myp        &
               ,scratch%vt3da             ,scratch%vt3db              &
               ,scratch%vt3dc             ,scratch%vt3dd              &
               ,scratch%vt3de             ,scratch%vt3df              &
               ,scratch%vt3dh             ,scratch%vt3di              &
               ,scratch%vt3dj             ,scratch%vt3dk              &
               ,mynum, &
               oneBasic, oneAveBasic)

       else

          call DeepCopyToBasicFields(oneBasic, oneAveBasic)
          call fa_preptc_adap(mzp,mxp,myp                               &
               ,scratch%vt3da            ,scratch%vt3db             &
               ,scratch%vt3dc            ,scratch%vt3dd             &
               ,scratch%vt3de            ,scratch%vt3df             &
               ,scratch%vt3dh            ,oneBasic%dn0    &
               ,oneBasic%dn0u ,oneBasic%dn0v   &
               ,grid_g(ngrid)%aru   ,grid_g(ngrid)%arv     &
               ,grid_g(ngrid)%arw   ,grid_g(ngrid)%volt    &
               ,grid_g(ngrid)%dxu   ,grid_g(ngrid)%dyv     &
               ,grid_g(ngrid)%dxt   ,grid_g(ngrid)%dyt     &
               ,zt,zm,dzm,vctr1,vctr2,jdim,mynum                          )
          call DeepCopyFromBasicFields(oneBasic, oneAveBasic)

       endif

       !- combine old adv for thermo + micro with mnt adv for tracers
       if( advmnt== 0) then
          i_scl=oneScalarTabSize  !- all scalars 
       elseif(advmnt == 2) then
          i_scl= oneScalarTabSize - NSPECIES_TRANSPORTED !- only theta_il+water+tke
       elseif(advmnt == 3) then
          i_scl=1  !- only theta_il
       endif

       do n=1,i_scl

          !- if RK or ABM3 schemes, THP/THC are not transported here
          if (dyncore_flag == 2) then
             if (oneScalarTab(n)%name == 'THC' .or. &
                  oneScalarTab(n)%name == 'THP') cycle
          endif

          if(ccatt == 1 .and. aerosol == 1) then
             if( n >= num_scalar_aer_1st ) then 
                if (if_adap == 0) then
                   call DeepCopyToBasicFields(oneBasic, oneAveBasic)
                   call fa_preptc_with_sedim(mzp,mxp,myp    &
                        ,scratch%vt3da             ,scratch%vt3db      &
                        ,scratch%vt3dc              ,scratch%vt3df &
                        ,scratch%vt3dk      &
                        ,oneBasic%dn0        ,grid_g(ngrid)%rtgt       &
                        ,grid_g(ngrid)%f13t        ,grid_g(ngrid)%f23t       &
                                !srf- aerosol section 
                        ,dtlt      &
                        ,N    & ! current scalar being transported
                        ,num_scalar_aer_1st    & ! 1st aerosol at scalar table
                        ,oneBasic%wp         & ! air vertical velocity (P time)
                        ,oneBasic%wc         & ! air vertical velocity (C time)
                        ,scratch%vt3dp             & ! ) ! to save horizontal contribution on the sigmaz velocity
                        ,nzpmax,hw4,dzm,dzt,dd_sedim(:,ngrid))   ! (DMK) deposicao seca (cod. limpo)
                   call DeepCopyFromBasicFields(oneBasic, oneAveBasic)
                else
                   print*,'sedim not yet prepared for shaved eta'
                   stop 3333
                endif
             endif
          endif
          !--(DMK-CCATT-FIM)-----------------------------------------------------

          call atob_long(mxyzp, oneScalarTab(n)%var_p_3D, scratch%scr1)

          if (if_adap == 0) then

             call fa_xc(mzp,mxp,myp,ia,iz,1,myp        &
                  ,oneScalarTab(n)%var_p_3D ,scratch%scr1   &
                  ,scratch%vt3da ,scratch%vt3dd  &
                  ,scratch%vt3dg ,scratch%vt3dh  &
                  ,scratch%vt3di ,mynum              )

             if (jdim .eq. 1)  &
                  call fa_yc(mzp,mxp,myp,ia,iz,ja,jz        &
                  ,oneScalarTab(n)%var_p_3D ,scratch%scr1    &
                  ,scratch%vt3db  ,scratch%vt3de   &
                  ,scratch%vt3dg  ,scratch%vt3dj   &
                  ,scratch%vt3di  ,jdim,mynum         )

             call fa_zc(mzp,mxp,myp,ia,iz,ja,jz        &
                  ,oneScalarTab(n)%var_p_3D ,scratch%scr1    &
                  ,scratch%vt3dc  ,scratch%vt3df   &
                  ,scratch%vt3dg  ,scratch%vt3dk   &
                  ,vctr1,vctr2,mynum                     )


             if (dumpLocal) then
                call MsgDump(h//" invokes advtndc")
             end if

             call advtndc(mzp,mxp,myp,ia,iz,ja,jz    &
                  ,oneScalarTab(n)%var_p_3D ,scratch%scr1  &
                  ,oneScalarTab(n)%var_t_1D ,dtlt,mynum        )

          else

             call fa_xc_adap(mzp,mxp,myp,ia,iz,1,myp         &
                  ,grid_g(ngrid)%lpw ,oneScalarTab(n)%var_p_3D            &
                  ,scratch%scr1      ,scratch%vt3da  &
                  ,scratch%vt3dd     ,scratch%vt3dg  &
                  ,scratch%vt3dh     ,mynum              )

             if (jdim .eq. 1)                                &
                  call fa_yc_adap(mzp,mxp,myp,ia,iz,ja,jz         &
                  ,grid_g(ngrid)%lpw ,oneScalarTab(n)%var_p_3D   &
                  ,scratch%scr1      ,scratch%vt3db   &
                  ,scratch%vt3de     ,scratch%vt3dg   &
                  ,scratch%vt3dh   &
                  ,jdim                    ,mynum              )

             call fa_zc_adap(mzp,mxp,myp,ia,iz,ja,jz         &
                  ,grid_g(ngrid)%lpw ,oneScalarTab(n)%var_p_3D  &
                  ,scratch%scr1      ,scratch%vt3dc   &
                  ,scratch%vt3df     ,scratch%vt3dg   &
                  ,scratch%vt3dh     ,vctr1              &
                  ,vctr2                   ,mynum              )

             call advtndc_adap(mzp,mxp,myp,ia,iz,ja,jz  &
                  ,grid_g(ngrid)%lpw ,oneScalarTab(n)%var_p_3D   &
                  ,scratch%scr1      ,oneScalarTab(n)%var_t_1D   &
                  ,dtlt                    ,mynum         )

          endif

       enddo

    endif
  end subroutine advectc

  !     *********************************************************************



  subroutine vel_advectc(m1,m2,m3,ia,iz,ja,jz,izu,jzv  &
       ,ut,vt,wt,flxu,flxv,flxw, oneBasic)
    !> @brief: vel_advectc
    !! @date:  18/Nov/2015
    !! @version:  5.2
    !! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
    !! 
    !! @copyright Under CC-GPL License by INPE/CPTEC
    !! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt

    integer,intent(in) :: m1,m2,m3,ia,iz,ja,jz,izu,jzv
    real, dimension(m1,m2,m3) :: ut,vt,wt,flxu,flxv,flxw
    type(BasicFields), pointer, intent(in) :: oneBasic


    integer :: j,i,k,jm,im
    real :: c1z,c1x,c1y

    ! Compute momentum fluxes flxu, flxv, flxw

    do j = 1,m3
       do i = 1,m2
          do k = 1,m1
             flxu(k,i,j) = oneBasic%uc(k,i,j) * oneBasic%dn0u(k,i,j) * grid_g(ngrid)%rtgu(i,j)  &
                  * grid_g(ngrid)%fmapui(i,j)
             flxv(k,i,j) = oneBasic%vc(k,i,j) * oneBasic%dn0v(k,i,j) * grid_g(ngrid)%rtgv(i,j)  &
                  * grid_g(ngrid)%fmapvi(i,j)
          enddo
       enddo
    enddo

    if(itopo.eq.0) then
       do j = 1,m3
          do i = 1,m2
             do k = 1,m1-1
                flxw(k,i,j) = oneBasic%wc(k,i,j)  &
                     * .5 * (oneBasic%dn0(k,i,j) + oneBasic%dn0(k+1,i,j))
             enddo
          enddo
       enddo
    else
       do j = 1,m3
          jm = max(j-1,1)
          do i = 1,m2
             im = max(i-1,1)
             do k = 1,m1-1
                flxw(k,i,j) = oneBasic%wc(k,i,j)  &
                     * .5 * (oneBasic%dn0(k,i,j) + oneBasic%dn0(k+1,i,j))  &
                     + hw4(k) * ((flxu(k,i,j) + flxu(k+1,i,j)  &
                     + flxu(k,im,j) + flxu(k+1,im,j)) * grid_g(ngrid)%f13t(i,j)  &
                     + (flxv(k,i,j) + flxv(k+1,i,j)  &
                     + flxv(k,i,jm) + flxv(k+1,i,jm)) * grid_g(ngrid)%f23t(i,j))
             enddo
          enddo
       enddo
    endif

    ! Compute advection contribution to U tendency

    do j = ja,jz
       do i = ia,izu
          c1z = .25 / grid_g(ngrid)%rtgu(i,j)
          c1x = c1z * grid_g(ngrid)%fmapu(i,j) * grid_g(ngrid)%dxu(i,j)
          c1y = c1z * grid_g(ngrid)%fmapu(i,j) * grid_g(ngrid)%dyu(i,j)

          do k = 2,m1-1
             ut(k,i,j) = ut(k,i,j) + c1x / oneBasic%dn0u(k,i,j) * (  &
                  (flxu(k,i,j) + flxu(k,i-1,j))  &
                  * (oneBasic%uc(k,i,j) + oneBasic%uc(k,i-1,j))  &
                  - (flxu(k,i,j) + flxu(k,i+1,j))  &
                  * (oneBasic%uc(k,i,j) + oneBasic%uc(k,i+1,j))  &
                  + (flxu(k,i+1,j) - flxu(k,i-1,j)) * 2.* oneBasic%uc(k,i,j) )
          enddo

          do k = 2,m1-1
             ut(k,i,j) = ut(k,i,j) + c1y / oneBasic%dn0u(k,i,j) * (  &
                  (flxv(k,i,j-jdim) + flxv(k,i+1,j-jdim))  &
                  * (oneBasic%uc(k,i,j) + oneBasic%uc(k,i,j-jdim))  &
                  - (flxv(k,i,j) + flxv(k,i+1,j))  &
                  * (oneBasic%uc(k,i,j) + oneBasic%uc(k,i,j+jdim))&
                  + (flxv(k,i,j) + flxv(k,i+1,j) - flxv(k,i,j-jdim)  &
                  - flxv(k,i+1,j-jdim)) * 2.* oneBasic%uc(k,i,j) )
          enddo

          do k = 2,m1-1
             ut(k,i,j) = ut(k,i,j) + c1z * dzt(k) / oneBasic%dn0u(k,i,j) * (  &
                  (flxw(k-1,i,j) + flxw(k-1,i+1,j))  &
                  * (oneBasic%uc(k,i,j) + oneBasic%uc(k-1,i,j))  &
                  - (flxw(k,i,j) + flxw(k,i+1,j))  &
                  * (oneBasic%uc(k,i,j) + oneBasic%uc(k+1,i,j))   &
                  + (flxw(k,i,j) + flxw(k,i+1,j) - flxw(k-1,i,j)  &
                  - flxw(k-1,i+1,j)) * 2.* oneBasic%uc(k,i,j) )
          enddo
       enddo
    enddo

    ! Compute advection contribution to V tendency

    do j = ja,jzv
       do i = ia,iz
          c1z = .25 / grid_g(ngrid)%rtgv(i,j)
          c1x = c1z * grid_g(ngrid)%fmapv(i,j) * grid_g(ngrid)%dxv(i,j)
          c1y = c1z * grid_g(ngrid)%fmapv(i,j) * grid_g(ngrid)%dyv(i,j)

          do k = 2,m1-1
             vt(k,i,j) = vt(k,i,j) + c1x / oneBasic%dn0v(k,i,j) * (  &
                  (flxu(k,i-1,j) + flxu(k,i-1,j+jdim))  &
                  * (oneBasic%vc(k,i,j) + oneBasic%vc(k,i-1,j))  &
                  - (flxu(k,i,j) + flxu(k,i,j+jdim))  &
                  * (oneBasic%vc(k,i,j) + oneBasic%vc(k,i+1,j))  &
                  + (flxu(k,i,j) + flxu(k,i,j+jdim) - flxu(k,i-1,j)  &
                  - flxu(k,i-1,j+jdim)) * 2.* oneBasic%vc(k,i,j) )
          enddo

          do k = 2,m1-1
             vt(k,i,j) = vt(k,i,j) + c1y / oneBasic%dn0v(k,i,j) * (  &
                  (flxv(k,i,j) + flxv(k,i,j-jdim))  &
                  * (oneBasic%vc(k,i,j) + oneBasic%vc(k,i,j-jdim))  &
                  - (flxv(k,i,j) + flxv(k,i,j+jdim))  &
                  * (oneBasic%vc(k,i,j) + oneBasic%vc(k,i,j+jdim))  &
                  + (flxv(k,i,j+jdim) - flxv(k,i,j-jdim))  &
                  * 2.* oneBasic%vc(k,i,j) )
          enddo

          do k = 2,m1-1
             vt(k,i,j) = vt(k,i,j) + c1z * dzt(k) / oneBasic%dn0v(k,i,j) * (  &
                  (flxw(k-1,i,j) + flxw(k-1,i,j+jdim))  &
                  * (oneBasic%vc(k,i,j) + oneBasic%vc(k-1,i,j))  &
                  - (flxw(k,i,j) + flxw(k,i,j+jdim))  &
                  * (oneBasic%vc(k,i,j) + oneBasic%vc(k+1,i,j))  &
                  + (flxw(k,i,j) + flxw(k,i,j+jdim) - flxw(k-1,i,j)  &
                  - flxw(k-1,i,j+jdim)) * 2.* oneBasic%vc(k,i,j) )
          enddo
       enddo
    enddo

    ! Compute advection contribution to W tendency
    do j = ja,jz
       do i = ia,iz
          c1z = .5 / grid_g(ngrid)%rtgt(i,j)
          c1x = c1z * grid_g(ngrid)%fmapt(i,j) * grid_g(ngrid)%dxt(i,j)
          c1y = c1z * grid_g(ngrid)%fmapt(i,j) * grid_g(ngrid)%dyt(i,j)

          do k = 2,m1-2
             wt(k,i,j) = wt(k,i,j)  &
                  + c1x / (oneBasic%dn0(k,i,j) + oneBasic%dn0(k+1,i,j)) * (  &
                  (flxu(k,i-1,j) + flxu(k+1,i-1,j))  &
                  * (oneBasic%wc(k,i,j) + oneBasic%wc(k,i-1,j))  &
                  - (flxu(k,i,j) + flxu(k+1,i,j))  &
                  * (oneBasic%wc(k,i,j) + oneBasic%wc(k,i+1,j))  &
                  + (flxu(k,i,j) + flxu(k+1,i,j) - flxu(k,i-1,j)  &
                  - flxu(k+1,i-1,j)) * 2.* oneBasic%wc(k,i,j) )
          enddo

          do k = 2,m1-2
             wt(k,i,j) = wt(k,i,j)  &
                  + c1y / (oneBasic%dn0(k,i,j) + oneBasic%dn0(k+1,i,j)) * (  &
                  (flxv(k,i,j-jdim) + flxv(k+1,i,j-jdim))  &
                  * (oneBasic%wc(k,i,j) + oneBasic%wc(k,i,j-jdim))  &
                  - (flxv(k,i,j) + flxv(k+1,i,j))  &
                  * (oneBasic%wc(k,i,j) + oneBasic%wc(k,i,j+jdim))  &
                  + (flxv(k,i,j) + flxv(k+1,i,j) - flxv(k,i,j-jdim)  &
                  - flxv(k+1,i,j-jdim)) * 2.* oneBasic%wc(k,i,j) )
          enddo

          do k = 2,m1-2
             wt(k,i,j) = wt(k,i,j)  &
                  + c1z * dzm(k) / (oneBasic%dn0(k,i,j) + oneBasic%dn0(k+1,i,j)) * (  &
                  (flxw(k,i,j) + flxw(k-1,i,j))  &
                  * (oneBasic%wc(k,i,j) + oneBasic%wc(k-1,i,j))  &
                  - (flxw(k,i,j) + flxw(k+1,i,j))  &
                  * (oneBasic%wc(k,i,j) + oneBasic%wc(k+1,i,j))   &
                  + (flxw(k+1,i,j) - flxw(k-1,i,j)) * 2.* oneBasic%wc(k,i,j) )
          enddo
       enddo
    enddo
  end subroutine vel_advectc




  !     *********************************************************************

  subroutine fa_preptc(m1,m2,m3,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df  &
       ,vt3dh,vt3di,vt3dj,vt3dk,mynum, oneBasic, oneAveBasic)
    !> @brief:  fa_preptc                                      
    !! @author:  unknow
    !! @date:  18/Nov/2015
    !! @version:  5.2
    !! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
    !! 
    !! @copyright Under CC-GPL License by INPE/CPTEC
    !! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt


    integer,intent(in) :: m1,m2,m3,mynum
    integer :: j,i,k,im,ip,jm,jp
    type(BasicFields), pointer, intent(in) :: oneBasic
    type(BasicFields), pointer, intent(in) :: oneAveBasic

    real :: c1,c2,c3,c4,rtgti

    real, dimension(m1,m2,m3) :: vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df  &
         ,vt3dh,vt3di,vt3dj,vt3dk

    ! VT3DA, VT3DB, and VT3DC are input as the velocity components (averaged
    ! between past and current time levels) times dtlt.

    ! Add contribution to VT3DC from horiz winds crossing sloping sigma surfaces,
    !    and include 1/rtgt factor in VT3DC
    ! Compute half Courant numbers: VT3DD, VT3DE, and VT3DF
    ! Compute weight at scalar point: VT3DH
    ! Compute advective weights for the linear term: VT3DI, VCTR1, and VCTR2

    call DeepCopyToBasicFields(oneBasic, oneAveBasic)
    do j = 1,m3
       jm = max(1,j-1)
       jp = min(m3,j+1)
       do i = 1,m2
          im = max(1,i-1)
          ip = min(m2,i+1)
          rtgti = 1. / grid_g(ngrid)%rtgt(i,j)
          c1 = .5 * grid_g(ngrid)%dxu(i,j)
          c2 = .5 * grid_g(ngrid)%dyv(i,j)
          c3 = grid_g(ngrid)%dxt(i,j) * grid_g(ngrid)%fmapt(i,j) * rtgti
          c4 = grid_g(ngrid)%dyt(i,j) * grid_g(ngrid)%fmapt(i,j) * rtgti

          do k = 1,m1-1
             vt3dc(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
                  + vt3da(k,im,j) + vt3da(k+1,im,j)) * grid_g(ngrid)%f13t(i,j)  &
                  + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
                  + vt3db(k+1,i,jm)) * grid_g(ngrid)%f23t(i,j)) * hw4(k)  &
                  + vt3dc(k,i,j) * rtgti
             vt3dd(k,i,j) = c1 * vt3da(k,i,j)
             vt3de(k,i,j) = c2 * vt3db(k,i,j)
             vt3df(k,i,j) = .5 * vt3dc(k,i,j) * dzm(k)
             vctr3(k) = 1. / oneBasic%dn0(k,i,j)
             vt3dh(k,i,j) = c3 * vctr3(k)
             vt3dj(k,i,j) = c4 * vctr3(k)
             vt3dk(k,i,j) = dzt(k) * vctr3(k)
          enddo

          ! temporary override
          vt3di(1,i,j) = .5
          vt3di(2,i,j) = .5
          vt3di(3,i,j) = .5
          vt3di(4,i,j) = .5
       enddo
    enddo

    do k = 1,m1-1
       vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
       vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
    enddo

    ! Convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
    ! into mass fluxes times dtlt.

    do j = 1,m3
       do i = 1,m2
          c1 = grid_g(ngrid)%fmapui(i,j) * grid_g(ngrid)%rtgu(i,j)
          c2 = grid_g(ngrid)%fmapvi(i,j) * grid_g(ngrid)%rtgv(i,j)
          do k = 1,m1-1
             vt3da(k,i,j) = vt3da(k,i,j) * c1 * oneBasic%dn0u(k,i,j)
             vt3db(k,i,j) = vt3db(k,i,j) * c2 * oneBasic%dn0v(k,i,j)
             vt3dc(k,i,j) = vt3dc(k,i,j) * .5  &
                  * (oneBasic%dn0(k,i,j) + oneBasic%dn0(k+1,i,j))
          enddo
       enddo
    enddo
    call DeepCopyFromBasicFields(oneBasic, oneAveBasic)
    
  end subroutine fa_preptc

  !     *********************************************************************

  subroutine fa_xc(m1,m2,m3,ia,iz,ja,jz  &
       ,scp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)
    !> @brief: Compute scalar flux times dtlt
    !! @date:  18/Nov/2015
    !! @version:  5.2
    !! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
    !! 
    !! @copyright Under CC-GPL License by INPE/CPTEC
    !! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt

    integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k,mynum

    real :: dfact
    real, dimension(m1,m2,m3) :: scp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di

    dfact = .5
    do j = 1,m3
       do i = 1,iz

          ! Compute scalar flux times dtlt [VT3DG]

          do k = 2,m1-1
             vt3dg(k,i,j) = vt3da(k,i,j)  &
                  * (vt3di(1,i,j) * scr1(k,i,j)  &
                  +  vt3di(2,i,j) * scr1(k,i+1,j)  &
                  +  vt3dd(k,i,j) * (scr1(k,i,j) - scr1(k,i+1,j)))
          enddo

          ! Modify fluxes to retain positive-definiteness on scalar quantities.
          !    If a flux will remove 1/2 quantity during a timestep,
          !    reduce to first order flux. This will remain positive-definite
          !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
          !    both fluxes are evacuating the box.

          do k = 2,m1-1
             if (vt3da(k,i,j) .gt. 0.) then
                if (vt3dg(k,i,j) * vt3dh(k,i,j) .gt.  &
                     dfact * scr1(k,i,j)) then
                   vt3dg(k,i,j) = vt3da(k,i,j) * scr1(k,i,j)
                endif
             elseif (vt3da(k,i,j) .lt. 0.) then
                if (-vt3dg(k,i,j) * vt3dh(k,i+1,j) .gt.  &
                     dfact * scr1(k,i+1,j)) then
                   vt3dg(k,i,j) = vt3da(k,i,j) * scr1(k,i+1,j)
                endif
             endif
          enddo
       enddo
    enddo

    ! Compute flux divergence

    do j = 1,m3
       do i = ia,iz
          do k = 2,m1-1
             scr1(k,i,j) = scr1(k,i,j)  &
                  + vt3dh(k,i,j) * (vt3dg(k,i-1,j) - vt3dg(k,i,j)  &
                  + scp(k,i,j) * (vt3da(k,i,j) - vt3da(k,i-1,j)))
          enddo
       enddo
    enddo
  end subroutine fa_xc

  !     *********************************************************************

  subroutine fa_yc(m1,m2,m3,ia,iz,ja,jz  &
       ,scp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)
    !> @brief: Compute scalar flux 
    !! @date:  18/Nov/2015
    !! @version:  5.2
    !! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
    !! 
    !! @copyright Under CC-GPL License by INPE/CPTEC
    !! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt

    integer :: m1,m2,m3,ia,iz,ja,jz,jdim,mynum,i,j,k

    real :: dfact
    real, dimension(m1,m2,m3) :: scp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di

    dfact = .5
    do j = 1,jz
       do i = ia,iz

          ! Compute scalar flux VT3DG

          do k = 2,m1-1
             vt3dg(k,i,j) = vt3db(k,i,j)  &
                  * (vt3di(3,i,j) * scr1(k,i,j)  &
                  +  vt3di(4,i,j) * scr1(k,i,j+jdim)  &
                  +  vt3de(k,i,j) * (scr1(k,i,j) - scr1(k,i,j+jdim)))
          enddo

          !      Modify fluxes to retain positive-definiteness on scalar quantities.
          !         If a flux will remove 1/2 quantity during a timestep,
          !         reduce to first order flux. This will remain positive-definite
          !         under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
          !         both fluxes are evacuating the box.

          do k = 2,m1-1
             if (vt3db(k,i,j) .gt. 0.) then
                if (vt3dg(k,i,j) * vt3dj(k,i,j) .gt.  &
                     dfact * scr1(k,i,j)) then
                   vt3dg(k,i,j) = vt3db(k,i,j) * scr1(k,i,j)
                endif
             elseif (vt3db(k,i,j) .lt. 0.) then
                if (-vt3dg(k,i,j) * vt3dj(k,i,j+jdim) .gt.  &
                     dfact * scr1(k,i,j+jdim)) then
                   vt3dg(k,i,j) = vt3db(k,i,j) * scr1(k,i,j+jdim)
                endif
             endif
          enddo
       enddo
    enddo

    ! Compute flux divergence

    do j = ja,jz
       do i = ia,iz
          do k = 2,m1-1
             scr1(k,i,j) = scr1(k,i,j)  &
                  + vt3dj(k,i,j) * (vt3dg(k,i,j-jdim) - vt3dg(k,i,j)  &
                  + scp(k,i,j) * (vt3db(k,i,j) - vt3db(k,i,j-jdim)))
          enddo
       enddo
    enddo
  end subroutine fa_yc

  !     *********************************************************************

  subroutine fa_zc(m1,m2,m3,ia,iz,ja,jz  &
       ,scp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)
    !> @brief: Compute scalar flux 
    !! @date:  18/Nov/2015
    !! @version:  5.2
    !! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
    !! 
    !! @copyright Under CC-GPL License by INPE/CPTEC
    !! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt

    integer :: m1,m2,m3,ia,iz,ja,jz,mynum,i,j,k

    real :: dfact
    real, dimension(m1,m2,m3) :: scp,scr1,vt3dc,vt3df,vt3dg,vt3dk
    real, dimension(*) :: vctr1,vctr2

    dfact = .5
    do j = ja,jz
       do i = ia,iz

          ! Compute scalar flux VT3DG

          do k = 1,m1-1
             vt3dg(k,i,j) = vt3dc(k,i,j)  &
                  * (vctr1(k) * scr1(k,i,j)  &
                  +  vctr2(k) * scr1(k+1,i,j)  &
                  +  vt3df(k,i,j) * (scr1(k,i,j) - scr1(k+1,i,j)))
          enddo

          ! Modify fluxes to retain positive-definiteness on scalar quantities.
          !    If a flux will remove 1/2 quantity during a timestep,
          !    reduce to first order flux. This will remain positive-definite
          !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
          !    both fluxes are evacuating the box.

          do k = 1,m1-1
             if (vt3dc(k,i,j) .gt. 0.) then
                if (vt3dg(k,i,j) * vt3dk(k,i,j) .gt.  &
                     dfact * scr1(k,i,j)) then
                   vt3dg(k,i,j) = vt3dc(k,i,j) * scr1(k,i,j)
                endif
             elseif (vt3dc(k,i,j) .lt. 0.) then
                if (-vt3dg(k,i,j) * vt3dk(k+1,i,j) .gt.  &
                     dfact * scr1(k+1,i,j)) then
                   vt3dg(k,i,j) = vt3dc(k,i,j) * scr1(k+1,i,j)
                endif
             endif
          enddo
       enddo
    enddo

    ! Compute flux divergence

    do j = ja,jz
       do i = ia,iz
          do k = 2,m1-1
             scr1(k,i,j) = scr1(k,i,j)  &
                  + vt3dk(k,i,j) * (vt3dg(k-1,i,j) - vt3dg(k,i,j)  &
                  + scp(k,i,j) * (vt3dc(k,i,j) - vt3dc(k-1,i,j)))
          enddo
       enddo
    enddo
  end subroutine fa_zc

  !     ****************************************************************

  subroutine advtndc(m1,m2,m3,ia,iz,ja,jz,scp,sca,sct,dtl,mynum)
    !> @brief: advtndc
    !! @date:  18/Nov/2015
    !! @version:  5.2
    !! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
    !! 
    !! @copyright Under CC-GPL License by INPE/CPTEC
    !! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt

    integer :: m1,m2,m3,ia,iz,ja,jz,mynum,i,j,k

    real :: dtl,dtli
    real, dimension(m1,m2,m3) :: scp,sca,sct

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(advtndc)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1-1
       write(str(2),"(i8)") ia
       write(str(3),"(i8)") iz
       write(str(4),"(i8)") ja
       write(str(5),"(i8)") jz
       call MsgDump(h//" updates (2:"//trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
    end if

    dtli = 1. / dtl
    do j = ja,jz
       do i = ia,iz
          do k = 2,m1-1
             sct(k,i,j) = sct(k,i,j) + (sca(k,i,j)-scp(k,i,j)) * dtli
          enddo
       enddo
    enddo
  end subroutine advtndc

end module ModRadvc
