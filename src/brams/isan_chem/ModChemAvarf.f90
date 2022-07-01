!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModChemAvarf

  use ModAVarF, only: &
       isnsig, &
       vshyd, &
       visurf, &
       varuv

  use isan_coms, only: &
       guess1st,  &
       hybbot,    &
       hybtop,    &
       nisn,      &
       nsigz,     &
       sigzwt,    &
       pi_s,      &
       pi_p,      &
       levth, &
       is_grids,  &
       rr_scr1, &
       rs_sfp,    &
       rs_sft,    &
       rs_snow,   &
       rs_sst,    &
       rs_slp

  use chem_isan_coms, only: &
       chem_is_grids,  &
       aer_is_grids,   &
       nspecies,       &
       pi_sc,  &
       ps_sc, &
       pi_aer_sc,   &
       ps_aer_sc, &
       nspecies_aer_in

  use mem_chem1, only: &
       chem_assim

  use mem_aer1, only: &
       aer_assim
  
  use mem_grid, only: &
       grid_g,   &
       nnxp,     &
       nnyp,     &
       nnzp,     &
       ztop,     &
       ztn

  use ModRbnd, only: &
       topset, &
       botset

  use rconstants, only: &
       cp,  &
       g,   &
       p00i,&
       rocp

  implicit none

  private

  public :: chem_makevarf
  public :: chem_varfile_nstfeed

contains



  
  subroutine chem_makevarf(ng)
    integer, intent(in) :: ng

    !---------------------------------------------------------------+
    !    Interpolate model grid from isentropic/sigmaz/surface data
    !---------------------------------------------------------------+

    !            Vertically interpolate isentropic data to sigma-z levels

    call isnsig(nnzp(ng),nnxp(ng),nnyp(ng) ,is_grids(ng)%rr_u  &
         ,is_grids(ng)%rr_v      ,is_grids(ng)%rr_t  &
         ,is_grids(ng)%rr_r      ,is_grids(ng)%rr_p  &
         ,grid_g(ng)%topt,ztn(1,ng),ztop)

    if(CHEM_ASSIM == 1 .and. nspecies>0 ) then
       call chem_isnsig(nnzp(ng),nnxp(ng),nnyp(ng),nspecies     &
            ,chem_is_grids(ng)%rr_sc        &
            ,grid_g(ng)%topt,ztn(:,ng),ztop)
    endif
    if(AER_ASSIM == 1 .and. nspecies_aer_in>0 ) then
       call aer_isnsig(nnzp(ng),nnxp(ng),nnyp(ng),nspecies_aer_in     &
            ,aer_is_grids(ng)%rr_sc        &
            ,grid_g(ng)%topt,ztn(:,ng),ztop)
    endif



    !            Compute Exner function on model sigma-z surfaces
    !              and change relative humidity to mixing ratio.

    call vshyd(nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_p  &
         ,is_grids(ng)%rr_t       ,is_grids(ng)%rr_r  &
         ,grid_g(ng)%topt,grid_g(ng)%rtgt,ztn(1,ng))

    !          Combine surface analysis with the upper air data.

    call visurf(nnzp(ng),nnxp(ng),nnyp(ng) ,is_grids(ng)%rr_u  &
         ,is_grids(ng)%rr_v        ,is_grids(ng)%rr_t  &
         ,is_grids(ng)%rr_r        ,is_grids(ng)%rr_p  &
         ,grid_g(ng)%topt,grid_g(ng)%rtgt,ztn(1,ng))

    is_grids(ng)%rr_slp (1:nnxp(ng),1:nnyp(ng)) =rs_slp (1:nnxp(ng),1:nnyp(ng))
    is_grids(ng)%rr_sfp (1:nnxp(ng),1:nnyp(ng)) =rs_sfp (1:nnxp(ng),1:nnyp(ng))
    is_grids(ng)%rr_sft (1:nnxp(ng),1:nnyp(ng)) =rs_sft (1:nnxp(ng),1:nnyp(ng))
    is_grids(ng)%rr_snow(1:nnxp(ng),1:nnyp(ng)) =rs_snow(1:nnxp(ng),1:nnyp(ng))
    is_grids(ng)%rr_sst (1:nnxp(ng),1:nnyp(ng)) =rs_sst (1:nnxp(ng),1:nnyp(ng))

    !          average the velocities to the correct points in the stagger
    !             and rotate for polar stereographic transformation.

    call varuv(nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_u  &
         ,is_grids(ng)%rr_v)

    return
  end subroutine chem_makevarf








  !     ****************************************************************

  subroutine chem_varfile_nstfeed(ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c &
       ,nbot,ntop)
    integer, intent(in) :: ifm
    integer, intent(in) :: icm
    integer, intent(in) :: n1f
    integer, intent(in) :: n2f
    integer, intent(in) :: n3f
    integer, intent(in) :: n1c
    integer, intent(in) :: n2c
    integer, intent(in) :: n3c
    integer, intent(in) :: nbot
    integer, intent(in) :: ntop

    integer :: nspc

    !     Feed back the finer mesh to the coarser mesh.

    call fdback(is_grids(icm)%rr_u   (1,1,1),is_grids(ifm)%rr_u   (1,1,1) &
         ,is_grids(icm)%rr_dn0u(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
         ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'u',rr_scr1(1))

    call fdback(is_grids(icm)%rr_v   (1,1,1),is_grids(ifm)%rr_v   (1,1,1) &
         ,is_grids(icm)%rr_dn0v(1,1,1),is_grids(ifm)%rr_dn0v(1,1,1) &
         ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'v',rr_scr1(1))

    call fdback(is_grids(icm)%rr_p  (1,1,1),is_grids(ifm)%rr_p  (1,1,1) &
         ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
         ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'p',rr_scr1(1))

    call fdback(is_grids(icm)%rr_t  (1,1,1),is_grids(ifm)%rr_t  (1,1,1) &
         ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
         ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'t',rr_scr1(1))

    call fdback(is_grids(icm)%rr_r  (1,1,1),is_grids(ifm)%rr_r  (1,1,1) &
         ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
         ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'t',rr_scr1(1))

    if(CHEM_ASSIM == 1 .and. nspecies>0 ) then
       do nspc=1,nspecies
          print*,'fdback for spc=',nspc
          call fdback(chem_is_grids(icm)%rr_sc  (1,1,1,nspc) &
               ,chem_is_grids(ifm)%rr_sc  (1,1,1,nspc) &
               ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
               ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'t',rr_scr1(1))
       enddo
    endif
    if(aer_ASSIM == 1 .and. nspecies_aer_in>0 ) then
       do nspc=1,nspecies_aer_in
          print*,'fdback for spc=',nspc
          call fdback(aer_is_grids(icm)%rr_sc  (1,1,1,nspc) &
               ,aer_is_grids(ifm)%rr_sc  (1,1,1,nspc) &
               ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
               ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'t',rr_scr1(1))
       enddo
    endif

    if(nbot == 1) then
       call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15 ,is_grids(icm)%rr_u,'U')
       call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15 ,is_grids(icm)%rr_v,'V')
       call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15 ,is_grids(icm)%rr_p,'P')
       call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15 ,is_grids(icm)%rr_t,'T')
       call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15 ,is_grids(icm)%rr_r,'T')
    endif

    if(ntop == 1) then
       call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15,is_grids(icm)%rr_u,'U')
       call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15,is_grids(icm)%rr_v,'V')
       call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15,is_grids(icm)%rr_p,'P')
       call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15,is_grids(icm)%rr_t,'T')
       call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
            ,15,is_grids(icm)%rr_r,'T')
    endif


    return
  end subroutine chem_varfile_nstfeed





  subroutine chem_isnsig(n1,n2,n3,nsp_dummy,sc,topt,zt,ztop)
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: nsp_dummy
    real, intent(inout) :: sc(:,:,:,:) ! (n1,n2,n3,nsp_dummy) 
    real, intent(in) :: topt(:,:) ! (n2,n3) 
    real, intent(in) :: zt(:) ! (n1)
    real, intent(in) :: ztop

    real, allocatable :: v1(:),v2(:),v6(:,:)
    integer :: i,j,k,nki,ki,kbeg,kend,nspc
    real :: rtg,wtsz

    allocate (v1(n1),v2(nisn),v6(nisn,nsp_dummy))

    do j=1,n3
       do i=1,n2

          rtg=(1.-topt(i,j)/ztop)
          do k=1,n1
             v1(k)=topt(i,j)+zt(k)*rtg
          enddo

          if (guess1st.eq.'PRESS') then
             do nspc=1,nsp_dummy
                nki=0
                do ki=1,nisn
                   if(pi_p(i,j,ki).lt.1e20 .and. pi_s(i,j,ki).lt.1e20  &
                        .and. pi_sc(i,j,ki,nspc).lt.1e20)  then
                      nki=nki+1
                      v2(nki)=(pi_s(i,j,ki)-cp*levth(ki)  &
                           *(pi_p(i,j,ki)*p00i)**rocp)/g
                      v6(nki,nspc)=pi_sc(i,j,ki,nspc)
                   endif
                enddo  ! loop ki
             enddo   ! loop nspc

             do nspc=1,nsp_dummy
                call htint(nki,v6(1,nspc),v2,n1,sc(1,i,j,nspc),v1)
             enddo  ! loop nspc

          endif


          kbeg=n1+1
          kend=n1+1
          do k=1,nsigz
             if(zt(k).gt.hybbot) then
                kbeg=k-1
                exit
             endif
          enddo

          do k=1,nsigz
             if(zt(k).gt.hybtop) then
                kend=k
                exit
             endif
          enddo

          do k=1,nsigz

             if(k.lt.kbeg) then
                wtsz=sigzwt
             elseif(k.gt.kend) then
                wtsz=0.
             elseif(k.ge.kbeg.and.k.le.kend) then
                wtsz= (zt(kend)-zt(k))  &
                     /(zt(kend)-zt(kbeg))  &
                     *sigzwt
             endif

             do nspc=1,nsp_dummy
                sc(k,i,j,nspc)=(1.-wtsz)*sc(k,i,j,nspc)+wtsz*ps_sc(i,j,k,nspc)
             enddo  ! loop nspc

          enddo

       enddo
    enddo

    deallocate (v1,v2,v6)

    return
  end subroutine chem_isnsig




  

  subroutine aer_isnsig(n1,n2,n3,nsp_dummy,sc,topt,zt,ztop)
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: nsp_dummy
    real, intent(inout) :: sc(:,:,:,:) ! (n1,n2,n3,nsp_dummy)
    real, intent(in) :: topt(:,:) ! (n2,n3)
    real, intent(in) :: zt(:) ! (n1)
    real, intent(in) :: ztop

    real, allocatable :: v1(:),v2(:),v6(:,:)
    integer :: i,j,k,nki,ki,kbeg,kend,nspc
    real :: rtg,wtsz

    allocate (v1(n1),v2(nisn),v6(nisn,nsp_dummy))

    do j=1,n3
       do i=1,n2

          rtg=(1.-topt(i,j)/ztop)
          do k=1,n1
             v1(k)=topt(i,j)+zt(k)*rtg
          enddo

          if (guess1st.eq.'PRESS') then
             do nspc=1,nsp_dummy
                nki=0
                do ki=1,nisn
                   if(pi_p(i,j,ki).lt.1e20 .and. pi_s(i,j,ki).lt.1e20  &
                        .and. pi_aer_sc(i,j,ki,nspc).lt.1e20)  then
                      nki=nki+1
                      v2(nki)=(pi_s(i,j,ki)-cp*levth(ki)  &
                           *(pi_p(i,j,ki)*p00i)**rocp)/g
                      v6(nki,nspc)=pi_aer_sc(i,j,ki,nspc)
                   endif
                enddo  ! loop ki
             enddo   ! loop nspc

             do nspc=1,nsp_dummy
                call htint(nki,v6(1,nspc),v2,n1,sc(1,i,j,nspc),v1)
             enddo  ! loop nspc

          endif


          kbeg=n1+1
          kend=n1+1
          do k=1,nsigz
             if(zt(k).gt.hybbot) then
                kbeg=k-1
                exit
             endif
          enddo

          do k=1,nsigz
             if(zt(k).gt.hybtop) then
                kend=k
                exit
             endif
          enddo

          do k=1,nsigz

             if(k.lt.kbeg) then
                wtsz=sigzwt
             elseif(k.gt.kend) then
                wtsz=0.
             elseif(k.ge.kbeg.and.k.le.kend) then
                wtsz= (zt(kend)-zt(k))  &
                     /(zt(kend)-zt(kbeg))  &
                     *sigzwt
             endif

             do nspc=1,nsp_dummy
                sc(k,i,j,nspc)=(1.-wtsz)*sc(k,i,j,nspc)+wtsz*ps_aer_sc(i,j,k,nspc)
             enddo  ! loop nspc

          enddo

       enddo
    enddo

    deallocate (v1,v2,v6)

    return
  end subroutine aer_isnsig
end module ModChemAvarf
