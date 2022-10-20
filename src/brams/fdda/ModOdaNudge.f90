!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModOdaNudge

  use ModOdaProcObs, only: &
       oda_proc_obs
  
  use ModOdaKrig, only: &
       krig
  
  use mem_oda, only: &
       oda_g, &
       ckrg, &
       ekobs, &
       frqoda, &
       rkobs, &
       tkobs, &
       tnudoda, &
       todabeg, &
       todaend, &
       ukobs, &
       vkobs, &
       wt_oda_rt, &
       wt_oda_th, &
       wt_oda_uv, &
       xkobs, &
       ykobs, &
       zkobs, &
       wt_oda_grid, &
       nstkrg, &
       roda_hgt, &
       rmaxkrg, &
       akrg, &
       roda_upa0, &
       roda_sfc0, &
       roda_sfce, &
       roda_upae, &
       roda_zfact, &
       caxkrg, &
       caykrg, &
       cazkrg
  
  use mem_scratch, only: &
       scratch
  
  use mem_tend, only: &
       tend
       
  use ModBasicFields, only: &
       BasicFields
  
  use io_params, only: &
       timstr
  
  use node_mod, only: &
       mynum, &
       nodemxp, &
       nodemyp, &
       nodei0, &
       nodej0, &
       nodeia, &
       nodeiz, &
       nodeja, &
       nodejz
  
  use mem_grid, only: &
       grid_g, &
       ngrids, &
       ngrid, &
       time, &
       ztn, &
       nnzp, &
       nnxp, &
       nnyp, &
       dtlong, &
       dtlongn, &
       xtn, &
       ytn, &
       ztn

  implicit none

  private

  public :: oda_nudge



contains

  subroutine oda_nudge(oneBasicFields)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    
    integer :: ncall=0, nobs, ng,i,j,k

    real, allocatable, dimension(:,:) :: plt

    ! Namelist variables:

    ! - Flag to turn on oda                     (IF_ODA)
    ! - Model start and end time for oda        (TODABEG,TODAEND)
    ! - Frequency of tendency updates           (FRQODA) (0. for every coarse grid dt)
    ! - For each grid:
    !     - Nudging timescale                   (TNUDODA)
    !     - Number of vertical levels to apply  (no)
    !     - Flag to turn on for this grid       (IG_ODA)


    ! Covariance stuff:

    ! For each grid:
    ! - Surface radii - 1) e-2  2) zero         (RODA_SFCE, RODA_SFC0) (SFC0 <= UPA0)
    ! - "mid-level" radii - same                (RODA_UPAE, RODA_UPA0)
    ! - "mid-level" radii height                (RODA_HGT)
    !         (linear change from surface to mid, constant above)
    !
    ! - vertical "100" factor                   (RODA_ZFACT)

    ! Obs processing:
    ! - Time interpolate limit (TIL)- if the future-past obs time 
    !    is > this limit, do not use to interpolate (ODA_SFC_TIL,ODA_UPA_TIL)
    !
    ! - Time extrapolate limit (TEL)- if past/future obs is greater than TIL,
    !    but less than TEL, use obs                 (ODA_SFC_TEL,ODA_UPA_TEL)



    if (ncall==0) then
       call oda_nudge_init(ngrids, nnzp(1), nodemxp(mynum,1), nodemyp(mynum,1))
       ncall = 1
    endif

    !print*,'time,frqoda:',mynum,time,frqoda,mod(time,frqoda),dtlongn(1)

    if (frqoda>0. .and. mod(time,frqoda)>=dtlongn(1) ) goto 20


    if ( (ngrid==1 .and. time>=todabeg .and. time<=todaend)  &
         .or. time==timstr ) then

       ! Compute new krigged fields and variances

       do ng=1,ngrids

          if (wt_oda_grid(ng)>0.0) then

             ! Process observations and call krigging 

             call oda_proc_obs(nnzp(ng), nodemxp(mynum,ng), nodemyp(mynum,ng), &
                  nodei0(mynum,ng), nodej0(mynum,ng),  &
                  oneBasicFields%pp, oneBasicFields%pi0,  &
                  scratch%scr1, ng, nobs)

             oda_g(ng)%uk(1:nnzp(ng),1:nodemxp(mynum,ng),1:nodemyp(mynum,ng)) = 0.
             oda_g(ng)%ukv(1:nnzp(ng),1:nodemxp(mynum,ng),1:nodemyp(mynum,ng)) = 0.
             call krig(nnzp(ng),nodemxp(mynum,ng),nodemyp(mynum,ng)  &
                  ,xtn(1+nodei0(mynum,ng):,ng),ytn(1+nodej0(mynum,ng):,ng),ztn(:,ng)  &
                  ,nobs,xkobs,ykobs,zkobs,ekobs,ukobs  &
                  ,ng,nnzp(ng),grid_g(ng)%topt  &
                  ,oda_g(ng)%uk,oda_g(ng)%ukv,1.)

             oda_g(ng)%vk(1:nnzp(ng),1:nodemxp(mynum,ng),1:nodemyp(mynum,ng))=0.
             oda_g(ng)%vkv(1:nnzp(ng),1:nodemxp(mynum,ng),1:nodemyp(mynum,ng))=0.
             call krig(nnzp(ng),nodemxp(mynum,ng),nodemyp(mynum,ng)  &
                  ,xtn(1+nodei0(mynum,ng):,ng),ytn(1+nodej0(mynum,ng):,ng),ztn(:,ng)  &
                  ,nobs,xkobs,ykobs,zkobs,ekobs,vkobs  &
                  ,ng,nnzp(ng),grid_g(ng)%topt  &
                  ,oda_g(ng)%vk,oda_g(ng)%vkv,1.)

             oda_g(ng)%tk(1:nnzp(ng),1:nodemxp(mynum,ng),1:nodemyp(mynum,ng))=0.
             oda_g(ng)%tkv(1:nnzp(ng),1:nodemxp(mynum,ng),1:nodemyp(mynum,ng))=0.
             call krig(nnzp(ng),nodemxp(mynum,ng),nodemyp(mynum,ng)  &
                  ,xtn(1+nodei0(mynum,ng):,ng),ytn(1+nodej0(mynum,ng):,ng),ztn(:,ng)  &
                  ,nobs,xkobs,ykobs,zkobs,ekobs,tkobs  &
                  ,ng,nnzp(ng),grid_g(ng)%topt  &
                  ,oda_g(ng)%tk,oda_g(ng)%tkv,1.)


             oda_g(ng)%rk(1:nnzp(ng),1:nodemxp(mynum,ng),1:nodemyp(mynum,ng))=0.
             oda_g(ng)%rkv(1:nnzp(ng),1:nodemxp(mynum,ng),1:nodemyp(mynum,ng))=0.
             call krig(nnzp(ng),nodemxp(mynum,ng),nodemyp(mynum,ng)  &
                  ,xtn(1+nodei0(mynum,ng):,ng),ytn(1+nodej0(mynum,ng):,ng),ztn(:,ng)  &
                  ,nobs,xkobs,ykobs,zkobs,ekobs,rkobs  &
                  ,ng,nnzp(ng),grid_g(ng)%topt  &
                  ,oda_g(ng)%rk,oda_g(ng)%rkv,1.)
          endif

       enddo

    endif

20  continue

    ! Compute and apply tendencies

    ng=ngrid
    if (wt_oda_grid(ng)>0.0 .and. time>=todabeg .and. time<=todaend) then
       if (allocated(plt)) deallocate(plt); allocate(plt(nnxp(ng),nnyp(ng)))

       call oda_tendency(nnzp(ng), nodemxp(mynum,ng), nodemyp(mynum,ng), &
            nodeia(mynum,ng), nodeiz(mynum,ng), &
            nodeja(mynum,ng), nodejz(mynum,ng), &
            oneBasicFields%up, tend%ut,  &
            oda_g(ng)%uk, oda_g(ng)%ukv, &
            wt_oda_uv*wt_oda_grid(ng)/tnudoda, &
            nodei0(mynum,ng), nodej0(mynum,ng))
       call oda_tendency(nnzp(ng), nodemxp(mynum,ng), nodemyp(mynum,ng), &
            nodeia(mynum,ng), nodeiz(mynum,ng), &
            nodeja(mynum,ng), nodejz(mynum,ng), &
            oneBasicFields%vp, tend%vt, &
            oda_g(ng)%vk, oda_g(ng)%vkv, &
            wt_oda_uv*wt_oda_grid(ng)/tnudoda, &
            nodei0(mynum,ng), nodej0(mynum,ng))
       call oda_tendency(nnzp(ng), nodemxp(mynum,ng), nodemyp(mynum,ng), &
            nodeia(mynum,ng), nodeiz(mynum,ng), &
            nodeja(mynum,ng), nodejz(mynum,ng), &
            oneBasicFields%theta, tend%tht, &
            oda_g(ng)%tk,oda_g(ng)%tkv,  &
            wt_oda_th*wt_oda_grid(ng)/tnudoda, &
            nodei0(mynum,ng), nodej0(mynum,ng))
       call oda_tendency(nnzp(ng), nodemxp(mynum,ng), nodemyp(mynum,ng), &
            nodeia(mynum,ng), nodeiz(mynum,ng), &
            nodeja(mynum,ng), nodejz(mynum,ng), &
            oneBasicFields%rtp, tend%rtt, &
            oda_g(ng)%rk, oda_g(ng)%rkv, &
            wt_oda_rt*wt_oda_grid(ng)/tnudoda, &
            nodei0(mynum,ng), nodej0(mynum,ng))


    endif
  end subroutine oda_nudge

  !========================================================================

  subroutine oda_nudge_init(ngg,m1m,m2m,m3m)
    integer :: ngg
    integer, dimension(ngg) :: m1m,m2m,m3m

    integer ::ng,k

    ! Turn namelist parameters into krigging routine parameters
    print*,'nud_init:'

    do ng=1,ngrids
       if (wt_oda_grid(ng) > 0.0) then
          nstkrg(ng)=1
          ckrg(1:nstkrg(ng),ng) = 1
          do k=1,nnzp(ng)
             if(ztn(k,ng) < roda_hgt(ng)) then
                rmaxkrg(k,ng)=   &
                     roda_sfc0(ng)+(roda_upa0(ng)-roda_sfc0(ng))  &
                     * ztn(k,ng)/roda_hgt(ng)
                ! Kriging routine needs the following to be negative...
                akrg(k,ng)   = -  &
                     (roda_sfce(ng)+(roda_upae(ng)-roda_sfce(ng))  &
                     * ztn(k,ng)/roda_hgt(ng)) 
             else
                rmaxkrg(k,ng)=roda_upa0(ng)
                akrg(k,ng) = -roda_upae(ng)
             endif
             print*,k,ztn(k,ng), roda_hgt(ng),rmaxkrg(k,ng),akrg(k,ng)
          enddo

          caxkrg(1,ng)= 1. ; caykrg(1,ng)= 0. ; cazkrg(1,ng)= 0.
          caxkrg(2,ng)= 0. ; caykrg(2,ng)= 1. ; cazkrg(2,ng)= 0.
          caxkrg(3,ng)= 0. ; caykrg(3,ng)= 0. ; cazkrg(3,ng)= roda_zfact(ng)

          ! Initialize analysis and variance fields so zero tendency will be computed

          oda_g(ng)%uk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
          oda_g(ng)%vk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
          oda_g(ng)%tk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
          oda_g(ng)%rk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
          oda_g(ng)%ukv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
          oda_g(ng)%vkv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
          oda_g(ng)%tkv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
          oda_g(ng)%rkv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.

       endif
    enddo
  end subroutine oda_nudge_init

  !========================================================================

  subroutine oda_tendency(n1,n2,n3,ia,iz,ja,jz,ap,at,ak,akv,tscalei,i0,j0)
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: i0
    integer, intent(in) :: j0
    real, intent(in) :: tscalei
    real, intent(in) :: ap(n1,n2,n3)
    real, intent(inout) :: at(n1,n2,n3)
    real, intent(in) :: ak(n1,n2,n3)
    real, intent(in) :: akv(n1,n2,n3)

    real :: ttmin,ttmax

    integer :: i,j,k,iiix,jjjx,kkkx,iiim,jjjm,kkkm
    real :: fnna

    ! Compute and apply tendencies based on krigged field/variance 

    !print*,'oda_tend:',ia,iz,ja,jz,n1,n2,n3

    ttmax=-1e20; ttmin=1e20
    do j=ja,jz
       do i=ia,iz
          do k=1,n1
             fnna=(1.0-akv(k,i,j)/2.0)*(ak(k,i,j)-ap(k,i,j))*tscalei
             at(k,i,j)=at(k,i,j)+fnna
             if(fnna > ttmax) then
                ttmax=fnna
                iiix=i+i0
                jjjx=j+j0
                kkkx=k
             endif
             if(fnna < ttmin) then
                ttmin=fnna
                iiim=i+i0
                jjjm=j+j0
                kkkm=k
             endif
          enddo
       enddo
    enddo
  end subroutine oda_tendency
end module ModOdaNudge
