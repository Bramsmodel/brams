!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModTurbKAdap

  implicit none

  private
  public :: strain_adap

contains


  subroutine strain_adap(m1,m2,m3,ia,iz,ja,jz  &
       ,ia_1,ja_1,iz1,jz1,jd,lpu_r,lpv_r,lpw_r  &
       ,up,vp,wp,vt3da,vt3db,vt3dc,vt3dd,vt3de  &
       ,vt3df,vt3dg,vt3dh,vt3di,vt3dn,scr2,idiffk  &
       ,dxm,dxt,dxu,dxv,dym,dyt,dyu,dyv,dzm,dzt)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ia_1
    integer, intent(in) :: ja_1
    integer, intent(in) :: iz1
    integer, intent(in) :: jz1
    integer, intent(in) :: jd
    real, intent(in) :: lpu_r(m2,m3)
    real, intent(in) :: lpv_r(m2,m3)
    real, intent(in) :: lpw_r(m2,m3)
    real, intent(in) :: up(m1,m2,m3)
    real, intent(in) :: vp(m1,m2,m3)
    real, intent(in) :: wp(m1,m2,m3)
    real, intent(out) :: vt3da(m1,m2,m3)
    real, intent(out) :: vt3db(m1,m2,m3)
    real, intent(out) :: vt3dc(m1,m2,m3)
    real, intent(out) :: vt3dd(m1,m2,m3)
    real, intent(out) :: vt3de(m1,m2,m3)
    real, intent(out) :: vt3df(m1,m2,m3)
    real, intent(out) :: vt3dg(m1,m2,m3)
    real, intent(out) :: vt3dh(m1,m2,m3)
    real, intent(out) :: vt3di(m1,m2,m3)
    real, intent(out) :: vt3dn(m1,m2,m3)
    real, intent(out) :: scr2(m1,m2,m3)
    integer, intent(in) :: idiffk
    real, intent(in) :: dxm(m2,m3)
    real, intent(in) :: dxt(m2,m3)
    real, intent(in) :: dxu(m2,m3)
    real, intent(in) :: dxv(m2,m3)
    real, intent(in) :: dym(m2,m3)
    real, intent(in) :: dyt(m2,m3)
    real, intent(in) :: dyu(m2,m3)
    real, intent(in) :: dyv(m2,m3)
    real, intent(in) :: dzm(:)
    real, intent(in) :: dzt(:)

    integer :: i
    integer :: j
    integer :: k
    integer :: ka
    integer :: lpu(m2,m3)
    integer :: lpv(m2,m3)
    integer :: lpw(m2,m3)
    

    !lfr: Solving a  problem with integer inside vtables
    lpu=int(lpu_r);lpv=int(lpv_r);lpw=int(lpw_r)

    ! du/dx

    do j = ja,jz
       do i = ia,iz1
          ka = max(lpu(i,j),lpu(i-1,j))
          do k = 1,ka-1
             vt3da(k,i,j) = 0.
          enddo
          do k = ka,m1
             vt3da(k,i,j) = (up(k,i,j) - up(k,i-1,j)) * dxt(i,j)
          enddo
       enddo
    enddo

    ! dv/dx

    do j = ja_1,jz
       do i = ia_1,iz
          ka = max(lpv(i+1,j),lpv(i,j))
          do k = 1,ka-1
             vt3db(k,i,j) = 0.
          enddo
          do k = ka,m1
             vt3db(k,i,j) = (vp(k,i+1,j) - vp(k,i,j)) * dxm(i,j)
          enddo
       enddo
    enddo

    ! dw/dx

    do j = ja,jz
       do i = ia_1,iz
          ka = max(lpw(i+1,j),lpw(i,j))
          do k = 1,ka-1
             vt3df(k,i,j) = 0.
          enddo
          do k = ka,m1
             vt3df(k,i,j) = (wp(k,i+1,j) - wp(k,i,j)) * dxu(i,j)
          enddo
       enddo
    enddo

    ! du/dy

    do j = ja_1,jz
       do i = ia_1,iz
          ka = max(lpu(i,j+1),lpu(i,j))
          do k = 1,ka-1
             vt3dn(k,i,j) = 0.
          enddo
          do k = ka,m1
             vt3dn(k,i,j) = (up(k,i,j+1) - up(k,i,j)) * dym(i,j)
          enddo
       enddo
    enddo

    ! dv/dy

    do j = ja,jz1
       do i = ia,iz
          ka = max(lpv(i,j),lpv(i,j-1))
          do k = 1,ka-1
             vt3dc(k,i,j) = 0.
          enddo
          do k = ka,m1
             vt3dc(k,i,j) = (vp(k,i,j) - vp(k,i,j-1)) * dyt(i,j)
          enddo
       enddo
    enddo

    ! dw/dy

    do j = ja_1,jz
       do i = ia,iz
          ka = max(lpw(i,j+1),lpw(i,j))
          do k = 1,ka-1
             vt3dg(k,i,j) = 0.
          enddo
          do k = ka,m1
             vt3dg(k,i,j) = (wp(k,i,j+1) - wp(k,i,j)) * dyv(i,j)
          enddo
       enddo
    enddo

    ! du/dz

    do j = ja,jz
       do i = ia_1,iz
          ka = lpu(i,j)
          do k = 1,ka-1
             vt3dd(k,i,j) = 0.
          enddo
          do k = ka,m1-1
             vt3dd(k,i,j) = (up(k+1,i,j) - up(k,i,j)) * dzm(k)
          enddo
       enddo
    enddo

    ! dv/dz

    do j = ja_1,jz
       do i = ia,iz
          ka = lpv(i,j)
          do k = 1,ka-1
             vt3de(k,i,j) = 0.
          enddo
          do k = ka,m1-1
             vt3de(k,i,j) = (vp(k+1,i,j) - vp(k,i,j)) * dzm(k)
          enddo
       enddo
    enddo

    ! dw/dz

    if (idiffk .ge. 3) then
       do j = ja,jz
          do i = ia,iz
             ka = lpw(i,j)
             do k = 1,ka
                scr2(k,i,j) = 0.
             enddo
             do k = ka+1,m1
                scr2(k,i,j) = (wp(k,i,j) - wp(k-1,i,j)) * dzt(k)
             enddo
          enddo
       enddo
    endif

    if (idiffk .le. 2) then
       do j = ja,jz
          do i = ia,iz
             do k = 2,m1-1
                vt3dh(k,i,j) =2. * (vt3da(k,i,j) * vt3da(k,i,j)  &
                     + vt3dc(k,i,j) * vt3dc(k,i,j))                &
                     + .0625 * (vt3db(k,i,j) + vt3db(k,i-1,j)      &
                     + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)         &
                     + vt3dn(k,i,j) + vt3dn(k,i-1,j)               &
                     + vt3dn(k,i,j-jd) + vt3dn(k,i-1,j-jd)) ** 2
                vt3di(k,i,j) = .0625 * ((vt3dd(k,i,j) + vt3dd(k-1,i,j)  &
                     + vt3dd(k,i-1,j) + vt3dd(k-1,i-1,j)) ** 2            &
                     + (vt3de(k,i,j) + vt3de(k-1,i,j)                     &
                     + vt3de(k,i,j-jd) + vt3de(k-1,i,j-jd)) ** 2)
             enddo
          enddo
       enddo
    else
       do j = ja,jz
          do i = ia,iz
             do k = 2,m1-1
                vt3da(k,i,j) = 2. * vt3da(k,i,j)
                vt3dc(k,i,j) = 2. * vt3dc(k,i,j)
                scr2(k,i,j) = 2. * scr2(k,i,j)
                vt3db(k,i,j) = vt3db(k,i,j) + vt3dn(k,i,j)
                vt3dn(k,i,j) = vt3db(k,i,j)
                vt3dd(k,i,j) = vt3dd(k,i,j) + vt3df(k,i,j)
                vt3de(k,i,j) = vt3de(k,i,j) + vt3dg(k,i,j)
                vt3di(k,i,j) = 0.333333  &
                     * (vt3da(k,i,j) + vt3dc(k,i,j) + scr2(k,i,j))
             enddo
          enddo

          do k = 2,m1-1
             vt3da(k,iz1,j) = 2. * vt3da(k,iz1,j)
             vt3db(k,ia_1,j) = vt3db(k,ia_1,j) + vt3dn(k,ia_1,j)
             vt3dn(k,ia_1,j) = vt3db(k,ia_1,j)
             vt3dd(k,ia_1,j) = vt3dd(k,ia_1,j) + vt3df(k,ia_1,j)
          enddo
       enddo

       do i = ia_1,iz
          do k = 2,m1-1
             vt3dc(k,i,jz1) = 2. * vt3dc(k,i,jz1)
             vt3db(k,i,ja_1) = vt3db(k,i,ja_1) + vt3dn(k,i,ja_1)
             vt3dn(k,i,ja_1) = vt3db(k,i,ja_1)
             vt3de(k,i,ja_1) = vt3de(k,i,ja_1) + vt3dg(k,i,ja_1)
          enddo
       enddo

       do j = ja,jz
          do i = ia,iz
             do k = 2,m1-1
                vt3dh(k,i,j) = .5 * (  &
                     (vt3da(k,i,j) - vt3di(k,i,j)) ** 2           &
                     + (vt3dc(k,i,j) - vt3di(k,i,j)) ** 2         &
                     + ( scr2(k,i,j) - vt3di(k,i,j)) ** 2)        &
                     + .0625 * ((vt3db(k,i,j) + vt3db(k,i-1,j)    &
                     + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)) ** 2  &
                     + (vt3dd(k,i,j) + vt3dd(k,i-1,j)             &
                     + vt3dd(k-1,i,j) + vt3dd(k-1,i-1,j)) ** 2    &
                     + (vt3de(k,i,j) + vt3de(k-1,i,j)             &
                     + vt3de(k,i,j-jd) + vt3de(k-1,i,j-jd)) ** 2)
                vt3di(k,i,j) = vt3dh(k,i,j)
             enddo
          enddo
       enddo
    endif

    return
  end subroutine strain_adap
end module ModTurbKAdap
