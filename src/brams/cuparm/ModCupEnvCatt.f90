module ModCupEnvCatt

  implicit none

  private

  public :: cup_kbcon_catt
  public :: get_zi

contains
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  subroutine cup_kbcon_catt(cap_inc,iloop, k22, kbcon, he_cup, hes_cup, &
       hkb, kzi, mix, mgmxp, mkx, mgmzp, istart, iend, ierr, kbmax, p_cup,&
       cap_max,j)

    !--- Determine the level of convective cloud base  - kbcon
    ! kbcon  = LFC of parcel from k22
    ! k22    = updraft originating level
    ! hkb    = moist static energy at originating level

    ! arguments
    real, intent(in) :: cap_inc(:) ! (mgmxp)
    integer, intent(in) :: iloop
    integer, intent(inout) :: k22(:) ! (mgmxp)
    integer, intent(inout) :: kbcon(:) ! (mgmxp)
    real, intent(in) :: he_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: hes_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: hkb(:) ! (mgmxp) **(JP)** unused
    integer, intent(in) :: kzi(:) ! (mgmxp) **(JP)** unused
    integer, intent(in) :: mix !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mkx !**(JP)** unused
    integer, intent(in) :: mgmzp
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(inout) :: ierr(:) ! (mgmxp)
    integer, intent(in) :: kbmax(:) ! (mgmxp)
    real, intent(in) :: p_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: cap_max(:) ! (mgmxp)
    integer, intent(in) :: j !**(JP)** unused

    ! local vars
    integer :: k
    integer :: kp
    integer :: i
    real :: pbcdif,plus,hkbpbl

    do i=istart,iend
       kbcon(i)=1
       if (ierr(i).ne.0) cycle
       kbcon(i)=k22(i)

       go to 32
31     continue
       kbcon(i)=kbcon(i)+1

       if (kbcon(i).gt.kbmax(i)+2) then
          if(iloop.lt.4)ierr(i)=3
          if(iloop.eq.4)ierr(i)=997

          cycle
       endif
32     continue

       if (he_cup(i,k22(i)).lt.hes_cup(i,kbcon(i))) then

          go to 31
       else
       endif

       !     cloud base pressure and max moist static energy pressure
       !     i.e., the depth (in mb) of the layer of negative buoyancy

       if (kbcon(i)-k22(i).eq.1) then
          cycle
       endif

       pbcdif = -p_cup(i,kbcon(i)) + p_cup(i,k22(i))
       plus=cap_max(i)-float(iloop-1)*cap_inc(i)
       plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))

       if (pbcdif.gt.plus) then
          k22(i)   = k22(i)+1
          kbcon(i) = k22(i)         
          go to 32
       endif
    enddo
  end subroutine cup_kbcon_catt


  !-------------------------------------------------------------------

  subroutine get_zi(mix,mgmxp,mkx,mgmzp,istart,iend,j,ierr,kzi,tkeg, &
       rcpg,z,ztop,tkmin)

    ! arguments
    integer, intent(in) :: mix !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mkx
    integer, intent(in) :: mgmzp
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(in) :: j !**(JP)** unused
    integer, intent(in) :: ierr(:) ! (mgmxp)
    integer, intent(inout) :: kzi(:) ! (mgmxp)
    real, intent(in) :: tkeg(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: rcpg(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: z(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: ztop(:) ! (mgmxp)
    real, intent(in) ::  tkmin

    ! local vars
    integer :: i
    integer :: k
    integer :: kzimax
    integer :: ktke_max
    real ::  tke_tmp
    real, parameter ::  rcpmin=1.e-5
    real, parameter ::  pblhmax=3000.

    do i=istart,iend
       kzi(i)  = 2

       if(ierr(i).eq.0)then
          tke_tmp = 0.
          ktke_max= 1
          !---  max level for kzi
          do k=1,mkx
             if(z(i,k).ge. pblhmax+ztop(i)) then
                kzimax = min(k,mgmzp-1)
                exit
             endif
          enddo
          !level of max tke  below kzimax and w/out clouds
          do  k=1,kzimax
             if(rcpg(i,k) .lt. rcpmin) then
                if( tkeg(i,k) .ge. tke_tmp) then
                   tke_tmp = tkeg(i,k)
                   cycle
                else
                   ktke_max= max(1,k-1)
                   exit
                endif
             endif
          enddo

          do k=ktke_max,kzimax+1
             if(rcpg(i,k) .lt. rcpmin) then
                if(tkeg(i,k) .gt. 1.1*tkmin)  then
                   kzi(i) = k
                   cycle
                endif
             else
                kzi(i) = k
                exit
             endif
          enddo
          kzi(i) = max(2     ,kzi(i))
          kzi(i) = min(kzimax,kzi(i))
       endif
    enddo
  end subroutine get_zi
end module ModCupEnvCatt
