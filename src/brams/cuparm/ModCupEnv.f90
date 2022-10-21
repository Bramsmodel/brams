module ModCupEnv

  implicit none
  
  private

  public :: cup_env
  public :: cup_env_clev
  public :: cup_ktop
  public :: maximi
  public :: minimi

contains



  !--------------------------------------------------------------------
  subroutine cup_env(j, ipr, jpr, z, qes, he, hes, t, q, p, z1,        &
       mix, mgmxp, mkx, mgmzp, istart, iend, psur, ierr, tcrit, itest)

    ! arguments:
    integer, intent(in) :: j, ipr, jpr, mix, mgmxp, mkx, mgmzp, istart, &
         iend, itest
    integer, intent(in) :: ierr(:) ! (mgmxp)
    real, intent(in)    :: tcrit
    real, intent(in)    :: t(:,:) ! (mgmxp,mgmzp)
    real, intent(in)    :: p(:,:) ! (mgmxp,mgmzp)
    real, intent(in)    :: z1(:) ! (mgmxp)
    real, intent(in)    :: psur(:) ! (mgmxp)
    real, intent(inout) :: qes(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: q(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: z(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: he(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: hes(:,:) ! (mgmxp,mgmzp)

    ! local variables:
    integer :: i, k, iph,  m
    real    :: tv(mgmxp,mgmzp)
    real    :: ae(2), be(2), ht(2), e, tvbar
    real    :: xl, cp

    xl=2.5e06
    cp=1004.
    ht(1)=xl/cp
    ht(2)=2.834e6/cp
    be(1)=.622*ht(1)/.286
    ae(1)=be(1)/273.+alog(610.71)
    be(2)=.622*ht(2)/.286
    ae(2)=be(2)/273.+alog(610.71)
    do k=1,mkx
       do i=istart,iend
          if(ierr(i).eq.0)then
             !sgb - iph is for phase, dependent on tcrit (water or ice)
             iph = 1
             if (t(i,k).le.tcrit) iph = 2
             e = exp(ae(iph)-be(iph)/t(i,k))
             qes(i,k) = .622*e/(100.*p(i,k)-e)
             if (qes(i,k) .le. 1.e-08)   qes(i,k)=1.e-08
             if (q(i,k)   .gt. qes(i,k)) q(i,k)=qes(i,k)
             tv(i,k) = t(i,k)+.608*q(i,k)*t(i,k)
          endif
       enddo
    enddo

    !--- z's are calculated with changed h's and q's and t's
    !--- if itest=2

    if (itest.ne.2) then
       do i=istart,iend
          if (ierr(i).eq.0) then
             z(i,1)=max(0.,z1(i))-(alog(p(i,1))-alog(psur(i)))*287.*tv(i,1)/9.81
          endif
       enddo

       !--- calculate heights
       do k=2,mkx
          do i=istart,iend
             if (ierr(i).eq.0) then
                tvbar =  .5*tv(i,k)+.5*tv(i,k-1)
                z(i,k) = z(i,k-1)-(alog(p(i,k))-alog(p(i,k-1)))*287.*tvbar/9.81
             endif
          enddo
       enddo
    else
       do k=1,mkx
          do i=istart,iend
             if (ierr(i).eq.0) then
                z(i,k) = (he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
                z(i,k) = max(1.e-3,z(i,k))
             endif
          enddo
       enddo
    endif

    !--- calculate moist static energy - he
    !--- saturated moist static energy - hes

    do k=1,mkx
       do i=istart,iend
          if (ierr(i).eq.0) then
             if (itest.eq.0) he(i,k) = 9.81*z(i,k)+1004.*t(i,k)+2.5e06*q(i,k)
             hes(i,k) = 9.81*z(i,k)+1004.*t(i,k)+2.5e06*qes(i,k)

             if (he(i,k).ge.hes(i,k)) he(i,k) = hes(i,k)
          endif
       enddo
    enddo

    return
  end subroutine cup_env

  !--------------------------------------------------------------------
  subroutine cup_env_clev(j, ipr, jpr, t, qes, q, he, hes, z, p, qes_cup,  &
       q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur,       &
       mix, mgmxp, mkx, mgmzp, istart, iend, ierr, z1)

    ! arguments
    integer, intent(in) :: j
    integer, intent(in) :: ipr !**(JP)** unused
    integer, intent(in) :: jpr !**(JP)** unused
    real, intent(in) :: t(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: qes(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: q(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: he(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: hes(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: z(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: p(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: qes_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: q_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: he_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: hes_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: z_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: p_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: gamma_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(inout) :: t_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: psur(:) ! (mgmxp)
    integer, intent(in) :: mix !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mkx
    integer, intent(in) :: mgmzp
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(in):: ierr(:) ! (mgmxp)
    real, intent(in) :: z1(:) ! (mgmxp)

    ! local vars
    integer :: i
    integer :: k
    integer :: m
    real :: xl, rv, cp

    xl=2.5e6
    rv=461.9
    cp=1004.
    do k=2,mkx
       do i=istart,iend
          if (ierr(i).eq.0)then
             qes_cup(i,k) = .5*(qes(i,k-1) + qes(i,k))
             q_cup(i,k)   = .5*(  q(i,k-1) +   q(i,k))
             hes_cup(i,k) = .5*(hes(i,k-1) + hes(i,k))
             he_cup(i,k)  = .5*( he(i,k-1) +  he(i,k))
             if (he_cup(i,k).gt.hes_cup(i,k)) he_cup(i,k) = hes_cup(i,k)

             z_cup(i,k) = .5*(z(i,k-1) + z(i,k))
             p_cup(i,k) = .5*(p(i,k-1) + p(i,k))
             t_cup(i,k) = .5*(t(i,k-1) + t(i,k))

             gamma_cup(i,k) =(xl/cp)*(xl/(rv*t_cup(i,k)*t_cup(i,k)))*qes_cup(i,k)
          endif
       enddo
    enddo
    do i=istart,iend
       if (ierr(i).eq.0)then
          qes_cup(i,1) = qes(i,1)
          q_cup(i,1)   = q(i,1)
          hes_cup(i,1) = hes(i,1)
          he_cup(i,1)  = he(i,1)

          !srf
          !        z_cup(i,1) = .5*( z(i,1) +   z1(i))
          !        p_cup(i,1) = .5*( p(i,1) + psur(i))
          !        t_cup(i,1) =      t(i,1)
          z_cup(i,1)   = z1(i)
          p_cup(i,1)   = psur(i)
          t_cup(i,1)   = t(i,1)
          !srf
          gamma_cup(i,1) = xl/cp*(xl/(rv*t_cup(i,1)*t_cup(i,1)))*qes_cup(i,1)
       endif
    enddo
  end subroutine cup_env_clev


  !--------------------------------------------------------------------
  subroutine cup_direction2(i, j, dir, id, mix, mjx, mgmxp, mgmyp,  &
       massflx, iresult, num, imass, nall, maxens3, massfld)

    ! arguments:
    integer, intent(in)  :: i, j, mix, mjx, mgmxp, mgmyp, num, imass, nall, maxens3
    integer, intent(out) :: iresult
    !srf      integer id(mix,mjx)
    !srf      real dir(mgmxp), massflx(mix,mjx)
    integer, intent(in)  ::  id(:,:) ! (mgmxp,mgmyp)
    real, intent(inout)  ::  dir(:) ! (mgmxp)
    real, intent(inout)  ::  massflx(:,:) ! (mgmxp,mgmyp)
    real, intent(out)    ::  massfld

    ! local variables:
    integer :: k
    integer :: ia, ib, ja, jb
    real    ::  diff

    if (imass.eq.1) then
       massfld = massflx(i,j)
    endif
    iresult=0
    !      return

    diff = 22.5
    if (dir(i).lt.22.5) dir(i)=360.+dir(i)
    if (id(i,j).eq.1)   iresult=1
    !      ja=max(2,j-1)
    !      ia=max(2,i-1)
    !      jb=min(mjx-1,j+1)
    !      ib=min(mix-1,i+1)
    ja=j-1
    ia=i-1
    jb=j+1
    ib=i+1
    if (dir(i).gt.90.-diff.and.dir(i).le.90.+diff) then

       !--- steering flow from the east
       if (id(ib,j).eq.1) then
          iresult=1
          if (imass.eq.1) then
             massfld = max(massflx(ib,j),massflx(i,j))
          endif
          return
       endif
    else if (dir(i).gt.135.-diff.and.dir(i).le.135.+diff)then

       !--- steering flow from the south-east
       if (id(ib,ja).eq.1) then
          iresult=1
          if (imass.eq.1) then
             massfld = max(massflx(ib,ja),massflx(i,j))
          endif
          return
       endif

       !--- steering flow from the south
    else if (dir(i).gt.180.-diff.and.dir(i).le.180.+diff) then
       if (id(i,ja).eq.1) then
          iresult=1
          if (imass.eq.1) then
             massfld = max(massflx(i,ja),massflx(i,j))
          endif
          return
       endif

       !--- steering flow from the south west
    else if (dir(i).gt.225.-diff.and.dir(i).le.225.+diff) then
       if (id(ia,ja).eq.1) then
          iresult=1
          if (imass.eq.1) then
             massfld = max(massflx(ia,ja),massflx(i,j))
          endif
          return
       endif

       !--- steering flow from the west
    else if (dir(i).gt.270.-diff.and.dir(i).le.270.+diff) then
       if (id(ia,j).eq.1) then
          iresult=1
          if (imass.eq.1)then
             massfld = max(massflx(ia,j),massflx(i,j))
          endif
          return
       endif

       !--- steering flow from the north-west
    else if (dir(i).gt.305.-diff.and.dir(i).le.305.+diff) then
       if (id(ia,jb).eq.1) then
          iresult=1
          if (imass.eq.1) then
             massfld = max(massflx(ia,jb),massflx(i,j))
          endif
          return
       endif

       !--- steering flow from the north
    else if (dir(i).gt.360.-diff.and.dir(i).le.360.+diff) then
       if (id(i,jb).eq.1) then
          iresult=1
          if (imass.eq.1) then
             massfld = max(massflx(i,jb),massflx(i,j))
          endif
          return
       endif

       !--- steering flow from the north-east
    else if (dir(i).gt.45.-diff.and.dir(i).le.45.+diff) then
       if (id(ib,jb).eq.1) then
          iresult=1
          if (imass.eq.1) then
             massfld = max(massflx(ib,jb),massflx(i,j))
          endif
          return
       endif
    endif

    !srf
    !      if(massfld.gt.0.) print*,'---- mass fld=',massfld        
    !srf

    return
  end subroutine cup_direction2

  !--------------------------------------------------------------------
  subroutine cup_ktop(ilo, dby, kbcon, ktop, mix, mgmxp, mkx, mgmzp,   &
       istart, iend, ierr)

    ! arguments
    integer, intent(in) :: ilo
    real, intent(inout) :: dby(:,:) ! (mgmxp,mgmzp)
    integer, intent(in) :: kbcon(:) ! (mgmxp)
    integer, intent(inout) :: ktop(:) ! (mgmxp)
    integer, intent(in) :: mix !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mkx
    integer, intent(in) :: mgmzp
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(inout) :: ierr(:) ! (mgmxp)
    ! local variables
    integer :: i, k

    do i=istart,iend
       !do 42 i=istart,iend
       ktop(i)=1
       if (ierr(i).eq.0) then

          do k=kbcon(i)+1,mkx-2
             if (dby(i,k).le.0.) then
                ktop(i) = k-1
                go to 41
             endif
          enddo
          if (ilo.eq.1) ierr(i)=5
          if (ilo.eq.2) ierr(i)=998
          go to 42
41        continue
          do k=ktop(i)+1,mkx
             dby(i,k)=0.
          enddo
       endif
42     continue
    enddo
  end subroutine cup_ktop


  !--------------------------------------------------------------------
  subroutine cup_kbcon(iloop, k22, kbcon, he_cup, hes_cup, mix, mgmxp,  &
       mkx, mgmzp, istart, iend, ierr, kbmax, p_cup, cap_max)

    !--- determine the level of convective cloud base  - kbcon

    ! arguments
    integer, intent(in) :: iloop
    integer, intent(inout) :: k22(:) ! (mgmxp)
    integer, intent(inout) :: kbcon(:) ! (mgmxp)
    real, intent(in) :: he_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: hes_cup(:,:) ! (mgmxp,mgmzp)
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
    ! local variables
    integer :: i
    real :: pbcdif

    do i=istart,iend
       !do 27 i=istart,iend
       kbcon(i)=1
       !if (ierr(i).ne.0) go to 27
       if (ierr(i).ne.0) cycle
       kbcon(i)=k22(i)

       !srf
       !       print*,'-----------------------------------------'
       !       print*,'1',i,ierr(i),k22(i),kbcon(i),cap_max(i)
       !       do k=mkx,k22(i),-1
       !       do k=mkx,1,-1
       !       print*,k,he_cup(i,k),hes_cup(i,k),p_cup(i,k)
       !       enddo
       !srf

       go to 32
31     continue
       kbcon(i)=kbcon(i)+1
       if (kbcon(i).gt.kbmax(i)+2) then
          if(iloop.eq.1)ierr(i)=3
          if(iloop.eq.2)ierr(i)=997
          !srf
          !       print*,'2',i,ierr(i),k22(i),kbcon(i),cap_max(i)
          !srf

          !go to 27
          cycle
       endif
32     continue
       if (he_cup(i,k22(i)).lt.hes_cup(i,kbcon(i))) go to 31

       !     cloud base pressure and max moist static energy pressure
       !     i.e., the depth (in mb) of the layer of negative buoyancy

       if (kbcon(i)-k22(i).eq.1) then
          !srf
          !       print*,'2',i,ierr(i),k22(i),kbcon(i),cap_max(i)
          !srf
          !go to 27
          cycle
       endif

       pbcdif = -p_cup(i,kbcon(i)) + p_cup(i,k22(i))
       if (pbcdif.gt.cap_max(i)) then
          k22(i)   = k22(i)+1
          kbcon(i) = k22(i)
          go to 32
       endif
       !27   continue
    enddo

    return
  end subroutine cup_kbcon



  !--------------------------------------------------------------------
  subroutine cup_kbcon_cin(iloop, k22, kbcon, he_cup, hes_cup, z, tmean, qes,  &
       mix, mgmxp, mkx, mgmzp, istart, iend, ierr, kbmax, p_cup, cap_max)

    ! arguments
    integer, intent(in) :: iloop
    integer, intent(inout) :: k22(:) ! (mgmxp)
    integer, intent(inout) :: kbcon(:) ! (mgmxp)
    real, intent(in) :: he_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: hes_cup(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: z(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: tmean(:,:) ! (mgmxp,mgmzp)
    real, intent(in) :: qes(:,:) ! (mgmxp,mgmzp)
    integer, intent(in) :: mix !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mkx !**(JP)** unused
    integer, intent(in) :: mgmzp
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(inout) :: ierr(:) ! (mgmxp)
    integer, intent(in) :: kbmax(:) ! (mgmxp)
    real, intent(in) :: p_cup(:,:) ! (mgmxp,mgmzp) **(JP)** unused
    real, intent(in) :: cap_max(:) ! (mgmxp)
    ! local variables
    integer :: i
    real :: pbcdif
    real :: cin
    real :: cin_max
    real :: dh
    real :: tprim
    real :: gamma
    real :: lovcp_p, xl, rv, cp
    
    xl=2.5e6
    rv=461.525
    cp=1004.
    lovcp_p = xl/cp
    !srf- -----------------------------------
    !
    !--- determine the level of convective cloud base  - kbcon
    !
    do i=istart,iend
       !do 27 i=istart,iend
       cin_max=-cap_max(i)

       kbcon(i)=1
       cin = 0.
       !if (ierr(i).ne.0) go to 27
       if (ierr(i).ne.0) cycle
       kbcon(i)=k22(i)
       go to 32
31     continue
       kbcon(i)=kbcon(i)+1
       if (kbcon(i).gt.kbmax(i)+2) then
          if (iloop.eq.1) ierr(i)=3
          if (iloop.eq.2) ierr(i)=997

          !srf
          !      print*,'kcon_cin ',i,k22(i),kbcon(i),ierr(i),cin
          !srf

          !go to 27
          cycle
       endif
32     continue
       dh = he_cup(i,k22(i)) - hes_cup(i,kbcon(i))
       if (dh.lt. 0.) then
          gamma = lovcp_p*(xl/(rv*(tmean(i,k22(i))**2)))*qes(i,k22(i))
          tprim = dh/(cp*(1.+gamma))

          cin   = cin + 9.8066*tprim* &
               (z(i,k22(i))-z(i,k22(i)-1))/tmean(i,k22(i))

          go to 31
       end if

       !     if negative energy in negatively buoyant layer
       !       exceeds convective inhibition (cin) threshold,
       !       then set k22 level one level up and see if that
       !       will work.

       if (cin.lt.cin_max) then
          k22(i)=k22(i)+1
          kbcon(i)=k22(i)
          go to 32
       endif

       !srf
       !      print*,'kcon_cin ',i,k22(i),kbcon(i),ierr(i),cin
       !srf

       !27   continue
    enddo

    return
  end subroutine cup_kbcon_cin



  !--------------------------------------------------------------------
  subroutine minimi(array, mxx, mgmxp, mzx, mgmzp, ks, kend,  &
       kt, istart, iend, ierr)

    ! arguments
    real, intent(in) :: array(:,:) ! (mgmxp,mgmzp)
    integer, intent(in) :: mxx !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mzx !**(JP)** unused
    integer, intent(in) :: mgmzp
    integer, intent(in) :: ks(:) ! (mgmxp)
    integer, intent(in) :: kend(:) ! (mgmxp)
    integer, intent(inout) :: kt(:) ! (mgmxp)
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(in) :: ierr(:) ! (mgmxp)
    ! local variables
    integer :: kstop
    integer :: i, k
    real :: x(mgmxp)
    
    do i=istart,iend
       !do 200 i=istart,iend
       kt(i)=ks(i)
       if (ierr(i).eq.0) then
          x(i)=array(i,ks(i))
          kstop=max(ks(i)+1,kend(i))

          do k=ks(i)+1,kstop
             !do 100 k=ks(i)+1,kstop
             if(array(i,k).lt.x(i)) then
                x(i)=array(i,k)
                kt(i)=k
             endif
             !100        continue
          enddo
       endif
       !200  continue
    enddo

    return
  end subroutine minimi

  !--------------------------------------------------------------------
  subroutine maximi(array, mxx, mgmxp, mzx, mgmzp, ks, ke, maxx, istart, &
       iend, ierr)

    ! arguments
    real, intent(in) :: array(:,:) ! (mgmxp,mgmzp)
    integer, intent(in) :: mxx !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mzx !**(JP)** unused
    integer, intent(in) :: mgmzp
    integer, intent(in) :: ks
    integer, intent(in) :: ke(:) ! (mgmxp)
    integer, intent(inout) :: maxx(:) ! (mgmxp)
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(in) :: ierr(:) ! (mgmxp)
    ! local variables
    integer :: i, k
    real :: x(mgmxp)
    real :: xar

    do i=istart,iend
       !do 200 i=istart,iend
       maxx(i)=ks
       if(ierr(i).eq.0)then
          x(i)=array(i,ks)

          do k=ks,ke(i)
             !do 100 k=ks,ke(i)
             xar=array(i,k)
             if(xar.ge.x(i)) then
                x(i)=xar
                maxx(i)=k
             endif
             !100        continue
          enddo
       endif
       !200  continue
    enddo

    return
  end subroutine maximi
end module ModCupEnv
