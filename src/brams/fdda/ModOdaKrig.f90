!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModOdaKrig

  use mem_oda, only: &
       nstkrg,   &
       akrg,     &
       caxkrg,   &
       caykrg,   &
       cazkrg,   &
       ckrg, &
       rmaxkrg

  implicit none

  private

  public :: krig
contains

  subroutine krig(n1,n2,n3,x,y,z,ndata,xd,yd,zd,ed,dvar,ngrid  &
       ,nnanzp,topt,var,varkrg,cazmod)
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3

    real, intent(in) :: x(:)
    real, intent(in) :: y(:)
    real, intent(in) :: z(:)

    integer, intent(in) :: ndata
    real, intent(in) :: xd(:)
    real, intent(in) :: yd(:)
    real, intent(in) :: zd(:)
    real, intent(in) :: ed(:)
    real, intent(in) :: dvar(:)

    integer, intent(in) :: ngrid
    integer, intent(in) :: nnanzp
    real, intent(in) :: topt(:,:)
    real, intent(out) :: var(:,:,:)
    real, intent(out) :: varkrg(:,:,:)
    real, intent(in) :: cazmod

    integer, parameter :: maxobs=7000,maxobs2=7000
    real, dimension(maxobs) :: xkp,ykp,zkp,dkp,c20,szg 
    real(kind=8), dimension(maxobs) :: alph,c2 
    real(kind=8), dimension(maxobs2) :: c1
    character(len=80) :: label
    integer :: p0,p1,p2,p3

    real :: sx,sy,sz,szgm,r,dmean
    integer :: i,j,k,nt,ntt,nxp,nyp,ii,n,m,nkp,ni,jk,ik,ktilt,l,ip,jp
    !
    !     Kriging optimal interpolation subroutine
    !
    !     Values at grid points are based on known data using Kriging.
    !     II, IND(II), NT, and NTT are indices for storing the Kriging
    !     configuration and weights to be used in Kriging the simulated
    !     data.


    dmean=0.0
    ii=1

    nt=0
    ntt=0
    nxp=n2
    nyp=n3
    do j=1,nyp
       do i=1,nxp
          do k=1,nnanzp
             ! ip=11; jp=8
             ! if(i==ip.and.j==jp)print*,'------------------'
             ! print*,'------------------',i,j,k
             var(k,i,j)=0.0
             varkrg(k,i,j)=0.0
             ii=ii+1
             !
             !     Calculate separation vector between grid point I,J,K and data
             !     points, retain those with magnitude less than RMAXKRG
             !
             n=0
             do m=1,ndata
                if( dvar(m) < -998. ) cycle
                sx=x(i)-xd(m)
                sy=y(j)-yd(m)
                sz=z(k)-zd(m)
                szg(m)=abs(topt(i,j)-ed(m))
                szg(m)=0.
                szgm=szg(m)
                r=(sx**2+sy**2+(sz*cazkrg(3,ngrid)*cazmod)**2)**0.5
                szgm=0.0
                if(r.gt.rmaxkrg(k,ngrid) ) cycle
                n=n+1
                nt=nt+1
                !
                !     Calculate covariance between data point and grid point
                !
                c2(n)=1.0-gam(sx,sy,sz,ngrid,k,szgm,cazmod)
                xkp(n)=xd(m)
                ykp(n)=yd(m)
                zkp(n)=zd(m)
                dkp(n)=dvar(m)
             enddo

             nkp=n+1
             if(nkp.eq.1) then
                var(k,i,j)=dmean
                varkrg(k,i,j)=2.0
                go to 901
             endif
             !
             !     Non-bias condition
             !
             c2(n+1)=1.0
             do ni=1,nkp
                c20(ni)=c2(ni)
             enddo
             !
             !     Calculate data-data covariance matrix (N+1 * N+1) for data
             !     configuration about grid point I,J,K
             !
             p0=0
             do jk=1,n
                do ik=1,jk
                   sx=xkp(jk)-xkp(ik)
                   sy=ykp(jk)-ykp(ik)
                   sz=zkp(jk)-zkp(ik)
                   p0=p0+1
                   szgm=abs(szg(jk)-szg(ik))
                   szgm=0.0
                   c1(p0)=1.0-gam(sx,sy,sz,ngrid,k,szgm,cazmod)
                   !if(i==ip.and.j==jp)print*,'c1def:',jk,ik,k,p0,c1(p0)
                enddo
             enddo
             p1=1+(n*(n+1)/2)
             p2=(n+1)*(n+2)/2-1
             !if(i==ip.and.j==jp)print*,'p12:',p1,p2
             !
             !     Last column in matrix is all 1'S; non-bias condition
             !
             do p0=p1,p2
                c1(p0)=1.0
             enddo
             p3=p2+1
             c1(p3)=0.0
             !
             !      Call matrix solver
             !
             call relms(alph,c1,c2,nkp,1,ktilt,i,j)
             if(ktilt /= 0) then
                var(k,i,j)=dmean
                varkrg(k,i,j)=2.0
                cycle
             endif

             !
             !      Calculate Kriged grid point
             !
             do l=1,n
                var(k,i,j)=var(k,i,j)+alph(l)*dkp(l)
                !!     if(i==ip.and.j==jp)print*,'vars:',l,var(k,i,j),alph(l),dkp(l),ktilt
                varkrg(k,i,j)=varkrg(k,i,j)+alph(l)*c20(l)
                ntt=ntt+1
             enddo
             varkrg(k,i,j)=1.0-alph(nkp)-varkrg(k,i,j)
             !    if(i==ip.and.j==jp)print*,'vars:',var(k,i,j),varkrg(k,i,j)
901          continue
          enddo
       enddo
    enddo

    return
  end subroutine krig
  !
  !     ****************************************************************
  !
  real function gam(hx,hy,hz,ngrid,k,hzg,cazmod)
    real, intent(in) :: hx
    real, intent(in) :: hy
    real, intent(in) :: hz
    integer, intent(in) :: ngrid
    integer, intent(in) :: k
    real, intent(in) :: hzg
    real, intent(in) :: cazmod


    real :: h,dxmod,dymod,dx,dy,dz
    integer :: is,ijs,ijs1,ijs2

    !
    !     Variogram calculation
    !
    gam=0.0
    h=(hx*hx+hy*hy+hz*hz)**0.5
    if(h .lt. 1.e-03) return
    do is=1,nstkrg(ngrid)
       ijs=is
       ijs1=ijs+nstkrg(ngrid)
       ijs2=ijs1+nstkrg(ngrid)
       dxmod=hzg/750.0*(-akrg(1,ngrid))
       dymod=dxmod
       dx=hx*caxkrg(ijs,ngrid)+hy*caxkrg(ijs1,ngrid)+  &
            hz*caxkrg(ijs2,ngrid)
       dx=abs(dx)+dxmod
       dy=hx*caykrg(ijs,ngrid)+hy*caykrg(ijs1,ngrid)+  &
            hz*caykrg(ijs2,ngrid)
       dy=abs(dy)+dymod
       dz=hx*cazkrg(ijs,ngrid)+hy*cazkrg(ijs1,ngrid)+  &
            hz*cazkrg(ijs2,ngrid)*cazmod
       h=(dx*dx+dy*dy+dz*dz)**0.5
       if(akrg(k,ngrid) < 0.0) then
          gam=gam+ckrg(is,ngrid)*(1.-exp(h/akrg(k,ngrid)))
       else if(akrg(k,ngrid) > 0.0) then
          if(h .ge. akrg(k,ngrid)) then
             gam=gam+ckrg(is,ngrid)
          else
             gam=gam+ckrg(is,ngrid)*(1.5*h/akrg(k,ngrid)-  &
                  0.5*h*h*h/(akrg(k,ngrid)*akrg(k,ngrid)*akrg(k,ngrid)))
          end if
       else
          gam=gam+ckrg(is,ngrid)
       end if
    end do
  end function gam



  subroutine relms(x,a,b,m,n,ktilt,ip,jp)
    real(kind=8), intent(inout) :: x(:)
    real(kind=8), intent(inout) :: a(:)
    real(kind=8), intent(inout) :: b(:)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(out) :: ktilt
    integer, intent(in) :: ip
    integer, intent(in) :: jp
    !
    !     Matrix solver for Kriging equations
    !
    real(kind=8) :: r,piv,tol
    real :: ak
    integer :: i,nm,m1,kk,k,lp,km1,ll,ij,ii,j,llb,in

    if (m .le. 0) then
       ktilt=-1
       return
    end if

    if (m .le. 1) then
       if(a(1) .eq. 0.) then
          ktilt=1
          return
       else
          do i=1,n
             x(i)=b(i)/a(1)
          end do
          ktilt=0
          return
       end if
    end if
    !
    !     Initialize
    !
    tol=0.0
    ktilt=0
    nm=n*m
    m1=m-1
    kk=0
    !
    !     Start triangulation
    !
    do k=1,m1
       kk=kk+k
       ak=a(kk)
       !if(ip==4.and.jp==11) print*,'relms1:',k,kk,ak,tol
       if(ak .eq. tol) then
          ktilt=k
          return
       end if
       piv=1./ak
       ii=kk
       lp=0
       km1=k-1
       do i=k,m1
          ll=ii
          ii=ii+i
          r=a(ii)*piv
          lp=lp+1
          ij=ii-km1
          do  j=i,m1
             ij=ij+j
             ll=ll+j
             a(ij)=a(ij)-r*a(ll)
          end do
          do llb=k,nm,m
             in=llb+lp
             b(in)=b(in)-r*b(llb)
          end do
       end do
    end do
    r=a(ij)
    if(r .eq. tol) then
       ktilt=m
       return
    end if
    !
    !     End triangulation
    !
    !     Start back solution
    !
    piv=1./r
    do llb=m,nm,m
       x(llb)=b(llb)*piv
    end do
    i=m
    kk=ij
    do ii=1,m1
       kk=kk-i
       piv=1./a(kk)
       i=i-1
       do llb=i,nm,m
          in=llb
          r=b(in)
          ij=kk
          do j=i,m1
             ij=ij+j
             in=in+1
             r=r-a(ij)*x(in)
          end do
          x(llb)=r*piv
       end do
    end do
    !
    !     End solution
    !
  end subroutine relms
end module ModOdaKrig
