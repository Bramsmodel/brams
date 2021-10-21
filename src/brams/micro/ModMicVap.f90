!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModMicVap

  use rconstants, only: &
       alvi,   &
       cp, &
       alvl,   &
       cp253i, &
       cpi,    &
       cpi4,   &
       cpor,   &
       p00

  use ModMicControl, only:&
       MicControl

  implicit none

  private

  public :: thrmstr
  public :: vapdiff
  public :: vapflux
  public :: psxfer
  public :: newtemp
  public :: diffprep

contains



  subroutine thrmstr(m1,k1,k2,lpw,pp,thp,theta,pi0,rtp,rv,i,j,oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,i,j,k,lcat,lpw
    real :: fracliq,tcoal,tairstr
    integer, dimension(10) :: k1,k2
    real, dimension(m1) :: pp,thp,theta,pi0,rtp,rv

    do k = lpw,m1
       oneMicControl%pitot(k) = pi0(k) + pp(k)
       oneMicControl%press(k) = p00 * (oneMicControl%pitot(k) * cpi) ** cpor
       oneMicControl%tair(k) = theta(k) * oneMicControl%pitot(k) * cpi
    enddo

    do k = 1,k1(10)-1
       theta(k) = thp(k)
       rv(k) = rtp(k)
    enddo

    do k = k2(10)+1,m1
       theta(k) = thp(k)
       rv(k) = rtp(k)
    enddo

    do k = k1(10),k2(10)
       oneMicControl%til(k) = thp(k) * oneMicControl%pitot(k) * cpi
       oneMicControl%rliq(k) = 0.
       oneMicControl%rice(k) = 0.
    enddo

    do lcat = 1,2
       do k = k1(lcat),k2(lcat)
          oneMicControl%rliq(k) = oneMicControl%rliq(k) + oneMicControl%rx(k,lcat)
       enddo
    enddo

    do lcat = 3,5
       do k = k1(lcat),k2(lcat)
          oneMicControl%rice(k) = oneMicControl%rice(k) + oneMicControl%rx(k,lcat)
       enddo
    enddo

    do lcat = 6,7
       do k = k1(lcat),k2(lcat)
          call qtc(oneMicControl%qx(k,lcat),tcoal,fracliq)
          oneMicControl%rliq(k) = oneMicControl%rliq(k) + oneMicControl%rx(k,lcat) * fracliq
          oneMicControl%rice(k) = oneMicControl%rice(k) + oneMicControl%rx(k,lcat) * (1. - fracliq)
       enddo
    enddo

    do k = k1(10),k2(10)
       oneMicControl%qhydm(k) = alvl * oneMicControl%rliq(k) + alvi * oneMicControl%rice(k)
       oneMicControl%rvstr(k) = rtp(k) - oneMicControl%rliq(k) - oneMicControl%rice(k)
       oneMicControl%sa(k,1) = oneMicControl%til(k) * oneMicControl%qhydm(k) / (1.e-12 + oneMicControl%rliq(k) + oneMicControl%rice(k))
    enddo

    do k = k1(10),k2(10)
       if (oneMicControl%tair(k) .gt. 253.) then
          tairstr = 0.5 * (oneMicControl%til(k)  &
               + sqrt(oneMicControl%til(k) * (oneMicControl%til(k) + cpi4 * oneMicControl%qhydm(k))))
          oneMicControl%sa(k,1) = oneMicControl%sa(k,1) * cpi / (2. * tairstr - oneMicControl%til(k))
       else
          tairstr = oneMicControl%til(k) * (1. + oneMicControl%qhydm(k) * cp253i)
          oneMicControl%sa(k,1) = oneMicControl%sa(k,1) * cp253i
       endif
       oneMicControl%tairstrc(k) = tairstr - 273.16
       !LFR
       if(tairstr<100.0) then
          print *,'tairstr: ',tairstr,k,i,j
          print *,'til  qhydm ='  ,oneMicControl%til(k) , oneMicControl%qhydm(k) 
          stop 'mic_vap.f90 routine thrmstr'
       end if
    enddo

    return
  end subroutine thrmstr

  !******************************************************************************

  subroutine diffprep(m1,lcat,k1,k2,rv,dn0,i,j,mynum,oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,lcat,k1,k2,i,j,k,mynum,if1,if4,if6,if8,lhcat
    real :: fre,scdei
    real, dimension(m1) :: rv,dn0

    if (lcat .le. 2) then
       if1 = 1
       if4 = 4
       if6 = 6
       if8 = 8
    else
       if1 = 2
       if4 = 5
       if6 = 7
       if8 = 9
    endif

    do k = k1,k2
       lhcat = oneMicControl%jhcat(k,lcat)

       if (oneMicControl%rx(k,lcat) .lt. 1.e-9) go to 229

       fre = oneMicControl%frefac1(lhcat) * oneMicControl%emb(k,lcat) ** oneMicControl%pwmasi(lhcat)  &
            + oneMicControl%rdynvsci(k) * oneMicControl%frefac2(lhcat) * oneMicControl%emb(k,lcat) ** oneMicControl%cdp1(lhcat)

       oneMicControl%sb(k,lcat) = oneMicControl%cx(k,lcat) * dn0(k) * fre * oneMicControl%pi4dt
       oneMicControl%su(k,lcat) = oneMicControl%vapdif(k) * oneMicControl%sb(k,lcat)
       oneMicControl%sd(k,lcat) = oneMicControl%sh(k,lcat) * oneMicControl%rx(k,lcat)
       oneMicControl%se(k,lcat) = oneMicControl%su(k,lcat) * oneMicControl%sa(k,if6) + oneMicControl%sb(k,lcat) * oneMicControl%thrmcon(k)
       oneMicControl%sf(k,lcat) = oneMicControl%su(k,lcat) * oneMicControl%sl(if1) - oneMicControl%sb(k,lcat) * oneMicControl%sa(k,2)
       oneMicControl%sg(k,lcat) = oneMicControl%su(k,lcat) * oneMicControl%sa(k,if8) + oneMicControl%sb(k,lcat) * oneMicControl%sa(k,3)  &
            + oneMicControl%sj(lcat) * oneMicControl%qr(k,lcat)
       !     + lambda_j [Joules/kg_air added by radiative heating this timestep]
       scdei = 1. / (oneMicControl%sc(if1) * oneMicControl%sd(k,lcat) + oneMicControl%se(k,lcat))
       oneMicControl%ss(k,lcat) = oneMicControl%sf(k,lcat) * scdei
       oneMicControl%sw(k,lcat) = (oneMicControl%sg(k,lcat) - oneMicControl%sk(if1) * oneMicControl%sd(k,lcat)) * scdei
       oneMicControl%ttest(k,lcat) = oneMicControl%ss(k,lcat) * rv(k) + oneMicControl%sw(k,lcat)

229    continue

    enddo

    if (lcat .ge. 3 .and. lcat .le. 5) then
       do k = k1,k2
          if (oneMicControl%rx(k,lcat) .lt. 1.e-9) go to 228
          if (oneMicControl%ttest(k,lcat) .ge. 0.) then
             oneMicControl%sm(k,lcat) = 0.
             oneMicControl%sh(k,lcat) = 1.
             oneMicControl%sd(k,lcat) = oneMicControl%sh(k,lcat) * oneMicControl%rx(k,lcat)
             scdei = 1. / (oneMicControl%sc(if1) * oneMicControl%sd(k,lcat) + oneMicControl%se(k,lcat))
             oneMicControl%ss(k,lcat) = oneMicControl%sf(k,lcat) * scdei
             oneMicControl%sw(k,lcat) = (oneMicControl%sg(k,lcat) - oneMicControl%sk(if1) * oneMicControl%sd(k,lcat)) * scdei
          else
             oneMicControl%sm(k,lcat) = 1.
          endif
228       continue
       enddo
    endif

    if (lcat .ge. 6) then
       do k = k1,k2
          if (oneMicControl%rx(k,lcat) .lt. 1.e-9) go to 227
          if (oneMicControl%ttest(k,lcat) .ge. 0.) then
             oneMicControl%sm(k,lcat) = 0.
          else
             oneMicControl%sm(k,lcat) = 1.
          endif
227       continue
       enddo
    endif

    do k = k1,k2
       if (oneMicControl%rx(k,lcat) .lt. 1.e-9) go to 226
       oneMicControl%sy(k,lcat) = oneMicControl%rvsrefp(k,if1) * oneMicControl%sm(k,lcat) * oneMicControl%sw(k,lcat) - oneMicControl%sa(k,if4)
       oneMicControl%sz(k,lcat) = 1. - oneMicControl%rvsrefp(k,if1) * oneMicControl%ss(k,lcat) * oneMicControl%sm(k,lcat)
       oneMicControl%sumuy(k) = oneMicControl%sumuy(k) + oneMicControl%su(k,lcat) * oneMicControl%sy(k,lcat)
       oneMicControl%sumuz(k) = oneMicControl%sumuz(k) + oneMicControl%su(k,lcat) * oneMicControl%sz(k,lcat)

226    continue
    enddo

    return
  end subroutine diffprep

  !******************************************************************************

  subroutine vapdiff (m1,kf1,kf2,rv,i,j,mynum,oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,kf1,kf2,i,j,k,mynum
    real, dimension(m1) :: rv

    do k = kf1,kf2
       rv(k) = (oneMicControl%rvstr(k) + oneMicControl%sumuy(k)) / (1.0 + oneMicControl%sumuz(k))
    enddo

    return
  end subroutine vapdiff

  !******************************************************************************

  subroutine vapflux(m1,lcat,i,j,mynum,k1,k2,dn0,rv,oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,lcat,i,j,k,mynum,k1,k2,if1,if4
    real :: rxx
    real, dimension(m1) :: dn0,rv

    rxx=0.0 !LFR

    if (lcat .le. 2) then
       if1 = 1
       if4 = 4
    else
       if1 = 2
       if4 = 5
    endif

    do k = k1,k2

       if (oneMicControl%rx(k,lcat) .lt. 1.e-9) go to 229

       oneMicControl%tx(k,lcat) = (oneMicControl%ss(k,lcat) * rv(k) + oneMicControl%sw(k,lcat)) * oneMicControl%sm(k,lcat)
       oneMicControl%vap(k,lcat) = oneMicControl%su(k,lcat) * (rv(k) + oneMicControl%sa(k,if4) - oneMicControl%rvsrefp(k,if1) * oneMicControl%tx(k,lcat))

       if (oneMicControl%vap(k,lcat) .gt. -oneMicControl%rx(k,lcat)) then

          rxx = oneMicControl%rx(k,lcat) + oneMicControl%vap(k,lcat)

          if (oneMicControl%sm(k,lcat) .gt. .5) then
             oneMicControl%qx(k,lcat) = oneMicControl%sc(if1) * oneMicControl%tx(k,lcat) + oneMicControl%sk(if1)
             oneMicControl%qr(k,lcat) = oneMicControl%qx(k,lcat) * rxx
          else
             oneMicControl%qx(k,lcat) = (rv(k) * oneMicControl%sf(k,lcat) + oneMicControl%sg(k,lcat)  &
                  - oneMicControl%tx(k,lcat) * oneMicControl%se(k,lcat)) / oneMicControl%sd(k,lcat)
             oneMicControl%qx(k,lcat) = min(350000.,max(-100000.,oneMicControl%qx(k,lcat)))
             oneMicControl%qr(k,lcat) = oneMicControl%qx(k,lcat) * rxx
          endif

       endif

       !bob Now also do the following section if pristine ice totally melts:
       ! evaporate it too.

       if ((lcat .eq. 3 .and. oneMicControl%qx(k,lcat) .gt. 330000.) .or.  &
            oneMicControl%vap(k,lcat) .le. -oneMicControl%rx(k,lcat)) then

          oneMicControl%sumuy(k) = oneMicControl%sumuy(k) - oneMicControl%su(k,lcat) * oneMicControl%sy(k,lcat)
          oneMicControl%sumuz(k) = oneMicControl%sumuz(k) - oneMicControl%su(k,lcat) * oneMicControl%sz(k,lcat)
          oneMicControl%sumvr(k) = oneMicControl%sumvr(k) + oneMicControl%rx(k,lcat)
          rv(k) = (oneMicControl%rvstr(k) + oneMicControl%sumuy(k) + oneMicControl%sumvr(k)) / (1.0 + oneMicControl%sumuz(k))

          oneMicControl%vap(k,lcat) = - oneMicControl%rx(k,lcat)
          oneMicControl%tx(k,lcat) = 0.
          oneMicControl%rx(k,lcat) = 0.
          oneMicControl%qx(k,lcat) = 0.
          oneMicControl%qr(k,lcat) = 0.
       else
          oneMicControl%rx(k,lcat) = rxx
       endif

229    continue

    enddo
    return
  end subroutine vapflux

  !******************************************************************************

  subroutine psxfer(m1,k1,k2,dn0,i,j,oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,k1,k2,i,j,k,lhcat,it
    real :: embx,dn,xlim,dvap,dqr,dnum
    real, dimension(m1) :: dn0

    do k = k1,k2

       if (oneMicControl%vap(k,3) .gt. 0. .or. oneMicControl%vap(k,4) .lt. 0.) then

          if (oneMicControl%vap(k,3) .gt. 0.) then
             lhcat = oneMicControl%jhcat(k,3)
             embx = max(1.e-9,oneMicControl%rx(k,3)) / max(1.e-3,oneMicControl%cx(k,3))
             dn = oneMicControl%dnfac(lhcat) * embx ** oneMicControl%pwmasi(lhcat)
             it = nint(dn * 1.e6)

             !srf
             it=min(it, 5000)   
             xlim = oneMicControl%gam(it,3) * oneMicControl%dps2 * &
                  (oneMicControl%dps / dn) ** (oneMicControl%gnu(3) - 1.)  &
                  / (oneMicControl%gamn1(3) * oneMicControl%pwmas(lhcat) * dn ** 2)

             dvap = min(oneMicControl%rx(k,3),  &
                  oneMicControl%vap(k,3) * (xlim + oneMicControl%gam(it,1) / oneMicControl%gamn1(3)))
             dqr = dvap * oneMicControl%qx(k,3)
             dnum = dvap * min(oneMicControl%dpsmi(lhcat),1./embx)
          else
             lhcat = oneMicControl%jhcat(k,4)
             embx = max(1.e-9,oneMicControl%rx(k,4)) / max(1.e-3,oneMicControl%cx(k,4))
             dn = oneMicControl%dnfac(lhcat) * embx ** oneMicControl%pwmasi(lhcat)
             it = nint(dn * 1.e6)
             !srf
             it=min(it, 5000)   
             xlim = oneMicControl%gam(it,3) * oneMicControl%dps2 * &
                  (oneMicControl%dps / dn) ** (oneMicControl%gnu(4) - 1.)  &
                  / (oneMicControl%gamn1(4) * oneMicControl%pwmas(lhcat) * dn ** 2)

             dvap = max(-oneMicControl%rx(k,4),oneMicControl%vap(k,4) * xlim)
             dqr = dvap * oneMicControl%qx(k,4)
             dnum = dvap * max(oneMicControl%dpsmi(lhcat),1./embx)
          endif

          oneMicControl%rx(k,3) = oneMicControl%rx(k,3) - dvap
          oneMicControl%cx(k,3) = oneMicControl%cx(k,3) - dnum
          oneMicControl%qr(k,3) = oneMicControl%qr(k,3) - dqr
          oneMicControl%rx(k,4) = oneMicControl%rx(k,4) + dvap
          oneMicControl%cx(k,4) = oneMicControl%cx(k,4) + dnum
          oneMicControl%qr(k,4) = oneMicControl%qr(k,4) + dqr

       endif
    enddo
    return
  end subroutine psxfer

  !******************************************************************************

  subroutine newtemp(m1,kf1,kf2,rv,theta,i,j,oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    real rslf,rsif

    integer :: m1,kf1,kf2,i,j,k
    real, dimension(m1) :: rv,theta

    do k = kf1,kf2
       oneMicControl%tairc(k) = oneMicControl%tairstrc(k) + oneMicControl%sa(k,1) * (oneMicControl%rvstr(k) - rv(k))
       oneMicControl%tair(k)  = oneMicControl%tairc(k) + 273.16
       theta(k) = oneMicControl%tair(k) * cp / oneMicControl%pitot(k)

       oneMicControl%rvlsair(k) = rslf(oneMicControl%press(k),oneMicControl%tair(k))
       oneMicControl%rvisair(k) = rsif (oneMicControl%press(k),oneMicControl%tair(k))
       !LFR
       if(theta(k)<100.0) then
          print *,theta(k)
       end if
    enddo

    return
  end subroutine newtemp
end module ModMicVap
