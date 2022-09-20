!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!  FOR ICNFLG=3, DEBATING WHETHER TO KEEP IT #/M4 OR CHANGE
!  PARM TO #/KG/M.  NEED TO DEFINE AVMIPSA, ETC. FOR ALL CATEGORIES.
!  MAY WANT TO DEFINE C1 TOO AND RENAME IT.
!  IMPORTANT ISSUE: k loop limits for the jnmb == 5 sections
!  need to consider collection efficiencies for different habits?
!  collection efficiency for hail too high.  big hail should not
!  coallesce.
module ModMicrophysicsMisc

  use ModMicroFields, only: &
       MicroFields

  use ModMicControl, only: &
       MicControl
  
  use ModBasicFields, only: &
       BasicFields

  use mem_grid,  only: &
       grid_g,         & ! INTENT(IN)
       ngrid             ! INTENT(IN)

  use rconstants, only : &
       cpi, &
       alvl,             & ! INTENT(IN)
       alvi, &             ! INTENT(IN)
       pi4,              & ! INTENT(IN)
       alvl,             & ! INTENT(IN)
       alvi,             & ! INTENT(IN)
       alli                ! INTENT(IN)



  implicit none

  include "MicConstants.h"
  
  private
  public :: negadj1
  public :: ae1kmic
  public :: each_column
  public :: each_call
  public :: enemb
  public :: x02
  public :: pc03
  public :: sedim
  public :: c03
  public :: range_check
  
contains

  subroutine each_call(m1,dtlt, oneMicControl)
    ! Arguments
    integer, intent(in) :: m1
    real, intent(in)    :: dtlt
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables
    integer :: lcat, k

    ! Initialize constants for vapor diffusion and, for fixed diameter cases, emb.

    oneMicControl%colf = .785 * dtlt
    oneMicControl%pi4dt = pi4 * dtlt
    oneMicControl%sl(1) = alvl
    oneMicControl%sl(2) = alvi
    oneMicControl%sc(1) = 4186.
    oneMicControl%sc(2) = 2093.
    oneMicControl%sj(1) = 0
    oneMicControl%sj(2) = 1
    oneMicControl%sj(3) = 0
    oneMicControl%sj(4) = 0
    oneMicControl%sj(5) = 0
    oneMicControl%sj(6) = 1
    oneMicControl%sj(7) = 1
    oneMicControl%sk(1) = alli
    oneMicControl%sk(2) = 0.

    do lcat = 1,7
       if (oneMicControl%jnmb(lcat) == 2) then
          do k = 2,m1-1
             oneMicControl%emb(k,lcat) = oneMicControl%cfmas(lcat) * oneMicControl%parm(lcat) ** oneMicControl%pwmas(lcat)
          enddo
       endif
       do k = 2,m1-1
          oneMicControl%jhcat(k,lcat) = lcat
       enddo
    enddo

    do k = 2,m1-1
       oneMicControl%sh(k,1) = 0.
       oneMicControl%sh(k,2) = 1.
       oneMicControl%sh(k,6) = 1.
       oneMicControl%sh(k,7) = 1.

       oneMicControl%sm(k,1) = 1.
       oneMicControl%sm(k,2) = 1.
    enddo

    return
  end subroutine each_call

  !******************************************************************************

  subroutine range_check(m1,k1,k2,k3,i,j,lpw_R, oneMicroFields, oneMicControl)
    ! Arguments
    integer, intent(in)                 :: m1, i, j
    integer, dimension(10), intent(out) :: k1, k2, k3
    real, intent(in) :: lpw_R
    type(MicroFields), pointer, intent(in)       :: oneMicroFields 
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables
    integer                             :: k, lcatt, lcat, l, jcat, lpw

    lpw=int(lpw_R)
    ! zero out microphysics scratch arrays for the present i,j column

    do lcat = 1,ncat
       do k = 2,m1-1
          oneMicControl%rx(k,lcat) = 0.
          oneMicControl%cx(k,lcat) = 0.
          oneMicControl%qr(k,lcat) = 0.
          oneMicControl%qx(k,lcat) = 0.
          oneMicControl%vap(k,lcat) = 0.
          oneMicControl%tx(k,lcat) = 0.
       enddo

       if (oneMicControl%jnmb(lcat) >= 3) then
          do k = 2,m1-1
             oneMicControl%emb(k,lcat) = 0.
          enddo
       endif

       do jcat = 1,ncat
          do k = 2,m1-1
             oneMicControl%rxfer(k,lcat,jcat) = 0.
             oneMicControl%qrxfer(k,lcat,jcat) = 0.
             oneMicControl%enxfer(k,lcat,jcat) = 0.
          enddo
       enddo
    enddo

    do l = 1,7
       k1(l) = lpw
       k2(l) = 1
    enddo

    ! fill scratch arrays for cloud water

    if (oneMicControl%jnmb(1) >= 1) then
       do k = lpw,m1-1
          if (oneMicroFields%rcp(k,i,j) >= 1.e-9) then
             k2(1) = k
             oneMicControl%rx(k,1) = oneMicroFields%rcp(k,i,j)
             if (oneMicControl%jnmb(1) >= 5) oneMicControl%cx(k,1) = oneMicroFields%ccp(k,i,j)

             !
             !Carrio 2012 -------------------------------------------
             ! ALL two moment options require this array!!!
             if (oneMicControl%jnmb(1) >=5) oneMicControl%cccnx(k) = oneMicroFields%cccnp(k,i,j)
             !-------------------------------------------------------
             !      

             !print*,"xccnp=",oneMicControl%cccnx(k), oneMicControl%cx(k,1)


             !--(DMK-CARRIO-OLD)------------------------------------------------------------
             !           if (oneMicControl%jnmb(1) == 7) oneMicControl%cccnx(k) = oneMicroFields%cccnp(k,i,j)
             !--(DMK-CARRIO-INI)------------------------------------------------------------

          else
             if (k2(1) == 1) k1(1) = k + 1
          endif
       enddo
    endif

    ! fill scratch arrays for rain

    if (oneMicControl%jnmb(2) >= 1) then
       do k = lpw,m1-1
          if (oneMicroFields%rrp(k,i,j) >= 1.e-9) then
             k2(2) = k
             oneMicControl%rx(k,2) = oneMicroFields%rrp(k,i,j)
             oneMicControl%qx(k,2) = oneMicroFields%q2(k,i,j)
             oneMicControl%qr(k,2) = oneMicControl%qx(k,2) * oneMicControl%rx(k,2)
             if (oneMicControl%jnmb(2) >= 5) oneMicControl%cx(k,2) = oneMicroFields%crp(k,i,j)
          else
             if (k2(2) == 1) k1(2) = k + 1
          endif
       enddo
    endif

    ! fill scratch arrays for pristine ice

    if (oneMicControl%jnmb(3) >= 1) then
       do k = lpw,m1-1
          if (oneMicroFields%rpp(k,i,j) >= 1.e-9) then
             k2(3) = k
             oneMicControl%rx(k,3) = oneMicroFields%rpp(k,i,j)
             oneMicControl%cx(k,3) = oneMicroFields%cpp(k,i,j)
             if (oneMicControl%jnmb(3) == 7) oneMicControl%cifnx(k) = oneMicroFields%cifnp(k,i,j)
          else
             if (k2(3) == 1) k1(3) = k + 1
          endif
       enddo
    endif

    ! fill scratch arrays for snow

    if (oneMicControl%jnmb(4) >= 1) then
       do k = lpw,m1-1
          if (oneMicroFields%rsp(k,i,j) >= 1.e-9) then
             k2(4) = k
             oneMicControl%rx(k,4) = oneMicroFields%rsp(k,i,j)
             if (oneMicControl%jnmb(4) >= 5) oneMicControl%cx(k,4) = oneMicroFields%csp(k,i,j)
          else
             if (k2(4) == 1) k1(4) = k + 1
          endif
       enddo
    endif

    ! fill scratch arrays for aggregates

    if (oneMicControl%jnmb(5) >= 1) then
       do k = lpw,m1-1
          if (oneMicroFields%rap(k,i,j) >= 1.e-9) then
             k2(5) = k
             oneMicControl%rx(k,5) = oneMicroFields%rap(k,i,j)
             if (oneMicControl%jnmb(5) >= 5) oneMicControl%cx(k,5) = oneMicroFields%cap(k,i,j)
          else
             if (k2(5) == 1) k1(5) = k + 1
          endif
       enddo
    endif

    ! fill scratch arrays for graupel

    if (oneMicControl%jnmb(6) >= 1) then
       do k = lpw,m1-1
          if (oneMicroFields%rgp(k,i,j) >= 1.e-9) then
             k2(6) = k
             oneMicControl%rx(k,6) = oneMicroFields%rgp(k,i,j)
             oneMicControl%qx(k,6) = oneMicroFields%q6(k,i,j)
             oneMicControl%qr(k,6) = oneMicControl%qx(k,6) * oneMicControl%rx(k,6)
             if (oneMicControl%jnmb(6) >= 5) oneMicControl%cx(k,6) = oneMicroFields%cgp(k,i,j)
          else
             if (k2(6) == 1) k1(6) = k + 1
          endif
       enddo
    endif

    ! fill scratch arrays for hail

    if (oneMicControl%jnmb(7) >= 1) then
       do k = lpw,m1-1
          if (oneMicroFields%rhp(k,i,j) >= 1.e-9) then
             k2(7) = k
             oneMicControl%rx(k,7) = oneMicroFields%rhp(k,i,j)
             oneMicControl%qx(k,7) = oneMicroFields%q7(k,i,j)
             oneMicControl%qr(k,7) = oneMicControl%qx(k,7) * oneMicControl%rx(k,7)
             if (oneMicControl%jnmb(7) >= 5) oneMicControl%cx(k,7) = oneMicroFields%chp(k,i,j)
          else
             if (k2(7) == 1) k1(7) = k + 1
          endif
       enddo
    endif

    k3(1) = k2(1)
    k3(3) = k2(3)

    k1(8) = min(k1(1),k1(2))
    k2(8) = max(k2(1),k2(2))
    k1(9) = min(k1(3),k1(4),k1(5),k1(6),k1(7))
    k2(9) = max(k2(3),k2(4),k2(5),k2(6),k2(7))
    k1(10) = min(k1(8),k1(9))
    k2(10) = max(k2(8),k2(9))

    return
  end subroutine range_check

  !******************************************************************************

  subroutine each_column(m1,k1,k2,i,j,lpw,rv,dn0, oneMicControl)
    ! Arguments
    integer, intent(in)                :: m1
    integer, intent(in)                :: i, j  ! Not used
    integer, intent(in)                :: lpw
    integer, dimension(10), intent(in) :: k1, k2
    real, dimension(m1), intent(in)    :: rv, dn0
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables
    integer :: k,nt,ns
    real :: ck1,ck2,ck3,elsref,elsrefp,dplinv,eisref,eisrefp,dpiinv,relhum

    real :: rslf,rsif,eslf,eslpf,esif,esipf

    data ck1,ck2,ck3/-4.818544e-3,1.407892e-4,-1.249986e-7/

    do k = lpw,m1-1
       oneMicControl%rvlsair(k) = rslf (oneMicControl%press(k),oneMicControl%tair(k))
       oneMicControl%rvisair(k) = rsif (oneMicControl%press(k),oneMicControl%tair(k))
       oneMicControl%dn0i(k) = 1. / dn0(k)
       oneMicControl%tairc(k)   = oneMicControl%tair(k) - 273.16
       oneMicControl%tx(k,1) = oneMicControl%tairc(k)
       oneMicControl%thrmcon(k) = ck1 + (ck2 + ck3 * oneMicControl%tair(k)) * oneMicControl%tair(k)
       oneMicControl%dynvisc(k) = .1718e-4 + .49e-7 * oneMicControl%tairc(k)

       ! Diagnose habit of pristine ice and snow

       nt = max(1,min(31,-nint(oneMicControl%tairc(k))))
       relhum = min(1.,rv(k) / oneMicControl%rvlsair(k))
       ns = max(1,nint(100. * relhum))
       oneMicControl%jhcat(k,3) = oneMicControl%jhabtab(nt,ns,1)
       oneMicControl%jhcat(k,4) = oneMicControl%jhabtab(nt,ns,2)

    enddo

    do k = k1(10),k2(10)
       oneMicControl%vapdif(k)     = 2.14 * (oneMicControl%tair(k) / 273.15) ** 1.94 / oneMicControl%press(k)
       oneMicControl%rdynvsci(k) = sqrt(1. / oneMicControl%dynvisc(k))
       oneMicControl%denfac(k) = sqrt(oneMicControl%dn0i(k))

       oneMicControl%colfacr(k) = oneMicControl%colf * oneMicControl%denfac(k) * dn0(k)
       oneMicControl%colfacr2(k) = 2. * oneMicControl%colfacr(k)
       oneMicControl%colfacc(k) = oneMicControl%colfacr(k) * dn0(k)
       oneMicControl%colfacc2(k) = 2. * oneMicControl%colfacc(k)

       oneMicControl%tref(k,1)   = oneMicControl%tairc(k) - min(25.,700. * (oneMicControl%rvlsair(k) - rv(k)))
       oneMicControl%sa(k,2) = oneMicControl%thrmcon(k) * oneMicControl%sa(k,1)
       oneMicControl%sa(k,3) = oneMicControl%thrmcon(k) * (oneMicControl%tairstrc(k) + oneMicControl%sa(k,1) * oneMicControl%rvstr(k))

       oneMicControl%sumuy(k) = 0.
       oneMicControl%sumuz(k) = 0.
       oneMicControl%sumvr(k) = 0.
    enddo

    do k = k1(8),k2(8)
       elsref       = eslf(oneMicControl%tref(k,1))
       elsrefp      = eslpf(oneMicControl%tref(k,1))
       dplinv       = 1. / (oneMicControl%press(k) - elsref)
       oneMicControl%rvsref (k,1) = .622 * elsref * dplinv
       oneMicControl%rvsrefp(k,1) = .622 * elsrefp * dplinv * (1. + elsref * dplinv)

       oneMicControl%sa(k,4) = oneMicControl%rvsrefp(k,1) * oneMicControl%tref(k,1) - oneMicControl%rvsref(k,1)
       oneMicControl%sa(k,6) = alvl * oneMicControl%rvsrefp(k,1)
       oneMicControl%sa(k,8) = alvl * oneMicControl%sa(k,4)
    enddo

    do k = k1(9),k2(9)
       oneMicControl%tref(k,2)    = min(0.,oneMicControl%tref(k,1))
       eisref       = esif (oneMicControl%tref(k,2))
       eisrefp      = esipf(oneMicControl%tref(k,2))
       dpiinv       = 1. / (oneMicControl%press(k) - eisref)
       oneMicControl%rvsref (k,2) = .622 * eisref * dpiinv
       oneMicControl%rvsrefp(k,2) = .622 * eisrefp * dpiinv * (1. + eisref * dpiinv)
       oneMicControl%rvs0(k)      = 379.4 / (oneMicControl%press(k) - 610.)

       oneMicControl%sa(k,5) = oneMicControl%rvsrefp(k,2) * oneMicControl%tref(k,2) - oneMicControl%rvsref(k,2)
       oneMicControl%sa(k,7) = alvi * oneMicControl%rvsrefp(k,2)
       oneMicControl%sa(k,9) = alvi * oneMicControl%sa(k,5)
       oneMicControl%sh(k,3) = 0.
       oneMicControl%sh(k,4) = 0.
       oneMicControl%sh(k,5) = 0.

    enddo

    return
  end subroutine each_column

  !******************************************************************************

  subroutine enemb(m1,k1,k2,lcat,jflag,dn0,i,j, oneMicControl)
    ! Arguments
    integer, intent(in)             :: m1, k1, k2, lcat, jflag
    integer, intent(in)             :: i, j ! Not used
    real, dimension(m1), intent(in) :: dn0
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables
    integer :: k,lhcat
    real :: embi,parmi,fracmass,cxloss


    if (oneMicControl%jnmb(lcat) == 2) then
       embi = 1. / oneMicControl%emb(2,lcat)
       do k = k1,k2
          oneMicControl%cx(k,lcat) = oneMicControl%rx(k,lcat) * embi
       enddo
    elseif (oneMicControl%jnmb(lcat) == 3) then
       do k = k1,k2
          lhcat = oneMicControl%jhcat(k,lcat)
          oneMicControl%emb(k,lcat) = oneMicControl%cfemb0(lhcat) * (dn0(k) * oneMicControl%rx(k,lcat)) ** oneMicControl%pwemb0(lhcat)
          oneMicControl%cx(k,lcat) = oneMicControl%cfen0(lhcat) * oneMicControl%dn0i(k)  &
               * (dn0(k) * oneMicControl%rx(k,lcat)) ** oneMicControl%pwen0(lhcat)
       enddo
    elseif (oneMicControl%jnmb(lcat) == 4) then
       parmi = 1. / oneMicControl%parm(lcat)
       do k = k1,k2
          oneMicControl%emb(k,lcat) = max(oneMicControl%emb0(lcat),min(oneMicControl%emb1(lcat),oneMicControl%rx(k,lcat) * parmi))
          oneMicControl%cx(k,lcat) = oneMicControl%rx(k,lcat) / oneMicControl%emb(k,lcat)
       enddo
    elseif (oneMicControl%jnmb(lcat) >= 5 .and. jflag == 1) then
       do k = k1,k2
          oneMicControl%emb(k,lcat) = max(oneMicControl%emb0(lcat),min(oneMicControl%emb1(lcat),oneMicControl%rx(k,lcat)  &
               / max(1.e-9,oneMicControl%cx(k,lcat))))
          oneMicControl%cx(k,lcat) = oneMicControl%rx(k,lcat) / oneMicControl%emb(k,lcat)
       enddo
    elseif (oneMicControl%jnmb(lcat) >= 5 .and. jflag == 2) then
       do k = k1,k2

          if (oneMicControl%rx(k,lcat) >= 1.e-9) then

             if (oneMicControl%vap(k,lcat) < 0.) then
                fracmass = min(1.,-oneMicControl%vap(k,lcat) / oneMicControl%rx(k,lcat))
                cxloss = oneMicControl%cx(k,lcat) * oneMicControl%enmlttab(int(200. * fracmass) + 1  &
                     ,oneMicControl%jhcat(k,lcat))
                oneMicControl%cx(k,lcat) = oneMicControl%cx(k,lcat) - cxloss
             endif
             oneMicControl%emb(k,lcat) = max(oneMicControl%emb0(lcat),min(oneMicControl%emb1(lcat),oneMicControl%rx(k,lcat)  &
                  / max(1.e-9,oneMicControl%cx(k,lcat))))
             oneMicControl%cx(k,lcat) = oneMicControl%rx(k,lcat) / oneMicControl%emb(k,lcat)

          endif

       enddo
    endif

    return
  end subroutine enemb

  !******************************************************************************

  subroutine x02(m1,k1,k2,lcat,dn0,i,j, oneMicControl)
    ! Arguments
    integer, intent(in)                :: m1, lcat
    integer, intent(in)                :: i, j  ! Not Used
    integer, dimension(10),intent(inout) :: k1, k2
    real, dimension(m1), intent(in)    :: dn0
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables
    integer :: k,jflag,lhcat,inc,idns
    real :: rinv,closs,rxinv,rmelt,fracliq,cmelt,tcoal,ricetor6,rshed,rmltshed  &
         ,qrmltshed,shedmass,fracmloss,dn


    k1(lcat) = k1(10)
    k2(lcat) = 1
    do k = k1(10),k2(10)
       if (oneMicControl%rx(k,lcat) >= 1.e-9) k2(lcat) = k
       if (k2(lcat) == 1 .and. oneMicControl%rx(k,lcat) < 1.e-9) k1(lcat) = k + 1
    enddo

    if (lcat == 2 .or. lcat >= 4) then
       jflag = 1

       call enemb(m1,k1(lcat),k2(lcat),lcat,jflag,dn0,i,j, oneMicControl)

       do k = k1(lcat),k2(lcat)

          if (oneMicControl%rx(k,lcat) >= 1.e-9) then

             lhcat = oneMicControl%jhcat(k,lcat)
             oneMicControl%vterm(k,lcat) = -oneMicControl%vtfac(lhcat) * oneMicControl%emb(k,lcat) ** oneMicControl%pwvtmasi(lhcat) &
                  * oneMicControl%denfac(k)

          endif

       enddo
    endif

    if (lcat == 2) then

       do k = k1(lcat),k2(lcat)

          if (oneMicControl%rx(k,lcat) >= 1.e-9) then

             rxinv = 1. / oneMicControl%rx(k,lcat)
             oneMicControl%qx(k,lcat) = oneMicControl%qr(k,lcat) * rxinv
             ! limit rain to under 48C and over -80C
             oneMicControl%qx(k,lcat) = max(0.,min(1.6*alli,oneMicControl%qx(k,lcat)))

          endif

       enddo

    elseif (lcat == 3) then

       do k = k1(lcat),k2(lcat)

          if (oneMicControl%rx(k,lcat) >= 1.e-9) then


             rinv = 1. / oneMicControl%rx(k,lcat)
             oneMicControl%qx(k,lcat) = oneMicControl%qr(k,lcat) * rinv

             call qtc(oneMicControl%qx(k,lcat),tcoal,fracliq)

             rmelt = oneMicControl%rx(k,lcat) * fracliq
             cmelt = oneMicControl%cx(k,lcat) * fracliq

             oneMicControl%rx(k,lcat) = oneMicControl%rx(k,lcat) - rmelt
             oneMicControl%rx(k,1) = oneMicControl%rx(k,1) + rmelt
             oneMicControl%cx(k,lcat) = oneMicControl%cx(k,lcat) - cmelt
             oneMicControl%cx(k,1) = oneMicControl%cx(k,1) + cmelt


          endif

       enddo
       !
       ! meyers - source for cloud aerosol number here?
       !
    elseif (lcat == 4 .or. lcat == 5) then

       do k = k1(lcat),k2(lcat)

          if (oneMicControl%rx(k,lcat) >= 1.e-9) then

             rinv = 1. / oneMicControl%rx(k,lcat)
             oneMicControl%qx(k,lcat) = oneMicControl%qr(k,lcat) * rinv
             call qtc(oneMicControl%qx(k,lcat),tcoal,fracliq)

             if (fracliq > 1.e-6) then
                rmelt = oneMicControl%rx(k,lcat) * fracliq

                ! change this??? move to rain instead ??? look at melting decisions in col2

                ricetor6 = min(oneMicControl%rx(k,lcat) - rmelt,rmelt)
                oneMicControl%rx(k,lcat) = oneMicControl%rx(k,lcat) - rmelt - ricetor6
                oneMicControl%rx(k,6) = oneMicControl%rx(k,6) + rmelt + ricetor6
                oneMicControl%qr(k,6) = oneMicControl%qr(k,6) + rmelt * alli
                oneMicControl%qx(k,lcat) = 0.

                ! keep the above the same with ricetor6
                ! meyers - use sa melt table here? yes
                !
                fracmloss = (rmelt + ricetor6) * rinv
                closs = oneMicControl%enmlttab(int(200. * fracmloss) + 1,oneMicControl%jhcat(k,lcat)) * oneMicControl%cx(k,lcat)
                oneMicControl%cx(k,lcat) = oneMicControl%cx(k,lcat) - closs
                oneMicControl%cx(k,6) = oneMicControl%cx(k,6) + closs
             endif


          endif

       enddo


    elseif (lcat == 6) then

       do k = k1(lcat),k2(lcat)

          if (oneMicControl%rx(k,lcat) >= 1.e-9) then

             rxinv = 1. / oneMicControl%rx(k,lcat)
             oneMicControl%qx(k,lcat) = oneMicControl%qr(k,lcat) * rxinv
             call qtc(oneMicControl%qx(k,lcat),tcoal,fracliq)

             if (fracliq > 0.95) then
                oneMicControl%rx(k,2) = oneMicControl%rx(k,2) + oneMicControl%rx(k,6)
                oneMicControl%qr(k,2) = oneMicControl%qr(k,2) + oneMicControl%rx(k,6) * alli
                oneMicControl%cx(k,2) = oneMicControl%cx(k,2) + oneMicControl%cx(k,6)
                oneMicControl%rx(k,6) = 0.
                oneMicControl%qr(k,6) = 0.
                oneMicControl%cx(k,6) = 0.
             endif

          endif

       enddo

    elseif (lcat == 7) then

       shedmass = 5.236e-7
       do k = k1(lcat),k2(lcat)

          if (oneMicControl%rx(k,lcat) >= 1.e-9) then

             rxinv = 1. / oneMicControl%rx(k,lcat)
             oneMicControl%qx(k,lcat) = oneMicControl%qr(k,lcat) * rxinv
             !c          oneMicControl%qx(k,lcat) = max(-50.,oneMicControl%qx(k,lcat))
             call qtc(oneMicControl%qx(k,lcat),tcoal,fracliq)

             if (fracliq > 0.95) then
                oneMicControl%rx(k,2) = oneMicControl%rx(k,2) + oneMicControl%rx(k,7)
                oneMicControl%qr(k,2) = oneMicControl%qr(k,2) + oneMicControl%rx(k,7) * alli
                oneMicControl%cx(k,2) = oneMicControl%cx(k,2) + oneMicControl%cx(k,7)
                oneMicControl%rx(k,7) = 0.
                oneMicControl%qr(k,7) = 0.
                oneMicControl%cx(k,7) = 0.
                !
                !  take out following IF statement?
                !

             elseif (fracliq > 0.3) then

                lhcat = oneMicControl%jhcat(k,lcat)
                inc = nint(200. * fracliq) + 1
                dn = oneMicControl%dnfac(lhcat) * oneMicControl%emb(k,lcat) ** oneMicControl%pwmasi(lhcat)
                idns = max(1,nint(1.e3 * dn * oneMicControl%gnu(lcat)))
                rshed = oneMicControl%rx(k,lcat) * oneMicControl%shedtab(inc,idns)
                !cc               rmltshed = oneMicControl%rx(k,lcat) * rmlttab(inc) + rshed
                rmltshed = rshed
                qrmltshed = rmltshed * alli

                oneMicControl%rx(k,2) = oneMicControl%rx(k,2) + rmltshed
                oneMicControl%qr(k,2) = oneMicControl%qr(k,2) + qrmltshed
                oneMicControl%rx(k,lcat) = oneMicControl%rx(k,lcat) - rmltshed
                oneMicControl%qr(k,lcat) = oneMicControl%qr(k,lcat) - qrmltshed
                !               closs = oneMicControl%cx(k,lcat) * oneMicControl%enmlttab(inc,lhcat)
                !               oneMicControl%cx(k,lcat) = oneMicControl%cx(k,lcat) - closs
                !               oneMicControl%cx(k,2) = oneMicControl%cx(k,2) + closs + rshed / shedmass
                oneMicControl%cx(k,2) = oneMicControl%cx(k,2) + rshed / shedmass
             endif

          endif

       enddo

    endif
    return
  end subroutine x02

  !******************************************************************************

  subroutine c03(m1,k1,k2,lcat,dn0,i,j, oneMicControl)
    ! Arguments
    integer, intent(in)             :: m1, k1, k2, lcat, i, j
    real, dimension(m1), intent(in) :: dn0
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables
    integer :: jflag

    jflag = 1
    if (oneMicControl%jnmb(lcat) >= 3) call enemb(m1,k1,k2,lcat,jflag,dn0,i,j, oneMicControl)

    return
  end subroutine c03

  !******************************************************************************

  subroutine pc03(m1,k1,k2,lcat,dn0,i,j, oneMicControl)
    ! Arguments
    integer, intent(in)             :: m1, k1, k2, lcat, i, j
    real, dimension(m1), intent(in) :: dn0
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables
    integer :: k, lhcat, jflag

    jflag = 1
    if (oneMicControl%jnmb(lcat) >= 3) call enemb(m1,k1,k2,lcat,jflag,dn0,i,j, oneMicControl)

    do k = k1,k2

       if (oneMicControl%rx(k,lcat) >= 1.e-9) then

          lhcat = oneMicControl%jhcat(k,lcat)
          oneMicControl%vterm(k,lcat) = -oneMicControl%vtfac(lhcat) * oneMicControl%emb(k,lcat) ** oneMicControl%pwvtmasi(lhcat) * oneMicControl%denfac(k)

       endif

    enddo

    return
  end subroutine pc03

  !******************************************************************************

  subroutine sedim(m1,lcat,ngr,nembfall,maxkfall,k1,k2,lpw,i,j  &
       ,rtp,thp,theta,dn0,alphasfc  &
       ,pcpg,qpcpg,dpcpg,dtlti,cnew,rnew,qrnew,pcpfillc,pcpfillr,sfcpcp, &
       dzt, if_adap, oneMicControl)
    ! Arguments
    integer, intent(in) :: m1, lcat, ngr, nembfall, maxkfall, k1, k2, lpw, i, j
    real, dimension(m1), intent(inout) :: rtp
    real, dimension(m1), intent(in)  :: thp, theta, dn0
    real, intent(in) :: alphasfc, dtlti
    real, intent(inout) :: pcpg, qpcpg, dpcpg
    real, dimension(m1), intent(out) :: cnew, rnew, qrnew
    real, dimension(m1,maxkfall,nembfall,nhcat), intent(in) :: pcpfillc, pcpfillr
    real, dimension(maxkfall,nembfall,nhcat), intent(in) :: sfcpcp
    real, dimension(m1), intent(in) :: dzt          ! From RAMS 6.0
    integer, intent(in) :: if_adap                  ! From RAMS 6.0
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables
    integer :: k,lhcat,iemb,iemb2,kkf,kk
    real :: colddn0,rolddn0,qrolddn0,dispemb,riemb,wt2,psfc,qnew

    oneMicControl%pcprx(lcat) = 0.

    !--(DMK-CARRIO-INI)-----------------------------------------------------
    do k = 2,k2  ! Old way
       !do k = 1, m1 ! From RAMS 6.0
       !--(DMK-CARRIO-OLD)-----------------------------------------------------
       !  !do k = 2,k2
       !  do k = 1, m1 ! From RAMS 6.0
       !--(DMK-CARRIO-END)-----------------------------------------------------
       rnew(k) = 0.
       cnew(k) = 0.
       qrnew(k) = 0.
    enddo

    !--(DMK-CARRIO-INI)-----------------------------------------------------
    !change_MP for chem
    do k =1,m1
       oneMicControl%rsedim(k) = 0.
    enddo
    !end change_MP
    !--(DMK-CARRIO-FIM)-----------------------------------------------------

    do k = k1,k2
       lhcat = oneMicControl%jhcat(k,lcat)

       if (oneMicControl%rx(k,lcat) > 1.e-9) then
          colddn0 = oneMicControl%cx(k,lcat) * dn0(k)
          rolddn0 = oneMicControl%rx(k,lcat) * dn0(k)
          qrolddn0 = oneMicControl%qx(k,lcat) * rolddn0

          dispemb = oneMicControl%ch1(lhcat)  &
               * (oneMicControl%emb(k,lcat)/oneMicControl%cfmas(lhcat)) ** oneMicControl%ch3(lhcat) * sqrt(oneMicControl%dn0i(k))
          riemb = 1. + oneMicControl%ch2(lhcat,ngr) * log10(dispemb / oneMicControl%dispemb0(lhcat,ngr))

          !Bob (10/24/00):  Now, limiting iemb to max of nembfall

          iemb = min(nint(riemb),nembfall)

          if (k <= maxkfall) then
             psfc = rolddn0 * sfcpcp(k,iemb,lhcat)
          endif

          do kkf = 1,min(maxkfall,k-1)
             kk = k + 1 - kkf

             cnew(kk)  = cnew(kk)   &
                  +  colddn0 * oneMicControl%dn0i(kk) * pcpfillc(k,kkf,iemb,lhcat)
             rnew(kk)  = rnew(kk)   &
                  +  rolddn0 * oneMicControl%dn0i(kk) * pcpfillr(k,kkf,iemb,lhcat)
             qrnew(kk) = qrnew(kk)  &
                  + qrolddn0 * oneMicControl%dn0i(kk) * pcpfillr(k,kkf,iemb,lhcat)

             !---------------------------------------------------------------
             ! adjustment for underground and partial grid cells for ada grid

             !if (ada_flag == 1) then
             !   do k = 2,lpw-1

             !must consider volt correction to pcpfillc and pcpfillr tables

             !---------------------------------------------------------------

          enddo

          if (k <= maxkfall) then
             qpcpg = qpcpg + psfc * oneMicControl%qx(k,lcat)
             oneMicControl%pcprx(lcat) = oneMicControl%pcprx(lcat) + psfc
          endif

       endif
    enddo

    ! From RAMS 6.0
    if (if_adap == 1) then
       do k = 2,lpw-1
          oneMicControl%pcprx(lcat) = oneMicControl%pcprx(lcat) + rnew(k) * dn0(k) / dzt(k)
          qpcpg = qpcpg + qrnew(k) * dn0(k) / dzt(k)

          cnew(k) = 0.
          rnew(k) = 0.
          qrnew(k) = 0.
       enddo
    endif

    pcpg = pcpg + oneMicControl%pcprx(lcat)
    oneMicControl%accpx(lcat) = oneMicControl%pcprx(lcat)
    dpcpg = dpcpg + oneMicControl%pcprx(lcat) * alphasfc
    oneMicControl%pcprx(lcat) = oneMicControl%pcprx(lcat) * dtlti

    !do k = 2,k2
    do k =lpw,k2 ! From RAMS 6.0
       rtp(k) = rtp(k) + rnew(k) - oneMicControl%rx(k,lcat)
       qnew = qrnew(k) / max(1.e-20, rnew(k))

       !         if (iqflag == 1) then
       oneMicControl%tairc(k) = oneMicControl%tairc(k) - thp(k) * thp(k)  &
            * (2820. * (rnew(k) - oneMicControl%rx(k,lcat))  &
            - cpi * (qrnew(k) - oneMicControl%qx(k,lcat) * oneMicControl%rx(k,lcat)))  &
            / (max(oneMicControl%tair(k), 253.) * theta(k))
       !         else
       !            oneMicControl%tairc(k) = oneMicControl%tairc(k) - thp(k) * thp(k) * 2820.
       !     +          * (rnew(k) - oneMicControl%rx(k,lcat)) / (max(oneMicControl%tair(k), 253.) * theta(k))
       !         endif

       !--(DMK-CARRIO-INI)------------------------------------------------------------------
       !change_MP rsedim is the fraction of rain removed or added by sedimentation
       if(lcat.eq.2 .and. oneMicControl%rx(k,2).gt.1.e-9)then
          oneMicControl%rsedim(k)=-(rnew(k)-oneMicControl%rx(k,lcat))/oneMicControl%rx(k,lcat)
       endif
       !end change_MP
       !--(DMK-CARRIO-FIM)------------------------------------------------------------------

       oneMicControl%rx(k,lcat) = rnew(k)
       oneMicControl%cx(k,lcat) = cnew(k)
       oneMicControl%qx(k,lcat) = qnew

       if (oneMicControl%rx(k,lcat) < 1.e-9) then
          oneMicControl%rx(k,lcat) = 0.
          oneMicControl%cx(k,lcat) = 0.
          oneMicControl%qx(k,lcat) = 0.
       endif

    enddo
    return
  end subroutine sedim

  !******************************************************************************

  subroutine negadj1(m1,m2,m3,oneBasicFields, oneMicControl, oneMicroFields)
    ! Arguments
    integer, intent(in) :: m1, m2, m3
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicControl), pointer, intent(in) :: oneMicControl
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    real :: vctr9(m1)

    if (oneMicControl%level == 0) return

    call adj1(m1,m2,m3,grid_g(ngrid)%lpw,oneBasicFields%rtp  &
         ,oneBasicFields%thp,oneMicroFields,vctr9,oneMicControl)

    return
  end subroutine negadj1

  !******************************************************************************

  subroutine adj1(m1,m2,m3,lpw_R,rtp,thp,oneMicroFields,vctr9, oneMicControl)
    ! Arguments
    integer, intent(in)                      :: m1, m2, m3
    real, dimension(m2,m3), intent(in)    :: lpw_R
    type(MicroFields), pointer, intent(in)         :: oneMicroFields
    real, dimension(m1,m2,m3), intent(inout) :: rtp, thp
    real, dimension(m1), intent(out)         :: vctr9
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables
    integer :: i,j,k,lcat,ka
    real :: frac
    integer, dimension(m2,m3)  :: lpw

    lpw=int(lpw_R)

    if (oneMicControl%level .eq. 0) return

    do lcat = 1,ncat
       do k = 1,m1
          oneMicControl%rx(k,lcat) = 0.
          oneMicControl%cx(k,lcat) = 0.
       enddo
    enddo

    do j = 1,m3
       do i = 1,m2

          !--(DMK-CARRIO-INI)-----------------------------------------------------
          !srf- rams60 mod
          ! Do this for all levels, regardless of ADAP
          ka = 1
          !--(DMK-CARRIO-OLD)-----------------------------------------------------
          !        ka = lpw(i,j)
          !--(DMK-CARRIO-END)-----------------------------------------------------


          if (oneMicControl%jnmb(1) > 0) call ae1kmic(ka,m1,oneMicControl%rx(:,1),oneMicroFields%rcp(:,i,j))
          if (oneMicControl%jnmb(2) > 0) call ae1kmic(ka,m1,oneMicControl%rx(:,2),oneMicroFields%rrp(:,i,j))
          if (oneMicControl%jnmb(3) > 0) call ae1kmic(ka,m1,oneMicControl%rx(:,3),oneMicroFields%rpp(:,i,j))
          if (oneMicControl%jnmb(4) > 0) call ae1kmic(ka,m1,oneMicControl%rx(:,4),oneMicroFields%rsp(:,i,j))
          if (oneMicControl%jnmb(5) > 0) call ae1kmic(ka,m1,oneMicControl%rx(:,5),oneMicroFields%rap(:,i,j))
          if (oneMicControl%jnmb(6) > 0) call ae1kmic(ka,m1,oneMicControl%rx(:,6),oneMicroFields%rgp(:,i,j))
          if (oneMicControl%jnmb(7) > 0) call ae1kmic(ka,m1,oneMicControl%rx(:,7),oneMicroFields%rhp(:,i,j))

          do lcat = 1,ncat
             do k = ka,m1
                if (oneMicControl%rx(k,lcat) < 1.e-9) oneMicControl%rx(k,lcat) = 0.
             enddo
          enddo

          do k = ka,m1
             rtp(k,i,j) = max(0.,rtp(k,i,j))
             vctr9(k) = 1.001 * (oneMicControl%rx(k,1)+ oneMicControl%rx(k,2) + oneMicControl%rx(k,3)  &
                  + oneMicControl%rx(k,4) + oneMicControl%rx(k,5) + oneMicControl%rx(k,6) + oneMicControl%rx(k,7))
          enddo

          do k = ka,m1
             if (vctr9(k) > rtp(k,i,j)) then
                frac = rtp(k,i,j) / (1.e-9 + vctr9(k))
                do lcat = 1,ncat
                   oneMicControl%rx(k,lcat) = oneMicControl%rx(k,lcat) * frac
                enddo
             endif
          enddo

          if (oneMicControl%jnmb(1) > 0) call ae1kmic(ka,m1,oneMicroFields%rcp(:,i,j),oneMicControl%rx(:,1))
          if (oneMicControl%jnmb(2) > 0) call ae1kmic(ka,m1,oneMicroFields%rrp(:,i,j),oneMicControl%rx(:,2))
          if (oneMicControl%jnmb(3) > 0) call ae1kmic(ka,m1,oneMicroFields%rpp(:,i,j),oneMicControl%rx(:,3))
          if (oneMicControl%jnmb(4) > 0) call ae1kmic(ka,m1,oneMicroFields%rsp(:,i,j),oneMicControl%rx(:,4))
          if (oneMicControl%jnmb(5) > 0) call ae1kmic(ka,m1,oneMicroFields%rap(:,i,j),oneMicControl%rx(:,5))
          if (oneMicControl%jnmb(6) > 0) call ae1kmic(ka,m1,oneMicroFields%rgp(:,i,j),oneMicControl%rx(:,6))
          if (oneMicControl%jnmb(7) > 0) call ae1kmic(ka,m1,oneMicroFields%rhp(:,i,j),oneMicControl%rx(:,7))

          if (oneMicControl%jnmb(1) >= 5)  &
               call ae1mic(ka,m1,oneMicroFields%ccp(:,i,j),oneMicroFields%rcp(:,i,j),oneMicControl%rx(1,1))
          if (oneMicControl%jnmb(2) >= 5)  &
               call ae1mic(ka,m1,oneMicroFields%crp(:,i,j),oneMicroFields%rrp(:,i,j),oneMicControl%rx(1,2))
          if (oneMicControl%jnmb(3) >= 5)  &
               call ae1mic(ka,m1,oneMicroFields%cpp(:,i,j),oneMicroFields%rpp(:,i,j),oneMicControl%rx(1,3))
          if (oneMicControl%jnmb(4) >= 5)  &
               call ae1mic(ka,m1,oneMicroFields%csp(:,i,j),oneMicroFields%rsp(:,i,j),oneMicControl%rx(1,4))
          if (oneMicControl%jnmb(5) >= 5)  &
               call ae1mic(ka,m1,oneMicroFields%cap(:,i,j),oneMicroFields%rap(:,i,j),oneMicControl%rx(1,5))
          if (oneMicControl%jnmb(6) >= 5)  &
               call ae1mic(ka,m1,oneMicroFields%cgp(:,i,j),oneMicroFields%rgp(:,i,j),oneMicControl%rx(1,6))
          if (oneMicControl%jnmb(7) >= 5)  &
               call ae1mic(ka,m1,oneMicroFields%chp(:,i,j),oneMicroFields%rhp(:,i,j),oneMicControl%rx(1,7))

          !
          !  Think about how thp should change here - should it be due to a change in
          !     rtp or to a change in the condensate?
          !
          !               vctr10(k) = rrp(k,i,j +rpp(k,i,j)  + rsp(k,i,j) + rap(k,i,j)
          !     +                   + rgp(k,i,j) + rhp(k,i,j)
          !               thp(k,i,j) = thp(k,i,j)
          !     +                    * (1. - aklv * (vctr8(k) - rtp(k,i,j))
          !c or +                    * (1. - aklv * (vctr10(k) - vctr9(k,i,j))
          !     +                    /(max(temp, 253.)))

       enddo
    enddo
    return
  end subroutine adj1

  !---------------------------------------------------------------------------

  subroutine ae1mic(ka,m1,c3,r3,r1)
    ! Arguments
    integer, intent(in)                :: m1,ka
    real, dimension(m1), intent(inout) :: c3
    real, dimension(m1), intent(in)    :: r3, r1

    ! Local Variables
    integer :: k

    do k = ka,m1
       c3(k) = c3(k) * r1(k) / (1.e-9 + r3(k))
       if (c3(k) < 0.) c3(k) = 0.
    enddo

    return
  end subroutine ae1mic

  !---------------------------------------------------------------------------

  subroutine ae1kmic(ka,kb,cr3,cr1)
    ! Arguments
    integer, intent(in) :: ka,kb
    real, dimension(kb), intent(out) :: cr3
    real, dimension(kb), intent(in)  :: cr1

    ! Local Variables
    integer :: k

    do k = ka,kb
       cr3(k) = cr1(k)
    enddo

    return
  end subroutine ae1kmic
end module ModMicrophysicsMisc
