!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModMicInit

  use rconstants, only: &
       cp !INTENT(IN)

  use mem_grid, only: &
       zt !INTENT(IN)

  use ModMicGamma, only: &
       gammln, &
       gammp, &
       gammq

  use ModMicTabs, only: &
       mkcoltb

  use dump, only: &
       dumpMessage

  use ModMicControl, only: &
       MicControl

  implicit none

  private

  include "files.h"
  include "constants.h"
  include "MicConstants.h"
  
  public :: micro_master
  public :: initqin
  public :: effective_radius
  public :: jnmbinit
  public :: micinit
contains



  subroutine micro_master(oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: lhcat,khcat,lcat,nd1,nd2,nip,ilcat,ilhcat,idum

    real, dimension(7,15) :: dstprms
    real, dimension(15,15) :: jpairr,jpairc
    character(len=80)            :: dataline
    character(len=f_name_length) :: cname
    character(len=8) :: c0
    character(len=*), parameter :: h="**(micro_master)**"

    REAL :: auxGnu,auxEmb0,auxEmb1,auxCfmas,auxPwmas &
         ,auxCfvt,auxPwvt !Auxiliary variables LFR

    data dstprms/ &
         !----------------------------------------------------------------------
         ! shape      cfmas   pwmas      cfvt    pwvt     dmb0      dmb1
         !----------------------------------------------------------------------
         .5,      524.,     3.,    3173.,     2.,   2.e-6,   40.e-6,  & !cloud
         .5,      524.,     3.,     149.,     .5,   .1e-3,    5.e-3,  & !rain
         .179,     110.8,   2.91,  5.769e5,   1.88,  15.e-6,  125.e-6,  & !pris col
         .179,  2.739e-3,   1.74,  188.146,   .933,   .1e-3,   10.e-3,  & !snow col
         .5,      .496,    2.4,    3.084,     .2,   .1e-3,   10.e-3,  & !aggreg
         .5,      157.,     3.,     93.3,     .5,   .1e-3,    5.e-3,  & !graup
         .5,      471.,     3.,     161.,     .5,   .8e-3,   10.e-3,  & !hail
         .0429,     .8854,    2.5,     316.,   1.01,      00,       00,  & !pris hex
         .3183,   .377e-2,     2.,     316.,   1.01,      00,       00,  & !pris den
         .1803,   1.23e-3,    1.8,  5.769e5,   1.88,      00,       00,  & !pris ndl
         .5,     .1001,  2.256,   3.19e4,   1.66,      00,       00,  & !pris ros
         .0429,     .8854,    2.5,    4.836,    .25,      00,       00,  & !snow hex
         .3183,   .377e-2,     2.,    4.836,    .25,      00,       00,  & !snow den
         .1803,   1.23e-3,    1.8,  188.146,   .933,      00,       00,  & !snow ndl
         .5,     .1001,  2.256,  1348.38,  1.241,      00,       00/    !snow ros

    data jpairr/  &
         0,  0,  0,  1,  2,  3,  4,  0,  0,  0,  0,  5,  6,  7,  8,  &
         0,  0,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,  &
         0, 22, 23, 24,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
         25, 26, 27, 28,  0,  0,  0, 29, 30, 31, 32,  0,  0,  0,  0,  &
         33, 34, 35, 36,  0,  0,  0, 37, 38, 39, 40, 41, 42, 43, 44,  &
         45, 46, 47, 48, 49,  0,  0, 50, 51, 52, 53, 54, 55, 56, 57,  &
         58, 59, 60, 61, 62, 63,  0, 64, 65, 66, 67, 68, 69, 70, 71,  &
         0, 72,  0, 73,  0,  0,  0, 74,  0,  0,  0, 75, 76, 77, 78,  &
         0, 79,  0, 80,  0,  0,  0,  0, 81,  0,  0, 82, 83, 84, 85,  &
         0, 86,  0, 87,  0,  0,  0,  0,  0, 88,  0, 89, 90, 91, 92,  &
         0, 93,  0, 94,  0,  0,  0,  0,  0,  0, 95, 96, 97, 98, 99,  &
         100,101,102,  0,  0,  0,  0,103,104,105,106,107,  0,  0,  0,  &
         108,109,110,  0,  0,  0,  0,111,112,113,114,  0,115,  0,  0,  &
         116,117,118,  0,  0,  0,  0,119,120,121,122,  0,  0,123,  0,  &
         124,125,126,  0,  0,  0,  0,127,128,129,130,  0,  0,  0,131/

    data jpairc/  &
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
         0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
         0,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
         4,  5,  6,  7,  0,  0,  0,  8,  9, 10, 11,  0,  0,  0,  0,  &
         12, 13, 14, 15, 16,  0,  0, 17, 18, 19, 20, 21, 22, 23, 24,  &
         25, 26, 27, 28, 29, 30,  0, 31, 32, 33, 34, 35, 36, 37, 38,  &
         39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,  &
         0, 54,  0,  0,  0,  0,  0, 55,  0,  0,  0,  0,  0,  0,  0,  &
         0, 56,  0,  0,  0,  0,  0,  0, 57,  0,  0,  0,  0,  0,  0,  &
         0, 58,  0,  0,  0,  0,  0,  0,  0, 59,  0,  0,  0,  0,  0,  &
         0, 60,  0,  0,  0,  0,  0,  0,  0,  0, 61,  0,  0,  0,  0,  &
         62, 63, 64,  0,  0,  0,  0, 65, 66, 67, 68, 69,  0,  0,  0,  &
         70, 71, 72,  0,  0,  0,  0, 73, 74, 75, 76,  0, 77,  0,  0,  &
         78, 79, 80,  0,  0,  0,  0, 81, 82, 83, 84,  0,  0, 85,  0,  &
         86, 87, 88,  0,  0,  0,  0, 89, 90, 91, 92,  0,  0,  0, 93/

    !  Define several parameters from above data list
    do lhcat=1,nhcat
       oneMicControl%var_shape(lhcat) = dstprms(1,lhcat)
       oneMicControl%cfmas(lhcat) = dstprms(2,lhcat)
       oneMicControl%pwmas(lhcat) = dstprms(3,lhcat)
       oneMicControl%cfvt (lhcat) = dstprms(4,lhcat)
       oneMicControl%pwvt (lhcat) = dstprms(5,lhcat)

       do khcat=1,nhcat
          oneMicControl%ipairc(lhcat,khcat) = jpairc(lhcat,khcat)
          oneMicControl%ipairr(lhcat,khcat) = jpairr(lhcat,khcat)
       enddo
    enddo

    do lcat=1,ncat
       oneMicControl%emb0 (lcat) = oneMicControl%cfmas(lcat) * dstprms(6,lcat) ** oneMicControl%pwmas(lcat)
       oneMicControl%emb1 (lcat) = oneMicControl%cfmas(lcat) * dstprms(7,lcat) ** oneMicControl%pwmas(lcat)
    enddo

    if (oneMicControl%level .ne. 3) return

    if(oneMicControl%mkcoltab.lt.0.or.oneMicControl%mkcoltab.gt.1)then
       write(c0,"(i8)") oneMicControl%mkcoltab
       !call fatal_error(h//'mkcoltab set to '//&
       !      trim(adjustl(c0))//'which is out of bounds')
       iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
            'mkcoltab set to '//&
            trim(adjustl(c0))//'which is out of bounds')
    endif

    cname=oneMicControl%coltabfn(1:len_trim(oneMicControl%coltabfn))

    if(oneMicControl%mkcoltab.eq.1)then

       !**(JP)** not worked yet

       !call fatal_error(h//" mkcoltab == 1 was not worked yet")
       iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
            " mkcoltab == 1 was not worked yet")
       ! Make collection table and write to file

       call mkcoltb(oneMicControl)
       open(91,file=cname(1:len_trim(cname)),form='formatted',status='unknown')
       rewind(91)
       write(91,181)
       do lcat = 1,ncat
          write(91,182)lcat,oneMicControl%gnu(lcat),oneMicControl%emb0(lcat),oneMicControl%emb1(lcat)
       enddo
       write(91,180)
       write(91,183)
       do lhcat = 1,nhcat
          write(91,182)lhcat,oneMicControl%cfmas(lhcat),oneMicControl%pwmas(lhcat)  &
               ,oneMicControl%cfvt(lhcat),oneMicControl%pwvt(lhcat)
       enddo
       write(91,180)
       do nip=1,npairc
          write(91,186)nip
          write(91,184)(nd2,(oneMicControl%coltabc(nd1,nd2,nip)  &
               ,nd1=1,nembc),nd2=1,nembc)
       enddo
       write(91,180)
       do nip=1,npairr
          write(91,187)nip
          write(91,184)(nd2,(oneMicControl%coltabr(nd1,nd2,nip)  &
               ,nd1=1,nembc),nd2=1,nembc)
       enddo

    else

       !**(JP)** badly worked (master should read and broadcast)

       !  Read collection table

       open(91,file=cname(1:len_trim(cname)),form='formatted',status='old')
       read(91,185)dataline
       do ilcat = 1,ncat
          read(91,182)lcat,auxGnu,auxEmb0,auxEmb1
          oneMicControl%gnu (lcat)= auxGnu
          oneMicControl%emb0(lcat)= auxEmb0
          oneMicControl%emb1(lcat)= auxEmb1
       enddo
       read(91,185)dataline
       read(91,185)dataline
       do ilhcat = 1,nhcat
          read(91,182)lhcat,auxCfmas,auxPwmas &
               ,auxCfvt,auxPwvt
          oneMicControl%cfmas(lhcat) = auxCfmas
          oneMicControl%pwmas(lhcat) = auxPwmas
          oneMicControl%cfvt (lhcat) = auxCfvt
          oneMicControl%pwvt (lhcat) = auxPwvt
       enddo
       read(91,185)dataline
       do nip=1,npairc
          read(91,185)dataline
          read(91,184)(idum,(oneMicControl%coltabc(nd1,nd2,nip)  &
               ,nd1=1,nembc),nd2=1,nembc)
       enddo
       read(91,185)dataline
       do nip=1,npairr
          read(91,185)dataline
          read(91,184)(idum,(oneMicControl%coltabr(nd1,nd2,nip)  &
               ,nd1=1,nembc),nd2=1,nembc)
       enddo

    endif

180 format(' ')
181 format(' lcat    gnu        emb0       emb1    ')
182 format(i4,7e11.4)
183 format(' lhcat  cfmas      pwmas       cfvt       pwvt')
184 format(i3,20f6.2)
185 format(a80)
186 format('ipairc',i4)
1186 format(6x,i4)
187 format('ipairr',i4)

    close(91)

    return
  end subroutine micro_master

  !******************************************************************************

  subroutine initqin(n1,n2,n3,q2,q6,q7,pi0,pp,theta,dn0,cccnp,cifnp,oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    ! Arguments
    integer,                   intent(in)  :: n1,n2,n3
    real, dimension(:,:,:), intent(out) :: q2, q6, q7, cifnp

    real, dimension(:,:,:), intent(out)  :: cccnp   ! Carrio 2012

    real, dimension(:,:,:), intent(in)  :: pi0, pp, theta, dn0

    ! Local Variables
    integer :: i,j,k

    ! Initialize Q2, Q6, Q7, CCN, IFN.

    do j = 1,n3
       do i = 1,n2
          do k = 1,n1
             oneMicControl%pitot(k) = pi0(k,i,j) + pp(k,i,j)
             oneMicControl%tair(k) = theta(k,i,j) * oneMicControl%pitot(k) / cp

             if(oneMicControl%irain  .ge. 1) q2(k,i,j) = oneMicControl%tair(k) - 193.16
             if(oneMicControl%igraup .ge. 1) q6(k,i,j) = 0.5 * min(0.,oneMicControl%tair(k) - 273.16)
             if(oneMicControl%ihail  .ge. 1) q7(k,i,j) = 0.5 * min(0.,oneMicControl%tair(k) - 273.16)
             !
             !
             !Carrio 2012 -------------------------------------------------
             !
             !---> icloud= 5: CCN is homogeneously initilaized by CPARM
             if(oneMicControl%icloud .eq. 5) cccnp(k,i,j) = oneMicControl%cparm
             !
             !---> icloud= 6: Same of 5, but decreases above 4km (reasonable!)
             if(oneMicControl%icloud .eq.6)then
                if(k<=2) cccnp(k,i,j)=oneMicControl%cparm
                if(k>2.and.zt(k)<=4000.) cccnp(k,i,j)=max(100.,oneMicControl%cparm * (1.-zt(k)/4000.))
                if(zt(k)>4000.) cccnp(k,i,j) = 100.
             endif
             !
             !---> oneMicControl%icloud= 7: YOU HAVE TO HARD-CODE IT HERE !!!!
             !---> to have a 3-D heterogeneous initialization of CCN
             if(oneMicControl%icloud .eq.7)then
                cccnp(k,i,j)=oneMicControl%cparm
                print*,'if you use ICLOUD = 7, you MUST set up a 3D field'
		print*,'to initialize CCN in mic_init.f90'
		stop
             endif
             !-------------------------------------------------------------
             !

             if (oneMicControl%ipris .eq. 7) cifnp(k,i,j) = 1.e5 * dn0(k,i,j) ** 5.4

          enddo
       enddo
    enddo
    return
  end subroutine initqin

  !******************************************************************************

  subroutine jnmbinit(oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    if (oneMicControl%level /= 3) then

       if (oneMicControl%level <= 1) then
          oneMicControl%jnmb(1) = 0
       else
          oneMicControl%jnmb(1) = 4
       endif

       oneMicControl%jnmb(2) = 0
       oneMicControl%jnmb(3) = 0
       oneMicControl%jnmb(4) = 0
       oneMicControl%jnmb(5) = 0
       oneMicControl%jnmb(6) = 0
       oneMicControl%jnmb(7) = 0

    else

       oneMicControl%jnmb(1) = oneMicControl%icloud
       oneMicControl%jnmb(2) = oneMicControl%irain
       oneMicControl%jnmb(3) = oneMicControl%ipris
       oneMicControl%jnmb(4) = oneMicControl%isnow
       oneMicControl%jnmb(5) = oneMicControl%iaggr
       oneMicControl%jnmb(6) = oneMicControl%igraup
       oneMicControl%jnmb(7) = oneMicControl%ihail

       if (oneMicControl%icloud .eq. 1) oneMicControl%jnmb(1) = 4
       if (oneMicControl%irain  .eq. 1) oneMicControl%jnmb(2) = 2
       if (oneMicControl%ipris  .ge. 1) oneMicControl%jnmb(3) = 5
       if (oneMicControl%isnow  .eq. 1) oneMicControl%jnmb(4) = 2
       if (oneMicControl%iaggr  .eq. 1) oneMicControl%jnmb(5) = 2
       if (oneMicControl%igraup .eq. 1) oneMicControl%jnmb(6) = 2
       if (oneMicControl%ihail  .eq. 1) oneMicControl%jnmb(7) = 2

       if (oneMicControl%irain == 5 .or. oneMicControl%isnow == 5 .or. oneMicControl%iaggr == 5 .or.  &
            oneMicControl%igraup == 5 .or. oneMicControl%ihail == 5) then

          if (oneMicControl%irain  >= 1) oneMicControl%jnmb(2) = 5
          if (oneMicControl%isnow  >= 1) oneMicControl%jnmb(4) = 5
          if (oneMicControl%iaggr  >= 1) oneMicControl%jnmb(5) = 5
          if (oneMicControl%igraup >= 1) oneMicControl%jnmb(6) = 5
          if (oneMicControl%ihail  >= 1) oneMicControl%jnmb(7) = 5

       endif

    endif
    return
  end subroutine jnmbinit

  !******************************************************************************

  subroutine micinit(oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    !Local Variables:
    integer :: lhcat,lcat,ia
    real :: cfmasi,c1,glg,glg1,glg2,glgm,glgc,glgmv,flngi,dpsi,embsip,dnsip
    real :: aux_loop1, aux_loop2

    ! Initialize arrays based on microphysics namelist parameters

    oneMicControl%parm(1) = oneMicControl%cparm
    oneMicControl%parm(2) = oneMicControl%rparm
    !     parm(3) = pparm   [obsolete]
    oneMicControl%parm(4) = oneMicControl%sparm
    oneMicControl%parm(5) = oneMicControl%aparm
    oneMicControl%parm(6) = oneMicControl%gparm
    oneMicControl%parm(7) = oneMicControl%hparm

    if (oneMicControl%icloud .le. 1) oneMicControl%parm(1) = .3e9
    if (oneMicControl%irain  .eq. 1) oneMicControl%parm(2) = .1e-2
    !     if (ipris  .eq. 1) parm(3) = .1e4     [obsolete]
    if (oneMicControl%isnow  .eq. 1) oneMicControl%parm(4) = .1e-2
    if (oneMicControl%iaggr  .eq. 1) oneMicControl%parm(5) = .1e-2
    if (oneMicControl%igraup .eq. 1) oneMicControl%parm(6) = .1e-2
    if (oneMicControl%ihail  .eq. 1) oneMicControl%parm(7) = .3e-2

    oneMicControl%dps = 125.e-6
    oneMicControl%dps2 = oneMicControl%dps ** 2
    oneMicControl%rictmin = 1.0001
    oneMicControl%rictmax = 0.9999 * float(nembc)

    do lhcat = 1,nhcat
       lcat = lhcat + (3 - lhcat) * (lhcat / 8) + lhcat / 12

       oneMicControl%cfden(lhcat) = oneMicControl%cfmas(lhcat) * 6.0 / 3.14159
       oneMicControl%pwden(lhcat) = oneMicControl%pwmas(lhcat) - 3.
       oneMicControl%emb0log(lcat) = log(oneMicControl%emb0(lcat))
       oneMicControl%emb1log(lcat) = log(oneMicControl%emb1(lcat))

       ! Define coefficients [vtfac, frefac1, frefac2] used for terminal velocity
       ! and Reynolds number

       cfmasi = 1. / oneMicControl%cfmas(lhcat)
       oneMicControl%pwmasi(lhcat) = 1. / oneMicControl%pwmas(lhcat)
       oneMicControl%pwen0(lhcat) = 1. / (oneMicControl%pwmas(lhcat) + 1.)
       oneMicControl%pwemb0(lhcat) = oneMicControl%pwmas(lhcat) / (oneMicControl%pwmas(lhcat) + 1.)
       c1 = 1.5 + .5 * oneMicControl%pwvt(lhcat)

       glg = gammln(oneMicControl%gnu(lcat))
       glg1 = gammln(oneMicControl%gnu(lcat) + 1.)
       glg2 = gammln(oneMicControl%gnu(lcat) + 2.)
       glgm = gammln(oneMicControl%gnu(lcat) + oneMicControl%pwmas(lhcat))
       glgc = gammln(oneMicControl%gnu(lcat) + c1)
       glgmv = gammln(oneMicControl%gnu(lcat) + oneMicControl%pwmas(lhcat) + oneMicControl%pwvt(lhcat))

       if (oneMicControl%jnmb(lcat) .eq. 3) then
          oneMicControl%cfemb0(lhcat) = oneMicControl%cfmas(lhcat) * exp(glgm - glg)  &
               ** oneMicControl%pwen0(lhcat) * (1. / oneMicControl%parm(lcat)) ** oneMicControl%pwemb0(lhcat)
          oneMicControl%cfen0(lhcat) = oneMicControl%parm(lcat) * (exp(glg - glgm) / oneMicControl%parm(lcat))  &
               ** oneMicControl%pwen0(lhcat)
       endif

       oneMicControl%dnfac(lhcat) = (cfmasi * exp(glg - glgm)) ** oneMicControl%pwmasi(lhcat)

       oneMicControl%vtfac(lhcat) = oneMicControl%cfvt(lhcat) * exp(glgmv - glgm)  &
            * (cfmasi * exp(glg - glgm)) ** (oneMicControl%pwvt(lhcat) *oneMicControl%pwmasi(lhcat))

       oneMicControl%frefac1(lhcat) = oneMicControl%var_shape(lhcat) * exp(glg1 - glg)  &
            * (cfmasi * exp(glg - glgm)) ** oneMicControl%pwmasi(lhcat)

       oneMicControl%frefac2(lhcat) = oneMicControl%var_shape(lhcat) * 0.229 * sqrt(oneMicControl%cfvt(lcat))  &
            * (cfmasi * exp(glg - glgm)) ** (oneMicControl%pwmasi(lhcat) * c1)  &
            * exp(glgc - glg)

       oneMicControl%sipfac(lhcat) = .785 * exp(glg2 - glg)  &
            * (cfmasi * exp(glg - glgm)) ** (2. * oneMicControl%pwmasi(lhcat))

       oneMicControl%cfmasft(lhcat) = oneMicControl%cfmas(lhcat) * exp(gammln  &
            (oneMicControl%gnu(lcat) + oneMicControl%pwmas(lhcat)) - gammln(oneMicControl%gnu(lcat)))

       oneMicControl%dict(lcat) = float(nembc-1) / (oneMicControl%emb1log(lcat) - oneMicControl%emb0log(lcat))

       oneMicControl%dpsmi(lhcat) = 1. / (oneMicControl%cfmas(lhcat) * oneMicControl%dps ** oneMicControl%pwmas(lhcat))
       if (lhcat .le. 4) oneMicControl%gamm(lhcat) = exp(glg)
       if (lhcat .le. 4) oneMicControl%gamn1(lhcat) = exp(glg1)

       ! gam1   :  the integral of the pristine distribution from dps to infty
       ! gam2   :  the integral of the snow dist. from 0 to dps
       ! gam3   :  values of the exponential exp(-dps/dn)

    enddo

    flngi = 1. / float(ngam)

    ! ALF
    aux_loop1 = oneMicControl%dps * 1.e6
    aux_loop2 = oneMicControl%emb1(1) * flngi

    do ia=1,ngam
       dpsi = aux_loop1 / float(ia)

       oneMicControl%gam(ia,1) = gammq(oneMicControl%gnu(3) + 1., dpsi)
       oneMicControl%gam(ia,2) = gammp(oneMicControl%gnu(4) + 1., dpsi)
       oneMicControl%gam(ia,3) = exp(-dpsi)

       oneMicControl%gaminc(ia,1)=gammq(oneMicControl%gnu(3),dpsi)
       oneMicControl%gaminc(ia,2)=gammp(oneMicControl%gnu(4),dpsi)

       !embsip = emb1(1) * float(ia) * flngi
       embsip = aux_loop2 * float(ia)
       dnsip = oneMicControl%dnfac(1) * embsip ** oneMicControl%pwmasi(1)
       oneMicControl%gamsip13(ia) = gammp(oneMicControl%gnu(1),13.e-6/dnsip)
       oneMicControl%gamsip24(ia) = gammq(oneMicControl%gnu(1),24.e-6/dnsip)
    enddo

    return
  end subroutine micinit
  !******************************************************************************
  subroutine effective_radius(n1,n2,n3,rei,rel)
    ! Arguments
    integer,                   intent(in)  :: n1,n2,n3
    real, dimension(:,:,:), intent(inout) ::rei,rel
    ! Local Variables
    integer :: i,j,k

    do j = 1,n3
       do i = 1,n2
          do k = 1,n1
             rei(k,i,j)=5.0!micrometers for RRTM
	     rel(k,i,j)=2.5!micrometers for RRTM
          enddo
       enddo
    enddo
    return
  end subroutine effective_radius
end module ModMicInit
