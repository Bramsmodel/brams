!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModMicColl

  use ModMicControl, only: &
       MicControl

  implicit none

  include "MicConstants.h"

  private

  public :: auto_accret
  public :: effxy
  public :: getict
  public :: cols
  public :: col3344
  public :: col3443
  public :: col1
  public :: col2
  public :: col3
  public :: colxfers

contains



  subroutine getict(k1,k2,lcat,i,j,mynum, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: k1,k2,lcat,i,j,k,mynum
    real :: rict,rictmm

    do k = k1,k2

       if (oneMicControl%rx(k,lcat) .ge. 1.e-9) then

          rict = oneMicControl%dict(lcat) * (log(oneMicControl%emb(k,lcat)) - &
               oneMicControl%emb0log(lcat)) + 1.
          rictmm = max(oneMicControl%rictmin,min(oneMicControl%rictmax,rict))
          oneMicControl%ict1(k,lcat) = int(rictmm)
          oneMicControl%ict2(k,lcat) = oneMicControl%ict1(k,lcat) + 1
          oneMicControl%wct2(k,lcat) = rictmm - float(oneMicControl%ict1(k,lcat))
          oneMicControl%wct1(k,lcat) = 1.0 - oneMicControl%wct2(k,lcat)

       endif

    enddo
    return
  end subroutine getict

  !******************************************************************************

  subroutine auto_accret(m1,k1,k2,dn0,dtlt,i,j, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,i,j,k,k1,k2,id1cc,id1cr,id1crn,ir2cr,id2cr,ir2rr,id2rr
    real :: dtlt3,dtlt6,dmb1cgs,dmb2cgs,r2cgs,en1cgs,ad1,ar2,d2minx,ad2  &
         ,bd1,br2,bd2,d2e,bd1cc,bd1cr,br2cr,bd2cr,br2rr,bd2rr,wd1cr,wr2dr  &
         ,wr2rr,wd2rr,tm1cc,tn1cc,tn2cc,tm1cr,tn1cr,tn2rr,en1cgs_2  &
         ,um1cc,un1cc,un2cc,um1cr,un1cr,un2rr,um2,un1,dtlt,cfmasi1,cfmasi2  &
         ,pwmasi1,pwmasi2,wr2cr

    real, dimension(m1) :: dn0

    dtlt3 = 1.e3 * dtlt
    dtlt6 = 1.e6 * dtlt
    cfmasi1 = 1. / oneMicControl%cfmas(1)
    cfmasi2 = 1. / oneMicControl%cfmas(2)
    pwmasi1 = 1. / oneMicControl%pwmas(1)
    pwmasi2 = 1. / oneMicControl%pwmas(2)

    do k = k1,k2
       if(oneMicControl%rx(k,1) .ge. 1.e-9) then

          ! This subroutine works in cgs units, so convert inputs from mks

          dmb1cgs = 100. * (oneMicControl%emb(k,1) * cfmasi1) ** pwmasi1
          dmb2cgs = 100. * (oneMicControl%emb(k,2) * cfmasi2) ** pwmasi2
          r2cgs = 1.e-3 * oneMicControl%rx(k,2) * dn0(k)
          en1cgs = 1.e-6 * oneMicControl%cx(k,1) * dn0(k)

          ad1 = max(oneMicControl%d1min,min(oneMicControl%d1max,dmb1cgs))
          ar2 = max(oneMicControl%r2min,min(oneMicControl%r2max,r2cgs))
          d2minx = max(oneMicControl%d2min,(r2cgs / (.1 * .5236)) ** pwmasi2)
          ad2 = max(d2minx,min(oneMicControl%d2max,dmb2cgs))

          bd1 = alog10(ad1/oneMicControl%d1min)
          br2 = alog10(ar2/oneMicControl%r2min)
          bd2 = alog10(ad2/d2minx)
          d2e =  alog10(oneMicControl%d2max / d2minx)

          bd1cc = float(nd1cc-1) * (ad1 - oneMicControl%d1min) / &
               (oneMicControl%d1max - oneMicControl%d1min) + 1.
          bd1cr = bd1 / oneMicControl%d1ecr + 1.
          br2cr = br2 / oneMicControl%r2ecr + 1.
          bd2cr = bd2 / d2e * float(nd2cr-1) + 1.
          br2rr = br2 / oneMicControl%r2err + 1.
          bd2rr = bd2 / d2e * float(nd2rr-1) + 1.

          !         id1cc  =  int(bd1cc)
          id1cc  =  nint(bd1cc)   
          id1cc  =  max(1,min(id1cc,nd1cc)) !srf bound exception

          id1cr  =  int(bd1cr)
          id1cr  =  max(1,min(id1cr,nd1cr-1)) !srf bound exception

          id1crn = nint(bd1cr)
          id1crn =  max(1,min(id1crn,nd1cr)) !srf bound exception

          ir2cr  =  int(br2cr)!r1tabcr
          ir2cr  =  max(1,min(ir2cr,nr2cr-1)) !srf bound exception

          id2cr  = nint(bd2cr)
          id2cr  =  max(1,min(id2cr,nd2cr)) !srf bound exception

          ir2rr  =  int(br2rr)!c2tabrr
          ir2rr  =  max(1,min(ir2rr,nr2rr-1)) !srf bound exception

          id2rr  =  int(bd2rr)
          id2rr  =  max(1,min(id2rr,nd2rr-1)) !srf bound exception

          wd1cr = bd1cr - float(id1cr)
          wr2cr = br2cr - float(ir2cr)
          wr2rr = br2rr - float(ir2rr)
          wd2rr = bd2rr - float(id2rr)

          !----------- from CSU version ----
          !To fix array bounds exceptions
          !   if(id1cc.gt.ndcc-1) then
          !     id1cc=ndcc-1
          !     wd1cc=1.-wd1cc
          !   endif
          !   if(id2dd.gt.ndcc-1) then
          !     id2dd=ndcc-1
          !     wd2dd=1.-wd2dd
          !   endif
          !   if(id1cr.gt.ndccr-1) then
          !     id1cr=ndccr-1
          !     wd1cr=1.-wd1cr
          !   endif
          !   if(ir3cr.gt.nrrcr-1) then
          !     ir3cr=nrrcr-1
          !     wr3cr=1.-wr3cr
          !   endif
          !   if(id1cd.gt.ndcd-1)  then
          !     id1cd=ndcd-1
          !     wd1cd=1.-wd1cd
          !   endif
          !   if(id2cd.gt.ndcd-1)  then
          !     id2cd=ndcd-1
          !     wd2cd=1.-wd2cd
          !   endif
          !   if(id2cr.gt.ndccr-1) then
          !     id2cr=ndccr-1
          !     wd2cr=1.-wd2cr
          !   endif
          !----------- from CSU version ----


          tm1cc = oneMicControl%r1tabcc(id1cc)

          tn1cc = oneMicControl%c1tabcc(id1cc)

          tn2cc = oneMicControl%c2tabcc(id1cc)

          !--(DMK-CARRIO-INI)----------------------------------------------------------------
          if(abs(wr2cr) < 1.e-5 .or. ir2cr > nr2cr) then

             tm1cr = (1.-wd1cr) * ((1.-wr2cr) * oneMicControl%r1tabcr(id1cr  ,ir2cr  ,id2cr)) & 
                  +     wd1cr  * ((1.-wr2cr) * oneMicControl%r1tabcr(id1cr+1,ir2cr  ,id2cr))

             tn1cr = (1.-wr2cr) * oneMicControl%c1tabcr(id1crn,ir2cr  ,id2cr)

          else
             !--(DMK-CARRIO-FIM)----------------------------------------------------------------

             tm1cr = (1.-wd1cr) * ((1.-wr2cr) * oneMicControl%r1tabcr(id1cr  ,ir2cr  ,id2cr)  &
                  +                   wr2cr  * oneMicControl%r1tabcr(id1cr  ,ir2cr+1,id2cr))  &
                  +     wd1cr  * ((1.-wr2cr) * oneMicControl%r1tabcr(id1cr+1,ir2cr  ,id2cr)  &
                  +                   wr2cr  * oneMicControl%r1tabcr(id1cr+1,ir2cr+1,id2cr))

             tn1cr =               (1.-wr2cr) * oneMicControl%c1tabcr(id1crn,ir2cr  ,id2cr)  &
                  +                   wr2cr  * oneMicControl%c1tabcr(id1crn,ir2cr+1,id2cr)
             !--(DMK-CARRIO-INI)----------------------------------------------------------------
          endif

          if(abs(wr2rr) < 1.e-5 .or. ir2rr > nr2rr) then

             tn2rr = (1.-wd2rr) * ((1.-wr2rr) * oneMicControl%c2tabrr(ir2rr  ,id2rr  ))  &
                  +     wd2rr  * ((1.-wr2rr) * oneMicControl%c2tabrr(ir2rr  ,id2rr+1))
          else	
             !--(DMK-CARRIO-FIM)----------------------------------------------------------------

             tn2rr = (1.-wd2rr) * ((1.-wr2rr) * oneMicControl%c2tabrr(ir2rr  ,id2rr  )  &
                  +     wr2rr  * oneMicControl%c2tabrr(ir2rr+1,id2rr  ))  &
                  +     wd2rr  * ((1.-wr2rr) * oneMicControl%c2tabrr(ir2rr  ,id2rr+1)  &
                  +     wr2rr  * oneMicControl%c2tabrr(ir2rr+1,id2rr+1))

             !--(DMK-CARRIO-INI)----------------------------------------------------------------
          endif
          !--(DMK-CARRIO-FIM)----------------------------------------------------------------

          en1cgs_2 = en1cgs ** 2

          um1cc = tm1cc * en1cgs_2 * dtlt3
          un1cc = tn1cc * en1cgs_2 * dtlt6
          un2cc = tn2cc * en1cgs_2 * dtlt6
          um1cr = 10. ** tm1cr * en1cgs * dtlt3
          un1cr = 10. ** tn1cr * en1cgs * dtlt6
          un2rr = 10. ** tn2rr * dtlt6

          ! The above values are amounts in kg/m^3 or #/m^3 converted in the
          ! present timestep, but must still be corrected for the effect of
          ! density on fall velocity.  Thus, they must be multiplied by
          ! (dn0i ** .5) which fall velocity is proportional to.  Also, since
          ! rxfer and enxfer are in units of kg/kg and #/kg, respectively, the
          ! above transfer amounts must also be multiplied by dn0i.  Together,
          ! these factors make (dn0i ** 1.5).

          um2 = min(oneMicControl%rx(k,1),(um1cc + um1cr) * oneMicControl%dn0i(k))
          un1 = min(oneMicControl%cx(k,1)*dn0(k),(un1cc + un1cr))

          oneMicControl%rxfer(k,1,2)  =  oneMicControl%rxfer(k,1,2) + um2
          oneMicControl%qrxfer(k,1,2) = oneMicControl%qrxfer(k,1,2) + um2 * oneMicControl%qx(k,1)
          oneMicControl%enxfer(k,1,1) = oneMicControl%enxfer(k,1,1) + un1 - un2cc
          oneMicControl%enxfer(k,1,2) = oneMicControl%enxfer(k,1,2) + un2cc

          ! no collis breakup yet - do not use next line but use col(2,2) in 3d micro

          !cc         enxfer(k,2,2) = enxfer(k,2,2) + un2rr

          ! aerosol loss here?

       endif
    enddo
    return
  end subroutine auto_accret

  !******************************************************************************

  subroutine effxy(m1,k1,k2,i,j, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,i,j,k,ncall7
    integer, dimension(10) :: k1,k2
    real :: dmr
    data ncall7/0/
    save

    !     1 = rp,rs,ra,rg,rh

    if (ncall7 .eq. 0 .and. oneMicControl%jnmb(2) .ge. 1 .and. oneMicControl%jnmb(3) .ge. 1) then
       ncall7 = 7
       do k = 2,m1-1
          oneMicControl%eff(k,1) = 1.0
       enddo
    endif

    !     2 = cs,ca

    if (oneMicControl%jnmb(2) .ge. 1 .or. oneMicControl%jnmb(3) .ge. 1) then
       do k = k1(1),k2(1)

          ! Rough fit from Pruppacher and Klett Fig. 14-14 p. 496:
          !  close to curve for 404 microns.  Replace with auto_accret eventually.

          if (oneMicControl%emb(k,1) .gt. 9.e-13) then
             oneMicControl%eff(k,2) = min(1.,30. * (oneMicControl%emb(k,1) - 9.e-13) ** .15)
          else
             oneMicControl%eff(k,2) = 0.
          endif
       enddo
    endif

    !     3 = rr

    if (oneMicControl%jnmb(2) .ge. 1) then
       do k = k1(2),k2(2)

          if (oneMicControl%rx(k,2) .ge. 1.e-9) then

             ! rain breakup (old)

             !            dmr = dn(k,2) * gnu2
             !            if (dmr .lt. .0006) then
             !               eff(k,3) = 1.0
             !            elseif (dmr .gt. .001446) then
             !               eff(k,3) = -5.0
             !            else
             !               eff(k,3) = exp(2300. * (dmr - .0006))
             !            endif

             ! rain breakup (new - temporary; eventually combine with autoconv/accret

             if (oneMicControl%emb(k,2) .lt. .113e-6) then
                oneMicControl%eff(k,3) = 1.0
             elseif (oneMicControl%emb(k,2) .gt. .158e-5) then
                oneMicControl%eff(k,3) = -5.0
             else
                oneMicControl%eff(k,3) = 2. - exp(.1326e7 * (oneMicControl%emb(k,2) - .113e-6))
             endif

          endif

       enddo

    endif

    !     4 = pp,ps,pa

    if (oneMicControl%jnmb(5) .ge. 1) then
       do k = k1(3),k2(3)
          if (abs(oneMicControl%tx(k,3)+14.) .le. 2.) then
             oneMicControl%eff(k,4) = 1.4
          else
             oneMicControl%eff(k,4) = min(0.2,10. ** (0.035 * oneMicControl%tx(k,3) - 0.7))
          endif

       enddo

       !     5 = ss,sa

       do k = k1(4),k2(4)
          if (abs(oneMicControl%tx(k,4)+14.) .le. 2.) then
             oneMicControl%eff(k,5) = 1.4
          else
             oneMicControl%eff(k,5) = min(0.2,10. ** (0.035 * oneMicControl%tx(k,4) - 0.7))
          endif
       enddo

       !     6 = aa

       do k = k1(5),k2(5)

          if (oneMicControl%rx(k,5) .ge. 1.e-9) then

             if (abs(oneMicControl%tx(k,5)+14.) .le. 2.) then
                oneMicControl%eff(k,6) = 1.4
             elseif (oneMicControl%tx(k,5) .ge. -1.) then
                oneMicControl%eff(k,6) = 1.
             else
                oneMicControl%eff(k,6) = min(0.2,10. ** (0.035 * oneMicControl%tx(k,5) - 0.7))
             endif

          endif

       enddo
    endif

    !     7 = pg,sg,ag,gg,gh

    if (oneMicControl%jnmb(6) .ge. 1) then
       do k = k1(6),k2(6)
          if (oneMicControl%qr(k,6) .gt. 0.) then
             oneMicControl%eff(k,7) = 1.0
          else
             oneMicControl%eff(k,7) = min(0.2,10. ** (0.035 * oneMicControl%tx(k,6) - 0.7))
          endif
       enddo
    endif

    !     8 = ph,sh,ah,gh

    if (oneMicControl%jnmb(7) .ge. 1) then
       do k = k1(7),k2(7)

          if (oneMicControl%rx(k,7) .ge. 1.e-9) then

             if (oneMicControl%qr(k,7) .gt. 0.) then
                oneMicControl%eff(k,8) = 1.0
             else
                oneMicControl%eff(k,8) = min(0.2,10. ** (0.035 * oneMicControl%tx(k,7) - 0.7))
             endif

          endif

       enddo
    endif

    !     9 = cg,ch

    if (oneMicControl%jnmb(2) .ge. 1 .or. oneMicControl%jnmb(3) .ge. 1) then
       do k = k1(1),k2(1)


          ! Rough fit from Pruppacher and Klett Fig. 14-11 p. 485:
          !  close to curves for 142 and 305 microns.  Replace with auto_accret eventually.

          if (oneMicControl%emb(k,1) .gt. 3.4e-14) then
             oneMicControl%eff(k,9) = min(1.,1426. * (oneMicControl%emb(k,1) - 3.4e-14) ** .28)
          else
             oneMicControl%eff(k,9) = 0.
          endif
       enddo
    endif

    !     10 = hh (trial)

    if (oneMicControl%jnmb(7) .ge. 1) then
       do k = k1(7),k2(7)
          oneMicControl%eff(k,10) = max(0.,.1 + .005 * oneMicControl%tx(k,7))
       enddo
    endif

    return
  end subroutine effxy

  !******************************************************************************

  !c      SUBROUTINE EFFAB(m1,MX,MY,EFF,DIAX,DIAY,TMPX,TMPY,TDEW)
  !c      DIMENSION EFF(*),PT(*),DIAX(*),DIAY(*),tmpx(*),tmpy(*),tdew(*)
  !
  !c      IF (MX.EQ.6.AND.MY.EQ.6) THEN
  !c        DO K=2,M1
  !c          EFF(K)=MIN(.2,10.**(.035*(MAX(TMPX(K),TMPY(K))
  !c     +          -273.16)-.7))
  !c        ENDDO
  !c      ELSE
  !c        DO K=2,M1
  !c          TMP=MAX(TMPX(K),TMPY(K))
  !c          EFF(K)=MAX(MIN(1.,(TMP-273.06)*1.E10),
  !c     +           MIN(10.**(.035*TMP-273.16)-.7),.2))
  !c         ENDDO
  !c         IF (MX.LE.5.AND.MY.LE.5) THEN
  !c           DO K=2,M1
  !c            EFF(K)=MAX(EFF(K),MIN(1.4,
  !c     +        (MAX(0.,(2.-ABS(TMPX(k)-259.16)))*(TDEW(k)-TMPX(k))
  !c     +        +MAX(0.,(2.-ABS(TMPY(k)-259.16)))*(TDEW(k)-TMPY(k)))
  !c     +        *1.E10))
  !c            ENDDO
  !c         ENDIF
  !c      ENDIF
  !c      RETURN
  !c      END

  !*********************************************************************************

  subroutine cols(m1,mx,mc1,k1,k2,i,j, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: ipc,m1,mx,mc1,k1,k2,i,j,k
    real :: colnum,tabval

    do k = k1,k2
       if(oneMicControl%rx(k,mx) .ge. 1.e-9) then
          ipc = oneMicControl%ipairc(oneMicControl%jhcat(k,mx),oneMicControl%jhcat(k,mx))

          tabval  &
               = oneMicControl%wct1(k,mx) ** 2 * &
               oneMicControl%coltabc(oneMicControl%ict1(k,mx),oneMicControl%ict1(k,mx),ipc)  &
               + 2. * oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,mx) * &
               oneMicControl%coltabc(oneMicControl%ict1(k,mx),oneMicControl%ict2(k,mx),ipc)  &
               + oneMicControl%wct2(k,mx) ** 2 * &
               oneMicControl%coltabc(oneMicControl%ict2(k,mx),oneMicControl%ict2(k,mx),ipc)

          colnum = oneMicControl%colfacc(k) * &
               oneMicControl%eff(k,mc1) * &
               oneMicControl%cx(k,mx) ** 2 * 10. ** (-tabval)
          oneMicControl%enxfer(k,mx,mx) = oneMicControl%enxfer(k,mx,mx) + &
               min(0.5 * oneMicControl%cx(k,mx),colnum)
       endif
    enddo
    return
  end subroutine cols

  !******************************************************************************

  subroutine col3344(m1,mx,mz,mc1,k1,k2,i,j, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,mx,mz,mc1,k1,k2,i,j,k,ip,ipc
    real :: c1,tabvalx,colamt,tabvaln,colnum

    do k = k1,k2
       if(oneMicControl%rx(k,mx) .ge. 1.e-9) then
          ip = oneMicControl%ipairr(oneMicControl%jhcat(k,mx),oneMicControl%jhcat(k,mx))
          ipc = oneMicControl%ipairc(oneMicControl%jhcat(k,mx),oneMicControl%jhcat(k,mx))
          c1 = oneMicControl%eff(k,mc1) * oneMicControl%cx(k,mx) ** 2

          tabvalx  = oneMicControl%wct1(k,mx) ** 2 &
               * oneMicControl%coltabr(oneMicControl%ict1(k,mx),oneMicControl%ict1(k,mx),ip)  &
               + 2. * oneMicControl%wct1(k,mx) &
               * oneMicControl%wct2(k,mx) * &
               oneMicControl%coltabr(oneMicControl%ict1(k,mx),oneMicControl%ict2(k,mx),ip)  &
               + oneMicControl%wct2(k,mx) ** 2 &
               * oneMicControl%coltabr(oneMicControl%ict2(k,mx),oneMicControl%ict2(k,mx),ip)

          colamt = min(oneMicControl%rx(k,mx),oneMicControl%colfacr2(k) * c1 * 10. ** (-tabvalx))
          oneMicControl%rxfer(k,mx,mz) = oneMicControl%rxfer(k,mx,mz) + colamt
          oneMicControl%qrxfer(k,mx,mz) = oneMicControl%qrxfer(k,mx,mz) + colamt * oneMicControl%qx(k,mx)

          if (oneMicControl%jnmb(mz) >= 5) then

             tabvaln  = oneMicControl%wct1(k,mx) ** 2 &
                  * oneMicControl%coltabc(oneMicControl%ict1(k,mx),oneMicControl%ict1(k,mx),ipc)  &
                  + 2. * oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,mx) * &
                  oneMicControl%coltabc(oneMicControl%ict1(k,mx),oneMicControl%ict2(k,mx),ipc)  &
                  + oneMicControl%wct2(k,mx) ** 2 &
                  * oneMicControl%coltabc(oneMicControl%ict2(k,mx),oneMicControl%ict2(k,mx),ipc)

             colnum = min(0.5 * oneMicControl%cx(k,mx),oneMicControl%colfacc2(k) * c1 * 10. ** (-tabvaln))
             oneMicControl%enxfer(k,mx,mz) = oneMicControl%enxfer(k,mx,mz) + colnum
             oneMicControl%enxfer(k,mx,mx) = oneMicControl%enxfer(k,mx,mx) + colnum

          endif
       endif
    enddo
    return
  end subroutine col3344

  !******************************************************************************

  subroutine col3443(m1,mx,my,mz,k1,k2,i,j, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,mx,my,mz,k1,k2,i,j,k,jhcatx,jhcaty,ipxy,ipyx,ipc
    real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum

    do k = k1,k2
       if(oneMicControl%rx(k,mx) .ge. 1.e-9 .and. oneMicControl%rx(k,my) .ge. 1.e-9) then
          jhcatx = oneMicControl%jhcat(k,mx)
          jhcaty = oneMicControl%jhcat(k,my)
          ipxy = oneMicControl%ipairr(jhcatx,jhcaty)
          ipyx = oneMicControl%ipairr(jhcaty,jhcatx)
          ipc  = oneMicControl%ipairc(jhcatx,jhcaty)
          c1 = oneMicControl%eff(k,4) * oneMicControl%cx(k,mx) * oneMicControl%cx(k,my)

          tabvalx  = oneMicControl%wct1(k,mx) * oneMicControl%wct1(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,mx),oneMicControl%ict1(k,my),ipxy)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct1(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,mx),oneMicControl%ict1(k,my),ipxy)  &
               + oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,mx),oneMicControl%ict2(k,my),ipxy)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct2(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,mx),oneMicControl%ict2(k,my),ipxy)

          rcx = min(oneMicControl%rx(k,mx),c1 * oneMicControl%colfacr(k) * 10. ** (-tabvalx))

          tabvaly = oneMicControl%wct1(k,my) * oneMicControl%wct1(k,mx) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,my),oneMicControl%ict1(k,mx),ipyx)  &
               + oneMicControl%wct2(k,my) * oneMicControl%wct1(k,mx) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,my),oneMicControl%ict1(k,mx),ipyx)  &
               + oneMicControl%wct1(k,my) * oneMicControl%wct2(k,mx) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,my),oneMicControl%ict2(k,mx),ipyx)  &
               + oneMicControl%wct2(k,my) * oneMicControl%wct2(k,mx) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,my),oneMicControl%ict2(k,mx),ipyx)
          
          rcy = min(oneMicControl%rx(k,my),c1 * oneMicControl%colfacr(k) * 10. ** (-tabvaly))

          oneMicControl%rxfer(k,mx,mz) = oneMicControl%rxfer(k,mx,mz) + rcx
          oneMicControl%qrxfer(k,mx,mz) = oneMicControl%qrxfer(k,mx,mz) + rcx * oneMicControl%qx(k,mx)

          oneMicControl%rxfer(k,my,mz) = oneMicControl%rxfer(k,my,mz) + rcy
          oneMicControl%qrxfer(k,my,mz) = oneMicControl%qrxfer(k,my,mz) + rcy * oneMicControl%qx(k,my)

          tabvaln = oneMicControl%wct1(k,mx) * oneMicControl%wct1(k,my) * &
               oneMicControl%coltabc (oneMicControl%ict1(k,mx),oneMicControl%ict1(k,my),ipc)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct1(k,my) * &
               oneMicControl%coltabc (oneMicControl%ict2(k,mx),oneMicControl%ict1(k,my),ipc)  &
               + oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,my) * &
               oneMicControl%coltabc (oneMicControl%ict1(k,mx),oneMicControl%ict2(k,my),ipc)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct2(k,my) * &
               oneMicControl%coltabc (oneMicControl%ict2(k,mx),oneMicControl%ict2(k,my),ipc)
          colnum = c1 * oneMicControl%colfacc(k) * 10. ** (-tabvaln)

          if (oneMicControl%cx(k,mx) .gt. oneMicControl%cx(k,my)) then
             oneMicControl%enxfer(k,my,mz) = min(oneMicControl%cx(k,my),colnum)
             oneMicControl%enxfer(k,mx,mx) = min(oneMicControl%cx(k,mx),colnum)
          else
             oneMicControl%enxfer(k,mx,mz) = min(oneMicControl%cx(k,mx),colnum)
             oneMicControl%enxfer(k,my,my) = min(oneMicControl%cx(k,my),colnum)
          endif

          ! also loss for aerosol

       endif
    enddo
    return
  end subroutine col3443

  !******************************************************************************

  subroutine col1(m1,mx,my,mz,mc4,k1,k2,i,j, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,mx,my,mz,mc4,k1,k2,i,j,k,ipxy,ipc
    real :: c1,tabvalx,rcx,tabvaln,colnum

    do k = k1,k2
       if(oneMicControl%rx(k,mx) .ge. 1.e-9 .and. oneMicControl%rx(k,my) .ge. 1.e-9) then
          ipxy = oneMicControl%ipairr(oneMicControl%jhcat(k,mx),oneMicControl%jhcat(k,my))
          ipc  = oneMicControl%ipairc(oneMicControl%jhcat(k,mx),oneMicControl%jhcat(k,my))
          c1 = oneMicControl%eff(k,mc4) * oneMicControl%cx(k,mx) * oneMicControl%cx(k,my)

          tabvalx = oneMicControl%wct1(k,mx) * oneMicControl%wct1(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,mx),oneMicControl%ict1(k,my),ipxy)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct1(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,mx),oneMicControl%ict1(k,my),ipxy)  &
               + oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,mx),oneMicControl%ict2(k,my),ipxy)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct2(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,mx),oneMicControl%ict2(k,my),ipxy)

          rcx = min(oneMicControl%rx(k,mx),c1 * oneMicControl%colfacr(k) * 10. ** (-tabvalx))

          oneMicControl%rxfer(k,mx,mz) = oneMicControl%rxfer(k,mx,mz) + rcx
          oneMicControl%qrxfer(k,mx,mz) = oneMicControl%qrxfer(k,mx,mz) + rcx * oneMicControl%qx(k,mx)

          if (oneMicControl%jnmb(mx) >= 5) then
             tabvaln = oneMicControl%wct1(k,mx) * oneMicControl%wct1(k,my) * &
                  oneMicControl%coltabc (oneMicControl%ict1(k,mx),oneMicControl%ict1(k,my),ipc)  &
                  + oneMicControl%wct2(k,mx) * oneMicControl%wct1(k,my) * &
                  oneMicControl%coltabc (oneMicControl%ict2(k,mx),oneMicControl%ict1(k,my),ipc)  &
                  + oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,my) * &
                  oneMicControl%coltabc (oneMicControl%ict1(k,mx),oneMicControl%ict2(k,my),ipc)  &
                  + oneMicControl%wct2(k,mx) * oneMicControl%wct2(k,my) * &
                  oneMicControl%coltabc (oneMicControl%ict2(k,mx),oneMicControl%ict2(k,my),ipc)

             colnum = c1 * oneMicControl%colfacc(k) * 10. ** (-tabvaln)
             oneMicControl%enxfer(k,mx,mx) = oneMicControl%enxfer(k,mx,mx) + &
                  min(colnum,oneMicControl%cx(k,mx))

             ! also loss for aerosol

          endif

       endif
    enddo
    return
  end subroutine col1

  !******************************************************************************

  subroutine col2(m1,mx,my,mz,mc2,k1,k2,dn0,dtlt,i,j, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,mx,my,mz,mc2,k1,k2,i,j,k,jhcatx,jhcaty,ipxy,ipyx,ipc,it
    real :: c1,c2,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum0,colnum,rcoal  &
         ,qrcx,qrcy,qrcoal,qcoal,fracliq,tcoal,coalliq,coalice,area,cn13,cn24  &
         ,sip,rsip,qrsip,rfinlz,xtoz,dtlt

    real, dimension(m1) :: dn0

    real, dimension(15) ::  alpha,beta
    !            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    data alpha /00.,00.,00., 1., 1., 1., 1.,00.,00.,00.,00., 1., 1., 1., 1./
    data beta  /00.,00.,00.,1.5,1.1,0.0,0.0,00.,00.,00.,00.,1.2,1.1,1.1,1.3/

    do k = k1,k2
       if(oneMicControl%rx(k,mx) .ge. 1.e-9 .and. oneMicControl%rx(k,my) .ge. 1.e-9) then
          jhcatx = oneMicControl%jhcat(k,mx)
          jhcaty = oneMicControl%jhcat(k,my)
          ipxy = oneMicControl%ipairr(jhcatx,jhcaty)
          ipyx = oneMicControl%ipairr(jhcaty,jhcatx)
          ipc  = oneMicControl%ipairc(jhcatx,jhcaty)
          c2 = oneMicControl%cx(k,mx) * oneMicControl%cx(k,my)
          c1 = oneMicControl%eff(k,mc2) * c2

          tabvalx = oneMicControl%wct1(k,mx) * oneMicControl%wct1(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,mx),oneMicControl%ict1(k,my),ipxy)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct1(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,mx),oneMicControl%ict1(k,my),ipxy)  &
               + oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,mx),oneMicControl%ict2(k,my),ipxy)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct2(k,my) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,mx),oneMicControl%ict2(k,my),ipxy)

          rcx = min(oneMicControl%rx(k,mx),c1 * oneMicControl%colfacr(k) * 10. ** (-tabvalx))

          tabvaly = oneMicControl%wct1(k,my) * oneMicControl%wct1(k,mx) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,my),oneMicControl%ict1(k,mx),ipyx)  &
               + oneMicControl%wct2(k,my) * oneMicControl%wct1(k,mx) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,my),oneMicControl%ict1(k,mx),ipyx)  &
               + oneMicControl%wct1(k,my) * oneMicControl%wct2(k,mx) * &
               oneMicControl%coltabr (oneMicControl%ict1(k,my),oneMicControl%ict2(k,mx),ipyx)  &
               + oneMicControl%wct2(k,my) * oneMicControl%wct2(k,mx) * &
               oneMicControl%coltabr (oneMicControl%ict2(k,my),oneMicControl%ict2(k,mx),ipyx)

          rcy = min(oneMicControl%rx(k,my),c1 * oneMicControl%colfacr(k) * 10. ** (-tabvaly))

          if (oneMicControl%jnmb(mx) >= 5 .or. oneMicControl%jnmb(my) >= 5) then

             tabvaln = oneMicControl%wct1(k,mx) * oneMicControl%wct1(k,my) * &
                  oneMicControl%coltabc (oneMicControl%ict1(k,mx),oneMicControl%ict1(k,my),ipc)  &
                  + oneMicControl%wct2(k,mx) * oneMicControl%wct1(k,my) * &
                  oneMicControl%coltabc (oneMicControl%ict2(k,mx),oneMicControl%ict1(k,my),ipc)  &
                  + oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,my) * &
                  oneMicControl%coltabc (oneMicControl%ict1(k,mx),oneMicControl%ict2(k,my),ipc)  &
                  + oneMicControl%wct2(k,mx) * oneMicControl%wct2(k,my) * &
                  oneMicControl%coltabc (oneMicControl%ict2(k,mx),oneMicControl%ict2(k,my),ipc)

             colnum0 = c2 * oneMicControl%colfacc(k) * 10. ** (-tabvaln)
             colnum = colnum0 * oneMicControl%eff(k,mc2)

          else

             ! Dummy initiation of variables colnum because of ifort compiler error detection:
             ! Run-Time Check Failure. The variable 'col2_$COLNUM' is being used without being defined
             ! Routine=col2, line=655, file=mic_coll.f90
             ! called by Routine=mcphys, line=337, file=mic_driv.f90
             colnum = 0. ! Check latter if correct

             ! Dummy initiation of variables colnum0 because of ifort compiler error detection:
             ! Run-Time Check Failure. The variable 'col2_$COLNUM0' is being used without being defined
             ! Routine=col2, line=624, file=mic_coll.f90
             ! called by Routine=mcphys, line=337, file=mic_driv.f90
             colnum0 = 0. ! Check latter if correct

          endif

          rcoal = rcx + rcy
          qrcx = rcx * oneMicControl%qx(k,mx)
          qrcy = rcy * oneMicControl%qx(k,my)
          qrcoal = qrcx + qrcy
          qcoal = qrcoal / (1.e-13 + rcoal)

          call qtc(qcoal,tcoal,fracliq)

          coalliq = rcoal * fracliq
          coalice = rcoal - coalliq

          ! secondary ice production: cn24 is the number fraction of collected cloud
          ! droplets larger than 24 microns and is obtained from an incomplete gamma
          ! function table.  cn13 is the fraction of collected cloud droplets
          ! smaller than 13 microns.  area is cross section area of collecting ice
          ! per m^3 of atmospheric volume.

          if (tcoal .gt. -8. .and. tcoal .lt. -3.) then

             area = oneMicControl%cx(k,my) * dn0(k) * oneMicControl%sipfac(jhcaty) * &
                  oneMicControl%emb(k,my)  &
                  ** (2.*oneMicControl%pwmasi(jhcaty))
             it = nint(oneMicControl%emb(k,mx) / oneMicControl%emb1(1) * 5000.)
             cn13 = colnum * oneMicControl%gamsip13(it) / (area * dtlt)
             cn24 = min(oneMicControl%cx(k,mx)*dn0(k),colnum0) * oneMicControl%gamsip24(it)
             sip = 9.1e-10 * cn24 * cn13 ** .93
             if (tcoal .lt. -5.) then
                sip = 0.33333 * (tcoal + 8.) * sip
             else
                sip = -0.5 * (tcoal + 3.) * sip
             endif

             rsip = sip * oneMicControl%emb0(3) * oneMicControl%dn0i(k)
             qrsip = qcoal * rsip

             rcoal = rcoal - rsip
             qrcoal = qrcoal - qrsip

             oneMicControl%enxfer(k,mx,3) = oneMicControl%enxfer(k,mx,3) + sip
             oneMicControl%rxfer(k,mx,3) = oneMicControl%rxfer(k,mx,3) + rsip
             oneMicControl%qrxfer(k,mx,3) = oneMicControl%qrxfer(k,mx,3) + qrsip

          endif

          ! ALWAYS NEED (ALPHA + BETA) .GE. 1 but in the (rare) case that
          ! fracliq may be a little larger than fracx due to collected
          ! liquid being above 0C, need (ALPHA + BETA) to be at least 1.1
          ! or 1.2, or need ALPHA itself to be at least 1.0.

          rfinlz = min(rcoal,  &
               alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

          xtoz = min(rcx,rfinlz)

          oneMicControl%rxfer(k,mx,mz) = oneMicControl%rxfer(k,mx,mz) + xtoz
          oneMicControl%rxfer(k,mx,my) = oneMicControl%rxfer(k,mx,my) + rcx - xtoz
          if (my .ne. mz) oneMicControl%rxfer(k,my,mz) = oneMicControl%rxfer(k,my,mz)  &
               + rfinlz - xtoz

          oneMicControl%qrxfer(k,mx,mz) = oneMicControl%qrxfer(k,mx,mz) + &
               oneMicControl%qx(k,mx) * xtoz
          oneMicControl%qrxfer(k,mx,my) = oneMicControl%qrxfer(k,mx,my) + &
               oneMicControl%qx(k,mx) * (rcx - xtoz)
          if (my .ne. mz) oneMicControl%qrxfer(k,my,mz) = oneMicControl%qrxfer(k,my,mz)  &
               + oneMicControl%qx(k,my) * (rfinlz - xtoz)

          oneMicControl%enxfer(k,mx,mx) = oneMicControl%enxfer(k,mx,mx) + &
               min(colnum,oneMicControl%cx(k,mx))
          if (my .ne. mz) oneMicControl%enxfer(k,my,mz) = oneMicControl%enxfer(k,my,mz)  &
               + (rfinlz - xtoz) * min(colnum,oneMicControl%cx(k,my)) / (1.e-20 + rcy)

          ! BUT NEED TO CHANGE THE ABOVE FOR 177 COLLECTION BECAUSE X = Y

          ! also include loss of aerosol

       endif
    enddo
    return
  end subroutine col2

  !******************************************************************************

  subroutine col3(m1,mx,my,mz,k1,k2,i,j, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,mx,my,mz,k1,k2,i,j,k,ipxy,ipyx,ipc,jhcaty
    real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum,colnumx,colnumy,coalnum  &
         ,rcoal,qrcx,qrcy,qrcoal,qcoal,fracliq,coalliq,coalice,xtoz  &
         ,rfinlz,tcoal,cfinlz
    real, dimension(15) :: alpha,beta

    !            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    data alpha /00.,00., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1./
    data beta  /00.,00., 2., 2., 2., 1., 0., 2., 2., 2., 2., 2., 2., 2., 2./

    do k = k1,k2
       if(oneMicControl%rx(k,mx) .ge. 1.e-9 .and. oneMicControl%rx(k,my) .ge. 1.e-9) then
          jhcaty = oneMicControl%jhcat(k,my)
          ipxy = oneMicControl%ipairr(oneMicControl%jhcat(k,mx),jhcaty)
          ipyx = oneMicControl%ipairr(jhcaty,oneMicControl%jhcat(k,mx))
          ipc  = oneMicControl%ipairc(oneMicControl%jhcat(k,mx),jhcaty)
          c1 = oneMicControl%eff(k,1) * oneMicControl%cx(k,mx) * oneMicControl%cx(k,my)

          tabvalx  &
               = oneMicControl%wct1(k,mx) * oneMicControl%wct1(k,my) * oneMicControl%coltabr (oneMicControl%ict1(k,mx),oneMicControl%ict1(k,my),ipxy)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct1(k,my) * oneMicControl%coltabr (oneMicControl%ict2(k,mx),oneMicControl%ict1(k,my),ipxy)  &
               + oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,my) * oneMicControl%coltabr (oneMicControl%ict1(k,mx),oneMicControl%ict2(k,my),ipxy)  &
               + oneMicControl%wct2(k,mx) * oneMicControl%wct2(k,my) * oneMicControl%coltabr (oneMicControl%ict2(k,mx),oneMicControl%ict2(k,my),ipxy)

          rcx = min(oneMicControl%rx(k,mx),c1 * oneMicControl%colfacr(k) * 10. ** (-tabvalx))

          tabvaly  &
               = oneMicControl%wct1(k,my) * oneMicControl%wct1(k,mx) * oneMicControl%coltabr (oneMicControl%ict1(k,my),oneMicControl%ict1(k,mx),ipyx)  &
               + oneMicControl%wct2(k,my) * oneMicControl%wct1(k,mx) * oneMicControl%coltabr (oneMicControl%ict2(k,my),oneMicControl%ict1(k,mx),ipyx)  &
               + oneMicControl%wct1(k,my) * oneMicControl%wct2(k,mx) * oneMicControl%coltabr (oneMicControl%ict1(k,my),oneMicControl%ict2(k,mx),ipyx)  &
               + oneMicControl%wct2(k,my) * oneMicControl%wct2(k,mx) * oneMicControl%coltabr (oneMicControl%ict2(k,my),oneMicControl%ict2(k,mx),ipyx)

          rcy = min(oneMicControl%rx(k,my),c1 * oneMicControl%colfacr(k) * 10. ** (-tabvaly))

          if (oneMicControl%jnmb(mx) >= 5) then
             tabvaln  &
                  = oneMicControl%wct1(k,mx) * oneMicControl%wct1(k,my) * oneMicControl%coltabc (oneMicControl%ict1(k,mx),oneMicControl%ict1(k,my),ipc)  &
                  + oneMicControl%wct2(k,mx) * oneMicControl%wct1(k,my) * oneMicControl%coltabc (oneMicControl%ict2(k,mx),oneMicControl%ict1(k,my),ipc)  &
                  + oneMicControl%wct1(k,mx) * oneMicControl%wct2(k,my) * oneMicControl%coltabc (oneMicControl%ict1(k,mx),oneMicControl%ict2(k,my),ipc)  &
                  + oneMicControl%wct2(k,mx) * oneMicControl%wct2(k,my) * oneMicControl%coltabc (oneMicControl%ict2(k,mx),oneMicControl%ict2(k,my),ipc)

             colnum = c1 * oneMicControl%colfacc(k) * 10. ** (-tabvaln)
             colnumx = min(oneMicControl%cx(k,mx),colnum)
             colnumy = min(oneMicControl%cx(k,my),colnum)
             coalnum = min(colnumx,colnumy)
          endif

          rcoal = rcx + rcy
          qrcx = rcx * oneMicControl%qx(k,mx)
          qrcy = rcy * oneMicControl%qx(k,my)
          qrcoal = qrcx + qrcy
          qcoal = qrcoal / (1.e-20 + rcoal)

          call qtc(qcoal,tcoal,fracliq)

          coalliq = rcoal * fracliq
          coalice = rcoal - coalliq

          if (fracliq .ge. .99) then

             oneMicControl%rxfer(k,my,mx) = oneMicControl%rxfer(k,my,mx) + rcy
             oneMicControl%qrxfer(k,my,mx) = oneMicControl%qrxfer(k,my,mx) + qrcy
             if (oneMicControl%jnmb(mx) >= 5)  &
                  oneMicControl%enxfer(k,my,my) = oneMicControl%enxfer(k,my,my) + colnumy
          else

             rfinlz = min(rcoal,  &
                  alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

             xtoz = min(rcx,rfinlz)

             oneMicControl%rxfer(k,mx,mz) = oneMicControl%rxfer(k,mx,mz) + xtoz
             oneMicControl%rxfer(k,mx,my) = oneMicControl%rxfer(k,mx,my) + rcx - xtoz
             if (my .ne. mz) oneMicControl%rxfer(k,my,mz) = oneMicControl%rxfer(k,my,mz)  &
                  + rfinlz - xtoz

             ! NEED TO USE QCOAL TO TRANSFER Q?

             oneMicControl%qrxfer(k,mx,mz) = oneMicControl%qrxfer(k,mx,mz) + oneMicControl%qx(k,mx) * xtoz
             oneMicControl%qrxfer(k,mx,my) = oneMicControl%qrxfer(k,mx,my) + oneMicControl%qx(k,mx) * (rcx - xtoz)
             if (my .ne. mz) oneMicControl%qrxfer(k,my,mz) = oneMicControl%qrxfer(k,my,mz)  &
                  + oneMicControl%qx(k,my) * (rfinlz - xtoz)

             if (oneMicControl%jnmb(mx) >= 5) then
                if (my .eq. mz) then
                   oneMicControl%enxfer(k,mx,mx) = oneMicControl%enxfer(k,mx,mx) + colnumx
                elseif (colnumy .ge. colnumx) then
                   cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)
                   oneMicControl%enxfer(k,mx,mz) = oneMicControl%enxfer(k,mx,mz) + cfinlz
                   oneMicControl%enxfer(k,mx,mx) = oneMicControl%enxfer(k,mx,mx) + colnumx - cfinlz
                   oneMicControl%enxfer(k,my,my) = oneMicControl%enxfer(k,my,my) + colnumy
                else
                   cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)
                   oneMicControl%enxfer(k,my,mz) = oneMicControl%enxfer(k,my,mz) + cfinlz
                   oneMicControl%enxfer(k,mx,mx) = oneMicControl%enxfer(k,mx,mx) + colnumx
                   oneMicControl%enxfer(k,my,my) = oneMicControl%enxfer(k,my,my) + colnumy - cfinlz
                endif
             endif

          endif
       endif
    enddo

    ! also include loss of aerosol

    return
  end subroutine col3

  !******************************************************************************

  subroutine colxfers(m1,k1,k2,i,j,rloss,enloss, oneMicControl)
    type(MicControl), pointer, intent(in) :: oneMicControl
    integer :: m1,i,j,k,lcat,kd1,kd2,jcat
    integer, dimension(10) :: k1,k2
    real, dimension(m1) :: rloss,enloss

    !  All rxfer values are nonnegative.

    !--(DMK-CARRIO-INI)------------------------------------------------------
    !change_MP
    do k = 1,m1 
       oneMicControl%xcoll(k)=0.
    enddo
    !end change_MP
    !--(DMK-CARRIO-FIM)------------------------------------------------------

    do lcat = 1,7
       if (oneMicControl%jnmb(lcat) .ge. 1) then
          kd1 = k1(lcat)
          kd2 = k2(lcat)

          do k = kd1,kd2
             rloss(k) = 0.
             enloss(k) = 0.
          enddo

          do jcat = 1,7
             ! change this to include enxfer of the same categories
             if (oneMicControl%jnmb(jcat) .ge. 1) then
                if (lcat .ne. jcat) then
                   do k = kd1,kd2
                      rloss(k) = rloss(k) + oneMicControl%rxfer(k,lcat,jcat)
                   enddo
                endif
                do k = kd1,kd2
                   enloss(k) = enloss(k) + oneMicControl%enxfer(k,lcat,jcat)
                enddo
             endif
          enddo

          do k = kd1,kd2
             rloss(k) = min(1.,oneMicControl%rx(k,lcat) / max(1.e-20,rloss(k)))
             enloss(k) = min(1.,oneMicControl%cx(k,lcat) / max(1.e-10,enloss(k)))
          enddo

          do jcat = 1,7
             if (oneMicControl%jnmb(jcat) .ge. 1) then
                if (lcat .ne. jcat) then
                   do k = kd1,kd2
                      oneMicControl%rxfer(k,lcat,jcat) = oneMicControl%rxfer(k,lcat,jcat)*rloss(k)
                      oneMicControl%qrxfer(k,lcat,jcat)=oneMicControl%qrxfer(k,lcat,jcat)*rloss(k)
                   enddo
                endif
                do k = kd1,kd2
                   oneMicControl%enxfer(k,lcat,jcat) = oneMicControl%enxfer(k,lcat,jcat)*enloss(k)
                enddo
             endif
          enddo
       endif
    enddo

    !--(DMK-CARRIO-INI)------------------------------------------------------
    !change_MP xcoll is the fraction of clouds transferred to rain
    if (oneMicControl%jnmb(1) .ge. 1 .and. oneMicControl%jnmb(2) .ge. 1) then
       kd1 = k1(1)
       kd2 = k2(1)
       do k= kd1,kd2
          if(oneMicControl%rx(k,1) .gt. 1.e-9)then
             oneMicControl%xcoll(k) = oneMicControl%rxfer(k,1,2)/oneMicControl%rx(k,1)
          endif
       enddo
    endif
    !end change_MP
    !--(DMK-CARRIO-FIM)------------------------------------------------------

    do lcat = 1,7

       if (oneMicControl%jnmb(lcat) .ge. 1) then

          kd1 = k1(lcat)
          kd2 = k2(lcat)

          do jcat = 1,7
             if (oneMicControl%jnmb(jcat) .ge. 1 .and. lcat .ne. jcat) then
                do k = kd1,kd2
                   oneMicControl%rx(k,lcat) = oneMicControl%rx(k,lcat) - oneMicControl%rxfer(k,lcat,jcat)
                   oneMicControl%rx(k,jcat) = oneMicControl%rx(k,jcat) + oneMicControl%rxfer(k,lcat,jcat)
                   oneMicControl%qr(k,lcat) = oneMicControl%qr(k,lcat) - oneMicControl%qrxfer(k,lcat,jcat)
                   oneMicControl%qr(k,jcat) = oneMicControl%qr(k,jcat) + oneMicControl%qrxfer(k,lcat,jcat)
                   oneMicControl%cx(k,lcat) = oneMicControl%cx(k,lcat) - oneMicControl%enxfer(k,lcat,jcat)
                   oneMicControl%cx(k,jcat) = oneMicControl%cx(k,jcat) + oneMicControl%enxfer(k,lcat,jcat)
                enddo
             endif
          enddo

          if (oneMicControl%jnmb(lcat) >= 5) then
             do k = kd1,kd2
                oneMicControl%cx(k,lcat) = oneMicControl%cx(k,lcat) - oneMicControl%enxfer(k,lcat,lcat)
             enddo
          endif

       endif
    enddo
    return
  end subroutine colxfers
end module ModMicColl
