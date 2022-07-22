!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModMicrophysicsDrive

  use ModMicControl, only: &
       MicControl
  
  use ModMicInit, only: &
       micinit
  
  use ModMicNuc, only: &
       cldnuc, &
       icenuc
  
  use ModMicTabs, only: &
       make_autotab, &
       haznuc, &
       tabmelt, &
       tabhab, &
       mksedim_tab, &
       homfrzcl

  use ModMicColl, only: &
       auto_accret, &
       effxy, &
       getict, &
       cols, &
       col3344, &
       col3443, &
       col1, &
       col2, &
       col3, &
       colxfers

  use ModMicrophysicsMisc, only: &
       each_column, &
       enemb, &
       x02, &
       pc03, &
       sedim, &
       c03, &
       each_call, &
       range_check, &
       ae1kmic

  use ModMicVap, only: &
       thrmstr, &
       vapdiff, &
       vapflux, &
       psxfer, &
       newtemp, &
       diffprep

  use ModMicroFields, only: &
       MicroFields
  
  use mem_micro, only:  &
       DeepCopyToMicroFields, &
       DeepCopyFromMicroFields, &
       micro_g, &
       micro_vars        ! INTENT(IN) ! Only a type structure

  use grid_dims, only: &
       maxgrds
  
  use mem_chemic, only: &
       chemic_vars, &        ! INTENT(IN) ! Only a type structure
       chemic_g

  use mem_radiate, only: &
       radiate_vars,     & ! INTENT(IN)
       iswrtyp,          & ! INTENT(IN)
       ilwrtyp,          & ! INTENT(IN)
       radfrq, &              ! INTENT(IN)
       radiate_g          ! INTENT(INOUT)

  use ModBasicFields, only : &
       BasicFields

  use mem_grid, only:   &
       ngrids,          & ! INTENT(IN)
       ngrid,           & ! INTENT(IN)
       zm,              & ! INTENT(IN)
       dzt,             & ! INTENT(IN)
       dtlt,            & ! INTENT(IN)
       jdim,            & ! INTENT(IN)
       maxnzp,          & ! INTENT(IN)
       time,            & ! INTENT(IN)
       zt,              & ! INTENT(IN)
       itime1,          & ! INTENT(IN)
       if_adap,         & ! INTENT(IN)
       grid_g,          & ! INTENT(IN)
       nnzp

  use node_mod, only :  &
       mzp,             & ! INTENT(IN)
       mxp,             & ! INTENT(IN)
       myp,             & ! INTENT(IN)
       ja,              & ! INTENT(IN)
       jz,              & ! INTENT(IN)
       ia,              & ! INTENT(IN)
       iz,              & ! INTENT(IN)
       mynum              ! INTENT(IN)

  use mem_chem1, only : &
       chemistry
  
  use mem_chem1aq, only: &
       chemistry_aq

  implicit none

  include "MicConstants.h"
  
  private
  public :: micro
  
contains





  subroutine micro(oneBasicFields, oneMicControl, oneMicroFields)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicControl), pointer, intent(in) :: oneMicControl
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    
    ! Local Variables:
    integer :: nembfall,maxkfall,ngr,lhcat,i,j
    integer, dimension(10)  :: k1,k2,k3
    integer, save :: ncall = 0
    integer, save, dimension(15)  :: ncall2g = 0
    real :: dtlti ! rshort, cosz, rlongup, rlong, albedt ! Not used
    type pcp_tab_type
       real, pointer, dimension(:,:,:,:) :: pcpfillc,pcpfillr
       real, pointer, dimension(:,:,:)   :: sfcpcp
    end type pcp_tab_type
    type (pcp_tab_type), save :: pcp_tab(maxgrds)

    character(len=*), parameter :: h="**(micro)**"
    
    if (oneMicControl%level .ne. 3) return

    nembfall = 20
    maxkfall = 4

    if(ncall == 0) then
       ncall = 1

       do ngr = 1,ngrids
          allocate (pcp_tab(ngr)%pcpfillc(nnzp(ngr),maxkfall,nembfall,nhcat))
          allocate (pcp_tab(ngr)%pcpfillr(nnzp(ngr),maxkfall,nembfall,nhcat))
          allocate (pcp_tab(ngr)%sfcpcp(maxkfall,nembfall,nhcat))
       enddo

       call micinit(oneMicControl)
       call make_autotab(oneMicControl)
       call haznuc(oneMicControl)
       call tabmelt(oneMicControl)
       call tabhab(oneMicControl)

       do lhcat = 1,nhcat
          oneMicControl%ch3(lhcat) = oneMicControl%pwvt(lhcat) * oneMicControl%pwmasi(lhcat)
          oneMicControl%cdp1(lhcat) = oneMicControl%pwmasi(lhcat) * (1.5 + .5 * oneMicControl%pwvt(lhcat))
          oneMicControl%pwvtmasi(lhcat) = oneMicControl%pwvt(lhcat) * oneMicControl%pwmasi(lhcat)
       enddo
    endif

    if (ncall2g(ngrid) .ne. 5) then
       ncall2g(ngrid) = 5

       call mksedim_tab(mzp,mxp,myp,ngrid,nembfall,maxkfall,zm,dzt  &
            ,pcp_tab(ngrid)%pcpfillc,pcp_tab(ngrid)%pcpfillr  &
            ,pcp_tab(ngrid)%sfcpcp,oneMicControl)

       do lhcat = 1,nhcat
          oneMicControl%ch2(lhcat,ngrid) = float(nembfall-1) &
               / log10(oneMicControl%dispemb1(lhcat,ngrid) / oneMicControl%dispemb0(lhcat,ngrid))
       enddo

       call homfrzcl(dtlt,ngrid,oneMicControl)
    endif

    call each_call(mzp,dtlt,oneMicControl)
    dtlti = 1. / dtlt
    ngr = ngrid

    do j = ja,jz
       do i = ia,iz

          call DeepCopyToMicroFields(oneMicroFields, h)
          call range_check(mzp,k1,k2,k3,i,j,grid_g(ngr)%lpw(i,j),oneMicroFields,oneMicControl)
          call DeepCopyFromMicroFields(oneMicroFields, h)

          call mcphys(mzp,k1,k2,k3,i,j,ngrid,jdim,maxnzp             &
               ,nembfall,maxkfall,mynum,dtlt,dtlti,time,zm,dzt                 &
               ,zt,itime1,radiate_g(ngr)                &
               ,oneBasicFields%thp(:,i,j), oneBasicFields%theta(:1,i,j)   &
               ,oneBasicFields%pp(:,i,j), oneBasicFields%rtp(:,i,j)   &
               ,oneBasicFields%rv(:,i,j), oneBasicFields%wp(:,i,j)   &
               ,oneBasicFields%dn0(:,i,j), oneBasicFields%pi0(:,i,j)   &
               ,grid_g(ngr)%rtgt(i,j), grid_g(ngr)%lpw(i,j)     &
               ,micro_g(ngr)%pcpg(i,j), micro_g(ngr)%qpcpg(i,j)     &
               ,micro_g(ngr)%dpcpg(i,j), pcp_tab(ngr)%pcpfillc &
               ,pcp_tab(ngr)%pcpfillr, pcp_tab(ngr)%sfcpcp    &
               ,grid_g(ngr)%glat(i,j), grid_g(ngr)%topt(i,j),if_adap, oneMicControl)

          call copyback(mzp,k1,k2,k3,grid_g(ngr)%lpw(i,j),i,j,micro_g(ngr), oneMicControl)

          !change_MP for chem
          if(chemistry > 0 .and. chemistry_aq >= 1)then
             call copy2(mzp,k1,k2,k3,grid_g(ngr)%lpw(i,j),i,j,chemic_g(ngr), oneMicControl)
          endif
          !end change_MP

       enddo
    enddo

    return
  end subroutine micro

  !     *****************************************************************

  subroutine mcphys(m1,k1,k2,k3,i,j,ngr,jdim,maxnzp  &
       ,nembfall,maxkfall,mynum,dtlt,dtlti,time,zm,dzt,zt  &
       ,itime1,radiate  &
       ,thp,theta,pp,rtp,rv,wp,dn0,pi0  &
       ,rtgt,lpw_R,pcpg,qpcpg,dpcpg  &
       ,pcpfillc,pcpfillr,sfcpcp  &
       ,glat,topt,if_adap, oneMicControl)
    ! Arguments:
    integer, intent(in) :: m1
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: ngr
    integer, intent(in) :: jdim
    integer, intent(in) :: maxnzp
    integer, intent(in) :: nembfall
    integer, intent(in) :: maxkfall
    integer, intent(in) :: mynum
    integer, intent(inout) :: k1(10)
    integer, intent(inout) :: k2(10)
    integer, intent(inout) :: k3(10)
    real, intent(in) :: dtlt
    real, intent(in) :: dtlti
    real, intent(in) :: time
    real, intent(in) :: zm(m1)
    real, intent(in) :: dzt(m1) ! Not used
    integer, intent(in) :: itime1 ! Not used
    type (radiate_vars), intent(inout)    :: radiate
    real, intent(inout) :: thp(m1)
    real, intent(inout) :: theta(m1)
    real, intent(in) :: pp(m1)
    real, intent(inout) :: rtp(m1)
    real, intent(out) :: rv(m1)
    real, intent(in) :: wp(m1)
    real, intent(in) :: dn0(m1)
    real, intent(in) :: pi0(m1)
    real, intent(in) :: rtgt
    real, intent(in) :: lpw_R
    real, intent(out) :: pcpg
    real, intent(out) :: qpcpg
    real, intent(out) :: dpcpg
    integer, intent(in) :: if_adap
    ! Variables needed for Harrington radiation scheme
    real, intent(in) :: pcpfillc(m1,maxkfall,nembfall,nhcat)
    real, intent(in) :: pcpfillr(m1,maxkfall,nembfall,nhcat)
    real, intent(in) :: sfcpcp(maxkfall,nembfall,nhcat)
    real, intent(in) :: glat
    real, intent(in) :: topt
    real, intent(in) :: zt(m1)
    type(MicControl), pointer, intent(in) :: oneMicControl

    ! Local Variables:
    integer :: lpw
    integer :: k,jflag,lcat,icv,icx,mc1,mc2,mc3,mc4  &
         ,mcat,k1cnuc,k2cnuc,k1pnuc,k2pnuc,lhcat
    integer, dimension(7)   :: mcats,mivap,mix02
    real, dimension(7)      :: dpcp0
    integer, dimension(9,4) :: mcat1
    integer, dimension(7,2) :: mcat2
    integer, dimension(4)   :: mcat33

    data mcats /0,3,0,0,6,7,10/
    data mcat1 /3,3,3,4,4,4,5,5,6  &
         ,5,6,7,5,6,7,6,7,7  &
         ,5,6,7,5,6,7,6,7,7  &
         ,4,7,8,5,7,8,7,8,8/
    data mcat2 /0,0,0,6,6,7,7  &
         ,0,0,0,2,2,9,9/
    data mcat33 /0,0,4,5/
    data mivap /1,3,4,5,2,6,7/
    data mix02 /3,1,4,5,6,7,2/
    data dpcp0 /.001,.001,.010,.010,.010,.003,.001/
    save

    lpw=int(lpw_R)
    call thrmstr(m1,k1,k2,lpw,pp(1),thp(1),theta(1),pi0(1),rtp(1),rv(1),i,j,oneMicControl)

    call each_column(m1,k1,k2,i,j,lpw,rv(1),dn0(1),oneMicControl)

    ! Diagnose hydrometeor mean mass emb, and if necessary, number concentration.

    jflag = 1

    do lcat = 1,7
       if (oneMicControl%jnmb(lcat) .ge. 1) then
          call enemb(m1,k1(lcat),k2(lcat),lcat,jflag,dn0(1),i,j,oneMicControl)
       endif
    enddo

    ! Evaluate radiative heating rates if using Harrington radiation scheme

    !srf- this radiation is not available in BRAMS 5.2 onwards
    !  if (iswrtyp .eq. 3 .or. ilwrtyp .eq. 3) then
    !     if (mod(time + .001,radfrq) .lt. dtlt .or. time .lt. .001) then
    !
    !        call radcalc3(m1,maxnzp,ncat,iswrtyp,ilwrtyp,if_adap,lpw  &
    !     ,glat,rtgt,topt  &
    !             ,radiate%albedt  (i,j) ,radiate%cosz  (i,j)  &
    !             ,radiate%rlongup (i,j) ,radiate%rshort(i,j)  &
    !             ,radiate%rlong   (i,j)  &
    !             ,zm,zt,rv(1),dn0(1),radiate%fthrd(1,i,j),i,j,time,ngr)!
    !
    !     endif
    !  endif

    do lcat = 1,7
       if (oneMicControl%jnmb(lcat) .ge. 1) then
          call diffprep(m1,lcat,k1(lcat),k2(lcat),rv(1),dn0(1),i,j,mynum,oneMicControl)
       endif
    enddo

    call vapdiff(m1,k1(10),k2(10),rv(1),i,j,mynum,oneMicControl)

    do icv = 1,7
       lcat = mivap(icv)

       if (oneMicControl%jnmb(lcat) .ge. 1) then
          call vapflux(m1,lcat,i,j,mynum,k1(lcat),k2(lcat),dn0(1),rv(1),oneMicControl)
       endif
    enddo

    if (oneMicControl%jnmb(4) .ge. 1) then
       call psxfer (m1,min(k1(3),k1(4)),max(k2(3),k2(4)),dn0(1),i,j,oneMicControl)
    endif

    jflag = 2
    do lcat = 1,7
       if (oneMicControl%jnmb(lcat) .ge. 1) then
          call enemb(m1,k1(lcat),k2(lcat),lcat,jflag,dn0(1),i,j,oneMicControl)
          call getict(k1(lcat),k2(lcat),lcat,i,j,mynum,oneMicControl)
       endif
    enddo

    call newtemp(m1,k1(10),k2(10),rv(1),theta(1),i,j,oneMicControl)

    if (oneMicControl%jnmb(2) .ge. 1) then
       call auto_accret(m1,k1(1),k2(1),dn0(1),dtlt,i,j, oneMicControl)
    endif

    call effxy(m1,k1,k2,i,j, oneMicControl)

    ! Self collection of rain, aggregates, graupel, hail:  number change only

    do lcat = 2,7

       if (lcat .eq. 3 .or. lcat .eq. 4) go to 29
       mc1 = mcats(lcat)
       if (oneMicControl%jnmb(lcat) >= 5) then
          call cols (m1,lcat,mc1,k1(lcat),k2(lcat),i,j, oneMicControl)
       endif
29     continue
    enddo

    ! Self collection of pristine ice, snow

    do lcat = 3,4
       mc1 = mcat33(lcat)
       if (oneMicControl%jnmb(lcat) .ge. 1 .and. oneMicControl%jnmb(5) .ge. 1) then
          call col3344 (m1,lcat,5,mc1,k1(lcat),k2(lcat),i,j, oneMicControl)
       endif
    enddo

    ! Collection between pristine ice and snow

    if (oneMicControl%jnmb(5) .ge. 1) then
       call col3443 (m1,3,4,5,max(k1(3),k1(4)),min(k2(3),k2(4)),i,j, oneMicControl)
    endif

    ! Ice-ice collisions

    do icx = 1,9
       mc1 = mcat1(icx,1)
       mc2 = mcat1(icx,2)
       mc3 = mcat1(icx,3)
       mc4 = mcat1(icx,4)

       if (oneMicControl%jnmb(mc1) .ge. 1 .and. oneMicControl%jnmb(mc3) .ge. 1) then
          call col1 (m1,mc1,mc2,mc3,mc4,max(k1(mc1),k1(mc2))  &
               ,min(k2(mc1),k2(mc2)),i,j, oneMicControl)
       endif
    enddo

    ! Ice-cloud collisions

    do lcat = 4,7
       mc1 = mcat2(lcat,1)
       mc2 = mcat2(lcat,2)

       if (oneMicControl%jnmb(lcat) .ge. 1 .and. oneMicControl%jnmb(mc1) .ge. 1) then
          call col2 (m1,1,lcat,mc1,mc2,max(k1(1),k1(lcat)),min(k2(1),k2(lcat))  &
               ,dn0(1),dtlt,i,j, oneMicControl)
       endif
    enddo

    ! Ice-rain collisions

    do lcat = 3,7
       if (oneMicControl%jnmb(lcat) .ge. 1 .and. oneMicControl%jnmb(7) .ge. 1) then
          call col3 (m1,2,lcat,7,max(k1(2),k1(lcat)),min(k2(2),k2(lcat)),i,j, oneMicControl)
       endif
    enddo

    call colxfers(m1,k1,k2,i,j,oneMicControl%scrmic1,oneMicControl%scrmic2, oneMicControl)

    do mcat = 1,7
       lcat = mix02(mcat)
       if (oneMicControl%jnmb(lcat) .ge. 1) call x02(m1,k1,k2,lcat,dn0(1),i,j,oneMicControl)
    enddo

    if (oneMicControl%jnmb(1) .ge. 1) then
       call cldnuc(m1,k1cnuc,k2cnuc,lpw,rv,wp,i,j,oneMicControl)
    endif

    k1(1) = min(k1(1),k1cnuc)
    k2(1) = max(k2(1),k2cnuc)
    k3(1) = max(k2(1),k3(1))

    if (oneMicControl%jnmb(1) .ge. 1) then
       call c03(m1,k1(1),k2(1),1,dn0(1),i,j,oneMicControl)
    endif

    if (oneMicControl%jnmb(3) .ge. 1) then
       call icenuc(m1,k1(1),k2(1),k1pnuc,k2pnuc,lpw,ngr,rv,dn0,dtlt,i,j,oneMicControl)
    endif

    k1(3) = min(k1(3),k1pnuc)
    k2(3) = max(k2(3),k2pnuc)
    k3(3) = max(k2(3),k3(3))

    do lcat = 3,1,-2
       if (oneMicControl%jnmb(lcat) .ge. 1) then
          call pc03(m1,k1(lcat),k2(lcat),lcat,dn0(1),i,j,oneMicControl)
       endif
    enddo

    !  Zero out precip arrays.

    pcpg = 0.
    qpcpg = 0.
    dpcpg = 0.

    ! tairc is used here to accumulate changes to thp from sedim

    do k = lpw,m1
       oneMicControl%tairc(k) = 0.
    enddo

    do lhcat = 2,nhcat
       oneMicControl%ch1(lhcat) = dtlt * oneMicControl%cfvt(lhcat) / rtgt
    enddo

    do lcat = 2,7
       if (oneMicControl%jnmb(lcat) .ge. 1) then
          ! New subroutine from RAMS 6.0
          !        call sedim (m1,lcat,ngr,nembfall,maxkfall  &
          !             ,k1(lcat),k2(lcat),lpw,i,j  &
          !             ,rtp(1),thp(1),theta(1),dn0(1),dpcp0(lcat)  &
          !             ,pcpg,qpcpg,dpcpg,dtlti,oneMicControl%scrmic1,oneMicControl%scrmic2,oneMicControl%scrmic3  &
          !             ,pcpfillc,pcpfillr,sfcpcp,dzt,if_adap)
          call sedim (m1,lcat,ngr,nembfall,maxkfall  &
               ,k1(lcat),k2(lcat),lpw,i,j  &
               ,rtp,thp,theta,dn0,dpcp0(lcat)  &
               ,pcpg,qpcpg,dpcpg,dtlti,&
               oneMicControl%scrmic1(1:m1),&
               oneMicControl%scrmic2(1:m1),&
               oneMicControl%scrmic3(1:m1)  &
               ,pcpfillc,pcpfillr,sfcpcp,dzt,if_adap,oneMicControl)
          !sedim            (m1,lcat,ngr,nembfall,maxkfall,
          !             k1,k2,lpw,i,j  &
          !             ,rtp,thp,theta,dn0,alphasfc  &
          !            ,pcpg,qpcpg,dpcpg,dtlti,cnew,rnew,qrnew,
          !              pcpfillc,pcpfillr,sfcpcp,dzt, if_adap)


       endif
    enddo

    do k = lpw,m1
       thp(k) = thp(k) + oneMicControl%tairc(k)
    enddo

    return
  end subroutine mcphys

  !******************************************************************************

  subroutine copyback(m1,k1,k2,k3,lpw_R,i,j,micro, oneMicControl)
    ! Arguments:
    integer, intent(in) :: m1
    type (micro_vars), intent(inout) :: micro
    integer, intent(in) :: k1(10)
    integer, intent(in) :: k2(10)
    integer, intent(in) :: k3(10)
    integer, intent(in) :: i
    integer, intent(in) :: j
    real, intent(in) :: lpw_R
    type(MicControl), pointer, intent(in) :: oneMicControl

    integer :: lpw
    lpw=int(lpw_R)


    if (oneMicControl%jnmb(1) >= 1) then
       call ae1kmic(lpw,k3(1),micro%rcp(:,i,j),oneMicControl%rx(:,1))
       if (oneMicControl%jnmb(1) >= 5) call ae1kmic(lpw,k3(1),micro%ccp(:,i,j),oneMicControl%cx(:,1))
    endif

    if (oneMicControl%jnmb(2) >= 1) then
       call ae1kmic(lpw,k2(10),micro%rrp(:,i,j),oneMicControl%rx(:,2))
       call ae1kmic(lpw,k2(10),micro%q2(:,i,j),oneMicControl%qx(:,2))
       micro%accpr(i,j) = micro%accpr(i,j) + oneMicControl%accpx(2)
       micro%pcprr(i,j) = oneMicControl%pcprx(2)
       if (oneMicControl%jnmb(2) >= 5) call ae1kmic(lpw,k3(1),micro%crp(:,i,j),oneMicControl%cx(:,2))
    endif

    if (oneMicControl%jnmb(3) >= 1) then
       call ae1kmic(lpw,k3(3),micro%rpp(:,i,j),oneMicControl%rx(:,3))
       micro%accpp(i,j) = micro%accpp(i,j) + oneMicControl%accpx(3)
       micro%pcprp(i,j) = oneMicControl%pcprx(3)
       if (oneMicControl%jnmb(3) >= 5) call ae1kmic(lpw,k3(3),micro%cpp(:,i,j),oneMicControl%cx(:,3))
    endif

    if (oneMicControl%jnmb(4) >= 1) then
       call ae1kmic(lpw,k2(10),micro%rsp(:,i,j),oneMicControl%rx(:,4))
       micro%accps(i,j) = micro%accps(i,j) + oneMicControl%accpx(4)
       micro%pcprs(i,j) = oneMicControl%pcprx(4)
       if (oneMicControl%jnmb(4) >= 5) call ae1kmic(lpw,k2(10),micro%csp(:,i,j),oneMicControl%cx(:,4))
    endif

    if (oneMicControl%jnmb(5) >= 1) then
       call ae1kmic(lpw,k2(10),micro%rap(:,i,j),oneMicControl%rx(:,5))
       micro%accpa(i,j) = micro%accpa(i,j) + oneMicControl%accpx(5)
       micro%pcpra(i,j) = oneMicControl%pcprx(5)
       if (oneMicControl%jnmb(5) >= 5) call ae1kmic(lpw,k2(10),micro%cap(:,i,j),oneMicControl%cx(:,5))
    endif

    if (oneMicControl%jnmb(6) >= 1) then
       call ae1kmic(lpw,k2(10),micro%rgp(:,i,j),oneMicControl%rx(:,6))
       call ae1kmic(lpw,k2(10),micro%q6(:,i,j),oneMicControl%qx(:,6))
       micro%accpg(i,j) = micro%accpg(i,j) + oneMicControl%accpx(6)
       micro%pcprg(i,j) = oneMicControl%pcprx(6)
       if (oneMicControl%jnmb(6) >= 5) call ae1kmic(lpw,k2(10),micro%cgp(:,i,j),oneMicControl%cx(:,6))
    endif

    if (oneMicControl%jnmb(7) >= 1) then
       call ae1kmic(lpw,k2(10),micro%rhp(:,i,j),oneMicControl%rx(:,7))
       call ae1kmic(lpw,k2(10),micro%q7(:,i,j),oneMicControl%qx(:,7))
       micro%accph(i,j) = micro%accph(i,j) + oneMicControl%accpx(7)
       micro%pcprh(i,j) = oneMicControl%pcprx(7)
       if (oneMicControl%jnmb(7) >= 5) call ae1kmic(lpw,k2(10),micro%chp(:,i,j),oneMicControl%cx(:,7))
    endif

    return
  end subroutine copyback


  !--(DMK-CCATT-INI)-----------------------------------------------------
  !******************************************************************************
  !change_MP
  !
  subroutine copy2(m1,k1,k2,k3,lpw_R,i,j,chemic, oneMicControl)
    ! Arguments:
    integer, intent(in)              :: m1
    type (chemic_vars), intent(inout) :: chemic 
    integer, dimension(10), intent(in) :: k1, k2, k3
    integer, intent(in) :: i,j
    real, intent(in) :: lpw_R
    type(MicControl), pointer, intent(in) :: oneMicControl
    
    integer :: lpw

    lpw=int(lpw_R)

    if (oneMicControl%jnmb(2) >= 1) then
       call ae1kmic(lpw,k2(10),chemic%coll(:,i,j),oneMicControl%xcoll)
       call ae1kmic(lpw,k2(10),chemic%sedimr(:,i,j),oneMicControl%rsedim)
    endif

    return
  end subroutine copy2

  !end change_MP
  !--(DMK-CCATT-FIM)-----------------------------------------------------

end module ModMicrophysicsDrive
