!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModNudAnalysis

  use ModNestFillDens, only: &
       fillscr, &
       fillvar
  
  use mem_tend, only: &
       tend

  use ModBasicFields, only: &
       BasicFields

  use mem_grid, only: &
       time, &
       ngrid, &
       jdim, &
       grid_g, &
       nnzp, &
       nnxp, &
       nnyp, &
       ztop, &
       zt, &
       dtlt, &
       nzp, &
       ztn, &
       nxtnest, &
       maxnxp, &
       maxnyp

  use mem_varinit, only: &
       nud_cond, &
       tcond_beg, &
       tcond_end, &
       varinit_g, &
       nud_type, &
       vtime1, &
       vtime2, &
       htime1, &
       htime2, &
       wt_nudge_grid, &
       wt_nudge_uv, &
       wt_nudge_th, &
       wt_nudge_pi, &
       wt_nudge_rt, &
       timeWindowIAU, &
       nudlat, &
       tnudcent, &
       tnudtop, &
       tnudlat, &
       znudtop, &
       ramp, &
       condtime1, &
       condtime2, &
       wt_nudgec_grid, &
       t_nudge_rc


  use node_mod, only: &
       ia, &
       iz, &
       ja, &
       jz, &
       mxp, &
       myp, &
       mzp, &
       ibcon, &
       mynum, &
       master_num, &
       mchnum, &
       i0,            &
       j0,            &
       nmachs,        &
       nodei0,        &
       nodej0,        &
       nodemxp,       &
       nodemyp

  use chem1_list, only: &
       nspecies, &
       fdda, &
       on, &
       spc_alloc

  use mem_chem1, only: &
       chem1_g,        &
       chem_assim

  use mem_scratch, only: &
       vctr1, &
       vctr2, &
       vctr3, &
       scratch, &
       vctr4, &
       vctr14, &
       vctr5, &
       vctr10, &
       vctr11, &
       vctr12, &
       vctr13

  use ModEvaluation, only: &
       comunicateStatistic, &
       evaluate

  use modIau, only: &
       applyIAU

  use dump,only: &
       dumpMessage

  implicit none

  private

  public :: datassim
  public :: varweight
  public :: varweight_chem
  public :: VariableWeight
  public :: VariableWeightChem
  public :: vfintrpf
  public :: VarfIntrp

contains

  subroutine datassim(oneBasicFields)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    
    integer :: il,ir,jl,jr

    !--(DMK-CCATT-INI)-----------------------------------------------------
    integer :: nspc
    !--(DMK-CCATT-END)-----------------------------------------------------

    !     Set bounds for nudging this sub-domain

    il=ia  ;  ir=iz  ;  jl=ja  ;  jr=jz

    !     West and east boundaries.
    if(iand(ibcon,1) /= 0)  il=1
    if(iand(ibcon,2) /= 0)  ir=mxp


    !     South and north boundaries.
    if(jdim == 1) then
       if(iand(ibcon,4) /= 0) jl=1
       if(iand(ibcon,8) /= 0) jr=myp
    endif


    ! Basic boundary and analysis nudging scheme

    call nudge(mzp,mxp,myp,il,ir,jl,jr,varinit_g(ngrid)%varwts   &
         ,varinit_g(ngrid)%varup ,varinit_g(ngrid)%varvp   &
         ,varinit_g(ngrid)%varpp ,varinit_g(ngrid)%vartp   &
         ,varinit_g(ngrid)%varrp  &
         ,varinit_g(ngrid)%varuf ,varinit_g(ngrid)%varvf   &
         ,varinit_g(ngrid)%varpf ,varinit_g(ngrid)%vartf   &
         ,varinit_g(ngrid)%varrf   &
         ,oneBasicFields%up    ,oneBasicFields%vp   &
         ,oneBasicFields%theta ,oneBasicFields%rtp   &
         ,oneBasicFields%pp   &
         ,tend%ut,tend%vt,tend%tht,tend%rtt,tend%pt, oneBasicFields)

    !--(DMK-CCATT-INI)-----------------------------------------------------
    if (chem_assim == 1) then
       do nspc=1,nspecies
          if(spc_alloc(fdda,nspc) == on) then
             call chem_nudge(mzp,mxp,myp,il,ir,jl,jr, nspc  &
                  ,varinit_g(ngrid)%varwts_chem   &
                  ,chem1_g(nspc,ngrid)%sc_pp  & !past
                  ,chem1_g(nspc,ngrid)%sc_pf  & !future
                  ,chem1_g(nspc,ngrid)%sc_p   & !current
 		  ,chem1_g(nspc,ngrid)%sc_t      ) !tendency
          endif
       enddo
    endif
    !--(DMK-CCATT-END)-----------------------------------------------------

    ! Condensate nudging scheme

    if (nud_cond == 1 .and. time >= tcond_beg .and. time <= tcond_end) &
         call nudge_cond(mzp,mxp,myp,il,ir,jl,jr,varinit_g(ngrid)%varwts   &
         ,varinit_g(ngrid)%varrph ,varinit_g(ngrid)%varcph  &
         ,varinit_g(ngrid)%varrfh ,varinit_g(ngrid)%varcfh   &
         ,oneBasicFields%rtp ,tend%rtt)

    return
  end subroutine datassim

  !     ******************************************************************

  subroutine nudge(m1,m2,m3,ia,iz,ja,jz,varwts  &
       ,varup,varvp,varpp,vartp,varrp  &
       ,varuf,varvf,varpf,vartf,varrf  &
       ,up,vp,theta,rtp,pp,ut,vt,tht,rtt,pt,oneBasicFields)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    
    integer :: m1,m2,m3,ia,iz,ja,jz
    real, dimension(m1,m2,m3) :: varup,varvp,vartp,varrp,varpp  &
         ,varuf,varvf,vartf,varrf,varpf  &
         ,varwts,up,vp,theta,rtp,pp  &
         ,ut,vt,tht,rtt,pt

    integer :: i,j,k,iCount,v
    real :: tfact, wt_uv, wt_th, wt_pi, wt_rt
    real :: rmse(m1,5),bias(m1,5),soma(m1,5),somaQ(m1,5),centFact

    !srf - special weights for pressure and/or UV (only for operations) 
    real, dimension(m1,m2,m3) :: varwts_for_operations_only
    if(  wt_nudge_pi < 0. .or. wt_nudge_uv < 0. ) varwts_for_operations_only=0.

    !         Linearly interpolate values in time, then nudge.

    !-- initialize and update (if desired) the nudging weights (lat, top, center)
    call VariableWeight(nnzp(1), nodemxp(mynum,1), nodemyp(mynum,1), nnxp(1),&
         nnyp(1), nodei0 (mynum,1), nodej0 (mynum,1),  &
         grid_g(1)%topt, grid_g(1)%rtgt, varinit_g(1)%varwts ,&
         varwts_for_operations_only)

    if (nud_type == 2 .or. nud_type == 4) then
       tfact=(time-vtime1)/(vtime2-vtime1)
    elseif (nud_type == 1) then
       tfact=(time-htime1)/(htime2-htime1)
    endif

    wt_uv=wt_nudge_grid(ngrid)* wt_nudge_uv
    wt_th=wt_nudge_grid(ngrid)* wt_nudge_th
    wt_pi=wt_nudge_grid(ngrid)* wt_nudge_pi
    wt_rt=wt_nudge_grid(ngrid)* wt_nudge_rt

    if(evaluate==1) then
       soma=0
       somaQ=0
       !iCount=0
    endif
    do j=ja,jz
       do i=ia,iz
          !iCount=iCount+1
          do k=1,m1
             vctr1(k)=varup(k,i,j)+(varuf(k,i,j)-varup(k,i,j))*tfact
             vctr2(k)=varvp(k,i,j)+(varvf(k,i,j)-varvp(k,i,j))*tfact
             vctr3(k)=vartp(k,i,j)+(vartf(k,i,j)-vartp(k,i,j))*tfact
             vctr4(k)=varrp(k,i,j)+(varrf(k,i,j)-varrp(k,i,j))*tfact
             vctr5(k)=varpp(k,i,j)+(varpf(k,i,j)-varpp(k,i,j))*tfact

             if(evaluate==1) then
                !Sum the quadratic error to statistic RMSE
                somaQ(k,1)=somaQ(k,1)+(oneBasicFields%up(k,i,j)-vctr1(k))**2
                somaQ(k,2)=somaQ(k,2)+(oneBasicFields%vp(k,i,j)-vctr2(k))**2
                somaQ(k,3)=somaQ(k,3)+(oneBasicFields%theta(k,i,j)-vctr3(k))**2
                somaQ(k,4)=somaQ(k,4)+(oneBasicFields%rtp(k,i,j)-vctr4(k))**2
                somaQ(k,5)=somaQ(k,5)+(oneBasicFields%pp(k,i,j)-vctr5(k))**2

                !Sum the absolute error to statistic
                soma(k,1)=soma(k,1)+(oneBasicFields%up(k,i,j)-vctr1(k))
                soma(k,2)=soma(k,2)+(oneBasicFields%vp(k,i,j)-vctr2(k))
                soma(k,3)=soma(k,3)+(oneBasicFields%theta(k,i,j)-vctr3(k))
                soma(k,4)=soma(k,4)+(oneBasicFields%rtp(k,i,j)-vctr4(k))
                soma(k,5)=soma(k,5)+(oneBasicFields%pp(k,i,j)-vctr5(k))
             endif

             vctr10(k)=(varwts(k,i,j)+varwts(k,min(m2,i+1),j)   )*.5* wt_uv
             vctr11(k)=(varwts(k,i,j)+varwts(k,i,min(m3,j+jdim)))*.5* wt_uv
             vctr12(k)=varwts(k,i,j)* wt_th
             vctr13(k)=varwts(k,i,j)* wt_pi
             vctr14(k)=varwts(k,i,j)* wt_rt
          enddo

          !-srf special weights for pressure (full domain)
          if(wt_nudge_pi < 0.) then   
             do k=1,m1
                vctr13(k)=varwts_for_operations_only(k,i,j)* abs (wt_pi)
             enddo
          endif
          if(wt_nudge_uv < 0.) then   
             do k=1,m1              
                vctr10(k)=(varwts_for_operations_only(k,i,j)+&
                     varwts_for_operations_only(k,min(m2,i+1),j)   )*.5* abs(wt_uv)
                vctr11(k)=(varwts_for_operations_only(k,i,j)+&
                     varwts_for_operations_only(k,i,min(m3,j+jdim)))*.5* abs(wt_uv)
             enddo
          endif
          !-srf end

          do k=1,m1
             ut (k,i,j) = ut(k,i,j) + vctr10(k)*(vctr1(k)-up   (k,i,j))
             vt (k,i,j) = vt(k,i,j) + vctr11(k)*(vctr2(k)-vp   (k,i,j))
             tht(k,i,j) =tht(k,i,j) + vctr12(k)*(vctr3(k)-theta(k,i,j))
             pt (k,i,j) = pt(k,i,j) + vctr13(k)*(vctr5(k)-pp   (k,i,j))
             rtt(k,i,j) =rtt(k,i,j) + vctr14(k)*(vctr4(k)-rtp  (k,i,j))
          enddo

       enddo
    enddo

    !Comunicate the local errors to master and let the master compute RMSE and BIAS
    if(evaluate==1)  call comunicateStatistic(soma,somaQ)

    return
  end subroutine nudge

  !     ******************************************************************

  subroutine VariableWeight(mzp, mxp, myp, nxp, nyp, i0, j0, &
       topt, rtgt, varwts,varwts_for_operations_only)
    integer, intent(in)  :: mzp
    integer, intent(in)  :: mxp
    integer, intent(in)  :: myp
    integer, intent(in)  :: nxp
    integer, intent(in)  :: nyp
    integer, intent(in)  :: i0
    integer, intent(in)  :: j0
    real,    intent(in)    :: topt(mxp,myp)
    real,    intent(in)    :: rtgt(mxp,myp)
    real,    intent(inout) :: varwts(mzp,mxp,myp)
    !srf - special weights for pressure (only for operations) 
    real,    intent(inout) :: varwts_for_operations_only(mzp,mxp,myp)

    integer :: i,j,k
    integer :: iGlobal
    integer :: jGlobal
    real :: tnudcenti,tnudtopi,tnudlati
    real :: tnudtopi_x,tnudlati_x
    real :: rown,rows,rowe,roww
    real :: zloc,wttop,wtlat,delzi,centFact,delPBL,wtbot
    character(len=*), parameter :: h="**(VariableWeight)**"
    real, parameter :: znudbot=1000.  ! meters above local terrain

    !         Get weights for large scale and model tendencies
    if (nudlat .le. 0) return

    ! 1) time = 0. do it anyway (needs for initialization)
    if(applyIAU > 0)  then 
       if((timeWindowIAU==0.0 .and. time>0.0) .or. (timeWindowIAU>0.01 .and. time>timeWindowIAU+RAMP)) return 
    endif

    !--- t     => t+ iau => nudging at center with centFact =1 
    !--- t+iau => t+ iau + ramp => nudging decays to zero in this interval with centFact 1 -> 0 


    if(timeWindowIAU>0.01 .and. time <= timeWindowIAU+RAMP .and. RAMP > 0.) then
       !centFact=-1*( (time-timeWindowIAU) /(timeWindowIAU+RAMP)) + 1
       centFact=-1*( (time-timeWindowIAU) /(RAMP))     + 1
    else
       centFact=0.0
    endif

    !-- this is for the case of requesting nudging at center of model domain 
    !-- for the total time integration
    if(applyIAU == 0)  centFact=1.0

    !if(mynum == 1) print*,"ND=",time,timeWindowIAU,timeWindowIAU+RAMP,centFact;call flush(6)

    tnudcenti=0.
    if(tnudcent.gt. .01) tnudcenti=centFact*(1./tnudcent) 
    tnudtopi=0.
    if(tnudtop.gt. .01) then 
       tnudtopi  =1./tnudtop-tnudcenti
       tnudtopi_x=1./tnudtop
    endif
    tnudlati=0.
    if(tnudlat.gt. .01) then 
       tnudlati  =1./tnudlat-tnudcenti
       tnudlati_x=1./tnudlat
    endif

    if(ztop.gt.znudtop) then
       delzi=1./(ztop-znudtop)
    elseif(tnudtop.gt. .01) then
       print*,'Incorrect specification of znudtop ! ! !'
       print*,' znudtop = ',znudtop
       print*,'    ztop = ',ztop
       stop 'varwt-znud'
    endif

    if(time <= dtlt + timeWindowIAU .and. tnudcent > .01 .and. mchnum == master_num) then
       !print*,"=> doing nudging in center:",time,centFact,tnudcenti
       write(*,100) "===> doing nudging in center domain              :",time,centFact,tnudcenti
       call flush(6)
    endif
    if(time >  dtlt + timeWindowIAU .and. tnudcent > .01 .and. mchnum == master_num) then
       !print*,"=> ramping the nudging in center to zero:",time,centFact,tnudcenti
       write(*,100) "===> ramping the nudging in center domain to zero:",time,centFact,tnudcenti
       call flush(6)
    endif

    !-- turn off nudging from surface to znudbot
    delPBL= 0.
    if(znudbot>0.01) delPBL= 1./znudbot

    do j=1,myp
       jGlobal = j + j0
       do i=1,mxp
          iGlobal = i + i0

          ! quadratic weight function for lateral boundaries

          rown=max(0.,float(jGlobal+nudlat-nyp))
          rows=max(0.,float(nudlat+1-jGlobal))
          rowe=max(0.,float(iGlobal+nudlat-nxp))
          roww=max(0.,float(nudlat+1-iGlobal))
          wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
               /float(nudlat*nudlat)

          ! linear weight function for top boundary

          do k=1,mzp
             zloc=zt(k)*rtgt(i,j)+topt(i,j)
             wttop=max(0.,(zloc-znudtop)*delzi)
             !-- turn off nudging from surface to znudbot
             if(znudbot>0.01) then 
                wtbot =min(1.,max(0.,(zloc-znudbot+topt(i,j))*delPBL))
             else
                wtbot =1.
             endif

             ! if( mchnum == master_num) then
             !    if(i==10 .and. j==10) print*,'wtbot=',k,wtbot,zloc-topt(i,j)
             !    call flush(6)
             ! endif

             ! full 3-D weight function

             !-srf orig varwts(k,i,j)=tnudcenti       + max(tnudlati*wtlat,tnudtopi*wttop)
             !-srf excluding nudging in the PBL (wtbot)
             varwts(k,i,j)=tnudcenti*wtbot + max(tnudlati*wtlat,tnudtopi*wttop)

             !-srf special weights for pressure (pressure will be nudged in the full domain)          
             varwts_for_operations_only(k,i,j) = tnudlati_x*wtbot + max(tnudlati_x*wtlat,tnudtopi_x*wttop)          
          end do
       end do
    end do

    !PRINT*,"MX",maxval(varwts),maxval(varwts_for_operations_only)
    !PRINT*,"MN",minval(varwts),minval(varwts_for_operations_only)

100 format(1x,A50,F10.1,F8.2,E12.3)

  end subroutine VariableWeight




  subroutine nudge_cond(m1,m2,m3,ia,iz,ja,jz,varwts  &
       ,varrph,varcph,varrfh,varcfh,rtp,rtt)
    integer :: m1,m2,m3,ia,iz,ja,jz
    real, dimension(m1,m2,m3) :: varrph,varcph,varrfh,varcfh  &
         ,varwts,rtp,rtt

    integer :: i,j,k
    real :: tfact, wt_rc



    ! tfact is temporal interpolation weight,
    ! wt_rc is timescale and grid-dependent weight

    tfact=(time-condtime1)/(condtime2-condtime1)   &
         * wt_nudgec_grid(ngrid)/t_nudge_rc

    wt_rc=wt_nudgec_grid(ngrid)/t_nudge_rc

    do j=ja,jz
       do i=ia,iz

          do k=1,m1
             vctr1(k)=varrph(k,i,j)+(varrfh(k,i,j)-varrph(k,i,j))*tfact
             vctr2(k)=varcph(k,i,j)+(varcfh(k,i,j)-varcph(k,i,j))*tfact

             vctr3(k)=varwts(k,i,j)*wt_rc
          enddo

          ! Only nudging total water where condensate exists...

          do k=1,m1
             if (vctr2(k) > 0.)  &
                  rtt(k,i,j)=rtt(k,i,j) + vctr3(k)*(vctr1(k)-rtp(k,i,j))
          enddo

       enddo
    enddo

    return
  end subroutine nudge_cond

  !     ****************************************************************

  subroutine varweight(n1,n2,n3,varwts,topt,rtgt)
    integer :: n1,n2,n3
    real :: varwts(n1,n2,n3),topt(n2,n3),rtgt(n2,n3)

    integer :: i,j,k
    real :: tnudcenti,tnudtopi,tnudlati,rown,rows,rowe,roww,zloc,wttop &
         ,wtlat,delzi

    !         Get weights for large scale and model tendencies
    if (nudlat .le. 0) return

    tnudcenti=0.
    if(tnudcent.gt. .01) tnudcenti=1./tnudcent
    tnudtopi=0.
    if(tnudtop.gt. .01) tnudtopi=1./tnudtop-tnudcenti
    tnudlati=0.
    if(tnudlat.gt. .01) tnudlati=1./tnudlat-tnudcenti

    if(ztop.gt.znudtop) then
       delzi=1./(ztop-znudtop)
    elseif(tnudtop.gt. .01) then
       print*,'Incorrect specification of znudtop ! ! !'
       print*,' znudtop = ',znudtop
       print*,'    ztop = ',ztop
       stop 'varwt-znud'
    endif


    do j=1,n3
       do i=1,n2

          !                       quadratic weight function for lateral boundaries

          rown=max(0.,float(j+nudlat-n3))
          rows=max(0.,float(nudlat+1-j))
          rowe=max(0.,float(i+nudlat-n2))
          roww=max(0.,float(nudlat+1-i))
          wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
               /float(nudlat*nudlat)

          !                       linear weight function for top boundary

          !--(DMK-CCATT-INI)-----------------------------------------------------
          !srf-rams60 mod
          do k=1,n1
             zloc=ztn(k,1)*rtgt(i,j)+topt(i,j)
             !srf-rams60 mod
             !--(DMK-original)------------------------------------------------------
             !      do k=1,nzp
             !         zloc=zt(k)*rtgt(i,j)+topt(i,j)
             !--(DMK-CCATT-END)-----------------------------------------------------

             wttop=max(0.,(zloc-znudtop)*delzi)

             !                       full 3-D weight function

             varwts(k,i,j)=tnudcenti  &
                  +max(tnudlati*wtlat,tnudtopi*wttop)

          enddo

       enddo
    enddo

    return
  end subroutine varweight


  !     *****************************************************************

  subroutine vfintrpf(ifm,ifflag)
    real :: scr1(maxnxp*maxnyp)
    real :: scr2(maxnxp*maxnyp)
    integer :: ifm,ifflag,icm

    icm = nxtnest(ifm)
    if (icm == 0) return

    !    Temporarily fill VT2DA with interpolated topography from coarser grid

    call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
         ,scr1,grid_g(icm)%topt)
    call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
         ,scr2,scratch%vt2da)

  end subroutine vfintrpf




  subroutine VarfIntrp(ifm,ifflag)
    include "constants.h"
    integer :: ifm,ifflag,icm
    real :: scr1(maxnxp*maxnyp)
    real :: scr2(maxnxp*maxnyp)
    character(len=*), parameter :: h="**(VarfIntrp)**"

    !--(DMK-CCATT-INI)-----------------------------------------------------
    integer :: nspc
    !--(DMK-CCATT-FIM)-----------------------------------------------------

    icm = nxtnest(ifm)
    if (icm == 0) return

    !    Temporarily fill VT2DA with interpolated topography from coarser grid
    !**(JP)** not converted

    !call fatal_error(h//"**(JP)** not converted")
    iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
         "**(JP)** not converted")

    call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
         ,scr1,grid_g(icm)%topt)
    call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
         ,scr2,scratch%vt2da)

  end subroutine VarfIntrp

  !--(DMK-CCATT-INI)-----------------------------------------------------------
  subroutine chem_nudge(m1,m2,m3,ia,iz,ja,jz,nspc &
       ,varwts_chem  &
       ,sc_pp  &
       ,sc_pf  &
       ,scp    &
       ,sct)
    integer :: m1,m2,m3,ia,iz,ja,jz,nspc
    real, dimension(m1,m2,m3) :: sc_pp,sc_pf,varwts_chem,scp,sct

    integer :: i,j,k
    real :: tfact,  wt_chem


    !return

    !         Linearly interpolate values in time, then nudge.

    if (nud_type == 2 .or. nud_type == 4) then
       tfact=(time-vtime1)/(vtime2-vtime1)
    elseif (nud_type == 1) then
       tfact=(time-htime1)/(htime2-htime1)
    endif

    ! - specific nudging  weight for specie 'nspc'
    !- for now the weight is 1 for all species.
    wt_chem = 0.083333 !=1/12 !valor original= 1.
    !- in the future, use different weights and pass them to the slaves via chem_mpass_init
    !!wt_chem=wt_nudge_grid(ngrid)* wt_nudge_chem(nspc)

    do j=ja,jz
       do i=ia,iz

          do k=1,m1
             vctr4(k)=sc_pp(k,i,j)+(sc_pf(k,i,j)-sc_pp(k,i,j))*tfact
             vctr14(k)=varwts_chem(k,i,j)* wt_chem
          enddo

          do k=1,m1
             sct(k,i,j)=sct(k,i,j) + vctr14(k)*(vctr4(k)-scp(k,i,j))
          enddo

       enddo
    enddo

    return
  end subroutine chem_nudge



  subroutine VariableWeightChem(mzp, mxp, myp, nxp, nyp, i0, j0, &
       topt, rtgt, varwts)
    integer, intent(in)  :: mzp
    integer, intent(in)  :: mxp
    integer, intent(in)  :: myp
    integer, intent(in)  :: nxp
    integer, intent(in)  :: nyp
    integer, intent(in)  :: i0
    integer, intent(in)  :: j0
    real,    intent(out) :: varwts(mzp,mxp,myp)
    real,    intent(in)  :: topt(mxp,myp)
    real,    intent(in)  :: rtgt(mxp,myp)

    integer :: i,j,k
    integer :: iGlobal
    integer :: jGlobal
    real :: tnudcenti,tnudtopi,tnudlati
    real :: rown,rows,rowe,roww
    real :: zloc,wttop,wtlat,delzi
    character(len=*), parameter :: h="**(VariableWeight)**"

    !         Get weights for large scale and model tendencies

    if (nudlat .le. 0) return

    tnudcenti=0.
    if(tnudcent.gt. .01) tnudcenti=1./tnudcent
    tnudtopi=0.
    if(tnudtop.gt. .01) tnudtopi=1./tnudtop-tnudcenti
    tnudlati=0.
    if(tnudlat.gt. .01) tnudlati=1./tnudlat-tnudcenti

    if(ztop.gt.znudtop) then
       delzi=1./(ztop-znudtop)
    elseif(tnudtop.gt. .01) then
       print*,'Incorrect specification of znudtop ! ! !'
       print*,' znudtop = ',znudtop
       print*,'    ztop = ',ztop
       stop 'varwt-znud'
    endif


    do j=1,myp
       jGlobal = j + j0
       do i=1,mxp
          iGlobal = i + i0

          ! quadratic weight function for lateral boundaries

          rown=max(0.,float(jGlobal+nudlat-nyp))
          rows=max(0.,float(nudlat+1-jGlobal))
          rowe=max(0.,float(iGlobal+nudlat-nxp))
          roww=max(0.,float(nudlat+1-iGlobal))
          wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
               /float(nudlat*nudlat)

          ! linear weight function for top boundary

          do k=1,mzp
             zloc=ztn(k,1)*rtgt(i,j)+topt(i,j)

             wttop=max(0.,(zloc-znudtop)*delzi)

             ! full 3-D weight function

             varwts(k,i,j)=tnudcenti  &
                  +max(tnudlati*wtlat,tnudtopi*wttop)
          end do
       end do
    end do
  end subroutine VariableWeightChem



  subroutine varweight_chem(n1,n2,n3,varwts_chem,topt,rtgt)
    integer :: n1,n2,n3
    real :: varwts_chem(n1,n2,n3),topt(n2,n3),rtgt(n2,n3)

    integer :: i,j,k
    real :: tnudcenti,tnudtopi,tnudlati,rown,rows,rowe,roww,zloc,wttop &
         ,wtlat,delzi

    !         Get weights for large scale and model tendencies

    if (nudlat .le. 0) return

    tnudcenti=0.
    !-srf
    !- no chemistry nudging in the inner of domain is allowed
    !if(tnudcent.gt. .01) tnudcenti=1./tnudcent

    tnudtopi=0.
    if(tnudtop.gt. .01) tnudtopi=1./tnudtop-tnudcenti
    tnudlati=0.
    if(tnudlat.gt. .01) tnudlati=1./tnudlat-tnudcenti



    if(ztop.gt.znudtop) then
       delzi=1./(ztop-znudtop)
    elseif(tnudtop.gt. .01) then
       print*,'Incorrect specification of znudtop ! ! !'
       print*,' znudtop = ',znudtop
       print*,'    ztop = ',ztop
       stop 'varwt-znud'
    endif


    do j=1,n3
       do i=1,n2

          !                       quadratic weight function for lateral boundaries

          rown=max(0.,float(j+nudlat-n3))
          rows=max(0.,float(nudlat+1-j))
          rowe=max(0.,float(i+nudlat-n2))
          roww=max(0.,float(nudlat+1-i))
          wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
               /float(nudlat*nudlat)

          !                       linear weight function for top boundary

          do k=1,nzp

             !srf-rams60
             zloc=ztn(k,1)*rtgt(i,j)+topt(i,j)
             !         zloc=zt(k)*rtgt(i,j)+topt(i,j)

             wttop=max(0.,(zloc-znudtop)*delzi)

             !                       full 3-D weight function

             varwts_chem(k,i,j)=tnudcenti  &
                  +max(tnudlati*wtlat,tnudtopi*wttop)

          enddo

       enddo
    enddo

    return
  end subroutine varweight_chem
  !srf-chem-end
  !--(DMK-CCATT-FIM)-----------------------------------------------------------
end module ModNudAnalysis
