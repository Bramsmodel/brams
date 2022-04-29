!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModRtimi
  use mem_grid, only: &
       ngrid, &
       dyncore_flag, &
       dtlt, &
       time, &
       eps, &
       dtlv, &
       icorflg, &
       jdim, &
       ngbegun, &
       nzp, &
       nxp, &
       nyp

  use mem_tend, only: &
       tend

  use ModScalarTable, only: &
       ScalarTable

  use node_mod, only: &
       mxp, &
       myp, &
       mzp, &
       nmachs

  use mem_grell, only: &
       cuforc_g, &
       cuforc_sh_g

  use shcu_vars_const, only: &
       NNSHCU

  use mem_cuparm, only: &
       confrq, &
       NNQPARM

  use mem_basic, only: &
       basic_g

  use mem_scratch, only: &
       scratch

  implicit none

  include "constants.h"

  private
  public :: tend0
  public :: hadvance
  public :: predtr

contains
  subroutine tend0(oneScalarTab, oneScalarTabSize)
    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(inout) :: oneScalarTabSize

    integer :: n
    integer(kind=i8) :: mxyzp

    !     This routine simply sets all tendency arrays to zero.

    !     First u,v tendencies

    mxyzp = mxp*myp*mzp
    tend%ut(1:mxyzp) = 0.
    tend%vt(1:mxyzp) = 0.
    tend%wt(1:mxyzp) = 0.
    tend%pt(1:mxyzp) = 0.
    !-srf if RK scheme - check if this is necessary -
    !if (dyncore_flag == 2 ) then
    !   tend%ut_rk (:) = 0.0
    !   tend%vt_rk (:) = 0.0
    !   tend%wt_rk (:) = 0.0
    !   tend%pt_rk (:) = 0.0
    !   tend%tht_rk(:) = 0.0
    !endif
    !IF(mod(time,confrq).lt.dtlt .or. time .lt. .01) then  
    !  if(NNQPARM(ngrid) >1) then
    !    cuforc_g   (ngrid)%lsfth(1:mzp,1:mxp,1:myp) = 0. 
    !    cuforc_g   (ngrid)%lsfrt(1:mzp,1:mxp,1:myp) = 0.
    !  endif 
    !  if(NNSHCU(ngrid) >1) then
    !    cuforc_sh_g(ngrid)%lsfth(1:mzp,1:mxp,1:myp) = 0. 
    !    cuforc_sh_g(ngrid)%lsfrt(1:mzp,1:mxp,1:myp) = 0.
    !  endif 
    !ENDIF 
    !     Now sclrr tendencies

    do n = 1,oneScalarTabSize        

       call azero_l(mxyzp, oneScalarTab(n)%var_t_1D)
    enddo

  end subroutine tend0

  !**************************************************************************

  subroutine hadvance(iac)

    integer, intent(in) :: iac

    integer(kind=i8) :: mxyzp
    logical :: RAW
    real :: vt3da(mzp,mxp,myp)

    !     It is here that the Asselin filter is applied.  For the velocities
    !     and pressure, this must be done in two stages, the first when
    !     IAC=1 and the second when IAC=2.
    RAW = .false.
    IF(dyncore_flag == 1) RAW = .true.

    mxyzp = mxp * myp * mzp
    eps = .2

    IF(RAW  .and. ngbegun(ngrid) == 0) then 
       !-srf: march-2012/revisited on Nov2015
       !- Implemented an upgrade of Robert-Asselin (RA) filter based on Williams (MWR, 2009)
       !- The new filter is called Robert-Asselin-Williams (RAW)
       !- The scratch arrays vt3di-j-k-l are used to store the displacment parameter
       !- between the 2 calls of HADVANCE in the subroutine timestep. These arrays
       !- must not be used for any other routine, within this procedure (before the 
       !- fist call of hadvance and after the last one). After that, they are free to
       !- be used for other procedure.
       !- For this filter, the arrays vt3d-i,-j,-k,-l must be allocated with size mxyzp.

       if(min(  size(scratch%vt3di), size(scratch%vt3dj) & 
            ,size(scratch%vt3dk), size(scratch%vt3dl) ) < mxyzp )&
            stop "===> size of vt3di or vt3dj or vt3dk or vt3dl is less than mxyzp in hadvance"  
    ENDIF

    !     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.

    call predict(&
         mxyzp, &
         basic_g(ngrid)%uc, &
         basic_g(ngrid)%up,  &
         tend%ut, &
         vt3da,&
         iac, dtlv,  &
         scratch%vt3di,&
         RAW)

    if (icorflg .eq. 1 .or. jdim .eq. 1) then
       call predict(&
            mxyzp, &
            basic_g(ngrid)%vc, &
            basic_g(ngrid)%vp, &
            tend%vt, &
            vt3da,&
            iac, dtlv, &
            scratch%vt3dj, &
            RAW)
    endif

    call predict(&
         mxyzp, &
         basic_g(ngrid)%wc, &
         basic_g(ngrid)%wp, &
         tend%wt, &
         vt3da,&
         iac, dtlv, &
         scratch%vt3dk, &
         RAW)

    call predict(mxyzp,&
         basic_g(ngrid)%pc, &
         basic_g(ngrid)%pp, &
         tend%pt, &
         vt3da, &
         iac, dtlv, &
         scratch%vt3dl, &
         RAW)


    return
  end subroutine hadvance


  !**************************************************************************


  !**************************************************************************

  subroutine predtr(oneScalarTab, oneScalarTabSize)
    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(inout) :: oneScalarTabSize
    integer :: n !mxyzp
    integer(kind=i8) :: mxyzp
    character(len=*), parameter :: h="**(predtr)**"

    !   -  Step thermodynamic variables from  t  to  t+1.

    mxyzp = mxp * myp * mzp

    do n = 1, oneScalarTabSize

       !- if RK scheme, THP/THC are not predicted here
       if (dyncore_flag == 2) then
          if (oneScalarTab(n)%name == 'THC' .or. &
               oneScalarTab(n)%name == 'THP') cycle
       endif

       if (associated(oneScalarTab(n)%var_p_1D)) then
          call update_long(mxyzp, oneScalarTab(n)%var_p_1D,  &
               oneScalarTab(n)%var_t_1D, dtlt)
       else if (associated(oneScalarTab(n)%var_p_2D)) then
          call update_long(mxyzp, oneScalarTab(n)%var_p_2D,  &
               oneScalarTab(n)%var_t_1D, dtlt)
       else if (associated(oneScalarTab(n)%var_p_3D)) then
          call update_long(mxyzp, oneScalarTab(n)%var_p_3D,  &
               oneScalarTab(n)%var_t_1D, dtlt)
       else
          call fatal_error(h//" no var_p_XD associated for "//trim(oneScalarTab(n)%name))
       end if
    enddo

  end subroutine predtr
end module ModRtimi

subroutine predict(&
     npts,&
     ac, &
     ap, &
     fa, &
     af, &
     iac,dtlp, &
     d, &
     RAW)

  use mem_grid, only: &
       ngrid, &
       ngbegun, &
       eps
  
  implicit none

  include "constants.h"

  integer(kind=i8), intent(in) :: npts
  real, intent(inout) :: ac(:)
  real, intent(inout) :: ap(:)
  real, intent(in) :: fa(:)
  real, intent(out) :: af(:)
  integer, intent(in) :: iac   !,npts, m
  real, intent(in) :: dtlp
  real, intent(inout):: d(:) ! the displacement parameter 
  logical, intent(in) :: RAW
  !     For IAC=3, this routine moves the arrays AC and AP forward by
  !     1 time level by adding in the prescribed tendency. It also
  !     applies the Asselin filter given by:

  !              {AC} = AC + EPS * (AP - 2 * AC + AF)

  !     where AP,AC,AF are the past, current and future time levels of A.
  !     All IAC=1 does is to perform the {AC} calculation without the AF
  !     term present.  IAC=2 completes the calculation of {AC} by adding
  !     the AF term only, and advances AC by filling it with input AP
  !     values which were already updated in ACOUSTC.
  !
  !-srf
  real :: epsu
  integer(kind=i8) :: m

  real   , parameter :: alpha = 0.5   !- use alpha = 0.5 for RAW  filter
  !- use alpha = 1.0 for RA   filter
  !-srf

  epsu = eps
  if (ngbegun(ngrid) .eq. 0) epsu = 0.5

  if(RAW) epsu=0.5*epsu  

  if (iac .eq. 1) then

     if(.not. RAW) then
        ac(1:npts) = ac(1:npts) + epsu * (ap(1:npts) -2. * ac(1:npts))

        !********srf: using RAW filter - start ****************
     elseif(RAW) then

        d(1:npts)  = epsu*(ap(1:npts) - 2.*ac(1:npts))    

     endif
     !********srf: using RAW filter - end   ****************

     return
  elseif (iac .eq. 2) then

     if(.not. RAW) Then
        af(1:npts) = ap(1:npts)                     ! AQUI ap eh o U^(N+1)
        ap(1:npts) = ac(1:npts) + epsu * af(1:npts) ! AQUI AP eh U^(N) filtrado

        !*******srf: using RAW filter - start ****************
     elseif(RAW) then

        do m = 1,npts
           !srf: aqui "ap"  ja e o "af" obtido na acoustic e buoyancy.      
           !     termino do calculo do parametro displacement (Williams 2009)
           d (m) = epsu*ap(m)+d(m)

           !srf: aqui o U(N+1)=af filtrado eh obtido
           af(m) = ap(m)+d(m)*(alpha-1.)

           !srf: aqui o U(N)=ac filtrado eh obtido e ja salvo na area de memoria
           !     do "ap" para o proximo timestep
           ap(m) = ac(m)+d(m)*alpha 

           !srf: aqui o af ou (U_n+1) se torna o "ac" para o proximo timestep
           ac(m) = af(m)        
        enddo
        !srf: aqui tudo que eh necessario para avancar em 2 * deltaT
        !     foi feito => return
        return
        !*******srf: using RAW filter - end ****************

     endif

  elseif (iac .eq. 3) then
     af = ap + dtlp * fa
     ap = ac + epsu * (ap - 2. * ac + af)
  endif

  ac(1:npts) = af(1:npts)   ! aqui af se torna ac para  o proximo timestep

end subroutine predict
