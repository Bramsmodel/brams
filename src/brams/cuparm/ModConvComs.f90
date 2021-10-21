!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModConvComs

  implicit none
  public
  
  integer, parameter :: nkp=100

  integer :: kcon
  integer :: klcl
  integer :: klfc
  integer :: ketl
  integer :: kct
  integer :: igo
  integer :: kmt
  integer :: icprtfl
  integer :: icpltfl
  real :: zmid
  real :: cdzmin
  real :: dzlow
  real :: dzhigh
  real :: plcl
  real :: tlcl
  real :: dzlcl
  real :: zlcl
  real :: garea
  real :: wconmin
  real :: contim
  real :: preff
  real :: envshr
  real :: supply
  real :: cptime
  real :: cprecip
  real :: ucon(nkp)
  real :: vcon(nkp)
  real :: wcon(nkp)
  real :: thtcon(nkp)
  real :: rvcon(nkp)
  real :: prcon(nkp)
  real :: picon(nkp)
  real :: tmpcon(nkp)
  real :: dncon(nkp)
  real :: zcon(nkp)
  real :: zzcon(nkp)
  real :: upe(nkp)
  real :: vpe(nkp)
  real :: wpe(nkp)
  real :: ze(nkp)
  real :: te(nkp)
  real :: the(nkp)
  real :: pe(nkp)
  real :: rte(nkp)
  real :: pke(nkp)
  real :: rhoe(nkp)
  real :: thve(nkp)
  real :: zc(nkp)
  real :: rve(nkp)
  real :: thee(nkp)
  real :: qvct1(nkp)
  real :: qvct2(nkp)
  real :: qvct3(nkp)
  real :: qvct4(nkp)
  real :: vheat(nkp)
  real :: vmois(nkp)
  real :: vmdry(nkp)
  real :: frcon(nkp)
  real :: ftcon(nkp)
  real :: tcon(nkp)
  real :: rcon(nkp)
  real :: theu(nkp)
  real :: rsu(nkp)
  real :: thu(nkp)
  real :: tu(nkp)
  real :: thd(nkp)
  real :: wtd(nkp)
  real :: thcon(nkp)
  real :: rtcon(nkp)
end module ModConvComs
