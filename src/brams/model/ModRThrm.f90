!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModRThrm

  use ModBasicFields, only: &
       BasicFields

  use ModMicControl, only: &
       MicControl
  
  use mem_grid, only: &
       ngrid !INTENT(IN)

  use ModMicroFields, only: &
       MicroFields

  use mem_scratch, only: &
       scratch, &
       vctr1,   &
       vctr2,   &
       vctr3

  use rconstants, only: &
       cpi,  & ! INTENT(IN)
       p00,  & ! INTENT(IN)
       alvl, & ! INTENT(IN)
       alvi, & ! INTENT(IN)
       cpi4, & ! INTENT(IN)
       cp253i, &  ! INTENT(IN)
       cp      ! INTENT(IN)

  implicit none

  private
  public :: thermo
  public :: thermo_boundary_driver
  public :: theta_thp_rk

contains




  subroutine thermo(mzp, mxp, myp, ia, iz, ja, jz, &
       oneBasic, oneMicVars, oneMicroFields)
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    type(BasicFields), pointer, intent(in) :: oneBasic
    type(MicControl), pointer, intent(in) :: oneMicVars
    type(MicroFields), pointer, intent(in) :: oneMicroFields

    character(len=*), parameter :: h="**(thermo)**"
    real :: vctr4(mzp)
    real :: vctr5(mzp)
    real :: vctr6(mzp)
    
    if (oneMicVars%level .le. 1) then

       call drythrm(mzp,mxp,myp,ia,iz,ja,jz  &
            ,oneBasic%thp ,oneBasic%theta   &
            ,oneBasic%rtp ,oneBasic%rv,oneMicVars%level)

    elseif (oneMicVars%level .eq. 2) then

       call satadjst(mzp,mxp,myp,ia,iz,ja,jz  &
            ,oneBasic%pp  ,scratch%scr1             &
            ,oneBasic%thp ,oneBasic%theta &
            ,scratch%vt3db      ,oneBasic%pi0   &
            ,oneBasic%rtp ,oneBasic%rv    &
            ,oneMicroFields%rcp )

    elseif (oneMicVars%level .eq. 3) then

       if(oneMicVars%mcphys_type == 0) then

          call wetthrm3(mzp,mxp,myp,ia,iz,ja,jz,oneMicVars%jnmb  &
               ,oneBasic%pi0 ,oneBasic%pp     &
               ,oneBasic%thp ,oneBasic%theta  &
               ,oneBasic%rtp ,oneBasic%rv     &
               ,oneMicroFields%rcp ,oneMicroFields%rrp    &
               ,oneMicroFields%rpp ,oneMicroFields%rsp    &
               ,oneMicroFields%rap ,oneMicroFields%rgp    &
               ,oneMicroFields%rhp ,oneMicroFields%q6     &
               ,oneMicroFields%q7 &
               ,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,ngrid,oneMicVars%mcphys_type)

       elseif(oneMicVars%mcphys_type >= 2 .or. oneMicVars%mcphys_type == 3 .or. oneMicVars%mcphys_type == 4) then

          !-srf for GThompson/GFDL uphysics
          call wetthrm3_GT(mzp,mxp,myp,ia,iz,ja,jz,oneMicVars%jnmb  &
               ,oneBasic%pi0 ,oneBasic%pp     &
               ,oneBasic%thp ,oneBasic%theta  &
               ,oneBasic%rtp ,oneBasic%rv     &
               ,oneMicroFields%rcp ,oneMicroFields%rrp    &
               ,oneMicroFields%rpp ,oneMicroFields%rsp    &
               ,oneMicroFields%rgp    &
               ,ngrid,oneMicVars%mcphys_type)
       endif


    else

       call fatal_error(h//' Thermo option not supported...LEVEL out of bounds')

    endif

    
  end subroutine thermo



  subroutine thermo_boundary_driver(time, dtlong, f_thermo_e, f_thermo_w, f_thermo_s, &
       f_thermo_n, mzp, mxp, myp, jdim, &
       oneBasic, oneMicVars, oneMicroFields)

    ! Arguments
    real, intent(in) :: time
    real, intent(in) :: dtlong
    logical, intent(in) :: f_thermo_n
    logical, intent(in) :: f_thermo_s
    logical, intent(in) :: f_thermo_e
    logical, intent(in) :: f_thermo_w
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: jdim
    type(BasicFields), pointer, intent(in) :: oneBasic
    type(MicControl), pointer, intent(in) :: oneMicVars
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    
    ! Local Variables
    ! real, parameter :: frq_thermo_bd = 100. !in seconds

    ! Checkig frequency of calls
    ! Calling at each Long Timestep
    if (mod(time,dtlong)/=0) then
       return
    end if

    if (f_thermo_e) then
       call thermo(mzp, mxp, myp, 1,   1,   1, myp, &
            oneBasic, oneMicVars, oneMicroFields)
    end if

    if (f_thermo_w) then
       call thermo(mzp, mxp, myp, mxp, mxp, 1, myp, &
            oneBasic, oneMicVars, oneMicroFields)
    end if

    if (jdim==1) then

       if (f_thermo_s) then
          call thermo(mzp, mxp, myp, 1, mxp, 1, 1, &
               oneBasic, oneMicVars, oneMicroFields)
       end if

       if (f_thermo_n) then
          call thermo(mzp, mxp, myp, 1, mxp, myp, myp, &
               oneBasic, oneMicVars, oneMicroFields)
       end if
    endif

  end subroutine thermo_boundary_driver

  !--------------------------------------------------------------------------------


  !
  !     ***************************************************************
  !
  subroutine drythrm(m1,m2,m3,ia,iz,ja,jz,thil,theta,rt,rv,level)

    ! This routine calculates theta and rv for the case where no condensate is
    ! allowed.

    ! Arguments:
    integer, intent(in)   :: m1,m2,m3,ia,iz,ja,jz,level
    real, intent(in)      :: thil(m1,m2,m3),rt(m1,m2,m3)
    real, intent(inout)   :: theta(m1,m2,m3), rv(m1,m2,m3)

    ! Local Variables:

    integer               :: i,j,k

    do j = ja,jz
       do i = ia,iz
          do k = 1,m1
             theta(k,i,j) = thil(k,i,j)
          enddo
          if (level .eq. 1) then
             do k = 1,m1
                rv(k,i,j) = rt(k,i,j)
             enddo
          endif
       enddo
    enddo
    return
  end subroutine drythrm
  !
  !     ***************************************************************
  !
  subroutine satadjst(m1,m2,m3,ia,iz,ja,jz  &
       ,pp,p,thil,theta,t,pi0,rtp,rv,rcp)

    ! This routine diagnoses theta, rv, and rcp using a saturation adjustment
    ! for the case when water is in the liquid phase only

    ! Arguments:
    integer, intent(in) :: m1, m2, m3, ia, iz, ja, jz
    real, intent(in)    :: pp(m1,m2,m3), thil(m1,m2,m3)  &
         ,pi0(m1,m2,m3), rtp(m1,m2,m3)
    real, intent(inout) ::  p(m1,m2,m3), t(m1,m2,m3), &
         rcp(m1,m2,m3), theta(m1,m2,m3), rv(m1,m2,m3)

    ! Local Variables:
    real ::  rvls
    real, external      :: rslf
    integer             :: i,j,k,iterate
    real                :: picpi,til,tt
    integer             :: n

    do j = ja,jz
       do i = ia,iz
          do k = 1,m1
             picpi = (pi0(k,i,j) + pp(k,i,j)) * cpi
             p(k,i,j) = p00 * picpi ** 3.498
             til = thil(k,i,j) * picpi
             t(k,i,j) = til

             do iterate = 1,20
                rvls = rslf(p(k,i,j),t(k,i,j))
                rcp(k,i,j) = max(rtp(k,i,j) - rvls, 0.)
                tt = 0.7 * t(k,i,j) + 0.3 * til  &
                     * (1. + alvl * rcp(k,i,j)  &
                     / (cp * max(t(k,i,j),253.)))
                if (abs(tt - t(k,i,j)) .le. 0.001) go to 1
                t(k,i,j) = tt
             enddo
1            continue
             rv(k,i,j) = rtp(k,i,j) - rcp(k,i,j)
             theta(k,i,j) = t(k,i,j) / picpi
          enddo
       enddo
    enddo
    return
  end subroutine satadjst
  !
  !     ***************************************************************
  !
  subroutine wetthrm3(m1,m2,m3,ia,iz,ja,jz,jnmb  &
       ,pi0,pp,thp,theta,rtp,rv,rcp,rrp,rpp,rsp,rap,rgp,rhp,q6,q7  &
       ,picpi,tair,til,rliq,rice,qhydm,ngrid,mcphys_type)

    ! This routine calculates theta and rv for "level 3 microphysics"
    ! given prognosed theta_il, cloud, rain, pristine ice, snow, aggregates,
    ! graupel, hail, q6, and q7.

    ! Arguments:
    integer, intent(in)  :: m1, m2, m3, ia, iz, ja, jz, jnmb(*), ngrid,mcphys_type
    real , intent(in)    :: pi0(m1,m2,m3), pp(m1,m2,m3), thp(m1,m2,m3)  &
         ,rtp(m1,m2,m3), rcp(m1,m2,m3), rrp(m1,m2,m3), rpp(m1,m2,m3),   &
         rsp(m1,m2,m3), rap(m1,m2,m3), rgp(m1,m2,m3), rhp(m1,m2,m3),    &
         q6(m1,m2,m3), q7(m1,m2,m3)
    real , intent(inout) :: picpi(*), tair(*), til(*), rliq(m1), rice(*), &
         qhydm(*), rv(m1,m2,m3), theta(m1,m2,m3)

    ! Local Variables:
    integer :: i, j, k
    real    :: tcoal, fracliq, tairstr

    do j = ja,jz
       do i = ia,iz

          do k = 1,m1
             picpi(k) = (pi0(k,i,j) + pp(k,i,j)) * cpi
             tair(k) = theta(k,i,j) * picpi(k)
             til(k) = thp(k,i,j) * picpi(k)
             rliq(k) = 0.
             rice(k) = 0.
          enddo
          if (jnmb(1) .ge. 1) then
             do k = 1,m1
                rliq(k) = rliq(k) + rcp(k,i,j)
             enddo
          endif

          if (jnmb(2) .ge. 1) then
             do k = 1,m1
                rliq(k) = rliq(k) + rrp(k,i,j)
             enddo
          endif

          if (jnmb(3) .ge. 1) then
             do k = 1,m1
                rice(k) = rice(k) + rpp(k,i,j)
             enddo
          endif

          if (jnmb(4) .ge. 1) then
             do k = 1,m1
                rice(k) = rice(k) + rsp(k,i,j)
             enddo
          endif

          if (jnmb(5) .ge. 1) then
             do k = 1,m1
                rice(k) = rice(k) + rap(k,i,j)
             enddo
          endif

          if (jnmb(6) .ge. 1) then
             if(mcphys_type == 0) then 
                do k = 1,m1
                   call qtc(q6(k,i,j),tcoal,fracliq)
                   rliq(k) = rliq(k) + rgp(k,i,j) * fracliq
                   rice(k) = rice(k) + rgp(k,i,j) * (1. - fracliq)
                enddo
             elseif(mcphys_type ==2 .or. mcphys_type ==3.or. mcphys_type ==4) then 
                !-for GThompson uphysics
                do k = 1,m1
                   rice(k) = rice(k) + rgp(k,i,j) 
                enddo
             endif
          endif

          if (jnmb(7) .ge. 1) then
             do k = 1,m1
                call qtc(q7(k,i,j),tcoal,fracliq)
                rliq(k) = rliq(k) + rhp(k,i,j) * fracliq
                rice(k) = rice(k) + rhp(k,i,j) * (1. - fracliq)
             enddo
          endif

          do k = 1,m1
             qhydm(k) = alvl * rliq(k) + alvi * rice(k)
             rv(k,i,j) = rtp(k,i,j) - rliq(k) - rice(k)
          enddo

          do k = 1,m1
             if (tair(k) .gt. 253.) then
                tairstr = 0.5 * (til(k)  &
                     + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
             else
                tairstr = til(k) * (1. + qhydm(k) * cp253i)
             endif
             theta(k,i,j) = tairstr / picpi(k)
          enddo

       enddo
    enddo
    return
  end subroutine wetthrm3
  !
  !     ***************************************************************
  !
  subroutine wetthrm3_GT(m1,m2,m3,ia,iz,ja,jz,jnmb  &
       ,pi0,pp,thp,theta,rtp,rv,rcp,rrp,rpp,rsp,rgp,ngrid,mcphys_type)

    ! This routine calculates theta and rv for "level 3 microphysics"
    ! given prognosed theta_il, cloud, rain, pristine ice, snow, graupel


    ! Arguments:
    integer, intent(in)  :: m1, m2, m3, ia, iz, ja, jz, jnmb(*), ngrid,mcphys_type
    real , intent(in)    :: pi0(m1,m2,m3), pp(m1,m2,m3), thp(m1,m2,m3)  &
         ,rtp(m1,m2,m3), rcp(m1,m2,m3), rrp(m1,m2,m3), rpp(m1,m2,m3),   &
         rsp(m1,m2,m3),  rgp(m1,m2,m3)
    real , intent(inout) ::  rv(m1,m2,m3), theta(m1,m2,m3)

    ! Local Variables:
    integer :: i, j, k
    real    :: tcoal, fracliq, tairstr
    real ,dimension(m1) :: picpi, tair, til, rliq, rice,  qhydm

    do j = ja,jz
       do i = ia,iz

          do k = 1,m1
             picpi(k) = (pi0(k,i,j) + pp(k,i,j)) * cpi
             tair(k) = theta(k,i,j) * picpi(k)
             til(k) = thp(k,i,j) * picpi(k)
             rliq(k) = 0.
             rice(k) = 0.
          enddo
          if (jnmb(1) .ge. 1) then
             do k = 1,m1
                rliq(k) = rliq(k) + rcp(k,i,j)
             enddo
          endif

          if (jnmb(2) .ge. 1) then
             do k = 1,m1
                rliq(k) = rliq(k) + rrp(k,i,j)
             enddo
          endif

          if (jnmb(3) .ge. 1) then
             do k = 1,m1
                rice(k) = rice(k) + rpp(k,i,j)
             enddo
          endif

          if (jnmb(4) .ge. 1) then
             do k = 1,m1
                rice(k) = rice(k) + rsp(k,i,j)
             enddo
          endif


          if (jnmb(6) .ge. 1) then
             do k = 1,m1
                rice(k) = rice(k) + rgp(k,i,j) 
             enddo
          endif

          do k = 1,m1
             qhydm(k) = alvl * rliq(k) + alvi * rice(k)
             rv(k,i,j) = rtp(k,i,j) - rliq(k) - rice(k)
          enddo

          do k = 1,m1
             if (tair(k) .gt. 253.) then
                tairstr = 0.5 * (til(k)  &
                     + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
             else
                tairstr = til(k) * (1. + qhydm(k) * cp253i)
             endif
             theta(k,i,j) = tairstr / picpi(k)
          enddo

       enddo
    enddo
    return
  end subroutine wetthrm3_GT





  
  subroutine theta_thp_rk(mzp,mxp,myp,ia,iz,ja,jz,action, &
       oneBasic, oneMicVars, oneMicroFields)
    !-this is only for RK scheme (uses thc and pc)
    ! Arguments:
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    character(len=*), intent(in) :: action
    type(BasicFields), pointer, intent(in) :: oneBasic
    type(MicControl), pointer, intent(in) :: oneMicVars
    type(MicroFields), pointer, intent(in) :: oneMicroFields

    character(len=*), parameter :: h="**(theta_thp_rk)**"
    real :: vctr4(mzp) 
    real :: vctr5(mzp) 
    real :: vctr6(mzp) 
    
    if (trim(action).ne. "get_thetail" .and. trim(action).ne."get_theta") then
       call fatal_error(h//" unknow action at theta_thp_rk routine")
    end if

    if (oneMicVars%level .le. 1) then
       if (trim(action)=="get_thetail") then
          call fatal_error(h//" not ready for option get_thetail")
       end if

       call drythrm(mzp,mxp,myp,ia,iz,ja,jz  &
            ,oneBasic%thc ,oneBasic%theta   &
            ,oneBasic%rtp ,oneBasic%rv,oneMicVars%level)

    else if (oneMicVars%level .eq. 2) then

       if (trim(action)=="get_thetail") then
          call fatal_error(h//" not ready for option get_thetail")
       end if

       call satadjst(mzp,mxp,myp,ia,iz,ja,jz  &
            ,oneBasic%pc  ,scratch%scr1             &
            ,oneBasic%thc ,oneBasic%theta &
            ,scratch%vt3db          ,oneBasic%pi0   &
            ,oneBasic%rtp ,oneBasic%rv    &
            ,oneMicroFields%rcp )

    else if (oneMicVars%level .eq. 3) then

       if(oneMicVars%mcphys_type == 0) then

          if (trim(action)=="get_thetail") then
             call fatal_error(h//" not ready for option get_thetail")
          end if

          call wetthrm3(mzp,mxp,myp,ia,iz,ja,jz,oneMicVars%jnmb  &
               ,oneBasic%pi0 ,oneBasic%pc     &
               ,oneBasic%thc ,oneBasic%theta  &
               ,oneBasic%rtp ,oneBasic%rv     &
               ,oneMicroFields%rcp ,oneMicroFields%rrp    &
               ,oneMicroFields%rpp ,oneMicroFields%rsp    &
               ,oneMicroFields%rap ,oneMicroFields%rgp    &
               ,oneMicroFields%rhp ,oneMicroFields%q6     &
               ,oneMicroFields%q7 &
               ,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,ngrid,oneMicVars%mcphys_type)

       else if(oneMicVars%mcphys_type == 2 .or. oneMicVars%mcphys_type == 3.or. oneMicVars%mcphys_type == 4) then

          !-srf for GThompson uphysics
          call theta_thp_GT(mzp,mxp,myp,ia,iz,ja,jz,oneMicVars%jnmb  &
               ,oneBasic%pi0 ,oneBasic%pc     &
               ,oneBasic%thc ,oneBasic%theta  &
               ,oneBasic%rtp ,oneBasic%rv     &
               ,oneMicroFields%rcp ,oneMicroFields%rrp    &
               ,oneMicroFields%rpp ,oneMicroFields%rsp    &
               ,oneMicroFields%rgp    &
               ,ngrid,oneMicVars%mcphys_type,action)
       endif


    else

       stop 'theta_thp option not supported...LEVEL out of bounds'

    endif
    
    
  end subroutine theta_thp_rk






  
  subroutine theta_thp_GT(m1,m2,m3,ia,iz,ja,jz,jnmb  &
       ,pi0,pp,thp,theta,rtp,rv,rcp,rrp,rpp,rsp,rgp,ngrid,mcphys_type &
       ,action)

    ! This routine calculates theta for "level 3 microphysics"
    ! given prognosed theta_il, cloud, rain, pristine ice, snow, graupel

    ! Arguments:
    integer, intent(in)  :: m1, m2, m3, ia, iz, ja, jz, jnmb(*), ngrid,mcphys_type
    real , intent(in)    :: pi0(m1,m2,m3), pp(m1,m2,m3)  &
         ,rtp(m1,m2,m3), rcp(m1,m2,m3), rrp(m1,m2,m3), rpp(m1,m2,m3),   &
         rsp(m1,m2,m3),  rgp(m1,m2,m3)
    real , intent(inout) ::  rv(m1,m2,m3), theta(m1,m2,m3),thp(m1,m2,m3)
    character*(*) :: action

    ! Local Variables:
    integer :: i, j, k
    real    :: tcoal, fracliq, tairstr
    real ,dimension(m1) :: picpi, tair, til, rliq, rice,  qhydm

    if(trim(action) == "get_theta") then

       do j = ja,jz
          do i = ia,iz

             do k = 1,m1
                picpi(k) = (pi0(k,i,j) + pp(k,i,j)) * cpi
                tair(k)  = theta(k,i,j) * picpi(k)
                til(k)   = thp(k,i,j) * picpi(k)

                rliq(k) = rcp(k,i,j) + rrp(k,i,j)
                rice(k) = rtp(k,i,j) - rv(k,i,j) - rliq(k)

                qhydm(k) = alvl * rliq(k) + alvi * rice(k)
             enddo

             do k = 1,m1
                if (tair(k) .gt. 253.) then
                   tairstr = 0.5 * (til(k)  &
                        + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
                else
                   tairstr = til(k) * (1. + qhydm(k) * cp253i)
                endif
                theta(k,i,j) = tairstr / picpi(k)
             enddo

          enddo
       enddo

    elseif(trim(action) == "get_thetail") then
       do j = ja,jz
          do i = ia,iz

             do k = 1,m1

                picpi(k) = (pi0(k,i,j) + pp(k,i,j)) * cpi
                tair (k) = theta(k,i,j)*picpi(k)
                rliq (k) = rcp(k,i,j) + rrp(k,i,j)
                rice (k) = rtp(k,i,j) - rv(k,i,j) - rliq(k)


                !- ice-liq potential temperature (Kelvin)
                thp(k,i,j)   =  theta(k,i,j)*(1. + alvl * rliq(k)/(cp * max(tair(k),253.))  &
                     + alvi * rice(k)/(cp * max(tair(k),253.)) ) **(-1.0)      
             enddo

          enddo
       enddo



    else
       stop "unknow action at theta_thp_GT subroutine"
    endif

    return
  end subroutine theta_thp_GT
end module ModRThrm
