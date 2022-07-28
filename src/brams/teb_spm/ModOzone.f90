!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModOzone
  use mem_grid, only: &
       grid_g,        & ! INTENT(IN)
       zt               ! INTENT(IN)

  use ModBasicFields, only: &
       BasicFields

  use ModGaspartFields, only: &
       GaspartFields

  use mem_radiate, only: &
       radiate_g

  use rconstants, only: &
       cpi, &
       cpor, &
       p00 ! INTENT(IN)

  implicit none

  private
  public :: ozone
  
  ! constantes referentes ao calculo de correcao de k em funcao de T, T=300K

  real, parameter :: A(15)=(/&
       1.185E-10, 0.0, 2.04E-05,2.642E+03, 2.055e+02,0.0,3.229E+05, &
       3.068E+02, 4.991E+03, 2.055E+01, 6.1E+03,4.11E+03, 2.936E+04,&
       4.213E+03, 5.578E+02 &
       /)


  real, parameter :: Ea(15)=(/&
       -1.05,  0.0,  0.0, 2.72, 4.91, 0.0, 0.0, 0.0, -0.54,&
       1.19, -1.0,-0.57,  0.0,  0.0, -1.57&
       /)

  real, parameter :: B(15)=(/&
       -2.00, 0.0, -4.80, -1.00, -1.00,0.0,-1.00,-1.00,-1.00,-1.00,-1.00,&
       -1.00, -1.00,-1.00, -1.00&
       /)

  real, parameter :: Rcal   = 0.0019872
  real, parameter :: RJOULE = 8.314
  real, parameter :: M      = 1.0e+06
  real, parameter :: o2     = 2.09e+05
  real, parameter :: Tref   = 300.0
  real, parameter :: Mih2o  = 18.01

  ! unidades das constantes
  ! A=fator de frequencia de choques (ppm-1min-1)
  ! R= constante dos gases (kcal/mol/K)
  ! B= coeficiente de correcao adicionado a equacao de 
  !    Arrhenius para correcao da T (adimensional)
  ! Ea= energia de ativacao (kca/mol)
  ! M=concentracao de ar (ppm)
  ! o2=concentracao do oxigenio constante (ppm)

contains



  subroutine ozone(mzp, mxp, myp, ia, iz, ja, jz, ng, deltat, &
       oneBasicFields, oneGaspartFields)
    ! Arguments:
    integer, intent(in) :: mzp, mxp, myp, ia, iz, ja, jz, ng
    real, intent(in)    :: deltat
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(GaspartFields), pointer, intent(in) :: oneGaspartFields

    call chemistry (mzp, mxp, myp, ia, iz, ja, jz,              &
         oneGaspartFields%pno,   oneGaspartFields%pno2, &
         oneGaspartFields%pco,   oneGaspartFields%pvoc, &
         oneGaspartFields%po3,   oneGaspartFields%pso2, &
         oneGaspartFields%pso4,                             &
         oneGaspartFields%prhco, oneGaspartFields%pho2, &
         oneGaspartFields%po3p,  oneGaspartFields%po1d, &
         oneGaspartFields%pho,   oneGaspartFields%proo, &
         oneBasicFields%theta,   oneBasicFields%dn0,    &
         oneBasicFields%pi0,     oneBasicFields%pp,     &
         oneBasicFields%rv,                                 &
         radiate_g(ng)%cosz,   &
         grid_g(ng)%rtgt,       grid_g(ng)%topma ,    &
         deltat, cpi, cpor, p00, zt,                            &
         oneGaspartFields%pnot,      oneGaspartFields%pno2t,    &
         oneGaspartFields%pcot,      oneGaspartFields%pvoct,    &
         oneGaspartFields%po3t,                                 &
         oneGaspartFields%pso4t,                                &
         oneGaspartFields%prhcot,    oneGaspartFields%pho2t,    &
         oneGaspartFields%po3pt,     oneGaspartFields%po1dt,    &
         oneGaspartFields%phot,      oneGaspartFields%proot     )


    !endif

  end subroutine ozone

  !##############################################################################

  subroutine chemistry(m1, m2, m3, ia, iz, ja, jz,    &
       noi, no2i, coi, vocsi, o3i, so2i, so4i, rchoi, &
       ho2i, O3Pi, O1Di, HOi, RO2i, theta, dn0, pi0,  &
       pp, rv, cosz, rtgt, topt, dtlt, cpi,           & !rshort
       cpor, p00, zt, not, no2t, cot, vocst, o3t,     &
       so4t, rchot, ho2t, O3Pt, O1Dt, HOt, RO2t       ) !so2t
    ! Arguments:
    integer, intent(in) :: m1, m2, m3, ia, iz, ja, jz
    ! chemichal elements and local altitude
    real, intent(in)    :: noi(m1,m2,m3)
    real, intent(in)    :: no2i(m1,m2,m3)
    real, intent(in)    :: coi(m1,m2,m3)
    real, intent(in)    :: vocsi(m1,m2,m3)
    real, intent(in)    :: o3i(m1,m2,m3)
    real, intent(in)    :: so2i(m1,m2,m3)
    real, intent(in)    :: so4i(m1,m2,m3)
    real, intent(in)    :: rchoi(m1,m2,m3)
    real, intent(in)    :: ho2i(m1,m2,m3)
    real, intent(in)    :: O3Pi(m1,m2,m3)
    real, intent(in)    :: O1Di(m1,m2,m3)
    real, intent(in)    :: HOi(m1,m2,m3)
    real, intent(in)    :: RO2i(m1,m2,m3)
    real, intent(inout) :: not(m1,m2,m3)
    real, intent(inout) :: no2t(m1,m2,m3)
    real, intent(inout) :: cot(m1,m2,m3)
    real, intent(inout) :: vocst(m1,m2,m3)
    real, intent(inout) :: o3t(m1,m2,m3)
    real, intent(inout) :: so4t(m1,m2,m3)
    real, intent(inout) :: rchot(m1,m2,m3)
    real, intent(inout) :: ho2t(m1,m2,m3)
    real, intent(inout) :: O3Pt(m1,m2,m3)
    real, intent(inout) :: O1Dt(m1,m2,m3)
    real, intent(inout) :: HOt(m1,m2,m3)
    real, intent(inout) :: RO2t(m1,m2,m3)
    ! 3-D model's variables used
    real, intent(in)    :: theta(m1,m2,m3)
    real, intent(in)    :: dn0(m1,m2,m3)
    real, intent(in)    :: pi0(m1,m2,m3)
    real, intent(in)    :: pp(m1,m2,m3)
    real, intent(in)    :: rv(m1,m2,m3)
    ! 2-D model's variables used
    real, intent(in)    :: cosz(m2,m3)
    real, intent(in)    :: rtgt(m2,m3)
    real, intent(in)    :: topt(m2,m3)
    !time-split determining variables
    real, intent(in)    :: dtlt
    ! model's "constants" used
    real, intent(in)    :: cpi
    real, intent(in)    :: cpor
    real, intent(in)    :: p00
    ! 1-D model's variables used
    real, intent(in)    :: zt(m1)

    ! Local Variables:
    integer :: i, j, ik,jk
    !variables used to keep time n-1 values
    real :: pno, pno2, po3p, po1d, po3, ph2o, pco, pho, pro2, prcho, pvocs, &
         pho2, pso2, pso4, hom
    !chemichal elements and local altitude
    real :: no(m1), no2(m1), co(m1), vocs(m1), o3(m1), rcho(m1), &
         ho2(m1), O3P(m1), O1D(m1), HO(m1), RO2(m1), h2o(m1), so2(m1), so4(m1)
    ! 3-D model's variables used
    real :: tempk(m1,m2,m3)
    ! model's "constants" used
    real, parameter :: avo=6.0221367E+23  !avogadro constant
    real            :: temp
    !velocities coefficients
    real            :: k(15)
    !photolysis rates
    real            :: j2, j6
    real, parameter :: kso2 = 1.0e-18 ! defining velocity constant for SO2
    !production and lost terms used in the implicit scheme
    real :: prod, plosc
    !time-split determining variables
    real :: dtll_factor, dtll, dtlti
    integer :: niter_ozo, inter


    !space loop
    do j=ja,jz
       do i=ia,iz

          !storing previous values for all gases
          do ik=2,m1

             !calculating absolute temperature using Exner function and Theta
             tempk(ik,i,j) = theta(ik,i,j)*pi0(ik,i,j)*cpi

          enddo
       enddo
    enddo

    !calculating time-split for chemical reactions (minimum 1 s)
    niter_ozo   = max(1, nint(dtlt/5.e-1))
    dtll_factor = 1./float(niter_ozo)
    dtll        = dtlt*dtll_factor

    do j=ja,jz
       do i=ia,iz
          do ik=1,m1

             no  (ik) = max(0., noi  (ik,i,j))
             no2 (ik) = max(0., no2i (ik,i,j))
             so2 (ik) = max(0., so2i (ik,i,j))
             so4 (ik) = max(0., so4i (ik,i,j))
             co  (ik) = max(0., coi  (ik,i,j))
             vocs(ik) = max(0., vocsi(ik,i,j))
             o3  (ik) = max(0., o3i  (ik,i,j))
             rcho(ik) = max(0., rchoi(ik,i,j))
             ho2 (ik) = max(0., ho2i (ik,i,j))
             o3p (ik) = max(0., o3pi (ik,i,j))
             o1d (ik) = max(0., o1di (ik,i,j))
             ho  (ik) = max(0., hoi  (ik,i,j))
             ro2 (ik) = max(0., ro2i (ik,i,j))

             h2o(ik)  = rv(ik,i,j)

          enddo

          call conv_ppm_rm(m1, o3  , 48.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, co  , 28.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, no  , 30.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, no2 , 46.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, vocs, 42.08, 28.97, 1.e6)
          call conv_ppm_rm(m1, rcho, 44.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, ho2 , 33.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, h2o , 18.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, O3P , 16.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, O1D , 16.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, HO  , 17.0 , 28.97, 1.e6)
          call conv_ppm_rm(m1, RO2 , 47.0 , 28.97, 1.e6)

          j2 = max(0.001,     42.92e-02*cosz(i,j))/60.
          j6 = max(0.0000001,  2.04E-03*cosz(i,j))/60.

          do ik=2,m1-1

             !calculating k's coefficients
             do jk=1,15

                k(jk) = (A(jk)*((tempk(ik,i,j)/Tref)**B(jk))* &
                     exp(-Ea(jk)/(Rcal*tempk(ik,i,j))))/60.

             enddo

             !initiate time-split for reactions in order to avoid numerical instalibity
             do inter=1,niter_ozo

                !keeping values from time n-1
                pno   = no  (ik)
                pno2  = no2 (ik)
                po3p  = o3p (ik)
                po1d  = o1d (ik)
                po3   = o3  (ik) 
                ph2o  = h2o (ik)
                pco   = co  (ik)
                pho   = ho  (ik)
                pro2  = ro2 (ik)
                prcho = rcho(ik)
                pvocs = vocs(ik)
                pho2  = ho2 (ik)  

                !updating concentrations

                !implicit integration of O3P
                o3p(ik) = (po3p + (j2*pno2*dtll))/(1. + (M*k(3)*(o2)*dtll))

                !implicit integration of O1D
                O1D(ik) = (pO1D + (j6* po3*dtll))/(1. + (k(7)*ph2o*dtll))

                !implicit integration of HO
                plosc   = (k(8)*pco) + (k(11)*pvocs) + (k(13)*prcho)
                prod    = (2.*k(7)*ph2o*pO1D) + (k(9)*pno*pho2) + (k(10)*pho2*po3)
                ho(ik)  = (pho + (prod*dtll))/(1. + (plosc*dtll))

                !implicit integration of ro2
                plosc   = (k(12)*pno) +(k(15)*pho2)
                prod    = (k(11)*pvocs*pho)
                ro2(ik) = (pro2 + (prod*dtll))/(1. + (plosc*dtll))

                !implicit integration of o3
                plosc   = (k(4)*pno) + (k(5)*pno2) + j6 + (k(10)*pho2)
                prod    = k(3)*o2*M*po3p
                o3(ik)  = (po3 + dtll*prod)/(1. + dtll*plosc)

                !implicit integration of no
                plosc   = (2.*k(1)*pno*o2) + (k(4)*po3)+(k(9)*pho2) + (k(12)*pro2)
                prod    = j2*pno2
                no(ik)  = (pno + dtll*prod)/(1. + dtll*plosc)

                !implicit integration of no2
                plosc   = j2 + (k(5)*po3)
                prod    = (2.*k(1)*pno*pno*o2) + (k(4)*pno*po3) + &
                     (k(9)*pno*pho2) +(k(12)*pro2*pno)
                no2(ik) = (pno2 + dtll*prod)/(1. + dtll*plosc)

                !implicit integration of vocs
                plosc   =  k(11)*pho 
                prod    = 0.
                vocs(ik)= (pvocs + dtll*prod)/(1. + dtll*plosc) 

                !implicit integration of co
                plosc   = k(8)*pho
                prod    = 0.
                co(ik)  = (pco + dtll*prod)/(1. + dtll*plosc) 

                !implicit integration of rcho
                plosc   = (k(13)*pho)
                prod    = k(12)*pro2*pno
                rcho(ik)= (prcho + dtll*prod)/(1. + dtll*plosc)

                !implicit integration of ho2
                plosc   = (k(9)*pno) + (k(10)*po3) + (2.*k(14)*pho2) + &
                     (k(15)*pro2)
                prod    = (k(8)*pho*pco) +(k(12)*pro2*pno)
                ho2(ik) = (pho2 + dtll*prod)/(1. + dtll*plosc)

             enddo !time
          enddo !levels

          !converting units from ppm to kg/kg
          call conv_ppm_rm(m1, o3  , 28.97, 48.0 , 1.e-6)
          call conv_ppm_rm(m1, co  , 28.97, 28.0 , 1.e-6)
          call conv_ppm_rm(m1, no  , 28.97, 30.0 , 1.e-6)
          call conv_ppm_rm(m1, no2 , 28.97, 46.0 , 1.e-6)
          call conv_ppm_rm(m1, vocs, 28.97, 42.08, 1.e-6)
          call conv_ppm_rm(m1, rcho, 28.97, 44.0 , 1.e-6)
          call conv_ppm_rm(m1, ho2 , 28.97, 33.0 , 1.e-6)
          call conv_ppm_rm(m1, h2o , 28.97, 18.0 , 1.e-6)
          call conv_ppm_rm(m1, O3P , 28.97, 16.0 , 1.e-6)
          call conv_ppm_rm(m1, O1D , 28.97, 16.0 , 1.e-6)
          call conv_ppm_rm(m1, HO  , 28.97, 17.0 , 1.e-6)
          call conv_ppm_rm(m1, RO2 , 28.97, 47.0 , 1.e-6)

          dtlti = 1./dtlt

          ! integration of so4
          do ik=1,m1

             !converting variables to units of s-1
             pso2    = (so2(ik)*dn0(ik,i,j)*avo)/(1e-3*64.)
             pso4    = (so4(ik)*dn0(ik,i,j)*avo)/(1E-3*96.)
             hom     = (ho(ik)*dn0(ik,i,j)*avo)/(1E-3*17.)
             temp    = pso4 + (kso2*hom*pso2)*dtlt 

             !converting back to kg/kg
             so4(ik) = (temp*(1e-3)*96.)/(dn0(ik,i,j)*avo) 

          enddo

          do ik=1,m1

             so4t (ik,i,j) = so4t (ik,i,j) + ((so4 (ik) - so4i (ik,i,j))*dtlti)
             o3t  (ik,i,j) = o3t  (ik,i,j) + ((o3  (ik) - o3i  (ik,i,j))*dtlti)
             cot  (ik,i,j) = cot  (ik,i,j) + ((co  (ik) - coi  (ik,i,j))*dtlti)
             not  (ik,i,j) = not  (ik,i,j) + ((no  (ik) - noi  (ik,i,j))*dtlti)
             no2t (ik,i,j) = no2t (ik,i,j) + ((no2 (ik) - no2i (ik,i,j))*dtlti)
             vocst(ik,i,j) = vocst(ik,i,j) + ((vocs(ik) - vocsi(ik,i,j))*dtlti)
             rchot(ik,i,j) = rchot(ik,i,j) + ((rcho(ik) - rchoi(ik,i,j))*dtlti)
             ho2t (ik,i,j) = ho2t (ik,i,j) + ((ho2 (ik) - ho2i (ik,i,j))*dtlti)
             o3pt (ik,i,j) = o3pt (ik,i,j) + ((o3p (ik) - o3pi (ik,i,j))*dtlti)
             o1dt (ik,i,j) = o1dt (ik,i,j) + ((o1d (ik) - o1di (ik,i,j))*dtlti)
             hot  (ik,i,j) = hot  (ik,i,j) + ((ho  (ik) - hoi  (ik,i,j))*dtlti)
             ro2t (ik,i,j) = ro2t (ik,i,j) + ((ro2 (ik) - ro2i (ik,i,j))*dtlti)

          enddo

       enddo
    enddo

  end subroutine chemistry

  ! ---------------------------------------------------------------------------

  subroutine conv_ppm_rm(n1, b, pden, pnum, expo)
    ! Arguments: 
    integer, intent(in) :: n1
    real, intent(inout) :: b(n1)
    real, intent(in)    :: pden
    real, intent(in)    :: pnum
    real, intent(in)    :: expo ! Not used
    ! Local variables:
    integer :: k
    real    :: fator


    fator = (pden/pnum)*expo

    do k=1,n1-1
       b(k) = b(k)*fator
    enddo

  end subroutine conv_ppm_rm
end module ModOzone
