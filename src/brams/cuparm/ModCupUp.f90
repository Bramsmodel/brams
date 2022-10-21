module ModCupUp

  implicit none
  
  private
  public :: cup_up_he
  public :: cup_up_nms
  public :: cup_up_moisture
  public :: cup_up_aa0
contains

  !----------------------------------------------------------------------
  subroutine cup_up_he(k22, hkb, z_cup, cd, entr2, he_cup, hc,               &
       mix, mgmxp, mkx, mgmzp, kbcon, ierr, istart, iend, dby, he, hes_cup)
    integer, intent(in) :: k22(:)
    real, intent(out) :: hkb(:)
    real, intent(in) :: z_cup(:,:)
    real, intent(in) :: cd(:,:)
    real, intent(in) :: entr2(:,:)
    real, intent(in) :: he_cup(:,:)
    real, intent(inout) :: hc(:,:) 
    integer, intent(inout) :: mix !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mkx
    integer, intent(in) :: mgmzp
    integer, intent(in) :: kbcon(:)
    integer, intent(in) :: ierr(:)
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    real, intent(out) :: dby(:,:)
    real, intent(in) :: he(:,:)
    real, intent(in) :: hes_cup(:,:)

    integer i, j, k
    real entr, dz


    ! hc = cloud moist static energy
    ! hkb = moist static energy at originating level
    ! he = moist static energy on model levels
    ! he_cup = moist static energy on model cloud levels
    ! hes_cup = saturation moist static energy on model cloud levels
    ! dby = buoancy term
    ! cd= detrainment function 
    ! z_cup = heights of model cloud levels 
    ! entr = entrainment rate
    !--- Moist static energy inside cloud

    do i=istart,iend
       if (ierr(i).eq.0.) then
          hkb(i)=he_cup(i,k22(i))
          do k=1,k22(i)
             hc(i,k)=he_cup(i,k)
             dby(i,k)=0.
          enddo
          do k=k22(i),kbcon(i)-1
             hc(i,k)=hkb(i)
             dby(i,k)=0.
          enddo
          k=kbcon(i)
          hc(i,k)=hkb(i)
          dby(i,kbcon(i))=hkb(i)-hes_cup(i,k)
       endif
    enddo
    do k=2,mkx-1
       do i=istart,iend
          if (k.gt.kbcon(i).and.ierr(i).eq.0.) then
             dz=z_cup(i,k)-z_cup(i,k-1)
             hc(i,k)=(hc(i,k-1)*(1.-.5*cd(i,k)*dz)+entr2(i,k)*      &
                  dz*he(i,k-1))/(1.+entr2(i,k)*dz-.5*cd(i,k)*dz)
             dby(i,k)=hc(i,k)-hes_cup(i,k)
          endif
       enddo
    enddo
  end subroutine cup_up_he


  !----------------------------------------------------------------------
  subroutine cup_up_moisture(ierr, z_cup, qc, qrc, pw, pwav,            &
       kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend, cd, dby,      &
       mentr_rate2, q, gamma_cup, zu, qes_cup, k22, qe_cup, &
       w_up, rho, ccn, trigg, iens, autoconv)
    integer, intent(in) :: ierr(:)
    real, intent(in) :: z_cup(:,:)
    real, intent(inout) :: qc(:,:)
    real, intent(out) :: qrc(:,:)
    real, intent(out) :: pw(:,:)
    real, intent(inout) :: pwav(:)
    integer, intent(in) :: kbcon(:)
    integer, intent(in) :: ktop(:)
    integer, intent(inout) :: mix !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mkx
    integer, intent(in) :: mgmzp
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    real, intent(in) :: cd(:,:)
    real, intent(in) :: dby(:,:)
    real, intent(in) :: mentr_rate2(:,:)
    real, intent(in) :: q(:,:)
    real, intent(in) :: gamma_cup(:,:)
    real, intent(in) :: zu(:,:)
    real, intent(in) :: qes_cup(:,:)
    integer, intent(in) :: k22(:)
    real, intent(in) :: qe_cup(:,:) 
    real, intent(in) :: w_up(:,:)
    real, intent(in) :: rho(:,:) ! g/cm^3
    real, intent(in) :: ccn(:)
    integer, intent(in) :: trigg
    integer, intent(in) :: iens
    integer, intent(in) :: autoconv


    integer i, k

    integer iall

    real radius, xl, dz, dh, qrch, c0, mentr_rate

    real qtc(mgmxp,mgmzp)
    real con, q1

    ! real, parameter :: anbase =  50. !*1.e+6 !berry-number at cloud base #/m^3(maritime)
    ! real, parameter :: ccn =1000.  !*1.e+6 !berry-number at cloud base #/m^3(continental)
    real, parameter :: bdispm = 0.366       !berry--size dispersion (maritime)
    real, parameter :: bdispc = 0.146       !berry--size dispersion (continental)
    !integer, parameter :: berry = 2     ! if berry = 2, berry parameterization on, else berry = 0 to 
    ! default autoconversion param.

    integer, parameter:: deep=1,shallow=2

    if(autoconv == 2 .and. trigg /= 2) stop "berry must be used with trigg = 2"

    iall=0
    xl=2.5e6
    c0=.002

    ! c0=.000   !para teste de conservacao (prec=0)
    ! cd= detrainment function 
    ! q = environmental q on model levels
    ! qe_cup = environmental q on model cloud levels
    ! qes_cup = saturation q on model cloud levels
    ! dby = buoancy term
    ! cd= detrainment function 
    ! zu = normalized updraft mass flux
    ! gamma_cup = gamma on model cloud levels
    ! mentr_rate = entrainment rate
    ! qc = cloud q (including liquid water) after entrainment
    ! qrch = saturation q in cloud
    ! qrc = liquid water content in cloud after rainout
    ! pw = condensate that will fall out at that level
    ! pwav = totan normalized integrated condensate (i1)
    ! c0 = conversion rate (cloud to rain)

    !--- no precip for small clouds
    if(iens.eq.shallow)  c0=0.

    !do k=1,mkx
    !     do i=istart,iend
    !        if(mentr_rate2(i,k).gt.0.)then
    !          radius = 0.2/(0.2/1200)
    !	 if(radius.lt.900.)
    !     !         if(radius.lt.900.)iall=0
    !         endif
    !     enddo
    !enddo


    do i=istart,iend
       pwav(i)=0.
    enddo
    do k=1,mkx
       do i=istart,iend
          pw(i,k) =0.
          !_srf     qc(i,k) =qes_cup(i,k)
          qc(i,k) =qe_cup(i,k)
          qrc(i,k)=0.
       enddo
    enddo
    do i=istart,iend
       if(ierr(i).eq.0.)then
          do k=k22(i),kbcon(i)-1
             qc(i,k)=qe_cup(i,k22(i))
          enddo
       endif
    enddo

    do k=2,mkx-1
       do i=istart,iend
          if (ierr(i).ne.0)  cycle
          if (k.lt.kbcon(i)) cycle
          if (k.gt.ktop(i))  cycle

          dz=z_cup(i,k)-z_cup(i,k-1)

          !--- 1. steady state plume equation, for what could
          !---    be in cloud without condensation

          qc(i,k)=(qc(i,k-1)*(1.-.5*cd(i,k)*dz)+mentr_rate2(i,k)*	     &
               dz*q(i,k-1))/(1.+mentr_rate2(i,k)*dz-.5*cd(i,k)*dz)

          !--- 3.condensation
          !   qtc  = total condensed water (var local)
          !--- saturation  in cloud, this is what is allowed to be in it

          qrch=qes_cup(i,k)+(1./xl)*(gamma_cup(i,k)/(1.+gamma_cup(i,k)))*dby(i,k)


          !------------------------- comeco do trecho modificado----------------------------
          if(autoconv.eq.2 .and. iens.eq.deep) then


             ! total condensed water (qtc)
             qtc(i,k) = max(0., qc(i,k)-qrch) ! kg[h2o]/kg[ar]= g[h2o]/g[ar]

             if(w_up(i,k) <= 0. .or. qtc(i,k) == 0. ) then
                con = 0.

             else

                !- rho = air density (g[ar]/cm^3)
                q1 = rho(i,k)*qtc(i,k)  ! g[h2o]/cm^3

                !- berry's formulation for autoconversion
                con = 1.e+6*q1*q1/(60.0*(5.0*rho(i,k) + 0.0366*ccn(i)/ &
                     ( 1.e6 * qtc(i,k) * bdispc)  ) ) ! kg[h2o]/ ( kg[ar] s)

                !- w_up is the vertical velocity of the updraft (m/s)
                con =      con/(0.75*min(10.,w_up(i,k))) ! kg[h2o]/ ( kg[ar] m)
             endif

             !- rain water production (pw)
             pw(i,k)=con*dz*zu(i,k)  !unit: kg[liq water]/kg[air]

             !- limit pw to the max allowed (the total condensed)
             pw(i,k)=min(pw(i,k),qtc(i,k))

             !- condensed water remained in cloud
             qrc(i,k) = qtc(i,k) - pw(i,k)

!!!!pw(i,k)=pw(i,k)*zu(i,k)  !unit: kg[liq water]/kg[air]
          else
             !--- 3.condensation
             !--- liquid water content in cloud after rainout

             qrc(i,k)=(qc(i,k)-qrch)/(1.+c0*dz*zu(i,k))
             if (qrc(i,k).lt.0.) then
                qrc(i,k)=0.
             endif

             pw(i,k)=c0*dz*qrc(i,k)*zu(i,k)  !unit: kg[liq water]/kg[air]
             !unit of c0 is m^-1
             !if (iall.eq.1) then
             !   qrc(i,k)=0.
             !   pw(i,k)=(qc(i,k)-qrch)*zu(i,k)
             !   if (pw(i,k).lt.0.) pw(i,k)=0.
             !endif
          endif
          !
          !print*,'cons=',i,k,j,pw(i,k)*1000.!,100.*(pw(i,k)/zu(i,k)+qrc(i,k)+qrch-qc(i,k))/(1.e-13+qc(i,k));call flush(6)


          !--- set next level for the in cloud total water
          qc(i,k)=qrc(i,k)+qrch
          !
          !--- integrated normalized ondensate
          !
          pwav(i)=pwav(i)+pw(i,k)

       enddo
    enddo
    return
  end subroutine cup_up_moisture

  !----------------------------------------------------------------------
  subroutine cup_up_nms(zu, z_cup, entr2, cd, kbcon, ktop,&
       mix, mgmxp, mkx, mgmzp, istart, iend, ierr, k22)
    real, intent(out) :: zu(:,:)
    real, intent(in) :: z_cup(:,:)
    real, intent(in) :: entr2(:,:)
    real, intent(in) :: cd(:,:)
    integer, intent(in) :: kbcon(:)
    integer, intent(in) :: ktop(:)
    integer, intent(inout) :: mix !**(JP)** unused
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mkx
    integer, intent(in) :: mgmzp
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(in) :: ierr(:)
    integer, intent(in) :: k22(:)

    integer i, k
    real dz



    do k=1,mkx
       do i=istart,iend
          zu(i,k)=0.
       enddo
    enddo
    do i=istart,iend
       if (ierr(i).eq.0) then
          do k=k22(i),kbcon(i)
             zu(i,k)=1.
          enddo
          do k=kbcon(i)+1,ktop(i)
             dz=z_cup(i,k)-z_cup(i,k-1)
             zu(i,k)=zu(i,k-1)*(1.+(entr2(i,k)-cd(i,k))*dz)
          enddo
       endif
    enddo
    return
  end subroutine cup_up_nms


  !----------------------------------------------------------------------
  subroutine cup_up_aa0(aa0, z, zu, dby, gamma_cup, t_cup,              &
       kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend, ierr)
    real, intent(out) :: aa0(:)
    real, intent(in) :: z(:,:)
    real, intent(in) :: zu(:,:)
    real, intent(in) :: dby(:,:)
    real, intent(in) :: gamma_cup(:,:)
    real, intent(in) :: t_cup(:,:)
    integer, intent(in) :: kbcon(:)
    integer, intent(in) :: ktop(:)
    integer, intent(in) :: mix
    integer, intent(in) :: mgmxp
    integer, intent(in) :: mkx
    integer, intent(in) :: mgmzp
    integer, intent(in) :: istart
    integer, intent(in) :: iend
    integer, intent(in) :: ierr(:)


    integer :: i
    integer :: k
    real :: dz
    real :: da


    do i=istart,iend
       aa0(i)=0.
    enddo
    do k=2,mkx-1
       !do 100 k=2,mkx-1
       do i=istart,iend
          !do 100 i=istart,iend
          !if (ierr(i).ne.0)  go to 100
          if (ierr(i).ne.0)  cycle
          !if (k.le.kbcon(i)) go to 100
          if (k.le.kbcon(i)) cycle
          !if (k.gt.ktop(i))  go to 100
          if (k.gt.ktop(i))  cycle
          dz = z(i,k)-z(i,k-1)
          da = zu(i,k)*dz*(9.81/(1004.*((t_cup(i,k)))))*dby(i,k-1)/  &
               (1.+gamma_cup(i,k))
          !if (k.eq.ktop(i).and.da.le.0.) go to 100
          if (k.eq.ktop(i).and.da.le.0.) cycle
          aa0(i)=aa0(i)+da
          if (aa0(i).lt.0.) aa0(i)=0.

          !100     continue
       enddo
    enddo
  end subroutine cup_up_aa0
end module ModCupUp
