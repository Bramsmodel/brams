module ModRrtmDriver

  use iso_fortran_env, only: &
       int64, &
       real64

  use ccatt_start, only: &
       ccatt

  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModNamelistFile, only: &
       NamelistFile
  
  use ModBasicFields, only: &
       BasicFields

  use ModMicControl, only: &
       MicControl
  
  use ModRadiateFields, only: &
       RadiateFields
  
  use mem_grid    , only: &
       ngrid, &
       time, &
       dtlt, &
       itime1, &
       nzg, &
       nzs, &
       npatch, &
       grid_g, &
       nnzp, &
       if_adap, &
       zm, &
       zt, &
       naddsc, &
       imonth1, &
       idate1, &
       iyear1, &
       centlat, &
       centlon, &
       ztop, &
       dzm , &
       ngrids, &
       grid_vars, &
       dzt, &
       deltaxn, &
       deltayn, &
       zmn

  use ModOptical, only: &
       optProp

  use mem_carma, only: &
       carma, &
       carma_aotMap

  use node_mod, only: &
       mynum ! INTENT(IN)

  use mem_rrtm, only: &
       aot_rrtm_lw, &
       aot_rrtm_sw, &
       sig, &
       delsig, &
       sigmid, &
       adjes, &
       adjust, &
       ch4ispresent, &
       ch4pos, &
       co2ispresent, &
       firsttime, &
       iceflglw, &
       iceflgsw, &
       icld, &
       idrv, &
       inflglw, &
       liqflglw, &
       nls, &
       no2ispresent, &
       no2pos, &
       o2ispresent, &
       o2pos, &
       co2pos, &
       inflgsw, &
       o3ispresent, &
       o3pos, &
       pi, &
       rch4ar, &
       rco2ar, &
       rh2oar, &
       rno2ar, &
       ro2ar, &
       ro3ar, &
       scon, &
       bndmax_sw, &
       bndmax_lw, &
       liqflgsw, &
       bndmin_sw, &
       bndmin_lw, &
       initRRTM, &
       flip, &
       irng, &
       permuteseed

  use mem_tend, only: &
       tend

  use rrtmg_sw_rad, only: &
       rrtmg_sw

  use rrtmg_lw_rad, only: &
       rrtmg_lw

  use mcica_subcol_gen_sw, only: &
       mcica_subcol_sw

  use mcica_subcol_gen_lw, only: &
       mcica_subcol_lw

  use grid_dims, only: &
       nzpmax

  use ModCuParmFields, only: &
       CuParmFields

  use ModMicroFields, only: &
       MicroFields

  use parkind, only : &
       im => kind_im, &
       rb => kind_rb

  use rconstants, only : &
       cp, &
       cpor, &
       p00, &
       stefan, &
       cpi, &
       pi180

  use parrrsw, only : &
       nbndsw, &
       ngptsw, &
       naerec, &
       jpband

  use mem_leaf, only: &
       leaf_g, &
       isfcl, &
       albedo ! intent(in)

  use moddateutils, only: &
       julday

  use mem_chem1, only: &
       chem1_g, &
       chemistry

  use ref_sounding, only: &
       pi01dn ! intent(in)

  use parrrtm, only : &
       nbndlw, &
       ngptlw, &
       maxxsec, &
       mxmol

  use rrtmg_lw_cldprop, only: &
       cldprop

  use rrtmg_sw_cldprop, only: &
       cldprop_sw

  use mem_tuv, only: &
       carma_tuv, &
       tuv2carma

  use mem_scratch1_grell, only: &
       xmb4d, &
       zup5d, &
       clwup5d, &
       up_massdetr5d

  use ModLeafComs, only: &
       rlonga_a, &
       rlonga_gs, &
       rlonggs_v, &
       rlonga_v, &
       rlongv_a, &
       rlongv_gs, &
       rlonggs_a, &
       rshort_g, &
       rshort_v, &
       snowfac, &
       vf, &
       tempk, &
       slcpd, &
       fracliq, &
       emisv, &
       emisg, &
       slmsts, &
       rshort_s

  implicit none

  private
  public :: rrtm_driver

  real(kind=real64), allocatable, dimension(:,:,:) :: ozone
  real(kind=real64) :: ozsig(18)
  integer, parameter :: nlm_getoz=18
  logical :: first_getoz,inter_getoz
  integer :: year_getoz,mon_getoz
  integer, parameter :: nl=37
  integer, parameter :: ns=4
  include "aerosol_setup.f90"

contains


  subroutine rrtm_driver(mzp, mxp, myp, ia, iz, ja, jz, mynum, &
       oneNamelistFile, oneBasicFields, oneMicVars, oneMicroFields, &
       oneRadiateFields, oneCuParmFields)

    integer, intent(in) :: mzp, mxp, myp, ia, iz, ja, jz, mynum
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicControl), pointer, intent(in) :: oneMicVars
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields
    type(CuParmFields), pointer, intent(in) :: oneCuParmFields

    real,dimension(mzp,mxp,myp) :: lwl,iwl
    real,dimension(mxp,myp) :: rain
    real :: patch_area,can_temp,veg_temp,leaf_class,veg_height
    real :: veg_fracarea,veg_albedo,sfcwater_nlev,rshort,rlong
    real :: albedt,rlongup,cosz,hrAngleLocal,cdec
    real, dimension(nzg) :: soil_energy,soil_water,soil_text
    real, dimension(nzs) :: sfcwater_energy,sfcwater_depth
    real :: maxCloud_fraction
    real :: solfac
    integer :: ip,i,j,jday
    integer :: icount = 0
    integer, parameter :: ngpt = 141

    character(len=16) :: str(10)
    character(len=*), parameter :: h="**(rrtm_driver)**"
    logical, parameter :: dumpLocal=.false.
    
    !- if not including radiation, return
    if ((oneNamelistFile%ilwrtyp + oneNamelistFile%iswrtyp)==0) return

    icount = icount + 1
    if (icount>ngpt) icount = 0

    !--- check if it is time to recompute radiative tendency and fluxes
    !
    !--- radiation calculation is updated only every radfrq seconds
    if ( (mod(time+.001, oneNamelistFile%radfrq) < dtlt .or. time<0.001)) then

       !--- set radiation tendency for theta to zero
       oneRadiateFields%fthrd(1:mzp,1:mxp,1:myp) = 0.0

       !- call routine to calculate cloud properties for RRTM/CARMA
       !--1st, set zero to the local arrays
       oneRadiateFields%cloud_fraction=0.0;  rain=0.0; lwl=0.0; iwl =0.0

       call  cloud_prop_rrtm(mzp, mxp, myp, ia, iz, ja, jz &
                                !-- output
            ,oneRadiateFields%cloud_fraction  &
            ,rain  &
            ,lwl  &
            ,iwl,  &
            oneNamelistFile, oneBasicFields, oneMicVars, oneMicroFields, &
            oneCuParmFields)

       !-srf tuning section for cloud fraction and other parameters for radiation
       if(oneNamelistFile%radtun /= 1.0) then
          oneRadiateFields%cloud_fraction=  min(1.,oneNamelistFile%radtun* oneRadiateFields%cloud_fraction)
          rain          =  oneNamelistFile%radtun*rain
          lwl           =  oneNamelistFile%radtun*lwl
          iwl           =  oneNamelistFile%radtun*iwl
       endif



       !- Compute solar zenith angle [cosz(i,j)] & solar constant factr [solfac].
       call zen_rtm(imonth1, idate1, iyear1, time, itime1, centlat, centlon, &
            oneNamelistFile%lonrad, pi180, ia, iz, ja, jz, jday, solfac, hranglelocal, cdec,&
            mynum, oneRadiateFields)

       !- compute patch-averaged surface albedo [albedt(i,j)] and up longwave
       !- radiative flux [rlongup(i,j)].
       !DSM ---- In case of JULES surface scheme, rlongup and albedt are already provided
       if (isfcl /= 5 .or. time<0.001) then
          oneRadiateFields%albedt  = 0.
          oneRadiateFields%rlongup = 0.


          if (dumpLocal) then
             call MsgDump(h//" copia soil_energy do leaf_g para var local")
             call MsgDump(h//" invoca sfcrad_rtm com a var local")
          end if
       
          do ip = 1,npatch
             do j = 1,jz
                do i = 1,iz

                   !By now just copy the memory to local variables
                   !in order to send to sfcrad
                   soil_energy    (1:nzg)=leaf_g(ngrid)%soil_energy(1:nzg,i,j,ip)
                   soil_water     (1:nzg)=leaf_g(ngrid)%soil_water (1:nzg,i,j,ip)
                   soil_text      (1:nzg)=leaf_g(ngrid)%soil_text  (1:nzg,i,j,ip)
                   sfcwater_energy(1:nzs)=leaf_g(ngrid)%sfcwater_energy(1:nzs,i,j,ip)
                   sfcwater_depth (1:nzs)=leaf_g(ngrid)%sfcwater_depth (1:nzs,i,j,ip)
                   patch_area            =leaf_g(ngrid)%patch_area   (i,j,ip)
                   can_temp              =leaf_g(ngrid)%can_temp     (i,j,ip)
                   veg_temp              =leaf_g(ngrid)%veg_temp     (i,j,ip)
                   leaf_class            =leaf_g(ngrid)%leaf_class   (i,j,ip)
                   veg_height            =leaf_g(ngrid)%veg_height   (i,j,ip)
                   veg_fracarea          =leaf_g(ngrid)%veg_fracarea (i,j,ip)
                   veg_albedo            =leaf_g(ngrid)%veg_albedo   (i,j,ip)
                   sfcwater_nlev         =leaf_g(ngrid)%sfcwater_nlev(i,j,ip)

                   rshort                =oneRadiateFields%rshort (i,j)
                   rlong                 =oneRadiateFields%rlong  (i,j)
                   albedt                =oneRadiateFields%albedt (i,j)
                   rlongup               =oneRadiateFields%rlongup(i,j)
                   cosz                  =oneRadiateFields%cosz   (i,j)

!!$                   if (patch_area < 0.0) then
!!$                      write(str(1),"(i3)") i
!!$                      write(str(2),"(i3)") j
!!$                      write(str(3),"(i3)") ip
!!$                      write(str(4),"(i3)") ngrid
!!$                      write(str(5),"(e15.7)") patch_area
!!$                      call MsgDump(h//" patch_area("//trim(adjustl(str(1)))//&
!!$                           ","//trim(adjustl(str(2)))//&
!!$                           ","//trim(adjustl(str(3)))//&
!!$                           ")="//trim(adjustl(str(5)))//&
!!$                           " para ngrid="//trim(adjustl(str(4)))//&
!!$                           " na linha 343")
!!$                      call fatal_error(h//" line 343")
!!$                   end if
                   !
                   call sfcrad_rtm(nzg, nzs, ip,             &
                        soil_energy,    soil_water,     &
                        soil_text  ,    sfcwater_energy,&
                        sfcwater_depth, patch_area,     &
                        can_temp,       veg_temp,       &
                        leaf_class,     veg_height,     &
                        veg_fracarea,   veg_albedo,     &
                        sfcwater_nlev,                  &
                        rshort, rlong,  albedt,         &
                        rlongup, cosz                   &
                        )

                   !Copy back albedo and rlongUP to memory
                   oneRadiateFields%albedt (i,j)=albedt
                   oneRadiateFields%rlongup(i,j)=rlongup

!!$                   if (rlongup < 0.0) then
!!$                      write(str(1),"(i3)") i
!!$                      write(str(2),"(i3)") j
!!$                      write(str(3),"(e15.7)") rlongup
!!$                      call MsgDump(h//" rlongup("//trim(adjustl(str(1)))//&
!!$                           ","//trim(adjustl(str(2)))//&
!!$                           "="//trim(adjustl(str(3))))
!!$                      call fatal_error(h//" line 357")
!!$                   end if
                end do
             end do
          end do
       endif



       !- RRTM Radiation
       call radrrtmdrv(ia,iz,ja,jz,mxp,myp,mzp,mynum&
            , oneRadiateFields%cloud_fraction         &
            ,rain                   &
            ,lwl                    &
            ,iwl                    &
            ,icount                 &
            ,ngpt,                  &
            oneBasicFields, oneMicroFields, oneRadiateFields)
    endif
    !--- apply radiative tendencies to model tendencies
    call tend_accum_rtm(mzp, mxp, myp, ia, iz, ja, jz, oneRadiateFields)

  end subroutine rrtm_driver

  ! ****************************************************************************

  subroutine tend_accum_rtm(m1,m2,m3,ia,iz,ja,jz, oneRadiateFields)
    integer, intent(in) :: m1, m2, m3, ia, iz, ja, jz
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields

    ! local variables:
    integer :: i, j, k,ipos

    ipos=0
    do j=1,m3
       do i=1,m2
          do k=1,m1
             ipos=ipos+1
             tend%tht(ipos) = tend%tht(ipos) + oneRadiateFields%fthrd(k,i,j)
          end do
       end do
    end do

  end subroutine tend_accum_rtm
  ! ****************************************************************************

  subroutine zen_rtm(imonth1, idate1, iyear1, time, itime1, centlat, centlon, &
       lonrad, pi180, ia, iz, ja, jz, jday, solfac, hrangle, cdec, mynum, &
       oneRadiateFields)

    ! arguments:
    integer, intent(in)  :: imonth1, idate1, iyear1, itime1,mynum
    real, intent(in)     :: time
    real, intent(in)     :: centlat(:), centlon(:)
    integer, intent(in)  :: lonrad
    real, intent(in)     :: pi180
    integer, intent(in)  :: ia, iz, ja, jz
    integer, intent(out) :: jday
    real, intent(out)    :: solfac
    real, intent(out)    :: hrangle
    real, intent(out)    :: cdec
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields
    ! local variables:
    integer :: i, j
    real    :: sdec, declin, d0, d02, dayhr, radlat, cslcsd, snlsnd, gglon, &
         dayhrr, eqt

    jday   = julday(imonth1, idate1, iyear1)
    jday   = jday + nint(time/86400.)
    !      sdec - sine of declination, cdec - cosine of declination
    declin = -23.5*cos(6.283/365.*(jday + 9))*pi180
    sdec   = sin(declin)
    cdec   = cos(declin)

    ! find the factor, solfac, to multiply the solar constant to correct
    ! for earth's varying distance to the sun.

    d0     = 6.2831853*float(jday-1)/365.
    d02    = d0*2.
    solfac = 1.000110 + 0.034221*cos(d0) + 0.001280*sin(d0) + &
         0.000719*cos(d02) + 0.000077*sin(d02)

    ! find the hour angle, then get cosine of zenith angle.

    !ner_i - including solar time equation ("eqt" must be defined, it is a new variable)
    eqt = (0.000075 + 0.001868*cos(d0) - 0.032077*sin(d0) - 0.014615*cos(d02) &
         - 0.040849*sin(d02))*1440/(2*3.141593)
    !ner_f - including solar time equation

    !-ner dayhr  = time/3600. + float(itime1/100) + float(mod(itime1,100))/60.
    dayhr = (time / 3600. + float(itime1/100) + float(mod(itime1,100)) / 60.) !&

    do j = ja,jz
       do i = ia,iz
          radlat = grid_g(ngrid)%glat(i,j)*pi180
          if (lonrad==0)      radlat = centlat(1)*pi180
          if (radlat==declin) radlat = radlat + 1.e-5
          cslcsd = cos(radlat)*cdec
          snlsnd = sin(radlat)*sdec
          gglon  = grid_g(ngrid)%glon(i,j)
          if (lonrad==0)      gglon = centlon(1)

          !ner_i new hour angle calculation
          hrangle=((dayhr+(gglon/15.)-(12.-eqt/60.))*15./1.)*3.141593/180.
          !          dayhrr    = mod(dayhr+gglon/15.+24., 24.)
          !          hrangl    = 15.*(dayhrr - 12.)*pi180

          oneRadiateFields%cosz(i,j) = snlsnd + cslcsd*cos(hrangle)

          !oneRadiateFields%cosz(i,j) = min(oneRadiateFields%cosz(i,j), 1.0)
          !oneRadiateFields%cosz(i,j) = max(oneRadiateFields%cosz(i,j),-1.0)

          oneRadiateFields%cosz(i,j) = max(oneRadiateFields%cosz(i,j),1.0e-10)

       end do
    end do
  end subroutine zen_rtm
  ! ****************************************************************************

  subroutine sfcrad_rtm(nzg, nzs, ip,                                            &
       soil_energy, soil_water, soil_text, sfcwater_energy, sfcwater_depth,  &
       patch_area, can_temp, veg_temp, leaf_class, veg_height, veg_fracarea, &
       veg_albedo, sfcwater_nlev, rshort, rlong, albedt, rlongup, cosz      &
       )

    ! Arguments:
    integer, intent(IN) :: nzg, nzs, ip
    real, intent(IN)    :: soil_energy(nzg)
    real, intent(IN)    :: soil_water(nzg)
    real, intent(IN)    :: soil_text(nzg)
    real, intent(IN)    :: sfcwater_energy(nzs)
    real, intent(IN)    :: sfcwater_depth(nzs)
    real, intent(IN)    :: patch_area
    real, intent(IN)    :: can_temp
    real, intent(IN)    :: veg_temp
    real, intent(IN)    :: leaf_class
    real, intent(IN)    :: veg_height
    real, intent(IN)    :: veg_fracarea
    real, intent(IN)    :: veg_albedo
    real, intent(IN)    :: sfcwater_nlev
    real, intent(IN)    :: rshort
    real, intent(IN)    :: rlong
    real, intent(INOUT) :: albedt
    real, intent(INOUT) :: rlongup
    real, intent(IN)    :: cosz

    ! Local Variables:
    integer :: k, nsoil, nveg, ksn
    real :: alb, vfc, fcpct, alg, rad, als, fractrans, absg, algs, emv, emgs, &
         gslong, vlong, alv
    real :: salvo
    
    real :: vctr32(nint(sfcwater_nlev)+10)

    character(len=16) :: str(10)
    character(len=*), parameter :: h="**(sfcrad_rtm)**"
    logical, parameter :: dumpLocal=.false.
    
    ! This routine is called by the radiation parameterization and by leaf.
    ! It computes net surface albedo plus radiative exchange between the
    ! atmosphere, vegetation, and the snow/ground given previously computed
    ! downward longwave and shortwave fluxes from the atmosphere.
    ! Also computed are functions of snowcover that are required for the above
    ! radiation calculations as well as other calculations in leaf.

    ! The shortwave parameterizations are only valid if the cosine of the zenith
    ! angle is greater than .03 .  Water albedo from Atwater and Bell (1981)

    ! alg, als, and alv are the albedos of the ground, snow, and vegetation
    ! (als needs a better formula based on age of the surface snow).

    ! absg and vctr32 are the actual fractions of shortwave incident on snow
    ! plus ground that get absorbed by the ground and each snow layer,
    ! respectively.  They currently use the variable fractrans, which is the
    ! fraction of light transmitted through each layer based on mass per square
    ! meter.  algs is the resultant albedo from snow plus ground.

    if (ip==1) then

!!$       if (rlongup < 0.0) then
!!$          write(str(3),"(e15.7)") rlongup
!!$          call MsgDump(h//" rlongup"//&
!!$               "="//trim(adjustl(str(3))))
!!$          call fatal_error(h//" na entrada, linha 545")
!!$       end if

       if (cosz>0.03) then
          alb    = min(max(-.0139 + .0467*tan(acos(cosz)), 0.03), 0.999)
          albedt = albedt + patch_area*alb
       endif

       call qtk(soil_energy(nzg), tempk(nzg), fracliq(nzg))
       salvo=rlongup
       rlongup = rlongup + patch_area*stefan*tempk(nzg)**4

!!$       if (rlongup < 0.0) then
!!$          write(str(1),"(e15.7)") salvo
!!$          write(str(2),"(e15.7)") patch_area
!!$          write(str(3),"(e15.7)") stefan
!!$          write(str(4),"(e15.7)") tempk(nzg)**4
!!$          write(str(5),"(e15.7)") rlongup
!!$          call MsgDump(h//" na linha 556, "//&
!!$               trim(adjustl(str(5)))//" = "//trim(adjustl(str(1)))//&
!!$               " + "//trim(adjustl(str(2)))//"*"//trim(adjustl(str(3)))//&
!!$               "*"//trim(adjustl(str(4))))
!!$          call fatal_error(h//" line 556")
!!$       end if
       
    elseif (isfcl==0) then

       albedt  = albedt + patch_area*albedo
       rlongup = rlongup + patch_area*stefan*can_temp**4
!!$       if (rlongup < 0.0) then
!!$          write(str(3),"(e15.7)") rlongup
!!$          call MsgDump(h//" rlongup"//&
!!$               "="//trim(adjustl(str(3))))
!!$          call fatal_error(h//" line 553")
!!$       end if

    else

       ! Diagnose soil temperature and liquid fraction


       do k=1,nzg
          nsoil = nint(soil_text(k))
          call qwtk(soil_energy(k), soil_water(k)*1.e3,  &
               slcpd(nsoil), tempk(k), fracliq(k))
       enddo

       ! Diagnose snow temperature and the influence of snow covering veg.

       nveg    = nint(leaf_class)
       ksn     = nint(sfcwater_nlev)
       snowfac = 0.
       do k=1,ksn
          snowfac = snowfac + sfcwater_depth(k)
          call qtk(sfcwater_energy(k), tempk(k+nzg), fracliq(k+nzg))
       enddo
       snowfac = min(0.99, snowfac/max(0.001, veg_height))

       vf  = veg_fracarea*(1. - snowfac)
       vfc = 1. - vf

       ! Shortwave radiation calculations

       !srf-25-02-2005
       nsoil = nint(soil_text(nzg))
       !srf

       fcpct = soil_water(nzg)/slmsts(nsoil)
       if (fcpct>0.5) then
          alg = 0.14
       else
          alg = 0.31 - 0.34*fcpct
       endif
       alv = veg_albedo

       rad = 1.
       if (ksn>0) then

          ! als = .14 (the wet soil value) for all-liquid

          als = 0.5 - 0.36*fracliq(ksn + nzg)
          rad = 1.  - als
       endif
       do k=ksn,1,-1
          fractrans = exp(-20.*sfcwater_depth(k))
          vctr32(k) = rad*(1. - fractrans)
          rad = rad * fractrans
       enddo
       absg = (1. - alg)*rad
       algs = 1. - absg
       do k=ksn,1,-1
          algs        = algs - vctr32(k)
          rshort_s(k) = rshort*vfc*vctr32(k)
       enddo
       rshort_g = rshort*vfc*absg
       rshort_v = rshort*vf*(1. - alv + vfc*algs)
       !  rshort_a = rshort*(vf*alv + vfc*vfc*algs)

       alb = vf*alv + vfc*vfc*algs
       albedt = albedt + patch_area*alb

       ! Longwave radiation calculations

       emv  = emisv(nveg)
       emgs = emisg(nsoil)
       if (ksn>0) emgs = 1.0
       gslong = emgs*stefan*tempk(ksn+nzg)**4
       vlong  = emv*stefan*veg_temp**4

       rlonga_v  = rlong*vf*(emv + vfc*(1. - emgs))
       rlonga_gs = rlong*vfc*emgs
       rlongv_gs = vlong*vf*emgs
       rlongv_a  = vlong*vf*(2. - emgs - vf + emgs*vf)
       rlonggs_v = gslong*vf*emv
       rlonggs_a = gslong*vfc
       rlonga_a  = rlong*(vf*(1. - emv) + vfc*vfc*(1. - emgs))
       rlongup = rlongup + patch_area*(rlongv_a + rlonggs_a + rlonga_a)

       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ! In case rlong is not computed, zero out all longwave fluxes other
       ! than rlongup.  [On the first timestep, radiative fluxes may not
       ! be available until microphysics is called, and zeroing these fluxes
       ! avoids the imbalance of having upward longwave without downward longwave
       ! in LEAF-2.  Also, this allows LEAF-2 to run without radiation for all
       ! timesteps, if desired for testing.]

       if (rlong<0.1) then
          rlonga_v  = 0.
          rlonga_gs = 0.
          rlongv_gs = 0.
          rlongv_a  = 0.
          rlonggs_v = 0.
          rlonggs_a = 0.
          rlonga_a  = 0.
       endif

    endif

  end subroutine sfcrad_rtm

  ! ****************************************************************************
  subroutine radrrtmdrv(ia,iz,ja,jz,mxp,myp,mzp,mynum &
       ,cloud_fraction           &
       ,rain                     &
       ,lwl                      &
       ,iwl                      &
       ,icount                   &
       ,ngpt,                    &
       oneBasicFields, oneMicroFields, oneRadiateFields)
    integer, intent(in) :: ia,iz,ja,jz
    integer, intent(in) :: icount
    integer, intent(in) :: ngpt

    integer, intent(in) :: mxp,myp,mzp,mynum
    real, intent(in), dimension(mzp,mxp,myp) :: cloud_fraction !cloud_fraction
    real, intent(in), dimension(    mxp,myp) :: rain !total rain water
    real, intent(in), dimension(mzp,mxp,myp) :: lwl !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
    real, intent(in), dimension(mzp,mxp,myp) :: iwl !total cloud ice water (kg/kg for carma and g/m2 for rrtm)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields

    integer,allocatable,dimension(:) :: ipos,jpos
    integer (kind=int64),allocatable,dimension(:)  :: imask
    logical, allocatable, dimension(:) :: aboveh

    integer :: nlay
    integer, parameter :: iplon=1

    real(kind=rb) :: taucloud_lw(mzp-1,nbndlw)

    real(kind=rb) :: taucloud_sw(mzp-1,jpband)
    real(kind=rb) :: tauc_lw(nbndlw,mzp-1)
    real(kind=rb) :: ssac_lw(nbndlw,mzp-1)
    real(kind=rb) :: asmc_lw(nbndlw,mzp-1)
    real(kind=rb) :: fsfc_lw(nbndlw,mzp-1)

    real(kind=rb) :: tauc_sw(nbndsw,mzp-1)
    real(kind=rb) :: ssac_sw(nbndsw,mzp-1)
    real(kind=rb) :: asmc_sw(nbndsw,mzp-1)
    real(kind=rb) :: fsfc_sw(nbndsw,mzp-1)



    real(kind=rb) :: taucldorig(mzp-1,jpband)

    real(kind=rb) :: ssacloud_sw(mzp-1,jpband)
    real(kind=rb) :: ssacloud_lw(mzp-1,jpband)

    real(kind=rb) :: asmcloud_sw(mzp-1,jpband)
    real(kind=rb) :: asmcloud_lw(mzp-1,jpband)

    real(kind=rb),parameter :: factorcor=1.0e-9_rb

    real(kind=rb), allocatable :: play(:,:)          ! layer pressures (hpa, mb)
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: plev(:,:)          ! interface pressures (hpa, mb)
    !    dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: tlay(:,:)          ! layer temperatures (k)
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: tlev(:,:)          ! interface temperatures (k)
    !    dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: tsfc(:)            ! surface temperature (k)
    !    dimensions: (ncol)
    real(kind=rb), allocatable :: h2ovmr(:,:)        ! h2o volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: o3vmr(:,:)         ! o3 volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: co2vmr(:,:)        ! co2 volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: ch4vmr(:,:)        ! methane volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: n2ovmr(:,:)        ! nitrous oxide volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: o2vmr(:,:)         ! oxygen volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: cfc11vmr(:,:)      ! cfc11 volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: cfc12vmr(:,:)      ! cfc12 volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: cfc22vmr(:,:)      ! cfc22 volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: ccl4vmr(:,:)       ! ccl4 volume mixing ratio
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: emis(:,:)          ! surface emissivity
    !    dimensions: (ncol,nbndlw)
    real(kind=rb), allocatable :: asdir(:)           ! uv/vis surface albedo direct rad
    !    dimensions: (ncol)
    real(kind=rb), allocatable :: aldir(:)           ! near-ir surface albedo direct rad
    !    dimensions: (ncol)
    real(kind=rb), allocatable :: asdif(:)           ! uv/vis surface albedo: diffuse rad
    !    dimensions: (ncol)
    real(kind=rb), allocatable :: aldif(:)           ! near-ir surface albedo: diffuse rad
    !    dimensions: (ncol)

    integer(kind=im) :: dyofyr          ! day of the year (used to get earth/sun
    !  distance if adjflx not provided)
    real(kind=rb), allocatable :: coszen(:)         ! cosine of solar zenith angle
    !    dimensions: (ncol)

    !    in upward flux as a function of
    !    surface temperature [0=off, 1=on]
    !    0: normal forward calculation
    !    1: normal forward calculation with
    !       duflx_dt and duflxc_dt output

!!!Cloud properties from mcica_sw with dimensions: (ngptsw,ncol,nlay)

    !LFR For version 5.0 of rrtmg
    real(kind=rb), allocatable :: alpha_sw(:,:) ! cloud fraction decorrelation length (m)
    real(kind=rb), allocatable :: cldfmcl_sw(:,:,:)    ! cloud fraction from mcica_sw
    real(kind=rb), allocatable :: taucmcl_sw(:,:,:)    ! in-cloud optical depth 
    real(kind=rb), allocatable :: ssacmcl_sw(:,:,:)    ! in-cloud single scattering albedo 
    real(kind=rb), allocatable :: asmcmcl_sw(:,:,:)    ! in-cloud asymmetry parameter 
    real(kind=rb), allocatable :: fsfcmcl_sw(:,:,:)    ! in-cloud forward scattering fraction 
    real(kind=rb), allocatable :: ciwpmcl_sw(:,:,:)    ! in-cloud ice water path (g/m2) 
    real(kind=rb), allocatable :: clwpmcl_sw(:,:,:)    ! in-cloud liquid water path (g/m2) 
    real(kind=rb), allocatable :: reicmcl_sw(:,:)      ! cloud ice effective radius (microns) 
    real(kind=rb), allocatable :: relqmcl_sw(:,:)      ! cloud water drop effective radius (microns)

!!!Cloud properties from mcica_lw with dimensions: (ngptlw,ncol,nlay)

    real(kind=rb), allocatable :: cldfmcl_lw(:,:,:)    ! cloud fraction from mcica_lw
    real(kind=rb), allocatable :: taucmcl_lw(:,:,:)    ! in-cloud optical depth
    real(kind=rb), allocatable :: ssacmcl_lw(:,:,:)    ! in-cloud single scattering albedo
    real(kind=rb), allocatable :: asmcmcl_lw(:,:,:)    ! in-cloud asymmetry parameter
    real(kind=rb), allocatable :: fsfcmcl_lw(:,:,:)    ! in-cloud forward scattering fraction
    real(kind=rb), allocatable :: ciwpmcl_lw(:,:,:)    ! in-cloud ice water path (g/m2)
    real(kind=rb), allocatable :: clwpmcl_lw(:,:,:)    ! in-cloud liquid water path (g/m2)

!!!No subcolumns for these properties yet
    real(kind=rb), allocatable :: reicmcl_lw(:,:)      ! cloud ice effective radius (microns) 
    real(kind=rb), allocatable :: relqmcl_lw(:,:)      ! cloud water drop effective radius (microns)

!!!Cloud properties rrtmg with dimension dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: cldf(:,:)            ! cloud fraction
    real(kind=rb), allocatable :: ciwp(:,:)            ! in-cloud ice water path (g/m2)
    real(kind=rb), allocatable :: clwp(:,:)            ! in-cloud liquid water path (g/m2)
    real(kind=rb), allocatable :: reic(:,:)            ! cloud ice effective radius (microns)
    ! specific definition of reicmcl depends on setting of iceflgsw:
    ! iceflgsw = 0: (inactive)
    !
    ! iceflgsw = 1: ice effective radius, r_ec, (ebert and curry, 1992),
    !               r_ec range is limited to 13.0 to 130.0 microns
    ! iceflgsw = 2: ice effective radius, r_k, (key, streamer ref. manual, 1996)
    !               r_k range is limited to 5.0 to 131.0 microns
    ! iceflgsw = 3: generalized effective size, dge, (fu, 1996),
    !               dge range is limited to 5.0 to 140.0 microns
    !               [dge = 0.9021 * r_ec - 7.0365] 
    real(kind=rb), allocatable :: relq(:,:)            ! cloud water drop effective radius (microns)

!!!Cloud properties rrtmg with dimensions : (nbndsw,ncol,nlay)

    real(kind=rb), allocatable :: taucl_sw(:,:,:)     ! in-cloud optical depth
    real(kind=rb), allocatable :: ssacl_sw(:,:,:)    ! in-cloud single scattering albedo
    real(kind=rb), allocatable :: asmcl_sw(:,:,:)    ! in-cloud asymmetry parameter
    real(kind=rb), allocatable :: fsfcl_sw(:,:,:)    ! in-cloud forward scattering fraction

!!!Cloud properties rrtmg with dimensions: (nbndlw,ncol,nlay)
    !LFR For version 5.0 of rrtmg
    real(kind=rb), allocatable :: alpha_lw(:,:) ! cloud fraction decorrelation length (m)

    real(kind=rb), allocatable :: taucl_lw(:,:,:)    ! in-cloud optical depth
    real(kind=rb), allocatable :: ssacl_lw(:,:,:)    ! in-cloud single scattering albedo
    real(kind=rb), allocatable :: asmcl_lw(:,:,:)    ! in-cloud asymmetry parameter
    real(kind=rb), allocatable :: fsfcl_lw(:,:,:)    ! in-cloud forward scattering fraction

!!!Aerosol properties rrtmg_sw with dimensions: (ncol,nlay,nbndsw)

    real(kind=rb), allocatable :: tauaer_sw(:,:,:)     ! aerosol optical depth (iaer=10 only) (non-delta scaled)
    real(kind=rb), allocatable :: ssaaer_sw(:,:,:)     ! aerosol single scattering albedo (iaer=10 only) (non-delta scaled)
    real(kind=rb), allocatable :: asmaer_sw(:,:,:)     ! aerosol asymmetry parameter (iaer=10 only) (non-delta scaled)


!!!Aerosol properties rrtmg_sw with dimensions: (ncol,nlay,nbndsw)

    real(kind=rb), allocatable :: tauaer_lw(:,:,:)     ! aerosol optical depth (iaer=10 only) (non-delta scaled)
    real(kind=rb), allocatable :: ssaaer_lw(:,:,:)     ! aerosol single scattering albedo (iaer=10 only) (non-delta scaled)
    real(kind=rb), allocatable :: asmaer_lw(:,:,:)     ! aerosol asymmetry parameter (iaer=10 only) (non-delta scaled)

!!!Aerosol optical depth at 0.55 micron
    real(kind=rb), allocatable :: ecaer    (:,:,:)     ! aerosol optical depth at 0.55 micron (iaer=6 only) (non-delta scaled)

    ! -- output -----

    real(kind=rb), allocatable :: swuflx(:,:)       ! total sky shortwave upward flux (w/m2)
    !    dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: swdflx(:,:)       ! total sky shortwave downward flux (w/m2)
    !    dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: swhr(:,:)         ! total sky shortwave radiative heating rate (k/d)
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: swuflxc(:,:)      ! clear sky shortwave upward flux (w/m2)
    !    dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: swdflxc(:,:)      ! clear sky shortwave downward flux (w/m2)
    !    dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: swhrc(:,:)        ! clear sky shortwave radiative heating rate (k/d)
    !    dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: uflx(:,:)         ! total sky longwave upward flux (w/m2)
    !        dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: dflx(:,:)         ! total sky longwave downward flux (w/m2)
    !        dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: hr(:,:)           ! total sky longwave radiative heating rate (k/d)
    !        dimensions: (ncol,nlay)
    real(kind=rb), allocatable :: uflxc(:,:)        ! clear sky longwave upward flux (w/m2)
    !        dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: dflxc(:,:)        ! clear sky longwave downward flux (w/m2)
    !        dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: hrc(:,:)           ! clear sky longwave radiative heating rate (k/d)
    !        dimensions: (ncol,nlay)

    real(kind=rb), allocatable :: duflx_dt(:,:)
    ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    dimensions: (ncol,nlay+1)
    real(kind=rb), allocatable :: duflxc_dt(:,:)
    ! change in clear sky upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !        dimensions: (ncol,nlay+1)
    integer :: nofcols,noc,np
    integer :: ilay,iSpc
    integer :: changeseed

    integer :: i,j,k,m
    real(kind=rb) :: picpi
    real(kind=real64) :: co2val,doy
    real(kind=rb),allocatable,dimension(:,:) :: p,t,q,fice,o3l,lmixr
    real(kind=rb),allocatable,dimension(:)   :: psurf,colrad
    character(len=2) :: cmzp
    integer :: itime(4)
    integer(kind=im) :: ncbandssw,ncbandslw

    character(len=16) :: str(10)
    character(len=*), parameter :: h="**(radrrtmdrv)**"
    
    integer(kind=im) :: iaer !LFR for the 5.0 version
    integer(kind=im), parameter :: isolvar=-1         ! Flag for solar variability method
    !#Solar variability scaling factors or indices (ISOLVAR=-1,1,2,3 only)
    !#            For ISOLVAR = 1:
    !#                    SOLVAR(1)    Facular (Mg) index amplitude scale factor
    !#                    SOLVAR(2)    Sunspot (SB) index amplitude scale factor
    !#
    !#            For ISOLVAR = 2:
    !#                    SOLVAR(1)    Facular (Mg) index as defined in the NRLSSI2 model;
    !#                                 used for modeling time-specific solar activity
    !#                    SOLVAR(2)    Sunspot (SB) index as defined in the NRLSSI2 model; 
    !#                                 used for modeling time-specific solar activity
    !#
    !#            For ISOLVAR = -1 or 3:
    !#                    SOLVAR(1:14) Band scale factors for modeling spectral variation of 
    !#                                 averaged solar cycle in each shortwave ban
    !real(kind=rb) :: indsolvar(2)
    !# Mg and SB index amplitude scale factors (isolvar=1), or
    !# Mg and SB indices as defined in the NRLSSI2 model (isolvar=2)
    !real(kind=rb) :: bndsolvar(nbndsw)
    !# Solar variability scale factors 
    !# for each shortwave band Dimensions: (nbndsw=14)
    !real(kind=rb) :: solcycfrac 
    !# Solar cycle fraction (0-1); fraction of the way through the mean 11-year
    !# cycle with 0.0 defined as the first day of year 1 and 1.0 defined as the
    !# last day of year 11 (ISOLVAR=1 only). Note that for the combined effect of
    !# the solar constant of 1360.85, and the mean facular brightening and sunspot 
    !# dimming components (without scaling), the minimum total solar irradiance of
    !# 1360.49 occurs at SOLCYCFRAC = 0.0265, and the maximum total solar irradiance 
    !# of 1361.34 occurs at SOLCYCFRAC = 0.3826. 

    ! **(JP)** required for binary reproducibitily whcn using GCC
    real :: inv86400=1.0/86400.0

    iaer=10 !LFR for the 5.0 version

    nofcols=(iz-ia+1)*(jz-ja+1)
    nlay=mzp-1

    !- integer
    allocate(ipos  (nofcols))     ;ipos =0 !integer
    allocate(jpos  (nofcols))     ;jpos =0 !integer
    allocate(imask (nofcols))     ;imask=0 !integer
    !- real
    allocate(play  (nofcols,nlay)  ) ;play  =0.0
    allocate(plev  (nofcols,nlay+1)) ;plev  =0.0
    allocate(tlay  (nofcols,nlay)  ) ;tlay  =0.0
    allocate(tlev  (nofcols,nlay+1)) ;tlev  =0.0
    !
    allocate(h2ovmr(nofcols,nlay))    ;h2ovmr=0.0
    allocate(o3vmr (nofcols,nlay))    ;o3vmr =0.0
    allocate(co2vmr(nofcols,nlay))    ;co2vmr=0.0
    allocate(ch4vmr(nofcols,nlay))    ;ch4vmr=0.0
    allocate(n2ovmr(nofcols,nlay))    ;n2ovmr=0.0
    allocate(o2vmr (nofcols,nlay))    ;o2vmr =0.0
    !
    allocate(cfc11vmr(nofcols,nlay))  ; cfc11vmr=0.0
    allocate(cfc12vmr(nofcols,nlay))  ; cfc12vmr=0.0
    allocate(cfc22vmr(nofcols,nlay))  ; cfc22vmr=0.0
    allocate(ccl4vmr (nofcols,nlay))  ; ccl4vmr =0.0
    allocate(emis    (nofcols,nbndlw)); emis    =0.0
    !

    allocate(alpha_sw(nofcols,nlay)); alpha_sw=0.0
    allocate(cldfmcl_sw(ngptsw,nofcols,nlay)) ;cldfmcl_sw=0.0
    allocate(taucmcl_sw(ngptsw,nofcols,nlay)) ;taucmcl_sw=0.0
    allocate(ssacmcl_sw(ngptsw,nofcols,nlay)) ;ssacmcl_sw=0.0
    allocate(asmcmcl_sw(ngptsw,nofcols,nlay)) ;asmcmcl_sw=0.0
    allocate(fsfcmcl_sw(ngptsw,nofcols,nlay)) ;fsfcmcl_sw=0.0
    allocate(ciwpmcl_sw(ngptsw,nofcols,nlay)) ;ciwpmcl_sw=0.0
    allocate(clwpmcl_sw(ngptsw,nofcols,nlay)) ;clwpmcl_sw=0.0

    allocate(reicmcl_sw(nofcols,nlay)) ;reicmcl_sw=0.0
    allocate(relqmcl_sw(nofcols,nlay)) ;relqmcl_sw=0.0

    allocate(cldfmcl_lw(ngptlw,nofcols,nlay)) ;cldfmcl_lw=0.0
    allocate(taucmcl_lw(ngptlw,nofcols,nlay)) ;taucmcl_lw=0.0
    allocate(ssacmcl_lw(ngptlw,nofcols,nlay)) ;ssacmcl_lw=0.0
    allocate(asmcmcl_lw(ngptlw,nofcols,nlay)) ;asmcmcl_lw=0.0
    allocate(fsfcmcl_lw(ngptlw,nofcols,nlay)) ;fsfcmcl_lw=0.0
    allocate(ciwpmcl_lw(ngptlw,nofcols,nlay)) ;ciwpmcl_lw=0.0
    allocate(clwpmcl_lw(ngptlw,nofcols,nlay)) ;clwpmcl_lw=0.0

    allocate(reicmcl_lw(nofcols,nlay)) ;reicmcl_lw=0.0
    allocate(relqmcl_lw(nofcols,nlay)) ;relqmcl_lw=0.0

    allocate(cldf(nofcols,nlay))      ;cldf=0.0
    allocate(ciwp(nofcols,nlay))      ;ciwp=0.0
    allocate(clwp(nofcols,nlay))      ;clwp=0.0
    allocate(reic(nofcols,nlay))      ;reic=0.0
    allocate(relq(nofcols,nlay))      ;relq=0.0

    allocate(taucl_sw(nbndsw,nofcols,nlay))    ;taucl_sw=0.0
    allocate(ssacl_sw(nbndsw,nofcols,nlay))    ;ssacl_sw=0.0
    allocate(asmcl_sw(nbndsw,nofcols,nlay))    ;asmcl_sw=0.0
    allocate(fsfcl_sw(nbndsw,nofcols,nlay))    ;fsfcl_sw=0.0

    allocate(alpha_lw(nofcols,nlay)); alpha_lw=0.0

    allocate(taucl_lw(nbndlw,nofcols,nlay))    ;taucl_lw=0.0
    allocate(ssacl_lw(nbndlw,nofcols,nlay))    ;ssacl_lw=0.0
    allocate(asmcl_lw(nbndlw,nofcols,nlay))    ;asmcl_lw=0.0
    allocate(fsfcl_lw(nbndlw,nofcols,nlay))    ;fsfcl_lw=0.0

    allocate(tauaer_sw (nofcols,nlay,nbndsw)) ;tauaer_sw =0.0
    allocate(ssaaer_sw (nofcols,nlay,nbndsw)) ;ssaaer_sw =0.0
    allocate(asmaer_sw (nofcols,nlay,nbndsw)) ;asmaer_sw =0.0
    allocate(tauaer_lw (nofcols,nlay,nbndlw)) ;tauaer_lw =0.0
    allocate(ssaaer_lw (nofcols,nlay,nbndlw)) ;ssaaer_lw =0.0
    allocate(asmaer_lw (nofcols,nlay,nbndlw)) ;asmaer_lw =0.0

    allocate(ecaer     (nofcols,nlay,naerec)) ;ecaer     =0.0

    allocate(swhr      (nofcols,nlay))        ;swhr      =0.0
    allocate(swhrc     (nofcols,nlay))        ;swhrc     =0.0
    allocate(hr        (nofcols,nlay))        ;hr        =0.0
    allocate(hrc       (nofcols,nlay))        ;hrc       =0.0

    !
    allocate(swuflx   (nofcols,nlay+2));  swuflx   =0.0 !srf set "nlay+2" for
    allocate(swdflx   (nofcols,nlay+2));  swdflx   =0.0 !    the extra layer
    allocate(swuflxc  (nofcols,nlay+2));  swuflxc  =0.0 !    at the top of
    allocate(swdflxc  (nofcols,nlay+2));  swdflxc  =0.0 !    the model.
    allocate(uflx     (nofcols,nlay+2));  uflx  =0.0 !    Otherwise, it can
    allocate(dflx     (nofcols,nlay+2));  dflx  =0.0 !    be "nlay+1".
    allocate(uflxc    (nofcols,nlay+2));  uflxc    =0.0
    allocate(dflxc    (nofcols,nlay+2));  dflxc    =0.0
    allocate(duflx_dt (nofcols,nlay+2));  duflx_dt =0.0
    allocate(duflxc_dt(nofcols,nlay+2));  duflxc_dt=0.0

    allocate(p     (nofcols,mzp)) ;  p     =0.0
    allocate(t     (nofcols,mzp)) ;  t     =0.0
    allocate(q     (nofcols,mzp)) ;  q     =0.0
    allocate(o3l   (nofcols,mzp)) ;  o3l   =0.0
    allocate(fice  (nofcols,nlay));  fice  =0.0
    allocate(lmixr (nofcols,nlay));  lmixr =0.0
    allocate(psurf (nofcols))     ;  psurf =0.0
    allocate(colrad(nofcols))     ;  colrad=0.0
    allocate(tsfc  (nofcols))     ;  tsfc  =0.0
    allocate(asdir (nofcols))     ;  asdir =0.0
    allocate(aldir (nofcols))     ;  aldir =0.0
    allocate(asdif (nofcols))     ;  asdif =0.0
    allocate(aldif (nofcols))     ;  aldif =0.0
    allocate(coszen(nofcols))     ;  coszen=0.0


    tauc_sw      =0.0
    ssac_sw      =0.0
    asmc_sw      =0.0
    fsfc_sw      =0.0
    tauc_lw      =0.0
    ssac_lw      =0.0
    asmc_lw      =0.0
    fsfc_lw      =0.0
    taucloud_lw  =0.0
    taucloud_sw=0.0
    taucldorig=0.0
    ssacloud_sw=0.0
    ssacloud_lw=0.0
    asmcloud_sw=0.0
    asmcloud_lw=0.0

    !-initialization of rrtm memory
    if(firsttime) call initrrtm()
    !
    !- day of year
    dyofyr=julday(imonth1, idate1, iyear1)

    aot_rrtm_sw=0.0
    aot_rrtm_lw=0.0

    if(firsttime .or. chemistry < 0) then

       tauaer_sw=0.0_rb
       ssaaer_sw=0.0_rb
       asmaer_sw=0.0_rb

       tauaer_lw=0.0_rb
       ssaaer_lw=0.0_rb
       asmaer_lw=0.0_rb

    end if

    !-
    !- column counter
    noc=0
    !- loop over the node domain
    do j=ja,jz  
       do i=ia,iz  
          noc=noc+1
          ipos(noc)=i
          jpos(noc)=j

          !- 2d input data
          coszen(noc)= oneRadiateFields%cosz(i,j)
          tsfc  (noc)=(oneRadiateFields%rlongup(i,j)/stefan)** 0.25e0_rb

          if(.not. firsttime) then
             asdir(noc)=oneRadiateFields%albedt(i,j)
             aldir(noc)=oneRadiateFields%albedt(i,j) !infrared
             asdif(noc)=oneRadiateFields%albedt(i,j)
             aldif(noc)=oneRadiateFields%albedt(i,j) !infrared
             emis(noc,:)=1.-oneRadiateFields%albedt(i,j)
          end if
          !
          !- surface pressure
          psurf(noc)=(oneBasicFields%pi0(1,i,j) + oneBasicFields%pp(1,i,j)+ &
               oneBasicFields%pi0(2,i,j) + oneBasicFields%pp(2,i,j))*cpi*0.50
          psurf(noc)=(p00*psurf(noc)**cpor)*adjust
          !land/ocean mask
          do np=1,npatch  
             if(leaf_g(ngrid)%leaf_class(i,j,np)==0) then
                imask(noc)=0
                cycle
             end if !ocean
             if(leaf_g(ngrid)%leaf_class(i,j,np)==2) then
                imask(noc)=13
                cycle
             end if!land ice
          end do
          colrad(noc)=pi*(1.00e0_rb-(grid_g(ngrid)%glat(i,j)+90.00e0_rb)/180.00e0_rb)

          !- 3d input data
          do k=2,mzp   
             picpi = (oneBasicFields%pi0(k,i,j) + oneBasicFields%pp(k,i,j)) * cpi
             p(noc,k-1) = (p00 * picpi ** cpor)*adjust
             t(noc,k-1) =oneBasicFields%theta(k,i,j) * picpi
             q(noc,k-1) =oneBasicFields%rv   (k,i,j)
             !
             cldf(noc,k-1)=cloud_fraction(k,i,j) !cloud_fraction
             clwp(noc,k-1)=lwl(k,i,j)            !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
             ciwp(noc,k-1)=iwl(k,i,j)            !total cloud ice water (kg/kg for carma and g/m2 for rrtm)

             relq(noc,k-1)=oneMicroFields%rel(k,i,j) ! effective radius liquid

             if(iceflgsw == 1) then

                reic(noc,k-1)=min(130.0,max(13.0,oneMicroFields%rei(k,i,j))) !effective radius ice

             elseif(iceflgsw == 2) then 

                reic(noc,k-1)=min(131.0,max(5.0,oneMicroFields%rei(k,i,j))) !effective radius ice

             elseif(iceflgsw == 3) then 

                reic(noc,k-1)=min(140.0,max(5.0,0.9021*oneMicroFields%rei(k,i,j)-7.0365)) ! generalized effective radius ice

             endif

          end do
          q(noc,mzp) = q(noc,mzp-1)
          p(noc,mzp)=p(noc,mzp-1)-2.*(p(noc,mzp-1)-p(noc,mzp-2))
          t(noc,mzp)=t(noc,mzp-1)+(t(noc,mzp-1)-t(noc,mzp-2))/&
               (p(noc,mzp-1)-p(noc,mzp-2))*&
               (p(noc,mzp  )-p(noc,mzp-1))

          !- these are values in the middle of the atmospheric layers
          do k=1,nlay  
             play(noc,k)=p(noc,k)
             tlay(noc,k)=t(noc,k)
          end do

          !- these are values in the interfaces of the atmos layers
          do k=2,nlay  
             plev(noc,k)=(p(noc,k)+p(noc,k-1))*0.5
             tlev(noc,k)=(t(noc,k)+t(noc,k-1))*0.5
          end do
          !-special treatment for level 1 (interface layer)
          plev(noc,1)=psurf(noc)
          tlev(noc,1)=0.5* ( oneBasicFields%theta(1,i,j)*&
               (oneBasicFields%pi0(1,i,j) + oneBasicFields%pp(1,i,j)) * cpi + &
               
               oneBasicFields%theta(2,i,j)*&
               (oneBasicFields%pi0(2,i,j) + oneBasicFields%pp(2,i,j)) * cpi   &
               )

          !-special treatment for level mzp (interface layer)
          plev(noc,mzp)=plev(noc,mzp-1)-2.*(plev(noc,mzp-1)-play(noc,nlay))

          tlev(noc,mzp)=tlev(noc,mzp-1)+(tlev(noc,mzp-1)-tlev(noc,mzp-2))/&
               (plev(noc,mzp-1)-plev(noc,mzp-2))*&
               (plev(noc,mzp  )-plev(noc,mzp-1))




          if(.not. (firsttime .or. chemistry < 0)) then
             do k=2,mzp  
                do np=1,nbndsw  
                   tauaer_sw(noc,k-1,np)=optProp(1,np)%tauaer(k,i,j)
                   aot_rrtm_sw(i,j,np) = aot_rrtm_sw(i,j,np)+tauaer_sw(noc,k-1,np)
                   !Total of AOT integrated in column for output
                   if(k==mzp) carma(ngrid)%aot(i,j,np)=aot_rrtm_sw(i,j,np)
                   ssaaer_sw(noc,k-1,np)=optProp(ngrid,np)%ssa(k,i,j)
                   asmaer_sw(noc,k-1,np)=optProp(ngrid,np)%asp(k,i,j)
                end do

                do np=1,nbndlw  
                   tauaer_lw(noc,k-1,np)=optProp(1,np+nbndsw)%tauaer(k,i,j)
                   aot_rrtm_lw(i,j,np) = aot_rrtm_lw(i,j,np)+tauaer_lw(noc,k-1,np)
                   !Total of AOT integrated in column for output
                   if(k==mzp) carma(ngrid)%aot(i,j,np+nbndsw)=aot_rrtm_lw(i,j,np)
                   ssaaer_lw(noc,k-1,np)=optProp(ngrid,np+nbndsw)%ssa(k,i,j)
                   asmaer_lw(noc,k-1,np)=optProp(ngrid,np+nbndsw)%asp(k,i,j)
                end do
             end do
          end if
       end do
    end do


    ecaer =0.0_rb

    do noc=1,nofcols  
       !calculating sig and delta sig
       do i=2,mzp  
          sig(i)=plev(noc,i)/ plev(noc,1) !plev(2,mzp)
       end do
       sig(1)=1.
       sig(mzp+1)=0.

       do i=1,mzp  
          if (sig(i)<0.1) exit
       end do
       nls=mzp-i+1
       do k=1,mzp   
          delsig(k)=sig(k)-sig(k+1)
       end do
       delsig(mzp+1)=delsig(mzp)
       do k=1,mzp   
          sigmid(k)=(sig(k)+sig(k+1))/2
       end do
       sigmid(mzp+1)=sigmid(mzp)
       call initgetoz_rrtm(365.2500_real64,mzp,sig)

       !kml calling getoz for both cases
       doy=real(dyofyr,real64)
       call getoz_rrtm (mzp,sigmid,colrad(noc),doy,o3l(noc,:))
    enddo

    if(o3ispresent) then

       !srf o3l(:,1)=o3l(:,2)    !kml eliminate o3 surface layer of the climatological profile

       do noc=1,nofcols  
          do k=2,mzp  
             !kml sum o3 climatological and model profiles
             o3vmr(noc,k-1)=chem1_g(o3pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*       &
                  ro3ar*oneBasicFields%dn0(k,ipos(noc),jpos(noc))*factorcor +&
                  o3l(noc,flip(k-1))*ro3ar
          end do
       end do
    else
       do noc=1,nofcols  
          do k=1,mzp-1   
             o3vmr(noc,k)=o3l(noc,flip(k))*ro3ar !!! *1.e+4  << check!
             !print*,"o3vmr=",noc,k,real(o3vmr(noc,k),4)
          end do
       end do
    end if

    if(co2ispresent) then
       do noc=1,nofcols
          do k=2,mzp
             co2vmr(noc,k-1)=chem1_g(co2pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*rco2ar*factorcor
          end do
       end do
    else
       co2vmr=379.e-6*rco2ar !fixed value from ipcc
    end if

    if(ch4ispresent) then
       do noc=1,nofcols
          do k=2,mzp
             ch4vmr(noc,k-1)=chem1_g(ch4pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*rch4ar*factorcor
          end do
       end do
    else
       ch4vmr=1774.e-9*rch4ar !fixed value from ipcc
    end if

    if(no2ispresent) then
       do noc=1,nofcols
          do k=2,mzp
             n2ovmr(noc,k-1)=chem1_g(no2pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*rno2ar*factorcor
          end do
       end do
    else
       n2ovmr=319.e-9 *rno2ar!fixed value from ipcc
    end if

    if(o2ispresent) then
       do noc=1,nofcols
          do k=2,mzp
             o2vmr(noc,k-1)=chem1_g(o2pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*ro2ar*factorcor
          end do
       end do
    else
       o2vmr=0.209488*ro2ar !fixed value from ipcc
    end if

    do noc=1,nofcols
       do k=2,mzp
          h2ovmr(noc,k-1)=q(noc,k)*rh2oar
       end do
    end do

    cfc11vmr=0.251e-9 !srf* factorcor
    cfc12vmr=0.538e-9 !srf* factorcor
    cfc22vmr=0.169e-9 !srf* factorcor
    ccl4vmr =0.093e-9 !srf* factorcor


    !- longwave section

    do noc=1,nofcols
       call cldprop(nlay, inflglw, iceflglw, liqflglw, cldf(noc,:), tauc_lw, &
            ciwp(noc,:), clwp(noc,:), &
            reic(noc,:), relq(noc,:), ncbandslw, taucloud_lw)
       do np=1,nbndlw
          taucl_lw(np,noc,:)=taucloud_lw(:,np)

       end do
    end do

    changeseed = permuteseed + icount 

    call mcica_subcol_lw(iplon, nofcols, nlay, icld, changeseed, irng, play, & 
         cldf, ciwp, clwp, reic, relq, taucl_lw,  &
         alpha_lw, &
         cldfmcl_lw, ciwpmcl_lw, clwpmcl_lw, reicmcl_lw, relqmcl_lw, &
         taucmcl_lw)

    !- longwave calculations
    call rrtmg_lw (&
         nofcols   , &!01
         nlay      , &!02
         icld      , &!03
         idrv      , &!04
         play      , &!05
         plev      , &!06
         tlay      , &!07
         tlev      , &!08
         tsfc      , &!09
         h2ovmr    , &!10
         o3vmr     , &!11
         co2vmr    , &!12
         ch4vmr    , &!13
         n2ovmr    , &!14
         o2vmr     , &!15
         cfc11vmr  , &!16
         cfc12vmr  , &!17
         cfc22vmr  , &!18
         ccl4vmr   , &!19
         emis      , &!20
         inflglw   , &!21
         iceflglw  , &!22
         liqflglw  , &!23
         cldfmcl_lw, &!24
         taucmcl_lw, &!25
         ciwpmcl_lw, &!26
         clwpmcl_lw, &!27
         reicmcl_lw, &!28
         relqmcl_lw, &!29
         tauaer_lw , &!30
         uflx      , &!33
         dflx      , &!34
         hr        , &!35
         uflxc     , &!36
         dflxc     , &!37
         hrc       , &!38
         duflx_dt  , &!39
         duflxc_dt   &!40
         )

    !--- shortwave section
    do noc=1,nofcols
       call cldprop_sw(nlay, inflgsw, iceflgsw, liqflgsw, cldf(noc,:) , &
            tauc_sw, ssac_sw, asmc_sw, fsfc_sw, ciwp(noc,:), &
            clwp(noc,:),reic(noc,:), relq(noc,:), &
            taucldorig, taucloud_sw, ssacloud_sw, asmcloud_sw)
       do np=1,nbndsw

          taucl_sw(np,noc,:)= taucloud_sw(:,np+15)
          ssacl_sw(np,noc,:)= ssacloud_sw(:,np+15)
          asmcl_sw(np,noc,:)= asmcloud_sw(:,np+15)

       end do
    end do

    changeseed = changeseed + ngpt

    call mcica_subcol_sw(iplon, nofcols, nlay, icld, changeseed, irng, play, & 
         cldf, ciwp, clwp, reic, relq, taucl_sw, ssacl_sw, asmcl_sw, fsfcl_sw, &
         alpha_sw, &
         cldfmcl_sw, ciwpmcl_sw, clwpmcl_sw, reicmcl_sw, relqmcl_sw, &
         taucmcl_sw, ssacmcl_sw, asmcmcl_sw, fsfcmcl_sw)

    changeseed = changeseed - ngpt

    !examples of dump 3d
    !call dumpvar(real(reicmcl),'rei','-001',nofcols,nlay)
    !call dumpvar(real(relqmcl),'rel','-001',nofcols,nlay)

    !- shortwave calculations
    call rrtmg_sw &
         (nofcols , &!01                         
         nlay    , &!02
         icld    , &!03
         iaer    , &! LFR for the 5.0 version
         play    , &!04
         plev    , &!05
         tlay    , &!06
         tlev    , &!07
         tsfc    , &!08 *
         h2ovmr  , &!09
         o3vmr   , &!10
         co2vmr  , &!11
         ch4vmr  , &!12
         n2ovmr  , &!13
         o2vmr   , &!14 *
         asdir   , &!15
         asdif   , &!16
         aldir   , &!17
         aldif   , &!18 *
         coszen  , &!19
         adjes   , &!20
         dyofyr  , &!21
         scon    , &!22 *
         isolvar , & !LFR for the 5.0 version
         inflgsw , &!23
         iceflgsw, &!24
         liqflgsw, &!25
         cldfmcl_sw, &!26 *
         taucmcl_sw, &!27
         ssacmcl_sw, &!28
         asmcmcl_sw, &!29
         fsfcmcl_sw, &!30 *
         ciwpmcl_sw, &!31
         clwpmcl_sw, &!32
         reicmcl_sw, &!33
         relqmcl_sw, &!34 *
         tauaer_sw, &!35
         ssaaer_sw, &!36
         asmaer_sw, &!37
         ecaer   , &!38 *
         swuflx  , &!40
         swdflx  , &!41
         swhr    , &!42 *
         swuflxc , &!43
         swdflxc , &!44
         swhrc     &!45 *
                                !bndsolvar, & !LFR for the 5.0 version
                                !indsolvar, & !LFR for the 5.0 version
                                !solcycfrac & !LFR for the 5.0 version
         )

    !-set initialization variable to false, so
    !- next timestep it will no be performed anymore.
    firsttime=.false.

    !- sending out theta tendency and surface radiative fluxes
    noc=0
    do j=ja,jz
       do i=ia,iz
          noc=noc+1
          oneRadiateFields%rlong (i,j)=   dflx(noc,1)
          oneRadiateFields%rshort(i,j)= swdflx(noc,1)
          if(swdflx(noc,1)<.5) oneRadiateFields%rshort(i,j)=0.0
          do k=2,mzp-1
             !- radiative tendency of temperature (Kelvin/sec)
             ! **(JP)** required for binary reproducibility when using GCC
             ! **(JP)** use next statement if gfortran version is 8-5
             oneRadiateFields%fthrd(k,i,j)=(swhr(noc,k-1)+hr(noc,k-1))/86400.
             ! **(JP)** use next statement if gfortran version is newer (required for 13.2)
             !!$oneRadiateFields%fthrd(k,i,j)=(swhr(noc,k-1)+hr(noc,k-1))*inv86400
          enddo
          oneRadiateFields%fthrd(1  ,i,j)= oneRadiateFields%fthrd(2    ,i,j)
          oneRadiateFields%fthrd(mzp,i,j)= oneRadiateFields%fthrd(mzp-1,i,j)
       end do
    end do

    deallocate(ipos  )
    deallocate(jpos  )
    deallocate(imask )
    deallocate(play  )
    deallocate(plev  )
    deallocate(tlay  )
    deallocate(tlev  )
    deallocate(h2ovmr)
    deallocate(o3vmr )
    deallocate(co2vmr)
    deallocate(ch4vmr)
    deallocate(n2ovmr)
    deallocate(o2vmr )
    deallocate(cfc11vmr)
    deallocate(cfc12vmr)
    deallocate(cfc22vmr)
    deallocate(ccl4vmr )
    deallocate(emis )

    deallocate(cldfmcl_sw)
    deallocate(taucmcl_sw)
    deallocate(ssacmcl_sw)
    deallocate(asmcmcl_sw)
    deallocate(fsfcmcl_sw)
    deallocate(ciwpmcl_sw)
    deallocate(clwpmcl_sw)

    deallocate(reicmcl_sw)
    deallocate(relqmcl_sw)

    deallocate(cldfmcl_lw)
    deallocate(taucmcl_lw)
    deallocate(ssacmcl_lw)
    deallocate(asmcmcl_lw)
    deallocate(fsfcmcl_lw)
    deallocate(ciwpmcl_lw)
    deallocate(clwpmcl_lw)

    deallocate(reicmcl_lw)
    deallocate(relqmcl_lw)

    deallocate(cldf)
    deallocate(ciwp)
    deallocate(clwp)
    deallocate(reic)
    deallocate(relq)

    deallocate(taucl_sw)
    deallocate(ssacl_sw)
    deallocate(asmcl_sw)
    deallocate(fsfcl_sw)

    deallocate(taucl_lw)
    deallocate(ssacl_lw)
    deallocate(asmcl_lw)
    deallocate(fsfcl_lw)

    deallocate(tauaer_sw) 
    deallocate(ssaaer_sw) 
    deallocate(asmaer_sw) 
    deallocate(tauaer_lw) 
    deallocate(ssaaer_lw) 
    deallocate(asmaer_lw) 

    deallocate(ecaer)

    deallocate(swhr      )
    deallocate(swhrc     )
    deallocate(hr        )
    deallocate(hrc       )

    deallocate(swuflx   )
    deallocate(swdflx   )
    deallocate(swuflxc  )
    deallocate(swdflxc  )
    deallocate(uflx     )
    deallocate(dflx     )
    deallocate(uflxc    )
    deallocate(dflxc    )
    deallocate(duflx_dt )
    deallocate(duflxc_dt)

    deallocate(p     )
    deallocate(t     )
    deallocate(q     )
    deallocate(o3l   )
    deallocate(fice  )
    deallocate(lmixr )
    deallocate(psurf )
    deallocate(colrad)
    deallocate(tsfc  )
    deallocate(asdir )
    deallocate(aldir )
    deallocate(asdif )
    deallocate(aldif )
    deallocate(coszen)
  end subroutine radrrtmdrv

  ! ****************************************************************************

  subroutine cloud_prop_rrtm(m1,m2,m3,ia, iz, ja, jz &
       , cloud_fraction  &
       , rain            &
       , lwl             &
       , iwl,            &
       oneNamelistFile, oneBasicFields, oneMicVars, oneMicroFields, &
       oneCuParmFields)

    integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicControl), pointer, intent(in) :: oneMicVars
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    type(CuParmFields), pointer, intent(in) :: oneCuParmFields
    
    real, intent(out), dimension(m1,m2,m3) :: cloud_fraction !cloud_fraction
    real, intent(out), dimension(m2,m3   ) :: rain !total rain water
    real, intent(out), dimension(m1,m2,m3) :: lwl !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
    real, intent(out), dimension(m1,m2,m3) :: iwl !total cloud ice water (kg/kg for carma and g/m2 for rrtm)

    !-- local variables
    integer :: i,j,k,k700,kp1,kdthdp

    real, parameter :: rhminl = .90              ! minimum rh for low stable clouds
    real, parameter :: rhminh = .80              ! minimum rh for high stable clouds
    real, parameter :: sh1 = 0.07 ,sh2= 500.0   ! parameters for shallow convection cloud fraction
    real, parameter :: dp1 = 0.14 ,dp2= 500.0   ! parameters for deep convection cloud fraction
    real, parameter :: premit= 750.e2              ! top pressure bound for mid level cloud
    real, parameter :: pnot = 1.e5       ! reference pressure
    real, parameter :: lapse = 6.5e-3    ! u.s. standard atmsophere lapse rate
    real, parameter :: premib = 750.e2   ! bottom pressure bound of middle cloud
    real, parameter :: pretop = 1.0e2    ! pressure bounding high cloud
    real, parameter :: abeta = 0.07
    real, parameter :: bbeta = -0.14
    real, parameter :: pi = 3.14159265358979323846
    real, parameter :: bx = 100.* (3./(4.*pi))**(1./3.)
    real, parameter :: r13 = 1./3.
    real, parameter :: r13bbeta = 1./3. - 0.14

    real, dimension(m1) :: press,rh,cldst,concld,rhu00,rpdeli
    real, dimension(m1,m2,m3) :: dummy_vec
    real, dimension(m2,m3   ) :: upmf,upmfsh
    real, dimension(m1,m2,m3) :: zup,zupshallow,clwup,clwupsh
    real :: ocnfrac,picpi,temp,dztri,strat,rhpert,shallowcu,deepcu,upmfshx,upmfx&
         ,rhwght,rhdif,rhlim,ps,thetas,dthdpmn,dthdp,dummy,bb,factor
    logical cldbnd          ! region below high cloud boundary
    real, external :: rs
    character(len=3) :: cm1
    !- tuning parameters to include direct coupling between cupar and radiation
    real, parameter :: tun_rad_shall=0.02, tun_rad_deep =0.02
    integer,parameter :: coupl_rad_cupar=1 ! 0 -no, 1-yes
    integer :: n
    ! set defaults for rhu00
    rhu00(:) = 2.0
    ! define rh perturbation in order to estimate rhdfda
    rhpert = 0.01

    !--- initialization of cuparm parameters
    !--- convective cloud fraction
    !clcn = 0.0
    !
    if(oneNamelistFile%nnqparm(ngrid) == 6 .or. &
         oneNamelistFile%nnqparm(ngrid) == 8) then
       upmf      (     1:m2,1:m3)= xmb4d(     1:m2,1:m3,1,ngrid) !- mass flux deep    convection
       upmfsh    (     1:m2,1:m3)= xmb4d(     1:m2,1:m3,2,ngrid) !- mass flux shallow convection
       zup       (1:m1,1:m2,1:m3)= zup5d(1:m1,1:m2,1:m3,1,ngrid) !- normalized mass flux
       zupshallow(1:m1,1:m2,1:m3)= zup5d(1:m1,1:m2,1:m3,2,ngrid) !- normalized mass flux

       if(coupl_rad_cupar == 1 ) then
          clwup     (1:m1,1:m2,1:m3)= tun_rad_deep *clwup5d(1:m1,1:m2,1:m3,1,ngrid)
          clwupsh   (1:m1,1:m2,1:m3)= tun_rad_shall*clwup5d(1:m1,1:m2,1:m3,2,ngrid)
       endif
       !--- convective cloud fraction (to implement this scheme, CLCN needs to be transported as a scalar)
    else
       upmf      (     1:m2,1:m3)=0.0
       upmfsh    (     1:m2,1:m3)=0.0
       zup       (1:m1,1:m2,1:m3)=0.0
       zupshallow(1:m1,1:m2,1:m3)=0.0
       if(coupl_rad_cupar == 1 ) then
          clwup     (1:m1,1:m2,1:m3)=0.0
          clwupsh   (1:m1,1:m2,1:m3)=0.0
       endif
    endif
    do j=ja,jz
       do i=ia,iz

          ! evaluate potential temperature and relative humidity
          do k=1,m1
             picpi = (oneBasicFields%pi0(k,i,j) + oneBasicFields%pp(k,i,j)) * cpi
             press(k) = p00 * picpi ** cpor
             temp = oneBasicFields%theta(k,i,j) * picpi
             rh(k) =min(1.,max(0.05,oneBasicFields%rv(k,i,j)/rs(press(k),temp)))
             !------
             !
             cloud_fraction(k,i,j)  = 0.
             cldst (k)      = 0.
             concld(k)      = 0.
          enddo
          ps    =0.5*(press(1)    +press(2))
          thetas=0.5*(oneBasicFields%theta(1,i,j)+oneBasicFields%theta(2,i,j))
          ocnfrac=0.
          if(leaf_g(ngrid)%patch_area(i,j,1)>0.99) ocnfrac=1. !flag < 1 para land
          !flag  =1 para water

          do k=1,m1-1
             rpdeli(k) = 1./(press(k+1) - press(k))
          end do
          !
          ! cloud mass flux in si units of kg/m2/s; should produce typical numbers of 20%
          ! shallow and deep convective cloudiness are evaluated separately (since processes
          ! are evaluated separately) and summed
          !
          do k=2,m1
             shallowcu=0.0; upmfshx=0.0
             if(upmfsh(i,j)>1.e-8) then
                upmfshx=zupshallow(k-1,i,j) * upmfsh(i,j)
                shallowcu = max(0.0,min(sh1*log(1.0+sh2*upmfshx),0.30))
             endif
             deepcu = 0.
             if(upmf(i,j) >1.e-8) then
                upmfx     = zup(k-1,i,j) * upmf(i,j)

                !-orig deepcu = max(0.0,min(dp1*log(1.0+dp2*(cmfmc(i,k+1)-cmfmc2(i,k+1))),0.60))
                ! check possibility of (upmfx-upmfshx) <0. => log (x<0)
                deepcu = max(0.0,min(dp1*log(1.0+dp2*max(0.,(upmfx-upmfshx))),0.60))
             endif
             concld(k) = min(shallowcu + deepcu,0.80)
             rh(k) = (rh(k) - concld(k))/(1.0 - concld(k))
          end do

          do k=2,m1
             kp1 = min(k + 1,m1)
             !
             cldbnd = press(k).ge.pretop
             if ( press(k).ge.premib ) then
                !==============================================================
                ! this is the low cloud (below premib) block
                !==============================================================
                ! enhance low cloud activation over land with no snow cover
                if ( ocnfrac < 0.999 ) then !.and. (snowh(i) <= 0.000001)) then
                   rhlim = rhminl - 0.10
                else
                   rhlim = rhminl
                endif
                !
                rhdif = (rh(k) - rhlim)/(1.0-rhlim)
                cloud_fraction(k,i,j) = min(0.999,(max(rhdif,0.0))**2)
             else if ( press(k).lt.premit ) then
                !==============================================================
                ! this is the high cloud (above premit) block
                !==============================================================
                !
                rhlim = rhminh
                !
                rhdif = (rh(k) - rhlim)/(1.0-rhlim)
                cloud_fraction(k,i,j) = min(0.999,(max(rhdif,0.0))**2)
             else
                !==============================================================
                ! this is the middle cloud block
                !==============================================================
                !
                !       linear rh threshold transition between thresholds for low & high cloud
                !
                rhwght = (premib-(max(press(k),premit)))/(premib-premit)

                if ( ocnfrac < 0.999 ) then !if (land(i) .and. (snowh(i) <= 0.000001)) then
                   rhlim = rhminh*rhwght + (rhminl - 0.10)*(1.0-rhwght)
                else
                   rhlim = rhminh*rhwght + rhminl*(1.0-rhwght)
                endif
                rhdif = (rh(k) - rhlim)/(1.0-rhlim)
                cloud_fraction(k,i,j) = min(0.999,(max(rhdif,0.0))**2)
             end if
             ! save rhlim to rhu00, it handles well by itself for low/high cloud
             !
             rhu00(k)=rhlim
             !==================================================================================
          end do


          !--- stratus
          ! find most stable level below 750 mb for evaluating stratus regimes
          ! nothing triggers unless a stability greater than this minimum threshold is found
          dthdpmn = -0.125
          kdthdp  = 0

          do k=2,m1
             if (press(k) >= premib .and. ocnfrac.gt. 0.01) then
                ! i think this is done so that dtheta/dp is in units of dg/mb (jjh)
                dthdp = 100.0*(oneBasicFields%theta(k,i,j) - oneBasicFields%theta(k-1,i,j))*rpdeli(k-1)
                if (dthdp < dthdpmn) then
                   dthdpmn = dthdp
                   kdthdp  = k     ! index of interface of max inversion
                endif
             endif
          enddo

          ! also check between the bottom layer and the surface
          ! only perform this check if the criteria were not met above

          if ( kdthdp .eq. 0 .and. ocnfrac.gt.0.01) then
             dthdp = 100.0 * (thetas - oneBasicFields%theta(m1,i,j)) / (ps-press(m1))
             if (dthdp < dthdpmn) then
                dthdpmn = dthdp
                kdthdp  = m1     ! index of interface of max inversion
             endif
          endif
          do k=2,m1-1
             k700=-99
             if(0.01*press(k) .le.  700.) then
                k700=k
                exit
             endif
          enddo
          if( k700 > m1 .or.   k700 < 1) then
             if (mynum == 1) then
                write (*,fmt='(a,i3.3,a,i3.3)') "wrong k700 at cloud_prop routine for i=",i,"; j=",j
                write (*,fmt='(a,i3.3,a,i8)') "m1=",m1,"; k700=",k700
                write (*,fmt='(a,45(e9.3,1x))') "press=",(press(k),k=1,m1)
                write (*,"(2(a,e15.7))") "p00=",p00,"; cpor=",cpor
                write (*,fmt='(a,45(e9.3,1x))') "oneBasicFields%pi0=",(oneBasicFields%pi0(k,i,j),k=1,m1)
                write (*,fmt='(a,45(e9.3,1x))') "oneBasicFields%pp=",(oneBasicFields%pp(k,i,j),k=1,m1)
             end if
             stop
          end if
          if (kdthdp /= 0) then
             k = kdthdp
             kp1 = min(k+1,m1)
             ! note: strat will be zero unless ocnfrac > 0.01
             strat = min(1.,max(0., ocnfrac * ((oneBasicFields%theta(k700,i,j)-thetas)*0.057-0.5573)))
             !
             ! assign the stratus to the layer just below max inversion
             ! the relative humidity changes so rapidly across the inversion
             ! that it is not safe to just look immediately below the inversion
             ! so limit the stratus cloud by rh in both layers below the inversion
             !
             cldst(k) = min(strat,max(rh(k),rh(kp1)))
          endif

          !
          ! aggregate cloud contributions (cldst should be zero everywhere except at level kdthdp(i))
          !
          do k=1,m1
             !
             !       which is greater; standard layered cloud amount or stratocumulus diagnosis
             !
             cloud_fraction(k,i,j) = max(cloud_fraction(k,i,j),cldst(k))
             !
             !       add in the contributions of convective cloud (determined separately and accounted
             !       for by modifications to the large-scale relative humidity.
             !
             cloud_fraction(k,i,j) = min(cloud_fraction(k,i,j)+concld(k), 1.0)
          end do


       enddo
    enddo


    !-------start calculation of cloud ....


    dummy_vec = 0.0
    ! if level == 1 do nothing
    if (oneMicVars%level==2) then
       lwl(1:m1,ia:iz,ja:jz) = oneMicroFields%rcp(1:m1,ia:iz,ja:jz)

       if (oneNamelistFile%nnqparm(ngrid)/=0) then
          rain(ia:iz,ja:jz)= oneCuParmFields%conprr(ia:iz,ja:jz)* 3600.
       endif

    elseif (oneMicVars%level>=3) then
       if (oneNamelistFile%nnqparm(ngrid)/=0) then
          rain(ia:iz,ja:jz) = oneCuParmFields%conprr(ia:iz,ja:jz) + &
               oneMicroFields%pcpg(ia:iz,ja:jz)
       else
          rain(ia:iz,ja:jz) = oneMicroFields%pcpg(ia:iz,ja:jz)
       endif
       rain(ia:iz,ja:jz) = rain(ia:iz,ja:jz)*3600.

       if (oneMicVars%icloud>0) then
          lwl(1:m1,ia:iz,ja:jz) = lwl(1:m1,ia:iz,ja:jz) + oneMicroFields%rcp(1:m1,ia:iz,ja:jz)
       end if
       if (oneMicVars%igraup>0) then
          if(oneMicVars%mcphys_type == 0) then
             do k=1,m1
                do i=ia,iz
                   do j=ja,jz
                      call qtc(oneMicroFields%q6(k,i,j), dummy,dummy_vec(k,i,j))
                   enddo
                enddo
             enddo
          elseif(oneMicVars%mcphys_type == 2 .or. oneMicVars%mcphys_type ==3 .or. oneMicVars%mcphys_type == 4) then  !srf -gthompson/gfdl microphysics - graupel only in ice phase
             dummy_vec=0.0
          endif
          lwl(1:m1,ia:iz,ja:jz) = dummy_vec(1:m1,ia:iz,ja:jz)*oneMicroFields%rgp(1:m1,ia:iz,ja:jz) &
               + lwl(1:m1,ia:iz,ja:jz) !kg/kg
          dummy_vec(:,:,:)      = 1. - dummy_vec(:,:,:)
          iwl(1:m1,ia:iz,ja:jz) = dummy_vec(1:m1,ia:iz,ja:jz)*oneMicroFields%rgp(1:m1,ia:iz,ja:jz) &
               + iwl(1:m1,ia:iz,ja:jz) !kg/kg
       endif
       if (oneMicVars%ihail>0) then
          dummy_vec = 0.
          do k=1,m1
             do i=ia,iz
                do j=ja,jz
                   call qtc(oneMicroFields%q7(k,i,j), dummy, dummy_vec(k,i,j))
                enddo
             enddo
          enddo

          lwl(1:m1,ia:iz,ja:jz) = dummy_vec(1:m1,ia:iz,ja:jz)*oneMicroFields%rhp(1:m1,ia:iz,ja:jz)&
               + lwl(1:m1,ia:iz,ja:jz)

          dummy_vec(:,:,:) = 1.0 - dummy_vec(:,:,:)

          iwl(1:m1,ia:iz,ja:jz) = dummy_vec(1:m1,ia:iz,ja:jz)*oneMicroFields%rhp(1:m1,ia:iz,ja:jz) &
               + iwl(1:m1,ia:iz,ja:jz)   !kg/kg
       endif
       if (oneMicVars%iaggr>0) &
            iwl(1:m1,ia:iz,ja:jz) = iwl(1:m1,ia:iz,ja:jz) + oneMicroFields%rap(1:m1,ia:iz,ja:jz)   !kg/kg

       if (oneMicVars%isnow>0) &
            iwl(1:m1,ia:iz,ja:jz) = iwl(1:m1,ia:iz,ja:jz) + oneMicroFields%rsp(1:m1,ia:iz,ja:jz)   !kg/kg

       if (oneMicVars%ipris>0) &
            iwl(1:m1,ia:iz,ja:jz) = iwl(1:m1,ia:iz,ja:jz) + oneMicroFields%rpp(1:m1,ia:iz,ja:jz)  !kg/kg
    endif
    !- making direct couplig between liq/ice water from cupar to radiation
    if(coupl_rad_cupar == 1 ) then
       do j = ja,jz
          do i = ia,iz
             do k = 1,m1
                temp = oneBasicFields%theta(k,i,j)*(oneBasicFields%pp(k,i,j)+oneBasicFields%pi0(k,i,j))*cpi ! air temp (kelvin)

                if(temp .gt. 253.)then
                   lwl(k,i,j)=lwl(k,i,j)+ clwup(k,i,j)+  clwupsh(k,i,j)
                else
                   iwl(k,i,j)=iwl(k,i,j)+ clwup(k,i,j)+  clwupsh(k,i,j)
                endif
             enddo
          end do
       enddo
    endif

    do j=ja,jz
       do i=ia,iz
          do k=1,m1
             lwl(k,i,j) = max(0.,lwl(k,i,j))  !kg/kg
             iwl(k,i,j) = max(0.,iwl(k,i,j))  !kg/kg
          enddo
       enddo
    enddo
    do j=ja,jz
       do i=ia,iz
          rain(i,j) = max(0.,rain(i,j))
       enddo
    enddo

    !- rrtmg radiation  -  continue calculation of cloud effective radius

    if (oneNamelistFile%ilwrtyp==6 .or. oneNamelistFile%iswrtyp==6) then

       if(oneMicVars%mcphys_type == 2 .or. oneMicVars%mcphys_type == 4 ) then
          !- rei and rel are calculated by gt microphysics (no-aer option)

          do j=ja,jz
             do i=ia,iz
                do k=1,m1

                   picpi = (oneBasicFields%pi0(k,i,j) + oneBasicFields%pp(k,i,j)) * cpi
                   press(k) = p00 * picpi ** cpor
                   temp = oneBasicFields%theta(k,i,j) * picpi

                   lwl(k,i,j) = oneBasicFields%dn0(k,i,j)*lwl(k,i,j)* 1.e+3  !g/m3
                   iwl(k,i,j) = oneBasicFields%dn0(k,i,j)*iwl(k,i,j)* 1.e+3  !g/m3
                   dztri= grid_g(ngrid)%rtgt(i,j)/dzt(k)

                   if(iwl(k,i,j)<1.0e-6 .or. temp>273.0) iwl(k,i,j)=0.0
                   lwl(k,i,j) = lwl(k,i,j)* dztri  !g/m2
                   iwl(k,i,j) = iwl(k,i,j)* dztri  !g/m2
                enddo
             enddo
          enddo

       else
          !- for all other options, rei and rel are calculated in the section below
          !
          do j=ja,jz
             do i=ia,iz
                do k=1,m1

                   !- liquid cloud effective radius -----
                   !- [liu&daum, 2000 and 2005. liu et al 2008]
                   !- cloud drop number concentration
                   dummy_vec(k,i,j) = oneMicroFields%ccp(k,i,j) * oneBasicFields%dn0(k,i,j) * 1.e-6  ! #/cm3
                   dummy_vec(k,i,j) = max(150.,dummy_vec(k,i,j))

                   !
                   !- liquid cloud effective radius ----- [liu&daum, 2000 and 2005]

                   lwl(k,i,j) = oneBasicFields%dn0(k,i,j)*lwl(k,i,j)* 1.e+3  !g/m3
                   iwl(k,i,j) = oneBasicFields%dn0(k,i,j)*iwl(k,i,j)* 1.e+3  !g/m3

                   !- u[lwl] = g/m3 , u[dummy_vec]= #/cm^3
                   !- the cte bx with the units above, provides rel in 10^-6 meter
                   !  oneMicroFields%rel(k,i,j)= bx *  ( lwl(k,i,j) /dummy_vec(k,i,j))**r13 &
                   !- the factor below is adimensional e provides correction
                   !- for dispersion of the cloud spectrum:
                   !- prefactor of liu et al (2008) - lwl must be in g/cm^3
                   !          *abeta*(1.e-6*lwl(k,i,j) /dummy_vec(k,i,j))**bbeta
                   !- expression to avoid nan when lwl = 0. in prefactor
                   oneMicroFields%rel(k,i,j)= bx *  ( lwl(k,i,j) /dummy_vec(k,i,j))**r13bbeta &
                        *abeta*6.92 !6.92=(1.e-6)**bbeta

                   ! rel is limited between 2.5 and 60 micrometers as
                   ! required by rrtm parameterization
                   oneMicroFields%rel(k,i,j) = max(2.5, min( 60.0, oneMicroFields%rel(k,i,j) ) )

                   !------ice cloud effective radius ----- [klaus wyser, 1998]

                   picpi = (oneBasicFields%pi0(k,i,j) + oneBasicFields%pp(k,i,j)) * cpi
                   press(k) = p00 * picpi ** cpor
                   temp = oneBasicFields%theta(k,i,j) * picpi
                   dztri= grid_g(ngrid)%rtgt(i,j)/dzt(k)
                   if(iwl(k,i,j)<1.0e-6 .or. temp>273.0) then
                      oneMicroFields%rei(k,i,j)=5.0
                      iwl(k,i,j)=0.0
                   else
                      bb = -2. + log10(iwl(k,i,j)/50.)*(1.e-3*(273.15-max(210.15,temp))**1.5)
                      oneMicroFields%rei(k,i,j) =377.4 + 203.3 * bb+ 37.91 * bb **2 + 2.3696 * bb **3
                   endif
                   lwl(k,i,j) = lwl(k,i,j)* dztri  !g/m2
                   iwl(k,i,j) = iwl(k,i,j)* dztri  !g/m2
                enddo
             enddo
          enddo
       endif

    endif
    return
    print*,"=============mcphys-radiation coupling ==================="
    print*,"max-min  cl_frac:  ",maxval(cloud_fraction(2:m1,ia:iz,ja:jz)),minval(cloud_fraction(2:m1,ia:iz,ja:jz)) !cloud_fraction
    print*,"max-min     rain:  ",maxval(rain(ia:iz,ja:jz)),minval(rain(ia:iz,ja:jz)) !total rain water
    print*,"max-min lwl: ",maxval(lwl(2:m1,ia:iz,ja:jz)),minval(lwl(2:m1,ia:iz,ja:jz)) !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
    print*,"max-min iwl: ",maxval(iwl(2:m1,ia:iz,ja:jz)),minval(iwl(2:m1,ia:iz,ja:jz))     !total cloud ice water (kg/kg for carma and g/m2 for rrtm)
    print*,"max-min oneMicroFields%rel: ",maxval(oneMicroFields%rel(2:m1,ia:iz,ja:jz)),minval(oneMicroFields%rel(2:m1,ia:iz,ja:jz))  !total cloud liquid water
    print*,"max-min rei: ",maxval(oneMicroFields%rei(2:m1,ia:iz,ja:jz)),minval(oneMicroFields%rei(2:m1,ia:iz,ja:jz))  !total cloud ice water

    call flush(6)

  end subroutine cloud_prop_rrtm
  


  ! ****************************************************************************


  subroutine InitGetoz_rrtm(yrl,kmax,sl)
    real(kind=real64),    intent(in   ) :: yrl
    integer, intent(in   ) :: kmax
    real(kind=real64),    intent(in   ) :: sl (kmax)
    integer                :: l, ll


    if(.not. allocated(ozone)) allocate(ozone(nlm_getoz,nl,ns))

    !
    !     four season climatological ozone data in nmc sigma layers
    !
    !     for seasonal variation
    !     season=1 - winter          season=2 - spring
    !     season=3 - summer          season=4 - fall
    !     unit of ozone mixing ratio is in ( 10**-4 g/g ).  the data is
    !     in 18 sigma layers from top to bottom.  for every layer, there
    !     are 37 latitudes at 5 degree interval from north pole to south
    !     pole.
    !     mrf86 18 layers
    !
    !
    !     1. winter
    !
    !     wint1(18,6)
    !
    ozone(1:18, 1:6, 1) = reshape( (/ &
         .068467e0_real64,.052815e0_real64,.035175e0_real64,.022334e0_real64,.013676e0_real64,.007363e0_real64, &
         .003633e0_real64,.001582e0_real64,.001111e0_real64,.000713e0_real64,.000517e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .069523e0_real64,.052249e0_real64,.034255e0_real64,.021379e0_real64,.012306e0_real64,.006727e0_real64, &
         .003415e0_real64,.001578e0_real64,.001072e0_real64,.000681e0_real64,.000517e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .070579e0_real64,.051684e0_real64,.033335e0_real64,.020423e0_real64,.010935e0_real64,.006091e0_real64, &
         .003197e0_real64,.001573e0_real64,.001034e0_real64,.000650e0_real64,.000517e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .074885e0_real64,.049987e0_real64,.030140e0_real64,.017894e0_real64,.009881e0_real64,.005543e0_real64, &
         .002907e0_real64,.001379e0_real64,.000961e0_real64,.000644e0_real64,.000512e0_real64,.000463e0_real64, &
         .000451e0_real64,.000408e0_real64,.000385e0_real64,.000361e0_real64,.000351e0_real64,.000349e0_real64, &
         .079190e0_real64,.048290e0_real64,.026945e0_real64,.015366e0_real64,.008826e0_real64,.004995e0_real64, &
         .002616e0_real64,.001184e0_real64,.000887e0_real64,.000637e0_real64,.000508e0_real64,.000486e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .082443e0_real64,.047591e0_real64,.025358e0_real64,.014294e0_real64,.008233e0_real64,.004664e0_real64, &
         .002430e0_real64,.001068e0_real64,.000851e0_real64,.000644e0_real64,.000508e0_real64,.000474e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64/), &
         (/18,6/))
    !
    !     wint2(18,6)
    !
    ozone(1:18, 7:12, 1) = reshape( (/ &
         .085695e0_real64,.046892e0_real64,.023772e0_real64,.013223e0_real64,.007640e0_real64,.004333e0_real64, &
         .002244e0_real64,.000951e0_real64,.000815e0_real64,.000650e0_real64,.000508e0_real64,.000463e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .089618e0_real64,.042869e0_real64,.019963e0_real64,.010502e0_real64,.005966e0_real64,.003525e0_real64, &
         .001936e0_real64,.000906e0_real64,.000769e0_real64,.000625e0_real64,.000508e0_real64,.000452e0_real64, &
         .000451e0_real64,.000408e0_real64,.000385e0_real64,.000361e0_real64,.000351e0_real64,.000349e0_real64, &
         .093540e0_real64,.038846e0_real64,.016155e0_real64,.007781e0_real64,.004292e0_real64,.002716e0_real64, &
         .001628e0_real64,.000862e0_real64,.000724e0_real64,.000600e0_real64,.000508e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .097097e0_real64,.034916e0_real64,.012983e0_real64,.006240e0_real64,.003666e0_real64,.002259e0_real64, &
         .001336e0_real64,.000730e0_real64,.000629e0_real64,.000549e0_real64,.000499e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .100654e0_real64,.030986e0_real64,.009812e0_real64,.004698e0_real64,.003041e0_real64,.001803e0_real64, &
         .001044e0_real64,.000599e0_real64,.000533e0_real64,.000499e0_real64,.000491e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .101724e0_real64,.026500e0_real64,.007228e0_real64,.003391e0_real64,.002058e0_real64,.001285e0_real64, &
         .000811e0_real64,.000531e0_real64,.000478e0_real64,.000449e0_real64,.000440e0_real64,.000421e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64/), &
         (/18,6/))
    !
    !     wint3(18,6)
    !
    ozone(1:18, 13:18, 1) = reshape( (/ &
         .102794e0_real64,.022015e0_real64,.004645e0_real64,.002084e0_real64,.001076e0_real64,.000767e0_real64, &
         .000577e0_real64,.000463e0_real64,.000423e0_real64,.000399e0_real64,.000389e0_real64,.000401e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .103456e0_real64,.018235e0_real64,.003195e0_real64,.001379e0_real64,.000771e0_real64,.000585e0_real64, &
         .000474e0_real64,.000411e0_real64,.000380e0_real64,.000362e0_real64,.000343e0_real64,.000348e0_real64, &
         .000346e0_real64,.000328e0_real64,.000317e0_real64,.000305e0_real64,.000302e0_real64,.000302e0_real64, &
         .104118e0_real64,.014455e0_real64,.001745e0_real64,.000674e0_real64,.000467e0_real64,.000403e0_real64, &
         .000370e0_real64,.000359e0_real64,.000337e0_real64,.000325e0_real64,.000296e0_real64,.000294e0_real64, &
         .000293e0_real64,.000302e0_real64,.000306e0_real64,.000302e0_real64,.000302e0_real64,.000302e0_real64, &
         .104106e0_real64,.012997e0_real64,.001479e0_real64,.000639e0_real64,.000468e0_real64,.000422e0_real64, &
         .000392e0_real64,.000372e0_real64,.000342e0_real64,.000325e0_real64,.000296e0_real64,.000294e0_real64, &
         .000293e0_real64,.000302e0_real64,.000306e0_real64,.000302e0_real64,.000302e0_real64,.000302e0_real64, &
         .104093e0_real64,.011539e0_real64,.001213e0_real64,.000604e0_real64,.000468e0_real64,.000442e0_real64, &
         .000414e0_real64,.000385e0_real64,.000347e0_real64,.000325e0_real64,.000296e0_real64,.000294e0_real64, &
         .000293e0_real64,.000302e0_real64,.000306e0_real64,.000302e0_real64,.000302e0_real64,.000302e0_real64, &
         .104087e0_real64,.010726e0_real64,.000971e0_real64,.000538e0_real64,.000440e0_real64,.000434e0_real64, &
         .000418e0_real64,.000397e0_real64,.000375e0_real64,.000343e0_real64,.000296e0_real64,.000294e0_real64, &
         .000293e0_real64,.000302e0_real64,.000306e0_real64,.000302e0_real64,.000302e0_real64,.000302e0_real64/), &
         (/18,6/))
    !
    !     wint4(18,6)
    !
    ozone(1:18, 19:24, 1) = reshape( (/ &
         .102665e0_real64,.010977e0_real64,.001237e0_real64,.000590e0_real64,.000498e0_real64,.000479e0_real64, &
         .000458e0_real64,.000436e0_real64,.000421e0_real64,.000387e0_real64,.000326e0_real64,.000298e0_real64, &
         .000246e0_real64,.000227e0_real64,.000211e0_real64,.000200e0_real64,.000194e0_real64,.000186e0_real64, &
         .100892e0_real64,.012873e0_real64,.001886e0_real64,.000785e0_real64,.000643e0_real64,.000568e0_real64, &
         .000519e0_real64,.000487e0_real64,.000471e0_real64,.000437e0_real64,.000368e0_real64,.000305e0_real64, &
         .000201e0_real64,.000151e0_real64,.000117e0_real64,.000098e0_real64,.000090e0_real64,.000093e0_real64, &
         .100534e0_real64,.013704e0_real64,.002028e0_real64,.000861e0_real64,.000701e0_real64,.000604e0_real64, &
         .000546e0_real64,.000513e0_real64,.000504e0_real64,.000462e0_real64,.000381e0_real64,.000307e0_real64, &
         .000201e0_real64,.000151e0_real64,.000117e0_real64,.000098e0_real64,.000090e0_real64,.000093e0_real64, &
         .100218e0_real64,.015035e0_real64,.002537e0_real64,.001037e0_real64,.000790e0_real64,.000726e0_real64, &
         .000673e0_real64,.000628e0_real64,.000579e0_real64,.000512e0_real64,.000440e0_real64,.000374e0_real64, &
         .000307e0_real64,.000253e0_real64,.000227e0_real64,.000208e0_real64,.000194e0_real64,.000186e0_real64, &
         .099903e0_real64,.016365e0_real64,.003045e0_real64,.001214e0_real64,.000879e0_real64,.000848e0_real64, &
         .000801e0_real64,.000744e0_real64,.000654e0_real64,.000562e0_real64,.000499e0_real64,.000441e0_real64, &
         .000410e0_real64,.000358e0_real64,.000342e0_real64,.000322e0_real64,.000302e0_real64,.000302e0_real64, &
         .099547e0_real64,.017725e0_real64,.003693e0_real64,.001578e0_real64,.001125e0_real64,.000985e0_real64, &
         .000879e0_real64,.000795e0_real64,.000712e0_real64,.000643e0_real64,.000584e0_real64,.000521e0_real64, &
         .000482e0_real64,.000384e0_real64,.000351e0_real64,.000322e0_real64,.000302e0_real64,.000302e0_real64/), &
         (/18,6/))
    !
    !     wint5(18,6)
    !
    ozone(1:18, 25:30, 1) = reshape( (/ &
         .099191e0_real64,.019085e0_real64,.004340e0_real64,.001943e0_real64,.001371e0_real64,.001122e0_real64, &
         .000957e0_real64,.000847e0_real64,.000770e0_real64,.000724e0_real64,.000669e0_real64,.000601e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .098107e0_real64,.020617e0_real64,.004758e0_real64,.002137e0_real64,.001516e0_real64,.001211e0_real64, &
         .000999e0_real64,.000848e0_real64,.000778e0_real64,.000730e0_real64,.000677e0_real64,.000603e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .097023e0_real64,.022148e0_real64,.005177e0_real64,.002332e0_real64,.001660e0_real64,.001300e0_real64, &
         .001041e0_real64,.000849e0_real64,.000786e0_real64,.000737e0_real64,.000686e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .093464e0_real64,.026177e0_real64,.008525e0_real64,.003892e0_real64,.002452e0_real64,.001609e0_real64, &
         .001116e0_real64,.000851e0_real64,.000809e0_real64,.000762e0_real64,.000690e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .089906e0_real64,.030206e0_real64,.011873e0_real64,.005453e0_real64,.003244e0_real64,.001918e0_real64, &
         .001192e0_real64,.000852e0_real64,.000832e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .080939e0_real64,.032414e0_real64,.014163e0_real64,.007241e0_real64,.004328e0_real64,.002522e0_real64, &
         .001481e0_real64,.000934e0_real64,.000861e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64/), &
         (/18,6/))
    !
    !     wint6(18,6)
    !
    ozone(1:18, 31:36, 1) = reshape( (/ &
         .071972e0_real64,.034622e0_real64,.016453e0_real64,.009029e0_real64,.005413e0_real64,.003127e0_real64, &
         .001770e0_real64,.001015e0_real64,.000890e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .069820e0_real64,.035028e0_real64,.016929e0_real64,.009389e0_real64,.005645e0_real64,.003260e0_real64, &
         .001843e0_real64,.001055e0_real64,.000905e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .067669e0_real64,.035434e0_real64,.017406e0_real64,.009749e0_real64,.005876e0_real64,.003393e0_real64, &
         .001916e0_real64,.001094e0_real64,.000920e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .065518e0_real64,.035975e0_real64,.017854e0_real64,.010100e0_real64,.006534e0_real64,.003985e0_real64, &
         .002321e0_real64,.001240e0_real64,.000966e0_real64,.000774e0_real64,.000640e0_real64,.000548e0_real64, &
         .000479e0_real64,.000384e0_real64,.000346e0_real64,.000316e0_real64,.000302e0_real64,.000302e0_real64, &
         .063367e0_real64,.036516e0_real64,.018302e0_real64,.010452e0_real64,.007192e0_real64,.004577e0_real64, &
         .002727e0_real64,.001387e0_real64,.001012e0_real64,.000762e0_real64,.000585e0_real64,.000490e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .061216e0_real64,.037359e0_real64,.019151e0_real64,.010633e0_real64,.006845e0_real64,.004382e0_real64, &
         .002691e0_real64,.001511e0_real64,.001061e0_real64,.000749e0_real64,.000568e0_real64,.000465e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64/), &
         (/18,6/))
    !
    !     wint7(18)
    !
    ozone(1:18, 37, 1) = (/ &
         .059066e0_real64,.038201e0_real64,.019999e0_real64,.010813e0_real64,.006498e0_real64,.004188e0_real64, &
         .002656e0_real64,.001636e0_real64,.001110e0_real64,.000737e0_real64,.000551e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64/)
    !
    !     2. spring
    !
    ozone(1:18, 1:6, 2) = reshape( (/ &
         .074229e0_real64,.050084e0_real64,.030930e0_real64,.018676e0_real64,.011965e0_real64,.008165e0_real64, &
         .005428e0_real64,.003399e0_real64,.002098e0_real64,.001138e0_real64,.000780e0_real64,.000632e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64, &
         .074927e0_real64,.049459e0_real64,.029215e0_real64,.018025e0_real64,.011754e0_real64,.007786e0_real64, &
         .004972e0_real64,.002926e0_real64,.001817e0_real64,.001025e0_real64,.000758e0_real64,.000632e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64, &
         .075625e0_real64,.048835e0_real64,.027500e0_real64,.017375e0_real64,.011544e0_real64,.007407e0_real64, &
         .004516e0_real64,.002453e0_real64,.001536e0_real64,.000912e0_real64,.000737e0_real64,.000632e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64, &
         .077409e0_real64,.048159e0_real64,.026661e0_real64,.016596e0_real64,.010962e0_real64,.006972e0_real64, &
         .004160e0_real64,.002132e0_real64,.001391e0_real64,.000868e0_real64,.000686e0_real64,.000601e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64, &
         .079194e0_real64,.047483e0_real64,.025822e0_real64,.015818e0_real64,.010380e0_real64,.006537e0_real64, &
         .003804e0_real64,.001811e0_real64,.001245e0_real64,.000825e0_real64,.000635e0_real64,.000570e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64, &
         .084591e0_real64,.046553e0_real64,.025037e0_real64,.015156e0_real64,.009841e0_real64,.006124e0_real64, &
         .003534e0_real64,.001693e0_real64,.001170e0_real64,.000793e0_real64,.000631e0_real64,.000537e0_real64, &
         .000551e0_real64,.000509e0_real64,.000486e0_real64,.000516e0_real64,.000548e0_real64,.000446e0_real64/), &
         (/18,6/))
    ozone(1:18, 7:12, 2) = reshape( (/ &
         .089988e0_real64,.045622e0_real64,.024253e0_real64,.014495e0_real64,.009303e0_real64,.005711e0_real64, &
         .003264e0_real64,.001574e0_real64,.001096e0_real64,.000762e0_real64,.000627e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64, &
         .092863e0_real64,.042419e0_real64,.020704e0_real64,.012034e0_real64,.007417e0_real64,.004504e0_real64, &
         .002590e0_real64,.001334e0_real64,.000977e0_real64,.000731e0_real64,.000622e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64, &
         .095737e0_real64,.039215e0_real64,.017155e0_real64,.009572e0_real64,.005532e0_real64,.003296e0_real64, &
         .001916e0_real64,.001094e0_real64,.000858e0_real64,.000699e0_real64,.000618e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64, &
         .097501e0_real64,.035382e0_real64,.014856e0_real64,.008207e0_real64,.004619e0_real64,.002720e0_real64, &
         .001610e0_real64,.001012e0_real64,.000829e0_real64,.000687e0_real64,.000610e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64, &
         .099264e0_real64,.031548e0_real64,.012557e0_real64,.006841e0_real64,.003705e0_real64,.002144e0_real64, &
         .001304e0_real64,.000930e0_real64,.000799e0_real64,.000675e0_real64,.000601e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64, &
         .101718e0_real64,.026523e0_real64,.008473e0_real64,.004382e0_real64,.002392e0_real64,.001505e0_real64, &
         .001036e0_real64,.000836e0_real64,.000727e0_real64,.000618e0_real64,.000550e0_real64,.000494e0_real64, &
         .000501e0_real64,.000479e0_real64,.000473e0_real64,.000509e0_real64,.000541e0_real64,.000445e0_real64/), &
         (/18,6/))
    ozone(1:18, 13:18, 2) = reshape( (/ &
         .104172e0_real64,.021499e0_real64,.004389e0_real64,.001922e0_real64,.001078e0_real64,.000865e0_real64, &
         .000767e0_real64,.000743e0_real64,.000654e0_real64,.000562e0_real64,.000499e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .104145e0_real64,.018082e0_real64,.003274e0_real64,.001493e0_real64,.000919e0_real64,.000762e0_real64, &
         .000678e0_real64,.000641e0_real64,.000584e0_real64,.000531e0_real64,.000495e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .104118e0_real64,.014665e0_real64,.002159e0_real64,.001063e0_real64,.000759e0_real64,.000659e0_real64, &
         .000589e0_real64,.000539e0_real64,.000514e0_real64,.000499e0_real64,.000491e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .107719e0_real64,.013052e0_real64,.001822e0_real64,.000953e0_real64,.000701e0_real64,.000604e0_real64, &
         .000551e0_real64,.000525e0_real64,.000509e0_real64,.000499e0_real64,.000491e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .111320e0_real64,.011439e0_real64,.001485e0_real64,.000843e0_real64,.000642e0_real64,.000549e0_real64, &
         .000512e0_real64,.000512e0_real64,.000504e0_real64,.000499e0_real64,.000491e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .112375e0_real64,.011255e0_real64,.001357e0_real64,.000744e0_real64,.000585e0_real64,.000533e0_real64, &
         .000512e0_real64,.000512e0_real64,.000504e0_real64,.000499e0_real64,.000491e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64/), &
         (/18,6/))
    ozone(1:18, 19:24, 2) = reshape( (/ &
         .109850e0_real64,.010424e0_real64,.001079e0_real64,.000567e0_real64,.000498e0_real64,.000479e0_real64, &
         .000463e0_real64,.000448e0_real64,.000418e0_real64,.000399e0_real64,.000389e0_real64,.000367e0_real64, &
         .000351e0_real64,.000328e0_real64,.000320e0_real64,.000337e0_real64,.000355e0_real64,.000304e0_real64, &
         .107002e0_real64,.009961e0_real64,.001025e0_real64,.000533e0_real64,.000497e0_real64,.000460e0_real64, &
         .000422e0_real64,.000385e0_real64,.000332e0_real64,.000300e0_real64,.000288e0_real64,.000249e0_real64, &
         .000202e0_real64,.000158e0_real64,.000132e0_real64,.000114e0_real64,.000104e0_real64,.000093e0_real64, &
         .107735e0_real64,.010146e0_real64,.001120e0_real64,.000576e0_real64,.000526e0_real64,.000477e0_real64, &
         .000430e0_real64,.000385e0_real64,.000332e0_real64,.000300e0_real64,.000288e0_real64,.000249e0_real64, &
         .000202e0_real64,.000158e0_real64,.000132e0_real64,.000114e0_real64,.000104e0_real64,.000093e0_real64, &
         .107021e0_real64,.012233e0_real64,.001533e0_real64,.000643e0_real64,.000556e0_real64,.000505e0_real64, &
         .000471e0_real64,.000448e0_real64,.000403e0_real64,.000362e0_real64,.000355e0_real64,.000296e0_real64, &
         .000251e0_real64,.000207e0_real64,.000180e0_real64,.000161e0_real64,.000152e0_real64,.000140e0_real64, &
         .106308e0_real64,.014320e0_real64,.001946e0_real64,.000709e0_real64,.000585e0_real64,.000533e0_real64, &
         .000512e0_real64,.000512e0_real64,.000473e0_real64,.000425e0_real64,.000423e0_real64,.000342e0_real64, &
         .000301e0_real64,.000257e0_real64,.000232e0_real64,.000212e0_real64,.000205e0_real64,.000209e0_real64, &
         .100592e0_real64,.015718e0_real64,.002411e0_real64,.001007e0_real64,.000802e0_real64,.000642e0_real64, &
         .000559e0_real64,.000526e0_real64,.000501e0_real64,.000474e0_real64,.000470e0_real64,.000439e0_real64, &
         .000430e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64/), &
         (/18,6/))
    ozone(1:18, 25:30, 2) = reshape( (/ &
         .094877e0_real64,.017116e0_real64,.002877e0_real64,.001305e0_real64,.001018e0_real64,.000751e0_real64, &
         .000606e0_real64,.000539e0_real64,.000529e0_real64,.000524e0_real64,.000516e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .094163e0_real64,.020198e0_real64,.004594e0_real64,.001772e0_real64,.001077e0_real64,.000806e0_real64, &
         .000649e0_real64,.000565e0_real64,.000547e0_real64,.000537e0_real64,.000521e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .093449e0_real64,.023279e0_real64,.006312e0_real64,.002240e0_real64,.001135e0_real64,.000862e0_real64, &
         .000692e0_real64,.000591e0_real64,.000564e0_real64,.000549e0_real64,.000525e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .089886e0_real64,.026029e0_real64,.008558e0_real64,.003312e0_real64,.001655e0_real64,.001124e0_real64, &
         .000807e0_real64,.000631e0_real64,.000602e0_real64,.000568e0_real64,.000525e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .086323e0_real64,.028778e0_real64,.010805e0_real64,.004383e0_real64,.002175e0_real64,.001386e0_real64, &
         .000923e0_real64,.000671e0_real64,.000640e0_real64,.000587e0_real64,.000525e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .082715e0_real64,.031096e0_real64,.013350e0_real64,.006131e0_real64,.003205e0_real64,.002043e0_real64, &
         .001304e0_real64,.000842e0_real64,.000734e0_real64,.000631e0_real64,.000555e0_real64,.000494e0_real64, &
         .000480e0_real64,.000408e0_real64,.000385e0_real64,.000361e0_real64,.000351e0_real64,.000349e0_real64/), &
         (/18,6/))
    ozone(1:18, 31:36, 2) = reshape( (/ &
         .079108e0_real64,.033415e0_real64,.015895e0_real64,.007878e0_real64,.004234e0_real64,.002700e0_real64, &
         .001686e0_real64,.001014e0_real64,.000829e0_real64,.000675e0_real64,.000584e0_real64,.000454e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .074807e0_real64,.034651e0_real64,.017056e0_real64,.008574e0_real64,.004769e0_real64,.002986e0_real64, &
         .001827e0_real64,.001079e0_real64,.000853e0_real64,.000675e0_real64,.000584e0_real64,.000454e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .070506e0_real64,.035887e0_real64,.018218e0_real64,.009270e0_real64,.005304e0_real64,.003271e0_real64, &
         .001969e0_real64,.001145e0_real64,.000878e0_real64,.000675e0_real64,.000584e0_real64,.000454e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .067669e0_real64,.037799e0_real64,.019680e0_real64,.009612e0_real64,.005481e0_real64,.003476e0_real64, &
         .002093e0_real64,.001123e0_real64,.000837e0_real64,.000631e0_real64,.000546e0_real64,.000447e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .064832e0_real64,.039712e0_real64,.021142e0_real64,.009954e0_real64,.005658e0_real64,.003681e0_real64, &
         .002218e0_real64,.001100e0_real64,.000796e0_real64,.000587e0_real64,.000508e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .063734e0_real64,.039842e0_real64,.022004e0_real64,.010859e0_real64,.005712e0_real64,.003589e0_real64, &
         .002155e0_real64,.001174e0_real64,.000856e0_real64,.000612e0_real64,.000508e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64/), &
         (/18,6/))
    ozone(1:18, 37, 2) = (/ &
         .062636e0_real64,.039972e0_real64,.022867e0_real64,.011765e0_real64,.005766e0_real64,.003498e0_real64, &
         .002092e0_real64,.001248e0_real64,.000917e0_real64,.000637e0_real64,.000508e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64/)
    !
    !     3. summer
    !
    ozone(1:18, 1:6, 3) = reshape( (/ &
         .059066e0_real64,.038201e0_real64,.019999e0_real64,.010813e0_real64,.006498e0_real64,.004188e0_real64, &
         .002656e0_real64,.001636e0_real64,.001110e0_real64,.000737e0_real64,.000551e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .061216e0_real64,.037359e0_real64,.019151e0_real64,.010633e0_real64,.006845e0_real64,.004382e0_real64, &
         .002691e0_real64,.001511e0_real64,.001061e0_real64,.000749e0_real64,.000568e0_real64,.000465e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .063367e0_real64,.036516e0_real64,.018302e0_real64,.010452e0_real64,.007192e0_real64,.004577e0_real64, &
         .002727e0_real64,.001387e0_real64,.001012e0_real64,.000762e0_real64,.000585e0_real64,.000490e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .065518e0_real64,.035975e0_real64,.017854e0_real64,.010100e0_real64,.006534e0_real64,.003985e0_real64, &
         .002321e0_real64,.001240e0_real64,.000966e0_real64,.000774e0_real64,.000640e0_real64,.000548e0_real64, &
         .000479e0_real64,.000384e0_real64,.000346e0_real64,.000316e0_real64,.000302e0_real64,.000302e0_real64, &
         .067669e0_real64,.035434e0_real64,.017406e0_real64,.009749e0_real64,.005876e0_real64,.003393e0_real64, &
         .001916e0_real64,.001094e0_real64,.000920e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .069820e0_real64,.035028e0_real64,.016929e0_real64,.009389e0_real64,.005645e0_real64,.003260e0_real64, &
         .001843e0_real64,.001055e0_real64,.000905e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64/), &
         (/18,6/))
    ozone(1:18, 7:12, 3) = reshape( (/ &
         .071972e0_real64,.034622e0_real64,.016453e0_real64,.009029e0_real64,.005413e0_real64,.003127e0_real64, &
         .001770e0_real64,.001015e0_real64,.000890e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .080939e0_real64,.032414e0_real64,.014163e0_real64,.007241e0_real64,.004328e0_real64,.002522e0_real64, &
         .001481e0_real64,.000934e0_real64,.000861e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .089906e0_real64,.030206e0_real64,.011873e0_real64,.005453e0_real64,.003244e0_real64,.001918e0_real64, &
         .001192e0_real64,.000852e0_real64,.000832e0_real64,.000787e0_real64,.000694e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .093464e0_real64,.026177e0_real64,.008525e0_real64,.003892e0_real64,.002452e0_real64,.001609e0_real64, &
         .001116e0_real64,.000851e0_real64,.000809e0_real64,.000762e0_real64,.000690e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .097023e0_real64,.022148e0_real64,.005177e0_real64,.002332e0_real64,.001660e0_real64,.001300e0_real64, &
         .001041e0_real64,.000849e0_real64,.000786e0_real64,.000737e0_real64,.000686e0_real64,.000606e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .098107e0_real64,.020617e0_real64,.004758e0_real64,.002137e0_real64,.001516e0_real64,.001211e0_real64, &
         .000999e0_real64,.000848e0_real64,.000778e0_real64,.000730e0_real64,.000677e0_real64,.000603e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64/), &
         (/18,6/))
    ozone(1:18, 13:18, 3) = reshape( (/ &
         .099191e0_real64,.019085e0_real64,.004340e0_real64,.001943e0_real64,.001371e0_real64,.001122e0_real64, &
         .000957e0_real64,.000847e0_real64,.000770e0_real64,.000724e0_real64,.000669e0_real64,.000601e0_real64, &
         .000557e0_real64,.000412e0_real64,.000362e0_real64,.000326e0_real64,.000309e0_real64,.000302e0_real64, &
         .099547e0_real64,.017725e0_real64,.003693e0_real64,.001578e0_real64,.001125e0_real64,.000985e0_real64, &
         .000879e0_real64,.000795e0_real64,.000712e0_real64,.000643e0_real64,.000584e0_real64,.000521e0_real64, &
         .000482e0_real64,.000384e0_real64,.000351e0_real64,.000322e0_real64,.000302e0_real64,.000302e0_real64, &
         .099903e0_real64,.016365e0_real64,.003045e0_real64,.001214e0_real64,.000879e0_real64,.000848e0_real64, &
         .000801e0_real64,.000744e0_real64,.000654e0_real64,.000562e0_real64,.000499e0_real64,.000441e0_real64, &
         .000410e0_real64,.000358e0_real64,.000342e0_real64,.000322e0_real64,.000302e0_real64,.000302e0_real64, &
         .100218e0_real64,.015035e0_real64,.002537e0_real64,.001037e0_real64,.000790e0_real64,.000726e0_real64, &
         .000673e0_real64,.000628e0_real64,.000579e0_real64,.000512e0_real64,.000440e0_real64,.000374e0_real64, &
         .000307e0_real64,.000253e0_real64,.000227e0_real64,.000208e0_real64,.000194e0_real64,.000186e0_real64, &
         .100534e0_real64,.013704e0_real64,.002028e0_real64,.000861e0_real64,.000701e0_real64,.000604e0_real64, &
         .000546e0_real64,.000513e0_real64,.000504e0_real64,.000462e0_real64,.000381e0_real64,.000307e0_real64, &
         .000201e0_real64,.000151e0_real64,.000117e0_real64,.000098e0_real64,.000090e0_real64,.000093e0_real64, &
         .100892e0_real64,.012873e0_real64,.001886e0_real64,.000785e0_real64,.000643e0_real64,.000568e0_real64, &
         .000519e0_real64,.000487e0_real64,.000471e0_real64,.000437e0_real64,.000368e0_real64,.000305e0_real64, &
         .000201e0_real64,.000151e0_real64,.000117e0_real64,.000098e0_real64,.000090e0_real64,.000093e0_real64/), &
         (/18,6/))
    ozone(1:18, 19:24, 3) = reshape( (/ &
         .102665e0_real64,.010977e0_real64,.001237e0_real64,.000590e0_real64,.000498e0_real64,.000479e0_real64, &
         .000458e0_real64,.000436e0_real64,.000421e0_real64,.000387e0_real64,.000326e0_real64,.000298e0_real64, &
         .000246e0_real64,.000227e0_real64,.000211e0_real64,.000200e0_real64,.000194e0_real64,.000186e0_real64, &
         .104087e0_real64,.010726e0_real64,.000971e0_real64,.000538e0_real64,.000440e0_real64,.000434e0_real64, &
         .000418e0_real64,.000397e0_real64,.000375e0_real64,.000343e0_real64,.000296e0_real64,.000294e0_real64, &
         .000293e0_real64,.000302e0_real64,.000306e0_real64,.000302e0_real64,.000302e0_real64,.000302e0_real64, &
         .104093e0_real64,.011539e0_real64,.001213e0_real64,.000604e0_real64,.000468e0_real64,.000442e0_real64, &
         .000414e0_real64,.000385e0_real64,.000347e0_real64,.000325e0_real64,.000296e0_real64,.000294e0_real64, &
         .000293e0_real64,.000302e0_real64,.000306e0_real64,.000302e0_real64,.000302e0_real64,.000302e0_real64, &
         .104106e0_real64,.012997e0_real64,.001479e0_real64,.000639e0_real64,.000468e0_real64,.000422e0_real64, &
         .000392e0_real64,.000372e0_real64,.000342e0_real64,.000325e0_real64,.000296e0_real64,.000294e0_real64, &
         .000293e0_real64,.000302e0_real64,.000306e0_real64,.000302e0_real64,.000302e0_real64,.000302e0_real64, &
         .104118e0_real64,.014455e0_real64,.001745e0_real64,.000674e0_real64,.000467e0_real64,.000403e0_real64, &
         .000370e0_real64,.000359e0_real64,.000337e0_real64,.000325e0_real64,.000296e0_real64,.000294e0_real64, &
         .000293e0_real64,.000302e0_real64,.000306e0_real64,.000302e0_real64,.000302e0_real64,.000302e0_real64, &
         .103456e0_real64,.018235e0_real64,.003195e0_real64,.001379e0_real64,.000771e0_real64,.000585e0_real64, &
         .000474e0_real64,.000411e0_real64,.000380e0_real64,.000362e0_real64,.000343e0_real64,.000348e0_real64, &
         .000346e0_real64,.000328e0_real64,.000317e0_real64,.000305e0_real64,.000302e0_real64,.000302e0_real64/), &
         (/18,6/))
    ozone(1:18, 25:30, 3) = reshape( (/ &
         .102794e0_real64,.022015e0_real64,.004645e0_real64,.002084e0_real64,.001076e0_real64,.000767e0_real64, &
         .000577e0_real64,.000463e0_real64,.000423e0_real64,.000399e0_real64,.000389e0_real64,.000401e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .101724e0_real64,.026500e0_real64,.007228e0_real64,.003391e0_real64,.002058e0_real64,.001285e0_real64, &
         .000811e0_real64,.000531e0_real64,.000478e0_real64,.000449e0_real64,.000440e0_real64,.000421e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .100654e0_real64,.030986e0_real64,.009812e0_real64,.004698e0_real64,.003041e0_real64,.001803e0_real64, &
         .001044e0_real64,.000599e0_real64,.000533e0_real64,.000499e0_real64,.000491e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .097097e0_real64,.034916e0_real64,.012983e0_real64,.006240e0_real64,.003666e0_real64,.002259e0_real64, &
         .001336e0_real64,.000730e0_real64,.000629e0_real64,.000549e0_real64,.000499e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .093540e0_real64,.038846e0_real64,.016155e0_real64,.007781e0_real64,.004292e0_real64,.002716e0_real64, &
         .001628e0_real64,.000862e0_real64,.000724e0_real64,.000600e0_real64,.000508e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .089618e0_real64,.042869e0_real64,.019963e0_real64,.010502e0_real64,.005966e0_real64,.003525e0_real64, &
         .001936e0_real64,.000906e0_real64,.000769e0_real64,.000625e0_real64,.000508e0_real64,.000452e0_real64, &
         .000451e0_real64,.000408e0_real64,.000385e0_real64,.000361e0_real64,.000351e0_real64,.000349e0_real64/), &
         (/18,6/))
    ozone(1:18, 31:36, 3) = reshape( (/ &
         .085695e0_real64,.046892e0_real64,.023772e0_real64,.013223e0_real64,.007640e0_real64,.004333e0_real64, &
         .002244e0_real64,.000951e0_real64,.000815e0_real64,.000650e0_real64,.000508e0_real64,.000463e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .082443e0_real64,.047591e0_real64,.025358e0_real64,.014294e0_real64,.008233e0_real64,.004664e0_real64, &
         .002430e0_real64,.001068e0_real64,.000851e0_real64,.000644e0_real64,.000508e0_real64,.000474e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .079190e0_real64,.048290e0_real64,.026945e0_real64,.015366e0_real64,.008826e0_real64,.004995e0_real64, &
         .002616e0_real64,.001184e0_real64,.000887e0_real64,.000637e0_real64,.000508e0_real64,.000486e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .074885e0_real64,.049987e0_real64,.030140e0_real64,.017894e0_real64,.009881e0_real64,.005543e0_real64, &
         .002907e0_real64,.001379e0_real64,.000961e0_real64,.000644e0_real64,.000512e0_real64,.000463e0_real64, &
         .000451e0_real64,.000408e0_real64,.000385e0_real64,.000361e0_real64,.000351e0_real64,.000349e0_real64, &
         .070579e0_real64,.051684e0_real64,.033335e0_real64,.020423e0_real64,.010935e0_real64,.006091e0_real64, &
         .003197e0_real64,.001573e0_real64,.001034e0_real64,.000650e0_real64,.000517e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .069523e0_real64,.052249e0_real64,.034255e0_real64,.021379e0_real64,.012306e0_real64,.006727e0_real64, &
         .003415e0_real64,.001578e0_real64,.001072e0_real64,.000681e0_real64,.000517e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64/), &
         (/18,6/))
    ozone(1:18, 37, 3) = (/ &
         .068467e0_real64,.052815e0_real64,.035175e0_real64,.022334e0_real64,.013676e0_real64,.007363e0_real64, &
         .003633e0_real64,.001582e0_real64,.001111e0_real64,.000713e0_real64,.000517e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64/)
    !
    !     4. fall
    !
    ozone(1:18, 1:6, 4) = reshape( (/ &
         .062636e0_real64,.039972e0_real64,.022867e0_real64,.011765e0_real64,.005766e0_real64,.003498e0_real64, &
         .002092e0_real64,.001248e0_real64,.000917e0_real64,.000637e0_real64,.000508e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .063734e0_real64,.039842e0_real64,.022004e0_real64,.010859e0_real64,.005712e0_real64,.003589e0_real64, &
         .002155e0_real64,.001174e0_real64,.000856e0_real64,.000612e0_real64,.000508e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .064832e0_real64,.039712e0_real64,.021142e0_real64,.009954e0_real64,.005658e0_real64,.003681e0_real64, &
         .002218e0_real64,.001100e0_real64,.000796e0_real64,.000587e0_real64,.000508e0_real64,.000441e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .067669e0_real64,.037799e0_real64,.019680e0_real64,.009612e0_real64,.005481e0_real64,.003476e0_real64, &
         .002093e0_real64,.001123e0_real64,.000837e0_real64,.000631e0_real64,.000546e0_real64,.000447e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .070506e0_real64,.035887e0_real64,.018218e0_real64,.009270e0_real64,.005304e0_real64,.003271e0_real64, &
         .001969e0_real64,.001145e0_real64,.000878e0_real64,.000675e0_real64,.000584e0_real64,.000454e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .074807e0_real64,.034651e0_real64,.017056e0_real64,.008574e0_real64,.004769e0_real64,.002986e0_real64, &
         .001827e0_real64,.001079e0_real64,.000853e0_real64,.000675e0_real64,.000584e0_real64,.000454e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64/), &
         (/18,6/))
    ozone(1:18, 7:12, 4) = reshape( (/ &
         .079108e0_real64,.033415e0_real64,.015895e0_real64,.007878e0_real64,.004234e0_real64,.002700e0_real64, &
         .001686e0_real64,.001014e0_real64,.000829e0_real64,.000675e0_real64,.000584e0_real64,.000454e0_real64, &
         .000401e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .082715e0_real64,.031096e0_real64,.013350e0_real64,.006131e0_real64,.003205e0_real64,.002043e0_real64, &
         .001304e0_real64,.000842e0_real64,.000734e0_real64,.000631e0_real64,.000555e0_real64,.000494e0_real64, &
         .000480e0_real64,.000408e0_real64,.000385e0_real64,.000361e0_real64,.000351e0_real64,.000349e0_real64, &
         .086323e0_real64,.028778e0_real64,.010805e0_real64,.004383e0_real64,.002175e0_real64,.001386e0_real64, &
         .000923e0_real64,.000671e0_real64,.000640e0_real64,.000587e0_real64,.000525e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .089886e0_real64,.026029e0_real64,.008558e0_real64,.003312e0_real64,.001655e0_real64,.001124e0_real64, &
         .000807e0_real64,.000631e0_real64,.000602e0_real64,.000568e0_real64,.000525e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .093449e0_real64,.023279e0_real64,.006312e0_real64,.002240e0_real64,.001135e0_real64,.000862e0_real64, &
         .000692e0_real64,.000591e0_real64,.000564e0_real64,.000549e0_real64,.000525e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .094163e0_real64,.020198e0_real64,.004594e0_real64,.001772e0_real64,.001077e0_real64,.000806e0_real64, &
         .000649e0_real64,.000565e0_real64,.000547e0_real64,.000537e0_real64,.000521e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64/), &
         (/18,6/))
    ozone(1:18, 13:18, 4) = reshape( (/ &
         .094877e0_real64,.017116e0_real64,.002877e0_real64,.001305e0_real64,.001018e0_real64,.000751e0_real64, &
         .000606e0_real64,.000539e0_real64,.000529e0_real64,.000524e0_real64,.000516e0_real64,.000535e0_real64, &
         .000558e0_real64,.000459e0_real64,.000436e0_real64,.000416e0_real64,.000406e0_real64,.000395e0_real64, &
         .100592e0_real64,.015718e0_real64,.002411e0_real64,.001007e0_real64,.000802e0_real64,.000642e0_real64, &
         .000559e0_real64,.000526e0_real64,.000501e0_real64,.000474e0_real64,.000470e0_real64,.000439e0_real64, &
         .000430e0_real64,.000358e0_real64,.000333e0_real64,.000311e0_real64,.000302e0_real64,.000302e0_real64, &
         .106308e0_real64,.014320e0_real64,.001946e0_real64,.000709e0_real64,.000585e0_real64,.000533e0_real64, &
         .000512e0_real64,.000512e0_real64,.000473e0_real64,.000425e0_real64,.000423e0_real64,.000342e0_real64, &
         .000301e0_real64,.000257e0_real64,.000232e0_real64,.000212e0_real64,.000205e0_real64,.000209e0_real64, &
         .107021e0_real64,.012233e0_real64,.001533e0_real64,.000643e0_real64,.000556e0_real64,.000505e0_real64, &
         .000471e0_real64,.000448e0_real64,.000403e0_real64,.000362e0_real64,.000355e0_real64,.000296e0_real64, &
         .000251e0_real64,.000207e0_real64,.000180e0_real64,.000161e0_real64,.000152e0_real64,.000140e0_real64, &
         .107735e0_real64,.010146e0_real64,.001120e0_real64,.000576e0_real64,.000526e0_real64,.000477e0_real64, &
         .000430e0_real64,.000385e0_real64,.000332e0_real64,.000300e0_real64,.000288e0_real64,.000249e0_real64, &
         .000202e0_real64,.000158e0_real64,.000132e0_real64,.000114e0_real64,.000104e0_real64,.000093e0_real64, &
         .107002e0_real64,.009961e0_real64,.001025e0_real64,.000533e0_real64,.000497e0_real64,.000460e0_real64, &
         .000422e0_real64,.000385e0_real64,.000332e0_real64,.000300e0_real64,.000288e0_real64,.000249e0_real64, &
         .000202e0_real64,.000158e0_real64,.000132e0_real64,.000114e0_real64,.000104e0_real64,.000093e0_real64/), &
         (/18,6/))
    ozone(1:18, 19:24, 4) = reshape( (/ &
         .109850e0_real64,.010424e0_real64,.001079e0_real64,.000567e0_real64,.000498e0_real64,.000479e0_real64, &
         .000463e0_real64,.000448e0_real64,.000418e0_real64,.000399e0_real64,.000389e0_real64,.000367e0_real64, &
         .000351e0_real64,.000328e0_real64,.000320e0_real64,.000337e0_real64,.000355e0_real64,.000304e0_real64, &
         .112375e0_real64,.011255e0_real64,.001357e0_real64,.000744e0_real64,.000585e0_real64,.000533e0_real64, &
         .000512e0_real64,.000512e0_real64,.000504e0_real64,.000499e0_real64,.000491e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .111320e0_real64,.011439e0_real64,.001485e0_real64,.000843e0_real64,.000642e0_real64,.000549e0_real64, &
         .000512e0_real64,.000512e0_real64,.000504e0_real64,.000499e0_real64,.000491e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .107719e0_real64,.013052e0_real64,.001822e0_real64,.000953e0_real64,.000701e0_real64,.000604e0_real64, &
         .000551e0_real64,.000525e0_real64,.000509e0_real64,.000499e0_real64,.000491e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .104118e0_real64,.014665e0_real64,.002159e0_real64,.001063e0_real64,.000759e0_real64,.000659e0_real64, &
         .000589e0_real64,.000539e0_real64,.000514e0_real64,.000499e0_real64,.000491e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .104145e0_real64,.018082e0_real64,.003274e0_real64,.001493e0_real64,.000919e0_real64,.000762e0_real64, &
         .000678e0_real64,.000641e0_real64,.000584e0_real64,.000531e0_real64,.000495e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64/), &
         (/18,6/))
    ozone(1:18, 25:30, 4) = reshape( (/ &
         .104172e0_real64,.021499e0_real64,.004389e0_real64,.001922e0_real64,.001078e0_real64,.000865e0_real64, &
         .000767e0_real64,.000743e0_real64,.000654e0_real64,.000562e0_real64,.000499e0_real64,.000486e0_real64, &
         .000501e0_real64,.000502e0_real64,.000509e0_real64,.000561e0_real64,.000607e0_real64,.000515e0_real64, &
         .101718e0_real64,.026523e0_real64,.008473e0_real64,.004382e0_real64,.002392e0_real64,.001505e0_real64, &
         .001036e0_real64,.000836e0_real64,.000727e0_real64,.000618e0_real64,.000550e0_real64,.000494e0_real64, &
         .000501e0_real64,.000479e0_real64,.000473e0_real64,.000509e0_real64,.000541e0_real64,.000445e0_real64, &
         .099264e0_real64,.031548e0_real64,.012557e0_real64,.006841e0_real64,.003705e0_real64,.002144e0_real64, &
         .001304e0_real64,.000930e0_real64,.000799e0_real64,.000675e0_real64,.000601e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64, &
         .097501e0_real64,.035382e0_real64,.014856e0_real64,.008207e0_real64,.004619e0_real64,.002720e0_real64, &
         .001610e0_real64,.001012e0_real64,.000829e0_real64,.000687e0_real64,.000610e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64, &
         .095737e0_real64,.039215e0_real64,.017155e0_real64,.009572e0_real64,.005532e0_real64,.003296e0_real64, &
         .001916e0_real64,.001094e0_real64,.000858e0_real64,.000699e0_real64,.000618e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64, &
         .092863e0_real64,.042419e0_real64,.020704e0_real64,.012034e0_real64,.007417e0_real64,.004504e0_real64, &
         .002590e0_real64,.001334e0_real64,.000977e0_real64,.000731e0_real64,.000622e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64/), &
         (/18,6/))
    ozone(1:18, 31:36, 4) = reshape( (/ &
         .089988e0_real64,.045622e0_real64,.024253e0_real64,.014495e0_real64,.009303e0_real64,.005711e0_real64, &
         .003264e0_real64,.001574e0_real64,.001096e0_real64,.000762e0_real64,.000627e0_real64,.000503e0_real64, &
         .000501e0_real64,.000459e0_real64,.000436e0_real64,.000460e0_real64,.000486e0_real64,.000398e0_real64, &
         .084591e0_real64,.046553e0_real64,.025037e0_real64,.015156e0_real64,.009841e0_real64,.006124e0_real64, &
         .003534e0_real64,.001693e0_real64,.001170e0_real64,.000793e0_real64,.000631e0_real64,.000537e0_real64, &
         .000551e0_real64,.000509e0_real64,.000486e0_real64,.000516e0_real64,.000548e0_real64,.000446e0_real64, &
         .079194e0_real64,.047483e0_real64,.025822e0_real64,.015818e0_real64,.010380e0_real64,.006537e0_real64, &
         .003804e0_real64,.001811e0_real64,.001245e0_real64,.000825e0_real64,.000635e0_real64,.000570e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64, &
         .077409e0_real64,.048159e0_real64,.026661e0_real64,.016596e0_real64,.010962e0_real64,.006972e0_real64, &
         .004160e0_real64,.002132e0_real64,.001391e0_real64,.000868e0_real64,.000686e0_real64,.000601e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64, &
         .075625e0_real64,.048835e0_real64,.027500e0_real64,.017375e0_real64,.011544e0_real64,.007407e0_real64, &
         .004516e0_real64,.002453e0_real64,.001536e0_real64,.000912e0_real64,.000737e0_real64,.000632e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64, &
         .074927e0_real64,.049459e0_real64,.029215e0_real64,.018025e0_real64,.011754e0_real64,.007786e0_real64, &
         .004972e0_real64,.002926e0_real64,.001817e0_real64,.001025e0_real64,.000758e0_real64,.000632e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64/), &
         (/18,6/))
    ozone(1:18, 37, 4) = (/ &
         .074229e0_real64,.050084e0_real64,.030930e0_real64,.018676e0_real64,.011965e0_real64,.008165e0_real64, &
         .005428e0_real64,.003399e0_real64,.002098e0_real64,.001138e0_real64,.000780e0_real64,.000632e0_real64, &
         .000603e0_real64,.000559e0_real64,.000538e0_real64,.000574e0_real64,.000614e0_real64,.000515e0_real64/)

    ozsig(:) = (/ &
         .020747_real64,.073986_real64,.124402_real64,.174576_real64,.224668_real64,.274735_real64, &
         .324767_real64,.374806_real64,.424818_real64,.497450_real64,.593540_real64,.688125_real64, &
         .777224_real64,.856317_real64,.920400_real64,.960480_real64,.981488_real64,.995004_real64/)


    first_getoz=.true.

    if(first_getoz)then
       mon_getoz=yrl/12.0_real64
       year_getoz=yrl
       if(nlm_getoz.ne.kmax)then
          inter_getoz=.true.
       else
          inter_getoz=.false.
          do l=1,nlm_getoz
             ll=nlm_getoz-l+1
             if(abs(ozsig(l)-sl(ll)).gt.0.0001_real64) inter_getoz=.true.
          end do
       endif
       first_getoz=.false.
    endif
  end subroutine InitGetoz_rrtm

  subroutine getoz_rrtm (kmax,sigmid,colrad,date,o3l)
    !
    ! input parameters and variables:
    !     ncols  =  number of atmospheric columns
    !     kmax   =  number of atmospheric layers
    !     sigmid =  sigma coordinate at middle of layer
    !     colrad =  colatitude of each column (0-3.14 from np to sp in radians)
    !     date   =  model julian date
    !
    ! tabulated data
    !     ozone  =  climatological ozone mixing ratio in 18 sigma layers
    !               and in 5 degree latitude interval
    !
    ! output variables:
    !     o3l   =  18 layers ozone mixing ratio in given lat and date
    !
    !==========================================================================
    ! :: kmax.....Number of grid points at vertical
    ! :: sigmid.......sigma coordinate at middle of layer
    ! :: pai......constant pi=3.1415926
    ! :: yrl......length of year in days
    !==========================================================================

    ! Input variables
    integer      ,    intent(IN   ) :: kmax
    real(KIND=real64),    intent(IN   ) :: sigmid (kmax)
    real(KIND=real64),    intent(IN   ) :: colrad
    real(KIND=real64),    intent(INOUT) :: date
    real(KIND=real64),parameter :: pai=3.1415926e00_real64

    ! Output variable
    real(KIND=real64),    intent(OUT) :: o3l(kmax)

    ! Local variables
    real(KIND=real64) :: a1   (nlm_getoz)
    real(KIND=real64) :: a2   (nlm_getoz)
    real(KIND=real64) :: a3   (nlm_getoz)
    real(KIND=real64) :: a4   (nlm_getoz)
    real(KIND=real64) :: b1   (nlm_getoz)
    real(KIND=real64) :: b2   (nlm_getoz)
    real(KIND=real64) :: b3   (nlm_getoz)
    real(KIND=real64) :: b4   (nlm_getoz)
    real(KIND=real64) :: do3a (nlm_getoz)
    real(KIND=real64) :: do3b (nlm_getoz)
    real(KIND=real64) :: ozo3l(nlm_getoz)

    real(KIND=real64), parameter :: rlag = 14.8125e0_real64

    integer :: l
    integer :: la
    integer :: ll
    integer :: kmx
    integer :: imon
    integer :: isea
    integer :: k
    integer :: i
    integer :: kk    (kmax)
    real(KIND=real64)    :: theta
    real(KIND=real64)    :: flat
    real(KIND=real64)    :: rang
    real(KIND=real64)    :: rsin1
    real(KIND=real64)    :: rcos1
    real(KIND=real64)    :: rcos2
    real(KIND=real64)    :: rate
    real(KIND=real64)    :: aa
    real(KIND=real64)    :: bb
    logical :: notfound(kMax)

    kmx=nlm_getoz
    !
    !     find closest place in the data according to input slat.
    !
    if(date.gt.year_getoz) date=date-year_getoz

    imon=date/mon_getoz + 1

    if(imon.lt.1)imon=1

    isea=imon/3 + 1

    if(isea.eq.5) isea=1
    if(isea.gt.5) then
       write(12,"('0 ERROR IN ISEA - TERMINATION IN SUBROUTINE GETOZ')")
       write(10,"('0 ERROR IN ISEA - TERMINATION IN SUBROUTINE GETOZ')")
       stop 9954
    end if
    theta = 90.0_real64-(180.0_real64/pai)*colrad ! colatitude -> latitude
    ! the 180 degrees are divided into 37 bands with 5deg each
    ! except for the first and last, which have 2.5 deg
    ! The centers of the bands are located at:
    !   90, 85, 80, ..., 5, 0, -5, ..., -85, -90 (37 latitudes)
    flat  = 0.2_real64*theta ! indexing the latitudes: goes from -18. to +18.
    ! find the latitude index before and after each latitude
    la    = 19.501e0_real64-flat !
    ll    = 19.001e0_real64-flat

    !
    !     find sin and cos coefficients for time interpolation.
    !
    rang=2.0e0_real64*pai*(date-rlag)/year_getoz
    rsin1=sin(rang)
    rcos1=cos(rang)
    rcos2=cos(2.0e0_real64*rang)
    rate=real(19-ll,real64)-flat
    !
    !     ozone interpolation in latitude and time
    !
    do k=1,kmx
       a1(k) =2.5e-1_real64*(ozone(k,la,1)+ozone(k,la,2)+ &
            ozone(k,la,3)+ozone(k,la,4))
       a2(k) =0.5e0_real64*(ozone(k,la,2)-ozone(k,la,4))
       a3(k) =0.5e0_real64*(ozone(k,la,1)-ozone(k,la,3))
       a4(k) =2.5e-1_real64*(ozone(k,la,1)+ozone(k,la,3)- &
            ozone(k,la,2)-ozone(k,la,4))
       b1(k) =2.5e-1_real64*(ozone(k,ll,1)+ozone(k,ll,2)+ &
            ozone(k,ll,3)+ozone(k,ll,4))
       b2(k) =0.5e0_real64*(ozone(k,ll,2)-ozone(k,ll,4))
       b3(k) =0.5e0_real64*(ozone(k,ll,1)-ozone(k,ll,3))
       b4(k) =2.5e-1_real64*(ozone(k,ll,1)+ozone(k,ll,3)- &
            ozone(k,ll,2)-ozone(k,ll,4))
       do3a(k)=a1(k)+rsin1*a2(k)+rcos1*a3(k)+rcos2*a4(k)
       do3b(k)=b1(k)+rsin1*b2(k)+rcos1*b3(k)+rcos2*b4(k)
       ozo3l(k)=do3a(k)+rate*(do3b(k)-do3a(k))
       ozo3l(k)=1.0e-04_real64*ozo3l(k)
    end do
    if(inter_getoz)then
       do l=1,kmax
          notfound(l) = sigmid(l) > ozsig(1)
          if (notfound(l)) then
             kk(l)=kmx
          else
             kk(l)=2
          end if
       end do
       do l=1,kmax
          if (notfound(l)) then
             do k=2,kmx
                if(sigmid(l).gt.ozsig(k-1).and.sigmid(l).le.ozsig(k))then
                   kk(l)=k
                   exit
                end if
             end do
          end if

       end do
    end if
    if(inter_getoz)then
       do l=1,kmax
          aa=(ozo3l(kk(l))-ozo3l(kk(l)-1))/(ozsig(kk(l))-ozsig(kk(l)-1))
          bb=ozo3l(kk(l)-1)-aa*ozsig(kk(l)-1)
          o3l(kmax+1-l)=bb+aa*sigmid(l)
       end do
    end if
    if(.not.inter_getoz)then
       do l=1,nlm_getoz
          o3l(l)=ozo3l(l)
       end do
    endif
  end subroutine getoz_rrtm
end module ModRrtmDriver

