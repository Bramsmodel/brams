module ModSeaSalt
  !# Module to determine the whitecap seasalt flux and mass over waves (surface)
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: The module is based on paper *Modeling sea-salt aerosol in a coupled climate and sectional
  !# microphysical model: mass, optical depth and number concentration*, Atmos. Chem. Phys., 11, 4587–
  !# 4610, 2011 - T. Fan and O. B. Toon, That determine the whitecap flux and mass of sea salt particles over (on surface) the
  !# sea waves.
  !# \[F_n=3.84 \text{x}10^{-6}  (ak_n  sst+bk_n) U_{10}^{3.41}\]
  !# and
  !# \[Mass_n=F_n*M_aer_n*1.0\text{x}10^{-3}\]
  !# Where: <br/>
  !# \[sst \text{ = sea surface temperature [k]}\]
  !# \[U_{10} \text{= wind at 10 m [m/s]}\]
  !# \[M_aer \text{= Aerosol particle mass in each bin}\]
  !# \[ak_n \text{=parameter factor for equation (a)}\]
  !# \[bk_n \text{=parameter factor for equation (b)}\]
  !# \[n \text{= especific bin of aerosol [1,2 or 3]S}\]
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
  !#
  !# **Date**: 2015Jul
  !# @endnote
  !#
  !# @changes
  !#
  !# +
  !# @endchanges
  !# @bug
  !# No active bugs reported now
  !# @endbug
  !#
  !# @todo
  !#  &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !# @endwarning
  !#
  !#--- ----------------------------------------------------------------------------------------
  use mem_leaf, only: &
       leaf_g           ! INTENT(IN)

  use mem_grid, only: &
       dtlt, &
       npatch,       &  ! INTENT(IN)
       deltaxn,      &  ! INTENT(IN)
       deltayn,      &  ! INTENT(IN)
       time,         &
       grid_g,       &
       dzt

  use io_params, only: &
       iupdsst,       & ! INTENT(IN)
       ssttime1,      & ! INTENT(IN)
       ssttime2         ! INTENT(IN)

  use ModBasicFields, only:  &
       BasicFields

  use node_mod, only: &
       mynum

  use mem_aer1, only: &
       aer1_g, &
       AEROSOL

  use aer1_list, only:   &
       marin       ,      &
       matrix_level,      &
       aer_name, &
       nmodes, &
       nspecies, &
       spc_alloc

  use ModAerClim, only: &
       SScam, &
       aerCam

  use ccatt_start, only: &
       ccatt

  use mem_chem1, only: &
       chemistry

  implicit none

  private

  public :: seaSaltDriver
  
  integer,parameter :: nbins=3
  real,parameter,dimension(nbins,0:4) :: c=reshape((/   &
       -2.881e+06,-6.743e+06, 2.181e+06,   &
       -3.003e+13, 1.183e+14,-4.165e+12,   &
       -2.867e+21,-8.148e+20, 3.132e+18,   &
       5.932e+28, 2.404e+27,-9.841e+23,   &
       -2.576e+35,-2.452e+33, 1.085e+29/), &
       (/3,5/))
  !# dimension of c param for nbins
  real,parameter,dimension(nbins,0:4) :: d=reshape((/   &
       7.609e+08, 2.279e+09,-5.800e+08,   &
       1.829e+16,-3.787e+16, 1.105e+15,   &
       6.791e+23, 2.528e+23,-8.297e+20,   &
       -1.616e+31,-7.310e+29, 2.601e+26,   &
       7.188e+37, 7.368e+35,-2.859e+31/), &
       (/3,5/))
  !# dimension of d param for nbins
  real, parameter,dimension(nbins) :: minDiam=(/0.020e-06,0.145e-06,0.419e-06/)
  !# Min diameter for aerosol particles in each bin
  real, parameter,dimension(nbins) :: maxDiam=(/0.145e-06,0.419e-06,2.800e-06/)
  !# Max diameter for aerosol particles in each bin
  real, parameter,dimension(nbins) :: Diam=(minDiam+maxDiam)/2
  !      ak(k)=c(k,4)*diam(k)**4+ &
  !            c(k,3)*diam(k)**3+ &
  !            c(k,2)*diam(k)**2+ &
  !            c(k,1)*diam(k)   + &
  !            c(k,0)
  !      bk(k)=d(k,4)*diam(k)**4+ &
  !            d(k,3)*diam(k)**3+ &
  !            d(k,2)*diam(k)**2+ &
  !            d(k,1)*diam(k)   + &
  !            d(k,0)
  !# Average diameter of aerosol particles
  real, parameter,dimension(nbins) :: ak=(/-3.496e+06, 2.264e+05, 1.433e+06/)
  !# ak factor applied on equation for each bin
  real, parameter,dimension(nbins) :: bk=(/ 1.148e+09,-3.034e+07,-3.815e+08/)
  !# bk factor applied on equation for each bin
  real, parameter :: DrySeaSaltDendity=2.17
  !# Dry seaSalt density in g/cm3
  real, parameter :: pi=3.1415926
  !# Pi constant
  real, parameter :: mp1=(4/3)*Pi*(diam(1)/2)**3*DrySeaSaltDendity*1.0e-06
  !# mass of seasalt for diam in bin 1
  real, parameter :: mp2=(4/3)*Pi*(diam(2)/2)**3*DrySeaSaltDendity*1.0e-06
  !# mass of seasalt for diam in bin 2
  real, parameter :: mp3=(4/3)*Pi*(diam(3)/2)**3*DrySeaSaltDendity*1.0e-06
  !# mass of seasalt for diam in bin 3

  real, parameter,dimension(nbins) :: massSSaer=(/mp1,mp2,mp3/)
  !MAss of aerosol for each bin
  logical, parameter :: dump=.false.
  !# Dump results in grads. Notice: Must be used just when running with one thread

contains

  subroutine aersaltFlux(flux,w10,sst)
    !# Obtain the seasalt flux over waves
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This subroutine solves the flux of seasalt aerosol in each specific bean
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2015Jul
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    real, intent(in) :: w10
    !# Wind at 10 m [m/s]
    real, intent(in) :: sst
    !# Sea surface temperature [K]
    real, intent(out) :: flux(3)
    !# Flux of aerosol #/m2s
    integer :: k

    do k=1,nbins

       flux(k)=3.84e-06*(ak(k)*sst+bk(k))*w10**3.41
       !print *, massSSaer(k)
    end do

  end subroutine aersaltFlux


  subroutine massSeaSalt(ia,iz,ja,jz,w10,sst,landMask,seaSaltMass)
    !# Obtain the mass of seasalt
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This subroutine determine the ammount of seasalt produced in each grid cell
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2015Jul
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(IN) :: ia
    !# First position (without halo) in x dimension
    integer, intent(IN) :: iz
    !# Last position (without halo) in x dimension
    integer, intent(IN) :: ja
    !# First position (without halo) in y dimension
    integer, intent(IN) :: jz
    !# Last position (without halo) in y dimension
    real,intent(IN) :: w10(ia:iz,ja:jz)
    !# wind at 10m [m/s]
    real,intent(IN) :: sst(ia:iz,ja:jz)
    !# Sea Surface Temperature [K]
    logical, intent(IN) :: landMask(ia:iz,ja:jz)
    !# Mask land surface (.false.=ocean)
    real,intent(out) :: seaSaltMass(nbins,ia:iz,ja:jz)
    !# Mass of SeaSalt for each grid cell [kg/m2s]

    integer :: i,j,k
    real,dimension(nbins) :: flux

    seaSaltMass=0.0

    do i=ia,iz
       do j=ja,jz
          !WRITE(*,FMT='("LFR: ",2(I2.2,1X),L1)'), i,j,landmask(i,j)
          if(landMask(i,j)) cycle
          call aersaltFlux(flux,w10(i,j),sst(i,j))
          do k=1,nbins
             seaSaltMass(k,i,j)=flux(k)*massSSaer(k)*1.0e-3 !Kg/m2S
          end do
       end do
    end do

  end subroutine massSeaSalt

  subroutine SeaSaltDriver(ia,iz,ja,jz,ngrid,m2,m3, oneBasicFields)
    !# Driver to mass of seasalt routines
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This subroutine adapt the BRAMS to use the seasalt routines
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2015Jul
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    integer, intent(IN) :: ia
    !# First position (without halo) in x dimension
    integer, intent(IN) :: iz
    !# Last position (without halo) in x dimension
    integer, intent(IN) :: ja
    !# First position (without halo) in y dimension
    integer, intent(IN) :: jz
    !# Last position (without halo) in y dimension
    integer, intent(IN) :: ngrid
    !# Current grid
    integer, intent(IN) :: m2
    !# Total of points in x dimension
    integer, intent(IN) :: m3
    !# total of points in y dimension
    real :: seaSaltMass(3,ia:iz,ja:jz)
    !# Mass of SeaSalt for each grid cell [kg/m2s]
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    
    real,parameter :: fcu =1.e+9 !=> ppbm

    integer :: np,i,j,nmd,ispc,imode
    logical :: landMask(ia:iz,ja:jz) ! .false. is ocean
    real :: dz
    real :: sst(ia:iz,ja:jz)
    real :: w10(ia:iz,ja:jz)
    real :: timefac_sst
    character(LEN=8) :: ctime
    character(LEN=6) :: clat,clon
    integer :: imonth

    if(AEROSOL == 0) return

    do i=ia,iz
       do j=ja,jz
          if(leaf_g(ngrid)%patch_area(i,j,1)==1) then
             landMask(i,j)=.false.
          else
             landMask(i,j)=.true.
          end if
          !WRITE(*,FMT='("LFR: ",2(I2.2,1X),2(F6.2,1X),L1)'), i,j,grid_g(ngrid)%glat(i,j),grid_g(ngrid)%glon(i,j),landmask(i,j)
          ! Time interpolation factor for updating SST

          if (iupdsst == 0) then
             timefac_sst = 0.
          else
             timefac_sst = (time - ssttime1(ngrid)) / (ssttime2(ngrid) - ssttime1(ngrid))
          endif
          sst(i,j) =   leaf_g(ngrid)%seatp(i,j) + &
               ( leaf_g(ngrid)%seatf(i,j) - leaf_g(ngrid)%seatp(i,j) ) * timefac_sst

          w10(i,j)=(sqrt(oneBasicFields%up(1,i,j)**2+oneBasicFields%vp(1,i,j)**2)+ &
               sqrt(oneBasicFields%up(2,i,j)**2+oneBasicFields%vp(2,i,j)**2))/2
       end do
    end do

    call massSeaSalt(ia,iz,ja,jz,w10,sst,landMask,seaSaltMass)

    if (aerosol==-1 .and. .not. (CCATT==1 .and. chemistry >= 1)) then
       do imonth=1,12
          do i=ia,iz
             do j=ja,jz
                do nmd=1,3
                   dz = grid_g(ngrid)%rtgt(i,j)/dzt(2)
                   aercam(SScam,imonth)%aer(1,i,j)=aercam(SScam,imonth)%aer(1,i,j) &
                        +fcu*seaSaltMass(nmd,i,j)/(dz)
                enddo
             enddo
          enddo
       enddo
       return
    endif

    if(AEROSOL == 1 .and. (CCATT==1 .and. chemistry >= 1)) then
       do i=ia,iz
          do j=ja,jz
             do nmd=1,3
                dz = grid_g(ngrid)%rtgt(i,j)/dzt(2)
                !units -> 10+9 kg/m3/sec
                aer1_g(nmd,marin,ngrid)%sc_src(1,i,j)=fcu*seaSaltMass(nmd,i,j)/(dz)!*oneBasicFields%dn0(2,i,j))
             end do
          end do
       end do
       !If dump is ON the code will create a grads file with output of sea salt aerosols
       !WARNING: Use it only if the model run with one processor (thread).
       if(dump) call writeSeaSalt(m2,m3,ia,iz,ja,jz,              &
            aer1_g(1,marin,ngrid)%sc_src(2,:,:),     &
            aer1_g(2,marin,ngrid)%sc_src(2,:,:),     &
            aer1_g(3,marin,ngrid)%sc_src(2,:,:),     &
            sst,                                     &
            w10,                                     &
            landMask,                                &
            grid_g(ngrid)%glat(1,1),                 &
            grid_g(ngrid)%glon(1,1),                 &
            time                                     &
            )
    elseif(AEROSOL == 2 .and. (CCATT==1 .and. chemistry >= 1)) then
       if(spc_alloc(1,7,5) == 0 .or. spc_alloc(1,8,5)== 0) &
            stop "memory not allocated for sea salt species: ssa_seas and ssc_seas"
       !
       ispc  = 5
       do i=ia,iz
          do j=ja,jz
             dz = grid_g(ngrid)%rtgt(i,j)/dzt(2)
             !ssa_seas Sea Salt accum. mode   !units -> 10+9 kg/m3/sec
             imode = 7
             aer1_g(imode,ispc,ngrid)%sc_src(1,i,j)=fcu*(seaSaltMass(1,i,j)+seaSaltMass(2,i,j))/dz
             !ssc_seas Sea Salt coarse mode   !units -> 10+9 kg/m3/sec
             imode = 8
             aer1_g(imode,ispc,ngrid)%sc_src(1,i,j)=fcu*(seaSaltMass(3,i,j))/dz
          end do
       end do
       !print*,"SEASALT",maxval(aer1_g(imode,ispc,ngrid)%sc_src(1,:,:)),&
       !                 minval(aer1_g(8,5,ngrid)%sc_src(1,:,:));call flush(6)
    endif
  end subroutine SeaSaltDriver

  subroutine writeSeaSalt(nx,ny,ia,iz,ja,jz,var1,var2,var3,sst,w10,mask,lat,lon,time)
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real,intent(in) :: var1(nx,ny)
    real,intent(in) :: var2(nx,ny)
    real,intent(in) :: var3(nx,ny)
    real,intent(in) :: lat
    real,intent(in) :: lon
    real,intent(in) :: time
    real,intent(in) :: sst(ia:iz,ja:jz)
    real,intent(in) :: w10(ia:iz,ja:jz)
    logical, intent(in) :: mask(ia:iz,ja:jz)

    integer,parameter :: nvar =6
    integer,save :: nt = 1
    integer,save :: it=0
    ! loop counter
    integer,save :: irec=0  ! Record number
    integer :: i
    integer :: j
    real, parameter :: timewrite=1200.0
    real :: var4(nx,ny),var5(nx,ny),var6(nx,ny)


    if (.not. (mod(time + .001,timewrite) .lt. dtlt .or. time .lt. 0.001)) return

    var4=0.0;var5=0.0;var6=0.0

    do i=ia,iz
       do j=ja,jz
          var4(i,j)=sst(i,j)
          var5(i,j)=w10(i,j)
          if(.not. mask(i,j)) var6(i,j)=1.0
       end do
    end do


    if(it==0) then
       open(30,file='Outfile.gra',action='write',form='unformatted',access='direct',&
            recl=4*nVar*nx*ny,status='replace')
       call writeSeaSaltCtl(nx,ny,lat,lon)
    end if

    !loop over timesteps
    it  = it+1
    irec = irec + 1
    write(30,rec=irec) Var1(:,:),Var2(:,:),Var3(:,:),Var4(:,:),Var5(:,:),Var6(:,:)

  end subroutine writeSeaSalt

  subroutine writeSeaSaltCtl(nx,ny,lat,lon)
    integer, intent(IN) :: nx,ny
    real, intent(IN) :: lat,lon
    open(31,FILE='Outfile.ctl',STATUS='replace')
    write (31,fmt='(A)') 'dset ./Outfile.gra'
    write (31,fmt='(A)') 'undef  -0.9990000E+34'
    write (31,fmt='(A)') 'title Sea Salt Aerosol'
    write (31,fmt='(A)') 'options little_endian'
    write (31,fmt='(A,I2.2,A,F8.4,1X,F8.4)') 'xdef ',nx,' linear ',lon,0.95
    write (31,fmt='(A,I2.2,A,F8.4,1X,F8.4)') 'ydef ',ny,' linear ',lat,0.89
    write (31,fmt='(A)') 'zdef 1 linear 1.0 1.0'
    write (31,fmt='(A)') 'tdef 5 linear 00:00Z01jan2012 1hr'
    write (31,fmt='(A)') 'vars 6'
    write (31,fmt='(A)') 'seas1       0 99    - Sea Salt 0.020e-06 a 0.145e-06  [ ]'
    write (31,fmt='(A)') 'seas2       0 99    - Sea Salt 0.145e-06 a 0.419e-06  [ ]'
    write (31,fmt='(A)') 'seas3       0 99    - Sea Salt 0.419e-06 a 2.800e-06  [ ]'
    write (31,fmt='(A)') 'sst         0 99    - Sea Surface Temperature         [K]'
    write (31,fmt='(A)') 'w10         0 99    - Wind 10m                        [m/s]'
    write (31,fmt='(A)') 'mask        0 99    - Mask land Ocean=1.0             [ ]'
    write (31,fmt='(A)') 'endvars'
  end subroutine writeSeaSaltCtl
end module ModSeaSalt

