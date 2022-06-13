module ModRbnd

  use ModTurbFields, only: &
       TurbFields
  
  use ModTurbKE, only: &
       tkeinit
  
  use ModMicrophysicsMisc, only: &
       negadj1

  use mem_grid, only: dtlt, &
       dtlv, &
       ibnd, &
       jbnd, &
       jdim, &
       lsflg,&
       dzm,  &
       dzt,  &
       grid_g, &
       ngrid,  &
       nxtnest,&
       cphas,  &
       if_adap,&
       nstbot, &
       nsttop, &
       nfpt, &
       distim, &
       nnzp,   &
       dts,    &
       icorflg,&
       nz,     &
       nzp,    &
       ztop,   &
       nnz,    &
       zmn,    &
       zt,     &
       naddsc

  use mem_scratch, only: vctr17, &
       vctr18, &
       scratch,&
       vctr2,  &
       vctr5 

  use mem_tend, only: tend       !tend%ut

  use ModBasicFields, only: &
       BasicFields
  
  use node_mod, only: ia,      &
       iz,      &
       ja,      &
       jz,      &
       mxp,     &
       myp,     &
       mzp,     &
       ibcon,   &
       mynum,   &
       nodemxp, &
       nodemyp

  use micphys, only: level
  use ref_sounding, only: u01dn, &
       v01dn

  use ccatt_start, only: ccatt
  use mem_chem1, only: chemistry,nspecies_transported

  use ModScalarTable, only: &
       ScalarTable

  implicit none

  include "constants.h"

  private
  public :: LatSetScalar
  public :: TopSetScalar
  public :: TopSet2Scalar
  public :: BotSetScalar
  public :: BotSetAdapScalar
  public :: KeepTracersNonnegScalar
  public :: latbnd 
  public :: latnormv
  public :: vpsets  
  public :: topset
  public :: botset
  public :: dumset
  public :: rayft 
  public :: rayf 
  public :: trsets
  public :: botset_adap
  public :: rayf_adap

contains

  subroutine LatSetScalar(m1,m2,m3,ia,iz,ja,jz,ibcon,vnam,ap,uc,vc,dxu,dxm,dyv,dym)

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon
    character(len=*) :: vnam
    real, pointer, intent(in) :: ap(:,:,:)
    ! pointer intent(in), values intent(inout)
    real, pointer, intent(in) :: uc(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: vc(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dxu(:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dxm(:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dyv(:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dym(:,:)
    ! pointer and values intent(in)

    integer :: i,j,k,lbw,lbe,lbs,lbn
    real :: thresh,dtlx,c1,dxr,dyr

    character(len=*), parameter :: h="**(LatSetScalar)**"

    if (.not. associated(ap)) then
       call fatal_error(h//" field ap not associated")
    end if

    if (iand(ibcon,1) .gt. 0) lbw = ia - 1
    if (iand(ibcon,2) .gt. 0) lbe = iz + 1
    if (iand(ibcon,4) .gt. 0) lbs = ja - 1
    if (iand(ibcon,8) .gt. 0) lbn = jz + 1

    thresh = 0.
    if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'W' .or. vnam .eq.'P') then
       dtlx = dtlv
    else
       dtlx = dtlt
    endif

    if (ibnd .ne. 4 .and. vnam .ne. 'U' .and. lsflg .ne. 3) then

       !     Western and Eastern boundaries for zero gradient option

       if (lsflg .eq. 0) then
          if (iand(ibcon,1) .gt. 0) then
             do j = 1,m3
                do k = 1,m1
                   ap(k,lbw,j) = ap(k,ia,j)
                enddo
             enddo
          endif
          if (iand(ibcon,2) .gt. 0) then
             do j = 1,m3
                do k = 1,m1
                   ap(k,lbe,j) = ap(k,iz,j)
                enddo
             enddo
          endif
       else

          !     Western boundary for lsflg = 1 or 2

          if (iand(ibcon,1) .gt. 0) then
             do j = 1,m3-1 !m3  !ALF
                if (vnam .eq. 'V') then
                   dxr = dxm(ia,j) / dxm(lbw,j)
                   c1 = .5 * dtlx * dxm(lbw,j)
                   do k = 1,m1
                      vctr17(k) = -c1 * (uc(k,lbw,j) + uc(k,lbw,j+jdim))
                   enddo
                elseif (vnam .eq. 'W') then
                   dxr = dxu(ia,j) / dxu(lbw,j)
                   c1 = .5 * dtlx * dxu(lbw,j)
                   do k = 1,m1-1 !m1  !ALF
                      vctr17(k) = -c1 * (uc(k,lbw,j) + uc(k+1,lbw,j))
                   enddo
                else
                   dxr = dxu(ia,j) / dxu(lbw,j)
                   c1 = dtlx * dxu(lbw,j)
                   do k = 1,m1
                      vctr17(k) = -c1 * uc(k,lbw,j)
                   enddo
                endif
                do k = 1,m1
                   vctr18(k) = ap(k,ia,j) + dxr * (ap(k,ia,j) - ap(k,ia+1,j))
                enddo
                do k = 1,m1
                   if (vctr17(k) .ge. thresh) then
                      ap(k,lbw,j) = vctr18(k)
                   elseif (lsflg .eq. 1) then
                      ap(k,lbw,j) = ap(k,ia,j)
                   endif
                enddo
             enddo
          endif

          !     Eastern Boundary for LSFLG = 1 or 2

          if (iand(ibcon,2) .gt. 0) then
             do j = 1,m3-1 !m3  !ALF
                if (vnam .eq. 'V') then
                   dxr = dxm(iz-1,j) / dxm(iz,j)
                   c1 = .5 * dtlx * dxm(iz,j)
                   do k = 1,m1
                      vctr17(k) = c1 * (uc(k,iz,j) + uc(k,iz,j+jdim))
                   enddo
                elseif (vnam .eq. 'W') then
                   dxr = dxu(iz-1,j) / dxu(iz,j)
                   c1 = .5 * dtlx * dxu(iz,j)
                   do k = 1,m1-1 !m1  !ALF
                      vctr17(k) = c1 * (uc(k,iz,j) + uc(k+1,iz,j))
                   enddo
                else
                   dxr = dxu(iz-1,j) / dxu(iz,j)
                   c1 = dtlx * dxu(iz,j)
                   do k = 1,m1
                      vctr17(k) = c1 * uc(k,iz,j)
                   enddo
                endif
                do k = 1,m1
                   vctr18(k) = ap(k,iz,j) + dxr * (ap(k,iz,j) - ap(k,iz-1,j))
                enddo
                do k = 1,m1
                   if (vctr17(k) .ge. thresh) then
                      ap(k,lbe,j) = vctr18(k)
                   elseif (lsflg .eq. 1) then
                      ap(k,lbe,j) = ap(k,iz,j)
                   endif
                enddo
             enddo
          endif
       endif
    endif

    if(jdim.eq.1.and.jbnd.ne.4.and.vnam.ne.'V'.and.lsflg.ne.3)then

       !     Southern and Northern boundaries for zero gradient option

       if (lsflg .eq. 0) then
          if (iand(ibcon,4) .gt. 0) then
             do i = 1,m2
                do k = 1,m1
                   ap(k,i,lbs) = ap(k,i,ja)
                enddo
             enddo
          endif
          if (iand(ibcon,8) .gt. 0) then
             do i = 1,m2
                do k = 1,m1
                   ap(k,i,lbn) = ap(k,i,jz)
                enddo
             enddo
          endif
       else

          !     Southern boundary for LSFLG = 1 or 2

          if (iand(ibcon,4) .gt. 0) then
             do i = 1,m2-1 !m2 !ALF
                if (vnam .eq. 'U') then
                   dyr = dym(i,ja) / dym(i,lbs)
                   c1 = .5 * dtlx * dym(i,lbs)
                   do k = 1,m1
                      vctr17(k) = -c1 * (vc(k,i,lbs) + vc(k,i+1,lbs))
                   enddo
                elseif (vnam .eq. 'W') then
                   dyr = dyv(i,ja) / dyv(i,lbs)
                   c1 = .5 * dtlx * dyv(i,lbs)
                   do k = 1,m1-1 !m1  !ALF
                      vctr17(k) = -c1 * (vc(k,i,lbs) + vc(k+1,i,lbs))
                   enddo
                else
                   dyr = dyv(i,ja) / dyv(i,lbs)
                   c1 = dtlx * dyv(i,lbs)

                   !--(DMK-CCATT-INI)-----------------------------------------------------
                   !srf - fix from rams 60
                   do k = 1,m1
                      !--(DMK-CCATT-OLD)-----------------------------------------------------
                      !              do k = 1,nz
                      !--(DMK-CCATT-FIM)-----------------------------------------------------         

                      vctr17(k) = -c1 * vc(k,i,lbs)
                   enddo
                endif
                do k = 1,m1
                   vctr18(k) = ap(k,i,ja) + dyr * (ap(k,i,ja) - ap(k,i,ja+1))
                enddo
                do k = 1,m1
                   if (vctr17(k) .ge. thresh) then
                      ap(k,i,lbs) = vctr18(k)
                   elseif (lsflg .eq. 1) then
                      ap(k,i,lbs) = ap(k,i,ja)
                   endif
                enddo
             enddo
          endif

          !     Northern Boundary for LSFLG = 1 or 2

          if (iand(ibcon,8) .gt. 0) then
             do i = 1,m2-1 !m2 !ALF
                if (vnam .eq. 'U') then
                   dyr = dym(i,jz-1) / dym(i,jz)
                   c1 = .5 * dtlx * dym(i,jz)
                   do k = 1,m1
                      vctr17(k) = c1 * (vc(k,i,jz) + vc(k,i+1,jz))
                   enddo
                elseif (vnam .eq. 'W') then
                   dyr = dyv(i,jz-1) / dyv(i,jz)
                   c1 = .5 * dtlx * dyv(i,jz)
                   do k = 1,m1-1 !m1  !ALF
                      vctr17(k) = c1 * (vc(k,i,jz) + vc(k+1,i,jz))
                   enddo
                else
                   dyr = dyv(i,jz-1) / dyv(i,jz)
                   c1 = dtlx * dyv(i,jz)
                   do k = 1,m1
                      vctr17(k) = c1 * vc(k,i,jz)
                   enddo
                endif
                do k = 1,m1
                   vctr18(k) = ap(k,i,jz) + dyr * (ap(k,i,jz) - ap(k,i,jz-1))
                enddo
                do k = 1,m1
                   if (vctr17(k) .ge. thresh) then
                      ap(k,i,lbn) = vctr18(k)
                   elseif (lsflg .eq. 1) then
                      ap(k,i,lbn) = ap(k,i,jz)
                   endif
                enddo
             enddo
          endif
       endif
    endif

    return
  end subroutine LatSetScalar


  subroutine TopSetScalar(m1,m2,m3,ia,iz,ja,jz,ibcon,ap,vnam)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon
    character(len=*) :: vnam
    real, pointer, intent(in) :: ap(:,:,:)
    ! pointer intent(in), values intent(inout)

    integer :: i,j
    real :: dzmr,dztr

    dzmr = dzm(m1-2) / dzm(m1-1)
    dztr = dzt(m1-2) / dzt(m1-1)

    !     Computation of all prognostic variables (other than W) at
    !       level NZP by extrapolation from below

    if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'P') then
       do j = 1,m3
          do i = 1,m2
             ap(m1,i,j) = ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j))
          enddo
       enddo
    endif
    if (vnam .eq. 'T') then
       do j = 1,m3
          do i = 1,m2
             ap(m1,i,j) = max(0.,ap(m1-1,i,j)+dzmr*(ap(m1-1,i,j)-ap(m1-2,i,j)))
          enddo
       enddo
    endif

    return
  end subroutine TopSetScalar

  subroutine TopSet2Scalar(m1,m2,m3,ia,iz,ja,jz,ibcon,ap,vnam)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon
    character(len=*) :: vnam
    real, pointer, intent(in) :: ap(:,:,:)
    ! pointer intent(in), values intent(inout)

    integer :: i, j
    !real :: dzmr,dztr

    !dzmr = dzm(m1-2) / dzm(m1-1)
    !dztr = dzt(m1-2) / dzt(m1-1)

    !     Computation of all prognostic variables (other than W) at
    !       level NZP by extrapolation from below

    !if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'P') then
    !  do j = 1,m3
    !      do i = 1,m2
    !         ap(m1,i,j) = ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j))
    !      enddo
    !   enddo
    !endif
    !if (vnam .eq. 'T') then
    do j = 1,m3
       do i = 1,m2
          ap(m1,i,j) = ap(m1-1,i,j)
       enddo
    enddo
    !endif

    return
  end subroutine TopSet2Scalar

  subroutine BotSetScalar(m1,m2,m3,ia,iz,ja,jz,ibcon,aa,vnam)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon
    character(len=*) :: vnam
    real, pointer, intent(in) :: aa(:,:,:)
    ! pointer intent(in), values intent(inout)

    integer :: i,j
    real :: dzmr

    if (vnam .eq. 'P') then
       dzmr = dzm(2) / dzm(1)
       do i = 1,m2
          do j = 1,m3
             aa(1,i,j) = aa(2,i,j) + (aa(2,i,j) - aa(3,i,j)) * dzmr
          enddo
       enddo
    else
       do i = 1,m2
          do j = 1,m3
             aa(1,i,j) = aa(2,i,j)
          enddo
       enddo
    endif

    return
  end subroutine BotSetScalar

  subroutine BotSetAdapScalar(m1,m2,m3,ia,iz,ja,jz,ibcon,lpx,aa,vnam)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon
    character(len=*) :: vnam
    integer, intent(in) :: lpx(m2,m3)
    real, pointer, intent(in) :: aa(:,:,:)
    ! pointer intent(in), values intent(inout)


    integer :: i,j,k,ka
    real :: dzmr

    if (vnam .eq. 'P') then
       dzmr = dzm(2) / dzm(1)
       do i = 1,m2
          do j = 1,m3
             ka = lpx(i,j)
             do k = ka-1,1,-1
                aa(k,i,j) = aa(k+1,i,j) + (aa(k+1,i,j) - aa(k+2,i,j)) * dzmr
             enddo
          enddo
       enddo
    else
       do i = 1,m2
          do j = 1,m3
             ka = lpx(i,j)
             do k = ka-1,1,-1
                aa(k,i,j) = aa(k+1,i,j)
             enddo
          enddo
       enddo
    endif

    if (vnam == 'U' .or. vnam == 'V' .or. vnam == 'W') then
       do i = 1,m2
          do j = 1,m3
             ka = lpx(i,j)
             do k = ka-1,1,-1
                aa(k,i,j) = 0.
             enddo
          enddo
       enddo
    endif

    return
  end subroutine BotSetAdapScalar

  subroutine KeepTracersNonnegScalar(mxyzp,scp)
    integer, intent(in) :: mxyzp
    real, pointer, intent(in) :: scp(:)
    ! pointer intent(in), values intent(inout)

    integer :: i 
    do i = 1,mxyzp
       scp(i) = max(0.,scp(i))
    enddo
  end subroutine KeepTracersNonnegScalar

  subroutine latbnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon,nxtnest,ngrid,ibnd,jbnd, &
       lpu,lpv,up,uc,ut,vp,vc,vt,dxt,dyt)
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon
    integer, intent(in) :: nxtnest(:)
    integer, intent(in) :: ngrid
    integer, intent(in) :: ibnd 
    integer, intent(in) :: jbnd 
    real,    pointer, intent(in)    :: lpu(:,:)
    real,    pointer, intent(in)    :: lpv(:,:)
    real,    pointer, intent(in)    :: up(:,:,:)
    real,    pointer, intent(in)    :: uc(:,:,:)
    real, dimension(mzp,mxp,myp), intent(in) :: ut
    real,    pointer, intent(in)    :: vp(:,:,:)
    real,    pointer, intent(in)    :: vc(:,:,:)
    real, dimension(mzp,mxp,myp), intent(in) :: vt
    real,    pointer, intent(in)    :: dxt(:,:)
    real,    pointer, intent(in)    :: dyt(:,:)

    !     This routine drives the computation of the radiative lateral
    !     boundary condition for normal velocity on the coarsest grid
    !     and the recomputation of all boundary tendencies on nested grids
    !     after the first nested grid timestep.

    if (nxtnest(ngrid) .eq. 0) then

       !         Radiative and/or mesoscale compensation region lateral
       !            boundary conditions.

       if (ibnd .le. 3 .or. jbnd .le. 3) then

          call latnormv(mzp,mxp,myp,ia,iz,ja,jz,ibcon                  &
               ,lpu ,lpv  &
               ,up  ,uc   &
               ,ut  ,vp   &
               ,vc  ,vt   &
               ,dxt ,dyt  )

       endif

    endif
    return
  end subroutine latbnd

  subroutine latnormv(m1,m2,m3,ia,iz,ja,jz,ibcon,lpu_R,lpv_R  &
       ,up,uc,ut,vp,vc,vt,dxt,dyt)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon

    real, pointer, intent(in) :: lpu_R(:,:) 
    real, pointer, intent(in) :: lpv_R(:,:)
    real, pointer, intent(in) :: dxt(:,:)
    real, pointer, intent(in) :: dyt(:,:)

    real, pointer, intent(in) :: up(:,:,:) 
    real, pointer, intent(in) :: uc(:,:,:) 
    real, dimension(m1,m2,m3) :: ut  !TO recebe pointer 1D
    real, pointer, intent(in) :: vp(:,:,:) 
    real, pointer, intent(in) :: vc(:,:,:) 
    real, dimension(m1,m2,m3) :: vt  !TO recebe pointer 1D

    integer :: i,j,k
    real :: dxl,dxr,cphx,cphy
    integer, dimension(m2,m3) :: lpu,lpv

    !     This routine ultimately updates tendencies at lateral boundaries
    !     after first diagnosing appropriate phase speeds.
    !
    !     IBND and JBND are flags for the radiative type in the X and Y
    !     direction. Their meaning is:
    !
    !        IBND=1......Klemp-Wilhelmson (1978) type; phase speed given
    !                    by CPHAS
    !        IBND=2......Klemp-Lilly (1980) type; doppler shifted phase
    !                    speed constant with height and diagnosed from
    !                    average of Orlanski speeds, i.e. function of
    !                    only (X,Y)
    !        IBND=3......Orlanski(1974) type; Phase speeds diagnosed
    !                    from local conditions and function of (x,y,z)

    !     Calculate the diagnostic phase
    !     speed. The Orlanski(1976) leapfrog method is to use three time
    !     levels of information, namely the T-2, T-1, and T level to
    !     evaluate the phase speed given by - du/dt / du/dx = u + C.
    !     If this is to be an Orlanski or Klemp-Lilly type boundary in the x
    !     direction then this following diagnostic procedure is necessary.

    !     If this is the first call to a routine, initialize the phase
    !       speed arrays if necessary.

    lpu=int(lpu_R);lpv=int(lpv_R)

    if (ibcon.eq.0) return

    !     first compute "X" boundaries.

    if (iand(ibcon,1) .ne. 0) then
       do j = 1,m3
          dxl = 1. / (dtlv * dxt(2,j))
          do k = lpu(1,j),m1

             cphx = min(0.,max(-dxl,(up(k,1,j)-cphas)))
             ut(k,1,j) = ut(k,1,j) - cphx * dxt(2,j)  &
                  * (up(k,2,j) + ut(k,2,j) * dtlv - up(k,1,j))

          enddo
       enddo
    endif

    if (iand(ibcon,2) .ne. 0) then
       do j = 1,m3
          dxr = 1. / (dtlv * dxt(m2-1,j))
          do k = lpu(m2-1,j),m1

             cphx = max(0.,min(dxr,(up(k,m2-1,j)+cphas)))
             ut(k,m2-1,j) = ut(k,m2-1,j) - cphx * dxt(m2-1,j)  &
                  * (up(k,m2-1,j) - (up(k,m2-2,j) + ut(k,m2-2,j) * dtlv))

          enddo
       enddo
    endif

    !     South and north boundaries.

    if (jdim .eq. 1) then

       if (iand(ibcon,4) .ne. 0) then
          do i = 1,m2
             dxl = 1. / (dtlv * dyt(i,2))
             do k = lpv(i,1),m1
                cphy = min(0.,max(-dxl,(vp(k,i,1)-cphas)))
                vt(k,i,1) = vt(k,i,1) - cphy * dyt(i,2)  &
                     * (vp(k,i,2) + vt(k,i,2) * dtlv - vp(k,i,1))
             enddo
          enddo
       endif

       if (iand(ibcon,8) .ne. 0) then
          do i = 1,m2
             dxr = 1. / (dtlv * dyt(i,m3-1))
             do k = lpv(i,m3-1),m1
                cphy = max(0.,min(dxr,(vp(k,i,m3-1)+cphas)))
                vt(k,i,m3-1) = vt(k,i,m3-1) - cphy * dyt(i,m3-1)  &
                     * (vp(k,i,m3-1) - (vp(k,i,m3-2) + vt(k,i,m3-2) * dtlv))
             enddo
          enddo
       endif

    endif
    return
  end subroutine latnormv

  subroutine vpsets(mzp,mxp,myp,ia,iz,ja,jz,ibcon,nstbot, &
       up,vp,wp,pp,uc,vc,wc,pc,dxu,dxm,dyv,dym,&
       lpu,lpv,lpw, oneBasicFields)
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon
    integer, intent(in) :: nstbot
    ! pointer intent(in), values intent(inout)
    real, pointer, intent(in) :: up(:,:,:) 
    real, pointer, intent(in) :: vp(:,:,:) 
    real, pointer, intent(in) :: wp(:,:,:) 
    real, pointer, intent(in) :: pp(:,:,:) 
    real, pointer, intent(in) :: uc(:,:,:) 
    real, pointer, intent(in) :: vc(:,:,:) 
    real, pointer, intent(in) :: wc(:,:,:) 
    real, pointer, intent(in) :: pc(:,:,:) 
    real, pointer, intent(in) :: dxu(:,:)
    real, pointer, intent(in) :: dxm(:,:)
    real, pointer, intent(in) :: dyv(:,:)
    real, pointer, intent(in) :: dym(:,:)
    real, pointer, intent(in) :: lpu(:,:)
    real, pointer, intent(in) :: lpv(:,:)
    real, pointer, intent(in) :: lpw(:,:)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    
    if (nxtnest(ngrid) .eq. 0) then
       call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'U'              &
            ,up ,up &
            ,vp ,dxu &
            ,dxm,dyv &
            ,dym )
       call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'V'              &
            ,vp ,up  &
            ,vp ,dxu  &
            ,dxm,dyv  &
            ,dym )
       call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'W'              &
            ,wp ,up  &
            ,vp ,dxu  &
            ,dxm,dyv  &
            ,dym )
       call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'P'              &
            ,pp ,up  &
            ,vp ,dxu  &
            ,dxm,dyv  &
            ,dym )
    endif
    if (nsttop .eq. 1) then
       call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
            ,up ,'U')
       call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
            ,vp ,'V')
    endif

    if (nstbot .eq. 1) then
       if (if_adap == 0) then
          call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
               ,up,'U')
          call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
               ,vp,'V')
       else
          call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,int(lpu)  &
               ,up,'U')
          call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,int(lpv)  &
               ,vp,'V')
       endif
    endif

    call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
         ,pp,'P')

    if (if_adap == 0) then
       call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
            ,pp,'P')
    else
       call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,int(lpw)  &
            ,pp,'P')
    endif

    call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
         ,oneBasicFields%wp,'W')

    call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
         ,oneBasicFields%up,'U')
    call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
         ,oneBasicFields%vp,'V')

    if (nxtnest(ngrid) .eq. 0) then
       call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'U'  &
            ,uc ,up  &
            ,vp ,dxu  &
            ,dxm,dyv  &
            ,dym )
       call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'V'  &
            ,vc ,up &
            ,vp ,dxu &
            ,dxm,dyv &
            ,dym )
       call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon ,'W' &
            ,wc ,up &
            ,vp ,dxu &
            ,dxm,dyv &
            ,dym )
       call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'P'  &
            ,pc ,up  &
            ,vp ,dxu  &
            ,dxm,dyv  &
            ,dym )
    endif
    if (nsttop .eq. 1) then
       call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
            ,uc ,'U')
       call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
            ,vc ,'V')
    endif

    if (nstbot .eq. 1) then
       if (if_adap == 0) then
          call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
               ,uc,'U')
          call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
               ,vc,'V')
       else
          call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,int(lpu)  &
               ,uc,'U')
          call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,int(lpv)  &
               ,vc,'V')
       endif
    endif


    call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
         ,pc,'P')

    if (if_adap == 0) then
       call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
            ,pc,'P')
    else
       call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,int(lpw)  &
            ,pc,'P')
    endif

    call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
         ,oneBasicFields%wc,'W')
    call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
         ,oneBasicFields%uc,'U')
    call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
         ,oneBasicFields%vc,'V')

    return
  end subroutine vpsets

  subroutine latset(m1,m2,m3,ia,iz,ja,jz,ibcon,vnam,ap,uc,vc,dxu,dxm,dyv,dym)
    integer,intent(in) :: m1,m2,m3,ia,iz,ja,jz,ibcon
    real,  pointer, intent(in) :: ap(:,:,:)
    real,  pointer, intent(in) :: uc(:,:,:)
    real,  pointer, intent(in) :: vc(:,:,:)
    real,  pointer, intent(in) :: dxu(:,:)
    real,  pointer, intent(in) :: dxm(:,:)
    real,  pointer, intent(in) :: dyv(:,:)
    real,  pointer, intent(in) :: dym(:,:)
    character(len=*) :: vnam

    integer :: i,j,k,lbw,lbe,lbs,lbn
    real :: thresh,dtlx,c1,dxr,dyr

    if (iand(ibcon,1) .gt. 0) lbw = ia - 1
    if (iand(ibcon,2) .gt. 0) lbe = iz + 1
    if (iand(ibcon,4) .gt. 0) lbs = ja - 1
    if (iand(ibcon,8) .gt. 0) lbn = jz + 1

!!$print *, "DEBUG-ALF:latset:m1,m2,m3,lbw,lbe,lbs,lbn=", &
!!$     m1,m2,m3,lbw,lbe,lbs,lbn
!!$call flush(8)

    thresh = 0.
    if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'W' .or. vnam .eq.'P') then
       dtlx = dtlv
    else
       dtlx = dtlt
    endif

    if (ibnd .ne. 4 .and. vnam .ne. 'U' .and. lsflg .ne. 3) then

       !     Western and Eastern boundaries for zero gradient option

       if (lsflg .eq. 0) then
          if (iand(ibcon,1) .gt. 0) then
             do j = 1,m3
                do k = 1,m1
                   ap(k,lbw,j) = ap(k,ia,j)
                enddo
             enddo
          endif
          if (iand(ibcon,2) .gt. 0) then
             do j = 1,m3
                do k = 1,m1
                   ap(k,lbe,j) = ap(k,iz,j)
                enddo
             enddo
          endif
       else

          !     Western boundary for lsflg = 1 or 2

          if (iand(ibcon,1) .gt. 0) then
             do j = 1,m3-1 !m3  !ALF
                if (vnam .eq. 'V') then
                   dxr = dxm(ia,j) / dxm(lbw,j)
                   c1 = .5 * dtlx * dxm(lbw,j)
                   do k = 1,m1
                      vctr17(k) = -c1 * (uc(k,lbw,j) + uc(k,lbw,j+jdim))
                   enddo
                elseif (vnam .eq. 'W') then
                   dxr = dxu(ia,j) / dxu(lbw,j)
                   c1 = .5 * dtlx * dxu(lbw,j)
                   do k = 1,m1-1 !m1  !ALF
                      vctr17(k) = -c1 * (uc(k,lbw,j) + uc(k+1,lbw,j))
                   enddo
                else
                   dxr = dxu(ia,j) / dxu(lbw,j)
                   c1 = dtlx * dxu(lbw,j)
                   do k = 1,m1
                      vctr17(k) = -c1 * uc(k,lbw,j)
                   enddo
                endif
                do k = 1,m1
                   vctr18(k) = ap(k,ia,j) + dxr * (ap(k,ia,j) - ap(k,ia+1,j))
                enddo
                do k = 1,m1
                   if (vctr17(k) .ge. thresh) then
                      ap(k,lbw,j) = vctr18(k)
                   elseif (lsflg .eq. 1) then
                      ap(k,lbw,j) = ap(k,ia,j)
                   endif
                enddo
             enddo
          endif

          !     Eastern Boundary for LSFLG = 1 or 2

          if (iand(ibcon,2) .gt. 0) then
             do j = 1,m3-1 !m3  !ALF
                if (vnam .eq. 'V') then
                   dxr = dxm(iz-1,j) / dxm(iz,j)
                   c1 = .5 * dtlx * dxm(iz,j)
                   do k = 1,m1
                      vctr17(k) = c1 * (uc(k,iz,j) + uc(k,iz,j+jdim))
                   enddo
                elseif (vnam .eq. 'W') then
                   dxr = dxu(iz-1,j) / dxu(iz,j)
                   c1 = .5 * dtlx * dxu(iz,j)
                   do k = 1,m1-1 !m1  !ALF
                      vctr17(k) = c1 * (uc(k,iz,j) + uc(k+1,iz,j))
                   enddo
                else
                   dxr = dxu(iz-1,j) / dxu(iz,j)
                   c1 = dtlx * dxu(iz,j)
                   do k = 1,m1
                      vctr17(k) = c1 * uc(k,iz,j)
                   enddo
                endif
                do k = 1,m1
                   vctr18(k) = ap(k,iz,j) + dxr * (ap(k,iz,j) - ap(k,iz-1,j))
                enddo
                do k = 1,m1
                   if (vctr17(k) .ge. thresh) then
                      ap(k,lbe,j) = vctr18(k)
                   elseif (lsflg .eq. 1) then
                      ap(k,lbe,j) = ap(k,iz,j)
                   endif
                enddo
             enddo
          endif
       endif
    endif

    if(jdim.eq.1.and.jbnd.ne.4.and.vnam.ne.'V'.and.lsflg.ne.3)then

       !     Southern and Northern boundaries for zero gradient option

       if (lsflg .eq. 0) then
          if (iand(ibcon,4) .gt. 0) then
             do i = 1,m2
                do k = 1,m1
                   ap(k,i,lbs) = ap(k,i,ja)
                enddo
             enddo
          endif
          if (iand(ibcon,8) .gt. 0) then
             do i = 1,m2
                do k = 1,m1
                   ap(k,i,lbn) = ap(k,i,jz)
                enddo
             enddo
          endif
       else

          !     Southern boundary for LSFLG = 1 or 2


          if (iand(ibcon,4) .gt. 0) then
             do i = 1,m2-1 !m2 !ALF
                if (vnam .eq. 'U') then
                   dyr = dym(i,ja) / dym(i,lbs)
                   c1 = .5 * dtlx * dym(i,lbs)
                   do k = 1,m1
                      vctr17(k) = -c1 * (vc(k,i,lbs) + vc(k,i+1,lbs))
                   enddo
                elseif (vnam .eq. 'W') then
                   dyr = dyv(i,ja) / dyv(i,lbs)
                   c1 = .5 * dtlx * dyv(i,lbs)
                   do k = 1,m1-1 !m1  !ALF
                      vctr17(k) = -c1 * (vc(k,i,lbs) + vc(k+1,i,lbs))
                   enddo
                else
                   dyr = dyv(i,ja) / dyv(i,lbs)
                   c1 = dtlx * dyv(i,lbs)

                   !--(DMK-CCATT-INI)-----------------------------------------------------
                   !srf - fix from rams 60
                   do k = 1,m1
                      !--(DMK-CCATT-OLD)-----------------------------------------------------
                      !                  do k = 1,nz
                      !--(DMK-CCATT-FIM)-----------------------------------------------------       

                      vctr17(k) = -c1 * vc(k,i,lbs)
                   enddo
                endif
                do k = 1,m1
                   vctr18(k) = ap(k,i,ja) + dyr * (ap(k,i,ja) - ap(k,i,ja+1))
                enddo
                do k = 1,m1
                   if (vctr17(k) .ge. thresh) then
                      ap(k,i,lbs) = vctr18(k)
                   elseif (lsflg .eq. 1) then
                      ap(k,i,lbs) = ap(k,i,ja)
                   endif
                enddo
             enddo
          endif

          !     Northern Boundary for LSFLG = 1 or 2

          if (iand(ibcon,8) .gt. 0) then
             do i = 1,m2-1 !m2 !ALF
                if (vnam .eq. 'U') then
                   dyr = dym(i,jz-1) / dym(i,jz)
                   c1 = .5 * dtlx * dym(i,jz)
                   do k = 1,m1
                      vctr17(k) = c1 * (vc(k,i,jz) + vc(k,i+1,jz))
                   enddo
                elseif (vnam .eq. 'W') then
                   dyr = dyv(i,jz-1) / dyv(i,jz)
                   c1 = .5 * dtlx * dyv(i,jz)
                   do k = 1,m1-1 !m1  !ALF
                      vctr17(k) = c1 * (vc(k,i,jz) + vc(k+1,i,jz))
                   enddo
                else
                   dyr = dyv(i,jz-1) / dyv(i,jz)
                   c1 = dtlx * dyv(i,jz)
                   do k = 1,m1
                      vctr17(k) = c1 * vc(k,i,jz)
                   enddo
                endif
                do k = 1,m1
                   vctr18(k) = ap(k,i,jz) + dyr * (ap(k,i,jz) - ap(k,i,jz-1))
                enddo
                do k = 1,m1
                   if (vctr17(k) .ge. thresh) then
                      ap(k,i,lbn) = vctr18(k)
                   elseif (lsflg .eq. 1) then
                      ap(k,i,lbn) = ap(k,i,jz)
                   endif
                enddo
             enddo
          endif
       endif
    endif

    return
  end subroutine latset


  subroutine topset(m1,m2,m3,ia,iz,ja,jz,ibcon,ap,vnam)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon

    real,    pointer, intent(in) :: ap(:,:,:)
    ! pointer intent(in), values intent(inout)
    character(len=*) :: vnam

    integer :: i,j
    real :: dzmr,dztr

    dzmr = dzm(m1-2) / dzm(m1-1)
    dztr = dzt(m1-2) / dzt(m1-1)

    !     Computation of all prognostic variables (other than W) at
    !       level NZP by extrapolation from below

    if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'P') then
       do j = 1,m3
          do i = 1,m2
             ap(m1,i,j) = ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j))
          enddo
       enddo
    endif
    if (vnam .eq. 'T') then
       do j = 1,m3
          do i = 1,m2
             ap(m1,i,j) = max(0.,ap(m1-1,i,j)+dzmr*(ap(m1-1,i,j)-ap(m1-2,i,j)))
          enddo
       enddo
    endif

    return
  end subroutine topset

  subroutine botset(m1,m2,m3,ia,iz,ja,jz,ibcon,aa,vnam)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon

    real, pointer, intent(in) :: aa(:,:,:)
    ! pointer intent(in), values intent(inout)
    character(len=*) :: vnam

    integer :: i,j
    real :: dzmr

    if (vnam .eq. 'P') then
       dzmr = dzm(2) / dzm(1)
       do i = 1,m2
          do j = 1,m3
             aa(1,i,j) = aa(2,i,j) + (aa(2,i,j) - aa(3,i,j)) * dzmr
          enddo
       enddo
    else
       do i = 1,m2
          do j = 1,m3
             aa(1,i,j) = aa(2,i,j)
          enddo
       enddo
    endif

    return
  end subroutine botset

  subroutine dumset(m1,m2,m3,ia,iz,ja,jz,ibcon,aa,vnam)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon

    real, pointer, intent(in) :: aa(:,:,:)
    ! pointer intent(in), values intent(inout)
    character(len=*) :: vnam

    integer :: i,j,k

    if (vnam .eq. 'U' .and. iand(ibcon,2) .gt. 0) then
       do j = 1,m3
          do k = 1,m1
             aa(k,m2,j) = aa(k,m2-1,j)
          enddo
       enddo
    elseif (vnam .eq. 'V' .and. iand(ibcon,8) .gt. 0) then
       do i = 1,m2
          do k = 1,m1
             aa(k,i,m3) = aa(k,i,m3-jdim)
          enddo
       enddo
    elseif (vnam .eq. 'W') then
       do j = 1,m3
          do i = 1,m2
             aa(m1,i,j) = aa(m1-1,i,j)
          enddo
       enddo
    endif

    return
  end subroutine dumset

  subroutine rayft(mxp,myp,mzp,mynum,ngrid,nnzp,if_adap,level,nodemyp,nodemxp,&
       vt3da, oneBasicFields)
    integer, intent(in)   :: mxp
    integer, intent(in)   :: myp
    integer, intent(in)   :: mzp
    integer, intent(in)   :: mynum
    integer, intent(in)   :: ngrid
    integer, intent(in)   :: nnzp(:)
    integer, intent(in)   :: if_adap
    integer, intent(in)   :: level  
    integer, intent(in)   :: nodemyp(:,:)
    integer, intent(in)   :: nodemxp(:,:)
    real   , pointer, intent(in)   :: vt3da(:)
    ! pointer intent(in), values intent(inout)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    
    integer :: ii,i,j,k

    integer(kind=i8) :: mxyzp, ind

    !     This routine is the rayleigh friction driver for the
    !     theta friction and is called from the long timestep.

    if (nfpt .eq. 0 .or. distim .eq. 0.) return

    mxyzp = mxp * myp * mzp

    !     First load past virtual theta into temporary.

    if (level .ge. 1) then
       ind = 0
       do j = 1,nodemyp(mynum,ngrid)
          do i = 1,nodemxp(mynum,ngrid)
             do k = 1,nnzp(ngrid)
                ind = ind + 1
                scratch%vt3da(ind) = oneBasicFields%theta(k,i,j)  &
                     * (1. + .61 * oneBasicFields%rv(k,i,j))
             enddo
          enddo
       enddo
    else
       call atob_long(mxyzp, oneBasicFields%theta, scratch%vt3da)
    endif

    !     Now get rayleigh friction tendency

    if (if_adap == 0) then

       call rayf(4,mzp,mxp,myp,ia,iz,ja,jz,ibcon                  &
            ,scratch%vt3da,oneBasicFields%th0  &
            ,tend%tht     ,grid_g(ngrid)%rtgt  &
            ,grid_g(ngrid)%topt)

    else

       call rayf_adap(4,mzp,mxp,myp,ia,iz,ja,jz,ibcon     &
            ,int(grid_g(ngrid)%lpw) ,scratch%vt3da  &
            ,oneBasicFields%th0 ,tend%tht      )

    endif
    return
  end subroutine rayft

  subroutine rayf(ifrom,m1,m2,m3,ia,iz,ja,jz,ibcon,var,th0,tht,rtgx,topx)
    integer, intent(in) :: ifrom
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: ibcon
    real, intent(inout), dimension(m1,m2,m3) :: var !TO recebe pointer 1D
    real, intent(inout), dimension(m1,m2,m3) :: tht !TO recebe pointer 1D

    real, pointer, intent(in) :: th0(:,:,:)
    ! pointer intent(in), values intent(in)
    real, pointer, intent(in) :: rtgx(:,:)
    ! pointer intent(in), values intent(in)
    real, pointer, intent(in) :: topx(:,:)
    ! pointer intent(in), values intent(in)


    real :: zmkf,c1,c2
    integer :: kf,i,j,k

    !     This routine calculates rayleigh friction terms velocity and theta_il

    if (nfpt .eq. 0 .or. distim .le. 0) return
    kf = nnz(1) - nfpt
    zmkf = zmn(kf,1)
    c1 = 1. / (distim * (ztop - zmkf))
    c2 = dts * c1
    goto(100,200,300,400) ifrom
100 continue

    !     u friction

    do j = ja,jz
       do i = ia,iz
          do k = 1,nzp
             vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
          enddo
          call htint(nzp,u01dn(1,ngrid),zt,nzp,vctr5,vctr2)
          do k = nz,2,-1
             if (vctr2(k) .le. zmkf) go to 10

             var(k,i,j) = var(k,i,j) + c2 * (vctr2(k) - zmkf)  &
                  * (vctr5(k) - var(k,i,j))

          enddo
10        continue
       enddo
    enddo
    return
200 continue

    !     V friction

    if (jdim .eq. 0 .and. icorflg .eq. 0) return
    do j = ja,jz
       do i = ia,iz
          do k = 1,nzp
             vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
          enddo
          call htint(nzp,v01dn(1,ngrid),zt,nzp,vctr5,vctr2)
          do k = nz,2,-1
             if (vctr2(k) .le. zmkf) go to 20
             var(k,i,j) = var(k,i,j) + c2 * (vctr2(k) - zmkf)  &
                  * (vctr5(k) - var(k,i,j))
          enddo
20        continue
       enddo
    enddo
    return
300 continue

    !     W friction

    do j = ja,jz
       do i = ia,iz
          do k = nz,2,-1
             vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
             if (vctr2(k) .le. zmkf) go to 30
             var(k,i,j) = var(k,i,j) - c2 * (vctr2(k) - zmkf) * var(k,i,j)
          enddo
30        continue
       enddo
    enddo
    return
400 continue

    !     THETA FRICTION

    do j = ja,jz
       do i = ia,iz
          do k = nz,2,-1
             vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
             if (vctr2(k) .le. zmkf) go to 40
             tht(k,i,j) = tht(k,i,j) + c1 * (vctr2(k) - zmkf)  &
                  * (th0(k,i,j) - var(k,i,j))
          enddo
40        continue
       enddo
    enddo
    return
  end subroutine rayf

  subroutine trsets(oneScalarTab, oneScalarTabSize, oneBasicFields, &
       oneTurbFields)
    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(inout) :: oneScalarTabSize
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields
    integer :: n, mxyzp
    character(len=*), parameter :: h="**(trsets)**"

    !     Apply lateral, top, and bottom boundary conditions.

    do n = 1,oneScalarTabSize
       if (nxtnest(ngrid) .eq. 0) then
          call LatSetScalar(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'TR'   &
               ,oneScalarTab(n)%var_p_3D,oneBasicFields%up  &
               ,oneBasicFields%vp ,grid_g(ngrid)%dxu   &
               ,grid_g(ngrid)%dxm ,grid_g(ngrid)%dyv   &
               ,grid_g(ngrid)%dym  )
       endif
       if (nsttop .eq. 1)  then

          if(n .le. (oneScalarTabSize - ( NADDSC + NSPECIES_TRANSPORTED) ))  then
             call TopSetScalar (mzp,mxp,myp,ia,iz,ja,jz,ibcon,oneScalarTab(n)%var_p_3D,'T')
          else
             call TopSet2Scalar(mzp,mxp,myp,ia,iz,ja,jz,ibcon,oneScalarTab(n)%var_p_3D,'T')
          endif
       endif
       if (nstbot .eq. 1)  then
          if (if_adap == 0) then
             call BotSetScalar(mzp,mxp,myp,ia,iz,ja,jz,ibcon,oneScalarTab(n)%var_p_3D,'T')
          else
             call BotSetAdapScalar(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
                  ,int(grid_g(ngrid)%lpw),oneScalarTab(n)%var_p_3D,'T')
          endif
       endif
    enddo

    !       Make sure all positive definite quantities remain such.

    call tkeinit(mzp,mxp,myp, oneTurbFields)

    call negadj1(mzp,mxp,myp, oneBasicFields)

    !--(DMK-CCATT-INI)-----------------------------------------------------
    !-srf for chem - aerosol quantities
    if( ccatt == 1 .and. CHEMISTRY >= 0) then
       mxyzp = mzp*mxp*myp
       do n = 1,oneScalarTabSize

          if(n .le. (oneScalarTabSize - ( NADDSC + NSPECIES_TRANSPORTED) )) cycle

          call KeepTracersNonnegScalar(mxyzp,oneScalarTab(n)%var_p_1D)

       enddo
    endif

  end subroutine trsets

  subroutine botset_adap(m1,m2,m3,ia,iz,ja,jz,ibcon,lpx,aa,vnam)
    integer, intent(in)    :: m1
    integer, intent(in)    :: m2
    integer, intent(in)    :: m3
    integer, intent(in)    :: ia
    integer, intent(in)    :: iz
    integer, intent(in)    :: ja
    integer, intent(in)    :: jz
    integer, intent(in)    :: ibcon
    integer, intent(in) :: lpx(m2,m3)
    real, intent(inout) :: aa(m1,m2,m3)
    character(len=*) :: vnam

    integer :: i,j,k,ka
    real :: dzmr

    if (vnam .eq. 'P') then
       dzmr = dzm(2) / dzm(1)
       do i = 1,m2
          do j = 1,m3
             ka = lpx(i,j)
             do k = ka-1,1,-1
                aa(k,i,j) = aa(k+1,i,j) + (aa(k+1,i,j) - aa(k+2,i,j)) * dzmr
             enddo
          enddo
       enddo
    else
       do i = 1,m2
          do j = 1,m3
             ka = lpx(i,j)
             do k = ka-1,1,-1
                aa(k,i,j) = aa(k+1,i,j)
             enddo
          enddo
       enddo
    endif

    if (vnam == 'U' .or. vnam == 'V' .or. vnam == 'W') then
       do i = 1,m2
          do j = 1,m3
             ka = lpx(i,j)
             do k = ka-1,1,-1
                aa(k,i,j) = 0.
             enddo
          enddo
       enddo
    endif

    return
  end subroutine botset_adap

  subroutine rayf_adap(ifrom,m1,m2,m3,ia,iz,ja,jz,ibcon,lpx,var,th0,tht)
    integer :: ifrom,m1,m2,m3,ia,iz,ja,jz,ibcon
    integer, dimension(m2,m3) :: lpx
    real, dimension(m1,m2,m3) :: var,th0,tht

    real :: zmkf,c1,c2
    integer :: kf,i,j,k

    !     This routine calculates rayleigh friction terms velocity and theta_il

    if (nfpt .eq. 0 .or. distim .le. 0) return
    kf = nnz(1) - nfpt
    zmkf = zmn(kf,1)
    c1 = 1. / (distim * (ztop - zmkf))
    c2 = dts * c1
    goto(100,200,300,400) ifrom
100 continue

    !     u friction
    do j = ja,jz
       do i = ia,iz
          do k = m1-1,lpx(i,j),-1
             if (zt(k) .le. zmkf) go to 10

             var(k,i,j) = var(k,i,j) + c2 * (zt(k) - zmkf)  &
                  * (u01dn(k,ngrid) - var(k,i,j))

          enddo
10        continue
       enddo
    enddo
    return
200 continue

    !     V friction

    if (jdim .eq. 0 .and. icorflg .eq. 0) return
    do j = ja,jz
       do i = ia,iz
          do k = m1-1,lpx(i,j),-1
             if (zt(k) .le. zmkf) go to 20
             var(k,i,j) = var(k,i,j) + c2 * (zt(k) - zmkf)  &
                  * (v01dn(1,ngrid) - var(k,i,j))
          enddo
20        continue
       enddo
    enddo
    return
300 continue

    !     W friction

    do j = ja,jz
       do i = ia,iz
          do k = m1-1,lpx(i,j),-1
             if (zt(k) .le. zmkf) go to 30
             var(k,i,j) = var(k,i,j) - c2 * (zt(k) - zmkf) * var(k,i,j)
          enddo
30        continue
       enddo
    enddo
    return
400 continue

    !     THETA FRICTION

    do j = ja,jz
       do i = ia,iz
          do k = m1-1,lpx(i,j),-1
             if (zt(k) .le. zmkf) go to 40
             tht(k,i,j) = tht(k,i,j) + c1 * (zt(k) - zmkf)  &
                  * (th0(k,i,j) - var(k,i,j))
          enddo
40        continue
       enddo
    enddo
    return
  end subroutine rayf_adap


end module ModRbnd

