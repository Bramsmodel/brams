!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module ModControlVars

  use ModNamelistFile, only: &
       NamelistFile

  implicit none

  include "files.h"

  private
  public :: ControlVars
  public :: CreateControlVars
  public :: DestroyControlVars


  type ControlVars

     ! namelist variables indexed by grid number
     ! copied to this grid scalar vars, to avoid
     ! indexing by grid number

     integer :: nnxp
     integer :: nnyp
     integer :: nnzp
     integer :: nxtnest
     integer :: nstratx
     integer :: nstraty
     integer :: nndtrat
     real :: centlat
     real :: centlon
     integer :: ninest
     integer :: njnest
     integer :: nknest
     integer :: nnsttop
     integer :: nnstbot
     real :: gridu
     real :: gridv
     integer :: ifusflg
     character(len=f_name_length) :: ifusfn
     real :: wt_nudge_grid
     real :: wt_nudgec_grid
     integer :: itoptflg
     integer :: isstflg
     integer :: ivegtflg
     integer :: isoilflg
     integer :: ndviflg
     integer :: nofilflg
     character(len=f_name_length) :: itoptfn
     character(len=f_name_length) :: isstfn
     character(len=f_name_length) :: ivegtfn
     character(len=f_name_length) :: isoilfn
     character(len=f_name_length) :: ndvifn
     integer :: itopsflg
     real :: toptenh
     real :: toptwvl
     integer :: iz0flg
     real :: z0max
     integer :: nnqparm
     integer :: nnshcu
     integer :: idiffk
     real :: csx
     real :: csz
     real :: xkhkm
     real :: zkhkm
     real :: akmin
     real :: lati
     real :: loni
     real :: latf
     real :: lonf
     integer :: zlevmax
  end type ControlVars


contains



  function CreateControlVars(oneNamelistFile, gridId) result(res)
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    integer :: gridId
    type(ControlVars), pointer :: res

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateControlVars)**"

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if

    ! namelist variables
    
    res%akmin = oneNamelistFile%akmin(gridId)
    res%centlat = oneNamelistFile%centlat(gridId)
    res%centlon = oneNamelistFile%centlon(gridId)
    res%csx = oneNamelistFile%csx(gridId)
    res%csz = oneNamelistFile%csz(gridId)
    res%gridu = oneNamelistFile%gridu(gridId)
    res%gridv = oneNamelistFile%gridv(gridId)
    res%idiffk = oneNamelistFile%idiffk(gridId)
    res%ifusflg = oneNamelistFile%ifusflg(gridId)
    res%ifusfn = oneNamelistFile%ifusfn(gridId)
    res%isoilflg = oneNamelistFile%isoilflg(gridId)
    res%isoilfn = oneNamelistFile%isoilfn(gridId)
    res%isstflg = oneNamelistFile%isstflg(gridId)
    res%isstfn = oneNamelistFile%isstfn(gridId)
    res%itopsflg = oneNamelistFile%itopsflg(gridId)
    res%itoptflg = oneNamelistFile%itoptflg(gridId)
    res%itoptfn = oneNamelistFile%itoptfn(gridId)
    res%ivegtflg = oneNamelistFile%ivegtflg(gridId)
    res%ivegtfn = oneNamelistFile%ivegtfn(gridId)
    res%iz0flg = oneNamelistFile%iz0flg(gridId)
    res%latf = oneNamelistFile%latf(gridId)
    res%lati = oneNamelistFile%lati(gridId)
    res%lonf = oneNamelistFile%lonf(gridId)
    res%loni = oneNamelistFile%loni(gridId)
    res%ndviflg = oneNamelistFile%ndviflg(gridId)
    res%ndvifn = oneNamelistFile%ndvifn(gridId)
    res%ninest = oneNamelistFile%ninest(gridId)
    res%njnest = oneNamelistFile%njnest(gridId)
    res%nknest = oneNamelistFile%nknest(gridId)
    res%nndtrat = oneNamelistFile%nndtrat(gridId)
    res%nnqparm = oneNamelistFile%nnqparm(gridId)
    res%nnshcu = oneNamelistFile%nnshcu(gridId)
    res%nnstbot = oneNamelistFile%nnstbot(gridId)
    res%nnsttop = oneNamelistFile%nnsttop(gridId)
    res%nnxp = oneNamelistFile%nnxp(gridId)
    res%nnyp = oneNamelistFile%nnyp(gridId)
    res%nnzp = oneNamelistFile%nnzp(gridId)
    res%nofilflg = oneNamelistFile%nofilflg(gridId)
    res%nstratx = oneNamelistFile%nstratx(gridId)
    res%nstraty = oneNamelistFile%nstraty(gridId)
    res%nxtnest = oneNamelistFile%nxtnest(gridId)
    res%toptenh = oneNamelistFile%toptenh(gridId)
    res%toptwvl = oneNamelistFile%toptwvl(gridId)
    res%wt_nudge_grid = oneNamelistFile%wt_nudge_grid(gridId)
    res%wt_nudgec_grid = oneNamelistFile%wt_nudgec_grid(gridId)
    res%xkhkm = oneNamelistFile%xkhkm(gridId)
    res%z0max = oneNamelistFile%z0max(gridId)
    res%zkhkm = oneNamelistFile%zkhkm(gridId)
    res%zlevmax = oneNamelistFile%zlevmax(gridId)
  end function CreateControlVars







  subroutine DestroyControlVars(oneControlVars)
    type(ControlVars), pointer, intent(inout) :: oneControlVars

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyControlVars)**"
    
    if (associated(oneControlVars)) then
       deallocate(oneControlVars, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneControlVars fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    nullify(oneControlVars)
  end subroutine DestroyControlVars

end module ModControlVars
