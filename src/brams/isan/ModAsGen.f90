!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module ModAsGen

  use mem_grid, only: ngrids 
  use isan_coms, only: nigrids

  implicit none

  private

  public :: opspec4

contains
  subroutine opspec4
    integer :: ifaterr,iwarerr,infoerr
    character(len=8) :: str(10)
    character, parameter :: h="**(opspec4)**"
    
    ! This routine checks the option specifications in the $MODEL_GRIDS
    !    $ISAN_ISENTROPIC namelists for consistency.

    ! don't allow nigrids <= ngrids

    if(nigrids.gt.ngrids)then
       write(str(1),"(i8)") nigrids
       write(str(2),"(i8)") ngrids
       call fatal_error(h//" nigrids ("//trim(adjustl(str(1)))//&
            ") must be <= ngrids ("//trim(adjustl(str(2)))//")")
       ifaterr=ifaterr+1
    endif

    ! make sure that timmax >= isan_inc

    !fisan_inc=int(isan_inc/100)*3600+float(mod(isan_inc,100))
    !if(fisan_inc.gt.timmax)then
    !  print*,' warning - timmax must be <= isan_inc'
    !  print*,'           resetting to timmax to (s) ',fisan_inc
    !  timmax=fisan_inc
    !  timstr=fisan_inc
    !  iwarerr=iwarerr+1
    !endif
  end subroutine opspec4
end module ModAsGen
