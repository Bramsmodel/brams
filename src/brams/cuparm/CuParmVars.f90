!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModCuParmVars

  use ModParallelEnvironment, only: &
       MsgDump

  use ModNamelistFile, only: &
       NamelistFile

  use grid_dims, only: &
       maxfiles          ! INTENT(IN)

  implicit none

  include "files.h"
  
  private
  public :: CuParmVars
  public :: CreateCuParmVars
  public :: DestroyCuParmVars
  public :: DumpCuParmVars
  
  type CuParmVars
     character(len=f_name_length) :: fnames_cu(maxfiles)
     character(len=14)  :: itotdate_cu(maxfiles)
     real :: cu_times(maxfiles)
     integer :: ncufiles
     integer :: ncufl
     real :: cutime1
     real :: cutime2
  end type CuParmVars
  

contains




  function CreateCuParmVars() result(res)
    type(CuParmVars), pointer :: res

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateCuParmVars)**"

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(str(1), "(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if

    res%fnames_cu=""
    res%itotdate_cu=""
    res%cu_times=0.0
    res%ncufiles=0
    res%ncufl=0
    res%cutime1=0.0
    res%cutime2=0.0
  end function CreateCuParmVars




  subroutine DestroyCuParmVars(oneCuParmVars)
    type(CuParmVars), pointer, intent(inout) :: oneCuParmVars

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyCuParmVars)**"

    if (associated(oneCuParmVars)) then
       deallocate(oneCuParmVars, stat=ierr)
       if (ierr /= 0) then
          write(str(1), "(i8)") ierr
          call fatal_error(h//" deallocate res fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    nullify(oneCuParmVars)
  end subroutine DestroyCuParmVars




  subroutine DumpCuParmVars(oneCuParmVars, name)
    type(CuParmVars), pointer, intent(in) :: oneCuParmVars
    character(len=*), intent(in) :: name

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpCuParmVars)**"

    if (associated(oneCuParmVars)) then
       call MsgDump(h//" oneCuParmVars from "//trim(name)//" is associated")
    else
       call MsgDump(h//" oneCuParmVars from "//trim(name)//" is not associated")
    end if
  end subroutine DumpCuParmVars

end module ModCuParmVars
