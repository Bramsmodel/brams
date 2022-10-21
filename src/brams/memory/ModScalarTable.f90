!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModScalarTable

  use ModParallelEnvironment, only: &
       MsgDump
  
  implicit none
  private
  public :: ScalarTable
  public :: CreateScalarTab
  public :: DestroyScalarTab
  public :: DumpScalarTab
  public :: InsertAtScalarTab

  include "constants.h"

  ! Maximum number of ScalarTab entries

  integer, parameter :: maxvars=1000

  ! Define data type for scalar variable table

  type ScalarTable
     character (len=16) :: name
     real, contiguous, pointer      :: var_p_1D(:) => null()
     real, contiguous, pointer      :: var_p_2D(:,:) => null()
     real, contiguous, pointer      :: var_p_3D(:,:,:) => null()
     real, contiguous, pointer      :: var_t_1D(:) => null()
     real, contiguous, pointer      :: a_var_p_1D(:) => null()
     real, contiguous, pointer      :: a_var_t_1D(:) => null()
  end type ScalarTable

  interface InsertAtScalarTab
     module procedure InsertAtScalarTab_1D
     module procedure InsertAtScalarTab_2D
     module procedure InsertAtScalarTab_3D
  end interface InsertAtScalarTab

contains



  function CreateScalarTab() result(ST)
    type(ScalarTable), pointer :: ST(:)

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateScalarTab)**"
    logical, parameter :: dumpLocal=.false.

    allocate(ST(maxvars), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate ScalarTab fails with stat="//&
            trim(adjustl(str(1))))
    end if
  end function CreateScalarTab
    


  subroutine DestroyScalarTab(oneScalarTab)
    type(ScalarTable), pointer, intent(inout) :: oneScalarTab(:)

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyScalarTab)**"
    
    if (associated(oneScalarTab)) then
       deallocate(oneScalarTab, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    nullify(oneScalarTab)
  end subroutine DestroyScalarTab


  

  function StringScalarTableEntry(oneScalarTable) result(str)
    type(ScalarTable), intent(in) :: oneScalarTable
    character(len=128) :: str

    integer :: nAssoc=0

    str="variable "//trim(oneScalarTable%name)

    if (associated(oneScalarTable%var_p_1D)) then
       nAssoc = nAssoc + 1
       str = trim(str)//", var_p_1D"
    end if

    if (associated(oneScalarTable%var_p_2D)) then
       nAssoc = nAssoc + 1
       str = trim(str)//", var_p_2D"
    end if

    if (associated(oneScalarTable%var_p_3D)) then
       nAssoc = nAssoc + 1
       str = trim(str)//", var_p_3D"
    end if

    if (associated(oneScalarTable%a_var_p_1D)) then
       nAssoc = nAssoc + 1
       str = trim(str)//", a_var_p_1D"
    end if

    if (associated(oneScalarTable%a_var_t_1D)) then
       nAssoc = nAssoc + 1
       str = trim(str)//", a_var_t_1D"
    end if

    if (nAssoc == 0) then
       str = trim(str)//" has none associated entires"
    else
       str = trim(str)//" associated entires"
    end if
  end function StringScalarTableEntry





  subroutine DumpScalarTab(oneScalarTab, oneScalarTabSize)
    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(in) :: oneScalarTabSize

    integer :: i
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpScalarTab)**"

    write(str(2),"(i8)") oneScalarTabSize
    if (associated(oneScalarTab)) then
       write(str(1),"(i8)") size(oneScalarTab)
       call MsgDump(h//" associated with size "//trim(adjustl(str(1)))//&
            ", currently with "//trim(adjustl(str(2)))//" entries; dumping entries:")
       do i = 1, oneScalarTabSize
          call MsgDump(h//StringScalarTableEntry(oneScalarTab(i)))
       end do
    else
       call MsgDump(h//" not associated with #entries="//trim(adjustl(str(2))))
    end if
  end subroutine DumpScalarTab

  


  
  subroutine InsertAtScalarTab_1D(varp, vart, tabstr, ScalarTab, ScalarTabSize)
    real, contiguous, pointer, intent(in) :: varp(:)
    real, contiguous, pointer, intent(in) :: vart(:)
    character (len=*), intent(in) :: tabstr
    type(ScalarTable), pointer, intent(in) :: ScalarTab(:)
    integer, intent(inout) :: ScalarTabSize

    integer :: ntok
    character (len=1), parameter :: toksep=':'
    character (len=32) :: tokens(10)

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertAtScalarTab_2D)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       write(str(1),"(i8)") size(ScalarTab)
       write(str(2),"(i8)") ScalarTabSize
       call MsgDump(h//" inserting field "//trim(tabstr)//&
            " at ScalarTab with size "//trim(adjustl(str(1)))//&
            "; last ScalarTab used position is "//&
            trim(adjustl(str(2))))
    end if

    ! break tabstr into tokens, looking for variable name
    call tokenize1(tabstr, tokens, ntok, toksep)

    ! advance entries counter
    ScalarTabSize = ScalarTabSize + 1
    if (ScalarTabSize > maxvars) then
       write(str(1),"(i8)") maxvars
       call fatal_error(h//" allocated size ("//trim(adjustl(str(1)))//&
            ") exceeded while inserting "//trim(adjustl(tokens(1))))
    end if
    
    ! insert name
    ScalarTab(ScalarTabSize)%name = tokens(1)

    ! insert field and tendency
    ScalarTab(ScalarTabSize)%var_p_1D => varp
    ScalarTab(ScalarTabSize)%var_t_1D => vart
    ScalarTab(ScalarTabSize)%a_var_p_1D => varp
    ScalarTab(ScalarTabSize)%a_var_t_1D => vart

    if (dumpLocal) then
       write(str(1),"(i8)") ScalarTabSize
       call MsgDump(h//" inserted at ScalarTab position "//&
            trim(adjustl(str(1)))//" "//&
            StringScalarTableEntry(ScalarTab(ScalarTabSize)))
    end if
  end subroutine InsertAtScalarTab_1D





  
  subroutine InsertAtScalarTab_2D(varp, vart, tabstr, ScalarTab, ScalarTabSize)
    real, contiguous, pointer, intent(in) :: varp(:,:)
    real, contiguous, pointer, intent(in) :: vart(:)
    character (len=*), intent(in) :: tabstr
    type(ScalarTable), pointer, intent(in) :: ScalarTab(:)
    integer, intent(inout) :: ScalarTabSize

    integer :: ntok
    character (len=1), parameter :: toksep=':'
    character (len=32) :: tokens(10)

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertAtScalarTab_2D)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       write(str(1),"(i8)") size(ScalarTab)
       write(str(2),"(i8)") ScalarTabSize
       call MsgDump(h//" inserting field "//trim(tabstr)//&
            " at ScalarTab with size "//trim(adjustl(str(1)))//&
            "; last ScalarTab used position is "//&
            trim(adjustl(str(2))))
    end if

    ! break tabstr into tokens, looking for variable name
    call tokenize1(tabstr, tokens, ntok, toksep)

    ! advance entries counter
    ScalarTabSize = ScalarTabSize + 1
    if (ScalarTabSize > maxvars) then
       write(str(1),"(i8)") maxvars
       call fatal_error(h//" allocated size ("//trim(adjustl(str(1)))//&
            ") exceeded while inserting "//trim(adjustl(tokens(1))))
    end if

    ! insert name
    ScalarTab(ScalarTabSize)%name = tokens(1)

    ! insert field and tendency
    ScalarTab(ScalarTabSize)%var_p_2D => varp
    ScalarTab(ScalarTabSize)%var_t_1D => vart
    ScalarTab(ScalarTabSize)%a_var_p_1D(1:size(varp)) => varp(:,:)
    ScalarTab(ScalarTabSize)%a_var_t_1D => vart

    if (dumpLocal) then
       write(str(1),"(i8)") ScalarTabSize
       call MsgDump(h//" inserted at ScalarTab position "//&
            trim(adjustl(str(1)))//" "//&
            StringScalarTableEntry(ScalarTab(ScalarTabSize)))
    end if
  end subroutine InsertAtScalarTab_2D





  
  subroutine InsertAtScalarTab_3D(varp, vart, tabstr, ScalarTab, ScalarTabSize)
    real, contiguous, pointer, intent(in) :: varp(:,:,:)
    real, contiguous, pointer, intent(in) :: vart(:)
    character (len=*), intent(in) :: tabstr
    type(ScalarTable), pointer, intent(in) :: ScalarTab(:)
    integer, intent(inout) :: ScalarTabSize

    integer :: ntok
    character (len=1), parameter :: toksep=':'
    character (len=32) :: tokens(10)

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertAtScalarTab_2D)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       write(str(1),"(i8)") size(ScalarTab)
       write(str(2),"(i8)") ScalarTabSize
       call MsgDump(h//" inserting field "//trim(tabstr)//&
            " at ScalarTab with size "//trim(adjustl(str(1)))//&
            "; last ScalarTab used position is "//&
            trim(adjustl(str(2))))
    end if

    ! break tabstr into tokens, looking for variable name
    call tokenize1(tabstr, tokens, ntok, toksep)

    ! advance entries counter
    ScalarTabSize = ScalarTabSize + 1
    if (ScalarTabSize > maxvars) then
       write(str(1),"(i8)") maxvars
       call fatal_error(h//" allocated size ("//trim(adjustl(str(1)))//&
            ") exceeded while inserting "//trim(adjustl(tokens(1))))
    end if

    ! insert name
    ScalarTab(ScalarTabSize)%name = tokens(1)

    ! insert field and tendency
    ScalarTab(ScalarTabSize)%var_p_3D => varp
    ScalarTab(ScalarTabSize)%var_t_1D => vart
    ScalarTab(ScalarTabSize)%a_var_p_1D(1:size(varp)) => varp(:,:,:)
    ScalarTab(ScalarTabSize)%a_var_t_1D => vart

    if (dumpLocal) then
       write(str(1),"(i8)") ScalarTabSize
       call MsgDump(h//" inserted at ScalarTab position "//&
            trim(adjustl(str(1)))//" "//&
            StringScalarTableEntry(ScalarTab(ScalarTabSize)))
    end if
  end subroutine InsertAtScalarTab_3D
end module ModScalarTable
