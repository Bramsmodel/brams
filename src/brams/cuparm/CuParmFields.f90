!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModCuParmFields

  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModNamelistFile, only: &
       NamelistFile

  use ModNodeDimensions, only: &
       NodeDimensions

  use ModVarTable, only: &
       VarTable, &
       InsertVarTable

  implicit none

  include "files.h"

  private

  public :: CuParmFields

  public :: hasAconpr

  public :: CreateCuParmFields
  public :: CreateEmptyCuParmFields
  public :: DestroyCuParmFields
  public :: DumpCuParmFields
  public :: InsertCuParmFieldsAtVarTable

  public :: CreateCuParmShFields
  public :: CreateEmptyCuParmShFields
  public :: DestroyCuParmShFields
  public :: DumpCuParmShFields
  public :: InsertCuParmShFieldsAtVarTable

  type CuParmFields
     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, contiguous :: thsrc(:,:,:) => null()
     real, pointer, contiguous :: rtsrc(:,:,:) => null()
     real, pointer, contiguous :: thsrcf(:,:,:) => null()
     real, pointer, contiguous :: rtsrcf(:,:,:) => null()
     real, pointer, contiguous :: thsrcp(:,:,:) => null()
     real, pointer, contiguous :: rtsrcp(:,:,:) => null()
     real, pointer, contiguous :: clsrc(:,:,:) => null() !srf-cloud/ice source term
     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, contiguous :: aconpr(:,:) => null()
     real, pointer, contiguous :: conprr(:,:) => null()
     real, pointer, contiguous :: conprrp(:,:) => null()
     real, pointer, contiguous :: conprrf(:,:) => null()
  end type CuParmFields



contains



  logical function hasAconpr(oneNamelistFile, ng)
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    integer, intent(IN) :: ng

    hasAconpr = oneNamelistFile%nnqparm(ng)>= 1 .or. oneNamelistFile%if_cuinv == 1

  end function hasAconpr



  function CreateCuParmFields(oneNamelistFile, oneNodeDimensions, gridId) result(res)
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    integer, intent(in) :: gridId
    type(CuParmFields), pointer :: res

    integer :: mzp
    integer :: mxp
    integer :: myp
    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateCuParmFields)**"

    ! Allocate arrays based on options (if necessary)

    if( hasAconpr(oneNamelistFile, gridId) )  then

       mzp=oneNodeDimensions%mzp
       mxp=oneNodeDimensions%mxp
       myp=oneNodeDimensions%myp

       allocate (res, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate res fails with stat="//&
               trim(adjustl(str(1))))
       end if

       if(oneNamelistFile%nnqparm(gridId) < 3) then
          !srf: for conv 3 and 4, special arrays
          !srf:  are allocated instead of these ones

          allocate (res%thsrc(mzp,mxp,myp), stat=ierr)
          if (ierr/=0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate thsrc fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%thsrc=0.0

          allocate (res%rtsrc(mzp,mxp,myp), stat=ierr)
          if (ierr/=0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate rtsrc fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%rtsrc=0.0

          allocate (res%clsrc(mzp,mxp,myp), stat=ierr)
          if (ierr/=0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate clsrc fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%clsrc=0.0 !srf-cloud/ice source term

       end if

       allocate (res%aconpr(mxp,myp), stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate aconpr fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%aconpr = 0.

       allocate (res%conprr(mxp,myp), stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate conprr fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%conprr = 0.


       if (oneNamelistFile%if_cuinv == 1) then

          allocate (res%thsrcp(mzp,mxp,myp), stat=ierr)
          if (ierr/=0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate thsrcp fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%thsrcp=0.0

          allocate (res%rtsrcp(mzp,mxp,myp), stat=ierr)
          if (ierr/=0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate rtsrcp fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%rtsrcp=0.0

          allocate (res%thsrcf(mzp,mxp,myp), stat=ierr)
          if (ierr/=0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate thsrcf fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%thsrcf=0.0

          allocate (res%rtsrcf(mzp,mxp,myp), stat=ierr)
          if (ierr/=0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate rtsrcf fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%rtsrcf=0.0

          allocate (res%conprrp(mxp,myp), stat=ierr)
          if (ierr/=0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate conprp fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%conprrp=0.0

          allocate (res%conprrf(mxp,myp), stat=ierr)
          if (ierr/=0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate conprrf fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res%conprrf=0.0

       end if
    else
       nullify(res)
    end if

  end function CreateCuParmFields




  function CreateEmptyCuParmFields() result(res)
    type(CuParmFields), pointer :: res

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateEmptyCuParmFields)**"

    allocate (res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if

    nullify(res%thsrc)
    nullify(res%rtsrc)
    nullify(res%clsrc)
    nullify(res%aconpr)
    nullify(res%conprr)
    nullify(res%thsrcp)
    nullify(res%rtsrcp)
    nullify(res%thsrcf)
    nullify(res%rtsrcf)
    nullify(res%conprrp)
    nullify(res%conprrf)
    
  end function CreateEmptyCuParmFields




  subroutine DestroyCuParmFields(oneCuParmFields)
    type(CuParmFields), pointer, intent(inout) :: oneCuParmFields

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyCuParmFields)**"

    if (.not. associated(oneCuParmFields)) then
       nullify(oneCuParmFields)
       return
    end if

    if (associated(oneCuParmFields%thsrc)) then
       deallocate (oneCuParmFields%thsrc, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate thsrc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%rtsrc)) then
       deallocate (oneCuParmFields%rtsrc, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rtsrc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%clsrc)) then
       deallocate (oneCuParmFields%clsrc, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate clsrc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%aconpr)) then
       deallocate (oneCuParmFields%aconpr, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate aconpr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%conprr)) then
       deallocate (oneCuParmFields%conprr, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate conprr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%thsrcp)) then
       deallocate (oneCuParmFields%thsrcp, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate thsrcp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%rtsrcp)) then
       deallocate (oneCuParmFields%rtsrcp, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rtsrcp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%thsrcf)) then
       deallocate (oneCuParmFields%thsrcf, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate thsrcf fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%rtsrcf)) then
       deallocate (oneCuParmFields%rtsrcf, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rtsrcf fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%conprrp)) then
       deallocate (oneCuParmFields%conprrp, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate conprp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%conprrf)) then
       deallocate (oneCuParmFields%conprrf, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate conprrf fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    deallocate (oneCuParmFields, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneCuParmFields fails with stat="//&
            trim(adjustl(str(1))))
    end if
    nullify(oneCuParmFields)
  end subroutine DestroyCuParmFields



  subroutine DumpCuParmFields(oneCuParmFields, name)
    type(CuParmFields), pointer, intent(in) :: oneCuParmFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpCuParmFields)**"

    if (associated(oneCuParmFields)) then
       call MsgDump(h//" "//trim(adjustl(name))//" is associated")
    else
       call MsgDump(h//" "//trim(adjustl(name))//" is not associated")
    end if
  end subroutine DumpCuParmFields




  subroutine InsertCuParmFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneCuParmFields, oneAveCuParmFields, imean)

    ! Arguments:
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type (CuParmFields), pointer, intent(in) :: oneCuParmFields
    type (CuParmFields), pointer, intent(in) :: oneAveCuParmFields
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(InsertCuParmFieldsAtVarTable)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" invoked with unasociated oneVarTable")
    end if

    if (.not. associated(oneCuParmFields)) then
       return
    else if (.not. associated(oneAveCuParmFields)) then
       call fatal_error(h//&
            " oneCuParmFields is associated but oneAveCuParmFields is not associated")
    end if

    ! Fill pointers to arrays into variable tables

    if (associated(oneCuParmFields%thsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%thsrc, 'THSRC :3:hist:anal:mpti:mpt3', &
            oneAveCuParmFields%thsrc, imean)
    end if

    if (associated(oneCuParmFields%rtsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%rtsrc, 'RTSRC :3:hist:anal:mpti:mpt3', &
            oneAveCuParmFields%rtsrc, imean)
    end if

    if (associated(oneCuParmFields%clsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%clsrc, 'CLSRC :3:hist:anal:mpti:mpt3', &
            oneAveCuParmFields%clsrc, imean)
    end if

    if (associated(oneCuParmFields%thsrcp)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%thsrcp, 'THSRCP :3:mpti:', &
            oneAveCuParmFields%thsrcp, imean)
    end if

    if (associated(oneCuParmFields%rtsrcp)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%rtsrcp, 'RTSRCP :3:mpti:', &
            oneAveCuParmFields%rtsrcp, imean)
    end if

    if (associated(oneCuParmFields%thsrcf)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%thsrcf, 'THSRCF :3:mpti:', &
            oneAveCuParmFields%thsrcf, imean)
    end if

    if (associated(oneCuParmFields%rtsrcf)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%rtsrcf, 'RTSRCF :3:mpti:', &
            oneAveCuParmFields%rtsrcf, imean)
    end if

    if (associated(oneCuParmFields%aconpr)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%aconpr, 'ACONPR :2:hist:anal:mpti:mpt3', &
            oneAveCuParmFields%aconpr, imean)
    end if

    if (associated(oneCuParmFields%conprr)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%conprr, 'CONPRR :2:hist:anal:mpt3', &
            oneAveCuParmFields%conprr, imean)
    end if

    if (associated(oneCuParmFields%conprrp)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%conprrp, 'CONPRRP :2:mpti', &
            oneAveCuParmFields%conprrp, imean)
    end if

    if (associated(oneCuParmFields%conprrf)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%conprrf, 'CONPRRF :2:mpti', &
            oneAveCuParmFields%conprrf, imean)
    end if

  end subroutine InsertCuParmFieldsAtVarTable





  function CreateCuParmShFields(oneNodeDimensions) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(CuParmFields), pointer :: res

    integer :: mzp
    integer :: mxp
    integer :: myp
    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateCuParmShFields)**"

    ! Allocate arrays for shallow cummulus feedback

    mzp=oneNodeDimensions%mzp
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp

    allocate (res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if

    allocate (res%thsrc(mzp,mxp,myp), stat=ierr)
    if (ierr/=0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate thsrc fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%thsrc=0.0

    allocate (res%rtsrc(mzp,mxp,myp), stat=ierr)
    if (ierr/=0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate rtsrc fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%rtsrc=0.0

    nullify(res%clsrc)
    nullify(res%aconpr)
    nullify(res%conprr)
    nullify(res%thsrcp)
    nullify(res%rtsrcp)
    nullify(res%thsrcf)
    nullify(res%rtsrcf)
    nullify(res%conprrp)
    nullify(res%conprrf)
  end function CreateCuParmShFields





  function CreateEmptyCuParmShFields() result(res)
    type(CuParmFields), pointer :: res

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateEmptyCuParmShFields)**"

    ! Allocate arrays for shallow cummulus feedback

    allocate (res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if

    nullify(res%thsrc)
    nullify(res%rtsrc)
    nullify(res%clsrc)
    nullify(res%aconpr)
    nullify(res%conprr)
    nullify(res%thsrcp)
    nullify(res%rtsrcp)
    nullify(res%thsrcf)
    nullify(res%rtsrcf)
    nullify(res%conprrp)
    nullify(res%conprrf)
  end function CreateEmptyCuParmShFields




  subroutine DestroyCuParmShFields(oneCuParmFields)
    type(CuParmFields), pointer, intent(inout) :: oneCuParmFields

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyCuParmShFields)**"

    if (.not. associated(oneCuParmFields)) then
       nullify(oneCuParmFields)
       return
    end if

    if (associated(oneCuParmFields%thsrc)) then
       deallocate (oneCuParmFields%thsrc, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate thsrc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneCuParmFields%rtsrc)) then
       deallocate (oneCuParmFields%rtsrc, stat=ierr)
       if (ierr/=0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rtsrc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    deallocate (oneCuParmFields, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneCuParmFields fails with stat="//&
            trim(adjustl(str(1))))
    end if
    nullify(oneCuParmFields)
  end subroutine DestroyCuParmShFields




  subroutine DumpCuParmShFields(oneCuParmFields, name)
    type(CuParmFields), pointer, intent(in) :: oneCuParmFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpCuParmShFields)**"

    if (associated(oneCuParmFields)) then
       call MsgDump(h//" "//trim(adjustl(name))//" is associated")
    else
       call MsgDump(h//" "//trim(adjustl(name))//" is not associated")
    end if
  end subroutine DumpCuParmShFields





  subroutine InsertCuParmShFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneCuParmFields, oneAveCuParmFields, imean)

    ! Arguments:
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type (CuParmFields), pointer, intent(in) :: oneCuParmFields
    type (CuParmFields), pointer, intent(in) :: oneAveCuParmFields
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(InsertCuParmShFieldsAtVarTable)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" invoked with unasociated oneVarTable")
    end if

    if (.not. associated(oneCuParmFields)) then
       return
    else if (.not. associated(oneAveCuParmFields)) then
       call fatal_error(h//&
            " oneCuParmFields is associated but oneAveCuParmFields is not associated")
    end if

    ! Fill pointers to arrays into variable tables

    if (associated(oneCuParmFields%thsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%thsrc, 'THSRC :3:hist:anal:mpti:mpt3', &
            oneAveCuParmFields%thsrc, imean)
    end if

    if (associated(oneCuParmFields%rtsrc)) then
       call InsertVarTable(oneVarTable, oneVarTableSize, &
            oneCuParmFields%rtsrc, 'RTSRC :3:hist:anal:mpti:mpt3', &
            oneAveCuParmFields%rtsrc, imean)
    end if

  end subroutine InsertCuParmShFieldsAtVarTable
end module ModCuParmFields
