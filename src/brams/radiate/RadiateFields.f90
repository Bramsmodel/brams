!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module ModRadiateFields

  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModNamelistFile, only: &
       NamelistFile

  use ModVarTable, only: &
       VarTable, &
       InsertVarTable

  use ModNodeDimensions, only: &
       NodeDimensions
  
  implicit none

  private

  public :: RadiateFields
  public :: CreateRadiateFields
  public :: CreateEmptyRadiateFields
  public :: DestroyRadiateFields
  public :: DumpRadiateFields
  public :: InsertRadiateFieldsAtVarTable

  type RadiateFields
     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, contiguous :: fthrd(:,:,:) => null()
     real, pointer, contiguous :: cloud_fraction(:,:,:) => null()
     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, contiguous :: rshort(:,:) => null()
     real, pointer, contiguous :: rlong(:,:) => null()
     real, pointer, contiguous :: rlongup(:,:) => null()
     real, pointer, contiguous :: albedt(:,:) => null()
     real, pointer, contiguous :: cosz(:,:) => null()
  end type RadiateFields


contains




  function CreateRadiateFields(oneNodeDimensions, oneNamelistFile) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(RadiateFields), pointer :: res

    integer :: mzp
    integer :: mxp
    integer :: myp
    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateRadiateFields)**"

    ! Allocate arrays based on options (if necessary)

    if (oneNamelistFile%ilwrtyp + oneNamelistFile%iswrtyp > 0)  then

       mzp=oneNodeDimensions%mzp
       mxp=oneNodeDimensions%mxp
       myp=oneNodeDimensions%myp
       
       allocate (res, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate res fails with stat="//&
               trim(adjustl(str(1))))
       end if

       allocate (res%fthrd(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate fthrd fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%fthrd=0.0

       allocate (res%cloud_fraction(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate cloud_fraction fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%cloud_fraction=0.0

       allocate (res%rshort(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate rshort fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%rshort=0.0

       allocate (res%rlong(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate rlong fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%rlong=0.0

       allocate (res%rlongup(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate rlongup fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%rlongup=0.0

       allocate (res%albedt(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate albedt fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%albedt=0.0

       allocate (res%cosz(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate cosz fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%cosz=0.0
    end if
  end function CreateRadiateFields




  function CreateEmptyRadiateFields(oneNamelistFile) result(res)
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(RadiateFields), pointer :: res

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateEmptyRadiateFields)**"

    ! Allocate arrays based on options (if necessary)

    if (oneNamelistFile%ilwrtyp + oneNamelistFile%iswrtyp > 0)  then

       allocate (res, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate res fails with stat="//&
               trim(adjustl(str(1))))
       end if

       nullify(res%fthrd)
       nullify(res%cloud_fraction)
       nullify(res%rshort)
       nullify(res%rlong)
       nullify(res%rlongup)
       nullify(res%albedt)
       nullify(res%cosz)

    end if
  end function CreateEmptyRadiateFields





  subroutine DestroyRadiateFields(oneRadiateFields)
    type(RadiateFields), pointer, intent(inout) :: oneRadiateFields

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyRadiateFields)**"

    ! Allocate arrays based on options (if necessary)

    if (associated(oneRadiateFields)) then

       if (associated(oneRadiateFields%fthrd)) then
          deallocate (oneRadiateFields%fthrd, stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" deallocate fthrd fails with stat="//&
                  trim(adjustl(str(1))))
          end if
       end if

       if (associated(oneRadiateFields%cloud_fraction)) then
          deallocate (oneRadiateFields%cloud_fraction, stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" deallocate cloud_fraction fails with stat="//&
                  trim(adjustl(str(1))))
          end if
       end if

       if (associated(oneRadiateFields%rshort)) then
          deallocate (oneRadiateFields%rshort, stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" deallocate rshort fails with stat="//&
                  trim(adjustl(str(1))))
          end if
       end if

       if (associated(oneRadiateFields%rlong)) then
          deallocate (oneRadiateFields%rlong, stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" deallocate rlong fails with stat="//&
                  trim(adjustl(str(1))))
          end if
       end if

       if (associated(oneRadiateFields%rlongup)) then
          deallocate (oneRadiateFields%rlongup, stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" deallocate rlongup fails with stat="//&
                  trim(adjustl(str(1))))
          end if
       end if

       if (associated(oneRadiateFields%albedt)) then
          deallocate (oneRadiateFields%albedt, stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" deallocate albedt fails with stat="//&
                  trim(adjustl(str(1))))
          end if
       end if

       if (associated(oneRadiateFields%cosz)) then
          deallocate (oneRadiateFields%cosz, stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" deallocate cosz fails with stat="//&
                  trim(adjustl(str(1))))
          end if
       end if

       deallocate (oneRadiateFields, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneRadiateFields fails with stat="//&
               trim(adjustl(str(1))))
       end if


    end if

    nullify (oneRadiateFields)

  end subroutine DestroyRadiateFields





  subroutine DumpRadiateFields(oneRadiateFields, name)
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpRadiateFields)**"

    if (associated(oneRadiateFields)) then
       call MsgDump(h//" "//name//" is associated")
    else
       call MsgDump(h//" "//name//" is not associated")
    end if
  end subroutine DumpRadiateFields





  subroutine InsertRadiateFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneRadiateFields, oneAveRadiateFields, imean)

    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields
    type(RadiateFields), pointer, intent(in) :: oneAveRadiateFields
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(InsertRadiateFieldsAtVarTable)**" 

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    end if

    if (associated(oneRadiateFields)) then

       if (.not. associated(oneAveRadiateFields)) then
          call fatal_error(h//" oneAveRadiateFields not associated")
       end if
       
       ! Fill pointers to arrays into variable tables

       if (associated(oneRadiateFields%cloud_fraction))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneRadiateFields%cloud_fraction, &
               'CLOUD_FRACTION :3:anal:mpti:mpt3', &
               oneAveRadiateFields%cloud_fraction, imean)
       end if

       if (associated(oneRadiateFields%fthrd))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneRadiateFields%fthrd, &
               'FTHRD :3:hist:anal:mpti:mpt3', &
               oneAveRadiateFields%fthrd, imean)
       end if

       if (associated(oneRadiateFields%rshort))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneRadiateFields%rshort, &
               'RSHORT :2:hist:anal:mpti:mpt3', &
               oneAveRadiateFields%rshort, imean)
       end if

       if (associated(oneRadiateFields%rlong))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneRadiateFields%rlong, &
               'RLONG :2:hist:anal:mpti:mpt3', &
               oneAveRadiateFields%rlong, imean)
       end if

       if (associated(oneRadiateFields%rlongup))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneRadiateFields%rlongup, &
               'RLONGUP :2:hist:anal:mpti:mpt3', &
               oneAveRadiateFields%rlongup, imean)
       end if

       if (associated(oneRadiateFields%albedt))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneRadiateFields%albedt, &
               'ALBEDT :2:hist:anal:mpti:mpt3', &
               oneAveRadiateFields%albedt, imean)
       end if

       if (associated(oneRadiateFields%cosz))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneRadiateFields%cosz, &
               'COSZ :2:hist:anal:mpt3', &
               oneAveRadiateFields%cosz, imean)
       end if
    end if
  end subroutine InsertRadiateFieldsAtVarTable

end module ModRadiateFields
