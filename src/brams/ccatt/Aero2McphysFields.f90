!###########################################################################
!  B - Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModAero2McphysFields

  use ModParallelEnvironment, only: &
       MsgDump

  use ModNodeDimensions, only: &
       NodeDimensions
  
  use ModVarTable, only: &
       VarTable, &
       InsertVarTable

  implicit none

  private
  public :: Aero2McPhysFields
  public :: CreateAero2McphysFields
  public :: CreateEmptyAero2McphysFields
  public :: DestroyAero2McphysFields
  public :: DumpAero2McphysFields
  public :: InsertAero2McphysFieldsAtVarTable

  !-kml/srf - for microphysics activation

  type Aero2McphysFields
     !Inputs to matrix
     real, pointer, contiguous :: kappa_eff(:,:,:) => null()
     real, pointer, contiguous :: diam_eff(:,:,:) => null()
     real, pointer, contiguous :: numb_water(:,:,:) => null()
     real, pointer, contiguous :: numb_ice(:,:,:) => null()
  end type Aero2McphysFields

contains




  function CreateAero2McphysFields(oneNodeDimensions, nmodes, mcphys_type) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    integer, intent(in) :: nmodes
    integer, intent(in) :: mcphys_type
    type (Aero2McphysFields), pointer, contiguous :: res(:)

    integer :: i
    integer :: mxp
    integer :: myp
    integer :: mzp
    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateAero2McphysFields)**"

    nullify(res)

    if(mcphys_type == 3) then
       allocate(res(nmodes), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(2),"(i8)") nmodes
          call fatal_error(h//" allocate res("//&
               trim(adjustl(str(2)))//") fails with stat="//&
               trim(adjustl(str(1))))
       end if

       mzp=oneNodeDimensions%mzp
       mxp=oneNodeDimensions%mxp
       myp=oneNodeDimensions%myp

       do i=1,nmodes
          allocate(res(i)%kappa_eff(mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             write(str(2),"(i8)") i
             call fatal_error(h//" allocate kappa_eff for i="//&
                  trim(adjustl(str(2)))//" fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res(i)%kappa_eff  = 0.

          allocate(res(i)%diam_eff(mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             write(str(2),"(i8)") i
             call fatal_error(h//" allocate diam_eff for i="//&
                  trim(adjustl(str(2)))//" fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res(i)%diam_eff   = 0.

          allocate(res(i)%numb_water(mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             write(str(2),"(i8)") i
             call fatal_error(h//" allocate numb_water for i="//&
                  trim(adjustl(str(2)))//" fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res(i)%numb_water = 0.

          allocate(res(i)%numb_ice(mzp,mxp,myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             write(str(2),"(i8)") i
             call fatal_error(h//" allocate numb_ice for i="//&
                  trim(adjustl(str(2)))//"fails with stat="//&
                  trim(adjustl(str(1))))
          end if
          res(i)%numb_ice   = 0.
       enddo
    endif
  end function CreateAero2McphysFields





  function CreateEmptyAero2McphysFields(nmodes, mcphys_type) result(res)
    type (Aero2McphysFields), pointer, contiguous :: res(:)
    integer, intent(in) :: nmodes
    integer, intent(in) :: mcphys_type

    integer :: i
    integer :: mxp
    integer :: myp
    integer :: mzp
    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateEmptyAero2McphysFields)**"

    nullify(res)

    if(mcphys_type == 3) then
       allocate(res(nmodes), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(2),"(i8)") nmodes
          call fatal_error(h//" allocate res("//&
               trim(adjustl(str(2)))//") fails with stat="//&
               trim(adjustl(str(1))))
       end if

       do i=1,nmodes
          nullify(res(i)%kappa_eff)
          nullify(res(i)%diam_eff)
          nullify(res(i)%numb_water)
          nullify(res(i)%numb_ice)
       end do
    end if
  end function CreateEmptyAero2McphysFields




  
  subroutine DestroyAero2McphysFields(oneAero2McphysFields)
    type (Aero2McphysFields), pointer, intent(inout) :: oneAero2McphysFields(:)

    integer :: i
    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyAero2McphysFields)**"

    if (associated(oneAero2McphysFields)) then

       do i=1,size(oneAero2McphysFields)
          
          if (associated(oneAero2McphysFields(i)%kappa_eff)) then
             deallocate(oneAero2McphysFields(i)%kappa_eff, stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                write(str(2),"(i8)") i
                call fatal_error(h//" deallocate kappa_eff for i="//&
                     trim(adjustl(str(2)))//" fails with stat="//&
                     trim(adjustl(str(1))))
             end if
          end if
          
          if (associated(oneAero2McphysFields(i)%diam_eff)) then
             deallocate(oneAero2McphysFields(i)%diam_eff, stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                write(str(2),"(i8)") i
                call fatal_error(h//" deallocate diam_eff for i="//&
                     trim(adjustl(str(2)))//" fails with stat="//&
                     trim(adjustl(str(1))))
             end if
          end if

          if (associated(oneAero2McphysFields(i)%numb_water)) then
             deallocate(oneAero2McphysFields(i)%numb_water, stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                write(str(2),"(i8)") i
                call fatal_error(h//" deallocate numb_water for i="//&
                     trim(adjustl(str(2)))//" fails with stat="//&
                     trim(adjustl(str(1))))
             end if
          end if

          if (associated(oneAero2McphysFields(i)%numb_ice)) then
             deallocate(oneAero2McphysFields(i)%numb_ice, stat=ierr)
             if (ierr /= 0) then
                write(str(1),"(i8)") ierr
                write(str(2),"(i8)") i
                call fatal_error(h//" deallocate numb_ice for i="//&
                     trim(adjustl(str(2)))//"fails with stat="//&
                     trim(adjustl(str(1))))
             end if
          end if
       end do

       deallocate(oneAero2McphysFields, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(2),"(i8)") size(oneAero2McphysFields)
          call fatal_error(h//" deallocate oneAero2McphysFields("//&
               trim(adjustl(str(2)))//") fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    nullify(oneAero2McphysFields)
  end subroutine DestroyAero2McphysFields





  subroutine DumpAero2McphysFields(oneAero2McphysFields, name)
    type (Aero2McphysFields), pointer, intent(in) :: oneAero2McphysFields(:)
    character(len=*), intent(in) :: name

    integer :: i
    integer :: ierr
    character(len=8) :: str(10)
    logical :: assThis(4)
    character(len=*), parameter :: h="**(DumpAero2McphysFields)**"

    if (associated(oneAero2McphysFields)) then
       write(str(1),"(i8)") size(oneAero2McphysFields)
       assThis(1) = associated(oneAero2McphysFields(1)%kappa_eff)
       assThis(2) = associated(oneAero2McphysFields(1)%diam_eff)
       assThis(3) = associated(oneAero2McphysFields(1)%numb_water)
       assThis(4) = associated(oneAero2McphysFields(1)%numb_ice)
       write(str(2),"(i8)") count(assThis)
       call MsgDump(h//" oneAero2McphysFields from "//trim(adjustl(name))//&
            " is associated with size "//trim(adjustl(str(1)))//&
            " and each element has "//trim(adjustl(str(2)))//&
            " components associated")
    else
       call MsgDump(h//" oneAero2McphysFields from "//trim(adjustl(name))//&
            " is not associated")
    end if
  end subroutine DumpAero2McphysFields







  subroutine InsertAero2McphysFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneAero2McphysFields, oneAveAero2McphysFields, imean)

    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(Aero2McphysFields), pointer, intent(in) :: oneAero2McphysFields(:)
    type(Aero2McphysFields), pointer, intent(in) :: oneAveAero2McphysFields(:)
    integer, intent(in) :: imean

    integer :: i
    character(len=*), parameter :: h="**(InsertAero2McphysFieldsAtVarTable)**" 

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    end if

    ! Fill pointers to arrays into variable tables

    if (.not. associated(oneAero2McphysFields)) then
       return
    end if

    do i=1,size(oneAero2McphysFields)
       if (associated(oneAero2McphysFields(i)%kappa_eff)) then
          call InsertVarTable(oneVarTable, oneVarTableSize, &
               oneAero2McphysFields(i)%kappa_eff, &
               'KAPPA :3:hist:anal:mpti:mpt3', &
               oneAveAero2McphysFields(i)%kappa_eff, imean)
       end if

       if (associated(oneAero2McphysFields(i)%diam_eff)) then
          call InsertVarTable(oneVarTable, oneVarTableSize, &
               oneAero2McphysFields(i)%diam_eff, &
               'DIAMT_AER :3:hist:anal:mpti:mpt3', &
               oneAveAero2McphysFields(i)%diam_eff, imean)
       end if

       if (associated(oneAero2McphysFields(i)%numb_water)) then
          call InsertVarTable(oneVarTable, oneVarTableSize, &
               oneAero2McphysFields(i)%numb_water, &
               'WATER_FAER :3:hist:anal:mpti:mpt3', &
               oneAveAero2McphysFields(i)%numb_water, imean)
       end if

       if (associated(oneAero2McphysFields(i)%numb_ice)) then
          call InsertVarTable(oneVarTable, oneVarTableSize, &
               oneAero2McphysFields(i)%numb_ice, &
               'ICE_FAER :3:hist:anal:mpti:mpt3', &
               oneAveAero2McphysFields(i)%numb_ice, imean)
       end if
    end do
  end subroutine InsertAero2McphysFieldsAtVarTable
end module ModAero2McphysFields
