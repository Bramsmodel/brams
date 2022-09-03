!!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModScalarFields

  ! Added scalar variables and tendencies

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
  public :: ScalarFields
  public :: CreateScalarFields
  public :: CreateEmptyScalarFields
  public :: DestroyScalarFields
  public :: DumpScalarFields
  public :: InsertScalarFieldsAtVarTable 

  type ScalarFields
     real, contiguous, pointer :: sclp(:,:,:) => null()
     real, contiguous, pointer :: drydep(:,:) => null()
     real, contiguous, pointer :: sclt(:) => null()
     real, contiguous, pointer :: wetdep(:,:) => null()
     real, contiguous, pointer :: srcsc(:,:,:) => null()
  end type ScalarFields


contains



  
  function CreateScalarFields(oneNodeDimensions, oneNamelistFile, gridId) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(ScalarFields), pointer :: res(:)
    integer, intent(in) :: gridId

    integer :: ierr
    integer :: nsc
    integer :: naddsc
    integer :: mzp
    integer :: mxp
    integer :: myp
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateScalarFields)**" 

    if (.not. associated(oneNodeDimensions)) then
       call fatal_error(h//" oneNodeDimensions not associated")
    else if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if

    naddsc = oneNamelistFile%naddsc
    if (naddsc == 0) then
       nullify(res)
       return
    end if

    mzp=oneNodeDimensions%mzp
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp
    
    allocate(res(naddsc), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       write(str(2),"(i8)") naddsc
       call fatal_error(h//" allocate res("//trim(adjustl(str(2)))//&
            " fails with stat="//trim(adjustl(str(1))))
    end if
    
    do nsc=1,naddsc
       allocate (res(nsc)%sclp(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(2),"(i8)") mzp
          write(str(3),"(i8)") mxp
          write(str(4),"(i8)") myp
          call fatal_error(h//" allocate sclp("//&
               trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//","//&
               trim(adjustl(str(4)))//")"//&
               " fails with stat="//trim(adjustl(str(1))))
       end if
       res(nsc)%sclp=0.0
       
       allocate (res(nsc)%drydep(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(3),"(i8)") mxp
          write(str(4),"(i8)") myp
          call fatal_error(h//" allocate drydep("//&
               trim(adjustl(str(3)))//","//&
               trim(adjustl(str(4)))//")"//&
               " fails with stat="//trim(adjustl(str(1))))
       end if
       res(nsc)%drydep=0.0

       ! single tendency for all grids
       
       if (gridId == 1) then
          allocate (res(nsc)%sclt(mzp*mxp*myp), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             write(str(2),"(i8)") mzp
             write(str(3),"(i8)") mxp
             write(str(4),"(i8)") myp
             call fatal_error(h//" allocate sclt("//&
                  trim(adjustl(str(2)))//"*"//&
                  trim(adjustl(str(3)))//"*"//&
                  trim(adjustl(str(4)))//")"//&
                  " fails with stat="//trim(adjustl(str(1))))
          end if
          res(nsc)%sclt=0.0
       else
          call fatal_error(h//" tendencies (sclt) for multiple grids "//&
               "share the first grid area; not implemented")
!!$          do nsc=1,naddsc
!!$             do ng=2,ngrs
!!$                scalar_g(nsc,ng)%sclt => scalar_g(nsc,1)%sclt
!!$             enddo
!!$          enddo
       end if
       
       allocate (res(nsc)%wetdep(mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(3),"(i8)") mxp
          write(str(4),"(i8)") myp
          call fatal_error(h//" allocate wetdep("//&
               trim(adjustl(str(3)))//","//&
               trim(adjustl(str(4)))//")"//&
               " fails with stat="//trim(adjustl(str(1))))
       end if
       res(nsc)%wetdep=0.0
       
       allocate (res(nsc)%srcsc(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          write(str(2),"(i8)") mzp
          write(str(3),"(i8)") mxp
          write(str(4),"(i8)") myp
          call fatal_error(h//" allocate srcsc("//&
               trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//","//&
               trim(adjustl(str(4)))//")"//&
               " fails with stat="//trim(adjustl(str(1))))
       end if
       res(nsc)%srcsc=0.0

    end do
  end function CreateScalarFields





  function CreateEmptyScalarFields(oneNamelistFile) result(res)
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(ScalarFields), pointer :: res(:)

    integer :: ierr
    integer :: nsc
    integer :: naddsc
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateEmptyScalarFields)**" 

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if

    naddsc = oneNamelistFile%naddsc
    if (naddsc == 0) then
       nullify(res)
       return
    end if

    allocate(res(naddsc), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       write(str(2),"(i8)") naddsc
       call fatal_error(h//" allocate res("//trim(adjustl(str(2)))//&
            " fails with stat="//trim(adjustl(str(1))))
    end if
    
    do nsc=1,naddsc
       nullify (res(nsc)%sclp)
       nullify (res(nsc)%drydep)
       nullify (res(nsc)%sclt)
       nullify (res(nsc)%wetdep)
       nullify (res(nsc)%srcsc)
    end do
  end function CreateEmptyScalarFields



  subroutine DestroyScalarFields(oneScalarFields)
    type(ScalarFields), pointer, intent(inout) :: oneScalarFields(:)

    integer :: ierr
    integer :: nsc
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyScalarFields)**"

    if (.not. associated(oneScalarFields)) then
       return
    end if
    
    do nsc=1,size(oneScalarFields)
       deallocate (oneScalarFields(nsc)%sclp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate sclpfails with stat="//&
               trim(adjustl(str(1))))
       end if
       
       deallocate (oneScalarFields(nsc)%drydep, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate drydep fails with stat="//&
               trim(adjustl(str(1))))
       end if
       
       deallocate (oneScalarFields(nsc)%wetdep, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate wetdep fails with stat="//&
               trim(adjustl(str(1))))
       end if
       
       deallocate (oneScalarFields(nsc)%srcsc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate srcsc fails with stat="//&
               trim(adjustl(str(1))))
       end if

    end do
    
    deallocate(oneScalarFields, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneScalarFields fails with stat="//&
            trim(adjustl(str(1))))
    end if

    nullify(oneScalarFields)
    
  end subroutine DestroyScalarFields





  subroutine DumpScalarFields(oneScalarFields, name)
    type(ScalarFields), pointer, intent(in) :: oneScalarFields(:)
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpScalarFields)**"

    if (associated(oneScalarFields)) then
       call MsgDump(h//" "//name//" associated")
    else
       call MsgDump(h//" "//name//" not associated")
    end if
  end subroutine DumpScalarFields




  


  subroutine InsertScalarFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneScalarFields, oneAveScalarFields, oneNamelistFile, imean)

    ! Fill pointers to arrays into variable tables

    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(ScalarFields), pointer, intent(in) :: oneScalarFields(:)
    type(ScalarFields), pointer, intent(in) :: oneAveScalarFields(:)
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    integer, intent(in) :: imean

    integer :: iSca
    character(len=15) :: sname
    character(len=8) :: str_recycle
    character(len=*), parameter :: h="**(InsertScalarFieldsAtVarTable)**" 

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    end if

    if (.not. associated(oneScalarFields)) then
       return
    end if

    str_recycle = ''
    if (oneNamelistFile%recycle_tracers == 1 .or. &
         oneNamelistFile%ioutput == 5) then
       str_recycle = ':recycle'
    endif

    do iSca=1,size(oneScalarFields)
       if (associated(oneScalarFields(iSca)%sclp)) then
          write(sname,'(a4,i3.3)') 'SCLP',iSca
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneScalarFields(iSca)%sclp, &
               trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle), &
               oneAveScalarFields(iSca)%sclp, imean)
       end if

       if (associated(oneScalarFields(iSca)%drydep)) then
          write(sname,'(a4,i3.3)') 'SCDD',iSca
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneScalarFields(iSca)%drydep, &
               trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1', &
               oneAveScalarFields(iSca)%drydep, imean)
       end if

       if (associated(oneScalarFields(iSca)%wetdep)) then
          write(sname,'(a6,i3.3)') 'wetdep',iSca
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneScalarFields(iSca)%wetdep, &
               trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1', &
               oneAveScalarFields(iSca)%wetdep, imean)
       end if

       if (associated(oneScalarFields(iSca)%srcsc)) then
          write(sname,'(a5,i3.3)') 'scrsc',iSca
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneScalarFields(iSca)%srcsc, &
               trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1', &
               oneAveScalarFields(iSca)%srcsc, imean)
       end if
    end do
  end subroutine InsertScalarFieldsAtVarTable

end module ModScalarFields
