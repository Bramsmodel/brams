!!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_scalar

  use ModNamelistFile, only: namelistFile

  ! Added scalar variables and tendencies

  type scalar_vars
     real, contiguous, pointer :: sclp(:,:,:) => null()
     real, contiguous, pointer :: drydep(:,:) => null()
     real, contiguous, pointer :: sclt(:) => null()
     real, contiguous, pointer :: wetdep(:,:) => null()
     real, contiguous, pointer :: srcsc(:,:,:) => null()
  end type scalar_vars

  ! scal_p allocated by (maxsclr,ngrids)

  type (scalar_vars), pointer :: scalar_g(:,:) => null()
  type (scalar_vars), pointer :: scalarm_g(:,:) => null()

  integer :: recycle_tracers ! from RAMSIN

contains
  !---------------------------------------------------------------

  subroutine alloc_scalar(scal,n1,n2,n3,naddsc)

    implicit none

    integer,intent(in) :: naddsc
    type(scalar_vars), intent(inout) :: scal(:)
    integer,intent(in) :: n1,n2,n3

    integer :: nsc

    ! print *,'Size of scal=' ,size(scal,1)
    ! Allocate arrays based on options (if necessary).
    do nsc=1,naddsc
       !print*,'escalar=',nsc,naddsc,n1,n2,n3
       allocate (scal(nsc)%sclp(n1,n2,n3))
       allocate (scal(nsc)%drydep(n2,n3))
       allocate (scal(nsc)%wetdep(n2,n3))
       allocate (scal(nsc)%srcsc(n1,n2,n3))
       
!--(DMK-LFR NEC-SX6)----------------------------------------------
       scal(nsc)%sclp   = 0.
       scal(nsc)%drydep = 0.
       scal(nsc)%wetdep = 0.
       scal(nsc)%srcsc  = 0.
!--(DMK-LFR NEC-SX6)----------------------------------------------       
       
    enddo

    return
  end subroutine alloc_scalar

  !---------------------------------------------------------------

  subroutine dealloc_scalar(scal)

    implicit none

    type(scalar_vars), pointer, intent(inout) :: scal(:,:)
    
    integer :: nsc
    integer :: ng

    !  Deallocate arrays

    if (associated(scal)) then
       do ng=1,size(scal,2)
          do nsc=1,size(scal,1)
             if (associated(scal(nsc,ng)%sclp))   deallocate (scal(nsc,ng)%sclp)
             if (associated(scal(nsc,ng)%drydep)) deallocate (scal(nsc,ng)%drydep)
             ! For CATT
             if (associated(scal(nsc,ng)%wetdep)) deallocate (scal(nsc,ng)%wetdep)
             if (associated(scal(nsc,ng)%srcsc)) deallocate (scal(nsc,ng)%srcsc)
          end do
       end do
       deallocate(scal)
       nullify(scal)
    end if
  end subroutine dealloc_scalar

  !---------------------------------------------------------------

  subroutine nullify_scalar(scal,naddsc)

    implicit none

    type (scalar_vars) :: scal(naddsc)

    integer :: naddsc
    integer :: nsc

    !  Deallocate arrays

    do nsc=1,naddsc
       if (associated(scal(nsc)%sclp))   nullify (scal(nsc)%sclp)
       if (associated(scal(nsc)%drydep)) nullify (scal(nsc)%drydep)
       ! For CATT
       if (associated(scal(nsc)%wetdep)) nullify (scal(nsc)%wetdep)
       if (associated(scal(nsc)%srcsc)) nullify (scal(nsc)%srcsc)
    enddo

    return
  end subroutine nullify_scalar

  !---------------------------------------------------------------

  subroutine filltab_scalar(oneVarTable, oneVarTableSize, &
       scal, scalm, na, imean)

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    use io_params, only : ioutput         ! INTENT(IN)

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(scalar_vars), pointer, intent(in) :: scal
    type(scalar_vars), pointer, intent(in) :: scalm
    integer, intent(in) :: na
    integer, intent(in) :: imean
    
    character(len=15) :: sname
    character(len=8) :: str_recycle
    character(len=*), parameter :: h="**(filltab_scalar)**" 

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(scal)) then
       call fatal_error(h//" scal not associated")
    end if

    str_recycle = ''
    if (RECYCLE_TRACERS == 1 .or. ioutput == 5) then
       str_recycle = ':recycle'
    endif

    ! Fill pointers to arrays into variable tables

    if (associated(scal%sclp)) then
       write(sname,'(a4,i3.3)') 'SCLP',na
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            scal%sclp, &
            trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle), &
            scalm%sclp, imean)
    end if

    if (associated(scal%drydep)) then
       write(sname,'(a4,i3.3)') 'SCDD',na
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            scal%drydep, &
            trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1', &
            scalm%drydep, imean)
    end if

    if (associated(scal%wetdep)) then
       write(sname,'(a6,i3.3)') 'wetdep',na
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            scal%wetdep, &
            trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1', &
            scalm%wetdep, imean)
    end if

    if (associated(scal%srcsc)) then
       write(sname,'(a5,i3.3)') 'scrsc',na
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            scal%srcsc, &
            trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1', &
            scalm%srcsc, imean)
    end if

  end subroutine filltab_scalar

  subroutine StoreNamelistFileAtMem_scalar(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    recycle_tracers = oneNamelistFile%recycle_tracers
  end subroutine StoreNamelistFileAtMem_scalar
end module mem_scalar
