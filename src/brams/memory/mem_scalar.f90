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
     real, contiguous, pointer :: sclp(:,:,:)
     real, contiguous, pointer :: drydep(:,:)
     real, contiguous, pointer :: sclt(:)
     real, contiguous, pointer :: wetdep(:,:)
     real, contiguous, pointer :: srcsc(:,:,:)
  end type scalar_vars

  ! scal_p allocated by (maxsclr,ngrids)

  type (scalar_vars), pointer :: scalar_g(:,:)
  type (scalar_vars), pointer :: scalarm_g(:,:)

  integer :: recycle_tracers ! from RAMSIN

contains
  !---------------------------------------------------------------

  subroutine alloc_scalar(scal,n1,n2,n3,naddsc)

    implicit none

    integer,intent(in) :: naddsc
    type (scalar_vars),dimension(naddsc) :: scal
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

  subroutine dealloc_scalar(scal,naddsc)

    implicit none

    type (scalar_vars) :: scal(naddsc)
    integer :: naddsc
    integer :: nsc

    !  Deallocate arrays

    do nsc=1,naddsc
       if (associated(scal(nsc)%sclp))   deallocate (scal(nsc)%sclp)
       if (associated(scal(nsc)%drydep)) deallocate (scal(nsc)%drydep)
       ! For CATT
       if (associated(scal(nsc)%wetdep)) deallocate (scal(nsc)%wetdep)
       if (associated(scal(nsc)%srcsc)) deallocate (scal(nsc)%srcsc)
    enddo

    return
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
       scal, scalm, na)

    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable
    
    use io_params, only : ioutput         ! INTENT(IN)

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(scalar_vars), pointer, intent(in) :: scal
    type(scalar_vars), pointer, intent(in) :: scalm
    integer, intent(in) :: na

    logical :: assAve
    logical :: assThis
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

    assAve=associated(scalm)
    
    ! Fill pointers to arrays into variable tables

    if (associated(scal%sclp)) then
       write(sname,'(a4,i3.3)') 'SCLP',na
       if (assAve) then
          assThis=associated(scalm%sclp)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               scal%sclp, &
               trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle), &
               scalm%sclp)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               scal%sclp, &
               trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))
       end if
    end if
    
    if (associated(scal%drydep)) then
       write(sname,'(a4,i3.3)') 'SCDD',na
       if (assAve) then
          assThis=associated(scalm%drydep)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               scal%drydep, &
               trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1', &
               scalm%drydep)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               scal%drydep, &
               trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
    
    if (associated(scal%wetdep)) then
       write(sname,'(a6,i3.3)') 'wetdep',na
       if (assAve) then
          assThis=associated(scalm%wetdep)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               scal%wetdep, &
               trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1', &
               scalm%wetdep)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               scal%wetdep, &
               trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
    
    if (associated(scal%srcsc)) then
       write(sname,'(a5,i3.3)') 'scrsc',na
       if (assAve) then
          assThis=associated(scalm%srcsc)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               scal%srcsc, &
               trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1', &
               scalm%srcsc)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               scal%srcsc, &
               trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1')
       endif
    end if

  end subroutine filltab_scalar

  subroutine StoreNamelistFileAtMem_scalar(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    recycle_tracers = oneNamelistFile%recycle_tracers
  end subroutine StoreNamelistFileAtMem_scalar
end module mem_scalar
