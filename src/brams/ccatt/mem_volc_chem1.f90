!###########################################################################
! CCATT- B - Regional Atmospheric Modeling System - RAMS
!###########################################################################

module mem_volc_chem1

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  use ModNamelistFile, only: namelistFile

  include "constants.h"
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  type volc_mean_vars   
     real, pointer, contiguous :: plum_heigth(:,:) => null()
     real, pointer, contiguous :: vent_elev(:,:) => null()
     real, pointer, contiguous :: duration(:,:) => null()
     real, pointer, contiguous :: begin_time(:,:) => null()
  end type volc_mean_vars  

  type(volc_mean_vars), allocatable, target :: volc_mean_g(:)
  type(volc_mean_vars), allocatable, target :: volc_meanm_g(:)
   
  integer:: volcanoes 

contains
  !---------------------------------------------------------------

  subroutine alloc_volc_chem1(volc_mean,n1,n2,n3)

    !use chem1_list, only : spc_alloc,spc_name,src,on
    implicit none

    integer,intent(in) :: n1,n2,n3
    type (volc_mean_vars) :: volc_mean
        
    allocate (volc_mean%plum_heigth(n2,n3));  volc_mean%plum_heigth(:,:)=0.
    allocate (volc_mean%vent_elev   (n2,n3)); volc_mean%vent_elev  (:,:)=0.
    allocate (volc_mean%duration   (n2,n3)); volc_mean%duration   (:,:)=0.
    !allocate (volc_mean%begin_time (n2,n3)); volc_mean%begin_time (:,:)=0.

  end subroutine alloc_volc_chem1

  !---------------------------------------------------------------
  !subroutine dealloc_volc_chem1(volc_mean)
  !
  ! implicit none
  !
  ! 
  !end subroutine dealloc_volc_chem1
  !---------------------------------------------------------------

  subroutine nullify_volc_chem1(volc_mean)

    implicit none
    type (volc_mean_vars) :: volc_mean

       if (associated(volc_mean%plum_heigth)) nullify(volc_mean%plum_heigth)
       if (associated(volc_mean%vent_elev  )) nullify(volc_mean%vent_elev   )
       if (associated(volc_mean%duration )) nullify(volc_mean%duration )
       !if (associated(volc_mean%begin_time )) nullify(volc_mean%begin_time )
     
  end subroutine nullify_volc_chem1

  !---------------------------------------------------------------

  subroutine filltab_volc_chem1(oneVarTable, oneVarTableSize, &
       volc_mean, volc_meanm)

    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(volc_mean_vars), pointer, intent(in) :: volc_mean
    type(volc_mean_vars), pointer, intent(in) :: volc_meanm

    logical :: assThis
    character(len=*), parameter :: h="**(filltab_volc_chem1)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(volc_mean)) then
       call fatal_error(h//" volc_mean not associated")
    end if
    
    ! Fill pointers to arrays into variable tables
    ! 2d var

    if (.not. associated(volc_meanm)) then
       assThis=.false.
    else
       assThis=associated(volc_meanm%plum_heigth)
    end if
    if (assThis) then
       call InsertAtVarTable(oneVarTable, oneVarTableSize, &
            volc_mean%plum_heigth, &
            'plum_heigth_volc :2:hist:anal:mpti:mpt3:mpt1', &
            volc_meanm%plum_heigth)
    else
       call InsertAtVarTable(oneVarTable, oneVarTableSize, &
            volc_mean%plum_heigth, &
            'plum_heigth_volc :2:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (.not. associated(volc_meanm)) then
       assThis=.false.
    else
       assThis=associated(volc_meanm%vent_elev)
    end if
    if (assThis) then
       call InsertAtVarTable(oneVarTable, oneVarTableSize, &
            volc_mean%vent_elev, &
            'vent_elev_volc :2:hist:anal:mpti:mpt3:mpt1', &
            volc_meanm%vent_elev)
    else
       call InsertAtVarTable(oneVarTable, oneVarTableSize, &
            volc_mean%vent_elev, &
            'vent_elev_volc :2:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (.not. associated(volc_meanm)) then
       assThis=.false.
    else
       assThis=associated(volc_meanm%duration)
    end if
    if (assThis) then
       call InsertAtVarTable(oneVarTable, oneVarTableSize, &
            volc_mean%duration, &
            'duration_volc :2:hist:anal:mpti:mpt3:mpt1', &
            volc_meanm%duration)
    else
       call InsertAtVarTable(oneVarTable, oneVarTableSize, &
            volc_mean%duration, &
            'duration_volc :2:hist:anal:mpti:mpt3:mpt1')
    end if
  end subroutine filltab_volc_chem1
  !---------------------------------------------------------------


!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  subroutine StoreNamelistFileAtMem_volcChem1(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    volcanoes = oneNamelistFile%volcanoes
  end subroutine StoreNamelistFileAtMem_volcChem1
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


end module mem_volc_chem1
