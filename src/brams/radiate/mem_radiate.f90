!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
!:DOC%BEGIN
!:DOC%TITLE Modulo Mem_Radiate

module mem_radiate

  use ModNamelistFile, only: namelistFile

  type radiate_vars
     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, contiguous :: fthrd(:,:,:)
     real, pointer, contiguous :: cldamnt(:,:,:)
     real, pointer, contiguous :: cluamnt(:,:,:)
     real, pointer, contiguous :: fthrd_sw(:,:,:)
     real, pointer, contiguous :: fthrd_lw(:,:,:)
     real, pointer, contiguous :: cloud_fraction(:,:,:)
     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, contiguous :: rshort(:,:)
     real, pointer, contiguous :: rlong(:,:)
     real, pointer, contiguous :: rlongup(:,:)
     real, pointer, contiguous :: albedt(:,:)
     real, pointer, contiguous :: cosz(:,:)
     real, pointer, contiguous :: rshortdif(:,:)
     real, pointer, contiguous :: sw_up_toa(:,:)
     real, pointer, contiguous :: lw_up_toa(:,:)
  end type radiate_vars

  type (radiate_vars), pointer :: radiate_g(:), radiatem_g(:)

  integer :: lonrad ! from RAMSIN
  integer :: ilwrtyp ! from RAMSIN
  integer :: iswrtyp ! from RAMSIN
  real    :: radfrq ! from RAMSIN
  real    :: radtun ! from RAMSIN
  integer :: ncall_i !Indica primeira chamada
  real    :: prsnz,prsnzp !Calculadas na primeira chamada

contains

  subroutine alloc_radiate(radiate,n1,n2,n3,ng)

    implicit none
    type (radiate_vars) :: radiate
    integer, intent(in) :: n1,n2,n3,ng

    ! Allocate arrays based on options (if necessary)

    if( ilwrtyp+iswrtyp > 0)  then
       allocate (radiate%fthrd(n1,n2,n3))
       allocate (radiate%cloud_fraction(n1,n2,n3))
       allocate (radiate%rshort(n2,n3))
       allocate (radiate%rlong(n2,n3))
       allocate (radiate%rlongup(n2,n3))
       allocate (radiate%albedt(n2,n3))
       allocate (radiate%cosz(n2,n3))


       !-20/10/2015: srf - not being used, actually
       !
       !--(DMK-CCATT-INI)-----------------------------------------------------------
       !NER
       !    allocate (radiate%sw_up_toa(n2,n3))
       !    allocate (radiate%lw_up_toa(n2,n3))
       !    allocate (radiate%rshortdif(n2,n3))
       !    allocate (radiate%cldamnt(n1,n2,n3))
       !    allocate (radiate%cluamnt(n1,n2,n3))  
       !    allocate (radiate%fthrd_sw(n1,n2,n3)) 
       !    allocate (radiate%fthrd_lw(n1,n2,n3)) 
       !--(DMK-CCATT-FIM)-----------------------------------------------------------

       !--(DMK-LFR NEC-SX6)----------------------------------------------
       radiate%fthrd   = 0.
       radiate%rshort  = 0.
       radiate%rlong   = 0.
       radiate%rlongup = 0.
       radiate%albedt  = 0.
       radiate%cosz    = 0.
       !-20/10/2015: srf - not being used, actually
       !--(DMK-CCATT-INI)-----------------------------------------------------------
       !    !NER
       !    radiate%cldamnt  = 0.
       !    radiate%cluamnt  = 0.
       !    radiate%fthrd_sw  = 0.
       !    radiate%fthrd_sw  = 0.
       !    radiate%sw_up_toa  = 0.
       !    radiate%lw_up_toa  = 0.
       !    radiate%rshortdif  = 0. 
       !--(DMK-CCATT-FIM)-----------------------------------------------------------    
       !--(DMK-LFR NEC-SX6)----------------------------------------------

    endif

    return
  end subroutine alloc_radiate


  subroutine nullify_radiate(radiate)

    implicit none
    type (radiate_vars) :: radiate


    if (associated(radiate%fthrd))    nullify (radiate%fthrd)
    if (associated(radiate%cloud_fraction)) nullify (radiate%cloud_fraction)
    if (associated(radiate%rshort))   nullify (radiate%rshort)
    if (associated(radiate%rlong))    nullify (radiate%rlong)
    if (associated(radiate%rlongup))  nullify (radiate%rlongup)
    if (associated(radiate%albedt))   nullify (radiate%albedt)
    if (associated(radiate%cosz))     nullify (radiate%cosz)

    !--(DMK-CCATT-INI)-----------------------------------------------------------
    !NER
    if (associated(radiate%cldamnt))  nullify (radiate%cldamnt)
    if (associated(radiate%cluamnt))  nullify (radiate%cluamnt)
    if (associated(radiate%rshortdif))nullify (radiate%rshortdif)
    if (associated(radiate%sw_up_toa))nullify (radiate%sw_up_toa)
    if (associated(radiate%lw_up_toa))nullify (radiate%lw_up_toa)
    if (associated(radiate%fthrd_sw)) nullify (radiate%fthrd_sw) 
    if (associated(radiate%fthrd_lw)) nullify (radiate%fthrd_lw)
    !--(DMK-CCATT-FIM)-----------------------------------------------------------

    return
  end subroutine nullify_radiate

  subroutine dealloc_radiate(radiate)

    implicit none
    type (radiate_vars) :: radiate


    if (associated(radiate%fthrd))    deallocate (radiate%fthrd)
    if (associated(radiate%cloud_fraction)) deallocate (radiate%cloud_fraction)
    if (associated(radiate%rshort))   deallocate (radiate%rshort)
    if (associated(radiate%rlong))    deallocate (radiate%rlong)
    if (associated(radiate%rlongup))  deallocate (radiate%rlongup)
    if (associated(radiate%albedt))   deallocate (radiate%albedt)
    if (associated(radiate%cosz))     deallocate (radiate%cosz)

    !--(DMK-CCATT-INI)-----------------------------------------------------------
    !NER
    if (associated(radiate%cldamnt))  deallocate (radiate%cldamnt)
    if (associated(radiate%cluamnt))  deallocate (radiate%cluamnt) 
    if (associated(radiate%rshortdif))deallocate (radiate%rshortdif)
    if (associated(radiate%sw_up_toa))deallocate (radiate%sw_up_toa)
    if (associated(radiate%lw_up_toa))deallocate (radiate%lw_up_toa)
    if (associated(radiate%fthrd_sw)) deallocate (radiate%fthrd_sw)
    if (associated(radiate%fthrd_lw)) deallocate (radiate%fthrd_lw)
    !--(DMK-CCATT-FIM)-----------------------------------------------------------

    return
  end subroutine dealloc_radiate


  subroutine filltab_radiate(oneVarTable, oneVarTableSize, &
       radiate, radiatem)

    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable
    
    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(radiate_vars), pointer, intent(in) :: radiate
    type(radiate_vars), pointer, intent(in) :: radiatem

    logical :: assAve

    assAve=associated(radiatem)

    ! Fill pointers to arrays into variable tables

    if (associated(radiate%cloud_fraction))  then
       if (assAve) then
          if (associated(radiatem%cloud_fraction)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%cloud_fraction, &
                  'CLOUD_FRACTION :3:anal:mpti:mpt3', &
                  radiatem%cloud_fraction)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%cloud_fraction, &
                  'CLOUD_FRACTION :3:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%cloud_fraction, &
               'CLOUD_FRACTION :3:anal:mpti:mpt3')
       end if
    end if
         
    if (associated(radiate%fthrd))  then
       if (assAve) then
          if (associated(radiatem%fthrd)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%fthrd, &
                  'FTHRD :3:hist:anal:mpti:mpt3', &
                  radiatem%fthrd)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%fthrd, &
                  'FTHRD :3:hist:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%fthrd, &
               'FTHRD :3:hist:anal:mpti:mpt3')
       end if
    end if
         
    if (associated(radiate%cldamnt))  then
       if (assAve) then
          if (associated(radiatem%cldamnt)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%cldamnt, &
                  'CLDAMNT :3:hist:anal:mpti:mpt3', &
                  radiatem%cldamnt)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%cldamnt, &
                  'CLDAMNT :3:hist:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%cldamnt, &
               'CLDAMNT :3:hist:anal:mpti:mpt3')
       end if
    end if
         
    if (associated(radiate%fthrd_sw))  then
       if (assAve) then
          if (associated(radiatem%fthrd_sw)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%fthrd_sw, &
                  'FTHRD_SW :3:hist:anal:mpti:mpt3', &
                  radiatem%fthrd_sw)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%fthrd_sw, &
                  'FTHRD_SW :3:hist:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%fthrd_sw, &
               'FTHRD_SW :3:hist:anal:mpti:mpt3')
       end if
    end if
         
    if (associated(radiate%fthrd_lw))  then
       if (assAve) then
          if (associated(radiatem%fthrd_lw)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%fthrd_lw, &
                  'FTHRD_LW :3:hist:anal:mpti:mpt3', &
                  radiatem%fthrd_lw)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%fthrd_lw, &
                  'FTHRD_LW :3:hist:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%fthrd_lw, &
               'FTHRD_LW :3:hist:anal:mpti:mpt3')
       end if
    end if
         
    if (associated(radiate%cluamnt))  then
       if (assAve) then
          if (associated(radiatem%cluamnt)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%cluamnt, &
                  'CLUAMNT :3:hist:anal:mpti:mpt3', &
                  radiatem%cluamnt)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%cluamnt, &
                  'CLUAMNT :3:hist:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%cluamnt, &
               'CLUAMNT :3:hist:anal:mpti:mpt3')
       end if
    end if
         
    if (associated(radiate%rshort))  then
       if (assAve) then
          if (associated(radiatem%rshort)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%rshort, &
                  'RSHORT :2:hist:anal:mpti:mpt3', &
                  radiatem%rshort)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%rshort, &
                  'RSHORT :2:hist:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%rshort, &
               'RSHORT :2:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(radiate%rlong))  then
       if (assAve) then
          if (associated(radiatem%rlong)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%rlong, &
                  'RLONG :2:hist:anal:mpti:mpt3', &
                  radiatem%rlong)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%rlong, &
                  'RLONG :2:hist:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%rlong, &
               'RLONG :2:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(radiate%rlongup))  then
       if (assAve) then
          if (associated(radiatem%rlongup)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%rlongup, &
                  'RLONGUP :2:hist:anal:mpti:mpt3', &
                  radiatem%rlongup)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%rlongup, &
                  'RLONGUP :2:hist:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%rlongup, &
               'RLONGUP :2:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(radiate%albedt))  then
       if (assAve) then
          if (associated(radiatem%albedt)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%albedt, &
                  'ALBEDT :2:hist:anal:mpti:mpt3', &
                  radiatem%albedt)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%albedt, &
                  'ALBEDT :2:hist:anal:mpti:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%albedt, &
               'ALBEDT :2:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(radiate%cosz))  then
       if (assAve) then
          if (associated(radiatem%cosz)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%cosz, &
                  'COSZ :2:hist:anal:mpt3', &
                  radiatem%cosz)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%cosz, &
                  'COSZ :2:hist:anal:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%cosz, &
               'COSZ :2:hist:anal:mpt3')
       end if
    end if
         
    if (associated(radiate%rshortdif))  then
       if (assAve) then
          if (associated(radiatem%rshortdif)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%rshortdif, &
                  'RSHORTDIF :2:hist:anal:mpt3', &
                  radiatem%rshortdif)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%rshortdif, &
                  'RSHORTDIF :2:hist:anal:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%rshortdif, &
               'RSHORTDIF :2:hist:anal:mpt3')
       end if
    end if
         
    if (associated(radiate%sw_up_toa))  then
       if (assAve) then
          if (associated(radiatem%sw_up_toa)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%sw_up_toa, &
                  'SW_UP_TOA :2:hist:anal:mpt3', &
                  radiatem%sw_up_toa)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%sw_up_toa, &
                  'SW_UP_TOA :2:hist:anal:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%sw_up_toa, &
               'SW_UP_TOA :2:hist:anal:mpt3')
       end if
    end if
         
    if (associated(radiate%lw_up_toa))  then
       if (assAve) then
          if (associated(radiatem%lw_up_toa)) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%lw_up_toa, &
                  'LW_UP_TOA :2:hist:anal:mpt3', &
                  radiatem%lw_up_toa)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  radiate%lw_up_toa, &
                  'LW_UP_TOA :2:hist:anal:mpt3')
          end if
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               radiate%lw_up_toa, &
               'LW_UP_TOA :2:hist:anal:mpt3')
       end if
    end if

  end subroutine filltab_radiate


  subroutine StoreNamelistFileAtMem_radiate(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    ilwrtyp = oneNamelistFile%ilwrtyp
    iswrtyp = oneNamelistFile%iswrtyp
    lonrad = oneNamelistFile%lonrad
    radfrq = oneNamelistFile%radfrq
    radtun = oneNamelistFile%radtun
  end subroutine StoreNamelistFileAtMem_radiate
end module mem_radiate
