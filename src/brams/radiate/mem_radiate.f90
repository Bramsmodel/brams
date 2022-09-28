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

  use ModRadiateFields, only: &
       RadiateFields
  
  type radiate_vars
     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, contiguous :: fthrd(:,:,:) => null()
     real, pointer, contiguous :: cldamnt(:,:,:) => null()
     real, pointer, contiguous :: cluamnt(:,:,:) => null()
     real, pointer, contiguous :: fthrd_sw(:,:,:) => null()
     real, pointer, contiguous :: fthrd_lw(:,:,:) => null()
     real, pointer, contiguous :: cloud_fraction(:,:,:) => null()
     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, contiguous :: rshort(:,:) => null()
     real, pointer, contiguous :: rlong(:,:) => null()
     real, pointer, contiguous :: rlongup(:,:) => null()
     real, pointer, contiguous :: albedt(:,:) => null()
     real, pointer, contiguous :: cosz(:,:) => null()
     real, pointer, contiguous :: rshortdif(:,:) => null()
     real, pointer, contiguous :: sw_up_toa(:,:) => null()
     real, pointer, contiguous :: lw_up_toa(:,:) => null()
  end type radiate_vars

  type (radiate_vars), pointer, contiguous :: radiate_g(:)
  type (radiate_vars), pointer, contiguous :: radiatem_g(:)

  integer :: lonrad ! from RAMSIN
  integer :: ilwrtyp ! from RAMSIN
  integer :: iswrtyp ! from RAMSIN
  real    :: radfrq ! from RAMSIN
  real    :: radtun ! from RAMSIN
  integer :: ncall_i !Indica primeira chamada

  character(len=128), private :: lastTo=""
  character(len=128), private :: lastFrom=""
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
       radiate, radiatem, ng, imean)

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(radiate_vars), pointer, intent(in) :: radiate(:)
    type(radiate_vars), pointer, intent(in) :: radiatem(:)
    integer, intent(in) :: ng
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(filltab_radiate)**"


    ! Fill pointers to arrays into variable tables

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(radiate)) then
       call fatal_error(h//" radiate not associated")
    else if (.not. associated(radiatem)) then
       call fatal_error(h//" radiatem not associated")
    end if

    if (associated(radiate(ng)%cloud_fraction))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%cloud_fraction, &
            'CLOUD_FRACTION :3:anal:mpti:mpt3', &
            radiatem(ng)%cloud_fraction, imean)
    end if

    if (associated(radiate(ng)%fthrd))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%fthrd, &
            'FTHRD :3:hist:anal:mpti:mpt3', &
            radiatem(ng)%fthrd, imean)
    end if

    if (associated(radiate(ng)%cldamnt))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%cldamnt, &
            'CLDAMNT :3:hist:anal:mpti:mpt3', &
            radiatem(ng)%cldamnt, imean)
    end if

    if (associated(radiate(ng)%fthrd_sw))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%fthrd_sw, &
            'FTHRD_SW :3:hist:anal:mpti:mpt3', &
            radiatem(ng)%fthrd_sw, imean)
    end if

    if (associated(radiate(ng)%fthrd_lw))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%fthrd_lw, &
            'FTHRD_LW :3:hist:anal:mpti:mpt3', &
            radiatem(ng)%fthrd_lw, imean)
    end if

    if (associated(radiate(ng)%cluamnt))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%cluamnt, &
            'CLUAMNT :3:hist:anal:mpti:mpt3', &
            radiatem(ng)%cluamnt, imean)
    end if

    if (associated(radiate(ng)%rshort))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%rshort, &
            'RSHORT :2:hist:anal:mpti:mpt3', &
            radiatem(ng)%rshort, imean)
    end if

    if (associated(radiate(ng)%rlong))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%rlong, &
            'RLONG :2:hist:anal:mpti:mpt3', &
            radiatem(ng)%rlong, imean)
    end if

    if (associated(radiate(ng)%rlongup))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%rlongup, &
            'RLONGUP :2:hist:anal:mpti:mpt3', &
            radiatem(ng)%rlongup, imean)
    end if

    if (associated(radiate(ng)%albedt))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%albedt, &
            'ALBEDT :2:hist:anal:mpti:mpt3', &
            radiatem(ng)%albedt, imean)
    end if

    if (associated(radiate(ng)%cosz))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%cosz, &
            'COSZ :2:hist:anal:mpt3', &
            radiatem(ng)%cosz, imean)
    end if

    if (associated(radiate(ng)%rshortdif))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%rshortdif, &
            'RSHORTDIF :2:hist:anal:mpt3', &
            radiatem(ng)%rshortdif, imean)
    end if

    if (associated(radiate(ng)%sw_up_toa))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%sw_up_toa, &
            'SW_UP_TOA :2:hist:anal:mpt3', &
            radiatem(ng)%sw_up_toa, imean)
    end if

    if (associated(radiate(ng)%lw_up_toa))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            radiate(ng)%lw_up_toa, &
            'LW_UP_TOA :2:hist:anal:mpt3', &
            radiatem(ng)%lw_up_toa, imean)
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




  subroutine DeepCopyToRadiateFields(oneRadiateFields, name)
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DeepCopyToRadiateFields)**"

    if (lastTo /= "") then
       call fatal_error(h//" wrong order; previously invoked by "//&
            trim(adjustl(lastTo))//" and now invoked by "//&
            trim(adjustl(name)))
    end if

    lastTo=name
    lastFrom=""
    
    oneRadiateFields%fthrd = radiate_g(1)%fthrd
    oneRadiateFields%cloud_fraction = radiate_g(1)%cloud_fraction
    oneRadiateFields%rshort = radiate_g(1)%rshort
    oneRadiateFields%rlong = radiate_g(1)%rlong
    oneRadiateFields%rlongup = radiate_g(1)%rlongup
    oneRadiateFields%albedt = radiate_g(1)%albedt
    oneRadiateFields%cosz = radiate_g(1)%cosz
  end subroutine DeepCopyToRadiateFields




  subroutine DeepCopyFromRadiateFields(oneRadiateFields, name)
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DeepCopyFromRadiateFields)**"

    if (lastFrom /= "") then
       call fatal_error(h//" wrong order; previously invoked by "//&
            trim(adjustl(lastFrom))//" and now invoked by "//&
            trim(adjustl(name)))
    end if

    lastFrom=name
    lastTo=""
    
    radiate_g(1)%fthrd = oneRadiateFields%fthrd
    radiate_g(1)%cloud_fraction = oneRadiateFields%cloud_fraction
    radiate_g(1)%rshort = oneRadiateFields%rshort
    radiate_g(1)%rlong = oneRadiateFields%rlong
    radiate_g(1)%rlongup = oneRadiateFields%rlongup
    radiate_g(1)%albedt = oneRadiateFields%albedt
    radiate_g(1)%cosz = oneRadiateFields%cosz
  end subroutine DeepCopyFromRadiateFields  
end module mem_radiate
