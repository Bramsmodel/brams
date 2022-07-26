!############################# Change Log ##################################
! 5.0.0
!
! Demerval S. Moreira - 27/Ago/2012
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModJulesFields

  use iso_fortran_env, only: &
       int64
  
  use ModNamelistFile, only: &
       NamelistFile

  use ModParallelEnvironment, only: &
       MsgDump

  use ModNodeDimensions, only: &
       NodeDimensions

  use ModVarTables, only: &
       InsertVTab
  
  implicit none

  private

  public :: JulesFields
  public :: CreateJulesFields
  public :: DestroyJulesFields
  public :: DumpJulesFields
  public :: InsertJulesFieldsAtVarTable
  
  type JulesFields
     real, pointer, contiguous :: gpp(:,:) => null()      
     real, pointer, contiguous :: resp_s(:,:) => null()   
     real, pointer, contiguous :: resp_p(:,:) => null()   
     real, pointer, contiguous :: npp(:,:) => null()      
     real, pointer, contiguous :: u10mj(:,:) => null()
     real, pointer, contiguous :: v10mj(:,:) => null()
     real, pointer, contiguous :: u10mj1hr(:,:,:) => null()
     real, pointer, contiguous :: v10mj1hr(:,:,:) => null()
     real, pointer, contiguous :: fracj(:,:,:) => null()
     real, pointer, contiguous :: t2mj(:,:) => null()
     real, pointer, contiguous :: t2mj_max(:,:) => null()
     real, pointer, contiguous :: t2mj_min(:,:) => null()
     real, pointer, contiguous :: rv2mj(:,:) => null()
     real, pointer, contiguous :: csj(:,:) => null()
     real, pointer, contiguous :: anthrop_heatj(:,:,:) => null()
     real, pointer, contiguous :: radnet_tilej(:,:,:) => null()
     real, pointer, contiguous :: ftl_tilej(:,:,:) => null()
     real, pointer, contiguous :: le_tilej(:,:,:) => null()
     real, pointer, contiguous :: htf_tilej(:,:,:) => null()
     real, pointer, contiguous :: snowdepthj(:,:,:) => null()
     real, pointer, contiguous :: ht_fluxj(:,:) => null()
     real, pointer, contiguous :: temp_surfj(:,:) => null()
  end type JulesFields

contains




  
  function CreateJulesFields(oneNodeDimensions, oneNamelistFile) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(JulesFields), pointer :: res

    integer :: mzp
    integer :: mxp
    integer :: myp
    integer :: nzg
    integer :: nzs
    integer :: npatch
    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateJulesFields)**"

    if (.not. associated(oneNodeDimensions)) then
       call fatal_error(h//" oneNodeDimensions not associated")
    else if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if
    
    mzp=oneNodeDimensions%mzp
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp
    nzg = oneNamelistFile%nzg
    nzs = oneNamelistFile%nzs
    npatch = oneNamelistFile%npatch

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if
    
    ! Allocate arrays based on options (if necessary)

    allocate (res%gpp(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate gpp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%resp_s(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate resp_s fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%resp_p(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate resp_p fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%npp(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate npp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%u10mj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate u10mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%v10mj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate v10mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%u10mj1hr(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate u10mj1hr fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%v10mj1hr(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate v10mj1hr fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%fracj(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate fracj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%t2mj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate t2mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%t2mj_max(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate t2mj_max fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%t2mj_min(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate t2mj_min fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%rv2mj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate rv2mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%csj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate csj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%anthrop_heatj(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate anthrop_heatj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%radnet_tilej(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate radnet_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%ftl_tilej(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate ftl_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%le_tilej(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate le_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%htf_tilej(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate htf_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%snowdepthj(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate snowdepthj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%ht_fluxj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate ht_fluxj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate (res%temp_surfj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate temp_surfj fails with stat="//&
            trim(adjustl(str(1))))
    end if

  end function CreateJulesFields


  

  subroutine DestroyJulesFields(oneJulesFields)
    type(JulesFields), pointer, intent(inout) :: oneJulesFields

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyJulesFields)**"

    if (.not. associated(oneJulesFields)) then
       return
    end if

    deallocate(oneJulesFields%gpp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate gpp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%resp_s, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate resp_s fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%resp_p, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate resp_p fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%npp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate npp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%u10mj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate u10mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%v10mj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate v10mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%u10mj1hr, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate u10mj1hr fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%v10mj1hr, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate v10mj1hr fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%fracj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate fracj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%t2mj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate t2mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%t2mj_max, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate t2mj_max fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%t2mj_min, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate t2mj_min fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%rv2mj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate rv2mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%csj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate csj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%anthrop_heatj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate anthrop_heatj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%radnet_tilej, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate radnet_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%ftl_tilej, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate ftl_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%le_tilej, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate le_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%htf_tilej, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate htf_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%snowdepthj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate snowdepthj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%ht_fluxj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate ht_fluxj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields%temp_surfj, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate temp_surfj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    deallocate(oneJulesFields, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneJulesFields fails with stat="//&
            trim(adjustl(str(1))))
    end if
    nullify(oneJulesFields)
  end subroutine DestroyJulesFields




  

  subroutine DumpJulesFields(oneJulesFields, name)
    type(JulesFields), pointer, intent(in) :: oneJulesFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpJulesFields)**"

    if (associated(oneJulesFields)) then
       call MsgDump(h//" "//name//" associated")
    else
       call MsgDump(h//" "//name//" not associated")
    end if
  end subroutine DumpJulesFields

  


  subroutine InsertJulesFieldsAtVarTable(oneJulesFields, oneAveJulesFields, &
       oneNamelistFile, gridId)
    type(JulesFields), pointer, intent(in) :: oneJulesFields
    type(JulesFields), pointer, intent(in) :: oneAveJulesFields
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    integer, intent(in) :: gridId

    integer :: imean
    integer(kind=int64) :: nptsXY
    integer(kind=int64) :: nptsXYZ
    integer(kind=int64) :: nptsXYP
    character(len=8) :: str_recycle
    character(len=*), parameter :: h="**(InsertJulesFieldsAtVarTable)**" 

    if (.not. associated(oneJulesFields)) then
       call fatal_error(h//" oneJulesFields not associated")
    else if (.not. associated(oneAveJulesFields)) then
       call fatal_error(h//" oneAveJulesFields not associated")
    else if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if

    ! Should average fields be stored at variable tables?

    if (oneNamelistFile%avgtim == 0) then
       imean=0 ! do not store
    else
       imean=1 ! store
    end if

    if (oneNamelistFile%ipastin == 1) then
       str_recycle = ':recycle'
    else
       str_recycle = ''
    endif

    ! Fill pointers to arrays into variable tables

    nptsXY=int(size(oneJulesFields%gpp,1),int64)* &
         int(size(oneJulesFields%gpp,2),int64)

    call InsertVTab (oneJulesFields%gpp,oneAveJulesFields%gpp  &           
         ,gridId, nptsXY, imean, 'GPP :2:anal:mpti:mpt3')  
	 
    call InsertVTab (oneJulesFields%resp_s,oneAveJulesFields%resp_s  &    
         ,gridId, nptsXY, imean, 'RESP_S :2:anal:mpti:mpt3') 
	 
    call InsertVTab (oneJulesFields%resp_p,oneAveJulesFields%resp_p  &     
         ,gridId, nptsXY, imean, 'RESP_P :2:anal:mpti:mpt3') 
	 
    call InsertVTab (oneJulesFields%npp,oneAveJulesFields%npp  &           
         ,gridId, nptsXY, imean, 'NPP :2:anal:mpti:mpt3')    

    call InsertVTab (oneJulesFields%u10mj,oneAveJulesFields%u10mj  &
         ,gridId, nptsXY, imean, 'U10MJ :2:anal:mpti:mpt3')

    call InsertVTab (oneJulesFields%v10mj,oneAveJulesFields%v10mj  &
         ,gridId, nptsXY, imean, 'V10MJ :2:anal:mpti:mpt3')

    call InsertVTab (oneJulesFields%t2mj,oneAveJulesFields%t2mj  &
         ,gridId, nptsXY, imean, 'T2MJ :2:anal:mpti:mpt3')

    call InsertVTab (oneJulesFields%t2mj_max,oneAveJulesFields%t2mj_max  &
         ,gridId, nptsXY, imean, 'T2MJ_MAX :2:anal:mpti:mpt3')

    call InsertVTab (oneJulesFields%t2mj_min,oneAveJulesFields%t2mj_min  &
         ,gridId, nptsXY, imean, 'T2MJ_MIN :2:anal:mpti:mpt3')

    call InsertVTab (oneJulesFields%rv2mj,oneAveJulesFields%rv2mj  &
         ,gridId, nptsXY, imean, 'RV2MJ :2:anal:mpti:mpt3')

    call InsertVTab (oneJulesFields%csj,oneAveJulesFields%csj  &
         ,gridId, nptsXY, imean, 'CSJ :2:anal:mpti:mpt3')

    call InsertVTab (oneJulesFields%ht_fluxj,oneAveJulesFields%ht_fluxj  &
         ,gridId, nptsXY, imean, 'ht_fluxj :2:anal:mpti:mpt3')

    call InsertVTab (oneJulesFields%temp_surfj,oneAveJulesFields%temp_surfj  &
         ,gridId, nptsXY, imean, 'temp_surfj :2:anal:mpti:mpt3')

    nptsXYZ=int(size(oneJulesFields%u10mj1hr,1),int64)* &
         int(size(oneJulesFields%u10mj1hr,2),int64)* &
         int(size(oneJulesFields%u10mj1hr,3),int64)

    call InsertVTab (oneJulesFields%u10mj1hr,oneAveJulesFields%u10mj1hr  &
         ,gridId, nptsXYZ, imean, 'U10MJ1hr :3:hist:anal:mpti:mpt3:mpt2')

    call InsertVTab (oneJulesFields%v10mj1hr,oneAveJulesFields%v10mj1hr  &
         ,gridId, nptsXYZ, imean, 'V10MJ1hr :3:hist:anal:mpti:mpt3:mpt2')

    call InsertVTab (oneJulesFields%fracj,oneAveJulesFields%fracj  &
         ,gridId, nptsXYZ, imean, 'fracj :3:hist:anal:mpti:mpt3:mpt2')

    nptsXYP=int(size(oneJulesFields%anthrop_heatj,1),int64)* &
         int(size(oneJulesFields%anthrop_heatj,2),int64)* &
         int(size(oneJulesFields%anthrop_heatj,3),int64)

    call InsertVTab (oneJulesFields%anthrop_heatj,oneAveJulesFields%anthrop_heatj  &
         ,gridId, nptsXYP, imean,  &
         'anthrop_heatj :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (oneJulesFields%radnet_tilej,oneAveJulesFields%radnet_tilej  &
         ,gridId, nptsXYP, imean,  &
         'radnet_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (oneJulesFields%ftl_tilej,oneAveJulesFields%ftl_tilej  &
         ,gridId, nptsXYP, imean,  &
         'ftl_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (oneJulesFields%le_tilej,oneAveJulesFields%le_tilej  &
         ,gridId, nptsXYP, imean,  &
         'le_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (oneJulesFields%htf_tilej,oneAveJulesFields%htf_tilej  &
         ,gridId, nptsXYP, imean,  &
         'htf_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (oneJulesFields%snowdepthj,oneAveJulesFields%snowdepthj  &
         ,gridId, nptsXYP, imean,  &
         'snowdepthj :6:hist:anal:mpti:mpt3'//trim(str_recycle))

  end subroutine InsertJulesFieldsAtVarTable

end module ModJulesFields
