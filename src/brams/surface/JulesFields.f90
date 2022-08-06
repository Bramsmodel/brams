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

  use ModVarTable, only: &
       VarTable, &
       InsertAtVarTable
  
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
    res%gpp = 0.0
    
    allocate (res%resp_s(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate resp_s fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%resp_s = 0.0
    
    allocate (res%resp_p(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate resp_p fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%resp_p = 0.0
    
    allocate (res%npp(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate npp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%npp = 0.0
    
    allocate (res%u10mj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate u10mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%u10mj = 0.0
    
    allocate (res%v10mj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate v10mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%v10mj = 0.0
    
    allocate (res%u10mj1hr(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate u10mj1hr fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%u10mj1hr = 0.0
    
    allocate (res%v10mj1hr(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate v10mj1hr fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%v10mj1hr = 0.0
    
    allocate (res%fracj(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate fracj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%fracj = 0.0
    
    allocate (res%t2mj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate t2mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%t2mj = 0.0
    
    allocate (res%t2mj_max(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate t2mj_max fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%t2mj_max = 0.0
    
    allocate (res%t2mj_min(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate t2mj_min fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%t2mj_min = 0.0
    
    allocate (res%rv2mj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate rv2mj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%rv2mj = 0.0
    
    allocate (res%csj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate csj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%csj = 0.0
    
    allocate (res%anthrop_heatj(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate anthrop_heatj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%anthrop_heatj = 0.0
    
    allocate (res%radnet_tilej(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate radnet_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%radnet_tilej = 0.0
    
    allocate (res%ftl_tilej(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate ftl_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%ftl_tilej = 0.0
    
    allocate (res%le_tilej(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate le_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%le_tilej = 0.0
    
    allocate (res%htf_tilej(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate htf_tilej fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%htf_tilej = 0.0
    
    allocate (res%snowdepthj(mxp,myp,npatch), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate snowdepthj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%snowdepthj = 0.0
    
    allocate (res%ht_fluxj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate ht_fluxj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%ht_fluxj = 0.0
    
    allocate (res%temp_surfj(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate temp_surfj fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%temp_surfj = 0.0
    
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

  


  subroutine InsertJulesFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneJulesFields, oneAveJulesFields, oneNamelistFile)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(JulesFields), pointer, intent(in) :: oneJulesFields
    type(JulesFields), pointer, intent(in) :: oneAveJulesFields
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile

    logical :: assAve
    logical :: assThis
    character(len=8) :: str_recycle
    character(len=*), parameter :: h="**(InsertJulesFieldsAtVarTable)**" 

    if (.not. associated(oneJulesFields)) then
       call fatal_error(h//" oneJulesFields not associated")
    else if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if


    assAve=associated(oneAveJulesFields)

    if (oneNamelistFile%ipastin == 1) then
       str_recycle = ':recycle'
    else
       str_recycle = ''
    endif

    ! Fill pointers to arrays into variable tables

    if (associated(oneJulesFields%gpp)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%gpp)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%gpp, &
               'GPP :2:anal:mpti:mpt3', &  
               oneAveJulesFields%gpp)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%gpp, &
               'GPP :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%resp_s)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%resp_s)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%resp_s, &
               'RESP_S :2:anal:mpti:mpt3', & 
               oneAveJulesFields%resp_s)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%resp_s, &
               'RESP_S :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%resp_p)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%resp_p)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%resp_p, &
               'RESP_P :2:anal:mpti:mpt3', & 
               oneAveJulesFields%resp_p)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%resp_p, &
               'RESP_P :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%npp)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%npp)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%npp, &
               'NPP :2:anal:mpti:mpt3', &    
               oneAveJulesFields%npp)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%npp, &
               'NPP :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%u10mj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%u10mj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%u10mj, &
               'U10MJ :2:anal:mpti:mpt3', &
               oneAveJulesFields%u10mj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%u10mj, &
               'U10MJ :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%v10mj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%v10mj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%v10mj, &
               'V10MJ :2:anal:mpti:mpt3', &
               oneAveJulesFields%v10mj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%v10mj, &
               'V10MJ :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%t2mj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%t2mj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%t2mj, &
               'T2MJ :2:anal:mpti:mpt3', &
               oneAveJulesFields%t2mj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%t2mj, &
               'T2MJ :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%t2mj_max)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%t2mj_max)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%t2mj_max, &
               'T2MJ_MAX :2:anal:mpti:mpt3', &
               oneAveJulesFields%t2mj_max)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%t2mj_max, &
               'T2MJ_MAX :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%t2mj_min)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%t2mj_min)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%t2mj_min, &
               'T2MJ_MIN :2:anal:mpti:mpt3', &
               oneAveJulesFields%t2mj_min)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%t2mj_min, &
               'T2MJ_MIN :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%rv2mj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%rv2mj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%rv2mj, &
               'RV2MJ :2:anal:mpti:mpt3', &
               oneAveJulesFields%rv2mj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%rv2mj, &
               'RV2MJ :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%csj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%csj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%csj, &
               'CSJ :2:anal:mpti:mpt3', &
               oneAveJulesFields%csj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%csj, &
               'CSJ :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%ht_fluxj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%ht_fluxj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%ht_fluxj, &
               'ht_fluxj :2:anal:mpti:mpt3', &
               oneAveJulesFields%ht_fluxj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%ht_fluxj, &
               'ht_fluxj :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%temp_surfj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%temp_surfj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%temp_surfj, &
               'temp_surfj :2:anal:mpti:mpt3', &
               oneAveJulesFields%temp_surfj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%temp_surfj, &
               'temp_surfj :2:anal:mpti:mpt3')
       end if
    end if
          

    if (associated(oneJulesFields%u10mj1hr)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%u10mj1hr)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%u10mj1hr, &
               'U10MJ1hr :3:hist:anal:mpti:mpt3:mpt2', &
               oneAveJulesFields%u10mj1hr)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%u10mj1hr, &
               'U10MJ1hr :3:hist:anal:mpti:mpt3:mpt2')
       end if
    end if
          

    if (associated(oneJulesFields%v10mj1hr)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%v10mj1hr)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%v10mj1hr, &
               'V10MJ1hr :3:hist:anal:mpti:mpt3:mpt2', &
               oneAveJulesFields%v10mj1hr)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%v10mj1hr, &
               'V10MJ1hr :3:hist:anal:mpti:mpt3:mpt2')
       end if
    end if
          

    if (associated(oneJulesFields%fracj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%fracj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%fracj, &
               'fracj :3:hist:anal:mpti:mpt3:mpt2', &
               oneAveJulesFields%fracj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%fracj, &
               'fracj :3:hist:anal:mpti:mpt3:mpt2')
       end if
    end if
          

    if (associated(oneJulesFields%anthrop_heatj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%anthrop_heatj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%anthrop_heatj, &
               'anthrop_heatj :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
               oneAveJulesFields%anthrop_heatj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%anthrop_heatj, &
               'anthrop_heatj :6:hist:anal:mpti:mpt3'//trim(str_recycle))
       end if
    end if
          

    if (associated(oneJulesFields%radnet_tilej)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%radnet_tilej)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%radnet_tilej, &
               'radnet_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
               oneAveJulesFields%radnet_tilej)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%radnet_tilej, &
               'radnet_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))
       end if
    end if
          

    if (associated(oneJulesFields%ftl_tilej)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%ftl_tilej)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%ftl_tilej, &
               'ftl_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
               oneAveJulesFields%ftl_tilej)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%ftl_tilej, &
               'ftl_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))
       end if
    end if
          

    if (associated(oneJulesFields%le_tilej)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%le_tilej)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%le_tilej, &
               'le_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
               oneAveJulesFields%le_tilej)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%le_tilej, &
               'le_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))
       end if
    end if
          

    if (associated(oneJulesFields%htf_tilej)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%htf_tilej)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%htf_tilej, &
               'htf_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
               oneAveJulesFields%htf_tilej)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%htf_tilej, &
               'htf_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))
       end if
    end if
          

    if (associated(oneJulesFields%snowdepthj)) then
       if (assAve) then
          assThis=associated(oneAveJulesFields%snowdepthj)
       else
          assThis=.false.
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%snowdepthj, &
               'snowdepthj :6:hist:anal:mpti:mpt3'//trim(str_recycle), &
               oneAveJulesFields%snowdepthj)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               oneJulesFields%snowdepthj, &
               'snowdepthj :6:hist:anal:mpti:mpt3'//trim(str_recycle))
       end if
    end if
          

  end subroutine InsertJulesFieldsAtVarTable

end module ModJulesFields
