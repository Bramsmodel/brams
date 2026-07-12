!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_varinit

  use ModNamelistFile, only: &
       namelistFile

  use grid_dims, only: &
       maxfiles, &
       maxgrds

  use chem1_list, only: &
       nspecies

  use mem_chem1, only: &
       chem_assim

  ! Memory for varfile, history, and condensate nudging

  type varinit_vars
     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, contiguous :: varup(:,:,:)
     real, pointer, contiguous :: varvp(:,:,:)
     real, pointer, contiguous :: varpp(:,:,:)
     real, pointer, contiguous :: vartp(:,:,:)
     real, pointer, contiguous :: varrp(:,:,:)
     real, pointer, contiguous :: varuf(:,:,:)
     real, pointer, contiguous :: varvf(:,:,:)
     real, pointer, contiguous :: varpf(:,:,:)
     real, pointer, contiguous :: vartf(:,:,:)
     real, pointer, contiguous :: varrf(:,:,:)
     real, pointer, contiguous :: varwts(:,:,:)
     real, pointer, contiguous :: varwts_chem(:,:,:)
     real, pointer, contiguous :: varrph(:,:,:)
     real, pointer, contiguous :: varrfh(:,:,:)
     real, pointer, contiguous :: varcph(:,:,:)
     real, pointer, contiguous :: varcfh(:,:,:)
  end type varinit_vars

  type (varinit_vars), pointer :: varinit_g(:) => null()


  integer, parameter :: maxnudfiles=maxfiles !500

  integer :: nud_type ! from RAMSIN
  integer :: nnudfiles
  integer :: nnudfl
  integer :: nudlat ! from RAMSIN

  include "files.h"

  character(len=f_name_length) :: fnames_nud(maxnudfiles)
  character(len=14)  :: itotdate_nud(maxnudfiles)
  real               :: nud_times(maxnudfiles)
  character(len=f_name_length) :: nud_hfile ! from RAMSIN
  real :: timeWindowIAU  ! from RAMSIN
  real :: ramp ! from RAMSIN
  real :: tnudlat ! from RAMSIN
  real :: tnudcent ! from RAMSIN
  real :: tnudtop ! from RAMSIN
  real :: znudtop ! from RAMSIN
  real :: wt_nudge_uv ! from RAMSIN
  real :: wt_nudge_th ! from RAMSIN
  real :: wt_nudge_pi ! from RAMSIN
  real :: wt_nudge_rt ! from RAMSIN
  real :: wt_nudge_grid(maxgrds) ! from RAMSIN
  real :: htime1
  real :: htime2

!--(DMK-CCATT-INI)-----------------------------------------------------
  real :: wt_nudge_chem(nspecies)
!--(DMK-CCATT-FIM)-----------------------------------------------------


  integer :: igrid_match(maxgrds)

  !----------------------------------------------------------------------------
  integer :: nud_cond ! from RAMSIN
  integer :: ncondfiles
  integer :: ncondfl

  character(len=f_name_length) :: fnames_cond(maxnudfiles)
  character(len=14)  :: itotdate_cond(maxnudfiles)
  real :: cond_times(maxnudfiles)
  character(len=f_name_length) :: cond_hfile ! from RAMSIN
  real :: tcond_beg  ! from RAMSIN
  real :: tcond_end  ! from RAMSIN
  real :: wt_nudgec_grid(maxgrds) ! from RAMSIN
  real :: t_nudge_rc  ! from RAMSIN
  real :: condtime1
  real :: condtime2


  !----------------------------------------------------------------------------
  character(len=f_name_length) :: fnames_varf(maxnudfiles)
  character(len=14)  :: itotdate_varf(maxnudfiles)
  real :: varf_times(maxnudfiles)

  character(len=f_name_length)       :: varfpfx ! from RAMSIN
  character(len=20)       :: varfpfx_vartype ! from RAMSIN
  ! Modif. by ALF

  real :: vtime1
  real :: vtime2
  real :: vwait1 ! from RAMSIN
  real :: vwaittot ! from RAMSIN

  integer :: nvarffiles
  integer :: nvarffl


  character(len=14)       :: lastdate_iv

  !----------------------------------------------------------------------------

contains

  subroutine alloc_varinit(varinit,n1,n2,n3)

    implicit none
    type (varinit_vars) :: varinit
    integer, intent(in) :: n1,n2,n3

    ! Allocate arrays based on options (if necessary)

    if( nud_type >= 1 ) then
       allocate (varinit%varup(n1,n2,n3))
       allocate (varinit%varvp(n1,n2,n3))
       allocate (varinit%varpp(n1,n2,n3))
       allocate (varinit%vartp(n1,n2,n3))
       allocate (varinit%varrp(n1,n2,n3))
       allocate (varinit%varuf(n1,n2,n3))
       allocate (varinit%varvf(n1,n2,n3))
       allocate (varinit%varpf(n1,n2,n3))
       allocate (varinit%vartf(n1,n2,n3))
       allocate (varinit%varrf(n1,n2,n3))                      
       allocate (varinit%varwts(n1,n2,n3))
       if(chem_assim == 1 ) allocate (varinit%varwts_chem(n1,n2,n3))
    endif

    if (nud_cond == 1) then
       allocate (varinit%varcph(n1,n2,n3))
       allocate (varinit%varcfh(n1,n2,n3))                      
       allocate (varinit%varrph(n1,n2,n3))
       allocate (varinit%varrfh(n1,n2,n3))                      
    endif
  end subroutine alloc_varinit


  subroutine nullify_varinit(varinit)

    implicit none
    type (varinit_vars) :: varinit


    if (associated(varinit%varup))     nullify (varinit%varup)
    if (associated(varinit%varvp))     nullify (varinit%varvp)
    if (associated(varinit%varpp))     nullify (varinit%varpp)
    if (associated(varinit%vartp))     nullify (varinit%vartp)
    if (associated(varinit%varrp))     nullify (varinit%varrp)
    if (associated(varinit%varuf))     nullify (varinit%varuf)
    if (associated(varinit%varvf))     nullify (varinit%varvf)
    if (associated(varinit%varpf))     nullify (varinit%varpf)
    if (associated(varinit%vartf))     nullify (varinit%vartf)
    if (associated(varinit%varrf))     nullify (varinit%varrf)
    if (associated(varinit%varwts))    nullify (varinit%varwts)
    if (associated(varinit%varwts_chem)) nullify (varinit%varwts_chem)
    if (associated(varinit%varcph))     nullify (varinit%varcph)
    if (associated(varinit%varcfh))     nullify (varinit%varcfh)
    if (associated(varinit%varrph))     nullify (varinit%varrph)
    if (associated(varinit%varrfh))     nullify (varinit%varrfh)
  end subroutine nullify_varinit

  subroutine dealloc_varinit(varinit)

    implicit none
    type (varinit_vars) :: varinit


    if (associated(varinit%varup))     deallocate (varinit%varup)
    if (associated(varinit%varvp))     deallocate (varinit%varvp)
    if (associated(varinit%varpp))     deallocate (varinit%varpp)
    if (associated(varinit%vartp))     deallocate (varinit%vartp)
    if (associated(varinit%varrp))     deallocate (varinit%varrp)
    if (associated(varinit%varuf))     deallocate (varinit%varuf)
    if (associated(varinit%varvf))     deallocate (varinit%varvf)
    if (associated(varinit%varpf))     deallocate (varinit%varpf)
    if (associated(varinit%vartf))     deallocate (varinit%vartf)
    if (associated(varinit%varrf))     deallocate (varinit%varrf)
    if (associated(varinit%varwts))    deallocate (varinit%varwts)

!--(DMK-CCATT-FIM)-----------------------------------------------------
    if (associated(varinit%varwts_chem)) deallocate(varinit%varwts_chem)
!--(DMK-CCATT-FIM)-----------------------------------------------------

    if (associated(varinit%varcph))     deallocate (varinit%varcph)
    if (associated(varinit%varcfh))     deallocate (varinit%varcfh)
    if (associated(varinit%varrph))     deallocate (varinit%varrph)
    if (associated(varinit%varrfh))     deallocate (varinit%varrfh)

    return
  end subroutine dealloc_varinit


  subroutine filltab_varinit(oneVarTable, oneVarTableSize, &
       varinit, imean, ng)

    ! Build VarTable entry with varinit_vars components

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(varinit_vars), pointer, intent(in) :: varinit(:)
    integer, intent(in) :: imean
    integer, intent(in) :: ng

    integer, parameter :: notIMean=0
    real, pointer, contiguous :: null3D(:,:,:) => null()
    character(len=*), parameter :: h="**(filltab_varinit)**"

    ! Fill pointers to arrays into variable tables

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(varinit)) then
       call fatal_error(h//" varinit not associated")
    end if

    if (associated(varinit(ng)%varup))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varup, &
            'VARUP :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%varvp))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varvp, &
            'VARVP :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%varpp))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varpp, &
            'VARPP :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%vartp))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%vartp, &
            'VARTP :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%varrp))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varrp, &
            'VARRP :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%varuf))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varuf, &
            'VARUF :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%varvf))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varvf, &
            'VARVF :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%varpf))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varpf, &
            'VARPF :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%vartf))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%vartf, &
            'VARTF :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%varrf))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varrf, &
            'VARRF :3:mpti', &
            null3D, &
            notIMean)
    end if

    if (associated(varinit(ng)%varwts))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varwts, &
            'VARWTS :3:mpti', &
            null3D, &
            notIMean)
    end if

    if(chem_assim == 1 .and. associated(varinit(ng)%varwts_chem)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            varinit(ng)%varwts_chem, &
            'VARWTS_CHEM :3:mpti', &
            null3D, &
            notIMean)
    end if


    if (nud_cond == 1) then
       if (associated(varinit(ng)%varcph))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               varinit(ng)%varcph, &
               'VARCPH :3:mpti', &
               null3D, &
               notIMean)
       end if

       if (associated(varinit(ng)%varcfh))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               varinit(ng)%varcfh, &
               'VARCFH :3:mpti', &
               null3D, &
               notIMean)
       end if

       if (associated(varinit(ng)%varrph))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               varinit(ng)%varrph, &
               'VARRPH :3:mpti', &
               null3D, &
               notIMean)
       end if

       if (associated(varinit(ng)%varrfh))  then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               varinit(ng)%varrfh, &
               'VARRFH :3:mpti', &
               null3D, &
               notIMean)
       end if
    endif
  end subroutine filltab_varinit


  subroutine StoreNamelistFileAtMem_varinit(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    cond_hfile = oneNamelistFile%cond_hfile
    nud_cond = oneNamelistFile%nud_cond
    nud_hfile = oneNamelistFile%nud_hfile
    nud_type = oneNamelistFile%nud_type
    nudlat = oneNamelistFile%nudlat
    t_nudge_rc = oneNamelistFile%t_nudge_rc
    tcond_beg = oneNamelistFile%tcond_beg
    tcond_end = oneNamelistFile%tcond_end
    tnudcent = oneNamelistFile%tnudcent
    timeWindowIAU = oneNamelistFile%timeWindowIAU
    ramp = oneNamelistFile%ramp
    tnudlat = oneNamelistFile%tnudlat
    tnudtop = oneNamelistFile%tnudtop
    varfpfx = oneNamelistFile%varfpfx
    varfpfx_vartype = oneNamelistFile%varfpfx_vartype
    vwait1 = oneNamelistFile%vwait1
    vwaittot = oneNamelistFile%vwaittot
    wt_nudge_grid = oneNamelistFile%wt_nudge_grid
    wt_nudge_pi = oneNamelistFile%wt_nudge_pi
    wt_nudge_rt = oneNamelistFile%wt_nudge_rt
    wt_nudge_th = oneNamelistFile%wt_nudge_th
    wt_nudge_uv = oneNamelistFile%wt_nudge_uv
    wt_nudgec_grid = oneNamelistFile%wt_nudgec_grid
    znudtop = oneNamelistFile%znudtop
  end subroutine StoreNamelistFileAtMem_varinit

end module mem_varinit
