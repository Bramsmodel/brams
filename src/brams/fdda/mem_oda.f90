!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_oda

  use ModNamelistFile, only: namelistFile

  use grid_dims, only: maxfiles !INTENT(IN)

  type oda_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)

     real, pointer, contiguous :: uk(:,:,:) => null()
     real, pointer, contiguous :: vk(:,:,:) => null()
     real, pointer, contiguous :: tk(:,:,:) => null()
     real, pointer, contiguous :: rk(:,:,:) => null()
     real, pointer, contiguous :: ukv(:,:,:) => null()
     real, pointer, contiguous :: vkv(:,:,:) => null()
     real, pointer, contiguous :: tkv(:,:,:) => null()
     real, pointer, contiguous :: rkv(:,:,:) => null()
  end type oda_vars

  type(oda_vars), allocatable, target :: oda_g(:)
  type(oda_vars), allocatable, target :: odam_g(:)

  integer, parameter :: maxodafiles=maxfiles !1000
  integer, parameter :: maxodasta=2000
  integer, parameter :: maxodagrids=10
  integer, parameter :: maxodanzp=200
  integer, parameter :: maxodatimes=3*maxodafiles

  include "files.h"

  character(len=f_name_length) :: fnames_upa(maxodafiles)
  character(len=f_name_length) :: fnames_sfc(maxodafiles)
  character(len=14)  :: itotdate_upa(maxodafiles)
  character(len=14)  :: itotdate_sfc(maxodafiles)
  character(len=8)   :: staid_sfc(maxodasta)
  character(len=8)   :: staid_upa(maxodasta)
  integer            :: ntimes_sfc(maxodasta)
  integer            :: ntimes_upa(maxodasta)
  integer            :: maxtimes_sfc
  integer            :: maxtimes_upa

  ! Namelist inputs

  character(len=128) :: oda_upaprefix ! from RAMSIN
  character(len=128) :: oda_sfcprefix ! from RAMSIN
  integer            :: if_oda ! from RAMSIN
  real               :: frqoda ! from RAMSIN
  real               :: todabeg ! from RAMSIN
  real               :: todaend ! from RAMSIN
  real               :: tnudoda ! from RAMSIN
  real               :: wt_oda_grid(maxodagrids) ! from RAMSIN
  real               :: oda_sfc_til ! from RAMSIN
  real               :: oda_sfc_tel ! from RAMSIN
  real               :: oda_upa_til ! from RAMSIN
  real               :: oda_upa_tel ! from RAMSIN
  real               :: wt_oda_uv ! from RAMSIN
  real               :: wt_oda_th ! from RAMSIN
  real               :: wt_oda_pi ! from RAMSIN
  real               :: wt_oda_rt ! from RAMSIN
  real               :: roda_sfce(maxodagrids) ! from RAMSIN
  real               :: roda_sfc0(maxodagrids) ! from RAMSIN
  real               :: roda_upae(maxodagrids) ! from RAMSIN
  real               :: roda_upa0(maxodagrids) ! from RAMSIN
  real               :: roda_zfact(maxodagrids) ! from RAMSIN
  real               :: roda_hgt(maxodagrids) ! from RAMSIN

  integer            :: nupafiles
  integer            :: nsfcfiles
  integer            :: num_oda_sfc
  integer            :: num_oda_upa

  ! Surface data 

  type oda_sfc_info_type
     character(len=8) :: id
     integer :: intid
     integer :: ntimes
     integer :: iactive(maxodagrids)
     real    :: xista(maxodagrids)
     real    :: xjsta(maxodagrids)
     real    :: xlat
     real    :: xlon
     real    :: xsta
     real    :: ysta
     real    :: stopo
  end type oda_sfc_info_type


  type oda_sfc_type
     real, pointer, contiguous :: temp(:)
     real, pointer, contiguous :: dewpt(:)
     real, pointer, contiguous :: us(:)
     real, pointer, contiguous :: vs(:)
     real, pointer, contiguous :: ps(:)
     real, pointer, contiguous :: u(:)
     real, pointer, contiguous :: v(:)
     real, pointer, contiguous :: time(:)
  end type oda_sfc_type

  type(oda_sfc_info_type), allocatable, target :: oda_sfc_info(:)
  type(oda_sfc_type)     , allocatable, target :: oda_sfc_obs(:)


  ! Upper air info

  type oda_upa_info_type
     character(len=8) :: id
     integer :: intid
     integer :: ntimes
     integer, dimension(maxodagrids) :: iactive
     real, dimension(maxodagrids) :: xista, xjsta
     real :: xlat,xlon,xsta,ysta,stopo
  end type oda_upa_info_type

  ! Upper air data
  type oda_upa_type
     character(len=14) :: ctotdate
     real, pointer, dimension(:,:) :: theta, rv, us, vs, u, v, zz, pi, zgeo
     real, pointer, dimension(:) :: time 
     integer, pointer, dimension(:) :: lp,lz
  end type oda_upa_type

  type(oda_upa_info_type), allocatable :: oda_upa_info(:)
  type(oda_upa_type)     , allocatable :: oda_upa_obs(:)

  integer, parameter :: maxupalevs=100

  ! Krigging routine info

  real :: rmaxkrg(maxodanzp,maxodagrids)  &
       ,ckrg(3,maxodagrids),akrg(maxodanzp,maxodagrids)  &
       ,caxkrg(9,maxodagrids),caykrg(9,maxodagrids),cazkrg(9,maxodagrids)
  integer :: nstkrg(maxodagrids)


  ! Filled obs arrays for an analysis time

  integer, parameter :: maxkobs=10000
  real, dimension(maxkobs) :: ukobs,vkobs,tkobs,rkobs,pkobs  &
       ,xkobs,ykobs,zkobs,ekobs  &
       ,ikobs,jkobs

contains

  subroutine alloc_oda(oda,n1,n2,n3,ng,proc_type)

    implicit none
    type (oda_vars) :: oda
    integer, intent(in) :: n1,n2,n3,ng,proc_type

    ! Allocate arrays based on options (if necessary)


    if( if_oda == 1 .and. proc_type /= 1) then
       allocate (oda%uk(n1,n2,n3))
       allocate (oda%vk(n1,n2,n3))
       allocate (oda%tk(n1,n2,n3))
       allocate (oda%rk(n1,n2,n3))
       allocate (oda%ukv(n1,n2,n3))
       allocate (oda%vkv(n1,n2,n3))
       allocate (oda%tkv(n1,n2,n3))
       allocate (oda%rkv(n1,n2,n3))
    endif

    return
  end subroutine alloc_oda


  subroutine nullify_oda(oda)

    implicit none
    type (oda_vars) :: oda

    if (associated(oda%uk))      nullify (oda%uk)
    if (associated(oda%vk))      nullify (oda%vk)
    if (associated(oda%tk))      nullify (oda%tk)
    if (associated(oda%rk))      nullify (oda%rk)
    if (associated(oda%ukv))     nullify (oda%ukv)
    if (associated(oda%vkv))     nullify (oda%vkv)
    if (associated(oda%tkv))     nullify (oda%tkv)
    if (associated(oda%rkv))     nullify (oda%rkv)

    return
  end subroutine nullify_oda

  subroutine dealloc_oda(oda)

    implicit none
    type (oda_vars) :: oda


    if (associated(oda%uk))      deallocate (oda%uk)
    if (associated(oda%vk))      deallocate (oda%vk)
    if (associated(oda%tk))      deallocate (oda%tk)
    if (associated(oda%rk))      deallocate (oda%rk)
    if (associated(oda%ukv))     deallocate (oda%ukv)
    if (associated(oda%vkv))     deallocate (oda%vkv)
    if (associated(oda%tkv))     deallocate (oda%tkv)
    if (associated(oda%rkv))     deallocate (oda%rkv)

    return
  end subroutine dealloc_oda


  subroutine filltab_oda(oneVarTable, oneVarTableSize, &
       oda, odam, imean)

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(oda_vars), pointer, intent(in) :: oda
    type(oda_vars), pointer, intent(in) :: odam
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(filltab_oda)**"

    if (.not. associated(oda)) then
       call fatal_error(h//" invoked with unassociated oda")
    else if (.not. associated(odam)) then
       call fatal_error(h//" invoked with unassociated odam")
    else if (.not. associated(oneVarTable)) then
       call fatal_error(h//" invoked with unasociated oneVarTable")
    end if

    ! Fill pointers to arrays into variable tables


    if (associated(oda%uk)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oda%uk, &
            'UKODA :3:', &
            odam%uk, imean)
    end if

    if (associated(oda%vk)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oda%vk, &
            'VKODA :3:', &
            odam%vk, imean)
    end if

    if (associated(oda%tk)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oda%tk, &
            'TKODA :3:', &
            odam%tk, imean)
    end if

    if (associated(oda%rk)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oda%rk, &
            'RKODA :3:', &
            odam%rk, imean)
    end if

    if (associated(oda%ukv)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oda%ukv, &
            'UVODA :3:', &
            odam%ukv, imean)
    end if

    if (associated(oda%vkv)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oda%vkv, &
            'VVODA :3:', &
            odam%vkv, imean)
    end if

    if (associated(oda%tkv)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oda%tkv, &
            'TVODA :3:', &
            odam%tkv, imean)
    end if

    if (associated(oda%rkv)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oda%rkv, &
            'RVODA :3:', &
            odam%rkv, imean)
    end if

  end subroutine filltab_oda

  subroutine StoreNamelistFileAtMem_oda(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    frqoda = oneNamelistFile%frqoda
    if_oda = oneNamelistFile%if_oda
    oda_sfc_tel = oneNamelistFile%oda_sfc_tel
    oda_sfc_til = oneNamelistFile%oda_sfc_til
    oda_sfcprefix = oneNamelistFile%oda_sfcprefix
    oda_upa_tel = oneNamelistFile%oda_upa_tel
    oda_upa_til = oneNamelistFile%oda_upa_til
    oda_upaprefix = oneNamelistFile%oda_upaprefix
    roda_hgt = oneNamelistFile%roda_hgt
    roda_sfc0 = oneNamelistFile%roda_sfc0
    roda_sfce = oneNamelistFile%roda_sfce
    roda_upa0 = oneNamelistFile%roda_upa0
    roda_upae = oneNamelistFile%roda_upae
    roda_zfact = oneNamelistFile%roda_zfact
    tnudoda = oneNamelistFile%tnudoda
    todabeg = oneNamelistFile%todabeg
    todaend = oneNamelistFile%todaend
    wt_oda_grid = oneNamelistFile%wt_oda_grid
    wt_oda_pi = oneNamelistFile%wt_oda_pi
    wt_oda_rt = oneNamelistFile%wt_oda_rt
    wt_oda_th = oneNamelistFile%wt_oda_th
    wt_oda_uv = oneNamelistFile%wt_oda_uv
  end subroutine StoreNamelistFileAtMem_oda
end module mem_oda
