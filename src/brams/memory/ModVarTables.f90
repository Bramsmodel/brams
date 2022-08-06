!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModVarTables

  use ModParallelEnvironment, only: &
       MsgDump

  use io_params, only: &
       nlite_vars, & ! INTENT(IN)
       lite_vars     ! INTENT(IN)

  use chem1_list, only:&
       chem_name=>spc_name,    &
       chem_alloc=>spc_alloc,  & 
       chem_on=>on,            &
       chem_fdda=>fdda,        &
       chem_transport=>transport, &
       chem_nspecies=>nspecies 

  use aer1_list, only: aer_name=>spc_name,     &
       aer_nspecies=>nspecies, &
       aer_alloc=>spc_alloc,     &
       aer_fdda=>fdda, &
       aer_transport=>transport, &
       aer_on=>on,               &
       aer_nmodes=>nmodes

  use ModVarTable, only: &
       VarTable
  
  implicit none
  private
  public :: maxvars
  public :: VarTableFields
  public :: vtab_r
  public :: nvgrids
  public :: num_var
!!$  public :: InsertVTab
  public :: GetVTabEntry
  public :: GetVTabSectionSize
  public :: VerifyVTabEntry
  public :: StringIndexing
  public :: ZeroVTab
  public :: lite_varset
  public :: DumpVTab
  public :: setInitial4Vtable
  public :: DeepCopyToVarTable
  public :: DeepCopyFromVarTable

  include "constants.h"
  character(len=32) :: nameTo=""
  character(len=32) :: nameFrom=""

  ! Maximum number of variables of all types (3d + 2d + leaf)

  integer, parameter :: maxvars=1000

  ! Define data type for main variable table

  type VarTableFields
     real, pointer      :: var_m_2D(:,:) => null()
     ! var_m_2D: pointer to full 2D field average
     real, pointer      :: var_m_3D(:,:,:) => null()
     ! var_m_3D: pointer to full 3D field average
     real, pointer      :: var_m_4D(:,:,:,:) => null()
     ! var_m_4D: pointer to full 4D field average
     real, pointer      :: var_p_2D(:,:) => null()
     ! var_p_2D: pointer to full real 2D field 
     real, pointer      :: var_p_3D(:,:,:) => null()
     ! var_p_3D: pointer to full real 3D field 
     real, pointer      :: var_p_4D(:,:,:,:) => null()
     ! var_p_4D: pointer to full real 4D field 
     integer, pointer   :: var_p_2D_I(:,:) => null()
     ! var_p_2D_I: pointer to full integer 2D field 
     integer, pointer   :: var_m_2D_I(:,:) => null()
     ! var_p_2D_I: pointer to full integer 2D field average
     integer            :: idim_type
     ! idim_type codes dimensionality of var_p:
     ! idim_type == 2 means (nmxp, nmyp)
     ! idim_type == 3 means (nmzp, nmxp, nmyp)
     ! idim_type == 4 means (nzg, nmxp, nmyp, npatch)
     ! idim_type == 5 means (nzs, nmxp, nmyp, npatch)
     ! idim_type == 6 means (nmxp, nmyp, npatch)
     ! idim_type == 7 means (nmxp, nmyp, nwave)
     integer(kind=i8)   :: npts
     ! npts is number of elements at field (product of dimensions)
     integer            :: ihist
     integer            :: ianal
     integer            :: imean
     integer            :: ilite
     integer            :: impti     ! commmunicate during initialization
     integer            :: impt1
     integer            :: impt2
     integer            :: impt3     ! communicate for output
     integer            :: imptd
     integer            :: irecycle
     character (len=16) :: name
     ! field name
  end type VarTableFields

  ! Main variable table allocated to (maxvars,maxgrds)

  type(VarTableFields), allocatable, target :: vtab_r(:,:)

  ! "nvgrids" is "ngrids", for convenience

  integer :: nvgrids

  ! number of variables for each grid, allocated to "ngrids"

  integer, allocatable :: num_var(:)


!!$  interface InsertVTab
!!$     module procedure InsertVTab_2D
!!$     module procedure InsertVTab_2D_I
!!$     module procedure InsertVTab_3D
!!$     module procedure InsertVTab_4D
!!$  end interface InsertVTab

  interface ZeroVTab
     module procedure zero_vtab_2D
     module procedure zero_vtab_2D_I
     module procedure zero_vtab_3D
     module procedure zero_vtab_4D
  end interface ZeroVTab


contains


  ! insert a variable into vtab_r


!!$  subroutine InsertTab_2D(var, varm, ng, npts, imean, tabstr)
!!$    real, target,      intent(in) :: var(:,:)
!!$    real, target,      intent(in) :: varm(:,:)
!!$    integer,           intent(in) :: ng
!!$    integer(kind=i8),  intent(in) :: npts
!!$    integer,           intent(in) :: imean
!!$    character (len=*), intent(in) :: tabstr
!!$
!!$    ! insert the old way
!!$
!!$    call vtables2(var(1,1), varm(1,1), ng, npts, imean, tabstr)
!!$
!!$    ! save full field and average field
!!$
!!$    vtab_r(num_var(ng),ng)%var_p_2D => var
!!$    vtab_r(num_var(ng),ng)%var_m_2D => varm
!!$  end subroutine InsertVTab_2D
!!$
!!$
!!$
!!$  subroutine InsertVTab_2D_I(var, varm, ng, npts, imean, tabstr)
!!$    integer, target,   intent(in) :: var(:,:)
!!$    integer, target,   intent(in) :: varm(:,:)
!!$    integer,           intent(in) :: ng
!!$    integer(kind=i8),  intent(in) :: npts
!!$    integer,           intent(in) :: imean
!!$    character (len=*), intent(in) :: tabstr
!!$
!!$    character (len=80) :: line
!!$    character (len=1) :: toksep=':', cdimen,ctype
!!$    character (len=32) :: tokens(10)
!!$    character (len=8) :: cname,ctab
!!$    integer :: ntok,nt,nv
!!$    
!!$    ! insert text of vtables2 and adapt to integer fields
!!$
!!$    call tokenize1(tabstr,tokens,ntok,toksep)
!!$
!!$    num_var(ng)=num_var(ng)+1
!!$    nv=num_var(ng)
!!$
!!$    vtab_r(nv,ng)%name=tokens(1)
!!$    vtab_r(nv,ng)%npts=npts
!!$    read(tokens(2),*) vtab_r(nv,ng)%idim_type
!!$
!!$    vtab_r(nv,ng)%ihist=0
!!$    vtab_r(nv,ng)%ianal=0
!!$    vtab_r(nv,ng)%imean=imean
!!$    vtab_r(nv,ng)%ilite=0
!!$    vtab_r(nv,ng)%impti=0
!!$    vtab_r(nv,ng)%impt1=0
!!$    vtab_r(nv,ng)%impt2=0
!!$    vtab_r(nv,ng)%impt3=0
!!$    vtab_r(nv,ng)%imptd=0
!!$    vtab_r(nv,ng)%irecycle=0
!!$
!!$    do nt=3,ntok
!!$       ctab=tokens(nt)
!!$
!!$       if(ctab == 'hist' ) then
!!$          vtab_r(nv,ng)%ihist=1
!!$       elseif(ctab == 'anal' ) then
!!$          vtab_r(nv,ng)%ianal=1
!!$       elseif(ctab == 'lite' ) then
!!$          vtab_r(nv,ng)%ilite=1
!!$       elseif(ctab == 'mpti' ) then
!!$          vtab_r(nv,ng)%impti=1
!!$       elseif(ctab == 'mpt1' ) then
!!$          vtab_r(nv,ng)%impt1=1
!!$       elseif(ctab == 'mpt2' ) then
!!$          vtab_r(nv,ng)%impt2=1
!!$       elseif(ctab == 'mpt3' ) then
!!$          vtab_r(nv,ng)%impt3=1
!!$       elseif(ctab == 'mptd' ) then
!!$          vtab_r(nv,ng)%imptd=1
!!$       elseif(ctab == 'recycle' ) then
!!$          vtab_r(nv,ng)%irecycle=1
!!$       else
!!$          print*, 'Illegal table specification for var:', tokens(1),ctab
!!$          stop 'bad var table'
!!$       endif
!!$
!!$    enddo
!!$
!!$    ! save full field and average field
!!$
!!$    vtab_r(num_var(ng),ng)%var_p_2D_I => var
!!$    vtab_r(num_var(ng),ng)%var_m_2D_I => varm
!!$  end subroutine InsertVTab_2D_I
!!$
!!$
!!$
!!$  subroutine InsertVTab_3D(var, varm, ng, npts, imean, tabstr)
!!$    real, target,      intent(in) :: var(:,:,:)
!!$    real, target,      intent(in) :: varm(:,:,:)
!!$    integer,           intent(in) :: ng
!!$    integer(kind=i8),  intent(in) :: npts
!!$    integer,           intent(in) :: imean
!!$    character (len=*), intent(in) :: tabstr
!!$    !write(*,*) 'LFR-DEBUG: 3D: ',tabstr
!!$
!!$    ! insert the old way
!!$
!!$    call vtables2(var(1,1,1), varm(1,1,1), ng, npts, imean, tabstr)
!!$
!!$    ! save full field and average field
!!$
!!$    vtab_r(num_var(ng),ng)%var_p_3D => var
!!$    vtab_r(num_var(ng),ng)%var_m_3D => varm
!!$  end subroutine InsertVTab_3D
!!$
!!$
!!$
!!$  subroutine InsertVTab_4D(var, varm, ng, npts, imean, tabstr)
!!$    real, target,      intent(in) :: var(:,:,:,:)
!!$    real, target,      intent(in) :: varm(:,:,:,:)
!!$    integer,           intent(in) :: ng
!!$    integer(kind=i8),  intent(in) :: npts
!!$    integer,           intent(in) :: imean
!!$    character (len=*), intent(in) :: tabstr
!!$
!!$    ! insert the old way
!!$
!!$    call vtables2(var(1,1,1,1), varm(1,1,1,1), ng, npts, imean, tabstr)
!!$
!!$    ! save full field and average field
!!$
!!$    vtab_r(num_var(ng),ng)%var_p_4D => var
!!$    vtab_r(num_var(ng),ng)%var_m_4D => varm
!!$  end subroutine InsertVTab_4D




  subroutine GetVTabEntry(tabstr, ng, vtabPtr)
    character (len=*), intent(in) :: tabstr
    integer,           intent(in) :: ng
    type(VarTableFields), pointer   :: vtabPtr

    integer :: ni

    vtabPtr => null()
    do ni = 1, num_var(ng)
       if (trim(vtab_r(ni,ng)%name) == trim(tabstr)) then
          vtabPtr => vtab_r(ni,ng)
          exit
       end if
    end do
  end subroutine GetVTabEntry



  subroutine vtables2(var, varm, ng, npts, imean, tabstr)
    real, target :: var,varm
    integer, intent(in) :: ng,imean !npts
    integer(kind=i8), intent(in) :: npts
    character (len=*), intent(in) :: tabstr

    character (len=80) ::line
    character (len=1) ::toksep=':', cdimen,ctype
    character (len=32) ::tokens(10)
    character (len=8) :: cname,ctab

    integer :: ntok,nt,nv

    call tokenize1(tabstr,tokens,ntok,toksep)
    !print *,'LFR->DBG: vtables2: ',npts,tokens(1); call flush(6)
    num_var(ng)=num_var(ng)+1
    nv=num_var(ng)

    vtab_r(nv,ng)%name=tokens(1)
    vtab_r(nv,ng)%npts=npts
    read(tokens(2),*) vtab_r(nv,ng)%idim_type
    !print*,'tab:',nv,ng,vtab_r(nv,ng)%name ,vtab_r(nv,ng)%npts

    vtab_r(nv,ng)%ihist=0
    vtab_r(nv,ng)%ianal=0
    vtab_r(nv,ng)%imean=imean
    vtab_r(nv,ng)%ilite=0
    vtab_r(nv,ng)%impti=0
    vtab_r(nv,ng)%impt1=0
    vtab_r(nv,ng)%impt2=0
    vtab_r(nv,ng)%impt3=0

    !--(DMK)------------------------------------------
    vtab_r(nv,ng)%imptd=0
    !--(DMK)------------------------------------------

    vtab_r(nv,ng)%irecycle=0

    do nt=3,ntok
       ctab=tokens(nt)

       if(ctab == 'hist' ) then
          vtab_r(nv,ng)%ihist=1
       elseif(ctab == 'anal' ) then
          vtab_r(nv,ng)%ianal=1
       elseif(ctab == 'lite' ) then
          vtab_r(nv,ng)%ilite=1
       elseif(ctab == 'mpti' ) then
          vtab_r(nv,ng)%impti=1
       elseif(ctab == 'mpt1' ) then
          vtab_r(nv,ng)%impt1=1
       elseif(ctab == 'mpt2' ) then
          vtab_r(nv,ng)%impt2=1
       elseif(ctab == 'mpt3' ) then
          vtab_r(nv,ng)%impt3=1

          !--(DMK)------------------------------------------
       elseif(ctab == 'mptd' ) then
          vtab_r(nv,ng)%imptd=1
          !--(DMK)------------------------------------------

       elseif(ctab == 'recycle' ) then
          vtab_r(nv,ng)%irecycle=1
       else
          print*, 'Illegal table specification for var:', tokens(1),ctab
          stop 'bad var table'
       endif

    enddo

    return
  end subroutine vtables2





  integer function GetVTabSectionSize(vTabPtr, &
       iStart, iEnd, jStart, jEnd)
    type(VarTableFields), pointer :: vTabPtr
    integer, intent(in) :: iStart
    integer, intent(in) :: iEnd
    integer, intent(in) :: jStart
    integer, intent(in) :: jEnd

    character(len=8) :: c0
    character(len=*), parameter :: h="**(GetVTabSectionSize)**"

    if (.not. associated(vTabPtr)) then
       call fatal_error(h//" vTabPtr not associated")
    else
       GetVTabSectionSize=(iEnd-iStart+1)*(jEnd-jStart+1)
       select case (vTabPtr%idim_type)
       case(2)
          ! idim_type == 2 means (nmxp, nmyp)
       case(3)
          ! idim_type == 3 means (nmzp, nmxp, nmyp)
          GetVTabSectionSize=GetVTabSectionSize*&
               size(vTabPtr%var_p_3D,1)
       case(4)
          ! idim_type == 4 means (nzg, nmxp, nmyp, npatch)
          GetVTabSectionSize=GetVTabSectionSize*&
               size(vTabPtr%var_p_4D,1)*&
               size(vTabPtr%var_p_4D,4)
       case(5)
          ! idim_type == 5 means (nzs, nmxp, nmyp, npatch)
          GetVTabSectionSize=GetVTabSectionSize*&
               size(vTabPtr%var_p_4D,1)*&
               size(vTabPtr%var_p_4D,4)
       case(6)
          ! idim_type == 6 means (nmxp, nmyp, npatch)
          GetVTabSectionSize=GetVTabSectionSize*&
               size(vTabPtr%var_p_3D,3)
       case(7)
          ! idim_type == 7 means (nmxp, nmyp, nwave)
          GetVTabSectionSize=GetVTabSectionSize*&
               size(vTabPtr%var_p_3D,3)
       case default
          write(c0,"(i8)") vTabPtr%idim_type
          call fatal_error(h//" unknown idim_type="//trim(adjustl(c0)))
       end select
    end if
  end function GetVTabSectionSize




  subroutine VerifyVTabEntry(vTabPtr)
    type(VarTableFields), pointer :: vTabPtr
    character(len=*), parameter :: h="**(VerifyVTabEntry)**"

    if (.not. associated(vTabPtr)) then
       call fatal_error(h//" null vTabPtr")
    else
       select case (vTabPtr%idim_type)
       case (2)
          if (.not. associated(vTabPtr%var_p_2D)) then
             call fatal_error(h//" vTabPtr%var_p_2D of field "//&
                  trim(adjustl(vTabPtr%name))//" is not associated")
          end if
       case (3)
          if (.not. associated(vTabPtr%var_p_3D)) then
             call fatal_error(h//" vTabPtr%var_p_3D of field "//&
                  trim(adjustl(vTabPtr%name))//" is not associated")
          end if
       case (4:5)
          if (.not. associated(vTabPtr%var_p_4D)) then
             call fatal_error(h//" vTabPtr%var_p_4D of field "//&
                  trim(adjustl(vTabPtr%name))//" is not associated")
          end if
       case (6:7)
          if (.not. associated(vTabPtr%var_p_3D)) then
             call fatal_error(h//" vTabPtr%var_p_3D of field "//&
                  trim(adjustl(vTabPtr%name))//" is not associated")
          end if
       case default
          call fatal_error(h//" vTabPtr%idim_type of field "//&
               trim(adjustl(vTabPtr%name))//" is outside range [2:7]")
       end select
    end if
  end subroutine VerifyVTabEntry


  subroutine StringIndexing(vTabPtr, &
       xStart, xEnd, yStart, yEnd, string)
    type(VarTableFields), pointer :: vTabPtr
    integer, intent(in) :: xStart
    integer, intent(in) :: xEnd
    integer, intent(in) :: yStart
    integer, intent(in) :: yEnd
    character(len=*), intent(out) :: string

    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(StringIndexing)**"


    if (.not. associated(vTabPtr)) then
       call fatal_error(h//" null vTabPtr")
    end if

    write(c0,"(i8)") xStart
    write(c1,"(i8)") xEnd
    write(c2,"(i8)") yStart
    write(c3,"(i8)") yEnd

    select case (vTabPtr%idim_type)
    case(2)
       string="("//&
            trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
            trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
    case(3)
       write(c4,"(i8)") size(vTabPtr%var_p_3D,1)
       string="(1:"//trim(adjustl(c4))//","//&
            trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
            trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
    case(4:5)
       write(c4,"(i8)") size(vTabPtr%var_p_4D,1)
       write(c5,"(i8)") size(vTabPtr%var_p_4D,4)
       string="(1:"//trim(adjustl(c4))//","//&
            trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
            trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
            "1:"//trim(adjustl(c5))//")"

    case(6:7)
       write(c4,"(i8)") size(vTabPtr%var_p_3D,3)
       string="("//&
            trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
            trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
            "1:"//trim(adjustl(c4))//")"

    case default
       write(c0,"(i8)") vTabPtr%idim_type
       call fatal_error(h//" field section "//trim(vTabPtr%name)//&
            " with unknown idim_type="//trim(adjustl(c0)))
    end select
  end subroutine StringIndexing


  subroutine zero_vtab_2D(var,nx,ny)
    integer, intent(in) :: nx,ny
    real, intent(inout) :: var(nx,ny)

    var=0.0

  end subroutine zero_vtab_2D

  subroutine zero_vtab_2D_I(var,nx,ny)
    integer, intent(in) :: nx,ny
    integer, intent(inout) :: var(nx,ny)

    var=0.0

  end subroutine zero_vtab_2D_I

  subroutine zero_vtab_3D(var,nx,ny,nz)
    integer,intent(in) :: nx,ny,nz
    real, intent(inout) :: var(nx,ny,nz)

    var=0.0

  end subroutine zero_vtab_3D

  subroutine zero_vtab_4D(var,nx,ny,nz,nk)
    integer,intent(in) :: nx,ny,nz,nk
    real, intent(inout) :: var(nx,ny,nz,nk)

    var=0.0

  end subroutine zero_vtab_4D



  subroutine lite_varset(proc_type)

    ! Arguments:
    integer, intent(in) :: proc_type

    ! Local variables:
    integer :: nv,ng,nvl,ifound


    ! Loop over each variable input in namelist "LITE_VARS" and set
    !   lite flag in ModVarTables

    do ng = 1,nvgrids   
       vtab_r(1:num_var(ng),ng)%ilite = 0
    enddo

    do nvl=1,nlite_vars
       ifound=0

       do ng=1,nvgrids

          do nv=1,num_var(ng)

             if (vtab_r(nv,ng)%name == lite_vars(nvl) ) then
                vtab_r(nv,ng)%ilite = 1
                ifound=1
             endif

          enddo

       enddo

       if (proc_type==0 .or. proc_type==1) then !Output only in Master Process
          if(ifound == 0) then
             print*,'!---------------------------------------------------------'
             print*,'! LITE_VARS variable does not exist in main variable table'
             print*,'!    variable name-->',lite_vars(nvl),'<--'
             print*,'!---------------------------------------------------------'
          else
             print*,'!---------------------------------------------------------'
             print*,'! LITE_VARS variable added--->',trim(lite_vars(nvl))
             print*,'!---------------------------------------------------------'
          endif
       endif

    enddo

    return
  end subroutine lite_varset

  !-------------------------------------------------------------------------


!!$subroutine GetVarFromMem (nxp, nyp, nzp, nzg, nzs, npatch, &
!!$     varName, itype, ngrd, arrayOut, sizeArray)
!!$
!!$  integer,            intent(in)    :: nxp  ! as at vartable
!!$  integer,            intent(in)    :: nyp  ! as at vartable
!!$  integer,            intent(in)    :: nzp  ! as at vartable
!!$  integer,            intent(in)    :: nzg  ! as at vartable
!!$  integer,            intent(in)    :: nzs  ! as at vartable
!!$  integer,            intent(in)    :: npatch  ! as at vartable
!!$  character(LEN=*),   intent(in)    :: varName
!!$  integer,            intent(in)    :: ngrd
!!$  integer,            intent(out)   :: itype
!!$  integer(kind=i8),   intent(in)    :: sizeArray
!!$  real,	              intent(inout) :: arrayOut(sizeArray)
!!$
!!$  character(len=16) :: c0, c1
!!$  character(len=*), parameter :: h="**(GetVarFromMem)**"
!!$  integer(kind=i8) :: ni
!!$  integer(kind=i8) :: npts
!!$  logical          :: found
!!$  character(len=len(varName)) :: varnIn, varnOut
!!$  real, pointer :: ptr ! points to a field at vartable
!!$  real :: scr1(sizeArray)
!!$
!!$  ! field name changes from vartable to analysis file
!!$  ! in two cases (PI and HKH); 
!!$  ! given output file field name, find vartable correspondent
!!$
!!$  if (trim(varName) == 'PI') then
!!$     varnIn = 'PP'
!!$  else if (trim(varName) == 'HKH') then
!!$     varnIn = 'HKM'
!!$  else
!!$     varnIn = varName
!!$  end if
!!$
!!$  ! search for vartable name at vartable
!!$  ! store result at array 
!!$
!!$  found = .false.
!!$  do ni = 1, num_var(ngrd)
!!$     if (trim(vtab_r(ni,ngrd)%name) == trim(varnIn)) then
!!$        itype = vtab_r(ni,ngrd)%idim_type
!!$        npts  = vtab_r(ni,ngrd)%npts
!!$        ptr => vtab_r(ni,ngrd)%var_p
!!$        if (npts > sizeArray) then
!!$           write(c0,"(i16)") npts
!!$           write(c1,"(i16)") sizeArray
!!$           call fatal_error(h//&
!!$                " array size for "//trim(varName)//&
!!$                " is "//trim(adjustl(c1))//&
!!$                ", smaller than "//trim(adjustl(c0))//" required")
!!$        end if
!!$        found = .true.
!!$        exit
!!$     end if
!!$  end do
!!$
!!$  ! halts if not there
!!$
!!$  if (.not. found) then
!!$     write(*,"(a)") h//" var "//trim(varnIn)//" not found in vtab_r; will dump vtab_r"
!!$     call DumpVTab(ngrd)
!!$     call fatal_error(h//" var "//trim(varnIn)//" not found in vtab_r")
!!$  end if
!!$
!!$  ! convert fields PP, HKM and VKH from vartables to analysis file
!!$  ! or store field at scr1
!!$
!!$  call PreProcForOutput(ngrd, varnIn, npts, ptr, scr1, varnOut)
!!$
!!$  ! verify output name
!!$
!!$  if (trim(varnOut) /= trim (varName)) then
!!$     call fatal_error(h//" fails computing "//trim(varName))
!!$  end if
!!$
!!$  ! move verticals from first to third dimension, if required
!!$
!!$  if (itype==3 .or. itype==4 .or. itype==5) then
!!$     call RearrangeForOutput(nxp, nyp, nzp, nzg, nzs, npatch, &
!!$          itype, scr1, arrayOut)
!!$  else
!!$     arrayOut = scr1
!!$  end if
!!$end subroutine GetVarFromMem

  subroutine DumpVTab(ngrd)
    integer, intent(in) :: ngrd

    integer :: i
    character(len=*), parameter :: h="**(DumpVTab)**"

    write(*,"(a,i2)") h//" dump of vtab_r names for grid ",ngrd
    write(*,"(a)") h//" name            idim_type           npts"
    do i = 1, num_var(ngrd)
       write(*,"(1x,a16,1x,i8,1x,i16)") vtab_r(i,ngrd)%name, &
            vtab_r(i,ngrd)%idim_type, vtab_r(i,ngrd)%npts
    end do
  end subroutine DumpVTab





  subroutine setInitial4Vtable(ng, chemistry, aerosol)
    integer, intent(in) :: ng
    integer, intent(in) :: chemistry
    integer, intent(in) :: aerosol
    integer :: ni, nspc, imode
    character(len=2) :: cmode

    do ni = 1, num_var(ng)

       vtab_r(ni,ng)%ianal=0	

       if (trim(vtab_r(ni,ng)%name) == 'TOPT'  .or. &
	    trim(vtab_r(ni,ng)%name) == 'UP'    .or. &
	    trim(vtab_r(ni,ng)%name) == 'VP'    .or. &
	    trim(vtab_r(ni,ng)%name) == 'THETA' .or. &
	    trim(vtab_r(ni,ng)%name) == 'PP'    .or. &
	    trim(vtab_r(ni,ng)%name) == 'RV') then
          vtab_r(ni,ng)%ianal=1
          cycle
       end if

       if(vtab_r(ni,ng)%irecycle == 1) vtab_r(ni,ng)%ianal=1


       if(CHEMISTRY >= 0) then 
          do nspc=1,chem_nspecies
             !print*, spc_alloc(fdda,nspc), on
             !print*, trim(vtab_r(ni,ng)%name), '>>', trim(spc_name(nspc))//'P'
             if(chem_alloc(chem_fdda,nspc) == chem_on .and. &
                  trim(vtab_r(ni,ng)%name) == trim(chem_name(nspc))//'P') then 
                vtab_r(ni,ng)%ianal=1
                cycle
             end if
          end do
       end if
       if(AEROSOL == 1 .and. CHEMISTRY >= 0) then
          do nspc=1,aer_nspecies
             do imode = 1, aer_nmodes
                write(cmode, '(BN, I2)')imode
                cmode = adjustl(cmode)
                !print*, trim(vtab_r(ni,ng)%name), trim(aer_name(nspc))//trim(cmode)//'P'
                if(aer_alloc(aer_fdda,imode,nspc) == 1  .and. &
                     trim(vtab_r(ni,ng)%name) == trim(aer_name(nspc))//trim(cmode)//'P') then
                   vtab_r(ni,ng)%ianal=1
                   cycle
                end if
             end do
          end do
       end if
    end do

  end subroutine setInitial4Vtable

  subroutine DeepCopyToVarTable(oneVarTable, oneVarTableSize, name)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(out) :: oneVarTableSize
    character(len=*), intent(in) :: name

    integer :: i
    character(len=*), parameter :: h="**(DeepCopyToVarTable)**"

    if (nameTo /= "") then
       call fatal_error(h//" consecutive calls; this one from "//trim(adjustl(name))//&
            " and previous one from "//trim(adjustl(nameTo)))
    else
       nameTo=name
       nameFrom=""
    end if
    
    do i = 1, num_var(1)
       oneVarTable(i)%var_p_2D => vtab_r(i,1)%var_p_2D
       oneVarTable(i)%var_p_2D_I => vtab_r(i,1)%var_p_2D_I
       oneVarTable(i)%var_p_3D => vtab_r(i,1)%var_p_3D
       oneVarTable(i)%var_p_4D => vtab_r(i,1)%var_p_4D
       oneVarTable(i)%var_m_2D => vtab_r(i,1)%var_m_2D
       oneVarTable(i)%var_m_2D_I => vtab_r(i,1)%var_m_2D_I
       oneVarTable(i)%var_m_3D => vtab_r(i,1)%var_m_3D
       oneVarTable(i)%var_m_4D => vtab_r(i,1)%var_m_4D
       oneVarTable(i)%idim_type = vtab_r(i,1)%idim_type
       oneVarTable(i)%npts = vtab_r(i,1)%npts
       oneVarTable(i)%ihist = vtab_r(i,1)%ihist
       oneVarTable(i)%ianal = vtab_r(i,1)%ianal
       oneVarTable(i)%imean = vtab_r(i,1)%imean
       oneVarTable(i)%ilite = vtab_r(i,1)%ilite
       oneVarTable(i)%impti = vtab_r(i,1)%impti
       oneVarTable(i)%impt1 = vtab_r(i,1)%impt1
       oneVarTable(i)%impt2 = vtab_r(i,1)%impt2
       oneVarTable(i)%impt3 = vtab_r(i,1)%impt3
       oneVarTable(i)%imptd = vtab_r(i,1)%imptd
       oneVarTable(i)%irecycle = vtab_r(i,1)%irecycle
       oneVarTable(i)%name = vtab_r(i,1)%name
       call MsgDump(h//" copyed "//trim(adjustl(oneVarTable(i)%name)))
    end do
    oneVarTableSize = num_var(1)
  end subroutine DeepCopyToVarTable





  subroutine DeepCopyFromVarTable(oneVarTable, oneVarTableSize, name)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(out) :: oneVarTableSize
    character(len=*), intent(in) :: name

    integer :: i
    character(len=*), parameter :: h="**(DeepCopyFromVarTable)**"

    if (nameFrom /= "") then
       call fatal_error(h//" consecutive calls; this one from "//trim(adjustl(name))//&
            " and previous one from "//trim(adjustl(nameFrom)))
    else
       nameFrom=name
       nameTo=""
    end if
    
    do i = 1, oneVarTableSize
       vtab_r(i,1)%var_p_2D => oneVarTable(i)%var_p_2D  
       vtab_r(i,1)%var_p_2D_I => oneVarTable(i)%var_p_2D_I  
       vtab_r(i,1)%var_p_3D => oneVarTable(i)%var_p_3D  
       vtab_r(i,1)%var_p_4D => oneVarTable(i)%var_p_4D  
       vtab_r(i,1)%var_m_2D => oneVarTable(i)%var_m_2D  
       vtab_r(i,1)%var_m_2D_I => oneVarTable(i)%var_m_2D_I  
       vtab_r(i,1)%var_m_3D => oneVarTable(i)%var_m_3D  
       vtab_r(i,1)%var_m_4D => oneVarTable(i)%var_m_4D  
       vtab_r(i,1)%idim_type = oneVarTable(i)%idim_type  
       vtab_r(i,1)%npts = oneVarTable(i)%npts  
       vtab_r(i,1)%ihist = oneVarTable(i)%ihist  
       vtab_r(i,1)%ianal = oneVarTable(i)%ianal  
       vtab_r(i,1)%imean = oneVarTable(i)%imean  
       vtab_r(i,1)%ilite = oneVarTable(i)%ilite  
       vtab_r(i,1)%impti = oneVarTable(i)%impti  
       vtab_r(i,1)%impt1 = oneVarTable(i)%impt1  
       vtab_r(i,1)%impt2 = oneVarTable(i)%impt2  
       vtab_r(i,1)%impt3 = oneVarTable(i)%impt3  
       vtab_r(i,1)%imptd = oneVarTable(i)%imptd  
       vtab_r(i,1)%irecycle = oneVarTable(i)%irecycle  
       vtab_r(i,1)%name = oneVarTable(i)%name  
       call MsgDump(h//" copyed "//trim(adjustl(vtab_r(i,1)%name)))
    end do
    num_var(1) = oneVarTableSize
  end subroutine DeepCopyFromVarTable  
end module ModVarTables
     


