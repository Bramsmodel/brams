!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module var_tables

  use ModParallelEnvironment, only: &
       MsgDump
  
  implicit none
  private
  public :: maxvars
  public :: var_tables_r
  public :: vtab_r
  public :: nvgrids
  public :: num_var
  public :: scalar_table
  public :: scalar_tab
  public :: num_scalar
  public :: InsertVTab
  public :: InsertScalarTab  
  public :: GetVTabEntry
  public :: GetVTabSectionSize
  public :: VerifyVTabEntry
  public :: StringIndexing
  public :: ZeroVTab

  include "constants.h"

  ! Maximum number of variables of all types (3d + 2d + leaf)

  integer, parameter :: maxvars=1000

  ! Define data type for main variable table

  type var_tables_r
     real, pointer      :: var_p
     ! var_p: scalar pointing to field first position
     real, pointer      :: var_m
     ! var_m: scalar pointing to field average first position
     real, pointer      :: var_p_2D(:,:) => null()
     ! var_p_2D: pointer to full real 2D field 
     real, pointer      :: var_p_3D(:,:,:) => null()
     ! var_p_3D: pointer to full real 3D field 
     real, pointer      :: var_p_4D(:,:,:,:) => null()
     ! var_p_4D: pointer to full real 4D field 
     integer, pointer   :: var_p_2D_I(:,:) => null()
     ! var_p_2D_I: pointer to full integer 2D field 
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
  end type var_tables_r

  ! Main variable table allocated to (maxvars,maxgrds)

  type(var_tables_r), allocatable, target :: vtab_r(:,:)

  ! "nvgrids" is "ngrids", for convenience

  integer :: nvgrids

  ! number of variables for each grid, allocated to "ngrids"

  integer, allocatable :: num_var(:)

  ! Define data type for scalar variable table

  type scalar_table
     real, pointer      :: var_p
     real, pointer      :: var_t
     character (len=16) :: name
     real, pointer      :: a_var_p(:)
     real, pointer      :: a_var_t(:)
     real, pointer      :: var_p_2D(:,:) => null()
     real, pointer      :: var_p_3D(:,:,:) => null()
  end type scalar_table

  ! Scalar variable table allocated to (maxsclr,maxgrds)

  type(scalar_table), allocatable :: scalar_tab(:,:)

  ! number of scalars for each grid, allocated to "ngrids"

  integer, allocatable :: num_scalar(:)


  interface InsertVTab
     module procedure InsertVTab_2D
     module procedure InsertVTab_2D_I
     module procedure InsertVTab_3D
     module procedure InsertVTab_4D
  end interface InsertVTab

  interface ZeroVTab
     module procedure zero_vtab_2D
     module procedure zero_vtab_2D_I
     module procedure zero_vtab_3D
     module procedure zero_vtab_4D
  end interface ZeroVTab

  !LFR em 06mai2020 
  interface
     subroutine vtables2_I(var, varm, ng, npts, imean, tabstr)
       !use var_tables
       !implicit none
       include "constants.h"
       integer, target :: var,varm
       integer, intent(in) :: ng,imean !npts
       integer(kind=i8), intent(in) :: npts
       character (len=*), intent(in) :: tabstr
     end subroutine vtables2_I
  end interface

  interface InsertScalarTab
     module procedure InsertScalarTab_2D
     module procedure InsertScalarTab_3D
  end interface InsertScalarTab

contains





  function StringScalarTableEntry(oneScalarTable) result(str)
    type(scalar_table), intent(in) :: oneScalarTable
    character(len=64) :: str

    character(len=*), parameter :: assoc =" assoc "
    character(len=*), parameter :: notAssoc =" not assoc "

    str="variable "//trim(oneScalarTable%name)//":"
    if (associated(oneScalarTable%a_var_p)) then
       str = trim(str)//" a_var_p"//assoc
    else
       str = trim(str)//" a_var_p"//notAssoc
    end if

    if (associated(oneScalarTable%a_var_t)) then
       str = trim(str)//", a_var_t"//assoc
    else
       str = trim(str)//", a_var_t"//notAssoc
    end if

    if (associated(oneScalarTable%var_p_2D)) then
       str = trim(str)//", var_p_2D"//assoc
    else
       str = trim(str)//", var_p_2D"//notAssoc
    end if

    if (associated(oneScalarTable%var_p_3D)) then
       str = trim(str)//", var_p_3D"//assoc
    else
       str = trim(str)//", var_p_3D"//notAssoc
    end if
  end function StringScalarTableEntry


  
  subroutine InsertScalarTab_2D(varp,vart,ng,tabstr,elements)
    real, target :: varp(:,:)
    real, target :: vart(:)
    integer, intent(in) :: ng
    character (len=*), intent(in) :: tabstr
    integer, intent(in) :: elements

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertScalarTab_2D)**"
    logical, parameter :: dumpLocal=.true.

    if (dumpLocal) then
       write(str(1),"(i8)") ng
       write(str(2),"(i8)") num_scalar(ng)
       call MsgDump(h//" at grid "//trim(adjustl(str(1)))//&
            " inserting field "//trim(tabstr)//"; last scalar_tab used position is "//&
            trim(adjustl(str(2))))
    end if
    
    ! insert the old way for vtables_scalar
    call vtables_scalar(varp(1,1),vart(1),ng,tabstr)

    ! insert the old way for vtables_scalar_new
    call vtables_scalar_new(varp,vart,ng,tabstr,elements)

    ! save full field vtables_scalar
    scalar_tab(num_scalar(ng),ng)%var_p_2D => varp

!!$    ! save  field vtables_scalar_new
!!$    scalar_tab(num_scalar(ng),ng)%a_var_t => vart(:)

    if (dumpLocal) then
       write(str(1),"(i8)") num_scalar(ng)
       call MsgDump(h//" inserted at scalar_tab position "//&
            trim(adjustl(str(1)))//" "//&
            StringScalarTableEntry(scalar_tab(num_scalar(ng),ng)))
    end if
  end subroutine InsertScalarTab_2D





  
  subroutine InsertScalarTab_3D(varp,vart,ng,tabstr,elements)
    real, target :: varp(:,:,:)
    real, target :: vart(:)
    integer, intent(in) :: ng
    character (len=*), intent(in) :: tabstr
    integer, intent(in) :: elements

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertScalarTab_2D)**"
    logical, parameter :: dumpLocal=.true.

    if (dumpLocal) then
       write(str(1),"(i8)") ng
       write(str(2),"(i8)") num_scalar(ng)
       call MsgDump(h//" at grid "//trim(adjustl(str(1)))//&
            " inserting field "//trim(tabstr)//"; last scalar_tab used position is "//&
            trim(adjustl(str(2))))
    end if
    
    ! insert the old way for vtables_scalar
    call vtables_scalar(varp(1,1,1),vart(1),ng,tabstr)

    ! insert the old way for vtables_scalar_new
    call vtables_scalar_new(varp,vart,ng,tabstr,elements)

    ! save full field vtables_scalar
    scalar_tab(num_scalar(ng),ng)%var_p_3D => varp

!!$    ! save  field vtables_scalar_new
!!$    scalar_tab(num_scalar(ng),ng)%a_var_t => vart(:)

    if (dumpLocal) then
       write(str(1),"(i8)") num_scalar(ng)
       call MsgDump(h//" inserted at scalar_tab position "//&
            trim(adjustl(str(1)))//" "//&
            StringScalarTableEntry(scalar_tab(num_scalar(ng),ng)))
    end if
  end subroutine InsertScalarTab_3D

  ! insert a variable into vtab_r


  subroutine InsertVTab_2D(var, varm, ng, npts, imean, tabstr)
    real, target,      intent(in) :: var(:,:)
    real, target,      intent(in) :: varm(:,:)
    integer,           intent(in) :: ng
    integer(kind=i8),  intent(in) :: npts
    integer,           intent(in) :: imean
    character (len=*), intent(in) :: tabstr

    ! insert the old way

    call vtables2(var(1,1), varm(1,1), ng, npts, imean, tabstr)

    ! save full field

    vtab_r(num_var(ng),ng)%var_p_2D => var
  end subroutine InsertVTab_2D



  subroutine InsertVTab_2D_I(var, varm, ng, npts, imean, tabstr)
    integer, target,   intent(in) :: var(:,:)
    integer, target,   intent(in) :: varm(:,:)
    integer,           intent(in) :: ng
    integer(kind=i8),  intent(in) :: npts
    integer,           intent(in) :: imean
    character (len=*), intent(in) :: tabstr

    ! insert the old way

    call vtables2_I(var(1,1), varm(1,1), ng, npts, imean, tabstr)

    ! save full field

    vtab_r(num_var(ng),ng)%var_p_2D_I => var
  end subroutine InsertVTab_2D_I



  subroutine InsertVTab_3D(var, varm, ng, npts, imean, tabstr)
    real, target,      intent(in) :: var(:,:,:)
    real, target,      intent(in) :: varm(:,:,:)
    integer,           intent(in) :: ng
    integer(kind=i8),  intent(in) :: npts
    integer,           intent(in) :: imean
    character (len=*), intent(in) :: tabstr
    !write(*,*) 'LFR-DEBUG: 3D: ',tabstr

    ! insert the old way

    call vtables2(var(1,1,1), varm(1,1,1), ng, npts, imean, tabstr)

    ! save full field

    vtab_r(num_var(ng),ng)%var_p_3D => var
  end subroutine InsertVTab_3D



  subroutine InsertVTab_4D(var, varm, ng, npts, imean, tabstr)
    real, target,      intent(in) :: var(:,:,:,:)
    real, target,      intent(in) :: varm(:,:,:,:)
    integer,           intent(in) :: ng
    integer(kind=i8),  intent(in) :: npts
    integer,           intent(in) :: imean
    character (len=*), intent(in) :: tabstr

    ! insert the old way

    call vtables2(var(1,1,1,1), varm(1,1,1,1), ng, npts, imean, tabstr)

    ! save full field

    vtab_r(num_var(ng),ng)%var_p_4D => var
  end subroutine InsertVTab_4D




  subroutine GetVTabEntry(tabstr, ng, vtabPtr)
    character (len=*), intent(in) :: tabstr
    integer,           intent(in) :: ng
    type(var_tables_r), pointer   :: vtabPtr

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

    vtab_r(nv,ng)%var_p => var
    vtab_r(nv,ng)%var_m => varm


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

  subroutine vtables_scalar(varp,vart,ng,tabstr)
    implicit none
    real, target :: varp,vart
    integer, intent(in) :: ng
    character (len=*), intent(in) :: tabstr

    character (len=1) ::toksep=':'
    character (len=32) ::tokens(10)

    integer :: ntok,nv,ns

    call tokenize1(tabstr,tokens,ntok,toksep)
    ! increase position counter
    num_scalar(ng)=num_scalar(ng)+1
    nv=num_scalar(ng)
    scalar_tab(nv,ng)%name = tokens(1)
    scalar_tab(nv,ng)%var_p => varp
    scalar_tab(nv,ng)%var_t => vart
  end subroutine vtables_scalar



  
  subroutine vtables_scalar_new(varp,vart,ng,tabstr,elements)
    implicit none
    integer :: elements !ALF
    real, target :: varp(elements), vart(elements)
    integer, intent(in) :: ng
    character (len=*), intent(in) :: tabstr

    character(len=1) :: toksep=':'
    character(len=32) :: tokens(10)

    integer :: ntok,nv,ns, i

    call tokenize1(tabstr,tokens,ntok,toksep)
    !    Fill in existing table slot
    nv=num_scalar(ng)
    scalar_tab(nv,ng)%a_var_p => varp(:)
    scalar_tab(nv,ng)%a_var_t => vart(:)
  end subroutine vtables_scalar_new



  
  integer function GetVTabSectionSize(vTabPtr, &
       iStart, iEnd, jStart, jEnd)
    type(var_tables_r), pointer :: vTabPtr
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
    type(var_tables_r), pointer :: vTabPtr
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
    type(var_tables_r), pointer :: vTabPtr
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


end module var_tables

subroutine vtables2_I(var, varm, ng, npts, imean, tabstr)
  use var_tables
  implicit none
  include "constants.h"
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

  num_var(ng)=num_var(ng)+1
  nv=num_var(ng)

  vtab_r(nv,ng)%var_p => var
  vtab_r(nv,ng)%var_m => varm


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
end subroutine vtables2_I
