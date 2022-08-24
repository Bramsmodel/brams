!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModVarTable

  use iso_fortran_env, only: &
       int64

  use ModParallelEnvironment, only: &
       MsgDump

  use io_params, only: &
       nlite_vars, & ! INTENT(IN)
       lite_vars     ! INTENT(IN)

  implicit none

  private

  public :: VarTable
  public :: CreateVarTable
  public :: DestroyVarTable
  public :: DumpVarTable
  public :: InsertAtVarTable
  public :: InsertVarTable
  public :: InsertVarTable_2D
  public :: FixVarTableForIOUTPUT5
  public :: MarkLiteVarsAtVarTable
  public :: Name2VarTableEntry

  ! Maximum number of variables of all types (3d + 2d + leaf)

  integer, parameter :: maxvars=1000

  ! Define data type for main variable table

  type VarTable
     real, pointer      :: var_p_2D(:,:) => null()
     ! var_p_2D: pointer to full real 2D field 
     integer, pointer   :: var_p_2D_I(:,:) => null()
     ! var_p_2D_I: pointer to full integer 2D field 
     real, pointer      :: var_p_3D(:,:,:) => null()
     ! var_p_3D: pointer to full real 3D field 
     real, pointer      :: var_p_4D(:,:,:,:) => null()
     ! var_p_4D: pointer to full real 4D field 
     real, pointer      :: var_m_2D(:,:) => null()
     ! var_m_2D: pointer to full 2D field average
     integer, pointer   :: var_m_2D_I(:,:) => null()
     ! var_p_2D_I: pointer to full integer 2D field average
     real, pointer      :: var_m_3D(:,:,:) => null()
     ! var_m_3D: pointer to full 3D field average
     real, pointer      :: var_m_4D(:,:,:,:) => null()
     ! var_m_4D: pointer to full 4D field average
     integer            :: idim_type
     ! idim_type codes dimensionality of var_p:
     ! idim_type == 2 means (nmxp, nmyp)
     ! idim_type == 3 means (nmzp, nmxp, nmyp)
     ! idim_type == 4 means (nzg, nmxp, nmyp, npatch)
     ! idim_type == 5 means (nzs, nmxp, nmyp, npatch)
     ! idim_type == 6 means (nmxp, nmyp, npatch)
     ! idim_type == 7 means (nmxp, nmyp, nwave)
     integer(kind=int64)   :: npts ! number of grid cells
     ! npts is number of elements at field (product of dimensions)
     integer :: ihist ! history write flag (=1 means to write)
     integer :: ianal ! analysis write flag (=1 means to write)
     integer :: imean ! mean value field  (=1 means is an average field)
     integer :: ilite ! lite write flat (=1 means to write on lite output)
     integer :: impti ! commmunicate during initialization
     integer :: impt1 ! ghost zone update demanded by timestep procedure
     integer :: impt2 ! unused at BRAMS6
     integer :: impt3 ! unused at BRAMS6
     integer :: imptd ! unused at BRAMS6
     integer :: irecycle ! recycle from previous run during initialization 
     character (len=16) :: name ! field name
  end type VarTable

  interface InsertAtVarTable
     module procedure InsertAtVarTable_2D
     module procedure InsertAtVarTable_2D_I
     module procedure InsertAtVarTable_3D
     module procedure InsertAtVarTable_4D
  end interface InsertAtVarTable

  interface InsertVarTable
     module procedure InsertVarTable_2D
     module procedure InsertVarTable_2D_I
     module procedure InsertVarTable_3D
     module procedure InsertVarTable_4D
  end interface InsertVarTable


contains





  function CreateVarTable() result(res)
    type(VarTable), pointer :: res(:)

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateVarTable)**"

    allocate(res(maxvars), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       write(str(2),"(i8)") maxvars
       call fatal_error(h//" allocate res("//&
            trim(adjustl(str(2)))//") fails with stat="//&
            trim(adjustl(str(1))))
    end if
  end function CreateVarTable




  subroutine DestroyVarTable(oneVarTable)
    type(VarTable), pointer, intent(inout) :: oneVarTable(:)

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyVarTable)**"

    if (associated(oneVarTable)) then
       deallocate(oneVarTable, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneVarTable fails with stat="//&
               trim(adjustl(str(1))))
       end if
       nullify(oneVarTable)
    end if
  end subroutine DestroyVarTable





  subroutine DumpVarTable(oneVarTable, oneVarTableSize, name)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(in) :: oneVarTableSize
    character(len=*), intent(in) :: name

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpVarTable)**"

    if (associated(oneVarTable)) then
       write(str(1),"(i8)") oneVarTableSize
       call MsgDump(h//" "//trim(adjustl(name))//" is associated with size="//&
            trim(adjustl(str(1))))
    else
       call MsgDump(h//" "//trim(adjustl(name))//" is not associated")
    end if
  end subroutine DumpVarTable


  ! insert a real 2D field into a VarTable


  subroutine InsertAtVarTable_2D(oneVarTable, oneVarTableSize, var, tabstr, varm)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    real, pointer, intent(in) :: var(:,:)
    character (len=*), intent(in) :: tabstr
    real, pointer, optional, intent(in) :: varm(:,:)

    integer(kind=int64) :: npts
    integer :: imean
    character(len=*), parameter :: h="**(InsertAtVarTable_2D)**"

    ! field size

    npts=int(size(var,1),int64) * &
         int(size(var,2),int64)

    ! average field is present?

    if (present(varm)) then
       imean=1
    else
       imean=0
    end if

    ! get new entry and fill char components

    call NewVarTableEntry(oneVarTable, oneVarTableSize, npts, imean, tabstr)

    ! save full field

    oneVarTable(oneVarTableSize)%var_p_2D => var

    ! if present, save average field

    if (present(varm)) then
       oneVarTable(oneVarTableSize)%var_m_2D => varm
    end if
  end subroutine InsertAtVarTable_2D


  ! insert an integer 2D field into a VarTable


  subroutine InsertAtVarTable_2D_I(oneVarTable, oneVarTableSize, var, tabstr, varm)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    integer, pointer, intent(in) :: var(:,:)
    character (len=*), intent(in) :: tabstr
    integer, pointer, optional, intent(in) :: varm(:,:)

    integer(kind=int64) :: npts
    integer :: imean
    character(len=*), parameter :: h="**(InsertAtVarTable_2D_I)**"

    ! field size

    npts=int(size(var,1),int64) * &
         int(size(var,2),int64)

    ! average field is present?

    if (present(varm)) then
       imean=1
    else
       imean=0
    end if

    ! get new entry and fill char components

    call NewVarTableEntry(oneVarTable, oneVarTableSize, npts, imean, tabstr)

    ! save full field

    oneVarTable(oneVarTableSize)%var_p_2D_I => var

    ! if present, save average field

    if (present(varm)) then
       oneVarTable(oneVarTableSize)%var_m_2D_I => varm
    end if
  end subroutine InsertAtVarTable_2D_I


  ! insert a real 3D field into a VarTable


  subroutine InsertAtVarTable_3D(oneVarTable, oneVarTableSize, var, tabstr, varm)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    real, pointer, intent(in) :: var(:,:,:)
    character (len=*), intent(in) :: tabstr
    real, pointer, optional, intent(in) :: varm(:,:,:)

    integer(kind=int64) :: npts
    integer :: imean
    character(len=*), parameter :: h="**(InsertAtVarTable_3D)**"

    ! field size

    npts=int(size(var,1),int64) * &
         int(size(var,2),int64) * &
         int(size(var,3),int64)

    ! average field is present?

    if (present(varm)) then
       imean=1
    else
       imean=0
    end if

    ! get new entry and fill char components

    call NewVarTableEntry(oneVarTable, oneVarTableSize, npts, imean, tabstr)

    ! save full field

    oneVarTable(oneVarTableSize)%var_p_3D => var

    ! if present, save average field

    if (present(varm)) then
       oneVarTable(oneVarTableSize)%var_m_3D => varm
    end if
  end subroutine InsertAtVarTable_3D


  ! insert a real 4D field into a VarTable


  subroutine InsertAtVarTable_4D(oneVarTable, oneVarTableSize, var, tabstr, varm)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    real, pointer, intent(in) :: var(:,:,:,:)
    character (len=*), intent(in) :: tabstr
    real, pointer, optional, intent(in) :: varm(:,:,:,:)

    integer(kind=int64) :: npts
    integer :: imean
    character(len=*), parameter :: h="**(InsertAtVarTable_4D)**"

    ! field size

    npts=int(size(var,1),int64) * &
         int(size(var,2),int64) * &
         int(size(var,3),int64) * &
         int(size(var,4),int64) 

    ! average field is present?

    if (present(varm)) then
       imean=1
    else
       imean=0
    end if

    ! get new entry and fill char components

    call NewVarTableEntry(oneVarTable, oneVarTableSize, npts, imean, tabstr)

    ! save full field

    oneVarTable(oneVarTableSize)%var_p_4D => var

    ! if present, save average field

    if (present(varm)) then
       oneVarTable(oneVarTableSize)%var_m_4D => varm
    end if
  end subroutine InsertAtVarTable_4D






  subroutine NewVarTableEntry(oneVarTable, oneVarTableSize, npts, imean, tabstr)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    integer(kind=int64), intent(in) :: npts
    integer, intent(in) :: imean
    character (len=*), intent(in) :: tabstr

    integer :: ntok
    integer :: nt
    character (len=1), parameter :: toksep=':'
    character (len=32) :: tokens(10)
    character (len=8) :: ctab
    character(len=8) :: str(10)
    character(len=256) :: strOut
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(NewVarTableEntry)**"

    ! get new entry

    oneVarTableSize = oneVarTableSize + 1
    if (oneVarTableSize > maxvars) then
       write(str(1),"(i8)") maxvars
       call fatal_error(h//" new entry exceeds allocated size of "//&
            trim(adjustl(str(1))))
    end if

    ! parse token string

    call tokenize1(tabstr,tokens,ntok,toksep)

    ! insert tokens at new entry;
    ! name, number of points, idim_type

    oneVarTable(oneVarTableSize)%name=tokens(1)
    oneVarTable(oneVarTableSize)%npts=npts
    read(tokens(2),*) oneVarTable(oneVarTableSize)%idim_type

    if (dumpLocal) then
       strOut=tokens(1)
    end if

    ! further tokens; set defaults and correct

    oneVarTable(oneVarTableSize)%ihist=0
    oneVarTable(oneVarTableSize)%ianal=0
    oneVarTable(oneVarTableSize)%imean=imean
    oneVarTable(oneVarTableSize)%ilite=0
    oneVarTable(oneVarTableSize)%impti=0
    oneVarTable(oneVarTableSize)%impt1=0
    oneVarTable(oneVarTableSize)%impt2=0
    oneVarTable(oneVarTableSize)%impt3=0
    oneVarTable(oneVarTableSize)%imptd=0
    oneVarTable(oneVarTableSize)%irecycle=0
    do nt=3,ntok
       ctab=tokens(nt)
       if(ctab == 'hist' ) then
          oneVarTable(oneVarTableSize)%ihist=1
          if (dumpLocal) strOut=trim(adjustl(strOut))//"; HIST"
       elseif(ctab == 'anal' ) then
          oneVarTable(oneVarTableSize)%ianal=1
          if (dumpLocal) strOut=trim(adjustl(strOut))//"; ANAL"
       elseif(ctab == 'lite' ) then
          oneVarTable(oneVarTableSize)%ilite=1
          if (dumpLocal) strOut=trim(adjustl(strOut))//"; LITE"
       elseif(ctab == 'mpti' ) then
          oneVarTable(oneVarTableSize)%impti=1
          if (dumpLocal) strOut=trim(adjustl(strOut))//"; MPTI"
       elseif(ctab == 'mpt1' ) then
          oneVarTable(oneVarTableSize)%impt1=1
          if (dumpLocal) strOut=trim(adjustl(strOut))//"; MPT1"
       elseif(ctab == 'mpt2' ) then
          oneVarTable(oneVarTableSize)%impt2=1
          if (dumpLocal) strOut=trim(adjustl(strOut))//"; MPT2"
       elseif(ctab == 'mpt3' ) then
          oneVarTable(oneVarTableSize)%impt3=1
          if (dumpLocal) strOut=trim(adjustl(strOut))//"; MPT3"
       elseif(ctab == 'mptd' ) then
          oneVarTable(oneVarTableSize)%imptd=1
          if (dumpLocal) strOut=trim(adjustl(strOut))//"; MPTD"
       elseif(ctab == 'recycle' ) then
          if (dumpLocal) strOut=trim(adjustl(strOut))//"; RECICLE"
          oneVarTable(oneVarTableSize)%irecycle=1
       else
          call fatal_error(h//" unknown token "//trim(adjustl(ctab)))
       endif
    end do
    if (dumpLocal) then
       call MsgDump(h//" "//trim(adjustl(strOut)))
    end if
  end subroutine NewVarTableEntry
!!$
!!$
!!$  subroutine StringIndexing(vTabPtr, &
!!$       xStart, xEnd, yStart, yEnd, string)
!!$    type(VarTableFields), pointer :: vTabPtr
!!$    integer, intent(in) :: xStart
!!$    integer, intent(in) :: xEnd
!!$    integer, intent(in) :: yStart
!!$    integer, intent(in) :: yEnd
!!$    character(len=*), intent(out) :: string
!!$
!!$    character(len=8) :: c0, c1, c2, c3, c4, c5
!!$    character(len=*), parameter :: h="**(StringIndexing)**"
!!$
!!$
!!$    if (.not. associated(vTabPtr)) then
!!$       call fatal_error(h//" null vTabPtr")
!!$    end if
!!$
!!$    write(c0,"(i8)") xStart
!!$    write(c1,"(i8)") xEnd
!!$    write(c2,"(i8)") yStart
!!$    write(c3,"(i8)") yEnd
!!$
!!$    select case (vTabPtr%idim_type)
!!$    case(2)
!!$       string="("//&
!!$            trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
!!$            trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
!!$    case(3)
!!$       write(c4,"(i8)") size(vTabPtr%var_p_3D,1)
!!$       string="(1:"//trim(adjustl(c4))//","//&
!!$            trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
!!$            trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
!!$    case(4:5)
!!$       write(c4,"(i8)") size(vTabPtr%var_p_4D,1)
!!$       write(c5,"(i8)") size(vTabPtr%var_p_4D,4)
!!$       string="(1:"//trim(adjustl(c4))//","//&
!!$            trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
!!$            trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
!!$            "1:"//trim(adjustl(c5))//")"
!!$
!!$    case(6:7)
!!$       write(c4,"(i8)") size(vTabPtr%var_p_3D,3)
!!$       string="("//&
!!$            trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
!!$            trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
!!$            "1:"//trim(adjustl(c4))//")"
!!$
!!$    case default
!!$       write(c0,"(i8)") vTabPtr%idim_type
!!$       call fatal_error(h//" field section "//trim(vTabPtr%name)//&
!!$            " with unknown idim_type="//trim(adjustl(c0)))
!!$    end select
!!$  end subroutine StringIndexing



  subroutine MarkLiteVarsAtVarTable(oneVarTable, oneVarTableSize)

    ! Arguments:
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(in) :: oneVarTableSize

    ! Local variables:
    integer :: nv
    integer :: nvl
    logical :: found

    ! Loop over each variable input in namelist "LITE_VARS" and set
    !   lite flag in ModVarTable

    if (nlite_vars > 0) then
       found=.false.
       do nvl=1,nlite_vars
          do nv = 1, oneVarTableSize
             ! find each lite var and set flag
             if (oneVarTable(nv)%name == lite_vars(nvl) ) then
                oneVarTable(nv)%ilite = 1
                found=.true.
                print*,'!---------------------------------------------------------'
                print*,'! LITE_VARS variable added--->',trim(lite_vars(nvl))
                print*,'!---------------------------------------------------------'
                exit
             end if
          end do
          if(.not. found) then
             print*,'!---------------------------------------------------------'
             print*,'! LITE_VARS variable does not exist in var table'
             print*,'!    variable name-->',lite_vars(nvl),'<--'
             print*,'!---------------------------------------------------------'
          end if
       end do
    end if
    
  end subroutine MarkLiteVarsAtVarTable


!!$  subroutine DumpVTab(ngrd)
!!$    integer, intent(in) :: ngrd
!!$
!!$    integer :: i
!!$    character(len=*), parameter :: h="**(DumpVTab)**"
!!$
!!$    write(*,"(a,i2)") h//" dump of vtab_r names for grid ",ngrd
!!$    write(*,"(a)") h//" name            idim_type           npts"
!!$    do i = 1, num_var(ngrd)
!!$       write(*,"(1x,a16,1x,i8,1x,i16)") vtab_r(i,ngrd)%name, &
!!$            vtab_r(i,ngrd)%idim_type, vtab_r(i,ngrd)%npts
!!$    end do
!!$  end subroutine DumpVTab




  subroutine FixVarTableForIOUTPUT5(oneVarTable, oneVarTableSize)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(in) :: oneVarTableSize

    integer :: ni
    integer :: nspc
    integer :: imode
    character(len=2) :: cmode
    character(len=*), parameter :: h="**(FixVarTableForIOUTPUT5)**"

    do ni = 1, oneVarTableSize

       oneVarTable(ni)%ianal=0	

       if (trim(oneVarTable(ni)%name) == 'TOPT'  .or. &
	    trim(oneVarTable(ni)%name) == 'UP'    .or. &
	    trim(oneVarTable(ni)%name) == 'VP'    .or. &
	    trim(oneVarTable(ni)%name) == 'THETA' .or. &
	    trim(oneVarTable(ni)%name) == 'PP'    .or. &
	    trim(oneVarTable(ni)%name) == 'RV') then
          oneVarTable(ni)%ianal=1
          cycle
       end if

       if (oneVarTable(ni)%irecycle == 1) then
          oneVarTable(ni)%ianal=1
       end if

    end do
  end subroutine FixVarTableForIOUTPUT5

  


  function Name2VarTableEntry(oneVarTable, oneVarTableSize, name) result(vtabPtr)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(in) :: oneVarTableSize
    character (len=*), intent(in) :: name
    type(VarTable), pointer   :: vtabPtr

    integer :: ni

    vtabPtr => null()
    do ni = 1, oneVarTableSize
       if (trim(oneVarTable(ni)%name) == trim(name)) then
          vtabPtr => oneVarTable(ni)
          exit
       end if
    end do
  end function Name2VarTableEntry




  
  subroutine InsertVarTable_2D_I(oneVarTable, oneVarTableSize, var, tabstr, varm, imean)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    integer, pointer, intent(in) :: var(:,:)
    character (len=*), intent(in) :: tabstr
    integer, pointer, intent(in) :: varm(:,:)
    integer, intent(in) :: imean

    integer(kind=int64) :: sizeVar
    integer(kind=int64) :: sizeVarm
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertVarTable_2D_I)**"

    ! field size

    sizeVar=int(size(var,1),int64) * &
         int(size(var,2),int64) 

    ! get new entry and fill VarTable token dependent integer components

    call NewVarTableEntry(oneVarTable, oneVarTableSize, sizeVar, imean, tabstr)

    ! imean and varm compatibility
    ! if imean==1, then varm should be allocated with the same size of var;
    ! if imean==0, then varm should not be allocated

    if (imean == 1) then
       if (.not. associated(varm)) then
          call fatal_error(h//" invoked for field "//&
               trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
               " to use unassociated average field")
       else
          sizeVarm=int(size(varm,1),int64) * &
               int(size(varm,2),int64)
          if (sizeVar /= sizeVarm) then
             write(str(1),"(i8)") size(var,1)
             write(str(2),"(i8)") size(var,2)
             write(str(5),"(i8)") size(varm,1)
             write(str(6),"(i8)") size(varm,2)
             call fatal_error(h//" invoked for field "//&
                  trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
                  " to use varm("//&
                  trim(adjustl(str(5)))//":"//&
                  trim(adjustl(str(6)))//")"//&
                  " with size incompatible to var("//&
                  trim(adjustl(str(1)))//":"//&
                  trim(adjustl(str(2)))//")")
          end if
       end if
    else if (imean == 0) then
       if (associated(varm)) then
          call fatal_error(h//" invoked for current field "//&
               trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
               " but average field is associated")
       end if
    else
       write(str(1),"(i8)") imean
       call fatal_error(h//" invoked for field "//&
            trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
            " with unknown imean="//trim(adjustl(str(1))))
    end if
    
    ! store full field

    oneVarTable(oneVarTableSize)%var_p_2D_I => var

    ! store average field if desired

    if (imean == 1) then
       oneVarTable(oneVarTableSize)%var_m_2D_I => varm
    end if
  end subroutine InsertVarTable_2D_I
  


  

  
  subroutine InsertVarTable_2D(oneVarTable, oneVarTableSize, var, tabstr, varm, imean)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    real, pointer, intent(in) :: var(:,:)
    character (len=*), intent(in) :: tabstr
    real, pointer, intent(in) :: varm(:,:)
    integer, intent(in) :: imean

    integer(kind=int64) :: sizeVar
    integer(kind=int64) :: sizeVarm
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertVarTable_2D)**"

    ! field size

    sizeVar=int(size(var,1),int64) * &
         int(size(var,2),int64) 

    ! get new entry and fill VarTable token dependent integer components

    call NewVarTableEntry(oneVarTable, oneVarTableSize, sizeVar, imean, tabstr)

    ! imean and varm compatibility
    ! if imean==1, then varm should be allocated with the same size of var;
    ! if imean==0, then varm should not be allocated

    if (imean == 1) then
       if (.not. associated(varm)) then
          call fatal_error(h//" invoked for field "//&
               trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
               " to use unassociated average field")
       else
          sizeVarm=int(size(varm,1),int64) * &
               int(size(varm,2),int64)
          if (sizeVar /= sizeVarm) then
             write(str(1),"(i8)") size(var,1)
             write(str(2),"(i8)") size(var,2)
             write(str(5),"(i8)") size(varm,1)
             write(str(6),"(i8)") size(varm,2)
             call fatal_error(h//" invoked for field "//&
                  trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
                  " to use varm("//&
                  trim(adjustl(str(5)))//":"//&
                  trim(adjustl(str(6)))//")"//&
                  " with size incompatible to var("//&
                  trim(adjustl(str(1)))//":"//&
                  trim(adjustl(str(2)))//")")
          end if
       end if
    else if (imean == 0) then
       if (associated(varm)) then
          call fatal_error(h//" invoked for current field "//&
               trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
               " but average field is associated")
       end if
    else
       write(str(1),"(i8)") imean
       call fatal_error(h//" invoked for field "//&
            trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
            " with unknown imean="//trim(adjustl(str(1))))
    end if
    
    ! store full field

    oneVarTable(oneVarTableSize)%var_p_2D => var

    ! store average field if desired

    if (imean == 1) then
       oneVarTable(oneVarTableSize)%var_m_2D => varm
    end if
  end subroutine InsertVarTable_2D
  
  

  subroutine InsertVarTable_3D(oneVarTable, oneVarTableSize, var, tabstr, varm, imean)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    real, pointer, intent(in) :: var(:,:,:)
    character (len=*), intent(in) :: tabstr
    real, pointer, intent(in) :: varm(:,:,:)
    integer, intent(in) :: imean

    integer(kind=int64) :: sizeVar
    integer(kind=int64) :: sizeVarm
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertVarTable_3D)**"

    ! field size

    sizeVar=int(size(var,1),int64) * &
         int(size(var,2),int64) * &
         int(size(var,3),int64)

    ! get new entry and fill VarTable token dependent integer components

    call NewVarTableEntry(oneVarTable, oneVarTableSize, sizeVar, imean, tabstr)

    ! imean and varm compatibility
    ! if imean==1, then varm should be allocated with the same size of var;
    ! if imean==0, then varm should not be allocated

    if (imean == 1) then
       if (.not. associated(varm)) then
          call fatal_error(h//" invoked for field "//&
               trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
               " to use unassociated average field")
       else
          sizeVarm=int(size(varm,1),int64) * &
               int(size(varm,2),int64) * &
               int(size(varm,3),int64)
          if (sizeVar /= sizeVarm) then
             write(str(1),"(i8)") size(var,1)
             write(str(2),"(i8)") size(var,2)
             write(str(3),"(i8)") size(var,3)
             write(str(4),"(i8)") size(varm,1)
             write(str(5),"(i8)") size(varm,2)
             write(str(6),"(i8)") size(varm,3)
             call fatal_error(h//" invoked for field "//&
                  trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
                  " to use varm("//&
                  trim(adjustl(str(4)))//":"//&
                  trim(adjustl(str(5)))//":"//&
                  trim(adjustl(str(6)))//")"//&
                  " with size incompatible to var("//&
                  trim(adjustl(str(1)))//":"//&
                  trim(adjustl(str(2)))//":"//&
                  trim(adjustl(str(3)))//")")
          end if
       end if
    else if (imean == 0) then
       if (associated(varm)) then
          call fatal_error(h//" invoked for current field "//&
               trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
               " but average field is associated")
       end if
    else
       write(str(1),"(i8)") imean
       call fatal_error(h//" invoked for field "//&
            trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
            " with unknown imean="//trim(adjustl(str(1))))
    end if
    
    ! store full field

    oneVarTable(oneVarTableSize)%var_p_3D => var

    ! store average field if desired

    if (imean == 1) then
       oneVarTable(oneVarTableSize)%var_m_3D => varm
    end if
  end subroutine InsertVarTable_3D
  


  subroutine InsertVarTable_4D(oneVarTable, oneVarTableSize, var, tabstr, varm, imean)
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    real, pointer, intent(in) :: var(:,:,:,:)
    character (len=*), intent(in) :: tabstr
    real, pointer, intent(in) :: varm(:,:,:,:)
    integer, intent(in) :: imean

    integer(kind=int64) :: sizeVar
    integer(kind=int64) :: sizeVarm
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertVarTable_4D)**"

    ! field size

    sizeVar=int(size(var,1),int64) * &
         int(size(var,2),int64) * &
         int(size(var,3),int64) * &
         int(size(var,4),int64)

    ! get new entry and fill VarTable token dependent integer components

    call NewVarTableEntry(oneVarTable, oneVarTableSize, sizeVar, imean, tabstr)

    ! imean and varm compatibility
    ! if imean==1, then varm should be allocated with the same size of var;
    ! if imean==0, then varm should not be allocated

    if (imean == 1) then
       if (.not. associated(varm)) then
          call fatal_error(h//" invoked for field "//&
               trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
               " to use unassociated average field")
       else
          sizeVarm=int(size(varm,1),int64) * &
               int(size(varm,2),int64) * &
               int(size(varm,3),int64) * &
               int(size(varm,4),int64)
          if (sizeVar /= sizeVarm) then
             write(str(1),"(i8)") size(var,1)
             write(str(2),"(i8)") size(var,2)
             write(str(3),"(i8)") size(var,3)
             write(str(4),"(i8)") size(var,4)
             write(str(5),"(i8)") size(varm,1)
             write(str(6),"(i8)") size(varm,2)
             write(str(7),"(i8)") size(varm,3)
             write(str(8),"(i8)") size(varm,4)
             call fatal_error(h//" invoked for field "//&
                  trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
                  " to use varm("//&
                  trim(adjustl(str(5)))//":"//&
                  trim(adjustl(str(6)))//":"//&
                  trim(adjustl(str(7)))//":"//&
                  trim(adjustl(str(8)))//")"//&
                  " with size incompatible to var("//&
                  trim(adjustl(str(1)))//":"//&
                  trim(adjustl(str(2)))//":"//&
                  trim(adjustl(str(3)))//":"//&
                  trim(adjustl(str(4)))//")")
          end if
       end if
    else if (imean == 0) then
       if (associated(varm)) then
          call fatal_error(h//" invoked for current field "//&
               trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
               " but average field is associated")
       end if
    else
       write(str(1),"(i8)") imean
       call fatal_error(h//" invoked for field "//&
            trim(adjustl(oneVarTable(oneVarTableSize)%name))//&
            " with unknown imean="//trim(adjustl(str(1))))
    end if
    
    ! store full field

    oneVarTable(oneVarTableSize)%var_p_4D => var

    ! store average field if desired

    if (imean == 1) then
       oneVarTable(oneVarTableSize)%var_m_4D => varm
    end if
  end subroutine InsertVarTable_4D
  
  
end module ModVarTable
