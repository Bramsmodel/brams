module ModFieldSection


  ! a list of field sections to be communicated
  ! among processes in a single message passing operation.
  ! List head and tail are stored elsewhere


  use ModParallelEnvironment, only: MsgDump

  implicit none

  private
  public :: FieldSection
  public :: SizeFieldSectionName
  public :: CreateFieldSection
  public :: StringFieldSection
  public :: DumpFieldSection
  public :: DestroyFieldSection
  public :: FieldSectionSize
  public :: FieldSectionName
  public :: FieldSectionData2Buffer
  public :: Buffer2FieldSectionData
  public :: UpdateFieldAdress

  integer, parameter :: SizeFieldSectionName=16

  type FieldSection
     ! one entry of a list of fields
     ! to be communicated to a single process
     ! in a single message passing operation.

     ! Data to communicate is the horizontal
     ! section [xStart:xEnd,yStart:yEnd] (local indices)
     ! of field_XXX (XXX=2D, 3D, 4D).

     ! If the field has more than 2 dimensions, then
     ! the remaining dimensions of each pair (x,y) of
     ! the section should be fully communicated.

     ! Component idim_type informs the remaining
     ! dimensions to be communicated, in a coded scheme.

     ! Component name has the field name to be communicated.

     ! Component fieldSectionSize is the size of the field
     ! to be communicated.

     real, pointer :: field_1D(:) => null()
     real, pointer :: field_2D(:,:) => null()
     real, pointer :: field_3D(:,:,:) => null()
     real, pointer :: field_4D(:,:,:,:) => null()
     ! field_XXX points to the array to extract
     ! the section to be communicated
     integer :: zStart = -1
     integer :: zEnd = -1
     integer :: xStart = -1
     integer :: xEnd = -1
     integer :: yStart = -1
     integer :: yEnd = -1
     ! the 2D section to be communicated is, in local indices,
     ! [xStart:xEnd,yStart:yEnd]
     ! zStart:zEnd simplifies coding for most 3D fields and
     ! accomodate sending a sub-section of the z dimension
     integer :: lbz = -1
     integer :: ubz = -1
     integer :: lbx = -1
     integer :: ubx = -1
     integer :: lby = -1
     integer :: uby = -1
     ! 1D fields are tipically a linearization of 3D arrays;
     ! as so, the bounds of the 3D array are unknown and should
     ! be stored
     integer :: idim_type = -1
     ! field dimensioning code, to know which other dimensions
     ! should be communicated:
     ! idim_type == 1 are linearized 3D fields,
     !   mapping 3D (lbz:ubz, lbx:ubx, lby:uby) into 1D
     ! idim_type == 2 means (nmxp, nmyp)
     !   no other dimensions to communicate
     ! idim_type == 3 means (nmzp, nmxp, nmyp)
     !   communicate first dimension for each (x,y)
     ! idim_type == 4 means (nzg, nmxp, nmyp, npatch)
     !   communicate first and last dimension for each (x,y)
     ! idim_type == 5 means (nzs, nmxp, nmyp, npatch)
     !   communicate first and last dimension for each (x,y)
     ! idim_type == 6 means (nmxp, nmyp, npatch)
     !   communicate  last dimension for each (x,y)
     ! idim_type == 7 means (nmxp, nmyp, nwave)
     !   communicate  last dimension for each (x,y)
     integer :: fieldSectionSize = -1
     ! number of data elements to communicate
     character(len=SizeFieldSectionName) :: name = ""
     ! field variable name
  end type FieldSection

  interface CreateFieldSection
     module procedure CreateFieldSection_1D
     module procedure CreateFieldSection_2D
     module procedure CreateFieldSection_3D
     module procedure CreateFieldSection_4D
     module procedure CreateFieldSection_Null
  end interface CreateFieldSection


  interface FieldSectionData2Buffer
     module procedure FieldSectionData2BufferFixedAdress
     module procedure FieldSectionData2BufferFixedAdress1D
  end interface FieldSectionData2Buffer


  interface Buffer2FieldSectionData
     module procedure Buffer2FieldSectionDataFixedAdress
     module procedure Buffer2FieldSectionDataFixedAdress1D
  end interface Buffer2FieldSectionData


  interface UpdateFieldAdress
     module procedure UpdateFieldAdress_0D
     module procedure UpdateFieldAdress_1D
     module procedure UpdateFieldAdress_2D
     module procedure UpdateFieldAdress_3D
     module procedure UpdateFieldAdress_4D
  end interface UpdateFieldAdress
contains





  function CreateFieldSection_1D(field, name, idim_type, &
       lbz, ubz, lbx, ubx, lby, uby, &
       zStart, zEnd, xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 1D field should be communicated

    real, pointer, intent(in) :: field(:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: lbz
    integer, intent(in) :: ubz
    integer, intent(in) :: lbx
    integer, intent(in) :: ubx
    integer, intent(in) :: lby
    integer, intent(in) :: uby
    integer, intent(in) :: zStart ! local index
    integer, intent(in) :: zEnd   ! local index
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateFieldSection_1D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(str(1))))
    end if
    oneFieldSection%field_1D => field
    oneFieldSection%lbz = lbz
    oneFieldSection%ubz = ubz
    oneFieldSection%lbx = lbx
    oneFieldSection%ubx = ubx
    oneFieldSection%lby = lby
    oneFieldSection%uby = uby
    oneFieldSection%zStart = zStart
    oneFieldSection%zEnd = zEnd
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type = idim_type
    if (idim_type == 1) then
       oneFieldSection%fieldSectionSize = &
            (zEnd - zStart +1) * &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1)
    else
       write(str(1),"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(str(1)))//&
            "; 1D fields are linearized 3D fields with idim_type=1")
    end if
    if (dumpLocal) then
       call DumpFieldSection(oneFieldSection, h//" created ")
    end if
  end function CreateFieldSection_1D





  function CreateFieldSection_2D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 2D field should be communicated

    real, pointer, intent(in) :: field(:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateFieldSection_2D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(str(1))))
    end if
    oneFieldSection%field_2d => field
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type = idim_type
    if (idim_type == 2) then
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1)
    else
       write(str(1),"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(str(1))))
    end if
    if (dumpLocal) then
       call DumpFieldSection(oneFieldSection, h//" created ")
    end if
  end function CreateFieldSection_2D





  function CreateFieldSection_3D(field, name, idim_type, &
       zStart, zEnd, xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 3D field should be communicated

    real, pointer, intent(in) :: field(:,:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: zStart ! local index
    integer, intent(in) :: zEnd   ! local index
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateFieldSection_3D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(str(1))))
    end if
    if (.not. associated(field)) then
       call fatal_error(h//" field "//trim(adjustl(name))//" not associated")
    end if
    oneFieldSection%field_3d => field
    oneFieldSection%lbz = lbound(field,1)
    oneFieldSection%ubz = ubound(field,1)
    oneFieldSection%lbx = lbound(field,2)
    oneFieldSection%ubx = ubound(field,2)
    oneFieldSection%lby = lbound(field,3)
    oneFieldSection%uby = ubound(field,3)
    oneFieldSection%zStart = lbound(field,1)
    oneFieldSection%zEnd = ubound(field,1)
    oneFieldSection%xStart = zStart
    oneFieldSection%xEnd = zEnd
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type=idim_type
    select case (idim_type)
    case (3)
       oneFieldSection%fieldSectionSize = &
            (zEnd - zStart +1) * &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1)
    case (6,7)
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,3)
    case default
       write(str(1),"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(str(1))))
    end select
    if (dumpLocal) then
       call DumpFieldSection(oneFieldSection, h//" created ")
    end if
  end function CreateFieldSection_3D





  function CreateFieldSection_4D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 4D field should be communicated

    real, pointer, intent(in) :: field(:,:,:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateFieldSection_4D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(str(1))))
    end if
    oneFieldSection%field_4d => field
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type=idim_type
    select case (idim_type)
    case (4, 5)
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,1) * &
            size(field,4)
    case default
       write(str(1),"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(str(1))))
    end select
    if (dumpLocal) then
       call DumpFieldSection(oneFieldSection, h//" created ")
    end if
  end function CreateFieldSection_4D




  function CreateFieldSection_Null(name, idim_type, &
       lbz, ubz, lbx, ubx, lby, uby, &
       zStart, zEnd, xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 4D field should be communicated

    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: lbz
    integer, intent(in) :: ubz
    integer, intent(in) :: lbx
    integer, intent(in) :: ubx
    integer, intent(in) :: lby
    integer, intent(in) :: uby
    integer, intent(in) :: zStart ! local index
    integer, intent(in) :: zEnd   ! local index
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateFieldSection_Null)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(str(1))))
    end if
    if (idim_type > 3) then
       write(str(1),"(i8)") idim_type
       call fatal_error(h//" not prepared for idim_type="//trim(adjustl(str(1)))//&
            " for field "//trim(adjustl(name)))
    end if
    if (idim_type == 2) then
       oneFieldSection%lbz = 1
       oneFieldSection%ubz = 1
    else
       oneFieldSection%lbz = lbz
       oneFieldSection%ubz = ubz
    end if
    oneFieldSection%lbx = lbx
    oneFieldSection%ubx = ubx
    oneFieldSection%lby = lby
    oneFieldSection%uby = uby
    oneFieldSection%zStart = zStart
    oneFieldSection%zEnd = zEnd
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type=idim_type
    oneFieldSection%fieldSectionSize = &
         (zEnd-zStart+1) * &
         (xEnd-xStart+1) * &
         (yEnd-yStart+1) 
    if (dumpLocal) then
       call DumpFieldSection(oneFieldSection, h//" created ")
    end if
  end function CreateFieldSection_Null



  function StringFieldSection(oneFieldSection) result(res)

    ! String with the fields of a type FieldSection variable

    type(FieldSection), pointer :: oneFieldSection
    character(len=256) :: res

    character(len=128) :: strSec
    character(len=128) :: strDim
    character(len=8) :: str(20)
    character(len=*), parameter :: h="**(StringFieldSection)**"

    if (.not. associated(oneFieldSection)) then
       res = " null FieldSection"
    else
       write(str(1),"(i8)") oneFieldSection%xStart
       write(str(2),"(i8)") oneFieldSection%xEnd
       write(str(3),"(i8)") oneFieldSection%yStart
       write(str(4),"(i8)") oneFieldSection%yEnd
       select case (oneFieldSection%idim_type)
       case(1)
          if (oneFieldSection%zStart /= -1) then
             write(str(5),"(i8)") oneFieldSection%zStart
             write(str(6),"(i8)") oneFieldSection%zEnd
          else
             str(5)="?"
             str(6)="?"
          end if
          strSec="1D map of ("//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//","//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")"
          strDim="(?:?)"
       case(2)
          strSec="("//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")"
          if (associated(oneFieldSection%field_2D)) then
             write(str(11),"(i8)") lbound(oneFieldSection%field_2D,1)
             write(str(12),"(i8)") ubound(oneFieldSection%field_2D,1)
             write(str(13),"(i8)") lbound(oneFieldSection%field_2D,2)
             write(str(14),"(i8)") ubound(oneFieldSection%field_2D,2)
             strDim="("//&
                  trim(adjustl(str(11)))//":"//trim(adjustl(str(12)))//","//&
                  trim(adjustl(str(13)))//":"//trim(adjustl(str(14)))//")"
          else
             strDim="(?:?,?:?)"
          end if
       case(3)
          if (oneFieldSection%zStart /= -1) then
             write(str(5),"(i8)") oneFieldSection%zStart
             write(str(6),"(i8)") oneFieldSection%zEnd
          else
             str(5)="?"
             str(6)="?"
          end if
          strSec="("//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//","//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")"
          if (associated(oneFieldSection%field_3D)) then
             write(str(11),"(i8)") lbound(oneFieldSection%field_3D,1)
             write(str(12),"(i8)") ubound(oneFieldSection%field_3D,1)
             write(str(13),"(i8)") lbound(oneFieldSection%field_3D,2)
             write(str(14),"(i8)") ubound(oneFieldSection%field_3D,2)
             write(str(15),"(i8)") lbound(oneFieldSection%field_3D,3)
             write(str(16),"(i8)") ubound(oneFieldSection%field_3D,3)
             strDim="("//&
                  trim(adjustl(str(11)))//":"//trim(adjustl(str(12)))//","//&
                  trim(adjustl(str(13)))//":"//trim(adjustl(str(14)))//","//&
                  trim(adjustl(str(16)))//":"//trim(adjustl(str(16)))//")"
          else
             strDim="(?:?,?:?,?:?)"
          end if
       case(4:5)
          if (associated(oneFieldSection%field_4D)) then
             write(str(5),"(i8)") lbound(oneFieldSection%field_4D,1)
             write(str(6),"(i8)") ubound(oneFieldSection%field_4D,1)
             write(str(7),"(i8)") lbound(oneFieldSection%field_4D,4)
             write(str(8),"(i8)") ubound(oneFieldSection%field_4D,4)
          else
             str(5)="?"
             str(6)="?"
             str(7)="?"
             str(8)="?"
          end if
          strSec="("//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//","//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
               trim(adjustl(str(7)))//":"//trim(adjustl(str(8)))//")"
          if (associated(oneFieldSection%field_4D)) then
             write(str(11),"(i8)") lbound(oneFieldSection%field_4D,1)
             write(str(12),"(i8)") ubound(oneFieldSection%field_4D,1)
             write(str(13),"(i8)") lbound(oneFieldSection%field_4D,2)
             write(str(14),"(i8)") ubound(oneFieldSection%field_4D,2)
             write(str(15),"(i8)") lbound(oneFieldSection%field_4D,3)
             write(str(16),"(i8)") ubound(oneFieldSection%field_4D,3)
             write(str(17),"(i8)") lbound(oneFieldSection%field_4D,4)
             write(str(18),"(i8)") ubound(oneFieldSection%field_4D,4)
             strDim="("//&
                  trim(adjustl(str(11)))//":"//trim(adjustl(str(12)))//","//&
                  trim(adjustl(str(13)))//":"//trim(adjustl(str(14)))//","//&
                  trim(adjustl(str(16)))//":"//trim(adjustl(str(16)))//","//&
                  trim(adjustl(str(17)))//":"//trim(adjustl(str(18)))//")"
          else
             strDim="(?:?,?:?,?:?,?,?)"
          end if
       case(6:7)
          if (oneFieldSection%zStart /= -1) then
             write(str(5),"(i8)") oneFieldSection%zStart
             write(str(6),"(i8)") oneFieldSection%zEnd
          else
             str(5)="?"
             str(6)="?"
          end if
          strSec="("//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
               trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")"
          if (associated(oneFieldSection%field_3D)) then
             write(str(11),"(i8)") lbound(oneFieldSection%field_3D,1)
             write(str(12),"(i8)") ubound(oneFieldSection%field_3D,1)
             write(str(13),"(i8)") lbound(oneFieldSection%field_3D,2)
             write(str(14),"(i8)") ubound(oneFieldSection%field_3D,2)
             write(str(15),"(i8)") lbound(oneFieldSection%field_3D,3)
             write(str(16),"(i8)") ubound(oneFieldSection%field_3D,3)
             strDim="("//&
                  trim(adjustl(str(11)))//":"//trim(adjustl(str(12)))//","//&
                  trim(adjustl(str(13)))//":"//trim(adjustl(str(14)))//","//&
                  trim(adjustl(str(16)))//":"//trim(adjustl(str(16)))//")"
          else
             strDim="(?:?,?:?,?:?)"
          end if
       case default
          write(str(1),"(i8)") oneFieldSection%idim_type
          call fatal_error(h//" field section "//trim(oneFieldSection%name)//&
               " with unknown idim_type="//trim(adjustl(str(1))))
       end select
       write(str(1),"(i8)") oneFieldSection%FieldSectionSize
       res = "field section "//trim(oneFieldSection%name)//&
            trim(strSec)//" of size "//trim(adjustl(str(1)))//&
            " dim "//trim(adjustl(strDim))
    end if
  end function StringFieldSection





  subroutine DumpFieldSection(oneFieldSection, strMsg)

    ! Dumps a variable of type FieldSection with own header
    ! or with header "strMsg", if present

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    character(len=*), intent(in), optional :: strMsg

    character(len=*), parameter :: h="**(DumpFieldSection)**"

    if (present(strMsg)) then
       call MsgDump(trim(strMsg)//" "//&
            trim(adjustl(StringFieldSection(oneFieldSection))))
    else
       call MsgDump(h//" "//&
            trim(adjustl(StringFieldSection(oneFieldSection))))
    end if
  end subroutine DumpFieldSection





  subroutine DestroyFieldSection(oneFieldSection)

    ! reclaims memory area and returns null pointer

    type(FieldSection), pointer, intent(inout) :: oneFieldSection

    integer :: ierr
    character(len=SizeFieldSectionName) :: name
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyFieldSection)**"
    logical, parameter :: dumpLocal=.false.

    if (associated(oneFieldSection)) then
       name = oneFieldSection%name

       ! avoid deallocating pointed area,
       ! just nullify field pointers since
       ! some other pointer may be pointing
       ! to the same field
       oneFieldSection%field_2D => null()
       oneFieldSection%field_3D => null()
       oneFieldSection%field_4D => null()
       deallocate(oneFieldSection, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate fails with stat="//&
               trim(adjustl(str(1))))
       end if
       if (dumpLocal) then
          call MsgDump(h//" named "//trim(adjustl(name)))
       end if
    end if
    nullify(oneFieldSection)
  end subroutine DestroyFieldSection




  function FieldSectionSize(oneFieldSection) result(nEle)

    ! returns the number of elements of a field section

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    integer :: nEle

    nEle=0
    if (associated(oneFieldSection)) then
       nEle=oneFieldSection%fieldSectionSize
    end if
  end function FieldSectionSize





  function FieldSectionName(oneFieldSection) result(name)

    ! name of a FieldSection; empty if not associated

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    character(len=SizeFieldSectionName) :: name

    name=""
    if (associated(oneFieldSection)) then
       name=oneFieldSection%name
    end if
  end function FieldSectionName





  subroutine FieldSectionData2BufferFixedAdress(oneFieldSection, &
       buf, bufStart, bufSize)

    ! copy field section values into a 1D buffer in array
    ! element order

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    ! field values of the field section to copy from
    real, intent(inout) :: buf(:)
    ! buffer to copy to
    integer, intent(inout) :: bufStart
    ! copy starts at buf(bufStart)
    integer, intent(in) :: bufSize
    ! buf maximum size

    integer :: finalPos
    integer :: posBuf
    integer :: x
    integer :: y
    integer :: z
    integer :: k
    integer :: lMax
    integer :: kMax
    character(len=8) :: str(10)
    character(len=8) :: buf0, bufn
    character(len=64) :: preStr, midStr, posStr
    character(len=*), parameter :: h="**(FieldSectionData2BufferFixedAdress)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" null oneFieldSection")
    end if

    finalPos=bufStart+FieldSectionSize(oneFieldSection)-1
    if (finalPos > bufSize) then
       write(str(1),"(i8)") FieldSectionSize(oneFieldSection)
       write(str(2),"(i8)") bufStart
       write(str(3),"(i8)") bufSize
       call fatal_error(h//" field "//&
            trim(adjustl(oneFieldSection%name))//&
            " with size "//trim(adjustl(str(1)))//&
            " is larger than buffer ("//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")")
    end if

    if (dumpLocal) then
       write(buf0,"(i8)") bufStart
       write(bufn,"(i8)") finalPos
       preStr=""
       write(str(1),"(i8)") oneFieldSection%xStart
       write(str(2),"(i8)") oneFieldSection%xEnd
       midStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//"," 
       write(str(1),"(i8)") oneFieldSection%yStart
       write(str(2),"(i8)") oneFieldSection%yEnd
       midStr=trim(adjustl(midStr))//trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))
       posStr=""
    end if

    posBuf=bufStart
    select case (oneFieldSection%idim_type)
    case(1)
       call fatal_error(h//" unknown kMax;"//&
            " use FieldSectionData2BufferFixedAdress1D interface")
    case(2)
       if (.not. associated(oneFieldSection%field_2D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
       end if
       do y=oneFieldSection%yStart, oneFieldSection%yEnd
          do x=oneFieldSection%xStart, oneFieldSection%xEnd
             buf(posBuf)=oneFieldSection%field_2D(x,y)
             posBuf=posBuf+1
          end do
       end do
    case(3)
       if (.not. associated(oneFieldSection%field_3D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
       end if
       kMax = size(oneFieldSection%field_3D,1)
       do y=oneFieldSection%yStart, oneFieldSection%yEnd
          do x=oneFieldSection%xStart, oneFieldSection%xEnd
             do z=oneFieldSection%zStart, oneFieldSection%zEnd
                buf(posBuf)=oneFieldSection%field_3D(z,x,y)
                posBuf=posBuf+1
             end do
          end do
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") oneFieldSection%zStart
          write(str(2),"(i8)") oneFieldSection%zEnd
          preStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//"," 
       end if

    case(4,5)
       if (.not. associated(oneFieldSection%field_4D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
       end if
       lMax = size(oneFieldSection%field_4D,1)
       kMax = size(oneFieldSection%field_4D,4)
       do k=lBound(oneFieldSection%field_4D,4), &
            uBound(oneFieldSection%field_4D,4)
          do y=oneFieldSection%yStart, oneFieldSection%yEnd
             do x=oneFieldSection%xStart, oneFieldSection%xEnd
                buf(posBuf:posBuf+lMax-1)=oneFieldSection%field_4D(:,x,y,k)
                posBuf=posBuf+lMax
             end do
          end do
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") lbound(oneFieldSection%field_4D,1)
          write(str(2),"(i8)") ubound(oneFieldSection%field_4D,1)
          preStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//"," 
          write(str(1),"(i8)") lbound(oneFieldSection%field_4D,4)
          write(str(2),"(i8)") ubound(oneFieldSection%field_4D,4)
          posStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))
       end if

    case(6,7)
       if (.not. associated(oneFieldSection%field_3D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
       end if
       do k=lBound(oneFieldSection%field_3D,3), &
            uBound(oneFieldSection%field_3D,3)
          do y=oneFieldSection%yStart, oneFieldSection%yEnd
             do x=oneFieldSection%xStart, oneFieldSection%xEnd
                buf(posBuf)=oneFieldSection%field_3D(x,y,k)
                posBuf=posBuf+1
             end do
          end do
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") lbound(oneFieldSection%field_3D,3)
          write(str(2),"(i8)") ubound(oneFieldSection%field_3D,3)
          posStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))
       end if

    end select
    bufStart=posBuf

    if (dumpLocal) then
       call MsgDump(h//" buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//&
            ") <- "//trim(adjustl(oneFieldSection%name))//"("//&
            trim(adjustl(preStr))//trim(adjustl(midStr))//trim(adjustl(posStr))//")")
    end if
  end subroutine FieldSectionData2BufferFixedAdress





  subroutine FieldSectionData2BufferFixedAdress1D(oneFieldSection, &
       nzp, nxp, nyp, &
       buf, bufStart, bufSize)

    ! copy field section values of a 3D pointer array that was
    ! created as a linearized 1D pointer array
    ! into a 1D buffer in array element order

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    ! field values of the field section to copy from
    real, intent(inout) :: buf(:)
    ! buffer to copy to
    integer, intent(inout) :: bufStart
    ! copy starts at buf(bufStart)
    integer, intent(in) :: bufSize
    ! buf maximum size
    integer, intent(in) :: nzp
    integer, intent(in) :: nxp
    integer, intent(in) :: nyp
    ! first dimension of 3D field

    integer :: finalPos
    integer :: x
    integer :: y
    integer :: k
    integer :: lMax
    character(len=8) :: str(10)
    character(len=8) :: buf0, bufn
    character(len=64) :: preStr, midStr
    character(len=*), parameter :: h="**(FieldSectionData2BufferFixedAdress1D)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" null oneFieldSection")
    end if

    finalPos=bufStart+FieldSectionSize(oneFieldSection)-1
    if (finalPos > bufSize) then
       write(str(1),"(i8)") FieldSectionSize(oneFieldSection)
       write(str(2),"(i8)") bufStart
       write(str(3),"(i8)") bufSize
       call fatal_error(h//" field "//&
            trim(adjustl(oneFieldSection%name))//&
            " with size "//trim(adjustl(str(1)))//&
            " is larger than buffer ("//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")")
    end if

    if (dumpLocal) then
       write(buf0,"(i8)") bufStart
       write(bufn,"(i8)") finalPos
       write(str(1),"(i8)") nzp
       preStr="1:"//trim(adjustl(str(1)))//","
       write(str(1),"(i8)") oneFieldSection%xStart
       write(str(2),"(i8)") oneFieldSection%xEnd
       midStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//"," 
       write(str(1),"(i8)") oneFieldSection%yStart
       write(str(2),"(i8)") oneFieldSection%yEnd
       midStr=trim(adjustl(midStr))//trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))
    end if

    select case (oneFieldSection%idim_type)
    case(1)
       if (.not. associated(oneFieldSection%field_1D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
       end if
       call HiddenFieldSectionData2Buffer(&
            oneFieldSection%field_1D, &
            oneFieldSection%name, &
            oneFieldSection%idim_type, &
            oneFieldSection%xStart, oneFieldSection%xEnd, &
            oneFieldSection%yStart,oneFieldSection% yEnd, &
            nzp, nxp, nyp, &
            bufStart, bufSize, buf)

    case default
       call fatal_error(h//" wrong FieldSectionData2Buffer interface call for field "//&
            trim(adjustl(oneFieldSection%name))//&
            "; remove nzp from the call")
    end select

    if (dumpLocal) then
       call MsgDump(h//" buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//&
            ") <- "//trim(adjustl(oneFieldSection%name))//"("//&
            trim(adjustl(preStr))//trim(adjustl(midStr))//") mapped 1D")
    end if
  end subroutine FieldSectionData2BufferFixedAdress1D




  subroutine Buffer2FieldSectionDataFixedAdress(oneFieldSection, &
       buf, bufStart, bufSize)

    ! copy field section values from a 1D buffer to
    ! the field at field section in array element order

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    ! field values of the field section to copy to
    real, intent(in) :: buf(:)
    ! buffer to copy from
    integer, intent(inout) :: bufStart
    ! copy starts at buf(bufStart)
    integer, intent(in) :: bufSize
    ! buf maximum size

    integer :: finalPos
    integer :: posBuf
    integer :: x
    integer :: y
    integer :: z
    integer :: k
    integer :: lMax
    integer :: kMax
    character(len=8) :: str(20)
    character(len=64) :: strDim
    character(len=8) :: buf0, bufn
    character(len=64) :: preStr, midStr, posStr
    character(len=*), parameter :: h="**(Buffer2FieldSectionDataFixedAdress)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" null oneFieldSection")
    end if

    finalPos=bufStart+FieldSectionSize(oneFieldSection)-1
    if (finalPos > bufSize) then
       write(str(1),"(i8)") FieldSectionSize(oneFieldSection)
       write(str(2),"(i8)") bufStart
       write(str(3),"(i8)") bufSize
       call fatal_error(h//" field "//&
            trim(adjustl(oneFieldSection%name))//&
            " with size "//trim(adjustl(str(1)))//&
            " is larger than buffer ("//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")")
    end if

    if (dumpLocal) then
       write(buf0,"(i8)") bufStart
       write(bufn,"(i8)") finalPos
       preStr=""
       write(str(1),"(i8)") oneFieldSection%xStart
       write(str(2),"(i8)") oneFieldSection%xEnd
       midStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//"," 
       write(str(1),"(i8)") oneFieldSection%yStart
       write(str(2),"(i8)") oneFieldSection%yEnd
       midStr=trim(adjustl(midStr))//trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))
       posStr=""
    end if

    posBuf=bufStart
    select case (oneFieldSection%idim_type)
    case (1)
       call fatal_error(h//" unknown kmax;"//&
            " use Buffer2FieldSectionDataFixedAdress1D interface")
    case(2)
       if (.not. associated(oneFieldSection%field_2D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
       end if
       if (dumpLocal) then
          write(str(11),"(i8)") lbound(oneFieldSection%field_2D,1)
          write(str(12),"(i8)") ubound(oneFieldSection%field_2D,1)
          write(str(13),"(i8)") lbound(oneFieldSection%field_2D,2)
          write(str(14),"(i8)") ubound(oneFieldSection%field_2D,2)
          strDim="("//&
               trim(adjustl(str(11)))//":"//trim(adjustl(str(12)))//","//&
               trim(adjustl(str(13)))//":"//trim(adjustl(str(14)))//")"
       end if
       do y=oneFieldSection%yStart, oneFieldSection%yEnd
          do x=oneFieldSection%xStart, oneFieldSection%xEnd
             oneFieldSection%field_2D(x,y) = buf(posBuf)
             posBuf=posBuf+1
          end do
       end do
    case(3)
       if (.not. associated(oneFieldSection%field_3D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
       end if
       if (dumpLocal) then
          write(str(11),"(i8)") lbound(oneFieldSection%field_3D,1)
          write(str(12),"(i8)") ubound(oneFieldSection%field_3D,1)
          write(str(13),"(i8)") lbound(oneFieldSection%field_3D,2)
          write(str(14),"(i8)") ubound(oneFieldSection%field_3D,2)
          write(str(15),"(i8)") lbound(oneFieldSection%field_3D,3)
          write(str(16),"(i8)") ubound(oneFieldSection%field_3D,3)
          strDim="("//&
               trim(adjustl(str(11)))//":"//trim(adjustl(str(12)))//","//&
               trim(adjustl(str(13)))//":"//trim(adjustl(str(14)))//","//&
               trim(adjustl(str(15)))//":"//trim(adjustl(str(16)))//")"
       end if
       kMax = size(oneFieldSection%field_3D,1)
       do y=oneFieldSection%yStart, oneFieldSection%yEnd
          do x=oneFieldSection%xStart, oneFieldSection%xEnd
             do z=oneFieldSection%zStart, oneFieldSection%zEnd
                oneFieldSection%field_3D(z,x,y) = buf(posBuf)
                posBuf=posBuf+1
             end do
          end do
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") lbound(oneFieldSection%field_3D,1)
          write(str(2),"(i8)") ubound(oneFieldSection%field_3D,1)
          preStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//"," 
       end if

    case(4,5)
       if (.not. associated(oneFieldSection%field_4D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
       end if
       if (dumpLocal) then
          write(str(11),"(i8)") lbound(oneFieldSection%field_4D,1)
          write(str(12),"(i8)") ubound(oneFieldSection%field_4D,1)
          write(str(13),"(i8)") lbound(oneFieldSection%field_4D,2)
          write(str(14),"(i8)") ubound(oneFieldSection%field_4D,2)
          write(str(15),"(i8)") lbound(oneFieldSection%field_4D,3)
          write(str(16),"(i8)") ubound(oneFieldSection%field_4D,3)
          write(str(17),"(i8)") lbound(oneFieldSection%field_4D,4)
          write(str(18),"(i8)") ubound(oneFieldSection%field_4D,4)
          strDim="("//&
               trim(adjustl(str(11)))//":"//trim(adjustl(str(12)))//","//&
               trim(adjustl(str(13)))//":"//trim(adjustl(str(14)))//","//&
               trim(adjustl(str(15)))//":"//trim(adjustl(str(16)))//","//&
               trim(adjustl(str(17)))//":"//trim(adjustl(str(18)))//")"
       end if
       lMax = size(oneFieldSection%field_4D,1)
       kMax = size(oneFieldSection%field_4D,4)
       do k=lBound(oneFieldSection%field_4D,4), &
            uBound(oneFieldSection%field_4D,4)
          do y=oneFieldSection%yStart, oneFieldSection%yEnd
             do x=oneFieldSection%xStart, oneFieldSection%xEnd
                oneFieldSection%field_4D(:,x,y,k) = buf(posBuf:posBuf+lMax-1)
                posBuf=posBuf+lMax
             end do
          end do
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") lbound(oneFieldSection%field_4D,1)
          write(str(2),"(i8)") ubound(oneFieldSection%field_4D,1)
          preStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//"," 
          write(str(1),"(i8)") lbound(oneFieldSection%field_4D,4)
          write(str(2),"(i8)") ubound(oneFieldSection%field_4D,4)
          posStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))
       end if

    case(6,7)
       if (.not. associated(oneFieldSection%field_3D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
       end if
       if (dumpLocal) then
          write(str(11),"(i8)") lbound(oneFieldSection%field_3D,1)
          write(str(12),"(i8)") ubound(oneFieldSection%field_3D,1)
          write(str(13),"(i8)") lbound(oneFieldSection%field_3D,2)
          write(str(14),"(i8)") ubound(oneFieldSection%field_3D,2)
          write(str(15),"(i8)") lbound(oneFieldSection%field_3D,3)
          write(str(16),"(i8)") ubound(oneFieldSection%field_3D,3)
          strDim="("//&
               trim(adjustl(str(11)))//":"//trim(adjustl(str(12)))//","//&
               trim(adjustl(str(13)))//":"//trim(adjustl(str(14)))//","//&
               trim(adjustl(str(15)))//":"//trim(adjustl(str(16)))//")"
       end if
       do k=lBound(oneFieldSection%field_3D,3), &
            uBound(oneFieldSection%field_3D,3)
          do y=oneFieldSection%yStart, oneFieldSection%yEnd
             do x=oneFieldSection%xStart, oneFieldSection%xEnd
                oneFieldSection%field_3D(x,y,k) = buf(posBuf)
                posBuf=posBuf+1
             end do
          end do
       end do

       if (dumpLocal) then
          write(str(1),"(i8)") lbound(oneFieldSection%field_3D,3)
          write(str(2),"(i8)") ubound(oneFieldSection%field_3D,3)
          posStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))
       end if

    end select
    bufStart=posBuf

    if (dumpLocal) then
       call MsgDump(h//" for "//trim(adjustl(oneFieldSection%name))//&
            " dimensioned "//trim(adjustl(strDim)))
       call MsgDump(h//" "//trim(adjustl(oneFieldSection%name))//"("//&
            trim(adjustl(preStr))//trim(adjustl(midStr))//trim(adjustl(posStr))//&
            ") <- buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//")")

    end if
  end subroutine Buffer2FieldSectionDataFixedAdress





  subroutine Buffer2FieldSectionDataFixedAdress1D(oneFieldSection, &
       nzp, nxp, nyp, buf, bufStart, bufSize)

    ! copy field section values from a 1D buffer to
    ! the field at field section in array element order

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    ! field values of the field section to copy to
    real, intent(in) :: buf(:)
    ! buffer to copy from
    integer, intent(inout) :: bufStart
    ! copy starts at buf(bufStart)
    integer, intent(in) :: bufSize
    ! buf maximum size
    integer, intent(in) :: nzp
    integer, intent(in) :: nxp
    integer, intent(in) :: nyp

    integer :: finalPos
    integer :: posBuf
    integer :: x
    integer :: y
    integer :: k
    integer :: lMax
    character(len=8) :: str(20)
    character(len=64) :: strDim
    character(len=8) :: buf0, bufn
    character(len=64) :: preStr, midStr, posStr
    character(len=*), parameter :: h="**(Buffer2FieldSectionDataFixedAdress1D)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" null oneFieldSection")
    end if

    finalPos=bufStart+FieldSectionSize(oneFieldSection)-1
    if (finalPos > bufSize) then
       write(str(1),"(i8)") FieldSectionSize(oneFieldSection)
       write(str(2),"(i8)") bufStart
       write(str(3),"(i8)") bufSize
       call fatal_error(h//" field "//&
            trim(adjustl(oneFieldSection%name))//&
            " with size "//trim(adjustl(str(1)))//&
            " is larger than buffer ("//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//")")
    end if

    if (dumpLocal) then
       write(buf0,"(i8)") bufStart
       write(bufn,"(i8)") finalPos
       preStr=""
       write(str(1),"(i8)") oneFieldSection%xStart
       write(str(2),"(i8)") oneFieldSection%xEnd
       midStr=trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//"," 
       write(str(1),"(i8)") oneFieldSection%yStart
       write(str(2),"(i8)") oneFieldSection%yEnd
       midStr=trim(adjustl(midStr))//trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))
       posStr=""
    end if

    if (dumpLocal) then
       call MsgDump(h//" will do "//trim(adjustl(oneFieldSection%name))//"("//&
            trim(adjustl(preStr))//trim(adjustl(midStr))//trim(adjustl(posStr))//&
            ") <- buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//")")

    end if

    select case (oneFieldSection%idim_type)
    case (1)
       call HiddenBuffer2FieldSectionData(&
            oneFieldSection%field_1D, &
            oneFieldSection%name, &
            oneFieldSection%idim_type, &
            oneFieldSection%xStart, oneFieldSection%xEnd, &
            oneFieldSection%yStart, oneFieldSection%yEnd, &
            nzp, nxp, nyp, &
            bufStart, bufSize, buf)
    case default
       call fatal_error(h//" this interface requires 1D field;"//&
            " instead, use Buffer2FieldSectionDataFixedAdress interface")
    end select

    if (dumpLocal) then
       call MsgDump(h//" "//trim(adjustl(oneFieldSection%name))//"("//&
            trim(adjustl(preStr))//trim(adjustl(midStr))//trim(adjustl(posStr))//&
            ") <- buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//")")

    end if
  end subroutine Buffer2FieldSectionDataFixedAdress1D




  subroutine UpdateFieldAdress_0D(oneFieldSection, field, fieldName)
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    real, pointer, intent(in) :: field
    character(len=*), intent(in) :: fieldName

    call HiddenUpdateFieldAdress_0D(oneFieldSection, field, fieldName)
  end subroutine UpdateFieldAdress_0D



  subroutine UpdateFieldAdress_1D(oneFieldSection, field, fieldName)
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    real, pointer, intent(in) :: field(:)
    character(len=*), intent(in) :: fieldName

    oneFieldSection%field_1D => field
    oneFieldSection%name = fieldName
  end subroutine UpdateFieldAdress_1D



  subroutine UpdateFieldAdress_2D(oneFieldSection, field, fieldName)
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    real, pointer, intent(in) :: field(:,:)
    character(len=*), intent(in) :: fieldName

    oneFieldSection%field_2D => field
    oneFieldSection%name = fieldName
  end subroutine UpdateFieldAdress_2D



  subroutine UpdateFieldAdress_3D(oneFieldSection, field, fieldName)
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    real, pointer, intent(in) :: field(:,:,:)
    character(len=*), intent(in) :: fieldName

    oneFieldSection%field_3D => field
    oneFieldSection%name = fieldName
  end subroutine UpdateFieldAdress_3D



  subroutine UpdateFieldAdress_4D(oneFieldSection, field, fieldName)
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    real, pointer, intent(in) :: field(:,:,:,:)
    character(len=*), intent(in) :: fieldName

    oneFieldSection%field_4D => field
    oneFieldSection%name = fieldName
  end subroutine UpdateFieldAdress_4D
end module ModFieldSection





subroutine HiddenFieldSectionData2Buffer(field, &
     fieldName, idim_type, xStart, xEnd, yStart, yEnd, &
     nzp, nxp, nyp, &
     bufStart, bufSize, buf)

  ! hiden interface to allow actual 1D or 3D field argument,
  ! due to invocations of FieldSectionData2BufferVariableAdress
  ! passing scalar_tab pointers and also tendency arrays, besides
  ! the usual 3D field arrays

  use ModParallelEnvironment, only: MsgDump

  implicit none
  integer, intent(in) :: nzp
  integer, intent(in) :: nxp
  integer, intent(in) :: nyp
  real, intent(in) :: field(nzp,nxp,nyp)
  character(len=*), intent(in) :: fieldName
  integer, intent(in) :: idim_type
  integer, intent(in) :: xStart
  integer, intent(in) :: xEnd
  integer, intent(in) :: yStart
  integer, intent(in) :: yEnd
  integer, intent(inout) :: bufStart
  integer, intent(in) :: bufSize
  real, intent(inout) :: buf(bufSize)

  integer :: x, y, z
  integer :: posBuf

  character(len=8) :: buf0
  character(len=8) :: bufn
  character(len=8) :: str(10)
  character(len=16) :: strReal
  integer :: sectionSize
  integer :: bufEnd

  character(len=*), parameter :: h="**(HiddenFieldSectionData2Buffer)**"
  logical, parameter :: dumpLocal=.false.

  sectionSize=nzp*(xEnd-xStart+1)*(yEnd-yStart+1)
  bufEnd=bufStart+sectionSize-1
  if (dumpLocal) then
     write(buf0,"(i8)") bufStart
     write(str(1),"(i8)") bufSize
     write(str(1),"(i8)") bufEnd
     call MsgDump(h//" with bufSize="//trim(adjustl(str(1)))//&
          ", bufStart="//trim(adjustl(buf0))//&
          ", bufStart+sectionSize="//trim(adjustl(str(2))))
  end if

  posBuf=bufStart

  if (dumpLocal) then
     write(strReal,"(e15.7)") buf(bufStart)
     call MsgDump(h//" buf(bufStart)="//trim(adjustl(strReal)))
  end if

  if (dumpLocal) then
     write(strReal,"(e15.7)") buf(bufEnd)
     call MsgDump(h//" buf(bufEnd)="//trim(adjustl(strReal)))
  end if

  if (dumpLocal) then
     write(str(1),"(i8)") xStart
     write(str(2),"(i8)") yStart
     write(strReal,"(e15.7)") field(1,xStart,yStart)
     call MsgDump(h//" field(1,"//&
          trim(adjustl(str(1)))//","//&
          trim(adjustl(str(2)))//")="//&
          trim(adjustl(strReal)))
  end if

  if (dumpLocal) then
     write(str(1),"(i8)") xEnd
     write(str(2),"(i8)") yEnd
     write(strReal,"(e15.7)") field(nzp,xEnd,yEnd)
     call MsgDump(h//" field(nzp,"//&
          trim(adjustl(str(1)))//","//&
          trim(adjustl(str(2)))//")="//&
          trim(adjustl(strReal)))
  end if

  select case (idim_type)
  case(1,3)
     do y=yStart, yEnd
        do x=xStart, xEnd
           buf(posBuf:posBuf+nzp-1) = field(1:nzp,x,y)
           posBuf=posBuf+nzp
        end do
     end do

  case default
     write(str(1),"(i8)") idim_type
     call fatal_error(h//" idim_type("//&
          trim(adjustl(str(1)))//") /= 3")
  end select

  bufStart=posBuf

  if (dumpLocal) then
     write(bufn,"(i8)") bufStart-1
     write(str(1),"(i8)") xStart
     write(str(2),"(i8)") xEnd
     write(str(3),"(i8)") yStart
     write(str(4),"(i8)") yEnd
     write(str(5),"(i8)") 1
     write(str(6),"(i8)") nzp

     call MsgDump(h//" buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//&
          ") <- "//trim(adjustl(fieldName))//"("//&
          trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//","//&
          trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
          trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")")
  end if
end subroutine HiddenFieldSectionData2Buffer





subroutine HiddenBuffer2FieldSectionData(field, &
     fieldName, idim_type, xStart, xEnd, yStart, yEnd, &
     nzp, nxp, nyp, &
     bufStart, bufSize, buf)

  ! hiden interface to allow actual 1D field argument to be
  ! treated as a 3D field argument due to invocations of
  ! Buffer2FieldSectionData passing tend% pointers to 1D arrays

  use ModParallelEnvironment, only: MsgDump

  implicit none
  integer, intent(in) :: nzp
  integer, intent(in) :: nxp
  integer, intent(in) :: nyp
  real, intent(inout) :: field(nzp,nxp,nyp)
  character(len=*), intent(in) :: fieldName
  integer, intent(in) :: idim_type
  integer, intent(in) :: xStart
  integer, intent(in) :: xEnd
  integer, intent(in) :: yStart
  integer, intent(in) :: yEnd
  integer, intent(inout) :: bufStart
  integer, intent(in) :: bufSize
  real, intent(in) :: buf(bufSize)

  integer :: x, y
  integer :: posBuf

  character(len=8) :: buf0
  character(len=8) :: bufn
  character(len=8) :: str(20)
  character(len=64) :: strDim

  character(len=*), parameter :: h="**(HiddenBuffer2FieldSectionData)**"
  logical, parameter :: dumpLocal=.false.

  if (dumpLocal) then
     write(buf0,"(i8)") bufStart
     write(str(1),"(i8)") bufSize
     call MsgDump(h//" with bufSize="//trim(adjustl(str(1))))
  end if

  posBuf=bufStart

  select case (idim_type)
  case(1)
     if (dumpLocal) then
        write(str(11),"(i8)") lbound(field,1)
        write(str(12),"(i8)") ubound(field,1)
        write(str(13),"(i8)") lbound(field,2)
        write(str(14),"(i8)") ubound(field,2)
        write(str(15),"(i8)") lbound(field,3)
        write(str(16),"(i8)") ubound(field,3)
        strDim="("//&
             trim(adjustl(str(11)))//":"//trim(adjustl(str(12)))//","//&
             trim(adjustl(str(13)))//":"//trim(adjustl(str(14)))//","//&
             trim(adjustl(str(15)))//":"//trim(adjustl(str(16)))//")"
     end if
     do y=yStart, yEnd
        do x=xStart, xEnd
           field(1:nzp,x,y) = buf(posBuf:posBuf+nzp-1)
           posBuf=posBuf+nzp
        end do
     end do

  case default
     write(str(1),"(i8)") idim_type
     call fatal_error(h//" idim_type("//&
          trim(adjustl(str(1)))//") /= 1")
  end select

  bufStart=posBuf

  if (dumpLocal) then
     write(bufn,"(i8)") bufStart-1
     write(str(1),"(i8)") xStart
     write(str(2),"(i8)") xEnd
     write(str(3),"(i8)") yStart
     write(str(4),"(i8)") yEnd
     write(str(5),"(i8)") 1
     write(str(6),"(i8)") nzp
     call MsgDump(h//" for "//trim(adjustl(fieldName))//&
          " declared "//trim(adjustl(strDim)))
     call MsgDump(h//trim(adjustl(fieldName))//"("//&
          trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//","//&
          trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
          trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//") <- "//&
          " buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn)))
  end if
end subroutine HiddenBuffer2FieldSectionData

subroutine HiddenUpdateFieldAdress_0D(oneFieldSection, field, fieldName)
  use ModFieldSection, only: &
       FieldSection
  implicit none
  type(FieldSection), pointer, intent(in) :: oneFieldSection
  real, target, intent(in) :: field(:,:,:)
  character(len=*), intent(in) :: fieldName

  oneFieldSection%field_3D => field
  oneFieldSection%name = fieldName
end subroutine HiddenUpdateFieldAdress_0D
