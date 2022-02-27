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
  public :: UpdateFieldAdress
  public :: FieldSectionSize
  public :: FieldSectionName
  public :: FieldSectionData2Buffer
  public :: FieldSectionData2BufferVariableAdressArr
  public :: FieldSectionData2BufferVariableAdressScalar
  public :: Buffer2FieldSectionData

  integer, parameter :: SizeFieldSectionName=16

  type FieldSection
     private
     ! one entry of a list of fields
     ! to be communicated to a single process
     ! in a single message passing operation.

     ! Data to communicate is the horizontal
     ! section [xStart:xEnd,yStart:yEnd] (local indices)
     ! of field_XXX (XXX=2D, 3D, 4D or I2D).

     ! If the field has more than 2 dimensions, then
     ! the remaining dimensions of each pair (x,y) of
     ! the section should be fully communicated.

     ! Component idim_type informs the remaining
     ! dimensions to be communicated, in a coded scheme.

     ! Component name has the field name to be communicated.

     ! Component fieldSectionSize is the size of the field
     ! to be communicated.

     real, pointer :: field_2D(:,:) => null()
     real, pointer :: field_3D(:,:,:) => null()
     real, pointer :: field_4D(:,:,:,:) => null()
     integer, pointer :: field_I2D(:,:) => null()
     ! field_XXX points to the array to extract
     ! the section to be communicated
     integer :: xStart = -1
     integer :: xEnd = -1
     integer :: yStart = -1
     integer :: yEnd = -1
     ! the 2D section to be communicated is, in local indices,
     ! [xStart:xEnd,yStart:yEnd]
     integer :: idim_type = -1
     ! field dimensioning code, to know which other dimensions
     ! should be communicated:
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
     module procedure CreateFieldSection_I2D
     module procedure CreateFieldSection_2D
     module procedure CreateFieldSection_3D
     module procedure CreateFieldSection_4D
     module procedure CreateFieldSection_Null
  end interface CreateFieldSection


  interface UpdateFieldAdress
     module procedure UpdateFieldAdress_2D
     module procedure UpdateFieldAdress_3D
     module procedure UpdateFieldAdress_4D
  end interface UpdateFieldAdress

  interface FieldSectionData2Buffer
     module procedure FieldSectionData2BufferFixedAdress
  end interface FieldSectionData2Buffer

  interface Buffer2FieldSectionData
     module procedure Buffer2FieldSectionDataFixedAdress
     module procedure Buffer2FieldSectionDataVariableAdress
  end interface Buffer2FieldSectionData
contains





  function CreateFieldSection_I2D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the integer 2D field should be communicated

    integer, pointer, intent(in) :: field(:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_I2D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(c0)))
    end if
    oneFieldSection%field_I2d => field
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
       write(c0,"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(c0)))
    end if
    if (dumpLocal) then
       call DumpFieldSection(oneFieldSection, h//" created ")
    end if
  end function CreateFieldSection_I2D





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
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_2D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(c0)))
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
       write(c0,"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(c0)))
    end if
    if (dumpLocal) then
       call DumpFieldSection(oneFieldSection, h//" created ")
    end if
  end function CreateFieldSection_2D





  function CreateFieldSection_3D(field, name, idim_type, &
       xStart, xEnd, yStart, yEnd) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 3D field should be communicated

    real, pointer, intent(in) :: field(:,:,:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_3D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(c0)))
    end if
    oneFieldSection%field_3d => field
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type=idim_type
    select case (idim_type)
    case (3)
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,1)
    case (6,7)
       oneFieldSection%fieldSectionSize = &
            (yEnd - yStart +1) * &
            (xEnd - xStart +1) * &
            size(field,3)
    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(c0)))
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
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_4D)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(c0)))
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
       write(c0,"(i8)") idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(c0)))
    end select
    if (dumpLocal) then
       call DumpFieldSection(oneFieldSection, h//" created ")
    end if
  end function CreateFieldSection_4D




  function CreateFieldSection_Null(name, idim_type, &
       xStart, xEnd, yStart, yEnd, fieldSectionSize) result(oneFieldSection)

    ! stores at oneFieldSection which elements of
    ! the real 4D field should be communicated

    character(len=*), intent(in) :: name
    integer, intent(in) :: idim_type
    integer, intent(in) :: xStart ! local index
    integer, intent(in) :: xEnd   ! local index
    integer, intent(in) :: yStart ! local index
    integer, intent(in) :: yEnd   ! local index
    integer, intent(in) :: fieldSectionSize   ! acount all dimensions
    type(FieldSection), pointer :: oneFieldSection

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreateFieldSection_Null)**"
    logical, parameter :: dumpLocal=.false.

    allocate(oneFieldSection, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate(oneFieldSection) fails with stat="//&
            trim(adjustl(c0)))
    end if
    oneFieldSection%xStart = xStart
    oneFieldSection%xEnd = xEnd
    oneFieldSection%yStart = yStart
    oneFieldSection%yEnd = yEnd
    oneFieldSection%name = name
    oneFieldSection%idim_type=idim_type
    oneFieldSection%fieldSectionSize = fieldSectionSize
    if (dumpLocal) then
       call DumpFieldSection(oneFieldSection, h//" created ")
    end if
  end function CreateFieldSection_Null



  function StringFieldSection(oneFieldSection) result(res)

    ! String with the fields of a type FieldSection variable

    type(FieldSection), pointer :: oneFieldSection
    character(len=256) :: res

    character(len=128) :: string
    character(len=8) :: c0, c1, c2, c3, c4, c5, c6, c7
    character(len=*), parameter :: h="**(StringFieldSection)**"

    if (.not. associated(oneFieldSection)) then
       res = " null FieldSection"
    else
       write(c0,"(i8)") oneFieldSection%xStart
       write(c1,"(i8)") oneFieldSection%xEnd
       write(c2,"(i8)") oneFieldSection%yStart
       write(c3,"(i8)") oneFieldSection%yEnd
       select case (oneFieldSection%idim_type)
       case(2)
          string="("//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
       case(3)
          if (associated(oneFieldSection%field_3D)) then
             write(c4,"(i8)") lbound(oneFieldSection%field_3D,1)
             write(c5,"(i8)") ubound(oneFieldSection%field_3D,1)
          else
             c4="Unknown"
             c5="Unknown"
          end if
          string="("//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//","//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
       case(4:5)
          if (associated(oneFieldSection%field_4D)) then
             write(c4,"(i8)") lbound(oneFieldSection%field_4D,1)
             write(c5,"(i8)") ubound(oneFieldSection%field_4D,1)
             write(c6,"(i8)") lbound(oneFieldSection%field_4D,4)
             write(c7,"(i8)") ubound(oneFieldSection%field_4D,4)
          else
             c4="Unknown"
             c5="Unknown"
             c6="Unknown"
             c7="Unknown"
          end if
          string="("//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//","//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c6))//":"//trim(adjustl(c7))//")"
       case(6:7)
          if (associated(oneFieldSection%field_3D)) then
             write(c4,"(i8)") lbound(oneFieldSection%field_3D,3)
             write(c5,"(i8)") ubound(oneFieldSection%field_3D,3)
          else
             c4="Unknown"
             c5="Unknown"
          end if
          string="("//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//")"
       case default
          write(c0,"(i8)") oneFieldSection%idim_type
          call fatal_error(h//" field section "//trim(oneFieldSection%name)//&
               " with unknown idim_type="//trim(adjustl(c0)))
       end select
       write(c0,"(i8)") oneFieldSection%FieldSectionSize
       res = "field section "//trim(oneFieldSection%name)//&
            trim(string)//" of size "//trim(adjustl(c0))
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
    character(len=8) :: c0
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
       oneFieldSection%field_I2D => null()
       deallocate(oneFieldSection, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate fails with stat="//&
               trim(adjustl(c0)))
       end if
       if (dumpLocal) then
          call MsgDump(h//" named "//trim(adjustl(name)))
       end if
    end if
    nullify(oneFieldSection)
  end subroutine DestroyFieldSection




  subroutine UpdateFieldAdress_2D(oneFieldSection, field)
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    real, pointer, intent(in) :: field(:,:)

    character(len=8) :: c0
    character(len=*), parameter :: h="**(UpdateFieldAdress_2D)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" oneFieldSection not associated")
    else if (oneFieldSection%idim_type /= 2) then
       write(c0,"(i8)") oneFieldSection%idim_type
       call fatal_error(h//" idim_type ("//trim(adjustl(c0))//&
            " incompatible with field_2D")
    end if
    oneFieldSection%field_2D => field
    if (dumpLocal) then
       call MsgDump(h//" successfully updated "//&
            trim(adjustl(oneFieldSection%name))//&
            " memory address")
    end if
  end subroutine UpdateFieldAdress_2D




  subroutine UpdateFieldAdress_3D(oneFieldSection, field)
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    real, pointer, intent(in) :: field(:,:,:)

    character(len=8) :: c0
    character(len=*), parameter :: h="**(UpdateFieldAdress_3D)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" oneFieldSection not associated")
    else if (&
         (oneFieldSection%idim_type /= 3) .and. &
         (oneFieldSection%idim_type /= 6) .and. &
         (oneFieldSection%idim_type /= 7) ) then
       write(c0,"(i8)") oneFieldSection%idim_type
       call fatal_error(h//" idim_type ("//trim(adjustl(c0))//&
            " incompatible with field_3D")
    end if
    oneFieldSection%field_3D => field
    if (dumpLocal) then
       call MsgDump(h//" successfully updated "//&
            trim(adjustl(oneFieldSection%name))//&
            " memory address")
    end if
  end subroutine UpdateFieldAdress_3D




  subroutine UpdateFieldAdress_4D(oneFieldSection, field)
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    real, pointer, intent(in) :: field(:,:,:,:)

    character(len=8) :: c0
    character(len=*), parameter :: h="**(UpdateFieldAdress_4D)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" oneFieldSection not associated")
    else if (&
         (oneFieldSection%idim_type /= 4) .and. &
         (oneFieldSection%idim_type /= 5) ) then
       write(c0,"(i8)") oneFieldSection%idim_type
       call fatal_error(h//" idim_type ("//trim(adjustl(c0))//&
            " incompatible with field_4D")
    end if
    oneFieldSection%field_4D => field
    if (dumpLocal) then
       call MsgDump(h//" successfully updated "//&
            trim(adjustl(oneFieldSection%name))//&
            " memory address")
    end if
  end subroutine UpdateFieldAdress_4D




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
    integer :: k
    integer :: lMax
    integer :: kMax
    character(len=8) :: c0, c1, c2
    character(len=8) :: buf0, bufn
    character(len=64) :: preStr, midStr, posStr
    character(len=*), parameter :: h="**(FieldSectionData2BufferFixedAdress)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" null oneFieldSection")
    end if

    finalPos=bufStart+FieldSectionSize(oneFieldSection)-1
    if (finalPos > bufSize) then
       write(c0,"(i8)") FieldSectionSize(oneFieldSection)
       write(c1,"(i8)") bufStart
       write(c2,"(i8)") bufSize
       call fatal_error(h//" field "//&
            trim(adjustl(oneFieldSection%name))//&
            " with size "//trim(adjustl(c0))//&
            " is larger than buffer ("//&
            trim(adjustl(c1))//":"//trim(adjustl(c2))//")")
    end if

    if (dumpLocal) then
       write(buf0,"(i8)") bufStart
       write(bufn,"(i8)") finalPos
       preStr=""
       write(c0,"(i8)") oneFieldSection%xStart
       write(c1,"(i8)") oneFieldSection%xEnd
       midStr=trim(adjustl(c0))//":"//trim(adjustl(c1))//"," 
       write(c0,"(i8)") oneFieldSection%yStart
       write(c1,"(i8)") oneFieldSection%yEnd
       midStr=trim(adjustl(midStr))//trim(adjustl(c0))//":"//trim(adjustl(c1))
       posStr=""
    end if

    posBuf=bufStart
    select case (oneFieldSection%idim_type)
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
             buf(posBuf:posBuf+kMax-1)=oneFieldSection%field_3D(:,x,y)
             posBuf=posBuf+kMax
          end do
       end do

       if (dumpLocal) then
          write(c0,"(i8)") lbound(oneFieldSection%field_3D,1)
          write(c1,"(i8)") ubound(oneFieldSection%field_3D,1)
          preStr=trim(adjustl(c0))//":"//trim(adjustl(c1))//"," 
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
          write(c0,"(i8)") lbound(oneFieldSection%field_4D,1)
          write(c1,"(i8)") ubound(oneFieldSection%field_4D,1)
          preStr=trim(adjustl(c0))//":"//trim(adjustl(c1))//"," 
          write(c0,"(i8)") lbound(oneFieldSection%field_4D,4)
          write(c1,"(i8)") ubound(oneFieldSection%field_4D,4)
          posStr=trim(adjustl(c0))//":"//trim(adjustl(c1))
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
          write(c0,"(i8)") lbound(oneFieldSection%field_3D,3)
          write(c1,"(i8)") ubound(oneFieldSection%field_3D,3)
          posStr=trim(adjustl(c0))//":"//trim(adjustl(c1))
       end if

    end select
    bufStart=posBuf

    if (dumpLocal) then
       call MsgDump(h//" buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//&
            ") <- "//trim(adjustl(oneFieldSection%name))//"("//&
            trim(adjustl(preStr))//trim(adjustl(midStr))//trim(adjustl(posStr))//")")
    end if
  end subroutine FieldSectionData2BufferFixedAdress






  subroutine FieldSectionData2BufferVariableAdressArr(field, kMax, oneFieldSection, &
       buf, bufStart, bufSize)

    ! copy field section values into a 1D buffer in array
    ! element order, for a 3D field with variable memory adress

    real, target, intent(in) :: field(:,:,:)
    integer, intent(in) :: kMax
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
    real, pointer :: fieldPtr(:,:,:)

    character(len=8) :: str(10)
    character(len=8) :: buf0, bufn
    character(len=*), parameter :: h="**(FieldSectionData2BufferVariableAdressArr)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" null oneFieldSection")
    end if

    fieldPtr => field
    call UpdateFieldAdress_3D(oneFieldSection, fieldPtr)
    
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
       write(str(1),"(i8)") bufSize
       call MsgDump(h//" with bufSize="//trim(adjustl(str(1))))
    end if

    posBuf = bufStart

    select case (oneFieldSection%idim_type)
    case(3)
       do y=oneFieldSection%yStart, oneFieldSection%yEnd
          do x=oneFieldSection%xStart, oneFieldSection%xEnd
             do z = 1, kMax
                buf(z-1+posBuf) = field(z,x,y)
             end do
             posBuf=posBuf+kMax
          end do
       end do

    case default
       write(str(1),"(i8)") oneFieldSection%idim_type
       call fatal_error(h//" oneFieldSection%idim_type("//&
            trim(adjustl(str(1)))//") /= 3")
    end select

    bufStart=posBuf

    if (dumpLocal) then
       write(bufn,"(i8)") bufStart-1
       write(str(1),"(i8)") oneFieldSection%xStart
       write(str(2),"(i8)") oneFieldSection%xEnd
       write(str(3),"(i8)") oneFieldSection%yStart
       write(str(4),"(i8)") oneFieldSection%yEnd
       write(str(5),"(i8)") 1
       write(str(6),"(i8)") kMax

       call MsgDump(h//" buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//&
            ") <- "//trim(adjustl(oneFieldSection%name))//"("//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//","//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")")
    end if
  end subroutine FieldSectionData2BufferVariableAdressArr




  subroutine FieldSectionData2BufferVariableAdressScalar(field, kMax, oneFieldSection, &
       buf, bufStart, bufSize)

    ! copy field section values into a 1D buffer in array
    ! element order, for a 3D field with variable memory adress

    real, intent(in) :: field
    integer, intent(in) :: kMax
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
    integer :: k
    integer :: lMax
    character(len=8) :: c0, c1, c2
    character(len=8) :: buf0, bufn
    character(len=64) :: preStr, midStr
    character(len=*), parameter :: h="**(FieldSectionData2BufferVariableAdressScalar)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" null oneFieldSection")
    end if

    finalPos=bufStart+FieldSectionSize(oneFieldSection)-1
    if (finalPos > bufSize) then
       write(c0,"(i8)") FieldSectionSize(oneFieldSection)
       write(c1,"(i8)") bufStart
       write(c2,"(i8)") bufSize
       call fatal_error(h//" field "//&
            trim(adjustl(oneFieldSection%name))//&
            " with size "//trim(adjustl(c0))//&
            " is larger than buffer ("//&
            trim(adjustl(c1))//":"//trim(adjustl(c2))//")")
    end if

    if (dumpLocal) then
       write(c1,"(i8)") bufStart
       write(c2,"(i8)") bufSize
       call MsgDump(h//" invokes HiddenFieldSection...;"//&
            " bufStart="//trim(adjustl(c1))//&
            " bufSize="//trim(adjustl(c2)))
    end if

    ! hide scalar field

    call HiddenFieldSectionData2BufferVariableAdress(field, &
         oneFieldSection%name, oneFieldSection%idim_type, &
         oneFieldSection%xStart, oneFieldSection%xEnd, &
         oneFieldSection%yStart, oneFieldSection%yEnd, kMax, &
         bufStart, bufSize, buf)
  end subroutine FieldSectionData2BufferVariableAdressScalar





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
    integer :: k
    integer :: lMax
    integer :: kMax
    character(len=8) :: c0, c1, c2
    character(len=8) :: buf0, bufn
    character(len=64) :: preStr, midStr, posStr
    character(len=*), parameter :: h="**(Buffer2FieldSectionDataFixedAdress)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" null oneFieldSection")
    end if

    finalPos=bufStart+FieldSectionSize(oneFieldSection)-1
    if (finalPos > bufSize) then
       write(c0,"(i8)") FieldSectionSize(oneFieldSection)
       write(c1,"(i8)") bufStart
       write(c2,"(i8)") bufSize
       call fatal_error(h//" field "//&
            trim(adjustl(oneFieldSection%name))//&
            " with size "//trim(adjustl(c0))//&
            " is larger than buffer ("//&
            trim(adjustl(c1))//":"//trim(adjustl(c2))//")")
    end if

    if (dumpLocal) then
       write(buf0,"(i8)") bufStart
       write(bufn,"(i8)") finalPos
       preStr=""
       write(c0,"(i8)") oneFieldSection%xStart
       write(c1,"(i8)") oneFieldSection%xEnd
       midStr=trim(adjustl(c0))//":"//trim(adjustl(c1))//"," 
       write(c0,"(i8)") oneFieldSection%yStart
       write(c1,"(i8)") oneFieldSection%yEnd
       midStr=trim(adjustl(midStr))//trim(adjustl(c0))//":"//trim(adjustl(c1))
       posStr=""
    end if

    posBuf=bufStart
    select case (oneFieldSection%idim_type)
    case(2)
       if (.not. associated(oneFieldSection%field_2D)) then
          call fatal_error(h//" field "//trim(adjustl(oneFieldSection%name))//&
               " not allocated")
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
       kMax = size(oneFieldSection%field_3D,1)
       do y=oneFieldSection%yStart, oneFieldSection%yEnd
          do x=oneFieldSection%xStart, oneFieldSection%xEnd
             oneFieldSection%field_3D(:,x,y) = buf(posBuf:posBuf+kMax-1)
             posBuf=posBuf+kMax
          end do
       end do

       if (dumpLocal) then
          write(c0,"(i8)") lbound(oneFieldSection%field_3D,1)
          write(c1,"(i8)") ubound(oneFieldSection%field_3D,1)
          preStr=trim(adjustl(c0))//":"//trim(adjustl(c1))//"," 
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
                oneFieldSection%field_4D(:,x,y,k) = buf(posBuf:posBuf+lMax-1)
                posBuf=posBuf+lMax
             end do
          end do
       end do

       if (dumpLocal) then
          write(c0,"(i8)") lbound(oneFieldSection%field_4D,1)
          write(c1,"(i8)") ubound(oneFieldSection%field_4D,1)
          preStr=trim(adjustl(c0))//":"//trim(adjustl(c1))//"," 
          write(c0,"(i8)") lbound(oneFieldSection%field_4D,4)
          write(c1,"(i8)") ubound(oneFieldSection%field_4D,4)
          posStr=trim(adjustl(c0))//":"//trim(adjustl(c1))
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
                oneFieldSection%field_3D(x,y,k) = buf(posBuf)
                posBuf=posBuf+1
             end do
          end do
       end do

       if (dumpLocal) then
          write(c0,"(i8)") lbound(oneFieldSection%field_3D,3)
          write(c1,"(i8)") ubound(oneFieldSection%field_3D,3)
          posStr=trim(adjustl(c0))//":"//trim(adjustl(c1))
       end if

    end select
    bufStart=posBuf

    if (dumpLocal) then
       call MsgDump(h//trim(adjustl(oneFieldSection%name))//"("//&
            trim(adjustl(preStr))//trim(adjustl(midStr))//trim(adjustl(posStr))//&
            ") <- buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//")")

    end if
  end subroutine Buffer2FieldSectionDataFixedAdress




  subroutine Buffer2FieldSectionDataVariableAdress(oneFieldSection, &
       field, msgStartZ, msgEndZ, &
       buf, bufStart, bufSize)

    ! copy field section values from a 1D buffer to
    ! the field at field section in array element order

    type(FieldSection), pointer, intent(in) :: oneFieldSection
    real, pointer, intent(in) :: field(:,:,:)
    integer, intent(in) :: msgStartZ
    integer, intent(in) :: msgEndZ
    ! indices of the field section to copy to;
    ! field first dimension is indexed lbField:ubField
    ! but corresponding send message has first dimension
    ! only at the range msgStartZ:msgEndZ;
    ! this happens whenever the send message operates
    ! on fields with ghost zone width of 1 and the
    ! received target field has a larger ghost zone width;
    ! such care is not necessary on the other dimensions,
    ! since the correct indices are taken from the field
    ! section xStart, xEnd, yStart and yEnd components
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
    integer :: zOffset
    character(len=8) :: c0, c1, c2
    character(len=8) :: buf0, bufn
    character(len=64) :: preStr, midStr
    character(len=*), parameter :: h="**(Buffer2FieldSectionDataVariableAdress)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(oneFieldSection)) then
       call fatal_error(h//" null oneFieldSection")
    else if (.not. associated(field)) then
       call fatal_error(h//" null field")
    end if

    call UpdateFieldAdress_3D(oneFieldSection, field)
    
    finalPos=bufStart+FieldSectionSize(oneFieldSection)-1
    if (finalPos > bufSize) then
       write(c0,"(i8)") FieldSectionSize(oneFieldSection)
       write(c1,"(i8)") bufStart
       write(c2,"(i8)") bufSize
       call fatal_error(h//" field "//&
            trim(adjustl(oneFieldSection%name))//&
            " with size "//trim(adjustl(c0))//&
            " is larger than buffer ("//&
            trim(adjustl(c1))//":"//trim(adjustl(c2))//")")
    end if

    kMax=msgEndZ-msgStartZ+1

    if (dumpLocal) then
       write(buf0,"(i8)") bufStart
       write(bufn,"(i8)") finalPos
       write(c0,"(i8)") msgStartZ
       write(c1,"(i8)") msgEndZ
       preStr=trim(adjustl(c0))//":"//trim(adjustl(c1))//"," 
       write(c0,"(i8)") oneFieldSection%xStart
       write(c1,"(i8)") oneFieldSection%xEnd
       midStr=trim(adjustl(c0))//":"//trim(adjustl(c1))//"," 
       write(c0,"(i8)") oneFieldSection%yStart
       write(c1,"(i8)") oneFieldSection%yEnd
       midStr=trim(adjustl(midStr))//trim(adjustl(c0))//":"//trim(adjustl(c1))
    end if

    posBuf=bufStart
    select case (oneFieldSection%idim_type)
    case(3)
       do y=oneFieldSection%yStart, oneFieldSection%yEnd
          do x=oneFieldSection%xStart, oneFieldSection%xEnd
             zOffset=posBuf-msgStartZ
             do z=msgStartZ, msgEndZ
                field(z,x,y) = buf(z+zOffset)
             end do
             posBuf=posBuf+kMax
          end do
       end do

    case default
       write(c0,"(i8)") oneFieldSection%idim_type
       call fatal_error(h//" incompatible idim_type="//&
            trim(adjustl(c0)))
    end select
    bufStart=posBuf

    if (dumpLocal) then
       call MsgDump(h//trim(adjustl(oneFieldSection%name))//"("//&
            trim(adjustl(preStr))//trim(adjustl(midStr))//&
            ") <- buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//")")

    end if
  end subroutine Buffer2FieldSectionDataVariableAdress
end module ModFieldSection





subroutine HiddenFieldSectionData2BufferVariableAdress(field, &
     fieldName, idim_type, xStart, xEnd, yStart, yEnd, kMax, &
     bufStart, bufSize, buf)

  ! hiden interface to allow actual 1D or 3D field argument,
  ! due to invocations of FieldSectionData2BufferVariableAdress
  ! passing scalar_tab pointers and also tendency arrays, besides
  ! the usual 3D field arrays

  use ModParallelEnvironment, only: MsgDump

  implicit none
  real, intent(in) :: field(:,:,:)
  character(len=*), intent(in) :: fieldName
  integer, intent(in) :: idim_type
  integer, intent(in) :: xStart
  integer, intent(in) :: xEnd
  integer, intent(in) :: yStart
  integer, intent(in) :: yEnd
  integer, intent(in) :: kMax
  integer, intent(inout) :: bufStart
  integer, intent(in) :: bufSize
  real, intent(inout) :: buf(bufSize)

  integer :: x, y, z
  integer :: posBuf

  character(len=8) :: buf0
  character(len=8) :: bufn
  character(len=8) :: str(10)

  character(len=*), parameter :: h="**(HiddenFieldSectionData2BufferVariableAdress)**"
  logical, parameter :: dumpLocal=.true.

  if (dumpLocal) then
     write(buf0,"(i8)") bufStart
     write(str(1),"(i8)") bufSize
     call MsgDump(h//" with bufSize="//trim(adjustl(str(1))))
  end if

  posBuf=bufStart

  select case (idim_type)
  case(3)
     do y=yStart, yEnd
        do x=xStart, xEnd
           buf(posBuf:posBuf+kMax-1) = field(1:kMax,x,y)
           posBuf=posBuf+kMax
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
     write(str(6),"(i8)") kMax

     call MsgDump(h//" buf("//trim(adjustl(buf0))//":"//trim(adjustl(bufn))//&
          ") <- "//trim(adjustl(fieldName))//"("//&
          trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//","//&
          trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
          trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")")
  end if
end subroutine HiddenFieldSectionData2BufferVariableAdress
