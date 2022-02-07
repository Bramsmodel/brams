module ModFieldSection


  ! a list of field sections to be communicated
  ! among processes in a single message passing operation.
  ! List head and tail are stored elsewhere


  use ModParallelEnvironment, only: MsgDump

  implicit none
!!$  include "constants.h"

  private
  public :: FieldSection
  public :: CreateFieldSection
  public :: StringFieldSection
  public :: DumpFieldSection
  public :: DestroyFieldSection
  public :: NextFieldSection
  public :: AppendFieldSection
  public :: FieldSectionSize
  public :: FieldSectionData2Buffer
  public :: Buffer2FieldSectionData

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
     character(len=16) :: name = ""
     ! field variable name
     type(FieldSection), pointer :: next => null()
     type(FieldSection), pointer :: previous => null()
     ! double linked list of FieldSection, since
     ! one communication may have multiple fields
     ! to communicate
  end type FieldSection

  interface CreateFieldSection
     module procedure CreateFieldSection_I2D
     module procedure CreateFieldSection_2D
     module procedure CreateFieldSection_3D
     module procedure CreateFieldSection_4D
  end interface CreateFieldSection
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
    logical, parameter :: dumpLocal=.true.

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
       call DumpFieldSection(oneFieldSection, h//" has entry ")
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
    logical, parameter :: dumpLocal=.true.

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
       call DumpFieldSection(oneFieldSection, h//" has entry ")
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
    logical, parameter :: dumpLocal=.true.

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
       call DumpFieldSection(oneFieldSection, h//" has entry ")
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
    logical, parameter :: dumpLocal=.true.

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
       call DumpFieldSection(oneFieldSection, h//" has entry ")
    end if
  end function CreateFieldSection_4D





  function StringFieldSection(oneFieldSection) result(res)

    ! String with the fields of a type FieldSection variable

    type(FieldSection), pointer :: oneFieldSection
    character(len=256) :: res

    character(len=128) :: string
    character(len=8) :: c0, c1, c2, c3, c4, c5
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
          write(c4,"(i8)") size(oneFieldSection%field_3D,1)
          string="(1:"//trim(adjustl(c4))//","//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//")"
       case(4:5)
          write(c4,"(i8)") size(oneFieldSection%field_4D,1)
          write(c5,"(i8)") size(oneFieldSection%field_4D,4)
          string="(1:"//trim(adjustl(c4))//","//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               "1:"//trim(adjustl(c5))//")"
       case(6:7)
          write(c4,"(i8)") size(oneFieldSection%field_3D,3)
          string="("//&
               trim(adjustl(c0))//":"//trim(adjustl(c1))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               "1:"//trim(adjustl(c4))//")"
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
    character(len=128) :: name
    type(FieldSection), pointer :: this
    type(FieldSection), pointer :: previous
    type(FieldSection), pointer :: next
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DestroyFieldSection)**"
    logical, parameter :: dumpLocal=.true.

    if (associated(oneFieldSection)) then
       name = oneFieldSection%name
       oneFieldSection%field_2D => null()
       oneFieldSection%field_3D => null()
       oneFieldSection%field_4D => null()
       oneFieldSection%field_I2D => null()
       previous => oneFieldSection%previous
       next => oneFieldSection%next
       if (associated(previous)) then
          previous%next => next
       end if
       if (associated(next)) then
          next%previous => previous
       end if
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




  
  function NextFieldSection(node) result(next)

    ! returns node following "node" at the list;
    ! if input "node" is empty, returns list head;
    ! if no more nodes in the list, returns null

    type(FieldSection), pointer :: node
    type(FieldSection), pointer :: next

    character(len=*), parameter :: h="**(NextFieldSection)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(node)) then
       call fatal_error(h//" invoked with not associated node")
    else
       next => node%next
    end if

    if (dumpLocal) then
       call MsgDump(h//" is "//trim(adjustl(StringFieldSection(next))))
    end if
  end function NextFieldSection




  
  subroutine AppendFieldSection(this, next)

    type(FieldSection), pointer, intent(in) :: this
    type(FieldSection), pointer, intent(in) :: next

    character(len=*), parameter :: h="**(AppendFieldSection)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(this)) then
       call fatal_error(h//" null this")
    end if
    if (.not. associated(next)) then
       call fatal_error(h//" null next")
    end if
    next%previous => this
    this%next => next
    if (dumpLocal) then
       call MsgDump(h//" "//&
            trim(adjustl(next%name))//" appended to "//&
            trim(adjustl(this%name)))
    end if
  end subroutine AppendFieldSection




  
  function FieldSectionSize(oneFieldSection) result(nEle)

    ! returns the number of elements of a field section
    
    type(FieldSection), pointer, intent(in) :: oneFieldSection
    integer :: nEle

    nEle=0
    if (associated(oneFieldSection)) then
       nEle=oneFieldSection%fieldSectionSize
    end if
  end function FieldSectionSize




  
  subroutine FieldSectionData2Buffer(oneFieldSection, &
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
    character(len=*), parameter :: h="**(FieldSectionData2Buffer)**"
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
       do y=oneFieldSection%yStart, oneFieldSection%yEnd
          do x=oneFieldSection%xStart, oneFieldSection%xEnd
             buf(posBuf)=oneFieldSection%field_2D(x,y)
             posBuf=posBuf+1
          end do
       end do
    case(3)
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
  end subroutine FieldSectionData2Buffer





  subroutine Buffer2FieldSectionData(oneFieldSection, &
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
    character(len=*), parameter :: h="**(Buffer2FieldSectionData)**"
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
       do y=oneFieldSection%yStart, oneFieldSection%yEnd
          do x=oneFieldSection%xStart, oneFieldSection%xEnd
             oneFieldSection%field_2D(x,y) = buf(posBuf)
             posBuf=posBuf+1
          end do
       end do
    case(3)
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
  end subroutine Buffer2FieldSectionData
end module ModFieldSection
