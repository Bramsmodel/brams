module ModDumpFields

  ! dumps one field on a file;
  ! compares one field with a file contents;
  ! handy to find discrepancies among two executions

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  implicit none

  private

  public :: timeStepNumber
  public :: CreateDumpFields
  public :: DumpField
  
  logical :: created = .false.

  integer, parameter :: maxFileName=256
  character(len=maxFileName) :: prefix = ""

  integer :: timeStepNumber=-1

  integer, external :: AvailableFileUnit
  
  interface DumpField
     module procedure DumpField_r1d
     module procedure DumpField_r2d
     module procedure DumpField_r3d
     module procedure DumpField_r4d
  end interface DumpField

contains




  subroutine CreateDumpFields(&
       inputParallelEnvironment, &
       inputPrefix)
    type(ParallelEnvironment), pointer :: inputParallelEnvironment
    character(len=*), intent(in) :: inputPrefix

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateDumpFields)**"
    logical, parameter :: dumpLocal=.false.

    if (len_trim(inputPrefix) > maxFileName) then
       write(str(1),"(i8)") len_trim(inputPrefix)
       write(str(2),"(i8)") maxFileName
       call fatal_error(h//" inputPrefix length ("//trim(adjustl(str(1)))//&
            ") > maxFileName ("//trim(adjustl(str(2)))//")")
    end if
    prefix = inputPrefix

    if (.not. associated(inputParallelEnvironment)) then
       call fatal_error(h//" inputParallelEnvironment not associated")
    end if

    if (inputParallelEnvironment%nmachs /= 1) then
       call fatal_error(h//" not ready for parallel executions")
    end if

    created = .true.

    if (dumpLocal) then
       call MsgDump(h//" defined prefix="//trim(prefix))
    end if
  end subroutine CreateDumpFields
       




  subroutine CheckEnvironment(h)
    character(len=*), intent(in) :: h

    if (.not. created) then
       call fatal_error(h//" ModDumpFields not created")
    end if
  end subroutine CheckEnvironment


  

  function BuildFileName(h, inputFileName) result(fileName)
    character(len=*), intent(in) :: h
    character(len=*), intent(in) :: inputFileName
    character(len=maxFileName) :: fileName

    integer :: lenFileName=-1
    character(len=MaxFileName) :: strTimeStepNumber=""

    character(len=8) :: str(10)
    character(len=*), parameter :: hLocal="**(BuildFileName)**"
    logical, parameter :: dumpLocal=.false.
    
    if ((timeStepNumber > 9999) .or. (timeStepNumber < 0)) then
       write(str(1),"(i8)") timeStepNumber
       call fatal_error(h//" timeStepNumber ("//trim(adjustl(str(1)))//&
            ") out of range [0:9999[")
    end if
    write(strTimeStepNumber,"('_',i4.4,'.bin')") timeStepNumber
    
    lenFileName = len_trim(prefix)+len_trim(inputFileName) + len_trim(strTimeStepNumber)
    if (lenFileName > maxFileName) then
       write(str(1),"(i8)") lenFileName
       write(str(2),"(i8)") maxFileName
       call fatal_error(h//" fileName length ("//trim(adjustl(str(1)))//&
            ") > maxFileName ("//trim(adjustl(str(2)))//")")
    end if
    fileName=trim(prefix)//trim(inputFileName)//trim(strTimeStepNumber)

    if (dumpLocal) then
       call MsgDump(h//hLocal//" defined fileName="//trim(fileName))
    end if
  end function BuildFileName


  
    
  subroutine DumpField_r1d(field, inputFileName)
    real, intent(in) :: field(:)
    character(len=*), intent(in) :: inputFileName

    character(len=maxFileName) :: fileName=""
    integer :: iun
    integer :: ios
    
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpField_r1d)**"
    logical, parameter :: dumpLocal=.false.

    call CheckEnvironment(h)
    
    fileName = BuildFileName(h, inputFileName)

    iun = AvailableFileUnit()

    open(unit=iun, file=trim(fileName), action="write", form="unformatted", status="replace", iostat=ios)
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" open "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    write(iun, iostat=ios) field
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" write "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    close(iun, iostat=ios)
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" close "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    if (dumpLocal) then
       call MsgDump(h//" wrote file "//trim(fileName))
    end if
  end subroutine DumpField_r1d





  subroutine DumpField_r2d(field, inputFileName)
    real, intent(in) :: field(:,:)
    character(len=*), intent(in) :: inputFileName

    character(len=maxFileName) :: fileName=""
    integer :: iun
    integer :: ios
    
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpField_r2d)**"
    logical, parameter :: dumpLocal=.false.

    call CheckEnvironment(h)
    
    fileName = BuildFileName(h, inputFileName)

    iun = AvailableFileUnit()

    open(unit=iun, file=trim(fileName), action="write", form="unformatted", status="replace", iostat=ios)
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" open "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    write(iun, iostat=ios) field
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" write "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    close(iun, iostat=ios)
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" close "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    if (dumpLocal) then
       call MsgDump(h//" wrote file "//trim(fileName))
    end if
  end subroutine DumpField_r2d





  subroutine DumpField_r3d(field, inputFileName)
    real, intent(in) :: field(:,:,:)
    character(len=*), intent(in) :: inputFileName

    character(len=maxFileName) :: fileName=""
    integer :: iun
    integer :: ios
    
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpField_r3d)**"
    logical, parameter :: dumpLocal=.false.

    call CheckEnvironment(h)
    
    fileName = BuildFileName(h, inputFileName)

    iun = AvailableFileUnit()

    open(unit=iun, file=trim(fileName), action="write", form="unformatted", status="replace", iostat=ios)
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" open "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    write(iun, iostat=ios) field
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" write "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    close(iun, iostat=ios)
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" close "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    if (dumpLocal) then
       call MsgDump(h//" wrote file "//trim(fileName))
    end if
  end subroutine DumpField_r3d





  subroutine DumpField_r4d(field, inputFileName)
    real, intent(in) :: field(:,:,:,:)
    character(len=*), intent(in) :: inputFileName

    character(len=maxFileName) :: fileName=""
    integer :: iun
    integer :: ios
    
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpField_r4d)**"
    logical, parameter :: dumpLocal=.false.

    call CheckEnvironment(h)
    
    fileName = BuildFileName(h, inputFileName)

    iun = AvailableFileUnit()

    open(unit=iun, file=trim(fileName), action="write", form="unformatted", status="replace", iostat=ios)
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" open "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    write(iun, iostat=ios) field
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" write "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    close(iun, iostat=ios)
    if (ios /= 0) then
       write(str(1),"(i8)") ios
       call fatal_error(h//" close "//trim(fileName)//" fails with iostat="//trim(adjustl(str(1))))
    end if

    if (dumpLocal) then
       call MsgDump(h//" wrote file "//trim(fileName))
    end if
  end subroutine DumpField_r4d
end module ModDumpFields
