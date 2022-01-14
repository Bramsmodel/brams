module InputTimeStamp

  implicit none

  private

  include "i8.h"

  type directory
     character(len=128)         :: fName
     integer                    :: MPISize
     integer                    :: count_rate
     integer                    :: count_max
     integer                    :: lastEventName
     character(len=12), pointer :: EventName(:)
     integer                    :: lastTS
  end type directory

  type TSFile
     character(len=128)  :: fName
     integer             :: MPIRank
     integer             :: lastTS
     integer(kind=i8), pointer     :: TS(:)
     integer, pointer    :: EventId(:)
     integer             :: eventMin
     integer             :: eventMax
     integer(kind=i8)    :: tsMax
  end type TSFile


  type Report
     character(len=256) :: fName
     integer            :: sizeNames
     integer, pointer   :: mapEventToPos(:)
  end type Report

  integer, parameter, public :: unitDir=59

  public :: directory
  public :: ReadDirectory
  public :: DestroyDirectory
  public :: TSFile
  public :: ReadTSFile
  public :: DumpTSFile
  public :: DestroyTSFile
  public :: Report
  public :: ReportOneTSFile
  public :: DestroyReport
contains


  subroutine ReadDirectory (machsize, dirName, dir)
    integer, intent(in) :: machsize
    character(len=*), intent(in) :: dirName
    type(directory), intent(out) :: dir

    integer :: name
    integer :: lastChar
    character(len=32) :: str1, str2,str3
    character(len=*), parameter :: h="**(ReadDirectory)**"

    dir%fName="ts.XXXXX.Directory"
    write(dir%fName(4:8),"(i5.5)") machsize
    lastChar = len_trim(dirName)
    if (dirName(lastChar:lastChar) == "/") then
       dir%fName= dirName(1:lastChar)//trim(dir%fName)
    else
       dir%fName= dirName(1:lastChar)//"/"//trim(dir%fName)
    end if

    dir%MPISize = machsize

    open (unitDir, file=trim(dir%fName), status="old", action="read")
    read (unitDir, "(i12)") dir%count_rate
    read (unitDir, "(i12)") dir%count_max
    read (unitDir, "(i12)") dir%lastEventName
    allocate(dir%EventName(dir%lastEventName))
    do name = 1, dir%lastEventName
       read (unitDir, "(a)") dir%EventName(name)
    end do
    read (unitDir, "(i12)") dir%lastTS
    close (unitDir)
    write(*, "(a)") h//" ingested file "//trim(adjustl(dir%fName))//" with" 
    write(str1,"(i32)") dir%count_rate
    write(str2,"(i32)") dir%count_max
    write(str3,"(i32)") dir%lastEventName
    write(*, "(a)") h//" count_rate="//trim(adjustl(str1))//", count_max="//trim(adjustl(str2))//&
         &" and "//trim(adjustl(str3))//" events:" 
    do name = 1, dir%lastEventName
       write(*, "(a)") h//" "//trim(dir%EventName(name))
    end do
  end subroutine ReadDirectory


  subroutine DestroyDirectory(dir)
    type(directory), intent(inout) :: dir
    deallocate(dir%EventName)
  end subroutine DestroyDirectory



  subroutine ReadTSFile (machsize, machrank, dirName, ts)
    integer, intent(in) :: machsize
    integer, intent(in) :: machrank
    character(len=*), intent(in) :: dirName
    type(TSFile), intent(out) :: ts
    integer :: i
    integer :: lastChar


    ts%fName = "tsXXXXXYYYYY.out"
    write(ts%fName(3: 7),"(i5.5)") machsize
    write(ts%fName(8:12),"(i5.5)") machrank
    lastChar = len_trim(dirName)
    if (dirName(lastChar:lastChar) == "/") then
       ts%fName= dirName(1:lastChar)//trim(ts%fName)
    else
       ts%fName= dirName(1:lastChar)//"/"//trim(ts%fName)
    end if

    ts%MPIRank=machrank

    ! File Containts

    open (unitDir, file=trim(ts%fName), status="old", form="unformatted", action="read")
    read (unitDir) ts%lastTS
    allocate(ts%TS(0:ts%lastTS))
    allocate(ts%EventId(0:ts%lastTS))
    do i = 0, ts%lastTS
       read (unitDir) ts%TS(i), ts%EventId(i)
    end do
    close (unitDir)

    ! extremes

    ts%eventMin = minval(ts%EventId(:))
    ts%eventMax = maxval(ts%EventId(:))
    ts%tsMax    = ts%TS(ts%lastTS)
  end subroutine ReadTSFile




  subroutine DumpTSFile(ts)
    type(TSFile), intent(in) :: ts

    integer :: i
    integer :: event
    integer, allocatable :: eventDist(:)
    character(len=*), parameter :: h="**(DumpTSFile)**"
    character(len=16) :: c0, c1

    write(c0,"(i16)") ts%lastTS
    write(c1,"(i16)") ts%tsMax
    write(*,"(a)") h//" file "//trim(ts%fName)//" has "//trim(adjustl(c0))//&
         " time stamps up to time "//trim(adjustl(c1))

    allocate (eventDist(ts%eventMin:ts%eventMax))
    eventDist(:) = 0
    do i = 0, ts%lastTS
       event = ts%EventId(i)
       eventDist(event) = eventDist(event) + 1
    end do

    write(*,"(a)",advance="no") h//" event:      "
    do event = ts%eventMin, ts%eventMax
       write(*, "(i6)", advance="no") event
    end do
    write(*,"(1x)")
    write(*,"(a)",advance="no") h//" occurences: "
    do event = ts%eventMin, ts%eventMax
       write(*, "(i6)", advance="no") eventDist(event)
    end do
    write(*,"(1x)")

    deallocate(eventDist)
  end subroutine DumpTSFile




  subroutine DestroyTSFile(ts)
    type(TSFile), intent(inout) :: ts
    deallocate(ts%TS)
    deallocate(ts%EventId)
  end subroutine DestroyTSFile









  subroutine ReportOneTSFile (dir, ts, rep)
    type(directory), intent(in) :: dir
    type(TSFile),    intent(in) :: ts
    type(Report),    intent(in) :: rep

    integer :: i
    integer :: iTS
    integer :: pos
    integer :: ierr
    integer(kind=i8), allocatable :: counters(:)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(ReportOneTSFile)**"

    ! open report file
    write(*, "(a,i6.6)") h//" for proc ", ts%MPIrank

    open(unit=unitDir, file=trim(rep%fName), status="old", action="write", position="append", iostat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(*,"(a)") h//"**(ERROR)** open file "//trim(rep%fName)//&
            &" to append returns iostat="//trim(adjustl(c0))
       stop
    end if

    ! process TS file

    allocate(counters(1:rep%sizeNames))
    counters = 0_i8
    do iTS = 1, ts%lastTS
       pos = rep%mapEventToPos(ts%EventId(iTS))
       counters(pos) = counters(pos) + ts%TS(iTS) - ts%TS(iTS-1)
    end do
    counters(rep%sizeNames) = ts%TS(iTS-1) - ts%TS(0)

    write(unitDir,"(i4.4)",advance="no") ts%MPIRank
    do i = 1, size(counters)
       write(unitDir,"(1x,f9.2)", advance="no") real(counters(i))/real(dir%count_rate)
    end do
    write(unitDir,"(1x)")

    ! finalize

    deallocate(counters)
    close(unitDir)
  end subroutine ReportOneTSFile




  subroutine DestroyReport(rep)
    type(Report),     intent(out) :: rep
    deallocate(rep%mapEventToPos)
  end subroutine DestroyReport
end module InputTimeStamp
