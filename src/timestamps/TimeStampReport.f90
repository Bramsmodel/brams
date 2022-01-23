program TimeStampReport

  use InputTimeStamp
  implicit none

  ! namelist 

  integer            :: machsize ! MPI size
  character(len=128) :: pathIn   ! data directory
  character(len=128) :: pathOut  ! output directory
  namelist /TSInput/ machsize, pathIn, pathOut

  integer :: mach
  integer :: ierr
  character(len=8)  :: c0, c1
  character(len=*), parameter :: h="**(TimeStampReport)**" 

  type(directory) :: dir
  type(TSFile) :: ts
  type(Report) :: rep

  ! read namelist and echo 

  read (*, NML=TSInput, iostat=ierr)
  if (ierr == 0) then
     write(*, "(a,i4.4)") h//" MPI rank is ",machsize
     write(*, "(a)") h//" Input Directory: "//trim(pathIn)
     write(*, "(a)") h//" Output Directory: "//trim(pathOut)
  else
     write(c0,"(i8)") ierr
     write(*,"(a)") "**(ERROR)** read namelist returns istat="&
          &//trim(adjustl(c0))
     stop
  end if

  ! read directory

  call ReadDirectory(machsize, pathIn, dir)

  ! read proc 0 trace file; write first line and proc 0 line; destroy trace file

  call ReadTSFile (machsize, 0, pathIn, ts)
  call CreateReport (dir, pathOut, ts, rep)
  call ReportOneTSFile (dir, ts, rep)
  call DestroyTSFile(ts)

  ! read one trace file, report and destroy

  do mach = 1, machsize-1
     call ReadTSFile (machsize, mach, pathIn, ts)
     call ReportOneTSFile (dir, ts, rep)
     call DestroyTSFile(ts)
  end do

  ! finalize

  call DestroyDirectory(dir)
  call DestroyReport(rep)
  stop
contains
  subroutine CreateReport (dir, dirName, ts, rep)
    use InputTimeStamp
    type(directory),  intent(in) :: dir
    character(len=*), intent(in) :: dirName
    type(TSFile),     intent(in) :: ts
    type(Report),     intent(out) :: rep

    integer :: iTS
    integer :: iTSFile
    integer :: lastChar
    integer :: i
    integer :: pos
    integer :: eventNumber
    integer :: ierr
    integer :: posNames
    integer :: posEventNames
    integer :: sizeNames
    logical, allocatable :: SyncEvents(:)
    character(len=9), allocatable :: names(:)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(CreateReport)**"

    ! report file name

    rep%fName="Report.XXXXX"
    write(rep%fName(8:12), "(i5.5)") dir%MPISize
    lastChar = len_trim(dirName)
    if (dirName(lastChar:lastChar) == "/") then
       rep%fName = dirName(1:lastChar)//trim(rep%fName)
    else
       rep%fName = dirName(1:lastChar)//"/"//trim(rep%fName)
    end if

    ! open report file

    open(unit=unitDir, file=trim(rep%fName), status="replace", action="write", iostat=ierr)
    if (ierr == 0) then
       write(*,"(a)") h//" Building file "//trim(rep%fName)
    else
       write(c0,"(i8)") ierr
       write(*,"(a)") h//"**(ERROR)** open file "//trim(rep%fName)//&
            &" returns iostat="//trim(adjustl(c0))
       stop
    end if

    ! which events are synchronized (requiring a new name)

    allocate (syncEvents(dir%lastEventName), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") dir%lastEventName
       write(*,"(a)") h//"**(ERROR)** allocate syncEvents("//trim(adjustl(c1))//&
            &") returns iostat="//trim(adjustl(c0))
       stop
    end if
    do i = 1, dir%lastEventName
       eventNumber = -i
       syncEvents(i) = any(ts%EventId(:) == eventNumber)
    end do

    ! event names on report's first line printing position

    rep%sizeNames = dir%lastEventName + count(syncEvents(:)) + 1  ! events, synchronizations, total
    allocate(names(rep%sizeNames), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") rep%sizeNames
       write(*,"(a)") h//"**(ERROR)** allocate names("//trim(adjustl(c1))//&
            &") returns iostat="//trim(adjustl(c0))
       stop
    end if

    ! map event id to printing position

    allocate(rep%mapEventToPos(-dir%lastEventName:dir%lastEventName), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") dir%lastEventName
       write(*,"(a)") h//"**(ERROR)** allocate mapEventToPos(-"//trim(adjustl(c1))//&
            &":"//trim(adjustl(c1))//") returns iostat="//trim(adjustl(c0))
       stop
    end if
    rep%mapEventToPos = 0

    ! position events (and possible synchronizations) for printing

    posNames=0
    do posEventNames = 1, dir%lastEventName  ! skip event 0 (baseline)
       posNames = posNames+1
       names(posNames) = dir%EventName(posEventNames)
       rep%mapEventToPos(posEventNames) = posNames
       if (syncEvents(posEventNames)) then
          posNames = posNames+1
          names(posNames) = trim(dir%EventName(posEventNames))//"Sync"
          rep%mapEventToPos(-posEventNames) = posNames
       end if
    end do
    posNames = posNames + 1
    if (posNames == rep%sizeNames) then
       names(posNames) = "Total"
    else
       write(c0,"(i8)") posNames
       write(c1,"(i8)") rep%sizeNames
       write(*,"(a)") h//"**(ERROR)** logical error: posNames ("//trim(adjustl(c0))//&
            &") and sizeNames ("//trim(adjustl(c1))//") differ"
       stop
    end if

    ! header line

    write(unitDir,"(a4)", advance="no") "proc"
    do i = 1, size(names)
       write(unitDir,"(1x,a9)", advance="no") adjustr(names(i))
    end do
    write(unitDir,"(1x)")

    close(unitDir)
    deallocate(names)
    deallocate(syncEvents)
  end subroutine CreateReport
end program TimeStampReport
