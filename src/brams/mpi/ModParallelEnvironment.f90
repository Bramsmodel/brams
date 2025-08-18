module ModParallelEnvironment

  ! ModParallelEnvironment: creates a type that stores info about the 
  !                         MPI parallel environment. Creates a dump file
  !                         for each process. Provides procedures for
  !                         dumping a msg at the dump file or at stdout.

  use iso_fortran_env, only: &
       output_unit

  implicit none
  
  private
  public :: ParallelEnvironment
  public :: CreateParallelEnvironment
  public :: DestroyParallelEnvironment
  public :: MsgDump
  public :: MsgOutput
  public :: Brams2MpiProcNbr
  public :: Mpi2BramsProcNbr
  public :: GetNumberOfProcesses
  public :: GetThisBramsProcessNumber
  public :: DumpUnit

  type ParallelEnvironment
     integer :: communicator 
     integer :: nmachs             ! number of processes (MPI size)
     integer :: mchnum             ! MPI number of this process
     integer :: master_num         ! MPI number of the master process
     integer :: myNum              ! Brams number of this process
  end type ParallelEnvironment

  include "mpif.h"
  integer :: DumpUnit
  character(len=6) :: strMchNum

  

contains



  !*** CreateParallelEnvironment: create and fill variable of this type



  subroutine CreateParallelEnvironment(nmachs, mchnum, master_num, &
       comm, oneParallelEnvironment)
    integer, intent(in) :: nmachs     ! number of processes (0 iff sequential run)
    integer, intent(in) :: master_num ! this process rank (0:nmachs-1); 0 on sequential runs
    integer, intent(in) :: mchnum     ! this process rank (0:nmachs-1); 0 on sequential runs
    integer, intent(in) :: comm       ! MPI communicator

    integer :: ierr
    logical :: op, ex
    character :: size, rank
    character(len=512) :: message
    character(len=16) :: dumpFName="Dump.XXXXX.YYYYY"
    type(ParallelEnvironment), pointer :: oneParallelEnvironment

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateParallelEnvironment)**" 
    logical, parameter :: dumpLocal=.false.
    
    ! creates variable and fill components

    allocate(oneParallelEnvironment)
    oneParallelEnvironment%nmachs = nmachs
    oneParallelEnvironment%mchnum = mchnum
    oneParallelEnvironment%master_num = master_num
    oneParallelEnvironment%communicator = comm
    oneParallelEnvironment%myNum = Mpi2BramsProcNbr(mchnum)

    write(dumpFName(6:10),"(i5.5)") mchnum
    write(strMchNum,"(a1,i5.5)") "_", mchnum
    write(dumpFName(12:16),"(i5.5)") nmachs
    open(newUnit=DumpUnit, file=dumpFName, action="write", status="replace", iostat=ierr, iomsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" open "//dumpFName//" fails with message "//trim(message))
    end if
    if (dumpLocal) then
       write(str(1),"(i8)") DumpUnit
       call MsgDump(h//" abriu arquivo "//dumpFName//" na unidade "//trim(adjustl(str(1))))
       call MsgOutput(h//strMchNum//" abriu arquivo "//dumpFName//" na unidade "//trim(adjustl(str(1))))
    end if
  end subroutine CreateParallelEnvironment



  ! DestroyParallelEnvironment: returns storage, destroy variable
  !                             closes dump file



  subroutine DestroyParallelEnvironment(oneParallelEnvironment)
    type(parallelEnvironment), pointer :: oneParallelEnvironment

    logical :: op
    if (associated(oneParallelEnvironment)) then
       deallocate(oneParallelEnvironment)

       ! only the first execution of this procedure closes dump file

       inquire(unit=DumpUnit, opened=op)
       if (op) then
          close(DumpUnit)
       end if
    end if
    nullify(oneParallelEnvironment)
  end subroutine DestroyParallelEnvironment



  ! MsgDump: Writes string at dump file



  subroutine MsgDump(msg, noAdvance)
    character(len=*), intent(in) :: msg
    logical, optional, intent(in) :: noAdvance

    integer :: mchnum
    integer :: nmachs
    integer :: ierr
    logical :: op
    character(len=512) :: message
    character(len=16) :: dumpFName="Dump.XXXXX.YYYYY"

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(MsgDump)**"
    logical, parameter :: dumpLocal=.false.

    if (present(noAdvance)) then
       if (noAdvance) then
          write(DumpUnit,"(a)",advance="no") trim(msg)
       else
          write(DumpUnit,"(a)") trim(msg)
       end if
    else
       write(DumpUnit,"(a)") trim(msg)
       if (dumpLocal) then
          write(str(1),"(i8)") DumpUnit
          call MsgOutput(h//strMchNum//" escreveu "//trim(msg)//" na unidade "//trim(adjustl(str(1))))
       end if
    end if
    flush(DumpUnit)
    if (dumpLocal) then
       write(str(1),"(i8)") DumpUnit
       call MsgOutput(h//strMchNum//" flush na unidade "//trim(adjustl(str(1))))
    end if
    if (present(noAdvance)) then
       if (noAdvance) then
          call MsgOutput(h//strMchNum//" "//trim(msg), noAdvance)
       else
          call MsgOutput(h//strMchNum//" "//trim(msg))
       end if
    else
       call MsgOutput(h//strMchNum//" "//trim(msg))
    end if
  end subroutine MsgDump



  ! MsgOutput: Writes string at output file



  subroutine MsgOutput(str, noAdvance)
    character(len=*), intent(in) :: str
    logical, optional, intent(in) :: noAdvance

    if (present(noAdvance)) then
       if (noAdvance) then
          write(output_unit,"(a)",advance="no") trim(str)
          flush(output_unit)
          return
       end if
    end if
    write(output_unit,"(a)") trim(str)
    flush(output_unit)
  end subroutine MsgOutput




  integer function Brams2MpiProcNbr(procNbr)
    integer, intent(in) :: procNbr
    Brams2MpiProcNbr = procNbr-1
  end function Brams2MpiProcNbr




  integer function Mpi2BramsProcNbr(procNbr)
    integer, intent(in) :: procNbr
    Mpi2BramsProcNbr = procNbr+1
  end function Mpi2BramsProcNbr




  integer function GetNumberOfProcesses(oneParallelEnvironment)
    type(parallelEnvironment), pointer :: oneParallelEnvironment

    character(len=*), parameter :: h="**(GetNumberOfProcesses)**"

    if (.not. associated(oneParallelEnvironment)) then
       call fatal_error(h//" oneParallelEnvironment not associated")
    else
       GetNumberOfProcesses = oneParallelEnvironment%nmachs
    end if
  end function GetNumberOfProcesses




  integer function GetThisBramsProcessNumber(oneParallelEnvironment)
    type(parallelEnvironment), pointer :: oneParallelEnvironment

    character(len=*), parameter :: h="**(GetThisBramsProcessNumber)**"

    if (.not. associated(oneParallelEnvironment)) then
       call fatal_error(h//" oneParallelEnvironment not associated")
    else
       GetThisBramsProcessNumber = oneParallelEnvironment%myNum
    end if
  end function GetThisBramsProcessNumber
end module ModParallelEnvironment

