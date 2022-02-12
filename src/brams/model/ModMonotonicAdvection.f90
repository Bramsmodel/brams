!----------------------------------------------------------------------!
! Optional advection scheme for CCATT-BRAMS/BRAMS models version 4.2+  !
! Based on Walcek, 2000 (JGR) and Walcek and Aleksic, 1998 (ATENV).    !
! The scheme is highly conservative, monotonic and keeps mass mixing   !
! ratio positive definite. 					       !
! Implemented by Saulo Freitas (saulo.freitas@cptec.inpe.br) @ Jun/2009!
! MPI/Paralelized by L. Flavio/J. Panneta                              !
!----------------------------------------------------------------------!

module ModMonotonicAdvection

  use ModParallelEnvironment, only: ParallelEnvironment
  use ModParallelEnvironment, only: MsgDump
  use ModGridDims, only: GridDims
  use ModDomainDecomp, only: DomainDecomp

  implicit none

  private
  public :: MonotonicAdvection
  public :: CreateMonotonicAdvection
  
  type MonotonicAdvection
     real,pointer :: u3d(:,:,:)
     real,pointer :: v3d(:,:,:)
     real,pointer :: w3d(:,:,:)
     real,pointer :: vc3d_in(:,:,:)
     real,pointer :: vc3d_out(:,:,:)
     real,pointer :: vc3d_x(:,:,:)
     real,pointer :: vc3d_y(:,:,:)
     real,pointer :: dd0_3d(:,:,:)
     real,pointer :: dd0_3du(:,:,:)
     real,pointer :: dd0_3dv(:,:,:)
     real,pointer :: dd0_3dw(:,:,:)
     real,pointer :: den0_3d(:,:,:)
     real,pointer :: den1_3d(:,:,:)
     real,pointer :: den2_3d(:,:,:)
     real,pointer :: den3_3d(:,:,:)
     real,pointer :: l_dxtW(:,:,:)
     real,pointer :: l_dytW(:,:,:)
     real,pointer :: dxtW(:,:)
     real,pointer :: dytW(:,:)
     real,pointer :: dztW(:)
  end type MonotonicAdvection

contains




  function CreateMonotonicAdvection(&
       oneParallelEnvironment, &
       oneGridDims, &
       oneLocalOwnMonoAdv)      result(oneMonoAdv)
    type(ParallelEnvironment), pointer, intent(in) :: oneParallelEnvironment
    type(GridDims), pointer, intent(in) :: oneGridDims
    type(DomainDecomp), pointer, intent(in) :: oneLocalOwnMonoAdv
    type(MonotonicAdvection), pointer :: oneMonoAdv

    integer :: myNum
    integer :: mzp
    integer :: mxpMonoAdv
    integer :: mypMonoAdv

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateMonoAdv)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(oneParallelEnvironment)) then
       call fatal_error(h//" oneParallelEnvironment not associated")
    end if
    if (.not. associated(oneGridDims)) then
       call fatal_error(h//" oneGridDims not associated")
    end if
    if (.not. associated(oneLocalOwnMonoAdv)) then
       call fatal_error(h//" oneLocalOwnMonoAdv not associated")
    end if

    myNum = oneParallelEnvironment%myNum
    mzp = oneGridDims%nnzp
    mxpMonoAdv = oneLocalOwnMonoAdv%nx(myNum)
    mypMonoAdv = oneLocalOwnMonoAdv%ny(myNum)

    allocate(oneMonoAdv, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%u3d(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%u3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%v3d(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%v3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%w3d(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%w3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%dd0_3d(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%dd0_3d  fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%dd0_3du(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%dd0_3du fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%dd0_3dv(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%dd0_3dv fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%dd0_3dw(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%dd0_3dw fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%den0_3d(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%den0_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%den1_3d(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%den1_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%den2_3d(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%den2_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%den3_3d(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%den3_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%l_dxtW(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%l_dxtW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%l_dytW(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%l_dytW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%dxtW(mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%dxtW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%dytW(mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%dytW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%dztW(mzp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%dztW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%vc3d_in(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%vc3d_in  fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneMonoAdv%vc3d_out(mzp,mxpMonoAdv,mypMonoAdv), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneMonoAdv%vc3d_out fails with stat="//&
            trim(adjustl(str(1))))
    end if

    if (dumpLocal) then
       write(str(2),"(i8)") mzp
       write(str(3),"(i8)") mxpMonoAdv
       write(str(4),"(i8)") mypMonoAdv
       call MsgDump(h//" allocated oneMonoAdv%u3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%v3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%w3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%dd0_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%dd0_3du("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%dd0_3dv("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%dd0_3dw("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%den0_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%den1_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%den2_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%den3_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%l_dxtW("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%l_dytW("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%dxtW("//&
            trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%dytW("//&
            trim(adjustl(str(3)))//")")
       call MsgDump(h//" allocated oneMonoAdv%dztW("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%vc3d_in("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneMonoAdv%vc3d_out("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
    end if
  end function CreateMonotonicAdvection
end module ModMonotonicAdvection
