!----------------------------------------------------------------------!
! Optional advection scheme for CCATT-BRAMS/BRAMS models version 4.2+  !
! Based on Walcek, 2000 (JGR) and Walcek and Aleksic, 1998 (ATENV).    !
! The scheme is highly conservative, monotonic and keeps mass mixing   !
! ratio positive definite. 					       !
! Implemented by Saulo Freitas (saulo.freitas@cptec.inpe.br) @ Jun/2009!
! MPI/Paralelized by L. Flavio/J. Panneta                              !
!----------------------------------------------------------------------!

module ModMonotonicAdvection

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModMessageSet, only: &
       UpdateFieldAdressAtAdvMnt, &
       PostSendRecvMsgs, &
       WaitSendRecvMsgs

  use ModDomainDecomp, only: &
       DomainDecomp

  use ModGrid, only: &
       Grid, &
       DumpGrid

  use ModScalarTable, only: &
       DeepCopyToScalarTab, &
       DeepCopyFromScalarTab       

  use ModNamelistFile, only: &
       NamelistFile

  use mem_grid, only:        &
       dtlt,   & !intent(in)
       time,   &
       ngrid,  & !intent(in)
       dzt,    & !intent(in)
       dztn,   & !intent(in)
       grid_g, & !intent(in)
       hw4   , & !intent(in)
       if_adap,& !intent(in)
       dyncore_flag  !intent(in)

  use mem_basic, only: &
       basic_g  !intent(in)

  use micphys, only: &
       level !intent(in)

  use rconstants, only: &
       cp,p00,cv,rgas,cpi   !intent(in)

  use mem_aer1, only: &
       aerosol,    &       !intent(in)
       num_scalar_aer_1st !intent(in)

  use mem_chem1, only: &
       nspecies_transported !intent(in)

  use module_dry_dep, only: &
       sedim_type,          &
       dd_sedim,            &
       naer_transported

  use var_tables, only : &
       scalar_tab, & ! (var_p = IN, var_t = INOUT)
       num_scalar   ! (in)

  use ccatt_start, only: &
       ccatt               ! (in)

  implicit none

  private

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
     real,pointer :: dxtW(:,:)
     real,pointer :: dytW(:,:)
     real,pointer :: dztW(:)
  end type MonotonicAdvection

  public :: MonotonicAdvection
  public :: CreateMonotonicAdvection
  public :: DestroyMonotonicAdvection
  public :: advmnt_driver  ! Subroutine
  public :: StoreNamelistFileAtAdvMnt ! Subroutine

  ! public names, set by StoreNamelistFileAtAdvMnt
  integer, public :: advmnt 
  integer, public :: GhostZoneLength 

  ! module private variables

  ! flow control flags
  logical, parameter :: use_true_density=.true.
  ! for theoretical experiments
  logical, parameter :: theor_wind=.false.

  ! constants
  real, parameter :: c1 = cv/rgas
  real, parameter :: c2 = p00/rgas

contains




  function CreateMonotonicAdvection(oneGrid) result(oneAdvMnt)
    type(Grid), pointer, intent(in) :: oneGrid
    type(MonotonicAdvection), pointer :: oneAdvMnt

    integer :: mzpAdvMnt
    integer :: mxpAdvMnt
    integer :: mypAdvMnt

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateMonotonicAdvection)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneGrid)) then
       call fatal_error(h//" oneGrid not associated")
    end if
    if (.not. associated(oneGrid%NodeDimsAdvMnt)) then
       call fatal_error(h//" oneGrid%NodeDimsAdvMnt not associated")
    end if

    mzpAdvMnt = oneGrid%NodeDimsAdvMnt%mzp
    mxpAdvMnt = oneGrid%NodeDimsAdvMnt%mxp
    mypAdvMnt = oneGrid%NodeDimsAdvMnt%myp

    allocate(oneAdvMnt, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%u3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%u3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%v3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%v3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%w3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%w3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dd0_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dd0_3d  fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dd0_3du(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dd0_3du fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dd0_3dv(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dd0_3dv fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dd0_3dw(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dd0_3dw fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%den0_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%den0_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%den1_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%den1_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%den2_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%den2_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%den3_3d(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%den3_3d fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dxtW(mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dxtW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dytW(mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dytW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%dztW(mzpAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%dztW fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%vc3d_in(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%vc3d_in  fails with stat="//&
            trim(adjustl(str(1))))
    end if
    allocate(oneAdvMnt%vc3d_out(mzpAdvMnt,mxpAdvMnt,mypAdvMnt), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneAdvMnt%vc3d_out fails with stat="//&
            trim(adjustl(str(1))))
    end if

    if (dumpLocal) then
       write(str(2),"(i8)") mzpAdvMnt
       write(str(3),"(i8)") mxpAdvMnt
       write(str(4),"(i8)") mypAdvMnt
       call MsgDump(h//" allocated oneAdvMnt%u3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%v3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%w3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dd0_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dd0_3du("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dd0_3dv("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dd0_3dw("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%den0_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%den1_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%den2_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%den3_3d("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dxtW("//&
            trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dytW("//&
            trim(adjustl(str(3)))//")")
       call MsgDump(h//" allocated oneAdvMnt%dztW("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%vc3d_in("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       call MsgDump(h//" allocated oneAdvMnt%vc3d_out("//&
            trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
    end if
  end function CreateMonotonicAdvection





  subroutine DestroyMonotonicAdvection(oneAdvMnt)
    type(MonotonicAdvection), pointer, intent(inout) :: oneAdvMnt

    integer :: ierr
    character(len=8) :: str(10)
    logical :: dumpLocal=.false.
    character(len=*), parameter :: h="**(DestroyMonotonicAdvection)**"

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    if (associated(oneAdvMnt)) then
       deallocate(oneAdvMnt%u3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%u3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%v3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%v3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%w3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%w3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dd0_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dd0_3d  fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dd0_3du, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dd0_3du fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dd0_3dv, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dd0_3dv fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dd0_3dw, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dd0_3dw fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%den0_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%den0_3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%den1_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%den1_3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%den2_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%den2_3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%den3_3d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%den3_3d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dxtW, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dxtW fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dytW, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dytW fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%dztW, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%dztW fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%vc3d_in, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%vc3d_in  fails with stat="//&
               trim(adjustl(str(1))))
       end if
       deallocate(oneAdvMnt%vc3d_out, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneAdvMnt%vc3d_out fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    nullify(oneAdvMnt)

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine DestroyMonotonicAdvection




  subroutine advmnt_driver(oneGrid, varn, &
       m1 ,m2 ,m3 ,ia,iz,ja,jz,izu,jzv,&
       i0,j0,nodemyp,nodemxp,nodemzp,mynum)

    type(Grid), pointer, intent(in) :: oneGrid
    integer , intent(in) :: m1
    integer , intent(in) :: m2
    integer , intent(in) :: m3
    integer , intent(in) :: ia
    integer , intent(in) :: iz
    integer , intent(in) :: ja
    integer , intent(in) :: jz
    integer , intent(in) :: izu
    integer , intent(in) :: jzv
    integer , intent(in) :: i0
    integer , intent(in) :: j0
    integer, intent(in) :: nodemxp(:,:)
    integer, intent(in) :: nodemyp(:,:)
    integer, intent(in) :: nodemzp(:,:)
    integer , intent(in) :: mynum
    character(len=*),intent(in) :: varn

    !--- local vars
    integer :: n
    integer :: ng
    integer :: i
    integer :: j
    integer :: k
    integer :: iExtern
    integer :: jExtern
    integer :: i_scl
    integer :: current_aer_ispc
    integer :: current_ndt_z
    integer, target :: ndt_z(naer_transported)
    integer, target :: ndtZ(naer_transported)
    logical  :: IsThisScalarAer =.false.

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(advmnt_driver)**"
    character(len=8) :: str(11)

    type(MonotonicAdvection), pointer :: oneAdvMnt

    integer :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer :: mxp
    ! x dimension of external fields 
    integer :: myp
    ! y dimension of external fields 
    integer :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection

    !- This scheme is not applied to advect  U, V, and W or shaved-eta

    if (varn .eq. 'V' .or. varn .eq. 'ALL') then
       call fatal_error(h//' not using mnt to advect u,v,w')
    end if
    if (if_adap /= 0) then
       call fatal_error(h//' MNT advection not ready for shaved eta')
    end if
    if (mynum == 0) then
       call fatal_error(h//' ADV MNT called with mynum = 0, try np = 2')
    end if

    ! dimension of external fields (regular ghost zone width)

    mzp=oneGrid%NodeDims%mzp
    mxp=oneGrid%NodeDims%mxp
    myp=oneGrid%NodeDims%myp

    ! dimension of Monotonic Advection fields (wide ghost zone width)

    mxpAdvMnt=oneGrid%NodeDimsAdvMnt%mxp
    mypAdvMnt=oneGrid%NodeDimsAdvMnt%myp

    ! index external = index Monotonic Advection + offset

    iOffset = oneGrid%NodeDimsAdvMnt%i0 - oneGrid%NodeDims%i0 
    i1ExternAtAdvMnt = 1 - iOffset
    iMxpExternAtAdvMnt = mxp - iOffset
    jOffset = oneGrid%NodeDimsAdvMnt%j0 - oneGrid%NodeDims%j0 
    j1ExternAtAdvMnt = 1 - jOffset
    jMypExternAtAdvMnt = myp - jOffset

    ! current grid
    ! necessary while there are global variables outside oneGrid

    ng = OneGrid%Id

    if (dumpLocal) then
       call MsgDump(h//" starts")
       write(str(1),"(i8)") mzp
       call MsgDump(h//"mzp="//trim(adjustl(str(1))))
       write(str(1),"(i8)") mxp
       call MsgDump(h//"mxp="//trim(adjustl(str(1))))
       write(str(1),"(i8)") myp
       call MsgDump(h//"myp="//trim(adjustl(str(1))))
       write(str(1),"(i8)") mxpAdvMnt
       call MsgDump(h//"mxpAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") mypAdvMnt
       call MsgDump(h//"mypAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") iOffset
       call MsgDump(h//"iOffset="//trim(adjustl(str(1))))
       write(str(1),"(i8)") i1ExternAtAdvMnt
       call MsgDump(h//"i1ExternAtAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") iMxpExternAtAdvMnt
       call MsgDump(h//"iMxpExternAtAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") jOffset
       call MsgDump(h//"jOffset="//trim(adjustl(str(1))))
       write(str(1),"(i8)") j1ExternAtAdvMnt
       call MsgDump(h//"j1ExternAtAdvMnt="//trim(adjustl(str(1))))
       write(str(1),"(i8)") jMypExternAtAdvMnt
       call MsgDump(h//"jMypExternAtAdvMnt="//trim(adjustl(str(1))))
    end if

    ! create local memory area for large GhostZoneWidth variables

    oneAdvMnt => CreateMonotonicAdvection(oneGrid)

    if(.not. use_true_density) then
       call InitializeDensities(mzp, mxp, myp, &
            mxpAdvMnt, mypAdvMnt, &
            iOffset, i1ExternAtAdvMnt, iMxpExternAtAdvMnt,  &
            jOffset, j1ExternAtAdvMnt, jMypExternAtAdvMnt,  &
            basic_g(ng)%dn0, basic_g(ng)%dn0u, basic_g(ng)%dn0v, &
            oneAdvMnt%dd0_3d, oneAdvMnt%dd0_3du, &
            oneAdvMnt%dd0_3dv, oneAdvMnt%dd0_3dw)
    end if

    ! update field addresses for message passing from recently allocated local memory area

    call UpdateFieldAdressAtAdvMnt(&
         oneGrid%AdvMntUVSendX, oneGrid%AdvMntUVRecvX, &
         oneGrid%AdvMntUVSendY, oneGrid%AdvMntUVRecvY, &
         oneGrid%AdvMntDxDySendX, oneGrid%AdvMntDxDyRecvX, &
         oneGrid%AdvMntDxDySendY, oneGrid%AdvMntDxDyRecvY, &
         oneGrid%AdvMntDd0SendX, oneGrid%AdvMntDd0RecvX, &
         oneGrid%AdvMntDd0SendY, oneGrid%AdvMntDd0RecvY, &
         oneGrid%AdvMntDenSendX, oneGrid%AdvMntDenRecvX, &
         oneGrid%AdvMntDenSendY, oneGrid%AdvMntDenRecvY, &
         oneGrid%AdvMntScaSendX, oneGrid%AdvMntScaRecvX, &
         oneGrid%AdvMntScaSendY, oneGrid%AdvMntScaRecvY, &
         oneAdvMnt%u3d, oneAdvMnt%v3d, &
         oneAdvMnt%dxtW, oneAdvMnt%dytW, &
         oneAdvMnt%dd0_3d, oneAdvMnt%dd0_3du, &
         oneAdvMnt%dd0_3dv, oneAdvMnt%dd0_3dw, &
         oneAdvMnt%den0_3d, oneAdvMnt%den1_3d, &
         oneAdvMnt%den2_3d, oneAdvMnt%den3_3d, &
         oneAdvMnt%vc3d_in, oneAdvMnt%vc3d_out)

    ! dxtW, dytW, dztW from external fields

    call InitializeGridSpacings(&
         mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
         iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
         jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
         grid_g(ng)%dxt, &
         grid_g(ng)%dyt, &
         grid_g(ng)%fmapt, &
         grid_g(ng)%rtgt, &
         dztn(:,ng), &
         oneAdvMnt%dxtW, &
         oneAdvMnt%dytW, &
         oneAdvMnt%dztW)

    ! start ghost zone update for dxtW, dytW, dztW

    call PostSendRecvMsgs(oneGrid%AdvMntDxDySendX, oneGrid%AdvMntDxDyRecvX)

    !- get actual air densities, if using them instead of basic state fields

    if(use_true_density) then
       ! dd0_3d, dd0_3du, dd0_3dv, dd0_3dw, from external fields
       call GetTrueDensities(&
            mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
            iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
            jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
            level,&
            basic_g(ng)%rtp, &
            basic_g(ng)%rv, &
            basic_g(ng)%pp, &
            basic_g(ng)%pi0, &
            basic_g(ng)%theta, &
            oneAdvMnt%dd0_3d, &
            oneAdvMnt%dd0_3du, &
            oneAdvMnt%dd0_3dv, &
            oneAdvMnt%dd0_3dw)
       if (.not. theor_wind) then
          call PostSendRecvMsgs(oneGrid%AdvMntDd0SendX, oneGrid%AdvMntDd0RecvX)
       end if
    end if

    !- prepare wind velocities including map factors

    ndtZ=0
    ! u3d, v3d, w3d from external fields
    call PrepareWinds(&
         ng, mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
         iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
         jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
         dtlt, &
         basic_g(ng)%uc, basic_g(ng)%up, &
         basic_g(ng)%vc, basic_g(ng)%vp, &
         basic_g(ng)%wc, basic_g(ng)%wp, &
         grid_g(ng)%fmapui, grid_g(ng)%fmapvi, &
         grid_g(ng)%rtgt, grid_g(ng)%rtgu, grid_g(ng)%rtgv, &
         grid_g(ng)%f13t, grid_g(ng)%f23t, &
         oneAdvMnt%u3d, oneAdvMnt%v3d, oneAdvMnt%w3d, &
         aerosol, naer_transported, &
         dd_sedim, dzt, ndtZ)
    if (.not. theor_wind) then
       call PostSendRecvMsgs(oneGrid%AdvMntUVSendX, oneGrid%AdvMntUVRecvX)
    end if

    if(theor_wind) then
       ! dd0_3d, dd0_3du, dd0_3dv, dd0_3dw, u3d, v3d, w3d from external fields
       call PrepareTheorWinds(mzp, mxp, myp,&
            iOffset, i1ExternAtAdvMnt, iMxpExternAtAdvMnt,  &
            jOffset, j1ExternAtAdvMnt, jMypExternAtAdvMnt,  &
            dtlt, time,  &
            oneAdvMnt%u3d, oneAdvMnt%v3d, oneAdvMnt%w3d, &
            oneAdvMnt%dd0_3d, oneAdvMnt%dd0_3du, &
            oneAdvMnt%dd0_3dv, oneAdvMnt%dd0_3dw)
       call PostSendRecvMsgs(oneGrid%AdvMntDd0SendX, oneGrid%AdvMntDd0RecvX)
       call PostSendRecvMsgs(oneGrid%AdvMntUVSendX, oneGrid%AdvMntUVRecvX)
    end if

    !- prepare Walcek's air densities

    ! den0_3d, den1_3d, den2_3d, den3_3d from previously computed
    ! u3d, v3d, w3d, dxtW, dytW, dztW, dd0_3d, dd0_3du, dd0_3dv, dd0_3dw
    call GetWalceksDensities(&
         mzp, dtlt, mxpAdvMnt, mypAdvMnt, &
         i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
         j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
         oneAdvMnt%u3d, oneAdvMnt%v3d, oneAdvMnt%w3d, &
         oneAdvMnt%dd0_3d, oneAdvMnt%dd0_3du, &
         oneAdvMnt%dd0_3dv, oneAdvMnt%dd0_3dw, &
         oneAdvMnt%dxtW, oneAdvMnt%dytW, oneAdvMnt%dztW, &
         oneAdvMnt%den0_3d, oneAdvMnt%den1_3d, &
         oneAdvMnt%den2_3d, oneAdvMnt%den3_3d)
    call PostSendRecvMsgs(oneGrid%AdvMntDenSendX, oneGrid%AdvMntDenRecvX)

    ! message passing to update ghost zones

    call WaitSendRecvMsgs(oneGrid%AdvMntUVSendX, oneGrid%AdvMntUVRecvX)
    call WaitSendRecvMsgs(oneGrid%AdvMntDxDySendX, oneGrid%AdvMntDxDyRecvX)
    call WaitSendRecvMsgs(oneGrid%AdvMntDd0SendX, oneGrid%AdvMntDd0RecvX)
    call WaitSendRecvMsgs(oneGrid%AdvMntDenSendX, oneGrid%AdvMntDenRecvX)
    call PostSendRecvMsgs(oneGrid%AdvMntUVSendY, oneGrid%AdvMntUVRecvY)
    call PostSendRecvMsgs(oneGrid%AdvMntDxDySendY, oneGrid%AdvMntDxDyRecvY)
    call PostSendRecvMsgs(oneGrid%AdvMntDd0SendY, oneGrid%AdvMntDd0RecvY)
    call PostSendRecvMsgs(oneGrid%AdvMntDenSendY, oneGrid%AdvMntDenRecvY)
    call WaitSendRecvMsgs(oneGrid%AdvMntUVSendY, oneGrid%AdvMntUVRecvY)
    call WaitSendRecvMsgs(oneGrid%AdvMntDxDySendY, oneGrid%AdvMntDxDyRecvY)
    call WaitSendRecvMsgs(oneGrid%AdvMntDd0SendY, oneGrid%AdvMntDd0RecvY)
    call WaitSendRecvMsgs(oneGrid%AdvMntDenSendY, oneGrid%AdvMntDenRecvY)

    !- ready to do advection, loop over all scalars

    if(advmnt == 1) then
       i_scl=1                                            !- all scalars
    else if(advmnt == 2) then
       i_scl=num_scalar(ng) - NSPECIES_TRANSPORTED +1  !- only chemical + aer species
    else if(advmnt == 3) then
       i_scl=2                                            !- all scalars, but not theta_il
    end if

    ! copy external scalar_tab into oneGrid

    call DeepCopyToScalarTab(oneGrid%ScalarTab, oneGrid%ScalarTabSize)

    !srf- do n=1,num_scalar(ng)     ! original
    do n=i_scl,num_scalar(ng)

       !- if RK or ABM3 scheme, THP/THC are not transported here

       if (dyncore_flag == 2) then
          if (oneGrid%ScalarTab(n)%name == 'THC' .or. &
               oneGrid%ScalarTab(n)%name == 'THP') cycle
       end if

       !srf - somente para gases e aerossois
       !     do n=num_scalar(ng) - NSPECIES_TRANSPORTED +1,num_scalar(ng)
       !      if (scalar_tab(n,ng)%name /= 'COP' .and. scalar_tab(n,ng)%name /= 'CH4P') cycle
       !          scalar_tab(n,ng)%name /= 'O3P'  ) cycle

       !- Aerosol sedimentation
       IsThisScalarAer  = .false.
       current_aer_ispc = 0
       current_ndt_z    = 1
       if(ccatt == 1 .and. aerosol > 0 .and. n >= num_scalar_aer_1st) then
          !srf-  We are going to include sedimentation of aerosols at
          !      vertical advection tendency. It is supposed that scalars
          !      with  N >= num_scalar_aer_1st are _all_ aerosols .
          !
          IsThisScalarAer=.true.
          current_aer_ispc = n - num_scalar_aer_1st + 1
          current_ndt_z    = ndt_z (current_aer_ispc)

       end if

       if (associated(oneGrid%ScalarTab(n)%var_p_3D)) then

          ! set oneAdvMnt%vc3d_in north border to zero
          do j = 1, j1ExternAtAdvMnt-1
             do i = 1, mxpAdvMnt
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = 0.0
                end do
             end do
          end do
          do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
             jExtern = j + jOffset
             ! set oneAdvMnt%vc3d_in west border to zero
             do i = 1, i1ExternAtAdvMnt-1
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = 0.0
                end do
             end do
             ! copy scalartab external field to the
             ! inner part of oneAdvMnt%vc3d_in
             do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
                iExtern = i + iOffset
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = oneGrid%ScalarTab(n)%var_p_3D(k,iExtern,jExtern)
                end do
             end do
             ! set oneAdvMnt%vc3d_in east border to zero
             do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = 0.0
                end do
             end do
          end do
          ! set oneAdvMnt%vc3d_in south border to zero
          do j = jMypExternAtAdvMnt+1, mypAdvMnt
             do i = 1, mxpAdvMnt
                do k = 1, mzp
                   oneAdvMnt%vc3d_in(k,i,j) = 0.0
                end do
             end do
          end do

          call AdvectMnt(oneAdvMnt, oneGrid, &
               ngrid, mzp, mxp, myp, mxpAdvMnt, mypAdvMnt,&
               ia, iz, ja, jz, n, dtlt, &
               current_aer_ispc, current_ndt_z, IsThisScalarAer)


          call AdvectTendency(mzp, mxp, &
               iOffset, jOffset, &
               ia, iz, ja, jz, dtlt, &
               scalarp3D=oneGrid%ScalarTab(n)%var_p_3D, &
               AdvMntField=oneAdvMnt%vc3d_out, &
               scalart1D=oneGrid%ScalarTab(n)%var_t_1D)

       end if


    end do

    call DeepCopyFromScalarTab(oneGrid%ScalarTab, oneGrid%ScalarTabSize)

    ! destroy local memory area for large GhostZoneWidth variables

    call DestroyMonotonicAdvection(oneAdvMnt)

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine advmnt_driver





  subroutine InitializeGridSpacings(&
       mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
       iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
       dxt, dyt, fmapt, rtgt, dztn, &
       dxtW, dytW, dztW)

    ! computes new grid spacing on x, y and z
    ! for Wal2cek monotonic advection

    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, pointer, intent(in) :: dxt(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: dyt(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: fmapt(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rtgt(:,:)
    ! external field, pointer and values are intent(in)
    real, intent(in) :: dztn(:)
    ! external field
    real, pointer, intent(in) :: dxtW(:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dytW(:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dztW(:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)

    ! local var

    integer :: iExtern
    integer :: jExtern
    integer :: i
    integer :: j
    integer :: k
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(InitializeGridSpacings)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") mzp
       write(str(2),"(i8)") i1ExternAtAdvMnt
       write(str(3),"(i8)") iMxpExternAtAdvMnt
       write(str(4),"(i8)") j1ExternAtAdvMnt
       write(str(5),"(i8)") jMypExternAtAdvMnt
       call MsgDump(h//" set values of"//&
            " dxtW, dytW both at ("//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")"//&
            " and dztW at (1:"//trim(adjustl(str(1)))//"), nullifying remaining positions")
       write(str(1),"(i8)") iOffset
       write(str(2),"(i8)") jOffset
       call MsgDump(h//" external fields index = monotonic fields index + ("//&
            trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//")")
    end if

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          dxtW(i,j) = 0.0
          dytW(i,j) = 0.0
       end do
    end do

    ! fill Monotonic Advection fields only where
    ! external fields are in range
    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jExtern = j + jOffset
       ! set Monotonic Advection west ghost zone fields to zero 
       do i = 1, i1ExternAtAdvMnt - 1
          dxtW(i,j) = 0.0
          dytW(i,j) = 0.0
       end do
       ! fill where both Monotonic Advection and external fields
       ! are in range
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          rtgti = 1. / rtgt(iExtern,jExtern)

          !- at init/rams_grid.f90:
          !     dxt(i,j)=fmapt(i,j)/(xmn(i,ngrid)-xmn(i-1,ngrid))
          !     dyt(i,j)=fmapt(i,j)/(ymn(j,ngrid)-ymn(j-1,ngrid))

          dxtW(i,j) = 1./(dxt(iExtern,jExtern) * fmapt(iExtern,jExtern) * rtgti)
          dytW(i,j) = 1./(dyt(iExtern,jExtern) * fmapt(iExtern,jExtern) * rtgti)
       end do
       ! set Monotonic Advection east ghost zone fields to zero 
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          dxtW(i,j) = 0.0
          dytW(i,j) = 0.0
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          dxtW(i,j) = 0.0
          dytW(i,j) = 0.0
       end do
    end do

    ! z ranges of Monotonic Advection and external fields are identical
    do k = 1,mzp
       !- at init/gridset.f90:
       !  dztn(k,ifm) = 1. / (zmn(k,ifm) - zmn(k-1,ifm))
       ! Por que o Jacobiano nao depende de Z, o dztw depende somente
       ! de z.
       !dztW(k,i,j) = 1./ ( dzt(k) * rtgti * fmapt(i,j)**2 )
       dztW(k) = 1./ ( dztn(k) ) !
    end do
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine InitializeGridSpacings




  subroutine GetTrueDensities(&
       mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
       iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
       level,&
       rtp, rv, pp, pi0, theta, &
       dd0_3d, dd0_3du, dd0_3dv, dd0_3dw)

    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    integer, intent(in) :: level
    real, pointer, intent(in) :: rtp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rv(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: pp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: pi0(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: theta(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: dd0_3d(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dd0_3du(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dd0_3dv(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: dd0_3dw(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)

    ! local var
    integer :: iExtern
    integer :: jExtern
    integer :: i
    integer :: j
    integer :: k
    integer :: i1
    integer :: j1
    real :: c3

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(GetTrueDensities)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") lbound(dd0_3d,1)
       write(str(2),"(i8)") ubound(dd0_3d,1)
       write(str(3),"(i8)") lbound(dd0_3d,2)
       write(str(4),"(i8)") ubound(dd0_3d,2)
       write(str(5),"(i8)") lbound(dd0_3d,3)
       write(str(6),"(i8)") ubound(dd0_3d,3)
       call MsgDump(h//" dd0_3d dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") lbound(dd0_3du,1)
       write(str(2),"(i8)") ubound(dd0_3du,1)
       write(str(3),"(i8)") lbound(dd0_3du,2)
       write(str(4),"(i8)") ubound(dd0_3du,2)
       write(str(5),"(i8)") lbound(dd0_3du,3)
       write(str(6),"(i8)") ubound(dd0_3du,3)
       call MsgDump(h//" dd0_3du dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") lbound(dd0_3dv,1)
       write(str(2),"(i8)") ubound(dd0_3dv,1)
       write(str(3),"(i8)") lbound(dd0_3dv,2)
       write(str(4),"(i8)") ubound(dd0_3dv,2)
       write(str(5),"(i8)") lbound(dd0_3dv,3)
       write(str(6),"(i8)") ubound(dd0_3dv,3)
       call MsgDump(h//" dd0_3dv dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") lbound(dd0_3dw,1)
       write(str(2),"(i8)") ubound(dd0_3dw,1)
       write(str(3),"(i8)") lbound(dd0_3dw,2)
       write(str(4),"(i8)") ubound(dd0_3dw,2)
       write(str(5),"(i8)") lbound(dd0_3dw,3)
       write(str(6),"(i8)") ubound(dd0_3dw,3)
       call MsgDump(h//" dd0_3dw dimensioned ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//")")
       write(str(1),"(i8)") mzp
       write(str(2),"(i8)") i1ExternAtAdvMnt
       write(str(3),"(i8)") iMxpExternAtAdvMnt
       write(str(4),"(i8)") j1ExternAtAdvMnt
       write(str(5),"(i8)") jMypExternAtAdvMnt
       call MsgDump(h//" set values of"//&
            " dd0_3d, dd0_3du, dd0_3dv, dd0_3dw, all at ("//&
            "1:"//trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")"//&
            ", nullifying remaining positions")
       write(str(1),"(i8)") iOffset
       write(str(2),"(i8)") jOffset
       call MsgDump(h//" external fields index = monotonic fields index + ("//&
            trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//")")
    end if

    ! dd0_3d computation

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
          end do
       end do
    end do

    !- true air density at points "T"

    c3 = c2 * (cpi**c1)
    if( level == 0 ) then
       do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
          jExtern = j + jOffset
          ! set Monotonic Advection west ghost zone fields to zero
          do i = 1, i1ExternAtAdvMnt-1
             do k = 1, mzp
                dd0_3d(k,i,j) = 0.0
             end do
          end do
          ! fill where both Monotonic Advection and external fields
          ! are in range
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                dd0_3d(k,i,j) = (c3/theta(k,iExtern,jExtern))*&
                     (pi0(k,iExtern,jExtern)+pp(k,iExtern,jExtern))**c1
             end do
          end do
          ! set Monotonic Advection east ghost zone fields to zero
          do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
             do k = 1, mzp
                dd0_3d(k,i,j) = 0.0
             end do
          end do
       end do
    else
       do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
          jExtern = j + jOffset
          ! set Monotonic Advection west ghost zone fields to zero
          do i = 1, i1ExternAtAdvMnt-1
             do k = 1, mzp
                dd0_3d(k,i,j) = 0.0
             end do
          end do
          ! fill where both Monotonic Advection and external fields
          ! are in range
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                dd0_3d(k,i,j) = (c3/theta(k,iExtern,jExtern))*&
                     (1. + rtp(k,iExtern,jExtern))/ &
                     (1. + 1.61*rv(k,iExtern,jExtern))*&
                     (pi0(k,iExtern,jExtern)+pp(k,iExtern,jExtern))**c1
             end do
          end do
          ! set Monotonic Advection east ghost zone fields to zero
          do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
             do k = 1, mzp
                dd0_3d(k,i,j) = 0.0
             end do
          end do
       end do
    end if

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
          end do
       end do
    end do


    ! use dd0_3d to compute dd0_3du and dd0_3dv

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
          end do
       end do
    end do

    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       j1 = min(j+1,jMypExternAtAdvMnt)
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt-1
          do k = 1, mzp
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
          end do
       end do
       !- true air density computation
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          i1 = min(i+1,iMxpExternAtAdvMnt)
          do k = 1,mzp
             dd0_3du(k,i,j) = .5 * (dd0_3d(k,i,j) + dd0_3d(k,i1,j))
             dd0_3dv(k,i,j) = .5 * (dd0_3d(k,i,j) + dd0_3d(k,i,j1))
          end do
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
          end do
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
          end do
       end do
    end do

    ! compute true air density for w

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do

    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt-1
          do k = 1, mzp
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
       ! compute dd0_3dw 
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          do k = 1,mzp-1
             dd0_3dw(k,i,j) = 0.5*(dd0_3d(k,i,j) + dd0_3d(k+1,i,j))
          end do
          dd0_3dw(mzp,i,j)=dd0_3dw(mzp-1,i,j)
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine GetTrueDensities



  subroutine PrepareWinds(&
       ng, mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
       iOffset, i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
       dtlt, uc, up, vc, vp, wc, wp, &
       fmapui, fmapvi, rtgt, rtgu, rtgv, f13t, f23t, &
       u3d, v3d, w3d, &
       aerosol, naer_transported, dd_sedim, dzt, ndt_z)

    integer, intent(in) :: ng ! grid number, should dissapear
    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, intent(in) :: dtlt

    real, pointer, intent(in) :: uc(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: up(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: vc(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: vp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: wc(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: wp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: fmapui(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: fmapvi(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rtgt(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rtgu(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: rtgv(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: f13t(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: f23t(:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, intent(in) :: u3d(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: v3d(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, intent(in) :: w3d(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    integer, intent(in) :: aerosol
    ! flag for aerosol computation
    integer, intent(in) :: naer_transported
    ! # transported aerosol
    type(sedim_type), intent(in) :: dd_sedim(:,:)
    real, intent(in) :: dzt(:)
    integer, intent(inout) :: ndt_z(:)
    ! aerosol sedimentation timestep control

    !- local var

    integer :: jm
    integer :: jp
    integer :: im
    integer :: ip
    integer :: ispc
    integer :: i
    integer :: j
    integer :: k
    integer :: iExtern
    integer :: jExtern
    real :: cx1
    real :: cx2
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(PrepareWinds)**"

    if (dumpLocal) then
       call MsgDump(h//" starts; computes u3d, v3d and w3d just at a section"//&
            " restricted to the original ghost zone of 1")
    end if

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             w3d(k,i,j) = 0.0
             u3d(k,i,j) = 0.0
             v3d(k,i,j) = 0.0
          end do
       end do
    end do

    ! u3d, u3d, and w3d are input as the velocity components (averaged
    ! between past and current time levels) times dtlt.
    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jExtern = j + jOffset
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt-1
          do k = 1, mzp
             w3d(k,i,j) = 0.0
             u3d(k,i,j) = 0.0
             v3d(k,i,j) = 0.0
          end do
       end do
       ! u3d, u3d, and w3d are input as the velocity components (averaged
       ! between past and current time levels) times dtlt.
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          do k = 1, mzp
             w3d(k,i,j) = ( wc(k,iExtern,jExtern) + wp(k,iExtern,jExtern) )*0.5
             u3d(k,i,j) = ( uc(k,iExtern,jExtern) + up(k,iExtern,jExtern) )*0.5
             v3d(k,i,j) = ( vc(k,iExtern,jExtern) + vp(k,iExtern,jExtern) )*0.5
          end do
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             w3d(k,i,j) = 0.0
             u3d(k,i,j) = 0.0
             v3d(k,i,j) = 0.0
          end do
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             w3d(k,i,j) = 0.0
             u3d(k,i,j) = 0.0
             v3d(k,i,j) = 0.0
          end do
       end do
    end do

    ! transform w3d from cartesian vertical velocity to sigma_z velocity

    ! Add contribution to w3d from horiz winds crossing sloping sigma surfaces,
    ! and include 1/rtgt factor in w3d
    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jm = max(j1ExternAtAdvMnt,j-1)
       jp = min(jMypExternAtAdvMnt,j+1)
       jExtern = j + jOffset
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          im = max(i1ExternAtAdvMnt,i-1)
          ip = min(iMxpExternAtAdvMnt,i+1)
          rtgti = 1. / rtgt(iExtern,jExtern)
          do k = 1,mzp-1
             w3d(k,i,j) = &
                  ( &
                  (u3d(k,i,j)+u3d(k+1,i,j)+u3d(k,im,j)+u3d(k+1,im,j)) * f13t(iExtern,jExtern) + &
                  (v3d(k,i,j)+v3d(k+1,i,j)+v3d(k,i,jm)+v3d(k+1,i,jm)) * f23t(iExtern,jExtern)  &
                  ) * hw4(k) + w3d(k,i,j) * rtgti
          end do
       end do
    end do

    ! include map factors on u and v

    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jExtern = j + jOffset
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          cx1 = fmapui(iExtern,jExtern) * rtgu(iExtern,jExtern)
          cx2 = fmapvi(iExtern,jExtern) * rtgv(iExtern,jExtern)
          do k = 1,mzp-1
             u3d(k,i,j) = u3d(k,i,j) * cx1
             v3d(k,i,j) = v3d(k,i,j) * cx2
          end do
       end do
    end do
    !-----------------------------------------
    !- control for aerosol sedimentation
    if(aerosol > 0 .and. naer_transported > 0) then
       ! very crude estimation of CFL violation and fix for the number of sub-timesteps
       ! for large particles
       do ispc=1,naer_transported
          ndt_z(ispc)=ceiling(maxval(abs(dd_sedim(ispc,ng)%v_sed_part))*dtlt*maxval(dzt(1:mzp)))
       end do
    end if
    !- end of aerosol sedimentation
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine PrepareWinds





  subroutine GetWalceksDensities(&
       mzp, dtlt, mxpAdvMnt, mypAdvMnt, &
       i1ExternAtAdvMnt,  iMxpExternAtAdvMnt,  &
       j1ExternAtAdvMnt,  jMypExternAtAdvMnt,  &
       u3d, v3d, w3d, &
       dd0_3d, dd0_3du, dd0_3dv, dd0_3dw, &
       dxtW, dytW, dztW, &
       den0_3d, den1_3d, den2_3d, den3_3d)

    integer, intent(in) :: mzp
    real, intent(in) :: dtlt
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, pointer, intent(in) :: u3d(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: v3d(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: w3d(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dd0_3d(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dd0_3du(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dd0_3dv(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dd0_3dw(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dztW(:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dxtW(:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dytW(:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: den0_3d(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: den1_3d(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: den2_3d(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: den3_3d(:,:,:)
    ! pointer intent(in), values intent(out)


    ! local var
    integer i,j,k

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(GetWalceksDensities)**"

    if (dumpLocal) then
       call MsgDump(h//" starts; computes den0_3d, den1_3d, den2_3d and den3_3d"//&
            " just at a section restricted to the original ghost zone of 1")
    end if


    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             den0_3d(k,i,j) = 0.0
             den1_3d(k,i,j) = 0.0
             den2_3d(k,i,j) = 0.0
             den3_3d(k,i,j) = 0.0
          end do
       end do
    end do

    do  j = jMypExternAtAdvMnt, j1ExternAtAdvMnt+1, -1
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt
          do k = 1, mzp
             den0_3d(k,i,j) = 0.0
             den1_3d(k,i,j) = 0.0
             den2_3d(k,i,j) = 0.0
             den3_3d(k,i,j) = 0.0
          end do
       end do

       do  i = i1ExternAtAdvMnt+1, iMxpExternAtAdvMnt
          den0_3d(1,i,j) = 0.0
          den1_3d(1,i,j) = 0.0
          den2_3d(1,i,j) = 0.0
          den3_3d(1,i,j) = 0.0
          do k = 2, mzp
             den0_3d(k,i,j)=dd0_3d(k,i,j)
             den1_3d(k,i,j)=den0_3d(k,i,j)- dtlt/dxtW(i,j)*&
                  (dd0_3du(k,i,j)*u3d(k,i,j)-dd0_3du(k,i-1,j)*u3d(k,i-1,j))
             den2_3d(k,i,j)=den1_3d(k,i,j)- dtlt/dytW(i,j)*&
                  (dd0_3dv(k,i,j)*v3d(k,i,j)-dd0_3dv(k,i,j-1)*v3d(k,i,j-1))
             den3_3d(k,i,j)=den2_3d(k,i,j)- dtlt/dztW(k)  *&
                  (dd0_3dw(k,i,j)*w3d(k,i,j)-dd0_3dw(k-1,i,j)*w3d(k-1,i,j))
          end do
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             den0_3d(k,i,j) = 0.0
             den1_3d(k,i,j) = 0.0
             den2_3d(k,i,j) = 0.0
             den3_3d(k,i,j) = 0.0
          end do
       end do
    end do

    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             den0_3d(k,i,j) = 0.0
             den1_3d(k,i,j) = 0.0
             den2_3d(k,i,j) = 0.0
             den3_3d(k,i,j) = 0.0
          end do
       end do
    end do

    !srf- BC for den3_3d
    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       do k = 1, mzp
          den3_3d(k,i1ExternAtAdvMnt,j)=den3_3d(k,i1ExternAtAdvMnt+1,j)
       end do
    end do

    do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
       do k = 1, mzp
          den3_3d(k,i,j1ExternAtAdvMnt)=den3_3d(k,i,j1ExternAtAdvMnt+1)
       end do
    end do
  end subroutine GetWalceksDensities
  







  subroutine AdvectTendency(mzp, mxp, &
       iOffset, jOffset, &
       ia, iz, ja, jz, dtl, &
       scalarp3D, AdvMntField, scalart1D)
    integer, intent(in) :: mzp ! external field dimension
    integer, intent(in) :: mxp ! external field dimension
    integer, intent(in) :: iOffset
    integer, intent(in) :: jOffset
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(in) :: dtl
    real, pointer, intent(in) :: scalarp3D(:,:,:)
    ! external field, regular ghost zone, values are intent(in)
    real, pointer, intent(in) :: AdvMntField(:,:,:)
    ! monotonic advection field, extended ghost zone, values are intent(in)
    real, pointer, intent(in) :: scalart1D(:)
    ! 3D external field with regular ghost zone mapped into a 1D field;
    ! values are intent(inout)

    integer :: i
    integer :: iAdvMnt
    integer :: j
    integer :: jAdvMnt
    integer :: k
    integer :: ind1D
    real :: dtli

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(AdvectTendency)**"
    character(len=8) :: str(10)
    if (dumpLocal) then
       write(str(1),"(i8)") mzp-1
       write(str(2),"(i8)") ia
       write(str(3),"(i8)") iz
       write(str(4),"(i8)") ja
       write(str(5),"(i8)") jz
       call MsgDump(h//" at (2:"//trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
    end if

    dtli = 1. / dtl

    ! indexing external field

    do j = ja,jz
       jAdvMnt=j-jOffset
       do i = ia,iz
          iAdvMnt=i-iOffset
          ind1D=(i-1+(j-1)*mxp)*mzp
          do k = 2,mzp-1
             scalart1D(k+ind1D) = scalart1D(k+ind1D) + &
                  (AdvMntField(k,iAdvMnt,jAdvMnt)-scalarp3D(k,i,j)) * dtli
          end do
       end do
    end do
  end subroutine AdvectTendency







  subroutine AdvectMnt(oneAdvMnt, oneGrid, & 
       ngrid, mzp, mxp, myp, mxpAdvMnt, mypAdvMnt, &
       ia, iz, ja, jz, n, dt, &
       current_aer_ispc, current_ndt_z, IsThisScalarAer)

    type(MonotonicAdvection), pointer, intent(in) :: oneAdvMnt
    type(Grid), pointer, intent(in) :: oneGrid
    integer, intent(in) :: ngrid
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mxpAdvMnt
    integer, intent(in) :: mypAdvMnt
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: n
    real   , intent(in) :: dt
    integer, intent(in) :: current_ndt_z
    integer, intent(in) :: current_aer_ispc
    logical, intent(in) :: IsThisScalarAer

    !- local var

    integer itz
    integer ibegin,iend,jbegin,jend
    integer :: yStartSouth
    integer :: yEndSouth
    integer :: yStartNorth
    integer :: yEndNorth
    !- type of sedimentation scheme (0= Walcek, 1=upwind)
    integer , parameter :: iupwind = 0
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(AdvectMnt)**"
    character(len=8) :: str(10)

    iBegin = oneGrid%NodeDimsAdvMnt%ia-1
    iEnd   = oneGrid%NodeDimsAdvMnt%iz+1
    jBegin = oneGrid%NodeDimsAdvMnt%ja-1
    jEnd   = oneGrid%NodeDimsAdvMnt%jz+1

    !--- update X borders preparing X advection
    if (dumpLocal) then
       call MsgDump(h//" starts; update borders of vc3d_in for x advection")
    end if
    call PostSendRecvMsgs(oneGrid%AdvMntScaSendX, oneGrid%AdvMntScaRecvX)
    call WaitSendRecvMsgs(oneGrid%AdvMntScaSendX, oneGrid%AdvMntScaRecvX)
    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3DX to advect vc3d_in at y border regions")
    end if

    ! initialize vc3d_out
    
    oneAdvMnt%vc3d_out = oneAdvMnt%vc3d_in 

    ! do X-advection on Y border exchange regions
    
    yStartSouth = oneGrid%NodeDimsAdvMnt%ja
    yEndSouth = yStartSouth + oneGrid%NodeDimsAdvMnt%GhostZoneWidth - 1
    yEndNorth = oneGrid%NodeDimsAdvMnt%jz
    yStartNorth = yEndNorth - oneGrid%NodeDimsAdvMnt%GhostZoneWidth + 1

    call Advec3DX(mzp, mxpAdvMnt, mypAdvMnt, &
         yStart=yStartSouth, &
         yEnd=yEndSouth, &
         q0=oneAdvMnt%vc3d_in, &
         u=oneAdvMnt%u3d, &
         den0=oneAdvMnt%den0_3d, &
         den1=oneAdvMnt%den1_3d, &
         dt=dt, &
         dxx=oneAdvMnt%dxtW, &
         dd0=oneAdvMnt%dd0_3du, &
         qn=oneAdvMnt%vc3d_out)

    call Advec3DX(mzp, mxpAdvMnt, mypAdvMnt, &
         yStart=yStartNorth, &
         yEnd=yEndNorth, &
         q0=oneAdvMnt%vc3d_in, &
         u=oneAdvMnt%u3d, &
         den0=oneAdvMnt%den0_3d, &
         den1=oneAdvMnt%den1_3d, &
         dt=dt, &
         dxx=oneAdvMnt%dxtW, &
         dd0=oneAdvMnt%dd0_3du, &
         qn=oneAdvMnt%vc3d_out)

    ! post send/recv, buffering Y border exchange regions
    
    if (dumpLocal) then
       call MsgDump(h//" post send/recv on Y border regions")
    end if
    call PostSendRecvMsgs(oneGrid%AdvMntScaSendY, oneGrid%AdvMntScaRecvY)

    ! do X-advection on remaining Y internal regions 

    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3DX to advect vc3d_in at remaining y regions")
    end if
    call Advec3DX(mzp, mxpAdvMnt, mypAdvMnt, &
         yStart=yEndSouth+1, &
         yEnd=yStartNorth-1, &
         q0=oneAdvMnt%vc3d_in, &
         u=oneAdvMnt%u3d, &
         den0=oneAdvMnt%den0_3d, &
         den1=oneAdvMnt%den1_3d, &
         dt=dt, &
         dxx=oneAdvMnt%dxtW, &
         dd0=oneAdvMnt%dd0_3du, &
         qn=oneAdvMnt%vc3d_out)

    ! wait on message passing and debuffer Y border exchange regions

    if (dumpLocal) then
       call MsgDump(h//" wait on message exchange to update borders of vc3d_out for y advection")
    end if
    call WaitSendRecvMsgs(oneGrid%AdvMntScaSendY, oneGrid%AdvMntScaRecvY)

    !--- do Y-advection

    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3DY to advect vc3d_out, storing result in vc3d_in")
    end if
    call Advec3DY(mzp, mxpAdvMnt, mypAdvMnt, &
         q0=oneAdvMnt%vc3d_out, &
         u=oneAdvMnt%v3d, &
         den0=oneAdvMnt%den1_3d, &
         den1=oneAdvMnt%den2_3d, &
         dt=dt, &
         dxx=oneAdvMnt%dytW, &
         dd0=oneAdvMnt%dd0_3dv, &
         qn=oneAdvMnt%vc3d_in)

    !--- do k-advection
    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3DZ to advect vc3d_in, storing result in vc3d_out")
    end if
    call Advec3DZ(mzp, mxpAdvMnt, mypAdvMnt, &
         q0=oneAdvMnt%vc3d_in, &
         u=oneAdvMnt%w3d, &
         den0=oneAdvMnt%den2_3d, &
         den1=oneAdvMnt%den3_3d, &
         dt=dt, &
         dxx=oneAdvMnt%dztW, &
         dd0=oneAdvMnt%dd0_3dw, &
         qn=oneAdvMnt%vc3d_out)


    !- aerosol section to include sedimentation
    !- the sedimentation process is done using pure cartesian coordinates
    !- so, all sedimentation velocities are treat as cartesian vertical velocities
    !- which are positive downwards.
    if (dumpLocal) then
       write(str(1),"(i8)") aerosol
       write(str(2),"(l)") IsThisScalarAer
       call MsgDump(h//" aerosol="//trim(adjustl(str(1)))//&
            "; IsThisScalarAer="//trim(adjustl(str(2))))
    end if
    if(aerosol > 0 .and. IsThisScalarAer) then

       !-srf introducing a time-splitting for aerosol sedimentation

       if (dumpLocal) then
          write(str(1),"(i8)") iupwind
          call MsgDump(h//" iupwind="//trim(adjustl(str(1))))
       end if
       if(iupwind == 0 ) then
          ! - Walcek method
          ! this routine works _only_ for mass concentration or density (kg/m3)
          ! converting mixing ratio (kg/kg) to density (kg/m3)
          oneAdvMnt%vc3d_in(:,:,:)=oneAdvMnt%vc3d_out(:,:,:) * oneAdvMnt%den0_3d(:,:,:)

          !- do time splitting for aerosols with large fall velocities
          do itz=1,current_ndt_z
             call Advec3DZSedim(mzp,mxp,myp,&
                  ia,iz,ja,jz,                        &
                  q0=oneAdvMnt%vc3d_in(:,iBegin:iEnd,jBegin:jEnd),	 &
                  u=dd_sedim(current_aer_ispc,ngrid)%v_sed_part,          & !fall velocity
                  dt=dt/float(current_ndt_z),                              & !subtimestep
                  dzt=dzt(1:mzp), &
                  rtgt=grid_g(ngrid)%rtgt,	                 &
                  qn=oneAdvMnt%vc3d_out(:,iBegin:iEnd,jBegin:jEnd))

             ! copy output to input array for the next sup-timestep
             if(itz < current_ndt_z) oneAdvMnt%vc3d_in(:,:,:)=oneAdvMnt%vc3d_out(:,:,:)

          end do
          ! converting back mass concentration to mixing ratio
          oneAdvMnt%vc3d_out(:,:,:)=&
               oneAdvMnt%vc3d_out(:,:,:)/&
               oneAdvMnt%den0_3d(:,:,:)

       else if(iupwind == 1 ) then
          ! - upwind method
          !- do time splitting for aerosols with large fall velocities
          do itz=1,current_ndt_z
             call Advec3DZSedimUpw(mzp,mxp,myp,&
                  ia,iz,ja,jz, &
                  u=dd_sedim(current_aer_ispc,ngrid)%v_sed_part, & !fall velocity
                  dt=dt/float(current_ndt_z), & !subtimestep
                  dzt=dzt(1:mzp), &
                  rtgt=grid_g(ngrid)%rtgt,	&
                  qn=oneAdvMnt%vc3d_out(:,iBegin:iEnd,jBegin:jEnd))

          end do
       end if
    end if
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine AdvectMnt






  subroutine Advec3DX(mzp, mxpAdvMnt, mypAdvMnt, &
       yStart, yEnd, &
       q0, u, den0, den1, dt, dxx, dd0, &
       qn)
    integer, intent(in) :: mzp
    integer, intent(in) :: mxpAdvMnt
    integer, intent(in) :: mypAdvMnt
    integer, intent(in) :: yStart
    integer, intent(in) :: yEnd
    real, pointer, intent(in) :: q0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: u(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: den0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: den1(:,:,:)
    ! pointer and values intent(in)
    real, intent(in) :: dt
    real, pointer, intent(in) :: dxx(:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dd0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: qn(:,:,:)
    ! pointer intent(in), values intent(out)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(mzp,mxpAdvMnt,mypAdvMnt)
    real :: vcmax(mzp,mxpAdvMnt,mypAdvMnt)
    real :: vcmin(mzp,mxpAdvMnt,mypAdvMnt)
    logical :: imxmn(mzp,mxpAdvMnt,mypAdvMnt)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(Advec3DX)**"

    if (dumpLocal) then
       write(str(1),"(i8)") yStart
       write(str(2),"(i8)") yEnd
       call MsgDump(h//" starts computing at y=("//trim(adjustl(str(1)))//&
            ":"//trim(adjustl(str(2)))//")")
    end if

    ! verify y limits
    
    if (yStart < 2) then
       write(str(1),"(i8)") yStart
       call fatal_error(h//" yStart ("//trim(adjustl(str(1)))//&
            ") should be >= 2")
    else if (yEnd >= mypAdvMnt) then
       write(str(1),"(i8)") yEnd
       write(str(2),"(i8)") mypAdvMnt
       call fatal_error(h//" yEnd ("//trim(adjustl(str(1)))//&
            ") should be < mypAdvMnt ("//trim(adjustl(str(2)))//&
            ")")
    end if
    
    imxmn=.false.

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do j=yStart,yEnd
       do k=2,mzp-1
          if(u(k,1,j)>=zr0) flux(k,1,j)= q0(k,1,j)*u(k,1,j)*dt*dd0(k,1,j)
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=yStart,yEnd
       do  i=2,mxpAdvMnt-1
          do k=2,mzp-1
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i-1,j),q0(k,i+1,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k,i-1,j),q0(k,i+1,j))+eps)        !       extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k,i,j  )< zr0) ck1= q0(k,i+1,j)
             if(u(k,i-1,j)>=zr0) ck2= q0(k,i-1,j)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    do j=yStart,yEnd
       do  i=2,mxpAdvMnt-1 ! ia,iz-1 or 1,iz-1
          do k=2,mzp-1
             if(u(k,i,j)<zr0) cycle
             if(u(k,i-1,j)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
                x1n= (1.-x1)*(q0(k,i+1,j)-q0(k,i-1,j))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k,i-1,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i+1,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i+1,j))  ), max(q0(k,i,j),q0(k,i+1,j)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j) = max(vcmin(k,i,j),min(vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i-1,j)
             end if                                                  !Eq-9a
          end do
       end do
    end do

    ! If periodic boundary conditions are assumed, it is necessary
    !   to recalculate the updated mixing ratio at cell 1 if there
    !   is inflow to that cell from the boundary between IDIM and 1
    !   Here these statements are commented out, but should be uncommented
    !   if this subroutine is needed for periodic boundary conditions,
    !   and then one of the calling arguements to the subroutine is IPERIOD
    !   which is set to "1" if you assume period boundary conditions
    !      IF(IPERIOD==1) THEN
    !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
    !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
    !      END IF
    !
    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    !
    do j=yStart,yEnd
       do k=2,mzp-1
          if(u(k,mxpAdvMnt-1,j)<zr0) flux(k,mxpAdvMnt-1,j)= &
               q0(k,mxpAdvMnt,j)*u(k,mxpAdvMnt-1,j)*dt*dd0(k,mxpAdvMnt-1,j)
       end do
    end do

    do j=yStart,yEnd
       do i=mxpAdvMnt-1,2,-1 !iz,ia,-1
          do k=2,mzp-1
             if(u(k,i-1,j)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + &
                     flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k,i-1,j))/dxx(i,j)     ! Courant number
                x1n= (1.-x1)*(q0(k,i-1,j)-q0(k,i+1,j))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k,i+1,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i-1,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i-1,j)) ), max(q0(k,i,j),q0(k,i-1,j)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j)= max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i-1,j))/den1(k,i,j) ))
                flux(k,i-1,j)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop y-z
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DX
  

  


  subroutine Advec3DY(mzp, mxp, myp, &
       q0, u, den0, den1, dt, dxx, dd0, &
       qn)
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    real, pointer, intent(in) :: q0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: u(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: den0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: den1(:,:,:)
    ! pointer and values intent(in)
    real, intent(in) :: dt
    real, pointer, intent(in) :: dxx(:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dd0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(out):: qn(:,:,:)
    ! pointer intent(in), values intent(out)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(mzp,mxp,myp)
    real :: vcmax(mzp,mxp,myp)
    real :: vcmin(mzp,mxp,myp)
    logical :: imxmn(mzp,mxp,myp)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3DY)**"

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    ! copy input field to output field
    qn= q0
    imxmn=.false.

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do i=2,mxp-1
       do k=2,mzp-1
          if(u(k,i,1)>=zr0) flux(k,i,1)= q0(k,i,1)*u(k,i,1)*dt*dd0(k,i,1)
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !	mixing ratio at t+dt. If these limits are ever violated,
    !	non-monotonic (oscillatory) behavior in solution results
    do i=2,mxp-1
       do  j=2,myp-1 ! ja,jz
          do k=2,mzp-1
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i,j-1),q0(k,i,j+1))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k,i,j-1),q0(k,i,j+1))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k,i,j  )< zr0) ck1= q0(k,i,j+1)
             if(u(k,i,j-1)>=zr0) ck2= q0(k,i,j-1)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    do i=2,mxp-1
       do  j=2,myp-1 ! ja,jz
          do k=2,mzp-1
             if(u(k,i,j)<zr0) cycle
             if(u(k,i,j-1)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
                x1n= (1.-x1)*(q0(k,i,j+1)-q0(k,i,j-1))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k,i,j-1)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i,j+1)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i,j+1))  ), max(q0(k,i,j),q0(k,i,j+1)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j)= max(  vcmin(k,i,j),   min(   vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i,j-1) !Eq-9a
             end if
          end do
       end do
    end do

    ! If periodic boundary conditions are assumed, it is necessary
    !   to recalculate the updated mixing ratio at cell 1 if there
    !   is inflow to that cell from the boundary between IDIM and 1
    !   Here these statements are commented out, but should be uncommented
    !   if this subroutine is needed for periodic boundary conditions,
    !   and then one of the calling arguements to the subroutine is IPERIOD
    !   which is set to "1" if you assume period boundary conditions
    !      IF(IPERIOD==1) THEN
    !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
    !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
    !      END IF
    !
    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    !
    do i=2,mxp-1
       do k=2,mzp-1
          if(u(k,i,myp-1)<zr0) flux(k,i,myp-1)= &
               q0(k,i,myp)*u(k,i,myp-1)*dt*dd0(k,i,myp-1)
       end do
    end do

    do i=2,mxp-1
       do j=myp-1,2,-1 !jz,ja,-1
          do k=2,mzp-1
             if(u(k,i,j-1)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + &
                     flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k,i,j-1))/dxx(i,j)     ! Courant number
                x1n= (1.-x1)*(q0(k,i,j-1)-q0(k,i,j+1))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k,i,j+1)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i,j-1)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i,j-1)) ), max(q0(k,i,j),q0(k,i,j-1)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j)= max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i,j-1))/den1(k,i,j) ))
                flux(k,i,j-1)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop x-z
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DY






  subroutine Advec3DZ(mzp, mxp, myp, &
       q0,&
       u,&
       den0,&
       den1,&
       dt,&
       dxx,&
       dd0,&
       qn)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    real, pointer, intent(in) :: q0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: u(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: den0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: den1(:,:,:)
    ! pointer and values intent(in)
    real, intent(in) :: dt
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dxx(:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dd0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in):: qn(:,:,:)
    ! pointer intent(in), values intent(out)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(mzp,mxp,myp)
    real :: vcmax(mzp,mxp,myp)
    real :: vcmin(mzp,mxp,myp)
    logical :: imxmn(mzp,mxp,myp)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3DZ)**"

    if (dumpLocal) then
       call MsgDump(h//" starts ")
    end if

    ! copy input field to output field
    qn = q0
    imxmn=.false.


    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=2,myp-1
       do i=2,mxp-1
          do k=2,mzp-1 
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
             if(u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do


    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do j=2,myp-1
       do i=2,mxp-1
          if(u(1,i,j)>=zr0) flux(1,i,j)= &
               q0(1,i,j)*u(1,i,j)*dt*dd0(1,i,j)
          do k=2,mzp-1
             if(u(k,i  ,j)<zr0) cycle
             if(u(k-1,i,j)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(k)               ! Courant number
                x1n= (1.-x1)*(q0(k+1,i,j)-q0(k-1,i,j))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k-1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k+1,i,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k+1,i,j))  ), max(q0(k,i,j),q0(k+1,i,j)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j)= max(  vcmin(k,i,j),   min(   vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(k)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k-1,i,j)
             end if                                                  !Eq-9a
          end do
       end do
    end do

    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    do j=2,myp-1
       do i=2,mxp-1
          if(u(mzp-1,i,j)<zr0) flux(mzp-1,i,j)=&
               q0(mzp,i,j)*u(mzp-1,i,j)*dt*dd0(mzp-1,i,j)
          do k=mzp-1,2,-1
             if(u(k-1,i,j)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k) + flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k-1,i,j))/dxx(k)     ! Courant number
                x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k+1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k-1,i,j)) ), max(q0(k,i,j),q0(k-1,i,j)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j) = max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k)-x1*cf1*dd0(k-1,i,j))/den1(k,i,j) ))
                flux(k-1,i,j)=dxx(k)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop y-x
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DZ



  subroutine Advec3DZSedim(mzp, mxp, myp, &
       ia, iz, ja, jz,&
       q0,&
       u,&
       dt,&
       dzt,&
       rtgt,&
       qn)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real   , intent(in) :: q0(mzp,mxp,myp)
    real   , intent(in) :: u(mzp,mxp,myp)
    real   , intent(in) :: dt
    real   , intent(in) :: dzt(mzp)
    real   , intent(in) :: rtgt(mxp,myp)
    real   , intent(out):: qn(mzp,mxp,myp)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(mzp,mxp,myp)
    real :: vcmax(mzp,mxp,myp)
    real :: vcmin(mzp,mxp,myp)
    logical :: imxmn(mzp,mxp,myp)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3DZSedim)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") mzp
       write(str(2),"(i8)") mxp
       write(str(3),"(i8)") myp
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" starts at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    qn = q0
    imxmn=.false.

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=ja,jz
       do i=ia,iz
          do  k=2,mzp-1 !
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(-u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
             if(-u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
             if(k==2) ck2= q0(k,i,j) !for sedim only
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    do j=ja,jz
       do i=ia,iz
          rtgti=1./rtgt(i,j)
          flux(mzp-1,i,j)=q0(mzp,i,j)*(-u(mzp-1,i,j))*dt
          do k=mzp-1,2,-1
             !srf       x1=  dt*ABS(u(k-1,i,j))/dxx(k)     ! Courant number
             x1=  dt*abs(u(k-1,i,j))*dzt(k)*rtgti     ! Courant number
             if(k==2) x1 = 0. ! no flux below sfc terrain,for sedim only
             x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
             cf= q0(k,i,j) + x1n                                       !Eq-4b
             if(imxmn(k+1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
             if(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
             cf1= min( max( cf, min(q0(k,i,j),q0(k-1,i,j)) ), max(q0(k,i,j),q0(k-1,i,j)) )
             if(k>2) then  !for sedim only
                qn(k,i,j) = max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                                !srf                 (q0(k,i,j)-flux(k,i,j)/dxx(k)      -x1*cf1) ))
                     (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1) ))
             else
                qn(k,i,j) = (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1)
             end if
             !srf	   flux(k-1,i,j)=dxx(k)             *(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
             flux(k-1,i,j)=(1./(dzt(k)*rtgti))*(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
          end do
       end do
    end do !- big loop y-x
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DZSedim




  subroutine Advec3DZSedimUpw(mzp, mxp, myp, &
       ia, iz, ja, jz, &
       u, &
       dt, &
       dzt, &
       rtgt,&
       qn)

    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real   , intent(in) :: u(mzp,mxp,myp)
    real   , intent(in) :: dt
    real   , intent(in) :: dzt(mzp)
    real   , intent(in) :: rtgt(mxp,myp)
    real   , intent(out):: qn(mzp,mxp,myp)

    integer :: i
    integer :: j
    integer :: k
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3DZSedimUpw)**"

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    !- big loop y-x
    do j=ja,jz
       do i=ia,iz
          rtgti=1./rtgt(i,j)
          !srf dxx = dz = rtgti/dzt
          !srf qn(mzp-1,i,j) = qn(mzp-1,i,j) / (1.0 - dt*u(mzp-1,i,j)/dxx(mzp-1)      )
          qn(mzp-1,i,j) = qn(mzp-1,i,j) / (1.0 + dt*u(mzp-1,i,j)*dzt(mzp-1)*rtgti)
          do k=mzp-2,2,-1 !
             !srf    qn(k,i,j)= 1.0/(1.0+dt*u(k,i,j)/dxx(k))&
             !srf               *( qn(k,i,j)+ dt*u(k,i,j) /dxx(k+1) * qn(k+1,i,j) )
             qn(k,i,j)= 1.0/(1.0 + dt*u(k,i,j)*dzt(k)*rtgti)&
                  *( qn(k,i,j) + dt*u(k+1,i,j)*dzt(k+1)*rtgti * qn(k+1,i,j) )
             !   tc(i,j,l,k) = 1.0/(1.0+dt_settl(k)*vd_cor/delz(i,j,l2))&
             !  	 *(tc(i,j,l,k) + dt_settl(k)*vd_cor /delz(i,j,l2-1) &
             !  	 * tc(i,j,l+1,k))
          end do
       end do
    end do !- big loop y-x
    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DZSedimUpw








  subroutine InitializeDensities(mzp, mxp, myp, &
       mxpAdvMnt, mypAdvMnt, &
       iOffset, i1ExternAtAdvMnt, iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt, jMypExternAtAdvMnt,  &
       dn0, dn0u, dn0v, &
       dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw )

    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: mxpAdvMnt
    ! x dimension of Monotonic Advection fields
    integer, intent(in) :: mypAdvMnt
    ! y dimension of Monotonic Advection fields
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, pointer, intent(in) :: dn0(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dn0u(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dn0v(:,:,:)
    ! pointer and values intent(in)
    real, pointer, intent(in) :: dd0_3d(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3du(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3dv(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3dw(:,:,:)
    ! pointer intent(in), values intent(out)

    ! local var
    integer :: i
    integer :: j
    integer :: k
    integer :: iExtern
    integer :: jExtern

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(InitializeDensities)**"

    ! set Monotonic Advection south ghost zone fields to zero
    do j = 1, j1ExternAtAdvMnt-1
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do


    do j = j1ExternAtAdvMnt, jMypExternAtAdvMnt
       jExtern = j + jOffset
       ! set Monotonic Advection west ghost zone fields to zero
       do i = 1, i1ExternAtAdvMnt-1
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
       ! fill where both Monotonic Advection and external fields
       ! are in range
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          do k = 1, mzp
             dd0_3d (k,i,j) = dn0 (k,iExtern,jExtern)
             dd0_3du(k,i,j) = dn0u(k,iExtern,jExtern)
             dd0_3dv(k,i,j) = dn0v(k,iExtern,jExtern)
          end do
          do k = 1,mzp-1
             dd0_3dw(k,i,j) = 0.5*&
                  (dn0(k,iExtern,jExtern) + dn0(k+1,iExtern,jExtern))
          end do
          dd0_3dw(mzp,i,j)=dd0_3dw(mzp-1,i,j)
       end do
       ! set Monotonic Advection east ghost zone fields to zero
       do i = iMxpExternAtAdvMnt+1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do


    ! set Monotonic Advection north ghost zone fields to zero
    do j = jMypExternAtAdvMnt+1, mypAdvMnt
       do i = 1, mxpAdvMnt
          do k = 1, mzp
             dd0_3d(k,i,j) = 0.0
             dd0_3du(k,i,j) = 0.0
             dd0_3dv(k,i,j) = 0.0
             dd0_3dw(k,i,j) = 0.0
          end do
       end do
    end do
  end subroutine InitializeDensities




  subroutine PrepareTheorWinds(mzp, mxp, myp,&
       iOffset, i1ExternAtAdvMnt, iMxpExternAtAdvMnt,  &
       jOffset, j1ExternAtAdvMnt, jMypExternAtAdvMnt,  &
       dtlt, time,  &
       u3d, v3d, w3d, &
       dd0_3d, dd0_3du, dd0_3dv, dd0_3dw)

    integer, intent(in) :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer, intent(in) :: mxp
    ! x dimension of external fields 
    integer, intent(in) :: myp
    ! y dimension of external fields 
    integer, intent(in) :: iOffset
    ! x index offset from external to Monotonic Advection 
    integer, intent(in) :: i1ExternAtAdvMnt
    ! first x position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: iMxpExternAtAdvMnt
    ! last x position of external fields (mxp) indexed Monotonic Advection
    integer, intent(in) :: jOffset
    ! y index offset from external to Monotonic Advection 
    integer, intent(in) :: j1ExternAtAdvMnt
    ! first y position of external fields (1) indexed Monotonic Advection
    integer, intent(in) :: jMypExternAtAdvMnt
    ! last y position of external fields (myp) indexed Monotonic Advection
    real, intent(in) :: dtlt
    real, intent(in) :: time
    real, pointer, intent(in) :: u3d(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: v3d(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: w3d(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3d(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3du(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3dv(:,:,:)
    ! pointer intent(in), values intent(out)
    real, pointer, intent(in) :: dd0_3dw(:,:,:)
    ! pointer intent(in), values intent(out)

    !- local var
    integer :: i
    integer :: j
    integer :: k
    integer :: iExtern
    integer :: jExtern
    real :: dtlto2
    real :: ai0s  =  25.0
    real :: umx   =   80.0
    real, parameter :: pii   =   3.141592653589793
    real    :: xa,ilop,ya
    real    :: periodo  =   6.*3600.
    real, parameter :: iwndty = 2 

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(PrepareTheorWinds)**"

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    dtlto2 =  10.!*dtlt


    if(iwndty==1) ai0s= 50.5
    ilop= ai0s-21.  ! needed for printouts
    ! Define wind fields (rotation or divergent) and initial mixing ratios
    ! Cone at (25,50) for rotating winds; Cone at (50,50) divergent winds
    do j = jMypExternAtAdvMnt, j1ExternAtAdvMnt, -1
       jExtern = j + jOffset
       do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
          iExtern = i + iOffset
          do k = 1, mzp
             dd0_3d (k,i,j)=1.
             dd0_3du(k,i,j)=1.
             dd0_3dv(k,i,j)=1.
             dd0_3dw(k,i,j)=1.

             u3d(k,i,j)=-2.*umx*(real(jExtern)-real(110)/2.-.5)/real(110)
             v3d(k,i,j)= 2.*umx*(real(iExtern)-real(100)/2.-.5)/real(100)
             w3d(k,i,j)= 0.
          end do
       end do
    end do

    if(iwndty==1) then
       do j = jMypExternAtAdvMnt, j1ExternAtAdvMnt, -1
          jExtern = j + jOffset
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                xa=pii/25.
                u3d(k,i,j)=umx*&
                     sin(xa*real(iExtern))*&
                     sin(xa*(real(jExtern)))
                v3d(k,i,j)=umx*&
                     cos(xa*(real(iExtern)-.5))*&
                     cos(xa*(real(jExtern)+.5))
             end do
          end do
       end do

    else if (iwndty==2) then
       do j = jMypExternAtAdvMnt, j1ExternAtAdvMnt, -1
          jExtern = j + jOffset
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                xa=pii/100. ! myp=mxp
                u3d(k,i,j)=umx*&
                     (sin(xa*real(iExtern)))**2*&
                     sin(2*xa*(real(jExtern)))*&
                     cos(pii*time/periodo)
                v3d(k,i,j)=-umx*&
                     (sin(xa*real(jExtern)))**2*&
                     sin(2*xa*(real(iExtern)))*&
                     cos(pii*time/periodo)
             end do
          end do
       end do

    else if (iwndty==3) then
       do j = jMypExternAtAdvMnt, j1ExternAtAdvMnt, -1
          jExtern = j + jOffset
          do i = i1ExternAtAdvMnt, iMxpExternAtAdvMnt
             iExtern = i + iOffset
             do k = 1, mzp
                xa=pii/100. ! myp=mxp
                ya=50.
                u3d(k,i,j)=-umx*&
                     (sin(xa*real(iExtern)))**2*&
                     sin(2.*xa*(real(jExtern)-ya))*&
                     cos(pii*time/periodo)
                v3d(k,i,j)=0.5*&
                     umx*&
                     sin(2.*xa*real(iExtern))*&
                     cos(xa*(real(jExtern)-ya))*&
                     cos(pii*time/periodo)
             end do
          end do
       end do
    end if
  end subroutine PrepareTheorWinds





  subroutine StoreNamelistFileAtAdvMnt(oneNamelistFile)

    ! import NameListFile values into module variables

    type(namelistFile), pointer :: oneNamelistFile

    advmnt = oneNamelistFile%advmnt
    GhostZoneLength=oneNamelistFile%GhostZoneLength
  end subroutine StoreNamelistFileAtAdvMnt
end module ModMonotonicAdvection
