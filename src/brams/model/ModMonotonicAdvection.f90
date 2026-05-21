!----------------------------------------------------------------------!
! Optional advection scheme for CCATT-BRAMS/BRAMS models version 4.2+  !
! Based on Walcek, 2000 (JGR) and Walcek and Aleksic, 1998 (ATENV).    !
! The scheme is highly conservative, monotonic and keeps mass mixing   !
! ratio positive definite. 					       !
! Implemented by Saulo Freitas (saulo.freitas@cptec.inpe.br) @ Jun/2009!
! MPI/Paralelized by L. Flavio/J. Panneta                              !
!----------------------------------------------------------------------!

module ModMonotonicAdvection

  use ModConvertDomainDecomp, only: &
       ConvertDomainDecomp, &
       SendRecvConvertDomainDecomp

  use ModMessageSet, only: &
       MessageSet, &
       PostSendRecvMsgs, &
       WaitSendRecvMsgs

  use ModScalarTable, only: &
       ScalarTable

  use ModBasicFields, only: &
       BasicFields

  use ModNodeDimensions, only: &
       NodeDimensions

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModDomainDecomp, only: &
       DomainDecomp

  use ModNamelistFile, only: &
       NamelistFile

  use ModMicControl, only: &
       MicControl

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

  use ccatt_start, only: &
       ccatt               ! (in)

  implicit none

  private

  type MonAdvBramsDD
     real, pointer, contiguous :: u3d(:,:,:)
     real, pointer, contiguous :: v3d(:,:,:)
     real, pointer, contiguous :: w3d(:,:,:)
     real, pointer, contiguous :: vc3d_in(:,:,:)
     real, pointer, contiguous :: vc3d_out(:,:,:)
     real, pointer, contiguous :: vc3d_x(:,:,:)
     real, pointer, contiguous :: vc3d_y(:,:,:)
     real, pointer, contiguous :: dd0_3d(:,:,:)
     real, pointer, contiguous :: dd0_3du(:,:,:)
     real, pointer, contiguous :: dd0_3dv(:,:,:)
     real, pointer, contiguous :: dd0_3dw(:,:,:)
     real, pointer, contiguous :: den0_3d(:,:,:)
     real, pointer, contiguous :: den1_3d(:,:,:)
     real, pointer, contiguous :: den2_3d(:,:,:)
     real, pointer, contiguous :: den3_3d(:,:,:)
     real, pointer, contiguous :: dxtW(:,:)
     real, pointer, contiguous :: dytW(:,:)
     real, pointer, contiguous :: dztW(:)
  end type MonAdvBramsDD


  type MonAdvOneDirDD
     real, pointer, contiguous :: u(:,:,:)
     real, pointer, contiguous :: q0(:,:,:)
     real, pointer, contiguous :: qn(:,:,:)
     real, pointer, contiguous :: dd0(:,:,:)
     real, pointer, contiguous :: den0(:,:,:)
     real, pointer, contiguous :: den1(:,:,:)
     real, pointer, contiguous :: dxx(:,:)
  end type MonAdvOneDirDD

  public :: advmnt_driver  ! Subroutine
  public :: StoreNamelistFileAtAdvMnt ! Subroutine

  ! public names, set by StoreNamelistFileAtAdvMnt
  integer, public :: advmnt 
  integer, public :: GhostZoneLength 

  ! module private variables

  ! constants
  real, parameter :: c1 = cv/rgas
  real, parameter :: c2 = p00/rgas

contains




  function CreateMonAdvBramsDD(oneNodeDimensions) result(oneMonAdvBramsDD)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(MonAdvBramsDD), pointer :: oneMonAdvBramsDD

    integer :: mzp
    integer :: mxp
    integer :: myp

    integer :: ierr
    character(len=512) :: message

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateMonAdvBramsDD)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneNodeDimensions)) then
       call fatal_error(h//" oneNodeDimensions not associated")
    end if

    mzp = oneNodeDimensions%mzp
    mxp = oneNodeDimensions%mxp
    myp = oneNodeDimensions%myp

    if (dumpLocal) then
       write(str(1),"(i8)") mzp
       write(str(2),"(i8)") mxp
       write(str(3),"(i8)") myp
       call MsgDump(h//" FullDirection="//oneNodeDimensions%FullDirection//" allocates "//&
            "all 3D at ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//"), 2D at ("//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//"), 1D at ("//&
            trim(adjustl(str(1)))//")")
    end if

    allocate(oneMonAdvBramsDD, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%u3d(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%u3d fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%v3d(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%v3d fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%w3d(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%w3d fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%dd0_3d(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%dd0_3d  fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%dd0_3du(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%dd0_3du fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%dd0_3dv(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%dd0_3dv fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%dd0_3dw(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%dd0_3dw fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%den0_3d(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%den0_3d fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%den1_3d(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%den1_3d fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%den2_3d(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%den2_3d fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%den3_3d(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%den3_3d fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%dxtW(mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%dxtW fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%dytW(mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%dytW fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%dztW(mzp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%dztW fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%vc3d_in(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%vc3d_in  fails with message "//trim(message))
    end if
    allocate(oneMonAdvBramsDD%vc3d_out(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvBramsDD%vc3d_out fails with message "//trim(message))
    end if
  end function CreateMonAdvBramsDD





  subroutine DestroyMonAdvBramsDD(oneMonAdvBramsDD)
    type(MonAdvBramsDD), pointer, intent(inout) :: oneMonAdvBramsDD

    integer :: ierr
    character(len=512) :: message

    logical :: dumpLocal=.false.
    character(len=*), parameter :: h="**(DestroyMonAdvBramsDD)**"

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    if (associated(oneMonAdvBramsDD)) then
       deallocate(oneMonAdvBramsDD%u3d, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%u3d fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%v3d, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%v3d fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%w3d, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%w3d fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%dd0_3d, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%dd0_3d  fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%dd0_3du, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%dd0_3du fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%dd0_3dv, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%dd0_3dv fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%dd0_3dw, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%dd0_3dw fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%den0_3d, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%den0_3d fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%den1_3d, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%den1_3d fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%den2_3d, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%den2_3d fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%den3_3d, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%den3_3d fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%dxtW, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%dxtW fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%dytW, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%dytW fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%dztW, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%dztW fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%vc3d_in, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%vc3d_in  fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD%vc3d_out, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD%vc3d_out fails with message "//trim(message))
       end if
       deallocate(oneMonAdvBramsDD, stat=ierr, errmsg=message)
       if (ierr /= 0) then
          call fatal_error(h//" deallocate oneMonAdvBramsDD fails with message "//trim(message))
       end if
    end if
    nullify(oneMonAdvBramsDD)

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine DestroyMonAdvBramsDD




  function CreateMonAdvOneDirDD(oneNodeDimensions, varName) result(oneMonAdvOneDirDD)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    character(len=*), intent(in) :: varName
    type(MonAdvOneDirDD), pointer :: oneMonAdvOneDirDD

    integer :: mzp
    integer :: mxp
    integer :: myp

    integer :: ierr
    character(len=512) :: message
    logical :: emptyFields

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateMonAdvOneDirDD)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneNodeDimensions)) then
       call fatal_error(h//" oneNodeDimensions not associated")
    end if

    ! in any case, create output variable of type MonAdvOneDirDD
    
    allocate(oneMonAdvOneDirDD, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneMonAdvOneDirDD fails with message "//trim(message))
    end if

    emptyFields = oneNodeDimensions%Empty

    ! treat empty and non-empty node dimensions
    
    if (emptyFields) then

       ! treat case output variable has empty fields
    
       nullify(oneMonAdvOneDirDD%u)
       nullify(oneMonAdvOneDirDD%q0)
       nullify(oneMonAdvOneDirDD%qn)
       nullify(oneMonAdvOneDirDD%dd0)
       nullify(oneMonAdvOneDirDD%den0)
       nullify(oneMonAdvOneDirDD%den1)
       nullify(oneMonAdvOneDirDD%dxx)
       if (dumpLocal) then
          call MsgDump(h//" all fields are empty for "//trim(varName))
       end if
    else

       ! treat case output variable has full fields
    
       mzp = oneNodeDimensions%mzp
       mxp = oneNodeDimensions%mxp
       myp = oneNodeDimensions%myp
       
       if (dumpLocal) then
          write(str(1),"(i8)") mzp
          write(str(2),"(i8)") mxp
          write(str(3),"(i8)") myp
          call MsgDump(h//" variable "//trim(varName)//": FullDirection="//oneNodeDimensions%FullDirection//" allocates "//&
               "all 3D at ("//&
               trim(adjustl(str(1)))//","//&
               trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//"), 2D at ("//&
               trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//"), 1D at ("//&
               trim(adjustl(str(1)))//")")
       end if
       
       allocate(oneMonAdvOneDirDD%u(mzp,mxp,myp), stat=ierr, errmsg=message)
       if (ierr/=0) then
          call fatal_error(h//" allocate u fails with message "//trim(message))
       end if
       allocate(oneMonAdvOneDirDD%q0(mzp,mxp,myp), stat=ierr, errmsg=message)
       if (ierr/=0) then
          call fatal_error(h//" allocate q0 fails with message "//trim(message))
       end if
       allocate(oneMonAdvOneDirDD%qn(mzp,mxp,myp), stat=ierr, errmsg=message)
       if (ierr/=0) then
          call fatal_error(h//" allocate qn fails with message "//trim(message))
       end if
       allocate(oneMonAdvOneDirDD%dd0(mzp,mxp,myp), stat=ierr, errmsg=message)
       if (ierr/=0) then
          call fatal_error(h//" allocate dd0 fails with message "//trim(message))
       end if
       allocate(oneMonAdvOneDirDD%den0(mzp,mxp,myp), stat=ierr, errmsg=message)
       if (ierr/=0) then
          call fatal_error(h//" allocate den0 fails with message "//trim(message))
       end if
       allocate(oneMonAdvOneDirDD%den1(mzp,mxp,myp), stat=ierr, errmsg=message)
       if (ierr/=0) then
          call fatal_error(h//" allocate den1 fails with message "//trim(message))
       end if
       allocate(oneMonAdvOneDirDD%dxx(mxp,myp), stat=ierr, errmsg=message)
       if (ierr/=0) then
          call fatal_error(h//" allocate dxx fails with message "//trim(message))
       end if
    end if
  end function CreateMonAdvOneDirDD





  subroutine DestroyMonAdvOneDirDD(oneMonAdvOneDirDD)
    type(MonAdvOneDirDD), pointer, intent(inout) :: oneMonAdvOneDirDD

    integer :: ierr
    character(len=512) :: message

    character(len=*), parameter :: h="**(DestroyMonAdvOneDirDD)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneMonAdvOneDirDD)) then
       return
    end if

    ! if one field is not associated, none are
    
    if (.not. associated(oneMonAdvOneDirDD%u)) then
       return
    end if

    
    deallocate(oneMonAdvOneDirDD%u, stat=ierr, errmsg=message)
    if (ierr/=0) then
       call fatal_error(h//" deallocate u fails with message "//trim(message))
    end if
    deallocate(oneMonAdvOneDirDD%q0, stat=ierr, errmsg=message)
    if (ierr/=0) then
       call fatal_error(h//" deallocate q0 fails with message "//trim(message))
    end if
    deallocate(oneMonAdvOneDirDD%qn, stat=ierr, errmsg=message)
    if (ierr/=0) then
       call fatal_error(h//" deallocate qn fails with message "//trim(message))
    end if
    deallocate(oneMonAdvOneDirDD%dd0, stat=ierr, errmsg=message)
    if (ierr/=0) then
       call fatal_error(h//" deallocate dd0 fails with message "//trim(message))
    end if
    deallocate(oneMonAdvOneDirDD%den0, stat=ierr, errmsg=message)
    if (ierr/=0) then
       call fatal_error(h//" deallocate den0 fails with message "//trim(message))
    end if
    deallocate(oneMonAdvOneDirDD%den1, stat=ierr, errmsg=message)
    if (ierr/=0) then
       call fatal_error(h//" deallocate den1 fails with message "//trim(message))
    end if
    deallocate(oneMonAdvOneDirDD%dxx, stat=ierr, errmsg=message)
    if (ierr/=0) then
       call fatal_error(h//" deallocate dxx fails with message "//trim(message))
    end if

    deallocate(oneMonAdvOneDirDD, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate oneMonAdvOneDirDD fails with message "//trim(message))
    end if
  end subroutine DestroyMonAdvOneDirDD





  subroutine advmnt_driver(&
       oneNodeDimensions, oneNodeDimensionsMonAdvX, oneNodeDimensionsMonAdvY, &
       oneBasicFields, &
       oneScalarTableSize, oneScalarTable, Id, varn, &
       oneMicVars, MonAdvInputSend, MonAdvInputRecv, &
       ConvertBramsToMonAdvX, ConvertBramsToMonAdvY, &
       ConvertMonAdvXToMonAdvY, ConvertBramsToBrams, &
       ConvertMonAdvYToBrams)

    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensionsMonAdvX
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensionsMonAdvY
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    integer, intent(in) :: oneScalarTableSize
    type(ScalarTable), pointer, intent(in) :: oneScalarTable(:)
    integer, intent(in) :: Id ! grid Id
    character(len=*),intent(in) :: varn
    type(MicControl), pointer, intent(in) :: oneMicVars
    type(MessageSet), pointer, intent(in) :: MonAdvInputSend
    type(MessageSet), pointer, intent(in) :: MonAdvInputRecv
    type(ConvertDomainDecomp), pointer, intent(in) :: ConvertBramsToMonAdvX
    type(ConvertDomainDecomp), pointer, intent(in) :: ConvertBramsToMonAdvY
    type(ConvertDomainDecomp), pointer, intent(in) :: ConvertMonAdvXToMonAdvY
    type(ConvertDomainDecomp), pointer, intent(in) :: ConvertBramsToBrams
    type(ConvertDomainDecomp), pointer, intent(in) :: ConvertMonAdvYToBrams

    !--- local vars
    integer :: scalarNbr
    integer :: ng
    integer :: i_scl
    integer :: current_aer_ispc
    integer :: current_ndt_z
    integer, target :: ndt_z(naer_transported)
    integer, target :: ndtZ(naer_transported)

    logical  :: IsThisScalarAer =.false.

    type(MonAdvBramsDD), pointer :: oneMonAdvBramsDD => null()
    type(MonAdvOneDirDD), pointer :: oneMonAdvXDD => null()
    type(MonAdvOneDirDD), pointer :: oneMonAdvYDD => null()

    integer :: mzp
    ! z dimension of external and Monotonic Advection fields 
    integer :: mxp
    ! x dimension of external fields 
    integer :: myp
    ! y dimension of external fields 
    character(len=64) :: scalarName

    integer :: ierr, x, y, z, ind
    character(len=512) :: message

    character(len=16) :: strL(11)
    character(len=*), parameter :: h="**(advmnt_driver)**"
    logical, parameter :: dumpLocal=.false.

    !- This scheme is not applied to advect  U, V, and W or shaved-eta

    if (varn .eq. 'V' .or. varn .eq. 'ALL') then
       call fatal_error(h//' not using mnt to advect u,v,w')
    end if
    if (if_adap /= 0) then
       call fatal_error(h//' MNT advection not ready for shaved eta')
    end if

    ! dimension of external fields (regular ghost zone width)

    mzp=oneNodeDimensions%mzp
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp

    ! current grid
    ! necessary while there are global variables outside oneGrid

    ng = Id

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    ! update borders of required input fields

    call PostSendRecvMsgs(MonAdvInputSend, MonAdvInputRecv)
    call WaitSendRecvMsgs(MonAdvInputSend, MonAdvInputRecv)

    ! create local memory area for Monotonic Advection variables

    if (dumpLocal) then
       call MsgDump(h//" call CreateMonAdvBramsDD")
    end if
    oneMonAdvBramsDD => CreateMonAdvBramsDD(oneNodeDimensions)

    if (dumpLocal) then
       call MsgDump(h//" call CreateMonAdvOneDirDD for X")
    end if
    oneMonAdvXDD => CreateMonAdvOneDirDD(oneNodeDimensionsMonAdvX,"oneMonAdvXDD")

    if (dumpLocal) then
       call MsgDump(h//" call CreateMonAdvOneDirDD for Y")
    end if
    oneMonAdvYDD => CreateMonAdvOneDirDD(oneNodeDimensionsMonAdvY,"oneMonAdvYDD")

    ! dxtW, dytW, dztW from external fields

    if (dumpLocal) then
       call MsgDump(h//" call InitializeGridSpacings")
    end if
    call InitializeGridSpacings(&
         oneNodeDimensions, &
         grid_g(ng)%dxt, &
         grid_g(ng)%dyt, &
         grid_g(ng)%fmapt, &
         grid_g(ng)%rtgt, &
         dztn(:,ng), &
         oneMonAdvBramsDD%dxtW, &
         oneMonAdvBramsDD%dytW, &
         oneMonAdvBramsDD%dztW)


    !- get actual air densities, if using them instead of basic state fields

    if (dumpLocal) then
       call MsgDump(h//" call GetTrueDensities")
    end if
    call GetTrueDensities(&
         oneNodeDimensions, &
         oneMicVars%level,&
         oneBasicFields%rtp, &
         oneBasicFields%rv, &
         oneBasicFields%pp, &
         oneBasicFields%pi0, &
         oneBasicFields%theta, &
         oneMonAdvBramsDD%dd0_3d, &
         oneMonAdvBramsDD%dd0_3du, &
         oneMonAdvBramsDD%dd0_3dv, &
         oneMonAdvBramsDD%dd0_3dw)

    ! update borders of dd0_X at Brams domain decomposition

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%dd0_3d, &
         nameSend="oneMonAdvBramsDD%dd0_3d", &
         fieldRecv=oneMonAdvBramsDD%dd0_3d, &
         nameRecv="oneMonAdvBramsDD%dd0_3d", &
         oneConvertDomainDecomp=ConvertBramsToBrams)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%dd0_3du, &
         nameSend="oneMonAdvBramsDD%dd0_3du", &
         fieldRecv=oneMonAdvBramsDD%dd0_3du, &
         nameRecv="oneMonAdvBramsDD%dd0_3du", &
         oneConvertDomainDecomp=ConvertBramsToBrams)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%dd0_3dv, &
         nameSend="oneMonAdvBramsDD%dd0_3dv", &
         fieldRecv=oneMonAdvBramsDD%dd0_3dv, &
         nameRecv="oneMonAdvBramsDD%dd0_3dv", &
         oneConvertDomainDecomp=ConvertBramsToBrams)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%dd0_3dw, &
         nameSend="oneMonAdvBramsDD%dd0_3dw", &
         fieldRecv=oneMonAdvBramsDD%dd0_3dw, &
         nameRecv="oneMonAdvBramsDD%dd0_3dw", &
         oneConvertDomainDecomp=ConvertBramsToBrams)

    !- prepare wind velocities including map factors

    if (dumpLocal) then
       call MsgDump(h//" call PrepareWinds")
    end if

    ndtZ=0
    call PrepareWinds(&
         ng, oneNodeDimensions, &
         dtlt, &
         oneBasicFields%uc, oneBasicFields%up, &
         oneBasicFields%vc, oneBasicFields%vp, &
         oneBasicFields%wc, oneBasicFields%wp, &
         grid_g(ng)%fmapui, grid_g(ng)%fmapvi, &
         grid_g(ng)%rtgt, grid_g(ng)%rtgu, grid_g(ng)%rtgv, &
         grid_g(ng)%f13t, grid_g(ng)%f23t, &
         oneMonAdvBramsDD%u3d, oneMonAdvBramsDD%v3d, oneMonAdvBramsDD%w3d, &
         aerosol, naer_transported, &
         dd_sedim, dzt, ndtZ)

    ! update borders of w3d at Brams domain decomposed fields

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%w3d, &
         nameSend="oneMonAdvBramsDD%w3d", &
         fieldRecv=oneMonAdvBramsDD%w3d, &
         nameRecv="oneMonAdvBramsDD%w3d", &
         oneConvertDomainDecomp=ConvertBramsToBrams)

    !- prepare Walcek's air densities

    if (dumpLocal) then
       call MsgDump(h//" call GetWalceksDensities")
    end if
    call GetWalceksDensities(&
         oneNodeDimensions, dtlt, &
         oneMonAdvBramsDD%u3d, oneMonAdvBramsDD%v3d, oneMonAdvBramsDD%w3d, &
         oneMonAdvBramsDD%dd0_3d, oneMonAdvBramsDD%dd0_3du, &
         oneMonAdvBramsDD%dd0_3dv, oneMonAdvBramsDD%dd0_3dw, &
         oneMonAdvBramsDD%dxtW, oneMonAdvBramsDD%dytW, oneMonAdvBramsDD%dztW, &
         oneMonAdvBramsDD%den0_3d, oneMonAdvBramsDD%den1_3d, &
         oneMonAdvBramsDD%den2_3d, oneMonAdvBramsDD%den3_3d)

    ! update borders of denX_3d at Brams domain decomposed fields

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%den0_3d, &
         nameSend="oneMonAdvBramsDD%den0_3d", &
         fieldRecv=oneMonAdvBramsDD%den0_3d, &
         nameRecv="oneMonAdvBramsDD%den0_3d", &
         oneConvertDomainDecomp=ConvertBramsToBrams)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%den1_3d, &
         nameSend="oneMonAdvBramsDD%den1_3d", &
         fieldRecv=oneMonAdvBramsDD%den1_3d, &
         nameRecv="oneMonAdvBramsDD%den1_3d", &
         oneConvertDomainDecomp=ConvertBramsToBrams)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%den2_3d, &
         nameSend="oneMonAdvBramsDD%den2_3d", &
         fieldRecv=oneMonAdvBramsDD%den2_3d, &
         nameRecv="oneMonAdvBramsDD%den2_3d", &
         oneConvertDomainDecomp=ConvertBramsToBrams)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%den3_3d, &
         nameSend="oneMonAdvBramsDD%den3_3d", &
         fieldRecv=oneMonAdvBramsDD%den3_3d, &
         nameRecv="oneMonAdvBramsDD%den3_3d", &
         oneConvertDomainDecomp=ConvertBramsToBrams)

    ! convert fields required by Advect3X from Brams domain decomposition
    ! to MonAdvX domain decomposition

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%u3d, &
         nameSend="oneMonAdvBramsDD%u3d", &
         fieldRecv=oneMonAdvXDD%u, &
         nameRecv="oneMonAdvXDD%u", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvX)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%dd0_3du, &
         nameSend="oneMonAdvBramsDD%dd0_3du", &
         fieldRecv=oneMonAdvXDD%dd0, &
         nameRecv="oneMonAdvXDD%dd0", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvX)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%den0_3d, &
         nameSend="oneMonAdvBramsDD%den0_3d", &
         fieldRecv=oneMonAdvXDD%den0, &
         nameRecv="oneMonAdvXDD%den0", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvX)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%den1_3d, &
         nameSend="oneMonAdvBramsDD%den1_3d", &
         fieldRecv=oneMonAdvXDD%den1, &
         nameRecv="oneMonAdvXDD%den1", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvX)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%dxtW, &
         nameSend="oneMonAdvBramsDD%dxtW", &
         fieldRecv=oneMonAdvXDD%dxx, &
         nameRecv="oneMonAdvXDD%dxx", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvX)

    ! convert fields required by Advect3Y from Brams domain decomposition
    ! to MonAdvY domain decomposition

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%v3d, &
         nameSend="oneMonAdvBramsDD%v3d", &
         fieldRecv=oneMonAdvYDD%u, &
         nameRecv="oneMonAdvYDD%u", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvY)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%dd0_3dv, &
         nameSend="oneMonAdvBramsDD%dd0_3dv", &
         fieldRecv=oneMonAdvYDD%dd0, &
         nameRecv="oneMonAdvYDD%dd0", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvY)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%den1_3d, &
         nameSend="oneMonAdvBramsDD%den1_3d", &
         fieldRecv=oneMonAdvYDD%den0, &
         nameRecv="oneMonAdvYDD%den0", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvY)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%den2_3d, &
         nameSend="oneMonAdvBramsDD%den2_3d", &
         fieldRecv=oneMonAdvYDD%den1, &
         nameRecv="oneMonAdvYDD%den1", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvY)

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvBramsDD%dytW, &
         nameSend="oneMonAdvBramsDD%dytW", &
         fieldRecv=oneMonAdvYDD%dxx, &
         nameRecv="oneMonAdvYDD%dxx", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvY)

    !- ready to do advection, loop over all scalars

    if(advmnt == 1) then
       i_scl=1                                            !- all scalars
    else if(advmnt == 2) then
       i_scl=oneScalarTableSize - NSPECIES_TRANSPORTED +1  !- only chemical + aer species
    else if(advmnt == 3) then
       i_scl=2                                            !- all scalars, but not theta_il
    end if

    !srf- do n=1,oneScalarTableSize     ! original
    SCALARS: do scalarNbr=i_scl,oneScalarTableSize

       !- if RK or ABM3 scheme, THP/THC are not transported here
       
       if (dyncore_flag == 2) then
          if (oneScalarTable(scalarNbr)%name == 'THC' .or. &
               oneScalarTable(scalarNbr)%name == 'THP') cycle
       end if
       scalarName=trim(oneScalarTable(scalarNbr)%name)

       !- Aerosol sedimentation
       IsThisScalarAer  = .false.
       current_aer_ispc = 0
       current_ndt_z    = 1
       if(ccatt == 1 .and. aerosol > 0 .and. scalarNbr >= num_scalar_aer_1st) then
          !srf-  We are going to include sedimentation of aerosols at
          !      vertical advection tendency. It is supposed that scalars
          !      with  N >= num_scalar_aer_1st are _all_ aerosols .
          !
          IsThisScalarAer=.true.
          current_aer_ispc = scalarNbr - num_scalar_aer_1st + 1
          current_ndt_z    = ndt_z (current_aer_ispc)
          
       end if

       if (associated(oneScalarTable(scalarNbr)%var_p_3D)) then

          if (dumpLocal) then
             write(strL(1),"(e15.7)") minval(oneScalarTable(scalarNbr)%var_p_3D)
             write(strL(2),"(e15.7)") maxval(oneScalarTable(scalarNbr)%var_p_3D)
             call MsgDump(h//" call AdvectMnt for variable "//trim(scalarName)//&
                  " range=["//trim(adjustl(strL(1)))//":"//trim(adjustl(strL(2)))//"]")
          end if
          call AdvectMnt(oneMonAdvBramsDD, oneMonAdvXDD, oneMonAdvYDD, &
               oneScalarTable, scalarNbr, dtlt, &
               current_aer_ispc, current_ndt_z, IsThisScalarAer, &
               oneNodeDimensions, oneNodeDimensionsMonAdvX, oneNodeDimensionsMonAdvY, &
               ConvertBramsToMonAdvX, ConvertMonAdvXToMonAdvY, ConvertMonAdvYToBrams)
          
          ! update borders of fully advected field at Brams domain decomposition
          
          call SendRecvConvertDomainDecomp(&
               fieldSend=oneMonAdvBramsDD%vc3d_out, &
               nameSend="oneMonAdvBramsDD%vc3d_out", &
               fieldRecv=oneMonAdvBramsDD%vc3d_out, &
               nameRecv="oneMonAdvBramsDD%vc3d_out", &
               oneConvertDomainDecomp=ConvertBramsToBrams)

          if (dumpLocal) then
             call MsgDump(h//" call AdvectTendency for variable "//trim(scalarName))
          end if
          
          call AdvectTendency(&
               scalarp3D=oneScalarTable(scalarNbr)%var_p_3D, &
               MonAdvField=oneMonAdvBramsDD%vc3d_out, &
               scalart1D=oneScalarTable(scalarNbr)%var_t_1D, &
               dtl=dtlt, &
               oneNodeDimensions=oneNodeDimensions)
          
       end if
    end do SCALARS

    ! destroy local memory area for large GhostZoneWidth variables

    if (dumpLocal) then
       call MsgDump(h//" call DestroyMonAdvBramsDD")
    end if
    call DestroyMonAdvBramsDD(oneMonAdvBramsDD)

    if (dumpLocal) then
       call MsgDump(h//" call DestroyMonAdvOneDirDD for X and Y")
    end if
    call DestroyMonAdvOneDirDD(oneMonAdvXDD)
    call DestroyMonAdvOneDirDD(oneMonAdvYDD)

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine advmnt_driver







  subroutine AdvectMnt(oneMonAdvBramsDD, oneMonAdvXDD, oneMonAdvYDD, &
       oneScalarTable, scalarNbr, dt, &
       current_aer_ispc, current_ndt_z, IsThisScalarAer, &
       oneNodeDimensions, oneNodeDimensionsMonAdvX, oneNodeDimensionsMonAdvY, &
       ConvertBramsToMonAdvX, ConvertMonAdvXToMonAdvY, ConvertMonAdvYToBrams)

    type(MonAdvBramsDD), pointer, intent(in) :: oneMonAdvBramsDD
    type(MonAdvOneDirDD), pointer, intent(inout) :: oneMonAdvXDD
    type(MonAdvOneDirDD), pointer, intent(inout) :: oneMonAdvYDD
    type(ScalarTable), pointer, intent(in) :: oneScalarTable(:)
    integer, intent(in) :: scalarNbr
    real   , intent(in) :: dt
    integer, intent(in) :: current_ndt_z
    integer, intent(in) :: current_aer_ispc
    logical, intent(in) :: IsThisScalarAer
    type(NodeDimensions), pointer, intent(in):: oneNodeDimensions
    type(NodeDimensions), pointer, intent(in):: oneNodeDimensionsMonAdvX
    type(NodeDimensions), pointer, intent(in):: oneNodeDimensionsMonAdvY
    type(ConvertDomainDecomp), pointer, intent(in) :: ConvertBramsToMonAdvX
    type(ConvertDomainDecomp), pointer, intent(in) :: ConvertMonAdvXToMonAdvY
    type(ConvertDomainDecomp), pointer, intent(in) :: ConvertMonAdvYToBrams

    integer :: mzp
    integer :: mxp
    integer :: myp
    logical :: emptyMonAdvX
    logical :: emptyMonAdvY
    
    !- type of sedimentation scheme (0= Walcek, 1=upwind)
    integer , parameter :: iupwind = 0
    character(len=128) :: fieldName

    character(len=*), parameter :: h="**(AdvectMnt)**"
    logical, parameter :: dumpLocal=.false.

    emptyMonAdvX = oneNodeDimensionsMonAdvX%Empty
    emptyMonAdvY = oneNodeDimensionsMonAdvY%Empty
    
    mzp=oneNodeDimensions%mzp
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp

    fieldName=oneScalarTable(scalarNbr)%name

    ! initialize vc3d_out for monotonic advection Z

    oneMonAdvBramsDD%vc3d_out = oneScalarTable(scalarNbr)%var_p_3D

    ! convert field to be advected from Brams domain decomposition
    ! to monotonic advection X domain decomposition

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneScalarTable(scalarNbr)%var_p_3D, &
         nameSend="oneScalarTable(scalarNbr)%var_p_3D", &
         fieldRecv=oneMonAdvXDD%q0, &
         nameRecv="oneMonAdvXDD%q0", &
         oneConvertDomainDecomp=ConvertBramsToMonAdvX)

    ! do X-advection 

    if (dumpLocal) then
       call MsgDump(h//" call Advec3DX")
    end if

    if (.not. emptyMonAdvX) then
       call Advec3DX(&
            q0=oneMonAdvXDD%q0, &
            u=oneMonAdvXDD%u, &
            den0=oneMonAdvXDD%den0, &
            den1=oneMonAdvXDD%den1, &
            dt=dt, &
            dxx=oneMonAdvXDD%dxx, &
            dd0=oneMonAdvXDD%dd0, &
            qn=oneMonAdvXDD%qn, &
            fieldName=trim(fieldName), &
            oneNodeDimensionsMonAdv=oneNodeDimensionsMonAdvX)
    end if

    ! convert field to be advected from monotonic advection X 
    ! domain decomposition to monotonic advection Y domain decomposition

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvXDD%qn, &
         nameSend="oneMonAdvXDD%qn", &
         fieldRecv=oneMonAdvYDD%q0, &
         nameRecv="oneMonAdvYDD%q0", &
         oneConvertDomainDecomp=ConvertMonAdvXToMonAdvY)

    !--- do Y-advection

    if (dumpLocal) then
       call MsgDump(h//" call Advec3DY")
    end if

    if (.not. emptyMonAdvY) then
       call Advec3DY(&
            q0=oneMonAdvYDD%q0, &
            u=oneMonAdvYDD%u, &
            den0=oneMonAdvYDD%den0, &
            den1=oneMonAdvYDD%den1, &
            dt=dt, &
            dxx=oneMonAdvYDD%dxx, &
            dd0=oneMonAdvYDD%dd0, &
            qn=oneMonAdvYDD%qn, &
            fieldName = trim(fieldName), &
            oneNodeDimensionsMonAdv=oneNodeDimensionsMonAdvY)
    end if

    ! convert field to be advected from monotonic advection Y
    ! domain decomposition to Brams domain decomposition

    call SendRecvConvertDomainDecomp(&
         fieldSend=oneMonAdvYDD%qn, &
         nameSend="oneMonAdvYDD%qn", &
         fieldRecv=oneMonAdvBramsDD%vc3d_in, &
         nameRecv="oneMonAdvBramsDD%vc3d_in", &
         oneConvertDomainDecomp=ConvertMonAdvYToBrams)

    !--- do k-advection
    if (dumpLocal) then
       call MsgDump(h//" call Advec3DZ")
    end if

    call Advec3DZ(&
         q0=oneMonAdvBramsDD%vc3d_in, &
         u=oneMonAdvBramsDD%w3d, &
         den0=oneMonAdvBramsDD%den2_3d, &
         den1=oneMonAdvBramsDD%den3_3d, &
         dt=dt, &
         dxx=oneMonAdvBramsDD%dztW, &
         dd0=oneMonAdvBramsDD%dd0_3dw, &
         qn=oneMonAdvBramsDD%vc3d_out, &
         fieldName = trim(fieldName), &
         oneNodeDimensionsMonAdv=oneNodeDimensions)


    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine AdvectMnt






  subroutine Advec3DX(&
       q0, u, den0, den1, dt, dxx, dd0, &
       qn, fieldName, oneNodeDimensionsMonAdv)
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
    character(len=*), intent(in) :: fieldName
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensionsMonAdv

    integer :: i
    integer :: j
    integer :: k
    integer :: yStart
    integer :: yEnd
    integer :: mzp
    integer :: mxp
    integer :: myp
    integer :: ierr
    real, allocatable :: flux(:,:,:)
    real, allocatable :: vcmax(:,:,:)
    real, allocatable :: vcmin(:,:,:)
    logical, allocatable :: imxmn(:,:,:)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    character(len=512) :: message

    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(Advec3DX)**"

    mzp = oneNodeDimensionsMonAdv%mzp
    mxp = oneNodeDimensionsMonAdv%mxp
    myp = oneNodeDimensionsMonAdv%myp

    allocate(flux(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate flux fails with message "//trim(message))
    end if

    allocate(vcmax(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate vcmax fails with message "//trim(message))
    end if

    allocate(vcmin(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate vcmin fails with message "//trim(message))
    end if

    allocate(imxmn(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate imxmn fails with message "//trim(message))
    end if

    ! avoid modifying Y boundary conditions

    if (oneNodeDimensionsMonAdv%borderSouth) then
       yStart=oneNodeDimensionsMonAdv%ja
    else
       yStart=1
    end if

    if (oneNodeDimensionsMonAdv%borderNorth) then
       yEnd=oneNodeDimensionsMonAdv%jz
    else
       yEnd=oneNodeDimensionsMonAdv%myp
    end if

    if (dumpLocal) then
       write(str(1),"(i8)") yStart
       write(str(2),"(i8)") yEnd
       call MsgDump(h//" starts computing at y=("//trim(adjustl(str(1)))//&
            ":"//trim(adjustl(str(2)))//")")
    end if

    imxmn=.false.
    flux = 0.0
    qn = q0

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
       do  i=2,mxp-1
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
       do  i=2,mxp-1 ! ia,iz-1 or 1,iz-1
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
          if(u(k,mxp-1,j)<zr0) flux(k,mxp-1,j)= &
               q0(k,mxp,j)*u(k,mxp-1,j)*dt*dd0(k,mxp-1,j)
       end do
    end do

    do j=yStart,yEnd
       do i=mxp-1,2,-1 !iz,ia,-1
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

    deallocate(flux, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate flux fails with message "//trim(message))
    end if

    deallocate(vcmax, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate vcmax fails with message "//trim(message))
    end if

    deallocate(vcmin, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate vcmin fails with message "//trim(message))
    end if

    deallocate(imxmn, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate imxmn fails with message "//trim(message))
    end if

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DX





  subroutine Advec3DY(&
       q0, u, den0, den1, dt, dxx, dd0, &
       qn, fieldName, oneNodeDimensionsMonAdv)
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
    real, pointer, intent(inout):: qn(:,:,:)
    ! pointer intent(in), values intent(out)
    character(len=*), intent(in) :: fieldName
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensionsMonAdv

    integer :: i
    integer :: j
    integer :: k
    integer :: mzp
    integer :: mxp
    integer :: myp
    integer :: ierr
    integer :: xStart
    integer :: xEnd
    real, allocatable :: flux(:,:,:)
    real, allocatable :: vcmax(:,:,:)
    real, allocatable :: vcmin(:,:,:)
    logical, allocatable :: imxmn(:,:,:)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    character(len=512) :: message

    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(Advec3DY)**"

    mzp = oneNodeDimensionsMonAdv%mzp
    mxp = oneNodeDimensionsMonAdv%mxp
    myp = oneNodeDimensionsMonAdv%myp

    allocate(flux(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate flux fails with message "//trim(message))
    end if

    allocate(vcmax(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate vcmax fails with message "//trim(message))
    end if

    allocate(vcmin(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate vcmin fails with message "//trim(message))
    end if

    allocate(imxmn(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate imxmn fails with message "//trim(message))
    end if

    if (oneNodeDimensionsMonAdv%borderWest) then
       xStart=2
    else
       xStart=1
    end if

    if (oneNodeDimensionsMonAdv%borderEast) then
       xEnd=mxp-1
    else
       xEnd=mxp
    end if

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    ! copy input field to output field

    qn= q0
    imxmn=.false.
    flux = 0.0

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain

    do i=xStart,xEnd
       do k=2,mzp-1
          if(u(k,i,1)>=zr0) flux(k,i,1)= q0(k,i,1)*u(k,i,1)*dt*dd0(k,i,1)
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !	mixing ratio at t+dt. If these limits are ever violated,
    !	non-monotonic (oscillatory) behavior in solution results

    do i=xStart,xEnd
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

    do i=xStart,xEnd
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

    do i=xStart,xEnd
       do k=2,mzp-1
          if(u(k,i,myp-1)<zr0) flux(k,i,myp-1)= &
               q0(k,i,myp)*u(k,i,myp-1)*dt*dd0(k,i,myp-1)
       end do
    end do

    do i=xStart,xEnd
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

    deallocate(flux, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate flux fails with message "//trim(message))
    end if

    deallocate(vcmax, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate vcmax fails with message "//trim(message))
    end if

    deallocate(vcmin, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate vcmin fails with message "//trim(message))
    end if

    deallocate(imxmn, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate imxmn fails with message "//trim(message))
    end if

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DY






  subroutine Advec3DZ(&
       q0,&
       u,&
       den0,&
       den1,&
       dt,&
       dxx,&
       dd0,&
       qn, &
       fieldName, &
       oneNodeDimensionsMonAdv)
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
    character(len=*), intent(in) :: fieldName
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensionsMonAdv

    integer :: i
    integer :: j
    integer :: k
    integer :: mzp
    integer :: mxp
    integer :: myp
    integer :: ierr
    real, allocatable :: flux(:,:,:)
    real, allocatable :: vcmax(:,:,:)
    real, allocatable :: vcmin(:,:,:)
    logical, allocatable :: imxmn(:,:,:)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    character(len=512) :: message

    logical, parameter :: dumpLocal=.false.
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(Advec3DZ)**"

    mzp = oneNodeDimensionsMonAdv%mzp
    mxp = oneNodeDimensionsMonAdv%mxp
    myp = oneNodeDimensionsMonAdv%myp

    allocate(flux(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate flux fails with message "//trim(message))
    end if

    allocate(vcmax(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate vcmax fails with message "//trim(message))
    end if

    allocate(vcmin(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate vcmin fails with message "//trim(message))
    end if

    allocate(imxmn(mzp,mxp,myp), stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" allocate imxmn fails with message "//trim(message))
    end if

    if (dumpLocal) then
       call MsgDump(h//" starts ")
    end if

    ! copy input field to output field

    qn = q0
    imxmn=.false.
    flux = 0.0

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

    deallocate(flux, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate flux fails with message "//trim(message))
    end if

    deallocate(vcmax, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate vcmax fails with message "//trim(message))
    end if

    deallocate(vcmin, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate vcmin fails with message "//trim(message))
    end if

    deallocate(imxmn, stat=ierr, errmsg=message)
    if (ierr /= 0) then
       call fatal_error(h//" deallocate imxmn fails with message "//trim(message))
    end if

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine Advec3DZ





  subroutine AdvectTendency(&
       scalarp3D, &
       MonAdvField, &
       scalart1D, &
       dtl, &
       oneNodeDimensions)
    real, pointer, intent(in) :: scalarp3D(:,:,:)
    ! external field, regular ghost zone, values are intent(in)
    real, pointer, intent(in) :: MonAdvField(:,:,:)
    ! monotonic advection field, extended ghost zone, values are intent(in)
    real, pointer, intent(in) :: scalart1D(:)
    ! 3D external field with regular ghost zone mapped into a 1D field;
    ! values are intent(inout)
    real, intent(in) :: dtl
    type(NodeDimensions), pointer, intent(in):: oneNodeDimensions

    integer :: i
    integer :: j
    integer :: k
    integer :: ind1D
    integer :: mzp
    integer :: mxp
    integer :: ia
    integer :: iz
    integer :: ja
    integer :: jz
    real :: dtli

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(AdvectTendency)**"
    character(len=8) :: str(10)

    mzp = oneNodeDimensions%mzp
    mxp = oneNodeDimensions%mxp
    ia = oneNodeDimensions%ia
    iz = oneNodeDimensions%iz
    ja = oneNodeDimensions%ja
    jz = oneNodeDimensions%jz

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
       do i = ia,iz
          ind1D=(i-1+(j-1)*mxp)*mzp
          do k = 2,mzp-1
             scalart1D(k+ind1D) = scalart1D(k+ind1D) + &
                  (MonAdvField(k,i,j)-scalarp3D(k,i,j)) * dtli
          end do
       end do
    end do
  end subroutine AdvectTendency






  subroutine InitializeGridSpacings(&
       oneNodeDimensions, &
       dxt, dyt, fmapt, rtgt, dztn, &
       dxtW, dytW, dztW)

    ! computes new grid spacing on x, y and z
    ! for Wal2cek monotonic advection

    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    ! pointer and values intent(in)
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

    integer :: i
    integer :: j
    integer :: k
    integer :: mzp
    integer :: mxp
    integer :: myp
    real :: rtgti

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InitializeGridSpacings)**"
    logical, parameter :: dumpLocal=.false.

    mzp = oneNodeDimensions%mzp
    mxp = oneNodeDimensions%mxp
    myp = oneNodeDimensions%myp

    if (dumpLocal) then
       write(str(1),"(i8)") mzp
       write(str(2),"(i8)") 1
       write(str(3),"(i8)") mxp
       write(str(4),"(i8)") 1
       write(str(5),"(i8)") myp
       call MsgDump(h//" set values of"//&
            " dxtW, dytW both at ("//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")"//&
            " and dztW at (1:"//trim(adjustl(str(1)))//"), nullifying remaining positions")
    end if

    ! fill Monotonic Advection fields 
    do j = 1, myp
       do i = 1, mxp
          rtgti = 1. / rtgt(i,j)

          !- at init/rams_grid.f90:
          !     dxt(i,j)=fmapt(i,j)/(xmn(i,ngrid)-xmn(i-1,ngrid))
          !     dyt(i,j)=fmapt(i,j)/(ymn(j,ngrid)-ymn(j-1,ngrid))

          dxtW(i,j) = 1./(dxt(i,j) * fmapt(i,j) * rtgti)
          dytW(i,j) = 1./(dyt(i,j) * fmapt(i,j) * rtgti)
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





  subroutine GetTrueDensities(oneNodeDimensions, &
       level, rtp, rv, pp, pi0, theta, &
       dd0_3d, dd0_3du, dd0_3dv, dd0_3dw)

    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    ! pointer and values are intent(in)
    integer, intent(in) :: level
    real, pointer, contiguous, intent(in) :: rtp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, contiguous, intent(in) :: rv(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, contiguous, intent(in) :: pp(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, contiguous, intent(in) :: pi0(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, contiguous, intent(in) :: theta(:,:,:)
    ! external field, pointer and values are intent(in)
    real, pointer, contiguous, intent(in) :: dd0_3d(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, contiguous, intent(in) :: dd0_3du(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, contiguous, intent(in) :: dd0_3dv(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)
    real, pointer, contiguous, intent(in) :: dd0_3dw(:,:,:)
    ! Monotonic Advection field, pointer is intent(in), values are intent(out)

    ! local var
    integer :: i
    integer :: j
    integer :: k
    integer :: mzp
    integer :: mxp
    integer :: myp
    integer :: i1
    integer :: j1
    real :: c3

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(GetTrueDensities)**"
    character(len=8) :: str(10)

    mzp = oneNodeDimensions%mzp
    mxp = oneNodeDimensions%mxp
    myp = oneNodeDimensions%myp

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
       write(str(2),"(i8)") 1
       write(str(3),"(i8)") mxp
       write(str(4),"(i8)") 1
       write(str(5),"(i8)") myp
       call MsgDump(h//" set values of"//&
            " dd0_3d, dd0_3du, dd0_3dv, dd0_3dw, all at ("//&
            "1:"//trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")")
    end if

    ! dd0_3d computation

    !- true air density at points "T"

    c3 = c2 * (cpi**c1)
    if( level == 0 ) then
       do j = 1, myp
          do i = 1, mxp
             do k = 1, mzp
                dd0_3d(k,i,j) = (c3/theta(k,i,j))*&
                     (pi0(k,i,j)+pp(k,i,j))**c1
             end do
          end do
       end do
    else
       do j = 1, myp
          do i = 1, mxp
             do k = 1, mzp
                dd0_3d(k,i,j) = (c3/theta(k,i,j))*&
                     (1. + rtp(k,i,j))/ &
                     (1. + 1.61*rv(k,i,j))*&
                     (pi0(k,i,j)+pp(k,i,j))**c1
             end do
          end do
       end do
    end if

    ! use dd0_3d to compute dd0_3du and dd0_3dv

    do j = 1, myp
       j1 = min(j+1,myp)
       !- true air density computation
       do i = 1, mxp
          i1 = min(i+1,mxp)
          do k = 1,mzp
             dd0_3du(k,i,j) = .5 * (dd0_3d(k,i,j) + dd0_3d(k,i1,j))
             dd0_3dv(k,i,j) = .5 * (dd0_3d(k,i,j) + dd0_3d(k,i,j1))
          end do
       end do
    end do

    ! use dd0_3d to compute dd0_3dw

    do j = 1, myp
       do i = 1, mxp
          do k = 1,mzp-1
             dd0_3dw(k,i,j) = 0.5*(dd0_3d(k,i,j) + dd0_3d(k+1,i,j))
          end do
          dd0_3dw(mzp,i,j)=dd0_3dw(mzp-1,i,j)
       end do
    end do

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine GetTrueDensities






  subroutine PrepareWinds(&
       ng, oneNodeDimensions, &
       dtlt, uc, up, vc, vp, wc, wp, &
       fmapui, fmapvi, rtgt, rtgu, rtgv, f13t, f23t, &
       u3d, v3d, w3d, &
       aerosol, naer_transported, dd_sedim, dzt, ndt_z)

    integer, intent(in) :: ng ! grid number, should dissapear
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
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
    integer :: mzp
    ! z dimension of Brams domain decomposed  fields 
    integer :: mxp
    ! x dimension of Brams domain decomposed  fields 
    integer :: myp
    ! y dimension of Brams domain decomposed fields 
    real :: cx1
    real :: cx2
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(PrepareWinds)**"

    ! dimension of Brams domain decomposed fields

    mzp=oneNodeDimensions%mzp
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp

    if (dumpLocal) then
       call MsgDump(h//" starts; computes u3d, v3d and w3d just at a section"//&
            " restricted to the original ghost zone of 1")
    end if

    ! u3d, u3d, and w3d are input as the velocity components (averaged
    ! between past and current time levels) times dtlt.
    do j = 1, myp
       ! u3d, u3d, and w3d are input as the velocity components (averaged
       ! between past and current time levels) times dtlt.
       do i = 1, mxp
          do k = 1, mzp
             w3d(k,i,j) = ( wc(k,i,j) + wp(k,i,j) )*0.5
             u3d(k,i,j) = ( uc(k,i,j) + up(k,i,j) )*0.5
             v3d(k,i,j) = ( vc(k,i,j) + vp(k,i,j) )*0.5
          end do
       end do
    end do

    ! transform w3d from cartesian vertical velocity to sigma_z velocity

    ! Add contribution to w3d from horiz winds crossing sloping sigma surfaces,
    ! and include 1/rtgt factor in w3d
    do j = 1, myp
       jm = max(1,j-1)
       jp = min(myp,j+1)
       do i = 1, mxp
          im = max(1,i-1)
          ip = min(mxp,i+1)
          rtgti = 1. / rtgt(i,j)
          do k = 1,mzp-1
             w3d(k,i,j) = &
                  ( &
                  (u3d(k,i,j)+u3d(k+1,i,j)+u3d(k,im,j)+u3d(k+1,im,j)) * f13t(i,j) + &
                  (v3d(k,i,j)+v3d(k+1,i,j)+v3d(k,i,jm)+v3d(k+1,i,jm)) * f23t(i,j)  &
                  ) * hw4(k) + w3d(k,i,j) * rtgti
          end do
       end do
    end do

    ! include map factors on u and v

    do j = 1, myp
       do i = 1, mxp
          cx1 = fmapui(i,j) * rtgu(i,j)
          cx2 = fmapvi(i,j) * rtgv(i,j)
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
       oneNodeDimensions, dtlt, &
       u3d, v3d, w3d, &
       dd0_3d, dd0_3du, dd0_3dv, dd0_3dw, &
       dxtW, dytW, dztW, &
       den0_3d, den1_3d, den2_3d, den3_3d)

    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    ! pointer and values intent(in)
    real, intent(in) :: dtlt
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
    integer :: i
    integer :: j
    integer :: k
    integer :: mzp
    integer :: mxp
    integer :: myp

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(GetWalceksDensities)**"

    if (dumpLocal) then
       call MsgDump(h//" starts; computes den0_3d, den1_3d, den2_3d and den3_3d"//&
            " just at a section restricted to the original ghost zone of 1")
    end if

    mzp = oneNodeDimensions%mzp
    mxp = oneNodeDimensions%mxp
    myp = oneNodeDimensions%myp

    ! set Monotonic Advection south ghost zone fields to zero
    do i = 1, mxp
       do k = 1, mzp
          den0_3d(k,i,1) = 0.0
          den1_3d(k,i,1) = 0.0
          den2_3d(k,i,1) = 0.0
          den3_3d(k,i,1) = 0.0
       end do
    end do

    do  j = myp, 2, -1
       ! set Monotonic Advection west ghost zone fields to zero
       do k = 1, mzp
          den0_3d(k,1,j) = 0.0
          den1_3d(k,1,j) = 0.0
          den2_3d(k,1,j) = 0.0
          den3_3d(k,1,j) = 0.0
       end do

       do  i = 2, mxp
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
    end do

    !srf- BC for den3_3d
    do j = 1, myp
       do k = 1, mzp
          den3_3d(k,1,j)=den3_3d(k,1+1,j)
       end do
    end do

    do i = 1, mxp
       do k = 1, mzp
          den3_3d(k,i,1)=den3_3d(k,i,1+1)
       end do
    end do
  end subroutine GetWalceksDensities





  subroutine StoreNamelistFileAtAdvMnt(oneNamelistFile)

    ! import NameListFile values into module variables

    type(namelistFile), pointer :: oneNamelistFile

    advmnt = oneNamelistFile%advmnt
    GhostZoneLength=oneNamelistFile%GhostZoneLength
  end subroutine StoreNamelistFileAtAdvMnt
end module ModMonotonicAdvection
