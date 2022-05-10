!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModBasicFields

  use ModNodeDimensions, only: &
       NodeDimensions

  use ModNamelistFile, only: &
       NamelistFile

  use ModParallelEnvironment, only: &
       MsgDump

  use iso_fortran_env, only: &
       real64

  use var_tables, only: &
       InsertVTab

  use mem_stilt, only: &
       iexev

  use mem_basic, only: &
       basic_g, &
       basicm_g
  
  implicit none

  private
  public :: BasicFields
  public :: CreateBasicFields
  public :: DestroyBasicFields
  public :: DumpBasicFields
  public :: InsertBasicFieldsAtVarTable 
  public :: DeepCopyToBasicFields
  public :: DeepCopyFromBasicFields
 
  type BasicFields

     ! Variables to be dimensioned by (nzp,nxp,nyp)

     real, pointer, contiguous :: up(:,:,:) => null()
     real, pointer, contiguous :: uc(:,:,:) => null()
     real, pointer, contiguous :: vp(:,:,:) => null()
     real, pointer, contiguous :: vc(:,:,:) => null()
     real, pointer, contiguous :: wp(:,:,:) => null()
     real, pointer, contiguous :: wc(:,:,:) => null()
     real, pointer, contiguous :: pp(:,:,:) => null()
     real, pointer, contiguous :: pc(:,:,:) => null()
     real, pointer, contiguous :: rv(:,:,:) => null()
     real, pointer, contiguous :: theta(:,:,:) => null()
     real, pointer, contiguous :: thp(:,:,:) => null()
     real, pointer, contiguous :: thc(:,:,:) => null()
     real, pointer, contiguous :: rtp(:,:,:) => null()
     real, pointer, contiguous :: pi0(:,:,:) => null()
     real, pointer, contiguous :: th0(:,:,:) => null()
     real, pointer, contiguous :: dn0(:,:,:) => null()
     real, pointer, contiguous :: dn0u(:,:,:) => null()
     real, pointer, contiguous :: dn0v(:,:,:) => null()

     ! Variables to be dimensioned by (nxp,nyp)

     real, pointer, contiguous :: fcoru(:,:) => null()
     real, pointer, contiguous :: fcorv(:,:) => null()
     real, pointer, contiguous :: cputime(:,:) => null()

  end type BasicFields

  logical :: updatingBasicFields=.false.
  
contains




  function CreateBasicFields(oneNodeDims, oneNamelistFile) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDims
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer :: res

    integer :: ierr
    integer :: mzp
    integer :: mxp
    integer :: myp
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateBasicFields)**"
    logical, parameter :: dumpLocal=.false.

    mzp=oneNodeDims%mzp
    mxp=oneNodeDims%mxp
    myp=oneNodeDims%myp

    if (dumpLocal) then
       write(str(1),"(i8)") mzp
       write(str(2),"(i8)") mxp
       write(str(3),"(i8)") myp
       call MsgDump(h//" starts with"//&
            " mzp="//trim(adjustl(str(1)))//&
            ", mxp="//trim(adjustl(str(2)))//&
            ", myp="//trim(adjustl(str(3))))
    end if

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if

    allocate(res%up(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate up fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%up=0.0

    allocate(res%uc(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate uc fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%uc=0.0

    allocate(res%vp(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate vp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%vp=0.0

    allocate(res%vc(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate vc fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%vc=0.0

    allocate(res%wp(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate wp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%wp=0.0

    allocate(res%wc(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate wc fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%wc=0.0

    allocate(res%pp(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate pp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%pp=0.0

    allocate(res%pc(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate pc fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%pc=0.0

    allocate(res%rv(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate rv fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%rv=0.0

    allocate(res%theta(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate theta fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%theta=0.0

    allocate(res%thp(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate thp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%thp=0.0

    if (oneNamelistFile%dyncore_flag == 2) then
       allocate(res%thc(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate thc fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%thc=0.0
    end if

    allocate(res%rtp(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate rtp fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%rtp=0.0

    allocate(res%pi0(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate pi0 fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%pi0=0.0

    allocate(res%th0(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate th0 fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%th0=0.0

    allocate(res%dn0(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate dn0 fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%dn0=0.0

    allocate(res%dn0u(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate dn0u fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%dn0u=0.0

    allocate(res%dn0v(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate dn0v fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%dn0v=0.0

    allocate(res%fcoru(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate fcoru fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%fcoru=0.0

    allocate(res%fcorv(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate fcorv fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%fcorv=0.0

    allocate(res%cputime(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate cputime fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%cputime=0.0

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end function CreateBasicFields





  subroutine DestroyBasicFields(oneBasicFields)
    type(BasicFields), pointer, intent(inout) :: oneBasicFields

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyBasicFields)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneBasicFields)) then
       if (dumpLocal) then
          call MsgDump(h//" null oneBasicFields")
       end if
       return
    end if

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    deallocate(oneBasicFields%up, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate up fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%uc, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate uc fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%vp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate vp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%vc, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate vc fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%wp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate wp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%wc, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate wc fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%pp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate pp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%pc, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate pc fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%rv, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate rv fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%theta, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate theta fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%thp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate thp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    if (associated(oneBasicFields%thc)) then
       deallocate(oneBasicFields%thc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate thc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    deallocate(oneBasicFields%rtp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate rtp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%pi0, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate pi0 fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%th0, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate th0 fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%dn0, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate dn0 fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%dn0u, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate dn0u fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%dn0v, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate dn0v fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%fcoru, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate fcoru fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%fcorv, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate fcorv fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields%cputime, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate cputime fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicFields, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneBasicFields fails with stat="//&
            trim(adjustl(str(1))))
    end if

    nullify(oneBasicFields)

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine DestroyBasicFields






  subroutine DumpBasicFields(oneBasicFields, name)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpBasicFields)**"

    if (associated(oneBasicFields)) then
       call MsgDump(h//" "//name//" associated")
    else
       call MsgDump(h//" "//name//" not associated")
    end if
  end subroutine DumpBasicFields


  


  subroutine InsertBasicFieldsAtVarTable(oneBasicFields, oneAveBasicFields, imean, gridId)
    type (BasicFields), pointer, intent(in) :: oneBasicFields
    type (BasicFields), pointer, intent(in) :: oneAveBasicFields
    integer, intent(in) :: imean
    integer, intent(in) :: gridId

    integer(kind=real64) :: npts
    character(len=*), parameter :: h="**(InsertBasicFieldsAtVarTable)**" 

    if (.not. associated(oneBasicFields)) then
       call fatal_error(h//" oneBasicFields not associated")
    else if (.not. associated(oneAveBasicFields)) then
       call fatal_error(h//" oneAveBasicFields not associated")
    end if

    ! Fill pointers to arrays into variable tables

    if (associated(oneBasicFields%up)) then
       npts=size(oneBasicFields%up)
       call InsertVTab (oneBasicFields%up,oneAveBasicFields%up  &
            ,gridId, npts, imean,  &
            'UP :3:hist:anal:mpti:mpt3:mpt2')
    else
       call fatal_error(h//" oneBasicFields%up not associated")
    end if

    if (associated(oneBasicFields%vp))  then
       call InsertVTab (oneBasicFields%vp,oneAveBasicFields%vp  &
            ,gridId, npts, imean,  &
            'VP :3:hist:anal:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicFields%wp))  then
       call InsertVTab (oneBasicFields%wp,oneAveBasicFields%wp  &
            ,gridId, npts, imean,  &
            'WP :3:hist:anal:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicFields%pp))  then
       call InsertVTab (oneBasicFields%pp,oneAveBasicFields%pp  &
            ,gridId, npts, imean,  &
            'PP :3:hist:anal:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicFields%uc))  then
       call InsertVTab (oneBasicFields%uc,oneAveBasicFields%uc  &
            ,gridId, npts, imean,  &
            'UC :3:hist:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicFields%vc))  then
       call InsertVTab (oneBasicFields%vc,oneAveBasicFields%vc  &
            ,gridId, npts, imean,  &
            'VC :3:hist:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicFields%wc))  then
       call InsertVTab (oneBasicFields%wc,oneAveBasicFields%wc  &
            ,gridId, npts, imean,  &
            'WC :3:hist:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicFields%pc))  then
       call InsertVTab (oneBasicFields%pc,oneAveBasicFields%pc  &
            ,gridId, npts, imean,  &
            'PC :3:hist:mpti:mpt3:mpt2')
    end if


    if (associated(oneBasicFields%thp)) then
       call InsertVTab (oneBasicFields%thp,oneAveBasicFields%thp  &
            ,gridId, npts, imean,  &
            'THP :3:hist:mpti:mpt3:mpt1')
    end if

    if (associated(oneBasicFields%rtp)) then
       call InsertVTab (oneBasicFields%rtp,oneAveBasicFields%rtp  &
            ,gridId, npts, imean,  &
            'RTP :3:hist:mpti:mpt3:mpt1')
    end if


    if(iexev == 2) then
       if (associated(oneBasicFields%theta)) then
          call InsertVTab (oneBasicFields%theta,oneAveBasicFields%theta  &
               ,gridId, npts, imean,  &
               'THETA :3:hist:anal:mpti:mpt3:mpt1')
       end if
    else
       if (associated(oneBasicFields%theta)) then
          call InsertVTab (oneBasicFields%theta,oneAveBasicFields%theta  &
               ,gridId, npts, imean,  &
               'THETA :3:hist:anal:mpti:mpt3')
       end if
    endif

    if (associated(oneBasicFields%thc)) then
       call InsertVTab (oneBasicFields%thc,oneAveBasicFields%thc  &
            ,gridId, npts, imean,  &
            'THC :3:hist:mpti:mpt3:mpt1')
    end if

    if(iexev == 2) then
       if (associated(oneBasicFields%rv)) then
          call InsertVTab (oneBasicFields%rv,oneAveBasicFields%rv  &
               ,gridId, npts, imean,  &
               'RV :3:hist:anal:mpti:mpt3:mpt1')
       end if
    else
       if (associated(oneBasicFields%rv)) then
          call InsertVTab (oneBasicFields%rv,oneAveBasicFields%rv  &
               ,gridId, npts, imean,  &
               'RV :3:hist:anal:mpti:mpt3')
       end if
    endif

    if (associated(oneBasicFields%pi0)) then
       call InsertVTab (oneBasicFields%pi0,oneAveBasicFields%pi0  &
            ,gridId, npts, imean,  &
            'PI0 :3:mpti')
    end if

    if (associated(oneBasicFields%th0)) then
       call InsertVTab (oneBasicFields%th0,oneAveBasicFields%th0  &
            ,gridId, npts, imean,  &
            'TH0 :3:mpti')
    end if

    if (associated(oneBasicFields%dn0)) then
       call InsertVTab (oneBasicFields%dn0,oneAveBasicFields%dn0  &
            ,gridId, npts, imean,  &
            'DN0 :3:mpti')
    end if

    if (associated(oneBasicFields%dn0u)) then
       call InsertVTab (oneBasicFields%dn0u,oneAveBasicFields%dn0u  &
            ,gridId, npts, imean,  &
            'DN0U :3:mpti')
    end if

    if (associated(oneBasicFields%dn0v)) then
       call InsertVTab (oneBasicFields%dn0v,oneAveBasicFields%dn0v  &
            ,gridId, npts, imean,  &
            'DN0V :3:mpti')
    end if


    if (associated(oneBasicFields%fcoru)) then
       npts = size(oneBasicFields%fcoru)
       call InsertVTab (oneBasicFields%fcoru,oneAveBasicFields%fcoru  &
            ,gridId, npts, imean,  &
            'FCORU :2:mpti')
    else
       call fatal_error(h//" oneBasicFields%fcoru not associated")
    end if

    if (associated(oneBasicFields%fcorv)) then
       call InsertVTab (oneBasicFields%fcorv,oneAveBasicFields%fcorv  &
            ,gridId, npts, imean,  &
            'FCORV :2:mpti')
    end if

    if (associated(oneBasicFields%cputime)) then
       call InsertVTab (oneBasicFields%cputime,oneAveBasicFields%cputime  &
            ,gridId, npts, imean,  &
            'CPUTIME :2:anal:mpti:mpt3')
    end if
  end subroutine InsertBasicFieldsAtVarTable



  

  subroutine DeepCopyToBasicFields(oneBasicFields, oneAveBasicFields)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(BasicFields), pointer, intent(in) :: oneAveBasicFields

    character(len=*), parameter :: h="**(DeepCopyToBasicFields)**"

    if (.not. associated(oneBasicFields)) then
       call fatal_error(h//" oneBasicFields not associated")
    else if (.not. associated(oneAveBasicFields)) then
       call fatal_error(h//" oneAveBasicFields not associated")
    end if

    if (updatingBasicFields) then
       call fatal_error(h//" invoked whenever BasicFields was beeing updated")
    else
       updatingBasicFields=.true.
    end if
    
    oneBasicFields%up => basic_g(1)%up
    oneAveBasicFields%up => basicm_g(1)%up
    
    oneBasicFields%uc => basic_g(1)%uc
    oneAveBasicFields%uc => basicm_g(1)%uc
    
    oneBasicFields%vp => basic_g(1)%vp
    oneAveBasicFields%vp => basicm_g(1)%vp
    
    oneBasicFields%vc => basic_g(1)%vc
    oneAveBasicFields%vc => basicm_g(1)%vc
    
    oneBasicFields%wp => basic_g(1)%wp
    oneAveBasicFields%wp => basicm_g(1)%wp
    
    oneBasicFields%wc => basic_g(1)%wc
    oneAveBasicFields%wc => basicm_g(1)%wc
    
    oneBasicFields%pp => basic_g(1)%pp
    oneAveBasicFields%pp => basicm_g(1)%pp
    
    oneBasicFields%pc => basic_g(1)%pc
    oneAveBasicFields%pc => basicm_g(1)%pc
    
    oneBasicFields%rv => basic_g(1)%rv
    oneAveBasicFields%rv => basicm_g(1)%rv
    
    oneBasicFields%theta => basic_g(1)%theta
    oneAveBasicFields%theta => basicm_g(1)%theta
    
    oneBasicFields%thp => basic_g(1)%thp
    oneAveBasicFields%thp => basicm_g(1)%thp
    
    oneBasicFields%thc => basic_g(1)%thc
    oneAveBasicFields%thc => basicm_g(1)%thc
    
    oneBasicFields%rtp => basic_g(1)%rtp
    oneAveBasicFields%rtp => basicm_g(1)%rtp
    
    oneBasicFields%pi0 => basic_g(1)%pi0
    oneAveBasicFields%pi0 => basicm_g(1)%pi0
    
    oneBasicFields%th0 => basic_g(1)%th0
    oneAveBasicFields%th0 => basicm_g(1)%th0
    
    oneBasicFields%dn0 => basic_g(1)%dn0
    oneAveBasicFields%dn0 => basicm_g(1)%dn0
    
    oneBasicFields%dn0u => basic_g(1)%dn0u
    oneAveBasicFields%dn0u => basicm_g(1)%dn0u
    
    oneBasicFields%dn0v => basic_g(1)%dn0v
    oneAveBasicFields%dn0v => basicm_g(1)%dn0v
    
    oneBasicFields%fcoru => basic_g(1)%fcoru
    oneAveBasicFields%fcoru => basicm_g(1)%fcoru
    
    oneBasicFields%fcorv => basic_g(1)%fcorv
    oneAveBasicFields%fcorv => basicm_g(1)%fcorv
    
    oneBasicFields%cputime => basic_g(1)%cputime
    oneAveBasicFields%cputime => basicm_g(1)%cputime
    
  end subroutine DeepCopyToBasicFields




  subroutine DeepCopyFromBasicFields(oneBasicFields, oneAveBasicFields)
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(BasicFields), pointer, intent(in) :: oneAveBasicFields

    character(len=*), parameter :: h="**(DeepCopyFromBasicFields)**"

    if (.not. associated(oneBasicFields)) then
       call fatal_error(h//" oneBasicFields not associated")
    else if (.not. associated(oneAveBasicFields)) then
       call fatal_error(h//" oneAveBasicFields not associated")
    end if

    if (.not. updatingBasicFields) then
       call fatal_error(h//" invoked whenever BasicFields was not beeing updated")
    else
       updatingBasicFields=.false.
    end if
    
    basic_g(1)%up => oneBasicFields%up
    basicm_g(1)%up => oneAveBasicFields%up

    basic_g(1)%uc => oneBasicFields%uc
    basicm_g(1)%uc => oneAveBasicFields%uc

    basic_g(1)%vp => oneBasicFields%vp
    basicm_g(1)%vp => oneAveBasicFields%vp

    basic_g(1)%vc => oneBasicFields%vc
    basicm_g(1)%vc => oneAveBasicFields%vc

    basic_g(1)%wp => oneBasicFields%wp
    basicm_g(1)%wp => oneAveBasicFields%wp

    basic_g(1)%wc => oneBasicFields%wc
    basicm_g(1)%wc => oneAveBasicFields%wc

    basic_g(1)%pp => oneBasicFields%pp
    basicm_g(1)%pp => oneAveBasicFields%pp

    basic_g(1)%pc => oneBasicFields%pc
    basicm_g(1)%pc => oneAveBasicFields%pc

    basic_g(1)%rv => oneBasicFields%rv
    basicm_g(1)%rv => oneAveBasicFields%rv

    basic_g(1)%theta => oneBasicFields%theta
    basicm_g(1)%theta => oneAveBasicFields%theta

    basic_g(1)%thp => oneBasicFields%thp
    basicm_g(1)%thp => oneAveBasicFields%thp

    basic_g(1)%thc => oneBasicFields%thc
    basicm_g(1)%thc => oneAveBasicFields%thc

    basic_g(1)%rtp => oneBasicFields%rtp
    basicm_g(1)%rtp => oneAveBasicFields%rtp

    basic_g(1)%pi0 => oneBasicFields%pi0
    basicm_g(1)%pi0 => oneAveBasicFields%pi0

    basic_g(1)%th0 => oneBasicFields%th0
    basicm_g(1)%th0 => oneAveBasicFields%th0

    basic_g(1)%dn0 => oneBasicFields%dn0
    basicm_g(1)%dn0 => oneAveBasicFields%dn0

    basic_g(1)%dn0u => oneBasicFields%dn0u
    basicm_g(1)%dn0u => oneAveBasicFields%dn0u

    basic_g(1)%dn0v => oneBasicFields%dn0v
    basicm_g(1)%dn0v => oneAveBasicFields%dn0v

    basic_g(1)%fcoru => oneBasicFields%fcoru
    basicm_g(1)%fcoru => oneAveBasicFields%fcoru

    basic_g(1)%fcorv => oneBasicFields%fcorv
    basicm_g(1)%fcorv => oneAveBasicFields%fcorv

    basic_g(1)%cputime => oneBasicFields%cputime
    basicm_g(1)%cputime => oneAveBasicFields%cputime

  end subroutine DeepCopyFromBasicFields
end module ModBasicFields
