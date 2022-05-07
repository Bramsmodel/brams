!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module ModBasicVars

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
  public :: BasicVars
  public :: CreateBasicVars
  public :: DestroyBasicVars
  public :: DumpBasicVars
  public :: InsertBasicVarsAtVarTable 
  public :: DeepCopyToBasicVars
  public :: DeepCopyFromBasicVars
 
  type BasicVars

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

  end type BasicVars

contains




  function CreateBasicVars(oneNodeDims, oneNamelistFile) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDims
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicVars), pointer :: res

    integer :: ierr
    integer :: mzp
    integer :: mxp
    integer :: myp
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateBasicVars)**"
    logical, parameter :: dumpLocal=.true.

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
  end function CreateBasicVars





  subroutine DestroyBasicVars(oneBasicVars)
    type(BasicVars), pointer, intent(inout) :: oneBasicVars

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyBasicVars)**"
    logical, parameter :: dumpLocal=.true.

    if (.not. associated(oneBasicVars)) then
       if (dumpLocal) then
          call MsgDump(h//" null oneBasicVars")
       end if
       return
    end if

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    deallocate(oneBasicVars%up, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate up fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%uc, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate uc fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%vp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate vp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%vc, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate vc fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%wp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate wp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%wc, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate wc fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%pp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate pp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%pc, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate pc fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%rv, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate rv fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%theta, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate theta fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%thp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate thp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    if (associated(oneBasicVars%thc)) then
       deallocate(oneBasicVars%thc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate thc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    deallocate(oneBasicVars%rtp, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate rtp fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%pi0, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate pi0 fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%th0, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate th0 fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%dn0, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate dn0 fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%dn0u, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate dn0u fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%dn0v, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate dn0v fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%fcoru, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate fcoru fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%fcorv, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate fcorv fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars%cputime, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate cputime fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneBasicVars, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneBasicVars fails with stat="//&
            trim(adjustl(str(1))))
    end if

    nullify(oneBasicVars)

    if (dumpLocal) then
       call MsgDump(h//" finishes")
    end if
  end subroutine DestroyBasicVars






  subroutine DumpBasicVars(oneBasicVars, name)
    type(BasicVars), pointer, intent(in) :: oneBasicVars
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpBasicVars)**"

    if (associated(oneBasicVars)) then
       call MsgDump(h//" "//name//" associated")
    else
       call MsgDump(h//" "//name//" not associated")
    end if
  end subroutine DumpBasicVars


  


  subroutine InsertBasicVarsAtVarTable(oneBasicVars, oneAveBasicVars, imean, gridId)
    type (BasicVars), pointer, intent(in) :: oneBasicVars
    type (BasicVars), pointer, intent(in) :: oneAveBasicVars
    integer, intent(in) :: imean
    integer, intent(in) :: gridId

    integer(kind=real64) :: npts
    character(len=*), parameter :: h="**(InsertBasicVarsAtVarTable)**" 

    if (.not. associated(oneBasicVars)) then
       call fatal_error(h//" oneBasicVars not associated")
    else if (.not. associated(oneAveBasicVars)) then
       call fatal_error(h//" oneAveBasicVars not associated")
    end if

    ! Fill pointers to arrays into variable tables

    if (associated(oneBasicVars%up)) then
       npts=size(oneBasicVars%up)
       call InsertVTab (oneBasicVars%up,oneAveBasicVars%up  &
            ,gridId, npts, imean,  &
            'UP :3:hist:anal:mpti:mpt3:mpt2')
    else
       call fatal_error(h//" oneBasicVars%up not associated")
    end if

    if (associated(oneBasicVars%vp))  then
       call InsertVTab (oneBasicVars%vp,oneAveBasicVars%vp  &
            ,gridId, npts, imean,  &
            'VP :3:hist:anal:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicVars%wp))  then
       call InsertVTab (oneBasicVars%wp,oneAveBasicVars%wp  &
            ,gridId, npts, imean,  &
            'WP :3:hist:anal:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicVars%pp))  then
       call InsertVTab (oneBasicVars%pp,oneAveBasicVars%pp  &
            ,gridId, npts, imean,  &
            'PP :3:hist:anal:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicVars%uc))  then
       call InsertVTab (oneBasicVars%uc,oneAveBasicVars%uc  &
            ,gridId, npts, imean,  &
            'UC :3:hist:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicVars%vc))  then
       call InsertVTab (oneBasicVars%vc,oneAveBasicVars%vc  &
            ,gridId, npts, imean,  &
            'VC :3:hist:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicVars%wc))  then
       call InsertVTab (oneBasicVars%wc,oneAveBasicVars%wc  &
            ,gridId, npts, imean,  &
            'WC :3:hist:mpti:mpt3:mpt2')
    end if

    if (associated(oneBasicVars%pc))  then
       call InsertVTab (oneBasicVars%pc,oneAveBasicVars%pc  &
            ,gridId, npts, imean,  &
            'PC :3:hist:mpti:mpt3:mpt2')
    end if


    if (associated(oneBasicVars%thp)) then
       call InsertVTab (oneBasicVars%thp,oneAveBasicVars%thp  &
            ,gridId, npts, imean,  &
            'THP :3:hist:mpti:mpt3:mpt1')
    end if

    if (associated(oneBasicVars%rtp)) then
       call InsertVTab (oneBasicVars%rtp,oneAveBasicVars%rtp  &
            ,gridId, npts, imean,  &
            'RTP :3:hist:mpti:mpt3:mpt1')
    end if


    if(iexev == 2) then
       if (associated(oneBasicVars%theta)) then
          call InsertVTab (oneBasicVars%theta,oneAveBasicVars%theta  &
               ,gridId, npts, imean,  &
               'THETA :3:hist:anal:mpti:mpt3:mpt1')
       end if
    else
       if (associated(oneBasicVars%theta)) then
          call InsertVTab (oneBasicVars%theta,oneAveBasicVars%theta  &
               ,gridId, npts, imean,  &
               'THETA :3:hist:anal:mpti:mpt3')
       end if
    endif

    if (associated(oneBasicVars%thc)) then
       call InsertVTab (oneBasicVars%thc,oneAveBasicVars%thc  &
            ,gridId, npts, imean,  &
            'THC :3:hist:mpti:mpt3:mpt1')
    end if

    if(iexev == 2) then
       if (associated(oneBasicVars%rv)) then
          call InsertVTab (oneBasicVars%rv,oneAveBasicVars%rv  &
               ,gridId, npts, imean,  &
               'RV :3:hist:anal:mpti:mpt3:mpt1')
       end if
    else
       if (associated(oneBasicVars%rv)) then
          call InsertVTab (oneBasicVars%rv,oneAveBasicVars%rv  &
               ,gridId, npts, imean,  &
               'RV :3:hist:anal:mpti:mpt3')
       end if
    endif

    if (associated(oneBasicVars%pi0)) then
       call InsertVTab (oneBasicVars%pi0,oneAveBasicVars%pi0  &
            ,gridId, npts, imean,  &
            'PI0 :3:mpti')
    end if

    if (associated(oneBasicVars%th0)) then
       call InsertVTab (oneBasicVars%th0,oneAveBasicVars%th0  &
            ,gridId, npts, imean,  &
            'TH0 :3:mpti')
    end if

    if (associated(oneBasicVars%dn0)) then
       call InsertVTab (oneBasicVars%dn0,oneAveBasicVars%dn0  &
            ,gridId, npts, imean,  &
            'DN0 :3:mpti')
    end if

    if (associated(oneBasicVars%dn0u)) then
       call InsertVTab (oneBasicVars%dn0u,oneAveBasicVars%dn0u  &
            ,gridId, npts, imean,  &
            'DN0U :3:mpti')
    end if

    if (associated(oneBasicVars%dn0v)) then
       call InsertVTab (oneBasicVars%dn0v,oneAveBasicVars%dn0v  &
            ,gridId, npts, imean,  &
            'DN0V :3:mpti')
    end if


    if (associated(oneBasicVars%fcoru)) then
       npts = size(oneBasicVars%fcoru)
       call InsertVTab (oneBasicVars%fcoru,oneAveBasicVars%fcoru  &
            ,gridId, npts, imean,  &
            'FCORU :2:mpti')
    else
       call fatal_error(h//" oneBasicVars%fcoru not associated")
    end if

    if (associated(oneBasicVars%fcorv)) then
       call InsertVTab (oneBasicVars%fcorv,oneAveBasicVars%fcorv  &
            ,gridId, npts, imean,  &
            'FCORV :2:mpti')
    end if

    if (associated(oneBasicVars%cputime)) then
       call InsertVTab (oneBasicVars%cputime,oneAveBasicVars%cputime  &
            ,gridId, npts, imean,  &
            'CPUTIME :2:anal:mpti:mpt3')
    end if
  end subroutine InsertBasicVarsAtVarTable



  

  subroutine DeepCopyToBasicVars(oneBasicVars, oneAveBasicVars)
    type(BasicVars), pointer, intent(in) :: oneBasicVars
    type(BasicVars), pointer, intent(in) :: oneAveBasicVars

    character(len=*), parameter :: h="**(DeepCopyToBasicVars)**"

    if (.not. associated(oneBasicVars)) then
       call fatal_error(h//" oneBasicVars not associated")
    else if (.not. associated(oneAveBasicVars)) then
       call fatal_error(h//" oneAveBasicVars not associated")
    end if
    
    oneBasicVars%up => basic_g(1)%up
    oneAveBasicVars%up => basicm_g(1)%up
    
    oneBasicVars%uc => basic_g(1)%uc
    oneAveBasicVars%uc => basicm_g(1)%uc
    
    oneBasicVars%vp => basic_g(1)%vp
    oneAveBasicVars%vp => basicm_g(1)%vp
    
    oneBasicVars%vc => basic_g(1)%vc
    oneAveBasicVars%vc => basicm_g(1)%vc
    
    oneBasicVars%wp => basic_g(1)%wp
    oneAveBasicVars%wp => basicm_g(1)%wp
    
    oneBasicVars%wc => basic_g(1)%wc
    oneAveBasicVars%wc => basicm_g(1)%wc
    
    oneBasicVars%pp => basic_g(1)%pp
    oneAveBasicVars%pp => basicm_g(1)%pp
    
    oneBasicVars%pc => basic_g(1)%pc
    oneAveBasicVars%pc => basicm_g(1)%pc
    
    oneBasicVars%rv => basic_g(1)%rv
    oneAveBasicVars%rv => basicm_g(1)%rv
    
    oneBasicVars%theta => basic_g(1)%theta
    oneAveBasicVars%theta => basicm_g(1)%theta
    
    oneBasicVars%thp => basic_g(1)%thp
    oneAveBasicVars%thp => basicm_g(1)%thp
    
    oneBasicVars%thc => basic_g(1)%thc
    oneAveBasicVars%thc => basicm_g(1)%thc
    
    oneBasicVars%rtp => basic_g(1)%rtp
    oneAveBasicVars%rtp => basicm_g(1)%rtp
    
    oneBasicVars%pi0 => basic_g(1)%pi0
    oneAveBasicVars%pi0 => basicm_g(1)%pi0
    
    oneBasicVars%th0 => basic_g(1)%th0
    oneAveBasicVars%th0 => basicm_g(1)%th0
    
    oneBasicVars%dn0 => basic_g(1)%dn0
    oneAveBasicVars%dn0 => basicm_g(1)%dn0
    
    oneBasicVars%dn0u => basic_g(1)%dn0u
    oneAveBasicVars%dn0u => basicm_g(1)%dn0u
    
    oneBasicVars%dn0v => basic_g(1)%dn0v
    oneAveBasicVars%dn0v => basicm_g(1)%dn0v
    
    oneBasicVars%fcoru => basic_g(1)%fcoru
    oneAveBasicVars%fcoru => basicm_g(1)%fcoru
    
    oneBasicVars%fcorv => basic_g(1)%fcorv
    oneAveBasicVars%fcorv => basicm_g(1)%fcorv
    
    oneBasicVars%cputime => basic_g(1)%cputime
    oneAveBasicVars%cputime => basicm_g(1)%cputime
    
  end subroutine DeepCopyToBasicVars




  subroutine DeepCopyFromBasicVars(oneBasicVars, oneAveBasicVars)
    type(BasicVars), pointer, intent(in) :: oneBasicVars
    type(BasicVars), pointer, intent(in) :: oneAveBasicVars

    character(len=*), parameter :: h="**(DeepCopyFromBasicVars)**"

    if (.not. associated(oneBasicVars)) then
       call fatal_error(h//" oneBasicVars not associated")
    else if (.not. associated(oneAveBasicVars)) then
       call fatal_error(h//" oneAveBasicVars not associated")
    end if

    basic_g(1)%up => oneBasicVars%up
    basicm_g(1)%up => oneAveBasicVars%up

    basic_g(1)%uc => oneBasicVars%uc
    basicm_g(1)%uc => oneAveBasicVars%uc

    basic_g(1)%vp => oneBasicVars%vp
    basicm_g(1)%vp => oneAveBasicVars%vp

    basic_g(1)%vc => oneBasicVars%vc
    basicm_g(1)%vc => oneAveBasicVars%vc

    basic_g(1)%wp => oneBasicVars%wp
    basicm_g(1)%wp => oneAveBasicVars%wp

    basic_g(1)%wc => oneBasicVars%wc
    basicm_g(1)%wc => oneAveBasicVars%wc

    basic_g(1)%pp => oneBasicVars%pp
    basicm_g(1)%pp => oneAveBasicVars%pp

    basic_g(1)%pc => oneBasicVars%pc
    basicm_g(1)%pc => oneAveBasicVars%pc

    basic_g(1)%rv => oneBasicVars%rv
    basicm_g(1)%rv => oneAveBasicVars%rv

    basic_g(1)%theta => oneBasicVars%theta
    basicm_g(1)%theta => oneAveBasicVars%theta

    basic_g(1)%thp => oneBasicVars%thp
    basicm_g(1)%thp => oneAveBasicVars%thp

    basic_g(1)%thc => oneBasicVars%thc
    basicm_g(1)%thc => oneAveBasicVars%thc

    basic_g(1)%rtp => oneBasicVars%rtp
    basicm_g(1)%rtp => oneAveBasicVars%rtp

    basic_g(1)%pi0 => oneBasicVars%pi0
    basicm_g(1)%pi0 => oneAveBasicVars%pi0

    basic_g(1)%th0 => oneBasicVars%th0
    basicm_g(1)%th0 => oneAveBasicVars%th0

    basic_g(1)%dn0 => oneBasicVars%dn0
    basicm_g(1)%dn0 => oneAveBasicVars%dn0

    basic_g(1)%dn0u => oneBasicVars%dn0u
    basicm_g(1)%dn0u => oneAveBasicVars%dn0u

    basic_g(1)%dn0v => oneBasicVars%dn0v
    basicm_g(1)%dn0v => oneAveBasicVars%dn0v

    basic_g(1)%fcoru => oneBasicVars%fcoru
    basicm_g(1)%fcoru => oneAveBasicVars%fcoru

    basic_g(1)%fcorv => oneBasicVars%fcorv
    basicm_g(1)%fcorv => oneAveBasicVars%fcorv

    basic_g(1)%cputime => oneBasicVars%cputime
    basicm_g(1)%cputime => oneAveBasicVars%cputime

  end subroutine DeepCopyFromBasicVars
end module ModBasicVars
