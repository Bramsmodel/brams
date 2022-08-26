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
       int64

  use ModVarTable, only: &
       VarTable, &
       InsertVarTable

  use mem_stilt, only: &
       iexev

  implicit none

  private
  public :: BasicFields
  public :: CreateBasicFields
  public :: DestroyBasicFields
  public :: DumpBasicFields
  public :: InsertBasicFieldsAtVarTable 
 
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
  character(len=64) :: fromProcedure=""
  
contains




  function CreateBasicFields(oneNodeDimensions, oneNamelistFile) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer :: res

    integer :: ierr
    integer :: mzp
    integer :: mxp
    integer :: myp
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateBasicFields)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneNodeDimensions)) then
       call fatal_error(h//" oneNodeDimensions not associated")
    else if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if
       
    mzp=oneNodeDimensions%mzp
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp

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

    if (associated(oneBasicFields%up)) then
       deallocate(oneBasicFields%up, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate up fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%uc)) then
       deallocate(oneBasicFields%uc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate uc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%vp)) then
       deallocate(oneBasicFields%vp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%vc)) then
       deallocate(oneBasicFields%vc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate vc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%wp)) then
       deallocate(oneBasicFields%wp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate wp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%wc)) then
       deallocate(oneBasicFields%wc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate wc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%pp)) then
       deallocate(oneBasicFields%pp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%pc)) then
       deallocate(oneBasicFields%pc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%rv)) then
       deallocate(oneBasicFields%rv, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rv fails with stat="//&
               trim(adjustl(str(1))))

       end if
    end if

    if (associated(oneBasicFields%theta)) then
       deallocate(oneBasicFields%theta, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate theta fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%thp)) then
       deallocate(oneBasicFields%thp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate thp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%thc)) then
       deallocate(oneBasicFields%thc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate thc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%rtp)) then
       deallocate(oneBasicFields%rtp, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate rtp fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%pi0)) then
       deallocate(oneBasicFields%pi0, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pi0 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%th0)) then
       deallocate(oneBasicFields%th0, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate th0 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%dn0)) then
       deallocate(oneBasicFields%dn0, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate dn0 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%dn0u)) then
       deallocate(oneBasicFields%dn0u, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate dn0u fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%dn0v)) then
       deallocate(oneBasicFields%dn0v, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate dn0v fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%fcoru)) then
       deallocate(oneBasicFields%fcoru, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate fcoru fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%fcorv)) then
       deallocate(oneBasicFields%fcorv, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate fcorv fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields%cputime)) then
       deallocate(oneBasicFields%cputime, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate cputime fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    if (associated(oneBasicFields)) then
       deallocate(oneBasicFields, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate oneBasicFields fails with stat="//&
               trim(adjustl(str(1))))
       end if
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


  


  subroutine InsertBasicFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneBasicFields, oneAveBasicFields, imean)

    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(BasicFields), pointer, intent(in) :: oneAveBasicFields
    integer, intent(in) :: imean

    character(len=*), parameter :: h="**(InsertBasicFieldsAtVarTable)**" 

    if (.not. associated(oneBasicFields)) then
       call fatal_error(h//" oneBasicFields not associated")
    else if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    end if

    ! Fill pointers to arrays into variable tables

    if (associated(oneBasicFields%up)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%up, &
            'UP :3:hist:anal:mpti:mpt3:mpt2', &
            oneAveBasicFields%up, imean)
    end if

    if (associated(oneBasicFields%vp))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%vp, &
            'VP :3:hist:anal:mpti:mpt3:mpt2', &
            oneAveBasicFields%vp, imean)
    end if


    if (associated(oneBasicFields%wp))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%wp, &
            'WP :3:hist:anal:mpti:mpt3:mpt2', &
            oneAveBasicFields%wp, imean)
    end if

    if (associated(oneBasicFields%pp))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%pp, &
            'PP :3:hist:anal:mpti:mpt3:mpt2', &
            oneAveBasicFields%pp, imean)
    end if

    if (associated(oneBasicFields%uc))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%uc, &
            'UC :3:hist:mpti:mpt3:mpt2', &
            oneAveBasicFields%uc, imean)
    end if

    if (associated(oneBasicFields%vc))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%vc, &
            'VC :3:hist:mpti:mpt3:mpt2', &
            oneAveBasicFields%vc, imean)
    end if

    if (associated(oneBasicFields%wc))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%wc, &
            'WC :3:hist:mpti:mpt3:mpt2', &
            oneAveBasicFields%wc, imean)
    end if

    if (associated(oneBasicFields%pc))  then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%pc, &
            'PC :3:hist:mpti:mpt3:mpt2', &
            oneAveBasicFields%pc, imean)
    end if


    if (associated(oneBasicFields%thp)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%thp, &
            'THP :3:hist:mpti:mpt3:mpt1', &
            oneAveBasicFields%thp, imean)
    end if

    if (associated(oneBasicFields%rtp)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%rtp, &
            'RTP :3:hist:mpti:mpt3:mpt1', &
            oneAveBasicFields%rtp, imean)
    end if


    if(iexev == 2) then
       if (associated(oneBasicFields%theta)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneBasicFields%theta, &
               'THETA :3:hist:anal:mpti:mpt3:mpt1', &
               oneAveBasicFields%theta, imean)
       end if
    else
       if (associated(oneBasicFields%theta)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneBasicFields%theta, &
               'THETA :3:hist:anal:mpti:mpt3', &
               oneAveBasicFields%theta, imean)
       end if
    endif

    if (associated(oneBasicFields%thc)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%thc, &
            'THC :3:hist:mpti:mpt3:mpt1', &
            oneAveBasicFields%thc, imean)
    end if

    if(iexev == 2) then
       if (associated(oneBasicFields%rv)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneBasicFields%rv, &
               'RV :3:hist:anal:mpti:mpt3:mpt1', &
               oneAveBasicFields%rv, imean)
       end if
    else
       if (associated(oneBasicFields%rv)) then
          call InsertVarTable (oneVarTable, oneVarTableSize, &
               oneBasicFields%rv, &
               'RV :3:hist:anal:mpti:mpt3', &
               oneAveBasicFields%rv, imean)
       end if
    endif

    if (associated(oneBasicFields%pi0)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%pi0, &
            'PI0 :3:mpti', &
            oneAveBasicFields%pi0, imean)
    end if

    if (associated(oneBasicFields%th0)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%th0, &
            'TH0 :3:mpti', &
            oneAveBasicFields%th0, imean)
    end if

    if (associated(oneBasicFields%dn0)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%dn0, &
            'DN0 :3:mpti', &
            oneAveBasicFields%dn0, imean)
    end if

    if (associated(oneBasicFields%dn0u)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%dn0u, &
            'DN0U :3:mpti', &
            oneAveBasicFields%dn0u, imean)
    end if

    if (associated(oneBasicFields%dn0v)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%dn0v, &
            'DN0V :3:mpti', &
            oneAveBasicFields%dn0v, imean)
    end if


    if (associated(oneBasicFields%fcoru)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%fcoru, &
            'FCORU :2:mpti', &
            oneAveBasicFields%fcoru, imean)
    else
       call fatal_error(h//" oneBasicFields%fcoru not associated")
    end if

    if (associated(oneBasicFields%fcorv)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%fcorv, &
            'FCORV :2:mpti', &
            oneAveBasicFields%fcorv, imean)
    end if

    if (associated(oneBasicFields%cputime)) then
       call InsertVarTable (oneVarTable, oneVarTableSize, &
            oneBasicFields%cputime, &
            'CPUTIME :2:anal:mpti:mpt3', &
            oneAveBasicFields%cputime, imean)
    end if
  end subroutine InsertBasicFieldsAtVarTable
end module ModBasicFields
