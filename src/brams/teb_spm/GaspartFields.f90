module ModGaspartFields

  use iso_fortran_env, only: &
       int64
  
  use ModVarTables, only: &
       InsertVtab

  use ModParallelEnvironment, only: &
       MsgDump

  use ModNamelistFile, only: &
       NamelistFile

  use ModNodeDimensions, only: &
       NodeDimensions
  
  implicit none

  private

  public :: GaspartFields
  public :: CreateGaspartFields
  public :: DestroyGaspartFields
  public :: DumpGaspartFields
  public :: InsertGaspartFieldsAtVarTable
  
  type GaspartFields

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, contiguous, pointer :: pco(:,:,:)
     real, contiguous, pointer :: pno(:,:,:)
     real, contiguous, pointer :: pno2(:,:,:)
     real, contiguous, pointer :: ppm25(:,:,:)
     real, contiguous, pointer :: pvoc(:,:,:)
     real, contiguous, pointer :: pso2(:,:,:)
     real, contiguous, pointer :: proo(:,:,:)
     real, contiguous, pointer :: pso4(:,:,:)
     real, contiguous, pointer :: paer(:,:,:)
     real, contiguous, pointer :: po3(:,:,:)
     real, contiguous, pointer :: prhco(:,:,:)
     real, contiguous, pointer :: pho2(:,:,:)
     real, contiguous, pointer :: po3p(:,:,:)
     real, contiguous, pointer :: po1d(:,:,:)
     real, contiguous, pointer :: pho(:,:,:)
     real, contiguous, pointer :: gasr(:,:,:)
     real, contiguous, pointer :: peoxid(:,:,:)

     ! variables to be dimensioned by (nzp,nxp)
     real, contiguous, pointer :: fusog(:,:)

     real, contiguous, pointer :: pcot(:)
     real, contiguous, pointer :: pnot(:)
     real, contiguous, pointer :: pno2t(:)
     real, contiguous, pointer :: ppm25t(:)
     real, contiguous, pointer :: pvoct(:)
     real, contiguous, pointer :: pso2t(:)
     real, contiguous, pointer :: pso4t(:)
     real, contiguous, pointer :: paert(:) 
     real, contiguous, pointer :: po3t(:)
     real, contiguous, pointer :: prhcot(:)
     real, contiguous, pointer :: pho2t(:)
     real, contiguous, pointer :: po3pt(:)
     real, contiguous, pointer :: po1dt(:)
     real, contiguous, pointer :: phot(:)
     real, contiguous, pointer :: proot(:)

  end type GaspartFields

contains



  function CreateGaspartFields(oneNodeDimensions, oneNamelistFile) result(res)
    type(NodeDimensions), pointer, intent(in) :: oneNodeDimensions
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(GaspartFields), pointer :: res

    integer :: ierr
    integer :: mzp
    integer :: mxp
    integer :: myp
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateGaspartFields)**"

    if (.not. associated(oneNodeDimensions)) then
       call fatal_error(h//" oneNodeDimensions not associated")
    else if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if

    mzp=oneNodeDimensions%mzp
    mxp=oneNodeDimensions%mxp
    myp=oneNodeDimensions%myp

    allocate(res, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate res fails with stat="//&
            trim(adjustl(str(1))))
    end if

    
    allocate(res%pco(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate pco fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%pco=0.0
    
    allocate(res%pno(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate pno fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%pno=0.0
    
    allocate(res%pno2(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate pno2 fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%pno2=0.0
    
    allocate(res%ppm25(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate ppm25 fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%ppm25=0.0
    
    allocate(res%pso2(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate pso2 fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%pso2=0.0
    
    allocate(res%pvoc(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate pvoc fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%pvoc=0.0
    
    allocate(res%gasr(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate gasr fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%gasr=0.0
    
    allocate(res%pso4(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate pso4 fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%pso4=0.0
    
    allocate(res%paer(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate paer fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%paer=0.0
    
    allocate(res%peoxid(mzp,mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate peoxid fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%peoxid=0.0
    
    allocate(res%fusog(mxp,myp), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate fusog fails with stat="//&
            trim(adjustl(str(1))))
    end if
    res%fusog=0.0
    
    if (oneNamelistFile%ichemi == 1) then

       allocate(res%po3(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate po3 fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%po3=0.0
       
       allocate(res%prhco(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate prhco fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%prhco=0.0
       
       allocate(res%pho2(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate pho2 fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%pho2=0.0
       
       allocate(res%po3p(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate po3p fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%po3p=0.0
       
       allocate(res%po1d(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate po1d fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%po1d=0.0
       
       allocate(res%pho(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate pho fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%pho=0.0
       
       allocate(res%proo(mzp,mxp,myp), stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate proo fails with stat="//&
               trim(adjustl(str(1))))
       end if
       res%proo=0.0
    endif

  end function CreateGaspartFields





  
  subroutine DestroyGaspartFields(oneGaspartFields)
    type(GaspartFields), pointer, intent(inout) :: oneGaspartFields

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyGaspartFields)**"

    if (.not. associated(oneGaspartFields)) then
       return
    end if

    if (associated(oneGaspartFields%pco)) then
       deallocate(oneGaspartFields%pco, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pco fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%pno)) then
       deallocate(oneGaspartFields%pno, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pno fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%pno2)) then
       deallocate(oneGaspartFields%pno2, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pno2 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%ppm25)) then
       deallocate(oneGaspartFields%ppm25, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate ppm25 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%pso2)) then
       deallocate(oneGaspartFields%pso2, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pso2 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%pvoc)) then
       deallocate(oneGaspartFields%pvoc, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pvoc fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%gasr)) then
       deallocate(oneGaspartFields%gasr, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate gasr fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%pso4)) then
       deallocate(oneGaspartFields%pso4, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pso4 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%paer)) then
       deallocate(oneGaspartFields%paer, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate paer fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%peoxid)) then
       deallocate(oneGaspartFields%peoxid, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate peoxid fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (associated(oneGaspartFields%fusog)) then
       deallocate(oneGaspartFields%fusog, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate fusog fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    

    if (associated(oneGaspartFields%po3)) then
       deallocate(oneGaspartFields%po3, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate po3 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
       
    if (associated(oneGaspartFields%prhco)) then
       deallocate(oneGaspartFields%prhco, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate prhco fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
       
    if (associated(oneGaspartFields%pho2)) then
       deallocate(oneGaspartFields%pho2, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pho2 fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
       
    if (associated(oneGaspartFields%po3p)) then
       deallocate(oneGaspartFields%po3p, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate po3p fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
       
    if (associated(oneGaspartFields%po1d)) then
       deallocate(oneGaspartFields%po1d, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate po1d fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
       
    if (associated(oneGaspartFields%pho)) then
       deallocate(oneGaspartFields%pho, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate pho fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
       
    if (associated(oneGaspartFields%proo)) then
       deallocate(oneGaspartFields%proo, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" deallocate proo fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
       
    deallocate(oneGaspartFields, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneGaspartFields fails with stat="//&
            trim(adjustl(str(1))))
    end if
    
    nullify(oneGaspartFields)
  end subroutine DestroyGaspartFields




  
  subroutine DumpGaspartFields(oneGaspartFields, name)
    type(GaspartFields), pointer, intent(inout) :: oneGaspartFields
    character(len=*), intent(in) :: name

    character(len=*), parameter :: h="**(DumpGaspartFields)**"

    if (associated(oneGaspartFields)) then
       call MsgDump(h//" "//name//" is associated")
    else
       call MsgDump(h//" "//name//" is not associated")
    end if
  end subroutine DumpGaspartFields
       


  


  subroutine InsertGaspartFieldsAtVarTable(oneGaspartFields, oneAveGaspartFields, &
       oneNamelistFile, gridId)
    type(GaspartFields), pointer, intent(in) :: oneGaspartFields
    type(GaspartFields), pointer, intent(in) :: oneAveGaspartFields
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    integer, intent(in) :: gridId

    integer :: imean
    integer(kind=int64) :: npts
    character(len=*), parameter :: h="**(InsertGaspartFieldsAtVarTable)**" 

    if (.not. associated(oneGaspartFields)) then
       call fatal_error(h//" oneGaspartFields not associated")
    else if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" oneNamelistFile not associated")
    end if

    ! Should average fields be stored at variable tables?

    if (oneNamelistFile%avgtim == 0) then
       imean=0 ! do not store
    else
       imean=1 ! store
    end if
    
    ! Fill pointers to arrays into variable tables


    if (associated(oneGaspartFields%fusog)) then
       npts=int(size(oneGaspartFields%fusog,1),int64) * &
            int(size(oneGaspartFields%fusog,2),int64)
       call InsertVTab (oneGaspartFields%FUSOG,oneGaspartFields%FUSOG &
            ,gridId, npts, imean,  &
            'FUSOG:2:hist:anal:mpti:mpt3:mpt1')
    end if


    if (associated(oneGaspartFields%pco))  then
       npts=int(size(oneGaspartFields%pco,1),int64) * &
            int(size(oneGaspartFields%pco,2),int64) * &
            int(size(oneGaspartFields%pco,3),int64)
       call InsertVTab (oneGaspartFields%PCO,oneGaspartFields%PCO&
            ,gridId, npts, imean,  &
            'PCO:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%pno))  then
       npts=int(size(oneGaspartFields%pno,1),int64) * &
            int(size(oneGaspartFields%pno,2),int64) * &
            int(size(oneGaspartFields%pno,3),int64)
       call InsertVTab (oneGaspartFields%PNO,oneGaspartFields%PNO&
            ,gridId, npts, imean,  &
            'PNO:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%pno2))  then
       npts=int(size(oneGaspartFields%pno2,1),int64) * &
            int(size(oneGaspartFields%pno2,2),int64) * &
            int(size(oneGaspartFields%pno2,3),int64)
       call InsertVTab (oneGaspartFields%PNO2,oneGaspartFields%PNO2&
            ,gridId, npts, imean,  &
            'PNO2:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%ppm25))  then
       npts=int(size(oneGaspartFields%ppm25,1),int64) * &
            int(size(oneGaspartFields%ppm25,2),int64) * &
            int(size(oneGaspartFields%ppm25,3),int64)
       call InsertVTab (oneGaspartFields%PPM25,oneGaspartFields%PPM25&
            ,gridId, npts, imean,  &
            'PPM25:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%pvoc))  then
       npts=int(size(oneGaspartFields%pvoc,1),int64) * &
            int(size(oneGaspartFields%pvoc,2),int64) * &
            int(size(oneGaspartFields%pvoc,3),int64)
       call InsertVTab (oneGaspartFields%PVOC,oneGaspartFields%PVOC&
            ,gridId, npts, imean,  &
            'PVOC:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%pso2))  then
       npts=int(size(oneGaspartFields%pso2,1),int64) * &
            int(size(oneGaspartFields%pso2,2),int64) * &
            int(size(oneGaspartFields%pso2,3),int64)
       call InsertVTab (oneGaspartFields%PSO2,oneGaspartFields%PSO2&
            ,gridId, npts, imean,  &
            'PSO2:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%pso4))  then
       npts=int(size(oneGaspartFields%pso4,1),int64) * &
            int(size(oneGaspartFields%pso4,2),int64) * &
            int(size(oneGaspartFields%pso4,3),int64)
       call InsertVTab (oneGaspartFields%PSO4,oneGaspartFields%PSO4&
            ,gridId, npts, imean,  &
            'PSO4:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%paer))  then
       npts=int(size(oneGaspartFields%paer,1),int64) * &
            int(size(oneGaspartFields%paer,2),int64) * &
            int(size(oneGaspartFields%paer,3),int64)
       call InsertVTab (oneGaspartFields%PAER,oneGaspartFields%PAER&
            ,gridId, npts, imean,  &
            'PAER:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%PEOXID))  then
       npts=int(size(oneGaspartFields%PEOXID,1),int64) * &
            int(size(oneGaspartFields%PEOXID,2),int64) * &
            int(size(oneGaspartFields%PEOXID,3),int64)
       call InsertVTab (oneGaspartFields%PEOXID,oneGaspartFields%PEOXID&
            ,gridId, npts, imean,  &
            'PEOXID:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%gasr))  then
       npts=int(size(oneGaspartFields%gasr,1),int64) * &
            int(size(oneGaspartFields%gasr,2),int64) * &
            int(size(oneGaspartFields%gasr,3),int64)
       call InsertVTab (oneGaspartFields%GASR,oneGaspartFields%GASR&
            ,gridId, npts, imean,  &
            'GASR:3:mpti:mpt3:mpt1')
    end if
    

    if (associated(oneGaspartFields%po3))  then
       npts=int(size(oneGaspartFields%po3,1),int64) * &
            int(size(oneGaspartFields%po3,2),int64) * &
            int(size(oneGaspartFields%po3,3),int64)
       call InsertVTab (oneGaspartFields%PO3,oneGaspartFields%PO3&
            ,gridId, npts, imean,  &
            'PO3:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%prhco))  then
       npts=int(size(oneGaspartFields%prhco,1),int64) * &
            int(size(oneGaspartFields%prhco,2),int64) * &
            int(size(oneGaspartFields%prhco,3),int64)
       call InsertVTab (oneGaspartFields%PRHCO,oneGaspartFields%PRHCO&
            ,gridId, npts, imean,  &
            'PRHCO:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%pho2))  then
       npts=int(size(oneGaspartFields%pho2,1),int64) * &
            int(size(oneGaspartFields%pho2,2),int64) * &
            int(size(oneGaspartFields%pho2,3),int64)
       call InsertVTab (oneGaspartFields%PHO2,oneGaspartFields%PHO2&
            ,gridId, npts, imean,  &
            'PHO2:3:hist:anal:mpti:mpt3:mpt1')
    end if
       
    if (associated(oneGaspartFields%po3p))  then
       npts=int(size(oneGaspartFields%po3p,1),int64) * &
            int(size(oneGaspartFields%po3p,2),int64) * &
            int(size(oneGaspartFields%po3p,3),int64)
       call InsertVTab (oneGaspartFields%PO3P,oneGaspartFields%PO3P&
            ,gridId, npts, imean,  &
            'PO3P:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%po1d))  then
       npts=int(size(oneGaspartFields%po1d,1),int64) * &
            int(size(oneGaspartFields%po1d,2),int64) * &
            int(size(oneGaspartFields%po1d,3),int64)
       call InsertVTab (oneGaspartFields%PO1D,oneGaspartFields%PO1D&
            ,gridId, npts, imean,  &
            'PO1D:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%pho))  then
       npts=int(size(oneGaspartFields%pho,1),int64) * &
            int(size(oneGaspartFields%pho,2),int64) * &
            int(size(oneGaspartFields%pho,3),int64)
       call InsertVTab (oneGaspartFields%PHO,oneGaspartFields%PHO&
            ,gridId, npts, imean,  &
            'PHO:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
    if (associated(oneGaspartFields%proo))  then
       npts=int(size(oneGaspartFields%proo,1),int64) * &
            int(size(oneGaspartFields%proo,2),int64) * &
            int(size(oneGaspartFields%proo,3),int64)
       call InsertVTab (oneGaspartFields%PROO,oneGaspartFields%PROO&
            ,gridId, npts, imean,  &
            'PROO:3:hist:anal:mpti:mpt3:mpt1')
    end if
    
  end subroutine InsertGaspartFieldsAtVarTable
end module ModGaspartFields
