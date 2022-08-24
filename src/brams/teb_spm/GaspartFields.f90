module ModGaspartFields

  use iso_fortran_env, only: &
       int64

  use ModVarTable, only: &
       VarTable, &
       InsertVarTable_2D, &
       InsertVarTable
  
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
     real, contiguous, pointer :: pco(:,:,:) => null()
     real, contiguous, pointer :: pno(:,:,:) => null()
     real, contiguous, pointer :: pno2(:,:,:) => null()
     real, contiguous, pointer :: ppm25(:,:,:) => null()
     real, contiguous, pointer :: pvoc(:,:,:) => null()
     real, contiguous, pointer :: pso2(:,:,:) => null()
     real, contiguous, pointer :: proo(:,:,:) => null()
     real, contiguous, pointer :: pso4(:,:,:) => null()
     real, contiguous, pointer :: paer(:,:,:) => null()
     real, contiguous, pointer :: po3(:,:,:) => null()
     real, contiguous, pointer :: prhco(:,:,:) => null()
     real, contiguous, pointer :: pho2(:,:,:) => null()
     real, contiguous, pointer :: po3p(:,:,:) => null()
     real, contiguous, pointer :: po1d(:,:,:) => null()
     real, contiguous, pointer :: pho(:,:,:) => null()
     real, contiguous, pointer :: gasr(:,:,:) => null()
     real, contiguous, pointer :: peoxid(:,:,:) => null()

     ! variables to be dimensioned by (nzp,nxp)
     real, contiguous, pointer :: fusog(:,:) => null()

     real, contiguous, pointer :: pcot(:) => null()
     real, contiguous, pointer :: pnot(:) => null()
     real, contiguous, pointer :: pno2t(:) => null()
     real, contiguous, pointer :: ppm25t(:) => null()
     real, contiguous, pointer :: pvoct(:) => null()
     real, contiguous, pointer :: pso2t(:) => null()
     real, contiguous, pointer :: pso4t(:) => null()
     real, contiguous, pointer :: paert(:)  => null()
     real, contiguous, pointer :: po3t(:) => null()
     real, contiguous, pointer :: prhcot(:) => null()
     real, contiguous, pointer :: pho2t(:) => null()
     real, contiguous, pointer :: po3pt(:) => null()
     real, contiguous, pointer :: po1dt(:) => null()
     real, contiguous, pointer :: phot(:) => null()
     real, contiguous, pointer :: proot(:) => null()

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
       


  


  subroutine InsertGaspartFieldsAtVarTable(oneVarTable, oneVarTableSize, &
       oneGaspartFields, oneAveGaspartFields)

    ! Build VarTable entry with GaspartFields components

    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(GaspartFields), pointer, intent(in) :: oneGaspartFields
    type(GaspartFields), pointer, intent(in) :: oneAveGaspartFields

    integer, parameter :: notImean=0
    character(len=*), parameter :: h="**(InsertGaspartFieldsAtVarTable)**" 

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(oneGaspartFields)) then
       call fatal_error(h//" oneGaspartFields not associated")
    end if

    if (associated(oneGaspartFields%fusog)) then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%fusog, &
                  'FUSOG:2:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%fusog, &
                  notImean)
    end if

    if (associated(oneGaspartFields%pco))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%pco, &
                  'PCO:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%pco, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%pno))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%pno, &
                  'PNO:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%pno, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%pno2))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%pno2, &
                  'PNO2:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%pno2, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%ppm25))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%ppm25, &
                  'PPM25:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%ppm25, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%pvoc))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%pvoc, &
                  'PVOC:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%pvoc, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%pso2))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%pso2, &
                  'PSO2:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%pso2, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%pso4))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%pso4, &
                  'PSO4:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%pso4, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%paer))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%paer, &
                  'PAER:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%paer, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%peoxid))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%peoxid, &
                  'PEOXID:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%peoxid, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%gasr))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%gasr, &
                  'GASR:3:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%gasr, &
                  notImean)
    end if

    if (associated(oneGaspartFields%po3))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%po3, &
                  'PO3:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%po3, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%prhco))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%prhco, &
                  'PRHCO:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%prhco, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%pho2))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%pho2, &
                  'PHO2:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%pho2, &
                  notImean)
    end if
       
    if (associated(oneGaspartFields%po3p))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%po3p, &
                  'PO3P:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%po3p, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%po1d))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%po1d, &
                  'PO1D:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%po1d, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%pho))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%pho, &
                  'PHO:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%pho, &
                  notImean)
    end if
    
    if (associated(oneGaspartFields%proo))  then
             call InsertVarTable (oneVarTable, oneVarTableSize, &
                  oneGaspartFields%proo, &
                  'PROO:3:hist:anal:mpti:mpt3:mpt1', &
                  oneAveGaspartFields%proo, &
                  notImean)
    end if
    
  end subroutine InsertGaspartFieldsAtVarTable


!!$  subroutine InsertGaspartFieldsAtVarTable(oneVarTable, oneVarTableSize, &
!!$       oneGaspartFields)
!!$
!!$    ! Build VarTable entry with GaspartFields components
!!$
!!$    type(VarTable), pointer, intent(in) :: oneVarTable(:)
!!$    integer, intent(inout) :: oneVarTableSize
!!$    type(GaspartFields), pointer, intent(in) :: oneGaspartFields
!!$
!!$    type(GaspartFields), contiguous, pointer :: noField_2D(:,:)
!!$    integer, parameter :: notImean=0
!!$    character(len=*), parameter :: h="**(InsertGaspartFieldsAtVarTable)**" 
!!$
!!$    if (.not. associated(oneVarTable)) then
!!$       call fatal_error(h//" oneVarTable not associated")
!!$    else if (.not. associated(oneGaspartFields)) then
!!$       call fatal_error(h//" oneGaspartFields not associated")
!!$    end if
!!$
!!$    nullify(noField_2D)
!!$    
!!$    if (associated(oneGaspartFields%fusog)) then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%fusog, &
!!$                  'FUSOG:2:hist:anal:mpti:mpt3:mpt1', &
!!$                  noField_2D, &
!!$                  notImean)
!!$    end if
!!$
!!$    if (associated(oneGaspartFields%pco))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%pco, &
!!$                  'PCO:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%pno))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%pno, &
!!$                  'PNO:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%pno2))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%pno2, &
!!$                  'PNO2:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%ppm25))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%ppm25, &
!!$                  'PPM25:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%pvoc))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%pvoc, &
!!$                  'PVOC:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%pso2))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%pso2, &
!!$                  'PSO2:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%pso4))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%pso4, &
!!$                  'PSO4:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%paer))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%paer, &
!!$                  'PAER:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%peoxid))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%peoxid, &
!!$                  'PEOXID:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%gasr))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%gasr, &
!!$                  'GASR:3:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$
!!$    if (associated(oneGaspartFields%po3))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%po3, &
!!$                  'PO3:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%prhco))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%prhco, &
!!$                  'PRHCO:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%pho2))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%pho2, &
!!$                  'PHO2:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$       
!!$    if (associated(oneGaspartFields%po3p))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%po3p, &
!!$                  'PO3P:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%po1d))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%po1d, &
!!$                  'PO1D:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%pho))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%pho, &
!!$                  'PHO:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$    if (associated(oneGaspartFields%proo))  then
!!$             call InsertVarTable (oneVarTable, oneVarTableSize, &
!!$                  oneGaspartFields%proo, &
!!$                  'PROO:3:hist:anal:mpti:mpt3:mpt1', &
!!$                  null(), &
!!$                  notImean)
!!$    end if
!!$    
!!$  end subroutine InsertGaspartFieldsAtVarTable  
end module ModGaspartFields
