module ModChem1Vars

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModNamelistFile, only: &
       NamelistFile
  
  use ModControlVars, only: &
       ControlVars

  use ModGridDims, only: &
       GridDims
  
  use chem1_list, only : &
       maxnspecies

  use ModChem1Constants, only: &
       antro, &
       bburn, &
       bioge, &
       geoge, &
       max_ntimes_src, &
       nsrc
  
  implicit none

  private
  public :: Chem1Vars
  public :: CreateChem1Vars
  public :: DestroyChem1Vars
  public :: DumpChem1Vars
  
  type Chem1Vars
     
     !- dimension of sources arrays (=1 for 2dim, =m1 for 3dim)
     integer, pointer, contiguous :: chem1_src_z_dim_g(:) => null()

     !- use of the prescribed diurnal cycle or linear interpolation 
     !- for the instantaneous emission rate
     !- for biomass burning, the diur_cycle must be always 1
     ! (the 2nd element of diur_cycle array)
     !- 1=on, 0=off (=> will use linterp.)
     ! diur_cycle(1)== antro; diur_cycle(2) == bburn
     ! diur_cycle(3)== bioge; diur_cycle(4) == geoge
     integer, pointer, contiguous :: diur_cycle(:)  => null()

     !- actual number used
     integer, pointer, contiguous :: ntimes_src(:) => null()

     integer :: recycle_tracers
     integer :: nspecies_transported
     integer :: nspecies_chem_transported
     integer :: nspecies_chem_no_transported
     ! determine if 4DDA will or  not be used
     integer :: chem_assim
     ! define if the chemistry will run and in this case the type of solver
     integer :: chemistry
     ! to control when chemistry will be called
     integer :: isplit
     ! determine the type of splliting method
     character(LEN=20) :: split_method
     ! timestep of chemistry integration
     real chem_timestep
     ! current n_dyn_chem_n
     integer :: n_dyn_chem
!!$     integer, dimension(maxnspecies) :: TRANSP_CHEM_INDEX,NO_TRANSP_CHEM_INDEX
     integer, pointer, contiguous :: transp_chem_index(:) => null()
     integer, pointer, contiguous :: no_transp_chem_index(:) => null()
!!$     integer :: N_DYN_CHEM_N(maxgrds)  ! number of dynamic timesteps per chem. timestep
     ! number of dynamic timesteps per chem. timestep
     integer :: n_dyn_chem_n     
  end type Chem1Vars

contains


  function CreateChem1Vars(oneParallelEnvironment, oneNamelistFile, &
       oneControlVars, oneGridDims, gridId) result(oneChem1Vars)
    type(ParallelEnvironment), pointer, intent(in) :: oneParallelEnvironment
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(ControlVars), pointer, intent(in) :: oneControlVars
    type(GridDims), pointer, intent(in) :: oneGridDims
    type(Chem1Vars), pointer :: oneChem1Vars
    integer, intent(in) :: gridId

    ! include code from subroutine define_n_dyn_chem

    integer :: ierr
    integer :: isrc
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateChem1Vars)**"

    allocate(oneChem1Vars, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" allocate oneChem1Vars fails with stat="//&
            trim(adjustl(str(1))))
    end if

    ! set scalar values
    oneChem1Vars%recycle_tracers = oneNamelistFile%recycle_tracers
    oneChem1Vars%nspecies_transported=0
    oneChem1Vars%nspecies_chem_transported=0
    oneChem1Vars%nspecies_chem_no_transported=0
    oneChem1Vars%chem_assim = oneNamelistFile%chem_assim
    oneChem1Vars%chemistry = oneNamelistFile%chemistry

    !- define split control to setup when chem will be called   
    if(trim(adjustl(oneNamelistFile%split_method))=='SYMMETRIC') then
       oneChem1Vars%isplit = 2
    else
       oneChem1Vars%isplit = 1
    end if

    oneChem1Vars%split_method = oneNamelistFile%split_method
    oneChem1Vars%chem_timestep = oneNamelistFile%chem_timestep

    !- n_dyn_chem for the coarser grid
    oneChem1Vars%n_dyn_chem_n = max(1,&
         nint(oneNamelistFile%chem_timestep/oneNamelistFile%dtlong)) 
    !- correct for a nested grid
    if (gridId > 1) then
       oneChem1Vars%n_dyn_chem_n = min(oneControlVars%nndtrat, &
            oneChem1Vars%n_dyn_chem_n)
    end if

    allocate(oneChem1Vars%chem1_src_z_dim_g(nsrc), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       write(str(2),"(i8)") nsrc
       call fatal_error(h//" allocate chem1_src_z_dim_g("//&
            trim(adjustl(str(2)))//") fails with stat="//&
            trim(adjustl(str(1))))
    end if
    !- determination of the dimension of Z-dir of source field array
    oneChem1Vars%chem1_src_z_dim_g(antro) = oneGridDims%nnzp    ! 2d 
    oneChem1Vars%chem1_src_z_dim_g(bburn) = oneGridDims%nnzp    ! 3d
    oneChem1Vars%chem1_src_z_dim_g(bioge) = oneGridDims%nnzp    ! 2d
    oneChem1Vars%chem1_src_z_dim_g(geoge) = oneGridDims%nnzp    ! 3d for volcanoes

    allocate(oneChem1Vars%diur_cycle(nsrc), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       write(str(2),"(i8)") nsrc
       call fatal_error(h//" allocate diur_cycle("//&
            trim(adjustl(str(2)))//") fails with stat="//&
            trim(adjustl(str(1))))
    end if
    oneChem1Vars%diur_cycle(:) = oneNamelistFile%diur_cycle(:)

    allocate(oneChem1Vars%ntimes_src(nsrc), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       write(str(2),"(i8)") nsrc
       call fatal_error(h//" allocate ntimes_src("//&
            trim(adjustl(str(2)))//") fails with stat="//&
            trim(adjustl(str(1))))
    end if
    !-determine the time dimension for src arrays
    do isrc=1,nsrc
       ! if   : ntimes_src(isrc) = 1, will use the prescribed diurnal cycle
       ! else : the array will be allocated with 2 dimensions to
       !        read 2 files and make linear interpolation
       oneChem1Vars%ntimes_src(isrc)= 2 - oneChem1Vars%diur_cycle(isrc)

       if (oneChem1Vars%ntimes_src(isrc) .gt. max_ntimes_src  .or. &
            oneChem1Vars%ntimes_src(isrc) .lt. 1  ) then
          call fatal_error(h//' ntimes_src > max_ntimes_src or < 1')
       end if

       if (oneChem1Vars%ntimes_src(bburn) > 1) then
          call fatal_error(h//" bburn must run with diurnal cycle ON")
       end if

       if(oneChem1Vars%ntimes_src(geoge) > 1) then
          call fatal_error(h//" geoge must run with diurnal cycle ON")
       end if
    end do

    allocate(oneChem1Vars%transp_chem_index(maxnspecies), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       write(str(2),"(i8)") maxnspecies
       call fatal_error(h//" allocate transp_chem_index("//&
            trim(adjustl(str(2)))//") fails with stat="//&
            trim(adjustl(str(1))))
    end if
    oneChem1Vars%transp_chem_index(:)=0
    
    allocate(oneChem1Vars%no_transp_chem_index(maxnspecies), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       write(str(2),"(i8)") maxnspecies
       call fatal_error(h//" allocate no_transp_chem_index("//&
            trim(adjustl(str(2)))//") fails with stat="//&
            trim(adjustl(str(1))))
    end if
    oneChem1Vars%no_transp_chem_index(:)=0

    !- get the number of dynamics cycles inside each chemistry cycle:
    !- observe that 'dtlong' (timestep of coarser grid) is used, instead
    !- of dtlongn(ngrid) , the grid dependent long timestep.

    
    if((oneParallelEnvironment%mynum == 0 .or. oneParallelEnvironment%mynum ==1) .and. &
         oneChem1Vars%chemistry > 0) then
       print*,'----------------------------------------------------------------'
       print*,' -- > chemistry splitting: N_DYN_CHEM=',oneChem1Vars%n_dyn_chem_n
       print*,'----------------------------------------------------------------'
    endif

  end function CreateChem1Vars



  subroutine DestroyChem1Vars(oneChem1Vars)
    type(Chem1Vars), pointer, intent(inout) :: oneChem1Vars

    integer :: ierr
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DestroyChem1Vars)**"
    
    if (.not. associated(oneChem1Vars)) then
       return
    end if
    
    deallocate(oneChem1Vars%chem1_src_z_dim_g, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate chem1_src_z_dim_g fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneChem1Vars%diur_cycle, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate diur_cycle fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneChem1Vars%ntimes_src, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate ntimes_src fails with stat="//&
            trim(adjustl(str(1))))
    end if

    deallocate(oneChem1Vars%transp_chem_index, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate transp_chem_index fails with stat="//&
            trim(adjustl(str(1))))
    end if
    
    deallocate(oneChem1Vars%no_transp_chem_index, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate no_transp_chem_index fails with stat="//&
            trim(adjustl(str(1))))
    end if
    
    deallocate(oneChem1Vars, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate oneChem1Vars fails with stat="//&
            trim(adjustl(str(1))))
    end if
    nullify(oneChem1Vars)
  end subroutine DestroyChem1Vars



  subroutine DumpChem1Vars(oneChem1Vars, name)
    type(Chem1Vars), pointer, intent(in) :: oneChem1Vars
    character(len=*), intent(in) :: name

    integer :: isrc
    character(len=8) :: str(10)
    character(len=16) :: strr
    character(len=*), parameter :: h="**(DumpChem1Vars)**"

    if (.not. associated(oneChem1Vars)) then
       call MsgDump(h//" at "//trim(adjustl(name))//" oneChem1Vars not associated")
       return
    else
       call MsgDump(h//" will dump data from "//trim(adjustl(name)))
    end if

    if (associated(oneChem1Vars%chem1_src_z_dim_g)) then
       do isrc=1,size(oneChem1Vars%chem1_src_z_dim_g)
          write(str(1),"(i8)") isrc
          write(str(2),"(i8)") oneChem1Vars%chem1_src_z_dim_g(isrc)
          call MsgDump(h//" oneChem1Vars%chem1_src_z_dim_g("//&
               trim(adjustl(str(1)))//") = "//trim(adjustl(str(2))))
       end do
    else
       call MsgDump(h//" oneChem1Vars%chem1_src_z_dim_g not associated")
    end if
    
    if (associated(oneChem1Vars%diur_cycle)) then
       do isrc=1,size(oneChem1Vars%diur_cycle)
          write(str(1),"(i8)") isrc
          write(str(2),"(i8)") oneChem1Vars%diur_cycle(isrc)
          call MsgDump(h//" oneChem1Vars%diur_cycle("//&
               trim(adjustl(str(1)))//") = "//trim(adjustl(str(2))))
       end do
    else
       call MsgDump(h//" oneChem1Vars%diur_cycle not associated")
    end if
    
    if (associated(oneChem1Vars%ntimes_src)) then
       do isrc=1,size(oneChem1Vars%ntimes_src)
          write(str(1),"(i8)") isrc
          write(str(2),"(i8)") oneChem1Vars%ntimes_src(isrc)
          call MsgDump(h//" oneChem1Vars%ntimes_src("//&
               trim(adjustl(str(1)))//") = "//trim(adjustl(str(2))))
       end do
    else
       call MsgDump(h//" oneChem1Vars%ntimes_src not associated")
    end if
    
    if (associated(oneChem1Vars%transp_chem_index)) then
       do isrc=1,size(oneChem1Vars%transp_chem_index)
          write(str(1),"(i8)") isrc
          write(str(2),"(i8)") oneChem1Vars%transp_chem_index(isrc)
          call MsgDump(h//" oneChem1Vars%transp_chem_index("//&
               trim(adjustl(str(1)))//") = "//trim(adjustl(str(2))))
       end do
    else
       call MsgDump(h//" oneChem1Vars%transp_chem_index not associated")
    end if
    
    if (associated(oneChem1Vars%no_transp_chem_index)) then
       do isrc=1,size(oneChem1Vars%no_transp_chem_index)
          write(str(1),"(i8)") isrc
          write(str(2),"(i8)") oneChem1Vars%no_transp_chem_index(isrc)
          call MsgDump(h//" oneChem1Vars%no_transp_chem_index("//&
               trim(adjustl(str(1)))//") = "//trim(adjustl(str(2))))
       end do
    else
       call MsgDump(h//" oneChem1Vars%no_transp_chem_index not associated")
    end if

    write(str(1),"(i8)") oneChem1Vars%recycle_tracers
    call MsgDump(h//" recycle_tracers="//trim(adjustl(str(1))))
    
    write(str(1),"(i8)") oneChem1Vars%nspecies_transported
    call MsgDump(h//" nspecies_transported="//trim(adjustl(str(1))))
    
    write(str(1),"(i8)") oneChem1Vars%nspecies_chem_transported
    call MsgDump(h//" nspecies_chem_transported="//trim(adjustl(str(1))))
    
    write(str(1),"(i8)") oneChem1Vars%nspecies_chem_no_transported
    call MsgDump(h//" nspecies_chem_no_transported="//trim(adjustl(str(1))))
    
    write(str(1),"(i8)") oneChem1Vars%chem_assim
    call MsgDump(h//" chem_assim="//trim(adjustl(str(1))))
    
    write(str(1),"(i8)") oneChem1Vars%chemistry
    call MsgDump(h//" chemistry="//trim(adjustl(str(1))))
    
    write(str(1),"(i8)") oneChem1Vars%isplit
    call MsgDump(h//" isplit="//trim(adjustl(str(1))))
    
    write(str(1),"(i8)") oneChem1Vars%n_dyn_chem
    call MsgDump(h//" n_dyn_chem="//trim(adjustl(str(1))))
    
    write(str(1),"(i8)") oneChem1Vars%n_dyn_chem_n
    call MsgDump(h//" n_dyn_chem_n     ="//trim(adjustl(str(1))))
    
    call MsgDump(h//" split_method="//trim(adjustl(oneChem1Vars%split_method)))

    write(strr,"(e15.7)") oneChem1Vars%chem_timestep
    call MsgDump(h//" chem_timestep="//trim(adjustl(strr)))
  end subroutine DumpChem1Vars
end module ModChem1Vars
