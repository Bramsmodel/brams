module ReadBcst

  use ModMPassFull, only: &
       mk_2_buff, &
       mk_3_buff, &
       ex_2_buff, &
       ex_2p_buff, &
       ex_3_buff, &
       ex_4_buff, &
       ex_buff_carma
  
  use ModControlVars, only: &
       ControlVars

  use ModTurbFields, only: &
       TurbFields
  
  use ModBasicFields, only: &
       BasicFields

  use mem_grid, only: &
       ngrids, nzs, nzg, npatch, &
       time, iyear1, imonth1, idate1, itime1, &
       GlobalSizes, runtype


  use node_mod, only:      &
       mchnum, mynum, master_num, &
       nmachs, nodemxp, nodemyp, &
       nodeia, nodeiz, nodeja, nodejz, nodeibcon, &
       nodei0, nodej0

  use an_header, only: &
       IOFileDS,       &
       FieldWriteStoreInfo

  use mem_aerad, only: &
       nwave

  use ModVarTables, only: &
       num_var,      &
       vtab_r

  use ParLib, only: &
       parf_bcast, &
       parf_minloc, &
       parf_allreduce_max, &
       parf_GatherAllChunks

  implicit none


  private
  public :: ReadStoreOwnChunk
  public :: ReadStoreFullFieldAndOwnChunk
  public :: Broadcast
  public :: ProcWithMin
  public :: MaxCFLOverall
  public :: LocalSizesAndDisp
  public :: PreProcAndGather
  public :: RearrangeAndDump
  public :: gatherData
  public :: storeOwnChunk_3D


  interface ReadStoreOwnChunk
     module procedure ReadStoreOwnChunk_2D, ReadStoreOwnChunk_3D
  end interface ReadStoreOwnChunk

  interface Broadcast
     module procedure Broadcast_I, Broadcast_I1D, &
          Broadcast_R, Broadcast_R1D, Broadcast_R2D, &
          Broadcast_C, Broadcast_C1D
  end interface Broadcast

  interface gatherData 
     module procedure gatherData2d, gatherData3d, gatherData4d
  end interface gatherData

  integer, parameter :: idim_type_min=2
  integer, parameter :: idim_type_max=7
  logical, parameter :: dumpLocal=.false.
  include "constants.h"
contains



  subroutine ReadStoreOwnChunk_2D(grid, fUnit, toStore, fieldName, oneControlVars)

    integer, intent(in) :: grid
    integer, intent(in) :: fUnit
    real, pointer :: toStore(:,:)
    character(len=*), intent(in) :: fieldName
    type(ControlVars), pointer, intent(in) :: oneControlVars


    integer :: ldimx, ldimy, lin
    integer :: ia, iz, ja, jz
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(ReadStoreOwnChunk_2D)**"

    ! check allocated memory
    if (runtype(1:9)=='MAKEVFILE') then
       ldimx = oneControlVars%nnxp
       ldimy = oneControlVars%nnyp
    else
       ldimx = nodemxp(mynum,grid)
       ldimy = nodemyp(mynum,grid)
    endif
    if (.not. associated(toStore)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if
    lin = size(toStore,1)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    lin = size(toStore,2)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ia = nodei0(mynum,grid)+1
    iz = nodei0(mynum,grid)+nodemxp(mynum,grid)
    ja = nodej0(mynum,grid)+1
    jz = nodej0(mynum,grid)+nodemyp(mynum,grid)

    ! allocate scratch area for full domain

    call AllocReadStoreOwnChunk_2D(fUnit, toStore, fieldName, &
         oneControlVars%nnxp, oneControlVars%nnyp, ldimx, ldimy, ia, iz, ja, jz, &
         mchnum, master_num, runtype)

  end subroutine ReadStoreOwnChunk_2D




  subroutine AllocReadStoreOwnChunk_2D(fUnit, toStore, fieldName, &
       nnxp, nnyp, ldimx, ldimy, ia, iz, ja, jz, &
       mchnum, master_num, runtype)

    integer, intent(in) :: fUnit
    real, pointer :: toStore(:,:)
    character(len=*), intent(in) :: fieldName
    integer, intent(in) :: nnxp
    integer, intent(in) :: nnyp
    integer, intent(in) :: ldimx
    integer, intent(in) :: ldimy
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mchnum
    integer, intent(in) :: master_num
    character(len=*), intent(in) :: runtype


    real :: fullGrid(nnxp,nnyp)
    integer :: lin
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(AllocReadStoreOwnChunk_2D)**"
    integer :: i,j,n

    fullGrid=0.0

    !'master process opens file and reads first data into full domain scratch

    if (mchnum == master_num) then
       call vfirec(fUnit,fullGrid(1,1),nnxp*nnyp,'LIN')
    end if

    ! broadcast full domain scratch; 
    !local chunk is extracted and stored at desired variable
    if (runtype(1:9)/='MAKEVFILE') then
       call parf_bcast(fullGrid, int(nnxp,i8), int(nnyp,i8), &
            master_num)
    endif

    call mk_2_buff(fullGrid, toStore, &
         nnxp, nnyp, ldimx, ldimy, ia, iz, ja, jz)

    if (dumpLocal) then
       write(*,"(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)") &
            h//trim(fieldName)//"(   1:",ldimx, ",   1:",ldimy, &
            ")=fullGrid(",ia,":",iz,",",ja,":",jz,")"
    end if
  end subroutine AllocReadStoreOwnChunk_2D


  subroutine ReadStoreOwnChunk_3D(grid, fUnit, toStore, nz, fieldName, oneControlVars)

    integer, intent(in) :: grid
    integer, intent(in) :: fUnit
    integer, intent(in) :: nz
    real, pointer :: toStore(:,:,:)
    character(len=*), intent(in) :: fieldName
    type(ControlVars), pointer, intent(in) :: oneControlVars


    integer :: ldimx, ldimy, lin
    integer :: ia, iz, ja, jz, ka, kz
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(ReadStoreOwnChunk_3D)**"

    ! check allocated memory

    if (.not. associated(toStore)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if
    lin = size(toStore,1)
    if (nz /= lin) then
       write(c0,"(i8)") nz
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" z dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ldimx = nodemxp(mynum,grid)
    lin = size(toStore,2)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ldimy = nodemyp(mynum,grid)
    lin = size(toStore,3)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ia = nodei0(mynum,grid)+1
    iz = nodei0(mynum,grid)+nodemxp(mynum,grid)
    ja = nodej0(mynum,grid)+1
    jz = nodej0(mynum,grid)+nodemyp(mynum,grid)

    call AllocReadStoreOwnChunk_3D(fUnit, toStore, fieldName, &
         nz, oneControlVars%nnxp, oneControlVars%nnyp, ldimx, ldimy, ia, iz, ja, jz, &
         mchnum, master_num, runtype)

    if (dumpLocal) then
       write(*,"(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)") &
            h//trim(fieldName)//"(   1: ",nz, ",   1:",ldimx, ",   1:",ldimy, &
            ")=fullGrid(   1: ",nz,",",ia,":",iz,",",ja,":",jz,")"
    end if
  end subroutine ReadStoreOwnChunk_3D





  subroutine AllocReadStoreOwnChunk_3D(fUnit, toStore, fieldName, &
       nz, nnxp, nnyp, ldimx, ldimy, ia, iz, ja, jz, &
       mchnum, master_num, runtype)

    integer, intent(in) :: fUnit
    real, pointer :: toStore(:,:,:)
    character(len=*), intent(in) :: fieldName
    integer, intent(in) :: nz
    integer, intent(in) :: nnxp
    integer, intent(in) :: nnyp
    integer, intent(in) :: ldimx
    integer, intent(in) :: ldimy
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mchnum
    integer, intent(in) :: master_num
    character(len=*), intent(in) :: runtype

    real :: fullGrid(nz,nnxp,nnyp)
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(AllocReadStoreOwnChunk_3D)**"


    ! master process opens file and reads first data into full domain scratch

    if (mchnum == master_num) then
       call vfirec(fUnit,fullGrid(1,1,1),nz*nnxp*nnyp,'LIN')
    end if

    ! broadcast full domain scratch; 
    ! local chunk is extracted and stored at desired variable

    call parf_bcast(fullGrid, int(nz,i8), int(nnxp,i8), &
         int(nnyp,i8), master_num)

    call mk_3_buff(fullGrid, toStore, &
         nz, nnxp, nnyp, nz, ldimx, ldimy, ia, iz, ja, jz)

    if (dumpLocal) then
       write(*,"(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)") &
            h//trim(fieldName)//"(   1: ",nz, ",   1:",ldimx, ",   1:",ldimy, &
            ")=fullGrid(   1: ",nz,",",ia,":",iz,",",ja,":",jz,")"
    end if
  end subroutine AllocReadStoreOwnChunk_3D


  subroutine ReadStoreFullFieldAndOwnChunk(grid, fUnit, full, toStore, fieldName, oneControlVars)
    integer, intent(in) :: grid
    integer, intent(in) :: fUnit
    real, pointer :: full(:,:)
    real, pointer :: toStore(:,:)
    character(len=*), intent(in) :: fieldName
    type(ControlVars), pointer, intent(in) :: oneControlVars


    integer :: ldimx, ldimy, lin
    integer :: ia, iz, ja, jz
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(ReadStoreFullFieldAndOwnChunk)**"

    ! check full field allocated memory

    !print *,'LFR-DBG Inside read Bor: ', grid,fUnit,size(full,1),size(full,2),size(toStore,1),size(toStore,2),fieldName
    !print *,'LFR-DBG Inside read Nodemxp: ',mynum,grid,nodemxp(mynum,grid); call flush(6)

    if (.not. associated(full)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if
    ldimx = oneControlVars%nnxp
    lin = size(full,1)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//" full field of "//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ldimy = oneControlVars%nnyp
    lin = size(full,2)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//" full field of "//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if

    ! check local chunk allocated memory

    if (.not. associated(toStore)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if
    ldimx = nodemxp(mynum,grid)
    lin = size(toStore,1)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ldimy = nodemyp(mynum,grid)
    lin = size(toStore,2)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ia = nodei0(mynum,grid)+1
    iz = nodei0(mynum,grid)+nodemxp(mynum,grid)
    ja = nodej0(mynum,grid)+1
    jz = nodej0(mynum,grid)+nodemyp(mynum,grid)

    ! master process opens file and reads first data into full domain
    !print *, 'LFR-DBG->',fieldName,'Reading... vfirec',mchnum 
    if (mchnum == master_num) then
       call vfirec(fUnit,full(1,1),oneControlVars%nnxp*oneControlVars%nnyp,'LIN')
    end if
    !print *,'FieldName, Max e min lido: ',fieldName,maxval(full),minval(full),mchnum,master_num
    ! broadcast full domain; 
    ! local chunk is extracted and stored at desired variable

    call parf_bcast(full, int(oneControlVars%nnxp,i8), int(oneControlVars%nnyp,i8), &
         master_num)

    call mk_2_buff(full, toStore, &
         oneControlVars%nnxp, oneControlVars%nnyp, ldimx, ldimy, ia, iz, ja, jz)

    if (dumpLocal) then
       write(*,"(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)") &
            h//trim(fieldName)//"(   1:",ldimx, ",   1:",ldimy, &
            ")=full(",ia,":",iz,",",ja,":",jz,")"
    end if
  end subroutine ReadStoreFullFieldAndOwnChunk




  subroutine Broadcast_I(dataBcst, root, dataName)
    integer, intent(inout) :: dataBcst
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    integer :: IntArr(1)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_I)**"

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(c1, "(i8)") dataBcst
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " before broadcast "//trim(dataName)//" with value "//&
            trim(adjustl(c1))
    endif

    IntArr(1) = dataBcst
    call parf_bcast(IntArr, 1_i8, root)
    dataBcst = IntArr(1)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(c1, "(i8)") dataBcst
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)//" with value "//&
            trim(adjustl(c1))
    end if
  end subroutine Broadcast_I






  subroutine Broadcast_I1D(dataBcst, root, dataName)
    integer, intent(inout) :: dataBcst(:)
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_I1D)**"

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(c1, "(i8)") size(dataBcst)
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " before broadcast "//trim(dataName)//" with size "//&
            trim(adjustl(c1))
       call flush(6)
    endif

    call parf_bcast(dataBcst, int(size(dataBcst),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_I1D






  subroutine Broadcast_R(dataBcst, root, dataName)
    real, intent(inout) :: dataBcst
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    real :: RealArr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_R)**"

    RealArr = dataBcst
    call parf_bcast(RealArr, root)
    dataBcst = RealArr

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(c1, "(f8.4)") dataBcst
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)//" with value "//&
            trim(adjustl(c1))
    end if
  end subroutine Broadcast_R






  subroutine Broadcast_R1D(dataBcst, root, dataName)
    real, intent(inout)          :: dataBcst(:)
    integer, intent(in)          :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_R1D)**"

    ! broadcast dataBcst

    call parf_bcast(dataBcst, int(size(dataBcst),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_R1D






  subroutine Broadcast_R2D(dataBcst, root, dataName)
    real, intent(inout)          :: dataBcst(:,:)
    integer, intent(in)          :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_R2D)**"

    ! broadcast dataBcst

    call parf_bcast(dataBcst, int(size(dataBcst,1),i8), &
         int(size(dataBcst,2),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_R2D






  subroutine Broadcast_C(dataBcst, root, dataName)
    character(len=*), intent(inout) :: dataBcst
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(Broadcast_C)**"

    call parf_bcast(dataBcst, int(len(dataBcst),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_C






  subroutine Broadcast_C1D(dataBcst, root, dataName)
    character(len=*), intent(inout) :: dataBcst(:)
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(Broadcast_C1D)**"

    call parf_bcast(dataBcst, int(len(dataBcst),i8), int(size(dataBcst),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_C1D





  subroutine ProcWithMin(val, rank)
    real,    intent(inout) :: val
    integer, intent(out)   :: rank

    integer :: ierr
    real :: bufin(2)
    real :: bufout(2)

    bufin(:) = (/ val, real(mchnum) /)

    call parf_minloc(bufin, bufout)

    val = bufout(1)
    rank = nint(bufout(2))
  end subroutine ProcWithMin






  subroutine MaxCFLOverall(cflxy, cflz)
    real, intent(inout) :: cflxy(:)
    real, intent(inout) :: cflz(:)

    integer :: sizeCflxy
    integer :: sizeCflz
    integer :: sizeVec
    integer :: ierr
    real, allocatable :: vecIn(:)
    real, allocatable :: vecOut(:)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(MaxCFLOverall)**"

    sizeCflxy = size(cflxy)
    sizeCflz  = size(cflz)

    if (sizeCflxy /= sizeCflz) then
       write(c0,"(i8)") sizeCflxy
       write(c1,"(i8)") sizeCflz
       call fatal_error(h//" sizes of cflxy ("//trim(adjustl(c0))//&
            ") and cflz ("//trim(adjustl(c1))//" disagree")
    end if

    sizeVec=sizeCflxy+sizeCflz

    allocate(vecOut(sizeVec), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") sizeVec
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate vecOut ("//trim(adjustl(c0))//&
            ") failed with stat="//trim(adjustl(c1)))
    end if

    allocate(vecIn(sizeVec), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") sizeVec
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate vecIn ("//trim(adjustl(c0))//&
            ") failed with stat="//trim(adjustl(c1)))
    end if

    vecIn(1:sizeCflxy) = cflxy
    vecIn(sizeCflxy+1:sizeCflxy+sizeCflz) = cflz

    call parf_allreduce_max(vecIn, vecOut, int(sizeVec,i8))

    cflxy = vecOut(1:sizeCflxy)
    cflz  = vecOut(sizeCflxy+1:sizeCflxy+sizeCflz)

    deallocate(vecOut, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate vecOut failed with stat="//trim(adjustl(c1)))
    end if

    deallocate(vecIn, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate vecIn failed with stat="//trim(adjustl(c1)))
    end if
  end subroutine MaxCFLOverall



  ! LocalSizesAndDisp : sizes and displacements required by MPI_GATHERV,
  !                      indexed by process id (1:nmachs) and by variable 
  !                      idim_type


  subroutine LocalSizesAndDisp (ngrid, il1, ir2, jb1, jt2, localSize, disp, oneControlVars)

    integer, intent(in ) :: ngrid            ! which grid to use
    integer, intent(out) :: il1(nmachs)
    integer, intent(out) :: ir2(nmachs)
    integer, intent(out) :: jb1(nmachs)
    integer, intent(out) :: jt2(nmachs)
    integer, intent(out) :: localSize(nmachs, idim_type_min:idim_type_max)
    integer, intent(out) :: disp(nmachs, idim_type_min:idim_type_max)
    type(ControlVars), pointer, intent(in) :: oneControlVars

    integer :: proc
    integer :: idim_type
    character(len=8) :: c0
    character(len=*), parameter :: h="**(LocalSizesAndDisp)**" 

    ! field size at each process and dimensionality

    do proc = 1, nmachs
       localSize(proc,2) = nodemxp(proc,ngrid)*nodemyp(proc,ngrid)
       localSize(proc,3) = oneControlVars%nnzp*nodemxp(proc,ngrid)*nodemyp(proc,ngrid)
       localSize(proc,4) = nzg*nodemxp(proc,ngrid)*nodemyp(proc,ngrid)*npatch
       localSize(proc,5) = nzs*nodemxp(proc,ngrid)*nodemyp(proc,ngrid)*npatch
       localSize(proc,6) = nodemxp(proc,ngrid)*nodemyp(proc,ngrid)*npatch
       localSize(proc,7) = nodemxp(proc,ngrid)*nodemyp(proc,ngrid)*nwave
    end do

    ! displacement from base address of gathered 1D array 
    ! where to start gathering field portion from each process

    do idim_type = idim_type_min, idim_type_max
       disp(1,idim_type) = 0   ! first process at base address
       do proc = 2, nmachs   ! next process at (proc-1) position + localSize(proc-1)
          disp(proc,idim_type) = &
               disp(proc-1,idim_type)+localSize(proc-1,idim_type)
       end do
    end do

    ! gathered 1D array will be unpacked; in this process,
    ! eliminate parallel ghost zones, marking start and end
    ! indices of internal zones

    do proc = 1, nmachs
       if (btest(nodeibcon(proc,ngrid),0)) then
          il1(proc) = nodeia(proc,ngrid) - 1
       else
          il1(proc) = nodeia(proc,ngrid)
       end if
       if (btest(nodeibcon(proc,ngrid),1)) then
          ir2(proc) = nodeiz(proc,ngrid) + 1
       else
          ir2(proc) = nodeiz(proc,ngrid)
       end if
       if (btest(nodeibcon(proc,ngrid),2)) then
          jb1(proc) = nodeja(proc,ngrid) - 1
       else
          jb1(proc) = nodeja(proc,ngrid)
       end if
       if (btest(nodeibcon(proc,ngrid),3)) then
          jt2(proc) = nodejz(proc,ngrid) + 1
       else
          jt2(proc) = nodejz(proc,ngrid)
       end if
    end do
  end subroutine LocalSizesAndDisp



  ! PreProcAndGather: returns field gathered from all processes (at master process)
  !                   pre-process field before gathering if required.



  subroutine PreProcAndGather(preProc, ngrid, idim_type, varn, &
       il1, ir2, jb1, jt2, localSize, disp, thisChunkSize, LocalChunk, &
       sizeGathered, gathered, sizeFullField, FullField, &
       oneControlVars, oneBasicFields, oneTurbFields)

    logical,          intent(in   ) :: preProc
    integer,          intent(in   ) :: ngrid
    integer,          intent(in   ) :: idim_type
    character(len=*), intent(inout) :: varn
    integer,          intent(in   ) :: il1(nmachs)
    integer,          intent(in   ) :: ir2(nmachs)
    integer,          intent(in   ) :: jb1(nmachs)
    integer,          intent(in   ) :: jt2(nmachs)
    integer,          intent(in   ) :: localSize(nmachs, idim_type_min:idim_type_max)
    integer,          intent(in   ) :: disp(nmachs, idim_type_min:idim_type_max)
    integer,          intent(in   ) :: thisChunkSize
    real,             intent(inout) :: LocalChunk(thisChunkSize)
    integer,          intent(in   ) :: sizeGathered
    real,             intent(out  ) :: gathered(sizeGathered)    !scratch
    integer,          intent(in   ) :: sizeFullField
    real,             intent(out  ) :: FullField(sizeFullField)
    type(ControlVars), pointer, intent(in) :: oneControlVars
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields

    character(len=len(varn)) :: varnOut
    integer :: ierr
    character(len=7)            :: cProc
    character(len=8)            :: c0, c1, c2
    character(len=*), parameter :: h="**(PreProcAndGather)**" 


    write(cProc,"(a5,i2)") " Proc",mchnum

    ! if requiring pre-processing

    if (preProc) then

       if (dumpLocal) then
          write(c0,"(i8)") thisChunkSize
          write(*,"(a)") h//cProc//"  with thisChunkSize="//trim(adjustl(c0))
       end if

       ! pre-process LocalChunk before gathering

       select case (varn)
       case ('PP')

          ! Output total Exner function
          
          call RAMS_aprep_p (thisChunkSize, LocalChunk,  &
               oneBasicFields%pi0, LocalChunk)
          varnOut='PI'
          
       case ('HKM')
          
          ! Convert to HKM to HKH (note that VKH is HKH for Deardorff)
          
          call RAMS_aprep_hkh (thisChunkSize, LocalChunk, &
               oneTurbFields%vkh, oneBasicFields%dn0,  &
               LocalChunk, oneControlVars%idiffk, &
               oneControlVars%xkhkm)
          varnOut='HKH'
          
       case ('VKH')
          
          ! Un-density weight VKH
          
          call RAMS_aprep_vkh (thisChunkSize, LocalChunk, &
               oneBasicFields%dn0, LocalChunk)
          varnOut='VKH'
          
       case default
          varnOut=varn
       end select

       varn = varnOut

    end if

    ! gather field scaterred over slaves; master gets full unpacked field

    call GatherOneFullField(ngrid, idim_type, thisChunkSize, &
         LocalChunk, localSize, disp, il1, ir2, jb1, jt2, &
         sizeGathered, gathered, sizeFullField, FullField, oneControlVars)
  end subroutine PreProcAndGather






  subroutine GatherOneFullField(ngrid, idim_type, thisChunkSize, &
       LocalChunk, localSize, disp, il1, ir2, jb1, jt2, &
       sizeGathered, gathered, sizeFullField, FullField, oneControlVars)
    integer, intent(in ) :: ngrid
    integer, intent(in ) :: idim_type
    integer, intent(in ) :: thisChunkSize
    real,    intent(in ) :: LocalChunk(thisChunkSize)
    integer, intent(in ) :: localSize(nmachs, idim_type_min:idim_type_max)
    integer, intent(in ) :: disp(nmachs, idim_type_min:idim_type_max)
    integer, intent(in ) :: il1(nmachs)
    integer, intent(in ) :: ir2(nmachs)
    integer, intent(in ) :: jb1(nmachs)
    integer, intent(in ) :: jt2(nmachs)
    integer, intent(in ) :: sizeGathered
    real,    intent(out) :: gathered(sizeGathered)    !scratch
    integer, intent(in ) :: sizeFullField
    real,    intent(out) :: FullField(sizeFullField)
    type(ControlVars), pointer, intent(in) :: oneControlVars

    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(GatherOneFullField)**"

    ! gather field (with ghost zones and wrong order) at master_num

    if (dumpLocal) then
       write(c0,"(i8)") thisChunkSize
       write(*,"(a)") h//" will gather with local size "//trim(adjustl(c0))
       call flush(6)
    end if
    call GatherAllChunks (LocalChunk, thisChunkSize, idim_type, &
         localSize, disp, gathered, sizeGathered)
    if (dumpLocal) then
       write(*,"(a)") h//" done gathering"
    end if

    ! master_num unpacks fields (removes unnecessary ghost zones and positions entries)

    if (mchnum == master_num) then
       if (dumpLocal) then
          write(*,"(a)") h//" master will RemoveGhost"
       end if
       call RemoveGhost (ngrid, idim_type, sizeGathered, gathered, &
            il1, ir2, jb1, jt2, disp, localSize, sizeFullField, FullField, oneControlVars)
       if (dumpLocal) then
          write(*,"(a)") h//" done RemoveGhost"
       end if
    end if
  end subroutine GatherOneFullField






  subroutine GatherAllChunks (LocalChunk, thisChunkSize, idim_type, &
       localSize, disp, gathered, sizeGathered)
    integer, intent(in ) :: thisChunkSize
    real,    intent(in ) :: LocalChunk(thisChunkSize)
    integer, intent(in ) :: idim_type
    integer, intent(in ) :: localSize(nmachs, idim_type_min:idim_type_max)
    integer, intent(in ) :: disp(nmachs, idim_type_min:idim_type_max)
    integer, intent(in ) :: sizeGathered
    real,    intent(out) :: gathered(sizeGathered)

    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(GatherAllChunks)**" 

    if (dumpLocal) then
       write(c0,"(i8)") thisChunkSize
       write(c1,"(i8)") sizeGathered
       write(*,"(a)") h//" thisChunkSize="//trim(adjustl(c0))//&
            "; sizeGathered="//trim(adjustl(c1))
       call flush(6)
    end if

    ! gather a field

    call parf_GatherAllChunks(LocalChunk, thisChunkSize, idim_type, &
         localSize, disp, gathered, sizeGathered, master_num, nmachs)
    if (dumpLocal) then
       write(*,"(a)") h//" done"
       call flush(6)
    end if

  end subroutine GatherAllChunks







  subroutine RemoveGhost (ngrid, idim_type, sizeGathered, gathered, &
       il1, ir2, jb1, jt2, disp, localSize, sizeFullField, FullField, oneControlVars)
    integer, intent(in)  :: ngrid
    integer, intent(in)  :: idim_type
    integer, intent(in)  :: sizeGathered
    real,    intent(in)  :: gathered(sizeGathered)
    integer, intent(in)  :: il1(nmachs)
    integer, intent(in)  :: ir2(nmachs)
    integer, intent(in)  :: jb1(nmachs)
    integer, intent(in)  :: jt2(nmachs)
    integer, intent(in)  :: disp(nmachs, idim_type_min:idim_type_max)
    integer, intent(in)  :: localSize(nmachs, idim_type_min:idim_type_max)
    integer, intent(in)  :: sizeFullField
    real,    intent(out) :: FullField(sizeFullField)
    type(ControlVars), pointer, intent(in) :: oneControlVars

    integer :: proc
    character(len=8) :: c0
    character(len=*), parameter :: h="**(RemoveGhost)**"

    ! unpack gathered field, removing unnecessary ghost zones and
    ! placing entries at correct field positions

    select case (idim_type)
    case (2)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_2_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  oneControlVars%nnxp, oneControlVars%nnyp, &
                  nodemxp(proc,ngrid), nodemyp(proc,ngrid), &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (3)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_3_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  oneControlVars%nnzp, oneControlVars%nnxp, oneControlVars%nnyp, &
                  oneControlVars%nnzp, nodemxp(proc,ngrid), nodemyp(proc,ngrid), &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (4)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_4_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  nzg, oneControlVars%nnxp, oneControlVars%nnyp, npatch, &
                  nzg, nodemxp(proc,ngrid), nodemyp(proc,ngrid), npatch, &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (5)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_4_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  nzs, oneControlVars%nnxp, oneControlVars%nnyp, npatch, &
                  nzs, nodemxp(proc,ngrid), nodemyp(proc,ngrid), npatch, &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (6)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_2p_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  oneControlVars%nnxp, oneControlVars%nnyp, npatch, &
                  nodemxp(proc,ngrid), nodemyp(proc,ngrid), npatch, &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (7)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_buff_carma(FullField, gathered(disp(proc,idim_type)+1), &
                  oneControlVars%nnxp, oneControlVars%nnyp, nwave, &
                  nodemxp(proc,ngrid), nodemyp(proc,ngrid), nwave, &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" invoked with unknown idim_type ("//&
            trim(adjustl(c0))//")")
    end select
  end subroutine RemoveGhost



  ! RearrangeAndDump: dump field into all required output files;
  !                   rearrange field prior to output if required



  subroutine RearrangeAndDump (ngrid, rearran, &
       histFlag, instFlag, liteFlag, meanFlag, &
       histFileDS, instFileDS, liteFileDS, meanFileDS, &
       varn, idim_type, sizeFullField, FullField, Rear, oneControlVars)

    integer,          intent(in   ) :: ngrid
    logical,          intent(in   ) :: rearran
    logical,          intent(in   ) :: histFlag
    logical,          intent(in   ) :: instFlag
    logical,          intent(in   ) :: liteFlag
    logical,          intent(in   ) :: meanFlag
    type(IOFileDS),   intent(inout) :: histFileDS
    type(IOFileDS),   intent(inout) :: instFileDS
    type(IOFileDS),   intent(inout) :: liteFileDS
    type(IOFileDS),   intent(inout) :: meanFileDS
    character(len=*), intent(in   ) :: varn
    integer,          intent(in   ) :: idim_type
    integer,          intent(in   ) :: sizeFullField
    real,             intent(in   ) :: FullField(sizeFullField)
    real,             intent(out  ) :: Rear(sizeFullField)
    type(ControlVars), pointer, intent(in) :: oneControlVars

    character(len=*), parameter :: h="**(RearrangeAndDump)**"

    if (dumpLocal) then
       write(*,"(4(a,l1))") h//" enter field "//trim(varn)//":"//&
            " histFlag=",histFlag, &
            " instFlag=",instFlag, &
            " liteFlag=",liteFlag, &
            " meanFlag=",meanFlag
    end if

    ! if field to be rearranged, rearrange and dump

    if (rearran) then
       call RearrangeForOutput(oneControlVars%nnxp, oneControlVars%nnyp, oneControlVars%nnzp, &
            nzg, nzs, npatch, idim_type, FullField, Rear)

       ! dump rearranged field

       if (instFlag) then
          call FieldWriteStoreInfo(Rear, sizeFullField, &
               varn, idim_type, ngrid, instFileDS)
       end if
       if (liteFlag) then
          call FieldWriteStoreInfo(Rear, sizeFullField, &
               varn, idim_type, ngrid, liteFileDS)
       end if
       if (meanFlag) then
          call FieldWriteStoreInfo(Rear, sizeFullField, &
               varn, idim_type, ngrid, meanFileDS)
       end if
       if (histFlag) then
          call FieldWriteStoreInfo(Rear, sizeFullField, &
               varn, idim_type, ngrid, histFileDS)
       end if
    else

       ! dump original field

       if (instFlag) then
          call FieldWriteStoreInfo(FullField, sizeFullField, &
               varn, idim_type, ngrid, instFileDS)
       end if
       if (liteFlag) then
          call FieldWriteStoreInfo(FullField, sizeFullField, &
               varn, idim_type, ngrid, liteFileDS)
       end if
       if (meanFlag) then
          call FieldWriteStoreInfo(FullField, sizeFullField, &
               varn, idim_type, ngrid, meanFileDS)
       end if
       if (histFlag) then
          call FieldWriteStoreInfo(FullField, sizeFullField, &
               varn, idim_type, ngrid, histFileDS)
       end if
    end if
  end subroutine RearrangeAndDump



  ! Recreating Global Information (Gathering data)
  subroutine gatherData2D(idim_type, varn, ifm, nnxp, nnyp, &
       nmachs, mchnum, mynum, master_num,                   &
       localData2D, globalData2D, &
       oneControlVars, oneBasicFields, oneTurbFields)

    implicit none
    include "constants.h"
    ! Arguments:
    integer, intent(IN)           :: idim_type, ifm, nnxp, nnyp, &
         nmachs, mchnum, mynum, master_num
    character(LEN=16), intent(IN) :: varn
    real, intent(IN)              :: localData2D(:,:)
    real, intent(OUT)             :: globalData2D(:,:)
    type(ControlVars), pointer, intent(in) :: oneControlVars
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields

    ! Local Variables:
    character(LEN=16)  :: localVarn
    integer            :: ierr
    integer, parameter :: idim_type_min = 2
    integer, parameter :: idim_type_max = 7
    integer            :: il1(nmachs)
    integer            :: ir2(nmachs)
    integer            :: jb1(nmachs)
    integer            :: jt2(nmachs)
    integer            :: localSize(nmachs,idim_type_min:idim_type_max)
    integer            :: disp(nmachs,idim_type_min:idim_type_max)
    integer            :: maxLocalSize
    integer            :: sizeGathered(idim_type_min:idim_type_max)
    integer            :: maxSizeGathered
    integer            :: sizeFullField(idim_type_min:idim_type_max)
    integer            :: maxsizeFullField
    integer            :: globalSize(idim_type_min:idim_type_max)
    real, allocatable  :: localChunk(:)
    real, allocatable  :: gathered(:)
    real, allocatable  :: fullField(:)

    ! Recreating Global information about Soil Water
    ! grid dependent, field independent constants for gather and unpacking
    ! as a function of idim_type
    call LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp, oneControlVars)
    maxLocalSize = maxval(localSize(mynum,:))
    allocate(localChunk(maxLocalSize), stat=ierr)
    if (ierr/=0) then
       call fatal_error("Error allocating localChunk (gatherData)")
    endif
    call CopyLocalChunk(localData2D(1,1), localChunk, &
         LocalSize(mynum,idim_type))
    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
    maxSizeGathered = maxval(sizeGathered)
    allocate(gathered(maxSizeGathered), stat=ierr)
    if (ierr/=0) then
       call fatal_error("Error allocating gathered (gatherData)")
    endif
    ! grid dependent field sizes as a function of idim_type
    call GlobalSizes(ifm, nmachs, nwave, globalSize)
    if (mchnum==master_num) then
       sizeFullField(:) = globalSize(:)
    else
       sizeFullField(:) = 1
    end if
    maxSizeFullField = maxval(sizeFullField)
    allocate(fullField(sizeFullField(idim_type)), stat=ierr)
    if (ierr/=0) then
       call fatal_error("Error allocating fullField (gatherData)")
    endif
    localVarn = trim(varn)
    call PreProcAndGather(.false., ifm, idim_type, localVarn, &
         il1, ir2, jb1, jt2, localSize, disp,                 &
         localSize(mynum,idim_type), LocalChunk,              &
         sizeGathered(idim_type), gathered,                   &
         sizeFullField(idim_type), fullField, &
         oneControlVars, oneBasicFields, oneTurbFields)

    if (mchnum==master_num) then
       globalData2D = reshape(fullField, (/nnxp, nnyp/))
    endif

    call parf_bcast(globalData2D, int(nnxp,i8), int(nnyp,i8), master_num)

    deallocate(fullField)
    deallocate(gathered)
    deallocate(localChunk)

  end subroutine gatherData2D



  ! Recreating Global Information (Gathering data)
  subroutine gatherData3D(idim_type, varn, ifm, nnzp, nnxp, nnyp, &
       nmachs, mchnum, mynum, master_num,                         &
       localData3D, globalData3D, &
       oneControlVars, oneBasicFields, oneTurbFields)

    implicit none
    include "constants.h"
    ! Arguments:
    integer, intent(IN)           :: idim_type, ifm, nnzp, nnxp, nnyp, &
         nmachs, mchnum, mynum, master_num
    character(LEN=16), intent(IN) :: varn
    real, intent(IN)              :: localData3D(:,:,:)
    real, intent(OUT)             :: globalData3D(:,:,:)
    type(ControlVars), pointer, intent(in) :: oneControlVars
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields

    ! Local Variables:
    character(LEN=16)  :: localVarn
    integer            :: ierr
    integer, parameter :: idim_type_min = 2
    integer, parameter :: idim_type_max = 7
    integer            :: il1(nmachs)
    integer            :: ir2(nmachs)
    integer            :: jb1(nmachs)
    integer            :: jt2(nmachs)
    integer            :: localSize(nmachs,idim_type_min:idim_type_max)
    integer            :: disp(nmachs,idim_type_min:idim_type_max)
    integer            :: maxLocalSize
    integer            :: sizeGathered(idim_type_min:idim_type_max)
    integer            :: maxSizeGathered
    integer            :: sizeFullField(idim_type_min:idim_type_max)
    integer            :: maxsizeFullField
    integer            :: globalSize(idim_type_min:idim_type_max)
    real, allocatable  :: localChunk(:)
    real, allocatable  :: gathered(:)
    real, allocatable  :: fullField(:)

    ! Recreating Global information about Soil Water
    ! grid dependent, field independent constants for gather and unpacking
    ! as a function of idim_type
    call LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp, oneControlVars)
    maxLocalSize = maxval(localSize(mynum,:))
    allocate(localChunk(maxLocalSize), stat=ierr)
    if (ierr/=0) then
       call fatal_error("Error allocating localChunk (gatherData)")
    endif
    call CopyLocalChunk(localData3D(1,1,1), localChunk, &
         LocalSize(mynum,idim_type))
    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
    maxSizeGathered = maxval(sizeGathered)
    allocate(gathered(maxSizeGathered), stat=ierr)
    if (ierr/=0) then
       call fatal_error("Error allocating gathered (gatherData)")
    endif
    ! grid dependent field sizes as a function of idim_type
    call GlobalSizes(ifm, nmachs, nwave, globalSize)
    if (mchnum==master_num) then
       sizeFullField(:) = globalSize(:)
    else
       sizeFullField(:) = 1
    end if
    maxSizeFullField = maxval(sizeFullField)
    allocate(fullField(sizeFullField(idim_type)), stat=ierr)
    if (ierr/=0) then
       call fatal_error("Error allocating fullField (gatherData)")
    endif
    localVarn = trim(varn)
    call PreProcAndGather(.false., ifm, idim_type, localVarn, &
         il1, ir2, jb1, jt2, localSize, disp,                 &
         localSize(mynum,idim_type), LocalChunk,              &
         sizeGathered(idim_type), gathered,                   &
         sizeFullField(idim_type), fullField, &
         oneControlVars, oneBasicFields, oneTurbFields)

    if (mchnum==master_num) then
       globalData3D = reshape(fullField, (/nnzp, nnxp, nnyp/))
    endif

    call parf_bcast(globalData3D, &
         int(nnzp,i8), int(nnxp,i8), int(nnyp,i8), master_num)

    deallocate(fullField)
    deallocate(gathered)
    deallocate(localChunk)

  end subroutine gatherData3D



  ! Recreating Global Information (Gathering data)
  subroutine gatherData4D(idim_type, varn, ifm, mzg, nnxp, nnyp, npat, &
       nmachs, mchnum, mynum, master_num,                              &
       localData4D, globalData4D, &
       oneControlVars, oneBasicFields, oneTurbFields)

    implicit none
    include "constants.h"
    ! Arguments:
    integer, intent(IN)           :: idim_type, ifm, mzg, nnxp, nnyp, npat, &
         nmachs, mchnum, mynum, master_num
    character(LEN=16), intent(IN) :: varn
    real, intent(IN)              :: localData4D(:,:,:,:)
    real, intent(OUT)             :: globalData4D(:,:,:,:)
    type(ControlVars), pointer, intent(in) :: oneControlVars
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields

    ! Local Variables:
    character(LEN=16)  :: localVarn
    integer            :: ierr
    integer, parameter :: idim_type_min = 2
    integer, parameter :: idim_type_max = 7
    integer            :: il1(nmachs)
    integer            :: ir2(nmachs)
    integer            :: jb1(nmachs)
    integer            :: jt2(nmachs)
    integer            :: localSize(nmachs,idim_type_min:idim_type_max)
    integer            :: disp(nmachs,idim_type_min:idim_type_max)
    integer            :: maxLocalSize
    integer            :: sizeGathered(idim_type_min:idim_type_max)
    integer            :: maxSizeGathered
    integer            :: sizeFullField(idim_type_min:idim_type_max)
    integer            :: maxsizeFullField
    integer            :: globalSize(idim_type_min:idim_type_max)
    real, allocatable  :: localChunk(:)
    real, allocatable  :: gathered(:)
    real, allocatable  :: fullField(:)

    ! Recreating Global information about Soil Water
    ! grid dependent, field independent constants for gather and unpacking
    ! as a function of idim_type
    call LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp, oneControlVars)
    maxLocalSize = maxval(localSize(mynum,:))
    allocate(localChunk(maxLocalSize), stat=ierr)
    if (ierr/=0) then
       call fatal_error("Error allocating localChunk (gatherData)")
    endif
    call CopyLocalChunk(localData4D(1,1,1,1), localChunk, &
         LocalSize(mynum,idim_type))
    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
    maxSizeGathered = maxval(sizeGathered)
    allocate(gathered(maxSizeGathered), stat=ierr)
    if (ierr/=0) then
       call fatal_error("Error allocating gathered (gatherData)")
    endif
    ! grid dependent field sizes as a function of idim_type
    call GlobalSizes(ifm, nmachs, nwave, globalSize)
    if (mchnum==master_num) then
       sizeFullField(:) = globalSize(:)
    else
       sizeFullField(:) = 1
    end if
    maxSizeFullField = maxval(sizeFullField)
    allocate(fullField(sizeFullField(idim_type)), stat=ierr)
    if (ierr/=0) then
       call fatal_error("Error allocating fullField (gatherData)")
    endif
    localVarn = trim(varn)
    call PreProcAndGather(.false., ifm, idim_type, localVarn, &
         il1, ir2, jb1, jt2, localSize, disp,                 &
         localSize(mynum,idim_type), LocalChunk,              &
         sizeGathered(idim_type), gathered,                   &
         sizeFullField(idim_type), fullField, &
         oneControlVars, oneBasicFields, oneTurbFields)

    if (mchnum==master_num) then
       globalData4D = reshape(fullField, (/mzg, nnxp, nnyp, npat/))
    endif

    call parf_bcast(globalData4D, &
         int(mzg,i8), int(nnxp,i8), int(nnyp,i8), int(npat,i8), master_num)

    deallocate(fullField)
    deallocate(gathered)
    deallocate(localChunk)

  end subroutine gatherData4D

  !--(DMK-CCATT-INI)-----------------------------------------------------
  !temporary function 
  subroutine storeOwnChunk_3D(grid, fullGrid, toStore, nz, nx, ny, fieldName, oneControlVars)

    include "constants.h"

    integer, intent(in) :: grid
    integer, intent(in) :: nz
    integer, intent(in):: nx
    integer, intent(in):: ny
    real,  pointer, dimension(:,:,:):: toStore
    real,  dimension(nz,nx,ny):: fullGrid
    character(len=*), intent(in):: fieldName
    type(ControlVars), pointer, intent(in) :: oneControlVars
    
    integer :: ldimx
    integer:: ldimy
    integer:: lin
    integer :: ia
    integer:: iz
    integer:: ja
    integer:: jz
    integer:: ka
    integer:: kz
    integer :: ierr
    character(len=8):: c0
    character(len=8):: c1

    character(len=*), parameter :: h="**(storeOwnChunk_3D)**"

    ! check allocated memory

    if (.not. associated(toStore)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if

    lin = size(toStore,1)
    if (nz /= lin) then
       write(c0,"(i8)") nz
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" z dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if

    ldimx = nodemxp(mynum,grid)
    lin = size(toStore,2)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if

    ldimy = nodemyp(mynum,grid)
    lin = size(toStore,3)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if

    ia = nodei0(mynum,grid)+1
    iz = nodei0(mynum,grid)+nodemxp(mynum,grid)
    ja = nodej0(mynum,grid)+1
    jz = nodej0(mynum,grid)+nodemyp(mynum,grid)

    call mk_3_buff(fullGrid, toStore, &
         nz, oneControlVars%nnxp, oneControlVars%nnyp, nz, ldimx, ldimy, ia, iz, ja, jz)



  end subroutine storeOwnChunk_3D
  !--(DMK-CCATT-FIM)-----------------------------------------------------

  subroutine rams_aprep_p (n1,a,b,c)
    integer :: n1
    real :: a(*),b(*),c(*)

    integer :: i

    do i=1,n1
       c(i)=a(i)+b(i)
    enddo

    return
  end subroutine rams_aprep_p



  !******************************************************************************

  subroutine rams_aprep_hkh(n1,hkm,vkh,dn0,scr1,idiffk,xkhkm)
    integer :: n1,idiffk
    real :: xkhkm
    real, dimension(*) :: hkm,vkh,dn0,scr1
    integer :: ind

    if (idiffk <= 3) then
       do ind = 1,n1
          scr1(ind) = hkm(ind) * xkhkm / dn0(ind)
       enddo
    elseif (idiffk >= 4) then
       do ind = 1,n1
          scr1(ind) = vkh(ind) / dn0(ind)
       enddo
    endif

    return
  end subroutine rams_aprep_hkh

  !******************************************************************************

  subroutine rams_aprep_vkh(n1,vkh,dn0,vt3dd)
    integer :: n1
    real :: vkh(*),dn0(*),vt3dd(*)
    integer :: ind

    do ind = 1,n1
       vt3dd(ind) = vkh(ind) / dn0(ind)
    enddo

    return
  end subroutine rams_aprep_vkh


end module ReadBcst
