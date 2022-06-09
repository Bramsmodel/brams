! # rotinas de escrita e leitura baseadas em trabalho de DMK.

module digitalFilter

  use ModControlVars, only: &
       ControlVars
  
  use ModBasicFields, only: &
       BasicFields
  
  use node_mod, only: &
       mchnum, &
       nmachs, &
       mynum, &
       master_num

  use mem_grid, only:    &
       iyear1, &
       imonth1, &
       idate1, &
       itime1, &
       time, &
       begtime,&
       istp, &
       dtlongn,&
       ngrids

  use grid_dims, only: &
       maxgrds

  use var_tables, only: &
       num_var, &
       vtab_r

  use ReadBcst, only: &
       LocalSizesAndDisp

  use ModNamelistFile, only: &
       namelistFile

  use ModDateUtils, only: &
       date_add_to

  use io_params,only: &
       hfilout

  implicit none

  private
  public :: StoreNamelistFileAtdigitalFilter
  public :: initDigitalFilter
  public :: applyDigitalFilter
  public :: fileNameDF
  public :: dfVars
  public :: applyDF
  public :: timeWindowDF
  public :: iposDF
  public :: frqanlDF

  include "files.h"

  type df_vars

     real, dimension(:,:,:), pointer :: UP
     real, dimension(:,:,:), pointer :: VP
     real, dimension(:,:,:), pointer :: WP
     real, dimension(:,:,:), pointer :: PP
     real, dimension(:,:,:), pointer :: UC
     real, dimension(:,:,:), pointer :: VC
     real, dimension(:,:,:), pointer :: WC
     real, dimension(:,:,:), pointer :: PC
     real, dimension(:,:,:), pointer :: THP
     real, dimension(:,:,:), pointer :: RTP
     real, dimension(:,:,:), pointer :: THETA
     real, dimension(:,:,:), pointer :: RV

  end type df_vars

  ! # ramsin <--
  logical:: applyDF
  real:: timeWindowDF

  ! # local
  character(len=f_name_length):: fileNameDF
  type(df_vars), pointer, dimension(:):: dfVars
  real, dimension(:), pointer:: dfHK
  integer:: timestepDF
  real:: writeTimeDF
  integer:: writeStpDF
  real:: frqanlDF
  integer:: iposDF


contains


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine StoreNamelistFileAtdigitalFilter(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile

    applyDF = oneNamelistFile%applyDigitalFilter
    timeWindowDF = oneNamelistFile%digitalFilterTimeWindow
  end subroutine StoreNamelistFileAtdigitalFilter



  subroutine initDigitalFilter(dfVars, ngrids, nz, nx, ny)

    type(df_vars), dimension(:), pointer, intent(out):: dfVars
    integer, intent(in):: ngrids
    integer, dimension(ngrids), intent(in):: nz
    integer, dimension(ngrids), intent(in):: nx
    integer, dimension(ngrids), intent(in):: ny


    integer :: ng


    if(applyDF)then

       allocate(dfvars(ngrids))

       writeTimeDF = 0
       writeStpDF = 0

       do ng = 1, ngrids

          allocate(dfVars(ng)%UP(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%VP(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%WP(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%PP(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%UC(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%VC(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%WC(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%PC(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%THP(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%RTP(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%THETA(nz(ng), nx(ng), ny(ng)))
          allocate(dfVars(ng)%RV(nz(ng), nx(ng), ny(ng)))

          dfVars(ng)%UP = 0.0
          dfVars(ng)%VP = 0.0
          dfVars(ng)%WP = 0.0
          dfVars(ng)%PP = 0.0
          dfVars(ng)%UC = 0.0
          dfVars(ng)%VC = 0.0
          dfVars(ng)%WC = 0.0
          dfVars(ng)%PC = 0.0
          dfVars(ng)%THP= 0.0
          dfVars(ng)%RTP = 0.0
          dfVars(ng)%THETA= 0.0
          dfVars(ng)%RV= 0.0
       end do
    end if

  end subroutine initDigitalFilter

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine destroyDigitalFilter(ngrids)

    integer, intent(in) :: ngrids

    integer :: ng


    deallocate(dfHK)

    do ng = 1, ngrids
       deallocate(dfVars(ng)%UP)
       deallocate(dfVars(ng)%VP)
       deallocate(dfVars(ng)%WP)
       deallocate(dfVars(ng)%PP)
       deallocate(dfVars(ng)%UC)
       deallocate(dfVars(ng)%VC)
       deallocate(dfVars(ng)%WC)
       deallocate(dfVars(ng)%PC)
       deallocate(dfVars(ng)%THP)
       deallocate(dfVars(ng)%RTP)
       deallocate(dfVars(ng)%THETA)
       deallocate(dfVars(ng)%RV)
    end do

    deallocate(dfvars)

  end subroutine destroyDigitalFilter

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine applyDigitalFilter(fileName, dfVars, oneBasicFields, oneControlVars)

    ! # parameter.
    character(len=*), intent(inout) :: fileName
    type(df_vars), dimension(:), pointer, intent(inout):: dfVars
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(ControlVars), pointer, intent(in) :: oneControlVars
    
    ! # local.
    integer:: ng
    integer:: idx, k, i
    integer :: n, nsteps,nhlf
    real:: tetac,timecut, TC,DT,thetac,wc,arg
    real, parameter:: pi = 3.1415927 !there is no such constant in rconstants.f90
    real, parameter        :: filter = 1

    if(.not. applyDF) return

    if(time .eq. 0.0)then

       if(filter == 1) then      

          timecut = timeWindowDF/2.0
          !timecut = 6.*3600.

          timestepDF = 0

          n = (timeWindowDF*0.5)/dtlongn(1)

          allocate(dfHK((n*2)+1)); dfHK=0.0

          !tetac= 2*pi*dtlongn(1)/(timeWindowDF/2.0)
          tetac= 2*pi*dtlongn(1)/(timecut)

          do idx = -n, n

             if (idx /= 0) then     !DSM  para nao dar divisao por zero
                dfHK(idx+n+1) = ( (sin(idx*tetac)) / (idx*pi)) * sin((idx*pi)/(n+1)) / ((idx*pi)/(n+1))
             else
                dfHK(n+1) = 2.0 * dtlongn(1) / (timecut)
             endif

          end do

          !do idx =1, n*2+1
          !  if(mchnum == master_num) print*,"DFHK=",idx,dfhk(idx),sum(dfhk)
          !enddo
          !dfHK(n+1) = 2.0 * dtlongn(1) / (timeWindowDF/2.0)

          ! checking dfhk normalization
          dfHK = dfHK/sum(dfhk)
          !IF(MYnum==1) print*,'SUM of dig filter = ',sum(dfhk)

       else

          TC=timeWindowDF/2.0
          dt=dtlongn(1) 
          thetac = TC ! cutoff
          nsteps=TC/DT+1
          nhlf = (nsteps+1)/2
          allocate(dfHK((nhlf*2)+1)); dfHK=0.0
          do k = 1, nhlf-1
             n   = k-nhlf
             arg = n*pi/nhlf            
             wc  = sin(arg)/arg ! Lanczos window
             dfHK(k) = wc*sin(n*2.0*pi*DT/tc)/(n*pi)
          end do
          dfHK(nhlf) = 2*DT/TC
          do i = nhlf+1, nsteps
             dfHK(i) = dfHK(nsteps-i+1)
          end do
          dfHK = dfHK/sum(dfHK)
          do idx =1, nhlf*2+1
             if(mchnum == master_num) print*,"DFHK=",idx,dfhk(idx),sum(dfhk)
          enddo
       endif

    end if

    !IF(MYnum==1) print*,' filtering = ',time,istp,dfHK(istp)

    do ng = 1, ngrids
       dfVars(ng)%UP=  dfVars(ng)%UP    + (dfHK(istp) * oneBasicFields%UP)
       dfVars(ng)%VP=  dfVars(ng)%VP    + (dfHK(istp) * oneBasicFields%VP)
       dfVars(ng)%WP=  dfVars(ng)%WP    + (dfHK(istp) * oneBasicFields%WP)
       dfVars(ng)%PP=  dfVars(ng)%PP    + (dfHK(istp) * oneBasicFields%PP)
       dfVars(ng)%UC=  dfVars(ng)%UC    + (dfHK(istp) * oneBasicFields%UC)
       dfVars(ng)%VC=  dfVars(ng)%VC    + (dfHK(istp) * oneBasicFields%VC)
       dfVars(ng)%WC=  dfVars(ng)%WC    + (dfHK(istp) * oneBasicFields%WC)
       dfVars(ng)%PC=  dfVars(ng)%PC    + (dfHK(istp) * oneBasicFields%PC)
       dfVars(ng)%THP=  dfVars(ng)%THP   + (dfHK(istp) * oneBasicFields%THP)
       dfVars(ng)%RTP=  dfVars(ng)%RTP   + (dfHK(istp) * oneBasicFields%RTP)
       dfVars(ng)%THETA=  dfVars(ng)%THETA + (dfHK(istp) * oneBasicFields%THETA)
       dfVars(ng)%RV=  dfVars(ng)%RV    + (dfHK(istp) * oneBasicFields%RV)
    end do

    if(time .ne. 0.0 .and. mod(begtime+dtlongn(1), (timeWindowDF/2.0)) .lt. dtlongn(1) .and. writeTimeDF .eq. 0)then
       writeTimeDF = begtime
       writeStpDF = istp
       if(MYNUM==1)then
          ng=1
          print*,'================================================================'
          print*, 'saving fields @time(h)=',time,begtime
          print*, 'max min U V = ', maxval(oneBasicFields%UP),minval(oneBasicFields%UP)&
               ,maxval(oneBasicFields%VP),minval(oneBasicFields%vP)
          call flush(6)
          print*,'================================================================'
       endif

       call saveNodeFields_digFilt(fileName, oneControlVars)
    end if

    if(time .ne. 0.0 .and. mod(begtime, timeWindowDF) .lt. dtlongn(1) .and. writeTimeDF .ne. 0)then

       call loadNodeFields_digFilt(fileName, oneControlVars)
       if(MYNUM==1)then
          ng=1
          print*,'================================================================'
          print*, 'loading fields @time(h)=',time,begtime
          print*, 'max min U V = ', maxval(oneBasicFields%UP),minval(oneBasicFields%UP)&
               ,maxval(oneBasicFields%VP),minval(oneBasicFields%vP)
          call flush(6)
       endif

       do ng = 1, ngrids
          oneBasicFields%UP       = dfVars(ng)%UP
          oneBasicFields%VP       = dfVars(ng)%VP
          oneBasicFields%WP       = dfVars(ng)%WP
          oneBasicFields%PP       = dfVars(ng)%PP
          oneBasicFields%UC       = dfVars(ng)%UC
          oneBasicFields%VC       = dfVars(ng)%VC
          oneBasicFields%WC       = dfVars(ng)%WC
          oneBasicFields%PC= dfVars(ng)%PC
          oneBasicFields%THP= dfVars(ng)%THP
          oneBasicFields%RTP = dfVars(ng)%RTP
          oneBasicFields%THETA= dfVars(ng)%THETA
          oneBasicFields%RV       = dfVars(ng)%RV
       end do

       if(MYNUM==1)then
          ng=1
          print*, 'filtering fields @time(h)=',time,begtime
          print*, 'max min U V = ', maxval(oneBasicFields%UP),minval(oneBasicFields%UP)&
               ,maxval(oneBasicFields%VP),minval(oneBasicFields%vP)
          print*,'================================================================'
          call flush(6)
       endif
       call destroyDigitalFilter(ngrids)

       begtime = writeTimeDF
       istp = writeStpDF - 1
       applyDF = .false.

    end if

  end subroutine applyDigitalFilter

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !--(DMK)-------------------------------------------------------------------

  subroutine saveNodeFields_digFilt(fileName, oneControlVars)

    character(len=*), intent(out) :: fileName
    type(ControlVars), pointer, intent(in) :: oneControlVars

    integer :: idim_type

    integer, parameter :: idim_type_min=2
    integer, parameter :: idim_type_max=7
    integer :: il1(nmachs)
    integer :: ir2(nmachs)
    integer :: jb1(nmachs)
    integer :: jt2(nmachs)
    integer :: localSize(nmachs,idim_type_min:idim_type_max)
    integer :: disp(nmachs,idim_type_min:idim_type_max)
    integer :: maxLocalSize

    real, allocatable :: LocalChunk(:)

    integer :: ng, nv

    character(len=7)            :: cProc
    character(len=16)           :: varn
    character(len=8)            :: c0, c1
    character(len=*), parameter :: h="**(saveNodeFields)**"
    logical, parameter          :: dumpLocal = .false.
    integer                     :: ierr

    write(cProc,"(a5,i2)") " Proc",mchnum

    do ng = 1, ngrids

       if (dumpLocal) then
          write(c0,"(i8)") ng
          write(*,"(a)") h//cProc//" starts processing grid "//trim(adjustl(c0))
          call flush(6)
       end if

       ! grid dependent, field independent constants for gather and unpacking
       ! as a function of idim_type
       call LocalSizesAndDisp(ng, il1, ir2, jb1, jt2, localSize, disp, oneControlVars)

       ! maximum sizes over all fields
       maxLocalSize = maxval(LocalSize(mynum,:))

       ! maximum space for pre-processed local field
       allocate(LocalChunk(maxLocalSize), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") maxLocalSize
          write(c1,"(i8)") ierr
          call fatal_error(h//" allocate LocalSize("// &
               trim(adjustl(c0))//") failed with stat="// &
               trim(adjustl(c1)))
       else if (dumpLocal) then
          write(c0,"(i8)") maxLocalSize
          write(*,"(a)") h//cProc//" allocated LocalSize("//trim(adjustl(c0))//")"
          call flush(6)
       end if

       call OpenNodeWrite_digFilt(259, time, iyear1, imonth1, idate1, itime1, maxgrds, ng, mchnum, nmachs, fileName)
       ! for all fields on this grid

       if (dumpLocal) then
          write(*, "(4(a,l1))") h//cProc//" File Open for IOUTPUT=3"
       endif

       do nv = 1, num_var(ng)

          ! dimensionality
          idim_type = vtab_r(nv,ng)%idim_type
          if (idim_type < idim_type_min .or. idim_type > idim_type_max) then
             write(c0,"(i8)") idim_type
             call fatal_error(h//" unknown idim_type="//trim(adjustl(c0)))
          end if

          ! case 1: output current field values (for hist, inst and lite output files)

          if (idim_type == 2) then
             call CopyLocalChunk(vtab_r(nv,ng)%var_p_2D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 3) then
             call CopyLocalChunk(vtab_r(nv,ng)%var_p_3D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 4) then
             call CopyLocalChunk(vtab_r(nv,ng)%var_p_4D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 5) then
             call CopyLocalChunk(vtab_r(nv,ng)%var_p_4D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 6) then
             call CopyLocalChunk(vtab_r(nv,ng)%var_p_3D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 7) then
             call CopyLocalChunk(vtab_r(nv,ng)%var_p_3D, LocalChunk, LocalSize(mynum,idim_type))
          end if
          call nodeWrite_digFilt(259, LocalChunk, LocalSize(mynum,idim_type))
       enddo

       deallocate(LocalChunk, stat=ierr)
       if (ierr /= 0) then
          write(c1,"(i8)") ierr
          call fatal_error(h//" deallocate LocalSize failed with stat="//trim(adjustl(c1)))
       end if

       call CloseNodeWrite_digFilt(259)

    end do ! # do ng = 1, ngrids

  end subroutine saveNodeFields_digFilt
  !--(DMK)-------------------------------------------------------------------
  subroutine loadNodeFields_digFilt(fileName, oneControlVars)

    character(len=*), intent(in) :: fileName
    type(ControlVars), pointer, intent(in) :: oneControlVars

    integer :: idim_type

    integer, parameter :: idim_type_min=2
    integer, parameter :: idim_type_max=7
    integer :: il1(nmachs)
    integer :: ir2(nmachs)
    integer :: jb1(nmachs)
    integer :: jt2(nmachs)
    integer :: localSize(nmachs,idim_type_min:idim_type_max)
    integer :: disp(nmachs,idim_type_min:idim_type_max)
    integer :: maxLocalSize

    real, allocatable :: LocalChunk(:)

    integer :: ng, nv

    character(len=7)            :: cProc
    character(len=16)           :: varn
    character(len=8)            :: c0, c1
    character(len=*), parameter :: h="**(saveNodeFields)**"
    logical, parameter          :: dumpLocal = .false.
    integer                     :: ierr

    write(cProc,"(a5,i2)") " Proc",mchnum

    do ng = 1, ngrids

       if (dumpLocal) then
          write(c0,"(i8)") ng
          write(*,"(a)") h//cProc//" starts processing grid "//trim(adjustl(c0))
          call flush(6)
       end if

       ! grid dependent, field independent constants for gather and unpacking
       ! as a function of idim_type
       call LocalSizesAndDisp(ng, il1, ir2, jb1, jt2, localSize, disp, oneControlVars)

       ! maximum sizes over all fields
       maxLocalSize = maxval(LocalSize(mynum,:))

       ! maximum space for pre-processed local field
       allocate(LocalChunk(maxLocalSize), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") maxLocalSize
          write(c1,"(i8)") ierr
          call fatal_error(h//" allocate LocalSize("// &
               trim(adjustl(c0))//") failed with stat="// &
               trim(adjustl(c1)))
       else if (dumpLocal) then
          write(c0,"(i8)") maxLocalSize
          write(*,"(a)") h//cProc//" allocated LocalSize("//trim(adjustl(c0))//")"
          call flush(6)
       end if

       open (unit=259, file=filename(1:len_trim(filename)), form='unformatted', iostat=ierr)
       ! for all fields on this grid

       if(dumpLocal) then
          write(*, "(4(a,l1))") h//cProc//" File Open for IOUTPUT=3"
       endif

       do nv = 1, num_var(ng)

          ! dimensionality
          idim_type = vtab_r(nv,ng)%idim_type
          if (idim_type < idim_type_min .or. idim_type > idim_type_max) then
             write(c0,"(i8)") idim_type
             call fatal_error(h//" unknown idim_type="//trim(adjustl(c0)))
          end if

          ! case 1: output current field values (for hist, inst and lite output files)

          call nodeRead_digFilt(259, LocalChunk)

          if (idim_type == 2) then
             call CopyLocalChunkReverse(vtab_r(nv,ng)%var_p_2D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 3) then
             call CopyLocalChunkReverse(vtab_r(nv,ng)%var_p_3D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 4) then
             call CopyLocalChunkReverse(vtab_r(nv,ng)%var_p_4D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 5) then
             call CopyLocalChunkReverse(vtab_r(nv,ng)%var_p_4D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 6) then
             call CopyLocalChunkReverse(vtab_r(nv,ng)%var_p_3D, LocalChunk, LocalSize(mynum,idim_type))
          else if (idim_type == 7) then
             call CopyLocalChunkReverse(vtab_r(nv,ng)%var_p_3D, LocalChunk, LocalSize(mynum,idim_type))
          end if
       enddo

       deallocate(LocalChunk, stat=ierr)
       if (ierr /= 0) then
          write(c1,"(i8)") ierr
          call fatal_error(h//" deallocate LocalSize failed with stat="//trim(adjustl(c1)))
       end if

       call CloseNodeWrite_digFilt(259)

    end do ! # do ng = 1, ngrids

  end subroutine loadNodeFields_digFilt
  !--(DMK)-------------------------------------------------------------------

  subroutine OpenNodeWrite_digFilt(iUnit, time, iyear1, imonth1, idate1, &
       itime1, maxgrds, ng, mchnum, nmachs, fileName)

    integer, intent(in) :: iUnit
    real,    intent(in) :: time
    integer, intent(in) :: iyear1
    integer, intent(in) :: imonth1
    integer, intent(in) :: idate1
    integer, intent(in) :: itime1
    integer, intent(in) :: maxgrds
    integer, intent(in) :: ng
    integer, intent(in) :: mchnum
    integer, intent(in) :: nmachs
    character(len=*), intent(out) :: fileName

    integer            :: ierr
    integer            :: oyear1
    integer            :: omonth1
    integer            :: odate1
    integer            :: otime1
    character(len=4)   :: yy
    character(len=2)   :: mm
    character(len=2)   :: dd
    character(len=6)   :: hh
    character(len=5)   :: pp
    character(len=5)   :: nn
    character(len=2)   :: gg
    character(len=8)   :: c0
!!$  !character(len=256), parameter :: lfsafilout='/tmp/massaru'
!!$  character(len=256), parameter :: lfsafilout='/tmp'
    character(len=7)            :: cProc
    character(len=*), parameter :: h="**(OpenNodeWrite)**"
    logical, parameter          :: dumpLocal = .false.

    write(cProc,"(a5,i2)") " Proc",mchnum


    call date_add_to(iyear1,imonth1,idate1,itime1,time,'s', &
         oyear1,omonth1,odate1,otime1)

    write (yy,'(i4.4)') oyear1
    write (mm,'(i2.2)') omonth1
    write (dd,'(i2.2)') odate1
    write (hh,'(i6.6)') otime1
    write (pp,'(i5.5)') mchnum
    write (nn,'(i5.5)') nmachs
    write (gg,'(i1)')   ng

    filename=trim(hfilout)//"-p"//pp//"-"//nn//"-A"//yy//"-"//mm//"-"//&
         dd//"-"//hh//"-g"//trim(gg)//"-bin-df.tmp"

    if (dumpLocal) then
       write(*, "(4(a,l1))") h//cProc//" filename="//filename
    endif

    open (iUnit, file=filename(1:len_trim(filename)), form='unformatted', &
         iostat=ierr)

    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error("opening "//trim(filename)//" failed with iostat "//&
            trim(adjustl(c0)))
    end if

  end subroutine OpenNodeWrite_digFilt
  !--(DMK)-------------------------------------------------------------------


  !--(DMK)-------------------------------------------------------------------
  subroutine CloseNodeWrite_digFilt(iUnit)
    integer, intent(in) :: iUnit

    close (iUnit)

    return
  end subroutine CloseNodeWrite_digFilt
  !--(DMK)-------------------------------------------------------------------



  !--(DMK)-------------------------------------------------------------------
  subroutine nodeWrite_digFilt(iUnit, field, LocalSize)
    integer, intent(in) :: iUnit
    real,    intent(in) :: field(LocalSize)
    integer, intent(in) :: LocalSize

    call writeField_digFilt(iUnit, field, LocalSize)

    return
  end subroutine nodeWrite_digFilt
  !--(DMK)-------------------------------------------------------------------



  !--(DMK)-------------------------------------------------------------------
  subroutine writeField_digFilt(iUnit, field, LocalSize)
    integer, intent(in) :: iUnit
    real,    intent(in) :: field(LocalSize)
    integer, intent(in) :: LocalSize

    write(iUnit) LocalSize
    write(iUnit) field(1:LocalSize)

    return
  end subroutine writeField_digFilt

  !--(DMK)-------------------------------------------------------------------
  subroutine nodeRead_digFilt(iUnit, field)
    integer, intent(in) :: iUnit
    real,    intent(inout) :: field(*)

    call readField_digFilt(iUnit, field)

    return
  end subroutine nodeRead_digFilt
  !--(DMK)-------------------------------------------------------------------



  !--(DMK)-------------------------------------------------------------------
  subroutine readField_digFilt(iUnit, field)
    integer, intent(in) :: iUnit
    real,    intent(inout) :: field(*)

    integer :: LocalSize

    read(iUnit) LocalSize
    read(iUnit) field(1:LocalSize)

    return
  end subroutine readField_digFilt

  !-----------

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module digitalFilter
