!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModMkSfcSst

  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModNestFillDens, only: &
       fillscr, &
       fillvar

  use ModRUser, only: &
       sstinit_user

  use mem_leaf, only: &
       seatmp

  use io_params,only: &
       isstflg,        &
       isstfn, &
       sstfpfx

  use mem_mksfc, only: &
       sfcfile_p, &
       mksfc_scr1, &
       mksfc_scr2, &
       nvsstf, &
       iyearvs, &
       imonthvs, &
       idatevs, &
       ihourvs, &
       mksfc_vt2da, &
       mksfc_vt2db, &
       scrx,           &
       vsstfil 

  use modgeodat, only: &
       geodat

  use grid_dims, only :&
       nxpmax, &
       nypmax

  use mem_grid, only: &
       platn, &
       plonn, &
       xtn, &
       ytn, &
       deltaxn, &
       deltayn, &
       nxtnest, &
       nnxp, &
       nnyp, &
       iyear1, &
       imonth1, &
       idate1, &
       ihour1

  implicit none

  include "files.h"
  include "UseVfm.h"

  integer, parameter :: fUnit=25
  private

  public :: sst_read_dataheader
  public :: sstnest
  public :: sst_write

contains

  subroutine sst_read_dataheader(ifm)
    integer :: ifm

    integer :: itime,nc
    character(len=f_name_length) :: flnm,line,line2
    character(len=1) :: dummy
    logical :: there
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(sst_read_dataheader)**"
    logical, parameter :: dumpLocal=.true.

    ! Read header file for all sstdata files (all times and locations).  The header
    ! file contains:
    ! first line:       geographic block size (degrees), file size, geographic
    !                   starting point for the dataset (south pole), and offsets
    ! second line:      number of data times (NDTIM)
    ! next NDTIM lines: file prefix, year, month, day, and hour for one data time
    ! last 3 lines:     comments describing the first, second, and following lines
    !                   above

    ! Construct header file name

    flnm=trim(isstfn(ifm))//'HEADER'

    inquire(file=flnm(1:len_trim(flnm)),exist=there)

    if (dumpLocal) then
       if (there) then
          call MsgDump(h//" found file "//trim(flnm))
       else
          call MsgDump(h//" did not find file "//trim(flnm))
       end if
    end if
          
    if (.not.there) then
       print*,'SSTDATA header file:', trim(flnm), ' for grid ',ifm,' not there.'
       stop 'sst_read_fileheader-1'
    endif

    ! Read this header file

    call rams_f_open(fUnit,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
    rewind fUnit

    ! read number of data times in dataset

    read(fUnit,*) dummy
    read(fUnit,*) nvsstf(ifm)

    if (nvsstf(ifm) <= 0) then
       print*, 'No SST input files found with specified prefix or incorrect header'
       close(fUnit)
       stop 'sst_read_fileheader-2'
    endif

    ! read prefix list and times 

    do itime = 1,nvsstf(ifm)
       read(fUnit,'(A80)') line
       call char_strip_var(line, flnm, line2)
       read (line2,*) iyearvs(itime,ifm), imonthvs(itime,ifm),  &
            idatevs(itime,ifm), ihourvs(itime,ifm)

       vsstfil(itime,ifm)=trim(isstfn(ifm))//trim(flnm)

       if (dumpLocal) then
          write(str(1),"(i8)") ihourvs(itime,ifm)
          write(str(2),"(i8)") idatevs(itime,ifm)
          write(str(3),"(i8)") imonthvs(itime,ifm)
          write(str(4),"(i8)") iyearvs(itime,ifm)
          call MsgDump(h//" inserted file "//trim(adjustl(vsstfil(itime,ifm)))//&
               " of "//trim(adjustl(str(1)))//":"//&
               trim(adjustl(str(2)))//"/"//&
               trim(adjustl(str(3)))//"/"//&
               trim(adjustl(str(4)))//" at mksfc file inventory")
       end if
    enddo

    close(fUnit)

    return
  end subroutine sst_read_dataheader

  !******************************************************************************

  subroutine sstnest(ifm, ivtime)
    ! Arguments:
    integer, intent(IN) :: ifm, ivtime
    ! Local Variables:
    integer :: icm

    icm = nxtnest(ifm)

    ! Initialize SEATP and SEATF in subroutine sstinit

    call sstinit(nnxp(ifm),nnyp(ifm),ifm, sfcfile_p(ifm)%seatf)

    if (icm>=1 .and. isstflg(ifm)==0) then

       ! Interpolate SEATF from coarser grid

       call fillscr(1, nxpmax, nypmax, 1, nnxp(icm), nnyp(icm), 1, 1,  &
            mksfc_scr1, sfcfile_p(icm)%seatf)
       call fillvar(1, nxpmax, nypmax, 1, nnxp(ifm), nnyp(ifm), 1, 1,  &
            mksfc_scr2, sfcfile_p(ifm)%seatf)

       nvsstf(ifm)                 = nvsstf(icm)
       iyearvs (1:nvsstf(ifm),ifm) = iyearvs (1:nvsstf(ifm),icm)
       imonthvs(1:nvsstf(ifm),ifm) = imonthvs(1:nvsstf(ifm),icm)
       idatevs (1:nvsstf(ifm),ifm) = idatevs (1:nvsstf(ifm),icm)
       ihourvs (1:nvsstf(ifm),ifm) = ihourvs (1:nvsstf(ifm),icm)

    elseif (isstflg(ifm)==1) then

       ! Interpolate SEATF from standard dataset

       call geodat(nnxp(ifm), nnyp(ifm), sfcfile_p(ifm)%seatf,     &
            isstfn(ifm)(1:len_trim(isstfn(ifm))), vsstfil(ivtime,ifm)(1:len_trim(vsstfil(ivtime,ifm))), mksfc_vt2da, mksfc_vt2db, &
            ifm, 'SST')

    else

       iyearvs(1,ifm) = iyear1; imonthvs(1,ifm) = imonth1
       idatevs(1,ifm) = idate1; ihourvs(1,ifm)  = ihour1        

    endif

    ! If desired, override current values of SEATF with user-defined
    ! changes to subroutine sstinit_user.

    call sstinit_user(nnxp(ifm), nnyp(ifm), ifm, sfcfile_p(ifm)%seatf)

  end subroutine sstnest

  ! ****************************************************************************

  subroutine sstinit(n2,n3,ifm,seatf)
    integer :: n2,n3,ifm,i,j
    real, dimension(n2,n3) :: seatf

    ! Fill the SEATF array with a default value of seatmp.  This 
    ! default is used only when a standard RAMS sst dataset is not used and when 
    ! no overrides to sea temperature are defined in subroutine sstinit_user 
    ! in the file ruser.f90.

    do j = 1,n3
       do i = 1,n2
          seatf(i,j) = seatmp
       enddo
    enddo
    return
  end subroutine sstinit


  !****************************************************************************

  subroutine sst_write(ifm,ivt)
    integer, intent(in) :: ifm
    integer, intent(in) :: ivt

    integer :: i,j
    integer :: ios
    
    real :: glatr,glonr
    character(len=f_name_length) :: flnm
    character(len=2) :: cgrid
    character(len=*), parameter :: h="**(sst_write)**"
    logical, parameter :: dumpLocal=.true.

    ! Write sst data to sst file for one grid and one time

    write(cgrid,'(a1,i1)') 'g',ifm

    if (useVfm) then
       call makefnam(flnm,sstfpfx,0.,iyearvs(ivt,ifm),imonthvs(ivt,ifm) &
            ,idatevs(ivt,ifm),ihourvs (ivt,ifm)*10000,'W',cgrid,'vfm')
    else
       call makefnam(flnm,sstfpfx,0.,iyearvs(ivt,ifm),imonthvs(ivt,ifm) &
            ,idatevs(ivt,ifm),ihourvs (ivt,ifm)*10000,'W',cgrid,'bin')
    end if
    if (dumpLocal) then
       call MsgDump(h//" build fileName "//trim(flnm))
    end if
       
    call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

    if (useVfm) then

       call rams_f_open(fUnit,flnm(1:len_trim(flnm)),'FORMATTED','REPLACE','WRITE',1)
       rewind fUnit
       
       write(fUnit,99) 999999,2
       write(fUnit,100) iyearvs(ivt,ifm),imonthvs(ivt,ifm) &
            ,idatevs(ivt,ifm),ihourvs (ivt,ifm)
       write(fUnit,101) nnxp(ifm),nnyp(ifm)
       write(fUnit,102) deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr
       
99     format(2i8)
100    format(1x,i4.4,2(1x,i2.2),1x,i4.4)
101    format(4i5)
102    format(6f16.5)
       
       call vforec(fUnit,sfcfile_p(ifm)%seatf,nnxp(ifm)*nnyp(ifm),24,scrx,'LIN')

       close(fUnit)


    else

       open(fUnit, action="write", file=trim(flnm), form="unformatted", iostat=ios)
       if (ios /= 0) then
          call fatal_error(h//" opening file "//trim(flnm))
       end if
       rewind fUnit
       
       write(fUnit) 999999,2
       write(fUnit) iyearvs(ivt,ifm),imonthvs(ivt,ifm) &
            ,idatevs(ivt,ifm),ihourvs (ivt,ifm)
       write(fUnit) nnxp(ifm),nnyp(ifm)
       write(fUnit) deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr
       
       write(fUnit) sfcfile_p(ifm)%seatf
       
       close(fUnit)
    end if

    if (dumpLocal) then
       call MsgDump(h//" wrote field seatf at file "//trim(flnm))
    end if
    return
  end subroutine sst_write
end module ModMkSfcSst
