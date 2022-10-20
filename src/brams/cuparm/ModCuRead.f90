!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModCuRead

  use ModNamelistFile, only: &
       NamelistFile

  use ModCuParmVars, only: &
       CuParmVars
  
  use ModCuParmFields, only: &
       CuParmFields

  use ModDateUtils, only: &
       date_abs_secs2,   &
       date_add_to,      &
       date_make_big

  use mem_grid, only: &
       nnzp, &
       nnxp, &
       nnyp, &
       nxtnest, &
       grid_g, &
       time, &
       timmax , &
       ngrids, &
       idate1,     &
       imonth1,    &
       itime1,     &
       iyear1,     &
       runtype

  use isan_coms, only: &
       isan_inc

  use grid_dims, only: &
       maxfiles          ! INTENT(IN)

  implicit none

  include "files.h"
  include "constants.h"

  private
  public :: cu_read

contains



  subroutine cu_read(initflag, oneNamelistFile, oneCuParmVars, oneCuParmFields)


    !------------------------------------------------------
    !  Read cumulus inversion tendencies
    !------------------------------------------------------
    integer, intent(in) :: initflag
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(CuParmVars), pointer, intent(in) :: oneCuParmVars
    type(CuParmFields), pointer, intent(in) :: oneCuParmFields

    character(len=14)  :: itotdate_start
    integer :: iyears,imonths,idates,ihours,nf,ifm

    if (initflag == 1) then   ! Initialization

       ! Inventory all cu inversion files. 
       call cu_file_inv (iyear1,imonth1,idate1,itime1, oneNamelistFile, oneCuParmVars)

       ! Find past time file

       if (runtype == 'HISTORY') then
          call date_add_to(iyear1,imonth1,idate1,itime1*100  &
               ,max(time,oneNamelistFile%tcu_beg),'s',iyears,imonths,idates,ihours)
          call date_make_big(iyears,imonths,idates,ihours  &
               ,itotdate_start)

       elseif (runtype == 'INITIAL') then
          call date_add_to(iyear1,imonth1,idate1,itime1*100  &
               ,oneNamelistFile%tcu_beg,'s',iyears,imonths,idates,ihours)
          call date_make_big(iyears,imonths,idates,ihours  &
               ,itotdate_start)
       endif

       do nf=1,oneCuParmVars%ncufiles
          if(itotdate_start >= oneCuParmVars%itotdate_cu(nf) .and.  &
               itotdate_start <  oneCuParmVars%itotdate_cu(nf+1) ) then
             oneCuParmVars%ncufl=nf
             exit
          endif
       enddo

       print*,'nud starting at history file:',oneCuParmVars%ncufl

       ! Read initial files.

       call cu_update(0,oneCuParmVars%ncufl,oneNamelistFile, oneCuParmVars, oneCuParmFields)


    elseif (initflag == 2) then   ! Runtime file increment

       if ( time >= oneCuParmVars%cu_times(oneCuParmVars%ncufl+1) ) oneCuParmVars%ncufl = oneCuParmVars%ncufl + 1

    endif


    ! Read new files.

    call cu_update(1,oneCuParmVars%ncufl+1,oneNamelistFile,oneCuParmVars, oneCuParmFields)


    oneCuParmVars%cutime1=oneCuParmVars%cu_times(oneCuParmVars%ncufl)
    oneCuParmVars%cutime2=oneCuParmVars%cu_times(oneCuParmVars%ncufl+1)

    return
  end subroutine cu_read


  subroutine cu_file_inv (iyear1,imonth1,idate1,itime1, oneNamelistFile, oneCuParmVars)
    integer, intent(in) :: iyear1
    integer, intent(in) :: imonth1
    integer, intent(in) :: idate1
    integer, intent(in) :: itime1
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(CuParmVars), pointer, intent(in) :: oneCuParmVars

    integer :: nc,nf,lnf,nhftot
    integer :: inyear,inmonth,indate,inhour


    character(len=f_name_length), dimension(maxfiles) :: fnames
    character(len=f_name_length) :: rams_filelist_arg
    character(len=14)  :: itotdate
    real(kind=8) :: secs_init,secs_cu

    logical there 
    integer :: localTime
    character(len=f_name_length) :: sVarName 
    integer :: iyears,imonths,idates,ihours
    integer :: indice,ng
    character(len=1)  :: cgrid

    ! Get abs seconds of run start

    call date_abs_secs2(iyear1,imonth1,idate1,itime1*100,secs_init)

    ! Go through and make inventory

    nhftot=-1

    nhftot = ((timmax/3600) / (isan_inc/100)) + 1    

    call date_add_to(iyear1,imonth1,idate1,itime1*100  &
         ,time,'s',iyears,imonths,idates,ihours)

    localTime = itime1
    indice = 1  


    do ng=1,ngrids

       write (cgrid,'(i1)') ng

       do nf=1,nhftot
          call date_add_to(iyears, imonths, idates, localTime*100,  &
               0., 's', iyears, imonths, idates, ihours)

          write(sVarName,100) trim(oneNamelistFile%cu_prefix), &
               iyears,'-',imonths,'-',idates,'-',ihours/100   
100       format(a,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)

          inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

          if (there) then
             fnames(indice) = trim(sVarName)//"g"//cgrid//".vfm"
             indice = indice + 1
          endif

          if(localTime .LE. 1800)then
             localTime = localTime + isan_inc
          else
             localTime = 000000
             localTime = localTime + isan_inc
          endif
       enddo

    enddo

    nhftot = indice - 1


    if(nhftot > maxfiles) then
       print*,'too many cu files'
       stop 'lots_of_cu_files'
    endif

    oneCuParmVars%ncufiles=0
    do nf=1,nhftot

       ! only save grid 1 files names and times. 

       if (index(fnames(nf),'-g1.') /= 0) then
          lnf=len_trim(fnames(nf))
          read(fnames(nf)(lnf-23:lnf-7),20) inyear,inmonth,indate,inhour
20        format(i4,1x,i2,1x,i2,1x,i6)

          call date_make_big(inyear,inmonth,indate,inhour,itotdate)

          oneCuParmVars%ncufiles=oneCuParmVars%ncufiles+1
          oneCuParmVars%fnames_cu(oneCuParmVars%ncufiles)=fnames(nf)
          oneCuParmVars%itotdate_cu(oneCuParmVars%ncufiles)=itotdate

          call date_abs_secs2(inyear,inmonth,indate,inhour,secs_cu)
          oneCuParmVars%cu_times(oneCuParmVars%ncufiles)=secs_cu - secs_init
       endif

    enddo

    call RAMS_dintsort(oneCuParmVars%ncufiles,oneCuParmVars%itotdate_cu,oneCuParmVars%fnames_cu)

    !  start printing section
    !--------------------------------------------------------------

    print*,' '
    print*,' '
    print*,' '
    print*,'-------------------------------------------------------------'
    print*,'-----------  Cumulus Tendency Input File Inventory -------------'
    print*,'-------------------------------------------------------------'
    do nf=1,oneCuParmVars%ncufiles
       print*,  oneCuParmVars%itotdate_cu(nf),'   ',oneCuParmVars%cu_times(nf)  &
            ,trim(oneCuParmVars%fnames_cu(nf))
    enddo
    print*,'------------------------------------------------------'

    !--------------------------------------------------------------


    return
  end subroutine cu_file_inv

  !******************************************************************************

  subroutine cu_update(iswap,ncu, oneNamelistFile, oneCuParmVars, oneCuParmFields)
    integer, intent(in) :: iswap
    integer, intent(in) :: ncu
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(CuParmVars), pointer, intent(in) :: oneCuParmVars
    type(CuParmFields), pointer, intent(in) :: oneCuParmFields

    include "files.h"

    character (len=f_name_length) :: cunamein
    character (len=1) :: cng
    integer :: ngr,nc,ifm,icm !,npts
    integer(kind=i8) :: npts

    integer,save :: iun=10
    logical :: there

    ! Put new fields into cu future arrays. If iswap == 1, 
    !     swap future into past first

    if (iswap == 1) then
       do ngr=1,ngrids
          oneCuParmFields%thsrcp(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
               oneCuParmFields%thsrcf(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
          oneCuParmFields%rtsrcp(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))=  &
               oneCuParmFields%rtsrcf(1:nnzp(ngr),1:nnxp(ngr),1:nnyp(ngr))
          oneCuParmFields%conprrp(1:nnxp(ngr),1:nnyp(ngr))=  &
               oneCuParmFields%conprrf(1:nnxp(ngr),1:nnyp(ngr))
       enddo
    endif


    ! Open the input file and read fields


    do ngr=1,ngrids

       ifm=ngr
       icm=nxtnest(ifm)

       print*,'ncu:',ncu,oneCuParmVars%fnames_cu(ncu)
       nc=len_trim(oneCuParmVars%fnames_cu(ncu))
       write(cng,'(i1)') ngr
       cunamein=oneCuParmVars%fnames_cu(ncu)(1:nc-5)//cng//'.vfm'
       print*,'ncu:',ncu,cunamein

       inquire (file=cunamein(1:len_trim(cunamein)), exist=there)

       if (oneNamelistFile%wt_cu_grid(ngr) > 0.) then
          if (there) then

             call rams_f_open(iun,cunamein(1:len_trim(cunamein)),'FORMATTED','OLD','READ',0)

             npts=nnzp(ngr)*nnxp(ngr)*nnyp(ngr)
             call vfirec(iun,oneCuParmFields%thsrcf,npts,'LIN')
             call vfirec(iun,oneCuParmFields%rtsrcf,npts,'LIN')

             npts=nnxp(ngr)*nnyp(ngr)
             call vfirec(iun,oneCuParmFields%conprrf,npts,'LIN')
          endif
       endif

    enddo

    ! Close the input file

    close(iun)

    return
  end subroutine cu_update


end module ModCuRead
