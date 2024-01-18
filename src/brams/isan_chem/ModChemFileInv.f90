!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModChemFileInv

  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModDateUtils, only: &
       date_add_to, &
       date_make_big

  use isan_coms, only: &
       igridfl, &
       iapr, &
       guess1st, &
       iarawi, &
       iasrfce, &
       isan_inc, &
       npdates, &
       iproc_dates, &
       i1st_flg, &
       iupa_flg, &
       isfc_flg, &
       iproc_flag, &
       iproc_names, &
       maxisfiles, &
       icfiletype, &
       icprefix

  use mem_grid, only: &
       time

  use dump, only: &
       dumpmessage  

  implicit none

  include "files.h" 
  include "constants.h"
  include "UseVfm.h"

  private

  public :: chem_isan_file_inv

contains




  subroutine chem_isan_file_inv (iyear1,imonth1,idate1,itime1,timmax,chem_assim, &
       chemistry)
    integer, intent(in) :: iyear1
    integer, intent(in) :: imonth1
    integer, intent(in) :: idate1
    integer, intent(in) :: itime1
    real, intent(in) :: timmax
    integer, intent(in) :: chem_assim
    integer, intent(in) :: chemistry

    ! Local variables ----------------------------------------------------------
    !
    !  nfgfiles = Number of First Guess FILES   (defined value = zero)
    !  nupfiles = Number of Upper Air FILES     (defined value = zero)
    !  nsffiles = Number of SurFace input FILES (defined value = zero)
    !
    ! -------------------------------------------------------------------------
    integer :: nfgfiles,nc,nf,lnf,nn,ndates,nupfiles,nsffiles,isan_err_flag
    integer :: inyear,inmonth,indate,inhour
    integer :: iyearf,imonthf,idatef,ihourf
    integer :: iyear2,imonth2,idate2,ihour2,ihour1
    real :: tinc
    character(len=14) :: itotdate_fg(maxisfiles),itotdate_up(maxisfiles)  &
         ,itotdate_sf(maxisfiles),itotdates(4*maxisfiles),idate_end
    character(len=f_name_length), dimension(maxisfiles) :: fnames_fg, &
         fnames_up, fnames_sf
    character(len=f_name_length) :: rams_filelist_arg
    logical :: there,fileExist
    integer :: localTime
    character(len=f_name_length) :: sVarName 
    integer :: iyears,imonths,idates,ihours
    integer :: indice
    character(len=3) :: suffix

    character(len=8) :: str(10)
    character(len=*),parameter :: h='**(chem_ISAN_file_inv)**'
    logical, parameter :: dumpLocal=.true.
    

    ! Define variables -------------------------------------------------------
    nfgfiles = 0
    nupfiles = 0
    nsffiles = 0
    ! ------------------------------------------------------------------------

    !          Go through first guess, upper air, surface input files
    !            and make inventory

    if(igridfl.ne.0) then
       if(iapr(1:1).ne.' '.and.iapr(1:1).ne.char(0) ) then

          print*,'----> performing the data inventory'

          nfgfiles=-1
          nc=len_trim(iapr)
          indice = 1

          nfgfiles = ((timmax/3600) / (isan_inc/100)) + 1    


          call date_add_to(iyear1,imonth1,idate1,itime1*100  &
               ,time,'s',iyears,imonths,idates,ihours)

          localTime = itime1
          indice = 1  

          do nf=1,nfgfiles

             call date_add_to(iyears, imonths, idates, localTime*100,  &
                  0., 's', iyears, imonths, idates, ihours)

             !--(DMK-CCATT-INI)----------------------------------------------------------------
             if(guess1st.eq.'PRESS')  then 
                if(CHEM_ASSIM == 1 .and. CHEMISTRY >= 0)  then 
                   write(sVarName,100) iapr(1:len_trim(iapr)),'-',iyears,'-',imonths,'-',idates,'-',ihours/100   
                   if (useVfm) then
                      sVarName = trim(sVarName) // '.vfm'
                   else
                      sVarName = trim(sVarName) // '.bin'
                   end if
                else 
                   write(sVarName,101) iapr(1:len_trim(iapr)),iyears,'-',imonths,'-',idates,'-',ihours/100   
                end if
             end if
100          format(a,a1,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)
101          format(a,   i4.4,a1,i2.2,a1,i2.2,a1,i4.4)


             if(guess1st.eq.'RAMS')   then 
                write(sVarName,102) iapr(1:len_trim(iapr)),'-A-',iyears,'-',imonths,'-',idates,'-',ihours   
                sVarName = sVarName // "-head.txt"
102             format(a,a3,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)
             end if

             inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

             if (there) then
                fnames_fg(indice) = trim(sVarName)
                indice = indice + 1
                if (dumpLocal) then
                   call MsgDump(h//" inserted "//trim(sVarName)//" at fnames_fg inventory")
                end if
             end if

             !--srf change LE 18000 to LE 2100 for  3/3 ivars
             if(localTime .le. 2100)then
                localTime = localTime + isan_inc
             else
                localTime = 000000
                localTime = localTime + isan_inc
             end if
          end do

          nfgfiles = indice - 1

          if(nfgfiles.gt.maxisfiles) then
             print*,'too many first guess files'
             call fatal_error('lots_of_first_guess')
          end if

          do nf=1,nfgfiles

             lnf=len_trim(fnames_fg(nf))

             suffix = trim(fnames_fg(nf)(len_trim(fnames_fg(nf))-2:len_trim(fnames_fg(nf))))

             if (guess1st.eq.'PRESS') then
                read(fnames_fg(nf)(lnf-18:lnf-4),20) inyear,inmonth,indate,inhour
             end if
             !--(DMK-CCATT-FIM)----------------------------------------------------------------

             ! form of a-A-2000-07-01-060000-head.txt
             if(guess1st.eq.'RAMS') then
                read(fnames_fg(nf)(lnf-25:lnf-10),20) inyear,inmonth,indate,inhour
             end if

20           format(i4,1x,i2,1x,i2,1x,i4)

             call date_make_big(inyear,inmonth,indate,inhour*100,itotdate_fg(nf))

          end do

          call RAMS_dintsort(nfgfiles,itotdate_fg,fnames_fg)
       end if
    end if


    if(igridfl.ne.3) then
       if(iarawi(1:1).ne.' '.and.iarawi(1:1).ne.char(0) ) then

          nupfiles=-1
          nc=len_trim(iarawi)

          nupfiles = ((timmax/3600) / (isan_inc/100)) + 1    

          call date_add_to(iyear1,imonth1,idate1,itime1*100  &
               ,time,'s',iyears,imonths,idates,ihours)

          localTime = itime1
          indice = 1  

          do nf=1,nupfiles

             call date_add_to(iyears, imonths, idates, localTime*100,  &
                  0., 's', iyears, imonths, idates, ihours)

             write(sVarName,400) iarawi(1:len_trim(iarawi)),iyears,'-',imonths,'-',idates,'-',ihours/100   
400          format(a,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)

             inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

             if (there) then
                fnames_up(indice) = trim(sVarName)
                indice = indice + 1
                if (dumpLocal) then
                   call MsgDump(h//" inserted "//trim(sVarName)//" at fnames_up inventory")
                end if
             end if

             if(localTime .le. 1800)then
                localTime = localTime + isan_inc
             else
                localTime = 000000
                localTime = localTime + isan_inc
             end if
          end do

          nupfiles = indice - 1

          if(nupfiles.gt.maxisfiles) then
             print*,'too many upper air files'
             call fatal_error('lots_of_upper_air')
          end if

          do nf=1,nupfiles
             lnf=len_trim(fnames_up(nf))
             read(fnames_up(nf)(lnf-14:lnf),20) inyear,inmonth,indate,inhour
             call date_make_big(inyear,inmonth,indate,inhour*100,itotdate_up(nf))
          end do
          call RAMS_dintsort(nupfiles,itotdate_up,fnames_up)

       end if


       if(iasrfce(1:1).ne.' '.and.iasrfce(1:1).ne.char(0) ) then

          nsffiles=-1
          nc=len_trim(iasrfce)
          nsffiles = ((timmax/3600) / (isan_inc/100)) + 1    

          call date_add_to(iyear1,imonth1,idate1,itime1*100  &
               ,time,'s',iyears,imonths,idates,ihours)

          localTime = itime1
          indice = 1  

          do nf=1,nsffiles

             call date_add_to(iyears, imonths, idates, localTime*100,  &
                  0., 's', iyears, imonths, idates, ihours)

             write(sVarName,300) iasrfce(1:len_trim(iasrfce)),iyears,'-',imonths,'-',idates,'-',ihours/100   
300          format(a,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)

             inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

             if (there) then
                fnames_sf(indice) = trim(sVarName)
                indice = indice + 1
                if (dumpLocal) then
                   call MsgDump(h//" inserted "//trim(sVarName)//" at fnames_sf inventory")
                end if
             end if

             if(localTime .le. 1800)then
                localTime = localTime + isan_inc
             else
                localTime = 000000
                localTime = localTime + isan_inc
             end if
          end do

          nsffiles = indice - 1

          if(nsffiles.gt.maxisfiles) then
             print*,'too many surface air files'
             call fatal_error('lots_of_surface')
          end if

          do nf=1,nsffiles
             lnf=len_trim(fnames_sf(nf))
             read(fnames_sf(nf)(lnf-14:lnf),20) inyear,inmonth,indate,inhour
             call date_make_big(inyear,inmonth,indate,inhour*100,itotdate_sf(nf))
          end do
          call RAMS_dintsort(nsffiles,itotdate_sf,fnames_sf)

       end if

    end if

    ! put dates in order, removing duplicates

    call RAMS_sort_dint3(nfgfiles,itotdate_fg  &
         ,nupfiles,itotdate_up,nsffiles,itotdate_sf  &
         ,ndates,itotdates)

    call RAMS_unique_dint(ndates,itotdates)

    !  start printing section
    !--------------------------------------------------------------

    print*,' '
    print*,' '
    print*,' '
    print*,'---------------------------------------------------------'
    write (*,fmt='(A,I1,A,I3,A)')'-- ISAN Input File Date Inventory for icFileType = ',icFileType,', with ',ndates,' dates.'
    print*,'---------------------------------------------------------'

    if (dumpLocal) then
       write(str(1),"(i8)") icFileType
       write(str(2),"(i8)") ndates
       call MsgDump(h//" ISAN Input File Date Inventory for icFileType="//trim(adjustl(str(1)))//&
            ", with "//trim(adjustl(str(2)))//" dates")
    end if
    
    do nn=1,ndates
       print*,'---- Date:',itotdates(nn)

       if(icFileType==1) then
          call makeGrib2fileName(trim(icPrefix),itime1,isan_inc,nn &
               ,iproc_names(nn,5),iproc_names(nn,6))
          print*,'---- First guess grib2 file:'  &
               ,iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))
          if (dumpLocal) then
             call MsgDump(h//'---- First guess grib2 file: '//trim(iproc_names(nn,5)))
          end if
       elseif(icFileType==2) then
          iproc_names(nn,5)=trim(icPrefix)
          print*,'---- First guess netCDF file:'  &
               ,iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))
          if (dumpLocal) then
             call MsgDump(h//'---- First guess netCDF file: '//trim(iproc_names(nn,5)))
          end if
       elseif(icFileType==3) then
          call makeGeosfileName(trim(icPrefix),isan_inc,nn &
               ,iproc_names(nn,5),iproc_names(nn,6),itotdates(nn),iyear1,imonth1,idate1,itime1)
          if (dumpLocal) then
             call MsgDump(h//'---- First guess Geos file: '//trim(iproc_names(nn,5))//&
                  " and "//trim(iproc_names(nn,6)))
          end if
       elseif(icFileType==4) then
          call makeGradsfileName(trim(icPrefix),isan_inc,nn &
               ,iproc_names(nn,5),iproc_names(nn,6),itotdates(nn),iyear1,imonth1,idate1,itime1)
          if (dumpLocal) then
             call MsgDump(h//'---- First guess grads file: '//trim(iproc_names(nn,5))//&
                  " and "//trim(iproc_names(nn,6)))
          end if
       else
          do nf=1,nfgfiles
             if(itotdates(nn).eq.itotdate_fg(nf)) then
                print*,'---- First guess file:'  &
                     ,fnames_fg(nf)(1:len_trim(fnames_fg(nf)))
                if (dumpLocal) then
                   call MsgDump(h//'---- First guess file: '//trim(fnames_fg(nf)))
                end if
             end if
          end do
       end if

       do nf=1,nupfiles
          if(itotdates(nn).eq.itotdate_up(nf)) then
             print*,'---- Upper air   file:'  &
                  ,fnames_up(nf)(1:len_trim(fnames_up(nf)))
             if (dumpLocal) then
                call MsgDump(h//'---- upper air file: '//trim(fnames_up(nf)))
             end if
          end if
       end do
       do nf=1,nsffiles
          if(itotdates(nn).eq.itotdate_sf(nf)) then
             print*,'---- Surface obs file:'  &
                  ,fnames_sf(nf)(1:len_trim(fnames_sf(nf)))
             if (dumpLocal) then
                call MsgDump(h//'---- surface obs file: '//trim(fnames_sf(nf)))
             end if
          end if
       end do
       print*,'------------------------------------------------------'
    end do
    !--------------------------------------------------------------

    ! Find dates we are going to process

    ! Find end date
    call date_add_to(iyear1,imonth1,idate1,itime1*100  &
         ,timmax,'s',iyearf,imonthf,idatef,ihourf)
    call date_make_big(iyearf,imonthf,idatef,ihourf,idate_end)

    ihour1=itime1*100
    tinc= (isan_inc/100) * 60.  + mod(isan_inc,100)
    npdates = 1
    call date_add_to (iyear1,imonth1,idate1,ihour1  &
         ,tinc*(npdates-1),'m',iyear2,imonth2,idate2,ihour2)
    call date_make_big(iyear2,imonth2,idate2,ihour2  &
         ,iproc_dates(npdates))
    do while (iproc_dates(npdates) .lt. idate_end)
       npdates = npdates + 1

       call date_add_to (iyear1,imonth1,idate1,ihour1  &
            ,tinc*(npdates-1),'m',iyear2,imonth2,idate2,ihour2)
       call date_make_big(iyear2,imonth2,idate2,ihour2  &
            ,iproc_dates(npdates))
    end do


    !   We have the dates to process. Find if files exist and set overall
    !     "go ahead" flag. Put filenames in iproc_names array if they will be used.

    !      iproc_flag array info:
    !      1) Process?     0=no, 1=yes
    !      2) First guess? 0=no file, 1=exists, 2=exists, don't use, 3=interpolate
    !      3) Upper air?   0=no file, 1=exists, 2=exists, don't use
    !      4) Surface?     0=no file, 1=exists, 2=exists, don't use

    isan_err_flag=0

    print*,' '
    print*,' '
    print*,' '
    print*,'---------------------------------------------------------'
    print*,'-----------  ISAN Processing Information    -------------'
    print*,'---------------------------------------------------------'
    print*,'----    Flags:  IGRIDFL =',igridfl,'  I1ST_FLG=',i1st_flg
    print*,'----    Flags:  IUPA_FLG=',iupa_flg,'  ISFC_FLG=',isfc_flg
    print*,'---------------------------------------------------------'

    do nn=1,npdates
       iproc_flag(nn,1)=0
       iproc_flag(nn,2)=0
       iproc_flag(nn,3)=0
       iproc_flag(nn,4)=0

       print*,'------------------------------------------------------'
       print*,'---- Date:', iproc_dates(nn)

       iproc_flag(nn,1)=1

       if(icFileType==1) then
          print *,itime1; call flush(6)
          call makeGrib2fileName(trim(icPrefix),itime1,isan_inc,nn &
               ,iproc_names(nn,5),iproc_names(nn,6))
          print*,'---- First guess grib2 (NCEP) file:'  &
               ,iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))
          if (dumpLocal) then
             call MsgDump(h//'---- First guess grib2 file: '//trim(iproc_names(nn,5)))
          end if
       elseif(icFileType==2) then
          iproc_names(nn,5)=trim(icPrefix)
          print*,'---- First guess nc4 (NASA) file:'  &
               ,iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))       
          if (dumpLocal) then
             call MsgDump(h//'---- First guess nc4 file: '//trim(iproc_names(nn,5)))
          end if
       elseif(icFileType==3) then
          call makeGeosfileName(trim(icPrefix),isan_inc,nn &
               ,iproc_names(nn,5),iproc_names(nn,6),iproc_dates(nn),iyear1,imonth1,idate1,itime1)
          if (dumpLocal) then
             call MsgDump(h//'---- First guess geos file: '//trim(iproc_names(nn,5))//&
                  " and "//trim(iproc_names(nn,6)))
          end if
       elseif(icFileType==4) then
          call makeGradsfileName(trim(icPrefix),isan_inc,nn &
               ,iproc_names(nn,5),iproc_names(nn,6),iproc_dates(nn),iyear1,imonth1,idate1,itime1)
          if (dumpLocal) then
             call MsgDump(h//'---- First guess grads file: '//trim(iproc_names(nn,5))//&
                  " and "//trim(iproc_names(nn,6)))
          end if
       end if

       if(icFileType==1) then

          inquire(file=iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))), exist=fileExist )

          if(.not. fileExist) then
             iErrNumber=dumpMessage(c_tty,c_yes,h,'468' &
                  ,c_fatal,'File '//iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))) &
                  //' not found. Please, verify and solve it!')
          end if
          write(*,fmt='(A)') '---- First guess file' &
               //iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))//' exists.'
          iproc_flag(nn,2)=1
       elseif(icFileType==2) then

          inquire(file=iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))), exist=fileExist )

          if(.not. fileExist)  iErrNumber=dumpMessage(c_tty,c_yes,h,'468' &
               ,c_fatal,'File '//iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))) &
               //' not found. Please, verify and solve it!')
          write(*,fmt='(A)') '---- First guess file' &
               //iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))//' exists.'
          iproc_flag(nn,2)=1
       elseif(icFileType==3) then

          inquire(file=iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))), exist=fileExist )

          if(.not. fileExist)  iErrNumber=dumpMessage(c_tty,c_yes,h,'468' &
               ,c_fatal,'File '//iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))) &
               //' not found. Please, verify and solve it!')
          write(*,fmt='(A)') '---- First guess file' &
               //iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))//' exists.'
          iproc_flag(nn,2)=1
       elseif(icFileType==4) then

          inquire(file=iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))), exist=fileExist )
          
          if(.not. fileExist)  iErrNumber=dumpMessage(c_tty,c_yes,h,'468' &
               ,c_fatal,'File '//iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))) &
               //' not found. Please, verify and solve it!')
          write(*,fmt='(A)') '---- First guess file' &
               //iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))//' exists.'
          iproc_flag(nn,2)=1
       else

          ! first guess
          iproc_flag(nn,2)=0
          do nf=1,nfgfiles
             if(iproc_dates(nn).eq.itotdate_fg(nf)) then
                print*,'---- First guess file exists.'
                if(igridfl.eq.0) then
                   iproc_flag(nn,2)=2
                else
                   iproc_flag(nn,2)=1
                end if
                iproc_names(nn,1)=fnames_fg(nf)

                if (dumpLocal) then
                   call MsgDump(h//" First guess file "//trim(iproc_names(nn,1))//" included")
                end if

                goto 71
             end if
          end do
          ! file doesn't exist
          if(i1st_flg.eq.1) then
             print*,'---- First guess file does not exist.'  &
                  ,' Will not process this time.'
             iproc_flag(nn,1)=0
             if (dumpLocal) then
                call MsgDump(h//" First guess file does not exist; will not process this time")
             end if
          elseif(i1st_flg.eq.2) then
             isan_err_flag=1
             print*,'---- First guess file does not exist.'  &
                  ,' Will STOP run.'
             iproc_flag(nn,1)=0
             if (dumpLocal) then
                call MsgDump(h//" First guess file does not exist; will stop run")
             end if
          elseif(i1st_flg.eq.3) then
             print*,'---- First guess file does not exist.'  &
                  ,' Will attempt interpolation.'
             iproc_flag(nn,2)=3
             if (dumpLocal) then
                call MsgDump(h//" First guess file does not exist; will attemp interpolation")
             end if
          end if

       end if
71     continue

       ! upper air
       iproc_flag(nn,3)=0
       do nf=1,nupfiles
          if(iproc_dates(nn).eq.itotdate_up(nf)) then
             print*,'---- Upper air file exists.'
             if(igridfl.eq.3) then
                iproc_flag(nn,3)=2
             else
                iproc_flag(nn,3)=1
             end if
             iproc_names(nn,2)=fnames_up(nf)
             if (dumpLocal) then
                call MsgDump(h//" Upper air file "//trim(adjustl(iproc_names(nn,2)))//" should exist")
             end if
             goto 72
          end if
       end do

       ! file doesn't exist
       if(iupa_flg.eq.1) then
          print*,'---- Upper air file does not exist.'  &
               ,' Will not process this time.'
          iproc_flag(nn,1)=0
          if (dumpLocal) then
             call MsgDump(h//" Upper air file does not exist; will not process this time")
          end if
       elseif(iupa_flg.eq.2) then
          isan_err_flag=1
          print*,'---- Upper air file does not exist.',' Will STOP run.'
          iproc_flag(nn,1)=0
          if (dumpLocal) then
             call MsgDump(h//" Upper air file does not exist; will stop run")
          end if
       elseif(iupa_flg.eq.3) then
          print*,'---- Upper air file does not exist.',' Will try running without.'
          iproc_flag(nn,3)=3
          if (dumpLocal) then
             call MsgDump(h//" Upper air file does not exist; will try running without")
          end if
       end if
72     continue

       ! surface
       iproc_flag(nn,4)=0
       do nf=1,nsffiles
          if(iproc_dates(nn).eq.itotdate_sf(nf)) then
             print*,'---- Surface obs file exists.'
             if(igridfl.eq.3) then
                iproc_flag(nn,4)=2
             else
                iproc_flag(nn,4)=1
             end if
             iproc_names(nn,3)=fnames_sf(nf)
             if (dumpLocal) then
                call MsgDump(h//" Surface obs file "//trim(adjustl(iproc_names(nn,3)))//" should exist")
             end if
             goto 73
          end if
       end do

       ! file doesn't exist
       if(isfc_flg.eq.1) then
          print*,'---- Surface file does not exist.',' Will not process this time.'
          iproc_flag(nn,1)=0
          if (dumpLocal) then
             call MsgDump(h//" Surface file does not exist; will not process this time")
          end if
       elseif(isfc_flg.eq.2) then
          isan_err_flag=1
          print*,'---- Surface file does not exist.',' Will STOP run.'
          iproc_flag(nn,1)=0
          if (dumpLocal) then
             call MsgDump(h//" Surface file does not exist; will stop run")
          end if
       elseif(isfc_flg.eq.3) then
          print*,'---- Surface file does not exist.',' Will try running without.'
          iproc_flag(nn,4)=3
          if (dumpLocal) then
             call MsgDump(h//" Surface file does not exist; will try running without")
          end if
       end if
73     continue

    end do
    print*,'---------------------------------------------------------'
    if (isan_err_flag.ne.0) then
       print*,'ISAN run stopping because of errors!'
       print*,'See previous output listing.'
       call fatal_error('isan_file_errors')
    end if
  end subroutine chem_isan_file_inv






  subroutine makeGrib2fileName(prefix,itime1,isan_inc,nn,gribFileName,invFileName)
    ! Input/Output variables
    character(len=*), intent(in) :: prefix 
    integer,intent(in)    :: itime1
    integer,intent(in)    :: isan_inc
    integer,intent(in)    :: nn
    character(len=*),intent(out) :: gribFileName
    character(len=*),intent(out) :: invFileName
    !#
    !Local variables
    character(len=3) :: ctime

    character(len=*),parameter :: h='**(makeGrib2fileName)**'
    logical, parameter :: dumpLocal=.true.

    !Code
    !print *,'In: ',itime1,nn,isan_inc
    write(ctime,fmt='(I3.3)') int(itime1/100)+(nn-1)*isan_inc/100
    gribFileName=trim(prefix)//ctime
    invFileName=trim(prefix)//ctime//'.inv'


  end subroutine makeGrib2fileName



  subroutine makeGradsfileName(prefix,isan_inc,nn,gradsFileName,invFileName,&
       iproc_date,iyear1,imonth1,idate1,itime1)

    ! Input/Output variables
    character(len=*), intent(in) :: prefix 
    integer,intent(in)    :: isan_inc
    integer,intent(in)    :: nn
    character(len=*),intent(out) :: gradsFileName
    character(len=*),intent(out) :: invFileName
    character(len=*), intent(in) :: iproc_date 
    integer,intent(in)    :: iyear1
    integer,intent(in)    :: imonth1
    integer,intent(in)    :: idate1
    integer,intent(in)    :: itime1

    !Local variables
    character(len=2) :: ctime,cdate,cmonth
    character(len=4) :: cyear!,cntime

    character(len=*),parameter :: h='**(makeGradsfileName)**'
    logical, parameter :: dumpLocal=.true.


    !Code
    !print *,'In: nn,isan_inc,iproc_date,iyear1,imonth1,idate1,itime1= ',nn,isan_inc,iproc_date,iyear1,imonth1,idate1,itime1
    write(ctime,fmt='(I2.2)') int(itime1/100)
    write(cyear,fmt='(I4.4)') iyear1
    write(cmonth,fmt='(I2.2)') imonth1
    write(cdate,fmt='(I2.2)') idate1
    gradsFileName=trim(prefix)//iproc_date(1:10)//'.gra'
    invFileName=trim(prefix)//ctime//'.inv'

  end subroutine makeGradsfileName






  subroutine makeGeosfileName(prefix,isan_inc,nn,geosFileName,invFileName,&
       iproc_date,iyear1,imonth1,idate1,itime1)

    ! Input/Output variables
    character(len=*), intent(in) :: prefix 
    integer,intent(in)    :: isan_inc
    integer,intent(in)    :: nn
    character(len=*),intent(out) :: geosFileName
    character(len=*),intent(out) :: invFileName
    character(len=*), intent(in) :: iproc_date 
    integer,intent(in)    :: iyear1
    integer,intent(in)    :: imonth1
    integer,intent(in)    :: idate1
    integer,intent(in)    :: itime1

    !Local variables
    character(len=2) :: ctime,cdate,cmonth
    character(len=4) :: cyear!,cntime

    character(len=*), parameter :: h='**(makeGeosfileName)**'
    logical, parameter :: dumpLocal=.true.
    
    !Code
    print *,'In: nn,isan_inc,iproc_date,iyear1,imonth1,idate1,itime1= ',nn,isan_inc,iproc_date,iyear1,imonth1,idate1,itime1
    ! GEOS.fp.fcst.inst3_3d_asm_Cp.20200213_00+20200214_0300.V01.nc4
    !write(cntime,fmt='(I4.4)') int(itime1)+(nn-1)*isan_inc
    write(ctime,fmt='(I2.2)') int(itime1/100)
    write(cyear,fmt='(I4.4)') iyear1
    write(cmonth,fmt='(I2.2)') imonth1
    write(cdate,fmt='(I2.2)') idate1
    geosFileName=trim(prefix)//cyear//cmonth//cdate//'_'//ctime//'+'//iproc_date(1:8)//'_'//iproc_date(9:12)//'.V01.nc4'
    invFileName=trim(prefix)//ctime//'.inv'

  end subroutine makeGeosfileName
end module ModChemFileInv
