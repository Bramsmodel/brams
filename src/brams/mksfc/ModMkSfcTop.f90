!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModMkSfcTop

  use ModParallelEnvironment, only: &
       MsgDump

  use ModControlVars, only: &
       ControlVars
  
  use mem_grid, only: &
       platn, &
       plonn, &
       xtn, &
       ytn, &
       deltaxn, &
       deltayn, &
       nnxp, &
       nnyp, &
       nnxyp, &
       grid_g, &
       oneGlobalGridData

  use io_params, only: &
       topfiles, &
       itoptflg, &
       itopsflg, &
       toptenh, &
       toptwvl, &
       iz0flg, &
       z0max, &
       z0fact

  use mem_mksfc, only: &
       sfcfile_p, &
       scrx

  use dump,only: &
       dumpMessage

  use node_mod, only:      &
       mchnum,             &
       master_num, &
       nodemxp, &
       mynum

  use ReadBcst, only: &
       ReadStoreFullFieldAndOwnChunk, &
       ReadStoreOwnChunk


  implicit none

  include "files.h"
  include "constants.h"
  include "UseVfm.h"
  
  integer, parameter :: fUnit=25

  private
  public :: top_check
  public :: toptinit
  public :: top_write
  public :: TopReadStoreOwnChunk
  public :: TRSFFieldAndOwnChunk


contains


  
  subroutine top_check(ifm,ierr)
    ! This subroutine checks for the existence of a surface file for
    ! grid number ifm, and if it exists, also checks for agreement of
    ! grid configuration between the file and the current model run.
    ! If the file does not exist or does not match grid configuration,
    ! the flag ifileok is returned with a value of 0.  If the file
    ! exists and is ok, ifileok is returned with a value of 1.
    integer, intent(in) :: ifm
    integer, intent(out) :: ierr

    integer :: lc,isfc_marker,isfc_ver,nsfx,nsfy  &
         ,nsitoptflg,nsitopsflg,nsiz0flg
    real ::  sfdx,sfdy,sfplat,sfplon,sflat,sflon,stoptenh,stoptwvl  &
         ,sz0max,sz0fact,glatr,glonr

    character(len=f_name_length) :: flnm
    character(len=2) :: cgrid
    logical there
    integer :: ios

    character(len=16) :: str(10)
    character(len=*), parameter :: h="**(top_check)**"
    logical, parameter :: dumpLocal=.true.

    lc=len_trim(topfiles)
    write(cgrid,'(a1,i1)') 'g',ifm
    if (useVfm) then
       flnm=trim(topfiles)//'-S-'//cgrid//'.vfm'
    else
       flnm=trim(topfiles)//'-S-'//cgrid//'.bin'
    end if
    
    inquire(file=flnm(1:len_trim(flnm)),exist=there)

    if (dumpLocal) then
       if (there) then
          call MsgDump(h//" found topo file "//trim(flnm))
       else
          call MsgDump(h//" did not find topo file "//trim(flnm))
       end if
    end if

    if(.not.there) then
       ierr = 1
       print*,'TOPfile:', trim(flnm), ' for grid ',ifm,' not there.'
       print*,'------------------------------------------------'
       return
    endif

    call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

    if (useVfm) then
       call rams_f_open(fUnit,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)

       read (fUnit,*) isfc_marker,isfc_ver
       read (fUnit,*) nsfx,nsfy  &
            ,sfdx,sfdy,sfplat,sfplon,sflat,sflon
       read (fUnit,*) nsitoptflg  &
            ,nsitopsflg,stoptenh,stoptwvl,nsiz0flg,sz0max,sz0fact
       close (fUnit)

    else
       open(fUnit, action="read", file=trim(flnm), form="unformatted", iostat=ios)
       if (ios /= 0) then
          call fatal_error(h//" failure opening file "//trim(flnm))
       end if
       read (fUnit) isfc_marker,isfc_ver
       read (fUnit) nsfx,nsfy  &
            ,sfdx,sfdy,sfplat,sfplon,sflat,sflon
       read (fUnit) nsitoptflg  &
            ,nsitopsflg,stoptenh,stoptwvl,nsiz0flg,sz0max,sz0fact
       close (fUnit)
    end if

    if (dumpLocal) then
       write(str(1),"(e16.8)") sflon
       call MsgDump(h//" read sflon="//trim(adjustl(str(1))))
    end if
    
    if (nsfx                       .ne. nnxp(ifm)     .or.  &
         nsfy                       .ne. nnyp(ifm)     .or.  &
         abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
         abs(sfdy-deltayn(ifm))     .gt. .001          .or.  &
         abs(sfplat-platn(ifm))     .gt. .001          .or.  &
         abs(sfplon-plonn(ifm))     .gt. .001          .or.  &
         abs(sflat-glatr)           .gt. .001          .or.  &
         abs(sflon-glonr)           .gt. .001          .or.  &
         nsitoptflg                 .ne. itoptflg(ifm) .or.  &
         nsitopsflg                 .ne. itopsflg(ifm) .or.  &
         abs(stoptenh-toptenh(ifm)) .gt. .001          .or.  &
         abs(stoptwvl-toptwvl(ifm)) .gt. .001          .or.  &
         nsiz0flg                   .ne. iz0flg(ifm)   .or.  &
         abs(sz0max-z0max(ifm))     .gt. .001          .or.  &
         abs(sz0fact-z0fact)        .gt. .00001) then

       ierr = 1

       print*,'SFCfile mismatch on grid:',ifm
       print*,'Values: model, file'
       print*,'-------------------'
       print*,'nnxp:',nnxp(ifm),nsfx
       print*,'nnyp:',nnyp(ifm),nsfy
       print*,'deltax:',deltaxn(ifm),sfdx
       print*,'deltay:',deltayn(ifm),sfdy
       print*,'platn:',platn(ifm),sfplat
       print*,'plonn:',plonn(ifm),sfplon
       print*,'SW lat:',glatr,sflat
       print*,'SW lon:',glonr,sflon
       print*,'itoptflg:',itoptflg(ifm),nsitoptflg
       print*,'itopsflg:',itopsflg(ifm),nsitopsflg
       print*,'toptenh:',toptenh(ifm),stoptenh
       print*,'toptwvl:',toptwvl(ifm),stoptwvl
       print*,'iz0flg:',iz0flg(ifm),nsiz0flg
       print*,'z0max:',z0max(ifm),sz0max
       print*,'z0fact:',z0fact,sz0fact
       print*,'-------------------'

       if (dumpLocal) then
          call MsgDump(h//" sizes mismatch at file "//trim(flnm))
       end if

    else

       ierr = 0
       if (dumpLocal) then
          call MsgDump(h//" sizes match at file "//trim(flnm))
       end if

    endif

  end subroutine top_check

  ! ****************************************************************************

  subroutine toptinit(n2,n3,ifm,topt,topzo)
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: ifm
    real, intent(out) :: topt(:,:)
    real, intent(out) :: topzo(:,:)

    integer :: i,j
    character(len=*), parameter :: h="**(toptinit)**"
    logical, parameter :: dumpLocal=.false.
    
    ! Fill the TOPT array with a default value of 0.  This default is used only
    ! when a standard RAMS topography dataset is not used and when no overrides
    ! to topography heights are defined in subroutine toptinit_user in the
    ! file ruser.f.

    do j = 1,n3
       do i = 1,n2
          topt(i,j) = 0.
          topzo(i,j) = .0001
       enddo
    enddo
    if (dumpLocal) then
       call MsgDump(h//" initializes topt and topzo arrays")
    end if
  end subroutine toptinit


  !*****************************************************************************

  subroutine top_write(ifm)
    integer, intent(in) :: ifm

    integer :: ip,k,i,j
    integer :: ios
    real :: glatr,glonr
    character(len=f_name_length) :: flnm
    character(len=2) :: cgrid
    character(len=*), parameter :: h="**(top_write)**"
    logical, parameter :: dumpLocal=.false.

    !     write surface characteristics, one file for each grid

    write(cgrid,'(a1,i1)') 'g',ifm

    if (useVfm) then
       flnm=trim(topfiles)//'-S-'//cgrid//'.vfm'
    else
       flnm=trim(topfiles)//'-S-'//cgrid//'.bin'
    end if

    if (dumpLocal) then
       call MsgDump(h//" will write file "//trim(flnm))
    end if
    
    call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

    if (useVfm) then
       call rams_f_open (fUnit,flnm(1:len_trim(flnm)),'FORMATTED','REPLACE','WRITE',1)
       rewind fUnit

       write(fUnit,99) 999999,3
99     format(2i8)

       write(fUnit,100) nnxp(ifm),nnyp(ifm)  &
            ,deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr
100    format(2i5,2f15.5,4f11.5)

       write(fUnit,101) itoptflg(ifm),itopsflg(ifm),toptenh(ifm),toptwvl(ifm)  &
            ,iz0flg(ifm),z0max(ifm),z0fact
101    format(2i5,2f11.5,i5,2f11.5)


       call vforec(fUnit,sfcfile_p(ifm)%topt,nnxyp(ifm),24,scrx,'LIN')
       call vforec(fUnit,sfcfile_p(ifm)%topzo,nnxyp(ifm),24,scrx,'LIN')
       
       close(fUnit)

    else
       open(fUnit, action="write", file=trim(flnm), form="unformatted", iostat=ios)
       if (ios /= 0) then
          call fatal_error(h//" failure opening file "//trim(flnm))
       end if
       rewind fUnit
       write(fUnit) 999999,3
       write(fUnit) nnxp(ifm),nnyp(ifm)  &
            ,deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr
       write(fUnit) itoptflg(ifm),itopsflg(ifm),toptenh(ifm),toptwvl(ifm)  &
            ,iz0flg(ifm),z0max(ifm),z0fact

       write(fUnit) sfcfile_p(ifm)%topt
       write(fUnit) sfcfile_p(ifm)%topzo
       
       close(fUnit)
    end if

    if (dumpLocal) then
       call MsgDump(h//" wrote file "//trim(flnm)//&
            " containing topt and topzo")
    end if
    
    return
  end subroutine top_write




  subroutine TopReadStoreOwnChunk(ifm, oneControlVars)
    integer, intent(in) :: ifm
    type(ControlVars), pointer, intent(in) :: oneControlVars

    character(len=f_name_length) :: flnm
    character(len=2) :: cgrid
    character(len=1) :: dummy
    logical :: there
    integer :: ierr
    integer :: ios
    character(len=8) :: c0
    character(len=*), parameter :: h="**(TopReadStoreOwnChunk)**"
    logical, parameter :: dumpLocal=.false.
    
    ! master opens file

    if (mchnum == master_num) then
       write(cgrid,'(a1,i1)') 'g',ifm

       if (useVfm) then
          flnm=trim(topfiles)//'-S-'//cgrid//'.vfm'
       else
          flnm=trim(topfiles)//'-S-'//cgrid//'.bin'
       end if

       inquire(file=flnm(1:len_trim(flnm)),exist=there)

       if (dumpLocal) then
          if (there) then
             call MsgDump(h//" found file "//trim(flnm))
          else
             call MsgDump(h//" did not find file "//trim(flnm))
          end if
       end if
       
       if(.not.there) then
          !call fatal_error(h//" file "//trim(flnm)//" not there")
          iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
               " file "//trim(flnm)//" not there")
       endif

       if (useVfm) then
          call rams_f_open(fUnit,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
          rewind fUnit

          ! Skip header records (already checked)
          
          read(fUnit,*) dummy
          read(fUnit,*) dummy
          read(fUnit,*) dummy
       else
          open(fUnit, action="read", file=trim(flnm), form="unformatted", iostat=ios)
          if (ios /= 0) then
             call fatal_error(h//" failure opening file "//trim(flnm))
          end if
          rewind fUnit
          ! Skip header records (already checked)
          
          read(fUnit) dummy
          read(fUnit) dummy
          read(fUnit) dummy
       end if
          
    end if

    if (dumpLocal) then
       call MsgDump(h//" opened file "//trim(flnm)//" and read preamble")
    end if

    ! master reads twice and broadcast full domain;
    ! local chunk is extracted and stored at desired variable

    call ReadStoreOwnChunk(ifm, fUnit, grid_g(ifm)%topta, "topta", oneControlVars)
    call ReadStoreOwnChunk(ifm, fUnit, grid_g(ifm)%topzo, "topzo", oneControlVars)

    ! master process close file

    if (mchnum == master_num) then
       close (fUnit)
    end if
  end subroutine TopReadStoreOwnChunk





  subroutine TRSFFieldAndOwnChunk(ifm, oneControlVars)
    integer, intent(in) :: ifm
    type(ControlVars), pointer, intent(in) :: oneControlVars

    character(len=f_name_length) :: flnm
    character(len=2) :: cgrid
    character(len=1) :: dummy
    logical :: there
    integer :: ierr
    integer :: ios
    character(len=8) :: c0
    character(len=*), parameter :: h="**(TRSFFieldAndOwnChunk)**"
    logical :: dumpLocal=.false.

    ! master opens file

    if (mchnum == master_num) then

       write(cgrid,'(a1,i1)') 'g',ifm

       if (useVfm) then
          flnm=trim(topfiles)//'-S-'//cgrid//'.vfm'
       else
          flnm=trim(topfiles)//'-S-'//cgrid//'.bin'
       end if

       inquire(file=flnm(1:len_trim(flnm)),exist=there)

       if (dumpLocal) then
          if (there) then
             call MsgDump(h//" found file "//trim(flnm))
          else
             call MsgDump(h//" did not find file "//trim(flnm))
          end if
       end if

       if(.not.there) then
          !call fatal_error(h//" file "//trim(flnm)//" not there")
          iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
               " file "//trim(flnm)//" not there")
       endif

       if (useVfm) then

          call rams_f_open(fUnit,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
          rewind fUnit

          ! Skip header records (already checked)
          
          read(fUnit,*) dummy
          read(fUnit,*) dummy
          read(fUnit,*) dummy

       else

          open(fUnit, action="read", file=trim(flnm), form="unformatted", iostat=ios)
          if (ios /= 0) then
             call fatal_error(h//" failure opening file "//trim(flnm))
          end if
          rewind fUnit
          ! Skip header records (already checked)
          
          read(fUnit) dummy
          read(fUnit) dummy
          read(fUnit) dummy
          
       end if

       if (dumpLocal) then
          call MsgDump(h//" opened file "//trim(flnm)//" and read preamble")
       end if
          
    end if

    ! master reads twice and broadcast full domain;
    ! local chunk is extracted and stored at desired variable

    call ReadStoreFullFieldAndOwnChunk(ifm, fUnit, &
         oneGlobalGridData(ifm)%global_topta, grid_g(ifm)%topta, "topta", &
         oneControlVars)

    call ReadStoreOwnChunk(ifm, fUnit, grid_g(ifm)%topzo, "topzo", &
         oneControlVars)

    ! master process close file

    if (mchnum == master_num) then
       close (fUnit)
    end if

    if (dumpLocal) then
       call MsgDump(h//" read topta and topzo ")
    end if

  end subroutine TRSFFieldAndOwnChunk

end module ModMkSfcTop
