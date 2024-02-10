!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module ModMkSfcSfc

  use ModParallelEnvironment, only: &
       MsgDump

  use mem_grid, only: &
       platn, &
       plonn, &
       xtn, &
       ytn, &
       deltaxn, &
       deltayn, &
       nnxp, &
       nnyp, &
       nzg, &
       npatch, &
       nnxyp

  use io_params, only: &
       sfcfiles, &
       ivegtflg, &
       isoilflg, &
       nofilflg

  use mem_mksfc, only: &
       sfcfile_p, &
       scrx

  use dump, only: &
       dumpMessage

  use mem_leaf, only: &
       leaf_g

  use node_mod, only:      &
       mchnum,             &
       mynum,              &
       master_num, &
       nodemxp, &
       nodemyp

  use ReadBcst, only: &
       ReadStoreOwnChunk

  use ModControlVars, only: &
       ControlVars

  implicit none

  include "files.h"
  include "constants.h"
  include "UseVfm.h"
  integer, parameter :: fUnit=25

  private
  public :: sfc_check
  public :: sfc_write
  public :: SfcReadStoreOwnChunk

contains




  subroutine sfc_check(ifm,ierr)

    ! This subroutine checks for the existence of a surface file for
    ! grid number ifm, and if it exists, also checks for agreement of
    ! grid configuration between the file and the current model run.
    ! If the file does not exist or does not match grid configuration,
    ! the flag ifileok is returned with a value of 0.  If the file
    ! exists and is ok, ifileok is returned with a value of 1.

    integer, intent(in) :: ifm
    integer, intent(out) :: ierr
    
    integer :: lc,isfc_marker,isfc_ver,nsfx,nsfy,nsfzg  &
         ,nsivegtflg,nsisoilflg,nsnofilflg,nspatch, ios
    real ::  sfdx,sfdy,sfplat,sfplon,sflat,sflon,glatr,glonr


    character(len=f_name_length) :: flnm
    character(len=2) :: cgrid
    logical there

    character(len=16) :: str(10)
    character(len=*), parameter :: h="**(sfc_check)**"
    logical, parameter :: dumpLocal=.false.

    flnm=""
    lc=len_trim(sfcfiles)
    write(cgrid,'(a1,i1)') 'g',ifm
    if (useVfm) then
       flnm=trim(sfcfiles)//'-S-'//cgrid//'.vfm'
    else
       flnm=trim(sfcfiles)//'-S-'//cgrid//'.bin'
    end if
    
    inquire(file=flnm(1:len_trim(flnm)),exist=there)

    if (dumpLocal) then
       if (there) then
          call MsgDump(h//" found file "//trim(adjustl(flnm))) 
       else
          call MsgDump(h//" did not find file "//trim(adjustl(flnm))) 
       end if
    end if

    if(.not.there) then
       ierr = 1
       print*,'SFCfile:', trim(flnm), ' for grid ',ifm,' not there.'
       print*,'------------------------------------------------'
       return
    endif

    call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

    if (useVfm) then

       call rams_f_open(fUnit,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
       read (fUnit,*) isfc_marker,isfc_ver
       read (fUnit,100) nsfx,nsfy,nsfzg,nspatch  &
            ,sfdx,sfdy,sfplat,sfplon,sflat,sflon
       read (fUnit,101) nsivegtflg,nsisoilflg,nsnofilflg
100    format(4i5,2f15.5,4f11.5)
101    format(5i5,2f11.5,i5,2f11.5)
       close (fUnit)
    else
       
       open(fUnit, action="read", file=trim(flnm), form="unformatted", iostat=ios)
       read (fUnit) isfc_marker,isfc_ver
       read (fUnit) nsfx,nsfy,nsfzg,nspatch  &
            ,sfdx,sfdy,sfplat,sfplon,sflat,sflon
       read (fUnit) nsivegtflg,nsisoilflg,nsnofilflg
       close (fUnit)
    end if

    if (dumpLocal) then
       write(str(1),"(e16.8)") sflon
       call MsgDump(h//" read sflon="//trim(adjustl(str(1))))
    end if

    if (nsfx                       .ne. nnxp(ifm)     .or.  &
         nsfy                       .ne. nnyp(ifm)     .or.  &
         nsfzg                      .ne. nzg           .or.  &
         nspatch                    .ne. npatch        .or.  &
         abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
         abs(sfdy-deltayn(ifm))     .gt. .001          .or.  &
         abs(sfplat-platn(ifm))     .gt. .001          .or.  &
         abs(sfplon-plonn(ifm))     .gt. .001          .or.  &
         abs(sflat-glatr)           .gt. .001          .or.  &
         abs(sflon-glonr)           .gt. .001          .or.  &
         nsivegtflg                 .ne. ivegtflg(ifm) .or.  &
         nsisoilflg                 .ne. isoilflg(ifm) .or.  &
         nsnofilflg                 .ne. nofilflg(ifm) ) then

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
       print*,'ivegtflg:',ivegtflg(ifm),nsivegtflg
       print*,'isoilflg:',isoilflg(ifm),nsisoilflg
       print*,'nofilflg:',nofilflg(ifm),nsnofilflg
       print*,'-------------------'

    else

       ierr = 0

    endif

    if (dumpLocal) then
       if (ierr == 0) then
          call MsgDump(h//" dimensions agree at file "//trim(adjustl(flnm)))
       else
          call MsgDump(h//" dimensions disagree at file "//trim(adjustl(flnm)))
       end if
    end if
    return
  end subroutine sfc_check

  !*****************************************************************************

  subroutine sfc_write(ifm)
    integer, intent(in) :: ifm

    integer :: ip,k,i,j, ios
!!$    integer :: ii
    real :: glatr,glonr
    character(len=f_name_length) :: flnm
    character(len=2) :: cgrid
    character(len=*), parameter :: h="**(sfc_write)**"
    character(len=8) :: str(10)
    character(len=16) :: lstr
    character(len=256) :: longStr
    logical, parameter :: dumpLocal=.false.

    !     write surface characteristics, one file for each grid

    flnm=""
    write(cgrid,'(a1,i1)') 'g',ifm
    
    if (useVfm) then
       flnm=trim(sfcfiles)//'-S-'//cgrid//'.vfm'
    else
       flnm=trim(sfcfiles)//'-S-'//cgrid//'.bin'
    end if
    
    if (dumpLocal) then
       call MsgDump(h//" will write file "//trim(adjustl(flnm)))
    end if

    call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

    print *, h//" glonr apos xy_ll=",glonr

    if (useVfm) then
       call rams_f_open (fUnit,flnm(1:len_trim(flnm)),'FORMATTED','REPLACE','WRITE',1)
       rewind fUnit
       
       write(fUnit,99) 999999,3
99     format(2i8)

       write(fUnit,100) nnxp(ifm),nnyp(ifm),nzg,npatch  &
            ,deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr
100    format(4i5,2f15.5,4f11.5)

       write(fUnit,101) ivegtflg(ifm),isoilflg(ifm),nofilflg(ifm)
101    format(3i5)


       do ip = 1,npatch
          call vforec(fUnit,sfcfile_p(ifm)%patch_area(1,1,ip),nnxyp(ifm),24,scrx,'LIN')
       enddo

       do ip = 1,npatch
          call vforec(fUnit,sfcfile_p(ifm)%leaf_class(1,1,ip),nnxyp(ifm),24,scrx,'LIN')
       enddo
          
       do ip = 1,npatch
          call vforec(fUnit,sfcfile_p(ifm)%soil_text(1,1,1,ip),nzg*nnxyp(ifm),24,scrx,'LIN')
       enddo

       close(fUnit)

    else

       open(fUnit, action="write", file=trim(flnm), form="unformatted", iostat=ios)
       if (ios /= 0) then
          call fatal_error(h//" opening file "//trim(flnm))
       end if
       
       write(fUnit) 999999,3

       write(fUnit) nnxp(ifm),nnyp(ifm),nzg,npatch  &
            ,deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr

       write(fUnit) ivegtflg(ifm),isoilflg(ifm),nofilflg(ifm)

       do ip = 1,npatch
          write(fUnit) sfcfile_p(ifm)%patch_area(:,:,ip)
       enddo

       do ip = 1,npatch
          write(fUnit) sfcfile_p(ifm)%leaf_class(:,:,ip)
       enddo
       
       do ip = 1,npatch
          write(fUnit) sfcfile_p(ifm)%soil_text(:,:,:,ip)
       enddo

       close(fUnit)
    end if

    if (dumpLocal) then
       call MsgDump(h//" wrote file "//trim(adjustl(flnm))//&
            " containing patch_area, leaf_class and soil_text")
    end if

    print *, h//" escreveu sflon (glonr)= ",glonr

  end subroutine sfc_write



  subroutine SfcReadStoreOwnChunk(ifm, oneControlVars)
    integer, intent(in) :: ifm
    type(ControlVars), pointer, intent(in) :: oneControlVars

    integer :: ipat, ios
    logical :: there
    character(len=f_name_length) :: flnm
    character(len=2) :: cgrid
    character(len=1) :: dummy
    character(len=2) :: cipat

    real, pointer :: p2D(:,:), p3D(:,:,:)

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(SfcReadStoreOwnChunk)**"
    logical, parameter :: dumpLocal=.false.

    ! master process opens file and skips headers
    flnm=""
    if (mchnum == master_num) then
       write(cgrid,'(a1,i1)') 'g',ifm
       if (useVfm) then
          flnm=trim(sfcfiles)//'-S-'//cgrid//'.vfm'
       else
          flnm=trim(sfcfiles)//'-S-'//cgrid//'.bin'
       end if
       
       inquire(file=flnm(1:len_trim(flnm)),exist=there)

       if (dumpLocal) then
          if (there) then
             call MsgDump(h//" found file "//flnm(1:len_trim(flnm)))
          else
             call MsgDump(h//" did not find file "//flnm(1:len_trim(flnm)))
          end if
       end if

       if(.not.there) then
          !call fatal_error(h//" file "//trim(flnm)//" not there")
          iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
               " file "//trim(flnm)//" not there")
       end if

       if (useVfm) then

          ! use VFM format
          
          call rams_f_open(fUnit,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
          rewind fUnit
          
          ! Skip header records (already been checked)
          read(fUnit,*) dummy
          read(fUnit,*) dummy
          read(fUnit,*) dummy
          
       else
          
          ! use binary format
          
          open(fUnit, action="read", file=trim(flnm), form="unformatted", iostat=ios)
          read(fUnit) dummy
          read(fUnit) dummy
          read(fUnit) dummy
       end if
    end if

    ! deals with patch_area
       
    do ipat = 1,npatch
       write(cipat,fmt='(I2.2)') ipat
       p2D => leaf_g(ifm)%patch_area(:,:,ipat)
       call ReadStoreOwnChunk(ifm, fUnit, p2D, "patch_area"//cipat, &
            oneControlVars)
    end do

    ! deals with leaf_class

    do ipat = 1,npatch
       p2D => leaf_g(ifm)%leaf_class(:,:,ipat)
       call ReadStoreOwnChunk(ifm, fUnit, p2D, "leaf_class", &
            oneControlVars)
    end do

    ! deals with soil_text

    do ipat = 1,npatch
       p3D => leaf_g(ifm)%soil_text(:,:,:,ipat)
       call ReadStoreOwnChunk(ifm, fUnit, p3D, nzg, "soil_text", &
            oneControlVars)
    end do


    if (mchnum == master_num) then
       close (fUnit)
    end if

  end subroutine SfcReadStoreOwnChunk
end module ModMkSfcSfc
