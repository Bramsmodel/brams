!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModRcio

  use iso_fortran_env, only: &
       real64

  use ModParallelEnvironment, only: &
       MsgDump

  use io_params, only: &
       ioutput, &
       iinput, &
       ntopsmth, &
       izflat

  use grid_dims, only: &
       nzgmax

  use mem_cuparm, only: &
       nnqparm, &
       wcldbs, &
       confrq

  use mem_grid, only: &
       iversion, &
       expnme, &
       if_adap, &
       ngrids, &
       nzg, &
       nzs, &
       naddsc, &
       time, &
       ztop, &
       polelat, &
       polelon, &
       dzrat, &
       dzmax, &
       nnxp, &
       nnyp, &
       nnzp, &
       nndtrat, &
       nstratx, &
       nstraty, &
       ngbegun, &
       nnacoust, &
       nxtnest, &
       nnsttop, &
       nnstbot, &
       ninest, &
       njnest, &
       nknest, &
       gridu, &
       gridv, &
       centlat, &
       centlon, &
       dimove, &
       djmove, &
       deltazn, &
       deltaxn, &
       deltayn, &
       platn, &
       plonn, &
       zz, &
       nestz1, &
       nestz2, &
       nstratz1, &
       nstratz2, &
       xmn, &
       xtn, &
       ymn, &
       ytn, &
       zmn, &
       ztn, &
       dzmn, &
       dztn, &
       itopo, &
       initial, &
       impl, &
       jdim, &
       iadvl, &
       iadvf, &
       lsflg, &
       ibnd, &
       jbnd, &
       icorflg, &
       vveldamp, &
       ihtran, &
       nfpt, &
       ideltat, &
       nacoust, &
       iflag, &
       iyear1, &
       imonth1, &
       idate1, &
       itime1, &
       npatch, &
       distim, &
       eps, &
       cphas, &
       sspct

  use mem_leaf, only: &
       isfcl, &
       nvegpat, &
       drtcon, &
       seatmp, &
       albedo, &
       dthcon, &
       slz

  use mem_radiate, only: &
       lonrad, &
       ilwrtyp, &
       iswrtyp, &
       radfrq

  use ModMicControl, only: &
       MicControl

  use ref_sounding, only: &
       maxsndg, &
       u01dn, &
       v01dn, &
       pi01dn, &
       th01dn, &
       dn01dn, &
       rt01dn, &
       iref, &
       jref, &
       nsndg, &
       topref, &
       us, &
       vs, &
       ts, &
       thds, &
       ps, &
       hs

  use ModLeafComs, only: &
       nstyp, &
       nvtyp, &
       nvtyp_teb, &
       kroot, &
       slden, &
       slcpd, &
       slbs, &
       slcond, &
       slcons, &
       slmsts, &
       slpots, &
       ssand, &
       sclay, &
       sorgan, &
       sporo, &
       soilcp, &
       slfc, &
       emisg, &
       emisv, &
       root, &
       cmin, &
       corg, &
       cwat, &
       cair, &
       cka, &
       ckw


  use mem_stilt, only: &
       iexev,          &
       imassflx

  use ModNamelistFile, only: &
       NamelistFile

  use an_header, only: &
       IOFileDS

  implicit none

  include "MicConstants.h"
  
  private

  public :: DumpIOHeadTable
  public :: cio

  interface cio
     module procedure cio_i_s
     module procedure cio_i_1d
     module procedure cio_c_s
     module procedure cio_f_s
     module procedure cio_f_1d
     module procedure cio_f_2d
     module procedure cio_f8_s
  end interface cio
contains


  integer function cio_i_s(iun,irw,cstr,ia)
    integer, intent(in) :: iun
    integer, intent(in) :: irw
    character(len=*), intent(in) :: cstr
    integer :: ia

    integer :: nn
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(cio_i_s)**"

    if (irw == 1) then
       call cio_pos_file (iun,cstr,cio_i_s)
       if (cio_i_s == 1) then
          write(str(1),"(i8)") iun
          call fatal_error(h//" string "//cstr//&
               " not found at unit "//trim(adjustl(str(1))))
       end if
       read(iun,*) nn
       read(iun,*) ia
    else if (irw == 2) then
       write(iun,"('__',a)") cstr
       write(iun,*) 1
       write(iun,"(i6)") ia
       cio_i_s=0
    else
       write(str(1),"(i8)") irw
       call fatal_error(h//" invoked with unknown irw="//&
            trim(adjustl(str(1))))
    end if
  end function cio_i_s




  integer function cio_i_1d(iun,irw,cstr,ia)
    integer, intent(in) :: iun
    integer, intent(in) :: irw
    character(len=*), intent(in) :: cstr
    integer :: ia(:)

    integer :: nn
    integer :: i
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(cio_i_1d)**"

    if (irw == 1) then
       call cio_pos_file (iun,cstr,cio_i_1d)
       if (cio_i_1d == 1) then
          write(str(1),"(i8)") iun
          call fatal_error(h//" string "//cstr//&
               " not found at unit "//trim(adjustl(str(1))))
       end if
       read(iun,*) nn
       read(iun,*) ia
    else if (irw == 2) then
       nn=size(ia)
       write(iun,"('__',a)") cstr
       write(iun,*) nn
       write(iun,"(i6)") (ia(i),i=1,nn)
       cio_i_1d=0
    else
       write(str(1),"(i8)") irw
       call fatal_error(h//" invoked with unknown irw="//&
            trim(adjustl(str(1))))
    end if
  end function cio_i_1d



  integer function cio_c_s(iun,irw,cstr,ia)
    integer, intent(in) :: iun
    integer, intent(in) :: irw
    character(len=*), intent(in) :: cstr
    character(len=*) :: ia

    integer :: nn
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(cio_c_s)**"

    if (irw == 1) then
       call cio_pos_file (iun,cstr,cio_c_s)
       if (cio_c_s == 1) then
          write(str(1),"(i8)") iun
          call fatal_error(h//" string "//cstr//&
               " not found at unit "//trim(adjustl(str(1))))
       end if
       read(iun,*) nn
       read(iun,"(a)") ia
    else if (irw == 2) then
       write(iun,"('__',a)") cstr
       write(iun,*) 1
       write(iun,"(a)") ia
       cio_c_s=0
    else
       write(str(1),"(i8)") irw
       call fatal_error(h//" invoked with unknown irw="//&
            trim(adjustl(str(1))))
    endif
  end function cio_c_s



  integer function cio_f_s(iun,irw,cstr,ia)
    integer, intent(in) :: iun
    integer, intent(in) :: irw
    character(len=*), intent(in) :: cstr
    real :: ia

    integer :: nn
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(cio_f_s)**"

    if (irw == 1) then
       call cio_pos_file (iun,cstr,cio_f_s)
       if (cio_f_s == 1) then
          write(str(1),"(i8)") iun
          call fatal_error(h//" string "//cstr//&
               " not found at unit "//trim(adjustl(str(1))))
       end if
       read(iun,*) nn
       read(iun,*) ia
    else if (irw == 2) then
       write(iun,"('__',a)") cstr
       write(iun,*) 1
       write(iun,"(e16.8)") ia
       cio_f_s=0
    else
       write(str(1),"(i8)") irw
       call fatal_error(h//" invoked with unknown irw="//&
            trim(adjustl(str(1))))
    endif
  end function cio_f_s




  integer function cio_f_1d(iun,irw,cstr,ia)
    integer, intent(in) :: iun
    integer, intent(in) :: irw
    character(len=*), intent(in) :: cstr
    real :: ia(:)

    integer :: nn
    integer :: i
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(cio_f_1d)**"

    nn=size(ia)
    if (irw == 1) then
       call cio_pos_file (iun,cstr,cio_f_1d)
       if (cio_f_1d == 1) then
          write(str(1),"(i8)") iun
          call fatal_error(h//" string "//cstr//&
               " not found at unit "//trim(adjustl(str(1))))
       end if
       read(iun,*) nn
       read(iun,*) (ia(i),i=1,nn)
    else if (irw == 2) then
       write(iun,"('__',a)") cstr
       write(iun,*) nn
       write(iun,"(e16.8)") (ia(i),i=1,nn)
       cio_f_1d=0
    else
       write(str(1),"(i8)") irw
       call fatal_error(h//" invoked with unknown irw="//&
            trim(adjustl(str(1))))
    endif
  end function cio_f_1d




  integer function cio_f_2d(iun,irw,cstr,ia)
    integer, intent(in) :: iun
    integer, intent(in) :: irw
    character(len=*), intent(in) :: cstr
    real :: ia(:,:)

    integer :: nn
    integer :: n1
    integer :: n2
    integer :: i
    integer :: j
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(cio_f_2d)**"

    nn=size(ia)
    n1=size(ia,1)
    n2=size(ia,2)
    if (irw == 1) then
       call cio_pos_file (iun,cstr,cio_f_2d)
       if (cio_f_2d == 1) then
          write(str(1),"(i8)") iun
          call fatal_error(h//" string "//cstr//&
               " not found at unit "//trim(adjustl(str(1))))
       end if
       read(iun,*) nn
       read(iun,*) ((ia(i,j),i=1,n1),j=1,n2)
    else if (irw == 2) then
       write(iun,"('__',a)") cstr
       write(iun,*) nn
       write(iun,"(e16.8)") ((ia(i,j),i=1,n1),j=1,n2)
       cio_f_2d=0
    else
       write(str(1),"(i8)") irw
       call fatal_error(h//" invoked with unknown irw="//&
            trim(adjustl(str(1))))
    endif
  end function cio_f_2d



  integer function cio_f8_s(iun,irw,cstr,ia)
    integer, intent(in) :: iun
    integer, intent(in) :: irw
    character(len=*), intent(in) :: cstr
    real(kind=real64) :: ia

    integer :: nn
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(cio_f8_s)**"

    if (irw == 1) then
       call cio_pos_file (iun,cstr,cio_f8_s)
       if (cio_f8_s == 1) then
          write(str(1),"(i8)") iun
          call fatal_error(h//" string "//cstr//&
               " not found at unit "//trim(adjustl(str(1))))
       end if
       read(iun,*) nn
       read(iun,*) ia
    else if (irw == 2) then
       write(iun,"('__',a)") cstr
       write(iun,*) 1
       write(iun,"(e24.16)") ia
       cio_f8_s=0
    else
       write(str(1),"(i8)") irw
       call fatal_error(h//" invoked with unknown irw="//&
            trim(adjustl(str(1))))
    endif
  end function cio_f8_s




  subroutine DumpIOHeadTable(oneIOFileDS, oneNamelistFile, oneMicControl)
    type(IOFileDS), intent(inout) :: oneIOFileDS
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(MicControl), pointer, intent(in) :: oneMicControl

    integer, parameter :: unitLow=10
    integer, parameter :: unitHigh=99
    integer :: iunit
    integer :: nv
    logical :: op
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DumpIOHeadTable)**"
    logical, parameter :: dumpLocal=.false.

    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    ! find available Fortran unit

    do iunit = unitLow, unitHigh
       inquire (unit=iunit, opened=op)
       if (.not. op) exit
    end do
    if (iunit <= unitHigh) then
       oneIOFileDS%unit = iunit
    else
       call fatal_error(h//" Fortran i/o units exausted")
    end if

    ! open file, write header information, close file

    call rams_f_open_u(oneIOFileDS%unit, oneIOFileDS%fHeadName,  &
         'FORMATTED','REPLACE','WRITE', "ASIS", oneIOFileDS%iclobber)
    if (dumpLocal) then
       write(c0,"(i8)") oneIOFileDS%unit
       write(*, "(a)") h//" opened new file "//trim(oneIOFileDS%fHeadName)//&
            " at unit "//trim(adjustl(c0))
    end if

    write(oneIOFileDS%unit,'(i6)') oneIOFileDS%ht%lastUsed
    do nv = 1, oneIOFileDS%ht%lastUsed
       write(oneIOFileDS%unit,fmt='(a16,1x,i12,i3,i3,1x,i9)') &
            oneIOFileDS%ht%f(nv)%string,                      &
            oneIOFileDS%ht%f(nv)%npointer,                    &
            oneIOFileDS%ht%f(nv)%idim_type,                   &
            oneIOFileDS%ht%f(nv)%ngrid,                       &
            oneIOFileDS%ht%f(nv)%nvalues
    end do

    call commio(oneIOFileDS%fId,'WRITE',oneIOFileDS%unit, oneNamelistFile, oneMicControl)
    close(oneIOFileDS%unit)
    if (dumpLocal) then
       write(c0,"(i8)") oneIOFileDS%unit
       write(*, "(a)") h//" dumped and close file "//trim(oneIOFileDS%fHeadName)//&
            " at unit "//trim(adjustl(c0))
    end if

    oneIOFileDS%unit = -1
  end subroutine DumpIOHeadTable




  subroutine commio (cfile,io,iun,oneNamelistFile,oneMicControl)
    character(len=*) :: cfile !**(JP)** not used
    character(len=*), intent(in) :: io
    integer, intent(in) :: iun
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(MicControl), pointer, intent(in) :: oneMicControl


    !     This routine reads or writes the history and analysis file common blocks.

    character(len=2) :: cng
    integer :: irw,ie,ng

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(commio)**"
    logical, parameter :: dumpLocal=.false.

    !**(JP)**
    ! variables originaly declared at mem_turb but not referenced elsewhere
    ! still dumped at XXX-head.txt to maintain reproducibility
    ! should be removed on future versions
    real, parameter :: brunt=0.0
    real, parameter :: rmin=0.0
    real, parameter :: rmax=0.0

    if (io == 'READ') then
       irw=1
    else if (io == 'WRITE') then
       irw=2
    else
       call fatal_error(h//" invoked with io not READ or WRITE, but **"//&
            io//"**")
    end if

    ie=cio(iun,irw,'iversion',iversion)
    ie=cio(iun,irw,'expnme',expnme)
    ie=cio(iun,irw,'ioutput',ioutput)
    ie=cio(iun,irw,'if_adap',if_adap)
    ie=cio(iun,irw,'ngrids',ngrids)
    ie=cio(iun,irw,'nzg',nzg)
    ie=cio(iun,irw,'nzs',nzs)
    ie=cio(iun,irw,'naddsc',naddsc)
    ie=cio(iun,irw,'time',time)
    ie=cio(iun,irw,'ztop',ztop)
    ie=cio(iun,irw,'polelat',polelat)
    ie=cio(iun,irw,'polelon',polelon)
    ie=cio(iun,irw,'dzrat',dzrat)
    ie=cio(iun,irw,'dzmax',dzmax)

    ie=cio(iun,irw,'nnxp',nnxp(1:ngrids))
    ie=cio(iun,irw,'nnyp',nnyp(1:ngrids))
    ie=cio(iun,irw,'nnzp',nnzp(1:ngrids))
    ie=cio(iun,irw,'nnqparm',nnqparm(1:ngrids))
    ie=cio(iun,irw,'nndtrat',nndtrat(1:ngrids))
    ie=cio(iun,irw,'nstratx',nstratx(1:ngrids))
    ie=cio(iun,irw,'nstraty',nstraty(1:ngrids))
    ie=cio(iun,irw,'ngbegun',ngbegun(1:ngrids))
    ie=cio(iun,irw,'nnacoust',nnacoust(1:ngrids))
    ie=cio(iun,irw,'nxtnest',nxtnest(1:ngrids))
    ie=cio(iun,irw,'nnsttop',nnsttop(1:ngrids))
    ie=cio(iun,irw,'nnstbot',nnstbot(1:ngrids))
    ie=cio(iun,irw,'ninest',ninest(1:ngrids))
    ie=cio(iun,irw,'njnest',njnest(1:ngrids))
    ie=cio(iun,irw,'nknest',nknest(1:ngrids))
    ie=cio(iun,irw,'idiffk',oneNamelistFile%idiffk(1:ngrids))
    ie=cio(iun,irw,'gridu',gridu(1:ngrids))
    ie=cio(iun,irw,'gridv',gridv(1:ngrids))
    ie=cio(iun,irw,'akmin',oneNamelistFile%akmin(1:ngrids))
    ie=cio(iun,irw,'csz',oneNamelistFile%csz(1:ngrids))
    ie=cio(iun,irw,'csx',oneNamelistFile%csx(1:ngrids))
    ie=cio(iun,irw,'xkhkm',oneNamelistFile%xkhkm(1:ngrids))
    ie=cio(iun,irw,'zkhkm',oneNamelistFile%zkhkm(1:ngrids))
    ie=cio(iun,irw,'centlat',centlat(1:ngrids))
    ie=cio(iun,irw,'centlon',centlon(1:ngrids))
    ie=cio(iun,irw,'dimove',dimove(1:ngrids))
    ie=cio(iun,irw,'djmove',djmove(1:ngrids))
    ie=cio(iun,irw,'deltazn',deltazn(1:ngrids))
    ie=cio(iun,irw,'deltaxn',deltaxn(1:ngrids))
    ie=cio(iun,irw,'deltayn',deltayn(1:ngrids))
    ie=cio(iun,irw,'platn',platn(1:ngrids))
    ie=cio(iun,irw,'plonn',plonn(1:ngrids))
    ie=cio(iun,irw,'zz',zz(1:nnzp(1)))

    ie=cio(iun,irw,'nestz1',nestz1)
    ie=cio(iun,irw,'nestz2',nestz2)
    ie=cio(iun,irw,'nstratz1',nstratz1(1:nnzp(1)))
    ie=cio(iun,irw,'nstratz2',nstratz2(1:nnzp(1)))

    do ng=1,ngrids
       write(cng,1) ng
1      format(i2.2)
       ie=cio(iun,irw,'xmn'//cng,xmn(1:nnxp(ng),ng))
       ie=cio(iun,irw,'xtn'//cng,xtn(1:nnxp(ng),ng))
       ie=cio(iun,irw,'ymn'//cng,ymn(1:nnyp(ng),ng))
       ie=cio(iun,irw,'ytn'//cng,ytn(1:nnyp(ng),ng))
       ie=cio(iun,irw,'zmn'//cng,zmn(1:nnzp(ng),ng))
       ie=cio(iun,irw,'ztn'//cng,ztn(1:nnzp(ng),ng))
       ie=cio(iun,irw,'dzmn'//cng,dzmn(1:nnzp(ng),ng))
       ie=cio(iun,irw,'dztn'//cng,dztn(1:nnzp(ng),ng))
       ie=cio(iun,irw,'u01dn'//cng,u01dn(1:nnzp(ng),ng))
       ie=cio(iun,irw,'v01dn'//cng,v01dn(1:nnzp(ng),ng))
       ie=cio(iun,irw,'pi01dn'//cng,pi01dn(1:nnzp(ng),ng))
       ie=cio(iun,irw,'th01dn'//cng,th01dn(1:nnzp(ng),ng))
       ie=cio(iun,irw,'dn01dn'//cng,dn01dn(1:nnzp(ng),ng))
       ie=cio(iun,irw,'rt01dn'//cng,rt01dn(1:nnzp(ng),ng))
    enddo

    ie=cio(iun,irw,'kroot',kroot(1:nvtyp+nvtyp_teb))

    ie=cio(iun,irw,'itopo',itopo)
    ie=cio(iun,irw,'initial',initial)
    ie=cio(iun,irw,'impl',impl)
    ie=cio(iun,irw,'iinput',iinput)
    ie=cio(iun,irw,'jdim',jdim)
    ie=cio(iun,irw,'iadvl',iadvl)
    ie=cio(iun,irw,'iadvf',iadvf)
    ie=cio(iun,irw,'lonrad',lonrad)
    ie=cio(iun,irw,'lsflg',lsflg)
    ie=cio(iun,irw,'ibnd',ibnd)
    ie=cio(iun,irw,'jbnd',jbnd)
    ie=cio(iun,irw,'icorflg',icorflg)
    ie=cio(iun,irw,'vveldamp',vveldamp)
    !--(DMK-CCATT)---------------------------------------------------------
    ![ML
    ie=cio(iun,irw,'iexev',iexev)
    ie=cio(iun,irw,'imassflx',imassflx)
    !ML]
    !--(DMK-CCATT-END)-----------------------------------------------------
    ie=cio(iun,irw,'ilwrtyp',ilwrtyp)
    ie=cio(iun,irw,'iswrtyp',iswrtyp)
    ie=cio(iun,irw,'iref',iref)
    ie=cio(iun,irw,'jref',jref)
    ie=cio(iun,irw,'ihtran',ihtran)
    ie=cio(iun,irw,'nfpt',nfpt)
    ie=cio(iun,irw,'nsndg',nsndg)
    ie=cio(iun,irw,'ideltat',ideltat)
    ie=cio(iun,irw,'nacoust',nacoust)
    ie=cio(iun,irw,'iflag',iflag)
    ie=cio(iun,irw,'ntopsmth',ntopsmth)
    ie=cio(iun,irw,'izflat',izflat)
    ie=cio(iun,irw,'iyear1',iyear1)
    ie=cio(iun,irw,'imonth1',imonth1)
    ie=cio(iun,irw,'idate1',idate1)
    ie=cio(iun,irw,'itime1',itime1)
    ie=cio(iun,irw,'isfcl',isfcl)
    ie=cio(iun,irw,'npatch',npatch)
    ie=cio(iun,irw,'nvegpat',nvegpat)
    ie=cio(iun,irw,'mcphys_type',oneMicControl%mcphys_type)
    ie=cio(iun,irw,'level',oneMicControl%level)
    ie=cio(iun,irw,'irain',oneMicControl%irain)
    ie=cio(iun,irw,'ipris',oneMicControl%ipris)
    ie=cio(iun,irw,'isnow',oneMicControl%isnow)
    ie=cio(iun,irw,'iaggr',oneMicControl%iaggr)
    ie=cio(iun,irw,'igraup',oneMicControl%igraup)
    ie=cio(iun,irw,'icloud',oneMicControl%icloud)
    ie=cio(iun,irw,'ihail',oneMicControl%ihail)

    ie=cio(iun,irw,'brunt',brunt)
    ie=cio(iun,irw,'wcldbs',wcldbs)
    ie=cio(iun,irw,'drtcon',drtcon)
    ie=cio(iun,irw,'rmin',rmin)
    ie=cio(iun,irw,'radfrq',radfrq)
    ie=cio(iun,irw,'distim',distim)
    ie=cio(iun,irw,'seatmp',seatmp)
    ie=cio(iun,irw,'confrq',confrq)
    ie=cio(iun,irw,'rmax',rmax)
    ie=cio(iun,irw,'eps',eps)
    ie=cio(iun,irw,'albedo',albedo)
    ie=cio(iun,irw,'dthcon',dthcon)
    ie=cio(iun,irw,'cphas',cphas)
    ie=cio(iun,irw,'topref',topref)
    ie=cio(iun,irw,'sspct',sspct)
    ie=cio(iun,irw,'rparm',oneMicControl%rparm)
    ie=cio(iun,irw,'pparm',oneMicControl%pparm)
    ie=cio(iun,irw,'sparm',oneMicControl%sparm)
    ie=cio(iun,irw,'aparm',oneMicControl%aparm)
    ie=cio(iun,irw,'gparm',oneMicControl%gparm)
    ie=cio(iun,irw,'cparm',oneMicControl%cparm)
    ie=cio(iun,irw,'hparm',oneMicControl%hparm)

    ie=cio(iun,irw,'cfmas',oneMicControl%cfmas(1:nhcat))
    ie=cio(iun,irw,'pwmas',oneMicControl%pwmas(1:nhcat))

    ie=cio(iun,irw,'us',us(1:maxsndg))
    ie=cio(iun,irw,'vs',vs(1:maxsndg))
    ie=cio(iun,irw,'ts',ts(1:maxsndg))
    ie=cio(iun,irw,'thds',thds(1:maxsndg))
    ie=cio(iun,irw,'ps',ps(1:maxsndg))
    ie=cio(iun,irw,'hs',hs(1:maxsndg))

    ie=cio(iun,irw,'slden',slden(1:nstyp))
    ie=cio(iun,irw,'slcpd',slcpd(1:nstyp))
    ie=cio(iun,irw,'slbs',slbs(1:nstyp))
    ie=cio(iun,irw,'slcond',slcond(1:nstyp))
    ie=cio(iun,irw,'slcons',slcons(1:nstyp))
    ie=cio(iun,irw,'slmsts',slmsts(1:nstyp))
    ie=cio(iun,irw,'slpots',slpots(1:nstyp))
    ie=cio(iun,irw,'ssand',ssand(1:nstyp))
    ie=cio(iun,irw,'sclay',sclay(1:nstyp))
    ie=cio(iun,irw,'sorgan',sorgan(1:nstyp))
    ie=cio(iun,irw,'sporo',sporo(1:nstyp))
    ie=cio(iun,irw,'soilcp',soilcp(1:nstyp))
    ie=cio(iun,irw,'slfc',slfc(1:nstyp))
    ie=cio(iun,irw,'emisg',emisg(1:nstyp))

    ie=cio(iun,irw,'emisv',emisv(1:nvtyp+nvtyp_teb))

    ie=cio(iun,irw,'root',root)
    ie=cio(iun,irw,'slz',slz(1:nzg))

    ie=cio(iun,irw,'cmin',cmin)
    ie=cio(iun,irw,'corg',corg)
    ie=cio(iun,irw,'cwat',cwat)
    ie=cio(iun,irw,'cair',cair)
    ie=cio(iun,irw,'cka',cka)
    ie=cio(iun,irw,'ckw',ckw)

    return
  end subroutine commio

  !---------------------------------------------------------

  subroutine cio_pos_file(iun,cstr,ierr)
    implicit none
    integer :: iun,ierr
    character(len=*) :: cstr
    character(len=256) :: line,csearch

    integer :: nl,nc,iend,ilen

    !      print*,'cio_pos:',iun,cstr

    iend=0
1   continue
    do nl=1,1000000
       read(iun,10,end=100) line
10     format(a)
       ilen=len(cstr)
       csearch='__'//cstr(1:ilen)
       nc=index(line,csearch(1:ilen+2) )
       !         print*,'cio_pos:',nl,nc,line
       if(nc.eq.1) then
          ierr=0
          !            print*,'---- Name found on header file:',cstr
          return
       endif
    enddo

100 continue
    if(iend.eq.1) then
       ierr=1
       print*,'---- Name NOT found on header file:',cstr
       rewind(iun)
       return
    endif
    rewind(iun)
    iend=1
    goto 1

  end subroutine cio_pos_file

  !---------------------------------------------------------

!!$integer function cio_i(iun,irw,cstr,ia,n)
!!$  implicit none
!!$  integer :: iun,irw,n
!!$  integer, dimension(n) :: ia
!!$  character(len=*) :: cstr
!!$  character(len=256) :: string
!!$  integer :: nn,i
!!$
!!$  if (irw.eq.1) then
!!$     call cio_pos_file (iun,cstr,cio_i)
!!$     if(cio_i.eq.1) return
!!$     read(iun,*) nn
!!$     read(iun,*) (ia(i),i=1,nn)
!!$  elseif(irw.eq.2) then
!!$     write(iun,20) cstr
!!$20   format('__',a)
!!$     write(iun,*) n
!!$     write(iun,11) (ia(i),i=1,n)
!!$11   format(i6)
!!$     cio_i=0
!!$  endif
!!$
!!$  return
!!$end function cio_i
!!$
!!$!---------------------------------------------------------
!!$
!!$integer function cio_f(iun,irw,cstr,ia,n)
!!$  implicit none
!!$  integer :: iun,irw,n
!!$  real, dimension(n) :: ia
!!$  character(len=*) :: cstr
!!$  character(len=256) :: string
!!$  integer :: nn,i
!!$
!!$  if (irw.eq.1) then
!!$     call cio_pos_file (iun,cstr,cio_f)
!!$     if(cio_f.eq.1) return
!!$     read(iun,*) nn
!!$     read(iun,*) (ia(i),i=1,nn)
!!$  elseif(irw.eq.2) then
!!$     write(iun,20) cstr
!!$20   format('__',a)
!!$     write(iun,*) n
!!$     write(iun,11) (ia(i),i=1,n)
!!$11   format(e16.8)
!!$     cio_f=0
!!$  endif
!!$
!!$  return
!!$end function cio_f
!!$
!!$!---------------------------------------------------------
!!$
!!$integer function cio_f8(iun,irw,cstr,ia,n)
!!$  implicit none
!!$  integer :: iun,irw,n
!!$  real(kind=8) :: ia(*)
!!$  character(len=*) :: cstr
!!$  character(len=256) :: string
!!$  integer :: nn,i
!!$
!!$  if (irw.eq.1) then
!!$     call cio_pos_file (iun,cstr,cio_f8)
!!$     if(cio_f8.eq.1) return
!!$     read(iun,*) nn
!!$     read(iun,*) (ia(i),i=1,nn)
!!$  elseif(irw.eq.2) then
!!$     write(iun,20) cstr
!!$20   format('__',a)
!!$     write(iun,*) n
!!$     write(iun,11) (ia(i),i=1,n)
!!$11   format(e24.16)
!!$     cio_f8=0
!!$  endif
!!$
!!$  return
!!$end function cio_f8
!!$
!!$!---------------------------------------------------------
!!$
!!$integer function cio_c(iun,irw,cstr,ia,n)
!!$  implicit none
!!$  integer :: iun,irw,n
!!$  character(len=*) :: ia(*)
!!$  character(len=*) :: cstr
!!$  character(len=256) :: string
!!$  integer :: nn,i
!!$
!!$  if (irw.eq.1) then
!!$     call cio_pos_file (iun,cstr,cio_c)
!!$     if(cio_c.eq.1) return
!!$     read(iun,*) nn
!!$     read(iun,10) (ia(i),i=1,nn)
!!$  elseif(irw.eq.2) then
!!$     write(iun,20) cstr
!!$20   format('__',a)
!!$     write(iun,*) n
!!$     write(iun,10) (ia(i),i=1,n)
!!$10   format(a)
!!$     cio_c=0
!!$  endif
!!$
!!$  return
!!$end function cio_c
end module ModRcio
