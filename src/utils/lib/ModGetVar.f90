!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModGetVar

  use dump, only: &
       dumpMessage

  use an_header,only:&
       anal_table,    &
       nvbtab

  implicit none

  include "files.h"
  include "constants.h"

  private

  public :: RAMS_getvar
  public :: FindFieldInAnalysisFile
  public :: GetFieldInAnalysisFile

  interface RAMS_getvar
     module procedure RAMS_getvar_1D
     module procedure RAMS_getvar_2D
     module procedure RAMS_getvar_3D
  end interface RAMS_getvar
contains



  integer function RAMS_getvar_1D(string, ngrd, a, b, flnm)
    ! Arguments:
    character(len=*), intent(in) :: string
    integer, intent(in)          :: ngrd
    real, intent(out)            :: a(:)
    real, intent(inout)          :: b(:)
    character(len=*), intent(in) :: flnm
    ! Local variables:
    integer                      :: itype !!, rams_c_pos
    character(LEN=1)             :: cgrid
    character(LEN=f_name_length) :: flng
    character(LEN=120)           :: errmsg
    logical                      :: there
    integer                      :: ni
    integer(kind=i8)             :: npts, iword

    character(len=*), parameter :: h="**(RAMS_getvar_1D)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       print*, h//' getvar:',string
    end if

    do ni=1,nvbtab

       if(string==anal_table(ni)%string .and. ngrd==anal_table(ni)%ngrid) then

          write(cgrid,'(i1)') ngrd
          flng=trim(flnm)//'-g'//cgrid//'.vfm'

          inquire(file=flng(1:len_trim(flng)),exist=there)
          if(.not.there) then
             errmsg='File not found - '//trim(flng)
             RAMS_getvar_1D = 1
             call error_mess(errmsg)
             return
          endif

          npts=anal_table(ni)%nvalues
          itype=anal_table(ni)%idim_type
          iword=anal_table(ni)%npointer

          !  print*,'gv:opening'
          call RAMS_c_open(flng(1:len_trim(flng))//char(0),'r'//char(0))
          !  print*,'gv:opened'
          call vfirecr(10,a,npts,'LIN',b,iword)
          !  print*,'gv:vfirecr'
          call RAMS_c_close()
          !  print*,'gv:closed'

          RAMS_getvar_1D=0
          if (dumpLocal) then
             print*,'getvar good:',string
          end if
          return

       endif
    enddo

    errmsg='Variable not available in this run - '//string
    call error_mess(errmsg)
    RAMS_getvar_1D=1

    return
  end function RAMS_getvar_1D
  


  integer function RAMS_getvar_2D(string, ngrd, a, b, flnm)
    ! Arguments:
    character(len=*), intent(in) :: string
    integer, intent(in)          :: ngrd
    real, intent(out)            :: a(:,:)
    real, intent(inout)          :: b(:)
    character(len=*), intent(in) :: flnm
    ! Local variables:
    integer                      :: itype !!, rams_c_pos
    character(LEN=1)             :: cgrid
    character(LEN=f_name_length) :: flng
    character(LEN=120)           :: errmsg
    logical                      :: there
    integer                      :: ni
    integer(kind=i8)             :: npts, iword

    character(len=*), parameter :: h="**(RAMS_getvar_2D)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       print*, h//' getvar:',string
    end if

    do ni=1,nvbtab

       if(string==anal_table(ni)%string .and. ngrd==anal_table(ni)%ngrid) then

          write(cgrid,'(i1)') ngrd
          flng=trim(flnm)//'-g'//cgrid//'.vfm'

          inquire(file=flng(1:len_trim(flng)),exist=there)
          if(.not.there) then
             errmsg='File not found - '//trim(flng)
             RAMS_getvar_2D = 1
             call error_mess(errmsg)
             return
          endif

          npts=anal_table(ni)%nvalues
          itype=anal_table(ni)%idim_type
          iword=anal_table(ni)%npointer

          !  print*,'gv:opening'
          call RAMS_c_open(flng(1:len_trim(flng))//char(0),'r'//char(0))
          !  print*,'gv:opened'
          call vfirecr(10,a,npts,'LIN',b,iword)
          !  print*,'gv:vfirecr'
          call RAMS_c_close()
          !  print*,'gv:closed'

          RAMS_getvar_2D=0
          if (dumpLocal) then
             print*,'getvar good:',string
          end if
          return

       endif
    enddo

    errmsg='Variable not available in this run - '//string
    call error_mess(errmsg)
    RAMS_getvar_2D=1

    return
  end function RAMS_getvar_2D
  


  integer function RAMS_getvar_3D(string, ngrd, a, b, flnm)
    ! Arguments:
    character(len=*), intent(in) :: string
    integer, intent(in)          :: ngrd
    real, intent(out)            :: a(:,:,:)
    real, intent(inout)          :: b(:)
    character(len=*), intent(in) :: flnm
    ! Local variables:
    integer                      :: itype !!, rams_c_pos
    character(LEN=1)             :: cgrid
    character(LEN=f_name_length) :: flng
    character(LEN=120)           :: errmsg
    logical                      :: there
    integer                      :: ni
    integer(kind=i8)             :: npts, iword

    character(len=*), parameter :: h="**(RAMS_getvar_3D)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       print*, h//' getvar:',string
    end if

    do ni=1,nvbtab

       if(string==anal_table(ni)%string .and. ngrd==anal_table(ni)%ngrid) then

          write(cgrid,'(i1)') ngrd
          flng=trim(flnm)//'-g'//cgrid//'.vfm'

          inquire(file=flng(1:len_trim(flng)),exist=there)
          if(.not.there) then
             errmsg='File not found - '//trim(flng)
             RAMS_getvar_3D = 1
             call error_mess(errmsg)
             return
          endif

          npts=anal_table(ni)%nvalues
          itype=anal_table(ni)%idim_type
          iword=anal_table(ni)%npointer

          !  print*,'gv:opening'
          call RAMS_c_open(flng(1:len_trim(flng))//char(0),'r'//char(0))
          !  print*,'gv:opened'
          call vfirecr(10,a,npts,'LIN',b,iword)
          !  print*,'gv:vfirecr'
          call RAMS_c_close()
          !  print*,'gv:closed'

          RAMS_getvar_3D=0
          if (dumpLocal) then
             print*,'getvar good:',string
          end if
          return

       endif
    enddo

    errmsg='Variable not available in this run - '//string
    call error_mess(errmsg)
    RAMS_getvar_3D=1

    return
  end function RAMS_getvar_3D
  


  

  subroutine FindFieldInAnalysisFile(string, ngrd, flnm, flng, npts, fPosition)
    ! Arguments:
    character(len=*), intent(in ) :: string
    integer,          intent(in ) :: ngrd
    character(len=*), intent(in ) :: flnm
    character(len=*), intent(out) :: flng
    integer(kind=i8), intent(out) :: npts
    integer(kind=i8), intent(out) :: fPosition

    ! Local variables:
    character(LEN=1)  :: cgrid
    logical           :: there
    integer           :: ni
    character(len=*), parameter  :: h="**(FindFieldInAnalysisFile)**"

    do ni=1,nvbtab

       if(string==anal_table(ni)%string .and. ngrd==anal_table(ni)%ngrid) then

          write(cgrid,'(i1)') ngrd
          flng=trim(flnm)//'-g'//cgrid//'.vfm'

          inquire(file=flng(1:len_trim(flng)),exist=there)
          if(.not.there) then
             !call fatal_error(h//'File not found - '//trim(flng)//&
             !    ' while searching for field '//trim(string)//&
             !    ' at grid '//cgrid)
             iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
                  'File not found - '//trim(flng)//&
                  ' while searching for field '//trim(string)//&
                  ' at grid '//cgrid)
          else
             npts=anal_table(ni)%nvalues
             fPosition=anal_table(ni)%npointer
             return
          end if
       end if
    end do

    !call fatal_error(h//'Variable '//trim(string)//&
    !     &' not available in the analysis file')
    iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
         'Variable '//trim(string)//&
         ' not available in the analysis file')

  end subroutine FindFieldInAnalysisFile



  subroutine GetFieldInAnalysisFile(flng, npts, fPosition, a, b)
    character(len=*),             intent(in   ) :: flng
    integer(kind=i8),             intent(in   ) :: npts
    integer(kind=i8),             intent(in   ) :: fPosition
    real,                         intent(inout) :: a(:)
    real,                         intent(inout) :: b(:)

    call RAMS_c_open(trim(flng)//char(0),'r'//char(0))

    call vfirecr(10,a,npts,'LIN',b,fPosition)

    call RAMS_c_close()
  end subroutine GetFieldInAnalysisFile
end module ModGetVar
