Module ModPostgridNetCDF
#ifdef cdf
  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModNamelistFile, only: namelistFile

  use ModBramsGrid, only: BramsGrid

  use ModPostUtils, only: UpperCase, &
                          LowerCase

  use ModPostTypes, ONLY: PostVarType, &
                          postgrid,   &
                          getPostVariable, &
                          all_post_variables

  use mem_grid, only: time


  !For netCDF
  include "constants.h"
  character(len=256) :: netCDFFileName
  logical :: netCDFFirstTime
  integer :: ncid,LatDimID,LonDimID,LevDimID,timDimId,SurDimID,SoiDimID
  character(len=256), allocatable :: netCdfFieldName(:)
  integer, allocatable :: netCdfNLevCode(:)
  integer :: VarDimId(1000)
  character(len=256) :: netCdfFieldDescription(1000)
  character(len=256) :: netCdfFieldUnits(1000)
  real,allocatable :: hoursFrom1900(:)

  public :: FillNetcdfVarControlFile
  public :: netCDFFileName
  public :: netCDFFirstTime
  public :: OpenNetCDFBinaryFile

contains

#define CHECK_NF90(x) call check_nf_err(x, __LINE__, __FILE__)

subroutine check_nf_err(iErr, line, filename)

    use netcdf, only: nf90_noerr, &
      nf90_strerror
    implicit none

    INTEGER, INTENT(IN) :: iErr
    INTEGER, INTENT(IN) :: line
    CHARACTER(LEN=*), INTENT(IN) :: filename

    if (iErr /= nf90_noerr) then
        print *, "NF90 ERROR: ", trim(nf90_strerror(iErr))
        print *, "At line ", line, " in file ", filename
        call exit(iErr)
    end if

end subroutine check_nf_err

subroutine OpenNetCDFBinaryFile(oneNamelistFile, onePostGrid, oneBramsGrid, igrid) 
    use netCDF, ONLY: nf90_create,&
                      nf90_write, &
                      nf90_clobber, &
                      nf90_netcdf4
    IMPLICIT NONE
    type(NamelistFile), pointer :: oneNamelistFile
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid  
    integer, intent(in) :: igrid

    integer :: ierr

    character(len=1) :: c0

    write(c0,fmt='(I1)') igrid

    !if(.not. netCDFFirstTime) return
    if (oneBramsGrid%mchnum /= oneBramsGrid%master_num) return

    CALL makefnam(netCDFFileName, trim(oneNamelistFile%gprefix), time, &
         oneNamelistFile%iyear1, oneNamelistFile%imonth1,  &
         oneNamelistFile%idate1, oneNamelistFile%itime1*100, &
          'A', 'g'//c0, 'nc ')

    iErrNumber = nf90_create(path = trim(netCDFFileName), cmode = or(nf90_clobber,nf90_netcdf4), ncid = ncid)
    CHECK_NF90(iErrNumber)
end subroutine OpenNetCDFBinaryFile

subroutine FillNetcdfVarControlFile(oneNamelistFile, onePostGrid, oneBramsGrid)
    use mem_grid, only: oneGlobalGridData, nnxp, nnyp, &
                        iyear1,imonth1,idate1,ihour1, &
                        timmax, timeunit,npatch
    use io_params, only: frqanl
    use netcdf, ONLY: nf90_def_var, &
                      nf90_def_dim, &
                      nf90_redef,   &
                      nf90_put_att, &
                      nf90_enddef,  &
                      nf90_put_var, &
                      nf90_float

    use dump, only: dumpMessage
    use ModDateUtils, ONLY: date_abs_secs2
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    implicit none

    type(NamelistFile), pointer :: oneNamelistFile
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostVarType) :: one_post_variable

    include "constants.h"

    integer :: i,ndims,nvars,cnt
    character(len=256) :: name,varname(30),atName(300),atValue(300)
    character(len = 16) :: varNameUpper,vName,vNameTmp

    integer, dimension(300) :: lenDim,nat
    integer :: vdim,levs,an,ntimes,id,ivp
    real :: slayer(oneBramsGrid%nzg)
    real(kind=r8) :: seconds
    real :: lon(nnxp(1)),lat(nnyp(1))
    character(len=8) :: cVar
    real :: nan
    integer :: all_post_idx

    integer, parameter :: deflate_level = 1
    character(len=*), parameter :: h='**(FillNetcdfVarControlFile)**'
    logical, parameter :: dumpLocal=.false.

    !if(.not. netCDFFirstTime) return
    if (oneBramsGrid%mchnum /= oneBramsGrid%master_num) return

    nan = IEEE_Value(nan, IEEE_QUIET_NAN)

    iErrNumber = nf90_def_dim(ncid, "longitude", nnxp(1), LonDimID)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_dim(ncid, "latitude", nnyp(1), LatDimID)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_dim(ncid, "level",onePostGrid%nVert, LevDimID)
    CHECK_NF90(iErrNumber)
    ntimes=1!timmax/frqanl
    iErrNumber = nf90_def_dim(ncid, "time", ntimes, TimDimID)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_dim(ncid, "soil_level", oneBramsGrid%nzg, SoiDimID)
    CHECK_NF90(iErrNumber)

    ntimes=1
    if (.not. allocated(hoursFrom1900)) allocate(hoursFrom1900(nTimes))

    ! Define the variable Longitude
    iErrNumber = nf90_def_var(ncid, "longitude", nf90_float, &
    (/LonDimId/), LonDimId)
    CHECK_NF90(iErrNumber)

    ! Define the variable Latitude
    iErrNumber = nf90_def_var(ncid, "latitude", nf90_float, &
    (/LatDimId/), LatDimId)
    CHECK_NF90(iErrNumber)

    ! Define the variable Level
    iErrNumber = nf90_def_var(ncid, "level", nf90_float, &
    (/LevDimId/), LevDimId)
    CHECK_NF90(iErrNumber)

    ! Define the variable Level
    iErrNumber = nf90_def_var(ncid, "time", nf90_float, &
    (/TimDimId/), TimDimId)
    CHECK_NF90(iErrNumber)

    ! Define the variable Level
    iErrNumber = nf90_def_var(ncid, "soil_level", nf90_float, &
    (/SoiDimId/), SoiDimId)
    CHECK_NF90(iErrNumber)

    iErrNumber = nf90_redef(ncid)
    CHECK_NF90(iErrNumber)

    iErrNumber = nf90_put_att(ncid, LonDimId, "units", "degrees_east")
    iErrNumber = nf90_put_att(ncid, LonDimId, "long name", "longitude")
    iErrNumber = nf90_put_att(ncid, LatDimId, "units", "degrees_north")
    iErrNumber = nf90_put_att(ncid, LatDimId, "long name", "latitude")
    iErrNumber = nf90_put_att(ncid, LevDimId, "units", "mBar")
    iErrNumber = nf90_put_att(ncid, LevDimId, "long name", "pressure level")
    iErrNumber = nf90_put_att(ncid, TimDimId, "units", "hours since 1900-01-01 00:00:00")
    iErrNumber = nf90_put_att(ncid, TimDimId, "long name", "time")   
    iErrNumber = nf90_put_att(ncid, TimDimId, "calendar", "gregorian")
    iErrNumber = nf90_put_att(ncid, soiDimId, "units", "m")
    iErrNumber = nf90_put_att(ncid, soiDimId, "long name", "soil level")

    cnt=0
    do ivp = 1, oneNamelistFile%nvp
      vName=oneNamelistFile%vp(ivp)
      varNameUpper = trim(UpperCase(vName))
      one_post_variable = getPostVariable(varNameUpper, all_post_idx)
      vName = trim(LowerCase(vName))

      if (len(trim(one_post_variable%fieldName)) .eq. 0) then
         write(logUnit, "(a)") "**(OnePostField)** Post field " // varName // " does not exists in list of variables."
         write(logUnit, "(a)") "Model will continue ..."
      else
         select case (one_post_variable%ivar_type)
         case (2) !Surface var
            cnt=cnt+1
            netCdfFieldDescription(cnt)=one_post_variable%fieldDescription
            netCdfFieldUnits(cnt)=one_post_variable%fieldUnits      
            iErrNumber = nf90_def_var(ncid,name=vName, xtype=nf90_float, &
                  dimids=(/LonDimId,LatDimId,TimDimId/),varid=VarDimId(cnt), &
                  shuffle=.true., deflate_level=deflate_level)
            CHECK_NF90(iErrNumber)
            all_post_variables(all_post_idx)%netcdfId(1) = VarDimId(cnt)

         case (3) !Athmos var
            cnt=cnt+1
            netCdfFieldDescription(cnt)=one_post_variable%fieldDescription
            netCdfFieldUnits(cnt)=one_post_variable%fieldUnits
            iErrNumber = nf90_def_var(ncid,name=vName, xtype=nf90_float, &
                  dimids=(/LonDimId,LatDimId,LevDimId,TimDimId/),varid=VarDimId(cnt), &
                  shuffle=.true., deflate_level=deflate_level)
            CHECK_NF90(iErrNumber)
            all_post_variables(all_post_idx)%netcdfId(1) = VarDimId(cnt)

         case (7) !Surface Soil var

            do i=1,npatch
              cnt=cnt+1
              write(cvar,fmt='(I8)') i 
              vNameTmp=trim(vName)//trim(adjustl(cvar))
              netCdfFieldDescription(cnt)=trim(one_post_variable%fieldDescription)//': patch #'//trim(cvar)
              netCdfFieldUnits(cnt)=one_post_variable%fieldUnits
              iErrNumber = nf90_def_var(ncid,name=vNameTmp, xtype=nf90_float, &
              dimids=(/LonDimId,LatDimId,TimDimId/),varid=VarDimId(cnt), &
                  shuffle=.true., deflate_level=deflate_level)
              CHECK_NF90(iErrNumber)
              all_post_variables(all_post_idx)%netcdfId(i) = VarDimId(cnt)

            enddo
         case (8) !

            do i=1,npatch
              cnt=cnt+1
              write(cvar,fmt='(I8)') i 
              vNameTmp=trim(vName)//trim(adjustl(cvar))
              netCdfFieldDescription(cnt)=trim(one_post_variable%fieldDescription)//': patch #'//trim(cvar)
              netCdfFieldUnits(cnt)=one_post_variable%fieldUnits
              iErrNumber = nf90_def_var(ncid,name=vNameTmp, xtype=nf90_float, &
                  dimids=(/LonDimId,LatDimId,SoiDimId,TimDimId/),varid=VarDimId(cnt), &
                  shuffle=.true., deflate_level=deflate_level)
              CHECK_NF90(iErrNumber)
              all_post_variables(all_post_idx)%netcdfId(i) = VarDimId(cnt)

            enddo
         end select

      end if
    enddo

    do i=1,cnt

       iErrNumber = nf90_put_att(ncid, VarDimId(i), "_FillValue", nan)
       CHECK_NF90(iErrNumber)
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "missing_value", nan)
       CHECK_NF90(iErrNumber)
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "units" &
                 ,trim(netCdfFieldUnits(i)))
       CHECK_NF90(iErrNumber)
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "long_name" &
                 , trim(netCdfFieldDescription(i)))
       CHECK_NF90(iErrNumber)
    enddo

    iErrNumber = nf90_enddef(ncid)
    CHECK_NF90(iErrNumber)

    !Fiiling lons, lats and pressure levels 
    do i=1,nnyp(1)
      lat(i)=oneGlobalGridData(1)%global_glat(1,nnyp(1)-i+1)
    end do

    !Making glon from 0 to 360 east direction
    do i=1,nnxp(1)
      !if(oneGlobalGridData(1)%global_glon(i,1)<0) &
        !lon(i)=360+oneGlobalGridData(1)%global_glon(i,1)
        lon(i)=oneGlobalGridData(1)%global_glon(i,1)
    enddo

    iErrNumber = nf90_put_var(ncid, LonDimId, lon)
    CHECK_NF90(iErrNumber)

    iErrNumber = nf90_put_var(ncid, LevDimId, onePostGrid%vertScaleValues)
    CHECK_NF90(iErrNumber)

    !Filling Soil layers - Just putting a number to represent it
    do i=1,oneBramsGrid%nzg
      slayer(i)=real(i)
    enddo
    iErrNumber = nf90_put_var(ncid, SoiDimId, slayer)
    CHECK_NF90(iErrNumber)
    
    !Filling Times
    !First convert the initial date to seconds
    if (dumpLocal) then
       call MsgDump(h//" invokes date_abs_secs2")
    end if
    call date_abs_secs2(iyear1,imonth1,idate1,ihour1,seconds)
    !Now put increments of frqanl in hours inside array
    do i=1,ntimes
      hoursFrom1900(i)=(((frqanl*(i-1))+seconds)/3600.)-24
    enddo
    !Finally put the array in TimDimId
    iErrNumber = nf90_put_var(ncid, TimDimId, hoursFrom1900)
    CHECK_NF90(iErrNumber)

    netCDFFirstTime=.false.

  end subroutine FillNetcdfVarControlFile

  subroutine netCdfPostField2D(fieldName,nLon,nLat,OutputArray,netcdfIdIndex)

    use netcdf, ONLY: nf90_put_var
    use dump, only: dumpMessage
    use mem_grid, only: time

    include "constants.h"
    integer, intent(in) :: nlon,nlat
    character(len=*), intent(in) :: fieldName
    real, intent(in) :: OutputArray(nlon, nlat)
    integer, intent(in) :: netcdfIdIndex

    integer :: id, nt
    character(len=256) :: varName(500),name
    real :: varArray(nlon, nlat)
    type(PostVarType) :: one_post_variable

    one_post_variable = getPostVariable(fieldName)
    call invertLats(OutputArray,varArray,nlon,nlat)

    iErrNumber = nf90_put_var(ncid,one_post_variable%netcdfId(netcdfIdIndex),varArray,start=(/1, 1, 1/))
    CHECK_NF90(iErrNumber)

  end subroutine netCdfPostField2D

  subroutine netCdfPostField3D(fieldName,nLon,nLat,ilev,OutputArray,netcdfIdIndex)

    use netcdf, ONLY: nf90_put_var
    use dump, only: dumpMessage
    use mem_grid, only: time

    include "constants.h"
    integer, intent(in) :: nlon,nlat,iLev
    character(len=*), intent(in) :: fieldName
    real, intent(in) :: OutputArray(nlon,nlat)
    integer, intent(in) :: netcdfIdIndex

    integer :: id,nt
    character(len=256) :: varName(500),name
    real :: varArray(nlon,nlat)
    type(PostVarType) :: one_post_variable

    iErrNumber = nf90_put_var(ncid, LatDimId, lat)
    CHECK_NF90(iErrNumber)


    one_post_variable = getPostVariable(fieldName)
    call invertLats(OutputArray,varArray,nlon,nlat)

    iErrNumber = nf90_put_var(ncid,one_post_variable%netcdfId(netcdfIdIndex),varArray,start=(/1,1,iLev,1/))
    CHECK_NF90(iErrNumber)

  end subroutine netCdfPostField3D

  subroutine invertLats(iArray,oArray,nLon,nLat)

    integer,intent(in) :: nLon,nLat
    real,intent(in) :: iArray(nLon,nLat)
    real,intent(out) :: oArray(nLon,nLat)

    integer :: i,j

    do i=1,nLon
      do j=1,nLat
        oArray(i,j)=iArray(i,nLat-j+1)
      enddo
    enddo

  end subroutine invertLats
#endif
end Module ModPostgridNetCDF

