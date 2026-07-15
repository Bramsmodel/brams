Module ModPostgridNetCDF
#ifdef cdf
  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModNamelistFile, only: namelistFile

  use ModBramsGrid, only: BramsGrid

  use ModPostUtils, only: UpperCase, &
                          LowerCase, &
                          undef

  use ModPostTypes, ONLY: PostVarType, &
                          postgrid,   &
                          getPostVariable, &
                          all_post_variables


  implicit none

  !For netCDF
  character(len=256) :: netCDFFileName
  logical :: netCDFFirstTime
  integer :: ncid,LatDimId,LonDimId,LevDimId,timDimId,SoiDimId
  integer :: LonVarId, LatVarId, LevVarId, TimVarId, SoiVarId
  integer :: VarDimId(1000)
  character(len=256) :: netCdfFieldDescription(1000)
  character(len=256) :: netCdfFieldUnits(1000)

  public :: FillNetcdfVarControlFile
  public :: netCDFFileName
  public :: netCDFFirstTime
  public :: OpenNetCDFBinaryFile
  public :: netCdfPostField2D
  public :: netCdfPostField3D

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
    use mem_grid, only: time

    IMPLICIT NONE
    type(NamelistFile), pointer :: oneNamelistFile
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid  
    integer, intent(in) :: igrid

    integer :: iErrNumber
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
    use mem_grid, only: oneGlobalGridData,npatch,time, &
                        iyear1,imonth1,idate1,ihour1
    use mem_leaf, only: slz
    use netcdf, only: nf90_def_var, &
                      nf90_def_dim, &
                      nf90_redef,   &
                      nf90_put_att, &
                      nf90_enddef,  &
                      nf90_put_var, &
                      nf90_float,   &
                      nf90_int,     &
                      nf90_set_fill,&
                      nf90_nofill

    use dump, only: dumpMessage
    use ModDateUtils, ONLY: date_abs_secs2
    implicit none

    type(NamelistFile), pointer :: oneNamelistFile
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostVarType) :: one_post_variable

    integer :: iErrNumber
    integer :: i,cnt,ivp
    character(len = 16) :: varNameUpper,vName,vNameTmp

    character(len=8) :: cVar
    character(len=24) :: datestring
    integer :: all_post_idx

    integer, parameter :: deflate_level = 1
    integer, parameter :: ntimes = 1
    character(len=*), parameter :: h='**(FillNetcdfVarControlFile)**'
    logical, parameter :: dumpLocal=.false.

    !if(.not. netCDFFirstTime) return
    if (oneBramsGrid%mchnum /= oneBramsGrid%master_num) return

    iErrNumber = nf90_set_fill(ncid, nf90_nofill, i)
    CHECK_NF90(iErrNumber)

    iErrNumber = nf90_def_dim(ncid, "longitude", onePostGrid%nLon, LonDimId)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_dim(ncid, "latitude", onePostGrid%nLat, LatDimId)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_dim(ncid, "pressure_level",onePostGrid%nVert, LevDimId)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_dim(ncid, "soil_level", oneBramsGrid%nzg, SoiDimId)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_dim(ncid, "time", ntimes, TimDimID)
    CHECK_NF90(iErrNumber)

    ! Define the variables of the dimensions
    iErrNumber = nf90_def_var(ncid, "longitude", nf90_float, (/LonDimId/), LonVarId)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_var(ncid, "latitude", nf90_float, (/LatDimId/), LatVarId)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_var(ncid, "pressure_level", nf90_float, (/LevDimId/), LevVarId)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_var(ncid, "soil_level", nf90_float, (/SoiDimId/), SoiVarId)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_def_var(ncid, "time", nf90_int, (/TimDimId/), TimVarId)
    CHECK_NF90(iErrNumber)

    iErrNumber = nf90_put_att(ncid, LonVarId, "units", "degrees_east")
    iErrNumber = nf90_put_att(ncid, LonVarId, "long_name", "longitude")
    iErrNumber = nf90_put_att(ncid, LatVarId, "units", "degrees_north")
    iErrNumber = nf90_put_att(ncid, LatVarId, "long_name", "latitude")
    iErrNumber = nf90_put_att(ncid, LevVarId, "units", "mb")
    iErrNumber = nf90_put_att(ncid, LevVarId, "long_name", "pressure level")

    write(datestring, '(I4.4, "-", I2.2, "-", I2.2, " ", I2.2, ":00:00")') iyear1, imonth1, idate1, ihour1

    iErrNumber = nf90_put_att(ncid, TimVarId, "units", "seconds since "//datestring)
    iErrNumber = nf90_put_att(ncid, TimVarId, "long_name", "time")
    iErrNumber = nf90_put_att(ncid, TimVarId, "calendar", "gregorian")
    iErrNumber = nf90_put_att(ncid, SoiVarId, "units", "m")
    iErrNumber = nf90_put_att(ncid, SoiVarId, "long_name", "soil level")

    cnt=0
    do ivp = 1, oneNamelistFile%nvp
      vName=oneNamelistFile%vp(ivp)
      varNameUpper = trim(UpperCase(vName))
      one_post_variable = getPostVariable(varNameUpper, all_post_idx)
      vName = trim(LowerCase(vName))

      if (len(trim(one_post_variable%fieldName)) .eq. 0) then
         print *,"**(OnePostField)** Post field " // vName // " does not exists in list of variables."
         print *,"Model will continue ..."
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

       iErrNumber = nf90_put_att(ncid, VarDimId(i), "_FillValue", undef)
       CHECK_NF90(iErrNumber)
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "missing_value", undef)
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

    iErrNumber = nf90_put_var(ncid, LatVarId, onePostGrid%lat)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_put_var(ncid, LonVarId, onePostGrid%lon)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_put_var(ncid, LevVarId, onePostGrid%vertScaleValues)
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_put_var(ncid, SoiVarId, slz(1:oneBramsGrid%nzg))
    CHECK_NF90(iErrNumber)
    iErrNumber = nf90_put_var(ncid, TimVarId, time)
    CHECK_NF90(iErrNumber)

    netCDFFirstTime=.false.

end subroutine FillNetcdfVarControlFile

subroutine netCdfPostField2D(fieldName,nLon,nLat,OutputArray,netcdfIdIndex)

    use netcdf, ONLY: nf90_put_var
    use dump, only: dumpMessage

    integer, intent(in) :: nlon,nlat
    character(len=*), intent(in) :: fieldName
    real, intent(in) :: OutputArray(nlon, nlat)
    integer, intent(in) :: netcdfIdIndex

    integer :: iErrNumber
    type(PostVarType) :: one_post_variable

    one_post_variable = getPostVariable(fieldName)

    iErrNumber = nf90_put_var(ncid,one_post_variable%netcdfId(netcdfIdIndex),OutputArray,start=(/1, 1, 1/))
    CHECK_NF90(iErrNumber)

end subroutine netCdfPostField2D

subroutine netCdfPostField3D(fieldName,nLon,nLat,ilev,OutputArray,netcdfIdIndex)

    use netcdf, ONLY: nf90_put_var
    use dump, only: dumpMessage

    integer, intent(in) :: nlon,nlat,iLev
    character(len=*), intent(in) :: fieldName
    real, intent(in) :: OutputArray(nlon,nlat)
    integer, intent(in) :: netcdfIdIndex

    integer :: iErrNumber
    type(PostVarType) :: one_post_variable

    one_post_variable = getPostVariable(fieldName)

    iErrNumber = nf90_put_var(ncid,one_post_variable%netcdfId(netcdfIdIndex),OutputArray,start=(/1,1,iLev,1/))
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

