Module ModPostOneFieldNetcdf
#ifdef cdf
  use ModPostGridNetCDF, only: &
      netCdfPostField2D, &
      netCdfPostField3D

  implicit none
  
  contains

  subroutine writeNetCdf2D(fieldName,nLon,nLat,OutputArray,netcdfIdIndex)
  
      integer, intent(in) :: nlon,nlat
      character(len=*), intent(in) :: fieldName
      real, intent(in) :: OutputArray(nlon, nlat)
      integer, intent(in) :: netcdfIdIndex
      
      call netCdfPostField2D(fieldName,nLon,nLat,OutputArray,netcdfIdIndex)
  
  end subroutine writeNetCdf2D
  
  subroutine writeNetCdf3D(fieldName,nLon,nLat,iLev,OutputArray,netcdfIdIndex)
    integer, intent(in) :: nlon,nlat,iLev
    character(len=*), intent(in) :: fieldName
    real, intent(in) :: OutputArray(nlon, nlat)
    integer, intent(in) :: netcdfIdIndex
  
    call netCdfPostField3D(fieldName,nLon,nLat,iLev,OutputArray,netcdfIdIndex)
  
  end subroutine writeNetCdf3D
#endif
end Module ModPostOneFieldNetcdf
