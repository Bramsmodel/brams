module radiation

  use ModNamelistFile, only: &
       NamelistFile
  
  use ModRrtmDriver  , only: &
       rrtm_driver

  use ModCarmaDriver , only: &
       carma_driver

  use ModBasicFields, only: &
       BasicFields

  use ModMicControl, only: &
       MicControl

  use ModMicroFields, only: &
       MicroFields
  
  use ModRadiateFields, only: &
       RadiateFields
  
  use ModCuParmFields, only: &
       CuParmFields

  use ModGridDims, only: &
       GridDims
  
  implicit none

  private

  public :: radiate

contains


  subroutine radiate(mzp, mxp, myp, ia, iz, ja, jz, mynum, &
       oneNamelistFile, oneBasicFields, oneMicVars, oneMicroFields, &
       oneRadiateFields, oneCuParmFields, oneGridDims)

    integer, intent(in) :: mzp, mxp, myp, ia, iz, ja, jz, mynum
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicControl), pointer, intent(in) :: oneMicVars
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields
    type(CuParmFields), pointer, intent(in) :: oneCuParmFields
    type(GridDims), pointer, intent(in) :: oneGridDims

    character(len=*), parameter :: h="**(radiate)**"
    
    if ((oneNamelistFile%ilwrtyp + oneNamelistFile%iswrtyp)==0) then

       return

    else if (oneNamelistFile%ilwrtyp==6 .and. oneNamelistFile%iswrtyp==6) then

       call rrtm_driver(mzp, mxp, myp, ia, iz, ja, jz, mynum, &
            oneNamelistFile, oneBasicFields, oneMicVars, oneMicroFields, &
            oneRadiateFields, oneCuParmFields, oneGridDims)

    else if (oneNamelistFile%ilwrtyp==4 .and. oneNamelistFile%iswrtyp==4) then

       call carma_driver(mzp,mxp,myp,ia,iz,ja,jz,mynum, &
            oneNamelistFile, oneBasicFields, oneMicVars, oneMicroFields, &
            oneRadiateFields, oneCuParmFields) !teste 2

    else
       call fatal_error(h//" unknown radiation scheme")
    end if

  end subroutine radiate

end module radiation
