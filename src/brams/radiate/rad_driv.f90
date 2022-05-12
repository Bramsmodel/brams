
MODULE RADIATION

CONTAINS

  SUBROUTINE radiate(mzp, mxp, myp, ia, iz, ja, jz, mynum, oneBasicFields)

    USE mem_radiate, only: ilwrtyp, iswrtyp
    USE rrtm_driv  , only: rrtm_driver
    USE ModCarmaDriver , only: carma_driver
    use ModBasicFields, only: &
         BasicFields

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: mzp, mxp, myp, ia, iz, ja, jz, mynum
    type(BasicFields), pointer, intent(in) :: oneBasicFields

    if &
       ((ilwrtyp + iswrtyp)==0) return ! teste

    if &
      ( ilwrtyp==6 .and. iswrtyp==6) &
      then

       call rrtm_driver(mzp, mxp, myp, ia, iz, ja, jz, mynum)

    elseif( ilwrtyp==4 .and. iswrtyp==4) then

       call carma_driver(mzp,mxp,myp,ia,iz,ja,jz,mynum, oneBasicFields) !teste 2

    else
       stop "unknown radiation scheme"
    endif

  END SUBROUTINE radiate

END MODULE RADIATION
