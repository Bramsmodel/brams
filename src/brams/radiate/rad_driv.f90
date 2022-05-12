
module RADIATION

contains

  
  subroutine radiate(mzp, mxp, myp, ia, iz, ja, jz, mynum, oneBasicFields)

    use mem_radiate, only: ilwrtyp, iswrtyp
    use ModRrtmDriver  , only: rrtm_driver
    use ModCarmaDriver , only: carma_driver
    use ModBasicFields, only: &
         BasicFields

    implicit none
    integer, intent(IN) :: mzp, mxp, myp, ia, iz, ja, jz, mynum
    type(BasicFields), pointer, intent(in) :: oneBasicFields

    if &
       ((ilwrtyp + iswrtyp)==0) return ! teste

    if &
      ( ilwrtyp==6 .and. iswrtyp==6) &
      then

       call rrtm_driver(mzp, mxp, myp, ia, iz, ja, jz, mynum, oneBasicFields)

    elseif( ilwrtyp==4 .and. iswrtyp==4) then

       call carma_driver(mzp,mxp,myp,ia,iz,ja,jz,mynum, oneBasicFields) !teste 2

    else
       stop "unknown radiation scheme"
    endif

  end subroutine radiate

end module RADIATION
