module radiation

  use mem_radiate, only: &
       ilwrtyp, &
       iswrtyp

  use ModRrtmDriver  , only: &
       rrtm_driver

  use ModCarmaDriver , only: &
       carma_driver

  use ModBasicFields, only: &
       BasicFields

  use ModMicControl, only: &
       MicControl

  implicit none

  private

  public :: radiate

contains


  subroutine radiate(mzp, mxp, myp, ia, iz, ja, jz, mynum, &
       oneBasicFields, oneMicControl)

    integer, intent(IN) :: mzp, mxp, myp, ia, iz, ja, jz, mynum
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicControl), pointer, intent(in) :: oneMicControl

    if &
         ((ilwrtyp + iswrtyp)==0) return ! teste

    if &
         ( ilwrtyp==6 .and. iswrtyp==6) &
         then

       call rrtm_driver(mzp, mxp, myp, ia, iz, ja, jz, mynum, &
            oneBasicFields, oneMicControl)

    elseif( ilwrtyp==4 .and. iswrtyp==4) then

       call carma_driver(mzp,mxp,myp,ia,iz,ja,jz,mynum, &
            oneBasicFields, oneMicControl) !teste 2

    else
       stop "unknown radiation scheme"
    endif

  end subroutine radiate

end module radiation
