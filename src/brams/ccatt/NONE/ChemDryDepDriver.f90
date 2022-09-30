module ChemDryDepDriver

  use ModNamelistFile, only: &
       NamelistFile

  use rconstants, only: &
       cpi,             &
       cpor,            &
       p00,             &
       g,               &
       vonk

  use mem_grid, only: &
       grid_g,        &
       jdim,          &
       dzt,           &
       zt,            &
       nzpmax,        &
       npatch,        &
       dtlt,          &
       imonth1,       &
       idate1,        &
       iyear1,        &
       ngrid

  use ModMicControl, only: &
       MicControl

  use mem_cuparm, only: &
       cuparm_g

  use ModBasicFields, only: &
       BasicFields

  use ModMicroFields, only: &
       MicroFields
  
  use ModTurbFields, only: &
       TurbFields

  use mem_leaf, only: &
       leaf_g

  use ModRadiateFields, only: &
       RadiateFields

  use mem_chem1, only: &
       chem1_g,        &
       chemistry

  use mem_aer1, only: &
       aerosol, &
       aer1_g

  use module_dry_dep, only: &
       dd_sedim,            &
       dry_dep                 ! Subroutine



  implicit none

  private


  public :: drydep_driver



contains

  !========================================================================
  subroutine drydep_driver(m1,m2,m3,ia,iz,ja,jz, &
       oneNamelistFile, &
       oneBasicFields, &
       oneTurbFields, &
       oneMicControl, &
       oneMicroFields, &
       oneRadiateFields)

    integer,              intent(IN)    :: m1
    integer,              intent(IN)    :: m2
    integer,              intent(IN)    :: m3
    integer,              intent(IN)    :: ia
    integer,              intent(IN)    :: iz
    integer,              intent(IN)    :: ja
    integer,              intent(IN)    :: jz
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields
    type(MicControl), pointer, intent(in) :: oneMicControl
    type(MicroFields), pointer, intent(in) :: oneMicroFields
    type(RadiateFields), pointer, intent(in) :: oneRadiateFields


    return
  end subroutine drydep_driver
  !========================================================================


end module ChemDryDepDriver
