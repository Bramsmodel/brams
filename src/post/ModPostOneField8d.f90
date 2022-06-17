module ModPostOneField8d

  use ModNamelistFile, only: &
       NamelistFile

  use ModTurbFields, only: &
       TurbFields

  use ModBasicFields, only: &
       BasicFields

  use ModOutputUtils, only: &
       GetVarFromMemToOutput

  use ModBramsGrid, only: &
       BramsGrid

  use ModPostTypes, only: &
       PostGrid

  use ModPostOneFieldUtils, only: &
       PrepareAndOutputGradsField

  use ModPostTypes, only: &
       PostVarType

  use ModPostUtils, only : &
       get_leaf_soil

  implicit none

  private

  public :: Brams2Post_8d

contains


  subroutine Brams2Post_8d (one_post_variable, oneBramsGrid, onePostGrid, &
       oneNamelistFile, oneBasicFields, oneTurbFields)
    type(PostVarType) :: one_post_variable
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields

    real :: ScrT4N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%nzg, oneBramsGrid%npatch)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%nzg, oneBramsGrid%npatch)

    select case (one_post_variable%fieldName)
    case ('SMOIST')
       call GetVarFromMemToOutput ('SOIL_WATER', oneBramsGrid%currGrid, ScrT4N01, &
            oneNamelistFile, oneBasicFields, oneTurbFields)
       call get_leaf_soil (ScrT4N01, OutputField)
    case ('SENERGY')
       call GetVarFromMemToOutput ('SOIL_ENERGY', oneBramsGrid%currGrid, ScrT4N01, &
            oneNamelistFile, oneBasicFields, oneTurbFields)
       call get_leaf_soil (ScrT4N01, OutputField)
    case ('SLTEX_P')
       call GetVarFromMemToOutput ('SOIL_TEXT', oneBramsGrid%currGrid, OutputField, &
            oneNamelistFile, oneBasicFields, oneTurbFields)
    case default
       write(*, "(a)") "**(OnePostField)** Post field 8d" // one_post_variable%fieldName // " not implemented!"
       stop
    end select

    ! most of cases run this. Some just returns before
    call PrepareAndOutputGradsField(one_post_variable, oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_8d


end module ModPostOneField8d
