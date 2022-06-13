module ModPostOneField

  use ModTurbFields, only: &
       TurbFields
  
  use ModBasicFields, only: &
       BasicFields

  use ModBramsGrid, only: &
       BramsGrid

  use ModPostTypes, only: &
       PostGrid, &
       PostVarType, &
       getPostVarible

  use ModPostOneField2d, only: &
       Brams2Post_2d

  use ModPostOneField3d, only: &
       Brams2Post_3d

  use ModPostOneField7d, only: &
       Brams2Post_7d

  use ModPostOneField8d, only: &
       Brams2Post_8d

  use ModPostOneFieldUtils, only: &
       initialize_all_post_variables, &
       finalize_all_post_variables, &
       add_post_variable

  use ModPostUtils, only: &
       UpperCase

  use ModNamelistFile, only: &
       namelistFile

  use dump, only: &
       dumpMessage

  use node_mod, only: &
       mchnum,        &
       master_num


  implicit none

  private

  public :: PostOneField
  public :: initialize_post_variables
  public :: finalize_post_variables


contains

  ! TODO - Read the variables from a file to exchange with Brabu
  subroutine initialize_post_variables(oneNamelistFile)
    type(NamelistFile), pointer :: oneNamelistFile
    type(PostVarType) :: postVar
    integer :: varfileUnit, err

    call initialize_all_post_variables()
    varfileUnit = 135
    open(varfileUnit, file = oneNamelistFile%csvFile, status = "old", action = "read")
    do
       read(varfileUnit, *, iostat = err) postVar%ivar_type, postVar%fieldName, postVar%fieldDescription, postVar%fieldUnits
       if (err /= 0) then
          exit
       end if
       call add_post_variable(postVar%ivar_type, postVar%fieldName, postVar%fieldDescription, postVar%fieldUnits)
    enddo
    close(varfileUnit)

  end subroutine initialize_post_variables



  
  subroutine finalize_post_variables()
    call finalize_all_post_variables()
  end subroutine finalize_post_variables


  

  subroutine PostOneField(varName, oneBramsGrid, onePostGrid, &
       oneNamelistFile, oneBasicFields, oneTurbFields)
    include "constants.h"
    character(len = *), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields

    type(PostVarType) :: one_post_variable
    character(len = 16) :: varNameUpper
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(PostOneField)**"
    integer :: err

    varNameUpper = trim(UpperCase(varName))
    one_post_variable = getPostVarible(varNameUpper)
    if (len(trim(one_post_variable%fieldName)) .eq. 0) then
       if (mchnum==master_num) &
            err=dumpMessage(c_tty,c_yes,'Post','',c_warning,'Post field "' // varName &
            // '" does not exists in list of variables. Model will continue!')
    else
       select case (one_post_variable%ivar_type)
       case (2)
          call Brams2Post_2d(one_post_variable, oneBramsGrid, onePostGrid, &
               oneNamelistFile, oneBasicFields, oneTurbFields)
       case (3)
          call Brams2Post_3d(one_post_variable, oneBramsGrid, onePostGrid, &
               oneNamelistFile, oneBasicFields, oneTurbFields)
       case (7)
          call Brams2Post_7d(one_post_variable, oneBramsGrid, onePostGrid, &
               oneNamelistFile, oneBasicFields, oneTurbFields)
       case (8)
          call Brams2Post_8d(one_post_variable, oneBramsGrid, onePostGrid, &
               oneNamelistFile, oneBasicFields, oneTurbFields)
       case default
          write(str(1),"(i8)") one_post_variable%ivar_type
          call fatal_error(h//" unknown ivar_type="//trim(adjustl(str(1))))
       end select
    end if

  end subroutine PostOneField

end module ModPostOneField
