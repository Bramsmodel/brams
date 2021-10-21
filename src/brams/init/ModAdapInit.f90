!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModAdapInit

  implicit none
  
  private

  public :: lpuvw_init

contains

  subroutine lpuvw_init(n2,n3,lpu_R,lpv_R,lpw_R)
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    real, intent(out) :: lpu_R(:,:) ! (n2,n3)
    real, intent(out) :: lpv_R(:,:) ! (n2,n3)
    real, intent(out) :: lpw_R(:,:) ! (n2,n3)

    integer :: i,j

    do j = 1,n3
       do i = 1,n2
          lpu_R(i,j) = 2.0
          lpv_R(i,j) = 2.0
          lpw_R(i,j) = 2.0
       enddo
    enddo
  end subroutine lpuvw_init

end module ModAdapInit
