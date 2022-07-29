!############################# Change Log ##################################
! 3.1.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003, 2004 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################

module mem_stilt

  use grid_dims, only: &
       maxgrds

  use io_params, only: frqanl

  !--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  use ModNamelistFile, only: namelistFile

  include "constants.h"
  !--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  type stilt_vars
     real, pointer, dimension(:,:,:) :: &
          thvlast,lnthvadv, lnthetav, lnthvtend, afxu, afxv, afxw &
          ,ltscale, sigw, ltscaleb, sigwb &
          ,tkepb, afxub, afxvb, afxwb &
          ,cfxup1, cfxdn1, dfxup1, efxup1 &
          ,dfxdn1, efxdn1, cfxup2, dfxup2 &
          ,efxup2                         & 
          !-srf : for the true air density
          ,dnp

     real, pointer, dimension(:,:) :: pblhgt,lmo
  end type stilt_vars

  type (stilt_vars), allocatable :: stilt_g(:), stiltm_g(:)
  integer                        :: iexev,imassflx   

  real                          :: frqmassave
  !----- These variables control the time when the averages should be reset. -------------!
  real   , dimension(maxgrds)   :: etime_adve
  real   , dimension(maxgrds)   :: etime_turb

contains

  subroutine alloc_stilt(idiffk,stilt,n1,n2,n3,ng)
    implicit none
    type (stilt_vars) :: stilt
    integer, intent(in) :: n1, n2, n3, ng,idiffk

    ! Allocate arrays based on options (if necessary)
    if (iexev == 2) then
       allocate (stilt%thvlast(n1,n2,n3))
       allocate (stilt%lnthvadv(n1,n2,n3))
       allocate (stilt%lnthetav(n1,n2,n3))
       allocate (stilt%lnthvtend(n1,n2,n3))
       !-srf : for the true air density
       allocate (stilt%dnp(n1,n2,n3))
    endif

    if (idiffk == 7) then
       allocate (stilt%ltscale(n1,n2,n3))
       allocate (stilt%sigw(n1,n2,n3))
       allocate (stilt%pblhgt(n2,n3))
       !
       allocate (stilt%lmo(n2,n3))
    end if

    if (imassflx == 1) then
       if (idiffk == 7) then
          allocate (stilt%ltscale(n1,n2,n3))
          allocate (stilt%sigw(n1,n2,n3))
          allocate (stilt%pblhgt(n2,n3))
       end if
       allocate (stilt%afxu(n1,n2,n3))
       allocate (stilt%afxv(n1,n2,n3))
       allocate (stilt%afxw(n1,n2,n3))
       allocate (stilt%afxub(n1,n2,n3))
       allocate (stilt%afxvb(n1,n2,n3))
       allocate (stilt%afxwb(n1,n2,n3))
       allocate (stilt%cfxup1(n1,n2,n3))
       allocate (stilt%cfxdn1(n1,n2,n3))
       allocate (stilt%dfxup1(n1,n2,n3))
       allocate (stilt%efxup1(n1,n2,n3))
       allocate (stilt%dfxdn1(n1,n2,n3))
       allocate (stilt%efxdn1(n1,n2,n3))
       allocate (stilt%cfxup2(n1,n2,n3))
       allocate (stilt%dfxup2(n1,n2,n3))
       allocate (stilt%efxup2(n1,n2,n3))
       if (idiffk /= 2 .and. idiffk /= 3) then
          allocate (stilt%tkepb(n1,n2,n3))
       end if
       if (idiffk == 7) then
          allocate (stilt%ltscaleb(n1,n2,n3))
          allocate (stilt%sigwb(n1,n2,n3))
       end if
    end if

    !srf- initialization
    frqmassave =frqanl
    etime_adve =0.0
    etime_turb =0.0

  end subroutine alloc_stilt


  subroutine nullify_stilt(stilt)
    implicit none
    type (stilt_vars)   :: stilt

    if (associated(stilt%thvlast ))  nullify (stilt%thvlast )
    if (associated(stilt%lnthvadv ))  nullify (stilt%lnthvadv )
    if (associated(stilt%lnthetav ))  nullify (stilt%lnthetav )
    if (associated(stilt%lnthvtend ))  nullify (stilt%lnthvtend )
    !-srf : for the true air density
    if (associated(stilt%dnp     ))  nullify (stilt%dnp     )

    if (associated(stilt%ltscale ))  nullify (stilt%ltscale )
    if (associated(stilt%sigw    ))  nullify (stilt%sigw    )
    if (associated(stilt%pblhgt  ))  nullify (stilt%pblhgt  )
    if (associated(stilt%lmo     ))  nullify (stilt%lmo     )
    if (associated(stilt%afxu    ))  nullify (stilt%afxu    )
    if (associated(stilt%afxv    ))  nullify (stilt%afxv    )
    if (associated(stilt%afxw    ))  nullify (stilt%afxw    )
    if (associated(stilt%ltscaleb))  nullify (stilt%ltscaleb)
    if (associated(stilt%sigwb   ))  nullify (stilt%sigwb   )
    if (associated(stilt%tkepb   ))  nullify (stilt%tkepb   ) 
    if (associated(stilt%afxub   ))  nullify (stilt%afxub   )
    if (associated(stilt%afxvb   ))  nullify (stilt%afxvb   )
    if (associated(stilt%afxwb   ))  nullify (stilt%afxwb   )
    if (associated(stilt%cfxup1  ))  nullify (stilt%cfxup1  )
    if (associated(stilt%cfxdn1  ))  nullify (stilt%cfxdn1  )
    if (associated(stilt%dfxup1  ))  nullify (stilt%dfxup1  )
    if (associated(stilt%efxup1  ))  nullify (stilt%efxup1  )
    if (associated(stilt%dfxdn1  ))  nullify (stilt%dfxdn1  )
    if (associated(stilt%efxdn1  ))  nullify (stilt%efxdn1  )
    if (associated(stilt%cfxup2  ))  nullify (stilt%cfxup2  )
    if (associated(stilt%dfxup2  ))  nullify (stilt%dfxup2  )
    if (associated(stilt%efxup2  ))  nullify (stilt%efxup2  )
    return
  end subroutine nullify_stilt

  subroutine dealloc_stilt(stilt)
    implicit none
    type (stilt_vars)   :: stilt

    if (associated(stilt%thvlast ))  deallocate (stilt%thvlast )
    if (associated(stilt%lnthvadv ))  deallocate (stilt%lnthvadv )
    if (associated(stilt%lnthetav ))  deallocate (stilt%lnthetav ) 
    if (associated(stilt%lnthvtend ))  deallocate (stilt%lnthvtend )
    !-srf : for the true air density
    if (associated(stilt%dnp     ))  deallocate (stilt%dnp     )

    if (associated(stilt%ltscale ))  deallocate (stilt%ltscale )
    if (associated(stilt%sigw    ))  deallocate (stilt%sigw    )
    if (associated(stilt%pblhgt  ))  deallocate (stilt%pblhgt  )
    if (associated(stilt%lmo     ))  deallocate (stilt%lmo     )
    if (associated(stilt%afxu    ))  deallocate (stilt%afxu    )
    if (associated(stilt%afxv    ))  deallocate (stilt%afxv    )
    if (associated(stilt%afxw    ))  deallocate (stilt%afxw    )
    if (associated(stilt%ltscaleb))  deallocate (stilt%ltscaleb)
    if (associated(stilt%sigwb   ))  deallocate (stilt%sigwb   )
    if (associated(stilt%tkepb   ))  deallocate (stilt%tkepb   )
    if (associated(stilt%afxub   ))  deallocate (stilt%afxub   )
    if (associated(stilt%afxvb   ))  deallocate (stilt%afxvb   )
    if (associated(stilt%afxwb   ))  deallocate (stilt%afxwb   )
    if (associated(stilt%cfxup1  ))  deallocate (stilt%cfxup1  )
    if (associated(stilt%cfxdn1  ))  deallocate (stilt%cfxdn1  )
    if (associated(stilt%dfxup1  ))  deallocate (stilt%dfxup1  )
    if (associated(stilt%efxup1  ))  deallocate (stilt%efxup1  )
    if (associated(stilt%dfxdn1  ))  deallocate (stilt%dfxdn1  )
    if (associated(stilt%efxdn1  ))  deallocate (stilt%efxdn1  )
    if (associated(stilt%cfxup2  ))  deallocate (stilt%cfxup2  )
    if (associated(stilt%dfxup2  ))  deallocate (stilt%dfxup2  )
    if (associated(stilt%efxup2  ))  deallocate (stilt%efxup2  )
    return
  end subroutine dealloc_stilt

  subroutine filltab_stilt(oneVarTable, oneVarTableSize, stilt, stiltm)

    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type (stilt_vars) :: stilt
    type (stilt_vars), optional :: stiltm
    
    if (associated(stilt%thvlast)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%thvlast, &
               'THVLAST :3:hist:mpti:mpt3:mpt1', &
               stiltm%thvlast)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%thvlast, &
               'THVLAST :3:hist:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%lnthvadv)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%lnthvadv, &
               'LNTHVADV :3:hist:mpti:mpt3:mpt1', &
               stiltm%lnthvadv)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%lnthvadv, &
               'LNTHVADV :3:hist:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%lnthetav)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%lnthetav, &
               'LNTHETAV :3:hist:mpti:mpt3:mpt1', &
               stiltm%lnthetav)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%lnthetav, &
               'LNTHETAV :3:hist:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%lnthvtend)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%lnthvtend, &
               'LNTHVTEND :3:hist:mpti:mpt3:mpt1', &
               stiltm%lnthvtend)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%lnthvtend, &
               'LNTHVTEND :3:hist:mpti:mpt3:mpt1')
       end if
    end if

    !-srf : for the true air density
    if (associated(stilt%dnp)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%dnp, &
               'DNP :3:hist:anal:mpti:mpt3', &
               stiltm%dnp)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%dnp, &
               'DNP :3:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(stilt%ltscale)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%ltscale, &
               'TL :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%ltscale)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%ltscale, &
               'TL :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
    
    if (associated(stilt%sigw)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%sigw, &
               'SIGW :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%sigw)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%sigw, &
               'SIGW :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
    
    if (associated(stilt%afxu)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxu, &
               'AFXU :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%afxu)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxu, &
               'AFXU :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
       
    if (associated(stilt%afxv)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxv, &
               'AFXV :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%afxv)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxv, &
               'AFXV :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
    
    if (associated(stilt%afxw)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxw, &
               'AFXW :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%afxw)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxw, &
               'AFXW :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
    
    if (associated(stilt%ltscaleb)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%ltscaleb, &
               'TLB :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%ltscaleb)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%ltscaleb, &
               'TLB :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
    
    if (associated(stilt%sigwb)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%sigwb, &
               'SIGWB :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%sigwb)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%sigwb, &
               'SIGWB :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%tkepb)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%tkepb, &
               'TKEPB :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%tkepb)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%tkepb, &
               'TKEPB :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
       
    if (associated(stilt%afxub)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxub, &
               'AFXUB :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%afxub)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxub, &
               'AFXUB :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
    
    if (associated(stilt%afxvb)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxvb, &
               'AFXVB :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%afxvb)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxvb, &
               'AFXVB :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
    
    if (associated(stilt%afxwb)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxwb, &
               'AFXWB :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%afxwb)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%afxwb, &
               'AFXWB :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
       
    if (associated(stilt%cfxup1)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%cfxup1, &
               'CFXUP1 :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%cfxup1)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%cfxup1, &
               'CFXUP1 :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%cfxdn1)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%cfxdn1, &
               'CFXDN1 :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%cfxdn1)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%cfxdn1, &
               'CFXDN1 :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%dfxup1)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%dfxup1, &
               'DFXUP1 :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%dfxup1)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%dfxup1, &
               'DFXUP1 :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%efxup1)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%efxup1, &
               'EFXUP1 :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%efxup1)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%efxup1, &
               'EFXUP1 :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%dfxdn1)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%dfxdn1, &
               'DFXDN1 :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%dfxdn1)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%dfxdn1, &
               'DFXDN1 :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%efxdn1)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%efxdn1, &
               'EFXDN1 :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%efxdn1)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%efxdn1, &
               'EFXDN1 :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%cfxup2)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%cfxup2, &
               'CFXUP2 :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%cfxup2)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%cfxup2, &
               'CFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%dfxup2)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%dfxup2, &
               'DFXUP2 :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%dfxup2)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%dfxup2, &
               'DFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%efxup2)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%efxup2, &
               'EFXUP2 :3:hist:anal:mpti:mpt3:mpt1', &
               stiltm%efxup2)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%efxup2, &
               'EFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%pblhgt)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%pblhgt, &
               'PBLHGT :2:hist:anal:mpti:mpt3:mpt1', &
               stiltm%pblhgt)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%pblhgt, &
               'PBLHGT :2:hist:anal:mpti:mpt3:mpt1')
       end if
    end if

    if (associated(stilt%lmo)) then
       if (present(stiltm)) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%lmo, &
               'LMO    :2:hist:anal:mpti:mpt3:mpt1', &
               stiltm%lmo)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               stilt%lmo, &
               'LMO    :2:hist:anal:mpti:mpt3:mpt1')
       end if
    end if
  end subroutine filltab_stilt

  !=======================================================================================!
  !=======================================================================================!
  !     This function finds the virtual temperature based on the temperature and mixing   !
  ! ratio. Two notes: 
  ! 1. It will use the condensation in case the total mixing ratio is provided.           !
  ! 2. This can be used for virtual potential temperature, just give potential tempera-   !
  !    ture instead of temperature. 
  !---------------------------------------------------------------------------------------!
  real function virtt(temp,rvap,rtot)

    use rconstants, only: ep

    implicit none

    !----- Arguments--------------------------------------------------------------------!
    real, intent(in)           :: temp     ! Temperature         [    K]
    real, intent(in)           :: rvap     ! Vapour mixing ratio [kg/kg]
    real, intent(in)           :: rtot     ! Total mixing ratio  [kg/kg]
    !------------------------------------------------------------------------------------!

    virtt = temp * (1. + rvap / ep) / (1. + rtot)

    return

  end function virtt
  !=======================================================================================!
  !=======================================================================================!
  subroutine zero_average_mass_adve(stilt)
    implicit none
    type (stilt_vars)   :: stilt

    if (associated(stilt%afxub   ))  stilt%afxub   = 0.0
    if (associated(stilt%afxvb   ))  stilt%afxvb   = 0.0
    if (associated(stilt%afxwb   ))  stilt%afxwb   = 0.0

    return
  end subroutine zero_average_mass_adve
  !=======================================================================================!  
  !==========================================================================================!
  !    This subroutine flushes the average variables to zero whenever necessary. This is     !
  ! done at the output subroutine, right after the analysis (lite and full) are saved.       !
  !------------------------------------------------------------------------------------------!
  subroutine zero_average_mass_turb(stilt)
    implicit none
    type (stilt_vars)   :: stilt

    if (associated(stilt%ltscaleb))  stilt%ltscaleb= 0.0
    if (associated(stilt%sigwb   ))  stilt%sigwb   = 0.0
    if (associated(stilt%tkepb   ))  stilt%tkepb   = 0.0

    return
  end subroutine zero_average_mass_turb

  !--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  subroutine StoreNamelistFileAtMem_stilt(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    iexev = oneNamelistFile%iexev
    imassflx = oneNamelistFile%imassflx
  end subroutine StoreNamelistFileAtMem_stilt
  !--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

end module mem_stilt
