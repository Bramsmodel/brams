! Module necessary to GRELL param.

module mem_grell

  type grell_vars

     ! Variables to be dimensioned by (m2,m3)
     real, pointer, contiguous :: UPMF(:,:)
     real, pointer, contiguous :: DNMF(:,:)
     real, pointer, contiguous :: XIACT_C(:,:)
     real, pointer, contiguous :: XIACT_P(:,:)
     real, pointer, contiguous :: XIERR(:,:)
     real, pointer, contiguous :: XKDT(:,:)
     real, pointer, contiguous :: XKTOP(:,:)
     real, pointer, contiguous :: XKBCON(:,:)
     real, pointer, contiguous :: XJMIN(:,:)
     real, pointer, contiguous :: XK22(:,:)

     ! Variables to be dimensioned by (m1,m2,m3)
     !real, pointer, dimension(:,:,:) :: &
     !lsfth, &
     !lsfrt

  end type grell_vars

  type (grell_vars), allocatable, target :: grell_g(:)
  type (grell_vars), allocatable, target :: grellm_g(:)
  type (grell_vars), allocatable, target :: grell_g_sh(:)
  type (grell_vars), allocatable, target :: grellm_g_sh(:)


  type cuforc_vars
     real, pointer, contiguous :: lsfth(:,:,:)
     real, pointer, contiguous :: lsfrt(:,:,:)
  end type cuforc_vars

  type (cuforc_vars), allocatable, target :: cuforc_g(:)
  type (cuforc_vars), allocatable, target :: cuforcm_g(:)
  type (cuforc_vars), allocatable, target :: cuforc_sh_g(:)
  type (cuforc_vars), allocatable, target :: cuforcm_sh_g(:)



contains

  subroutine alloc_grell(grell, m1, m2, m3, ng)

    use mem_cuparm, only : nnqparm  ! INTENT(IN)

    implicit none
    ! Arguments
    type (grell_vars), intent(INOUT) :: grell
    integer, intent(IN)              :: m1, m2, m3, ng
    ! Local Variables:
    integer :: ierr
    character(len=*), parameter :: h="**(alloc_grell)**"

    ! Allocate arrays based on options (if necessary)
    if (abs(nnqparm(ng))==2)  then
       allocate (grell%UPMF     (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%UPMF")
       allocate (grell%DNMF     (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%DNMF")
       allocate (grell%XIACT_C  (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XIACT_C")
       allocate (grell%XIACT_P  (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XIACT_P")
       allocate (grell%XIERR    (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XIERR")
       allocate (grell%XKDT     (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XKDT")
       allocate (grell%XKTOP    (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XKTOP")
       allocate (grell%XKBCON   (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XKBCON")
       allocate (grell%XJMIN    (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XJMIN")
       allocate (grell%XK22     (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XK22")


       !--(DMK-LFR NEC-SX6)----------------------------------------------
       grell%upmf = 0.
       grell%dnmf = 0.
       grell%xiact_c = 0.
       grell%xiact_p = 0.
       grell%xierr = 0.
       grell%xkdt = 0.
       grell%xktop = 0.
       grell%xkbcon = 0.
       grell%xjmin = 0.
       grell%xk22 = 0.
       !--(DMK-LFR NEC-SX6)----------------------------------------------

    endif

  end subroutine alloc_grell
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine alloc_grell_sh(grell, m1, m2, m3, ng)

    use mem_cuparm, only : nnqparm  ! INTENT(IN)
    use shcu_vars_const, only : nnshcu   ! INTENT(IN)

    implicit none
    type (grell_vars) :: grell
    integer, intent(in) :: m1, m2, m3, ng

    ! Allocate arrays based on options (if necessary)

    allocate (grell%UPMF     (m2, m3));grell%UPMF    =0.0
    allocate (grell%DNMF     (m2, m3));grell%DNMF    =0.0
    allocate (grell%XIACT_C  (m2, m3));grell%XIACT_C =0.0
    allocate (grell%XIACT_P  (m2, m3));grell%XIACT_P =0.0
    allocate (grell%XIERR    (m2, m3));grell%XIERR   =0.0
    allocate (grell%XKDT     (m2, m3));grell%XKDT    =0.0
    allocate (grell%XKTOP    (m2, m3));grell%XKTOP   =0.0
    allocate (grell%XKBCON   (m2, m3));grell%XKBCON  =0.0
    allocate (grell%XJMIN    (m2, m3));grell%XJMIN   =0.0
    allocate (grell%XK22     (m2, m3));grell%XK22=0.0

    return
  end subroutine alloc_grell_sh
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine alloc_cu_forcings(cuforc, m1, m2, m3, ng)

    !    use mem_cuparm, only : nnqparm 
    !    use shcu_vars_const, only : nnshcu   

    implicit none
    type (cuforc_vars) :: cuforc
    integer, intent(in) :: m1, m2, m3, ng

    allocate (cuforc%lsfth(m1, m2, m3));cuforc%lsfth=0.0
    allocate (cuforc%lsfrt(m1, m2, m3));cuforc%lsfrt=0.0

  end subroutine alloc_cu_forcings
  !---------------------------------------------------------------
  subroutine nullify_cuforc(cuforc)

    implicit none
    type (cuforc_vars) :: cuforc

    if (associated(cuforc%lsfth))   nullify (cuforc%lsfth)
    if (associated(cuforc%lsfrt))   nullify (cuforc%lsfrt)

  end subroutine nullify_cuforc
  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine filltab_cuforc_sh(oneVarTable, oneVarTableSize, &
       cuforc, cuforcm)

    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable
    
    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(cuforc_vars), pointer, intent(in) :: cuforc
    type(cuforc_vars), pointer, intent(in) :: cuforcm

    logical :: assThis
    character(len=*), parameter :: h="**(filltab_cuforc_sh)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(cuforc)) then
       call fatal_error(h//" cuforc not associated")
    end if
    
    if (associated(cuforc%lsfth)) then
       if (.not. associated(cuforcm)) then
          assThis=.false.
       else
          assThis=associated(cuforcm%lsfth)
       end if
       if (assThis) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuforc%lsfth, &
               'LSFTH_SH :3:hist:anal:mpti:mpt3', &
               cuforcm%lsfth)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuforc%lsfth, &
               'LSFTH_SH :3:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(cuforc%lsfrt)) then
       if (.not. associated(cuforcm)) then
          assThis=.false.
       else
          assThis=associated(cuforcm%lsfrt)
       end if
       if (assThis) then
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuforc%lsfrt, &
               'LSFRT_SH :3:hist:anal:mpti:mpt3', &
               cuforcm%lsfrt)
       else
          call InsertAtVarTable(oneVarTable, oneVarTableSize, &
               cuforc%lsfrt, &
               'LSFRT_SH :3:hist:anal:mpti:mpt3')
       end if
    end if
  end subroutine filltab_cuforc_sh

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine filltab_cuforc(oneVarTable, oneVarTableSize, &
       cuforc, cuforcm)

    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(cuforc_vars), pointer, intent(in) :: cuforc
    type(cuforc_vars), pointer, intent(in) :: cuforcm

    logical :: assThis
    character(len=*), parameter :: h="**(filltab_cuforc)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(cuforc)) then
       call fatal_error(h//" cuforc not associated")
    end if
    
    if (associated(cuforc%lsfth)) then
       if (.not. associated(cuforcm)) then
          assThis=.false.
       else
          assThis=associated(cuforcm%lsfth)
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               cuforc%lsfth, &
               'LSFTH :3:hist:anal:mpti:mpt3', &
               cuforcm%lsfth)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               cuforc%lsfth, &
               'LSFTH :3:hist:anal:mpti:mpt3')
       end if
    end if

    if (associated(cuforc%lsfrt)) then
       if (.not. associated(cuforcm)) then
          assThis=.false.
       else
          assThis=associated(cuforcm%lsfrt)
       end if
       if (assThis) then
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               cuforc%lsfrt, &
               'LSFRT :3:hist:anal:mpti:mpt3', &
               cuforcm%lsfrt)
       else
          call InsertAtVarTable (oneVarTable, oneVarTableSize, &
               cuforc%lsfrt, &
               'LSFRT :3:hist:anal:mpti:mpt3')
       end if
    end if
  end subroutine filltab_cuforc

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine nullify_grell(grell)

    implicit none
    type (grell_vars), intent(INOUT) :: grell

    nullify (grell%UPMF)
    nullify (grell%DNMF)
    nullify (grell%XIACT_C)
    nullify (grell%XIACT_P)
    nullify (grell%XIERR)
    nullify (grell%XKDT)
    nullify (grell%XKTOP)
    nullify (grell%XKBCON)
    nullify (grell%XJMIN)
    nullify (grell%XK22)


  end subroutine nullify_grell

  ! *************************************************************************

  subroutine dealloc_grell(grell)

    implicit none
    type (grell_vars), intent(INOUT) :: grell

    if (associated(grell%UPMF))    deallocate (grell%UPMF)
    if (associated(grell%DNMF))    deallocate (grell%DNMF)
    if (associated(grell%XIACT_C)) deallocate (grell%XIACT_C)
    if (associated(grell%XIACT_P)) deallocate (grell%XIACT_P)
    if (associated(grell%XIERR))   deallocate (grell%XIERR)
    if (associated(grell%XKDT))    deallocate (grell%XKDT)
    if (associated(grell%XKTOP))   deallocate (grell%XKTOP)
    if (associated(grell%XKBCON))  deallocate (grell%XKBCON)
    if (associated(grell%XJMIN))   deallocate (grell%XJMIN)
    if (associated(grell%XK22))    deallocate (grell%XK22)

  end subroutine dealloc_grell

  ! *************************************************************************

  subroutine filltab_grell(oneVarTable, oneVarTableSize, &
       grell, grellm, nnqparm)

    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable
    
    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(grell_vars), pointer, intent(in) :: grell
    type(grell_vars), pointer, intent(in) :: grellm
    integer, intent(in) :: nnqparm

    ! Local Variables:
    logical :: assAve
    logical :: assThis
    character(len=*), parameter :: h="**(filltab_grell)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(grell)) then
       call fatal_error(h//" grell not associated")
    end if
    
    ! Fill pointers to arrays into variable tables

    if (nnqparm == 2 )  then
       assAve=associated(grellm)

       if (assAve) then
          assThis=associated(grellm%UPMF)
       else
          assThis=.false.
       end if
       if (associated(grell%UPMF)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%UPMF, &
                  'UPMF :2:hist:anal:mpti:mpt3', &
                  grellm%UPMF)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%UPMF, &
                  'UPMF :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (assAve) then
          assThis=associated(grellm%DNMF)
       else
          assThis=.false.
       end if
       if (associated(grell%DNMF)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%DNMF, &
                  'DNMF :2:hist:anal:mpti:mpt3', &
                  grellm%DNMF)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%DNMF, &
                  'DNMF :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (assAve) then
          assThis=associated(grellm%XIACT_C)
       else
          assThis=.false.
       end if
       if (associated(grell%XIACT_C)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XIACT_C, &
                  'XIACT_C :2:hist:anal:mpti:mpt3', &
                  grellm%XIACT_C)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XIACT_C, &
                  'XIACT_C :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (assAve) then
          assThis=associated(grellm%XIACT_P)
       else
          assThis=.false.
       end if
       if (associated(grell%XIACT_P)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XIACT_P, &
                  'XIACT_P :2:hist:anal:mpti:mpt3', &
                  grellm%XIACT_P)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XIACT_P, &
                  'XIACT_P :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (assAve) then
          assThis=associated(grellm%XIERR)
       else
          assThis=.false.
       end if
       if (associated(grell%XIERR)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XIERR, &
                  'XIERR :2:hist:anal:mpti:mpt3', &
                  grellm%XIERR)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XIERR, &
                  'XIERR :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (assAve) then
          assThis=associated(grellm%XKDT)
       else
          assThis=.false.
       end if
       if (associated(grell%XKDT)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XKDT, &
                  'XKDT :2:hist:anal:mpti:mpt3', &
                  grellm%XKDT)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XKDT, &
                  'XKDT :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (assAve) then
          assThis=associated(grellm%XKTOP)
       else
          assThis=.false.
       end if
       if (associated(grell%XKTOP)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XKTOP, &
                  'XKTOP :2:hist:anal:mpti:mpt3', &
                  grellm%XKTOP)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XKTOP, &
                  'XKTOP :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (assAve) then
          assThis=associated(grellm%XKBCON)
       else
          assThis=.false.
       end if
       if (associated(grell%XKBCON)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XKBCON, &
                  'XKBCON :2:hist:anal:mpti:mpt3', &
                  grellm%XKBCON)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XKBCON, &
                  'XKBCON :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (assAve) then
          assThis=associated(grellm%XJMIN)
       else
          assThis=.false.
       end if
       if (associated(grell%XJMIN)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XJMIN, &
                  'XJMIN :2:hist:anal:mpti:mpt3', &
                  grellm%XJMIN)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XJMIN, &
                  'XJMIN :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (assAve) then
          assThis=associated(grellm%XK22)
       else
          assThis=.false.
       end if
       if (associated(grell%XK22)) then
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XK22, &
                  'XK22 :2:hist:anal:mpti:mpt3', &
                  grellm%XK22)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell%XK22, &
                  'XK22 :2:hist:anal:mpti:mpt3')
          end if
       end if
    endif

  end subroutine filltab_grell

  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine filltab_grell_sh(oneVarTable, oneVarTableSize, &
       grell_sh, grellm_sh, nnshcu)

    use ModVarTable, only: &
         VarTable, &
         InsertAtVarTable
    
    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(grell_vars), pointer, intent(in) :: grell_sh
    type(grell_vars), pointer, intent(in) :: grellm_sh
    integer, intent(in) :: nnshcu

    ! Local Variables:
    logical :: assAve
    logical :: assThis
    character(len=*), parameter :: h="**(filltab_grell_sh)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    else if (.not. associated(grell_sh)) then
       call fatal_error(h//" grell_sh not associated")
    end if
    
    ! Fill pointers to arrays into variable tables

    if (nnshcu==1 .or. nnshcu==2) then
       assAve=associated(grellm_sh)
       
       if (associated(grell_sh%UPMF)) then
          if (assAve) then
             assThis=associated(grellm_sh%UPMF)
          else
             assThis=.false.
          end if
          if (assThis) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  grell_sh%UPMF, &
                  'UPMFSH :2:hist:anal:mpti:mpt3', &
                  grellm_sh%UPMF)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  grell_sh%UPMF, &
                  'UPMFSH :2:hist:anal:mpti:mpt3')
          end if
       end if
       
       if (associated(grell_sh%XIERR)) then
          if (assAve) then
             assThis=associated(grellm_sh%XIERR)
          else
             assThis=.false.
          end if
          if (assThis) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  grell_sh%XIERR, &
                  'XIERRSH :2:hist:anal:mpti:mpt3', &
                  grellm_sh%XIERR)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  grell_sh%XIERR, &
                  'XIERRSH :2:hist:anal:mpti:mpt3')
          end if
       end if
       
       if (associated(grell_sh%XKTOP)) then
          if (assAve) then
             assThis=associated(grellm_sh%XKTOP)
          else
             assThis=.false.
          end if
          if (assThis) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  grell_sh%XKTOP, &
                  'XKTOPSH :2:hist:anal:mpti:mpt3', &
                  grellm_sh%XKTOP)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  grell_sh%XKTOP, &
                  'XKTOPSH :2:hist:anal:mpti:mpt3')
          end if
       end if

       if (associated(grell_sh%XKBCON)) then
          if (assAve) then
             assThis=associated(grellm_sh%XKBCON)
          else
             assThis=.false.
          end if
          if (assThis) then
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  grell_sh%XKBCON, &
                  'XKBCONSH :2:hist:mpti:mpt3', &
                  grellm_sh%XKBCON)
          else
             call InsertAtVarTable (oneVarTable, oneVarTableSize, &
                  grell_sh%XKBCON, &
                  'XKBCONSH :2:hist:mpti:mpt3')
          end if
       end if

       if (associated(grell_sh%XK22)) then
          if (assAve) then
             assThis=associated(grellm_sh%XK22)
          else
             assThis=.false.
          end if
          if (assThis) then
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell_sh%XK22, &
                  'XK22SH :2:hist:mpti:mpt3', &
                  grellm_sh%XK22)
          else
             call InsertAtVarTable(oneVarTable, oneVarTableSize, &
                  grell_sh%XK22, &
                  'XK22SH :2:hist:mpti:mpt3')
          end if
       end if

    endif
  end subroutine filltab_grell_sh

end module mem_grell
