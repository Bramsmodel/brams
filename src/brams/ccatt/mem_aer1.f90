!###########################################################################
!  B - Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_aer1

  use grid_dims, only : &
       maxgrds, &
       maxsclr

  use ModNamelistFile, only: &
       namelistFile

  use ModScalarTable, only: &
       ScalarTable, &
       InsertAtScalarTab

  use aer1_list, only : &
       nspecies_aer=> nspecies, &
       nmodes, spc_alloc, spc_name, src, ddp, wdp, fdda, offline, on ,off, &
       mode_alloc, transport, aer_name, mode_name, numb_mod_alloc, numb_alloc, &
       numb_name, inorg_spc_name, inorg_alloc, inorg_mod_alloc 


  use ModVarTable, only: &
       VarTable, &
       InsertVarTable

  use io_params, only : ioutput         ! INTENT(IN)

  use mem_chem1, only: chem_assim,RECYCLE_TRACERS,&
       nspecies_transported ! this is first calculated at chemistry
  ! "filltab_tend_chem1" routine

  use chem1_list, only: &
       chem_name=>spc_name,    &
       chem_alloc=>spc_alloc,  & 
       chem_on=>on,            &
       chem_fdda=>fdda,        &
       chem_nspecies=>nspecies

  include "constants.h"

  type aer1_vars
     !--- All families
     real, pointer, contiguous :: sc_p(:,:,:)
     real, pointer, contiguous :: sc_pp(:,:,:)
     real, pointer, contiguous :: sc_pf(:,:,:)
     real, pointer, contiguous :: sc_src(:,:,:)
     real, pointer, contiguous :: sc_dd(:,:)
     real, pointer, contiguous :: sc_wd(:,:)
     real, pointer, contiguous :: sc_t(:)
     !-----------
  end type aer1_vars

  type (aer1_vars), pointer :: aer1_g(:,:,:) => null() !mass
  type (aer1_vars), pointer :: aer1m_g(:,:,:) => null() !mass
  type (aer1_vars), pointer :: aer2_g(:,:) => null() !number
  type (aer1_vars), pointer :: aer2m_g(:,:) => null()  !number
  type (aer1_vars), pointer :: aer1_inorg_g(:,:) => null() !mass inorg
  type (aer1_vars), pointer :: aer1m_inorg_g(:,:) => null() !mass inorg


  !- dimension of sources arrays (=1 for 2dim, =m1 for 3dim)
  integer :: aer1_src_z_dim_g(nspecies_aer*nmodes,maxgrds)
  integer :: aer2_src_z_dim_g(nmodes,maxgrds)
  integer :: AER_ASSIM



  !  character(LEN=20),dimension(nsrc),parameter :: src_name= &
  !! '12345678901234567890'
  !(/                      &
  !  'antro               '&
  !, 'bburn               '&
  !, 'bioge               '/)!

  !  integer     :: &
  !   antro   = 01, & ! anthropogenic sources
  !   bburn   = 02, & ! biomass burning sources
  !   bioge   = 03    ! biogenic sources ! must be equal to "nsrc"

  integer :: AEROSOL
  integer :: ind_mode_sedim(maxsclr+nspecies_aer*nmodes)
  integer :: ind_aer_sedim (maxsclr+nspecies_aer*nmodes)
  integer :: num_scalar_aer_1st
  real    :: AER_TIMESTEP         ! timestep of matrixiaintegration
contains
  !---------------------------------------------------------------

  subroutine alloc_aer1(aer1,nvert_src,n1,n2,n3,nmodes,nspecies)

    implicit none

    integer,intent(in) :: n1,n2,n3,nspecies,nmodes
    integer :: ispc,isrc,imode

    type (aer1_vars)    ,dimension(     nmodes,nspecies) :: aer1
    !    type (aer1_src_vars),dimension(nsrc,nmodes,nspecies) :: aer1_src
    integer,dimension(nspecies)    :: nvert_src

    !print*,'----------------------------------------------------------------'
    !print*,' memory allocation for aerosol species:'

    do ispc=1,nspecies

       do imode=1,nmodes

          !- this control is only for AEROSOL MODEL = 1
          IF(AEROSOL == 1 ) THEN
             !- for aerosols all allocated must be transported
             if(spc_alloc(transport,imode,ispc) /= 1 .and. mode_alloc(imode,ispc) == 1) then
                print*,"================================================"
                print*,'aerosol=',ispc,'mode=', imode
                print*,'all allocated must be transported, change aer1-list'
                print*,'transport =',spc_alloc(transport,imode,ispc)
                print*,'allocation=',mode_alloc(imode,ispc),aer_name(imode,ispc)
                print*,"================================================"
                stop 3334
             endif
          ENDIF

          !1st test: if the mode does not exist, cycle
          if(mode_alloc(imode,ispc) /= 1 ) cycle
          !
          !- allocate memory for the past time tracer mixing ratio
          allocate (aer1(imode,ispc)%sc_p  (n1,n2,n3)) ;aer1(imode,ispc)%sc_p = 0.


          !- allocate memory for the dry deposition of the tracer mixing ratio
          if(spc_alloc(ddp,imode,ispc) == on) then
             allocate (aer1(imode,ispc)%sc_dd    (n2,n3)) ;aer1(imode,ispc)%sc_dd= 0.
          endif

          !- allocate memory for the wet deposition of the tracer mixing ratio
          if(spc_alloc(wdp,imode,ispc) == on) then
             allocate (aer1(imode,ispc)%sc_wd    (n2,n3)) ;aer1(imode,ispc)%sc_wd= 0.
          endif

          !- allocate memory for the past/future time tracer mixing ratio
          !- (4D data assimilation)
          if(chem_assim == on) then
             if(spc_alloc(fdda,imode,ispc) == on) then
                allocate (aer1(imode,ispc)%sc_pp (n1,n2,n3),aer1(imode,ispc)%sc_pf (n1,n2,n3))
                aer1(imode,ispc)%sc_pp =0.; aer1(imode,ispc)%sc_pf = 0.
             endif
          endif

          !- allocate memory for the sources (3d=bburn and 2d= urban/bioge/marin/sdust)
          if(spc_alloc(src,imode,ispc) == on) then
             allocate (aer1(imode,ispc)%sc_src  (nvert_src(ispc),n2,n3))
             aer1(imode,ispc)%sc_src = 0.
          endif

       enddo
    enddo
    return
  end subroutine alloc_aer1

  !--------------------------------------------------------------------------

  subroutine dealloc_aer1(aer1,nmodes,nspecies)

    implicit none

    integer,intent(in) :: nspecies,nmodes
    type (aer1_vars)    ,dimension(nmodes,nspecies) :: aer1
    integer :: ispc,imode

    !  Deallocate arrays
    do ispc=1,nspecies
       do imode=1,nmodes

          if (associated(aer1(imode,ispc)%sc_p ) )    deallocate(aer1(imode,ispc)%sc_p  )
          if (associated(aer1(imode,ispc)%sc_src))    deallocate(aer1(imode,ispc)%sc_src)
          if (associated(aer1(imode,ispc)%sc_dd) )    deallocate(aer1(imode,ispc)%sc_dd )
          if (associated(aer1(imode,ispc)%sc_wd) )    deallocate(aer1(imode,ispc)%sc_wd )
          if (associated(aer1(imode,ispc)%sc_pp) )    deallocate(aer1(imode,ispc)%sc_pp )
          if (associated(aer1(imode,ispc)%sc_pf) )    deallocate(aer1(imode,ispc)%sc_pf )

       enddo
    enddo

    return
  end subroutine dealloc_aer1

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine nullify_aer1(aer1,nmodes,nspecies)

    implicit none

    integer,intent(in) ::nspecies,nmodes
    type (aer1_vars),dimension(nmodes,nspecies) :: aer1
    integer :: ispc,imode

    do ispc=1,nspecies
       do imode=1,nmodes

          if (associated(aer1(imode,ispc)%sc_p ) )    nullify (aer1(imode,ispc)%sc_p  )
          if (associated(aer1(imode,ispc)%sc_src))    nullify (aer1(imode,ispc)%sc_src)
          if (associated(aer1(imode,ispc)%sc_dd) )    nullify (aer1(imode,ispc)%sc_dd )
          if (associated(aer1(imode,ispc)%sc_wd) )    nullify (aer1(imode,ispc)%sc_wd )
          if (associated(aer1(imode,ispc)%sc_pp) )    nullify (aer1(imode,ispc)%sc_pp )
          if (associated(aer1(imode,ispc)%sc_pf) )    nullify (aer1(imode,ispc)%sc_pf )

       enddo
    enddo

    return
  end subroutine nullify_aer1

  !---------------------------------------------------------------
  subroutine filltab_aer1(oneVarTable, oneVarTableSize, &
       aer1, aer1m, nvert_src, imean)

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type (aer1_vars), pointer, intent(in) :: aer1(:,:)
    type (aer1_vars), pointer, intent(in) :: aer1m(:,:)
    integer, intent(in) :: nvert_src(:)
    integer, intent(in) :: imean

    integer :: ispc
    integer :: imode
    character(len=8) :: register_name
    character(len=8) :: str_recycle
    character(len=1) :: str_src_dim


    str_recycle = ''; str_src_dim = ''
    if (RECYCLE_TRACERS == 1 .or. RECYCLE_TRACERS == 2 .or. ioutput == 5) then
       str_recycle = ':recycle'
    endif

    !- Fill pointers to arrays into variable tables
    do ispc = 1, size(aer1,2)
       do imode = 1, size(aer1,1)

          if (aerosol==1) then
             register_name=trim(spc_name(ispc))//trim(mode_name(imode))
          else if (aerosol==2) then
             register_name=trim(aer_name(imode,ispc))
          end if

          if (associated(aer1(imode,ispc)%sc_p)) then

             !---- tracer mixing ratio (dimension 3d)
             call InsertVarTable(oneVarTable, oneVarTableSize, &
                  aer1(imode,ispc)%sc_p, &
                  trim(register_name)//'P :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle), &
                  aer1m(imode,ispc)%sc_p, imean)

             !---- sources (3 and 2 dimension)
             if(spc_alloc(src,imode,ispc) == on) then
                if (nvert_src(ispc) ==  1) then
                   str_src_dim = '2'  ! for 2d sources
                else
                   str_src_dim = '3'  ! for 3d sources
                end if

                call InsertVarTable(oneVarTable, oneVarTableSize, &
                     aer1(imode,ispc)%sc_src,                       &
                     trim(register_name)//'_SRC :'//trim(str_src_dim)//':hist:anal:mpti:mpt3:mpt1', &
                     aer1m(imode,ispc)%sc_src, imean)
             end if

             !---- dry and wet deposition (dimension 2d)
             if(spc_alloc(ddp,imode,ispc) == on) then
                call InsertVarTable(oneVarTable, oneVarTableSize, &
                     aer1(imode,ispc)%sc_dd, &
                     trim(register_name) //'DD :2:hist:anal:mpti:mpt3', &
                     aer1m(imode,ispc)%sc_dd, imean)
             end if

             if(spc_alloc(wdp,imode,ispc) == on) then
                call InsertVarTable(oneVarTable, oneVarTableSize, &
                     aer1(imode,ispc)%sc_wd, &
                     trim(register_name) //'WD :2:hist:anal:mpti:mpt3', &
                     aer1m(imode,ispc)%sc_wd, imean)
             end if
             !----  data assimilation (dimension 3d)
             if(chem_assim == on) then
                if(spc_alloc(fdda,imode,ispc) == on) then
                   call InsertVarTable(oneVarTable, oneVarTableSize, &
                        aer1(imode,ispc)%sc_pp, &
                        trim(register_name) //'PP :3:mpti', &
                        aer1m(imode,ispc)%sc_pp, imean)
                   call InsertVarTable(oneVarTable, oneVarTableSize, &
                        aer1(imode,ispc)%sc_pf, &
                        trim(register_name) //'PF :3:mpti', &
                        aer1m(imode,ispc)%sc_pf, imean)
                end if
             end if
          end if
       end do
    end do
  end subroutine filltab_aer1


  !---------------------------------------------------------------

  subroutine alloc_tend_aer1(nmzp,nmxp,nmyp,ngrs,nmodes,nspecies,proc_type)
    implicit none
    integer,intent(in)                   :: ngrs,proc_type,nspecies,nmodes
    integer,intent(in), dimension (ngrs) :: nmzp,nmxp,nmyp
    integer :: ng,ntpts,ispc,imode

    !         Find the maximum number of grid points needed for any grid.

    if(proc_type==1) then
       ntpts=1
    else
       ntpts=0
       do ng=1,ngrs
          ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
       enddo
    endif

!!!!!  WE ARE ONLY CHECKING GRID 1 !!!!!!!!!
!!!!!    All grids must have same scalars defined !!!!!!!
    do ispc=1,nspecies
       do imode=1,nmodes

          !1st test: if the mode does not exist, cycle
          if(mode_alloc(imode,ispc) /= 1 ) cycle

          if (associated(aer1_g(imode,ispc,1)%sc_p)) allocate (aer1_g(imode,ispc,1)%sc_t(ntpts))

          do ng=2,ngrs
             aer1_g(imode,ispc,ng)%sc_t => aer1_g(imode,ispc,1)%sc_t
          enddo

       enddo
    enddo

  end subroutine alloc_tend_aer1
  !---------------------------------------------------------------

  subroutine nullify_tend_aer1(nmodes,nspecies)

    implicit none
    integer,intent(in) :: nspecies,nmodes
    integer ::ispc,imode

    do ispc=1,nspecies
       do imode=1,nmodes
          if (associated(aer1_g(imode,ispc,1)%sc_t)) nullify (aer1_g(imode,ispc,1)%sc_t)
       enddo
    enddo


  end subroutine nullify_tend_aer1

  !---------------------------------------------------------------
  subroutine dealloc_tend_aer1(nmodes,nspecies)
    implicit none
    integer,intent(in) ::  nspecies,nmodes
    integer ::ispc,imode

    do ispc=1,nspecies
       do imode=1,nmodes
          if (associated(aer1_g(imode,ispc,1)%sc_t)) deallocate (aer1_g(imode,ispc,1)%sc_t)
       enddo
    enddo

  end  subroutine dealloc_tend_aer1

  !---------------------------------------------------------------

  subroutine filltab_tend_aer1(oneScalarTab, oneScalarTabSize, &
       nmodes, nspecies, ng)
    implicit none

    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(inout) :: oneScalarTabSize
    integer,intent(in) :: nmodes,nspecies,ng
    integer ::ispc,imode
    integer :: elements
    character(len=8) :: register_name

    num_scalar_aer_1st = 0 !- defines which num_scalar corresponds to the 1st aerosol

    do ispc=1,nspecies
       do imode=1,nmodes

          if(aerosol==1) register_name=trim(spc_name(ispc))//trim(mode_name(imode))
          if(aerosol==2) register_name=trim(aer_name(imode,ispc))

          ! Fill pointers to scalar arrays into scalar tables

          if (associated(aer1_g(imode,ispc,ng)%sc_t)) then
             elements = size(aer1_g(imode,ispc,ng)%sc_t)
             call InsertAtScalarTab(aer1_g(imode,ispc,ng)%sc_p, &
                  aer1_g(imode,ispc,ng)%sc_t, &
                  trim(register_name) //'P', &
                  oneScalarTab, oneScalarTabSize)
             !- total number of transported species (CHEM + AER)
             nspecies_transported = nspecies_transported + 1

             !-save the aerosol identity for sedimentation/advection routine
             if(num_scalar_aer_1st == 0) num_scalar_aer_1st = oneScalarTabSize

          endif

       enddo
    enddo

  end subroutine filltab_tend_aer1

  !--------------------------------------------------------------------------
  subroutine StoreNamelistFileAtMem_aer1(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    aerosol = oneNamelistFile%aerosol
    aer_assim = oneNamelistFile%aer_assim
    aer_timestep=oneNamelistFile%aer_timestep
  end subroutine StoreNamelistFileAtMem_aer1

  !- LFR: For matrix included aer2
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  subroutine alloc_aer2(aer2,n1,n2,n3,nmodes)

    implicit none

    integer,intent(in) :: n1,n2,n3,nmodes
    integer :: imode

    type (aer1_vars)    ,dimension(nmodes) :: aer2

    do imode=1,nmodes
       !1st test: if the mode does not exist, cycle
       if(numb_mod_alloc(imode) /= 1 ) cycle

       !- allocate memory for the past time tracer number concentration
       allocate (aer2(imode)%sc_p  (n1,n2,n3)) ;aer2(imode)%sc_p = 0.

       !- allocate memory for the dry deposition 
       if(numb_alloc(ddp,imode) == on) then
          allocate (aer2(imode)%sc_dd    (n2,n3)) ;aer2(imode)%sc_dd= 0.
       endif
       !- allocate memory for the wet deposition 
       if(numb_alloc(wdp,imode) == on) then
          allocate (aer2(imode)%sc_wd    (n2,n3)) ;aer2(imode)%sc_wd= 0.
       endif
       !- allocate memory for the past/future time
       !- (4D data assimilation)
       if(chem_assim == on) then
          if(numb_alloc(fdda,imode) == on) then
             allocate (aer2(imode)%sc_pp (n1,n2,n3),aer2(imode)%sc_pf (n1,n2,n3))
             aer2(imode)%sc_pp =0.; aer2(imode)%sc_pf = 0.
          endif
       endif
       !- allocate memory for the sources (all 3-d
       if(numb_alloc(src,imode) == on) then
          allocate (aer2(imode)%sc_src  (n1,n2,n3))
          aer2(imode)%sc_src = 0.
       endif
    enddo
  end subroutine alloc_aer2




  !--------------------------------------------------------------------------

  subroutine dealloc_aer2(aer2,nmodes)

    implicit none

    integer,intent(in) :: nmodes
    type (aer1_vars)    ,dimension(nmodes) :: aer2

  end subroutine dealloc_aer2

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine nullify_aer2(aer2,nmodes)

    implicit none

    integer,intent(in) :: nmodes
    type (aer1_vars),dimension(nmodes) :: aer2
    integer :: ispc

    do ispc=1,nmodes
       if (associated(aer2(ispc)%sc_p ) )    nullify (aer2(ispc)%sc_p  )
       if (associated(aer2(ispc)%sc_src))    nullify (aer2(ispc)%sc_src)
       if (associated(aer2(ispc)%sc_dd) )    nullify (aer2(ispc)%sc_dd )
       if (associated(aer2(ispc)%sc_wd) )    nullify (aer2(ispc)%sc_wd )
       if (associated(aer2(ispc)%sc_pp) )    nullify (aer2(ispc)%sc_pp )
       if (associated(aer2(ispc)%sc_pf) )    nullify (aer2(ispc)%sc_pf )
    enddo
  end subroutine nullify_aer2




  subroutine filltab_aer2(oneVarTable, oneVarTableSize, &
       aer2, aer2m, nvert_src, imean)

    use io_params, only : ioutput         ! INTENT(IN)

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(aer1_vars), pointer, intent(in) :: aer2(:)
    type(aer1_vars), pointer, intent(in) :: aer2m(:)
    integer, intent(in) :: nvert_src(:)
    integer, intent(in) :: imean

    integer :: imode,i

    character(len=8) :: str_recycle
    character(len=1) :: str_src_dim
    str_recycle = ''; str_src_dim = ''
    if (RECYCLE_TRACERS == 1 .or. ioutput == 5) then
       str_recycle = ':recycle'
    end if

    !- Fill pointers to arrays into variable tables
    do imode = 1, size(nvert_src)
       if (associated(aer2(imode)%sc_p)) then

          !---- number concentration (dimension 3d)
          call InsertVarTable(oneVarTable, oneVarTableSize, &
               aer2(imode)%sc_p, &
               trim(numb_name(imode))//'P :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle), &
               aer2m(imode)%sc_p, imean)

          if(numb_alloc(src,imode) == on) then
             !- for now we are using allways 3-d
             str_src_dim = '3' 
             call InsertVarTable(oneVarTable, oneVarTableSize, &
                  aer2 (imode)%sc_src, &
                  trim(numb_name(imode))//'_SRC :'//trim(str_src_dim)//':hist:anal:mpti:mpt3:mpt1', &
                  aer2m(imode)%sc_src, imean)
          endif

          !---- dry and wet deposition (dimension 2d)
          if(numb_alloc(ddp,imode) == on) then
             call InsertVarTable(oneVarTable, oneVarTableSize, &
                  aer2(imode)%sc_dd, &
                  trim(numb_name(imode))//'DD :2:hist:anal:mpti:mpt3', &
                  aer2m(imode)%sc_dd, imean)
          end if

          if(numb_alloc(wdp,imode) == on) then
             call InsertVarTable(oneVarTable, oneVarTableSize, &
                  aer2(imode)%sc_wd, &
                  trim(numb_name(imode))//'WD :2:hist:anal:mpti:mpt3', &
                  aer2m(imode)%sc_wd, imean)
          end if
          !----  data assimilation (dimension 3d)
          if(chem_assim == on) then
             if(numb_alloc(fdda,imode) == on) then
                call InsertVarTable(oneVarTable, oneVarTableSize, &
                     aer2(imode)%sc_pp, &
                     trim(numb_name(imode))//'PP :3:mpti', &
                     aer2m(imode)%sc_pp, imean)
                call InsertVarTable(oneVarTable, oneVarTableSize, &
                     aer2(imode)%sc_pf, &
                     trim(numb_name(imode))//'PF :3:mpti', &
                     aer2m(imode)%sc_pf, imean)
             end if
          end if
       end if
    end do
  end subroutine filltab_aer2





  !---------------------------------------------------------------

  subroutine alloc_tend_aer2(nmzp,nmxp,nmyp,ngrs,nmodes,proc_type)
    implicit none
    integer,intent(in)                   :: ngrs,proc_type,nmodes
    integer,intent(in), dimension (ngrs) :: nmzp,nmxp,nmyp
    integer :: ng,ntpts,ispc,imode
    if(proc_type==1) then
       ntpts=1
    else
       ntpts=0
       do ng=1,ngrs
          ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
       enddo
    endif
!!!!!  WE ARE ONLY CHECKING GRID 1 !!!!!!!!!
!!!!!    All grids must have same scalars defined !!!!!!!
    do imode=1,nmodes

       !1st test: if the mode does not exist, cycle
       if(numb_mod_alloc(imode) /= 1 ) cycle

       if (associated(aer2_g(imode,1)%sc_p)) allocate (aer2_g(imode,1)%sc_t(ntpts))

       do ng=2,ngrs
          aer2_g(imode,ng)%sc_t => aer2_g(imode,1)%sc_t
       enddo

    enddo
  end subroutine alloc_tend_aer2
  !---------------------------------------------------------------

  subroutine nullify_tend_aer2(nmodes)

    implicit none
    integer,intent(in) :: nmodes
    integer ::imode

    do imode=1,nmodes
       if (associated(aer2_g(imode,1)%sc_t)) nullify (aer2_g(imode,1)%sc_t)
    enddo

  end subroutine nullify_tend_aer2

  !---------------------------------------------------------------
  subroutine dealloc_tend_aer2(nmodes)
    implicit none
    integer,intent(in) ::  nmodes

  end  subroutine dealloc_tend_aer2

  !---------------------------------------------------------------

  subroutine filltab_tend_aer2(oneScalarTab, oneScalarTabSize, &
       nmodes,ng)
    use mem_chem1 , only: nspecies_transported ! this is first calculated at chemistry
    ! "filltab_tend_chem1" routine
    implicit none

    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(inout) :: oneScalarTabSize
    integer,intent(in) :: nmodes,ng
    integer :: imode
    integer :: elements
    ! Fill pointers to scalar arrays into scalar tables

    do imode=1,nmodes
       if (associated(aer2_g(imode,ng)%sc_t)) then
          elements = size(aer2_g(imode,ng)%sc_t)
          call InsertAtScalarTab(aer2_g(imode,ng)%sc_p, aer2_g(imode,ng)%sc_t, trim(numb_name(imode))//'P',&
               oneScalarTab, oneScalarTabSize)
          !- total number of transported species (CHEM + AER)
          nspecies_transported = nspecies_transported + 1

       endif
    enddo
  end subroutine filltab_tend_aer2
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  subroutine alloc_aer1_inorg(aer1,n1,n2,n3,ninorg)

    implicit none

    integer,intent(in) :: n1,n2,n3,ninorg
    integer :: ispc,isrc

    type (aer1_vars)    ,dimension(ninorg) :: aer1

    !print*,'----------------------------------------------------------------'
    !print*,' memory allocation for aerosol inorg species:'

    do ispc=1,ninorg

       !- allocate memory for the past time tracer mixing ratio
       allocate (aer1(ispc)%sc_p  (n1,n2,n3)) ;aer1(ispc)%sc_p = 0.


       !- allocate memory for the dry deposition of the tracer mixing ratio
       if(inorg_alloc(ddp,ispc) == on) then
          allocate (aer1(ispc)%sc_dd    (n2,n3)) ;aer1(ispc)%sc_dd= 0.
       endif

       !- allocate memory for the wet deposition of the tracer mixing ratio
       if(inorg_alloc(wdp,ispc) == on) then
          allocate (aer1(ispc)%sc_wd    (n2,n3)) ;aer1(ispc)%sc_wd= 0.
       endif

       !- allocate memory for the past/future time tracer mixing ratio
       !- (4D data assimilation)
       if(chem_assim == on) then
          if(inorg_alloc(fdda,ispc) == on) then
             allocate (aer1(ispc)%sc_pp (n1,n2,n3),aer1(ispc)%sc_pf (n1,n2,n3))
             aer1(ispc)%sc_pp =0.; aer1(ispc)%sc_pf = 0.
          endif
       endif

       !- allocate memory for the sources (3d=bburn and 2d= urban/bioge/marin/sdust)
       if(inorg_alloc(src,ispc) == on) then
          allocate (aer1(ispc)%sc_src  (n1,n2,n3))
          aer1(ispc)%sc_src = 0.
       endif

    enddo

    return
  end subroutine alloc_aer1_inorg

  !--------------------------------------------------------------------------

  subroutine dealloc_aer1_inorg(aer1,ninorg)

    implicit none

    integer,intent(in) :: ninorg
    type (aer1_vars)    ,dimension(ninorg) :: aer1
    integer :: ispc

    !  Deallocate arrays
    do ispc=1,ninorg

       if (associated(aer1(ispc)%sc_p ) )    deallocate(aer1(ispc)%sc_p  )
       if (associated(aer1(ispc)%sc_src))    deallocate(aer1(ispc)%sc_src)
       if (associated(aer1(ispc)%sc_dd) )    deallocate(aer1(ispc)%sc_dd )
       if (associated(aer1(ispc)%sc_wd) )    deallocate(aer1(ispc)%sc_wd )
       if (associated(aer1(ispc)%sc_pp) )    deallocate(aer1(ispc)%sc_pp )
       if (associated(aer1(ispc)%sc_pf) )    deallocate(aer1(ispc)%sc_pf )
    enddo

    return
  end subroutine dealloc_aer1_inorg

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine nullify_aer1_inorg(aer1,ninorg)

    implicit none

    integer,intent(in) ::ninorg
    type (aer1_vars),dimension(ninorg) :: aer1
    integer :: ispc

    do ispc=1,ninorg


       if (associated(aer1(ispc)%sc_p ) )    nullify (aer1(ispc)%sc_p  )
       if (associated(aer1(ispc)%sc_src))    nullify (aer1(ispc)%sc_src)
       if (associated(aer1(ispc)%sc_dd) )    nullify (aer1(ispc)%sc_dd )
       if (associated(aer1(ispc)%sc_wd) )    nullify (aer1(ispc)%sc_wd )
       if (associated(aer1(ispc)%sc_pp) )    nullify (aer1(ispc)%sc_pp )
       if (associated(aer1(ispc)%sc_pf) )    nullify (aer1(ispc)%sc_pf )

    enddo

    return
  end subroutine nullify_aer1_inorg

  !---------------------------------------------------------------

  subroutine filltab_aer1_inorg(oneVarTable, oneVarTableSize, &
       aer1, aer1m, imean)

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    use io_params, only : ioutput         ! INTENT(IN)

    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(aer1_vars), pointer, intent(in) :: aer1(:)
    type(aer1_vars), pointer, intent(in) :: aer1m(:)
    integer, intent(in) :: imean

    integer :: ispc,imode
    character(len=8) :: str_recycle
    character(len=1) :: str_src_dim

    str_recycle = ''; str_src_dim = ''
    if (RECYCLE_TRACERS == 1 .or. ioutput == 5) then
       str_recycle = ':recycle'
    endif

    !- Fill pointers to arrays into variable tables
    do ispc = 1, size(aer1)

       if (associated(aer1(ispc)%sc_p)) then

          !---- tracer mixing ratio (dimension 3d)
          call InsertVarTable(oneVarTable, oneVarTableSize, &
               aer1(ispc)%sc_p, &
               trim(inorg_spc_name(ispc))//'P :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle), &
               aer1m(ispc)%sc_p, imean)

          !---- sources (3 and 2 dimension)
          if(inorg_alloc(src,ispc) == on) then
             str_src_dim = '3'  ! for 3d sources
             call InsertVarTable(oneVarTable, oneVarTableSize, &
                  aer1(ispc)%sc_src, &
                  trim(inorg_spc_name(ispc))//'_SRC :'//trim(str_src_dim)//':hist:anal:mpti:mpt3:mpt1', &
                  aer1m(ispc)%sc_src, imean)
          endif

          !---- dry and wet deposition (dimension 2d)
          if(inorg_alloc(ddp,ispc) == on) then
             call InsertVarTable(oneVarTable, oneVarTableSize, &
                  aer1(ispc)%sc_dd, &
                  trim(inorg_spc_name(ispc))//'DD :2:hist:anal:mpti:mpt3', &
                  aer1m(ispc)%sc_dd, imean)
          end if

          if(inorg_alloc(wdp,ispc) == on) then
             call InsertVarTable(oneVarTable, oneVarTableSize, &
                  aer1(ispc)%sc_wd, &
                  trim(inorg_spc_name(ispc))//'WD :2:hist:anal:mpti:mpt3', &
                  aer1m(ispc)%sc_wd, imean)
          end if

          !----  data assimilation (dimension 3d)
          if(chem_assim == on) then
             if(inorg_alloc(fdda,ispc) == on) then
                call InsertVarTable(oneVarTable, oneVarTableSize, &
                     aer1(ispc)%sc_pp, &
                     trim(inorg_spc_name(ispc))//'PP :3:mpti', &
                     aer1m(ispc)%sc_pp, imean)
                call InsertVarTable(oneVarTable, oneVarTableSize, &
                     aer1(ispc)%sc_pf, &
                     trim(inorg_spc_name(ispc))//'PF :3:mpti', &
                     aer1m(ispc)%sc_pf, imean)
             end if
          end if
       end if
    end do
  end subroutine filltab_aer1_inorg


  !---------------------------------------------------------------

  subroutine alloc_tend_aer1_inorg(nmzp,nmxp,nmyp,ngrs,ninorg,proc_type)
    implicit none
    integer,intent(in)                   :: ngrs,proc_type,ninorg
    integer,intent(in), dimension (ngrs) :: nmzp,nmxp,nmyp
    integer :: ng,ntpts,ispc

    !         Find the maximum number of grid points needed for any grid.

    if(proc_type==1) then
       ntpts=1
    else
       ntpts=0
       do ng=1,ngrs
          ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
       enddo
    endif

!!!!!  WE ARE ONLY CHECKING GRID 1 !!!!!!!!!
!!!!!    All grids must have same scalars defined !!!!!!!
    do ispc=1,ninorg


       !1st test: if the mode does not exist, cycle
       if(inorg_mod_alloc(ispc) /= 1 ) cycle

       if (associated(aer1_inorg_g(ispc,1)%sc_p)) allocate (aer1_inorg_g(ispc,1)%sc_t(ntpts))

       do ng=2,ngrs
          aer1_inorg_g(ispc,ng)%sc_t => aer1_inorg_g(ispc,1)%sc_t
       enddo

    enddo

  end subroutine alloc_tend_aer1_inorg
  !---------------------------------------------------------------

  subroutine nullify_tend_aer1_inorg(ninorg)

    implicit none
    integer,intent(in) ::ninorg
    integer ::ispc

    do ispc=1,ninorg

       if (associated(aer1_inorg_g(ispc,1)%sc_t)) nullify (aer1_inorg_g(ispc,1)%sc_t)
    enddo


  end subroutine nullify_tend_aer1_inorg

  !---------------------------------------------------------------
  subroutine dealloc_tend_aer1_inorg(ninorg)
    implicit none
    integer,intent(in) :: ninorg
    integer ::ispc

    do ispc=1,ninorg

       if (associated(aer1_inorg_g(ispc,1)%sc_t)) deallocate (aer1_inorg_g(ispc,1)%sc_t)
    enddo

  end  subroutine dealloc_tend_aer1_inorg

  !---------------------------------------------------------------

  subroutine filltab_tend_aer1_inorg(oneScalarTab, oneScalarTabSize, &
       ninorg,ng)
    use mem_chem1, only: nspecies_transported ! this is first calculated at chemistry
    ! "filltab_tend_chem1" routine
    implicit none

    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(inout) :: oneScalarTabSize
    integer,intent(in) :: ninorg,ng
    integer ::ispc
    integer :: elements

    do ispc=1,ninorg
       ! Fill pointers to scalar arrays into scalar tables

       if (associated(aer1_inorg_g(ispc,ng)%sc_t)) then
          elements = size(aer1_inorg_g(ispc,ng)%sc_t)
          call InsertAtScalarTab(&
               aer1_inorg_g(ispc,ng)%sc_p, &
               aer1_inorg_g(ispc,ng)%sc_t, &
               trim(inorg_spc_name(ispc))//'P',&
               oneScalarTab, oneScalarTabSize)
          !- total number of transported species (CHEM + AER)
          nspecies_transported = nspecies_transported + 1

       endif

    enddo

  end subroutine filltab_tend_aer1_inorg


  !=============================================================================================
  subroutine dumpAer(fileName)
    !# Dum aerossois to a file
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: Dump aerossois to a file
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 30 September 2020 (Wednesday)
    !# @endnote
    !#
    !# @changes
    !# &#9744; <br/>
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#

    !Use area
    use dump, only: &
         dumpMessage

    use node_mod, only: &
         mzp, mxp, myp,  & ! INTENT(IN)
         izu, jzv,       & ! INTENT(IN)
         mynum,          & ! INTENT(IN)
         ibcon,          & ! INTENT(IN)
         nmachs,         &   ! INTENT(IN)
         mchnum, &
         master_num, &
         nodemxp, &
         nodemyp, &
         nodemzp, &
         nodei0, &
         nodej0, &
         ixb,ixe,iyb,iye

    use mem_grid, only: &
         ngrids, nnxp, nnyp, nzg, npatch, nzs, nnzp, &
         oneGlobalGridData, &
         iyear1,imonth1,idate1,ihour1,itime1 ! from RAMSIN

    implicit none

    include "constants.h"
    character(len=*),parameter :: sourceName='mem_chem1.f90' !Name of source code
    character(len=*),parameter :: procedureName='**dumpAer**' !Name of this procedure
    !
    !Local Parameters

    !Input/Output variables
    character(len=*), intent(in) :: fileName

    !Local variables
    integer :: recordLen,irec,nspc
    integer :: ispc,nm,k
    character(len=15) :: tdef
    real :: dlat,dlon

    !Code

    if(nmachs>1) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
         ,c_warning,'Funciona apenas para rodadas de 1 processador!... sorry...')

    write(tDef,fmt='(I2.2,":00z",I2.2,A3,I4.4)')  ihour1,idate1,month_Name(imonth1),iyear1

    nspc=0
    irec=1
    recordLen=4*mxp*myp
    open(unit=33,file=fileName//'.gra',action='WRITE',status='REPLACE' &
         ,form='UNFORMATTED',access='DIRECT',recl=recordLen)
    do ispc=1,nspecies_aer
       do nm=1,nmodes
          if(mode_alloc(nm,ispc) /= 1) cycle
          nspc=nspc+1
          do k=1,mzp
             write(33,rec=irec) aer1_g(nm,ispc,1)%sc_p(k,:,:)
             irec=irec+1
          enddo
       enddo
    enddo
    do ispc=1,nspecies_aer
       do nm=1,nmodes
          if(mode_alloc(nm,ispc) /= 1) cycle
          nspc=nspc+1
          do k=1,mzp
             write(33,rec=irec) aer1_g(nm,ispc,1)%sc_pf(k,:,:)
             irec=irec+1
          enddo
       enddo
    enddo
    close(33)

    dlon=oneGlobalGridData(1)%global_glon(2,1)-oneGlobalGridData(1)%global_glon(1,1)
    dlat=oneGlobalGridData(1)%global_glat(1,2)-oneGlobalGridData(1)%global_glat(1,1)

    open(33, file=fileName//'.ctl', action='write', status='replace')

    !writing the name of grads file
    write(33,*) 'dset ^'//fileName//'.gra'
    !writing others infos to ctl
    write(33,*) 'undef -0.9990000E+34'
    write(33,*) 'title Aer Test File'
    write(33,*) 'xdef ',mxp,' linear ',oneGlobalGridData(1)%global_glon(1,1),dlon
    write(33,*) 'ydef ',myp,' linear ',oneGlobalGridData(1)%global_glat(1,1),dlat
    write(33,*) 'zdef ',mzp,'levels',(k,k=1,mzp)
    write(33,*) 'tdef 1 linear '//tDef//' 1mo'
    write(33,*) 'vars ',nspc*2
    do ispc=1,nspecies_aer
       do nm=1,nmodes
          if(mode_alloc(nm,ispc) /= 1) cycle
          write(33,*) trim(aer_name(nm,ispc)),mzp,'99 ',trim(aer_name(nm,ispc))
       enddo
    enddo
    do ispc=1,nspecies_aer
       do nm=1,nmodes
          if(mode_alloc(nm,ispc) /= 1) cycle
          write(33,*) trim(aer_name(nm,ispc))//'_pf',mzp,'99 ',trim(aer_name(nm,ispc))//'_pf'
       enddo
    enddo
    write(33,*) 'endvars'
    close(33)

  end subroutine dumpAer


  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  subroutine define_aer1_src_zdim(aer1_src_z_dim,n1)

    implicit none
    integer, intent(out) :: aer1_src_z_dim(nspecies_aer)
    integer, intent(in) :: n1

    !- determination of the dimension of Z-dir of source field array
    aer1_src_z_dim = n1     ! 3d
    !  aer1_src_z_dim(bburn) = n1     ! 3d
    !  aer1_src_z_dim(urban) = n1     ! 2d
    !  aer1_src_z_dim(bioge) = n1     ! 3d
    !  aer1_src_z_dim(marin) = n1     ! 3d
    !  aer1_src_z_dim(v_ash) = n1     ! 3d
    return
  end subroutine define_aer1_src_zdim



  subroutine FixChemAerVarTableForIOUTPUT5(oneVarTable, oneVarTableSize, chemistry, aerosol)
    implicit none
    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(in) :: oneVarTableSize
    integer, intent(in) :: chemistry
    integer, intent(in) :: aerosol
    integer :: ni, nspc, imode
    character(len=2) :: cmode

    do ni = 1, oneVarTableSize

       if (trim(oneVarTable(ni)%name) == 'TOPT'  .or. &
	    trim(oneVarTable(ni)%name) == 'UP'    .or. &
	    trim(oneVarTable(ni)%name) == 'VP'    .or. &
	    trim(oneVarTable(ni)%name) == 'THETA' .or. &
	    trim(oneVarTable(ni)%name) == 'PP'    .or. &
	    trim(oneVarTable(ni)%name) == 'RV') then
          cycle
       end if


       if(chemistry >= 0) then 
          do nspc=1,chem_nspecies
             if(chem_alloc(chem_fdda,nspc) == chem_on .and. &
                  trim(oneVarTable(ni)%name) == trim(chem_name(nspc))//'P') then 
                oneVarTable(ni)%ianal=1
                cycle
             end if
          end do
       end if

       if(aerosol == 1 .and. chemistry >= 0) then
          do nspc=1, nspecies_aer
             do imode = 1, nmodes
                write(cmode, '(BN, I2)')imode
                cmode = adjustl(cmode)
                if(spc_alloc(fdda,imode,nspc) == 1  .and. &
                     trim(oneVarTable(ni)%name) == trim(aer_name(imode,nspc))//trim(cmode)//'P') then
                   oneVarTable(ni)%ianal=1
                   cycle
                end if
             end do
          end do
       end if
    end do
  end subroutine FixChemAerVarTableForIOUTPUT5

end module mem_aer1

