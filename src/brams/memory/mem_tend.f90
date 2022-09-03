!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module mem_tend

    use ModScalarTable, only:  &
         ScalarTable, &
         InsertAtScalarTab 

    use ModBasicFields, only: &
         BasicFields
    
    use ModMicroFields, only: &
         MicroFields
    
    use ModTurbFields, only: &
         TurbFields

    use ModTurbFields, only: &
         TurbFields

    use ModScalarFields, only: &
         ScalarFields
    
    use mem_grid,   only: &
         dyncore_flag

    use ModGaspartFields, only: &
         GaspartFields
    
    use teb_spm_start, only: &
         TEB_SPM
    
    use mem_emiss, only: &
         ichemi, &
         isource

    use ModNamelistFile, only: &
         namelistFile

    implicit none

    private

    public :: tend_vars
    public :: tend
    public :: alloc_tend
    public :: nullify_tend
    public :: dealloc_tend
    public :: filltab_tend
    
  type tend_vars
     real, contiguous, pointer :: ut(:)
     real, contiguous, pointer :: vt(:)
     real, contiguous, pointer :: wt(:)
     real, contiguous, pointer :: pt(:)
     real, contiguous, pointer :: tht(:)
     real, contiguous, pointer :: rtt(:)
     real, contiguous, pointer :: rct(:)
     real, contiguous, pointer :: rrt(:)
     real, contiguous, pointer :: rpt(:)
     real, contiguous, pointer :: rst(:)
     real, contiguous, pointer :: rat(:)
     real, contiguous, pointer :: rgt(:)
     real, contiguous, pointer :: rht(:)
     real, contiguous, pointer :: cct(:)
     real, contiguous, pointer :: crt(:)
     real, contiguous, pointer :: cpt(:)
     real, contiguous, pointer :: cst(:)
     real, contiguous, pointer :: cat(:)
     real, contiguous, pointer :: cgt(:)
     real, contiguous, pointer :: cht(:)
     real, contiguous, pointer :: cccnt(:)
     real, contiguous, pointer :: cifnt(:)
     real, contiguous, pointer :: tket(:)
     real, contiguous, pointer :: epst(:)
     real, contiguous, pointer :: rdt(:)
     real, contiguous, pointer :: cdt(:)
     real, contiguous, pointer :: gccnt(:)
     real, contiguous, pointer :: cccmt(:)
     real, contiguous, pointer :: gccmt(:)
     real, contiguous, pointer :: cnm1t(:)
     real, contiguous, pointer :: cnm2t(:)
     real, contiguous, pointer :: cnm3t(:)
     real, contiguous, pointer :: cnm8t(:)
     real, contiguous, pointer :: md1nt(:)
     real, contiguous, pointer :: md2nt(:)
     real, contiguous, pointer :: salt_filmt(:)
     real, contiguous, pointer :: salt_jett(:)
     real, contiguous, pointer :: salt_spmt(:)
     real, contiguous, pointer :: ut_rk(:)
     real, contiguous, pointer :: vt_rk(:)
     real, contiguous, pointer :: wt_rk(:)
     real, contiguous, pointer :: pt_rk(:)
     real, contiguous, pointer :: tht_rk(:)
     real, contiguous, pointer :: ut_past(:)
     real, contiguous, pointer :: vt_past(:)
     real, contiguous, pointer :: wt_past(:)
     real, contiguous, pointer :: pt_past(:)
     real, contiguous, pointer :: tht_past(:)
     real, contiguous, pointer :: cldfrt(:)
  end type tend_vars

  type (tend_vars) :: tend

contains
  !---------------------------------------------------------------

  subroutine alloc_tend(nmzp,nmxp,nmyp,ngrs,naddsc,proc_type,&
       oneBasicFields, oneTurbFields, oneMicroFields)
    ! Arguments:
    integer, intent(in) :: nmzp(:), nmxp(:), nmyp(:)
    integer, intent(in) :: ngrs, proc_type, naddsc
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields
    type(MicroFields), pointer, intent(in) :: oneMicroFields

    ! Local Variables:
    integer :: ng, ntpts, nsc

    character(len=*), parameter :: h="**(alloc_tend)**"
    !         Find the maximum number of grid points needed for any grid.

    if(proc_type==1) then
       ntpts=1
    else
       ntpts=0
       do ng=1,ngrs
          ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
       enddo
    endif

    ! Allocate arrays based on options (if necessary).
    !   Do not need these arrays if it is a master process in a parallel run,
    !      so just allocate to 1 word.

!!!!!  WE ARE ONLY CHECKING GRID 1 !!!!!!!!!
!!!!!    All grids must have same scalars defined !!!!!!!

    if (associated(oneBasicFields%up))      then
       allocate (tend%ut(ntpts))
       tend%ut = 0.
       if( dyncore_flag==2 .or. dyncore_flag==3) then
          allocate (tend%ut_rk(ntpts))
          tend%ut_rk = 0.
       endif
       if( dyncore_flag==3 ) then
          allocate (tend%ut_past(ntpts))
          tend%ut_past = 0.
       endif
    endif
    if (associated(oneBasicFields%vp))      then
       allocate (tend%vt(ntpts))
       tend%vt = 0.
       if( dyncore_flag==2 .or. dyncore_flag==3) then
          allocate (tend%vt_rk(ntpts))
          tend%vt_rk = 0.
       endif
       if( dyncore_flag==3 ) then
          allocate (tend%vt_past(ntpts))
          tend%vt_past = 0.
       endif
    endif
    if (associated(oneBasicFields%wp))      then
       allocate (tend%wt(ntpts))
       tend%wt = 0.
       if( dyncore_flag==2 .or. dyncore_flag==3) then
          allocate (tend%wt_rk(ntpts))
          tend%wt_rk = 0.
       endif
       if( dyncore_flag==3 ) then
          allocate (tend%wt_past(ntpts))
          tend%wt_past = 0.
       endif
    endif
    if (associated(oneBasicFields%pp))      then
       allocate (tend%pt(ntpts))
       tend%pt = 0.
       if( dyncore_flag==2 .or. dyncore_flag==3) then
          allocate (tend%pt_rk(ntpts))
          tend%pt_rk = 0.
       endif
       if( dyncore_flag==3 ) then
          allocate (tend%pt_past(ntpts))
          tend%pt_past = 0.
       endif
    endif

    if (associated(oneBasicFields%thp))     then
       allocate (tend%tht(ntpts))
       tend%tht = 0.
       if( dyncore_flag==2 .or. dyncore_flag==3) then
          allocate (tend%tht_rk(ntpts))
          tend%tht_rk = 0.
       endif
       if( dyncore_flag==3 ) then
          allocate (tend%tht_past(ntpts))
          tend%tht_past = 0.
       endif
    endif
    if (associated(oneBasicFields%rtp))     then
       allocate (tend%rtt(ntpts))
       tend%rtt = 0.
    endif
    if (associated(oneMicroFields%rcp))     then
       allocate (tend%rct(ntpts))
       tend%rct = 0.
    endif
    if (associated(oneMicroFields%rrp))     then
       allocate (tend%rrt(ntpts))
       tend%rrt = 0.
    endif
    if (associated(oneMicroFields%rpp))     then
       allocate (tend%rpt(ntpts))
       tend%rpt = 0.
    endif
    if (associated(oneMicroFields%rsp))     then
       allocate (tend%rst(ntpts))
       tend%rst = 0.
    endif
    if (associated(oneMicroFields%rap))     then
       allocate (tend%rat(ntpts))
       tend%rat = 0.
    endif
    if (associated(oneMicroFields%rgp))     then
       allocate (tend%rgt(ntpts))
       tend%rgt = 0.
    endif
    if (associated(oneMicroFields%rhp))     then
       allocate (tend%rht(ntpts))
       tend%rht = 0.
    endif
    if (associated(oneMicroFields%ccp))     then
       allocate (tend%cct(ntpts))
       tend%cct = 0.
    endif
    if (associated(oneMicroFields%crp))     then
       allocate (tend%crt(ntpts))
       tend%crt = 0.
    endif
    if (associated(oneMicroFields%cpp))     then
       allocate (tend%cpt(ntpts))
       tend%cpt = 0.
    endif
    if (associated(oneMicroFields%csp))     then
       allocate (tend%cst(ntpts))
       tend%cst = 0.
    endif
    if (associated(oneMicroFields%cap))     then
       allocate (tend%cat(ntpts))
       tend%cat = 0.
    endif
    if (associated(oneMicroFields%cgp))     then
       allocate (tend%cgt(ntpts))
       tend%cgt = 0.
    endif
    if (associated(oneMicroFields%chp))     then
       allocate (tend%cht(ntpts))
       tend%cht = 0.
    endif
    if (associated(oneMicroFields%cccnp))   then
       allocate (tend%cccnt(ntpts))
       tend%cccnt = 0.
    endif
    if (associated(oneMicroFields%cifnp))   then
       allocate (tend%cifnt(ntpts))
       tend%cifnt = 0.
    endif
    if (associated(oneMicroFields%cldfr))   then
       allocate (tend%cldfrt(ntpts))
       tend%cldfrt = 0.
    endif
    if (associated(oneTurbFields%tkep))     then
       allocate (tend%tket(ntpts))
       tend%tket = 0.
    endif
    if (associated(oneTurbFields%epsp))     then
       allocate (tend%epst(ntpts))
       tend%epst = 0.
    endif
    !
    !-2015- for 2M microphysics (from G. Camponogara)
    if (associated(oneMicroFields%rdp))     then
       allocate (tend%rdt(ntpts))  ;tend%rdt=0.0
    endif

    if (associated(oneMicroFields%cdp))     then
       allocate (tend%cdt(ntpts))  ;tend%cdt=0.0
    endif

    if (associated(oneMicroFields%gccnp))   then
       allocate (tend%gccnt(ntpts))  ;tend%gccnt=0.0
    endif

    if (associated(oneMicroFields%cccmp))   then
       allocate (tend%cccmt(ntpts))  ;tend%cccmt=0.0
    endif

    if (associated(oneMicroFields%gccmp))   then
       allocate (tend%gccmt(ntpts))  ;tend%gccmt=0.0
    endif

    if (associated(oneMicroFields%cnm1p))   then
       allocate (tend%cnm1t(ntpts))  ;tend%cnm1t=0.0
    endif

    if (associated(oneMicroFields%cnm2p))   then
       allocate (tend%cnm2t(ntpts))  ;tend%cnm2t=0.0
    endif
    if (associated(oneMicroFields%cnm3p))   then
       allocate (tend%cnm3t(ntpts))  ;tend%cnm3t=0.0
    endif

    if (associated(oneMicroFields%cnm8p))   then
       allocate (tend%cnm8t(ntpts))  ;tend%cnm8t=0.0
    endif

    if (associated(oneMicroFields%md1np))   then
       allocate (tend%md1nt(ntpts))  ;tend%md1nt=0.0
    endif

    if (associated(oneMicroFields%md2np))   then
       allocate (tend%md2nt(ntpts))  ;tend%md2nt=0.0
    endif

    if (associated(oneMicroFields%salt_filmp)) then
       allocate (tend%salt_filmt(ntpts))  ;tend%salt_filmt=0.0
    endif
    if (associated(oneMicroFields%salt_jetp))  then
       allocate (tend%salt_jett(ntpts))  ;tend%salt_jett=0.0
    endif

    if (associated(oneMicroFields%salt_spmp))  then
       allocate (tend%salt_spmt(ntpts))  ;tend%salt_spmt=0.0
    endif
    !GC - 2M microphysics

    !--(DMK-LFR NEC-SX6)----------------------------------------------

    ! TEB_SPM
    if (TEB_SPM==1) then
       if(isource==1)then
          call fatal_error(h//" allocation of gaspart not converted; "//&
               "seems wrong whenever TEB_SPM and isource are both 1")
!!$          if (associated(gaspart_g(1)%pno).and.  &
!!$               (.not.associated(gaspart_g(1)%pnot)))  &
!!$               allocate (gaspart_g(1)%pnot(ntpts))
!!$
!!$          if (associated(gaspart_g(1)%pno2).and.  &
!!$               (.not.associated(gaspart_g(1)%pno2t)))  &
!!$               allocate (gaspart_g(1)%pno2t(ntpts))
!!$
!!$          if (associated(gaspart_g(1)%ppm25).and.  &
!!$               (.not.associated(gaspart_g(1)%ppm25t)))  &
!!$               allocate (gaspart_g(1)%ppm25t(ntpts))
!!$
!!$          if (associated(gaspart_g(1)%pco).and.  &
!!$               (.not.associated(gaspart_g(1)%pcot)))  &
!!$               allocate (gaspart_g(1)%pcot(ntpts))
!!$
!!$          if (associated(gaspart_g(1)%pso2).and.  &
!!$               (.not.associated(gaspart_g(1)%pso2t)))  &
!!$               allocate (gaspart_g(1)%pso2t(ntpts))
!!$
!!$          if (associated(gaspart_g(1)%pso4).and.  &
!!$               (.not.associated(gaspart_g(1)%pso4t)))  &
!!$               allocate (gaspart_g(1)%pso4t(ntpts))
!!$
!!$          if (associated(gaspart_g(1)%paer).and.  &
!!$               (.not.associated(gaspart_g(1)%paert)))  &
!!$               allocate (gaspart_g(1)%paert(ntpts))
!!$
!!$          if (associated(gaspart_g(1)%pvoc).and.  &
!!$               (.not.associated(gaspart_g(1)%pvoct)))  &
!!$               allocate (gaspart_g(1)%pvoct(ntpts))
!!$
!!$          if(ichemi==1)then
!!$
!!$             if (associated(gaspart_g(1)%po3).and.  &
!!$                  (.not.associated(gaspart_g(1)%po3t)))  &
!!$                  allocate (gaspart_g(1)%po3t(ntpts))
!!$
!!$             if (associated(gaspart_g(1)%prhco).and.  &
!!$                  (.not.associated(gaspart_g(1)%prhcot)))  &
!!$                  allocate (gaspart_g(1)%prhcot(ntpts))
!!$
!!$             if (associated(gaspart_g(1)%pho2).and.  &
!!$                  (.not.associated(gaspart_g(1)%pho2t)))  &
!!$                  allocate (gaspart_g(1)%pho2t(ntpts))
!!$
!!$             if (associated(gaspart_g(1)%po3p).and.  &
!!$                  (.not.associated(gaspart_g(1)%po3pt)))  &
!!$                  allocate (gaspart_g(1)%po3pt(ntpts))
!!$
!!$             if (associated(gaspart_g(1)%po1d).and.  &
!!$                  (.not.associated(gaspart_g(1)%po1dt)))  &
!!$                  allocate (gaspart_g(1)%po1dt(ntpts))
!!$
!!$             if (associated(gaspart_g(1)%pho).and.  &
!!$                  (.not.associated(gaspart_g(1)%phot)))  &
!!$                  allocate (gaspart_g(1)%phot(ntpts))
!!$
!!$             if (associated(gaspart_g(1)%proo).and.  &
!!$                  (.not.associated(gaspart_g(1)%proot)))  &
!!$                  allocate (gaspart_g(1)%proot(ntpts))
!!$          endif
!!$
!!$          do ng=2,ngrs
!!$             gaspart_g(ng)%pnot   => gaspart_g(1)%pnot
!!$             gaspart_g(ng)%pno2t  => gaspart_g(1)%pno2t
!!$             gaspart_g(ng)%ppm25t => gaspart_g(1)%ppm25t
!!$             gaspart_g(ng)%pcot   => gaspart_g(1)%pcot
!!$             gaspart_g(ng)%pso2t  => gaspart_g(1)%pso2t
!!$             gaspart_g(ng)%pso4t  => gaspart_g(1)%pso4t
!!$             gaspart_g(ng)%paert  => gaspart_g(1)%paert
!!$             gaspart_g(ng)%pvoct  => gaspart_g(1)%pvoct
!!$             if(ichemi==1)then
!!$                gaspart_g(ng)%po3t   => gaspart_g(1)%po3t
!!$                gaspart_g(ng)%prhcot => gaspart_g(1)%prhcot
!!$                gaspart_g(ng)%pho2t  => gaspart_g(1)%pho2t
!!$                gaspart_g(ng)%po3pt  => gaspart_g(1)%po3pt
!!$                gaspart_g(ng)%po1dt  => gaspart_g(1)%po1dt
!!$                gaspart_g(ng)%phot   => gaspart_g(1)%phot
!!$                gaspart_g(ng)%proot  => gaspart_g(1)%proot
!!$             endif
!!$          enddo
!!$
       endif

    endif
    !

!!$    do nsc=1,naddsc
!!$       if (associated(scalar_g(nsc,1)%sclp).and.  &
!!$            (.not.associated(scalar_g(nsc,1)%sclt)))  &
!!$            allocate (scalar_g(nsc,1)%sclt(ntpts))
!!$       do ng=2,ngrs
!!$          scalar_g(nsc,ng)%sclt => scalar_g(nsc,1)%sclt
!!$       enddo
!!$    enddo

  end subroutine alloc_tend

  !---------------------------------------------------------------

  subroutine nullify_tend(naddsc)
    ! Arguments:
    integer, intent(in) :: naddsc

    ! Local Variables:
    integer :: nsc
    character(len=*), parameter :: h="**(nullify_tend)**"

    ! Deallocate all tendency arrays

    if (associated(tend%ut))   nullify (tend%ut)
    if (associated(tend%vt))   nullify (tend%vt)
    if (associated(tend%wt))   nullify (tend%wt)
    if (associated(tend%pt))   nullify (tend%pt)
    if (associated(tend%tht))  nullify (tend%tht)
    if (associated(tend%rtt))  nullify (tend%rtt)
    if (associated(tend%rct))  nullify (tend%rct)
    if (associated(tend%rrt))  nullify (tend%rrt)
    if (associated(tend%rpt))  nullify (tend%rpt)
    if (associated(tend%rst))  nullify (tend%rst)
    if (associated(tend%rat))  nullify (tend%rat)
    if (associated(tend%rgt))  nullify (tend%rgt)
    if (associated(tend%rht))  nullify (tend%rht)
    if (associated(tend%cct))  nullify (tend%cct)
    if (associated(tend%crt))  nullify (tend%crt)
    if (associated(tend%cpt))  nullify (tend%cpt)
    if (associated(tend%cst))  nullify (tend%cst)
    if (associated(tend%cat))  nullify (tend%cat)
    if (associated(tend%cgt))  nullify (tend%cgt)
    if (associated(tend%cht))  nullify (tend%cht)
    if (associated(tend%cccnt))nullify (tend%cccnt)
    if (associated(tend%cifnt))nullify (tend%cifnt)
    if (associated(tend%cldfrt))nullify (tend%cldfrt)
    if (associated(tend%tket)) nullify (tend%tket)
    if (associated(tend%epst)) nullify (tend%epst)

    !-2015- for 2M microphysics (from G. Camponogara)
    if (associated(tend%rdt  )) nullify(tend%rdt)
    if (associated(tend%cdt  )) nullify(tend%cdt)
    if (associated(tend%gccnt)) nullify(tend%gccnt)
    if (associated(tend%cccmt)) nullify(tend%cccmt)
    if (associated(tend%gccmt)) nullify(tend%gccmt)
    if (associated(tend%cnm1t)) nullify(tend%cnm1t)
    if (associated(tend%cnm2t)) nullify(tend%cnm2t)
    if (associated(tend%cnm3t)) nullify(tend%cnm3t)
    if (associated(tend%cnm8t)) nullify(tend%cnm8t)
    if (associated(tend%md1nt)) nullify(tend%md1nt)
    if (associated(tend%md2nt)) nullify(tend%md2nt)
    if (associated(tend%salt_filmt ))  nullify(tend%salt_filmt)
    if (associated(tend%salt_jett  ))  nullify(tend%salt_jett)
    if (associated(tend%salt_spmt  ))  nullify(tend%salt_spmt)
    !-2015- for 2M microphysics (from G. Camponogara)
    !- for RK/ABM3 method
    if (associated(tend%ut_rk))   nullify (tend%ut_rk)
    if (associated(tend%vt_rk))   nullify (tend%vt_rk)
    if (associated(tend%wt_rk))   nullify (tend%wt_rk)
    if (associated(tend%pt_rk))   nullify (tend%pt_rk)
    if (associated(tend%tht_rk))  nullify (tend%tht_rk)
    !- for ABM3 method
    if (associated(tend%ut_past))   nullify (tend%ut_past)
    if (associated(tend%vt_past))   nullify (tend%vt_past)
    if (associated(tend%wt_past))   nullify (tend%wt_past)
    if (associated(tend%pt_past))   nullify (tend%pt_past)
    if (associated(tend%tht_past))  nullify (tend%tht_past)

    ! TEB_SPM
    if (TEB_SPM==1) then
       if(isource==1)then
          call fatal_error(h//" nullify gaspart not converted; "//&
               "seems wrong whenever TEB_SPM and isource are both 1")
!!$          if (associated(gaspart_g(1)%pnot  )) nullify (gaspart_g(1)%pnot  )
!!$          if (associated(gaspart_g(1)%pno2t )) nullify (gaspart_g(1)%pno2t )
!!$          if (associated(gaspart_g(1)%ppm25t)) nullify (gaspart_g(1)%ppm25t)
!!$          if (associated(gaspart_g(1)%pcot  )) nullify (gaspart_g(1)%pcot  )
!!$          if (associated(gaspart_g(1)%pso2t )) nullify (gaspart_g(1)%pso2t )
!!$          if (associated(gaspart_g(1)%pso4t )) nullify (gaspart_g(1)%pso4t )
!!$          if (associated(gaspart_g(1)%paert )) nullify (gaspart_g(1)%paert )
!!$          if (associated(gaspart_g(1)%pvoct )) nullify (gaspart_g(1)%pvoct )
!!$
!!$          if(ichemi==1)then
!!$             if (associated(gaspart_g(1)%po3t  )) nullify (gaspart_g(1)%po3t  )
!!$             if (associated(gaspart_g(1)%prhcot)) nullify (gaspart_g(1)%prhcot)
!!$             if (associated(gaspart_g(1)%pho2t )) nullify (gaspart_g(1)%pho2t )
!!$             if (associated(gaspart_g(1)%po3pt )) nullify (gaspart_g(1)%po3pt )
!!$             if (associated(gaspart_g(1)%po1dt )) nullify (gaspart_g(1)%po1dt )
!!$             if (associated(gaspart_g(1)%phot  )) nullify (gaspart_g(1)%phot  )
!!$             if (associated(gaspart_g(1)%proot )) nullify (gaspart_g(1)%proot )
!!$          endif
       endif
    endif
    !

!!$    do nsc=1,naddsc
!!$       if (associated(scalar_g(nsc,1)%sclt)) nullify (scalar_g(nsc,1)%sclt)
!!$    enddo

    return
  end subroutine nullify_tend
  !---------------------------------------------------------------

  subroutine dealloc_tend(naddsc)
    ! Arguments:
    integer, intent(in) :: naddsc

    ! Local Variables:
    integer :: nsc
    character(len=*), parameter :: h="**(dealloc_tend)**"
    
    ! Deallocate all tendency arrays

    if (associated(tend%ut))   deallocate (tend%ut)
    if (associated(tend%vt))   deallocate (tend%vt)
    if (associated(tend%wt))   deallocate (tend%wt)
    if (associated(tend%pt))   deallocate (tend%pt)
    if (associated(tend%tht))  deallocate (tend%tht)
    if (associated(tend%rtt))  deallocate (tend%rtt)
    if (associated(tend%rct))  deallocate (tend%rct)
    if (associated(tend%rrt))  deallocate (tend%rrt)
    if (associated(tend%rpt))  deallocate (tend%rpt)
    if (associated(tend%rst))  deallocate (tend%rst)
    if (associated(tend%rat))  deallocate (tend%rat)
    if (associated(tend%rgt))  deallocate (tend%rgt)
    if (associated(tend%rht))  deallocate (tend%rht)
    if (associated(tend%cct))  deallocate (tend%cct)
    if (associated(tend%crt))  deallocate (tend%crt)
    if (associated(tend%cpt))  deallocate (tend%cpt)
    if (associated(tend%cst))  deallocate (tend%cst)
    if (associated(tend%cat))  deallocate (tend%cat)
    if (associated(tend%cgt))  deallocate (tend%cgt)
    if (associated(tend%cht))  deallocate (tend%cht)
    if (associated(tend%cccnt))deallocate (tend%cccnt)
    if (associated(tend%cifnt))deallocate (tend%cifnt)
    if (associated(tend%cldfrt))deallocate (tend%cldfrt)
    if (associated(tend%tket)) deallocate (tend%tket)
    if (associated(tend%epst)) deallocate (tend%epst)

    !-2015- for 2M microphysics (from G. Camponogara)
    if (associated(tend%rdt  )) deallocate(tend%rdt)
    if (associated(tend%cdt  )) deallocate(tend%cdt)
    if (associated(tend%gccnt)) deallocate(tend%gccnt)
    if (associated(tend%cccmt)) deallocate(tend%cccmt)
    if (associated(tend%gccmt)) deallocate(tend%gccmt)
    if (associated(tend%cnm1t)) deallocate(tend%cnm1t)
    if (associated(tend%cnm2t)) deallocate(tend%cnm2t)
    if (associated(tend%cnm3t)) deallocate(tend%cnm3t)
    if (associated(tend%cnm8t)) deallocate(tend%cnm8t)
    if (associated(tend%md1nt)) deallocate(tend%md1nt)
    if (associated(tend%md2nt)) deallocate(tend%md2nt)
    if (associated(tend%salt_filmt ))  deallocate(tend%salt_filmt)
    if (associated(tend%salt_jett  ))  deallocate(tend%salt_jett)
    if (associated(tend%salt_spmt  ))  deallocate(tend%salt_spmt)
    !-2015- for 2M microphysics (from G. Camponogara)

    ! TEB_SPM
    if (TEB_SPM==1) then
       if(isource==1)then
          call fatal_error(h//" deallocation of gaspart not converted; "//&
               "seems wrong whenever TEB_SPM and isource are both 1")
!!$          if (associated(gaspart_g(1)%pnot  )) deallocate (gaspart_g(1)%pnot  )
!!$          if (associated(gaspart_g(1)%pno2t )) deallocate (gaspart_g(1)%pno2t )
!!$          if (associated(gaspart_g(1)%ppm25t)) deallocate (gaspart_g(1)%ppm25t)
!!$          if (associated(gaspart_g(1)%pcot  )) deallocate (gaspart_g(1)%pcot  )
!!$          if (associated(gaspart_g(1)%pso2t )) deallocate (gaspart_g(1)%pso2t )
!!$          if (associated(gaspart_g(1)%pso4t )) deallocate (gaspart_g(1)%pso4t )
!!$          if (associated(gaspart_g(1)%paert )) deallocate (gaspart_g(1)%paert )
!!$          if (associated(gaspart_g(1)%pvoct )) deallocate (gaspart_g(1)%pvoct )
!!$
!!$          if(ichemi==1)then
!!$             if (associated(gaspart_g(1)%po3t  )) deallocate (gaspart_g(1)%po3t)
!!$             if (associated(gaspart_g(1)%prhcot)) deallocate (gaspart_g(1)%prhcot)
!!$             if (associated(gaspart_g(1)%pho2t )) deallocate (gaspart_g(1)%pho2t)
!!$             if (associated(gaspart_g(1)%po3pt )) deallocate (gaspart_g(1)%po3pt)
!!$             if (associated(gaspart_g(1)%po1dt )) deallocate (gaspart_g(1)%po1dt)
!!$             if (associated(gaspart_g(1)%phot  )) deallocate (gaspart_g(1)%phot)
!!$             if (associated(gaspart_g(1)%proot )) deallocate (gaspart_g(1)%proot)
!!$          endif

       endif
    endif

!!$    do nsc=1,naddsc
!!$       if (associated(scalar_g(nsc,1)%sclt)) deallocate (scalar_g(nsc,1)%sclt)
!!$    enddo
    if (associated(tend%ut_rk))   deallocate (tend%ut_rk)
    if (associated(tend%vt_rk))   deallocate (tend%vt_rk)
    if (associated(tend%wt_rk))   deallocate (tend%wt_rk)
    if (associated(tend%pt_rk))   deallocate (tend%pt_rk)
    if (associated(tend%tht_rk))  deallocate (tend%tht_rk)

    if (associated(tend%ut_past))   deallocate (tend%ut_past)
    if (associated(tend%vt_past))   deallocate (tend%vt_past)
    if (associated(tend%wt_past))   deallocate (tend%wt_past)
    if (associated(tend%pt_past))   deallocate (tend%pt_past)
    if (associated(tend%tht_past))  deallocate (tend%tht_past)

    return
  end subroutine dealloc_tend

  !---------------------------------------------------------------

  subroutine filltab_tend(oneScalarTab, oneScalarTabSize, &
       oneBasicFields, oneMicroFields, oneTurbFields, oneGaspartFields, &
       oneScalarFields, naddsc, ng)
    ! Arguments:
    type(ScalarTable), pointer, intent(in) :: oneScalarTab(:)
    integer, intent(inout) :: oneScalarTabSize
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(MicroFields), intent(in) :: oneMicroFields 
    type(TurbFields), pointer, intent(in) :: oneTurbFields
    type(GaspartFields), pointer, intent(in) :: oneGaspartFields
    type(ScalarFields), pointer, intent(in) :: oneScalarFields(:)
    integer, intent(in)             :: naddsc, ng

    ! Local Variables:
    integer :: nsc
    character (len=7) :: sname
    character(len=*), parameter :: h="**(filltab_tend)**"

    ! Fill pointers to scalar arrays into scalar tables

    !if ( ( dyncore_flag==0) .or. (dyncore_flag==1) ) then
    !- not need to include THC in scalar table, including only THP
    !- will works better for leapfrog and RK schemes
    if (associated(tend%tht)) then
       call InsertAtScalarTab(oneBasicFields%thp,tend%tht, 'THP', &
            oneScalarTab, oneScalarTabSize)     
    endif

    if (associated(tend%rtt)) then
       call InsertAtScalarTab(oneBasicFields%rtp,tend%rtt, 'RTP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%rct)) then
       call InsertAtScalarTab(oneMicroFields%rcp,tend%rct, 'RCP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%rrt)) then
       call InsertAtScalarTab(oneMicroFields%rrp,tend%rrt, 'RRP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%rpt)) then
       call InsertAtScalarTab(oneMicroFields%rpp,tend%rpt, 'RPP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%rst)) then
       call InsertAtScalarTab(oneMicroFields%rsp,tend%rst, 'RSP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%rat)) then
       call InsertAtScalarTab(oneMicroFields%rap,tend%rat, 'RAP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%rgt)) then
       call InsertAtScalarTab(oneMicroFields%rgp,tend%rgt, 'RGP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%rht)) then
       call InsertAtScalarTab(oneMicroFields%rhp,tend%rht, 'RHP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%cct)) then
       call InsertAtScalarTab(oneMicroFields%ccp,tend%cct, 'CCP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%crt)) then
       call InsertAtScalarTab(oneMicroFields%crp,tend%crt, 'CRP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%cpt)) then
       call InsertAtScalarTab(oneMicroFields%cpp,tend%cpt, 'CPP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%cst)) then
       call InsertAtScalarTab(oneMicroFields%csp,tend%cst, 'CSP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%cat)) then
       call InsertAtScalarTab(oneMicroFields%cap,tend%cat, 'CAP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%cgt)) then
       call InsertAtScalarTab(oneMicroFields%cgp,tend%cgt, 'CGP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%cht)) then
       call InsertAtScalarTab(oneMicroFields%chp,tend%cht, 'CHP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%cccnt)) then
       call InsertAtScalarTab(oneMicroFields%cccnp,tend%cccnt, 'CCCNP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%cifnt)) then
       call InsertAtScalarTab(oneMicroFields%cifnp,tend%cifnt, 'CIFNP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if (associated(tend%cldfrt)) then
       call InsertAtScalarTab(oneMicroFields%cldfr,tend%cldfrt, 'CLDFR', &
            oneScalarTab, oneScalarTabSize)
    endif

    if( associated(tend%tket)) then
       call InsertAtScalarTab(oneTurbFields%tkep,tend%tket, 'TKEP', &
            oneScalarTab, oneScalarTabSize)
    endif

    if( associated(tend%epst)) then
       call InsertAtScalarTab(oneTurbFields%epsp,tend%epst, 'EPSP', &
            oneScalarTab, oneScalarTabSize)
    endif



    !-2015- for 2M microphysics (from G. Camponogara)
    if (associated(tend%rdt  )) then
       call InsertAtScalarTab(oneMicroFields%rdp,tend%rdt, 'RDP', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%cdt  )) then
       call InsertAtScalarTab(oneMicroFields%cdp,tend%cdt, 'CDP', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%gccnt)) then
       call InsertAtScalarTab(oneMicroFields%gccnp,tend%gccnt, 'GCCNP', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%cccmt))  then
       call InsertAtScalarTab(oneMicroFields%cccmp,tend%cccmt, 'CCCMP', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%gccmt))  then
       call InsertAtScalarTab(oneMicroFields%gccmp,tend%gccmt, 'GCCMP', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%cnm1t))  then
       call InsertAtScalarTab(oneMicroFields%cnm1p,tend%cnm1t, 'CNM1P', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%cnm2t))  then
       call InsertAtScalarTab(oneMicroFields%cnm2p,tend%cnm2t, 'CNM2P', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%cnm3t))  then
       call InsertAtScalarTab(oneMicroFields%cnm3p,tend%cnm3t, 'CNM3P', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%cnm8t))  then
       call InsertAtScalarTab(oneMicroFields%cnm8p,tend%cnm8t, 'CNM8P', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%md1nt))  then
       call InsertAtScalarTab(oneMicroFields%md1np,tend%md1nt, 'MD1NP', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%md2nt))  then
       call InsertAtScalarTab(oneMicroFields%md2np,tend%md2nt, 'MD2NP', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%salt_filmt))  then
       call InsertAtScalarTab(oneMicroFields%salt_filmp,tend%salt_filmt, 'SALT_FILMP', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%salt_jett))  then
       call InsertAtScalarTab(oneMicroFields%salt_jetp,tend%salt_jett, 'SALT_JETP', &
            oneScalarTab, oneScalarTabSize)
    endif
    if (associated(tend%salt_spmt))  then
       call InsertAtScalarTab(oneMicroFields%salt_spmp,tend%salt_spmt, 'SALT_SPMP', &
            oneScalarTab, oneScalarTabSize)
    endif
    !-2015- for 2M microphysics (from G. Camponogara)


    ! TEB_SPM
    if (TEB_SPM==1) then
       if (.not. associated(oneGaspartFields)) then
          call fatal_error(h//" oneGaspartFields not associated, but TEB_SPM==1")
       end if
       if(isource==1)then
          
          if (associated(oneGaspartFields%pnot)) then
             call InsertAtScalarTab(oneGaspartFields%pno,oneGaspartFields%pnot, 'PNO', &
                  oneScalarTab, oneScalarTabSize)
          endif

          if (associated(oneGaspartFields%pno2t)) then
             call InsertAtScalarTab(oneGaspartFields%pno2,oneGaspartFields%pno2t, 'PNO2', &
                  oneScalarTab, oneScalarTabSize)
          endif

          if (associated(oneGaspartFields%ppm25t)) then
             call InsertAtScalarTab(oneGaspartFields%ppm25,oneGaspartFields%ppm25t, 'PM25', &
                  oneScalarTab, oneScalarTabSize)
          endif

          if (associated(oneGaspartFields%pcot)) then
             call InsertAtScalarTab(oneGaspartFields%pco,oneGaspartFields%pcot, 'PCO', &
                  oneScalarTab, oneScalarTabSize)
          endif

          if (associated(oneGaspartFields%pso2t)) then
             call InsertAtScalarTab(oneGaspartFields%pso2,oneGaspartFields%pso2t, 'PSO2', &
                  oneScalarTab, oneScalarTabSize)
          endif

          if (associated(oneGaspartFields%pso4t)) then
             call InsertAtScalarTab(oneGaspartFields%pso4,oneGaspartFields%pso4t, 'PSO4', &
                  oneScalarTab, oneScalarTabSize)
          endif

          if (associated(oneGaspartFields%paert)) then
             call InsertAtScalarTab(oneGaspartFields%paer,oneGaspartFields%paert, 'PAER', &
                  oneScalarTab, oneScalarTabSize)
          endif

          if (associated(oneGaspartFields%pvoct)) then
             call InsertAtScalarTab(oneGaspartFields%pvoc,oneGaspartFields%pvoct, 'PVOC', &
                  oneScalarTab, oneScalarTabSize)
          endif

          if(ichemi==1) then
             if (associated(oneGaspartFields%po3t)) then
                call InsertAtScalarTab(oneGaspartFields%po3,oneGaspartFields%po3t, 'PO3', &
                     oneScalarTab, oneScalarTabSize)
             endif

             if (associated(oneGaspartFields%prhcot)) then
                call InsertAtScalarTab(oneGaspartFields%prhco,oneGaspartFields%prhcot, 'PRHCO', &
                     oneScalarTab, oneScalarTabSize)
             endif

             if (associated(oneGaspartFields%pho2t)) then
                call InsertAtScalarTab(oneGaspartFields%pho2,oneGaspartFields%pho2t, 'PHO2', &
                     oneScalarTab, oneScalarTabSize)
             endif

             if (associated(oneGaspartFields%po3pt)) then
                call InsertAtScalarTab(oneGaspartFields%po3p,oneGaspartFields%po3pt, 'PO3P', &
                     oneScalarTab, oneScalarTabSize)
             endif

             if (associated(oneGaspartFields%po1dt)) then
                call InsertAtScalarTab(oneGaspartFields%po1d,oneGaspartFields%po1dt, 'PO1D', &
                     oneScalarTab, oneScalarTabSize)
             endif

             if (associated(oneGaspartFields%phot)) then
                call InsertAtScalarTab(oneGaspartFields%pho,oneGaspartFields%phot, 'PHO', &
                     oneScalarTab, oneScalarTabSize)
             endif

             if (associated(oneGaspartFields%proot)) then
                call InsertAtScalarTab(oneGaspartFields%proo,oneGaspartFields%proot, 'PROO', &
                     oneScalarTab, oneScalarTabSize)
             endif
          endif

       endif

    endif
    !

    do nsc=1,naddsc
       write(sname,'(a4,i3.3)') 'SCLP',nsc
       if (associated(oneScalarFields(nsc)%sclt)) then
          call InsertAtScalarTab(oneScalarFields(nsc)%sclp,oneScalarFields(nsc)%sclt, sname, &
               oneScalarTab, oneScalarTabSize)
       endif
    enddo

  end subroutine filltab_tend

end module mem_tend
