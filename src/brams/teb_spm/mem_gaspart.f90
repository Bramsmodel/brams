module mem_gaspart

  type gaspart_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, contiguous, pointer :: PCO(:,:,:)
     real, contiguous, pointer :: PNO(:,:,:)
     real, contiguous, pointer :: PNO2(:,:,:)
     real, contiguous, pointer :: PPM25(:,:,:)
     real, contiguous, pointer :: PVOC(:,:,:)
     real, contiguous, pointer :: PSO2(:,:,:)
     real, contiguous, pointer :: PROO(:,:,:)
     real, contiguous, pointer :: PSO4(:,:,:)
     real, contiguous, pointer :: PAER(:,:,:)
     real, contiguous, pointer :: PO3(:,:,:)
     real, contiguous, pointer :: PRHCO(:,:,:)
     real, contiguous, pointer :: PHO2(:,:,:)
     real, contiguous, pointer :: PO3P(:,:,:)
     real, contiguous, pointer :: PO1D(:,:,:)
     real, contiguous, pointer :: PHO(:,:,:)
     real, contiguous, pointer :: GASR(:,:,:)
     real, contiguous, pointer :: PEOXID(:,:,:)

     ! Variables to be dimensioned by (nzp,nxp)
     real, contiguous, pointer :: fusog(:,:)

     real, contiguous, pointer :: PCOT(:)
     real, contiguous, pointer :: PNOT(:)
     real, contiguous, pointer :: PNO2T(:)
     real, contiguous, pointer :: PPM25T(:)
     real, contiguous, pointer :: PVOCT(:)
     real, contiguous, pointer :: PSO2T(:)
     real, contiguous, pointer :: PSO4T(:)
     real, contiguous, pointer :: PAERT(:) 
     real, contiguous, pointer :: PO3T(:)
     real, contiguous, pointer :: PRHCOT(:)
     real, contiguous, pointer :: PHO2T(:)
     real, contiguous, pointer :: PO3PT(:)
     real, contiguous, pointer :: PO1DT(:)
     real, contiguous, pointer :: PHOT(:)
     real, contiguous, pointer :: PROOT(:)

  end type gaspart_vars

  type (gaspart_vars), allocatable, target :: gaspart_g(:), gaspartm_g(:)

contains

  subroutine alloc_gaspart(gaspart, n1, n2, n3)

    use mem_emiss, only: ichemi

    implicit none
    ! Arguments:
    type (gaspart_vars), intent(INOUT) :: gaspart
    integer, intent(in) :: n1, n2, n3

    allocate (gaspart%pco(n1,n2,n3))
    allocate (gaspart%pno(n1,n2,n3))
    allocate (gaspart%pno2(n1,n2,n3))
    allocate (gaspart%ppm25(n1,n2,n3))
    allocate (gaspart%pso2(n1,n2,n3))
    allocate (gaspart%pvoc(n1,n2,n3))
    allocate (gaspart%gasr(n1,n2,n3))
    allocate (gaspart%pso4(n1,n2,n3))
    allocate (gaspart%paer(n1,n2,n3))
    allocate (gaspart%PEOXID(n1,n2,n3))
    allocate (gaspart%fusog(n2,n3))

    if (ichemi==1) then
       allocate (gaspart%po3   (n1,n2,n3))
       allocate (gaspart%prhco (n1,n2,n3))
       allocate (gaspart%pho2  (n1,n2,n3))
       allocate (gaspart%po3p  (n1,n2,n3))
       allocate (gaspart%po1d  (n1,n2,n3))
       allocate (gaspart%pho   (n1,n2,n3))
       allocate (gaspart%proo  (n1,n2,n3))
    endif

  end subroutine alloc_gaspart



  subroutine zero_gaspart(gaspart, n1, n2, n3)

    use mem_emiss, only: ichemi

    implicit none
    ! Arguments:
    type (gaspart_vars), intent(OUT) :: gaspart
    integer, intent(in) :: n1, n2, n3

    gaspart%pco    = 0.0
    gaspart%pno    = 0.0
    gaspart%pno2   = 0.0
    gaspart%ppm25  = 0.0
    gaspart%pso2   = 0.0
    gaspart%pvoc   = 0.0
    gaspart%gasr   = 0.0
    gaspart%pso4   = 0.0
    gaspart%paer   = 0.0
    gaspart%PEOXID = 0.0
    gaspart%fusog  = 0.0

    if (ichemi==1) then
       gaspart%po3   = 0.0
       gaspart%prhco = 0.0
       gaspart%pho2  = 0.0
       gaspart%po3p  = 0.0
       gaspart%po1d  = 0.0
       gaspart%pho   = 0.0
       gaspart%proo  = 0.0
    endif

  end subroutine zero_gaspart


  subroutine nullify_gaspart(gaspart)
    use mem_emiss, only: ichemi

    implicit none
    ! Arguments:
    type (gaspart_vars), intent(INOUT) :: gaspart

    if (associated(gaspart%pco))    nullify (gaspart%pco)
    if (associated(gaspart%pno))    nullify (gaspart%pno)
    if (associated(gaspart%pno2))   nullify (gaspart%pno2)
    if (associated(gaspart%ppm25))  nullify (gaspart%ppm25)
    if (associated(gaspart%pvoc))   nullify (gaspart%pvoc)
    if (associated(gaspart%pso2))   nullify (gaspart%pso2)
    if (associated(gaspart%pso4))   nullify (gaspart%pso4)
    if (associated(gaspart%paer))   nullify (gaspart%paer)
    if (associated(gaspart%PEOXID))   nullify (gaspart%PEOXID)
    if (associated(gaspart%gasr))   nullify (gaspart%gasr)
    if (associated(gaspart%fusog))   nullify (gaspart%fusog)

    if (ichemi==1) then
       if (associated(gaspart%po3  ))  nullify (gaspart%po3  )
       if (associated(gaspart%prhco))  nullify (gaspart%prhco)
       if (associated(gaspart%pho2 ))  nullify (gaspart%pho2 )
       if (associated(gaspart%po3p ))  nullify (gaspart%po3p )
       if (associated(gaspart%po1d ))  nullify (gaspart%po1d )
       if (associated(gaspart%pho  ))  nullify (gaspart%pho  )
       if (associated(gaspart%proo ))  nullify (gaspart%proo )
    endif

  end subroutine nullify_gaspart



  subroutine dealloc_gaspart(gaspart)

    use mem_emiss, only: ichemi

    implicit none

    ! Arguments
    type (gaspart_vars), intent(INOUT) :: gaspart

    if (associated(gaspart%pco))    deallocate (gaspart%pco)
    if (associated(gaspart%pno))    deallocate (gaspart%pno)
    if (associated(gaspart%pno2))   deallocate (gaspart%pno2)
    if (associated(gaspart%ppm25))  deallocate (gaspart%ppm25)
    if (associated(gaspart%pvoc))   deallocate (gaspart%pvoc)
    if (associated(gaspart%pso2))   deallocate (gaspart%pso2)
    if (associated(gaspart%pso4))   deallocate (gaspart%pso4)
    if (associated(gaspart%paer))   deallocate (gaspart%paer)
    if (associated(gaspart%PEOXID))   deallocate (gaspart%PEOXID)
    if (associated(gaspart%gasr))   deallocate (gaspart%gasr)
    if (associated(gaspart%fusog))   deallocate (gaspart%fusog)

    if (ichemi==1) then
       if (associated(gaspart%po3  ))  deallocate (gaspart%po3  )
       if (associated(gaspart%prhco))  deallocate (gaspart%prhco)
       if (associated(gaspart%pho2 ))  deallocate (gaspart%pho2 )
       if (associated(gaspart%po3p ))  deallocate (gaspart%po3p )
       if (associated(gaspart%po1d ))  deallocate (gaspart%po1d )
       if (associated(gaspart%pho  ))  deallocate (gaspart%pho  )
       if (associated(gaspart%proo ))  deallocate (gaspart%proo )
    endif

  end subroutine dealloc_gaspart


  subroutine filltab_gaspart(gaspart, gaspartm, imean, n1, n2, n3, ng)

    use mem_emiss, only: ichemi
    use ModVarTables, only: &
         InsertVtab

    implicit none
    include "constants.h"
    ! Arguments:
    type (gaspart_vars), intent(IN) :: gaspart, gaspartm
    integer, intent(IN) :: imean, n1, n2, n3, ng
    ! Local Variables:
    integer(kind=i8)  :: npts
    real, pointer :: var, varm

    ! Fill pointers to arrays into variable tables

    npts=n2*n3

    if (associated(gaspart%fusog))  &
         call InsertVTab (gaspart%FUSOG,gaspartm%FUSOG &
         ,ng, npts, imean,  &
         'FUSOG:2:hist:anal:mpti:mpt3:mpt1')

    npts=n1*n2*n3

    if (associated(gaspart%pco))  &
         call InsertVTab (gaspart%PCO,gaspartm%PCO&
         ,ng, npts, imean,  &
         'PCO:3:hist:anal:mpti:mpt3:mpt1')
    if (associated(gaspart%pno))  &
         call InsertVTab (gaspart%PNO,gaspartm%PNO&
         ,ng, npts, imean,  &
         'PNO:3:hist:anal:mpti:mpt3:mpt1')
    if (associated(gaspart%pno2))  &
         call InsertVTab (gaspart%PNO2,gaspartm%PNO2&
         ,ng, npts, imean,  &
         'PNO2:3:hist:anal:mpti:mpt3:mpt1')
    if (associated(gaspart%ppm25))  &
         call InsertVTab (gaspart%PPM25,gaspartm%PPM25&
         ,ng, npts, imean,  &
         'PPM25:3:hist:anal:mpti:mpt3:mpt1')
    if (associated(gaspart%pvoc))  &
         call InsertVTab (gaspart%PVOC,gaspartm%PVOC&
         ,ng, npts, imean,  &
         'PVOC:3:hist:anal:mpti:mpt3:mpt1')
    if (associated(gaspart%pso2))  &
         call InsertVTab (gaspart%PSO2,gaspartm%PSO2&
         ,ng, npts, imean,  &
         'PSO2:3:hist:anal:mpti:mpt3:mpt1')
    if (associated(gaspart%pso4))  &
         call InsertVTab (gaspart%PSO4,gaspartm%PSO4&
         ,ng, npts, imean,  &
         'PSO4:3:hist:anal:mpti:mpt3:mpt1')
    if (associated(gaspart%paer))  &
         call InsertVTab (gaspart%PAER,gaspartm%PAER&
         ,ng, npts, imean,  &
         'PAER:3:hist:anal:mpti:mpt3:mpt1')
    if (associated(gaspart%PEOXID))  &
         call InsertVTab (gaspart%PEOXID,gaspartm%PEOXID&
         ,ng, npts, imean,  &
         'PEOXID:3:hist:anal:mpti:mpt3:mpt1')
    if (associated(gaspart%gasr))  &
         call InsertVTab (gaspart%GASR,gaspartm%GASR&
         ,ng, npts, imean,  &
         'GASR:3:mpti:mpt3:mpt1')

    if (ichemi==1) then
       if (associated(gaspart%po3))  &
            call InsertVTab (gaspart%PO3,gaspartm%PO3&
            ,ng, npts, imean,  &
            'PO3:3:hist:anal:mpti:mpt3:mpt1')
       if (associated(gaspart%prhco))  &
            call InsertVTab (gaspart%PRHCO,gaspartm%PRHCO&
            ,ng, npts, imean,  &
            'PRHCO:3:hist:anal:mpti:mpt3:mpt1')
       if (associated(gaspart%pho2))  &
            call InsertVTab (gaspart%PHO2,gaspartm%PHO2&
            ,ng, npts, imean,  &
            'PHO2:3:hist:anal:mpti:mpt3:mpt1')
       if (associated(gaspart%po3p))  &
            call InsertVTab (gaspart%PO3P,gaspartm%PO3P&
            ,ng, npts, imean,  &
            'PO3P:3:hist:anal:mpti:mpt3:mpt1')
       if (associated(gaspart%po1d))  &
            call InsertVTab (gaspart%PO1D,gaspartm%PO1D&
            ,ng, npts, imean,  &
            'PO1D:3:hist:anal:mpti:mpt3:mpt1')
       if (associated(gaspart%pho))  &
            call InsertVTab (gaspart%PHO,gaspartm%PHO&
            ,ng, npts, imean,  &
            'PHO:3:hist:anal:mpti:mpt3:mpt1')
       if (associated(gaspart%proo))  &
            call InsertVTab (gaspart%PROO,gaspartm%PROO&
            ,ng, npts, imean,  &
            'PROO:3:hist:anal:mpti:mpt3:mpt1')
    endif

  end subroutine filltab_gaspart

end module mem_gaspart
