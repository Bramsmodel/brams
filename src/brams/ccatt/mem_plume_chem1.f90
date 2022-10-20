!###########################################################################
!  B - Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_plume_chem1

  use ModNamelistFile, only: namelistFile

  include "constants.h"

  type plume_vars   
     real, pointer, contiguous :: fct(:,:)
  end type plume_vars
  
  type(plume_vars), allocatable, target :: plume_g (:,:,:)
  type(plume_vars), allocatable, target :: plumem_g (:,:,:)

  type plume_mean_vars   
     real, pointer, contiguous :: fire_size(:,:)
     real, pointer, contiguous :: flam_frac(:,:)
  end type plume_mean_vars

  type(plume_mean_vars), allocatable, target :: plume_mean_g(:,:)
  type(plume_mean_vars), allocatable, target :: plume_meanm_g(:,:)


  integer, parameter :: nveg_agreg      = 4
  integer, parameter :: tropical_forest = 1
  integer, parameter :: boreal_forest   = 2
  integer, parameter :: savannah        = 3
  integer, parameter :: grassland       = 4 ! must be equal to nveg_agreg
  character(len=20), parameter :: veg_name(nveg_agreg) = (/ &
                               'Tropical-Forest', &
                               'Boreal-Forest  ', &
                               'Savanna        ', &
                               'Grassland      ' /)
  character(len=20), parameter :: spc_suf(nveg_agreg) = (/ &
                               'agtf' , &  ! trop forest
                               'agef' , &  ! extratrop forest
                               'agsv' , &  ! savanna
                               'aggr'   /) ! grassland
  real                         :: prfrq 	    
  integer                      :: plumerise 


 !-- this section is for the FRP methodology

  type plume_fre_vars   
     real, pointer, contiguous :: pvar(:,:)
  end type plume_fre_vars
  
  type(plume_fre_vars), allocatable, target :: plume_fre_g (:,:)
  type(plume_fre_vars), allocatable, target :: plumem_fre_g (:,:)
 
  integer, parameter ::      &
           iflam_frac  =1  &
          ,imean_frp   =2  &
          ,istd_frp    =3  &
          ,imean_size  =4  &
          ,istd_size   =5  

  character(len=10), parameter :: fre_var_name(5) = (/ &
                               'flam_frac ' , &  ! 
                               'mean_frp  ' , &  ! 
                               'std_frp   ' , &  ! 
                               'mean_size ' , &  ! 
                               'std_size  '   /) ! 




contains
  !---------------------------------------------------------------

  subroutine alloc_plume_chem1(plume,plume_mean,plume_fre,n1,n2,n3,nspecies )
                              
    use chem1_list, only : spc_alloc,spc_name,src,on
    implicit none

    integer,intent(in) :: n1,n2,n3,nspecies
    integer :: ispc,iv
    
    type (plume_vars)      ,dimension(nveg_agreg,nspecies) :: plume
    type (plume_mean_vars) ,dimension(nveg_agreg         ) :: plume_mean

    type (plume_fre_vars)  ,dimension(5                  ) :: plume_fre

    integer:: imean_plume
    imean_plume = 1 !change this at emis_flam_smold routine also
         
    !print*,'----------------------------------------------------------------'
    !print*,' memory allocation for plumerise sources:'
    if(PLUMERISE == 1) then
      if(imean_plume /= on) then
    
    	  do ispc=1,nspecies
    	   !print*,'spc=',spc_name(ispc),'size=',n1,n2,n3

    	   if(spc_alloc(src,ispc) == on) then 

    	      do iv=1,nveg_agreg
    		  allocate (plume(iv,ispc)%fct(n2,n3))
    		  plume(iv,ispc)%fct(:,:) = 0.
    	      enddo
    	   endif
          enddo
       else 
 
          do iv=1,nveg_agreg
            allocate (plume_mean(iv)%flam_frac(n2,n3))
            plume_mean(iv)%flam_frac(:,:)=0.
          enddo
      endif


      do iv=1,nveg_agreg
        allocate (plume_mean(iv)%fire_size(n2,n3))
        plume_mean(iv)%fire_size(:,:)=0.
      enddo
    endif
    !- for FRP methodology
    if(PLUMERISE == 2)then
      do iv=1,5
         allocate (plume_fre(iv)%pvar(n2,n3))
         plume_fre(iv)%pvar(:,:)=0.
      enddo
    
    endif

    return
  end subroutine alloc_plume_chem1

  !---------------------------------------------------------------

  subroutine dealloc_plume_chem1(plume,plume_mean,plume_fre,nspecies)

   implicit none

   integer,intent(in) :: nspecies
   type (plume_vars)    ,dimension(nveg_agreg,nspecies) :: plume
   type (plume_mean_vars) ,dimension(nveg_agreg         ) :: plume_mean
   type (plume_fre_vars)  ,dimension(5                  ) :: plume_fre
   integer :: ispc,iv
  
    !  Deallocate arrays
    do ispc=1,nspecies
     do iv=1,nveg_agreg
       if (associated(plume(iv,ispc)%fct)) deallocate(plume(iv,ispc)%fct)
     enddo
    enddo
    do iv=1,nveg_agreg
       if (associated(plume_mean(iv)%fire_size)) deallocate(plume_mean(iv)%fire_size)
       if (associated(plume_mean(iv)%flam_frac)) deallocate(plume_mean(iv)%flam_frac)
    enddo
     
    do iv=1,5
       if (associated(plume_fre(iv)%pvar)) deallocate(plume_fre(iv)%pvar)
    enddo
    
     
 end subroutine dealloc_plume_chem1

  !---------------------------------------------------------------

  subroutine nullify_plume_chem1(plume,plume_mean,plume_fre,nspecies)

    implicit none

    integer,intent(in) :: nspecies
    type (plume_vars)    ,dimension(nveg_agreg,nspecies) :: plume
    type (plume_mean_vars) ,dimension(nveg_agreg         ) :: plume_mean
    type (plume_fre_vars)  ,dimension(5                  ) :: plume_fre
    
    integer :: ispc,iv
 
    do ispc=1,nspecies

      do iv=1,nveg_agreg
        if (associated(plume(iv,ispc)%fct)) nullify(plume(iv,ispc)%fct)
      enddo
    
    enddo
    do iv=1,nveg_agreg
       if (associated(plume_mean(iv)%fire_size)) nullify(plume_mean(iv)%fire_size)
       if (associated(plume_mean(iv)%flam_frac)) nullify(plume_mean(iv)%flam_frac)
    enddo
 
    !- this is for FRP method
    do iv=1,5
       if (associated(plume_fre(iv)%pvar)) nullify(plume_fre(iv)%pvar)
    enddo

    return
  end subroutine nullify_plume_chem1

  !---------------------------------------------------------------

  subroutine filltab_plume_chem1(oneVarTable, oneVarTableSize, &
       plume, plumem, plume_mean, plume_meanm, plume_fre, plumem_fre, &
       nspecies, ng, imean)

    use chem1_list, only: &
         spc_alloc, &
         spc_name, &
         src, &
         on

    use mem_chem1, only: &
         chem1_g

    use ModVarTable, only: &
         VarTable, &
         InsertVarTable

    implicit none

    type(VarTable), pointer, intent(in) :: oneVarTable(:)
    integer, intent(inout) :: oneVarTableSize
    type(plume_vars), pointer, intent(in) :: plume(:,:)
    type(plume_vars), pointer, intent(in) :: plumem(:,:)
    type(plume_mean_vars), pointer, intent(in) :: plume_mean(:)
    type(plume_mean_vars), pointer, intent(in) :: plume_meanm(:)
    type(plume_fre_vars), pointer, intent(in) :: plume_fre(:)
    type(plume_fre_vars), pointer, intent(in) :: plumem_fre(:)
    integer, intent(in) :: nspecies
    integer, intent(in) :: ng
    integer, intent(in) :: imean

    logical :: assThis
    integer :: ispc,iv  
    integer :: imean_plume
    character(len=*), parameter :: h="**(filltab_plume_chem1)**"

    if (.not. associated(oneVarTable)) then
       call fatal_error(h//" oneVarTable not associated")
    end if


    imean_plume = 1 !change this at emis_flam_smold routine also


    ! Fill pointers to arrays into variable tables


    if(PLUMERISE == 1) then
       if(imean_plume /= on) then
          do ispc=1,nspecies
             if (associated(chem1_g(ispc,ng)%sc_p)) then
                !------- sources 
                if(spc_alloc(src,ispc) == on) then
                   if (.not. associated(plume)) then
                      call fatal_error(h//" plume not associated but plumerise=1 and "//&
                           "imean_plume/=on and chem1_g as well as spc_alloc are associated")
                   end if
                   do iv=1,nveg_agreg
                      call InsertVarTable(oneVarTable, oneVarTableSize, &
                           plume(iv,ispc)%fct,       &
                           trim(spc_name(ispc))//'_'//trim(spc_suf(iv))//&
                           ' :2:hist:anal:mpti:mpt3:mpt1', &
                           plumem(iv,ispc)%fct, imean)
                   end do
                end if
             end if
          end do
       else
          if (.not. associated(plume_mean)) then
             call fatal_error(h//" plume_mean not associated but plumerise=1 and imean_plume=on")
          end if
          !---  flam frac mean
          do iv=1,nveg_agreg
             call InsertVarTable(oneVarTable, oneVarTableSize, &
                  plume_mean(iv)%flam_frac,  &
                  'flam_frac'//'_'//trim(spc_suf(iv))//&
                  ' :2:hist:anal:mpti:mpt3:mpt1', &
                  plume_meanm(iv)%flam_frac, imean)
          enddo
       endif
       !--- fire size   
       if (.not. associated(plume_mean)) then
          call fatal_error(h//" plume_mean not associated but plumerise=1")
       end if
       do iv=1,nveg_agreg
          call InsertVarTable(oneVarTable, oneVarTableSize, &
               plume_mean(iv)%fire_size,  &
               'firesize'//'_'//trim(spc_suf(iv))//&
               ' :2:hist:anal:mpti:mpt3:mpt1', &
               plume_meanm(iv)%fire_size, imean)
       end do
    else if(PLUMERISE == 2)then
       if (.not. associated(plume_fre)) then
          call fatal_error(h//" plume_fre not associated but plumerise=2")
       end if
       do iv=1,5
          call InsertVarTable(oneVarTable, oneVarTableSize, &
               plume_fre(iv)%pvar, &
               fre_var_name(iv)//' :2:hist:anal:mpti:mpt3:mpt1', &
               plumem_fre(iv)%pvar, imean)
       end do
    end if
  end subroutine filltab_plume_chem1

  !---------------------------------------------------------------

  subroutine StoreNamelistFileAtMem_plumeChem1(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    plumerise  = oneNamelistFile%plumerise
    prfrq = oneNamelistFile%prfrq
  end subroutine StoreNamelistFileAtMem_plumeChem1


  !---------------------------------------------------------------
end module mem_plume_chem1
