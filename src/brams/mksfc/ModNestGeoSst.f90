!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModNestGeoSst

  use ModParallelEnvironment, only: &
       MsgDump
  
  use ModNestFillDens, only: &
       fillscr, &
       fillvar

  use ModNestFeed, only: &
       fdback, &
       fdback
  
  use ModLeaf3Init, only: &
       sfcinit_nofile, &
       sfcinit_file
  
  use ModRUser, only: &
       toptinit_user, &
       sfcinit_file_user, &
       sfcinit_nofile_user
       
  use ModLanduseInput, only: &
       landuse_opqr
  
  use ModGeodat, only: &
       geodat, &
       geodat_var, &
       geodat_var_opt
  
  use ModMkSfcTop, only: &
       toptinit
  
  use ModControlVars, only: &
       ControlVars
  
  use ModInitHis, only: &
       patch_land_average, &
       patch_land_unaverage
  
  use mem_mksfc, only: &
       sfcfile_p, &
       vndvifil, &
       mksfc_scr1, &
       mksfc_scr2, &
       mksfc_vt2da, &
       mksfc_vt2db

  use mem_grid, only: &
       nxtnest, &
       nnxp, &
       nnyp, &
       nzg, &
       npatch, &
       ipm, &
       jpm, &
       platn, &
       plonn, &
       maxnzp, &
       maxnxp, &
       maxnyp, &
       grid_g, &
       nnzp, &
       nzs, &
       maxnxp, &
       maxnyp, &
       ngridsh, &
       jdim, &
       nstratx, &
       nstraty

  use io_params, only: &
       ivegtflg, &
       ivegtfn, &
       isoilflg, &
       isoilfn, &
       ndviflg, &
       ndvifn, &
       ifusflg, &
       ifusfn, &
       nofilflg, &
       itoptflg, &
       itoptfn, &
       toptwvl

  use mem_leaf, only: &
       leaf_g, &
       nvegpat

  use dump, only: &
       dumpMessage

  use ModBasicFields, only: &
       BasicFields

  use ModTurbFields, only: &
       TurbFields

  use node_mod, only: &
       nodemxp, nodemyp, &
       mchnum,        &
       master_num, &
       mynum

  use memSoilMoisture, only: &
       SOIL_MOIST  ! INTENT(IN)

  use ModSoilMoisture, only: &
       soilMoistureInit ! Subroutine

  use ccatt_start, only: &
       ccatt

  use grid_dims, only: &
       nxpmax, &
       nypmax

  implicit none

  private

  public :: toptnest
  public :: geonest_file
  public :: GeonestNoFile

contains

  subroutine toptnest(ngra,ngrb)
    include "constants.h"
    ! Arguments:
    integer, intent(in) :: ngra, ngrb
    ! Local Variables:
    integer :: ifm,icm,ipat,i,j,k,indfm,ivtime,nc1
    integer :: i1, i2, j1, j2
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(toptnest)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if
    
    do ifm=ngra,ngrb
       icm = nxtnest(ifm)
       ! Initialize TOPOGRAPHY in toptinit.

       call toptinit(nnxp(ifm), nnyp(ifm), ifm,   &
            sfcfile_p(ifm)%topt, sfcfile_p(ifm)%topzo)

       if (icm>=1 .and. itoptflg(ifm)==0) then

          ! Interpolate TOPO from coarser grid:

          call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
               ,mksfc_scr1,sfcfile_p(icm)%topt)
!!$        call eintp(mksfc_scr1,mksfc_scr2,1,maxnxp,maxnyp  &
!!$             ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
          call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
               ,mksfc_scr2,sfcfile_p(ifm)%topt)
          if (dumpLocal) then
             call MsgDump(h//" interpolate topt from coarser grid")
          end if
          
          ! Interpolate TOPO ZO from coarser grid:
          call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
               ,mksfc_scr1,sfcfile_p(icm)%topzo)
!!$        call eintp(mksfc_scr1,mksfc_scr2,1,maxnxp,maxnyp  &
!!$             ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
          call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
               ,mksfc_scr2,sfcfile_p(ifm)%topzo)
          if (dumpLocal) then
             call MsgDump(h//" interpolate topzo from coarser grid")
          end if

       elseif (itoptflg(ifm)==1) then

          !-srf: para topografia com suavizacao variavel
          if ( TOPTWVL(ifm)<0.) then
             ! Interpolate TOPO from standard dataset:
             !SRF-OPT   call geodat_var    (nnxp(ifm), nnyp(ifm), sfcfile_p(ifm)%topt(1,1),  &
             call geodat_var_OPT(nnxp(ifm), nnyp(ifm), sfcfile_p(ifm)%topt,  &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  mksfc_vt2da, mksfc_vt2db, ifm, 'TOP')
             ! Interpolate TOPO ZO from standard dataset:
             !SRF-OPT   call geodat_var    (nnxp(ifm), nnyp(ifm), sfcfile_p(ifm)%topzo(1,1),  &
             call geodat_var_OPT(nnxp(ifm), nnyp(ifm), sfcfile_p(ifm)%topzo,  &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  mksfc_vt2da, mksfc_vt2db, ifm, 'ZOT')
             if (dumpLocal) then
                call MsgDump(h//" smooth interpolation of topt and topzo from standard data set")
             end if
          else
             ! Interpolate TOPO from standard dataset:
             call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topt,  &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  mksfc_vt2da,mksfc_vt2db,ifm,'TOP')
             ! Interpolate TOPO ZO from standard dataset:
             call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topzo,  &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  mksfc_vt2da,mksfc_vt2db,ifm,'ZOT')
             if (dumpLocal) then
                call MsgDump(h//" interpolate topt and topzo from standard data set")
             end if
          endif

       elseif (itoptflg(ifm)==3) then

          if (TOPTWVL(ifm)<0.) then
             ! Interpolate TOPO from dted dataset:
             call geodat_var(nnxp(ifm), nnyp(ifm), sfcfile_p(ifm)%topt,  &
                  itoptfn(ifm), itoptfn(ifm), mksfc_vt2da, mksfc_vt2db, ifm, &
                  'TOD')
             ! Interpolate TOPO ZO from dted dataset:
             call geodat_var(nnxp(ifm), nnyp(ifm), sfcfile_p(ifm)%topzo,  &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  mksfc_vt2da, mksfc_vt2db, ifm, 'ZOD')
             if (dumpLocal) then
                call MsgDump(h//" smooth interpolation of topt and topzo from dted data set")
             end if
          else
             ! Interpolate TOPO from dted dataset:
             call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topt,  &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  mksfc_vt2da,mksfc_vt2db,ifm,'TOD')
             ! Interpolate TOPO ZO from dted dataset:
             call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topzo,  &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  itoptfn(ifm)(1:len_trim(itoptfn(ifm))), &
                  mksfc_vt2da,mksfc_vt2db,ifm,'ZOD')
             if (dumpLocal) then
                call MsgDump(h//" interpolate topt and topzo from dted data set")
             end if
          endif

       endif

       ! If desired, override current values of TOPOGRAPHY in ruser.f subroutine.

       call toptinit_user(nnxp(ifm),nnyp(ifm),ifm  &
            ,sfcfile_p(ifm)%topt ,sfcfile_p(ifm)%topzo)

    enddo

    if (ngra .eq. ngrb) then
       if (dumpLocal) then
          call MsgDump(h//" finishes")
       end if
       return
    end if

    ! In case topography data have been independently reassigned on any grid,
    ! average fine mesh topography sequentially to the coarser grids.

    do ifm = ngrb,ngra,-1
       if (nxtnest(ifm) .gt. ngridsh .and. ifm .ge. 2) then
          icm = nxtnest(ifm)
          call fdback(sfcfile_p(icm)%topt,sfcfile_p(ifm)%topt  &
               ,mksfc_vt2da,mksfc_scr2,1,nnxp(icm),nnyp(icm)  &
               ,1,nnxp(ifm),nnyp(ifm),ifm,'terr',mksfc_vt2db)

          if (dumpLocal) then
             write(str(1),"(i8)") icm
             write(str(2),"(i8)") ifm
             call MsgDump(h//" topt of coarser grid "//trim(adjustl(str(1)))//&
                  " rebuild due to reassignment at finer grid "//trim(adjustl(str(2))))
          end if

       endif
    enddo

    ! In case terrain heights have been independently reassigned on
    ! any grid, interpolate coarse grid terrain heights to a temporary
    ! fine mesh array.  Fill the fine mesh boundary terrain heights
    ! from the temporary array.

    do ifm = ngra,ngrb
       icm = nxtnest(ifm)
       if (icm .ge. 1) then
          call fillscr(1,nxpmax,nypmax,1,nnxp(icm),nnyp(icm),1,1  &
               ,mksfc_scr1,sfcfile_p(icm)%topt)

!!$        call eintp(mksfc_scr1,mksfc_scr2,1,nxpmax,nypmax  &
!!$             ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
          call fillvar(1,nxpmax,nypmax,1,nnxp(ifm),nnyp(ifm),1,1  &
               ,mksfc_scr2,mksfc_scr1)

          nc1 = jdim * (nstraty(ifm) + 1)
!!$        call ae2(nnxp(ifm),nnyp(ifm),2+nstratx(ifm)  &
!!$             ,nnxp(ifm)-1-nstratx(ifm),1+nc1,nnyp(ifm)-nc1  &
!!$             ,mksfc_scr1,sfcfile_p(ifm)%topt(1,1))
          i1 = 2+nstratx(ifm)
          i2 = nnxp(ifm)-1-nstratx(ifm)
          j1 = 1+nc1
          j2 = nnyp(ifm)-nc1
          mksfc_scr1(i1:i2,j1:j2) = sfcfile_p(ifm)%topt(i1:i2,j1:j2)

!!$        call ae1(nnxp(ifm)*nnyp(ifm),sfcfile_p(ifm)%topt(1,1),mksfc_scr1)
          call ae1_l(int(nnxp(ifm)*nnyp(ifm),i8), &
               sfcfile_p(ifm)%topt, mksfc_scr1)

          if (dumpLocal) then
             write(str(1),"(i8)") icm
             write(str(2),"(i8)") ifm
             call MsgDump(h//" fill boundary of topt at finer grid "//trim(adjustl(str(2)))//&
                  " from coarser grid "//trim(adjustl(str(1))))
          end if

       endif
    enddo

  end subroutine toptnest

  !*************************************************************************

  subroutine geonest_file(ifm)
    integer, intent(IN) :: ifm

    integer :: icm,ipat,i,j,k,indfm,ivtime,nc1,mynum
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(geonest_file)**"
    logical, parameter :: dumpLocal=.false.

    icm = nxtnest(ifm)

    if (dumpLocal) then
       write(str(1),"(i8)") ifm
       call MsgDump(h//" initialize patch_area, leaf_class, soil_text of grid "//&
            trim(adjustl(str(1))))
    end if
    
    ! Initialize PATCH AREA, LANDUSE CLASS, and SOIL TEXTURAL CLASS
    ! in subroutine sfcinit.

    call sfcinit_file(nnxp(ifm),nnyp(ifm),nzg,npatch,ifm      &
         ,sfcfile_p(ifm)%patch_area   &
         ,sfcfile_p(ifm)%leaf_class   &
         ,sfcfile_p(ifm)%soil_text)


    print*, ' '
    print*,'====================================================='
    print*,'Starting landuse data input on grid ',ifm
    print*,'====================================================='

    if (icm .ge. 1 .and. ivegtflg(ifm) .eq. 0) then

       ! Assign PATCH AREAS and PATCH CLASSES from coarser grid:

       if (dumpLocal) then
          write(str(1),"(i8)") ifm
          write(str(2),"(i8)") icm
          call MsgDump(h//" copy patch_area and leaf_class from coarser grid "//&
               trim(adjustl(str(2)))//" into finer grid "//trim(adjustl(str(1))))
       end if

       do ipat = 1,npatch
          do j = 1,nnyp(ifm)
             do i = 1,nnxp(ifm)

                sfcfile_p(ifm)%patch_area(i,j,ipat) =  &
                     sfcfile_p(icm)%patch_area(ipm(i,ifm),jpm(j,ifm),ipat)
                sfcfile_p(ifm)%leaf_class(i,j,ipat) =  &
                     sfcfile_p(icm)%leaf_class(ipm(i,ifm),jpm(j,ifm),ipat)

             enddo
          enddo
       enddo


    elseif (ivegtflg(ifm) .eq. 1) then

       ! Assign PATCH AREAS and PATCH CLASSES from standard dataset:

       if (dumpLocal) then
          write(str(1),"(i8)") ifm
          write(str(2),"(i8)") icm
          call MsgDump(h//" copy patch_area and leaf_class of finer grid "//&
               trim(adjustl(str(1)))//" from standard dataset")
       end if

       call landuse_opqr(nnxp(ifm),nnyp(ifm),nzg,npatch,nvegpat  &
            ,ivegtflg(ifm),ivegtfn(ifm),isoilflg(ifm),isoilfn(ifm) &
            ,ndviflg(ifm),ndvifn(ifm),vndvifil(1,ifm)  &
            ,'veg',platn(ifm),plonn(ifm)        &
            ,sfcfile_p(ifm)%soil_text  &
            ,sfcfile_p(ifm)%patch_area   &
            ,sfcfile_p(ifm)%leaf_class   &
            ,sfcfile_p(ifm)%veg_ndvif)


    endif

    if (icm .ge. 1 .and. isoilflg(ifm) .eq. 0) then

       ! Assign SOIL TEXTURE CLASS from coarser grid

       if (dumpLocal) then
          write(str(1),"(i8)") ifm
          write(str(2),"(i8)") icm
          call MsgDump(h//" copy soil_text from coarser grid "//&
               trim(adjustl(str(2)))//" into finer grid "//trim(adjustl(str(1))))
       end if

       do ipat = 2,npatch
          do k = 1,nzg
             do j = 1,nnyp(ifm)
                do i = 1,nnxp(ifm)
                   sfcfile_p(ifm)%soil_text(k,i,j,ipat) =  &
                        sfcfile_p(icm)%soil_text(k,ipm(i,ifm),jpm(j,ifm),ipat)
                enddo
             enddo
          enddo
       enddo

    elseif (isoilflg(ifm) .eq. 1) then

       ! Assign SOIL TEXTURE CLASS from standard dataset:

       if (dumpLocal) then
          write(str(1),"(i8)") ifm
          call MsgDump(h//" copy soil_text, patch_area, leaf_class and veg_ndvif of finer grid "//&
               trim(adjustl(str(1)))//" from standard dataset")
       end if

       call landuse_opqr(nnxp(ifm),nnyp(ifm),nzg,npatch,nvegpat  &
            ,ivegtflg(ifm),ivegtfn(ifm),isoilflg(ifm),isoilfn(ifm) &
            ,ndviflg(ifm),ndvifn(ifm),vndvifil(1,ifm)  &
            ,'soil',platn(ifm),plonn(ifm)        &
            ,sfcfile_p(ifm)%soil_text  &
            ,sfcfile_p(ifm)%patch_area   &
            ,sfcfile_p(ifm)%leaf_class   &
            ,sfcfile_p(ifm)%veg_ndvif)

    endif

    ! If desired, override current values of PATCH AREA, PATCH CLASS,
    ! LEAF-2 VEGETATION CLASS, SOIL TEXTURAL CLASS, and/or
    ! NDVI in ruser.f subroutines.

    if (dumpLocal) then
       write(str(1),"(i8)") ifm
       call MsgDump(h//" overwrite soil_text, patch_area and leaf_class of finer grid "//&
            trim(adjustl(str(1)))//" from user given dataset")
    end if

    call sfcinit_file_user(nnxp(ifm),nnyp(ifm),nzg,npatch,ifm &
         ,sfcfile_p(ifm)%patch_area      &
         ,sfcfile_p(ifm)%leaf_class      &
         ,sfcfile_p(ifm)%soil_text    )

    ! As a final initialization step, eliminate any land patch area that is less
    ! than 1% of the total grid cell area.  Set its area to zero, and compensate
    ! by enlarging areas of remaining patches.

    if (dumpLocal) then
       write(str(1),"(i8)") ifm
       call MsgDump(h//" eliminate too small  patch_area cells of finer grid "//&
            trim(adjustl(str(1))))
    end if

    call patch_minsize(nnxp(ifm),nnyp(ifm),npatch  &
         ,sfcfile_p(ifm)%patch_area)


    return
  end subroutine geonest_file

  !******************************************************************************

  subroutine patch_interp(icm,ifm,nc1,nc2,nc3,nc4,nf1,nf2,nf3,nf4 &
       ,ac,af,pareac,pareaf,avgc,avgf,slabc,slabf)
    integer :: icm,ifm,nc1,nc2,nc3,nc4,nf1,nf2,nf3,nf4
    real :: ac(nc1,nc2,nc3,nc4), af(nf1,nf2,nf3,nf4)
    real :: pareac(nc2,nc3,nc4), pareaf(nf2,nf3,nf4)
    real :: avgc(nc1,nc2,nc3),   avgf(nf1,nf2,nf3)
    real :: slabc(nc2,nc3), slabf(nf2,nf3)

    integer :: k,i,j
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(patch_interp)**"
    logical, parameter :: dumpLocal=.false.
    
    ! Average coarse grid field over all land patches

    if (dumpLocal) then
       write(str(1),"(i8)") icm
       write(str(2),"(i8)") ifm
       call MsgDump(h//" interpolate finer grid "//trim(adjustl(str(2)))// &
            " field with coarser grid "//trim(adjustl(str(1)))//&
            " field values average over patches")
    end if
    
    call patch_land_average(nc1,nc2,nc3,nc4  &
         ,pareac,ac,avgc)

    ! Interpolate patch-averaged to fine grid

    do k=1,nc1    ! nc1 and nf1 are the same

       do j=1,nc3
          do i=1,nc2
             slabc(i,j)=avgc(k,i,j)
          enddo
       enddo

!!$     call fmint2d(icm,ifm,'t',slabc,slabf)

       do j=1,nf3
          do i=1,nf2
             avgf(k,i,j)=slabf(i,j)
          enddo
       enddo

    enddo

    ! Fill fine grid field back into all land patches

    call patch_land_unaverage(nf1,nf2,nf3,nf4,avgf,af)

    return
  end subroutine patch_interp


  !     *****************************************************************

!!$subroutine fmint5(var1,var2,dn0xc,dn0xf,vt2da,ifm,icm,vpnt,idwt)
!!$
!!$  use mem_scratch, only: &
!!$       scratch
!!$
!!$  use mem_grid, only: &
!!$       maxnzp, maxnxp, maxnyp, nnzp, nnxp, nnyp, grid_g
!!$
!!$  implicit none
!!$  integer :: ifm,icm,idwt
!!$
!!$  real, dimension(*) :: var1,var2,vt2da,dn0xc,dn0xf
!!$  character(len=*) :: vpnt
!!$
!!$  if (icm .eq. 0) return
!!$
!!$  call fillscr(maxnzp,maxnxp,maxnyp,nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$       ,1,nnzp(icm),scratch%scr1(1),var1)
!!$
!!$
!!$  if (idwt .eq. 1) then
!!$     call dnswt2(maxnzp,maxnxp,maxnyp,nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$          ,scratch%scr1(1),dn0xc,vpnt,1)
!!$  endif
!!$
!!$  call eintp(scratch%scr1(1),scratch%scr2(1),maxnzp,maxnxp,maxnyp  &
!!$       ,nnzp(ifm),nnxp(ifm),nnyp(ifm),ifm,3,vpnt,0,0)
!!$
!!$  call fillvar(maxnzp,maxnxp,maxnyp,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
!!$       ,1,nnzp(ifm),scratch%scr2(1),var2)
!!$
!!$  if (idwt .eq. 1) call dnswt2(nnzp(ifm),nnxp(ifm),nnyp(ifm),nnzp(ifm)  &
!!$       ,nnxp(ifm),nnyp(ifm),var2,dn0xf,vpnt,2)
!!$
!!$  call rtgintrp(nnzp(ifm),nnxp(ifm),nnyp(ifm),var2,vt2da  &
!!$       ,grid_g(ifm)%topt(1,1),ifm,vpnt)
!!$
!!$  return
!!$end subroutine fmint5


  !******************************************************************************

  subroutine patch_minsize(n2,n3,npat,patch_area)
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: npat
    real, intent(inout) :: patch_area(n2,n3,npat)

    integer :: i,j,ipat,jpat
    real :: orig_size
    character(len=*), parameter :: h="**(patch_minsize)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" enlarge too small patch areas")
    end if
    
    do j = 1,n3
       do i = 1,n2
          do ipat = 2,npat
             if (patch_area(i,j,ipat) .gt. 0. .and.  &
                  patch_area(i,j,ipat) .lt. .01) then

                orig_size = patch_area(i,j,ipat)
                patch_area(i,j,ipat) = 0.

                do jpat = 1,npat
                   if (jpat .ne. ipat) then
                      patch_area(i,j,jpat) = patch_area(i,j,jpat)  &
                           / (1. - orig_size)
                   endif
                enddo
             endif
          enddo
       enddo
    enddo

    return
  end subroutine patch_minsize


  !*************************************************************************

  subroutine GeonestNofile(ngra, ngrb, &
       oneControlVars, oneBasicFields, oneTurbFields)
    include "constants.h"
    integer, intent(IN) :: ngra, ngrb
    type(ControlVars), pointer, intent(in) :: oneControlVars
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields

    integer :: isiz,ifm,icm,ipat,i,j,k,indfm,ivtime,nc1,ic,jc

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(GeonestNofile)**"
    logical, parameter :: dumpLocal=.true.
    
    real :: vt2da(maxnxp * maxnyp)
    real :: vt2db(maxnxp * maxnyp)
    real :: vt3db(maxnzp * maxnxp * maxnyp)
    real :: vt3da(maxnzp * maxnxp * maxnyp)

    ! Initialization/interpolation of leaf-2 variables for which standard RAMS
    ! datasets never exist.

    isiz = maxnxp * maxnyp

    do ifm = ngra,ngrb
       icm = nxtnest(ifm)

       if (dumpLocal) then
          write(str(1),"(i8)") ifm
          call MsgDump(h//" assign default values for leaf_g variables of grid "//trim(adjustl(str(1))))
       end if

       ! First, fill NOFILE LEAF-2 variables with default values in SFCINIT.

       if (dumpLocal) then
          call MsgDump(h//" invoca sfcinit_nofile")
       end if
       
       call sfcinit_nofile(nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm),nzg,nzs,npatch,ifm    &
            ,oneBasicFields%theta              &
            ,oneBasicFields%pi0                &
            ,oneBasicFields%pp                 &
            ,oneBasicFields%rv                 &
            ,leaf_g(ifm)%seatp, leaf_g(ifm)%seatf      &
            ,leaf_g(ifm)%soil_water        &
            ,leaf_g(ifm)%soil_energy       &
            ,leaf_g(ifm)%soil_text         &
            ,leaf_g(ifm)%sfcwater_mass     &
            ,leaf_g(ifm)%sfcwater_energy   &
            ,leaf_g(ifm)%sfcwater_depth    &
            ,leaf_g(ifm)%ustar             &
            ,leaf_g(ifm)%tstar             &
            ,leaf_g(ifm)%rstar             &
            ,leaf_g(ifm)%veg_fracarea      &
            ,leaf_g(ifm)%veg_lai           &
            ,leaf_g(ifm)%veg_tai           &
            ,leaf_g(ifm)%veg_rough         &
            ,leaf_g(ifm)%veg_height        &
            ,leaf_g(ifm)%veg_albedo        &
            ,leaf_g(ifm)%patch_area        &
            ,leaf_g(ifm)%patch_rough       &
            ,leaf_g(ifm)%patch_wetind      &
            ,leaf_g(ifm)%leaf_class        &
            ,leaf_g(ifm)%soil_rough        &
            ,leaf_g(ifm)%sfcwater_nlev     &
            ,leaf_g(ifm)%stom_resist       &
            ,leaf_g(ifm)%ground_rsat, leaf_g(ifm)%ground_rvap    &
            ,leaf_g(ifm)%veg_water, leaf_g(ifm)%veg_temp    &
            ,leaf_g(ifm)%can_rvap, leaf_g(ifm)%can_temp     &
            ,leaf_g(ifm)%veg_ndvip          &
            ,leaf_g(ifm)%veg_ndvic          &
            ,leaf_g(ifm)%veg_ndvif          &
            ,leaf_g(ifm)%snow_mass, leaf_g(ifm)%snow_depth      &
            ,grid_g(ifm)%glat  &
            ,grid_g(ifm)%glon        ,grid_g(ifm)%topzo         &
            ,grid_g(ifm)%lpw        )

       ! Assignment section for NOFILE leaf-2 variables

       if (icm > 0) then

          !**(JP)** multiple grids not converted yet

          !call fatal_error(h//"**(JP)** multiple grids not converted yet")
          iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
               "**(JP)** multiple grids not converted yet")

          if( nofilflg(ifm) == 0) then

             ! Assign values from coarse grid cells and patches


             if (dumpLocal) then
                call MsgDump(h//" Assign values from coarse grid cells and patches")
             end if
       
             do ipat = 1,npatch
                do j = 1,nnyp(ifm)
                   do i = 1,nnxp(ifm)
                      ic = ipm(i,ifm)
                      jc = jpm(j,ifm)

                      do k = 1,nzg
                         leaf_g(ifm)%soil_water           (k,i,j,ipat) = &
                              leaf_g(icm)%soil_water      (k,ic,jc,ipat)
                         leaf_g(ifm)%soil_energy          (k,i,j,ipat) = &
                              leaf_g(icm)%soil_energy     (k,ic,jc,ipat)
                      enddo

                      do k = 1,nzs
                         leaf_g(ifm)%sfcwater_mass        (k,i,j,ipat) = &
                              leaf_g(icm)%sfcwater_mass   (k,ic,jc,ipat)
                         leaf_g(ifm)%sfcwater_energy      (k,i,j,ipat) = &
                              leaf_g(icm)%sfcwater_energy (k,ic,jc,ipat)
                         leaf_g(ifm)%sfcwater_depth       (k,i,j,ipat) = &
                              leaf_g(icm)%sfcwater_depth  (k,ic,jc,ipat)
                      enddo

                      leaf_g(ifm)%veg_fracarea         (i,j,ipat) = &
                           leaf_g(icm)%veg_fracarea    (ic,jc,ipat)
                      leaf_g(ifm)%veg_lai              (i,j,ipat) = &
                           leaf_g(icm)%veg_lai         (ic,jc,ipat)
                      leaf_g(ifm)%veg_tai              (i,j,ipat) = &
                           leaf_g(icm)%veg_tai         (ic,jc,ipat)
                      leaf_g(ifm)%veg_rough            (i,j,ipat) = &
                           leaf_g(icm)%veg_rough       (ic,jc,ipat)
                      leaf_g(ifm)%veg_height           (i,j,ipat) = &
                           leaf_g(icm)%veg_height      (ic,jc,ipat)
                      leaf_g(ifm)%veg_albedo           (i,j,ipat) = &
                           leaf_g(icm)%veg_albedo      (ic,jc,ipat)
                      leaf_g(ifm)%patch_rough          (i,j,ipat) = &
                           leaf_g(icm)%patch_rough     (ic,jc,ipat)
                      leaf_g(ifm)%patch_wetind         (i,j,ipat) = &
                           leaf_g(icm)%patch_wetind    (ic,jc,ipat)
                      leaf_g(ifm)%soil_rough           (i,j,ipat) = &
                           leaf_g(icm)%soil_rough      (ic,jc,ipat)
                      leaf_g(ifm)%sfcwater_nlev        (i,j,ipat) = &
                           leaf_g(icm)%sfcwater_nlev   (ic,jc,ipat)
                      leaf_g(ifm)%stom_resist          (i,j,ipat) = &
                           leaf_g(icm)%stom_resist     (ic,jc,ipat)
                      leaf_g(ifm)%ground_rsat          (i,j,ipat) = &
                           leaf_g(icm)%ground_rsat     (ic,jc,ipat)
                      leaf_g(ifm)%ground_rvap          (i,j,ipat) = &
                           leaf_g(icm)%ground_rvap     (ic,jc,ipat)
                      leaf_g(ifm)%veg_water            (i,j,ipat) = &
                           leaf_g(icm)%veg_water       (ic,jc,ipat)
                      leaf_g(ifm)%veg_temp             (i,j,ipat) = &
                           leaf_g(icm)%veg_temp        (ic,jc,ipat)
                      leaf_g(ifm)%can_rvap             (i,j,ipat) = &
                           leaf_g(icm)%can_rvap        (ic,jc,ipat)
                      leaf_g(ifm)%can_temp             (i,j,ipat) = &
                           leaf_g(icm)%can_temp        (ic,jc,ipat)
                      leaf_g(ifm)%veg_ndvic            (i,j,ipat) = &
                           leaf_g(icm)%veg_ndvic       (ic,jc,ipat)

                   enddo
                enddo
             enddo

          elseif(nofilflg(ifm) == 1) then

             ! Interpolate from coarse grid.
             ! We can interpolate water patch directly.
             !   For land patches, do this by first averaging all
             !   coarse grid land patches, interpolate, then assign back to all
             !   fine grid land patches.

             if (dumpLocal) then
                call MsgDump(h//" Interpolate from coarse grid")
             end if
             
             call patch_interp(icm,ifm  &
                  ,nzg,nnxp(icm),nnyp(icm),npatch,nzg,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%soil_water,leaf_g(ifm)%soil_water &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )

             call patch_interp(icm,ifm  &
                  ,nzg,nnxp(icm),nnyp(icm),npatch,nzg,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%soil_energy,leaf_g(ifm)%soil_energy &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )

             call patch_interp(icm,ifm  &
                  ,nzs,nnxp(icm),nnyp(icm),npatch,nzs,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%sfcwater_mass,leaf_g(ifm)%sfcwater_mass &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,nzs,nnxp(icm),nnyp(icm),npatch,nzs,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%sfcwater_energy,leaf_g(ifm)%sfcwater_energy &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,nzs,nnxp(icm),nnyp(icm),npatch,nzs,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%sfcwater_depth,leaf_g(ifm)%sfcwater_depth &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )

             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_fracarea,leaf_g(ifm)%veg_fracarea &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_lai,leaf_g(ifm)%veg_lai &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_tai,leaf_g(ifm)%veg_tai &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_rough,leaf_g(ifm)%veg_rough &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_height,leaf_g(ifm)%veg_height &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_albedo,leaf_g(ifm)%veg_albedo &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%patch_rough,leaf_g(ifm)%patch_rough &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_fracarea,leaf_g(ifm)%veg_fracarea &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%patch_wetind,leaf_g(ifm)%patch_wetind &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%soil_rough,leaf_g(ifm)%soil_rough &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%sfcwater_nlev,leaf_g(ifm)%sfcwater_nlev &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%stom_resist,leaf_g(ifm)%stom_resist &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%ground_rsat,leaf_g(ifm)%ground_rsat &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%ground_rvap,leaf_g(ifm)%ground_rvap &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_water,leaf_g(ifm)%veg_water &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_temp,leaf_g(ifm)%veg_temp &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%can_rvap,leaf_g(ifm)%can_rvap &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%can_temp,leaf_g(ifm)%can_temp &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )
             call patch_interp(icm,ifm  &
                  ,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch  &
                  ,leaf_g(icm)%veg_ndvic,leaf_g(ifm)%veg_ndvic &
                  ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area &
                  ,vt3da,vt3db,vt2da,vt2db )

          endif

       endif

       ! Heterogeneous Soil Moisture Initialization

       if (dumpLocal) then
          write(str(1),"(i8)") ifm
          call MsgDump(h//" assign default values for leaf_g soil variables of grid "//trim(adjustl(str(1))))
       end if

       if ((SOIL_MOIST == 'i').or.(SOIL_MOIST == 'I').or.  &
            (SOIL_MOIST == 'a').or.(SOIL_MOIST == 'A')) then

          if (dumpLocal) then
             call MsgDump(h//" invoca soilMoistureInit")
          end if
          
          call soilMoistureInit(nnzp(ifm), nodemxp(mynum,ifm),         &
               nodemyp(mynum,ifm), nzg, nzs, npatch, ifm,              &
               oneBasicFields%theta, oneBasicFields%pi0, oneBasicFields%pp,  &
               leaf_g(ifm)%soil_water, leaf_g(ifm)%soil_energy,        &
               leaf_g(ifm)%soil_text,                                  &
               grid_g(ifm)%glat , grid_g(ifm)%glon, grid_g(ifm)%lpw,   &
               leaf_g(ifm)%seatp, leaf_g(ifm)%seatf, &
               oneControlVars, oneBasicFields, oneTurbFields)

          !-moved to initOneProc
          !        call change_soil_moisture_init(nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm)    &
          !             ,nzg,nzs,npatch,ifm       &
          !             ,oneBasicFields%theta (1,1,1) ,oneBasicFields%pi0 (1,1,1)  &
          !             ,oneBasicFields%pp    (1,1,1)  &
          !             ,leaf_g(ifm)%soil_water   (1,1,1,1)  &
          !             ,leaf_g(ifm)%soil_energy     (1,1,1,1)  &
          !             ,leaf_g(ifm)%soil_text   (1,1,1,1)  &
          !             ,leaf_g(ifm)%sfcwater_mass   (1,1,1,1)  &
          !             ,leaf_g(ifm)%sfcwater_energy (1,1,1,1)  &
          !             ,leaf_g(ifm)%sfcwater_depth  (1,1,1,1)  &
          !             ,grid_g(ifm)%glat   (1,1)       &
          !             ,grid_g(ifm)%glon (1,1)      &
          !             ,grid_g(ifm)%lpw (1,1)    &
          !             ,leaf_g(ifm)%leaf_class(1,1,1)  )

       endif

       ! Override any of the above variable assignments by user-specified changes
       ! to subroutine sfcinit_nofile_user.

       if (dumpLocal) then
          write(str(1),"(i8)") ifm
          call MsgDump(h//" overwrite assigned values of leaf_g variables of grid "//&
               trim(adjustl(str(1)))//" by user selected values")
       end if

       if (dumpLocal) then
          call MsgDump(h//" invoca sfcinit_nofile_user")
       end if
       
       call sfcinit_nofile_user(nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm)      &
            ,nzg,nzs,npatch,ifm              &
            ,oneBasicFields%theta ,oneBasicFields%pi0   &
            ,oneBasicFields%pp    ,oneBasicFields%rv    &
            ,leaf_g(ifm)%soil_water        &
            ,leaf_g(ifm)%soil_energy       &
            ,leaf_g(ifm)%soil_text         &
            ,leaf_g(ifm)%sfcwater_mass     &
            ,leaf_g(ifm)%sfcwater_energy   &
            ,leaf_g(ifm)%sfcwater_depth    &
            ,leaf_g(ifm)%ustar             &
            ,leaf_g(ifm)%tstar             &
            ,leaf_g(ifm)%rstar             &
            ,leaf_g(ifm)%veg_fracarea      &
            ,leaf_g(ifm)%veg_lai           &
            ,leaf_g(ifm)%veg_tai           &
            ,leaf_g(ifm)%veg_rough         &
            ,leaf_g(ifm)%veg_height        &
            ,leaf_g(ifm)%veg_albedo        &
            ,leaf_g(ifm)%patch_area        &
            ,leaf_g(ifm)%patch_rough       &
            ,leaf_g(ifm)%patch_wetind      &
            ,leaf_g(ifm)%leaf_class        &
            ,leaf_g(ifm)%soil_rough        &
            ,leaf_g(ifm)%sfcwater_nlev     &
            ,leaf_g(ifm)%stom_resist       &
            ,leaf_g(ifm)%ground_rsat,leaf_g(ifm)%ground_rvap    &
            ,leaf_g(ifm)%veg_water         &
            ,leaf_g(ifm)%veg_temp          &
            ,leaf_g(ifm)%can_rvap          &
            ,leaf_g(ifm)%can_temp          &
            ,leaf_g(ifm)%veg_ndvip         &
            ,leaf_g(ifm)%veg_ndvic         &
            ,leaf_g(ifm)%veg_ndvif         &
            ,leaf_g(ifm)%snow_mass, leaf_g(ifm)%snow_depth      &
            ,grid_g(ifm)%glat          &
            ,grid_g(ifm)%glon, grid_g(ifm)%topzo, grid_g(ifm)%lpw)

    enddo
  end subroutine GeonestNofile
end module ModNestGeoSst
