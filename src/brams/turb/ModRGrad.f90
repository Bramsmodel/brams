!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module ModRGrad


  use mem_grid, only : jdim   &  !INTENT(IN)
       , dzt    &  !INTENT(IN)
       , hw     &  !INTENT(IN)
       , ht     &  !INTENT(IN)
       , dzm    &  !INTENT(IN)
       , ngrid  &  !INTENT(IN)
       , grid_g    !INTENT(IN)


  implicit none

  private
  public :: grad
  public :: divcart
  
contains


  
  subroutine grad(m1,m2,m3,ia,iz,ja,jz  &
       ,vc3da,vc3db,dir,gpnt)
    integer, intent(IN) :: m1   &
         , m2   &
         , m3   &
         , ia   &
         , iz   &
         , ja   &
         , jz

    real, intent(IN)    :: vc3da(:,:,:)    

    real, intent(INOUT) :: vc3db(:,:,:)

    character(len=*), intent(IN) :: dir,gpnt

    character(len=6) :: optyp

    optyp='GRADNT'

    call rams_grad(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,DIR,GPNT,optyp)

  end subroutine grad

  !------------------------------------------------------------------------------------

  subroutine divcart(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt)
    integer, intent(IN) :: m1   &
         , m2   &
         , m3   &
         , ia   &
         , iz   &
         , ja   &
         , jz

    real, intent(IN)    :: vc3da(:,:,:)    

    real, intent(INOUT) :: vc3db(m1,m2,m3)

    character(len=*), intent(IN) :: dir,gpnt

    character(len=6) :: optyp

    optyp='DIVCRT'

    call rams_grad(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,DIR,GPNT,optyp)

    return
  end subroutine divcart

  !------------------------------------------------------------------------------------

!!$  subroutine divstar(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt)
!!$    integer, intent(IN) :: m1    &
!!$         , m2    &
!!$         , m3    &
!!$         , ia    &
!!$         , iz    &
!!$         , ja    &
!!$         , jz
!!$
!!$    real, intent(IN)    :: vc3da(m1,m2,m3)
!!$
!!$    real, intent(INOUT) :: vc3db(m1,m2,m3)
!!$
!!$    character(len=*), intent(IN) :: dir,gpnt
!!$
!!$    character(len=6) :: optyp
!!$
!!$    optyp='DIVSTR'
!!$
!!$    call rams_grad(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,DIR,GPNT,optyp)
!!$
!!$    return
!!$  end subroutine divstar

  !------------------------------------------------------------------------------------

  subroutine rams_grad(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,optyp)
    integer, intent(IN) :: m1    &
         , m2    &
         , m3    &
         , ia    &
         , iz    &
         , ja    &
         , jz

    real, intent(IN)    :: vc3da(:,:,:)    

    real, intent(INOUT) :: vc3db(m1,m2,m3)

    character(len=*), intent(IN) :: dir,gpnt

    character(len=6), intent(IN) :: optyp

    integer :: jaa,jzz
    real :: vctr2(m1)
    real :: vctr1(m1)

    jaa=ja
    jzz=jz
    if(jdim.eq.0) then
       jaa=1
       jzz=1
    endif

    if(DIR.eq.'XDIR')then
       if(GPNT.eq.'UPNT')then
          call gradxu(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGU  &
               ,GRID_G(NGRID)%RTGT,GRID_G(NGRID)%DXT,DZT  &
               ,GRID_G(NGRID)%FMAPUI,GRID_G(NGRID)%FMAPT  &
               ,GRID_G(NGRID)%F13T  &
               ,HW,VCTR2,'T',JDIM)
       elseif(GPNT.eq.'VPNT')then
          call gradxt(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGV  &
               ,GRID_G(NGRID)%RTGM,GRID_G(NGRID)%DXM,DZT  &
               ,GRID_G(NGRID)%FMAPVI,GRID_G(NGRID)%FMAPM  &
               ,GRID_G(NGRID)%F13M  &
               ,HW,VCTR2,'T',JDIM)
       elseif(GPNT.eq.'WPNT')then
          call gradxt(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGT  &
               ,GRID_G(NGRID)%RTGU,GRID_G(NGRID)%DXU,DZM  &
               ,GRID_G(NGRID)%FMAPTI,GRID_G(NGRID)%FMAPU  &
               ,GRID_G(NGRID)%F13U  &
               ,HT,VCTR2,'W',JDIM)
       elseif(GPNT.eq.'TPNT')then
          call gradxt(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGT  &
               ,GRID_G(NGRID)%RTGU,GRID_G(NGRID)%DXU,DZT  &
               ,GRID_G(NGRID)%FMAPTI,GRID_G(NGRID)%FMAPU  &
               ,GRID_G(NGRID)%F13U  &
               ,HW,VCTR2,'T',JDIM)
       elseif(GPNT.eq.'NPNT')then
          call gradxt(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGV  &
               ,GRID_G(NGRID)%RTGM,GRID_G(NGRID)%DXM,DZM  &
               ,GRID_G(NGRID)%FMAPVI,GRID_G(NGRID)%FMAPM  &
               ,GRID_G(NGRID)%F13M  &
               ,HT,VCTR2,'W',JDIM)
       elseif(GPNT.eq.'OPNT')then
          call gradxu(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGU  &
               ,GRID_G(NGRID)%RTGT,GRID_G(NGRID)%DXT,DZM  &
               ,GRID_G(NGRID)%FMAPUI,GRID_G(NGRID)%FMAPT  &
               ,GRID_G(NGRID)%F13T  &
               ,HT,VCTR2,'W',JDIM)
       elseif(GPNT.eq.'PPNT')then
          call gradxu(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGM  &
               ,GRID_G(NGRID)%RTGV,GRID_G(NGRID)%DXV,DZT  &
               ,GRID_G(NGRID)%FMAPMI,GRID_G(NGRID)%FMAPV  &
               ,GRID_G(NGRID)%F13V  &
               ,HW,VCTR2,'T',JDIM)
       elseif(GPNT.eq.'MPNT')then
          call gradxu(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGM  &
               ,GRID_G(NGRID)%RTGV,GRID_G(NGRID)%DXV,DZM  &
               ,GRID_G(NGRID)%FMAPMI,GRID_G(NGRID)%FMAPV  &
               ,GRID_G(NGRID)%F13V  &
               ,HT,VCTR2,'W',JDIM)
       endif
    elseif(DIR.eq.'YDIR')then
       if(GPNT.eq.'UPNT')then
          call gradyt(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGU  &
               ,GRID_G(NGRID)%RTGM,GRID_G(NGRID)%DYM,DZT  &
               ,GRID_G(NGRID)%FMAPUI,GRID_G(NGRID)%FMAPM  &
               ,GRID_G(NGRID)%F23M  &
               ,HW,VCTR2,'T',JDIM)
       elseif(GPNT.eq.'VPNT')then
          call gradyv(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGV  &
               ,GRID_G(NGRID)%RTGT,GRID_G(NGRID)%DYT,DZT  &
               ,GRID_G(NGRID)%FMAPVI,GRID_G(NGRID)%FMAPT  &
               ,GRID_G(NGRID)%F23T  &
               ,HW,VCTR2,'T',JDIM)
       elseif(GPNT.eq.'WPNT')then
          call gradyt(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGT  &
               ,GRID_G(NGRID)%RTGV,GRID_G(NGRID)%DYV,DZM  &
               ,GRID_G(NGRID)%FMAPTI,GRID_G(NGRID)%FMAPV  &
               ,GRID_G(NGRID)%F23V  &
               ,HT,VCTR2,'W',JDIM)
       elseif(GPNT.eq.'TPNT')then
          call gradyt(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGT  &
               ,GRID_G(NGRID)%RTGV,GRID_G(NGRID)%DYV,DZT  &
               ,GRID_G(NGRID)%FMAPTI,GRID_G(NGRID)%FMAPV  &
               ,GRID_G(NGRID)%F23V  &
               ,HW,VCTR2,'T',JDIM)
       elseif(GPNT.eq.'NPNT')then
          call gradyv(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGV  &
               ,GRID_G(NGRID)%RTGT,GRID_G(NGRID)%DYT,DZM  &
               ,GRID_G(NGRID)%FMAPVI,GRID_G(NGRID)%FMAPT  &
               ,GRID_G(NGRID)%F23T  &
               ,HT,VCTR2,'W',JDIM)
       elseif(GPNT.eq.'OPNT')then
          call gradyt(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGU  &
               ,GRID_G(NGRID)%RTGM,GRID_G(NGRID)%DYM,DZM  &
               ,GRID_G(NGRID)%FMAPUI,GRID_G(NGRID)%FMAPM  &
               ,GRID_G(NGRID)%F23M  &
               ,HT,VCTR2,'W',JDIM)
       elseif(GPNT.eq.'PPNT')then
          call gradyv(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGM &
               ,GRID_G(NGRID)%RTGU,GRID_G(NGRID)%DYU,DZT  &
               ,GRID_G(NGRID)%FMAPMI,GRID_G(NGRID)%FMAPU  &
               ,GRID_G(NGRID)%F23U  &
               ,HW,VCTR2,'T',JDIM)
       elseif(GPNT.eq.'MPNT')then
          call gradyv(m1,m2,m3,ia,iz,jaa,jzz  &
               ,OPTYP,vc3da,vc3db,VCTR1,GRID_G(NGRID)%RTGM  &
               ,GRID_G(NGRID)%RTGU,GRID_G(NGRID)%DYU,DZM  &
               ,GRID_G(NGRID)%FMAPMI,GRID_G(NGRID)%FMAPU  &
               ,GRID_G(NGRID)%F23U  &
               ,HT,VCTR2,'W',JDIM)
       endif
    elseif(DIR.eq.'ZDIR')then
       if(GPNT.eq.'UPNT')then
          call gradzt(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
               ,GRID_G(NGRID)%RTGU,DZM)
       elseif(GPNT.eq.'VPNT')then
          call gradzt(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
               ,GRID_G(NGRID)%RTGV,DZM)
       elseif(GPNT.eq.'WPNT')then
          call gradzw(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
               ,GRID_G(NGRID)%RTGT,DZT)
       elseif(GPNT.eq.'TPNT')then
          call gradzt(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
               ,GRID_G(NGRID)%RTGT,DZM)
       elseif(GPNT.eq.'NPNT')then
          call gradzw(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
               ,GRID_G(NGRID)%RTGV,DZT)
       elseif(GPNT.eq.'OPNT')then
          call gradzw(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
               ,GRID_G(NGRID)%RTGU,DZT)
       elseif(GPNT.eq.'PPNT')then
          call gradzt(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
               ,GRID_G(NGRID)%RTGM,DZM)
       elseif(GPNT.eq.'MPNT')then
          call gradzw(m1,m2,m3,ia,iz,jaa,jzz,vc3da,vc3db  &
               ,GRID_G(NGRID)%RTGM,DZT)
       endif
    endif

    return
  end subroutine rams_grad

  !     ******************************************************************
  !
  !     This is a general subroutine which computes any component of the
  !     gradient or divergence of vc3da and stores it in vc3db.

  subroutine gradxu(m1,m2,m3,ia,iz,ja,jz  &
       ,optyp,vc3da,vc3db,vc1da,rtge,rtgc  &
       ,dx,dz,fmapi,fmap,fq,hq,hq4,lev,jd)
    integer, intent(IN) :: m1  &
         ,m2  &
         ,m3  &
         ,ia  &
         ,iz  &
         ,ja  &
         ,jz  &
         ,jd

    real, intent(IN) :: vc3da(:,:,:)       
    real, intent(IN) ::    &
          rtge(m2,m3)       &
         , rtgc(m2,m3)       &
         , dx(m2,m3)         &
         , fmap(m2,m3)       &
         , fmapi(m2,m3)      &
         , dz(m1)            &
         , fq(m2,m3)         &
         , hq(*)

    real, intent(INOUT)  :: vc3db(m1,m2,m3)   &
         , vc1da(*)          &
         , hq4(*)

    character(len=*), intent(IN) :: optyp,lev

    integer :: i,j,k

    if(OPTYP.eq.'GRADNT')then
       do J=ja,jz
          do I=ia,iz
             do K=1,m1
                vc3db(K,I,J)=(vc3da(K,I,J)*RTGE(I,J)  &
                     -vc3da(K,I-1,J)*RTGE(I-1,J))  &
                     *DX(I,J)/RTGC(I,J)
             enddo
          enddo
       enddo
    else
       do J=ja,jz
          do I=ia,iz
             do K=1,m1
                vc3db(K,I,J)=(vc3da(K,I,J)*RTGE(I,J)  &
                     *FMAPI(I,J)  &
                     -vc3da(K,I-1,J)*RTGE(I-1,J)  &
                     *FMAPI(I-1,J))  &
                     *DX(I,J)/RTGC(I,J)*FMAP(I,J)
             enddo
          enddo
       enddo
    endif

    if(OPTYP.ne.'DIVSTR')then
       if(LEV.eq.'W')then
          do K=1,m1
             HQ4(K)=0.25*HQ(K)
          enddo
       else
          do K=2,m1
             HQ4(K)=0.25*HQ(K-1)
          enddo
       endif

       do j=ja,jz
          do i=ia,iz
             do k=2,m1
                vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
                     +vc3da(k,i-1,j)+vc3da(k-1,i-1,j))
             enddo
             !srf- Antoon - suggested modification
             if(OPTYP .ne. 'GRADNT') VC1DA(2)=0. 
             !srf	 
             do k=2,m1-1
                vc3db(k,i,j)=vc3db(k,i,j)  &
                     +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
             enddo
          enddo
       enddo
       do j=ja,jz
          do i=ia,iz
             vc3db(1,i,j)=vc3db(2,i,j)
          enddo
       enddo
       if(lev.eq.'W')then
          do j=ja,jz
             do i=ia,iz
                vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
             end do
          end do
       else if(lev.eq.'T')then
          do j=ja,jz
             do i=ia,iz
                vc3db(m1,i,j)=vc3db(m1-1,i,j)
             end do
          end do
       end if
       ! **(JP)** fim modificacao
    endif

    return
  end subroutine gradxu

  subroutine gradxt(m1,m2,m3,ia,iz,ja,jz  &
       ,OPTYP,vc3da,vc3db,VC1DA,RTGE,RTGC  &
       ,DX,DZ,FMAPI,FMAP,FQ,HQ,HQ4,LEV,JD)
    integer, intent(IN) :: m1   &
         , m2   &
         , m3   &
         , ia   &
         , iz   & 
         , ja   &
         , jz   &
         , jd

    real, intent(IN)    :: vc3da(:,:,:)      
    real, intent(IN)    ::  &
          RTGE(m2,m3)      &
         , RTGC(m2,m3)      &
         , DX(m2,m3)        &
         , FMAP(m2,m3)      &
         , FMAPI(m2,m3)     &
         , DZ(*)            &
         , FQ(m2,m3)        &
         , HQ(*)

    real, intent(INOUT) :: VC1DA(*)         &
         , vc3db(m1,m2,m3)  &
         , HQ4(*)           

    character(len=*), intent(IN) :: OPTYP,LEV


    integer :: i,j,k

    if(OPTYP.eq.'GRADNT')then
       do J=ja,jz
          do I=ia,iz
             do K=1,m1
                vc3db(K,I,J)=(vc3da(K,I+1,J)*RTGE(I+1,J)  &
                     -vc3da(K,I,J)*RTGE(I,J))  &
                     *DX(I,J)/RTGC(I,J)
             enddo
          enddo
       enddo
    else
       do J=ja,jz
          do I=ia,iz
             do K=1,m1
                vc3db(K,I,J)=(vc3da(K,I+1,J)*RTGE(I+1,J)  &
                     *FMAPI(I+1,J)  &
                     -vc3da(K,I,J)*RTGE(I,J)  &
                     *FMAPI(I,J))  &
                     *DX(I,J)/RTGC(I,J)*FMAP(I,J)
             enddo
          enddo
       enddo
    endif

    if(OPTYP.ne.'DIVSTR')then
       if(LEV.eq.'W')then
          do K=1,m1
             HQ4(K)=0.25*HQ(K)
          enddo
       else
          do K=2,m1
             HQ4(K)=0.25*HQ(K-1)
          enddo
       endif

       do j=ja,jz
          do i=ia,iz
             do k=2,m1
                vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
                     +vc3da(k,i+1,j)+vc3da(k-1,i+1,j))
             enddo
             !srf- Antoon - suggested modification
             if(OPTYP .ne. 'GRADNT') VC1DA(2)=0.
             !srf
             do k=2,m1-1
                vc3db(k,i,j)=vc3db(k,i,j)  &
                     +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
             enddo
          enddo
       enddo
       do j=ja,jz
          do i=ia,iz
             vc3db(1,i,j)=vc3db(2,i,j)
          enddo
       enddo
       if (lev.eq.'W') then
          do j=ja,jz
             do i=ia,iz
                vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
             enddo
          enddo
       else if(lev.eq.'T') then
          do j=ja,jz
             do i=ia,iz
                vc3db(m1,i,j)=vc3db(m1-1,i,j)
             enddo
          enddo
       end if
       ! **(JP)** fim modificacao
    endif

    return
  end subroutine gradxt

  !

  subroutine gradyv(m1,m2,m3,ia,iz,ja,jz  &
       ,OPTYP,vc3da,vc3db,VC1DA,RTGE,RTGC  &
       ,DY,DZ,FMAPI,FMAP,FQ,HQ,HQ4,LEV,JD)
    integer, intent(IN) :: m1   &
         , m2   &
         , m3   &
         , ia   &
         , iz   &
         , ja   &
         , jz   &
         , jd

    real, intent(IN)  :: vc3da(:,:,:)       
    real, intent(IN)  ::    &
          RTGE(m2,m3)       &
         , RTGC(m2,m3)       &
         , DY(m2,m3)         &
         , FMAP(m2,m3)       &
         , FMAPI(m2,m3)      &
         , DZ(*)             &
         , FQ(m2,m3)         &
         , HQ(*)

    real, intent(INOUT) :: vc3db(m1,m2,m3)   &
         , HQ4(*)            &
         , VC1DA(*)


    character(len=*) :: OPTYP,LEV

    integer :: i,j,k

    if(OPTYP.eq.'GRADNT')then
       do J=ja,jz
          do I=ia,iz
             do K=1,m1
                vc3db(K,I,J)=(vc3da(K,I,J)*RTGE(I,J)  &
                     -vc3da(K,I,J-jd)*RTGE(I,J-jd))  &
                     *DY(I,J)/RTGC(I,J)
             enddo
          enddo
       enddo
    else
       do J=ja,jz
          do I=ia,iz
             do K=1,m1
                vc3db(K,I,J)=(vc3da(K,I,J)*RTGE(I,J)  &
                     *FMAPI(I,J)  &
                     -vc3da(K,I,J-jd)*RTGE(I,J-jd)  &
                     *FMAPI(I,J-jd))  &
                     *DY(I,J)/RTGC(I,J)*FMAP(I,J)
             enddo
          enddo
       enddo
    endif

    if(OPTYP.ne.'DIVSTR')then
       if(LEV.eq.'W')then
          do K=1,m1
             HQ4(K)=0.25*HQ(K)
          enddo
       else
          do K=2,m1
             HQ4(K)=0.25*HQ(K-1)
          enddo
       endif

       do j=ja,jz
          do i=ia,iz
             do k=2,m1
                vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
                     +vc3da(k,i,j-jd)+vc3da(k-1,i,j-jd))
             enddo
             !srf- antoon mod suggested
             if(OPTYP .ne. 'GRADNT') VC1DA(2)=0.
             !srf
             do k=2,m1-1
                vc3db(k,i,j)=vc3db(k,i,j)  &
                     +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
             enddo
          enddo
       enddo
       do j=ja,jz
          do i=ia,iz
             vc3db(1,i,j)=vc3db(2,i,j)
          enddo
       enddo
       if (lev.eq.'W') then
          do j=ja,jz
             do i=ia,iz
                vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
             enddo
          enddo
       else if (lev.eq.'T') then
          do j=ja,jz
             do i=ia,iz
                vc3db(m1,i,j)=vc3db(m1-1,i,j)
             enddo
          enddo
       end if
       ! **(JP)** fim modificacao
    endif

    return
  end subroutine gradyv

  subroutine gradyt(m1,m2,m3,ia,iz,ja,jz  &
       ,OPTYP,vc3da,vc3db,VC1DA,RTGE,RTGC  &
       ,DY,DZ,FMAPI,FMAP,FQ,HQ,HQ4,LEV,JD)
    integer, intent(IN) :: m1    &
         , m2    &
         , m3    &
         , ia    &
         , iz    &
         , ja    &
         , jz    &
         , jd

    real, intent(IN) :: vc3da(:,:,:)       
    real, intent(IN) ::    &
          RTGE(m2,m3)       &
         , RTGC(m2,m3)       &
         , DY(m2,m3)         &
         , FMAP(m2,m3)       &
         , FMAPI(m2,m3)      &
         , DZ(*)             &
         , FQ(m2,m3)         &
         , HQ(*)

    real, intent(INOUT) :: vc3db(m1,m2,m3)   &
         , HQ4(*)            &
         , VC1DA(*)          

    character(len=*), intent(IN) :: OPTYP, LEV

    integer :: i,j,k

    if(OPTYP.eq.'GRADNT')then
       do J=ja,jz
          do I=ia,iz
             do K=1,m1
                vc3db(K,I,J)=(vc3da(K,I,J+jd)*RTGE(I,J+jd)  &
                     -vc3da(K,I,J)*RTGE(I,J))  &
                     *DY(I,J)/RTGC(I,J)
             enddo
          enddo
       enddo
    else
       do J=ja,jz
          do I=ia,iz
             do K=1,m1
                vc3db(K,I,J)=(vc3da(K,I,J+jd)*RTGE(I,J+jd)  &
                     *FMAPI(I,J+jd)  &
                     -vc3da(K,I,J)*RTGE(I,J)  &
                     *FMAPI(I,J))  &
                     *DY(I,J)/RTGC(I,J)*FMAP(I,J)
             enddo
          enddo
       enddo
    endif

    if(OPTYP.ne.'DIVSTR')then
       if(LEV.eq.'W')then
          do K=1,m1
             HQ4(K)=0.25*HQ(K)
          enddo
       else
          do K=2,m1
             HQ4(K)=0.25*HQ(K-1)
          enddo
       endif

       do j=ja,jz
          do i=ia,iz
             do k=2,m1
                vc1da(k)=hq4(k)*(vc3da(k,i,j)+vc3da(k-1,i,j)  &
                     +vc3da(k,i,j+jd)+vc3da(k-1,i,j+jd))
             enddo
             !srf- antoon sug. modification
             if(OPTYP .ne. 'GRADNT') VC1DA(2)=0.
             !srf
             do k=2,m1-1
                vc3db(k,i,j)=vc3db(k,i,j)  &
                     +fq(i,j)*dz(k)*(vc1da(k+1)-vc1da(k))
             enddo
          enddo
       enddo
       do j=ja,jz
          do i=ia,iz
             vc3db(1,i,j)=vc3db(2,i,j)
          enddo
       enddo
       if (lev.eq.'W') then
          do j=ja,jz
             do i=ia,iz
                vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
             enddo
          enddo
       else if (lev.eq.'T') then
          do j=ja,jz
             do i=ia,iz
                vc3db(m1,i,j)=vc3db(m1-1,i,j)
             enddo
          enddo
       end if
       ! **(JP)** fim modificacao
    endif

    return
  end subroutine gradyt

  subroutine gradzw(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,RTGC,DZ)
    integer, intent(IN) :: m1  &
         , m2  &
         , m3  &
         , ia  &
         , iz  &
         , ja  &
         , jz

    real, intent(IN)    :: vc3da(:,:,:)      
    real, intent(IN)    ::   &
          RTGC(m2,m3)      &
         , DZ(*)

    real, intent(INOUT) :: vc3db(m1,m2,m3) 

    integer :: i,j,k

    do J=ja,jz
       do I=ia,iz
          do K=2,m1
             vc3db(K,I,J)=(vc3da(K,I,J)-vc3da(K-1,I,J))*DZ(K)  &
                  /RTGC(I,J)
          enddo
       enddo
    enddo
    return
  end subroutine gradzw

  subroutine gradzt(m1,m2,m3,ia,iz,ja,jz,vc3da,vc3db,RTGC,DZ)
    integer :: m1,m2,m3,ia,iz,ja,jz

    real, intent(IN)    :: vc3da(:,:,:)     
    real, intent(IN)    :: &
          RTGC(m2,m3)     &
         , DZ(*)

    real, intent(INOUT) :: vc3db(m1,m2,m3)

    integer :: i,j,k

    do J=ja,jz
       do I=ia,iz
          do K=1,m1-1
             vc3db(K,I,J)=(vc3da(K+1,I,J)-vc3da(K,I,J))*DZ(K)  &
                  /RTGC(I,J)
          enddo
       enddo
    enddo
    return
  end subroutine gradzt
end module ModRGrad
