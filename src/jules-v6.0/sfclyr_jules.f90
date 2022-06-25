!############################# Change Log ##################################
! Interface entre BRAMS e JULES
! 
module ModSfcLyrJules

  use ModLeaf3Init, only: &
       sfcdata
  
  use ModBasicFields, only: &
       BasicFields

  use io_params, only : &
       frqanl, &
       afilout, &
       frqhis, &
       hfilin

  use node_mod, only: &
       MYNUM

  use leaf_coms, only: &
       slcpd, &
       slmsts, &
       gzotheta 

  use mem_leaf, only: &
       slz, &
       leaf_g

  use mem_jules, only: &
       jules_g

  use mem_grid, only: &
       npatch, &
       nzg, &
       grid_g, &
       time, &
       dtlong, &
       iyear1, &
       imonth1, &
       idate1 , &
       istp, &
       itime1, &
       timmax, &
       zt, &
       ngrid, &
       runtype                   

  use rconstants, only: &
       cpi, &
       cp, &
       alvl, &
       vonk, &
       g, &
       p00, &
       cpor

  use mem_radiate, only: &
       radiate_g 

  use ModTurbFields, only: &
       TurbFields

  use mem_cuparm, only: &
       cuparm_g, &
       nnqparm

  use micphys, only: &
       level

  use mem_micro, only: &
       micro_g

  use chem1_list, only: &
       CO2, &
       chemical_mechanism

  use mem_chem1, only: &
       chem1_g, &
       chemistry

  !--- Modulos do JULES ---
  use io_constants, only: &
       max_file_name_len

  use mem_brams_jules, only: &
       mynumB, &
       glatB, &
       glonB, &
       nxB, &
       nyB, &
       land_fracB, &
       precipB, &
       swdownB, &
       lwdownB, &
       diff_radB, &
       tempB, &
       upsB, &
       vpsB, &
       pstarB, &
       qB, &
       fracB, &
       timestepB, &
       main_run_startB, &
       main_run_endB, &
       ntimestepB, &
       sm_levelsB, &
       dzsoilB, &
       sthuB, &
       tsoilB, &
       tstarB, &
       output_periodB, &
       dir_run_idB, &
       dump_periodB, &
       runtypeB, &
       hfilinB

  use gridbox_mean_mod, only: &
       surftiles_to_gbm

  use jules_surface_types_mod, only: &
       npft, &
       ntype

  use csigma, only: &
       sbcon

  use model_time_mod, only: &
       end_of_run

  use ancil_info, only: &
       land_pts, &
       row_length

  use jules_fields_mod, only: &
       progs, &
       psparms, &
       ainfo, &
       trifctltype

  use sf_diags_mod, only: &
       sf_diag

  use datetime_mod, only: &
       datetime, &
       datetime_to_string, &
       datetime_advance, &
       datetime_diff

  use gridmean_fluxes, only: &
       fqw_1_ij,         &  ! latent heat flux 
       ftl_1_ij,         &  ! fluxo de calor sensivel do JULES
       taux_1_ij,        &  !   W'ly component of surface wind stress (N/sq m)
       tauy_1_ij            !   S'ly component of surface wind stress (N/sq m)

  use fluxes, only : &
       emis_surft, &
       land_albedo_ij

  implicit none

  private

  public :: sfclyr_jules

contains


  subroutine sfclyr_jules(mzp,mxp,myp,iaI,izI,jaI,jzI,jdim,julesFile,&
       oneBasicFields, oneTurbFields)

    !--- Modulos do BRAMS ---

    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: iaI
    integer, intent(in) :: izI
    integer, intent(in) :: jaI
    integer, intent(in) :: jzI
    integer, intent(in) :: jdim
    character(len=*), intent(in) :: julesFile
    type(BasicFields), pointer, intent(in) :: oneBasicFields
    type(TurbFields), pointer, intent(in) :: oneTurbFields

    integer               :: nsoil, fase,ia,iz,ja,jz,hh,mm
    integer, parameter :: fat_dtlong=1  ! > 1 para nao executar o JULES em todos os timestep do BRAMS


    integer            :: ng,i ,j, k2,l,n,ip,k,ntimestep, itime1_seg

    logical              :: there, start=.true.

    real               :: hcpi, zoverl,cx,wtol,psin,piv &
         ,tempk,fracliq,zts2,ths2,ustar2,tstar2,rstar2

    real, dimension(:), allocatable :: rlongupJ

    character (LEN=2)  :: tB_str    

    real, dimension(mxp,myp) :: rvs2,ups2,vps2,pp2,temp2,pcpgl

    logical, save , allocatable :: run_jules(:)

    real, allocatable :: dens2(:,:)

    type(datetime) :: dti,dtf  ! Placeholder for a datetime in a calculation

    data wtol/1e-20/

    character(LEN=max_file_name_len) :: nml_dir  ! Directory containing namelists

    character(LEN=256) :: aux

    nml_dir=trim(julesFile)

    inquire(FILE=trim(nml_dir)//'/drive.nml',EXIST=there)

    if (.not. there) then
       print*;print*;print*;print*, 'Not found:  '//trim(nml_dir)//'/drive.nml'
       print*;print*, 'Check JULESIN variable in RAMSIN'
       print*;stop
    endif

    ia=iaI-1
    iz=izI+1
    ja=jaI-1
    jz=jzI+1

    mynumB=MYNUM
    ng=ngrid
    nxB=iz-ia+1
    nyB=jz-ja+1

    if (.not. allocated(swdownB)) then
       allocate( swdownB(nxB,nyB),lwdownB(nxB,nyB),diff_radB(nxB,nyB),precipB(nxB,nyB), &
            tempB(nxB,nyB),upsB(nxB,nyB),vpsB(nxB,nyB),pstarB(nxB,nyB),qB(nxB,nyB), &
            tstarB(nxB,nyB) )
    endif

    !--- Precipitacao total ---
    pcpgl(:,:)=0.
    if (nnqparm(ng) > 0 .and. level >= 3) then
       pcpgl(:,:)=cuparm_g(ng)%conprr(:,:) + micro_g(ng)%pcpg(:,:)
    elseif(nnqparm(ng) == 0 .and. level >= 3) then
       pcpgl(:,:)=micro_g(ng)%pcpg(:,:)
    endif

    rvs2(:,:) = oneBasicFields%rv(2,:,:)
    ups2(:,:) = oneBasicFields%up(2,:,:)
    vps2(:,:) = oneBasicFields%vp(2,:,:)

    pp2(:,:)=p00 * ((oneBasicFields%pi0(2,:,:) + oneBasicFields%pp(2,:,:)) * cpi) ** cpor

    hcpi = .5 * cpi
    do j=ja,jz
       do i=ia,iz
          k2=int(grid_g(ng)%lpw(i,j))
          piv = hcpi * (oneBasicFields%pi0(k2-1,i,j) + oneBasicFields%pi0(k2,i,j)   &
               + oneBasicFields%pp(k2-1,i,j) + oneBasicFields%pp(k2,i,j))
          temp2(i,j) = oneBasicFields%theta(k2,i,j) * piv   != airtemp = tstar_tile
       enddo
    enddo

    runtypeB=runtype
    hfilinB=hfilin

    swdownB(:,:)   =radiate_g(ng)%rshort(ia:iz, ja:jz) 
    lwdownB(:,:)   =radiate_g(ng)%rlong(ia:iz, ja:jz)

    !--- TMP ate resolver o problema de borda da radiacao ---{
    swdownB(ia,:)=swdownB(ia+1,:)
    swdownB(iz,:)=swdownB(iz-1,:)
    swdownB(:,ja)=swdownB(:,ja+1)
    swdownB(:,jz)=swdownB(:,jz-1)

    lwdownB(ia,:)=lwdownB(ia+1,:)
    lwdownB(iz,:)=lwdownB(iz-1,:)
    lwdownB(:,ja)=lwdownB(:,ja+1)
    lwdownB(:,jz)=lwdownB(:,jz-1)
    !------------------------------------------------}

    precipB(:,:)   =pcpgl(ia:iz, ja:jz)
    diff_radB(:,:) =0.0 !radiate_g(ng)%rshortdif(ia:iz, ja:jz)  !TMP
    tempB(:,:)     =temp2(ia:iz, ja:jz)
    upsB(:,:)      =ups2(ia:iz, ja:jz)
    vpsB(:,:)      =vps2(ia:iz, ja:jz)
    pstarB(:,:)    =pp2(ia:iz, ja:jz)
    qB(:,:)        =rvs2(ia:iz, ja:jz)/(1+rvs2(ia:iz, ja:jz))

    !print*, ia,iz,ja,jz,swdownB(10,1),swdownB(10,18),swdownB(10,2), lwdownB(10,1),lwdownB(10,18),lwdownB(10,2)
    !stop 'em /home/demerval/JULES-BRAMS/JULES5.8-BRAMS5.4/source/BRAMS5.4/src/jules/sfclyr_jules.f90'

    if (start) then
       start=.false.

       dump_periodB=nint(frqhis)
       output_periodB=nint(frqanl)
       dir_run_idB=trim(afilout)

       sm_levelsB=nzg
       if (.not. allocated(dzsoilB)) allocate(dzsoilB(sm_levelsB))
       if (.not. allocated(sthuB)) allocate(sthuB(sm_levelsB,nxB,nyB))
       if (.not. allocated(tsoilB)) allocate(tsoilB(sm_levelsB,nxB,nyB))

       dzsoilB(1)=-1.0*slz(sm_levelsB)
       do k=2,sm_levelsB
          dzsoilB(k)=-1.0*slz(sm_levelsB-k+1) + slz(sm_levelsB-k+2)
       enddo

       call sfcdata

       !--- ACOPLANDO - BRAMS p/ JULES ------{
       do j=ja,jz
          do i=ia,iz
             do k=1,sm_levelsB
                nsoil = nint(leaf_g(ng)%soil_text(k,i,j,2))
                sthuB(sm_levelsB-k+1,i,j)=leaf_g(ng)%soil_water(k,i,j,2)/slmsts(nsoil)

                call qwtk(leaf_g(ng)%soil_energy(k,i,j,2),     &
                     leaf_g(ng)%soil_water(k,i,j,2)*1.e3, slcpd(nsoil),tempk,fracliq)
                tsoilB(sm_levelsB-k+1,i,j)=tempk
             enddo
          enddo
       enddo

       tstarB(:,:)=temp2(ia:iz, ja:jz)
       !-------------------------------------} 

       allocate(glatB(nxB,nyB),glonB(nxB,nyB),land_fracB(nxB,nyB))
       glatB=grid_g(ng)%glat(ia:iz,ja:jz)
       glonB=grid_g(ng)%glon(ia:iz,ja:jz)

       !WHERE(leaf_g(ng)%patch_area(ia:iz,ja:jz,1)>0.99)
       !   land_fracB=0.0
       !ELSEWHERE
       land_fracB=1.0  ! calculando em todos os pontos, pois senao dah problema no history, pois 
       ! o patch_area estah mudando ao longo do tempo.
       !END WHERE

       !write(44,*)land_fracB
       !print*,ia,iz,ja,jz,nxB,nyB,time
       !stop 'em /home/demerval/JULES-BRAMS/JULES5.8-BRAMS5.4/source/BRAMS5.4/src/jules/sfclyr_jules.f90'



       timestepB=nint(dtlong*fat_dtlong)  !DSM foi colocado em sfclyr_jules.f90 um fator multiplicando o dtlong para nao chamar o JULES a cada timeStep do BRAMS
       dti%year=iyear1
       dti%month=imonth1
       dti%day=idate1
       itime1_seg=(itime1/100)*3600 + mod(itime1,100)*60
       dti%time=itime1_seg

       dtf = datetime_advance(dti, nint(timmax))

       !print*, dti%time,timmax

       if (trim(runtypeB)=='HISTORY') then
          aux=trim(hfilinB(index(hfilinB,'-head.txt',BACK = .true.)-17:index(hfilinB,'-head.txt',BACK = .true.)-1))
          read(aux(1:4),*) dti%year
          read(aux(6:7),*) dti%month
          read(aux(9:10),*) dti%day

          read(aux(12:13),*) hh
          read(aux(14:15),*) mm
          dti%time=hh*3600+mm*60
       endif

       main_run_startB=datetime_to_string(dti)

       main_run_endB=datetime_to_string(dtf)

       ntimestep=datetime_diff(dtf, dti)/dtlong+1
       allocate (run_jules(ntimestep))

       !print*,'b===>',dti
       !print*,dtf
       !print*,fat_dtlong,ntimestep,fat_dtlong
       !stop 'aasfhd'

       ntimestepB=1
       run_jules=.false.
       run_jules(1)=.true.
       do i=1+fat_dtlong,ntimestep,fat_dtlong
          run_jules(i)=.true.
          ntimestepB=ntimestepB+1
       enddo

       if(mynum==1) then
          open(unit=66,file='jules.log',status='old',position='append',action='write')
          write(unit=66,fmt='(A)') '---- Inicio da FASE-1 ---'
          close(unit=66)
       endif
       fase=1; call jules_subroutine(nml_dir,fase)  ! faz leitura do namelist (le o ntype)

       !--- Converte o vegetacao do BRAMS (Leaf3) para a do JULES ---{
       allocate(fracB(nxB,nyB,ntype))
       call frac_from_leaf(ntype,nxB,nyB,npatch   &
            ,leaf_g(ng)%patch_area(ia:iz,ja:jz,:)  &
            ,leaf_g(ng)%leaf_class(ia:iz,ja:jz,:)  &
            ,fracB)
       !-------------------------------------------------------------}

       if(mynum==1) then
          open(unit=66,file='jules.log',status='old',position='append',action='write')
          write(unit=66,fmt='(A)') '---- Inicio da FASE-2 ---'
          close(unit=66)
       endif
       fase=2; call jules_subroutine(nml_dir,fase)  ! finaliza a inicializacao do JULES

       if (allocated(dzsoilB)) deallocate(dzsoilB)
       if (allocated(sthuB)) deallocate(sthuB)
       if (allocated(tsoilB)) deallocate(tsoilB)

    endif

    if (.not. run_jules(istp)) return

    if(mynum==1 .and. time==0) then
       open(unit=66,file='jules.log',status='old',position='append',action='write')
       write(unit=66,fmt='(A)') '---- Inicio da FASE-3 ---'
       close(unit=66)
    endif
    !--- ACOPLANDO - BRAMS p/ JULES ------{
    !DSM    IF ( (chemical_mechanism .eq. 'RELACS_TUV' .or. chemical_mechanism .eq. 'CO2') &
    !DSM          .and. CHEMISTRY == 0) THEN
    !DSM       aerotype%co2_3d_ij(1:nxB, 1:nyB) = chem1_g(CO2,ng)%sc_p (2, ia:iz, ja:jz) * 1.E-9
    !DSM    ELSE
    !DSM       aerotype%co2_3d_ij = 384.
    !DSM    END IF
    !-------------------------------------}

    fase=3; call jules_subroutine(nml_dir,fase) !--- para cada timstep do BRAMS

    !--- ACOPLANDO - JULES p/ BRAMS ------{
    !--- Acoplando o albedt ---{
    radiate_g(ng)%albedt(:,:)=(land_albedo_ij(:,:,1)+land_albedo_ij(:,:,2)+land_albedo_ij(:,:,3)+land_albedo_ij(:,:,4))/4

    allocate(rlongupJ(land_pts))
    rlongupJ(:)=sbcon * surftiles_to_gbm(emis_surft * progs%tstar_surft**4, ainfo)
    do l=1,land_pts
       j = ( ainfo%land_index(l)-1 ) / row_length + 1
       i = ainfo%land_index(l) - ( j-1 ) * row_length
       if (i<1 .or. i>nxB .or. j<1 .or. j>nyB .or. i+ia-1>iz .or. j+ja-1>jz) then
          print*, "ERRO... conversao incorreta de l para i,j -> l,i,j=",l,i,j
          stop
       endif
       !--- Acoplando o rlongup ---{
       radiate_g(ng)%rlongup(i+ia-1,j+ja-1) = rlongupJ(l)

       !--- Acoplando a umidade do solo ---
       do k=1,sm_levelsB
          nsoil = nint(leaf_g(ng)%soil_text(k,i,j,2))
          leaf_g(ng)%soil_water(k,i+ia-1,j+ja-1,2:npatch)=psparms%sthu_soilt(l,1,sm_levelsB-k+1)*slmsts(nsoil) + &
               psparms%sthf_soilt(l,1,sm_levelsB-k+1)*slmsts(nsoil)
       enddo

       !--- Acoplando o lai ---
       leaf_g(ng)%veg_lai(i+ia-1,j+ja-1,2:npatch) = sum(progs%lai_pft(l,:) * ainfo%frac_surft(l,1:npft))

    enddo
    deallocate(rlongupJ)

    !--- Acoplando fluxo de calor latente ---
    oneTurbFields%sflux_r(ia:iz,ja:jz) = sf_diag%latent_heat(:,:)/alvl

    !--- Acoplando fluxo de calor sensivel ---
    oneTurbFields%sflux_t(ia:iz,ja:jz) = ftl_1_ij(:,:)/cp

    allocate(dens2(nxB,nyB))
    !--- Acoplando fluxo sflux_u e sflux_v ---
    dens2(:,:) = (oneBasicFields%dn0(1,ia:iz,ja:jz) + oneBasicFields%dn0(2,ia:iz,ja:jz)) * .5
    oneTurbFields%sflux_u(ia:iz,ja:jz) = -1*taux_1_ij(:,:)/dens2(:,:)
    oneTurbFields%sflux_v(ia:iz,ja:jz) = -1*tauy_1_ij(:,:)/dens2(:,:)

    do j=ja,jz
       do i=ia,iz
          !--- Acoplando ustar, tstar, rstar ---
          ustar2 = max(0.000001,sqrt( &
               sqrt((oneTurbFields%sflux_u(i+ia-1,j+ja-1))**2 + &
               (oneTurbFields%sflux_v(i+ia-1,j+ja-1))**2)))
          tstar2 = ftl_1_ij(i,j)/(dens2(i,j) * cp * ustar2)
          rstar2 = sf_diag%latent_heat(i,j)/(dens2(i,j) * alvl * ustar2)
          leaf_g(ng)%ustar(i+ia-1,j+ja-1,1:npatch)=ustar2
          leaf_g(ng)%tstar(i+ia-1,j+ja-1,1:npatch)=tstar2
          leaf_g(ng)%rstar(i+ia-1,j+ja-1,1:npatch)=rstar2

          !--- Acoplando fluxo sflux_w ---
          zts2 = zt(2) * grid_g(ng)%rtgt(i+ia-1,j+ja-1)
          ths2 = oneBasicFields%theta(2,i+ia-1,j+ja-1)
          gzotheta = g * zts2 / ths2
          zoverl = gzotheta * vonk * tstar2 / (ustar2 * ustar2)
          if (zoverl < 0.) then
             cx = zoverl * sqrt(sqrt(1. - 15. * zoverl))
          else
             cx = zoverl / (1.0 + 4.7 * zoverl)
          endif
          psin = sqrt((1.-2.86 * cx) / (1. + cx * (-5.39 + cx * 6.998 )))
          ip=2
          oneTurbFields%sflux_w(i+ia-1,j+ja-1) = (0.27 * max(6.25 * (1. - cx) * &
               psin,wtol)- 1.18 * cx * psin) * ustar2 * ustar2 * leaf_g(ng)%patch_area(i+ia-1,j+ja-1,ip)    
       enddo
    enddo
    deallocate(dens2)
    !-------------------------------------}

    !--- ESCREVENDO AS VARIAVEIS DO JULES NO OUTPUT DO BRAMS ---{
    jules_g(ng)%u10mj(ia:iz,ja:jz)=sf_diag%u10m(:,:)
    jules_g(ng)%v10mj(ia:iz,ja:jz)=sf_diag%v10m(:,:)
    jules_g(ng)%t2mj(ia:iz,ja:jz)=sf_diag%t1p5m(:,:)
    jules_g(ng)%rv2mj(ia:iz,ja:jz)=sf_diag%q1p5m(:,:)/(1-sf_diag%q1p5m(:,:))


    do l=1,land_pts
       j = ( ainfo%land_index(l)-1 ) / row_length + 1
       i = ainfo%land_index(l) - ( j-1 ) * row_length


       if (i<1 .or. i>nxB .or. j<1 .or. j>nyB .or. i+ia-1>iz .or. j+ja-1>jz) then
          print*, "ERRO... conversao incorreta de l para i,j -> l,i,j=",l,i,j
          stop
       endif

       jules_g(ng)%gpp(i+ia-1,j+ja-1)=trifctltype%gpp_gb(l)
       jules_g(ng)%resp_p(i+ia-1,j+ja-1)=trifctltype%resp_p_gb(l)
       jules_g(ng)%npp(i+ia-1,j+ja-1)=trifctltype%npp_gb(l)
       jules_g(ng)%resp_s(i+ia-1,j+ja-1)=sum(sum(trifctltype%resp_s_soilt(l,1,:,:), 2), 1)
    enddo
    !-----------------------------------------------------------}

    !IF ( end_of_run ) THEN  !esta funcao eh do JULES - pode tambem utilizar o ultimo timestep do BRAMS
    !LFR
    if(time>=timmax) then
       if(mynum==1) then 
          open(unit=66,file='jules.log',status='old',position='append',action='write')
          write(unit=66,fmt='(A)') '---- Inicio da FASE-4 ---'
          close(unit=66)
       endif


       fase=4; call jules_subroutine(nml_dir,fase) ! apos o ultimo timestep do BRAMS
    endif

  end subroutine sfclyr_jules
end module ModSfcLyrJules

!{DSM
!--------------------------------------------------------------------------------------------------
!--- Encontra a fracao de vegetacao a partir do mapa da fracao de vegetacao definida pelo leaf3 ---
!--------------------------------------------------------------------------------------------------
subroutine frac_from_leaf( ntype,nx,ny,npatch,patch_area,leaf_class,frac )
  implicit none
  integer, intent(IN)  :: ntype,nx,ny,npatch
  real, intent(IN)     :: patch_area(nx,ny,npatch), leaf_class(nx,ny,npatch)
  real, intent(OUT)    :: frac(nx,ny,ntype)

  integer              :: j,i,tJ,n
  character (LEN=80)   :: veg(ntype)
  character (LEN=2)    :: tB_str          

  !--- Convertendo o tipo de vegetacao do BRAMS para o JULES ---
  call brams2jules(veg,ntype)

  if (ntype == 9 ) then
     frac=0.

     do j=1,ny
        do i=1,nx
           if (i > nx .or. j > ny) then
              print*, 'ERRO!!! i > nx ou j > ny - i, nx, j, ny =',i, nx, j, ny
              stop
           endif

           do n=2,npatch
              write(tB_str,'(i2.2)') nint(leaf_class(i,j,n))

              !--- Encontrando o tipo correspondente ao JULES ---
              do tJ=1,ntype+1    !--- o +1 eh apenas para checar se foi encontrado (condicao abaixo)
                 if (index(veg(tJ),tB_str)/=0) exit
              enddo

              !--- Checando se encontrou um indice valido para a vegetacao do JULES ---
              if (tJ>ntype) then
                 print*, 'ERRO!!! Nao foi encontrado uma correspondencia entre BRAMS e JULES'
                 stop
              endif

              frac(i,j,tJ)=frac(i,j,tJ) + max(0.,patch_area(i,j,n))

           enddo  !--- DO n=1,npatch

           frac(i,j,7)=frac(i,j,7)+max(0.,patch_area(i,j,1))

           where ( frac(i,j,:) <= 1.00004E-06 ) frac(i,j,:) = 1.00004E-06  ! valor minimo que o JULES trabalha

           frac(:,:,9)=0.0 !--- Mas ice fraction deve iniciar com zero.

           n=1
           do while (sum(frac(i,j,:)) > 1.0)
              if (frac(i,j,n)>0.01) frac(i,j,n)=frac(i,j,n)-0.01                             !delta=0.01
              n=n+1
              if (n>8) n=1     
           enddo

           !--- Para garantir que a fracao total seja igual a 1 ---
           frac(i,j,1)=1.-( frac(i,j,2)+frac(i,j,3)+frac(i,j,4)     &
                +frac(i,j,5)+frac(i,j,6)+frac(i,j,7)+frac(i,j,8)+frac(i,j,9) )

        enddo  !i=1,nx
     enddo  !j=1,ny
  elseif (ntype == 13 ) then
     frac=0.

     do j=1,ny
        do i=1,nx
           if (i > nx .or. j > ny) then
              print*, 'ERRO!!! i > nx ou j > ny - i, nx, j, ny =',i, nx, j, ny
              stop
           endif

           do n=2,npatch
              write(tB_str,'(i2.2)') nint(leaf_class(i,j,n))

              !--- Encontrando o tipo correspondente ao JULES ---
              do tJ=1,ntype+1    !--- o +1 eh apenas para checar se foi encontrado (condicao abaixo)
                 if (index(veg(tJ),tB_str)/=0) exit
              enddo

              !--- Checando se encontrou um indice valido para a vegetacao do JULES ---
              if (tJ>ntype) then
                 print*, 'ERRO!!! Nao foi encontrado uma correspondencia entre BRAMS e JULES'
                 stop
              endif

              frac(i,j,tJ)=frac(i,j,tJ) + max(0.,patch_area(i,j,n))

           enddo  !--- DO n=1,npatch

           frac(i,j,11)=frac(i,j,11)+max(0.,patch_area(i,j,1))

           where ( frac(i,j,:) <= 1.00004E-06 ) frac(i,j,:) = 1.00004E-06  ! valor minimo que o JULES trabalha

           frac(:,:,13)=0.0 !--- Mas ice fraction deve iniciar com zero.

           n=1
           do while (sum(frac(i,j,:)) > 1.0)
              if (frac(i,j,n)>0.01) frac(i,j,n)=frac(i,j,n)-0.01                             !delta=0.01
              n=n+1
              if (n>12) n=1     
           enddo

           !--- Para garantir que a fracao total seja igual a 1 ---
           frac(i,j,1)=1.-( frac(i,j,2)+frac(i,j,3)+frac(i,j,4)     &
                +frac(i,j,5)+frac(i,j,6)+frac(i,j,7)+frac(i,j,8)  &
                +frac(i,j,9)+frac(i,j,10)+frac(i,j,11)+frac(i,j,12)+frac(i,j,13) )

        enddo  !i=1,nx
     enddo  !j=1,ny
  else
     print*, "ntype=",ntype
     print*, "ATENCAO... ntype <> 9 ou de 13, Deve-se ajustar a subrotina brams2jules em init_frac2brams.f90"
     stop  
  endif

end subroutine frac_from_leaf

!------------------------------------------------------------------------------
!--- Convertendo o tipo de vegetacao do BRAMS para o JULES ---
!-------------------------------------------------------------------------------
subroutine brams2jules(veg,ntype)
  implicit none
  integer, intent(in)              :: ntype
  character (len=80), intent(out)  :: veg(ntype)

  if (ntype == 9 ) then
     veg(1)='06 07 20'       !tJ=1 => BT=broadleaf trees
     veg(2)='04 05 14'       !tJ=2 => NT=needleleaf trees
     veg(3)='15 08 16 17'    !tJ=3 => C3G=C3 (temperate) grass
     veg(4)='09'             !tJ=4 => C4G=C4 (tropical) grass
     veg(5)='11 12 13 18'    !tJ=5 => shrub
     veg(6)='19 21'          !tJ=6 => urban
     veg(7)='00 01'          !tJ=7 => lake=inland water
     veg(8)='03 10'          !tJ=8 => soil=bare soil
     veg(9)='02'             !tJ=9 => ice
  elseif (ntype == 13 ) then
     veg(1)='07 20'         !tJ=1 => BT="Broadleaf"  tropical 
     veg(2)='14'            !tJ=2 => Broadleaf temperada
     veg(3)='06'            !tJ=3 => Broadleaf deciduas
     veg(4)='04'            !tJ=4 => Needle-leaf" "evergreen
     veg(5)='05'            !tJ=5 => "Needle-leaf" deciduas
     veg(6)='15 08 11 16 17'!tJ=6 => Gramineas - C3
     veg(7)='09'            !tJ=7 => Gramineas - C4
     veg(8)='12 18'         !tJ=8 => shrub - Cerrado "evergreen"
     veg(9)='13'            !tJ=9 => shrub - Cerrado deciduo
     veg(10)='19 21'        !tJ=10 => urban
     veg(11)='00 01'        !tJ=11 => lake=inland water
     veg(12)='03 10'        !tJ=12 => bare soil - Solo nu
     veg(13)='02'           !tJ=13 => ice - Gelo
  else
     print*, "ntype=",ntype
     print*, "ATENCAO... ntype <> 9 ou de 13, Deve-se ajustar a subrotina brams2jules em init_frac2brams.f90"
     stop  
  endif
end subroutine brams2jules

!DSM}
