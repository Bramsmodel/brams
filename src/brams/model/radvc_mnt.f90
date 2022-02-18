!----------------------------------------------------------------------!
! Optional advection scheme for CCATT-BRAMS/BRAMS models version 4.2+  !
! Based on Walcek, 2000 (JGR) and Walcek and Aleksic, 1998 (ATENV).    !
! The scheme is highly conservative, monotonic and keeps mass mixing   !
! ratio positive definite. 					       !
! Implemented by Saulo Freitas (saulo.freitas@cptec.inpe.br) @ Jun/2009!
! MPI/Paralelized by L. Flavio/J. Panneta                              !
!----------------------------------------------------------------------!

module monotonic_adv

  use ModParallelEnvironment, only: MsgDump

  use node_mod, only:           &
       ibcon,    &  !intent(in)
       mynum,    &  !intent(in)
       i0,   &  !intent(in)
       j0,   &  !intent(in)
       nodemyp,  &  !intent(in)
       nodemxp,  &  !intent(in)
       nodemzp      !intent(in)

  use grid_dims, only: &
       nzpmax !intent(in)
       
  use mem_grid, only:        &
       dtlt,   & !intent(in)
       time,   &
       ngrids, & !intent(in)
       ngrid,  & !intent(in)
       dzt,    & !intent(in)
       dztn,   & !intent(in)
       grid_g, & !intent(in)
       grid_g, & !intent(in)
       naddsc, & !intent(in)
       hw4   , & !intent(in)
       if_adap,& !intent(in)
       dyncore_flag  !intent(in)

  use mem_basic, only: basic_g  !intent(in)

  use micphys    ,  only: level !intent(in)

  use rconstants ,  only: cp,p00,cv,rgas,cpi   !intent(in)

  use mem_aer1, only: &
       aerosol,    &       !intent(in)
       num_scalar_aer_1st !intent(in)

  use mem_chem1, only: &
       nspecies_transported !intent(in)

  use module_dry_dep, only: &
       dd_sedim,            &
       naer_transported,    &
       sedim_type

  use mem_scratch, only  : scratch  ! only scr1, inout

  use var_tables, only : scalar_tab & ! (var_p = IN, var_t = INOUT)
       ,num_scalar   ! (in)

  use advMessageMod, only: SendMessageI
  use advMessageMod, only: RecvMessageI
  use advMessageMod, only: SendMessageJ
  use advMessageMod, only: RecvMessageJ
  use advMessageMod, only: newM2
  use advMessageMod, only: newM3
  use advMessageMod, only: newIa
  use advMessageMod, only: newIz
  use advMessageMod, only: newJa
  use advMessageMod, only: newJz
  use advMessageMod, only: nRecvI
  use advMessageMod, only: nRecvJ
  use advMessageMod, only: nSendI
  use advMessageMod, only: nSendJ
  use advMessageMod, only: totalrecvi
  use advMessageMod, only: totalsendi
  use advMessageMod, only: totalrecvj
  use advMessageMod, only: totalsendj

  use ModComm, only: i0LGZ
  use ModComm, only: j0LGZ


  use ParLib, only: parf_send_noblock_real
  use ParLib, only: parf_get_noblock_real
  use ParLib, only: parf_wait_any_nostatus
  use ParLib, only: parf_wait_all_nostatus


  use ModNamelistFile, only: NamelistFile

  use ccatt_start, only: &
       ccatt               ! (in)

  implicit none

  private
  public :: advmnt_driver  ! Subroutine
  public :: StoreNamelistFileAtRadvc_mnt ! Subroutine

  ! public names, set by StoreNamelistFileAtRadvc_mnt
  integer, public :: advmnt 
  integer, public :: GhostZoneLength 

  ! module private variables

  ! flow control flags
  integer, parameter :: ON=1,OFF=0
  integer, parameter :: use_true_density  = 1 ! 0= OFF, 1=ON

  ! for theoretical experiments
  integer, parameter :: theor_wind = 0        ! 0= OFF, 1=ON

  ! constants
  real, parameter :: c1 = cv/rgas
  real, parameter :: c2 = p00/rgas

  ! all fields with enlarged ghost zones
  type advmnt_vars
     real,pointer :: u3d(:,:,:)
     real,pointer :: v3d(:,:,:)
     real,pointer :: w3d(:,:,:)
     real,pointer :: vc3d_in(:,:,:)
     real,pointer :: vc3d_out(:,:,:)
     real,pointer :: vc3d_x(:,:,:)
     real,pointer :: vc3d_y(:,:,:)
     real,pointer :: dd0_3d(:,:,:)
     real,pointer :: dd0_3du(:,:,:)
     real,pointer :: dd0_3dv(:,:,:)
     real,pointer :: dd0_3dw(:,:,:)
     real,pointer :: den0_3d(:,:,:)
     real,pointer :: den1_3d(:,:,:)
     real,pointer :: den2_3d(:,:,:)
     real,pointer :: den3_3d(:,:,:)
     real,pointer :: l_dxtW(:,:,:)
     real,pointer :: l_dytW(:,:,:)
     real,pointer :: dxtW(:,:)
     real,pointer :: dytW(:,:)
     real,pointer :: dztW(:)
  end type advmnt_vars

  ! single variable containing all enlarged ghost zone fields
  ! for all grids
  type(advmnt_vars), allocatable :: advmnt_g(:)

  ! advmnt_g initialization flag
  ! advmnt_g should be initialized only once
  integer :: mnt_adv_jnitialized=0 ! 0=not initialized; 1=initialized

  integer :: nSend_i
  integer :: nSend_j
  integer :: nRecv_i
  integer :: nRecv_j
  integer, parameter :: bigdump=1
  real, allocatable :: buffcomm(:,:,:)
  integer :: nRec_i
  integer :: nSnd_i
  integer :: nRec_j
  integer :: nSnd_j
  integer :: bufSendTotalLength_i
  integer :: bufSendTotalLength_j
  integer :: bufREcvTotalLength_i
  integer :: bufRecvTotalLength_j
  real, allocatable :: bufRecv(:)
  real, allocatable :: bufSend(:)

contains





  subroutine advmnt_driver(varn,m1 ,m2 ,m3 ,ia,iz,ja,jz,izu,jzv,mynum)

    integer , intent(in) :: m1
    integer , intent(in) :: m2
    integer , intent(in) :: m3
    integer , intent(in) :: ia
    integer , intent(in) :: iz
    integer , intent(in) :: ja
    integer , intent(in) :: jz
    integer , intent(in) :: izu
    integer , intent(in) :: jzv
    integer , intent(in) :: mynum
    character(len=*),intent(in) :: varn

    !--- local vars
    integer :: n
    integer :: ng
    integer :: mxyzp
    integer :: i
    integer :: j
    integer :: procfile
    integer :: ibegin
    integer :: iend
    integer :: jbegin
    integer :: jend
    integer :: i_scl
    integer :: sori
    integer :: sorj
    integer :: sosi
    integer :: sosj
    integer :: current_aer_ispc
    integer :: current_ndt_z
    integer :: ndt_z(naer_transported)
    real, pointer :: scalarp
    real, pointer :: scalart
    logical  :: IsThisScalarAer =.false.

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(advmnt_driver)**"
    character(len=8) :: str(11)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       call MsgDump(h//" starts; MPI rank fields are dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
       write(str(1),"(i8)") 1+i0
       write(str(2),"(i8)") m2+i0
       write(str(3),"(i8)") 1+j0
       write(str(4),"(i8)") m3+j0
       call MsgDump(h//" fields global indices ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")")
       write(str(1),"(i8)") ia+i0
       write(str(2),"(i8)") iz+i0
       write(str(3),"(i8)") ja+j0
       write(str(4),"(i8)") jz+j0
       call MsgDump(h//" own (no ghost) fields global indices ("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")")
    end if

    if(mnt_adv_jnitialized == OFF) then
       if (mynum == 0) stop 'ADV MNT called with mynum = 0, try np = 2'
       if (dumpLocal) then
          call MsgDump(h//" enter initialization")
       end if

       ! allocates advmnt_g and its fields for all grids
       if (dumpLocal) then
          write(str(1),"(i8)") newM2
          write(str(2),"(i8)") newM3
          call MsgDump(h//" invokes initialize_advmnt passing"//&
               " newM2="//trim(adjustl(str(1)))//&
               " newM3="//trim(adjustl(str(2))))
       end if
       call initialize_advmnt(ngrids, &
            nodemzp(mynum,:), &
            newM2, &
            newM3)

       ! allocates a single buffer for message passing sends
       ! and a single buffer for message passing receives
       sori=maxval(totalrecvi)
       sorj=maxval(totalrecvj)
       sosi=maxval(totalsendi)
       sosj=maxval(totalsendj)
       allocate(bufRecv(max(sori,sorj)))
       allocate(bufSend(max(sosi,sosj)))

       do ng=1,ngrids
          iBegin=newIa(ng)-1
          iEnd=newIz(ng)+1
          jBegin=newJa(ng)-1
          jEnd=newJz(ng)+1
          if (dumpLocal) then
             call MsgDump(h//" invokes initialize_grid_spacings")
          end if
          call initialize_grid_spacings(ng, &
               nodemzp(mynum,ng), &
               nodemxp(mynum,ng), &
               nodemyp(mynum,ng), &
               grid_g(ng)%dxt, &
               grid_g(ng)%dyt, &
               grid_g(ng)%fmapt, &
               grid_g(ng)%rtgt, &
               advmnt_g(ng)%dxtW(iBegin:iEnd,jBegin:jEnd), &
               advmnt_g(ng)%dytW(iBegin:iEnd,jBegin:jEnd), &
               advmnt_g(ng)%dztW )
       end do

       if(use_true_density == OFF) then
          do ng=1,ngrids
             iBegin=newIa(ng)-1
             iEnd=newIz(ng)+1
             jBegin=newJa(ng)-1
             jEnd=newJz(ng)+1
             if (dumpLocal) then
                call MsgDump(h//" invokes initialize_densities")
             end if
             call initialize_densities(nodemzp(mynum,ng),&
                  nodemxp(mynum,ng),nodemyp(mynum,ng) &
                  , basic_g(ng)%dn0     &
                  , basic_g(ng)%dn0u    &
                  , basic_g(ng)%dn0v    &
                  ,advmnt_g(ng)%dd0_3d (1:m1,iBegin:iEnd,jBegin:jEnd)  &
                  ,advmnt_g(ng)%dd0_3du(1:m1,iBegin:iEnd,jBegin:jEnd) &
                  ,advmnt_g(ng)%dd0_3dv(1:m1,iBegin:iEnd,jBegin:jEnd) &
                  ,advmnt_g(ng)%dd0_3dw(1:m1,iBegin:iEnd,jBegin:jEnd) )
          end do
       end if

       mnt_adv_jnitialized= ON
    end if

    mxyzp=m1*m2*m3

    !- This scheme is not applied to advect  U, V, and W
    if (varn .eq. 'V' .or. varn .eq. 'ALL') then
       stop 'not using mnt to advect u,v,w'
    end if
    if (if_adap /= 0) then
       stop 'MNT advection not ready for shaved eta'
    end if
    !
    ndt_z =0 ! integer initialization
    !
    !- Advect  scalars
    iBegin=newIa(ngrid)-1
    iEnd  =newIz(ngrid)+1
    jBegin=newJa(ngrid)-1
    jEnd  =newJz(ngrid)+1

    !- get actual air densities, if using them instead of basic state fields
    if(use_true_density == ON) then
       if (dumpLocal) then
          write(str(1),"(i8)") m1
          write(str(2),"(i8)") m2
          write(str(3),"(i8)") m3
          write(str(4),"(i8)") iBegin
          write(str(5),"(i8)") iEnd
          write(str(6),"(i8)") jBegin
          write(str(7),"(i8)") jEnd
          call MsgDump(h//" invokes get_true_densities")
          call MsgDump(h//" passing input arrays of dim (1:"//&
               trim(adjustl(str(1)))//",1:"//&
               trim(adjustl(str(2)))//",1:"//trim(adjustl(str(3)))//")")
          call MsgDump(h//" and section of output arrays dim"//&
               " (1:"//trim(adjustl(str(1)))//","//&
               trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
               trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")")
       end if
       call get_true_densities(m1,m2,m3,level, &
            basic_g(ngrid)%rtp,     &
            basic_g(ngrid)%rv,      &
            basic_g(ngrid)%pp,      &
            basic_g(ngrid)%pi0,     &
            basic_g(ngrid)%theta,   &
            advmnt_g(ngrid)%dd0_3d (1:m1,iBegin:iEnd,jBegin:jEnd), &
            advmnt_g(ngrid)%dd0_3du(1:m1,iBegin:iEnd,jBegin:jEnd), &
            advmnt_g(ngrid)%dd0_3dv(1:m1,iBegin:iEnd,jBegin:jEnd), &
            advmnt_g(ngrid)%dd0_3dw(1:m1,iBegin:iEnd,jBegin:jEnd)  )

    end if

    !- prepare wind velocities including map factors
    if (dumpLocal) then
          write(str(1),"(i8)") ia
          write(str(2),"(i8)") iz
          write(str(3),"(i8)") ja
          write(str(4),"(i8)") jz
          write(str(5),"(i8)") m1
          write(str(6),"(i8)") m2
          write(str(7),"(i8)") m3
          write(str(8),"(i8)") iBegin
          write(str(9),"(i8)") iEnd
          write(str(10),"(i8)") jBegin
          write(str(11),"(i8)") jEnd
          call MsgDump(h//" invokes prepare_winds on output fields (1:"//&
               trim(adjustl(str(5)))//","//&
               trim(adjustl(str(8)))//":"//trim(adjustl(str(9)))//","//&
               trim(adjustl(str(10)))//":"//trim(adjustl(str(11)))//")"//&
               " using"//&
               " m1="//trim(adjustl(str(5)))//&
               " m2="//trim(adjustl(str(6)))//&
               " m3="//trim(adjustl(str(7)))//&
               " ia="//trim(adjustl(str(1)))//&
               " iz="//trim(adjustl(str(2)))//&
               " ja="//trim(adjustl(str(3)))//&
               " jz="//trim(adjustl(str(4))))
    end if
    call prepare_winds(dtlt,m1,m2,m3,ia,iz,ja,jz     &
         ,basic_g(ngrid)%uc  &
         ,basic_g(ngrid)%up  &
         ,basic_g(ngrid)%vc  &
         ,basic_g(ngrid)%vp  &
         ,basic_g(ngrid)%wc  &
         ,basic_g(ngrid)%wp  &
                                !
         ,grid_g(ngrid)%fmapui &
         ,grid_g(ngrid)%fmapvi &
         ,grid_g(ngrid)%rtgt   &
         ,grid_g(ngrid)%rtgu   &
         ,grid_g(ngrid)%rtgv   &
         ,grid_g(ngrid)%f13t   &
         ,grid_g(ngrid)%f23t   &
                                !
         ,advmnt_g(ngrid)%u3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
         ,advmnt_g(ngrid)%v3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
         ,advmnt_g(ngrid)%w3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
         ,ndt_z                )


    if(theor_wind == on) then
       if (dumpLocal) then
          call MsgDump (h//" invokes prepare_theor_winds")
       end if

       call prepare_theor_winds(dtlt,m1,m2,m3,ia,iz,ja,jz,time     &
            ,advmnt_g(ngrid)%u3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
            ,advmnt_g(ngrid)%v3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
            ,advmnt_g(ngrid)%w3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
            ,grid_g(ngrid)%dxt    &
            ,grid_g(ngrid)%dyt    &
            ,advmnt_g(ngrid)%dd0_3d (1:m1,iBegin:iEnd,jBegin:jEnd)  &
            ,advmnt_g(ngrid)%dd0_3du(1:m1,iBegin:iEnd,jBegin:jEnd) &
            ,advmnt_g(ngrid)%dd0_3dv(1:m1,iBegin:iEnd,jBegin:jEnd) &
            ,advmnt_g(ngrid)%dd0_3dw(1:m1,iBegin:iEnd,jBegin:jEnd) )
    end if
    !- prepare Walcek's air densities
   if (dumpLocal) then
          write(str(1),"(i8)") ia
          write(str(2),"(i8)") iz
          write(str(3),"(i8)") ja
          write(str(4),"(i8)") jz
          write(str(5),"(i8)") m1
          write(str(6),"(i8)") m2
          write(str(7),"(i8)") m3
          write(str(8),"(i8)") iBegin
          write(str(9),"(i8)") iEnd
          write(str(10),"(i8)") jBegin
          write(str(11),"(i8)") jEnd
          call MsgDump(h//" invokes get_Walceks_densities on fields (1:"//&
               trim(adjustl(str(5)))//","//&
               trim(adjustl(str(8)))//":"//trim(adjustl(str(9)))//","//&
               trim(adjustl(str(10)))//":"//trim(adjustl(str(11)))//")"//&
               " using"//&
               " m1="//trim(adjustl(str(5)))//&
               " m2="//trim(adjustl(str(6)))//&
               " m3="//trim(adjustl(str(7)))//&
               " ia="//trim(adjustl(str(1)))//&
               " iz="//trim(adjustl(str(2)))//&
               " ja="//trim(adjustl(str(3)))//&
               " jz="//trim(adjustl(str(4))))
    end if
    call get_Walceks_densities(dtlt,m1,m2,m3 &
         ,advmnt_g(ngrid)%u3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
         ,advmnt_g(ngrid)%v3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
         ,advmnt_g(ngrid)%w3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
         ,advmnt_g(ngrid)%dd0_3d (1:m1,iBegin:iEnd,jBegin:jEnd)  &
         ,advmnt_g(ngrid)%dd0_3du(1:m1,iBegin:iEnd,jBegin:jEnd) &
         ,advmnt_g(ngrid)%dd0_3dv(1:m1,iBegin:iEnd,jBegin:jEnd) &
         ,advmnt_g(ngrid)%dd0_3dw(1:m1,iBegin:iEnd,jBegin:jEnd) &
         ,advmnt_g(ngrid)%den0_3d(1:m1,iBegin:iEnd,jBegin:jEnd) &
         ,advmnt_g(ngrid)%den1_3d(1:m1,iBegin:iEnd,jBegin:jEnd) &
         ,advmnt_g(ngrid)%den2_3d(1:m1,iBegin:iEnd,jBegin:jEnd) &
         ,advmnt_g(ngrid)%den3_3d(1:m1,iBegin:iEnd,jBegin:jEnd) &
         ,advmnt_g(ngrid)%dxtW(iBegin:iEnd,jBegin:jEnd) &
         ,advmnt_g(ngrid)%dytW(iBegin:iEnd,jBegin:jEnd) &
         ,advmnt_g(ngrid)%dztW &
         ,grid_g(ngrid)%dxt    &
         ,grid_g(ngrid)%dyt    )


    if (dumpLocal) then
       call MsgDump (h//" invokes InitialFieldsUpdate exchanging borders on x")
    end if
    call InitialFieldsUpdate(ngrids,m1,m2,m3,newm2(ngrid),newm3(ngrid),ngrid,mynum, &
         nrecvI(ngrid),RecvMessageI(ngrid)%proc,RecvMessageI(ngrid)%tag, &
         RecvMessageI(ngrid)%ia,RecvMessageI(ngrid)%iz,&
         RecvMessageI(ngrid)%ja,RecvMessageI(ngrid)%jz, &
         RecvMessageI(ngrid)%start,RecvMessageI(ngrid)%mSize,TotalRecvI(ngrid), &
         nSendI(ngrid),sendMessageI(ngrid)%proc,sendMessageI(ngrid)%tag, &
         sendMessageI(ngrid)%ia,sendMessageI(ngrid)%iz,&
         sendMessageI(ngrid)%ja,sendMessageI(ngrid)%jz, &
         sendMessageI(ngrid)%start,sendMessageI(ngrid)%mSize,TotalSendI(ngrid))

    if (dumpLocal) then
       call MsgDump (h//" invokes InitialFieldsUpdate exchanging borders on y")
    end if
    call InitialFieldsUpdate(ngrids,m1,m2,m3,newm2(ngrid),newm3(ngrid),ngrid,mynum, &
         nrecvJ(ngrid),RecvMessageJ(ngrid)%proc,RecvMessageJ(ngrid)%tag, &
         RecvMessageJ(ngrid)%ia,RecvMessageJ(ngrid)%iz,&
         RecvMessageJ(ngrid)%ja,RecvMessageJ(ngrid)%jz, &
         RecvMessageJ(ngrid)%start,RecvMessageJ(ngrid)%mSize,TotalRecvJ(ngrid), &
         nSendJ(ngrid),sendMessageJ(ngrid)%proc,sendMessageJ(ngrid)%tag, &
         sendMessageJ(ngrid)%ia,sendMessageJ(ngrid)%iz,&
         sendMessageJ(ngrid)%ja,sendMessageJ(ngrid)%jz, &
         sendMessageJ(ngrid)%start,sendMessageJ(ngrid)%mSize,TotalSendJ(ngrid))


    !- ready to do advection, loop over all scalars
    if(advmnt == 1) then
       i_scl=1                                            !- all scalars
    else if(advmnt == 2) then
       i_scl=num_scalar(ngrid) - NSPECIES_TRANSPORTED +1  !- only chemical + aer species
    else if(advmnt == 3) then
       i_scl=2                                            !- all scalars, but not theta_il
    end if

    !srf- do n=1,num_scalar(ngrid)     ! original
    do n=i_scl,num_scalar(ngrid)

       !- if RK or ABM3 scheme, THP/THC are not transported here
       if (dyncore_flag == 2) then
          if (scalar_tab(n,ngrid)%name == 'THC' .or. &
               scalar_tab(n,ngrid)%name == 'THP') cycle
       end if

       !srf - somente para gases e aerossois
       !     do n=num_scalar(ngrid) - NSPECIES_TRANSPORTED +1,num_scalar(ngrid)
       !      if (scalar_tab(n,ngrid)%name /= 'COP' .and. scalar_tab(n,ngrid)%name /= 'CH4P') cycle
       !          scalar_tab(n,ngrid)%name /= 'O3P'  ) cycle

       !- Aerosol sedimentation
       IsThisScalarAer  = .false.
       current_aer_ispc = 0
       current_ndt_z    = 1
       if(ccatt == 1 .and. aerosol > 0 .and. n >= num_scalar_aer_1st) then
          !srf-  We are going to include sedimentation of aerosols at
          !      vertical advection tendency. It is supposed that scalars
          !      with  N >= num_scalar_aer_1st are _all_ aerosols .
          !
          IsThisScalarAer=.true.
          current_aer_ispc = n - num_scalar_aer_1st + 1
          current_ndt_z    = ndt_z (current_aer_ispc)

       end if

       scalarp => scalar_tab(n,ngrid)%var_p
       scalart => scalar_tab(n,ngrid)%var_t
       if (dumpLocal) then
          call MsgDump (h//" advects scalar "//&
               trim(adjustl(scalar_tab(n,ngrid)%name)))
       end if

       if (dumpLocal) then
          write(str(1),"(i8)") m1
          write(str(2),"(i8)") iBegin
          write(str(3),"(i8)") iEnd
          write(str(4),"(i8)") jBegin
          write(str(5),"(i8)") jEnd
          call MsgDump (h//" atob copies "//&
               trim(adjustl(scalar_tab(n,ngrid)%name))//&
               " with ghost zone 1 into vc3d_in(1:"//&
               trim(adjustl(str(1)))//","//&
               trim(adjustl(str(2)))//":"//trim(adjustl(str(3)))//","//&
               trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//")"//&
               " with ghost zone 3")
       end if
       call atob(mxyzp,scalarp,advmnt_g(ngrid)%vc3d_in(1:m1,iBegin:iEnd,jBegin:jEnd))

       if (dumpLocal) then
          write(str(1),"(i8)") ia
          write(str(2),"(i8)") iz
          write(str(3),"(i8)") ja
          write(str(4),"(i8)") jz
          write(str(5),"(i8)") m1
          write(str(6),"(i8)") m2
          write(str(7),"(i8)") m3
          call MsgDump(h//" invokes advect_mnt with"//&
               " m1="//trim(adjustl(str(5)))//&
               " m2="//trim(adjustl(str(6)))//&
               " m3="//trim(adjustl(str(7)))//&
               " ia="//trim(adjustl(str(1)))//&
               " iz="//trim(adjustl(str(2)))//&
               " ja="//trim(adjustl(str(3)))//&
               " jz="//trim(adjustl(str(4))))
       end if
       call advect_mnt(ngrid,m1,m2,m3,ia,iz,ja,jz,dtlt,mynum,n, &
            current_aer_ispc,current_ndt_z,IsThisScalarAer)

       if (dumpLocal) then
          write(str(1),"(i8)") ia
          write(str(2),"(i8)") iz
          write(str(3),"(i8)") ja
          write(str(4),"(i8)") jz
          write(str(5),"(i8)") m1
          write(str(6),"(i8)") m2
          write(str(7),"(i8)") m3
          write(str(8),"(i8)") iBegin
          write(str(9),"(i8)") iEnd
          write(str(10),"(i8)") jBegin
          write(str(11),"(i8)") jEnd
          call MsgDump(h//" invokes advtndc to increment scalart with vc3d_out-scalarp on (1:"//&
               trim(adjustl(str(5)))//","//&
               trim(adjustl(str(8)))//":"//trim(adjustl(str(9)))//","//&
               trim(adjustl(str(10)))//":"//trim(adjustl(str(11)))//")")
!!$               " using"//&
!!$               " m1="//trim(adjustl(str(5)))//&
!!$               " m2="//trim(adjustl(str(6)))//&
!!$               " m3="//trim(adjustl(str(7)))//&
!!$               " ia="//trim(adjustl(str(1)))//&
!!$               " iz="//trim(adjustl(str(2)))//&
!!$               " ja="//trim(adjustl(str(3)))//&
!!$               " jz="//trim(adjustl(str(4))))
       end if
       call advtndc(m1,m2,m3,ia,iz,ja,jz        &
            ,scalarp  ,advmnt_g(ngrid)%vc3d_out(1:m1,iBegin:iEnd,jBegin:jEnd)  &
            ,scalart  ,dtlt,mynum        )
    end do
  end subroutine advmnt_driver





  subroutine initialize_advmnt(ngrids, mmzp, newM2, newM3)

    ! allocates module variable advmnt_g of type advmnt_vars
    ! containing all fields with enlarged ghost zones

    ! allocates each field of advmnt_g for each grid and initializes to zero

    integer, intent(in) :: ngrids
    integer, intent(in) :: mmzp(ngrids)
    integer, intent(in) :: newM2(ngrids)
    integer, intent(in) :: newM3(ngrids)

    integer :: ng
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(initialize_advmnt)**"
    character(len=8) :: str(10)

    if(allocated(advmnt_g)) then
       print *,'Error in initialize_advmnt, sub: radvc_mnt: advmnt_g already allocated!'
       stop 1000
    end if

    allocate (advmnt_g(ngrids))
    do ng=1,ngrids
       allocate(advmnt_g(ng)%u3d    (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%u3d=0.
       allocate(advmnt_g(ng)%v3d    (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%v3d=0.
       allocate(advmnt_g(ng)%w3d    (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%w3d=0.
       allocate(advmnt_g(ng)%dd0_3d (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%dd0_3d =0.
       allocate(advmnt_g(ng)%dd0_3du(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%dd0_3du=0.
       allocate(advmnt_g(ng)%dd0_3dv(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%dd0_3dv=0.
       allocate(advmnt_g(ng)%dd0_3dw(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%dd0_3dw=0.
       allocate(advmnt_g(ng)%den0_3d(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%den0_3d=0.
       allocate(advmnt_g(ng)%den1_3d(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%den1_3d=0.
       allocate(advmnt_g(ng)%den2_3d(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%den2_3d=0.
       allocate(advmnt_g(ng)%den3_3d(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%den3_3d=0.
       allocate(advmnt_g(ng)%l_dxtW (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%l_dxtW=0.
       allocate(advmnt_g(ng)%l_dytW (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%l_dytW=0.
       allocate(advmnt_g(ng)%dxtW   (         newM2(ng),newM3(ng))); advmnt_g(ng)%dxtW=0.
       allocate(advmnt_g(ng)%dytW   (         newM2(ng),newM3(ng))); advmnt_g(ng)%dytW=0.
       allocate(advmnt_g(ng)%dztW   (mmzp(ng)                    )); advmnt_g(ng)%dztW=0.
       allocate(advmnt_g(ng)%vc3d_in (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%vc3d_in =0.
       allocate(advmnt_g(ng)%vc3d_out(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%vc3d_out=0.

       if (dumpLocal) then
          write(str(1),"(i8)") ng
          write(str(2),"(i8)") mmzp(ng)
          write(str(3),"(i8)") newM2(ng)
          write(str(4),"(i8)") newM3(ng)
          call MsgDump(h//" grid "//trim(adjustl(str(1)))//" allocates advmnt_g fields:")
          call MsgDump(h//" allocates u3d("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates v3d("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates w3d("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates dd0_3d("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates dd0_3du("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates dd0_3dv("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates dd0_3dw("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates den0_3d("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates den1_3d("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates den2_3d("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates den3_3d("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates l_dxtW("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates l_dytW("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates dxtW("//&
               trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates dytW("//&
               trim(adjustl(str(3)))//")")
          call MsgDump(h//" allocates dztW("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates vc3d_in("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
          call MsgDump(h//" allocates vc3d_out("//&
               trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//","//trim(adjustl(str(4)))//")")
       end if

    end do
  end subroutine initialize_advmnt





  subroutine Deinitialize_advmnt(ngrids, mmzp,mmxp,mmyp)

    integer , intent(in) :: ngrids
    integer , intent(in) :: mmxp(ngrids)
    integer , intent(in) :: mmyp(ngrids)
    integer , intent(in) :: mmzp(ngrids)

    integer :: ng

    do ng=1,ngrids

       deallocate(advmnt_g(ng)%u3d)
       deallocate(advmnt_g(ng)%v3d)
       deallocate(advmnt_g(ng)%w3d)

       deallocate(advmnt_g(ng)%dd0_3d )
       deallocate(advmnt_g(ng)%dd0_3du)
       deallocate(advmnt_g(ng)%dd0_3dv)
       deallocate(advmnt_g(ng)%dd0_3dw)

       deallocate(advmnt_g(ng)%den0_3d)
       deallocate(advmnt_g(ng)%den1_3d)
       deallocate(advmnt_g(ng)%den2_3d)
       deallocate(advmnt_g(ng)%den3_3d)

       deallocate(advmnt_g(ng)%dxtW)
       deallocate(advmnt_g(ng)%dytW)
       deallocate(advmnt_g(ng)%dztW)

       deallocate(advmnt_g(ng)%l_dxtW)
       deallocate(advmnt_g(ng)%l_dytW)

       deallocate(advmnt_g(ng)%vc3d_in )
       deallocate(advmnt_g(ng)%vc3d_out)
       !deallocate(advmnt_g(ng)%vc3d_x)
       !deallocate(advmnt_g(ng)%vc3d_y)

    end do
    deallocate (advmnt_g)
  end subroutine Deinitialize_advmnt





  subroutine initialize_densities(m1,m2,m3,dn0,dn0u,dn0v &
       ,dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw )

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    real, intent(in) :: dn0(m1,m2,m3)
    real, intent(in) :: dn0u(m1,m2,m3)
    real, intent(in) :: dn0v(m1,m2,m3)
    real, intent(out) :: dd0_3d(m1,m2,m3)
    real, intent(out) :: dd0_3du(m1,m2,m3)
    real, intent(out) :: dd0_3dv(m1,m2,m3)
    real, intent(out) :: dd0_3dw(m1,m2,m3)

    ! local var
    integer i,j,k

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(initialize_densities)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       call MsgDump(h//" set values of"//&
            " dd0_3d("//trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//")"//&
            ", dd0_3du("//trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//")"//&
            ", dd0_3dv("//trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//")"//&
            ", dd0_3dw("//trim(adjustl(str(1)))//","//trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//")")
    end if

    dd0_3d (:,:,:)=  dn0 (:,:,:)
    dd0_3du(:,:,:)=  dn0u(:,:,:)
    dd0_3dv(:,:,:)=  dn0v(:,:,:)
    do j = 1,m3
       do i = 1,m2
          do k = 1,m1-1
             dd0_3dw(k,i,j) = 0.5*(dn0(k,i,j) +dn0(k+1,i,j))
          end do
          dd0_3dw(m1,i,j)=dd0_3dw(m1-1,i,j)
       end do
    end do
  end subroutine initialize_densities





  subroutine initialize_grid_spacings(ng,m1,m2,m3,&
       dxt,dyt,fmapt,rtgt, &
       dxtW,dytW,dztW)

    ! computes new grid spacing on x, y and z
    ! for Walcek monotonic advection

    integer, intent(in) :: ng
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    real, intent(in) :: dxt(m2,m3)
    real, intent(in) :: dyt(m2,m3)
    real, intent(in) :: fmapt(m2,m3)
    real, intent(in) :: rtgt(m2,m3)
    real, intent(out) :: dxtW(m2,m3)
    real, intent(out) :: dytW(m2,m3)
    real, intent(out) :: dztW(m1)

    ! local var
    integer i,j,k
    real rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(initialize_grid_spacings)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       call MsgDump(h//" set values of"//&
            " dxtW("//trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//")"//&
            ", dytW("//trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//")"//&
            ", dztW("//trim(adjustl(str(1)))//")")
    end if

    do j = 1,m3
       do i = 1,m2
          rtgti = 1. / rtgt(i,j)

          !- at init/rams_grid.f90:
          !     dxt(i,j)=fmapt(i,j)/(xmn(i,ngrid)-xmn(i-1,ngrid))
          !     dyt(i,j)=fmapt(i,j)/(ymn(j,ngrid)-ymn(j-1,ngrid))

          dxtW(i,j) = 1./(dxt(i,j) * fmapt(i,j) * rtgti)
          dytW(i,j) = 1./(dyt(i,j) * fmapt(i,j) * rtgti)
       end do
    end do
    do k = 1,m1
       !- at init/gridset.f90:
       !  dztn(k,ifm) = 1. / (zmn(k,ifm) - zmn(k-1,ifm))
       ! Por que o Jacobiano nao depende de Z, o dztw depende somente
       ! de z.
       !dztW(k,i,j) = 1./ ( dzt(k) * rtgti * fmapt(i,j)**2 )
       dztW(k)	 = 1./ ( dztn(k,ng) ) !
    end do
  end subroutine initialize_grid_spacings





  subroutine get_true_densities(m1,m2,m3,level,rtp,rv,pp,pi0,theta &
       ,dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw )

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: level
    real, intent(in) :: rtp(m1,m2,m3)
    real, intent(in) :: rv(m1,m2,m3)
    real, intent(in) :: pp(m1,m2,m3)
    real, intent(in) :: pi0(m1,m2,m3)
    real, intent(in) :: theta(m1,m2,m3)
    real, intent(out) :: dd0_3d(:,:,:)
    real, intent(out) :: dd0_3du(m1,m2,m3)
    real, intent(out) :: dd0_3dv(m1,m2,m3)
    real, intent(out) :: dd0_3dw(m1,m2,m3)

    ! local var
    integer i,j,k
    real c3

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(get_true_densities)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") lbound(dd0_3d,1)
       write(str(5),"(i8)") ubound(dd0_3d,1)
       write(str(6),"(i8)") lbound(dd0_3d,2)
       write(str(7),"(i8)") ubound(dd0_3d,2)
       write(str(8),"(i8)") lbound(dd0_3d,3)
       write(str(9),"(i8)") ubound(dd0_3d,3)
       write(str(10),"(i8)") level
       call MsgDump(h//" declared ouput array"//&
            " dd0_3d("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//","//&
            trim(adjustl(str(8)))//":"//trim(adjustl(str(9)))//")")
       call MsgDump(h//" declared ouput array"//&
            " dd0_3du("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
       call MsgDump(h//" declared ouput array"//&
            " dd0_3dv("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
       call MsgDump(h//" declared ouput array"//&
            " dd0_3dw("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
       call MsgDump(h//" level="//trim(adjustl(str(10))))
    end if

    c3 = c2 * (cpi**c1)

    !- true air density at points "T"

    if( level == 0 ) then
       dd0_3d(:,:,:) = (c3/theta(:,:,:))*(pi0(:,:,:)+pp(:,:,:))**c1
    else
       do j = 1,m3
          do i = 1,m2
             do k = 1,m1
                dd0_3d(k,i,j) = (c3/theta(k,i,j))* (1. + rtp(k,i,j))/ &
                     (1. + 1.61*rv(k,i,j))*(pi0(k,i,j)+pp(k,i,j))**c1
             end do
          end do
       end do
    end if

    !- true air density at points "U", "V" and "W":

    call fill_dn0uv(m1,m2,m3,dd0_3d,dd0_3du,dd0_3dv)

    do j = 1,m3
       do i = 1,m2
          do k = 1,m1-1
             dd0_3dw(k,i,j) = 0.5*(dd0_3d(k,i,j) + dd0_3d(k+1,i,j))
          end do
          dd0_3dw(m1,i,j)=dd0_3dw(m1-1,i,j)
       end do
    end do
  end subroutine get_true_densities





  subroutine prepare_winds(dtlt,m1,m2,m3,ia,iz,ja,jz &
       ,uc,up,vc,vp,wc,wp &
       ,fmapui &
       ,fmapvi &
       ,rtgt   &
       ,rtgu   &
       ,rtgv   &
       ,f13t   &
       ,f23t   &
       ,u3d,v3d,w3d &
       ,ndt_z                )

    integer , intent(in) :: m1,m2,m3,ia,iz,ja,jz
    real    , intent(in) :: dtlt
    real,dimension(m1,m2,m3),intent(in) :: uc,up,vc,vp,wc,wp
    real,dimension(m2,m3)   ,intent(in) :: rtgt,rtgu,rtgv,fmapui,fmapvi,f13t,f23t

    real,dimension(m1,m2,m3),intent(out)::u3d,v3d,w3d

    !- aerosol sedimentation
    integer, dimension(naer_transported) , intent(inout) :: ndt_z

    !- local var
    !real   dtlto2
    integer jm,jp,im,ip , ispc
    integer i,j,k
    real :: cx1,cx2,rtgti,dum(m1)

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(prepare_winds)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       call MsgDump(h//" computes u3d, v3d and w3d just at a section"//&
            " restricted to the original ghost zone of 1")
    end if
    ! dtlto2 = .5


    ! u3d, u3d, and w3d are input as the velocity components (averaged
    ! between past and current time levels) times dtlt.
    do j=1,m3
       do i = 1,m2
          do k = 1,m1

             w3d(k,i,j) = ( wc(k,i,j) + wp(k,i,j) )*0.5
             u3d(k,i,j) = ( uc(k,i,j) + up(k,i,j) )*0.5
             v3d(k,i,j) = ( vc(k,i,j) + vp(k,i,j) )*0.5

          end do
       end do
    end do
    ! after this point w3d is the cartesian vertical velocity


    !return ! for pure cartesian coordinates

    ! here w3d is the cartesian vertical velocity

    ! Add contribution to w3d from horiz winds crossing sloping sigma surfaces,
    ! and include 1/rtgt factor in w3d
    do j = 1,m3
       jm = max(1,j-1)
       jp = min(m3,j+1)
       do i = 1,m2
          im = max(1,i-1)
          ip = min(m2,i+1)
          rtgti = 1. / rtgt(i,j)

          do k = 1,m1-1
             w3d(k,i,j) = &
                  ( &
                  (u3d(k,i,j)+u3d(k+1,i,j)+u3d(k,im,j)+u3d(k+1,im,j)) * f13t(i,j) + &
                  (v3d(k,i,j)+v3d(k+1,i,j)+v3d(k,i,jm)+v3d(k+1,i,jm)) * f23t(i,j)  &
                  ) * hw4(k) + w3d(k,i,j) * rtgti
          end do
       end do
    end do
    !- after this point w3d is the sigma_z velocity

    !- including map factors on U,V:
    do j = 1,m3
       do i = 1,m2
          cx1 = fmapui(i,j) * rtgu(i,j)
          cx2 = fmapvi(i,j) * rtgv(i,j)
          do k = 1,m1-1
             u3d(k,i,j) = u3d(k,i,j) * cx1
             v3d(k,i,j) = v3d(k,i,j) * cx2
          end do
       end do
    end do
    !-----------------------------------------
    !- control for aerosol sedimentation
    if(aerosol > 0 .and. naer_transported > 0) then
       ! very crude estimation of CFL violation and fix for the number of sub-timesteps
       ! for large particles
       do ispc=1,naer_transported
          ndt_z(ispc)=ceiling(maxval(abs(dd_sedim(ispc,ngrid)%v_sed_part))*dtlt*maxval(dzt(1:m1)))
       end do
    end if
    !- end of aerosol sedimentation
  end subroutine prepare_winds





  subroutine get_Walceks_densities(dt,m1,m2,m3,u3d,v3d,w3d &
       ,dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw &
       ,den0_3d,den1_3d,den2_3d,den3_3d &
       ,dxtW,dytW,dztW,dxt,dyt)

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    real, intent(in) :: dt
    real, intent(in) :: dztW(m1)
    real, intent(in) :: dxtW(m2,m3)
    real, intent(in) :: dytW(m2,m3)
    real, intent(in) :: dxt(m2,m3)
    real, intent(in) :: dyt(m2,m3)
    real, intent(in) :: u3d(m1,m2,m3)
    real, intent(in) :: v3d(m1,m2,m3)
    real, intent(in) :: w3d(m1,m2,m3)
    real, intent(in) :: dd0_3d(m1,m2,m3)
    real, intent(in) :: dd0_3du(m1,m2,m3)
    real, intent(in) :: dd0_3dv(m1,m2,m3)
    real, intent(in) :: dd0_3dw(m1,m2,m3)
    real, intent(out) :: den0_3d(m1,m2,m3)
    real, intent(out) :: den1_3d(m1,m2,m3)
    real, intent(out) :: den2_3d(m1,m2,m3)
    real, intent(out) :: den3_3d(m1,m2,m3)


    ! local var
    integer i,j,k

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(get_Walceks_densities)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       call MsgDump(h//" computes den0_3d, den1_3d, den2_3d and den3_3d"//&
            " just at a section restricted to the original ghost zone of 1")
    end if

    do  j=m3,2,-1
       do  i=2,m2
          do k = 2,m1
             den0_3d(k,i,j)=dd0_3d(k,i,j)
             den1_3d(k,i,j)=den0_3d(k,i,j)- dt/dxtW(i,j)*&
                  (dd0_3du(k,i,j)*u3d(k,i,j)-dd0_3du(k,i-1,j)*u3d(k,i-1,j))
             den2_3d(k,i,j)=den1_3d(k,i,j)- dt/dytW(i,j)*&
                  (dd0_3dv(k,i,j)*v3d(k,i,j)-dd0_3dv(k,i,j-1)*v3d(k,i,j-1))
             den3_3d(k,i,j)=den2_3d(k,i,j)- dt/dztW(k)  *&
                  (dd0_3dw(k,i,j)*w3d(k,i,j)-dd0_3dw(k-1,i,j)*w3d(k-1,i,j))
          end do
       end do
    end do
    !srf- BC for den3_3d
    den3_3d(:,1,:)=den3_3d(:,2,:)
    den3_3d(:,:,1)=den3_3d(:,:,2)
  end subroutine get_Walceks_densities





  subroutine advect_mnt(ngrid,m1,m2,m3,ia,iz,ja,jz,dt,mynum,n,&
       current_aer_ispc,current_ndt_z,IsThisScalarAer)

    integer , intent(in) :: m1,ngrid
    integer , intent(in) :: m2
    integer , intent(in) :: m3
    integer , intent(in) :: ia
    integer , intent(in) :: iz
    integer , intent(in) :: ja
    integer , intent(in) :: jz,n
    integer , intent(in) :: mynum
    real    , intent(in) :: dt
    integer , intent(in) :: current_ndt_z,current_aer_ispc
    logical , intent(in) :: IsThisScalarAer
    !- local var
    !REAL,DIMENSION(m1)               :: dxx
    !REAL,DIMENSION(m2,m3)            :: dxy
    real masscon,initialmass,vol
    integer nrec,itz
    integer ibegin,iend,jbegin,jend
    !- type of sedimentation scheme (0= Walcek, 1=upwind)
    integer , parameter :: iupwind = 0
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(advect_mnt)**"
    character(len=8) :: str(10)
    
    iBegin= newIa(ngrid)-1
    iEnd  = newIz(ngrid)+1
    jBegin= newJa(ngrid)-1
    jEnd  = newJz(ngrid)+1

    !--- do X-advection
    if (dumpLocal) then
       call MsgDump(h//" update borders of vc3d_in for x advection")
    end if
    call UpdateBorders(m1, newm2(ngrid), newm3(ngrid),advmnt_g(ngrid)%vc3d_in, &
         nrecvI(ngrid), RecvMessageI(ngrid)%proc, RecvMessageI(ngrid)%tag, &
         RecvMessageI(ngrid)%ia, RecvMessageI(ngrid)%iz, &
         RecvMessageI(ngrid)%ja, RecvMessageI(ngrid)%jz, &
         RecvMessageI(ngrid)%start, RecvMessageI(ngrid)%mSize, TotalRecvI(ngrid), &
         nSendI(ngrid), sendMessageI(ngrid)%proc, sendMessageI(ngrid)%tag, &
         sendMessageI(ngrid)%ia, sendMessageI(ngrid)%iz, &
         sendMessageI(ngrid)%ja, sendMessageI(ngrid)%jz, &
         sendMessageI(ngrid)%start, sendMessageI(ngrid)%mSize, TotalSendI(ngrid))

    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3d_X to advect vc3d_in, storing result in vc3d_out")
    end if
    call Advec3d_X(m1,newM2(ngrid),newM3(ngrid),2,newM2(ngrid)-1,2,newM3(ngrid)-1 &
         ,advmnt_g(ngrid)%vc3d_in                             &
         ,advmnt_g(ngrid)%u3d,advmnt_g(ngrid)%den0_3d         &
         ,advmnt_g(ngrid)%den1_3d,dt,advmnt_g(ngrid)%dxtW     &
         ,advmnt_g(ngrid)%dd0_3du                             &
         ,advmnt_g(ngrid)%vc3d_out    ,mynum )

    !--- do Y-advection

    if (dumpLocal) then
       call MsgDump(h//" update borders of vc3d_out for y advection")
    end if
    call UpdateBorders(m1, newm2(ngrid), newm3(ngrid),advmnt_g(ngrid)%vc3d_out, &
         nrecvJ(ngrid), RecvMessageJ(ngrid)%proc, RecvMessageJ(ngrid)%tag, &
         RecvMessageJ(ngrid)%ia, RecvMessageJ(ngrid)%iz, &
         RecvMessageJ(ngrid)%ja, RecvMessageJ(ngrid)%jz, &
         RecvMessageJ(ngrid)%start, RecvMessageJ(ngrid)%mSize, TotalRecvJ(ngrid), &
         nSendJ(ngrid), sendMessageJ(ngrid)%proc, sendMessageJ(ngrid)%tag, &
         sendMessageJ(ngrid)%ia, sendMessageJ(ngrid)%iz, &
         sendMessageJ(ngrid)%ja, sendMessageJ(ngrid)%jz, &
         sendMessageJ(ngrid)%start, sendMessageJ(ngrid)%mSize, TotalSendJ(ngrid))

    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3d_Y to advect vc3d_out, storing result in vc3d_in")
    end if
    call Advec3d_Y(m1,newM2(ngrid),newM3(ngrid),2,newM2(ngrid)-1,2,newM3(ngrid)-1 &
         ,advmnt_g(ngrid)%vc3d_out                                        &
         ,advmnt_g(ngrid)%v3d,advmnt_g(ngrid)%den1_3d                     &
         ,advmnt_g(ngrid)%den2_3d,dt,advmnt_g(ngrid)%dytW                 &
         ,advmnt_g(ngrid)%dd0_3dv                                         &
         ,advmnt_g(ngrid)%vc3d_in  ,mynum )

    !--- do k-advection
    if (dumpLocal) then
       call MsgDump(h//" invoke Advec3d_Z to advect vc3d_in, storing result in vc3d_out")
    end if
    call Advec3d_Z(m1,newM2(ngrid),newM3(ngrid),ibegin,iend,jbegin,jend &
         ,advmnt_g(ngrid)%vc3d_in                                   &
         ,advmnt_g(ngrid)%w3d,advmnt_g(ngrid)%den2_3d               &
         ,advmnt_g(ngrid)%den3_3d,dt,advmnt_g(ngrid)%dztW           &
         ,advmnt_g(ngrid)%dd0_3dw                                   &
         ,advmnt_g(ngrid)%vc3d_out  ,mynum )


    !- aerosol section to include sedimentation
    !- the sedimentation process is done using pure cartesian coordinates
    !- so, all sedimentation velocities are treat as cartesian vertical velocities
    !- which are positive downwards.
    if (dumpLocal) then
       write(str(1),"(i8)") aerosol
       write(str(2),"(l)") IsThisScalarAer
       call MsgDump(h//" aerosol="//trim(adjustl(str(1)))//&
            "; IsThisScalarAer="//trim(adjustl(str(2))))
    end if
    if(aerosol > 0 .and. IsThisScalarAer) then

       !-srf introducing a time-splitting for aerosol sedimentation

       if (dumpLocal) then
          write(str(1),"(i8)") iupwind
          call MsgDump(h//" iupwind="//trim(adjustl(str(1))))
       end if
       if(iupwind == 0 ) then
          ! - Walcek method
          ! this routine works _only_ for mass concentration or density (kg/m3)
          ! converting mixing ratio (kg/kg) to density (kg/m3)
          advmnt_g(ngrid)%vc3d_in(:,:,:)=advmnt_g(ngrid)%vc3d_out(:,:,:) * advmnt_g(ngrid)%den0_3d(:,:,:)

          !- do time splitting for aerosols with large fall velocities
          do itz=1,current_ndt_z
             call Advec3d_Z_sedim(m1,m2,m3,ia,iz,ja,jz                        &
                  ,advmnt_g(ngrid)%vc3d_in(:,iBegin:iEnd,jBegin:jEnd)	 &
                  ,dd_sedim(current_aer_ispc,ngrid)%v_sed_part          & !fall velocity
                  ,dt/float(current_ndt_z)                              & !subtimestep
                  ,dzt(1:m1),grid_g(ngrid)%rtgt	                 &
                  ,advmnt_g(ngrid)%vc3d_out(:,iBegin:iEnd,jBegin:jEnd)  &
                  ,mynum )

             ! copy output to input array for the next sup-timestep
             if(itz < current_ndt_z) advmnt_g(ngrid)%vc3d_in(:,:,:)=advmnt_g(ngrid)%vc3d_out(:,:,:)

          end do
          ! converting back mass concentration to mixing ratio
          advmnt_g(ngrid)%vc3d_out(:,:,:)=&
               advmnt_g(ngrid)%vc3d_out(:,:,:)/&
               advmnt_g(ngrid)%den0_3d(:,:,:)

       else if(iupwind == 1 ) then
          ! - upwind method
          !- do time splitting for aerosols with large fall velocities
          do itz=1,current_ndt_z
             call Advec3d_Z_sedim_upw(m1,m2,m3,ia,iz,ja,jz                          &
                  ,dd_sedim(current_aer_ispc,ngrid)%v_sed_part          & !fall velocity
                  ,dt/float(current_ndt_z)                              & !subtimestep
                  ,dzt(1:m1),grid_g(ngrid)%rtgt	                                 &
                  ,advmnt_g(ngrid)%vc3d_out(:,iBegin:iEnd,jBegin:jEnd)  &
                  ,mynum )

          end do
       end if
    end if
  end subroutine advect_mnt





  subroutine prepare_theor_winds(dtlt,m1,m2,m3,ia,iz,ja,jz,time &
       ,u3d,v3d,w3d&
       ,dxt,dyt &
       ,dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw )


    integer , intent(in) :: m1,m2,m3,ia,iz,ja,jz
    real, intent(in) :: dtlt,time
    real, dimension(m2,m3)   , intent(in) :: dxt,dyt

    real,dimension(m1,m2,m3),intent(out)::u3d,v3d,w3d     &
         ,dd0_3d ,dd0_3du &
         ,dd0_3dv,dd0_3dw

    !- local var
    real   :: dtlto2
    integer :: i,j,k
    real  :: ai0s  =  25.0
    real  :: aj0s  =  50.0
    real  :: umx   =   80.0
    real,parameter :: pii   =   3.141592653589793
    real    :: umax  =   0.0
    real    :: anrev,curnt,rx,xa,ilop,iwndty,nrec,ya
    real    :: periodo  =   6.*3600.

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(prepare_theor_winds)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       call MsgDump(h//" is invoked")
    end if

    dtlto2 =  10.!*dtlt

    !WRITE(6,*) ' Wind fields?  0-rotating; or  1-divergent winds'
    !iwndty = 0  ! 0-rotating
    !iwndty = 1  ! 1-divergent winds
    iwndty = 2  ! 1-divergent winds

    if(iwndty==1) ai0s= 50.5
    ilop= ai0s-21.  ! needed for printouts
    ! Define wind fields (rotation or divergent) and initial mixing ratios
    !  Cone at (25,50) for rotating winds; Cone at (50,50) divergent winds
    do k = 1,m1
       do  j=m3,1,-1
          do  i=1,m2

             dd0_3d (k,i,j)=1.
             dd0_3du(k,i,j)=1.
             dd0_3dv(k,i,j)=1.
             dd0_3dw(k,i,j)=1.

             u3d(k,i,j)= -2.*umx*(real(j)-real(110)/2.-.5)/real(110)
             v3d(k,i,j)=  2.*umx*(real(i)-real(100)/2.-.5)/real(100)
             w3d(k,i,j)= 0.

             if(iwndty==1) then
                xa=pii/25.
                if(J>0) u3d(k,i,j)=umx*sin(xa*real(i))*sin(xa*(real(j)))
                if(I>0) v3d(k,i,j)=umx*cos(xa*(real(i)-.5))*cos(XA*(real(j)+.5))
             end if


             if(iwndty==2) then
                xa=pii/100. ! m3=m2
                if(J>0) u3d(k,i,j)=  umx* (sin(xa*real(i)))**2 *sin(2*xa*(real(j)))*cos(pii*time/periodo)
                if(I>0) v3d(k,i,j)=- umx* (sin(xa*real(j)))**2 *sin(2*xa*(real(i)))*cos(pii*time/periodo)
             end if


             if(iwndty==3) then
                xa=pii/100. ! m3=m2
                ya=50.
                if(J>0) u3d(k,i,j)= -   umx* (sin(    xa*real(i)))**2 * sin(2.*xa*(real(j)-ya)) *cos(pii*time/periodo)
                if(I>0) v3d(k,i,j)= 0.5*umx* (sin(2.* xa*real(i)))    * cos(   xa*(real(j)-ya)) *cos(pii*time/periodo)
             end if

             umax= max(abs(u3d(k,i,j)),abs(v3d(k,i,j)),umax)
             rx= sqrt((real(i)-ai0s)**2.+(real(j)-aj0s)**2.)

          end do !i
       end do !j
    end do !k
  end subroutine prepare_theor_winds





  subroutine UpdateBorders(m1, m2, m3, field, &
       nRecv, procRecv_ext, tagRecv, iaRecv, izRecv, jaRecv, jzRecv, &
       bufRecvStart, bufRecvLength, bufRecvTotalLength, &
       nSend, procSend_ext, tagSend, iaSend, izSend, jaSend, jzSend, &
       bufSendStart, bufSendLength, bufSendTotalLength)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    real,    intent(inout) :: field(m1,m2,m3)
    integer, intent(in) :: nRecv
    integer, intent(in) :: procRecv_ext(nRecv)
    integer, intent(in) :: tagRecv(nRecv)
    integer, intent(in) :: iaRecv(nRecv)
    integer, intent(in) :: izRecv(nRecv)
    integer, intent(in) :: jaRecv(nRecv)
    integer, intent(in) :: jzRecv(nRecv)
    integer, intent(in) :: bufRecvStart(nRecv)
    integer, intent(in) :: bufRecvLength(nRecv)
    integer, intent(in) :: bufRecvTotalLength
    integer, intent(in) :: nSend
    integer, intent(in) :: procSend_ext(nSend)
    integer, intent(in) :: tagSend(nSend)
    integer, intent(in) :: iaSend(nSend)
    integer, intent(in) :: izSend(nSend)
    integer, intent(in) :: jaSend(nSend)
    integer, intent(in) :: jzSend(nSend)
    integer, intent(in) :: bufSendStart(nSend)
    integer, intent(in) :: bufSendLength(nSend)
    integer, intent(in) :: bufSendTotalLength

    integer :: procRecv(nRecv)
    integer :: procSend(nSend)

    integer :: i, j, iCnt, i1, i2, cnt, iRecv, iSend, ierr, iRecS
    integer :: reqRecv(nRecv)
    integer :: reqSend(nSend), recNum

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(UpdateBorders)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       call MsgDump(h//" starts with"//&
            " m1="//trim(adjustl(str(1)))//&
            " m2="//trim(adjustl(str(2)))//&
            " m3="//trim(adjustl(str(3))))
    end if

    procRecv=procRecv_ext-1
    procSend=procSend_ext-1
    do iRecv = 1, nRecv
       if (dumpLocal) then
          write(str(1),"(i8)") bufRecvLength(iRecv)
          write(str(2),"(i8)") procRecv(iRecv)
          write(str(3),"(i8)") iRecv
          write(str(4),"(i8)") tagRecv(iRecv)
          call MsgDump(h//" dispatch recv #"//trim(adjustl(str(3)))//&
               " with tag "//trim(adjustl(str(4)))//&
               " of "//trim(adjustl(str(1)))//&
               " reals from MPI proc "//trim(adjustl(str(2))))
       end if
       call parf_get_noblock_real(bufRecv(bufRecvStart(iRecv)), &
            bufRecvLength(iRecv),procRecv(iRecv) , &
            tagRecv(iRecv), reqRecv(iRecv))
    end do

    do iSend = 1, nSend
       i1 = bufSendStart(iSend)
       iCnt = bufSendStart(iSend)
       i2 = bufSendLength(iSend)
       do j = jaSend(iSend), jzSend(iSend)
          do i = iaSend(iSend), izSend(iSend)
             bufSend(iCnt:iCnt+m1-1) = field(1:m1,i,j)
             iCnt = iCnt+m1
          end do
       end do
       if (dumpLocal) then
          write(str(1),"(i8)") iaSend(iSend)+i0LGZ
          write(str(2),"(i8)") izSend(iSend)+i0LGZ
          write(str(3),"(i8)") jaSend(iSend)+j0LGZ
          write(str(4),"(i8)") jzSend(iSend)+j0LGZ
          write(str(5),"(i8)") i2
          write(str(6),"(i8)") procSend(iSend)
          write(str(7),"(i8)") tagSend(iSend)
          write(str(8),"(i8)") m1
          call MsgDump(h//" dispatch send of global section (1:"//&
               trim(adjustl(str(8)))//","//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")"//&
               " with "//trim(adjustl(str(5)))//&
               " reals to MPI proc "//trim(adjustl(str(6)))//&
               " with tag "//trim(adjustl(str(7))))
       end if
       call parf_send_noblock_real(bufSend(i1), i2, procSend(iSend), &
            tagSend(iSend), reqSend(iSend))
    end do

    !     do cnt = 1, nRecv
    do iRecV = 1, nRecv
       call parf_wait_any_nostatus(nRecv,reqRecv,recNum)
       i1 = bufRecvStart(recNum)
       iCnt = bufRecvStart(recNum)
       i2 = bufRecvLength(recNum)
       do j = jaRecv(recNum), jzRecv(recNum)
          do i = iaRecv(recNum), izRecv(recNum)
             field(1:m1,i,j) = bufRecv(iCnt:iCnt+m1-1)
             iCnt = iCnt+m1
          end do
       end do
       if (dumpLocal) then
          write(str(1),"(i8)") iaRecv(recNum)+i0LGZ
          write(str(2),"(i8)") izRecv(recNum)+i0LGZ
          write(str(3),"(i8)") jaRecv(recNum)+j0LGZ
          write(str(4),"(i8)") jzRecv(recNum)+j0LGZ
          write(str(5),"(i8)") i2
          write(str(6),"(i8)") recNum
          write(str(7),"(i8)") tagRecv(recNum)
          write(str(8),"(i8)") m1
          write(str(9),"(i8)") procRecv(recNum)
          call MsgDump(h//" recv #"//trim(adjustl(str(6)))//&
               " from MPI proc "//trim(adjustl(str(9)))//&
               " of tag "//trim(adjustl(str(7)))//&
               " with "//trim(adjustl(str(5)))//" reals and stores at global section (1:"//&
               trim(adjustl(str(8)))//","//&
               trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
               trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//")")
       end if
    end do
    do iRecS=1,nSend
       call parf_wait_all_nostatus(nSend,reqSend)
    end do
  end subroutine UpdateBorders





  subroutine InitialFieldsUpdate(ngrids,m1,m2,m3,Nm2,Nm3,ng,mynum, &
       nRec, procRecv, tagRecv, iaRecv, izRecv, jaRecv, jzRecv, &
       bufRecvStart, bufRecvLength, bufRecvTotalLength, &
       nSnd, procSend, tagSend, iaSend, izSend, jaSend, jzSend, &
       bufSendStart, bufSendLength, bufSendTotalLength)

    integer, intent(in) :: ngrids
    integer, intent(in) :: m1
    integer, intent(in) :: m2,nm2
    integer, intent(in) :: m3,nm3,ng,mynum
    integer, intent(in) :: nRec
    integer, intent(in) :: procRecv(nRec)
    integer, intent(in) :: tagRecv(nRec)
    integer, intent(in) :: iaRecv(nRec)
    integer, intent(in) :: izRecv(nRec)
    integer, intent(in) :: jaRecv(nRec)
    integer, intent(in) :: jzRecv(nRec)
    integer, intent(in) :: bufRecvStart(nRec)
    integer, intent(in) :: bufRecvLength(nRec)
    integer, intent(in) :: bufRecvTotalLength
    integer, intent(in) :: nSnd
    integer, intent(in) :: procSend(nSnd)
    integer, intent(in) :: tagSend(nSnd)
    integer, intent(in) :: iaSend(nSnd)
    integer, intent(in) :: izSend(nSnd)
    integer, intent(in) :: jaSend(nSnd)
    integer, intent(in) :: jzSend(nSnd)
    integer, intent(in) :: bufSendStart(nSnd)
    integer, intent(in) :: bufSendLength(nSnd)
    integer, intent(in) :: bufSendTotalLength

    integer :: i,j,k
    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(InitialFieldsUpdate)**"
    character(len=8) :: str(10)

    integer, save :: iupdate_dxy=0

    if(bufSendTotalLength==0 .or. bufRecvTotalLength==0) return


    do i=1,nm2
       do j=1,nm3
          do k=1,1!m1
             advmnt_g(ng)%l_dxtW(k,i,j)=advmnt_g(ng)%dxtW(i,j)
             advmnt_g(ng)%l_dytW(k,i,j)=advmnt_g(ng)%dytW(i,j)
          end do
       end do
    end do

    if (dumpLocal) then
       call MsgDump(h//" update borders of u3d")
    end if

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%u3d, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of v3d")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%v3d, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of dd0_3d")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%dd0_3d, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of dd0_3du")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%dd0_3du, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of dd0_3dv")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%dd0_3dv, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of dd0_3dw")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%dd0_3dw, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of den0_3d")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%den0_3d, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of den1_3d")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%den1_3d, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of den2_3d")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%den2_3d, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of den3_3d")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%den3_3d, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of l_dxtW")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%l_dxtW, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    if (dumpLocal) then
       call MsgDump(h//" update borders of l_dytW")
    end if
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%l_dytW, &
         nRec, procRecv, tagRecv, &
         iaRecv, izRecv, jaRecv, jzRecv, &
         bufRecvStart, bufRecvLength, bufRecvTotalLength, &
         nSnd, procSend, tagSend, &
         iaSend, izSend, jaSend, jzSend, &
         bufSendStart, bufSendLength, bufSendTotalLength)

    do i=1,nm2
       do j=1,nm3
          advmnt_g(ng)%dxtW(i,j)=advmnt_g(ng)%l_dxtW(1,i,j)
          advmnt_g(ng)%dytW(i,j)=advmnt_g(ng)%l_dytW(1,i,j)
       end do
    end do
  end subroutine InitialFieldsUpdate





  subroutine Advec3d_X(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,den0,&
       den1,dt,dxx,&
       dd0,&
       qn,mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: den0(m1,m2,m3)
    real   , intent(in) :: den1(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dxx(m2,m3)
    real   , intent(in) :: dd0(m1,m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    integer :: ii
    integer :: ji
    integer :: ii0
    integer :: ji0
    integer :: ie
    integer :: je
    integer :: ie0
    integer :: je0
    integer :: ipos
    integer :: iia
    integer :: iiz
    integer :: nvar
    integer :: nf
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_X)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    ! copy input field to output field
    qn=q0
    imxmn=.false.

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do j=ja,jz
       do k=2,m1-1
          if(u(k,1,j)>=zr0) flux(k,1,j)= q0(k,1,j)*u(k,1,j)*dt*dd0(k,1,j)
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
    do j=ja,jz
       do  i=2,m2-1 ! ia,iz-1 or 1,iz-1
          do k=2,m1-1
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i-1,j),q0(k,i+1,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k,i-1,j),q0(k,i+1,j))+eps)        !       extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k,i,j  )< zr0) ck1= q0(k,i+1,j)
             if(u(k,i-1,j)>=zr0) ck2= q0(k,i-1,j)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    do j=ja,jz
       do  i=2,m2-1 ! ia,iz-1 or 1,iz-1
          do k=2,m1-1
             if(u(k,i,j)<zr0) cycle
             if(u(k,i-1,j)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
                x1n= (1.-x1)*(q0(k,i+1,j)-q0(k,i-1,j))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k,i-1,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i+1,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i+1,j))  ), max(q0(k,i,j),q0(k,i+1,j)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j) = max(vcmin(k,i,j),min(vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i-1,j)
             end if                                                  !Eq-9a
          end do
       end do
    end do

    ! If periodic boundary conditions are assumed, it is necessary
    !   to recalculate the updated mixing ratio at cell 1 if there
    !   is inflow to that cell from the boundary between IDIM and 1
    !   Here these statements are commented out, but should be uncommented
    !   if this subroutine is needed for periodic boundary conditions,
    !   and then one of the calling arguements to the subroutine is IPERIOD
    !   which is set to "1" if you assume period boundary conditions
    !      IF(IPERIOD==1) THEN
    !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
    !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
    !      END IF
    !
    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    !
    do j=ja,jz
       do k=2,m1-1
          if(u(k,m2-1,j)<zr0) flux(k,m2-1,j)= &
               q0(k,m2,j)*u(k,m2-1,j)*dt*dd0(k,m2-1,j)
       end do
    end do

    do j=ja,jz
       do i=m2-1,2,-1 !iz,ia,-1
          do k=2,m1-1
             if(u(k,i-1,j)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + &
                     flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k,i-1,j))/dxx(i,j)     ! Courant number
                x1n= (1.-x1)*(q0(k,i-1,j)-q0(k,i+1,j))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k,i+1,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i-1,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i-1,j)) ), max(q0(k,i,j),q0(k,i-1,j)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j)= max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i-1,j))/den1(k,i,j) ))
                flux(k,i-1,j)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop y-z
  end subroutine Advec3d_X





  subroutine Advec3d_Y(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,den0,&
       den1,dt,dxx,&
       dd0,&
       qn,mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: den0(m1,m2,m3)
    real   , intent(in) :: den1(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dxx(m2,m3)
    real   , intent(in) :: dd0(m1,m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    integer :: ii
    integer :: ji
    integer :: ii0
    integer :: ji0
    integer :: ie
    integer :: je
    integer :: ie0
    integer :: je0
    integer :: ipos
    integer :: iia
    integer :: iiz
    integer :: nvar
    integer :: nf
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_Y)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    ! copy input field to output field
    qn= q0
    imxmn=.false.

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do i=ia,iz
       do k=2,m1-1
          if(u(k,i,1)>=zr0) flux(k,i,1)= q0(k,i,1)*u(k,i,1)*dt*dd0(k,i,1)
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !	mixing ratio at t+dt. If these limits are ever violated,
    !	non-monotonic (oscillatory) behavior in solution results
    do i=ia,iz
       do  j=2,m3-1 ! ja,jz
          do k=2,m1-1
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i,j-1),q0(k,i,j+1))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k,i,j-1),q0(k,i,j+1))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k,i,j  )< zr0) ck1= q0(k,i,j+1)
             if(u(k,i,j-1)>=zr0) ck2= q0(k,i,j-1)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    ! Identify local max and min, specify mixing ratio limits at new time
    do i=ia,iz
       do  j=2,m3-1 ! ja,jz
          do k=2,m1-1
             if(u(k,i,j)<zr0) cycle
             if(u(k,i,j-1)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
                x1n= (1.-x1)*(q0(k,i,j+1)-q0(k,i,j-1))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k,i,j-1)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i,j+1)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i,j+1))  ), max(q0(k,i,j),q0(k,i,j+1)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j)= max(  vcmin(k,i,j),   min(   vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i,j-1) !Eq-9a
             end if
          end do
       end do
    end do

    ! If periodic boundary conditions are assumed, it is necessary
    !   to recalculate the updated mixing ratio at cell 1 if there
    !   is inflow to that cell from the boundary between IDIM and 1
    !   Here these statements are commented out, but should be uncommented
    !   if this subroutine is needed for periodic boundary conditions,
    !   and then one of the calling arguements to the subroutine is IPERIOD
    !   which is set to "1" if you assume period boundary conditions
    !      IF(IPERIOD==1) THEN
    !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
    !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
    !      END IF
    !
    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    !
    do i=ia,iz
       do k=2,m1-1
          if(u(k,i,m3-1)<zr0) flux(k,i,m3-1)= &
               q0(k,i,m3)*u(k,i,m3-1)*dt*dd0(k,i,m3-1)
       end do
    end do

    do i=ia,iz
       do j=m3-1,2,-1 !jz,ja,-1
          do k=2,m1-1
             if(u(k,i,j-1)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + &
                     flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k,i,j-1))/dxx(i,j)     ! Courant number
                x1n= (1.-x1)*(q0(k,i,j-1)-q0(k,i,j+1))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k,i,j+1)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k,i,j-1)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k,i,j-1)) ), max(q0(k,i,j),q0(k,i,j-1)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j)= max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i,j-1))/den1(k,i,j) ))
                flux(k,i,j-1)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop x-z
  end subroutine Advec3d_Y





  subroutine Advec3d_Z(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,den0,&
       den1,dt,dxx,&
       dd0,&
       qn,mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: den0(m1,m2,m3)
    real   , intent(in) :: den1(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dxx(m1)
    real   , intent(in) :: dd0(m1,m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_Z)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    ! copy input field to output field
    qn = q0
    imxmn=.false.


    ! Identify local max and min, specify mixing ratio limits at new time
    !  VCMAX and VCMIN are the absolute physical limits to the
    !     mixing ratio at t+dt. If these limits are ever violated,
    !     non-monotonic (oscillatory) behavior in solution results
   do j=ja,jz
       do i=ia,iz
          do k=2,m1-1 
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
             if(u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do


    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    do j=ja,jz
       do i=ia,iz
          if(u(1,i,j)>=zr0) flux(1,i,j)= &
               q0(1,i,j)*u(1,i,j)*dt*dd0(1,i,j)
          do k=2,m1-1
             if(u(k,i  ,j)<zr0) cycle
             if(u(k-1,i,j)<zr0) then
                flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
             else                              !      use upstream
                x1= dt*u(k,i,j)/dxx(k)               ! Courant number
                x1n= (1.-x1)*(q0(k+1,i,j)-q0(k-1,i,j))/4.

                ! First, estimate mixing ratio in outflowing fluid (Cf)
                cf= q0(k,i,j) + x1n                                       !Eq-4a

                !   Check to see if there is a peak (min) upwind and/or
                !    downwind of cell face
                if(imxmn(k-1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k+1,i,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
                !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"

                !   Limit Cf to be between mixing ratio on either side of edge
                !      where flux is being calculated
                cf1= min( max( cf, min(q0(k,i,j),q0(k+1,i,j))  ), max(q0(k,i,j),q0(k+1,i,j)) )

                !   Calculate mixing ratio at new time, but limit to physically
                !    reasonable values
                qn(k,i,j)= max(  vcmin(k,i,j),   min(   vcmax(k,i,j),          &   !eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))

                !   Re-calculate OUTFLOWING flux before moving on to next cell
                !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
                !    is encountered.
                flux(k,i,j)= dxx(k)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k-1,i,j)
             end if                                                  !Eq-9a
          end do
       end do
    end do

    ! Update mixing ratios and limit Fluxes going DOWN where u<0
    !  The logic of this loop through the grid line is identical
    !  to the "DO 10" Loop above, only you start at the highest I
    !  edge and work backwards to I=1
    do j=ja,jz
       do i=ia,iz
          if(u(m1-1,i,j)<zr0) flux(m1-1,i,j)=&
               q0(m1,i,j)*u(m1-1,i,j)*dt*dd0(m1-1,i,j)
          do k=m1-1,2,-1
             if(u(k-1,i,j)>=zr0) then           ! Inflow-only cell
                if(u(k,i,j)<zr0) qn(k,i,j)=  max(  vcmin(k,i,j),   min(   vcmax(k,i,j),&
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k) + flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))
             else
                x1=  dt*abs(u(k-1,i,j))/dxx(k)     ! Courant number
                x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
                cf= q0(k,i,j) + x1n                                       !Eq-4b
                if(imxmn(k+1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
                if(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
                cf1= min( max( cf, min(q0(k,i,j),q0(k-1,i,j)) ), max(q0(k,i,j),q0(k-1,i,j)) )
                if(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
                qn(k,i,j) = max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k)-x1*cf1*dd0(k-1,i,j))/den1(k,i,j) ))
                flux(k-1,i,j)=dxx(k)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
             end if
          end do
       end do
    end do !- big loop y-x
  end subroutine Advec3d_Z





  subroutine Advec3d_Z_sedim(m1,m2,m3,ia,iz,ja,jz,&
       q0,&
       u,&
       dt,&
       dzt,rtgt,&
       qn,&
       mynum)
    !-------------------------
    ! This subroutine calculates change in mixing ratio (Q0) during time
    !  step DT due to advection along a grid IDIM in length. Mixing ratios
    !  from host code (C) are loaded into Q0 array, which is updated to QN.
    !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
    !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
    !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
    !
    ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
    ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
    ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
    ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
    ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
    !                 Q0 defined along 0 - IDIM+1 cells:
    !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
    !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
    !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
    !   lower BC |             <---   Q0 grid   --->         | upper BC
    !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
    !
    !  Input to this subroutine, provided in common /sub/, and the calling
    !  arguments to this subroutine:
    !     IDIM - #of grid cells being updated
    !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
    !                 additional boundary value mixing ratios padded into the
    !                 0th and IDIM+1 cell locations
    !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
    !                each grid cell in the array, units consistent with DX, DT
    !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in Calling code
    !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
    !                 multi-dimensional calculations, as noted in calling code
    !     DT-         time step- units consistent with U
    !     DXX(IDIM)-  Grid cell length along advection direction, Units
    !                   consistent with DT and U
    !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
    !                  (remains constant for all dimensions at the initial
    !                  fluid density of the 1st dimension of a 2-3 D calculation
    !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
    !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
    !               fluid density at the beginning of the 1st dimensional
    !               advection step of a 2 or 3 D advection calculation done one
    !               step at a time
    !
    !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
    !

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: q0(m1,m2,m3)
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dzt(m1)
    real   , intent(in) :: rtgt(m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    real :: flux(m1,m2,m3)
    real :: vcmax(m1,m2,m3)
    real :: vcmin(m1,m2,m3)
    logical :: imxmn(m1,m2,m3)
    real, parameter :: zr0=0.0
    real, parameter :: EPS=1.e-6
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_Z_sedim)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    qn = q0
    imxmn=.false.

    ! Identify local max and min, specify mixing ratio limits at new time
          !  VCMAX and VCMIN are the absolute physical limits to the
          !     mixing ratio at t+dt. If these limits are ever violated,
          !     non-monotonic (oscillatory) behavior in solution results
    do j=ja,jz
       do i=ia,iz
          do  k=2,m1-1 !
             imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
             ck1= q0(k,i,j)
             ck2= q0(k,i,j)
             if(-u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
             if(-u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
             if(k==2) ck2= q0(k,i,j) !for sedim only
             vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
             vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
          end do
       end do
    end do

    do j=ja,jz
       do i=ia,iz
          rtgti=1./rtgt(i,j)
          flux(m1-1,i,j)=q0(m1,i,j)*(-u(m1-1,i,j))*dt
          do k=m1-1,2,-1
             !srf       x1=  dt*ABS(u(k-1,i,j))/dxx(k)     ! Courant number
             x1=  dt*abs(u(k-1,i,j))*dzt(k)*rtgti     ! Courant number
             if(k==2) x1 = 0. ! no flux below sfc terrain,for sedim only
             x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
             cf= q0(k,i,j) + x1n                                       !Eq-4b
             if(imxmn(k+1,i,j)) cf= q0(k,i,j) +max(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
             if(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
             cf1= min( max( cf, min(q0(k,i,j),q0(k-1,i,j)) ), max(q0(k,i,j),q0(k-1,i,j)) )
             if(k>2) then  !for sedim only
                qn(k,i,j) = max(  vcmin(k,i,j),  min(   vcmax(k,i,j), 	  &   !Eq-3&8
                                !srf                 (q0(k,i,j)-flux(k,i,j)/dxx(k)      -x1*cf1) ))
                     (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1) ))
             else
                qn(k,i,j) = (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1)
             end if
             !srf	   flux(k-1,i,j)=dxx(k)             *(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
             flux(k-1,i,j)=(1./(dzt(k)*rtgti))*(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
          end do
       end do
    end do !- big loop y-x
  end subroutine Advec3d_Z_sedim





  subroutine Advec3d_Z_sedim_upw(m1,m2,m3, ia,iz,ja,jz,u,dt,dzt,rtgt,qn,mynum)

    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum
    real   , intent(in) :: u(m1,m2,m3)
    real   , intent(in) :: dt
    real   , intent(in) :: dzt(m1)
    real   , intent(in) :: rtgt(m2,m3)
    real   , intent(out):: qn(m1,m2,m3)

    integer :: i
    integer :: j
    integer :: k
    real :: cf
    real :: cf1
    real :: ck1
    real :: ck2
    real :: x1
    real :: x1n
    real :: rtgti

    logical, parameter :: dumpLocal=.false.
    character(len=*), parameter :: h="**(Advec3d_Z_sedim_upw)**"
    character(len=8) :: str(10)

    if (dumpLocal) then
       write(str(1),"(i8)") m1
       write(str(2),"(i8)") m2
       write(str(3),"(i8)") m3
       write(str(4),"(i8)") ia
       write(str(5),"(i8)") iz
       write(str(6),"(i8)") ja
       write(str(7),"(i8)") jz
       call MsgDump(h//" at surface area ("//&
            trim(adjustl(str(4)))//":"//trim(adjustl(str(5)))//","//&
            trim(adjustl(str(6)))//":"//trim(adjustl(str(7)))//")"//&
            " of fields dimensioned ("//&
            trim(adjustl(str(1)))//","//&
            trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//")")
    end if

    !- big loop y-x
    do j=ja,jz
       do i=ia,iz
          rtgti=1./rtgt(i,j)
          !srf dxx = dz = rtgti/dzt
          !srf qn(m1-1,i,j) = qn(m1-1,i,j) / (1.0 - dt*u(m1-1,i,j)/dxx(m1-1)      )
          qn(m1-1,i,j) = qn(m1-1,i,j) / (1.0 + dt*u(m1-1,i,j)*dzt(m1-1)*rtgti)
          do k=m1-2,2,-1 !
             !srf    qn(k,i,j)= 1.0/(1.0+dt*u(k,i,j)/dxx(k))&
             !srf               *( qn(k,i,j)+ dt*u(k,i,j) /dxx(k+1) * qn(k+1,i,j) )
             qn(k,i,j)= 1.0/(1.0 + dt*u(k,i,j)*dzt(k)*rtgti)&
                  *( qn(k,i,j) + dt*u(k+1,i,j)*dzt(k+1)*rtgti * qn(k+1,i,j) )
             !   tc(i,j,l,k) = 1.0/(1.0+dt_settl(k)*vd_cor/delz(i,j,l2))&
             !  	 *(tc(i,j,l,k) + dt_settl(k)*vd_cor /delz(i,j,l2-1) &
             !  	 * tc(i,j,l+1,k))
          end do
       end do
    end do !- big loop y-x
  end subroutine Advec3d_Z_sedim_upw





  subroutine CheckBorders(m1, m2, m3, field, &
       nRecv, procRecv, tagRecv, iaRecv, izRecv, jaRecv, jzRecv, &
       bufRecvStart, bufRecvLength, bufRecvTotalLength, &
       nSend, procSend, tagSend, iaSend, izSend, jaSend, jzSend, &
       bufSendStart, bufSendLength, bufSendTotalLength,mynum,op,ie)
    integer, intent(in) :: m1,mynum,op,ie
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    real,    intent(inout) :: field(m1,m2,m3)
    integer, intent(in) :: nRecv
    integer, intent(in) :: procRecv(nRecv)
    integer, intent(in) :: tagRecv(nRecv)
    integer, intent(in) :: iaRecv(nRecv)
    integer, intent(in) :: izRecv(nRecv)
    integer, intent(in) :: jaRecv(nRecv)
    integer, intent(in) :: jzRecv(nRecv)
    integer, intent(in) :: bufRecvStart(nRecv)
    integer, intent(in) :: bufRecvLength(nRecv)
    integer, intent(in) :: bufRecvTotalLength
    integer, intent(in) :: nSend
    integer, intent(in) :: procSend(nSend)
    integer, intent(in) :: tagSend(nSend)
    integer, intent(in) :: iaSend(nSend)
    integer, intent(in) :: izSend(nSend)
    integer, intent(in) :: jaSend(nSend)
    integer, intent(in) :: jzSend(nSend)
    integer, intent(in) :: bufSendStart(nSend)
    integer, intent(in) :: bufSendLength(nSend)
    integer, intent(in) :: bufSendTotalLength

    integer :: fout,i
    character :: opc

    fout=80+mynum
    opc='Y'
    if(op==1) opc='X'

    if(ie==0) then
       write (fout,'("Borders updated, direction ",A)') opc
       return
    end if


    write (fout,'(" Updating borders, direction  ",A)') opc ; call flush(fout)
    write (fout,'("nRecv: ",I3.3," nSend: ",I3.3)') nRecv,nSend; call flush(fout)
    write (fout,'("TotRecv: ",I6.6," TotSend: ",I6.6)') bufRecvTotalLength,bufSendTotalLength; call flush(fout)
    write (fout,'(A)') '---------------------------------- Send ---------------------------'; call flush(fout)
    write (fout,'(7(A3,1X),2(A,1X))') 'nSn','prc','tag','ia ','iz ','ja ','jz ','Start','Length'; call flush(fout)
    do i=1,nRecv
       write (fout,'(7(I3.3,1X),2(I6.6,1X))') i,procRecv(i),tagRecv(i),iaRecv(i),izRecv(i),jaRecv(i),jzRecv(i), &
            bufRecvStart(i),bufRecvLength(i); call flush(fout)
    end do
    write (fout,'(A)') '------------------------------- Receive  --------------------------'; call flush(fout)
    write (fout,'(7(A3,1X),2(A,1X))') 'nRv','prc','tag','ia ','iz ','ja ','jz ','Start','Length'
    do i=1,nRecv
       write (fout,'(7(I3.3,1X),2(I6.6,1X))')	i,procSend(i),tagSend(i),iaSend(i),izSend(i),jaSend(i),jzSend(i), &
            bufSendStart(i),bufSendLength(i); call flush(fout)
    end do
  end subroutine CheckBorders





  subroutine StoreNamelistFileAtRadvc_mnt(oneNamelistFile)

    ! import NameListFile values into module variables

    type(namelistFile), pointer :: oneNamelistFile

    advmnt = oneNamelistFile%advmnt
    GhostZoneLength=oneNamelistFile%GhostZoneLength
  end subroutine StoreNamelistFileAtRadvc_mnt
end module monotonic_adv
