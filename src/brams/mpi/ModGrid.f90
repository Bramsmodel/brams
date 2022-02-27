module ModGrid

  ! ModGrid: 

!!$  use ParLib, only: &
!!$       parf_barrier, &
!!$       parf_exit_mpi

  use ModNamelistFile, only: &
       NamelistFile

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModGridDims, only: &
       GridDims, &
       CreateGridDims, &
       DumpGridDims, &
       DestroyGridDims

  use ModDomainDecomp, only: &
       DomainDecomp, &
       CreateGlobalOwn, &
       CreateGlobalOwnWithBC, &
       CreateGlobalWithGhost, &
       CreateLocalOwn, &
       DumpDomainDecomp, &
       DestroyDomainDecomp

  use ModNeighbourNodes, only: &
       NeighbourNodes, &
       CreateNeighbourNodes, &
       DumpNeighbourNodes, &
       DestroyNeighbourNodes

  use ModMessageSet, only: &
       MessageSet, &
       DumpMessageSet, &
       CreateAcousticMessageSet, &
       DestroyAcousticMessageSet, &
       CreateDn0MessageSet, &
       DestroyDn0MessageSet, &
       CreateG3DMessageSet, &
       DestroyG3DMessageSet, &
       CreateSelectedGhostZoneMessageSet, &
       DestroySelectedGhostZoneMessageSet, &
       CreateAllGhostZoneMessageSet, &
       DestroyAllGhostZoneMessageSet, &
       CreateAcouDampOneMessageSet, &
       DestroyAcouDampOneMessageSet, &
       CreateWideGhostZoneMessageSet, &
       DestroyWideGhostZoneMessageSet

  use ModMonotonicAdvection, only: &
       MonotonicAdvection, &
       CreateMonotonicAdvection, &
       DestroyMonotonicAdvection


  ! JP: temporariamente usa variaveis globais enquanto
  !     var_tables nao for inclusa no tipo Grid

  use var_tables, only: &
       num_var, &
       vtab_r

  use meteogramType, only: &
       PolygonContainer

  implicit none

  private
  public :: Grid
  public :: CreateGrid
  public :: InsertMessageSetAtOneGrid
  public :: DestroyGrid
  public :: DumpGrid


  type Grid
     integer :: Id
     ! Id: grid number on Namelist
     type(NamelistFile), pointer :: Ramsin => null()
     ! Ramsin: this grid namelist file
     type(ParallelEnvironment), pointer :: ParEnv => null()
     ! ParEnv: mpi size, rank and communicator for this run
     type(GridDims), pointer :: GridSize => null()
     ! GridSize: this grid dimensions as defined by namelist
     type(DomainDecomp), pointer :: GlobalOwn => null()
     ! GlobalOwn: global indices of this grid domain
     !            decomposition (domain partition) owned by
     !            each rank - Ghost Zone not included
     type(DomainDecomp), pointer :: GlobalOwnWithBC => null()
     ! GlobalOwnWithBC: global indices of this grid domain
     !            decomposition (domain partition) owned by
     !            each rank including Boundary Conditions
     type(DomainDecomp), pointer :: GlobalWithGhost => null()
     ! GlobalWithGhost: global indices of this grid domain
     !                  decomposition at each rank, including
     !                  the owned points and a ghost zone of length one.
     !                  Not a domain partition, due to ghost
     !                  zone inclusion.
     type(DomainDecomp), pointer :: LocalOwn => null()
     ! LocalOwn: local indices of this grid domain
     !           decomposition owned by each rank. 
     !           Convertion of GlobalOwn to local indices
     type(NeighbourNodes), pointer :: Neigh => null()
     ! Neigh: list of BRAMS process numbers that are neighbours
     !        of this node for usual ghost zone update operations
     type(MessageSet), pointer :: AcouSendU => null()
     type(MessageSet), pointer :: AcouRecvU => null()
     ! AcouSend/RecvU: Ghost Zone update at acoust_new and acoust_adap
     type(MessageSet), pointer :: AcouSendV => null()
     type(MessageSet), pointer :: AcouRecvV => null()
     ! AcouSend/RecvV: Ghost Zone update at acoust_new and acoust_adap
     type(MessageSet), pointer :: AcouSendPNorth => null()
     type(MessageSet), pointer :: AcouRecvPNorth => null()
     ! AcouSend/RecvPNorth: Ghost Zone update at acoust_new and acoust_adap
     type(MessageSet), pointer :: AcouSendPEast => null()
     type(MessageSet), pointer :: AcouRecvPEast => null()
     ! AcouSend/RecvPEast: Ghost Zone update at acoust_new and acoust_adap
     type(MessageSet), pointer :: AcouSendUV => null()
     type(MessageSet), pointer :: AcouRecvUV => null()
     ! AcouSend/RecvUV: Ghost Zone update at acoust_new and acoust_adap
     type(MessageSet), pointer :: AcouSendWP => null()
     type(MessageSet), pointer :: AcouRecvWP => null()
     ! AcouSend/RecvWP: Ghost Zone update at acoust_new and acoust_adap
     type(MessageSet), pointer :: SendDn0u => null()
     type(MessageSet), pointer :: RecvDn0u => null()
     type(MessageSet), pointer :: SendDn0v => null()
     type(MessageSet), pointer :: RecvDn0v => null()
     ! Send/RecvDn0u/v: Ghost Zone update at FillDn0uv
     type(MessageSet), pointer :: SendG3D => null()
     type(MessageSet), pointer :: RecvG3D => null()
     ! Send/RecvG3D: Ghost Zone update at cuparm_grell3_catt
     type(MessageSet), pointer :: SelectedGhostZoneSend => null()
     type(MessageSet), pointer :: SelectedGhostZoneRecv => null()
     ! SelectedGhostZoneSend/RecvG3D: Ghost Zone update at timestep and timestep_rk
     type(MessageSet), pointer :: AllGhostZoneSend => null()
     type(MessageSet), pointer :: AllGhostZoneRecv => null()
     ! AllGhostZoneSend/RecvG3D: Ghost Zone update at PostProcess
     ! type(MessageSet) contains all information required for
     ! ghost zone update. See description at ModMessageSet 
     type(PolygonContainer), pointer :: meteoPolygons => null()

     type(MessageSet), pointer :: AcouDampDivSend => null()
     type(MessageSet), pointer :: AcouDampDivRecv => null()
     type(MessageSet), pointer :: AcouDampPPSend => null()
     type(MessageSet), pointer :: AcouDampPPRecv => null()
     type(MessageSet), pointer :: AcouDampAlphaSend => null()
     type(MessageSet), pointer :: AcouDampAlphaRecv => null()
     type(MessageSet), pointer :: AcouDampThtSend => null()
     type(MessageSet), pointer :: AcouDampThtRecv => null()
     ! AcouDampSend/Recv: Ghost Zone update of a single field
     ! on Runge Kutta Dynamics, acoust_new and init_div_damping_coef.
     ! Fields to update are local variables to these procedures,
     ! allocated and deallocated at each call. As so, field
     ! memory address vary with procedure invocation
     type(MessageSet), pointer :: WideGhostZoneSend => null()
     type(MessageSet), pointer :: WideGhostZoneRecv => null()
     ! WideGhostZoneSend/Recv: Ghost Zone update of four fields
     ! with large ghost zones on advect_ws.
     ! Fields to update are local variables to these procedures,
     ! allocated and deallocated at each call. As so, field
     ! memory address vary with procedure invocation
  end type Grid



contains



  ! CreateGrid: create and fill variable of this type



  subroutine CreateGrid(gridId, oneNamelistFile, oneParallelEnvironment, oneGrid)
    integer, intent(in) :: gridId
    type(NamelistFile), pointer :: oneNamelistFile
    type(ParallelEnvironment), pointer :: oneParallelEnvironment
    type(Grid), pointer :: oneGrid

    logical :: largeGhostZone
    ! if will use extended ghost zone width

    character(len=16) :: c0, c1
    character(len=*), parameter :: h="**(CreateGrid)**"
    logical, parameter :: dumpLocal=.true.

    ! correctness of input arguments

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneNamelistFile")
    else if (.not. associated(oneParallelEnvironment)) then
       call fatal_error(h//" invoked with null oneParallelEnvironment")
    else if (associated(oneGrid)) then
       call fatal_error(h//" invoked with already associated oneGrid")
    end if

    ! create a variable of type grid and fill entries

    allocate(oneGrid)

    ! stores input arguments

    oneGrid%id = gridId
    oneGrid%Ramsin => oneNamelistFile
    oneGrid%ParEnv => oneParallelEnvironment

    ! run options

    largeGhostZone = oneNamelistFile%dyncore_flag==2 .and. oneNamelistFile%advmnt==1
    largeGhostZone = largeGhostZone .or. (&
         (oneNamelistFile%dyncore_flag==0 .or. oneNamelistFile%dyncore_flag==1) .and. &
         oneNamelistFile%advmnt>=1)

    ! store GridDims extracted from OneNamelistFile 

    oneGrid%GridSize => CreateGridDims(gridId, &
         oneNamelistFile)

    ! compute domain decomposition, obtaining
    ! cells owned by each rank and store at GlobalOwn

    oneGrid%GlobalOwn => CreateGlobalOwn(&
         GridSize=oneGrid%GridSize, &
         ParEnv=oneGrid%ParEnv, &
         varName="GlobalOwn" &
         )


    oneGrid%GlobalOwnWithBC => CreateGlobalOwnWithBC(&
         GridSize=oneGrid%GridSize, &
         ParEnv=oneGrid%ParEnv, &
         GlobalOwn=oneGrid%GlobalOwn &
         )

    ! insert original ghost zone of widht 1
    ! at GlobalOwn and store at GlobalWithGhost

    oneGrid%GlobalWithGhost => CreateGlobalWithGhost(&
         GridSize=oneGrid%GridSize, &
         ParEnv=oneGrid%ParEnv, &
         GlobalOwn=oneGrid%GlobalOwn, &
         GhostZoneWidth=1, &
         varName="GlobalWithGhost" &
         )

    ! convert global indices from GlobalWithGhost
    ! into local indices stored at LocalOwn

    oneGrid%LocalOwn => CreateLocalOwn(&
         ParEnv=oneGrid%ParEnv, &
         GlobalWithGhost=oneGrid%GlobalWithGhost, &
         GlobalOwn=oneGrid%GlobalOwn, &
         varName="LocalOwn" &
         )

    ! neighbour nodes for original ghost zone update operations

    oneGrid%Neigh => CreateNeighbourNodes(&
         ParEnv=oneGrid%ParEnv, &
         GlobalOwn=oneGrid%GlobalOwn, &
         GlobalWithGhost=oneGrid%GlobalWithGhost, &
         varName="oneGrid%Neigh" &
         )

    if (dumpLocal) then
       call MsgDump(h//" dumping OneGrid at the end")
       call DumpGrid(OneGrid)
    end if
  end subroutine CreateGrid





  subroutine InsertMessageSetAtOneGrid(oneGrid)
    type(Grid), pointer :: oneGrid

    character(len=16) :: str(10)
    character(len=*), parameter :: h="**(InsertMessageSetAtOneGrid)**"
    logical, parameter :: dumpLocal=.false.

    integer, parameter :: TagU=25
    integer, parameter :: TagV=26
    integer, parameter :: TagPNorth=27
    integer, parameter :: TagPEast=28
    integer, parameter :: TagUV=29
    integer, parameter :: TagWP=30
    integer, parameter :: TagDn0u=31
    integer, parameter :: TagDn0v=32
    integer, parameter :: TagG3D=33
    integer, parameter :: TagSelectedGhostZone=34
    integer, parameter :: TagAllGhostZone=35
    integer, parameter :: TagAcouDampDiv=36
    integer, parameter :: TagAcouDampPP=37
    integer, parameter :: TagAcouDampAlpha=38
    integer, parameter :: TagAcouDampTht=39
    integer, parameter :: TagWideGhostZone=40


    ! Field pointer for fields not yet allocated
    ! not yet allocated; CreateAcouDampOneMessageSet
    ! takes bounds from field pointer.
    ! Field address will be replaced by
    ! a call to UpdateFieldAdress, once the
    ! desired field is allocated
    integer :: ierr
    integer :: myNum
    integer :: lbx, ubx
    integer :: lby, uby
    integer :: lbz, ubz
    real, pointer :: boundsRegularGhostZone(:,:,:)

    if (.not. associated(oneGrid)) then
       call fatal_error(h//" invoked with null grid")
    end if

    call CreateAcousticMessageSet(oneGrid%Id, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, &
         oneGrid%GlobalOwnWithBC, &
         oneGrid%GlobalWithGhost, &
         oneGrid%AcouSendU, oneGrid%AcouRecvU, TagU, &
         oneGrid%AcouSendV, oneGrid%AcouRecvV, TagV, &
         oneGrid%AcouSendPNorth, oneGrid%AcouRecvPNorth, TagPNorth, &
         oneGrid%AcouSendPEast, oneGrid%AcouRecvPEast, TagPEast, &
         oneGrid%AcouSendUV, oneGrid%AcouRecvUV, TagUV, &
         oneGrid%AcouSendWP, oneGrid%AcouRecvWP, TagWP)

    call CreateDn0MessageSet(oneGrid%Id, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
         oneGrid%SendDn0u, oneGrid%RecvDn0u, TagDn0u, &
         oneGrid%SendDn0v, oneGrid%RecvDn0v, TagDn0v)

    call CreateG3DMessageSet(oneGrid%Id, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhost, &
         oneGrid%Ramsin, &
         oneGrid%SendG3D, oneGrid%RecvG3D, TagG3D)

    ! temporariamente, num_var e vtab_r sao variaveis globais,
    ! enquanto nao inclusas no tipo Grid

    call CreateSelectedGhostZoneMessageSet(&
         oneGrid%Id, num_var, vtab_r, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhost, &
         oneGrid%SelectedGhostZoneSend, &
         oneGrid%SelectedGhostZoneRecv, &
         TagSelectedGhostZone)

    call CreateAllGhostZoneMessageSet(&
         oneGrid%Id, num_var, vtab_r, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhost, &
         oneGrid%AllGhostZoneSend, &
         oneGrid%AllGhostZoneRecv, &
         TagAllGhostZone)

    ! allocate field with desired bounds
    ! desired bonds assure that field section will
    ! be correctly computed for fields that were
    ! not yet allocated, such as local fields at procedures
    myNum=oneGrid%ParEnv%myNum
    lbx=1
    ubx=oneGrid%LocalOwn%nx(myNum)
    lby=1
    uby=oneGrid%LocalOwn%ny(myNum)
    lbz=1
    ubz=oneGrid%GridSize%nnzp
    allocate(boundsRegularGhostZone(lbz:ubz,lbx:ubx,lby:uby), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") lbz
       write(str(2),"(i8)") ubz
       write(str(3),"(i8)") lbx
       write(str(4),"(i8)") ubx
       write(str(5),"(i8)") lby
       write(str(6),"(i8)") uby
       write(str(7),"(i8)") ierr
       call fatal_error(h//" allocate boundsRegularGhostZone("//&
            trim(adjustl(str(1)))//":"//trim(adjustl(str(2)))//","//&
            trim(adjustl(str(3)))//":"//trim(adjustl(str(4)))//","//&
            trim(adjustl(str(5)))//":"//trim(adjustl(str(6)))//&
            ") fails with stat="//trim(adjustl(str(7))))
    end if

    ! use desired bounds fields to create AcouDamp Message Sets;
    ! correct field memory address by invoking UpdateFieldAdress
    ! prior to use the Message Sets

    call CreateAcouDampOneMessageSet(&
         boundsRegularGhostZone, "Div", 3,  &
         oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
         TagAcouDampDiv, "AcouDampDivSend", "AcouDampDivRecv", &
         oneGrid%AcouDampDivSend, oneGrid%AcouDampDivRecv)

    call CreateAcouDampOneMessageSet(&
         boundsRegularGhostZone, "PP", 3,  &
         oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
         TagAcouDampPP, "AcouDampPPSend", "AcouDampPPRecv", &
         oneGrid%AcouDampPPSend, oneGrid%AcouDampPPRecv)

    call CreateAcouDampOneMessageSet(&
         boundsRegularGhostZone, "Alpha", 3,  &
         oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
         TagAcouDampAlpha, "AcouDampAlphaSend", "AcouDampAlphaRecv", &
         oneGrid%AcouDampAlphaSend, oneGrid%AcouDampAlphaRecv)

    call CreateAcouDampOneMessageSet(&
         boundsRegularGhostZone, "Tht", 3,  &
         oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
         TagAcouDampTht, "AcouDampThtSend", "AcouDampThtRecv", &
         oneGrid%AcouDampThtSend, oneGrid%AcouDampThtRecv)

    deallocate(boundsRegularGhostZone, stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") ierr
       call fatal_error(h//" deallocate boundsRegularGhostZone"//&
            " fails with stat="//trim(adjustl(str(1))))
    end if

    call CreateWideGhostZoneMessageSet(&
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhost, oneGrid%LocalOwn, &
         TagWideGhostZone, oneGrid%WideGhostZoneSend, oneGrid%WideGhostZoneRecv)

    if (dumpLocal) then
       call MsgDump(h//" dumping oneGrid")
       call DumpGrid(OneGrid)
    end if
  end subroutine InsertMessageSetAtOneGrid



  ! DestroyGrid: deallocate area of a variable of type grid



  subroutine DestroyGrid(oneGrid)
    type(Grid), pointer :: oneGrid

    character(len=*), parameter :: h="**(DestroyGrid)**"

    if (associated(oneGrid)) then
       call DestroyGridDims(oneGrid%GridSize)
       call DestroyDomainDecomp(oneGrid%GlobalOwn)
       call DestroyDomainDecomp(oneGrid%GlobalWithGhost)
       call DestroyDomainDecomp(oneGrid%LocalOwn)
       call DestroyNeighbourNodes(oneGrid%Neigh)
       call DestroyAcousticMessageSet(&
            oneGrid%AcouSendU, oneGrid%AcouRecvU, &
            oneGrid%AcouSendV, oneGrid%AcouRecvV, &
            oneGrid%AcouSendPNorth, oneGrid%AcouRecvPNorth, &
            oneGrid%AcouSendPEast, oneGrid%AcouRecvPEast, &
            oneGrid%AcouSendUV, oneGrid%AcouRecvUV, &
            oneGrid%AcouSendWP, oneGrid%AcouRecvWP)
       call DestroyDn0MessageSet( &
            oneGrid%SendDn0u, oneGrid%RecvDn0u, &
            oneGrid%SendDn0v, oneGrid%RecvDn0v)
       call DestroyG3DMessageSet( &
            oneGrid%SendG3D, oneGrid%RecvG3D)
       call DestroySelectedGhostZoneMessageSet( &
            oneGrid%SelectedGhostZoneSend, &
            oneGrid%SelectedGhostZoneRecv)
       call DestroyAllGhostZoneMessageSet( &
            oneGrid%AllGhostZoneSend, &
            oneGrid%AllGhostZoneRecv)
       call DestroyAcouDampOneMessageSet( &
            oneGrid%AcouDampDivSend, oneGrid%AcouDampDivRecv)
       call DestroyAcouDampOneMessageSet( &
            oneGrid%AcouDampPPSend, oneGrid%AcouDampPPRecv)
       call DestroyAcouDampOneMessageSet( &
            oneGrid%AcouDampAlphaSend, oneGrid%AcouDampAlphaRecv)
       call DestroyAcouDampOneMessageSet( &
            oneGrid%AcouDampThtSend, oneGrid%AcouDampThtRecv)
       call DestroyWideGhostZoneMessageSet(&
            oneGrid%WideGhostZoneSend, oneGrid%WideGhostZoneRecv)
       deallocate(oneGrid)
    end if
    nullify(oneGrid)
  end subroutine DestroyGrid



  ! DumpGrid:


  subroutine DumpGrid(oneGrid)
    type(Grid), pointer :: oneGrid

    character(len=8) :: c0
    character(len=*), parameter :: h="**(DumpGrid)**"

    if (.not. associated(oneGrid)) then
       call fatal_error(h//" invoked with null oneGrid")
    else if (.not. associated(oneGrid%Ramsin)) then
       call fatal_error(h//" invoked with null oneGrid%Ramsin")
    else if (.not. associated(oneGrid%ParEnv)) then
       call fatal_error(h//" invoked with null oneGrid%ParEnv")
    end if

    write(c0,"(i8)") oneGrid%Id
    call MsgDump(h//" for grid "//trim(adjustl(c0)))

    call MsgDump(h//" dumping component GridSize")
    call DumpGridDims(oneGrid%GridSize)

    call MsgDump(h//" dumping domain decomposed components")
    call DumpDomainDecomp(oneGrid%GlobalOwn, "GlobalOwn")
    call DumpDomainDecomp(oneGrid%GlobalOwnWithBC, "GlobalOwnWithBC")
    call DumpDomainDecomp(oneGrid%GlobalWithGhost, "GlobalWithGhost")
    call DumpDomainDecomp(oneGrid%LocalOwn, "LocalOwn")

    call MsgDump(h//" dumping neighborhood components")
    call DumpNeighbourNodes(oneGrid%Neigh,"oneGrid%Neigh")

    call MsgDump(h//" dumping message set components")
    call MsgDump(h//" dumping AcouSendU")
    call DumpMessageSet(oneGrid%AcouSendU)
    call MsgDump(h//" dumping AcouRecvU")
    call DumpMessageSet(oneGrid%AcouRecvU)
    call MsgDump(h//" dumping AcouSendV")
    call DumpMessageSet(oneGrid%AcouSendV)
    call MsgDump(h//" dumping AcouRecvV")
    call DumpMessageSet(oneGrid%AcouRecvV)
    call MsgDump(h//" dumping AcouSendPNorth")
    call DumpMessageSet(oneGrid%AcouSendPNorth)
    call MsgDump(h//" dumping AcouRecvPNorth")
    call DumpMessageSet(oneGrid%AcouRecvPNorth)
    call MsgDump(h//" dumping AcouSendPEast")
    call DumpMessageSet(oneGrid%AcouSendPEast)
    call MsgDump(h//" dumping AcouRecvPEast")
    call DumpMessageSet(oneGrid%AcouRecvPEast)
    call MsgDump(h//" dumping AcouSendUV")
    call DumpMessageSet(oneGrid%AcouSendUV)
    call MsgDump(h//" dumping AcouRecvUV")
    call DumpMessageSet(oneGrid%AcouRecvUV)
    call MsgDump(h//" dumping AcouSendWP")
    call DumpMessageSet(oneGrid%AcouSendWP)
    call MsgDump(h//" dumping AcouRecvWP")
    call DumpMessageSet(oneGrid%AcouRecvWP)
    call MsgDump(h//" dumping SendDn0u")
    call DumpMessageSet(oneGrid%SendDn0u)
    call MsgDump(h//" dumping RecvDn0u")
    call DumpMessageSet(oneGrid%RecvDn0u)
    call MsgDump(h//" dumping SendDn0v")
    call DumpMessageSet(oneGrid%SendDn0v)
    call MsgDump(h//" dumping RecvDn0v")
    call DumpMessageSet(oneGrid%RecvDn0v)
    call MsgDump(h//" dumping SendG3D")
    call DumpMessageSet(oneGrid%SendG3D)
    call MsgDump(h//" dumping RecvG3D")
    call DumpMessageSet(oneGrid%RecvG3D)
    call MsgDump(h//" dumping SelectedGhostZoneSend")
    call DumpMessageSet(oneGrid%SelectedGhostZoneSend)
    call MsgDump(h//" dumping SelectedGhostZoneRecv")
    call DumpMessageSet(oneGrid%SelectedGhostZoneRecv)
    call MsgDump(h//" dumping AllGhostZoneSend")
    call DumpMessageSet(oneGrid%AllGhostZoneSend)
    call MsgDump(h//" dumping AllGhostZoneRecv")
    call DumpMessageSet(oneGrid%AllGhostZoneRecv)
    call MsgDump(h//" dumping AcouDampDivSend")
    call DumpMessageSet(oneGrid%AcouDampDivSend)
    call MsgDump(h//" dumping AcouDampDivRecv")
    call DumpMessageSet(oneGrid%AcouDampDivRecv)
    call MsgDump(h//" dumping AcouDampPPSend")
    call DumpMessageSet(oneGrid%AcouDampPPSend)
    call MsgDump(h//" dumping AcouDampPPRecv")
    call DumpMessageSet(oneGrid%AcouDampPPRecv)
    call MsgDump(h//" dumping AcouDampAlphaSend")
    call DumpMessageSet(oneGrid%AcouDampAlphaSend)
    call MsgDump(h//" dumping AcouDampAlphaRecv")
    call DumpMessageSet(oneGrid%AcouDampAlphaRecv)
    call MsgDump(h//" dumping AcouDampThtSend")
    call DumpMessageSet(oneGrid%AcouDampThtSend)
    call MsgDump(h//" dumping AcouDampThtRecv")
    call DumpMessageSet(oneGrid%AcouDampThtRecv)

    call MsgDump(h//" dumping WideGhostZoneSend")
    call DumpMessageSet(oneGrid%WideGhostZoneSend)
    call MsgDump(h//" dumping WideGhostZoneRecv")
    call DumpMessageSet(oneGrid%WideGhostZoneRecv)
  end subroutine DumpGrid
end module ModGrid
