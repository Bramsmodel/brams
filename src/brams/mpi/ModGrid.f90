module ModGrid

  ! ModGrid: 

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
       CreateAcoustNewMessageSet, &
       DestroyAcoustNewMessageSet, &
       CreateWideGhostZoneMessageSet, &
       DestroyWideGhostZoneMessageSet, &
       CreateAdvMntMessageSet, &
       DestroyAdvMntMessageSet

  use ModNodeDimensions, only: &
       NodeDimensions, &
       CreateNodeDimensions, &
       DestroyNodeDimensions, &
       DumpNodeDimensions

  ! JP: temporariamente usa variaveis globais enquanto
  !     var_tables nao for inclusa no tipo Grid

  use var_tables, only: &
       num_var, &
       vtab_r, &
       num_scalar, &
       scalar_table, &
       scalar_tab

  use mem_tend, only: &
       tend

  use meteogramType, only: &
       PolygonContainer

  implicit none

  private
  public :: Grid
  public :: CreateGrid
  public :: InsertMessageSetAtOneGrid
  public :: DestroyGrid
  public :: DumpGrid
  public :: InsertScalarTabAtOneGrid
  public :: DeepCopyToScalarTabAtOneGrid
  public :: DeepCopyFromScalarTabAtOneGrid

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
     type(DomainDecomp), pointer :: GlobalWithGhostAdvMnt => null()
     ! GlobalWithGhostAdvMnt: global indices of this grid domain
     !                        decomposition at each rank, including
     !                        the owned points and a ghost zone of 
     !                        parametrized width, used at MonotonicAdvection. 
     !                        Not a domain partition, due to ghost zone inclusion.
     type(DomainDecomp), pointer :: LocalOwnAdvMnt => null()
     ! LocalOwnAdvMnt: local indices of this grid domain
     !                 decomposition owned by each rank,
     !                 use at MonotonicAdvection.
     !                 Convertion of GlobalWithGhostAdvMnt to local indices
     type(NodeDimensions), pointer :: NodeDims => null()
     ! NodeDims: indices and dimensions of this process
     ! domain decomposed sub-domain
     type(NodeDimensions), pointer :: NodeDimsAdvMnt => null()
     ! NodeDimsAdvMnt: indices and dimensions of this process
     ! domain decomposed sub-domain for use inside MonotonicAdvection
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

     type(MessageSet), pointer :: AcoustNewDivSend => null()
     type(MessageSet), pointer :: AcoustNewDivRecv => null()
     type(MessageSet), pointer :: AcoustNewPPSend => null()
     type(MessageSet), pointer :: AcoustNewPPRecv => null()
     type(MessageSet), pointer :: AcoustNewAlphaSend => null()
     type(MessageSet), pointer :: AcoustNewAlphaRecv => null()
     type(MessageSet), pointer :: AcoustNewThtSend => null()
     type(MessageSet), pointer :: AcoustNewThtRecv => null()
     ! AcoustNewSend/Recv: Ghost Zone update of a single field
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

     type(MessageSet), pointer :: AdvMntUVSendX => null()
     type(MessageSet), pointer :: AdvMntUVRecvX => null()
     type(MessageSet), pointer :: AdvMntUVSendY => null()
     type(MessageSet), pointer :: AdvMntUVRecvY => null()

     type(MessageSet), pointer :: AdvMntDxDySendX => null()
     type(MessageSet), pointer :: AdvMntDxDyRecvX => null()
     type(MessageSet), pointer :: AdvMntDxDySendY => null()
     type(MessageSet), pointer :: AdvMntDxDyRecvY => null()

     type(MessageSet), pointer :: AdvMntDd0SendX => null()
     type(MessageSet), pointer :: AdvMntDd0RecvX => null()
     type(MessageSet), pointer :: AdvMntDd0SendY => null()
     type(MessageSet), pointer :: AdvMntDd0RecvY => null()

     type(MessageSet), pointer :: AdvMntDenSendX => null()
     type(MessageSet), pointer :: AdvMntDenRecvX => null()
     type(MessageSet), pointer :: AdvMntDenSendY => null()
     type(MessageSet), pointer :: AdvMntDenRecvY => null()

     type(MessageSet), pointer :: AdvMntScaSendX => null()
     type(MessageSet), pointer :: AdvMntScaRecvX => null()
     type(MessageSet), pointer :: AdvMntScaSendY => null()
     type(MessageSet), pointer :: AdvMntScaRecvY => null()

     type(scalar_table), pointer :: ScalarTab(:) => null()
     integer :: ScalarTabSize=0
  end type Grid



contains



  ! CreateGrid: create and fill variable of this type



  subroutine CreateGrid(gridId, oneNamelistFile, oneParallelEnvironment, oneGrid)
    integer, intent(in) :: gridId
    type(NamelistFile), pointer :: oneNamelistFile
    type(ParallelEnvironment), pointer :: oneParallelEnvironment
    type(Grid), pointer :: oneGrid

    character(len=16) :: c0, c1
    character(len=*), parameter :: h="**(CreateGrid)**"
    logical, parameter :: dumpLocal=.false.

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

    ! include boundary conditions (no ghost zone)

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

    ! this node dimensions and indexing limits

    oneGrid%NodeDims => CreateNodeDimensions(&
         GridSize=oneGrid%GridSize, &
         ParEnv=oneGrid%ParEnv, &
         LocalOwn=oneGrid%LocalOwn, &
         GlobalOwn=oneGrid%GlobalOwn, &
         verticalGhostZoneWidth=0, &
         surfaceGhostZoneWidth=1, &
         varName="NodeDims" &
         )

    ! for MonotonicAdvection, insert ghost zone of parametrized widht
    ! at GlobalOwn and store at GlobalWithGhostAdvMnt

    oneGrid%GlobalWithGhostAdvMnt => CreateGlobalWithGhost(&
         GridSize=oneGrid%GridSize, &
         ParEnv=oneGrid%ParEnv, &
         GlobalOwn=oneGrid%GlobalOwn, &
         GhostZoneWidth=oneNamelistFile%ghostzonelength, &
         varName="GlobalWithGhostAdvMnt" &
         )

    ! convert global indices from GlobalWithGhostAdvMnt
    ! into local indices stored at LocalOwnAdvMnt

    oneGrid%LocalOwnAdvMnt => CreateLocalOwn(&
         ParEnv=oneGrid%ParEnv, &
         GlobalWithGhost=oneGrid%GlobalWithGhostAdvMnt, &
         GlobalOwn=oneGrid%GlobalOwn, &
         varName="LocalOwnAdvMnt" &
         )

    ! this node dimensions and indexing limits

    oneGrid%NodeDimsAdvMnt => CreateNodeDimensions(&
         GridSize=oneGrid%GridSize, &
         ParEnv=oneGrid%ParEnv, &
         LocalOwn=oneGrid%LocalOwnAdvMnt, &
         GlobalOwn=oneGrid%GlobalOwn, &
         verticalGhostZoneWidth=0, &
         surfaceGhostZoneWidth=oneNamelistFile%ghostzonelength, &
         varName="NodeDimsAdvMnt" &
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
    integer, parameter :: TagAcoustNewDiv=36
    integer, parameter :: TagAcoustNewPP=37
    integer, parameter :: TagAcoustNewAlpha=38
    integer, parameter :: TagAcoustNewTht=39
    integer, parameter :: TagWideGhostZone=40
    integer, parameter :: TagAdvMntUVX=41
    integer, parameter :: TagAdvMntUVY=42
    integer, parameter :: TagAdvMntDxDyX=43
    integer, parameter :: TagAdvMntDxDyY=44
    integer, parameter :: TagAdvMntDd0X=45
    integer, parameter :: TagAdvMntDd0Y=46
    integer, parameter :: TagAdvMntDenX=47
    integer, parameter :: TagAdvMntDenY=48
    integer, parameter :: TagAdvMntScaX=49
    integer, parameter :: TagAdvMntScaY=50

    ! Field pointer for fields not yet allocated
    ! not yet allocated; CreateAcoustNewMessageSet
    ! takes bounds from field pointer.
    ! Field address will be replaced at
    ! PostSendRecvMsgs, since when this procedure
    ! is invoked, field ought to be allocated
    integer :: ierr
    integer :: lbx, ubx
    integer :: lby, uby
    integer :: lbz, ubz

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

    ! use desired bounds fields to create AcoustNew Message Sets;
    ! correct field memory address by invoking UpdateFieldAdress
    ! prior to use the Message Sets

    call CreateAcoustNewMessageSet(&
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, oneGrid%NodeDims, &
         TagAcoustNewDiv, oneGrid%AcoustNewDivSend, oneGrid%AcoustNewDivRecv, &
         TagAcoustNewPP, oneGrid%AcoustNewPPSend, oneGrid%AcoustNewPPRecv, &
         TagAcoustNewAlpha, oneGrid%AcoustNewAlphaSend, oneGrid%AcoustNewAlphaRecv, &
         TagAcoustNewTht, oneGrid%AcoustNewThtSend, oneGrid%AcoustNewThtRecv, &
         tend%tht_rk)


    call CreateWideGhostZoneMessageSet(&
         oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhost, &
         oneGrid%NodeDims%mzp, 1, oneGrid%NodeDims%mzp, &
         TagWideGhostZone, oneGrid%WideGhostZoneSend, oneGrid%WideGhostZoneRecv)

    call CreateAdvMntMessageSet(&
         oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhostAdvMnt, &
         oneGrid%NodeDims, oneGrid%NodeDimsAdvMnt, &
         TagAdvMntUVX, oneGrid%AdvMntUVSendX, oneGrid%AdvMntUVRecvX, &
         TagAdvMntUVY, oneGrid%AdvMntUVSendY, oneGrid%AdvMntUVRecvY, &
         TagAdvMntDxDyX, oneGrid%AdvMntDxDySendX, oneGrid%AdvMntDxDyRecvX, &
         TagAdvMntDxDyY, oneGrid%AdvMntDxDySendY, oneGrid%AdvMntDxDyRecvY, &
         TagAdvMntDd0X, oneGrid%AdvMntDd0SendX, oneGrid%AdvMntDd0RecvX, &
         TagAdvMntDd0Y, oneGrid%AdvMntDd0SendY, oneGrid%AdvMntDd0RecvY, &
         TagAdvMntDenX, oneGrid%AdvMntDenSendX, oneGrid%AdvMntDenRecvX, &
         TagAdvMntDenY, oneGrid%AdvMntDenSendY, oneGrid%AdvMntDenRecvY, &
         TagAdvMntScaX, oneGrid%AdvMntScaSendX, oneGrid%AdvMntScaRecvX, &
         TagAdvMntScaY, oneGrid%AdvMntScaSendY, oneGrid%AdvMntScaRecvY)


    if (dumpLocal) then
       call MsgDump(h//" dumping oneGrid")
       call DumpGrid(OneGrid)
       call MsgDump(h//" done dumping oneGrid")
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
       call DestroyDomainDecomp(oneGrid%GlobalWithGhostAdvMnt)
       call DestroyDomainDecomp(oneGrid%LocalOwnAdvMnt)
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
       call DestroyAcoustNewMessageSet( &
            oneGrid%AcoustNewDivSend, oneGrid%AcoustNewDivRecv, &
            oneGrid%AcoustNewPPSend, oneGrid%AcoustNewPPRecv, &
            oneGrid%AcoustNewAlphaSend, oneGrid%AcoustNewAlphaRecv, &
            oneGrid%AcoustNewThtSend, oneGrid%AcoustNewThtRecv)
       call DestroyWideGhostZoneMessageSet(&
            oneGrid%WideGhostZoneSend, oneGrid%WideGhostZoneRecv)
       call DestroyNodeDimensions(oneGrid%NodeDims)
       call DestroyNodeDimensions(oneGrid%NodeDimsAdvMnt)
       call DestroyAdvMntMessageSet(&
            oneGrid%AdvMntUVSendX, oneGrid%AdvMntUVRecvX, &
            oneGrid%AdvMntUVSendY, oneGrid%AdvMntUVRecvY, &
            oneGrid%AdvMntDxDySendX, oneGrid%AdvMntDxDyRecvX, &
            oneGrid%AdvMntDxDySendY, oneGrid%AdvMntDxDyRecvY, &
            oneGrid%AdvMntDd0SendX, oneGrid%AdvMntDd0RecvX, &
            oneGrid%AdvMntDd0SendY, oneGrid%AdvMntDd0RecvY, &
            oneGrid%AdvMntDenSendX, oneGrid%AdvMntDenRecvX, &
            oneGrid%AdvMntDenSendY, oneGrid%AdvMntDenRecvY, &
            oneGrid%AdvMntScaSendX, oneGrid%AdvMntScaRecvX, &
            oneGrid%AdvMntScaSendY, oneGrid%AdvMntScaRecvY)

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
    call DumpDomainDecomp(oneGrid%GlobalWithGhostAdvMnt, "GlobalWithGhostAdvMnt")
    call DumpDomainDecomp(oneGrid%LocalOwnAdvMnt, "LocalOwnAdvMnt")

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
    call MsgDump(h//" dumping AcoustNewDivSend")
    call DumpMessageSet(oneGrid%AcoustNewDivSend)
    call MsgDump(h//" dumping AcoustNewDivRecv")
    call DumpMessageSet(oneGrid%AcoustNewDivRecv)
    call MsgDump(h//" dumping AcoustNewPPSend")
    call DumpMessageSet(oneGrid%AcoustNewPPSend)
    call MsgDump(h//" dumping AcoustNewPPRecv")
    call DumpMessageSet(oneGrid%AcoustNewPPRecv)
    call MsgDump(h//" dumping AcoustNewAlphaSend")
    call DumpMessageSet(oneGrid%AcoustNewAlphaSend)
    call MsgDump(h//" dumping AcoustNewAlphaRecv")
    call DumpMessageSet(oneGrid%AcoustNewAlphaRecv)
    call MsgDump(h//" dumping AcoustNewThtSend")
    call DumpMessageSet(oneGrid%AcoustNewThtSend)
    call MsgDump(h//" dumping AcoustNewThtRecv")
    call DumpMessageSet(oneGrid%AcoustNewThtRecv)
    call MsgDump(h//" dumping WideGhostZoneSend")
    call DumpMessageSet(oneGrid%WideGhostZoneSend)
    call MsgDump(h//" dumping WideGhostZoneRecv")
    call DumpMessageSet(oneGrid%WideGhostZoneRecv)
    call DumpNodeDimensions(oneGrid%NodeDims, "NodeDims")
    call DumpNodeDimensions(oneGrid%NodeDimsAdvMnt, "NodeDimsAdvMnt")
    call MsgDump(h//" dumping AdvMntUVSendX")
    call DumpMessageSet(oneGrid%AdvMntUVSendX)
    call MsgDump(h//" dumping AdvMntUVRecvX")
    call DumpMessageSet(oneGrid%AdvMntUVRecvX)
    call MsgDump(h//" dumping AdvMntUVSendY")
    call DumpMessageSet(oneGrid%AdvMntUVSendY)
    call MsgDump(h//" dumping AdvMntUVRecvY")
    call DumpMessageSet(oneGrid%AdvMntUVRecvY)
    call MsgDump(h//" dumping AdvMntDxDySendX")
    call DumpMessageSet(oneGrid%AdvMntDxDySendX)
    call MsgDump(h//" dumping AdvMntDxDyRecvX")
    call DumpMessageSet(oneGrid%AdvMntDxDyRecvX)
    call MsgDump(h//" dumping AdvMntDxDySendY")
    call DumpMessageSet(oneGrid%AdvMntDxDySendY)
    call MsgDump(h//" dumping AdvMntDxDyRecvY")
    call DumpMessageSet(oneGrid%AdvMntDxDyRecvY)
    call MsgDump(h//" dumping AdvMntDd0SendX")
    call DumpMessageSet(oneGrid%AdvMntDd0SendX)
    call MsgDump(h//" dumping AdvMntDd0RecvX")
    call DumpMessageSet(oneGrid%AdvMntDd0RecvX)
    call MsgDump(h//" dumping AdvMntDd0SendY")
    call DumpMessageSet(oneGrid%AdvMntDd0SendY)
    call MsgDump(h//" dumping AdvMntDd0RecvY")
    call DumpMessageSet(oneGrid%AdvMntDd0RecvY)
    call MsgDump(h//" dumping AdvMntDenSendX")
    call DumpMessageSet(oneGrid%AdvMntDenSendX)
    call MsgDump(h//" dumping AdvMntDenRecvX")
    call DumpMessageSet(oneGrid%AdvMntDenRecvX)
    call MsgDump(h//" dumping AdvMntDenSendY")
    call DumpMessageSet(oneGrid%AdvMntDenSendY)
    call MsgDump(h//" dumping AdvMntDenRecvY")
    call DumpMessageSet(oneGrid%AdvMntDenRecvY)
    call MsgDump(h//" dumping AdvMntScaSendX")
    call DumpMessageSet(oneGrid%AdvMntScaSendX)
    call MsgDump(h//" dumping AdvMntScaRecvX")
    call DumpMessageSet(oneGrid%AdvMntScaRecvX)
    call MsgDump(h//" dumping AdvMntScaSendY")
    call DumpMessageSet(oneGrid%AdvMntScaSendY)
    call MsgDump(h//" dumping AdvMntScaRecvY")
    call DumpMessageSet(oneGrid%AdvMntScaRecvY)
    call MsgDump(h//" finishes")
  end subroutine DumpGrid


  subroutine InsertScalarTabAtOneGrid(oneGrid)
    type(Grid), pointer, intent(in) :: oneGrid

    integer :: ng
    integer :: ierr
    integer :: iEle
    integer :: nEle

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(InsertScalarTabAtOneGrid)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    ng=OneGrid%Id
    nEle=num_scalar(ng)

    allocate(oneGrid%ScalarTab(nEle), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") nEle
       write(str(2),"(i8)") ierr
       call fatal_error(h//" allocate ScalarTab("//&
            trim(adjustl(str(1)))//") fails with stat="//&
            trim(adjustl(str(2))))
    end if

    oneGrid%ScalarTabSize=nEle

    do iEle = 1, nEle
       oneGrid%ScalarTab(iEle)%name = scalar_tab(iEle,ng)%name
       oneGrid%ScalarTab(iEle)%a_var_p => scalar_tab(iEle,ng)%a_var_p
       oneGrid%ScalarTab(iEle)%a_var_t => scalar_tab(iEle,ng)%a_var_t
       oneGrid%ScalarTab(iEle)%var_p_2D => scalar_tab(iEle,ng)%var_p_2D
       oneGrid%ScalarTab(iEle)%var_p_3D => scalar_tab(iEle,ng)%var_p_3D
    end do

    if (dumpLocal) then
       write(str(1),"(i8)") nEle
       call MsgDump(h//" finishes building ScalarTab with "//&
            trim(adjustl(str(1)))//" entries")
    end if

  end subroutine InsertScalarTabAtOneGrid



  subroutine DeepCopyToScalarTabAtOneGrid(oneGrid)
    type(Grid), pointer, intent(in) :: oneGrid

    integer :: ng
    integer :: ierr
    integer :: iEle
    integer :: nEle
    integer :: dim1
    integer :: dim2
    integer :: dim3

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DeepCopyToScalarTabAtOneGrid)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    call DestroyScalarTabAtOneGrid(oneGrid)
    ng=OneGrid%Id
    nEle=num_scalar(ng)

    allocate(oneGrid%ScalarTab(nEle), stat=ierr)
    if (ierr /= 0) then
       write(str(1),"(i8)") nEle
       write(str(2),"(i8)") ierr
       call fatal_error(h//" allocate ScalarTab("//&
            trim(adjustl(str(1)))//") fails with stat="//&
            trim(adjustl(str(2))))
    end if

    oneGrid%ScalarTabSize=nEle

    do iEle = 1, nEle
       oneGrid%ScalarTab(iEle)%name = scalar_tab(iEle,ng)%name

       if (associated(scalar_tab(iEle,ng)%a_var_p)) then
          dim1 = size(scalar_tab(iEle,ng)%a_var_p)
          allocate(oneGrid%ScalarTab(iEle)%a_var_p(dim1), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") nEle
             write(str(2),"(i8)") ierr
             call fatal_error(h//" allocate a_var_p"//&
                  trim(adjustl(str(1)))//") fails with stat="//&
                  trim(adjustl(str(2))))
          end if
          oneGrid%ScalarTab(iEle)%a_var_p = scalar_tab(iEle,ng)%a_var_p
       else
          nullify(oneGrid%ScalarTab(iEle)%a_var_p)
       end if

       if (associated(scalar_tab(iEle,ng)%a_var_t)) then
          dim1 = size(scalar_tab(iEle,ng)%a_var_t)
          allocate(oneGrid%ScalarTab(iEle)%a_var_t(dim1), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") nEle
             write(str(2),"(i8)") ierr
             call fatal_error(h//" allocate a_var_t"//&
                  trim(adjustl(str(1)))//") fails with stat="//&
                  trim(adjustl(str(2))))
          end if
          oneGrid%ScalarTab(iEle)%a_var_t = scalar_tab(iEle,ng)%a_var_t
       else
          nullify(oneGrid%ScalarTab(iEle)%a_var_t)
       end if

       if (associated(scalar_tab(iEle,ng)%a_var_p_3D)) then
          dim1 = size(scalar_tab(iEle,ng)%a_var_p_3D,1)
          dim2 = size(scalar_tab(iEle,ng)%a_var_p_3D,2)
          dim3 = size(scalar_tab(iEle,ng)%a_var_p_3D,3)
          allocate(oneGrid%ScalarTab(iEle)%a_var_p_3D(dim1,dim2,dim3), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             write(str(2),"(i8)") dim1
             write(str(3),"(i8)") dim2
             write(str(4),"(i8)") dim3
             call fatal_error(h//" allocate a_var_p_3D("//&
                  trim(adjustl(str(2)))//","//&
                  trim(adjustl(str(3)))//","//&
                  trim(adjustl(str(4)))//")"//&
                  " fails with stat="//trim(adjustl(str(1))))
          end if
          oneGrid%ScalarTab(iEle)%a_var_p_3D = scalar_tab(iEle,ng)%a_var_p_3D
       else
          nullify(oneGrid%ScalarTab(iEle)%a_var_p_3D)
       end if

       if (associated(scalar_tab(iEle,ng)%a_var_t_3D)) then
          dim1 = size(scalar_tab(iEle,ng)%a_var_t_3D,1)
          dim2 = size(scalar_tab(iEle,ng)%a_var_t_3D,2)
          dim3 = size(scalar_tab(iEle,ng)%a_var_t_3D,3)
          allocate(oneGrid%ScalarTab(iEle)%a_var_t_3D(dim1,dim2,dim3), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             write(str(2),"(i8)") dim1
             write(str(3),"(i8)") dim2
             write(str(4),"(i8)") dim3
             call fatal_error(h//" allocate a_var_t_3D("//&
                  trim(adjustl(str(2)))//","//&
                  trim(adjustl(str(3)))//","//&
                  trim(adjustl(str(4)))//")"//&
                  " fails with stat="//trim(adjustl(str(1))))
          end if
          oneGrid%ScalarTab(iEle)%a_var_t_3D = scalar_tab(iEle,ng)%a_var_t_3D
       else
          nullify(oneGrid%ScalarTab(iEle)%a_var_t_3D)
       end if

       if (associated(scalar_tab(iEle,ng)%var_p_1D)) then
          dim1 = size(scalar_tab(iEle,ng)%var_p_1D,1)
          allocate(oneGrid%ScalarTab(iEle)%var_p_1D(dim1), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") nEle
             write(str(2),"(i8)") ierr
             call fatal_error(h//" allocate var_p_1D"//&
                  trim(adjustl(str(1)))//") fails with stat="//&
                  trim(adjustl(str(2))))
          end if
          oneGrid%ScalarTab(iEle)%var_p_1D = scalar_tab(iEle,ng)%var_p_1D
       else
          nullify(oneGrid%ScalarTab(iEle)%var_p_1D)
       end if

       if (associated(scalar_tab(iEle,ng)%var_p_2D)) then
          dim1 = size(scalar_tab(iEle,ng)%var_p_2D,1)
          dim2 = size(scalar_tab(iEle,ng)%var_p_2D,2)
          allocate(oneGrid%ScalarTab(iEle)%var_p_2D(dim1,dim2), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") nEle
             write(str(2),"(i8)") ierr
             call fatal_error(h//" allocate var_p_2D"//&
                  trim(adjustl(str(1)))//") fails with stat="//&
                  trim(adjustl(str(2))))
          end if
          oneGrid%ScalarTab(iEle)%var_p_2D = scalar_tab(iEle,ng)%var_p_2D
       else
          nullify(oneGrid%ScalarTab(iEle)%var_p_2D)
       end if

       if (associated(scalar_tab(iEle,ng)%var_p_3D)) then
          dim1 = size(scalar_tab(iEle,ng)%var_p_3D,1)
          dim2 = size(scalar_tab(iEle,ng)%var_p_3D,2)
          dim3 = size(scalar_tab(iEle,ng)%var_p_3D,3)
          allocate(oneGrid%ScalarTab(iEle)%var_p_3D(dim1,dim2,dim3), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") nEle
             write(str(2),"(i8)") ierr
             call fatal_error(h//" allocate var_p_3D"//&
                  trim(adjustl(str(1)))//") fails with stat="//&
                  trim(adjustl(str(2))))
          end if
          oneGrid%ScalarTab(iEle)%var_p_3D = scalar_tab(iEle,ng)%var_p_3D
       else
          nullify(oneGrid%ScalarTab(iEle)%var_p_3D)
       end if

       if (associated(scalar_tab(iEle,ng)%var_t_1D)) then
          dim1 = size(scalar_tab(iEle,ng)%var_t_1D,1)
          allocate(oneGrid%ScalarTab(iEle)%var_t_1D(dim1), stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") nEle
             write(str(2),"(i8)") ierr
             call fatal_error(h//" allocate var_t_1D"//&
                  trim(adjustl(str(1)))//") fails with stat="//&
                  trim(adjustl(str(2))))
          end if
          oneGrid%ScalarTab(iEle)%var_t_1D = scalar_tab(iEle,ng)%var_t_1D
       else
          nullify(oneGrid%ScalarTab(iEle)%var_t_1D)
       end if
    end do

    if (dumpLocal) then
       write(str(1),"(i8)") nEle
       call MsgDump(h//" finishes building ScalarTab with "//&
            trim(adjustl(str(1)))//" entries")
    end if

  end subroutine DeepCopyToScalarTabAtOneGrid





  subroutine DeepCopyFromScalarTabAtOneGrid(oneGrid)
    type(Grid), pointer, intent(in) :: oneGrid

    integer :: ng
    integer :: ierr
    integer :: iEle
    integer :: nEle

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DeepCopyFromScalarTabAtOneGrid)**"
    logical, parameter :: dumpLocal=.false.

    if (dumpLocal) then
       call MsgDump(h//" starts")
    end if

    ng=OneGrid%Id
    nEle=num_scalar(ng)

    do iEle = 1, nEle

       if (associated(scalar_tab(iEle,ng)%a_var_p)) then
          scalar_tab(iEle,ng)%a_var_p = oneGrid%ScalarTab(iEle)%a_var_p
       end if

       if (associated(scalar_tab(iEle,ng)%a_var_t)) then
          scalar_tab(iEle,ng)%a_var_t = oneGrid%ScalarTab(iEle)%a_var_t
       end if

       if (associated(scalar_tab(iEle,ng)%a_var_p_3D)) then
          scalar_tab(iEle,ng)%a_var_p_3D = oneGrid%ScalarTab(iEle)%a_var_p_3D
       end if

       if (associated(scalar_tab(iEle,ng)%a_var_t_3D)) then
          scalar_tab(iEle,ng)%a_var_t_3D = oneGrid%ScalarTab(iEle)%a_var_t_3D
       end if

       if (associated(scalar_tab(iEle,ng)%var_p_1D)) then
          scalar_tab(iEle,ng)%var_p_1D = oneGrid%ScalarTab(iEle)%var_p_1D
       end if

       if (associated(scalar_tab(iEle,ng)%var_p_2D)) then
          scalar_tab(iEle,ng)%var_p_2D = oneGrid%ScalarTab(iEle)%var_p_2D
       end if

       if (associated(scalar_tab(iEle,ng)%var_p_3D)) then
          scalar_tab(iEle,ng)%var_p_3D = oneGrid%ScalarTab(iEle)%var_p_3D
       end if

       if (associated(scalar_tab(iEle,ng)%var_t_1D)) then
          scalar_tab(iEle,ng)%var_t_1D = oneGrid%ScalarTab(iEle)%var_t_1D
       end if
    end do

    if (dumpLocal) then
       write(str(1),"(i8)") nEle
       call MsgDump(h//" finishes building ScalarTab with "//&
            trim(adjustl(str(1)))//" entries")
    end if
  end subroutine DeepCopyFromScalarTabAtOneGrid




  subroutine DestroyScalarTabAtOneGrid(oneGrid)
    type(Grid), pointer, intent(in) :: oneGrid

    integer :: iEle

    if (associated(oneGrid%ScalarTab)) then
       do iEle = 1, oneGrid%ScalarTabSize
          if (associated(oneGrid%ScalarTab(iEle)%a_var_p)) then
             deallocate(oneGrid%ScalarTab(iEle)%a_var_p)
             nullify(oneGrid%ScalarTab(iEle)%a_var_p)
          end if
          if (associated(oneGrid%ScalarTab(iEle)%a_var_t)) then
             deallocate(oneGrid%ScalarTab(iEle)%a_var_t)
             nullify(oneGrid%ScalarTab(iEle)%a_var_t)
          end if
          if (associated(oneGrid%ScalarTab(iEle)%var_p_2D)) then
             deallocate(oneGrid%ScalarTab(iEle)%var_p_2D)
             nullify(oneGrid%ScalarTab(iEle)%var_p_2D)
          end if
          if (associated(oneGrid%ScalarTab(iEle)%var_p_3D)) then
             deallocate(oneGrid%ScalarTab(iEle)%var_p_3D)
             nullify(oneGrid%ScalarTab(iEle)%var_p_3D)
          end if
       end do
       deallocate(oneGrid%ScalarTab)
       nullify(oneGrid%ScalarTab)
    end if
  end subroutine DestroyScalarTabAtOneGrid
end module ModGrid
