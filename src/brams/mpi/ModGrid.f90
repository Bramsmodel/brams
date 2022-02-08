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
       DestroyAllGhostZoneMessageSet


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
     type(DomainDecomp), pointer :: GlobalWithGhost => null()
     ! GlobalWithGhost: global indices of this grid domain
     !                  decomposition at each rank, including
     !                  the owned points and a ghost zone of length one.
     !                  Not a domain partition, due to ghost
     !                  zone inclusion..
     type(DomainDecomp), pointer :: LocalOwn => null()
     ! LocalOwn: local indices of this grid domain
     !           decomposition owned by each rank. 
     !           Convertion of GlobalOwn to local indices
     type(NeighbourNodes), pointer :: Neigh => null()
     ! Neigh: list of BRAMS process numbers that are neighbours
     !        of this node for usual ghost zone update operations
     type(DomainDecomp), pointer :: GlobalWithGhostAdvectc_rk => null()
     ! GlobalWithGhostAdvectc_rk: global indices of this grid domain
     !                            decomposition considering the ghost zone
     !                            width required by procedure advectc_rk 
     !                            for all ranks
     type(DomainDecomp), pointer :: LocalOwnAdvectc_rk => null()
     ! LocalOwnAdvectc_rk: local indices of this grid domain
     !                     decomposition considering the ghost zone
     !                     width required by procedure advectc_rk 
     !                     for all ranks.
     !                     Convertion of GlobalWithGhostAdvectc_rk
     !                     to local indices
     type(NeighbourNodes), pointer :: NeighAdvectc_rk => null()
     ! NeighAdvectc_rk: list of BRAMS process numbers that are neighbours
     !                  of this node for Advectc_rk ghost zone updates
     type(MessageSet), pointer :: AcouSendU
     type(MessageSet), pointer :: AcouRecvU
     type(MessageSet), pointer :: AcouSendV
     type(MessageSet), pointer :: AcouRecvV
     type(MessageSet), pointer :: AcouSendP
     type(MessageSet), pointer :: AcouRecvP
     type(MessageSet), pointer :: AcouSendUV
     type(MessageSet), pointer :: AcouRecvUV
     type(MessageSet), pointer :: AcouSendWP
     type(MessageSet), pointer :: AcouRecvWP
     type(MessageSet), pointer :: SendDn0u
     type(MessageSet), pointer :: RecvDn0u
     type(MessageSet), pointer :: SendDn0v
     type(MessageSet), pointer :: RecvDn0v
     type(MessageSet), pointer :: SendG3D
     type(MessageSet), pointer :: RecvG3D
     type(MessageSet), pointer :: SelectedGhostZoneSend
     type(MessageSet), pointer :: SelectedGhostZoneRecv
     type(MessageSet), pointer :: AllGhostZoneSend
     type(MessageSet), pointer :: AllGhostZoneRecv
     type(PolygonContainer), pointer :: meteoPolygons
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

    ! insert original ghost zone of required by
    ! procedure advectc_rk at GlobalOwn and store
    ! at GlobalWithGhostAdvectc_rk
    
    oneGrid%GlobalWithGhostAdvectc_rk => CreateGlobalWithGhost(&
         GridSize=oneGrid%GridSize, &
         ParEnv=oneGrid%ParEnv, &
         GlobalOwn=oneGrid%GlobalOwn, &
         GhostZoneWidth=3, &
         varName="GlobalWithGhostAdvectc_rk" &
         )

    ! neighbour nodes for original ghost zone update operations
    
    oneGrid%Neigh => CreateNeighbourNodes(&
         ParEnv=oneGrid%ParEnv, &
         GlobalOwn=oneGrid%GlobalOwn, &
         GlobalWithGhost=oneGrid%GlobalWithGhost, &
         varName="oneGrid%Neigh" &
         )


    ! convert global indices from GlobalWithGhostAdvectc_rk
    ! into local indices stored at LocalOwnAdvectc_rk

    oneGrid%LocalOwnAdvectc_rk => CreateLocalOwn(&
         ParEnv=oneGrid%ParEnv, &
         GlobalWithGhost=oneGrid%GlobalWithGhostAdvectc_rk, &
         GlobalOwn=oneGrid%GlobalOwn, &
         varName="LocalOwnAdvectc_rk" &
         )

    ! neighbour nodes for Advectc_rk ghost zone update operations

    oneGrid%NeighAdvectc_rk => CreateNeighbourNodes(&
         ParEnv=oneGrid%ParEnv, &
         GlobalOwn=oneGrid%GlobalOwn, &
         GlobalWithGhost=oneGrid%GlobalWithGhostAdvectc_rk, &
         varName="oneGrid%NeighAdvectc_rk" &
         )
	 
    oneGrid%AcouSendU => null()
    oneGrid%AcouRecvU => null()
    oneGrid%AcouSendV => null()
    oneGrid%AcouRecvV => null()
    oneGrid%AcouSendP => null()
    oneGrid%AcouRecvP => null()
    oneGrid%AcouSendUV => null()
    oneGrid%AcouRecvUV => null()
    oneGrid%AcouSendWP => null()
    oneGrid%AcouRecvWP => null()
    oneGrid%SendDn0u => null()
    oneGrid%RecvDn0u => null()
    oneGrid%SendDn0v => null()
    oneGrid%RecvDn0v => null()
    oneGrid%SendG3D => null()
    oneGrid%RecvG3D => null()
    oneGrid%SelectedGhostZoneSend => null()
    oneGrid%SelectedGhostZoneRecv => null()
    oneGrid%AllGhostZoneSend => null()
    oneGrid%AllGhostZoneRecv => null()
    oneGrid%meteoPolygons => null()
    
  end subroutine CreateGrid





  subroutine InsertMessageSetAtOneGrid(oneGrid)
    type(Grid), pointer :: oneGrid

    character(len=16) :: c0, c1
    character(len=*), parameter :: h="**(InsertMessageSetAtOneGrid)**"
    logical, parameter :: dumpLocal=.false.

    if (.not. associated(oneGrid)) then
       call fatal_error(h//" invoked with null grid")
    end if

    call CreateAcousticMessageSet(oneGrid%Id, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, &
         oneGrid%GlobalWithGhost, &
         oneGrid%AcouSendU, oneGrid%AcouRecvU, &
         oneGrid%AcouSendV, oneGrid%AcouRecvV,&
         oneGrid%AcouSendP, oneGrid%AcouRecvP, &
         oneGrid%AcouSendUV, oneGrid%AcouRecvUV, &
         oneGrid%AcouSendWP, oneGrid%AcouRecvWP)

    call CreateDn0MessageSet(oneGrid%Id, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
         oneGrid%SendDn0u, oneGrid%RecvDn0u, &
         oneGrid%SendDn0v, oneGrid%RecvDn0v)

    call CreateG3DMessageSet(oneGrid%Id, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
         oneGrid%Ramsin, &
         oneGrid%SendG3D, oneGrid%RecvG3D)

    ! temporariamente, num_var e vtab_r sao variaveis globais,
    ! enquanto nao inclusas no tipo Grid

    call CreateSelectedGhostZoneMessageSet(&
       oneGrid%Id, num_var, vtab_r, &
       oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
       oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
       oneGrid%SelectedGhostZoneSend, oneGrid%SelectedGhostZoneRecv)

    call CreateAllGhostZoneMessageSet(&
       oneGrid%Id, num_var, vtab_r, &
       oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
       oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
       oneGrid%AllGhostZoneSend, oneGrid%AllGhostZoneRecv)

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
            oneGrid%AcouSendP, oneGrid%AcouRecvP, &
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
    call DumpGridDims(oneGrid%GridSize, h)

    call MsgDump(h//" dumping domain decomposed components")
    call DumpDomainDecomp(oneGrid%GlobalOwn, "GlobalOwn")
    call DumpDomainDecomp(oneGrid%GlobalWithGhost, "GlobalWithGhost")
    call DumpDomainDecomp(oneGrid%LocalOwn, "LocalOwn")
    call DumpDomainDecomp(oneGrid%GlobalWithGhostAdvectc_rk, "GlobalWithGhostAdvectc_rk")
    call DumpDomainDecomp(oneGrid%LocalOwnAdvectc_rk, "LocalOwnAdvectc_rk")

    call MsgDump(h//" dumping neighborhood")
    call DumpNeighbourNodes(oneGrid%NeighAdvectc_rk,"oneGrid%NeighAdvectc_rk")

    call MsgDump(h//" dumping AcouSendU")
    call DumpMessageSet(oneGrid%AcouSendU)
    call MsgDump(h//" dumping AcouRecvU")
    call DumpMessageSet(oneGrid%AcouRecvU)
    call MsgDump(h//" dumping AcouSendV")
    call DumpMessageSet(oneGrid%AcouSendV)
    call MsgDump(h//" dumping AcouRecvV")
    call DumpMessageSet(oneGrid%AcouRecvV)
    call MsgDump(h//" dumping AcouSendP")
    call DumpMessageSet(oneGrid%AcouSendP)
    call MsgDump(h//" dumping AcouRecvP")
    call DumpMessageSet(oneGrid%AcouRecvP)
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
  end subroutine DumpGrid
end module ModGrid
