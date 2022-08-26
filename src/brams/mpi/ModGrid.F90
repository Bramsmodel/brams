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

  use ModScalarTable, only: &
       ScalarTable, &
       CreateScalarTab, &
       DestroyScalarTab, &
       DumpScalarTab

  use ModVarTable, only: &
       VarTable, &
       CreateVarTable, &
       DestroyVarTable, &
       DumpVarTable
  
  use ModBasicFields, only: &
       BasicFields, &
       CreateBasicFields, &
       DestroyBasicFields, &
       DumpBasicFields
  
  use ModTurbFields, only: &
       TurbFields, &
       CreateTurbFields, &
       DestroyTurbFields, &
       DumpTurbFields

  use ModControlVars, only: &
       ControlVars, &
       CreateControlVars, &
       DestroyControlVars

  use ModMicControl, only: &
       MicControl, &
       CreateMicControl, &
       DestroyMicControl, &
       DumpMicControl

  use ModMicroFields, only: &
       MicroFields, &
       CreateMicroFields, &
       DestroyMicroFields, &
       DumpMicroFields

  use ModShcuFields, only: &
       ShcuFields, &
       CreateShcuFields, &
       DestroyShcuFields, &
       DumpShcuFields

  use ModGaspartFields, only: &
       GaspartFields, &
       CreateGaspartFields, &
       DestroyGaspartFields, &
       DumpGaspartFields
  
  use meteogramType, only: &
       PolygonContainer

#ifdef JULES
  use ModJulesFields, only: &
       JulesFields, &
       CreateJulesFields, &
       DestroyJulesFields, &
       DumpJulesFields
#endif  
  
  use mem_tend, only: &
       tend

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
     type(NamelistFile), pointer :: oneNamelistFile => null()
     ! oneNamelistFile: full namelist file
     type(ControlVars), pointer :: oneControlVars => null()
     ! oneControlVars: all variables used to control flow; include all oneNamelistFile
     !          variables, replacing arrays indexed by grid number by
     !          scalars for this grid number, avoiding grid number
     !          reference; include all other variables spreaded by
     !          the code that control if statements
     type(ParallelEnvironment), pointer :: oneParallelEnvironment => null()
     ! oneParallelEnvironment: mpi size, rank and communicator for this run
     type(GridDims), pointer :: oneGridDims => null()
     ! oneGridDims: this grid dimensions as defined by namelist
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
     type(NodeDimensions), pointer :: oneNodeDimensions => null()
     ! oneNodeDimensions: indices and dimensions of this process
     ! domain decomposed sub-domain
     type(NodeDimensions), pointer :: oneNodeDimensionsAdvMnt => null()
     ! oneNodeDimensionsAdvMnt: indices and dimensions of this process
     ! domain decomposed sub-domain for use inside MonotonicAdvection
     type(BasicFields), pointer :: oneBasicFields => null()
     type(BasicFields), pointer :: oneAveBasicFields => null()

     type(TurbFields), pointer :: oneTurbFields => null()
     type(TurbFields), pointer :: oneAveTurbFields => null()

     type(MicroFields), pointer :: oneMicroFields => null()
     type(MicroFields), pointer :: oneAveMicroFields => null()

     type(ScalarTable), pointer :: oneScalarTable(:) => null()
     integer :: oneScalarTableSize=0

     type(VarTable), pointer :: oneVarTable(:) => null()
     integer :: oneVarTableSize=0
     
#ifdef JULES
     type(JulesFields), pointer :: oneJulesFields => null()
     type(JulesFields), pointer :: oneAveJulesFields => null()
#endif

     type(ShcuFields), pointer :: oneShcuFields => null()
     type(ShcuFields), pointer :: oneAveShcuFields => null()

     type(GaspartFields), pointer :: oneGaspartFields => null()
     type(GaspartFields), pointer :: oneAveGaspartFields => null()

     type(MicControl), pointer :: oneMicVars => null()
     
     ! AllGhostZoneSend/RecvG3D: Ghost Zone update at PostProcess
     ! type(MessageSet) contains all information required for
     ! ghost zone update. See description at ModMessageSet 
     type(PolygonContainer), pointer :: meteoPolygons => null()

     type(NeighbourNodes), pointer :: oneNeighbourNodes => null()
     ! oneNeighbourNodes: list of BRAMS process numbers that are neighbours
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
  end type Grid



contains



  ! CreateGrid: create and fill variable of this type



  function CreateGrid(gridId, oneNamelistFile, oneParallelEnvironment) result(oneGrid)
    integer, intent(in) :: gridId
    type(NamelistFile), pointer, intent(in) :: oneNamelistFile
    type(ParallelEnvironment), pointer, intent(in) :: oneParallelEnvironment
    type(Grid), pointer :: oneGrid

    integer :: ierr
    logical :: createAve
    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(CreateGrid)**"
    logical, parameter :: dumpLocal=.false.

    ! correctness of input arguments

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneNamelistFile")
    else if (.not. associated(oneParallelEnvironment)) then
       call fatal_error(h//" invoked with null oneParallelEnvironment")
    end if

    ! define if average fields should be created or not,
    ! according to namelist

    createAve = oneNamelistFile%avgtim > 0.0 .and. &
         (oneNamelistFile%frqmean /= 0.0 .or. &
         oneNamelistFile%frqboth /= 0.0)
    
    ! create a variable of type grid and fill entries

    allocate(oneGrid, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneGrid fails with stat="//&
            trim(adjustl(str(1))))
    end if

    ! stores input arguments

    oneGrid%id = gridId
    oneGrid%oneNamelistFile => oneNamelistFile
    oneGrid%oneParallelEnvironment => oneParallelEnvironment

    ! store GridDims extracted from OneNamelistFile 

    oneGrid%oneGridDims => CreateGridDims(gridId, &
         oneNamelistFile)

    oneGrid%oneControlVars => CreateControlVars(&
         oneGrid%oneNamelistFile, gridId)
    
    ! compute domain decomposition, obtaining
    ! cells owned by each rank and store at GlobalOwn

    oneGrid%GlobalOwn => CreateGlobalOwn(&
         oneGrid%oneGridDims, &
         oneGrid%oneParallelEnvironment, &
         "GlobalOwn" &
         )

    ! include boundary conditions (no ghost zone)

    oneGrid%GlobalOwnWithBC => CreateGlobalOwnWithBC(&
         oneGrid%oneGridDims, &
         oneGrid%oneParallelEnvironment, &
         oneGrid%GlobalOwn &
         )

    ! insert original ghost zone of widht 1
    ! at GlobalOwn and store at GlobalWithGhost

    oneGrid%GlobalWithGhost => CreateGlobalWithGhost(&
         oneGrid%oneGridDims, &
         oneGrid%oneParallelEnvironment, &
         oneGrid%GlobalOwn, &
         1, &
         "GlobalWithGhost" &
         )

    ! convert global indices from GlobalWithGhost
    ! into local indices stored at LocalOwn

    oneGrid%LocalOwn => CreateLocalOwn(&
         oneGrid%oneParallelEnvironment, &
         oneGrid%GlobalWithGhost, &
         oneGrid%GlobalOwn, &
         "LocalOwn" &
         )

    ! neighbour nodes for original ghost zone update operations

    oneGrid%oneNeighbourNodes => CreateNeighbourNodes(&
         ParEnv=oneGrid%oneParallelEnvironment, &
         GlobalOwn=oneGrid%GlobalOwn, &
         GlobalWithGhost=oneGrid%GlobalWithGhost, &
         varName="oneGrid%oneNeighbourNodes" &
         )

    ! this node dimensions and indexing limits

    oneGrid%oneNodeDimensions => CreateNodeDimensions(&
         GridSize=oneGrid%oneGridDims, &
         ParEnv=oneGrid%oneParallelEnvironment, &
         LocalOwn=oneGrid%LocalOwn, &
         GlobalOwn=oneGrid%GlobalOwn, &
         verticalGhostZoneWidth=0, &
         surfaceGhostZoneWidth=1, &
         varName="oneNodeDimensions" &
         )

    ! for MonotonicAdvection, insert ghost zone of parametrized widht
    ! at GlobalOwn and store at GlobalWithGhostAdvMnt

    oneGrid%GlobalWithGhostAdvMnt => CreateGlobalWithGhost(&
         GridSize=oneGrid%oneGridDims, &
         ParEnv=oneGrid%oneParallelEnvironment, &
         GlobalOwn=oneGrid%GlobalOwn, &
         GhostZoneWidth=oneNamelistFile%ghostzonelength, &
         varName="GlobalWithGhostAdvMnt" &
         )

    ! convert global indices from GlobalWithGhostAdvMnt
    ! into local indices stored at LocalOwnAdvMnt

    oneGrid%LocalOwnAdvMnt => CreateLocalOwn(&
         ParEnv=oneGrid%oneParallelEnvironment, &
         GlobalWithGhost=oneGrid%GlobalWithGhostAdvMnt, &
         GlobalOwn=oneGrid%GlobalOwn, &
         varName="LocalOwnAdvMnt" &
         )

    ! this node dimensions and indexing limits

    oneGrid%oneNodeDimensionsAdvMnt => CreateNodeDimensions(&
         GridSize=oneGrid%oneGridDims, &
         ParEnv=oneGrid%oneParallelEnvironment, &
         LocalOwn=oneGrid%LocalOwnAdvMnt, &
         GlobalOwn=oneGrid%GlobalOwn, &
         verticalGhostZoneWidth=0, &
         surfaceGhostZoneWidth=oneNamelistFile%ghostzonelength, &
         varName="oneNodeDimensionsAdvMnt" &
         )

    ! this node Basic Fields

    oneGrid%oneBasicFields => CreateBasicFields(&
         oneGrid%oneNodeDimensions, &
         oneGrid%oneNamelistFile)
    if (createAve) then
       oneGrid%oneAveBasicFields => CreateBasicFields(&
            oneGrid%oneNodeDimensions, &
            oneGrid%oneNamelistFile)
    else
       ! oneAveBasicFields is created with null components
       allocate(oneGrid%oneAveBasicFields, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate oneGrid%oneAveBasicFields fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    ! this node Turb Fields

    oneGrid%oneTurbFields => CreateTurbFields(&
         oneGrid%oneNodeDimensions, &
         oneGrid%oneNamelistFile, &
         gridId)
    if (createAve) then
       oneGrid%oneAveTurbFields => CreateTurbFields(&
            oneGrid%oneNodeDimensions, &
            oneGrid%oneNamelistFile, &
            gridId)
    else
       ! oneAveTurbFields is created with null components
       allocate(oneGrid%oneAveTurbFields, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate oneGrid%oneAveTurbFields fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if

    ! this node MicControl

    oneGrid%oneMicVars => CreateMicControl(oneGrid%oneNamelistFile)

    ! this node Micro Fields

    oneGrid%oneMicroFields => CreateMicroFields(&
         oneGrid%Id, &
         oneGrid%oneNamelistFile, &
         oneGrid%oneNodeDimensions, &
         oneGrid%oneMicVars)
    if (createAve) then
       oneGrid%oneAveMicroFields => CreateMicroFields(&
            oneGrid%Id, &
            oneGrid%oneNamelistFile, &
            oneGrid%oneNodeDimensions, &
            oneGrid%oneMicVars)
    end if

    ! this node Scalar Table

    oneGrid%oneScalarTable => CreateScalarTab()
    oneGrid%oneScalarTableSize = 0

    ! this node Var Table

    oneGrid%oneVarTable => CreateVarTable()
    oneGrid%oneVarTableSize = 0

    ! this node Jules Fields

#ifdef JULES
    if (oneGrid%oneNamelistFile%isfcl == 5) then
       oneGrid%oneJulesFields => CreateJulesFields(&
            oneGrid%oneNodeDimensions, &
            oneGrid%oneNamelistFile)
       if (createAve) then
          oneGrid%oneAveJulesFields => CreateJulesFields(&
               oneGrid%oneNodeDimensions, &
               oneGrid%oneNamelistFile)
       else
          ! oneAveJulesFields is created with null components
          allocate(oneGrid%oneAveJulesFields, stat=ierr)
          if (ierr /= 0) then
             write(str(1),"(i8)") ierr
             call fatal_error(h//" allocate oneGrid%oneAveJulesFields fails with stat="//&
                  trim(adjustl(str(1))))
          end if
       end if
    end if
#endif

    ! this node Shcu Fields
    
    oneGrid%oneShcuFields => CreateShcuFields(&
         oneGrid%oneNodeDimensions, &
         oneGrid%oneControlVars)
    if (createAve) then
       oneGrid%oneAveShcuFields => CreateShcuFields(&
            oneGrid%oneNodeDimensions, &
            oneGrid%oneControlVars)
    end if

    ! this node Gaspart Fields

    oneGrid%oneGaspartFields => CreateGaspartFields(&
         oneGrid%oneNodeDimensions, &
         oneGrid%oneNamelistFile)
    if (createAve) then
       oneGrid%oneAveGaspartFields => CreateGaspartFields(&
            oneGrid%oneNodeDimensions, &
            oneGrid%oneNamelistFile)
    else
       ! oneAveGaspartFields is created with null components
       allocate(oneGrid%oneAveGaspartFields, stat=ierr)
       if (ierr /= 0) then
          write(str(1),"(i8)") ierr
          call fatal_error(h//" allocate oneGrid%oneAveGaspartFields fails with stat="//&
               trim(adjustl(str(1))))
       end if
    end if
    
    if (dumpLocal) then
       call MsgDump(h//" dumping OneGrid at the end of CreateGrid")
       call DumpGrid(OneGrid)
    end if
  end function CreateGrid




  subroutine InsertMessageSetAtOneGrid(oneGrid)
    type(Grid), pointer :: oneGrid

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

    if (.not. associated(oneGrid)) then
       call fatal_error(h//" invoked with null grid")
    end if

    call CreateAcousticMessageSet(&
         oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
         oneGrid%Id, &
         oneGrid%oneGridDims, oneGrid%oneParallelEnvironment, oneGrid%oneNeighbourNodes, &
         oneGrid%GlobalOwn, &
         oneGrid%GlobalOwnWithBC, &
         oneGrid%GlobalWithGhost, &
         oneGrid%AcouSendU, oneGrid%AcouRecvU, TagU, &
         oneGrid%AcouSendV, oneGrid%AcouRecvV, TagV, &
         oneGrid%AcouSendPNorth, oneGrid%AcouRecvPNorth, TagPNorth, &
         oneGrid%AcouSendPEast, oneGrid%AcouRecvPEast, TagPEast, &
         oneGrid%AcouSendUV, oneGrid%AcouRecvUV, TagUV, &
         oneGrid%AcouSendWP, oneGrid%AcouRecvWP, TagWP)

    call CreateDn0MessageSet(&
         oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
         oneGrid%Id, &
         oneGrid%oneGridDims, oneGrid%oneParallelEnvironment, oneGrid%oneNeighbourNodes, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, &
         oneGrid%SendDn0u, oneGrid%RecvDn0u, TagDn0u, &
         oneGrid%SendDn0v, oneGrid%RecvDn0v, TagDn0v)

    call CreateG3DMessageSet(&
         oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
         oneGrid%Id, &
         oneGrid%oneGridDims, oneGrid%oneParallelEnvironment, oneGrid%oneNeighbourNodes, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhost, &
         oneGrid%oneNamelistFile, &
         oneGrid%SendG3D, oneGrid%RecvG3D, TagG3D)


    call CreateSelectedGhostZoneMessageSet(&
         oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
         oneGrid%Id, &
         oneGrid%oneGridDims, oneGrid%oneParallelEnvironment, oneGrid%oneNeighbourNodes, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhost, &
         oneGrid%SelectedGhostZoneSend, &
         oneGrid%SelectedGhostZoneRecv, &
         TagSelectedGhostZone)

    call CreateAllGhostZoneMessageSet(&
         oneGrid%oneVarTable, oneGrid%oneVarTableSize, &
         oneGrid%Id, &
         oneGrid%oneGridDims, oneGrid%oneParallelEnvironment, oneGrid%oneNeighbourNodes, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhost, &
         oneGrid%AllGhostZoneSend, &
         oneGrid%AllGhostZoneRecv, &
         TagAllGhostZone)

    ! use desired bounds fields to create AcoustNew Message Sets;
    ! correct field memory address by invoking UpdateFieldAdress
    ! prior to use the Message Sets

    call CreateAcoustNewMessageSet(&
         oneGrid%oneGridDims, oneGrid%oneParallelEnvironment, oneGrid%oneNeighbourNodes, &
         oneGrid%GlobalOwn, oneGrid%GlobalWithGhost, oneGrid%oneNodeDimensions, &
         TagAcoustNewDiv, oneGrid%AcoustNewDivSend, oneGrid%AcoustNewDivRecv, &
         TagAcoustNewPP, oneGrid%AcoustNewPPSend, oneGrid%AcoustNewPPRecv, &
         TagAcoustNewAlpha, oneGrid%AcoustNewAlphaSend, oneGrid%AcoustNewAlphaRecv, &
         TagAcoustNewTht, oneGrid%AcoustNewThtSend, oneGrid%AcoustNewThtRecv, &
         tend%tht_rk)


    call CreateWideGhostZoneMessageSet(&
         oneGrid%oneParallelEnvironment, oneGrid%oneNeighbourNodes, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhost, oneGrid%oneNodeDimensions, &
         1, oneGrid%oneNodeDimensions%mzp, &
         TagWideGhostZone, oneGrid%WideGhostZoneSend, oneGrid%WideGhostZoneRecv)

    call CreateAdvMntMessageSet(&
         oneGrid%oneParallelEnvironment, oneGrid%oneNeighbourNodes, &
         oneGrid%GlobalOwnWithBC, oneGrid%GlobalWithGhostAdvMnt, &
         oneGrid%oneNodeDimensions, oneGrid%oneNodeDimensionsAdvMnt, &
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
       call MsgDump(h//" dumping oneGrid at the end of InsertMessageSetAtOneGrid")
       call DumpGrid(OneGrid)
    end if
  end subroutine InsertMessageSetAtOneGrid



  ! DestroyGrid: deallocate area of a variable of type grid



  subroutine DestroyGrid(oneGrid)
    type(Grid), pointer :: oneGrid

    character(len=*), parameter :: h="**(DestroyGrid)**"

    if (associated(oneGrid)) then
       call DestroyGridDims(oneGrid%oneGridDims)
       call DestroyControlVars(oneGrid%oneControlVars)
       call DestroyDomainDecomp(oneGrid%GlobalOwn)
       call DestroyDomainDecomp(oneGrid%GlobalOwnWithBC)
       call DestroyDomainDecomp(oneGrid%GlobalWithGhost)
       call DestroyDomainDecomp(oneGrid%LocalOwn)
       call DestroyNeighbourNodes(oneGrid%oneNeighbourNodes)
       call DestroyNodeDimensions(oneGrid%oneNodeDimensions)
       call DestroyDomainDecomp(oneGrid%GlobalWithGhostAdvMnt)
       call DestroyDomainDecomp(oneGrid%LocalOwnAdvMnt)
       call DestroyNodeDimensions(oneGrid%oneNodeDimensionsAdvMnt)
       call DestroyBasicFields(oneGrid%oneBasicFields)
       call DestroyBasicFields(oneGrid%oneAveBasicFields)
       call DestroyTurbFields(oneGrid%oneTurbFields)
       call DestroyTurbFields(oneGrid%oneAveTurbFields)
       call DestroyMicroFields(oneGrid%oneMicroFields)
       call DestroyMicroFields(oneGrid%oneAveMicroFields)
       call DestroyScalarTab(oneGrid%oneScalarTable)
       oneGrid%oneScalarTableSize=0
       call DestroyVarTable(oneGrid%oneVarTable)
       oneGrid%oneVarTableSize=0
       call DestroyMicControl(oneGrid%oneMicVars)
#ifdef JULES
       call DestroyJulesFields(oneGrid%oneJulesFields)
       call DestroyJulesFields(oneGrid%oneAveJulesFields)
#endif
       call DestroyShcuFields(oneGrid%oneShcuFields)
       call DestroyShcuFields(oneGrid%oneAveShcuFields)
       call DestroyGaspartFields(oneGrid%oneGaspartFields)
       call DestroyGaspartFields(oneGrid%oneAveGaspartFields)
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

    character(len=8) :: str(10)
    character(len=*), parameter :: h="**(DumpGrid)**"

    if (.not. associated(oneGrid)) then
       call fatal_error(h//" invoked with null oneGrid")
    else if (.not. associated(oneGrid%oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneGrid%oneNamelistFile")
    else if (.not. associated(oneGrid%oneParallelEnvironment)) then
       call fatal_error(h//" invoked with null oneGrid%oneParallelEnvironment")
    end if

    write(str(1),"(i8)") oneGrid%Id
    call MsgDump(h//" for grid "//trim(adjustl(str(1))))

    call MsgDump(h//" dumping component oneGridDims")
    call DumpGridDims(oneGrid%oneGridDims)

    call MsgDump(h//" dumping domain decomposed components")
    call DumpDomainDecomp(oneGrid%GlobalOwn, "GlobalOwn")
    call DumpDomainDecomp(oneGrid%GlobalOwnWithBC, "GlobalOwnWithBC")
    call DumpDomainDecomp(oneGrid%GlobalWithGhost, "GlobalWithGhost")
    call DumpDomainDecomp(oneGrid%LocalOwn, "LocalOwn")
    call DumpDomainDecomp(oneGrid%GlobalWithGhostAdvMnt, "GlobalWithGhostAdvMnt")
    call DumpDomainDecomp(oneGrid%LocalOwnAdvMnt, "LocalOwnAdvMnt")
    
    call MsgDump(h//" dumping neighborhood components")
    call DumpNeighbourNodes(oneGrid%oneNeighbourNodes,"oneGrid%oneNeighbourNodes")

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
    call DumpNodeDimensions(oneGrid%oneNodeDimensions, "oneNodeDimensions")
    call DumpNodeDimensions(oneGrid%oneNodeDimensionsAdvMnt, "oneNodeDimensionsAdvMnt")
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


    call DumpBasicFields(oneGrid%oneBasicFields, "oneGrid%oneBasicFields")
    call DumpBasicFields(oneGrid%oneAveBasicFields, "oneGrid%oneAveBasicFields")
    call DumpTurbFields(oneGrid%oneTurbFields, "oneGrid%oneTurbFields")
    call DumpTurbFields(oneGrid%oneAveTurbFields, "oneGrid%oneAveTurbFields")
    call DumpMicroFields(oneGrid%oneMicroFields, "oneGrid%oneMicroFields")
    call DumpMicroFields(oneGrid%oneAveMicroFields, "oneGrid%oneAveMicroFields")
    call DumpVarTable(oneGrid%oneVarTable, oneGrid%oneVarTableSize, "oneGrid%oneVarTable")
    call MsgDump(h//" dumping MicControl")
    call DumpMicControl(oneGrid%oneMicVars)
#ifdef JULES
    call DumpJulesFields(oneGrid%oneJulesFields, "oneGrid%oneJulesFields")
    call DumpJulesFields(oneGrid%oneAveJulesFields, "oneGrid%oneAveJulesFields")
#endif
    call DumpShcuFields(oneGrid%oneShcuFields, "oneGrid%oneShcuFields")
    call DumpShcuFields(oneGrid%oneAveShcuFields, "oneGrid%oneAveShcuFields")
    call DumpGaspartFields(oneGrid%oneGaspartFields, "oneGrid%oneGaspartFields")
    call DumpGaspartFields(oneGrid%oneAveGaspartFields, "oneGrid%oneAveGaspartFields")
    call MsgDump(h//" finishes")
  end subroutine DumpGrid
end module ModGrid
