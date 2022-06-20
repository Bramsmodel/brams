mem_stilt.o : $(STILT)/mem_stilt.f90 var_tables.o io_params.o grid_dims.o \
	ModNamelistFile.o rconstants.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfFile.o : $(FDDA)/ModVarfFile.f90 ModBasicFields.o mem_chem1.o node_mod.o \
	mem_leaf.o ModNudAnalysis.o ModMessageSet.o mem_aer1.o micphys.o parlibf.o \
	ModRcio.o ModDateUtils.o rconstants.o mem_scratch.o ModGrid.o mem_varinit.o \
	mem_grid.o ReadBcst.o ModControlVars.o aer1_list.o ModGridTree.o ref_sounding.o \
	isan_coms.o chem1_list.o ModRamsReadHeader.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBuffering.o : $(MPI)/ModBuffering.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNodeDimensions.o : $(MPI)/ModNodeDimensions.f90 ModDomainDecomp.o \
	ModGridDims.o ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageData.o : $(MPI)/ModMessageData.f90 ModFieldSectionList.o \
	ModFieldSection.o ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageSet.o : $(MPI)/ModMessageSet.f90 ModMessageData.o ModFieldSection.o \
	var_tables.o ModFieldSectionList.o ModParallelEnvironment.o mem_grid.o \
	parlibf.o ModNeighbourNodes.o ModNodeDimensions.o ModNamelistFile.o \
	ModGridDims.o ModDomainDecomp.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridDims.o : $(MPI)/ModGridDims.f90 ModNamelistFile.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGrid.o : $(MPI)/ModGrid.f90 ModTurbFields.o var_tables.o mem_tend.o \
	ModBasicFields.o ModParallelEnvironment.o ModNeighbourNodes.o \
	ModNodeDimensions.o ModNamelistFile.o ModGridDims.o ModMessageSet.o \
	ModScalarTable.o ModDomainDecomp.o meteogramType.o ModControlVars.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridTree.o : $(MPI)/ModGridTree.f90 ModGrid.o ModNamelistFile.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNeighbourNodes.o : $(MPI)/ModNeighbourNodes.f90 ModDomainDecomp.o \
	ModGridDims.o ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModFieldSection.o : $(MPI)/ModFieldSection.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModFieldSectionList.o : $(MPI)/ModFieldSectionList.f90 ModFieldSection.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDomainDecomp.o : $(MPI)/ModDomainDecomp.f90 ModGridDims.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBasicFields.o : $(MEMORY)/ModBasicFields.f90 var_tables.o mem_stilt.o \
	ModParallelEnvironment.o ModNodeDimensions.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ReadBcst.o : $(MPI)/ReadBcst.f90 ModTurbFields.o var_tables.o an_header.o \
	ModBasicFields.o mem_grid.o parlibf.o node_mod.o mem_aerad.o ModControlVars.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mpi_io_engine-5d.o : $(IO)/mpi_io_engine-5d.f90 ModParaInit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

an_header.o : $(UTILS_MODS)/an_header.f90 $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMemory.o : $(UTILS_LIB)/ModMemory.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

getvar.o : $(UTILS_LIB)/getvar.f90 dump.o an_header.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

sfclyr_jules.o : $(JULES_DIR)/sfclyr_jules.f90 mem_radiate.o mem_cuparm.o \
	ModBasicFields.o io_params.o mem_chem1.o node_mod.o io_constants.o \
	model_time_mod.o mem_leaf.o csigma_mod.o gridbox_mean_mod.o micphys.o \
	leaf_coms.o mem_brams_jules.o rconstants.o datetime_mod.o mem_jules.o \
	gridmean_fluxes.o jules_surface_types_mod.o mem_grid.o mem_micro.o ancil_info.o \
	fluxes.o ModTurbFields.o jules_fields_mod.o chem1_list.o sf_diags_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDateUtils.o : $(UTILS_MODS)/ModDateUtils.f90 dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cup_dn.o : $(CUPARM)/cup_dn.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ccatt_start.o : $(CCATT)/ccatt_start.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

machine_arq.o : $(MODEL)/machine_arq.F90 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(CMACH) $(<F:.F90=.F90)
	rm -f $(<F:.f90=.f90)

mem_grid_dim_defs.o : $(MEMORY)/mem_grid_dim_defs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

teb_spm_start.o : $(TEB_SPM)/teb_spm_start.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

grid_dims.o : $(MEMORY)/grid_dims.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

io_params.o : $(IO)/io_params.f90 ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ref_sounding.o : $(MODEL)/ref_sounding.f90 ModNamelistFile.o grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

micphys.o : $(MICRO)/micphys.f90 ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModScalarTable.o : $(MEMORY)/ModScalarTable.f90 ModParallelEnvironment.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

var_tables.o : $(MEMORY)/var_tables.f90 ModParallelEnvironment.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_cuparm.o : $(CUPARM)/mem_cuparm.f90 ModNamelistFile.o var_tables.o \
	grid_dims.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_turb_scalar.o : $(TURB)/mem_turb_scalar.f90 var_tables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_precision.o : $(RADIATE)/mem_precision.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tuv.o : $(TUV)/mem_tuv.f90 mem_globrad.o var_tables.o mem_stilt.o \
	ModTuv2.7.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

tuvParameter.o : $(TUV)/tuvParameter.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModTuv2.7.o : $(TUV)/ModTuv2.7.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuvDriver2.7.o : $(TUV)/ModTuvDriver2.7.f90 mem_globrad.o mem_tuv.o \
	mem_radiate.o ModBasicFields.o mem_grid.o mem_chem1.o ref_sounding.o \
	chem_fastjx_driv.o node_mod.o chem1_list.o tuvParameter.o mem_rrtm.o mem_leaf.o \
	mem_aerad.o extra.o rconstants.o mem_carma.o ModTuv2.7.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aerad.o : $(RADIATE)/mem_aerad.f90 mem_grid_dim_defs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globrad.o : $(RADIATE)/mem_globrad.f90 mem_precision.o mem_radiate.o \
	parlibf.o ModNamelistFile.o mem_aerad.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globaer.o : $(RADIATE)/mem_globaer.f90 mem_precision.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_carma.o : $(RADIATE)/mem_carma.f90 mem_globrad.o ModTurbFields.o \
	var_tables.o ModBasicFields.o io_params.o mem_grid.o grid_dims.o parlibf.o \
	node_mod.o mem_scalar.o ModNamelistFile.o mem_aerad.o ReadBcst.o \
	ModControlVars.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_leaf.o : $(SURFACE)/mem_leaf.f90 teb_spm_start.o var_tables.o io_params.o \
	grid_dims.o ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_jules.o : $(SURFACE)/mem_jules.f90 io_params.o ModNamelistFile.o \
	var_tables.o grid_dims.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_micro.o : $(MICRO)/mem_micro.f90 mem_radiate.o mem_cuparm.o var_tables.o \
	micphys.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_radiate.o : $(RADIATE)/mem_radiate.f90 ModNamelistFile.o var_tables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scalar.o : $(MEMORY)/mem_scalar.f90 var_tables.o ModNamelistFile.o \
	io_params.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch.o : $(MEMORY)/mem_scratch.f90 mem_radiate.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch1_brams.o : $(MEMORY)/mem_scratch1_brams.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbFields.o : $(TURB)/ModTurbFields.f90 var_tables.o \
	ModParallelEnvironment.o ModNodeDimensions.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tend.o : $(MEMORY)/mem_tend.f90 teb_spm_start.o ModTurbFields.o mem_emiss.o \
	ModBasicFields.o mem_micro.o mem_grid.o mem_scalar.o ModNamelistFile.o \
	ModScalarTable.o mem_gaspart.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_varinit.o : $(MEMORY)/mem_varinit.f90 var_tables.o grid_dims.o mem_chem1.o \
	chem1_list.o ModNamelistFile.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid.o : $(MEMORY)/mem_grid.f90 ModNamelistFile.o var_tables.o grid_dims.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_nestb.o : $(NESTING)/mem_nestb.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_mksfc.o : $(MKSFC)/mem_mksfc.f90 teb_spm_start.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_micro_optij.o : $(MICRO)/mem_micro_optij.f90 micphys.o mem_radiate.o \
	mem_grid.o mem_micro.o grid_dims.o node_mod.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rconstants.o : $(MEMORY)/rconstants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

node_mod.o : $(MPI)/node_mod.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

conv_coms.o : $(CUPARM)/conv_coms.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ke_coms.o : $(TURB)/ke_coms.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

leaf_coms.o : $(SURFACE)/leaf_coms.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_oda.o : $(FDDA)/mem_oda.f90 ModNamelistFile.o var_tables.o grid_dims.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

obs_input.o : $(FDDA)/obs_input.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

grid_struct.o : $(MEMORY)/grid_struct.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_opt_scratch.o : $(TURB)/mem_opt_scratch.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

Phys_const.o : $(CUPARM)/Phys_const.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cup_up.o : $(CUPARM)/cup_up.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cup_output_vars.o : $(CUPARM)/cup_output_vars.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

shcu_vars_const.o : $(CUPARM)/shcu_vars_const.f90 conv_coms.o ModNamelistFile.o \
	grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_shcu.o : $(CUPARM)/mem_shcu.f90 var_tables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell_param2.o : $(CUPARM)/mem_grell_param2.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell.o : $(CUPARM)/mem_grell.f90 mem_cuparm.o var_tables.o \
	shcu_vars_const.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tconv.o : $(CCATT)/mem_tconv.f90 chem1_list.o aer1_list.o mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch1_grell.o : $(CUPARM)/mem_scratch1_grell.f90 ccatt_start.o dump.o \
	mem_grell_param2.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch2_grell.o : $(CUPARM)/mem_scratch2_grell.f90 node_mod.o \
	mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch2_grell_sh.o : $(CUPARM)/mem_scratch2_grell_sh.f90 node_mod.o \
	mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch3_grell.o : $(CUPARM)/mem_scratch3_grell.f90 mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch3_grell_sh.o : $(CUPARM)/mem_scratch3_grell_sh.f90 mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_cutrans.o : $(CUPARM)/mem_cutrans.f90 mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

isan_coms.o : $(ISAN_MODS)/isan_coms.f90 ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

oda_read.o : $(FDDA)/oda_read.f90 isan_coms.o mem_oda.o ModDateUtils.o \
	mem_grid.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

oda_krig.o : $(FDDA)/oda_krig.f90 mem_oda.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaNudge.o : $(FDDA)/ModOdaNudge.f90 mem_tend.o ModBasicFields.o io_params.o \
	mem_grid.o mem_oda.o node_mod.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

oda_proc_obs.o : $(FDDA)/oda_proc_obs.f90 mem_oda.o rconstants.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

oda_sta_count.o : $(FDDA)/oda_sta_count.f90 mem_oda.o obs_input.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

oda_sta_input.o : $(FDDA)/oda_sta_input.f90 mem_oda.o obs_input.o ModDateUtils.o \
	mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

read_ralph.o : $(FDDA)/read_ralph.f90 obs_input.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cu_read.o : $(CUPARM)/cu_read.f90 isan_coms.o mem_cuparm.o ModDateUtils.o \
	mem_grid.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNamelistFile.o : $(INIT)/ModNamelistFile.f90 dump.o modPrintInitial.o \
	ModParallelEnvironment.o parlibf.o grid_dims.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModControlVars.o : $(INIT)/ModControlVars.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParallelEnvironment.o : $(MPI)/ModParallelEnvironment.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOutputUtils.o : $(IO)/ModOutputUtils.f90 dump.o ModTurbFields.o var_tables.o \
	ModBasicFields.o ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

alloc.o : $(MEMORY)/alloc.F90 mem_globrad.o mem_tend.o mem_cuparm.o \
	mem_radiate.o mem_leaf.o mem_aerad.o mem_globaer.o mem_shcu.o mem_gaspart.o \
	mem_scratch.o mem_jules.o mem_micro.o mem_grid.o mem_varinit.o mem_oda.o \
	mem_scalar.o mem_teb_common.o teb_spm_start.o var_tables.o mem_teb.o \
	mem_opt_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMemAlloc.o : $(MEMORY)/ModMemAlloc.F90 mem_globrad.o chem1aq_list.o \
	mem_stilt.o mem_tend.o mem_emiss.o mem_cuparm.o io_params.o mem_radiate.o \
	ModBasicFields.o ModCuParGrell3.o mem_turb_scalar.o node_mod.o mem_chem1.o \
	digitalFilter.o parrrsw.o mem_leaf.o mem_micro_optij.o ModOptical.o mem_aerad.o \
	mem_aer1.o mem_grell.o mem_scratch3_grell.o machine_arq.o mem_chemic.o \
	mem_globaer.o micphys.o carma_fastjx.o mem_shcu.o ModParallelEnvironment.o \
	mem_chem1aq.o shcu_vars_const.o mem_grid_dim_defs.o mem_grell_param2.o \
	mem_gaspart.o mem_carma.o mem_scratch.o mem_volc_chem1.o mem_teb_vars_const.o \
	mem_jules.o chem_dry_dep.o ModGrid.o mem_grid.o grid_dims.o mem_micro.o \
	mem_varinit.o mem_oda.o mem_scalar.o mem_scratch3_grell_sh.o mem_teb_common.o \
	mem_scratch1_grell.o extra.o mem_plume_chem1.o mem_nestb.o ModTuv2.7.o \
	teb_spm_start.o aer1_list.o modIau.o ModTurbFields.o var_tables.o ccatt_start.o \
	mem_tuv.o mem_teb.o mem_scratch1_brams.o ModEvaluation.o mem_scratch2_grell.o \
	chem1_list.o mem_opt_scratch.o chem_sources.o mem_scratch2_grell_sh.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestep.o : $(MODEL)/ModTimestep.F90 mem_stilt.o ModTimeStamp.o mem_tend.o \
	rad_driv.o ModRConv.o ModMicrophysicsMisc.o mem_emiss.o ModRadvc.o \
	ModBasicFields.o mem_cuparm.o ModMicThompsonDriver.o ModSeaSalt.o \
	ModCuParGrell3.o node_mod.o mem_chem1.o ModAcoust.o mem_radiate.o mem_leaf.o \
	ModNudAnalysis.o ModRtimi.o ModMessageSet.o mem_aer1.o ModLeaf3.o ModRShCuPar.o \
	ModRamsMicrophysics2M.o ModOptical.o ModOdaNudge.o ModOzone.o \
	ChemDryDepDriver.o machine_arq.o micphys.o ModTurbK.o ModGasPart.o \
	shcu_vars_const.o ModRexev.o ModRThrm.o ModChemistryDriver.o rconstants.o \
	mem_scratch.o digitalFilter.o ModGrid.o mem_varinit.o grid_dims.o \
	ModRConvGrellCatt.o mem_grid.o mem_oda.o mem_scalar.o ModRbnd.o \
	ModMicrophysicsDrive.o ModCoriolis.o ModMicGfdlDriver.o mem_plume_chem1.o \
	ModMatrixDriver.o teb_spm_start.o ccatt_start.o ModWindFarm.o ModDiffuse.o \
	ModMonotonicAdvection.o sfclyr_jules.o ModUrbanCanopy.o chem_sources.o \
	ChemSourcesDriver.o ModRstilt.o $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestepRK.o : $(MODEL)/ModTimestepRK.F90 mem_stilt.o ModTimeStamp.o \
	mem_tend.o utilsMod.o ModRConv.o ModMicrophysicsMisc.o mem_emiss.o ModRadvc.o \
	mem_cuparm.o mem_radiate.o ModMicThompsonDriver.o ModSeaSalt.o ModCuParGrell3.o \
	node_mod.o mem_chem1.o ModAcoust.o digitalFilter.o mem_leaf.o ModNudAnalysis.o \
	ModRtimi.o ModMessageSet.o mem_aer1.o ModLeaf3.o ModRShCuPar.o \
	ModRamsMicrophysics2M.o ModOptical.o ModOdaNudge.o ModOzone.o \
	ChemDryDepDriver.o machine_arq.o micphys.o ModTurbK.o ModParallelEnvironment.o \
	ModTimestep.o ModGasPart.o shcu_vars_const.o ModRexev.o ModRThrm.o \
	ModChemistryDriver.o rconstants.o mem_scratch.o ModGrid.o ModAerClim.o \
	mem_varinit.o grid_dims.o mem_grid.o mem_oda.o mem_scalar.o ModRbnd.o \
	ModMicrophysicsDrive.o ModCoriolis.o ModMicGfdlDriver.o mem_plume_chem1.o \
	ModMatrixDriver.o modIau.o teb_spm_start.o ModRadvcRK.o ccatt_start.o \
	ModWindFarm.o ModLeaf3OceanOnly.o ModDiffuse.o ModMonotonicAdvection.o \
	sfclyr_jules.o ModUrbanCanopy.o rad_driv.o chem_sources.o ChemSourcesDriver.o \
	ModRstilt.o $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOneProc.o : $(MODEL)/ModOneProc.F90 ModInitHis.o mem_emiss.o domain_decomp.o \
	ModRamsMicrophysics2M.o ModNestGeoSst.o ModSched.o ModParallelEnvironment.o \
	shcu_vars_const.o mem_carma.o soilMoisture.o mem_volc_chem1.o mem_micro.o \
	mem_grid.o ModMkSfcDriver.o extra.o ReadBcst.o teb_spm_start.o dump.o \
	var_tables.o ccatt_start.o tuvParameter.o chem1_list.o mem_globrad.o \
	ModBasicFields.o ModCuParGrell3.o digitalFilter.o ModMkSfcTop.o ModRinit.o \
	ModTimestep.o ModTuvDriver2.7.o parlibf.o mem_chem1aq.o mem_grell_param2.o \
	mem_scratch.o mem_oda.o mem_teb_common.o ModCoriolis.o mem_plume_chem1.o \
	ModTuv2.7.o modIau.o ModVarfFile.o ModEvaluation.o ModMonotonicAdvection.o \
	ModNamelistFile.o ModRhhi.o local_proc.o mem_stilt.o ModTimeStamp.o \
	ModMicrophysicsMisc.o mem_cuparm.o mem_radiate.o mem_chem1.o node_mod.o \
	ModMkSfcFuso.o ModNdviRead.o mem_aer1.o machine_arq.o ModRio.o ModGasPart.o \
	ModRThrm.o ModChemistryDriver.o mem_teb_vars_const.o ModSstRead.o ModGrid.o \
	ModAerClim.o ModDomainDecomp.o ModWindFarm.o isan_coms.o chem_sources.o \
	ModInitMicThompson.o io_params.o mem_leaf.o ModMkSfcSfc.o ModRecycle.o \
	micphys.o ModTimestepRK.o ModParaInit.o ModChemAsgen.o mem_gaspart.o \
	chem_dry_dep.o mem_varinit.o grid_dims.o mem_scalar.o dam.o ModRanlavg.o \
	ModMemAlloc.o aer1_list.o ModGridTree.o meteogram.o ModPostGridNetCDF.o \
	mem_teb.o ref_sounding.o ModUrbanCanopy.o ModOpspec.o memSoilMoisture.o \
	ModPostProcess.o $(UTILS_INCS)/tsNames.h $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsReadHeader.o : $(IO)/ModRamsReadHeader.f90 an_header.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rams_grid.o : $(INIT)/rams_grid.f90 node_mod.o dump.o rconstants.o mem_grid.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

gridset.o : $(INIT)/gridset.f90 rconstants.o mem_grid.o grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

adap_init.o : $(INIT)/adap_init.f90 mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRio.o : $(IO)/ModRio.f90 ModTurbFields.o var_tables.o an_header.o \
	ModBasicFields.o io_params.o mem_grid.o grid_dims.o ref_sounding.o ModRcio.o \
	node_mod.o parlibf.o mem_chem1.o ModParallelEnvironment.o ModNamelistFile.o \
	mem_aerad.o mpi_io_engine-5d.o ModDateUtils.o ReadBcst.o ModControlVars.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h $(UTILS_INCS)/interface.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mpass_dtl.o : $(MPI)/mpass_dtl.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rnest_par.o : $(MPI)/rnest_par.f90 mem_grid.o grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

domain_decomp.o : $(INIT)/domain_decomp.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParaInit.o : $(MPI)/ModParaInit.f90 dump.o var_tables.o mem_grid.o \
	grid_dims.o node_mod.o ModScalarTable.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rnode.o : $(MODEL)/rnode.f90 var_tables.o mem_grid.o grid_dims.o parlibf.o \
	node_mod.o mem_leaf.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

nud_read.o : $(FDDA)/nud_read.f90 mem_varinit.o mem_grid.o mem_chem1.o \
	isan_coms.o ModNudAnalysis.o ModDateUtils.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

nud_update.o : $(FDDA)/nud_update.f90 var_tables.o ModInitHis.o an_header.o \
	grid_struct.o mem_varinit.o mem_grid.o ModRcio.o mem_chem1.o chem1_list.o \
	mem_aerad.o rconstants.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cond_read.o : $(FDDA)/cond_read.f90 mem_varinit.o mem_grid.o isan_coms.o \
	ModNudAnalysis.o ModDateUtils.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cond_update.o : $(FDDA)/cond_update.f90 var_tables.o ModInitHis.o an_header.o \
	grid_struct.o mem_grid.o mem_varinit.o ModRcio.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudAnalysis.o : $(FDDA)/ModNudAnalysis.f90 dump.o mem_tend.o ModBasicFields.o \
	mem_scratch.o mem_varinit.o mem_grid.o mem_chem1.o ModEvaluation.o node_mod.o \
	chem1_list.o modIau.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCoriolis.o : $(MODEL)/ModCoriolis.f90 mem_tend.o ModBasicFields.o \
	ref_sounding.o mem_grid.o parlibf.o ModBuffering.o node_mod.o rconstants.o \
	mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGeodat.o : $(MKSFC)/ModGeodat.f90 teb_spm_start.o io_params.o mem_grid.o \
	mem_leaf.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLanduseInput.o : $(MKSFC)/ModLanduseInput.f90 mem_mksfc.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3.o : $(SURFACE)/ModLeaf3.f90 teb_spm_start.o ModTurbFields.o micphys.o \
	ccatt_start.o mem_teb.o ModBasicFields.o mem_cuparm.o mem_micro.o io_params.o \
	mem_grid.o mem_radiate.o node_mod.o leaf_coms.o mem_leaf.o mem_teb_common.o \
	rconstants.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3OceanOnly.o : $(SURFACE)/ModLeaf3OceanOnly.f90 ModTurbFields.o micphys.o \
	ccatt_start.o mem_cuparm.o ModBasicFields.o mem_micro.o mem_grid.o \
	mem_radiate.o io_params.o node_mod.o ModCuParGrell3.o mem_leaf.o leaf_coms.o \
	ModLeaf3.o ConvPar_GF_GEOS5.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

leaf3_hyd.o : $(SURFACE)/leaf3_hyd.f90 mem_grid.o mem_leaf.o leaf_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

leaf3_init.o : $(SURFACE)/leaf3_init.f90 teb_spm_start.o io_params.o mem_grid.o \
	grid_dims.o mem_leaf.o leaf_coms.o ModLeaf3.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

leaf3_teb.o : $(SURFACE)/leaf3_teb.f90 mem_emiss.o ModGasPart.o \
	mem_teb_vars_const.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGasPart.o : $(TEB_SPM)/ModGasPart.f90 var_tables.o mem_teb_vars_const.o \
	an_header.o mem_emiss.o ModBasicFields.o mem_grid.o grid_dims.o ModRcio.o \
	parlibf.o node_mod.o mem_leaf.o mem_gaspart.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOzone.o : $(TEB_SPM)/ModOzone.f90 mem_radiate.o ModBasicFields.o mem_grid.o \
	mem_gaspart.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSstRead.o : $(MKSFC)/ModSstRead.f90 mem_mksfc.o io_params.o mem_grid.o \
	grid_dims.o node_mod.o mem_leaf.o ModMkSfcSst.o ModDateUtils.o ReadBcst.o \
	ModControlVars.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNdviRead.o : $(MKSFC)/ModNdviRead.f90 mem_mksfc.o io_params.o mem_grid.o \
	grid_dims.o node_mod.o mem_leaf.o ModMkSfcNdvi.o ModDateUtils.o ReadBcst.o \
	ModControlVars.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcDriver.o : $(MKSFC)/ModMkSfcDriver.f90 teb_spm_start.o mem_mksfc.o \
	ModNestGeoSst.o ModSstRead.o io_params.o mem_grid.o grid_dims.o ModMkSfcTop.o \
	node_mod.o ModMkSfcFuso.o ModNdviRead.o ModMkSfcNdvi.o ModMkSfcSst.o \
	ModMkSfcSfc.o ReadBcst.o ModLanduseInput.o ModControlVars.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSfc.o : $(MKSFC)/ModMkSfcSfc.f90 dump.o mem_mksfc.o io_params.o \
	mem_grid.o node_mod.o mem_leaf.o ReadBcst.o ModControlVars.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcTop.o : $(MKSFC)/ModMkSfcTop.f90 dump.o mem_mksfc.o io_params.o \
	mem_grid.o node_mod.o ReadBcst.o ModControlVars.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSst.o : $(MKSFC)/ModMkSfcSst.f90 mem_mksfc.o io_params.o mem_grid.o \
	grid_dims.o mem_leaf.o ModGeodat.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcNdvi.o : $(MKSFC)/ModMkSfcNdvi.f90 dump.o mem_mksfc.o io_params.o \
	mem_grid.o mem_leaf.o ModLanduseInput.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcFuso.o : $(MKSFC)/ModMkSfcFuso.f90 mem_teb_vars_const.o mem_mksfc.o \
	mem_teb.o mem_emiss.o io_params.o mem_grid.o node_mod.o mem_gaspart.o \
	ReadBcst.o ModControlVars.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mpass_full.o : $(MPI)/mpass_full.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestGeoSst.o : $(MKSFC)/ModNestGeoSst.f90 dump.o ModTurbFields.o \
	ModLanduseInput.o ModInitHis.o mem_mksfc.o ccatt_start.o ModBasicFields.o \
	io_params.o mem_grid.o grid_dims.o ModMkSfcTop.o node_mod.o mem_leaf.o \
	soilMoisture.o ModGeodat.o memSoilMoisture.o mem_scratch.o ModControlVars.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

nest_filldens.o : $(NESTING)/nest_filldens.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

nest_intrp.o : $(NESTING)/nest_intrp.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

nest_feed.o : $(NESTING)/nest_feed.f90 mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAcoust.o : $(MODEL)/ModAcoust.f90 raco_adap.o mem_tend.o micphys.o ModGrid.o \
	ModBasicFields.o ref_sounding.o mem_grid.o ModParallelEnvironment.o node_mod.o \
	ModMessageSet.o rconstants.o mem_scratch.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mod_GhostBlockPartition.o : $(MODEL)/mod_GhostBlockPartition.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mod_GhostBlock.o : $(MODEL)/mod_GhostBlock.f90 mod_GhostBlockPartition.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvc.o : $(MODEL)/ModRadvc.f90 mem_tend.o chem_dry_dep.o ccatt_start.o \
	ModBasicFields.o ModParallelEnvironment.o mem_grid.o grid_dims.o \
	ModMonotonicAdvection.o mem_chem1.o ModNamelistFile.o mem_aer1.o \
	ModScalarTable.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcRK.o : $(MODEL)/ModRadvcRK.f90 mem_stilt.o mem_tend.o ModGrid.o \
	mem_grid.o grid_dims.o mem_chem1.o ModParallelEnvironment.o node_mod.o \
	ModMessageSet.o ModRexev.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRanlavg.o : $(IO)/ModRanlavg.f90 var_tables.o ModBasicFields.o io_params.o \
	mem_grid.o grid_dims.o ModRThrm.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/interface.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRbnd.o : $(BC)/ModRbnd.f90 ModTurbKE.o ModTurbFields.o mem_tend.o micphys.o \
	ccatt_start.o ModMicrophysicsMisc.o ModBasicFields.o ref_sounding.o mem_grid.o \
	mem_chem1.o node_mod.o ModScalarTable.o mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConv.o : $(CUPARM)/ModRConv.f90 mem_tend.o conv_coms.o mem_cuparm.o \
	ModBasicFields.o mem_grid.o node_mod.o rconstants.o mem_scratch.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRinit.o : $(INIT)/ModRinit.f90 ModTurbKE.o ModTurbFields.o micphys.o \
	ModBasicFields.o ref_sounding.o mem_varinit.o mem_grid.o io_params.o \
	mem_micro.o node_mod.o ModRbnd.o rconstants.o mem_scratch.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRtimi.o : $(MODEL)/ModRtimi.f90 mem_tend.o mem_cuparm.o ModBasicFields.o \
	mem_grid.o node_mod.o shcu_vars_const.o ModScalarTable.o mem_grell.o \
	mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ruser.o : $(SURFACE)/ruser.f90 ccatt_start.o io_params.o mem_grid.o node_mod.o \
	leaf_coms.o mem_leaf.o memSoilMoisture.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

vtab_fill.o : $(MEMORY)/vtab_fill.f90 aer1_list.o var_tables.o io_params.o \
	mem_chem1.o chem1_list.o mem_aer1.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

extra.o : $(MEMORY)/extra.f90 dump.o ModNamelistFile.o var_tables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

aer1_list.o : $(AEROSOL)/aer1_list_$(AERLEVEL).f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)
	@ln -fs aer1_list_$(AERLEVEL).o aer1_list.o

mem_aer1.o : $(CCATT)/mem_aer1.f90 aer1_list.o dump.o var_tables.o io_params.o \
	mem_grid.o grid_dims.o mem_chem1.o node_mod.o ModNamelistFile.o \
	ModScalarTable.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem1_list.o : $(MODEL_CHEM)/chem1_list.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem1aq_list.o : $(MODEL_CHEM)/chem1aq_list.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1aq.o : $(CCATT)/mem_chem1aq.f90 chem1aq_list.o var_tables.o grid_dims.o \
	mem_chem1.o ModNamelistFile.o ModScalarTable.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1.o : $(CCATT)/mem_chem1.f90 var_tables.o io_params.o grid_dims.o \
	chem1_list.o ModNamelistFile.o ModScalarTable.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_plume_chem1.o : $(CCATT)/mem_plume_chem1.f90 chem1_list.o ModNamelistFile.o \
	var_tables.o mem_chem1.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chemic.o : $(CCATT)/mem_chemic.f90 micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_volc_chem1.o : $(CCATT)/mem_volc_chem1.f90 ModNamelistFile.o var_tables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_sources.o : $(CCATT)/chem_sources.f90 aer1_list.o mem_volc_chem1.o \
	io_params.o mem_grid.o mem_chem1.o parlibf.o ModNamelistFile.o mem_aer1.o \
	mem_plume_chem1.o ModDateUtils.o ReadBcst.o ModControlVars.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_plumerise_scalar.o : $(CCATT)/chem_plumerise_scalar.f90 node_mod.o \
	mem_plume_chem1.o mem_chem1.o mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemSourcesDriver.o : $(CCATT)/ChemSourcesDriver.f90 aer1_list.o \
	mem_volc_chem1.o mem_stilt.o chem_plumerise_scalar.o io_params.o mem_chem1.o \
	chem1_list.o mem_leaf.o mem_aer1.o chem_sources.o mem_plume_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_dry_dep.o : $(MODEL_CHEM)/chem_dry_dep.f90 aer1_list.o mem_chem1.o \
	chem1_list.o mem_aer1.o extra.o ModDateUtils.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemDryDepDriver.o : $(MODEL_CHEM)/ChemDryDepDriver.f90 ModTurbFields.o \
	chem_dry_dep.o micphys.o mem_cuparm.o ModBasicFields.o mem_grid.o grid_dims.o \
	mem_micro.o mem_radiate.o mem_chem1.o mem_leaf.o mem_aer1.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbK.o : $(TURB)/ModTurbK.f90 mem_stilt.o ke_coms.o mem_tend.o mem_cuparm.o \
	ModTurbDiff.o ModBasicFields.o mem_chem1.o mem_turb_scalar.o node_mod.o \
	mem_leaf.o ModScalarTable.o mem_grell.o ModTurbKE.o ModRGrad.o micphys.o \
	ModTurbDiffAdap.o ModTurbKAdap.o rconstants.o mem_scratch.o mem_grid.o \
	grid_dims.o mem_micro.o ModTurbFields.o ccatt_start.o ModMonotonicAdvection.o \
	ModNamelistFile.o ModRstilt.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffuse.o : $(TURB)/ModDiffuse.f90 ModTurbKE.o ModTurbFields.o ke_coms.o \
	mem_tend.o micphys.o ModTurbK.o ModBasicFields.o ModTurbDiff.o mem_grid.o \
	mem_micro.o node_mod.o ModNamelistFile.o mem_leaf.o mem_opt_scratch.o \
	ModScalarTable.o ModDiffSclr.o mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKE.o : $(TURB)/ModTurbKE.f90 ModTurbFields.o ke_coms.o mem_grid.o \
	rconstants.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiff.o : $(TURB)/ModTurbDiff.f90 ModRGrad.o mem_cuparm.o mem_grid.o \
	mem_opt_scratch.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffSclr.o : $(TURB)/ModDiffSclr.f90 ModTurbDiff.o mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

urban.o : $(SURFACE)/urban.f90 mem_teb_vars_const.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModUrbanCanopy.o : $(SURFACE)/ModUrbanCanopy.f90 ModTurbFields.o mem_tend.o \
	ModBasicFields.o mem_grid.o node_mod.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsDrive.o : $(MICRO)/ModMicrophysicsDrive.f90 micphys.o \
	mem_chemic.o ModMicrophysicsMisc.o ModBasicFields.o mem_radiate.o mem_micro.o \
	mem_grid.o mem_chem1.o node_mod.o mem_chem1aq.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsMicrophysics2M.o : $(MICRO)/ModRamsMicrophysics2M.f90 dump.o micphys.o \
	ModBasicFields.o mem_micro.o mem_grid.o grid_dims.o node_mod.o mem_leaf.o \
	rconstants.o mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSched.o : $(MODEL)/ModSched.f90 dump.o local_proc.o ModBasicFields.o \
	ref_sounding.o mem_grid.o parlibf.o mem_varinit.o isan_coms.o node_mod.o \
	io_params.o ModNamelistFile.o mem_radiate.o mem_cuparm.o shcu_vars_const.o \
	ReadBcst.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

local_proc.o : $(MODEL)/local_proc.F90 dump.o mem_stilt.o io_params.o mem_grid.o \
	grid_dims.o ref_sounding.o node_mod.o rconstants.o ReadBcst.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mic_init.o : $(MICRO)/mic_init.f90 dump.o rconstants.o mem_grid.o micphys.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsMisc.o : $(MICRO)/ModMicrophysicsMisc.f90 micphys.o \
	ModBasicFields.o mem_micro.o mem_grid.o rconstants.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mic_nuc.o : $(MICRO)/mic_nuc.f90 micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mic_vap.o : $(MICRO)/mic_vap.f90 rconstants.o micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mic_gamma.o : $(MICRO)/mic_gamma.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mic_tabs.o : $(MICRO)/mic_tabs.f90 micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mic_coll.o : $(MICRO)/mic_coll.f90 micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_carma.o : $(RADIATE)/rad_carma.F90 mem_aerad.o aer1_list.o mem_globrad.o \
	machine_arq.o mem_globaer.o ccatt_start.o mem_tuv.o carma_fastjx.o \
	mem_radiate.o mem_grid.o grid_dims.o mem_chem1.o chem1_list.o node_mod.o \
	mem_leaf.o mem_aer1.o ModDateUtils.o mem_carma.o rconstants.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOpspec.o : $(IO)/ModOpspec.f90 mem_globrad.o chem1aq_list.o mem_stilt.o \
	mem_emiss.o mem_radiate.o io_params.o mem_cuparm.o mem_chem1.o mem_leaf.o \
	mem_aer1.o micphys.o mem_chem1aq.o shcu_vars_const.o mem_grell_param2.o \
	mem_grid.o grid_dims.o mem_varinit.o modIau.o teb_spm_start.o aer1_list.o \
	ccatt_start.o chem1_list.o ModNamelistFile.o chem_sources.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_driv.o : $(RADIATE)/rad_driv.f90 mem_radiate.o ModCarmaDriver.o \
	ModBasicFields.o ModRrtmDriver.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCarmaDriver.o : $(RADIATE)/ModCarmaDriver.f90 teb_spm_start.o mem_tend.o \
	micphys.o mem_radiate.o rad_carma.o mem_micro.o ModBasicFields.o mem_cuparm.o \
	grid_dims.o node_mod.o mem_grid.o mem_leaf.o ModLeaf3.o mem_teb_common.o \
	mem_scratch1_grell.o ModDateUtils.o mem_carma.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRrtmDriver.o : $(RADIATE)/ModRrtmDriver.f90 mem_tend.o mem_cuparm.o \
	ModBasicFields.o mem_radiate.o mem_chem1.o node_mod.o parrrsw.o mem_leaf.o \
	parrrtm.o ModOptical.o rrtmg_lw_rad.o rrtmg_lw_cldprop.o micphys.o leaf_coms.o \
	mem_grell_param2.o ModDateUtils.o rconstants.o mem_carma.o rrtmg_sw_rad.o \
	rrtmg_sw_cldprop.o mem_rrtm.o parkind.o mem_grid.o grid_dims.o mem_micro.o \
	mcica_subcol_gen_sw.o mem_scratch1_grell.o mcica_subcol_gen_lw.o \
	teb_spm_start.o ccatt_start.o mem_tuv.o ref_sounding.o \
	$(UTILS_INCS)/aerosol_setup.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRShCuPar.o : $(CUPARM)/ModRShCuPar.f90 ModTurbFields.o mem_tend.o ModRConv.o \
	conv_coms.o ModBasicFields.o mem_shcu.o mem_grid.o mem_micro.o node_mod.o \
	shcu_vars_const.o mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConvGrellCatt.o : $(CUPARM)/ModRConvGrellCatt.f90 mem_stilt.o mem_tend.o \
	mem_radiate.o io_params.o mem_cuparm.o ModCuParGrell3.o node_mod.o mem_leaf.o \
	mem_grell.o micphys.o mem_grell_param2.o rconstants.o mem_scratch.o ModGrid.o \
	mem_grid.o mem_micro.o mem_scalar.o mem_scratch1_grell.o ModChemConvTransp.o \
	ccatt_start.o ModRstilt.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_cu_g3.o : $(CUPARM)/module_cu_g3.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_cu_gf.o : $(CUPARM)/module_cu_gf.f90 node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_cu_gf_v5.1.o : $(CUPARM)/module_cu_gf_v5.1.f90 module_gate.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_cu_gd_fim.o : $(CUPARM)/module_cu_gd_fim.f90 module_gate.o Phys_const.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

MAPL_Constants.o : $(CUPARM)/MAPL_Constants.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

Henrys_Law_cts.o : $(CUPARM)/Henrys_Law_cts.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ConvPar_GF_GEOS5.o : $(CUPARM)/ConvPar_GF_GEOS5.F90 module_gate.o \
	MAPL_Constants.o Henrys_Law_cts.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuParGrell3.o : $(CUPARM)/ModCuParGrell3.F90 mem_stilt.o mem_tend.o \
	mem_cuparm.o io_params.o mem_radiate.o mem_chem1.o ModBasicFields.o ModRadvc.o \
	node_mod.o mem_leaf.o ModMessageSet.o mem_grell.o micphys.o mem_grell_param2.o \
	rconstants.o mem_carma.o mem_scratch.o mem_jules.o ModGrid.o mem_grid.o \
	mem_varinit.o mem_micro.o module_cu_gf.o grid_dims.o mem_scratch1_grell.o \
	ConvPar_GF_GEOS5.o ModChemConvTransp.o Phys_const.o var_tables.o ccatt_start.o \
	module_cu_gf_v5.1.o module_cu_g3.o ModNamelistFile.o ModRstilt.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

kbcon_ecmwf.o : $(CUPARM)/kbcon_ecmwf.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cup_grell_catt_deep.o : $(CUPARM)/cup_grell_catt_deep.f90 Phys_const.o \
	mem_scratch3_grell.o kbcon_ecmwf.o ccatt_start.o mem_grid.o mem_varinit.o \
	mem_scratch2_grell.o node_mod.o mem_grell_param2.o mem_carma.o \
	cup_output_vars.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemConvTransp.o : $(CCATT)/ModChemConvTransp.f90 Phys_const.o aer1_list.o \
	mem_cuparm.o mem_grid.o mem_chem1.o node_mod.o chem1_list.o mem_tconv.o \
	mem_aer1.o mem_scratch1_grell.o mem_grell_param2.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cup_grell_catt_shallow.o : $(CUPARM)/cup_grell_catt_shallow.f90 Phys_const.o \
	mem_grid.o mem_varinit.o mem_scratch3_grell_sh.o node_mod.o mem_grell_param2.o \
	cup_output_vars.o mem_scratch2_grell_sh.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cup_env.o : $(CUPARM)/cup_env.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cup_env_catt.o : $(CUPARM)/cup_env_catt.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

upcase.o : $(CUPARM)/upcase.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

memSoilMoisture.o : $(SOIL_MOISTURE)/memSoilMoisture.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

soilMoisture.o : $(SOIL_MOISTURE)/soilMoisture.F90 ModTurbFields.o \
	ModBasicFields.o io_params.o mem_grid.o parlibf.o node_mod.o ModNamelistFile.o \
	leaf_coms.o mem_leaf.o mem_aerad.o memSoilMoisture.o ReadBcst.o rconstants.o \
	ModControlVars.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRcio.o : $(IO)/ModRcio.f90 mem_stilt.o an_header.o micphys.o mem_cuparm.o \
	io_params.o ModParallelEnvironment.o grid_dims.o mem_grid.o mem_radiate.o \
	ref_sounding.o ModNamelistFile.o mem_leaf.o leaf_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRGrad.o : $(TURB)/ModRGrad.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRhhi.o : $(INIT)/ModRhhi.f90 micphys.o ModBasicFields.o ref_sounding.o \
	mem_grid.o ModRinit.o rconstants.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitHis.o : $(IO)/ModInitHis.f90 var_tables.o micphys.o an_header.o \
	ModBasicFields.o ref_sounding.o mem_grid.o io_params.o ModRcio.o mem_chem1.o \
	mem_varinit.o chem1_list.o ModRamsReadHeader.o mem_leaf.o ModRinit.o \
	mem_aerad.o leaf_coms.o ModLeaf3.o rconstants.o mem_scratch.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRecycle.o : $(IO)/ModRecycle.f90 aer1_list.o dump.o var_tables.o an_header.o \
	ModDateUtils.o io_params.o mem_grid.o mem_chem1.o node_mod.o chem1_list.o \
	ModRamsReadHeader.o mem_aerad.o mem_aer1.o ReadBcst.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRThrm.o : $(MODEL)/ModRThrm.f90 micphys.o ModBasicFields.o mem_grid.o \
	mem_micro.o rconstants.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

varf_update.o : $(FDDA)/varf_update.f90 ModInitHis.o ref_sounding.o mem_grid.o \
	rconstants.o mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

radvc_adap.o : $(MODEL)/radvc_adap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKAdap.o : $(TURB)/ModTurbKAdap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiffAdap.o : $(TURB)/ModTurbDiffAdap.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

raco_adap.o : $(MODEL)/raco_adap.f90 ModGrid.o mem_grid.o node_mod.o ModRbnd.o \
	ModMessageSet.o rconstants.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_teb.o : $(TEB_SPM)/mem_teb.f90 var_tables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_teb_common.o : $(TEB_SPM)/mem_teb_common.f90 var_tables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_teb_vars_const.o : $(TEB_SPM)/mem_teb_vars_const.f90 ModNamelistFile.o \
	grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_gaspart.o : $(TEB_SPM)/mem_gaspart.f90 mem_emiss.o var_tables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_emiss.o : $(TEB_SPM)/mem_emiss.f90 ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

aobj.o : $(ISAN)/aobj.f90 isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

asgen.o : $(ISAN)/asgen.f90 isan_coms.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

asti2.o : $(ISAN)/asti2.f90 isan_coms.o ModDateUtils.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

asti.o : $(ISAN)/asti.f90 isan_coms.o rconstants.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

astp.o : $(ISAN)/astp.f90 isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

avarf.o : $(ISAN)/avarf.f90 isan_coms.o ModRbnd.o rconstants.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

file_inv.o : $(ISAN)/file_inv.f90 isan_coms.o mem_grid.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

first_rams.o : $(ISAN)/first_rams.f90 an_header.o mem_grid.o grid_dims.o \
	ModRcio.o isan_coms.o rconstants.o mem_scratch.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

isan_io.o : $(ISAN)/isan_io.f90 isan_coms.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

refstate.o : $(ISAN)/refstate.f90 rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

v_interps.o : $(ISAN)/v_interps.f90 isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_isan_coms.o : $(ISAN_CHEM)/chem_isan_coms.f90 isan_coms.o chem1_list.o \
	aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_aobj.o : $(ISAN_CHEM)/chem_aobj.f90 isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsgen.o : $(ISAN_CHEM)/ModChemAsgen.F90 aer1_list.o dump.o \
	chem_isan_coms.o io_params.o mem_grid.o grid_dims.o mem_chem1.o isan_coms.o \
	ModMkSfcTop.o chem1_list.o node_mod.o mem_aer1.o ModDateUtils.o \
	ModControlVars.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_asti2.o : $(ISAN_CHEM)/chem_asti2.f90 isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_asti.o : $(ISAN_CHEM)/chem_asti.f90 chem_isan_coms.o mem_grid.o mem_chem1.o \
	isan_coms.o mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_astp.o : $(ISAN_CHEM)/chem_astp.F90 dump.o chem_isan_coms.o ModDateUtils.o \
	mem_varinit.o mem_chem1.o isan_coms.o chem1_list.o rconstants.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_avarf.o : $(ISAN_CHEM)/chem_avarf.f90 chem_isan_coms.o mem_grid.o \
	mem_chem1.o isan_coms.o ModRbnd.o mem_aer1.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_file_inv.o : $(ISAN_CHEM)/chem_file_inv.f90 isan_coms.o dump.o \
	ModDateUtils.o mem_grid.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_first_rams.o : $(ISAN_CHEM)/chem_first_rams.f90 an_header.o mem_grid.o \
	grid_dims.o ModRcio.o isan_coms.o rconstants.o mem_scratch.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_isan_io.o : $(ISAN_CHEM)/chem_isan_io.f90 isan_coms.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_refstate.o : $(ISAN_CHEM)/chem_refstate.f90 ccatt_start.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_v_interps.o : $(ISAN_CHEM)/chem_v_interps.f90 isan_coms.o dump.o \
	rconstants.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

carma_fastjx.o : $(CCATT)/carma_fastjx.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_fastjx_data.o : $(CCATT)/chem_fastjx_data.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_fastjx57.o : $(CCATT)/chem_fastjx57.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_fastjx_driv.o : $(CCATT)/chem_fastjx_driv.f90 chem1_list.o rconstants.o \
	chem_fastjx57.o chem_fastjx_data.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_utils.o : $(CCATT)/chem_spack_utils.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_spack.o : $(CCATT)/mem_spack.f90 chem_spack_utils.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_uv_att.o : $(CCATT)/chem_uv_att.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_solve_sparse.o : $(CCATT)/chem_spack_solve_sparse.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_lu.o : $(CCATT)/chem_spack_lu.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_jacdchemdc.o : $(MODEL_CHEM)/chem_spack_jacdchemdc.f90 \
	chem_spack_dratedc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_dratedc.o : $(MODEL_CHEM)/chem_spack_dratedc.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_kinetic.o : $(MODEL_CHEM)/chem_spack_kinetic.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_fexprod.o : $(MODEL_CHEM)/chem_spack_fexprod.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_fexloss.o : $(MODEL_CHEM)/chem_spack_fexloss.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rates.o : $(MODEL_CHEM)/chem_spack_rates.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_fexchem.o : $(MODEL_CHEM)/chem_spack_fexchem.f90 chem_spack_rates.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros.o : $(CCATT)/chem_spack_ros.f90 chem_spack_solve_sparse.o \
	mem_chem1.o mem_spack.o chem_spack_fexchem.o chem_spack_jacdchemdc.o \
	chem_spack_kinetic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros_dyndt.o : $(CCATT)/chem_spack_ros_dyndt.f90 \
	chem_spack_solve_sparse.o chem_spack_ros.o mem_chem1.o mem_spack.o \
	chem_spack_fexchem.o chem_spack_jacdchemdc.o chem_spack_kinetic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rodas3_dyndt.o : $(CCATT)/chem_spack_rodas3_dyndt.f90 extra.o \
	chem_spack_ros.o mem_grid.o mem_chem1.o mem_spack.o chem_spack_fexchem.o \
	chem_spack_jacdchemdc.o chem_spack_kinetic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_qssa.o : $(CCATT)/chem_spack_qssa.f90 chem_spack_fexprod.o \
	chem_spack_rates.o mem_chem1.o chem_spack_dratedc.o chem_spack_kinetic.o \
	chem_spack_fexloss.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_trans_gasaq.o : $(MODEL_CHEM)/chem_trans_gasaq.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_trans_liq.o : $(CCATT)/chem_trans_liq.f90 mem_chem1aq.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_orage.o : $(CCATT)/chem_orage.f90 mem_scratch1_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemistryDriver.o : $(CCATT)/ModChemistryDriver.f90 mem_globrad.o \
	chem_spack_solve_sparse.o mem_stilt.o chem1aq_list.o chem_spack_utils.o \
	ModBasicFields.o mem_radiate.o mem_cuparm.o mem_chem1.o node_mod.o \
	chem_spack_qssa.o mem_aerad.o mem_aer1.o parrrtm.o chem_trans_liq.o micphys.o \
	mem_chemic.o carma_fastjx.o chem_spack_ros_dyndt.o ModTuvDriver2.7.o \
	mem_spack.o mem_chem1aq.o chem_spack_rodas3_dyndt.o mem_grell_param2.o \
	rconstants.o mem_carma.o mem_scratch.o chem_spack_ros.o chem_trans_gasaq.o \
	mem_rrtm.o mem_grid.o mem_micro.o grid_dims.o chem_uv_att.o \
	mem_scratch1_grell.o extra.o chem_orage.o aer1_list.o chem_fastjx_driv.o \
	chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimeStamp.o : $(MODEL)/ModTimeStamp.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBramsGrid.o : $(POST_SRC)/ModBramsGrid.f90 ModPostUtils.o ref_sounding.o \
	ModParallelEnvironment.o mem_grid.o node_mod.o ModNamelistFile.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostGrid.o : $(POST_SRC)/ModPostGrid.F90 ModTurbFields.o ModPostUtils.o \
	ModBasicFields.o io_params.o ModParallelEnvironment.o mem_grid.o parlibf.o \
	ModNamelistFile.o ModBramsGrid.o ModOutputUtils.o ModPostOneFieldNetCDF.o \
	ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostGridNetCDF.o : $(POST_SRC)/ModPostGridNetCDF.F90 dump.o ModPostUtils.o \
	io_params.o mem_grid.o ModNamelistFile.o ModBramsGrid.o ModDateUtils.o \
	ModPostTypes.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneFieldNetCDF.o : $(POST_SRC)/ModPostOneFieldNetCDF.F90 \
	ModPostGridNetCDF.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostTypes.o : $(POST_SRC)/ModPostTypes.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneFieldUtils.o : $(POST_SRC)/ModPostOneFieldUtils.f90 ModPostGrid.o \
	ModBramsGrid.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField.o : $(POST_SRC)/ModPostOneField.f90 ModPostOneField2d.o dump.o \
	ModTurbFields.o ModPostUtils.o ModPostOneField7d.o ModBasicFields.o \
	ModPostOneFieldUtils.o ModPostOneField8d.o ModPostOneField3d.o node_mod.o \
	ModNamelistFile.o ModBramsGrid.o ModPostTypes.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField2d.o : $(POST_SRC)/ModPostOneField2d.f90 dump.o ModTurbFields.o \
	ModPostUtils.o micphys.o ModBasicFields.o io_params.o mem_grid.o ModPostGrid.o \
	ModPostOneFieldUtils.o mem_radiate.o node_mod.o mem_cuparm.o ModNamelistFile.o \
	ModBramsGrid.o ModOutputUtils.o mem_aerad.o ModPostTypes.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField3d.o : $(POST_SRC)/ModPostOneField3d.f90 ModTurbFields.o \
	ModPostUtils.o micphys.o ModBasicFields.o mem_grid.o mem_varinit.o \
	ModPostGrid.o ModPostOneFieldUtils.o node_mod.o ModNamelistFile.o \
	ModBramsGrid.o ModOutputUtils.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField7d.o : $(POST_SRC)/ModPostOneField7d.f90 ModTurbFields.o \
	ModPostUtils.o ModBasicFields.o ModPostOneFieldUtils.o ModNamelistFile.o \
	ModBramsGrid.o ModOutputUtils.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField8d.o : $(POST_SRC)/ModPostOneField8d.f90 ModTurbFields.o \
	ModPostUtils.o ModBasicFields.o ModPostOneFieldUtils.o ModNamelistFile.o \
	ModBramsGrid.o ModOutputUtils.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostProcess.o : $(POST_SRC)/ModPostProcess.F90 ModGridTree.o ModTurbFields.o \
	ModPostGridNetCDF.o ModGrid.o ModBasicFields.o io_params.o \
	ModParallelEnvironment.o ModPostGrid.o ModPostOneField.o ModNamelistFile.o \
	ModBramsGrid.o ModMessageSet.o ModPostTypes.o $(UTILS_INCS)/tsNames.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostUtils.o : $(POST_SRC)/ModPostUtils.f90 dump.o ModParallelEnvironment.o \
	mem_leaf.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h \
	$(POST_INCS)/post_rconfig.h $(POST_INCS)/post_rconstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRexev.o : $(STILT)/ModRexev.f90 mem_stilt.o mem_tend.o micphys.o \
	ModBasicFields.o ModRadvc.o mem_grid.o mem_micro.o rconstants.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRstilt.o : $(STILT)/ModRstilt.f90 ModTurbFields.o mem_stilt.o mem_cuparm.o \
	ModBasicFields.o mem_grid.o grid_dims.o ModMonotonicAdvection.o \
	mem_scratch1_grell.o mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

turb_constants.o : $(STILT)/turb_constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

tkenn.o : $(STILT)/tkenn.f90 mem_stilt.o mem_grid.o rconstants.o \
	turb_constants.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

digitalFilter.o : $(MODEL)/digitalFilter.f90 var_tables.o ModBasicFields.o \
	io_params.o mem_grid.o grid_dims.o node_mod.o ModNamelistFile.o ModDateUtils.o \
	ReadBcst.o ModControlVars.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMonotonicAdvection.o : $(MODEL)/ModMonotonicAdvection.f90 chem_dry_dep.o \
	micphys.o ccatt_start.o ModGrid.o ModParallelEnvironment.o mem_grid.o \
	mem_chem1.o ModNamelistFile.o ModMessageSet.o mem_aer1.o ModDomainDecomp.o \
	rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSeaSalt.o : $(CCATT)/ModSeaSalt.f90 aer1_list.o ccatt_start.o ModAerClim.o \
	ModBasicFields.o io_params.o mem_grid.o mem_chem1.o node_mod.o mem_leaf.o \
	mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

meteogramType.o : $(IO)/meteogramType.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

meteogram.o : $(IO)/meteogram.f90 var_tables.o ModPostUtils.o satPolyColision.o \
	mem_grid.o node_mod.o ModNamelistFile.o meteogramType.o $(UTILS_INCS)/files.h \
	$(POST_INCS)/post_rconstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

parkind.o : $(RRTMG_SW_MOD)/parkind.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

parrrsw.o : $(RRTMG_SW_MOD)/parrrsw.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrsw_aer.o : $(RRTMG_SW_MOD)/rrsw_aer.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_cld.o : $(RRTMG_SW_MOD)/rrsw_cld.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_con.o : $(RRTMG_SW_MOD)/rrsw_con.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg16.o : $(RRTMG_SW_MOD)/rrsw_kg16.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg17.o : $(RRTMG_SW_MOD)/rrsw_kg17.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg18.o : $(RRTMG_SW_MOD)/rrsw_kg18.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg19.o : $(RRTMG_SW_MOD)/rrsw_kg19.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg20.o : $(RRTMG_SW_MOD)/rrsw_kg20.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg21.o : $(RRTMG_SW_MOD)/rrsw_kg21.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg22.o : $(RRTMG_SW_MOD)/rrsw_kg22.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg23.o : $(RRTMG_SW_MOD)/rrsw_kg23.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg24.o : $(RRTMG_SW_MOD)/rrsw_kg24.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg25.o : $(RRTMG_SW_MOD)/rrsw_kg25.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg26.o : $(RRTMG_SW_MOD)/rrsw_kg26.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg27.o : $(RRTMG_SW_MOD)/rrsw_kg27.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg28.o : $(RRTMG_SW_MOD)/rrsw_kg28.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg29.o : $(RRTMG_SW_MOD)/rrsw_kg29.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_ref.o : $(RRTMG_SW_MOD)/rrsw_ref.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_tbl.o : $(RRTMG_SW_MOD)/rrsw_tbl.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_vsn.o : $(RRTMG_SW_MOD)/rrsw_vsn.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_wvn.o : $(RRTMG_SW_MOD)/rrsw_wvn.f90 parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_random_numbers.o : $(RRTMG_SW_SRC)/mcica_random_numbers.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_sw.o : $(RRTMG_SW_SRC)/mcica_subcol_gen_sw.f90 rrsw_vsn.o \
	rrsw_wvn.o mcica_random_numbers.o parkind.o parrrsw.o rrsw_con.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprmc.f90 rrsw_vsn.o parkind.o \
	parrrsw.o rrsw_cld.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprop.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprop.f90 rrsw_vsn.o parkind.o \
	parrrsw.o rrsw_cld.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_init.o : $(RRTMG_SW_SRC)/rrtmg_sw_init.f90 rrsw_kg17.o rrsw_kg24.o \
	parrrsw.o rrsw_kg20.o rrsw_vsn.o rrsw_aer.o parkind.o rrsw_kg18.o rrsw_kg29.o \
	rrsw_kg19.o rrsw_kg27.o rrsw_kg28.o rrsw_kg21.o rrsw_tbl.o rrsw_kg25.o \
	rrsw_wvn.o rrsw_con.o rrsw_kg23.o rrsw_kg22.o rrsw_cld.o rrtmg_sw_setcoef.o \
	rrsw_kg26.o rrsw_kg16.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_k_g.o : $(RRTMG_SW_SRC)/rrtmg_sw_k_g.f90 rrsw_vsn.o rrsw_kg23.o \
	rrsw_kg22.o rrsw_kg17.o parkind.o rrsw_kg24.o rrsw_kg18.o rrsw_kg20.o \
	rrsw_kg29.o rrsw_kg19.o rrsw_kg27.o rrsw_kg28.o rrsw_kg21.o rrsw_kg26.o \
	rrsw_kg25.o rrsw_kg16.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_rad.o : $(RRTMG_SW_SRC)/rrtmg_sw_rad.f90 rrtmg_sw_spcvmc.o rrsw_aer.o \
	rrtmg_sw_cldprmc.o parkind.o mcica_subcol_gen_sw.o parrrsw.o rrsw_con.o \
	rrtmg_sw_setcoef.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_reftra.o : $(RRTMG_SW_SRC)/rrtmg_sw_reftra.f90 rrsw_vsn.o parkind.o \
	rrsw_tbl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_setcoef.o : $(RRTMG_SW_SRC)/rrtmg_sw_setcoef.f90 rrsw_ref.o rrsw_vsn.o \
	parrrsw.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvmc.f90 rrsw_vsn.o \
	rrtmg_sw_taumol.o rrtmg_sw_vrtqdr.o parkind.o parrrsw.o rrtmg_sw_reftra.o \
	rrsw_tbl.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvrt.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvrt.f90 rrsw_vsn.o \
	rrtmg_sw_taumol.o rrtmg_sw_vrtqdr.o parkind.o parrrsw.o rrtmg_sw_reftra.o \
	rrsw_tbl.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_taumol.o : $(RRTMG_SW_SRC)/rrtmg_sw_taumol.f90 rrsw_vsn.o rrsw_kg23.o \
	rrsw_wvn.o rrsw_kg22.o parkind.o rrsw_kg17.o rrsw_kg24.o rrsw_kg18.o parrrsw.o \
	rrsw_con.o rrsw_kg20.o rrsw_kg19.o rrsw_kg27.o rrsw_kg28.o rrsw_kg21.o \
	rrsw_kg29.o rrsw_kg26.o rrsw_kg25.o rrsw_kg16.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_vrtqdr.o : $(RRTMG_SW_SRC)/rrtmg_sw_vrtqdr.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

parrrtm.o : $(RRTMG_LW_MOD)/parrrtm.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_cld.o : $(RRTMG_LW_MOD)/rrlw_cld.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_con.o : $(RRTMG_LW_MOD)/rrlw_con.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg01.o : $(RRTMG_LW_MOD)/rrlw_kg01.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg02.o : $(RRTMG_LW_MOD)/rrlw_kg02.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg03.o : $(RRTMG_LW_MOD)/rrlw_kg03.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg04.o : $(RRTMG_LW_MOD)/rrlw_kg04.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg05.o : $(RRTMG_LW_MOD)/rrlw_kg05.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg06.o : $(RRTMG_LW_MOD)/rrlw_kg06.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg07.o : $(RRTMG_LW_MOD)/rrlw_kg07.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg08.o : $(RRTMG_LW_MOD)/rrlw_kg08.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg09.o : $(RRTMG_LW_MOD)/rrlw_kg09.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg10.o : $(RRTMG_LW_MOD)/rrlw_kg10.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg11.o : $(RRTMG_LW_MOD)/rrlw_kg11.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg12.o : $(RRTMG_LW_MOD)/rrlw_kg12.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg13.o : $(RRTMG_LW_MOD)/rrlw_kg13.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg14.o : $(RRTMG_LW_MOD)/rrlw_kg14.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg15.o : $(RRTMG_LW_MOD)/rrlw_kg15.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg16.o : $(RRTMG_LW_MOD)/rrlw_kg16.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_ncpar.o : $(RRTMG_LW_MOD)/rrlw_ncpar.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_ref.o : $(RRTMG_LW_MOD)/rrlw_ref.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_tbl.o : $(RRTMG_LW_MOD)/rrlw_tbl.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_vsn.o : $(RRTMG_LW_MOD)/rrlw_vsn.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_wvn.o : $(RRTMG_LW_MOD)/rrlw_wvn.f90 parkind.o parrrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_lw.o : $(RRTMG_LW_SRC)/mcica_subcol_gen_lw.f90 \
	mcica_random_numbers.o parkind.o rrlw_con.o rrlw_wvn.o parrrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprmc.f90 rrlw_cld.o rrlw_vsn.o \
	parkind.o rrlw_wvn.o parrrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprop.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprop.f90 rrlw_vsn.o parkind.o \
	parrrtm.o rrlw_cld.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_init.o : $(RRTMG_LW_SRC)/rrtmg_lw_init.f90 rrlw_kg15.o rrlw_kg13.o \
	rrlw_kg07.o rrlw_tbl.o rrlw_cld.o parrrtm.o rrlw_kg12.o rrlw_kg08.o rrlw_vsn.o \
	rrlw_wvn.o rrlw_kg14.o rrlw_kg05.o rrlw_kg02.o rrtmg_lw_setcoef.o rrlw_kg16.o \
	parkind.o rrlw_kg06.o rrlw_kg01.o rrlw_kg03.o rrlw_kg04.o rrlw_kg11.o \
	rrlw_con.o rrlw_kg10.o rrlw_kg09.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_k_g.o : $(RRTMG_LW_SRC)/rrtmg_lw_k_g.f90 rrlw_kg14.o rrlw_kg15.o \
	rrlw_kg04.o rrlw_kg16.o rrlw_kg07.o rrlw_kg06.o rrlw_kg08.o rrlw_kg13.o \
	rrlw_vsn.o rrlw_kg11.o parkind.o rrlw_kg12.o rrlw_kg10.o rrlw_kg09.o \
	rrlw_kg01.o rrlw_kg03.o rrlw_kg05.o rrlw_kg02.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rad.o : $(RRTMG_LW_SRC)/rrtmg_lw_rad.f90 rrtmg_lw_setcoef.o \
	rrtmg_lw_taumol.o rrtmg_lw_cldprmc.o parkind.o rrlw_con.o rrlw_wvn.o parrrtm.o \
	rrtmg_lw_rtrnmc.o mcica_subcol_gen_lw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrn.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrn.f90 rrlw_tbl.o rrlw_vsn.o \
	parkind.o rrlw_con.o rrlw_wvn.o parrrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmc.f90 rrlw_tbl.o rrlw_vsn.o \
	parkind.o rrlw_con.o rrlw_wvn.o parrrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmr.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmr.f90 rrlw_tbl.o rrlw_vsn.o \
	parkind.o rrlw_con.o rrlw_wvn.o parrrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_setcoef.o : $(RRTMG_LW_SRC)/rrtmg_lw_setcoef.f90 rrlw_ref.o rrlw_vsn.o \
	parkind.o rrlw_wvn.o parrrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_taumol.o : $(RRTMG_LW_SRC)/rrtmg_lw_taumol.f90 rrlw_kg15.o rrlw_kg13.o \
	rrlw_kg07.o parrrtm.o rrlw_kg12.o rrlw_kg08.o rrlw_vsn.o rrlw_wvn.o rrlw_kg14.o \
	rrlw_kg05.o rrlw_kg02.o rrlw_ref.o rrlw_kg16.o parkind.o rrlw_kg06.o \
	rrlw_kg01.o rrlw_kg03.o rrlw_kg04.o rrlw_kg11.o rrlw_con.o rrlw_kg10.o \
	rrlw_kg09.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_rrtm.o : $(RADIATE)/mem_rrtm.f90 rrtmg_lw_init.o parkind.o mem_grid.o \
	mem_chem1.o chem1_list.o parrrsw.o rrtmg_sw_init.o node_mod.o parrrtm.o \
	rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

isrpia.o : $(MATRIX)/isrpia.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

actv.o : $(MATRIX)/actv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

coag.o : $(MATRIX)/coag.f90 memMatrix.o setup.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

depv.o : $(MATRIX)/depv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

matrix.o : $(MATRIX)/matrix.f90 actv.o setup.o subs.o depv.o memMatrix.o npf.o \
	coag.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

npf.o : $(MATRIX)/npf.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

quad.o : $(MATRIX)/quad.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

setup.o : $(MATRIX)/setup.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

subs.o : $(MATRIX)/subs.f90 memMatrix.o setup.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

thermo_isorr.o : $(MATRIX)/thermo_isorr.f90 memMatrix.o issoropia.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

solut.o : $(MATRIX)/solut.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

issoropia.o : $(MATRIX)/issoropia.f90 isrpia.o solut.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModMatrixDriver.o : $(MATRIX)/ModMatrixDriver.F90 mem_aer1.o aer1_list.o \
	ModTurbFields.o micphys.o ModBasicFields.o mem_radiate.o mem_grid.o mem_chem1.o \
	mem_micro.o chem1_list.o node_mod.o ModParticle.o mem_leaf.o subs.o setup.o \
	isrpia.o rconstants.o memMatrix.o npf.o coag.o 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ModParticle.o : $(MATRIX)/ModParticle.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

memMatrix.o : $(MATRIX)/memMatrix.f90 aer1_list.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModMicThompsonDriver.o : $(MICRO)/ModMicThompsonDriver.f90 micphys.o \
	ModBasicFields.o io_params.o mem_micro.o mem_grid.o mem_radiate.o node_mod.o \
	mem_leaf.o module_mp_thompson.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_mp_thompson.o : $(MICRO)/module_mp_thompson.f90 module_mp_radar.o \
	node_mod.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicGfdlDriver.o : $(MICRO)/ModMicGfdlDriver.f90 mem_radiate.o \
	ModBasicFields.o mem_micro.o io_params.o mem_grid.o gfdl_cloud_microphys.o \
	node_mod.o mem_leaf.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

gfdl_cloud_microphys.o : $(MICRO)/gfdl_cloud_microphys.F90 module_mp_radar.o \
	node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_mp_radar.o : $(MICRO)/module_mp_radar.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_wind_fitch.o : $(WIND_FARM)/module_wind_fitch.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModWindFarm.o : $(WIND_FARM)/ModWindFarm.f90 ModTurbFields.o mem_tend.o \
	ModBasicFields.o io_params.o mem_grid.o node_mod.o ModNamelistFile.o \
	module_wind_fitch.o ModDateUtils.o rconstants.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

debugTools.o : $(UTILS_TOOLS)/debugTools.f90 node_mod.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOptical.o : $(RADIATE)/ModOptical.f90 aer1_list.o dump.o ModTurbFields.o \
	var_tables.o ccatt_start.o ModBasicFields.o mem_radiate.o mem_grid.o parlibf.o \
	node_mod.o ModNamelistFile.o mem_leaf.o mem_aer1.o ReadBcst.o ModControlVars.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

dam.o : $(ENERGY)/dam.f90 dump.o ModNamelistFile.o ModDateUtils.o mem_grid.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAerClim.o : $(AERCLIM)/ModAerClim.f90 dump.o ModTurbFields.o ModBasicFields.o \
	mem_grid.o parlibf.o node_mod.o ReadBcst.o ModControlVars.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitMicThompson.o : $(MICRO)/ModInitMicThompson.f90 dump.o ModBasicFields.o \
	mem_micro.o mem_grid.o parlibf.o node_mod.o generic.o ModDateUtils.o ReadBcst.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModEvaluation.o : $(EVAL)/ModEvaluation.f90 node_mod.o ModNamelistFile.o \
	mem_grid.o parlibf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modIau.o : $(MODEL)/modIau.f90 dump.o mem_tend.o mem_varinit.o mem_grid.o \
	parlibf.o node_mod.o ModNamelistFile.o ReadBcst.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modPrintInitial.o : $(INIT)/modPrintInitial.F90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

include jules_depend_model.mk

