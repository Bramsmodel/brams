actv.o : $(MATRIX)/actv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModAdapInit.o : $(INIT)/ModAdapInit.f90 mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

aer1_list.o : $(AEROSOL)/aer1_list_$(AERLEVEL).f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)
	@ln -fs aer1_list_$(AERLEVEL).o aer1_list.o

utils_f.o : $(UTILS_LIB)/utils_f.f90 dump.o ModDateUtils.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

an_header.o : $(UTILS_MODS)/an_header.f90 $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsGen.o : $(ISAN)/ModAsGen.f90 mem_grid.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsTi.o : $(ISAN)/ModAsTi.f90 isan_coms.o rconstants.o mem_grid.o \
	ModChemAObj.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsTp.o : $(ISAN)/ModAsTp.f90 rconstants.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAVarF.o : $(ISAN)/ModAVarF.f90 ModRbnd.o rconstants.o mem_grid.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

carma_fastjx.o : $(CCATT)/carma_fastjx.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ccatt_start.o : $(CCATT)/ccatt_start.f90 ModNamelistFile.o 
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

ModChemAObj.o : $(ISAN_CHEM)/ModChemAObj.f90 isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsti.o : $(ISAN_CHEM)/ModChemAsti.f90 mem_grid.o ModChemFirstRams.o \
	ModAsTi.o chem_isan_coms.o ModChemAsti2.o ModChemVInterps.o mem_aer1.o \
	ModChemAObj.o isan_coms.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsti2.o : $(ISAN_CHEM)/ModChemAsti2.f90 ModDateUtils.o rconstants.o \
	isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAstp.o : $(ISAN_CHEM)/ModChemAstp.F90 ModAsTp.o mem_chem1.o rconstants.o \
	dump.o chem1_list.o chem_isan_coms.o ModDateUtils.o mem_varinit.o isan_coms.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAvarf.o : $(ISAN_CHEM)/ModChemAvarf.f90 rconstants.o mem_grid.o ModRbnd.o \
	ModNestFeed.o chem_isan_coms.o ModAVarF.o mem_aer1.o isan_coms.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_dry_dep.o : $(MODEL_CHEM)/chem_dry_dep.f90 aer1_list.o chem1_list.o \
	ModDateUtils.o mem_aer1.o extra.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_fastjx57.o : $(CCATT)/chem_fastjx57.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_fastjx_data.o : $(CCATT)/chem_fastjx_data.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_fastjx_driv.o : $(CCATT)/chem_fastjx_driv.f90 chem_fastjx_data.o \
	rconstants.o chem_fastjx57.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFileInv.o : $(ISAN_CHEM)/ModChemFileInv.f90 dump.o ModDateUtils.o \
	mem_grid.o isan_coms.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFirstRams.o : $(ISAN_CHEM)/ModChemFirstRams.f90 grid_dims.o rconstants.o \
	mem_grid.o mem_scratch.o ModNestFillDens.o ModGetVar.o an_header.o \
	ModChemRefState.o ModRcio.o isan_coms.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_isan_coms.o : $(ISAN_CHEM)/chem_isan_coms.f90 aer1_list.o isan_coms.o \
	chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemIsanIo.o : $(ISAN_CHEM)/ModChemIsanIo.f90 isan_coms.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_plumerise_scalar.o : $(CCATT)/chem_plumerise_scalar.f90 mem_plume_chem1.o \
	node_mod.o mem_aer1.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemRefState.o : $(ISAN_CHEM)/ModChemRefState.f90 rconstants.o ccatt_start.o \
	ModNestFillDens.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_sources.o : $(CCATT)/chem_sources.f90 mem_grid.o parlibf.o \
	mem_plume_chem1.o mem_volc_chem1.o ModControlVars.o aer1_list.o \
	ModNamelistFile.o ModDateUtils.o mem_aer1.o ReadBcst.o io_params.o mem_chem1.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_dratedc.o : $(MODEL_CHEM)/chem_spack_dratedc.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_fexchem.o : $(MODEL_CHEM)/chem_spack_fexchem.f90 chem_spack_rates.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_fexloss.o : $(MODEL_CHEM)/chem_spack_fexloss.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_fexprod.o : $(MODEL_CHEM)/chem_spack_fexprod.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_jacdchemdc.o : $(MODEL_CHEM)/chem_spack_jacdchemdc.f90 \
	chem_spack_dratedc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_kinetic.o : $(MODEL_CHEM)/chem_spack_kinetic.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_lu.o : $(CCATT)/chem_spack_lu.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_qssa.o : $(CCATT)/chem_spack_qssa.f90 chem_spack_fexprod.o \
	chem_spack_kinetic.o chem_spack_rates.o chem_spack_dratedc.o \
	chem_spack_fexloss.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rates.o : $(MODEL_CHEM)/chem_spack_rates.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rodas3_dyndt.o : $(CCATT)/chem_spack_rodas3_dyndt.f90 \
	chem_spack_jacdchemdc.o mem_grid.o mem_spack.o chem_spack_kinetic.o extra.o \
	chem_spack_fexchem.o chem_spack_ros.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros.o : $(CCATT)/chem_spack_ros.f90 chem_spack_jacdchemdc.o \
	mem_spack.o chem_spack_kinetic.o chem_spack_fexchem.o chem_spack_solve_sparse.o \
	mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros_dyndt.o : $(CCATT)/chem_spack_ros_dyndt.f90 \
	chem_spack_jacdchemdc.o mem_spack.o chem_spack_kinetic.o chem_spack_fexchem.o \
	chem_spack_solve_sparse.o chem_spack_ros.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_solve_sparse.o : $(CCATT)/chem_spack_solve_sparse.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_utils.o : $(CCATT)/chem_spack_utils.f90 
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

chem_uv_att.o : $(CCATT)/chem_uv_att.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemVInterps.o : $(ISAN_CHEM)/ModChemVInterps.f90 dump.o rconstants.o \
	isan_coms.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemDryDepDriver.o : $(MODEL_CHEM)/ChemDryDepDriver.f90 mem_leaf.o grid_dims.o \
	rconstants.o mem_grid.o RadiateFields.o ModNamelistFile.o ModBasicFields.o \
	ModMicroFields.o chem_dry_dep.o mem_aer1.o ModMicControl.o CuParmFields.o \
	ModTurbFields.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemSourcesDriver.o : $(CCATT)/ChemSourcesDriver.f90 mem_leaf.o \
	mem_plume_chem1.o aer1_list.o mem_volc_chem1.o chem1_list.o chem_sources.o \
	chem_plumerise_scalar.o mem_stilt.o mem_aer1.o io_params.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

coag.o : $(MATRIX)/coag.f90 memMatrix.o setup.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModCondRead.o : $(FDDA)/ModCondRead.f90 mem_grid.o ModCondUpdate.o VarTable.o \
	ModRamsGrid.o mem_varinit.o ModDateUtils.o ModNudAnalysis.o isan_coms.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCondUpdate.o : $(FDDA)/ModCondUpdate.f90 mem_grid.o VarTable.o an_header.o \
	mem_varinit.o ModInitHis.o grid_struct.o ModRcio.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModConvComs.o : $(CUPARM)/ModConvComs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ConvPar_GF_GEOS5.o : $(CUPARM)/ConvPar_GF_GEOS5.F90 Henrys_Law_cts.o \
	MAPL_Constants.o module_gate.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuRead.o : $(CUPARM)/ModCuRead.f90 grid_dims.o mem_grid.o CuParmFields.o \
	ModNamelistFile.o ModDateUtils.o CuParmVars.o isan_coms.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

dam.o : $(ENERGY)/dam.f90 dump.o ModDateUtils.o mem_grid.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

depv.o : $(MATRIX)/depv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

digitalFilter.o : $(MODEL)/digitalFilter.f90 utilsMod.o grid_dims.o mem_grid.o \
	ModControlVars.o VarTable.o ModNamelistFile.o ModBasicFields.o ModDateUtils.o \
	node_mod.o ReadBcst.o io_params.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

domain_decomp.o : $(INIT)/domain_decomp.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

extra.o : $(MEMORY)/extra.f90 dump.o VarTable.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGetVar.o : $(UTILS_LIB)/ModGetVar.f90 dump.o an_header.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

gfdl_cloud_microphys.o : $(MICRO)/gfdl_cloud_microphys.F90 node_mod.o \
	module_mp_radar.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

grid_dims.o : $(MEMORY)/grid_dims.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

grid_struct.o : $(MEMORY)/grid_struct.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridSet.o : $(INIT)/ModGridSet.f90 grid_dims.o rconstants.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

Henrys_Law_cts.o : $(CUPARM)/Henrys_Law_cts.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

io_params.o : $(IO)/io_params.f90 grid_dims.o ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

isan_coms.o : $(ISAN_MODS)/isan_coms.f90 grid_dims.o ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

isrpia.o : $(MATRIX)/isrpia.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

issoropia.o : $(MATRIX)/issoropia.f90 isrpia.o solut.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ke_coms.o : $(TURB)/ke_coms.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Hyd.o : $(SURFACE)/ModLeaf3Hyd.f90 mem_leaf.o mem_grid.o ModLeafComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Init.o : $(SURFACE)/ModLeaf3Init.f90 mem_leaf.o teb_spm_start.o \
	grid_dims.o rconstants.o mem_grid.o ModLeaf3.o io_params.o ModLeafComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Teb.o : $(SURFACE)/ModLeaf3Teb.f90 ModGasPart.o mem_teb_vars_const.o \
	ModUrban.o mem_emiss.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeafComs.o : $(SURFACE)/ModLeafComs.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

local_proc.o : $(MODEL)/local_proc.F90 ref_sounding.o grid_dims.o rconstants.o \
	mem_grid.o dump.o mem_stilt.o node_mod.o ReadBcst.o io_params.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

machine_arq.o : $(MODEL)/machine_arq.F90 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(CMACH) $(<F:.F90=.F90)
	rm -f $(<F:.f90=.f90)

MAPL_Constants.o : $(CUPARM)/MAPL_Constants.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_random_numbers.o : $(RRTMG_SW_SRC)/mcica_random_numbers.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_lw.o : $(RRTMG_LW_SRC)/mcica_subcol_gen_lw.f90 parrrtm.o \
	rrlw_con.o mcica_random_numbers.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_sw.o : $(RRTMG_SW_SRC)/mcica_subcol_gen_sw.f90 rrsw_wvn.o \
	rrsw_con.o parrrsw.o mcica_random_numbers.o parkind.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aer1.o : $(CCATT)/mem_aer1.f90 grid_dims.o mem_grid.o dump.o aer1_list.o \
	VarTable.o ModNamelistFile.o chem1_list.o ModScalarTable.o node_mod.o \
	io_params.o mem_chem1.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

Aero2McphysFields.o : $(CCATT)/Aero2McphysFields.f90 ModNodeDimensions.o \
	VarTable.o ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aerad.o : $(RADIATE)/mem_aerad.f90 mem_grid_dim_defs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_carma.o : $(RADIATE)/mem_carma.f90 grid_dims.o mem_grid.o mem_aerad.o \
	parlibf.o ModControlVars.o VarTable.o ModRamsGrid.o ModNamelistFile.o \
	mem_globrad.o ModSoilMoisture.o ModBasicFields.o node_mod.o ModMPassFull.o \
	ReadBcst.o io_params.o ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1.o : $(CCATT)/mem_chem1.f90 grid_dims.o VarTable.o ModNamelistFile.o \
	chem1_list.o ModScalarTable.o io_params.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1aq.o : $(CCATT)/mem_chem1aq.f90 grid_dims.o VarTable.o \
	ModNamelistFile.o ModScalarTable.o chem1aq_list.o mem_chem1.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chemic.o : $(CCATT)/mem_chemic.f90 ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

CuParmVars.o : $(CUPARM)/CuParmVars.f90 grid_dims.o ModParallelEnvironment.o \
	ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

CuParmFields.o : $(CUPARM)/CuParmFields.f90 ModNodeDimensions.o VarTable.o \
	ModParallelEnvironment.o ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_emiss.o : $(TEB_SPM)/mem_emiss.f90 ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

GaspartFields.o : $(TEB_SPM)/GaspartFields.f90 ModParallelEnvironment.o \
	VarTable.o ModNamelistFile.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globaer.o : $(RADIATE)/mem_globaer.f90 mem_precision.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globrad.o : $(RADIATE)/mem_globrad.f90 mem_aerad.o mem_precision.o \
	ModNamelistFile.o parlibf.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell.o : $(CUPARM)/mem_grell.f90 VarTable.o ModNamelistFile.o \
	shcu_vars_const.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell_param2.o : $(CUPARM)/mem_grell_param2.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid.o : $(MEMORY)/mem_grid.f90 grid_dims.o ModNamelistFile.o VarTable.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid_dim_defs.o : $(MEMORY)/mem_grid_dim_defs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

JulesFields.o : $(SURFACE)/JulesFields.f90 ModParallelEnvironment.o VarTable.o \
	ModNamelistFile.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_leaf.o : $(SURFACE)/mem_leaf.f90 teb_spm_start.o grid_dims.o VarTable.o \
	ModNamelistFile.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicroFields.o : $(MICRO)/ModMicroFields.f90 ModParallelEnvironment.o \
	VarTable.o ModNamelistFile.o ModNodeDimensions.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_mksfc.o : $(MKSFC)/mem_mksfc.f90 teb_spm_start.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_nestb.o : $(NESTING)/mem_nestb.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_oda.o : $(FDDA)/mem_oda.f90 grid_dims.o ModNamelistFile.o VarTable.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_opt_scratch.o : $(TURB)/mem_opt_scratch.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_plume_chem1.o : $(CCATT)/mem_plume_chem1.f90 VarTable.o chem1_list.o \
	ModNamelistFile.o mem_chem1.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_precision.o : $(RADIATE)/mem_precision.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

RadiateFields.o : $(RADIATE)/RadiateFields.f90 ModNodeDimensions.o VarTable.o \
	ModParallelEnvironment.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_rrtm.o : $(RADIATE)/mem_rrtm.f90 rconstants.o mem_grid.o rrtmg_lw_init.o \
	parrrtm.o chem1_list.o parrrsw.o node_mod.o parkind.o rrtmg_sw_init.o \
	mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ScalarFields.o : $(MEMORY)/ScalarFields.f90 ModNodeDimensions.o VarTable.o \
	ModParallelEnvironment.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch.o : $(MEMORY)/mem_scratch.f90 ModNamelistFile.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch1_brams.o : $(MEMORY)/mem_scratch1_brams.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch1_grell.o : $(CUPARM)/mem_scratch1_grell.f90 dump.o \
	mem_grell_param2.o ccatt_start.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ShcuFields.o : $(CUPARM)/ShcuFields.f90 ModParallelEnvironment.o \
	ModControlVars.o VarTable.o ModNamelistFile.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_spack.o : $(CCATT)/mem_spack.f90 chem_spack_utils.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_stilt.o : $(STILT)/mem_stilt.f90 grid_dims.o rconstants.o VarTable.o \
	ModNamelistFile.o io_params.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tconv.o : $(CCATT)/mem_tconv.f90 aer1_list.o mem_aer1.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_teb.o : $(TEB_SPM)/mem_teb.f90 VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_teb_common.o : $(TEB_SPM)/mem_teb_common.f90 VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_teb_vars_const.o : $(TEB_SPM)/mem_teb_vars_const.f90 grid_dims.o \
	ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tend.o : $(MEMORY)/mem_tend.f90 GaspartFields.o teb_spm_start.o mem_grid.o \
	ScalarFields.o ModNamelistFile.o ModMicroFields.o ModBasicFields.o \
	ModScalarTable.o mem_emiss.o ModTurbFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_turb_scalar.o : $(TURB)/mem_turb_scalar.f90 VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tuv.o : $(TUV)/mem_tuv.f90 ModTuv2.7.o VarTable.o mem_globrad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mem_varinit.o : $(MEMORY)/mem_varinit.f90 grid_dims.o VarTable.o \
	ModNamelistFile.o chem1_list.o mem_chem1.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_volc_chem1.o : $(CCATT)/mem_volc_chem1.f90 VarTable.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

memMatrix.o : $(MATRIX)/memMatrix.f90 aer1_list.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

memSoilMoisture.o : $(SOIL_MOISTURE)/memSoilMoisture.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

meteogram.o : $(IO)/meteogram.f90 satPolyColision.o ModPostUtils.o mem_grid.o \
	VarTable.o ModNamelistFile.o node_mod.o ModMPassDtl.o meteogramType.o \
	$(POST_INCS)/post_rconstants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

meteogramType.o : $(IO)/meteogramType.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicColl.o : $(MICRO)/ModMicColl.f90 ModMicControl.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicGamma.o : $(MICRO)/ModMicGamma.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicInit.o : $(MICRO)/ModMicInit.f90 rconstants.o ModMicGamma.o mem_grid.o \
	dump.o ModMicTabs.o ModMicControl.o $(MICRO)/MicConstants.h \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicNuc.o : $(MICRO)/ModMicNuc.f90 ModMicControl.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicTabs.o : $(MICRO)/ModMicTabs.f90 ModMicControl.o ModMicGamma.o \
	$(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicVap.o : $(MICRO)/ModMicVap.f90 ModMicControl.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicControl.o : $(MICRO)/ModMicControl.f90 grid_dims.o \
	ModParallelEnvironment.o ModNamelistFile.o $(MICRO)/MicConstants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mod_GhostBlock.o : $(MODEL)/mod_GhostBlock.f90 mod_GhostBlockPartition.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mod_GhostBlockPartition.o : $(MODEL)/mod_GhostBlockPartition.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAcoust.o : $(MODEL)/ModAcoust.f90 ref_sounding.o mem_tend.o ModAcoustAdap.o \
	rconstants.o mem_grid.o mem_scratch.o ModParallelEnvironment.o ModBasicFields.o \
	node_mod.o ModGrid.o ModMicControl.o ModMessageSet.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAerClim.o : $(AERCLIM)/ModAerClim.f90 mem_grid.o parlibf.o dump.o \
	ModControlVars.o ModBasicFields.o ModSoilMoisture.o node_mod.o ReadBcst.o \
	ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBasicFields.o : $(MEMORY)/ModBasicFields.f90 ModParallelEnvironment.o \
	VarTable.o ModNamelistFile.o mem_stilt.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBramsGrid.o : $(POST_SRC)/ModBramsGrid.f90 ref_sounding.o ModPostUtils.o \
	ModParallelEnvironment.o mem_grid.o mem_aerad.o ModNamelistFile.o node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBuffering.o : $(MPI)/ModBuffering.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCarmaDriver.o : $(RADIATE)/ModCarmaDriver.f90 mem_leaf.o teb_spm_start.o \
	mem_tend.o rconstants.o RadiateFields.o rad_carma.o grid_dims.o mem_grid.o \
	mem_scratch1_grell.o ModLeaf3.o mem_teb_common.o ModNamelistFile.o \
	ModMicroFields.o mem_carma.o ModBasicFields.o ModDateUtils.o node_mod.o \
	ModMicControl.o CuParmFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsgen.o : $(ISAN_CHEM)/ModChemAsgen.F90 grid_dims.o ModMkSfcTop.o \
	node_mod.o ModChemRefState.o ModAsGen.o aer1_list.o ModChemAsti.o mem_aer1.o \
	mem_chem1.o mem_grid.o dump.o ModChemIsanIo.o ModChemAstp.o ModRamsGrid.o \
	chem_isan_coms.o ModDateUtils.o isan_coms.o ModControlVars.o chem1_list.o \
	ModChemFileInv.o ModChemAvarf.o io_params.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemConvTransp.o : $(CCATT)/ModChemConvTransp.f90 mem_scratch1_grell.o \
	mem_grid.o mem_scratch.o aer1_list.o mem_tconv.o chem1_list.o node_mod.o \
	mem_aer1.o Phys_const.o mem_grell_param2.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemistryDriver.o : $(CCATT)/ModChemistryDriver.f90 grid_dims.o \
	mem_scratch1_grell.o chem_trans_liq.o chem_fastjx_driv.o mem_chem1aq.o \
	mem_chemic.o chem_spack_utils.o node_mod.o aer1_list.o mem_globrad.o \
	mem_carma.o mem_aer1.o extra.o chem1aq_list.o chem_spack_ros.o \
	chem_trans_gasaq.o mem_chem1.o rconstants.o mem_grid.o mem_spack.o \
	chem_spack_rodas3_dyndt.o parrrtm.o ModNamelistFile.o carma_fastjx.o \
	chem_spack_qssa.o RadiateFields.o chem_uv_att.o chem_spack_ros_dyndt.o \
	mem_aerad.o mem_rrtm.o chem1_list.o ModBasicFields.o ModMicroFields.o \
	mem_stilt.o chem_spack_solve_sparse.o ModTuvDriver2.7.o mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModControlVars.o : $(INIT)/ModControlVars.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCoriolis.o : $(MODEL)/ModCoriolis.f90 ref_sounding.o mem_tend.o rconstants.o \
	mem_grid.o mem_scratch.o parlibf.o ModBasicFields.o ModBuffering.o node_mod.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuParGrell3.o : $(CUPARM)/ModCuParGrell3.F90 grid_dims.o mem_scratch1_grell.o \
	ModRstilt.o node_mod.o Phys_const.o ModRadvc.o module_cu_g3.o ModMessageSet.o \
	mem_scratch.o mem_carma.o mem_varinit.o module_cu_gf.o mem_grell.o \
	CuParmFields.o mem_chem1.o mem_leaf.o ModChemConvTransp.o rconstants.o \
	mem_grid.o ModNamelistFile.o ModMicControl.o mem_tend.o module_cu_gf_v5.1.o \
	RadiateFields.o VarTable.o ModMicroFields.o ModBasicFields.o mem_stilt.o \
	ModGrid.o ccatt_start.o io_params.o mem_grell_param2.o ConvPar_GF_GEOS5.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDateUtils.o : $(UTILS_MODS)/ModDateUtils.f90 dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffSclr.o : $(TURB)/ModDiffSclr.f90 ModTurbDiff.o mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffuse.o : $(TURB)/ModDiffuse.f90 mem_leaf.o mem_tend.o mem_grid.o \
	mem_scratch.o ke_coms.o ModTurbK.o ModDiffSclr.o ModNamelistFile.o \
	ModBasicFields.o ModTurbDiff.o ModMicroFields.o mem_opt_scratch.o \
	ModScalarTable.o node_mod.o ModTurbKE.o ModMicControl.o ModTurbFields.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDomainDecomp.o : $(MPI)/ModDomainDecomp.f90 ModParallelEnvironment.o \
	ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModEvaluation.o : $(EVAL)/ModEvaluation.f90 node_mod.o mem_grid.o \
	ModNamelistFile.o parlibf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModFieldSection.o : $(MPI)/ModFieldSection.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModFieldSectionList.o : $(MPI)/ModFieldSectionList.f90 ModParallelEnvironment.o \
	ModFieldSection.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGasPart.o : $(TEB_SPM)/ModGasPart.f90 GaspartFields.o mem_leaf.o grid_dims.o \
	mem_grid.o parlibf.o VarTable.o ModBasicFields.o an_header.o node_mod.o \
	mem_emiss.o mem_teb_vars_const.o ModRcio.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGeodat.o : $(MKSFC)/ModGeodat.f90 mem_leaf.o io_params.o teb_spm_start.o \
	mem_grid.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGrid.o : $(MPI)/ModGrid.F90 ModTurbFields.o ModMessageSet.o \
	ModParallelEnvironment.o ScalarFields.o aer1_list.o JulesFields.o \
	CuParmFields.o ModDomainDecomp.o ModGridDims.o GaspartFields.o \
	ModNamelistFile.o ShcuFields.o ModNeighbourNodes.o ModMicControl.o \
	meteogramType.o mem_tend.o RadiateFields.o VarTable.o ModControlVars.o \
	ModBasicFields.o ModMicroFields.o ModScalarTable.o Aero2McphysFields.o \
	ModNodeDimensions.o CuParmVars.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridDims.o : $(MPI)/ModGridDims.f90 ModParallelEnvironment.o \
	ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridTree.o : $(MPI)/ModGridTree.f90 ModParallelEnvironment.o ModGrid.o \
	ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modIau.o : $(MODEL)/modIau.f90 mem_tend.o mem_grid.o parlibf.o dump.o \
	ModNamelistFile.o mem_varinit.o node_mod.o ModMPassFull.o ReadBcst.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitHis.o : $(IO)/ModInitHis.f90 an_header.o ref_sounding.o mem_scratch.o \
	ModRamsReadHeader.o mem_varinit.o ModRcio.o ModLeafComs.o mem_chem1.o \
	mem_leaf.o rconstants.o mem_grid.o ModGetVar.o ModRamsGrid.o ModRinit.o \
	ModMicControl.o mem_aerad.o VarTable.o ModLeaf3.o ModBasicFields.o chem1_list.o \
	io_params.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitMicThompson.o : $(MICRO)/ModInitMicThompson.f90 mem_grid.o parlibf.o \
	dump.o ModMicroFields.o ModBasicFields.o ModDateUtils.o node_mod.o ReadBcst.o \
	generic.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLanduseInput.o : $(MKSFC)/ModLanduseInput.f90 ModLeaf3Init.o mem_mksfc.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3.o : $(SURFACE)/ModLeaf3.f90 ModLeaf3Teb.o teb_spm_start.o mem_leaf.o \
	rconstants.o mem_grid.o RadiateFields.o ModLeafComs.o mem_teb_common.o \
	ModNamelistFile.o ModBasicFields.o ModMicroFields.o node_mod.o ccatt_start.o \
	ModLeaf3Hyd.o mem_teb.o io_params.o ModMicControl.o CuParmFields.o \
	ModTurbFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3OceanOnly.o : $(SURFACE)/ModLeaf3OceanOnly.f90 mem_leaf.o \
	ModTurbFields.o rconstants.o mem_grid.o RadiateFields.o ModLeaf3.o \
	ModNamelistFile.o ModBasicFields.o node_mod.o ccatt_start.o ModCuParGrell3.o \
	io_params.o CuParmFields.o ModLeafComs.o ConvPar_GF_GEOS5.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMatrixDriver.o : $(MATRIX)/ModMatrixDriver.F90 mem_leaf.o memMatrix.o npf.o \
	rconstants.o mem_grid.o isrpia.o subs.o aer1_list.o ModParticle.o chem1_list.o \
	ModBasicFields.o ModMicroFields.o coag.o node_mod.o mem_aer1.o \
	Aero2McphysFields.o ModMicControl.o ModTurbFields.o setup.o mem_chem1.o 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ModMemAlloc.o : $(MEMORY)/ModMemAlloc.F90 grid_dims.o mem_scratch1_grell.o \
	ModTuv2.7.o mem_chem1aq.o mem_chemic.o mem_nestb.o node_mod.o ModCuParGrell3.o \
	ModTurbFields.o ModParallelEnvironment.o mem_scratch.o ScalarFields.o \
	aer1_list.o mem_volc_chem1.o mem_globrad.o mem_tuv.o mem_carma.o mem_varinit.o \
	parrrsw.o mem_teb.o mem_scratch1_brams.o shcu_vars_const.o JulesFields.o \
	mem_grell.o extra.o CuParmFields.o ModLeafComs.o mem_grid_dim_defs.o \
	chem1aq_list.o mem_leaf.o GaspartFields.o mem_chem1.o mem_globaer.o mem_grid.o \
	mem_plume_chem1.o mem_teb_common.o mem_grell_param2.o mem_aer1.o ShcuFields.o \
	carma_fastjx.o mem_emiss.o mem_teb_vars_const.o modIau.o mem_oda.o \
	teb_spm_start.o mem_tend.o RadiateFields.o mem_aerad.o VarTable.o \
	digitalFilter.o ModEvaluation.o ModMicroFields.o mem_opt_scratch.o \
	ModBasicFields.o chem1_list.o chem_sources.o ModGrid.o ccatt_start.o \
	chem_dry_dep.o mem_stilt.o ModOptical.o io_params.o Aero2McphysFields.o \
	mem_turb_scalar.o machine_arq.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMemory.o : $(UTILS_LIB)/ModMemory.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageData.o : $(MPI)/ModMessageData.f90 ModFieldSectionList.o \
	ModParallelEnvironment.o ModFieldSection.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageSet.o : $(MPI)/ModMessageSet.f90 ModParallelEnvironment.o mem_grid.o \
	parlibf.o VarTable.o ModMessageData.o ModNamelistFile.o ModFieldSection.o \
	ModNodeDimensions.o ModNeighbourNodes.o ModFieldSectionList.o ModDomainDecomp.o \
	ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicGfdlDriver.o : $(MICRO)/ModMicGfdlDriver.f90 mem_leaf.o rconstants.o \
	mem_grid.o gfdl_cloud_microphys.o ModNamelistFile.o ModMicroFields.o \
	ModBasicFields.o node_mod.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsDrive.o : $(MICRO)/ModMicrophysicsDrive.f90 ModMicNuc.o \
	grid_dims.o mem_grid.o ModMicColl.o mem_chem1aq.o mem_chemic.o ModMicInit.o \
	ModMicTabs.o ModMicroFields.o ModBasicFields.o ModMicVap.o node_mod.o \
	ModMicrophysicsMisc.o ModMicControl.o mem_chem1.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsMisc.o : $(MICRO)/ModMicrophysicsMisc.f90 rconstants.o mem_grid.o \
	ModMicroFields.o ModBasicFields.o ModMicControl.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicThompsonDriver.o : $(MICRO)/ModMicThompsonDriver.f90 mem_leaf.o \
	rconstants.o mem_grid.o module_mp_thompson.o ModNamelistFile.o ModBasicFields.o \
	ModMicroFields.o node_mod.o io_params.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcDriver.o : $(MKSFC)/ModMkSfcDriver.f90 ModMkSfcSst.o teb_spm_start.o \
	grid_dims.o mem_grid.o ModMkSfcTop.o ModControlVars.o ModSstRead.o \
	ModMkSfcFuso.o ModNdviRead.o mem_mksfc.o ModMkSfcSfc.o node_mod.o \
	ModMkSfcNdvi.o ReadBcst.o io_params.o ModLanduseInput.o ModNestGeoSst.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcFuso.o : $(MKSFC)/ModMkSfcFuso.f90 GaspartFields.o io_params.o \
	mem_grid.o ModControlVars.o mem_mksfc.o node_mod.o mem_teb.o mem_emiss.o \
	mem_teb_vars_const.o ReadBcst.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcNdvi.o : $(MKSFC)/ModMkSfcNdvi.f90 mem_leaf.o mem_grid.o ModRUser.o \
	dump.o mem_mksfc.o io_params.o ModLanduseInput.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSfc.o : $(MKSFC)/ModMkSfcSfc.f90 mem_leaf.o mem_grid.o dump.o \
	ModControlVars.o mem_mksfc.o node_mod.o ReadBcst.o io_params.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSst.o : $(MKSFC)/ModMkSfcSst.f90 mem_leaf.o grid_dims.o mem_grid.o \
	ModGeodat.o ModRUser.o mem_mksfc.o io_params.o ModNestFillDens.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcTop.o : $(MKSFC)/ModMkSfcTop.f90 mem_grid.o dump.o ModControlVars.o \
	mem_mksfc.o node_mod.o ReadBcst.o io_params.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMonotonicAdvection.o : $(MODEL)/ModMonotonicAdvection.f90 rconstants.o \
	ModParallelEnvironment.o mem_grid.o ModNamelistFile.o ModGrid.o chem_dry_dep.o \
	mem_aer1.o ccatt_start.o ModMicControl.o ModDomainDecomp.o ModMessageSet.o \
	mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNamelistFile.o : $(INIT)/ModNamelistFile.f90 grid_dims.o \
	ModParallelEnvironment.o parlibf.o dump.o modPrintInitial.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNdviRead.o : $(MKSFC)/ModNdviRead.f90 mem_leaf.o grid_dims.o mem_grid.o \
	ModControlVars.o mem_mksfc.o ModDateUtils.o node_mod.o ModMkSfcNdvi.o \
	ReadBcst.o io_params.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNeighbourNodes.o : $(MPI)/ModNeighbourNodes.f90 ModParallelEnvironment.o \
	ModDomainDecomp.o ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestGeoSst.o : $(MKSFC)/ModNestGeoSst.f90 grid_dims.o ModGeodat.o \
	ModMkSfcTop.o ModSoilMoisture.o node_mod.o memSoilMoisture.o ModLanduseInput.o \
	ModLeaf3Init.o ModTurbFields.o mem_scratch.o ModInitHis.o mem_leaf.o mem_grid.o \
	dump.o ModNestFeed.o ModRUser.o ModControlVars.o ModBasicFields.o mem_mksfc.o \
	ccatt_start.o io_params.o ModNestFillDens.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNodeDimensions.o : $(MPI)/ModNodeDimensions.f90 ModParallelEnvironment.o \
	ModDomainDecomp.o ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudAnalysis.o : $(FDDA)/ModNudAnalysis.f90 mem_tend.o mem_chem1.o mem_grid.o \
	mem_scratch.o dump.o ModEvaluation.o ModBasicFields.o chem1_list.o \
	mem_varinit.o node_mod.o ModNestFillDens.o modIau.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaNudge.o : $(FDDA)/ModOdaNudge.f90 mem_oda.o mem_tend.o mem_grid.o \
	mem_scratch.o ModOdaKrig.o ModBasicFields.o ModOdaProcObs.o node_mod.o \
	io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOneProc.o : $(MODEL)/ModOneProc.F90 grid_dims.o ModMkSfcTop.o ModOpspec.o \
	ModSoilMoisture.o ModRnode.o ModLeaf3Init.o ModChemistryDriver.o \
	mem_volc_chem1.o ModSstRead.o aer1_list.o mem_varinit.o ModCuRead.o \
	local_proc.o mem_aer1.o ModPostGridNetCDF.o ModParaInit.o ModRio.o dump.o \
	ModNestIntrp.o ModNamelistFile.o ModMPassDtl.o ModNestGeoSst.o ModRThrm.o \
	teb_spm_start.o ModBasicFields.o mem_stilt.o ModMicrophysicsMisc.o io_params.o \
	ModTuvDriver2.7.o machine_arq.o ModTuv2.7.o ModGasPart.o ModCondRead.o \
	ModMkSfcSfc.o node_mod.o ModTimestep.o domain_decomp.o ModAerClim.o mem_carma.o \
	ModInitHis.o ModDomainDecomp.o ModGridTree.o ModMkSfcDriver.o ModNudRead.o \
	mem_chem1.o mem_plume_chem1.o mem_teb_common.o ModOdaRead.o \
	mem_teb_vars_const.o mem_oda.o ModRUser.o parlibf.o digitalFilter.o \
	ModMicInit.o meteogram.o ccatt_start.o ModRanlavg.o ModPostProcess.o \
	ref_sounding.o ModUrbanCanopy.o mem_scratch.o ModMonotonicAdvection.o extra.o \
	ModMemAlloc.o ModChemAsgen.o mem_emiss.o ReadBcst.o ModRinit.o modIau.o dam.o \
	ModRecycle.o ModNdviRead.o chem_sources.o CuParmVars.o mem_grell_param2.o \
	ModWindFarm.o ModVarfFile.o mem_chem1aq.o ModMkSfcFuso.o ModVarfUpdate.o \
	tuvParameter.o ModRhhi.o ModCuParGrell3.o memSoilMoisture.o \
	ModParallelEnvironment.o ModInitMicThompson.o mem_globrad.o ModSched.o \
	mem_teb.o shcu_vars_const.o mem_leaf.o mem_grid.o ModRamsGrid.o isan_coms.o \
	ModLeaf3Teb.o ModTimeStamp.o VarTable.o ModEvaluation.o chem1_list.o \
	ModTimestepRK.o chem_dry_dep.o ModGrid.o ModCoriolis.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/tsNames.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOpspec.o : $(IO)/ModOpspec.f90 grid_dims.o mem_chem1aq.o aer1_list.o \
	mem_globrad.o mem_varinit.o mem_aer1.o shcu_vars_const.o chem1aq_list.o \
	mem_chem1.o mem_leaf.o mem_grid.o ModNamelistFile.o mem_emiss.o ModMicControl.o \
	modIau.o teb_spm_start.o chem1_list.o chem_sources.o mem_stilt.o \
	modPrintInitial.o ccatt_start.o io_params.o mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOptical.o : $(RADIATE)/ModOptical.f90 mem_leaf.o mem_grid.o parlibf.o dump.o \
	aer1_list.o ModControlVars.o VarTable.o ModRamsGrid.o ModNamelistFile.o \
	ModBasicFields.o ModSoilMoisture.o node_mod.o ccatt_start.o ModMPassFull.o \
	mem_aer1.o ModTurbFields.o ReadBcst.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOutputUtils.o : $(IO)/ModOutputUtils.f90 dump.o VarTable.o ModNamelistFile.o \
	ModBasicFields.o ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOzone.o : $(TEB_SPM)/ModOzone.f90 GaspartFields.o rconstants.o mem_grid.o \
	RadiateFields.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParaInit.o : $(MPI)/ModParaInit.f90 grid_dims.o mem_grid.o dump.o VarTable.o \
	ModScalarTable.o node_mod.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParallelEnvironment.o : $(MPI)/ModParallelEnvironment.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParticle.o : $(MATRIX)/ModParticle.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModPostGrid.o : $(POST_SRC)/ModPostGrid.F90 ModPostUtils.o \
	ModParallelEnvironment.o mem_grid.o parlibf.o VarTable.o ModOutputUtils.o \
	ModNamelistFile.o ModBasicFields.o ModPostOneFieldNetCDF.o ModBramsGrid.o \
	io_params.o ModTurbFields.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostGridNetCDF.o : $(POST_SRC)/ModPostGridNetCDF.F90 ModPostUtils.o \
	mem_grid.o dump.o ModNamelistFile.o ModDateUtils.o ModBramsGrid.o io_params.o \
	ModPostTypes.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField.o : $(POST_SRC)/ModPostOneField.f90 ModPostUtils.o \
	ModPostOneField3d.o dump.o VarTable.o ModNamelistFile.o ModBasicFields.o \
	ModPostOneField8d.o node_mod.o ModPostOneField7d.o ModBramsGrid.o \
	ModMicControl.o ModPostOneField2d.o ModPostOneFieldUtils.o ModTurbFields.o \
	ModPostTypes.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField2d.o : $(POST_SRC)/ModPostOneField2d.f90 ModPostUtils.o \
	mem_grid.o mem_aerad.o dump.o ModPostOneFieldUtils.o VarTable.o \
	ModOutputUtils.o ModNamelistFile.o ModBasicFields.o ModPostGrid.o node_mod.o \
	ModBramsGrid.o io_params.o ModMicControl.o CuParmFields.o ModTurbFields.o \
	ModPostTypes.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField3d.o : $(POST_SRC)/ModPostOneField3d.f90 ModPostUtils.o \
	mem_grid.o VarTable.o ModOutputUtils.o ModNamelistFile.o ModBasicFields.o \
	mem_varinit.o ModPostGrid.o node_mod.o ModBramsGrid.o ModMicControl.o \
	ModPostOneFieldUtils.o ModTurbFields.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField7d.o : $(POST_SRC)/ModPostOneField7d.f90 ModPostUtils.o \
	VarTable.o ModOutputUtils.o ModNamelistFile.o ModBasicFields.o ModBramsGrid.o \
	ModPostOneFieldUtils.o ModTurbFields.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField8d.o : $(POST_SRC)/ModPostOneField8d.f90 ModPostUtils.o \
	VarTable.o ModOutputUtils.o ModNamelistFile.o ModBasicFields.o ModBramsGrid.o \
	ModPostOneFieldUtils.o ModTurbFields.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneFieldNetCDF.o : $(POST_SRC)/ModPostOneFieldNetCDF.F90 \
	ModPostGridNetCDF.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneFieldUtils.o : $(POST_SRC)/ModPostOneFieldUtils.f90 ModPostGrid.o \
	ModPostTypes.o ModBramsGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostProcess.o : $(POST_SRC)/ModPostProcess.F90 ModPostOneField.o \
	ModParallelEnvironment.o ModPostTypes.o VarTable.o ModNamelistFile.o \
	ModBasicFields.o ModPostGridNetCDF.o ModPostGrid.o ModGrid.o ModBramsGrid.o \
	io_params.o ModGridTree.o ModTurbFields.o ModMessageSet.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostTypes.o : $(POST_SRC)/ModPostTypes.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostUtils.o : $(POST_SRC)/ModPostUtils.f90 dump.o ModParallelEnvironment.o \
	mem_leaf.o $(POST_INCS)/post_rconstants.h $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h $(POST_INCS)/post_rconfig.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modPrintInitial.o : $(INIT)/modPrintInitial.F90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvc.o : $(MODEL)/ModRadvc.f90 grid_dims.o ModParallelEnvironment.o \
	mem_grid.o mem_scratch.o mem_tend.o ModNamelistFile.o ModBasicFields.o \
	ModScalarTable.o chem_dry_dep.o ccatt_start.o mem_aer1.o \
	ModMonotonicAdvection.o ModRadvcAdap.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcRK.o : $(MODEL)/ModRadvcRK.f90 mem_tend.o grid_dims.o \
	ModParallelEnvironment.o mem_grid.o ModRexev.o mem_stilt.o ModGrid.o node_mod.o \
	ModMessageSet.o mem_chem1.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsReadHeader.o : $(IO)/ModRamsReadHeader.f90 an_header.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRanlavg.o : $(IO)/ModRanlavg.f90 grid_dims.o mem_grid.o VarTable.o \
	ModBasicFields.o ModMicroFields.o io_params.o ModMicControl.o ModRThrm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRbnd.o : $(BC)/ModRbnd.f90 ref_sounding.o mem_tend.o mem_grid.o mem_scratch.o \
	ModMicroFields.o ModBasicFields.o ModScalarTable.o node_mod.o ccatt_start.o \
	ModMicrophysicsMisc.o ModTurbKE.o ModMicControl.o ModTurbFields.o mem_chem1.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRcio.o : $(IO)/ModRcio.f90 mem_leaf.o grid_dims.o ModParallelEnvironment.o \
	mem_grid.o ref_sounding.o ModNamelistFile.o an_header.o mem_stilt.o io_params.o \
	ModMicControl.o ModLeafComs.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRecycle.o : $(IO)/ModRecycle.f90 mem_grid.o mem_aerad.o ModGetVar.o dump.o \
	aer1_list.o VarTable.o ModRamsReadHeader.o chem1_list.o an_header.o \
	ModDateUtils.o node_mod.o ModMPassFull.o mem_aer1.o io_params.o ReadBcst.o \
	mem_chem1.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRexev.o : $(STILT)/ModRexev.f90 mem_tend.o rconstants.o mem_grid.o \
	mem_scratch.o ModBasicFields.o mem_stilt.o ModRadvc.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRGrad.o : $(TURB)/ModRGrad.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRhhi.o : $(INIT)/ModRhhi.f90 ref_sounding.o grid_dims.o rconstants.o \
	mem_grid.o mem_scratch.o ModRamsGrid.o ModRinit.o ModBasicFields.o \
	ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRinit.o : $(INIT)/ModRinit.f90 ref_sounding.o ModTurbFields.o rconstants.o \
	mem_grid.o mem_scratch.o ModRbnd.o ModBasicFields.o ModMicroFields.o \
	mem_varinit.o node_mod.o io_params.o ModMicControl.o ModTurbKE.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRio.o : $(IO)/ModRio.f90 grid_dims.o an_header.o node_mod.o ModTurbFields.o \
	utilsMod.o ref_sounding.o ModParallelEnvironment.o ModMPassFull.o ModRcio.o \
	mem_chem1.o mem_grid.o mpi_io_engine-5d.o ModNamelistFile.o ModDateUtils.o \
	ReadBcst.o ModMicControl.o mem_aerad.o parlibf.o ModControlVars.o VarTable.o \
	ModBasicFields.o io_params.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/interface.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRrtmDriver.o : $(RADIATE)/ModRrtmDriver.f90 grid_dims.o mem_scratch1_grell.o \
	node_mod.o ref_sounding.o mem_tuv.o mem_carma.o parrrsw.o mcica_subcol_gen_lw.o \
	parkind.o CuParmFields.o ModLeafComs.o mem_chem1.o mem_leaf.o rconstants.o \
	mem_grid.o rrtmg_sw_cldprop.o parrrtm.o mcica_subcol_gen_sw.o ModNamelistFile.o \
	ModDateUtils.o rrtmg_sw_rad.o rrtmg_lw_rad.o ModMicControl.o mem_tend.o \
	teb_spm_start.o RadiateFields.o rrtmg_lw_cldprop.o mem_rrtm.o ModBasicFields.o \
	ModOptical.o ModMicroFields.o ccatt_start.o mem_grell_param2.o \
	$(UTILS_INCS)/aerosol_setup.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRShCuPar.o : $(CUPARM)/ModRShCuPar.f90 mem_tend.o shcu_vars_const.o \
	rconstants.o mem_grid.o ModMicroFields.o ModBasicFields.o ShcuFields.o \
	node_mod.o ModConvComs.o ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRstilt.o : $(STILT)/ModRstilt.f90 grid_dims.o mem_grid.o mem_scratch1_grell.o \
	mem_scratch.o ModBasicFields.o mem_stilt.o ModMonotonicAdvection.o \
	ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRThrm.o : $(MODEL)/ModRThrm.f90 rconstants.o mem_grid.o mem_scratch.o \
	ModBasicFields.o ModMicroFields.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRtimi.o : $(MODEL)/ModRtimi.f90 mem_tend.o mem_grid.o mem_scratch.o \
	ModBasicFields.o ModScalarTable.o node_mod.o shcu_vars_const.o mem_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModScalarTable.o : $(MEMORY)/ModScalarTable.f90 ModParallelEnvironment.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSched.o : $(MODEL)/ModSched.f90 ref_sounding.o ModParallelEnvironment.o \
	mem_grid.o parlibf.o dump.o ModNamelistFile.o ModBasicFields.o mem_varinit.o \
	local_proc.o node_mod.o shcu_vars_const.o io_params.o isan_coms.o ReadBcst.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSeaSalt.o : $(CCATT)/ModSeaSalt.f90 mem_leaf.o mem_grid.o aer1_list.o \
	ModAerClim.o ModBasicFields.o node_mod.o ccatt_start.o mem_aer1.o io_params.o \
	mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSstRead.o : $(MKSFC)/ModSstRead.f90 ModMkSfcSst.o mem_leaf.o grid_dims.o \
	mem_grid.o ModControlVars.o mem_mksfc.o ModDateUtils.o node_mod.o ReadBcst.o \
	io_params.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimeStamp.o : $(MODEL)/ModTimeStamp.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestep.o : $(MODEL)/ModTimestep.F90 grid_dims.o ModMicGfdlDriver.o \
	ModRstilt.o ModRtimi.o ModGasPart.o node_mod.o rad_driv.o ModMatrixDriver.o \
	ModCuParGrell3.o ModRadvc.o ModMessageSet.o ModUrbanCanopy.o ModOdaNudge.o \
	ModChemistryDriver.o mem_scratch.o ModRbnd.o ModOzone.o ModMicrophysicsDrive.o \
	mem_varinit.o mem_aer1.o shcu_vars_const.o ModMonotonicAdvection.o mem_chem1.o \
	mem_leaf.o rconstants.o mem_grid.o mem_plume_chem1.o ModSeaSalt.o ModTurbK.o \
	ChemDryDepDriver.o ModMicThompsonDriver.o mem_emiss.o ModRThrm.o \
	ModNudAnalysis.o ChemSourcesDriver.o ModCoriolis.o mem_oda.o teb_spm_start.o \
	mem_tend.o ModTimeStamp.o ModRShCuPar.o ModDiffuse.o ModLeaf3.o digitalFilter.o \
	ModRexev.o ModBasicFields.o chem_sources.o ModOptical.o mem_stilt.o ModGrid.o \
	ccatt_start.o ModMicrophysicsMisc.o ModAcoust.o sfclyr_jules.o ModWindFarm.o \
	machine_arq.o $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestepRK.o : $(MODEL)/ModTimestepRK.F90 grid_dims.o ModMicGfdlDriver.o \
	ModRstilt.o ModRtimi.o ModRadvcRK.o ModGasPart.o ModLeaf3OceanOnly.o node_mod.o \
	rad_driv.o ModMatrixDriver.o ModCuParGrell3.o ModRadvc.o ModMessageSet.o \
	utilsMod.o ModUrbanCanopy.o ModParallelEnvironment.o ModOdaNudge.o \
	ModChemistryDriver.o mem_scratch.o ModTimestep.o ModRbnd.o ModOzone.o \
	ModAerClim.o ModMicrophysicsDrive.o mem_varinit.o mem_aer1.o shcu_vars_const.o \
	ModMonotonicAdvection.o mem_chem1.o mem_leaf.o rconstants.o mem_grid.o \
	mem_plume_chem1.o ModSeaSalt.o ModTurbK.o ChemDryDepDriver.o \
	ModMicThompsonDriver.o mem_emiss.o ModRThrm.o ModNudAnalysis.o \
	ChemSourcesDriver.o ModCoriolis.o modIau.o mem_oda.o teb_spm_start.o mem_tend.o \
	ModTimeStamp.o ModRShCuPar.o ModDiffuse.o ModLeaf3.o digitalFilter.o ModRexev.o \
	chem_sources.o ModOptical.o mem_stilt.o ModGrid.o ccatt_start.o \
	ModMicrophysicsMisc.o ModAcoust.o sfclyr_jules.o ModWindFarm.o machine_arq.o \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiff.o : $(TURB)/ModTurbDiff.f90 mem_grid.o ModRGrad.o mem_scratch.o \
	mem_opt_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiffAdap.o : $(TURB)/ModTurbDiffAdap.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbFields.o : $(TURB)/ModTurbFields.f90 ModParallelEnvironment.o VarTable.o \
	ModNamelistFile.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbK.o : $(TURB)/ModTurbK.f90 grid_dims.o ModRstilt.o ModTKenn.o \
	ModTurbDiff.o node_mod.o ModTurbFields.o mem_scratch.o ModMonotonicAdvection.o \
	mem_grell.o mem_chem1.o mem_leaf.o rconstants.o mem_grid.o ModNamelistFile.o \
	ModMicControl.o mem_tend.o ModTurbDiffAdap.o ke_coms.o ModMicroFields.o \
	ModBasicFields.o ModScalarTable.o mem_stilt.o ModRGrad.o ccatt_start.o \
	ModTurbKAdap.o ModTurbKE.o mem_turb_scalar.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKAdap.o : $(TURB)/ModTurbKAdap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKE.o : $(TURB)/ModTurbKE.f90 rconstants.o mem_grid.o mem_scratch.o \
	ke_coms.o ModTurbFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuv2.7.o : $(TUV)/ModTuv2.7.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuvDriver2.7.o : $(TUV)/ModTuvDriver2.7.f90 ref_sounding.o mem_leaf.o \
	rconstants.o mem_grid.o RadiateFields.o mem_aerad.o ModTuv2.7.o mem_rrtm.o \
	mem_globrad.o ModNamelistFile.o ModBasicFields.o mem_carma.o chem1_list.o \
	mem_tuv.o node_mod.o tuvParameter.o extra.o chem_fastjx_driv.o mem_chem1.o 
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

module_mp_radar.o : $(MICRO)/module_mp_radar.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_mp_thompson.o : $(MICRO)/module_mp_thompson.f90 node_mod.o \
	module_mp_radar.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_wind_fitch.o : $(WIND_FARM)/module_wind_fitch.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModUrbanCanopy.o : $(SURFACE)/ModUrbanCanopy.f90 mem_tend.o mem_grid.o \
	ModBasicFields.o node_mod.o ModTurbFields.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfFile.o : $(FDDA)/ModVarfFile.f90 ModVarfUpdate.o node_mod.o \
	ModMessageSet.o ref_sounding.o mem_scratch.o aer1_list.o ModRamsReadHeader.o \
	mem_varinit.o mem_aer1.o ModRcio.o ModGridTree.o mem_chem1.o mem_leaf.o \
	rconstants.o mem_grid.o ModGetVar.o ModRamsGrid.o ModDateUtils.o ReadBcst.o \
	ModMicControl.o ModNudAnalysis.o isan_coms.o parlibf.o ModControlVars.o \
	ModBasicFields.o chem1_list.o ModGrid.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModWindFarm.o : $(WIND_FARM)/ModWindFarm.f90 mem_tend.o rconstants.o mem_grid.o \
	ModNamelistFile.o ModBasicFields.o ModDateUtils.o node_mod.o \
	module_wind_fitch.o io_params.o ModTurbFields.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMPassDtl.o : $(MPI)/ModMPassDtl.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMPassFull.o : $(MPI)/ModMPassFull.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mpi_io_engine-5d.o : $(IO)/mpi_io_engine-5d.f90 ModParaInit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestFeed.o : $(NESTING)/ModNestFeed.f90 mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestFillDens.o : $(NESTING)/ModNestFillDens.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestIntrp.o : $(NESTING)/ModNestIntrp.f90 ref_sounding.o grid_dims.o \
	rconstants.o mem_grid.o mem_scratch.o mem_nestb.o ModBasicFields.o ModRinit.o \
	ModNestFillDens.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

node_mod.o : $(MPI)/node_mod.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

npf.o : $(MATRIX)/npf.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModNudRead.o : $(FDDA)/ModNudRead.f90 mem_grid.o ModNudUpdate.o VarTable.o \
	ModRamsGrid.o mem_varinit.o ModDateUtils.o ModNudAnalysis.o isan_coms.o \
	mem_chem1.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudUpdate.o : $(FDDA)/ModNudUpdate.f90 mem_grid.o mem_aerad.o VarTable.o \
	chem1_list.o an_header.o mem_varinit.o ModInitHis.o grid_struct.o ModRcio.o \
	mem_chem1.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

obs_input.o : $(FDDA)/obs_input.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaKrig.o : $(FDDA)/ModOdaKrig.f90 mem_oda.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaProcObs.o : $(FDDA)/ModOdaProcObs.f90 mem_oda.o rconstants.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaRead.o : $(FDDA)/ModOdaRead.f90 mem_oda.o mem_grid.o ModDateUtils.o \
	ModOdaStaInput.o isan_coms.o ModOdaStaCount.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaCount.o : $(FDDA)/ModOdaStaCount.f90 mem_oda.o ModReadRalph.o \
	obs_input.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaInput.o : $(FDDA)/ModOdaStaInput.f90 mem_oda.o mem_grid.o \
	ModReadRalph.o ModDateUtils.o obs_input.o ModOdaStaCount.o 
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

parrrtm.o : $(RRTMG_LW_MOD)/parrrtm.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

Phys_const.o : $(CUPARM)/Phys_const.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAcoustAdap.o : $(MODEL)/ModAcoustAdap.f90 rconstants.o mem_grid.o \
	mem_scratch.o ModRbnd.o node_mod.o ModGrid.o ModMessageSet.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_carma.o : $(RADIATE)/rad_carma.F90 mem_leaf.o grid_dims.o rconstants.o \
	mem_grid.o mem_aerad.o aer1_list.o ModNamelistFile.o mem_globrad.o mem_carma.o \
	mem_tuv.o ModDateUtils.o carma_fastjx.o chem1_list.o ccatt_start.o mem_aer1.o \
	node_mod.o mem_globaer.o machine_arq.o mem_chem1.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_driv.o : $(RADIATE)/rad_driv.f90 ModCarmaDriver.o RadiateFields.o \
	ModNamelistFile.o ModMicroFields.o ModBasicFields.o ModRrtmDriver.o \
	ModMicControl.o CuParmFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcAdap.o : $(MODEL)/ModRadvcAdap.f90 ModAdapInit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsGrid.o : $(INIT)/ModRamsGrid.f90 rconstants.o mem_grid.o dump.o \
	ModGridSet.o node_mod.o ModAdapInit.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rconstants.o : $(MEMORY)/rconstants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModReadRalph.o : $(FDDA)/ModReadRalph.f90 rconstants.o obs_input.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ReadBcst.o : $(MPI)/ReadBcst.f90 utilsMod.o mem_grid.o mem_aerad.o parlibf.o \
	ModControlVars.o ModBasicFields.o an_header.o node_mod.o ModMPassFull.o \
	ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ref_sounding.o : $(MODEL)/ref_sounding.f90 grid_dims.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRnode.o : $(MODEL)/ModRnode.f90 mem_leaf.o grid_dims.o mem_grid.o parlibf.o \
	VarTable.o node_mod.o $(UTILS_INCS)/constants.h 
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

rrlw_wvn.o : $(RRTMG_LW_MOD)/rrlw_wvn.f90 parrrtm.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
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

rrtmg_lw_cldprmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprmc.f90 rrlw_vsn.o parrrtm.o \
	rrlw_cld.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprop.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprop.f90 parrrtm.o rrlw_cld.o \
	rrlw_vsn.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_init.o : $(RRTMG_LW_SRC)/rrtmg_lw_init.f90 rrlw_kg04.o rrlw_kg12.o \
	rrlw_cld.o rrlw_kg01.o rrlw_kg15.o rrlw_kg07.o rrlw_kg13.o rrlw_kg09.o \
	rrlw_kg05.o rrlw_tbl.o parkind.o rrtmg_lw_k_g.o rrlw_kg02.o rrlw_kg06.o \
	rrtmg_lw_setcoef.o rrlw_kg10.o parrrtm.o rrlw_con.o rrlw_kg11.o rrlw_kg03.o \
	rrlw_vsn.o rrlw_kg16.o rrlw_kg08.o rrlw_wvn.o rrlw_kg14.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_k_g.o : $(RRTMG_LW_SRC)/rrtmg_lw_k_g.f90 parkind.o rrlw_kg05.o \
	rrlw_kg03.o rrlw_kg11.o rrlw_vsn.o rrlw_kg04.o rrlw_kg10.o rrlw_kg12.o \
	rrlw_kg16.o rrlw_kg01.o rrlw_kg06.o rrlw_kg08.o rrlw_kg07.o rrlw_kg13.o \
	rrlw_kg09.o rrlw_kg15.o rrlw_kg02.o rrlw_kg14.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rad.o : $(RRTMG_LW_SRC)/rrtmg_lw_rad.f90 rrtmg_lw_setcoef.o parrrtm.o \
	mcica_subcol_gen_lw.o rrlw_con.o rrtmg_lw_rtrnmc.o rrtmg_lw_taumol.o \
	rrtmg_lw_cldprmc.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrn.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrn.f90 rrlw_vsn.o parrrtm.o \
	rrlw_tbl.o rrlw_con.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmc.f90 rrlw_vsn.o parrrtm.o \
	rrlw_tbl.o rrlw_con.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmr.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmr.f90 rrlw_vsn.o parrrtm.o \
	rrlw_tbl.o rrlw_con.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_setcoef.o : $(RRTMG_LW_SRC)/rrtmg_lw_setcoef.f90 rrlw_vsn.o parrrtm.o \
	rrlw_ref.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_taumol.o : $(RRTMG_LW_SRC)/rrtmg_lw_taumol.f90 rrlw_kg04.o rrlw_kg12.o \
	rrlw_kg01.o rrlw_kg15.o rrlw_kg07.o rrlw_kg13.o rrlw_kg09.o rrlw_kg05.o \
	parkind.o rrlw_kg02.o rrlw_kg06.o rrlw_kg10.o parrrtm.o rrlw_con.o rrlw_kg11.o \
	rrlw_kg03.o rrlw_vsn.o rrlw_kg16.o rrlw_kg08.o rrlw_ref.o rrlw_wvn.o \
	rrlw_kg14.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprmc.f90 rrsw_wvn.o rrsw_cld.o \
	parrrsw.o parkind.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprop.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprop.f90 rrsw_wvn.o rrsw_cld.o \
	parrrsw.o parkind.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_init.o : $(RRTMG_SW_SRC)/rrtmg_sw_init.f90 rrtmg_sw_k_g.o rrsw_wvn.o \
	rrsw_kg19.o rrsw_kg22.o rrsw_con.o rrsw_aer.o rrsw_kg24.o rrsw_kg27.o \
	rrsw_kg23.o parrrsw.o rrsw_kg28.o rrsw_kg29.o parkind.o rrsw_kg25.o rrsw_vsn.o \
	rrsw_kg17.o rrsw_cld.o rrtmg_sw_setcoef.o rrsw_kg26.o rrsw_kg20.o rrsw_kg18.o \
	rrsw_tbl.o rrsw_kg16.o rrsw_kg21.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_k_g.o : $(RRTMG_SW_SRC)/rrtmg_sw_k_g.f90 rrsw_kg19.o rrsw_kg22.o \
	rrsw_kg24.o rrsw_kg17.o rrsw_kg27.o rrsw_kg23.o rrsw_kg28.o rrsw_kg16.o \
	rrsw_kg29.o rrsw_kg25.o rrsw_kg21.o parkind.o rrsw_kg26.o rrsw_kg20.o \
	rrsw_vsn.o rrsw_kg18.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_rad.o : $(RRTMG_SW_SRC)/rrtmg_sw_rad.f90 rrsw_wvn.o rrtmg_sw_spcvmc.o \
	rrsw_con.o rrsw_aer.o mcica_subcol_gen_sw.o rrtmg_sw_cldprmc.o parrrsw.o \
	rrtmg_sw_setcoef.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_reftra.o : $(RRTMG_SW_SRC)/rrtmg_sw_reftra.f90 rrsw_vsn.o rrsw_tbl.o \
	parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_setcoef.o : $(RRTMG_SW_SRC)/rrtmg_sw_setcoef.f90 parrrsw.o rrsw_ref.o \
	rrsw_vsn.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvmc.f90 rrtmg_sw_reftra.o \
	rrsw_wvn.o rrtmg_sw_taumol.o parrrsw.o parkind.o rrtmg_sw_vrtqdr.o rrsw_vsn.o \
	rrsw_tbl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvrt.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvrt.f90 rrtmg_sw_reftra.o \
	rrsw_wvn.o rrtmg_sw_taumol.o parrrsw.o parkind.o rrtmg_sw_vrtqdr.o rrsw_vsn.o \
	rrsw_tbl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_taumol.o : $(RRTMG_SW_SRC)/rrtmg_sw_taumol.f90 rrsw_wvn.o rrsw_kg19.o \
	rrsw_kg22.o rrsw_con.o rrsw_kg24.o rrsw_kg17.o rrsw_kg27.o rrsw_kg23.o \
	parrrsw.o rrsw_kg16.o rrsw_kg28.o rrsw_kg29.o rrsw_kg25.o rrsw_kg21.o parkind.o \
	rrsw_kg26.o rrsw_kg20.o rrsw_vsn.o rrsw_kg18.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_vrtqdr.o : $(RRTMG_SW_SRC)/rrtmg_sw_vrtqdr.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRUser.o : $(SURFACE)/ModRUser.f90 mem_leaf.o rconstants.o mem_grid.o \
	node_mod.o ccatt_start.o io_params.o memSoilMoisture.o ModLeafComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

setup.o : $(MATRIX)/setup.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

sfclyr_jules.o : $(JULES_DIR)/sfclyr_jules.f90 gridbox_mean_mod.o node_mod.o \
	datetime_mod.o ancil_info.o jules_fields_mod.o model_time_mod.o ModLeaf3Init.o \
	ModTurbFields.o fluxes.o sf_diags_mod.o JulesFields.o CuParmFields.o \
	ModLeafComs.o mem_chem1.o mem_leaf.o csigma_mod.o rconstants.o mem_grid.o \
	ModNamelistFile.o jules_surface_types_mod.o ModMicControl.o gridmean_fluxes.o \
	RadiateFields.o ModBasicFields.o ModMicroFields.o chem1_list.o io_constants.o \
	mem_brams_jules.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

shcu_vars_const.o : $(CUPARM)/shcu_vars_const.f90 grid_dims.o ModNamelistFile.o \
	ModConvComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSoilMoisture.o : $(SOIL_MOISTURE)/ModSoilMoisture.F90 mem_leaf.o \
	ModTurbFields.o rconstants.o mem_grid.o mem_aerad.o parlibf.o ModControlVars.o \
	ModNamelistFile.o ModBasicFields.o node_mod.o ModMPassFull.o ReadBcst.o \
	io_params.o memSoilMoisture.o ModLeafComs.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

solut.o : $(MATRIX)/solut.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

subs.o : $(MATRIX)/subs.f90 memMatrix.o setup.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

teb_spm_start.o : $(TEB_SPM)/teb_spm_start.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTKenn.o : $(STILT)/ModTKenn.f90 rconstants.o mem_grid.o mem_scratch.o \
	turb_constants.o mem_stilt.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

turb_constants.o : $(STILT)/turb_constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

tuvParameter.o : $(TUV)/tuvParameter.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModUrban.o : $(SURFACE)/ModUrban.f90 mem_teb_vars_const.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

VarTable.o : $(MEMORY)/VarTable.f90 io_params.o ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfUpdate.o : $(FDDA)/ModVarfUpdate.f90 ref_sounding.o rconstants.o \
	mem_grid.o mem_scratch.o ModInitHis.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

include jules_depend_model.mk

