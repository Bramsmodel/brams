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

utils_f.o : $(UTILS_LIB)/utils_f.f90 ModDateUtils.o dump.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

an_header.o : $(UTILS_MODS)/an_header.f90 $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsGen.o : $(ISAN)/ModAsGen.f90 mem_grid.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsTi.o : $(ISAN)/ModAsTi.f90 mem_grid.o ModChemAObj.o rconstants.o \
	isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsTp.o : $(ISAN)/ModAsTp.f90 rconstants.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAVarF.o : $(ISAN)/ModAVarF.f90 ModRbnd.o mem_grid.o rconstants.o isan_coms.o 
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

ModChemAsti.o : $(ISAN_CHEM)/ModChemAsti.f90 mem_grid.o chem_isan_coms.o \
	isan_coms.o ModChemAsti2.o mem_chem1.o ModChemVInterps.o ModChemAObj.o \
	mem_aer1.o ModChemFirstRams.o ModAsTi.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsti2.o : $(ISAN_CHEM)/ModChemAsti2.f90 ModDateUtils.o rconstants.o \
	isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAstp.o : $(ISAN_CHEM)/ModChemAstp.F90 ModAsTp.o chem_isan_coms.o \
	isan_coms.o mem_varinit.o mem_chem1.o chem1_list.o rconstants.o ModDateUtils.o \
	dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAvarf.o : $(ISAN_CHEM)/ModChemAvarf.f90 mem_grid.o ModRbnd.o \
	chem_isan_coms.o isan_coms.o mem_chem1.o ModNestFeed.o rconstants.o ModAVarF.o \
	mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_dry_dep.o : $(MODEL_CHEM)/chem_dry_dep.f90 extra.o mem_chem1.o mem_aer1.o \
	chem1_list.o ModDateUtils.o aer1_list.o 
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

chem_fastjx_driv.o : $(CCATT)/chem_fastjx_driv.f90 rconstants.o \
	chem_fastjx_data.o chem1_list.o chem_fastjx57.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFileInv.o : $(ISAN_CHEM)/ModChemFileInv.f90 mem_grid.o ModDateUtils.o \
	dump.o isan_coms.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFirstRams.o : $(ISAN_CHEM)/ModChemFirstRams.f90 ModChemRefState.o \
	mem_grid.o isan_coms.o grid_dims.o ModNestFillDens.o rconstants.o ModGetVar.o \
	an_header.o ModRcio.o mem_scratch.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_isan_coms.o : $(ISAN_CHEM)/chem_isan_coms.f90 chem1_list.o aer1_list.o \
	isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemIsanIo.o : $(ISAN_CHEM)/ModChemIsanIo.f90 isan_coms.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_plumerise_scalar.o : $(CCATT)/chem_plumerise_scalar.f90 mem_aer1.o \
	node_mod.o mem_plume_chem1.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemRefState.o : $(ISAN_CHEM)/ModChemRefState.f90 ModNestFillDens.o \
	rconstants.o ccatt_start.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_sources.o : $(CCATT)/chem_sources.f90 mem_grid.o ModNamelistFile.o \
	io_params.o mem_volc_chem1.o ModControlVars.o mem_chem1.o mem_plume_chem1.o \
	mem_aer1.o ReadBcst.o ModDateUtils.o parlibf.o aer1_list.o \
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

chem_spack_qssa.o : $(CCATT)/chem_spack_qssa.f90 chem_spack_dratedc.o \
	chem_spack_kinetic.o mem_chem1.o chem_spack_fexprod.o chem_spack_fexloss.o \
	chem_spack_rates.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rates.o : $(MODEL_CHEM)/chem_spack_rates.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rodas3_dyndt.o : $(CCATT)/chem_spack_rodas3_dyndt.f90 mem_grid.o \
	extra.o chem_spack_ros.o chem_spack_kinetic.o chem_spack_fexchem.o \
	chem_spack_jacdchemdc.o mem_spack.o mem_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros.o : $(CCATT)/chem_spack_ros.f90 chem_spack_kinetic.o \
	chem_spack_fexchem.o chem_spack_jacdchemdc.o mem_spack.o mem_chem1.o \
	chem_spack_solve_sparse.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros_dyndt.o : $(CCATT)/chem_spack_ros_dyndt.f90 chem_spack_ros.o \
	chem_spack_kinetic.o chem_spack_fexchem.o chem_spack_jacdchemdc.o mem_spack.o \
	mem_chem1.o chem_spack_solve_sparse.o 
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

ModChemVInterps.o : $(ISAN_CHEM)/ModChemVInterps.f90 rconstants.o dump.o \
	isan_coms.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemDryDepDriver.o : $(MODEL_CHEM)/ChemDryDepDriver.f90 mem_grid.o mem_leaf.o \
	RadiateFields.o ModNamelistFile.o ModTurbFields.o CuParmFields.o mem_chem1.o \
	grid_dims.o rconstants.o ModMicroFields.o mem_aer1.o ModMicControl.o \
	chem_dry_dep.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemSourcesDriver.o : $(CCATT)/ChemSourcesDriver.f90 mem_leaf.o mem_stilt.o \
	io_params.o chem_plumerise_scalar.o mem_volc_chem1.o mem_chem1.o \
	mem_plume_chem1.o chem1_list.o mem_aer1.o aer1_list.o chem_sources.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

coag.o : $(MATRIX)/coag.f90 setup.o memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModCondRead.o : $(FDDA)/ModCondRead.f90 mem_grid.o isan_coms.o VarTable.o \
	mem_varinit.o ModCondUpdate.o ModNudAnalysis.o ModDateUtils.o ModRamsGrid.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCondUpdate.o : $(FDDA)/ModCondUpdate.f90 mem_grid.o VarTable.o mem_varinit.o \
	grid_struct.o ModInitHis.o an_header.o ModRcio.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModConvComs.o : $(CUPARM)/ModConvComs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ConvPar_GF_GEOS5.o : $(CUPARM)/ConvPar_GF_GEOS5.F90 module_gate.o \
	Henrys_Law_cts.o MAPL_Constants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuRead.o : $(CUPARM)/ModCuRead.f90 mem_grid.o ModNamelistFile.o isan_coms.o \
	CuParmFields.o CuParmVars.o grid_dims.o ModDateUtils.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

dam.o : $(ENERGY)/dam.f90 mem_grid.o ModDateUtils.o ModNamelistFile.o dump.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

depv.o : $(MATRIX)/depv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

digitalFilter.o : $(MODEL)/digitalFilter.f90 mem_grid.o ModNamelistFile.o \
	io_params.o VarTable.o ModControlVars.o node_mod.o grid_dims.o ReadBcst.o \
	ModDateUtils.o utilsMod.o ModBasicFields.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

domain_decomp.o : $(INIT)/domain_decomp.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

extra.o : $(MEMORY)/extra.f90 VarTable.o ModNamelistFile.o dump.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGetVar.o : $(UTILS_LIB)/ModGetVar.f90 dump.o an_header.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
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

ModGridSet.o : $(INIT)/ModGridSet.f90 mem_grid.o grid_dims.o rconstants.o 
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

ModLeaf3Hyd.o : $(SURFACE)/ModLeaf3Hyd.f90 mem_grid.o mem_leaf.o ModLeafComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Init.o : $(SURFACE)/ModLeaf3Init.f90 mem_grid.o mem_leaf.o io_params.o \
	teb_spm_start.o grid_dims.o rconstants.o ModLeaf3.o ModLeafComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Teb.o : $(SURFACE)/ModLeaf3Teb.f90 mem_teb_vars_const.o mem_emiss.o \
	ModUrban.o ModGasPart.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeafComs.o : $(SURFACE)/ModLeafComs.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

local_proc.o : $(MODEL)/local_proc.F90 mem_grid.o mem_stilt.o io_params.o \
	ref_sounding.o grid_dims.o rconstants.o node_mod.o ReadBcst.o dump.o \
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

mcica_subcol_gen_lw.o : $(RRTMG_LW_SRC)/mcica_subcol_gen_lw.f90 parkind.o \
	rrlw_con.o parrrtm.o rrlw_wvn.o mcica_random_numbers.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_sw.o : $(RRTMG_SW_SRC)/mcica_subcol_gen_sw.f90 parrrsw.o \
	parkind.o rrsw_con.o rrsw_wvn.o mcica_random_numbers.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aer1.o : $(CCATT)/mem_aer1.f90 mem_grid.o ModNamelistFile.o io_params.o \
	ModScalarTable.o VarTable.o mem_chem1.o grid_dims.o chem1_list.o node_mod.o \
	dump.o aer1_list.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

Aero2McphysFields.o : $(CCATT)/Aero2McphysFields.f90 ModParallelEnvironment.o \
	VarTable.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aerad.o : $(RADIATE)/mem_aerad.f90 mem_grid_dim_defs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_carma.o : $(RADIATE)/mem_carma.f90 mem_grid.o ModNamelistFile.o io_params.o \
	mem_aerad.o mem_globrad.o ModMPassFull.o VarTable.o ModControlVars.o \
	ModTurbFields.o grid_dims.o node_mod.o ModSoilMoisture.o ReadBcst.o parlibf.o \
	ModRamsGrid.o ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1.o : $(CCATT)/mem_chem1.f90 ModNamelistFile.o io_params.o \
	ModScalarTable.o VarTable.o grid_dims.o chem1_list.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1aq.o : $(CCATT)/mem_chem1aq.f90 ModNamelistFile.o ModScalarTable.o \
	VarTable.o mem_chem1.o grid_dims.o chem1aq_list.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chemic.o : $(CCATT)/mem_chemic.f90 ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

CuParmVars.o : $(CUPARM)/CuParmVars.f90 ModParallelEnvironment.o grid_dims.o \
	ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

CuParmFields.o : $(CUPARM)/CuParmFields.f90 ModParallelEnvironment.o \
	ModNamelistFile.o VarTable.o ModNodeDimensions.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_emiss.o : $(TEB_SPM)/mem_emiss.f90 ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

GaspartFields.o : $(TEB_SPM)/GaspartFields.f90 ModNamelistFile.o VarTable.o \
	ModParallelEnvironment.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globaer.o : $(RADIATE)/mem_globaer.f90 mem_precision.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globrad.o : $(RADIATE)/mem_globrad.f90 mem_precision.o parlibf.o \
	ModNamelistFile.o mem_aerad.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell.o : $(CUPARM)/mem_grell.f90 shcu_vars_const.o ModNamelistFile.o \
	VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell_param2.o : $(CUPARM)/mem_grell_param2.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid.o : $(MEMORY)/mem_grid.f90 VarTable.o grid_dims.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid_dim_defs.o : $(MEMORY)/mem_grid_dim_defs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

JulesFields.o : $(SURFACE)/JulesFields.f90 ModNamelistFile.o \
	ModParallelEnvironment.o VarTable.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_leaf.o : $(SURFACE)/mem_leaf.f90 ModNamelistFile.o io_params.o VarTable.o \
	teb_spm_start.o grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicroFields.o : $(MICRO)/ModMicroFields.f90 ModNamelistFile.o \
	ModParallelEnvironment.o VarTable.o ModMicControl.o ModNodeDimensions.o 
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

mem_oda.o : $(FDDA)/mem_oda.f90 VarTable.o grid_dims.o ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_opt_scratch.o : $(TURB)/mem_opt_scratch.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_plume_chem1.o : $(CCATT)/mem_plume_chem1.f90 VarTable.o ModNamelistFile.o \
	chem1_list.o mem_chem1.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_precision.o : $(RADIATE)/mem_precision.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

RadiateFields.o : $(RADIATE)/RadiateFields.f90 ModParallelEnvironment.o \
	ModNamelistFile.o VarTable.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_rrtm.o : $(RADIATE)/mem_rrtm.f90 mem_grid.o parrrsw.o parkind.o parrrtm.o \
	mem_chem1.o rrtmg_lw_init.o node_mod.o rconstants.o chem1_list.o \
	rrtmg_sw_init.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ScalarFields.o : $(MEMORY)/ScalarFields.f90 ModParallelEnvironment.o \
	ModNamelistFile.o VarTable.o ModNodeDimensions.o 
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

ShcuFields.o : $(CUPARM)/ShcuFields.f90 ModNamelistFile.o \
	ModParallelEnvironment.o ModControlVars.o VarTable.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_spack.o : $(CCATT)/mem_spack.f90 chem_spack_utils.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_stilt.o : $(STILT)/mem_stilt.f90 ModNamelistFile.o io_params.o VarTable.o \
	grid_dims.o rconstants.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tconv.o : $(CCATT)/mem_tconv.f90 mem_aer1.o chem1_list.o aer1_list.o 
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

mem_tend.o : $(MEMORY)/mem_tend.f90 mem_grid.o ModTurbFields.o ScalarFields.o \
	ModNamelistFile.o ModScalarTable.o teb_spm_start.o mem_emiss.o GaspartFields.o \
	ModMicroFields.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_turb_scalar.o : $(TURB)/mem_turb_scalar.f90 VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tuv.o : $(TUV)/mem_tuv.f90 VarTable.o ModTuv2.7.o mem_globrad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mem_varinit.o : $(MEMORY)/mem_varinit.f90 ModNamelistFile.o VarTable.o \
	mem_chem1.o grid_dims.o chem1_list.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_volc_chem1.o : $(CCATT)/mem_volc_chem1.f90 VarTable.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

memMatrix.o : $(MATRIX)/memMatrix.f90 ModNamelistFile.o aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

memSoilMoisture.o : $(SOIL_MOISTURE)/memSoilMoisture.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

meteogram.o : $(IO)/meteogram.f90 mem_grid.o ModNamelistFile.o VarTable.o \
	ModMPassDtl.o meteogramType.o satPolyColision.o ModPostUtils.o node_mod.o \
	$(UTILS_INCS)/files.h $(POST_INCS)/post_rconstants.h 
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

ModMicInit.o : $(MICRO)/ModMicInit.f90 mem_grid.o ModMicGamma.o rconstants.o \
	ModMicTabs.o ModMicControl.o dump.o $(MICRO)/MicConstants.h \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
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

ModMicVap.o : $(MICRO)/ModMicVap.f90 rconstants.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicControl.o : $(MICRO)/ModMicControl.f90 ModParallelEnvironment.o \
	grid_dims.o ModNamelistFile.o $(MICRO)/MicConstants.h $(UTILS_INCS)/files.h 
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

ModAcoust.o : $(MODEL)/ModAcoust.f90 mem_grid.o mem_scratch.o ModBasicFields.o \
	ref_sounding.o mem_tend.o ModParallelEnvironment.o ModGrid.o node_mod.o \
	rconstants.o ModMessageSet.o ModMicControl.o ModAcoustAdap.o \
	$(UTILS_INCS)/tsNames.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAerClim.o : $(AERCLIM)/ModAerClim.f90 mem_grid.o ModTurbFields.o \
	ModControlVars.o node_mod.o ModSoilMoisture.o ReadBcst.o parlibf.o dump.o \
	ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBasicFields.o : $(MEMORY)/ModBasicFields.f90 mem_stilt.o ModNamelistFile.o \
	ModParallelEnvironment.o VarTable.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBramsGrid.o : $(POST_SRC)/ModBramsGrid.f90 mem_grid.o ModNamelistFile.o \
	mem_aerad.o ref_sounding.o ModParallelEnvironment.o ModPostUtils.o node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBuffering.o : $(MPI)/ModBuffering.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCarmaDriver.o : $(RADIATE)/ModCarmaDriver.f90 mem_leaf.o mem_grid.o \
	ModDateUtils.o ModNamelistFile.o RadiateFields.o rad_carma.o teb_spm_start.o \
	mem_carma.o mem_tend.o mem_scratch1_grell.o node_mod.o mem_teb_common.o \
	rconstants.o ModLeaf3.o grid_dims.o ModMicroFields.o CuParmFields.o \
	ModMicControl.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsgen.o : $(ISAN_CHEM)/ModChemAsgen.F90 ModChemFileInv.o chem_isan_coms.o \
	io_params.o ModChemAvarf.o node_mod.o isan_coms.o ModControlVars.o grid_dims.o \
	ModChemIsanIo.o dump.o ModChemAstp.o mem_grid.o ModAsGen.o mem_chem1.o \
	ModRamsGrid.o ModChemRefState.o ModDateUtils.o ModMkSfcTop.o chem1_list.o \
	mem_aer1.o ModChemAsti.o aer1_list.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemConvTransp.o : $(CCATT)/ModChemConvTransp.f90 mem_grid.o mem_tconv.o \
	mem_scratch1_grell.o Phys_const.o mem_chem1.o node_mod.o chem1_list.o \
	mem_grell_param2.o mem_aer1.o aer1_list.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemistryDriver.o : $(CCATT)/ModChemistryDriver.f90 chem_spack_rodas3_dyndt.o \
	RadiateFields.o chem_trans_liq.o mem_carma.o carma_fastjx.o chem_spack_utils.o \
	chem_trans_gasaq.o node_mod.o rconstants.o mem_aer1.o mem_rrtm.o \
	ModNamelistFile.o mem_aerad.o chem_spack_ros.o mem_chemic.o chem_fastjx_driv.o \
	grid_dims.o chem_spack_qssa.o ModBasicFields.o mem_grid.o mem_stilt.o \
	mem_globrad.o extra.o mem_scratch1_grell.o chem_uv_att.o mem_chem1.o \
	chem_spack_solve_sparse.o chem_spack_ros_dyndt.o ModTuvDriver2.7.o mem_spack.o \
	parrrtm.o chem1_list.o ModMicroFields.o mem_chem1aq.o aer1_list.o \
	chem1aq_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModControlVars.o : $(INIT)/ModControlVars.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCoriolis.o : $(MODEL)/ModCoriolis.f90 mem_grid.o mem_scratch.o ref_sounding.o \
	mem_tend.o node_mod.o rconstants.o parlibf.o ModBuffering.o ModBasicFields.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuParGrell3.o : $(CUPARM)/ModCuParGrell3.F90 mem_leaf.o RadiateFields.o \
	io_params.o module_cu_g3.o mem_varinit.o mem_carma.o ModGrid.o node_mod.o \
	rconstants.o mem_grell_param2.o ModMicControl.o ModNamelistFile.o VarTable.o \
	ModRadvc.o grid_dims.o ModRstilt.o ModChemConvTransp.o CuParmFields.o \
	ModBasicFields.o mem_stilt.o mem_grid.o mem_tend.o mem_grell.o \
	mem_scratch1_grell.o ccatt_start.o mem_chem1.o module_cu_gf_v5.1.o \
	ConvPar_GF_GEOS5.o ModMessageSet.o Phys_const.o ModMicroFields.o mem_scratch.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDateUtils.o : $(UTILS_MODS)/ModDateUtils.f90 dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffSclr.o : $(TURB)/ModDiffSclr.f90 mem_grid.o ModTurbDiff.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffuse.o : $(TURB)/ModDiffuse.f90 mem_grid.o mem_leaf.o ModTurbFields.o \
	ModNamelistFile.o ModBasicFields.o ModScalarTable.o mem_tend.o ModDiffSclr.o \
	mem_opt_scratch.o ke_coms.o ModTurbKE.o node_mod.o ModTurbDiff.o \
	ModMicroFields.o ModMicControl.o ModTurbK.o mem_scratch.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDomainDecomp.o : $(MPI)/ModDomainDecomp.f90 ModParallelEnvironment.o \
	ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModEvaluation.o : $(EVAL)/ModEvaluation.f90 mem_grid.o node_mod.o parlibf.o \
	ModNamelistFile.o 
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

ModGasPart.o : $(TEB_SPM)/ModGasPart.f90 mem_grid.o mem_leaf.o \
	mem_teb_vars_const.o VarTable.o mem_emiss.o grid_dims.o node_mod.o \
	GaspartFields.o parlibf.o an_header.o ModRcio.o ModBasicFields.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGeodat.o : $(MKSFC)/ModGeodat.f90 mem_grid.o teb_spm_start.o io_params.o \
	mem_leaf.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGrid.o : $(MPI)/ModGrid.F90 RadiateFields.o ModTurbFields.o JulesFields.o \
	meteogramType.o ModNeighbourNodes.o ModMicControl.o ModNodeDimensions.o \
	ModGridDims.o ModNamelistFile.o ModParallelEnvironment.o VarTable.o \
	ModControlVars.o Aero2McphysFields.o ShcuFields.o ModDomainDecomp.o \
	CuParmFields.o ModBasicFields.o ModScalarTable.o mem_tend.o CuParmVars.o \
	ModMessageSet.o ScalarFields.o GaspartFields.o ModMicroFields.o aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridDims.o : $(MPI)/ModGridDims.f90 ModParallelEnvironment.o \
	ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridTree.o : $(MPI)/ModGridTree.f90 ModParallelEnvironment.o \
	ModNamelistFile.o ModGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modIau.o : $(MODEL)/modIau.f90 mem_grid.o ModNamelistFile.o ModMPassFull.o \
	mem_tend.o mem_varinit.o node_mod.o ReadBcst.o parlibf.o dump.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitHis.o : $(IO)/ModInitHis.f90 mem_leaf.o io_params.o mem_varinit.o \
	ModRamsReadHeader.o rconstants.o ModMicControl.o mem_aerad.o VarTable.o \
	ModRinit.o ModGetVar.o ModBasicFields.o mem_grid.o ref_sounding.o mem_chem1.o \
	ModLeaf3.o ModRamsGrid.o ModLeafComs.o chem1_list.o an_header.o ModRcio.o \
	mem_scratch.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitMicThompson.o : $(MICRO)/ModInitMicThompson.f90 mem_grid.o node_mod.o \
	ReadBcst.o generic.o ModMicroFields.o ModDateUtils.o parlibf.o dump.o \
	ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLanduseInput.o : $(MKSFC)/ModLanduseInput.f90 mem_mksfc.o ModLeaf3Init.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3.o : $(SURFACE)/ModLeaf3.f90 mem_grid.o mem_leaf.o mem_teb.o \
	ModNamelistFile.o io_params.o ModTurbFields.o RadiateFields.o teb_spm_start.o \
	ModLeaf3Hyd.o ccatt_start.o ModLeafComs.o node_mod.o rconstants.o \
	mem_teb_common.o ModLeaf3Teb.o ModMicroFields.o CuParmFields.o ModMicControl.o \
	ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3OceanOnly.o : $(SURFACE)/ModLeaf3OceanOnly.f90 mem_leaf.o mem_grid.o \
	RadiateFields.o ModNamelistFile.o ModTurbFields.o io_params.o ModBasicFields.o \
	ccatt_start.o ModCuParGrell3.o node_mod.o rconstants.o ConvPar_GF_GEOS5.o \
	ModLeaf3.o CuParmFields.o ModLeafComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMatrixDriver.o : $(MATRIX)/ModMatrixDriver.F90 mem_grid.o mem_leaf.o \
	ModTurbFields.o memMatrix.o ModParticle.o setup.o Aero2McphysFields.o \
	mem_chem1.o coag.o node_mod.o chem1_list.o rconstants.o isrpia.o \
	ModMicroFields.o mem_aer1.o ModMicControl.o subs.o npf.o aer1_list.o \
	ModBasicFields.o 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ModMemAlloc.o : $(MEMORY)/ModMemAlloc.F90 mem_leaf.o ModTurbFields.o \
	RadiateFields.o io_params.o mem_tuv.o JulesFields.o mem_teb_vars_const.o \
	mem_varinit.o mem_carma.o mem_scratch.o carma_fastjx.o ModGrid.o \
	ModCuParGrell3.o node_mod.o mem_grell_param2.o mem_aerad.o digitalFilter.o \
	VarTable.o ModParallelEnvironment.o mem_chemic.o mem_oda.o Aero2McphysFields.o \
	mem_plume_chem1.o grid_dims.o mem_teb_common.o ShcuFields.o mem_globaer.o \
	machine_arq.o modIau.o mem_chem1aq.o CuParmFields.o ModBasicFields.o mem_grid.o \
	chem_sources.o mem_stilt.o mem_volc_chem1.o mem_globrad.o extra.o mem_tend.o \
	mem_grell.o mem_opt_scratch.o mem_scratch1_grell.o mem_emiss.o ccatt_start.o \
	mem_chem1.o ModOptical.o mem_scratch1_brams.o chem_dry_dep.o shcu_vars_const.o \
	mem_teb.o chem1aq_list.o ScalarFields.o parrrsw.o teb_spm_start.o \
	mem_grid_dim_defs.o ModTuv2.7.o ModEvaluation.o ModLeafComs.o chem1_list.o \
	mem_nestb.o GaspartFields.o ModMicroFields.o mem_aer1.o aer1_list.o \
	mem_turb_scalar.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMemory.o : $(UTILS_LIB)/ModMemory.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageData.o : $(MPI)/ModMessageData.f90 ModParallelEnvironment.o \
	ModFieldSectionList.o ModFieldSection.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageSet.o : $(MPI)/ModMessageSet.f90 mem_grid.o ModGridDims.o \
	ModNamelistFile.o ModParallelEnvironment.o ModFieldSectionList.o VarTable.o \
	ModMessageData.o ModNeighbourNodes.o ModDomainDecomp.o parlibf.o \
	ModFieldSection.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicGfdlDriver.o : $(MICRO)/ModMicGfdlDriver.f90 mem_grid.o mem_leaf.o \
	ModNamelistFile.o gfdl_cloud_microphys.o io_params.o node_mod.o rconstants.o \
	ModMicroFields.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsDrive.o : $(MICRO)/ModMicrophysicsDrive.f90 mem_grid.o \
	mem_chemic.o ModMicrophysicsMisc.o ModMicNuc.o mem_chem1.o grid_dims.o \
	ModMicTabs.o ModMicColl.o ModMicVap.o node_mod.o ModMicroFields.o mem_chem1aq.o \
	ModMicControl.o ModMicInit.o ModBasicFields.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsMisc.o : $(MICRO)/ModMicrophysicsMisc.f90 mem_grid.o rconstants.o \
	ModMicroFields.o ModMicControl.o ModBasicFields.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicThompsonDriver.o : $(MICRO)/ModMicThompsonDriver.f90 mem_grid.o mem_leaf.o \
	ModNamelistFile.o io_params.o module_mp_thompson.o node_mod.o rconstants.o \
	ModMicroFields.o ModMicControl.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcDriver.o : $(MKSFC)/ModMkSfcDriver.f90 mem_grid.o io_params.o \
	ModMkSfcSst.o ModNestGeoSst.o ModControlVars.o ModLanduseInput.o \
	teb_spm_start.o ModNdviRead.o ModMkSfcFuso.o grid_dims.o ModMkSfcTop.o \
	node_mod.o ModMkSfcSfc.o ReadBcst.o ModMkSfcNdvi.o mem_mksfc.o ModSstRead.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcFuso.o : $(MKSFC)/ModMkSfcFuso.f90 mem_grid.o mem_teb.o io_params.o \
	mem_teb_vars_const.o ModControlVars.o mem_emiss.o node_mod.o GaspartFields.o \
	ReadBcst.o mem_mksfc.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcNdvi.o : $(MKSFC)/ModMkSfcNdvi.f90 mem_grid.o mem_leaf.o io_params.o \
	ModRUser.o ModLanduseInput.o dump.o mem_mksfc.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSfc.o : $(MKSFC)/ModMkSfcSfc.f90 mem_grid.o mem_leaf.o io_params.o \
	ModControlVars.o node_mod.o ReadBcst.o dump.o mem_mksfc.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSst.o : $(MKSFC)/ModMkSfcSst.f90 mem_leaf.o ModGeodat.o mem_grid.o \
	io_params.o ModRUser.o grid_dims.o ModNestFillDens.o mem_mksfc.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcTop.o : $(MKSFC)/ModMkSfcTop.f90 mem_grid.o io_params.o ModControlVars.o \
	node_mod.o ReadBcst.o dump.o mem_mksfc.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMonotonicAdvection.o : $(MODEL)/ModMonotonicAdvection.f90 mem_grid.o \
	ModNamelistFile.o ModParallelEnvironment.o ccatt_start.o ModGrid.o mem_chem1.o \
	rconstants.o ModDomainDecomp.o mem_aer1.o ModMicControl.o chem_dry_dep.o \
	ModMessageSet.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNamelistFile.o : $(INIT)/ModNamelistFile.f90 ModParallelEnvironment.o \
	modPrintInitial.o grid_dims.o parlibf.o dump.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNdviRead.o : $(MKSFC)/ModNdviRead.f90 mem_grid.o mem_leaf.o io_params.o \
	ModControlVars.o node_mod.o grid_dims.o ReadBcst.o ModDateUtils.o \
	ModMkSfcNdvi.o mem_mksfc.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNeighbourNodes.o : $(MPI)/ModNeighbourNodes.f90 ModParallelEnvironment.o \
	ModDomainDecomp.o ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestGeoSst.o : $(MKSFC)/ModNestGeoSst.f90 mem_leaf.o memSoilMoisture.o \
	ModTurbFields.o io_params.o mem_scratch.o node_mod.o ModRUser.o \
	ModControlVars.o grid_dims.o ModNestFillDens.o ModSoilMoisture.o dump.o \
	ModBasicFields.o mem_grid.o ModGeodat.o ccatt_start.o ModNestFeed.o mem_mksfc.o \
	ModLanduseInput.o ModMkSfcTop.o ModInitHis.o ModLeaf3Init.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNodeDimensions.o : $(MPI)/ModNodeDimensions.f90 ModParallelEnvironment.o \
	ModDomainDecomp.o ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudAnalysis.o : $(FDDA)/ModNudAnalysis.f90 mem_grid.o ModBasicFields.o \
	mem_tend.o mem_varinit.o ModEvaluation.o mem_chem1.o node_mod.o \
	ModNestFillDens.o chem1_list.o modIau.o dump.o mem_scratch.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaNudge.o : $(FDDA)/ModOdaNudge.f90 mem_grid.o ModOdaKrig.o mem_scratch.o \
	io_params.o mem_tend.o mem_oda.o node_mod.o ModOdaProcObs.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOneProc.o : $(MODEL)/ModOneProc.F90 ModMkSfcDriver.o ModUrbanCanopy.o \
	node_mod.o ModOpspec.o ModTimestep.o digitalFilter.o ModRinit.o mem_oda.o \
	ModVarfUpdate.o mem_globrad.o mem_volc_chem1.o ref_sounding.o mem_chem1.o \
	ModRhhi.o parlibf.o meteogram.o mem_teb.o ModOdaRead.o tuvParameter.o \
	ModCuRead.o mem_aer1.o aer1_list.o mem_leaf.o ModRnode.o io_params.o \
	ModMicInit.o isan_coms.o VarTable.o ModMicrophysicsMisc.o dam.o ModCondRead.o \
	grid_dims.o chem_sources.o mem_stilt.o extra.o ModRanlavg.o ModNestIntrp.o \
	mem_emiss.o ModNudRead.o ModTimestepRK.o ModMemAlloc.o ModMkSfcTop.o \
	ModInitHis.o mem_chem1aq.o ModRThrm.o ModParaInit.o mem_scratch.o \
	ModNestGeoSst.o ModGrid.o mem_grell_param2.o ModRUser.o \
	ModParallelEnvironment.o machine_arq.o mem_teb_common.o ModSoilMoisture.o \
	modIau.o ModDomainDecomp.o dump.o ModSstRead.o ModBasicFields.o mem_grid.o \
	local_proc.o ModInitMicThompson.o ModRamsGrid.o shcu_vars_const.o ModGasPart.o \
	ModMonotonicAdvection.o ModTuvDriver2.7.o ModMPassDtl.o ModSched.o \
	ModEvaluation.o ModMkSfcFuso.o chem1_list.o ModLeaf3Teb.o memSoilMoisture.o \
	ModAerClim.o mem_teb_vars_const.o mem_varinit.o mem_carma.o ModChemAsgen.o \
	ModCuParGrell3.o ModRecycle.o ModMkSfcSfc.o ModRio.o ModCoriolis.o \
	ModWindFarm.o ModNamelistFile.o ModPostGridNetCDF.o mem_plume_chem1.o \
	domain_decomp.o CuParmVars.o ccatt_start.o ModNdviRead.o ModTimeStamp.o \
	ModChemistryDriver.o ModPostProcess.o chem_dry_dep.o ModGridTree.o \
	teb_spm_start.o ModTuv2.7.o ReadBcst.o ModVarfFile.o ModLeaf3Init.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/tsNames.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOpspec.o : $(IO)/ModOpspec.f90 mem_leaf.o io_params.o mem_varinit.o \
	mem_grell_param2.o ModMicControl.o ModNamelistFile.o modPrintInitial.o \
	grid_dims.o modIau.o mem_chem1aq.o chem_sources.o mem_grid.o mem_stilt.o \
	mem_globrad.o mem_emiss.o ccatt_start.o mem_chem1.o shcu_vars_const.o \
	teb_spm_start.o chem1_list.o mem_aer1.o aer1_list.o chem1aq_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOptical.o : $(RADIATE)/ModOptical.f90 mem_leaf.o mem_grid.o ModNamelistFile.o \
	ModTurbFields.o ModMPassFull.o VarTable.o ModControlVars.o ccatt_start.o \
	node_mod.o ModSoilMoisture.o ReadBcst.o ModRamsGrid.o mem_aer1.o parlibf.o \
	dump.o aer1_list.o ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOutputUtils.o : $(IO)/ModOutputUtils.f90 ModTurbFields.o ModNamelistFile.o \
	VarTable.o dump.o ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOzone.o : $(TEB_SPM)/ModOzone.f90 mem_grid.o RadiateFields.o rconstants.o \
	GaspartFields.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParaInit.o : $(MPI)/ModParaInit.f90 mem_grid.o ModScalarTable.o VarTable.o \
	grid_dims.o node_mod.o dump.o $(UTILS_INCS)/constants.h 
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

ModPostGrid.o : $(POST_SRC)/ModPostGrid.F90 mem_grid.o ModOutputUtils.o \
	ModNamelistFile.o ModTurbFields.o io_params.o ModBramsGrid.o \
	ModParallelEnvironment.o VarTable.o ModPostOneFieldNetCDF.o ModPostUtils.o \
	parlibf.o ModPostTypes.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostGridNetCDF.o : $(POST_SRC)/ModPostGridNetCDF.F90 mem_grid.o \
	ModNamelistFile.o io_params.o ModBramsGrid.o ModPostUtils.o ModDateUtils.o \
	ModPostTypes.o dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField.o : $(POST_SRC)/ModPostOneField.f90 ModTurbFields.o \
	ModNamelistFile.o ModBramsGrid.o VarTable.o ModPostOneField2d.o ModPostUtils.o \
	node_mod.o ModPostOneFieldUtils.o ModPostOneField8d.o ModPostOneField3d.o \
	ModPostOneField7d.o dump.o ModPostTypes.o ModMicControl.o ModBasicFields.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField2d.o : $(POST_SRC)/ModPostOneField2d.f90 mem_grid.o \
	ModOutputUtils.o ModNamelistFile.o ModTurbFields.o io_params.o mem_aerad.o \
	ModBramsGrid.o VarTable.o ModPostUtils.o node_mod.o ModPostOneFieldUtils.o \
	dump.o CuParmFields.o ModPostTypes.o ModMicControl.o ModPostGrid.o \
	ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField3d.o : $(POST_SRC)/ModPostOneField3d.f90 mem_grid.o \
	ModOutputUtils.o ModNamelistFile.o ModTurbFields.o ModBramsGrid.o VarTable.o \
	mem_varinit.o ModPostUtils.o node_mod.o ModPostOneFieldUtils.o ModPostTypes.o \
	ModMicControl.o ModPostGrid.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField7d.o : $(POST_SRC)/ModPostOneField7d.f90 ModOutputUtils.o \
	ModTurbFields.o ModNamelistFile.o ModBramsGrid.o VarTable.o ModPostUtils.o \
	ModPostOneFieldUtils.o ModPostTypes.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField8d.o : $(POST_SRC)/ModPostOneField8d.f90 ModOutputUtils.o \
	ModTurbFields.o ModNamelistFile.o ModBramsGrid.o VarTable.o ModPostUtils.o \
	ModPostOneFieldUtils.o ModPostTypes.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneFieldNetCDF.o : $(POST_SRC)/ModPostOneFieldNetCDF.F90 \
	ModPostGridNetCDF.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneFieldUtils.o : $(POST_SRC)/ModPostOneFieldUtils.f90 ModPostTypes.o \
	ModPostGrid.o ModBramsGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostProcess.o : $(POST_SRC)/ModPostProcess.F90 ModTurbFields.o io_params.o \
	ModNamelistFile.o ModBramsGrid.o ModParallelEnvironment.o VarTable.o \
	ModPostGrid.o ModPostGridNetCDF.o ModBasicFields.o ModPostOneField.o ModGrid.o \
	ModPostTypes.o ModGridTree.o ModMessageSet.o $(UTILS_INCS)/tsNames.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostTypes.o : $(POST_SRC)/ModPostTypes.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostUtils.o : $(POST_SRC)/ModPostUtils.f90 ModParallelEnvironment.o dump.o \
	mem_leaf.o $(UTILS_INCS)/files.h $(POST_INCS)/post_rconstants.h \
	$(UTILS_INCS)/constants.h $(POST_INCS)/post_rconfig.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modPrintInitial.o : $(INIT)/modPrintInitial.F90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvc.o : $(MODEL)/ModRadvc.f90 mem_grid.o ModNamelistFile.o \
	ModMonotonicAdvection.o ModBasicFields.o ModScalarTable.o \
	ModParallelEnvironment.o mem_tend.o ModRadvcAdap.o ccatt_start.o mem_chem1.o \
	grid_dims.o mem_aer1.o chem_dry_dep.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcRK.o : $(MODEL)/ModRadvcRK.f90 mem_grid.o mem_stilt.o \
	ModParallelEnvironment.o mem_tend.o mem_chem1.o ModGrid.o grid_dims.o \
	node_mod.o ModRexev.o ModMessageSet.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsReadHeader.o : $(IO)/ModRamsReadHeader.f90 an_header.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRanlavg.o : $(IO)/ModRanlavg.f90 mem_grid.o io_params.o VarTable.o \
	grid_dims.o ModMicroFields.o ModMicControl.o ModRThrm.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRbnd.o : $(BC)/ModRbnd.f90 mem_grid.o ModTurbFields.o mem_scratch.o \
	ModScalarTable.o ref_sounding.o mem_tend.o ModMicrophysicsMisc.o \
	ModMicroFields.o ccatt_start.o mem_chem1.o node_mod.o ModTurbKE.o \
	ModMicControl.o ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRcio.o : $(IO)/ModRcio.f90 mem_grid.o mem_leaf.o mem_stilt.o \
	ModNamelistFile.o io_params.o ref_sounding.o ModParallelEnvironment.o \
	grid_dims.o ModMicControl.o an_header.o ModLeafComs.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRecycle.o : $(IO)/ModRecycle.f90 mem_grid.o io_params.o mem_aerad.o \
	ModMPassFull.o VarTable.o ModRamsReadHeader.o mem_chem1.o mem_aer1.o node_mod.o \
	chem1_list.o ReadBcst.o ModGetVar.o ModDateUtils.o dump.o an_header.o \
	aer1_list.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRexev.o : $(STILT)/ModRexev.f90 mem_grid.o mem_stilt.o ModBasicFields.o \
	mem_tend.o ModRadvc.o rconstants.o ModMicControl.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRGrad.o : $(TURB)/ModRGrad.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRhhi.o : $(INIT)/ModRhhi.f90 mem_grid.o ModBasicFields.o ref_sounding.o \
	ModRinit.o grid_dims.o rconstants.o ModMicControl.o ModRamsGrid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRinit.o : $(INIT)/ModRinit.f90 mem_grid.o ModRbnd.o ModTurbFields.o \
	io_params.o mem_scratch.o ref_sounding.o mem_varinit.o ModTurbKE.o node_mod.o \
	rconstants.o ModMicroFields.o ModMicControl.o ModBasicFields.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRio.o : $(IO)/ModRio.f90 ModTurbFields.o io_params.o node_mod.o \
	ModMicControl.o ModNamelistFile.o mem_aerad.o VarTable.o ModControlVars.o \
	ModParallelEnvironment.o grid_dims.o utilsMod.o ModBasicFields.o mem_grid.o \
	ModMPassFull.o ref_sounding.o mpi_io_engine-5d.o mem_chem1.o parlibf.o \
	ReadBcst.o ModDateUtils.o an_header.o ModRcio.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/interface.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRrtmDriver.o : $(RADIATE)/ModRrtmDriver.f90 mem_leaf.o rrtmg_sw_cldprop.o \
	RadiateFields.o mem_tuv.o parkind.o mem_carma.o mcica_subcol_gen_lw.o \
	node_mod.o rconstants.o mem_rrtm.o ModMicControl.o ModNamelistFile.o \
	mcica_subcol_gen_sw.o grid_dims.o CuParmFields.o rrtmg_lw_rad.o \
	ModBasicFields.o mem_grid.o ref_sounding.o mem_tend.o mem_scratch1_grell.o \
	ccatt_start.o mem_chem1.o ModOptical.o rrtmg_sw_rad.o parrrsw.o teb_spm_start.o \
	parrrtm.o rrtmg_lw_cldprop.o ModMicroFields.o ModDateUtils.o ModLeafComs.o \
	$(UTILS_INCS)/aerosol_setup.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRShCuPar.o : $(CUPARM)/ModRShCuPar.f90 mem_grid.o shcu_vars_const.o \
	ModTurbFields.o ModConvComs.o mem_tend.o node_mod.o rconstants.o ShcuFields.o \
	ModMicroFields.o ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRstilt.o : $(STILT)/ModRstilt.f90 mem_stilt.o mem_grid.o ModTurbFields.o \
	ModMonotonicAdvection.o mem_scratch.o mem_scratch1_grell.o grid_dims.o \
	ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRThrm.o : $(MODEL)/ModRThrm.f90 mem_grid.o mem_scratch.o rconstants.o \
	ModMicroFields.o ModMicControl.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRtimi.o : $(MODEL)/ModRtimi.f90 mem_grid.o shcu_vars_const.o ModBasicFields.o \
	ModScalarTable.o mem_tend.o node_mod.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModScalarTable.o : $(MEMORY)/ModScalarTable.f90 ModParallelEnvironment.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSched.o : $(MODEL)/ModSched.f90 mem_grid.o local_proc.o shcu_vars_const.o \
	ModNamelistFile.o io_params.o ref_sounding.o ModParallelEnvironment.o \
	isan_coms.o mem_varinit.o node_mod.o ReadBcst.o parlibf.o dump.o \
	ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSeaSalt.o : $(CCATT)/ModSeaSalt.f90 mem_leaf.o mem_grid.o io_params.o \
	ModAerClim.o ccatt_start.o mem_chem1.o node_mod.o mem_aer1.o aer1_list.o \
	ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSstRead.o : $(MKSFC)/ModSstRead.f90 mem_grid.o mem_leaf.o io_params.o \
	ModMkSfcSst.o ModControlVars.o node_mod.o grid_dims.o ReadBcst.o ModDateUtils.o \
	mem_mksfc.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimeStamp.o : $(MODEL)/ModTimeStamp.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestep.o : $(MODEL)/ModTimestep.F90 mem_leaf.o ModRShCuPar.o ModRtimi.o \
	ModMatrixDriver.o mem_varinit.o ModOdaNudge.o ModUrbanCanopy.o ModCuParGrell3.o \
	ModGrid.o node_mod.o rconstants.o ModMicThompsonDriver.o ModDiffuse.o \
	ModCoriolis.o ModAcoust.o ModWindFarm.o ModSeaSalt.o ModRexev.o digitalFilter.o \
	ModMicrophysicsMisc.o mem_oda.o ModRadvc.o ModNudAnalysis.o grid_dims.o \
	machine_arq.o ModRstilt.o mem_plume_chem1.o ChemDryDepDriver.o \
	ChemSourcesDriver.o ModBasicFields.o mem_grid.o mem_stilt.o chem_sources.o \
	sfclyr_jules.o mem_tend.o mem_emiss.o ccatt_start.o mem_chem1.o ModOzone.o \
	ModLeaf3.o ModChemistryDriver.o ModTimeStamp.o ModOptical.o ModTurbK.o \
	ModMessageSet.o shcu_vars_const.o ModMicrophysicsDrive.o ModGasPart.o \
	ModMonotonicAdvection.o ModRbnd.o teb_spm_start.o ModMicGfdlDriver.o rad_driv.o \
	mem_aer1.o ModRThrm.o mem_scratch.o $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestepRK.o : $(MODEL)/ModTimestepRK.F90 mem_leaf.o ModRShCuPar.o ModRtimi.o \
	ModAerClim.o ModMatrixDriver.o mem_varinit.o ModOdaNudge.o ModUrbanCanopy.o \
	ModCuParGrell3.o ModGrid.o node_mod.o rconstants.o ModMicThompsonDriver.o \
	ModDiffuse.o ModCoriolis.o ModAcoust.o ModWindFarm.o ModSeaSalt.o ModRexev.o \
	ModTimestep.o digitalFilter.o ModParallelEnvironment.o ModMicrophysicsMisc.o \
	mem_oda.o ModRadvc.o ModNudAnalysis.o grid_dims.o machine_arq.o ModRstilt.o \
	mem_plume_chem1.o modIau.o ChemDryDepDriver.o ChemSourcesDriver.o utilsMod.o \
	chem_sources.o mem_grid.o mem_stilt.o sfclyr_jules.o mem_tend.o mem_emiss.o \
	ccatt_start.o mem_chem1.o ModOzone.o ModLeaf3.o ModChemistryDriver.o \
	ModTimeStamp.o ModOptical.o ModTurbK.o ModMessageSet.o ModRbnd.o \
	ModMicrophysicsDrive.o ModGasPart.o shcu_vars_const.o ModMonotonicAdvection.o \
	ModRadvcRK.o teb_spm_start.o ModMicGfdlDriver.o ModLeaf3OceanOnly.o rad_driv.o \
	mem_aer1.o ModRThrm.o mem_scratch.o $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiff.o : $(TURB)/ModTurbDiff.f90 mem_grid.o mem_opt_scratch.o ModRGrad.o \
	mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiffAdap.o : $(TURB)/ModTurbDiffAdap.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbFields.o : $(TURB)/ModTurbFields.f90 ModNamelistFile.o \
	ModParallelEnvironment.o VarTable.o ModNodeDimensions.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbK.o : $(TURB)/ModTurbK.f90 mem_leaf.o ModTurbFields.o ModMicroFields.o \
	ModRGrad.o node_mod.o rconstants.o ModMicControl.o ModTKenn.o ModNamelistFile.o \
	grid_dims.o ModRstilt.o ModBasicFields.o mem_grid.o mem_stilt.o ModTurbKAdap.o \
	ModScalarTable.o mem_turb_scalar.o mem_tend.o mem_grell.o ke_coms.o mem_chem1.o \
	ccatt_start.o ModMonotonicAdvection.o ModTurbDiff.o ModTurbDiffAdap.o \
	ModTurbKE.o mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKAdap.o : $(TURB)/ModTurbKAdap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKE.o : $(TURB)/ModTurbKE.f90 mem_grid.o ModTurbFields.o ke_coms.o \
	rconstants.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuv2.7.o : $(TUV)/ModTuv2.7.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuvDriver2.7.o : $(TUV)/ModTuvDriver2.7.f90 mem_grid.o mem_leaf.o \
	RadiateFields.o ModNamelistFile.o mem_aerad.o mem_globrad.o ref_sounding.o \
	mem_tuv.o mem_carma.o extra.o ModTuv2.7.o tuvParameter.o mem_chem1.o \
	chem_fastjx_driv.o node_mod.o rconstants.o chem1_list.o mem_rrtm.o \
	ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_cu_g3.o : $(CUPARM)/module_cu_g3.f90 
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

ModUrbanCanopy.o : $(SURFACE)/ModUrbanCanopy.f90 mem_grid.o ModTurbFields.o \
	mem_tend.o node_mod.o ModBasicFields.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfFile.o : $(FDDA)/ModVarfFile.f90 mem_leaf.o mem_varinit.o ModRcio.o \
	ModRamsReadHeader.o ModGrid.o node_mod.o rconstants.o ModMicControl.o \
	isan_coms.o ModControlVars.o ModNudAnalysis.o ModVarfUpdate.o ModGetVar.o \
	ModBasicFields.o mem_grid.o ref_sounding.o mem_chem1.o parlibf.o ModRamsGrid.o \
	ModGridTree.o ModMessageSet.o ModDateUtils.o chem1_list.o ReadBcst.o mem_aer1.o \
	aer1_list.o mem_scratch.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModWindFarm.o : $(WIND_FARM)/ModWindFarm.f90 mem_grid.o ModTurbFields.o \
	io_params.o ModNamelistFile.o mem_tend.o node_mod.o rconstants.o \
	module_wind_fitch.o ModDateUtils.o ModBasicFields.o $(UTILS_INCS)/files.h 
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

ModNestIntrp.o : $(NESTING)/ModNestIntrp.f90 mem_grid.o ModBasicFields.o \
	ref_sounding.o ModRinit.o grid_dims.o ModNestFillDens.o rconstants.o \
	mem_nestb.o mem_scratch.o 
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

ModNudRead.o : $(FDDA)/ModNudRead.f90 mem_grid.o isan_coms.o VarTable.o \
	mem_varinit.o mem_chem1.o ModNudAnalysis.o ModNudUpdate.o ModDateUtils.o \
	ModRamsGrid.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudUpdate.o : $(FDDA)/ModNudUpdate.f90 mem_grid.o mem_aerad.o VarTable.o \
	mem_varinit.o grid_struct.o mem_chem1.o chem1_list.o ModInitHis.o an_header.o \
	ModRcio.o $(UTILS_INCS)/files.h 
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

ModOdaProcObs.o : $(FDDA)/ModOdaProcObs.f90 mem_grid.o rconstants.o mem_oda.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaRead.o : $(FDDA)/ModOdaRead.f90 mem_grid.o isan_coms.o ModOdaStaCount.o \
	mem_oda.o ModOdaStaInput.o ModDateUtils.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaCount.o : $(FDDA)/ModOdaStaCount.f90 mem_grid.o ModReadRalph.o \
	mem_oda.o obs_input.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaInput.o : $(FDDA)/ModOdaStaInput.f90 obs_input.o mem_grid.o \
	ModReadRalph.o ModOdaStaCount.o mem_oda.o ModDateUtils.o 
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

ModAcoustAdap.o : $(MODEL)/ModAcoustAdap.f90 mem_grid.o ModRbnd.o ModGrid.o \
	node_mod.o rconstants.o ModMessageSet.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_carma.o : $(RADIATE)/rad_carma.F90 mem_grid.o mem_leaf.o ModNamelistFile.o \
	mem_tuv.o mem_aerad.o mem_globrad.o mem_carma.o ccatt_start.o carma_fastjx.o \
	mem_chem1.o mem_aer1.o grid_dims.o rconstants.o mem_globaer.o chem1_list.o \
	machine_arq.o node_mod.o ModDateUtils.o aer1_list.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_driv.o : $(RADIATE)/rad_driv.f90 RadiateFields.o ModNamelistFile.o \
	ModCarmaDriver.o ModRrtmDriver.o ModMicroFields.o CuParmFields.o \
	ModMicControl.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcAdap.o : $(MODEL)/ModRadvcAdap.f90 ModAdapInit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsGrid.o : $(INIT)/ModRamsGrid.f90 mem_grid.o ModAdapInit.o node_mod.o \
	rconstants.o ModGridSet.o dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rconstants.o : $(MEMORY)/rconstants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModReadRalph.o : $(FDDA)/ModReadRalph.f90 obs_input.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ReadBcst.o : $(MPI)/ReadBcst.f90 mem_grid.o ModTurbFields.o mem_aerad.o \
	ModMPassFull.o ModControlVars.o utilsMod.o node_mod.o parlibf.o an_header.o \
	ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ref_sounding.o : $(MODEL)/ref_sounding.f90 grid_dims.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRnode.o : $(MODEL)/ModRnode.f90 mem_grid.o mem_leaf.o VarTable.o grid_dims.o \
	node_mod.o parlibf.o $(UTILS_INCS)/constants.h 
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

rrtmg_lw_cldprmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprmc.f90 parkind.o parrrtm.o \
	rrlw_wvn.o rrlw_cld.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprop.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprop.f90 parkind.o rrlw_cld.o \
	rrlw_vsn.o parrrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_init.o : $(RRTMG_LW_SRC)/rrtmg_lw_init.f90 rrlw_kg14.o parkind.o \
	rrlw_kg05.o rrlw_kg06.o rrlw_kg08.o rrlw_kg07.o rrlw_kg03.o rrtmg_lw_setcoef.o \
	rrtmg_lw_k_g.o rrlw_kg11.o rrlw_kg04.o rrlw_kg15.o rrlw_wvn.o rrlw_kg10.o \
	rrlw_cld.o rrlw_vsn.o rrlw_kg01.o rrlw_kg16.o rrlw_kg12.o rrlw_kg02.o \
	rrlw_tbl.o rrlw_con.o parrrtm.o rrlw_kg09.o rrlw_kg13.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_k_g.o : $(RRTMG_LW_SRC)/rrtmg_lw_k_g.f90 rrlw_kg01.o rrlw_kg14.o \
	parkind.o rrlw_kg05.o rrlw_kg16.o rrlw_kg11.o rrlw_kg12.o rrlw_kg04.o \
	rrlw_kg09.o rrlw_kg15.o rrlw_kg06.o rrlw_kg10.o rrlw_kg08.o rrlw_kg07.o \
	rrlw_kg03.o rrlw_kg13.o rrlw_vsn.o rrlw_kg02.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rad.o : $(RRTMG_LW_SRC)/rrtmg_lw_rad.f90 rrtmg_lw_cldprmc.o \
	rrtmg_lw_setcoef.o parkind.o rrlw_con.o mcica_subcol_gen_lw.o parrrtm.o \
	rrlw_wvn.o rrtmg_lw_taumol.o rrtmg_lw_rtrnmc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrn.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrn.f90 parkind.o rrlw_con.o \
	rrlw_tbl.o parrrtm.o rrlw_wvn.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmc.f90 parkind.o rrlw_con.o \
	rrlw_tbl.o parrrtm.o rrlw_wvn.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmr.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmr.f90 parkind.o rrlw_con.o \
	rrlw_tbl.o parrrtm.o rrlw_wvn.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_setcoef.o : $(RRTMG_LW_SRC)/rrtmg_lw_setcoef.f90 parkind.o parrrtm.o \
	rrlw_wvn.o rrlw_ref.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_taumol.o : $(RRTMG_LW_SRC)/rrtmg_lw_taumol.f90 rrlw_kg14.o parkind.o \
	rrlw_kg05.o rrlw_kg06.o rrlw_kg08.o rrlw_kg07.o rrlw_kg03.o rrlw_kg11.o \
	rrlw_kg04.o rrlw_kg15.o rrlw_wvn.o rrlw_kg10.o rrlw_vsn.o rrlw_kg01.o \
	rrlw_kg16.o rrlw_kg12.o rrlw_kg02.o rrlw_con.o parrrtm.o rrlw_kg09.o rrlw_ref.o \
	rrlw_kg13.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprmc.f90 parrrsw.o parkind.o \
	rrsw_wvn.o rrsw_vsn.o rrsw_cld.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprop.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprop.f90 parrrsw.o parkind.o \
	rrsw_wvn.o rrsw_vsn.o rrsw_cld.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_init.o : $(RRTMG_SW_SRC)/rrtmg_sw_init.f90 rrsw_kg16.o parkind.o \
	rrsw_kg23.o rrsw_wvn.o rrsw_kg29.o rrsw_kg17.o rrsw_kg27.o rrsw_kg25.o \
	rrsw_kg22.o rrsw_kg18.o rrsw_kg26.o rrsw_cld.o rrsw_aer.o rrsw_kg24.o \
	rrsw_kg28.o rrtmg_sw_setcoef.o parrrsw.o rrtmg_sw_k_g.o rrsw_kg21.o rrsw_tbl.o \
	rrsw_kg20.o rrsw_con.o rrsw_kg19.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_k_g.o : $(RRTMG_SW_SRC)/rrtmg_sw_k_g.f90 rrsw_kg27.o rrsw_kg25.o \
	parkind.o rrsw_kg16.o rrsw_kg21.o rrsw_kg24.o rrsw_kg20.o rrsw_kg28.o \
	rrsw_kg22.o rrsw_kg23.o rrsw_kg19.o rrsw_kg29.o rrsw_kg18.o rrsw_kg26.o \
	rrsw_kg17.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_rad.o : $(RRTMG_SW_SRC)/rrtmg_sw_rad.f90 parrrsw.o \
	mcica_subcol_gen_sw.o parkind.o rrsw_aer.o rrtmg_sw_cldprmc.o \
	rrtmg_sw_setcoef.o rrsw_con.o rrsw_wvn.o rrtmg_sw_spcvmc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_reftra.o : $(RRTMG_SW_SRC)/rrtmg_sw_reftra.f90 rrsw_vsn.o rrsw_tbl.o \
	parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_setcoef.o : $(RRTMG_SW_SRC)/rrtmg_sw_setcoef.f90 rrsw_vsn.o parrrsw.o \
	rrsw_ref.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvmc.f90 rrtmg_sw_reftra.o \
	parrrsw.o parkind.o rrtmg_sw_taumol.o rrsw_tbl.o rrtmg_sw_vrtqdr.o rrsw_wvn.o \
	rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvrt.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvrt.f90 rrtmg_sw_reftra.o \
	parrrsw.o parkind.o rrtmg_sw_taumol.o rrsw_tbl.o rrtmg_sw_vrtqdr.o rrsw_wvn.o \
	rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_taumol.o : $(RRTMG_SW_SRC)/rrtmg_sw_taumol.f90 rrsw_kg27.o parrrsw.o \
	rrsw_kg25.o rrsw_kg16.o parkind.o rrsw_kg21.o rrsw_kg24.o rrsw_kg20.o \
	rrsw_kg28.o rrsw_kg22.o rrsw_kg23.o rrsw_con.o rrsw_wvn.o rrsw_kg19.o \
	rrsw_kg29.o rrsw_kg18.o rrsw_kg26.o rrsw_kg17.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_vrtqdr.o : $(RRTMG_SW_SRC)/rrtmg_sw_vrtqdr.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRUser.o : $(SURFACE)/ModRUser.f90 mem_grid.o memSoilMoisture.o mem_leaf.o \
	io_params.o ccatt_start.o node_mod.o rconstants.o ModLeafComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

setup.o : $(MATRIX)/setup.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

sfclyr_jules.o : $(JULES_DIR)/sfclyr_jules.f90 mem_leaf.o RadiateFields.o \
	io_params.o model_time_mod.o ModTurbFields.o JulesFields.o node_mod.o \
	rconstants.o ModMicControl.o jules_fields_mod.o jules_surface_types_mod.o \
	ModNamelistFile.o csigma_mod.o CuParmFields.o sf_diags_mod.o ModBasicFields.o \
	mem_grid.o ModLeaf3Init.o mem_chem1.o gridmean_fluxes.o datetime_mod.o \
	gridbox_mean_mod.o ancil_info.o io_constants.o chem1_list.o mem_brams_jules.o \
	ModMicroFields.o fluxes.o ModLeafComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

shcu_vars_const.o : $(CUPARM)/shcu_vars_const.f90 grid_dims.o ModNamelistFile.o \
	ModConvComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSoilMoisture.o : $(SOIL_MOISTURE)/ModSoilMoisture.F90 memSoilMoisture.o \
	mem_grid.o mem_leaf.o ModNamelistFile.o ModTurbFields.o io_params.o mem_aerad.o \
	ModMPassFull.o ModControlVars.o ModBasicFields.o node_mod.o rconstants.o \
	ReadBcst.o parlibf.o ModLeafComs.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

solut.o : $(MATRIX)/solut.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

subs.o : $(MATRIX)/subs.f90 setup.o memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

teb_spm_start.o : $(TEB_SPM)/teb_spm_start.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTKenn.o : $(STILT)/ModTKenn.f90 mem_grid.o mem_stilt.o rconstants.o \
	turb_constants.o mem_scratch.o 
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

VarTable.o : $(MEMORY)/VarTable.f90 ModParallelEnvironment.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfUpdate.o : $(FDDA)/ModVarfUpdate.f90 mem_grid.o ref_sounding.o \
	rconstants.o ModInitHis.o mem_scratch.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

include jules_depend_model.mk

