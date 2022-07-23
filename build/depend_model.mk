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

utils_f.o : $(UTILS_LIB)/utils_f.f90 dump.o ModDateUtils.o $(UTILS_INCS)/files.h \
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

ModAsTi.o : $(ISAN)/ModAsTi.f90 isan_coms.o rconstants.o ModChemAObj.o \
	mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsTp.o : $(ISAN)/ModAsTp.f90 isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAVarF.o : $(ISAN)/ModAVarF.f90 ModRbnd.o isan_coms.o rconstants.o mem_grid.o 
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

ModChemAsti.o : $(ISAN_CHEM)/ModChemAsti.f90 mem_chem1.o ModAsTi.o \
	ModChemAsti2.o ModChemVInterps.o mem_aer1.o chem_isan_coms.o mem_grid.o \
	ModChemFirstRams.o ModChemAObj.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsti2.o : $(ISAN_CHEM)/ModChemAsti2.f90 isan_coms.o ModDateUtils.o \
	rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAstp.o : $(ISAN_CHEM)/ModChemAstp.F90 mem_chem1.o chem_isan_coms.o dump.o \
	rconstants.o ModDateUtils.o chem1_list.o isan_coms.o ModAsTp.o mem_varinit.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAvarf.o : $(ISAN_CHEM)/ModChemAvarf.f90 mem_chem1.o ModAVarF.o mem_aer1.o \
	mem_grid.o chem_isan_coms.o ModRbnd.o rconstants.o isan_coms.o ModNestFeed.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_dry_dep.o : $(MODEL_CHEM)/chem_dry_dep.f90 mem_chem1.o mem_aer1.o \
	ModDateUtils.o chem1_list.o extra.o aer1_list.o 
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

chem_fastjx_driv.o : $(CCATT)/chem_fastjx_driv.f90 chem1_list.o \
	chem_fastjx_data.o chem_fastjx57.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFileInv.o : $(ISAN_CHEM)/ModChemFileInv.f90 isan_coms.o ModDateUtils.o \
	dump.o mem_grid.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFirstRams.o : $(ISAN_CHEM)/ModChemFirstRams.f90 grid_dims.o an_header.o \
	mem_scratch.o ModRcio.o ModChemRefState.o mem_grid.o ModGetVar.o rconstants.o \
	isan_coms.o ModNestFillDens.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_isan_coms.o : $(ISAN_CHEM)/chem_isan_coms.f90 chem1_list.o isan_coms.o \
	aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemIsanIo.o : $(ISAN_CHEM)/ModChemIsanIo.f90 isan_coms.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_orage.o : $(CCATT)/chem_orage.f90 mem_scratch1_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_plumerise_scalar.o : $(CCATT)/chem_plumerise_scalar.f90 mem_aer1.o \
	node_mod.o mem_chem1.o mem_plume_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemRefState.o : $(ISAN_CHEM)/ModChemRefState.f90 ccatt_start.o \
	ModNestFillDens.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_sources.o : $(CCATT)/chem_sources.f90 ModNamelistFile.o ModControlVars.o \
	mem_chem1.o ReadBcst.o mem_plume_chem1.o mem_aer1.o mem_grid.o io_params.o \
	ModDateUtils.o parlibf.o aer1_list.o mem_volc_chem1.o $(UTILS_INCS)/constants.h 
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

chem_spack_qssa.o : $(CCATT)/chem_spack_qssa.f90 chem_spack_fexloss.o \
	chem_spack_dratedc.o chem_spack_kinetic.o mem_chem1.o chem_spack_fexprod.o \
	chem_spack_rates.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rates.o : $(MODEL_CHEM)/chem_spack_rates.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rodas3_dyndt.o : $(CCATT)/chem_spack_rodas3_dyndt.f90 \
	chem_spack_kinetic.o mem_chem1.o chem_spack_jacdchemdc.o mem_grid.o mem_spack.o \
	chem_spack_fexchem.o extra.o chem_spack_ros.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros.o : $(CCATT)/chem_spack_ros.f90 chem_spack_kinetic.o mem_chem1.o \
	chem_spack_jacdchemdc.o mem_spack.o chem_spack_fexchem.o \
	chem_spack_solve_sparse.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros_dyndt.o : $(CCATT)/chem_spack_ros_dyndt.f90 chem_spack_kinetic.o \
	mem_chem1.o chem_spack_jacdchemdc.o mem_spack.o chem_spack_fexchem.o \
	chem_spack_ros.o chem_spack_solve_sparse.o 
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

ModChemVInterps.o : $(ISAN_CHEM)/ModChemVInterps.f90 isan_coms.o rconstants.o \
	dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemDryDepDriver.o : $(MODEL_CHEM)/ChemDryDepDriver.f90 grid_dims.o \
	ModMicControl.o mem_chem1.o ModMicroFields.o mem_aer1.o mem_grid.o \
	ModTurbFields.o mem_cuparm.o rconstants.o mem_leaf.o chem_dry_dep.o \
	mem_radiate.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemSourcesDriver.o : $(CCATT)/ChemSourcesDriver.f90 mem_chem1.o mem_stilt.o \
	mem_plume_chem1.o mem_aer1.o io_params.o chem_sources.o mem_leaf.o chem1_list.o \
	aer1_list.o chem_plumerise_scalar.o mem_volc_chem1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

coag.o : $(MATRIX)/coag.f90 memMatrix.o setup.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModCondRead.o : $(FDDA)/ModCondRead.f90 ModRamsGrid.o mem_grid.o ModCondUpdate.o \
	ModDateUtils.o ModNudAnalysis.o isan_coms.o mem_varinit.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCondUpdate.o : $(FDDA)/ModCondUpdate.f90 an_header.o ModRcio.o ModInitHis.o \
	mem_grid.o grid_struct.o ModVarTables.o mem_varinit.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModConvComs.o : $(CUPARM)/ModConvComs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ConvPar_GF_GEOS5.o : $(CUPARM)/ConvPar_GF_GEOS5.F90 MAPL_Constants.o \
	Henrys_Law_cts.o module_gate.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuRead.o : $(CUPARM)/ModCuRead.f90 mem_grid.o ModDateUtils.o isan_coms.o \
	mem_cuparm.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCupDn.o : $(CUPARM)/ModCupDn.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCupEnv.o : $(CUPARM)/ModCupEnv.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCupEnvCatt.o : $(CUPARM)/ModCupEnvCatt.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCupGrellCattDeep.o : $(CUPARM)/ModCupGrellCattDeep.f90 node_mod.o ModCupEnv.o \
	mem_carma.o Phys_const.o ModCupDn.o mem_grell_param2.o ccatt_start.o mem_grid.o \
	ModCupEnvCatt.o cup_output_vars.o kbcon_ecmwf.o ModCupUp.o mem_scratch2_grell.o \
	mem_scratch3_grell.o mem_varinit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCupGrellCattShallow.o : $(CUPARM)/ModCupGrellCattShallow.f90 node_mod.o \
	ModCupEnv.o Phys_const.o mem_grell_param2.o mem_scratch3_grell_sh.o mem_grid.o \
	ModCupEnvCatt.o cup_output_vars.o ModCupUp.o mem_scratch2_grell_sh.o \
	mem_varinit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cup_output_vars.o : $(CUPARM)/cup_output_vars.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCupUp.o : $(CUPARM)/ModCupUp.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

dam.o : $(ENERGY)/dam.f90 ModNamelistFile.o mem_grid.o dump.o ModDateUtils.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

depv.o : $(MATRIX)/depv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

digitalFilter.o : $(MODEL)/digitalFilter.f90 ModNamelistFile.o grid_dims.o \
	node_mod.o ModControlVars.o ReadBcst.o mem_grid.o io_params.o ModDateUtils.o \
	utilsMod.o ModBasicFields.o ModVarTables.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

domain_decomp.o : $(INIT)/domain_decomp.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

extra.o : $(MEMORY)/extra.f90 ModNamelistFile.o ModVarTables.o dump.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGetVar.o : $(UTILS_LIB)/ModGetVar.f90 an_header.o dump.o \
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

ModGridSet.o : $(INIT)/ModGridSet.f90 grid_dims.o mem_grid.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

Henrys_Law_cts.o : $(CUPARM)/Henrys_Law_cts.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

io_params.o : $(IO)/io_params.f90 ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

isan_coms.o : $(ISAN_MODS)/isan_coms.f90 ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

isrpia.o : $(MATRIX)/isrpia.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

issoropia.o : $(MATRIX)/issoropia.f90 solut.o isrpia.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

kbcon_ecmwf.o : $(CUPARM)/kbcon_ecmwf.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ke_coms.o : $(TURB)/ke_coms.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Hyd.o : $(SURFACE)/ModLeaf3Hyd.f90 ModLeafComs.o mem_grid.o mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Init.o : $(SURFACE)/ModLeaf3Init.f90 teb_spm_start.o ModLeafComs.o \
	grid_dims.o ModLeaf3.o mem_grid.o io_params.o rconstants.o mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Teb.o : $(SURFACE)/ModLeaf3Teb.f90 ModUrban.o mem_emiss.o ModGasPart.o \
	mem_teb_vars_const.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeafComs.o : $(SURFACE)/ModLeafComs.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

local_proc.o : $(MODEL)/local_proc.F90 grid_dims.o node_mod.o ReadBcst.o \
	mem_stilt.o mem_grid.o io_params.o rconstants.o dump.o ref_sounding.o \
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

mcica_subcol_gen_lw.o : $(RRTMG_LW_SRC)/mcica_subcol_gen_lw.f90 \
	mcica_random_numbers.o parrrtm.o parkind.o rrlw_con.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_sw.o : $(RRTMG_SW_SRC)/mcica_subcol_gen_sw.f90 parrrsw.o \
	mcica_random_numbers.o parkind.o rrsw_vsn.o rrsw_con.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aer1.o : $(CCATT)/mem_aer1.f90 ModNamelistFile.o grid_dims.o node_mod.o \
	mem_chem1.o mem_grid.o io_params.o dump.o ModScalarTable.o aer1_list.o \
	ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aerad.o : $(RADIATE)/mem_aerad.f90 mem_grid_dim_defs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_carma.o : $(RADIATE)/mem_carma.f90 ModNamelistFile.o grid_dims.o node_mod.o \
	ModControlVars.o ReadBcst.o ModRamsGrid.o ModMPassFull.o mem_grid.o \
	ModSoilMoisture.o mem_scalar.o io_params.o ModTurbFields.o parlibf.o \
	mem_aerad.o ModBasicFields.o ModVarTables.o mem_globrad.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1.o : $(CCATT)/mem_chem1.f90 ModNamelistFile.o grid_dims.o io_params.o \
	ModScalarTable.o chem1_list.o ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1aq.o : $(CCATT)/mem_chem1aq.f90 ModNamelistFile.o grid_dims.o \
	chem1aq_list.o mem_chem1.o ModScalarTable.o ModVarTables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chemic.o : $(CCATT)/mem_chemic.f90 ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_cuparm.o : $(CUPARM)/mem_cuparm.f90 ModNamelistFile.o grid_dims.o \
	ModVarTables.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_cutrans.o : $(CUPARM)/mem_cutrans.f90 mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_emiss.o : $(TEB_SPM)/mem_emiss.f90 ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_gaspart.o : $(TEB_SPM)/mem_gaspart.f90 mem_emiss.o ModVarTables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globaer.o : $(RADIATE)/mem_globaer.f90 mem_precision.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globrad.o : $(RADIATE)/mem_globrad.f90 ModNamelistFile.o mem_radiate.o \
	mem_precision.o parlibf.o mem_aerad.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell.o : $(CUPARM)/mem_grell.f90 ModVarTables.o shcu_vars_const.o \
	mem_cuparm.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell_param2.o : $(CUPARM)/mem_grell_param2.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid.o : $(MEMORY)/mem_grid.f90 ModNamelistFile.o grid_dims.o ModVarTables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid_dim_defs.o : $(MEMORY)/mem_grid_dim_defs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_jules.o : $(SURFACE)/mem_jules.f90 ModNamelistFile.o grid_dims.o \
	ModVarTables.o io_params.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_leaf.o : $(SURFACE)/mem_leaf.f90 ModNamelistFile.o grid_dims.o \
	teb_spm_start.o io_params.o ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_micro.o : $(MICRO)/mem_micro.f90 ModMicControl.o ModMicroFields.o \
	mem_cuparm.o mem_radiate.o ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicroFields.o : $(MICRO)/ModMicroFields.f90 ModNamelistFile.o \
	ModParallelEnvironment.o ModNodeDimensions.o ModMicControl.o 
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

mem_oda.o : $(FDDA)/mem_oda.f90 ModNamelistFile.o grid_dims.o ModVarTables.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_opt_scratch.o : $(TURB)/mem_opt_scratch.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_plume_chem1.o : $(CCATT)/mem_plume_chem1.f90 ModNamelistFile.o \
	ModVarTables.o chem1_list.o mem_chem1.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_precision.o : $(RADIATE)/mem_precision.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_radiate.o : $(RADIATE)/mem_radiate.f90 ModNamelistFile.o ModVarTables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_rrtm.o : $(RADIATE)/mem_rrtm.f90 node_mod.o rrtmg_sw_init.o parrrsw.o \
	mem_chem1.o parrrtm.o parkind.o mem_grid.o rrtmg_lw_init.o rconstants.o \
	chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scalar.o : $(MEMORY)/mem_scalar.f90 ModNamelistFile.o ModVarTables.o \
	io_params.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch.o : $(MEMORY)/mem_scratch.f90 mem_aerad.o mem_radiate.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch1_brams.o : $(MEMORY)/mem_scratch1_brams.f90 
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

mem_shcu.o : $(CUPARM)/mem_shcu.f90 ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_spack.o : $(CCATT)/mem_spack.f90 chem1_list.o chem_spack_utils.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_stilt.o : $(STILT)/mem_stilt.f90 ModNamelistFile.o grid_dims.o io_params.o \
	rconstants.o ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tconv.o : $(CCATT)/mem_tconv.f90 chem1_list.o aer1_list.o mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_teb.o : $(TEB_SPM)/mem_teb.f90 ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_teb_common.o : $(TEB_SPM)/mem_teb_common.f90 ModVarTables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_teb_vars_const.o : $(TEB_SPM)/mem_teb_vars_const.f90 ModNamelistFile.o \
	grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tend.o : $(MEMORY)/mem_tend.f90 teb_spm_start.o mem_emiss.o \
	ModNamelistFile.o mem_grid.o ModTurbFields.o mem_scalar.o mem_gaspart.o \
	ModScalarTable.o ModBasicFields.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_turb_scalar.o : $(TURB)/mem_turb_scalar.f90 ModVarTables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tuv.o : $(TUV)/mem_tuv.f90 ModTuv2.7.o mem_globrad.o mem_stilt.o \
	ModVarTables.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mem_varinit.o : $(MEMORY)/mem_varinit.f90 ModNamelistFile.o grid_dims.o \
	mem_chem1.o chem1_list.o ModVarTables.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_volc_chem1.o : $(CCATT)/mem_volc_chem1.f90 ModNamelistFile.o ModVarTables.o \
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

meteogram.o : $(IO)/meteogram.f90 ModNamelistFile.o ModPostUtils.o node_mod.o \
	mem_grid.o ModMPassDtl.o meteogramType.o ModVarTables.o satPolyColision.o \
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

ModMicInit.o : $(MICRO)/ModMicInit.f90 ModMicControl.o mem_grid.o ModMicTabs.o \
	dump.o rconstants.o ModMicGamma.o $(MICRO)/MicConstants.h $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
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

ModMicControl.o : $(MICRO)/ModMicControl.f90 ModNamelistFile.o grid_dims.o \
	ModParallelEnvironment.o $(UTILS_INCS)/files.h $(MICRO)/MicConstants.h 
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

ModAcoust.o : $(MODEL)/ModAcoust.f90 mem_scratch.o node_mod.o ModMicControl.o \
	ModGrid.o ModMessageSet.o mem_grid.o mem_tend.o ModParallelEnvironment.o \
	ModBasicFields.o rconstants.o ModAcoustAdap.o ref_sounding.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAerClim.o : $(AERCLIM)/ModAerClim.f90 node_mod.o ModControlVars.o ReadBcst.o \
	mem_grid.o ModTurbFields.o ModSoilMoisture.o dump.o parlibf.o ModBasicFields.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBasicFields.o : $(MEMORY)/ModBasicFields.f90 ModNamelistFile.o mem_stilt.o \
	ModNodeDimensions.o ModParallelEnvironment.o ModVarTables.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBramsGrid.o : $(POST_SRC)/ModBramsGrid.f90 ModNamelistFile.o ModPostUtils.o \
	node_mod.o mem_grid.o ModParallelEnvironment.o mem_aerad.o ref_sounding.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBuffering.o : $(MPI)/ModBuffering.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCarmaDriver.o : $(RADIATE)/ModCarmaDriver.f90 teb_spm_start.o grid_dims.o \
	node_mod.o ModLeaf3.o mem_carma.o ModMicControl.o ModMicroFields.o mem_grid.o \
	mem_tend.o rad_carma.o mem_teb_common.o mem_cuparm.o ModDateUtils.o mem_leaf.o \
	rconstants.o mem_radiate.o ModBasicFields.o mem_scratch1_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsgen.o : $(ISAN_CHEM)/ModChemAsgen.F90 ModChemAvarf.o dump.o node_mod.o \
	ModChemAsti.o ModChemAstp.o mem_aer1.o ModChemFileInv.o ModDateUtils.o \
	ModChemIsanIo.o aer1_list.o ModMkSfcTop.o grid_dims.o mem_chem1.o \
	ModChemRefState.o chem_isan_coms.o io_params.o ModControlVars.o ModAsGen.o \
	ModRamsGrid.o mem_grid.o chem1_list.o isan_coms.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemConvTransp.o : $(CCATT)/ModChemConvTransp.f90 mem_scratch.o node_mod.o \
	mem_chem1.o Phys_const.o mem_grell_param2.o mem_aer1.o mem_grid.o mem_tconv.o \
	mem_cuparm.o chem1_list.o mem_scratch1_grell.o aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemistryDriver.o : $(CCATT)/ModChemistryDriver.f90 mem_cuparm.o \
	chem_fastjx_driv.o mem_radiate.o mem_aerad.o chem_uv_att.o ModBasicFields.o \
	mem_chemic.o chem_spack_rodas3_dyndt.o node_mod.o carma_fastjx.o mem_carma.o \
	ModMicroFields.o mem_aer1.o ModTuvDriver2.7.o mem_chem1aq.o mem_spack.o \
	chem_trans_liq.o chem_orage.o mem_scratch1_grell.o extra.o aer1_list.o \
	chem_spack_solve_sparse.o mem_globrad.o grid_dims.o mem_scratch.o mem_chem1.o \
	parrrtm.o mem_grell_param2.o chem_spack_qssa.o chem_trans_gasaq.o \
	chem_spack_ros.o mem_rrtm.o chem1aq_list.o mem_stilt.o mem_grid.o \
	chem_spack_ros_dyndt.o chem_spack_utils.o rconstants.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModControlVars.o : $(INIT)/ModControlVars.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCoriolis.o : $(MODEL)/ModCoriolis.f90 node_mod.o mem_scratch.o mem_grid.o \
	ModBuffering.o ModBasicFields.o rconstants.o parlibf.o mem_tend.o \
	ref_sounding.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuParGrell3.o : $(CUPARM)/ModCuParGrell3.F90 ModNamelistFile.o mem_grell.o \
	Phys_const.o ModMessageSet.o ModGrid.o mem_cuparm.o ModRadvc.o module_cu_gf.o \
	mem_radiate.o ModBasicFields.o ModVarTables.o node_mod.o mem_jules.o \
	mem_carma.o ModMicroFields.o module_cu_g3.o ModChemConvTransp.o \
	module_cu_gf_v5.1.o mem_leaf.o mem_tend.o mem_scratch1_grell.o ModRstilt.o \
	grid_dims.o mem_scratch.o mem_chem1.o mem_grell_param2.o ccatt_start.o \
	io_params.o ModMicControl.o mem_stilt.o mem_grid.o ConvPar_GF_GEOS5.o \
	rconstants.o mem_varinit.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDateUtils.o : $(UTILS_MODS)/ModDateUtils.f90 dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffSclr.o : $(TURB)/ModDiffSclr.f90 ModTurbDiff.o mem_scratch.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffuse.o : $(TURB)/ModDiffuse.f90 ModNamelistFile.o mem_scratch.o node_mod.o \
	ModMicControl.o ModMicroFields.o mem_grid.o ModTurbFields.o ModTurbKE.o \
	ModTurbK.o ModTurbDiff.o ModBasicFields.o ModScalarTable.o mem_leaf.o ke_coms.o \
	ModDiffSclr.o mem_opt_scratch.o mem_tend.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDomainDecomp.o : $(MPI)/ModDomainDecomp.f90 ModParallelEnvironment.o \
	ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModEvaluation.o : $(EVAL)/ModEvaluation.f90 ModNamelistFile.o node_mod.o \
	mem_grid.o parlibf.o 
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

ModGasPart.o : $(TEB_SPM)/ModGasPart.f90 mem_emiss.o grid_dims.o node_mod.o \
	an_header.o ModRcio.o mem_grid.o mem_gaspart.o mem_leaf.o mem_teb_vars_const.o \
	parlibf.o ModBasicFields.o ModVarTables.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGeodat.o : $(MKSFC)/ModGeodat.f90 teb_spm_start.o mem_grid.o mem_leaf.o \
	io_params.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGrid.o : $(MPI)/ModGrid.f90 ModNamelistFile.o ModControlVars.o \
	ModMicControl.o ModMessageSet.o ModMicroFields.o ModTurbFields.o \
	ModNodeDimensions.o ModParallelEnvironment.o ModBasicFields.o meteogramType.o \
	ModScalarTable.o ModGridDims.o mem_tend.o ModVarTables.o ModDomainDecomp.o \
	ModNeighbourNodes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridDims.o : $(MPI)/ModGridDims.f90 ModNamelistFile.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridTree.o : $(MPI)/ModGridTree.f90 ModNamelistFile.o \
	ModParallelEnvironment.o ModGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modIau.o : $(MODEL)/modIau.f90 ModNamelistFile.o node_mod.o ReadBcst.o \
	ModMPassFull.o mem_grid.o dump.o parlibf.o mem_tend.o mem_varinit.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitHis.o : $(IO)/ModInitHis.f90 ModLeafComs.o ModRinit.o mem_aerad.o \
	ModRamsReadHeader.o ModBasicFields.o ModVarTables.o ref_sounding.o ModRcio.o \
	mem_leaf.o ModLeaf3.o an_header.o mem_scratch.o mem_chem1.o io_params.o \
	ModGetVar.o ModMicControl.o ModRamsGrid.o mem_grid.o rconstants.o chem1_list.o \
	mem_varinit.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitMicThompson.o : $(MICRO)/ModInitMicThompson.f90 node_mod.o generic.o \
	ReadBcst.o ModMicroFields.o mem_grid.o dump.o ModDateUtils.o parlibf.o \
	ModBasicFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLanduseInput.o : $(MKSFC)/ModLanduseInput.f90 mem_mksfc.o ModLeaf3Init.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3.o : $(SURFACE)/ModLeaf3.f90 ModLeafComs.o teb_spm_start.o mem_scratch.o \
	node_mod.o ModMicControl.o ModMicroFields.o ccatt_start.o mem_grid.o \
	ModTurbFields.o io_params.o mem_teb_common.o ModLeaf3Hyd.o mem_cuparm.o \
	rconstants.o mem_leaf.o mem_radiate.o ModLeaf3Teb.o ModBasicFields.o mem_teb.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3OceanOnly.o : $(SURFACE)/ModLeaf3OceanOnly.f90 ModLeafComs.o \
	ModCuParGrell3.o node_mod.o ModLeaf3.o ccatt_start.o mem_grid.o ModTurbFields.o \
	io_params.o ConvPar_GF_GEOS5.o mem_cuparm.o rconstants.o mem_leaf.o \
	mem_radiate.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMatrixDriver.o : $(MATRIX)/ModMatrixDriver.F90 node_mod.o ModParticle.o \
	ModMicControl.o mem_chem1.o setup.o npf.o mem_aer1.o ModMicroFields.o \
	memMatrix.o mem_grid.o ModTurbFields.o isrpia.o coag.o subs.o rconstants.o \
	mem_leaf.o mem_radiate.o chem1_list.o ModBasicFields.o aer1_list.o 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ModMemAlloc.o : $(MEMORY)/ModMemAlloc.F90 ModLeafComs.o mem_grell.o ModGrid.o \
	mem_scratch3_grell_sh.o mem_cuparm.o mem_radiate.o mem_aerad.o ModBasicFields.o \
	mem_chemic.o ModVarTables.o mem_scratch3_grell.o ModCuParGrell3.o node_mod.o \
	machine_arq.o mem_jules.o mem_carma.o mem_plume_chem1.o carma_fastjx.o \
	mem_aer1.o mem_grid_dim_defs.o mem_chem1aq.o mem_teb_common.o shcu_vars_const.o \
	ModParallelEnvironment.o chem_sources.o mem_tuv.o mem_leaf.o ModEvaluation.o \
	mem_teb_vars_const.o mem_oda.o mem_teb.o mem_tend.o mem_scratch1_brams.o \
	mem_scratch1_grell.o extra.o mem_turb_scalar.o aer1_list.o mem_globrad.o \
	grid_dims.o mem_shcu.o mem_scratch.o mem_chem1.o mem_grell_param2.o \
	ccatt_start.o ModTurbFields.o mem_gaspart.o io_params.o mem_nestb.o \
	mem_globaer.o mem_opt_scratch.o ModTuv2.7.o mem_scratch2_grell.o \
	mem_scratch2_grell_sh.o mem_micro.o teb_spm_start.o mem_emiss.o chem1aq_list.o \
	digitalFilter.o modIau.o parrrsw.o mem_stilt.o ModOptical.o mem_grid.o \
	mem_scalar.o chem_dry_dep.o chem1_list.o mem_volc_chem1.o mem_varinit.o 
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

ModMessageSet.o : $(MPI)/ModMessageSet.f90 ModNamelistFile.o ModMessageData.o \
	ModDomainDecomp.o mem_grid.o ModNodeDimensions.o ModFieldSection.o \
	ModParallelEnvironment.o ModGridDims.o parlibf.o ModVarTables.o \
	ModFieldSectionList.o ModNeighbourNodes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicGfdlDriver.o : $(MICRO)/ModMicGfdlDriver.f90 gfdl_cloud_microphys.o \
	node_mod.o ModMicroFields.o mem_grid.o io_params.o rconstants.o mem_leaf.o \
	mem_radiate.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsDrive.o : $(MICRO)/ModMicrophysicsDrive.f90 grid_dims.o \
	node_mod.o ModMicControl.o mem_chem1.o ModMicroFields.o mem_grid.o \
	ModMicrophysicsMisc.o mem_chem1aq.o ModMicColl.o ModMicNuc.o ModMicTabs.o \
	ModMicVap.o mem_radiate.o ModBasicFields.o mem_chemic.o ModMicInit.o \
	$(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsMisc.o : $(MICRO)/ModMicrophysicsMisc.f90 mem_scratch.o \
	ModMicControl.o ModMicroFields.o mem_grid.o rconstants.o ModBasicFields.o \
	$(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicThompsonDriver.o : $(MICRO)/ModMicThompsonDriver.f90 node_mod.o \
	module_mp_thompson.o ModMicControl.o ModMicroFields.o mem_grid.o io_params.o \
	rconstants.o mem_leaf.o mem_radiate.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcDriver.o : $(MKSFC)/ModMkSfcDriver.f90 teb_spm_start.o grid_dims.o \
	node_mod.o ModMkSfcNdvi.o ModControlVars.o ModNestGeoSst.o mem_mksfc.o \
	ReadBcst.o mem_grid.o ModNdviRead.o io_params.o ModMkSfcSfc.o ModMkSfcFuso.o \
	ModLanduseInput.o ModMkSfcSst.o ModMkSfcTop.o ModSstRead.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcFuso.o : $(MKSFC)/ModMkSfcFuso.f90 mem_emiss.o node_mod.o \
	ModControlVars.o ReadBcst.o mem_mksfc.o mem_grid.o mem_gaspart.o io_params.o \
	mem_teb_vars_const.o mem_teb.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcNdvi.o : $(MKSFC)/ModMkSfcNdvi.f90 ModRUser.o mem_mksfc.o mem_grid.o \
	io_params.o dump.o mem_leaf.o ModLanduseInput.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSfc.o : $(MKSFC)/ModMkSfcSfc.f90 node_mod.o ModControlVars.o ReadBcst.o \
	mem_mksfc.o mem_grid.o io_params.o dump.o mem_leaf.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSst.o : $(MKSFC)/ModMkSfcSst.f90 grid_dims.o ModRUser.o mem_mksfc.o \
	mem_grid.o io_params.o ModGeodat.o mem_leaf.o ModNestFillDens.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcTop.o : $(MKSFC)/ModMkSfcTop.f90 node_mod.o ModControlVars.o ReadBcst.o \
	mem_mksfc.o mem_grid.o io_params.o dump.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMonotonicAdvection.o : $(MODEL)/ModMonotonicAdvection.f90 ModNamelistFile.o \
	ModMicControl.o mem_chem1.o ModMessageSet.o ModGrid.o mem_aer1.o mem_grid.o \
	ccatt_start.o ModParallelEnvironment.o rconstants.o chem_dry_dep.o \
	ModDomainDecomp.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNamelistFile.o : $(INIT)/ModNamelistFile.f90 grid_dims.o modPrintInitial.o \
	ModParallelEnvironment.o dump.o parlibf.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNdviRead.o : $(MKSFC)/ModNdviRead.f90 grid_dims.o node_mod.o ModMkSfcNdvi.o \
	ModControlVars.o ReadBcst.o mem_mksfc.o mem_grid.o io_params.o ModDateUtils.o \
	mem_leaf.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNeighbourNodes.o : $(MPI)/ModNeighbourNodes.f90 ModParallelEnvironment.o \
	ModDomainDecomp.o ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestGeoSst.o : $(MKSFC)/ModNestGeoSst.f90 dump.o ModBasicFields.o node_mod.o \
	memSoilMoisture.o mem_leaf.o ModLanduseInput.o ModMkSfcTop.o grid_dims.o \
	mem_scratch.o ccatt_start.o ModTurbFields.o io_params.o ModGeodat.o \
	ModNestFillDens.o ModNestFeed.o ModLeaf3Init.o ModRUser.o ModControlVars.o \
	mem_mksfc.o ModInitHis.o mem_grid.o ModSoilMoisture.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNodeDimensions.o : $(MPI)/ModNodeDimensions.f90 ModParallelEnvironment.o \
	ModDomainDecomp.o ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudAnalysis.o : $(FDDA)/ModNudAnalysis.f90 mem_scratch.o node_mod.o \
	mem_chem1.o mem_grid.o mem_tend.o dump.o ModEvaluation.o chem1_list.o \
	ModBasicFields.o ModNestFillDens.o mem_varinit.o modIau.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaNudge.o : $(FDDA)/ModOdaNudge.f90 node_mod.o mem_scratch.o ModOdaKrig.o \
	mem_grid.o mem_tend.o io_params.o ModOdaProcObs.o mem_oda.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOneProc.o : $(MODEL)/ModOneProc.F90 ModCuRead.o ModLeaf3Teb.o mem_carma.o \
	ModSched.o ModParallelEnvironment.o mem_teb.o mem_globrad.o grid_dims.o \
	mem_chem1.o domain_decomp.o ccatt_start.o ModRecycle.o ModNudRead.o mem_micro.o \
	teb_spm_start.o ModRio.o local_proc.o chem_dry_dep.o chem1_list.o ModMkSfcSfc.o \
	mem_varinit.o ModAerClim.o ModOpspec.o ModMicrophysicsMisc.o mem_cuparm.o \
	ModMkSfcFuso.o ModBasicFields.o ref_sounding.o mem_chem1aq.o \
	ModPostGridNetCDF.o mem_teb_vars_const.o extra.o ModParaInit.o ModNestIntrp.o \
	ModOdaRead.o mem_gaspart.o io_params.o ModVarfFile.o ModGridTree.o \
	ModMPassDtl.o ModWindFarm.o ModTuv2.7.o ModLeaf3Init.o mem_emiss.o \
	tuvParameter.o ModInitHis.o ModRamsGrid.o ModSoilMoisture.o mem_scalar.o \
	ModMonotonicAdvection.o ModChemAsgen.o ModCuParGrell3.o mem_plume_chem1.o \
	mem_aer1.o memSoilMoisture.o ModTimestep.o ModRThrm.o ModTimestepRK.o \
	shcu_vars_const.o meteogram.o aer1_list.o ModMkSfcTop.o ModRnode.o \
	ModCoriolis.o ModVarfUpdate.o ModChemistryDriver.o modIau.o \
	ModRamsMicrophysics2M.o ModNestGeoSst.o mem_stilt.o ModInitMicThompson.o \
	mem_grid.o dam.o isan_coms.o mem_volc_chem1.o ModMkSfcDriver.o \
	ModDomainDecomp.o ModSstRead.o ModNamelistFile.o ModMicInit.o ModCondRead.o \
	ModRinit.o ModGrid.o dump.o mem_radiate.o parlibf.o ModVarTables.o ModGasPart.o \
	node_mod.o machine_arq.o ReadBcst.o ModTuvDriver2.7.o ModNdviRead.o \
	mem_teb_common.o chem_sources.o mem_leaf.o ModEvaluation.o mem_oda.o \
	mem_scratch.o ModMemAlloc.o mem_grell_param2.o ModTimeStamp.o ModUrbanCanopy.o \
	ModRanlavg.o digitalFilter.o ModRUser.o ModRhhi.o ModPostProcess.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOpspec.o : $(IO)/ModOpspec.f90 ModNamelistFile.o mem_cuparm.o mem_radiate.o \
	mem_aer1.o mem_chem1aq.o shcu_vars_const.o chem_sources.o mem_leaf.o \
	aer1_list.o mem_globrad.o grid_dims.o mem_chem1.o mem_grell_param2.o \
	ccatt_start.o io_params.o modIau.o teb_spm_start.o mem_emiss.o chem1aq_list.o \
	ModMicControl.o mem_stilt.o mem_grid.o chem1_list.o mem_varinit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOptical.o : $(RADIATE)/ModOptical.f90 ModNamelistFile.o mem_radiate.o \
	node_mod.o ModControlVars.o ReadBcst.o ModRamsGrid.o ModMPassFull.o \
	ModTurbFields.o ModSoilMoisture.o ccatt_start.o mem_aer1.o mem_grid.o dump.o \
	mem_leaf.o parlibf.o ModBasicFields.o aer1_list.o ModVarTables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOutputUtils.o : $(IO)/ModOutputUtils.f90 ModNamelistFile.o ModTurbFields.o \
	dump.o ModBasicFields.o ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOzone.o : $(TEB_SPM)/ModOzone.f90 mem_grid.o mem_gaspart.o rconstants.o \
	mem_radiate.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParaInit.o : $(MPI)/ModParaInit.f90 grid_dims.o node_mod.o mem_grid.o dump.o \
	ModScalarTable.o ModVarTables.o $(UTILS_INCS)/constants.h 
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

ModPostGrid.o : $(POST_SRC)/ModPostGrid.F90 ModNamelistFile.o ModPostUtils.o \
	ModPostTypes.o ModPostOneFieldNetCDF.o mem_grid.o ModTurbFields.o io_params.o \
	ModParallelEnvironment.o ModOutputUtils.o parlibf.o ModBasicFields.o \
	ModBramsGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostGridNetCDF.o : $(POST_SRC)/ModPostGridNetCDF.F90 ModNamelistFile.o \
	ModPostUtils.o ModPostTypes.o mem_grid.o io_params.o dump.o ModDateUtils.o \
	ModBramsGrid.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField.o : $(POST_SRC)/ModPostOneField.f90 ModNamelistFile.o \
	ModPostUtils.o node_mod.o ModMicControl.o ModPostTypes.o ModTurbFields.o \
	ModPostOneField2d.o ModPostOneFieldUtils.o ModBasicFields.o dump.o \
	ModPostOneField3d.o ModPostOneField8d.o ModBramsGrid.o ModPostOneField7d.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField2d.o : $(POST_SRC)/ModPostOneField2d.f90 ModNamelistFile.o \
	ModPostUtils.o node_mod.o ModMicControl.o ModPostTypes.o mem_grid.o \
	ModTurbFields.o io_params.o mem_cuparm.o ModOutputUtils.o \
	ModPostOneFieldUtils.o dump.o mem_radiate.o mem_aerad.o ModBasicFields.o \
	ModBramsGrid.o ModPostGrid.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField3d.o : $(POST_SRC)/ModPostOneField3d.f90 ModNamelistFile.o \
	ModPostUtils.o node_mod.o ModMicControl.o ModPostTypes.o mem_grid.o \
	ModTurbFields.o ModOutputUtils.o ModPostOneFieldUtils.o ModBasicFields.o \
	ModBramsGrid.o ModPostGrid.o mem_varinit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField7d.o : $(POST_SRC)/ModPostOneField7d.f90 ModNamelistFile.o \
	ModPostUtils.o ModPostTypes.o ModTurbFields.o ModOutputUtils.o \
	ModPostOneFieldUtils.o ModBasicFields.o ModBramsGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField8d.o : $(POST_SRC)/ModPostOneField8d.f90 ModNamelistFile.o \
	ModPostUtils.o ModPostTypes.o ModTurbFields.o ModOutputUtils.o \
	ModPostOneFieldUtils.o ModBasicFields.o ModBramsGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneFieldNetCDF.o : $(POST_SRC)/ModPostOneFieldNetCDF.F90 \
	ModPostGridNetCDF.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneFieldUtils.o : $(POST_SRC)/ModPostOneFieldUtils.f90 ModBramsGrid.o \
	ModPostTypes.o ModPostGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostProcess.o : $(POST_SRC)/ModPostProcess.F90 ModNamelistFile.o \
	ModPostTypes.o ModMessageSet.o ModGrid.o ModTurbFields.o io_params.o \
	ModPostOneField.o ModParallelEnvironment.o ModGridTree.o ModPostGridNetCDF.o \
	ModBasicFields.o ModBramsGrid.o ModPostGrid.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostTypes.o : $(POST_SRC)/ModPostTypes.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostUtils.o : $(POST_SRC)/ModPostUtils.f90 ModParallelEnvironment.o dump.o \
	mem_leaf.o $(POST_INCS)/post_rconfig.h $(POST_INCS)/post_rconstants.h \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modPrintInitial.o : $(INIT)/modPrintInitial.F90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvc.o : $(MODEL)/ModRadvc.f90 ModNamelistFile.o ModMonotonicAdvection.o \
	mem_scratch.o grid_dims.o mem_chem1.o ccatt_start.o mem_grid.o mem_aer1.o \
	ModParallelEnvironment.o ModBasicFields.o ModScalarTable.o chem_dry_dep.o \
	mem_tend.o ModRadvcAdap.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcRK.o : $(MODEL)/ModRadvcRK.f90 grid_dims.o node_mod.o mem_chem1.o \
	mem_stilt.o ModRexev.o ModGrid.o ModMessageSet.o mem_grid.o \
	ModParallelEnvironment.o mem_tend.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsMicrophysics2M.o : $(MICRO)/ModRamsMicrophysics2M.f90 grid_dims.o \
	node_mod.o mem_scratch.o ModMicControl.o ModMicroFields.o mem_grid.o \
	rconstants.o mem_leaf.o dump.o ModBasicFields.o ModMicGamma.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsReadHeader.o : $(IO)/ModRamsReadHeader.f90 an_header.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRanlavg.o : $(IO)/ModRanlavg.f90 grid_dims.o ModMicControl.o mem_grid.o \
	io_params.o ModRThrm.o ModBasicFields.o ModVarTables.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRbnd.o : $(BC)/ModRbnd.f90 node_mod.o mem_scratch.o ModMicControl.o \
	mem_chem1.o ModMicroFields.o ccatt_start.o mem_grid.o ModTurbFields.o \
	ModMicrophysicsMisc.o ModTurbKE.o mem_tend.o ModScalarTable.o ModBasicFields.o \
	ref_sounding.o mem_micro.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRcio.o : $(IO)/ModRcio.f90 ModLeafComs.o grid_dims.o an_header.o \
	ModNamelistFile.o ModMicControl.o mem_stilt.o mem_grid.o io_params.o \
	mem_cuparm.o ModParallelEnvironment.o mem_leaf.o mem_radiate.o ref_sounding.o \
	$(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConv.o : $(CUPARM)/ModRConv.f90 node_mod.o mem_scratch.o mem_grid.o \
	mem_tend.o mem_cuparm.o ModBasicFields.o rconstants.o ModConvComs.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConvGrellCatt.o : $(CUPARM)/ModRConvGrellCatt.f90 mem_grell.o ModGrid.o \
	mem_cuparm.o mem_radiate.o ModCupGrellCattShallow.o ModCuParGrell3.o node_mod.o \
	ModChemConvTransp.o mem_leaf.o mem_tend.o mem_scratch1_grell.o \
	ModCupGrellCattDeep.o ModRstilt.o mem_scratch.o mem_grell_param2.o \
	ccatt_start.o io_params.o ModMicControl.o mem_stilt.o mem_grid.o mem_scalar.o \
	rconstants.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRecycle.o : $(IO)/ModRecycle.f90 node_mod.o an_header.o mem_chem1.o \
	ReadBcst.o mem_aer1.o ModMPassFull.o mem_grid.o io_params.o ModGetVar.o \
	ModDateUtils.o dump.o mem_aerad.o chem1_list.o ModRamsReadHeader.o \
	ModVarTables.o aer1_list.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRexev.o : $(STILT)/ModRexev.f90 mem_scratch.o ModMicControl.o mem_stilt.o \
	mem_grid.o ModBasicFields.o rconstants.o ModRadvc.o mem_tend.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRGrad.o : $(TURB)/ModRGrad.f90 mem_scratch.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRhhi.o : $(INIT)/ModRhhi.f90 grid_dims.o ModRinit.o mem_scratch.o \
	ModMicControl.o ModRamsGrid.o mem_grid.o rconstants.o ModBasicFields.o \
	ref_sounding.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRinit.o : $(INIT)/ModRinit.f90 node_mod.o mem_scratch.o ModMicControl.o \
	ModMicroFields.o ModRbnd.o mem_grid.o ModTurbFields.o ModTurbKE.o io_params.o \
	rconstants.o ModBasicFields.o ref_sounding.o mem_varinit.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRio.o : $(IO)/ModRio.f90 ModNamelistFile.o ModMPassFull.o mpi_io_engine-5d.o \
	parlibf.o utilsMod.o mem_aerad.o ModBasicFields.o ModVarTables.o ref_sounding.o \
	node_mod.o ModRcio.o ReadBcst.o ModParallelEnvironment.o ModDateUtils.o \
	grid_dims.o an_header.o mem_chem1.o ModTurbFields.o io_params.o \
	ModControlVars.o ModMicControl.o mem_grid.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/interface.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRrtmDriver.o : $(RADIATE)/ModRrtmDriver.f90 ModLeafComs.o \
	mcica_subcol_gen_lw.o mcica_subcol_gen_sw.o mem_cuparm.o mem_radiate.o \
	ModBasicFields.o ref_sounding.o node_mod.o mem_carma.o ModMicroFields.o \
	mem_tuv.o ModDateUtils.o mem_leaf.o rrtmg_sw_rad.o mem_tend.o \
	mem_scratch1_grell.o grid_dims.o rrtmg_lw_cldprop.o mem_chem1.o parrrtm.o \
	mem_grell_param2.o ccatt_start.o parkind.o mem_rrtm.o teb_spm_start.o parrrsw.o \
	ModMicControl.o ModOptical.o mem_grid.o rconstants.o rrtmg_lw_rad.o \
	rrtmg_sw_cldprop.o $(UTILS_INCS)/aerosol_setup.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRShCuPar.o : $(CUPARM)/ModRShCuPar.f90 mem_shcu.o node_mod.o mem_scratch.o \
	ModMicroFields.o mem_grid.o ModTurbFields.o mem_tend.o shcu_vars_const.o \
	ModBasicFields.o ModConvComs.o ModRConv.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRstilt.o : $(STILT)/ModRstilt.f90 grid_dims.o mem_scratch.o \
	ModMonotonicAdvection.o mem_stilt.o mem_grid.o ModTurbFields.o mem_cuparm.o \
	ModBasicFields.o mem_scratch1_grell.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRThrm.o : $(MODEL)/ModRThrm.f90 mem_scratch.o ModMicControl.o mem_grid.o \
	rconstants.o ModBasicFields.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRtimi.o : $(MODEL)/ModRtimi.f90 mem_scratch.o node_mod.o mem_grell.o \
	mem_grid.o shcu_vars_const.o mem_cuparm.o ModBasicFields.o ModScalarTable.o \
	mem_tend.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModScalarTable.o : $(MEMORY)/ModScalarTable.f90 ModParallelEnvironment.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSched.o : $(MODEL)/ModSched.f90 ModNamelistFile.o node_mod.o ReadBcst.o \
	local_proc.o mem_grid.o io_params.o shcu_vars_const.o mem_cuparm.o parlibf.o \
	dump.o mem_radiate.o ModBasicFields.o isan_coms.o ref_sounding.o mem_varinit.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSeaSalt.o : $(CCATT)/ModSeaSalt.f90 ModAerClim.o node_mod.o mem_chem1.o \
	mem_aer1.o ccatt_start.o mem_grid.o io_params.o mem_leaf.o ModBasicFields.o \
	aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSstRead.o : $(MKSFC)/ModSstRead.f90 grid_dims.o node_mod.o ModControlVars.o \
	ReadBcst.o mem_mksfc.o mem_grid.o io_params.o ModDateUtils.o mem_leaf.o \
	ModMkSfcSst.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimeStamp.o : $(MODEL)/ModTimeStamp.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestep.o : $(MODEL)/ModTimestep.F90 ModMonotonicAdvection.o ModRtimi.o \
	ModMessageSet.o ModGrid.o ModRexev.o ModMicrophysicsMisc.o mem_cuparm.o \
	mem_radiate.o ModRadvc.o sfclyr_jules.o ModDiffuse.o ModBasicFields.o \
	ModMicThompsonDriver.o ModGasPart.o ModCuParGrell3.o node_mod.o machine_arq.o \
	ModMatrixDriver.o ModOdaNudge.o mem_plume_chem1.o mem_aer1.o ModSeaSalt.o \
	shcu_vars_const.o ModRThrm.o ModAcoust.o chem_sources.o mem_leaf.o \
	ModNudAnalysis.o mem_oda.o mem_tend.o ModCoriolis.o ModRConvGrellCatt.o \
	ModMicGfdlDriver.o ModLeaf3.o ModRstilt.o grid_dims.o mem_scratch.o mem_chem1.o \
	ChemDryDepDriver.o ModRShCuPar.o ModRbnd.o ModMicrophysicsDrive.o ccatt_start.o \
	ModChemistryDriver.o ModTimeStamp.o ModUrbanCanopy.o ModTurbK.o rad_driv.o \
	ModWindFarm.o ModOzone.o ModRConv.o digitalFilter.o mem_micro.o teb_spm_start.o \
	mem_emiss.o ModRamsMicrophysics2M.o mem_stilt.o ModOptical.o \
	ChemSourcesDriver.o mem_grid.o mem_scalar.o rconstants.o mem_varinit.o \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestepRK.o : $(MODEL)/ModTimestepRK.F90 ModAerClim.o \
	ModMonotonicAdvection.o ModRtimi.o ModMessageSet.o ModRexev.o ModGrid.o \
	ModMicrophysicsMisc.o mem_cuparm.o mem_radiate.o utilsMod.o ModRadvc.o \
	ModDiffuse.o sfclyr_jules.o ModMicThompsonDriver.o ModGasPart.o \
	ModCuParGrell3.o node_mod.o machine_arq.o ModMatrixDriver.o ModOdaNudge.o \
	mem_plume_chem1.o ModRadvcRK.o mem_aer1.o ModSeaSalt.o ModTimestep.o ModRThrm.o \
	ModParallelEnvironment.o ModAcoust.o shcu_vars_const.o chem_sources.o \
	mem_leaf.o ModNudAnalysis.o mem_oda.o mem_tend.o ModCoriolis.o \
	ModMicGfdlDriver.o ModRstilt.o ModLeaf3.o grid_dims.o ModLeaf3OceanOnly.o \
	mem_scratch.o mem_chem1.o ChemDryDepDriver.o ModRShCuPar.o ModRbnd.o \
	ModMicrophysicsDrive.o ccatt_start.o ModChemistryDriver.o ModTimeStamp.o \
	ModUrbanCanopy.o ModTurbK.o rad_driv.o ModWindFarm.o ModOzone.o ModRConv.o \
	digitalFilter.o mem_micro.o teb_spm_start.o mem_emiss.o modIau.o \
	ModRamsMicrophysics2M.o mem_stilt.o ModOptical.o ChemSourcesDriver.o mem_grid.o \
	mem_scalar.o rconstants.o mem_varinit.o $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiff.o : $(TURB)/ModTurbDiff.f90 mem_scratch.o mem_grid.o mem_cuparm.o \
	mem_opt_scratch.o ModRGrad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiffAdap.o : $(TURB)/ModTurbDiffAdap.f90 mem_scratch.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbFields.o : $(TURB)/ModTurbFields.f90 ModNamelistFile.o \
	ModNodeDimensions.o ModParallelEnvironment.o ModVarTables.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbK.o : $(TURB)/ModTurbK.f90 ModNamelistFile.o ModMonotonicAdvection.o \
	mem_grell.o ModTurbKAdap.o ModTurbKE.o mem_cuparm.o ke_coms.o ModBasicFields.o \
	node_mod.o ModMicroFields.o mem_leaf.o mem_tend.o mem_turb_scalar.o ModRGrad.o \
	ModRstilt.o grid_dims.o mem_scratch.o mem_chem1.o ccatt_start.o ModTurbFields.o \
	ModTurbDiffAdap.o ModMicControl.o mem_stilt.o mem_grid.o ModTurbDiff.o \
	rconstants.o ModScalarTable.o ModTKenn.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKAdap.o : $(TURB)/ModTurbKAdap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKE.o : $(TURB)/ModTurbKE.f90 mem_scratch.o mem_grid.o ModTurbFields.o \
	rconstants.o ke_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuv2.7.o : $(TUV)/ModTuv2.7.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuvDriver2.7.o : $(TUV)/ModTuvDriver2.7.f90 tuvParameter.o mem_rrtm.o \
	node_mod.o mem_carma.o mem_chem1.o ModTuv2.7.o mem_grid.o chem1_list.o \
	chem_fastjx_driv.o mem_tuv.o rconstants.o mem_leaf.o mem_radiate.o mem_aerad.o \
	ModBasicFields.o extra.o ref_sounding.o mem_globrad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_cu_g3.o : $(CUPARM)/module_cu_g3.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_cu_gd_fim.o : $(CUPARM)/module_cu_gd_fim.f90 Phys_const.o module_gate.o 
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

ModUrbanCanopy.o : $(SURFACE)/ModUrbanCanopy.f90 node_mod.o mem_grid.o \
	ModTurbFields.o ModBasicFields.o mem_tend.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfFile.o : $(FDDA)/ModVarfFile.f90 ModMessageSet.o ModGrid.o parlibf.o \
	ModBasicFields.o ModRamsReadHeader.o ref_sounding.o node_mod.o ModRcio.o \
	ReadBcst.o mem_aer1.o ModDateUtils.o mem_leaf.o ModNudAnalysis.o aer1_list.o \
	mem_scratch.o mem_chem1.o ModVarfUpdate.o ModGetVar.o ModGridTree.o \
	ModControlVars.o ModMicControl.o ModRamsGrid.o mem_grid.o rconstants.o \
	chem1_list.o isan_coms.o mem_varinit.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModWindFarm.o : $(WIND_FARM)/ModWindFarm.f90 ModNamelistFile.o node_mod.o \
	mem_grid.o ModTurbFields.o io_params.o ModBasicFields.o rconstants.o \
	ModDateUtils.o module_wind_fitch.o mem_tend.o $(UTILS_INCS)/files.h 
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

ModNestIntrp.o : $(NESTING)/ModNestIntrp.f90 grid_dims.o ModRinit.o \
	mem_scratch.o mem_grid.o mem_nestb.o rconstants.o ModBasicFields.o \
	ModNestFillDens.o ref_sounding.o 
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

ModNudRead.o : $(FDDA)/ModNudRead.f90 mem_chem1.o ModRamsGrid.o mem_grid.o \
	ModDateUtils.o ModNudAnalysis.o ModNudUpdate.o isan_coms.o mem_varinit.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudUpdate.o : $(FDDA)/ModNudUpdate.f90 an_header.o mem_chem1.o ModRcio.o \
	ModInitHis.o mem_grid.o grid_struct.o mem_aerad.o chem1_list.o ModVarTables.o \
	mem_varinit.o $(UTILS_INCS)/files.h 
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

ModOdaRead.o : $(FDDA)/ModOdaRead.f90 ModOdaStaCount.o mem_grid.o ModDateUtils.o \
	mem_oda.o ModOdaStaInput.o isan_coms.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaCount.o : $(FDDA)/ModOdaStaCount.f90 ModReadRalph.o mem_grid.o \
	mem_oda.o obs_input.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaInput.o : $(FDDA)/ModOdaStaInput.f90 ModReadRalph.o ModOdaStaCount.o \
	mem_grid.o obs_input.o ModDateUtils.o mem_oda.o 
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

ModAcoustAdap.o : $(MODEL)/ModAcoustAdap.f90 mem_scratch.o node_mod.o ModGrid.o \
	ModRbnd.o ModMessageSet.o mem_grid.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_carma.o : $(RADIATE)/rad_carma.F90 grid_dims.o node_mod.o machine_arq.o \
	mem_carma.o carma_fastjx.o mem_chem1.o mem_aer1.o ccatt_start.o mem_grid.o \
	mem_aerad.o mem_tuv.o rconstants.o mem_leaf.o ModDateUtils.o mem_globaer.o \
	mem_radiate.o chem1_list.o aer1_list.o mem_globrad.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_driv.o : $(RADIATE)/rad_driv.f90 ModMicControl.o ModMicroFields.o \
	ModCarmaDriver.o mem_radiate.o ModBasicFields.o ModRrtmDriver.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcAdap.o : $(MODEL)/ModRadvcAdap.f90 ModAdapInit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsGrid.o : $(INIT)/ModRamsGrid.f90 node_mod.o ModAdapInit.o mem_grid.o \
	dump.o ModGridSet.o rconstants.o $(UTILS_INCS)/constants.h 
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

ReadBcst.o : $(MPI)/ReadBcst.f90 node_mod.o ModControlVars.o an_header.o \
	ModMPassFull.o mem_grid.o ModTurbFields.o utilsMod.o parlibf.o mem_aerad.o \
	ModBasicFields.o ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ref_sounding.o : $(MODEL)/ref_sounding.f90 ModNamelistFile.o grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRnode.o : $(MODEL)/ModRnode.f90 grid_dims.o node_mod.o mem_grid.o mem_leaf.o \
	parlibf.o ModVarTables.o $(UTILS_INCS)/constants.h 
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

rrsw_aer.o : $(RRTMG_SW_MOD)/rrsw_aer.f90 parkind.o parrrsw.o 
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

rrsw_kg16.o : $(RRTMG_SW_MOD)/rrsw_kg16.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg17.o : $(RRTMG_SW_MOD)/rrsw_kg17.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg18.o : $(RRTMG_SW_MOD)/rrsw_kg18.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg19.o : $(RRTMG_SW_MOD)/rrsw_kg19.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg20.o : $(RRTMG_SW_MOD)/rrsw_kg20.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg21.o : $(RRTMG_SW_MOD)/rrsw_kg21.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg22.o : $(RRTMG_SW_MOD)/rrsw_kg22.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg23.o : $(RRTMG_SW_MOD)/rrsw_kg23.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg24.o : $(RRTMG_SW_MOD)/rrsw_kg24.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg25.o : $(RRTMG_SW_MOD)/rrsw_kg25.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg26.o : $(RRTMG_SW_MOD)/rrsw_kg26.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg27.o : $(RRTMG_SW_MOD)/rrsw_kg27.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg28.o : $(RRTMG_SW_MOD)/rrsw_kg28.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg29.o : $(RRTMG_SW_MOD)/rrsw_kg29.f90 parkind.o parrrsw.o 
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

rrsw_wvn.o : $(RRTMG_SW_MOD)/rrsw_wvn.f90 parkind.o parrrsw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprmc.f90 rrlw_cld.o parrrtm.o \
	parkind.o rrlw_wvn.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprop.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprop.f90 parkind.o rrlw_cld.o \
	parrrtm.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_init.o : $(RRTMG_LW_SRC)/rrtmg_lw_init.f90 rrlw_kg16.o rrlw_kg08.o \
	rrlw_kg05.o rrlw_kg12.o rrlw_kg06.o rrlw_wvn.o rrlw_tbl.o rrlw_cld.o \
	rrlw_kg15.o rrlw_kg13.o rrlw_kg09.o rrlw_kg11.o rrlw_kg07.o rrtmg_lw_setcoef.o \
	rrlw_kg01.o rrlw_kg03.o rrtmg_lw_k_g.o parrrtm.o parkind.o rrlw_kg14.o \
	rrlw_kg10.o rrlw_con.o rrlw_kg02.o rrlw_kg04.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_k_g.o : $(RRTMG_LW_SRC)/rrtmg_lw_k_g.f90 rrlw_kg16.o rrlw_kg01.o \
	rrlw_kg03.o rrlw_kg08.o rrlw_kg05.o parkind.o rrlw_kg15.o rrlw_kg13.o \
	rrlw_kg09.o rrlw_kg11.o rrlw_kg07.o rrlw_kg12.o rrlw_kg06.o rrlw_kg14.o \
	rrlw_kg10.o rrlw_kg02.o rrlw_kg04.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rad.o : $(RRTMG_LW_SRC)/rrtmg_lw_rad.f90 parrrtm.o \
	mcica_subcol_gen_lw.o rrtmg_lw_taumol.o parkind.o rrtmg_lw_rtrnmc.o \
	rrtmg_lw_cldprmc.o rrlw_con.o rrtmg_lw_setcoef.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrn.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrn.f90 parrrtm.o parkind.o \
	rrlw_wvn.o rrlw_con.o rrlw_tbl.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmc.f90 parrrtm.o parkind.o \
	rrlw_wvn.o rrlw_con.o rrlw_tbl.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmr.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmr.f90 parrrtm.o parkind.o \
	rrlw_wvn.o rrlw_con.o rrlw_tbl.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_setcoef.o : $(RRTMG_LW_SRC)/rrtmg_lw_setcoef.f90 rrlw_ref.o parrrtm.o \
	parkind.o rrlw_wvn.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_taumol.o : $(RRTMG_LW_SRC)/rrtmg_lw_taumol.f90 rrlw_kg16.o rrlw_kg08.o \
	rrlw_kg05.o rrlw_kg12.o rrlw_kg06.o rrlw_wvn.o rrlw_kg15.o rrlw_kg13.o \
	rrlw_kg09.o rrlw_kg11.o rrlw_kg07.o rrlw_ref.o rrlw_kg01.o rrlw_kg03.o \
	parrrtm.o parkind.o rrlw_kg14.o rrlw_kg10.o rrlw_con.o rrlw_kg02.o rrlw_kg04.o \
	rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprmc.f90 rrsw_cld.o parrrsw.o \
	parkind.o rrsw_vsn.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprop.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprop.f90 rrsw_cld.o parrrsw.o \
	parkind.o rrsw_vsn.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_init.o : $(RRTMG_SW_SRC)/rrtmg_sw_init.f90 rrtmg_sw_k_g.o rrsw_kg22.o \
	rrsw_kg24.o rrsw_tbl.o rrsw_kg21.o rrsw_wvn.o rrsw_kg28.o rrsw_aer.o rrsw_vsn.o \
	rrsw_con.o rrsw_kg17.o rrtmg_sw_setcoef.o parkind.o rrsw_kg19.o rrsw_kg23.o \
	rrsw_kg20.o rrsw_kg29.o rrsw_kg25.o rrsw_kg18.o rrsw_kg16.o rrsw_cld.o \
	parrrsw.o rrsw_kg27.o rrsw_kg26.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_k_g.o : $(RRTMG_SW_SRC)/rrtmg_sw_k_g.f90 rrsw_kg28.o parkind.o \
	rrsw_kg22.o rrsw_kg19.o rrsw_kg23.o rrsw_kg20.o rrsw_kg24.o rrsw_kg25.o \
	rrsw_kg29.o rrsw_kg27.o rrsw_kg18.o rrsw_vsn.o rrsw_kg26.o rrsw_kg16.o \
	rrsw_kg21.o rrsw_kg17.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_rad.o : $(RRTMG_SW_SRC)/rrtmg_sw_rad.f90 rrsw_aer.o parrrsw.o \
	rrtmg_sw_setcoef.o parkind.o mcica_subcol_gen_sw.o rrtmg_sw_spcvmc.o rrsw_con.o \
	rrtmg_sw_cldprmc.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_reftra.o : $(RRTMG_SW_SRC)/rrtmg_sw_reftra.f90 parkind.o rrsw_vsn.o \
	rrsw_tbl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_setcoef.o : $(RRTMG_SW_SRC)/rrtmg_sw_setcoef.f90 parkind.o rrsw_ref.o \
	parrrsw.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvmc.f90 parrrsw.o parkind.o \
	rrtmg_sw_vrtqdr.o rrtmg_sw_reftra.o rrsw_tbl.o rrsw_vsn.o rrsw_wvn.o \
	rrtmg_sw_taumol.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvrt.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvrt.f90 parrrsw.o parkind.o \
	rrtmg_sw_vrtqdr.o rrtmg_sw_reftra.o rrsw_tbl.o rrsw_vsn.o rrsw_wvn.o \
	rrtmg_sw_taumol.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_taumol.o : $(RRTMG_SW_SRC)/rrtmg_sw_taumol.f90 rrsw_kg28.o parrrsw.o \
	parkind.o rrsw_kg22.o rrsw_kg19.o rrsw_kg23.o rrsw_kg20.o rrsw_kg24.o \
	rrsw_kg25.o rrsw_kg29.o rrsw_kg27.o rrsw_vsn.o rrsw_kg18.o rrsw_kg26.o \
	rrsw_con.o rrsw_kg16.o rrsw_kg21.o rrsw_wvn.o rrsw_kg17.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_vrtqdr.o : $(RRTMG_SW_SRC)/rrtmg_sw_vrtqdr.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRUser.o : $(SURFACE)/ModRUser.f90 ModLeafComs.o node_mod.o ccatt_start.o \
	mem_grid.o memSoilMoisture.o io_params.o rconstants.o mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

setup.o : $(MATRIX)/setup.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

sfclyr_jules.o : $(JULES_DIR)/sfclyr_jules.f90 ModLeafComs.o io_constants.o \
	mem_cuparm.o mem_radiate.o ModBasicFields.o node_mod.o mem_jules.o \
	gridbox_mean_mod.o fluxes.o ancil_info.o csigma_mod.o mem_leaf.o \
	jules_surface_types_mod.o datetime_mod.o mem_chem1.o jules_fields_mod.o \
	ModTurbFields.o io_params.o ModLeaf3Init.o mem_micro.o ModMicControl.o \
	mem_grid.o mem_brams_jules.o model_time_mod.o gridmean_fluxes.o rconstants.o \
	chem1_list.o sf_diags_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

shcu_vars_const.o : $(CUPARM)/shcu_vars_const.f90 ModNamelistFile.o \
	ModConvComs.o grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSoilMoisture.o : $(SOIL_MOISTURE)/ModSoilMoisture.F90 ModNamelistFile.o \
	ModLeafComs.o node_mod.o ModControlVars.o ReadBcst.o ModMPassFull.o mem_grid.o \
	ModTurbFields.o memSoilMoisture.o io_params.o rconstants.o mem_leaf.o parlibf.o \
	mem_aerad.o ModBasicFields.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
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

ModTKenn.o : $(STILT)/ModTKenn.f90 mem_scratch.o mem_stilt.o mem_grid.o \
	turb_constants.o rconstants.o 
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

ModVarTables.o : $(MEMORY)/ModVarTables.f90 ModParallelEnvironment.o aer1_list.o \
	chem1_list.o io_params.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfUpdate.o : $(FDDA)/ModVarfUpdate.f90 mem_scratch.o ModInitHis.o \
	mem_grid.o rconstants.o ref_sounding.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

include jules_depend_model.mk

