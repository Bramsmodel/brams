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

ModAsTi.o : $(ISAN)/ModAsTi.f90 isan_coms.o rconstants.o mem_grid.o \
	ModChemAObj.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsTp.o : $(ISAN)/ModAsTp.f90 rconstants.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAVarF.o : $(ISAN)/ModAVarF.f90 isan_coms.o rconstants.o mem_grid.o ModRbnd.o 
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
	isan_coms.o mem_aer1.o ModChemFirstRams.o mem_chem1.o ModAsTi.o ModChemAsti2.o \
	ModChemVInterps.o ModChemAObj.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsti2.o : $(ISAN_CHEM)/ModChemAsti2.f90 ModDateUtils.o rconstants.o \
	isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAstp.o : $(ISAN_CHEM)/ModChemAstp.F90 dump.o chem1_list.o \
	chem_isan_coms.o ModDateUtils.o isan_coms.o mem_varinit.o mem_chem1.o \
	rconstants.o ModAsTp.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAvarf.o : $(ISAN_CHEM)/ModChemAvarf.f90 mem_grid.o chem_isan_coms.o \
	isan_coms.o ModRbnd.o mem_aer1.o mem_chem1.o rconstants.o ModAVarF.o \
	ModNestFeed.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_dry_dep.o : $(MODEL_CHEM)/chem_dry_dep.f90 chem1_list.o ModDateUtils.o \
	mem_aer1.o mem_chem1.o extra.o aer1_list.o 
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

chem_fastjx_driv.o : $(CCATT)/chem_fastjx_driv.f90 rconstants.o chem1_list.o \
	chem_fastjx57.o chem_fastjx_data.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFileInv.o : $(ISAN_CHEM)/ModChemFileInv.f90 ModDateUtils.o dump.o \
	isan_coms.o mem_grid.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFirstRams.o : $(ISAN_CHEM)/ModChemFirstRams.f90 mem_grid.o ModRcio.o \
	ModGetVar.o isan_coms.o ModNestFillDens.o rconstants.o mem_scratch.o \
	ModChemRefState.o grid_dims.o an_header.o $(UTILS_INCS)/files.h \
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

chem_orage.o : $(CCATT)/chem_orage.f90 mem_scratch1_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_plumerise_scalar.o : $(CCATT)/chem_plumerise_scalar.f90 mem_chem1.o \
	mem_plume_chem1.o node_mod.o mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemRefState.o : $(ISAN_CHEM)/ModChemRefState.f90 rconstants.o ccatt_start.o \
	ModNestFillDens.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_sources.o : $(CCATT)/chem_sources.f90 mem_grid.o parlibf.o ModDateUtils.o \
	mem_plume_chem1.o mem_volc_chem1.o mem_aer1.o mem_chem1.o ReadBcst.o \
	ModControlVars.o aer1_list.o io_params.o ModNamelistFile.o \
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

chem_spack_qssa.o : $(CCATT)/chem_spack_qssa.f90 chem_spack_fexloss.o \
	chem_spack_fexprod.o mem_chem1.o chem_spack_kinetic.o chem_spack_rates.o \
	chem_spack_dratedc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rates.o : $(MODEL_CHEM)/chem_spack_rates.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rodas3_dyndt.o : $(CCATT)/chem_spack_rodas3_dyndt.f90 \
	chem_spack_jacdchemdc.o mem_grid.o mem_spack.o mem_chem1.o chem_spack_kinetic.o \
	chem_spack_fexchem.o extra.o chem_spack_ros.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros.o : $(CCATT)/chem_spack_ros.f90 chem_spack_jacdchemdc.o \
	mem_spack.o mem_chem1.o chem_spack_solve_sparse.o chem_spack_kinetic.o \
	chem_spack_fexchem.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros_dyndt.o : $(CCATT)/chem_spack_ros_dyndt.f90 \
	chem_spack_jacdchemdc.o mem_spack.o mem_chem1.o chem_spack_solve_sparse.o \
	chem_spack_kinetic.o chem_spack_fexchem.o chem_spack_ros.o 
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

chem_trans_liq.o : $(CCATT)/chem_trans_liq.f90 mem_chem1.o mem_chem1aq.o 
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
	ModBasicFields.o ModMicroFields.o mem_aer1.o chem_dry_dep.o mem_radiate.o \
	rconstants.o mem_cuparm.o ModTurbFields.o mem_chem1.o grid_dims.o \
	ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemSourcesDriver.o : $(CCATT)/ChemSourcesDriver.f90 mem_stilt.o chem1_list.o \
	mem_leaf.o mem_plume_chem1.o mem_aer1.o mem_volc_chem1.o mem_chem1.o \
	chem_plumerise_scalar.o chem_sources.o aer1_list.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

coag.o : $(MATRIX)/coag.f90 memMatrix.o setup.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModCondRead.o : $(FDDA)/ModCondRead.f90 mem_grid.o ModDateUtils.o isan_coms.o \
	ModRamsGrid.o mem_varinit.o ModNudAnalysis.o ModCondUpdate.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCondUpdate.o : $(FDDA)/ModCondUpdate.f90 ModRcio.o mem_grid.o ModVarTables.o \
	ModInitHis.o grid_struct.o mem_varinit.o an_header.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModConvComs.o : $(CUPARM)/ModConvComs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ConvPar_GF_GEOS5.o : $(CUPARM)/ConvPar_GF_GEOS5.F90 module_gate.o \
	MAPL_Constants.o Henrys_Law_cts.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuRead.o : $(CUPARM)/ModCuRead.f90 ModDateUtils.o mem_cuparm.o mem_grid.o \
	isan_coms.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
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

ModCupGrellCattDeep.o : $(CUPARM)/ModCupGrellCattDeep.f90 kbcon_ecmwf.o \
	mem_carma.o Phys_const.o ModCupUp.o mem_grid.o cup_output_vars.o ModCupDn.o \
	ModCupEnv.o mem_varinit.o ccatt_start.o mem_grell_param2.o node_mod.o \
	mem_scratch2_grell.o ModCupEnvCatt.o mem_scratch3_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCupGrellCattShallow.o : $(CUPARM)/ModCupGrellCattShallow.f90 \
	mem_scratch3_grell_sh.o ModCupUp.o mem_grid.o cup_output_vars.o ModCupEnv.o \
	mem_varinit.o mem_scratch2_grell_sh.o mem_grell_param2.o node_mod.o \
	Phys_const.o ModCupEnvCatt.o 
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

dam.o : $(ENERGY)/dam.f90 ModDateUtils.o dump.o mem_grid.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

depv.o : $(MATRIX)/depv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

digitalFilter.o : $(MODEL)/digitalFilter.f90 mem_grid.o ModVarTables.o \
	ModDateUtils.o ModBasicFields.o utilsMod.o ReadBcst.o node_mod.o \
	ModControlVars.o grid_dims.o io_params.o ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
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
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

gfdl_cloud_microphys.o : $(MICRO)/gfdl_cloud_microphys.F90 module_mp_radar.o \
	node_mod.o 
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

ModGridSet.o : $(INIT)/ModGridSet.f90 rconstants.o grid_dims.o mem_grid.o 
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

ModLeaf3Hyd.o : $(SURFACE)/ModLeaf3Hyd.f90 mem_leaf.o ModLeafComs.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Init.o : $(SURFACE)/ModLeaf3Init.f90 ModLeaf3.o mem_grid.o \
	teb_spm_start.o mem_leaf.o rconstants.o ModLeafComs.o grid_dims.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Teb.o : $(SURFACE)/ModLeaf3Teb.f90 mem_teb_vars_const.o ModUrban.o \
	ModGasPart.o mem_emiss.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeafComs.o : $(SURFACE)/ModLeafComs.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

local_proc.o : $(MODEL)/local_proc.F90 mem_stilt.o dump.o mem_grid.o \
	rconstants.o ReadBcst.o node_mod.o ref_sounding.o grid_dims.o io_params.o \
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
	mcica_random_numbers.o rrlw_wvn.o rrlw_con.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_sw.o : $(RRTMG_SW_SRC)/mcica_subcol_gen_sw.f90 \
	mcica_random_numbers.o rrsw_vsn.o parrrsw.o rrsw_con.o parkind.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aer1.o : $(CCATT)/mem_aer1.f90 dump.o VarTable.o mem_grid.o mem_chem1.o \
	grid_dims.o node_mod.o ModScalarTable.o aer1_list.o io_params.o \
	ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aerad.o : $(RADIATE)/mem_aerad.f90 mem_grid_dim_defs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_carma.o : $(RADIATE)/mem_carma.f90 mem_aerad.o mem_grid.o parlibf.o \
	ModBasicFields.o ModNamelistFile.o ModRamsGrid.o grid_dims.o ModTurbFields.o \
	ReadBcst.o node_mod.o mem_scalar.o ModMPassFull.o mem_globrad.o \
	ModControlVars.o VarTable.o io_params.o ModSoilMoisture.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1.o : $(CCATT)/mem_chem1.f90 chem1_list.o grid_dims.o ModScalarTable.o \
	VarTable.o io_params.o ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1aq.o : $(CCATT)/mem_chem1aq.f90 ModVarTables.o mem_chem1.o \
	chem1aq_list.o ModScalarTable.o grid_dims.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chemic.o : $(CCATT)/mem_chemic.f90 ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_cuparm.o : $(CUPARM)/mem_cuparm.f90 VarTable.o grid_dims.o ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
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

GaspartFields.o : $(TEB_SPM)/GaspartFields.f90 ModNodeDimensions.o \
	ModParallelEnvironment.o VarTable.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globaer.o : $(RADIATE)/mem_globaer.f90 mem_precision.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globrad.o : $(RADIATE)/mem_globrad.f90 mem_precision.o mem_aerad.o parlibf.o \
	mem_radiate.o ModNamelistFile.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell.o : $(CUPARM)/mem_grell.f90 mem_cuparm.o shcu_vars_const.o \
	ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell_param2.o : $(CUPARM)/mem_grell_param2.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid.o : $(MEMORY)/mem_grid.f90 grid_dims.o VarTable.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid_dim_defs.o : $(MEMORY)/mem_grid_dim_defs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

JulesFields.o : $(SURFACE)/JulesFields.f90 ModNodeDimensions.o \
	ModParallelEnvironment.o VarTable.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_leaf.o : $(SURFACE)/mem_leaf.f90 teb_spm_start.o grid_dims.o VarTable.o \
	io_params.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicroFields.o : $(MICRO)/ModMicroFields.f90 ModNodeDimensions.o \
	ModParallelEnvironment.o VarTable.o ModNamelistFile.o ModMicControl.o 
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

mem_oda.o : $(FDDA)/mem_oda.f90 grid_dims.o ModNamelistFile.o ModVarTables.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_opt_scratch.o : $(TURB)/mem_opt_scratch.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_plume_chem1.o : $(CCATT)/mem_plume_chem1.f90 mem_chem1.o chem1_list.o \
	ModNamelistFile.o ModVarTables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_precision.o : $(RADIATE)/mem_precision.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_radiate.o : $(RADIATE)/mem_radiate.f90 VarTable.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_rrtm.o : $(RADIATE)/mem_rrtm.f90 parrrtm.o chem1_list.o mem_grid.o parrrsw.o \
	mem_chem1.o rconstants.o node_mod.o rrtmg_sw_init.o parkind.o rrtmg_lw_init.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scalar.o : $(MEMORY)/mem_scalar.f90 VarTable.o io_params.o ModNamelistFile.o 
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

mem_scratch1_grell.o : $(CUPARM)/mem_scratch1_grell.f90 ccatt_start.o dump.o \
	mem_grell_param2.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch2_grell.o : $(CUPARM)/mem_scratch2_grell.f90 mem_grell_param2.o \
	node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch2_grell_sh.o : $(CUPARM)/mem_scratch2_grell_sh.f90 mem_grell_param2.o \
	node_mod.o 
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

ShcuFields.o : $(CUPARM)/ShcuFields.f90 ModNodeDimensions.o \
	ModParallelEnvironment.o ModControlVars.o VarTable.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_spack.o : $(CCATT)/mem_spack.f90 chem1_list.o chem_spack_utils.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_stilt.o : $(STILT)/mem_stilt.f90 rconstants.o grid_dims.o VarTable.o \
	io_params.o ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tconv.o : $(CCATT)/mem_tconv.f90 chem1_list.o aer1_list.o mem_aer1.o 
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

mem_tend.o : $(MEMORY)/mem_tend.f90 mem_grid.o teb_spm_start.o ModBasicFields.o \
	ModMicroFields.o ModTurbFields.o mem_scalar.o mem_emiss.o ModScalarTable.o \
	ModNamelistFile.o GaspartFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_turb_scalar.o : $(TURB)/mem_turb_scalar.f90 VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tuv.o : $(TUV)/mem_tuv.f90 mem_stilt.o ModTuv2.7.o mem_globrad.o \
	ModVarTables.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mem_varinit.o : $(MEMORY)/mem_varinit.f90 chem1_list.o mem_chem1.o grid_dims.o \
	VarTable.o ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_volc_chem1.o : $(CCATT)/mem_volc_chem1.f90 ModNamelistFile.o ModVarTables.o \
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

meteogram.o : $(IO)/meteogram.f90 ModMPassDtl.o mem_grid.o ModVarTables.o \
	meteogramType.o satPolyColision.o node_mod.o ModPostUtils.o ModNamelistFile.o \
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

ModMicInit.o : $(MICRO)/ModMicInit.f90 dump.o mem_grid.o ModMicTabs.o \
	rconstants.o ModMicGamma.o ModMicControl.o $(UTILS_INCS)/files.h \
	$(MICRO)/MicConstants.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicNuc.o : $(MICRO)/ModMicNuc.f90 ModMicControl.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicTabs.o : $(MICRO)/ModMicTabs.f90 ModMicGamma.o ModMicControl.o \
	$(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicVap.o : $(MICRO)/ModMicVap.f90 rconstants.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicControl.o : $(MICRO)/ModMicControl.f90 ModParallelEnvironment.o \
	grid_dims.o ModNamelistFile.o $(UTILS_INCS)/files.h $(MICRO)/MicConstants.h 
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

ModAcoust.o : $(MODEL)/ModAcoust.f90 ModMessageSet.o mem_grid.o ModBasicFields.o \
	ModParallelEnvironment.o ModMicControl.o ModAcoustAdap.o rconstants.o \
	mem_scratch.o node_mod.o ref_sounding.o mem_tend.o ModGrid.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAerClim.o : $(AERCLIM)/ModAerClim.f90 dump.o mem_grid.o parlibf.o \
	ModBasicFields.o ModTurbFields.o ReadBcst.o node_mod.o ModControlVars.o \
	ModSoilMoisture.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBasicFields.o : $(MEMORY)/ModBasicFields.f90 mem_stilt.o ModNodeDimensions.o \
	ModParallelEnvironment.o VarTable.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBramsGrid.o : $(POST_SRC)/ModBramsGrid.f90 mem_aerad.o mem_grid.o \
	ModParallelEnvironment.o node_mod.o ModPostUtils.o ref_sounding.o \
	ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBuffering.o : $(MPI)/ModBuffering.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCarmaDriver.o : $(RADIATE)/ModCarmaDriver.f90 mem_carma.o mem_teb_common.o \
	mem_tend.o ModLeaf3.o mem_grid.o ModDateUtils.o teb_spm_start.o mem_leaf.o \
	ModBasicFields.o ModMicroFields.o mem_radiate.o rconstants.o mem_cuparm.o \
	node_mod.o rad_carma.o mem_scratch1_grell.o grid_dims.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsgen.o : $(ISAN_CHEM)/ModChemAsgen.F90 ModChemAstp.o ModMkSfcTop.o \
	ModDateUtils.o ModRamsGrid.o ModAsGen.o grid_dims.o chem1_list.o \
	ModChemIsanIo.o dump.o mem_grid.o isan_coms.o mem_aer1.o node_mod.o \
	ModChemRefState.o ModChemFileInv.o ModChemAsti.o ModChemAvarf.o \
	chem_isan_coms.o mem_chem1.o ModControlVars.o aer1_list.o io_params.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemConvTransp.o : $(CCATT)/ModChemConvTransp.f90 chem1_list.o mem_grid.o \
	mem_aer1.o mem_chem1.o mem_scratch.o mem_grell_param2.o mem_cuparm.o node_mod.o \
	Phys_const.o mem_scratch1_grell.o aer1_list.o mem_tconv.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemistryDriver.o : $(CCATT)/ModChemistryDriver.f90 mem_stilt.o mem_aerad.o \
	chem_spack_qssa.o rconstants.o chem_uv_att.o carma_fastjx.o mem_scratch.o \
	chem_fastjx_driv.o mem_rrtm.o chem_trans_liq.o chem_spack_rodas3_dyndt.o \
	extra.o mem_chemic.o chem_spack_ros.o chem_orage.o chem_trans_gasaq.o \
	ModTuvDriver2.7.o mem_spack.o ModBasicFields.o chem_spack_utils.o grid_dims.o \
	parrrtm.o chem1_list.o mem_grid.o mem_aer1.o node_mod.o mem_scratch1_grell.o \
	mem_carma.o chem_spack_ros_dyndt.o mem_chem1aq.o ModMicroFields.o mem_chem1.o \
	mem_radiate.o mem_cuparm.o mem_grell_param2.o chem_spack_solve_sparse.o \
	chem1aq_list.o aer1_list.o mem_globrad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModControlVars.o : $(INIT)/ModControlVars.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCoriolis.o : $(MODEL)/ModCoriolis.f90 mem_grid.o parlibf.o ModBasicFields.o \
	rconstants.o mem_scratch.o node_mod.o ModBuffering.o ref_sounding.o mem_tend.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuParGrell3.o : $(CUPARM)/ModCuParGrell3.F90 mem_stilt.o mem_varinit.o \
	rconstants.o mem_scratch.o ModRadvc.o module_cu_gf.o ModMicControl.o \
	ModVarTables.o ModBasicFields.o ModChemConvTransp.o module_cu_gf_v5.1.o \
	grid_dims.o mem_grell.o ModMessageSet.o mem_grid.o mem_leaf.o ccatt_start.o \
	node_mod.o mem_scratch1_grell.o ModNamelistFile.o mem_carma.o module_cu_g3.o \
	ModMicroFields.o mem_radiate.o mem_chem1.o mem_cuparm.o mem_grell_param2.o \
	ConvPar_GF_GEOS5.o Phys_const.o mem_tend.o ModRstilt.o io_params.o ModGrid.o \
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

ModDiffuse.o : $(TURB)/ModDiffuse.f90 mem_tend.o mem_grid.o mem_leaf.o \
	ModBasicFields.o ModMicroFields.o ModTurbK.o ModDiffSclr.o mem_scratch.o \
	ModTurbFields.o node_mod.o mem_opt_scratch.o ModTurbDiff.o ke_coms.o \
	ModTurbKE.o ModScalarTable.o ModNamelistFile.o ModMicControl.o \
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

ModGasPart.o : $(TEB_SPM)/ModGasPart.f90 ModRcio.o ModVarTables.o mem_grid.o \
	parlibf.o mem_leaf.o ModBasicFields.o node_mod.o mem_emiss.o \
	mem_teb_vars_const.o grid_dims.o an_header.o GaspartFields.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGeodat.o : $(MKSFC)/ModGeodat.f90 teb_spm_start.o mem_leaf.o io_params.o \
	mem_grid.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGrid.o : $(MPI)/ModGrid.F90 ModNodeDimensions.o ModMessageSet.o \
	ModVarTables.o ModNeighbourNodes.o ModParallelEnvironment.o ModBasicFields.o \
	ModScalarTable.o ModMicControl.o ModMicroFields.o meteogramType.o ShcuFields.o \
	ModControlVars.o ModDomainDecomp.o ModTurbFields.o ModGridDims.o JulesFields.o \
	mem_tend.o VarTable.o ModNamelistFile.o GaspartFields.o 
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

modIau.o : $(MODEL)/modIau.f90 dump.o mem_grid.o parlibf.o mem_varinit.o \
	ReadBcst.o node_mod.o mem_tend.o ModNamelistFile.o ModMPassFull.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitHis.o : $(IO)/ModInitHis.f90 mem_aerad.o ModLeaf3.o mem_varinit.o \
	rconstants.o mem_scratch.o an_header.o ModMicControl.o ModVarTables.o \
	ModBasicFields.o ModRamsReadHeader.o ModRamsGrid.o ModRinit.o chem1_list.o \
	ModRcio.o mem_grid.o mem_leaf.o ModGetVar.o ref_sounding.o mem_chem1.o \
	ModLeafComs.o io_params.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitMicThompson.o : $(MICRO)/ModInitMicThompson.f90 dump.o generic.o \
	mem_grid.o parlibf.o ModDateUtils.o ModBasicFields.o ModMicroFields.o \
	ReadBcst.o node_mod.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLanduseInput.o : $(MKSFC)/ModLanduseInput.f90 ModLeaf3Init.o mem_mksfc.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3.o : $(SURFACE)/ModLeaf3.f90 mem_teb_common.o mem_grid.o teb_spm_start.o \
	ModBasicFields.o mem_leaf.o ModMicroFields.o ModLeaf3Teb.o ModLeaf3Hyd.o \
	rconstants.o ModLeafComs.o ccatt_start.o mem_cuparm.o node_mod.o mem_scratch.o \
	ModTurbFields.o mem_radiate.o mem_teb.o io_params.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3OceanOnly.o : $(SURFACE)/ModLeaf3OceanOnly.f90 ModLeaf3.o mem_grid.o \
	mem_leaf.o ModBasicFields.o ModCuParGrell3.o mem_radiate.o rconstants.o \
	mem_cuparm.o ModTurbFields.o ModLeafComs.o node_mod.o ccatt_start.o \
	ConvPar_GF_GEOS5.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMatrixDriver.o : $(MATRIX)/ModMatrixDriver.F90 setup.o chem1_list.o \
	mem_grid.o subs.o mem_leaf.o ModBasicFields.o ModMicroFields.o mem_aer1.o \
	mem_chem1.o rconstants.o mem_radiate.o ModTurbFields.o isrpia.o node_mod.o \
	npf.o memMatrix.o coag.o aer1_list.o ModParticle.o ModMicControl.o 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ModMemAlloc.o : $(MEMORY)/ModMemAlloc.F90 mem_stilt.o mem_nestb.o \
	mem_teb_common.o mem_aerad.o ModGrid.o mem_volc_chem1.o mem_varinit.o mem_teb.o \
	mem_scratch.o mem_scratch2_grell_sh.o carma_fastjx.o mem_scalar.o \
	mem_opt_scratch.o mem_emiss.o extra.o mem_chemic.o mem_scratch3_grell_sh.o \
	ModEvaluation.o machine_arq.o ModVarTables.o ModParallelEnvironment.o \
	ModBasicFields.o chem_dry_dep.o ShcuFields.o mem_turb_scalar.o \
	shcu_vars_const.o chem_sources.o mem_globrad.o grid_dims.o GaspartFields.o \
	mem_scratch3_grell.o mem_grell.o chem1_list.o digitalFilter.o mem_grid.o \
	mem_leaf.o parrrsw.o ModCuParGrell3.o mem_aer1.o ccatt_start.o mem_oda.o \
	ModTurbFields.o mem_tuv.o node_mod.o mem_scratch2_grell.o mem_teb_vars_const.o \
	mem_scratch1_grell.o mem_carma.o mem_chem1aq.o teb_spm_start.o mem_globaer.o \
	ModMicroFields.o mem_plume_chem1.o mem_radiate.o mem_chem1.o ModLeafComs.o \
	mem_cuparm.o mem_grell_param2.o mem_scratch1_brams.o ModTuv2.7.o ModOptical.o \
	chem1aq_list.o JulesFields.o mem_tend.o mem_grid_dim_defs.o aer1_list.o \
	io_params.o modIau.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMemory.o : $(UTILS_LIB)/ModMemory.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageData.o : $(MPI)/ModMessageData.f90 ModParallelEnvironment.o \
	ModFieldSection.o ModFieldSectionList.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageSet.o : $(MPI)/ModMessageSet.f90 ModNodeDimensions.o mem_grid.o \
	ModVarTables.o parlibf.o ModNeighbourNodes.o ModParallelEnvironment.o \
	ModFieldSectionList.o ModMessageData.o ModDomainDecomp.o ModFieldSection.o \
	ModGridDims.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicGfdlDriver.o : $(MICRO)/ModMicGfdlDriver.f90 mem_grid.o mem_leaf.o \
	ModBasicFields.o gfdl_cloud_microphys.o ModMicroFields.o mem_radiate.o \
	rconstants.o node_mod.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsDrive.o : $(MICRO)/ModMicrophysicsDrive.f90 ModMicNuc.o \
	mem_grid.o mem_chem1aq.o ModMicVap.o ModMicrophysicsMisc.o ModBasicFields.o \
	ModMicColl.o ModMicroFields.o ModMicTabs.o mem_radiate.o mem_chem1.o node_mod.o \
	ModMicInit.o grid_dims.o mem_chemic.o ModMicControl.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsMisc.o : $(MICRO)/ModMicrophysicsMisc.f90 mem_grid.o \
	ModBasicFields.o ModMicroFields.o rconstants.o mem_scratch.o ModMicControl.o \
	$(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicThompsonDriver.o : $(MICRO)/ModMicThompsonDriver.f90 mem_grid.o \
	module_mp_thompson.o mem_leaf.o ModBasicFields.o ModMicroFields.o mem_radiate.o \
	rconstants.o node_mod.o io_params.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcDriver.o : $(MKSFC)/ModMkSfcDriver.f90 ModLanduseInput.o ModMkSfcNdvi.o \
	ModMkSfcTop.o mem_grid.o mem_mksfc.o teb_spm_start.o ModSstRead.o \
	ModNestGeoSst.o ModNdviRead.o ModMkSfcFuso.o ReadBcst.o node_mod.o \
	ModMkSfcSfc.o ModControlVars.o ModMkSfcSst.o grid_dims.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcFuso.o : $(MKSFC)/ModMkSfcFuso.f90 mem_teb_vars_const.o mem_grid.o \
	mem_mksfc.o mem_teb.o ReadBcst.o node_mod.o mem_emiss.o ModControlVars.o \
	io_params.o GaspartFields.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcNdvi.o : $(MKSFC)/ModMkSfcNdvi.f90 ModLanduseInput.o dump.o mem_grid.o \
	mem_mksfc.o mem_leaf.o ModRUser.o io_params.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSfc.o : $(MKSFC)/ModMkSfcSfc.f90 dump.o mem_grid.o mem_mksfc.o \
	mem_leaf.o ReadBcst.o node_mod.o ModControlVars.o io_params.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSst.o : $(MKSFC)/ModMkSfcSst.f90 mem_grid.o mem_mksfc.o ModGeodat.o \
	mem_leaf.o ModRUser.o ModNestFillDens.o grid_dims.o io_params.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcTop.o : $(MKSFC)/ModMkSfcTop.f90 dump.o mem_grid.o mem_mksfc.o \
	ReadBcst.o node_mod.o ModControlVars.o io_params.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMonotonicAdvection.o : $(MODEL)/ModMonotonicAdvection.f90 ModGrid.o \
	mem_grid.o ModMessageSet.o ModParallelEnvironment.o mem_aer1.o chem_dry_dep.o \
	mem_chem1.o rconstants.o ccatt_start.o ModDomainDecomp.o ModNamelistFile.o \
	ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNamelistFile.o : $(INIT)/ModNamelistFile.f90 dump.o parlibf.o \
	ModParallelEnvironment.o modPrintInitial.o grid_dims.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNdviRead.o : $(MKSFC)/ModNdviRead.f90 ModMkSfcNdvi.o mem_mksfc.o mem_grid.o \
	ModDateUtils.o mem_leaf.o ReadBcst.o node_mod.o ModControlVars.o grid_dims.o \
	io_params.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNeighbourNodes.o : $(MPI)/ModNeighbourNodes.f90 ModParallelEnvironment.o \
	ModDomainDecomp.o ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestGeoSst.o : $(MKSFC)/ModNestGeoSst.f90 ModLanduseInput.o ModGeodat.o \
	mem_scratch.o ModNestFeed.o ModMkSfcTop.o ModInitHis.o ModBasicFields.o \
	grid_dims.o ModSoilMoisture.o dump.o mem_grid.o mem_leaf.o ModLeaf3Init.o \
	ccatt_start.o ModControlVars.o ModTurbFields.o node_mod.o mem_mksfc.o \
	ModRUser.o ModNestFillDens.o memSoilMoisture.o io_params.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNodeDimensions.o : $(MPI)/ModNodeDimensions.f90 ModParallelEnvironment.o \
	ModDomainDecomp.o ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudAnalysis.o : $(FDDA)/ModNudAnalysis.f90 dump.o chem1_list.o \
	ModEvaluation.o mem_grid.o ModBasicFields.o ModNestFillDens.o mem_varinit.o \
	mem_chem1.o mem_scratch.o node_mod.o mem_tend.o modIau.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaNudge.o : $(FDDA)/ModOdaNudge.f90 ModOdaProcObs.o mem_grid.o \
	ModBasicFields.o mem_oda.o mem_scratch.o node_mod.o mem_tend.o io_params.o \
	ModOdaKrig.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOneProc.o : $(MODEL)/ModOneProc.F90 ModSched.o ModAerClim.o ModTimeStamp.o \
	ModMemAlloc.o ModMicInit.o meteogram.o ModTuvDriver2.7.o ModEvaluation.o \
	ModMkSfcTop.o ModInitHis.o ModRinit.o ModNdviRead.o ModOpspec.o ref_sounding.o \
	ModRanlavg.o ModMkSfcDriver.o ModRThrm.o ModCoriolis.o mem_volc_chem1.o extra.o \
	machine_arq.o ModBasicFields.o ModRamsGrid.o ModChemAsgen.o grid_dims.o \
	ModSoilMoisture.o ModVarfFile.o mem_aer1.o ModUrbanCanopy.o ccatt_start.o \
	ModDomainDecomp.o node_mod.o ModCondRead.o ModInitMicThompson.o \
	ModNamelistFile.o mem_carma.o ModPostProcess.o mem_chem1aq.o teb_spm_start.o \
	ModSstRead.o ModLeaf3Teb.o mem_chem1.o ModChemistryDriver.o ModMkSfcSfc.o \
	ModRio.o memSoilMoisture.o aer1_list.o io_params.o mem_globrad.o modIau.o \
	ModRnode.o ModGrid.o ModMonotonicAdvection.o ModPostGridNetCDF.o \
	ModNestGeoSst.o mem_scalar.o dam.o ModCuRead.o ModParallelEnvironment.o \
	ModMkSfcFuso.o ModOdaRead.o ModRecycle.o ModMPassDtl.o chem1_list.o \
	digitalFilter.o tuvParameter.o mem_leaf.o ModCuParGrell3.o ModParaInit.o \
	mem_oda.o mem_teb_vars_const.o ModTimestepRK.o domain_decomp.o \
	mem_plume_chem1.o mem_cuparm.o ModWindFarm.o ModRhhi.o mem_stilt.o \
	mem_teb_common.o parlibf.o ModTimestep.o mem_varinit.o mem_teb.o mem_scratch.o \
	ReadBcst.o mem_emiss.o ModVarTables.o ModMicrophysicsMisc.o ModGridTree.o \
	chem_dry_dep.o shcu_vars_const.o chem_sources.o ModGasPart.o local_proc.o \
	dump.o mem_grid.o ModNudRead.o ModLeaf3Init.o ModVarfUpdate.o isan_coms.o \
	ModNestIntrp.o ModRUser.o ModRamsMicrophysics2M.o mem_radiate.o \
	mem_grell_param2.o ModTuv2.7.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOpspec.o : $(IO)/ModOpspec.f90 mem_stilt.o io_params.o mem_varinit.o \
	mem_emiss.o ModMicControl.o chem_sources.o shcu_vars_const.o grid_dims.o \
	chem1_list.o mem_grid.o mem_leaf.o mem_aer1.o ccatt_start.o ModNamelistFile.o \
	mem_chem1aq.o teb_spm_start.o mem_chem1.o mem_radiate.o mem_cuparm.o \
	mem_grell_param2.o chem1aq_list.o aer1_list.o mem_globrad.o modIau.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOptical.o : $(RADIATE)/ModOptical.f90 dump.o mem_grid.o ModVarTables.o \
	parlibf.o mem_leaf.o ModBasicFields.o ModRamsGrid.o mem_aer1.o mem_radiate.o \
	ccatt_start.o ModTurbFields.o ReadBcst.o node_mod.o ModMPassFull.o \
	ModControlVars.o aer1_list.o ModNamelistFile.o ModSoilMoisture.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOutputUtils.o : $(IO)/ModOutputUtils.f90 dump.o ModVarTables.o \
	ModBasicFields.o ModTurbFields.o ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOzone.o : $(TEB_SPM)/ModOzone.f90 mem_grid.o ModBasicFields.o mem_radiate.o \
	rconstants.o GaspartFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParaInit.o : $(MPI)/ModParaInit.f90 dump.o mem_grid.o ModVarTables.o \
	node_mod.o ModScalarTable.o grid_dims.o $(UTILS_INCS)/constants.h 
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

ModPostGrid.o : $(POST_SRC)/ModPostGrid.F90 mem_grid.o parlibf.o ModBramsGrid.o \
	ModBasicFields.o ModParallelEnvironment.o ModOutputUtils.o io_params.o \
	ModTurbFields.o ModPostUtils.o ModPostOneFieldNetCDF.o ModPostTypes.o \
	ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostGridNetCDF.o : $(POST_SRC)/ModPostGridNetCDF.F90 dump.o mem_grid.o \
	ModBramsGrid.o ModDateUtils.o ModNamelistFile.o ModPostTypes.o ModPostUtils.o \
	io_params.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField.o : $(POST_SRC)/ModPostOneField.f90 ModPostOneField7d.o dump.o \
	ModBramsGrid.o ModBasicFields.o ModPostOneFieldUtils.o ModPostTypes.o \
	ModTurbFields.o ModPostOneField8d.o node_mod.o ModPostOneField2d.o \
	ModPostUtils.o ModPostOneField3d.o ModNamelistFile.o ModMicControl.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField2d.o : $(POST_SRC)/ModPostOneField2d.f90 dump.o mem_aerad.o \
	mem_grid.o ModBramsGrid.o ModBasicFields.o ModOutputUtils.o \
	ModPostOneFieldUtils.o io_params.o mem_radiate.o mem_cuparm.o ModTurbFields.o \
	node_mod.o ModPostUtils.o ModPostGrid.o ModPostTypes.o ModNamelistFile.o \
	ModMicControl.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField3d.o : $(POST_SRC)/ModPostOneField3d.f90 mem_grid.o \
	ModBramsGrid.o ModBasicFields.o ModOutputUtils.o ModPostOneFieldUtils.o \
	mem_varinit.o ModTurbFields.o node_mod.o ModPostUtils.o ModPostGrid.o \
	ModPostTypes.o ModNamelistFile.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField7d.o : $(POST_SRC)/ModPostOneField7d.f90 ModBramsGrid.o \
	ModBasicFields.o ModOutputUtils.o ModPostOneFieldUtils.o ModTurbFields.o \
	ModPostUtils.o ModPostTypes.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField8d.o : $(POST_SRC)/ModPostOneField8d.f90 ModBramsGrid.o \
	ModBasicFields.o ModOutputUtils.o ModPostOneFieldUtils.o ModTurbFields.o \
	ModPostUtils.o ModPostTypes.o ModNamelistFile.o 
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

ModPostProcess.o : $(POST_SRC)/ModPostProcess.F90 ModMessageSet.o \
	ModPostOneField.o ModBramsGrid.o ModPostGridNetCDF.o ModBasicFields.o \
	ModParallelEnvironment.o ModGridTree.o ModPostTypes.o ModNamelistFile.o \
	ModTurbFields.o ModPostGrid.o io_params.o ModGrid.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostTypes.o : $(POST_SRC)/ModPostTypes.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostUtils.o : $(POST_SRC)/ModPostUtils.f90 mem_leaf.o \
	ModParallelEnvironment.o dump.o $(UTILS_INCS)/files.h \
	$(POST_INCS)/post_rconstants.h $(POST_INCS)/post_rconfig.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modPrintInitial.o : $(INIT)/modPrintInitial.F90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvc.o : $(MODEL)/ModRadvc.f90 ModRadvcAdap.o ModMonotonicAdvection.o \
	mem_grid.o ModBasicFields.o ModParallelEnvironment.o ModScalarTable.o \
	mem_aer1.o chem_dry_dep.o mem_chem1.o ccatt_start.o mem_scratch.o mem_tend.o \
	grid_dims.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcRK.o : $(MODEL)/ModRadvcRK.f90 mem_stilt.o ModMessageSet.o mem_grid.o \
	ModRexev.o ModParallelEnvironment.o mem_chem1.o node_mod.o mem_tend.o \
	grid_dims.o ModGrid.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsMicrophysics2M.o : $(MICRO)/ModRamsMicrophysics2M.f90 dump.o mem_grid.o \
	mem_leaf.o ModBasicFields.o ModMicroFields.o rconstants.o mem_scratch.o \
	node_mod.o ModMicGamma.o grid_dims.o ModMicControl.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsReadHeader.o : $(IO)/ModRamsReadHeader.f90 an_header.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRanlavg.o : $(IO)/ModRanlavg.f90 mem_grid.o ModVarTables.o ModBasicFields.o \
	ModRThrm.o ModMicroFields.o grid_dims.o io_params.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRbnd.o : $(BC)/ModRbnd.f90 mem_tend.o mem_grid.o ModBasicFields.o \
	ModMicrophysicsMisc.o ModMicroFields.o mem_chem1.o ccatt_start.o mem_scratch.o \
	ModTurbFields.o node_mod.o ref_sounding.o ModTurbKE.o ModScalarTable.o \
	ModMicControl.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRcio.o : $(IO)/ModRcio.f90 mem_stilt.o mem_grid.o mem_leaf.o \
	ModParallelEnvironment.o ModNamelistFile.o mem_radiate.o mem_cuparm.o \
	ModLeafComs.o ref_sounding.o grid_dims.o io_params.o an_header.o \
	ModMicControl.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConv.o : $(CUPARM)/ModRConv.f90 mem_grid.o ModBasicFields.o ModConvComs.o \
	rconstants.o mem_cuparm.o mem_scratch.o node_mod.o mem_tend.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConvGrellCatt.o : $(CUPARM)/ModRConvGrellCatt.f90 mem_stilt.o rconstants.o \
	mem_scratch.o mem_scalar.o ModMicControl.o ModChemConvTransp.o \
	ModCupGrellCattShallow.o mem_grell.o mem_grid.o mem_leaf.o ModCuParGrell3.o \
	ccatt_start.o node_mod.o ModCupGrellCattDeep.o mem_scratch1_grell.o \
	mem_radiate.o mem_cuparm.o mem_grell_param2.o mem_tend.o ModRstilt.o \
	io_params.o ModGrid.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRecycle.o : $(IO)/ModRecycle.f90 dump.o chem1_list.o mem_aerad.o mem_grid.o \
	ModVarTables.o ModDateUtils.o ModRamsReadHeader.o ModGetVar.o mem_aer1.o \
	mem_chem1.o ReadBcst.o node_mod.o ModMPassFull.o aer1_list.o io_params.o \
	an_header.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRexev.o : $(STILT)/ModRexev.f90 mem_stilt.o mem_grid.o ModBasicFields.o \
	rconstants.o mem_scratch.o ModRadvc.o mem_tend.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRGrad.o : $(TURB)/ModRGrad.f90 mem_scratch.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRhhi.o : $(INIT)/ModRhhi.f90 mem_grid.o ModBasicFields.o ModRamsGrid.o \
	ModRinit.o mem_scratch.o rconstants.o ref_sounding.o grid_dims.o \
	ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRinit.o : $(INIT)/ModRinit.f90 mem_grid.o ModBasicFields.o ModMicroFields.o \
	ModRbnd.o mem_varinit.o rconstants.o mem_scratch.o ModTurbFields.o node_mod.o \
	ref_sounding.o ModTurbKE.o io_params.o ModMicControl.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRio.o : $(IO)/ModRio.f90 mem_aerad.o parlibf.o ReadBcst.o an_header.o \
	ModMicControl.o ModVarTables.o ModDateUtils.o ModParallelEnvironment.o \
	ModBasicFields.o grid_dims.o utilsMod.o mem_grid.o ModRcio.o ModTurbFields.o \
	node_mod.o ref_sounding.o mpi_io_engine-5d.o ModNamelistFile.o mem_chem1.o \
	ModControlVars.o io_params.o ModMPassFull.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/interface.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRrtmDriver.o : $(RADIATE)/ModRrtmDriver.f90 rconstants.o rrtmg_sw_rad.o \
	mem_rrtm.o rrtmg_lw_rad.o ModMicControl.o ModDateUtils.o ModBasicFields.o \
	parkind.o grid_dims.o parrrtm.o rrtmg_lw_cldprop.o mem_grid.o mem_leaf.o \
	parrrsw.o ccatt_start.o mem_tuv.o node_mod.o mcica_subcol_gen_sw.o \
	ref_sounding.o mem_scratch1_grell.o rrtmg_sw_cldprop.o mem_carma.o \
	mcica_subcol_gen_lw.o teb_spm_start.o ModMicroFields.o mem_radiate.o \
	mem_chem1.o mem_cuparm.o ModOptical.o mem_grell_param2.o ModLeafComs.o \
	mem_tend.o $(UTILS_INCS)/aerosol_setup.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRShCuPar.o : $(CUPARM)/ModRShCuPar.f90 ModRConv.o mem_grid.o ModBasicFields.o \
	ModMicroFields.o ModConvComs.o ShcuFields.o mem_scratch.o ModTurbFields.o \
	node_mod.o shcu_vars_const.o mem_tend.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRstilt.o : $(STILT)/ModRstilt.f90 mem_stilt.o ModMonotonicAdvection.o \
	mem_grid.o ModBasicFields.o mem_cuparm.o mem_scratch.o ModTurbFields.o \
	mem_scratch1_grell.o grid_dims.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRThrm.o : $(MODEL)/ModRThrm.f90 mem_grid.o ModBasicFields.o ModMicroFields.o \
	rconstants.o mem_scratch.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRtimi.o : $(MODEL)/ModRtimi.f90 mem_grell.o mem_tend.o mem_grid.o \
	ModBasicFields.o mem_cuparm.o mem_scratch.o node_mod.o shcu_vars_const.o \
	ModScalarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModScalarTable.o : $(MEMORY)/ModScalarTable.f90 ModParallelEnvironment.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSched.o : $(MODEL)/ModSched.f90 dump.o mem_grid.o parlibf.o ModBasicFields.o \
	ModNamelistFile.o isan_coms.o mem_varinit.o mem_radiate.o mem_cuparm.o \
	ReadBcst.o node_mod.o shcu_vars_const.o ref_sounding.o io_params.o local_proc.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSeaSalt.o : $(CCATT)/ModSeaSalt.f90 mem_grid.o mem_leaf.o ModBasicFields.o \
	ModAerClim.o mem_aer1.o mem_chem1.o ccatt_start.o node_mod.o aer1_list.o \
	io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSstRead.o : $(MKSFC)/ModSstRead.f90 mem_mksfc.o mem_grid.o ModDateUtils.o \
	mem_leaf.o ReadBcst.o node_mod.o ModControlVars.o ModMkSfcSst.o grid_dims.o \
	io_params.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimeStamp.o : $(MODEL)/ModTimeStamp.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestep.o : $(MODEL)/ModTimestep.F90 mem_stilt.o sfclyr_jules.o \
	ModMonotonicAdvection.o ModDiffuse.o ModLeaf3.o rad_driv.o ModMatrixDriver.o \
	ModMicrophysicsDrive.o mem_varinit.o ModTimeStamp.o rconstants.o ModRadvc.o \
	mem_scratch.o mem_scalar.o ModRShCuPar.o mem_emiss.o ModRConv.o machine_arq.o \
	ModBasicFields.o ModMicrophysicsMisc.o ModTurbK.o ModRtimi.o shcu_vars_const.o \
	chem_sources.o ModOdaNudge.o ModGasPart.o ModMicGfdlDriver.o ModOzone.o \
	grid_dims.o digitalFilter.o ModMessageSet.o mem_grid.o mem_leaf.o ModSeaSalt.o \
	ModCuParGrell3.o mem_aer1.o ModUrbanCanopy.o ccatt_start.o mem_oda.o \
	ChemDryDepDriver.o node_mod.o ModAcoust.o ChemSourcesDriver.o \
	ModMicThompsonDriver.o ModRexev.o ModRConvGrellCatt.o teb_spm_start.o \
	mem_plume_chem1.o ModRThrm.o ModRamsMicrophysics2M.o ModRbnd.o mem_chem1.o \
	mem_radiate.o mem_cuparm.o ModNudAnalysis.o ModChemistryDriver.o ModWindFarm.o \
	ModOptical.o mem_tend.o ModRstilt.o ModCoriolis.o ModGrid.o \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestepRK.o : $(MODEL)/ModTimestepRK.F90 mem_stilt.o sfclyr_jules.o \
	ModMonotonicAdvection.o ModDiffuse.o ModLeaf3.o rad_driv.o ModGrid.o \
	ModMatrixDriver.o ModAerClim.o ModMicrophysicsDrive.o ModTimestep.o \
	mem_varinit.o ModTimeStamp.o rconstants.o ModRadvc.o mem_scratch.o mem_scalar.o \
	ModRShCuPar.o mem_emiss.o ModRConv.o machine_arq.o ModParallelEnvironment.o \
	ModMicrophysicsMisc.o ModTurbK.o ModRtimi.o shcu_vars_const.o chem_sources.o \
	ModOdaNudge.o ModGasPart.o ModMicGfdlDriver.o ModOzone.o grid_dims.o utilsMod.o \
	digitalFilter.o ModMessageSet.o mem_grid.o mem_leaf.o ModSeaSalt.o \
	ModCuParGrell3.o mem_aer1.o ModUrbanCanopy.o ccatt_start.o mem_oda.o \
	ChemDryDepDriver.o ModLeaf3OceanOnly.o node_mod.o ModAcoust.o \
	ChemSourcesDriver.o ModMicThompsonDriver.o ModRexev.o teb_spm_start.o \
	mem_plume_chem1.o ModRThrm.o ModRbnd.o ModRamsMicrophysics2M.o mem_chem1.o \
	mem_radiate.o mem_cuparm.o ModNudAnalysis.o ModChemistryDriver.o ModWindFarm.o \
	ModOptical.o ModRadvcRK.o mem_tend.o ModRstilt.o ModCoriolis.o modIau.o \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiff.o : $(TURB)/ModTurbDiff.f90 mem_grid.o ModRGrad.o mem_scratch.o \
	mem_cuparm.o mem_opt_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiffAdap.o : $(TURB)/ModTurbDiffAdap.f90 mem_scratch.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbFields.o : $(TURB)/ModTurbFields.f90 ModNodeDimensions.o \
	ModParallelEnvironment.o VarTable.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbK.o : $(TURB)/ModTurbK.f90 mem_stilt.o ModMonotonicAdvection.o \
	rconstants.o mem_scratch.o ke_coms.o ModMicControl.o ModTKenn.o \
	ModBasicFields.o ModTurbKAdap.o mem_turb_scalar.o grid_dims.o mem_grell.o \
	mem_tend.o mem_grid.o mem_leaf.o ccatt_start.o ModTurbFields.o node_mod.o \
	ModTurbDiffAdap.o ModScalarTable.o ModTurbKE.o ModNamelistFile.o \
	ModMicroFields.o ModRGrad.o mem_chem1.o mem_cuparm.o ModTurbDiff.o ModRstilt.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKAdap.o : $(TURB)/ModTurbKAdap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKE.o : $(TURB)/ModTurbKE.f90 mem_grid.o rconstants.o mem_scratch.o \
	ModTurbFields.o ke_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuv2.7.o : $(TUV)/ModTuv2.7.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuvDriver2.7.o : $(TUV)/ModTuvDriver2.7.f90 mem_carma.o chem1_list.o \
	mem_aerad.o mem_grid.o tuvParameter.o mem_leaf.o ModBasicFields.o mem_radiate.o \
	rconstants.o mem_chem1.o ModTuv2.7.o mem_tuv.o node_mod.o chem_fastjx_driv.o \
	ref_sounding.o mem_rrtm.o extra.o mem_globrad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_cu_g3.o : $(CUPARM)/module_cu_g3.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_cu_gd_fim.o : $(CUPARM)/module_cu_gd_fim.f90 module_gate.o Phys_const.o 
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

module_mp_thompson.o : $(MICRO)/module_mp_thompson.f90 module_mp_radar.o \
	node_mod.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

module_wind_fitch.o : $(WIND_FARM)/module_wind_fitch.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModUrbanCanopy.o : $(SURFACE)/ModUrbanCanopy.f90 mem_grid.o ModBasicFields.o \
	ModTurbFields.o node_mod.o mem_tend.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfFile.o : $(FDDA)/ModVarfFile.f90 parlibf.o mem_varinit.o rconstants.o \
	mem_scratch.o ReadBcst.o ModMicControl.o ModDateUtils.o ModBasicFields.o \
	ModRamsReadHeader.o ModGridTree.o ModRamsGrid.o chem1_list.o ModRcio.o \
	mem_grid.o ModMessageSet.o mem_leaf.o ModGetVar.o ModVarfUpdate.o isan_coms.o \
	mem_aer1.o node_mod.o ref_sounding.o mem_chem1.o ModNudAnalysis.o \
	ModControlVars.o aer1_list.o ModGrid.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModWindFarm.o : $(WIND_FARM)/ModWindFarm.f90 mem_grid.o ModDateUtils.o \
	ModBasicFields.o rconstants.o ModTurbFields.o node_mod.o mem_tend.o \
	module_wind_fitch.o io_params.o ModNamelistFile.o $(UTILS_INCS)/files.h 
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

ModNestIntrp.o : $(NESTING)/ModNestIntrp.f90 mem_nestb.o mem_grid.o \
	ModBasicFields.o ModNestFillDens.o ModRinit.o mem_scratch.o rconstants.o \
	ref_sounding.o grid_dims.o 
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

ModNudRead.o : $(FDDA)/ModNudRead.f90 mem_grid.o ModDateUtils.o ModNudUpdate.o \
	ModRamsGrid.o mem_varinit.o isan_coms.o mem_chem1.o ModNudAnalysis.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudUpdate.o : $(FDDA)/ModNudUpdate.f90 chem1_list.o mem_aerad.o mem_grid.o \
	ModRcio.o ModVarTables.o ModInitHis.o grid_struct.o mem_varinit.o mem_chem1.o \
	an_header.o $(UTILS_INCS)/files.h 
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

ModOdaProcObs.o : $(FDDA)/ModOdaProcObs.f90 rconstants.o mem_oda.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaRead.o : $(FDDA)/ModOdaRead.f90 mem_grid.o ModDateUtils.o ModOdaStaCount.o \
	isan_coms.o mem_oda.o ModOdaStaInput.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaCount.o : $(FDDA)/ModOdaStaCount.f90 mem_oda.o ModReadRalph.o \
	mem_grid.o obs_input.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaInput.o : $(FDDA)/ModOdaStaInput.f90 ModReadRalph.o mem_grid.o \
	ModDateUtils.o ModOdaStaCount.o mem_oda.o obs_input.o 
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

ModAcoustAdap.o : $(MODEL)/ModAcoustAdap.f90 mem_grid.o ModMessageSet.o \
	ModRbnd.o rconstants.o mem_scratch.o node_mod.o ModGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_carma.o : $(RADIATE)/rad_carma.F90 mem_carma.o chem1_list.o mem_aerad.o \
	machine_arq.o mem_grid.o mem_leaf.o ModDateUtils.o mem_globaer.o mem_aer1.o \
	mem_radiate.o rconstants.o ccatt_start.o carma_fastjx.o mem_tuv.o mem_chem1.o \
	node_mod.o grid_dims.o aer1_list.o mem_globrad.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_driv.o : $(RADIATE)/rad_driv.f90 ModBasicFields.o ModMicroFields.o \
	ModRrtmDriver.o mem_radiate.o ModCarmaDriver.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcAdap.o : $(MODEL)/ModRadvcAdap.f90 ModAdapInit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsGrid.o : $(INIT)/ModRamsGrid.f90 dump.o mem_grid.o ModGridSet.o \
	rconstants.o node_mod.o ModAdapInit.o $(UTILS_INCS)/constants.h 
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

ReadBcst.o : $(MPI)/ReadBcst.f90 mem_aerad.o mem_grid.o ModVarTables.o parlibf.o \
	ModBasicFields.o ModTurbFields.o node_mod.o ModMPassFull.o ModControlVars.o \
	utilsMod.o an_header.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ref_sounding.o : $(MODEL)/ref_sounding.f90 grid_dims.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRnode.o : $(MODEL)/ModRnode.f90 mem_grid.o ModVarTables.o parlibf.o \
	mem_leaf.o node_mod.o grid_dims.o $(UTILS_INCS)/constants.h 
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

rrtmg_lw_cldprmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprmc.f90 parrrtm.o rrlw_wvn.o \
	parkind.o rrlw_cld.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprop.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprop.f90 parrrtm.o parkind.o \
	rrlw_cld.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_init.o : $(RRTMG_LW_SRC)/rrtmg_lw_init.f90 rrlw_tbl.o rrlw_kg01.o \
	rrlw_cld.o rrtmg_lw_k_g.o rrlw_kg04.o rrlw_kg08.o rrlw_kg16.o rrlw_kg11.o \
	rrlw_kg09.o parkind.o rrlw_con.o parrrtm.o rrlw_kg15.o rrlw_kg12.o rrlw_kg14.o \
	rrlw_kg13.o rrlw_kg06.o rrtmg_lw_setcoef.o rrlw_vsn.o rrlw_wvn.o rrlw_kg07.o \
	rrlw_kg02.o rrlw_kg03.o rrlw_kg05.o rrlw_kg10.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_k_g.o : $(RRTMG_LW_SRC)/rrtmg_lw_k_g.f90 rrlw_kg15.o rrlw_kg12.o \
	rrlw_kg08.o rrlw_kg14.o rrlw_kg04.o rrlw_kg16.o rrlw_kg01.o rrlw_kg13.o \
	rrlw_kg07.o rrlw_kg11.o rrlw_kg02.o rrlw_kg03.o parkind.o rrlw_kg09.o \
	rrlw_kg05.o rrlw_kg10.o rrlw_kg06.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rad.o : $(RRTMG_LW_SRC)/rrtmg_lw_rad.f90 mcica_subcol_gen_lw.o \
	parrrtm.o rrtmg_lw_taumol.o rrlw_wvn.o rrtmg_lw_cldprmc.o parkind.o \
	rrtmg_lw_rtrnmc.o rrlw_con.o rrtmg_lw_setcoef.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrn.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrn.f90 parrrtm.o rrlw_wvn.o \
	rrlw_con.o rrlw_tbl.o parkind.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmc.f90 parrrtm.o rrlw_wvn.o \
	rrlw_con.o rrlw_tbl.o parkind.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmr.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmr.f90 parrrtm.o rrlw_wvn.o \
	rrlw_con.o rrlw_tbl.o parkind.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_setcoef.o : $(RRTMG_LW_SRC)/rrtmg_lw_setcoef.f90 parrrtm.o rrlw_ref.o \
	rrlw_wvn.o parkind.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_taumol.o : $(RRTMG_LW_SRC)/rrtmg_lw_taumol.f90 rrlw_kg01.o rrlw_kg08.o \
	rrlw_kg16.o rrlw_kg04.o rrlw_kg11.o rrlw_kg09.o rrlw_con.o parkind.o parrrtm.o \
	rrlw_kg15.o rrlw_kg12.o rrlw_ref.o rrlw_kg14.o rrlw_kg13.o rrlw_kg06.o \
	rrlw_vsn.o rrlw_wvn.o rrlw_kg07.o rrlw_kg02.o rrlw_kg03.o rrlw_kg05.o \
	rrlw_kg10.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprmc.f90 rrsw_vsn.o parrrsw.o \
	parkind.o rrsw_cld.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprop.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprop.f90 rrsw_vsn.o parrrsw.o \
	parkind.o rrsw_cld.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_init.o : $(RRTMG_SW_SRC)/rrtmg_sw_init.f90 rrsw_kg21.o rrsw_vsn.o \
	rrsw_kg29.o rrsw_tbl.o rrsw_kg22.o rrsw_kg24.o rrsw_kg16.o rrtmg_sw_setcoef.o \
	rrsw_kg20.o rrsw_kg19.o rrsw_kg27.o parkind.o rrsw_cld.o rrsw_aer.o \
	rrtmg_sw_k_g.o rrsw_kg23.o parrrsw.o rrsw_con.o rrsw_kg17.o rrsw_kg28.o \
	rrsw_kg25.o rrsw_wvn.o rrsw_kg26.o rrsw_kg18.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_k_g.o : $(RRTMG_SW_SRC)/rrtmg_sw_k_g.f90 rrsw_kg21.o rrsw_kg19.o \
	rrsw_kg23.o rrsw_kg26.o rrsw_vsn.o rrsw_kg29.o rrsw_kg22.o rrsw_kg20.o \
	rrsw_kg18.o rrsw_kg17.o rrsw_kg28.o rrsw_kg24.o rrsw_kg27.o parkind.o \
	rrsw_kg16.o rrsw_kg25.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_rad.o : $(RRTMG_SW_SRC)/rrtmg_sw_rad.f90 rrtmg_sw_spcvmc.o parrrsw.o \
	rrsw_con.o rrtmg_sw_setcoef.o mcica_subcol_gen_sw.o rrtmg_sw_cldprmc.o \
	parkind.o rrsw_aer.o rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_reftra.o : $(RRTMG_SW_SRC)/rrtmg_sw_reftra.f90 rrsw_vsn.o parkind.o \
	rrsw_tbl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_setcoef.o : $(RRTMG_SW_SRC)/rrtmg_sw_setcoef.f90 rrsw_vsn.o parrrsw.o \
	rrsw_ref.o parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvmc.f90 rrsw_vsn.o parrrsw.o \
	rrsw_tbl.o parkind.o rrtmg_sw_reftra.o rrtmg_sw_vrtqdr.o rrtmg_sw_taumol.o \
	rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvrt.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvrt.f90 rrsw_vsn.o parrrsw.o \
	rrsw_tbl.o parkind.o rrtmg_sw_reftra.o rrtmg_sw_vrtqdr.o rrtmg_sw_taumol.o \
	rrsw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_taumol.o : $(RRTMG_SW_SRC)/rrtmg_sw_taumol.f90 rrsw_kg21.o rrsw_kg19.o \
	rrsw_kg23.o rrsw_kg26.o rrsw_vsn.o parrrsw.o rrsw_con.o rrsw_wvn.o rrsw_kg29.o \
	rrsw_kg22.o rrsw_kg20.o rrsw_kg18.o rrsw_kg17.o rrsw_kg28.o rrsw_kg24.o \
	rrsw_kg27.o parkind.o rrsw_kg16.o rrsw_kg25.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_vrtqdr.o : $(RRTMG_SW_SRC)/rrtmg_sw_vrtqdr.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRUser.o : $(SURFACE)/ModRUser.f90 mem_grid.o mem_leaf.o ccatt_start.o \
	ModLeafComs.o rconstants.o node_mod.o memSoilMoisture.o io_params.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

setup.o : $(MATRIX)/setup.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

sfclyr_jules.o : $(JULES_DIR)/sfclyr_jules.f90 sf_diags_mod.o io_params.o \
	rconstants.o ModMicControl.o gridmean_fluxes.o ModBasicFields.o \
	model_time_mod.o datetime_mod.o io_constants.o chem1_list.o mem_grid.o \
	mem_leaf.o ModLeaf3Init.o ModTurbFields.o node_mod.o gridbox_mean_mod.o \
	jules_surface_types_mod.o csigma_mod.o ModMicroFields.o fluxes.o mem_radiate.o \
	mem_chem1.o ModLeafComs.o mem_cuparm.o ancil_info.o JulesFields.o \
	mem_brams_jules.o jules_fields_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

shcu_vars_const.o : $(CUPARM)/shcu_vars_const.f90 ModConvComs.o grid_dims.o \
	ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSoilMoisture.o : $(SOIL_MOISTURE)/ModSoilMoisture.F90 mem_aerad.o mem_grid.o \
	parlibf.o mem_leaf.o ModBasicFields.o ModNamelistFile.o rconstants.o \
	ModLeafComs.o ModTurbFields.o ReadBcst.o node_mod.o ModControlVars.o \
	io_params.o ModMPassFull.o memSoilMoisture.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
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

ModTKenn.o : $(STILT)/ModTKenn.f90 mem_stilt.o turb_constants.o mem_grid.o \
	rconstants.o mem_scratch.o 
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

ModVarTables.o : $(MEMORY)/ModVarTables.f90 chem1_list.o \
	ModParallelEnvironment.o aer1_list.o VarTable.o io_params.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

VarTable.o : $(MEMORY)/VarTable.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfUpdate.o : $(FDDA)/ModVarfUpdate.f90 mem_grid.o ModInitHis.o rconstants.o \
	mem_scratch.o ref_sounding.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

include jules_depend_model.mk

