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

ModAsTi.o : $(ISAN)/ModAsTi.f90 mem_grid.o rconstants.o isan_coms.o \
	ModChemAObj.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsTp.o : $(ISAN)/ModAsTp.f90 rconstants.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAVarF.o : $(ISAN)/ModAVarF.f90 mem_grid.o rconstants.o isan_coms.o ModRbnd.o 
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

ModChemAsti.o : $(ISAN_CHEM)/ModChemAsti.f90 ModChemVInterps.o mem_chem1.o \
	ModChemAObj.o isan_coms.o mem_aer1.o chem_isan_coms.o ModChemAsti2.o ModAsTi.o \
	mem_grid.o ModChemFirstRams.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsti2.o : $(ISAN_CHEM)/ModChemAsti2.f90 ModDateUtils.o rconstants.o \
	isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAstp.o : $(ISAN_CHEM)/ModChemAstp.F90 mem_varinit.o mem_chem1.o dump.o \
	rconstants.o ModDateUtils.o isan_coms.o chem1_list.o chem_isan_coms.o ModAsTp.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAvarf.o : $(ISAN_CHEM)/ModChemAvarf.f90 mem_chem1.o ModRbnd.o \
	rconstants.o mem_aer1.o isan_coms.o ModNestFeed.o chem_isan_coms.o ModAVarF.o \
	mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_dry_dep.o : $(MODEL_CHEM)/chem_dry_dep.f90 mem_chem1.o ModDateUtils.o \
	mem_aer1.o extra.o chem1_list.o aer1_list.o 
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

chem_fastjx_driv.o : $(CCATT)/chem_fastjx_driv.f90 chem1_list.o rconstants.o \
	chem_fastjx_data.o chem_fastjx57.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFileInv.o : $(ISAN_CHEM)/ModChemFileInv.f90 ModDateUtils.o mem_grid.o \
	isan_coms.o dump.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemFirstRams.o : $(ISAN_CHEM)/ModChemFirstRams.f90 ModNestFillDens.o \
	an_header.o ModChemRefState.o rconstants.o isan_coms.o grid_dims.o \
	mem_scratch.o ModGetVar.o mem_grid.o ModRcio.o $(UTILS_INCS)/files.h \
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

chem_plumerise_scalar.o : $(CCATT)/chem_plumerise_scalar.f90 node_mod.o \
	mem_plume_chem1.o mem_chem1.o mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemRefState.o : $(ISAN_CHEM)/ModChemRefState.f90 ModNestFillDens.o \
	rconstants.o ccatt_start.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_sources.o : $(CCATT)/chem_sources.f90 ReadBcst.o ModControlVars.o \
	mem_chem1.o ModNamelistFile.o ModDateUtils.o mem_volc_chem1.o mem_aer1.o \
	parlibf.o io_params.o mem_plume_chem1.o mem_grid.o aer1_list.o \
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
	mem_chem1.o chem_spack_fexloss.o chem_spack_kinetic.o chem_spack_rates.o \
	chem_spack_fexprod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rates.o : $(MODEL_CHEM)/chem_spack_rates.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rodas3_dyndt.o : $(CCATT)/chem_spack_rodas3_dyndt.f90 mem_chem1.o \
	mem_spack.o chem_spack_fexchem.o extra.o chem_spack_kinetic.o chem_spack_ros.o \
	mem_grid.o chem_spack_jacdchemdc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros.o : $(CCATT)/chem_spack_ros.f90 chem_spack_solve_sparse.o \
	mem_chem1.o mem_spack.o chem_spack_fexchem.o chem_spack_kinetic.o \
	chem_spack_jacdchemdc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros_dyndt.o : $(CCATT)/chem_spack_ros_dyndt.f90 \
	chem_spack_solve_sparse.o mem_chem1.o mem_spack.o chem_spack_fexchem.o \
	chem_spack_kinetic.o chem_spack_ros.o chem_spack_jacdchemdc.o 
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

ModChemVInterps.o : $(ISAN_CHEM)/ModChemVInterps.f90 rconstants.o isan_coms.o \
	dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemDryDepDriver.o : $(MODEL_CHEM)/ChemDryDepDriver.f90 chem_dry_dep.o \
	mem_chem1.o ModNamelistFile.o ModMicroFields.o rconstants.o mem_aer1.o \
	mem_leaf.o CuParmFields.o grid_dims.o RadiateFields.o ModBasicFields.o \
	mem_grid.o ModTurbFields.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemSourcesDriver.o : $(CCATT)/ChemSourcesDriver.f90 mem_chem1.o \
	chem_plumerise_scalar.o mem_stilt.o mem_volc_chem1.o mem_aer1.o mem_leaf.o \
	chem1_list.o io_params.o mem_plume_chem1.o chem_sources.o aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

coag.o : $(MATRIX)/coag.f90 setup.o memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModCondRead.o : $(FDDA)/ModCondRead.f90 mem_varinit.o ModRamsGrid.o \
	ModDateUtils.o isan_coms.o ModNudAnalysis.o VarTable.o mem_grid.o \
	ModCondUpdate.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCondUpdate.o : $(FDDA)/ModCondUpdate.f90 grid_struct.o mem_varinit.o \
	an_header.o ModInitHis.o mem_grid.o VarTable.o ModRcio.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
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

ModCuRead.o : $(CUPARM)/ModCuRead.f90 ModNamelistFile.o ModDateUtils.o \
	isan_coms.o CuParmFields.o grid_dims.o CuParmVars.o mem_grid.o \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
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

ModCupGrellCattDeep.o : $(CUPARM)/ModCupGrellCattDeep.f90 ModCupUp.o \
	mem_grell_param2.o ModCupEnv.o mem_varinit.o mem_carma.o ModCupEnvCatt.o \
	Phys_const.o node_mod.o cup_output_vars.o mem_scratch3_grell.o ModCupDn.o \
	kbcon_ecmwf.o mem_grid.o ccatt_start.o mem_scratch2_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCupGrellCattShallow.o : $(CUPARM)/ModCupGrellCattShallow.f90 ModCupUp.o \
	mem_grell_param2.o ModCupEnv.o mem_varinit.o mem_scratch3_grell_sh.o \
	ModCupEnvCatt.o Phys_const.o node_mod.o cup_output_vars.o \
	mem_scratch2_grell_sh.o mem_grid.o 
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

dam.o : $(ENERGY)/dam.f90 ModDateUtils.o ModNamelistFile.o mem_grid.o dump.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

depv.o : $(MATRIX)/depv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

digitalFilter.o : $(MODEL)/digitalFilter.f90 ReadBcst.o ModControlVars.o \
	ModNamelistFile.o utilsMod.o ModDateUtils.o node_mod.o grid_dims.o io_params.o \
	ModBasicFields.o mem_grid.o VarTable.o $(UTILS_INCS)/files.h 
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

ModGridSet.o : $(INIT)/ModGridSet.f90 grid_dims.o mem_grid.o rconstants.o 
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

ModLeaf3Init.o : $(SURFACE)/ModLeaf3Init.f90 rconstants.o mem_leaf.o \
	ModLeafComs.o io_params.o grid_dims.o teb_spm_start.o mem_grid.o ModLeaf3.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Teb.o : $(SURFACE)/ModLeaf3Teb.f90 ModGasPart.o mem_teb_vars_const.o \
	mem_emiss.o ModUrban.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeafComs.o : $(SURFACE)/ModLeafComs.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

local_proc.o : $(MODEL)/local_proc.F90 ReadBcst.o dump.o mem_stilt.o \
	rconstants.o node_mod.o grid_dims.o io_params.o ref_sounding.o mem_grid.o \
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
	mcica_random_numbers.o rrlw_con.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_sw.o : $(RRTMG_SW_SRC)/mcica_subcol_gen_sw.f90 \
	mcica_random_numbers.o rrsw_wvn.o rrsw_con.o parrrsw.o parkind.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aer1.o : $(CCATT)/mem_aer1.f90 mem_chem1.o ModNamelistFile.o dump.o \
	ModScalarTable.o node_mod.o chem1_list.o grid_dims.o io_params.o VarTable.o \
	mem_grid.o aer1_list.o $(UTILS_INCS)/constants.h 
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

mem_carma.o : $(RADIATE)/mem_carma.f90 ReadBcst.o ModControlVars.o mem_aerad.o \
	ModNamelistFile.o ModRamsGrid.o ModSoilMoisture.o mem_globrad.o node_mod.o \
	grid_dims.o io_params.o parlibf.o ModMPassFull.o ModBasicFields.o mem_grid.o \
	VarTable.o ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1.o : $(CCATT)/mem_chem1.f90 ModNamelistFile.o ModScalarTable.o \
	chem1_list.o grid_dims.o io_params.o VarTable.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1aq.o : $(CCATT)/mem_chem1aq.f90 mem_chem1.o ModNamelistFile.o \
	chem1aq_list.o ModScalarTable.o grid_dims.o VarTable.o \
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

CuParmFields.o : $(CUPARM)/CuParmFields.f90 ModNodeDimensions.o \
	ModParallelEnvironment.o VarTable.o ModNamelistFile.o $(UTILS_INCS)/files.h 
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

GaspartFields.o : $(TEB_SPM)/GaspartFields.f90 ModNamelistFile.o \
	ModNodeDimensions.o ModParallelEnvironment.o VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globaer.o : $(RADIATE)/mem_globaer.f90 mem_precision.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globrad.o : $(RADIATE)/mem_globrad.f90 parlibf.o mem_aerad.o mem_precision.o \
	ModNamelistFile.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
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

mem_grid.o : $(MEMORY)/mem_grid.f90 grid_dims.o VarTable.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid_dim_defs.o : $(MEMORY)/mem_grid_dim_defs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

JulesFields.o : $(SURFACE)/JulesFields.f90 ModNamelistFile.o ModNodeDimensions.o \
	ModParallelEnvironment.o VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_leaf.o : $(SURFACE)/mem_leaf.f90 ModNamelistFile.o grid_dims.o io_params.o \
	teb_spm_start.o VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicroFields.o : $(MICRO)/ModMicroFields.f90 ModNamelistFile.o \
	ModNodeDimensions.o ModParallelEnvironment.o VarTable.o ModMicControl.o 
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

mem_oda.o : $(FDDA)/mem_oda.f90 grid_dims.o VarTable.o ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_opt_scratch.o : $(TURB)/mem_opt_scratch.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_plume_chem1.o : $(CCATT)/mem_plume_chem1.f90 ModNamelistFile.o chem1_list.o \
	VarTable.o mem_chem1.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_precision.o : $(RADIATE)/mem_precision.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

RadiateFields.o : $(RADIATE)/RadiateFields.f90 ModNodeDimensions.o \
	ModParallelEnvironment.o VarTable.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_rrtm.o : $(RADIATE)/mem_rrtm.f90 rrtmg_sw_init.o mem_chem1.o parrrtm.o \
	rconstants.o node_mod.o chem1_list.o parrrsw.o parkind.o rrtmg_lw_init.o \
	mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ScalarFields.o : $(MEMORY)/ScalarFields.f90 ModNodeDimensions.o \
	ModParallelEnvironment.o VarTable.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch.o : $(MEMORY)/mem_scratch.f90 mem_aerad.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch1_brams.o : $(MEMORY)/mem_scratch1_brams.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scratch1_grell.o : $(CUPARM)/mem_scratch1_grell.f90 mem_grell_param2.o \
	ccatt_start.o dump.o $(UTILS_INCS)/constants.h 
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

ShcuFields.o : $(CUPARM)/ShcuFields.f90 ModControlVars.o ModNamelistFile.o \
	ModNodeDimensions.o ModParallelEnvironment.o VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_spack.o : $(CCATT)/mem_spack.f90 chem1_list.o chem_spack_utils.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_stilt.o : $(STILT)/mem_stilt.f90 ModNamelistFile.o rconstants.o grid_dims.o \
	io_params.o VarTable.o $(UTILS_INCS)/constants.h 
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

mem_tend.o : $(MEMORY)/mem_tend.f90 ScalarFields.o ModNamelistFile.o \
	ModScalarTable.o ModMicroFields.o GaspartFields.o ModBasicFields.o \
	teb_spm_start.o mem_grid.o mem_emiss.o ModTurbFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_turb_scalar.o : $(TURB)/mem_turb_scalar.f90 VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tuv.o : $(TUV)/mem_tuv.f90 ModTuv2.7.o mem_globrad.o VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mem_varinit.o : $(MEMORY)/mem_varinit.f90 mem_chem1.o ModNamelistFile.o \
	chem1_list.o grid_dims.o VarTable.o $(UTILS_INCS)/files.h 
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

meteogram.o : $(IO)/meteogram.f90 ModNamelistFile.o meteogramType.o \
	ModPostUtils.o node_mod.o ModMPassDtl.o satPolyColision.o mem_grid.o VarTable.o \
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

ModMicInit.o : $(MICRO)/ModMicInit.f90 dump.o rconstants.o ModMicGamma.o \
	mem_grid.o ModMicTabs.o ModMicControl.o $(MICRO)/MicConstants.h \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
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

ModMicControl.o : $(MICRO)/ModMicControl.f90 grid_dims.o ModNamelistFile.o \
	ModParallelEnvironment.o $(MICRO)/MicConstants.h $(UTILS_INCS)/files.h 
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

ModAcoust.o : $(MODEL)/ModAcoust.f90 rconstants.o ModParallelEnvironment.o \
	node_mod.o mem_scratch.o mem_tend.o ModMessageSet.o ref_sounding.o \
	ModBasicFields.o ModAcoustAdap.o mem_grid.o ModGrid.o ModMicControl.o \
	$(UTILS_INCS)/tsNames.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAerClim.o : $(AERCLIM)/ModAerClim.f90 ReadBcst.o ModControlVars.o dump.o \
	ModSoilMoisture.o node_mod.o parlibf.o ModBasicFields.o mem_grid.o \
	ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBasicFields.o : $(MEMORY)/ModBasicFields.f90 ModNamelistFile.o \
	ModNodeDimensions.o mem_stilt.o ModParallelEnvironment.o VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBramsGrid.o : $(POST_SRC)/ModBramsGrid.f90 mem_aerad.o ModNamelistFile.o \
	ModPostUtils.o ModParallelEnvironment.o node_mod.o ref_sounding.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBuffering.o : $(MPI)/ModBuffering.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCarmaDriver.o : $(RADIATE)/ModCarmaDriver.f90 rad_carma.o ModNamelistFile.o \
	mem_carma.o ModDateUtils.o ModMicroFields.o rconstants.o mem_leaf.o node_mod.o \
	CuParmFields.o RadiateFields.o mem_tend.o grid_dims.o mem_scratch1_grell.o \
	ModBasicFields.o teb_spm_start.o mem_grid.o mem_teb_common.o ModLeaf3.o \
	ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsgen.o : $(ISAN_CHEM)/ModChemAsgen.F90 dump.o ModMkSfcTop.o \
	ModRamsGrid.o ModDateUtils.o node_mod.o chem1_list.o grid_dims.o io_params.o \
	ModChemAstp.o mem_aer1.o mem_grid.o ModChemFileInv.o isan_coms.o ModAsGen.o \
	ModChemIsanIo.o chem_isan_coms.o ModChemAsti.o ModControlVars.o mem_chem1.o \
	ModChemRefState.o ModChemAvarf.o aer1_list.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemConvTransp.o : $(CCATT)/ModChemConvTransp.f90 mem_grell_param2.o \
	mem_chem1.o mem_aer1.o Phys_const.o node_mod.o chem1_list.o mem_scratch.o \
	mem_scratch1_grell.o mem_tconv.o mem_grid.o aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemistryDriver.o : $(CCATT)/ModChemistryDriver.f90 mem_grell_param2.o \
	chem_spack_solve_sparse.o chem_orage.o mem_stilt.o rconstants.o node_mod.o \
	chem1_list.o grid_dims.o chem_spack_rodas3_dyndt.o mem_scratch1_grell.o \
	mem_chemic.o chem_fastjx_driv.o mem_aerad.o chem1aq_list.o mem_chem1aq.o \
	ModMicroFields.o mem_spack.o mem_aer1.o chem_uv_att.o carma_fastjx.o mem_grid.o \
	chem_spack_ros_dyndt.o chem_trans_liq.o ModNamelistFile.o parrrtm.o \
	mem_globrad.o extra.o RadiateFields.o ModBasicFields.o chem_spack_utils.o \
	mem_rrtm.o ModTuvDriver2.7.o mem_chem1.o mem_carma.o chem_spack_qssa.o \
	chem_trans_gasaq.o chem_spack_ros.o aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModControlVars.o : $(INIT)/ModControlVars.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCoriolis.o : $(MODEL)/ModCoriolis.f90 rconstants.o node_mod.o mem_scratch.o \
	parlibf.o mem_tend.o ref_sounding.o ModBasicFields.o mem_grid.o ModBuffering.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuParGrell3.o : $(CUPARM)/ModCuParGrell3.F90 mem_grell_param2.o \
	module_cu_g3.o mem_stilt.o rconstants.o node_mod.o Phys_const.o grid_dims.o \
	io_params.o mem_scratch1_grell.o module_cu_gf.o ModGrid.o mem_varinit.o \
	ModMicroFields.o mem_grid.o ModNamelistFile.o CuParmFields.o RadiateFields.o \
	mem_tend.o ModMessageSet.o mem_grell.o ModBasicFields.o ccatt_start.o \
	ModMicControl.o ConvPar_GF_GEOS5.o mem_carma.o mem_chem1.o module_cu_gf_v5.1.o \
	mem_leaf.o ModChemConvTransp.o mem_scratch.o ModRstilt.o ModRadvc.o VarTable.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDateUtils.o : $(UTILS_MODS)/ModDateUtils.f90 dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffSclr.o : $(TURB)/ModDiffSclr.f90 mem_scratch.o mem_grid.o ModTurbDiff.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffuse.o : $(TURB)/ModDiffuse.f90 ModNamelistFile.o ModScalarTable.o \
	ModTurbKE.o ModMicroFields.o ModDiffSclr.o mem_leaf.o node_mod.o \
	mem_opt_scratch.o mem_scratch.o mem_tend.o ModTurbDiff.o ke_coms.o ModTurbK.o \
	ModBasicFields.o mem_grid.o ModTurbFields.o ModMicControl.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDomainDecomp.o : $(MPI)/ModDomainDecomp.f90 ModParallelEnvironment.o \
	ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModEvaluation.o : $(EVAL)/ModEvaluation.f90 parlibf.o mem_grid.o \
	ModNamelistFile.o node_mod.o 
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

ModGasPart.o : $(TEB_SPM)/ModGasPart.f90 mem_teb_vars_const.o an_header.o \
	ModRcio.o mem_leaf.o node_mod.o grid_dims.o parlibf.o GaspartFields.o \
	ModBasicFields.o mem_grid.o VarTable.o mem_emiss.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGeodat.o : $(MKSFC)/ModGeodat.f90 teb_spm_start.o io_params.o mem_grid.o \
	mem_leaf.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGrid.o : $(MPI)/ModGrid.F90 Aero2McphysFields.o CuParmVars.o GaspartFields.o \
	ModGridDims.o ShcuFields.o ModMicroFields.o ModParallelEnvironment.o \
	ModTurbFields.o ModNamelistFile.o ModDomainDecomp.o meteogramType.o \
	CuParmFields.o RadiateFields.o mem_tend.o ModMessageSet.o ModBasicFields.o \
	ModMicControl.o ModControlVars.o ScalarFields.o ModScalarTable.o \
	ModNodeDimensions.o ModNeighbourNodes.o JulesFields.o VarTable.o aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridDims.o : $(MPI)/ModGridDims.f90 ModNamelistFile.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGridTree.o : $(MPI)/ModGridTree.f90 ModNamelistFile.o ModGrid.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modIau.o : $(MODEL)/modIau.f90 ReadBcst.o mem_varinit.o ModNamelistFile.o dump.o \
	node_mod.o parlibf.o mem_tend.o ModMPassFull.o mem_grid.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitHis.o : $(IO)/ModInitHis.f90 ModRamsGrid.o rconstants.o chem1_list.o \
	io_params.o ModRamsReadHeader.o ModLeaf3.o mem_aerad.o mem_varinit.o \
	ModLeafComs.o ModGetVar.o mem_grid.o ModRcio.o ModRinit.o ModBasicFields.o \
	ModMicControl.o mem_chem1.o an_header.o mem_leaf.o mem_scratch.o ref_sounding.o \
	VarTable.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitMicThompson.o : $(MICRO)/ModInitMicThompson.f90 ReadBcst.o generic.o \
	dump.o ModDateUtils.o ModMicroFields.o node_mod.o parlibf.o ModBasicFields.o \
	mem_grid.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLanduseInput.o : $(MKSFC)/ModLanduseInput.f90 ModLeaf3Init.o mem_mksfc.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3.o : $(SURFACE)/ModLeaf3.f90 ModNamelistFile.o ModLeaf3Teb.o \
	ModMicroFields.o rconstants.o node_mod.o CuParmFields.o mem_leaf.o \
	ModLeafComs.o io_params.o RadiateFields.o ModLeaf3Hyd.o ModBasicFields.o \
	mem_teb.o teb_spm_start.o ModMicControl.o mem_grid.o ccatt_start.o \
	ModTurbFields.o mem_teb_common.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3OceanOnly.o : $(SURFACE)/ModLeaf3OceanOnly.f90 ConvPar_GF_GEOS5.o \
	ModNamelistFile.o rconstants.o mem_leaf.o CuParmFields.o node_mod.o \
	RadiateFields.o io_params.o ModLeafComs.o ModTurbFields.o ModBasicFields.o \
	mem_grid.o ccatt_start.o ModLeaf3.o ModCuParGrell3.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMatrixDriver.o : $(MATRIX)/ModMatrixDriver.F90 Aero2McphysFields.o \
	mem_chem1.o ModParticle.o isrpia.o ModMicroFields.o rconstants.o mem_aer1.o \
	memMatrix.o mem_leaf.o node_mod.o chem1_list.o coag.o npf.o setup.o subs.o \
	ModBasicFields.o mem_grid.o aer1_list.o ModTurbFields.o ModMicControl.o 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ModMemAlloc.o : $(MEMORY)/ModMemAlloc.F90 Aero2McphysFields.o mem_grell_param2.o \
	mem_globaer.o mem_tuv.o mem_stilt.o node_mod.o chem1_list.o mem_oda.o \
	grid_dims.o io_params.o GaspartFields.o mem_scratch1_grell.o \
	mem_scratch2_grell_sh.o mem_turb_scalar.o mem_plume_chem1.o mem_nestb.o \
	chem_sources.o ModGrid.o mem_emiss.o mem_chemic.o mem_teb_common.o \
	mem_varinit.o mem_aerad.o chem1aq_list.o digitalFilter.o mem_chem1aq.o \
	ShcuFields.o ModEvaluation.o ModMicroFields.o mem_scratch1_brams.o \
	ModParallelEnvironment.o mem_aer1.o mem_opt_scratch.o ModLeafComs.o \
	carma_fastjx.o ModOptical.o mem_grid.o VarTable.o ModTurbFields.o \
	mem_teb_vars_const.o machine_arq.o mem_globrad.o mem_volc_chem1.o \
	CuParmFields.o extra.o RadiateFields.o mem_tend.o mem_grell.o ModBasicFields.o \
	modIau.o teb_spm_start.o shcu_vars_const.o ccatt_start.o mem_scratch3_grell.o \
	mem_grid_dim_defs.o chem_dry_dep.o ScalarFields.o mem_carma.o mem_chem1.o \
	mem_scratch3_grell_sh.o mem_leaf.o mem_scratch.o parrrsw.o JulesFields.o \
	mem_teb.o ModTuv2.7.o aer1_list.o ModCuParGrell3.o mem_scratch2_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMemory.o : $(UTILS_LIB)/ModMemory.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageData.o : $(MPI)/ModMessageData.f90 ModFieldSection.o \
	ModParallelEnvironment.o ModFieldSectionList.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMessageSet.o : $(MPI)/ModMessageSet.f90 ModNamelistFile.o ModMessageData.o \
	ModDomainDecomp.o ModNodeDimensions.o ModFieldSection.o \
	ModParallelEnvironment.o ModNeighbourNodes.o parlibf.o ModGridDims.o \
	ModFieldSectionList.o mem_grid.o VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicGfdlDriver.o : $(MICRO)/ModMicGfdlDriver.f90 ModNamelistFile.o \
	ModMicroFields.o rconstants.o mem_leaf.o node_mod.o io_params.o \
	ModBasicFields.o mem_grid.o gfdl_cloud_microphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsDrive.o : $(MICRO)/ModMicrophysicsDrive.f90 ModMicInit.o \
	mem_chem1.o mem_chem1aq.o ModMicroFields.o ModMicrophysicsMisc.o node_mod.o \
	ModMicVap.o grid_dims.o ModMicColl.o mem_grid.o ModBasicFields.o ModMicNuc.o \
	ModMicTabs.o mem_chemic.o ModMicControl.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsMisc.o : $(MICRO)/ModMicrophysicsMisc.f90 ModMicroFields.o \
	rconstants.o ModBasicFields.o mem_grid.o ModMicControl.o \
	$(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicThompsonDriver.o : $(MICRO)/ModMicThompsonDriver.f90 ModNamelistFile.o \
	ModMicroFields.o rconstants.o mem_leaf.o node_mod.o io_params.o \
	ModBasicFields.o module_mp_thompson.o mem_grid.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcDriver.o : $(MKSFC)/ModMkSfcDriver.f90 ModSstRead.o ReadBcst.o \
	ModNdviRead.o ModControlVars.o ModNestGeoSst.o ModMkSfcSfc.o ModMkSfcFuso.o \
	ModMkSfcTop.o mem_mksfc.o node_mod.o grid_dims.o io_params.o ModLanduseInput.o \
	teb_spm_start.o mem_grid.o ModMkSfcNdvi.o ModMkSfcSst.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcFuso.o : $(MKSFC)/ModMkSfcFuso.f90 ReadBcst.o ModControlVars.o \
	mem_teb_vars_const.o mem_mksfc.o node_mod.o io_params.o GaspartFields.o \
	mem_teb.o mem_grid.o mem_emiss.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcNdvi.o : $(MKSFC)/ModMkSfcNdvi.f90 dump.o ModRUser.o mem_mksfc.o \
	mem_leaf.o io_params.o ModLanduseInput.o mem_grid.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSfc.o : $(MKSFC)/ModMkSfcSfc.f90 ReadBcst.o ModControlVars.o dump.o \
	mem_mksfc.o mem_leaf.o node_mod.o io_params.o mem_grid.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSst.o : $(MKSFC)/ModMkSfcSst.f90 ModNestFillDens.o ModRUser.o \
	mem_mksfc.o mem_leaf.o grid_dims.o io_params.o ModGeodat.o mem_grid.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcTop.o : $(MKSFC)/ModMkSfcTop.f90 ReadBcst.o ModControlVars.o dump.o \
	mem_mksfc.o node_mod.o io_params.o mem_grid.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMonotonicAdvection.o : $(MODEL)/ModMonotonicAdvection.f90 chem_dry_dep.o \
	mem_chem1.o ModNamelistFile.o ModDomainDecomp.o rconstants.o \
	ModParallelEnvironment.o mem_aer1.o ModMessageSet.o ModGrid.o mem_grid.o \
	ccatt_start.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNamelistFile.o : $(INIT)/ModNamelistFile.f90 dump.o ModParallelEnvironment.o \
	modPrintInitial.o grid_dims.o parlibf.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNdviRead.o : $(MKSFC)/ModNdviRead.f90 ReadBcst.o ModControlVars.o mem_mksfc.o \
	ModDateUtils.o mem_leaf.o node_mod.o grid_dims.o io_params.o mem_grid.o \
	ModMkSfcNdvi.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNeighbourNodes.o : $(MPI)/ModNeighbourNodes.f90 ModParallelEnvironment.o \
	ModGridDims.o ModDomainDecomp.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestGeoSst.o : $(MKSFC)/ModNestGeoSst.f90 dump.o ModMkSfcTop.o node_mod.o \
	ModLeaf3Init.o grid_dims.o io_params.o ModLanduseInput.o mem_grid.o \
	ModTurbFields.o ModRUser.o mem_mksfc.o memSoilMoisture.o ModNestFeed.o \
	ModBasicFields.o ccatt_start.o ModNestFillDens.o ModControlVars.o \
	ModSoilMoisture.o mem_leaf.o mem_scratch.o ModInitHis.o ModGeodat.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNodeDimensions.o : $(MPI)/ModNodeDimensions.f90 ModParallelEnvironment.o \
	ModGridDims.o ModDomainDecomp.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudAnalysis.o : $(FDDA)/ModNudAnalysis.f90 ModNestFillDens.o mem_varinit.o \
	mem_chem1.o dump.o ModEvaluation.o node_mod.o chem1_list.o mem_scratch.o \
	mem_tend.o ModBasicFields.o mem_grid.o modIau.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaNudge.o : $(FDDA)/ModOdaNudge.f90 ModOdaProcObs.o node_mod.o mem_oda.o \
	mem_tend.o mem_scratch.o ModOdaKrig.o io_params.o ModBasicFields.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOneProc.o : $(MODEL)/ModOneProc.F90 ModRThrm.o mem_stilt.o ModGasPart.o \
	mem_oda.o CuParmVars.o io_params.o mem_plume_chem1.o ModMkSfcDriver.o \
	ModPostGridNetCDF.o mem_teb_common.o ModNestIntrp.o ModRanlavg.o ModRnode.o \
	ModRUser.o ModRinit.o extra.o parlibf.o ModOpspec.o chem_dry_dep.o \
	ModTuvDriver2.7.o mem_carma.o dump.o ModRamsGrid.o chem1_list.o tuvParameter.o \
	chem_sources.o ModNestGeoSst.o mem_varinit.o ModCuRead.o ModChemistryDriver.o \
	mem_aer1.o dam.o ModRecycle.o ModMicInit.o mem_teb_vars_const.o ModMkSfcSfc.o \
	machine_arq.o mem_globrad.o mem_volc_chem1.o ModUrbanCanopy.o ModMPassDtl.o \
	ModBasicFields.o ModSstRead.o mem_chem1.o ModTimestep.o mem_leaf.o ModRhhi.o \
	ModGridTree.o ref_sounding.o mem_grell_param2.o ModMkSfcTop.o ModMkSfcFuso.o \
	ModLeaf3Init.o ModTimestepRK.o mem_emiss.o ReadBcst.o mem_chem1aq.o \
	ModEvaluation.o ModParallelEnvironment.o ModTimeStamp.o ModNamelistFile.o \
	ModDomainDecomp.o ModNudRead.o ModMicrophysicsMisc.o memSoilMoisture.o \
	isan_coms.o ModRio.o teb_spm_start.o ModCondRead.o mem_scratch.o \
	ModVarfUpdate.o mem_teb.o aer1_list.o modIau.o VarTable.o \
	ModRamsMicrophysics2M.o ModWindFarm.o local_proc.o node_mod.o grid_dims.o \
	ModGrid.o digitalFilter.o ModOdaRead.o ModPostProcess.o mem_grid.o ModSched.o \
	ModMemAlloc.o ModChemAsgen.o ModMonotonicAdvection.o meteogram.o ModVarfFile.o \
	ModAerClim.o ModCoriolis.o ccatt_start.o ModInitMicThompson.o ModParaInit.o \
	ModNdviRead.o ModSoilMoisture.o ModLeaf3Teb.o ModInitHis.o ModTuv2.7.o \
	shcu_vars_const.o ModCuParGrell3.o domain_decomp.o $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/tsNames.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOpspec.o : $(IO)/ModOpspec.f90 mem_grell_param2.o mem_stilt.o chem1_list.o \
	grid_dims.o io_params.o chem_sources.o mem_emiss.o mem_varinit.o chem1aq_list.o \
	mem_chem1aq.o mem_aer1.o mem_grid.o ModNamelistFile.o mem_globrad.o aer1_list.o \
	teb_spm_start.o ccatt_start.o shcu_vars_const.o ModMicControl.o mem_chem1.o \
	mem_leaf.o modIau.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOptical.o : $(RADIATE)/ModOptical.f90 ReadBcst.o ModControlVars.o \
	ModNamelistFile.o dump.o ModRamsGrid.o ModSoilMoisture.o mem_aer1.o mem_leaf.o \
	node_mod.o parlibf.o ModMPassFull.o ModBasicFields.o VarTable.o mem_grid.o \
	ccatt_start.o ModTurbFields.o aer1_list.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOutputUtils.o : $(IO)/ModOutputUtils.f90 ModNamelistFile.o dump.o \
	ModBasicFields.o VarTable.o ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOzone.o : $(TEB_SPM)/ModOzone.f90 rconstants.o RadiateFields.o \
	GaspartFields.o ModBasicFields.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParaInit.o : $(MPI)/ModParaInit.f90 dump.o ModScalarTable.o node_mod.o \
	grid_dims.o mem_grid.o VarTable.o $(UTILS_INCS)/constants.h 
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

ModPostGrid.o : $(POST_SRC)/ModPostGrid.F90 ModOutputUtils.o ModNamelistFile.o \
	ModBramsGrid.o ModPostTypes.o ModPostUtils.o ModParallelEnvironment.o parlibf.o \
	io_params.o ModBasicFields.o mem_grid.o VarTable.o ModPostOneFieldNetCDF.o \
	ModTurbFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostGridNetCDF.o : $(POST_SRC)/ModPostGridNetCDF.F90 ModNamelistFile.o dump.o \
	ModPostTypes.o ModDateUtils.o ModPostUtils.o io_params.o mem_grid.o \
	ModBramsGrid.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField.o : $(POST_SRC)/ModPostOneField.f90 ModPostOneField2d.o \
	ModNamelistFile.o dump.o ModPostOneField7d.o ModPostTypes.o ModPostOneField3d.o \
	ModPostUtils.o ModPostOneField8d.o node_mod.o ModPostOneFieldUtils.o \
	ModBasicFields.o VarTable.o ModBramsGrid.o ModTurbFields.o ModMicControl.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField2d.o : $(POST_SRC)/ModPostOneField2d.f90 mem_aerad.o \
	ModNamelistFile.o dump.o ModPostTypes.o ModPostUtils.o node_mod.o \
	CuParmFields.o io_params.o ModPostGrid.o mem_grid.o ModPostOneFieldUtils.o \
	ModBasicFields.o ModOutputUtils.o VarTable.o ModBramsGrid.o ModTurbFields.o \
	ModMicControl.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField3d.o : $(POST_SRC)/ModPostOneField3d.f90 mem_varinit.o \
	ModNamelistFile.o ModPostTypes.o ModPostUtils.o node_mod.o ModPostGrid.o \
	mem_grid.o ModPostOneFieldUtils.o ModBasicFields.o ModOutputUtils.o VarTable.o \
	ModBramsGrid.o ModTurbFields.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField7d.o : $(POST_SRC)/ModPostOneField7d.f90 ModNamelistFile.o \
	ModPostTypes.o ModPostOneFieldUtils.o ModPostUtils.o ModBasicFields.o \
	ModOutputUtils.o VarTable.o ModBramsGrid.o ModTurbFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField8d.o : $(POST_SRC)/ModPostOneField8d.f90 ModNamelistFile.o \
	ModPostTypes.o ModPostOneFieldUtils.o ModPostUtils.o ModBasicFields.o \
	ModOutputUtils.o VarTable.o ModBramsGrid.o ModTurbFields.o 
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

ModPostProcess.o : $(POST_SRC)/ModPostProcess.F90 ModNamelistFile.o \
	ModPostOneField.o ModPostTypes.o ModParallelEnvironment.o ModGridTree.o \
	ModPostGrid.o ModMessageSet.o io_params.o ModBasicFields.o VarTable.o \
	ModPostGridNetCDF.o ModGrid.o ModBramsGrid.o ModTurbFields.o \
	$(UTILS_INCS)/tsNames.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostTypes.o : $(POST_SRC)/ModPostTypes.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostUtils.o : $(POST_SRC)/ModPostUtils.f90 ModParallelEnvironment.o dump.o \
	mem_leaf.o $(POST_INCS)/post_rconfig.h $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/constants.h $(POST_INCS)/post_rconstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modPrintInitial.o : $(INIT)/modPrintInitial.F90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvc.o : $(MODEL)/ModRadvc.f90 chem_dry_dep.o ModMonotonicAdvection.o \
	mem_chem1.o ModNamelistFile.o ModScalarTable.o ModParallelEnvironment.o \
	mem_aer1.o mem_scratch.o grid_dims.o mem_tend.o ModRadvcAdap.o ModBasicFields.o \
	mem_grid.o ccatt_start.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcRK.o : $(MODEL)/ModRadvcRK.f90 mem_chem1.o mem_stilt.o \
	ModParallelEnvironment.o node_mod.o grid_dims.o mem_tend.o ModMessageSet.o \
	mem_grid.o ModGrid.o ModRexev.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsMicrophysics2M.o : $(MICRO)/ModRamsMicrophysics2M.f90 dump.o \
	ModMicroFields.o rconstants.o mem_leaf.o node_mod.o grid_dims.o \
	ModBasicFields.o ModMicGamma.o mem_grid.o ModMicControl.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsReadHeader.o : $(IO)/ModRamsReadHeader.f90 an_header.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRanlavg.o : $(IO)/ModRanlavg.f90 ModRThrm.o ModMicroFields.o grid_dims.o \
	io_params.o ModBasicFields.o mem_grid.o VarTable.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRbnd.o : $(BC)/ModRbnd.f90 mem_chem1.o ModScalarTable.o ModTurbKE.o \
	ModMicroFields.o ModMicrophysicsMisc.o node_mod.o mem_scratch.o mem_tend.o \
	ref_sounding.o ModBasicFields.o mem_grid.o ccatt_start.o ModTurbFields.o \
	ModMicControl.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRcio.o : $(IO)/ModRcio.f90 ModNamelistFile.o an_header.o mem_stilt.o \
	ModParallelEnvironment.o mem_leaf.o grid_dims.o io_params.o ModLeafComs.o \
	ref_sounding.o mem_grid.o ModMicControl.o $(MICRO)/MicConstants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConv.o : $(CUPARM)/ModRConv.f90 ModNamelistFile.o rconstants.o node_mod.o \
	CuParmFields.o CuParmVars.o mem_tend.o ModConvComs.o ModBasicFields.o \
	mem_grid.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConvGrellCatt.o : $(CUPARM)/ModRConvGrellCatt.f90 mem_grell_param2.o \
	ModCupGrellCattShallow.o ModCupGrellCattDeep.o rconstants.o mem_stilt.o \
	node_mod.o mem_leaf.o ModChemConvTransp.o mem_scratch.o mem_tend.o ModRstilt.o \
	io_params.o mem_scratch1_grell.o mem_grell.o ModGrid.o mem_grid.o ccatt_start.o \
	ModCuParGrell3.o ModMicControl.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRecycle.o : $(IO)/ModRecycle.f90 ReadBcst.o mem_aerad.o mem_chem1.o \
	an_header.o dump.o ModDateUtils.o mem_aer1.o node_mod.o chem1_list.o \
	io_params.o ModGetVar.o ModMPassFull.o VarTable.o ModRamsReadHeader.o \
	mem_grid.o aer1_list.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRexev.o : $(STILT)/ModRexev.f90 mem_stilt.o rconstants.o mem_scratch.o \
	mem_tend.o ModBasicFields.o ModRadvc.o mem_grid.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRGrad.o : $(TURB)/ModRGrad.f90 mem_scratch.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRhhi.o : $(INIT)/ModRhhi.f90 ModRamsGrid.o ModRinit.o rconstants.o \
	grid_dims.o mem_scratch.o ref_sounding.o ModBasicFields.o mem_grid.o \
	ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRinit.o : $(INIT)/ModRinit.f90 mem_varinit.o ModRbnd.o ModTurbKE.o \
	ModMicroFields.o rconstants.o node_mod.o mem_scratch.o io_params.o \
	ref_sounding.o ModBasicFields.o mem_grid.o ModTurbFields.o ModMicControl.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRio.o : $(IO)/ModRio.f90 ModDateUtils.o mpi_io_engine-5d.o node_mod.o \
	grid_dims.o io_params.o ReadBcst.o mem_aerad.o ModParallelEnvironment.o \
	mem_grid.o ModRcio.o ModTurbFields.o ModNamelistFile.o utilsMod.o parlibf.o \
	ModMPassFull.o ModBasicFields.o ModMicControl.o ModControlVars.o mem_chem1.o \
	an_header.o ref_sounding.o VarTable.o $(UTILS_INCS)/interface.h \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRrtmDriver.o : $(RADIATE)/ModRrtmDriver.f90 mcica_subcol_gen_sw.o \
	mem_grell_param2.o mem_tuv.o rconstants.o mcica_subcol_gen_lw.o rrtmg_sw_rad.o \
	ModDateUtils.o node_mod.o grid_dims.o parkind.o parrrsw.o mem_scratch1_grell.o \
	ModMicroFields.o ModLeafComs.o rrtmg_lw_rad.o ModOptical.o mem_grid.o \
	ModNamelistFile.o parrrtm.o CuParmFields.o RadiateFields.o mem_tend.o \
	ModBasicFields.o teb_spm_start.o ccatt_start.o mem_rrtm.o ModMicControl.o \
	rrtmg_lw_cldprop.o mem_chem1.o mem_carma.o mem_leaf.o rrtmg_sw_cldprop.o \
	ref_sounding.o $(UTILS_INCS)/aerosol_setup.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRShCuPar.o : $(CUPARM)/ModRShCuPar.f90 ShcuFields.o ModMicroFields.o \
	node_mod.o ModRConv.o mem_tend.o ModConvComs.o ModBasicFields.o mem_grid.o \
	shcu_vars_const.o ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRstilt.o : $(STILT)/ModRstilt.f90 ModMonotonicAdvection.o mem_stilt.o \
	mem_scratch.o grid_dims.o mem_scratch1_grell.o ModBasicFields.o mem_grid.o \
	ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRThrm.o : $(MODEL)/ModRThrm.f90 ModMicroFields.o rconstants.o mem_scratch.o \
	ModBasicFields.o mem_grid.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRtimi.o : $(MODEL)/ModRtimi.f90 ModScalarTable.o node_mod.o mem_scratch.o \
	mem_tend.o mem_grell.o ModBasicFields.o mem_grid.o shcu_vars_const.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModScalarTable.o : $(MEMORY)/ModScalarTable.f90 ModParallelEnvironment.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSched.o : $(MODEL)/ModSched.f90 ReadBcst.o mem_varinit.o ModNamelistFile.o \
	dump.o local_proc.o ModParallelEnvironment.o isan_coms.o node_mod.o parlibf.o \
	io_params.o ref_sounding.o ModBasicFields.o mem_grid.o shcu_vars_const.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSeaSalt.o : $(CCATT)/ModSeaSalt.f90 mem_chem1.o mem_aer1.o mem_leaf.o \
	node_mod.o io_params.o ModAerClim.o ModBasicFields.o aer1_list.o mem_grid.o \
	ccatt_start.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSstRead.o : $(MKSFC)/ModSstRead.f90 ReadBcst.o ModControlVars.o mem_mksfc.o \
	ModDateUtils.o mem_leaf.o node_mod.o grid_dims.o io_params.o mem_grid.o \
	ModMkSfcSst.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimeStamp.o : $(MODEL)/ModTimeStamp.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestep.o : $(MODEL)/ModTimestep.F90 ChemDryDepDriver.o \
	ModRamsMicrophysics2M.o ModRShCuPar.o ModWindFarm.o ModRThrm.o mem_stilt.o \
	sfclyr_jules.o ModGasPart.o rconstants.o ChemSourcesDriver.o node_mod.o \
	grid_dims.o mem_oda.o mem_plume_chem1.o ModMatrixDriver.o ModAcoust.o \
	chem_sources.o ModGrid.o mem_emiss.o ModLeaf3.o ModMicrophysicsDrive.o \
	mem_varinit.o digitalFilter.o ModChemistryDriver.o ModRbnd.o ModOdaNudge.o \
	ModMicThompsonDriver.o mem_aer1.o ModTimeStamp.o ModRConv.o ModOptical.o \
	mem_grid.o rad_driv.o ModRConvGrellCatt.o ModMonotonicAdvection.o \
	ModMicGfdlDriver.o machine_arq.o ModOzone.o ModMicrophysicsMisc.o \
	ModUrbanCanopy.o mem_tend.o ModSeaSalt.o ModCoriolis.o ModMessageSet.o \
	ModRtimi.o ModTurbK.o ModBasicFields.o teb_spm_start.o ccatt_start.o \
	mem_chem1.o mem_leaf.o mem_scratch.o ModRstilt.o ModNudAnalysis.o ModDiffuse.o \
	ModCuParGrell3.o ModRadvc.o shcu_vars_const.o ModRexev.o \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestepRK.o : $(MODEL)/ModTimestepRK.F90 ChemDryDepDriver.o \
	ModRamsMicrophysics2M.o ModRShCuPar.o ModWindFarm.o ModRThrm.o mem_stilt.o \
	sfclyr_jules.o ModGasPart.o rconstants.o ChemSourcesDriver.o node_mod.o \
	grid_dims.o mem_oda.o mem_plume_chem1.o ModMatrixDriver.o ModAcoust.o \
	chem_sources.o ModGrid.o mem_emiss.o ModLeaf3.o ModMicrophysicsDrive.o \
	mem_varinit.o digitalFilter.o ModRbnd.o ModChemistryDriver.o ModOdaNudge.o \
	ModMicThompsonDriver.o ModParallelEnvironment.o mem_aer1.o ModTimeStamp.o \
	ModRConv.o ModOptical.o mem_grid.o rad_driv.o ModMonotonicAdvection.o \
	ModMicGfdlDriver.o machine_arq.o ModOzone.o utilsMod.o ModRadvcRK.o \
	ModMicrophysicsMisc.o ModLeaf3OceanOnly.o ModUrbanCanopy.o mem_tend.o \
	ModAerClim.o ModSeaSalt.o ModCoriolis.o ModMessageSet.o ModRtimi.o ModTurbK.o \
	teb_spm_start.o ccatt_start.o shcu_vars_const.o mem_chem1.o ModTimestep.o \
	mem_leaf.o mem_scratch.o ModRstilt.o ModNudAnalysis.o ModDiffuse.o \
	ModCuParGrell3.o ModRadvc.o modIau.o ModRexev.o $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiff.o : $(TURB)/ModTurbDiff.f90 ModRGrad.o mem_scratch.o mem_grid.o \
	mem_opt_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiffAdap.o : $(TURB)/ModTurbDiffAdap.f90 mem_scratch.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbFields.o : $(TURB)/ModTurbFields.f90 ModNamelistFile.o \
	ModNodeDimensions.o ModParallelEnvironment.o VarTable.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbK.o : $(TURB)/ModTurbK.f90 mem_stilt.o rconstants.o node_mod.o \
	grid_dims.o mem_turb_scalar.o ke_coms.o ModTKenn.o ModMicroFields.o \
	ModTurbDiffAdap.o ModTurbDiff.o mem_grid.o ModTurbFields.o \
	ModMonotonicAdvection.o ModNamelistFile.o ModTurbKE.o ModRGrad.o ModTurbKAdap.o \
	mem_tend.o mem_grell.o ModBasicFields.o ccatt_start.o ModMicControl.o \
	mem_chem1.o ModScalarTable.o mem_leaf.o mem_scratch.o ModRstilt.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKAdap.o : $(TURB)/ModTurbKAdap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKE.o : $(TURB)/ModTurbKE.f90 rconstants.o mem_scratch.o ke_coms.o \
	mem_grid.o ModTurbFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuv2.7.o : $(TUV)/ModTuv2.7.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuvDriver2.7.o : $(TUV)/ModTuvDriver2.7.f90 mem_aerad.o mem_chem1.o \
	ModNamelistFile.o mem_carma.o mem_tuv.o mem_globrad.o rconstants.o node_mod.o \
	mem_leaf.o chem1_list.o RadiateFields.o tuvParameter.o extra.o \
	chem_fastjx_driv.o ref_sounding.o ModBasicFields.o ModTuv2.7.o mem_grid.o \
	mem_rrtm.o 
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

ModUrbanCanopy.o : $(SURFACE)/ModUrbanCanopy.f90 node_mod.o mem_tend.o \
	ModBasicFields.o mem_grid.o ModTurbFields.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfFile.o : $(FDDA)/ModVarfFile.f90 ModRamsGrid.o ModDateUtils.o \
	rconstants.o node_mod.o chem1_list.o ModRamsReadHeader.o ModGrid.o ReadBcst.o \
	mem_varinit.o mem_aer1.o ModGetVar.o mem_grid.o ModRcio.o isan_coms.o parlibf.o \
	ModMessageSet.o ModBasicFields.o ModMicControl.o ModControlVars.o mem_chem1.o \
	mem_leaf.o mem_scratch.o ModGridTree.o ModVarfUpdate.o ModNudAnalysis.o \
	ref_sounding.o aer1_list.o $(UTILS_INCS)/files.h $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModWindFarm.o : $(WIND_FARM)/ModWindFarm.f90 ModNamelistFile.o ModDateUtils.o \
	rconstants.o node_mod.o mem_tend.o io_params.o ModTurbFields.o ModBasicFields.o \
	mem_grid.o module_wind_fitch.o $(UTILS_INCS)/files.h 
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

ModNestIntrp.o : $(NESTING)/ModNestIntrp.f90 ModNestFillDens.o ModRinit.o \
	rconstants.o mem_scratch.o grid_dims.o ref_sounding.o ModBasicFields.o \
	mem_nestb.o mem_grid.o 
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

ModNudRead.o : $(FDDA)/ModNudRead.f90 mem_varinit.o mem_chem1.o ModRamsGrid.o \
	ModNudUpdate.o ModDateUtils.o isan_coms.o ModNudAnalysis.o mem_grid.o \
	VarTable.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudUpdate.o : $(FDDA)/ModNudUpdate.f90 grid_struct.o mem_varinit.o \
	mem_chem1.o mem_aerad.o an_header.o chem1_list.o ModInitHis.o mem_grid.o \
	VarTable.o ModRcio.o $(UTILS_INCS)/files.h 
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

ModOdaRead.o : $(FDDA)/ModOdaRead.f90 ModDateUtils.o isan_coms.o \
	ModOdaStaCount.o ModOdaStaInput.o mem_oda.o mem_grid.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaCount.o : $(FDDA)/ModOdaStaCount.f90 mem_oda.o mem_grid.o \
	ModReadRalph.o obs_input.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaInput.o : $(FDDA)/ModOdaStaInput.f90 obs_input.o ModDateUtils.o \
	ModOdaStaCount.o mem_oda.o ModReadRalph.o mem_grid.o 
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

ModAcoustAdap.o : $(MODEL)/ModAcoustAdap.f90 ModRbnd.o rconstants.o node_mod.o \
	mem_scratch.o ModMessageSet.o mem_grid.o ModGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_carma.o : $(RADIATE)/rad_carma.F90 mem_aerad.o mem_chem1.o ModNamelistFile.o \
	mem_carma.o machine_arq.o mem_globaer.o mem_tuv.o mem_globrad.o rconstants.o \
	ModDateUtils.o mem_aer1.o mem_leaf.o node_mod.o chem1_list.o grid_dims.o \
	carma_fastjx.o aer1_list.o mem_grid.o ccatt_start.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_driv.o : $(RADIATE)/rad_driv.f90 ModNamelistFile.o ModMicroFields.o \
	CuParmFields.o RadiateFields.o ModBasicFields.o ModCarmaDriver.o \
	ModRrtmDriver.o ModMicControl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcAdap.o : $(MODEL)/ModRadvcAdap.f90 ModAdapInit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsGrid.o : $(INIT)/ModRamsGrid.f90 dump.o ModGridSet.o rconstants.o \
	ModAdapInit.o node_mod.o mem_grid.o $(UTILS_INCS)/constants.h 
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

ReadBcst.o : $(MPI)/ReadBcst.f90 ModControlVars.o mem_aerad.o an_header.o \
	utilsMod.o node_mod.o parlibf.o ModMPassFull.o ModBasicFields.o mem_grid.o \
	ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ref_sounding.o : $(MODEL)/ref_sounding.f90 grid_dims.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRnode.o : $(MODEL)/ModRnode.f90 mem_leaf.o node_mod.o parlibf.o grid_dims.o \
	mem_grid.o VarTable.o $(UTILS_INCS)/constants.h 
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

rrtmg_lw_cldprmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprmc.f90 parrrtm.o rrlw_vsn.o \
	parkind.o rrlw_wvn.o rrlw_cld.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprop.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprop.f90 rrlw_cld.o parkind.o \
	rrlw_vsn.o parrrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_init.o : $(RRTMG_LW_SRC)/rrtmg_lw_init.f90 rrlw_kg14.o rrlw_con.o \
	rrlw_kg16.o parkind.o rrlw_kg07.o rrlw_cld.o rrlw_kg15.o rrlw_kg02.o rrlw_tbl.o \
	rrlw_vsn.o rrlw_kg04.o rrlw_kg05.o rrlw_kg12.o rrlw_kg09.o rrlw_kg10.o \
	parrrtm.o rrlw_kg11.o rrlw_kg03.o rrtmg_lw_setcoef.o rrtmg_lw_k_g.o rrlw_kg06.o \
	rrlw_kg08.o rrlw_kg13.o rrlw_wvn.o rrlw_kg01.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_k_g.o : $(RRTMG_LW_SRC)/rrtmg_lw_k_g.f90 rrlw_kg15.o rrlw_kg02.o \
	rrlw_kg09.o rrlw_kg10.o rrlw_kg06.o rrlw_kg14.o rrlw_kg11.o rrlw_vsn.o \
	rrlw_kg08.o rrlw_kg04.o rrlw_kg16.o rrlw_kg03.o rrlw_kg13.o parkind.o \
	rrlw_kg05.o rrlw_kg01.o rrlw_kg07.o rrlw_kg12.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rad.o : $(RRTMG_LW_SRC)/rrtmg_lw_rad.f90 rrtmg_lw_cldprmc.o parrrtm.o \
	mcica_subcol_gen_lw.o rrlw_con.o parkind.o rrtmg_lw_taumol.o rrtmg_lw_setcoef.o \
	rrlw_wvn.o rrtmg_lw_rtrnmc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrn.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrn.f90 parrrtm.o rrlw_tbl.o \
	rrlw_con.o rrlw_vsn.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmc.f90 parrrtm.o rrlw_tbl.o \
	rrlw_con.o rrlw_vsn.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmr.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmr.f90 parrrtm.o rrlw_tbl.o \
	rrlw_con.o rrlw_vsn.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_setcoef.o : $(RRTMG_LW_SRC)/rrtmg_lw_setcoef.f90 rrlw_ref.o parrrtm.o \
	rrlw_vsn.o parkind.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_taumol.o : $(RRTMG_LW_SRC)/rrtmg_lw_taumol.f90 rrlw_kg14.o rrlw_con.o \
	rrlw_kg16.o parkind.o rrlw_kg07.o rrlw_kg15.o rrlw_kg02.o rrlw_ref.o rrlw_vsn.o \
	rrlw_kg04.o rrlw_kg05.o rrlw_kg12.o rrlw_kg09.o rrlw_kg10.o parrrtm.o \
	rrlw_kg11.o rrlw_kg03.o rrlw_kg06.o rrlw_kg08.o rrlw_kg13.o rrlw_wvn.o \
	rrlw_kg01.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprmc.f90 rrsw_wvn.o parrrsw.o \
	parkind.o rrsw_cld.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprop.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprop.f90 rrsw_wvn.o parrrsw.o \
	parkind.o rrsw_cld.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_init.o : $(RRTMG_SW_SRC)/rrtmg_sw_init.f90 rrsw_kg17.o rrsw_con.o \
	parkind.o parrrsw.o rrsw_kg16.o rrsw_kg20.o rrsw_vsn.o rrsw_kg28.o rrsw_kg26.o \
	rrtmg_sw_setcoef.o rrsw_kg25.o rrsw_wvn.o rrsw_cld.o rrsw_kg24.o rrsw_aer.o \
	rrsw_tbl.o rrsw_kg23.o rrsw_kg27.o rrtmg_sw_k_g.o rrsw_kg29.o rrsw_kg19.o \
	rrsw_kg22.o rrsw_kg18.o rrsw_kg21.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_k_g.o : $(RRTMG_SW_SRC)/rrtmg_sw_k_g.f90 rrsw_kg17.o rrsw_kg23.o \
	rrsw_kg25.o rrsw_kg27.o rrsw_kg22.o rrsw_kg18.o parkind.o rrsw_kg21.o \
	rrsw_kg16.o rrsw_kg20.o rrsw_kg24.o rrsw_vsn.o rrsw_kg29.o rrsw_kg28.o \
	rrsw_kg19.o rrsw_kg26.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_rad.o : $(RRTMG_SW_SRC)/rrtmg_sw_rad.f90 mcica_subcol_gen_sw.o \
	rrtmg_sw_cldprmc.o rrtmg_sw_setcoef.o rrsw_wvn.o rrsw_con.o parkind.o parrrsw.o \
	rrsw_aer.o rrtmg_sw_spcvmc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_reftra.o : $(RRTMG_SW_SRC)/rrtmg_sw_reftra.f90 rrsw_vsn.o rrsw_tbl.o \
	parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_setcoef.o : $(RRTMG_SW_SRC)/rrtmg_sw_setcoef.f90 rrsw_vsn.o parkind.o \
	parrrsw.o rrsw_ref.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvmc.f90 rrsw_tbl.o \
	rrtmg_sw_vrtqdr.o rrtmg_sw_reftra.o rrsw_wvn.o parkind.o parrrsw.o \
	rrtmg_sw_taumol.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvrt.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvrt.f90 rrsw_tbl.o \
	rrtmg_sw_vrtqdr.o rrtmg_sw_reftra.o rrsw_wvn.o parkind.o parrrsw.o \
	rrtmg_sw_taumol.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_taumol.o : $(RRTMG_SW_SRC)/rrtmg_sw_taumol.f90 rrsw_kg17.o rrsw_kg23.o \
	rrsw_kg25.o rrsw_kg27.o rrsw_wvn.o rrsw_kg22.o rrsw_kg18.o rrsw_con.o parkind.o \
	parrrsw.o rrsw_kg16.o rrsw_kg20.o rrsw_kg21.o rrsw_kg24.o rrsw_vsn.o \
	rrsw_kg29.o rrsw_kg28.o rrsw_kg19.o rrsw_kg26.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_vrtqdr.o : $(RRTMG_SW_SRC)/rrtmg_sw_vrtqdr.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRUser.o : $(SURFACE)/ModRUser.f90 rconstants.o memSoilMoisture.o node_mod.o \
	mem_leaf.o ModLeafComs.o io_params.o mem_grid.o ccatt_start.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

setup.o : $(MATRIX)/setup.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

sfclyr_jules.o : $(JULES_DIR)/sfclyr_jules.f90 rconstants.o node_mod.o \
	ModLeaf3Init.o chem1_list.o io_params.o jules_fields_mod.o ModMicroFields.o \
	ModLeafComs.o mem_grid.o ModTurbFields.o csigma_mod.o ModNamelistFile.o \
	sf_diags_mod.o CuParmFields.o RadiateFields.o io_constants.o ModBasicFields.o \
	gridmean_fluxes.o datetime_mod.o ModMicControl.o jules_surface_types_mod.o \
	mem_chem1.o model_time_mod.o mem_leaf.o mem_brams_jules.o gridbox_mean_mod.o \
	JulesFields.o ancil_info.o fluxes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

shcu_vars_const.o : $(CUPARM)/shcu_vars_const.f90 grid_dims.o ModNamelistFile.o \
	ModConvComs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSoilMoisture.o : $(SOIL_MOISTURE)/ModSoilMoisture.F90 ReadBcst.o \
	ModControlVars.o mem_aerad.o ModNamelistFile.o rconstants.o memSoilMoisture.o \
	mem_leaf.o node_mod.o ModLeafComs.o io_params.o parlibf.o ModMPassFull.o \
	ModBasicFields.o mem_grid.o ModTurbFields.o $(UTILS_INCS)/files.h \
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

ModTKenn.o : $(STILT)/ModTKenn.f90 mem_stilt.o rconstants.o mem_scratch.o \
	mem_grid.o turb_constants.o 
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

ModVarfUpdate.o : $(FDDA)/ModVarfUpdate.f90 rconstants.o mem_scratch.o \
	ModInitHis.o ref_sounding.o mem_grid.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

include jules_depend_model.mk

