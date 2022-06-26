actv.o : $(MATRIX)/actv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

adap_init.o : $(INIT)/adap_init.f90 mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

aer1_list.o : $(AEROSOL)/aer1_list_$(AERLEVEL).f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)
	@ln -fs aer1_list_$(AERLEVEL).o aer1_list.o

alloc.o : $(MEMORY)/alloc.F90 mem_grid.o mem_gaspart.o mem_globrad.o \
	mem_radiate.o mem_micro.o mem_oda.o teb_spm_start.o var_tables.o mem_scalar.o \
	mem_teb.o mem_jules.o mem_aerad.o mem_leaf.o mem_teb_common.o mem_varinit.o \
	mem_shcu.o mem_cuparm.o mem_scratch.o mem_globaer.o mem_opt_scratch.o \
	mem_tend.o 
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

ModAsTi.o : $(ISAN)/ModAsTi.f90 mem_grid.o rconstants.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAsTp.o : $(ISAN)/ModAsTp.f90 isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAVarF.o : $(ISAN)/ModAVarF.f90 mem_grid.o ModRbnd.o rconstants.o isan_coms.o 
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

chem_aobj.o : $(ISAN_CHEM)/chem_aobj.f90 isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_asti.o : $(ISAN_CHEM)/chem_asti.f90 mem_grid.o mem_chem1.o ModAsTi.o \
	chem_isan_coms.o mem_aer1.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_asti2.o : $(ISAN_CHEM)/chem_asti2.f90 ModDateUtils.o rconstants.o \
	isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_astp.o : $(ISAN_CHEM)/chem_astp.F90 ModDateUtils.o ModAsTp.o mem_chem1.o \
	rconstants.o mem_varinit.o chem_isan_coms.o dump.o isan_coms.o chem1_list.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_avarf.o : $(ISAN_CHEM)/chem_avarf.f90 mem_grid.o mem_chem1.o ModAVarF.o \
	ModRbnd.o rconstants.o chem_isan_coms.o mem_aer1.o isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_dry_dep.o : $(MODEL_CHEM)/chem_dry_dep.f90 ModDateUtils.o mem_chem1.o \
	aer1_list.o extra.o mem_aer1.o chem1_list.o 
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

chem_fastjx_driv.o : $(CCATT)/chem_fastjx_driv.f90 chem1_list.o chem_fastjx57.o \
	rconstants.o chem_fastjx_data.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_file_inv.o : $(ISAN_CHEM)/chem_file_inv.f90 mem_grid.o ModDateUtils.o \
	dump.o isan_coms.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_first_rams.o : $(ISAN_CHEM)/chem_first_rams.f90 mem_scratch.o ModRcio.o \
	grid_dims.o mem_grid.o an_header.o rconstants.o isan_coms.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_isan_coms.o : $(ISAN_CHEM)/chem_isan_coms.f90 isan_coms.o chem1_list.o \
	aer1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_isan_io.o : $(ISAN_CHEM)/chem_isan_io.f90 isan_coms.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_orage.o : $(CCATT)/chem_orage.f90 mem_scratch1_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_plumerise_scalar.o : $(CCATT)/chem_plumerise_scalar.f90 mem_aer1.o \
	mem_chem1.o mem_plume_chem1.o node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_refstate.o : $(ISAN_CHEM)/chem_refstate.f90 ccatt_start.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_sources.o : $(CCATT)/chem_sources.f90 ModDateUtils.o mem_plume_chem1.o \
	io_params.o parlibf.o mem_chem1.o mem_grid.o aer1_list.o ReadBcst.o mem_aer1.o \
	ModControlVars.o ModNamelistFile.o mem_volc_chem1.o $(UTILS_INCS)/constants.h 
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
	chem_spack_fexloss.o mem_chem1.o chem_spack_rates.o chem_spack_fexprod.o \
	chem_spack_kinetic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rates.o : $(MODEL_CHEM)/chem_spack_rates.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_rodas3_dyndt.o : $(CCATT)/chem_spack_rodas3_dyndt.f90 \
	chem_spack_fexchem.o mem_grid.o mem_chem1.o extra.o chem_spack_kinetic.o \
	mem_spack.o chem_spack_ros.o chem_spack_jacdchemdc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros.o : $(CCATT)/chem_spack_ros.f90 chem_spack_fexchem.o mem_chem1.o \
	chem_spack_kinetic.o mem_spack.o chem_spack_solve_sparse.o \
	chem_spack_jacdchemdc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

chem_spack_ros_dyndt.o : $(CCATT)/chem_spack_ros_dyndt.f90 chem_spack_fexchem.o \
	chem_spack_ros.o mem_chem1.o chem_spack_kinetic.o mem_spack.o \
	chem_spack_solve_sparse.o chem_spack_jacdchemdc.o 
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

chem_v_interps.o : $(ISAN_CHEM)/chem_v_interps.f90 dump.o rconstants.o \
	isan_coms.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemDryDepDriver.o : $(MODEL_CHEM)/ChemDryDepDriver.f90 grid_dims.o mem_grid.o \
	ModBasicFields.o mem_chem1.o chem_dry_dep.o mem_cuparm.o rconstants.o micphys.o \
	mem_leaf.o mem_aer1.o mem_radiate.o ModTurbFields.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ChemSourcesDriver.o : $(CCATT)/ChemSourcesDriver.f90 chem_plumerise_scalar.o \
	mem_plume_chem1.o io_params.o chem_sources.o mem_chem1.o mem_stilt.o \
	aer1_list.o mem_leaf.o mem_aer1.o mem_volc_chem1.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

coag.o : $(MATRIX)/coag.f90 memMatrix.o setup.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModCondRead.o : $(FDDA)/ModCondRead.f90 ModDateUtils.o ModNudAnalysis.o \
	mem_grid.o ModCondUpdate.o mem_varinit.o isan_coms.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCondUpdate.o : $(FDDA)/ModCondUpdate.f90 var_tables.o ModInitHis.o ModRcio.o \
	mem_grid.o an_header.o grid_struct.o mem_varinit.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

conv_coms.o : $(CUPARM)/conv_coms.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ConvPar_GF_GEOS5.o : $(CUPARM)/ConvPar_GF_GEOS5.F90 module_gate.o \
	Henrys_Law_cts.o MAPL_Constants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuRead.o : $(CUPARM)/ModCuRead.f90 mem_grid.o ModDateUtils.o mem_cuparm.o \
	isan_coms.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

cup_dn.o : $(CUPARM)/cup_dn.f90 
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

ModCupGrellCattDeep.o : $(CUPARM)/ModCupGrellCattDeep.f90 node_mod.o \
	mem_scratch2_grell.o mem_scratch3_grell.o mem_carma.o ccatt_start.o mem_grid.o \
	Phys_const.o kbcon_ecmwf.o mem_grell_param2.o mem_varinit.o cup_output_vars.o \
	ModCupUp.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCupGrellCattShallow.o : $(CUPARM)/ModCupGrellCattShallow.f90 node_mod.o \
	mem_scratch3_grell_sh.o mem_grid.o Phys_const.o mem_scratch2_grell_sh.o \
	mem_grell_param2.o mem_varinit.o cup_output_vars.o ModCupUp.o 
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

dam.o : $(ENERGY)/dam.f90 mem_grid.o ModDateUtils.o dump.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

debugTools.o : $(UTILS_TOOLS)/debugTools.f90 node_mod.o mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

depv.o : $(MATRIX)/depv.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

digitalFilter.o : $(MODEL)/digitalFilter.f90 node_mod.o ModDateUtils.o \
	var_tables.o io_params.o grid_dims.o mem_grid.o ModBasicFields.o ReadBcst.o \
	ModControlVars.o ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

domain_decomp.o : $(INIT)/domain_decomp.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

extra.o : $(MEMORY)/extra.f90 var_tables.o dump.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

getvar.o : $(UTILS_LIB)/getvar.f90 an_header.o dump.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
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

gridset.o : $(INIT)/gridset.f90 mem_grid.o grid_dims.o rconstants.o 
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

ModLeaf3Hyd.o : $(SURFACE)/ModLeaf3Hyd.f90 mem_grid.o leaf_coms.o mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Init.o : $(SURFACE)/ModLeaf3Init.f90 teb_spm_start.o io_params.o \
	grid_dims.o mem_grid.o leaf_coms.o ModLeaf3.o rconstants.o mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3Teb.o : $(SURFACE)/ModLeaf3Teb.f90 ModGasPart.o mem_teb_vars_const.o \
	ModUrban.o mem_emiss.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

leaf_coms.o : $(SURFACE)/leaf_coms.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

local_proc.o : $(MODEL)/local_proc.F90 node_mod.o io_params.o grid_dims.o \
	mem_grid.o mem_stilt.o rconstants.o ref_sounding.o dump.o ReadBcst.o \
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
	rrlw_con.o parkind.o mcica_random_numbers.o rrlw_wvn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_sw.o : $(RRTMG_SW_SRC)/mcica_subcol_gen_sw.f90 rrsw_wvn.o \
	parrrsw.o parkind.o rrsw_con.o mcica_random_numbers.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aer1.o : $(CCATT)/mem_aer1.f90 node_mod.o var_tables.o io_params.o \
	grid_dims.o ModScalarTable.o mem_grid.o mem_chem1.o aer1_list.o dump.o \
	ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aerad.o : $(RADIATE)/mem_aerad.f90 mem_grid_dim_defs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_carma.o : $(RADIATE)/mem_carma.f90 node_mod.o var_tables.o mem_scalar.o \
	grid_dims.o io_params.o mem_grid.o ModBasicFields.o parlibf.o mem_aerad.o \
	ModMPassFull.o ReadBcst.o mem_globrad.o ModControlVars.o ModTurbFields.o \
	ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1.o : $(CCATT)/mem_chem1.f90 var_tables.o io_params.o grid_dims.o \
	ModScalarTable.o ModNamelistFile.o chem1_list.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chem1aq.o : $(CCATT)/mem_chem1aq.f90 var_tables.o grid_dims.o \
	ModScalarTable.o mem_chem1.o chem1aq_list.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_chemic.o : $(CCATT)/mem_chemic.f90 micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_cuparm.o : $(CUPARM)/mem_cuparm.f90 var_tables.o ModNamelistFile.o \
	grid_dims.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
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

mem_gaspart.o : $(TEB_SPM)/mem_gaspart.f90 var_tables.o mem_emiss.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globaer.o : $(RADIATE)/mem_globaer.f90 mem_precision.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_globrad.o : $(RADIATE)/mem_globrad.f90 mem_precision.o parlibf.o mem_aerad.o \
	mem_radiate.o ModNamelistFile.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell.o : $(CUPARM)/mem_grell.f90 shcu_vars_const.o var_tables.o \
	mem_cuparm.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grell_param2.o : $(CUPARM)/mem_grell_param2.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid.o : $(MEMORY)/mem_grid.f90 var_tables.o ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_grid_dim_defs.o : $(MEMORY)/mem_grid_dim_defs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_jules.o : $(SURFACE)/mem_jules.f90 var_tables.o ModNamelistFile.o \
	grid_dims.o io_params.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_leaf.o : $(SURFACE)/mem_leaf.f90 teb_spm_start.o var_tables.o io_params.o \
	grid_dims.o ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_micro.o : $(MICRO)/mem_micro.f90 mem_radiate.o var_tables.o mem_cuparm.o \
	micphys.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_micro_optij.o : $(MICRO)/mem_micro_optij.f90 node_mod.o grid_dims.o \
	mem_grid.o rconstants.o micphys.o mem_radiate.o mem_micro.o 
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

mem_oda.o : $(FDDA)/mem_oda.f90 var_tables.o ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_opt_scratch.o : $(TURB)/mem_opt_scratch.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_plume_chem1.o : $(CCATT)/mem_plume_chem1.f90 mem_chem1.o var_tables.o \
	ModNamelistFile.o chem1_list.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_precision.o : $(RADIATE)/mem_precision.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_radiate.o : $(RADIATE)/mem_radiate.f90 var_tables.o ModNamelistFile.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_rrtm.o : $(RADIATE)/mem_rrtm.f90 node_mod.o rrtmg_lw_init.o mem_grid.o \
	mem_chem1.o parrrsw.o parrrtm.o rconstants.o parkind.o rrtmg_sw_init.o \
	chem1_list.o 
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

mem_scratch1_grell.o : $(CUPARM)/mem_scratch1_grell.f90 dump.o ccatt_start.o \
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

mem_shcu.o : $(CUPARM)/mem_shcu.f90 var_tables.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_spack.o : $(CCATT)/mem_spack.f90 chem1_list.o chem_spack_utils.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_stilt.o : $(STILT)/mem_stilt.f90 var_tables.o io_params.o grid_dims.o \
	rconstants.o ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tconv.o : $(CCATT)/mem_tconv.f90 mem_aer1.o chem1_list.o aer1_list.o 
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

mem_tend.o : $(MEMORY)/mem_tend.f90 teb_spm_start.o mem_scalar.o \
	ModScalarTable.o mem_grid.o ModBasicFields.o mem_gaspart.o mem_emiss.o \
	ModTurbFields.o ModNamelistFile.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_turb_scalar.o : $(TURB)/mem_turb_scalar.f90 var_tables.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_tuv.o : $(TUV)/mem_tuv.f90 ModTuv2.7.o var_tables.o mem_stilt.o \
	mem_globrad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

mem_varinit.o : $(MEMORY)/mem_varinit.f90 var_tables.o grid_dims.o mem_chem1.o \
	ModNamelistFile.o chem1_list.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_volc_chem1.o : $(CCATT)/mem_volc_chem1.f90 var_tables.o ModNamelistFile.o \
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

meteogram.o : $(IO)/meteogram.f90 node_mod.o ModPostUtils.o satPolyColision.o \
	meteogramType.o ModMPassDtl.o var_tables.o mem_grid.o ModNamelistFile.o \
	$(POST_INCS)/post_rconstants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

meteogramType.o : $(IO)/meteogramType.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicColl.o : $(MICRO)/ModMicColl.f90 micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicGamma.o : $(MICRO)/ModMicGamma.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicInit.o : $(MICRO)/ModMicInit.f90 mem_grid.o ModMicTabs.o rconstants.o \
	micphys.o dump.o ModMicGamma.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicNuc.o : $(MICRO)/ModMicNuc.f90 micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicTabs.o : $(MICRO)/ModMicTabs.f90 ModMicGamma.o micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicVap.o : $(MICRO)/ModMicVap.f90 rconstants.o micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

micphys.o : $(MICRO)/micphys.f90 ModNamelistFile.o grid_dims.o \
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

ModAcoust.o : $(MODEL)/ModAcoust.f90 node_mod.o mem_scratch.o raco_adap.o \
	mem_grid.o ModBasicFields.o rconstants.o micphys.o ref_sounding.o ModGrid.o \
	ModParallelEnvironment.o ModMessageSet.o mem_tend.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModAerClim.o : $(AERCLIM)/ModAerClim.f90 node_mod.o mem_grid.o ModBasicFields.o \
	parlibf.o dump.o ReadBcst.o ModControlVars.o ModTurbFields.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBasicFields.o : $(MEMORY)/ModBasicFields.f90 var_tables.o mem_stilt.o \
	ModNodeDimensions.o ModParallelEnvironment.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBramsGrid.o : $(POST_SRC)/ModBramsGrid.f90 node_mod.o ModPostUtils.o \
	mem_grid.o mem_aerad.o ref_sounding.o ModParallelEnvironment.o \
	ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModBuffering.o : $(MPI)/ModBuffering.f90 ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCarmaDriver.o : $(RADIATE)/ModCarmaDriver.f90 node_mod.o ModDateUtils.o \
	teb_spm_start.o grid_dims.o mem_carma.o mem_grid.o ModLeaf3.o ModBasicFields.o \
	mem_cuparm.o rconstants.o micphys.o mem_scratch1_grell.o mem_leaf.o \
	mem_teb_common.o mem_tend.o mem_radiate.o mem_micro.o rad_carma.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemAsgen.o : $(ISAN_CHEM)/ModChemAsgen.F90 node_mod.o ModDateUtils.o \
	io_params.o grid_dims.o mem_grid.o mem_chem1.o aer1_list.o isan_coms.o \
	ModAsGen.o chem_isan_coms.o dump.o mem_aer1.o ModControlVars.o ModMkSfcTop.o \
	chem1_list.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemConvTransp.o : $(CCATT)/ModChemConvTransp.f90 node_mod.o mem_scratch.o \
	mem_grid.o mem_chem1.o Phys_const.o mem_cuparm.o aer1_list.o \
	mem_scratch1_grell.o mem_grell_param2.o mem_tconv.o mem_aer1.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModChemistryDriver.o : $(CCATT)/ModChemistryDriver.f90 chem_spack_ros_dyndt.o \
	mem_carma.o mem_grid.o mem_chem1aq.o parrrtm.o aer1_list.o chem_orage.o \
	chem_spack_rodas3_dyndt.o mem_globrad.o mem_radiate.o chem_spack_ros.o \
	ModTuvDriver2.7.o carma_fastjx.o mem_micro.o ModBasicFields.o \
	chem_fastjx_driv.o mem_stilt.o rconstants.o micphys.o chem_trans_liq.o \
	mem_aerad.o mem_spack.o mem_aer1.o grid_dims.o mem_rrtm.o chem_trans_gasaq.o \
	mem_cuparm.o extra.o mem_scratch1_grell.o mem_grell_param2.o \
	chem_spack_solve_sparse.o chem_spack_utils.o node_mod.o mem_scratch.o \
	mem_chem1.o chem1aq_list.o chem_spack_qssa.o mem_chemic.o chem_uv_att.o \
	chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModControlVars.o : $(INIT)/ModControlVars.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCoriolis.o : $(MODEL)/ModCoriolis.f90 node_mod.o mem_scratch.o mem_grid.o \
	parlibf.o ModBasicFields.o ModBuffering.o rconstants.o ref_sounding.o \
	mem_tend.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModCuParGrell3.o : $(CUPARM)/ModCuParGrell3.F90 module_cu_gf_v5.1.o ModRadvc.o \
	mem_carma.o mem_grid.o Phys_const.o ModRstilt.o mem_radiate.o ModMessageSet.o \
	mem_micro.o ConvPar_GF_GEOS5.o var_tables.o ccatt_start.o mem_jules.o \
	ModBasicFields.o mem_stilt.o ModChemConvTransp.o rconstants.o micphys.o \
	module_cu_gf.o mem_leaf.o mem_varinit.o mem_grell.o io_params.o grid_dims.o \
	mem_cuparm.o module_cu_g3.o mem_scratch1_grell.o mem_grell_param2.o \
	ModNamelistFile.o node_mod.o mem_scratch.o mem_chem1.o ModGrid.o mem_tend.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDateUtils.o : $(UTILS_MODS)/ModDateUtils.f90 dump.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffSclr.o : $(TURB)/ModDiffSclr.f90 mem_grid.o mem_scratch.o ModTurbDiff.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDiffuse.o : $(TURB)/ModDiffuse.f90 node_mod.o ModDiffSclr.o mem_scratch.o \
	ke_coms.o ModTurbDiff.o ModScalarTable.o mem_grid.o ModBasicFields.o \
	mem_opt_scratch.o ModTurbKE.o micphys.o mem_leaf.o mem_tend.o ModTurbK.o \
	ModTurbFields.o ModNamelistFile.o mem_micro.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModDomainDecomp.o : $(MPI)/ModDomainDecomp.f90 ModGridDims.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModEvaluation.o : $(EVAL)/ModEvaluation.f90 parlibf.o node_mod.o \
	ModNamelistFile.o mem_grid.o 
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

ModGasPart.o : $(TEB_SPM)/ModGasPart.f90 node_mod.o var_tables.o ModRcio.o \
	grid_dims.o mem_teb_vars_const.o mem_grid.o parlibf.o an_header.o \
	ModBasicFields.o mem_gaspart.o mem_emiss.o mem_leaf.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGeodat.o : $(MKSFC)/ModGeodat.f90 mem_grid.o teb_spm_start.o io_params.o \
	mem_leaf.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModGrid.o : $(MPI)/ModGrid.f90 var_tables.o meteogramType.o ModScalarTable.o \
	ModBasicFields.o ModGridDims.o ModNeighbourNodes.o ModNodeDimensions.o \
	ModDomainDecomp.o ModMessageSet.o ModParallelEnvironment.o ModControlVars.o \
	ModTurbFields.o ModNamelistFile.o mem_tend.o 
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

modIau.o : $(MODEL)/modIau.f90 node_mod.o mem_grid.o parlibf.o ModMPassFull.o \
	dump.o ReadBcst.o mem_varinit.o ModNamelistFile.o mem_tend.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitHis.o : $(IO)/ModInitHis.f90 var_tables.o mem_scratch.o ModRcio.o \
	io_params.o mem_grid.o ModBasicFields.o an_header.o mem_chem1.o ModLeaf3.o \
	leaf_coms.o ModRinit.o rconstants.o ref_sounding.o micphys.o mem_aerad.o \
	mem_leaf.o mem_varinit.o chem1_list.o ModRamsReadHeader.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModInitMicThompson.o : $(MICRO)/ModInitMicThompson.f90 node_mod.o ModDateUtils.o \
	mem_grid.o ModBasicFields.o parlibf.o generic.o dump.o ReadBcst.o mem_micro.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLanduseInput.o : $(MKSFC)/ModLanduseInput.f90 mem_mksfc.o ModLeaf3Init.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3.o : $(SURFACE)/ModLeaf3.f90 node_mod.o teb_spm_start.o mem_scratch.o \
	io_params.o ccatt_start.o mem_teb.o ModLeaf3Teb.o leaf_coms.o ModBasicFields.o \
	mem_grid.o ModLeaf3Hyd.o mem_cuparm.o rconstants.o micphys.o mem_leaf.o \
	mem_teb_common.o mem_radiate.o ModTurbFields.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModLeaf3OceanOnly.o : $(SURFACE)/ModLeaf3OceanOnly.f90 node_mod.o \
	ConvPar_GF_GEOS5.o ModCuParGrell3.o io_params.o ccatt_start.o mem_grid.o \
	ModBasicFields.o leaf_coms.o ModLeaf3.o mem_cuparm.o rconstants.o micphys.o \
	mem_leaf.o mem_radiate.o ModTurbFields.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMatrixDriver.o : $(MATRIX)/ModMatrixDriver.F90 node_mod.o chem1_list.o \
	setup.o coag.o mem_grid.o mem_chem1.o ModBasicFields.o isrpia.o aer1_list.o \
	rconstants.o micphys.o npf.o memMatrix.o mem_leaf.o subs.o mem_aer1.o \
	mem_radiate.o ModTurbFields.o ModParticle.o mem_micro.o 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

ModMemAlloc.o : $(MEMORY)/ModMemAlloc.F90 mem_scratch2_grell.o mem_carma.o \
	mem_grid.o parrrsw.o mem_chem1aq.o aer1_list.o mem_grid_dim_defs.o ModOptical.o \
	mem_gaspart.o mem_tend.o mem_globrad.o mem_radiate.o carma_fastjx.o mem_oda.o \
	mem_volc_chem1.o mem_micro.o teb_spm_start.o var_tables.o ModCuParGrell3.o \
	mem_scalar.o mem_teb_vars_const.o mem_teb.o ccatt_start.o mem_jules.o \
	ModBasicFields.o chem_sources.o digitalFilter.o mem_stilt.o mem_aerad.o \
	mem_scratch1_brams.o micphys.o mem_scratch2_grell_sh.o mem_leaf.o modIau.o \
	mem_teb_common.o mem_varinit.o mem_aer1.o ModTurbFields.o mem_shcu.o \
	mem_grell.o mem_micro_optij.o io_params.o grid_dims.o mem_turb_scalar.o \
	chem_dry_dep.o mem_cuparm.o mem_tuv.o extra.o mem_scratch1_grell.o \
	mem_grell_param2.o shcu_vars_const.o node_mod.o ModTuv2.7.o mem_scratch.o \
	mem_plume_chem1.o mem_scratch3_grell.o mem_scratch3_grell_sh.o mem_chem1.o \
	mem_globaer.o chem1aq_list.o mem_opt_scratch.o mem_chemic.o ModGrid.o \
	ModParallelEnvironment.o mem_emiss.o machine_arq.o ModEvaluation.o mem_nestb.o \
	chem1_list.o 
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

ModMessageSet.o : $(MPI)/ModMessageSet.f90 var_tables.o ModMessageData.o \
	mem_grid.o parlibf.o ModNeighbourNodes.o ModGridDims.o ModFieldSection.o \
	ModNodeDimensions.o ModFieldSectionList.o ModDomainDecomp.o \
	ModParallelEnvironment.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicGfdlDriver.o : $(MICRO)/ModMicGfdlDriver.f90 node_mod.o io_params.o \
	mem_grid.o ModBasicFields.o rconstants.o mem_leaf.o mem_radiate.o mem_micro.o \
	gfdl_cloud_microphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsDrive.o : $(MICRO)/ModMicrophysicsDrive.f90 node_mod.o \
	ModMicVap.o mem_grid.o ModBasicFields.o ModMicTabs.o ModMicColl.o ModMicNuc.o \
	mem_chem1.o mem_chemic.o micphys.o mem_chem1aq.o ModMicrophysicsMisc.o \
	mem_radiate.o ModMicInit.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicrophysicsMisc.o : $(MICRO)/ModMicrophysicsMisc.f90 mem_scratch.o \
	mem_grid.o ModBasicFields.o rconstants.o micphys.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMicThompsonDriver.o : $(MICRO)/ModMicThompsonDriver.f90 node_mod.o \
	module_mp_thompson.o io_params.o mem_grid.o ModBasicFields.o rconstants.o \
	micphys.o mem_leaf.o mem_radiate.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcDriver.o : $(MKSFC)/ModMkSfcDriver.f90 node_mod.o teb_spm_start.o \
	ModNestGeoSst.o ModMkSfcSfc.o ModMkSfcFuso.o grid_dims.o io_params.o mem_grid.o \
	ModMkSfcSst.o ModNdviRead.o ModMkSfcNdvi.o ReadBcst.o mem_mksfc.o \
	ModControlVars.o ModSstRead.o ModMkSfcTop.o ModLanduseInput.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcFuso.o : $(MKSFC)/ModMkSfcFuso.f90 node_mod.o io_params.o mem_teb.o \
	mem_teb_vars_const.o mem_grid.o mem_gaspart.o mem_emiss.o ReadBcst.o \
	mem_mksfc.o ModControlVars.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcNdvi.o : $(MKSFC)/ModMkSfcNdvi.f90 io_params.o ModRUser.o mem_grid.o \
	dump.o mem_leaf.o mem_mksfc.o ModLanduseInput.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSfc.o : $(MKSFC)/ModMkSfcSfc.f90 node_mod.o io_params.o mem_grid.o \
	dump.o mem_leaf.o ReadBcst.o mem_mksfc.o ModControlVars.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcSst.o : $(MKSFC)/ModMkSfcSst.f90 io_params.o grid_dims.o ModRUser.o \
	mem_grid.o ModGeodat.o mem_leaf.o mem_mksfc.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMkSfcTop.o : $(MKSFC)/ModMkSfcTop.f90 node_mod.o io_params.o mem_grid.o \
	dump.o ReadBcst.o mem_mksfc.o ModControlVars.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModMonotonicAdvection.o : $(MODEL)/ModMonotonicAdvection.f90 ccatt_start.o \
	mem_grid.o mem_chem1.o chem_dry_dep.o rconstants.o ModDomainDecomp.o micphys.o \
	ModGrid.o ModParallelEnvironment.o mem_aer1.o ModMessageSet.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNamelistFile.o : $(INIT)/ModNamelistFile.f90 grid_dims.o parlibf.o \
	modPrintInitial.o dump.o ModParallelEnvironment.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNdviRead.o : $(MKSFC)/ModNdviRead.f90 node_mod.o ModDateUtils.o io_params.o \
	grid_dims.o mem_grid.o ModMkSfcNdvi.o mem_leaf.o ReadBcst.o mem_mksfc.o \
	ModControlVars.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNeighbourNodes.o : $(MPI)/ModNeighbourNodes.f90 ModGridDims.o \
	ModParallelEnvironment.o ModDomainDecomp.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNestGeoSst.o : $(MKSFC)/ModNestGeoSst.f90 node_mod.o mem_scratch.o \
	ModInitHis.o io_params.o ccatt_start.o ModRUser.o grid_dims.o mem_grid.o \
	ModBasicFields.o ModMkSfcTop.o ModGeodat.o soilMoisture.o ModLeaf3Init.o \
	mem_leaf.o dump.o memSoilMoisture.o mem_mksfc.o ModControlVars.o \
	ModTurbFields.o ModLanduseInput.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNodeDimensions.o : $(MPI)/ModNodeDimensions.f90 ModGridDims.o \
	ModParallelEnvironment.o ModDomainDecomp.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudAnalysis.o : $(FDDA)/ModNudAnalysis.f90 node_mod.o mem_scratch.o \
	mem_grid.o ModBasicFields.o mem_chem1.o dump.o modIau.o mem_tend.o \
	mem_varinit.o ModEvaluation.o chem1_list.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaNudge.o : $(FDDA)/ModOdaNudge.f90 node_mod.o mem_scratch.o io_params.o \
	ModOdaKrig.o ModBasicFields.o mem_grid.o ModOdaProcObs.o mem_oda.o mem_tend.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOneProc.o : $(MODEL)/ModOneProc.F90 ModMkSfcFuso.o mem_chem1aq.o ModAerClim.o \
	aer1_list.o ModSched.o mem_gaspart.o meteogram.o mem_radiate.o \
	ModInitMicThompson.o digitalFilter.o ModCuRead.o mem_scalar.o ModWindFarm.o \
	ModBasicFields.o ModMkSfcTop.o isan_coms.o ModInitHis.o io_params.o dump.o \
	ModMicrophysicsMisc.o ModRhhi.o mem_chem1.o ModParaInit.o ref_sounding.o \
	ModGrid.o ModSstRead.o ModCondRead.o ModRio.o dam.o ModTimestepRK.o \
	ModVarfUpdate.o ModRamsMicrophysics2M.o mem_oda.o var_tables.o ModMkSfcDriver.o \
	ModOpspec.o mem_stilt.o micphys.o mem_aer1.o ModMkSfcSfc.o parlibf.o \
	chem_dry_dep.o mem_cuparm.o ModDomainDecomp.o extra.o ModChemAsgen.o \
	mem_grell_param2.o shcu_vars_const.o ModMicInit.o ModRThrm.o domain_decomp.o \
	mem_emiss.o ModParallelEnvironment.o ReadBcst.o machine_arq.o ModCoriolis.o \
	ModRUser.o mem_grid.o ModGridTree.o ModMemAlloc.o ModChemistryDriver.o \
	mem_teb_vars_const.o ModCuParGrell3.o ModNudRead.o mem_leaf.o mem_teb_common.o \
	mem_varinit.o local_proc.o ModTimestep.o memSoilMoisture.o ModTuv2.7.o \
	mem_scratch.o ModLeaf3Init.o ModRanlavg.o ModEvaluation.o ModPostProcess.o \
	chem1_list.o mem_carma.o ModTimeStamp.o mem_globrad.o ModTuvDriver2.7.o \
	mem_volc_chem1.o mem_micro.o ModNdviRead.o teb_spm_start.o ModMPassDtl.o \
	mem_teb.o chem_sources.o ccatt_start.o ModPostGridNetCDF.o modIau.o \
	ModGasPart.o ModRecycle.o tuvParameter.o ModNestGeoSst.o grid_dims.o \
	ModMonotonicAdvection.o ModRinit.o ModVarfFile.o ModNamelistFile.o node_mod.o \
	ModUrbanCanopy.o mem_plume_chem1.o ModLeaf3Teb.o ModOdaRead.o soilMoisture.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOpspec.o : $(IO)/ModOpspec.f90 mem_grid.o mem_chem1aq.o aer1_list.o \
	mem_globrad.o mem_radiate.o teb_spm_start.o ccatt_start.o chem_sources.o \
	mem_stilt.o micphys.o mem_leaf.o modIau.o mem_varinit.o mem_aer1.o io_params.o \
	grid_dims.o mem_cuparm.o mem_grell_param2.o shcu_vars_const.o ModNamelistFile.o \
	mem_chem1.o chem1aq_list.o mem_emiss.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOptical.o : $(RADIATE)/ModOptical.f90 node_mod.o var_tables.o ccatt_start.o \
	mem_grid.o ModBasicFields.o parlibf.o ModNamelistFile.o aer1_list.o \
	ModMPassFull.o dump.o mem_leaf.o ReadBcst.o mem_aer1.o mem_radiate.o \
	ModTurbFields.o ModControlVars.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOutputUtils.o : $(IO)/ModOutputUtils.f90 var_tables.o ModBasicFields.o dump.o \
	ModTurbFields.o ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOzone.o : $(TEB_SPM)/ModOzone.f90 mem_grid.o ModBasicFields.o rconstants.o \
	mem_gaspart.o mem_radiate.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModParaInit.o : $(MPI)/ModParaInit.f90 node_mod.o var_tables.o grid_dims.o \
	ModScalarTable.o mem_grid.o dump.o $(UTILS_INCS)/constants.h 
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

ModPostGrid.o : $(POST_SRC)/ModPostGrid.F90 ModOutputUtils.o ModPostUtils.o \
	ModPostOneFieldNetCDF.o io_params.o mem_grid.o ModBasicFields.o parlibf.o \
	ModNamelistFile.o ModParallelEnvironment.o ModBramsGrid.o ModTurbFields.o \
	ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostGridNetCDF.o : $(POST_SRC)/ModPostGridNetCDF.F90 ModPostUtils.o \
	ModDateUtils.o io_params.o ModPostTypes.o mem_grid.o dump.o ModBramsGrid.o \
	ModNamelistFile.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField.o : $(POST_SRC)/ModPostOneField.f90 ModPostOneField3d.o \
	ModPostUtils.o node_mod.o ModPostTypes.o ModPostOneField2d.o ModBasicFields.o \
	ModPostOneField7d.o dump.o ModBramsGrid.o ModPostOneFieldUtils.o \
	ModTurbFields.o ModNamelistFile.o ModPostOneField8d.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField2d.o : $(POST_SRC)/ModPostOneField2d.f90 ModOutputUtils.o \
	ModPostUtils.o node_mod.o io_params.o mem_grid.o ModBasicFields.o \
	ModNamelistFile.o ModPostGrid.o mem_aerad.o micphys.o mem_cuparm.o dump.o \
	ModBramsGrid.o ModPostOneFieldUtils.o mem_radiate.o ModTurbFields.o \
	ModPostTypes.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField3d.o : $(POST_SRC)/ModPostOneField3d.f90 node_mod.o \
	ModPostUtils.o ModOutputUtils.o mem_grid.o ModBasicFields.o ModNamelistFile.o \
	ModPostGrid.o micphys.o ModBramsGrid.o ModPostOneFieldUtils.o mem_varinit.o \
	ModTurbFields.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField7d.o : $(POST_SRC)/ModPostOneField7d.f90 ModOutputUtils.o \
	ModPostUtils.o ModBasicFields.o ModNamelistFile.o ModBramsGrid.o \
	ModPostOneFieldUtils.o ModTurbFields.o ModPostTypes.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostOneField8d.o : $(POST_SRC)/ModPostOneField8d.f90 ModOutputUtils.o \
	ModPostUtils.o ModBasicFields.o ModNamelistFile.o ModBramsGrid.o \
	ModPostOneFieldUtils.o ModTurbFields.o ModPostTypes.o 
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

ModPostProcess.o : $(POST_SRC)/ModPostProcess.F90 io_params.o \
	ModPostGridNetCDF.o ModBasicFields.o ModPostOneField.o ModNamelistFile.o \
	ModPostGrid.o ModGridTree.o ModGrid.o ModParallelEnvironment.o ModBramsGrid.o \
	ModMessageSet.o ModTurbFields.o ModPostTypes.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostTypes.o : $(POST_SRC)/ModPostTypes.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModPostUtils.o : $(POST_SRC)/ModPostUtils.f90 dump.o ModParallelEnvironment.o \
	mem_leaf.o $(UTILS_INCS)/constants.h $(POST_INCS)/post_rconfig.h \
	$(POST_INCS)/post_rconstants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

modPrintInitial.o : $(INIT)/modPrintInitial.F90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvc.o : $(MODEL)/ModRadvc.f90 mem_scratch.o grid_dims.o ccatt_start.o \
	ModScalarTable.o mem_grid.o ModBasicFields.o ModMonotonicAdvection.o \
	mem_chem1.o chem_dry_dep.o ModParallelEnvironment.o mem_aer1.o \
	ModNamelistFile.o mem_tend.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRadvcRK.o : $(MODEL)/ModRadvcRK.f90 node_mod.o grid_dims.o ModRexev.o \
	mem_grid.o mem_chem1.o mem_stilt.o ModGrid.o ModParallelEnvironment.o \
	ModMessageSet.o mem_tend.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsMicrophysics2M.o : $(MICRO)/ModRamsMicrophysics2M.f90 node_mod.o \
	mem_scratch.o grid_dims.o mem_grid.o ModBasicFields.o rconstants.o micphys.o \
	dump.o mem_leaf.o ModMicGamma.o mem_micro.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRamsReadHeader.o : $(IO)/ModRamsReadHeader.f90 an_header.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRanlavg.o : $(IO)/ModRanlavg.f90 var_tables.o io_params.o grid_dims.o \
	ModRThrm.o mem_grid.o ModBasicFields.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/interface.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRbnd.o : $(BC)/ModRbnd.f90 node_mod.o mem_scratch.o ccatt_start.o \
	ModScalarTable.o mem_grid.o ModBasicFields.o mem_chem1.o ModTurbKE.o micphys.o \
	ref_sounding.o ModMicrophysicsMisc.o ModTurbFields.o mem_tend.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRcio.o : $(IO)/ModRcio.f90 io_params.o grid_dims.o mem_grid.o leaf_coms.o \
	an_header.o mem_stilt.o mem_cuparm.o micphys.o ref_sounding.o \
	ModParallelEnvironment.o mem_leaf.o mem_radiate.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConv.o : $(CUPARM)/ModRConv.f90 node_mod.o mem_scratch.o conv_coms.o \
	mem_grid.o ModBasicFields.o mem_cuparm.o rconstants.o mem_tend.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRConvGrellCatt.o : $(CUPARM)/ModRConvGrellCatt.f90 mem_grid.o \
	ModCupGrellCattShallow.o ModRstilt.o mem_radiate.o ModCupGrellCattDeep.o \
	mem_micro.o ModCuParGrell3.o mem_scalar.o ccatt_start.o ModChemConvTransp.o \
	mem_stilt.o rconstants.o micphys.o mem_leaf.o mem_grell.o io_params.o \
	mem_cuparm.o mem_scratch1_grell.o mem_grell_param2.o node_mod.o mem_scratch.o \
	ModGrid.o mem_tend.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRecycle.o : $(IO)/ModRecycle.f90 node_mod.o ModDateUtils.o var_tables.o \
	io_params.o mem_grid.o mem_chem1.o an_header.o aer1_list.o mem_aerad.o \
	ModMPassFull.o dump.o ReadBcst.o mem_aer1.o chem1_list.o ModRamsReadHeader.o \
	$(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRexev.o : $(STILT)/ModRexev.f90 mem_scratch.o ModRadvc.o mem_grid.o \
	ModBasicFields.o mem_stilt.o rconstants.o micphys.o mem_micro.o mem_tend.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRGrad.o : $(TURB)/ModRGrad.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRhhi.o : $(INIT)/ModRhhi.f90 mem_scratch.o mem_grid.o ModBasicFields.o \
	ModRinit.o ref_sounding.o rconstants.o micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRinit.o : $(INIT)/ModRinit.f90 node_mod.o mem_scratch.o io_params.o \
	mem_grid.o ModBasicFields.o ModRbnd.o ModTurbKE.o ref_sounding.o rconstants.o \
	micphys.o mem_varinit.o ModTurbFields.o mem_micro.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRio.o : $(IO)/ModRio.f90 node_mod.o ModDateUtils.o var_tables.o \
	mpi_io_engine-5d.o io_params.o grid_dims.o ModRcio.o mem_grid.o \
	ModBasicFields.o an_header.o parlibf.o mem_chem1.o mem_aerad.o ref_sounding.o \
	ModMPassFull.o ModParallelEnvironment.o ReadBcst.o ModControlVars.o \
	ModTurbFields.o ModNamelistFile.o $(UTILS_INCS)/constants.h \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/interface.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRrtmDriver.o : $(RADIATE)/ModRrtmDriver.f90 ModDateUtils.o mem_carma.o \
	mem_grid.o rrtmg_sw_cldprop.o parrrsw.o parrrtm.o parkind.o ModOptical.o \
	mem_radiate.o mem_micro.o teb_spm_start.o ccatt_start.o leaf_coms.o \
	ModBasicFields.o rconstants.o micphys.o mcica_subcol_gen_sw.o mem_leaf.o \
	rrtmg_sw_rad.o mcica_subcol_gen_lw.o mem_rrtm.o grid_dims.o rrtmg_lw_cldprop.o \
	mem_cuparm.o mem_tuv.o mem_scratch1_grell.o mem_grell_param2.o node_mod.o \
	mem_chem1.o rrtmg_lw_rad.o ref_sounding.o mem_tend.o \
	$(UTILS_INCS)/aerosol_setup.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRShCuPar.o : $(CUPARM)/ModRShCuPar.f90 node_mod.o mem_scratch.o conv_coms.o \
	mem_grid.o ModBasicFields.o ModRConv.o mem_micro.o shcu_vars_const.o \
	ModTurbFields.o mem_shcu.o mem_tend.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRstilt.o : $(STILT)/ModRstilt.f90 mem_scratch.o grid_dims.o mem_grid.o \
	ModBasicFields.o ModMonotonicAdvection.o mem_stilt.o mem_cuparm.o \
	mem_scratch1_grell.o ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRThrm.o : $(MODEL)/ModRThrm.f90 mem_scratch.o mem_grid.o ModBasicFields.o \
	rconstants.o micphys.o mem_micro.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRtimi.o : $(MODEL)/ModRtimi.f90 node_mod.o mem_scratch.o ModScalarTable.o \
	mem_grid.o ModBasicFields.o mem_cuparm.o mem_tend.o shcu_vars_const.o \
	mem_grell.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModScalarTable.o : $(MEMORY)/ModScalarTable.f90 ModParallelEnvironment.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSched.o : $(MODEL)/ModSched.f90 node_mod.o io_params.o mem_grid.o parlibf.o \
	ModBasicFields.o local_proc.o mem_cuparm.o ref_sounding.o mem_varinit.o dump.o \
	ReadBcst.o shcu_vars_const.o mem_radiate.o ModNamelistFile.o isan_coms.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSeaSalt.o : $(CCATT)/ModSeaSalt.f90 node_mod.o io_params.o ccatt_start.o \
	mem_grid.o ModBasicFields.o ModAerClim.o mem_chem1.o aer1_list.o mem_leaf.o \
	mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModSstRead.o : $(MKSFC)/ModSstRead.f90 node_mod.o ModDateUtils.o io_params.o \
	grid_dims.o mem_grid.o ModMkSfcSst.o mem_leaf.o ReadBcst.o mem_mksfc.o \
	ModControlVars.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimeStamp.o : $(MODEL)/ModTimeStamp.f90 $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestep.o : $(MODEL)/ModTimestep.F90 ModCoriolis.o ModRadvc.o \
	ModNudAnalysis.o mem_grid.o ChemDryDepDriver.o sfclyr_jules.o \
	ModMicGfdlDriver.o ModRstilt.o ModTimeStamp.o ModOptical.o mem_radiate.o \
	ModMessageSet.o ModRamsMicrophysics2M.o rad_driv.o mem_oda.o ModSeaSalt.o \
	teb_spm_start.o ModDiffuse.o ModCuParGrell3.o mem_scalar.o ccatt_start.o \
	ModRexev.o chem_sources.o ModChemistryDriver.o ModMicrophysicsDrive.o \
	ModLeaf3.o ModAcoust.o ModBasicFields.o mem_stilt.o rconstants.o micphys.o \
	ModWindFarm.o mem_leaf.o ModRConvGrellCatt.o ModGasPart.o ModRtimi.o \
	mem_varinit.o mem_aer1.o ModOzone.o ModMicThompsonDriver.o grid_dims.o \
	ModRbnd.o ModMonotonicAdvection.o digitalFilter.o mem_cuparm.o ModRShCuPar.o \
	ModRConv.o ModMicrophysicsMisc.o shcu_vars_const.o ModMatrixDriver.o node_mod.o \
	mem_scratch.o ModUrbanCanopy.o mem_plume_chem1.o ModRThrm.o ModTurbK.o \
	mem_chem1.o ChemSourcesDriver.o ModOdaNudge.o ModGrid.o mem_emiss.o \
	machine_arq.o mem_tend.o $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTimestepRK.o : $(MODEL)/ModTimestepRK.F90 ModCoriolis.o utilsMod.o ModRadvc.o \
	ModNudAnalysis.o mem_grid.o ChemDryDepDriver.o ModAerClim.o sfclyr_jules.o \
	ModMicGfdlDriver.o ModRstilt.o ModTimeStamp.o ModOptical.o mem_radiate.o \
	ModMessageSet.o ModRamsMicrophysics2M.o rad_driv.o mem_oda.o ModSeaSalt.o \
	teb_spm_start.o ModDiffuse.o ModCuParGrell3.o mem_scalar.o ccatt_start.o \
	ModRexev.o chem_sources.o ModChemistryDriver.o ModMicrophysicsDrive.o \
	ModLeaf3.o ModAcoust.o mem_stilt.o ModWindFarm.o rconstants.o micphys.o \
	mem_leaf.o ModRtimi.o ModGasPart.o modIau.o mem_varinit.o mem_aer1.o ModOzone.o \
	ModMicThompsonDriver.o grid_dims.o ModRbnd.o ModMonotonicAdvection.o \
	digitalFilter.o mem_cuparm.o ModRShCuPar.o ModRConv.o ModTimestep.o \
	ModRadvcRK.o ModMicrophysicsMisc.o shcu_vars_const.o ModMatrixDriver.o \
	ModLeaf3OceanOnly.o node_mod.o mem_scratch.o ModUrbanCanopy.o mem_plume_chem1.o \
	ModRThrm.o ModTurbK.o mem_chem1.o ChemSourcesDriver.o ModOdaNudge.o ModGrid.o \
	ModParallelEnvironment.o mem_emiss.o machine_arq.o mem_tend.o \
	$(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiff.o : $(TURB)/ModTurbDiff.f90 mem_scratch.o mem_grid.o mem_cuparm.o \
	mem_opt_scratch.o ModRGrad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbDiffAdap.o : $(TURB)/ModTurbDiffAdap.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbFields.o : $(TURB)/ModTurbFields.f90 var_tables.o ModNodeDimensions.o \
	ModParallelEnvironment.o ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbK.o : $(TURB)/ModTurbK.f90 mem_grid.o ModRstilt.o mem_micro.o ke_coms.o \
	ccatt_start.o ModBasicFields.o mem_stilt.o ModTurbKE.o micphys.o rconstants.o \
	mem_leaf.o ModTurbFields.o mem_grell.o ModTurbDiff.o grid_dims.o ModTKenn.o \
	mem_turb_scalar.o ModTurbKAdap.o ModMonotonicAdvection.o mem_cuparm.o \
	ModNamelistFile.o node_mod.o mem_scratch.o ModScalarTable.o mem_chem1.o \
	ModRGrad.o ModTurbDiffAdap.o mem_tend.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKAdap.o : $(TURB)/ModTurbKAdap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTurbKE.o : $(TURB)/ModTurbKE.f90 mem_scratch.o ke_coms.o mem_grid.o \
	rconstants.o ModTurbFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuv2.7.o : $(TUV)/ModTuv2.7.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModTuvDriver2.7.o : $(TUV)/ModTuvDriver2.7.f90 node_mod.o ModTuv2.7.o \
	tuvParameter.o mem_rrtm.o mem_carma.o mem_grid.o ModBasicFields.o mem_chem1.o \
	chem_fastjx_driv.o mem_tuv.o mem_aerad.o rconstants.o ref_sounding.o extra.o \
	mem_leaf.o mem_globrad.o mem_radiate.o chem1_list.o 
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
	ModBasicFields.o ModTurbFields.o mem_tend.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfFile.o : $(FDDA)/ModVarfFile.f90 ModDateUtils.o ModNudAnalysis.o \
	mem_grid.o aer1_list.o ModGridTree.o ModVarfUpdate.o ModMessageSet.o \
	ModBasicFields.o rconstants.o micphys.o mem_leaf.o mem_aer1.o mem_varinit.o \
	isan_coms.o parlibf.o ModControlVars.o node_mod.o mem_scratch.o ModRcio.o \
	mem_chem1.o ref_sounding.o ModGrid.o ReadBcst.o chem1_list.o \
	ModRamsReadHeader.o $(UTILS_INCS)/constants.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModWindFarm.o : $(WIND_FARM)/ModWindFarm.f90 node_mod.o ModDateUtils.o \
	io_params.o module_wind_fitch.o mem_grid.o ModBasicFields.o rconstants.o \
	ModTurbFields.o ModNamelistFile.o mem_tend.o $(UTILS_INCS)/files.h 
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

nest_feed.o : $(NESTING)/nest_feed.f90 mem_grid.o 
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

node_mod.o : $(MPI)/node_mod.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

npf.o : $(MATRIX)/npf.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModNudRead.o : $(FDDA)/ModNudRead.f90 ModDateUtils.o ModNudAnalysis.o mem_grid.o \
	mem_chem1.o mem_varinit.o ModNudUpdate.o isan_coms.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModNudUpdate.o : $(FDDA)/ModNudUpdate.f90 var_tables.o ModInitHis.o ModRcio.o \
	mem_grid.o an_header.o mem_chem1.o mem_aerad.o grid_struct.o mem_varinit.o \
	chem1_list.o $(UTILS_INCS)/files.h 
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

ModOdaRead.o : $(FDDA)/ModOdaRead.f90 ModDateUtils.o mem_grid.o isan_coms.o \
	ModOdaStaInput.o mem_oda.o ModOdaStaCount.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaCount.o : $(FDDA)/ModOdaStaCount.f90 mem_grid.o obs_input.o mem_oda.o \
	ModReadRalph.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModOdaStaInput.o : $(FDDA)/ModOdaStaInput.f90 ModDateUtils.o mem_grid.o \
	obs_input.o ModReadRalph.o mem_oda.o ModOdaStaCount.o 
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

raco_adap.o : $(MODEL)/raco_adap.f90 node_mod.o mem_scratch.o mem_grid.o \
	ModRbnd.o rconstants.o ModGrid.o ModMessageSet.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_carma.o : $(RADIATE)/rad_carma.F90 node_mod.o ModDateUtils.o grid_dims.o \
	mem_carma.o ccatt_start.o mem_grid.o mem_chem1.o mem_globaer.o mem_tuv.o \
	rconstants.o mem_aerad.o aer1_list.o mem_aer1.o mem_leaf.o mem_globrad.o \
	machine_arq.o mem_radiate.o carma_fastjx.o chem1_list.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rad_driv.o : $(RADIATE)/rad_driv.f90 mem_radiate.o ModCarmaDriver.o \
	ModRrtmDriver.o ModBasicFields.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

radvc_adap.o : $(MODEL)/radvc_adap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rams_grid.o : $(INIT)/rams_grid.f90 mem_grid.o node_mod.o dump.o rconstants.o \
	$(UTILS_INCS)/constants.h 
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

ReadBcst.o : $(MPI)/ReadBcst.f90 node_mod.o var_tables.o mem_grid.o an_header.o \
	ModBasicFields.o parlibf.o mem_aerad.o ModMPassFull.o ModControlVars.o \
	ModTurbFields.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ref_sounding.o : $(MODEL)/ref_sounding.f90 ModNamelistFile.o grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rnode.o : $(MODEL)/rnode.f90 node_mod.o var_tables.o grid_dims.o mem_grid.o \
	parlibf.o mem_leaf.o $(UTILS_INCS)/constants.h 
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

rrtmg_lw_init.o : $(RRTMG_LW_SRC)/rrtmg_lw_init.f90 parrrtm.o rrlw_con.o \
	rrlw_kg10.o parkind.o rrlw_kg12.o rrlw_wvn.o rrlw_kg03.o rrlw_kg07.o \
	rrlw_kg13.o rrlw_kg11.o rrlw_kg06.o rrlw_kg09.o rrlw_kg15.o rrlw_vsn.o \
	rrtmg_lw_setcoef.o rrlw_kg02.o rrlw_kg14.o rrlw_cld.o rrlw_tbl.o rrlw_kg05.o \
	rrlw_kg08.o rrlw_kg16.o rrlw_kg04.o rrlw_kg01.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_k_g.o : $(RRTMG_LW_SRC)/rrtmg_lw_k_g.f90 rrlw_kg03.o rrlw_kg05.o \
	rrlw_kg07.o rrlw_kg08.o rrlw_kg13.o rrlw_kg04.o rrlw_kg16.o rrlw_kg10.o \
	rrlw_kg02.o parkind.o rrlw_kg01.o rrlw_kg11.o rrlw_kg06.o rrlw_kg14.o \
	rrlw_kg12.o rrlw_kg09.o rrlw_kg15.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rad.o : $(RRTMG_LW_SRC)/rrtmg_lw_rad.f90 rrtmg_lw_cldprmc.o \
	rrtmg_lw_setcoef.o parrrtm.o rrlw_con.o parkind.o rrtmg_lw_rtrnmc.o \
	rrtmg_lw_taumol.o rrlw_wvn.o mcica_subcol_gen_lw.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrn.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrn.f90 rrlw_tbl.o parrrtm.o \
	rrlw_con.o parkind.o rrlw_wvn.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmc.f90 rrlw_tbl.o parrrtm.o \
	rrlw_con.o parkind.o rrlw_wvn.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmr.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmr.f90 rrlw_tbl.o parrrtm.o \
	rrlw_con.o parkind.o rrlw_wvn.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_setcoef.o : $(RRTMG_LW_SRC)/rrtmg_lw_setcoef.f90 parrrtm.o parkind.o \
	rrlw_ref.o rrlw_wvn.o rrlw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_taumol.o : $(RRTMG_LW_SRC)/rrtmg_lw_taumol.f90 parrrtm.o rrlw_con.o \
	rrlw_kg10.o parkind.o rrlw_ref.o rrlw_kg12.o rrlw_wvn.o rrlw_kg03.o rrlw_kg07.o \
	rrlw_kg13.o rrlw_kg11.o rrlw_kg06.o rrlw_kg09.o rrlw_kg15.o rrlw_vsn.o \
	rrlw_kg02.o rrlw_kg14.o rrlw_kg05.o rrlw_kg08.o rrlw_kg16.o rrlw_kg04.o \
	rrlw_kg01.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprmc.f90 rrsw_cld.o rrsw_wvn.o \
	parrrsw.o parkind.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprop.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprop.f90 rrsw_cld.o rrsw_wvn.o \
	parrrsw.o parkind.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_init.o : $(RRTMG_SW_SRC)/rrtmg_sw_init.f90 rrsw_kg29.o rrsw_wvn.o \
	parrrsw.o rrsw_kg22.o parkind.o rrsw_kg24.o rrsw_kg27.o rrsw_con.o rrsw_kg18.o \
	rrsw_kg21.o rrsw_tbl.o rrsw_kg19.o rrsw_kg23.o rrsw_aer.o rrsw_kg26.o \
	rrsw_vsn.o rrsw_kg20.o rrsw_kg17.o rrsw_cld.o rrtmg_sw_setcoef.o rrsw_kg25.o \
	rrsw_kg16.o rrsw_kg28.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_k_g.o : $(RRTMG_SW_SRC)/rrtmg_sw_k_g.f90 rrsw_kg24.o rrsw_kg27.o \
	rrsw_kg29.o rrsw_kg22.o rrsw_kg28.o rrsw_kg20.o parkind.o rrsw_kg21.o \
	rrsw_kg18.o rrsw_kg25.o rrsw_kg19.o rrsw_kg23.o rrsw_kg26.o rrsw_kg16.o \
	rrsw_vsn.o rrsw_kg17.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_rad.o : $(RRTMG_SW_SRC)/rrtmg_sw_rad.f90 rrsw_aer.o rrsw_wvn.o \
	parrrsw.o rrtmg_sw_spcvmc.o rrtmg_sw_setcoef.o parkind.o rrsw_con.o \
	mcica_subcol_gen_sw.o rrtmg_sw_cldprmc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_reftra.o : $(RRTMG_SW_SRC)/rrtmg_sw_reftra.f90 parkind.o rrsw_vsn.o \
	rrsw_tbl.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_setcoef.o : $(RRTMG_SW_SRC)/rrtmg_sw_setcoef.f90 parkind.o parrrsw.o \
	rrsw_vsn.o rrsw_ref.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvmc.f90 rrsw_wvn.o \
	rrtmg_sw_reftra.o parrrsw.o rrtmg_sw_taumol.o parkind.o rrsw_tbl.o \
	rrtmg_sw_vrtqdr.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvrt.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvrt.f90 rrsw_wvn.o \
	rrtmg_sw_reftra.o parrrsw.o rrtmg_sw_taumol.o parkind.o rrsw_tbl.o \
	rrtmg_sw_vrtqdr.o rrsw_vsn.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_taumol.o : $(RRTMG_SW_SRC)/rrtmg_sw_taumol.f90 rrsw_kg24.o rrsw_kg27.o \
	rrsw_kg29.o rrsw_wvn.o parrrsw.o rrsw_kg17.o rrsw_kg22.o parkind.o rrsw_con.o \
	rrsw_kg18.o rrsw_kg21.o rrsw_kg19.o rrsw_kg23.o rrsw_kg25.o rrsw_kg26.o \
	rrsw_vsn.o rrsw_kg16.o rrsw_kg20.o rrsw_kg28.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_vrtqdr.o : $(RRTMG_SW_SRC)/rrtmg_sw_vrtqdr.f90 parkind.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModRUser.o : $(SURFACE)/ModRUser.f90 node_mod.o io_params.o ccatt_start.o \
	mem_grid.o leaf_coms.o rconstants.o memSoilMoisture.o mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

setup.o : $(MATRIX)/setup.f90 memMatrix.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

sfclyr_jules.o : $(JULES_DIR)/sfclyr_jules.f90 gridmean_fluxes.o \
	jules_surface_types_mod.o mem_grid.o gridbox_mean_mod.o mem_radiate.o \
	csigma_mod.o mem_micro.o leaf_coms.o ModBasicFields.o mem_jules.o ancil_info.o \
	rconstants.o micphys.o datetime_mod.o mem_leaf.o ModTurbFields.o io_constants.o \
	sf_diags_mod.o mem_brams_jules.o io_params.o mem_cuparm.o model_time_mod.o \
	node_mod.o jules_fields_mod.o mem_chem1.o ModLeaf3Init.o fluxes.o chem1_list.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

shcu_vars_const.o : $(CUPARM)/shcu_vars_const.f90 ModNamelistFile.o grid_dims.o \
	conv_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

soilMoisture.o : $(SOIL_MOISTURE)/soilMoisture.F90 node_mod.o io_params.o \
	mem_grid.o leaf_coms.o ModBasicFields.o parlibf.o rconstants.o mem_aerad.o \
	ModMPassFull.o memSoilMoisture.o mem_leaf.o ReadBcst.o ModControlVars.o \
	ModTurbFields.o ModNamelistFile.o $(UTILS_INCS)/constants.h \
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

ModTKenn.o : $(STILT)/ModTKenn.f90 turb_constants.o mem_scratch.o mem_grid.o \
	mem_stilt.o rconstants.o 
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

var_tables.o : $(MEMORY)/var_tables.f90 ModParallelEnvironment.o \
	$(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

ModVarfUpdate.o : $(FDDA)/ModVarfUpdate.f90 mem_scratch.o ModInitHis.o \
	mem_grid.o rconstants.o ref_sounding.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

vtab_fill.o : $(MEMORY)/vtab_fill.f90 var_tables.o io_params.o mem_chem1.o \
	aer1_list.o mem_aer1.o chem1_list.o $(UTILS_INCS)/constants.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

include jules_depend_model.mk

