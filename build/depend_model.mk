actv.o : $(MATRIX)/actv.f90 memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

adap_init.o : $(INIT)/adap_init.f90 mem_leaf.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

advSendMod.o : $(ADVC)/advSendMod.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

aer1_list.o : $(AEROSOL)/aer1_list_$(AERLEVEL).f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 
	@ln -fs aer1_list_$(AERLEVEL).o aer1_list.o

alloc.o : $(MEMORY)/alloc.F90 mem_jules.o mem_globrad.o teb_spm_start.o \
	mem_teb_common.o mem_opt_scratch.o mem_aerad.o mem_shcu.o mem_globaer.o \
	mem_teb.o mem_all.o mem_gaspart.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

an_header.o : $(UTILS_MODS)/an_header.f90 $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

aobj.o : $(ISAN)/aobj.f90 isan_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

asgen.o : $(ISAN)/asgen.f90 mem_grid.o ModDateUtils.o isan_coms.o io_params.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

asti.o : $(ISAN)/asti.f90 mem_grid.o isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

asti2.o : $(ISAN)/asti2.f90 ModDateUtils.o isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

astp.o : $(ISAN)/astp.f90 isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

avarf.o : $(ISAN)/avarf.f90 mem_grid.o isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

BoundaryMod.o : $(ADVC)/BoundaryMod.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

carma_driver.o : $(RADIATE)/carma_driver.f90 mem_grid.o mem_radiate.o node_mod.o \
	mem_carma.o mem_leaf.o mem_micro.o rad_carma.o teb_spm_start.o mem_teb_common.o \
	mem_scratch1_grell.o mem_cuparm.o rconstants.o micphys.o mem_tend.o mem_basic.o \
	ModDateUtils.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

carma_fastjx.o : $(CCATT)/carma_fastjx.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ccatt_start.o : $(CCATT)/ccatt_start.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem1_list.o : $(MODEL_CHEM)/chem1_list.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem1aq_list.o : $(MODEL_CHEM)/chem1aq_list.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_aobj.o : $(ISAN_CHEM)/chem_aobj.f90 isan_coms.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_asgen.o : $(ISAN_CHEM)/chem_asgen.F90 mem_grid.o aer1_list.o node_mod.o \
	dump.o chem_isan_coms.o mem_chem1.o chem1_list.o io_params.o isan_coms.o \
	ModDateUtils.o mem_aer1.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/files.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_asti.o : $(ISAN_CHEM)/chem_asti.f90 mem_grid.o chem_isan_coms.o mem_chem1.o \
	isan_coms.o mem_aer1.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_asti2.o : $(ISAN_CHEM)/chem_asti2.f90 ModDateUtils.o isan_coms.o \
	rconstants.o chem_isan_coms.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_astp.o : $(ISAN_CHEM)/chem_astp.F90 dump.o chem_isan_coms.o mem_varinit.o \
	mem_chem1.o chem1_list.o isan_coms.o rconstants.o ModDateUtils.o \
	$(UTILS_INCS)/constants.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_avarf.o : $(ISAN_CHEM)/chem_avarf.f90 mem_grid.o chem_isan_coms.o \
	mem_chem1.o isan_coms.o rconstants.o mem_aer1.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_conv_transp.o : $(CCATT)/chem_conv_transp.f90 mem_grid.o aer1_list.o \
	node_mod.o mem_scratch.o mem_chem1.o chem1_list.o mem_scratch1_grell.o \
	Phys_const.o mem_cuparm.o mem_grell_param2.o mem_basic.o mem_tconv.o mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_dry_dep.o : $(MODEL_CHEM)/chem_dry_dep.f90 aer1_list.o mem_chem1.o \
	chem1_list.o extra.o ModDateUtils.o mem_aer1.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_fastjx57.o : $(CCATT)/chem_fastjx57.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_fastjx_data.o : $(CCATT)/chem_fastjx_data.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_fastjx_driv.o : $(CCATT)/chem_fastjx_driv.f90 rconstants.o chem_fastjx57.o \
	chem_fastjx_data.o chem1_list.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_file_inv.o : $(ISAN_CHEM)/chem_file_inv.f90 mem_grid.o ModDateUtils.o \
	dump.o isan_coms.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/files.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_first_rams.o : $(ISAN_CHEM)/chem_first_rams.f90 mem_grid.o mem_scratch.o \
	isan_coms.o an_header.o rconstants.o $(UTILS_INCS)/i8.h $(UTILS_INCS)/files.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_isan_coms.o : $(ISAN_CHEM)/chem_isan_coms.f90 aer1_list.o isan_coms.o \
	chem1_list.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_isan_io.o : $(ISAN_CHEM)/chem_isan_io.f90 isan_coms.o $(UTILS_INCS)/i8.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_orage.o : $(CCATT)/chem_orage.f90 mem_scratch1_grell.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_plumerise_scalar.o : $(CCATT)/chem_plumerise_scalar.f90 mem_plume_chem1.o \
	node_mod.o mem_chem1.o mem_aer1.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_refstate.o : $(ISAN_CHEM)/chem_refstate.f90 ccatt_start.o rconstants.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_sources.o : $(CCATT)/chem_sources.f90 ModNamelistFile.o aer1_list.o \
	mem_grid.o ReadBcst.o mem_plume_chem1.o mem_chem1.o parlibf.o io_params.o \
	mem_volc_chem1.o ModDateUtils.o mem_aer1.o $(UTILS_INCS)/ranks.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_dratedc.o : $(MODEL_CHEM)/chem_spack_dratedc.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_fexchem.o : $(MODEL_CHEM)/chem_spack_fexchem.f90 chem_spack_rates.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_fexloss.o : $(MODEL_CHEM)/chem_spack_fexloss.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_fexprod.o : $(MODEL_CHEM)/chem_spack_fexprod.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_jacdchemdc.o : $(MODEL_CHEM)/chem_spack_jacdchemdc.f90 \
	chem_spack_dratedc.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_kinetic.o : $(MODEL_CHEM)/chem_spack_kinetic.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_lu.o : $(CCATT)/chem_spack_lu.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_qssa.o : $(CCATT)/chem_spack_qssa.f90 chem_spack_rates.o mem_chem1.o \
	chem_spack_dratedc.o chem_spack_fexloss.o chem_spack_kinetic.o \
	chem_spack_fexprod.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_rates.o : $(MODEL_CHEM)/chem_spack_rates.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_rodas3_dyndt.o : $(CCATT)/chem_spack_rodas3_dyndt.f90 mem_grid.o \
	mem_spack.o mem_chem1.o extra.o chem_spack_jacdchemdc.o chem_spack_ros.o \
	chem_spack_fexchem.o chem_spack_kinetic.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_ros.o : $(CCATT)/chem_spack_ros.f90 mem_spack.o \
	chem_spack_solve_sparse.o mem_chem1.o chem_spack_jacdchemdc.o \
	chem_spack_fexchem.o chem_spack_kinetic.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_ros_dyndt.o : $(CCATT)/chem_spack_ros_dyndt.f90 mem_spack.o \
	chem_spack_solve_sparse.o mem_chem1.o chem_spack_ros.o chem_spack_jacdchemdc.o \
	chem_spack_fexchem.o chem_spack_kinetic.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_solve_sparse.o : $(CCATT)/chem_spack_solve_sparse.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_spack_utils.o : $(CCATT)/chem_spack_utils.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_trans_gasaq.o : $(MODEL_CHEM)/chem_trans_gasaq.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_trans_liq.o : $(CCATT)/chem_trans_liq.f90 mem_chem1aq.o mem_chem1.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_uv_att.o : $(CCATT)/chem_uv_att.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chem_v_interps.o : $(ISAN_CHEM)/chem_v_interps.f90 isan_coms.o dump.o \
	rconstants.o $(UTILS_INCS)/constants.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ChemDryDepDriver.o : $(MODEL_CHEM)/ChemDryDepDriver.f90 mem_grid.o mem_radiate.o \
	mem_leaf.o mem_micro.o mem_chem1.o mem_turb.o mem_cuparm.o chem_dry_dep.o \
	rconstants.o micphys.o mem_basic.o mem_aer1.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

chemistry.o : $(CCATT)/chemistry.f90 node_mod.o mem_scratch.o parrrtm.o \
	grid_dims.o extra.o chem_spack_ros.o mem_basic.o mem_rrtm.o ModTuvDriver2.7.o \
	mem_radiate.o aer1_list.o carma_fastjx.o mem_carma.o mem_micro.o \
	chem_spack_ros_dyndt.o chem_uv_att.o chem1_list.o chem_fastjx_driv.o \
	chem1aq_list.o mem_cuparm.o mem_chem1aq.o rconstants.o mem_aer1.o mem_grid.o \
	mem_spack.o chem_orage.o mem_aerad.o micphys.o mem_stilt.o chem_trans_liq.o \
	chem_spack_solve_sparse.o chem_trans_gasaq.o mem_globrad.o chem_spack_qssa.o \
	mem_chem1.o mem_scratch1_grell.o chem_spack_rodas3_dyndt.o mem_grell_param2.o \
	mem_chemic.o chem_spack_utils.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ChemSourcesDriver.o : $(CCATT)/ChemSourcesDriver.f90 aer1_list.o chem_sources.o \
	mem_leaf.o mem_plume_chem1.o mem_chem1.o chem1_list.o chem_plumerise_scalar.o \
	io_params.o mem_stilt.o mem_volc_chem1.o mem_aer1.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

coag.o : $(MATRIX)/coag.f90 setup.o memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

cond_read.o : $(FDDA)/cond_read.f90 mem_grid.o ModDateUtils.o isan_coms.o \
	mem_varinit.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cond_update.o : $(FDDA)/cond_update.f90 mem_grid.o grid_struct.o var_tables.o \
	an_header.o rconstants.o mem_varinit.o $(UTILS_INCS)/i8.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

conv_coms.o : $(CUPARM)/conv_coms.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ConvPar_GF_GEOS5.o : $(CUPARM)/ConvPar_GF_GEOS5.F90 module_gate.o \
	MAPL_Constants.o Henrys_Law_cts.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

coriolis.o : $(MODEL)/coriolis.f90 mem_grid.o ref_sounding.o node_mod.o \
	mem_scratch.o parlibf.o rconstants.o mem_tend.o mem_basic.o ModBuffering.o \
	$(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cu_read.o : $(CUPARM)/cu_read.f90 mem_grid.o mem_cuparm.o isan_coms.o \
	mem_basic.o ModDateUtils.o $(UTILS_INCS)/i8.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cup_dn.o : $(CUPARM)/cup_dn.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cup_env.o : $(CUPARM)/cup_env.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cup_env_catt.o : $(CUPARM)/cup_env_catt.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cup_grell3.o : $(CUPARM)/cup_grell3.F90 ModNamelistFile.o module_cu_g3.o \
	node_mod.o mem_grell.o ModMessageSet.o mem_scratch.o Phys_const.o mem_basic.o \
	mem_radiate.o mem_carma.o mem_leaf.o mem_micro.o mem_cuparm.o var_tables.o \
	rconstants.o ccatt_start.o mem_grid.o ModGrid.o io_params.o micphys.o \
	mem_stilt.o module_cu_gf.o mem_jules.o module_cu_gf_v5.1.o mem_chem1.o \
	mem_turb.o mem_scratch1_grell.o ConvPar_GF_GEOS5.o mem_grell_param2.o \
	mem_tend.o mem_varinit.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cup_grell_catt_deep.o : $(CUPARM)/cup_grell_catt_deep.f90 mem_grid.o node_mod.o \
	mem_carma.o Phys_const.o mem_scratch3_grell.o cup_output_vars.o ccatt_start.o \
	mem_grell_param2.o kbcon_ecmwf.o mem_varinit.o mem_scratch2_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cup_grell_catt_shallow.o : $(CUPARM)/cup_grell_catt_shallow.f90 mem_grid.o \
	mem_scratch3_grell_sh.o node_mod.o mem_scratch2_grell_sh.o Phys_const.o \
	cup_output_vars.o mem_grell_param2.o mem_varinit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cup_output_vars.o : $(CUPARM)/cup_output_vars.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

cup_up.o : $(CUPARM)/cup_up.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

dam.o : $(ENERGY)/dam.f90 ModNamelistFile.o ModDateUtils.o dump.o mem_grid.o \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

debugTools.o : $(UTILS_TOOLS)/debugTools.f90 mem_grid.o node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

depv.o : $(MATRIX)/depv.f90 memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

diam.o : $(MATRIX)/diam.f90 memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

dicrete.o : $(MATRIX)/dicrete.f90 coag.o memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

diffsclr.o : $(TURB)/diffsclr.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

diffuse.o : $(TURB)/diffuse.f90 mem_grid.o node_mod.o mem_leaf.o mem_micro.o \
	mem_scratch.o mem_opt_scratch.o mem_turb.o ke_coms.o var_tables.o micphys.o \
	mem_tend.o mem_basic.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

digitalFilter.o : $(MODEL)/digitalFilter.f90 ModNamelistFile.o mem_grid.o \
	node_mod.o ReadBcst.o io_params.o grid_dims.o var_tables.o mem_basic.o \
	ModDateUtils.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

domain_decomp.o : $(INIT)/domain_decomp.f90 ModNamelistFile.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

errorMod.o : $(ADVC)/errorMod.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

extra.o : $(MEMORY)/extra.f90 ModNamelistFile.o var_tables.o dump.o \
	$(UTILS_INCS)/i8.h $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

file_inv.o : $(ISAN)/file_inv.f90 mem_grid.o ModDateUtils.o isan_coms.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

first_rams.o : $(ISAN)/first_rams.f90 mem_grid.o mem_scratch.o isan_coms.o \
	an_header.o rconstants.o $(UTILS_INCS)/i8.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

gaspart.o : $(TEB_SPM)/gaspart.f90 mem_grid.o mem_teb_vars_const.o node_mod.o \
	mem_leaf.o parlibf.o var_tables.o an_header.o mem_emiss.o mem_basic.o \
	mem_gaspart.o $(UTILS_INCS)/i8.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

geodat.o : $(MKSFC)/geodat.f90 mem_grid.o io_params.o mem_leaf.o teb_spm_start.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

gfdl_cloud_microphys.o : $(MICRO)/gfdl_cloud_microphys.F90 node_mod.o \
	module_mp_radar.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

grid_dims.o : $(MEMORY)/grid_dims.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

grid_struct.o : $(MEMORY)/grid_struct.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

GridMod.o : $(ADVC)/GridMod.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

gridset.o : $(INIT)/gridset.f90 mem_grid.o rconstants.o grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

Henrys_Law_cts.o : $(CUPARM)/Henrys_Law_cts.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

init_mtx.o : $(MATRIX)/init_mtx.f90 setup.o dicrete.o memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

InitAdvect.o : $(ADVC)/InitAdvect.f90 ProcessorMod.o GridMod.o dump.o MapMod.o \
	BoundaryMod.o advSendMod.o errorMod.o $(UTILS_INCS)/constants.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

initComm.o : $(COMM_SPC)/initComm.f90 mem_grid.o parlibf.o ReadBcst.o \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

inithis.o : $(IO)/inithis.f90 mem_grid.o ref_sounding.o mem_leaf.o leaf_coms.o \
	mem_scratch.o mem_chem1.o io_params.o chem1_list.o mem_aerad.o var_tables.o \
	an_header.o rconstants.o micphys.o mem_basic.o mem_varinit.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

initMicThompson.o : $(MICRO)/initMicThompson.f90 mem_grid.o dump.o node_mod.o \
	ReadBcst.o mem_micro.o parlibf.o generic.o mem_basic.o ModDateUtils.o \
	$(UTILS_INCS)/constants.f90 $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

io_params.o : $(IO)/io_params.f90 ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

isan_coms.o : $(ISAN_MODS)/isan_coms.f90 ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

isan_io.o : $(ISAN)/isan_io.f90 isan_coms.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

isofwd.o : $(MATRIX)/isofwd.f90 solut.o solut.o solut.o issoropia.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

isorev.o : $(MATRIX)/isorev.f90 solut.o issoropia.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

isrpia.o : $(MATRIX)/isrpia.f90 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

issoropia.o : $(MATRIX)/issoropia.f90 isrpia.o solut.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

kbcon_ecmwf.o : $(CUPARM)/kbcon_ecmwf.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ke_coms.o : $(TURB)/ke_coms.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

landuse_input.o : $(MKSFC)/landuse_input.f90 mem_mksfc.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

leaf3.o : $(SURFACE)/leaf3.f90 mem_grid.o mem_radiate.o node_mod.o mem_leaf.o \
	teb_spm_start.o mem_teb_common.o mem_micro.o leaf_coms.o mem_scratch.o \
	mem_turb.o io_params.o mem_cuparm.o rconstants.o micphys.o mem_teb.o \
	mem_basic.o ccatt_start.o mem_all.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

leaf3_hyd.o : $(SURFACE)/leaf3_hyd.f90 mem_grid.o mem_leaf.o leaf_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

leaf3_init.o : $(SURFACE)/leaf3_init.f90 mem_grid.o mem_leaf.o leaf_coms.o \
	teb_spm_start.o io_params.o grid_dims.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

leaf3_ocean_only.o : $(SURFACE)/leaf3_ocean_only.f90 mem_grid.o mem_radiate.o \
	node_mod.o mem_leaf.o mem_micro.o leaf_coms.o mem_turb.o io_params.o \
	ConvPar_GF_GEOS5.o mem_cuparm.o rconstants.o micphys.o mem_basic.o \
	ccatt_start.o mem_all.o cup_grell3.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

leaf3_teb.o : $(SURFACE)/leaf3_teb.f90 mem_teb_vars_const.o mem_emiss.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

leaf_coms.o : $(SURFACE)/leaf_coms.f90 grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

local_proc.o : $(MODEL)/local_proc.F90 mem_grid.o ref_sounding.o node_mod.o \
	dump.o ReadBcst.o io_params.o grid_dims.o rconstants.o mem_stilt.o \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

machine_arq.o : $(MODEL)/machine_arq.F90 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(CMACH) $(<F:.F90=.F90)
	@rm -f $(<F:.f90=.f90) 

MAPL_Constants.o : $(CUPARM)/MAPL_Constants.F90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

MapMod.o : $(ADVC)/MapMod.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

matrix.o : $(MATRIX)/matrix.f90 subs.o depv.o actv.o setup.o npf.o memMatrix.o \
	coag.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

MatrixDriver.o : $(MATRIX)/MatrixDriver.F90 mem_grid.o aer1_list.o subs.o \
	node_mod.o mem_radiate.o mem_leaf.o mem_micro.o setup.o mem_chem1.o \
	chem1_list.o mem_turb.o memMatrix.o npf.o mem_aer1.o rconstants.o ModParticle.o \
	micphys.o mem_basic.o isrpia.o coag.o 
	cp -f  $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER) $(<F:.F90=.F90)
	rm -f $(<F:.F90=.F90)

mcica_random_numbers.o : $(RRTMG_SW_SRC)/mcica_random_numbers.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_lw.o : $(RRTMG_LW_SRC)/mcica_subcol_gen_lw.f90 parkind.o \
	mcica_random_numbers.o parrrtm.o rrlw_vsn.o rrlw_con.o rrlw_wvn.o \
	mcica_random_numbers.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mcica_subcol_gen_sw.o : $(RRTMG_SW_SRC)/mcica_subcol_gen_sw.f90 parkind.o \
	rrsw_con.o mcica_random_numbers.o rrsw_vsn.o parrrsw.o rrsw_wvn.o \
	mcica_random_numbers.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_aer1.o : $(CCATT)/mem_aer1.f90 ModNamelistFile.o aer1_list.o mem_grid.o \
	dump.o node_mod.o mem_chem1.o io_params.o grid_dims.o var_tables.o \
	$(UTILS_INCS)/i8.h $(UTILS_INCS)/constants.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_aerad.o : $(RADIATE)/mem_aerad.f90 mem_grid_dim_defs.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_all.o : $(MEMORY)/mem_all.f90 mem_grid.o mem_radiate.o ref_sounding.o \
	mem_leaf.o mem_micro.o mem_scratch.o mem_basic.o mem_turb.o mem_scalar.o \
	io_params.o mem_cuparm.o var_tables.o mem_oda.o mem_scratch1_brams.o mem_tend.o \
	mem_nestb.o micphys.o mem_varinit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_basic.o : $(MEMORY)/mem_basic.f90 ModNamelistFile.o var_tables.o mem_grid.o \
	mem_stilt.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_carma.o : $(RADIATE)/mem_carma.f90 ModNamelistFile.o mem_grid.o node_mod.o \
	mem_globrad.o ReadBcst.o parlibf.o io_params.o mem_scalar.o grid_dims.o \
	mem_aerad.o var_tables.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_chem1.o : $(CCATT)/mem_chem1.f90 ModNamelistFile.o io_params.o chem1_list.o \
	grid_dims.o var_tables.o $(UTILS_INCS)/i8.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_chem1aq.o : $(CCATT)/mem_chem1aq.f90 ModNamelistFile.o mem_chem1.o \
	chem1aq_list.o grid_dims.o var_tables.o $(UTILS_INCS)/i8.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_chemic.o : $(CCATT)/mem_chemic.f90 micphys.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_cuparm.o : $(CUPARM)/mem_cuparm.f90 ModNamelistFile.o var_tables.o \
	grid_dims.o $(UTILS_INCS)/i8.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_cutrans.o : $(CUPARM)/mem_cutrans.f90 mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_emiss.o : $(TEB_SPM)/mem_emiss.f90 ModNamelistFile.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_gaspart.o : $(TEB_SPM)/mem_gaspart.f90 var_tables.o mem_emiss.o \
	$(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_globaer.o : $(RADIATE)/mem_globaer.f90 mem_precision.o mem_aerad.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_globrad.o : $(RADIATE)/mem_globrad.f90 ModNamelistFile.o mem_radiate.o \
	node_mod.o parlibf.o mem_aerad.o mem_precision.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_grell.o : $(CUPARM)/mem_grell.f90 mem_cuparm.o var_tables.o \
	shcu_vars_const.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_grell_param2.o : $(CUPARM)/mem_grell_param2.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_grid.o : $(MEMORY)/mem_grid.f90 ModNamelistFile.o var_tables.o grid_dims.o \
	$(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_grid_dim_defs.o : $(MEMORY)/mem_grid_dim_defs.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_jules.o : $(SURFACE)/mem_jules.f90 ModNamelistFile.o var_tables.o \
	io_params.o grid_dims.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_leaf.o : $(SURFACE)/mem_leaf.f90 ModNamelistFile.o teb_spm_start.o \
	io_params.o grid_dims.o var_tables.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_micro.o : $(MICRO)/mem_micro.f90 mem_cuparm.o mem_radiate.o var_tables.o \
	micphys.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_micro_optij.o : $(MICRO)/mem_micro_optij.f90 mem_grid.o mem_radiate.o \
	node_mod.o mem_micro.o grid_dims.o rconstants.o micphys.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_mksfc.o : $(MKSFC)/mem_mksfc.f90 teb_spm_start.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_nestb.o : $(NESTING)/mem_nestb.f90 var_tables.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_oda.o : $(FDDA)/mem_oda.f90 ModNamelistFile.o var_tables.o grid_dims.o \
	$(UTILS_INCS)/i8.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_opt_scratch.o : $(TURB)/mem_opt_scratch.f90 $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_plume_chem1.o : $(CCATT)/mem_plume_chem1.f90 ModNamelistFile.o var_tables.o \
	mem_chem1.o chem1_list.o $(UTILS_INCS)/i8.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_precision.o : $(RADIATE)/mem_precision.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_radiate.o : $(RADIATE)/mem_radiate.f90 ModNamelistFile.o var_tables.o \
	$(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_rrtm.o : $(RADIATE)/mem_rrtm.f90 mem_grid.o node_mod.o parkind.o mem_chem1.o \
	chem1_list.o parrrtm.o rrtmg_lw_init.o parrrsw.o rconstants.o rrtmg_sw_init.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

mem_scalar.o : $(MEMORY)/mem_scalar.f90 ModNamelistFile.o var_tables.o \
	io_params.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_scratch.o : $(MEMORY)/mem_scratch.f90 mem_aerad.o mem_radiate.o node_mod.o \
	grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_scratch1_brams.o : $(MEMORY)/mem_scratch1_brams.f90 var_tables.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_scratch1_grell.o : $(CUPARM)/mem_scratch1_grell.f90 ccatt_start.o dump.o \
	mem_grell_param2.o $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_scratch2_grell.o : $(CUPARM)/mem_scratch2_grell.f90 mem_grell_param2.o \
	node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_scratch2_grell_sh.o : $(CUPARM)/mem_scratch2_grell_sh.f90 mem_grell_param2.o \
	node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_scratch3_grell.o : $(CUPARM)/mem_scratch3_grell.f90 mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_scratch3_grell_sh.o : $(CUPARM)/mem_scratch3_grell_sh.f90 mem_grell_param2.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_shcu.o : $(CUPARM)/mem_shcu.f90 var_tables.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_spack.o : $(CCATT)/mem_spack.f90 chem1_list.o chem_spack_utils.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_stilt.o : $(STILT)/mem_stilt.f90 ModNamelistFile.o io_params.o grid_dims.o \
	var_tables.o rconstants.o $(UTILS_INCS)/i8.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_tconv.o : $(CCATT)/mem_tconv.f90 aer1_list.o chem1_list.o mem_aer1.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_teb.o : $(TEB_SPM)/mem_teb.f90 var_tables.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_teb_common.o : $(TEB_SPM)/mem_teb_common.f90 var_tables.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_teb_vars_const.o : $(TEB_SPM)/mem_teb_vars_const.f90 ModNamelistFile.o \
	grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_tend.o : $(MEMORY)/mem_tend.f90 mem_grid.o ModNamelistFile.o teb_spm_start.o \
	mem_micro.o mem_turb.o mem_scalar.o mem_emiss.o mem_basic.o mem_gaspart.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_turb.o : $(TURB)/mem_turb.f90 ModNamelistFile.o var_tables.o mem_stilt.o \
	grid_dims.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_turb_scalar.o : $(TURB)/mem_turb_scalar.f90 var_tables.o grid_dims.o \
	$(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_tuv.o : $(TUV)/mem_tuv.f90 var_tables.o mem_globrad.o ModTuv2.7.o \
	mem_stilt.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND)  $(<F:.f90=.f90)
	@rm -f $(<F:.f90=.f90) 

mem_varinit.o : $(MEMORY)/mem_varinit.f90 ModNamelistFile.o mem_chem1.o \
	chem1_list.o grid_dims.o var_tables.o $(UTILS_INCS)/i8.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mem_volc_chem1.o : $(CCATT)/mem_volc_chem1.f90 ModNamelistFile.o var_tables.o \
	$(UTILS_INCS)/i8.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

memMatrix.o : $(MATRIX)/memMatrix.f90 ModNamelistFile.o aer1_list.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

memSoilMoisture.o : $(SOIL_MOISTURE)/memSoilMoisture.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

meteogram.o : $(IO)/meteogram.f90 ModNamelistFile.o mem_grid.o node_mod.o \
	ModPostUtils.o var_tables.o satPolyColision.o meteogramType.o \
	$(POST_INCS)/post_rconstants.h $(UTILS_INCS)/files.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

meteogramType.o : $(IO)/meteogramType.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_coll.o : $(MICRO)/mic_coll.f90 rconstants.o micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_driv.o : $(MICRO)/mic_driv.f90 mem_grid.o mem_radiate.o node_mod.o \
	mem_micro.o mem_chem1.o mem_chem1aq.o micphys.o mem_chemic.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_driv_new.o : $(MICRO)/mic_driv_new.f90 mem_grid.o mem_radiate.o node_mod.o \
	mem_micro_optij.o mem_micro.o mem_chem1.o mem_chem1aq.o micphys.o mem_chemic.o \
	mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_gamma.o : $(MICRO)/mic_gamma.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_gfdl_driver.o : $(MICRO)/mic_gfdl_driver.f90 mem_grid.o mem_radiate.o \
	node_mod.o mem_leaf.o mem_micro.o gfdl_cloud_microphys.o io_params.o \
	rconstants.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_init.o : $(MICRO)/mic_init.f90 mem_grid.o dump.o rconstants.o micphys.o \
	$(UTILS_INCS)/constants.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_misc.o : $(MICRO)/mic_misc.f90 mem_grid.o mem_micro.o mem_scratch.o \
	rconstants.o micphys.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_nuc.o : $(MICRO)/mic_nuc.f90 rconstants.o micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_tabs.o : $(MICRO)/mic_tabs.f90 micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_thompson_driver.o : $(MICRO)/mic_thompson_driver.f90 mem_grid.o \
	mem_radiate.o node_mod.o module_mp_thompson.o mem_leaf.o mem_micro.o \
	io_params.o rconstants.o micphys.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mic_vap.o : $(MICRO)/mic_vap.f90 rconstants.o micphys.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

micphys.o : $(MICRO)/micphys.f90 ModNamelistFile.o grid_dims.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mksfc_driver.o : $(MKSFC)/mksfc_driver.f90 mem_grid.o node_mod.o ReadBcst.o \
	teb_spm_start.o io_params.o grid_dims.o mem_mksfc.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mksfc_fuso.o : $(MKSFC)/mksfc_fuso.f90 mem_grid.o mem_teb_vars_const.o \
	node_mod.o ReadBcst.o io_params.o mem_emiss.o mem_teb.o mem_mksfc.o \
	mem_gaspart.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mksfc_ndvi.o : $(MKSFC)/mksfc_ndvi.f90 mem_grid.o dump.o mem_leaf.o io_params.o \
	mem_mksfc.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mksfc_sfc.o : $(MKSFC)/mksfc_sfc.f90 mem_grid.o dump.o node_mod.o mem_leaf.o \
	ReadBcst.o io_params.o mem_mksfc.o $(UTILS_INCS)/constants.f90 \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mksfc_sst.o : $(MKSFC)/mksfc_sst.f90 mem_grid.o mem_mksfc.o io_params.o \
	mem_leaf.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mksfc_top.o : $(MKSFC)/mksfc_top.f90 mem_grid.o dump.o node_mod.o ReadBcst.o \
	io_params.o mem_mksfc.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mod_advect_kit.o : $(MODEL)/mod_advect_kit.f90 mem_grid.o node_mod.o \
	var_tables.o mod_GhostBlock.o mem_tend.o mem_basic.o mod_GhostBlockPartition.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mod_aer.o : $(AERCLIM)/mod_aer.f90 mem_grid.o dump.o node_mod.o ReadBcst.o \
	parlibf.o mem_basic.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mod_GhostBlock.o : $(MODEL)/mod_GhostBlock.f90 mod_GhostBlockPartition.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mod_GhostBlockPartition.o : $(MODEL)/mod_GhostBlockPartition.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mod_ozone.o : $(TEB_SPM)/mod_ozone.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModBramsGrid.o : $(POST_SRC)/ModBramsGrid.f90 mem_grid.o ModNamelistFile.o \
	ref_sounding.o node_mod.o ModPostUtils.o mem_aerad.o ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModBuffering.o : $(MPI)/ModBuffering.f90 ModParallelEnvironment.o 
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

ModDateUtils.o : $(UTILS_MODS)/ModDateUtils.f90 dump.o \
	$(UTILS_INCS)/constants.f90 $(UTILS_INCS)/ranks.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModDomainDecomp.o : $(MPI)/ModDomainDecomp.f90 ModGridDims.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModEvaluation.o : $(EVAL)/ModEvaluation.f90 mem_grid.o parlibf.o node_mod.o \
	ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModFieldSectionList.o : $(MPI)/ModFieldSectionList.f90 var_tables.o \
	ModParallelEnvironment.o ModDomainDecomp.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModGrid.o : $(MPI)/ModGrid.f90 ModNamelistFile.o ModMessageSet.o ModGridDims.o \
	ModDomainDecomp.o ModNeighbourNodes.o ModParallelEnvironment.o var_tables.o \
	meteogramType.o ModMessagePassing.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModGridDims.o : $(MPI)/ModGridDims.f90 ModNamelistFile.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModGridTree.o : $(MPI)/ModGridTree.f90 ModNamelistFile.o ModGrid.o \
	ModParallelEnvironment.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

modIau.o : $(MODEL)/modIau.f90 ModNamelistFile.o mem_grid.o dump.o node_mod.o \
	ReadBcst.o parlibf.o mem_tend.o mem_varinit.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/constants.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModMessageData.o : $(MPI)/ModMessageData.f90 ModDomainDecomp.o \
	ModNeighbourNodes.o ModParallelEnvironment.o var_tables.o ModFieldSectionList.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModMessagePassing.o : $(MPI)/ModMessagePassing.f90 ModNamelistFile.o mem_grid.o \
	ModMessageSet.o ModDomainDecomp.o ModNeighbourNodes.o ModParallelEnvironment.o \
	var_tables.o ModGridDims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModMessageSet.o : $(MPI)/ModMessageSet.f90 parlibf.o ModDomainDecomp.o \
	ModNeighbourNodes.o ModParallelEnvironment.o var_tables.o ModBuffering.o \
	ModMessageData.o ModFieldSectionList.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModNamelistFile.o : $(INIT)/ModNamelistFile.f90 dump.o parlibf.o grid_dims.o \
	ModParallelEnvironment.o modPrintInitial.o $(UTILS_INCS)/constants.f90 \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModNeighbourNodes.o : $(MPI)/ModNeighbourNodes.f90 ModGridDims.o \
	ModParallelEnvironment.o ModDomainDecomp.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModOneProc.o : $(MODEL)/ModOneProc.F90 ModNamelistFile.o node_mod.o \
	soilMoisture.o ModEvaluation.o initComm.o mem_teb_common.o mem_scratch.o \
	parlibf.o grid_dims.o ModVarfFile.o extra.o mem_basic.o ModTuvDriver2.7.o \
	mem_radiate.o radvc_mnt.o aer1_list.o module_wind_farm.o mem_leaf.o mem_micro.o \
	mem_carma.o mem_plume_chem1.o chem1_list.o mem_scalar.o tuvParameter.o \
	mem_cuparm.o var_tables.o isan_coms.o rtimh_rk.o mem_chem1aq.o \
	ModParallelEnvironment.o ModGridTree.o mem_volc_chem1.o \
	module_rams_microphysics_2M.o ccatt_start.o mem_aer1.o cup_grell3.o mem_grid.o \
	mod_aer.o dump.o ModPostProcess.o rtimh_abm.o ReadBcst.o teb_spm_start.o \
	ModGrid.o digitalFilter.o rtimh.o io_params.o ModTuv2.7.o dam.o chem_dry_dep.o \
	memSoilMoisture.o mem_emiss.o domain_decomp.o micphys.o mem_teb.o \
	ModPostGridNetCDF.o mem_stilt.o shcu_vars_const.o InitAdvect.o meteogram.o \
	mem_gaspart.o mem_teb_vars_const.o initMicThompson.o ref_sounding.o \
	mem_globrad.o chem_sources.o machine_arq.o modIau.o mem_chem1.o mem_turb.o \
	mod_advect_kit.o ModTimeStamp.o mem_grell_param2.o mem_oda.o local_proc.o \
	mem_varinit.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/files.h \
	$(UTILS_INCS)/i8.h $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModOutputUtils.o : $(IO)/ModOutputUtils.f90 mem_basic.o var_tables.o dump.o \
	mem_turb.o $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModParallelEnvironment.o : $(MPI)/ModParallelEnvironment.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModParticle.o : $(MATRIX)/ModParticle.f90 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

ModPostGrid.o : $(POST_SRC)/ModPostGrid.F90 ModNamelistFile.o mem_grid.o \
	ModOutputUtils.o ModPostTypes.o ModPostUtils.o parlibf.o io_params.o \
	ModBramsGrid.o ModParallelEnvironment.o ModPostOneFieldNetCDF.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostGridNetCDF.o : $(POST_SRC)/ModPostGridNetCDF.F90 ModNamelistFile.o \
	mem_grid.o ModPostTypes.o ModPostUtils.o dump.o io_params.o ModBramsGrid.o \
	ModDateUtils.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/ranks.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostOneField.o : $(POST_SRC)/ModPostOneField.f90 ModNamelistFile.o \
	ModPostTypes.o ModPostUtils.o ModPostOneField7d.o ModPostOneFieldUtils.o dump.o \
	node_mod.o ModPostOneField8d.o ModBramsGrid.o ModPostOneField3d.o \
	ModPostOneField2d.o $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostOneField2d.o : $(POST_SRC)/ModPostOneField2d.f90 mem_grid.o mem_radiate.o \
	ModOutputUtils.o ModPostTypes.o ModPostGrid.o ModPostOneFieldUtils.o \
	ModPostUtils.o dump.o node_mod.o io_params.o ModBramsGrid.o mem_aerad.o \
	mem_cuparm.o micphys.o $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostOneField3d.o : $(POST_SRC)/ModPostOneField3d.f90 mem_grid.o \
	ModOutputUtils.o ModPostTypes.o node_mod.o ModPostGrid.o ModPostOneFieldUtils.o \
	ModPostUtils.o ModBramsGrid.o micphys.o mem_varinit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostOneField7d.o : $(POST_SRC)/ModPostOneField7d.f90 ModOutputUtils.o \
	ModPostTypes.o ModPostUtils.o ModPostOneFieldUtils.o ModBramsGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostOneField8d.o : $(POST_SRC)/ModPostOneField8d.f90 ModOutputUtils.o \
	ModPostTypes.o ModPostUtils.o ModPostOneFieldUtils.o ModBramsGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostOneFieldNetCDF.o : $(POST_SRC)/ModPostOneFieldNetCDF.F90 \
	ModPostGridNetCDF.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostOneFieldUtils.o : $(POST_SRC)/ModPostOneFieldUtils.f90 ModOutputUtils.o \
	ModPostTypes.o ModPostGrid.o ModBramsGrid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostProcess.o : $(POST_SRC)/ModPostProcess.F90 ModNamelistFile.o \
	ModPostTypes.o ModPostGrid.o ModMessageSet.o ModGrid.o io_params.o \
	ModBramsGrid.o ModParallelEnvironment.o ModGridTree.o ModPostGridNetCDF.o \
	ModPostOneField.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostTypes.o : $(POST_SRC)/ModPostTypes.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModPostUtils.o : $(POST_SRC)/ModPostUtils.f90 ModParallelEnvironment.o dump.o \
	mem_leaf.o $(POST_INCS)/post_rconfig.h $(POST_INCS)/post_rconstants.h \
	$(UTILS_INCS)/constants.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

modPrintInitial.o : $(INIT)/modPrintInitial.F90 dump.o \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

modsched.o : $(MODEL)/modsched.f90 mem_grid.o mem_radiate.o ref_sounding.o \
	dump.o node_mod.o ReadBcst.o mem_varinit.o parlibf.o io_params.o mem_turb.o \
	mem_cuparm.o isan_coms.o mem_basic.o local_proc.o shcu_vars_const.o \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModTimeStamp.o : $(MODEL)/ModTimeStamp.f90 $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModTuv2.7.o : $(TUV)/ModTuv2.7.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModTuvDriver2.7.o : $(TUV)/ModTuvDriver2.7.f90 mem_grid.o mem_radiate.o \
	ref_sounding.o node_mod.o mem_globrad.o mem_carma.o mem_leaf.o mem_tuv.o \
	mem_chem1.o ModTuv2.7.o chem1_list.o chem_fastjx_driv.o mem_aerad.o \
	tuvParameter.o extra.o mem_all.o rconstants.o mem_basic.o mem_rrtm.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_cu_g3.o : $(CUPARM)/module_cu_g3.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_cu_gd_fim.o : $(CUPARM)/module_cu_gd_fim.f90 module_gate.o Phys_const.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_cu_gf.o : $(CUPARM)/module_cu_gf.f90 node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_cu_gf_v5.1.o : $(CUPARM)/module_cu_gf_v5.1.f90 module_gate.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_mp_radar.o : $(MICRO)/module_mp_radar.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_mp_thompson.o : $(MICRO)/module_mp_thompson.f90 module_mp_radar.o \
	node_mod.o $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_rams_microphysics_2M.o : $(MICRO)/module_rams_microphysics_2M.f90 \
	mem_grid.o dump.o node_mod.o mem_leaf.o mem_micro.o mem_scratch.o grid_dims.o \
	rconstants.o micphys.o mem_basic.o $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_wind_farm.o : $(WIND_FARM)/module_wind_farm.f90 ModNamelistFile.o \
	mem_grid.o module_wind_fitch.o node_mod.o io_params.o mem_turb.o rconstants.o \
	mem_tend.o mem_basic.o ModDateUtils.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_wind_fitch.o : $(WIND_FARM)/module_wind_fitch.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ModVarfFile.o : $(FDDA)/ModVarfFile.f90 mem_grid.o aer1_list.o ref_sounding.o \
	node_mod.o ModDateUtils.o ModMessageSet.o ReadBcst.o mem_leaf.o mem_scratch.o \
	ModGrid.o parlibf.o mem_chem1.o chem1_list.o isan_coms.o rconstants.o \
	ModGridTree.o micphys.o mem_basic.o mem_varinit.o mem_aer1.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mpass_dtl.o : $(MPI)/mpass_dtl.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mpass_feed.o : $(MPI)/mpass_feed.f90 mem_grid.o node_mod.o parlibf.o grid_dims.o \
	var_tables.o mem_scratch1_brams.o mem_basic.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mpass_full.o : $(MPI)/mpass_full.f90 mem_grid.o node_mod.o mem_turb.o \
	io_params.o mem_aerad.o var_tables.o an_header.o mem_basic.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mpass_nest.o : $(MPI)/mpass_nest.f90 mem_grid.o node_mod.o dump.o mem_nestb.o \
	parlibf.o grid_dims.o var_tables.o mem_basic.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

mpi_io_engine-5d.o : $(IO)/mpi_io_engine-5d.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ndvi_read.o : $(MKSFC)/ndvi_read.f90 mem_grid.o node_mod.o mem_leaf.o ReadBcst.o \
	io_params.o ModDateUtils.o mem_mksfc.o $(UTILS_INCS)/i8.h $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

nest_drivers.o : $(NESTING)/nest_drivers.f90 mem_grid.o node_mod.o mem_nestb.o \
	mem_scratch.o var_tables.o mem_tend.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

nest_feed.o : $(NESTING)/nest_feed.f90 mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

nest_filldens.o : $(NESTING)/nest_filldens.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

nest_geosst.o : $(MKSFC)/nest_geosst.f90 mem_grid.o dump.o node_mod.o \
	soilMoisture.o mem_leaf.o mem_scratch.o io_params.o memSoilMoisture.o \
	mem_basic.o ccatt_start.o mem_mksfc.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

nest_intrp.o : $(NESTING)/nest_intrp.f90 mem_grid.o ref_sounding.o mem_nestb.o \
	mem_scratch.o grid_dims.o rconstants.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

node_mod.o : $(MPI)/node_mod.f90 ModNamelistFile.o mem_grid.o ModDomainDecomp.o \
	grid_dims.o ModParallelEnvironment.o ModGridTree.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

npf.o : $(MATRIX)/npf.f90 memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

nud_analysis.o : $(FDDA)/nud_analysis.f90 mem_grid.o node_mod.o dump.o \
	ModEvaluation.o mem_scratch.o modIau.o mem_chem1.o chem1_list.o mem_tend.o \
	mem_basic.o mem_varinit.o $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

nud_read.o : $(FDDA)/nud_read.f90 mem_grid.o mem_varinit.o mem_chem1.o \
	isan_coms.o ModDateUtils.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

nud_update.o : $(FDDA)/nud_update.f90 mem_grid.o mem_chem1.o grid_struct.o \
	chem1_list.o mem_aerad.o var_tables.o an_header.o rconstants.o mem_basic.o \
	mem_varinit.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

obs_input.o : $(FDDA)/obs_input.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

oda_krig.o : $(FDDA)/oda_krig.f90 mem_oda.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

oda_nudge.o : $(FDDA)/oda_nudge.f90 mem_grid.o node_mod.o mem_scratch.o \
	io_params.o mem_oda.o mem_tend.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

oda_proc_obs.o : $(FDDA)/oda_proc_obs.f90 mem_grid.o mem_oda.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

oda_read.o : $(FDDA)/oda_read.f90 mem_grid.o mem_oda.o ModDateUtils.o \
	isan_coms.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

oda_sta_count.o : $(FDDA)/oda_sta_count.f90 mem_grid.o mem_oda.o obs_input.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

oda_sta_input.o : $(FDDA)/oda_sta_input.f90 mem_grid.o mem_oda.o ModDateUtils.o \
	obs_input.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

opspec.o : $(IO)/opspec.f90 grid_dims.o mem_radiate.o aer1_list.o mem_leaf.o \
	chem1_list.o chem1aq_list.o mem_cuparm.o mem_chem1aq.o ccatt_start.o mem_aer1.o \
	mem_grid.o teb_spm_start.o io_params.o mem_emiss.o micphys.o mem_stilt.o \
	shcu_vars_const.o mem_globrad.o chem_sources.o modIau.o mem_chem1.o mem_turb.o \
	mem_grell_param2.o mem_varinit.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

optical.o : $(RADIATE)/optical.f90 mem_grid.o aer1_list.o mem_radiate.o dump.o \
	node_mod.o mem_leaf.o ReadBcst.o parlibf.o var_tables.o mem_basic.o \
	ccatt_start.o mem_aer1.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ozone.o : $(TEB_SPM)/ozone.f90 mem_grid.o mem_radiate.o mod_ozone.o rconstants.o \
	mem_basic.o mem_gaspart.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

para_init.o : $(MPI)/para_init.f90 mem_grid.o node_mod.o dump.o grid_dims.o \
	var_tables.o domain_decomp.o mem_basic.o $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

parkind.o : $(RRTMG_SW_MOD)/parkind.f90 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

parrrsw.o : $(RRTMG_SW_MOD)/parrrsw.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

parrrtm.o : $(RRTMG_LW_MOD)/parrrtm.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

Phys_const.o : $(CUPARM)/Phys_const.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ProcessorMod.o : $(ADVC)/ProcessorMod.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

quad.o : $(MATRIX)/quad.f90 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

raco.o : $(MODEL)/raco.f90 mem_grid.o ModNamelistFile.o ref_sounding.o \
	node_mod.o dump.o ModMessageSet.o initComm.o mem_scratch.o ModGrid.o \
	ModParallelEnvironment.o rconstants.o micphys.o mem_tend.o mem_basic.o \
	raco_adap.o $(UTILS_INCS)/constants.f90 $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF) 
	@rm -f $(<F:.f90=.f90) 

raco_adap.o : $(MODEL)/raco_adap.f90 mem_grid.o node_mod.o ModMessageSet.o \
	mem_scratch.o ModGrid.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rad_carma.o : $(RADIATE)/rad_carma.F90 mem_grid.o mem_radiate.o aer1_list.o \
	mem_globrad.o carma_fastjx.o mem_carma.o mem_leaf.o node_mod.o mem_scratch.o \
	mem_tuv.o machine_arq.o mem_chem1.o chem1_list.o mem_aerad.o ccatt_start.o \
	mem_globaer.o rconstants.o ModDateUtils.o mem_aer1.o \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(AER)  $(<F:.F90=.F90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rad_driv.o : $(RADIATE)/rad_driv.f90 mem_radiate.o rtm_driver.o carma_driver.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

radvc.o : $(MODEL)/radvc.f90 mem_grid.o ModNamelistFile.o radvc_mnt.o node_mod.o \
	mem_scratch.o mem_chem1.o grid_dims.o var_tables.o chem_dry_dep.o mem_tend.o \
	mem_basic.o ccatt_start.o mem_aer1.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

radvc_adap.o : $(MODEL)/radvc_adap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

radvc_mnt.o : $(MODEL)/radvc_mnt.f90 mem_grid.o ModNamelistFile.o node_mod.o \
	mem_scratch.o parlibf.o mem_chem1.o var_tables.o chem_dry_dep.o rconstants.o \
	micphys.o mem_aer1.o mem_basic.o ccatt_start.o advSendMod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

radvc_new.o : $(MODEL)/radvc_new.f90 node_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

radvc_rk.o : $(MODEL)/radvc_rk.f90 mem_grid.o node_mod.o initComm.o mem_chem1.o \
	grid_dims.o var_tables.o mem_stilt.o mem_tend.o mem_basic.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rams_grid.o : $(INIT)/rams_grid.f90 mem_grid.o dump.o node_mod.o rconstants.o \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rams_mem_alloc.o : $(MEMORY)/rams_mem_alloc.F90 mem_scratch3_grell_sh.o \
	node_mod.o mem_grell.o ModEvaluation.o mem_scratch2_grell_sh.o mem_scratch.o \
	mem_teb_common.o mem_opt_scratch.o mem_tuv.o grid_dims.o extra.o parrrsw.o \
	mem_shcu.o mem_basic.o mem_radiate.o aer1_list.o carma_fastjx.o mem_carma.o \
	mem_leaf.o mem_micro.o mem_chemic.o mem_plume_chem1.o mem_turb_scalar.o \
	mem_scalar.o chem1_list.o mem_cuparm.o var_tables.o chem1aq_list.o \
	mem_chem1aq.o mem_volc_chem1.o ccatt_start.o optical.o mem_aer1.o cup_grell3.o \
	mem_grid.o teb_spm_start.o digitalFilter.o io_params.o ModTuv2.7.o mem_aerad.o \
	chem_dry_dep.o mem_emiss.o mem_stilt.o mem_teb.o mem_nestb.o micphys.o \
	shcu_vars_const.o mem_gaspart.o mem_grid_dim_defs.o mem_teb_vars_const.o \
	mem_jules.o mem_globrad.o mem_micro_optij.o chem_sources.o machine_arq.o \
	modIau.o mem_chem1.o mem_turb.o mem_scratch1_grell.o mem_scratch3_grell.o \
	mem_oda.o mem_grell_param2.o mem_globaer.o mem_scratch1_brams.o mem_tend.o \
	mem_varinit.o mem_scratch2_grell.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rams_read_header.o : $(IO)/rams_read_header.f90 an_header.o \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ranlavg.o : $(IO)/ranlavg.f90 mem_grid.o var_tables.o io_params.o grid_dims.o \
	$(UTILS_INCS)/i8.h $(UTILS_INCS)/interface.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rbnd.o : $(BC)/rbnd.f90 mem_grid.o ref_sounding.o node_mod.o mem_scratch.o \
	mem_chem1.o mem_turb.o var_tables.o micphys.o mem_tend.o mem_basic.o \
	ccatt_start.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rbnd_adap.o : $(BC)/rbnd_adap.f90 mem_grid.o ref_sounding.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rcio.o : $(IO)/rcio.f90 mem_all.o mem_stilt.o leaf_coms.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rconstants.o : $(MEMORY)/rconstants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rconv.o : $(CUPARM)/rconv.f90 mem_grid.o node_mod.o conv_coms.o mem_scratch.o \
	mem_cuparm.o rconstants.o mem_tend.o mem_basic.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rconv_grell_catt.o : $(CUPARM)/rconv_grell_catt.f90 mem_grid.o mem_radiate.o \
	node_mod.o mem_grell.o mem_leaf.o mem_micro.o mem_scratch.o mem_turb.o \
	mem_scalar.o io_params.o mem_cuparm.o mem_scratch1_grell.o mem_grell_param2.o \
	rconstants.o micphys.o mem_tend.o mem_basic.o mem_stilt.o ccatt_start.o \
	$(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

read_ralph.o : $(FDDA)/read_ralph.f90 obs_input.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ReadBcst.o : $(MPI)/ReadBcst.f90 mem_grid.o node_mod.o parlibf.o mem_turb.o \
	mem_aerad.o var_tables.o an_header.o mem_basic.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

recycle.o : $(IO)/recycle.f90 mem_grid.o aer1_list.o dump.o node_mod.o \
	ReadBcst.o parlibf.o mem_chem1.o io_params.o chem1_list.o mem_aerad.o \
	var_tables.o an_header.o ModDateUtils.o mem_aer1.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/constants.f90 $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ref_sounding.o : $(MODEL)/ref_sounding.f90 ModNamelistFile.o grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

refstate.o : $(ISAN)/refstate.f90 rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rexev.o : $(STILT)/rexev.f90 mem_grid.o mem_micro.o mem_scratch.o rconstants.o \
	mem_stilt.o mem_tend.o mem_basic.o micphys.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rgrad.o : $(TURB)/rgrad.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rhhi.o : $(INIT)/rhhi.f90 mem_grid.o ref_sounding.o mem_scratch.o rconstants.o \
	micphys.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rinit.o : $(INIT)/rinit.f90 mem_grid.o ref_sounding.o node_mod.o mem_micro.o \
	mem_scratch.o io_params.o rconstants.o micphys.o mem_basic.o mem_varinit.o \
	$(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rio.o : $(IO)/rio.f90 mem_grid.o ref_sounding.o node_mod.o ReadBcst.o parlibf.o \
	mem_chem1.o io_params.o mem_turb.o grid_dims.o mem_aerad.o var_tables.o \
	mpi_io_engine-5d.o an_header.o mem_basic.o ModDateUtils.o $(UTILS_INCS)/i8.h \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/interface.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rnest_par.o : $(MPI)/rnest_par.f90 mem_grid.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rnode.o : $(MODEL)/rnode.f90 mem_grid.o node_mod.o mem_leaf.o parlibf.o \
	grid_dims.o var_tables.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rrlw_cld.o : $(RRTMG_LW_MOD)/rrlw_cld.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_con.o : $(RRTMG_LW_MOD)/rrlw_con.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg01.o : $(RRTMG_LW_MOD)/rrlw_kg01.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg02.o : $(RRTMG_LW_MOD)/rrlw_kg02.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg03.o : $(RRTMG_LW_MOD)/rrlw_kg03.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg04.o : $(RRTMG_LW_MOD)/rrlw_kg04.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg05.o : $(RRTMG_LW_MOD)/rrlw_kg05.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg06.o : $(RRTMG_LW_MOD)/rrlw_kg06.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg07.o : $(RRTMG_LW_MOD)/rrlw_kg07.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg08.o : $(RRTMG_LW_MOD)/rrlw_kg08.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg09.o : $(RRTMG_LW_MOD)/rrlw_kg09.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg10.o : $(RRTMG_LW_MOD)/rrlw_kg10.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg11.o : $(RRTMG_LW_MOD)/rrlw_kg11.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg12.o : $(RRTMG_LW_MOD)/rrlw_kg12.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg13.o : $(RRTMG_LW_MOD)/rrlw_kg13.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg14.o : $(RRTMG_LW_MOD)/rrlw_kg14.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg15.o : $(RRTMG_LW_MOD)/rrlw_kg15.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_kg16.o : $(RRTMG_LW_MOD)/rrlw_kg16.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_ncpar.o : $(RRTMG_LW_MOD)/rrlw_ncpar.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_ref.o : $(RRTMG_LW_MOD)/rrlw_ref.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_tbl.o : $(RRTMG_LW_MOD)/rrlw_tbl.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_vsn.o : $(RRTMG_LW_MOD)/rrlw_vsn.f90 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrlw_wvn.o : $(RRTMG_LW_MOD)/rrlw_wvn.f90 parkind.o parrrtm.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_aer.o : $(RRTMG_SW_MOD)/rrsw_aer.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_cld.o : $(RRTMG_SW_MOD)/rrsw_cld.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_con.o : $(RRTMG_SW_MOD)/rrsw_con.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg16.o : $(RRTMG_SW_MOD)/rrsw_kg16.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg17.o : $(RRTMG_SW_MOD)/rrsw_kg17.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg18.o : $(RRTMG_SW_MOD)/rrsw_kg18.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg19.o : $(RRTMG_SW_MOD)/rrsw_kg19.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg20.o : $(RRTMG_SW_MOD)/rrsw_kg20.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg21.o : $(RRTMG_SW_MOD)/rrsw_kg21.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg22.o : $(RRTMG_SW_MOD)/rrsw_kg22.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg23.o : $(RRTMG_SW_MOD)/rrsw_kg23.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg24.o : $(RRTMG_SW_MOD)/rrsw_kg24.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg25.o : $(RRTMG_SW_MOD)/rrsw_kg25.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg26.o : $(RRTMG_SW_MOD)/rrsw_kg26.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg27.o : $(RRTMG_SW_MOD)/rrsw_kg27.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg28.o : $(RRTMG_SW_MOD)/rrsw_kg28.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_kg29.o : $(RRTMG_SW_MOD)/rrsw_kg29.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_ref.o : $(RRTMG_SW_MOD)/rrsw_ref.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_tbl.o : $(RRTMG_SW_MOD)/rrsw_tbl.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_vsn.o : $(RRTMG_SW_MOD)/rrsw_vsn.f90 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrsw_wvn.o : $(RRTMG_SW_MOD)/rrsw_wvn.f90 parrrsw.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprmc.f90 rrlw_cld.o parkind.o \
	parrrtm.o rrlw_vsn.o rrlw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_cldprop.o : $(RRTMG_LW_SRC)/rrtmg_lw_cldprop.f90 rrlw_vsn.o rrlw_cld.o \
	parrrtm.o parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_init.o : $(RRTMG_LW_SRC)/rrtmg_lw_init.f90 rrlw_kg01.o rrlw_kg14.o \
	rrlw_kg06.o rrlw_kg16.o rrtmg_lw_setcoef.o parrrtm.o rrlw_vsn.o rrlw_con.o \
	rrlw_kg15.o rrlw_kg04.o rrlw_kg08.o rrlw_kg10.o parkind.o rrlw_kg09.o \
	rrlw_kg05.o rrlw_kg12.o rrlw_tbl.o rrlw_kg02.o rrlw_wvn.o rrlw_cld.o \
	rrlw_kg03.o rrlw_kg11.o rrlw_kg13.o rrlw_kg07.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_k_g.o : $(RRTMG_LW_SRC)/rrtmg_lw_k_g.f90 rrlw_kg01.o rrlw_kg06.o \
	rrlw_kg14.o rrlw_kg15.o rrlw_kg10.o rrlw_kg16.o parkind.o rrlw_kg09.o \
	rrlw_kg04.o rrlw_kg05.o rrlw_kg03.o rrlw_kg12.o rrlw_vsn.o rrlw_kg11.o \
	rrlw_kg08.o rrlw_kg13.o rrlw_kg07.o rrlw_kg02.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rad.o : $(RRTMG_LW_SRC)/rrtmg_lw_rad.f90 rrtmg_lw_rtrnmc.o \
	rrtmg_lw_setcoef.o parkind.o mcica_subcol_gen_lw.o rrtmg_lw_taumol.o parrrtm.o \
	rrlw_vsn.o rrtmg_lw_cldprmc.o rrlw_con.o rrlw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrn.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrn.f90 parkind.o parrrtm.o \
	rrlw_vsn.o rrlw_tbl.o rrlw_con.o rrlw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmc.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmc.f90 parkind.o parrrtm.o \
	rrlw_vsn.o rrlw_tbl.o rrlw_con.o rrlw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_rtrnmr.o : $(RRTMG_LW_SRC)/rrtmg_lw_rtrnmr.f90 parkind.o parrrtm.o \
	rrlw_vsn.o rrlw_tbl.o rrlw_con.o rrlw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_setcoef.o : $(RRTMG_LW_SRC)/rrtmg_lw_setcoef.f90 parkind.o parrrtm.o \
	rrlw_vsn.o rrlw_wvn.o rrlw_ref.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_lw_taumol.o : $(RRTMG_LW_SRC)/rrtmg_lw_taumol.f90 rrlw_kg01.o rrlw_kg14.o \
	rrlw_kg06.o rrlw_kg16.o parrrtm.o rrlw_vsn.o rrlw_con.o rrlw_kg15.o rrlw_kg04.o \
	rrlw_kg08.o rrlw_kg10.o parkind.o rrlw_kg09.o rrlw_kg05.o rrlw_kg12.o \
	rrlw_kg02.o rrlw_wvn.o rrlw_kg03.o rrlw_kg11.o rrlw_kg13.o rrlw_kg07.o \
	rrlw_ref.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprmc.f90 rrsw_cld.o parkind.o \
	rrsw_vsn.o parrrsw.o rrsw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_cldprop.o : $(RRTMG_SW_SRC)/rrtmg_sw_cldprop.f90 rrsw_cld.o parkind.o \
	rrsw_vsn.o parrrsw.o rrsw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_init.o : $(RRTMG_SW_SRC)/rrtmg_sw_init.f90 rrsw_cld.o rrsw_kg25.o \
	rrsw_kg17.o parrrsw.o rrsw_kg19.o rrsw_kg29.o rrsw_aer.o rrsw_kg20.o \
	rrsw_kg22.o rrsw_con.o rrsw_kg21.o rrsw_tbl.o parkind.o rrsw_kg28.o rrsw_kg18.o \
	rrsw_wvn.o rrsw_kg16.o rrsw_kg26.o rrsw_kg24.o rrsw_vsn.o rrtmg_sw_setcoef.o \
	rrsw_kg23.o rrsw_kg27.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_k_g.o : $(RRTMG_SW_SRC)/rrtmg_sw_k_g.f90 rrsw_kg22.o rrsw_kg16.o \
	parkind.o rrsw_kg25.o rrsw_kg18.o rrsw_kg26.o rrsw_kg21.o rrsw_kg24.o \
	rrsw_kg17.o rrsw_vsn.o rrsw_kg28.o rrsw_kg19.o rrsw_kg29.o rrsw_kg23.o \
	rrsw_kg27.o rrsw_kg20.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_rad.o : $(RRTMG_SW_SRC)/rrtmg_sw_rad.f90 parkind.o rrsw_con.o \
	mcica_subcol_gen_sw.o rrsw_vsn.o parrrsw.o rrtmg_sw_setcoef.o \
	rrtmg_sw_cldprmc.o rrtmg_sw_spcvmc.o rrsw_aer.o rrsw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_reftra.o : $(RRTMG_SW_SRC)/rrtmg_sw_reftra.f90 rrsw_tbl.o parkind.o \
	rrsw_vsn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_setcoef.o : $(RRTMG_SW_SRC)/rrtmg_sw_setcoef.f90 rrsw_ref.o parrrsw.o \
	parkind.o rrsw_vsn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvmc.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvmc.f90 parkind.o \
	rrtmg_sw_taumol.o rrsw_vsn.o rrtmg_sw_vrtqdr.o parrrsw.o rrtmg_sw_reftra.o \
	rrsw_tbl.o rrsw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_spcvrt.o : $(RRTMG_SW_SRC)/rrtmg_sw_spcvrt.f90 parkind.o \
	rrtmg_sw_taumol.o rrsw_vsn.o rrtmg_sw_vrtqdr.o parrrsw.o rrtmg_sw_reftra.o \
	rrsw_tbl.o rrsw_wvn.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_taumol.o : $(RRTMG_SW_SRC)/rrtmg_sw_taumol.f90 rrsw_kg22.o rrsw_kg16.o \
	parkind.o rrsw_con.o rrsw_kg25.o rrsw_kg18.o rrsw_kg26.o rrsw_kg21.o \
	rrsw_kg24.o rrsw_vsn.o rrsw_kg17.o rrsw_kg28.o parrrsw.o rrsw_kg19.o \
	rrsw_kg29.o rrsw_kg23.o rrsw_kg27.o rrsw_wvn.o rrsw_kg20.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rrtmg_sw_vrtqdr.o : $(RRTMG_SW_SRC)/rrtmg_sw_vrtqdr.f90 parkind.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	rm -f $(<F:.f90=.f90)

rshcupar.o : $(CUPARM)/rshcupar.f90 mem_grid.o node_mod.o conv_coms.o \
	mem_micro.o mem_scratch.o mem_turb.o mem_shcu.o mem_tend.o mem_basic.o \
	shcu_vars_const.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rstilt.o : $(STILT)/rstilt.f90 mem_grid.o radvc_mnt.o mem_scratch.o mem_turb.o \
	mem_scratch1_grell.o mem_cuparm.o mem_stilt.o mem_basic.o $(UTILS_INCS)/i8.h 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rthrm.o : $(MODEL)/rthrm.f90 mem_grid.o mem_micro.o mem_scratch.o rconstants.o \
	micphys.o mem_basic.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rtimh.o : $(MODEL)/rtimh.F90 ModNamelistFile.o node_mod.o ModMessageSet.o \
	ChemDryDepDriver.o chemistry.o mem_basic.o rtm_driver.o mem_radiate.o \
	radvc_mnt.o module_wind_farm.o mem_leaf.o mem_plume_chem1.o mem_scalar.o \
	mem_cuparm.o rconstants.o ChemSourcesDriver.o rad_driv.o \
	module_rams_microphysics_2M.o ccatt_start.o optical.o mem_aer1.o cup_grell3.o \
	mem_grid.o teb_spm_start.o ModGrid.o digitalFilter.o mem_emiss.o micphys.o \
	mem_stilt.o shcu_vars_const.o MatrixDriver.o raco.o chem_sources.o \
	machine_arq.o mod_advect_kit.o mem_turb.o mem_chem1.o ModTimeStamp.o mem_oda.o \
	mem_tend.o mem_varinit.o $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rtimh_abm.o : $(MODEL)/rtimh_abm.F90 ModNamelistFile.o node_mod.o \
	ModMessageSet.o ChemDryDepDriver.o chemistry.o mem_basic.o rtm_driver.o \
	mem_radiate.o radvc_mnt.o module_wind_farm.o mem_leaf.o mem_plume_chem1.o \
	mem_scalar.o mem_cuparm.o rconstants.o ChemSourcesDriver.o rad_driv.o \
	module_rams_microphysics_2M.o ccatt_start.o optical.o mem_aer1.o cup_grell3.o \
	mem_grid.o teb_spm_start.o ModGrid.o digitalFilter.o mem_emiss.o micphys.o \
	mem_stilt.o shcu_vars_const.o MatrixDriver.o raco.o chem_sources.o \
	machine_arq.o mod_advect_kit.o mem_turb.o mem_chem1.o ModTimeStamp.o mem_oda.o \
	mem_tend.o mem_varinit.o $(UTILS_INCS)/i8.h $(UTILS_INCS)/tsNames.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rtimh_rk.o : $(MODEL)/rtimh_rk.F90 ModNamelistFile.o node_mod.o ModMessageSet.o \
	initComm.o ChemDryDepDriver.o mem_scratch.o chemistry.o mem_basic.o \
	rtm_driver.o mem_radiate.o radvc_mnt.o module_wind_farm.o mem_leaf.o \
	mem_plume_chem1.o mem_scalar.o mem_cuparm.o rconstants.o ChemSourcesDriver.o \
	rad_driv.o module_rams_microphysics_2M.o ccatt_start.o optical.o mem_aer1.o \
	cup_grell3.o mem_grid.o mod_aer.o teb_spm_start.o ModGrid.o digitalFilter.o \
	leaf3_ocean_only.o mem_emiss.o micphys.o mem_stilt.o shcu_vars_const.o \
	MatrixDriver.o raco.o chem_sources.o machine_arq.o modIau.o mod_advect_kit.o \
	mem_turb.o mem_chem1.o ModTimeStamp.o mem_oda.o mem_tend.o mem_varinit.o \
	$(UTILS_INCS)/constants.f90 $(UTILS_INCS)/tsNames.h $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rtimi.o : $(MODEL)/rtimi.f90 mem_grid.o node_mod.o mem_grell.o mem_scratch.o \
	mem_cuparm.o var_tables.o mem_tend.o mem_basic.o shcu_vars_const.o \
	$(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

rtm_driver.o : $(RADIATE)/rtm_driver.f90 node_mod.o mcica_subcol_gen_lw.o \
	mem_tuv.o parrrtm.o parrrsw.o mem_basic.o mem_rrtm.o mem_radiate.o \
	rrtmg_sw_cldprop.o mem_carma.o mem_leaf.o mem_micro.o leaf_coms.o mem_cuparm.o \
	rconstants.o ccatt_start.o optical.o ModDateUtils.o mem_grid.o \
	rrtmg_lw_cldprop.o parkind.o teb_spm_start.o mcica_subcol_gen_sw.o micphys.o \
	ref_sounding.o rrtmg_sw_rad.o mem_chem1.o mem_scratch1_grell.o rrtmg_lw_rad.o \
	mem_grell_param2.o mem_tend.o $(UTILS_INCS)/aerosol_setup.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

ruser.o : $(SURFACE)/ruser.f90 mem_grid.o node_mod.o mem_leaf.o leaf_coms.o \
	io_params.o memSoilMoisture.o rconstants.o ccatt_start.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

seasalt.o : $(CCATT)/seasalt.f90 mem_grid.o aer1_list.o node_mod.o mem_leaf.o \
	mem_chem1.o io_params.o mem_basic.o ccatt_start.o mem_aer1.o mod_aer.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

setup.o : $(MATRIX)/setup.f90 memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

sfclyr_jules.o : $(JULES_DIR)/sfclyr_jules.f90 node_mod.o datetime_mod.o \
	io_constants.o mem_basic.o model_time_mod.o jules_fields_mod.o mem_radiate.o \
	mem_leaf.o leaf_coms.o mem_micro.o csigma_mod.o ancil_info.o chem1_list.o \
	mem_cuparm.o rconstants.o mem_grid.o gridbox_mean_mod.o io_params.o \
	gridmean_fluxes.o micphys.o mem_brams_jules.o sf_diags_mod.o mem_jules.o \
	mem_chem1.o mem_turb.o fluxes.o jules_surface_types_mod.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

shcu_vars_const.o : $(CUPARM)/shcu_vars_const.f90 ModNamelistFile.o conv_coms.o \
	grid_dims.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

soilMoisture.o : $(SOIL_MOISTURE)/soilMoisture.F90 ModNamelistFile.o mem_grid.o \
	dump.o node_mod.o mem_leaf.o leaf_coms.o ReadBcst.o parlibf.o io_params.o \
	mem_aerad.o memSoilMoisture.o rconstants.o $(UTILS_INCS)/constants.f90 \
	$(UTILS_INCS)/files.h $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

solut.o : $(MATRIX)/solut.f90 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

sst_read.o : $(MKSFC)/sst_read.f90 mem_grid.o node_mod.o mem_leaf.o ReadBcst.o \
	io_params.o grid_dims.o ModDateUtils.o mem_mksfc.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

subs.o : $(MATRIX)/subs.f90 setup.o memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

teb_spm_start.o : $(TEB_SPM)/teb_spm_start.f90 ModNamelistFile.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

thermo_isorr.o : $(MATRIX)/thermo_isorr.f90 issoropia.o memMatrix.o 
	cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90)
	rm -f $(<F:.f90=.f90)

tkenn.o : $(STILT)/tkenn.f90 mem_grid.o mem_scratch.o rconstants.o mem_stilt.o \
	turb_constants.o 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

turb_constants.o : $(STILT)/turb_constants.f90 
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

turb_diff.o : $(TURB)/turb_diff.f90 mem_grid.o mem_scratch.o mem_opt_scratch.o \
	mem_turb.o mem_cuparm.o var_tables.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

turb_diff_adap.o : $(TURB)/turb_diff_adap.f90 mem_grid.o mem_scratch.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

turb_k.o : $(TURB)/turb_k.f90 mem_grid.o radvc_mnt.o node_mod.o mem_grell.o \
	mem_leaf.o mem_micro.o mem_scratch.o mem_chem1.o mem_turb.o mem_turb_scalar.o \
	ke_coms.o var_tables.o mem_cuparm.o rconstants.o micphys.o mem_tend.o \
	mem_basic.o mem_stilt.o ccatt_start.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

turb_k_adap.o : $(TURB)/turb_k_adap.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

turb_ke.o : $(TURB)/turb_ke.f90 mem_grid.o mem_scratch.o mem_turb.o ke_coms.o \
	rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

tuvParameter.o : $(TUV)/tuvParameter.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND)  $(<F:.f90=.f90)
	@rm -f $(<F:.f90=.f90) 

upcase.o : $(CUPARM)/upcase.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

urban.o : $(SURFACE)/urban.f90 mem_teb_vars_const.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

urban_canopy.o : $(SURFACE)/urban_canopy.f90 mem_grid.o node_mod.o mem_turb.o \
	mem_tend.o mem_basic.o $(UTILS_INCS)/files.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

v_interps.o : $(ISAN)/v_interps.f90 isan_coms.o rconstants.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

var_tables.o : $(MEMORY)/var_tables.f90 $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

varf_update.o : $(FDDA)/varf_update.f90 mem_grid.o ref_sounding.o mem_scratch.o \
	rconstants.o $(UTILS_INCS)/constants.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

vtab_fill.o : $(MEMORY)/vtab_fill.f90 aer1_list.o mem_chem1.o io_params.o \
	chem1_list.o var_tables.o mem_aer1.o $(UTILS_INCS)/i8.h 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

include jules_depend_model.mk

