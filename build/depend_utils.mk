## Make object rules ##

an_header.o  : $(UTILS_MODS)/an_header.f90 $(UTILS_INCS)/i8.h
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

ModDateUtils.o  : $(UTILS_MODS)/ModDateUtils.f90 dump.o $(UTILS_INCS)/ranks.h
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

charutils.o  : $(UTILS_LIB)/charutils.f90
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

dateutils.o  : $(UTILS_LIB)/dateutils.f90
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

error_mess.o : $(UTILS_LIB)/error_mess.f90
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

getvar.o     : $(UTILS_LIB)/getvar.f90  an_header.o dump.o
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

interp_lib.o : $(UTILS_LIB)/interp_lib.f90
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

hdf5_utils.o  : $(UTILS_LIB)/hdf5_utils.F90
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

htint-opt.o  : $(UTILS_LIB)/htint-opt.f90
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

map_proj.o   : $(UTILS_LIB)/map_proj.f90
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

utilsMod.o : $(UTILS_LIB)/utilsMod.f90 dump.o 
	@cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

ModMemory.o   : $(UTILS_LIB)/ModMemory.f90
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

numutils.o   : $(UTILS_LIB)/numutils.f90 $(UTILS_INCS)/i8.h
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

polarst.o    : $(UTILS_LIB)/polarst.f90
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

therm_lib.o  : $(UTILS_LIB)/therm_lib.f90
	 @cp -f $< $(<F:.f90=.f90)
#	 $(F_COMMAND) -pi $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

utils_f.o    : $(UTILS_LIB)/utils_f.f90 ModDateUtils.o dump.o
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

vformat.o    : $(UTILS_LIB)/vformat.f90 dump.o
	 @cp -f $< $(<F:.f90=.f90)
	 $(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	 @rm -f $(<F:.f90=.f90) 

filelist.o   : $(UTILS_LIB)/filelist.F90
	  @cp -f $< $(<F:.F90=.F90)
	  $(F_COMMAND) -D$(CMACH) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	  @rm -f $(<F:.f90=.f90) 

rsys.o       : $(UTILS_LIB)/rsys.F90
	  @cp -f $< $(<F:.F90=.F90)
	  $(F_COMMAND) -D$(CMACH) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	  @rm -f $(<F:.f90=.f90) 

parlibf.o : $(UTILS_LIB)/parlibf.F90 dump.o
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

satPolyColision.o : $(UTILS_LIB)/satPolyColision.f90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

dump.o : $(UTILS_DUMP)/dump.F90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

module_gate.o : $(UTILS_MDTOOLS)/module_gate.f90
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

dted.o       : $(UTILS_LIB)/dted.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

parlib.o     : $(UTILS_LIB)/parlib.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

tmpname.o    : $(UTILS_LIB)/tmpname.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

proc.o    : $(UTILS_LIB)/proc.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

utils_c.o    : $(UTILS_LIB)/utils_c.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

eenviron.o   : $(EFF)/eenviron.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

spAllocate.o : $(UTILS_LIB)/spAllocate.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

spBuild.o : $(UTILS_LIB)/spBuild.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

spFactor.o : $(UTILS_LIB)/spFactor.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

spOutput.o : $(UTILS_LIB)/spOutput.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

spSolve.o : $(UTILS_LIB)/spSolve.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

spUtils.o : $(UTILS_LIB)/spUtils.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

spFortran.o : $(UTILS_LIB)/spFortran.c
	  $(C_COMMAND) $< $(EXTRAFLAGSC)

gammaFunction.o:  $(UTILS_LIB)/gammaFunction.f90
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) -D$(CMACH) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

generic.o:  $(UTILS_LIB)/generic.f90
	@cp -f $< $(<F:.F90=.F90)
	$(F_COMMAND) $(<F:.F90=.F90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90) 

