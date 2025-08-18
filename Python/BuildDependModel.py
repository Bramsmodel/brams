#!/usr/bin/env python3
import os

from CommonDefs import homeDir, ObjSrcActionAtDependModel, FindFortranFiles, \
    ExpandDic, ModuleName, UseSet, IncSet, RemoveModuleNameFromUseSet, \
    NewEntryAtDepend, BuildDicObjName, BuildDicModName

############################
# MAIN SCRIPT
############################


basePath=homeDir+"/"
# dictionary DicObjName
# indexed by obj names to build at depend_model (first entry at each triplet of ObjSrcActionAtDependModel)
# each entry is a list with:
#   src name without path expansion
#   full src path
#   list of module names declared at src file
#   set of module names used by src file
#   set of include names included at src file
#   action to be performed at Makefile
DicObjName=BuildDicObjName(ObjSrcActionAtDependModel, ExpandDic)
print("There are",len(DicObjName),"entries @ DicObjName (dictionary of object files at depend_model)")



# dictionary DicIncName
# indexed by include name at source file
# each entry contains the corresponding name at depend_model
DicIncName={
    "interface.h" : "$(UTILS_INCS)/interface.h",
    "usevfm.h" : "$(UTILS_INCS)/UseVfm.h",
    "post_rconfig.h" : "$(POST_INCS)/post_rconfig.h",
    "aerosol_setup.f90" : "$(UTILS_INCS)/aerosol_setup.f90",
    "mpif.h" : None,
    "tsnames.h" : "$(UTILS_INCS)/tsNames.h",
    "ranks.h" : "$(UTILS_INCS)/ranks.h",
    "constants.f90" : "$(UTILS_INCS)/constants.f90",
    "i8.h" : "$(UTILS_INCS)/i8.h",
    "netcdf.inc" : None,
    "files.h" : "$(UTILS_INCS)/files.h",
    "constants.h" : "$(UTILS_INCS)/constants.h",
    "micconstants.h" : "$(MICRO)/MicConstants.h",
    "post_rconstants.h" : "$(POST_INCS)/post_rconstants.h" }



# dictionary DicModName
# indexed by module name
# each entry contains the corresponding obj name
DicModName=BuildDicModName(DicObjName)
print("There are",len(DicModName),"entries @ DicModName (dictionary of module names defined by source files at depend_model)")



# include module names defined by software tools previously compiled
DicModName["iso_fortran_env"] = None
DicModName["netcdf"] = None
DicModName["hdf5_utils"] = None
DicModName["wgrib2mpi"] = None
DicModName["wgrib2api"] = None


# complete DicModName and DicObjName with files that will not be compiled
# include jules object model_time_mod, cited at depend_model but not target
modName="model_time_mod"
objName="model_time_mod.o"
srcNameToExpand="$(JULES_03)"
epilogSrcName="/var/model_time_mod.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object csigma, cited at depend_model but not target
modName="csigma"
objName="csigma_mod.o"
srcNameToExpand="$(JULES_21)"
epilogSrcName="/csigma_mod.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include utils object genericfunctions, cited at depend_model but not target
modName="genericfunctions"
objName="generic.o"
srcNameToExpand="$(UTILS_LIB)"
epilogSrcName="/generic.f90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object fluxes, cited at depend_model but not target
modName="fluxes"
objName="fluxes.o"
srcNameToExpand="$(JULES_02)"
epilogSrcName="/fluxes.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object gridbox_mean_mod, cited at depend_model but not target
modName="gridbox_mean_mod"
objName="gridbox_mean_mod.o"
srcNameToExpand="$(JULES_30)"
epilogSrcName="/gridbox_mean_mod.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include utils object module_gate, cited at depend_model but not target
modName="module_gate"
objName="module_gate.o"
srcNameToExpand="$(UTILS_MDTOOLS)"
epilogSrcName="/module_gate.f90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object ancil_info, cited at depend_model but not target
modName="ancil_info"
objName="ancil_info.o"
srcNameToExpand="$(JULES_02)"
epilogSrcName="/ancil_info.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object jules_fields_mod, cited at depend_model but not target
modName="jules_fields_mod"
objName="jules_fields_mod.o"
srcNameToExpand="$(JULES_03)"
epilogSrcName="/jules_fields_mod.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object sf_diags_mod, cited at depend_model but not target
modName="sf_diags_mod"
objName="sf_diags_mod.o"
srcNameToExpand="$(JULES_26)"
epilogSrcName="/sf_diags_mod.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object datetime_mod, cited at depend_model but not target
modName="datetime_mod"
objName="datetime_mod.o"
srcNameToExpand="$(JULES_29)"
epilogSrcName="/datetime_mod.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include utils object dump, cited at depend_model but not target
modName="dump"
objName="dump.o"
srcNameToExpand="$(UTILS_DUMP)"
epilogSrcName="/dump.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object io_constants, cited at depend_model but not target
modName="io_constants"
objName="io_constants.o"
srcNameToExpand="$(JULES_18)"
epilogSrcName="/io_constants.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object jules_surface_types_mod, cited at depend_model but not target
modName="jules_surface_types_mod"
objName="jules_surface_types_mod.o"
srcNameToExpand="$(JULES_02)"
epilogSrcName="/jules_surface_types_mod.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object gridmean_fluxes, cited at depend_model but not target
modName="gridmean_fluxes"
objName="gridmean_fluxes.o"
srcNameToExpand="$(JULES_03)"
epilogSrcName="/var/gridmean_fluxes.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include utils object parlib, cited at depend_model but not target
modName="parlib"
objName="parlibf.o"
srcNameToExpand="$(UTILS_LIB)"
epilogSrcName="/parlibf.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include jules object mem_brams_jules, cited at depend_model but not target
modName="mem_brams_jules"
objName="mem_brams_jules.o"
srcNameToExpand="$(JULES_27)"
epilogSrcName="/mem_brams_jules.F90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include utils object utilsmod, cited at depend_model but not target
modName="utilsmod"
objName="utilsMod.o"
srcNameToExpand="$(UTILS_LIB)"
epilogSrcName="/utilsMod.f90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]

# include utils object satpolycolision, cited at depend_model but not target
modName="satpolycolision"
objName="satPolyColision.o"
srcNameToExpand="$(UTILS_LIB)"
epilogSrcName="/satPolyColision.f90"
expandedSrcName=basePath+ExpandDic[srcNameToExpand]+epilogSrcName
srcName=srcNameToExpand+epilogSrcName
DicModName[modName]=objName
DicObjName[objName]=[srcName,expandedSrcName,[],set(),set(),None]



# if present at useSet, remove own module name
# to avoid circular dependencies
for modName,objName in DicModName.items():
    if objName == None:
        continue
    useSet=DicObjName[objName][3]
    if useSet == None:
        continue
    useSet=RemoveModuleNameFromUseSet(useSet,modName)
    DicObjName[objName][3]=useSet


    
# verify if all module names @ use set are known
OutSet=set()
for key, value in DicObjName.items():
    useSet=value[3]
    for modName in useSet:
        if modName not in DicModName.keys():
            OutSet.add(modName)
            print("module **"+modName+"** used by file **"+value[0]+" unknown to DicModName")
if len(OutSet) != 0:
    print(OutSet)
    raise ValueError("these are modules used by target files @ depend_model.mk not present at DicModName")



# build new file depend_model.mk
fNameOut=basePath+"build/depend_model.mk"
fOut=open(fNameOut,"w")
for objName,tpl in DicObjName.items():
    srcName=tpl[0]
    useSet=tpl[3]
    incSet=tpl[4]
    action=tpl[5]
    if action != None:
        strOut=NewEntryAtDepend(
            objName,
            srcName,
            useSet,
            DicModName,
            incSet,
            DicIncName,
            action)
        fOut.write(strOut+"\n")
fOut.write("include jules_depend_model.mk\n\n")
fOut.close()

print("created file",fNameOut)
