import os
import sys


############################
# INPUT DATA (MAY BE EDITED)
############################

# cloned directory name
pythonDir=os.getcwd()
homeDir=pythonDir[:-7]
buildDir=homeDir+"/build"

# desired git branch name
branch="Walcek"

# import list ObjSrcActionAtDependModel from file
# pythonDir+"/Depend_"+"_"+branch
# File contains a list of entries to be inserted at depend_model
# each entry is a triple (objName, srcName to be expanded, action to be performed)
sys.path.insert(1,pythonDir)
modName="Depend_"+branch
mod=__import__(modName)
ObjSrcActionAtDependModel=mod.ObjSrcActionAtDependModel

modName="Utils_"+branch
mod=__import__(modName)
ObjSrcActionAtDependUtils=mod.ObjSrcActionAtDependUtils



# expand string to replace source file names $() at ObjSrcActionAtDependModel
ExpandDic={
    "$(UTILS_INCS)" : "src/utils/include",
    "$(UTILS_MODS)" : "src/utils/lib/modules",
    "$(POST_SRC)" : "src/post",
    "$(UTILS_DUMP)" : "src/utils/dump",
    "$(UTILS_MDTOOLS)" : "src/utils/model-devel-tools/src",
    "$(BC)" : "src/brams/bc",
    "$(MODEL)" : "src/brams/model",
    "$(ADVC)" : "src/brams/advect",
    "$(CCATT)" : "src/brams/ccatt",
    "$(AEROSOL)" : "src/brams/ccatt",
    "$(MODEL_CHEM)" : "src/brams/ccatt/RELACS_TUV", # varies with configure
    "$(ISAN_CHEM)" : "src/brams/isan_chem",
    "$(TUV)" : "src/brams/ccatt/TUV",
    "$(CUPARM)" : "src/brams/cuparm",
    "$(FDDA)" : "src/brams/fdda",
    "$(INIT)" : "src/brams/init",
    "$(IO)" : "src/brams/io",
    "$(ISAN)" : "src/brams/isan",
    "$(MEMORY)" : "src/brams/memory",
    "$(JULES_DIR)" : "src/jules-v6.0",
    "$(MICRO)" : "src/brams/micro",
    "$(MKSFC)" : "src/brams/mksfc",
    "$(MPI)" : "src/brams/mpi",
    "$(NESTING)" : "src/brams/nesting",
    "$(RADIATE)" : "src/brams/radiate",
    "$(RRTMG_SW_SRC)" : "src/rrtmg/RRTMG_SW/src",
    "$(RRTMG_SW_MOD)" : "src/rrtmg/RRTMG_SW/modules",
    "$(RRTMG_LW_SRC)" : "src/rrtmg/RRTMG_LW/src",
    "$(RRTMG_LW_MOD)" : "src/rrtmg/RRTMG_LW/modules",
    "$(SIB)" : "src/brams/sib",
    "$(SOIL_MOISTURE)" : "src/brams/soil_moisture",
    "$(SURFACE)" : "src/brams/surface",
    "$(TEB_SPM)" : "src/brams/teb_spm",
    "$(TURB)" : "src/brams/turb",
    "$(STILT)" : "src/brams/stilt",
    "$(MATRIX)" : "src/brams/matrix",
    "$(UTILS_LIB)" : "src/utils/lib",
    "$(UTILS_TOOLS)" : "src/utils/tools",
    "$(ISAN)" : "src/brams/isan",
    "$(ISAN_MODS)" : "src/brams/isan",
    "$(WIND_FARM)" : "src/brams/wind_farm",
    "$(COMM_SPC)" : "src/brams/comm",
    "$(ENERGY)" : "src/brams/energy",
    "$(AERCLIM)" : "src/brams/aerosol",
    "$(EVAL)" : "src/brams/evaluate",
    "$(PRE-BRAMS)" : "src/pre-brams",
    "$(JULES_02)" : "src/jules-v6.0/src/control/shared",
    "$(JULES_03)" : "src/jules-v6.0/src/control/standalone",
    "$(JULES_18)" : "src/jules-v6.0/src/params/standalone",
    "$(JULES_21)" : "src/jules-v6.0/src/science/params",
    "$(JULES_26)" : "src/jules-v6.0/src/science/surface",
    "$(JULES_27)" : "src/jules-v6.0/src/science/vegetation", 
    "$(JULES_29)" : "src/jules-v6.0/src/util",
    "$(JULES_30)" : "src/jules-v6.0/src/util/shared",
    "$(AERLEVEL)" : "SIMPLE" # varies with configure
}


# types already worked upon
NewTypes=[ 'grid', 
	   'namelistfile', 
	   'controlvars', 
	   'parallelenvironment', 
	   'griddims', 
	   'domaindecomp', 
	   'nodedimensions', 
           'scalartable', 
	   'basicfields',  # conclido em 18/05/2022
	   'turbfields',  # 03/06/2022
           'miccontrol',  # 19/07/2022
           'microfields',  # 23/07/2022
           'julesfields',  # 26/07/2022
           'shcufields',  # 27/07/2022
           'gaspartfields',  # 27/07/2022
           'vartable',  # 16/08/2022
           'scalarfields',  # 03/09/2022
           'aero2mcphysfields',  # 09/09/2022
           'radiatefields',  # 28/09/2022
           'cuparmvars',  # 01/10/2022
           'cuparmfields',  # 16/10/2022
           'polygoncontainer', 
	   'neighbournodes', 
	   'messagedata', 
	   'ghostblockpartition', 
	   'fieldsectionlist', 
	   'messageset', 
	   'fieldsection', 
	   'fieldsectionnode', 
	   'gridtree' 
           ]

typesDecl=["integer", "real", "logical", "character", "type", "double precision", "complex", "external"]
lenTypesDecl=[7, 4, 7, 9, 4, 16, 7, 8]

#typeNameWord1=["integer", "real", "logical", "character", "type", "double", "complex", "external"]
typeNameWord1=["integer", "real", "logical", "character", "type", "double", "complex"]
lenTypeNameWord1=[7, 4, 7, 9, 4, 6, 7, 8]

############################
# FUNCTION DECLARATIONS
############################





def AllFilesAtDir(sourcePath, subPath, FileNames):
    """recursively find all files at sourcePath+subPath"""
    with os.scandir(sourcePath+subPath) as it:
        for entry in it:
            if entry.is_dir():
                FileNames=AllFilesAtDir(sourcePath, subPath+entry.name+"/", FileNames)
            else:
                FileNames.append(sourcePath+subPath+entry.name)
    return(FileNames)



def WC(fileNames):
    """count number of files and number of file lines at file list"""
    fNameWc="TotalLines"
    os.system("rm -rf "+fNameWc)
    os.system("touch "+fNameWc)
    for f in fileNames:
        os.system("wc -l "+f+">> "+fNameWc)
    fWc=open(fNameWc,"r")
    lineCnt=0
    fileCnt=0
    for line in fWc:
        fileCnt+=1
        lineSplit=line.split()
        lineCnt+=int(lineSplit[0])
    fWc.close()
    os.system("rm -rf "+fNameWc)
    return(fileCnt,lineCnt)




def FindFortranFiles(sourcePath, subPath, FileNames):
    """recursively find all files with suffix .f90 or .F90 at sourcePath+subPath"""
    with os.scandir(sourcePath+subPath) as it:
        for entry in it:
            if entry.is_dir():
                FileNames=FindFortranFiles(sourcePath, subPath+entry.name+"/", FileNames)
            else:
                nameSplit=entry.name.split(".")
                if len(nameSplit) > 1:
                    suffix=nameSplit[len(nameSplit)-1]
                    if (suffix=="f90") or (suffix=="F90"):
                        FileNames.append(sourcePath+subPath+entry.name)
    return(FileNames)




def ModuleName(fName):
    """all module names @ file module statements"""
    fIn=open(fName,"r")
    ret=[]
    for line in fIn:
        lineLower=line[:-1].lower().split()
        if (len(lineLower)==2) and (lineLower[0] == "module"):
            modName=lineLower[1]
            ret.append(modName)
    fIn.close()
    return(ret)



def UseSet(fName):
    """set of module names used @ file"""
    Ret=set()
    fIn=open(fName,"r")
    for line in fIn:
        lineLower=line.lower()
        if "use" in lineLower:
            lineSplit=lineLower[:-1].split()
            if len(lineSplit) > 1:
                if lineSplit[0] == "use":
                    if "," in line:
                        modUsed=lineSplit[1].split(",")[0]
                    else:
                        modUsed=lineSplit[1]
                    Ret.add(modUsed)
    return(Ret)





def UseWithoutOnlySet(fName):
    """set of module names used @ file without only"""
    Ret=set()
    fIn=open(fName,"r")
    for line in fIn:
        lineLower=line.lower()
        if "use" in lineLower:
            lineSplit=lineLower[:-1].split()
            if len(lineSplit) > 1:
                if lineSplit[0] == "use":
                    if "," not in line:
                        modUsed=lineSplit[1].split(",")[0]
                        Ret.add(modUsed)
    return(Ret)


def IncSet(fName):
    """all include file names @ file include statements"""
    Ret=set()
    fIn=open(fName,"r")
    for line in fIn:
        lineLower=line.lower()
        if "include" in lineLower:
            lineSplit=lineLower[:-1].split()
            if len(lineSplit) > 1:
                if lineSplit[0] == "include":
                    incUsed=lineSplit[1].strip('"').strip("'")
                    Ret.add(incUsed)
    return(Ret)



def RemoveModuleNameFromUseSet(useSet,modName):
    """remove modName from useSet"""
    res=set()
    for mod in useSet:
        if mod != modName:
            res.add(mod)
    return(res)
            


        
def NewEntryAtDepend(objName, srcName, useSet, DicModName,
                     incSet, DicIncName, action):
    """build one entry of depend file"""

    # objName depends upon source file
    strOut=objName+" : "+srcName+" "
    lenStart=len(strOut)

    # insert use dependency, replacing module name by objName
    for use in useSet:
        if use not in DicModName.keys():
            raise ValueError(use+" not in DicModName")
        elif DicModName[use] != None:
            strOut+=DicModName[use]+" "

    # insert include dependency, replacing file name by objName
    for inc in incSet:
        if inc not in DicIncName.keys():
            raise ValueError(inc+" not in DicIncName")
        elif DicIncName[inc] != None:
            strOut+=DicIncName[inc]+" "

    # break strOut in lines of maximum length
    maxLineLength=80
    strSplit=strOut.split()
    line=strSplit[0]+" "+strSplit[1]+" "+strSplit[2]+" "
    lineLen=len(line)
    for name in strSplit[3:]:
        if lineLen+len(name) > maxLineLength:
            line+="\\\n\t"+name+" "
            lineLen=len("\t"+name+" ")
        else:
            line+=name+" "
            lineLen+=len(name+" ")
    return(line+"\n"+action)



def BuildDicObjName(ObjSrcActionAtDependModel, ExpandDic):
    """build dictionary indexed by object name composed by tuples

    Tuples are composed by source name without expansion, source
    name expanded (full path), list of declared modules, set of
    used modules, set of included files and action at makefile
"""
    DicObjName={}
    for tpl in ObjSrcActionAtDependModel:

        # object (target) name
        objName=tpl[0]
        if objName in DicObjName.keys():
            raise ValueError(objName+" repicated at DicObjName")

        # corresponding source file with $(XXX)
        srcName=tpl[1]
    
        # expand $() at source file name
        srcNameSplit=srcName.split("/")
        if srcNameSplit[0] not in ExpandDic:
            raise ValueError ("for "+srcName+", "+srcNameSplit[0]+" not in ExpandDic")
        expandedSrcName = srcName.replace(srcNameSplit[0],ExpandDic[srcNameSplit[0]])
        # expand twice $() at srcName
        if "$" in expandedSrcName:
            indStart=expandedSrcName.find("$")
            indEnd=expandedSrcName.find(")",indStart)
            toReplace=expandedSrcName[indStart:indEnd+1]
            if toReplace not in ExpandDic.keys():
                raise ValueError ("for "+expandedSrcName+", "+toReplace+" not in ExpandDic")
            replaceStr=ExpandDic[toReplace]
            expandedSrcName=expandedSrcName.replace(toReplace,replaceStr)
        expandedSrcName=homeDir+"/"+expandedSrcName

        # action to be performed
        action=tpl[2]
        if action == 'None':
            action=None

        # modules declared at source file
        modNameList=ModuleName(expandedSrcName)

        # set of module names cited at source file use statements
        useSet=UseSet(expandedSrcName)

        # set of include names cited at source file include statements
        incSet=IncSet(expandedSrcName)

        # Dictionary entry
        DicObjName[objName]=[
            srcName,
            expandedSrcName,
            modNameList,
            useSet,
            incSet,
            action
        ]
    return(DicObjName)



def BuildDicModName(DicObjName):
    DicModName={}
    for key, val in DicObjName.items():
        modNameList=val[2]
        for modName in modNameList:
            if modName in DicModName.keys():
                raise ValueError("module "+modName+" replicated at DicModName")
            DicModName[modName]=key
    return(DicModName)




def AllEntriesAtDepend(fNameDepend):
    """ingest depend_model file into tuples"""
    fDep=open(fNameDepend,"r")
    Res=[]
    searchObjName=True
    consumeDepList=False
    getAction=False
    for line in fDep:
        if searchObjName:
            if ":" in line:
                lineSplit=line.split(":")
                if len(lineSplit) != 2:
                    print(line)
                    raise ValueError("line with more than one :")
                else:  # single ":" only at a starting line of a new object
                    objName=lineSplit[0].strip() # object name prior to ":"
                    srcName=lineSplit[1].strip().split()[0]
                    searchObjName=False
                    consumeDepList=lineSplit[1][-2]=="\\"
                    getAction=not consumeDepList
                    action=""
        elif consumeDepList:
            if line[-2] != "\\": # dependency list ends if not a continuation line
                consumeDepList=False
                getAction=True
        elif getAction:
            if len(line[:-1].strip()) == 0: # blank line finishes action
                newEntry=objName,srcName,action
                Res.append(newEntry)
                searchObjName=True
                consumeDepList=False
                getAction=False
            else:
                action+=line
        else:
            raise RuntimeError("programming error at input line **"+line+\
                               "**: none of searchObjName, consumeDepList, getAction is true")
    if getAction:
        newEntry=objName,srcName,action
        Res.append(newEntry)
    return(Res)




def BuildDicTypes(DicObjName):
    """Dictionary indexed by object name, entries set of type definition and set of type declarations"""
    DicTypes={}
    for objName, tplDicObjName in DicObjName.items():
        srcLocalPath=tplDicObjName[1]
        defSet=set()
        useSet=set()
        fSrc=open(srcLocalPath,"r")
        for line in fSrc:
            lineLow=line.lower()[:-1]
            lineLowSplit=lineLow.split()
            lll = len(lineLowSplit)
            if (lll < 2):
                continue
            if "type" in lineLowSplit[0]:
                if (lineLowSplit[0] == "type") and (lineLowSplit[1][0] != "("):
                    typeName=lineLowSplit[1].strip()
                    defSet.add(typeName)
                elif lineLowSplit[0][0:4] == "type":
                    if "(" not in lineLow:
                        typeName=lineLowSplit[-1].strip()
                        defSet.add(typeName)
                    else:
                        parStart=lineLow.find("(")
                        parEnd=lineLow.find(")")
                        typeName=lineLow[parStart+1:parEnd].strip()
                        useSet.add(typeName)
        fSrc.close()
        DicTypes[objName]=(defSet,useSet)
    return(DicTypes)



def StripCommentsLowerNewLine(line):
    """remove comments, newline, leading and trailing blanks, force down case"""
    ll=line.lower()[:-1].strip()
    if "!" not in ll:
        return(ll)
    elif len(ll) == 0:
        return("")
    else:
        ind=ll.find("!")
        if ind == 0:
            return("")
        else:
            return(ll[:ind])





        
def ModuleTextAtFile(fName):
    """split file text into lists of modules texts and outside modules texts"""
    f=open(fName,"r")
    moduleList=[]
    moduleText=[]
    outsideList=[]
    outsideText=[]
    inModule=False
    for line in f:
        stripLine=StripCommentsLowerNewLine(line)
        if len(stripLine) > 0:
            splitLine=stripLine.split()
            if (not inModule):
                if "module" in splitLine[0]:
                    inModule=True
                    moduleText.append(stripLine)
                    if len(outsideText) != 0:
                        outsideList.append(outsideText)
                    outsideText=[]
                else:
                    outsideText.append(stripLine)
            else:
                if len(splitLine) != 0:
                    moduleText.append(stripLine)
                    if (splitLine[0]=="end") and (splitLine[1]=="module"):
                        inModule=False
                        moduleList.append(moduleText)
                        moduleText=[]
                        outsideText=[]
    else:
        if len(outsideText) != 0:
            outsideList.append(outsideText)
    f.close()
    return(moduleList,outsideList)





def RemoveContinuation(lineList):
    """colapse Fortran continuation lines into a single line"""
    noCont=[]
    textCont=""
    for line in lineList:
        line=line.strip()
        if line[-1] == "&":
            if "&" in line[:-1]:
                textCont += line[:-1].strip("&").strip("\t")
            else:
                textCont += line[:-1]
        elif textCont == "":
            noCont.append(line.strip("&").strip("\t"))
        else:
            noCont.append(textCont+line.strip("&").strip("\t"))
            textCont=""
    return(noCont)


def SplitProcList(someList):
    """Break procedure list into list of procedures"""
    typeListOneWord=["integer",
                     "real",
                     "logical",
                     "complex",
                     "character",
                     "pure",
                     "elemental",
                     "recursive"]
    inSubr=False
    inFunc=False
    procList=[]
    procText=[]
    noCont=RemoveContinuation(someList)
    for line in noCont:
        lineSplit=line.split()
        if inSubr:
            if len(lineSplit)>1:
                if (lineSplit[0]=="end") and (lineSplit[1]=="subroutine"):
                    inSubr=False
                    procText.append(line)
                    procList.append(procText)
                    procText=[]
                elif line.split()[0] == "subroutine":
                    continue
#                    print("!!!! subroutine starts within a subroutine !!!")
        elif line.split()[0] == "subroutine":
            inSubr=True
        elif "function" in line:
            if len(lineSplit) > 1:
                if lineSplit[0] == "function":
                    inFunc=True
                elif (lineSplit[0]=="end") and (lineSplit[1]=="function"):
                    if not inFunc:
                        print("function end outside function")
                        print(line)
                        raise ValueError("function end outside function")
                    procText.append(line)
                    procList.append(procText)
                    procText=[]
                    inFunc=False
                elif lineSplit[1]=="function":
                    if lineSplit[0] in typeListOneWord:
                        inFunc=True
                    elif "(" in lineSplit[0]:
                        lineOpen=lineSplit[0].split("(")
                        if lineOpen[0] in typeListOneWord:
                            inFunc=True
                    elif "*" in lineSplit[0]:
                        lineOpen=lineSplit[0].split("*")
                        if lineOpen[0] in typeListOneWord:
                            inFunc=True
                elif len(lineSplit) > 2:
                    if lineSplit[2]=="function":
                        if (lineSplit[0]=="double") and (lineSplit[1]=="precision"):
                            inFunc=True
                        elif lineSplit[0]=="elemental":
                            inFunc=True
                        elif lineSplit[0]=="pure":
                            inFunc=True
        if inFunc or inSubr:
            procText.append(line)
    return(procList)


def NameArgListFromInterface(inter):
    """ split procedure first line into name and list of arguments"""
    
    argList=[]
    
    openPar="(" in inter
    closePar=")" in inter

    # case name
    if (not openPar) and (not closePar):
        name=inter.strip()
    
    # inconsistent interface 
    elif (openPar and not closePar) or \
         (closePar and not openPar):
        raise ValueError("inconsistent interface at **"+inter+"**")

    # get name and arg list of interface with "(" and ")"
    else:
        indOpen=inter.index("(")
        indClose=inter.index(")")

        name=inter[:indOpen].strip()
        
        for arg in inter[indOpen+1:indClose].split(","):

            # remove case name()
            if len(arg) != 0:
                argList.append(arg.strip())
    return(name,argList)


def isDeclaration(line):
    """is this line a declaration?"""
    if "::" in line:
        lineBreak=line.split("::")
    elif " " in line:
        lineBreak=line.split(" ")
    else:
#        print("isDeclaration: not decl **"+line+"**")
        return(False)
    for n,d in enumerate(typesDecl):
        if len(lineBreak[0]) >= lenTypesDecl[n]:
            if d == lineBreak[0][0:lenTypesDecl[n]]:
#              print("isDeclaration: is **"+line+"**")
              return(True)
    else:
#        print("isDeclaration: not decl **"+line+"**")
        return(False)

def SplitDeclaration(line):
    varString=""
    attString=""

    print("Split Declaration: linha entrada=",line)
    if not isDeclaration(line):
        return(False,attString,varString)
              
    if "::" in line:
        lineBreak=line.split("::")
        for n,d in enumerate(typesDecl):
            if len(lineBreak[0]) >= lenTypesDecl[n]:
                if d == lineBreak[0][0:lenTypesDecl[n]]:
                    attString=lineBreak[0].expandtabs().strip().replace(" ","")
                    varString=lineBreak[1].expandtabs().replace(" ","")
                    print("SplitDeclaration sep ::",
                          "attributos **",attString,"**",
                          "nomes **",varString,"**")
                    return(True,attString,varString)
    elif " " in line:
        lineBreak=line.strip("\t").split()
#        print(lineBreak)
        for n,d in enumerate(typesDecl):
            if len(lineBreak[0]) >= lenTypesDecl[n]:
                if d == lineBreak[0][0:lenTypesDecl[n]]:
                    attString=lineBreak[0]
                    indStart=1
                    if lineBreak[1]=="*(*)":
                        attString+=lineBreak[1]
                        indStart=2
#                        print("ACHOU",attString)
                    elif lineBreak[1][0]=="(":
                        attString+=lineBreak[1]
                        indStart=2
                    for part in lineBreak[indStart:]:
                        varString+=part
                    print("SplitDeclaration sep branco",
                          "attributos",attString,
                          "nomes",varString)
                    return(True,attString,varString)
    return(False,attString,varString)

def ProcInterface(proc):
    """Split procedure interface into name and argList"""
    if len(proc) == 0:
        raise ValueError("proc has length 0")
    if "subroutine" in proc[0]:
        # name and argument list from subroutine declaration
        subInterface=proc[0].split(" ",1)[1].replace(" ","")
        name,argList=NameArgListFromInterface(subInterface)
    elif "function" in proc[0]:
        # name and argument list from function declaration
        subInterface=proc[0].split("function",1)[-1].replace(" ","")
        name,argList=NameArgListFromInterface(subInterface)
    else:
        raise ValueError("no procedure declaration")
    return(name,argList)


#def ProcInterfaceDeclarations(proc):
#    """Split procedure text into procedure name, argument and declarations lists"""
#    declList=[]
#    if len(proc) == 0:
#        raise ValueError("empty procedure text")
#    if "subroutine" in proc[0]:
#        # name and argument list from subroutine declaration
#        subInterface=proc[0].split(" ",1)[1].replace(" ","")
#        name,argList=NameArgListFromInterface(subInterface)
#        for line in proc[1:-1]:
#            if isDeclaration(line):
#                declList.append(line)
#    elif "function" in proc[0]:
#        # name and argument list from function declaration
#        subInterface=proc[0].split("function",1)[-1].replace(" ","")
#        name,argList=NameArgListFromInterface(subInterface)
#        for line in proc[1:-1]:
#            if isDeclaration(line):
#                declList.append(line)
#    else:
#        raise ValueError("no procedure declaration")
#    print("ProcInterfaceDeclaration:")
#    print("procedure",proc)
#    print("name **"+name+"**")
#    print("argList=",argList)
#    print("declList=",declList)
#    return(name,argList,declList)


#########################################################

def Tokenize(line):
    """Break line in tokens

    Tokens are names, numbers, expression, strings, (, ), ::, =, single : and *
    """

    tokenList=[]
    posTokenStart=0
    blankSeq=False
    inString=False
    lineStrip=line.replace("\t"," ").strip()
    for n, c in enumerate(lineStrip):

        # consuming strings
        if (c=="'") or (c=='"'):
            if inString:
                tokenCand=lineStrip[posTokenStart:n+1]
                if len(tokenCand) > 0:
                    tokenList.append(tokenCand)
                tokenList.append(c)
                posTokenStart=n+1
                inString=False
            else:
                tokenCand=lineStrip[posTokenStart:n]
                if len(tokenCand) > 0:
                    tokenList.append(tokenCand)
                posTokenStart=n
                inString=True
                
        # consuming blank sequence
        elif inString:
            continue
        
        elif c==" ":
            if not blankSeq:
                blankSeq=True
                tokenCand=lineStrip[posTokenStart:n]
                if len(tokenCand) > 0:
                    tokenList.append(tokenCand)
            posTokenStart=n+1
        else:
            blankSeq=False
            
            if c=="(":
                tokenCand=lineStrip[posTokenStart:n]
                if len(tokenCand) > 0:
                    tokenList.append(tokenCand)
                tokenList.append(c)
                posTokenStart=n+1
            elif c==")":
                tokenCand=lineStrip[posTokenStart:n]
                if len(tokenCand) > 0:
                    tokenList.append(tokenCand)
                tokenList.append(c)
                posTokenStart=n+1
            elif c=="=":
                tokenCand=lineStrip[posTokenStart:n]
                if len(tokenCand) > 0:
                    tokenList.append(tokenCand)
                tokenList.append(c)
                posTokenStart=n+1
            elif c=="," :
                tokenCand=lineStrip[posTokenStart:n]
                if len(tokenCand) > 0:
                    tokenList.append(tokenCand)
                posTokenStart=n+1
            elif c==":":
                if lineStrip[n-1] == ":":
                    tokenList.append("::")
                    posTokenStart=n+1
                elif len(lineStrip) == n+1:
                    tokenCand=lineStrip[posTokenStart:n]
                    if len(tokenCand) > 0:
                        tokenList.append(tokenCand)
                    posTokenStart=n
                elif lineStrip[n+1] == ":":
                    tokenCand=lineStrip[posTokenStart:n]
                    if len(tokenCand) > 0:
                        tokenList.append(tokenCand)
                    posTokenStart=n
                else:
                    tokenCand=lineStrip[posTokenStart:n]
                    if len(tokenCand) > 0:
                        tokenList.append(tokenCand)
                    tokenList.append(c)
                    posTokenStart=n+1
            elif c=="*":
                tokenCand=lineStrip[posTokenStart:n]
                if len(tokenCand) > 0:
                    tokenList.append(tokenCand)
                tokenList.append(c)
                posTokenStart=n+1
    else:
        tokenCand=lineStrip[posTokenStart:]
        if len(tokenCand) > 0:
            tokenList.append(tokenCand)
            
    return(tokenList)



def TypeName(tokenList):
    # type name
    if len(tokenList) < 2:
        return("")
    if (tokenList[0]=="double") and \
       (tokenList[1]=="precision"):
        return("double precision")
    elif tokenList[0] in typeNameWord1:
        return(tokenList[0])
    else:
        return("")



def isDec(tokenList):
    """is this line a declaration?"""
    res=(TypeName(tokenList) != "") or \
        (tokenList[0]=="external")
    return(res)


def KindName(tokenList,nextToken):
    kindName=""
    cntPar=0
    isPar=tokenList[nextToken] == "("
    isAst=tokenList[nextToken] == "*"
    if isAst:
        if len(tokenList) > nextToken+1:
            if tokenList[nextToken+1].isdecimal():
                kindName="*"+tokenList[nextToken+1]
                nextToken+=2
                return(nextToken,kindName)
        else:
            print(tokenList)
            raise ValueError("unknown text after *")
    if isPar or isAst:
        for nbr, token in enumerate(tokenList[nextToken:]):
            kindName+=token
            if token == "(":
                cntPar+=1
            elif token == ")":
                cntPar-=1
                if cntPar==0:
                    nextToken += nbr+1
                    break
        else:
            print(tokenList)
            raise ValueError("unclosed open par")            
    return(nextToken,kindName)

    


def IntentName(tokenList,nextToken):
    intentName=""
    cntPar=0
    if tokenList[nextToken] == "intent":
        nextToken+=1
        for nbr, token in enumerate(tokenList[nextToken:]):
            intentName+=token
            if token == "(":
                cntPar+=1
            elif token == ")":
                cntPar-=1
                if cntPar==0:
                    nextToken += nbr+1
                    break
        else:
            print(tokenList)
            raise ValueError("unclosed open par")            
    return(nextToken,intentName)


def DimensionBounds(tokenList,nextToken):
    dimensionBounds=""
    cntPar=0
    for nbr, token in enumerate(tokenList[nextToken:]):
        dimensionBounds+=token
        if token == "(":
            cntPar+=1
        elif token == ")":
            cntPar-=1
            if cntPar==0:
                nextToken += nbr+1
                break
    else:
        print(tokenList)
        raise ValueError("unclosed open par")            
    return(nextToken,dimensionBounds)
    


def DimensionName(tokenList,nextToken):
    dimensionName=""
    cntPar=0
    if tokenList[nextToken] == "dimension":
        nextToken,dimensionName=DimensionBounds(tokenList,nextToken+1)
    return(nextToken,dimensionName)

    


def VarNames(tokenList):
    dicNames={}
    currName=""
    currDim=""
    inDim=False
    cntPar=0
    for nbr, tok in enumerate(tokenList):
        if inDim:
            if tok == ")":
                cntPar-=1
                if cntPar == 0:
                    inDim=False
                    if currDim[-1]==",":
                        currDim=currDim[:-1]+tok
                    else:
                        currDim+=tok
                else:
                    currDim+=tok+","
            else:
                currDim+=tok+","
        elif tok == "(":
            cntPar=1
            currDim=tok
            inDim=True
        elif currName == "":
            currName=tok
        else:
            dicNames[currName]=currDim
            currDim=""
            currName=tok
    dicNames[currName]=currDim
    return(dicNames)



def ParseDeclList(declList):
    AfterKind=["save","parameter","intent","dimension","external","pointer","allocatable","optional"]
    # declDic indexed by declared variable name
    # entry is a tuple with:
    #   type
    #   kind
    #   save
    #   intent
    #   dimension
    #   optional
    #   pointer
    #   allocatable
    declDic={}
    for tokenList in declList:
        # type name
        typeName=TypeName(tokenList)
        if len(typeName) > 0:
            if typeName == "double precision":
                nextToken=2
            else:
                nextToken=1
        else:
            nextToken=0 # for line "external ..."
        nextToken,kindName=KindName(tokenList,nextToken)
        saveName=""
        parameterName=""
        intentName=""
        dimensionName=""
        externalName=""
        optionalName=""
        pointerName=""
        allocatableName=""
        while tokenList[nextToken] in AfterKind:
            if tokenList[nextToken] == "save":
                saveName="save"
                nextToken+=1
            elif tokenList[nextToken] == "parameter":
                parameterName="parameter"
                nextToken+=1
            elif tokenList[nextToken] == "intent":
                nextToken,intentName=IntentName(tokenList,nextToken)
            elif tokenList[nextToken] == "dimension":
                nextToken,dimensionName=DimensionName(tokenList,nextToken)
            elif tokenList[nextToken] == "external":
                externalName="external"
                nextToken+=1
            elif tokenList[nextToken] == "optional":
                optionalName="optional"
                nextToken+=1
            elif tokenList[nextToken] == "pointer":
                pointerName="pointer"
                nextToken+=1
            elif tokenList[nextToken] == "allocatable":
                allocatableName="allocatable"
                nextToken+=1
            if nextToken >= len(tokenList):
                raise ValueError("search exceeded tokenList")
        if tokenList[nextToken] == "::":
            nextToken+=1
        if (parameterName == ""):
            dicNames=VarNames(tokenList[nextToken:])
            for key,val in dicNames.items():
                if val !="":
                    declDic[key]=(typeName,\
                                  kindName,\
                                  saveName,\
                                  intentName,\
                                  val,\
                                  optionalName,\
                                  pointerName,\
                                  allocatableName)
                else:
                    declDic[key]=(typeName,\
                                  kindName,\
                                  saveName,\
                                  intentName,\
                                  dimensionName,\
                                  optionalName,\
                                  pointerName,\
                                  allocatableName)
    return(declDic)


def HasIntent(declDicEntry):
    return(declDicEntry[3]!="")

def IsArray(declDicEntry):
    return(declDicEntry[4]!="")

def IsAssumedShape(declDicEntry):
    charAssum=["(",")",",",":"]
    if declDicEntry[4]=="":
        return(False)
    elif declDicEntry[4][0] != "(":
        raise ValueError("unknown dimension")
    else:
        for c in declDicEntry[4]:
            if c not in charAssum:
                return(False)
        return(True)


            

def ProcInterfaceDeclarations(proc):
    """Split procedure text into procedure name, argument and declarations lists"""
    declList=[]
    if len(proc) == 0:
        raise ValueError("empty procedure text")
    if "subroutine" in proc[0]:
        # name and argument list from subroutine declaration
        subInterface=proc[0].split(" ",1)[1].replace(" ","")
        name,argList=NameArgListFromInterface(subInterface)
        for line in proc[1:]:
            tokenList=Tokenize(line)
            if isDec(tokenList):
                declList.append(tokenList)
    elif "function" in proc[0]:
        # name and argument list from function declaration
        subInterface=proc[0].split("function",1)[-1].replace(" ","")
        name,argList=NameArgListFromInterface(subInterface)
        for line in proc[1:]:
            tokenList=Tokenize(line)
            if isDec(tokenList):
                declList.append(tokenList)
    else:
        raise ValueError("no procedure declaration")

    declDic=ParseDeclList(declList)
    argDic={}
    undeclaredList=[]
    for arg in argList:
        if arg in declDic.keys():
            argDic[arg]=declDic[arg]
        else:
            undeclaredList.append(arg)
            
    return(name,argDic,undeclaredList)
        

