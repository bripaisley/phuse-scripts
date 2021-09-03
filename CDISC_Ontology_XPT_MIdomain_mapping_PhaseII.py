###########################################################################
###                        Loading Packages                             ###
###########################################################################
import os
import sys
import io
import requests
import warnings
import time
import glob
import re
#import wget
import pandas as pd
from functools import reduce
from collections import defaultdict
import chardet

#import csv
#import pathlib
pd.options.mode.chained_assignment = None

###########################################################################
###                     Helper Functions                                ###
###########################################################################
def get_actual_filename(name):
    name = "%s[%s]" % (name[:-1], name[-1])
    return glob.glob(name)[0]

def Domain_test(mydir):
    fMI = None
    fTS = None
    fEX = None
    fDM = None
    fDS = None

    try:
        MIfile = os.path.join(mydir, "mi.xpt")
        fMI = get_actual_filename(MIfile)
    except:
        pass

    try:
        TSfile = os.path.join(mydir, "ts.xpt")
        fTS = get_actual_filename(TSfile)
    except:
        pass

    try:
        EXfile = os.path.join(mydir, "ex.xpt")
        fEX = get_actual_filename(EXfile)
    except:
        pass

    try:
        DSfile = os.path.join(mydir, "ds.xpt")
        fDS = get_actual_filename(DSfile)
    except:
        pass

    try:
        DMfile = os.path.join(mydir, "dm.xpt")
        fDM = get_actual_filename(DMfile)
    except:
        pass

    return fMI, fTS, fEX, fDM, fDS

def Add_Term(columnName, row, term):
    for index in row.index:
        valueOut = ""
        val = row[columnName]
        if(term in row.at[columnName]):
            valueOut = val
        elif(val and not val.isspace()):
            strings = [row[columnName],term]
            valueOut = '; '.join(strings)
        else:
            valueOut = term
    return valueOut

def Search_Term(rowSubset, terms):
    booleanOut = ""
    booleanOut = any(term in str(rowSubset) for term in terms)
    return booleanOut

def rowFindingSEND(row):

    row = row.replace(to_replace =[",",";"], value=" ", regex=True).astype(str) #multiple white-space to 1 space "\\"
    row = row.replace(to_replace =" +", value=" ", regex=True).astype(str) #multiple white-space to 1 space "\\"
    rowFindingSEND = []
    rowtemp = row[['MIORRES_raw', 'MISTRESC_raw']].copy()
    rowFindingSEND = rowtemp.values.tolist()
    return rowFindingSEND

def codetermExtract(dfSENDCodelist, codeVal):

    codeTerm = dfSENDCodelist.loc[dfSENDCodelist["Code"] == codeVal, "CDISC Submission Value"].item()
    return codeTerm

def invert_dict(d):
    reversed_d = dict()
    for k, v in d.items():
        if isinstance(v, list): #if value is a list of items
            for item in v:
                reversed_d[item] = k
        else: #if value is a single item
            reversed_d[v] = k

    return reversed_d

###########################################################################
###                     Domain Data Extraction                          ###
###########################################################################

def MI_dataframe(fMI):
    with open(fMI, 'rb') as f:
        result = chardet.detect(f.read())  # or readline if the file is large
    dfMItemp = pd.read_sas(fMI, encoding=result['encoding']).drop_duplicates()

    dfMI = dfMItemp[dfMItemp['MITESTCD'] == "GHISTXQL"].drop_duplicates() #General Histopathologic Exam, Qual
    if(len(dfMI) != len(dfMItemp)):
        MIdiff = len(dfMItemp)-len(dfMI)
        warnings.warn(str(MIdiff) + " rows are not General Histopathologic Exam, Qual and will be excluded")

    dfuANIMAL = dfMI.reindex(columns=["STUDYID", "USUBJID"]).drop_duplicates()

    return dfMI, dfuANIMAL

def EX_dataframe(fEX, dfuANIMAL):
    with open(fEX, 'rb') as f:
        result = chardet.detect(f.read())  # or readline if the file is large
    dfEXtemp = pd.read_sas(fEX, encoding=result['encoding']).drop_duplicates()

    dfEXtemp2 = dfEXtemp.reindex(columns=["STUDYID", "USUBJID", "EXROUTE"])
    dfRoute = pd.merge(dfuANIMAL, dfEXtemp2, how='left', left_on=["STUDYID", "USUBJID"], right_on=["STUDYID", "USUBJID"]).drop_duplicates()

    return dfRoute

def TS_dataframe(fTS, dfuANIMAL):
    with open(fTS, 'rb') as f:
        result = chardet.detect(f.read())  # or readline if the file is large
    dfTStemp = pd.read_sas(fTS, encoding=result['encoding']).drop_duplicates()

    #Returns values present, no error if any are missing. Ex: if ROUTE missing, then other still outputs.
    dfTStemp2 = dfTStemp[dfTStemp["TSPARMCD"].isin(["STUDYID", "ROUTE"])].drop_duplicates()
    dfTStemp3 = dfTStemp2.reindex(columns=["STUDYID", "TSVAL"])
    dfRoute = pd.merge(dfuANIMAL, dfTStemp3, how='left', left_on=["STUDYID"], right_on=["STUDYID"]).drop_duplicates()
    dfRoute = dfRoute.rename(columns={'TSVAL': 'EXROUTE'})  ##RENAME ROUTE TO EXROUTE

    dfTStemp2= dfTStemp[dfTStemp["TSPARMCD"].isin(["STUDYID", "SPECIES"])].drop_duplicates()
    dfTStemp3 = dfTStemp2.reindex(columns=["STUDYID", "TSVAL"])
    dfSpecies = pd.merge(dfuANIMAL, dfTStemp3, how='left', left_on=["STUDYID"], right_on=["STUDYID"]).drop_duplicates()
    dfSpecies = dfSpecies.rename(columns={'TSVAL': 'SPECIES'})  ##RENAME ROUTE TO EXROUTE

    dfTStemp2= dfTStemp[dfTStemp["TSPARMCD"].isin(["STUDYID", "STRAIN"])].drop_duplicates()
    dfTStemp3 = dfTStemp2.reindex(columns=["STUDYID", "TSVAL"])
    dfStrain = pd.merge(dfuANIMAL, dfTStemp3, how='left', left_on=["STUDYID"], right_on=["STUDYID"]).drop_duplicates()
    dfStrain = dfStrain.rename(columns={'TSVAL': 'STRAIN'})  ##RENAME ROUTE TO EXROUTE

    return dfRoute, dfSpecies, dfStrain

def DS_dataframe(fDS, dfuANIMAL):
    with open(fDS, 'rb') as f:
        result = chardet.detect(f.read())  # or readline if the file is large
    dfDStemp = pd.read_sas(fDS, encoding=result['encoding']).drop_duplicates()

    dfDStemp2 = dfDStemp.reindex(columns=["STUDYID", "USUBJID", "DSDECOD"])
    dfDisposition = pd.merge(dfuANIMAL, dfDStemp2, how='left', left_on=["STUDYID", "USUBJID"], right_on=["STUDYID", "USUBJID"]).drop_duplicates()

    return dfDisposition

def DM_dataframe(fDM, dfuANIMAL):
    with open(fDM, 'rb') as f:
        result = chardet.detect(f.read())  # or readline if the file is large
    dfDMtemp = pd.read_sas(fDM, encoding=result['encoding']).drop_duplicates()
    dfDMtemp2 = dfDMtemp.reindex(columns=["STUDYID", "USUBJID", "SEX"])
    dfGender = pd.merge(dfuANIMAL, dfDMtemp2, how='left', left_on=["STUDYID", "USUBJID"], right_on=["STUDYID", "USUBJID"]).drop_duplicates()

    return dfGender

def Raw_dataframe(fMI, fTS, fEX, fDM, fDS):

    dfMI, dfuANIMAL = MI_dataframe(fMI)

    ##Generate Route(EX/TS) and Species(TS) domain dataframe
    if fEX == None and fTS==None:
        warnings.warn("No EX.xpt or TS.xpt file detected. Route cannot be determined.")  # if no EX or TS domain, print warning and proceed
        dfSpecies = dfuANIMAL.copy()
        dfSpecies["SPECIES"] = ""
        dfStrain = dfuANIMAL.copy()
        dfStrain["STRAIN"] = ""
        dfRoute = dfuANIMAL.copy()
        dfSpecies["EXROUTE"] = ""

    elif fEX==None:
        warnings.warn("No EX.xpt file detected. Route will be determined from TS domain.")  # if no EX domain, print warning and proceed
        dfRoute, dfSpecies, dfStrain = TS_dataframe(fTS, dfuANIMAL)

    elif fTS==None:
        warnings.warn("No TS.xpt file detected. Species cannot be determined.")  # if no TS domain, print warning and proceed

        dfRoute = EX_dataframe(fEX, dfuANIMAL)
        dfSpecies = dfuANIMAL.copy()
        dfSpecies["SPECIES"] = ""
        dfStrain = dfuANIMAL.copy()
        dfStrain["STRAIN"] = ""

    else:
        dfRoute, dfSpecies, dfStrain = TS_dataframe(fTS, dfuANIMAL)
        dfRoute = EX_dataframe(fEX, dfuANIMAL)

    ##Generate DS(Disposition) domain dataframe
    if fDS == None:
        warnings.warn("No DS.xpt file detected.")  # if no DS domain, print warning and proceed

        ##additional logic to address missing DS domain
        dfDisposition = dfuANIMAL.copy()
        dfDisposition["DSDECOD"] = ""

    else:
        dfDisposition = DS_dataframe(fDS, dfuANIMAL)

    ##Generate DM(Gender) domain dataframe
    if fDM == None:
        warnings.warn("No DM.xpt file detected.")  # if no DM domain, print warning and proceed

        ##additional logic to address missing DM domain
        dfGender = dfuANIMAL.copy()
        dfGender["SEX"] = ""

    else:
        dfGender = DM_dataframe(fDM, dfuANIMAL)

    dfList = [dfuANIMAL, dfSpecies, dfStrain, dfGender, dfRoute, dfDisposition]
    dfAnimalInfo = reduce(lambda  left,right: pd.merge(left, right, on=["STUDYID", "USUBJID"], how='left'), dfList)
    df = pd.merge(dfMI, dfAnimalInfo, how='left', left_on=["STUDYID", "USUBJID"], right_on=["STUDYID", "USUBJID"]).drop_duplicates()

    return(df)

###########################################################################
###                     Domain Data Extraction                          ###
###########################################################################

def dfClean(df):

    df = df.fillna("")
    df = df.apply(lambda x: x.astype(str).str.upper())
    df = df.replace(to_replace =[",","-","=","/","\[","]","\(","\)",";",":","\.","\\\\"], value=" ", regex=True).astype(str) #multiple white-space to 1 space "\\"
    df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    df = df.replace(to_replace =" +", value=" ", regex=True).astype(str) #multiple white-space to 1 space "\\"

    return df

def dfInitialization(df):
    df = df.fillna("")
    df = df.add_suffix('_raw')
    df["id"] = df.index + 1
    dftemp = df.copy()
    dftemp = dfClean(dftemp)

    dftemp["SEND_Domain"] = "MI"

    ##Animal Details columns
    dftemp["SEX"] = "" #SEND: Sex
    dftemp["SPECIES"] = "" #SEND: Species
    dftemp["STRAIN"] = ""  # SEND: Strain
    dftemp["EXROUTE"] = ""  # SEND: Route of Administration
    dftemp["DSDECOD"] = ""  # SEND: Animal Disposition

    # Specimen columns
    dftemp["MISPEC"] = "" #SEND: Specimen
    dftemp["MIANTREG"] = "" #SEND: Anatomical Location
    dftemp["Sample_Type_Group"] = "" #Non-standardized
    dftemp["Anatomical_Region_of_Specimen"] = "" #Non-standardized
    dftemp["Celltype"] = "" #Non-standardized
    dftemp["Sub-Cellular_Locator"] = "" #Non-standardized

    # Result temp columns
    dftemp["Non-Neoplastic_Finding_Type"] = "" #SEND: Non-Neoplastic Finding Type
    dftemp["Neoplasm_Type"] = "" #SEND: Neoplasm Type
    dftemp["Within_Normal_Limits_Results"] = "" #SEND: Within Normal Limits Results
    dftemp["Null_Flavor_Reason"] = "" #SEND: Null Flavor Reason

    # Result final columns
    dftemp["MISTRESC"] = "" #SEND: Result Value
    dftemp["MIRESCAT"] = ""  # SEND: Result Type
    dftemp["MISTAT"] = "" #SEND: Null Value
    dftemp["MIREASND"] = ""  # SEND: Null Reason
    dftemp["MISPCCND"] = "" #SEND: sample condition (ex: autolysis)
    dftemp["MISPCUFL"] = "" #SEND: specimen usability (N=not usable, NULL=usable)

    # Modifier columns
    dftemp["MISEV"] = "" #SEND: Severity
    dftemp["MILAT"] = "" #SEND: Laterality
    dftemp["MIDIR"] = "" #SEND: Specimen Directionality
    dftemp["MIDISTR"] = "" #SEND: Distribution of finding
    dftemp["MICHRON"] = "" #SEND: Chronicity
    dftemp["MIMETHOD"] = "" #SEND: Method used for staining
    dftemp["Directionality_of_Change"] = "" #Non-standardized
    dftemp["Extra_finding_terms"] = "" #Non-standardized

    dftemp["CODE_CONFLICT_FLAG"] = ""
    dftemp = dftemp.fillna("")

    return df, dftemp

###########################################################################
###              Gender, Species, Severity, Route mapping               ###
###########################################################################

def Gender_Mapping(dataframe, dfSENDCodelist):
    # Create dictionary
    #Gender_Codelist = dfSENDCodelist.set_index('Code').to_dict()['CDISC Submission Value']
    Gender_dict = dict(zip(dfSENDCodelist['CDISC Submission Value'], dfSENDCodelist['CDISC Submission Value']))

    Gender_dict[codetermExtract(dfSENDCodelist, "C20197")] = [Gender_dict[codetermExtract(dfSENDCodelist, "C20197")], "MALE", "M"]
    Gender_dict[codetermExtract(dfSENDCodelist, "C16576")] = [Gender_dict[codetermExtract(dfSENDCodelist, "C16576")], "FEMALE", "F"]
    Gender_dict[codetermExtract(dfSENDCodelist, "C17998")] = [Gender_dict[codetermExtract(dfSENDCodelist, "C17998")], "COMBINED", "UNKNOWN", "U", ""]

    reverse_dict = invert_dict(Gender_dict)
    dataframe['SEX'] = dataframe['SEX_raw'].map(reverse_dict)

    return dataframe

def Species_Mapping(dataframe, dfSENDCodelist, dfSTRAINCodelist):
    # Create dictionary
    Species_dict = dict(zip(dfSENDCodelist['CDISC Submission Value'], dfSENDCodelist['CDISC Submission Value']))
    Strain_dict = dict(zip(dfSTRAINCodelist['CDISC Submission Value'], dfSTRAINCodelist['CDISC Submission Value']))

    reverse_dict = invert_dict(Species_dict)

    dataframe['SPECIES'] = dataframe['SPECIES_raw'].map(reverse_dict)

    Strain_dict["BEAGLE"] = "BEAGLE"
    Strain_dict["HARTLEY MALE"] = "HARTLEY"
    Strain_dict["GOLDEN SYRIAN"] = "SYRIAN"
    Strain_dict["CYNOMOLGUS"] = "CYNOMOLGUS"
    Strain_dict["BALB C"] = "BALB/C"
    Strain_dict["C3H HENCRL"] = "C3H/He"
    Strain_dict["C57BL 6"] = "C57BL/6"
    Strain_dict["CB6F1 TGN RASH2"] = "CB6F1-TgN (RasH2)"
    Strain_dict["CD 1"] = "CD1(ICR)"
    Strain_dict["CD 1 ICR"] = "CD1(ICR)"
    Strain_dict["KNOCKOUT"] = "KNOCKOUT" #NON-STANDARD
    Strain_dict["PDAPP"] = "PDAPP" #NON-STANDARD
    Strain_dict["SPRAGUE DAWLEY"] = "SPRAGUE-DAWLEY"
    Strain_dict["TRANSGENIC"] = "TRANSGENIC" #NON-STANDARD
    Strain_dict["GOTTINGEN"] = "GOTTINGEN"
    Strain_dict["NEW ZEALAND WHITE"] = "NEW ZEALAND"
    Strain_dict["CD IGS"] = "SPRAGUE-DAWLEY"
    Strain_dict["CD IGS SPRAGUE DAWLEY"] = "SPRAGUE-DAWLEY"
    Strain_dict["FISCHER 344"] = "FISCHER 344"
    Strain_dict["LONG EVANS"] = "LONG EVANS"
    Strain_dict["SPRAGUE DAWLEY"] = "SPRAGUE-DAWLEY"
    Strain_dict["UNKNOWN"] = ""
    Strain_dict[""] = ""

    reverse_dict = invert_dict(Strain_dict)

    dataframe['STRAIN'] = dataframe['STRAIN_raw'].map(reverse_dict)

    return dataframe

def Severity_Mapping(dataframe, dfSENDCodelist):
    # Create dictionary
    Severity_dict = dict(zip(dfSENDCodelist['CDISC Submission Value'], dfSENDCodelist['CDISC Submission Value']))

    Severity_dict[codetermExtract(dfSENDCodelist, "C147501")] = [Severity_dict[codetermExtract(dfSENDCodelist, "C147501")],"1 OF 5", "MINIMAL"]
    Severity_dict[codetermExtract(dfSENDCodelist, "C147504")] = [Severity_dict[codetermExtract(dfSENDCodelist, "C147504")],"2 OF 5", "SLIGHT", "MILD"]
    Severity_dict[codetermExtract(dfSENDCodelist, "C147507")] = [Severity_dict[codetermExtract(dfSENDCodelist, "C147507")],"3 OF 5", "MODERATE"]
    Severity_dict[codetermExtract(dfSENDCodelist, "C147509")] = [Severity_dict[codetermExtract(dfSENDCodelist, "C147509")],"4 OF 5", "MARKED"]
    Severity_dict[codetermExtract(dfSENDCodelist, "C147510")] = [Severity_dict[codetermExtract(dfSENDCodelist, "C147510")],"5 OF 5", "SEVERE"]
    Severity_dict["NO GRADE ASSIGNED"] = ["PRESENT", "PRESENT NO GRADE ASSIGNED", "NO GRADE ASSIGNED"]  # extensible

    reverse_dict = invert_dict(Severity_dict)

    dataframe['MISEV'] = dataframe['MISEV_raw'].map(reverse_dict)

    return dataframe

def Route_Mapping(dataframe, dfSENDCodelist):
    # Create dictionary
    Route_dict = dict(zip(dfSENDCodelist['CDISC Submission Value'], dfSENDCodelist['CDISC Submission Value']))

    Route_dict[codetermExtract(dfSENDCodelist, "C78373")] = [Route_dict[codetermExtract(dfSENDCodelist, "C78373")],"WATER", "CAPSULE"] #Is "CAPSULE" correct??
    Route_dict[codetermExtract(dfSENDCodelist, "C78374")] = [Route_dict[codetermExtract(dfSENDCodelist, "C78374")],"GAVAGE", "OROGASTRIC", "NASOGASTRIC"]
    Route_dict[codetermExtract(dfSENDCodelist, "C38288")] = [Route_dict[codetermExtract(dfSENDCodelist, "C38288")],"ORAL"]
    Route_dict[codetermExtract(dfSENDCodelist, "C38299")] = [Route_dict[codetermExtract(dfSENDCodelist, "C38299")],"SUBCUTANEOUS"]
    Route_dict[codetermExtract(dfSENDCodelist, "C38311")] = [Route_dict[codetermExtract(dfSENDCodelist, "C38311")],"SEE PROTOCOL", ""]
    Route_dict[codetermExtract(dfSENDCodelist, "C48623")] = [Route_dict[codetermExtract(dfSENDCodelist, "C48623")],"NOT DOSED", "NOT APPLICABLE"]
    Route_dict[codetermExtract(dfSENDCodelist, "C38279")] = [Route_dict[codetermExtract(dfSENDCodelist, "C38279")],"INFUSION", "INTRAVENOUS INFUSION"]
    Route_dict[codetermExtract(dfSENDCodelist, "C38276")] = [Route_dict[codetermExtract(dfSENDCodelist, "C38276")],"INTRAVENOUS"]
    Route_dict[codetermExtract(dfSENDCodelist, "C38274")] = [Route_dict[codetermExtract(dfSENDCodelist, "C38274")],"INTRAVENOUS BOLUS", "INTRAVENOUS BOLUS INFUSION"]
    Route_dict[codetermExtract(dfSENDCodelist, "C38216")] = [Route_dict[codetermExtract(dfSENDCodelist, "C38216")],"INHALATION"]
    Route_dict[codetermExtract(dfSENDCodelist, "C102399")] = [Route_dict[codetermExtract(dfSENDCodelist, "C102399")],"INTRAJEJUNAL"]
    Route_dict[codetermExtract(dfSENDCodelist, "C28161")] = [Route_dict[codetermExtract(dfSENDCodelist, "C28161")],"INTRAMUSCULAR"]
    Route_dict[codetermExtract(dfSENDCodelist, "C38258")] = [Route_dict[codetermExtract(dfSENDCodelist, "C38258")],"INTRAPERITONEAL"]

    reverse_dict = invert_dict(Route_dict)

    dataframe['EXROUTE'] = dataframe['EXROUTE_raw'].map(reverse_dict)

    return dataframe

def Disposition_Mapping(dataframe, dfSENDCodelist):
    # Create dictionary
    Disposition_dict = dict(zip(dfSENDCodelist['CDISC Submission Value'], dfSENDCodelist['CDISC Submission Value']))


    Disposition_dict[codetermExtract(dfSENDCodelist, "C90351")] = [Disposition_dict[codetermExtract(dfSENDCodelist, "C90351")], "ACCIDENTAL DEATH", "ACCIDENTAL", "UNPLANNED KILLED TERMINAL"]
    Disposition_dict[codetermExtract(dfSENDCodelist, "C90425")] = [Disposition_dict[codetermExtract(dfSENDCodelist, "C90425")], "KILLED CLINICAL", "KILLED MORIBUND"]
    Disposition_dict[codetermExtract(dfSENDCodelist, "C90436")] = [Disposition_dict[codetermExtract(dfSENDCodelist, "C90436")], "KILLED INTERIM", "INTERIM SACRIFICE"]
    Disposition_dict[codetermExtract(dfSENDCodelist, "C123635")] = [Disposition_dict[codetermExtract(dfSENDCodelist, "C123635")], "NON MORIBUND SACRIFICE"]
    Disposition_dict[codetermExtract(dfSENDCodelist, "C90445")] = [Disposition_dict[codetermExtract(dfSENDCodelist, "C90445")], "RECOVERY SACRIFICE"]
    Disposition_dict[codetermExtract(dfSENDCodelist, "C90465")] = [Disposition_dict[codetermExtract(dfSENDCodelist, "C90465")], "KILLED TERMINAL", "TERMINAL SACRIFICE"]
    Disposition_dict[codetermExtract(dfSENDCodelist, "C90387")] = [Disposition_dict[codetermExtract(dfSENDCodelist, "C90387")], "DIED"]
    Disposition_dict[codetermExtract(dfSENDCodelist, "C90447")] = [Disposition_dict[codetermExtract(dfSENDCodelist, "C90447")], "TREATMENT GROUP ABORTED"]
    Disposition_dict[codetermExtract(dfSENDCodelist, "C96372")] = [Disposition_dict[codetermExtract(dfSENDCodelist, "C96372")], "OTHER SEE TEXT COMMENTS", ""]

    reverse_dict = invert_dict(Disposition_dict)

    dataframe['DSDECOD'] = dataframe['DSDECOD_raw'].map(reverse_dict)

    return dataframe

###########################################################################
###                        Specimen functions                           ###
###########################################################################

#Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12744")) #

##Core Specimen function
def Specimen_Mapping(row, dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist):

    rowtemp = row[['MISPEC_raw', 'MIANTREG_raw']].copy()
    rowSpecimen = rowtemp.values.tolist()

    if (Search_Term(rowSpecimen, ["ADRENAL"])):
        row = Adrenal_Gland_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["NERVE", "GANGLI"])):
        row = Nervous_System_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["SPINAL"])):
        row = Spinal_Cord_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["MARROW", "BUFFY"])):
        row = Hematopoietic_System_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["JOINT", "ANKLE", "KNEE"])):
        row = Joint_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["BONE", "VERTEBRA", "SKULL", "TIBIA", "FEMUR", "MAXILLA", "MANDIBLE"])):
        row = Skeletal_System_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    ##NUCLEI may be an issue
    elif (Search_Term(rowSpecimen, ["BRAIN", "MENINGE", "SUBSTANTIA NIGRA", "NUCLEI", "COLLICULUS", "CINGULATE", "COERULEUS", "PIRIFORM", "PITUITARY"])):
        row = Brain_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist)
    elif (Search_Term(rowSpecimen, ["HEART", "VALVE", "PERICARD"])):
        row = Heart_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["VEIN", "TUNICA", "VESSEL", "VASCULAR", "ARTERY", "MESENTER", "AORT"])):
        row = Vascular_System_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["STOMACH", "PYLORUS"])):
        row = Stomach_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["INTESTIN","ANUS","RECTUM","PEYER"])):
        row = Intestine_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["KIDNEY","URINARY","URETER","URETHRA"])):
        row = Kidney_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["LIVER","GALLBLADDER","BILE"])):
        row = Hepatobiliary_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["LUNG","PLEURA"])):
        row = Lung_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["LYMPH"])):
        row = Lymphoid_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["MUSCLE","DIAPHRAGM"])):
        row = Musculature_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["PANCREAS"])):
        row = Pancreas_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["SPLEEN", "SPLENIC"])):
        row = Spleen_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["THYMUS"])):
        row = Thymus_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["THYROID"])):
        row = Thyroid_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["MAMMARY","VAGINA","OVARY","UTERUS", "UTERINE", "CERVIX", "OVIDUCT","REPRODUCTIVE","ESTRUS"])):
        row = Female_Repro_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    ##"EFFERENT" may be an issue
    elif (Search_Term(rowSpecimen, ["TESTIS","VAS DEFERENS","EPIDIDYM","SPERMAT","PROSTATE","PENIS","SCROTUM","SEMINAL","EFFERENT"])):
        row = Male_Repro_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
        ##"GLAND" may be an issue
    elif (Search_Term(rowSpecimen, ["GLAND"])):
        row = Gland_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
        ##"SITE" may be an issue
    elif (Search_Term(rowSpecimen, ["SITE"])):
        row = Injection_Application_Site_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
        ##"LIP" may be an issue
    elif (Search_Term(rowSpecimen, ["SKIN","LIP","DERMIS","HINDPAW"])):
        row = Integument_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["EYE","ORBIT","CONJUNCTIVA","LACRIMAL"])):
        row = Eye_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist)
        ##"FAT" may be an issue
    elif (Search_Term(rowSpecimen, ["FAT","ADIPOSE"])):
        row = Adipose_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
        ##"ORAL" may be an issue
    elif (Search_Term(rowSpecimen, ["SALIVARY", "BODY CAVITY ORAL", "GINGIVA", "TONGUE", "TOOTH"])):
        row = Mouth_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["NASAL", "NOSE", "TURBINATE", "OLFACTORY", "NASOLACRIMAL"])):
        row = Nose_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    ##"MULTIPLE" and "CAVITY" may be an issue
    elif (Search_Term(rowSpecimen, ["ENDOTHELIUM","MULTIPLE","WHOLE BODY","MEDIASTINUM","PERITONE","CAVITY"])):
        row = Whole_Body_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)

    return row

##Specimen-specific functions
def Adrenal_Gland_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # ADRENAL GLAND #UBERON:0002369

    row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12666") #adrenal gland
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12666") #adrenal gland

    if (Search_Term(rowSpecimen, ["CORTEX"])):
        row['Anatomical_Region_of_Specimen'] = "CORTEX"
        if ((Search_Term(rowSpecimen, ["MEDULLA"])) and not (Search_Term(rowSpecimen, ["EXTRAMEDULLARY"]))):
            if (Search_Term(rowSpecimen, ["JUNCTION"])):
                row['Anatomical_Region_of_Specimen'] = "CORTICOMEDULLARY JUNCTION"
            else:
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MEDULLA")
    elif ((Search_Term(rowSpecimen, ["MEDULLA"])) and not (Search_Term(rowSpecimen, ["EXTRAMEDULLARY"]))):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MEDULLA")
    elif (Search_Term(rowSpecimen, ["CAPSULE"])):
        row['Anatomical_Region_of_Specimen'] = "CAPSULE"

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["X ZONE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CORTEX")
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "X-ZONE")

    if (Search_Term(rowSpecimen, ["FASCICUL", "FASICULATA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CORTEX")
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ZONA FASCICULATA")

    if (Search_Term(rowSpecimen, ["GLOMERULOSA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CORTEX")
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ZONA GLOMERULOSA")

    if (Search_Term(rowSpecimen, ["RETICULARIS"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CORTEX")
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ZONA RETICULARIS")

    if (Search_Term(rowSpecimen, ["INTERMEDIA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CORTEX")
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ZONA INTERMEDIA")

    ##CELLTYPES
    if (Search_Term(row, ["CHROMAFFIN CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "CHROMAFFIN CELL")
    if (Search_Term(row, ["ISLET CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "ISLET CELL")
    if (Search_Term(row, ["NEUROENDOCRINE"])):
        row['Celltype'] = Add_Term('Celltype', row, "NEUROENDOCRINE CELL")

    return row

def Nervous_System_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # NERVOUS SYSTEM #UBERON:0001016
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12768")  # PERIPHERAL NERVE

    if (Search_Term(rowSpecimen, ["GANGLI"])):
        if (Search_Term(rowSpecimen, ["CERVICOTHORACIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92211")  # GANGLION, CERVICOTHORACIC
        elif (Search_Term(rowSpecimen, ["CERVICAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C98713")  # GANGLION, CERVICAL
        elif (Search_Term(rowSpecimen, ["DORSAL ROOT"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12462")  # GANGLION, DORSAL ROOT
            if (Search_Term(rowSpecimen, ["LUMBAR"])):
                row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12744") # LUMBAR VERTEBRA
            elif (Search_Term(rowSpecimen, ["CERVICAL"])):
                row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12693") # CERVICAL VERTEBRA
            elif (Search_Term(rowSpecimen, ["THORACIC"])):
                row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12798") # THORACIC VERTEBRA
            elif (Search_Term(rowSpecimen, ["SACRAL"])):
                row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12853") # SACRAL VERTEBRA
        elif (Search_Term(rowSpecimen, ["TRIGEMINAL"])):
            if (Search_Term(rowSpecimen, ["NERVE"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92214")  # GANGLION, TRIGEMINAL/NERVE, TRIGEMINAL
            else:
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C62642")  # GANGLION, TRIGEMINAL
        elif (Search_Term(rowSpecimen, ["BASAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12447")  # BRAIN, BASAL GANGLIA
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12719")  # GANGLION
        if (Search_Term(rowSpecimen, ["SYMPATHETIC"])):
            if (Search_Term(rowSpecimen, ["PARASYMPATHETIC"])):
                row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C52557")  # PARASYMPATHETIC GANGLIA
            else:
                row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12467")  # SYMPATHETIC GANGLIA
        if (Search_Term(rowSpecimen, ["PARAVERTEBRAL"])):
            row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C97925")  # PARAVERTEBRAL GANGLIA
    elif (Search_Term(rowSpecimen, ["NERVE"])):
        if (Search_Term(rowSpecimen, ["BRACHIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12682")  # NERVE, BRACHIAL PLEXUS
            #row['UBERON_Code'] = "UBERON:0003433"  # arm nerve
        elif (Search_Term(rowSpecimen, ["OPTIC"])):
            if (Search_Term(rowSpecimen, ["DISC", "NERVE HEAD"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12760")  # OPTIC DISC
            else:
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12761")  # NERVE, OPTIC
                # row['UBERON_Code'] = "UBERON:0000941"  # cranial nerve II
        elif (Search_Term(rowSpecimen, ["CRANIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12700")  # NERVE, CRANIAL
        elif (Search_Term(rowSpecimen, ["SCIATIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52810")  # NERVE, SCIATIC
        elif (Search_Term(rowSpecimen, ["SURAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77675")  # NERVE, SURAL
        elif (Search_Term(rowSpecimen, ["TIBIA"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52809")  # NERVE, TIBIAL
        elif (Search_Term(rowSpecimen, ["ULNA"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52807")  # NERVE, ULNAR
        elif (Search_Term(rowSpecimen, ["TRIGEMINAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12806")  # NERVE, TRIGEMINAL
        elif (Search_Term(rowSpecimen, ["FACIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12714")  # NERVE, FACIAL
        elif (Search_Term(rowSpecimen, ["VAGUS","VAGAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12812")  # NERVE, VAGUS
        elif (Search_Term(rowSpecimen, ["MEDIAN"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52815")  # NERVE, MEDIAN
        elif (Search_Term(rowSpecimen, ["ROOT"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C54024")  # NERVE ROOT
        elif (Search_Term(rowSpecimen, ["SAPHENOUS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C147511")  # NERVE, SAPHENOUS
        elif (Search_Term(rowSpecimen, ["PERONEAL"])):
            if (Search_Term(rowSpecimen, ["COMMON"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92601")  # NERVE, PERONEAL, COMMON
            elif (Search_Term(rowSpecimen, ["DEEP"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92602")  # NERVE, PERONEAL, DEEP
            elif (Search_Term(rowSpecimen, ["SUPERFICIAL"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92603")  # NERVE, PERONEAL, SUPERFICIAL
            else:
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52814")  # NERVE, PERONEAL
        elif (Search_Term(rowSpecimen, ["PLANTAR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77674")  # NERVE, PLANTAR
        elif (Search_Term(rowSpecimen, ["RADIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52812")  # NERVE, RADIAL
        elif (Search_Term(rowSpecimen, ["SPINAL"])):
            if (Search_Term(rowSpecimen, ["SPINAL ACCESSORY"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12700")  # NERVE, CRANIAL
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32041"))  # SPINAL ACCESSORY NERVE
            else:
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12792")  # NERVE, SPINAL
        elif (Search_Term(rowSpecimen, ["PERIPHERAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12768")  # NERVE, PERIPHERAL
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12466")  # NERVE
        if (Search_Term(rowSpecimen, ["ABDUCENS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12700")  # NERVE, CRANIAL
            row['MIANTREG'] = Add_Term('MIANTREG', row,
                                       codetermExtract(dfAnatLocCodelist, "C12665"))  # ABDUCENS NERVE
        if (Search_Term(rowSpecimen, ["GLOSSOPHARYNGEAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12700")  # NERVE, CRANIAL
            row['MIANTREG'] = Add_Term('MIANTREG', row,
                                       codetermExtract(dfAnatLocCodelist, "C12723"))  # GLOSSOPHARYNGEAL NERVE
        if (Search_Term(rowSpecimen, ["HYPOGLOSSAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12700")  # NERVE, CRANIAL
            row['MIANTREG'] = Add_Term('MIANTREG', row,
                                       codetermExtract(dfAnatLocCodelist, "C12732"))  # HYPOGLOSSAL NERVE
        if (Search_Term(rowSpecimen, ["OCULOMOTOR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12700")  # NERVE, CRANIAL
            row['MIANTREG'] = Add_Term('MIANTREG', row,
                                       codetermExtract(dfAnatLocCodelist, "C12758"))  # OCULOMOTOR NERVE
        if (Search_Term(rowSpecimen, ["OLFACTORY"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12700")  # NERVE, CRANIAL
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12759"))  # OLFACTORY NERVE
        if (Search_Term(rowSpecimen, ["TROCHLEAR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12700")  # NERVE, CRANIAL
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12808"))  # TROCHLEAR NERVE
        if (Search_Term(rowSpecimen, ["VESTIBULOCOCHLEAR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12700")  # NERVE, CRANIAL
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12996"))  # VESTIBULOCOCHLEAR NERVE

    if (Search_Term(rowSpecimen, ["EPINEURI"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "EPINEURIUM")
    if (Search_Term(rowSpecimen, ["PERINEUR"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PERINEURAL")
    if (Search_Term(row, ["ADJACENT"])): #or row
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ADJACENT TISSUE")
    if (Search_Term(rowSpecimen, ["RETINAL NERVE FIBER LAYER"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32953"))  # RETINAL NERVE FIBER LAYER

    ##CELLTYPES
    if (Search_Term(row, ["NEURON", "NEURODEGEN"])):
        row['Celltype'] = Add_Term('Celltype', row, "NEURON")

    return row

def Spinal_Cord_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # SPINAL CORD #UBERON:0002240
    # spinal column--vertebra
    row['Sample_Type_Group'] = codetermExtract(dfAnatLocCodelist, "C12438")  # CENTRAL NERVOUS SYSTEM
    if (Search_Term(rowSpecimen, ["LUMBAR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12895")  # SPINAL CORD, LUMBAR
    elif (Search_Term(rowSpecimen, ["CERVICAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12892")  # SPINAL CORD, CERVICAL
    elif (Search_Term(rowSpecimen, ["THORACIC"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12894")  # SPINAL CORD, THORACIC
    elif (Search_Term(rowSpecimen, ["SACRAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12896")  # SPINAL CORD, SACRAL
    elif (Search_Term(rowSpecimen, ["COLUMN"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12998")  # SPINAL COLUMN
    elif (Search_Term(rowSpecimen, ["FLUID"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12692")  # FLUID, CEREBROSPINAL
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12464")  # SPINAL CORD

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["DURA MATER"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "DURA MATER")
    if (Search_Term(rowSpecimen, ["CENTRAL CANAL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "CENTRAL CANAL")
    if (Search_Term(rowSpecimen, ["GRAY MATTER"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GRAY MATTER")
    if (Search_Term(rowSpecimen, ["WHITE MATTER"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "WHITE MATTER")
    if (Search_Term(rowSpecimen, ["FUNICUL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "FUNICULUS")
    if (Search_Term(rowSpecimen, ["MENINGES"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12348")  # MENINGES

    ##CELLTYPES
    if (Search_Term(row, ["ASTROCYTE"])):
        row['Celltype'] = Add_Term('Celltype', row, "ASTROCYTE")
    if (Search_Term(row, ["NEURON", "NEURODEGEN"])):
        row['Celltype'] = Add_Term('Celltype', row, "NEURON")
    if (Search_Term(row, ["MICROGLIOSIS", "MICROGLIA"])):
        row['Celltype'] = Add_Term('Celltype', row, "MICROGLIA")
    if (Search_Term(row, ["GLIOSIS", "GLIAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "GLIAL CELL")

    return row

def Hematopoietic_System_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # HEMATOPOIETIC SYSTEM #UBERON:0002390
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12431") # BONE MARROW

    if (Search_Term(rowSpecimen, ["BUFFY COAT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C84507") # BUFFY COAT
        row['Sample_Type_Group'] = "OTHER" # blood
    else:
        if (Search_Term(rowSpecimen, ["FEMUR"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77686")) # BONE MARROW, FEMUR
        elif (Search_Term(rowSpecimen, ["STERNUM"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77690")) # BONE MARROW, STERNUM
        elif (Search_Term(rowSpecimen, ["TIBIA"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77691")) # BONE MARROW, TIBIA
        elif (Search_Term(rowSpecimen, ["HUMERUS"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77687")) # BONE MARROW, HUMERUS
        elif (Search_Term(rowSpecimen, ["SCAPULA"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77689")) # BONE MARROW, SCAPULA
        elif (Search_Term(rowSpecimen, ["VERTEBR"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77692")) # BONE MARROW, VERTEBRUM
        elif (Search_Term(rowSpecimen, ["RIB"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77688")) # BONE MARROW, RIB
        else:
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12431")) # BONE MARROW

    ###EXTRA TERMS
    if (Search_Term(row, ["METAPHYSIS"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "METAPHYSIS")
    if (Search_Term(row, ["ERYTHROID"])):
        row['Celltype'] = Add_Term('Celltype', row, "ERYTHROCYTE")
        if (Search_Term(row, ["RATIO"])):
            pass
        else:
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ERYTHROID")
    if (Search_Term(row, ["MYELOID"])):
        row['Celltype'] = Add_Term('Celltype', row, "MYELOID CELL")
        if (Search_Term(row, ["LEFT SHIFT"])):
            row['MISTRESC'] = Add_Term('MISTRESC', row, "MYELOID:ERYTHROID RATIO, INCREASED") # EXTENSIBLE
        elif (Search_Term(row, ["RATIO"])):
            if (("INCREAS" in row.at['MIORRES']) or ("INCREAS" in row.at['MISTRESC'])):
                row['MISTRESC'] = Add_Term('MISTRESC', row, "MYELOID:ERYTHROID RATIO, INCREASED") # EXTENSIBLE
            elif (("DECREAS" in row.at['MIORRES']) or ("DECREAS" in row.at['MISTRESC'])):
                row['MISTRESC'] = Add_Term('MISTRESC', row, "MYELOID:ERYTHROID RATIO, DECREASED") # EXTENSIBLE
            else:
                row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MISTRESC_non-match")
        else:
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MYELOID")
    # Non_neoplastic_finding_dict["INCREASED ME RATIO"] = "MYELOID:ERYTHROID RATIO, INCREASED"

    ##CELLTYPES
    if (Search_Term(row, ["MYELOGENIC PRECURSOR", "MYELOID PRECURSOR"])):
        row['Celltype'] = Add_Term('Celltype', row, "MYELOID PRECURSOR CELL")
    if (Search_Term(row, ["MYELOPOIE"])):
        row['Celltype'] = Add_Term('Celltype', row, "MYELOID CELL")
    if (Search_Term(row, ["MEGAKARYOCYT", "MEGAKERYOCYT"])):
        row['Celltype'] = Add_Term('Celltype', row, "MEGAKARYOCYTE")

    return row

def Joint_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # JOINT UBERON:0000982
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C13044") # JOINT

    if (Search_Term(rowSpecimen, ["CARPUS","CARPAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32264") # JOINT, CARPUS
        #row['UBERON_Code'] = "UBERON:0011132"  # intercarpal joint
    elif (Search_Term(rowSpecimen, ["ELBOW"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32497") # JOINT, ELBOW
        if (Search_Term(rowSpecimen, ["FLEXOR"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C163513")) # ELBOW FLEXOR MUSCLES
        elif (Search_Term(rowSpecimen, ["EXTENSOR"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C163512")) # ELBOW ENTENSOR MUSCLES
    elif (Search_Term(rowSpecimen, ["FEMOROTIBIA", "PATELLA", "KNEE", "TIBIOFEM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32898") # JOINT, FEMOROTIBIAL
    elif (Search_Term(rowSpecimen, ["TARSUS","TARSAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33735") # JOINT, TARSUS
        #row['UBERON_Code'] = "UBERON:0008447"  # intertarsal joint
    elif (Search_Term(rowSpecimen, ["SCAPUL", "HUMER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C111308") # JOINT, SCAPULOHUMERAL
    elif (Search_Term(rowSpecimen, ["HIP", "COXOFEM" ])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32742") # JOINT, HIP
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C13044") # JOINT

    ##EXTRA TERMS
    if (Search_Term(rowSpecimen, ["BICEPS"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C32200")) # MUSCLE, BICEPS BRACHII
    if (Search_Term(rowSpecimen, ["TRICEPS"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C90604")) # MUSCLE, TRICEPS BRACHII
    if (Search_Term(rowSpecimen, ["PALM"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33252")) # PALM
    if (Search_Term(rowSpecimen, ["PLANTAR", "SOLE"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33326")) # SOLE
    if (Search_Term(rowSpecimen, ["FOOTPAD"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C92654")) # FOOTPAD
    if (Search_Term(rowSpecimen, ["ANKLE"])):
        if (Search_Term(rowSpecimen, ["FLEXOR TENDOR"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C161389"))  # ANKLE JOINT ANTERIOR FLEXOR TENDONS
        elif (Search_Term(rowSpecimen, ["EXTENSOR TENDON"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C161390"))  # ANKLE JOINT ANTERIOR EXTENSOR TENDONS

    return row

def Skeletal_System_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # SKELETAL SYSTEM #UBERON:0001434
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12366") # BONE
    if (Search_Term(rowSpecimen, ["MARROW", "ERYTHROID", "MYELOID", "MEDULLA", "POIETIC", "POIESIS"])):
        row = Hematopoietic_System_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
    elif (Search_Term(rowSpecimen, ["VERTEBRA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33868") # BONE, VERTEBRA
        if (Search_Term(rowSpecimen, ["LUMBAR"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12744")) # LUMBAR VERTEBRA
        elif (Search_Term(rowSpecimen, ["CERVICAL"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12693")) # CERVICAL VERTEBRA
        elif (Search_Term(rowSpecimen, ["THORACIC"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12798")) # THORACIC VERTEBRA
        elif (Search_Term(rowSpecimen, ["SACRAL"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12853")) # SACRAL VERTEBRA
    elif (Search_Term(rowSpecimen, ["CALVARIUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C81188") # BONE, CALVARIUM
    elif (Search_Term(rowSpecimen, ["CARPAL", "CARPUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12688") #BONE, CARPAL
    elif (Search_Term(rowSpecimen, ["CLAVICLE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12695") # BONE, CLAVICLE
    elif (Search_Term(rowSpecimen, ["CONDYLE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C83002") #BONE, CONDYLE
    elif (Search_Term(rowSpecimen, ["FEMUR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12717") # BONE, FEMUR
    elif (Search_Term(rowSpecimen, ["FIBULA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12718") # BONE, FIBULA
    elif (Search_Term(rowSpecimen, ["HUMERUS", "HUMERAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12731") # BONE, HUMERUS
    elif (Search_Term(rowSpecimen, ["MANDIB", "LOWER JAW", "INFERIOR MAXILLA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12290") # BONE, MANDIBLE
    elif (Search_Term(rowSpecimen, ["MAXILLA", "UPPER JAW"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C26470") # BONE, MAXILLA
    elif (Search_Term(rowSpecimen, ["PATELLA", "KNEE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33282") # BONE, PATELLA
    elif (Search_Term(rowSpecimen, ["PELVIS", "PELVIC"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33287") # BONE, PELVIS
    elif (Search_Term(rowSpecimen, ["RADIUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12777") # BONE, RADIUS
    elif (Search_Term(rowSpecimen, ["SCAPULA", "SHOULDER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12783") # BONE, SCAPULA
    elif (Search_Term(rowSpecimen, ["SKULL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12789") # BONE, SKULL
        if (Search_Term(row, ["SCALLOPED REVERSAL LINE"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "CORTICAL BONE")
    elif (Search_Term(rowSpecimen, ["STERNUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12793") # BONE, STERNUM
    elif (Search_Term(rowSpecimen, ["TARSAL","TARSUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12796") # BONE, TARSUS
    elif (Search_Term(rowSpecimen, ["TIBIA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12800") # BONE, TIBIA
    elif (Search_Term(rowSpecimen, ["ULNA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12809") # BONE, ULNA
    elif (Search_Term(rowSpecimen, ["RIB"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12782") # BONE, RIB
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12366") # BONE

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["PHYSIS", "PHYSEAL", "GROWTH PLATE", "EPIPHYSEAL PLATE"])):
        if (Search_Term(rowSpecimen, ["METAPHYSIS", "METAPHYSEAL"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "METAPHYSIS")
            if (Search_Term(rowSpecimen, ["EPIPHYSIS"])):
                row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C32529")) # EPIPHYSIS
        elif (Search_Term(rowSpecimen, ["EPIPHYSIS"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C32529")) # EPIPHYSIS
        elif (Search_Term(rowSpecimen, ["DIAPHYSIS"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "DIAPHYSIS")
        else:
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GROWTH PLATE")
    if (Search_Term(rowSpecimen, ["ENDOCHONDRAL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GROWTH PLATE")
        # row['Anatomical_Region_of_Specimen']= Add_Term('Anatomical_Region_of_Specimen', row, "ENDOCHONDRAL")
    if (Search_Term(rowSpecimen, ["ENDOSTE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ENDOSTEUM")
    if (Search_Term(rowSpecimen, ["HYPERTROPHIC ZONE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "HYPERTROPHIC ZONE")
    if (Search_Term(row, ["WOVEN BONE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "WOVEN BONE")
    if (Search_Term(rowSpecimen, ["PRIMARY SPONGIOSA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "METAPHYSIS")
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PRIMARY SPONGIOSA")
    if (Search_Term(rowSpecimen, ["ARTICULAR"])):
        if (Search_Term(rowSpecimen, ["CARTILAGE"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ARTICULAR CARTILAGE")
        elif (Search_Term(rowSpecimen, ["PERIARTICULAR"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PERIARTICULAR")
        else:
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C13044")) # JOINT
    if (Search_Term(rowSpecimen, ["PERIOSTEUM"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PERIOSTEUM")
    if (Search_Term(rowSpecimen, ["OSTEOID"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "OSTEOID")
    if (Search_Term(rowSpecimen, ["CHONDROID"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CHONDROID")
    if (Search_Term(rowSpecimen, ["T11"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33722")) # T11 VERTEBRA
    elif (Search_Term(rowSpecimen, ["L1"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32899")) # L1 VERTEBRA
    elif (Search_Term(rowSpecimen, ["L2"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32900")) # L2 VERTEBRA
    if (Search_Term(row, ["INTERVERTEBRAL DIS"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C49571")) # INTERVERTEBRAL DISC
    if (Search_Term(rowSpecimen, ["TRABECULA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "TRABECULA")
        #row['UBERON_Code'] = "UBERON:0002483"  # trabecular bone tissue
    if (Search_Term(rowSpecimen, ["VERTEBRAL BODY"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "VERTEBRAL BODY")

    ##CELLTYPES
    if (Search_Term(row, ["CHONDROCYTE"])):
        row['Celltype'] = Add_Term('Celltype', row, "CHONDROCYTE")
    if (Search_Term(row, ["OSTEOBLAST"])):
        row['Celltype'] = Add_Term('Celltype', row, "OSTEOBLAST")
    if (Search_Term(row, ["OSTEOCLAST"])):
        row['Celltype'] = Add_Term('Celltype', row, "OSTEOCLAST")
    if (Search_Term(row, ["MYELOGENIC PRECURSOR", "MYELOID PRECURSOR"])):
        row['Celltype'] = Add_Term('Celltype', row, "MYELOID PRECURSOR CELL")
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12431") # BONE MARROW
    if (Search_Term(row, ["MYELOPOIE"])):
        row['Celltype'] = Add_Term('Celltype', row, "MYELOID CELL")
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12431") # BONE MARROW
    if (Search_Term(row, ["MEGAKARYOCYT", "MEGAKERYOCYT"])):
        row['Celltype'] = Add_Term('Celltype', row, "MEGAKARYOCYTE")
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12431") # BONE MARROW

    return row

def Brain_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist):

    # BRAIN #UBERON:0000955
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12439") # BRAIN

    if (Search_Term(rowSpecimen, ["OBEX"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92592") # BRAIN, OBEX
    elif (Search_Term(rowSpecimen, ["OLFACTORY BULB"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C28401") # BRAIN, OLFACTORY BULB
    elif (Search_Term(rowSpecimen, ["DEEP CEREBELLAR NUCLE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12445") # BRAIN, CEREBELLUM
        row['Anatomical_Region_of_Specimen'] = "DEEP CEREBELLAR NUCLEI"
    elif (Search_Term(rowSpecimen, ["LATERAL SEPTAL NUCLE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12447") # BRAIN, BASAL GANGLIA
        row['Anatomical_Region_of_Specimen'] = "LATERAL SEPTAL NUCLEI"
    elif (Search_Term(rowSpecimen, ["DORSAL MOTOR NUCLE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12442") # BRAIN, MEDULLA OBLONGATA
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12934")) # DORSAL MOTOR NUCLEUS
    elif (Search_Term(rowSpecimen, ["ROSTRAL COLLICULUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12510") # BRAIN, MIDBRAIN
        row['Anatomical_Region_of_Specimen'] = "SUPERIOR COLLICULUS"
    elif (Search_Term(rowSpecimen, ["LOCUS COERULEUS", "LOCUS CERULEUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12511") # BRAIN, PONS
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C97333")) # LOCUS CERULEUS
    elif (Search_Term(rowSpecimen, ["PONS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12511") # BRAIN, PONS
        if (Search_Term(rowSpecimen, ["VAROLII"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12511")) # PONS VAROLII
    elif (Search_Term(rowSpecimen, ["SUBSTANTIA NIGRA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12453") # BRAIN, SUBSTANTIA NIGRA
    elif (Search_Term(rowSpecimen, ["CORPUS CALLOSUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12446") # BRAIN, CORPUS CALLOSUM
        if (Search_Term(rowSpecimen, ["BODY"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32216")) # CORPUS CALLOSUM, BODY
        elif (Search_Term(rowSpecimen, ["GENU"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32675")) # CORPUS CALLOSUM, GENU
        elif (Search_Term(rowSpecimen, ["SPLENIUM"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33610")) # CORPUS CALLOSUM, SPLENIUM
    elif (Search_Term(rowSpecimen, ["BASAL GANGLIA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12447") # BRAIN, BASAL GANGLIA
    elif (Search_Term(rowSpecimen, ["VISUAL CORTEX"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C97340") # BRAIN, VISUAL CORTEX
        if (Search_Term(rowSpecimen, ["PRIMARY"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C97340")) # PRIMARY VISUAL CORTEX
    elif (Search_Term(rowSpecimen, ["AMYGDAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12440") # BRAIN, AMYGDALOID BODY
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12440")) # AMYGDALA
    elif (Search_Term(rowSpecimen, ["CEREBR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12351") # BRAIN, CEREBRUM
        #row['UBERON_Code'] = "UBERON:0001893"  # telencephalon
        if (Search_Term(rowSpecimen, ["CORTEX", "CORTICAL"])):
            if (Search_Term(rowSpecimen, ["SUBCORTEX", "SUBCORTICAL"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C98712"))  # CEREBRAL SUBCORTEX
            else:
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12443")) # CEREBRAL CORTEX
        elif (Search_Term(rowSpecimen, ["HIPPOCAMPUS"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12444")) # HIPPOCAMPUS
        elif (Search_Term(rowSpecimen, ["HEMISPHERE"])):
            if (Search_Term(rowSpecimen, ["LEFT"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32955")) # CEREBRAL HEMISPHERE, LEFT
            elif (Search_Term(rowSpecimen, ["RIGHT"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33472"))  # CEREBRAL HEMISPHERE, RIGHT
            elif (Search_Term(row, ["LEFT"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32955")) # CEREBRAL HEMISPHERE, LEFT
            elif (Search_Term(row, ["RIGHT"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33472"))  # CEREBRAL HEMISPHERE, RIGHT
        elif (Search_Term(rowSpecimen, ["ARTERY"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12691")) # CEREBRAL ARTERY
        elif (Search_Term(rowSpecimen, ["VEIN"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C53037")) # CEREBRAL VEIN
        if (Search_Term(rowSpecimen, ["EPENDYM"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12694")) # BRAIN, CHOROID PLEXUS
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "EPENDYMA")
        if (Search_Term(rowSpecimen, ["PIRIFORM CORTEX"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PIRIFORM CORTEX")
        if (Search_Term(rowSpecimen, ["CHOROID"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12694")) # CHOROID PLEXUS
    elif (Search_Term(rowSpecimen, ["HIPPOCAMPUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12351") # BRAIN, CEREBRUM
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12444")) # HIPPOCAMPUS
    elif (Search_Term(rowSpecimen, ["PIRIFORM CORTEX"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12351") # BRAIN, CEREBRUM
        row['Anatomical_Region_of_Specimen'] = "PIRIFORM CORTEX"
    elif (Search_Term(rowSpecimen, ["CINGULATE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12351") # BRAIN, CEREBRUM
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C52713")) # CINGULATE CORTEX
    elif (Search_Term(rowSpecimen, ["BRAIN STEM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12441") # BRAIN, BRAIN STEM
        if (Search_Term(rowSpecimen, ["CHOROID"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12694")) # CHOROID PLEXUS
        if ((Search_Term(rowSpecimen, ["MEDULLA"])) and not (Search_Term(rowSpecimen, ["EXTRAMEDULLARY"]))):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MEDULLA")
        if (Search_Term(rowSpecimen, ["PYRAMIDAL TRACT"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C142308"))  # PYRAMIDAL TRACTS, BRAINSTEM
            # row['UBERON_Code'] = "UBERON:0022272"  # corticobulbar tract
            if ((Search_Term(rowSpecimen, ["MEDULLA"])) and not (Search_Term(rowSpecimen, ["EXTRAMEDULLARY"]))):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MEDULLA")
    elif (Search_Term(rowSpecimen, ["INTERNAL CAPSULE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12439") # BRAIN
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "INTERNAL CAPSULE")
        if (Search_Term(rowSpecimen, ["PYRAMIDAL TRACT"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C142309"))  # PYRAMIDAL TRACTS, INTERNAL CAPSULE
    elif (Search_Term(rowSpecimen, ["PYRAMIDAL TRACT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12441") # BRAIN, BRAIN STEM
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C142308")) # PYRAMIDAL TRACTS, BRAINSTEM
        #row['UBERON_Code'] = "UBERON:0022272"  # corticobulbar tract
        if ((Search_Term(rowSpecimen, ["MEDULLA"])) and not (Search_Term(rowSpecimen, ["EXTRAMEDULLARY"]))):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MEDULLA")
    elif (Search_Term(rowSpecimen, ["CEREBELLUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12445") # BRAIN, CEREBELLUM
        if (Search_Term(rowSpecimen, ["CHOROID"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12694")) # CHOROID PLEXUS
    elif (Search_Term(rowSpecimen, ["CHOROID"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12694") # BRAIN, CHOROID PLEXUS
    elif (Search_Term(rowSpecimen, ["EPENDYM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12694") # BRAIN, CHOROID PLEXUS
        row['Anatomical_Region_of_Specimen'] = "EPENDYMA"
    elif ((Search_Term(rowSpecimen, ["MEDULLA"])) and not (Search_Term(rowSpecimen, ["EXTRAMEDULLARY"]))):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12442") # BRAIN, MEDULLA OBLONGATA
    elif (Search_Term(rowSpecimen, ["HYPOTHALAMUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12458") # BRAIN, HYPOTHALAMUS
    elif (Search_Term(rowSpecimen, ["THALAMUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12459") # BRAIN, THALAMUS
    elif (Search_Term(rowSpecimen, ["MIDBRAIN"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12510") # BRAIN, MIDBRAIN
    elif (Search_Term(rowSpecimen, ["FOREBRAIN"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C40185") # BRAIN, FOREBRAIN
    elif (Search_Term(rowSpecimen, ["GENICULATE NUCLE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12459") # BRAIN, THALAMUS
        row['Anatomical_Region_of_Specimen'] = "LATERAL GENICULATE NUCLEUS"
    elif (Search_Term(rowSpecimen, ["PITUITARY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12399") # GLAND, PITUITARY
        row['Sample_Type_Group'] = "GLAND" #EXTENSIBLE
        if (Search_Term(rowSpecimen, ["CRANIOPHARYNGEAL"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "RATHKE'S POUCH")
        elif (Search_Term(rowSpecimen, ["RATHKE"])):
            if (Search_Term(rowSpecimen, ["POUCH"])):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "RATHKE'S POUCH")
            elif (Search_Term(rowSpecimen, ["CLEFT"])):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "RATHKE'S CLEFT")
    elif (Search_Term(rowSpecimen, ["MENINGES"])):
        if (Search_Term(rowSpecimen, ["BRAIN"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12439") # BRAIN
            if (Search_Term(rowSpecimen, ["LEPTO"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32979"))  # LEPTOMENINGES
            else:
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12348")) # MENINGES
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12348") # MENINGES
            if (Search_Term(rowSpecimen, ["LEPTO"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32979"))  # LEPTOMENINGES
        if (Search_Term(rowSpecimen, ["LEPTO"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32979"))  # LEPTOMENINGES
    elif (Search_Term(rowSpecimen, ["COCHLEAR"])):
        if (Search_Term(rowSpecimen, ["VESTIBULOCOCHLEAR NUCLEUS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12439") # BRAIN
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12996")) # VESTIBULOCOCHLEAR NERVE
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "VESTIBULOCOCHLEAR NUCLEUS")
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12837") # BRAIN, COCHLEAR NUCLEI
    elif (Search_Term(rowSpecimen, ["VENTRIC"])):
        if (Search_Term(rowSpecimen, ["THIRD"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12827")  # BRAIN, THIRD VENTRICLE
        elif (Search_Term(rowSpecimen, ["FOURTH"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12356")  # BRAIN VENTRICLE
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "FOURTH VENTRICLE")
        elif (Search_Term(rowSpecimen, ["LATERAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12834")  # BRAIN VENTRICLE, LATERAL
        elif (Search_Term(rowSpecimen, ["PERIVENTRICUL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C132390")  # BRAIN, PERIVENTRICULAR REGION
        else:
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12356")) # BRAIN VENTRICLE
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12439") # BRAIN

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["CAUDATE"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12451")) # CAUDATE NUCLEUS
    if (Search_Term(rowSpecimen, ["PUTAMEN"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12452")) # PUTAMEN
    if (Search_Term(rowSpecimen, ["MOTOR CORTEX"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C97339")) # MOTOR CORTEX
    if (Search_Term(rowSpecimen, ["SENSORIMOTOR CORTEX"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C154777")) # SENSORIMOTOR CORTEX
    if (Search_Term(rowSpecimen, ["ACCUMBENS"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C52733")) # NUCLEUS ACCUMBENS
    if (Search_Term(rowSpecimen, ["NUCLEUS OF DIAGONAL BAND"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C97342")) # NUCLEUS OF DIAGONAL BAND
    if (Search_Term(rowSpecimen, ["SUPRATENTORIAL BRAIN"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12512")) # SUPRATENTORIAL BRAIN
    if (Search_Term(rowSpecimen, ["FRONTAL LOBE"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12352")) # FRONTAL LOBE
    if (Search_Term(rowSpecimen, ["INFRATENTORIAL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12509")) # INFRATENTORIAL BRAIN
    if (Search_Term(rowSpecimen, ["MOLECULAR LAYER"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MOLECULAR LAYER")
    if (Search_Term(rowSpecimen, ["ADENOHYPOPHYSIS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12399") # GLAND, PITUITARY
        row['Sample_Type_Group'] = "GLAND" # EXTENSIBLE
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ADENOHYPOPHYSIS")
    if (Search_Term(rowSpecimen, ["GRAY MATTER"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GRAY MATTER")
    if (Search_Term(rowSpecimen, ["WHITE MATTER"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "WHITE MATTER")
    if (Search_Term(rowSpecimen, ["GRANULAR LAYER"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GRANULAR LAYER")
    if (Search_Term(rowSpecimen, ["PITUITARY STALK"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PITUITARY STALK")
    if (Search_Term(rowSpecimen, ["MENINGES", "MENINGITIS", "MENINGEAL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12348")) # MENINGES
    if (Search_Term(rowSpecimen, ["VESTIBULOCOCHLEAR NUCLEUS"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12996")) # VESTIBULOCOCHLEAR NERVE
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "VESTIBULOCOCHLEAR NUCLEUS")
    elif (Search_Term(rowSpecimen, ["COCHLEAR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12837") # BRAIN, COCHLEAR NUCLEI
    if (Search_Term(rowSpecimen, ["NEUROPIL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "NEUROPIL")
    if (Search_Term(rowSpecimen, ["OPTIC TRACT"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "OPTIC TRACT")
    if (Search_Term(rowSpecimen, ["PINEAL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12398")) # PINEAL GLAND
    if (Search_Term(rowSpecimen, ["PARS"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12399")) # GLAND, PITUITARY
        row['Sample_Type_Group'] = "GLAND" # EXTENSIBLE
        if (Search_Term(rowSpecimen, ["NERVOSA"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "NEUROHYPOPHYSIS")
        elif (Search_Term(rowSpecimen, ["INTERMEDIA"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ADENOHYPOPHYSIS")
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PARS INTERMEDIA")
        elif (Search_Term(rowSpecimen, ["DISTALIS"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ADENOHYPOPHYSIS")
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PARS DISTALIS")
    if (Search_Term(rowSpecimen, ["PARENCHYMA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PARENCHYMA")
    if ((Search_Term(row, ["FLUORO"])) and (Search_Term(row, ["JADE"]))):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12439") # BRAIN
        row['Extra_finding_terms'] = Add_Term('Extra_finding_terms', row, "STAIN: FLUORO-JADE B")
        if (row['MISTRESC_raw'] == "DEGENERATION"):
            row['MISTRESC'] = codetermExtract(dfNonNeoplasmCodelist, "C50774") # DEGENERATION
        else:
            row['MISTRESC'] = "SPECIAL STAIN OR METHOD, EXAMINED" #EXTENSIBLE

    ##CELLTYPES
    if (Search_Term(row, ["ACIDOPHIL"])):  ##PITUITARY GLAND
        row['Celltype'] = Add_Term('Celltype', row, "ACIDOPHIL")
    if (Search_Term(row, ["CHROMOPHOBE CELL"])):  ##PITUITARY GLAND
        row['Celltype'] = Add_Term('Celltype', row, "CHROMOPHOBE CELL")
    if (Search_Term(row, ["ASTROCYT"])):
        row['Celltype'] = Add_Term('Celltype', row, "ASTROCYTE")
    if (Search_Term(row, ["NEURON", "NEURODEGEN"])):
        if (Search_Term(row, ["GRANULE CELL"])):
            row['Celltype'] = Add_Term('Celltype', row, "GRANULE CELL")
        else:
            row['Celltype'] = Add_Term('Celltype', row, "NEURON")
    if (Search_Term(row, ["MICROGLIOSIS", "MICROGLIA"])):
        row['Celltype'] = Add_Term('Celltype', row, "MICROGLIA")
    elif (Search_Term(row, ["GLIOSIS", "GLIAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "GLIAL CELL")
    if (Search_Term(row, ["GRANULE CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "GRANULE CELL")
    if (Search_Term(row, ["MONOCYT"])):
        row['Celltype'] = Add_Term('Celltype', row, "MONOCYTE")
    if (Search_Term(row, ["PURKINJE CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "PURKINJE CELL")

    return row

def Heart_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # HEART #UBERON:0000948
    row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12727") # HEART
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12727") # HEART

    if (Search_Term(rowSpecimen, ["VALV"])):
        if (Search_Term(rowSpecimen, ["AORT"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12727")) # VALVE, AORTIC
            if (Search_Term(rowSpecimen, ["ANNULUS"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C130167")) # AORTIC VALVE ANNULUS
            elif (Search_Term(rowSpecimen, ["CUSP"])):
                if (Search_Term(rowSpecimen, ["LEFT"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C127638"))  # AORTIC VALVE, LEFT CORONARY CUSP
                elif (Search_Term(rowSpecimen, ["RIGHT"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C127640"))  # AORTIC VALVE, RIGHT CORONARY CUSP
                elif (Search_Term(rowSpecimen, ["POSTERIOR","NON CORONARY"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C127639"))  # AORTIC VALVE, NON-CORONARY CUSP
        elif (Search_Term(rowSpecimen, ["MITRAL"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12753")) # VALVE, LEFT ATRIOVENTRICULAR
            if (Search_Term(rowSpecimen, ["ANNULUS"])):
                if (Search_Term(rowSpecimen, ["ANTERIOR", "ANTERO"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C127668"))  # MITRAL VALVE, ANTERIOR ANNULUS
                elif (Search_Term(rowSpecimen, ["POSTERIOR", "POSTERO"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C127670"))  # MITRAL VALVE, POSTERIOR ANNULUS
                else:
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C127306")) # MITRAL VALVE ANNULUS
            elif (Search_Term(rowSpecimen, ["CUSP"])):
                if (Search_Term(rowSpecimen, ["ANTERIOR", "ANTERO"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C127669"))  # MITRAL VALVE, ANTERIOR CUSP
                elif (Search_Term(rowSpecimen, ["POSTERIOR", "POSTERO"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C127671"))  # MITRAL VALVE, POSTERIOR CUSP
                else:
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12753"))  # MITRAL VALVE
            else:
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12753")) # MITRAL VALVE
        elif (Search_Term(rowSpecimen, ["TRICUSPID"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12805")) # VALVE, RIGHT ATRIOVENTRICULAR
            if (Search_Term(rowSpecimen, ["ANNULUS"])):
                if (Search_Term(rowSpecimen, ["POSTERIOR", "POSTERO"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C130169"))  # TRICUSPID VALVE, POSTERIOR ANNULUS
                else:
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C130047")) # TRICUSPID VALVE ANNULUS
            elif (Search_Term(rowSpecimen, ["CUSP"])):
                if (Search_Term(rowSpecimen, ["ANTERIOR", "ANTERO"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C32799"))  # TRICUSPID VALVE, ANTERIOR CUSP
                elif (Search_Term(rowSpecimen, ["POSTERIOR", "POSTERO"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C33055"))  # TRICUSPID VALVE, POSTERIOR CUSP
                elif (Search_Term(rowSpecimen, ["POSTERIOR", "POSTERO"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C33534"))  # TRICUSPID VALVE, SEPTAL CUSP
                else:
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12805")) # TRICUSPID VALVE
            else:
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12805")) # TRICUSPID VALVE
        elif (Search_Term(rowSpecimen, ["ATRIOVENTRICULAR"])):
            if (Search_Term(rowSpecimen, ["LEFT"])):
                row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12753"))  # VALVE, LEFT ATRIOVENTRICULAR
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12753"))  # MITRAL VALVE
            elif (Search_Term(rowSpecimen, ["RIGHT"])):
                row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12805"))  # VALVE, RIGHT ATRIOVENTRICULAR
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12805"))  # TRICUSPID VALVE
            else:
                row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12753"))  # VALVE, LEFT ATRIOVENTRICULAR
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12753"))  # MITRAL VALVE
                row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12805"))  # VALVE, RIGHT ATRIOVENTRICULAR
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12805"))  # TRICUSPID VALVE
                row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MISPEC")
                row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MIANTREG")
        elif (Search_Term(rowSpecimen, ["LUNG", "PULMON"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12775"))  # VALVE, PULMONARY
            if (Search_Term(rowSpecimen, ["CUSP"])):
                if (Search_Term(rowSpecimen, ["ANTERIOR", "ANTERO"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C127673"))  # PULMONARY VALVE, ANTERIOR CUSP
                elif (Search_Term(rowSpecimen, ["LEFT"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C127674"))  # PULMONARY VALVE, LEFT CUSP
                elif (Search_Term(rowSpecimen, ["RIGHT"])):
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C127675"))  # PULMONARY VALVE, RIGHT CUSP
                else:
                    row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12775")) # PULMONARY VALVE
        else:
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12729"))  # VALVE, CARDIAC
    elif (Search_Term(rowSpecimen, ["PERICARD"])):
        if (Search_Term(rowSpecimen, ["CAVITY"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C38662")  # BODY CAVITY, PERICARDIAL
        else:
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C13005"))  # PERICARDIUM

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["SMOOTH MUSCLE ACTIN"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12437"))  # MUSCLE, SMOOTH
        row['Extra_finding_terms'] = "STAIN: SMOOTH MUSCLE ACTIN"
        if (Search_Term(rowSpecimen, ["POSITIVE"])):
            row['MISTRESC'] = "SPECIAL STAIN OR METHOD, POSITIVE" # EXTENSIBLE
        elif (Search_Term(rowSpecimen, ["NEGAITIVE","NULL","NO VISIBLE LESIONS"])):
            row['MISTRESC'] = "SPECIAL STAIN OR METHOD, NEGATIVE" # EXTENSIBLE
        else:
            row['MISTRESC'] = "SPECIAL STAIN OR METHOD, EXAMINED" # EXTENSIBLE
    if (Search_Term(rowSpecimen, ["ATRIUM", "ATRIAL"])):
        if ("RIGHT" in row.at['MISPEC']):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12868"))  # HEART, RIGHT ATRIUM
        elif ("LEFT" in row.at['Non Grade Modifiers']):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12869"))  # HEART, LEFT ATRIUM
        else:
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12728"))  # HEART, ATRIUM
    if (Search_Term(rowSpecimen, ["SEPTUM"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C49485"))  # HEART, SEPTUM
    if ((Search_Term(rowSpecimen, ["VENTRIC"])) and not (Search_Term(rowSpecimen, ["ATRIOVENTRICULAR"]))):
        if (Search_Term(rowSpecimen, ["INTERVENTRICULAR", "INTRAVENTRICULAR"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C32874"))  # INTERVENTRICULAR SEPTUM
        else:
            if (Search_Term(rowSpecimen, ["LEFT"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12871"))  # HEART, LEFT VENTRICLE
            elif (Search_Term(rowSpecimen, ["RIGHT"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12870"))  # HEART, RIGHT VENTRICLE
            else:
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12730"))  # HEART, VENTRICLE
    if (Search_Term(rowSpecimen, ["ATRIOVENTRICULAR"])):
        if (Search_Term(rowSpecimen, ["LEFT"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12753"))  # VALVE, LEFT ATRIOVENTRICULAR
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12753"))  # MITRAL VALVE
        elif (Search_Term(rowSpecimen, ["RIGHT"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12805"))  # VALVE, RIGHT ATRIOVENTRICULAR
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C12805"))  # TRICUSPID VALVE
        else:
            row['MISPEC'] = Add_Term('MISPEC', row, "VALVE, LEFT ATRIOVENTRICULAR")
            row['MIANTREG'] = Add_Term('MIANTREG', row, "MITRAL VALVE")
            row['MISPEC'] = Add_Term('MISPEC', row, "VALVE, RIGHT ATRIOVENTRICULAR")
            row['MIANTREG'] = Add_Term('MIANTREG', row, "TRICUSPID VALVE")
            row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MISPEC")
            row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MIANTREG")
            #row['UBERON_Code'] = "UBERON:0002133"  # atrioventricular valve
    if (Search_Term(rowSpecimen, ["BASE", "BASILAR"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C48589"))  # HEART, BASE
        #row['UBERON_Code'] = "UBERON:0035213"  # basal zone of heart
    if (Search_Term(rowSpecimen, ["APEX", "APICAL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32126"))  # HEART, APEX
    if (Search_Term(rowSpecimen, ["EPICARD"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C13164"))  # EPICARDIUM
    if (Search_Term(rowSpecimen, ["ENDOCARD"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C13004"))  # ENDOCARDIUM
    if (Search_Term(rowSpecimen, ["PERICARD", "SEROSA"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C13005"))  # PERICARDIUM
    if (Search_Term(rowSpecimen, ["MYOCARD", "CARDIOMYOP"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MYOCARDIUM")
    if (Search_Term(rowSpecimen, ["MYOFIBER"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MYOCARDIUM")
        row['Celltype'] = Add_Term('Celltype', row, "CARDIOMYOCYTE")
    if (Search_Term(rowSpecimen, ["CHORDA", "TENDIN"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CHORDAE TENDINEAE")
    if ((Search_Term(rowSpecimen, ["AORT"])) and not (Search_Term(rowSpecimen, ["VALV"]))):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12669") # ARTERY, AORTA
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12679") #  VESSEL, BLOOD
        if (Search_Term(rowSpecimen, ["ROOT"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C127641"))  # ASCENDING AORTA, AORTIC ROOT
        elif (Search_Term(rowSpecimen, ["AORTIC ARCH"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12669") # ARTERY, AORTA
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfSpecimenCodelist, "C32123"))  # AORTIC ARCH
        else:
            if (Search_Term(rowSpecimen, ["ASCENDING"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12669") # ARTERY, AORTA
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32150"))  # ASCENDING AORTA
            else:
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12669") # ARTERY, AORTA
    elif (Search_Term(rowSpecimen, ["ARTER"])):
        if (Search_Term(rowSpecimen, ["LUNG", "PULMON"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12774"))  # PULMONARY ARTERY BRANCH
        else:
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12843"))  # ARTERY, CORONARY
        if (Search_Term(rowSpecimen, ["ARTERIOL"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ARTERIOLE")
    if ((Search_Term(rowSpecimen, ["PAPILLARY"])) and (Search_Term(rowSpecimen, ["MUSCLE"]))):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PAPILLARY MUSCLE")
    elif (Search_Term(rowSpecimen, ["PAPILLARY"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PAPILLARY")
    if (Search_Term(rowSpecimen, ["FAT", "ADIPO"])):
        if (Search_Term(rowSpecimen, ["BROWN"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C32235"))  # ADIPOSE TISSUE, BROWN
        elif (Search_Term(rowSpecimen, ["WHITE"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C33889"))  #  ADIPOSE TISSUE, WHITE
        else:
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12472"))  # ADIPOSE TISSUE

    ##CELLTYPES
    if (Search_Term(rowSpecimen, ["CARDIOMYOCYTE"])):
        row['Celltype'] = Add_Term('Celltype', row, "CARDIOMYOCYTE")

    return row

def Vascular_System_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # VASCULAR SYSTEM #UBERON:0007798
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12679") # VESSEL, BLOOD

    if (Search_Term(rowSpecimen, ["AORT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12669")  # ARTERY, AORTA
        if (Search_Term(rowSpecimen, ["ASCENDING"])):
            if (Search_Term(rowSpecimen, ["ROOT"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C127641"))  # ASCENDING AORTA, AORTIC ROOT
            elif (Search_Term(rowSpecimen, ["SINOTUBUL"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C127642"))  # ASCENDING AORTA, SINOTUBULAR JUNCTION
            elif (Search_Term(rowSpecimen, ["VALSALVA"])):
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33557"))  # ASCENDING AORTA, SINUS OF VALSALVA
            else:
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32150"))  # ASCENDING AORTA

    elif (Search_Term(rowSpecimen, ["VEIN"])):
        if (Search_Term(rowSpecimen, ["PAMPINIFORM PLEXUS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12814") # VEIN
            row['Anatomical_Region_of_Specimen'] = "PAMPINIFORM PLEXUS"
        elif (Search_Term(rowSpecimen, ["AURICULAR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77673") # VEIN, AURICULAR
        elif (Search_Term(rowSpecimen, ["BRACHIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12883") # VEIN, BRACHIAL
        elif (Search_Term(rowSpecimen, ["TAIL", "CAUDAL"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C92598")) # VEIN, CAUDAL
            row['Sample_Type_Group'] = "INJECTION/APPLICATION SITE"
        elif (Search_Term(rowSpecimen, ["CEPHALIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32286") # VEIN, CEPHALIC
        elif (Search_Term(rowSpecimen, ["FEMORAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12716") # VEIN, FEMORAL
        elif (Search_Term(rowSpecimen, ["JUGULAR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12738") # VEIN, JUGULAR
        elif (Search_Term(rowSpecimen, ["MESENTERIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53055") # VEIN, MESENTERIC
        elif (Search_Term(rowSpecimen, ["PORTAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33343") # VEIN, PORTAL
        elif (Search_Term(rowSpecimen, ["PULMONARY"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12776") # VEIN, PULMONARY
        elif (Search_Term(rowSpecimen, ["RENAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33462") # VEIN, RENAL
        elif (Search_Term(rowSpecimen, ["SAPHENA","SAPHENOUS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33511") # VEIN, SAPHENA
        elif (Search_Term(rowSpecimen, ["VENA CAVA"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12817") # VEIN, VENA CAVA
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12814") # VEIN

    elif (Search_Term(rowSpecimen, ["ARTERY"])):
        if (Search_Term(rowSpecimen, ["AURICULAR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52849") # ARTERY, AURICULAR
        elif (Search_Term(rowSpecimen, ["BRACHIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12681") # ARTERY, BRACHIAL
        elif (Search_Term(rowSpecimen, ["CAROTID"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12687") # ARTERY, CAROTID
        elif (Search_Term(rowSpecimen, ["CORONARY"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12843") # ARTERY, CORONARY
        elif (Search_Term(rowSpecimen, ["FEMORAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12715") # ARTERY, FEMORAL
        elif (Search_Term(rowSpecimen, ["ILIAC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12733") # ARTERY, ILIAC
        elif (Search_Term(rowSpecimen, ["INTERNAL THORACIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52941") # ARTERY, INTERNAL THORACIC
        elif (Search_Term(rowSpecimen, ["MESENTERIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52975") # ARTERY, MESENTERIC
        elif (Search_Term(rowSpecimen, ["PULMONARY"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12774") # ARTERY, PULMONARY
        elif (Search_Term(rowSpecimen, ["RENAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12778") # ARTERY, RENAL
        elif (Search_Term(rowSpecimen, ["SPINAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33587") # ARTERY, SPINAL
        elif (Search_Term(rowSpecimen, ["SUBCLAVIAN"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33643") # ARTERY, SUBCLAVIAN
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12372") # ARTERY


        if (Search_Term(rowSpecimen, ["CORONARY"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12727") # HEART
            row['MISPEC'] = Add_Term('MISPEC', row, "ARTERY, CORONARY")
        elif (Search_Term(rowSpecimen, ["FEMORAL"])):
            row['MISPEC'] = Add_Term('MISPEC', row, "ARTERY, FEMORAL")
        elif (Search_Term(rowSpecimen, ["MESENTER"])):
            row['MISPEC'] = Add_Term('MISPEC', row, "ARTERY, MESENTERIC")
        elif (Search_Term(rowSpecimen, ["RENAL"])):
            row['MISPEC'] = Add_Term('MISPEC', row, "ARTERY, RENAL")
            if (Search_Term(row, ["ARTERIOLE"])):
                row['Anatomical_Region_of_Specimen'] = "ARTERIOLE"
        elif (Search_Term(rowSpecimen, ["CAROTID"])):
            row['MISPEC'] = Add_Term('MISPEC', row, "ARTERY, CAROTID")
        elif (Search_Term(rowSpecimen, ["BRACHIAL"])):
            row['MISPEC'] = Add_Term('MISPEC', row, "ARTERY, BRACHIAL")
        elif (Search_Term(rowSpecimen, ["PULMONARY"])):
            row['MISPEC'] = Add_Term('MISPEC', row, "ARTERY, PULMONARY")
        else:
            row['MISPEC'] = Add_Term('MISPEC', row, "ARTERY")

    elif (Search_Term(rowSpecimen, ["TUNICA MEDIA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12679") # VESSEL, BLOOD
        row['Anatomical_Region_of_Specimen'] = "TUNICA MEDIA"
    elif (Search_Term(rowSpecimen, ["MESENTER"])):
        if (Search_Term(rowSpecimen, ["LYMPH NODE"])):
            row = Lymphoid_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist)
        elif (Search_Term(rowSpecimen, ["ARTERIOLE"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C52975")) # ARTERY, MESENTERIC
            row['Anatomical_Region_of_Specimen'] = "ARTERIOLE"
        elif (Search_Term(rowSpecimen, ["ARTER"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C52975")) # ARTERY, MESENTERIC
        elif (Search_Term(rowSpecimen, ["VEIN"])):
            row['MISPEC'] = row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53055") # VEIN, MESENTERIC
        elif (Search_Term(rowSpecimen, ["FAT", "ADIPO"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C33103")) # MESENTERY
            if (Search_Term(rowSpecimen, ["BROWN"])):
                row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C32235"))  # ADIPOSE TISSUE, BROWN
            elif (Search_Term(rowSpecimen, ["WHITE"])):
                row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C33889"))  #  ADIPOSE TISSUE, WHITE
            else:
                row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12472"))  # ADIPOSE TISSUE
        elif (Search_Term(rowSpecimen, ["VASCULAR", "VESSEL"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C33103")) # MESENTERY
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12679")) # VESSEL, BLOOD
        else:
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C33103")) # MESENTERY
            row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C25444") # BODY CAVITY
    else:
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12679")) # VESSEL, BLOOD

    return row

def Stomach_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12391") # STOMACH

    if (Search_Term(rowSpecimen, ["BRUNNER"])):
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12386")  # SMALL INTESTINE
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12263") # SMALL INTESTINE, DUODENUM
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C13010")) # GLAND, BRUNNER'S
    elif (Search_Term(rowSpecimen, ["PYLORUS", "PYLORIC"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12260") # STOMACH, PYLORUS
    elif (Search_Term(rowSpecimen, ["CARDIA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12256") # STOMACH, CARDIA
    elif (Search_Term(rowSpecimen, ["FUND"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12257") # STOMACH, FUNDUS
    elif (Search_Term(rowSpecimen, ["NONGLANDULAR","NON GLANDULAR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77662") #  STOMACH, NONGLANDULAR
    elif (Search_Term(rowSpecimen, ["GLANDULAR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77661") #  STOMACH, GLANDULAR
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12391") # STOMACH

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["ANTRUM"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12259"))  # ANTRUM PYLORI
    if (Search_Term(rowSpecimen, ["GASTROESOPHAGEAL JUNCTION"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32668"))  # GASTROESOPHAGEAL JUNCTION
    if (Search_Term(rowSpecimen, ["GASTRIC PIT"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GASTRIC PIT")
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GASTRIC GLAND")
    if (Search_Term(rowSpecimen, ["GASTRIC GLAND"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GASTRIC GLAND")
    if (Search_Term(rowSpecimen, ["PYLORIC GLAND", "PYLORUS GLAND"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PYLORIC GLAND")
    if (Search_Term(rowSpecimen, ["LIMITING RIDGE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "LIMITING RIDGE")
    if (Search_Term(rowSpecimen, ["BASAL LAYER"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "BASAL LAYER")
    if (Search_Term(rowSpecimen, ["MYENTERIC PLEXUS"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MYENTERIC PLEXUS")
    if (Search_Term(rowSpecimen, ["CRYPT"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CRYPT")

    ##CELLTYPES
    if (Search_Term(rowSpecimen, ["CHIEF CELL"])):
        row['Celltype'] = "CHIEF CELL"
    if (Search_Term(rowSpecimen, ["GOBLET CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "GOBLET CELL")
    if (Search_Term(rowSpecimen, ["MUCOUS NECK CELL"])):
        row['Celltype'] = "MUCOUS NECK CELL"
    if (Search_Term(rowSpecimen, ["PARIETAL CELL"])):
        row['Celltype'] = "PARIETAL CELL"

    return row

def Intestine_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    if (Search_Term(rowSpecimen, ["RECTUM"])):
        if (Search_Term(rowSpecimen, ["ANUS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92217")  # LARGE INTESTINE, RECTUM/LARGE INTESTINE, ANUS
            row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12379") # LARGE INTESTINE
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12390")  # LARGE INTESTINE, RECTUM
            row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12379") # LARGE INTESTINE
        if (Search_Term(rowSpecimen, ["GLAND"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77620"))  # GLAND, PERIANAL
    elif (Search_Term(rowSpecimen, ["ANUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C43362") # LARGE INTESTINE, ANUS
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12379") # LARGE INTESTINE
        if (Search_Term(rowSpecimen, ["GLAND"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77620"))  # GLAND, PERIANAL
    elif (Search_Term(rowSpecimen, ["PEYER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12387") # SMALL INTESTINE, ILEUM
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12771")) # PEYER'S PATCH
        row['Sample_Type_Group'] = "LYMPHOID TISSUE" # EXTENSIBLE
    elif (Search_Term(rowSpecimen, ["TRACT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12736") # INTESTINE
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C34082"))  # GASTROINTESTINAL TRACT
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12736") # INTESTINE
    elif (Search_Term(rowSpecimen, ["DUODENUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12263") # SMALL INTESTINE, DUODENUM
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12386") # SMALL INTESTINE
    elif (Search_Term(rowSpecimen, ["ILEUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12387") # SMALL INTESTINE, ILEUM
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12386") # SMALL INTESTINE
    elif (Search_Term(rowSpecimen, ["JEJUNUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12388") # "SMALL INTESTINE, JEJUNUM
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12386") # SMALL INTESTINE
    elif (Search_Term(rowSpecimen, ["SACCULUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C88024") # SMALL INTESTINE, SACCULUS ROTUNDUS
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12386") # SMALL INTESTINE
    elif (Search_Term(rowSpecimen, ["SMALL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12386") # SMALL INTESTINE
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12386") # SMALL INTESTINE
    elif (Search_Term(rowSpecimen, ["APPENDIX"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12380") # LARGE INTESTINE, APPENDIX
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12379") # LARGE INTESTINE
    elif (Search_Term(rowSpecimen, ["CECUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12381") # LARGE INTESTINE, CECUM
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12379") # LARGE INTESTINE
    elif (Search_Term(rowSpecimen, ["COLON"])):
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12379") # LARGE INTESTINE
        if (Search_Term(rowSpecimen, ["COLONY", "COLONIES"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12736") # INTESTINE
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12382") # LARGE INTESTINE, COLON
    elif (Search_Term(rowSpecimen, ["LARGE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12379") # LARGE INTESTINE
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12379") # LARGE INTESTINE
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12736") # INTESTINE
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12736") # INTESTINE

        ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["LYMPHOID TISSU", "GALT"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12936"))  # GUT-ASSOCIATED LYMPHOID TISSUE
        row['Sample_Type_Group'] = "LYMPHOID TISSUE" # EXTENSIBLE
    if (Search_Term(rowSpecimen, ["VILL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "VILLI")
    if (Search_Term(rowSpecimen, ["CRYPT"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CRYPT")
    if (Search_Term(rowSpecimen, ["MYENTERIC PLEX"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MYENTERIC PLEXUS")
        row['Celltype'] = Add_Term('Celltype', row, "NEURON")
    if (Search_Term(rowSpecimen, ["LACTEAL", "GALACTOCELE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "LACTEAL")
    if (Search_Term(rowSpecimen, ["BRUNNER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12263") # SMALL INTESTINE, DUODENUM
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C13010")) # "GLAND, BRUNNER'S
    if (Search_Term(rowSpecimen, ["CENTRAL VEIN"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12814"))  # VEIN
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CENTRAL VEIN")

    ##CELLTYPES
    if (Search_Term(rowSpecimen, ["ENTEROCYTE"])):
        row['Celltype'] = Add_Term('Celltype', row, "ENTEROCYTE")
    if (Search_Term(rowSpecimen, ["GOBLET CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "GOBLET CELL")
    if (Search_Term(rowSpecimen, ["PANETH CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "PANETH CELL")

    return row

def Kidney_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    if (Search_Term(rowSpecimen, ["URETER"])):
        row['Sample_Type_Group'] = codetermExtract(dfAnatLocCodelist, "C12413")  # URINARY SYSTEM
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12416")  # URETER
    elif (Search_Term(rowSpecimen, ["URETHRA"])):
        row['Sample_Type_Group'] = codetermExtract(dfAnatLocCodelist, "C12413")  # URINARY SYSTEM
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12417")  # URETHRA
    elif (Search_Term(rowSpecimen, ["URINARY"])):
        row['Sample_Type_Group'] = codetermExtract(dfAnatLocCodelist, "C12413")  # URINARY SYSTEM
        if (Search_Term(rowSpecimen, ["TRACT"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12414")  # URINARY BLADDER
            row['Anatomical_Region_of_Specimen'] = "URINARY TRACT"
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12414")  # URINARY BLADDER
    elif (Search_Term(rowSpecimen, ["KIDNEY"])):
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12415")  # KIDNEY
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12415")  # KIDNEY

        ###EXTRA TERMS
        if (Search_Term(rowSpecimen, ["PELVI"])):
            row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12887")) # RENAL PELVIS
        if (Search_Term(rowSpecimen, ["CORTEX", "CORTICAL"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12739"))  # KIDNEY, CORTEX
        if (Search_Term(rowSpecimen, ["GLOMERUL"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C13250"))  # GLOMERULUS
        if (Search_Term(rowSpecimen, ["MEDULLA"])):
            if (Search_Term(rowSpecimen, ["CORTICO"])):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CORTICOMEDULLARY JUNCTION")
            else:
                row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12740"))  # KIDNEY, MEDULLA
        if (Search_Term(rowSpecimen, ["HILUS", "HILUM"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32740"))  # KIDNEY, HILUM
        if (Search_Term(rowSpecimen, ["STRIPE"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12740"))  # KIDNEY, MEDULLA
            if ((Search_Term(rowSpecimen, ["INNER STRIPE"])) and (Search_Term(rowSpecimen, ["OUTER STRIPE"]))):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      "INNER STRIPE")
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      "OUTER STRIPE")
            elif (Search_Term(rowSpecimen, ["INNER STRIPE"])):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      "INNER STRIPE")
            elif (Search_Term(rowSpecimen, ["OUTER STRIPE"])):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      "OUTER STRIPE")
        if (Search_Term(rowSpecimen, ["HENLE"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "LOOP OF HENLE")
        if (Search_Term(rowSpecimen, ["TUBUL"])):
            if (Search_Term(rowSpecimen, ["DISTAL"])):
                if (Search_Term(rowSpecimen, ["CONVOLUTED"])):
                    row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                          "DISTAL CONVOLUTED TUBULE")
                else:
                    row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                          "DISTAL TUBULE")
            elif (Search_Term(rowSpecimen, ["PROXIMAL"])):
                if (Search_Term(rowSpecimen, ["CONVOLUTED"])):
                    row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                          "PROXIMAL CONVOLUTED TUBULE")
                else:
                    row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                          "PROXIMAL TUBULE")
            else:
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      "TUBULE")
        if (Search_Term(rowSpecimen, ["THIN SEGMENT"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "LOOP OF HENLE")
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "THIN SEGMENT")
        if (Search_Term(rowSpecimen, ["SUBCAPSUL"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "SUBCAPSULE")
        elif (Search_Term(rowSpecimen, ["EXTRACAPSUL"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "EXTRACAPSULAR")
        elif (Search_Term(rowSpecimen, ["CAPSUL"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CAPSULE")
        if (Search_Term(rowSpecimen, ["MESANGIUM"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "MESANGIUM")
        if ((Search_Term(rowSpecimen, ["PAPILLA"])) and not (Search_Term(rowSpecimen, ["PAPILLARY DUCT"]))):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12740"))  # KIDNEY, MEDULLA
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PAPILLA")
        if (Search_Term(rowSpecimen, ["ARTER"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12778"))  # RENAL ARTERY
            if (Search_Term(rowSpecimen, ["ARTERIOLE"])):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      "ARTERIOLE")
    if (Search_Term(rowSpecimen, ["UROTHELI", "URINARY EPITHELIUM"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "UROTHELIUM")
    if (Search_Term(rowSpecimen, ["COLLECTING DUCT", "PAPILLARY DUCT"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "COLLECTING DUCT")

    ##CELLTYPES
    if (Search_Term(rowSpecimen, ["TUBULE CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "RENAL TUBULE CELL")
    if (Search_Term(rowSpecimen, ["JUXTAGLOMERULAR CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "JUXTAGLOMERULAR CELL")
    if (Search_Term(rowSpecimen, ["TRANSITIONAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "TRANSITIONAL CELL")

    return row

def Hepatobiliary_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # HEPATOBILIARY SYSTEM #UBERON:0002423
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12392")  # LIVER

    if (Search_Term(rowSpecimen, ["BILE DUCT", "BILIARY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12392")  # LIVER
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12376"))  # BILE DUCT
    elif (Search_Term(rowSpecimen, ["BETA CATENIN STAIN"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12392")  # LIVER
        row['Extra_finding_terms'] = "STAIN: BETA CATENIN"
    elif (Search_Term(rowSpecimen, ["GALLBLADDER"])):
        if (Search_Term(rowSpecimen, ["LIVER"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77669")  # LIVER/GALLBLADDER
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12377")  # GALLBLADDER
    elif (Search_Term(rowSpecimen, ["LIVER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12392")  # LIVER
    elif (Search_Term(rowSpecimen, ["BILE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C13192")  # BILE
    else:
        row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MISPEC")

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["PARENCHYMA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PARENCHYMA")
    if (Search_Term(rowSpecimen, ["LOBE"])):
        if (Search_Term(rowSpecimen, ["CAUDAL"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33000"))  # LIVER, CAUDATE LOBE
        elif (Search_Term(rowSpecimen, ["QUADRATE"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C112404"))  # LIVER, QUADRATE LOBE
        elif (Search_Term(rowSpecimen, ["LEFT"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32965"))  # LIVER, LEFT LOBE
        elif (Search_Term(rowSpecimen, ["RIGHT"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33481"))  # LIVER, RIGHT LOBE
    if (Search_Term(rowSpecimen, ["CANALIC"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "BILE CANALICULUS")
    if (Search_Term(rowSpecimen, ["PORTAL"])):
        if (Search_Term(rowSpecimen, ["PERIPORTAL"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "PERIPORTAL")
        elif (Search_Term(rowSpecimen, ["PORTAL TRIAD", "PORTAL TRACT"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "PORTAL TRIAD")
        elif (Search_Term(rowSpecimen, ["PORTAL VEIN", "PORTAL VESSEL"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33343")) # PORTAL VEIN
        else:
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PORTAL")
    if (Search_Term(rowSpecimen, ["CENTRILOBULAR"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "CENTRILOBULAR")
    if (Search_Term(rowSpecimen, ["MIDZONAL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MIDZONAL")
    if (Search_Term(rowSpecimen, ["CENTRAL VEIN"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12814")) # VEIN
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "CENTRAL VEIN")
    if (Search_Term(rowSpecimen, ["LYMPH NODE"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77640"))  # LYMPH NODE, HEPATIC
    if (Search_Term(rowSpecimen, ["SUBCAPSULAR SINUS"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77640"))  # LYMPH NODE, HEPATIC
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "SUBCAPSULAR SINUS")

    ##CELLTYPES
    if (Search_Term(rowSpecimen, ["HEPATOCELLULAR", "HEPATOCYTE"])):
        row['Celltype'] = Add_Term('Celltype', row, "HEPATOCYTE")
    if ((Search_Term(rowSpecimen, ["CLEAR CELL"])) and not (Search_Term(rowSpecimen, ["MONONUCLEAR CELL"]))):
        row['Celltype'] = Add_Term('Celltype', row, "CLEAR CELL")
    if (Search_Term(rowSpecimen, ["BILIARY CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "BILIARY CELL")
    if (Search_Term(rowSpecimen, ["PERISINUSOIDAL CELL", "ITO CELL", "STELLATE CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "ITO CELL")
    elif (Search_Term(rowSpecimen, ["SINUSOIDAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "SINUSOIDAL CELL")
    if (Search_Term(rowSpecimen, ["KUPFFER"])):
        row['Celltype'] = Add_Term('Celltype', row, "KUPFFER CELL")
    if (Search_Term(rowSpecimen, ["OVAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "OVAL CELL")

    return row

def Lung_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C92218")  # LUNG/BRONCHUS

    if (Search_Term(rowSpecimen, ["LUNG"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92218")  # LUNG/BRONCHUS
    elif (Search_Term(rowSpecimen, ["PLEURA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12469")  # PLEURA
    else:
        row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MISPEC")

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["MYOCARD", "CARDIOMYOP"])):
        row = Heart_SPECIMEN_Mapping(row, rowSpecimen)
    if (Search_Term(rowSpecimen, ["ARTER"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12774")) # ARTERY, PULMONARY
    if (Search_Term(rowSpecimen, ["ALVEOL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12986"))  # ALVEOLUS
    if (Search_Term(rowSpecimen, ["BRONCHIO", "BRONCHO"])):
        if (Search_Term(rowSpecimen, ["TERMINAL"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "TERMINAL BRONCHIOLE")
        else:
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "BRONCHIOLE")
    elif (Search_Term(rowSpecimen, ["BRONCHUS"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12683"))  # BRONCHUS
    elif (Search_Term(rowSpecimen, ["BRONCHI"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "BRONCHI")
    if ((Search_Term(rowSpecimen, ["PLEURA"])) and not ("PLEURA" in row.at['Sample Type'])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12469"))  # PLEURA
    elif (Search_Term(rowSpecimen, ["SUBSEROSA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "SUBSEROSA")
    elif ((Search_Term(rowSpecimen, ["SEROSA"])) and not (Search_Term(rowSpecimen, ["PLEURA"]))):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12469"))  # PLEURA
    if (Search_Term(rowSpecimen, ["BALT","BRONCHUS ASSOCIATED LYMPHOID"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C32234")) # BRONCHUS-ASSOCIATED LYMPHOID TISSUE
    if (Search_Term(rowSpecimen, ["PERICARD"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C13005"))  # PERICARDIUM

    ##CELLTYPES
    if ((Search_Term(rowSpecimen, ["C CELL", "PARAFOLLICULAR CELL"])) and not (
    Search_Term(rowSpecimen, ["IC CELL"]))):
        row['Celltype'] = "PARAFOLLICULAR CELL"
    if (Search_Term(rowSpecimen, ["CLARA CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "CLARA CELL")
    if (Search_Term(rowSpecimen, ["KUPFFER"])):
        row['Celltype'] = Add_Term('Celltype', row, "KUPFFER CELL")
    if (Search_Term(rowSpecimen, ["PNEUMOCYTE TYPE II", "TYPE II PNEUMOCYTE", "PNEUMOCYTE TYPE 2", "TYPE 2 PNEUMOCYTE"])):
        row['Celltype'] = Add_Term('Celltype', row, "TYPE II PNEUMOCYTE")
    elif (Search_Term(rowSpecimen, ["PNEUMOCYTE TYPE I", "TYPE I PNEUMOCYTE", "PNEUMOCYTE TYPE 1", "TYPE 1 PNEUMOCYTE"])):
        row['Celltype'] = Add_Term('Celltype', row, "TYPE I PNEUMOCYTE")
    elif (Search_Term(rowSpecimen, ["PNEUMOCYTE"])):
        row['Celltype'] = Add_Term('Celltype', row, "PNEUMOCYTE")

    return row

def Lymphoid_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = "LYMPHOID SYSTEM" # EXTENSIBLE

    if (Search_Term(rowSpecimen, ["LYMPH NODE"])):
        if (Search_Term(rowSpecimen, ["AXILLA"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33073")  # LYMPH NODE, AXILLARY
        elif (Search_Term(rowSpecimen, ["BRACHIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92221")  # LYMPH NODE, BRACHIAL
        elif (Search_Term(rowSpecimen, ["BRONCHIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32232")  # LYMPH NODE, BRONCHIAL
        elif (Search_Term(rowSpecimen, ["CERVICAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32298")  # LYMPH NODE, CERVICAL
        elif (Search_Term(rowSpecimen, ["HEPATIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77640")  # LYMPH NODE, HEPATIC
        elif (Search_Term(rowSpecimen, ["ILEOCECOCOLIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77653")  # LYMPH NODE, ILEOCECOCOLIC
        elif (Search_Term(rowSpecimen, ["ILIAC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32761")  # LYMPH NODE, ILIAC
        elif (Search_Term(rowSpecimen, ["INGUINAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32801")  # LYMPH NODE, INGUINAL
        elif (Search_Term(rowSpecimen, ["INTERCOSTAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77652")  # LYMPH NODE, INTERCOSTAL
        elif (Search_Term(rowSpecimen, ["LUMBAR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77643")  # LYMPH NODE, LUMBAR
        elif (Search_Term(rowSpecimen, ["MAMMARY GLAND"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32853")  # LYMPH NODE, MAMMARY GLAND
        elif (Search_Term(rowSpecimen, ["MANDIBULAR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77650")  # LYMPH NODE, MANDIBULAR
        elif (Search_Term(rowSpecimen, ["MEDIASTINAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33073")  # LYMPH NODE, MEDIASTINAL
        elif (Search_Term(rowSpecimen, ["MESENTER"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77641")  # LYMPH NODE, MESENTERIC
        elif (Search_Term(rowSpecimen, ["PANCREATIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77642")  # LYMPH NODE, PANCREATIC
        elif (Search_Term(rowSpecimen, ["PAROTID"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33278")  # LYMPH NODE, PAROTID
        elif (Search_Term(rowSpecimen, ["POPLITEAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53146")  # LYMPH NODE, POPLITEAL
        elif (Search_Term(rowSpecimen, ["PORTAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77645")  # LYMPH NODE, PORTAL
        elif (Search_Term(rowSpecimen, ["REGIONAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C49018")  # LYMPH NODE, REGIONAL
        elif (Search_Term(rowSpecimen, ["RENAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77646")  # LYMPH NODE, RENAL
        elif (Search_Term(rowSpecimen, ["PHARYNGEAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77649")  # LYMPH NODE, RETROPHARYNGEAL
        elif (Search_Term(rowSpecimen, ["SACRAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77647")  # LYMPH NODE, SACRAL
        elif (Search_Term(rowSpecimen, ["SUBILIAC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92594")  # LYMPH NODE, SUBILIAC
        elif (Search_Term(rowSpecimen, ["SUBLINGUAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92434")  # LYMPH NODE, SUBLINGUAL
        elif (Search_Term(rowSpecimen, ["TRACHEOBRONCHIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77651")  # LYMPH NODE, TRACHEOBRONCHIAL
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12745")  # LYMPH NODE
    elif (Search_Term(rowSpecimen, ["GUT ASSOCIATED", "GALT"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12936"))  # GUT-ASSOCIATED LYMPHOID TISSUE
        row['Sample_Type_Group'] = "LYMPHOID TISSUE" # EXTENSIBLE
    else:
        if (Search_Term(rowSpecimen, ["HEMOLYMPHORETICULAR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C41168")  # HEMOLYMPHORETICULAR TISSUE
        elif (Search_Term(rowSpecimen, ["BRONCH"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32234")  # BRONCHUS-ASSOCIATED LYMPHOID
        elif (Search_Term(rowSpecimen, ["NASAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77659")  # NASAL-ASSOCIATED LYMPHOID TISSUE
        else:
            row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MISPEC")

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["VASCULA", "VESSEL"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C33038")) # VESSEL, LYMPHATIC
    if ((Search_Term(rowSpecimen, ["MEDULLARY"])) and not (Search_Term(rowSpecimen, ["EXTRAMEDULLARY"]))):
        if (Search_Term(rowSpecimen, ["MEDULLARY CORD"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "MEDULLARY CORD")
        elif (Search_Term(rowSpecimen, ["SINUS"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "MEDULLARY SINUS")
        else:
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "MEDULLARY")
    elif (Search_Term(rowSpecimen, ["MEDULLA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MEDULLA")
    if (Search_Term(rowSpecimen, ["SUBCAPSULAR SINUS"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "SUBCAPSULAR SINUS")
    if (Search_Term(rowSpecimen, ["PARACORT"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PARACORTEX")

    return row

def Musculature_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # MUSCULATURE #UBERON:0001015
    row['Sample_Type_Group'] = codetermExtract(dfAnatLocCodelist, "C13056")  # MUSCLE

    if (Search_Term(rowSpecimen, ["DELTOID"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C13050")  # MUSCLE, SKELETAL
        row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C13056")  # DELTOID MUSCLE
    elif (Search_Term(rowSpecimen, ["DIAPHRAGM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12702")  # MUSCLE, DIAPHRAGM
    elif (Search_Term(rowSpecimen, ["VASTUS"])):
        if (Search_Term(rowSpecimen, ["INTERMED"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C117876")  # MUSCLE, VASTUS INTERMEDIUS
        elif (Search_Term(rowSpecimen, ["LATERAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53073")  # MUSCLE, VASTUS LATERALIS
        elif (Search_Term(rowSpecimen, ["MEDIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C117736")  # MUSCLE, VASTUS MEDIALIS
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C13050")  # MUSCLE, SKELETAL
            row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MISPEC")
    elif (Search_Term(rowSpecimen, ["BICEPS FEMORIS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53147")  # MUSCLE, BICEPS FEMORIS
    elif (Search_Term(rowSpecimen, ["BICEPS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32200")  # MUSCLE, BICEPS BRACHII
    elif (Search_Term(rowSpecimen, ["GASTROCNEMIUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32666")  # MUSCLE, GASTROCNEMIUS
    elif (Search_Term(rowSpecimen, ["PECTORAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33286")  # MUSCLE, PECTORALIS
    elif (Search_Term(rowSpecimen, ["PLANTARIS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C117979")  # MUSCLE, PLANTARIS
    elif (Search_Term(rowSpecimen, ["PSOAS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33422")  # MUSCLE, PSOAS
    elif (Search_Term(rowSpecimen, ["QUADRICEPS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33441")  # MUSCLE, QUADRICEPS FEMORIS
    elif (Search_Term(rowSpecimen, ["RECTUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53175")  # MUSCLE, RECTUS FEMORIS
    elif (Search_Term(rowSpecimen, ["SOLEUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53075")  # MUSCLE, SOLEUS
    elif (Search_Term(rowSpecimen, ["TRICEPS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C90604")  # MUSCLE, TRICEPS BRACHII
    elif (Search_Term(rowSpecimen, ["SKELETAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C13050")  # MUSCLE, SKELETAL
    elif (Search_Term(rowSpecimen, ["SMOOTH"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12437")  # MUSCLE, SMOOTH
    elif (Search_Term(rowSpecimen, ["ABDOMIN"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32040")  # MUSCLE, ABDOMINAL
    elif (Search_Term(rowSpecimen, ["ADDUCTOR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53039")  # MUSCLE, ADDUCTOR
    elif (Search_Term(rowSpecimen, ["BULBO"])):
        if (Search_Term(rowSpecimen, ["LEVATOR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C112430")  # MUSCLE, LEVATOR ANI/BULBOSPONGIOSUS
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C112234")  # MUSCLE, BULBOSPONGIOSUS
    elif (Search_Term(rowSpecimen, ["LEVATOR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32984")  # MUSCLE, LEVATOR ANI
    elif (Search_Term(rowSpecimen, ["SPINAE", "SACROSPINAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52902")  # MUSCLE, ERECTOR SPINAE
    elif (Search_Term(rowSpecimen, ["DIGIT", "EXTENSOR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C52918")  # MUSCLE, EXTENSOR DIGITORUM LONGUS
    elif (Search_Term(rowSpecimen, ["OCULAR", "OCULO"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33199")  # MUSCLE, EXTRAOCULAR
    elif (Search_Term(rowSpecimen, ["GLUTEUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C78205")  # MUSCLE, GLUTEUS
    elif (Search_Term(rowSpecimen, ["INTERCOSTAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32824")  # MUSCLE, INTERCOSTAL
    elif (Search_Term(rowSpecimen, ["LATISSIMUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33150")  # MUSCLE, LATISSIMUS
    elif (Search_Term(rowSpecimen, ["MASSETER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C13074")  # MUSCLE, MASSETER
    elif (Search_Term(rowSpecimen, ["SEMITENDINOSUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53176")  # MUSCLE, SEMITENDINOSUS
    elif (Search_Term(rowSpecimen, ["STERNOCEPHALICUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C117980")  # MUSCLE, STERNOCEPHALICUS
    elif (Search_Term(rowSpecimen, ["TIBIALIS", "CRANIALIS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C53079")  # MUSCLE, TIBIALIS CRANIALIS
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C13050")  # MUSCLE, SKELETAL
        row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MISPEC")

    ###EXTRA TERMS
    if (Search_Term(row, ["MYOFIBER"])):
        row['Celltype'] = Add_Term('Celltype', row, "MYOCYTE")
    if (Search_Term(row, ["SARCOLEMMA"])):
        row['Celltype'] = Add_Term('Celltype', row, "MYOCYTE")
        row['Sub-Cellular_Locator'] = Add_Term('Sub-Cellular_Locator', row, "SARCOLEMMA")
    if (Search_Term(row, ["SARCOPLASM"])):
        row['Celltype'] = Add_Term('Celltype', row, "MYOCYTE")
        row['Sub-Cellular_Locator'] = Add_Term('Sub-Cellular_Locator', row, "SARCOPLASM")

    return row

def Pancreas_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12393")  # PANCREAS
    row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12393")  # PANCREAS

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["EXOCRINE"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32546"))  # PANCREAS, EXOCRINE
        if (Search_Term(row, ["EXOCRINE CELL"])):
            row['Celltype'] = Add_Term('Celltype', row, "EXOCRINE CELL")
    if (Search_Term(rowSpecimen, ["HEAD"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12269"))  # PANCREAS, HEAD
    if (Search_Term(rowSpecimen, ["TAIL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12271"))  # PANCREAS, TAIL

    ##CELLTYPES
    if (Search_Term(row, ["ACINAR", "ACINI", "ACINUS"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32546"))  # PANCREAS, EXOCRINE
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ACINAR")
        if (Search_Term(row, ["ACINAR CELL"])):
            row['Celltype'] = Add_Term('Celltype', row, "ACINAR CELL")
    if (Search_Term(rowSpecimen, ["DUCT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12272")  # DUCT, PANCREATIC
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32546"))  # PANCREAS, EXOCRINE
        if (Search_Term(rowSpecimen, ["DUCTAL CELL"])):
            row['Celltype'] = Add_Term('Celltype', row, "DUCTAL CELL")
    if (Search_Term(rowSpecimen, ["ISLET", "ENDOCRINE", "INSULAR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12608")  # PANCREAS, ENDOCRINE
        if (Search_Term(rowSpecimen, ["PERI ISLET", "PERI INSULAR", "PERIISLET", "PERIINSULAR"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PERIISLET")
        else:
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ISLET")
        if (Search_Term(row, ["ISLET CELL", "ENDOCRINE CELL"])):
            row['Celltype'] = Add_Term('Celltype', row, "ISLET CELL")
    if (Search_Term(rowSpecimen, ["SMALL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "SMALL CELL")

    return row

def Spleen_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12432")  # SPLEEN
    row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12432")  # SPLEEN

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["HILUM", "HILAR"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33601")) # SPLEEN, HILUM
    if (Search_Term(rowSpecimen, ["MARGINAL ZONE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "MARGINAL ZONE")
    if (Search_Term(rowSpecimen, ["MANTLE ZONE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MANTLE ZONE")
        #row['UBERON_Code'] = "UBERON:0010421"  # spleen B cell corona
    if (Search_Term(rowSpecimen, ["PULP"])):
        if (Search_Term(rowSpecimen, ["RED PULP"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "RED PULP")
        elif (Search_Term(rowSpecimen, ["WHITE PULP"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  "WHITE PULP")
    if (Search_Term(rowSpecimen, ["PERIARTERIOLAR LYMPHIOD SHEATH", "PALS"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "PERIARTERIOLAR LYMPHOID SHEATH (PALS)")
    if (Search_Term(rowSpecimen, ["ARTERY"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33597"))  # SPLENIC ARTERY
    if (Search_Term(rowSpecimen, ["VEIN"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33608"))  # SPLENIC VEIN


    return row

def Thymus_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12433") # THYMUS
    row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12433") # THYMUS

    if (Search_Term(rowSpecimen, ["CORTICOMEDULLARY JUNCTION"])):
        row['Anatomical_Region_of_Specimen'] = "CORTICOMEDULLARY JUNCTION"
    if (Search_Term(rowSpecimen, ["CORTEX", "CORTICAL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CORTEX")
    if (Search_Term(rowSpecimen, ["MEDULLA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MEDULLA")
    if (Search_Term(rowSpecimen, ["SUBCAPSUL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "SUBCAPSULE")
    elif (Search_Term(rowSpecimen, ["EXTRACAPSUL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "EXTRACAPSULAR")
    elif (Search_Term(rowSpecimen, ["CAPSUL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CAPSULE")

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["HASSAL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MEDULLA")
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "HASSAL'S CORPUSCLES")
    if (Search_Term(rowSpecimen, ["VASCUL", "VESSEL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12679"))  # BLOOD VESSEL
    if (Search_Term(rowSpecimen, ["PLEUR"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12469"))  # PLEURA

    return row

    ###############LESS COMMON####################

def Mouth_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):

    if (Search_Term(rowSpecimen, ["SALIVARY"])):
        row['Sample_Type_Group'] = "GLAND" # EXTENSIBLE
        if (Search_Term(rowSpecimen, ["PAROTID"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12427") # GLAND, SALIVARY, PAROTID
        elif (Search_Term(rowSpecimen, ["MANDIB", "SUBMAXILLA"])):
            if (Search_Term(rowSpecimen, ["SUBLINGUAL"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92215")  # GLAND, SALIVARY, SUBMANDIBULAR/GLAND, SALIVARY, SUBLINGUAL
            else:
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12233") # GLAND, SALIVARY, SUBMANDIBULAR
            if ((Search_Term(row, ["GRANULAR DUCT", "SECRETORY GRANULE"])) and (("RAT" in row.at['SPECIES']) or ("MOUSE" in row.at['SPECIES']))):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GRANULAR DUCT")
        elif (Search_Term(rowSpecimen, ["SUBLINGUAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12234") # GLAND, SALIVARY, SUBLINGUAL
        elif (Search_Term(rowSpecimen, ["MUCOUS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33141") # GLAND, SALIVARY, MUCOUS
        elif (Search_Term(rowSpecimen, ["SEROUS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33539") # GLAND, SALIVARY, SEROUS
        elif (Search_Term(rowSpecimen, ["ZYGOMATIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77624") # GLAND, SALIVARY, ZYGOMATIC
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12426") # GLAND, SALIVARY
    elif (Search_Term(rowSpecimen, ["TONGUE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12422") # TONGUE
        row['Sample_Type_Group'] = "OTHER" # EXTENSIBLE
    elif (Search_Term(rowSpecimen, ["TOOTH"])):
        row['Sample_Type_Group'] = "OTHER" # EXTENSIBLE
        if (Search_Term(rowSpecimen, ["INCISOR"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32769") # TOOTH, INCISOR
        elif (Search_Term(rowSpecimen, ["CANINE"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32258") # TOOTH, CANINE
        elif (Search_Term(rowSpecimen, ["MOLAR"])):
            if (Search_Term(rowSpecimen, ["PREMOLAR"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32201")  # TOOTH, PREMOLAR
            else:
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33136") # TOOTH, MOLAR
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12506") # TOOTH
    elif (Search_Term(rowSpecimen, ["GINGIVA"])):
        row['Sample_Type_Group'] = "OTHER" # EXTENSIBLE
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32677") # GINGIVA
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12421") # BODY CAVITY, ORAL
        row['Sample_Type_Group'] = "OTHER" # EXTENSIBLE

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["PERIODONTAL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32677")) # GINGIVA
    elif (Search_Term(rowSpecimen, ["LIGAMENT"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PERIODONTAL LIGAMENT")
    if (Search_Term(rowSpecimen, ["CEMENTUM"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CEMENTUM")
    if (Search_Term(rowSpecimen, ["PULP"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PULP")

    ##CELLTYPES
    if (Search_Term(rowSpecimen, ["AMELOBLAST"])):
        row['Celltype'] = Add_Term('Celltype', row, "AMELOBLAST")

    return row

def Thyroid_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):

    if ((Search_Term(rowSpecimen, ["PARATHYROID"])) and (Search_Term(rowSpecimen, ["GLAND THYROID"]))):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77667")  # GLAND, THYROID/GLAND, PARATHYROID
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12765")  # GLAND, PARATHYROID
    elif (Search_Term(rowSpecimen, ["PARATHYROID"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12765")  # GLAND, PARATHYROID
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12765")  # GLAND, PARATHYROID
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12400")  # GLAND, THYROID
        row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12400")  # GLAND, THYROID

    ##EXTRA TERMS
    if (Search_Term(row, ["FOLLICLE", "FOLLICULAR"])):
        if (Search_Term(row, ["FOLLICULAR CELL"])):
            if (Search_Term(row, ["PARAFOLLICULAR CELL"])):
                row['Celltype'] = "PARAFOLLICULAR CELL"
            else:
                row['Celltype'] = "FOLLICULAR CELL"
        else:
            if (Search_Term(row, ["LYMPHOID"])):
                row['MISTRESC'] = Add_Term('MISTRESC', row, "LYMPHOID FOLLICLE") #EXTENSIBLE
            else:
                row['Anatomical_Region_of_Specimen'] = "FOLLICLE"

    if (Search_Term(rowSpecimen, ["SUBCAPSUL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "SUBCAPSULE")
    elif (Search_Term(rowSpecimen, ["EXTRACAPSUL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                              "EXTRACAPSULAR")
    elif (Search_Term(rowSpecimen, ["CAPSULE"])):
        row['Anatomical_Region_of_Specimen'] = "CAPSULE"

    ##CELLTYPES
    if ((Search_Term(row, ["C CELL", "PARAFOLLICULAR CELL"])) and not (Search_Term(row, ["TIC CELL","ILIC CELL"]))):
        row['Celltype'] = "PARAFOLLICULAR CELL"
    if (Search_Term(rowSpecimen, ["CHIEF CELL"])):
        row['Celltype'] = "CHIEF CELL"

    return row

def Gland_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = "GLAND"  # EXTENSIBLE

    if (Search_Term(rowSpecimen, ["SALIVARY"])):
        row['Sample_Type_Group'] = "GLAND" # EXTENSIBLE
        if (Search_Term(rowSpecimen, ["PAROTID"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12427") # GLAND, SALIVARY, PAROTID
        elif (Search_Term(rowSpecimen, ["MANDIB", "SUBMAXILLA"])):
            if (Search_Term(rowSpecimen, ["SUBLINGUAL"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92215")  # GLAND, SALIVARY, SUBMANDIBULAR/GLAND, SALIVARY, SUBLINGUAL
            else:
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12233") # GLAND, SALIVARY, SUBMANDIBULAR
            if ((Search_Term(row, ["GRANULAR DUCT", "SECRETORY GRANULE"])) and (("RAT" in row.at['SPECIES']) or ("MOUSE" in row.at['SPECIES']))):
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GRANULAR DUCT")
        elif (Search_Term(rowSpecimen, ["SUBLINGUAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12234") # GLAND, SALIVARY, SUBLINGUAL
        elif (Search_Term(rowSpecimen, ["MUCOUS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33141") # GLAND, SALIVARY, MUCOUS
        elif (Search_Term(rowSpecimen, ["SEROUS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33539") # GLAND, SALIVARY, SEROUS
        elif (Search_Term(rowSpecimen, ["ZYGOMATIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77624") # GLAND, SALIVARY, ZYGOMATIC
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12426") # GLAND, SALIVARY
    elif (Search_Term(rowSpecimen, ["PREPUTIAL"])):
        if (Search_Term(rowSpecimen, ["CLITORAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C117978")  # GLAND, PREPUTIAL/GLAND, CLITORAL
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C79432")  # GLAND, PREPUTIAL
    elif (Search_Term(rowSpecimen, ["THIRD EYELID"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77616")  # GLAND OF THE THIRD EYELID
    elif (Search_Term(rowSpecimen, ["AMPULLARY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77955")  # GLAND, AMPULLARY
    elif (Search_Term(rowSpecimen, ["ANAL SAC"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C125895")  # GLAND, ANAL SAC
    elif (Search_Term(rowSpecimen, ["BRUNNER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C13010")  # GLAND, BRUNNER'S
    elif (Search_Term(rowSpecimen, ["BULBOURETHRAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32395")  # GLAND, BULBOURETHRAL
    elif (Search_Term(rowSpecimen, ["CIRCUMANAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77610")  # GLAND, CIRCUMANAL
    elif (Search_Term(rowSpecimen, ["CLITORAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77617")  # GLAND, CLITORAL
    elif (Search_Term(rowSpecimen, ["COAGULATING"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77618")  # GLAND, COAGULATING
    elif (Search_Term(rowSpecimen, ["ENDOMETRI"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33842")  # GLAND, ENDOMETRIAL
    elif (Search_Term(rowSpecimen, ["HARDERIAN"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77619")  # GLAND, HARDERIAN
    elif (Search_Term(rowSpecimen, ["LACRIMAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12346")  # GLAND, LACRIMAL
    elif (Search_Term(rowSpecimen, ["MEIBOMIAN"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33075")  # GLAND, MEIBOMIAN
    elif (Search_Term(rowSpecimen, ["PERIANAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77620")  # GLAND, PERIANAL
    elif (Search_Term(rowSpecimen, ["PINEAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12398")  # GLAND, PINEAL
    elif (Search_Term(rowSpecimen, ["SEBACEOUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33519")  # GLAND, SEBACEOUS
    elif (Search_Term(rowSpecimen, ["ZYMBAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77954")  # GLAND, ZYMBAL
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C70699")  # BIOSPECIMEN
        row['Anatomical_Region_of_Specimen'] = "GLAND"  # EXTENSIBLE

    if (Search_Term(rowSpecimen, ["ACINAR CELL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ACINAR")
        row['Celltype'] = Add_Term('Celltype', row, "ACINAR CELL")

    return row

def Eye_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist):
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12401")  # EYE

    if (Search_Term(rowSpecimen, ["EYELID"])):
        row['Sample_Type_Group'] = "OTHER" # EXTENSIBLE
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12713")  # EYELID
    elif (Search_Term(rowSpecimen, ["LACRIMAL"])):
        row['Sample_Type_Group'] = "GLAND"  # EXTENSIBLE
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12346")  # GLAND, LACRIMAL
    elif (Search_Term(rowSpecimen, ["ANTERIOR CHAMBER", "ANTERIOR COMPART"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12667")  # EYE, ANTERIOR CHAMBER
    elif (Search_Term(rowSpecimen, ["ANTERIOR SEGMENT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12668")  # EYE, ANTERIOR SEGMENT
    elif (Search_Term(rowSpecimen, ["AQUEOUS HUMOR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C13190")  # EYE, AQUEOUS HUMOR
    elif (Search_Term(rowSpecimen, ["POSTERIOR CHAMBER", "POSTERIOR COMPART"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12900")  # EYE, POSTERIOR CHAMBER
    elif (Search_Term(rowSpecimen, ["CILIARY BODY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12345")  # EYE, CILIARY BODY
    elif (Search_Term(rowSpecimen, ["CORNEA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12342")  # EYE, CORNEA
    elif (Search_Term(rowSpecimen, ["CHOROID"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12344")  # EYE, CHOROID
    elif (Search_Term(rowSpecimen, ["IRIS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12737")  # EYE, IRIS
    elif (Search_Term(row, ["RETINA", "ROSETTE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C49328")  # EYE, RETINA
        if ((row.str.contains("EPITHEL", "RPE", regex=True, na=False).any())):
            row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C33470")  # RETINAL PIGMENTED EPITHELIAL LAYER
            row['Celltype'] = Add_Term('Celltype', row, "RETINAL PIGMENTED EPITHELIAL CELL")
    elif (Search_Term(row, ["RPE", "PIGMENTED EPITHELIAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C49328")  # EYE, RETINA
        row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C33470")  # RETINAL PIGMENTED EPITHELIAL LAYER
        row['Celltype'] = Add_Term('Celltype', row, "RETINAL PIGMENTED EPITHELIAL CELL")
    elif (Search_Term(rowSpecimen, ["SCLERA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12784")  # EYE, SCLERA
    elif (Search_Term(rowSpecimen, ["VITREOUS", "POSTREMAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33884")  # EYE, VITREOUS
        if (Search_Term(rowSpecimen, ["CHAMBER"])):
            row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C33885")  # EYE, VITREOUS CHAMBER
    elif (Search_Term(rowSpecimen, ["LENS", "LENTICULAR"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12743")  # EYE, LENS
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12401")  # EYE

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["OUTER NUCLEAR LAYER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C49328")  # EYE, RETINA
        row['Anatomical_Region_of_Specimen'] = "OUTER NUCLEAR LAYER"
    if (Search_Term(rowSpecimen, ["PHOTORECEPTOR LAYER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C49328")  # EYE, RETINA
        row['Anatomical_Region_of_Specimen'] = "PHOTORECEPTOR LAYER"
    if (Search_Term(rowSpecimen, ["OCULAR MUSCLE"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33199"))  # EXTRAOCULAR MUSCLE
    if (Search_Term(rowSpecimen, ["CONJUNCTIVA"])):
        if (Search_Term(rowSpecimen, ["PALPEBRAL"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12901"))  # PALPEBRAL CONJUNCTIVA
        elif (Search_Term(rowSpecimen, ["BULBA"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12902"))  # BULBAR CONJUNCTIVA
        else:
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12341"))  # CONJUNCTIVA
    if (Search_Term(rowSpecimen, ["BULBA"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12902"))  # BULBAR CONJUNCTIVA
    if (Search_Term(rowSpecimen, ["PALPEBRAL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12901"))  # PALPEBRAL CONJUNCTIVA
    if (Search_Term(row, ["GANGLION"])):
        if (Search_Term(rowSpecimen, ["CELL LAYER"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "GANGLION CELL LAYER")
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C49328")  # EYE, RETINA
            row['Celltype'] = Add_Term('Celltype', row, "RETINAL GANGLION CELL")
    if (Search_Term(rowSpecimen, ["HARDERIAN"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77619")  # GLAND, HARDERIAN
    if (Search_Term(row, ["UVEITIS"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12811"))  # UVEA
        row['MISTRESC'] = Add_Term('MISTRESC', row, codetermExtract(dfNonNeoplasmCodelist, "C3137"))  # INFLAMMATION")
    if (Search_Term(rowSpecimen, ["UVEA"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12811"))  # UVEA
    if (Search_Term(row, ["NUCLEAR LAYER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C49328")  # EYE, RETINA
        if (Search_Term(row, ["INNER"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "INNER NUCLEAR LAYER")
        if (Search_Term(row, ["OUTER"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "OUTER NUCLEAR LAYER")
    if (Search_Term(row, ["PHOTORECEPTOR"])):
        row['Celltype'] = Add_Term('Celltype', row, "PHOTORECEPTOR CELL")
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C49328")  # EYE, RETINA
        if (Search_Term(row, ["INNER"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "INNER NUCLEAR LAYER")
        if (Search_Term(row, ["OUTER"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "OUTER NUCLEAR LAYER")
    if (Search_Term(row, ["LENS FIBER", "CATARACT"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "LENS FIBER")
    if (Search_Term(row, ["OPTIC NERVE"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12761")) # NERVE, OPTIC
    if (Search_Term(row, ["OPTIC DISK", "OPTIC DISC"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12760")) # OPTIC DISC
    if (Search_Term(row, ["TAPETAL", "TAPETUM"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "TAPETUM")
    if (Search_Term(rowSpecimen, ["LIMBUS"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "LIMBUS")
    if (Search_Term(rowSpecimen, ["MUSCLE"])):
        if (Search_Term(rowSpecimen, ["CILIARY MUSCLE"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CILIARY MUSCLE")
        else:
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33199"))  # EXTRAOCULAR MUSCLE
    if (Search_Term(rowSpecimen, ["HYALOID"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "HYALOID ARTERY")

    ##CELLTYPES
    if (Search_Term(row, ["GLIOSIS", "GLIAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "GLIAL CELL")
    if (Search_Term(row, ["RODS"])):
        row['Celltype'] = Add_Term('Celltype', row, "PHOTORECEPTOR CELL")
        row['Celltype'] = Add_Term('Celltype', row, "RODS")
    if (Search_Term(row, ["CONES"])):
        row['Celltype'] = Add_Term('Celltype', row, "PHOTORECEPTOR CELL")
        row['Celltype'] = Add_Term('Celltype', row, "CONES")

    return row

def Injection_Application_Site_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = "INJECTION/APPLICATION SITE"  #EXTENSIBLE

    if (Search_Term(rowSpecimen, ["TAIL VEIN", "CAUDAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92598")  # VEIN, CAUDAL
    elif (Search_Term(rowSpecimen, ["SURGICAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77683")  # SITE, SURGICAL
    elif (Search_Term(rowSpecimen, ["IMPLANT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77678")  # SITE, IMPLANTATION
    elif (Search_Term(rowSpecimen, ["INFUSION"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77679")  # SITE, INFUSION
    elif (Search_Term(rowSpecimen, ["CATHETER"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92596")  # SITE, CATHETER
    elif (Search_Term(rowSpecimen, ["EXTERIORIZATION"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77685")  # SITE, EXTERIORIZATION
    elif (Search_Term(rowSpecimen, ["INJURY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77681")  # SITE, INJURY
    elif (Search_Term(rowSpecimen, ["INJECTION"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77680")  # SITE, INJECTION
    elif (Search_Term(rowSpecimen, ["BIOPSY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77677")  # SITE, BIOPSY
    elif (Search_Term(rowSpecimen, ["MICROCHIP"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77682")  # SITE, MICROCHIP
    elif (Search_Term(rowSpecimen, ["SUBCUTANEOUS PORT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C147512")  # SITE, SUBCUTANEOUS PORT
    elif (Search_Term(rowSpecimen, ["TATTOO"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77684")  # SITE, TATTOO
    elif (Search_Term(rowSpecimen, ["APPLICATION"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77676")  # SITE, APPLICATION
    elif (Search_Term(rowSpecimen, ["UNCERTAIN PRIMARY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C48322")  # SITE, UNCERTAIN PRIMARY
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77676")  # SITE, APPLICATION
        row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MISPEC")

        ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["CEPHALIC"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32286"))  # CEPHALIC VEIN
    if (Search_Term(rowSpecimen, ["SAPHEN"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33511"))  # SAPHENOUS VEIN
    if (Search_Term(rowSpecimen, ["CRANIAL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C77638"))  # CRANIAL CAVITY
    '''
    if (Search_Term(rowSpecimen, ["LOWER"])):
        row['MIDIR'] = Add_Term('MIDIR', row, codetermExtract(dfCodelist, "C25309"))  # LOWER
    if (Search_Term(rowSpecimen, ["UPPER"])):
        row['MIDIR'] = Add_Term('MIDIR', row, codetermExtract(dfCodelist, "C25355"))  # UPPER
    if (Search_Term(rowSpecimen, ["DORSAL"])):
        row['MIDIR'] = Add_Term('MIDIR', row, codetermExtract(dfCodelist, "C45874"))  # DORSAL
    if (Search_Term(rowSpecimen, ["POSTERIOR"])):
        row['MIDIR'] = Add_Term('MIDIR', row, codetermExtract(dfCodelist, "C25622"))  # POSTERIOR
    '''
    if (Search_Term(rowSpecimen, ["LUMBAR"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12744"))  # LUMBAR REGION
    if (Search_Term(rowSpecimen, ["SCAPULA"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12783"))  # SCAPULA
    if (Search_Term(row, ["BLOOD PRESSURE DEVICE"])):
        row['Extra_finding_terms'] = Add_Term('Extra_finding_terms', row, "BLOOD PRESSURE DEVICE")
    if (Search_Term(rowSpecimen, ["PELVIS"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12767"))  # PELVIS
    if (Search_Term(rowSpecimen, ["THIGH"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33763"))  # THIGH

    '''
    ##EXTRACTING SITE INFO
    if ("ADMINISTRATION SITE" in row.at['Sample Type']):  # info in Sample Type for ADMINISTRATION SITE
        my_string = row.at['Sample Type']
        new_string = my_string.split("SITE", 1)
        if (len(new_string) > 1):
            SITE_new_string = "SITE: " + new_string[1]
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  SITE_new_string)
    elif ("SITE" in row.at['Sample Type']):  # info in Locator column for all other SITE
        my_string = row.at['Locator']
        if ("DTH" in row.at['Locator']):
            SITE_new_string = "SITE: " + my_string
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  SITE_new_string)
        elif ("KLH" in row.at['Locator']):
            SITE_new_string = "SITE: KLH "
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  SITE_new_string)
        elif re.match("PBS [0-9]+", row.at['Locator']) is not None:
            SITE_new_string = "SITE: " + my_string
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  SITE_new_string)
        elif re.match("SC[0-9]+", row.at['Locator']) is not None:
            SITE_new_string = "SITE: " + my_string
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                  SITE_new_string)
        elif re.match("SITE [0-9]+", row.at['Locator']) is not None:
            new_string = my_string.split("SITE ", 1)
            if (len(new_string) > 1):
                SITE_new_string = "SITE: " + new_string[1]
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      SITE_new_string)
        elif re.match("SUBCUTANEOUS [0-9]+", row.at['Locator']) is not None:
            new_string = my_string.split("SUBCUTANEOUS ", 1)
            if (len(new_string) > 1):
                SITE_new_string = "SITE: SC " + new_string[1]
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      SITE_new_string)
        elif ("SUBCUTANEOUS SITE" in row.at['Locator']):
            new_string = my_string.split("SITE ", 1)
            if (len(new_string) > 1):
                SITE_new_string = "SITE: SC " + new_string[1]
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      SITE_new_string)
        elif re.match("SITE [0-9]+", row.at['ALL_result_comments']) is not None:
            new_string = my_string.split("SITE ", 1)
            if (len(new_string) > 1):
                SITE_new_string = "SITE: " + new_string[1]
                row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row,
                                                                      SITE_new_string)
    '''
    return row

def Male_Repro_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = codetermExtract(dfAnatLocCodelist, "C12722")  # MALE REPRODUCTIVE SYSTEM

    if (Search_Term(rowSpecimen, ["SPERMATIC CORD"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12412")  # TESTIS
        row['Anatomical_Region_of_Specimen'] = "SPERMATIC CORD"
    elif (Search_Term(row, ["SPERMATOGENIC CYCLE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12412")  # TESTIS
        row['Extra_finding_terms'] = "SPERMATOGENIC CYCLE ASSESSMENT"
    elif (Search_Term(rowSpecimen, ["SEMINAL VESICLE", "SEMINAL SAC"])):
        if (Search_Term(rowSpecimen, ["COAGULATING"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92216")  # GLAND, SEMINAL VESICLE/GLAND, COAGULATING
        elif (Search_Term(rowSpecimen, ["PROSTATE"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77670")  # GLAND, PROSTATE/GLAND, SEMINAL VESICLE
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12787")  # GLAND, SEMINAL VESICLE"
    elif (Search_Term(row, ["PENIS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12409")  # PENIS
        if (Search_Term(rowSpecimen, ["CORPUS", "BODY"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12325"))  # PENIS, BODY
        elif (Search_Term(rowSpecimen, ["GLANS"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12324"))  # PENIS, GLANS
        elif (Search_Term(rowSpecimen, ["RADIX"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C124350"))  # PENIS, RADIX
    elif (Search_Term(row, ["SCROTUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12785")  # SCROTUM
    elif (Search_Term(row, ["DEFERENS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12412")  # TESTIS
        row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12813") # VAS DEFERENS
    elif (Search_Term(rowSpecimen, ["TESTIS"])):
        if (Search_Term(rowSpecimen, ["EPIDIDYM"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77668")  # TESTIS/EPIDIDYMIS
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12412")  # TESTIS
    elif (Search_Term(row, ["EPIDIDYM"])):
        if (Search_Term(rowSpecimen, ["CAUDA", "TAIL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33732")  # EPIDIDYMIS, CAUDA
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12328")  # EPIDIDYMIS
        if (Search_Term(rowSpecimen, ["CAPUT", "HEAD"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CAPUT")
        if (Search_Term(rowSpecimen, ["CORPUS", "BODY"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CORPUS")
    elif (Search_Term(rowSpecimen, ["PROSTATE"])):
        if (Search_Term(row, ["DORSOLATERAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77622")  # GLAND, PROSTATE DORSOLATERAL
        elif (Search_Term(row, ["VENTRAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77623")  # GLAND, PROSTATE VENTRAL
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12410")  # GLAND, PROSTATE
    elif (Search_Term(row, ["EFFERENT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32492")  # DUCT, EFFERENT
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C70699")  # BIOSPECIMEN
        row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MIANTREG")
        codetermExtract(dfAnatLocCodelist, "C12722")  # MALE REPRODUCTIVE SYSTEM

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["SEMINIFEROUS"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "SEMINIFEROUS TUBULE")
    if (Search_Term(rowSpecimen, ["SEX CORD"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "SEX CORD")
    if (Search_Term(rowSpecimen, ["RETE TEST"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "RETE TESTIS")

    ##CELLTYPES
    if (Search_Term(rowSpecimen, ["SERTOLI"])):
        row['Celltype'] = Add_Term('Celltype', row, "SERTOLI CELL")
    if (Search_Term(rowSpecimen, ["INTERSTITIAL CELL", "LEYDIG"])):
        row['Celltype'] = Add_Term('Celltype', row, "LEYDIG CELL")
    if (Search_Term(rowSpecimen, ["GERM CELL", "GERMINAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "GERM CELL")
    if (Search_Term(rowSpecimen, ["SPERM"])):
        if ((Search_Term(rowSpecimen, ["SPERMATID"])) and (Search_Term(rowSpecimen, ["GIANT CELL"]))):
            row['Celltype'] = Add_Term('Celltype', row, "SPERMATID GIANT CELL")
        elif (Search_Term(rowSpecimen, ["SPERMATID"])):
            row['Celltype'] = Add_Term('Celltype', row, "SPERMATID")
        elif (Search_Term(rowSpecimen, ["SPERMATOGONIA"])):
            row['Celltype'] = Add_Term('Celltype', row, "SPERMATOGONIA")
        else:
            row['Celltype'] = Add_Term('Celltype', row, "SPERM")
    if (Search_Term(rowSpecimen, ["ACINAR CELL"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ACINAR")
        row['Celltype'] = Add_Term('Celltype', row, "ACINAR CELL")

    return row

def Female_Repro_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    # FEMALE REPRO #UBERON:0000474
    row['Sample_Type_Group'] = codetermExtract(dfAnatLocCodelist, "C12402")  # FEMALE REPRODUCTIVE SYSTEM

    if (Search_Term(rowSpecimen, ["REPRODUCTIVE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C70699")  # BIOSPECIMEN
        if (row['SEX'] == "F"):
            row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12402")  # FEMALE REPRODUCTIVE SYSTEM
        elif (row['SEX'] == "M"):
            row['Sample_Type_Group'] = codetermExtract(dfAnatLocCodelist, "C12722")  # MALE REPRODUCTIVE SYSTEM
            row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12722")  # MALE REPRODUCTIVE SYSTEM
        else:
            row['Sample_Type_Group'] = "REPRODUCTIVE SYSTEM"
            row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MIANTREG")
    elif (Search_Term(rowSpecimen, ["MAMMARY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12367")  # GLAND, MAMMARY
        row['Sample_Type_Group'] = "GLAND"  # EXTENSIBLE
        if (row['SEX'] == "M"):
            row['Sample_Type_Group'] = codetermExtract(dfAnatLocCodelist, "C12722")  # MALE REPRODUCTIVE SYSTEM
    elif (Search_Term(rowSpecimen, ["FALLOPIAN", "OVIDUCT"])):
        if (Search_Term(rowSpecimen, ["OVARY"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92595")  # OVARY/OVIDUCT
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12403")  # OVIDUCT
    elif (Search_Term(rowSpecimen, ["CERVIX"])):
        if (Search_Term(rowSpecimen, ["UTERUS", "UTERINE"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92436")  # UTERUS/CERVIX
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12311")  # CERVIX
        if ((Search_Term(rowSpecimen, ["PLUG"])) and (Search_Term(rowSpecimen, ["MUC"]))):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MUCOUS PLUG")
    elif (Search_Term(rowSpecimen, ["UTERUS", "UTERINE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12405")  # UTERUS
    elif (Search_Term(rowSpecimen, ["OVARY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12404")  # OVARY
        if (Search_Term(rowSpecimen, ["LUTEUM", "LUTEA", "CORPUS"])):
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C26465")) # CORPUS LUTEUM
    elif (Search_Term(rowSpecimen, ["VAGINA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12407")  # VAGINA
    elif (Search_Term(rowSpecimen, ["ESTRUS"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77665")  # WHOLE ANIMAL
        row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12402")  # FEMALE REPRODUCTIVE SYSTEM
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77665")  # WHOLE ANIMAL
        row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12402")  # FEMALE REPRODUCTIVE SYSTEM
        row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MIANTREG")

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["FOLLIC"])):
        if (Search_Term(rowSpecimen, ["ANTRAL"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "ANTRAL FOLLICLE")
        elif (Search_Term(rowSpecimen, ["TERTIARY"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "TERTIARY FOLLICLE")
        else:
            row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33244"))  # OVARIAN FOLLICLE
    if (Search_Term(row, ["ATRESIA"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33244"))  # OVARIAN FOLLICLE
    if (Search_Term(rowSpecimen, ["MYOMETRI", "MYOMETRAL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12314"))  # MYOMETRIUM
    if (Search_Term(rowSpecimen, ["ENDOMETRI"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12313"))  # ENDOMETRIUM
    if (Search_Term(row, ["RETE OVARII"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "RETE OVARII")
    if (Search_Term(row, ["RETE TUBULE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "RETE TUBULE")
    if (Search_Term(rowSpecimen, ["SEX CORD"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "SEX CORD")
    if (Search_Term(rowSpecimen, ["BURSA"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "BURSA")
    if (Search_Term(rowSpecimen, ["INTERSTITIAL GLAND"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "INTERSTITIAL GLAND")
    if (Search_Term(rowSpecimen, ["PAROVARI", "MESONEPHRIC"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MESONEPHRIC DUCT")
    if (Search_Term(rowSpecimen, ["PAPILLA"])):
        if (Search_Term(rowSpecimen, ["PAPILLARY"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PAPILLARY")
        else:
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "PAPILLA")

    ##CELLTYPES
    if (Search_Term(rowSpecimen, ["SERTOLI"])):
        row['Celltype'] = Add_Term('Celltype', row, "SERTOLIFORM CELL")
    if (Search_Term(rowSpecimen, ["ENDOMETRIAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "ENDOMETRIAL CELL")
    if (Search_Term(rowSpecimen, ["GERM CELL", "GERMINAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "GERM CELL")
    if (Search_Term(rowSpecimen, ["GRANULOSA"])):
        row['Celltype'] = Add_Term('Celltype', row, "FOLLICULAR CELL")

    return row

def Integument_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92441")  # SKIN/SUBCUTIS
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C92441")  # SKIN/SUBCUTIS

    if (Search_Term(rowSpecimen, ["LIP"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12220")  # LIP
    elif (Search_Term(rowSpecimen, ["ARM"])):
        row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C32141")  # ARM
    elif (Search_Term(rowSpecimen, ["HEAD"])):
        row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12419")  # HEAD
    elif (Search_Term(rowSpecimen, ["TAIL"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77663"))  # TAIL
    elif (Search_Term(rowSpecimen, ["NOSE"])):
        row['MIANTREG'] = codetermExtract(dfAnatLocCodelist, "C12756")  # NOSE
    elif (Search_Term(rowSpecimen, ["PLANTAR"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32622"))  # FOOT
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33326"))  # SOLE
    elif (Search_Term(rowSpecimen, ["FOOT"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32622"))  # FOOT
    elif (Search_Term(rowSpecimen, ["HAND"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32712"))  # HAND
    elif (Search_Term(rowSpecimen, ["LEG"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C32974"))  # LEG
    elif (Search_Term(rowSpecimen, ["PAW"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C77660"))  # PAW
    elif (Search_Term(row, ["SCROTUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12785")  # SCROTUM
    elif (Search_Term(rowSpecimen, ["ISCHIAL CALLOSITY", "ISCHIAL PAD"])):
        row['Anatomical_Region_of_Specimen'] = "ISCHIAL CALLOSITY"
    elif (Search_Term(rowSpecimen, ["NASAL PLANUM"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12756"))  # NOSE
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "NASAL PLANUM")
    elif (Search_Term(rowSpecimen, ["PINNA"])):
        row['MISPEC'] = Add_Term('MISPEC', row, codetermExtract(dfSpecimenCodelist, "C12292"))  # PINNA
    elif (Search_Term(rowSpecimen, ["PALM"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33252"))  # PALM


    ###EXTRA TERMS
    '''
    if (Search_Term(rowSpecimen, ["LEFT"])):
        row['MILAT'] = "LEFT"
    if (Search_Term(rowSpecimen, ["RIGHT"])):
        row['MILAT'] = "RIGHT"
    if (Search_Term(rowSpecimen, ["LOWER"])):
        row['MIDIR'] = Add_Term('MIDIR', row, codetermExtract(dfCodelist, "C25309"))  # LOWER
    if (Search_Term(rowSpecimen, ["UPPER"])):
        row['MIDIR'] = Add_Term('MIDIR', row, codetermExtract(dfCodelist, "C25355"))  # UPPER
    if (Search_Term(rowSpecimen, ["DORSAL"])):
        row['MIDIR'] = Add_Term('MIDIR', row, codetermExtract(dfCodelist, "C45874"))  # DORSAL
    if (Search_Term(rowSpecimen, ["POSTERIOR"])):
        row['MIDIR'] = Add_Term('MIDIR', row, codetermExtract(dfCodelist, "C25622"))  # POSTERIOR
    if (Search_Term(rowSpecimen, ["VENTRAL", "VENTRUM"])):
        row['MIDIR'] = Add_Term('MIDIR', row, codetermExtract(dfCodelist, "C45875"))  # VENTRAL
    '''
    if (Search_Term(rowSpecimen, ["LUMBAR"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C45875"))  # LUMBAR REGION
    if (Search_Term(rowSpecimen, ["SCAPULA"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12783"))  # SCAPULA
    if (Search_Term(rowSpecimen, ["MUZZLE"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "MUZZLE")
    if (Search_Term(rowSpecimen, ["DERMIS", "DERMAL", "DERMOID"])):
        if (Search_Term(rowSpecimen, ["EPIDERMIS", "EPIDERMAL", "EPIDERMOID"])):
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "EPIDERMIS")
        else:
            row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "DERMIS")
    #if (Search_Term(rowSpecimen, ["HIND"])): # or ("HIND" in row.at['Specimen Directionality'])):
    #    row['MIDIR'] = Add_Term('MIDIR', row, "HIND")

    ##CELLTYPES
    if (Search_Term(rowSpecimen, ["GERM CELL", "GERMINAL CELL"])):
        row['Celltype'] = Add_Term('Celltype', row, "GERM CELL")

    return row

def Adipose_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C12472")  # ADIPOSE TISSUE

    if (Search_Term(rowSpecimen, ["BROWN"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C32235")  # ADIPOSE TISSUE, BROWN
    elif (Search_Term(rowSpecimen, ["WHITE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33889")  # ADIPOSE TISSUE, WHITE
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12472")  # ADIPOSE TISSUE

    ###EXTRA TERMS
    if (Search_Term(rowSpecimen, ["AXILLA"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12674"))  # AXILLA

    elif (Search_Term(rowSpecimen, ["INGUINAL"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12726"))  # INGUINAL REGION
    elif (Search_Term(rowSpecimen, ["SUBCUT"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33645"))  # SUBCUTIS
    elif (Search_Term(rowSpecimen, ["SCAPULA"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12783"))  # SCAPULA
    elif (Search_Term(rowSpecimen, ["OMENT"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33209"))  # OMENTUM
    else:
        row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MIANTREG")

    return row

def Nose_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = "OTHER"  # EXTENSIBLE

    if (Search_Term(rowSpecimen, ["TURBINATE", "CONCHA"])):
        if (Search_Term(rowSpecimen, ["CONCHA"])):
            if (Search_Term(rowSpecimen, ["DORSAL"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C139163")  # NASAL TURBINATE, DORSAL CONCHA
            elif (Search_Term(rowSpecimen, ["ETHMOIDAL"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C139162")  # NASAL TURBINATE, ETHMOIDAL CONCHA
            elif (Search_Term(rowSpecimen, ["MIDDLE"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C139164")  # NASAL TURBINATE, MIDDLE CONCHA
            elif (Search_Term(rowSpecimen, ["VENTRAL"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C139165")  # NASAL TURBINATE, VENTRAL CONCHA
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C49594")  # NASAL TURBINATE
    elif (Search_Term(rowSpecimen, ["NOSE"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12756")  # NOSE
    elif (Search_Term(rowSpecimen, ["MUCOSA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33205")  # MUCOSA, NASAL
    elif (Search_Term(row, ["OLFACTORY"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C98765")  # OLFACTORY REGION
    elif (Search_Term(rowSpecimen, ["NASOLACRIMAL DUCT"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C33161")  # DUCT, NASOLACRIMAL
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12424")  # BODY CAVITY, NASAL

    ##EXTRA TERMS
    if (Search_Term(rowSpecimen, ["TRABECULA"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C33157"))  # NASAL BONE
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "TRABECULA")
    if (Search_Term(row, ["TOOTH"])):
        row['MIANTREG'] = Add_Term('MIANTREG', row, codetermExtract(dfAnatLocCodelist, "C12506"))  # TOOTH
    if (Search_Term(row, ["CEMENTUM"])):
        row['Anatomical_Region_of_Specimen'] = Add_Term('Anatomical_Region_of_Specimen', row, "CEMENTUM")

    return row

def Whole_Body_SPECIMEN_Mapping(row, rowSpecimen, dfSpecimenCodelist, dfAnatLocCodelist):
    row['Sample_Type_Group'] = codetermExtract(dfSpecimenCodelist, "C77665")  # WHOLE ANIMAL

    if (Search_Term(rowSpecimen, ["ENDOTHELIUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77665")  # WHOLE ANIMAL
        row['Anatomical_Region_of_Specimen'] = "ENDOTHELIUM"
    elif (Search_Term(rowSpecimen, ["MEDIASTINUM"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12748")  # MEDIASTINUM
    elif (Search_Term(rowSpecimen, ["PERITONEUM", "SEROSA"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12770")  # PERITONEUM
    elif (Search_Term(rowSpecimen, ["WHOLE BODY", "WHOLE ANIMAL"])):
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77665")  # WHOLE ANIMAL
    elif (Search_Term(rowSpecimen, ["CAVITY"])):
        if (Search_Term(rowSpecimen, ["THORACIC"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12905")  # BODY CAVITY, THORACIC
        elif (Search_Term(rowSpecimen, ["CRANIAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77638")  # BODY CAVITY, CRANIAL
        elif (Search_Term(rowSpecimen, ["PERITONEAL"])):
            if (Search_Term(rowSpecimen, ["EXTRAPERITONEAL"])):
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C92208")  # BODY CAVITY, EXTRAPERITONEAL
            else:
                row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12769")  # BODY CAVITY, PERITONEAL
        elif (Search_Term(rowSpecimen, ["PELVIC" , "PELVIS"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12767")  # BODY CAVITY, PELVIC
        elif (Search_Term(rowSpecimen, ["ABDOMINAL"])):
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C12664")  # BODY CAVITY, ABDOMINAL
        else:
            row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C25444")  # BODY CAVITY
    else:
        row['MISPEC'] = codetermExtract(dfSpecimenCodelist, "C77665")  # WHOLE ANIMAL
        row['CODE_CONFLICT_FLAG'] = Add_Term('CODE_CONFLICT_FLAG', row, "MIANTREG")

    return row

##########################################################################
###                         Main function                               ###
###########################################################################
def main():
    ##In console: wget "https://evs.nci.nih.gov/ftp1/CDISC/SEND/SEND%20Terminology.txt"
    path_to_SENDcodelist = "C:\\Users\c143390\Desktop\SEND Terminology.txt"
    #dfSENDCodelist = pd.read_csv(path_to_SENDcodelist, delimiter="\t").fillna("")
    with open(path_to_SENDcodelist, 'rb') as f:
        result = chardet.detect(f.read())  # or readline if the file is large
    dfSENDCodelist = pd.read_csv(path_to_SENDcodelist, encoding=result['encoding'], delimiter="\t").drop_duplicates()
    # dfSENDCodelisttemp = dfClean(dfSENDCodelist)

    ##Directory to .xpt files
    #has MI
    mydir = os.path.abspath(r"\\ix1invivo-p\ivdr\InVivoFileFormatAdapter\Mappings\SEND\2021-02-03-11-00-02\done\8413493")
    #no MI
    #mydir = os.path.abspath(r"\\ix1invivo-p\ivdr\InVivoFileFormatAdapter\Mappings\SEND\2021-06-22-11-00-03\done\8461379")
    #outputFile = 'filepathXPT_test2.csv'

    ##Check which domains have .xpt files in directory
    fMI, fTS, fEX, fDM, fDS = Domain_test(mydir)

    ##Generate MI domain dataframe
    if fMI==None:
        sys.exit("No MI.xpt file detected.") #if no MI domain, then exit with error
    else:
        dfMIRaw = Raw_dataframe(fMI, fTS, fEX, fDM, fDS)

    dfMIRaw, dftemp = dfInitialization(dfMIRaw)
    dfSEXCodelist = dfSENDCodelist[dfSENDCodelist["Codelist Name"].isin(["Sex"])].drop_duplicates()
    dftemp = Gender_Mapping(dftemp, dfSEXCodelist)

    dfSTRAINCodelist = dfSENDCodelist[dfSENDCodelist["Codelist Name"].isin(["Strain/Substrain"])].drop_duplicates()
    dfSPECIESCodelist = dfSENDCodelist[dfSENDCodelist["Codelist Name"].isin(["Species"])].drop_duplicates()
    dftemp = Species_Mapping(dftemp, dfSPECIESCodelist, dfSTRAINCodelist)

    dfSEVERITYCodelist = dfSENDCodelist[dfSENDCodelist["Codelist Name"].isin(["SEND Severity"])].drop_duplicates()
    dftemp = Severity_Mapping(dftemp, dfSEVERITYCodelist)

    dfROUTECodelist = dfSENDCodelist[dfSENDCodelist["Codelist Name"].isin(["Route of Administration Response"])].drop_duplicates()
    dftemp = Route_Mapping(dftemp, dfROUTECodelist)

    dfDISPOSITIONCodelist = dfSENDCodelist[dfSENDCodelist["Codelist Name"].isin(["Standardized Disposition Term"])].drop_duplicates()
    dftemp = Disposition_Mapping(dftemp, dfDISPOSITIONCodelist)
    dftemp = dftemp.fillna("")

    dfSpecimenCodelist = dfSENDCodelist[dfSENDCodelist["Codelist Name"].isin(["Specimen"])].drop_duplicates()
    dfAnatLocCodelist = dfSENDCodelist[dfSENDCodelist["Codelist Name"].isin(["Anatomical Location"])].drop_duplicates()
    dfNonNeoplasmCodelist = dfSENDCodelist[dfSENDCodelist["Codelist Name"].isin(["Non-Neoplastic Finding Type"])].drop_duplicates()

    for index in dftemp.index:
        row = Specimen_Mapping(dftemp.iloc[index], dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist)
        dftemp.iloc[index] = row
        dftemp.iloc[index] = dftemp.iloc[index].fillna("")
        row = Specimen_Mapping(dftemp.iloc[index], dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist)
        val1 = dftemp.at[index, 'MISPEC']
        val2 = dftemp.at[index, 'MIANTREG']
        val3 = dftemp.at[index, 'MISPEC_raw']
        val4 = dftemp.at[index, 'MIANTREG_raw']
        val5 = dftemp.at[index, 'Anatomical_Region_of_Specimen']

        if(val1 and not val1.isspace()):
            pass
        else:
            if (val3 and not val3.isspace()):
                try:
                    dfSpecimenCodelisttemp = dfClean(dfSpecimenCodelist)

                    row['MISPEC'] = dfSpecimenCodelist.loc[dfSpecimenCodelisttemp["CDISC Submission Value"] == val3, "CDISC Submission Value"].item()
                    dftemp.iloc[index] = row
                    dftemp.iloc[index] = dftemp.iloc[index].fillna("")
                    row = Specimen_Mapping(dftemp.iloc[index], dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist)
                except:
                    row['MISPEC'] = "BIOSPECIMEN"
                    row['Sample_Type_Group'] = "OTHER"
                    row['CODE_CONFLICT_FLAG'] = "MISPEC: " + val3
                    dftemp.iloc[index] = row
                    dftemp.iloc[index] = dftemp.iloc[index].fillna("")
                    row = Specimen_Mapping(dftemp.iloc[index], dfSpecimenCodelist, dfAnatLocCodelist,
                                           dfNonNeoplasmCodelist)
            else:
                row['MISPEC'] = "BIOSPECIMEN"
                row['Sample_Type_Group'] = "OTHER"
                row['CODE_CONFLICT_FLAG'] = "MISPEC: NULL"
                dftemp.iloc[index] = row
                dftemp.iloc[index] = dftemp.iloc[index].fillna("")
                row = Specimen_Mapping(dftemp.iloc[index], dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist)

        if(val2 and not val2.isspace()):
            pass
        else:
            if (val4 and not val4.isspace()):
                try:
                    dfAnatLocCodelisttemp = dfClean(dfAnatLocCodelist)

                    row['MIANTREG'] = dfAnatLocCodelist.loc[dfAnatLocCodelisttemp["CDISC Submission Value"] == val4, "CDISC Submission Value"].item()
                    dftemp.iloc[index] = row
                    dftemp.iloc[index] = dftemp.iloc[index].fillna("")
                    row = Specimen_Mapping(dftemp.iloc[index], dfSpecimenCodelist, dfAnatLocCodelist, dfNonNeoplasmCodelist)
                except:
                    if (val5 and not val5.isspace()):
                        pass
                    else:
                        row['CODE_CONFLICT_FLAG'] = "MIANTREG: " + val4
                        dftemp.iloc[index] = row
                        dftemp.iloc[index] = dftemp.iloc[index].fillna("")


    # Creating Output File
    dftemp.to_csv("//lrlhps/users/c143390/Ontology/HistopathTerminology_test_5Aug2021.csv", index=False, header=True)

if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
