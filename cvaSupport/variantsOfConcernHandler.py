import json
import gzip
import typing
try:
    import fileHandling
    import mutationDataMerge
except ImportError:
    from . import fileHandling
    from . import mutationDataMerge

variantsOfConcernFile = "/opt/vep/ronavep/references/variantsOfConcern.json"


def loadVariantsOfConcern(filePath:str = variantsOfConcernFile):
    if fileHandling.isGzipped(filePath):
        jsonFile = gzip.open(filePath, 'rt')
    else:
        jsonFile = open(filePath, 'r')
    data = json.load(jsonFile)
    jsonFile.close()
    return data


def applyVariantsOfConcern(mutantDataList:typing.List[mutationDataMerge.CombinedMutantData], variantsTable:dict=None):
    if not variantsTable:
        variantsTable = loadVariantsOfConcern()
    strainTable = variantsTable["strains"]
    molecularTable = variantsTable["molecular"]
    strainObservationTable = {}
    variantToStrainTable = {}
    allVariantsOfConcern = set()
    for strain in strainTable:
        strainObservationTable[strain] = {}
        for gene, variantList in strainTable[strain]["mutations"].items():
            for variant in variantList:
                if not variant in variantToStrainTable:
                    variantToStrainTable[variant] = []
                variantToStrainTable[variant].append(strain)
                strainObservationTable[strain][variant] = 0
                allVariantsOfConcern.add(variant)
    for variant in molecularTable:
        allVariantsOfConcern.add(variant)
    for variantRecord in mutantDataList:
        proteinChange = variantRecord.proteinChange
        if not proteinChange:
            continue
        if not proteinChange in allVariantsOfConcern:
            continue
        if proteinChange in variantToStrainTable:
            concern = "Strain associated: %s" %(", ".join(variantToStrainTable[proteinChange]))
            variantRecord.addAlert(concern)
        if proteinChange in molecularTable:
            variantRecord.addAlert(molecularTable[proteinChange]['note'])
        for strain, mutations in strainObservationTable.items():
            if proteinChange in mutations:
                strainObservationTable[strain][proteinChange] = variantRecord.percentAlt * 100
    return strainObservationTable


