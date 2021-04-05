import os
import gzip
import typing
try:
    import fileHandling
    import mutationDataMerge
except ImportError:
    from . import fileHandling
    from . import mutationDataMerge


problematicSitesFile = "/opt/vep/ronavep/references/problematicSiteFilter.vcf"


def getProblematicSitesTable(problematicSitesFilePath:str=problematicSitesFile):
    if not os.path.isfile(problematicSitesFilePath):
        raise FileNotFoundError("Unable to find problematic sites list at %s" %problematicSitesFilePath)
    if fileHandling.isGzipped(problematicSitesFilePath):
        vcfFileHandle = gzip.open(problematicSitesFilePath, 'rt')
    else:
        vcfFileHandle = open(problematicSitesFilePath, 'r')
    positionTable = {}
    for line in vcfFileHandle:
        if not line.strip():
            continue
        if line.startswith("#"):
            continue
        contig, position, identifier, ref, alt, qual, filterValue, info = line.split("\t")
        position = int(position)
        positionTable[position] = filterValue
    return positionTable


def applySiteWarnings(mutantDataList:typing.List[mutationDataMerge.CombinedMutantData], biasWarningInvLog:float=2.0, problemSiteTable:dict=None):
    if not problemSiteTable:
        problemSiteTable = getProblematicSitesTable()
    for variantRecord in mutantDataList:
        position = variantRecord.locus[1]
        if position in problemSiteTable:
            flag = "ProblemSite:%s" %problemSiteTable[position]
            variantRecord.addFlag(flag)
        if variantRecord.strandBiasInvLog and variantRecord.strandBiasInvLog >= biasWarningInvLog:
            variantRecord.addFlag("StrandBias")
        if variantRecord.mateBiasInvLog and variantRecord.mateBiasInvLog >= biasWarningInvLog:
            variantRecord.addFlag("MateBias")