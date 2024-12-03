import os
import cvaSupport
import operator
import typing
import json
import scipy.stats
import math
import re

workingFolderEnv = os.environ.setdefault("WORKINGFOLDER", "/data")
if not os.path.isdir(workingFolderEnv):
    raise NotADirectoryError("Unable to find working directory at %s" %workingFolderEnv)
inputFolderEnv = os.environ.setdefault("INPUTFOLDER", os.path.join(workingFolderEnv, "filteredVCF"))
if not os.path.isdir(inputFolderEnv):
    raise NotADirectoryError("Unable to find input folder at %s" %inputFolderEnv)
# stringentFilteredVCFFolderEnv = os.environ.setdefault("STRINGENTVCFFOLDER", os.path.join(workingFolderEnv, "alignmentArtifactFilteredVCF"))
# if not os.path.isdir(stringentFilteredVCFFolderEnv):
#     print("WARNING: Stringent filtered VCF folder was not found at %s. High confidence variants will not be identified." %stringentFilteredVCFFolderEnv)
vepIntermediatesFolderEnv = os.environ.setdefault("VEPINTERMEDIATESFOLDER", os.path.join(workingFolderEnv, "vepOutputs"))
if not os.path.isdir(vepIntermediatesFolderEnv):
    os.mkdir(vepIntermediatesFolderEnv)
resultsFolderEnv = os.environ.setdefault("RESULTSFOLDER", os.path.join(workingFolderEnv, "results"))
if not os.path.isdir(resultsFolderEnv):
    os.mkdir(resultsFolderEnv)
freyjaOutputFolderEnv = os.environ.setdefault("FREYJAOUTPUTFOLDER", os.path.join(workingFolderEnv, "freyjaOutput"))
if not os.path.isdir(freyjaOutputFolderEnv):
    os.mkdir(freyjaOutputFolderEnv)


def getVCFList(folder:str=inputFolderEnv):
    folderContents = os.listdir(folder)
    folderFilesRaw = [os.path.join(folder, item) for item in folderContents]
    folderFilesFiltered = []
    for item in folderFilesRaw:
        if not os.path.isfile(item):
            continue
        if item.endswith(".vcf") or item.endswith(".vcf.gz"):
            folderFilesFiltered.append(item)
    return folderFilesFiltered


def runVEP(vcfPath:str):
    if not os.path.isfile(vcfPath):
        raise FileNotFoundError("Unable to find a VCF at %s" %vcfPath)
    vepOutput = cvaSupport.vepRunner.runVEP(vcfPath, vepIntermediatesFolderEnv)
    if vepOutput:
        return vepOutput
    else:
        raise RuntimeError("VEP appears to have had a failed run on sample %s" %vcfPath)


def makeVCFJoiningTable(vcfPath:str, returnSampleID:bool=False):
    if not os.path.isfile(vcfPath):
        raise FileNotFoundError("Unable to find a VCF at %s" %vcfPath)
    processedRecords = cvaSupport.vcfHandler.processVCF(vcfPath)
    sampleID = processedRecords[0].sampleID
    joiningTable = {}
    for mutation in processedRecords:
        joiningTable[mutation.vepIdentifier] = mutation
    if returnSampleID:
        return joiningTable, sampleID
    return joiningTable


def makeVEPJoiningTable(vepOutputPath:str):
    if not os.path.isfile(vepOutputPath):
        raise ValueError("Unable to find VEP output file at %s" %vepOutputPath)
    mutationEffects = cvaSupport.vepHandler.processVEPFile(vepOutputPath)
    joiningTable = {}
    for mutation in mutationEffects:
        joiningTable[mutation.identifier] = mutation
    return joiningTable


# def findStringentFilteredVCF(originalVCFPath:str, stringentFilteredVCFFolder:str=stringentFilteredVCFFolderEnv):
#     def baseName(filePath:str):
#         return os.path.split(filePath)[1].split(".")[0]
#     if not os.path.isdir(stringentFilteredVCFFolder):
#         return None
#     vcfList = getVCFList(stringentFilteredVCFFolder)
#     originalVCFBaseName = baseName(originalVCFPath)
#     candidates = []
#     for vcf in vcfList:
#         if baseName(vcf) == originalVCFBaseName and vcf.endswith("alignmentArtifactFilter.vcf"):
#             candidates.append(vcf)
#     if len(candidates) == 1:
#         return candidates[0]
#     elif not candidates:
#         return None
#     else:
#         raise RuntimeError("Found multiple potential candidates for stringent filter for original VCF %s. Candidates were %s" %(originalVCFPath, candidates))


def calculateAndApplyConfidenceScoreToMergedMutation(mergedMutation:cvaSupport.mutationDataMerge.CombinedMutantData, highConfidenceMutationList:typing.List[str]):
    confidenceLossFlags = ["weak_evidence", "orientation", "mate_bias", "strand_bias", "ProblemSite:mask", "ProblemSite:caution"]
    confidence = 2
    if mergedMutation.standardDNANotation in highConfidenceMutationList:
        if "ProblemSite:mask" in mergedMutation.flags:
            confidence = 3
        elif "ProblemSite:caution" in mergedMutation.flags:
            confidence = 2
        else:
            for badFlag in confidenceLossFlags:
                badFlagged = False
                if badFlag in mergedMutation.flags:
                    confidence = 3
                    badFlagged = True
                if not badFlagged:
                    confidence = 1
    else:
        if mergedMutation.variantCall.isIndel:
            if "clustered_events" in mergedMutation.flags:
                confidence = 3
        for badFlag in confidenceLossFlags:
            if badFlag in mergedMutation.flags:
                confidence = 3
    mergedMutation.confidenceLevel = confidence
    return mergedMutation


# def applyConfidenceScoresToMergedMutationList(vcfPath:str, mergedMutationList:typing.List[cvaSupport.mutationDataMerge.CombinedMutantData], stringentFilteredVCFFolder:str=stringentFilteredVCFFolderEnv):
#     stringentFilterVCF = findStringentFilteredVCF(vcfPath, stringentFilteredVCFFolder)
#     if not stringentFilterVCF:
#         stringentVCFTable = {}
#         print("No stringent filter VCF found to match %s. No high-confidence variants will be identified." %vcfPath)
#     else:
#         stringentVCFTable = makeVCFJoiningTable(stringentFilterVCF)
#         freyjaModVCF = makeFreyjaVCFMods(stringentFilterVCF)
#     highConfidenceList = list(stringentVCFTable.keys())
#     for mergedMutation in mergedMutationList:
#         calculateAndApplyConfidenceScoreToMergedMutation(mergedMutation, highConfidenceList)
#     return mergedMutationList


def makeFreyjaVCFMods(vcfPath:str, outputFolder:str=freyjaOutputFolderEnv):
    outputVCFPath = os.path.join(outputFolder, os.path.split(vcfPath[:-3])[1] + "freyjaMod.vcf")
    cvaSupport.freyjaVCFModder.processVCF(vcfPath, outputVCFPath)
    return outputVCFPath


def makeResultsTables():
    vcfList = getVCFList()
    results = {}
    strainObservationsTable = {}
    for vcfPath in vcfList:
        vepOutput = runVEP(vcfPath)
        freyjaModVCF = makeFreyjaVCFMods(vcfPath)
        vcfTable, sampleID = makeVCFJoiningTable(vcfPath, returnSampleID=True)
        print("Analyzing %s" %sampleID)
        vepTable = makeVEPJoiningTable(vepOutput)
        mergedResults = cvaSupport.mutationDataMerge.mergeVCFandVEPTables(vcfTable, vepTable)
        results[sampleID] = mergedResults
        results[sampleID].sort(key=operator.attrgetter("locus"))
        cvaSupport.problematicSites.applySiteWarnings(results[sampleID])
        strainObservations = cvaSupport.variantsOfConcernHandler.applyVariantsOfConcern(results[sampleID])
        #applyConfidenceScoresToMergedMutationList(vcfPath, results[sampleID])
        strainObservationsTable[sampleID] = strainObservations
        print("%s analysis completed." %sampleID)
    return results, strainObservationsTable


def writeBetaTable(mutationRecordTable:typing.Dict[str, typing.List[cvaSupport.mutationDataMerge.CombinedMutantData]]):
    outputFileName = "betaTable.txt"
    outputFilePath = os.path.join(resultsFolderEnv, outputFileName)
    outputFile = open(outputFilePath, 'w')
    columns = "\t".join(["#Sample", "Alpha", "Beta", "LowerLimit", "Scale", "PartialVariants", "PartialVariantsUnflagged"])
    print(columns, file=outputFile)
    for sampleID, variantRecords in mutationRecordTable.items():
        percentAlts = []
        partialVariants = 0
        partialVariantsUnflagged = 0
        for variantRecord in variantRecords:
            variantPercent = variantRecord.percentAlt
            if math.isnan(variantPercent):
                continue
            if variantPercent < 0:
                variantPercent = 0
            if 0.10 < variantPercent < 0.90:
                partialVariants += 1
                if not variantRecord.flags:
                    partialVariantsUnflagged += 1
            percentAlts.append(variantPercent)
        alpha, beta, lowerLimit, scale = scipy.stats.beta.fit(percentAlts)
        outputList = [sampleID, alpha, beta, lowerLimit, scale, partialVariants, partialVariantsUnflagged]
        outputList = [str(item) for item in outputList]
        outputString = "\t".join(outputList)
        print(outputString, file=outputFile)
    outputFile.close()


def writeVariantTables(mutationRecordTable:typing.Dict[str, typing.List[cvaSupport.mutationDataMerge.CombinedMutantData]]):
    columns = "\t".join(["#Sequence", "Consequence", "Protein", "PercentPrevalence", "ReadDepth", "InvLogMateBias", "InvLogStrandBias", "Flag/Filter", "Alert", "Confidence"])
    for sampleID, variantRecords in mutationRecordTable.items():
        sanitizedSampleID = re.sub("\W", "_", sampleID)
        outputFileName = "%s.variants.txt" %sanitizedSampleID
        outputFilePath = os.path.join(resultsFolderEnv, outputFileName)
        outputFile = open(outputFilePath, 'w')
        print(columns, file=outputFile)
        for variantRecord in variantRecords:
            print(variantRecord, file=outputFile)
        outputFile.close()


def writeStrainObservationTables(strainObservationTables:typing.Dict[str, typing.Dict[str, typing.Dict[str, float]]]):
    for sampleID, strainObservation in strainObservationTables.items():
        sanitizedSampleID = re.sub("\W", "_", sampleID)
        outputFileName = "%s.strainObservations.json" %sanitizedSampleID
        outputFilePath = os.path.join(resultsFolderEnv, outputFileName)
        outputFile = open(outputFilePath, 'w')
        json.dump(strainObservation, outputFile, indent=4)
        outputFile.close()


def writeVarShare(mutationRecordTable:typing.Dict[str, typing.List[cvaSupport.mutationDataMerge.CombinedMutantData]]):
    for sampleID, variantRecords in mutationRecordTable.items():
        sanitizedSampleID = re.sub("\W", "_", sampleID)
        outputFileName = "%s.variants.json" %sanitizedSampleID
        outputFilePath = os.path.join(resultsFolderEnv, outputFileName)
        finalOutput = {sampleID: {}}
        variantDict = finalOutput[sampleID]
        for variantRecord in variantRecords:
            variantID = variantRecord.vepDNANotation
            variantDict[variantID] = variantRecord.sharingDict
        outputFile = open(outputFilePath, 'w')
        json.dump(finalOutput, outputFile, indent=4)
        outputFile.close()


def main():
    results, strainObservationsTable = makeResultsTables()
    writeBetaTable(results)
    writeVariantTables(results)
    writeStrainObservationTables(strainObservationsTable)
    writeVarShare(results)
    print("DONE")


if __name__ == "__main__":
    main()
