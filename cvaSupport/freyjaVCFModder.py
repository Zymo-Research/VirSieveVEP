FORMATLINESTOADD = [
    '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Alternate allele depth">',
    '##FORMAT=<ID=ALT_FREQ,Number=1,Type=Float,Description="Alternate allele frequency">'
]

def processVCFLine(vcfLine:str) -> str:
    vcfLine = vcfLine.strip()
    if vcfLine.startswith("#"):
        return vcfLine
    else:
        return modifyVCFLine(vcfLine)


def modifyVCFLine(vcfLine:str) -> str:
    lineElements = vcfLine.split("\t")
    formatFields = lineElements[-2]
    dataFields = lineElements[-1]
    formatFieldList = formatFields.split(":")
    dataFieldList = dataFields.split(":")
    fieldTable = {}
    for field, data in zip(formatFieldList, dataFieldList):
        fieldTable[field] = data
    allelicDepths = fieldTable["AD"].split(",")
    allelicFrequencies = fieldTable["AF"].split(",")
    variantAllelicDepths = allelicDepths[1:]
    newFormatFieldList = formatFieldList.copy()
    newFormatFieldList.append("ALT_DP")
    newFormatFieldList.append("ALT_FREQ")
    fieldTable["ALT_DP"] = ",".join(variantAllelicDepths)
    fieldTable["ALT_FREQ"] = ",".join(allelicFrequencies)
    newDataFieldsList = []
    for element in newFormatFieldList:
        newDataFieldsList.append(fieldTable[element])
    newDataFields = ":".join(newDataFieldsList)
    newFormatFields = ":".join(newFormatFieldList)
    newLineElements = lineElements.copy()
    newLineElements[-2] = newFormatFields
    newLineElements[-1] = newDataFields
    newLine = "\t".join(newLineElements)
    return newLine


def processVCF(vcfPathInputPath:str, vcfOutputPath:str):
    addedFormatLines = False
    vcf = open(vcfPathInputPath, 'r')
    vcfOutput = open(vcfOutputPath, 'w')
    for line in vcf:
        newLine = processVCFLine(line)
        if newLine.startswith("#CHROM"):
            if FORMATLINESTOADD and not addedFormatLines:
                for line in FORMATLINESTOADD:
                    print(line, file=vcfOutput)
                addedFormatLines = True
        print(newLine, file=vcfOutput)
    vcf.close()
    vcfOutput.close()


if __name__ == "__main__":
    processVCF("C:\\Users\\mweinstein\\miWastewater\\alignmentArtifactFilteredVCF\\Cov-NEB-210503.orientationFilter.alignmentArtifactFilter.vcf", "C:\\Users\\mweinstein\\miWastewater\\alignmentArtifactFilteredVCF\\Cov-NEB-210503.orientationFilter.alignmentArtifactFilter.freyjaMod.vcf")