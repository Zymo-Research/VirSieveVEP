try:
    import fileHandling
except ImportError:
    from . import fileHandling
import gzip
import os


missenseVariants = [
    "missense_variant"
]

nonsenseVariants = [
    "frameshift_variant",
    "frameshift_variant,start_lost",
    "stop_gained"
]

inFrameIndelVariants = [
    "inframe_deletion",
    "inframe_insertion"
]

synonymousVariants = [
    "synonymous_variant"
]

outsideGeneVariants = [
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intergenic_variant",
    "-"
]

codingVariants = missenseVariants + nonsenseVariants + inFrameIndelVariants


class VariantEffect:

    def __init__(self, identifier:str, location:[str, int], allele:str, gene:str, feature:str, featureType:str, consequence:str, cDNAPosition:str, cdsPosition:str, proteinPosition:str, aminoAcid:str, codons:str, existingVariation:str, notes:str):
        self.identifier = identifier
        identifierSplit = identifier.split("_")
        if not len(identifierSplit) >= 3:
            raise ValueError("VEP identifiers should be at least three fields")
        variation = identifierSplit[-1]
        self.position = int(identifierSplit[-2])
        self.contig = "_".join(identifierSplit[:-2])
        self.ref, self.alt = variation.split("/")[:2] #rare cases with 3 items, but item 3 appears not annotated. Need to dig more on this.
        self.location = str(location)
        self.allele = allele
        self.gene = gene
        self.feature = feature
        self.featureType = featureType
        self.consequence = consequence
        self.cDNAPosition = cDNAPosition
        if self.cDNAPosition == "-":
            self.cDNAPosition = None
        self.cdsPosition = cdsPosition
        if self.cdsPosition == "-":
            self.cdsPosition = None
        self.proteinPosition = proteinPosition
        if self.proteinPosition == "-":
            self.proteinPosition = None
        self.aminoAcid = aminoAcid
        if self.aminoAcid == "-":
            self.aminoAcid = None
        self.codons = codons
        if self.codons == "-":
            self.codons = None
        self.existingVariation = existingVariation
        self.notes = notes
        self.isMissense = consequence in missenseVariants
        self.isNonsense = consequence in nonsenseVariants
        self.isCoding = consequence in codingVariants
        self.isMakingMutantProtein = consequence in codingVariants and consequence not in nonsenseVariants
        self.isSynonymous = consequence in synonymousVariants
        self.isInsideGene = consequence not in outsideGeneVariants

    @classmethod
    def fromVEPLine(cls, vepLine:str):
        vepSplit = vepLine.split("\t")
        if not len(vepSplit) == 14:
            raise ValueError("VEP output lines should have exactly 14 elements")
        vepSplit = [item.strip() for item in vepSplit]
        return cls(*vepSplit)

    @classmethod
    def variantOfNoEffect(cls, identifier:str, position:[str, int], alt:str):
        return cls(identifier, position, alt, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "")

    @property
    def proteinChangeNotation(self):
        if not self.aminoAcid:
            return ""
        if not self.isCoding:
            return ""
        reference, alt = self.aminoAcid.split("/")
        if reference == "-":
            return "%s:ins%s %s" %(self.gene, self.proteinPosition, alt)
        elif alt == "-":
            return "%s:del%s %s" %(self.gene, self.proteinPosition, reference)
        else:
            return "%s:%s%s%s" %(self.gene, reference, self.proteinPosition, alt)

    def __str__(self):
        proteinChange = self.proteinChangeNotation
        if proteinChange:
            proteinChange = " " + proteinChange
        return "%s%s" %(self.identifier, proteinChange)
    
    
def processVEPFile(vepFilePath:str):
    vepLineCollection = []
    if not os.path.isfile(vepFilePath):
        raise FileNotFoundError("Unable to find VCF at %s" % vepFilePath)
    if fileHandling.isGzipped(vepFilePath):
        vepFileHandle = gzip.open(vepFilePath, 'rt')
    else:
        vepFileHandle = open(vepFilePath, 'r')
    for line in vepFileHandle:
        if line.startswith("#"):
            continue
        if not line.strip():
            continue
        vepLineCollection.append(VariantEffect.fromVEPLine(line))
    vepFileHandle.close()
    return vepLineCollection


if __name__ == "__main__":
    variantEffects = processVEPFile("/opt/vep/ronavep/references/in2442-23.vep.txt")
    print("something")