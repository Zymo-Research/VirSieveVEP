import math
import typing
try:
    import vcfHandler
    import vepHandler
except ImportError:
    from . import vcfHandler
    from . import vepHandler


class CombinedMutantData:

    def __init__(self, vcf:vcfHandler.VariantRecord, vep:vepHandler.VariantEffect=None):
        self.variantCall = vcf
        self.locus = (vcf.contig, vcf.position)
        if vep is None:
            identifier = vcf.vepIdentifier
            location = vcf.position
            alt = vcf.alt
            vep = vepHandler.VariantEffect.variantOfNoEffect(identifier, location, alt)
        self.variantEffect = vep
        self.vepDNANotation = vep.identifier
        self.standardDNANotation = vcf.vepIdentifier
        self.consequence = vep.consequence
        self.proteinChange = vep.proteinChangeNotation
        self.isMissense = vep.isMissense
        self.isNonsense = vep.isNonsense
        self.isCoding = vep.isCoding
        self.isMakingMutantProtein = vep.isMakingMutantProtein
        self.isSynonymous = vep.isSynonymous
        self.isInsideGene = vep.isInsideGene
        self.refDepth = vcf.refDepth
        self.altDepth = vcf.altDepth
        self.totalDepth = self.refDepth + self.altDepth
        self.confidenceLevel = None
        try:
            self.percentAlt = self.altDepth / (self.refDepth + self.altDepth)
        except ZeroDivisionError:
            self.percentAlt = float("NaN")
        if not vcf.mateBias:
            self.mateBiasInvLog = None
        else:
            try:
                self.mateBiasInvLog = -math.log10(vcf.mateBias.pValue)
            except ValueError:
                self.mateBiasInvLog = None
        if not vcf.strandBias:
            self.strandBiasInvLog = None
        else:
            self.strandBiasInvLog = -math.log10(vcf.strandBias.pValue)
        self.flags = []
        if vcf.filter:
            for filterValue in vcf.filter:
                if not filterValue == "PASS":
                    self.flags.append(filterValue)
        self.alerts = []

    def addFlag(self, flag:str):
        if not flag in self.flags:
            self.flags.append(flag)

    def addAlert(self, alert:str):
        if not alert in self.alerts:
            self.alerts.append(alert)

    @property
    def sharingDict(self):
        if self.variantCall.strandBias is None:
            strandBiasDict = None
        else:
            strandBiasDict = self.variantCall.strandBias.toDict
        if self.variantCall.mateBias is None:
            mateBiasDict = None
        else:
            mateBiasDict = self.variantCall.mateBias.toDict
        dictionary = {
            "nucleic acid change": self.vepDNANotation,
            "consequence": self.consequence,
            "protein change": self.proteinChange,
            "total depth": self.totalDepth,
            "alt percent": self.percentAlt,
            "strand counts": strandBiasDict,
            "mate counts": mateBiasDict
        }
        return dictionary

    def __hash__(self):
        return hash(self.vepDNANotation)

    def __str__(self):
        proteinChange = self.proteinChange
        if not proteinChange:
            proteinChange = "N/A"
        if self.mateBiasInvLog is None:
            mateBias = "None"
        else:
            mateBias = round(self.mateBiasInvLog, 2)
        if self.strandBiasInvLog is None:
            strandBias = "None"
        else:
            strandBias = round(self.strandBiasInvLog, 2)
        outputStr = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(self.standardDNANotation, self.consequence, proteinChange, round(self.percentAlt * 100, 2), self.totalDepth, mateBias, strandBias , "|".join(self.flags), "|".join(self.alerts), self.confidenceLevel)
        return outputStr


def mergeVCFandVEPTables(vcfTable:typing.Dict[str, vcfHandler.VariantRecord], vepTable:typing.Dict[str, vepHandler.VariantEffect]):
    mergeList = []
    for vcfIdentifier, vcfRecord in vcfTable.items():
        vepRecord = vepTable.setdefault(vcfIdentifier, None)
        mergeList.append(CombinedMutantData(vcfRecord, vepRecord))
    return mergeList
