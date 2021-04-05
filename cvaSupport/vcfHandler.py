import vcf
import vcf.model

try:
    import fileHandling
except ImportError:
    from . import fileHandling
import gzip
import os
import io
import scipy.stats
import json


class BiasTable:

    def __init__(self, biasTable: list):
        if not len(biasTable) == 4:
            raise ValueError("Bias tables must be 4 items long")
        for item in biasTable:
            if not type(item) == int:
                raise ValueError("All items in bias table must be integers")
        self.forwardRefCount, self.reverseRefCount, self.forwardAltCount, self.reverseAltCount = biasTable
        self.oddsRatio, self.pValue = scipy.stats.fisher_exact(
            [[self.forwardRefCount, self.reverseRefCount], [self.forwardAltCount, self.reverseAltCount]])

    @property
    def totalRefReads(self):
        return self.forwardRefCount + self.reverseRefCount

    @property
    def totalAltReads(self):
        return self.forwardAltCount + self.reverseAltCount

    @property
    def toDict(self):
        return {
            "ref": {
                "fwd": self.forwardRefCount,
                "rev": self.reverseRefCount
            },
            "alt": {
                "fwd": self.forwardAltCount,
                "rev": self.reverseAltCount
            }
        }

    @classmethod
    def fromDict(cls, dictionary: dict):
        biasTableList = [
            dictionary["ref"]["fwd"],
            dictionary["ref"]["rev"],
            dictionary["alt"]["fwd"],
            dictionary["alt"]["rev"]
        ]
        return cls(biasTableList)

    def __str__(self):
        return "%s/%s:%s/%s %s" % (
        self.forwardRefCount, self.reverseRefCount, self.forwardAltCount, self.reverseAltCount, self.pValue)


class VariantRecord:

    def __init__(self, contig: str, position: int, ref: str, alt: str, filter: list, isDeletion: bool, isIndel: bool,
                 isSNV: bool, mateBiasTable: list = None, strandBiasTable: list = None, sampleID:str=""):
        self.contig = contig
        self.position = position
        self.ref = ref
        self.alt = alt
        self.filter = filter
        self.isDeletion = isDeletion
        self.isIndel = isIndel
        self.isSNV = isSNV
        if strandBiasTable:
            self.strandBias = BiasTable(strandBiasTable)
        else:
            self.strandBias = None
        if mateBiasTable:
            self.mateBias = BiasTable(mateBiasTable)
        else:
            self.mateBias = None
        self.sampleID = sampleID

    @property
    def refDepth(self):
        if self.strandBias:
            return self.strandBias.totalRefReads
        elif self.mateBias:
            return self.mateBias.totalRefReads
        else:
            return None

    @property
    def altDepth(self):
        if self.strandBias:
            return self.strandBias.totalAltReads
        elif self.mateBias:
            return self.mateBias.totalAltReads
        else:
            return None

    @classmethod
    def fromVCFRecord(cls, record: vcf.model._Record, altAlleleNumber: int = 0, sampleNumber: int = 0):
        try:
            mateBias = record.samples[sampleNumber].data.MB
        except AttributeError:
            mateBias = None
        try:
            strandBias = record.samples[sampleNumber].data.SB
        except AttributeError:
            strandBias = None
        return cls(record.CHROM, record.POS, record.REF, record.ALT[altAlleleNumber], record.FILTER, record.is_deletion,
                   record.is_indel, record.is_snp, strandBias,
                   mateBias, record.samples[sampleNumber].sample)

    def toDict(self):
        dictionary = {
            "contig": self.contig,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
            "filter": self.filter,
            "isDeletion": self.isDeletion,
            "isIndel": self.isIndel,
            "isSNV": self.isSNV,
            "strandBias": None,
            "mateBias": None
        }
        if not self.strandBias is None:
            dictionary["strandBias"] = self.strandBias.toDict()
        if not self.mateBias is None:
            dictionary["mateBias"] = self.mateBias.toDict()
        return dictionary

    def toJSON(self):
        return json.dumps(self.toDict())

    @property
    def vepIdentifier(self):
        if self.isSNV:
            return "%s_%s_%s/%s" % (self.contig, self.position, self.ref, self.alt)
        elif self.isIndel:
            if self.isDeletion:
                return "%s_%s_%s/%s" % (self.contig, self.position + 1, self.ref[1:], "-")
            else:
                return "%s_%s_%s/%s" % (self.contig, self.position + 1, "-", str(self.alt)[1:])
        else:
            raise RuntimeError("There should be no way to get to this value, please investigate")

    @property
    def standardMutationIdentifier(self):
        if self.isSNV:
            return "%s_%s_%s/%s" % (self.contig, self.position, self.ref, self.alt)
        elif self.isIndel:
            if self.isDeletion:
                return "%s_%s:del%s" % (self.contig, self.position + 1, self.ref[1:])
            else:
                return "%s_%s:ins%s" % (self.contig, self.position + 1, str(self.alt)[1:])
        else:
            raise RuntimeError("There should be no way to get to this value, please investigate")

    def __str__(self):
        return self.standardMutationIdentifier

def readVCF(vcfPath: str):
    if not os.path.isfile(vcfPath):
        raise FileNotFoundError("Unable to find VCF at %s" % vcfPath)
    if fileHandling.isGzipped(vcfPath):
        vcfFileHandle = gzip.open(vcfPath, 'rt')
    else:
        vcfFileHandle = open(vcfPath, 'r')
    vcfStream = io.StringIO()
    vcfStream.write(vcfFileHandle.read())
    vcfFileHandle.close()
    vcfStream.seek(0)
    vcfHandle = vcf.Reader(vcfStream)
    recordList = [record for record in vcfHandle]
    return vcfHandle, recordList


def processVCFRecordList(recordList: list):
    processedRecordCollection = []
    for record in recordList:
        processedRecordCollection.append(VariantRecord.fromVCFRecord(record))
    return processedRecordCollection


def processVCF(vcfPath:str):
    if not os.path.isfile(vcfPath):
        raise FileNotFoundError("Unable to find VCF at %s" %vcfPath)
    vcfHandle, recordList = readVCF(vcfPath)
    return processVCFRecordList(recordList)



if __name__ == "__main__":
    vcfHandle, recordList = readVCF("/opt/vep/ronavep/references/in2442-23.hard-filtered.vcf.gz")
    processedRecords = processVCFRecordList(recordList)
    print("something")