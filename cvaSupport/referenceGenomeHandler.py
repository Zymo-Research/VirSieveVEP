import os
from Bio import SeqIO, Seq
import io

referenceGenomePath = os.path.join(os.path.split(__file__)[0], "references", "MN908947.gb")


class ReferenceGenome:

    def __init__(self, gbkPath:str=referenceGenomePath):
        self.gbkPath = gbkPath
        if not os.path.isfile(gbkPath):
            raise FileNotFoundError("Unable to find reference file %s" %gbkPath)
        file = open(gbkPath, 'r')
        self.__gbkStream = io.StringIO()
        self.__gbkStream.write(file.read())
        file.close()
        self.__gbkStream.seek(0)
        self.__gbHandle = SeqIO.parse(self.__gbkStream, 'genbank')
        recordList = [record for record in self.__gbHandle]
        if not recordList:
            raise RuntimeError("Reference genome %s returned no records" % self.gbkPath)
        self.gbRecord = recordList[0]
        if len(recordList) > 1:
            print("WARNING: %s appears to have multiple records" % gbkPath)
        self.makeAttributes()

    def makeAttributes(self):
        self.geneTable = self.getGeneTable()
        self.geneTableByLocation = self.getGeneTableByLocation()
        self.sequenceTable = self.getSequenceTableByGene()
        self.translations = self.getTranslationTableByGene()

    def getGeneTable(self):
        def getGeneList():
            genes = {}
            for feature in self.gbRecord.features:
                if feature.type == "gene":
                    genes[feature.qualifiers["gene"][0]] = feature
            return genes
        geneList = getGeneList()
        return geneList

    def getGeneTableByLocation(self):
        geneTableByLocation = {}
        for geneID, geneData in self.geneTable.items():
            geneStart = geneData.location.nofuzzy_start
            geneEnd = geneData.location.nofuzzy_end
            geneLocationTuple = (geneStart, geneEnd)
            geneTableByLocation[geneLocationTuple] = geneID
        return geneTableByLocation

    def getSequenceTableByGene(self):
        sequenceTable = {}
        for coordinates, geneID in self.geneTableByLocation.items():
            start, end = coordinates
            sequence = self.gbRecord.seq[start:end]
            sequenceTable[geneID] = str(sequence)
        return sequenceTable

    def getTranslationTableByGene(self):
        translationTable = {}
        for geneID, sequence in self.sequenceTable.items():
            translationTable[geneID] = Seq.translate(sequence)
        return translationTable

if __name__ == "__main__":
    reference = ReferenceGenome()
    print("something")


