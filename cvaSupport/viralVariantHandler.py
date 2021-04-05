import json


class StrainIdentifier:

    def __init__(self, identifier: str, commonName: str, location: str, aliases: list, source: str, proteinChanges: list):
        self.identifier = identifier
        self.commonName = commonName
        self.location = location
        self.aliases = aliases
        self.source = source
        self.proteinChanges = self.makeMutationDictionary(proteinChanges)

    @staticmethod
    def makeMutationDictionary(mutations: list):
        if type(mutations) == dict:
            return mutations
        mutationDict = {}
        for mutation in mutations:
            gene, variant = mutation.split(":")
            gene = gene.strip()
            variant = mutation.strip()
            if not gene in mutationDict:
                mutationDict[gene] = []
            mutationDict[gene].append(variant)
        return mutationDict

    @property
    def dictionary(self):
        dictionary = {
            "identifier": self.identifier,
            "common name": self.commonName,
            "location": self.location,
            "aliases": self.aliases,
            "mutations": self.proteinChanges,
            "source": self.source
        }
        return dictionary

    @property
    def json(self):
        return json.dumps(self.dictionary, indent=4)

    @classmethod
    def fromDict(cls, dictionary: dict):
        return cls(
            dictionary["identifier"],
            dictionary["common name"],
            dictionary["location"],
            dictionary["aliases"],
            dictionary["mutations"],
            dictionary["source"]
        )

    @classmethod
    def fromJSON(cls, jsonString: str):
        dictionary = json.loads(jsonString)
        return cls.fromDict(dictionary)


class ProteinMutationIdentifier:

    def __init__(self, gene:str, ref:str, position:[str, int], alt:str, note:str="", source:str=""):
        self.gene = gene
        self.position = str(position)
        self.ref = ref
        if not ref:
            self.ref = "-"
        self.alt = alt
        if not alt:
            self.alt = "-"
        self.note = note
        self.source = source

    @property
    def proteinChangeNotation(self):
        if self.ref == "-":
            return "%s:ins%s %s" % (self.gene, self.position, self.alt)
        elif self.alt == "-":
            return "%s:del%s %s" % (self.gene, self.position, self.ref)
        else:
            return "%s:%s%s%s" % (self.gene, self.ref, self.position, self.alt)

    @property
    def dictionary(self):
        dictionary = {
            "gene": self.gene,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
            "note": self.note,
            "source": self.source
        }
        return dictionary

    @property
    def json(self):
        return json.dumps(self.dictionary, indent=4)

    @classmethod
    def fromDict(cls, dictionary:dict):
        return cls(**dictionary)

    @classmethod
    def fromJSON(cls, jsonString:str):
        return cls.fromDict(json.loads(jsonString))


class NucleicAcidMutationIdentifier:

    def __init__(self, contig:str, position:[str, int], ref:str, alt:str, note:str="", source:str=""):
        self.contig = contig
        self.position = str(position)
        self.ref = ref
        self.alt = alt
        self.note = note
        self.source = source

    @property
    def vepIdentifier(self):
        return "%s_%s_%s/%s" % (self.contig, self.position, self.ref, self.alt)

    @property
    def dictionary(self):
        dictionary = {
            "contig": self.contig,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
            "note": self.note,
            "source": self.source
        }
        return dictionary

    @property
    def json(self):
        return json.dumps(self.dictionary, indent=4)

    @classmethod
    def fromDict(cls, dictionary: dict):
        return cls(**dictionary)

    @classmethod
    def fromJSON(cls, jsonString: str):
        return cls.fromDict(json.loads(jsonString))
