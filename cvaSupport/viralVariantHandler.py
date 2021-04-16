import json
import typing


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

    @classmethod
    def fromString(cls, string, note:str="", source:str=""):
        string = string.strip()
        try:
            gene, change = string.split(":")
            change = change.strip()
        except ValueError:
            raise ValueError("Protein change notation should be Gene:change with only one colon character. That was not the number detected in %s." %string)
        if "ins" in change:
            ref = "-"
            change = change.replace("ins", "")
            try:
                pos, alt = change.split()
            except ValueError:
                raise ValueError("Protein insertion notation should be gene:ins[position] [sequence].  That notation was not parsed in %s" %string)
        elif "del" in change:
            alt = "-"
            change = change.replace("del", "")
            try:
                pos, ref = change.split()
            except ValueError:
                raise ValueError("Protein deletion notation should be gene:ins[position] [sequence].  That notation was not parsed in %s" %string)
        else:
            ref = change[0]
            alt = change[-1]
            pos = change[1:-1]
            try:
                pos = int(pos)
            except ValueError:
                raise ValueError("Protein substitution notation should be gene:[ref][position][alt] (such as S:D614G), but that format was not parsed from %s" %string)
        return cls(gene, ref, pos, alt, note, source)

    def __str__(self):
        return "%s_%s_%s_%s" %(self.gene, self.position, self.ref, self.alt)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.gene == other.gene and self.position == other.position and self.ref == other.ref and self.alt == other.alt

    def __gt__(self, other):
        if self == other:
            return False
        if self.gene > other.gene:
            return True
        elif self.gene < other.gene:
            return False
        else:
            if self.position > other.position:
                return True
            elif self.position < other.position:
                return False
            else:
                if self.ref > other.ref:
                    return True
                elif self.ref < other.ref:
                    return False
                else:
                    if self.alt > other.alt:
                        return True
                    elif self.alt < other.alt:
                        return False
                    else:
                        raise ValueError("Should never be able to reach this point.  Did so comparing greater than for %s and %s" %(self, other))

    def __lt__(self, other):
        return not (self > other or self == other)

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return self > other or self == other



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

    @classmethod
    def fromVEPString(cls, vepString:str, note:str="", source:str=""):
        identifierSplit = vepString.split("_")
        if not len(identifierSplit) >= 3:
            raise ValueError("VEP identifiers should be at least three fields")
        variation = identifierSplit[-1]
        position = int(identifierSplit[-2])
        contig = "_".join(identifierSplit[:-2])
        ref, alt = variation.split("/")[:2]
        return cls(contig, position, ref, alt, note, source)

    def __str__(self):
        return "%s_%s_%s_%s" % (self.contig, self.position, self.ref, self.alt)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.contig == other.contig and self.position == other.position and self.ref == other.ref and self.alt == other.alt

    def __gt__(self, other):
        if self == other:
            return False
        if self.contig > other.gene:
            return True
        elif self.contig < other.gene:
            return False
        else:
            if self.position > other.position:
                return True
            elif self.position < other.position:
                return False
            else:
                if self.ref > other.ref:
                    return True
                elif self.ref < other.ref:
                    return False
                else:
                    if self.alt > other.alt:
                        return True
                    elif self.alt < other.alt:
                        return False
                    else:
                        raise ValueError(
                            "Should never be able to reach this point.  Did so comparing greater than for %s and %s" % (
                            self, other))

    def __lt__(self, other):
        return not (self > other or self == other)

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return self > other or self == other


class MutationIdentifier:

    def __init__(self, proteinChange:[str, ProteinMutationIdentifier]=None, nucleicAcidChange:[str, NucleicAcidMutationIdentifier]=None, notes:str="", source:str=""):
        if not proteinChange and not nucleicAcidChange:
            raise ValueError("Unable to create a mutation data set with no mutation data associated.")
        if type(proteinChange) == str:
            proteinChange = ProteinMutationIdentifier.fromString(proteinChange)
        self.proteinChange = proteinChange
        if type(nucleicAcidChange) == str:
            nucleicAcidChange = NucleicAcidMutationIdentifier.fromVEPString(nucleicAcidChange)
        self.nucleicAcidChange = nucleicAcidChange
        self.notes = notes
        self.source = source

    @property
    def dictionary(self):
        proteinDict = None
        if self.proteinChange:
            proteinDict = self.proteinChange.dictionary
        nucleicAcidDict = None
        if self.nucleicAcidChange:
            nucleicAcidDict = self.nucleicAcidChange.dictionary
        returnDict = {
            "protein": proteinDict,
            "nucleic acid": nucleicAcidDict,
            "notes": self.notes,
            "source": self.source
        }
        return returnDict

    @classmethod
    def fromDict(cls, dictionary:dict):
        if not dictionary["protein"]:
            proteinChange = None
        else:
            proteinChange = ProteinMutationIdentifier.fromDict(dictionary["protein"])
        if not dictionary["nucleic acid"]:
            nucleicAcidChange = None
        else:
            nucleicAcidChange = NucleicAcidMutationIdentifier.fromDict(dictionary["nucleic acid"])
        source = dictionary["source"]
        notes = dictionary["notes"]
        return cls(proteinChange, nucleicAcidChange, notes, source)

    def __eq__(self, other):
        if self.nucleicAcidChange and other.nucleicAcidChange:
            return self.nucleicAcidChange == other.nucleicAcidChange
        elif self.proteinChange and other.proteinChange:
            return self.proteinChange == other.proteinChange


class StrainIdentifier:

    def __init__(self, identifier: str, commonName: str, location: str, aliases: list, mutations:typing.List[MutationIdentifier], source:str="", notes:str=""):
        self.identifier = identifier
        self.commonName = commonName
        self.location = location
        self.aliases = aliases
        self.source = source
        self.notes = notes
        self.mutationList = mutations
        self.nucleotideMutations, self.proteinMutations = self.makeMutationDictionaries(mutations)

    @staticmethod
    def makeMutationDictionaries(mutations: typing.List[MutationIdentifier]):
        nucleotideMutations = {}
        proteinMutations = {}
        for mutation in mutations:
            if mutation.proteinChange:
                gene = mutation.proteinChange.gene
                if not gene in proteinMutations:
                    proteinMutations[gene] = []
                proteinMutations[gene].append(mutation)
            if mutation.nucleicAcidChange:
                contig = mutation.nucleicAcidChange.contig
                if not contig in nucleotideMutations:
                    nucleotideMutations[contig] = {}
                position = mutation.nucleicAcidChange.position
                if not position in nucleotideMutations[contig]:
                    nucleotideMutations[contig][position] = []
                nucleotideMutations[contig][position].append(mutation)
        return nucleotideMutations, proteinMutations

    @property
    def dictionary(self):
        mutationList = [mutation.dictionary for mutation in self.mutationList]
        dictionary = {
            "identifier": self.identifier,
            "common name": self.commonName,
            "location": self.location,
            "aliases": self.aliases,
            "mutations": mutationList,
            "source": self.source,
            "notes": self.notes
        }
        return dictionary

    @property
    def json(self):
        return json.dumps(self.dictionary, indent=4)

    @classmethod
    def fromDict(cls, dictionary: dict):
        mutationList = []
        for mutationDict in dictionary["mutations"]:
            mutationList.append(MutationIdentifier.fromDict(mutationDict))
        return cls(
            dictionary["identifier"],
            dictionary["common name"],
            dictionary["location"],
            dictionary["aliases"],
            mutationList,
            dictionary["source"],
            dictionary["notes"]
        )

    @classmethod
    def fromJSON(cls, jsonString: str):
        dictionary = json.loads(jsonString)
        return cls.fromDict(dictionary)
