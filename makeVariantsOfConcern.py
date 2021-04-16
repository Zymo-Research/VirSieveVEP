import os
from cvaSupport.viralVariantHandler import StrainIdentifier, MutationIdentifier, ProteinMutationIdentifier, NucleicAcidMutationIdentifier
import json

contig = "NC_045512.2"

mutationList_B_1_351 = [
    MutationIdentifier(None, NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_174_G/T')),
    MutationIdentifier(None, NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_241_C/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf1ab:T265I'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_1059_C/T')),
    MutationIdentifier(None, NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_2692_A/T')),
    MutationIdentifier(None, NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_3037_C/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf1ab:K1655N'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_5230_G/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf1ab:H2799Y'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_8660_C/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf1ab:S2900L'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_8964_C/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf1ab:K3353R'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_10323_A/G')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf1ab:D4527Y'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_13843_G/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf1ab:P314L'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_14408_C/T')),
    MutationIdentifier(None, NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_14925_C/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf1ab:T5912I'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_17999_C/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('S:L18F'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_21614_C/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('S:D80A'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_21801_A/C')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('S:D215G'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_22206_A/G')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('S:R246I'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_22299_G/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('S:K417N'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_22813_G/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('S:E484K'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_23012_G/A')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('S:N501Y'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_23063_A/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('S:D614G'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_23403_A/G')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('S:A701V'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_23664_C or G?/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf3a:Q57H'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_25563_G/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('orf3a:S171L'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_25904_C/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('E:P71L'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_26456_C/T')),
    MutationIdentifier(None, NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_26645_C/T')),
    MutationIdentifier(None, NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_27504_T/C')),
    MutationIdentifier(None, NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_28253_C/T')),
    MutationIdentifier(ProteinMutationIdentifier.fromString('N:T205I'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_28887_C/T'))
]

mutationList_B_1_1_7 = [
MutationIdentifier(ProteinMutationIdentifier.fromString('Orf1ab:T1001I'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_3267_C/T')),
MutationIdentifier(ProteinMutationIdentifier.fromString('Orf1ab:A1708D'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_5388_C/alt')),
MutationIdentifier(ProteinMutationIdentifier.fromString('Orf1ab:I2230T'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_6954_T/C')),
MutationIdentifier(ProteinMutationIdentifier.fromString('Orf1ab:SGF 3675-3677 deletion'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_11288_/del')),
MutationIdentifier(ProteinMutationIdentifier.fromString('S:HV 69-70 deletion'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_21765_/del')),
MutationIdentifier(ProteinMutationIdentifier.fromString('S:Y144 deletion'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_21991_/del')),
MutationIdentifier(ProteinMutationIdentifier.fromString('S:N501Y'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_23063_A/T')),
MutationIdentifier(ProteinMutationIdentifier.fromString('S:A570D'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_23271_C/A')),
MutationIdentifier(ProteinMutationIdentifier.fromString('S:P681H'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_23604_C/A')),
MutationIdentifier(ProteinMutationIdentifier.fromString('S:T716I'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_23709_C/T')),
MutationIdentifier(ProteinMutationIdentifier.fromString('S:S982A'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_24506_T/G')),
MutationIdentifier(ProteinMutationIdentifier.fromString('S:D1118H'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_24914_G/C')),
MutationIdentifier(ProteinMutationIdentifier.fromString('Orf8:Q27stop'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_27972_C/T')),
MutationIdentifier(ProteinMutationIdentifier.fromString('Orf8:R52I'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_28048_G/T')),
MutationIdentifier(ProteinMutationIdentifier.fromString('Orf8:Y73C'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_28111_A/G')),
MutationIdentifier(ProteinMutationIdentifier.fromString('N:D3L'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_28280_G/C')),
MutationIdentifier(ProteinMutationIdentifier.fromString('N:D3L'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_28281_A/T')),
MutationIdentifier(ProteinMutationIdentifier.fromString('N:D3L'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_28282_T/A')),
MutationIdentifier(ProteinMutationIdentifier.fromString('N:S235F'), NucleicAcidMutationIdentifier.fromVEPString('NC_045512.2_28977_C/T'))

]

outputFolder = os.path.join(os.path.split(__file__)[0], "references")

b_1_1_7 = StrainIdentifier(
    identifier = "B.1.1.7",
    commonName = "UK",
    location = "UK",
    aliases = ["UK", "20I/501Y.V1"],
    source = "https://outbreak.info",
    proteinChanges= "ORF1ab:T1001I, ORF1ab:A1708D, ORF1ab:I2230T, ORF1ab:del3675-3677 SGF, S:del69-70 HV, S:del144 Y, S:N501Y, S:A570D, S:D614G, S:P681H, S:T716I, S:S982A, S:D1118H, ORF8:Q27*, ORF8:R52I, ORF8:Y73C, N:D3L, N:S235F".split(", ")
)

b_1_351 = StrainIdentifier(
    identifier = "B.1.351",
    commonName = "South Africa",
    location = "South Africa",
    aliases = ["South Africa", "20H/501Y.V2"],
    source="https://outbreak.info",
    proteinChanges= "ORF1ab:K1655N, E:P71L, N:T205I, S:K417N, S:E484K, S:N501Y, S:D614G, S:A701V".split(", ")
)

p_1 = StrainIdentifier(
    identifier = "P.1",
    commonName = "Brazil",
    location = "Brazil",
    aliases = ["Brazil", "Japan", "20J/501Y.V3"],
    source = "https://outbreak.info",
    proteinChanges= "N:P80R, ORF8:E92K, ORF3a:G174C, S:L18F, S:T20N, S:P26S, S:D138Y, S:R190S, S:K417T, S:E484K, S:N501Y, S:D614G, S:H655Y, S:T1027I, ORF1ab:E5662D, ORF1ab:K1795Q, ORF1ab:S1188L, ORF1ab:I760T, ORF1ab:del3675-3677 SGF, ORF1ab:F681L".split(", ")
)

b_1_427 = StrainIdentifier(
    identifier = "B.1.427",
    commonName = "West Coast",
    location = "California",
    aliases = ["CA VUI1", "CAL.20C"],
    source = "https://outbreak.info",
    proteinChanges= "S:D614G, ORF1ab:T265I, N:T205I, S:L452R, ORF1ab:S3158T, ORF3a:Q57H".split(", ") #ORF1b:D1183Y, ORF1b:P976L
)

b_1_429 = StrainIdentifier(
    identifier = "B.1.429",
    commonName="West Coast",
    location="California",
    aliases=["CA VUI1", "CAL.20C"],
    source = "https://outbreak.info",
    proteinChanges= "ORF1ab:T265I, ORF1ab:I4205V, S:S13I, S:W152C, S:L452R, S:D614G, ORF3a:Q57H, N:T205I".split(", ") #ORF1b:P314L, ORF1b:D1183Y
)

b_1_526 = StrainIdentifier(
    identifier = "B.1.526",
    commonName="New York",
    location="New York",
    aliases=[],
    source = "https://outbreak.info",
    proteinChanges= "S:D253G, S:D614G, ORF1ab:L3201P, ORF1ab:T265I, ORF3a:Q57H, S:T95I, ORF3a:P42L, ORF8:T11I".split(", ") #ORF1b:Q1011H,  ORF1b:P314L,
)

d614g = ProteinMutationIdentifier(
    gene = "S",
    ref = "D",
    position = 614,
    alt = "G",
    note = "May increase contagion and/or pathogenicity",
    source = "https://www.lanl.gov/updates/sars-cov-2-variant.php"
)

variantsOfConcern = {
    "strains": {
        "B.1.1.7": b_1_1_7.dictionary,
        "B.1.351": b_1_351.dictionary,
        "P.1": p_1.dictionary,
        "B.1.429": b_1_429.dictionary,
        "B.1.427": b_1_427.dictionary,
        "B.1.526": b_1_526.dictionary
    },
    "molecular": {
        "S:D614G": d614g.dictionary
    }
}

if __name__ == "__main__":
    outputFilePath = os.path.join(outputFolder, "variantsOfConcern.json")
    outputFile = open(outputFilePath, 'w')
    json.dump(variantsOfConcern, outputFile, indent=4)
    outputFile.close()