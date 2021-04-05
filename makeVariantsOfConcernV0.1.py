import os
from cvaSupport.viralVariantHandler import StrainIdentifier, ProteinMutationIdentifier, NucleicAcidMutationIdentifier
import json


outputFolder = os.path.join(os.path.split(__file__)[0], "references")

b_1_1_7 = StrainIdentifier(
    identifier = "B.1.1.7",
    commonName = "UK",
    location = "UK",
    aliases = ["UK", "20I/501Y.V1"],
    source = "https://www.cdc.gov/mmwr/volumes/70/wr/mm7003e2.htm#T1_down",
    proteinChanges=
        ["ORF1ab:" + variant for variant in "T1001I, A1708D, I2230T, del3675-3677 SGF".split(", ")] +
        ["S:" + variant for variant in "del69-70 HV, del144 Y, N501Y, A570D, D614G, P681H, T761I, S982A, D1118H".split(", ")] +
        ["ORF8:" + variant for variant in "Q27*, R52I, Y73C".split(", ")] +
        ["N:" + variant for variant in "D3L, S235F".split(", ")]
)

b_1_351 = StrainIdentifier(
    identifier = "B.1.351",
    commonName = "South Africa",
    location = "South Africa",
    aliases = ["South Africa", "20H/501Y.V2"],
    source="https://www.cdc.gov/mmwr/volumes/70/wr/mm7003e2.htm#T1_down",
    proteinChanges=
        ["ORF1ab:" + variant for variant in "K1655N".split(", ")] +
        ["S:" + variant for variant in "K417N, E484K, N501Y, D614G, A701V".split(", ")] +
        ["E:" + variant for variant in "P71L".split(", ")] +
        ["N:" + variant for variant in "T205I".split(", ")]
)

p_1 = StrainIdentifier(
    identifier = "P.1",
    commonName = "Brazil",
    location = "Brazil",
    aliases = ["Brazil", "Japan", "20J/501Y.V3"],
    source = "https://www.cdc.gov/mmwr/volumes/70/wr/mm7003e2.htm#T1_down",
    proteinChanges=
        ["ORF1ab:" + variant for variant in "F681L, I760T, S1188L, K1795Q, del3675-3677 SGF, E5662D".split(", ")] +
        ["S:" + variant for variant in "L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, D614G, H655Y, T1027I".split(", ")] +
        ["ORF3a:" + variant for variant in "C174G".split(", ")] +
        ["ORF8:" + variant for variant in "E92K".split(", ")] +
        ["ORF9:" + variant for variant in "Q77E".split(", ")] +
        ["ORF14:" + variant for variant in "V49L".split(", ")] +
        ["N:" + variant for variant in "P80R".split(", ")]
)

b_1_429 = StrainIdentifier(
    identifier = "B.1.429",
    commonName = "California",
    location = "California",
    aliases = ["California", "CAL.20C"],
    source = "https://www.cdc.gov/mmwr/volumes/70/wr/mm7003e2.htm#T1_down",
    proteinChanges=
        ["ORF1ab:" + variant for variant in "I4205V,D1183Y".split(", ")] +
        ["S:" + variant for variant in "S13I, W152C, L452R".split(", ")]
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
        "B.1.429": b_1_429.dictionary
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