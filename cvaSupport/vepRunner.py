import os
from shlex import quote as shlex_quote

vepPath = "/opt/vep/src/ensembl-vep/vep"
gtfPath = "/opt/vep/ronavep/references/Sars_cov_2.ASM985889v3.101.gtf.named.bgzipSort.gz"
synonymsPath = "/opt/vep/ronavep/references/sars_cov_2_ASM985889v3_chr_synonyms.txt"
fastaPath = "/opt/vep/ronavep/references/Sars_cov_2.ASM985889v3.dna_sm.toplevel.fa.gz"
geneRadius = 1


def runVEP(vcfPath:str, outputFolder:str=""):
    vcfName = os.path.split(vcfPath)[1]
    if vcfName.lower().endswith(".vcf.gz"):
        vcfName = ".".join(vcfName.split(".")[:-2])
    if vcfName.lower().endswith(".vcf"):
        vcfName = ".".join(vcfName.split(".")[:-1])
    outputFileName = vcfName + ".vep.txt"
    outputFilePath = os.path.join(outputFolder, outputFileName)
    if not os.path.isdir(outputFolder):
        os.mkdir(outputFolder)
    vepCommand = "%s -i %s  -gtf %s  -fasta %s  -synonyms %s  -distance %s  -o %s  --force_overwrite" %(vepPath, vcfPath, gtfPath, fastaPath, synonymsPath, geneRadius, outputFilePath)
    print("RUN: %s" %vepCommand)
    exitStatus = os.system(vepCommand)
    if exitStatus == 0:
        print("VEP Successful")
    else:
        print("WARNING: VEP EXITED WITH A NON-ZERO STATUS OF %s" %exitStatus)
    if os.path.isfile(outputFilePath):
        return outputFilePath
    else:
        return ""