#! /bin/bash

gzip -d Sars_cov_2.ASM985889v3.101.primary_assembly.MN908947.3.gff3.gz
grep -v "#" Sars_cov_2.ASM985889v3.101.primary_assembly.MN908947.3.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > Sars_cov_2.ASM985889v3.101.primary_assembly.MN908947.3.gff3.bgzipSort.gz
tabix -p gff Sars_cov_2.ASM985889v3.101.primary_assembly.MN908947.3.gff3.bgzipSort.gz
gzip -d Sars_cov_2.ASM985889v3.101.gtf.gz
grep -v "#" Sars_cov_2.ASM985889v3.101.gtf | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > Sars_cov_2.ASM985889v3.101.gtf.gz.bgzipSort.gz
tabix -p gff Sars_cov_2.ASM985889v3.101.gtf.gz.bgzipSort.gz
grep -v "#" Sars_cov_2.ASM985889v3.101.named.gtf | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > Sars_cov_2.ASM985889v3.101.gtf.named.bgzipSort.gz
tabix -p gff Sars_cov_2.ASM985889v3.101.gtf.named.bgzipSort.gz
gzip -d Sars_cov_2.ASM985889v3.dna_sm.toplevel.fa.gz
bgzip Sars_cov_2.ASM985889v3.dna_sm.toplevel.fa