#!/bin/bash
#BSUB -J trim_l35_polyA
#BSUB -q bigmem
#BSUB -P c_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o trim_l35_polyA%J.out
#BSUB -e trim_l35_polyA%J.err
#BSUB -u nxa945@miami.edu
#BSUB -N


todata="/scratch/projects/c_transcriptomics/NOVA_SCTLD/"

cd "/nethome/nxa945/trimmed_l35_polyA"

for sample in ${todata}/data/trimmed_l35/trimmed_l35_data/*.gz ;

do \
trim_galore ${sample} \
--polyA ; \

done
