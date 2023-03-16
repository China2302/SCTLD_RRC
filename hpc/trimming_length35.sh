#!/bin/bash
#BSUB -J trim_l35
#BSUB -q bigmem
#BSUB -P c_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o trim_l35%J.out
#BSUB -e trim_l35%J.err
#BSUB -u nxa945@miami.edu
#BSUB -N


ln -s ${todata}/data/trimmed_l35 trimmed_l35

cd "/nethome/nxa945/trimmed_l35"


todata="/scratch/projects/c_transcriptomics/NOVA_SCTLD/"

for sample in ${todata}/raw_data/*.gz ;

do \
trim_galore ${sample} \
--gzip \
--fastqc \
--fastqc_args "--outdir ${todata}/data/trimmed_l35/" \
--illumina \
--cores 4 \
--three_prime_clip_R1 12 \
--nextseq 30 \
--length 35 \
-o ${todata}/data/trimmed_l35/ ; \

done

multiqc ${todata}/data/trimmed_l35/ \
--outdir ${todata}/data/trimmed_l35/
