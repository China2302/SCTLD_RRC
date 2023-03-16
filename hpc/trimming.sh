#!/bin/bash
#BSUB -J trim_all
#BSUB -q bigmem
#BSUB -P c_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o trim_all%J.out
#BSUB -e trim_all%J.err
#BSUB -u nxa945@miami.edu
#BSUB -N



todata="/scratch/projects/c_transcriptomics/NOVA_SCTLD/"

for sample in ${todata}/raw_data/*.gz ;

do \
trim_galore ${sample} \
--gzip \
--fastqc \
--fastqc_args "--outdir ${todata}/data/trimmed/" \
--illumina \
--cores 4 \
--three_prime_clip_R1 12 \
--nextseq 30 \
--length 20 \
-o ${todata}/data/trimmed/ ; \

done

multiqc ${todata}/data/trimmed/ \
--outdir ${todata}/data/trimmed/
