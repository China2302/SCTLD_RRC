#!/bin/bash
#BSUB -J star_index
#BSUB -q general
#BSUB -P c_transcriptomics
#BSUB -n 8
#BSUB -o /scratch/projects/c_transcriptomics/Genomes/Ofav/star_index%J.out
#BSUB -e /scratch/projects/c_transcriptomics/Genomes/Ofav/star_index%J.err
#BSUB -u nxa945@miami.edu
#BSUB -N



todata="/scratch/projects/c_transcriptomics"

STAR \
--runThreadN 16 NumberOfThreads \
--runMode genomeGenerate \
--genomeDir ${todata}/Genomes/Ofav/Ofav_index \
--genomeFastaFiles ${todata}/Genomes/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.fna \
--sjdbGTFfile /${todata}/Genomes/Ofav/genomic.gtf \
--sjdbOverhang ReadLength-100 \
--genomeSAindexNbases 13
