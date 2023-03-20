#!/bin/bash
#BSUB -J star_align_l35_polyA
#BSUB -q bigmem
#BSUB -P c_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o star_align_l35_polyA%J.out
#BSUB -e star_align_l35_polyA%J.err
#BSUB -u nxa945@miami.edu
#BSUB -N


todata="/scratch/projects/c_transcriptomics/"


cd "/nethome/nxa945/forStar_l35_polyA"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
STAR \
--runMode alignReads \
--genomeDir ${todata}/Genomes/Ofav/Ofav_index \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile  ${todata}/Genomes/Ofav/genomic.gtf \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 16 \
--readFilesCommand gunzip -c \
--readFilesIn ${sample} \
--outFilterMultimapNmax 20 \
--quantMode TranscriptomeSAM GeneCounts \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${todata}/NOVA_SCTLD/data/alligned_l35_polyA/${sample} ; \

done

multiqc ${todata}/NOVA_SCTLD/data/alligned_l35_polyA/ \
--interactive \
--outdir ${todata}/NOVA_SCTLD/data/alligned_l35_polyA/
