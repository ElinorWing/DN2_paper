#!/bin/bash
#
# Execute in the current working directory
#$ -cwd
#
#set environment
#
#max runtime of 48 hours
#$ -V
#$ -l h_vmem=16G
#$ -pe sharedmem 8
#$ -l h_rt=47:59:59
#
# Task range. Tasks need to go from 1 to the number of files:
#$ -t 1-3

# load modules and conda env
module load igmm/apps/samtools/1.6
module load anaconda
source activate velocity


# example from eddie wiki
# SGE_TASK_ID will go from 1 to the number of files when we submit an array job
sample_name=`sed -n ${SGE_TASK_ID}p < filename_list.txt`
echo Processing sample: ${sample_name} on $HOSTNAME


#run velocyto
velocyto run -vvv -o /home/s2126362/scratch/velocyto/output -b /home/s2126362/scratch/barcodes/${sample_name}_barcodes.tsv.gz /home/s2126362/scratch/bam/${sample_name}_sample_alignments.bam /home/s2126362/scratch/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf


#deactivate velocity environment
source deactivate

