#! /usr/bin/env bash

#BSUB -J cellranger
#BSUB -o logs/vdj_%J.out
#BSUB -e logs/vdj_%J.err
#BSUB -R "select[mem>30] rusage[mem=30]"
#BSUB -q rna

module load cellranger

cellranger vdj \
    --id test \
    --reference=/beevol/home/rbilab/ref/cellranger/mouse/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
    --fastqs=/beevol/home/wellskri/Analysis/Mia_Smith/210825_b_cell_analysis/results_vdj/fastqs \
    --sample=healthy_bcells_all_VDJ \
    --chain=IG \
    --jobmode=lsf.template \
    --maxjobs=10