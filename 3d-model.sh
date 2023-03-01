#!/bin/bash

#SBATCH --account=
#SBATCH --nodes=
#SBATCH --ntasks=
#SBATCH --mem=
#SBATCH --qos=
#SBATCH --partition=
#SBATCH --time=14-23:59:59
#SBATCH --array=1-40

############################
#SBATCH --job-name=chrom3d-model
#SBATCH --output=%j-chrom3d-model.out
############################


SAMPLE="samplename"
SAMPLE_SLURM="${SLURM_ARRAY_TASK_ID}/${SAMPLE}"

echo "Run Chrom3D to generate 3D genome models"
echo "random seed ${SLURM_ARRAY_TASK_ID}"
Chrom3D -s $SLURM_ARRAY_TASK_ID -y 0.15 -r 5.0 -n 2000000 --nucleus \
    -o ${SAMPLE_SLURM}_inter_intra_chr.diploid.cmm \
    ${SAMPLE}_inter_intra_chr.diploid.gtrack

echo "get coordinates"
grep "marker" ${SAMPLE_SLURM}_inter_intra_chr.diploid.cmm \
    | grep "_A" \
	| awk -F'"' '{print $20"\t"$4"\t"$6"\t"$8}' \
    | awk -F'_A:' '{print $1"\t"$2}' \
    | sed 's/-/\t/' > ${SAMPLE_SLURM}_coords.bed

echo "sort coordinated bed"
sort -k1,1 -k2,2n ${SAMPLE_SLURM}_coords.bed > ${SAMPLE_SLURM}_coords_sorted.bed

echo "DONE"