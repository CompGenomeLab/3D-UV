#!/bin/bash

#SBATCH --account=
#SBATCH --nodes=
#SBATCH --ntasks=
#SBATCH --mem=
#SBATCH --qos=
#SBATCH --partition=
#SBATCH --time=13:49:59

############################
#SBATCH --job-name=chrom3d-preprocess
#SBATCH --output=%j-chrom3d-preprocess.out
############################


PPS="chrom3D/preprocess_scripts" # from chrom 3D pipeline
GENOME_SIZE="hg38.genome"


SAMPLE="samplename"

matrices="pathtohicmatrices"

matrix_50kb="${matrices}/${SAMPLE}.50kb.hicpro"
matrix_1mb="${matrices}/${SAMPLE}.1mb.hicpro"
abs_50kb="${matrices}/50kb_abs_chr.grch38.bed"
abs_1mb="${matrices}/1mb_abs_chr.grch38.bed"

tad="${SAMPLE}/tads_sorted.bed3"

cp $tad tad.bed

echo "convert hic-pro output"

echo "create intra-chromosomal contact matrices"
python ${PPS}/conv_hicpro_mat.py \
    $matrix_50kb \
    $abs_50kb \
    > ${SAMPLE}_50000.intermediate.bedpe

mkdir -p intra_chr_RAWobserved
bash ${PPS}/make_intrachr_rawObserved.sh \
    ${GENOME_SIZE} \
    ${SAMPLE}_50000.intermediate.bedpe

echo "create inter-chromosomal contact matrices"
python ${PPS}/conv_hicpro_mat.py \
    $matrix_1mb \
    $abs_1mb \
    > ${SAMPLE}_1000000.intermediate.bedpe


mkdir -p inter_chr_RAWobserved
bash ${PPS}/make_interchr_rawObserved.sh \
    ${GENOME_SIZE} \
    ${SAMPLE}_1000000.intermediate.bedpe


echo "TAD to domains, input is merged arrowhead out"
bash ${PPS}/arrowhead_to_domains.sh \
    tad.bed \
    ${GENOME_SIZE}


echo "Concatenate all the .domains to use in a later step"
cat *.chr*.domains > ${SAMPLE}_domainlist.domains


echo "Compute intra-chromosomal interaction counts between TADs"
mkdir -p intrachr_bedpe
bash ${PPS}/intrachr_NCHG_input_auto.sh \
    tad \
    ${GENOME_SIZE} \
    50kb

echo "Concatenate all intra-chromosomal interaction counts"
cat intrachr_bedpe/chr*.bedpe > intrachr_bedpe/${SAMPLE}_50kb.domain.RAW.bedpe


echo "Remove domains that contain centromeres from the BEDPE file"
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen | \
    pairToBed \
    -a intrachr_bedpe/${SAMPLE}_50kb.domain.RAW.bedpe \
    -b stdin -type neither > intrachr_bedpe/${SAMPLE}_50kb.domain.RAW.no_cen.bedpe


echo "Calculate the P-value and odds ratio for each pair of TADs"
${PPS}/NCHG_hic/NCHG -m 50000 \
    -p intrachr_bedpe/${SAMPLE}_50kb.domain.RAW.no_cen.bedpe \
    > ${SAMPLE}_50kb.domain.RAW.no_cen.NCHG.out


echo "Calculate FDR and filter significant interactions"
python ${PPS}/NCHG_fdr_oddratio_calc.py \
    ${SAMPLE}_50kb.domain.RAW.no_cen.NCHG.out fdr_bh 2 0.01 > ${SAMPLE}_50kb.domain.RAW.no_cen.NCHG.sig


echo "Create GTrack using significant interactions"
bash ${PPS}/make_gtrack.sh \
    ${SAMPLE}_50kb.domain.RAW.no_cen.NCHG.sig \
    ${SAMPLE}_domainlist.domains \
    ${SAMPLE}_intra_chromosome.gtrack


# # bedtools version needs to be 2.27.1 for this script to work
echo "Prepare inter-chromosomal Hi-C interaction counts"
bash ${PPS}/interchr_NCHG_input_auto.sh \
    ${GENOME_SIZE} \
    ${BLACKLIST} \
    1mb > ${SAMPLE}_1mb_inter.bedpe

echo "Call significant inter-chromosomal interactions"
${PPS}/NCHG_hic/NCHG \
    -i -p ${SAMPLE}_1mb_inter.bedpe > ${SAMPLE}_1mb_inter_chr.NCHG.out

echo "Calculate FDR and filter significant interactions"
python ${PPS}/NCHG_fdr_oddratio_calc.py \
    ${SAMPLE}_1mb_inter_chr.NCHG.out fdr_bh 2 0.01 > ${SAMPLE}_1mb_inter_chr.NCHG.sig

echo "Add significant inter-chromosomal interaction information to the GTrack"
bash ${PPS}/add_inter_chrom_beads_wo_lads.sh \
    ${SAMPLE}_intra_chromosome.gtrack \
    ${SAMPLE}_1mb_inter_chr.NCHG.sig \
    ${SAMPLE}_inter_intra_chr.gtrack

echo "Modify the Model Setup File to make a diploid model"
python ${PPS}/make_diploid_gtrack.py \
    ${SAMPLE}_inter_intra_chr.gtrack > ${SAMPLE}_inter_intra_chr.diploid.gtrack

echo "GTRACK sanity check"
grep -v '^#' ${SAMPLE}_inter_intra_chr.diploid.gtrack | cut -f 1 | sort | uniq -c

echo "DONE"