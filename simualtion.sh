
# simulation with reference genome
boquila input_reads.fa \
    --bed output_reads.bed \
    --fasta \
    --ref grch38_genome.fa \
    --regions Grch38.ron \
    --seed 7 \
    > output_reads.fa

# simulation with input DNA sequencing
boquila input_reads.fa \
    --fasta \
    --inseq input_DNA_sequencing.fq \
    --seed 7 \
    > output_reads.fa