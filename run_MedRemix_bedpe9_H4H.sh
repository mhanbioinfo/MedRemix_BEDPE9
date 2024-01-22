#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH --job-name=MedRemixBEDPE_bedpe9
#SBATCH --time=1-00:00:00
#SBATCH --account=pughlab
#SBATCH --mem=30G
#SBATCH --partition=himem
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=./logs_slurm/%j-%x.out
#SBATCH --error=./logs_slurm/%j-%x.err

module load R/4.2.1

echo "Job started at "$(date) 
time1=$(date +%s)


## set these paths
BEDPE9_FPATH="./GSM5067062_2020_6544_human_bedpe.1000lines.bed.gz"
# BEDPE9_FPATH="./GSM5067062_2020_6544_human_bedpe.bed.gz"
OUT_DIR="./outputs"

## set these resource paths
SRC_DIR="./scripts"
GENOME_DIR="/cluster/projects/pughlab/references/cfMeDIPseq_refs/medremix_ref/bsgenome/BSgenome.Hsapiens.UCSC.hg38"


## script
FNAME="$(basename ${BEDPE9_FPATH} .bed.gz)"
mkdir -p ${OUT_DIR}

## convert BEDPE9 to format accepted by MedRemix
bash ${SRC_DIR}/bedpe9_to_bedpe4medremix.sh \
${BEDPE9_FPATH} \
${OUT_DIR}/${FNAME}.bedpe4medremix

## bin_stat each chromosome
mkdir -p ${OUT_DIR}/bin_stats
chr_arr=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
# chr_arr=("chr1" "chr2")
for CHR in "${chr_arr[@]}"; do
    echo ${CHR}
    
    Rscript ${SRC_DIR}/bin_stats_bedpe.R \
    -i ${OUT_DIR}/${FNAME}.bedpe4medremix \
    -g ${GENOME_DIR} \
    -c ${CHR} \
    -o "${OUT_DIR}/bin_stats/bedpe_bin_stats_${FNAME}_${CHR}.tsv"

done

## merge bin stats across all chromosomes
Rscript ${SRC_DIR}/row_bind_tables.R \
--glob "${OUT_DIR}/bin_stats/bedpe_bin_stats_${FNAME}_chr*.tsv" \
-o ${OUT_DIR}/bedpe_bin_stats_${FNAME}.feather \
--in-tsv \
--out-feather \
--omit-paths

## run MedRemix
Rscript ${SRC_DIR}/cfmedip_nbglm.R \
-i ${OUT_DIR}/bedpe_bin_stats_${FNAME}.feather \
-o ${OUT_DIR}/bedpe_${FNAME}_fit_nbglm.tsv \
--modelout ${OUT_DIR}/bedpe_${FNAME}_fit_nbglm_model.Rds



time2=$(date +%s)
echo "Job ended at "$(date)
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
