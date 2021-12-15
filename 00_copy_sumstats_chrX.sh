#!/bin/bash
#-------------------------------------------------------------------------------
# SLURM command to run
#-------------------------------------------------------------------------------
# sbatch --time=24:00:00 --error errandout/copy_sumstats.e --output errandout/copy_sumstats.o 00_copy_sumstats.sh

#-------------------------------------------------------------------------------
# Description
#-------------------------------------------------------------------------------
# This script copies all of the summary statistic files from the pgc DAC account
# into a local sumstats directory (originally created to be run on pgca1pts on
# the LISA supercomputer system)

# change this to be the directory you are running the gwas in (sumstats will be
# copied to the "sumstats" folder inside this working directory)
WORKING_DIR=/home/pgca1pts/freeze3_gwas
cd $WORKING_DIR

# don't change this
ss_dir=${WORKING_DIR}/sumstats
mkdir -p ${ss_dir}
DAC_DIR=/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts



#-------------------------------------------------------------------------------
# VETS & UKBB
#-------------------------------------------------------------------------------
# VETS & UKBB - run these two studies with BOLT LMM
# VETS - ptsd_qt_vetsa_may12_2021_related_filtered.imputed.stats.gz
# UKBB - pts_ukbb_may13_2021_unrelated.bgen.stats.gz

#-------------------------------------------------------------------------------
# VETS
#-------------------------------------------------------------------------------
# Study only has males - use same file for both
echo VETS
dir=/home/pgca1pts/freeze3_bolt_data
zcat ${dir}/ptsd_qt_vetsa_may12_2021_related_filtered.imputed.stats.gz \
  | awk 'BEGIN{OFS="\t"}{if (NR == 1){print}else if($2 == 23 && $7 >= 0.01 && $7 <= 0.99){$2="X"; print}}' \
  | gzip > ${ss_dir}/VETS_broad_chrX_eur.txt.gz
# male gwas
cp ${ss_dir}/VETS_broad_chrX_eur.txt.gz \
  ${ss_dir}/VETS_males_chrX_eur.txt.gz

#-------------------------------------------------------------------------------
# UKBB
#-------------------------------------------------------------------------------
# Study only has males - use same file for both
echo UKBB
dir=/home/pgca1pts/freeze3_bolt_data
# base gwas
cp ${dir}/pts_ukbb_may13_2021_unrelated_chrX.bgen.stats.gz ${TMPDIR}/UKBB_broad_chrX_eur.txt.gz
zcat ${TMPDIR}/UKBB_broad_chrX_eur.txt.gz | awk 'BEGIN{OFS="\t"}{if(NR == 1){print}else if($7 >= 0.01 && $7 <= 0.99){$2="X"; print}}' \
  | gzip > ${ss_dir}/UKBB_broad_chrX_eur.txt.gz
# male gwas
# female gwas


#-------------------------------------------------------------------------------
# GROUP 3
#-------------------------------------------------------------------------------
# GROUP 3 sumstats - freeze 2 summary stats


# No MIRE, INTR, DAMI


#-------------------------------------------------------------------------------
# QIMR
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
echo QIMR
# base gwas
zcat ${ss_dir}/QIMR_broad_eur.txt.gz | awk '{if(NR == 1 || $1 == "X"){print}}' \
  | gzip > ${ss_dir}/QIMR_broad_chrX_eur.txt.gz
# male gwas
zcat ${ss_dir}/QIMR_males_eur.txt.gz | awk '{if(NR == 1 || $1 == "X"){print}}' \
  | gzip > ${ss_dir}/QIMR_males_chrX_eur.txt.gz
# female gwas
zcat ${ss_dir}/QIMR_females_eur.txt.gz | awk '{if(NR == 1 || $1 == "X"){print}}' \
  | gzip > ${ss_dir}/QIMR_females_chrX_eur.txt.gz


# No NCPT, TRAC

#-------------------------------------------------------------------------------
# GROUP 4
#-------------------------------------------------------------------------------
# GROUP 4 sumstats - freeze 3 summary stats


#-------------------------------------------------------------------------------
# WTCS
#-------------------------------------------------------------------------------
# Annotate RS IDs
# Filter to minor allele frequency >= 0.01
echo WTCS
# base gwas
zcat ${ss_dir}/WTCS_broad_eur.txt.gz | awk '{if(NR == 1 || $3 == "X"){print}}' \
  | gzip > ${ss_dir}/WTCS_broad_chrX_eur.txt.gz
# male gwas
zcat ${ss_dir}/WTCS_males_eur.txt.gz | awk '{if(NR == 1 || $1 == "X"){print}}' \
  | gzip > ${ss_dir}/WTCS_males_chrX_eur.txt.gz
# female gwas
zcat ${ss_dir}/WTCS_females_eur.txt.gz | awk '{if(NR == 1 || $1 == "X"){print}}' \
  | gzip > ${ss_dir}/WTCS_females_chrX_eur.txt.gz


# No AGDS


#-------------------------------------------------------------------------------
# CANA
#-------------------------------------------------------------------------------
# Already filtered to minor allele frequency >= 0.01
echo CANA
# base gwas
zcat ${ss_dir}/CANA_broad_eur.txt.gz | awk '{if(NR == 1 || $1 == "X"){print}}' \
  | gzip > ${ss_dir}/CANA_broad_chrX_eur.txt.gz
# male gwas - no sex specific
#cp ${dir}/clsa_1-23xm2_covars_c4_20210512m.PSD_DCTOFF_COM.assoc.logistic.add.gz \
#  ${ss_dir}/CANA_males_eur.txt.gz
# female gwas - no sex specific
#cp ${dir}/clsa_1-23xm2_covars_c4_20210512f.PSD_DCTOFF_COM.assoc.logistic.add.gz \
#  ${ss_dir}/CANA_females_eur.txt.gz


# No QIM2, RCOG, WTCM


#-------------------------------------------------------------------------------
# DAI2
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
# Change 23 to X
echo DAI2
dir=${DAC_DIR}/wave3/summary_stats/74_dai2
# base gwas
cp ${dir}/daner_iPSYCH2015_PTSDbroad_chrX_HRC_MAF01.gz \
  ${TMPDIR}/DAI2_broad_chrX_eur.txt.gz
zcat ${TMPDIR}/DAI2_broad_chrX_eur.txt.gz | awk 'BEGIN{OFS=" "}{if($1 == 23){$1 = "X"}else{$1 = $1} print}' \
  | gzip > ${ss_dir}/DAI2_broad_chrX_eur.txt.gz
# male gwas
cp ${dir}/daner_iPSYCH2015_PTSDbroad_males_chrX_HRC_MAF01.gz \
  ${TMPDIR}/DAI2_males_chrX_eur.txt.gz
zcat ${TMPDIR}/DAI2_males_chrX_eur.txt.gz | awk 'BEGIN{OFS=" "}{if($1 == 23){$1 = "X"}else{$1 = $1} print}' \
  | gzip > ${ss_dir}/DAI2_males_chrX_eur.txt.gz
# female gwas
cp ${dir}/daner_iPSYCH2015_PTSDbroad_females_chrX_HRC_MAF01.gz \
  ${TMPDIR}/DAI2_females_chrX_eur.txt.gz
zcat ${TMPDIR}/DAI2_females_chrX_eur.txt.gz | awk 'BEGIN{OFS=" "}{if($1 == 23){$1 = "X"}else{$1 = $1} print}' \
  | gzip > ${ss_dir}/DAI2_females_chrX_eur.txt.gz
# clean TMPDIR
rm ${TMPDIR}/DAI2*


# No BIOV, MGBB


#-------------------------------------------------------------------------------
# HUNT
#-------------------------------------------------------------------------------
# Filter to minor allele frequency >= 0.01
# Annotate RS IDs
echo HUNT
# base gwas
zcat ${ss_dir}/HUNT_broad_eur.txt.gz | awk '{if(NR == 1 || $3 == "X"){print}}' \
  | gzip > ${ss_dir}/HUNT_broad_chrX_eur.txt.gz

#-------------------------------------------------------------------------------
# SWED
#-------------------------------------------------------------------------------
# Already filtered to minor allele frequency >= 0.01
echo SWED
# base gwas
zcat ${ss_dir}/SWED_broad_eur.txt.gz | awk '{if(NR == 1 || $1 == "X"){print}}' \
  | gzip > ${ss_dir}/SWED_broad_chrX_eur.txt.gz

#-------------------------------------------------------------------------------
# FING
#-------------------------------------------------------------------------------
# Annotate RS IDs, liftover to hg19
# Filter to minor allele frequency >= 0.01
echo FING
# base gwas
zcat ${ss_dir}/FING_broad_eur.txt.gz | awk '{if(NR == 1 || $2 == "X"){print}}' \
  | gzip > ${ss_dir}/FING_broad_chrX_eur.txt.gz
# male gwas
zcat ${ss_dir}/FING_males_eur.txt.gz | awk '{if(NR == 1 || $2 == "X"){print}}' \
  | gzip > ${ss_dir}/FING_males_chrX_eur.txt.gz
# female gwas
zcat ${ss_dir}/FING_females_eur.txt.gz | awk '{if(NR == 1 || $2 == "X"){print}}' \
  | gzip > ${ss_dir}/FING_females_chrX_eur.txt.gz


#-------------------------------------------------------------------------------
# UKB2
#-------------------------------------------------------------------------------
# Copy X chromsome data
dir=${DAC_DIR}/wave3/summary_stats/91_ukb2
cp ${dir}/PGC_Return_2021_11_23.zip ${TMPDIR}
cd $TMPDIR
unzip PGC_Return_2021_11_23.zip
cd Broad_Data_Analysis
# base gwas
bzip2 -d PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_With_X_WG.txt_Broad.regenie.bz2
awk 'BEGIN{OFS=" "}{if(NR == 1){print}else if($1 == 23){$1="X"; print}}' PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_With_X_WG.txt_Broad.regenie \
  | gzip > ${ss_dir}/UKB2_broad_chrX_eur.txt.gz
# male gwas
bzip2 -d PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_Male_With_X_WG.txt_Broad.regenie.bz2
awk 'BEGIN{OFS=" "}{if(NR == 1){print}else if($1 == 23){$1="X"; print}}' PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_Male_With_X_WG.txt_Broad.regenie \
  | gzip > ${ss_dir}/UKB2_males_chrX_eur.txt.gz
# female gwas
bzip2 -d PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_Female_With_X_WG.txt_Broad.regenie.bz2
awk 'BEGIN{OFS=" "}{if(NR == 1){print}else if($1 == 23){$1="X"; print}}' PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_Female_With_X_WG.txt_Broad.regenie \
  | gzip > ${ss_dir}/UKB2_females_chrX_eur.txt.gz
# clean TMPDIR
rm -r ${TMPDIR}/Broad_Data_Analysis
# cd back into working directory
cd $WORKING_DIR


# No BIOM, ESBB


#-------------------------------------------------------------------------------
# MAYO
#-------------------------------------------------------------------------------
# Annotate RS IDs
dir=${DAC_DIR}/wave3/summary_stats/99_mayo
# base gwas
zcat ${ss_dir}/MAYO_broad_eur.txt.gz | awk '{if(NR == 1){$3="ID";print}else if($1 == 23){$1="X";$3=$1":"$2":"$5":"$6;print}}' \
  | gzip > ${TMPDIR}/MAYO_broad_chrX_eur.txt.gz
# male gwas
zcat ${ss_dir}/MAYO_males_eur.txt.gz | awk '{if(NR == 1){$3="ID";print}else if($1 == 23){$1="X";$3=$1":"$2":"$5":"$6;print}}' \
  | gzip > ${TMPDIR}/MAYO_males_chrX_eur.txt.gz
# female gwas
zcat ${ss_dir}/MAYO_females_eur.txt.gz | awk '{if(NR == 1){$3="ID";print}else if($1 == 23){$1="X";$3=$1":"$2":"$5":"$6;print}}' \
  | gzip > ${TMPDIR}/MAYO_females_chrX_eur.txt.gz

# Annotate RS IDs
for run in "broad" "males" "females"; do
  echo annotating $run
  if [ $run == "broad" ]; then
    base=MAYO_broad_chrX_eur
  elif [ $run == "males" ]; then
    base=MAYO_males_chrX_eur
  else
    base=MAYO_females_chrX_eur
  fi
  cat vcf_header1b.txt vcf_header2.txt <(zcat ${TMPDIR}/${base}.txt.gz | awk '{OFS="\t"}{if (NR>1) print $1,$2, ".", $5,$6 , "100", "PASS", "MVP="$3}' \
    | LC_ALL=C sort -g -k1r,1 -k2,2 ) > ${TMPDIR}/${base}.vcf
  bcftools view ${TMPDIR}/${base}.vcf -Oz -o ${TMPDIR}/${base}.vcf.gz
  tabix -f ${TMPDIR}/${base}.vcf.gz
  bcftools annotate -r X -a /home/pgca1pts/All_20180423_hg19.vcf.gz -c INFO ${TMPDIR}/${base}.vcf.gz -o ${TMPDIR}/${base}.vcf.annotated
  # header start with # so it stays at top after sorting
  echo -e "#SNPID\tSNP" > ${TMPDIR}/${base}_snpheader.txt
  tail -n+82 ${TMPDIR}/${base}.vcf.annotated  | grep RS \
    | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' \
    | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat ${TMPDIR}/${base}_snpheader.txt -  | LC_ALL=C sort -k1b,1 \
    > ${TMPDIR}/${base}.vcf.annotated.success
  tail -n+82 ${TMPDIR}/${base}.vcf.annotated | grep -v RS \
    | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}' \
    | sed 's/MVP=//g' | cat ${TMPDIR}/${base}_snpheader.txt - | LC_ALL=C sort -k1b,1 \
    > ${TMPDIR}/${base}.vcf.annotated.failed
  LC_ALL=C join -1 1 -2 3  ${TMPDIR}/${base}.vcf.annotated.success \
    <(zcat ${TMPDIR}/${base}.txt.gz | sed s/ID/\#SNPID/g | LC_ALL=C sort -k3b,3 ) \
    > ${TMPDIR}/${base}.txt.success
  LC_ALL=C join -1 1 -2 3  ${TMPDIR}/${base}.vcf.annotated.failed \
    <(zcat ${TMPDIR}/${base}.txt.gz | sed s/ID/\#SNPID/g | LC_ALL=C sort -k3b,3 ) \
    > ${TMPDIR}/${base}.txt.failed
  echo "annotation success:"
  wc -l ${TMPDIR}/${base}.txt.success
  echo "annotation failed:"
  wc -l ${TMPDIR}/${base}.txt.failed
  cat ${TMPDIR}/${base}.txt.success ${TMPDIR}/${base}.txt.failed \
    | LC_ALL=C sort -g -k 3  | awk '{if(NR == 2  || ($5 >= 0.01 && $5 <= 0.99)) print $2,$1,$3,$4,$6,$7,$5,$8,$9,$10,$11,$12,$13}' \
    | grep -v NA | sed s/\#SNPID/SNPID/g | gzip > ${ss_dir}/${base}.txt.gz
done
# clean TMPDIR
rm ${TMPDIR}/MAYO*

#-------------------------------------------------------------------------------
# GROUP 5
#-------------------------------------------------------------------------------
# GROUP 5 sumstats - MVP summary stats


# No MVP
