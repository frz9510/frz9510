#!/bin/bash
cd /lustre03/project/6033529/ARSACS/scripts
source /lustre03/project/6033529/IBD_denovo/scripts/ibdEnv/bin/activate

#patients_file=/lustre03/project/6033529/IBD_denovo/data/HPT_patients_clean_sureHPT_45-47mb.txt
patients_file=/lustre03/project/6033529/IBD_denovo/data/HG1_ind2keep.txt #HGA
#patients_file=/lustre03/project/6033529/IBD_denovo/data/HG2_ind2keep.txt #HGB
#patients_file=/lustre03/project/6033529/IBD_denovo/data/HG3_ind2keep.txt #HGU
hpt_file=/lustre03/project/6033529/DM1/results/genotypes/haplotypes/haps_patients_ctrls_chr19_45206933-47625828.txt
blocks_file=/lustre03/project/6033529/DM1/results/genotypes/haplotypes/dm1_CaG_sag_chr19_45206933-47625828_onlyhet_haps.haps.4GAMblocks
outdir=../results/tests

python DM1_haploblock_graph_FFF.py \
--patients ${patients_file} \
--hpt_file ${hpt_file} \
--out_dir ../results/tests \
--blocks_file ${blocks_file} \
--red_line 46273462 \
--prefix FINAL_HGU_FFF_ \
--main_color '#993F00' \
#--single_hpt H1872 \
#--percent ../HPT_patients_clean_sureHPT_45-47mb_percent.txt