#! /bin/sh

#SBATCH --account=def-girardsi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=0-12:00


module load StdEnv/2023 gcc/12.3
module load r/4.3
module load r-bundle-bioconductor/3.18
module load scipy-stack/2023b
module load haploview/4.2
#java -jar $EBROOTHAPLOVIEW/Haploview.jar

data_dir=/lustre09/project/6033529/ARSACS/results/IBD2/
working_dir=/lustre09/project/6033529/ARSACS/results/haplotypes2/
cd ${working_dir}

my_chr=13
## Your positiions for hapltypes around ARSACS gene if you want to send on the command line
start=19121958
end=115106996


python /lustre09/project/6033529/ARSACS/scripts/make_phased_vcf_as_GenHunterTDT_FFF2.py ${my_chr} ${start} ${end}

awk '{ 
    printf "%s %s", $1, $1; 
    for (i=2; i<=NF; i++) printf " %s", $i; 
    printf "\n"; 
}' ARSACS_CaG_sag_chr${my_chr}_${start}-${end}_onlyhet_haps_tmp.haps > ARSACS_CaG_sag_chr${my_chr}_${start}-${end}_onlyhet_haps_tmp2.haps


zgrep -v "##" ARSACS_CaG_Sag_commonsnps_geno0.05_noATGC_refinedIBD.calculated_chr13.vcf.gz | \
awk -v s=19121958 -v e=115106996 '{ if ($2 >= s && $2 <= e) print $0 }' | head


#export JAVA_TOOL_OPTIONS="-Xmx4g"#
#export _JAVA_OPTIONS="-Xmx4g"#
#java -jar $EBROOTHAPLOVIEW/Haploview.jar -nogui -haps ${working_dir}ARSACS_CaG_sag_chr${my_chr}_${start}-${end}_onlyhet_haps_tmp2.haps -info ${working_dir}ARSACS_CaG_sag_chr${my_chr}_${start}-${end}_onlyhet_haps.info -blockoutput ALL -dprime -compressedpng -hapthresh 0






## Have to define files in the python script
#python /lustre03/project/6033529/ARSACS/scripts/make_phased_vcf_as_haps-Copy1.py ${data_dir}ARSACS_CaG_Sag_commonsnps_geno0.05_noATGC_refinedIBD.calculated_chr13.vcf.gz

#ARSACS_CaG_Sag_commonsnps_geno0.05_noATGC_refinedIBD.calculated_chr21.vcf.gz
## Also have to define files output from python to the R scripts
#Rscript /lustre03/project/6033529/ARSACS/scripts/haplotype_network-Copy1.r 
