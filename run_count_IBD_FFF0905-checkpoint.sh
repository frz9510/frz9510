#! /bin/sh

#SBATCH --account=def-girardsi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-10:00

module load StdEnv/2020
module load plink/2.00-10252019-avx2
module load scipy-stack/2023a
module load r/4.1.2

mydir=$PWD
cd ${mydir}

data_dir='/lustre03/project/6033529/ARSACS/data/genotypes'
plink_file='ARSACS_CaG_Sag_commonsnps_geno0.05'
length=2

for ichr in {1..22}

ibd_gz_file='ARSACS_CaG_Sag_commonsnps_geno0.05_refinedIBD.calculated_chr'${ichr}'_window40_scale10_lod3_length1.ibd.gz'

do

	ibd_file_tmp=${ibd_gz_file}_chr${ichr}_window40_scale10_lod3_length1.ibd.gz
	ibd_file=${ibd_gz_file}_chr${ichr}_window40_scale10_lod3_length${length}.ibd
    
	zcat ${ibd_file_tmp} | awk -v len=${length} '$9>len'  > ${ibd_file}


	python grep_patients_in_ibd_file_CM.py ${ibd_file} ${plink_file}.fam
	plink2 --bfile ${plink_file} --chr ${ichr} --make-just-bim --out ${ibd_gz_file}_chr${ichr}

	python IBDsharing_byPosition.py ${ibd_file}.patients ${ichr} ${ibd_gz_file}_chr${ichr}.bim
	python IBDsharing_byPosition.py ${ibd_file}.patientsVSctrls ${ichr} ${ibd_gz_file}_chr${ichr}.bim
	python IBDsharing_byPosition.py ${ibd_file}.pop.ctrls ${ichr} ${ibd_gz_file}_chr${ichr}.bim
done


