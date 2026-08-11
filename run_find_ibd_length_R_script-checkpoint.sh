#! /bin/sh

#SBATCH --account=def-girardsi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=120G
#SBATCH --time=0-06:00


module load StdEnv/2023 gcc/12.3
module load r/4.3
module load r-bundle-bioconductor/3.18
module load scipy-stack/2023b
module load haploview/4.2


Rscript ~/projects/rrg-girardsi/ARSACS/scripts/find_ibd_length_in_hgs_FFF.r