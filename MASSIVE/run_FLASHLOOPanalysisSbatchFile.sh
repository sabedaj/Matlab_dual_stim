#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=sabrina.meikle@monash.edu
#SBATCH --job-name=Loop_dataAnalysis
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --array=1-7
#SBATCH --mem-per-cpu=30000
module load matlab
matlab -nodisplay -nojvm -nosplash < loopSubdirectMASSIVE_flash.m
