#!/bin/bash

#SBATCH --job-name=ants_preprocessing
#SBATCH --output=ants_preprocessing.out
#SBATCH --error=ants_preprocessing.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --partition=sandyb
#SBATCH --mem-per-cpu=2000

# Note: at least 2000 mem-per-cpu (2 GB) is required
IN_PATH=/scratch/midway/elayden/NTD  		# insert full path for data folder
N_SUBJECTS=3 					# insert number of subjects
N_RUNS=7 					# insert number of runs
ANTSPATH=/software/ANTs-2.1-el6-x86_64/bin/ 	# insert full path to ANTs binaries on your Midway account
SEGMENT=1					# boolean declaring whether to segment normalized anatomical
MAX_JOBS=6					# declare max jobs to submit at once for parallel processing (should == nodes * ntasks-per-node)

module load ANTs
bash ants_preprocessing.sh -i $IN_PATH -n $N_SUBJECTS -r $N_RUNS -a $ANTSPATH -s $SEGMENT -j $MAX_JOBS

