#!/bin/bash

#SBATCH --job-name=ants_preprocessing
#SBATCH --output=ants_preprocessing.out
#SBATCH --error=ants_preprocessing.err
#SBATCH --nodes=1
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000

# Note: at least 2000 mem-per-cpu (2 GB) is required
IN_PATH=/project2/bermanm/NTDants		# insert full path for data folder
N_SUBJECTS=46 					# insert number of subjects
N_RUNS=4 					# insert number of runs
SEGMENT=1					# boolean declaring whether to segment normalized anatomical
MAX_JOBS=16					# declare max jobs to submit at once for parallel processing (should == nodes * ntasks-per-node)
ANTSPATH=/software/ANTs-2.1-el6-x86_64/bin/ 	# insert full path to ANTs binaries

bash ants_preprocessing_single.sh -i $IN_PATH -n $N_SUBJECTS -r $N_RUNS -a $ANTSPATH -s $SEGMENT -j $MAX_JOBS

