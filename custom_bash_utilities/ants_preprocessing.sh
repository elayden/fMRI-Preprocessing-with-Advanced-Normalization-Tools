#!/bin/bash

# ants_preprocessing.sh

function Usage {
    cat <<USAGE

Usage:

`basename $0` -i IN_PATH -n N_SUBJECTS -r N_RUNS -a ANTSPATH -s SEGMENT -j MAX_JOBS

Compulsory arguments:

     -i IN_PATH: 	Input folder containing data for all subjects; Note: cannot contain spaces

     -n N_SUBJECTS:    	Integer specifying number of subjects/anatomicals.

     -r N_RUNS:		Integer specifying number of functional runs

     -a ANTSPATH:       Full path to ANTs binaries; Note: cannot contain spaces

     -s SEGMENT:	boolean: segment normalized anatomicals?

     -j MAX_JOBS: 	Integer specifying the max number of jobs to run in parallel (nodes * cores)

--------------------------------------------------------------------------------------
ants_preprocessing.sh was created by:
--------------------------------------------------------------------------------------
Elliot Layden
Environmental Neuroscience Laboratory
The University of Chicago
--------------------------------------------------------------------------------------

USAGE
    exit 1
}

# Function for limiting # jobs submitted to specified MAX_JOBS:
function limitJobs {
	jobslist=($(jobs -p))
	while (( ${#jobslist[*]} >= $1 ))
	do
    		sleep 1
    		jobslist=($(jobs -p))
	done
}

################################
# Get Inputs:
################################

# Provide output for Help
if [[ "$1" == "-h" ]]; then
	Usage >&2
fi

# Read in command line arguments
while getopts "h:i:n:r:a:s:j:" OPT
  do
  case $OPT in
      h) #help
   echo "$USAGE"
   exit 0
   ;;
       i) # Input Folder
   	IN_PATH=$OPTARG
   ;;
       n) # N Subjects
      	N_SUBJECTS=$OPTARG
   ;;
       r) # N Runs
      	N_RUNS=$OPTARG
   ;;
       a)
   	ANTSPATH=$OPTARG
   ;;
       s)
   	SEGMENT=$OPTARG
   ;;
       j)
	MAX_JOBS=$OPTARG
   ;;
      \?) # getopts issues an error message
      echo "$USAGE" >&2
      exit 1
      ;;
  esac
done

####################################
# Perform some preliminary items:
####################################

# Change directory to path specified:
cd $IN_PATH

# Prefix for Motion Correction Files:
motion_out="moco"

# Set ANTS Environment Variable (required):
export ANTSPATH=$ANTSPATH  

# Leave One Core Available For Main Thread:
MAX_JOBS=$[$MAX_JOBS-1]

# Initialize subjects and trials arrays for full parallel:
N_COUNT=$[$N_SUBJECTS * $N_RUNS - 1] # minus 1 because this is used as the max for zero-based indexing
declare -a subjects 
declare -a runs 
count=0
if (( $N_SUBJECTS < 100 )); then
	for i in $(seq -f "%02g" 1 $N_SUBJECTS)
	do
		for j in $(seq -f "%02g" 1 $N_RUNS)
		do
			subjects[$count]=$i
			runs[$count]=$j
			count=$[$count + 1]
		done
	done
elif (( $N_SUBJECTS<1000 )); then
	for i in $(seq -f "%03g" 1 $N_SUBJECTS)
	do
		for j in $(seq -f "%02g" 1 $N_RUNS)
		do
			subjects[count]=$i
			runs[count]=$j
			count=$[$count + 1]
		done
	done
else
	for i in $(seq -f "%04g" 1 $N_SUBJECTS)
	do
		for j in $(seq -f "%02g" 1 $N_RUNS)
		do
			subjects[count]=$i
			runs[count]=$j
			count=$[$count + 1]
		done
	done
fi

# Check Arrays:
# echo ${subjects[*]}
# echo ${runs[*]}
# exit

####################################
# FUNCTIONAL BIAS CORRECTION
####################################
for iter in $(seq -f "%01g" 0 $N_COUNT)
do
	# Limit Number of Jobs Submitted to MAX_JOBS:
	limitJobs $MAX_JOBS

	i=${subjects[$iter]}
	j=${runs[$iter]}

	(in=sub-${i}_run-${j}_bold

	# Make sure file exists:
	if [ -f ${in}.nii ]; then
		$ANTSPATH/N4BiasFieldCorrection -d 4 -i ${in}.nii -o ${in}.nii
	else
		echo ${in}.nii "does not exist"
	fi
	       
	)&
done
wait

####################################
# STRUCTURAL BIAS CORRECTION
####################################
for i in $(seq -f "%02g" 1 $N_SUBJECTS)
do
	# Limit Number of Jobs Submitted to MAX_JOBS:
	limitJobs $MAX_JOBS

	# Make sure file exists
	(if [ ! -f sub-${i}_anat.nii ]; then
		echo sub-${i}_anat.nii "does not exist"
    		continue
	fi
	$ANTSPATH/N4BiasFieldCorrection -d 3 -i sub-${i}_anat.nii -o b_sub-${i}_anat.nii)&
done
wait

####################################
# FUNCTIONAL MOTION CORRECTION
####################################
for iter in $(seq -f "%01g" 0 $N_COUNT)
do
	# Limit Number of Jobs Submitted to MAX_JOBS:
	limitJobs $MAX_JOBS

	(
	i=${subjects[$iter]}
	j=${runs[$iter]}

	in=sub-${i}_run-${j}_bold

	# Make sure file exists
	if [ -f ${in}.nii ]; then

	# average the time series
	${ANTSPATH}/antsMotionCorr -d 3 -a ${in}.nii -o ${in}_avg.nii.gz
		    
	# do affine motion correction
	${ANTSPATH}/antsMotionCorr -d 3 -o [${motion_out}_${in}_, ${motion_out}_${in}_Warped.nii.gz, ${motion_out}_${in}_avg.nii.gz] -m MI[${in}_avg.nii.gz, ${in}.nii, 1, 32, Regular, 0.2] -t Rigid[0.05] -i 20x10 -u 1 -e 1 -s 1x0 -f 2x1 -n 10 -w 1
		
	# Compute Voxel-wise Displacements:
	${ANTSPATH}/antsMotionCorrStats -x ${in}_avg.nii.gz -m ${motion_out}_${in}_MOCOparams.csv -o voxelWiseDisplace_${in}.nii.gz -f 1 -s displaceMap_${in}.nii -d ${motion_out}_${in}_Warp.nii.gz

	else
		echo ${in}.nii "does not exist"
	fi
	)&
done
wait

#####################################
# COREGISTER AVG MOCO FUNCT & ANAT
#####################################
for iter in $(seq -f "%01g" 0 $N_COUNT)
do
	# Limit Number of Jobs Submitted to MAX_JOBS:
	limitJobs $MAX_JOBS

	(
	i=${subjects[$iter]}
	j=${runs[$iter]}

	${ANTSPATH}/antsRegistrationSyN.sh -d 3 -n 6 -f b_sub-${i}_anat.nii -m moco_sub-${i}_run-${j}_bold_avg.nii.gz -t a -o coreg_sub-${i}_run-${j}_
	)&
done
wait

####################################
# BUILD GROUP-WISE TEMPLATE
# note: initialize to MNI Template
# May take several hours
####################################
${ANTSPATH}/antsMultivariateTemplateConstruction.sh -d 3 -o groupWise_ -c 2 -i 4 -j $MAX_JOBS -r 0 -y 0 -z mni_temp.nii -t GR b_sub*
wait
# -r 1 (rigid-body registration of all inputs) is not supposed to help if an initial template -z is specified
# -y 0 prevents template drift

####################################
# Register groupWise to ICBM-152
# Fairly quick, ~30 min.
####################################
${ANTSPATH}/antsRegistrationSyN.sh -d 3 -f mni_temp.nii -m groupWise_template0.nii.gz -n $MAX_JOBS -t s -o groupWise_MNI_
wait

####################################
# APPLY WARPS TO STRUCTS IN 1 STEP
####################################
for i in $(seq -f "%02g" 1 $N_SUBJECTS)
do
	# Limit Number of Jobs Submitted to MAX_JOBS:
	limitJobs $MAX_JOBS

	(# Collapse all transforms to 1 displacement field
	${ANTSPATH}/antsApplyTransforms -d 3 -o [CollapsedWarp_${i}.nii.gz,1] \
		-t groupWise_MNI_1Warp.nii.gz \
		-t groupWise_MNI_0GenericAffine.mat \
		-t groupWise_b_sub-${i}_anat.nii$((10#$i-1))Warp.nii.gz \
		-t groupWise_b_sub-${i}_anat.nii$((10#$i-1))Affine.txt \
		-r groupWise_MNI_Warped.nii.gz

	# Apply displacement field
	${ANTSPATH}/antsApplyTransforms -d 3 -e 0 -n BSpline[3] \
		-i b_sub-${i}_anat.nii \
		-r groupWise_MNI_Warped.nii.gz \
		-t CollapsedWarp_${i}.nii.gz \
		-o groupWise_MNI_${i}_anat.nii

	# Delete displacement field to save space
	rm CollapsedWarp_${i}.nii.gz)&

done
wait

#####################################
# Register mni_temp.nii to rmni_temp.nii
####################################
${ANTSPATH}/antsRegistrationSyN.sh -d 3 -f rmni_temp.nii -m mni_temp.nii -t r -o downsamp_temp_
wait

#####################################
# 1. Collapse Funct2Template Transforms
# 2. Replicate these to 4D along with Template
# 3. Apply to Original Functs
#####################################
for iter in $(seq -f "%01g" 0 $N_COUNT)
do
	# Limit Number of Jobs Submitted to MAX_JOBS:
	limitJobs $MAX_JOBS

	(
	i=${subjects[$iter]}
	j=${runs[$iter]}

	# Get header information from funct:
	fmri=sub-${i}_run-${j}_bold.nii
	hislice=`${ANTSPATH}/PrintHeader $fmri | grep Dimens | cut -d ',' -f 4 | cut -d ']' -f 1`
	tr=`${ANTSPATH}/PrintHeader $fmri | grep "Voxel Spac" | cut -d ',' -f 4 | cut -d ']' -f 1`

	# 1. Collapse all transforms to 1 displacement field
	${ANTSPATH}/antsApplyTransforms -d 3 -o [CollapsedWarp_sub-${i}_run-${j}_bold.nii.gz,1] \
		-t downsamp_temp_0GenericAffine.mat \
		-t groupWise_MNI_1Warp.nii.gz \
		-t groupWise_MNI_0GenericAffine.mat \
		-t groupWise_b_sub-${i}_anat.nii$((10#$i-1))Warp.nii.gz \
		-t groupWise_b_sub-${i}_anat.nii$((10#$i-1))Affine.txt \
		-t coreg_sub-${i}_run-${j}_0GenericAffine.mat \
		-r rmni_temp.nii

	# dimension, output, ReplicateDisplacement VectorFieldName TimeDims TimeSpacing TimeOrigin
	${ANTSPATH}/ImageMath 3 rep_CollapsedWarp_sub-${i}_run-${j}_bold.nii.gz ReplicateDisplacement \
  		CollapsedWarp_sub-${i}_run-${j}_bold.nii.gz $hislice $tr 0

	${ANTSPATH}/ImageMath 3 rep_rmni_temp_${i}_${j}.nii.gz ReplicateImage \
       		rmni_temp.nii $hislice $tr 0
	)&
done
wait

for iter in $(seq -f "%01g" 0 $N_COUNT)
do
	# Limit Number of Jobs Submitted to MAX_JOBS:
	limitJobs $[$MAX_JOBS/2]  # limit to 1/2 available due to RAM usage

	(
	i=${subjects[$iter]}
	j=${runs[$iter]}
	fmri=sub-${i}_run-${j}_bold.nii

	# Get header information from funct:
	fmri=sub-${i}_run-${j}_bold.nii

	# Apply displacement field to original bold data
	${ANTSPATH}/antsApplyTransforms -d 4 -o w_$fmri \
  		-t rep_CollapsedWarp_sub-${i}_run-${j}_bold.nii.gz \
		-t ${motion_out}_sub-${i}_run-${j}_bold_Warp.nii.gz  \
  		-r rep_rmni_temp_${i}_${j}.nii.gz \
		-i sub-${i}_run-${j}_bold.nii

	# Delete displacement fields to save space
	rm rep_rmni_temp_${i}_${j}.nii.gz
       	rm rep_CollapsedWarp_sub-${i}_run-${j}_bold.nii.gz
	rm CollapsedWarp_sub-${i}_run-${j}_bold.nii.gz
	)&
done
wait

################################
# Segment normalized structs if requested
###############################
if [ $SEGMENT ]; then
    for i in $(seq -f "%02g" 1 $N_SUBJECTS)
	do
	# Limit Number of Jobs Submitted to MAX_JOBS:
	limitJobs $MAX_JOBS

	(${ANTSPATH}/Atropos -d 3 \
	     -a groupWise_MNI_${i}_anat.nii \
	     -i PriorProbabilityImages[6, mni_tissue_%02d.nii, .6] \
	     -x mni_mask.nii \
	     -o [labelledTissues_${i}.nii,pTissue_%02d_s_${i}.nii])&
	done
	wait
fi
