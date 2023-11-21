#!/bin/bash

# Automatically creates .fsf files for running 1st-level block-wise GLM based on template
#
# Make executable:
# chmod a+x feat_create_GLM_Block.sh # make executable
#
# Execute:
# ./feat_create_GLM_Block.sh # run
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

GLMID=1

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

# Create directory if it doesn't exist yet 
mkdir -p ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Block_Scripts

# Set subject ID, loop
for (( subject=1 ; subject<=36 ; subject++ )); do

	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "Start subject ${subjectID}"

	# Set block ID, loop
	for (( block=1 ; block<=6 ; block++ ));	do
		blockID=`zeropad $block 1` # block ID with 1 digit

		echo "Start subject ${subjectID} block ${blockID}"

		# Get number of volumes for this block: use fslinfo, first five rows, last row of these, print number of volumes to variable npts
		npts=`fslinfo ${rootDir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/AROMA/denoised_func_data_nonaggr.nii.gz | head -n 5 | tail -n 1 | awk '{print $2}'`

		# Specify template, number of time points, number of blocks, output variable
		cat ${rootDir}/Analyses/fMRI_Scripts/FEAT_Templates/feat_template_GLM${GLMID}_block.fsf | sed s/TIMEPOINTS/$npts/g | sed s/SUBJECT/$subjectID/g | sed s/BLOCK/$blockID/g > ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Block_Scripts/feat_GLM${GLMID}_Block_Sub${subjectID}_Bl${blockID}.fsf

	done # end block loop

done # end subject loop
