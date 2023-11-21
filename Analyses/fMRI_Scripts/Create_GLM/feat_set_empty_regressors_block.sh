#!/bin/bash

# Sets regressors in .fsf scripts to empty if specified as empty in emptyregressors.txt file in specific folder of respective GLM of respective block of respective subject.
# 
# Make executable:
# chmod a+x feat_set_empty_regressors_block.sh # make executable
#
# Execute:
# ./feat_set_empty_regressors_block.sh # run
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

GLMID=1
nRegressors=10

rootDir=/project/3017042.02

# Create directory if it doesn't exist yet 
mkdir -p ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Block_Scripts

# set subject ID, loop
for (( subject=1 ; subject<=36 ; subject++ )); do

	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "start subject ${subjectID}"

	# set block ID, loop
	for (( block=1 ; block<=6 ; block++ ));	do

		blockID=`zeropad $block 1` # block ID with 1 digit
		echo "start subject ${subjectID} block ${blockID}"

		file="${rootDir}/Log/fMRI/sub-${subjectID}/GLM${GLMID}/timings_block${blockID}/emptyregressors.txt" # set file name
		emptyregressors=$(cat "$file") # store in variable

		# Check for each regressor whether 1
		for (( regressorID=1 ; regressorID<=${nRegressors} ; regressorID++ )); do

			position=$(($regressorID  * 2 - 2)) # given the commas, we need positions 0, 2, 4, ...
			status=${emptyregressors:$position:1} # from this position, only 1 character
			# echo "Regressor No. ${regressorID} has value $status"

			if (($status == 1)); then
				# set regressor status from 3 to 10 (empty)
				echo "Regressor No. $regressorID is empty, set to status empty (10)"
				sed -i "s/set fmri(shape$regressorID) 3/set fmri(shape$regressorID) 10/g" ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Block_Scripts/feat_GLM${GLMID}_Block_Sub${subjectID}_Bl${blockID}.fsf
			fi
		done # regressor-loop
	done # block-loop
done # subject-loop
