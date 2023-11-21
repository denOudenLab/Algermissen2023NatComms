#!/bin/bash

# Run .fsf file per block (1st level GLM in 3-level approach).
# 
# Make executable:
# chmod a+x feat_run_GLM_Block_1stlevel.sh 
#
# Execute:
# ./feat_run_GLM_Block_1stlevel.sh
# running GLM with 4 contrasts on each block takes ~ 20 minutes and ~ 5 GB
# running GLM with 6 contrasts on each block takes ~ 100 minutes and ~ 5 GB
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure
submitDir=${rootDir}/Analyses/fMRI_Scripts/Run_GLM # where fsl_sub_DCCN.sh is

GLMID=1

#for subject in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32 33 34 35 36; do
for (( subject=1 ; subject<=36 ; subject++ )); do # loop through these subjects

	subjectID=`zeropad $subject 3`; # 3 digit subject ID
#	subjectID=001

	for (( block=1 ; block<=6 ; block++ )); do # loop through these blocks

		blockID=`zeropad $block 1`; # 1 digit block ID

		# use fsl_sub_DCCN.sh (in my home directory), run on -T, minutes, storage space, feat command, file to run (with complete path)
		$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_Sub${subjectID}_Bl${blockID}" -T 150 -u 10gb feat /project/3017042.02/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Block_Scripts/feat_GLM${GLMID}_Block_Sub${subjectID}_Bl${blockID}.fsf;

	done 
done
