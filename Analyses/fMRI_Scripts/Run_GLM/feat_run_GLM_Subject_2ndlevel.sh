#!/bin/bash

# Run .fsf file per subject (2nd level GLM in 3-level approach).
#
# Make executable:
# chmod a+x feat_run_GLM_Subject_2ndlevel.sh # make executable
#
# Execute:
# ./feat_run_GLM_Subject_2ndlevel.sh # run
# running GLM with 4 contrasts on each block takes ~ 3 minutes and ~ 7 GB
# running GLM with 6 contrasts on each block takes ~ 4 minutes and ~ 1 GB
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure
submitDir=${rootDir}/Analyses/fMRI_Scripts/Run_GLM # where fsl_sub_DCCN.sh is

GLMID=1
#for subject in 26 27 28 29 30 31 32 33 34 35 36; do 
for (( subject=1 ; subject<=36 ; subject++ )); do # loop throught these subjects

	subjectID=`zeropad $subject 3`; # 3 digit subject ID

	# use fsl_sub_DCCN.sh, run on -T, minutes, storage space, feat command, file to run (with complete path)

	$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_Sub${subjectID}" -T 60 -u 1gb feat $rootDir/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Subject_Scripts/feat_GLM${GLMID}_Subject_Sub${subjectID}.fsf;

done
