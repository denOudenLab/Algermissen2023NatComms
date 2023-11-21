#!/bin/bash

# Automatically creates .fsf files for running 2nd-level subject-wise GLM based on GLM-specific template.
#
# Make executable:
# chmod a+x feat_create_GLM_Subject.sh # make executable
# 
# Execute:
# ./feat_create_GLM_Subject.sh # run
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

GLMID=1 # GLM number

# Create directory if it doesn't exist yet 
mkdir -p ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Subject_Scripts

# set subject ID, loop
for (( subject=1 ; subject<=36 ; subject++ )); do

	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "start subject $subjectID"

	# Specify template, subject number:
	cat ${rootDir}/Analyses/fMRI_Scripts/FEAT_Templates/feat_template_GLM${GLMID}_subject.fsf | sed s/SUBJECT/$subjectID/g  > ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Subject_Scripts/feat_GLM${GLMID}_Subject_Sub${subjectID}.fsf

done # end subject loop
