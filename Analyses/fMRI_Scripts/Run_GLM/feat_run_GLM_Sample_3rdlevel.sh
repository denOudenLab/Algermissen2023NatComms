#!/bin/bash

# Run .fsf file across subjects per GLM (3nd level GLM in 3-level approach).
#
# Make executable:
# chmod a+x feat_run_GLM_Sample_3rdlevel.sh
# 
# Execute:
# ./feat_run_GLM_Sample_3rdlevel.sh # run
#use fsl_sub_DCCN.sh, run on -T, minutes, storage space, feat command, file to run (with complete path)
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

GLMID=1

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure
submitDir=${rootDir}/Analyses/fMRI_Scripts/Run_GLM # where fsl_sub_DCCN.sh is

for cope in 1 2 3 4 5; do
#for (( cope=1 ; cope<=5 ; cope++ )); do # loop throught these copes

	copeID=`zeropad $cope 1`; # 1 digit cope ID

	# 1) Run contrasts with cluster correction, z > 3.1, p < 0.05
	$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_Cope${copeID}_pos" -T 120 -u 10gb feat /project/3017042.02/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_sample_pos_cope${copeID}.fsf;
	$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_Cope${copeID}_neg" -T 120 -u 10gb feat /project/3017042.02/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_sample_neg_cope${copeID}.fsf;

	# 2) Run with excluding 29 subjects:
#	$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_Cope${copeID}_pos" -T 120 -u 10gb feat /project/3017042.02/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}without7_sample_pos_cope${copeID}.fsf;
#	$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_Cope${copeID}_pos" -T 120 -u 10gb feat /project/3017042.02/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}without7_sample_neg_cope${copeID}.fsf;

	# 3) Run with loop over thresholds:
#	for thresh in {11..51}; do
#		$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_Cope${copeID}_pos_${thresh}" -T 120 -u 10gb feat /project/3017042.02/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_sample_pos_cope${copeID}_thresh${thresh}.fsf;
#		$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_Cope${copeID}_neg_${thresh}" -T 120 -u 10gb feat /project/3017042.02/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_sample_neg_cope${copeID}_thresh${thresh}.fsf;

#	done
done
