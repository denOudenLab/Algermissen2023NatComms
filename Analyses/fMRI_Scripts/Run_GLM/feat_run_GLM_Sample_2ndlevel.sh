#!/bin/bash

# Runs second-level GLM one sample level based on .fsf file (see respective "FEAT_GLM${GLMID}_Sample_Scripts" folder)
#
# Make executable:
# chmod a+x feat_run_GLM_Sample_2ndlevel.sh
#
# ./feat_run_GLM_Sample_2ndlevel.sh # run
#use fsl_sub_DCCN.sh (in my home directory), run on -T, minutes, storage space, feat command, file to run (with complete path)
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure
submitDir=${rootDir}/Analyses/fMRI_Scripts/Run_GLM # where fsl_sub_DCCN.sh is

GLMID=1

# NOTE: Because this is second-level, you can fit all contrasts in one single command. No need to fit separate GLMs per contrast

# Regular GLMs:
$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_pos" -T 240 -u 10gb feat ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_Sample_Scripts_pos.fsf;
$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_neg" -T 240 -u 10gb feat ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_Sample_Scripts_neg.fsf;

# -----------------------------------------------------------------------------------
# Control with only 29 subjects (included in TAfT):
$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_pos" -T 240 -u 10gb feat ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_without7_Sample_Scripts_pos.fsf;
$submitDir/fsl_sub_DCCN.sh -N "GLM${GLMID}_neg" -T 240 -u 10gb feat ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_without7_Sample_Scripts_neg.fsf;

# END
