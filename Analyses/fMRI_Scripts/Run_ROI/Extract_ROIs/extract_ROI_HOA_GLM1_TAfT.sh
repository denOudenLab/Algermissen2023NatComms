#!/bin/bash

# Create ROIs as conjunctions of an anatomical mask (from Harvard-Oxford Atlas; potentially combine multiple smaller regions to one bigger one) and selected clusters above zThreshold in GLM1.
# Used for creating volume-by-volume time-series regressors for TAfT and for extract mean PE estimates from selected regions.
#
# Involves the following steps:
# 1) Extract masks from NIFTI of Harvard-Oxford Atlas provided by FSL
# 2) Binarize masks based on probabiliy pThreshold
# 3) Combine smaller masks into bigger masks (like striatum, vmPFC)
# 4) Binarize again
# 5) Binarize pThresholded z-map from GLM1 (select respective contrasts; also do for GLM with SVC)
# 6) Multiply anatomical masks from HOA with pThresholded z-maps of respective contrast
#
# Make executable:
# chmod a+x extract_ROI_HOA_GLM1_TAfT.sh # make executable
#
# Execute:
# ./extract_ROI_HOA_GLM1_TAfT.sh # just run locally, even just single jobs 
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

# Inventory cortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-cort-maxprob-thr25-2mm.xml
# Inventory subcortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-sub-maxprob-thr25-1mm.xml

# Set directories:
rootDir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
anatDir=$rootDir/Log/fMRI/fMRI_Masks/masksHarvardOxford
targetDir=$rootDir/Log/fMRI/fMRI_Masks/masksTAfTCueLocked

fslDir=/opt/fsl/6.0.0 # FSL's directory--needs to be adapted to users' own directory structure
fslBinDir=$fslDir/bin
fslDataDir=$fslDir/data

pThreshold=10 # extract masks based on probabilities > 10%
zThreshold=3.1 # binarize thresholded z-masks at 3.1

# ---------------------------------------------------------------------------------------------------------------- #
## A) Cortex: CingulateAnteriorAction, LeftMotorHand, RightMotorHand:

echo "Extract and cortical masks"

## 1) Extract anatomical masks:
# ACC:
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatDir/CingulateAnterior 28 1
# JBL (used for mean PE extraction only):
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatDir/JuxtapositionalLobule 25 1
# Motor cortex:
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatDir/PrecentralGyrus 6 1
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatDir/PostcentralGyrus 16 1
# vmPFC:
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatDir/FrontalPole 0 1
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatDir/FrontalMedial 24 1
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatDir/Paracingulate 27 1

## 2) Binarize:
# ACC:
$fslBinDir/fslmaths $anatDir/CingulateAnterior -thr $pThreshold -bin $anatDir/CingulateAnterior_bin
# JBL:
$fslBinDir/fslmaths $anatDir/JuxtapositionalLobule -thr $pThreshold -bin $anatDir/JuxtapositionalLobule_bin
# Motor cortex:
$fslBinDir/fslmaths $anatDir/PrecentralGyrus -thr $pThreshold -bin $anatDir/PrecentralGyrus_bin
$fslBinDir/fslmaths $anatDir/PoscentralGyrus -thr $pThreshold -bin $anatDir/PostcentralGyrus_bin
# vmPFC:
$fslBinDir/fslmaths $anatDir/FrontalPole -thr $pThreshold -bin $anatDir/FrontalPole_bin
$fslBinDir/fslmaths $anatDir/FrontalMedial -thr $pThreshold -bin $anatDir/FrontalMedial_bin
$fslBinDir/fslmaths $anatDir/Paracingulate -thr $pThreshold -bin $anatDir/Paracingulate_bin

## 3) Combine to bigger masks:
# Motor cortex:
$fslBinDir/fslmaths $anatDir/PrecentralGyrus_bin -add $anatDir/PostcentralGyrus_bin $anatDir/MotorCortex
# vMPFC:
$fslBinDir/fslmaths $anatDir/FrontalPole_bin -add $anatDir/FrontalMedial_bin -add $anatDir/Paracingulate_bin $anatDir/vmPFC

## 4) Binarize again:
$fslBinDir/fslmaths $anatDir/MotorCortex -thr $pThreshold -bin $anatDir/MotorCortex_bin
$fslBinDir/fslmaths $anatDir/vmPFC -thr $pThreshold -bin $anatDir/vmPFC_bin

# ---------------------------------------------------------------------------------------------------------------- #
## B) Striatum:

echo "Extract and binarize striatal masks"

# 1) Extract all subparts:
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatDir/AccumbensLeft 10 1
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatDir/AccumbensRight 20 1
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatDir/CaudateLeft 4 1
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatDir/CaudateRight 15 1
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatDir/PutamenLeft 5 1
$fslBinDir/fslroi $fslDataDir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatDir/PutamenRight 16 1

# 2) Binarize:
$fslBinDir/fslmaths $anatDir/AccumbensLeft -thr $pThreshold -bin $anatDir/AccumbensLeft_bin
$fslBinDir/fslmaths $anatDir/AccumbensRight -thr $pThreshold -bin $anatDir/AccumbensRight_bin
$fslBinDir/fslmaths $anatDir/CaudateLeft -thr $pThreshold -bin $anatDir/CaudateLeft_bin
$fslBinDir/fslmaths $anatDir/CaudateRight -thr $pThreshold -bin $anatDir/CaudateRight_bin
$fslBinDir/fslmaths $anatDir/PutamenLeft -thr $pThreshold -bin $anatDir/PutamenLeft_bin
$fslBinDir/fslmaths $anatDir/PutamenRight -thr $pThreshold -bin $anatDir/PutamenRight_bin

# 3) Combine to entire bilateral striatum and accumbens masks:
$fslBinDir/fslmaths $anatDir/AccumbensLeft_bin -add $anatDir/AccumbensRight_bin -add $anatDir/CaudateLeft_bin -add $anatDir/CaudateRight_bin -add $anatDir/PutamenLeft_bin -add $anatDir/PutamenRight_bin $anatDir/Striatum
$fslBinDir/fslmaths $anatDir/AccumbensLeft_bin -add $anatDir/AccumbensRight_bin $anatDir/Accumbens

# 4) Binarize again:
$fslBinDir/fslmaths $anatDir/Striatum -bin $anatDir/Striatum_bin 
$fslBinDir/fslmaths $anatDir/Accumbens -bin $anatDir/Accumbens_bin

# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# 5) Binarize group-level z-maps from GLM1:
echo "Binarize group-level masks at $zThreshold"

# a) Valence positive: to be multiplied with vmPFC mask:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope1.feat/thresh_zstat1 -thr $zThreshold -bin ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope1.feat/thresh_zstat1_bin

# Also for GLM with SVC with striatal mask to get posterior putamen:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos_striatum.gfeat/cope1.feat/thresh_zstat1 -thr $zThreshold -bin ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos_striatum.gfeat/cope1.feat/thresh_zstat1_bin

# b) Valence negative: to be multiplied with ACC mask:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope1.feat/thresh_zstat1 -thr $zThreshold -bin ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope1.feat/thresh_zstat1_bin

# Also for GLM with SVC with striatal mask to get medial caudate:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_neg_striatum.gfeat/cope1.feat/thresh_zstat1 -thr $zThreshold -bin ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_neg_striatum.gfeat/cope1.feat/thresh_zstat1_bin

# c) Action positive: to be multiplied with striatum and ACC masks:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope2.feat/thresh_zstat1 -thr $zThreshold -bin ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope2.feat/thresh_zstat1_bin

# d) Conflict: to be multipled with JBL (only used for mean PE extract):
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GL1_FEAT_Combined_pos_ACC_JBL.gfeat/cope3.feat/thresh_zstat1 -thr $zThreshold -bin ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos_ACC_JBL.gfeat/cope3.feat/thresh_zstat1_bin

# e) Hand left: to be multiplied with right motor cortex mask:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope4.feat/thresh_zstat1 -thr $zThreshold -bin ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope4.feat/thresh_zstat1_bin

# f) Hand right: to be multiplied with left motor cortex mask:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope4.feat/thresh_zstat1 -thr $zThreshold -bin ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope4.feat/thresh_zstat1_bin

# ---------------------------------------------------------------------------------------------------------------- #
# 6) Multiply anatomical and group-level mask:
echo "Multiply anatomical and group masks"

# ----------------------------------------------------------------------------- #
## a) Valence contrast:

# vmPFCValence (adjust manually, see at the bottom of this script):
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope1.feat/thresh_zstat1_bin -mul $anatDir/vmPFC_bin -bin $targetDir/GLM1vmPFCValence.nii.gz

# ACCValence:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope1.feat/thresh_zstat1_bin -mul $anatDir/CingulateAnterior_bin -bin $targetDir/GLM1CingulateAnteriorValence.nii.gz

# LeftPutamenValence:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos_striatum.gfeat/cope1.feat/thresh_zstat1_bin -mul $anatDir/Striatum_bin -bin $targetDir/GLM1LeftPutamenValence.nii.gz

# MedialCaudateValence:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_neg_striatum.gfeat/cope1.feat/thresh_zstat1_bin -mul $anatDir/Striatum_bin -bin $targetDir/GLM1MedialCaudateValence.nii.gz

# ----------------------------------------------------------------------------- #
## b) Action contrast:

# StriatumAction:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope2.feat/thresh_zstat1_bin -mul $anatDir/Striatum_bin -bin $targetDir/GLM1StriatumAction.nii.gz

# ACCAction:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope2.feat/thresh_zstat1_bin -mul $anatDir/CingulateAnterior_bin -bin $targetDir/GLM1CingulateAnteriorAction.nii.gz

# ----------------------------------------------------------------------------- #
## c) Conflict contrast:

$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos_ACC_JBL.gfeat/cope3.feat/thresh_zstat1_bin -mul $anatDir/JuxtapositionalLobule_bin -bin $targetDir/GLM1JBLIncongruency.nii.gz

# ----------------------------------------------------------------------------- #
## d) Left/ right hand contrast:

# RightMotorHand:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope4.feat/thresh_zstat1_bin -mul $anatDir/MotorCortex_bin -bin $targetDir/GLM1RightMotorHand.nii.gz

# LeftMotorHand:
$fslBinDir/fslmaths ${rootDir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope4.feat/thresh_zstat1_bin -mul $anatDir/MotorCortex_bin -bin $targetDir/GLM1LeftMotorHand.nii.gz

echo "Finished :-)"

# e) Adjust manually in FSLeyes --> GLM1vmPFCValenceMan.nii.gz
# gunzip ${rootDir}/Log/fMRI/fMRI_Masks/masksTAfTCueLocked/GLM1vmPFCValenceMan.nii.gz
# END
