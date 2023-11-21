#!/bin/bash

# Create an ROI based on a region from HarvardOxford Atlas and a group-level result (cope).
#
# Make executable:
# chmod a+x extract_ROI_HOA_GLM2_SVC.sh # make executable
#
# Execute:
# ./extract_ROI_HOA_GLM2_SVC.sh # just run locally, even just single jobs 
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

# Inventory cortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-cort-maxprob-thr25-2mm.xml
# Inventory subcortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-sub-maxprob-thr25-1mm.xml

# Set directories:
rootDir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
hoMaskDir=$rootDir/Log/fMRI/fMRI_Masks/masksHarvardOxford

fslDir=/opt/fsl/6.0.0 # FSL's directory--needs to be adapted to users' own directory
fslBinDir=$fslDir/bin
fslDataDir=$fslDir/data

pThreshold=10 # extract masks based on probabilities > 10%

GLMID=2 # GLM from which to extract group-level masks

# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
## I) Extract anatomical masks:

echo "Extract, binarize, and combine cortical regions"

# 1) ACC:
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz ${hoMaskDir}/CingulateAnterior 28 1
${fslBinDir}/fslmaths ${hoMaskDir}/CingulateAnterior -thr pThreshold -bin ${hoMaskDir}/CingulateAnterior_bin

# ---------------------------------------------------------------------- #
# 2) vmPFC:

# a) Extract all subparts:
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz ${hoMaskDir}/FrontalPole 0 1
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz ${hoMaskDir}/FrontalMedial 24 1
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz ${hoMaskDir}/Paracingulate 27 1

# b) Binarize:
${fslBinDir}/fslmaths ${hoMaskDir}/FrontalPole -thr pThreshold -bin ${hoMaskDir}/FrontalPole_bin
${fslBinDir}/fslmaths ${hoMaskDir}/FrontalMedial -thr pThreshold -bin ${hoMaskDir}/FrontalMedial_bin
${fslBinDir}/fslmaths ${hoMaskDir}/Paracingulate -thr pThreshold -bin ${hoMaskDir}/Paracingulate_bin

# c) Combine:
${fslBinDir}/fslmaths ${hoMaskDir}/FrontalPole_bin -add ${hoMaskDir}/FrontalMedial_bin -add ${hoMaskDir}/Paracingulate_bin ${hoMaskDir}/vmPFC

# d) Binarize again:
${fslBinDir}/fslmaths ${hoMaskDir}/vmPFC -bin ${hoMaskDir}/vmPFC_bin

# ---------------------------------------------------------------------- #
# 3) Striatum:

# a) Extract all subparts:
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz ${hoMaskDir}/AccumbensLeft 10 1
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz ${hoMaskDir}/AccumbensRight 20 1
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz ${hoMaskDir}/CaudateLeft 4 1
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz ${hoMaskDir}/CaudateRight 15 1
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz ${hoMaskDir}/PutamenLeft 5 1
${fslBinDir}/fslroi ${fslDataDir}/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz ${hoMaskDir}/PutamenRight 16 1

# b) Binarize:
${fslBinDir}/fslmaths ${hoMaskDir}/AccumbensLeft -thr pThreshold -bin ${hoMaskDir}/AccumbensLeft_bin
${fslBinDir}/fslmaths ${hoMaskDir}/AccumbensRight -thr pThreshold -bin ${hoMaskDir}/AccumbensRight_bin
${fslBinDir}/fslmaths ${hoMaskDir}/CaudateLeft -thr pThreshold -bin ${hoMaskDir}/CaudateLeft_bin
${fslBinDir}/fslmaths ${hoMaskDir}/CaudateRight -thr pThreshold -bin ${hoMaskDir}/CaudateRight_bin
${fslBinDir}/fslmaths ${hoMaskDir}/PutamenLeft -thr pThreshold -bin ${hoMaskDir}/PutamenLeft_bin
${fslBinDir}/fslmaths ${hoMaskDir}/PutamenRight -thr pThreshold -bin ${hoMaskDir}/PutamenRight_bin

# c) Combine:
${fslBinDir}/fslmaths ${hoMaskDir}/AccumbensLeft_bin -add ${hoMaskDir}/AccumbensRight_bin -add ${hoMaskDir}/CaudateLeft_bin -add ${hoMaskDir}/CaudateRight_bin -add ${hoMaskDir}/PutamenLeft_bin -add ${hoMaskDir}/PutamenRight_bin ${hoMaskDir}/Striatum
${fslBinDir}/fslmaths ${hoMaskDir}/AccumbensLeft_bin -add ${hoMaskDir}/AccumbensRight_bin ${hoMaskDir}/Accumbens

# d) Binarize again:
${fslBinDir}/fslmaths ${hoMaskDir}/Striatum -bin ${hoMaskDir}/Striatum_bin 
${fslBinDir}/fslmaths ${hoMaskDir}/Accumbens -bin ${hoMaskDir}/Accumbens_bin

# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
# II) Binarize group-level z-maps:

echo "Binarize group-level masks"

# ----------------------------- #
## Contrast 1: Action
copeID=1
signID=neg
copeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_${signID}.gfeat/cope${copeID}.feat
echo "Binarize functional map for contrast $copeID $signID"

${fslBinDir}/fslmaths ${copeDir}/thresh_zstat1 -thr 3.1 -bin ${copeDir}/thresh_zstat1_bin

# ----------------------------- #
## Contrast 2: Valence
copeID=2
signID=pos
copeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_${signID}.gfeat/cope${copeID}.feat
echo "Binarize functional map for contrast $copeID $signID"

${fslBinDir}/fslmaths ${copeDir}/thresh_zstat1 -thr 3.1 -bin ${copeDir}/thresh_zstat1_bin

# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
# III) Multiply anatomical and group-level mask:

echo "Combine anatomical and group masks"

taftMaskDir=${rootDir}/Log/fMRI/fMRI_Masks/masksTAfTOutcomeLocked

# ----------------------------------------------------------- #
## Contrast 1: Action

copeID=1
signID=neg
copeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_${signID}.gfeat/cope${copeID}.feat
echo "Multiple anatomical x functional for contrast $copeID $signID"

# ACC:
echo "ACC"
${fslBinDir}/fslmaths ${copeDir}/thresh_zstat1_bin -mul ${hoMaskDir}/ACC_bin -bin ${copeDir}/GLM${GLMID}ACCAction.nii.gz
cp ${copeDir}/GLM${GLMID}ACCAction.nii.gz ${taftMaskDir}/GLM${GLMID}ACCAction.nii.gz

# ----------------------------------------------------------- #
## Contrast 2: Valence

copeID=2
signID=pos
copeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_${signID}.gfeat/cope${copeID}.feat
echo "Multiple anatomical x functional for contrast $copeID $signID"

# Striatum:
echo "Striatum"
${fslBinDir}/fslmaths ${copeDir}/thresh_zstat1_bin -mul ${hoMaskDir}/Striatum_bin -bin ${copeDir}/GLM${GLMID}StriatumValence.nii.gz
cp ${copeDir}/GLM${GLMID}StriatumValence.nii.gz ${taftMaskDir}/GLM${GLMID}StriatumValence.nii.gz

# vmPFC:
echo "vmPFC"
${fslBinDir}/fslmaths ${copeDir}/thresh_zstat1_bin -mul ${hoMaskDir}/vmPFC_bin -bin ${copeDir}/GLM${GLMID}vmPFCValence.nii.gz
cp ${copeDir}/GLM${GLMID}vmPFCValence.nii.gz ${taftMaskDir}/GLM${GLMID}vmPFCValence.nii.gz

# END
