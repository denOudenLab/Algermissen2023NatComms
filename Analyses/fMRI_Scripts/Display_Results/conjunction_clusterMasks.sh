#!/bin/bash

# Multiply two cluster masks to see overlap (like conjunction)
# use fslmaths with option -mul
# uses cluster_mask_zstat1.nii.gz which is automatically created by FSL
# binarize at the end, so everything > 1 becomes 1

# chmod a+x conjunction_clusterMasks.sh # make executable
# ./conjunction_clusterMasks.sh # run

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

## Settings:
GLMID=1
signID=pos
firstCope=4
secondCope=5
newCope=6 # ${secondCope+1}

## Directories:
firstCopeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined/cope${firstCope}_${signID}.gfeat/cope1.feat
secondCopeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined/cope${secondCope}_${signID}.gfeat/cope1.feat
targetDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_all

echo "Combine for GLM $GLMID: Cope $firstCope and $secondCope"

# Store in second folder:
# thresh_zstat1
fslmaths ${firstCopeDir}/thresh_zstat1.nii.gz -mul ${secondCopeDir}/thresh_zstat1.nii.gz -bin ${secondCopeDir}/thresh_zstat1_${firstCope}_${signID}_${secondCope}.nii.gz

# cluster_mask_zstat1
newName=cluster_mask_${firstCope}_${signID}_${secondCope}_${signID}.nii.gz
fslmaths ${firstCopeDir}/cluster_mask_zstat1.nii.gz -mul ${secondCopeDir}/cluster_mask_zstat1.nii.gz -bin ${secondCopeDir}/${newName}

echo "New file ${newName} lies in ${secondCopeDir}/"

# Copy and unzip:
cp ${secondCopeDir}/cluster_mask_${firstCope}_${signID}_${secondCope}_${signID}.nii.gz ${targetDir}/thresh_zstat1_cope${newCope}_all.nii.gz
gunzip ${targetDir}/thresh_zstat1_cope${newCope}_all.nii.gz
echo "Copied as cope1_cope${newCope}_all.nii to ${targetDir}"

# Copy underlying cope and thres:
cp ${targetDir}/cope1_cope${firstCope}_pos.nii ${targetDir}/cope1_cope${newCope}_pos.nii
cp ${targetDir}/tstat1_cope${firstCope}_pos.nii ${targetDir}/tstat1_cope${newCope}_pos.nii
cp ${targetDir}/zstat1_cope${firstCope}_pos.nii ${targetDir}/zstat1_cope${newCope}_pos.nii

echo "Done :-)"

