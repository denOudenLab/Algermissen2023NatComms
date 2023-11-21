#!/bin/bash

# Gather thresh_zstat1, cope1, tstat1, and zstat1 files of positive and negative group-level contrast of given GLM,
# unzip all, combine and unzip thresholded masks
# 
# Make executable:
# chmod a+x gather_unzip_contrasts_2level.sh
# 
# Execute:
# ./gather_unzip_contrasts_2level.sh # run
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

GLMID=2  # GLM name plus suffices, like "1_without7"
maxCope=6 # number of copes in this GLM

# Adjust number of copes in for-loop, both for retrieval and for thresholded-map combination!
# if additional postfixes: adjust pos and neg and start of signID for-loop below

# ------------------------------------------------------------------- #
# Define directory where to copy niftis too:
rootDir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
targetDir=$rootDir/Log/fMRI/GLM${GLMID}_FEAT_Combined_all # where files go to

fslDir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory structure

# ------------------------------------------------------------------- #
# Make directory if it doesn't exist yet:
mkdir -p $targetDir
echo "Target directory is $targetDir"

# 1) Copy all relevant files into one common directory:

for ((copeID=1; copeID<=$maxCope; copeID++)); do

	for signID in pos neg; do

		echo "copy cope $copeID $signID"
	
		copeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_${signID}.gfeat/cope${copeID}.feat

		# ------------------------------------------------------------------- #
		# Copy & rename thresholded map after contrast:
		cp $copeDir/thresh_zstat1.nii.gz $targetDir/thresh_zstat1_cope${copeID}_${signID}.nii.gz

		# ------------------------------------------------------------------- #
		# Copy & rename cope, tstat, zstat:
		cp $copeDir/stats/cope1.nii.gz $targetDir/cope1_cope${copeID}_${signID}.nii.gz
		cp $copeDir/stats/tstat1.nii.gz $targetDir/tstat1_cope${copeID}_${signID}.nii.gz
		cp $copeDir/stats/zstat1.nii.gz $targetDir/zstat1_cope${copeID}_${signID}.nii.gz

	done
done

# ------------------------------------------------------------------- #
# 2) Unzip all nii.gz files:

echo "Unzip all"

for file in $targetDir/*.nii.gz; do

	gunzip $file

done

# ------------------------------------------------------------------- #
# 3) Combine thresholded maps (and unzip again):

echo "Combine thresholded maps"
cd $targetDir # cd to output directory

for ((copeID=1; copeID<=$maxCope; copeID++)); do

	$fslDir/fslmaths $targetDir/thresh_zstat1_cope${copeID}_pos.nii -add $targetDir/thresh_zstat1_cope${copeID}_neg.nii $targetDir/thresh_zstat1_cope${copeID}_all.nii
	gunzip $targetDir/thresh_zstat1_cope${copeID}_all.nii

done 

echo "Finished :-)"
# END
