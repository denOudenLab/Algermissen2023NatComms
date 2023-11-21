#!/bin/bash

# Gather necessary files for running sliceDisplay for all contrasts of 1 GLM (with 3 levels)
# use thresh_zstat1 to create maps of significant clusters
# use cope1, tstat1, zstat1 as "background"
# chmod a+x gather_unzip_contrasts_3level.sh # make executable
# ./gather_unzip_contrasts_3level.sh # run

GLMID=1
maxCope=5
suffix= #

# Define directory where to copy niftis too:
rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure
targetDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_all${suffix}

# Make directory if it doesn't exist yet:
mkdir -p $targetDir
echo "Target directory is $targetDir"

# 1) Copy all relevant files into one common directory:

# for copeID in {1..10}; do
for copeID in {4..5}; do
# for ((copeID=1; copeID<=$maxCope; copeID++)); do # cope number

	for signID in pos neg; do # sign positive or negative

		copeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined${suffix}/cope${copeID}_${signID}.gfeat/cope1.feat
		echo "copy cope $copeID $signID"

		# Rename thresholded map after contrast:
		cp ${copeDir}/thresh_zstat1.nii.gz $targetDir/thresh_zstat1_cope${copeID}_${signID}${suffix}.nii.gz

		# CD to stats folder, copy and rename cope/tstat1/zstat1 after contrast:
		cp ${copeDir}/stats/cope1.nii.gz $targetDir/cope1_cope${copeID}_${signID}${suffix}.nii.gz
		cp ${copeDir}/stats/tstat1.nii.gz $targetDir/tstat1_cope${copeID}_${signID}${suffix}.nii.gz
		cp ${copeDir}/stats/zstat1.nii.gz $targetDir/zstat1_cope${copeID}_${signID}${suffix}.nii.gz

	done
done

# 2) Unzip all nii.gz files:

echo "Unzip all"
cd $targetDir # cd to output directory

for file in *.nii.gz; do

	gunzip $file

done

# 3) Combine thresholded maps:

echo "Combine thresholded maps"
cd $targetDir # cd to output directory

#for copeID in 5; do
for copeID in {4..5}; do
#for ((copeID=1; copeID<=$maxCope; copeID++)); do

	fslmaths thresh_zstat1_cope${copeID}_pos.nii -add thresh_zstat1_cope${copeID}_neg.nii thresh_zstat1_cope${copeID}_all${suffix}.nii
	gunzip thresh_zstat1_cope${copeID}_all${suffix}.nii

done 
# END
