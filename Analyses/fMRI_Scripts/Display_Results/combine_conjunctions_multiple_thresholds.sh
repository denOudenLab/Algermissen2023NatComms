#!/bin/bash

# For FEAT GLM outputs computed at different cluster-forming thresholds:
# Threshold 2 copes at given threshold,
# combine both cluster

# chmod a+x combine_conjunctions_multiple_thresholds.sh # make executable
# ./combine_conjunctions_multiple_thresholds.sh # run

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

GLMID=1 # 1without7
firstCope=4
secondCope=5
signID=pos

startThresh=11
stopThresh=51
clusterThresh=31

# ------------------------------------------------------------------------- #
# Loop over different thresholds:
# for signID in pos neg; do

iCount=0 # initialize count

baseDir=$rootdir/Log/fMRI/GLM${GLMID}_FEAT_Combined

# ------------------------------------------------------------------------- #
# Delete any old files:

rm -f $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii
rm -f $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}_bin$startThresh.nii
rm -f $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}_bin$stopThresh.nii


# ------------------------------------------------------------------------- #
#for threshID in ((21..31)); do
for (( threshID=${startThresh} ; threshID<=${stopThresh} ; threshID++ )); do # 

	echo "Start threshold $threshID"

	# Define folders:
	firstCopeDir=${baseDir}/cope${firstCope}_${signID}_thresh${threshID}.gfeat/cope1.feat
	secondCopeDir=${baseDir}/cope${secondCope}_${signID}_thresh${threshID}.gfeat/cope1.feat

	# Determine floatThresh:
	# https://stackoverflow.com/questions/12722095/how-do-i-use-floating-point-division-in-bash/35402635#35402635

	floatThresh=${threshID:0:-1}.${threshID: -1} # till last but 1 entry, dot, last entry

	echo "GLM $GLMID 1st cope $firstCope: Threshold at $floatThresh" # echo "$(($threshID/10))"
	fslmaths $firstCopeDir/thresh_zstat1.nii.gz -thr $floatThresh -bin $firstCopeDir/thresh_zstat1_bin.nii.gz
	
	echo "GLM $GLMID 2nd cope $secondCope: Threshold at $floatThresh"
	fslmaths $secondCopeDir/thresh_zstat1.nii.gz -thr $floatThresh -bin $secondCopeDir/thresh_zstat1_bin.nii.gz

	echo "GLM $GLMID: Combine copes $firstCope and $secondCope"
	fslmaths $firstCopeDir/thresh_zstat1_bin.nii.gz -mul $secondCopeDir/thresh_zstat1_bin.nii.gz $secondCopeDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz

	if (($threshID == $startThresh)); then
	# ------------------------------------------------------------------------- #
	# Initialize map:	
		echo "Initialize map at threshID $threshID"
		fslmaths $secondCopeDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz -mul $floatThresh $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz

	else
	# ------------------------------------------------------------------------- #
	# Add map at next threshold:
		echo "Add threshID $threshID"

		# 1) Multiply current map with binarized map of next threshold to create mask with values of current map:
		echo "1) Create mask"
		fslmaths $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz -mul $secondCopeDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz $baseDir/mask.nii.gz 

		# 2) Subtract this mask (so space where to add map of next threshold becomes zero): 
		echo "2) Subtract mask"
		fslmaths $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz -sub $baseDir/mask.nii.gz $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz

		# 3) Muliply map of next threshold with threshold itself:
		echo "3) Multiply with threshold"
		fslmaths $secondCopeDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz -mul $floatThresh $baseDir/tmp.nii.gz

		# 4) Add map of next threshold back to current map (fill space previous set to zero):
		echo "3) Add back to map"
	  	fslmaths $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz -add $baseDir/tmp.nii.gz $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz # add to previous
	fi

done

# ----------------------------------------------------------------------------------------- #
# Threshold:

# Lowest treshold for opacity:
floatThresh=${startThresh:0:-1}.${startThresh: -1}
echo "Threshold final map at $floatThresh"
fslmaths $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz -thr $floatThresh -bin $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}_bin$startThresh.nii.gz

# Cluster-forming treshold for black edges for significant clusters:
floatThresh=${clusterThresh:0:-1}.${clusterThresh: -1}
echo "Threshold final map at $floatThresh"
fslmaths $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz -thr $floatThresh -bin $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}_bin$clusterThresh.nii.gz

# ----------------------------------------------------------------------------------------- #
# Unzip:

echo "Unzip thresholded maps"
gunzip -f $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}.nii.gz
gunzip -f $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}_bin$startThresh.nii.gz
gunzip -f $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}_bin$stopThresh.nii.gz
gunzip -f $baseDir/thresh_zstat1_copes_${firstCope}_${secondCope}_${signID}_bin$clusterThresh.nii.gz

# Delete unnecessary files:
rm $baseDir/tmp.nii.gz
rm $baseDir/mask.nii.gz

echo "Finished :-)"

# END
