#!/bin/bash

# Gather cope and zstat from first cope, combine 
# use thresh_zstat1 to create maps of significant clusters
# use cope1, tstat1, zstat1 as "background"
# chmod a+x gather_conjunction_3level.sh # make executable
# ./gather_conjunction_3level.sh # run

## Define root directory:
rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

# ------------------------------------------------------------------- #
## GLM and cope settings:
GLMID=1

maxCope=5
newCope=$(($maxCope + 1))

copeIDA=4
copeIDB=5

signID=pos

suffix= # empty

# ------------------------------------------------------------------- #
## Define directory where to copy niftis too:
targetDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_all${suffix}

# Make directory if it doesn't exist yet:
mkdir -p ${targetDir}
echo "Target directory is ${targetDir}"

# ------------------------------------------------------------------- #
# 1) Copy and rename cope and z_stat of first contrast:

copeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined${suffix}/cope${copeIDA}_${signID}.gfeat/cope1.feat
cp ${copeDir}/stats/cope1.nii.gz ${targetDir}/cope1_cope${newCope}_${signID}${suffix}.nii.gz
cp ${copeDir}/stats/tstat1.nii.gz ${targetDir}/tstat1_cope${newCope}_${signID}${suffix}.nii.gz
cp ${copeDir}/stats/zstat1.nii.gz ${targetDir}/zstat1_cope${newCope}_${signID}${suffix}.nii.gz

# ------------------------------------------------------------------- #
# 2) Unzip files:
echo "Unzip all"
cd ${targetDir} # cd to output directory
for file in *${newCope}_${signID}${suffix}.nii.gz; do
	gunzip $file
done

# ------------------------------------------------------------------- #
# 3) Copy thresholded masks:

echo "Copy files $copeIDA"
copeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined${suffix}/cope${copeIDA}_${signID}.gfeat/cope1.feat
cp ${copeDir}/thresh_zstat1.nii.gz ${targetDir}/thresh_zstat1_cope${newCope}A.nii.gz

echo "Copy files $copeIDB"
copeDir=${rootDir}/Log/fMRI/GLM${GLMID}_FEAT_Combined${suffix}/cope${copeIDB}_${signID}.gfeat/cope1.feat
cp ${copeDir}/thresh_zstat1.nii.gz ${targetDir}/thresh_zstat1_cope${newCope}B.nii.gz

# ------------------------------------------------------------------- #
# 4) Combine & unzip thresholded masks:
echo "Combine thresholded masks"
fslmaths ${targetDir}/thresh_zstat1_cope${newCope}A.nii.gz -mul ${targetDir}/thresh_zstat1_cope${newCope}B.nii.gz ${targetDir}/thresh_zstat1_cope${newCope}_all${suffix}.nii.gz
gunzip ${targetDir}/thresh_zstat1_cope${newCope}_all${suffix}.nii.gz

# Delete temporary files:
rm -rf ${targetDir}/thresh_zstat1_cope${newCope}A.nii.gz
rm -rf ${targetDir}/thresh_zstat1_cope${newCope}B.nii.gz

# END
