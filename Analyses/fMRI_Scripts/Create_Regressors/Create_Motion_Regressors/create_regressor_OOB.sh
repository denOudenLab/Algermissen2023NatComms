#!/bin/bash

# Extract volume-by-volume out-of-brain (OOB) signal as suggested by Maarten Mennes.
# Saved as OOB_noise.txt in AROMA folder of specific block of specific subject.
# 
# Make executable:
# chmod a+x create_regressor_OOB.sh
# 
# Submit as job to the cluster:
# qsub -N "create_regressor_OOB" -l walltime=2:00:00,mem=2gb create_regressor_OOB.sh
# Much faster than CSF: takes ~ 20 seconds per block, takes ~ 1.8 GB
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
fslDir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory structure

for (( subject=1 ; subject<=36; subject++ )); do 

	subjectID=`zeropad $subject 3`;
	subDir=${rootDir}/Log/fMRI/sub-${subjectID}

	for (( block=1 ; block<=6 ; block++ )); do 
		blockID=`zeropad $block 1`; 

		echo "Start Subject ${subjectID} Block ${blockID}"
		blockDir=${subDir}/FEAT_Block${blockID}.feat

		# Create new temporary folder to-be-deleted at the end:
		mkdir -p ${blockDir}/Nuisance_OOB

		# 1) Create out-of-brain mask based on mean functional image:
		# Maarten: "just take the mean image, it does not have to fit perfectly with the brain. Rather, you want to make sure there is absolutely no brain tissue that might overlap with your mask. That is why I would include even more dilation, so that you are definitely away from the brain's edge" 
		# a) Take mean_func as example of functional scan (already brain-extracted),		
		# b) Apply mean dilation of non-zero voxels several (here: 3) times (to make brain bigger/ out-of-brain smaller, i.e. been on the safe side to include no brain), 
		# c) Binarize (so voxels within brain becomes 1 and out-of-brain becomes 0);
		# d) invert mask (mul -1, add 1), so brain becomes 0 and out-of-brain becomes 1; 
		# "you could leave/add -dilM to make the mask larger/smaller according to your liking. I found one -dilM to be too close to the brain"		
		$fslDir/fslmaths ${blockDir}/mean_func.nii.gz -dilM -dilM -dilM -bin -mul -1 -add 1 ${blockDir}/Nuisance_OOB/OOBmask.nii.gz

		# 2) Extract time series of mean OOB values:
		# Input raw nifti data (not brain-extracted yet) and mask, get OOB_noise.txt as output: one column with mean OOB activity per volume
		$fslDir/fslmeants -i ${subDir}/NIFTI/sub-${subjectID}_block${blockID} -o ${blockDir}/Nuisance_OOB/OOB_noise.txt -m ${blockDir}/Nuisance_OOB/OOBmask.nii.gz

		# 3) Copy output into AROMA folder (usable across GLMs):
		cp ${blockDir}/Nuisance_OOB/OOB_noise.txt ${blockDir}/AROMA/OOB_noise.txt

		# 4) Optional: remove folder again to save space:
		rm -rf ${blockDir}/Nuisance_OOB
		echo "Finished Subject ${subjectID} Block ${blockID}"
	done;
done
