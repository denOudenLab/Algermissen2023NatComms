#!/bin/bash

# Loop over subjects, extract summary statistics of parameter estimates within given mask using featquery.
# Steps involved:
# 1) Run invert_highres2standard_warp_subject.sh first to invert FNIRT registration matrix.
# 2) Brings given mask from standard to functional space for each subject.
# 3) Uses featquery to get summary statistics from parameter estimates within given mask for each subject.
# 
# Make executable:
# chmod a+x run_ROI_PEs_2ndlevel.sh # make executable
#
# Execute:
# ./run_ROI_PEs_2ndlevel.sh # run
# takes about 7-8 minutes per subject, so 4-5 h in total
# qsub -N "HippocampusValence" -l walltime=5:00:00,mem=10gb /project/3017042.02/Analyses/fMRI_Scripts/Run_ROI/run_ROI_PEs_2ndlevel.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

# Follow Andrew Jahn's video under: https://www.youtube.com/watch?v=SNVt7smHLm8
# Perform on single-subject level, thus loop over subjects
# Extract mean parameter estimates and analyse in R or Matlab or whatever

# Options:
# -p convert PE/COPE values into % 
# -s create time-series plots
# -b popup results in browser when finished (rather disable)

# Settings:

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

${fslDir}=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory

# GLM to extract from:
GLMIDextract=2

# a) Mask from 2nd-level GLM:
GLMIDmask=2

# Select mask for which to extract mean PEs: 
maskID=${GLMIDmask}StriatumValence
#maskID=${GLMIDmask}vmPFCValence
#maskID=${GLMIDmask}ACCAction

## File of mask:
maskFile=${rootDir}/Log/fMRI/fMRI_Masks/masksTAfTOutcomeLocked/${maskID}.nii.gz # functional masks from outcome-locked GLMs
#maskFile=${rootDir}/Log/fMRI/fMRI_Masks/masksHarvardOxford/${maskID}_bin.nii.gz # anatomical masks

# set subject ID, loop
for (( subject=1 ; subject<=36; subject++ )); do
# for subject in 36; do

	# Extract for each subject; mind: registration to MNI space already applied at first level; so no inverse registration needed
	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "Subject ${subjectID}: Start"

	subDir=${rootDir}/Log/fMRI/sub-${subjectID}/GLM${GLMIDextract}/GLM${GLMIDextract}_sub${subjectID}.feat

	# 1) Invert FNIRT registration matrix per subject:
	# run invert_highres2standard_warp_subject.sh first !!!!!

	# 2a) Bring mask from standard into functional space:
	echo "Subject ${subjectID}: Bring mask from standard to functinoal space"
	${fslDir}/applywarp -i ${maskFile} -r ${subDir}/example_func -o ${subDir}/${maskID}_func -w ${subDir}/reg/standard2highres_warp --postmat=${subDir}/reg/highres2example_func.mat

	# 3) Apply functional-space mask to PEs: (specify ho many PEs wanted!!!!!)
	echo "Subject ${subjectID}: Extract PEs"

	# 3a) First 4 PEs (GLM1): 
#	${fslDir}/featquery 1 ${subDir}/ 4 stats/pe1 stats/pe2 stats/pe3 stats/pe4 featquery_${maskID} -a 9 -p ${subDir}/${maskID}_func

	# 3b) 4 PEs later (GLM2):
#	${fslDir}/featquery 1 ${subDir}/ 4 stats/pe14 stats/pe15 stats/pe16 stats/pe17 featquery_${maskID} -a 9 -p ${subDir}/${maskID}_func

	# 3c) 8 PEs (GLM8E):
	${fslDir}/featquery 1 ${subDir}/ 8 stats/pe1 stats/pe2 stats/pe3 stats/pe4 stats/pe5 stats/pe6 stats/pe7 stats/pe8 featquery_${maskID} -a 9 -p ${subDir}/${maskID}_func

	echo "Subject ${subjectID}: Finished"

done # end subject loop
