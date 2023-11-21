#!/bin/bash

# Automatically creates .fsf files for running 3rd-level GLM on group-level based on template; vary cluster-forming threshold.
#
# Make executable:
# chmod a+x feat_create_GLM_Sample_threshold.sh # make executable
#
# Execute:
# ./feat_create_GLM_Sample_threshold.sh # run
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

rootDir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

GLMID=1
sign="pos"


# Create directory if it doesn't exist yet 
if (($sign == "pos")); then
	signVal=1
else
	signVal=-1
fi
echo "SignVal is $signVal"

# set subject ID, loop
for cope in 4 5; do

	echo "Start cope $cope"

	# Make new directory if non-existent:
	mkdir -p ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/

	for thresh in {11..51}; do
			
		floatThresh=${thresh:0:-1}.${thresh: -1} # add digit

		echo "Start threshold $floatThresh"

		# Specify template, fill in COPE, SIGN, sign multiplier (VAL), THRESH (string), FLOAT (with dot)
		cat ${rootDir}/Analyses/fMRI_Scripts/FEAT_Templates/feat_template_GLM${GLMID}_sample.fsf | sed s/COPE/$cope/g | sed s/SIGN/$sign/g | sed s/VAL/$signVal/g | sed s/THRESH/$thresh/g | sed s/FLOAT/$floatThresh/g > ${rootDir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_sample_${sign}_cope${cope}_thresh${thresh}.fsf

	done

done # end subject loop
