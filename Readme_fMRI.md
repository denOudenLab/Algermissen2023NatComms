# Readme fMRI

Analysis scripts for analyzing fMRI data and some relevant intermediate files (concatenated pre-processed blocks, regressors per GLM, masks used for ROI analyses) and results (maps from 3 GLMs). 

For fMRI **raw data** itself, see the separate collection under https://doi.org/10.34973/pezs-pw62.

Mind that in most files, you have to adjust the **root directory** to your own folder structure before running them.

## Scripts ##

We provide analyses scripts under *rootDir/Analyses/fMRI_Scripts*. 

- Folder with **template** .fsf files are:
	- *rootDir/Analyses/fMRI_Scripts/FEAT_Templates* --- contains .fsf template files which allow to create individual .fsf files for (each block of) each subject.	
 	- *rootDir/Analyses/fMRI_Scripts/FEAT_GLM_Sample_Scripts* --- contains .fsf files for 2nd-level/ 3rd level GLMs.
- For **submitting .fsf jobs** on a high performancecomputing cluster, use Run_GLM/fsl_sub_DCCN.sh (see in respective .sh files).
- For **pre-processing**, run scripts in the following order:
	- Run *rootDir/Analyses/fMRI_Scripts/Create_FEAT/feat_create_preAroma.sh* to create preprocessing files for each block for each subject
	- Run *rootDir/Analyses/fMRI_Scripts/Run_Feat/feat_run_preAroma.sh* to run pre-processing files for each block for each subject
	- Run *rootDir/Analyses/fMRI_Scripts/Run_ICA_AROMA/ICA_Aroma_run.sh* to run ICA-AROMA for each block for each subject
	- Run *rootDir/Analyses/fMRI_Scripts/Run_Feat/concat_flirt_blocks.sh* to concatenate blocks for each subject.

- For **first-level (subject level GLMs)**, run scripts in the following order:
	- **Nuisance regressors**: Run scripts in *rootDir/Analyses/fMRI_Scripts/Create_Regressors/Create_Motion_Regressors* to:
		- *create_regressor_CSF.sh*: Create mean CSF signal per subject per block.
		- *create_regressor_OOB.sh*: Create mean OOB signal per subject per block.
		- *Combine_Confound_EVs_Subject.m*: Combine CSF and OOB per subject across blocks, add realignment parameters, add spike regressor when relative displacement above cutoff, add intercept per block.
	- **Task regressors**: use *rootDir/Analyses/fMRI_Scripts/Create_Regressors/Create_Motion_Regressors/EEGfMRIPav_GLMX.m* to create tasks regressor per GLM X for each subject:
		- *GLM1*: GLM with 4 action x cue valence regressors, response (hand), error, standard update term, difference to biased update term, invalid trials. First level is block. Reported in main text (Fig. 3B-C).
		- *GLM1B*: similar to GLM1, but with standard and difference update term from M8 (dual perseveration model).
		- *GLM1C*: similar to GLM1, but with standard and difference update term from M9 (neutral outcomes reinterpretation model).
		- *GLM2*: GLM with 8 action x outcome valence regressors (at the time of outcome), 2 response (hand) regressors (at the time of responses), error, outcome onset, invalid trials. First level is subject. Reported in the main text (Fig. 3A and 6C-D).
		- *GLM3*: Similar to GLM1, but with trial-by-trial EEG regressors (frontal theta, midfrontal beta, midfrontal alpha) added. Contains versions 3A, 3B, 3C (select at start of file). First level is subject. Reported in the main text (Fig. 5).
			- Requires you to first run *rootDir/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/EEGfMRIPav_OutcomeLocked_8_TF_grouplevel_test.m* with contrastType = *'Valence'* and contrastType = 'Action' to save  clusters as masks, and then *rootDir/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/EEGfMRIPav_OutcomeLocked_8_extract_singleTrial.m* on respective masks to create trial-by-trial summary metric as BOLD regressor. Masks are also provided in folder rootDir/Log/EEG/OutcomeLockedResults/TF_mask.
	- **Fitting first (and second-) -level GLM**: 
	    - Run *rootDir/Analyses/fMRI_Scripts/Create_GLM/feat_create_GLM_Block.sh* or *rootDir/Analyses/fMRI_Scripts/Create_GLM/feat_create_GLM_Subject_1stlevel.sh* to create .fsf file for 1st-level GLM for each subject.
	    - Run *rootDir/Analyses/fMRI_Scripts/Create_GLM/feat_create_GLM_Subject.sh* to create .fsf file for 2nd-level GLM for each subject.
	    - Run *rootDir/Analyses/fMRI_Scripts/Create_GLM/feat_set_empty_regressors_block.sh* or *rootDir/Analyses/fMRI_Scripts/Create_GLM/feat_set_empty_regressors_subject.sh* to insert empty regressor settings into GLM .fsf files.
	    - Run *rootDir/Analyses/fMRI_Scripts/Run_GLM/feat_run_GLM_Block.sh* or *rootDir/Analyses/fMRI_Scripts/Run_GLM/feat_run_GLM_Subject_1stlevel.sh* to run 1st-level GLM for each subject.
	    - Run *rootDir/Analyses/fMRI_Scripts/Run_GLM/feat_run_GLM_Subject.sh* to run 2nd-level GLM for each subject.
- For **highest level (group-level) GLMs**, run scripts in the following order:
	- To create files with different cluster-forming thresholds (for plotting of GLM1), run *rootDir/Analyses/fMRI_Scripts/Create_GLM/feat_create_GLM_Sample_threshold.sh*.
	- Run *rootDir/Analyses/fMRI_Scripts/Run_GLM/Run_GLM/feat_run_GLM_Sample_3rdlevel.sh* or *rootDir/Analyses/fMRI_Scripts/Run_GLM/Run_GLM/feat_run_GLM_Sample_2ndlevel.sh* to run group-level GLM across subjects.
- For creating **dual coded fMRI images** with the slice display toolbox, run scripts in the following order:
	- Run *rootDir/Analyses/fMRI_Scripts/Display_Results/gather_unzip_contrasts_3level.sh* or *rootDir/Analyses/fMRI_Scripts/Display_Results/gather_unzip_contrasts_2level.sh* to assemble necessary files from all contrasts from a given GLM.
	- Run *rootDir/Analyses/fMRI_Scripts/Display_Results/conjunction_clusterMasks.sh* once to create conjunction of relevant contrasts 4 and 5 in GLM1.
	- Run *rootDir/Analyses/fMRI_Scripts/Display_Results/combine_conjunctions_multiple_thresholds.sh* to combine maps from multiple thresholds used for Fig. 3C.
	- Adjust directories to your own directory structure/ add paths to relevant toolboxes in *rootDir/Analyses/fMRI_Scripts/Display_Results/slice_display_set_dirs.m*.
	- Run *rootDir/Analyses/fMRI_Scripts/Display_Results/EEGfMRIPav_sliceDisplay_run.m* to create plots for each GLM for relevant slices.
- For creating **bar plots** of mean parameter estimates per regressor in selected ROI, 
	- Masks are provided under *rootDir/Log/fMRI/fMRI_Masks/masksTAfTOutcomeLocked* or can be alternatively created via executing *rootDir/Analyses/fMRI_Scripts/Run_ROI/Extract_ROIs/extract_ROI_HOA_GLM2_SVC.sh*.
	- Run *rootDir/Analyses/fMRI_Scripts/Run_ROI/invert_highres2standard_warp_subject.sh* to invert the native-space-to-MNI registration per subject.
	- Run *rootDir/Analyses/fMRI_Scripts/Run_ROI/run_ROI_PEs_2ndlevel.sh* to extract summary statistics on PE within mask.
	- Run *rootDir/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs/preprocess_PEs_outcomeLocked.R* to extract, rearrange, and save mean PE per regressor.
	- Run *rootDir/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs/plot_PEs_ROI_outcome_loop.m* (which calls *plot_PEs_ROI_outcome.m)* to create bar plots of parameter estimates per regressor.
		- Mind to adjust root directory in both *rootDir/Analyses/fMRI_Scripts/plot_PEs_ROI_outcome_run.m* and *rootDir/Analyses/fMRI_Scripts/plot_PEs_ROI_outcome_outcomeSorted.m*.
- For creating files to be used as regressors in **fMRI-informed EEG analyses (TAfT)**, run scripts in the following order: 
	- Masks are provided under *rootDir/Log/fMRI/fMRI_Masks/masksTAfTOutcomeLocked* or can be alternatively created via executing *rootDir/Analyses/fMRI_Scripts/Run_ROI/Extract_ROIs/extract_ROI_HOA_GLM1_TAfT.sh*.
	- Run *rootDir/Analyses/fMRI_Scripts/Run_ROI/invert_highres2standard_warp_block.sh* to invert the native-space-to-MNI registration per block.
	- Run *rootDir/Analyses/fMRI_Scripts/Run_ROI/run_ROI_rawdata_TAfT.sh* to extract mean time series from selected ROIs for fMRI-informed EEG analyses (will get stored in AROMA folder of each block of each subject).

## Data

We provide the following data:
- Preprocessed data  (after FEAT preprocessing) in *rootDir/Log/fMRI/sub-XXX*:
    - *postAROMA.nii.gz* --- preprocessed data from blocks concatenated for use in GLM.
    - *SubXXX_Nvolumes.txt* --- number of volumes per block.
- Regressors for GLM in GLM1: for each block: 
    - timings_blockX --- folder with all necessary regressors for block-level GLM, including:
        - *emptyregressors.txt* --- vector indicating whether regressor in design matrix is empty (1) or not (0)
        - *tConfoundEVs.txt* --- matrix with confound EVs (1 column, 1 value per volume): 6 realignment parameters, CSF, OOB, motion spikes.
- Regressors for GLMs in GLM2, GLM3A, GLM3B, GLM3C: 
    - timings_regressors --- folder with all necessary regressors for subject-level GLM, including:
        - *emptyregressors.txt* --- vector indicating whether regressor in design matrix is empty (1) or not (0)
        - *tConfoundEVs.txt* --- matrix with confound EVs (1 column, 1 value per volume): 6 realignment parameters, CSF, OOB, motion spikes.
- 2nd-level GLM results in *rootDir/Log/fMRI/*
    - *GLM1*: standard update term and difference to biased update term, reported in main text (Fig. 3 B-C):
        - GLM1_FEAT_Combined: contains separate folders for each contrast (4-5), both positive and negative, whole-brain correction.
    - *GLM2*: separate regressors for each action and outcome valence, actions regressors both outcome-locked anf response-locked, reported in main text (Fig. 3A and 6C-D):
        - GLM2_FEAT_Combined_neg and GLM2_FEAT_Combined_pos: whole-brain correction.
    - *GLM3*: like GLM1 plus EEG power regrsessor added, reported in main text (Fig. 5E-G):
        - GLM3_FEAT_Combined_neg and GLM3_FEAT_Combined_pos: whole-brain correction.
- Masks used in rootDir/Log/fMRI/fMRI_Masks:
    - *masksHarvardOxford*: anatomical masks of striatum and ACC/JBL used for GLMs with small volume correction (SVC).
    - *masksTAfTOutcomeLocked*: conjunctions of anatomical & functional masks used for (a) extraction of mean parameter estimates in selected ROIs) and (b) extraction of mean volume-by-volume time series for regressors in fMRI-informed EEG analyses (TAfT).

## Special notes on fMRI data
- For subject 004, no fieldmap was collected due to time constraints.
- For subjects 015 and 025, (structural-to-standard space) registration failed (likely due to a problem with the T1), so they are excluded from all analyses involving fMRI.
- For subject 025, no DTI data was collected due to time constraints.
- For subject 032, DTI and T1 scan were acquired in a separate session.

END OF FILE.
