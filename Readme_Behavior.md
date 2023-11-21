# README behavior

Analysis scripts for analyzing behavioral data and some model outputs.

**Code** for the entire paper will be **maintained** under https://github.com/johalgermissen/Algermissen2021CerCor, with a permanent copy of the code at the time of publication under https://github.com/denoudenlab/Algermissen2021CerCor.

For behavioral **raw data** itself, see the separate collection under https://doi.org/10.34973/pezs-pw62.

Note to set the **root directory** in respective script to your own folder structure before running them.

## Scripts
Analyses scripts are provided under *rootDir/Analyses/Behavior_Scripts*.
- Steps to recreate behavioral results:
	- *EEGfMRIPav_extract_rawdata.m*: Run to convert .mat files under *rootDir/Log/Behavior/Data_beh_mat* into .csv files under rootDir/Log/Behavior/Data_beh_csv
	- *EEGfMRIPav_1_Read_Initial.R*: Read all csv files, recode variables, concatenate subjects and store them as EEGfMRIPav_all.csv under Log/Behavior/Behavior_concatenated/EEGfMRIPav_all.csv
	- *EEGfMRIPav_2_Preprocess_Automated.R*: When opening R, refresh variable coding (e.g. ordered factors) and create more convenience variables.
	- *EEGfMRIPav_3_MixedModels.R*: Fit mixed-effects models testing effects of cue valence and required action on choices and RTs.
	- *EEGfMRIPav_4_StayBehavior.R*: Fit mixed-effects models testing effects of outcome valence, outcome salience, and performed action on stay behavior.
	- *EEGfMRIPav_5_StayBehavior_Neural.R*: Fit mixed-effects models testing how trial-by-trial BOLD amplitude and EEG power predict stay behavior.
	- *EEGfMRIPav_6_EEG_Behavior_outcomelocked.R*: Fit mixed-effects models testing whether trial-by-trial EEG power changes over time and whether effects of behavioral variables interact with time.

## Data
Linear mixed-effects models fitted and reported in this paper as *.rds* files under *rootDir/Log/Behavior/Behavior_Models/*:
    - *mod_pGo_cor_valence_action.rds*: lme4 model with proportion Go responses as function of cue valence and required action.
    - *LRT_mod_pGo_cor_valence_action.rds*: respective model fitted with afex likelihood ration tests.
    - *mod_RT_nna_valence_action.rds*: lme4 model with reaction times (RTs) as function of cue valence and required action.
    - *LRT_mod_RT_nna_valence_action.rds*: respective model fitted with afex likelihood ration tests.
    - *mod_stay_action_salience_valence.rds: lme4 model with proportion stay behavior as a function of outcome valence, outcome salience, and performed action. 
	- also fitted for certain subsets of data.
    - *LRT_mod_stay_action_salience_valence.rds*: respective model fitted with afex likelihood ration tests.
    -  Every model exists once more with *\_withoutInvalid* and the end, which is the model fitted on only those 30 people that were included in the fMRI-informed EEG analyses (TAfT) reported in Supplementary Material S01.

END OF FILE.
