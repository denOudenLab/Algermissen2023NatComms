# README behavior

Analysis scripts for computational modeling of behavior using the CBM toolbox (https://github.com/payampiray/cbm).

For behavioral **raw data** itself, see the separate collection under https://doi.org/10.34973/pezs-pw62.

Note to set the **root directory** in respective script to your own folder structure before running them.

## Scripts
Analyses scripts are provided under *rootDir/Analyses/Behavior_Scripts*.
- Steps to recreate behavioral results:
	- *EEGfMRIPav_cbm_prepareData.m*: Run to convert .mat files under *rootDir/Log/Behavior/Data_beh_mat* into one central .mat input file under *rootDir/Log/Behavior/Modelling_Stan/CBM_Results*.
	- *EEGfMRIPav_cbm_fit.m*: Fit computational models, some simple evaluations.
	- *sim_subj.m*: Run once to create inputs for one example subject (stimuli, required actions, feedback validity) for model simulations.
	- *EEGfMRIPav_cbm_sim.m*: Run model simulations or one-step-ahead predictions for given model based on given parameter type.
	- *loop_sim.m*: Loop that specifies inputs to and calls *EEGfMRIPav_cbm_sim*.

## Data
Results from fitting and simulations are under *rootDir/Log/Behavior/Modelling_CBM/*:
    - *HBI_Results*: results from fitting using hierarchical Bayesian inference. Each file name starts with "hbi_mod" and then comprises the model identifiers included in the fitting process.
    - *LAP_Results*: results from fitting using LaPlace approximation. Each file name starts with "hbi_mod" and then comprises the model identifiers.
    - *Simulations*: results from 1000 model simulations either using parameters fitted via hierarchical Bayesian inference (hbi) or LaPlace approximation (lap).
    - *EEGfMRIPav_cbm_inputData.mat*: raw data used for fitting. 

END OF FILE.
