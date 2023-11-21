function taft_postprocess_TF_main()

% taft_postprocess_TF_main()
% 
% Interactive script.
% Main script to load a previously created TAfT object, align data
% (channels) across subjects, do one-sample-t-test across subjects and plot
% as TF plot and as topoplot. 
% Use this script to select fMRI and behavioral regressors as well as
% subset of trials; all other settings are to be set in
% taft_postprocess_load_job.m.
% 
% INPUTS (interactive script):
% EEGdomain 	= string, set to 'TF' (data in TF domain) in this script.
% ROIs2use      = cell, one or more string(s) specifying the file names
%               of fMRI ROIs data files to include in design matrix. For any volume-by-volume regressor (also nuisance regressor).
% behav2use       cell, one or more string(s) specifying the behavioral variables 
%               to include in design matrix (optional). For any trial-by-trial regressor.
% selTrials     = string, sub-selection of trials to use (optional).
%
% OUTPUTS:
% Will print results to console and create plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT/

% ----------------------------------------------------------------------- %
%% Load EEG-fMRI regression weights, specify directories:

EEGdomain           = 'TF';

% ROIs2use            = {}; % empty
ROIs2use            = {'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1LeftMotorConj', ...
    'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1LeftITGConj', 'GLM1V1Conj'}; % --> GOES INTO PAPER

if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end

% 2) Behavioral regressors:
% behav2use           = {}; % none
behav2use           = {'Updatestd', 'Updatedif'}; % --> GOES INTO PAPER
% behav2use           = {'Updatestd', 'Updatedif', 'fbrel'}; % Updates + outcome sign

% 3) Perform only on selected trials:
selTrials           = 'all';

[job, dirs, betas]  = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
%% TF plot: One single ROI:

% job.invalidSubs = []; % all people
% job.invalidSubs = [11 12 23 30]; % for Updatestd and Updatedif
job.invalidSubs = [11 12 15 23 25 26 30]; % for fMRI regressors --> TAfT outliers

rng(20190822) % set random number generator for constant p-values

iROI            = 1; % iROI to be tested/ plotted

[sortBetas,~]   = taft_postprocess_TF_selectData(job, betas, iROI);
sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:

selChans        = {'Fz', 'FCz', 'Cz'}; % Jenn's a-priori selection
[~,~]           = taft_postprocess_TF_TFplot(job, dirs, sortBetas, ...
    iROI, selChans, 2, 1000, 2, false); % job, dirs, sortBetas, iROI, selChans, thresh,nP,zlim, isSave
pause(2)
close gcf

% ----------------------------------------------------------------------- %
%% TF plot: Loop over all ROIs:

dirs.TFplot = '/project/3017042.02/Log/OutcomeLockedPaperPlots';

% job.invalidSubs = []; % all people
% job.invalidSubs = [11 12 23 30]; % for Updatestd and Updatedif
job.invalidSubs = [11 12 15 23 25 26 30]; % for fMRI regressors

for iROI = 1:size(betas{1}.powspctrm, 1)
% for iROI = 1:2

    rng(20190822); % set seed for reproducibility of permutation test.

    [sortBetas,~]   = taft_postprocess_TF_selectData(job, betas, iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
    allChanROIs     = {{'Fz', 'FCz', 'Cz'}};

    for iChan = 1:length(allChanROIs) % iROI = 1;
        selChans    = allChanROIs{iChan};
        job.yScale  = 'log';
%         job.yScale  = 'lin';
        [~,~]  = taft_postprocess_TF_TFplot(job, dirs, sortBetas, ...
            iROI, selChans, 2, 10000, 2, true); % iROI, selChans, thresh,nP
%         pause(2)        
        close gcf
        pause(3);   
    end
end

% ----------------------------------------------------------------------- %
%% Topoplot: for all ROIs, for all frequency bands:

% job.invalidSubs = []; % for Updatestd/Updatedif
% job.invalidSubs = [11 12 23 30]; % for Updatestd and Updatedif
job.invalidSubs = [11 12 15 23 25 26 30]; % for fMRI regressors

dirs.topoplot = '/project/3017042.02/Log/OutcomeLockedPaperPlots';

% Select ROI:
for iROI = 1:size(betas{1}.powspctrm, 1) % iROI = 1;
    [~, Tvalues]        = taft_postprocess_TF_selectData(job, betas, iROI);
    taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, ...
        [1 4 8 13], [4 8 13 30], 0.0, 0.9, 0.1, 3, false); % ROI, startFreq, endFreq
    pause(2)
%     close all
end
fprintf('Finished :-) \n')


% ----------------------------------------------------------------------- %
%% Topoplot: one ROI, separate command for each frequency band:

startTime = 0; endTime = 0.9; step = 0.1; % stimulus-locked
iROI = 1;
[~, Tvalues] = taft_postprocess_TF_selectData(job, betas, iROI);

% Each band:
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, 1,4, startTime, endTime, step, 3, false); % iROI, startFreq, endFreq
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, 4, 8, startTime, endTime, step, 3, false); % iROI, startFreq, endFreq
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, 8, 13, startTime, endTime, step, 3, false); % iROI, startFreq, endFreq
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, 13, 30, startTime, endTime, step, 3, false); % iROI, startFreq, endFreq
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, [1 4 8 13], [4 8 13 30], 0, 1.1, 0.1, 3, false); % ROI, startFreq, endFreq
pause(5)
close all

% Custom frequency bands:
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, 6, 10, startTime, endTime, step, 3, false); % low alpha
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, 14, 28, startTime, endTime, step, 3, false); % focused beta

% Single topoplot:
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, 7, 22, 0.100, 0.100, 0.800, 3, false); % striatum cluster
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, 1, 6, 0.175, 0.175, 0.325, 3, false); % PCC cluster
taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, 6, 16, 0.100, 0.100, 0.300, 3, false); % ACC cluster

% ----------------------------------------------------------------------- %
%% Topoplot: Loop through ROIs, all frequency bands:

for iROI = 1:7 % iROI = 7;

    [~, Tvalues] = taft_postprocess_TF_selectData(job, betas, iROI);
    taft_postprocess_TF_topoplot(job, dirs, Tvalues, iROI, ...
        [1 4 8 13], [4 8 13 30], 0, 1.1, 0.1, 3, false) % ROI, startFreq, endFreq, startTime, endTime, steps, zlim, isSave
    close all

end

end % END OF FUNCTION.