function taft_postprocess_time_main()

% taft_postprocess_time_main()
% 
% Interactive script.
% Main script to load a previously created TAfT object, align data
% (channels) across subjects, do one-sample-t-test across subjects and plot
% as ER (line-) plot and as topoplot. 
% Use this script to select fMRI and behavioral regressors as well as
% subset of trials; all other settings are to be set in
% taft_postprocess_load_job.m.
% 
% INPUTS (interactive script):
% EEGdomain 	= string, set to 'time' (data in time domain) in this script.
% ROIs2use      = cell, one or more string(s) specifying the file names
%       of fMRI ROIs data files to include in design matrix. For any volume-by-volume regressor (also nuisance regressor).
% behav2use       cell, one or more string(s) specifying the behavioral variables 
%       to include in design matrix (optional). For any trial-by-trial regressor.
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

EEGdomain   = 'time';

% ROIs2use    = {}; % empty
ROIs2use    = {'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1LeftMotorConj', 'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1LeftITGConj', 'GLM1V1Conj'}; % --> GOES INTO PAPER

if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end

% 2) Behavioral regressors:
% behav2use   = {}; % none
behav2use    = {'Updatestd', 'Updatedif'}; % --> GOES INTO PAPER
% behav2use    = {'Updatestd', 'Updatedif', 'fbrel'}; 

% 3) Perform only on selected trials
selTrials   = 'all';

% selTrials   = 'GoTrials';
% selTrials   = 'NoGoTrials';
% selTrials   = 'PositiveTrials';
% selTrials   = 'NegativeTrials';
% selTrials   = 'RewardTrials';
% selTrials   = 'PunishmentTrials';
% selTrials   = 'NeutralTrials';
% selTrials   = 'NoRewardTrials';
% selTrials   = 'NoPunishmentTrials';

% selTrials   = 'GoRewTrials';
% selTrials   = 'GoNoRewTrials';
% selTrials   = 'GoNoPunTrials';
% selTrials   = 'GoPunTrials';
% selTrials   = 'NoGoRewTrials';
% selTrials   = 'NoGoNoRewTrials'; % sub 024 only 1 trial
% selTrials   = 'NoGoNoPunTrials'; % sub 030 only 1 trial
% selTrials   = 'NoGoPunTrials';

[job, dirs, betas] = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
%% ERplot: One single ROI:

% job.invalidSubs = []; % for Updatestd/Updatedif
job.invalidSubs = [11 12 23 30]; % for Updatestd and Updatedif
% job.invalidSubs = [11 12 15 23 25 26 30]; % for fMRI regressors

rng(20190823) % set random number generator for constant p-values

iROI            = 1;
[sortBetas,~]   = taft_postprocess_time_selectData(job, betas, iROI);
sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
selChans        = {'Fz', 'FCz', 'Cz'}; % Jenn's a-priori selection
% selChans        = {'Oz', 'O1', 'O2', 'O3', 'O4', 'POz', 'PO1', 'PO2'}; % for Updatestd and Updatedif
% selChans        = {'O2', 'O4', 'PO4'}; % for Updatestd and Updatedif
% selChans        = {'FCz', 'Cz'}; % Jenn's a-priori selection--more focused

% Without saving:
[corrp,tg]      = taft_postprocess_time_ERplot(job,dirs, sortBetas, iROI, selChans,2.0, 1000,false);
pause(2)
close gcf

% ----------------------------------------------------------------------- %
%% Identify single clusters:

timeVec     = sortBetas{1}.time;
timeVec     = timeVec(timeVec >= 0 & timeVec <= 0.7);
timeVec     = timeVec(squeeze(corrp) < 0.05);
tVal        = squeeze(tg);
tVal        = tVal(squeeze(corrp) < 0.05);

fprintf('Identified clusters:\n')
startTime   = timeVec(1); lastTime = timeVec(1);
for iTime = 1:length(timeVec) % iTime = 2;
    thisTime = timeVec(iTime);
    if round(thisTime - lastTime,3) > .0010 || iTime == length(timeVec)
        stopTime = lastTime; % new end of cluster
        if tVal(iTime - 1) < 0 
            signID = 'negative'; 
        else
            signID = 'positive';
        end
        fprintf('Significant %s cluster from %.03f - %.03f sec.\n', signID, startTime, stopTime);
        startTime = thisTime; % new start of cluster
    end
    lastTime = timeVec(iTime); % overwrite timing of last bin
end

% ----------------------------------------------------------------------- %
%% ERplot: Loop over all ROIs:

dirs.TFplot = '/project/3017042.02/Log/OutcomeLockedPaperPlots';

% job.invalidSubs = []; % for Updatestd/Updatedif
job.invalidSubs = [11 12 23 30]; % for Updatestd and Updatedif
% job.invalidSubs = [11 12 15 23 25 26 30]; % for fMRI regressors

for iROI = 1:size(betas{1}.avg, 1)
    rng(20190822) % set random number generator for constant p-values
    [sortBetas, ~]  = taft_postprocess_time_selectData(job, betas, iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
    allChanROIs     = {{'Fz', 'FCz', 'Cz'}}; % Jenn's a-priori selection
    for iChan = 1:length(allChanROIs) % iChan = 1;
        selChans    = allChanROIs{iChan};
        [~, ~]      = taft_postprocess_time_ERplot(job, dirs, sortBetas, iROI, selChans, 2.0, 1000, true);
        pause(1)
        close gcf
    end
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% TOPOPLOT: Multiple time domain topoplot:

dirs.topoplot = '/project/3017042.02/Log/OutcomeLockedPaperPlots';

for iROI = 1:size(betas{1}.avg, 1)
    [~, Tvalues]    = taft_postprocess_time_selectData(job, betas, iROI);
    taft_postprocess_time_topoplot(job, dirs, Tvalues, iROI, 0, 0.70, 0.05, 3, false);
    pause(2)
    close gcf
end

end % END OF FUNCTION.