function taft_preprocess_main()

% taft_preprocess_main()
% 
% Interactive script.
% Main script to create a new TAfT object by loading EEG data (TF or time
% domain), fMRI data (volume-by-volume, upsample, epoch into
% trial-by-trial), and behavioral data and perform multiple linear
% regression for each channel/ (frequency)/ time/ point across trials per
% subject. 
% Use this script to select fMRI and behavioral regressors as well as
% subset of trials; all other settings are to be set in
% taft_preprocess_initialize_job.m.
% 
% INPUTS (interactive script):
% ROIs2use      = cell, one or more string(s) specifying the file names
%               of fMRI ROIs data files to include in design matrix. For any volume-by-volume regressor (also nuisance regressor).
% behav2use       cell, one or more string(s) specifying the behavioral variables 
%               to include in design matrix (optional). For any trial-by-trial regressor.
% selTrials     = string, sub-selection of trials to use (optional).
%
% OUTPUTS:
% Will save cell with Fieldtrip object per subject to disk.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT/

%% General settings:

% 0) EEG domain:
EEGdomain = 'TF';
% EEGdomain = 'time';

% ----------------------------------------------------------------------- %
% 1) ROI to use:
% 7 regions in main paper:
ROIs2use = {'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1LeftMotorConj', ...
    'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1LeftITGConj', 'GLM1V1Conj'}; % --> goes into paper


if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end
% ----------------------------------------------------------------------- %
% 2) Behavioral regressors:
% behav2use   = {}; % none

behav2use   = {'Updatestd', 'Updatedif'}; % --> goes into paper
% behav2use   = {'Updatestd', 'Updatedif', 'fbrel'}; % Updates + main effect outcome valence

% ----------------------------------------------------------------------- %
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


% ----------------------------------------------------------------------- %
%% Initialize final output object:

nSub            = 36;
allSubs         = 1:nSub;
invalidSubs     = [];
% invalidSubs     = 24; % for NoGoNoRewTrials
% invalidSubs     = 30; % for NoGoNoPunTrials
validSubs       = setdiff(allSubs, invalidSubs);

% Initialize final output object:
betas = cell(1,nSub); % final output for each channel/time/frequency 

% ----------------------------------------------------------------------- %
%% Loop over subjects, create separate Fieldtrip object per subject, store in cell, save:

for iSub = validSubs % loop through ROIs
    
    %% Initialize job:
    job                         = taft_preprocess_initialize_job(EEGdomain, iSub, ROIs2use, behav2use, selTrials);
    
    %% Load and reshape MEEG data:
    [Y, job, betas{iSub}]       = taft_preprocess_load_EEG(job);
    
    %% Load, fit, combine, and reshape fMRI-beta data of different ROIs and blocks:

    fprintf('Subject %03d: Load fMRI data\n', iSub);
    X                           = taft_preprocess_load_fMRI(job);
    
    %% FitGLM for each channel/timing/ frequency combination:

    if strcmp(job.EEGdomain, 'time') % data in time-frequency domain
        betas{iSub}.avg         = taft_preprocess_combine_EEG_fMRI(job, X, Y);
    elseif strcmp(job.EEGdomain, 'TF') % data in time-frequency domain
        betas{iSub}.powspctrm   = taft_preprocess_combine_EEG_fMRI(job, X, Y);
    end
    fprintf('Subject %03d: Finished\n', iSub);
    clear X Y
    
end

% ----------------------------------------------------------------------- %
%% Fill up for invalid subjects:

if ~ isempty(invalidSubs)
    betas{invalidSubs}.powspctrm = nan(size(betas{validSubs(1)}.powspctrm));
    betas{invalidSubs}.label     = betas{validSubs(1)}.label;
end

% ----------------------------------------------------------------------- %
%% Save:

fprintf('Save data as %s\n', job.outputFile);
save(job.outputFile, 'betas', 'job');
fprintf('Finished :-)\n');

end % END OF FUNCTION.