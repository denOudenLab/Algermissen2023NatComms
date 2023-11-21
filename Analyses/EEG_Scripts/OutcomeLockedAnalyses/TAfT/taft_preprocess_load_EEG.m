function [Y, job,betas] = taft_preprocess_load_EEG(job)

% [Y, job,betas] = taft_preprocess_load_EEG(job)
%
% Loads EEG data, reformat into matrix Y. 
% Also update job for which (good) trials are kept.
%
% INPUTS:
% job           = structure with settings to be specified via taft_preprocess_initialize_job
% .EEGdomain 	= string, either 'TF' (data in TF domain) or 'time' (data in time domain).
% .EEGfile      = string, full file path of subject-specific EEG data set to load.
% .selTrialsFile= string, name of file containing indices of selected trials to load. 
% .subID        = numeric scalar, subject ID.
%  
% OUTPUTS:
% Y             = a 3D (trials x channels x time) or 4D (trials x channels x frequency x time) matrix
% job           = structure with settings for creating TAfT object, new field added:
% .goodTrlIdxEEG= numeric vector, indices of trials to include; to be applied to EEG data.
% .goodTrlIdxfMRI= numeric vector, indices of trials to include; to be applied to fMRI data.
% betas         = trial-averaged Fieldtrip object, used as template to insert beta-map later while forwarding info on timing, channels, etc.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT/

fprintf('Load EEG data\n');

% ----------------------------------------------------------------------- %
%% Load EEG file:

% Load file:
inputEEG    = load(job.EEGfile);

% Extract:
subEEG      = inputEEG.(char(job.EEGname));
    
% Remove electrode positions if existing field:
if isfield(subEEG, 'elec')
    subEEG = rmfield(subEEG, 'elec');
end

% ----------------------------------------------------------------------- %
%% Retrieve good trials:

% job.goodTrlIdx      = subEEG.trialinfo; % invalid feedback trials not epoched; counts till 640 - #invalid

% ---------------------------------- %
% 1) Load trial indices for valid outcomes:
behav               = taft_preprocess_load_behavior(job);
validTrialsfMRI     = behav.valid; % all behavioral trials, counts till 640
validTrialsEEG      = (1:length(validTrialsfMRI))'; % count till 640 - # invalid trials

% ---------------------------------- %
% 2) Load Boolean indices of trials rejected in EEG, compare:
tmp                 = load(job.rejTrialsfile); % invalid trials not epoched; counts till 640 - #invalid
rejectedTrials      = tmp.rejectedtrials;
if length(validTrialsEEG) ~= length(rejectedTrials)
    error('Number of valid trials and length of vector of rejected trials do not match'); 
end

% ---------------------------------- %
% 3) Apply rejection to valid trials:
% EEG:
job.goodTrlIdxEEG  = validTrialsEEG(~rejectedTrials);
%[subEEG.trialinfo job.goodTrlIdxEEG] % should be identical
% mean(subEEG.trialinfo == job.goodTrlIdxEEG) % should be 1

% Number:
if length(job.goodTrlIdxEEG) ~= length(subEEG.trialinfo)
    error('Number of valid non-rejected trials and trials in EEG object do not match'); 
end
% Identity:
if sum(job.goodTrlIdxEEG == subEEG.trialinfo) ~= length(job.goodTrlIdxEEG)
    warning('Trial indices of behavior and EEG do no match; likely some invalid trials omitted; rely on behavior'); 
end

% fMRI:
% rejectedTrials contains logical indices, i.e. although EEG and fMRI data
% gives different names (numbers) to trials, the relative positions of
% to-be-rejected trials in both vectors are the same
job.goodTrlIdxfMRI  = validTrialsfMRI(~rejectedTrials);
if length(job.goodTrlIdxEEG) ~= length(job.goodTrlIdxfMRI)
    error('Number of valid non-rejected trials in fMRI and EEG differ'); 
end

% ---------------------------------- %
% 4) Check if other trials provided for selection:
if isfield(job, 'selTrialsFile')
    
    % Load selected behavioral trials:
    fprintf('Found file indicating selected %s: Load and adjust good trials\n', job.selTrials);
    selTrl              = load(job.selTrialsFile); % load selected trials
    selTrlfMRI          = selTrl.selTrials; % extract
    
    % ------------------------------ %
    % fMRI:
    % trial numbers are indices on full vector of 1:640
    % So just intersection between good fMRI trials and selected trials
    job.goodTrlIdxfMRI  = intersect(job.goodTrlIdxfMRI, selTrlfMRI); % intersection with fMRI trials
    
    % ------------------------------ %
    % EEG:
    % Translate to EEG metric:
    selTrlLog           = ismember(validTrialsfMRI, selTrlfMRI); % convert to logical indices on valid trials
    selTrlEEG          = validTrialsEEG(selTrlLog); % apply to EEG valid trials
    job.goodTrlIdxEEG  = intersect(job.goodTrlIdxEEG, selTrlEEG); % intersection with EEG trials

    % ------------------------------ %
    % Check again:
    if length(job.goodTrlIdxEEG) ~= length(job.goodTrlIdxfMRI)
        error('Number of selected valid non-rejected trials in fMRI and EEG differ'); 
    end
    
    % Re-select EEG data:
    cfg             = [];
    cfg.trials      = find(ismember(subEEG.trialinfo, job.goodTrlIdxEEG));
    subEEG         = ft_selectdata(cfg, subEEG); % select time samples
    
end

% ----------------------------------------------------------------------- %
%% Select EEG data:

% Select time:
if job.toi ~= [-inf inf] 
    cfg             = [];
    cfg.latency     = job.toi;
    subEEG          = ft_selectdata(cfg, subEEG); % select time samples
end

% Select frequencies:
if isfield(subEEG, 'freq') && any(job.foi ~= [-inf inf])
    cfg             = [];
    cfg.frequency   = job.foi;
    subEEG          = ft_selectdata(cfg, subEEG); % select frequencies
end    

% ----------------------------------------------------------------------- %
%% Create 3- or 4-D EEG object (channels, samples, (frequencies), trials): 

fprintf('Reshape EEG data\n');
Y                       = []; % initialize empty object

if strcmp(job.EEGdomain, 'time')

    cfg                 = [];
    betas               = ft_timelockanalysis(cfg, subEEG); % average over trials to obtain settings for final betascome object
    Y                   = zeros(size(betas.label, 1), size(betas.time, 2), size(subEEG.trialinfo, 1)); % initialize

    for iTrial = 1:size(subEEG.trialinfo, 1) 
        Y(:, :, iTrial)         = subEEG.trial{iTrial}; % trial becomes last dimension
    end
    
elseif strcmp(job.EEGdomain, 'TF')

    cfg                 = [];
    betas               = ft_freqdescriptives(cfg, subEEG); % average over trials to obtain settings for final betascome object
    Y                   = zeros(size(betas.label, 1), size(betas.freq, 2), ...
        size(betas.time, 2), size(subEEG.trialinfo, 1)); % initialize

    for iTrial = 1:size(subEEG.trialinfo, 1) 
            Y(:, :, :, iTrial)  = squeeze(subEEG.powspctrm(iTrial, :, :, :)); % trial becomes last dimension
    end
    
else
    error('Unknown EEGdomain data format');
end

end % END OF FUNCTION.