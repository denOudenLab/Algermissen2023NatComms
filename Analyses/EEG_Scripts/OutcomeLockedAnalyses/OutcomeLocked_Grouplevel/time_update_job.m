function [job] = time_update_job(job)

% [job] = time_update_job(job)
%
% Complete settings for job on time-domain data to create.
% (data sets to create; plotting settings).
% 
% INPUTS:
% job                   = cell, created via time_update_job.m, needs at 
% least fields:
%   .nSub               = integer, number of subjects.
%   .sub2exclude        = numeric vector, numbers of subjects to exclude.
%   baselineSettings  = string, type of baseline correction, 'trend'
%   (regression-based), 'ft' (with Fieldtrip's ft_timelockbaseline).
%   .baselineTimings    = vector of 1 or 2, timing of baseline correction
%   for which to load data.
%   .responseSettings   = string, type of response setting for which
%   conditions are split, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo).
%   .outcomeSettings    = string, type of outcome coding for which 
%   conditions are split, 'abs' (positive/neutral/negative), 'rel' (positive/negative), 'all' (reward/no reward/ no punishment/ punishment).
%   .chanArea           = string, area of channels to select, 'frontal',
%   'completefrontal', 'midfrontal', 'midoccipital', 'rightoccipital', 'occipital'.
%   .band               = string, frequency band to select, 'delta', 'theta', 'thetadelta', 'lowalpha', 'verylowalpha', 'alpha', 'beta', 'middlebeta', 'broad'.
%   .contrastType       = string, contrast to be used: 'Preferred', 'Action',
%   'GoPreferred', 'SalientPreferred', 'SalientAction'.
%   .outERPcor          = Boolean, read data corrected for outcome-locked 
% 	ERP (true) or not (false), default false.
%   .ROI2use 	    = string, split data by ROI (instead of by conditions:'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1V1Conj'), optional.
% 
% OUTPUTS:
% job                   = cell, with the following fields added:
%   .inputFileName      = name of file to load: 
%   .layout             = string, cap to use for topoplots.
%   .timetiming         = numeric vector of 2 elements, timing for TF
%   plots.
%   .sigTime            = numeric vector of 2 elements, timing for
%   topoplots.
%   .channels           = vector of strings, selected channel names.
%   .chanName           = string, name of selected channel area
%   capitalized.
%   .invalidSubs        = numeric vector, subjects to exclude.
%   .validSubs          = numeric vector, subjects to include.
%   .nValidSubs         = numeric, number of subjects to include.
%   .lockName           = string, name of event-locking, fully written out.
%   .nCond              = numeric, number of conditions.
%   .condOrder          = numeric vector, order of conditions to plot.
%   .colMat             = matrix, colors of lines to plot per condition.
%   .lineStyle          = vector of strings, line types to plot per
%   condition.
%   .condNames          = vector of strings, condition labels to use in
%   plots.
%   .plotName           = string, important settings to be included in plot
%   names.
%   .twoLineColor       = 2 x 3 matrix, line colors to use for twoLinePlot.
%   .twoLineLabels      = vector of 2 strings, labels for conditions.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

fprintf('Update and complete job settings\n')

% ----------------------------------------------------------------------- %
%% File to load:

% Default input file:
job.inputFileName       = sprintf('EEGfMRIPav_time_baseCor_%s_baseTime_%03dms_%03dms', ...
    job.baselineSettings, job.baselineTimings(1)*-1000, job.baselineTimings(end)*-1000);

% If split by BOLD:
if isfield(job, 'ROI2use')
    job.inputFileName   = sprintf('%s_BOLD_%s', job.inputFileName, strjoin(job.ROI2use, ''));
    job.contrastType    = 'BOLD';
else
    % Standard input file:
    job.inputFileName   = [job.inputFileName sprintf('_resp_%s_fb_%s', ...
        job.responseSettings, job.outcomeSettings)];
end

% Add final ending:
job.inputFileName       = sprintf('%s.mat', job.inputFileName); % add final ending

% ----------------------------------------------------------------------- %
%% Plotting options:

% Cap layout:
job.layout      = 'easycapM11.mat'; %% Timing of time plot:

% Timing for topoplot: retrieve significant time from job if available
if isfield(job, 'sigTime') % delete in case previously loaded
    job         = rmfield(job, 'sigTime');
end

% ----------------------------------------------------------------------- %
%% Time of significant differences:

if isfield(job, 'ROI2use') 

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels   = {'High BOLD', 'Low BOLD'};

    if strcmp(job.ROI2use, 'GLM1StriatumConj')
        job.sigTime     = [0.158 0.194 0.318 0.381]; % two clusters
        
    elseif strcmp(job.ROI2use, 'GLM1ACCConjMan')
        job.sigTime     = [];

    elseif strcmp(job.ROI2use, 'GLM1vmPFCConjMan')
        job.sigTime     = [0.350  0.413];
        
    elseif strcmp(job.ROI2use, 'GLM1PCCConj')
        job.sigTime     = [];
        
    elseif strcmp(job.ROI2use, 'GLM1V1Conj')
        job.sigTime     = [0.297 0.356];

    else
        warning('unknown sigTime')
        job.sigTime     = [];
    end
    
elseif strcmp(job.contrastType, 'Preferred')

    job.twoLineColor    = [0 0.6 0.2; 0.8 0 0]; % green, red
    job.twoLineLabels   = {'Positive', 'Negative'};
        
        if strcmp(job.outcomeSettings, 'rel')
            job.sigTime = [0.246 0.294 0.344 0.414];
            
        elseif strcmp(job.outcomeSettings, 'all')
            job.sigTime = [0.248 0.292 0.341 0.412 ];
        else
            warning('unknown sigTime')
            job.sigTime = [];
        end
        
elseif strcmp(job.contrastType, 'Action')

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels   = {'Go', 'NoGo'};

elseif strcmp(job.contrastType, 'GoPreferred')

    job.twoLineColor    = [0 0.6 0.2; 0.8 0 0]; % green, red
    job.twoLineLabels   = {'preferred Go', 'non-preferred Go'};

elseif strcmp(job.contrastType, 'SalientPreferred')

    job.twoLineColor    = [0 0.6 0.2; 0.8 0 0]; % green, red
    job.twoLineLabels   = {'Reward', 'Punishment'};

    job.sigTime = [0.127 0.181 0.313 0.424];

elseif strcmp(job.contrastType, 'SalientAction')

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels   = {'Salient Go', 'Salient NoGo'};

end

% ----------------------------------------------------------------------- %
%% Channels and labels:

if strcmp(job.chanArea, 'frontopolar')
    job.channels = {'Fpz', 'Fp1', 'Fp2'};
    job.chanName = 'Frontopolar';
    
elseif strcmp(job.chanArea, 'frontal')
    job.channels = {'Fpz', 'AFz', 'Fz'};
    job.chanName = 'Frontal';
    
elseif strcmp(job.chanArea, 'midfrontal')
    job.channels = {'Fz', 'FCz', 'Cz'};
    job.chanName = 'Midfrontal';
    
elseif strcmp(job.chanArea, 'midoccipital')
    job.channels = {'Pz', 'POz', 'Oz'};
    job.chanName = 'Midoccipital';
    
elseif strcmp(job.chanArea, 'rightoccipital')
    job.channels = {'P4', 'PO4', 'O2'};
    job.chanName = 'Right occipital';
    
elseif strcmp(job.chanArea, 'occipital')
    job.channels = {'Oz', 'O1', 'O2', 'POz', 'PO3', 'PO4'};
    job.chanName = 'Occipital';
    
else
    fprintf('No channels set\n')
end

% ----------------------------------------------------------------------- %
%% Invalid subjects:

% Update which subjects to drop:
job.invalidSubs = []; 
% job.invalidSubs = [1 11 15 19 21 25 26]; % outliers in TAfT and fMRI

% If any subjects given manually, via job.sub2exclude
if ~isempty(job.sub2exclude)
    fprintf('Will exclude subjects %s\n', num2str(job.sub2exclude));
    job.invalidSubs = job.sub2exclude;
end

% Remove repetitions:
job.invalidSubs = unique(job.invalidSubs); % remove repetitions

% Infer valid subjects:
job.validSubs   = setdiff(1:job.nSub, job.invalidSubs);
job.nValidSubs  = length(job.validSubs);

% ----------------------------------------------------------------------- %
%% Complete settings based on response setting:

if isfield(job, 'ROI2use')
        job.nCond       = 2;
        job.condOrder   = 1:2;
        job.colMat      = [0 0 1; 1 0 0];
        job.lineStyle   = {'-', '-'};    
        job.condNames   = {'low BOLD', 'high BOLD'};
        
elseif strcmp(job.responseSettings, 'Go')
    
    if strcmp(job.outcomeSettings, 'rel')  
        job.nCond       = 4;
        job.condOrder   = 1:4;
        job.colMat      = [0 0.6 0.2; 0.8 0 0; 0 0.6 0.2; 0.8 0 0];
        job.lineStyle   = {'-', '-', '--', '--'};
        job.condNames   = {'Go Preferred', 'Go Non-preferred', 'NoGo Preferred', 'NoGo Non-preferred'};
        
    elseif strcmp(job.outcomeSettings, 'abs')
        job.nCond       = 6;
        job.condOrder   = 1:6;
        job.colMat      = [0 0.6 0.2; 0.7 0.7 0.7; 0.8 0 0; 0 0.6 0.2; 0.7 0.7 0.7; 0.8 0 0];
        job.lineStyle   = {'-', '-', '-', '--', '--', '--'};
        job.condNames   = {'Go Reward', 'Go Neutral', 'Go Punishment', 'NoGo Reward', 'NoGo Neutral', 'NoGo Punishment'};
        
    elseif strcmp(job.outcomeSettings, 'all')
        job.nCond       = 8;
        job.condOrder   = 1:8;
        job.colMat      = [0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0; 0.8 0 0; 0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0; 0.8 0 0];
        job.lineStyle   = {'-', '-', '-', '-', '--', '--', '--', '--'};
        job.condNames   = {'Go Reward', 'Go No Reward', 'Go No Punishment', 'Go Punishment', 'NoGo Reward', 'NoGo No Reward', 'NoGo No Punishment', 'NoGo Punishment'};
        
    else
        error('Unknown outcome setting');
    end
    
elseif strcmp(job.responseSettings, 'Hand')
    
    if strcmp(job.outcomeSettings, 'rel')
        job.nCond       = 6;
        job.condOrder   = 1:6;
        job.colMat      = [0 0.6 0.2; 0.8 0 0; 0 0.6 0.2; 0.8 0 0; 0 0.6 0.2; 0.8 0 0];
        job.lineStyle   = {'-', '-', ':', ':', '--', '--'};
        job.condNames   = {...
            'LeftGo Preferred', 'LeftGo Non-preferred', ...
            'RightGo Preferred', 'RightGo Non-preferred', ...
            'NoGo Preferred', 'NoGo Non-preferred'};

    elseif strcmp(job.outcomeSettings, 'abs')
        job.nCond       = 9;
        job.condOrder   = 1:9;
        job.colMat      = [0 0.6 0.2; 0.7 0.7 0.7; 0.8 0 0; 0 0.6 0.2; 0.7 0.7 0.7; 0.8 0 0; 0 0.6 0.2; 0.7 0.7 0.7; 0.8 0 0];
        job.lineStyle   = {'-', '-', '-', ':', ':', ':', '--', '--', '--'};
        job.condNames   = {...
            'LeftGo Reward', 'LeftGo Neutral', 'LeftGo Punishment', ...
            'RightGo Reward', 'RightGo Neutral', 'RightGo Punishment', ...
            'NoGo Reward', 'NoGo Neutral', 'NoGo Punishment'};
        
    elseif strcmp(job.outcomeSettings, 'all')
        job.nCond       = 12;
        job.condOrder   = 1:12;
        job.colMat      = [0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0; 0.8 0 0; 0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0; 0.8 0 0; 0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0; 0.8 0 0];
        job.lineStyle   = {'-', '-', '-', '-', ':', ':', ':', ':', '--', '--', '--', '--'};
        job.condNames   = {...
            'LeftGo Reward', 'LeftGo No Reward', 'LeftGo No Punishment', 'LeftGo Punishment', ...
            'RightGo Reward', 'RightGo No Reward', 'RightGo No Punishment', 'RightGo Punishment', ...
            'NoGo Reward', 'NoGo No Reward', 'NoGo No Punishment', 'NoGo Punishment'};
        
    else
        error('Unknown outcome setting');
    end
    
elseif strcmp(job.responseSettings, 'none')
    
    if strcmp(job.outcomeSettings, 'rel')  
        job.nCond       = 2;
        job.condOrder   = 1:2;
        job.colMat      = [0 0.6 0.2; 0.8 0 0];
        job.lineStyle   = {'-', '-'};
        job.condNames   = {'Preferred', 'Non-preferred'};

        elseif strcmp(job.outcomeSettings, 'abs')
        job.nCond       = 3;
        job.condOrder   = 1:3;
        job.colMat      = [0 0.6 0.2; 0.7 0.7 0.7; 0.8 0 0];
        job.lineStyle   = {'-', '-', '-'};
        job.condNames   = {'Reward', 'Neutral', 'Punishment'};
        
    elseif strcmp(job.outcomeSettings, 'all')
        job.nCond       = 4;
        job.condOrder   = 1:4;
        job.colMat      = [0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0; 0.8 0 0];
        job.lineStyle   = {'-', '-', '-', '-'};
        job.condNames   = {'Reward', 'No Reward', 'No Punishment', 'Punishment'};
        
    else
        error('Unknown outcome setting');
    end
    
else
    error('Unknown response setting');
end

% ----------------------------------------------------------------------- %
%% Handle for saving:

% Distinguish BOLD bins vs. contrast:
if isfield(job, 'ROI2use') % % In case split up per BOLD
    job.plotName = sprintf('BOLD_%s', strjoin(job.ROI2use, ''));

else % otherwise contrast, response, outcome settings
    job.plotName = sprintf('%s_%s_%s', ...
        job.contrastType, job.responseSettings, job.outcomeSettings);
end

% Add channels and baseline:
job.plotName    = sprintf('%s_%s_baseCor_%s_baseTime_%d_%d', ...
    job.plotName, job.chanArea, job.baselineSettings, ...
    job.baselineTimings(1)*1000, job.baselineTimings(end)*1000);

end % END OF FUNCTION.