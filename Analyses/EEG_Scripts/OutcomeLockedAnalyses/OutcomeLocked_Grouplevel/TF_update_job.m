function [job] = TF_update_job(job)

% [job] = TF_update_job(job)
%
% Complete settings for job on TF-domain data to create.
% (data sets to create; plotting settings).
% 
% INPUTS:
% job                   = cell, created via TF_update_job.m, needs at least
% fields:
%   .nSub               = integer, number of subjects.
%   .sub2exclude        = numeric vector, numbers of subjects to exclude.
%   .TFtype             = type of TF decomposition, 'hanning' or 'morlet'.
%   .nFreq              = numeric, number of frequency bins decomposed.
%   (default: 15).
%   .baselineSettings   = string, type of baseline correction, 'trend'
%   (regression-based), 'all' (grand mean oer subject), 'condition'
%   (grand mean separately per condition), 'trial' (per trial, own
%   implementation), 'ft_trial' (with Fieldtrip's ft_baselinecorrect).
%   .baselineTimings    = vector of 1 or 2, timing of baseline correction
%   for which to load data.
%   .responseSettings   = string, type of response setting for which
%   conditions are split, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo).
%   .outcomeSettings    = string, type of outcome coding for which 
%   conditions are split, 'abs' (positive/neutral/negative), 'rel' (positive/negative), 'all' (reward/no reward/ no punishment/ punishment).
%   .chanArea           = string, area of channels to select, 'frontal',
%   'completefrontal', 'midfrontal'.
%   .band               = string, frequency band to select, 'delta', 'theta', 'thetadelta', 'lowalpha', 'verylowalpha', 'alpha', 'beta', 'middlebeta', 'broad'.
%   .contrastType       = string, contrast to be used: 'Preferred', 'Action',
%   'GoPreferred', 'SalientPreferred', 'SalientAction'.
%   .outERPcor         = Boolean, read data corrected for outcome-locked 
% ERP (true) or not (false), default false.
%   .outERPresponseSettings   = string, type of response setting for which
%   conditions are split for ERP correction, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo), optional (default: 'Go').
%   .outERPoutcomeSettings    = string, type of outcome coding for which 
%   conditions are split for ERP correction, 'abs' (positive/neutral/negative), 'rel' (positive/negative), 'all' (reward/no reward/ no punishment/ punishment), optional (default: 'all').
% 
% OUTPUTS:
% job                   = cell, with the following fields added:
%   .inputFileName      = name of file to load: 
%   .layout             = string, cap to use for topoplots.
%   .TFtiming           = numeric vector of 2 elements, timing for TF
%   plots.
%   .sigTime            = numeric vector of 2 elements, timing for
%   topoplots.
%   .freq               = numeric vector of 2 elements, range of selected
%   frequencies.
%   .bandName           = string, name of selected frequency band
%   capitalized.
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
%   .plotName           = string, important settings to be included in plot
%   names.
%   .twoLineColor       = 2 x 3 matrix, line colors to use for twoLinePlot.
%   .twoLineLabels      = vector of 2 strings, labels for conditions.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

fprintf('Update and complete job settings\n')

% ----------------------------------------------------------------------- %
%% File to load:

% Standard input file:
job.inputFileName = sprintf('EEGfMRIPav_TF_baseCor_%s_baseTime_%03dms_%03dms_resp_%s_fb_%s_hanning%d', ...
    job.baselineSettings, job.baselineTimings(1)*-1000, job.baselineTimings(end)*-1000, ...
    job.responseSettings, job.outcomeSettings, job.nFreqs);

% If ERP-corrected:
if job.outERPcor
    outcomeSettingsERPcor = sprintf('_outERPcor_%s_%s', job.outERPresponseSettings, job.outERPoutcomeSettings);
    job.inputFileName = sprintf('%s%s', job.inputFileName, outcomeSettingsERPcor);
end

% Add final matlab ending:
job.inputFileName = sprintf('%s.mat', job.inputFileName); % add final ending

% ----------------------------------------------------------------------- %
%% Plotting options:  

% Cap layout:
job.layout   = 'easycapM11.mat'; %

% Timing of TF plot:
job.TFtiming = [-0.25 1.0];

% ----------------------------------------------------------------------- %
%% Two line plots: color, labels, time of significant differences:

% Timing for topoplot:
if isfield(job, 'sigTime') % save in case it is given as input; return back in again later
    sigTime = job.sigTime;
end
job.sigTime = []; % initialize/delete

if strcmp(job.contrastType, 'Preferred')

    job.twoLineColor    = [0 0.6 0.2; 0.8 0 0]; % green, red
    job.twoLineLabels   = {'Positive', 'Negative'};

    if strcmp(job.band, 'theta')
        job.sigTime     = [0.2500 0.4750]; % 
    elseif strcmp(job.band, 'beta') || strcmp(job.band, 'middlebeta')
        job.sigTime     = [0.3000 1.1750]; % 
    else
        job.sigTime     = [];
    end
        
elseif strcmp(job.contrastType, 'Action')

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels   = {'Go', 'NoGo'};

    if strcmp(job.band, 'theta')
        job.sigTime     = [0.000 0.2000]; % 
    elseif strcmp(job.band, 'beta')
        job.sigTime     = [0.200 0.4500]; % 
    else
        job.sigTime     = [];
    end
        
elseif strcmp(job.contrastType, 'GoPreferred')

    job.twoLineColor    = [0 0.6 0.2; 0.8 0 0]; % green, red
    job.twoLineLabels   = {'preferred Go', 'non-preferred Go'};

    if strcmp(job.band, 'theta')
        job.sigTime     = [0.2250 0.5000]; % 
    elseif strcmp(job.band, 'beta')
        job.sigTime     = [0.4500 0.8500]; % 
    else
        job.sigTime     = [];
    end
                        
elseif strcmp(job.contrastType, 'SalientPreferred')

    job.twoLineColor    = [0 0.6 0.2; 0.8 0 0]; % green, red
    job.twoLineLabels   = {'Reward', 'Punishment'};

    if strcmp(job.band, 'theta')
        job.sigTime     = [0.2250 0.5000]; % 
    else
        job.sigTime     = [0.1250 0.3750 0.3750 0.5750 0.5750 0.8250];
    end
        
elseif strcmp(job.contrastType, 'SalientAction')

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels   = {'Salient Go', 'Salient NoGo'};

    if strcmp(job.band, 'theta')
        job.sigTime     = [0 0.375]; % 
    elseif strcmp(job.band, 'lowalpha')
        job.sigTime     = [0 0.375]; % 
    elseif strcmp(job.band, 'beta')
        job.sigTime     = [0.200 0.425];
    end
                
else
    warning('unknown sigTime')
    job.sigTime = [];
end

% ----------------------------------------------------------------------- %
%% Frequency and band name of interest:

if strcmp(job.band, 'delta')
    job.freq        = [1 4]; % delta range
    job.bandName    = 'Delta';
    
elseif strcmp(job.band, 'theta')
    job.freq        = [4 8]; % theta range
    job.bandName    = 'Theta';
    
elseif strcmp(job.band, 'thetadelta')
    job.freq        = [1 8]; % theta range
    job.bandName    = 'Theta';

elseif strcmp(job.band, 'alpha')
    job.freq        = [8 13]; % alpha range
    job.bandName    = 'Alpha';
    
elseif strcmp(job.band, 'lowalpha')
    job.freq        = [6 10]; % alpha range
    job.bandName    = 'Low Alpha';
    
elseif strcmp(job.band, 'verylowalpha')
    job.freq        = [7 9]; % alpha range
    job.bandName    = 'Low Alpha';
    
elseif strcmp(job.band, 'beta')
    job.freq        = [13 30]; % beta range
    job.bandName    = 'Beta';
    
elseif strcmp(job.band, 'middlebeta')
    job.freq        = [14 28]; % beta range
    job.bandName    = 'Beta';
    
elseif strcmp(job.band, 'broad')
    job.freq        = [1 33]; % broadband range
    job.bandName    = 'Broadband';
else
    error('No frequency band set\n');
end

% ----------------------------------------------------------------------- %
%% Channels and labels:

if strcmp(job.chanArea, 'frontal')
    job.channels = {'AF3', 'AF4', 'Fz', 'F1', 'F2', 'FCz', 'FC1', 'FC2'};
    job.chanName = 'Frontal';

elseif strcmp(job.chanArea, 'completefrontal') % for complete frontal ROIs
    job.channels = {'Fpz', 'Fp1', 'Fp2', 'AFz', 'AF3', 'AF4', 'AF7', 'AF8', ...
        'Fz', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', ...
        'FCz', 'FC1', 'FC2', 'FC3', 'FC4', 'FC5', 'FC6', 'FC7', 'FC8'};
    job.chanName = 'Frontal';

elseif strcmp(job.chanArea, 'midfrontal')
    job.channels = {'Fz', 'FCz', 'Cz'};
    job.chanName = 'Midfrontal';

else
    warning('No channels set\n');
end

% Sort channels alphebetically:
job.channels = sort(job.channels);

% ----------------------------------------------------------------------- %
%% Invalid subjects:

% Update which subjects to drop:
job.invalidSubs = []; 
% job.invalidSubs = [11 12 18 24 30 34]; % some cell sizes < 20
% job.invalidSubs = [11 12 15 23 25 26 30]; % outliers in TAfT and fMRI

% If any subjects given manually, via job.sub2exclude
if ~isempty(job.sub2exclude)
    fprintf('Will exclude manually specified subjects %s\n', num2str(job.sub2exclude));
    job.invalidSubs = job.sub2exclude;
end

% Remove repetitions:
job.invalidSubs = unique(job.invalidSubs); % remove repetitions

% Infer valid subjects:
job.validSubs   = setdiff(1:job.nSub, job.invalidSubs);
job.nValidSubs  = length(job.validSubs);

% ----------------------------------------------------------------------- %
%% Condition names and labels:

% https://nl.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/24497/versions/2/screenshot.PNG
% for colors, see https://github.com/jkitchin/matlab-cmu/blob/master/%2Bcmu/colors.m

if strcmp(job.responseSettings, 'Go')
    
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
        job.condNames   = {'LeftGo Preferred', 'LeftGo Non-preferred', 'RightGo Preferred', 'RightGo Non-preferred', 'NoGo Preferred', 'NoGo Non-preferred'};
        
    elseif strcmp(job.outcomeSettings, 'abs')
        
        job.nCond       = 9;
        job.condOrder   = 1:9;
        job.colMat      = [0 0.6 0.2; 0.7 0.7 0.7; 0.8 0 0; 0 0.6 0.2; 0.7 0.7 0.7; 0.8 0 0; 0 0.6 0.2; 0.7 0.7 0.7; 0.8 0 0];
        job.lineStyle   = {'-', '-', '-', ':', ':', ':', '--', '--', '--'};
        job.condNames   = {'LeftGo Reward', 'LeftGo Neutral', 'LeftGo Punishment', 'RightGo Reward', 'RightGo Neutral', 'RightGo Punishment', 'NoGo Reward', 'NoGo Neutral', 'NoGo Punishment'};
        
    elseif strcmp(job.outcomeSettings, 'all')
        
        job.nCond       = 12;
        job.condOrder   = 1:12;
        job.colMat      = [0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0; 0.8 0 0; 0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0; 0.8 0 0; 0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0; 0.8 0 0];
        job.lineStyle   = {'-', '-', '-', '-', ':', ':', ':', ':', '--', '--', '--', '--'};
        job.condNames   = {'LeftGo Reward', 'LeftGo No Reward', 'LeftGo No Punishment', 'LeftGo Punishment', 'RightGo Reward', 'RightGo No Reward', 'RightGo No Punishment', 'RightGo Punishment', 'NoGo Reward', 'NoGo No Reward', 'NoGo No Punishment', 'NoGo Punishment'};
        
    else
        error('Unknown outcome setting');
    end
    
elseif strcmp(job.responseSettings, 'none')
    
    if strcmp(job.outcomeSettings, 'all')
        
        job.nCond       = 4;
        job.condOrder   = 1:4;
        job.colMat      = [0 0.6 0.2; 0.65 0.73 0.67; 0.76 0.61 0.62; 0.8 0 0]; % 165 185 170 and 195 155 157
        job.lineStyle   = {'-', '-', '-', '-'};
        job.condNames   = {'Reward', 'Neutral (no reward)', 'Neutral (no punishment)', 'Punishment'};
        
    else
        error('Unknown outcome setting');
    end
    
else
    fprintf('Unknown response setting');
end

% ----------------------------------------------------------------------- %
%% Handle for saving:

% Baseline name: contrast, response, outcome, channels, baseline type and
% timing, frequency band
job.plotName = sprintf('%s_%s_%s_%s_baseCor_%s_baseTime_%d_%d_%s', ...
    job.contrastType, job.responseSettings, job.outcomeSettings, job.chanArea, ...
    job.baselineSettings, job.baselineTimings(1)*1000, job.baselineTimings(end)*1000, ...
    job.bandName);

% If ERP-corrected:
if job.outERPcor
    job.plotName = sprintf('%s%s', job.plotName, outcomeSettingsERPcor);
end

end % END OF FUNCTION.