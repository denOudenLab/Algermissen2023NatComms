function EEGfMRIPav_OutcomeLocked_8_TF_create(rootDir, subVec, baselineSettings, baselineTimings, ...
    responseSettings, outcomeSettings, ...
    TFtype, nFreq, outERPcor, outERPresponseSettings, outERPoutcomeSettings)

% EEGfMRIPav_OutcomeLocked_8_TF_create(rootDir, subVec, baselineSettings, baselineTimings, ...
%     responseSettings, outcomeSettings ,...
%     TFtype, nFreq, outERPcor, outERPresponseSettings, outERPoutcomeSettings)
% 
% Load EEG, rejected trials, behavior; apply baseline correction; sort
% trials into specified conditions, aggregate within conditions, save.
% 
% INPUTS:
% rootDir           = string, root directory of project.
% subVec            = vector of integers, numbers of subjects to include
% (default: 1:36).
% baselineSettings  = string, type of baseline correction, 'trend'
% (regression-based), 'all' (grand mean oer subject), 'condition'
% (grand mean separately per condition), 'trial' (per trial, own
% implementation), 'ft_trial' (with Fieldtrip's ft_freqbaseline).
% baselineTimings   = vector of 1 or 2, timing of baseline correction
% which to use.
% responseSettings   = string, type of response setting for which
%   conditions are split, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo).
% outcomeSettings    = string, type of outcome coding for which 
%   conditions are split, 'abs' (positive/neutral/negative), 'rel' (positive/negative), 'all' (reward/no reward/ no punishment/ punishment).
% TFtype            = type of TF decomposition, 'hanning' or 'morlet'.
% nFreq             = numeric, number of frequency bins decomposed
% (default: 15).
% outERPcor         = Boolean, read data corrected for outcome-locked 
% ERP (true) or not (false), default false.
% outERPresponseSettings   = string, type of response setting for which
%   conditions are split for ERP correction, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo), optional (default: 'Go').
% outERPoutcomeSettings    = string, type of outcome coding for which 
%   conditions are split for ERP correction, 'abs' (positive/neutral/negative), 'rel' (positive/negative), 'all' (reward/no reward/ no punishment/ punishment), optional (default: 'all').
%
% OUTPUTS:
% saves data to disk under dirs.TFgroup.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

% ----------------------------------------------------------------------- %
%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootDir', 'var')
    rootDir = grouplevel_set_rootDir(); % '/project/3017042.02';
    fprintf('rootDir unspecified, assume %s\n', rootDir);
end

if ~exist('subVec', 'var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n', strjoin(string(subVec), ', '));
end

if ~exist('baselineSettings', 'var')
    baselineSettings = 'trend'; % trend, all, condition, trial, ft_trial
    fprintf('baselineSettings unspecified, assume %s\n', baselineSettings);
end

if ~exist('baselineTimings', 'var')
    baselineTimings = [-0.250 -0.050]; % 
    fprintf('baselineTimings unspecified, assume %s\n', strjoin(string(baselineTimings), ' - '));
end

if ~exist('responseSettings', 'var')
    responseSettings = 'Go'; % Go Hand none
    fprintf('responseSettings unspecified, assume %s\n', responseSettings);
end

if ~exist('outcomeSettings', 'var')
    outcomeSettings = 'all'; % abs rel all
    fprintf('outcomeSettings unspecified, assume %s\n', outcomeSettings);
end

if ~exist('TFtype', 'var')
    TFtype = 'hanning'; % hanning, morlet
    fprintf('TFtype unspecified, assume %s\n', TFtype);
end

if ~exist('nFreq', 'var')
    nFreq = 3; 
    fprintf('nFreq unspecified, assume %d\n', nFreq);
end

if ~exist('outERPcor', 'var')
    outERPcor = false;
end

if exist('outERPcor', 'var') && ~exist('outERPresponseSettings', 'var')
    outERPresponseSettings = 'Go';
end

if exist('outERPcor', 'var') && ~exist('outERPoutcomeSettings', 'var')
    outERPresponseSettings = 'all';
end
    
% ----------------------------------------------------------------------- %
%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootDir, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));

% Set directories and parameters:
dirs        = set_dirs(rootDir);
par         = set_par();

% Add TF directory based on TF decomposition settings:
dirs.TF     = fullfile(dirs.results, sprintf('TF_subjectlevel_%s%d', TFtype, nFreq));

if outERPcor
    outERPSettings  = sprintf('_outERPcor_%s_%s', outERPresponseSettings, outERPoutcomeSettings);
    dirs.TF         = [dirs.TF outERPSettings];
end

% ----------------------------------------------------------------------- %
%% Automatically fill in downstream settings:

fprintf('Retrieve condition information\n');

nCond = 1; % initialize

if strcmp(responseSettings, 'Go')
    allResp = [1 0]; % Go/NoGo.
    nCond   = nCond * 2;
elseif strcmp(responseSettings, 'Hand')
    allResp = [101 97 0]; % left/right/nogo.
    nCond   = nCond * 3;
elseif strcmp(responseSettings, 'none')
    allResp = 1; % just anything
    nCond   = nCond * 3;
else
    error('Unknown response setting');
end

if strcmp(outcomeSettings, 'abs')
    allFb = [1 0 -1]; % reward/neutral/punishment
    nCond = nCond * 3;
elseif strcmp(outcomeSettings, 'rel')
    allFb = [1 2]; % preferred/non-preferred
    nCond = nCond * 2;
elseif strcmp(outcomeSettings, 'all')
    allFb = 1:4; % reward/no reward/no punishment/punishment
    nCond = nCond * 4;
else
    error('Unknown response setting');
end

% ----------------------------------------------------------------------- %
%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.TF, sprintf('*_outcomelocked.mat')));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n', nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList, 8, 10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames, subVec));
nSub        = length(selIndices);

% ----------------------------------------------------------------------- %
%% Create final file name:

outputFile = sprintf('EEGfMRIPav_TF_baseCor_%s_baseTime_%03dms_%03dms_resp_%s_fb_%s_hanning%d', ...
    baselineSettings, baselineTimings(1)*-1000, baselineTimings(end)*-1000, responseSettings, outcomeSettings, nFreq);

if outERPcor
    outputFile = [outputFile outERPSettings];
end

% Add directory and .mat ending:
outputFile = fullfile(dirs.TFgroup, sprintf('%s.mat', outputFile));

fprintf('Selected file is %s\n', outputFile);

% ----------------------------------------------------------------------- %
%% Loop over subjects:

if exist(outputFile, 'file')
    
    warning('File %s already exists, not created again', outputFile);
    
else

    TFall = cell(nSub, 1);
    nTrialCond = nan(nSub, nCond);
    
    for iSub = selIndices % %iSub = 1;

        % --------------------------------------------------------------- %
        %% 1) Load EEG data:
        
        fprintf('Subject %03d: Start loading TF data ...\n', iSub);
        load(fullfile(dirs.TF, fileList{iSub}));
        fprintf('Subject %03d: Finished loading TF data\n', iSub);
        
        % --------------------------------------------------------------- %
        %% 2) Load rejected trials:

        fprintf('Subject %03d: Load rejected trials\n', iSub);
        filenameRejTrials = fullfile(dirs.rejTrials, sprintf('MRIpav_%03d_rejTrials.mat', iSub));
        load(filenameRejTrials);

        % --------------------------------------------------------------- %
        %% 3) Load and recode behavior:
        
        fprintf('Subject %03d: Load behavioral data\n', iSub);
        behav = load_recode_behav(rootDir, iSub, rejectedtrials); 
               
        % --------------------------------------------------------------- %
        %% 4) Baseline-correction per subject (and transformation into decibel): 
        
        % Retrieve baseline timings:
        timeIdx = dsearchn(freq.time', baselineTimings'); % select indices for baseline time window
        if length(timeIdx) > 1; timeIdx = timeIdx(1):timeIdx(end); end % extend to interval if start and end given
        fprintf('Subject %03d: Compute baseline, setting %s\n', iSub, baselineSettings);

        % --------------------------------------------------------------- %
        % Perform baseline correction based on given setting:
        if strcmp(baselineSettings, 'all')
            baseline        = nanmean(nanmean(freq.powspctrm(:, :, :, timeIdx), 1), 4); % first average over trials, then over time bins
            freq.powspctrm  = 10*log10(freq.powspctrm ./ ...
                repmat(baseline, size(freq.powspctrm, 1), 1, 1, size(freq.powspctrm, 4)));

        % --------------------------------------------------------------- %
        elseif strcmp(baselineSettings, 'condition')
            
           for iResp = allResp
                for iVal = allFb % win/avoid cue.
                   iCondi = iCondi + 1;
                   fprintf('Subject %03d condition %d: Compute baseline subject per val and %s for iResp %d, iVal %d\n', iSub, iCondi, responseSettings, iResp, iVal)
                   trials = find(behav.val == iVal & behav.resp == iResp); % based on exact response
                    % Compute baseline averaged over trials:
                    baseline     = nanmean(nanmean(baselineData.powspctrm(trials, :, :, timeIdx), 1), 4); % first average over time bins, then over trials     
                    % Correct each trial by that baseline:
                    freq.powspctrm(trials, :, :, :)    = 10*log10(freq.powspctrm(trials, :, :, :) ./ ...
                        repmat(baseline, length(trials), 1, 1, size(freq.powspctrm, 4)));
                end % end of iVal
            end % end of iResp

        % --------------------------------------------------------------- %
        elseif strcmp(baselineSettings, 'trial')
            fprintf('Subject %03d: Compute baseline per trial\n', iSub);
            % Average over specified time, then repeat across time dimension
            baseline   = nanmean(freq.powspctrm(:, :, :, timeIdx), 4); % keep trial, channel, freq bin, average over time bin
            freq.powspctrm      = 10*log10(freq.powspctrm ./ ...
                repmat(baseline, 1, 1, 1, size(freq.powspctrm, 4)));

        % --------------------------------------------------------------- %
        elseif strcmp(baselineSettings, 'ft_trial')
            fprintf('Subject %03d: Compute baseline per trial using Fieldtrip\n', iSub);
            % no adaptation possible for response-locked data
            cfg                 = [];
            cfg.baseline        = par.TF.baselinetime;
            cfg.baselinetype    = 'db'; % use db; try out relative
            freq                = ft_freqbaseline(cfg, freq);   
     
        % --------------------------------------------------------------- %
        elseif strcmp(baselineSettings, 'trend')
            fprintf('Subject %03d: Compute baseline per trial by detrending\n', iSub);
            freq.powspctrm      = regress_baseline(freq.powspctrm, freq.trialinfo, timeIdx(1), timeIdx(end), 'division');
            freq.powspctrm      = real(10*log10(freq.powspctrm)); % decibel conversion
	else
 	    error('Unknown baseline setting');
        end % end baseline correction
        
        % --------------------------------------------------------------- %
        %% 5) Select trials into conditions:

        cfg = []; iCondi = 0; % iCondi is stimulus number (Response x Valence, i.e. 6)
        for iResp = allResp % iResp = 1; loop over Go/NoGo or GoLeft/GoRight/NoGo
            for iVal = allFb % iVal = 1; loop over preferred/non-preferred or reward/neutral/punishment
               iCondi = iCondi + 1;

               % -------------------------------------------------------- %
               % 1) Select trial numbers of condition:
               fprintf('Subject %03d condition %d: Start TF extraction per %s and %s for iResp = %d, iVal = %d\n', ...
                   iSub, iCondi, responseSettings, outcomeSettings, iResp, iVal);

               % -------------------------------------------------------- %
               % 1a) Select based on response:
               if strcmp(responseSettings, 'Go')
                   respTrials = find(behav.isgo == iResp); % based on Go/NoGo
               elseif strcmp(responseSettings, 'Hand')
                   respTrials = find(behav.resp == iResp); % based on exact response
               elseif strcmp(responseSettings, 'none')
                   respTrials = behav.trlIdx; % include all trials
               else
                   error('Invalid response specification: %s', responseSettings);
               end % 

               % -------------------------------------------------------- %
               % 1b) Select based on feedback:
               if strcmp(outcomeSettings, 'abs')
                   fbTrials = find(behav.fb.abs == iVal); % reward/neutral/punishment
               elseif strcmp(outcomeSettings, 'rel')
                   fbTrials = find(behav.fb.rel == iVal); % preferred/non-preferred
               elseif strcmp(outcomeSettings, 'all')
                   fbTrials = find(behav.fb.all == iVal); % reward/no-reward/no-punishment/punishment
               else
                   error('Invalid outcome specification: %s', outcomeSettings);
               end % end of trial selection

               % -------------------------------------------------------- %
               % 1c) Form intersection:
               cfg.trials               = intersect(respTrials, fbTrials);
               fprintf('Found %d trials \n',length(cfg.trials));
               nTrialCond(iSub, iCondi) = length(cfg.trials);

               % -------------------------------------------------------- %
               % 2) Average over selected trials:
               if ~isempty(cfg.trials) % if any trials in this condition
                   TFall{iSub}.ValxAct{iCondi} = ft_freqdescriptives(cfg, freq); % define for this condition

               % -------------------------------------------------------- %
               else % otherwise fill in NaNs
                   fprintf('Subject %03d condition %d: Found no trials in condition iResp = %d, iVal=%d, store NaNs !!!!!!!!!!!!!!!!!!!!!!\n', ...
                       iSub, iCondi, iResp, iVal);
                   TFall{iSub}.ValxAct{iCondi}.powspctrm = nan(size(freq.powspctrm, 2), ...
                       size(freq.powspctrm,3), size(freq.powspctrm, 4));

               end % end of check for empty trials
            end % end of for loop for iVal
        end % end of for loop for iResp

	% --------------------------------------------------------------- %
        % Store subject info and condition labels.

        fprintf('Store subject %03d\n', iSub);
 
        clear freq seq par baselineData baseline

        % Update effective subject number
        iSub = iSub + 1;

    end % end of for loop for subject  
 
    % ------------------------------------------------------------------- %
    %% Save data for all subjects:
    
    fprintf('Save file under %s\n', outputFile);
    save(outputFile, 'TFall', 'nTrialCond', '-v7.3');
    fprintf('Finished :-)\n');

end % end if ~exist(outputFile)

fprintf('Done :-)\n');

end % END OF FUNCTION.