function EEGfMRIPav_OutcomeLocked_8_time_create(rootDir, subVec, baselineSettings, baselineTimings, ...
    responseSettings, outcomeSettings, ROI2use)

% EEGfMRIPav_OutcomeLocked_8_time_create(rootDir, subVec, baselineSettings, baselineTimings, ...
%     responseSettings, outcomeSettings, ROI2use)
% 
% Load EEG, rejected trials, behavior; apply baseline correction; sort
% trials into specified conditions, aggregate within conditions, save.
% 
% INPUTS:
% rootDir           = string, root directory of project.
% subVec            = vector of integers, numbers of subjects to include
% (default: 1:36).
% baselineSettings  = string, type of baseline correction, 'trend'
% (regression-based), 'ft' (with Fieldtrip's ft_timelockbaseline).
% baselineTimings   = vector of 1 or 2, timing of baseline correction
% which to use.
% responseSettings  = string, type of response setting for which
%   conditions are split, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo).
% outcomeSettings   = string, type of outcome coding for which 
%   conditions are split, 'abs' (positive/neutral/negative), 'rel' (positive/negative), 'all' (reward/no reward/ no punishment/ punishment).
% ROI2use           = string, split data by ROI (instead of by conditions), optional.
%
% OUTPUTS:
% saves data to disk under dirs.timegroup.
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

% Potential ROIs to use from GLM1:
% ROI2use = {'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1LeftMotorConj', 'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1LeftITGConj', 'GLM1V1Conj'};

% ----------------------------------------------------------------------- %
%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootDir, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));

% Set directories and parameters:
dirs        = set_dirs(rootDir);
par         = set_par();

% ----------------------------------------------------------------------- %
%% Automatically fill in downstream settings:

nCond = 1; % initialize
if strcmp(responseSettings, 'Go')
    allResp = [1 0]; % Go/NoGo.
    nCond   = nCond * 2;
elseif strcmp(responseSettings, 'Hand')
    allResp = [101 97 0]; % left/right/nogo.
    nCond   = nCond * 3;
elseif strcmp(responseSettings, 'none')
    allResp = 1; % just anything
    nCond   = nCond * 1;
else
    fprintf('Unknown response setting')
end

if strcmp(outcomeSettings, 'abs')
    allFb   = [1 0 -1]; % reward/neutral/punishment
    nCond   = nCond * 3;
elseif strcmp(outcomeSettings, 'rel')
    allFb   = [1 2]; % preferred/non-preferred
    nCond   = nCond * 2;
elseif strcmp(outcomeSettings, 'all')
    allFb   = 1:4; % reward/no reward/no punishment/punishment
    nCond   = nCond * 4;
else
    fprintf('Unknown response setting')
end

% Overwrite if ROI specified:
if exist('ROI2use', 'var')
    nCond   = 2;
end

% ----------------------------------------------------------------------- %
%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.finalPP, '*_finalPP.mat'));
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

% Create filename:
outputFile = sprintf('EEGfMRIPav_time_baseCor_%s_baseTime_%03dms_%03dms', ...
    baselineSettings, baselineTimings(1)*-1000, baselineTimings(end)*-1000);

if exist('ROIs2use', 'var')
    outputFile = [outputFile sprintf('_BOLD_%s', strjoin(ROIs2use, ''))]; % add ROI name (response and outcome irrelevant)
else
    outputFile = sprintf('%s_resp_%s_fb_%s', outputFile, responseSettings, outcomeSettings); % otherwise add resposne and outcome settings
end

% Add directory and .mat ending:
outputFile = fullfile(dirs.timegroup, sprintf('%s.mat', outputFile));

fprintf('Selected file is %s\n', outputFile)

% ----------------------------------------------------------------------- %
%% Loop over subjects:

if exist(outputFile, 'file')
    
    warning('File %s already exists, not created again', outputFile);
    
else

    ERPdata     = cell(nSub, 1);
    nTrialCond  = nan(nSub, nCond);
    
    for iSub = selIndices % %iSub = 1;

        % --------------------------------------------------------------- %
        %% 1) Load EEG data:

        fprintf('Subject %03d: Load time-domain EEG data\n', iSub);
        load(fullfile(dirs.finalPP, fileList{iSub}))
        scd = rmfield(scd, 'elec'); % remove elec field
        fprintf('Subject %03d: Finished loading\n', iSub);
               
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
        %% 4) Baseline correction:
        
        timeIdx = dsearchn(scd.time{1}', baselineTimings'); % select indices for baseline time window
        if length(timeIdx) > 1; timeIdx = timeIdx(1):timeIdx(end); end % extend to interval if start and end given
        
        if strcmp(baselineSettings, 'ft')
            fprintf('Subject %03d: Subtract baseline per trial using Fieldtrip\n', iSub)
            cfg         = [];
            cfg.baseline= baselineTimings; % -250 - 50 ms
            scd         = ft_timelockbaseline(cfg, scd);
        elseif strcmp(baselineSettings, 'trend')
            fprintf('Subject %03d: Subtract baseline per trial by detrending\n', iSub)
            scd.trial   = regress_baseline(scd.trial, scd.trialinfo, timeIdx(1), timeIdx(end), 'subtraction');
            scd.time    = scd.time{1}; % all trials with identical timing
            scd.dimord  = 'rpt_chan_time';
        end
        
        % --------------------------------------------------------------- %
        %% 5) Load trial-by-trial BOLD amplitude betas:
           
        if exist('ROIs2use', 'var')

            fprintf('Add path to TAfT\n');
            addpath(fullfile(rootDir, '/Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT'));
            fprintf('Add path to SPM\n');
            addpath /home/common/matlab/spm12;

            % Initialize job:
            job.subID       = iSub;
            job.RT          = 1.4; % original repetition time
            job.hpass       = 128; % high-pass filter about 128 sec.
            job.ups         = 10; % gives resolution of 140 ms, should be sufficient (Hauser Hunt et al. used 185 ms)
            job.ons_unit    = 1; % Boolean that onsets are given in seconds
            job.trialdur    = 8; % Trial durations are 4700 - 6650 ms (8 sec. as Hauser, Hunt et al. 2015)
            job.demeanX     = true;
            job.demeanY     = true;
            job.add_intercept = false; % add intercept for GLMs
            job.HRFtype     = 'trial';
            job.save        = false;

            % Files per block:
            dirs.fMRI = fullfile(rootDir, '/Log/fMRI')';    
            for iBlock = 1:6 % loop through blocks iBlock = 2;
                job.ROIs(1).ROIname = char(ROIs2use); %  name of ROI just carried forward
                job.ROIs(1).ROIdef(iBlock).rawfMRIfile  = {fullfile(dirs.fMRI, ...
                    sprintf('sub-%03d/FEAT_Block%d.feat/AROMA/%s.txt', iSub, iBlock, job.ROIs(1).ROIname))}; % where to load volumes from
                job.ROIs(1).ROIdef(iBlock).nuisance     = {fullfile(dirs.fMRI, ...
                    sprintf('sub-%03d/FEAT_Block%d.feat/AROMA/NuisanceRegressors.txt', iSub, iBlock))}; % where to load volumes from
                job.ROIs(1).ROIdef(iBlock).onsets       = {fullfile(dirs.fMRI, ...
                    sprintf('sub-%03d/FEAT_Block%d.feat/AROMA/outcomeOnsets.txt', iSub, iBlock))}; % trial onsets
                job.ROIs(1).ROIdef(iBlock).betafMRIfile = {fullfile(dirs.fMRI, ...
                    sprintf('fMRI_TAfT_Betas/betas_ROI_%s_sub%03d_block%d.txt', job.ROIs(1).ROIname, iSub, iBlock))};
                job.ROIs(1).ROIdef(iBlock).fitType      = 'HRF';
            end

            % Upsample ROI-data:
            fprintf('Subject %03d Block %d: Upsample data\n', iSub, iBlock);
            out         = taft_loop_upsample_fit(job); % dat is of format trials x time bins
            BOLDdata    = [out.ROIs{:}.ROIdef{:}]; % concatenate blocks

            % Select valid trials:
            X           = BOLDdata(scd.trialinfo);      

            % Determine quantiles:
            lowQ        = quantile(X, 1/3); % lowest third cutoff
            highQ       = quantile(X, 2/3); % highest third cutoff
            Xlow        = X <= lowQ; % lowest third indices
            Xhigh       = X >= highQ; % highest third indices
            Xsel        = Xlow + 2*Xhigh; % low becomes 1, high becomes 2
        end

        % --------------------------------------------------------------- %
	%% 6) Select trials per condition:

        if exist('ROIs2use', 'var')
            for iCondi = 1:nCond

                cfg                                 = [];
                cfg.trials                          = find(Xsel == iCondi);
                nTrialCond(iSub, iCondi)            = length(cfg.trials);
                fprintf('Subject %03d condition %d of %s: Found %d trials \n', ...
                    iSub, iCondi, char(ROIs2use), length(cfg.trials));

                   if ~isempty(cfg.trials) % if any trials in this condition
                       ERPdata{iSub, iCondi}        = ft_timelockanalysis(cfg, scd);

                   else % otherwise fill in NaNs
                       fprintf('Subject %03d condition %d: Found no trials, store NaNs !!!!!!!!!!!!!!!!!!!!!!\n', ...
                           iSub, iCondi);
                       ERPdata{iSub, iCondi}.avg    = nan(length(scd.label), length(scd.time{1}));

                   end % end of check for empty trials
            end % end of iCondi

        else % use response and outcome settings

            for iResp = allResp % iResp = 1; loop over Go/NoGo or GoLeft/GoRight/NoGo
                for iVal = allFb % iVal = 1; loop over preferred/non-preferred or reward/neutral/punishment
                   iCondi = iCondi + 1;

                   % 1) Select trial numbers of condition:
                   fprintf('Subject %03d condition %d: Start trial extraction per %s and %s for iResp = %d, iVal = %d\n', ...
                       iSub, iCondi, responseSettings, outcomeSettings, iResp, iVal);

                   % Select based on response:
                   if strcmp(responseSettings, 'Go')
                       respTrials = find(behav.isgo == iResp); % based on Go/NoGo
                   elseif strcmp(responseSettings, 'Hand')
                       respTrials = find(behav.resp == iResp); % based on exact response
                   elseif strcmp(responseSettings, 'none')
                       respTrials = behav.trlIdx; % include all trials
                   else
                       error('Invalid response specification: %s', responseSettings);
                   end % 

                   % Select based on feedback:
                   if strcmp(outcomeSettings, 'abs')
                       fbTrials = find(behav.fb.abs == iVal); % reward/neutral/punishment
                   elseif strcmp(outcomeSettings, 'rel')
                       fbTrials = find(behav.fb.rel == iVal); % preferred/non-preferred
                   elseif strcmp(outcomeSettings, 'all')
                       fbTrials = find(behav.fb.all == iVal); % reward/no-reward/no-punishment/punishment
                   else
                       error('Invalid outcome specification: %s', outcomeSettings);
                   end % end of trial selection

                   % Form intersection:
                   cfg                              = [];
                   cfg.trials                       = intersect(respTrials, fbTrials);
                   fprintf('Found %d trials \n', length(cfg.trials));
                   nTrialCond(iSub, iCondi)         = length(cfg.trials);
                   % Test with selectdata first:

                   % 2) Average over selected trials:
                   if ~isempty(cfg.trials) % if any trials in this condition
                       ERPdata{iSub, iCondi}        = ft_timelockanalysis(cfg, scd);
                   else % otherwise fill in NaNs
                       fprintf('Subject %03d condition %d: Found no trials in condition iResp = %d, iVal=%d, store NaNs !!!!!!!!!!!!!!!!!!!!!!\n', iSub, iCondi, iResp, iVal);
                       ERPdata{iSub, iCondi}.avg    = nan(length(scd.label),length(scd.time{1}));
                   end % end of check for empty trials
                end % end of for loop for iVal
            end % end of for loop for iResp
        end % end of case distinction ROI or not
        fprintf('Subject %03d: Finished\n', iSub);
        clear scd behav behav
    end % end loop iSub
    
    % ------------------------------------------------------------------- %
    %% Save data for all subjects:
    
    fprintf('Save file under %s\n', outputFile);
    save(outputFile, 'ERPdata', 'nTrialCond', '-v7.3');
    fprintf('Finished :-)\n')

end % end if ~exist(outputFile)

fprintf('Done :-)\n');

end % END OF FUNCTION.