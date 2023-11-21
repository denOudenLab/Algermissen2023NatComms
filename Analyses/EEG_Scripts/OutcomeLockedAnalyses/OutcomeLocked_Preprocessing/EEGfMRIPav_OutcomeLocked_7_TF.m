function EEGfMRIPav_OutcomeLocked_7_TF(rootDir, subVec, stimERPcor, responseSettings, outcomeSettings)

% EEGfMRIPav_OutcomeLocked_7_TF(rootDir, subVec, stimERPcor, responseSettings, outcomeSettings)
% 
% Reject trials, compute surface LaPlacian filter (and at the same time
% interpolate bad channels). Final steps of pre-processing.
% 
% INPUTS:
% rootDir           = string, root directory of project.
% subVec            = vector of integers, numbers of subjects to process.
% stimERPcor        = Boolean, whether to trial-wise subtract condition-wise 
% ERP before TF decompositions or not (default: false).
% responseSettings  = string, for condition-wise ERP correction: split conditions by Go/NoGo ('Go') or left Go/right Go/NoGo ('Hand') (default: Go).
% outcomeSettings   = string, for condition-wise ERP correction: split conditions by absolute feedback reward/neutral/punishment ('abs') or relative feedback preferred/non-preferred ('rel') or exact feedback reward/no reward/no punishment/punishment ('all') (default: all).
%
% OUTPUTS:
% saves data to disk under dirs.TF.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Preprocessing/

% ----------------------------------------------------------------------- %
%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootDir', 'var')
    rootDir = preprocessing_set_rootDir(); % '/project/3017042.02';
    fprintf('rootDir unspecified, assume %s\n',rootDir)
end

if ~exist('subVec', 'var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n', strjoin(string(subVec), ', '));
end

if ~exist('stimERPcor', 'var')
    stimERPcor = false;
    fprintf('stimERPcor unspecified, perform no condition-wise ERP correction');
end

if stimERPcor && ~exist('responseSettings', 'var')
    responseSettings = 'Go';
    fprintf('responseSettings for condition-wise ERP correction unspecified, assume %s\n',responseSettings);
end

if stimERPcor && ~exist('outcomeSettings', 'var')
    outcomeSettings = 'all';
    fprintf('outcomeSettings for condition-wise ERP correction unspecified, assume %s\n',outcomeSettings);
end
    
% ----------------------------------------------------------------------- %
%% Settings for ERP correction:

if stimERPcor
    fprintf('Correct for stimulus-locked ERP with response setting %s, outcome setting %s\n',...
        responseSettings, outcomeSettings);

    if strcmp(responseSettings, 'Go')
        allResp = [1 0]; % Go/NoGo.
        
    elseif strcmp(responseSettings, 'Hand')
        allResp = [101 97 0]; % left/right/nogo.
    
    else
        fprintf('Unknown response setting')
    end
    
    % Initialize outcome settings:
    if strcmp(outcomeSettings, 'abs')
        allFb = [1 0 -1]; % reward/neutral/punishment
    
    elseif strcmp(outcomeSettings, 'rel')
        allFb = [1 2]; % preferred/non-preferred
    
    elseif strcmp(outcomeSettings, 'all')
        allFb = 1:4; % reward/no reward/no punishment/punishment
    
    else
        fprintf('Unknown response setting')
    end
end

% ----------------------------------------------------------------------- %
%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootDir, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));

% Set directories and parameters:
dirs        = set_dirs(rootDir);
par         = set_par();

dirs.TF     = fullfile(dirs.results, sprintf('TF_subjectlevel_%s%d', par.TF.TFtype, par.TF.nFreqs));

if stimERPcor 
    dirs.TF = sprintf('%s_outERPcor_%s_%s',dirs.TF, responseSettings, outcomeSettings); % overwrite target file name
end

if ~exist(dirs.TF, 'dir'); mkdir(dirs.TF); end

fprintf('Perform TF decomposition from %.03f to %.03f sec. \n',...
    par.TF.toi4tf(1), par.TF.toi4tf(end));

% ----------------------------------------------------------------------- %
%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.finalPP, '*.mat'));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n', nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList, 8, 10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames, subVec));

% ----------------------------------------------------------------------- %
%% Loop over subjects:

for iSub = selIndices % %iSub = 1;
    
    % ------------------------------------------------------------------- %
    %% Create output file name:
    
    outputFile      = fullfile(dirs.TF, sprintf('MRIpav_%03d_outcomelocked.mat', iSub));
    fprintf('EEG output file is %s\n', outputFile);
    
    % Check if output file already exists:
    if exist(outputFile, 'file')
        warning('Subject %03d: File %s already exists, skip subject\n', iSub, outputFile);
        continue;
    end

    % ------------------------------------------------------------------- %
    %% Create input file name:
    
    inputFile       = fullfile(dirs.finalPP, fileList{iSub}); % go to subject directory
    fprintf('EEG input file is %s\n', inputFile);
       
    % ------------------------------------------------------------------- %
    %% Load data from previous step:
    
    fprintf('Subject %03d: Start loading time-domain data\n', iSub);
    tmp     = load(inputFile);
    scd     = tmp.scd;
    hist    = tmp.hist;
    fprintf('Subject %03d: Finished loading time-domain data\n', iSub);

    % Filter out non-EEG channels:
    scd.elec.elecpos    = scd.elec.elecpos(par.chan.EEG, :); 
   
    % ------------------------------------------------------------------- %
    %% Load rejected trials:
    
    fprintf('Subject %03d: Load rejected trials\n', iSub)
    
    filenameRejTrials = fullfile(dirs.rejTrials, sprintf('MRIpav_%03d_rejTrials.mat', iSub));
    load(filenameRejTrials);

    % ------------------------------------------------------------------- %
    %% Load and store relevant behavioral data:
    
    fprintf('Subject %03d: Load behavior\n', iSub)

    behav = load_recode_behav(rootDir, iSub, rejectedtrials); 
        
    % ------------------------------------------------------------------- %
    %% Subtract stimulus-locked ERP:
    
    if stimERPcor
        
        % --------------------------------------------------------------- %
        % 1) Convert trials stored as cells into matrix:
        cfg = [];
        cfg.keeptrials = 'yes';
        scd = ft_timelockanalysis(cfg, scd);
        
        % --------------------------------------------------------------- %
        % 2) Baseline correction:
        % a) With Fieldtrip:
%         cfg = [];
%         cfg.baseline = par.TF.baselineTimings; % -250 - 50 ms
%         base_scd = ft_timelockbaseline(cfg, scd);

        % b) Regression based:
        startIdx    = dsearchn(scd.time', par.TF.baselinetime(1)); % sample index for start of baseline period
        stopIdx     = dsearchn(scd.time', par.TF.baselinetime(end)); % sample index for stop of baseline period
        scd.trial   = regress_baseline(scd.trial, scd.trialinfo, startIdx, stopIdx, 'subtraction');

%        scd_old     = scd; % save for later comparison
        
        % --------------------------------------------------------------- %
        % 3) Compute and subtract ERP per condition:
        iCondi = 0; % iCondi is stimulus number (response x outcome, i.e. 6)

        for iResp = allResp % iResp = 1; loop over Go/NoGo or GoLeft/GoRight/NoGo

            for iVal = allFb % iVal = 1; loop over preferred/non-preferred or reward/neutral/punishment
                    
                iCondi = iCondi + 1;
                   fprintf('Subject %03d: Correct for stimulus-locked ERP in iResp = %d, iVal = %d\n', iSub, iResp, iVal);

                % ----------------------------------------- %
                % a) Select trials based on response:
                if strcmp(responseSettings, 'Go')
                   respTrials = find(behav.isgo == iResp); % based on Go/NoGo
                elseif strcmp(responseSettings, 'Hand')
                   respTrials = find(behav.resp == iResp); % based on exact response
                elseif strcmp(responseSettings, 'none')
                   respTrials = behav.trlidx; % all trials
                else
                   error('Invalid response specification: %s', responseSettings);
                end
                   
                % ----------------------------------------- %
                % b) Select trials based on feedback:
                if strcmp(outcomeSettings, 'abs')
                   fbTrials = find(behav.fb.abs == iVal); % reward/neutral/punishment
                elseif strcmp(outcomeSettings, 'rel')
                   fbTrials = find(behav.fb.rel == iVal); % preferred/non-preferred
                elseif strcmp(outcomeSettings, 'all')
                   fbTrials = find(behav.fb.all == iVal); % reward/no-reward/no-punishment/punishment
                elseif strcmp(outcomeSettings, 'none')
                   fbTrials = behav.trlidx; % all trials
                else
                   error('Invalid outcome specification: %s', outcomeSettings);
                end 
                   
                % ----------------------------------------- %
                % c) Form intersection:
                selTrials    = intersect(respTrials, fbTrials);
                fprintf('Found %d trials\n', length(selTrials));

                % ----------------------------------------- %
                % d) Compute ERP for selected trials, subtract from selected trials:
                scd_cor = subtract_ERP(scd, scd, selTrials, selTrials, 'Fieldtrip', 0, 1.3); % entire data with Fieldtrip
                   
                % ----------------------------------------- %
                % e) Insert selectively corrected trials into original object:
                scd.trial(selTrials, :, :) = scd_cor.trial(selTrials, :, :);

            end % end of for loop for iVal
        end % end of for loop for iResp

       % 4) Do baseline-correction again:
       startIdx  = dsearchn(scd.time', parTF.baselinetime(1));
       stopIdx   = dsearchn(scd.time', parTF.baselinetime(end));
       scd.trial = regress_baseline(scd.trial, scd.trialinfo, startIdx, stopIdx, 'subtraction');

       % Clear objects created intermediately:
       clear respTrials fbTrials startIdx stopIdx
                
    end % end stimERPcor
    
    % ------------------------------------------------------------------- %
    %% Perform Morlet wavelet convolution.
    
    if strcmp(par.TF.TFtype, 'morlet')

        fprintf('Subject %03d: Start TF decomposition using Morlet wavelets\n', iSub);
        
        cfg                 = [];
        cfg.method          = 'wavelet';
        cfg.output          = 'pow'; %'fourier';% to also keep the phase.
        cfg.keeptrials      = 'yes';
%         cfg.foilim          = [parTF.minFreq parTF.maxFreq]; % either frequency min and max
        cfg.foi             = par.TF.freq4tf;      % or directly frequencies of interest.
        cfg.width           = par.TF.width4tf;     % number of cycles.
        cfg.toi             = par.TF.toi4tf;       % times (s) to centre the analysis window.
        cfg.channel         = 'all';
        freq                = ft_freqanalysis(cfg, scd);
        
    % ------------------------------------------------------------------- %
    %% Perform multi-taper-method convolution with Hanning tapers:
    
    elseif strmcp(par.TF.TFtype, 'hanning')
        
        fprintf('Subject %03d: Start TF decomposition using Hanning tapers\n', iSub);
        
        cfg                 = [];
        cfg.method          = 'mtmconvol'; % mtmfft mtmconvol
        cfg.taper           = 'hanning';
%         cfg.foilim          = [parTF.minFreq parTF.maxFreq]; % either frequency min and max        
        cfg.foi             = par.TF.freq4tf; % or directly frequencies of interest.
        cfg.t_ftimwin       = repmat(par.TF.ftimwin,1,length(cfg.foi)); % length of time window for each frequency
        cfg.toi             = par.TF.toi4tf; % time bins for which to compute frequency
        cfg.pad             = par.TF.pad; % pad up to 8 sec. to get well-behaved frequency bins
        cfg.output          = 'pow'; %'fourier' to also keep the phase.
        cfg.keeptrials      = 'yes';
        cfg.channel         = 'all';
        freq                = ft_freqanalysis(cfg, scd);

    end
    
    % Check whether input frequencies are output frequencies:
    if(sum(cfg.foi == freq.freq) ~= length(freq.freq))
        warning('Input and output frequencies do not match')
        fprintf(' Input frequencies: %s\n', strjoin(string(cfg.foi)));
        fprintf('Output frequencies: %s\n', strjoin(string(freq.freq)));
    end

    % ------------------------------------------------------------------- %
    %% Save data.
    
    fprintf('Subject %03d: Start saving data\n', iSub);
    
    hist.TF.dirs       = dirs;
    hist.TF.par        = par;

    save(outputFile, 'freq', 'hist', '-v7.3')
    clear scd freq par hist behav 
    
    fprintf('Subject %03d: Finished saving data :-)\n', iSub); 
        
end % iSub-loop.

fprintf('Done :-)\n');

end % END OF FUNCTION.