function TF_extract_singleTrial(rootDir, subVec, maskType, band, channels, contrast)

% TF_extract_singleTrial(rootDir, subVec, maskType, band, channels, contrast, GLMID)
% 
% Compute trial-by-trial mean power within mask defined previously in
% two-condition permutation test (from EEGfMRIPav_grouplevel_TF_permutation.m).
% 
% INPUTS:
% rootDir           = string, root directory of project
% subVec            = vector of integers, numbers of subjects to process
% maskType 	    = string, whether mask comes from t-test based on EEG power ('ttest') or from fMRI-EEG correlations in TAfT ('TAfT').
% band              = string, band for which mask was created (featured in
% mask name), either 'alpha' (default) or 'theta'
% channels          = string, list of channels (without separators) for
% which mask was created (featured in mask name), default 'CzFCzCz'
% contrast          = string, contrast from which mask was created
% (featured in mask name), either 'Preferred' (default) or 'Action' or 'ActionLong' or a regressor used in TAfT.
%
% OUTPUTS:
% none, save vectors with trial-by-trial EEG power as .csv to disk.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

% ----------------------------------------------------------------------- %
%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootDir', 'var')
    rootDir = grouplevel_set_rootDir();
    fprintf('rootDir unspecified, assume %s\n', rootDir)
end

if ~exist('iSub', 'var')
    subVec  = 1:36; 
    fprintf('subVec unspecified, assume %s\n', strjoin(string(subVec), ', '));
end

if ~exist('maskType', 'var')
    band    = 'ttest'; % ttest, TAFT
end

if ~exist('contrast', 'var')
    contrast = 'Preferred'; % Preferred Action ActionLong GLM1StriatumConj
end

if ~exist('channels', 'var')
    channels = 'CzFCzFz'; % CzFCzFz
end

if ~exist('band', 'var')
    band    = 'beta'; % beta, lowalpha, theta, deltatheta
end

if ~exist('GLMID', 'var')
    GLMID   = '3A';
end

nSub        = length(subVec);

% ----------------------------------------------------------------------- %
%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootDir, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));

% Set directories and parameters:
dirs        = set_dirs(rootDir);
par         = set_par();

% Add TF directory based on TF settings:
dirs.TF     = fullfile(dirs.results, sprintf('TF_subjectlevel_%s%d', par.TF.TFtype, par.TF.nFreqs));

% fMRI (output):
dirs.fMRI   = fullfile(rootDir, 'Log', 'fMRI'); % subject-specific directory created later

% ----------------------------------------------------------------------- %
%% ALTERNATIVE A: Select mask from 2D conditions permutation test:

if strcmp(maskType, 'ttest')
    maskName        = sprintf('Mask_%s_%s_%s.mat', contrast, channels, band);
fprintf('Selected mask is %s\n', maskName);

% Masks used for EEG regressor in fMRI GLM:
% maskName    = sprintf('%s_%s_%s', job.contrastType, strjoin(sort(job.channels), ''), job.bandName);
% maskName    = 'Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta'; targetName  = 'PreferredFrontalDeltatheta'; GLMID = '3A';
% maskName    = 'Preferred_CzFCzFz_beta'; targetName  = 'PreferredMidfrontalBeta'; GLMID = '3B';
% maskName    = 'Action_CzFCzFz_lowalpha'; targetName  = 'ActionMidfrontalAlpha'; GLMID = '3C';

% Used for single-trial:
% maskName    = 'Preferred_CzFCzFz_beta';
% maskName    = 'Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta';
% maskName    = 'Action_CzFCzFz_lowalpha';

% ----------------------------------------------------------------------- %
%% ALTERNATIVE B: SELECT FROM TAFT:

%     maskType        = 'TAfT';
%     maskName        = 'Updatebias_CzFCzFz';

elseif strcmp(maskType, 'TAfT')
    
    maskName        = sprintf('Mask_%s_%s.mat', contrast, channels);
    fprintf('Selected mask is %s\n', maskName);
    
else
    error('Unknown mask type');
end
% ----------------------------------------------------------------------- %
%% Load mask:
    
fprintf('Load mask %s \n', maskName);   
if strcmp(maskType, 'ttest')
    dirs.mask           = fullfile(dirs.results, 'TF_Mask');
    load(fullfile(dirs.mask, sprintf('Mask_%s.mat', maskName)));
elseif strcmp(maskType, 'TAfT')
    dirs.mask           = fullfile(dirs.results, 'TAfT_Betas/TAfT_Masks');
    load(fullfile(dirs.mask, sprintf('Mask_%s.mat', maskName)));
else
    error('Unknown mask type');
end

% Plot:
contourf(-1.0:0.025:1.3, 1:33, squeeze(nanmean(freqMask, 1)), ...
    40,  'linestyle', 'none');
title(sprintf('Mask %s', maskName));
pause(3)
close gcf

% ----------------------------------------------------------------------- %
%% Locate single-subject EEG files:

tmp         = fullfile(dirs.TF, sprintf('*_outcomelocked.mat');
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n', nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList, 8, 10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames, subVec));

% ----------------------------------------------------------------------- %
%% Create target directory:

if ~exist('targetName', 'var'); targetName = maskName; fprintf('Use maskName as target file name\n'); end

% GLM 3A: % Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta becomes PreferredFrontalDeltatheta
% GLM 3B: Preferred_CzFCzFz_beta becomes PreferredMidfrontalBeta
% GLM 3C: Action_CzFCzFz_lowalpha becomes StrongWeakMidfrontalLowalpha

% Directory with all trial-by-trial EEG power files:
dirs.singleTrial    = fullfile(dirs.results, 'TF_singleTrial');
if ~exist(dirs.singleTrial, 'dir'); mkdir(dirs.singleTrial); end

% Initialize target directory:
dirs.target         = fullfile(dirs.singleTrial, maskName);
if ~exist(dirs.target, 'dir'); mkdir(dirs.target); end

% ----------------------------------------------------------------------- %
%% Loop over subjects: 

allImputeMean       = nan(nSub, par.nCond + 1); % initialize mean per condition (+ invalid trials) per subject
allImputeCount      = nan(nSub, par.nCond + 1); % initialize count of NaN per condition (+ invalid trials) per subject

for iSub = selIndices

    % ------------------------------------------------------------------- %
    %% Create output file name:
      
    outputFile  = fullfile(dir.target, sprintf('EEG_%s.txt', targetFile));
    fprintf('EEG output file is %s\n', outputFile);
    
    % Check if output file already exists:
    if exist(outputFile, 'file')
        warning('Subject %03d: File %s already exists, skip subject\n', iSub, outputFile);
        return;
    end

    % ------------------------------------------------------------------- %
    %% Create input file name:
    
    inputFile   = fullfile(dirs.TF, fileList{iSub}); % go to subject directory
    fprintf('EEG input file is %s\n', inputFile);

    % ------------------------------------------------------------------- %
    %% 1) Load EEG:
    
    fprintf('Subject %03d: Start loading TF data\n', iSub);
    tmp     = load(inputFile);
    freq    = tmp.freq;
    hist    = tmp.hist;
    fprintf('Subject %03d: Finished loading TF data\n', iSub);

    % ------------------------------------------------------------------- %
    %% 2) Load rejected trials: 
    
    fprintf('Subject %03d: Load rejected trials\n', iSub);
    
    filenameRejTrials = fullfile(dirs.rejTrials, sprintf('MRIpav_%03d_rejTrials.mat', iSub));
    load(filenameRejTrials);

    % ------------------------------------------------------------------- %
    %% 3) Load and store relevant behavioral data.
    
    fprintf('Subject %03d: Load behavior\n', iSub);

    behav = load_recode_behav(rootDir, iSub, rejectedtrials); 
    
    % ------------------------------------------------------------------- %
    %% 4) Loop over non-rejected valid trials and extract summary measure of mask:
    
    fprintf('Extract trial-by-trial mean TF power in mask for non-rejected valid trials\n');
    
    nNonrejTrial    = length(freq.trialinfo); % number non-rejected trials
    powSumRaw       = nan(nNonrejTrial, 1); % initialize vector for non-rejected trials
    for iTrial = 1:nNonrejTrial % iTrial = 1;
        trlPow              = squeeze(freq.powspctrm(iTrial, :, :, :)); % power for this trial
        if sum(size(trlPow) ~= size(freqMask)) > 0; error('Data and mask have different dimensions'); end
        trlPowMask          = trlPow .* freqMask; % multiply power of this trial with mask
        powSumRaw(iTrial)   = sum(trlPowMask(:)); % sum across dimensions
    end
%     plot(powSumRaw)
    
    % ------------------------------------------------------------------- %
    %% 5) Compute mean per action x outcome:
    
%   nNonrejTrial % number non-rejected trials in EEG
%   behav % non-rejected trials in behavior (selected above via load_recode_behav())
    if(nNonrejTrial ~= length(behav.go)); error('EEG and behavior have different numbers of trials'); end

    fprintf('Subject %03d: Compute condition means \n', iSub);
    imputeVec   = nan(par.nCond, 1);
    iCond = 0;
    for iResp = [1 0] % iResp = 1; loop over Go/NoGo
        for iVal = [1 2 3 4] % iVal = 1;
            iCond               = iCond + 1;
            selTrials = behav.go == iResp & behav.fb.all == iVal; % based on Go/NoGo and outcome
            imputeVec(iCond)    = mean(powSumRaw(selTrials));
        end
    end
    imputeVecMean               = mean(imputeVec); % mean of all conditions
    
    % ------------------------------------------------------------------- %
    %% 6) Interpolate rejected/ invalid trials:
    
    fprintf('Subject %03d: Interpolate trials \n', iSub);
    
    powSumImp   = nan(par.nTrial, 1); % initialize 640 x 1 vector
    valTrlCount = 1; % initialize valid trial count
    imputeCount = zeros(par.nCond + 1, 1); % count # imputations per condition + imputations of grand mean
    
    for iTrial  = 1:par.nTrial
        if ismember(iTrial, behav.recoded.trlIdx) % if valid, non-rejected trial: insert
            powSumImp(iTrial)   = powSumRaw(valTrlCount);
            valTrlCount         = valTrlCount + 1; % increase valid trial count
            
        else % retrieve condition of trial, impute:
            iCond               = 4 * (1-behav.go(iTrial)) + behav.fb.all(iTrial); % determine condition of this trials
            
            if isnan(iCond) % if outcome undetermined (invalid trials):
                powSumImp(iTrial)           = imputeVecMean; % impute overall mean
                imputeCount(par.nCond+1)    = imputeCount(par.nCond+1) + 1; % increase imputation count for invalid trials
                
            else % if outcome determined:
                powSumImp(iTrial)           = imputeVec(iCond); % retrieve value used for imputation
                imputeCount(iCond)          = imputeCount(iCond) + 1; % increase imputation count
            end
            
        end
    end
    
    allImputeCount(iSub, :)  = imputeCount;
    fprintf('Subject %03d: Imputed %s trials \n', iSub, strjoin(string(imputeCount)));
    
    % ------------------------------------------------------------------- %
    %% Check if still NaN left:

    if sum(isnan(powSumImp)) > 0; error('Found remaining NaN in vector'); end
    
    % ------------------------------------------------------------------- %
    %% 7) Demean (for fMRI):

    powSumImp      = powSumImp - mean(powSumImp);
    
    % ------------------------------------------------------------------- %
    %% 8) Save as csv file:
    
    fileName    = fullfile(dirs.target, sprintf('TF_singleTrial_%s_sub%03d.csv', maskName, iSub));
    csvwrite(fileName, powSumImp);
    
    clear freq tmp behav rejectedtrials trlPow powSum
    
end
fprintf('Finished :-)\n');

% ----------------------------------------------------------------------- %
%% Evaluate:

% Imputation per subject per condiition:
% allImputeCount

% Overall imputation per condition:
% mean(allImputeCount, 1)  %
% Beta: 18.7222   22.5000   15.1944   20.0833   11.8611    6.8889   12.8611    9.6111    1.1944

% Overall imputation per subject:
% mean(allImputeCount, 2)'

end % END OF FUNCTION.