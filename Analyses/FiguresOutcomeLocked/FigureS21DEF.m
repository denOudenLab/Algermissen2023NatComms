% FigureS21DEF.m

% Plots for Supplementary Figure S21, panels D, E, and F.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresOutcomeLocked

% clear all; close all; clc

% Set root directory:
cd(fileparts(which('figures_set_rootDir.m')));
dirs.root           = figures_set_rootDir();

% ----------------------------------------------------------------------- %
%% Figure S21D: Go/NoGo per trial per subject

% Add helper functions:
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers/')); % for load_recode_behav()

% Load complete directories:
dirs            = set_dirs(dirs.root);

% Fixed settings:
nSub            = 36;
nTrial          = 640;

% ----------------------------------------------------------------------- %
% Loop over subject, extract Go/NoGo responses:

pGo     = nan(nSub, nTrial);

for iSub = 1:nSub % iSub = 1;

    fprintf('Start subject %03d\n', iSub);

    % Load behavior:
    behav                   = load_recode_behav(dirs.root, iSub);
    
    % Store per subject:
    pGo(iSub, behav.trlIdx) = behav.go;
    
end

% ----------------------------------------------------------------------- %
% Plot with imagesc plot:

% Settings:
fontSize    = 32; % 32
set(0, 'DefaultFigureColormap', 1 - gray(2)); % invert gray color matrix

% Figure:
close all
f = figure('Position', [100 100 1200 800], 'Color', 'white'); hold on
imagesc(1:nTrial, 1:nSub, pGo);

% Vertical lines:
plot([nTrial/2+1 nTrial/2+1], get(gca, 'ylim'), '-b', 'LineWidth', 3);

% Plot settings:
set(gca, 'xlim', [1 nTrial], 'ylim', [1 nSub], 'clim', [0 1], ...
    'xtick', 0:100:700, 'ytick', 1:5:40,...
    'fontsize', fontSize, 'Linewidth', 3);

% Labels:
xlabel('Trial number', 'fontsize', fontSize);
ylabel('Participant number', 'fontsize', fontSize);

% Save:
saveas(gcf, sprintf('/project/3017042.02/Log/OutcomeLockedPaperPlots/FigS21D.png'))
pause(1)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21D.csv'));
csvwrite(fullFileName, pGo);

% ----------------------------------------------------------------------- %
%% Figure S21E: Upsampled BOLD and plot time course per action per outcome:

% Add paths:
addpath '/home/common/matlab/spm12' % for spm_pinv
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers')); % for boundedline()
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT'));

% Select ROI:
ROIs2use    = {'GLM1ACCConjMan'};
iROI        = 1;

% Data dimensions:
nSub        = 36;
nBlock      = 6;
nCond       = 4; % action x PE valence

% Initialize dirs for TAfT:
dirs.log    = fullfile(dirs.root, 'Log');
dirs.behav  = fullfile(dirs.log, 'Behavior/Data_beh_mat');
dirs.EEG    = fullfile(dirs.log, 'EEG/OutcomeLockedResults');    
dirs.fMRI   = fullfile(dirs.log, 'fMRI');    

% Initialize settings for TAfT:
job.hpass       = 128;
job.ups         = 10;
job.ons_unit    = 1;
job.TR          = 1.4;
job.trialdur    = 14; % 16 too long

% ----------------------------------------------------------------------- %
% Loop over subjects, extract fMRI and behavior:

subfMRI     = cell(nSub, 1); % initialize
subBehav    = cell(nSub, 1); % initialize

for iSub = 1:nSub
    
    data        = []; % initialize object to store data in
    
    % Load fMRI:   
    for iBlock = 1:nBlock % for each block % iBlock = 6
        
        % Name of ROIs just carried forward:
        job.ROIs(iROI).ROIname                      = ROIs2use{iROI}; 
        
        % Where BOLD time series, onsets and nuisance regressors are expected to be:
        fmriBlockDir                                = fullfile(dirs.fMRI, ...
            sprintf('sub-%03d/FEAT_Block%d.feat/AROMA', iSub, iBlock));
        
        % BOLD time series extracted from ROI:
        job.ROIs(iROI).ROIdef(iBlock).rawfMRIfile   = {fullfile(fmriBlockDir, ...
            sprintf('%s.txt', job.ROIs(iROI).ROIname))}; % where to load volumes from

        % Trial onsets (text file):
        job.ROIs(iROI).ROIdef(iBlock).onsets        = {fullfile(fmriBlockDir, ...
            'outcomeOnsets.txt')}; % trial onsets    
        
        % Upsample ROI-data:
        fprintf('Subject %03d Block %d: Upsample data\n', iSub, iBlock);
        [~, blockData, ~, ~]    = taft_preprocess_filter_upsample_epoch(job, iROI, iBlock); % dat is of format trials x time bins
        data                    = [data; blockData]; % append data for this block
        % data is of format # trials x # time points
        
    end
    % Save entire fMRI data per subject:
    subfMRI{iSub}           = data;
    
    % Load and store behavior per subject:
    job.behavFile           = fullfile(dirs.root, '/Log/Behavior/Data_beh_mat', ...
    sprintf('3017042.02_emmvdij_%03d_001_results.mat', iSub)); 
    subBehav{iSub}          = taft_preprocess_load_behavior(job);
    
end

% ----------------------------------------------------------------------- %
% Select trials based on behavior:

SubCondMean     = nan(nSub, size(subfMRI{1}, 2), nCond);
SubCondTrial    = nan(nSub, nCond); % initialize

for iSub = 1:nSub
    data    = subfMRI{iSub};
    iCond   = 1;
    for iGo = [1 0]
        for iVal = [1 2]
            selIdx                      = find(...
                subBehav{iSub}.isgo == iGo & ...
                subBehav{iSub}.fb.rel == iVal);
            SubCondMean(iSub, :, iCond)     = nanmean(data(selIdx, :), 1);
            SubCondTrial(iSub, iCond)       = length(selIdx);
            iCond = iCond + 1;
        end
    end
    clear data
end

% ----------------------------------------------------------------------- %
% Average over subjects:

% Overall mean per subject (for SEs):
SubMean                 = nan(nSub, size(SubCondMean, 2));
for iSub                = 1:nSub
    SubMean(iSub, :)    = nanmean(SubCondMean(iSub, :, :), 3);
end

% Grand mean across subjects (for SEs):
invalidSubs         = [15 25]; % subjects with bad co-registrations to start with
validSubs           = setdiff(1:nSub, invalidSubs);
GrandMean           = nanmean(SubMean(validSubs, :), 1);

% ----------------------------------------------------------------------- %
% Start figure:

% Plotting settings:
colMat      = [0 113 116; 201 61 33] ./ 255; % red-green-blind-friendly colors
colMat      = repmat(colMat, nCond / size(colMat, 1), 1);
fontSize    = 32;

p           = cell(nCond, 1);
figure('Position', [100 100 1000 800], 'Color', 'white'); hold on
for iCond = 1:nCond
    p{iCond} = boundedline((1:size(SubCondMean, 2)) * job.TR/job.ups, ...
        squeeze(nanmean(SubCondMean(validSubs, :, iCond))), ...
        squeeze(nCond / (nCond-1) * nanstd(squeeze(SubCondMean(validSubs, :, iCond)) - ...
        SubMean(validSubs, :) + repmat(GrandMean, length(validSubs), 1)) ./ ...
        sqrt(length(validSubs))), ...
        'cmap', colMat(iCond, :), 'alpha');
    if iCond > 2
        set(p{iCond}, 'linestyle', '--');
    end
    set(p{iCond}, 'linewidth', 4);
end
set(gca, 'xlim', [0 job.trialdur], 'ylim', [-0.15 0.15], ...
    'xtick',0:4:40, 'ytick',-1:.1:1, ...
    'FontSize', fontSize, 'FontName', 'Arial', 'Linewidth', 4);
ylabel('BOLD signal (a.u.)', 'fontsize', fontSize);
xlabel('Time (s)', 'fontsize', fontSize);

% ----------------------------------------------------------------------- %
% Save:

fprintf('Save\n')
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS21E.png'));
pause(3);
close gcf
fprintf('Finished :-)\n')

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21E_GoPositive.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 1));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21E_GoNegative.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 2));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21E_NoGoPositive.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 3));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21E_NoGoNegative.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 4));

% ----------------------------------------------------------------------- %
%% Figure S21F: Upsampled BOLD and plot time course per action per outcome split by next action:

% Add paths:
addpath '/home/common/matlab/spm12' % for spm_pinv
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers')); % for boundedline()
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT'));

% Select ROI:
ROIs2use    = {'GLM1ACCConjMan'};
iROI        = 1;

% Data dimensions:
nSub        = 36;
nBlock      = 6;
nCond       = 8; % action x PE valence x previous action !!!!!!!!!!!!!!!!!

% Initialize dirs for TAfT:
dirs.log    = fullfile(dirs.root, 'Log');
dirs.behav  = fullfile(dirs.log, 'Behavior/Data_beh_mat');
dirs.EEG    = fullfile(dirs.log, 'EEG/OutcomeLockedResults');    
dirs.fMRI   = fullfile(dirs.log, 'fMRI');    

% Initialize settings for TAfT:
job.hpass       = 128;
job.ups         = 10;
job.ons_unit    = 1;
job.TR          = 1.4;
job.trialdur    = 14; % 16 too long

% ----------------------------------------------------------------------- %
% Loop over subjects, extract fMRI and behavior:

subfMRI     = cell(nSub, 1); % initialize
subBehav    = cell(nSub, 1); % initialize

for iSub = 1:nSub
    
    data        = []; % initialize object to store data in
    
    % Load fMRI:   
    for iBlock = 1:nBlock % for each block % iBlock = 6
        
        % Name of ROIs just carried forward:
        job.ROIs(iROI).ROIname                      = ROIs2use{iROI}; 
        
        % Where BOLD time series, onsets and nuisance regressors are expected to be:
        fmriBlockDir                                = fullfile(dirs.fMRI, ...
            sprintf('sub-%03d/FEAT_Block%d.feat/AROMA', iSub, iBlock));
        
        % BOLD time series extracted from ROI:
        job.ROIs(iROI).ROIdef(iBlock).rawfMRIfile   = {fullfile(fmriBlockDir, ...
            sprintf('%s.txt', job.ROIs(iROI).ROIname))}; % where to load volumes from

        % Trial onsets (text file):
        job.ROIs(iROI).ROIdef(iBlock).onsets        = {fullfile(fmriBlockDir, ...
            'outcomeOnsets.txt')}; % trial onsets    
        
        % Upsample ROI-data:
        fprintf('Subject %03d Block %d: Upsample data\n', iSub, iBlock);
        [~, blockData, ~, ~]    = taft_preprocess_filter_upsample_epoch(job, iROI, iBlock); % dat is of format trials x time bins
        data                    = [data; blockData]; % append data for this block
        % data is of format # trials x # time points
        
    end
    % Save entire fMRI data per subject:
    subfMRI{iSub}           = data;
    
    % Load and store behavior per subject:
    job.behavFile           = fullfile(dirs.root, '/Log/Behavior/Data_beh_mat', ...
    sprintf('3017042.02_emmvdij_%03d_001_results.mat', iSub)); 
    subBehav{iSub}          = taft_preprocess_load_behavior(job);
    
end

% ----------------------------------------------------------------------- %
% Select trials based on behavior:

SubCondMean     = nan(nSub, size(subfMRI{1}, 2), nCond);
SubCondTrial    = nan(nSub, nCond); % initialize

for iSub = 1:nSub
    data    = subfMRI{iSub}; % retrieve data for this subject
    iCond   = 1; % initialize condition index
    for iNextGo = [1 0]
        for iGo = [1 0]
            for iVal = [1 2]
                selIdx  = find(...
                    subBehav{iSub}.isnextgo == iNextGo & ...
                    subBehav{iSub}.isgo == iGo & ...
                    subBehav{iSub}.fb.rel == iVal);
                SubCondMean(iSub, :, iCond)   = nanmean(data(selIdx, :), 1);
                SubCondTrial(iSub, iCond)    = length(selIdx);
                iCond = iCond + 1;
            end
        end
    end
    clear data
end

% ----------------------------------------------------------------------- %
% Average over subjects:

% Overall mean per subject (for SEs):
SubMean                 = nan(nSub, size(SubCondMean, 2));
for iSub                = 1:nSub
    SubMean(iSub, :)    = nanmean(SubCondMean(iSub, :, :), 3);
end

% Grand mean across subjects (for SEs):
invalidSubs         = [15 25]; % subjects with bad co-registrations to start with
validSubs           = setdiff(1:nSub, invalidSubs);
GrandMean           = nanmean(SubMean(validSubs, :), 1);

% ----------------------------------------------------------------------- %
% Start figure:

% Fixed settings:
colMat      = [0 113 116; 201 61 33] ./ 255; % red-green-blind-friendly colors
colMat      = repmat(colMat, nCond / size(colMat, 1), 1);
fontSize    = 32;

p           = cell(nCond, 1);
figure('Position', [100 100 1200 800], 'Color', 'white'); hold on

% Next action is Go: left plot
subplot(1, 2, 1);
for iCond = 1:nCond/2
    p{iCond} = boundedline((1:size(SubCondMean, 2)) * job.TR/job.ups, ...
        squeeze(nanmean(SubCondMean(validSubs, :, iCond))), ...
        squeeze(nCond / (nCond-1) * nanstd(squeeze(SubCondMean(validSubs, :, iCond)) - ...
        SubMean(validSubs, :) + repmat(GrandMean, length(validSubs), 1)) ./ ...
        sqrt(length(validSubs))), ...
        'cmap', colMat(iCond, :), 'alpha');
    if iCond > 2
        set(p{iCond}, 'linestyle', '--');
    end
    set(p{iCond}, 'linewidth', 4);
end

% Plot settings:
set(gca, 'xlim', [0 job.trialdur], 'ylim', [-0.15 0.15], ...
    'xtick',0:4:40, 'ytick',-1:.1:1, ...
    'FontSize', fontSize, 'FontName', 'Arial', 'Linewidth', 4);
ylabel('BOLD signal (a.u.)', 'fontsize',fontSize);
xlabel('Time (s)', 'fontsize', fontSize);

% Next action is NoGo: right plot
subplot(1, 2, 2)
for iCond = (nCond/2+1):nCond
    p{iCond} = boundedline((1:size(SubCondMean, 2)) * job.TR/job.ups, ...
        squeeze(nanmean(SubCondMean(validSubs, :, iCond))), ...
        squeeze(nCond / (nCond-1) * nanstd(squeeze(SubCondMean(validSubs, :, iCond)) - ...
        SubMean(validSubs, :) + repmat(GrandMean, length(validSubs), 1)) ./ ...
        sqrt(length(validSubs))), 'cmap', colMat(iCond - nCond/2, :), 'alpha');
    if iCond > 6
        set(p{iCond}, 'linestyle', '--');
    end
    set(p{iCond}, 'linewidth', 4);
end

% Plot settings:
set(gca, 'xlim', [0 job.trialdur], 'ylim', [-0.15 0.15], ...
    'xtick',0:4:40, 'ytick',-5:2:5, 'FontSize', ...
    fontSize, 'FontName', 'Arial', 'Linewidth', 4);
xlabel('Time (s)', 'fontsize', fontSize);

% ----------------------------------------------------------------------- %
% Save:

fprintf('Save\n')
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS21F.png'));
pause(3);
close gcf
fprintf('Finished :-)\n')

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21F_GoPositiveNextGo.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 1));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21F_GoNegativeNextGo.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 2));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21F_NoGoPositiveNextGo.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 3));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21F_NoGoNegativeNextGo.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 4));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21F_GoPositiveNextNoGo.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 5));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21F_GoNegativeNextNoGo.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 6));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21F_NoGoPositiveNextNoGo.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 7));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS21F_NoGoNegativeNextNoGo.csv'));
csvwrite(fullFileName, SubCondMean(validSubs, :, 8));

% END OF FILE.