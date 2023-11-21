% FigureS04.m

% Plots for Figure S04 in manuscript.
% The same approach is used for Supplementary Figure S04.
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

% Add Fieldtrip:
addpath /home/common/matlab/fieldtrip;
ft_defaults;

% ----------------------------------------------------------------------- %
%% FigS04A: TF left:

clear job data

% ----------------------------------------------------------------------- %
%%  Set contrast:

% Needs code from https://github.com/johalgermissen/Algermissen2021CerCor/,
% add under rootDir/Analyses/EEG_Scripts/CueLockedAnalyses.

addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel'));

job                     = []; % initialize empty job

dirs.EEG                = fullfile(dirs.root, 'Log', 'EEG', 'CueLockedResults');
dirs.TFgroup            = fullfile(dirs.EEG, 'TF_grouplevel');
job.dirs                = dirs; % add directories

% Data settings:
job.nSub                = 36; % necessary for validSubs before loading data
job.sub2exclude         = [11 12 15 23 25 26 30]; % exclusion in TAfT-outcome-locked

job.lockSettings        = 'resplocked'; % stimlocked resplocked

job.baselineSettings    = 'trend'; % 'all' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand'  
job.actionSettings      = 'exeAct'; % 'reqAct' or 'exeAct'
job.accSettings         = 'correct'; % 'bothAcc' or 'correct' or 'incorrect' or 'noAcc' or 'reqAct'

% Frequency settings:
job.TFtype              = 'hanning'; % hanning morlet
job.nFreq               = 15; % 15

% Channel settings:
job.chanArea            = 'midfrontal'; % 'midfrontal' or 'frontal' or 'central' or 'leftmotor' or 'rightmotor'

% Band settings:
job.band                = 'theta'; % 'theta' or 'alpha' or 'beta' or 'broad' 

% Contrast of interest:
job.contrastType        = 'Go'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'CongruentAccuracy' or 'IncongruentAccuracy'

% ERP-corrected:
job.stimERPcor          = false;
job.respERPcor          = false;

% Load and prepare data:
job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job,data); % prepare generic objects
job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job,data); % prepare objects specific to contrast

% ----------------------------------------------------------------------- %
%% Create TF plot:

set(0, 'DefaultFigureColormap', redblue(64)); nContour = 10;

% Settings:
zlim        = 0.5;
fontSize    = 32;
yScale      = 'log';

% Start figure:
figure('Position', [100 100 1200 800], 'Color', 'white'); hold on

% Contour plot:
contourf(data.mu.time, data.mu.freq, data.TF2plot, nContour, 'linestyle', 'none'); hold on

% Settings:
strcmp(job.lockSettings, 'resplocked')
    set(gca, 'xlim', [-0.5 0.9], 'ylim', [1 33], 'clim', [-1*zlim 1*zlim], 'yscale', yScale, ... % yscale log
        'xtick', -1:0.5:1, 'xtickLabel', {'-1, 000', '-500', 'Resp', '500', '1, 000'}, ...
        'ytick', [2 4 8 16 32], 'yscale', yScale, ...
        'fontsize', fontSize, 'Linewidth', 3) % -.25 1.3
plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', 3);

% Labels:
xlabel('Time (ms)', 'fontsize', fontSize, 'fontweight', 'bold');
ylabel('Frequency (Hz)', 'fontsize', fontSize, 'fontweight', 'bold');

% Color bar:
colorbar('Ticks', (-1*zlim):(zlim/2):zlim, 'Fontsize', fontSize);

% Save:
figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS04ALeft'));
saveas(gcf, [figName '.png'])
pause(1);
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04ALeft.csv'));
csvwrite(fullFileName, data.TF2plot);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% FigS04A: TFplot right:

clear job data

% ----------------------------------------------------------------------- %
%%  Set contrast:

% Add paths to helper files:
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel'));

job                     = [];

dirs.EEG                = fullfile(dirs.root, 'Log', 'EEG', 'OutcomeLockedResults');
dirs.TFgroup            = fullfile(dirs.EEG, 'TF_grouplevel');
job.dirs                = dirs;

job.nSub                = 36; % necessary for validSubs
job.sub2exclude         = [11 12 15 23 25 26 30]; % outliers in TAfT and fMRI
job.nFreqs              = 33; % 33 or 17
job.baselineSettings    = 'trend'; % 'grandAvg' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'all'; % 'rel' or 'abs' or 'all'

job.outERPcor           = false;

% ----------------------------------------------------- %
% Channels:
job.chanArea            = 'midfrontal'; % 'midfrontal'
% Band settings:
job.band                = 'lowalpha'; % 'theta' or 'alpha' or 'beta' or 'broad' % lowalpha middlebeta
% Contrast of interest:
job.contrastType        = 'Action'; % 'Preferred' or 'Action' or 


if ~exist('data', 'var')
    job         = TF_update_job(job);
    [job, data] = TF_load_data(job);
    [job, data] = TF_prepare_generic_data(job, data);
end

job             = TF_update_job(job);
[~, data]       = TF_prepare_contrast_data(job, data);

% ----------------------------------------------------------------------- %
%% Create TF plot:

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)) ; cmap = 'redblue'; nContour = 10; % default 10

% Settings:
xLim        = [-0.9 0.8]; 
zlim        = 0.5;
fontSize    = 32;
yScale      = 'log';

% Plot:
figure('Position', [100 100 1439 800], 'Color', 'white'); hold on

% Contour plot: 
contourf(data.mu.time, data.mu.freq, data.TF2plot, nContour, 'linestyle', 'none'); hold on

% Vertical lines:
plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', 3);
plot([0.7 0.7], get(gca, 'ylim'), ':k', 'LineWidth', 3);

% Settings:
set(gca, 'xlim', xLim, 'ylim', [1 33], 'clim', [-1*zlim 1*zlim], ...
    'xtick', [-1 -0.50 -0.25 0 0.25 0.5 0.7 1], ...
    'xtickLabel', {'-1, 000', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
    'ytick', [2 4 8 16 32], 'yscale', yScale, ...
    'fontsize', fontSize, 'Linewidth', 3);

% Labels:
ylabel('Frequency (Hz)', 'fontsize', fontSize, 'fontweight', 'bold');
xlabel('Time (ms)', 'fontsize', fontSize, 'fontweight', 'bold');

% Color bar:
colorbar('Ticks', (-1*zlim):(zlim/2):zlim, 'Fontsize', fontSize);

% Save:
figName = fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS04ARight');
saveas(gcf, [figName '.png']);

% Close:
pause(3)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04ARight.csv'));
csvwrite(fullFileName, data.TF2plot);

% ----------------------------------------------------------------------- %
%% Fig5B: Line plot:

clear job data

% ----------------------------------------------------------------------- %
%%  Set contrast:

job.dirs                = dirs;
job.nSub                = 36; % necessary for validSubs
job.sub2exclude         = [11 12 15 23 25 26 30]; % outliers in TAfT and fMRI
job.nFreqs              = 33; % 33 or 17
job.baselineSettings    = 'trend'; % 'grandAvg' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'all'; % 'rel' or 'abs' or 'all'

job.outERPcor           = false;

% ----------------------------------------------------- %
% Channels:
job.chanArea            = 'midfrontal'; % 'midfrontal'
% Band settings:
job.band                = 'lowalpha'; % 'theta' or 'alpha' or 'beta' or 'broad' % lowalpha middlebeta
% Contrast of interest:
job.contrastType        = 'Action'; % 'Preferred' or 'Action' or 


if ~exist('data', 'var')
    job         = TF_update_job(job);
    [job, data] = TF_load_data(job);
    [job, data] = TF_prepare_generic_data(job, data);
end

job             = TF_update_job(job);
[job, data]     = TF_prepare_contrast_data(job, data);

% ----------------------------------------------------------------------- %
%% Create lineplot:

% New color map for red-green blind:
job.colMat  = [0 113 116; 87 196 173; 240 174 102; 201 61 33] ./ 255;
job.colMat  = repmat(job.colMat, job.nCond / size(job.colMat, 1), 1);

% Baselines:
iBaseline   = find(round(data.mu.time, 3) == 0); % find baseline via time in sec
baseline    = min(squeeze(nanmean(data.SubCondTime(job.validSubs, :, iBaseline)))); % still extract only valid subjects
until       = find(round(data.mu.time, 3) == 1); % end for outcome-locked

% Other settings:
xLim        = [-0.9 0.8];
lineWidth   = 5;
fontSize    = 32;
transp      = 0.10;

% Start plot:
figure('Position', [100 100 1200 800], 'color', 'white'); hold on

% Loop over conditions, make bounded line plots:
p   = cell(job.nCond, 1);
for iCond = 1:job.nCond

    p{iCond} = boundedline(data.mu.time(1:until), ...
        squeeze(nanmean(data.SubCondTime(job.validSubs, iCond, 1:until))) - baseline, ...
        job.nCond / (job.nCond-1) * ...
        squeeze(nanstd(squeeze(data.SubCondTime(job.validSubs, iCond, 1:until)) - ...
        data.SubTime(job.validSubs, 1:until) + ...
        repmat(data.GrandTime(1:until), length(job.validSubs), 1)))' ./ ...
        sqrt(length(job.validSubs)), ...
        'cmap', job.colMat(iCond, :) , 'alpha', 'transparency', transp); % , 'linewidth', 2); hold on
    set(p{iCond}, 'linestyle', job.lineStyle{iCond});
    set(p{iCond}, 'Linewidth', lineWidth);

end

% Axis labels and limits:
xlabel('Time (ms)', 'fontweight', 'bold', 'fontsize', fontSize);
yLim = get(gca, 'ylim');
yMinLim = yLim(1); 
yMaxLim = yLim(2); 

ylabel('Power (dB)', 'fontsize', fontSize, 'fontweight', 'bold')
yMinLim = yMinLim - 0.1;
yMaxLim = yMaxLim + 0.1;

set(gca, 'xlim', xLim, 'ylim', [yMinLim yMaxLim], ...
    'xtick', [-1 -0.75 -0.50 -0.25 0 0.25 0.5 0.7 1 1.25], ...
    'xtickLabel', {'-1, 000', '-750', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1, 000'}, ...
    'fontsize', fontSize, 'Linewidth', 3);
set(gca, 'TickLength', [0 0]);

% Gray patch:
yLim = get(gca, 'ylim');
patch([0 .7 .7 0], [yLim(1) yLim(1) yLim(2) yLim(2)], [0.8 0.8 0.8], 'facealpha', 0.2, 'EdgeColor', 'none');

% Save:
figName = fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS04B');
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data file:
% See job.condNames
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04B_GoReward.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 1, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04B_GoNoReward.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 2, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04B_GoNoPunishment.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 3, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04B_GoPunishment.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 4, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04B_NoGoReward.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 5, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04B_NoGoNoReward.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 6, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04B_NoGoNoPunishment.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 7, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04B_NoGoPunishment.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 8, 1:until));

% ----------------------------------------------------------------------- %
%% FigS04C: Brain slices:

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.GLMID           = '2without7';

% ---------------------------------------------- %
% Sums of both hands during response:
job.iCope = 4; job.iView = 'coronal'; job.iSlice = 10; job.zLim = [0 5]; job.cLim = 150; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS04CLeft.png'));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04CLeft_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04CLeft_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% ---------------------------------------------- %
% Difference of both hands during response:
job.iCope = 5; job.iView = 'coronal'; job.iSlice = -30; job.zLim = [0 5]; job.cLim = 150; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS04CRight.png'));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04CRight_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04CRight_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% FigS04D: Brain slices:

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.GLMID           = '2without7';

% ---------------------------------------------- %
% Action (Go - NoGo) during outcome:
job.iCope = 1; job.iView = 'coronal'; job.iSlice = 10; job.zLim = [0 5]; job.cLim = 500; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS04DLeft.png'));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04DLeft_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04DLeft_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

job.iCope = 1; job.iView = 'coronal'; job.iSlice = -30; job.zLim = [0 5]; job.cLim = 500; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS04DRight.png'));
close gcf


% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04DRight_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS04DRight_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% END OF FILE.