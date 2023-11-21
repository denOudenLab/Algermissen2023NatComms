% FigureS17.m

% Plots for Supplementary Figure S17.
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

addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT/'));

% Add Fieldtrip:
fprintf('Add Fieldtrip\n');
addpath /home/common/matlab/fieldtrip;
ft_defaults;

% Remove SPM (if necessary):
fprintf('Remove SPM\n');
rmpath /home/common/matlab/spm12;

% ----------------------------------------------------------------------- %
%% Figure S17A:

% This figure displays the 7 regions in the conjunction of PE_STD and
% PE_DIF in GLM1. 

% ----------------------------------------------------------------------- %
%% General settings for TAfT:

% 0) EEG domain:
EEGdomain = 'TF';

% ----------------------------------------------------------------------- %
% 1) ROI to use:

ROIs2use    = {'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1LeftMotorConj', ...
    'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1LeftITGConj', 'GLM1V1Conj'};
ROInames    = {'Striatum', 'ACC', 'LeftMotor', 'vmPFC', 'PCC', 'LeftITG', 'V1'};
ROIColors   = [190 190 88; 255 0 0; 255 161 123; 51 43 170; 1 176 240; 74 203 145; 152 252 0] / 255;
nROI        = length(ROIs2use);

% ----------------------------------------------------------------------- %
% 2) Behavioral regressors:

behav2use = {'Updatestd', 'Updatedif'}; % --> goes into paper

% ----------------------------------------------------------------------- %
% 3) Perform only on selected trials

selTrials = 'all';

% ----------------------------------------------------------------------- %
% Data selection and plotting settings:

iSub        = 1;
iBlock      = 1;

fontSize    = 18;
LWD         = 4;
markerSize  = 12;

% ----------------------------------------------------------------------- %
% Initialize TAfT job:

job         = taft_preprocess_initialize_job(EEGdomain, iSub, ROIs2use, behav2use, selTrials);
    
% ----------------------------------------------------------------------- %
%% Figure S17B: Load volume-by-volume raw data of one block per ROI, plot:

% Set number of "dots" in time line:
xMax        = 40;
xMaxTime    = xMax * job.ups;

% Add lines (panel C) or not (panel B):
% addLine     = false;
addLine     = true;

% Initialize matrix for saving:
saveMat     = nan(47, nROI); 

for iROI = 1:nROI

    % job.ROIs(iROI).ROIdef(iBlock).rawfMRIfile
    [job, data, data_ups, i_ons] = taft_preprocess_filter_upsample_epoch(job, iROI, iBlock);

    close all

    figure('Position', [100 100 1200 150], 'Color', 'white'); hold on

    xVec        = 0:job.TR_ups:(job.TR_ups*(length(data_ups) - 1));
    pointIdx    = 1:job.ups:(length(data_ups) / job.ups  - 1);
    yVec        = data_ups;

    % Plot line:
    plot(xVec(1:xMaxTime), yVec(1:xMaxTime), '-', 'color', ROIColors(iROI, :), 'linewidth', LWD);

    % Plot points:
    plot(xVec(pointIdx), yVec(pointIdx), 'o', 'color', ROIColors(iROI, :), 'linewidth', LWD, 'markersize', markerSize);
    saveMat(:, iROI)    = yVec(pointIdx); % save
    
    if addLine
        timeVec = job.trialdur:job.trialdur:(job.trialdur*xMax);
        yLim    = get(gca, 'ylim');
        for iTime = 1:xMax
            plot([timeVec(iTime) timeVec(iTime)], yLim, 'k--', 'linewidth', 2);
        end
    end

    % Plotting settings:
    set(gca, 'xlim', [0 xMaxTime*job.TR_ups]);
    set(gca, 'linewidth', LWD, 'fontsize', fontSize);
    xlabel('Time (in sec.)');
    ylabel('BOLD');

    % Save:
    figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS17B_%s', ROInames{iROI}));
    saveas(gcf, [figName '.png']);
    pause(1);
    close gcf

end

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS17B.csv');
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S17C, top row: 

% recycle Figure #2 from Panel B.

% ----------------------------------------------------------------------- %
%% Figure S17C, middle row: Load SPM HRF, plot:

close all

hrf     = spm_hrf(job.TR_ups); 

figure('Position', [100 100 200 200]); hold on
plot(1:length(hrf), hrf, 'k-', 'linewidth', LWD);
set(gca, 'xtick', [], 'ytick', []);
set(gca, 'linewidth', LWD);

figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS17C_HRF');
saveas(gcf, [figName '.png']);
pause(1);
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS17C_HRF.csv');
csvwrite(fullFileName, hrf);

% ----------------------------------------------------------------------- %
%% Figure S17C, bottom row: Load trial-by-trial HRF amplitude, plot points horizontally:

xMax        = 7; % number data points
saveMat     = nan(xMax, nROI); % initialize for saving

for iROI = 1:nROI

    [job, data, data_ups, i_ons]    = taft_preprocess_filter_upsample_epoch(job, iROI, iBlock);

    amplitudes                      = taft_preprocess_fit_HRF_trial(job, data);

    close all

    figure('Position', [100 100 1200 150]); hold on

    xVec        = (1:length(amplitudes)) - 0.5;

    % Plot points:
    plot(xVec(1:xMax), amplitudes(1:xMax), ...
        'o', 'color', ROIColors(iROI, :), 'linewidth', LWD, 'markersize', markerSize);
    saveMat(:, iROI)    = amplitudes(1:xMax); % save

    set(gca, 'xlim', [0 xMax]);
    yLim = get(gca, 'ylim');
    yLim(1) = yLim(1) - 10;
    set(gca, 'ylim', yLim);
    set(gca, 'ytick', []);
    set(gca, 'linewidth', LWD, 'fontsize', fontSize);
    xlabel('Trial');
    ylabel('beta');

    figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS17C_bottom_%s', ROInames{iROI}));
    saveas(gcf, [figName '.png']);
    pause(1);
    close gcf
    
end

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS17C_bottom.csv');
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S17D: Load trial-by-trial HRF amplitude, plot lines vertically:

xMax        = 50; % number data points
saveMat     = nan(xMax, nROI); % initialize for saving

for iROI    = 1:7

    [job, data, data_ups, i_ons]    = taft_preprocess_filter_upsample_epoch(job, iROI, iBlock);

    amplitudes                      = taft_preprocess_fit_HRF_trial(job, data);

    close all

    figure('Position', [100 100 100 600], 'Color', 'white'); hold on

    xVec        = (1:length(amplitudes)) - 0.5;
    yVec        = amplitudes;

    % Plot points:
    plot(amplitudes(1:xMax), xVec(1:xMax), ...
        '-', 'color', ROIColors(iROI, :), 'linewidth', LWD);
    saveMat(:, iROI)    = amplitudes(1:xMax); % save

    set(gca, 'visible', 'off')

    figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS17D_%s', ROInames{iROI}));
    saveas(gcf, [figName '.png']);
    %     pause(1);
    close gcf
    
end

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS17D.csv');
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Load EEG data:

fprintf('Start loading EEG data ... \n');
inputEEG        = load(job.EEGfile);
fprintf('Finished loading EEG data! :-)\n');

% Average over channels:

cfg             = [];
cfg.channels    = 'FCz';
cfg.avgoverchan = 'yes'; 
EEGData         = ft_selectdata(cfg, inputEEG.freq); % extract frequency data
fprintf('Selected channel %s\n', cfg.channels);

% ----------------------------------------------------------------------- %
%% Figure S17E: Plot EEG per trial:

addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)) ; % default 10

nTrial  = 5; % number trials
selFreq = 1:15;  % frequency range

close all

for iTrial = 1:nTrial
    
    figure('Position', [100 100 200 200]); hold on
    contourf(EEGData.time, EEGData.freq(selFreq), squeeze(EEGData.powspctrm(iTrial, :, selFreq, :)));

    set(gca, 'xtick', [], 'ytick', [], 'linewidth', 3);
    xlabel('Time');
    ylabel('Frequency');
    title(sprintf('FCz, Trial %d', iTrial), 'fontweight', 'normal');
    
    % Save:
    figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS17E_trial%d', iTrial));
    saveas(gcf, [figName '.png']);
    pause(1);
    close gcf

    % Save source data file:
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS17E_trial%d.csv', iTrial));
    csvwrite(fullFileName, squeeze(EEGData.powspctrm(iTrial, :, selFreq, :)));
    
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Load TAfT:

EEGdomain       = 'TF';
ROIs2use        = {'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1LeftMotorConj', ...
    'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1LeftITGConj', 'GLM1V1Conj'};
behav2use       = {'Updatestd', 'Updatedif'};
selTrials       = 'all';

[job, dirs, betas] = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
%% FigureS17F, top row: Plot individual subjects' TF data:

selFreq = 1:15;     % frequency range
selChan = 'FCz';    % find(strcmp(betas{2}.label, 'FCz')) % selected channel
iROI    = 2;
nSub    = 3;

for iSub = 1:nSub
    
    iChan = find(strcmp(betas{2}.label, selChan));
    
    figure('Position', [100 100 200 200]); hold on
    contourf(betas{iSub}.time, betas{iSub}.freq(selFreq), ...
        squeeze(betas{iSub}.powspctrm(iROI, iChan, selFreq, :)));
    set(gca, 'xtick', [], 'ytick', [], 'yscale', 'log', 'linewidth', 3);
%     set(gca, 'visible', 'off')
    xlabel('Time');
    ylabel('Frequency');
    title(sprintf('%s, PPN %d', selChan, iSub), 'fontweight', 'normal');
    
    % Save:
    figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS17F_sub%03d', iSub));
    saveas(gcf, [figName '.png']);
    pause(1);
    close gcf

    % Save source data file:
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS17F_sub%d.csv', iSub));
    csvwrite(fullFileName, squeeze(betas{iSub}.powspctrm(iROI, iChan, selFreq, :)));

end

% ----------------------------------------------------------------------- %
%% Panel F, bottom row: Plot group-level TF data:

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)); nContour = 10; % default 10

job.invalidSubs = [11 12 15 23 25 26 30]; % for fMRI regressors --> TAfT outliers

rng(20190822) % set random number generator for constant p-values

iROI            = 2; % iROI to be tested/ plotted

[sortBetas,~]   = taft_postprocess_TF_selectData(job, betas, iROI);
sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:

selChans        = {'Fz', 'FCz', 'Cz'}; % Jenn's a-priori selection
[tg, corrp]     = taft_postprocess_TF_TFplot(job, dirs, sortBetas, iROI, selChans, 2, 1000, 3, false); % job, dirs, sortBetas, iROI, selChans,thresh,nP,zlim, isSave
close gcf

% Settings:
timeIdx         = size(tg, 4); % maximum time
zlim            = 3;
lineWidth       = 3;
fontSize        = 32; % 44; % uni: 32; smaller: 24 (3/4)
pCrit           = 0.05;

% ----------------------------------------------------------------------- %
% Start figure:
figure('Position', [100 100 1000 800], 'color', 'white'); hold on % narrower for Figure 5 TAfT

% Make contour plot:
contourf(sortBetas{1}.time(1:timeIdx), sortBetas{1}.freq, squeeze(real(tg)), ...
    nContour, 'linestyle', 'none');

% ----------------------------------------------------------------------- %
% Add vertical lines:

plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', 3);
plot([0.7 0.7], get(gca, 'ylim'), ':k', 'LineWidth', 3);

% ----------------------------------------------------------------------- %
% Settings:

set(gca, 'xlim', [-0.25 0.8], 'ylim', [1.5 33], 'clim', [-1*zlim  1*zlim], ...
    'xtick', [-1 -0.5 -0.25 0 0.25 0.5 0.7 1], ...
    'xtickLabel', {'-1', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
    'ytick', [2 4 8 16 32], 'yscale', 'log', ... 
    'fontsize', fontSize, 'Linewidth', lineWidth);

% Labels:
xlabel('Time (ms)', 'FontSize', fontSize, 'FontName', 'Arial', 'fontweight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', fontSize, 'FontName', 'Arial', 'fontweight', 'bold');

% Color bar:
colorbar('Ticks', (-1*zlim):(zlim/2):zlim, 'Fontsize', fontSize);

% ----------------------------------------------------------------------- %
% Add highlights based on p-values:

if sum(corrp(:) < pCrit) > 0 % if anything significant
   isSig = squeeze(double(corrp < pCrit)); % map significant clusters
   % see https://nl.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a#answer_211204

   pause(1) % pause to allow for plotting both plots into one:           
   hold on % on top of old plot

   % Extra contour:
   [~, hContour]  = contourf(sortBetas{1}.time(1:timeIdx), sortBetas{1}.freq, squeeze(isSig), 1);
   hContour.LineWidth = 5;
   drawnow;  % this is important, to ensure that FacePrims is ready in the next line!

   hFills = hContour.FacePrims;  % array of TriangleStrip objects
   [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
   for idx = 1:numel(hFills)
       hFills(idx).ColorData(4) = 1;   % default=255
       pause(1)
   end

   hold off

end
    
% ----------------------------------------------------------------------- %
% Save:
figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS17F_group');
saveas(gcf, [figName '.png']);
pause(1);
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS17F_group.csv');
csvwrite(fullFileName, squeeze(real(tg)));

% END OF FILE.