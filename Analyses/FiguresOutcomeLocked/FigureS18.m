% FigureS17.m

% Plots for Supplementary Figure S18.
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
dirs.EEG            = fullfile(dirs.root, 'Log', 'EEG', 'OutcomeLockedResults');
dirs.TFgroup        = fullfile(dirs.EEG, 'TF_grouplevel');

% Add paths to helper files:
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel'));
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT'));

% Add Fieldtrip:
addpath /home/common/matlab/fieldtrip;
ft_defaults;

% ----------------------------------------------------------------------- %
%% Figure S18A top row: TF plot:

clear job data

% ----------------------------------------------------------------------- %
%% Prepare contrast:

job.dirs                = dirs;
job.nSub                = 36; % necessary for validSubs
job.sub2exclude         = [11 12 23 30]; % exclude 11 12 23 30 because very noisy during trial rejection
job.nFreqs              = 33; % 33
job.baselineSettings    = 'trend';
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'rel'; % 'rel' or 'abs' or 'all'

job.outERPcor           = false;

% ----------------------------------------------------- %
% Channels of interest:
job.chanArea            = 'midfrontal'; % 'midfrontal'
% Band settings:
job.band                = 'theta'; % 'theta' or 'alpha' or 'beta' or 'broad' % lowalpha middlebeta
% Contrast of interest:
job.contrastType        = 'Preferred'; 

if ~exist('data', 'var')
    job         = TF_update_job(job);
    [job, data] = TF_load_data(job);
    [job, data] = TF_prepare_generic_data(job, data);
end

job             = TF_update_job(job);
[~, data]       = TF_prepare_contrast_data(job, data);

% ----------------------------------------------------------------------- %
% Plot:

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)); nContour = 10; % default 10

% Plotting settings:
zlim        = 0.5;
fontSize    = 32;
yScale      = 'log';

% ----------------------------------------------------------------------- %
%% Create TF plot:

figure('Position', [100 100 1200 800], 'color', 'white'); hold on

% Contour plot: 
contourf(data.mu.time, data.mu.freq, data.TF2plot, ...
    nContour, 'linestyle', 'none'); hold on

% Vertical lines:
plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', 3);
plot([0.7 0.7], get(gca, 'ylim'), ':k', 'LineWidth', 3);

% Settings:
set(gca, 'xlim', [-0.25 0.8], 'ylim', [1 33], 'clim', [-1*zlim 1*zlim], ...
    'xtick', [-1 -0.50 -0.25 0 0.25 0.5 0.7 1], ...
    'xtickLabel', {'-1, 000', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
    'ytick', [2 4 8 16 32], 'yscale', yScale, ...
    'fontsize', fontSize, 'Linewidth', 3); % -.25 1.3

% Labels:
ylabel('Frequency (Hz)', 'fontsize', fontSize, 'fontweight', 'bold');
xlabel('Time (ms)', 'fontsize', fontSize, 'fontweight', 'bold');

% Color bar:
colorbar('Ticks', (-1*zlim):(zlim/2):zlim, 'Fontsize', fontSize);

% Save:
figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS18A_top');
saveas(gcf, [figName '.png']);

% Close:
pause(1);
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS18A_top.csv'));
csvwrite(fullFileName, data.TF2plot);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Load TAfT:

EEGdomain       = 'TF';
ROIs2use        = {'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1LeftMotorConj', ...
    'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1LeftITGConj', 'GLM1V1Conj'}; % --> GOES INTO PAPER
behav2use       = {'Updatestd', 'Updatedif'};
selTrials       = 'all';

[job, dirs, betas] = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);


% ----------------------------------------------------------------------- %
%% Load EEG data:

fprintf('Start loading EEG data ... \n');
inputEEG        = load(job.EEGfile);
fprintf('Finished loading EEG data! :-)\n');

% % Extract frequency data:
EEGData         = inputEEG.freq;

% ----------------------------------------------------------------------- %
%% Load EEG mask for panel A bottom:

maskDir = fullfile(dirs.root, 'Log/EEG/OutcomeLockedResults/TF_Mask');

% Theta:
maskName = 'Mask_Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta'; colName = 'blue'; saveName = 'deltatheta';

% Beta:
% maskName = 'Mask_Preferred_CzFCzFz_beta.mat'; colName = 'red'; saveName = 'beta';

fprintf('Load %s \n', maskName);
load(fullfile(maskDir, maskName));

% ----------------------------------------------------------------------- %
%% Panel A: Plot TF mask in 3D space (channels/ frequencies/ time bins):

% Copy to shorter object name:
V       = freqMask;

% Delete channels without any data:
chanIdx = find(squeeze(mean(squeeze(mean(V ~= 0, 3)), 2)) > 0);
V       = squeeze(V(chanIdx, :, :));

% Permute dimensions:
[x,y,z] = meshgrid(1:length(inputEEG.freq.time), 1:length(chanIdx), inputEEG.freq.freq); V = permute(V , [1 3 2]);

% Start plot:
close all

p       = patch(isosurface(x, y, z, V));
isonormals(x, y, z, V, p);

p.FaceColor = colName;
p.EdgeColor = 'none';
set(gca, 'linewidth', 3);
set(gca, 'fontsize', 16);
daspect([1 1 1])
view(3)
camlight; lighting phong

% Labels:
xlabel('Time samples');
ylabel('Channel');
zlabel('Frequency');

% Save:
figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS18A_bottom_%s', saveName));
saveas(gcf, [figName '.png']);
pause(3);
close gcf

% Save source data file:
for iChan = 1:size(V, 1)
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS18A_bottom_%s_chan%s.csv', ...
        saveName, inputEEG.freq.label{chanIdx(iChan)}));
    csvwrite(fullFileName, squeeze(V(iChan, :, :)));
end

% ----------------------------------------------------------------------- %
%% Panel B: Plot EEG per trials:

% Add color bar:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)) ; cmap = 'redblue'; nContour = 10; % default 10

% Settings:
nTrial      = 2; % select trials
selFreq     = 1:15; % select frequencies
selChans    = {'FCz', 'Cz'}; % select channels

close all
for chanName = selChans
    
    iChan = find(ismember(EEGData.label, chanName));

    for iTrial = 1:nTrial

        figure('Position', [100 100 200 200]); hold on
        contourf(EEGData.time, EEGData.freq(selFreq), squeeze(EEGData.powspctrm(iTrial, iChan, selFreq, :)));
        set(gca, 'xtick', [], 'ytick', [], 'linewidth', 3);
        xlabel('Time');
        ylabel('Frequency');
        title(sprintf('%s, Trial %d', chanName{:}, iTrial), 'fontweight', 'normal');

        % Save:
        figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS18B_chan%s_trial%03d', chanName{:}, iTrial));
        saveas(gcf, [figName '.png']);
        pause(1);
        close gcf

        % Save source data file:
        fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS18B_chan%s_trial%d.csv', chanName{:}, iTrial));
        csvwrite(fullFileName, squeeze(EEGData.powspctrm(iTrial, iChan, selFreq, :)));

    end
    
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Load EEG predictor:

% Theta:
% dirs.regressors  = fullfile(dirs.root, 'Log/fMRI/sub-001/GLM3A/timings_regressors/');
% EEGname = 'tPreferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta.txt'; saveName = 'deltatheta';
% colName = 'blue';

% Beta:
dirs.regressors  = fullfile(dirs.root, 'Log/fMRI/sub-001/GLM3B/timings_regressors/');
EEGname = 'tPreferred_CzFCzFz_beta.txt'; saveName = 'beta';
colName = 'red';

% Load regressor:
tmp = load(fullfile(dirs.regressors, EEGname));

% Extract 3rd column:
EEGpred = tmp(:, 3);

% ----------------------------------------------------------------------- %
%% Panel C: Plot EEG trial-by-trial time course:

xMax        = 10;
fontSize    = 18;
LWD         = 4;
markerSize  = 12;

close all

figure('Position', [100 100 100 600], 'Color', 'white'); hold on

xVec        = (1:length(EEGpred)) - 0.5;
yVec        = EEGpred;

% Plot points:
plot(yVec(1:xMax), xVec(1:xMax), ...
    '-o', 'color', colName, 'linewidth', LWD, 'markerSize', 10);

set(gca, 'visible', 'off')

% Save figure:
figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS18C_%s', ...
    saveName));
saveas(gcf, [figName '.png']);
pause(1);
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS18C_%s.csv', saveName));
csvwrite(fullFileName, yVec(1:xMax));

% END OF FILE.