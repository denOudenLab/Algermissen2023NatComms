function EEGfMRIPav_plot_behavior()

% EEGfMRIPav_plot_behavior()
% 
% Reads data by calling EEGfMRIPav_aggr_data(), then creates
% 1. Line plot of p(Go) per condition over time.
% 2. Bar plot of p(Go) per condition with individual subject means as scattered data points.
% 3. Bar plot of RT per condition with individual subject means as scattered data points.
%
% Mind setting root directory and target directory for saving.
% 
% INPUT:
% none, calls EEGfMRIPav_aggr_data().
%
% Output:
% Creates and saves plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% ----------------------------------------------------------------------- %
%% Set directories:

fprintf('Set root directory\n')
rootDir     = '/project/3017042.02';
targetDir   = fullfile(rootDir, 'Log/CueLockedPaperPlots/');
if ~exist(targetDir, 'dir'); mkdir(targetDir); end

% ----------------------------------------------------------------------- %
%% Add additional paths:

fprintf('Set path to Matlab plots folder\n')
addpath(fullfile(rootDir, 'Analyses/Behavior_Scripts/Matlab_Plots'))

% ----------------------------------------------------------------------- %
%% Retrieve recoded and pre-processed behavioral data:

fprintf('Retrieve recoded and pre-processed behavioral data\n')
out = EEGfMRIPav_aggr_data(); % execute to get out

% ----------------------------------------------------------------------- %
%% Extract relevant fields into separate objects:

fprintf('Extract relevant fields into separate objects\n')
pGoCond         = out.pGoCond;
pCorrectCond    = out.pCorrectCond;
RTCond          = out.RTCond;
RTAcc           = out.RTAcc;
nSub            = size(pGoCond, 1);
nCond           = size(pGoCond, 2);
nCondRT         = size(RTAcc, 2);
nRep            = size(pGoCond, 3);

% ----------------------------------------------------------------------- %
%% Subject selection:

fprintf('Select valid subjects \n')

invalidSubs     = []; % outliers in TAfT and fMRI
% invalidSubs     = [1 11 15 19 21 25 26]; % outliers in TAfT and fMRI
fprintf('Exclude subjects %s\n', num2str(invalidSubs));
validSubs       = setdiff(1:nSub, invalidSubs);
nSubValid       = length(validSubs);

% ----------------------------------------------------------------------- %
%% Figure 1: Learning curve (line plot) per condition:

% Average over conditions per subject:
fprintf('Create mean per subject across conditions\n');
subMean     = squeeze(nanmean(pGoCond(validSubs, :, :), 2)); % average across conditions

% Grand mean over subjects:
fprintf('Create grand mean across subjects\n');
grandMean   = squeeze(nanmean(subMean, 1)); % average across subjects

% Average per condition:
fprintf('Create mean/ SE per condition\n');
condMean    = nan(nCond, nRep);
condSE      = nan(nCond, nRep);
for iCond = 1:nCond
    condMean(iCond, :) = squeeze(nanmean(pGoCond(validSubs, iCond, :), 1));
    condSE(iCond, :)   = nCond / (nCond - 1) * nanstd(squeeze(pGoCond(validSubs, iCond, :)) - ...
        subMean + repmat(grandMean, nSubValid, 1)) ./ sqrt(nSubValid);
end

% Plot settings:
fprintf('Set global plotting settings\n');
% colMat      = [0 .6 .2; .8 0 0;0 .6 .2; .8 0 0]; % old colors
colMat      = [0 113 116; 201 61 33; 0 113 116; 201 61 33] ./ 255; % new colors: dark green, light green, orange, red
FTS    = 38;
xWidth      = 1200;

% ----------------------------------- %
% Start figure:
p           = cell(nCond, 1);
fprintf('Create line plot per condition over time\n');

figure('Position', [100 100 xWidth 800]); 
for iCond = 1:nCond
    p{iCond} = boundedline(1:nRep, condMean(iCond, :), condSE(iCond, :), ...
        'cmap', colMat(iCond, :), 'alpha');
    set(p{iCond}, 'Linewidth', 4)
    if iCond > nCond/2
        set(p{iCond}, 'Linestyle', '--');
    end
end

% Add plot features:
set(gca, 'xlim', [0 nRep], 'ylim', [0 1], 'Linewidth', 3);
set(gca, 'xtick', 0:10:nRep, 'ytick', 0:.2:1, 'FontSize', FTS, 'Linewidth', 3);
ylabel('p(Go)', 'FontSize', FTS, 'FontName', 'Arial');
xlabel('Trial', 'FontSize', FTS, 'FontName', 'Arial');
box off

% Save:
fprintf('Save plot under %s\n', targetDir);
if isempty(invalidSubs)
    saveas(gcf, fullfile(targetDir, 'Behavior_CueConditions_%d.png',xWidth));
else
    saveas(gcf, fullfile(targetDir, 'Behavior_CueConditions_withoutInvalid_%d.png',xWidth));
end
close gcf

% ----------------------------------------------------------------------- %
%% Figure 2. p(Go) bar plot:

fprintf('Create bar plot per condition\n');

% Average across trials (grand means):
fprintf('Average across trials\n');
pGoCondMean         = squeeze(nanmean(pGoCond(validSubs, :, :), 3)); % Go responses
pGoCorrectCondMean  = squeeze(nanmean(pCorrectCond(validSubs, :, :), 3)); % correct responses
subMean             = squeeze(nanmean(pGoCondMean, 2)); % average across conditions
grandMean           = squeeze(nanmean(subMean, 1)); % average across subjects

% Mean per condition:
fprintf('Create mean/ SE per condition \n');
condMean            = nan(nCond, 1);
condMeanCorrect     = nan(nCond, 1);
condSE              = nan(nCond, 1);
for iCond = 1:nCond
    condMean(iCond)         = squeeze(nanmean(pGoCondMean(:, iCond), 1));
    condSE(iCond)           = nCond / (nCond - 1) * nanstd(squeeze(pGoCondMean(:, iCond)) - ...
        subMean + repmat(grandMean, nSubValid, 1)) ./ sqrt(nSubValid);
    condMeanCorrect(iCond)  = squeeze(nanmean(pGoCorrectCondMean(:, iCond), 1));
end

% General plot settings:
% colMat    = [0 .6 .2; .8 0 0; 0 .6 .2; .8 0 0]; % old colors
colMat    = [0 113 116; 201 61 33; 0 113 116; 201 61 33] ./ 255; % new colors: dark green, light green, orange, red
colMatCor = [0.43 0.75 0.39; 0.95 0.44 0.17]; % 110 192 101; 243 112 43
xLoc      = [1 2 3.5 4.5];
LWD = 4; capSize = 12; FTS = 36;

% --------------------------------------------- %
% Start figure:
close all
figure('Position', [100 100 800 800]); hold on

% a) Plot bars with errorbars:
for iCond = 1:nCond
    bar(xLoc(iCond), condMean(iCond), 0.75, 'FaceColor', colMat(iCond, :)); % bar plot
    errorbar(xLoc(iCond), condMean(iCond), condSE(iCond), ...
        'k', 'linestyle', 'none', 'linewidth', LWD, 'Capsize', capSize); % error bars
end

% b) Plot bars with correct Gos:
for iCond = 1:2
    bar(xLoc(iCond), condMeanCorrect(iCond), ...
        0.75, 'FaceColor', colMatCor(iCond, :)); % set(p1, 'FaceAlpha',0.5);
end

% c) Points:
for iCond = 1:nCond
    s = scatter(repmat(xLoc(iCond), 1, nSubValid), pGoCondMean(:, iCond)', ...
        [], 'k', 'jitter', 'on', 'jitterAmount', 0.15); hold on % was 0.05
    set(s, 'MarkerEdgeColor', [0.4 0.4 0.4], 'linewidth', 3); % was 1 
end

% Add plot features:
set(gca, 'xlim', [0.5 5.5], 'ylim', [0 1], 'xtick', 0:10:nRep,...
    'xtick', [1.5 4], 'xticklabel',{'Go', 'NoGo'}, 'ytick', 0:.2:1,...
    'FontSize', FTS, 'FontName', 'Arial', 'FontWeight', 'normal', 'Linewidth', 4);
ylabel('p(Go)', 'FontSize', FTS, 'FontName', 'Arial');
xlabel('Required Action', 'FontSize', FTS, 'FontName', 'Arial');
box off

% Save:
fprintf('Save plot under %s\n', targetDir);
if isempty(invalidSubs)
    saveas(gcf, fullfile(targetDir, 'Behavior_pGoBars.png'));
else
    saveas(gcf, fullfile(targetDir, 'Behavior_pGoBars_withoutInvalid.png'));
end
close gcf

% ----------------------------------------------------------------------- %
%% Figure 3. Bar plot reaction times:

% Select subjects:
RTCondMean          = RTAcc(validSubs, :); % select subjects

% Average per subject across conditions:
subMean             = squeeze(nanmean(RTCondMean, 2)); % average across conditions

% Average over subjects (grand average):
grandMean           = squeeze(nanmean(subMean, 1)); % average across subjects

% Compute mean/SE per condition per subject:
condMean            = nan(nCondRT, 1);
condMeanCorrect     = nan(nCondRT, 1);
condSE              = nan(nCondRT, 1);
for iCond = 1:nCondRT
    condMean(iCond)         = squeeze(nanmean(RTCondMean(:, iCond), 1));
    condSE(iCond)           = nCond / (nCond - 1) * nanstd(squeeze(RTCondMean(:, iCond)) - ...
        subMean + repmat(grandMean, nSubValid, 1)) ./ sqrt(nSubValid);
end

% General plot settings:
yLimMax         = round(max(RTCondMean(:)), 2)+.01;
% colMat          = [0.43 .75 .39; .95 .44 .17; 0 .6 .2; .8 0 0; 0 .6 .2; % .8 0 0]; old colors
colMat          = [0 113 116; 201 61 33; 0 113 116; 201 61 33] ./ 255; % new colors: dark green, light green, orange, red
xLoc            = [1 2 3.5 4.5 6 7];
FTS             = 38;

% --------------------------------------------- %
% Start figure:

figure('Position', [100 100 1000 800]); hold on

% Loop over conditions, bar, error bar, scatters:
for iCond = 1:nCondRT % loop over conditions to create bar and scatter plots
    bar(xLoc(iCond), condMean(iCond), 0.75, ...
        'FaceColor', colMat(iCond, :)); % bar plot
    errorbar(xLoc(iCond), condMean(iCond), condSE(iCond), ...
        'k', 'linestyle', 'none', 'linewidth', LWD, 'Capsize', capSize); % error bars
    s = scatter(repmat(xLoc(iCond), 1, nSubValid), RTCondMean(:, iCond)', ...
        [], 'k', 'jitter', 'on', 'jitterAmount', 0.15); hold on % was 0.05
    set(s, 'MarkerEdgeColor', [0.4 0.4 0.4], 'linewidth', 3); % was 1 
end

% Add plot features:
set(gca, 'xlim', [xLoc(1)-0.5 xLoc(end)+0.5], 'ylim', [0 yLimMax], ...
    'xtick', [1.5 4 5.25 6.5], 'xticklabel', {'correct Go', '\newline(other Go)', 'incorrect Go', '\newline(NoGo)'}, ...
    'ytick', 0:0.2:2,...
    'FontSize', FTS, 'FontName', 'Arial', 'FontWeight', 'normal', 'Linewidth', 4);
% xlabel('Performed (required) action', 'FontSize', FTS, 'FontName', 'Arial', 'Color', [0 0 0]);
ylabel('Reaction time', 'FontSize', FTS, 'FontName', 'Arial', 'Color', [0 0 0]);
box off

% --------------------------------------------- %
% Save:
fprintf('Save plot under %s\n', targetDir);
if isempty(invalidSubs)
    saveas(gcf, fullfile(targetDir, 'Behavior_RTBars.png'))
else
    saveas(gcf, fullfile(targetDir, 'Behavior_RTBars_withoutInvalid.png'))
end
close gcf

fprintf('Finished :-)\n');

end % END OF FUNCTION.
