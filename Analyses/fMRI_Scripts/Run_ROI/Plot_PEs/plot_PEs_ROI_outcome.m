function output = plot_PEs_ROI_outcome(ROIname, ROItitle, GLMID, condOrderName, isPoint, isLegend, yLim)

% Create bar plots for outcome-locked parameter estimates per condition
% extracted from ROI.
% Plot in order of extraction (i.e. Valence, then Action).
% THIS SCRIPT MIGHT ONCE HAVE BEEN JUST A TEMPLATE FOR THE OTHER SCRIPTS???
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs

% ROIname = 'GLM2vmPFC'; ROItitle = ROIname; GLMID = '2'; isPoint = true; yLim = [-0.2 0.2];

if nargin < 2
    ROItitle = ROIname;
end

if nargin < 3
    GLMID   = '2';
end

if nargin < 4
    condOrderName = 'actOut';
end

if nargin < 5
    isPoint = false;
end

if nargin < 6
    isLegend = false;
end

% yLim: see below

% ----------------------------------------------------------------------- %
%% Set directories:

dirs = plot_PE_set_dirs(GLMID);

% ----------------------------------------------------------------------- %
%% Load data:

fileName    = sprintf('%s.txt', ROIname);
fprintf('Load data file %s\n', fileName);
data        = load(fullfile(dirs.data, fileName));
data        = data'; % transpose so subjects in rows, conditions in columns
fprintf('Data %s: Found %d subjects and %d variables\n', ...
    ROIname, size(data, 1), size(data, 2));

% ----------------------------------------------------------------------- %
%% Select data:

nSub        = size(data, 1);
validSubs   = 1:nSub; % 15 and 25 already excluded because failed registration
% validSubs   = [1:10 13:21 23 25:27 29:34]; % drop 11 12 15 23 25 26 30, i.e. 11 12 22 24 28
nSubValid   = length(validSubs);

fprintf('Select %d subjects\n', nSubValid);

data = data(validSubs, :);

% ----------------------------------------------------------------------- %
%% Aggregate and create standard errors:

nSub        = size(data, 1);
nCond       = size(data, 2);
condMean    = nan(nCond, 1);
condSE      = nan(nCond, 1);
subMean     = nanmean(data, 2); % average over conditions --> one value per subject
grandMean   = nanmean(subMean); % average over subjects --> grand mean (1 value)

% Correct SE by subtracting overall mean per subject and adding back grand
% average; correct for # conditions (Cousineau-Morey correction):
for iCond = 1:nCond
    condMean(iCond) = nanmean(data(:, iCond));
    condSE(iCond) = nCond / (nCond - 1) * nanstd(data(:, iCond) - ...
        subMean + repmat(grandMean, nSub, 1)) ./ sqrt(nSub);
end

% Y-axis limits:
if ~exist('yLim', 'var')
    yOffset     = 0.05;
    if isPoint 
        yLim = [min(data(:)) - yOffset, max(data(:)) + yOffset];
    else
        yLim = [min(condMean-condSE) - yOffset, max(condMean+condSE) + yOffset];
    end
end
fprintf('Minimum   : %.03f, Maximum   : %.03f\n', min(data(:)), max(data(:)));
fprintf('y-axis min: %.03f, y-axis max: %.03f\n', min(yLim(1)), yLim(end));

% ----------------------------------------------------------------------- %
%% Fixed plotting parameters:

% a) x-axis ticks:
xOffset = 0.075;
xLoc    = [0.5+xOffset 1.5-xOffset 2.5+xOffset 3.5-xOffset ...
    5.0+xOffset 6.0-xOffset 7.0+xOffset 8.0-xOffset]; % always 1-2 a bit closter together
%xLoc     = [0.5 1.5 2.5 3.5 5 6 7 8];

% b) Color matrix:
% colMat  = [0 0.6 0.2;1 0.63 0.48; 0.49 0.99 0;0.8 0 0]; % old
colMat  = [0 113 116; 87 196 173; 240 174 102; 201 61 33] ./ 255;
colMat  = repmat(colMat, nCond / size(colMat, 1), 1);

% c) Other plotting settings:
LWD     = 3; FTS = 35; opacValue = 0.4; capSize = 12; sz = 24; yTick = -1:.05:1; 

% d) Legend names:
legendNames = {'Rewarded', 'Non-rewarded', 'Non-punished', 'Punished'};

% d) Condition order:
% Original order is:
% GoReward	GoNoReward	GoNoPunishment	GoPunishment	NoGoReward	NoGoNoReward	NoGoNoPunishment	NoGoPunishment
if strcmp(condOrderName, 'actOut')
    condOrder   = 1:nCond; % in original order (just different x-axis ticks)
    xLabel      = 'Performed action';
    xTicks      = [mean(xLoc(2:3)) mean(xLoc(6:7))];
    xTickNames  = {'Go', 'NoGo'};
    legendIdx   = [1 2 3 4];

elseif strcmp(condOrderName, 'actValSal')
    % Order should become:
    % GoReward GoNoPunishment GoNoReward GoPunishment NoGoReward NoGoNoPunishment NoGoNoReward NoGoPunishment 
    condOrder   = [1 3 2 4 5 7 6 8];
    xLabel      = 'Performed action';
    xTicks      = [mean(xLoc(2:3)) mean(xLoc(6:7))];
    xTickNames  = {'Go', 'NoGo'};
    legendIdx   = [1 3 2 4];

elseif strcmp(condOrderName, 'salActVal')
    % Order should become:
    % GoReward GoPunishment NoGoReward NoGoPunishment GoNoPunishment GoNoReward NoGoNoPunishment NoGoNoReward 
    condOrder   = [1 4 5 8 3 2 7 6];    
    xLabel      = 'Performed action';
    xTicks      = [mean(xLoc(1:2)) mean(xLoc(3:4)) mean(xLoc(5:6)) mean(xLoc(7:8))];
    xTickNames  = {'Go', 'NoGo', 'Go', 'NoGo'};
    legendIdx   = [1 6 5 2];

elseif strcmp(condOrderName, 'valActSal')
    % Order should become:
    % GoReward GoNoPunishment NoGoReward NoGoNoPunishment GoPunishment GoNoReward NoGoPunishment NoGoNoReward
    condOrder   = [1 3 5 7 4 2 8 6];
    xLabel      = 'Performed action';
    xTickNames  = {'Go', 'NoGo', 'Go', 'NoGo'};
    xTicks      = [mean(xLoc(1:2)) mean(xLoc(3:4)) mean(xLoc(5:6)) mean(xLoc(7:8))];
    legendIdx   = [1 6 2 5];
    
else 
    error('Unknown condOrderName %s', condOrderName);
end

% Re-sort data for output:
output  = data(:, condOrder);

% ----------------------------------------------------------------------- %
%% Start figure:

p = cell(nCond, 1);
figure('Position', [0 0 1200 800], 'Color', 'white'); hold on

% Bar plots:
for iCond = 1:nCond
    condIdx = condOrder(iCond);
    p{iCond} = bar(xLoc(iCond), condMean(condIdx), ...
        0.75, 'FaceColor', colMat(condIdx, :)); % bar plot
end

% Error bars:
for iCond = 1:nCond
    condIdx = condOrder(iCond);
    errorbar(xLoc(iCond), condMean(condIdx), condSE(condIdx), ...
        'k', 'linestyle', 'none', 'linewidth', LWD, 'Capsize', capSize); % error bars
end

% c) Scatter plots:
if isPoint
    for iCond = 1:nCond
        condIdx = condOrder(iCond);
            s = scatter(repmat(xLoc(condIdx), nSub, 1), data(:, condIdx), [], ...
                'k', 'jitter', 'on', 'jitterAmount', 0.15); hold on % was 0.05
            set(s, 'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerFaceAlpha', opacValue, ...
                'MarkerEdgeAlpha', opacValue, 'linewidth', 3); % was 1 
    end
    yTick = -2:.1:2;
end

% Settings:
set(gca, 'FontSize', FTS-10, 'FontName', 'Arial', 'Linewidth', LWD, ...
    'xlim', [0 xLoc(nCond) + 0.6], ...
    'xtick', xTicks, ...
    'xticklabel', xTickNames, ...
    'ytick', yTick);

if exist('yLim', 'var') % adjust y axis if applicable
    set(gca, 'ylim', yLim);
end

% Axis labels:
xlabel(xLabel, ...
    'FontSize', FTS, 'FontName', 'Arial', 'Color', [0 0 0]);
ylabel('Parameter estimates (a.u.)', ...
    'FontSize', FTS, 'FontName', 'Arial', 'Color', [0 0 0]);
title(sprintf('%s', ROItitle), ...
    'FontSize', FTS, 'FontName', 'Arial', 'Color', [0 0 0]);

% Add legend:
if isLegend
    legend([p{legendIdx}], legendNames); % legend box off
end

% ----------------------------------------------------------------------- %
%% Save:

figName     = sprintf('Matlab_Barplot_%s_%s', ROIname, condOrderName);
if isPoint; figName = sprintf('%s_withPoints', figName); end
if nSubValid < 34; figName = [figName '_without7']; end
fprintf('Save plot as %s\n', figName);
saveas(gcf, fullfile(dirs.plot, [figName '.png']));
% pause(3)
% close gcf

% END OF FILE.