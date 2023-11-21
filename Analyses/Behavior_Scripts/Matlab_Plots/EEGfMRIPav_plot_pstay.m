function EEGfMRIPav_plot_pstay(datatype)

% EEGfMRIPav_plot_pstay()
% 
% Reads data by calling EEGfMRIPav_aggr_data(), then creates bar plot of
% p(Stay) as function of condition specified under datatype input.
%
% Mind setting root directory and target directory for saving.
% 
% INPUT:
% datatype          = string, either 'OutVal' (split only by outcome
% valence) or 'SalActOutVal' (split by action and outcome)
%
% Output:
% Creates and saves plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% ----------------------------------------------------------------------- %
%% Complete settings:

if nargin < 1
    datatype = 'SalActOutVal';
end

% ----------------------------------------------------------------------- %
%% Set directories:

fprintf('Set root directory\n')
rootDir     = '/project/3017042.02';
targetDir   = fullfile(rootDir, 'Log/OutcomeLockedPaperPlots/');
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

stayOutVal          = out.stayOutVal;
staySalActOutVal    = out.staySalActOutVal;

% round(staySalActOutVal, 2)
nSub                = size(staySalActOutVal, 1);
nCond               = size(staySalActOutVal, 2);

% ----------------------------------------------------------------------- %
%% Subject selection:

fprintf('Select valid subjects \n')
invalidSubs     = []; % outliers in TAfT and fMRI
fprintf('Exclude subjects %s\n', num2str(invalidSubs));
validSubs       = setdiff(1:nSub, invalidSubs);
nSubValid       = length(validSubs);

% ----------------------------------------------------------------------- %
%% General plot settings:

fprintf('Set general plotting settings\n');
nCond     = size(out.staySalActOutVal, 2);
condNames = {'Go Reward', 'Go No Reward', 'Go No Punishment', 'Go Punishment', ...
    'NoGo Reward', 'NoGo No Reward', 'NoGo No Punishment', 'NoGo Punishment'}; % input order
% colAll    = [0 0.6 0.2; 1 0.63 0.48; 0.49 0.99 0; 0.8 0 0]; % old colors
colAll    = [0 113 116; 201 61 33; 0 113 116; 201 61 33] ./ 255; % new colors: dark green, light green, orange, red
xLoc      = nan(nCond, 1); % position of bars on x-axis

% Plotting locations:
for iCond = 1:nCond
    if iCond <= nCond/2
        xLoc(iCond) = iCond; % first half
    else
        xLoc(iCond) = iCond + 0.5; % second half
    end
end

xMax      = xLoc(nCond)+0.5;
LWD = 4; capSize = 12; FTS = 36;

% ----------------------------------------------------------------------- %
%% Select data:

colMat   = nan(nCond, 3);
if strcmp(datatype, 'OutVal') % outdated   
    data  = stayOutVal;
    xLabel      = 'Outcome';
    xTick       = [xLoc(1) xLoc(2)];
    xTickLabel  = {'Positive', 'Negative'};
    colMat     	= colAll(1:nCond, :); % just as many colors as conditions
elseif strcmp(datatype, 'SalActOutVal')    
    data  	= staySalActOutVal;
    condOrder   = [1 4 5 8 3 2 7 6]; % for splitting up - i would suggest to split up for 'neutral' and 'nonNeutral' outcomes I think, as the former is where you'd expect to see the biased learning. 
    colMat      = colAll([1:nCond/2 1:nCond/2], :); % just as many colors as conditions
    
    % x-labels for action only:
    xTick       = [(xLoc(1)+xLoc(2))/2 (xLoc(3)+xLoc(4))/2 (xLoc(5)+xLoc(6))/2 (xLoc(7)+xLoc(8))/2];
    xTickLabel  = {'Go', 'NoGo', 'Go', 'NoGo'};
    xLabel      = 'Performed action';
    
elseif strcmp(datatype, 'ActCueValOutVal') % Option 2: Action - Cue Valence - Outcome Valence
    data        = pStay;
    condOrder   = [1 2 3 4 5 6 7 8]; % GoRew, GoNoRew, GoNoPun, GoPun, NoGoRew, NoGoNoRew, NoGoNoPun, NoGoPun
    colMat      = colAll([1:nCond/2 1:nCond/2], :); % just as many colors as conditions
    
    % x-labels for action only:
    xTick       = [(xLoc(2)+xLoc(3))/2 (xLoc(6)+xLoc(7))/2];
    xTickLabel  = {'Go', 'NoGo'};
    xLabel      = 'Performed action';
    
elseif strcmp(datatype, 'CueValOutValAct') % Option 3: Cue Valence - Outcome Valence - Action
    data        = pStay;
    condOrder   = [1 5 2 6 3 7 4 8]; % GoRew, NoGoRew, GoNoRew, NoGoNoRew, GoNoPun, NoGoNoPun, GoPun, NoGoNoPun
    colMat      = colAll([1:nCond/2 1:nCond/2], :); % just as many colors as conditions
    
    % x-labels for action only:
    xTick       = [(xLoc(2)+xLoc(3))/2 (xLoc(6)+xLoc(7))/2];
    xTickLabel  = {'Win ', 'Avoid'};
    xLabel      = 'Cue valence';
    
elseif strcmp(datatype, 'SalOutValAct') % Option 4: Salience - Outcome Valence - Action
    data        = pStay;
    condOrder   = [1 5 4 8 2 6 3 7]; % GoRew, NoGoRew, GoPun, NoGoPun, GoNoRew, NoGoNoRew, GoNoPun, NoGoNoPun
    colMat      = colAll([1:nCond/2 1:nCond/2], :); % just as many colors as conditions
    
    % x-labels for salience only:
    xTick       = [(xLoc(2)+xLoc(3))/2 (xLoc(6)+xLoc(7))/2];
    xTickLabel  = {'Salient', 'Neutral'};
    xLabel      = 'Salience';
    
else
    error('Unknown data type')
end

% ----------------------------------------------------------------------- %
%% Correct SEs for between-subjects variability:

fprintf('Compute condition means and SEs correcting for between-subjects variability\n');
nCond               = size(data, 2);
subMean             = squeeze(nanmean(data, 2)); % average across conditions
grandMean           = squeeze(nanmean(subMean, 1)); % average across subjects
condMean            = nan(nCond, 1);
condSE              = nan(nCond, 1);

fprintf('\n')
for iCond = 1:nCond
    condMean(iCond)         = squeeze(nanmean(data(validSubs, iCond), 1)); % average over subjects
    fprintf('Condition %s: M = %.02f\n', condNames{iCond}, condMean(iCond));
    condSE(iCond)           = nCond / (nCond - 1) * nanstd(squeeze(data(validSubs, iCond)) - ...
        subMean + repmat(grandMean, nSubValid, 1)) ./ sqrt(nSubValid);
end

% ----------------------------------------------------------------------- %
%% Start figure:

fprintf('Start figure\n');

p = cell(nCond, 1);
figure('Position', [100 100 800 800]); hold on

% ----------------------------------------------------------------------- %
% a) Plot bars with errorbars:
for iCond = 1:nCond
    condIdx = condOrder(iCond);
    p{iCond} = bar(xLoc(iCond), condMean(condIdx), ...
        0.75, 'FaceColor', colMat(condIdx, :)); % bar plot
    errorbar(xLoc(iCond), condMean(condIdx), condSE(condIdx), ...
        'k', 'linestyle', 'none', 'linewidth', LWD, 'Capsize', capSize); % error bars
end

% ----------------------------------------------------------------------- %
% b) Individual data points:
for iCond = 1:nCond
    condIdx = condOrder(iCond);
    s = scatter(repmat(xLoc(iCond), 1, nSubValid),data(:, condIdx)', ...
        [], 'k', 'jitter', 'on', 'jitterAmount', 0.15); hold on % was 0.05
    set(s, 'MarkerEdgeColor', [0.4 0.4 0.4], 'linewidth', 3); % was 1 
end
plot([0 xMax], ones(2, 1)*1/3, 'k--', 'linewidth', LWD); % line at chance

% ----------------------------------------------------------------------- %
% Add plot features:
set(gca, 'xlim', [.5 xMax], 'ylim', [0 1], ...
    'xtick',xTick, 'xticklabel',xTickLabel, 'ytick', 0:.2:1, ...
    'FontSize', FTS, 'FontName', 'Arial', 'FontWeight', 'normal', 'Linewidth', 4);
xlabel(xLabel, 'FontSize', FTS, 'FontName', 'Arial');
ylabel('p(Stay)', 'FontSize', FTS, 'FontName', 'Arial');
if strcmp(datatype, 'SalActOutVal')
    legend([p{:}],{'Positive', 'Negative'}); box off
end

fprintf('Save plot under %s\n',targetDir);
saveas(gcf, fullfile(targetDir, sprintf('Behavior_Stay_%s.png',datatype)));
pause(3)
close gcf

fprintf('Finished :-)\n');

end % END OF FUNCTION.
