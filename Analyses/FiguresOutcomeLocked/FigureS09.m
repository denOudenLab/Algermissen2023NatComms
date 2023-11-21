% FigureS09.m

% Plots for Supplementary Figure S09.
% Will call taft_save_PEs() to load PEs for selected subject (sub 07).
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

% Scripts:
dirs.scripts    = fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT');

% Plots:
dirs.plot       = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');

% ----------------------------------------------------------------------- %
%% Load PEs for selected subject:

% Select subject; genereate PEs:
iSub    = 7;
out     = taft_save_PEs(iSub, false);
nCond   = 4;
trlMax  = 212; % maximal number trials for any condition (checked for sub 07)

% ----------------------------------------------------------------------- %
%% Figure S09A: Plot PE_STD and PE_BIAS per action x valence over trial repetitions:

FTS = 16; LWD = 2; MKS = 24;
colMat    = [87 196 173; 201 61 33; 87 196 173; 201 61 33] ./ 255; % green, red, green red

% Colour condition:
out.cond    = 2 * (1 - out.isgo) + out.valence;
condNames   = {'Go responses to Win cues', 'Go responses to Avoid cues', ...
    'NoGo responses to Win cues', 'NoGo responses to Avoid cues'};

% Start figure:
close all
p   = cell(1, 2);

figure('Position', [100 100 1000 800], 'Color', 'white'); hold on

% Y-axis: STD & BIAS:
yVec1   = out.PEstd; yName1 = 'PE_{STD}'; yLim = 3.5; yLabel = 'Prediction error term';
yVec2   = out.PEbias; yName2 = 'PE_{BIAS}';

% Initialize for saving:
saveMat1    = nan(trlMax, nCond);
saveMat2    = nan(trlMax, nCond);

% Separately per stimulus:
for iCond = 1:nCond
    subplot(2, 2, iCond); hold on
    
    iRow    = ceil(iCond/2);
    iCol    = mod(iCond - 1, 2) + 1;
    selIdx  = out.cond == iCond;
    xVec    = 1:sum(selIdx);
    p{1}    = plot(xVec, yVec1(selIdx), 'k.', 'LineWidth', 2, 'MarkerSize', MKS);
    p{2}    = plot(xVec, yVec2(selIdx), '.', 'Color', colMat(iCond, :), 'LineWidth', 2, 'MarkerSize', MKS * 3/4); % MKS * 2/4
    
    saveMat1(1:sum(selIdx), iCond)  = yVec1(selIdx);
    saveMat2(1:sum(selIdx), iCond)  = yVec2(selIdx);
    
    % Line at zero:
    plot([1 sum(selIdx)], zeros(1, 2), 'k--', 'LineWidth', 2);

    % Plot settings:
    xlim([min(xVec) max(xVec)]); ylim([-1*yLim 1*yLim]);
    if iRow == 2; xlabel('Trial number'); end
    if iCol == 1; ylabel(yLabel); end
    set(gca, 'FontSize', FTS, 'LineWidth', LWD);

    title(sprintf('%s', condNames{iCond}), 'FontWeight', 'normal');
    
end

% Save:
figName = 'FigS09A.png';
saveas(gcf, fullfile(dirs.plot, figName));
fprintf('Figure saved as %s\n', figName);
pause(2); close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS09A_PE_STD.csv');
csvwrite(fullFileName, saveMat1);
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS09A_PE_BIAS.csv');
csvwrite(fullFileName, saveMat2);

% ----------------------------------------------------------------------- %
%% Figure S09B: Plot PE_STD and PE_DIF per action x valence over trial repetitions:

FTS = 16; LWD = 2; MKS = 24;
colMat    = [87 196 173; 201 61 33; 87 196 173; 201 61 33] ./ 255; % green, red, green red

% Colour condition:
out.cond    = 2 * (1 - out.isgo) + out.valence;
condNames   = {'Go responses to Win cues', 'Go responses to Avoid cues', ...
    'NoGo responses to Win cues', 'NoGo responses to Avoid cues'};

% Start figure:
close all
p   = cell(1, 2);

figure('Position', [100 100 1000 800], 'Color', 'white'); hold on

% Y-axis: STD & DIF:
yVec1   = out.PEstd; yName1 = 'PE_{STD}'; col1 = [0.64 0.08 0.18]; yLim = 3.5; yLabel = 'Prediction error term';
yVec2   = out.PEdif; yName2 = 'PE_{DIF}'; col2 = [0.49 0.99 0];

% Initialize for saving:
saveMat1    = nan(trlMax, nCond);
saveMat2    = nan(trlMax, nCond);

% Separately per stimulus:
for iCond = 1:nCond
    subplot(2, 2, iCond); hold on
    
    iRow    = ceil(iCond/2);
    iCol    = mod(iCond - 1, 2) + 1;
    selIdx  = out.cond == iCond;
    xVec    = 1:sum(selIdx);
    p{1}    = plot(xVec, yVec1(selIdx), 'k.', 'LineWidth', 2, 'MarkerSize', MKS);
    p{2}    = plot(xVec, yVec2(selIdx), '.', 'Color', colMat(iCond, :), 'LineWidth', 2, 'MarkerSize', MKS * 3/4); % MKS * 2/4
    
    saveMat1(1:sum(selIdx), iCond)  = yVec1(selIdx);
    saveMat2(1:sum(selIdx), iCond)  = yVec2(selIdx);
    
    % Line at zero:
    plot([1 sum(selIdx)], zeros(1, 2), 'k--', 'LineWidth', 2);

    % Plot settings:
    xlim([min(xVec) max(xVec)]); ylim([-1*yLim 1*yLim]);
    if iRow == 2; xlabel('Trial number'); end
    if iCol == 1; ylabel(yLabel); end
    set(gca, 'FontSize', FTS, 'LineWidth', LWD);

    title(sprintf('%s', condNames{iCond}), 'FontWeight', 'normal');
    
end

% Save:
figName = 'FigS09B.png';
saveas(gcf, fullfile(dirs.plot, figName));
fprintf('Figure saved as %s\n', figName);
pause(2); close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS09B_PE_STD.csv');
csvwrite(fullFileName, saveMat1);
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS09B_PE_DIF.csv');
csvwrite(fullFileName, saveMat2);

% END OF FILE.