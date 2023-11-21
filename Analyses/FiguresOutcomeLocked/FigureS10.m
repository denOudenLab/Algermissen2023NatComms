% FigureS10.m

% Plots for Supplementary Figure S10.
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
nTrial  = length(out.PEdif);
nRep    = max(out.stimRep);
nStim   = max(out.stim);

% ----------------------------------------------------------------------- %
%% Figure S10A: Time course PE_STD and PE_BIAS:

FTS = 24; LWD = 2; MKS = 32;

% X-axis:
xVec    = 1:nTrial;

% Y-axis: STD & BIAS:
yVec1   = out.PEstd; yName1 = 'PE_{STD}'; col1 = [0.64 0.08 0.18]; yLim = 3.5; yLabel = 'Prediction error term';
yVec2   = out.PEbias; yName2 = 'PE_{BIAS}'; col2 = [0.30 0.75 0.93];

% Start figure:
close all
p   = cell(1, 2);

figure('Position', [100 100 1800 800], 'Color', 'white'); hold on

% Across stimuli:
p{1} = plot(xVec, yVec1, '.', 'Color', col1, 'LineWidth', 2, 'MarkerSize', MKS);
p{2} = plot(xVec, yVec2, '.', 'Color', col2, 'LineWidth', 1, 'MarkerSize', MKS/4*3);

% Vertical lines:
plot([320.5 320.5], [-1*yLim 1*yLim], 'k--', 'LineWidth', LWD);

% Plot settings:
xlim([min(xVec) max(xVec)]); ylim([-1*yLim 1*yLim]);
xlabel('Trial number');
ylabel(yLabel);
set(gca, 'FontSize', FTS, 'LineWidth', LWD);

% Legend:
legend([p{:}], {yName1, yName2}, 'Location', 'northeast', 'Orientation', 'horizontal'); % legend boxoff;

% Save:
figName = 'FigS10A.png';
saveas(gcf, fullfile(dirs.plot, figName));
fprintf('Figure saved as %s\n', figName);
pause(2); 
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS10A.csv');
csvwrite(fullFileName, [yVec1 yVec2]);

% ----------------------------------------------------------------------- %
%% Figure S10B: Correlation PE_STD and PE_BIAS:

FTS = 24; LWD = 2; MKS = 24;

% STD & BIAS:
yVec1   = out.PEstd; yName1 = 'PE_{STD}'; yLim = 3;
yVec2   = out.PEbias; yName2 = 'PE_{BIAS}';

% Plot:
close all
figure('Position', [100 100 800 800], 'Color', 'white'); hold on

% Plot dots:
plot(yVec1, yVec2, 'k.',  'LineWidth', 1, 'MarkerSize', MKS);

% Plotting settings:
axis equal;
xlim([-1*yLim 1*yLim]); 
ylim([-1*yLim 1*yLim]);
xlabel(yName1); ylabel(yName2);
set(gca, 'FontSize', FTS, 'LineWidth', LWD);

% Title:
fprintf('%s and %s: r = %.02f\n', yName1, yName2, corr(yVec1, yVec2))

% Save:
figName = 'FigS10B.png';
saveas(gcf, fullfile(dirs.plot, figName));
pause(2); 
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS10B.csv');
csvwrite(fullFileName, [yVec1 yVec2]);

% ----------------------------------------------------------------------- %
%% Figure S10C: Time course PE_STD and PE_DIF:

FTS = 24; LWD = 2; MKS = 32;

% X-axis:
xVec    = 1:nTrial;

% Y-axis: STD & DIF:
yVec1   = out.PEstd; yName1 = 'PE_{STD}'; col1 = [0.64 0.08 0.18]; yLim = 3.5; yLabel = 'Prediction error term';
yVec2   = out.PEdif; yName2 = 'PE_{DIF}'; col2 = [0.49 0.99 0];

% Start figure:
close all
p   = cell(1, 2);

figure('Position', [100 100 1800 800], 'Color', 'white'); hold on

% Across stimuli:
p{1} = plot(xVec, yVec1, '.', 'Color', col1, 'LineWidth', 2, 'MarkerSize', MKS);
p{2} = plot(xVec, yVec2, '.', 'Color', col2, 'LineWidth', 1, 'MarkerSize', MKS/4*3);

% Vertical lines:
plot([320.5 320.5], [-1*yLim 1*yLim], 'k--', 'LineWidth', LWD);

% Plot settings:
xlim([min(xVec) max(xVec)]); ylim([-1*yLim 1*yLim]);
xlabel('Trial number');
ylabel(yLabel);
set(gca, 'FontSize', FTS, 'LineWidth', LWD);

% Legend:
legend([p{:}], {yName1, yName2}, 'Location', 'northeast', 'Orientation', 'horizontal'); % legend boxoff;

% Save:
figName = 'FigS10C.png';
saveas(gcf, fullfile(dirs.plot, figName));
fprintf('Figure saved as %s\n', figName);
pause(2); 
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS10C.csv');
csvwrite(fullFileName, [yVec1 yVec2]);

% ----------------------------------------------------------------------- %
%% Figure S10D: Correlation PE_STD and PE_DIF:

FTS = 24; LWD = 2; MKS = 24;

% STD & DIF:
yVec1   = out.PEstd; yName1 = 'PE_{STD}'; yLim = 3;
yVec2   = out.PEdif; yName2 = 'PE_{DIF}';

% Plot:
close all
figure('Position', [100 100 800 800], 'Color', 'white'); hold on

% Plot dots:
plot(yVec1, yVec2, 'k.',  'LineWidth', 1, 'MarkerSize', MKS);

% Plotting settings:
axis equal;
xlim([-1*yLim 1*yLim]); 
ylim([-1*yLim 1*yLim]);
xlabel(yName1); ylabel(yName2);
set(gca, 'FontSize', FTS, 'LineWidth', LWD);

% Title:
fprintf('%s and %s: r = %.02f\n', yName1, yName2, corr(yVec1, yVec2))

% Save:
figName = 'FigS10D.png';
saveas(gcf, fullfile(dirs.plot, figName));
pause(2); 
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS10D.csv');
csvwrite(fullFileName, [yVec1 yVec2]);

% END OF FILE.