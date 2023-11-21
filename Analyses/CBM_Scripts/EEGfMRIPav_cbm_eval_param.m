% EEGfMRIPav_eval_param.m

% This is an interactive script---execute it step-by-step.
% It fits a serious of computational reinforcement learning models using the CBM toolbox.
% and evaluates them.
%
% EXPLANATION OF SETTINGS:
% dirs.root               = string, root directory of project.
%
% OUTPUTS:
% no outputs, just plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd dirs.root/Analyses/CBM_Scripts/

% clear all; close all; clc

% ----------------------------------------------------------------------- %
%% 00a) Set directories:

% Directories:
dirs.root   = '/project/3017042.02';
dir.results = fullfile(dirs.root, '/Log/Behavior/Modelling_CBM');
dir.lap     = fullfile(dir.results, 'LAP_Results');
dir.hbi     = fullfile(dir.results, 'HBI_Results');

% Add paths:
addpath('/home/action/johalg/cbm-master/codes'); % CBM toolbox
addpath(fullfile(dirs.root, '/Analyses/CBM_Scripts/models')); % models

% ----------------------------------------------------------------------- %
%% 00b) Select model:

nMod    = 9; % number of models defined by code

selMod  = 1; % model to evaluate
fprintf('Selected model is no. %02d\n', selMod);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 01) LOAD AND SAVE:
%% 01a) Alternative A: Load model fitted with LAP:

parType     = 'lap'; % use only for LAP

% Possible output names:
fname_mod   = cell(nMod, 1);
for iMod = 1:nMod
    fname_mod{iMod} = fullfile(dir.lap, sprintf('lap_mod%02d.mat', iMod));
end

% Load model:
fprintf('Load model %d fit with LAP\n', selMod);
fname       = load(fname_mod{selMod});
cbm         = fname.cbm;
subParam    = cbm.output.parameters;
nSub        = size(subParam, 1);
nParam      = size(subParam, 2);

% ----------------------------------------------------------------------- %
%% 01b) Alternative B: Load model fitted with HBI:

parType     = 'hbi'; % use only for hbi

% Load:
fname_hbi   = fullfile(dir.hbi, sprintf('hbi_mod_%02d.mat', selMod)); % all models
fprintf('Load model %d fit with HBI\n', selMod);
load(fname_hbi);

% ----------------------------- %
% Extract parameters:
fprintf('Extract model %d fit with HBI\n', selMod);
groupParam  = cbm.output.group_mean{:};
subParam    = cbm.output.parameters{:};
nSub        = size(subParam, 1);
nParam      = size(subParam, 2);

% ----------------------------------------------------------------------- %
%% 01c) Transform group and subject parameters appropriately:

fprintf('Transform parameters of model %d\n', selMod);

for iParam = 1:nParam
    if iParam == 1 || iParam == 6 && selMod == 9 % inverse temperature (always first parameter) or neural outcome revaluation
        groupParam(iParam)      = exp(groupParam(iParam));
        subParam(:, iParam)     = exp(subParam(:, iParam)); 
        fprintf('Exponentially transform parameter %d\n', iParam); 
    end
    
    if iParam == 2 || iParam == 5 && selMod == 6 % learning rates
        groupParam(iParam)      = 1 ./ (1 + exp(-groupParam(:, iParam)));
        subParam(:, iParam)     = 1 ./ (1 + exp(-subParam(:, iParam))); 
        fprintf('Sigmoid transform parameter %d\n', iParam); 
    end
end

% ----------------------------------------------------------------------- %
%% 01d) Save (either LAP or HBI):

fprintf('Save transformed group and subject level parameters under %s, model %02d\n', ...
    parType, selMod);

fileName = fullfile(eval(sprintf("dir.%s",parType)), sprintf('CBM_%s_M%02d_groupParameters.csv', ...
    parType, selMod));
csvwrite(fileName, groupParam);

fileName = fullfile(eval(sprintf("dir.%s",parType)), sprintf('CBM_%s_M%02d_subjectParameters.csv', ...
    parType, selMod));
csvwrite(fileName, subParam);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 02) GROUP-LEVEL PARAMETERS.
%% 02a) Mean per proup level parameter:

fprintf('\nGroup-level parameters for model %d:\n', selMod);

for iParam = 1:size(groupParam, 2)
    x = groupParam(:, iParam);
    fprintf('Parameter %d: M = %.02f\n', iParam, mean(x));
end

% ----------------------------------------------------------------------- %
%% 02b) Epsilon plus/ minus kappa: group level:

fprintf('\nGroup-level parameters: epsilon +/- kappa\n')
eps     = cbm.output.group_mean{:}(2);
kappa   = cbm.output.group_mean{:}(5);
epsPlusKappa    = 1 ./ (1 + exp(-(eps + kappa)));
epsBase         = 1 ./ (1 + exp(-(eps)));
epsMinusKappa   = 1 ./ (1 + exp(-(eps - kappa)));

fprintf('Epsilon minus kappa: %.03f\n', epsMinusKappa);
fprintf('Epsilon 0:           %.03f\n', epsBase);
fprintf('Epsilon plus kappa:  %.03f\n', epsPlusKappa);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 03) SUBJECT-LEVEL PARAMETERS.
%% 03a) Mean, SD, range per subject-level parameter:

fprintf('\nSubject-level parameters for model %d\n', selMod);

for iParam = 1:size(subParam, 2)
    x = subParam(:, iParam);
    fprintf('Parameter %d: M = %.03f, SD = %.03f, range %.03f - %.03f \n', ...
        iParam, mean(x), std(x), min(x), max(x));
end

% ----------------------------------------------------------------------- %
%% 03b) Epsilon plus/ minus kappa: subject level:

fprintf('\nSubject-level parameters: epsilon +/- kappa\n')
epsPlusKappa    = nan(nSub, 1);
epsBase         = nan(nSub, 1);
epsMinusKappa   = nan(nSub, 1);

for iSub = 1:nSub
    if strcmp(parType, 'hbi')
        eps     = cbm.output.parameters{:}(iSub, 2); % untransformed
        kappa   = cbm.output.parameters{:}(iSub, selMod); % untransformed
    elseif strcmp(parType, 'lap')
        eps     = cbm.output.parameters(iSub, 2); % untransformed
        kappa   = cbm.output.parameters(iSub, selMod); % untransformed
    else
        error('Unknown parType')
    end
    epsMinusKappa(iSub)     = 1 ./ (1 + exp(-(eps - kappa)));
    epsBase(iSub)           = 1 ./ (1 + exp(-(eps)));
    epsPlusKappa(iSub)      = 1 ./ (1 + exp(-(eps + kappa)));
end

fprintf('Epsilon minus kappa:   M = %.03f, SD = %.03f, range %.03f - %.03f \n', ...
    mean(epsMinusKappa), std(epsMinusKappa), min(epsMinusKappa), max(epsMinusKappa));
fprintf('Epsilon base:          M = %.03f, SD = %.03f, range %.03f - %.03f \n', ...
    mean(epsBase), std(epsBase), min(epsBase), max(epsBase));
fprintf('Epsilon plus kappa:    M = %.03f, SD = %.03f, range %.03f - %.03f \n', ... 
    mean(epsPlusKappa), std(epsPlusKappa), min(epsPlusKappa), max(epsPlusKappa));

% ----------------------------------------------------------------------- %
%% 03c) Find extreme values:

fprintf('Find subjects with extreme values for model %d\n', selMod);
for iParam = 1:size(subParam, 2)
    x       = subParam(:, iParam);
    xmin    = min(x);
    xminidx = find(x == xmin);
    xmax    = max(x);
    xmaxidx = find(x == xmax);
    fprintf('Parameter %d: min = %.05f for subjects %s, max = %.02f for subjects %s \n', ...
        iParam, xmin, num2str(xminidx), xmax, num2str(xmaxidx));
end

% ----------------------------------------------------------------------- %
%% 03d) Sign of parameter:

fprintf('Determine percentage negative parameters for model %d\n', selMod);
for iParam = 1:size(subParam, 2)
    x       = subParam(:, iParam);
    xNeg    = sum(x < 0); % how many negative
    xIdx    = find(x < 0)'; % who negative
    fprintf('Parameter %d: negative for %d subjects: %s \n', ...
        iParam, xNeg, num2str(xIdx));
end

% format long
% [[1:nSub]' x] % how negative
% format short

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 04) PLOTS.
%% 04a) Plot parameters with cbm function:

% (looks rather ugly, might use scientific notation for y-axis)
% fname_hbi = fullfile(dir.hbi, sprintf('hbi_mod%02d-%02d.mat', 1,nMod));
model_names = {'Q', 'Q + Go', 'Q + Go + \pi', 'Q + Go + \kappa', 'Q + Go + \pi + \kappa'}; % M1-5
param_names = {'\rho', '\epsilon', 'Go', '\pi', '\kappa'}; % M5
transform   = {'exp', 'sigmoid', 'none', 'none', 'none'}; % M5
cbm_hbi_plot(fname_hbi, model_names, param_names, transform);

% ----------------------------------------------------------------------- %
%% 04b) Density plots of single parameters:

f_hbi   = load(fname_hbi);
cbm     = f_hbi.cbm;
selMod    = 1; % specify model to look at
% Click through densityplots of each parameter:
for iParam = 1:size(cbm.output.parameters{selMod}, 2) % iParam = 2
    x = cbm.output.parameters{selMod}(:, iParam);
    ksdensity(x); hold on
    plot(x, 0.1 + randn(length(x), 1)/50, 'b.');
    title(sprintf('Model %02d: Parameter %d', selMod, iParam));
    w = waitforbuttonpress;
    close gcf
end
% https://stackoverflow.com/questions/5757387/how-to-draw-probability-density-function-in-matlab?rq=1

% ----------------------------------------------------------------------- %
%% 04c) Correlations/ bivariate distributions:

% Select parameters:
iParam1 = 3; % specify first parameter
iParam2 = 4; % specify second parameter

% Retrieve parameter values:
x       = subParam(:, iParam1);
y       = subParam(:, iParam2);
p       = polyfit(x', y', 1);
yhat    = polyval(p, x);

% Plot:
plot(x, y, 'b.'); hold on;
plot(x, yhat, 'r-');
axis([min(x) max(x) min(y) max(y)]);

% Print correlation to console:
corVal  = corr(x, y);
fprintf('Model %02d: Parameters %d and %d, cor = %0.2f\n', ...
    selMod, iParam1, iParam2, corVal);
title(sprintf('Model %02d: Parameters %d and %d, cor = %0.2f', ...
    selMod, iParam1, iParam2, corVal));

% ----------------------------------------------------------------------- %
%% 04d) Bar plots with dots:

% Settings:
paramNames  = {{'$\rho$', '$\epsilon$'},...
    {'$\rho$', '$\epsilon$', '{\it Go}'},...
    {'$\rho$', '$\epsilon$', '{\it Go}', '$\pi$'},...
    {'$\rho$', '$\epsilon$', '{\it Go}', '$\kappa$'},...
    {'$\rho$', '$\epsilon$', '{\it Go}', '$\pi$', '$\kappa$'}};
ScatterMatrix = subParam';
colMat = [0 0 1; 1 0 0; 0 01 0; 0 1 1; 1 0 1;1 0 0; 1 0 0]; % color
posMat = 1:1:(0.5 + nParam);

% Make figure:
addpath(fullfile(dirs.root, '/Analyses/Behavior_Scripts/Matlab_Plots/'); % add function barScatter
figure('units', 'normalized', 'outerposition', [0 0 1 1]); hold on
barScatter(ScatterMatrix, [], [], true, true, colMat, posMat);

% Add plot features:
set(gca, 'xlim', [0 (nParam+1)], 'ylim', [-5 5],...
    'xtick', 1:nParam, 'xticklabel', paramNames{selMod},...
    'FontSize', 32, 'FontName', 'Arial', 'Linewidth', 4, 'TickLabelInterpreter', 'latex');
xlabel('Parameter', 'FontSize',32, 'FontName', 'Arial');
ylabel('Parameter estimates', 'FontSize',32, 'FontName', 'Arial');
box off; hold off

% Print to console:
fprintf('Model %02d:\n', selMod)
for ii = 1:size(ScatterMatrix, 1)
    paramVec = ScatterMatrix(ii, :);
    fprintf('Parameter %02d: M = %.02f, SD = %.02f, range %.02f - %.02f\n', ...
        ii, nanmean(paramVec), nanstd(paramVec), min(paramVec), max(paramVec));
end

% Save:
saveas(gcf, fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/Parameters_%s_mod%02d.png', ...
    parType, selMod)));
pause(2)
close gcf

% END OF FILE.
