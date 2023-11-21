% EEGfMRIPav_cbm_parameter_recovery.m

% This is an interactive script---execute it step-by-step.
% Its performs parameter recovery for a selected model.
% The following steps have to be executed beforehand:
% - Fit model using EEGfMRIPav_cbm_fit to obtain parameters fitted to
% actual data.
% - Simulate new data using EEGfMRIPav_cbm_sim given parameters fitted.
% Mind to set the root directory dirs.root.
%
% EXPLANATION OF SETTINGS:
% dirs.root     = string, root directory of project.
% iMod          = numeric integer, ID of selected model for which to
% perform model recovery.
% simType       = scalar string, type of simulations, should be 'modSim'
% for probabilistic simulations of actual responses ('osap' will only give
% response probabilities).
% parType       = scalar string, type of parameters used for simulations,
% either LaPlace approximation ('LAP') or hierarchical Bayesian inference
% (HBI).
% nIter         = scalar integer, number of simulations (will be averaged
% to get 'recovered' parameter).
% nSub          = scalar integer, number of subjects.
%
% OUTPUTS:
% no outputs, just plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/CBM_Scripts/

% clear all; close all; clc

% ----------------------------------------------------------------------- %
%% 00a) Complete settings:

% Other settings:
iMod        = 5;
simType     = 'modSim'; % modSim osap
parType     = 'hbi'; % lap hbi
nIter       = 100; % number simulations to be recovered
nSub        = 36; % number subjects

fprintf('Perform model recovery for M%02d, simType %s, parType %s \n', ...
    iMod, simType, parType);

% ----------------------------------------------------------------------- %
%% 00b) Set directories:

dirs.root    = '/project/3017042.02'; % root directory, needs to be adjusted to local folder structure.
dirs.scripts = fullfile(dirs.root, 'Analyses/CBM_Scripts');
dirs.results = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');
dirs.lap     = fullfile(dirs.results, 'LAP_Results');
dirs.hbi     = fullfile(dirs.results, 'HBI_Results');
dirs.sim     = fullfile(dirs.results, 'Simulations');

% Create new directories if necessary:
dirs.paramRecovery = fullfile(dirs.results, 'parameterRecovery');
if ~exist(dirs.paramRecovery, 'dir'); mkdir(dirs.paramRecovery); end

% ----------------------------------------------------------------------- %
%% 00c) Add paths:

fprintf('Add paths\n');

addpath('/home/action/johalg/cbm-master/codes'); % add CBM
addpath(fullfile(dirs.root, 'Analyses/CBM_Scripts/models')); % models
addpath(fullfile(dirs.root, 'Analyses/CBM_Scripts')); % for sim_subj
% cd /project/3017042.02/Analyses/CBM_Scripts

addpath '/home/common/matlab/fieldtrip/qsub'; % qsubcellfun in Fieldtrip

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 01a) Load true parameters:

% Must have been generated with EEGfMRIPav_cbm_fit.

fprintf('Load ground true parameters for M%02d, parType %s\n', ...
    iMod, parType);

fprintf('Load parameters based on %s\n',parType)
if strcmp(parType, 'lap')
    lap_name_mod    = fullfile(dirs.lap,sprintf('lap_mod%02d.mat', iMod));
    fname           = load(lap_name_mod);
    cbm             = fname.cbm;
    trueParam       = cbm.output.parameters;
elseif strcmp(parType, 'hbi')
    hbi_name_mod    = fullfile(dirs.hbi,sprintf('hbi_mod_%02d.mat', iMod)); % hbi output object
    fname           = load(hbi_name_mod);
    cbm             = fname.cbm;
    trueParam       = cbm.output.parameters{:};
end

% Initialize number parameters:
nParam = size(trueParam, 2);
% trueParam

% ----------------------------------------------------------------------- %
%% 01b) Load model simulations:

% Must have been generated with EEGfMRIPav_cbm_sim.

fprintf('Load modSim simulations for M%02d, parType %s\n', ...
    iMod, parType);

simFile     = fullfile(dirs.sim, sprintf('%s_mod%02d_%s_iter%04d.mat', ...
    simType, iMod, parType, nIter));

tmp         = load(simFile);
simulations = tmp.sim;

% simulations.action % nSub x nIter x nTrial

% ----------------------------------------------------------------------- %
%% 01c) Prepare fitting simulated data:

% Load data used for model simulations:
subj = sim_subj; % sim_subj() is a function in Scripts

% Priors:
fprintf('Initialize priors\n')
% var: 2 for learning rate, otherwise 3
priors{1} = struct('mean', [2 0], 'variance', [3 2]); % note dimension of 'mean'
priors{2} = struct('mean', [2 0 0], 'variance', [3 2 3]); % note dimension of 'mean'
priors{3} = struct('mean', [2 0 0 0], 'variance', [3 2 3 3]); % note dimension of 'mean'
priors{4} = struct('mean', [2 0 0 0], 'variance', [3 2 3 2]); % note dimension of 'mean'
priors{5} = struct('mean', [2 0 0 0 0], 'variance', [3 2 3 3 2]); % note dimension of 'mean'
priors{6} = struct('mean', [2 0 0 0 0], 'variance', [3 2 3 3 2]); % note dimension of 'mean'
priors{7} = struct('mean', [2 0 0 0 0 0], 'variance', [3 2 3 3 2 3]); % note dimension of 'mean'
priors{8} = struct('mean', [2 0 0 0 0 0 0], 'variance', [3 2 3 3 2 3 3]); % note dimension of 'mean'
priors{9} = struct('mean', [2 0 0 0 0 0 0], 'variance', [3 2 3 3 2 3 3]); % note dimension of 'mean'

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 02a) Start Laplace approximation (LAP) via qsubeval jobs for each iteration:

% For each iteration:
% - add simulated actions to data
% - create file name 

rng(70)

% for iIter = [50 56 69 70 72 99]
for iIter = 1:nIter % iIter = 1;
    
    % ------------------------------------------------------------------- %
    % Prepare data for each subject for given iteration:
    simData = cell(nSub, 1);
    for iSub = 1:nSub % iSub = 1;
        simData{iSub} = subj;
        simData{iSub}.actions = squeeze(simulations.action(iSub, iIter, :));
        simData{iSub}.outcome = squeeze(simulations.outcome(iSub, iIter, :));
    end
    
    % ------------------------------------------------------------------- %
    % Output name:
    fname_mod = fullfile(dirs.paramRecovery, sprintf('lap_mod%02d_iter%04d.mat', iMod, iIter));
    
    % ------------------------------------------------------------------- %
    % Fit with Laplace approximation sequentially:
    cbm_lap(simData, eval(sprintf('@EEGfMRIPav_cbm_mod%02d', iMod)), priors{iMod}, fname_mod);

    % ------------------------------------------------------------------- %
    % Fit with Laplace approximation in parallel with qsubfeval on the
    % Donders HCP cluster:
%     req_mem   = 1024^3 * 5; % (bytes * KB * MB) * GB
%     req_etime = 60 * 60; % seconds * minutes % 30 minutes to be on the safe side
%     qsubfeval(@cbm_lap, simData, eval(sprintf('@EEGfMRIPav_cbm_mod%02d', iMod)), ...
%         priors{iMod}, fname_mod, 'memreq',  req_mem,  'timreq',  req_etime);

end

% ----------------------------------------------------------------------- %
%% 02b) Read in LAP parameters fitted to simulated data:

fprintf('Load fitted parameters for lap_mod%02d_iter%04d.mat \n', iMod, iIter);

parType     = 'lap';
fittedParam = nan(nSub, nParam, nIter);

for iIter = 1:nIter % iIter = 1;
    fname_lap = fullfile(dirs.paramRecovery,sprintf('lap_mod%02d_iter%04d.mat', iMod, iIter));
    tmp = load(fname_lap);
    fittedParam(:, :, iIter) = tmp.cbm.output.parameters;
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 03a) Fit with Hierarchical Bayesian Inference (HBI) via qsubeval jobs for each iteration:

if strcmp(parType, 'hbi')

    rng(70)

    % Model code:
    models = {eval(sprintf('@EEGfMRIPav_cbm_mod%02d', iMod))};

    for iIter = 1:nIter % iIter = 1;
%     for iIter = [11 41 50 79] % iIter = 1;
        % ------------------------------------------------------------------- %
        % Prepare data for each subject for given iteration:
        simData = cell(nSub, 1);
        for iSub = 1:nSub % iSub = 1;
            simData{iSub} = subj;
            simData{iSub}.actions = squeeze(simulations.action(iSub, iIter, :));
            simData{iSub}.outcome = squeeze(simulations.outcome(iSub, iIter, :));
        end
        % ------------------------------------------------------------------- %
        % LAP fitted object:
        fcbm_maps = {fullfile(dirs.paramRecovery,sprintf('lap_mod%02d_iter%04d.mat', iMod, iIter))};
        % ------------------------------------------------------------------- %
        % HBI output file:
        fname_hbi = {fullfile(dirs.paramRecovery,sprintf('hbi_mod%02d_iter%04d.mat', iMod, iIter))};
        % ------------------------------------------------------------------- %
        % Fit with HBI sequentially:
        cbm_hbi(simData, models, fcbm_maps, fname_hbi);
        % ------------------------------------------------------------------- %
        % Fit with HBI in parallel with qsubfeval:
%         req_mem   = 1024^3 * 5; % (bytes * KB * MB) * GB
%         req_etime = 60 * 45; % seconds * minutes % 10-15 minutes
%         qsubfeval(@cbm_hbi, simData, models, fcbm_maps, fname_hbi, 'memreq',  req_mem,  'timreq',  req_etime);

    end
    
else
    
    warning('parType is %s, skip HBI fitting\n', parType);
    
end
% ----------------------------------------------------------------------- %
%% 03b) Read in HBI parameters fitted to simulated data:

fprintf('Load fitted parameters for hbi_mod%02d_iterXXXX.mat \n', iMod);

parType     = 'hbi';
fittedParam = nan(nSub, nParam, nIter);

for iIter = 1:nIter % iIter = 1;
    fname_hbi = fullfile(dirs.paramRecovery,sprintf('hbi_mod%02d_iter%04d.mat', iMod, iIter));
    tmp = load(fname_hbi);
    fittedParam(:, :, iIter) = tmp.cbm.output.parameters{:};
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 04a) Descriptives and averaging:

% Descriptives of parameters:
fittedParamAvgSub   = squeeze(mean(fittedParam, 1)); % average over subjects

% Average over iterations:
fittedParamAvgIter  = squeeze(mean(fittedParam, 3)); % average over iterations

% Descriptives of parameters:
fprintf('Mean parameter averaged over subjects: \n');
for iParam = 1:nParam
    fprintf('Parameter %d: M = %.02f, SD = %.02f, Min = %.02f, Max = %.02f\n', ....
        iParam, mean(fittedParamAvgSub(iParam, :)), std(fittedParamAvgSub(iParam, :)), ...
        min(fittedParamAvgSub(iParam, :)), max(fittedParamAvgSub(iParam, :)));
end

% ----------------------------------------------------------------------- %
%% 04b) Transform paramaters:

rhoTrue         = exp(trueParam(:, 1));
rhoFitted       = exp(fittedParamAvgIter(:, 1));

epsilonTrue     = 1 ./ (1 + exp(-1*trueParam(:, 2))); 
epsilonFitted   = 1 ./ (1 + exp(-1*fittedParamAvgIter(:, 2))); 

% Store in objects:
trueParam(:, nParam + 1) = rhoTrue;
trueParam(:, nParam + 2) = epsilonTrue;
fittedParamAvgIter(:, nParam + 1) = rhoFitted; % transformed rho
fittedParamAvgIter(:, nParam + 2) = epsilonFitted; % transformed epsilon

% ----------------------------------------------------------------------- %
%% 04c) Detect outliers:

iParam = 5; % select parameter

trueParam(:, iParam)
plot(trueParam(:, iParam))

% M05: subject 29 is outlier
trueParam(29, 5)
fittedParamAvgIter(29, 5)

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 05a) Select subjects:

fprintf('Select subjects\n');

invalidSubs     = []; % no outliers, all subjects
% invalidSubs     = [29]; % outlier in M5
% invalidSubs     = [3 10]; % outlier in M7
% invalidSubs     = [11 12 15 23 25 26 30]; % outliers in TAfT and fMRI
fprintf('Exclude subjects %s\n',num2str(invalidSubs));
validSubs       = setdiff(1:nSub, invalidSubs);

% ----------------------------------------------------------------------- %
%% 05b) Correlate true vs. fitted parameters:

for iParam = 1:nParam
    fprintf('Correlation parameter %d: Pearson r = %.03f,  Spearman r = %.03f\n', ....
        iParam, ...
        corr(trueParam(validSubs, iParam), fittedParamAvgIter(validSubs, iParam), 'Type', 'Pearson'), ...
        corr(trueParam(validSubs, iParam), fittedParamAvgIter(validSubs, iParam), 'Type', 'Spearman'))
end

% Save source data files:
fullFileName    = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS05_trueParam.csv');
csvwrite(fullFileName, trueParam(validSubs, :));
fullFileName    = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS05_fittedParam.csv');
csvwrite(fullFileName, fittedParamAvgIter(validSubs, :));

% ----------------------------------------------------------------------- %
%% 05c) Plot true vs. fitted parameters:

% Index parameter names:
param_names{1} = {'\rho', '\epsilon'}; % M1
param_names{2} = {'\rho', '\epsilon', 'Go'}; % M2
param_names{3} = {'\rho', '\epsilon', 'Go', '\pi'}; % M3
param_names{4} = {'\rho', '\epsilon', 'Go', '\kappa'}; % M4
param_names{5} = {'\rho', '\epsilon', 'Go', '\pi', '\kappa'}; % M5
param_names{6} = {'\rho', '\epsilon', 'Go', '\pi', '\epsilon_{SalGo}'}; % M6
param_names{7} = {'\rho', '\epsilon', 'Go', '\pi', '\kappa', '\phi'}; % M7
param_names{8} = {'\rho', '\epsilon', 'Go', '\pi', '\kappa', '\phi_{WIN}', '\phi_{AVOID}'}; % M8
param_names{9} = {'\rho', '\epsilon', 'Go', '\pi', '\kappa', '\eta', '\phi'}; % M9

% Consider excluding subject 29 from regression.

% Choose parameter:
iParam          = 5; % 1-5 regular untransformed parameters

% Plotting settings:
fontSize        = 36;
lineWidth       = 4;
markerSize      = 35; % 25

% Select variables:
x = trueParam(:, iParam);
y = fittedParamAvgIter(:, iParam);

% Start plot:
close all
figure('Position', [100 100 900 800]); hold on

% Plot scatterplot:
plot(x, y, 'k.', 'linewidth', lineWidth, 'markerSize', markerSize)

% Add regression line:
b = polyfit(x(validSubs), y(validSubs), 1); % fit regression
Yhat = polyval(b, x(validSubs)); % predict
plot(x(validSubs), Yhat, 'r-', 'linewidth', lineWidth)

% Add identity line:
plot([-100 100], [-100 100], 'k--', 'linewidth', 1)

% Axis limits:
minLim = min([x; y]);
maxLim = max([x; y]);
axisLim = [minLim - 0.1*maxLim maxLim + 0.1*maxLim];
xlim(axisLim);
ylim(axisLim);

% Settings:
set(gca, 'fontsize', fontSize, 'linewidth', lineWidth);
xlabel('True parameter');
ylabel('Fitted parameter');

% Save:
figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/paramRecovery_%s_mod%02d_param%02d', ...
    parType, iMod, iParam));
if ~isempty(invalidSubs); figName = sprintf('%s_without%s', figName, strjoin(string(invalidSubs), '_')); end
saveas(gcf, [figName '.png']);

% Close:
pause(1)
close gcf


% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% 06a) Compute biased learning rates:

biasepsTrue     = nan(nSub, 2);
biasepsFitted   = nan(nSub, 2);

for iSub = 1:nSub
    if epsilonTrue(iSub) < .5 % If default learning rate below 0.5
      biasepsTrue(iSub, 2) = 1./(1+exp(-(trueParam(iSub, 2) - trueParam(iSub, 5)))); % Negative bias (Punishment after NoGo): subtract untransformed bias from untransformed epsilon, then transform
      biasepsTrue(iSub, 1) = 2.*epsilonTrue(iSub) - biasepsTrue(iSub, 2);              % Positive bias (Reward after Go): take difference transformed epsilon and transformed negative bias, add to transformed epsilon (= 2*transformed epsilon - transformed negative bias)
     else           % If default learning rate above 0.5
      biasepsTrue(iSub, 1) = 1./(1+exp(-(trueParam(iSub, 2) + trueParam(iSub, 5)))); % Positive bias (Reward after Go): add untransformed bias to untransformed epsilon, then transform
      biasepsTrue(iSub, 2) = 2.*epsilonTrue(iSub) - biasepsTrue(iSub, 1);              % Negative bias (Punishment after NoGo): take difference transformed epsilon and transformed positive bias, substract from transformed epsilon (= 2*transformed epsilon - transformed positive bias)
    end
    if epsilonFitted(iSub) < .5 % If default learning rate below 0.5
      biasepsFitted(iSub, 2) = 1./(1+exp(-(fittedParamAvgIter(iSub, 2) - fittedParamAvgIter(iSub, 5)))); % Negative bias (Punishment after NoGo): subtract untransformed bias from untransformed epsilon, then transform
      biasepsFitted(iSub, 1) = 2.*epsilonFitted(iSub) - biasepsFitted(iSub, 2);              % Positive bias (Reward after Go): take difference transformed epsilon and transformed negative bias, add to transformed epsilon (= 2*transformed epsilon - transformed negative bias)
    else
      biasepsFitted(iSub, 1) = 1./(1+exp(-(fittedParamAvgIter(iSub, 2) + fittedParamAvgIter(iSub, 5)))); % Positive bias (Reward after Go): add untransformed bias to untransformed epsilon, then transform
      biasepsFitted(iSub, 2) = 2.*epsilonFitted(iSub) - biasepsFitted(iSub, 1);              % Negative bias (Punishment after NoGo): take difference transformed epsilon and transformed positive bias, substract from transformed epsilon (= 2*transformed epsilon - transformed positive bias)
    end
end

% ----------------------------------------------------------------------- %
%% 06b) Select valid subjects:

invalidSubs     = []; % no outliers, all subjects
% invalidSubs     = [29]; % outlier in M5
fprintf('Exclude subjects %s\n',num2str(invalidSubs));
validSubs       = setdiff(1:nSub, invalidSubs);

% ----------------------------------------------------------------------- %
%% 06c) Plot biased learning rates:

% Exclude subject 29 from regression:

% Choose parameter:
iParam          = 2; % 1 or 2 for rewGo or punNoGo

% Plotting settings:
fontSize        = 36;
lineWidth       = 4;
markerSize      = 35; % 25

% Select variables:
x   = biasepsTrue(:, iParam);
y   = biasepsFitted(:, iParam);

% Start plot:
close all
% figure('Position', [100 100 900 600]); hold on
figure('Position', [100 100 900 800]); hold on

% Plot scatterplot:
plot(x, y, 'k.', 'linewidth', lineWidth, 'markerSize', markerSize)

% Add regression line:
b = polyfit(x(validSubs), y(validSubs), 1); % fit regression
Yhat = polyval(b, x(validSubs)); % predict
plot(x(validSubs), Yhat, 'r-', 'linewidth', lineWidth)

% Add identity line:
plot([-100 100], [-100 100], 'k--', 'linewidth', 1)

% Axis limits:
minLim = min([x; y]);
maxLim = max([x; y]);
axisLim = [minLim - 0.1*maxLim maxLim + 0.1*maxLim];
xlim(axisLim);
ylim(axisLim);

% Settings:
set(gca, 'fontsize', fontSize, 'linewidth', lineWidth);
xlabel('True parameter');
ylabel('Fitted parameter');

% Save:
figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/paramRecovery_%s_mod%02d_biaseps%02d', ...
    parType, iMod, iParam));
if ~isempty(invalidSubs); figName = sprintf('%s_without%s', figName, strjoin(string(invalidSubs), '_')); end
saveas(gcf, [figName '.png']);
% saveas(gcf, [figName '.eps'], 'epsc');

fprintf('Correlation parameter %d: r = %.03f\n', ....
    iParam, corr(biasepsTrue(validSubs, iParam), biasepsFitted(validSubs, iParam)))

% Close:
pause(1)
close gcf

% END OF FILE.
