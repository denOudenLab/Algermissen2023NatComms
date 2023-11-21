% EEGfMRIPav_cbm_fit.m

% This is an interactive script---execute it step-by-step.
% It fits a series of computational reinforcement learning models using the
% CBM toolbox and evaluates them.
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
% cd /project/3017042.02/Analyses/CBM_Scripts/

% clear all; close all; clc

% ----------------------------------------------------------------------- %
%% 00a) Settings:

nMod        = 9; % number of models to fit
dirs        = [];
dirs.root   = '/project/3017042.02';

fprintf('Fit %d models\n', nMod);

% ----------------------------------------------------------------------- %
%% 00b) Add paths:

fprintf('Add paths\n');

% Add paths:
addpath('/home/action/johalg/cbm-master/codes'); % CBM toolbox
addpath(fullfile(dirs.root, '/Analyses/CBM_Scripts/models')); % models

% ----------------------------------------------------------------------- %
%% 00c) Directories:

fprintf('Initialize directories\n');

dir.results = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');

dir.lap     = fullfile(dir.results, 'LAP_Results');
if ~exist(dir.lap, 'dir'); mkdir(dir.lap); end

dir.hbi     = fullfile(dir.results, 'HBI_Results');
if ~exist(dir.hbi, 'dir'); mkdir(dir.hbi); end

% ----------------------------------------------------------------------- %
%% 00d) Load and extract data, priors, output name:

fprintf('Load data\n');
inputFile   = fullfile(dir.results, 'EEGfMRIPav_cbm_inputData.mat');    
fdata       = load(inputFile); % only contains .data
data        = fdata.data; % extract .data
nSub        = length(data);

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

% Output names:
fprintf('Initialize output file names\n')
fname_mod = struct([]);
for iMod = 1:nMod
    fname_mod{iMod} = fullfile(dir.lap, sprintf('lap_mod%02d.mat', iMod));
end

% ----------------------------------------------------------------------- %
%% 00e) Check models in dry run:

fprintf('Test models (dry run)\n')
subj1 = data{1};

% a) Random parameter values:
for iMod = 1:nMod
    parameters = randn(1, 7);
    F1 = eval(sprintf('EEGfMRIPav_cbm_mod%02d(parameters, subj1)', iMod));
    fprintf('Model %02d: loglik = %f\n', iMod, F1);
end

% b) Extreme parameter values:
% for iMod = 1:nMod
%     parameters = [-10 10 -10 10 10 10 10];
%     F1 = eval(sprintf('EEGfMRIPav_cbm_mod%02d(parameters, subj1)', iMod));
%     fprintf('Model %02d: loglik = %f\n', iMod,F1);
% end

% c) Manually:
% parameters = randn(1, 2);
% F1 = EEGfMRIPav_cbm_mod01(parameters, subj1)

% parameters = [-10 10];
% F1 = EEGfMRIPav_cbm_mod01(parameters, subj1)

% ----------------------------------------------------------------------- %
%% 1a) LaPlace approximation (cbm_lap):

% All models are fit non-hierarchically.

% Fit for each model separately:
for iMod = 1:nMod
    fprintf('Fit model %02d with LaPlace approximation\n', iMod)
    % Format data, model, prior, output file
    cbm_lap(data, eval(sprintf('@EEGfMRIPav_cbm_mod%02d', iMod)), priors{iMod}, fname_mod{iMod});
end

% ----------------------------------------------------------------------- %
%% 1b) Evaluate LAP fit:

% Output names:
fprintf('Initialize output file names\n')
fname_mod = struct([]);
for iMod = 1:nMod
    fname_mod{iMod} = fullfile(dir.lap, sprintf('lap_mod%02d.mat', iMod));
end

% Retrieve log model-evidence per model per subject:
logModEvi = nan(nSub, nMod);
for iMod = 1:nMod % nMod = 5;
    fname = load(fname_mod{iMod});
    cbm   = fname.cbm;
    logModEvi(:, iMod) = cbm.output.log_evidence;
end

% Mean per model:
fprintf('Mean log-model evidence per model:\n');
mean(logModEvi, 1)
round(mean(logModEvi, 1))

% Sum per model:
fprintf('Summed log-model evidence per model:\n');
round(sum(logModEvi, 1))

% Plot log model-evidence per subject:
figure('units', 'normalized', 'outerposition', [0 0 1 1]); hold on
legendVec = [];
for iMod = 1:nMod
    plot(1:nSub, logModEvi(:, iMod), 'linewidth', 2);
    legendVec = [legendVec; sprintf('Model %02d', iMod)];
end
xlabel('Subjects'); ylabel('Log model evidence');
legend(legendVec);

% t-test between 2 models:
iMod1 = 1; iMod2 = 5;
fprintf('Test log model-evidence of models %d and %d against each other:\n', iMod1, iMod2);
[~, p, ~, STATS] = ttest(logModEvi(:, iMod1), logModEvi(:, iMod2));
fprintf('t(%d) = %.3f, p = %.03f\n', STATS.df, STATS.tstat, p);

% ---------------------------------------- %
% Model frequency:
modRange    = 1:nMod;
modFreq     = nan(nSub, 1);

% Per subject:
for iSub = 1:nSub % iSub = 1
    [~,I]   = max(logModEvi(iSub, modRange));
    modFreq(iSub) = I;
end
% Per model:
fprintf('Model frequency per model:\n');
for iMod = 1:length(modRange)
    idx     = find(modFreq == iMod);
    fprintf('Model M%02d is best for %02d subjects: %s\n', ...
        modRange(iMod), length(idx), mat2str(idx));
end

% ----------------------------------------------------------------------- %
%% 2a) Prepare and fit Hierarchical Bayesian inference (cbm_hbi) per SINGLE model:

% 1st input: data for all subjects:
inputFile = fullfile(dir.results, 'EEGfMRIPav_cbm_inputData.mat');    
fdata = load(inputFile); % only contains .data
data  = fdata.data; % extract .data

% 2nd input: priors--see above

% Other inputs: set within loop:
for iMod = 1:nMod
    fprintf('Fit model %02d with HBI (singular model)\n', iMod)
    % 3rd input: models
    models      = {eval(sprintf('@EEGfMRIPav_cbm_mod%02d', iMod))};
    fcbm_maps   = {fullfile(dir.lap, sprintf('lap_mod%02d.mat', iMod))};
    % 4th input: lap output file
    fname_hbi   = {fullfile(dir.hbi, sprintf('hbi_mod_%02d.mat', iMod))};
    % Fit:
    cbm_hbi(data, models, fcbm_maps, fname_hbi);
end

% ----------------------------------------------------------------------- %
%% 3a) Prepare Hierarchical Bayesian inference (cbm_hbi) across models:

% 1st input: data for all subjects:
inputFile   = fullfile(dir.results, 'EEGfMRIPav_cbm_inputData.mat');    
fdata       = load(inputFile); % only contains .data
data        = fdata.data; % extract .data

fprintf('Select subjects\n');

% Select subjects:
invalidSubs     = []; % outliers in TAfT and fMRI
% invalidSubs     = [11 12 15 23 25 26 30]; % outliers in TAfT and fMRI
fprintf('Exclude subjects %s\n', num2str(invalidSubs));
validSubs       = setdiff(1:nSub, invalidSubs);
nSubValid       = length(validSubs);
if ~isempty(invalidSubs)
    data = data(validSubs); 
end

% 2nd input: a cell input containing function handle to models

% 3rd input: another cell input containing file-address to files saved by cbm_lap
% All models:
modVec = 1:nMod;

fprintf('Prepare HBI for comparing models %s\n', strjoin({num2str(modVec)}, ' '));

% Create names of input files:
iCount      = 0;
models      = cell(length(modVec), 1);
fcbm_maps   = cell(length(modVec), 1);
for iMod = modVec
    iCount = iCount + 1;
    models{iCount} = eval(sprintf('@EEGfMRIPav_cbm_mod%02d', iMod));
    fcbm_maps{iCount} = fullfile(dir.lap, sprintf('lap_mod%02d.mat', iMod));
end

% 4th input: a file address for saving the output
% All models named explicitly:
fname_hbi = fullfile(dir.hbi, sprintf('hbi_mod%s', num2str(modVec, '_%02d'))); % single model
fprintf('Output file will be %s\n',fname_hbi)

if ~isempty(invalidSubs); fname_hbi = [fname_hbi '_without7']; end
fname_hbi = [fname_hbi '.mat'];

% Fit:
fprintf('Fit models %s with HBI \n', num2str(modVec, '_%02d'))
cbm_hbi(data,models, fcbm_maps, fname_hbi);
fprintf('Finished fitting :-)\n')

% ----------------------------------------------------------------------- %
%% 3b) Additionally run protected exceedance probability (including null model):

fprintf('Re-run models %s with HBI including null model\n', num2str(modVec, '_%02d'))
cbm_hbi_null(data,fname_hbi);
fprintf('Finished fitting including null model:-)\n')

% Evaluate:
f_hbi   = load(fname_hbi);
cbm     = f_hbi.cbm;
mf      = cbm.output.model_frequency;
xp      = cbm.output.exceedance_prob;
pxp     = cbm.output.protected_exceedance_prob;
xp_dif  = xp - pxp; % no difference
fprintf('Model frequency                 : %s\n', num2str(mf, '%.02f '));
fprintf('Exceedance probability          : %s\n', num2str(xp, '%.02f '));
fprintf('Protected Exceedance probability: %s\n', num2str(pxp, '%.02f '));
fprintf('Difference                      : %s\n', num2str(xp_dif, '%.02f '));
fprintf('Finished :-)\n')

% ----------------------------------------------------------------------- %
%% 3c) Evaluate HBI fit:

% Create output name:
modVec      = 1:6;
fname_hbi   = fullfile(dir.hbi, sprintf('hbi_mod%s.mat', num2str(modVec, '_%02d'))); % all models

% ----------------------------- %
% Load model:
f_hbi       = load(fname_hbi);
cbm         = f_hbi.cbm;
nSub        = size(cbm.output.parameters{1}, 1);

% ----------------------------- %
% Model frequency:
fprintf('Model frequency\n')
cbm.output.model_frequency

% ----------------------------- %
% Responsibility per subject:
fprintf('Maximal model responsibility\n')
responsibility  = cbm.output.responsibility;

% Per subject:
subResp     = nan(nSub, 1);
for iSub = 1:nSub
    subMax          = max(responsibility(iSub, :));
    subResp(iSub)   = find(responsibility(iSub, :) == subMax);
end
tabulate(subResp)

% Per model:
for iMod = 1:length(modVec)
    idx     = find(subResp == iMod);
    fprintf('Model M%02d is most responsible for %02d subjects: %s\n', ...
        modVec(iMod), length(idx), mat2str(idx));
end

% ----------------------------- %
% Exceedance probability:
fprintf('Exceedance probability:\n')
cbm.output.exceedance_prob

% ----------------------------- %
% Protected exceedance probability:
fprintf('Protected exceedance probability:\n')
cbm.output.protected_exceedance_prob

% END OF FILE.
