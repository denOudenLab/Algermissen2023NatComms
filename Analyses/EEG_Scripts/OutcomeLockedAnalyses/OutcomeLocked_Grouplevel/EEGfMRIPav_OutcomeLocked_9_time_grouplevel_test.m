% EEGfMRIPav_OutcomeLocked_9_time_grouplevel_test

% This is an interactive script---execute it step-by-step.
% Set the appropriate job settings for loading TF data aggregated across
% subjects, then perform tests.
% - cluster-based permutation test with Fieldtrip.
% - evaluate stat output object.
% - perform t-test per electrode and averaged over electrodes.
% 
% EXPLANATION OF SETTINGS:
% rootDir               = string, root directory of project.
% job                   = cell, created via time_update_job.m, needs at least
% fields:
%   .nSub               = integer, number of subjects.
%   .dirs               = cell, directories
%   .sub2exclude        = numeric vector, numbers of subjects to exclude.
%   baselineSettings  = string, type of baseline correction, 'trend'
%   (regression-based), 'ft' (with Fieldtrip's ft_timelockbaseline).
%   .baselineTimings    = vector of 1 or 2, timing of baseline correction
%   for which to load data.
%   .responseSettings   = string, type of response setting for which
%   conditions are split, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo).
%   .outcomeSettings    = string, type of outcome coding for which 
%   conditions are split, 'abs' (positive/neutral/negative), 'rel' (positive/negative), 'all' (reward/no reward/ no punishment/ punishment).
%   .chanArea           = string, area of channels to select, 'frontal',
%   'completefrontal', 'midfrontal', 'midoccipital', 'rightoccipital', 'occipital'.
%   .band               = string, frequency band to select, 'delta', 'theta', 'thetadelta', 'lowalpha', 'verylowalpha', 'alpha', 'beta', 'middlebeta', 'broad'.
%   .contrastType       = string, contrast to be used: 'Preferred', 'Action',
%   'GoPreferred', 'SalientPreferred', 'SalientAction'.
%   .ROI2use 	        = string, split data by ROI (instead of by conditions:'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1V1Conj'), optional.
% OUTPUTS:
% no outputs, just plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

% ----------------------------------------------------------------------- %
%% Set directories:

% Set root directory:
rootDir     = grouplevel_set_rootDir(); % '/project/3017042.02';

% Set directories and parameters:
dirs        = set_dirs(rootDir);
par         = set_par();

% ----------------------------------------------------------------------- %
%% Initialize job:

job         = []; % initialize empty job

job.nSub    = 36; % necessary for validSubs
job.dirs    = dirs; % add directories

% Data settings:
job.sub2exclude         = [11 12 23 30]; % 03 11 12 23 30 35

job.baselineSettings    = 'trend'; % 'all' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'
job.outcomeSettings     = 'rel'; % 'abs' or 'rel' or 'all'

% Channel settings:
job.chanArea            = 'midfrontal'; % 'frontal', 'completefrontal', 'midfrontal', 'midoccipital', 'rightoccipital', 'occipital'.

% Contrast of interest:
job.contrastType        = 'Preferred'; % 'Preferred' or 'Action' or 'GoPreferred' or 'SalientPreferred' or 'SalientAction'. 

% ----------------------------------------------------------------------- %
% Use ROI:
% job.ROI2use = {'GLM1StriatumConj'}; 
% job.ROI2use = {'GLM1vmPFCConjMan'}; 
% job.ROI2use = {'GLM1V1Conj'};

% ----------------------------------------------------------------------- %
% Load and prepare data:

job         = time_update_job(job); % initialize job

[job, data] = time_load_data(job); % load data
[job, data] = time_prepare_generic_data(job, data); % prepare generic objects

job         = time_update_job(job); % update job
[job, data] = time_prepare_contrast_data(job, data); % prepare objects specific to contrast

% ----------------------------------------------------------------------- %
%% Permutation test:

% Initialize neighbors (needed only of not averaged over channels):
rng(70); % set seed for reproducibility of permutation test.

% Retrieve neighbours:
cfg_neighb          = [];
cfg_neighb.method   = 'distance';
cfg_neighb.elecfile = 'easycap-M1.txt';
cfg.neighbours      = ft_prepare_neighbours(cfg_neighb);
cfg.channel         = job.channels; 

% Average over channels, but not timebins:
cfg.avgoverchan     = 'yes';
% cfg.channel         = {'Fz', 'FCz', 'Cz'};
% cfg.channel       = {'F1', 'F3', 'FCz', 'FC1', 'FC3', 'FC5', 'Cz', 'C1', 'C3'}; % vmPFC time domain
% cfg.channel       = {'Fz', 'F1', 'F2', 'F3', 'F4', 'FCz', 'FC1', 'FC2', 'FC3', 'FC4', 'Cz', 'C1', 'C2', 'C3', 'C4'}; % only midfrontal electrodes
% cfg.channel       = {'AFz', 'AF3', 'AF4', 'Fz', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'FCz', 'FC1', 'FC2', 'FC3', 'FC4', 'FC5', 'FC6', 'Cz', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6'}; % all frontal electrodes
cfg.channel         = job.channels;
cfg.latency         = [0.0 0.7];
cfg.avgovertime     = 'no';

% Other settings:
cfg.parameter                       = 'avg';
cfg.method                          = 'montecarlo';
cfg.statistic                       = 'ft_statfun_depsamplesT';
cfg.alpha                           = 0.05;
cfg.correctm                        = 'cluster';
cfg.correcttail                     = 'prob';
cfg.numrandomization                = 1000;
cfg.design(1, 1:2*job.nValidSubs)   = [ones(1, job.nValidSubs) 2*ones(1, job.nValidSubs)];
cfg.design(2, 1:2*job.nValidSubs)   = [1:job.nValidSubs 1:job.nValidSubs];
cfg.ivar                            = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                            = 2; % the 2nd row in cfg.design contains the subject number
fprintf('Perform permutation test over channels %s, time range %.03f - %.03f\n', ...
    strjoin(string(cfg.channel), ', '), cfg.latency(1), cfg.latency(end));
[stat]                  = ft_timelockstatistics(cfg, data.time1{job.validSubs}, data.time2{job.validSubs});

% ----------------------------------------------------------------------- %
%% Evaluate output:

evaluate_stat(stat, 0.70);

% ----------------------------------------------------------------------- %
%% T-test with Fieldtrip:

rng(70); % set seed for reproducibility of permutation test.

cfg             = [];
cfg.channel     = job.channels; % default 'all'
cfg.avgoverchan = 'yes';
cfg.latency     = [0.2 0.3];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
cfg.design(1, 1:2*job.nValidSubs)   = [ones(1, job.nValidSubs) 2*ones(1, job.nValidSubs)];
cfg.design(2, 1:2*job.nValidSubs)   = [1:job.nValidSubs 1:job.nValidSubs];
cfg.ivar                            = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                            = 2; % the 2nd row in cfg.design contains the subject number
stat = ft_timelockstatistics(cfg, data.time1{job.validSubs}, data.time2{job.validSubs}); 
fprintf('T-test of %s vs. %s in range %.03f - %.03f seconds, channels %s:\nt(%d) = %.03f, p = %.03f\n', ...
    job.twoLineLabels{1}, job.twoLineLabels{2}, cfg.latency(1), cfg.latency(end), strjoin(string(cfg.channel), '/'), ...
    stat.df, stat.stat, stat.prob);

% ----------------------------------------------------------------------- %
%% T-test in selected time range for selected channels:

% Select channels to average over:
selChans    = {'Fz', 'FCz', 'Cz'};

% Select timerange to average over:
selTime     = [0.2 0.3];

selTimeIdx  = dsearchn(data.mu.time', selTime');
selTimeIdx  = selTimeIdx(1):selTimeIdx(2);

% Extract mean signal in defined time/channel range for each subject:
difVec      = nan(size(data.ERPdata, 1), 1); % Initialize empty vector
for iSub = job.validSubs
    selChanIdx = find(ismember(data.ERPdata{iSub, 1}.label, selChans)); % determine indices of channels for each subject
    % Define contrast per hand:
    % Average over time (2) and channels (1)
    % Preferred if responseSettings = Go and outcomeSettings = rel (only 4 conditions):
    difVec(iSub) = mean(mean( ...
        data.ERPdata{iSub, 1}.avg(selChanIdx, selTimeIdx) + ...
        data.ERPdata{iSub, 3}.avg(selChanIdx, selTimeIdx) - ...
        data.ERPdata{iSub, 2}.avg(selChanIdx, selTimeIdx) - ...
        data.ERPdata{iSub, 4}.avg(selChanIdx, selTimeIdx), 2), 1);
end

% Perform t-test:
[~, p, ~, stats] = ttest(difVec); % only electrode [1 2] show the effect significantly.
fprintf('T-test of %s vs. %s in range %.03f - %.03f seconds, channels %s:\nt(%d) = %.03f, p = %.03f\n', ...
    job.twoLineLabels{1}, job.twoLineLabels{2}, selTime(1), selTime(end), strjoin(string(selChans), '/'), ...
    stats.df, stats.tstat, p);

% END OF FILE.