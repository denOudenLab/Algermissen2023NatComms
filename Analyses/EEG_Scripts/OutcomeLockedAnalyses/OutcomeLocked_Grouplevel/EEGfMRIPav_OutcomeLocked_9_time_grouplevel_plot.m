% EEGfMRIPav_OutcomeLocked_9_time_grouplevel_plot

% This is an interactive script---execute it step-by-step.
% Set the appropriate job settings for loading TF data aggregated across
% subjects, then create plots.
% - two-line plot contrasting two conditions of cotnrast
% - multi-line plot contrasting all conditions within data set
% - Topoplot contrasting two conditions of contrast
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

% ------------------------------------------------------------------ %
%% Set directories:

% Set root directory:
rootDir = grouplevel_set_rootDir(); % '/project/3017042.02';

% Set directories and parameters:
dirs        = set_dirs(rootDir);
par         = set_par();

% ------------------------------------------------------------------ %
%% Initialize job:

job         = []; % initialize empty job

job.nSub    = 36; % necessary for validSubs
job.dirs    = dirs; % add directories

% Data settings:
job.sub2exclude         = [11 12 23 30]; % exclude 11 12 23 30 because very noisy during trial rejection
% job.sub2exclude         = [11 12 18 23 24 30 34]; % exclude subjects with < 10 trials in some cells (all NoGo cells)

job.baselineSettings    = 'trend'; % 'all' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'all'; % 'abs' or 'rel' or 'all'

% Channel settings:
job.chanArea            = 'midfrontal'; % 'frontal', 'completefrontal', 'midfrontal', 'midoccipital', 'rightoccipital', 'occipital'.

% Contrast of interest:
job.contrastType        = 'Preferred'; % 'Preferred' or 'Action' or 'GoPreferred' or 'SalientPreferred' or 'SalientAction'. 

% ----------------------------------------------------------------------- %
% Use ROI:
job.ROI2use = {'GLM1StriatumConj'}; 
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
%% 1) TWOLINEPLOT ON MAIN CONTRAST:

% job = rmfield(job, 'sigTime'); % remove sigTime if inappropriate
twoLinePlot(job, data);

% set(gca, 'xtick',0:0.05:1, 'xtickLabel',0:0.05:1, 'fontsize', 12)

saveas(gcf,fullfile(dirs.timeplot, sprintf('twoLinePlot_%s.png', job.plotName)));
pause(3)
close gcf

% ----------------------------------------------------------------------- %
%% 2) MULTILINEPLOT FOR ALL CUE CONDITIONS:

% job = rmfield(job, 'sigTime');

multiLinePlot(job, data);

% set(gca, 'ylim', [-0.025 0.035]); % for 'all'
% set(gca, 'xtick', 0:0.05:1, 'xtickLabel', 0:0.05:1, 'fontsize', 12); % higher resolution x-axis

saveas(gcf,fullfile(dirs.timeplot, sprintf('multiLineplot_%s.png', job.plotName)));
close gcf

% ----------------------------------------------------------------------- %
%% 3A) SINGLE TOPO PLOTS:

% zlim:
zlim = 0.015; % task effects
% zlim = 0.010; % BOLD effects
% zlim = 0.005; % BOLD effects even weaker

% Timing:
% Preferred vs. non-preferred:
% Positive cluster 1: 0.248 - 0.292 sec., p = 0.042
% Negative cluster 1: 0.341 - 0.412 sec., p = 0.006
% job.sigTime = [0.248 0.292]; % N1/FRN
job.sigTime = [0.341 0.412]; % N2/P3

singleTopoPlot(job, data,zlim);

% Save:
saveas(gcf,fullfile(dirs.timeplot, sprintf('Topoplot_%s_TOI_%03d_%03dms_zlim%.03f.png', ...
    job.plotName, job.sigTime(1)*1000, job.sigTime(end)*1000, zlim)));
pause(3)
close gcf

% ----------------------------------------------------------------------- %
%% 3B) MULTIPLE TOPO PLOTS OVER TIME:

zlim = 0.015; % task effects
% zlim = 0.010; % BOLD effects
% zlim = 0.005; % BOLD effects (stronger)

% Settings for grid:
% nRows = 2; startTime = 0; endTime = 0.9; steps = 0.1; % 2 rows (for theta)
% nRows = 3; startTime = 0; endTime = 1.4; steps = 0.1; % 3 rows (for beta)
nRows = 2; startTime = 0; endTime = 0.45; steps = 0.05; % 2 rows (for BOLD)

% Plot:
multiTopoPlot(job, data, zlim, nRows, startTime, endTime, steps); % 

% Save:
saveas(gcf,fullfile(dirs.timeplot, sprintf('multiTopoPlot_%s_TOISequence_%d-%dms_zlim%d.jpg',...
    job.plotName, 1000*startTime, 1000*(endTime+steps), zlim)));
pause(3)
close gcf

% END OF FILE.