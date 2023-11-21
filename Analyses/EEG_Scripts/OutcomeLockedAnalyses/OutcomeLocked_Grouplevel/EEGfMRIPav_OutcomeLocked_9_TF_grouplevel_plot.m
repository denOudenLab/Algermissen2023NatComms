% EEGfMRIPav_OutcomeLocked_9_TF_grouplevel_plot

% This is an interactive script---execute it step-by-step.
% Set the appropriate job settings for loading TF data aggregated across
% subjects, then create plots.
% - two-line plot contrasting two conditions of cotnrast
% - multi-line plot contrasting all conditions within data set
% - Topoplot contrasting two conditions of contrast
% - TF plot contrasting two conditions of contrast
% - these plots on the condition averaged signal
% - these plots for each subject separately
% 
% EXPLANATION OF SETTINGS:
% rootDir               = string, root directory of project.
% job                   = cell, created via TF_update_job.m, needs at least
% fields:
%   .nSub               = integer, number of subjects.
%   .dirs               = cell, directories
%   .sub2exclude        = numeric vector, numbers of subjects to exclude.
%   .lockSettings       = string, type of event-locking, 'stimlocked' or
%   'resplocked'.
%   baselineSettings  = string, type of baseline correction, 'trend'
%   (regression-based), 'ft' (with Fieldtrip's ft_timelockbaseline).
%   .baselineTimings    = vector of 1 or 2, timing of baseline correction
%   for which to load data.
%   .responseSettings   = string, type of response setting for which
%   conditions are split, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo).
%   .outcomeSettings    = string, type of outcome coding for which 
%   conditions are split, 'abs' (positive/neutral/negative), 'rel' (positive/negative), 'all' (reward/no reward/ no punishment/ punishment).
%   .TFtype             = type of TF decomposition, 'hanning' or 'morlet'.
%   .nFreq              = numeric, number of frequency bins decomposed
%   (default: 15).
%   .chanArea           = string, area of channels to select, 'frontal',
%   'completefrontal', 'midfrontal', 'midoccipital', 'rightoccipital', 'occipital'.
%   .band               = string, frequency band to select, 'delta', 'theta', 'thetadelta', 'lowalpha', 'verylowalpha', 'alpha', 'beta', 'middlebeta', 'broad'.
%   .outERPcor         = Boolean, read data corrected for outcome-locked 
% ERP (true) or not (false), default false.
%   .outERPresponseSettings   = string, type of response setting for which
%   conditions are split for ERP correction, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo), optional (default: 'Go').
%   .outERPoutcomeSettings    = string, type of outcome coding for which 
%   conditions are split for ERP correction, 'abs' (positive/neutral/negative), 'rel' (positive/negative), 'all' (reward/no reward/ no punishment/ punishment), optional (default: 'all').
%
% OUTPUTS:
% no outputs, just plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

% ----------------------------------------------------------------------- %
%% Set directories:

rootDir     = grouplevel_set_rootDir(); % '/project/3017042.02';

% Set directories and parameters:
dirs        = set_dirs(rootDir);
par         = set_par();

% ----------------------------------------------------------------------- %
%% Initialize job:

job         = []; % initialize empty job

job.nSub    = 36; % necessary for validSubs
job.dirs    = dirs; % add directories

% Subjects to exclude:
job.sub2exclude         = [11 12 23 30]; % exclude 11 12 23 30 because very noisy during trial rejection
% job.sub2exclude         = [11 12 15 23 25 26 30]; % exclusion in TAfT
% job.sub2exclude         = [11 12 18 24 30 34]; % cell sizes < 10 (all NoGo cells, mostly NoGo & chance for reward)

% Number frequencies:
job.nFreqs              = 33; % 33

% Baseline:
job.baselineSettings    = 'trend'; % 'all' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];
% job.baselineTimings     = 0;

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'rel'; % 'rel' or 'abs' or 'all'

job.TFtype              = 'hanning'; % hanning morlet
job.nFreq               = 15; % 64 15 6

% Channel settings:
job.chanArea            = 'midfrontal'; % 'frontal', 'completefrontal', 'midfrontal', 'midoccipital', 'rightoccipital', 'occipital'.

% Band settings:
job.band                = 'theta'; % 'delta' or 'theta' or 'thetadelta' or 'lowalpha' or 'verylowalpha' or 'alpha' or 'beta' or 'middlebeta' or 'broad'.

% Contrast of interest:
job.contrastType        = 'Preferred'; % 'Preferred' or 'Action' or 'GoPreferred' or 'SalientPreferred' or 'SalientAction'. 

% ----------------------------------------------------------------------- %
% ERP-corrected:
job.outERPcor               = false;
job.outERPresponseSettings  = 'Go'; 
job.outERPoutcomeSettings   = 'all';

% ----------------------------------------------------------------------- %
% Load and prepare data:

job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job, data); % prepare generic objects

job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job, data); % prepare objects specific to contrast

% ----------------------------------------------------------------------- %
%% 01A) LINEPLOT ON MAIN CONTRAST:

% Plot:
twoLinePlot(job, data);

% Save:
saveas(gcf,fullfile(dirs.TFplot, sprintf('twoLinePlot_%s.png', ...
    job.plotName)));
pause(3)
close gcf

% ----------------------------------------------------------------------- %
%% 01B) MULTILINEPLOT FOR ALL CUE CONDITIONS:

% Plot:
multiLinePlot(job, data);

% set(gca, 'ylim', [-0.2 1.2]); % wider

saveas(gcf,fullfile(dirs.TFplot, sprintf('multiLinePlot_%s.png', ...
    job.plotName)));
pause(3)
close gcf

% ----------------------------------------------------------------------- %
%% 01C) SINGLE TOPO PLOTS:

% uses job.sigTime for timing
job.sigTime = [0.225 0.425];
% zlim = 0.25;

% Settings:
zlim = 0.5;
% zlim = 1;
% zlim = 2;

% Plot:
singleTopoPlot(job, data, zlim);

% Save:
saveas(gcf,fullfile(dirs.TFplot, sprintf('Topoplot_%s_TOI_%03d_%03dms_zlim%d.png', ...
    job.plotName, job.sigTime(1)*1000, job.sigTime(end)*1000, zlim)))
pause(3)
close gcf

% ----------------------------------------------------------------------- %
%% 01D) MULTIPLE TOPOPLOTS OVER TIME:

% Color scale limits:
zlim = 0.5; % 

nRows = 2; startTime = 0; endTime = 0.9; steps = 0.1; % 2 rows (for theta)
% nRows = 3; startTime = 0; endTime = 1.4; steps = 0.1; % 3 rows (for beta)
% nRows = 2; startTime = 0; endTime = 0.3; steps = 0.1; % 1 rows (for lowalpha)
% nRows = 2; startTime = -1; endTime = -0.1; steps = 0.1; % extraearly

% Plot:
multiTopoPlot(job, data, zlim, nRows, startTime, endTime, steps); % stim-locked 

% Save:
saveas(gcf,fullfile(dirs.TFplot, sprintf('multiTopoPlot_%s_TOI_%03d_%03dms_zlim%d.png', ...
    job.plotName, 1000*startTime, 1000*(endTime+steps),zlim)));
pause(3)
close gcf

% Manual:
% multiTopoPlot(job, data, 1, 2) % stim-locked (zlim = 2 for Go)
% multiTopoPlot(job, data, 2, 2) % stim-locked (zlim = 2 for Go)
% multiTopoPlot(job, data, 1,3,-1,0.4,0.1) % response-locked
% multiTopoPlot(job, data, 2,3,-1,0.4,0.1) % response-locked (zlim = 2 for Go)

% ----------------------------------------------------------------------- %
%% 01E) TF CONTRAST PLOT:

% Color scale limits:
zlim = 0.5;

% Y axis scaling:
% yscale = 'log'; % log, lin
% yscale = 'lin'; % log, lin

% job.TFtiming = [-0.250 1.000];
% job.TFtiming = [-1.000 1.000];

% Plot:
TFplot(job, data, zlim);

% Save:
saveas(gcf, fullfile(dirs.TFplot, sprintf('TFplot_%s_freq_1-15Hz_zlim%d.png', ...
    job.plotName, zlim)));
close gcf

% END OF FILE.