% EEGfMRIPav_OutcomeLocked_9_TF_grouplevel_test

% This is an interactive script---execute it step-by-step.
% Set the appropriate job settings for loading TF data aggregated across
% subjects, then perform tests.
% - cluster-based permutation test with Fieldtrip
% - evaluate stat output object
% - plot stat output object as topoplot or TF plot
% - perform t-test per electrode and averaged over electrodes.
% - export mean of time/frequency/channel range per subject per condition.
% - bar-plot mean of time/frequency/channel range per subject per condition.
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
% By J.C.Swart, 2016/ J. Algermissen, 2018-2023.
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
job.outcomeSettings     = 'all'; % 'rel' or 'abs' or 'all'

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
job.outERPcor           = false;
% job.outERPresponseSettings = 'Go'; 
% job.outERPoutcomeSettings  = 'all';

% ----------------------------------------------------------------------- %
% Load and prepare data:

job             = TF_update_job(job); % initialize job
[job, data]     = TF_load_data(job); % load data
[job, data]     = TF_prepare_generic_data(job, data); % prepare generic objects

job             = TF_update_job(job); % update job
[job, data]     = TF_prepare_contrast_data(job, data); % prepare objects specific to contrast

% ----------------------------------------------------------------------- %
%% PERMUTATION TEST:

% 1) Set whether averaging over channels (yes or no)
% 2) Set whether averaging over frequencies (yes or no)
% 3) Set band

% When determining timing: plotting = 0; otherwise = 1
% plotting = 0;

% ----------------------------------------------------------------------- %
%% 1) Averaging over channels and/ or frequencies:

% Option A: Average over nothing, retrieve both separate channels and 
% frequencies that are part of cluster above threshold:
cfg                 = [];
cfg.avgoverchan     = 'no'; % 
cfg.avgoverfreq     = 'no'; % 
% fprintf('Average neither over channels nor frequencies\n');

% Option B: average over frequencies, only retrieve separate channels that 
% are part of cluster above threshold (--> topoplot):
cfg                 = [];
cfg.avgoverchan     = 'no'; % 
cfg.avgoverfreq     = 'yes'; %
% fprintf('Average over frequencies, but not channels\n');

% Option C: average over channels, only retrieve separate frequencies that 
% are part of cluster above threshold (--> TF plot):
cfg                 = [];
cfg.avgoverchan     = 'yes'; % 
cfg.avgoverfreq     = 'no'; %
% fprintf('Average over channels, but not frequencies\n');

% Option D: average over both channels and frequencies, only obtain p-value
cfg                 = [];
cfg.avgoverchan     = 'yes'; % 
cfg.avgoverfreq     = 'yes'; % 
% fprintf('Average over both channels and frequencies\n');

% ----------------------------------------------------------------------- %
%% 2) Frequencies band to select:

Band = 'theta'; % delta theta thetadelta alpha lowalpha verylowalpha beta middlebeta broad
fprintf('Selected frequency band is %s\n', Band);

% ----------------------------------------------------------------------- %
%% 3) Run permutation test:

rng(70); % set seed for reproducibility of permutation test.

% Initialize neighbours:
cfg_neighb          = [];
cfg_neighb.method   = 'distance';
cfg_neighb.elecfile = 'easycap-M1.txt';
cfg.neighbours      = ft_prepare_neighbours(cfg_neighb);

% Select channels:
cfg.channel         = job.channels; %  {'Fz', 'FCz', 'Cz'};
if length(cfg.channel) == 1 % if single channel, must average:
    cfg.avgoverchan = 'yes';
end

% Select frequency band:
if strcmp(Band, 'broad')
    cfg.frequency       = [1.5 33]; % 
elseif strcmp(Band, 'delta')
    cfg.frequency       = [1.5 4]; % 
elseif strcmp(Band, 'theta')
    cfg.frequency       = [4 8]; %
elseif strcmp(Band, 'deltatheta')
    cfg.frequency       = [1.5 8]; % 
elseif strcmp(Band, 'lowalpha')
    cfg.frequency       = [6 10]; % 
elseif strcmp(Band, 'verylowalpha')
    cfg.frequency       = [7 9]; % 
elseif strcmp(Band, 'alpha')
    cfg.frequency       = [8 13]; % 
elseif strcmp(Band, 'beta')
    cfg.frequency       = [13 30]; % 
else
    error('Band %s not specified', Band);
end

% if plotting==0
%     cfg.avgoverfreq     = 'yes'; % yes; always yes for single channel
% else
%     cfg.avgoverfreq     = 'no'; % for plotting: don't average over frequencies
% end

% Time:
cfg.latency             = [-1 0.7]; % actual test
% cfg.latency             = [0 1.0]; % plotting
% cfg.latency             = [0 1.3]; % extended test for beta

% Permutation test settings:
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05; % threshold at t > 2
% cfg.clusteralpha        = 0.001; % threshold at t > 3.1
cfg.clusterstatistic    = 'maxsum'; % maxsum, maxsize
cfg.minnbchan           = 1;
cfg.tail                = 0; % 0 for two-sided, 1 for one-sided positive, -1 for one-sided negative; see http://www.fieldtriptoolbox.org/reference/ft_statistics_montecarlo/
cfg.clustertail         = cfg.tail; % must correspond to cfg.tail
% cfg.correcttail         = 'alpha';
cfg.alpha               = 0.05; % cut-off for thresholding t-values
cfg.numrandomization    = 500; % 500 by default, can be set higher
subj                    = length(job.validSubs); % number of subjects
design                  = zeros(2, 2*subj); % initialize design matrix
design(1, :)            = repmat(1:subj, [1, 2]); % subject numbers twice (congruent/ incongruent)
design(2, :)            = [ones(1, subj) 2*ones(1, subj)]; % first 1 for each subject (congruent), then 2 for each subject (incongruent)
cfg.design              = design;
cfg.uvar                = 1; % unit variables: (1 for) within 1 subject, not (2 for) across subjects
cfg.ivar                = 2; % independent variables: 2 conditions
fprintf('>>> Perform permutation test for channels %s, %.1f - %.1f Hz, %.1f - %.1f sec\n', ...
    strjoin(job.channels, '/'), cfg.frequency(1), cfg.frequency(end), cfg.latency(1), cfg.latency(end));
% General purpose statement:
[stat]                  = ft_freqstatistics(cfg, data.TF1{job.validSubs}, data.TF2{job.validSubs});

% ----------------------------------------------------------------------- %
%% 4) Evaluate stat:

evaluate_stat(stat, 0.05);

% ----------------------------------------------------------------------- %
%% 5) Plot STATS output as TF plot with Fieldtrip (averaged across channels):
% Needs stat being averaged over channels, but not (!) over frequencies
% Will ignore limited frequency range; always broadband.

zlim                = 2;

% Settings for topoplot:
cfg                 = [];
cfg.xlim            = [-0 1.0]; % -0.25 1.3
cfg.ylim            = [1.5 33];
cfg.zlim            = [zlim*-1 zlim*1];
cfg.maskparameter   = 'mask';
cfg.maskalpha       = 0.05; % 0.025;
cfg.parameter       = 'stat';
ft_singleplotTFR(cfg, stat);

% Save:
saveas(gcf, fullfile(dirs.TFplot, sprintf('TFplot_PermutationCluster_%s_%s_zlim%d.png', ...
    job.plotName, strjoin(sort(job.channels), ''), zlim)));
pause(3)
close gcf

% ----------------------------------------------------------------------- %
%% 6) Plot STATS output as topoplot with Fieldtrip (not averaged across channels):

% Topoplot with clusterplot (needs non-averaged channels):
% needs time or frequency be a singleton dimension (but not channels!).

zlim            = 6;

cfg             = [];
cfg.alpha       = 0.05;
cfg.parameter   = 'stat';
cfg.zlim        = [-1*zlim 1*zlim];
cfg.layout      = job.layout; 
% cfg.freq        = [1 8];
ft_clusterplot(cfg, stat);
saveas(gcf, fullfile(dirs.TFplot, sprintf('Topoplot_permutation_%s_Band_%s_zlim%d.png',...
    job.plotName, Band, zlim)));
pause(3)
close gcf

% ----------------------------------------------------------------------- %
%% T-TEST separately per electrode & averaged over electrodes:

% Which electrodes contribute to this effect? Post-hoc alpha = .05/3 = .017

% Select time:
% ip          = 1; % cluster index
% z           = squeeze(any(any(stat.posclusterslabelmat(:, :, :)==ip, 1), 2))>0;
% selTime     = stat.time(z); % time
% selTime     = MFtimeNeg;
selTime     = job.sigTime;
selTimeIdx  = dsearchn(data.mu.time', selTime'); % indices of start and end time
selTimeIdx  = selTimeIdx(1):selTimeIdx(end); % all indices in between

% Select frequencies:
selFreq     = job.freq';
selFreqIdx  = dsearchn(data.mu.freq', selFreq); % indices of lowest/highest frequency in that range
selFreqIdx  = selFreqIdx(1):selFreqIdx(end); % all indices in between

% Compute differences per subject per channel in selected time and frequency range:
fprintf('Extract mean TF power for %.03f - %0.3f sec., %d-%d Hz:\n', ...
    selTime(1), selTime(end), selFreq(1), selFreq(end))
difData     = nan(job.nSub, length(job.channels));
for iSub = 1:job.nSub % iSub = 1;
    selChanIdx          = find(ismember(data.TF1{iSub}.label, job.channels)); % determine indices of channels per subject
    difData(iSub, :)    = mean(mean(data.TF1{iSub}.powspctrm(selChanIdx, selFreqIdx, selTimeIdx), 3), 2) - ...
        mean(mean(data.TF2{iSub}.powspctrm(selChanIdx, selFreqIdx, selTimeIdx), 3), 2);
end

% Loop over channels, compute t-test:
fprintf('t-test per channel:\n')
for iChan = 1:length(job.channels)
    difVec              = difData(:, iChan);
    [~, P, ~, STATS]    = ttest(difVec);
    fprintf('Channel %s: t(%d) = %.02f, p = %.03f, d = %.02f\n', ...
        job.channels{iChan}, STATS.df, STATS.tstat, P/2, mean(difVec) / std(difVec)); % divide p by 2 because one-sided
end

% ----------------------------------------------------------------------- %
% Overall t-test:

difVec                  = mean(difData, 2); % average over electrodes
[H, P, CI, STATS]       = ttest(difVec);
fprintf('Overall: t(%d) = %.02f, p = %.03f, d = %.02f\n', ...
    STATS.df, STATS.tstat, P/2, mean(difVec) / std(difVec)); % divide p by 2 because one-sided

% ----------------------------------------------------------------------- %
%% BAR PLOT data in selected time/frequency/channel range per subject per condition:

% x-labels need to be adjusted based on particular contrast.

% Determine timing:
selTime         = job.sigTime; % retrieving timing of ROI
if isempty(selTime); error('No selected time'); end

selTimeIdx      = dsearchn(data.mu.time', selTime'); % retrieve edge indices in matrix
selTimeIdx      = selTimeIdx(1):selTimeIdx(end); % insert intermediate indices

% Baseline:
baseTimeIdx     = dsearchn(data.mu.time', 0); % retrieve index at t = 0
baseline        = -1*nanmean(nanmean(nanmean(data.SubCondTime(:, :, baseTimeIdx)))); % retrieve baseline data

% Baseline correction:
subCondMean     = baseline + squeeze(nanmean(data.SubCondTime(:, :, selTimeIdx), 3));
job.nValidSubs  = size(subCondMean, 1);

% Cousineau-Morey SE: average over time, subtract subject overall mean, add sample overall mean
nCond           = job.nCond;
subMean         = nanmean(subCondMean, 2);
grandMean       = nanmean(nanmean(nanmean(data.SubCondTime)));
condMean        = nan(job.nCond, 1); condSE = nan(job.nCond, 1);
for iCond = 1:nCond
    condMean(iCond) = nanmean(subCondMean(:, iCond));
    condSE(iCond)   = nCond/(nCond-1) * nanstd(subCondMean(:, iCond) - ...
        subMean + repmat(grandMean, job.nValidSubs, 1)) ./ sqrt(job.nValidSubs);
end

% ----------------------------------------------------------------------- %
%% Version 1: Action on x-axis, outcome as color:
% Settings:
colMat      = repmat([0 .6 .2; 1 0.63 0.48; 0.49 0.99 0; .8 0 0], job.nCond/4, 1);
xLoc        = [1 2 3 4 5.5 6.5 7.5 8.5];

% yLim        = [-0.6 0.6]; yTick = -1:.2:1;

lineWidth   = 3;
fontSize    = 30;

points 	    = false;

figure('Position', [100 100 1200 800]); hold on
p   = cell(job.nCond, 1);
for iCond = 1:job.nCond
    p{iCond} = bar(xLoc(iCond), condMean(iCond), 0.75, 'FaceColor', colMat(iCond, :));
end
for iCond = 1:job.nCond
    errorbar(xLoc(iCond), condMean(iCond), condSE(iCond), 'k', 'linestyle', 'none', 'Linewidth', lineWidth);
    if points
        s = scatter(repmat(xLoc(iCond), 1, job.nValidSubs), subCondMean(:, iCond)', ...
            [], 'k', 'jitter', 'on', 'jitterAmount', 0.15); hold on % was 0.05
        set(s, 'MarkerEdgeColor', [0.4 0.4 0.4], 'linewidth', 3); % was 1 
        yLim = [round(min(subCondMean(:)), 2) - 0.01 round(max(subCondMean(:)), 2) + 0.01];
        yTick = -10:1:10;
    end
end
xlabel('Performed action', 'FontSize', fontSize, 'FontName', 'Arial', 'Color', [0 0 0]);
ylabel(sprintf('%s Power', job.bandName), 'FontSize', fontSize, 'FontName', 'Arial', 'Color', [0 0 0]);
set(gca, 'FontName', 'Arial', 'FontSize', fontSize, 'ylim', [0 1.5], 'xlim', [0 10],...
    'xtick', [2.5 6.5], 'xticklabel', {'Go', 'NoGo'}, 'ytick',-2:.2:2, 'XColor', [0 0 0], 'YColor', [0 0 0])
legend('Rewarded', 'Non-rewarded', 'Non-punished', 'Punished'); box off
title(sprintf('Action x obtained outcome \n over channels %s, %d-%d ms., %d-%d Hz', ...
    strjoin(job.channels, '/'), 1000*selTime(1), 1000*selTime(end), ...
    job.freq(1), job.freq(end)), 'FontName', 'Arial', ...
    'FontSize', fontSize, 'FontWeight', 'normal', 'Color', [0 0 0]);

if points
    set(gca, 'ylim', [min(subCondMean(:)) max(subCondMean(:))],...
        'yTick',-5:.5:5);
end

% Adjust y-limits:
% set(gca, 'ylim', [-0.5 1.5]); % For theta or lowalpha
% set(gca, 'ylim', [-0.5 1.2]); % For alpha
% set(gca, 'ylim', [-0.3 0.3]); % For small effects, e.g. beta

% Save:
saveas(gcf, fullfile(dirs.TFplot, sprintf('Barplot_%s.png', job.plotName)));
pause(3)
close gcf

% ----------------------------------------------------------------------- %
%% Version 2: Outcome on x-axis, action as color:

condOrder   = [1 5 2 6 3 7 4 8];

% Settings:
colMat      = [1 0 0; 1 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 0 1; 0 0 1];
% colMat      = [0 .6 .2; 0 .6 .2; 0 .6 .2; 0 .6 .2; .8 0 0; .8 0 0; .8 0 0; .8 0 0];
xLoc        = [1 2 3 4 5.5 6.5 7.5 8.5];

% yLim        = [-0.6 0.6]; yTick = -1:.2:1;

lineWidth   = 3;
fontSize    = 30;

points 	    = false;

figure('Position', [100 100 1200 800]); hold on

p   = cell(job.nCond, 1);

iCount      = 0;
for iCond = condOrder
    iCount  = iCount + 1;
    bar(xLoc(iCount), condMean(iCond), 0.75, 'FaceColor', colMat(iCond, :));
end

iCount      = 0;
for iCond = condOrder
    iCount  = iCount + 1;
    errorbar(xLoc(iCount), condMean(iCond), condSE(iCond), ...
        'k', 'linestyle', 'none', 'Linewidth', lineWidth);
    if points
        s = scatter(repmat(xLoc(iCount), 1, job.nValidSubs), subCondMean(:, iCond)', ...
            [], 'k', 'jitter', 'on', 'jitterAmount', 0.15); hold on % was 0.05
        set(s, 'MarkerEdgeColor', [0.4 0.4 0.4], 'linewidth', 3); % was 1 
        yLim = [round(min(subCondMean(:)), 2) - 0.01 round(max(subCondMean(:)), 2) + 0.01];
        yTick = -10:1:10;
    end
end

xlabel('Outcome', 'FontSize', fontSize, 'FontName', 'Arial', 'Color', [0 0 0]);
ylabel(sprintf('%s Power', job.bandName), 'FontSize', fontSize, 'FontName', 'Arial', 'Color', [0 0 0]);

set(gca, 'FontName', 'Arial', 'FontSize', fontSize-10, 'ylim', [-0.5 1.5], 'xlim', [0 10],...
    'xtick', [1.5 3.5 6 8], 'xticklabel', {'Rewarded', 'Non-rewarded', 'Non-punished', 'Punished'}, 'ytick', -2:.2:2,...
    'XColor', [0 0 0], 'YColor', [0 0 0]);

if points
    set(gca, 'ylim', [min(subCondMean(:)) max(subCondMean(:))], ...
        'ytick',-5:.5:5);
end
% Adjust y-limits:
% set(gca, 'ylim', [-0.5 1.5]); % For theta or lowalpha
% set(gca, 'ylim', [-0.5 1.2]); % For alpha
% set(gca, 'ylim', [-0.3 0.3]); % For small effects, e.g. beta

legend('Go', 'NoGo'); box off
title(sprintf('Obtained outcome x action \n over channels %s, %d-%d ms., %d-%d Hz', ...
    strjoin(job.channels, '/'), 1000*selTime(1), 1000*selTime(end), ...
    job.freq(1), job.freq(end)), ...
    'FontName', 'Arial', 'FontSize', fontSize, 'FontWeight', ...
    'normal', 'Color', [0 0 0]);

% Save:
saveas(gcf, fullfile(dirs.TFplot, sprintf('Barplot_%s_fb_act.png', job.plotName)));
% pause(3)
% close gcf

% END OF FUNCTION.