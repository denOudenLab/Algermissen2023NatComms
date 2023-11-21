% FigureS20.m

% Plots for Supplementary Figure S20.
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

% Add paths to helper files:
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel'));
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT'));

% Add Fieldtrip:
addpath /home/common/matlab/fieldtrip;
ft_defaults;

% ----------------------------------------------------------------------- %
%% Load saved regression coefficients for 7 fMRI ROIs:

EEGdomain           = 'time';

% 1) ROIs:
ROIs2use            = {'GLM1StriatumConj', 'GLM1ACCConjMan', ...
    'GLM1LeftMotorConj', 'GLM1vmPFCConjMan', 'GLM1PCCConj', ...
    'GLM1LeftITGConj', 'GLM1V1Conj'}; 
if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end

% 2) Behavioral regressors:
behav2use           = {'Updatestd', 'Updatedif'};

% 3) Perform only on selected trials:
selTrials           = 'all';

[job, dirs, betas]  = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
%% Figure S20A: Multiple lines with highlighted patches in one plot: striatum, dACC, pgACC, PCC

% a) Channels:
selChans    = {'Fz', 'FCz', 'Cz'}; % Jenn's a-priori selection

% b) ROIs:
allROIs     = [1 2 4 5]; % Primary value regions: striatum, ACC, vmPFC, PCC

% c) ROI legend names:
legendNames = {'Striatum', 'ACC', 'vmPFC', 'PCC'};

% d) Colors:
colMat      = [0.75 0.75 0; 1 0 0; 0 0 1;0 0.75 0.75]; % 4 colors: Striatum, ACC, vmPFC, PCC

% e) General settings:
thresh      = 2.0;
nP          = 1000;
ylim        = 5;
lineWidth   = 4;
fontSize    = 32;

% Start figure:
p           = cell(length(allROIs), 1); % initialize p
saveMat     = nan(701, length(allROIs));
close gcf
figure('Position', [100 100 1200 800], 'Color', 'white'); 

iCondi      =  0;
for iROI = allROIs
    
    % Compute:
    iCondi          = iCondi + 1;
    [sortBetas,~]   = taft_postprocess_time_selectData(job, betas, iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
    c               = [];
    cfg             = [];
    cfg.channel     = selChans;
    cfg.avgoverchan = 'yes';
    cfg.latency     = [0 0.7]; % time where ERPs occur
    cfg.avgovertime = 'no';
    
    fprintf('Regressor %s: Average over channels %s, reshape into 3D\n', char(job.regNames(iROI)), strjoin(selChans, '/'));
    chanBetas       = cell(length(sortBetas), 1);
    for iSub = 1:length(sortBetas) % iSub = 1;
        chanBetas{iSub} = ft_selectdata(cfg, sortBetas{iSub}); % average over selected channels
        c(iSub, :, :)     = chanBetas{iSub}.avg; % bring into 4-D, with subject as first dimension
    end
    
    [corrp, tg]  = clustertf(c, thresh, nP);
    
    % Plot:
    endTime     = cfg.latency(end);
    timeVec     = sortBetas{1}.time(1:find(sortBetas{1}.time == endTime));
    Tvec        = squeeze(tg);
    pvec        = squeeze(corrp);
    p{iCondi}   = plot(chanBetas{1}.time, Tvec, 'color', colMat(iCondi, :), 'linewidth', lineWidth); hold on
    saveMat(:, iCondi) = Tvec; % save
    
    if any(pvec < .05 & Tvec < -2) % negative
        sigTimeIdx  = (Tvec < -2 & pvec < 0.05); % negative
        sigTime     = find(sigTimeIdx==1);
        fprintf('Significant negative cluster from %.03f - %.03f sec.\n', timeVec(sigTime(1)), timeVec(sigTime(end)));
        patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))], [repmat(-2, 1, sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')], colMat(iCondi, :), 'EdgeColor', 'none'); hold on
    end
    
    if any(pvec < .05 & Tvec > 2) % positive
        sigTimeIdx  = (Tvec > 2) & (pvec < 0.05); % positive
        sigTime     = find(sigTimeIdx==1);
        fprintf('Significant positive cluster from %.03f - %.03f sec.\n', timeVec(sigTime(1)), timeVec(sigTime(end)));
        patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))], [repmat(2, 1, sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')], colMat(iCondi, :), 'EdgeColor', 'none'); hold on
    end
end

% Horizontal lines:
nTime   = length(chanBetas{iSub}.time);
plot(chanBetas{iSub}.time, repelem(2, nTime), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 2);
plot(chanBetas{iSub}.time, repelem(-2, nTime), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 2);

% Final settings:
set(gca, 'xlim', cfg.latency, ...
    'xtick', [-1 -0.5 -0.25 0 0.25 0.5 0.7 1], ...
    'xtickLabel', {'-1', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
    'ylim', [-ylim ylim], 'ytick',-10:2:10, ...
    'fontsize', fontSize, 'LineWidth', lineWidth);
xlabel('Time (ms)', 'fontweight', 'bold', 'fontsize', fontSize); 
ylabel('T-values', 'fontweight', 'bold', 'fontsize', fontSize); 
title(sprintf('T-values over %s', strjoin(selChans, '/')));
box off;

% Legend:
legend([p{:}], legendNames, 'linewidth', lineWidth, 'fontsize', fontSize, 'Interpreter', 'tex'); legend boxoff

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS20A.png'));
pause(3);
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS20A.csv');
csvwrite(fullFileName, saveMat);

fprintf('Done :-)\n');

% ----------------------------------------------------------------------- %
%% Figure S20B: Multiple lines with highlighted patches in one plot: left M1, left ITG, V1

% a) Channels:
selChans    = {'Fz', 'FCz', 'Cz'}; % Jenn's a-priori selection

% b) ROIs:
allROIs     = [3 6 7]; % rest: M1, left ITG, V1

% c) ROI legend names:
legendNames = {'Left M1', 'Left ITG', 'V1'};

% d) Colors:
colMat      = [1 0.66 0.4; 0.33 0.87 0.62; 0.30 0.79 0.33]; % 3 colors rest: left M1, left ITG, V

% e) General settings:
thresh      = 2.0;
nP          = 1000;
ylim        = 5;
lineWidth   = 4;
fontSize    = 32;

% Start figure:
p           = cell(length(allROIs), 1); % initialize p
saveMat     = nan(701, length(allROIs));
close gcf
figure('Position', [100 100 1200 800], 'Color', 'white'); 

iCondi      =  0;
for iROI = allROIs
    
    % Compute:
    iCondi          = iCondi + 1;
    [sortBetas,~]   = taft_postprocess_time_selectData(job, betas, iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
    c               = [];
    cfg             = [];
    cfg.channel     = selChans;
    cfg.avgoverchan = 'yes';
    cfg.latency     = [0 0.7]; % time where ERPs occur
    cfg.avgovertime = 'no';
    
    fprintf('Regressor %s: Average over channels %s, reshape into 3D\n', char(job.regNames(iROI)), strjoin(selChans, '/'));
    chanBetas       = cell(length(sortBetas), 1);
    for iSub = 1:length(sortBetas) % iSub = 1;
        chanBetas{iSub} = ft_selectdata(cfg, sortBetas{iSub}); % average over selected channels
        c(iSub, :, :)     = chanBetas{iSub}.avg; % bring into 4-D, with subject as first dimension
    end
    
    [corrp, tg]  = clustertf(c, thresh, nP);
    
    % Plot:
    endTime     = cfg.latency(end);
    timeVec     = sortBetas{1}.time(1:find(sortBetas{1}.time == endTime));
    Tvec        = squeeze(tg);
    pvec        = squeeze(corrp);
    p{iCondi}   = plot(chanBetas{1}.time, Tvec, 'color', colMat(iCondi, :), 'linewidth', lineWidth); hold on
    saveMat(:, iCondi) = Tvec; % save
    
    if any(pvec < .05 & Tvec < -2) % negative
        sigTimeIdx  = (Tvec < -2 & pvec < 0.05); % negative
        sigTime     = find(sigTimeIdx==1);
        fprintf('Significant negative cluster from %.03f - %.03f sec.\n', timeVec(sigTime(1)), timeVec(sigTime(end)));
        patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))], [repmat(-2, 1, sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')], colMat(iCondi, :), 'EdgeColor', 'none'); hold on
    end
    
    if any(pvec < .05 & Tvec > 2) % positive
        sigTimeIdx  = (Tvec > 2) & (pvec < 0.05); % positive
        sigTime     = find(sigTimeIdx==1);
        fprintf('Significant positive cluster from %.03f - %.03f sec.\n', timeVec(sigTime(1)), timeVec(sigTime(end)));
        patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))], [repmat(2, 1, sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')], colMat(iCondi, :), 'EdgeColor', 'none'); hold on
    end
end

% Horizontal lines:
nTime   = length(chanBetas{iSub}.time);
plot(chanBetas{iSub}.time, repelem(2, nTime), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 2);
plot(chanBetas{iSub}.time, repelem(-2, nTime), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 2);

% Final settings:
set(gca, 'xlim', cfg.latency, ...
    'xtick', [-1 -0.5 -0.25 0 0.25 0.5 0.7 1], ...
    'xtickLabel', {'-1', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
    'ylim', [-ylim ylim], 'ytick',-10:2:10, ...
    'fontsize', fontSize, 'LineWidth', lineWidth);
xlabel('Time (ms)', 'fontweight', 'bold', 'fontsize', fontSize); 
ylabel('T-values', 'fontweight', 'bold', 'fontsize', fontSize); 
title(sprintf('T-values over %s', strjoin(selChans, '/')));
box off;

% Legend:
legend([p{:}], legendNames, 'linewidth', lineWidth, 'fontsize', fontSize, 'Interpreter', 'tex'); legend boxoff

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS20B.png'));
pause(3);
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS20B.csv');
csvwrite(fullFileName, saveMat);

fprintf('Done :-)\n');

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Figure S20C: Line plot striatum:

clear job data

% ----------------------------------------------------------------------- %
%% Prepare contrast:

dirs.EEG                = fullfile(dirs.root, 'Log', 'EEG', 'OutcomeLockedResults');
dirs.timegroup          = fullfile(dirs.EEG, 'time_grouplevel');
job.dirs                = dirs;
job.nSub                = 36; % necessary for validSubs
job.sub2exclude         = [11 12 23 30]; % exclude 11 12 23 30 because very noisy during trial rejection
job.nFreqs              = 33; % 33
job.baselineSettings    = 'trend';
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'all'; % 'rel' or 'abs' or 'all'

% ----------------------------------------------------- %
% Channels of interest:
job.chanArea            = 'midfrontal'; % 'midfrontal'
% Contrast of interest:
job.contrastType        = 'Preferred'; 

% ----------------------------------------------------- %
% Use ROI (overwrites contrast):

job.ROI2use = {'GLM1StriatumConj'}; 

if ~exist('data', 'var')
    job         = time_update_job(job);
    [job, data] = time_load_data(job);
    [job, data] = time_prepare_generic_data(job, data);
end

job             = time_update_job(job);
[~, data]       = time_prepare_contrast_data(job, data);

% ----------------------------------------------------------------------- %
%% Create line plot:

% Fixed settings:
fontSize            = 32;
lineWidth1          = 5; % actual lines in plot
lineWidth2          = 3; % all other lines

% Determine baseline (i.e. lowest line at time point zero):
timeIdx             = dsearchn(data.mu.time', [-0.25 0.80]');
timeIdx             = timeIdx(1):timeIdx(end);
iBaseline           = find(round(data.mu.time, 3) == 0); % find baseline via time in sec
baseline            = min(nanmean(data.mat1(:, iBaseline)), nanmean(data.mat2(:, iBaseline)));

% Overwrite colors:
job.twoLineColor    = [1 0 0; 0 0 1];

% Start figure:
p   = cell(2, 1);
figure('Position', [100 100 1200 800], 'Color', 'white'); hold on % DCC screen

% Lines:
% Use Cousineau 2005 method for uncorrected-within-subjects CIs 
% (only correction, no z; i.e. approx SE with only 2 conditions):
p{1} = boundedline(data.mu.time(timeIdx), ...
    nanmean(data.mat1(:, timeIdx)) - baseline, ...
    2 / (2-1) * nanstd(data.mat1(:, timeIdx) - ...
    data.matMean(:, timeIdx) + ... 
    repmat(data.matGrandMean(timeIdx), job.nValidSubs, 1)) ./ ...
    sqrt(job.nValidSubs), ...
    'cmap', job.twoLineColor(1, :), 'alpha', 'transparency', 0.2); %
set(p{1}, 'linewidth', lineWidth1); % adjust linewidth (not possible within boundedline)
p{2} = boundedline(data.mu.time(timeIdx), ...
    nanmean(data.mat2(:, timeIdx)) - baseline, ...
    2 / (2-1) * nanstd(data.mat2(:, timeIdx) - ...
    data.matMean(:, timeIdx) + ...
    repmat(data.matGrandMean(timeIdx), job.nValidSubs, 1)) ./ ...
    sqrt(job.nValidSubs), ...
    'cmap', job.twoLineColor(2, :), 'alpha', 'transparency', 0.2); %
set(p{2}, 'linewidth', lineWidth1); % adjust linewidth (not possible within boundedline)

% Extra vertical lines:
plot([0 0], get(gca,'ylim'), ':k', 'LineWidth', 3); % outcome onset
plot([.7 .7], get(gca,'ylim'), ':k', 'LineWidth', 3); % outcome offset

% x-axis label:
xlabel('Time (ms)', 'fontsize', fontSize, 'fontweight', 'bold', 'linewidth', lineWidth2);

% Retrieve y-axis limits, change later:
yLim    = get(gca, 'ylim'); % retrieve automatic settings, change later
yMinLim = yLim(1); 
yMaxLim = yLim(2); 

% Add x-axis settings and labels:
set(gca, 'xlim', [-0.25 0.8], ...
    'ylim', [yMinLim yMaxLim], ...
    'xtick', [-1.0 -0.75 -0.5 -0.25 0 0.25 0.5 0.70], ...
    'xtickLabel',{'-1000', '-750', '-500', '-250', 'FbOn', '250', '500', 'FbOff'}, ...
    'fontsize', fontSize, 'linewidth', lineWidth2);

% Add title y-label:
ylabel('Amplitude (A/cm²)', 'fontweight', 'bold', 'fontsize', fontSize); 

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS20C.png'));
pause(3)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20C_High.csv'));
csvwrite(fullFileName, data.mat1(:, timeIdx) - baseline);
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20C_Low.csv'));
csvwrite(fullFileName, data.mat2(:, timeIdx)- baseline);


% ----------------------------------------------------------------------- %
%% Figure S20D: Line plot pgACC:

clear job data

% ----------------------------------------------------------------------- %
%% Prepare contrast:

dirs.EEG                = fullfile(dirs.root, 'Log', 'EEG', 'OutcomeLockedResults');
dirs.timegroup          = fullfile(dirs.EEG, 'time_grouplevel');
job.dirs                = dirs;
job.nSub                = 36; % necessary for validSubs
job.sub2exclude         = [11 12 23 30]; % exclude 11 12 23 30 because very noisy during trial rejection
job.nFreqs              = 33; % 33
job.baselineSettings    = 'trend';
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'all'; % 'rel' or 'abs' or 'all'

% ----------------------------------------------------- %
% Channels of interest:
job.chanArea            = 'midfrontal'; % 'midfrontal'
% Contrast of interest:
job.contrastType        = 'Preferred'; 

% ----------------------------------------------------- %
% Use ROI (overwrites contrast):

job.ROI2use = {'GLM1vmPFCConjMan'}; 


if ~exist('data', 'var')
    job         = time_update_job(job);
    [job, data] = time_load_data(job);
    [job, data] = time_prepare_generic_data(job, data);
end

job             = time_update_job(job);
[~, data]       = time_prepare_contrast_data(job, data);

% ----------------------------------------------------------------------- %
%% Create line plot:

% Fixed settings:
fontSize            = 32;
lineWidth1          = 5; % actual lines in plot
lineWidth2          = 3; % all other lines

% Determine baseline (i.e. lowest line at time point zero):
timeIdx             = dsearchn(data.mu.time', [-0.25 0.80]');
timeIdx             = timeIdx(1):timeIdx(end);
iBaseline           = find(round(data.mu.time, 3) == 0); % find baseline via time in sec
baseline            = min(nanmean(data.mat1(:, iBaseline)), nanmean(data.mat2(:, iBaseline)));

% Overwrite colors:
job.twoLineColor    = [1 0 0; 0 0 1];

% Start figure:
p   = cell(2, 1);
figure('Position', [100 100 1200 800], 'Color', 'white'); hold on % DCC screen

% Lines:
% Use Cousineau 2005 method for uncorrected-within-subjects CIs 
% (only correction, no z; i.e. approx SE with only 2 conditions):
p{1} = boundedline(data.mu.time(timeIdx), ...
    nanmean(data.mat1(:, timeIdx)) - baseline, ...
    2 / (2-1) * nanstd(data.mat1(:, timeIdx) - ...
    data.matMean(:, timeIdx) + ... 
    repmat(data.matGrandMean(timeIdx), job.nValidSubs, 1)) ./ ...
    sqrt(job.nValidSubs), ...
    'cmap', job.twoLineColor(1, :), 'alpha', 'transparency', 0.2); %
set(p{1}, 'linewidth', lineWidth1); % adjust linewidth (not possible within boundedline)
p{2} = boundedline(data.mu.time(timeIdx), ...
    nanmean(data.mat2(:, timeIdx)) - baseline, ...
    2 / (2-1) * nanstd(data.mat2(:, timeIdx) - ...
    data.matMean(:, timeIdx) + ...
    repmat(data.matGrandMean(timeIdx), job.nValidSubs, 1)) ./ ...
    sqrt(job.nValidSubs), ...
    'cmap', job.twoLineColor(2, :), 'alpha', 'transparency', 0.2); %
set(p{2}, 'linewidth', lineWidth1); % adjust linewidth (not possible within boundedline)

% Extra vertical lines:
plot([0 0], get(gca,'ylim'), ':k', 'LineWidth', 3); % outcome onset
plot([.7 .7], get(gca,'ylim'), ':k', 'LineWidth', 3); % outcome offset

% x-axis label:
xlabel('Time (ms)', 'fontsize', fontSize, 'fontweight', 'bold', 'linewidth', lineWidth2);

% Retrieve y-axis limits, change later:
yLim    = get(gca, 'ylim'); % retrieve automatic settings, change later
yMinLim = yLim(1); 
yMaxLim = yLim(2); 

% Add x-axis settings and labels:
set(gca, 'xlim', [-0.25 0.8], ...
    'ylim', [yMinLim yMaxLim], ...
    'xtick', [-1.0 -0.75 -0.5 -0.25 0 0.25 0.5 0.70], ...
    'xtickLabel',{'-1000', '-750', '-500', '-250', 'FbOn', '250', '500', 'FbOff'}, ...
    'fontsize', fontSize, 'linewidth', lineWidth2);

% Add title y-label:
ylabel('Amplitude (A/cm²)', 'fontweight', 'bold', 'fontsize', fontSize); 

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS20D.png'));
pause(3)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20D_High.csv'));
csvwrite(fullFileName, data.mat1(:, timeIdx) - baseline);
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20D_Low.csv'));
csvwrite(fullFileName, data.mat2(:, timeIdx)- baseline);

% ----------------------------------------------------------------------- %
%% Figure S20E: Line plot V1:

clear job data

% ----------------------------------------------------------------------- %
%% Prepare contrast:

dirs.EEG                = fullfile(dirs.root, 'Log', 'EEG', 'OutcomeLockedResults');
dirs.timegroup          = fullfile(dirs.EEG, 'time_grouplevel');
job.dirs                = dirs;
job.nSub                = 36; % necessary for validSubs
job.sub2exclude         = [11 12 23 30]; % exclude 11 12 23 30 because very noisy during trial rejection
job.nFreqs              = 33; % 33
job.baselineSettings    = 'trend';
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'all'; % 'rel' or 'abs' or 'all'

% ----------------------------------------------------- %
% Channels of interest:
job.chanArea            = 'midfrontal'; % 'midfrontal'
% Contrast of interest:
job.contrastType        = 'Preferred'; 

% ----------------------------------------------------- %
% Use ROI (overwrites contrast):
 
job.ROI2use = {'GLM1V1Conj'};

if ~exist('data', 'var')
    job         = time_update_job(job);
    [job, data] = time_load_data(job);
    [job, data] = time_prepare_generic_data(job, data);
end

job             = time_update_job(job);
[~, data]       = time_prepare_contrast_data(job, data);

% ----------------------------------------------------------------------- %
%% Create plot:

% Fixed settings:
fontSize            = 32;
lineWidth1          = 5; % actual lines in plot
lineWidth2          = 3; % all other lines

% Determine baseline (i.e. lowest line at time point zero):
timeIdx             = dsearchn(data.mu.time', [-0.25 0.80]');
timeIdx             = timeIdx(1):timeIdx(end);
iBaseline           = find(round(data.mu.time, 3) == 0); % find baseline via time in sec
baseline            = min(nanmean(data.mat1(:, iBaseline)), nanmean(data.mat2(:, iBaseline)));

% Overwrite colors:
job.twoLineColor    = [1 0 0; 0 0 1];

% Start figure:
p   = cell(2, 1);
figure('Position', [100 100 1200 800], 'Color', 'white'); hold on % DCC screen

% Lines:
% Use Cousineau 2005 method for uncorrected-within-subjects CIs 
% (only correction, no z; i.e. approx SE with only 2 conditions):
p{1} = boundedline(data.mu.time(timeIdx), ...
    nanmean(data.mat1(:, timeIdx)) - baseline, ...
    2 / (2-1) * nanstd(data.mat1(:, timeIdx) - ...
    data.matMean(:, timeIdx) + ... 
    repmat(data.matGrandMean(timeIdx), job.nValidSubs, 1)) ./ ...
    sqrt(job.nValidSubs), ...
    'cmap', job.twoLineColor(1, :), 'alpha', 'transparency', 0.2); %
set(p{1}, 'linewidth', lineWidth1); % adjust linewidth (not possible within boundedline)
p{2} = boundedline(data.mu.time(timeIdx), ...
    nanmean(data.mat2(:, timeIdx)) - baseline, ...
    2 / (2-1) * nanstd(data.mat2(:, timeIdx) - ...
    data.matMean(:, timeIdx) + ...
    repmat(data.matGrandMean(timeIdx), job.nValidSubs, 1)) ./ ...
    sqrt(job.nValidSubs), ...
    'cmap', job.twoLineColor(2, :), 'alpha', 'transparency', 0.2); %
set(p{2}, 'linewidth', lineWidth1); % adjust linewidth (not possible within boundedline)

% Extra vertical lines:
plot([0 0], get(gca,'ylim'), ':k', 'LineWidth', 3); % outcome onset
plot([.7 .7], get(gca,'ylim'), ':k', 'LineWidth', 3); % outcome offset

% x-axis label:
xlabel('Time (ms)', 'fontsize', fontSize, 'fontweight', 'bold', 'linewidth', lineWidth2);

% Retrieve y-axis limits, change later:
yLim    = get(gca, 'ylim'); % retrieve automatic settings, change later
yMinLim = yLim(1); 
yMaxLim = yLim(2); 

% Add x-axis settings and labels:
set(gca, 'xlim', [-0.25 0.8], ...
    'ylim', [yMinLim yMaxLim], ...
    'xtick', [-1.0 -0.75 -0.5 -0.25 0 0.25 0.5 0.70], ...
    'xtickLabel',{'-1000', '-750', '-500', '-250', 'FbOn', '250', '500', 'FbOff'}, ...
    'fontsize', fontSize, 'linewidth', lineWidth2);

% Add title y-label:
ylabel('Amplitude (A/cm²)', 'fontweight', 'bold', 'fontsize', fontSize); 

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS20E.png'));
pause(3)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20E_High.csv'));
csvwrite(fullFileName, data.mat1(:, timeIdx) - baseline);
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20E_Low.csv'));
csvwrite(fullFileName, data.mat2(:, timeIdx)- baseline);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Load again saved regression coefficients for 7 fMRI ROIs:

EEGdomain           = 'time';

% 1) ROIs:
ROIs2use            = {'GLM1StriatumConj', 'GLM1ACCConjMan', ...
    'GLM1LeftMotorConj', 'GLM1vmPFCConjMan', 'GLM1PCCConj', ...
    'GLM1LeftITGConj', 'GLM1V1Conj'}; 
if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end

% 2) Behavioral regressors:
behav2use           = {'Updatestd', 'Updatedif'};

% 3) Perform only on selected trials:
selTrials           = 'all';

[job, dirs, betas]  = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
%% Figures S20F-L: Multiple topoplots:

addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64));

% General settings:
zlim            = 4;
nCols           = 4;

job.invalidSubs = [11 12 15 23 25 26 30]; % bad co-registrations & outliers for outcome-locked pre-processing & inspection (add 26)

% Timing settings:
startTime       = 0.0;
endTime         = 0.70;
steps           = 0.10;

% Create vectors:
startTimeVec    = startTime:steps:endTime;
endTimeVec      = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);
nPlot           = length(startTimeVec);
nRows           = ceil(nPlot/nCols); % number rows

% Select ROI:
for iROI = 1:length(ROIs2use)
        
    % Prepare data:
    [~, Tvalues]    = taft_postprocess_time_selectData(job, betas, iROI);

    % Start Plot:
    saveMat         = nan(64, nPlot);
    figure('Position', [100 100 1200 800], 'Color', 'white'); 

    for iPlot = 1:nPlot
        subplot(nRows, nCols, iPlot)

        cfg             = []; 
        cfg.figure      = gcf; 
        cfg.zlim        = [-1*zlim 1*zlim]; 
        cfg.marker      = 'on';
        cfg.style       = 'straight';
        cfg.layout      = 'easycapM11.mat'; 
        cfg.comment     = 'no';
        cfg.xlim        = [startTimeVec(iPlot) endTimeVec(iPlot)];
        cfg.colorbar    = 'no'; % want 'no', i.e. do it yourself % --> add externally
        cfg.style       = 'fill';
        ft_topoplotER(cfg, Tvalues);
        title(sprintf('%.1f-%.1f s', startTimeVec(iPlot), endTimeVec(iPlot)), 'fontsize', 28); % 1 digit: 28

        % Set linewidth:
        p = gca;
        for i = 7:11
            t = p.Children(i); 
            t.LineWidth = lineWidth;
        end

        % Save:
        timeIdx             = dsearchn(Tvalues.time', [startTimeVec(iPlot) endTimeVec(iPlot)]'); 
        timeIdx             = timeIdx(1):timeIdx(end);
        saveMat(:, iPlot)   = squeeze(mean(Tvalues.avg(:, timeIdx), 2));
    
    end

    % Save:
    if strcmp(ROIs2use{iROI}, 'GLM1StriatumConj') 
        figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS20F'));
        fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20F.csv'));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1ACCConjMan') 
        figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS20G'));
        fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20G.csv'));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1vmPFCConjMan') 
        figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS20H'));
        fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20H.csv'));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1PCCConj') 
        figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS20I'));
        fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20I.csv'));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1LeftMotorConj') 
        figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS20J'));
        fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20J.csv'));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1LeftITGConj') 
        figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS20K'));
        fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20K.csv'));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1V1Conj') 
        figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS20L'));
        fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS20L.csv'));
    end
    saveas(gcf, [figName '.png'])    
    saveMat         = horzcat(Tvalues.label, num2cell(saveMat));
    writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);

    % Close:
    pause(2)
    close gcf
    
end % end iROI
    
fprintf('Done :-)\n')

% END OF FILE.