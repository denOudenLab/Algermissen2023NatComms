% FigureS14.m

% Plots for Supplementary Figure S14.
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
dirs.EEG            = fullfile(dirs.root, 'Log', 'EEG', 'OutcomeLockedResults');
dirs.timegroup      = fullfile(dirs.EEG, 'time_grouplevel');

% Add paths to helper files:
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel'));

% Add Fieldtrip:
addpath /home/common/matlab/fieldtrip;
ft_defaults;

% ----------------------------------------------------------------------- %
%% Figure S14A:

clear job data

% ----------------------------------------------------------------------- %
%% Prepare contrast:

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
    'xtick', [-0.25 0 0.25 0.5 0.70], ...
    'xtickLabel',{'-250', 'FbOn', '250', '500', 'FbOff'}, ...
    'fontsize', fontSize, 'linewidth', lineWidth2);

% Add title y-label:
ylabel('Amplitude (A/cm²)', 'fontweight', 'bold', 'fontsize', fontSize); 

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS14A.png'));
pause(3)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14A_Favorable.csv'));
csvwrite(fullFileName, data.mat1(:, timeIdx) - baseline);
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14A_Unfavorable.csv'));
csvwrite(fullFileName, data.mat2(:, timeIdx)- baseline);

% ----------------------------------------------------------------------- %
%% Figure S14B:

% same data and before, only update job

% ----------------------------------------------------------------------- %
%% Update channels:

% Channels of interest:
job.chanArea    = 'rightoccipital';

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
    'xtick', [-0.25 0 0.25 0.5 0.70], ...
    'xtickLabel',{'-250', 'FbOn', '250', '500', 'FbOff'}, ...
    'fontsize', fontSize, 'linewidth', lineWidth2);

% Add title y-label:
ylabel('Amplitude (A/cm²)', 'fontweight', 'bold', 'fontsize', fontSize); 

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS14B.png'));
pause(3)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14B_Favorable.csv'));
csvwrite(fullFileName, data.mat1(:, timeIdx) - baseline);
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14B_Unfavorable.csv'));
csvwrite(fullFileName, data.mat2(:, timeIdx)- baseline);

% ----------------------------------------------------------------------- %
%% Figure S14C: Multiple topoplots

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64));

% Settings:
zlim            = 0.010;
nCols           = 4;
lineWidth       = 3;

% Timing settings:
startTime       = 0.0;
endTime         = 0.7;
steps           = 0.1;
startTimeVec    = startTime:steps:endTime;
endTimeVec      = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);
nPlot           = length(endTimeVec);
nRows           = ceil(nPlot/5);
saveMat         = nan(64, nPlot);

% Start plot:
close gcf
figure('Position', [100 100 1200 800], 'Color', 'white'); 
for iPlot = 1:nPlot
    subplot(nRows, nCols, iPlot)
    cfg             = []; 
    cfg.figure      = gcf; 
    cfg.zlim        = [-1*zlim 1*zlim]; 
    cfg.marker      ='on'; 
    cfg.style       ='straight';
    cfg.layout      = 'easycapM11.mat'; 
    cfg.comment     = 'no'; 
    cfg.xlim        = [startTimeVec(iPlot) endTimeVec(iPlot)];
    cfg.colorbar    = 'no'; % want 'no', i.e. do it yourself % --> add externally
    cfg.style       = 'fill';
    ft_topoplotER(cfg, data.topo2plot);
    title(sprintf('%.1f-%.1f s', startTimeVec(iPlot), endTimeVec(iPlot)), 'fontsize', 28);
    
    % Set linewidth:
    p = gca;
    for i = 7:11
        t = p.Children(i); 
        t.LineWidth = lineWidth;
    end
    
    % Save:
    timeIdx             = dsearchn(data.topo2plot.time', [startTimeVec(iPlot) endTimeVec(iPlot)]'); 
    timeIdx             = timeIdx(1):timeIdx(end);
    saveMat(:, iPlot)   = squeeze(mean(data.topo2plot.avg(:, timeIdx), 2));

end

% Save:
figName = fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS14C');
saveas(gcf, [figName '.png']);

% Close:
pause(2)
close gcf
fprintf('Done :-)\n');

% Save source data file:
saveMat         = horzcat(data.topo2plot.label, num2cell(saveMat));
fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14C.csv'));
writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);

% ----------------------------------------------------------------------- %
%% Figure S14D:

% same data and before, only update job

% ----------------------------------------------------------------------- %
%% Update contrast and channels:

% Contrast of interest:
job.contrastType        = 'Action'; 

% Channels of interest:
job.chanArea    = 'midfrontal';

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
    'xtick', [-0.25 0 0.25 0.5 0.70], ...
    'xtickLabel',{'-250', 'FbOn', '250', '500', 'FbOff'}, ...
    'fontsize', fontSize, 'linewidth', lineWidth2);

% Add title y-label:
ylabel('Amplitude (A/cm²)', 'fontweight', 'bold', 'fontsize', fontSize); 

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS14D.png'));
pause(3)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14D_Go.csv'));
csvwrite(fullFileName, data.mat1(:, timeIdx) - baseline);
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14D_NoGo.csv'));
csvwrite(fullFileName, data.mat2(:, timeIdx)- baseline);

% ----------------------------------------------------------------------- %
%% Figure S14E:

% same data and before, only update job

% ----------------------------------------------------------------------- %
%% Update channels:

% Channels of interest:
job.chanArea    = 'rightoccipital';

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
    'xtick', [-0.25 0 0.25 0.5 0.70], ...
    'xtickLabel',{'-250', 'FbOn', '250', '500', 'FbOff'}, ...
    'fontsize', fontSize, 'linewidth', lineWidth2);

% Add title y-label:
ylabel('Amplitude (A/cm²)', 'fontweight', 'bold', 'fontsize', fontSize); 

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS14E.png'));
pause(3)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14E_Go.csv'));
csvwrite(fullFileName, data.mat1(:, timeIdx) - baseline);
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14E_NoGo.csv'));
csvwrite(fullFileName, data.mat2(:, timeIdx)- baseline);

% ----------------------------------------------------------------------- %
%% Figure S14F: Multiple topoplots

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)) ;

% Settings:
zlim            = 0.010;
nCols           = 4;
lineWidth       = 3;

% Timing settings:
startTime       = 0.0;
endTime         = 0.7;
steps           = 0.1;
startTimeVec    = startTime:steps:endTime;
endTimeVec      = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);
nPlot           = length(endTimeVec);
nRows           = ceil(nPlot/5);
saveMat         = nan(64, nPlot);

% Start plot:
close gcf
figure('Position', [100 100 1200 800], 'Color', 'white'); 
for iPlot = 1:nPlot
    subplot(nRows, nCols, iPlot)
    cfg             = []; 
    cfg.figure      = gcf; 
    cfg.zlim        = [-1*zlim 1*zlim]; 
    cfg.marker      ='on'; 
    cfg.style       ='straight';
    cfg.layout      = 'easycapM11.mat'; 
    cfg.comment     = 'no'; 
    cfg.xlim        = [startTimeVec(iPlot) endTimeVec(iPlot)];
    cfg.colorbar    = 'no'; % want 'no', i.e. do it yourself % --> add externally
    cfg.style       = 'fill';
    ft_topoplotER(cfg, data.topo2plot);
    title(sprintf('%.1f-%.1f s', startTimeVec(iPlot), endTimeVec(iPlot)), 'fontsize', 28);
    
    % Set linewidth:
    p = gca;
    for i = 7:11
        t = p.Children(i); 
        t.LineWidth = lineWidth;
    end
    
    % Save:
    timeIdx             = dsearchn(data.topo2plot.time', [startTimeVec(iPlot) endTimeVec(iPlot)]'); 
    timeIdx             = timeIdx(1):timeIdx(end);
    saveMat(:, iPlot)   = squeeze(mean(data.topo2plot.avg(:, timeIdx), 2));

end

% Save:
figName = fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS14F');
saveas(gcf, [figName '.png']);

% Close:
pause(2)
close gcf
fprintf('Done :-)\n');

% Save source data file:
saveMat         = horzcat(data.topo2plot.label, num2cell(saveMat));
fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS14F.csv'));
writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);

% END OF FILE.