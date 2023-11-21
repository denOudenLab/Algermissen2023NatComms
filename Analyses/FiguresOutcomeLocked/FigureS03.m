% FigureS03.m

% Plots for Supplementary Figure S03.
% The same approach is used for Figure 4 and Supplementary Figure S13.
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
dirs.TFgroup        = fullfile(dirs.EEG, 'TF_grouplevel');

% Add paths to helper files:
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel'));
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT'));

% Add Fieldtrip:
addpath /home/common/matlab/fieldtrip;
ft_defaults;

% ----------------------------------------------------------------------- %
%% Figure S03A: TF plot:

clear job data

% ----------------------------------------------------------------------- %
%% Prepare contrast:

job.dirs                = dirs;
job.nSub                = 36; % necessary for validSubs
job.sub2exclude         = [11 12 15 23 25 26 30]; % outliers in TAfT and fMRI
job.nFreqs              = 33; % 33
job.baselineSettings    = 'trend';
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'rel'; % 'rel' or 'abs' or 'all'

job.outERPcor           = false;

% ----------------------------------------------------- %
% Channels of interest:
job.chanArea            = 'midfrontal'; % 'midfrontal'
% Band settings:
job.band                = 'theta'; % 'theta' or 'alpha' or 'beta' or 'broad' % lowalpha middlebeta
% Contrast of interest:
job.contrastType        = 'Preferred'; 

if ~exist('data', 'var')
    job         = TF_update_job(job);
    [job, data] = TF_load_data(job);
    [job, data] = TF_prepare_generic_data(job, data);
end

job             = TF_update_job(job);
[~, data]       = TF_prepare_contrast_data(job, data);

% ----------------------------------------------------------------------- %
% Plot:

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)); nContour = 10; % default 10

% Plotting settings:
zlim        = 0.5;
fontSize    = 32;
yScale      = 'log';

% ----------------------------------------------------------------------- %
%% Create TF plot:

figure('Position', [100 100 1200 800], 'color', 'white'); hold on

% Contour plot: 
contourf(data.mu.time, data.mu.freq, data.TF2plot, ...
    nContour, 'linestyle', 'none'); hold on

% Vertical lines:
plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', 3);
plot([0.7 0.7], get(gca, 'ylim'), ':k', 'LineWidth', 3);

% Settings:
set(gca, 'xlim', [-0.25 0.8], 'ylim', [1 33], 'clim', [-1*zlim 1*zlim], ...
    'xtick', [-1 -0.50 -0.25 0 0.25 0.5 0.7 1], ...
    'xtickLabel', {'-1, 000', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
    'ytick', [2 4 8 16 32], 'yscale', yScale, ...
    'fontsize', fontSize, 'Linewidth', 3); % -.25 1.3

% Labels:
ylabel('Frequency (Hz)', 'fontsize', fontSize, 'fontweight', 'bold');
xlabel('Time (ms)', 'fontsize', fontSize, 'fontweight', 'bold');

% Color bar:
colorbar('Ticks', (-1*zlim):(zlim/2):zlim, 'Fontsize', fontSize);

% Save:
figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03A');
saveas(gcf, [figName '.png']);

% Close:
pause(1);
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03A.csv'));
csvwrite(fullFileName, data.TF2plot);

% ----------------------------------------------------------------------- %
%% Figure S03B: Line plot theta power:

clear job data

% ----------------------------------------------------------------------- %
%% Prepare contrast:

job.dirs                = dirs;
job.nSub                = 36; % necessary for validSubs
job.sub2exclude         = [11 12 15 23 25 26 30]; % outliers in TAfT and fMRI
job.nFreqs              = 33; % 33
job.baselineSettings    = 'trend';
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'none'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'all'; % 'rel' or 'abs' or 'all'

job.outERPcor           = false;

% ----------------------------------------------------- %
% Channels of interest:
job.chanArea            = 'midfrontal'; % 'midfrontal'
% Band settings:
job.band                = 'theta'; % 'theta' or 'alpha' or 'beta' or 'broad' % lowalpha middlebeta
% Contrast of interest:
job.contrastType        = 'Preferred'; 

if ~exist('data', 'var')
    job         = TF_update_job(job);
    [job, data] = TF_load_data(job);
    [job, data] = TF_prepare_generic_data(job, data);
end

job             = TF_update_job(job);
[job, data]     = TF_prepare_contrast_data(job, data);

% ----------------------------------------------------------------------- %
%% Create lineplot:

% New color map for red-green blind:
job.colMat  = [0 113 116; 87 196 173; 240 174 102; 201 61 33] ./ 255;
job.colMat  = repmat(job.colMat, job.nCond / size(job.colMat, 1), 1);

% Baselines:
iBaseline   = find(round(data.mu.time, 3) == 0); % find baseline via time in sec
baseline    = min(squeeze(nanmean(data.SubCondTime(job.validSubs, :, iBaseline)))); % still extract only valid subjects
until       = find(round(data.mu.time, 3) == 1); % end for outcome-locked

% Other settings:
xLim        = [-0.25 0.8];
lineWidth   = 5;
fontSize    = 32;
transp      = 0.10;

% Start plot:
figure('Position', [100 100 1200 800], 'color', 'white'); hold on

% Loop over conditions, make bounded line plots:
p   = cell(job.nCond, 1);
for iCond = 1:job.nCond

    p{iCond} = boundedline(data.mu.time(1:until), ...
        squeeze(nanmean(data.SubCondTime(job.validSubs, iCond, 1:until))) - baseline, ...
        job.nCond / (job.nCond-1) * ...
        squeeze(nanstd(squeeze(data.SubCondTime(job.validSubs, iCond, 1:until)) - ...
        data.SubTime(job.validSubs, 1:until) + ...
        repmat(data.GrandTime(1:until), length(job.validSubs), 1)))' ./ ...
        sqrt(length(job.validSubs)), ...
        'cmap', job.colMat(iCond, :) , 'alpha', 'transparency', transp); % , 'linewidth', 2); hold on
    set(p{iCond}, 'linestyle', job.lineStyle{iCond})
    set(p{iCond}, 'Linewidth', lineWidth)

end

% Axis labels and limits:
xlabel('Time (ms)', 'fontweight', 'bold', 'fontsize', fontSize)
yLim = get(gca, 'ylim'); % yLim = [0 2.5];
yMinLim = yLim(1); 
yMaxLim = yLim(2); 
    
if strcmp(job.contrastType, 'Preferred') && strcmp(job.band, 'theta')
    yMinLim = -0.2; % adjust a bit: use -0.2
    yMaxLim = 1.4; % adjust a bit: use 1.3

elseif strcmp(job.contrastType, 'Preferred') && strcmp(job.band, 'beta')
    yMinLim = -0.3; % adjust a bit: use -0.2
    yMaxLim = 0.3; % adjust a bit: use 1.3
end

set(gca, 'xlim', xLim, 'ylim', [yMinLim yMaxLim], ...
    'xtick', [-1 -0.75 -0.50 -0.25 0 0.25 0.5 0.7 1 1.25], ...
    'xtickLabel', {'-1, 000', '-750', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1, 000'}, ...
    'fontsize', fontSize, 'Linewidth', 3);
set(gca, 'TickLength', [0 0]);

% Gray patch:
yLim = get(gca, 'ylim');
patch([0 .7 .7 0], [yLim(1) yLim(1) yLim(2) yLim(2)], [0.8 0.8 0.8], 'facealpha', 0.2, 'EdgeColor', 'none');

% Save:
figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03B');
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data file:
% See job.condNames
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03B_Reward.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 1, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03B_NoReward.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 2, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03B_NoPunishment.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 3, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03B_Punishment.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 4, 1:until));

% ----------------------------------------------------------------------- %
%% Create topoplot:

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)) ;

% Other settings:
lineWidth   = 3;

% Sigtime:
sigTime     = job.sigTime;

% Further settings:
if isfield(job, 'freq')
    zlim = 0.5;
else
    zlim = 0.010;
end

% Start figure:
figure('Position', [100 100 1200 800]); hold on

% Set config object for Fieldtrip:
cfg                     = []; 
cfg.figure              = gcf; 
cfg.zlim                = [-1*zlim 1*zlim]; 
cfg.marker              = 'on'; 
cfg.layout              = 'easycapM11.mat'; 
cfg.comment             = 'no'; cfg.xlim = sigTime;
cfg.style               = 'fill';
cfg.highlight           = 'on'; 
cfg.highlightchannel    = job.channels; 
cfg.highlightsymbol     = 'o'; 
cfg.highlightsize       = 12; 
cfg.highlightfontsize   = 12;
cfg.colorbar            = 'yes';

% Make plot:
if isfield(job, 'freq')
    cfg.ylim = job.freq; 
    ft_topoplotTFR(cfg, data.topo2plot);
else
    ft_topoplotER(cfg, data.topo2plot);
end

% Set linewidth around it:
p = gca;
for i = 12:16
    t = p.Children(i); 
    t.LineWidth = lineWidth;
end

% Save:
figName = fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS03B_topo');
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data file:
fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03B_topo.csv'));
timeIdx         = dsearchn(data.topo2plot.time', sigTime'); timeIdx = timeIdx(1):timeIdx(end);
freqIdx         = dsearchn(data.topo2plot.freq', job.freq'); freqIdx = freqIdx(1):freqIdx(end);
saveMat         = squeeze(mean(mean(data.topo2plot.powspctrm(:, freqIdx, timeIdx), 3), 2));
saveMat         = horzcat(data.topo2plot.label, num2cell(saveMat));
writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);
    
% ----------------------------------------------------------------------- %
%% Figure S03C: Line plot beta power:

clear job data

% ----------------------------------------------------------------------- %
%% Prepare contrast:

job.dirs                = dirs;
job.nSub                = 36; % necessary for validSubs
job.sub2exclude         = [11 12 15 23 25 26 30]; % outliers in TAfT and fMRI
job.nFreqs              = 33; % 33
job.baselineSettings    = 'trend';
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'none'; % 'Go' or 'Hand' or 'none'  
job.outcomeSettings     = 'all'; % 'rel' or 'abs' or 'all'

job.outERPcor           = false;

% ----------------------------------------------------- %
% Channels of interest:
job.chanArea            = 'midfrontal'; % 'midfrontal'
% Band settings:
job.band                = 'beta'; % 'theta' or 'alpha' or 'beta' or 'broad' % lowalpha middlebeta
% Contrast of interest:
job.contrastType        = 'Preferred'; 

if ~exist('data', 'var')
    job         = TF_update_job(job);
    [job, data] = TF_load_data(job);
    [job, data] = TF_prepare_generic_data(job, data);
end

job             = TF_update_job(job);
[job, data]     = TF_prepare_contrast_data(job, data);

% ----------------------------------------------------------------------- %
%% Create lineplot:

% New color map for red-green blind:
job.colMat  = [0 113 116; 87 196 173; 240 174 102; 201 61 33] ./ 255;
job.colMat  = repmat(job.colMat, job.nCond / size(job.colMat, 1), 1);

% Baselines:
iBaseline   = find(round(data.mu.time, 3) == 0); % find baseline via time in sec
baseline    = min(squeeze(nanmean(data.SubCondTime(job.validSubs, :, iBaseline)))); % still extract only valid subjects
until       = find(round(data.mu.time, 3) == 1); % end for outcome-locked

% Other settings:
xLim        = [-0.25 0.8];
lineWidth   = 5;
fontSize    = 32;
transp      = 0.10;

% Start plot:
figure('Position', [100 100 1200 800], 'color', 'white'); hold on

% Loop over conditions, make bounded line plots:
p   = cell(job.nCond, 1);
for iCond = 1:job.nCond

    p{iCond} = boundedline(data.mu.time(1:until), ...
        squeeze(nanmean(data.SubCondTime(job.validSubs, iCond, 1:until))) - baseline, ...
        job.nCond / (job.nCond-1) * ...
        squeeze(nanstd(squeeze(data.SubCondTime(job.validSubs, iCond, 1:until)) - ...
        data.SubTime(job.validSubs, 1:until) + ...
        repmat(data.GrandTime(1:until), length(job.validSubs), 1)))' ./ ...
        sqrt(length(job.validSubs)), ...
        'cmap', job.colMat(iCond, :) , 'alpha', 'transparency', transp); % , 'linewidth', 2); hold on
    set(p{iCond}, 'linestyle', job.lineStyle{iCond})
    set(p{iCond}, 'Linewidth', lineWidth)

end

% Axis labels and limits:
xlabel('Time (ms)', 'fontweight', 'bold', 'fontsize', fontSize)
yLim = get(gca, 'ylim');
yMinLim = yLim(1); 
yMaxLim = yLim(2); 


if strcmp(job.contrastType, 'Preferred') && strcmp(job.band, 'theta')
    yMinLim = -0.2; % adjust a bit: use -0.2
    yMaxLim = 1.4; % adjust a bit: use 1.3

elseif strcmp(job.contrastType, 'Preferred') && strcmp(job.band, 'beta')
    yMinLim = -0.3; % adjust a bit: use -0.2
    yMaxLim = 0.3; % adjust a bit: use 1.3
end

set(gca, 'xlim', xLim, 'ylim', [yMinLim yMaxLim], ...
    'xtick', [-1 -0.75 -0.50 -0.25 0 0.25 0.5 0.7 1 1.25], ...
    'xtickLabel', {'-1, 000', '-750', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1, 000'}, ...
    'fontsize', fontSize, 'Linewidth', 3);
set(gca, 'TickLength', [0 0]);

% Gray patch:
yLim    = get(gca, 'ylim');
patch([0 .7 .7 0], [yLim(1) yLim(1) yLim(2) yLim(2)], [0.8 0.8 0.8], 'facealpha', 0.2, 'EdgeColor', 'none');

% Save:
figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03C');
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data file:
% See job.condNames
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03C_Reward.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 1, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03C_NoReward.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 2, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03C_NoPunishment.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 3, 1:until));
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03C_Punishment.csv'));
csvwrite(fullFileName, data.SubCondTime(job.validSubs, 4, 1:until));

% ----------------------------------------------------------------------- %
%% Create topoplot:

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)) ;

% Other settings:
lineWidth   = 3;

% Sigtime:
sigTime     = job.sigTime;

% Further settings:
if isfield(job, 'freq')
    zlim = 0.5;
else
    zlim = 0.010;
end

% Start figure:
figure('Position', [100 100 1200 800]); hold on

% Set config object for Fieldtrip:
cfg                     = []; 
cfg.figure              = gcf; 
cfg.zlim                = [-1*zlim 1*zlim]; 
cfg.marker              = 'on'; 
cfg.layout              = 'easycapM11.mat'; 
cfg.comment             = 'no'; cfg.xlim = sigTime;
cfg.style               = 'fill';
cfg.highlight           = 'on'; 
cfg.highlightchannel    = job.channels; 
cfg.highlightsymbol     = 'o'; 
cfg.highlightsize       = 12; 
cfg.highlightfontsize   = 12;
cfg.colorbar            = 'yes';

% Make plot:
if isfield(job, 'freq')
    cfg.ylim = job.freq; 
    ft_topoplotTFR(cfg, data.topo2plot);
else
    ft_topoplotER(cfg, data.topo2plot);
end

% Set linewidth around it:
p = gca;
for i = 12:16
    t = p.Children(i); 
    t.LineWidth = lineWidth;
end

% Save:
figName = fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS03C_topo');
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data file:
fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03C_topo.csv'));
timeIdx         = dsearchn(data.topo2plot.time', sigTime'); timeIdx = timeIdx(1):timeIdx(end);
freqIdx         = dsearchn(data.topo2plot.freq', job.freq'); freqIdx = freqIdx(1):freqIdx(end);
saveMat         = squeeze(mean(mean(data.topo2plot.powspctrm(:, freqIdx, timeIdx), 3), 2));
saveMat         = horzcat(data.topo2plot.label, num2cell(saveMat));
writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Figure S03D: EEG correlates PE_BIAS:

EEGdomain           = 'TF';

% 1) ROIs:
ROIs2use            = {''}; % empty
if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end

% 2) Behavioral regressors:
behav2use           = {'Updatebias', 'fbrel'}; % Updates + outcome sign

% 3) Perform only on selected trials:
selTrials           = 'all';

[job, dirs, betas]  = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
%% Create TF plot:

% NOTE: Within within taft_postprocess_TF_TFplot, set latency to
% cfg.latency = [-0.25 0.8]; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64));

% Exclude subjects:
job.invalidSubs = [11 12 15 23 25 26 30]; % bad co-registrations & outliers for outcome-locked pre-processing & inspection (add 26)

% Select channel:
iROI            = 1;

% Select ROI:
selChans        = {'Fz', 'FCz', 'Cz'};
        
% Set seed:
rng(20190822) % set random number generator for constant p-values

% Preprocess betas:
[sortBetas,~]   = taft_postprocess_TF_selectData(job, betas, iROI);
sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:

% Plot to retrieve t-values and p-values:
[tg, corrp]     = taft_postprocess_TF_TFplot(job, dirs, sortBetas, iROI, selChans, 2, 1000, 2, false); % iROI, selChans, thresh, nP
close gcf

% Settings:
timeIdx         = size(tg, 4); % maximum time
zlim            = 3;
lineWidth       = 3;
fontSize        = 32; % 44; % uni: 32; smaller: 24 (3/4)
pCrit           = 0.05;

% Exception in pCrit because cluster not-significant when tested
% -0.25 - 0.8 sec. across all frequencies:
if strcmp(job.regNames(iROI), 'Updatestd') && isfield(job, 'suffix'); pCrit = .10; end
if strcmp(job.regNames(iROI), 'Updatebias') && ~isfield(job, 'suffix'); pCrit = .17; end

% ---------------------------------------- %
% Start figure:

figure('Position', [100 100 1200 800], 'color', 'white'); hold on % same as TFplot for standard EEG

% ---------------------------------------- %
% Make contour plot:
contourf(sortBetas{1}.time(1:timeIdx), sortBetas{1}.freq, squeeze(real(tg)), ...
    nContour, 'linestyle', 'none'); % 1.3 sec at index 53

% ---------------------------------------- %
% Add vertical lines:

plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', 3);
plot([0.7 0.7], get(gca, 'ylim'), ':k', 'LineWidth', 3);

% ---------------------------------------- %
% Settings:

set(gca, 'xlim', xLim, 'ylim', [1.5 33], 'clim', [-1*zlim  1*zlim], ...
    'xtick', [-1 -0.5 -0.25 0 0.25 0.5 0.7 1], 'xtickLabel', {'-1', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
    'ytick', [2 4 8 16 32], 'yscale', 'log', ... % 'ytick', [2 4 8 16 32], ... 
    'fontsize', fontSize, 'Linewidth', lineWidth); % -.25 1.3

% Labels:
xlabel('Time (ms)', 'FontSize', fontSize, 'FontName', 'Arial', 'fontweight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', fontSize, 'FontName', 'Arial', 'fontweight', 'bold');

% Color bar:
colorbar('Ticks', (-1*zlim):(zlim/2):zlim, 'Fontsize', fontSize);

% ---------------------------------------- %
% Add highlights based on p-values:

if sum(corrp(:) < pCrit) > 0 % if anything significant
   isSig = squeeze(double(corrp < pCrit)); % map significant clusters
    if strcmp(job.regNames{iROI}, 'Updatebias')
             isSig(:, 1:18)     = 0; % don't display clusters not significant in longer time window
             isSig(8:end, :)    = 0; % don't display clusters not significant in longer time window
             isSig(:, 34:end)   = 0; % don't display clusters not significant in longer time window
    end
    % see https://nl.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a#answer_211204

    pause(1) % pause to allow for plotting both plots into one:           
    hold on % on top of old plot

    % Extra contour:
    [~, hContour]  = contourf(sortBetas{1}.time(1:timeIdx), sortBetas{1}.freq, squeeze(isSig), 1);
    hContour.LineWidth = 5;
    drawnow;  % this is important, to ensure that FacePrims is ready in the next line!

    hFills = hContour.FacePrims;  % array of TriangleStrip objects
    [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
    for idx = 1:numel(hFills)
        hFills(idx).ColorData(4) = 1;   % default=255
        pause(1)
    end

    hold off

end

% -------------------------------- %
% Save figure:
fig = gcf;
figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03D');
saveas(gcf, [figName '.png'])

% Close:
pause(3);
close gcf

% Save source data file:
if strcmp(job.regNames(iROI), 'Updatebias')
    fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03D.csv');
    csvwrite(fullFileName, squeeze(real(tg)));
end
fprintf('Done :-)\n')

% ----------------------------------------------------------------------- %
%% Figure S03D: Topoplot:

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)) ; % default 10

% Color range:
zlim        = 3;

% Line width:
lineWidth   = 3;

% ------------------------------------------------- %
% Select ROI:

iROI = 1; iTime = [0.225 0.475]; iFreq = [1 5]; 

% ------------------------------------------------- %
% Compute T-values across subjects:

[sortBetas, Tvalues] = taft_postprocess_TF_selectData(job, betas, iROI);

% Plot:
figure('units', 'normalized', 'outerposition', [0 0 1 1]); hold on % fullscreen--no risk of implicitly resizing any axes
cfg         = []; 
cfg.figure  = gcf; 
cfg.ylim    = iFreq; 
cfg.zlim    = [-1*zlim 1*zlim]; 
cfg.marker  ='on'; 
cfg.layout  = 'easycapM11.mat'; 
cfg.comment = 'no'; 
cfg.xlim    = [iTime(1) iTime(end)];
cfg.style   = 'fill';
ft_topoplotTFR(cfg, Tvalues);
title(sprintf('%s, %.2f to %.2f sec', job.regNames{iROI}, iTime(1), iTime(end)), 'fontsize', 32);

% Set linewidth:
p = gca;
for i = 7:11
    t = p.Children(i); 
    t.LineWidth = lineWidth;
end

% Save:
figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03D_topo.png');
saveas(gcf, figName);

% Close:
pause(2);
close gcf

% Save source data file:
fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS03D_topo.csv'));
timeIdx         = dsearchn(Tvalues.time', iTime'); timeIdx = timeIdx(1):timeIdx(end);
freqIdx         = dsearchn(Tvalues.freq', iFreq'); freqIdx = freqIdx(1):freqIdx(end);
saveMat         = squeeze(mean(mean(Tvalues.powspctrm(:, freqIdx, timeIdx), 3), 2));
saveMat         = horzcat(Tvalues.label, num2cell(saveMat));
writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);

% ----------------------------------------------------------------------- %
%% Figure S03E & S03F: EEG correlates PE_STD and PE_DIF:

EEGdomain           = 'TF';

% 1) ROIs:
ROIs2use            = {''}; % empty
if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end

% 2) Behavioral regressors:
behav2use           = {'Updatestd', 'Updatedif', 'fbrel'}; % Updates + outcome sign

% 3) Perform only on selected trials:
selTrials           = 'all';

[job, dirs, betas]  = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
%% Plot:

% NOTE: Within within taft_postprocess_TF_TFplot, set latency to
% cfg.latency = [-0.25 0.8]; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)); cmap = 'redblue'; nContour = 10; % default 10


% Exclude subjects:
job.invalidSubs         = [11 12 15 23 25 26 30]; % bad co-registrations & outliers for outcome-locked pre-processing & inspection (add 26)

% Selected channels:
selChans                = {'Fz', 'FCz', 'Cz'};

for iROI = [1 2] % iROI = 9; % UpdateStd and UpdateDif: Fig. 4
        
    % Set seed:
    rng(20190822) % set random number generator for constant p-values

    % Preprocess betas:
    [sortBetas,~]   = taft_postprocess_TF_selectData(job, betas, iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:

    % Plot to retrieve t-values and p-values:
    [tg, corrp]     = taft_postprocess_TF_TFplot(job, dirs, sortBetas, iROI, selChans, 2, 1000, 2, false); % iROI, selChans, thresh, nP
    close gcf

    % Settings:
    timeIdx         = size(tg, 4); % maximum time
    xLim            = [-0.25 0.8];
    zlim            = 3;
    lineWidth       = 3;
    fontSize        = 32; % 44; % uni: 32; smaller: 24 (3/4)
    pCrit           = 0.05;

    % Exception in pCrit because cluster not-significant when tested
    % -0.25 - 0.8 sec. across all frequencies:
    if contains(job.regNames(iROI), 'PCCConj'); pCrit = 0.20; end
    if strcmp(job.regNames(iROI), 'Updatestd') && isfield(job, 'suffix'); pCrit = .10; end
    if strcmp(job.regNames(iROI), 'Updatebias') && ~isfield(job, 'suffix'); pCrit = .17; end

    % ------------------------------------------------------------------- %
    % Start figure:

    figure('Position', [100 100 1200 800], 'color', 'white'); hold on % same as TFplot for standard EEG
    
    % ------------------------------------------------------------------- %
    % Make contour plot:
    contourf(sortBetas{1}.time(1:timeIdx), sortBetas{1}.freq, squeeze(real(tg)), ...
        nContour, 'linestyle', 'none'); % 1.3 sec at index 53
    
    % ------------------------------------------------------------------- %
    % Add vertical lines:

    plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', 3);
    plot([0.7 0.7], get(gca, 'ylim'), ':k', 'LineWidth', 3);

    % ------------------------------------------------------------------- %
    % Settings:
    
    set(gca, 'xlim', xLim, 'ylim', [1.5 33], 'clim', [-1*zlim  1*zlim], ...
        'xtick', [-1 -0.5 -0.25 0 0.25 0.5 0.7 1], 'xtickLabel', {'-1', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
        'ytick', [2 4 8 16 32], 'yscale', 'log', ... % 'ytick', [2 4 8 16 32], ... 
        'fontsize', fontSize, 'Linewidth', lineWidth) % -.25 1.3

    % Labels:
    xlabel('Time (ms)', 'FontSize', fontSize, 'FontName', 'Arial', 'fontweight', 'bold');
    ylabel('Frequency (Hz)', 'FontSize', fontSize, 'FontName', 'Arial', 'fontweight', 'bold');

    % Color bar:
    colorbar('Ticks', (-1*zlim):(zlim/2):zlim, 'Fontsize', fontSize)

    % ------------------------------------------------------------------- %
    % Add highlights based on p-values:

    if sum(corrp(:) < pCrit) > 0 % if anything significant
       isSig = squeeze(double(corrp < pCrit)); % map significant clusters
        % see https://nl.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a#answer_211204

        pause(1) % pause to allow for plotting both plots into one:           
        hold on % on top of old plot

        % Extra contour:
        [~, hContour]  = contourf(sortBetas{1}.time(1:timeIdx), sortBetas{1}.freq, squeeze(isSig), 1);
        hContour.LineWidth = 5;
        drawnow;  % this is important, to ensure that FacePrims is ready in the next line!

        hFills = hContour.FacePrims;  % array of TriangleStrip objects
        [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
        for idx = 1:numel(hFills)
            hFills(idx).ColorData(4) = 1;   % default=255
            pause(1)
        end

        hold off

    end

    % ------------------------------------------------------------------- %
    % Save and close:

    % Save source data file:
    if strcmp(job.regNames(iROI), 'Updatestd')
        figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03E');
        saveas(gcf, [figName '.png']);
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03E.csv');
        csvwrite(fullFileName, squeeze(real(tg)));
    end
    % Save source data file:
    if strcmp(job.regNames(iROI), 'Updatedif')
        figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03F');
        saveas(gcf, [figName '.png']);
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS03F.csv');
        csvwrite(fullFileName, squeeze(real(tg)));
    end

    % Close:
    pause(3);
    close gcf

end % end iROI
fprintf('Done :-)\n')

% END OF FILE.