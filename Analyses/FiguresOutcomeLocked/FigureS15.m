% FigureS15.m

% Plots for Supplementary Figure S15.
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

% Settings:
lineWidth1  = 3; % actual lines in plot
lineWidth2  = 1; % all other lines
fontSize1   = 24; % 24
fontSize2   = 12; 
transp      = 0.10;

% Overwrite colMat to be color-blind friendly:
job.colMat  = [0 113 116; 240 174 102; 87 196 173; 201 61 33; ...
    0 113 116; 240 174 102; 87 196 173; 201 61 33] ./ 255;

% Baselines:
iBaseline   = round(data.mu.time, 3) == 0; % find baseline via time in sec
baseline    = min(squeeze(nanmean(data.SubCondTime(job.validSubs, :, iBaseline)))); % still extract only valid subjects

% Time range:
startIdx    = find(round(data.mu.time, 3) == -0.25); % start for outcome-locked
stopIdx     = find(round(data.mu.time, 3) == 1); % end for outcome-locked
timeIdx     = startIdx:stopIdx;

% Initialize for saving:
saveMat     = nan(length(timeIdx), job.nCond);

% Start plot:
p = cell(job.nCond, 1); % delete p
figure('Position', [100 100 1600 800], 'Color', 'white'); hold on

for iCond = 1:job.nCond
    
    p{iCond} = boundedline(data.mu.time(timeIdx), ...
        squeeze(nanmean(data.SubCondTime(job.validSubs, iCond, timeIdx))) - baseline, ...
        job.nCond / (job.nCond-1) * squeeze(nanstd(squeeze(data.SubCondTime(job.validSubs, iCond, timeIdx)) - ...
        data.SubTime(job.validSubs, timeIdx)+repmat(data.GrandTime(timeIdx), length(job.validSubs), 1)))' ./ ...
        sqrt(length(job.validSubs)), ...
        'cmap', job.colMat(iCond, :), 'alpha', 'transparency', transp); % , 'linewidth', 2); hold on
    set(p{iCond}, 'linestyle', job.lineStyle{iCond});
    set(p{iCond}, 'Linewidth', lineWidth1);
    
    % Save:
    saveMat(:, iCond) = squeeze(nanmean(data.SubCondTime(job.validSubs, iCond, timeIdx))) - baseline;

end

% Axis labels:
xlabel('Time (ms)', 'fontsize', fontSize1);
ylabel('Amplitude (A/cm²)', 'fontsize', fontSize1); 

% Add x-axis labels, vertical lines:
set(gca, 'xlim', [-0.25 1.0], 'ylim', [-0.03 0.030], ...
    'xtick', [-1.0 -0.75 -0.5 -0.25 0 0.25 0.5 0.70 1 1.25], ...
    'xtickLabel', {'-1000', '-750', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1, 000'}, ...
    'fontsize', fontSize1, 'Linewidth', lineWidth2) %, 'ytick',-5:0.1:4)

% Extra vertical lines:
plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', lineWidth2); % outcome onset
plot([.7 .7], get(gca, 'ylim'), ':k', 'LineWidth', lineWidth2); % outcome offset

% Add legend:
legend([p{:}], job.condNames, 'fontSize', fontSize2); legend boxoff;

% Save:
saveas(gcf, fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS15.png'));
pause(3)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS15.csv'));
csvwrite(fullFileName, saveMat);

% END OF FILE.