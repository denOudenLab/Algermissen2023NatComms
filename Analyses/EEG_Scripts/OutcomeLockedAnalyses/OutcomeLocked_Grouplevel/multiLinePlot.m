function multiLinePlot(job, data, errorbar)

% multiLinePlot(job, data, errorbar)
%
% Creates lineplot for all conditions specified in data.
% Works for both ERP and TFR data.
%
% INPUTS:
% job                   = cell with relevant fields for plotting settings, 
% required fields:
%   .nCond              = numeric, number of conditions in data.
%   .validSubs          = numeric vector, subjects to include.
%   .nValSubs           = numeric, number of subjects to include.
%   .condNames          = vector of strings, labels for conditions.
%   .colMat             = matrix, colors of lines to plot per condition.
%   .freq               = numeric vector of 2 elements, range of selected
%   frequencies (required if input data is TFR).
%   .channels           = vector of strings, selected channel names.
%   .accSettings        = string, type of accuracy setting for which 
%   conditions are split, 'correct', 'incorrect', 'bothAcc' (both
%   accuracies separately), 'noAcc' (no distinction), 'reqAct' (required
%   action) (for determining optimal y-axis limits).
%   .contrastType       = string, contrast to be used, either 
%   'Congruency', 'Go', 'GoLeftRight', 'GoValence', 'GoAccuracy', 
%   'Accuracy', 'CongruentAccuracy' 'IncongruentAccuracy' (for determining 
%   optimal y-axis limits).
%   .sigTime            = numeric vector of 2 elements, timing for
%   significant difference to highlight (optional).
% data                  = cell with the following fields:
%   .mu                 = grand average across both conditions and
%   subjects.
%   .SubCondTime        = per subject per condition, over time, averaged 
%   over frequencies and channels.
%   .SubTime            = SubCondTime averaged over conditions.
%   .GrandTime          = SubTime averaged over subjects.
% errorbar              = add errorbars (true) or not (false) to plot
% (default: true):
%
% OUTPUTS:
% none, just plotting.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

% ----------------------------------------------------------------------- %
%% Complete settings:

if nargin < 3
    errorbar = true;
end

% ----------------------------------------------------------------------- %
%% Close any open figure windows:

close all

% ----------------------------------------------------------------------- %
%% Settings:

lineWidth1  = 5; % actual lines in plot
lineWidth2  = 3; % all other lines
fontSize    = 32; % 24
transp      = 0.10;

% ----------------------------------------------------------------------- %
%% Until when data are completely available:

until       = find(round(data.mu.time, 3) == 1); % end for outcome-locked

% ----------------------------------------------------------------------- %
%% Determine baseline (i.e. lowest line at time point zero):

iBaseline   = round(data.mu.time, 3) == 0; % find baseline via time in sec
baseline    = min(squeeze(nanmean(data.SubCondTime(job.validSubs, :, iBaseline)))); % still extract only valid subjects

% ----------------------------------------------------------------------- %
%% Start plot:

p = cell(nCond, 1); % delete p

% figure('units', 'normalized', 'outerposition', [0 0 1 1]); hold on % fullscreen
figure('Position', [100 100 1200 800]); hold on

for iCond = 1:job.nCond
       
    % With error bars:
    if errorbar
        p{iCond} = boundedline(data.mu.time(1:until),...
            squeeze(nanmean(data.SubCondTime(job.validSubs, iCond, 1:until))) - baseline, ...
            job.nCond / (job.nCond - 1) * squeeze(nanstd(squeeze(data.SubCondTime(job.validSubs, iCond, 1:until)) - ...
            data.SubTime(job.validSubs, 1:until) + repmat(data.GrandTime(1:until), job.nValidSubs, 1)))' ./ ...
            sqrt(job.nValidSubs), ...
            'cmap', job.colMat(iCond, :), 'alpha', 'transparency', transp);
    else
        
%     Without error bars:
        p{iCond} = plot(data.mu.time(1:until),...
            squeeze(nanmean(data.SubCondTime(job.validSubs, iCond, 1:until))) - baseline,...
            'cmap', job.colMat(iCond, :), 'alpha', 'transparency', transp);
    end
    
    set(p{iCond}, 'linestyle', job.lineStyle{iCond}); % adjust style (solid, dashed)
    set(p{iCond}, 'Linewidth', lineWidth1); % adjust line width

end

% ----------------------------------------------------------------------- %
%% General axis settings:

% Labels:
xlabel('Time (s)', 'fontsize', fontSize)

set(gca, 'linewidth', lineWidth2);
set(gca, 'fontsize', fontSize);
yLim    = get(gca, 'ylim'); % retrieve current y axis limits, overwrite later
yMinLim = yLim(1); 
yMaxLim = yLim(2); 

if isfield(job, 'freq')
    title(sprintf('Mean time-frequency power %.1d-%.1d Hz over %s', ...
        job.freq(1), job.freq(2), strjoin(job.channels, '/')), 'fontsize', fontSize);
    ylabel('Power (dB)', 'fontsize', fontSize, 'fontweight', 'bold');
    yMinLim = yMinLim - 0.1; % adjust a bit
    yMaxLim = yMaxLim + 0.1; % adjust a bit

else
    title(sprintf('Voltage over %s', strjoin(job.channels, '/')), 'fontsize', fontSize);
    ylabel('Amplitude (A/cm²)', 'fontsize', fontSize); 
end

% ----------------------------------------------------------------------- %
%% Add x-axis labels, vertical lines:

set(gca, 'xlim', [-0.25 1.0], 'ylim', [yMinLim yMaxLim],...
    'xtick', [-1.0 -0.75 -0.5 -0.25 0 0.25 0.5 0.70 1 1.25], ...
    'xtickLabel', {'-1000', '-750', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1,000'}, ...
    'fontsize', fontSize, 'Linewidth', lineWidth2) %, 'ytick',-5:0.1:4)

% ----------------------------------------------------------------------- %
%% Extra vertical lines:

plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', lineWidth2); % outcome onset
plot([.7 .7], get(gca, 'ylim'), ':k', 'LineWidth', lineWidth2); % outcome offset

% ----------------------------------------------------------------------- %
%% Add area of significance:

legendNames = job.condNames;
if isfield(job, 'sigTime') && ~isempty(job.sigTime) 
    p{end+1} = plot([job.sigTime(1) job.sigTime(end)], [yMaxLim-0.1 yMaxLim-0.1], 'k-', 'linewidth', lineWidth1);
    legendNames{end+1} = 'Significant difference';
end

% ----------------------------------------------------------------------- %
%% Add legend:

legend([p{:}], legendNames); legend boxoff;

end % END OF FUNCTION.