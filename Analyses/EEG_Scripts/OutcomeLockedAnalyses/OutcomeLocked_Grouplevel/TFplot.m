function TFplot(job, data, zlim, yScale)

% TFplot(job, data, zlim, yScale)
%
% Creates time-frequency plot (2D heat map) contasting two conditions.
%
% INPUTS:
% job                   = cell with relevant fields for plotting settings, 
% required fields:
%   .TFtiming           = vector of 2 scalars, start and end timing of TF
%   plot (determined usually by lockSettings).
%   .channels           = vector of strings, selected channel names.
%   .chanName           = string, name of selected channel area
%   capitalized.
%   .twoLineLabels      = vector of 2 strings, labels for conditions.
% data                  = cell with the following fields:
%   .TF2plot            = averaged over channels/ conditions/ subjects, for
%   time-frequency plots.
% zlim                  = numeric, positive/ negative limit for for color 
% scale (default: 1).
% yScale        	= string, spacing of the y-axis, either 'log' or 'lin' (default: 'log')
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
    zlim = 1;
end

if nargin < 4
    yScale = 'log';
    fprintf('Unspecified y-axis scale, use %s scaling\n', yScale);
end

% ----------------------------------------------------------------------- %
%% Close any open figure window:

close all

% ----------------------------------------------------------------------- %
%% Fixed settings:

fontSize    = 24;
lineWidth   = 3;

% ----------------------------------------------------------------------- %
%% Plot figure:

figure('name', sprintf('%s vs. %s: Time-frequency power over %s', ...
    job.twoLineLabels{1}, job.twoLineLabels{end}, strjoin(job.channels, '/')), ...
    'Position', [100 100 1200 800]); hold on % fullscreen

% Contour plot:
contourf(data.mu.time, data.mu.freq, data.TF2plot, 40, 'linestyle', 'none'); hold on

% ----------------------------------------------------------------------- %
%% Vertical lines:

plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', lineWidth); % cue onset
plot([0.7 0.7], get(gca, 'ylim'), ':k', 'LineWidth', lineWidth); % cue offset

% ----------------------------------------------------------------------- %
%% Settings:

set(gca, 'xlim', [job.TFtiming], 'ylim', [1 33], 'clim', [-1*zlim 1*zlim], ...
    'xtick', [-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.7 1], ...
    'xtickLabel', {'-1,000', '-750', '-500', '-250', 'OutOn', '250', '500', 'OutOff', '1,000'}, ... 
    'ytick', [2 4 8 16 32], 'yscale', yScale, ...
    'fontsize', fontSize, 'Linewidth', lineWidth); %

box off

% ----------------------------------------------------------------------- %
%% Labels:

xlabel('Time (ms)', 'fontsize', fontSize);
ylabel('Frequency (Hz)', 'fontsize', fontSize);

% ----------------------------------------------------------------------- %
%% Title:

% Title:
title(sprintf('%s vs. %s: \n Time-frequency power over %s (%s-spaced)', ...
    job.twoLineLabels{1}, job.twoLineLabels{end}, ...
    strjoin(job.channels, '/'), yScale), 'fontsize', fontSize);

% ----------------------------------------------------------------------- %
%% Colorbar:

colorbar

end % END OF FUNCTION.