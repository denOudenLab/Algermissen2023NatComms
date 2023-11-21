function twoLinePlot(job, data, errorbar, within)

% twoLinePlot(job, data, errorbar, within)
%
% Creates lineplot for two conditions specified as contrast in job.
% Works for both ERP and TFR data.
%
% INPUTS:
% job                   = cell with relevant fields for plotting settings, 
% required fields:
%   .validSubs          = numeric vector, subjects to include.
%   .nValidSubs         = numeric, number of subjects to include.
%   .twoLineLabels      = vector of 2 strings, labels for conditions.
%   .freq               = numeric vector of 2 elements, range of selected
%   frequencies.
%   .channels           = vector of strings, selected channel names.
%   .sigTime            = numeric vector of 2 elements, timing for
%   significant difference to highlight (optional).
% data                  = cell with the following fields:
%   .mat1 and mat2      = averaged over channels/ frequencies/ conditions.
%   .matMean            = mean of mat1 and mat2 (over all conditions).
%   .matGrandMean       = matMean averaged over subjects.
% errorbar              = add errorbars (true) or not (false) to plot
% (default: true):
% within                = perform Cousineau-Morey correct on error bars
% (true) or not (default: true)
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

if nargin < 4
    within = true;
end

% ----------------------------------------------------------------------- %
%% Close any open figure windows:

close all

% ----------------------------------------------------------------------- %
%% Fixed settings:

fontSize        = 32; % 24; % 12 
fontSizeTitle   = 15; % keep overall title at 15
lineWidth1      = 5; % actual lines in plot
lineWidth2      = 3; % all other lines

% ----------------------------------------------------------------------- %
%% Determine baseline (i.e. lowest line at time point zero):

until           = find(round(data.mu.time, 3) == 1); % for outcome-locked
iBaseline       = find(round(data.mu.time, 3) == 0); % find baseline via time in sec
baseline        = min(nanmean(data.mat1(:, iBaseline)), nanmean(data.mat2(:, iBaseline)));

% ----------------------------------------------------------------------- %
%% Plot:

p   = cell(2, 1);

figure('Position', [100 100 1200 800]); hold on % DCC screen

% Lines:
% Use Cousineau 2005 method for uncorrected-within-subjects CIs 
% (only correction, no z; i.e. approx SE with only 2 conditions):

if errorbar % with error bars:

    if within % within-subjects (use Cousineau/Morey method): subtract individual mean, add group level mean, correct for #conds bz x/(x-1)
        p{1} = boundedline(data.mu.time(1:until), ...
            nanmean(data.mat1(:, 1:until)) - baseline, ...
            2 / (2-1) * nanstd(data.mat1(:, 1:until) - ...
            data.matMean(:, 1:until) + ... 
            repmat(data.matGrandMean(1:until), job.nValidSubs, 1)) ./ ...
            sqrt(job.nValidSubs), ...
            'cmap', job.twoLineColor(1, :), 'alpha', 'transparency', 0.2); %
        set(p{1}, 'linewidth', lineWidth1); % adjust linewidth (not possible within boundedline)
        p{2} = boundedline(data.mu.time(1:until), ...
            nanmean(data.mat2(:, 1:until)) - baseline, ...
            2 / (2-1) * nanstd(data.mat2(:, 1:until) - ...
            data.matMean(:, 1:until) + ...
            repmat(data.matGrandMean(1:until), job.nValidSubs, 1)) ./ ...
            sqrt(job.nValidSubs), ...
            'cmap', job.twoLineColor(2, :), 'alpha', 'transparency', 0.2); %
        set(p{2}, 'linewidth', lineWidth1); % adjust linewidth (not possible within boundedline)

    else % between-subjects (standard way of compute SEs):
        p{1} = boundedline(data.mu.time(1:until), ...
            nanmean(data.mat1(:, 1:until)), ...
            nanstd(data.mat1(:, 1:until)) ./ sqrt(length(job.validSubs)), ...
            'cmap', job.twoLineColor(1, :), 'alpha', 'transparency', 0.2);
        p{2} = boundedline(data.mu.time(1:until), ...
            nanmean(data.mat2(:, 1:until)), ...
            nanstd(data.mat2(:, 1:until)) ./ sqrt(length(job.validSubs)), ...
            'cmap', job.twoLineColor(2, :), 'alpha', 'transparency', 0.2);
    end
    
else % without error bars:
    p{1} = plot(data.mu.time, squeeze(nanmean(data.mat1)), ...
        '-', 'linewidth', lineWidth1, 'Color', job.twoLineColor(1, :)); % data.mat1: blue
    p{2} = plot(data.mu.time, squeeze(nanmean(data.mat2)), ...
        '-', 'linewidth', lineWidth1, 'Color', job.twoLineColor(2, :)); % data.mat2: red
end

% ----------------------------------------------------------------------- %
%% Extra vertical lines:

plot([0 0], get(gca,'ylim'), ':k', 'LineWidth', 3); % outcome onset
plot([.7 .7], get(gca,'ylim'), ':k', 'LineWidth', 3); % outcome offset

% ----------------------------------------------------------------------- %
%% x-axis label:

xlabel('Time (ms)', 'fontsize', fontSize, 'fontweight', 'bold', 'linewidth', lineWidth2);

% ----------------------------------------------------------------------- %
%% Retrieve y-axis limits, change later:

yLim    = get(gca, 'ylim'); % retrieve automatic settings, change later
yMinLim = yLim(1); 
yMaxLim = yLim(2); 

% ----------------------------------------------------------------------- %
%% Add x-axis settings and labels, add vertical lines:

set(gca, 'xlim', [-0.25 1.0], ...
    'ylim', [yMinLim yMaxLim], ...
    'xtick', [-1.0 -0.75 -0.5 -0.25 0 0.25 0.5 0.70 1 1.25], ...
    'xtickLabel',{'-1000', '-750', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1, 000'}, ...
    'fontsize', fontSize, 'linewidth', lineWidth2);

% ----------------------------------------------------------------------- %
%% Add title and y-label:

if isfield(job, 'freq')
    
    title(sprintf('%s vs. %s: Time-frequency power %d-%d Hz over %s', ...
        job.twoLineLabels{1}, job.twoLineLabels{end}, ...
        job.freq(1), job.freq(2), strjoin(job.channels, '/')), ...
        'fontsize', fontSizeTitle);
    ylabel('Power (dB)', 'fontsize', fontSize, 'fontweight', 'bold', 'fontsize', fontSize);
    yMinLim = yMinLim - 0.1; % adjust a bit
    yMaxLim = yMaxLim + 0.1; % adjust a bit
    set(gca, 'ylim', [yMinLim yMaxLim]); % set anew
    boost   = 0.05; % boost: negative offset from yMaxLim to draw significance lines
    
else    
    if isfield(job, 'ROIs2use')
        title(sprintf('Split per %s BOLD: Voltage over %s', ...
            job.ROIname, strjoin(job.channels, '/')), ...
            'fontsize', fontSizeTitle)
    
    else
        title(sprintf('%s vs. %s: Voltage over %s', ...
            job.twoLineLabels{1}, job.twoLineLabels{end}, strjoin(job.channels, '/')), ...
            'fontsize', fontSizeTitle)
    end
    
    ylabel('Amplitude (A/cm²)', 'fontweight', 'bold', 'fontsize', fontSize); 
    boost   = 0.0015; % boost: negative offset from yMaxLim to draw significance lines

end

% ----------------------------------------------------------------------- %
%% Add area of significance:

legendNames = job.twoLineLabels; % initialize labels

if isfield(job, 'sigTime') && ~isempty(job.sigTime) % add line of significance
    
    nSig = length(job.sigTime)/2; % job.sigTime allows for pairs of number from where to where difference is significant
    
    for iSig = 1:nSig % loop over pairs
        
        fprintf('Found significant timing at %.03f - %.03f sec\n', ...
            job.sigTime(iSig * 2 - 1), job.sigTime(iSig * 2));
        p{end+1} = plot([job.sigTime(iSig*2-1) job.sigTime(iSig*2)], ...
            [yMaxLim-boost yMaxLim-boost], 'k', 'linewidth', lineWidth1);

    end
    
    legendNames{end+1} = 'Significant difference'; % update if necessary
    
end

% ----------------------------------------------------------------------- %
%% Add legend:

legend([p{:}], legendNames, 'fontsize', fontSize); legend boxoff

end % END OF FUNCTION.