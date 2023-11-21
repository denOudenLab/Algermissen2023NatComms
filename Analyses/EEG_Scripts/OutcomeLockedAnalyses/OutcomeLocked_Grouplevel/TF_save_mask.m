function TF_save_mask(job, data, stat)

% TF_save_mask(job, data, stat)
%
% Save mask (T-values above threshold) for given contrast as yielded by
% cluster-based permutation test in Fieldtrip.
% DO NOT average over channels or frequencies;
% run permutation test with cfg.avgoverchan = 'no'; cfg.avgoverfreq = 'no';
% one-sided (only get all-negative or all-positive mask);
% band-limited (only focus on one cluster).
%
% INPUTS:
% job                     = cell, needs at least fields:
%   .dirs.data            = string, directory where all EEG data are.
%   .contrastType         = string, contrast to be used: 'Preferred', 'Action', 'GoPreferred', 'SalientPreferred', 'SalientAction'.
% data                    = cell, need at least the following fields:
%   .mu                   = grand average across both conditions and subjects.
% stat                    = output object of ft_freqstatistics, with following fields:
%   .label                = vector strings, channel labels of channels
% included in permutation test.
%   .mask                 = matrix, 1 if within significant cluster, 0 outside.
%   .posclusterslabelmat  = matrix, integer > 0 if part of respective cluster number, 0 outside.
%   .negclusterslabelmat  = matrix, integer > 0 if part of respective cluster number, 0 outside.
%   .stat                 = matrix of t-values, unthresholded.
%
% OUTPUTS:
% none, save 3D mask and plot to disk.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

% ----------------------------------------------------------------------- %
%% 1) PREPARE MASK:

% Run permutation test with cfg.avgoverchan = 'no'; cfg.avgoverfreq = 'no'; % 
% one-sided (only get all-negative or all-positive mask);
% band-limited (only focus on one cluster)

% Outputs from permutation test used here:
% stat.mask % 1 if within significant cluster, 0 outside
% stat.posclusterslabelmat % > 0 if part of respective cluster number, 0 outside
% stat.negclusterslabelmat % > 0 if part of respective cluster number, 0 outside
% stat.stat % t-values unthresholded

fprintf('Create mask for contrast %s, channels %s, frequency band %s\n', job.contrastType, strjoin(stat.label, ''), Band);

% Threshold t-values based on significance:
statMask    = stat.mask .* stat.stat; % mask [significance: 0/1] * t-values (continuous) * only negative cluster

% Automatically select channel, time and frequency indices of permutation test settings (not result!) performed:
maskChanIdx = find(ismember(data.mu.label, stat.label)); %
maskTimeIdx = dsearchn(data.mu.time', stat.time'); % retrieve indices of timing of stat
maskFreqIdx = dsearchn(data.mu.freq', stat.freq'); % retrieve indices of frequencies of stat

% Create empty mask, insert t-values:
freqMask    = zeros(size(data.mu.powspctrm)); % create 0s as default
freqMask(maskChanIdx, maskFreqIdx, maskTimeIdx) = statMask; % add statMask to it
% freqMask(maskChanIdx, maskFreqIdx, maskTimeIdx) % Display values of mask 

% ----------------------------------------------------------------------- %
%% Impose additional restriction based on time:

% For theta/ thetadelta: delete positive cluster after 600 ms

% selTime = 0.6; delTime = 'before';
% selTime = 0.6; delTime = 'after';

if exist('selTime', 'var')

    fprintf('Delete clusters %s %0.1f sec.\n', delTime, selTime);
    
    selTimeIdx = dsearchn(data.mu.time', 0.6); 
    
    if strcmp(delTime, 'after')
        freqMask(:, :, selTimeIdx:end) = 0; % delete from selTimeIdx onwards
    elseif strcmp(delTime, 'before')
        freqMask(:, :, 1:selTimeIdx) = 0; % delete till selTimeIdx 
    else
        error('Unknown whether delete before or after cutoff')

    end
end

% If negative: flip
if sum(freqMask(:)) < 0; freqMask = freqMask*-1; end

% Plot TF power (averaged over channels):
contourf(data.mu.time, data.mu.freq, squeeze(nanmean(freqMask, 1)), 40, 'linestyle', 'none');
pause(2);
close gcf

% ----------------------------------------------------------------------- %
%% Save:

dirs.mask = fullfile(job.dirs.results, 'TF_Mask');
if ~exist(dirs.mask, 'dir'); mkdir(dirs.mask); end
maskName = sprintf('Mask_%sLong_%s_%s.mat', job.contrastType, strjoin(sort(stat.label), ''), Band);
save(fullfile(dirs.mask, maskName), 'freqMask');

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Load and plot and save:

contrast    = 'ActionLong'; % Preferred Action ActionLong
channels    = 'CzFCzFz'; % CzFCzFz AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz
band        = 'beta'; % theta lowalpha alpha beta deltathetaalpha deltatheta broad

% Load:
dirs.mask = fullfile(job.dirs.data, 'TF_Mask');
load(fullfile(dirs.mask, sprintf('Mask_%s_%s_%s.mat', contrast, channels, band)));

nFreq   = size(freqMask, 2);
timeVec = -1:0.025:1.3;
freqVec = 1:nFreq;

% Plot:
fontSize = 24; lineWidth = 3;
figure('Position',[100 100 800 800]); 
contourf(timeVec, freqVec, squeeze(nanmean(freqMask, 1)), 40, 'linestyle', 'none');
set(gca, 'fontsize', fontSize, 'Linewidth', lineWidth) %
xlabel('Time (ms)', 'fontsize', fontSize, 'fontweight', 'bold');
ylabel('Frequency (Hz)', 'fontsize', fontSize, 'fontweight', 'bold');
title(sprintf('Mask Contrast %s, channels %s, band %s', contrast, channels, band));

% Save:
saveas(gcf, fullfile(dirs.mask, sprintf('Mask_%s_%s_%s.jpg', contrast, channels, band)));
pause(2)
close gcf

end % END OF FUNCTION.