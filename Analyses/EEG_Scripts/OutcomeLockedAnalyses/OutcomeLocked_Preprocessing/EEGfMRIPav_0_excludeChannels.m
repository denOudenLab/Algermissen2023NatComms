function EEGfMRIPav_0_excludeChannels(rootDir, iSub)

% EEGfMRIPav_o_excludeChannels(rootDir, iSub)
% 
% Interactive script.
% Load complete data of one subject, inspect with
% EEGlab, manually note down channels to be rejected.
% 
% INPUTS:
% rootDir           = string, root directory of project
% iSub              = scalar integer, subject number of data to inspect
%
% OUTPUTS:
% none, just plotting channels
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% -------------------------------------------------------------------------
% NOTE: RUNS BEST IN MATLAB 2013B!!! NEWER VERSION MIGHT NOT SUPPORT THIS
% EEGLAB CODE!!!

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Preprocessing/

% ----------------------------------------------------------------------- %
%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootDir', 'var')
    rootDir = '/project/3017042.02';
    fprintf('rootDir unspecified, assume %s\n', rootDir);
end

if ~exist('iSub', 'var')
    iSub = 1; 
    fprintf('iSub unspecified, assume %d\n', iSub);
end

% ----------------------------------------------------------------------- %
%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootDir, 'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs    = set_dirs(rootDir);
par     = set_par();

% ----------------------------------------------------------------------- %
%% Add EEGlab:

fprintf('Add EEGlab to path\n');
dirs.eeglab         = fullfile('~/matlabScripts/eeglab12_0_1_0b');
addpath(dirs.eeglab);

% ----------------------------------------------------------------------- %
%% Detect input directories:

dirs.sub    = fullfile(dirs.raw, sprintf('3017042.02_emmvdij_001_%03d', iSub));
EEGfiles    = dir(fullfile(dirs.sub, '*.eeg')); % read all eeg files for this subject
EEGfiles    = {EEGfiles.name};
nBlock      = length(EEGfiles);

fprintf('Subject %03d: Found %d files (blocks)\n', iSub, nBlock);

% ----------------------------------------------------------------------- %
%% Loop over blocks, read in with Fieldtrip:

% Initialize empty objects:
blockData   = cell(nBlock, 1);
trialCount = 0;                     % count trials to explicitly set trial numbers

for iBlock = 1:nBlock

    fprintf('Subject %03d: Load block %d with Fieldtrip \n', iSub, iBlock);
    inputFile               = EEGfiles{iBlock}; % select respective block

    % Define trials:
    cfg                     = [];
    cfg.dataset             = inputFile; % full path of file
    cfg.trialdef            = [];
    cfg.trialdef.eventvalue = par.epoch.eventCode; % triggers to expect
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.prestim    = abs(par.epoch.epochtime(1)); % time to select before trigger
    cfg.trialdef.poststim   = abs(par.epoch.epochtime(2)); % time to select after trigger
    cfg.event               = ft_read_event(dirs.fileArrayBlock);
    cfg                     = ft_definetrial(cfg);

    % Read in epochs:
    blockData{iBlock}       = ft_preprocessing(cfg); 
    
    % Adjust trial numbers (trialinfo):
    blockLength                 = length(blockData{iBlock}.trialinfo);
    blockData{iBlock}.trialinfo = ((trialCount+1):(trialCount+blockLength))';
    trialCount                  = trialCount + blockLength;

    % Adjust exact timing (sampleinfo):
    % Add 1,000,000 to sampleinfo with increasing block number:
    blockData{iBlock}.sampleinfo = blockData{iBlock}.sampleinfo + 1000000*(iBlock-1);
    
end

% ----------------------------------------------------------------------- %
%% Concatenate blocks:

fprintf('Subject %03d: Concatenate blocks \n', iSub);
cfg = [];
data = ft_appenddata(cfg, blockData{1:nBlocks}); % concatenate blocks

% ----------------------------------------------------------------------- %
%% Initialize EEGlab:

fprintf('Prepare EEGlab \n')

EEG             = [];

% First all settings that stay empty:
EEG.setname     = ''; % give name
EEG.filename    = '';
EEG.filepath    = '';
EEG.subject     = '';
EEG.group       = '';
EEG.condition   = '';
EEG.session     = [];
EEG.comments    = 'Original file: '; 
EEG.nbchan      = []; % #
EEG.trials      = []; % #
EEG.pnts        = []; % # timepoints in trial
EEG.srate       = []; %hz
EEG.xmin        = [];
EEG.xmax        = [];
EEG.times       = []; % 1xpnts array with timings
EEG.data        = []; % chan x time x trial data
EEG.icaact      = [];
EEG.icawinv     = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.chanlocs    = []; 
EEG.urchanlocs  = [];
EEG.chaninfo    = []; 
EEG.ref         = ''; 
EEG.event       = [];
EEG.urevent     = [];
EEG.eventdescription    = {};
EEG.epoch       = [];
EEG.epochdescription    = {};
EEG.reject      = [];
EEG.stats       = [];
EEG.specdata    = [];
EEG.specicaact  = [];
EEG.splinefile  = '';
EEG.icasplinefile       = '';
EEG.dipfit      = [];
EEG.history     = '';
EEG.saved       = 'no';
EEG.etc         = [];

% Prepare data:
EEG.nbchan      = par.chan.nEEGchan; % number of EEG channels (minus ECG, HR, RESP, so 64 left)
EEG.eloc_file   = par.chan.EEG; % exact channels to plot
EEG.srate       = data.fsample; %Sampling rate: 1000 Hz
EEG.trials      = numel(data.trial); % number of trials
EEG.pnts        = size(data.trial{1}, 2); % number timepoints in trial
EEG.times       = data.time{1}; % 1xpnts array with timings
EEG.xmin        = data.time{1}(1); % -1.7 seconds before trial onset 
EEG.xmax        = data.time{1}(end); % 2.8 seconds after trial osnet

% Transfer data from Fieldtrip format to EEGlab format:
for iTrial=1:EEG.trials
    EEG.data(:, :, iTrial) = data.trial{iTrial}(EEG.eloc_file, :); %
end
    
% plot in eeglab for trialrejection.
eeglab redraw
pop_eegplot(EEG, 1, 1, 0);

end % END OF FUNCTION.