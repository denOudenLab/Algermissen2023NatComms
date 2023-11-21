function [subj] = sim_subj()

% sim_subj()
% 
% Initialize data for one "standard" subject for model simulations.
%
% INPUTS:
% none.
%
% OUTPUTS:
% subj          = structure with the following fields:
% .stimuli      = vector of integers, stimuli seen (1-16).
% .reqactions   = vector of integers, required action given stimulus (1 = Go, 0 = NoGo).
% .feedback 	= vector of integers, feedback validity (1 = valid, 0 = invalid). 
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/CBM_Scripts/

% clear all; close all; clc

% ----------------------------------------------------------------------- %
%% Settings:

nStim   = 16;   % number of stimuli.
nResp   = 3;    % number of response options.
nTrial  = 640;  % number of trials.
fbProb  = 0.8;  % probability of valid feedback.

% ----------------------------------------------------------------------- %
%% Stimulus IDs:

% {1 = Go2Win-left, 2 = Go2Win-right, 3 = Go2Avoid-left, 4 = Go2Avoid-right, 
% 5 = NoGo2Win, 6 = NoGo2Win, 7 = noGo2Avoid1, 8 = noGo2Avoid2}.

rewID         = [1, 2, 5, 6, 9, 10, 13, 14]; % reward stimuli
goID          = [1:4, 9:12]; % Go stimuli
goLeftID      = [1, 3, 9, 11]; % Go-left stimuli
goRightID     = [2, 4, 10, 12]; % Go-right stimuli

% ----------------------------------------------------------------------- %
%% Create stimulus sequence (vector with all stimuli in presented order):

nTrialPerStim = nTrial / nStim; % number of trials per stimulus: 40
stimSeq       = repelem(1:nStim, nTrialPerStim); % vector with all stimulus presented nTrialPerStim times, NOT random!

% ----------------------------------------------------------------------- %
%% Categorize trial type:

% valenceSeq    = ismember(stimSeq, rewID); % win = 1, pun = 0. % actually later not needed
% goSeq         = ismember(stimSeq, goID); % 1 = go, 0 = nogo. % actually later not needed
actionSeq     = 3 - 2*ismember(stimSeq, goID); % 1 = go, 3 = nogo.
actionSeq(ismember(stimSeq, goRightID)) = 2; % 1=go-left, 2=go-right, 3 = nogo. 
% tabulate(actionSeq)

% ----------------------------------------------------------------------- %
%% Sample feedback validity:

feedback = nan(1, nTrial);
for iTrial = 1:nTrial
    feedback(iTrial) = binornd(1, fbProb);
end

% ----------------------------------------------------------------------- %
%% Save into subj object:

subj.stimuli    = stimSeq';
subj.reqactions = actionSeq';
subj.feedback   = feedback';

% END
