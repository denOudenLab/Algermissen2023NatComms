% EEGfMRIPav_cbm_prepareData.m
% 
% Execute this script to prepare the behavioral data for an appropriate format for the CBM toolbox.
% Mind adjusting the root directory.
%
% INPUTS:
% none.
%
% OUTPUTS:
% Saves all data as EEGfMRIPav_cbm_inputData.mat to indicated directory.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/CBM_Scripts/

% clear all; close all; clc

% ----------------------------------------------------------------------- %
%% Initialize directories:

dirs        = [];
dirs.root   = '/project/3017042.02'; % adjust to local folder structure

dirs.behav  = fullfile(dirs.root, 'Log/Behavior/Data_beh_mat');
dirs.target = fullfile(dirs.root, 'Log/Behavior/Modelling_Stan/CBM_Results');
if ~exist(dirs.target, 'dir'); mkdir(dirs.target); end

% ----------------------------------------------------------------------- %
%% Settings:

fprintf('Initialize settings\n')

% Fixed input settings:
nSessions       = 2; % number sessions, here 2
nCueConditions  = 4; % number cue conditions, here 4: Go2Win, Go2Avoid, NoGo2Win, NoGo2Avoid
nExemplars      = 2; % number of cues in each cue condition, here 2 because Go conditions need a left Go and a right Go cue
nPresentations  = 40; % number presentations of each cue, here 40

% Downstream settings:
nCues           = nCueConditions * nExemplars * nSessions; % total number of used cues, here 16
nTrials         = nCueConditions * nExemplars * nPresentations; % number trials per session, should be 320

% ----------------------------------------------------------------------- %
%% Retrieve subject IDs:

% Extract subject numbers:
subFind = dir(fullfile(dirs.behav, '*results*')); % extract all file names within directory
subFind = {subFind.name}; % keep only column "name"
subList = nan(numel(subFind), 1);

% Retrieve names of subject fiels:
for iSub = 1:numel(subFind) % 1 till number of elements in subFind
    subList(iSub) = str2double(string(extractBetween(subFind{iSub}, 'dij_', '_001_results'))); 
    % extracts the actual "sub" number of each subject --> check which ones found
end

sID     = sort(subList); % create a final subjectlist
nSub    = numel(sID); % number of subjects

% ----------------------------------------------------------------------- %
%% Load rawData:

fprintf('Load raw data\n')

for iSub = 1:nSub % for each subject
    
    fprintf('Load data for subject %03d\n', iSub)
    load(fullfile(rawDataPath_mat, subFind{iSub})); % load rawData for that subject
    
    for iSes = 1:2 % for each session
        
        rawData.sub(iSub, iSes, :)      = repelem(sID(iSub),nTrials); % subject number
        rawData.ses(iSub, iSes, :)      = repelem(iSes,nTrials); % session number
        rawData.trial(iSub, iSes, :)    = 1:nTrials; % trial number
        rawData.stim(iSub, iSes, :)     = prep.seq.learn.stim{iSes}; % stimulus presented
        rawData.corresp(iSub, iSes, :)  = prep.seq.learn.resp{iSes}; % correct response
        rawData.RT(iSub, iSes, :)       = results.learn{iSes}.RT; % RT per trial, NaN for nogo
        rawData.acc(iSub, iSes, :)      = results.learn{iSes}.acc; % 1 for correct, 0 for incorrect
        rawData.resp(iSub, iSes, :)     = results.learn{iSes}.response; % pressed response key (0 for nogo, 37 = left, 39 == right)
        rawData.out(iSub, iSes, :)      = results.learn{iSes}.outcome; % 1 for reward, 0 for neutral, -1 for punishment
        rawData.go(iSub, iSes, :)       = results.learn{iSes}.go; % 1 for (any) go, 0 for nogo
        
    end
end % takes only a few seconds

% ----------------------------------------------------------------------- %
%% Recode variables:

% Overview:
% 1. Code all responses with RT > 1.3035 as NoGos 
% 2. Recode Go so that only absence of any response (0) counts as NoGo, everything else as Go
% 4. Code response given on side of box (left/right), not precise button
% 4. Code accuracy given on side of box (left/right), not precise button
% Check correct responses:
% tabulate(rawData.corresp(:)) % 0 for NoGo, 101 for left Go, 97 for right Go
% Check possible responses:
% tabulate(rawData.resp(:)) % also 65, 69, 98, 102
% 5. Delete RTs < 0.2
% 6. Delete RTs > 1.3035
% 7. Outcome NaN becomes neutral outcome

fprintf('Recode raw data\n')

for iSub = 1:nSub % for each subject:
    fprintf('Recode data for subject %03d\n', iSub)
    
   % a) Let stimuli of second session range from 9 to 16:
    rawData.stim(iSub, 2, :) = rawData.stim(iSub, 2, :) + 8;
    
    for iSes = 1:nSessions % for each session
        
        for iTrial = 1:nTrials % for each trial
            
            % 1. Recode Resp: set to NoGo if RT > 1.3035:
            if rawData.RT(iSub, iSes, iTrial) > 1.3035
                rawData.resp(iSub, iSes, iTrial)        = 0;
            end
            
            % 2. Recode Go: NoGo only if really no response, otherwise Go:
            if rawData.resp(iSub, iSes, iTrial) == 0
                rawData.go(iSub, iSes, iTrial)          = 0;
            else
                rawData.go(iSub, iSes, iTrial)          = 1;  
            end
            
            % 3. Recode true handedness response: take any left or any right button press as corresponding to left/ right correct response:
            if ismember(rawData.resp(iSub, iSes, iTrial), [69, 70, 71, 72, 101, 102, 103, 104]) % correct left Go response
                rawData.resp(iSub, iSes, iTrial)        = 101;
                rawData.handresp(iSub, iSes, iTrial)    = 1;
            elseif ismember(rawData.resp(iSub, iSes, iTrial), [65, 66, 67, 68, 97, 98, 99, 100]) % correct right Go response
                rawData.resp(iSub, iSes, iTrial)        = 97;
                rawData.handresp(iSub, iSes, iTrial)    = 2;
            elseif rawData.resp(iSub, iSes, iTrial) == 0 % correct NoGo response
                rawData.resp(iSub, iSes, iTrial)        = 0;
                rawData.handresp(iSub, iSes, iTrial)    = 3;
            else
                error('Unknown response for Sub %03d, Ses %d, Trial %d', iSub, iSes, iTrial);
            end
            
            % 4. Recode accuracy: take any left or any right button press as corresponding to left/ right correct response:
            if rawData.corresp(iSub, iSes, iTrial) == rawData.corresp(iSub, iSes, iTrial) % correct left Go response
                rawData.acc(iSub, iSes, iTrial) = 1;
            else % otherwise incorrect
                rawData.acc(iSub, iSes, iTrial) = 0; 
            end
            
            % 5. and 6. Recode RTs: delete RT < 0.2 or RT > 1.3035:
            if rawData.RT(iSub, iSes, iTrial) < 0.2 || rawData.RT(iSub, iSes, iTrial) > 1.3035 
                rawData.RT(iSub, iSes, iTrial) = NaN;
            end
            
            % 7. Outcomes NaN recoded to neutral outcome:
            if isnan(rawData.out(iSub, iSes, iTrial))
                rawData.out(iSub, iSes, iTrial) = 0;
            end
        end
    end
end % takes only a few seconds

% Check:
% sum(rawData.RT(:) < 0.2)
% sum(rawData.RT(:) > 1.3035)
% sum(isnan(rawData.out(:)))
% Statistics:
% tabulate(rawData.stim(:))
% tabulate(rawData.resp(:))
% tabulate(rawData.handresp(:))
% tabulate(rawData.out(:))

% ----------------------------------------------------------------------- %
%% Reshape into structure:

fprintf('Reshape data\n');

% Initialize empty 36 x 1 structure:
data = struct([]);

% Fill structure:
for iSub = 1:nSub
    
    data{iSub}.stimuli = [squeeze(rawData.stim(iSub, 1, :)); squeeze(rawData.stim(iSub, 2, :))];
    data{iSub}.actions = [squeeze(rawData.handresp(iSub, 1, :)); squeeze(rawData.handresp(iSub, 2, :))];
    data{iSub}.outcome = [squeeze(rawData.out(iSub, 1, :)); squeeze(rawData.out(iSub, 2, :))];
    
end

% ----------------------------------------------------------------------- %
%% Save as one big rawData structure:

fprintf('Save data\n')
outputFile = fullfile(dirs.target, 'EEGfMRIPav_cbm_inputData.mat');    
save(outputFile, 'data');

% END OF FILE.
