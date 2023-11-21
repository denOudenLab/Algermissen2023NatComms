function out = load_recode_behav(rootDir, iSub, rejectedtrials)

% out = load_recode_seq(rootDir, iSub, rejectedtrials) 
%
% Load behavioral data; 
% recode response and accuracy variables (uninstructed keys to
% instructed keys; delete Gos with too long RTs, accuracy based on recoded
% keys presses);
% filter out rejected trials (rejected in EEG trial rejection).
%
% INPUT:
% rootDir           = string, root directory of project (needed to find
% behavioral data)
% iSub              = scalar integer, number of subject which to retrieve
% data from
% rejectedtrials    = vector of length 1 x number of valid (!) trials containing
% whether trials rejected (1) or not (0) (note that during EEG epoching, invalid trials are not even epoched because different trigger).
%
% OUTPUT:
% out               = pre-processed behavioral data (length: number of valid, non-rejected trials), with following fields:
% .trlIdx           = numeric, indices of trials kept after EEG rejection
% .stim             = numeric, identifiers of stimulus identify (1-16)
% .goCue            = numeric, identifiers of Go cues
% .reqaction        = numeric, index whether trial is Go cue (1) or not (0)
% .winCue           = numeric, identifiers of Win cues
% .val              = numeric, index whether trial is Win cue (1) or not
% (0)
% .go               = numeric, Go (1) or NoGo (1) performed; old variable
% .ISI              = numeric, inter-stimulus interval (between cue and outcome), in seconds.
% .ITI              = numeric, inter-trial interval (between cue and outcome), in seconds.
% .
% from raw data; rather use recoded isgo variable
% .resp             = numeric, left Go (101) or right Go (97) or NoGo (0)
% performed (recoded: uninstructed keys to instructed keys, remove too long
% RTs).
% .isgo             = numeric, Go (1) or NoGo (1) performed (recoded: 
% uninstructed keys to instructed keys, remove too long RTs).
% .reqresp          = numeric, left Go (101) or right Go (97) or NoGo (0)
% required.
% .reqgo            = numeric, Go (1) or NoGo (1) required.
% .accuracy         = numeric, correct (1) or incorrect (0) based on
% reqresp and resp.
% .RT              	= numeric, RTs, NaN if too fast (< 0.2 sec.) or too
% slow (> 1.3035) or NoGo.
% .outcome 	    = numeric, absolute outcome obtained, either reward (1), neutral (0), punishment (-1).
% .fb.abs 	    = numeric, same as .outcome.
% .fb.rel 	    = numeric, relative outcome obtained, either positive (1) or negative (2).
% .fb.all 	    = numeric, exact outcome obtained, either reward (1), no-reward (2), no-punishment (3), punishment (4).
% .splithalf        = numeric, whether first (1) or second (2) session.
% 
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% ----------------------------------------------------------------------- %
%% Complete input settings:

if ~exist('rootDir', 'var')
    rootDir = helpers_set_rootDir(); % '/project/3017042.02';
    fprintf('rootDir unspecified, assume %s\n', rootDir);
end

if ~exist('iSub', 'var')
    iSub    = 1;
    fprintf('iSub unspecified, assume %d\n', iSub);
end

if exist('rejectedtrials', 'var')
    isRej   = true; % fill in later once number valid trials is loaded 
else 
    isRej   = false; 
    warning('No rejected trials provided, keep all valid trials');
end

% ----------------------------------------------------------------------- %
%% Load raw behavioral data:

fprintf('Subject %03d: Load seqioral data\n', iSub);
load(fullfile(rootDir, sprintf('Log/Behavior/Data_beh_mat/3017042.02_emmvdij_%03d_001_results.mat', iSub)));
% loads prep with task settings
% loads results with results
% loads tm with timings

% ----------------------------------------------------------------------- %
%% Extract behavioral data, combine both sessions:

fprintf('Subject %03d: Combine both sessions of seqioral data\n', iSub);

seq.stim        = [prep.seq.learn.stim{1}; prep.seq.learn.stim{2}];
seq.resp        = [results.learn{1}.response; results.learn{2}.response];
seq.RT          = [results.learn{1}.RT; results.learn{2}.RT];
seq.outcome     = [results.learn{1}.outcome; results.learn{2}.outcome];
seq.accuracy    = [results.learn{1}.acc; results.learn{2}.acc];
seq.go          = [results.learn{1}.go; results.learn{2}.go];
seq.splithalf   = [ones(160, 1); 2*ones(160, 1); ones(160, 1); 2*ones(160, 1)];

% ----------------------------------------------------------------------- %
%% Extract event timings:

tCue    = [tm.learn{1}.stim; tm.learn{2}.stim];
tResp   = [tm.learn{1}.response; tm.learn{2}.response];
tSFI    = [tm.learn{1}.sfi; tm.learn{2}.sfi];
tFb     = [tm.learn{1}.outcome; tm.learn{2}.outcome];
tWarn   = [tm.learn{1}.warning; tm.learn{2}.warning];
tFb(isnan(tFb)) = tWarn; % insert warnings into fb timing

% Extract ISI and ITI:

% tSFI - tCue % mostly 1.3166
seq.ISI = tFb - tCue - 1.3166;
% fprintf('ISI statistics: min = %.03f, max = %.03f, mean = %.03f\n', min(seq.ISI), max(seq.ISI), mean(seq.ISI)); 
seq.ITI = [prep.seq.learn.ITI{1}; prep.seq.learn.ITI{2}];
% fprintf('ITI statistics: min = %.03f, max = %.03f, mean = %.03f\n', min(seq.ITI), max(seq.ITI), mean(seq.ITI)); 

% ----------------------------------------------------------------------- %
%% Indices of valid trials (no invalid outcomes):

fprintf('Subject %03d: Found %d trials with invalid trials---to-be-dropped\n', iSub, sum(isnan(seq.outcome)));
trlIdx = find(~isnan(seq.outcome)); % trial indices of valid trials

% trlIdx has length of number of valid trials. 
% rejectedtrials contains whether trial is rejected (1) or not (0).
% Now apply rejectedtrials (of same length because invalid trials not even 
% epoched in EEG) on old trlIdx to make it potentially even shorter:

% In case rejected trials provided:
if isRej

    % Check if matching sizes:
    if length(rejectedtrials) ~= length(trlIdx); error('Rejected trials and trlIdx contain different numbers of trials'); end

    % Now reject marked trials:
    fprintf('Subject %03d: Found %d rejected trials---to-be-dropped\n', iSub, sum(rejectedtrials));
    trlIdx = trlIdx(~rejectedtrials);
    % trlIdx contains indices of trials (trial number relative to 1:640).
    % Can now be applied to behavior (or fMRI).

end

% trlIdx contains indices of trials (trial number relative to 1:640).
% Can now be applied to behavior (or fMRI).

% ----------------------------------------------------------------------- %
%% Transfer relevant timings:

out.ISI     = seq.ISI(trlIdx);
out.ITI     = seq.ITI(trlIdx);

% ----------------------------------------------------------------------- %
%% Recode response variables:

% Initialize empty variables:
nTrial          = length(seq.resp); 
out.reqresp     = nan(nTrial, 1); % initialize in same format 
out.resp        = nan(nTrial, 1); % initialize in same format 
out.reqgo       = nan(nTrial, 1); % initialize in same format 
out.isgo        = nan(nTrial, 1); % initialize in same format 
out.accuracy    = nan(nTrial, 1); % initialize in same format 
out.RT          = nan(nTrial, 1); % initialize in same format 

for iTrial = 1:nTrial % k = 1

    % ------------------------------------------------------------------- %
    %% a) Recode recorded response (also for uninstructed key presses):
    
    if ismember(seq.resp(iTrial), [69, 70, 71, 72, 101, 102, 103, 104]) % left Go response
            out.resp(iTrial) = 101; % actually instructed left Go response is 101
            out.isgo(iTrial) = 1; % any Go
            
    elseif ismember(seq.resp(iTrial), [65, 66, 67, 68, 97, 98, 99, 100]) % right Go response
            out.resp(iTrial) = 97; % actually instructed right Go response is 97
            out.isgo(iTrial) = 1; % any Go
            
    else % NoGo response
            out.resp(iTrial) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
            out.isgo(iTrial) = 0; % any NoGo
    end
    
    % ------------------------------------------------------------------- %
    %% b) Recode required response (based on stimulus number):
    
    if ismember(seq.stim(iTrial), [1, 3]) % left Go response
            out.reqresp(iTrial) = 101; % actually instructed left Go response is 101
            out.reqgo(iTrial)   = 1;
    elseif ismember(seq.stim(iTrial), [2, 4]) % right Go response
            out.reqresp(iTrial) = 97; % actually instructed right Go response is 97
            out.reqgo(iTrial)   = 1;
    else % NoGo response
            out.reqresp(iTrial) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
            out.reqgo(iTrial)   = 0;
    end
    
    % ------------------------------------------------------------------- %
    %% c) Recode RTs:
    
    if out.isgo(iTrial) == 0
        out.RT(iTrial)      = NaN; % NoGo response
        
    elseif seq.RT(iTrial) < 0.2 % if too short: delete
        out.RT(iTrial)      = NaN;
    
    elseif seq.RT(iTrial) > 1.3035 % if too long: delete, set to NoGo
        out.RT(iTrial)      = NaN;
        out.resp(iTrial)    = 0; % late response = NoGo
        out.isgo(iTrial)    = 0; % late response = NoGo
    
    else
        out.RT(iTrial)      = round(seq.RT(iTrial), 3); % round to
        % ms be consistent with EEG sampling rate (in ms)
    end
    
    % ------------------------------------------------------------------- %
    %% d) Recode accuracy based on recorded responses:
    
    if out.reqresp(iTrial) ~= out.resp(iTrial) || seq.RT(iTrial) > 1.3035 % wrong response (key) or after time
        out.accuracy(iTrial) = 0; % 0 if incorrect
    
    else
        out.accuracy(iTrial) = 1; % 1 if correct
    end
    
end % end of for-loop over trials iTrial

% ----------------------------------------------------------------------- %
%% Report number of recoded variables:

fprintf('Subject %03d: Number of recoded responses: %d\n', iSub, sum(seq.resp ~= out.resp)); % # recoded trials
fprintf('Subject %03d: Number of recoded accuracies: %d\n', iSub, sum(seq.accuracy ~= out.accuracy)); % # recoded trials

% ----------------------------------------------------------------------- %
%% Select data only for non-rejected trials:

fprintf('Subject %03d: Remove rejected trials from behavioral data\n', iSub);

out.trlIdx          = trlIdx;

% Stimuli:
out.stim            = seq.stim(trlIdx);
out.goCue           = 1:4;
out.reqaction       = ~ismember(out.stim,out.goCue) + 1; % 1 if yes, 2 if no
out.winCue          = [1 2 5 6];
out.val             = ~ismember(out.stim,out.winCue) + 1; % 1 if yes, 2 if no

% Response (already recoded):
out.go              = seq.go(trlIdx); % old go variable from raw input, carry forward, but don't use
out.resp            = out.resp(trlIdx);
out.isgo            = out.isgo(trlIdx);
out.reqresp         = out.reqresp(trlIdx);
out.reqgo           = out.reqgo(trlIdx);
out.accuracy        = out.accuracy(trlIdx);
out.RT              = out.RT(trlIdx);

% Outcome:
out.outcome             = seq.outcome(trlIdx); % what subject see: -1, 0, 1
out.fb.abs              = seq.outcome(trlIdx); % what subject see: -1, 0, 1
out.fb.rel              = 1 - out.outcome; % 1 becomes good, 2 becomes bad
out.fb.rel(out.val==1)  = out.fb.rel(out.val == 1)+1; % Win cues: one up (0 becomes 1, 1 becomes 2)
out.fb.all              = out.fb.rel+2*(out.val == 2); % 1-4: reward, no-reward, no-punishment, punishment
out.splithalf           = seq.splithalf(trlIdx);

% Session:
out.splithalf       = seq.splithalf(trlIdx);

end % END OF FUNCTION.