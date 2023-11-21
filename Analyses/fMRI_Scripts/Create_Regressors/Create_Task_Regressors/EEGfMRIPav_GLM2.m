function EEGfMRIPav_GLM2()

% EEGfMRIPav_GLM2
% 
% Create regressors for subject-level GLM.
% Outcome x action conditions (time of outcome), 
% (hand of) response and error (time of response), 
% outcome onset, invalid responses.
% Mind to select respective EEG regressor & GLM number below.
% Use three-column (onset, duration, factor level) format.
% Saves each regressor as separate .txt file in respective
% timings_regressor folder of respective GLM of respective subject.
%
% Regressors of interest are:
% Actual Action (Go/ NoGo) x Outcome valence (Reward/ NoReward/ NoPunishment/ Punishment)
% 1) GoReward
% 2) GoNoReward
% 3) GoNoPunishment
% 4) GoPunishment
% 5) NoGoReward
% 6) NoGoNoReward
% 7) NoGoNoPunishment
% 8) NoGoPunishment
% 
% Response:
% 9) Left Hand response
% 10) Right Hand response
% 11) Performance: Error (correct/incorrect)
% 
% Outcome:
% 12) Outcome Onset
% 13) Invalid trials (uninstructed button press, no outcome).
%
% Mind setting root directory dirs.root
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% clear all; clc

%% Set directories:

rootDir     = '/project/3017042.02'; % root directory--needs to be adapted to users' folder structure

% Settings:
GLMID       = '2';
nSubjects   = 36;
nTrials     = 640;
nRegressors = 13; % for specification of empty regressors

% Directory where raw data are located:
dirs.behav  = fullfile(rootDir, 'Log/Behavior');
dirs.raw    = fullfile(rootDir, 'Data_beh_mat');

%% Detect which files are in behavioral directory:

input           = dir(dirs.raw); % store details of all files there
allFiles        = {input.name}; % extract names and store as cell array
allFiles        = allFiles'; % tranpose
allFiles{1, 1}  = []; % delete '.'
allFiles{2, 1}  = []; % delete '..'   
dirs.subs       = allFiles(~cellfun('isempty', allFiles)); % remove deleted entries
% --> array with all raw data files

%% MR regressors per block:

validTrialsAll = nan(nSubjects, nTrials);

for iSub = 1:nSubjects % length(dirs.subs) % specify for which subjects to create regressors
    
    fprintf('Subject %03d: Start \n', iSub)
 
    % Load event and timing data:
    filename = strcat('3017042.02_emmvdij_',sprintf('%03d', iSub), '_001_results.mat'); % raw data of specific subject (mind that subject ID can be 2 or 3 digits)
    cd(dirs.raw); % go to folder with result files of all subjects.
    load(filename);
    
    % Set directory where regressors will be printed:
    dirs.timings = fullfile(rootDir, 'Log', 'fMRI',...
            sprintf('sub-%03d', iSub),sprintf('GLM%s', GLMID), 'timings_regressors');
    % Create directory if it doesn't exist yet: 
    if ~exist(dirs.timings , 'dir'); mkdir(dirs.timings); end
      
    % ------------------------------------------------------------------- %
    %% Retrieve variables of interest:
    
    % Just concatenate both sessions:
    stimIDAll      = [prep.seq.learn.stim{1}; prep.seq.learn.stim{2}]; % cue condition (1-8)
    respAll        = [results.learn{1}.response; results.learn{2}.response]; % button ID pressed: 69, 70, 71, 72, 101, 102, 103, 104 for left; 65,66,67,68,97,98,99, 100 for right
    respNumAll     = NaN(640, 1); % responses recoded to 1, 2,3
    motorNumAll    = NaN(640, 1); % responses recoded to 1,-1,0
    RTAll          = [results.learn{1}.RT; results.learn{2}.RT]; % RT for determining whether response was too late (if too late, always wrong)
    correspAll     = [prep.seq.learn.resp{1}; prep.seq.learn.resp{2}]; % correct response
    outcomeAll     = [results.learn{1}.outcome; results.learn{2}.outcome]; % outcome obtained: -1, 0, 1
    incorrectAll   = repelem(NaN,length(stimIDAll))'; % incorrect responses; created after responses recoded

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Recode button presses, accuracy, RTs:
    % Recode uninstructed button presses on correct button box as instructed button presses
    % Recode accuracy based on these recoded responses; too late responses always as incorrect
    
    for k = 1:length(stimIDAll)
        if ismember(respAll(k), [69, 70, 71, 72, 101, 102, 103, 104]) % left Go response
                respAll(k) = prep.par.key.left; % actually instructed left Go response is 101
                respNumAll(k) = 1; % recode left Go numerically to 1
                motorNumAll(k) = 1; % left hand
        elseif ismember(respAll(k), [65, 66, 67, 68, 97, 98, 99, 100]) % right Go response
                respAll(k) = prep.par.key.right; % actually instructed right Go response is 97
                respNumAll(k) = 2; % recode right Go numerically to 2
                motorNumAll(k) = -1; % right hand
        else % NoGo response
                respAll(k) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
                respNumAll(k) = 3; % recode NoGo numerically to 3
                motorNumAll(k) = 0; % no hand
        end
        % Accuracy based on recoded responses
        if respAll(k) ~= correspAll(k) || RTAll(k) > 1.3035 % wrong key or too late (1.3035 latest RT) --> incorrect
            incorrectAll(k) = 1; % 1 if incorrect
        else
            incorrectAll(k) = 0; % 0 if correct
        end
    end
    % Don't recode too short or too long RTs because motor response still given
    % Don't recode trials with too long RTs as NoGos (and accuracy respectively) because people might have thought they made a Go response
        
    % Motor responses:
    leftHandAll     = respNumAll == 1;
    rightHandAll    = respNumAll == 2;
    
    % Initialize other task settings:
    nStim          = 16;
    nTrial         = 640;
    stimAll     = [stimIDAll(1:320); 8 + stimIDAll(321:640)];            % add 8 to stimulusID in second session of task
    cueDetected = zeros(1, nStim);                                       % valenced outcome for cue already experienced or not
    cueDetectAll= nan(nTrial, 1);                                        % cue valence known or not per trial
    validOutcome = nan(nTrial, 1);                                       % store trials where outcome was NaN to later delete row in regressor
    for iTrial = 1:nTrial        
        % Retrieve stimulus ID
        stim = stimAll(iTrial);
        % Store whether cue valence experienced in this trial or not
        cueDetectAll(iTrial) = cueDetected(stim);        
        % Retrieve actually received outcome:
        r = outcomeAll(iTrial);
        % "Unmute" stimulus if valence seen for the first time
        if abs(r) == 1; cueDetected(stim) = 1; end
        % Retrieve whether valid or invalid outcome received
        if isnan(r) == 1; validOutcome(iTrial) = 0; else; validOutcome(iTrial) = 1; end

    end
    % Store valid trials (just for check afterwards):
    validTrialsAll(iSub,:) = validOutcome;
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Create timing for regressors:

    % Strategy for concatenated blocks:
    % a) Retrieve number of volumes per block, compute effective fMRITime in GLM as TE * position of volume
    % b) Extract when in fMRITime each blocks starts,
    % c) extract for when block started in behavioral time (set to 0),
    % d) extract timing of trials from behavioral data files and correct for when blocks started in behavioral time,
    % e) add when block started in fMRITime    
    
    % a) Retrieve block lengths in volumes, create fMRI time:
    dirs.sub            = fullfile(rootDir, sprintf('Log/fMRI/sub-%03d', iSub));
    cd(dirs.sub)
    nVolumes            = load(sprintf('Sub%03d_Nvolumes.txt', iSub)); % number of volumes for all 6 blocks
    TE                  = 1.4; % echo time
    fMRITime            = (0:(sum(nVolumes)-1))*TE; % when each volume starts
    % how many volumes have occurred till end of each block, minus the block
    % itself (so rather at the start of each block), +1 to start at index 1:
    fMRIBlockStart      = cumsum(nVolumes)-nVolumes+1; 
    
    % b) Extract when each block started in fMRI time:
    trialBlockLength    = [110 110 100 110 110 100];
    trialTime0          = nan(0, 0); % time of first volume per block
    for i = 1:length(trialBlockLength)
        trialTime0   = [trialTime0; repmat(fMRITime(fMRIBlockStart(i)), 1, trialBlockLength(i))']; % repeat for each volume in block   
    end    
       
    % c) Extract when blocks started in behavioral time:
    % extract timing at start of blocks, repeat 110 (or 100 times), 
    % put together into t0 that will be subtracted from behavioral time:
    % retrieve timing of start of first 3 blocks, repeat for each of the 110 trials:
    t0part1 = reshape(repmat(tm.learn{1}.stim([1 111 221]), 1, 110)', [330 1]); % one long vector with 330 elements
    t0part1 = t0part1(1:320); % drop last 10 entries as block 3 was 10 trials shorter
    % retrieve timing of start of second 3 blocks, repeat for each of the 110 trials:
    t0part2 = reshape(repmat(tm.learn{2}.stim([1 111 221]), 1, 110)', [330 1]);
    t0part2 = t0part2(1:320); % drop last 10 entries as block 6 was 10 trials shorter
    t0      = [t0part1; t0part2] - 10; % concatenate; now subtract 10 seconds that were used for getting MRI signal steady-state (dummy scans)
   
    % d) Extract trial timing in behavioral time, correct for when block started:
    % 1) Retrieve timings of cues corrected for t0 (onset of block):
    tCueAllBlock        = [tm.learn{1}.stim; tm.learn{2}.stim] - t0; % take cue onset; subtract start of each block
    % 2) Retrieve timings of responses corrected for t0 (onset of block):
    tRespAllBlock       = [tm.learn{1}.response; tm.learn{2}.response] - t0; % note that in case of nogo, the cue offset time was recorded.
    % 3) Retrieve timing of outcomes (warning), correct for NaNs when people got error message, corrected for t0:
    tFBPrel        = [tm.learn{1}.outcome; tm.learn{2}.outcome]; % take onset of learning, so far uncorrected
    tWarningAll    = [tm.learn{1}.warning; tm.learn{2}.warning]; % take timing of warning messages, so far uncorrected
    tFBPrel(isnan(tFBPrel)) = tWarningAll; % replace NaNs with timing of warnings
    tFBAllBlock    = tFBPrel - t0; % subtract start of each block --> now corrected
    % sum(isnan(tFBAllBlock)) % check whether NaN due to uninstructed key presses (NAs left)?
    
    % e) Add start of blocks in fMRI time:
    tCueAll         = tCueAllBlock + trialTime0;
    tRespAll        = tRespAllBlock + trialTime0;
    tFBAll          = tFBAllBlock + trialTime0;
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Identify which trial has which event (indices):
    
    % 1) Cues:
    % Split up trials into those with Win vs. Avoid cues and Go vs. NoGo responses:
    % mind cueDetectAll: don't model cues for trials where outcome not yet experienced
    GoRewardTrials          = (ismember(respNumAll, [1 2]) & outcomeAll == 1  & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    GoNoRewardTrials        = (ismember(respNumAll, [1 2]) & outcomeAll == 0 & ismember(stimIDAll, [1 2 5 6]) & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    GoNoPunishmentTrials    = (ismember(respNumAll, [1 2]) & outcomeAll == 0 & ismember(stimIDAll, [3 4 7 8]) & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    GoPunishmentTrials      = (ismember(respNumAll, [1 2]) & outcomeAll == -1  & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    NoGoRewardTrials        = (respNumAll == 3 & outcomeAll == 1  & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    NoGoNoRewardTrials      = (respNumAll == 3 & outcomeAll == 0 & ismember(stimIDAll, [1 2 5 6]) & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    NoGoNoPunishmentTrials  = (respNumAll == 3 & outcomeAll == 0 & ismember(stimIDAll, [3 4 7 8]) & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    NoGoPunishmentTrials    = (respNumAll == 3 & outcomeAll == -1  & cueDetectAll == 1); % indicator whether trial was Win trial with Go response

    % 2) Responses:
    % leftHandAll
    % rightHandAll
    ErrorTrials = incorrectAll == 1;
    
    % 3) Feedback:
    OutcomeTrials = validOutcome == 1; 
    InvalidTrials = validOutcome == 0;  

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %%  Create three-column regressors:
    % 3-columns format:
    % Onset, Duration (1), Value (1 for factorial regressors)

    emptyRegressorMatrix = zeros(1, nRegressors); % initialize empty regressors (default: non-empty, so zero)

    % RESPONSE X FEEDBACK:
    % 1) Go response with reward: Regressor 1
    if sum(GoRewardTrials) == 0
        tGoRew                  = [0 0 0];
        emptyRegressorMatrix(1) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d Block %d has empty regressor tGoRew (1)!\n', iSub, iBlock)
    else
        tGoRew                  = [tFBAll(GoRewardTrials) ones(sum(GoRewardTrials), 2)];
    end
    
    % 2) Go response with no reward: Regressor 2
    if sum(GoNoRewardTrials) == 0
        tGoNoRew                = [0 0 0];
        emptyRegressorMatrix(2) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d Block %d has empty regressor tGoNoRew (2)! \n', iSub, iBlock)
    else
        tGoNoRew                = [tFBAll(GoNoRewardTrials) ones(sum(GoNoRewardTrials), 2)];
    end
    
    % 3) Go response with no punishment: Regressor 3
    if sum(GoNoPunishmentTrials) == 0
        tGoNoPun                = [0 0 0];
        emptyRegressorMatrix(3) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tGoNoPun (3)! \n', iSub)
    else
        tGoNoPun                = [tFBAll(GoNoPunishmentTrials) ones(sum(GoNoPunishmentTrials), 2)];
    end 

    % 4) Go response with punishment: Regressor 4
    if sum(GoPunishmentTrials) == 0
        tGoPun                  = [0 0 0];
        emptyRegressorMatrix(4) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressor tGoPun (4)! \n', iSub)
    else
        tGoPun                  = [tFBAll(GoPunishmentTrials) ones(sum(GoPunishmentTrials), 2)];
    end
    
    % 5) Go response with reward: Regressor 5
    if sum(NoGoRewardTrials) == 0
        tNoGoRew                = [0 0 0];
        emptyRegressorMatrix(5) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressor tNoGoRew (5)!\n', iSub)
    else
        tNoGoRew                = [tFBAll(NoGoRewardTrials) ones(sum(NoGoRewardTrials), 2)];
    end
    
    % 6) Go response with no reward: Regressor 6
    if sum(NoGoNoRewardTrials) == 0
        tNoGoNoRew              = [0 0 0];
        emptyRegressorMatrix(6) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressor tNoGoNoRew (6)! \n', iSub)
    else
        tNoGoNoRew              = [tFBAll(NoGoNoRewardTrials) ones(sum(NoGoNoRewardTrials), 2)];
    end
    
    % 7) Go response with no punishment: Regressor 7
    if sum(NoGoNoPunishmentTrials) == 0
        tNoGoNoPun              = [0 0 0];
        emptyRegressorMatrix(7) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tNoGoNoPun (7)! \n', iSub)
    else
        tNoGoNoPun              = [tFBAll(NoGoNoPunishmentTrials) ones(sum(NoGoNoPunishmentTrials), 2)];
    end 

    % 8) Go response with punishment: Regressor 8
    if sum(NoGoPunishmentTrials) == 0
        tNoGoPun                = [0 0 0];
        emptyRegressorMatrix(8) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressor tNoGoPun (8)! \n', iSub)
    else
        tNoGoPun                = [tFBAll(NoGoPunishmentTrials) ones(sum(NoGoPunishmentTrials), 2)];
    end
    
    % ------------------------------------------------------------------- %
    % 9-10) Motor response. use response-timing:
    
    if sum(leftHandAll == 1) == 0 % if in fact no errors        
        tLeftHand = [0 0 0];
        fprintf('Note: Subject %d has empty regressor leftHand (9)! \n', iSub)
        emptyRegressorMatrix(9) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
    else
        tLeftHand   = [tRespAll(leftHandAll) ones(sum(leftHandAll), 2)]; % time of response
    end
    
   if sum(rightHandAll == 1) == 0 % if in fact no errors        
        tRightHand = [0 0 0];
        fprintf('Note: Subject %d has empty regressor rightHand (10)! \n', iSub)
        emptyRegressorMatrix(10) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
   else
        tRightHand  = [tRespAll(rightHandAll) ones(sum(rightHandAll), 2)]; % time of response
    end

    % 11) Error trials. Note that we take cue onset as timing here (also NoGo responses)
    if sum(ErrorTrials == 1) == 0 % if in fact no errors        
        tError = [0 0 0];
        fprintf('Note: Subject %d has empty regressor Error (11)! \n', iSub)
        emptyRegressorMatrix(10) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
    else
        tError  = [tCueAll(ErrorTrials) ones(sum(ErrorTrials), 2)]; % sum command checks how many incorrect responses
    end
    
    % ------------------------------------------------------------------- %
    % 12) Outcome onset.
    tOutcome       = [tFBAll(OutcomeTrials) ones(sum(OutcomeTrials), 2)];

    % ------------------------------------------------------------------- %
    % 13) Invalid trials. Uninstructed key pressed, thus no outcome, but error message: Regressor 16
    if sum(InvalidTrials) == 0 % if in fact no uninstructed key presses
        tInvalid                = [0 0 0];
        emptyRegressorMatrix(13) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressor Invalid (13)! \n', iSub)
    else
        tInvalid                = [tFBAll(InvalidTrials) ones(sum(InvalidTrials), 2)]; % sum command checks how many incorrect responses
%            fprintf('Note: Subject %d Block %d has NONEMPTY regressor Invalid! \n', iSub, iBlock)
    end

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Save regressor matrices:
    
    % Action x outcome valence:
    save(fullfile(dirs.timings, 'tGoRew.txt'), 'tGoRew', '-ascii');
    save(fullfile(dirs.timings, 'tGoNoRew.txt'), 'tGoNoRew', '-ascii');
    save(fullfile(dirs.timings, 'tGoNoPun.txt'), 'tGoNoPun', '-ascii');
    save(fullfile(dirs.timings, 'tGoPun.txt'), 'tGoPun', '-ascii');
    save(fullfile(dirs.timings, 'tNoGoRew.txt'), 'tNoGoRew', '-ascii');
    save(fullfile(dirs.timings, 'tNoGoNoRew.txt'), 'tNoGoNoRew', '-ascii');
    save(fullfile(dirs.timings, 'tNoGoNoPun.txt'), 'tNoGoNoPun', '-ascii');
    save(fullfile(dirs.timings, 'tNoGoPun.txt'), 'tNoGoPun', '-ascii');

    % Errors:
    save(fullfile(dirs.timings, 'tLeftHand.txt'), 'tLeftHand', '-ascii');
    save(fullfile(dirs.timings, 'tRightHand.txt'), 'tRightHand', '-ascii');
    save(fullfile(dirs.timings, 'tError.txt'), 'tError', '-ascii');

    % Outcomes:
    save(fullfile(dirs.timings, 'tOutcome.txt'), 'tOutcome', '-ascii');
    save(fullfile(dirs.timings, 'tInvalid.txt'), 'tInvalid', '-ascii');
    
    % Save empty regressors:
    csvwrite(fullfile(dirs.timings, 'emptyregressors.txt'),emptyRegressorMatrix);

end % end iSub-loop.

% Check whether any regressor empty per file:
% emptyRegressorMatrix
% 
% Subject 001: Start 
% Note: Subject 1 has empty regressor Invalid (13)! 
% Subject 002: Start 
% Note: Subject 2 has empty regressor Invalid (13)! 
% Subject 003: Start 
% Subject 004: Start 
% Note: Subject 4 has empty regressor Invalid (13)! 
% Subject 005: Start 
% Note: Subject 5 has empty regressor Invalid (13)! 
% Subject 006: Start 
% Note: Subject 6 has empty regressor Invalid (13)! 
% Subject 007: Start 
% Note: Subject 7 has empty regressor Invalid (13)! 
% Subject 008: Start 
% Subject 009: Start 
% Note: Subject 9 has empty regressor Invalid (13)! 
% Subject 010: Start 
% Subject 011: Start 
% Subject 012: Start 
% Note: Subject 12 has empty regressor Invalid (13)! 
% Subject 013: Start 
% Note: Subject 13 has empty regressor Invalid (13)! 
% Subject 014: Start 
% Note: Subject 14 has empty regressor Invalid (13)! 
% Subject 015: Start 
% Note: Subject 15 has empty regressor Invalid (13)! 
% Subject 016: Start 
% Subject 017: Start 
% Note: Subject 17 has empty regressor Invalid (13)! 
% Subject 018: Start 
% Note: Subject 18 has empty regressor Invalid (13)! 
% Subject 019: Start 
% Note: Subject 19 has empty regressor Invalid (13)! 
% Subject 020: Start 
% Note: Subject 20 has empty regressor Invalid (13)! 
% Subject 021: Start 
% Note: Subject 21 has empty regressor Invalid (13)! 
% Subject 022: Start 
% Note: Subject 22 has empty regressor Invalid (13)! 
% Subject 023: Start 
% Note: Subject 23 has empty regressor Invalid (13)! 
% Subject 024: Start 
% Note: Subject 24 has empty regressor tNoGoNoRew (6)! 
% Note: Subject 24 has empty regressor Invalid (13)! 
% Subject 025: Start 
% Subject 026: Start 
% Subject 027: Start 
% Note: Subject 27 has empty regressor Invalid (13)! 
% Subject 028: Start 
% Note: Subject 28 has empty regressor Invalid (13)! 
% Subject 029: Start 
% Subject 030: Start 
% Subject 031: Start 
% Subject 032: Start 
% Note: Subject 32 has empty regressor Invalid (13)! 
% Subject 033: Start 
% Note: Subject 33 has empty regressor Invalid (13)! 
% Subject 034: Start 
% Note: Subject 34 has empty regressor Invalid (13)! 
% Subject 035: Start 
% Subject 036: Start 
% Note: Subject 36 has empty regressor Invalid (13)!

end % end of function.
