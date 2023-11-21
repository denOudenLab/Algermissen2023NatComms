function EEGfMRIPav_GLM3()

% EEGfMRIPav_GLM3
% 
% Create regressors for block-level GLM.
% Cue x action conditions, response, error, outcome onset, update and
% difference to biased update, invalid responses, trial-by-trial EEG power.
% Mind to select respective EEG regressor & GLM number below.
% Use three-column (onset, duration, factor level) format.
% Saves each regressor as separate .txt file in respective
% timings_regressor folder of respective GLM of respective subject.
%
% Use group-level HBI parameters from M5 (biased learning).
%
% Regressors of interest are:
% Valence (Win/Avoid) x Actual Action (Go/NoGo)
% (Go/NoGo):
% 1) Win_ActGo
% 2) Win_ActNoGo
% 3) Avoid_ActGo
% 4) Avoid_ActNoGo
%
% Covariates of no interest are:
% 5) Handedness of response: hand (left, NoGo, right)
% 6) Performance: Error (correct/incorrect)
% 
% Outcome:
% 7) Outcome Onset
% 8) Standard PE update term
% 9) Difference term to biased PE update
% 10) Invalid trials (uninstructed button press, no outcome).
%
% EEG:
% 11) Trial-by-trial EEG power.
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

% GLMID       = '3A'; EEGPredName = 'Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F38FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta';
% GLMID       = '3B'; EEGPredName = 'Preferred_CzFCzFz_beta';
% GLMID       = '3C'; EEGPredName = 'Action_CzFCzFz_lowalpha';
% only loaded at the very end when 3-column matrices are created and saved!

nSub        = 36;
nRegressors = 11; % for specification of empty regressors

% Directory where raw data are located:
dirs.raw            = fullfile(rootDir, 'Log/Behavior/Data_beh_mat');
dirs.singleTrial    = fullfile(rootDir, 'Log/EEG/OutcomeLockedResults/TF_singleTrial');

%% Detect which files are in behavioral directory:

cd(dirs.raw); % go to folder with result files of all subjects.
input = dir; % store details of all files there
allfiles = {input.name}; % extract names and store as cell array
allfiles = allfiles'; % tranpose
allfiles{1, 1} = []; % delete '.'
allfiles{2, 1} = []; % delete '..'   
dirs.subs = allfiles(~cellfun('isempty', allfiles)); % remove deleted entries
% --> array with all raw data files

%% Initialize group-level model parameters:

% 5 means and 5 standard deviations
modelParameters = [1.0889; -1.7971; 0.1224; 0.5457; 2.4423; 1.5299; 2.4356;
    0.5339; 0.9706; 2.0508];

% Initialize parameters
% Rho and epsilon still untransformed:
rho = exp(modelParameters(1));
epsilon = exp(modelParameters(2))/(1 + exp(modelParameters(2)));

% gobias and pibias not necessary as subjects' actual choices and
% outcomes used for computations; thus no action weights computed

% Biased epsilon:
if epsilon > 0.5
    epsilon_plus_kappa = exp(modelParameters(2) + modelParameters(5)) / (1 + exp(modelParameters(2) + modelParameters(5)));
    epsilon_minus_kappa = epsilon - (epsilon_plus_kappa - epsilon); 
else
   epsilon_minus_kappa = exp(modelParameters(2) - modelParameters(5)) / (1 + exp(modelParameters(2) - modelParameters(5)));
   epsilon_plus_kappa = epsilon + (epsilon - epsilon_minus_kappa);         
end

%% MR regressors per block:

for iSub = 1:nSub % length(dirs.subs) % specify for which subjects to create regressors
    
    fprintf('Starting with subject %d\n', iSub)
 
    % Load event and timing data:
    filename = strcat('3017042.02_emmvdij_', sprintf('%03d', iSub), '_001_results.mat'); % raw data of specific subject (mind that subject ID can be 2 or 3 digits)
    load(fullfile(dirs.raw,filename));
    
    % Set directory where regressors will be printed:
    dirs.timings = fullfile(rootDir, 'Log', 'fMRI',...
                sprintf('sub-%03d', iSub), sprintf('GLM%s',GLMID), 'timings_regressors');
    if ~exist(dirs.timings, 'dir'); mkdir(dirs.timings); end
            
    % ------------------------------------------------------------------- %
    %% Retrieve variables of interest:
    
    % Just concatenate both sessions:
    stimIDAll      = [prep.seq.learn.stim{1}; prep.seq.learn.stim{2}]; % cue condition (1-8)
    respAll        = [results.learn{1}.response; results.learn{2}.response]; % button ID pressed: 69,70,71,72, 101, 102, 103, 104 for left; 65,66,67,68,97,98,99, 100 for right
    respNumAll     = NaN(640, 1); % responses recoded to 1, 2, 3
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

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Simulate from computational model:
    
    % Initialize other task settings:
    nStim          = 16;
    nResp          = 3;
    nTrial         = 640;
    
    % Initialize variables:
    % Define initial value of each stimulus (separate matrix for each stimulus):
    % States:
    stimAll     = [stimIDAll(1:320); 8 + stimIDAll(321:640)];            % add 8 to stimulusID in second session of task
    PavValues   = repmat([0.5,0.5,-0.5,-0.5], [1,4]);                     % cue values, relevant for initializing action values
    stimcount   = zeros(1, nStim);                                       % count number of occurences of stimuli
    cueDetected = zeros(1, nStim);                                       % valenced outcome for cue already experienced or not
    cueDetectAll= nan(nTrial, 1);                                        % cue valence known or not per trial

    % Actions:
    %av0         = zeros(1, nResp, nStim);                                   % initial action values to zero
    av0         = reshape(repelem(rho*PavValues, nResp), 1, nResp, nStim);   % initial action values
    avwithbias          = av0;                                               % initialise current action values av with bias
    avwithoutbias       = av0;                                               % initialise current action values av without bias
    PEwithbias          = nan(nTrial, 1);                                    % initialise PE with bias
    Updatewithbias      = nan(nTrial, 1);                                    % initialise PE with bias
    PEwithoutbias       = nan(nTrial, 1);                                    % initialize PE without bias
    Updatewithoutbias   = nan(nTrial, 1);                                    % initialize PE without bias
    validOutcome        = nan(nTrial, 1);                                    % store trials where outcome was NaN to later delete row in regressor

    % Loop over trials, simulate action values and PEs:
    for iTrial = 1:nTrial
        
        % 1) Retrieve stimulus ID:
        stim = stimAll(iTrial);
        % Number of occurrences for this stimulus:
        stimcount(stim) = stimcount(stim) + 1;
        % Store whether cue valence experienced in this trial or not
        cueDetectAll(iTrial) = cueDetected(stim);
        
        % 2) Make choice:
        % Retrieve actually taken action:
        c = respNumAll(iTrial);
        % Retrieve actually received outcome:
        r = outcomeAll(iTrial);

        % 3) Learn based on outcome:
        % Determine learning rate considering Pavlovian bias on instrumental learning:
        if ismember(c, [1, 2]) && r == 1
            effepsilon = epsilon_plus_kappa;
        elseif c == 3 && r == -1
            effepsilon = epsilon_minus_kappa;
        else
            effepsilon = epsilon;
        end
        
        % Check whether valid key pressed and outcome ever obtained, otherwise no learning:
        if isnan(r) == 1 % no updating if outcome is NaN, i.e. subject got feedback that uninstructed key was pressed
            
            validOutcome(iTrial)        = 0;
            PEwithoutbias(iTrial)       = 0;
            PEwithbias(iTrial)          = 0;            
            Updatewithoutbias(iTrial)   = 0;
            Updatewithbias(iTrial)      = 0;  
            
        else % else: normal Rescorla-Wagner model updating:
            
            validOutcome(iTrial) = 1;
            
            % Update action values without bias (epsilon):
            PEwithoutbias(iTrial)       = rho * r - avwithoutbias(:, c, stim); % compute prediction error as reward minus expected value to be updated
            Updatewithoutbias(iTrial)   = epsilon * PEwithoutbias(iTrial);
            avwithoutbias(:, c, stim)   = avwithoutbias(:, c, stim) + Updatewithoutbias(iTrial); % update respective action value
            
            % Update action values with bias (effepsilon):
            PEwithbias(iTrial)          = rho * r - avwithbias(:, c, stim); % compute prediction error as reward minus expected value to be updated
            Updatewithbias(iTrial)      = effepsilon * PEwithbias(iTrial);
            avwithbias(:, c, stim)      = avwithbias(:, c, stim) + Updatewithbias(iTrial); % update respective action value
            
            % "Unmute" stimulus if valence seen for the first time:
            if abs(r) == 1; cueDetected(stim) = 1; end
        end
    end % end trial loop
    
    % Compute difference between PEs and updates:
    Updatedif   = Updatewithbias - Updatewithoutbias;

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
    fMRITime            = (0:(sum(nVolumes)-1)) * TE; % when each volume starts
    % how many volumes have occurred till end of each block, minus the block
    % itself (so rather at the start of each block), +1 to start at index 1:
    fMRIBlockStart      = cumsum(nVolumes) - nVolumes + 1; 
    
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
    % 1) retrieve timings of cues corrected for t0 (onset of block):
    tCueAllBlock        = [tm.learn{1}.stim; tm.learn{2}.stim] - t0; % take cue onset; subtract start of each block
    
    % 2) retrieve timings of responses corrected for t0 (onset of block):
%     tRespAllBlock       = [tm.learn{1}.response; tm.learn{2}.response] - t0; % note that in case of nogo, the cue offset time was recorded.
    
    % 3) retrieve timing of outcomes (warning), correct for NaNs when people got error message, corrected for t0:
    tFBPrel        = [tm.learn{1}.outcome; tm.learn{2}.outcome]; % take onset of learning, so far uncorrected
    tWarningAll    = [tm.learn{1}.warning; tm.learn{2}.warning]; % take timing of warning messages, so far uncorrected
    tFBPrel(isnan(tFBPrel)) = tWarningAll; % replace NaNs with timing of warnings
    tFBAllBlock    = tFBPrel - t0; % subtract start of each block --> now corrected
    % sum(isnan(tFBAllBlock)) % check whether NaN due to uninstructed key presses (NAs left)?
    
    % e) Add start of blocks in fMRI time:
    tCueAll         = tCueAllBlock + trialTime0;
%     tRespAll        = tRespAllBlock+trialTime0;
    tFBAll          = tFBAllBlock + trialTime0;
    
    % Check by looking at block transitions:
%     i=105:115
%     [tCueAllBlock(i) tCueAll(i)] % transition block 1 to 2 in block time and in total fMRI time
%     (fMRIBlockStart(2)-1)*TE+10 % supposed start of block 2 in total fMRI time

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Identify which trial has which event (indices):
    
    % 1) Cues:
    % Split up trials into those with Win vs. Avoid cues and Go vs. NoGo responses:
    % mind cueDetectAll: don't model cues for trials where outcome not yet experienced
    Win_ActGo_Trials         = (ismember(stimIDAll, [1 2 5 6]) & ismember(respNumAll, [1 2]) & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    Win_ActNoGo_Trials       = (ismember(stimIDAll, [1 2 5 6]) & ismember(respNumAll, 3) & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    Avoid_ActGo_Trials       = (ismember(stimIDAll, [3 4 7 8]) & ismember(respNumAll, [1 2]) & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    Avoid_ActNoGo_Trials     = (ismember(stimIDAll, [3 4 7 8]) & ismember(respNumAll, 3) & cueDetectAll == 1); % indicator whether trial was Win trial with Go response
    
    % 2) Responses:
    % handedness: motorNumAll
    ErrorTrials = incorrectAll == 1;
    
    % 3) Outcome:
    OutcomeTrials = validOutcome == 1; 
    InvalidTrials = validOutcome == 0;  
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %%  Create three-column regressors:
    % 3-columns format:
    % Onset, Duration (1), Value (1 for factorial regressors)

    emptyRegressorMatrix = zeros(1, nRegressors); % initialize empty regressors (default: non-empty, so zero)
            
    % 1) Win ActGo response:
    if sum(Win_ActGo_Trials) == 0
        tWin_ActGo = [0 0 0];
        emptyRegressorMatrix(1) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tWin_ActGo! \n', iSub)
    else
        tWin_ActGo  = [tCueAll(Win_ActGo_Trials  == 1) ones(sum(Win_ActGo_Trials), 2)];
    end
    
    % 2) Avoid ActGo response:
    if sum(Avoid_ActGo_Trials) == 0
        tAvoid_ActGo = [0 0 0];
        emptyRegressorMatrix(2) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tAvoid_ActGo! \n', iSub)
    else
        tAvoid_ActGo  = [tCueAll(Avoid_ActGo_Trials  == 1) ones(sum(Avoid_ActGo_Trials), 2)];
    end
    
    % 3) Win ActNoGo response:
    if sum(Win_ActNoGo_Trials) == 0
        tWin_ActNoGo = [0 0 0];
        emptyRegressorMatrix(3) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tWin_ActNoGo! \n', iSub)
    else
        tWin_ActNoGo = [tCueAll(Win_ActNoGo_Trials  == 1) ones(sum(Win_ActNoGo_Trials), 2)];
    end
    
    % 4) Avoid ActNoGo response:
    if sum(Avoid_ActNoGo_Trials) == 0
        tAvoid_ActNoGo = [0 0 0];
        emptyRegressorMatrix(4) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tAvoid_ActNoGo! \n', iSub)
    else
        tAvoid_ActNoGo = [tCueAll(Avoid_ActNoGo_Trials  == 1) ones(sum(Avoid_ActNoGo_Trials), 2)];
    end

    % Check for NAs:
    if sum(isnan(tWin_ActGo(:, 3))) > 0; fprintf("Sub%d tWin_ActGo has NAs!", iSub); end
    if sum(isnan(tAvoid_ActGo(:, 3))) > 0; fprintf("Sub%d tAvoid_ActGo has NAs!", iSub); end
    if sum(isnan(tWin_ActNoGo(:, 3))) > 0; fprintf("Sub%d tWin_ActNoGo has NAs!", iSub); end
    if sum(isnan(tAvoid_ActNoGo(:, 3))) > 0; fprintf("Sub%d tAvoid_ActNoGo has NAs!", iSub); end
    
    % 5) Motor response. Note that we take cue onset as timing here: Regressor 5
    tHand       = [tCueAll ones(length(tCueAll), 1) motorNumAll]; % sum command checks how many incorrect responses

    % 6) Error trials. Note that we take cue onset as timing here.
    if sum(ErrorTrials == 1) == 0 % if in fact no errors        
        tError = [0 0 0];
        fprintf('Note: Subject %d has empty regressor Error! \n', iSub)
        emptyRegressorMatrix(6) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
    else
        tError  = [tCueAll(ErrorTrials) ones(sum(ErrorTrials), 2)]; % sum command checks how many incorrect responses
    end
        
    % 7) Outcome onset.
    tOutcomeOnset       = [tFBAll(OutcomeTrials) ones(sum(OutcomeTrials), 2)];

    % 8-9) Update:
    if sum(isnan(Updatewithoutbias(OutcomeTrials))) > 0
        tUpdatestd     = [0 0 0];
        emptyRegressorMatrix(8) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d Block %d has empty regressor tUpdatestd (8)!\n', iSub, iBlock)
    else
        tUpdatestd          = [tFBAll(OutcomeTrials) ones(sum(OutcomeTrials), 1) Updatewithoutbias(OutcomeTrials)];
    end        
    tUpdatestd(:, 3)     = tUpdatestd(:, 3) - mean(tUpdatestd(:, 3)); % demean

    if sum(isnan(Updatedif(OutcomeTrials))) > 0
        tUpdatedif         = [0 0 0];
        emptyRegressorMatrix(9) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d Block %d has empty regressor tUpdatedif (9)!\n', iSub, iBlock)
    else
        tUpdatedif         = [tFBAll(OutcomeTrials) ones(sum(OutcomeTrials), 1) Updatedif(OutcomeTrials)];
    end        
    tUpdatedif(:, 3)    = tUpdatedif(:, 3) - mean(tUpdatedif(:, 3)); % demean

    % 10) Invalid trials. Uninstructed key pressed, thus no outcome, but error message.
    if sum(InvalidTrials) == 0 % if in fact no uninstructed key presses
        tInvalid = [0 0 0];
        fprintf('Note: Subject %d has empty regressor Invalid! \n', iSub)
        emptyRegressorMatrix(10) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
    else
        tInvalid = [tFBAll(InvalidTrials) ones(sum(InvalidTrials), 2)]; % sum command checks how many incorrect responses
    end
    
    % 11) EEG Predictor:
    Trial_EEG   = load(fullfile(dirs.singleTrial, sprintf('%s/TF_singleTrial_%s_sub%03d.csv',EEGPredName,EEGPredName, iSub)));
    if round(mean(Trial_EEG),5)~=0; error('EEG predictor not centered!'); end
    tEEGPred     =  [tFBAll ones(length(tFBAll), 1) Trial_EEG];
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Save regressor matrices:
    
    % Onsets:
    save(fullfile(dirs.timings, 'tWin_ActGo.txt'), 'tWin_ActGo', '-ascii');
    save(fullfile(dirs.timings, 'tWin_ActNoGo.txt'), 'tWin_ActNoGo', '-ascii');
    save(fullfile(dirs.timings, 'tAvoid_ActGo.txt'), 'tAvoid_ActGo', '-ascii');
    save(fullfile(dirs.timings, 'tAvoid_ActNoGo.txt'), 'tAvoid_ActNoGo', '-ascii');

    % Handedness, errors, and outcomes:
    save(fullfile(dirs.timings, 'tHand.txt'), 'tHand', '-ascii');
    save(fullfile(dirs.timings, 'tError.txt'), 'tError', '-ascii');
    save(fullfile(dirs.timings, 'tOutcomeOnset.txt'), 'tOutcomeOnset', '-ascii');
    save(fullfile(dirs.timings, 'tUpdatestd.txt'), 'tUpdatestd', '-ascii');
    save(fullfile(dirs.timings, 'tUpdatedif.txt'), 'tUpdatedif', '-ascii');
    save(fullfile(dirs.timings, 'tInvalid.txt'), 'tInvalid', '-ascii');

    % EEG:
    save(fullfile(dirs.timings, ['t' EEGPredName '.txt']), 'tEEGPred', '-ascii');

    % Save empty regressors:
    csvwrite(fullfile(dirs.timings, 'emptyregressors.txt'),emptyRegressorMatrix);
    
end % end iSub-loop.

% Check whether any regressor empty per file:
% emptyRegressorMatrix
% Starting with subject 1
% Note: Subject 1 has empty regressor Invalid! 
% Starting with subject 2
% Note: Subject 2 has empty regressor Invalid! 
% Starting with subject 3
% Starting with subject 4
% Note: Subject 4 has empty regressor Invalid! 
% Starting with subject 5
% Note: Subject 5 has empty regressor Invalid! 
% Starting with subject 6
% Note: Subject 6 has empty regressor Invalid! 
% Starting with subject 7
% Note: Subject 7 has empty regressor Invalid! 
% Starting with subject 8
% Starting with subject 9
% Note: Subject 9 has empty regressor Invalid! 
% Starting with subject 10
% Starting with subject 11
% Starting with subject 12
% Note: Subject 12 has empty regressor Invalid! 
% Starting with subject 13
% Note: Subject 13 has empty regressor Invalid! 
% Starting with subject 14
% Note: Subject 14 has empty regressor Invalid! 
% Starting with subject 15
% Note: Subject 15 has empty regressor Invalid! 
% Starting with subject 16
% Starting with subject 17
% Note: Subject 17 has empty regressor Invalid! 
% Starting with subject 18
% Note: Subject 18 has empty regressor Invalid! 
% Starting with subject 19
% Note: Subject 19 has empty regressor Invalid! 
% Starting with subject 20
% Note: Subject 20 has empty regressor Invalid! 
% Starting with subject 21
% Note: Subject 21 has empty regressor Invalid! 
% Starting with subject 22
% Note: Subject 22 has empty regressor Invalid! 
% Starting with subject 23
% Note: Subject 23 has empty regressor Invalid! 
% Starting with subject 24
% Note: Subject 24 has empty regressor Invalid! 
% Starting with subject 25
% Starting with subject 26
% Starting with subject 27
% Note: Subject 27 has empty regressor Invalid! 
% Starting with subject 28
% Note: Subject 28 has empty regressor Invalid! 
% Starting with subject 29
% Starting with subject 30
% Starting with subject 31
% Starting with subject 32
% Note: Subject 32 has empty regressor Invalid! 
% Starting with subject 33
% Note: Subject 33 has empty regressor Invalid! 
% Starting with subject 34
% Note: Subject 34 has empty regressor Invalid! 
% Starting with subject 35
% Starting with subject 36
% Note: Subject 36 has empty regressor Invalid! 
% Go back to directory where we came from

end % end of function.
