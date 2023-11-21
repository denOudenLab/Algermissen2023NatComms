function EEGfMRIPav_GLM1B()

% EEGfMRIPav_GLM1B
% 
% Create regressors for block-level GLM.
% Cue x action conditions, response, error, outcome onset, update and
% difference to biased update, invalid responses.
% Use three-column (onset, duration, factor level) format.
% Saves each regressor as separate .txt file in respective
% timings_regressor folder of respective GLM of respective subject.
% 
% Use group-level HBI parameters from M8 (with separate perseveration
% parameters for Win and Avoid cues).
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
% Mind setting root directory dirs.root
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% clear all; clc

% ----------------------------------------------------------------------- %
%% Set directories:

rootDir     = '/project/3017042.02'; % root directory--needs to be adapted to users' folder structure

% This is where we are
% cd /project/3017042.02/Analyses/fMRI_Scripts/Create_Regressors/Create_Task_Regressors

% Settings:
GLMID       = '1B';
nSubjects   = 36;
nBlocks     = 6;
nRegressors = 10; % for specification of empty regressors

% Directory where raw data are located:
dirs.behav  = fullfile(rootDir, 'Log/Behavior');
dirs.raw    = fullfile(dirs.behav, 'Data_beh_mat');

% ----------------------------------------------------------------------- %
%% Detect which files are in behavioral directory:
 
cd(dirs.raw); % go to folder with result files of all subjects.
input           = dir; % store details of all files there
allfiles        = {input.name}; % extract names and store as cell array
allfiles        = allfiles'; % tranpose
allfiles{1, 1}  = []; % delete '.'
allfiles{2, 1}  = []; % delete '..'   
dirs.subs       = allfiles(~cellfun('isempty', allfiles)); % remove deleted entries
% --> array with all raw data files

% ----------------------------------------------------------------------- %
%% Initialize group-level model parameters:

% Group-level HBI parameters:
% dirs.cbm    = fullfile(dirs.behav, 'Modelling_CBM');
% tmp         = load(fullfile(dirs.cbm, '/HBI_Results/hbi_mod_8.mat'));
% modHBI      = tmp.cbm.output.group_mean{:}; % group-level

% 5 means:
modelParameters = [1.2104   -2.5158    0.0970    0.2287    2.1696    1.8973    0.7661];

% Initialize parameters:

% Rho and epsilon still untransformed:
rho     = exp(modelParameters(1));
epsilon = exp(modelParameters(2)) / (1 + exp(modelParameters(2)));

% gobias and pibias not necessary as subjects' actual choices and
% outcomes used for computations; thus no action weights computed

% Biased epsilon:
if epsilon > 0.5
    epsilon_plus_kappa  = exp(modelParameters(2) + modelParameters(5)) / (1 + exp(modelParameters(2) + modelParameters(5)));
    epsilon_minus_kappa = epsilon - (epsilon_plus_kappa - epsilon); 
else
   epsilon_minus_kappa  = exp(modelParameters(2) - modelParameters(5)) / (1 + exp(modelParameters(2) - modelParameters(5)));
   epsilon_plus_kappa   = epsilon + (epsilon - epsilon_minus_kappa);         
end

% ----------------------------------------------------------------------- %
%% MR regressors per block:

for iSub = 1:nSubjects % length(dirs.subs) % specify for which subjects to create regressors
    
    fprintf('Starting with subject %d\n', iSub)
 
    % Load event and timing data:
    filename = strcat('3017042.02_emmvdij_', sprintf('%03d', iSub), '_001_results.mat'); % raw data of specific subject (mind that subject ID can be 2 or 3 digits)
    cd(dirs.raw); % go to folder with result files of all subjects.
    load(filename);
    
    % Set directory where regressors will be printed:
    dirs.timings = fullfile(rootDir, 'Log', 'fMRI',...
                sprintf('sub-%03d', iSub), sprintf('GLM%s', GLMID), 'timings_regressors');
    
    % Retrieve data for regressors.
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
                respAll(k)      = prep.par.key.left; % actually instructed left Go response is 101
                respNumAll(k)   = 1; % recode left Go numerically to 1
                motorNumAll(k)  = 1; % left hand
        elseif ismember(respAll(k), [65, 66, 67, 68, 97, 98, 99, 100]) % right Go response
                respAll(k)      = prep.par.key.right; % actually instructed right Go response is 97
                respNumAll(k)   = 2; % recode right Go numerically to 2
                motorNumAll(k)  = -1; % right hand
        else % NoGo response
                respAll(k)      = 0; % no response = NoGo; implicitly also sets NAs to NoGo
                respNumAll(k)   = 3; % recode NoGo numerically to 3
                motorNumAll(k)  = 0; % no hand
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
    stimAll     = [stimIDAll(1:320);8+stimIDAll(321:640)];                  % add 8 to stimulusID in second session of task
    PavValues   = repmat([0.5,0.5,-0.5,-0.5], [1,4]);                       % cue values, relevant for initializing action values
    stimcount   = zeros(1, nStim);                                          % count number of occurences of stimuli
    cueDetected = zeros(1, nStim);                                          % valenced outcome for cue already experienced or not
    cueDetectAll= nan(nTrial, 1);                                           % cue valence known or not per trial

    % Actions:
    %av0         = zeros(1, nResp, nStim);                                  % initial action values to zero
    av0         = reshape(repelem(rho*PavValues, nResp), 1, nResp, nStim);  % initial action values
    avwithbias          = av0;                                              % initialise current action values av with bias
    avwithoutbias       = av0;                                              % initialise current action values av without bias
    PEwithbias          = nan(nTrial, 1);                                   % initialise PE with bias
    Updatewithbias      = nan(nTrial, 1);                                   % initialise PE with bias
    PEwithoutbias       = nan(nTrial, 1);                                   % initialize PE without bias
    Updatewithoutbias   = nan(nTrial, 1);                                   % initialize PE without bias
    validOutcome        = nan(nTrial, 1);                                   % store trials where outcome was NaN to later delete row in regressor

    % Loop over trials, simulate action values and PEs:
    for iTrial = 1:nTrial
        
        % 1) Retrieve stimulus ID:
        stim = stimAll(iTrial);
        % Number of occurrences for this stimulus:
        stimcount(stim) = stimcount(stim)+1;
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
        elseif c == 3 && r ==  -1
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
            
        else % normal Rescorla-Wagner model updating
            
            validOutcome(iTrial) = 1;
            
            % Update action values without bias (epsilon):
            PEwithoutbias(iTrial)       = rho * r - avwithoutbias(:, c, stim); % compute prediction error as reward minus expected value to be updated
            Updatewithoutbias(iTrial)   = epsilon * PEwithoutbias(iTrial);
            avwithoutbias(:, c, stim)   = avwithoutbias(:, c, stim) + Updatewithoutbias(iTrial); % update respective action value
            
            % Update action values with bias (effepsilon):
            PEwithbias(iTrial)          = rho * r - avwithbias(:, c, stim); % compute prediction error as reward minus expected value to be updated
            Updatewithbias(iTrial)      = effepsilon * PEwithbias(iTrial);
            avwithbias(:, c, stim)      = avwithbias(:, c, stim) + Updatewithbias(iTrial); % update respective action value

            % "Unmute" stimulus if valence seen for the first time
            if abs(r) == 1; cueDetected(stim) = 1; end
        end
    end % end trial loop
    
    % Compute difference between PEs and updates:
    Updatedif   = Updatewithbias - Updatewithoutbias;

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Create timing for regressors:
    
    % a) Extract when blocks started in behavioral time:
    % Break at trial 110 and 220, where trial 1-110 are block one, 111-220 block two, i.e.
    % Retrieve timing of start of first 3 blocks, repeate for each of the 110 trials:
    t0part1 = reshape(repmat(tm.learn{1}.stim([1 111 221]), 1, 110)', [330 1]); % one long vector with 330 elements
    t0part1 = t0part1(1:320); % drop last 10 entries (actually first 10, see below)
    % Retrieve timing of start of second 3 blocks, repeate for each of the 110 trials:
    t0part2 = reshape(repmat(tm.learn{2}.stim([1 111 221]), 1, 110)', [330 1]);
    t0part2 = t0part2(1:320); % drop last 10 entries (actually first 10, see below)
    t0      = [t0part1; t0part2] - 10; % concatenate; now subtract 10 seconds that were used for getting MRI signal steady-state 
    
    % b) Extract trial timing in behavioral time, correct for when block started:
    % 1) retrieve timings of events corrected for t0 (onset of block):
    tCueAll        = [tm.learn{1}.stim; tm.learn{2}.stim] - t0; % take cue onset; subtract start of each block
%     tRespAll       = [tm.learn{1}.response; tm.learn{2}.response] - t0; % note that in case of nogo, the cue offset time was recorded.
   
    % 2) Retrieve timing of learning, correct for NaNs when people got error message:
    tFBPrel        = [tm.learn{1}.outcome; tm.learn{2}.outcome]; % take onset of learning, so far uncorrected
    tWarningAll    = [tm.learn{1}.warning; tm.learn{2}.warning]; % take timing of warning messages
    tFBPrel(isnan(tFBPrel)) = tWarningAll; % replace NaNs with timing of warnings
    tFBAll         = tFBPrel - t0; % subtract start of each block
    
    % 3) create timingfiles per block (array with 6 entries)
    trials2use = {1:110, 111:220, 221:320, 321:430,431:540,541:640};
    

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Extract behavioral data for each block:

    for iBlock = 1:nBlocks % iBlock = 1
        
        % Indices of trials in this block:
        idx         = trials2use{iBlock}; % retrieve trials for this block
        nBlockTrials= length(idx);
      
        % Behavioral variables in this block:
        stimID      = stimIDAll(idx); % extract cue conditions (1-8)
        resp        = respNumAll(idx); % extract responses (recoded above: 1, 2 , 3)
        motor       = motorNumAll(idx); % extract responses (recoded above: 1, 2 , 3)
        
        incorrect   = incorrectAll(idx); % extract accuracy (1 for incorrect, 0 for correct)
        valid       = validOutcome(idx); % extract in which trials people pressed an uninstructed key and thus got no outcome, but an error message

        cuedetect   = cueDetectAll(idx); % regressor for cue only if valence known
        
        % Saving directory for the current block: 
        dirs.timings = fullfile(rootDir, 'Log', 'fMRI',...
                sprintf('sub-%03d', iSub), sprintf('GLM%s', GLMID), sprintf('timings_block%d', iBlock));
        % Create directory if it doesn't exist yet: 
        if ~exist(dirs.timings , 'dir'); mkdir(dirs.timings); end

        % --------------------------------------------------------------- %
        %% Identify which trial has which event (indices):

        % 1) Cues:
        % Split up trials into those with Win vs. Avoid cues and Go vs. NoGo responses:
        % mind cueDetectAll: don't model cues for trials where outcome not yet experienced
        WinGoTrials     = (ismember(stimID, [1 2 5 6]) & ismember(resp, [1 2]) & cuedetect == 1); % indicator whether trial was Win trial with Go response
        AvoidGoTrials   = (ismember(stimID, [3 4 7 8]) & ismember(resp, [1 2]) & cuedetect == 1); % indicator whether trial was Win trial with Go response
        WinNoGoTrials   = (ismember(stimID, [1 2 5 6]) & resp == 3 & cuedetect == 1); % indicator whether trial was Win trial with NoGo response
        AvoidNoGoTrials = (ismember(stimID, [3 4 7 8]) & resp == 3 & cuedetect == 1); % indicator whether trial was Avoid trial with NoGo response

        % 2) Responses:
        ErrorTrials     = incorrect == 1;
        
        % 3) Feedback:
        InvalidTrials   = valid == 0;  
        
        emptyRegressorMatrix = zeros(1, nRegressors); % initialize empty regressors (default: non-empty, so zero)

        % --------------------------------------------------------------- %
        % --------------------------------------------------------------- %
        % --------------------------------------------------------------- %
        %%  Create three-column regressors:
        % 3-columns format:
        % Onset, Duration (1), Value (1 for factorial regressors)

        % CUE AND RESPONSES:
        % 1 Cue condition and actual action:
        
        % 1a) Win Go response: Regressor 1
        if sum(WinGoTrials) == 0
            tWinGo            = [0 0 0];
            emptyRegressorMatrix(1) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
            fprintf('Note: Subject %d Block %d has empty regressor tWinGo (1)!\n', iSub, iBlock)
        else
            tWinGo            = [tCueAll(idx(WinGoTrials)) ones(sum(WinGoTrials), 2)];
        end
        
        % 1b) Avoid Go response: Regressor 2
        if sum(AvoidGoTrials) == 0
            tAvoidGo            = [0 0 0];
            emptyRegressorMatrix(2) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
            fprintf('Note: Subject %d Block %d has empty regressor tAvoidGo (2)! \n', iSub, iBlock)
        else
            tAvoidGo            = [tCueAll(idx(AvoidGoTrials)) ones(sum(AvoidGoTrials), 2)];
        end
        
        % 1c) Win NoGo response: Regressor 3
        if sum(WinNoGoTrials) == 0
            tWinNoGo            = [0 0 0];
            emptyRegressorMatrix(3) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
            fprintf('Note: Subject %d Block %d has empty regressors tWinNoGo (3)! \n', iSub, iBlock)
        else
            tWinNoGo            = [tCueAll(idx(WinNoGoTrials)) ones(sum(WinNoGoTrials), 2)];
        end 

        % 1d) Avoid NoGo response: Regressor 4
        if sum(AvoidNoGoTrials) == 0
            tAvoidNoGo          = [0 0 0];
            emptyRegressorMatrix(4) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
            fprintf('Note: Subject %d Block %d has empty regressor tAvoidNoGo (4)! \n', iSub, iBlock)
        else
            tAvoidNoGo          = [tCueAll(idx(AvoidNoGoTrials)) ones(sum(AvoidNoGoTrials), 2)];
        end
        
        % --------------------------------------------------------------- %
        % 2) RESPONSE
        % 2a) Motor response. Note that we take cue onset as timing here: Regressor 5
        tMotor      = [tCueAll(idx) ones(sum(nBlockTrials), 1) motor]; % sum command checks how many incorrect responses

        % 2b) Error trials. Note that we take cue onset as timing here: Regressor 6
        if sum(ErrorTrials == 1) == 0 % if in fact no errors        
            tError = [0 0 0];
            emptyRegressorMatrix(6) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
            fprintf('Note: Subject %d Block %d has empty regressor Error (6)! \n', iSub, iBlock)
        else
            tError  = [tCueAll(idx(ErrorTrials)) ones(sum(ErrorTrials), 2)]; % sum command checks how many incorrect responses
        end
        
        % --------------------------------------------------------------- %
        % 3) OUTCOME:

        % 3a) Outcome onset
        tOutcome                = [tFBAll(idx) ones(nBlockTrials, 2)];

        % 3b) Update of standard RW model:
        if sum(isnan(Updatewithoutbias(idx))) > 0
            tUpdatestd     = [0 0 0];
            emptyRegressorMatrix(8) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
            fprintf('Note: Subject %d Block %d has empty regressor tUpdatestd (8)!\n', iSub, iBlock)
        else
            tUpdatestd          = [tFBAll(idx) ones(nBlockTrials, 1) Updatewithoutbias(idx)];
        end        
                
        % 3c) Update difference between biased and standard RW model:
        if sum(isnan(Updatedif(idx))) > 0
            tUpdatedif         = [0 0 0];
            emptyRegressorMatrix(9) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
            fprintf('Note: Subject %d Block %d has empty regressor tUpdatedif (9)!\n', iSub, iBlock)
        else
            tUpdatedif         = [tFBAll(idx) ones(nBlockTrials, 1) Updatedif(idx)];
        end        

        % --------------------------------------------------------------- %
        % 4) Invalid trials. Uninstructed key pressed, thus no outcome, but error message: Regressor 16
        if sum(InvalidTrials) == 0 % if in fact no uninstructed key presses
            tInvalid = [0 0 0];
            emptyRegressorMatrix(10) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
            fprintf('Note: Subject %d Block %d has empty regressor Invalid (10)! \n', iSub, iBlock)
        else
            tInvalid = [tFBAll(idx(InvalidTrials)) ones(sum(InvalidTrials), 2)]; % sum command checks how many incorrect responses
%            fprintf('Note: Subject %d Block %d has NONEMPTY regressor Invalid! \n', iSub, iBlock)
        end

        % --------------------------------------------------------------- %
        % Mean correct parametric modulators:
        tUpdatestd(:, 3)    = tUpdatestd(:, 3) - mean(tUpdatestd(:, 3));
        tUpdatedif(:, 3)    = tUpdatedif(:, 3) - mean(tUpdatedif(:, 3));

        % --------------------------------------------------------------- %
        % --------------------------------------------------------------- %
        % --------------------------------------------------------------- %
        %% Save regressor matrices:
        
        % Onsets:
        save(fullfile(dirs.timings, 'tWinGo.txt'), 'tWinGo', '-ascii');
        save(fullfile(dirs.timings, 'tAvoidGo.txt'), 'tAvoidGo', '-ascii');
        save(fullfile(dirs.timings, 'tWinNoGo.txt'), 'tWinNoGo', '-ascii');
        save(fullfile(dirs.timings, 'tAvoidNoGo.txt'), 'tAvoidNoGo', '-ascii');
        
        % Responses:
        save(fullfile(dirs.timings, 'tMotor.txt'), 'tMotor', '-ascii');
        save(fullfile(dirs.timings, 'tError.txt'), 'tError', '-ascii');
        
        % Outcomes:
        save(fullfile(dirs.timings, 'tOutcome.txt'), 'tOutcome', '-ascii');
        save(fullfile(dirs.timings, 'tUpdatestd.txt'), 'tUpdatestd', '-ascii');
        save(fullfile(dirs.timings, 'tUpdatedif.txt'), 'tUpdatedif', '-ascii');
        save(fullfile(dirs.timings, 'tInvalid.txt'), 'tInvalid', '-ascii');
        
        % Save empty regressors:
        csvwrite(fullfile(dirs.timings, 'emptyregressors.txt'), emptyRegressorMatrix);
        
    end % end iBlock-loop.
end % end iSub-loop.

end % end of function.
