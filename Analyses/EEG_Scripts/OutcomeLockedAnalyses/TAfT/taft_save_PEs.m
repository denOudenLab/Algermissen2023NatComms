function out = taft_save_PEs(iSub, saveToDisk)

% taft_save_PEs(iSub)
%
% Load behavioral data, compute standard and biased PEs, save as .csv
% files.
% 
% INPUTS:
% iSub              = scalar integer, ID of subject to load.
% saveToDisk        = Boolean, save to disk under dirs.target (optional;
% default: true).
%
% OUTPUTS:
% .csv file with standard/ biased/ difference PEs and Updates saved to disk
% under dirs.root/Log/Behavior/Behavior_PEs.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% for iSub = 1:36
%     taft_save_PEs(iSub)
% end

% clear all
% iSub = 3;

% ----------------------------------------------------------------------- %
%% Complete settings:

fprintf('Start subject %03d\n', iSub);


if nargin < 2
    saveToDisk = true;
end

% ----------------------------------------------------------------------- %
%% Determine behavioral file name per subject and load:

% Root directory:
dirs.root   	= taft_set_rootDir(); % /project/3017042.02

% Input directory:
behavFile = fullfile(dirs.root, 'Log/Behavior/Data_beh_mat', ...
    sprintf('3017042.02_emmvdij_%03d_001_results.mat', iSub));    

% Load behavior:
load(behavFile)

% Output directory:
dirs.target     = fullfile(dirs.root, 'Log/Behavior/Behavior_PEs');
if ~exist(dirs.target, 'dir'); mkdir(dirs.target); end

% ----------------------------------------------------------------------- %
%% Combine data from both sessions (blocks 1-3 and 4-6) into one vector:

seq.stim        = [prep.seq.learn.stim{1}; prep.seq.learn.stim{2}];
seq.resp        = [results.learn{1}.response; results.learn{2}.response];
seq.RT          = [results.learn{1}.RT; results.learn{2}.RT];
seq.outcome     = [results.learn{1}.outcome; results.learn{2}.outcome];
seq.accuracy    = [results.learn{1}.acc; results.learn{2}.acc];
seq.go          = [results.learn{1}.go; results.learn{2}.go];
seq.splithalf   = [ones(160, 1); 2*ones(160, 1);ones(160, 1); 2*ones(160, 1)];

% ----------------------------------------------------------------------- %
%% Carry forward stimulus properties: 

clear out % initialize out

out.stim        = seq.stim;
out.goCue       = 1:4;
out.reqaction   = ~ismember(out.stim,out.goCue) + 1; % 1 if yes, 2 if no
out.winCue      = [1 2 5 6];
out.valence     = ~ismember(out.stim,out.winCue) + 1; % 1 if yes, 2 if no
out.iswin       = ismember(out.stim,out.winCue); % 1 is Win, 0 if Avoid

% ----------------------------------------------------------------------- %
%% Initialize objects to recode: 

out.reqresp     = nan(size(out.stim)); % initialize in same format 
out.resp        = nan(size(out.stim)); % initialize in same format 
out.reqgo       = nan(size(out.stim)); % initialize in same format 
out.isgo        = nan(size(out.stim)); % initialize in same format 
out.isnextgo    = nan(size(out.stim)); % initialize in same format 
out.accuracy    = nan(size(out.stim)); % initialize in same format 
out.RT          = nan(size(out.stim)); % initialize in same format 
out.c           = nan(size(out.stim)); % initialize in same format 

% Whether trials are valid or not:
out.outcome     = seq.outcome; % outcome
out.valid       = find(~isnan(seq.outcome));

% ----------------------------------------------------------------------- %
%% Loop over trials:

nTrials = length(seq.resp);

for iTrial = 1:nTrials % k = 1

    % ------------------------------------------------------------------- %
    %% a) Recode recorded response (also for uninstructed key presses):
    
    if ismember(seq.resp(iTrial), [69, 70, 71, 72, 101, 102, 103, 104]) % left Go response
            out.resp(iTrial) = 101; % actually instructed left Go response is 101
            out.isgo(iTrial) = 1; % any Go
            out.c(iTrial)    = 1; % recode left Go numerically to 1
            
    elseif ismember(seq.resp(iTrial), [65, 66, 67, 68, 97, 98, 99, 100]) % right Go response
            out.resp(iTrial) = 97; % actually instructed right Go response is 97
            out.isgo(iTrial) = 1; % any Go
            out.c(iTrial)    = 2; % recode right Go numerically to 2
   
    else % NoGo response
            out.resp(iTrial) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
            out.isgo(iTrial) = 0; % any NoGo
            out.c(iTrial)    = 3; % recode NoGo numerically to 3
    end
    
    % ------------------------------------------------------------------- %
    %% b) Recode required response (based on stimulus number):
    
    if ismember(out.stim(iTrial), [1, 3]) % left Go response
            out.reqresp(iTrial) = 101; % actually instructed left Go response is 101
            out.reqgo(iTrial)   = 1;
            out.RT(iTrial)	= seq.RT(iTrial);
            
    elseif ismember(out.stim(iTrial), [2, 4]) % right Go response
            out.reqresp(iTrial) = 97; % actually instructed right Go response is 97
            out.reqgo(iTrial)   = 1;
            out.RT(iTrial) 	= seq.RT(iTrial);
            
    else % NoGo response
            out.reqresp(iTrial) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
            out.reqgo(iTrial)   = 0;
            out.RT(iTrial) 	= NaN;
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
        out.isgo(iTrial)    = 0; % 
        
    else
        out.RT(iTrial)      = round(seq.RT(iTrial),3); % round to
        % ms be consistent with EEG sampling rate (in ms)
    end
    
    % ------------------------------------------------------------------- %
    %% d) Recode accuracy based on recorded responses:
    
    if out.reqresp(iTrial) ~= out.resp(iTrial) || seq.RT(iTrial) > 1.3035 % wrong response (key) or after time
        out.accuracy(iTrial) = 0; % 0 if incorrect
        
    else
        out.accuracy(iTrial) = 1; % 1 if correct
    end
    
end % end of iTrial loop.

% ----------------------------------------------------------------------- %
%% Executed action x valence:

out.isW2G   = out.iswin == 1 & out.isgo == 1;
out.isW2NG  = out.iswin == 1 & out.isgo == 0;
out.isA2G   = out.iswin == 0 & out.isgo == 1;
out.isA2NG  = out.iswin == 0 & out.isgo == 0;

% ----------------------------------------------------------------------- %
%% Outcome:

out.fb.abs                  = out.outcome;
out.fb.rel                  = 1 - out.outcome; % 1 becomes good, 2 becomes bad
out.fb.rel(out.valence==1)  = out.fb.rel(out.valence == 1) + 1; % Win cues: one up (0 becomes 1, 1 becomes 2)
out.fb.all                  = out.fb.rel + 2*(out.valence == 2); % 1-4: reward, no-reward, no-punishment, punishment
% [out.fb.rel out.fb.abs out.fb.all]

out.fbrel                   = 2 - out.fb.rel; % 1 = reward, 0 = punishment
out.salience                = ismember(out.fb.all, [1 4]); % 1 = salient, 0 = not salient

% ----------------------------------------------------------------------- %
%% Executed action x outcome valence:

out.isRew       = out.fb.all == 1;
out.isNeu       = ismember(out.fb.all, [2 3]);
out.isPun       = out.fb.all == 4;

out.isGoRew     = out.isgo == 1 & out.fb.all == 1;
out.isGoNoRew   = out.isgo == 1 & out.fb.all == 2;
out.isGoNoPun   = out.isgo == 1 & out.fb.all == 3;
out.isGoPun     = out.isgo == 1 & out.fb.all == 4;
out.isNoGoRew   = out.isgo == 0 & out.fb.all == 1;
out.isNoGoNoRew = out.isgo == 0 & out.fb.all == 2;
out.isNoGoNoPun = out.isgo == 0 & out.fb.all == 3;
out.isNoGoPun   = out.isgo == 0 & out.fb.all == 4;
% 'isGoRew', 'isGoNoRew', 'isGoNoPun', 'isGoPun', 'isNoGoRew', 'isNoGoNoRew', 'isNoGoNoPun', 'isNoGoPun'
% sum(out.isW2G + out.isW2NG + out.isA2G + out.isA2NG)

% ----------------------------------------------------------------------- %
%% Load and prepare model parameter:

% Initialize parameters:
modelParameters = [1.0889; -1.7971; 0.1224; 0.5457; 2.4423; 1.5299; 2.4356;
    0.5339; 0.9706; 2.0508];

% Rho and epsilon still untransformed
rho     = exp(modelParameters(1));
epsilon = exp(modelParameters(2)) / (1 + exp(modelParameters(2)));
% gobias and pibias not necessary as subjects' actual choices and
% outcomes used for computations; thus no action weights computed

% Biased epsilon:
if epsilon > 0.5
    epsilon_plus_kappa  = exp(modelParameters(2) + modelParameters(5)) / ...
    (1 + exp(modelParameters(2) + modelParameters(5)));
    epsilon_minus_kappa = epsilon - (epsilon_plus_kappa - epsilon); 
else
   epsilon_minus_kappa  = exp(modelParameters(2)-modelParameters(5)) / ...
   (1+exp(modelParameters(2)-modelParameters(5)));
   epsilon_plus_kappa   = epsilon + (epsilon - epsilon_minus_kappa);         
end

% ----------------------------------------------------------------------- %
%% Simulate model:

% Initialize other task settings:
nStim           = 16;
nResp           = 3;
nTrial          = 640;

% Initialize variables:
% Define initial value of each stimulus (separate matrix for each stimulus):

% States:
stimAll         = [out.stim(1:320); 8 + out.stim(321:640)];     % add 8 to stimulusID in second session of task
PavValues       = repmat([0.5, 0.5, -0.5, -0.5], [1, 4]);       % cue values, relevant for initializing action values
stimcount       = zeros(1, nStim);                              % count number of occurences of stimuli
stimRep         = nan(nTrial, 1);                               % cue valence known or not per trial
cueDetected     = zeros(1, nStim);                              % valenced outcome for cue already experienced or not
cueDetectAll    = nan(nTrial, 1);                               % cue valence known or not per trial

% Actions:
%av0         = zeros(1,nResp,nStim);                                    % initial action values to zero
av0                 = reshape(repelem(rho * PavValues, nResp), ....
    1, nResp, nStim);  % initial action values
Outcome             = nan(nTrial, 1);                                   % initialise outcome
avwithoutbias       = av0;                                              % initialise current action values av without bias
avwithbias          = av0;                                              % initialise current action values av with bias

% Final target objects:
EVwithoutbias       = nan(nTrial, 1);                                   % initialize EV without bias
PEwithoutbias       = nan(nTrial, 1);                                   % initialize PE without bias
Updatewithoutbias   = nan(nTrial, 1);                                   % initialize update without bias
EVwithbias          = nan(nTrial, 1);                                   % initialise EV with bias
PEwithbias          = nan(nTrial, 1);                                   % initialise PE with bias
Updatewithbias      = nan(nTrial, 1);                                   % initialise update with bias

% Loop over trials:
for iTrial = 1:nTrial

    % 1) Retrieve stimulus ID
    stim = stimAll(iTrial);
    % Number of occurrences for this stimulus:
    stimcount(stim) = stimcount(stim) + 1;
    stimRep(iTrial) = stimcount(stim); % save
    % Store whether cue valence experienced in this trial or not
    cueDetectAll(iTrial) = cueDetected(stim);

    % 2) Make choice:
    % Retrieve actually taken action:
    c = out.c(iTrial);
    % Retrieve actually received outcome:
    r = out.outcome(iTrial);

    % 3) Learn based on outcome
    % Determine learning rate considering Pavlovian bias on instrumental learning:
    if ismember(c, [1, 2]) && r == 1
        effepsilon = epsilon_plus_kappa;
    elseif c == 3 && r == -1
        effepsilon = epsilon_minus_kappa;
    else
        effepsilon = epsilon;
    end

    % check whether valid key pressed and outcome ever obtained, otherwise no learning:
    if isnan(r) == 1 % no updating if outcome is NaN, i.e. subject got feedback that uninstructed key was pressed
        
        Outcome(iTrial)             = 0;
        EVwithoutbias(iTrial)       = avwithoutbias(:, c, stim); % carry forward Q-values
        EVwithbias(iTrial)          = avwithbias(:, c, stim); % carry forward Q-values
        PEwithoutbias(iTrial)       = 0;
        PEwithbias(iTrial)          = 0;            
        Updatewithoutbias(iTrial)   = 0;
        Updatewithbias(iTrial)      = 0;    
        
    else % else: normal Rescorla-Wagner model updating
        Outcome(iTrial)             = rho * r; % compute prediction error as reward minus expected value to be updated
 
        % update action values without bias (epsilon):
        EVwithoutbias(iTrial)       = avwithoutbias(:, c, stim);
        PEwithoutbias(iTrial)       = Outcome(iTrial) - EVwithoutbias(iTrial); % compute prediction error as reward minus expected value to be updated
        Updatewithoutbias(iTrial)   = epsilon * PEwithoutbias(iTrial);
        avwithoutbias(:, c, stim)   = EVwithoutbias(iTrial) + Updatewithoutbias(iTrial); % update respective action value

        % update action values with bias (effepsilon):
        EVwithbias(iTrial)          = avwithbias(:, c, stim);
        PEwithbias(iTrial)          = Outcome(iTrial) - EVwithbias(iTrial); % compute prediction error as reward minus expected value to be updated
        Updatewithbias(iTrial)      = effepsilon * PEwithbias(iTrial);
        avwithbias(:, c, stim)      = EVwithbias(iTrial) + Updatewithbias(iTrial); % update respective action value
        
        % "Unmute" stimulus if valence seen for the first time
        if abs(r) == 1; cueDetected(stim) = 1; end
    end
end % end trial loop

% Compute difference between PEs and updates:
EVdif           = EVwithbias - EVwithoutbias; % std + dif = bias
PEdif           = PEwithbias - PEwithoutbias; % std + dif = bias
Updatedif       = Updatewithbias - Updatewithoutbias; % std + dif = bias

% ----------------------------------------------------------------------- %
%% Store in output file:

% stimuli:
out.stimAll     = stimAll;
out.stimRep     = stimRep;

% Outcomes:
out.Outcome     = Outcome;

% Q-values:
out.EVstd       = EVwithoutbias;
out.EVbias      = EVwithbias;
out.EVdif       = EVdif;

% PEs:
out.PEstd       = PEwithoutbias;
out.PEbias      = PEwithbias;
out.PEdif       = PEdif;

% Updates:
out.Updatestd   = Updatewithoutbias;
out.Updatebias  = Updatewithbias;
out.Updatedif   = Updatedif; 

% Absolute prediction errors:
out.Updatestdabs   = abs(Updatewithoutbias);
out.Updatebiasabs  = abs(Updatewithbias);

% ----------------------------------------------------------------------- %
%% Check for next response:

out.isnextgo = nan(nTrial, 1);
for iTrial = 1:nTrial % iTrial = 1
    if ~ismember(iTrial, [110 220 320 430 540 640]) % if not end of block (leave initialized NA)
        out.isnextgo(iTrial) = out.isgo(iTrial+1);
    end
end

% ----------------------------------------------------------------------- %
%% Save:

if saveToDisk
    
    fprintf('Save PEs and Updates under %s\n', dirs.target);
    outputFile = fullfile(dirs.target, sprintf('PEstd/PEstd_sub%03d.csv', iSub)); 
    csvwrite(outputFile, out.PEstd)
    outputFile = fullfile(dirs.target, sprintf('PEbias/PEbias_sub%03d.csv', iSub)); 
    csvwrite(outputFile, out.PEbias)
    outputFile = fullfile(dirs.target, sprintf('PEdif/PEdif_sub%03d.csv', iSub)); 
    csvwrite(outputFile, out.PEdif)

    outputFile = fullfile(dirs.target, sprintf('Updatestd/Updatestd_sub%03d.csv', iSub)); 
    csvwrite(outputFile, out.Updatestd)
    outputFile = fullfile(dirs.target, sprintf('Updatebias/Updatebias_sub%03d.csv', iSub)); 
    csvwrite(outputFile, out.Updatebias)
    outputFile = fullfile(dirs.target, sprintf('Updatedif/Updatedif_sub%03d.csv', iSub)); 
    csvwrite(outputFile, out.Updatedif)
    
end

fprintf('Finished subject %03d\n', iSub);
% clear all

end % END OF FUNCTION.