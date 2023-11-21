function taft_runonce_select_trials()

% taft_runonce_select_trials()
%
% Load behavior, retrieves indices of trials on which certain relevant
% behaviors/ conditions happen (Go/ NoGo, Win/ Avoid, congruent/
% incongruent), save indices as separate files under TAfT_Betas/selTrials.
% Mind setting the root directory in taft_set_rootDir().
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT/

% ----------------------------------------------------------------------- %
%% Settings:

nSub            = 36; % number subjects

% ----------------------------------------------------------------------- %
%% Directories:

dirs.root       = taft_set_rootDir(); % /project/3017042.02

dirs.target     = fullfile(dirs.root, 'Log/EEG/OutcomeLockedResults/TAfT_Betas/selTrials');
% Create target directory if not existing yet:
if ~exist(dirs.target, 'dir'); mkdir(dirs.target); end

% ----------------------------------------------------------------------- %
%% Loop through subjects:

for iSub = 1:nSub
    
    fprintf('Start subject %03d\n', iSub);
    
    % ------------------------------------------------------------------- %
    %% Load behavior:    
    
    job.behavFile   = fullfile(dirs.root, '/Log/Behavior/Data_beh_mat', ...
    sprintf('3017042.02_emmvdij_%03d_001_results.mat', iSub));    

    out             = taft_preprocess_load_behavior(job);    
    
    % ------------------------------------------------------------------- %
    %% Select indices of where relevant behavior/ condition occurs, save:
    
    % Go trials:
    fprintf('Go actions\n')
    selTrials   = find(out.isgo == 1);
    save(fullfile(dirs.target, sprintf('selTrials_GoTrials_Sub%03d.mat', iSub)), 'selTrials');
    
    % NoGo trials:
    fprintf('NoGo actions\n')
    selTrials   = find(out.isgo == 0);
    save(fullfile(dirs.target, sprintf('selTrials_NoGoTrials_Sub%03d.mat', iSub)), 'selTrials');
    
    % Reward trials:
    fprintf('Reward outcomes\n')
    selTrials   = find(out.fb.abs == 1);
    save(fullfile(dirs.target, sprintf('selTrials_RewardTrials_Sub%03d.mat', iSub)), 'selTrials');

    fprintf('No Reward outcomes\n')
    selTrials   = find(out.fb.all == 2);
    save(fullfile(dirs.target, sprintf('selTrials_NoRewardTrials_Sub%03d.mat', iSub)), 'selTrials');

    % Neutral trials:
    fprintf('Neutral outcomes\n')
    selTrials   = find(out.fb.abs == 0);
    save(fullfile(dirs.target, sprintf('selTrials_NeutralTrials_Sub%03d.mat', iSub)), 'selTrials');
    
    % Punishment trials:
    fprintf('Punishment outcomes\n')
    selTrials   = find(out.fb.abs == -1);
    save(fullfile(dirs.target, sprintf('selTrials_PunishmentTrials_Sub%03d.mat', iSub)), 'selTrials');

    % No Punishemnt trials:
    fprintf('NoPunishment outcomes\n')
    selTrials   = find(out.fb.all == 3);
    save(fullfile(dirs.target, sprintf('selTrials_NoPunishmentTrials_Sub%03d.mat', iSub)), 'selTrials');

    % Positive trials:
    fprintf('Positive outcomes\n')
    selTrials   = find(out.fb.rel == 1);
    save(fullfile(dirs.target, sprintf('selTrials_PositiveTrials_Sub%03d.mat', iSub)), 'selTrials');
    
    % Negative trials:
    fprintf('Negative outcomes\n')
    selTrials   = find(out.fb.rel == 2);
    save(fullfile(dirs.target, sprintf('selTrials_NegativeTrials_Sub%03d.mat', iSub)), 'selTrials');

    % Go trials with rewards:
    fprintf('Go actions with reward\n');
    selTrials   = find(out.isgo == 1 & out.fb.all == 1);
    save(fullfile(dirs.target, sprintf('selTrials_GoRewTrials_Sub%03d.mat', iSub)), 'selTrials');

    % Go trials with no rewards:
    fprintf('Go actions with no reward\n');
    selTrials   = find(out.isgo == 1 & out.fb.all == 2);
    save(fullfile(dirs.target, sprintf('selTrials_GoNoRewTrials_Sub%03d.mat', iSub)), 'selTrials');

    % Go trials with no punishments:
    fprintf('Go actions with no punishment\n');
    selTrials   = find(out.isgo == 1 & out.fb.all == 3);
    save(fullfile(dirs.target, sprintf('selTrials_GoNoPunTrials_Sub%03d.mat', iSub)), 'selTrials');

    % Go trials with punishments:
    fprintf('Go actions with punishment\n');
    selTrials   = find(out.isgo == 1 & out.fb.all == 4);
    save(fullfile(dirs.target, sprintf('selTrials_GoPunTrials_Sub%03d.mat', iSub)), 'selTrials');
    
    % NoGo trials with rewards:
    fprintf('NoGo actions with reward\n');
    selTrials   = find(out.isgo == 0 & out.fb.all == 1);
    save(fullfile(dirs.target, sprintf('selTrials_NoGoRewTrials_Sub%03d.mat', iSub)), 'selTrials');

    % NoGo trials with no rewards:
    fprintf('NoGo actions with no reward\n');
    selTrials   = find(out.isgo == 0 & out.fb.all == 2);
    save(fullfile(dirs.target, sprintf('selTrials_NoGoNoRewTrials_Sub%03d.mat', iSub)), 'selTrials');

    % NoGo trials with no punishments:
    fprintf('NoGo actions with no punishment\n');
    selTrials   = find(out.isgo == 0 & out.fb.all == 3);
    save(fullfile(dirs.target, sprintf('selTrials_NoGoNoPunTrials_Sub%03d.mat', iSub)), 'selTrials');

    % NoGo trials with punishments:
    fprintf('NoGo actions with punishment\n');
    selTrials   = find(out.isgo == 0 & out.fb.all == 4);
    save(fullfile(dirs.target, sprintf('selTrials_NoGoPunTrials_Sub%03d.mat', iSub)), 'selTrials');
 
end % end iSub.

end % END OF FUNCTION.