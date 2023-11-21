function taft_save_BOLD_HRF_per_trial(ROI2use)

% taft_save_BOLD_HRF_per_trial(ROI2use)
%
% For selected ROI, loop over subjects, initialize TAfT job, upsample
% volume-by-volume data, epoch to trial-by-trial data, save to disk.
% Mind adjusting the root directory.
%
% INPUTs:
% ROI2use       = string, 1 ROI to upsample, epoch, and save.
%
% OUTPUTS:
% none, save to disk.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here: 
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT

% ----------------------------------------------------------------------- %
% Select ROI:
% ROIs2use = {'GLM1StriatumConj', 'GLM1ACCConjMan', 'GLM1LeftMotorConj', 'GLM1vmPFCConjMan', 'GLM1PCCConj', 'GLM1LeftITGConj', 'GLM1V1Conj'}; % --> GOES INTO PAPER

% ROIs2use = {'GLM1StriatumConj'}; % 
% ROIs2use = {'GLM1ACCConjMan'}; % 
% ROIs2use = {'GLM1LeftMotorConj'}; % 
% ROIs2use = {'GLM1vmPFCConjMan'}; % 
% ROIs2use = {'GLM1PCCConj'}; % 
% ROIs2use = {'GLM1LeftITGConj'}; % 
% ROIs2use = {'GLM1V1Conj'}; % 

if length(ROI2use) > 1; error('Error; more than 1 ROI string specified'); end

% ----------------------------------------------------------------------- %
%% Settings:

nSub                = 36;
nBlock              = 6;
nTrial              = 640;

% ----------------------------------------------------------------------- %
%% Set directories:

dirs.project        = taft_set_rootDir();
dirs.EEG 		    = fullfile(dirs.project, 'Log', 'EEG', 'OutcomeLockedResults');
dirs.TAfT           = fullfile(dirs.EEG, 'TAfT_Betas');
dirs.save           = fullfile(dirs.TAfT, sprintf('TAfT_BOLD_%s', ROI2use{:}));
if ~exist(dirs.save, 'dir'); mkdir(dirs.save); end

% ----------------------------------------------------------------------- %
%% Loop over subjects, extract fMRI, save:

for iSub = 1:nSub % iSub = 1;
    
    % ------------------------------------------------------------------- %
    %% 1) Initialize job:
    
    job                 = taft_preprocess_initialize_job('TF', iSub, ROI2use);
    
    % ------------------------------------------------------------------- %
    %% 2) Select trials:
    
    job.goodTrlIdxfMRI  = 1:nTrial; % select all (!) trials

    % ------------------------------------------------------------------- %
    %% 3) Load, fit, combine, and reshape fMRI-beta data of different ROIs and blocks:
    
    fprintf('Subject %03d: Load fMRI data\n', iSub);
    X                   = taft_preprocess_load_fMRI(job);
    if(length(X) ~= nTrial); error('Length of X differs from number of trials %d', nTrial); end

    % ------------------------------------------------------------------- %
    %% 4) Save as vector in csv file:
    
    % Create file name:
    fileName    = sprintf(fullfile(dirs.save, sprintf('TAfT_BOLD_%s_sub%03d.csv', ...
    strjoin(ROI2use), iSub)));

    % Save as csv:
    csvwrite(fileName, X);
    
end

fprintf('Finished :-)\n');

end % END OF FUNCTION.