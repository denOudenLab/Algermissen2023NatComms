function EEGfMRIPav_OutcomeLocked_4_rejICs(rootDir, subVec)

% EEGfMRIPav_OutcomeLocked_4_rejICs(rootDir, subVec)
% 
% Reject identified components (after ICA and visual inspection of plots),
% save data.
% 
% INPUTS:
% rootDir           = string, root directory of project
% subVec            = vector of integers, numbers of subjects to process
%
% OUTPUTS:
% saves data after component rejection to disk under dirs.postICA.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Preprocessing/

% ----------------------------------------------------------------------- %
%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootDir', 'var')
    rootDir = preprocessing_set_rootDir(); % '/project/3017042.02';
    fprintf('rootDir unspecified, assume %s\n',rootDir);
end

if ~exist('subVec', 'var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n', strjoin(string(subVec), ', '));
end

% ----------------------------------------------------------------------- %
%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootDir, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Helpers'));

% Set directories and parameters:
dirs    = set_dirs(rootDir);
par     = set_par();

% ----------------------------------------------------------------------- %
%% Overview of components to reject:

comp2remove     = { % sID, [IC to remove]
    1, [1 2 3 4 5 6 7 8 9 10 12 13 16 20 21 22 23 24 25 27 29 30 31 32 33 34 35 37 38 39 40 41 42 43 44 45 46 47 51 52 54 55];
    2, [1 2 3 4 6 7 8 11 12 14 15 16 20 22 23 30 31 33 36 41 43 44 46 47 51 52 56 58];
    3, [1 3 4 5 6 7 8 9 11 13 15 18 19 21 23 24 27 28 31 33 36 37 39 44 45 46 47 48 49 50 51 52 53 59 62];
    4, [1 2 3 4 5 6 7 8 9 10 11 12 13 15 17 21 22 23 24 25 26 31 32 35 36 37 38 40 41 42 44 45 47 53 54];
    5, [1 2 3 5 10 11 12 15 17 19 22 23 25 29 30 31 32 33 36 37 40 41 42 43 44 45 48 49 50 51 57 59];
    6, [1 2 3 6 7 8 9 13 18 21 24 32 33 34 35 36 37 40 41 44 45 53 54 55 56 58 59 60];
    7, [1 2 3 4 5 6 8 9 10 11 12 16 17 18 20 21 23 25 26 29 31 32 33 37 38 39 41 43 44 47 51 55 56 59 62];
    8, [1 2 4 5 6 8 14 15 16 20 25 27 28 30 35 36 37 39 41 42 43 45 46 51 55];
    9, [1 2 3 4 5 8 10 11 13 17 19 20 22 23 24 27 29 31 34 36 41 42 43 45 49 52 54 57 58];
    10, [1 2 3 4 5 6 10 13 22 26 27 28 29 30 31 32 35 37 39 43 44 45 48 52 53 58 60 61];
    11, [1 2 3 4 5 6 7 8 9 11 12 13 14 15 18 19 20 21 22 23 24 25 26 27 30 32 34 38 42 44 45 46 47 48 51 52 53 54 56 57 58 59 60 61 62];
    12, [1 2 3 5 6 11 16 19 20 21 22 23 24 28 29 34 35 36 39 41 42 44 45 48 53 56];
    13, [1 2 3 4 5 6 7 8 9 13 15 19 21 23 24 27 28 30 33 36 41 45 46 47 48 49 50 51 52 53 55 56];
    14, [1 2 3 4 5 6 9 10 12 16 17 19 20 24 28 29 31 35 43 47 49 50 52 53 54 55 59];
    15, [1 2 3 4 5 8 9 10 15 16 26 28 29 38 39 40 41 44 45 46 47 48 49 51 52 56 60 62];
    16, [1 2 3 4 5 6 7 11 13 14 15 16 17 18 19 22 28 30 31 33 35 36 37 38 39 40 42 44 45 46 50 52 58 61 62];
    17, [1 3 4 5 6 7 8 9 12 21 23 25 26 27 31 32 35 36 37 38 39 41 42 43 44 58];
    18, [1 2 3 4 5 8 9 10 17 18 21 22 27 28 29 32 34 35 37 38 40 49 50 52];
    19, [1 2 3 4 6 7 8 9 10 12 14 15 17 20 21 22 23 24 27 30 33 34 36 38 39 40 43 44 45 47 49 51];
    20, [1 2 3 4 6 9 12 14 15 17 19 21 23 24 28 29 30 31 32 34 35 36 37 38 40 41 43 46 47 48 49 51 56 57];
    21, [1 2 3 4 5 6 7 9 10 15 16 17 20 21 23 25 27 28 29 30 31 34 35 36 37 38 39 41 42 43 45 46 49 50];
    22, [1 2 3 4 8 10 16 18 20 22 23 25 27 28 29 31 32 33 34 35 36 37 40 44 45 46 50 51 53];
    23, [1 3 4 5 7 9 10 12 13 14 16 18 19 22 23 24 26 27 28 29 30 31 33 37 38 39 40 41 42 49 50 51 52 55 56 60 61];
    24, [1 2 3 4 7 11 14 15 25 26 28 29 30 31 32 33 34 36 37 38 39 44 45 47 50 51 52 53 54 56 58];
    25, [1 2 3 4 5 6 7 8 9 11 12 13 14 17 19 22 23 24 25 26 27 28 31 33 34 35 36 37 38 39 40 41 42 43 44 46 47 48 49 51 53 54 56 57 58];
    26, [1 2 3 4 6 8 19 23 25 28 29 30 32 33 36 37 38 39 41 43 44 46 47 48 49 50 51 55 56 58];
    27, [1 2 3 4 5 6 8 9 14 18 19 23 24 25 28 31 33 34 35 36 38 40 41 42 43 44 45 46 47 48 49 50 52 53 54 55 56 57 58 59];
    28, [1 2 3 4 5 6 7 8 14 16 17 23 27 29 31 32 33 34 35 36 37 40 42 43 44 45 46 47 55 56];
    29, [1 2 3 4 5 6 8 10 18 19 22 25 26 27 29 30 33 35 36 37 39 40 41 42 46 47 48 49 50 51 52 53 54 56];
    30, [1 2 4 5 7 8 9 10 11 14 15 16 17 18 19 20 21 22 23 25 27 28 29 30 31 32 33 34 35 36 37 38 39 42 45 46 48 53];
    31, [1 2 3 4 5 6 10 15 19 20 21 22 23 24 25 26 28 31 32 35 37 38 39 40 43 46 48 49 51 52 53 54];
    32, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 17 18 20 23 24 28 32 34 48 49 53 56 57 58 59];
    33, [1 2 3 4 5 6 7 9 10 11 12 13 14 16 21 22 26 27 28 29 31 32 33 34 35 37 38 39 40 41 42 43 44 45 46 49 51 52 56];
    34, [1 2 3 4 5 6 7 8 14 20 25 26 29 30 32 33 34 37 38 39 42 45 44 46 48 49 51 52 53 55];
    35, [1 2 3 5 6 7 8 9 10 11 15 16 21 25 27 29 30 31 35 36 38 39 40 41 42 43 45 46 47 49 50 51 52 55 56];
    36, [1 2 3 4 5 6 7 9 10 12 15 16 20 22 24 25 26 28 32 33 36 37 39 40 41 42 43 44 45 46 47 48 49 50 51 52 57 58];
    };

% ----------------------------------------------------------------------- %
%% Check how many got removed:

nSub    = 36;
nComp   = nan(nSub, 1); % initialize
for iSub = 1:nSub   
    nComp(iSub) = length(comp2remove{iSub, 2});    
end

fprintf('Number of rejected components: M = %.03f, SD = %.03f, range %d - %d \n', ...
    mean(nComp), std(nComp), min(nComp), max(nComp));
% Number of rejected components: M = 32.694, SD = 5.350, range 24 - 45 

fprintf('Number of rejected components: %s\n',num2str(sort(nComp, 'descend')', ' %d')); % pretty continuous at high end, makes not much sense to set arbitrary cut-off...

% % Reject invalid subjects:
% invalidSubs = [11 12 23 30];
% nComp(invalidSubs)' % 45    26    37    38
% nComp_keep = nComp(setdiff(1:size(comp2remove, 1), invalidSubs));
% fprintf('Number of rejected trials: M = %.03f, SD = %.03f, range %d - %d \n',mean(nComp_keep), std(nComp_keep),min(nComp_keep),max(nComp_keep));
% Number of rejected trials: M = 32.219, SD = 4.924, range 24 - 45 
% % --> Not much difference...

% ----------------------------------------------------------------------- %
%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.ICA, '*.mat'));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n', nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList, 8, 10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames, subVec));

% ----------------------------------------------------------------------- %
%% Loop over subjects.

for iSub = selIndices % iSub = 1;
    
    % ------------------------------------------------------------------- %
    %% Create output file name:
    
    outputFile      = fullfile(dirs.postICA, sprintf('MRIpav_%03d_postIC.mat', iSub));
    fprintf('EEG output file is %s\n', outputFile);
    
    % Check if output file already exists:
    if exist(outputFile, 'file')
        warning('Subject %03d: File %s already exists, skip subject\n', iSub, outputFile);
        continue;
    end

    % ------------------------------------------------------------------- %
    %% Create input file name:
    
    inputFile       = fullfile(dirs.ICA, fileList{iSub}); % go to subject directory
    fprintf('EEG input file is %s\n', inputFile);
    
    % ------------------------------------------------------------------- %
    %% Load data:
    
    fprintf('Subject %03d: Loading data ... (can take a few seconds)\n', iSub);
    load(inputFile) % 
    fprintf('Subject %03d: Loading completed\n', iSub);
        
    % ------------------------------------------------------------------- %
    %% Remove selected ICs:
    
    rejComps = comp2remove{[comp2remove{:, 1}] == iSub, 2};
    fprintf('Subject %03d: Reject components %s \n', iSub, strjoin(string(rejComps), ', '));

    cfg                     = [];
    cfg.component           = rejComps;
    data                    = ft_rejectcomponent(cfg, comp);
   
    % ------------------------------------------------------------------- %
    %% Save data: 

    fprintf('Subject %03d: Start saving data\n', iSub);
    
    hist.postICA.dirs       = dirs;
    hist.postICA.par        = par;
    
    save(outputFile, 'data', 'hist', '-v7.3');
    clear comp data hist

    fprintf('Subject %03d: Finished saving data\n', iSub);
    
end % end iSub-loop.

fprintf('Done :-)\n');

end % END OF FUNCTION.