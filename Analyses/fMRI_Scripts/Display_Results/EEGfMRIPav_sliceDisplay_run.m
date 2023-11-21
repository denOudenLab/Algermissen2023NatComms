% EEGfMRIPav_sliceDisplay_run.m

% Directories:
dirs        = []; 
dirs.root   = '/project/3017042.02'; % root directory--needs to be adapted to users' folder structure

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));


% ----------------------------------------------------------------------- %
%% Outcome-locked: GLM2: Contrasts Action, HandSum, HandDif

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.GLMID           = '2without7'; % 2 2without7

% Action @ outcome:
job.iCope = 1; job.iView = 'coronal'; job.iSlice = 10; job.zLim = [0 5]; job.cLim = 500; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

job.iCope = 1; job.iView = 'coronal'; job.iSlice = -30; job.zLim = [0 5]; job.cLim = 500; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% Same slices as below:
job.iCope = 1; job.iView = 'sagittal'; job.iSlice = 4; job.zLim = [0 5];job.cLim = 500;  % 500
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iCope = 1; job.iView = 'coronal'; job.iSlice = 6; job.zLim = [0 5]; job.cLim = 500; % 500
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iCope = 1; job.iView = 'axial'; job.iSlice = -10;  job.zLim = [0 5]; job.cLim = 500; % 500
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% ---------------------------------------------- %
% Hand Sum @ response:
job.iCope = 4; job.iView = 'coronal'; job.iSlice = 10; job.zLim = [0 5]; job.cLim = 150; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% ---------------------------------------------- %
% Hand Dif @ response:
job.iCope = 5; job.iView = 'coronal'; job.iSlice = -30; job.zLim = [0 5]; job.cLim = 150; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% ----------------------------------------------------------------------- %
%% Outcome-locked: GLM2: Contrast Valence

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.GLMID           = '2without7'; % 2 2without7

% Sagittal:
job.iCope = 2; job.iView = 'sagittal'; job.iSlice = 0; job.zLim = [0 3];job.cLim = 150; % 100
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iCope = 2; job.iView = 'sagittal'; job.iSlice = 4; job.zLim = [0 3];job.cLim = 150;  % 100
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iCope = 2; job.iView = 'sagittal'; job.iSlice = 6; job.zLim = [0 3];job.cLim = 150;  % 100
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% Coronal:
job.iCope = 2; job.iView = 'coronal'; job.iSlice = 6; job.zLim = [0 3]; job.cLim = 150; % 100
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% Axial:
job.iCope = 2; job.iView = 'axial'; job.iSlice = -10;  job.zLim = [0 3]; job.cLim = 150; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% ----------------------------------------------------------------------- %
%% Outcome-locked: Model-based Updatestd and Updatedif
% calls EEGfMRIPav_sliceDisplay_cope, makes red-to-blue plots

dirs.save = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.GLMID           = '1'; % 1 1without7 1B 1C

job.zLim    = [0 3.0]; % 

% ------------------------------------------------ %
% UpdateSTD:

% Sagittal:
job.iCope = 4; job.iView = 'sagittal'; job.iSlice = 0; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iCope = 4; job.iView = 'sagittal'; job.iSlice = 4; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iCope = 4; job.iView = 'sagittal'; job.iSlice = 6; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% Coronal:
job.iCope = 4; job.iView = 'coronal'; job.iSlice = 6; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% Axial:
job.iCope = 4; job.iView = 'axial'; job.iSlice = -10; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% ------------------------------------------------ %
% UpdateDIF:

% Sagittal:
job.iCope = 5; job.iView = 'sagittal'; job.iSlice = 0; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iCope = 5; job.iView = 'sagittal'; job.iSlice = 4; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iCope = 5; job.iView = 'sagittal'; job.iSlice = 6; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% Coronal:
job.iCope = 5; job.iView = 'coronal'; job.iSlice = 6; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% Axial:
job.iCope = 5; job.iView = 'axial'; job.iSlice = -10; job.cLim = 90; % 150
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% ----------------------------------------------------------------------- %
%% Outcome-locked: Conjunction Updatestd & Updatedif (SHORT)

% calls EEGfMRIPav_sliceDisplay_conjunction, makes hot plots

dirs.save   = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'conjunction';
job.GLMID           = '1'; % 1 1without7 1B 1C
job.firstCope       = 4; % 
job.secondCope      = 5; % 
job.sign            = 'pos'; % 

% OLD settings:
job.cLim            = [0.0 1.0]; % keep rather narrow

% job.cLim            = [0.0 1.0]; % keep rather narrow
% job.zLim            = [0 5.0]; % 

% Sagittal:
job.iView = 'sagittal'; job.iSlice = 4; % compromise ACC and vmPFC (no striatum)
EEGfMRIPav_sliceDisplay_conjunction(dirs, job); close gcf;


% Coronal:

job.iView = 'coronal'; job.iSlice = 6; % striatum and ACC
EEGfMRIPav_sliceDisplay_conjunction(dirs, job); close gcf;

% Axial:
job.iView = 'axial'; job.iSlice = -10; % striatum
EEGfMRIPav_sliceDisplay_conjunction(dirs, job); close gcf;

fprintf('Finished with conjunction plots:-) \n');

% ----------------------------------------------------------------------- %
%% GLMs 3A-3B: Using EEG as regressor:

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.iCope           = 6; 

% Theta:
job.GLMID           = '3A';
job.cLim            = 4;
job.iView = 'sagittal'; job.iSlice = 0;
EEGfMRIPav_sliceDisplay_cope(dirs, job);

% Beta:
job.GLMID           = '3B';
job.cLim            = 250;
job.iView = 'sagittal'; job.iSlice = 0;
EEGfMRIPav_sliceDisplay_cope(dirs, job);

% Alpha:
job.GLMID           = '3C';
job.cLim            = 35;
job.iView = 'sagittal'; job.iSlice = 0;
EEGfMRIPav_sliceDisplay_cope(dirs, job);

% ----------------------------------------------------------------------- %
% Other settings for GLMs with EEG as regressor:

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.iCope           = 6; 

% Theta:
job.GLMID           = '3A';
job.cLim            = 4;

% Beta:
% job.GLMID           = '3B';
% job.cLim            = 250;

% Alpha:
% job.GLMID           = '3C';
% job.cLim            = 35;

job.iView = 'sagittal'; job.iSlice = 0;
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iView = 'coronal'; job.iSlice = 10;  
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iView = 'coronal'; job.iSlice = 0;  
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iView = 'coronal'; job.iSlice = -60;  
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

job.iView = 'axial'; job.iSlice = -10; 
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;
job.iView = 'axial'; job.iSlice = 40; 
EEGfMRIPav_sliceDisplay_cope(dirs, job); close gcf;

% END OF FILE.
