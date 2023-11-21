% Figure03.m

% Plots for Figure 3A-C in manuscript.
% The same approach is for Supplementary Figures S02 and S08.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresOutcomeLocked

% clear all; close all; clc

% Set root directory:
cd(fileparts(which('figures_set_rootDir.m')));
dirs.root           = figures_set_rootDir();

% ----------------------------------------------------------------------- %
%% Figure 3A: Bar plots:

% Add additional paths for scripts:
addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Run_ROI/Plot_PEs'));

% For sub-selection of subjects,
% manipulate validSubs inside plot_PEs_ROI_outcome.m!!!!!!!!!!!!!!!!!!!!!!!

% Use results from GLM2.
% Run preprocess_PEs_outcomeLocked.R first to extract PEs into .csv files.
% See also in plot_PEs_ROI_outcome_loop.m.

GLMID       = '2';
roi2use     = {'GLM2vmPFCValence'; 'GLM2StriatumValence'};
roiTitle    = {'vmPFC'; 'Striatum'};

% Loop over ROIs:
for iPlot = 1:length(roi2use)
    
    % Sorted by actions, then by outcomes:
    data = plot_PEs_ROI_outcome(char(roi2use(iPlot)), char(roiTitle(iPlot)), GLMID, 'actOut', true, false, [-0.4 0.4]); % with points, without legend
    
    % Save plot:
    saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_%s.png', roiTitle{iPlot})));
    close gcf
    
    % Save source data:
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_%s.csv', roiTitle{iPlot}));
    csvwrite(fullFileName, data);
    
end

% ----------------------------------------------------------------------- %
%% Figure 3A: Brain slices displaying outcome valence:

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.GLMID           = '2';

job.zLim            = [0 3.0]; % 

% Sagittal:
job.iCope = 2; job.iView = 'sagittal'; job.iSlice = 4; job.cLim = 150; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
job.iCope = 2; job.iView = 'coronal'; job.iSlice = 6; job.cLim = 150; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
job.iCope = 2; job.iView = 'axial'; job.iSlice = -10; job.cLim = 150; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3A_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure 3B: Brain slices displaying only contrast PE_STD:

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

dirs.save = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.GLMID           = '1'; % 1, 1without7, 1B, 1C

job.zLim            = [0 3.0]; %

% Sagittal:
job.iCope = 4; job.iView = 'sagittal'; job.iSlice = 4; job.cLim = 90; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3B_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3B_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3B_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
job.iCope = 4; job.iView = 'coronal'; job.iSlice = 6; job.cLim = 90; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3B_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3B_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3B_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
job.iCope = 4; job.iView = 'axial'; job.iSlice = -10; job.cLim = 90; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3B_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3B_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3B_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure 3C: Brain slices displaying conjunction at different levels of cluster-forming threshold:

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

% GLM1: use this for GLM1, GLM1without7, GLM1B, GLM1C
dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'conjunction';
job.GLMID           = '1'; % 1, 1without7, 1B, 1C
job.firstCope       = 4;
job.secondCope      = 5;
job.sign            = 'pos';

job.cLim            = [0.0 1.0]; % keep rather narrow

% Sagittal:
job.iView = 'sagittal'; job.iSlice = 4; % compromise ACC and vmPFC (no striatum)
EEGfMRIPav_sliceDisplay_conjunction(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3C_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign))); 
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3C_thresh_zstat1_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
job.iView = 'coronal'; job.iSlice = 6; % striatum and ACC
EEGfMRIPav_sliceDisplay_conjunction(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3C_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign))); 
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3C_thresh_zstat1_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
job.iView = 'axial'; job.iSlice = -10; % striatum
EEGfMRIPav_sliceDisplay_conjunction(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3C_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign))); 
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig3C_thresh_zstat1_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% END OF FILE.