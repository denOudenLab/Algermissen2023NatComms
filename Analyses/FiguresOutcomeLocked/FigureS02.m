% FigureS02.m

% Plots for Supplementary Figure S02.
% The same approach is used for Figure 3 and for Supplementary Figures S08.
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
%% Figure S02A: bar plots:

% Add additional paths for scripts:
addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Run_ROI/Plot_PEs'));

% For sub-selection of subjects,
% manipulate validSubs inside plot_PEs_ROI_outcome.m!!!!!!!!!!!!!!!!!!!!!!!

% Use results from GLM2.
% Run preprocess_PEs_outcomeLocked.R first to extract PEs into .csv files.
% See also in plot_PEs_ROI_outcome_loop.m.

GLMID       = '2';
ROI2use     = {'GLM2vmPFCValence'; 'GLM2StriatumValence'};
ROITitle    = {'vmPFC'; 'Striatum'};

% Loop over ROIs:
for iPlot = 1:length(ROI2use)
    
    % Sorted by actions, then by outcomes:
    data = plot_PEs_ROI_outcome(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, 'actOut', true, false, [-0.4 0.4]); % with points, without legend
    
    % Save plot:
    saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_%s.png', ROITitle{iPlot})));
    close gcf
    
    % Save source data:
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_%s.csv', ROITitle{iPlot}));
    csvwrite(fullFileName, data);
    
end

% ----------------------------------------------------------------------- %
%% Figure S02A: Brain slices displaying outcome valence:

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.GLMID           = '2without7';

job.zLim            = [0 3.0]; % 

% Sagittal:
job.iCope = 2; job.iView = 'sagittal'; job.iSlice = 4; job.cLim = 150; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
job.iCope = 2; job.iView = 'coronal'; job.iSlice = 6; job.cLim = 150; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
job.iCope = 2; job.iView = 'axial'; job.iSlice = -10; job.cLim = 150; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02A_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S02B: Brain slices displaying only contrast PE_STD:

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

dirs.save = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.GLMID           = '1without7'; %

job.zLim            = [0 3.0]; %

% Sagittal:
job.iCope = 4; job.iView = 'sagittal'; job.iSlice = 4; job.cLim = 90; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02B_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02B_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02B_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
job.iCope = 4; job.iView = 'coronal'; job.iSlice = 6; job.cLim = 90; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02B_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02B_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02B_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
job.iCope = 4; job.iView = 'axial'; job.iSlice = -10; job.cLim = 90; job.sign = 'pos';
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02B_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02B_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, job.sign)));
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02B_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S02C: Brain slices displaying conjunction at different levels of cluster-forming threshold:

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'conjunction';
job.GLMID           = '1without7';
job.firstCope       = 4;
job.secondCope      = 5;
job.sign            = 'pos';

job.cLim            = [0.0 1.0]; % keep rather narrow

% Sagittal:
job.iView = 'sagittal'; job.iSlice = 4; job.sign = 'pos'; % compromise ACC and vmPFC (no striatum)
EEGfMRIPav_sliceDisplay_conjunction(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02C_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign))); 
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02C_thresh_zstat1_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
job.iView = 'coronal'; job.iSlice = 6; job.sign = 'pos'; % striatum and ACC
EEGfMRIPav_sliceDisplay_conjunction(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02C_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign))); 
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02C_thresh_zstat1_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
job.iView = 'axial'; job.iSlice = -10; job.sign = 'pos'; % striatum
EEGfMRIPav_sliceDisplay_conjunction(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02C_%s%d.png', ...
    job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign))); 
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS02C_thresh_zstat1_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% END OF FILE.