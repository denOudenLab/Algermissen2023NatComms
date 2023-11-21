% FigureS08.m

% Plots for Supplementary Figure S08.
% The same approach is used for Figure 3 and for Supplementary Figure S02.
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
%% Figure S08A-C: 

job                 = [];

% Select model to plot:
job.GLMID           = '1'; plotLetter = 'A';
% job.GLMID           = '1B'; plotLetter = 'B';
% job.GLMID           = '1C'; plotLetter = 'C';

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

% Further settings:
dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job.lockSettings    = 'outcomelocked';
job.type            = 'conjunction';
job.firstCope       = 4;
job.secondCope      = 5;
job.sign            = 'pos';

job.cLim            = [0.0 1.0]; % keep rather narrow

% Sagittal:
job.iView = 'sagittal'; job.iSlice = 4; % compromise ACC and vmPFC (no striatum)
EEGfMRIPav_sliceDisplay_conjunction(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS08%s_%s%d.png', ...
    plotLetter, job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign))); 
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS08%s_thresh_zstat1_%s%d.csv', ...
    plotLetter, job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
job.iView = 'coronal'; job.iSlice = 6; % striatum and ACC
EEGfMRIPav_sliceDisplay_conjunction(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS08%s_%s%d.png', ...
    plotLetter, job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign))); 
sliceIdx        = floor(size(data, 2)/2) + 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS08%s_thresh_zstat1_%s%d.csv', ...
    plotLetter, job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
job.iView = 'axial'; job.iSlice = -10; % striatum
EEGfMRIPav_sliceDisplay_conjunction(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS08%s_%s%d.png', ...
    plotLetter, job.iView, job.iSlice)));
close gcf

% Save source data file:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign))); 
sliceIdx        = floor(size(data, 3)/2) - 9 + job.iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS08%s_thresh_zstat1_%s%d.csv', ...
    plotLetter, job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% END OF FILE.