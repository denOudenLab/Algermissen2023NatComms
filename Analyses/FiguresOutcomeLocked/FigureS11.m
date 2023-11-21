% FigureS11.m

% Create source data files for Supplementary Figure S11.
% Plots were created as screenshots from FSLeyes.
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

dirs.mask       = fullfile(dirs.root, 'Log/fMRI/fMRI_Masks/masksTAfTOutcomeLocked');
cd(dirs.mask);

% ----------------------------------------------------------------------- %
%% Figure S11A:

% Specify:
maskFile        = 'GLM2vmPFCValence.nii';
% Unzip:
gunzip(fullfile(dirs.mask, [maskFile '.gz']));
% Load:
data            = niftiread(fullfile(dirs.mask, maskFile));

% Sagittal:
iView           = 'sagittal';
iSlice          = 0;
sliceIdx        = floor(size(data, 1)/2) - iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11A_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
iView           = 'coronal';
iSlice          = 42;
sliceIdx        = floor(size(data, 2)/2) + 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11A_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
iView           = 'axial';
iSlice          = -8;
sliceIdx        = floor(size(data, 3)/2) - 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11A_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S11B:

% Specify:
maskFile        = 'GLM2StriatumValence.nii';
% Unzip:
gunzip(fullfile(dirs.mask, [maskFile '.gz']));
% Load:
data            = niftiread(fullfile(dirs.mask, maskFile));

% Sagittal:
iView           = 'sagittal';
iSlice          = 10;
sliceIdx        = floor(size(data, 1)/2) - iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11B_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
iView           = 'coronal';
iSlice          = 12;
sliceIdx        = floor(size(data, 2)/2) + 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11B_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
iView           = 'axial';
iSlice          = -8;
sliceIdx        = floor(size(data, 3)/2) - 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11B_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S11C:

% Specify:
maskFile        = 'GLM1vmPFCConjMan.nii';
% Unzip:
gunzip(fullfile(dirs.mask, [maskFile '.gz']));
% Load:
data            = niftiread(fullfile(dirs.mask, maskFile));

% Sagittal:
iView           = 'sagittal';
iSlice          = 10;
sliceIdx        = floor(size(data, 1)/2) - iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11C_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
iView           = 'coronal';
iSlice          = 46;
sliceIdx        = floor(size(data, 2)/2) + 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11C_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
iView           = 'axial';
iSlice          = 8;
sliceIdx        = floor(size(data, 3)/2) - 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11C_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S11D:

% Specify:
maskFile        = 'GLM1StriatumConj.nii';
% Unzip:
gunzip(fullfile(dirs.mask, [maskFile '.gz']));
% Load:
data            = niftiread(fullfile(dirs.mask, maskFile));

% Sagittal:
iView           = 'sagittal';
iSlice          = 10;
sliceIdx        = floor(size(data, 1)/2) - iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11D_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
iView           = 'coronal';
iSlice          = 12;
sliceIdx        = floor(size(data, 2)/2) + 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11D_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
iView           = 'axial';
iSlice          = 0;
sliceIdx        = floor(size(data, 3)/2) - 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS11D_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% END OF FILE.