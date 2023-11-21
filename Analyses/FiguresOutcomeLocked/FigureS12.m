% FigureS12.m

% Create source data files for Supplementary Figure S12.
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
%% Figure S12A:

% Specify:
maskFile        = 'GLM1ACCUpdateMan.nii';
% Unzip:
gunzip(fullfile(dirs.mask, [maskFile '.gz']));
% Load:
data            = niftiread(fullfile(dirs.mask, maskFile));

% Sagittal:
iView           = 'sagittal';
iSlice          = 0;
sliceIdx        = floor(size(data, 1)/2) - iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12A_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
iView           = 'coronal';
iSlice          = 18;
sliceIdx        = floor(size(data, 2)/2) + 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12A_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
iView           = 'axial';
iSlice          = 26;
sliceIdx        = floor(size(data, 3)/2) - 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12A_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S12B:

% Specify:
maskFile        = 'GLM1PCCConj.nii';
% Unzip:
gunzip(fullfile(dirs.mask, [maskFile '.gz']));
% Load:
data            = niftiread(fullfile(dirs.mask, maskFile));

% Sagittal:
iView           = 'sagittal';
iSlice          = -2;
sliceIdx        = floor(size(data, 1)/2) - iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12B_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
iView           = 'coronal';
iSlice          = -26;
sliceIdx        = floor(size(data, 2)/2) + 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12B_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
iView           = 'axial';
iSlice          = 38;
sliceIdx        = floor(size(data, 3)/2) - 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12B_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S12C:

% Specify:
maskFile        = 'GLM1LeftMotorConj.nii';
% Unzip:
gunzip(fullfile(dirs.mask, [maskFile '.gz']));
% Load:
data            = niftiread(fullfile(dirs.mask, maskFile));

% Sagittal:
iView           = 'sagittal';
iSlice          = -50;
sliceIdx        = floor(size(data, 1)/2) - iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12C_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
iView           = 'coronal';
iSlice          = -26;
sliceIdx        = floor(size(data, 2)/2) + 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12C_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
iView           = 'axial';
iSlice          = 34;
sliceIdx        = floor(size(data, 3)/2) - 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12C_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S12D:

% Specify:
maskFile        = 'GLM1LeftITGConj.nii';
% Unzip:
gunzip(fullfile(dirs.mask, [maskFile '.gz']));
% Load:
data            = niftiread(fullfile(dirs.mask, maskFile));

% Sagittal:
iView           = 'sagittal';
iSlice          = -50;
sliceIdx        = floor(size(data, 1)/2) - iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12D_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
iView           = 'coronal';
iSlice          = -52;
sliceIdx        = floor(size(data, 2)/2) + 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12D_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
iView           = 'axial';
iSlice          = -16;
sliceIdx        = floor(size(data, 3)/2) - 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12D_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% ----------------------------------------------------------------------- %
%% Figure S12E:

% Specify:
maskFile        = 'GLM1V1Conj.nii';
% Unzip:
gunzip(fullfile(dirs.mask, [maskFile '.gz']));
% Load:
data            = niftiread(fullfile(dirs.mask, maskFile));

% Sagittal:
iView           = 'sagittal';
iSlice          = 22;
sliceIdx        = floor(size(data, 1)/2) - iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12E_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Coronal:
iView           = 'coronal';
iSlice          = -82;
sliceIdx        = floor(size(data, 2)/2) + 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, sliceIdx, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12E_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% Axial:
iView           = 'axial';
iSlice          = -14;
sliceIdx        = floor(size(data, 3)/2) - 9 + iSlice/2; % convert into MNI coordinates
saveMat         = squeeze(data(:, :, sliceIdx)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS12E_%s%d.csv', iView, iSlice));
csvwrite(fullFileName, saveMat);

% END OF FILE.