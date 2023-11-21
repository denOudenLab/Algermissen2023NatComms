% EEGfMRIPav_GLM2_labels.m

% Create labels for each contrast of GLM2.
% For GLM2without7, just copy and rename .mat output file.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/fMRI_Scripts/Display_Results/
% clear all; close all; clc

GLMlabels = {
                1, 'Action (Go - NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                2, 'Valence (Positive - Negative)', '\beta_{Positive} - \beta_{Negative} (a.u.)', 'pos'; 
                2, 'Salience (Salient - Neutral)', '\beta_{Salient} - \beta_{Neutral} (a.u.)', 'pos'; 
                4, 'Hand Sum (Left + Right)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                5, 'Hand Difference (Left - Right)', '\beta_{Left} - \beta_{Right} (a.u.)', 'pos'; 
                };
rootDir     = '/project/3017042.02';
fileName    = fullfile(rootDir, 'Analyses/fMRI_Scripts/Display_Results/Labels/GLM2labels.mat');
save(fileName, 'GLMlabels');

% END OF FILE.
