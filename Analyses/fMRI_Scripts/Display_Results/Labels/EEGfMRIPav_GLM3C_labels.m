% EEGfMRIPav_GLM3C_labels.m

% Create labels for each contrast of GLM3C.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/fMRI_Scripts/Display_Results/
% clear all; close all; clc

GLMlabels = { % sID, [filter to select]
                1, 'Executed action (Go . NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                2, 'Cue Valence (Win - Avoid)', '\beta_{Win} - \beta_{Avoid} (a.u.)', 'pos'; 
                3, 'Motor (Left - Right)', '\beta_{Left} - \beta_{Right} (a.u.)', 'pos'; 
                4, 'Standard update', '\beta_{Update_{standard}} (a.u.)', 'pos';
                5, 'Biased minus standard update', '\beta_{Update_{dif}} (a.u.)', 'pos';
                6, 'Midfrontal alpha', '\beta_{Midfrontal alpha} (a.u.)', 'pos';
                };
rootDir     = '/project/3017042.02';
fileName    = fullfile(rootDir, 'Analyses/fMRI_Scripts/Display_Results/Labels/GLM3Clabels.mat');
save(fileName, 'GLMlabels');

% END OF FILE.
