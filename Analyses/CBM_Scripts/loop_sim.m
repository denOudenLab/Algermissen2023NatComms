function loop_sim()

% loop_sim()
%
% Wrapper to call EEGfMRIPav_cbm_sim for model simulations and one-step-ahead predictions.
%
% INPUTS:
% None, set settings for simulations interactively.
%
% OUTPUTS:
% None, EEGfMRIPav_cbm_sim saves to disk.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/CBM_Scripts/

% ----------------------------------------------------------------------- %
%% Initialize root directory, add path:

rootDir = '/project/3017042.02';
addpath(fullfile(rootDir, '/Analyses/Stan_Scripts/CBM_Scripts');

% ----------------------------------------------------------------------- %
%% Initialize job settings:

nMod        = 5;
job.simType = 'modSim'; % modSim, osap
job.parType = 'lap'; % lap, hbi

% ----------------------------------------------------------------------- %
%% Loop over models:

for iMod = 1:nMod
    	
    fprintf('\n\n>>> Start simulating model %02d\n', iMod);
    EEGfMRIPav_cbm_sim(job);
end

end % END OF FUNCTION.
