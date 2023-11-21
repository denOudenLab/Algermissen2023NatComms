function [simulations, job] = EEGfMRIPav_cbm_sim(job)

% [simulations, job] = EEGfMRIPav_cbm_sim(nMod)
%
% Perform model simulations or one-step-ahead predictions for given model for given type of parameters.
%
% INPUTS:
% job 	      = structure with the following fields:
% .simType    = string, type of simulation, either 'modSim' (model simulations) or 'osap' (one-step-ahead predictions).
% .parType    = string, type of input parameters, either 'lap' (LaPlace approximation) or 'hbi' (Hierarchical Bayesian inference).
% .nMod       = integer, model number to be simulated.
%
% OUTPUTS:
% Save to disk.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/CBM_Scripts/

% clear all; close all; clc

% ----------------------------------------------------------------------- %
%% Complete input fields:

if ~isfield(job, 'simType')
    job.simType = 'modSim'; % modSim, osap
end

if ~isfield(job, 'parType')
    job.parType = 'lap'; % lap, hbi
end

if ~isfield(job, 'nMod')
    nMod = 1;
end

fprintf('Simulate for model %02d\n', nMod);
fprintf('Simulate data using method %s based on %s fits\n', ...
    job.simType, job.parType);

% ----------------------------------------------------------------------- %
%% Directories:

dirs.root    = '/project/3017042.02'; % root directory, needs to be adjusted to local folder structure.
dirs.scripts = fullfile(dirs.root, 'Analyses/CBM_Scripts');
dirs.results = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');
dirs.lap     = fullfile(dirs.results, 'LAP_Results');
dirs.hbi     = fullfile(dirs.results, 'HBI_Results');
dirs.sim     = fullfile(dirs.results, 'Simulations');
if ~exist(dirs.sim, 'dir'); mkdir(dirs.sim); end

% Add additional paths:
addpath(fullfile(dirs.scripts, 'osaps'));
addpath(fullfile(dirs.scripts, 'modelSimulations'));

% ----------------------------------------------------------------------- %
%% Complete downstream settings:

if strcmp(job.simType, 'osap')
    job.nIter       = 1; % 1 for osap
elseif strcmp(job.simType, 'modSim')
    job.nIter       = 100; % more for modSim
else
    error('Unknown simulation type')
end

% Object name:
job.hbi_suffix  = sprintf('%02d', nMod);

% Complete name of HBI file:
job.hbi_name_mod    = fullfile(dirs.hbi, sprintf('hbi_mod_%s.mat', job.hbi_suffix)); % hbi output object

% Number of HBI model to inspect:
job.hbi_nmod        = job.nMod;

% ----------------------------------------------------------------------- %
%% Define output file name:

outputfile = fullfile(dirs.sim,sprintf('%s_mod%02d_%s_iter%04d.mat', ...
    job.simType, job.nMod, job.parType, job.nIter));

% ----------------------------------------------------------------------- %
%% Start simulation:

if ~exist(outputfile, 'file')

    % ------------------------------------------------------------------- %
    %% Load data:

    % Use original subject data --> original stimulus order
    fprintf('Load data\n');
    inputFile   = fullfile(dirs.results, 'EEGfMRIPav_cbm_inputData.mat');    
    fdata       = load(inputFile); % only contains .data
    data        = fdata.data; % extract .data
    job.nSub    = length(data);

    % ------------------------------------------------------------------- %
    %% Load parameters:

    fprintf('Load parameters based on %s\n', job.parType)
    if strcmp(job.parType, 'lap')
        job.lap_name_mod    = fullfile(dirs.lap,sprintf('lap_mod%02d.mat', job.nMod));
        fname               = load(job.lap_name_mod);
        cbm                 = fname.cbm;
        allParam            = cbm.output.parameters;
    elseif strcmp(job.parType, 'hbi')
        fname               = load(job.hbi_name_mod);
        cbm                 = fname.cbm;
        if length(cbm.output.parameters) == 1 % if only 1 model
            allParam            = cbm.output.parameters{:};
        else
            allParam            = cbm.output.parameters{job.hbi_nmod}; % all models
        end
    end

    % ------------------------------------------------------------------- %
    %% Simulate data:

    fprintf('Simulate %s based on %s parameters for model %02d with %d iterations\n', ...
        job.simType, job.parType, job.nMod, job.nIter)

    for iSub = 1:job.nSub % iSub = 1;

        % Extract subject data:
        fprintf('Start subject %03d\n', iSub)
        if strcmp(job.simType, 'osap') 
            subj = data{iSub}; % retrieve subject data
        elseif strcmp(job.simType, 'modSim') 
            subj = sim_subj; % contains reqactions and feedback to compute feedback validity
        else
            error('Unknown simulation type')
        end

        % Retrieve parameters:
        parameters = allParam(iSub, :);

        % Save per subject:
        sim.parameters{iSub}   = parameters; % add parameters
        sim.subj{iSub}         = subj; % add subject data

        % Iterate:
        for iIter = 1:job.nIter % iIter = 1;

            % Progress bar:
            if strcmp(job.simType, 'modSim') && job.nIter >= 100
                fraction_done = iIter/job.nIter;
                waitbar(fraction_done)
            end

            % Simulate:
            out = eval(sprintf('EEGfMRIPav_cbm_mod%02d_%s(parameters,subj)', job.nMod, job.simType));

            % Save results:
            sim.p(iSub, iIter, :, :)        = out.p;
            sim.pGo(iSub, iIter, :)         = sum(out.p(:,[1 2]), 2); % pGo (average of both Gos)
            sim.PE(iSub, iIter, :)          = out.PE ./ exp(parameters(1)); % normalize PEs to range [-1 1]
            sim.EV(iSub, iIter, :, :, :)    = out.EV;
            sim.lik(iSub, iIter, :)         = out.lik;
            if strcmp(job.simType, 'modSim') 
                sim.action(iSub, iIter, :, :)   = out.action;
                sim.stay(iSub, iIter, :)        = out.stay;
                sim.outcome(iSub, iIter, :)     = out.outcome;
            end

        end % end iIter

    end % end iSub
    fprintf('Finished simulation :-)\n')

    % ------------------------------------------------------------------- %
    %% Save:

    fprintf('Save outputs\n')
    
    job.dirs = dirs;
    save(outputfile, 'sim', 'job', '-v7.3');

    fprintf('Done :-)\n')

% ----------------------------------------------------------------------- %
%% Otherwise load:
else

    fprintf('File %s already exists; \nload file\n', outputfile);
    load(outputfile);
    fprintf('File loaded :-) \n');

end % end of if

end % END OF FUNCTION.