function [job, data] = time_prepare_generic_data(job, data)

% [job, data] = time_prepare_generic_data(job, data)
%
% Add aggregated time-domain data sets per subject across conditions, per 
% condition across subjects, and both.
% Creates data sets that are only dependent on the conditions used for
% aggregation; NOT on contrast/ band/ channels etc.
% Those data sets are then used in time_prepare_contrast_data.m
% 
% INPUTS:
% job               = cell, created via time_update_job.m, needs at least fields 
%   .nSub           = integer, number of subjects
%   .nCond          = integer, number of conditions
%   .validSubs      = numeric vector, subject numbers of valid subjects to
%   be included in analyses
% data              = cell, need at least the following fields:
% 	ERPdata{iSub, iCond} = aggregated data per subject per condition,
% 	created via EEGfMRIPav_OutcomeLocked_8_time_grouplevel_create.m
%
% OUTPUTS:
% job               = cell, no fields added
% data              = cell, with the following fields added:
%   .tmpSubCond     = average data per subject per condition in 2D cell
%   format.
%   .Sub            = grand average per subject across conditions.
%   .Cond           = grand average per condition across subjects.
%   .mu             = grand average across both conditions and subjects.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

fprintf('Prepare generic data sets\n');

% ----------------------------------------------------------------------- %
%% Average within subjects over conditions (grand average per subject):

data.Sub = cell(job.nSub, 1);
for iSub = job.validSubs % 1:job.nSub; % iSub = 1;
    % Average over conditions (grand average per subject):
    cfg                         = [];
    cfg.parameter               = 'avg';
    data.Sub{iSub}              = ft_timelockgrandaverage(cfg, data.ERPdata{iSub, :}); % average over all conditions
    nanVec                      = isnan(data.Sub{iSub}.avg);
    powVec                      = length(data.Sub{iSub}.avg(:));
    fprintf('Subject %d, percent NaNs: %d percent\n', ...
        iSub, sum(nanVec(:)) / powVec * 100);
end

% ----------------------------------------------------------------------- %
%% Average within conditions over subjects (grand average per condition):
% Use only valid subjects:

data.Cond                       = cell(job.nCond, 1);
for iCondi = 1:job.nCond
    fprintf('Condition %d \n', iCondi);
    cfg                         = [];
    cfg.parameter               = 'avg';
    data.Cond{iCondi}           = ft_timelockgrandaverage(cfg, ...
        data.ERPdata{job.validSubs, iCondi}); % average over subjects
end

% ------------------------------------------------------------------ %
%% Overall average of subjects and conditions (reference for time, job.channels, etc.):
% Use output of previous step, thus use only valid subjects:

cfg                             = [];
cfg.parameter                   = 'avg';
data.mu                         = ft_timelockgrandaverage(cfg, data.Cond{:}); % average over conditions

end % END OF FUNCTION.