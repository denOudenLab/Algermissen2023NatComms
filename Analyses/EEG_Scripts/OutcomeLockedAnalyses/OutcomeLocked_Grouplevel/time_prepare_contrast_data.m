function [job, data] = time_prepare_contrast_data(job, data)

% [job, data] = time_prepare_generic_data(job, data)
%
% Add aggregated time-domain data sets for selected contrast, channels 
% to job object previously set up via time_update_job.m
% Preliminary data sets created via time_prepare_contrast_data.m
% 
% INPUTS:
% job               = cell, created via time_update_job.m, needs at least
% fields:
%   .nSub           = integer, number of subjects.
%   .nCond          = integer, number of conditions.
%   .validSubs      = numeric vector, subject numbers of valid subjects to
%   be included in analyses.
%   .channels       = vector of strings, selected channels.
%   .contrastType   = string, contrast to be used: 'Preferred', 'Action',
%   'GoPreferred', 'SalientPreferred', 'SalientAction'.
%   .responseSettings   = string, type of response setting for which
%   conditions are split, 'Go' (Go/ NoGo), 'Hand' (Left Go/ Right Go/
%   NoGo), 'none' (averaged over Go/NoGo).
%   .outcomeSettings    = string, type of outcome coding for which 
%   conditions are split, 'abs' (positive/neutral/negative), 'rel'
%   (positive/negative), 'all' (reward/no reward/ no punishment/
%   punishment).
% data              = cell, need at least the following fields:
%   .ERPdata{iSub, iCond} = aggregated data per subject per condition,
% 	created via EEGfMRIPav_OutcomeLocked_8_time_create.m
%   .Cond           = grand average per condition across subjects.
%   .mu             = grand average across both conditions and subjects.
%
% OUTPUTS:
% job               = cell, with the following fields added:
%   .chanIdx        = numeric vector, indices of selected channels.
% data              = cell, with the following fields added:
%   .SubCondTime    = per subject per condition, over time, averaged over
%   channels.
%   .time1 and time2= averaged over conditions, for permutation tests.
%   .mat1 and mat2  = averaged over channels/ conditions, for two-line 
%   plots.
%   .matMean        = mean of mat1 and mat2 (over all conditions).
%   .matGrandMean   = matMean averaged over subjects.
%   .SubTime        = SubCondTime averaged over conditions.
%   .GrandTime      = SubTime averaged over subjects.
%   .topo2plot      = averaged over conditions/ subjects, for topoplots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/OutcomeLockedAnalyses/OutcomeLocked_Grouplevel/

fprintf('Prepare data for contrast %s, channels %s\n', job.contrastType, strjoin(job.channels, '/'));

% ----------------------------------------------------------------------- %
%% Retrieve indices of selected channels:

% Channel indices:
job.chanIdx = find(ismember(data.mu.label, job.channels)); % determine indices of channels per subject

% ----------------------------------------------------------------------- %
%% Keep subjects and conditions and time, average over selected channels:

% a) Average within subject and within subject/condition over channels:
data.SubCondTime = zeros(job.nSub, job.nCond, length(data.mu.time));

for iSub = job.validSubs % 1:parAll.nSub
    for iCondi = 1:job.nCond
        fprintf('Subject %d, Condition %d: Average over channels %s \n', iSub, iCondi, strjoin(job.channels, '/'));
        subChanIdx = ismember(data.ERPdata{iSub, iCondi}.label, job.channels); % determine indices of channels per subject
        data.SubCondTime(iSub, iCondi, :) = nanmean(data.ERPdata{iSub, iCondi}.avg(subChanIdx, :), 1); % average over channels
    end
end

% ----------------------------------------------------------------------- %
%% Set and initialize data objects :

clear data.time1 data.time2 data.mat1 data.mat2 data.topoplot

% Initialize objects:

% a) Permutation test:
data.time1 = cell(job.nSub, 1);
data.time2 = cell(job.nSub, 1);

% b) 2-line plot:
data.mat1 = nan(1, length(data.mu.time)); 
data.mat2 = nan(1, length(data.mu.time)); 

% c) Topoplot:
data.topo2plot = data.Cond{2}; 

% ----------------------------------------------------------------------- %
%% Compute objects WITH LOOPING over subjects (i.e. permutation test):

for iSub = 1:job.nSub
    
    if strcmp(job.contrastType, 'BOLD')

        data.time1{iSub} = data.ERPdata{iSub, 1}; % low BOLD.
        data.time2{iSub} = data.ERPdata{iSub, 2}; % high BOLD.
        
    elseif strcmp(job.contrastType, 'Preferred')
        if strcmp(job.responseSettings, 'Go')
            if strcmp(job.outcomeSettings, 'rel')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize preferred.
                data.time1{iSub}.avg = real(data.ERPdata{iSub, 1}.avg + data.ERPdata{iSub, 3}.avg) ./ 2; % mean of GoPreferred and NoGoPreferred
                data.time2{iSub} = data.ERPdata{iSub, 2}; % initialize non-preferred.
                data.time2{iSub}.avg = real(data.ERPdata{iSub, 2}.avg + data.ERPdata{iSub, 4}.avg) ./ 2; % mean of GoNonPreferred and NoGoNonPreferred

            elseif strcmp(job.outcomeSettings, 'abs')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize preferred.
                data.time1{iSub}.avg = real(data.ERPdata{iSub, 1}.avg + data.ERPdata{iSub, 4}.avg) ./ 2; % mean of GoPreferred and NoGoPreferred
                data.time2{iSub} = data.ERPdata{iSub, 2}; % initialize non-preferred.
                data.time2{iSub}.avg = real(data.ERPdata{iSub, 3}.avg + data.ERPdata{iSub, 6}.avg) ./ 2; % mean of GoNonPreferred and NoGoNonPreferred

            elseif strcmp(job.outcomeSettings, 'all') 
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize preferred.
                data.time1{iSub}.avg = real(...
                    data.ERPdata{iSub, 1}.avg + ...
                    data.ERPdata{iSub, 3}.avg + ...
                    data.ERPdata{iSub, 5}.avg + ...
                    data.ERPdata{iSub, 7}.avg) ./ 4; % mean of GoReward and GoNoPunishment and NoGoReward and NoGoNoPunishment
                data.time2{iSub} = data.ERPdata{iSub, 2}; % initialize non-preferred.
                data.time2{iSub}.avg = real(...
                    data.ERPdata{iSub, 2}.avg + ...
                    data.ERPdata{iSub, 4}.avg + ...
                    data.ERPdata{iSub, 6}.avg + ...
                    data.ERPdata{iSub, 8}.avg) ./ 4; % mean of GoNoReward and GoPunishment and NoGoNoReward and NoGoPunishment
            else
                error('Contrast preferred not available for outcome setting %s', job.outcomeSettings);
            end

        elseif strcmp(job.responseSettings, 'Hand')
            error('Contrast not completed yet');
        else
            error('Unknown response setting');
        end
        
    elseif strcmp(job.contrastType, 'Action')

        if strcmp(job.responseSettings, 'Go')

            if strcmp(job.outcomeSettings, 'rel')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize Go.
                data.time1{iSub}.avg = real(data.ERPdata{iSub, 1}.avg + data.ERPdata{iSub, 2}.avg) ./ 2; % mean of Go.
                data.time2{iSub} = data.ERPdata{iSub, 3}; % initialize NoGo.
                data.time2{iSub}.avg = real(data.ERPdata{iSub, 3}.avg + data.ERPdata{iSub, 4}.avg) ./ 2; % mean of NoGo.

            elseif strcmp(job.outcomeSettings, 'abs')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize Go.
                data.time1{iSub}.avg = real(...
                    data.ERPdata{iSub, 1}.avg + ...
                    data.ERPdata{iSub, 2}.avg + ...
                    data.ERPdata{iSub, 3}.avg) ./ 3; % mean of Go.
                data.time2{iSub} = data.ERPdata{iSub, 4}; % initialize NoGo.
                data.time2{iSub}.avg = real(...
                    data.ERPdata{iSub, 4}.avg + ...
                    data.ERPdata{iSub, 5}.avg + ...
                    data.ERPdata{iSub, 6}.avg) ./ 3; % mean of NoGo.

            elseif strcmp(job.outcomeSettings, 'all')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize Go.
                data.time1{iSub}.avg = real(...
                    data.ERPdata{iSub, 1}.avg + ...
                    data.ERPdata{iSub, 2}.avg + ...
                    data.ERPdata{iSub, 3}.avg + data.ERPdata{iSub, 4}.avg) ./ 4; % mean of Go.
                data.time2{iSub} = data.ERPdata{iSub, 5}; % initialize NoGo.
                data.time2{iSub}.avg = real(...
                    data.ERPdata{iSub, 5}.avg + ...
                    data.ERPdata{iSub, 6}.avg + ...
                    data.ERPdata{iSub, 7}.avg + ...
                    data.ERPdata{iSub, 8}.avg) ./ 4; % mean of NoGo.
            else
                error('Unknown outcome setting');
            end

        elseif strcmp(job.responseSettings, 'Hand')

            if strcmp(job.outcomeSettings, 'rel')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize Go.
                data.time1{iSub}.avg = real(...
                    data.ERPdata{iSub, 1}.avg + ...
                    data.ERPdata{iSub, 2}.avg + ...
                    data.ERPdata{iSub, 3}.avg + ...
                    data.ERPdata{iSub, 4}.avg) ./ 4; % mean of Go.
                data.time2{iSub} = data.ERPdata{iSub, 5}; % initialize NoGo.
                data.time1{iSub}.avg = real(...
                    data.ERPdata{iSub, 5}.avg + ...
                    data.ERPdata{iSub, 6}.avg) ./ 2; % mean of NoGo.

            elseif strcmp(job.outcomeSettings, 'abs')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize Go.
                data.time1{iSub}.avg = real(...
                    data.ERPdata{iSub, 1}.avg + ...
                    data.ERPdata{iSub, 2}.avg + ...
                    data.ERPdata{iSub, 3}.avg + ...
                    data.ERPdata{iSub, 4}.avg + ...
                    data.ERPdata{iSub, 5}.avg + ...
                    data.ERPdata{iSub, 6}.avg) ./ 6; % mean of Go.
                data.time2{iSub} = data.ERPdata{iSub, 7}; % initialize NoGo.
                data.time1{iSub}.avg = real(...
                    data.ERPdata{iSub, 7}.avg + ...
                    data.ERPdata{iSub, 8}.avg + ...
                    data.ERPdata{iSub, 9}.avg) ./ 3; % mean of NoGo.

            elseif strcmp(job.outcomeSettings, 'all')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize Go.
                data.time1{iSub}.avg = real(...
                    data.ERPdata{iSub, 1}.avg + ...
                    data.ERPdata{iSub, 2}.avg + ...
                    data.ERPdata{iSub, 3}.avg + ...
                    data.ERPdata{iSub, 4}.avg + ...
                    data.ERPdata{iSub, 5}.avg + ...
                    data.ERPdata{iSub, 6}.avg + ...
                    data.ERPdata{iSub, 7}.avg + ...
                    data.ERPdata{iSub, 8}.avg) ./ 8; % mean of Go.
                data.time2{iSub} = data.ERPdata{iSub, 9}; % initialize NoGo.
                data.time1{iSub}.avg = real(...
                    data.ERPdata{iSub, 9}.avg + ...
                    data.ERPdata{iSub, 10}.avg + ...
                    data.ERPdata{iSub, 11}.avg + ...
                    data.ERPdata{iSub, 12}.avg) ./ 4; % mean of NoGo.

            else
                error('Unknown outcome setting');
            end
        else
            error('Unknown response setting');
        end

    % ------------------------------------------------------------------- % 
    % Follow-up contrasts:

    elseif strcmp(job.contrastType, 'GoPreferred')
        if strcmp(job.responseSettings, 'Go')

            if strcmp(job.responseSettings, 'rel')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize GoPreferred.                
                data.time2{iSub} = data.ERPdata{iSub, 2}; % initialize GoNonpreferred.

            elseif strcmp(job.outcomeSettings, 'abs')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize GoReward. 
                data.time2{iSub} = data.ERPdata{iSub, 3}; % initialize GoPunishment.

            elseif strcmp(job.outcomeSettings, 'all')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize GoPreferred.
                data.time1{iSub}.avg = real(data.ERPdata{iSub, 1}.avg + data.ERPdata{iSub, 3}.avg) ./ 2; % Mean of GoReward and GoNoPunishment. 
                data.time2{iSub} = data.ERPdata{iSub, 2}; % initialize GoNonpreferred.
                data.time2{iSub}.avg = real(data.ERPdata{iSub, 2}.avg + data.ERPdata{iSub, 4}.avg) ./ 2; % Mean of GoNoReward and GoPunishment.

            else
                error('Unknown outcome setting');
            end                
        elseif strcmp(job.responseSettings, 'Hand')
            error('Contrast not completed yet');
        else
            error('Unknown response setting');
        end

    elseif strcmp(job.contrastType, 'SalientPreferred')
        if strcmp(job.responseSettings, 'Go')

            if strcmp(job.responseSettings, 'rel')
                error('Contrast %s not implemented for outcome setting %s', job.contrastType, job.outcomeSettings);
            elseif strcmp(job.outcomeSettings, 'abs')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize true Reward.
                data.time1{iSub}.avg = real(data.ERPdata{iSub, 1}.avg + data.ERPdata{iSub, 4}.avg) ./ 2; % Mean of GoReward and NoGoReward.
                data.time2{iSub} = data.ERPdata{iSub, 3}; % initialize true Punishment.
                data.time2{iSub}.avg = real(data.ERPdata{iSub, 3}.avg + data.ERPdata{iSub, 6}.avg) ./ 2; % Mean of GoPunishment and NoGoPunishment.
            elseif strcmp(job.outcomeSettings, 'all')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize true Reward.
                data.time1{iSub}.avg = real(data.ERPdata{iSub, 1}.avg + data.ERPdata{iSub, 5}.avg) ./ 2; % Mean of GoReward and NoGoReward.
                data.time2{iSub} = data.ERPdata{iSub, 4}; % initialize true Punishment.
                data.time2{iSub}.avg = real(data.ERPdata{iSub, 4}.avg + data.ERPdata{iSub, 8}.avg) ./ 2; % Mean of GoPunishment and NoGoPunishment.
            else
                error('Unknown outcome setting');
            end                
        elseif strcmp(job.responseSettings, 'Hand')
            error('Contrast not completed yet');
        else
            error('Unknown response setting');
        end

    elseif strcmp(job.contrastType, 'SalientAction')
        if strcmp(job.responseSettings, 'Go')

            if strcmp(job.responseSettings, 'rel')
                
                error('Contrast %s not implemented for outcome setting %s', job.contrastType, job.outcomeSettings);

            elseif strcmp(job.outcomeSettings, 'abs')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize Go.
                data.time1{iSub}.avg = real(data.ERPdata{iSub, 1}.avg + data.ERPdata{iSub, 3}.avg) ./ 2; % Mean of GoReward and Go Punishment.
                data.time2{iSub} = data.ERPdata{iSub, 4}; % initialize NoGo.
                data.time2{iSub}.avg = real(data.ERPdata{iSub, 4}.avg + data.ERPdata{iSub, 6}.avg) ./ 2; % Mean of NoGoReward and NoGoPunishment.

            elseif strcmp(job.outcomeSettings, 'all')
                data.time1{iSub} = data.ERPdata{iSub, 1}; % initialize Go.
                data.time1{iSub}.avg = real(data.ERPdata{iSub, 1}.avg + data.ERPdata{iSub, 4}.avg) ./ 2; % Mean of GoReward and Go Punishment.
                data.time2{iSub} = data.ERPdata{iSub, 5}; % initialize NoGo.
                data.time2{iSub}.avg = real(data.ERPdata{iSub, 5}.avg + data.ERPdata{iSub, 8}.avg) ./ 2; % Mean of NoGoReward and NoGoPunishment.

            else
                error('Unknown outcome setting');
            end                
        elseif strcmp(job.responseSettings, 'Hand')
            error('Contrast not completed yet');
        else
            error('Unknown response setting');
        end

    else
        error('Unknown contrast setting');
    end
end

% ----------------------------------------------------------------------- %
%% Compute objects WITHOUT LOOPING over subjects

if strcmp(job.contrastType, 'BOLD')

    if strcmp(job.ROI2use, 'GLM1StriatumConj')
        job.ROIname = 'striatal';
    elseif strcmp(job.ROI2use, 'GLM1ACCConjMan')
        job.ROIname = 'ACC';
    elseif strcmp(job.ROI2use, 'GLM1vmPFCConjMan')
        job.ROIname = 'vmPFC';
    elseif strcmp(job.ROI2use, 'GLM1PCCConj')
        job.ROIname = 'PCC';
    elseif strcmp(job.ROI2use, 'GLM1V1Conj')
        job.ROIname = 'V1';
    else
        job.ROIname = 'High vs. low BOLD';        
        warning('Unknown region, use default labels');
    end

    % a) 2-line plot
    data.mat1 = squeeze(data.SubCondTime(job.validSubs, 1, :)); % low BOLD.
    data.mat2 = squeeze(data.SubCondTime(job.validSubs, 2, :)); % high BOLD.
    % b) Topoplot: high - low
    data.topo2plot.avg = data.Cond{2}.avg - data.Cond{1}.avg;

elseif strcmp(job.contrastType, 'Preferred')

    if strcmp(job.responseSettings, 'Go') 

        if strcmp(job.outcomeSettings, 'rel') 
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((data.SubCondTime(job.validSubs, 1, :) + data.SubCondTime(job.validSubs, 3, :)) ./ 2); % preferred
            data.mat2 = squeeze((data.SubCondTime(job.validSubs, 2, :) + data.SubCondTime(job.validSubs, 4, :)) ./ 2); % non-preferred
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{3}.avg) ./ 2 - (...
                data.Cond{2}.avg + ...
                data.Cond{4}.avg) ./ 2; % nonpreferred - preferred

        elseif strcmp(job.outcomeSettings, 'abs') 
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((data.SubCondTime(job.validSubs, 1, :) + data.SubCondTime(job.validSubs, 4, :)) ./ 2); % preferred
            data.mat2 = squeeze((data.SubCondTime(job.validSubs, 3, :) + data.SubCondTime(job.validSubs, 6, :)) ./ 2); % non-preferred
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{4}.avg) ./ 2 - (...
                data.Cond{3}.avg + ...
                data.Cond{6}.avg) ./ 2; % nonpreferred - preferred

        elseif strcmp(job.outcomeSettings, 'all') 
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((...
                data.SubCondTime(job.validSubs, 1, :) + ...
                data.SubCondTime(job.validSubs, 3, :) + ...
                data.SubCondTime(job.validSubs, 5, :) + ...
                data.SubCondTime(job.validSubs, 7, :)) ./ 4); % preferred
            data.mat2 = squeeze((...
                data.SubCondTime(job.validSubs, 2, :) + ...
                data.SubCondTime(job.validSubs, 4, :) + ...
                data.SubCondTime(job.validSubs, 6, :) + ...
                data.SubCondTime(job.validSubs, 8, :)) ./ 4); % nonpreferred
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{3}.avg + ...
                data.Cond{5}.avg + ...
                data.Cond{7}.avg) ./ 4 - (...
                data.Cond{2}.avg + ...
                data.Cond{4}.avg + ...
                data.Cond{6}.avg + ...
                data.Cond{8}.avg) ./ 4; % preferred - nonpreferred

        else
            error('Unknown outcome setting');
        end            
    elseif strcmp(job.responseSettings, 'Hand')
        if strcmp(job.outcomeSettings, 'rel')
            
            error('Contrast not completed yet');
            
        elseif strcmp(job.outcomeSettings, 'all') 
            
           error('Contrast not completed yet');
        else
            error('Unknown outcome setting');
        end
    else
        error('Inknown response setting');
    end
    
elseif strcmp(job.contrastType, 'Action')

    if strcmp(job.responseSettings, 'Go')

        if strcmp(job.outcomeSettings, 'rel')
            % a) 2-line plot
            data.mat1 = squeeze((data.SubCondTime(job.validSubs, 1, :) + data.SubCondTime(job.validSubs, 2, :)) ./ 2); % Go.
            data.mat2 = squeeze((data.SubCondTime(job.validSubs, 3, :) + data.SubCondTime(job.validSubs, 4, :)) ./ 2); % NoGo.
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{2}.avg) ./ 2 - (...
                data.Cond{3}.avg + ...
                data.Cond{4}.avg) ./ 2;

        elseif strcmp(job.outcomeSettings, 'abs')
            % a) 2-line plot
            data.mat1 = squeeze((...
                data.SubCondTime(job.validSubs, 1, :) + ...
                data.SubCondTime(job.validSubs, 2, :) + ...
                data.SubCondTime(job.validSubs, 3, :)) ./ 3); % Go.
            data.mat2 = squeeze((...
                data.SubCondTime(job.validSubs, 4, :) + ...
                data.SubCondTime(job.validSubs, 5, :) + ...
                data.SubCondTime(job.validSubs, 6, :)) ./ 3); % NoGo.
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{2}.avg + ...
                data.Cond{3}.avg) ./ 3 - (...
                data.Cond{4}.avg + ...
                data.Cond{5}.avg + ...
                data.Cond{6}.avg) ./ 3;

    elseif strcmp(job.outcomeSettings, 'all')
            % a) 2-line plot
            data.mat1 = squeeze((...
                data.SubCondTime(job.validSubs, 1, :) + ...
                data.SubCondTime(job.validSubs, 2, :) + ...
                data.SubCondTime(job.validSubs, 3, :) + ...
                data.SubCondTime(job.validSubs, 4, :)) ./ 4); % Go.
            data.mat2 = squeeze((...
                data.SubCondTime(job.validSubs, 5, :) + ...
                data.SubCondTime(job.validSubs, 6, :) + ...
                data.SubCondTime(job.validSubs, 7, :) + ...
                data.SubCondTime(job.validSubs, 8, :)) ./ 4); % NoGo.
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{2}.avg + ...
                data.Cond{3}.avg + ...
                data.Cond{4}.avg) ./ 4 - (...
                data.Cond{5}.avg + ...
                data.Cond{6}.avg + ...
                data.Cond{7}.avg + ...
                data.Cond{8}.avg) ./ 4;
        else
            error('Unknown outcome setting')
        end

    elseif strcmp(job.responseSettings, 'Hand')

        if strcmp(job.outcomeSettings, 'rel')
            % a) 2-line plot:
            data.mat1 = squeeze((...
                data.SubCondTime(job.validSubs, 1, :) + ...
                data.SubCondTime(job.validSubs, 2, :) + ...
                data.SubCondTime(job.validSubs, 3, :) + ...
                data.SubCondTime(job.validSubs, 4, :)) ./ 4); 
            data.mat2 = squeeze((...
                data.SubCondTime(job.validSubs, 5, :) + ...
                data.SubCondTime(job.validSubs, 6, :)) ./ 2); 
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{2}.avg + ...
                data.Cond{3}.avg + ...
                data.Cond{4}.avg) ./ 4 - (...
                data.Cond{5}.avg + ...
                data.Cond{6}.avg) ./ 2;

        elseif strcmp(job.outcomeSettings, 'abs')
            % a) 2-line plot:
            data.mat1 = squeeze((...
                data.SubCondTime(job.validSubs, 1, :) + ...
                data.SubCondTime(job.validSubs, 2, :) + ...
                data.SubCondTime(job.validSubs, 3, :) + ...
                data.SubCondTime(job.validSubs, 4, :) + ...
                data.SubCondTime(job.validSubs, 5, :) + ...
                data.SubCondTime(job.validSubs, 6, :)) ./ 6); 
            data.mat2 = squeeze((...
                data.SubCondTime(job.validSubs, 7, :) + ...
                data.SubCondTime(job.validSubs, 8, :) + ...
                data.SubCondTime(job.validSubs, 9, :)) ./ 3); 
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{2}.avg + ...
                data.Cond{3}.avg + ...
                data.Cond{4}.avg + ...
                data.Cond{5}.avg + ...
                data.Cond{6}.avg) ./ 6 - (...
                data.Cond{7}.avg + ...
                data.Cond{8}.avg + ...
                data.Cond{9}.avg) ./ 3;

        elseif strcmp(job.outcomeSettings, 'all')
            % a) 2-line plot:
            data.mat1 = squeeze((...
                data.SubCondTime(job.validSubs, 1, :) + ...
                data.SubCondTime(job.validSubs, 2, :) + ...
                data.SubCondTime(job.validSubs, 3, :) + ...
                data.SubCondTime(job.validSubs, 4, :) + ...
                data.SubCondTime(job.validSubs, 5, :) + ...
                data.SubCondTime(job.validSubs, 6, :) + ...
                data.SubCondTime(job.validSubs, 7, :) + ...
                data.SubCondTime(job.validSubs, 8, :)) ./ 8); 
            data.mat2 = squeeze((...
                data.SubCondTime(job.validSubs, 9, :) + ...
                data.SubCondTime(job.validSubs, 10, :) + ...
                data.SubCondTime(job.validSubs, 11, :) + ...
                data.SubCondTime(job.validSubs, 12, :)) ./ 4); 
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{2}.avg + ...
                data.Cond{3}.avg + ...
                data.Cond{4}.avg + ...
                data.Cond{5}.avg + ...
                data.Cond{6}.avg + ...
                data.Cond{7}.avg + ...
                data.Cond{8}.avg) ./ 8 - (...
                data.Cond{9}.avg + ...
                data.Cond{10}.avg + ...
                data.Cond{11}.avg + ...
                data.Cond{12}.avg) ./ 4;

        else
            error('Unknown outcome setting');
        end
    else
        error('Unknown response setting');
    end
    
% ----------------------------------------------------------------------- %
% Follow-up contrasts:

elseif strcmp(job.contrastType, 'GoPreferred')
    if strcmp(job.responseSettings, 'Go')

        if strcmp(job.outcomeSettings, 'rel')

            % a) 2-line plot
            data.mat1 = squeeze(data.SubCondTime(job.validSubs, 1, :)); % preferred Go.
            data.mat2 = squeeze(data.SubCondTime(job.validSubs, 2, :)); % nonpreferred Go.
            % b) Topoplot:
            data.topo2plot.avg = data.Cond{1}.avg - data.Cond{2}.avg;

        elseif strcmp(job.outcomeSettings, 'abs')

            % a) 2-line plot
            data.mat1 = squeeze(data.SubCondTime(job.validSubs, 1, :)); % GoReward.
            data.mat2 = squeeze(data.SubCondTime(job.validSubs, 3, :)); % GoPunishment.
            % b) Topoplot:
            data.topo2plot.avg = data.Cond{1}.avg - data.Cond{3}.avg;

        elseif strcmp(job.outcomeSettings, 'all')

            % a) 2-line plot
            data.mat1 = squeeze((data.SubCondTime(job.validSubs, 1, :) + data.SubCondTime(job.validSubs, 3, :)) ./ 2); % preferred Go.
            data.mat2 = squeeze((data.SubCondTime(job.validSubs, 2, :) + data.SubCondTime(job.validSubs, 4, :)) ./ 2); % non-preferred Go.
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{3}.avg) ./ 2 - (...
                data.Cond{2}.avg + ...
                data.Cond{4}.avg) ./ 2;

        else
            error('Unknown outcome setting');
        end
    elseif strcmp(job.responseSettings, 'Hand')
        error('Contrast GoPreferred under response setting Hand not yet implemented');
    else
        error('Unknown response setting');
    end

elseif strcmp(job.contrastType, 'SalientPreferred')

    if strcmp(job.responseSettings, 'Go')

        if strcmp(job.outcomeSettings, 'rel')
            error('No contrast %s under response setting %s', job.contrastType, job.responseSettings);

        elseif strcmp(job.outcomeSettings, 'abs')
            % a) 2-line plot
            data.mat1 = squeeze((data.SubCondTime(job.validSubs, 1, :) + data.SubCondTime(job.validSubs, 4, :)) ./ 2); % true reward.
            data.mat2 = squeeze((data.SubCondTime(job.validSubs, 3, :) + data.SubCondTime(job.validSubs, 6, :)) ./ 2); % true punishment.
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{4}.avg) ./ 2 - (...
                data.Cond{3}.avg + ...
                data.Cond{6}.avg) ./ 2;

        elseif strcmp(job.outcomeSettings, 'all')
            % a) 2-line plot
            data.mat1 = squeeze((data.SubCondTime(job.validSubs, 1, :) + data.SubCondTime(job.validSubs, 5, :)) ./ 2); % true reward.
            data.mat2 = squeeze((data.SubCondTime(job.validSubs, 4, :) + data.SubCondTime(job.validSubs, 8, :)) ./ 2); % true punishment.
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{5}.avg) ./ 2 - (...
                data.Cond{4}.avg + ...
                data.Cond{8}.avg) ./ 2;
        else
            error('Unknown outcome setting');
        end

    elseif strcmp(job.responseSettings, 'Hand')
            error('No contrast %s under response setting %s', job.contrastType, job.responseSettings);
    else
        error('Unknown response setting');
    end

elseif strcmp(job.contrastType, 'SalientAction')

    if strcmp(job.responseSettings, 'Go')

        if strcmp(job.outcomeSettings, 'rel')
            error('No contrast %s under response setting %s', job.contrastType, job.responseSettings);

        elseif strcmp(job.outcomeSettings, 'abs')
            % a) 2-line plot:
            data.mat1 = squeeze((data.SubCondTime(job.validSubs, 1, :) + data.SubCondTime(job.validSubs, 3, :)) ./ 2); % Go.
            data.mat2 = squeeze((data.SubCondTime(job.validSubs, 4, :) + data.SubCondTime(job.validSubs, 6, :)) ./ 2); % NoGo.
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{3}.avg) ./ 2 - (...
                data.Cond{4}.avg + ...
                data.Cond{6}.avg) ./ 2;

        elseif strcmp(job.outcomeSettings, 'all')
            % a) 2-line plot:
            data.mat1 = squeeze((data.SubCondTime(job.validSubs, 1, :) + data.SubCondTime(job.validSubs, 4, :)) ./ 2); % Go.
            data.mat2 = squeeze((data.SubCondTime(job.validSubs, 5, :) + data.SubCondTime(job.validSubs, 8, :)) ./ 2); % NoGo.
            % b) Topoplot:
            data.topo2plot.avg = (...
                data.Cond{1}.avg + ...
                data.Cond{4}.avg) ./ 2 - (...
                data.Cond{5}.avg + ...
                data.Cond{8}.avg) ./ 2;
        else
            error('Unknown outcome setting');
        end
    elseif strcmp(job.responseSettings, 'Hand')
            error('No contrast %s under response setting %s', job.contrastType, job.responseSettings);
    else
        error('Unknown response setting');
    end

else
    error('Unknown contrast type');
end

% ------------------------------------------------------------------ %
%% Error bounds for 2-line plots:

data.matDiff        = data.mat1 - data.mat2; % for difference-CIs
data.matMean        = (data.mat1 + data.mat2) ./ 2; % for Cousineau method
data.matGrandMean   = nanmean(data.matMean); % for Cousineau method
data.SubTime        = squeeze(nanmean(data.SubCondTime, 2)); % average over conditions for Cousineau method
data.GrandTime      = squeeze(nanmean(data.SubTime(job.validSubs, :), 1)); % average over subjects for Cousineau method

fprintf('Finished preparing data for contrast %s, channels %s\n', job.contrastType, strjoin(job.channels, '/'));

end % END OF FUNCTION.