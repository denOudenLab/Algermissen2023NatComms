function job = slice_display_set_job(job)

% job = slice_display_set_job(job)
% 
% Complete any job settings for slice display given job name that might not
% have been specified yet.
% 
% INPUT:
% 
% job               = structure with options (necessary and optional
% settings).
%
% NECESSARY INPUTS:
% .lockSettings     = string, either 'stimlocked' or 'outcomelocked'.
% .iView            = string, either 'sagittal' or 'coronal' or 'axial'.
% .iSlice           = scalar integer, value of slice to display.
%
% OPTIONAL INPUTS:
% .type             = string, either 'standard' or 'conjunction'.
% .GLMID            = string, name of GLM.
% .levels           = scalar integer, number of levels in GLM (depends on
% lockSettings).
% .iCope            = scalar integer, number of cope under type 'standard'.
% .firstCope        = scalar integer, number of first cope under type
% 'conjunction'.
% .secondCope       = scalar integer, number of second cope under type
% 'conjunction'.
% .sign             = string, either 'pos' or 'neg'.
% .cLim             = vector of two floats, c-value range (color
% intensity).
% .zLim             = vector of two floats, z-value range (opacity).
% .isSave           = Boolean, save plot or not.
% 
% OUTPUTS:
% job               = same object with settings filled in.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% ----------------------------------------------------------------------- %
%% 1) Required settings:

if ~isfield(job, 'lockSettings')
    error('No lockSettings specified');
end

if ~isfield(job, 'type')
    error('No type specified');
end

if strcmp(job.type, 'standard')
    if ~isfield(job, 'iCope')
        error('No iCope specified');
    end
end

if ~isfield(job, 'iView')
    error('No iView specified');
end

if ~isfield(job, 'iSlice')
    error('No iSlice specified');
end

% ----------------------------------------------------------------------- %
%% 2) Optional settings: further contrast details

% ------------------------------------------ %
% GLM ID if missing:

if ~isfield(job, 'GLMID')
    if strcmp(job.lockSettings, 'stimlocked')
        error('Default GLMID for this lockSetting not set');
    elseif strcmp(job.lockSettings, 'outcomelocked')
        job.GLMID       = '1'; % 
    else
        error('Unknown lock setting')
    end
    fprintf('No GLMID specified, set default to %s\n', job.GLMID)
end

% ------------------------------------------ %
% Number of levels if missing:

if ~isfield(job, 'levels')
    if contains(job.GLMID, '2') || contains(job.GLMID, '3')
        job.levels       = 2; % 
    elseif contains(job.GLMID, '1')
        job.levels       = 3; % 
    else
        error('Unknown GLM')
    end
    fprintf('No number of levels specified, set default to %s\n', job.levels)
end

% ------------------------------------------ %
% Lock setting if missing:

if strcmp(job.lockSettings, 'outcomelocked') && strcmp(job.type, 'conjunction')

    if ~isfield(job, 'firstCope')
        job.firstCope   = 4; % 
        fprintf('No firstCope specified, set default to %d\n', job.firstCope);
    end

    if ~isfield(job, 'secondCope')
        job.secondCope  = 5; % 
        fprintf('No secondCope specified, set default to %d\n', job.secondCope);
    end
end

% ------------------------------------------ %
% Sign of contrast if missing:
if ~isfield(job, 'sign')
    job.sign       = 'pos'; % 
    fprintf('No sign specified, set default to %s\n', job.sign);
end

% ----------------------------------------------------------------------- %
%% Optional settings: c and z range

% ------------------------------------------ %
% c range if missing:
if ~isfield(job, 'cLim')
    if strcmp(job.lockSettings, 'stimlocked')
        error('cLim for this lockSetting not set');
        
    elseif strcmp(job.lockSettings, 'outcomelocked')
        if strcmp(job.type, 'conjunction')
            job.cLim    = [0 1]; % 
            
        else % standard
            job.cLim    = [0 5]; %   
        end
    else
        error('Unknown lock setting %s', job.lockSettings);
    end
    
    fprintf('No cLim specified, set default to %0.1f-%0.1f\n', job.cLim(1), job.cLim(2));
end

% ------------------------------------------ %
% z range if missing:
if ~isfield(job, 'zLim')
    
    if strcmp(job.lockSettings, 'stimlocked')
        error('zLim for this lockSetting not set');
        
    elseif strcmp(job.lockSettings, 'outcomelocked')
        if strcmp(job.type, 'conjunction')
            job.zLim    = [1 5]; % 
        else
            job.zLim    = [0 5]; %   
        end
    else
        error('Unknown lock setting %s', job.lockSettings);
    end
    
    fprintf('No zLim specified, set default to %0.1f-%0.1f\n', job.zLim(1), job.zLim(2));
    
end

% ----------------------------------------------------------------------- %
%% Save or not:

if ~isfield(job, 'isSave')
    job.isSave = true;
    fprintf('No isSave specified, set default to %s\n', job.isSave);
end

end % END OF FUNCTION.
