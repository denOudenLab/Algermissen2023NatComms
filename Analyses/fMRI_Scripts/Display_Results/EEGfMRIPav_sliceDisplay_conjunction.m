function EEGfMRIPav_sliceDisplay_conjunction(dirs, job)

% EEGfMRIPav_sliceDisplay_conjunction(dirs, job)
%
% Plots conjunction of 2 contrasts at multiple cluster-forming thresholds
% using slice display (load files into layers using different file name
% format). Code zstat with both color and opacity.
%
% Before using: mind to set directories/paths in slice_display_set_dirs to respective toolboxes needed!
%
% INPUTS:
% dirs              = structure with various directories (none required,
% defaults set automatically).
% job               = structure with various settings:
% NECESSARY INPUTS:
% .lockSettings     = string, either 'stimlocked' or 'outcomelocked'.
% .type             = string, either 'standard' or 'conjunction'.
% .iView            = string, either 'sagittal' or 'coronal' or 'axial'.
% .iSlice           = scalar integer, value of slice to display.
%
% OPTIONAL INPUTS:
% .type             = string, either 'standard' or 'conjunction'.
% .GLMID            = string, name of GLM.
% .Cope             = number of cope under type 'standard'.
% .firstCope        = scalar integer, number of first cope under type
% 'conjunction'.
% .secondCope       = scalar integer, number of second cope under type
% 'conjunction'.
% .sign             = string, either 'pos' or 'neg'.
% .cLim             = vector of two floats, c-value range (color
% intensity).
% .zLim             = vector of two floats, z-value range (opacity).
%
% OUTPUTS:
% generate plot, save if specified.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% ----------------------------------------------------------------------- %
%% 1) Set directory and job:

dirs    = slice_display_set_dirs(dirs); % complete directories
job     = slice_display_set_job(job); % complete job settings

% Load color maps:
fprintf('Load slice display color maps\n');
load(fullfile(dirs.sliceDisplay, 'colormaps.mat')); % Get custom colormaps from slice display

% ----------------------------------------------------------------------- %
%% 2) Determine data directory:

% Directory with cope/zstat/thresh files:
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined', job.GLMID));

% ----------------------------------------------------------------------- %
%% 3) Create plot:

fprintf('Create plot for GLM %s, contrasts %d and %d, %s view, slice %d \n', ...
    job.GLMID, job.firstCope, job.secondCope, job.iView, job.iSlice);
            
% ----------------------------------------------------------------------- %
% Step 1: Initialize empty layers and settings variables:
layers                              = sd_config_layers('init', {'truecolor', 'dual', 'contour'});
settings                            = sd_config_settings('init');

% ----------------------------------------------------------------------- %
% Step 2: Define layers:

% ----------------------------------------------------------------------- %
% Layer 1: Anatomical map:

layers(1).color.file                = fullfile(dirs.sliceDisplay, 'MNI152_T1_2mm_brain.nii');
layers(1).color.map                 = gray(256);

% ----------------------------------------------------------------------- %
% Layer 2: Dual-coded layer (contrast estimates color-coded; z-statistics opacity-coded):

% Layer 2 color settings:
% Load zstat combining two copes (conjunction).
layers(2).color.file                = fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s.nii', ...
    job.firstCope, job.secondCope, job.sign)); % t-map/ z-map input: zstat1.nii
% layers(2).color.map                 = CyBuGyRdYl; % one of the custom color-maps in slice_display
layers(2).color.map                 = hot(256); % one of the custom color-maps in slice_display

layers(2).color.label               = 'Conjunction'; % '\beta_{left} - \beta_{right} (a.u.)'; % title for legend
layers(2).color.range               = job.zLim; % y-axis range (opacity) of legend (z/t-values)

% --------------------------------------------- %
% Layer 2 opacity settings:
% Load zstat combining two copes (conjunction).
layers(2).opacity.file              = fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s_bin21.nii', ...
    job.firstCope, job.secondCope, job.sign)); % t-map/ z-map input: zstat1.nii
layers(2).opacity.label             = '| z |'; % y-axis label (opacity) of legend
layers(2).opacity.range             = job.cLim; % y-axis range (opacity) of legend (z/t-values)

% ----------------------------------------------------------------------- %
% Layer 3: Contour of significantly activated voxels:
% Thresholded map: thresh_zstat1.nii.gz.

layers(3).color.file                = fullfile(dirs.data, sprintf('thresh_zstat1_copes_%d_%d_%s_bin31.nii', ...
    job.firstCope, job.secondCope, job.sign)); % t-map/ z-map input: zstat1.nii
layers(3).color.map                 = [0 0 0]; % white is [1 1 1]; gray is [0.3 0.3 0.3]; # black is [0 0 0];
layers(3).color.line_width          = 2; % default 3
layers(3).color.line_style          = '-';

% ----------------------------------------------------------------------- %
% Image settings:

settings.fig_specs.n.slice_column   = 1; % # columms for display
settings.fig_specs.title            = 'Conjunction'; % 'left - right'; % caption
settings.slice.orientation          = char(job.iView); % axial sagittal coronal
settings.slice.disp_slices          = job.iSlice; % single slice

% ----------------------------------------------------------------------- %
% Display the layers:

[settings, p]   = sd_display(layers, settings);
frame_h         = get(handle(gcf), 'JavaFrame');
            
%% 4) Save figure:

fileName = sprintf('SD_GLM%s_Cope%02d_%02d_%s_slice%d_z_%d_%d', job.GLMID, job.firstCope, job.secondCope, ...
    settings.slice.orientation, int8(settings.slice.disp_slices(1)), job.zLim(1)*10, job.zLim(end)*10);
if job.isSave
    saveas(gcf, fullfile(dirs.save, [fileName '.png']));
end
% pause(1);
% close gcf

end % END OF FUNCTION.
