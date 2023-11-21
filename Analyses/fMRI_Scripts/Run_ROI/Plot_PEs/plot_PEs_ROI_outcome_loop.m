%% Be here:

% Overview script to create bar plots for all kinds of ROIs.
% Johannes Algermissen, 2020.
%
% We are here: cd /project/3017042.02/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs

% Overview script for running all kinds of 

% Set directories:
dirs.root       = '/project/3017042.02';
dirs.scripts    = fullfile(dirs.root, 'Analyses/fMRI_Scripts/Run_ROI/Plot_PEs');
addpath(dirs.scripts);

% ======================================================================= %
%% Settings for GLM2: long version

GLMID       = '2';
ROI2use     = {'GLM2vmPFCValence'; 'GLM2StriatumValence'};
ROITitle    = {'vmPFC'; 'Striatum'};

% Manipulate validSubs inside plot_PEs_ROI_outcome_actValSal

% Loop over ROIs:
for iPlot = 1:length(ROI2use)

    % Without points, with legend:
    plot_PEs_ROI_outcome(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, 'actOut', false, false); close gcf;
    % With points, without legend --> Figure 2:
    plot_PEs_ROI_outcome(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, 'actOut', true, false, [-0.4 0.4]); close gcf;
    
    % First grouped by action, then outcome valence next to each other:
    plot_PEs_ROI_outcome_actValSal(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, false, true) % leave yLim open
    plot_PEs_ROI_outcome_actValSal(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, true, false, [-0.35 0.4]); close gcf; % fix yLim for vmPFC and striatum

    % First grouped by action, then outcome valence next to each other:
    plot_PEs_ROI_outcome_salActVal(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, false, true); close gcf; % leave yLim open
    plot_PEs_ROI_outcome_salActVal(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, true, false); close gcf; % leave yLim open
    
    % First grouped by valence, then action, then salience:
    plot_PEs_ROI_outcome_valActSal(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, false, true); close gcf; % leave yLim open
    plot_PEs_ROI_outcome_valActSal(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, true, false); close gcf; % leave yLim open
end

% ======================================================================= %
%% Settings for GLM2: short version

GLMID       = '2';
ROI2use     = {'GLM2vmPFCValence'; 'GLM2StriatumValence'};
ROITitle    = {'vmPFC'; 'Striatum'};

% Manipulate validSubs inside plot_PEs_ROI_outcome_actValSal

% Loop over ROIs:
for iPlot = 1:length(ROI2use)

    % With points, without legend --> Figure 2:
    plot_PEs_ROI_outcome(char(ROI2use(iPlot)), char(ROITitle(iPlot)), GLMID, 'actOut', true, false, [-0.4 0.4]); close gcf;
    
end


% END OF FILE.
