% FigureS19.m

% Plots for Supplementary Figure S19.
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

% Add paths to helper files:
addpath(fullfile(dirs.root, 'Analyses/EEG_Scripts/OutcomeLockedAnalyses/TAfT'));

% Add Fieldtrip:
addpath /home/common/matlab/fieldtrip;
ft_defaults;

% ----------------------------------------------------------------------- %
%% Load saved regression coefficients:

EEGdomain           = 'TF';

% 1) ROIs:
ROIs2use            = {'GLM1StriatumConj', 'GLM1ACCConjMan', ...
    'GLM1LeftMotorConj', 'GLM1vmPFCConjMan', 'GLM1PCCConj', ...
    'GLM1LeftITGConj', 'GLM1V1Conj'}; 
if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end

% 2) Behavioral regressors:
behav2use           = {'Updatestd', 'Updatedif'}; % Updates + outcome sign

% 3) Perform only on selected trials:
selTrials           = 'all';

[job, dirs, betas]  = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
%% Figure S19A-D: TF plots:

% NOTE: Within within taft_postprocess_TF_TFplot, set latency to
% cfg.latency = [-0.25 0.8]; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)); nContour = 10; % default 10

% Exclude subjects:
job.invalidSubs     = [11 12 15 23 25 26 30]; % bad co-registrations & outliers for outcome-locked pre-processing & inspection (add 26)

% Selected channels:
selChans            = {'Fz', 'FCz', 'Cz'};

% Loop over ROIs:
for iROI = 1:7
        
    % Set seed:
    rng(20190822) % set random number generator for constant p-values

    % Preprocess betas:
    [sortBetas,~]   = taft_postprocess_TF_selectData(job, betas, iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:

    % Plot to retrieve t-values and p-values:
    [tg, corrp]     = taft_postprocess_TF_TFplot(job, dirs, sortBetas, iROI, selChans, 2, 1000, 2, false); % iROI, selChans, thresh, nP
    close gcf

    % Settings:
    timeIdx         = size(tg, 4); % maximum time
    zlim            = 3;
    lineWidth       = 3;
    fontSize        = 32; % 44; % uni: 32; smaller: 24 (3/4)
    pCrit           = 0.05;

    % Exception in pCrit because cluster not-significant when tested
    % -0.25 - 0.8 sec. across all frequencies:
    if contains(job.regNames(iROI), 'PCCConj'); pCrit = .20; end

    % ------------------------------------------------------------------- %
    % Start figure:
    figure('Position', [100 100 1000 800], 'color', 'white'); hold on % narrower for Figure 5 TAfT

    % Make contour plot:
    contourf(sortBetas{1}.time(1:timeIdx), sortBetas{1}.freq, squeeze(real(tg)), ...
        nContour, 'linestyle', 'none');

    % ------------------------------------------------------------------- %
    % Add vertical lines:

    plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', 3);
    plot([0.7 0.7], get(gca, 'ylim'), ':k', 'LineWidth', 3);

    % ------------------------------------------------------------------- %
    % Settings:

    set(gca, 'xlim', [-0.25 0.8], 'ylim', [1.5 33], 'clim', [-1*zlim  1*zlim], ...
        'xtick', [-1 -0.5 -0.25 0 0.25 0.5 0.7 1], ...
        'xtickLabel', {'-1', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
        'ytick', [2 4 8 16 32], 'yscale', 'log', ... 
        'fontsize', fontSize, 'Linewidth', lineWidth);

    % Labels:
    xlabel('Time (ms)', 'FontSize', fontSize, 'FontName', 'Arial', 'fontweight', 'bold');
    ylabel('Frequency (Hz)', 'FontSize', fontSize, 'FontName', 'Arial', 'fontweight', 'bold');

    % Color bar:
    colorbar('Ticks', (-1*zlim):(zlim/2):zlim, 'Fontsize', fontSize);

    % ------------------------------------------------------------------- %
    % Add highlights based on p-values:

    if sum(corrp(:) < pCrit) > 0 % if anything significant
       isSig = squeeze(double(corrp < pCrit)); % map significant clusters
       % see https://nl.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a#answer_211204

       pause(1) % pause to allow for plotting both plots into one:           
       hold on % on top of old plot

       % Extra contour:
       [~, hContour]  = contourf(sortBetas{1}.time(1:timeIdx), sortBetas{1}.freq, squeeze(isSig), 1);
       hContour.LineWidth = 5;
       drawnow;  % this is important, to ensure that FacePrims is ready in the next line!

       hFills = hContour.FacePrims;  % array of TriangleStrip objects
       [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
       for idx = 1:numel(hFills)
           hFills(idx).ColorData(4) = 1;   % default=255
           pause(1)
       end

       hold off

    end

    % ------------------------------------------------------------------- %
    % Save as:
    fig = gcf;
    if strcmp(ROIs2use{iROI}, 'GLM1LeftMotorConj')
        figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS19A');
        saveas(gcf, [figName '.png']);
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS19A.csv');
        csvwrite(fullFileName, squeeze(real(tg)));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1vmPFCConjMan')
        figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS19B');
        saveas(gcf, [figName '.png']);
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS19B.csv');
        csvwrite(fullFileName, squeeze(real(tg)));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1LeftITGConj')
        figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS19C');
        saveas(gcf, [figName '.png']);
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS19C.csv');
        csvwrite(fullFileName, squeeze(real(tg)));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1V1Conj')
        figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS19D');
        saveas(gcf, [figName '.png']);
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS19D.csv');
        csvwrite(fullFileName, squeeze(real(tg)));
    end
    close gcf

end % end iROI
fprintf('Done :-)\n');

% ----------------------------------------------------------------------- %
%% Figures S19E: Multiple topoplots: Left motor cortex 13-30 Hz

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64));

% General settings:
zlim            = 4;
nCols           = 4;
lineWidth       = 3;

% Timing settings:
startTime       = 0.0;
endTime         = 0.7;
steps           = 0.1;
startTimeVec    = startTime:steps:endTime;
endTimeVec      = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);
nPlot           = length(endTimeVec);
nRows           = ceil(nPlot/nCols); % number rows
saveMat         = nan(64, nPlot);

% Frequency settings:
startFreq       = 13;
endFreq         = 30;

% Select ROIs:
iROI            = 3;
                    
% Prepare data:
[~, Tvalues]    = taft_postprocess_TF_selectData(job, betas, iROI);

% Start plot:
figure('Position', [100 100 1200 800], 'Color', 'white'); % default
for iPlot = 1:length(endTimeVec)

    subplot(nRows, nCols, iPlot);

    % Set config object for Fieldtrip:
    cfg         = []; 
    cfg.xlim    = [startTimeVec(iPlot) endTimeVec(iPlot)]; 
    cfg.ylim    = [startFreq endFreq]; 
    cfg.zlim    = [-1*zlim 1*zlim]; 
    cfg.marker  = 'on'; 
    cfg.style   = 'fill';
    cfg.layout  = 'easycapM11.mat'; cfg.comment = 'no'; 
    cfg.figure  = gcf;
    cfg.colorbar= 'no'; % want 'no', i.e. do it yourself % --> add externally

    % Plot:
    ft_topoplotTFR(cfg, Tvalues);
    if length(endTimeVec) > 1 % subheading if more than one plot
        title(sprintf('%.1f-%.1f s', startTimeVec(iPlot), endTimeVec(iPlot)), 'fontsize', 28); % 1 digit: 28
    end % end iTime
    
    % Set linewidth:
    p = gca;
    for i = 7:11
        t = p.Children(i); 
        t.LineWidth = lineWidth;
    end
    
    % Save:
    timeIdx             = dsearchn(Tvalues.time', [startTimeVec(iPlot) endTimeVec(iPlot)]'); 
    timeIdx             = timeIdx(1):timeIdx(end);
    freqIdx             = dsearchn(Tvalues.freq', [startFreq endFreq]'); 
    freqIdx             = freqIdx(1):freqIdx(end);
    saveMat(:, iPlot)   = squeeze(mean(mean(Tvalues.powspctrm(:, freqIdx, timeIdx), 3), 2));

end % end iPlot
        
% Save:
figName = fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS19E');
saveas(gcf, [figName '.png'])
pause(2)
close gcf

% Save source data file:
saveMat         = horzcat(Tvalues.label, num2cell(saveMat));
fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS19E.csv'));
writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);

% ----------------------------------------------------------------------- %
%% Figures S19F: Multiple topoplots: Primary visual cortex 8-13 Hz

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64));

% General settings:
zlim            = 4;
nCols           = 4;
lineWidth       = 3;

% Timing settings:
startTime       = 0.0;
endTime         = 0.7;
steps           = 0.1;
startTimeVec    = startTime:steps:endTime;
endTimeVec      = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);
nPlot           = length(endTimeVec);
nRows           = ceil(nPlot/nCols); % number rows
saveMat         = nan(64, nPlot);

% Frequency settings:
startFreq       = 8;
endFreq         = 13;

% Select ROIs:
iROI            = 7;
                    
% Prepare data:
[~, Tvalues]    = taft_postprocess_TF_selectData(job, betas, iROI);

% Start plot:
figure('Position', [100 100 1200 800]); % default
for iPlot = 1:length(endTimeVec)

    subplot(nRows, nCols, iPlot);

    % Set config object for Fieldtrip:
    cfg         = []; 
    cfg.xlim    = [startTimeVec(iPlot) endTimeVec(iPlot)]; 
    cfg.ylim    = [startFreq endFreq]; 
    cfg.zlim    = [-1*zlim 1*zlim]; 
    cfg.marker  = 'on'; 
    cfg.style   = 'fill';
    cfg.layout  = 'easycapM11.mat'; cfg.comment = 'no'; 
    cfg.figure  = gcf;
    cfg.colorbar= 'no'; % want 'no', i.e. do it yourself % --> add externally

    % Plot:
    ft_topoplotTFR(cfg, Tvalues);
    if length(endTimeVec) > 1 % subheading if more than one plot
        title(sprintf('%.1f-%.1f s', startTimeVec(iPlot), endTimeVec(iPlot)), 'fontsize', 28); % 1 digit: 28
    end % end iTime
    
    % Set linewidth:
    p = gca;
    for i = 7:11
        t = p.Children(i); 
        t.LineWidth = lineWidth;
    end
    
    % Save:
    timeIdx             = dsearchn(Tvalues.time', [startTimeVec(iPlot) endTimeVec(iPlot)]'); 
    timeIdx             = timeIdx(1):timeIdx(end);
    freqIdx             = dsearchn(Tvalues.freq', [startFreq endFreq]'); 
    freqIdx             = freqIdx(1):freqIdx(end);
    saveMat(:, iPlot)   = squeeze(mean(mean(Tvalues.powspctrm(:, freqIdx, timeIdx), 3), 2));

end % end iPlot
        
% Save:
figName = fullfile(dirs.root, '/Log/OutcomeLockedPaperPlots/FigS19F');
saveas(gcf, [figName '.png'])
pause(2)
close gcf

% Save source data file:
saveMat         = horzcat(Tvalues.label, num2cell(saveMat));
fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS19F.csv'));
writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);
                
% END OF FILE.