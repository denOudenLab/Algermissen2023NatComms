% Figure05.m

% Plots for Figure 5A-G in manuscript.
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
behav2use           = {'Updatestd', 'Updatedif'};

% 3) Perform only on selected trials:
selTrials           = 'all';

[job, dirs, betas]  = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
%% Figure 5A-C: TF plots:

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
    if strcmp(ROIs2use{iROI}, 'GLM1ACCConjMan')
        figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig5A');
        saveas(gcf, [figName '.png']);
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig5A.csv');
        csvwrite(fullFileName, squeeze(real(tg)));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1PCCConj')
        figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig5B');
        saveas(gcf, [figName '.png']);
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig5B.csv');
        csvwrite(fullFileName, squeeze(real(tg)));
    end
    if strcmp(ROIs2use{iROI}, 'GLM1StriatumConj')
        figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig5C');
        saveas(gcf, [figName '.png']);
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig5C.csv');
        csvwrite(fullFileName, squeeze(real(tg)));
    end
    close gcf

end % end iROI
fprintf('Done :-)\n');

% ----------------------------------------------------------------------- %
%% Figure 5A-C: Topoplots:

% Manually added color maps:
addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64)) ; cmap = 'redblue'; nContour = 10; % default 10

% Color range:
zlim = 3;

% ------------------------------------------------- %
% Select ROI:

% Fig 5A:
% iROI = 2; iTime = [0.100 0.575]; iFreq = [7 9]; saveLetter = 'A'; % ACC

% Fig 5B:
% iROI = 5; iTime = [0.2 0.5]; iFreq = [2 8]; saveLetter = 'B';  % PCC

% Fig 5C:
% iROI = 1; iTime = [0.100 0.800]; iFreq = [16 23]; saveLetter = 'C';  % Striatum

% ------------------------------------------------- %
% Compute T-values across subjects:

[sortBetas, Tvalues] = taft_postprocess_TF_selectData(job, betas, iROI);

% ------------------------------------------------- %
% Plot:
figure('units', 'normalized', 'outerposition', [0 0 1 1]); hold on % fullscreen--no risk of implicitly resizing any axes
cfg         = []; 
cfg.figure  = gcf; 
cfg.ylim    = iFreq; 
cfg.zlim    = [-1*zlim 1*zlim]; 
cfg.marker  ='on'; 
cfg.layout  = 'easycapM11.mat'; 
cfg.comment = 'no'; 
cfg.xlim    = [iTime(1) iTime(end)];
cfg.style   = 'fill';
ft_topoplotTFR(cfg, Tvalues);
title(sprintf('%s, %.2f to %.2f sec', job.regNames{iROI}, iTime(1), iTime(end)), 'fontsize', 32);

% Set linewidth:
p = gca;
for i = 7:11
    t = p.Children(i); 
    t.LineWidth = lineWidth;
end

% ------------------------------------------------- %
% Save:
figName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5%s_topo.png', saveLetter));
saveas(gcf, figName);

% Close:
pause(2);
close gcf

% Save source data file:
fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5%s_topo.csv', saveLetter));
timeIdx         = dsearchn(Tvalues.time', iTime'); timeIdx = timeIdx(1):timeIdx(end);
freqIdx         = dsearchn(Tvalues.freq', iFreq'); freqIdx = freqIdx(1):freqIdx(end);
saveMat         = squeeze(mean(mean(Tvalues.powspctrm(:, freqIdx, timeIdx), 3), 2));
saveMat         = horzcat(Tvalues.label, num2cell(saveMat));
writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);

% ----------------------------------------------------------------------- %
%% Plot 5D: Line plot:

% ----------------------------------------------------------------------- %
% Compute sigLine:

% NOTE: Within within taft_postprocess_TF_TFplot, set latency to
% cfg.latency = [-0.25 0.8]; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Initialize:
sigLine = zeros(size(betas{1}.powspctrm, 1), 43); % hard-coded based on latency

for iROI = [1 2 5] % iROI = 5;
        
    selChans        = {'Fz', 'FCz', 'Cz'};
    [sortBetas,~]   = taft_postprocess_TF_selectData(job, betas, iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
    [tg, corrp]     = taft_postprocess_TF_TFplot(job, dirs, sortBetas, iROI, selChans, 2, 1000, 2, false);
    close all
    
    % p-value settings:
    pCrit       = 0.05;
    
    % Adjust threshold because significant only in testing window (0-0.7 sec.):
    if iROI == 2; pCrit = .10; end % lower p-threshold because first beta cluster otherwise not visible
    if iROI == 5; pCrit = .25; end % lower p-threshold because intermediate delta cluster otherwise not visible
    
    % Extract absolute sum of t-values above threshold:
    sigLine(iROI, :) = abs(mean(squeeze(tg .* double(corrp < pCrit)), 1));
    if iROI == 5; sigLine(iROI, 32:end) = 0; end % crop late theta cluster likely not outcome-induced any more

    % Normalize by peak: 
    sigLine(iROI, :) = sigLine(iROI, :) / max(sigLine(iROI, :));
    
end
fprintf('Finished :-)\n');

% ----------------------------------------------------------------------- %
% Plot sigLine:

% Fixed settings:
lineWidth   = 8;
fontSize    = 32; % 44
legendNames = {'ACC', 'PCC', 'Striatum'};
colMat      = [0.75 0.75 0; 1 0 0; 0 0 0; 0 0 0; 0 0.75 0.75];

% Timing:
startTime   = -0.25;
stopTime    = 0.8;
startIdx    = find(round(sortBetas{1}.time, 3) == startTime);
stopIdx     = find(round(sortBetas{1}.time, 3) == stopTime);

% Start figure:
close all
figure('Position', [100 100 1500 800], 'Color', 'white'); hold on
iCond = 0;
p   = cell(3, 1);
for iROI = [2 5 1] % iROI = 5;
    iCond = iCond + 1;
    p{iCond} = plot(sortBetas{1}.time(startIdx:stopIdx), sigLine(iROI, startIdx:stopIdx)*100, ...
        'color', colMat(iROI, :), 'linewidth', lineWidth);
end

% Add vertical line at zero:
plot(sortBetas{1}.time, zeros(length(sortBetas{1}.time)), 'k-', 'linewidth', lineWidth);

% Add horizontal lines:
plot([0 0], get(gca, 'ylim'), ':k', 'LineWidth', 3);
plot([0.7 0.7], get(gca, 'ylim'), ':k', 'LineWidth', 3);

% Settings:
set(gca, 'xlim', [-0.25 0.8], ...
    'xtick', [0 0.25 0.5 0.7 1], ...
    'xtickLabel', {'FbOn', '250', '500', 'FbOff', '1,000'}, ... 
    'fontsize', fontSize, 'Linewidth', lineWidth);
xlabel('Time (ms)', 'FontSize', fontSize, 'FontName', 'Arial', 'fontweight', 'bold');
ylabel('% of maximum', 'fontweight', 'bold', 'fontsize', fontSize); 

% Save figure:
saveas(gcf, fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig5D.png')); % automatic ROI name
pause(3); 
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig5D.csv');
csvwrite(fullFileName, sigLine([2 5 1], startIdx:stopIdx)'*100);
 
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Figure 5E-G: Using EEG as regressors for fMRI:

addpath(fullfile(dirs.root, 'Analyses/fMRI_Scripts/Display_Results'));

dirs.save           = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots');
job                 = [];
job.lockSettings    = 'outcomelocked';
job.type            = 'standard';
job.iCope           = 6; 

% Alpha:
job.GLMID       = '3C';
job.cLim        = 35;
job.iView = 'sagittal'; job.iSlice = 0;
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5E_GLM%s_%s%d.png', job.GLMID, job.iView, job.iSlice)));
close gcf

% Save source data file:
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, 'pos')));
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5E_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, 'pos')));
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5E_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Theta:
job.GLMID       = '3A';
job.cLim        = 4;
job.iView = 'sagittal'; job.iSlice = 0;
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5F_GLM%s_%s%d.png', job.GLMID, job.iView, job.iSlice)));
close gcf

% Save source data file:
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, 'pos')));
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5F_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, 'pos')));
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5F_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% Beta:
job.GLMID       = '3B';
job.cLim        = 250;
job.iView = 'sagittal'; job.iSlice = 0;
EEGfMRIPav_sliceDisplay_cope(dirs, job);
saveas(gcf, fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5G_GLM%s_%s%d.png', job.GLMID, job.iView, job.iSlice)));
close gcf

% Save source data file:
sliceIdx        = floor(size(data, 1)/2) - job.iSlice/2; % convert into MNI coordinates
dirs.data       = fullfile(dirs.root, sprintf('Log/fMRI/GLM%s_FEAT_Combined_all', job.GLMID));
data            = niftiread(fullfile(dirs.data, sprintf('cope1_cope%d_%s.nii', job.iCope, 'pos')));
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5G_cope_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);
data            = niftiread(fullfile(dirs.data, sprintf('zstat1_cope%d_%s.nii', job.iCope, 'pos')));
saveMat         = squeeze(data(sliceIdx, :, :)); % imagesc(imrotate(saveMat, 90)); pbaspect([1 1 1]);
fullFileName    = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig5G_zstat_%s%d.csv', job.iView, job.iSlice));
csvwrite(fullFileName, saveMat);

% END OF FILE.