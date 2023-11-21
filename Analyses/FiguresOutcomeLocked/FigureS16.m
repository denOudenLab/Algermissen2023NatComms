% FigureS16.m

% Plots for Supplementary Figure S16.
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
%% Figure S16A: Line plot with highlighted patches:

% ----------------------------------------------------------------------- %
% Load saved regression coefficients for Updatebias:

EEGdomain           = 'time';

% 1) ROIs:
ROIs2use            = {''}; % empty
if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end

% 2) Behavioral regressors:
behav2use           = {'Updatebias', 'fbrel'}; % Updates + outcome sign

% 3) Perform only on selected trials:
selTrials           = 'all';

[job, dirs, betas]  = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
% Create figure:

job.invalidSubs = [11 12 23 30]; 

thresh      = 2.0;
nP          = 1000;
ylim        = 7;
lineWidth   = 4;
fontSize    = 32;

% Select channels:
selChans    = {'Fz', 'FCz', 'Cz'};

% Select ROI:
iROI        = 1;
        
% Set seed:
rng(20190822) % set random number generator for constant p-values

% Prepare data:
[sortBetas,~]   = taft_postprocess_time_selectData(job, betas, iROI);
sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:

% Start Plot:
close gcf

figure('Position', [100 100 1200 800], 'Color', 'white'); 

% Perform permutation test:
cfg             = [];
cfg.channel     = selChans;
cfg.avgoverchan = 'yes';
cfg.latency     = [0 0.7]; % time where ERPs occur
cfg.avgovertime = 'no';
fprintf('Regressor %s: Average over channels %s, reshape into 3D\n', ...
    char(job.regNames(iROI)), strjoin(selChans, '/'));

c               = [];
chanBetas       = cell(length(sortBetas), 1);
for iSub = 1:length(sortBetas) % iSub = 1;
    chanBetas{iSub}     = ft_selectdata(cfg, sortBetas{iSub}); % average over selected channels
    c(iSub, :, :)       = chanBetas{iSub}.avg; % bring into 4-D, with subject as first dimension
end
[corrp, tg]  = clustertf(c, thresh, nP);

% Extract timing and t-value and p-value info:
startTime   = cfg.latency(1);
endTime     = cfg.latency(end);
timeVec     = sortBetas{1}.time(find(sortBetas{1}.time==startTime):find(sortBetas{1}.time==endTime));
Tvec        = squeeze(tg);
pvec        = squeeze(corrp);

% Plot:
plot(chanBetas{1}.time, Tvec, 'k-', 'linewidth', lineWidth); hold on

% Add negative patches:
if any(pvec < .05 & Tvec < -2)
    sigTimeIdx  = (Tvec < -2 & pvec < 0.05); % negative
    sigTime     = find(sigTimeIdx==1);
    fprintf('Significant negative cluster from %.03f - %.03f sec.\n', timeVec(sigTime(1)), timeVec(sigTime(end)));
    patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))], [repmat(-2, 1, sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')], 'b', 'EdgeColor', 'none'); hold on
end

% Add positive patches:
if any(pvec < .05 & Tvec > 2)
    sigTimeIdx  = (Tvec > 2) & (pvec < 0.05); % positive
    sigTime     = find(sigTimeIdx==1);
    fprintf('Significant positive cluster from %.03f - %.03f sec.\n', timeVec(sigTime(1)), timeVec(sigTime(end)));
    patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))], [repmat(2, 1, sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')], 'b', 'EdgeColor', 'none'); hold on
end

% Horizontal lines:
nTime   = length(chanBetas{iSub}.time);
plot(chanBetas{iSub}.time, repelem(2, nTime), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 2);
plot(chanBetas{iSub}.time, repelem(-2, nTime), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 2);

% Final settings:
set(gca, 'xlim', cfg.latency, ....
    'xtick', [-1 -0.5 -0.25 0 0.25 0.5 0.7 1], ...
    'xtickLabel', {'-1', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
    'ylim', [-ylim ylim], 'ytick',-10:2:10, ...
    'fontsize', fontSize, 'LineWidth', lineWidth); % hard-coded so far
xlabel('Time (ms)', 'fontweight', 'bold', 'fontsize', fontSize); 
ylabel('T-values', 'fontweight', 'bold', 'fontsize', fontSize); 
box off

% Save:
figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS16A'));
saveas(gcf, [figName '.png'])
pause(1)
close gcf

% Save source data file:
fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS16A.csv');
csvwrite(fullFileName, Tvec);

fprintf('Done :-)\n')

% ----------------------------------------------------------------------- %
%% Figure S16A: Topoplots:

addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64));

% General settings:
zlim            = 4;
nCols           = 4;
lineWidth       = 3;

job.invalidSubs = [11 12 23 30]; % for Updatestd and Updatedif

% Timing settings:
startTime       = 0.0;
endTime         = 0.70;
steps           = 0.10;

% Create vectors:
startTimeVec    = startTime:steps:endTime;
endTimeVec      = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);
nPlot           = length(endTimeVec);
nRows           = ceil(nPlot/nCols); % number rows
saveMat         = nan(64, nPlot);

% Select ROI:
iROI            = 1;
        
% Prepare data:
[~, Tvalues] = taft_postprocess_time_selectData(job, betas, iROI);

% Start Plot:
figure('Position', [100 100 1200 800], 'Color', 'white'); 

for iPlot = 1:nPlot
    subplot(nRows, nCols, iPlot)

    cfg             = []; 
    cfg.figure      = gcf; 
    cfg.zlim        = [-1*zlim 1*zlim]; 
    cfg.marker      = 'on';
    cfg.style       = 'straight';
    cfg.layout      = 'easycapM11.mat'; 
    cfg.comment     = 'no';
    cfg.xlim        = [startTimeVec(iPlot) endTimeVec(iPlot)];
    cfg.colorbar    = 'no'; % want 'no', i.e. do it yourself % --> add externally
    cfg.style       = 'fill';
    ft_topoplotER(cfg, Tvalues);
    title(sprintf('%.1f-%.1f s', startTimeVec(iPlot), endTimeVec(iPlot)), 'fontsize', 28); % 1 digit: 28

    % Set linewidth:
    p = gca;
    for i = 7:11
        t = p.Children(i); 
        t.LineWidth = lineWidth;
    end
    
    % Save:
    timeIdx             = dsearchn(Tvalues.time', [startTimeVec(iPlot) endTimeVec(iPlot)]'); 
    timeIdx             = timeIdx(1):timeIdx(end);
    saveMat(:, iPlot)   = squeeze(mean(Tvalues.avg(:, timeIdx), 2));

end

% Save:
figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS16A_topo'));
saveas(gcf, [figName '.png'])

% Close:
pause(2)
close gcf
fprintf('Done :-)\n')

% Save source data file:
saveMat         = horzcat(Tvalues.label, num2cell(saveMat));
fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS16A_topo.csv'));
writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);

% ----------------------------------------------------------------------- %
%% Figures S16B-C: Line plots with highlighted patches:

% ----------------------------------------------------------------------- %
% Load saved regression coefficients for Updatebias:

EEGdomain           = 'time';

% 1) ROIs:
ROIs2use            = {''}; % empty
if ~exist('ROIs2use', 'var'); ROIs2use = {''}; end

% 2) Behavioral regressors:
behav2use           = {'Updatestd', 'Updatedif', 'fbrel'}; % Updates + outcome sign

% 3) Perform only on selected trials:
selTrials           = 'all';

[job, dirs, betas]  = taft_postprocess_load_job(EEGdomain, ROIs2use, behav2use, selTrials);

% ----------------------------------------------------------------------- %
% Create figure:

job.invalidSubs = [11 12 23 30]; 

thresh      = 2.0;
nP          = 1000;
ylim        = 7;
lineWidth   = 4;
fontSize    = 32;

% Select channels:
selChans        = {'Fz', 'FCz', 'Cz'};

% Select ROI:
for iROI = 1:2

    % Set seed:
    rng(20190822) % set random number generator for constant p-values

    % Prepare data:
    [sortBetas,~]   = taft_postprocess_time_selectData(job, betas, iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:

    % Start Plot:
    close gcf

    figure('Position', [100 100 1200 800], 'Color', 'white'); 

    % Perform permutation test:
    cfg             = [];
    cfg.channel     = selChans;
    cfg.avgoverchan = 'yes';
    cfg.latency     = [0 0.7]; % time where ERPs occur
    cfg.avgovertime = 'no';
    fprintf('Regressor %s: Average over channels %s, reshape into 3D\n', ...
        char(job.regNames(iROI)), strjoin(selChans, '/'));

    c               = [];
    chanBetas       = cell(length(sortBetas), 1);
    for iSub = 1:length(sortBetas) % iSub = 1;
        chanBetas{iSub}     = ft_selectdata(cfg, sortBetas{iSub}); % average over selected channels
        c(iSub, :, :)       = chanBetas{iSub}.avg; % bring into 4-D, with subject as first dimension
    end
    [corrp, tg]  = clustertf(c, thresh, nP);

    % Extract timing and t-value and p-value info:
    startTime   = cfg.latency(1);
    endTime     = cfg.latency(end);
    timeVec     = sortBetas{1}.time(find(sortBetas{1}.time==startTime):find(sortBetas{1}.time==endTime));
    Tvec        = squeeze(tg);
    pvec        = squeeze(corrp);

    % Plot:
    plot(chanBetas{1}.time, Tvec, 'k-', 'linewidth', lineWidth); hold on

    % Add negative patches:
    if any(pvec < .05 & Tvec < -2)
        sigTimeIdx  = (Tvec < -2 & pvec < 0.05); % negative
        sigTime     = find(sigTimeIdx==1);
        fprintf('Significant negative cluster from %.03f - %.03f sec.\n', timeVec(sigTime(1)), timeVec(sigTime(end)));
        patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))], [repmat(-2, 1, sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')], 'b', 'EdgeColor', 'none'); hold on
    end

    % Add positive patches:
    if any(pvec < .05 & Tvec > 2)
        sigTimeIdx  = (Tvec > 2) & (pvec < 0.05); % positive
        sigTime     = find(sigTimeIdx==1);
        fprintf('Significant positive cluster from %.03f - %.03f sec.\n', timeVec(sigTime(1)), timeVec(sigTime(end)));
        patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))], [repmat(2, 1, sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')], 'b', 'EdgeColor', 'none'); hold on
    end

    % Horizontal lines:
    nTime   = length(chanBetas{iSub}.time);
    plot(chanBetas{iSub}.time, repelem(2, nTime), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 2);
    plot(chanBetas{iSub}.time, repelem(-2, nTime), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 2);

    % Final settings:
    set(gca, 'xlim', cfg.latency, ....
        'xtick', [-1 -0.5 -0.25 0 0.25 0.5 0.7 1], ...
        'xtickLabel', {'-1', '-500', '-250', 'FbOn', '250', '500', 'FbOff', '1'}, ...
        'ylim', [-ylim ylim], 'ytick',-10:2:10, ...
        'fontsize', fontSize, 'LineWidth', lineWidth); % hard-coded so far
    box off
    xlabel('Time (ms)', 'fontweight', 'bold', 'fontsize', fontSize); 
    ylabel('T-values', 'fontweight', 'bold', 'fontsize', fontSize); 

    % Save:
    if strcmp(job.regNames(iROI), 'Updatestd')
        figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS16B'));
        saveas(gcf, [figName '.png'])
        % Save source data file:
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS16B.csv');
        csvwrite(fullFileName, Tvec);
    end
    if strcmp(job.regNames(iROI), 'Updatedif')
        figName = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS16C'));
        saveas(gcf, [figName '.png'])
        % Save source data file:
        fullFileName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/FigS16C.csv');
        csvwrite(fullFileName, Tvec);
    end
    pause(3);
    close gcf
        
end % end iROI
fprintf('Done :-)\n')

% ----------------------------------------------------------------------- %
%% Figures S16B-C: Topoplots:

addpath(fullfile(dirs.root, 'Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', redblue(64));

% General settings:
zlim            = 4;
nCols           = 4;
lineWidth       = 3;

job.invalidSubs = [11 12 23 30]; % for Updatestd and Updatedif

% Timing settings:
startTime       = 0.0;
endTime         = 0.70;
steps           = 0.10;

% Create vectors:
startTimeVec    = startTime:steps:endTime;
endTimeVec      = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);
nRows           = ceil(length(startTimeVec)/nCols); % number rows

for iROI = 1:2
        
    % Initialize for saving:
    saveMat         = nan(64, nPlot);

    % Prepare data:
    [~, Tvalues] = taft_postprocess_time_selectData(job, betas, iROI);
    
    % Start Plot:
    figure('Position', [100 100 1200 800], 'Color', 'white'); 
    
    for iPlot = 1:length(endTimeVec)
        subplot(nRows, nCols, iPlot)
              
        cfg             = []; 
        cfg.figure      = gcf; 
        cfg.zlim        = [-1*zlim 1*zlim]; 
        cfg.marker      = 'on';
        cfg.style       = 'straight';
        cfg.layout      = 'easycapM11.mat'; 
        cfg.comment     = 'no';
        cfg.xlim        = [startTimeVec(iPlot) endTimeVec(iPlot)];
        cfg.colorbar    = 'yes'; % want 'no', i.e. do it yourself % --> add externally
        cfg.style       = 'fill';
        ft_topoplotER(cfg, Tvalues);
        title(sprintf('%.1f-%.1f s', startTimeVec(iPlot), endTimeVec(iPlot)), 'fontsize', 28); % 1 digit: 28
        
        % Set linewidth:
        p = gca;
        for i = 7:11
            t = p.Children(i); 
            t.LineWidth = lineWidth;
        end

        % Save:
        timeIdx             = dsearchn(Tvalues.time', [startTimeVec(iPlot) endTimeVec(iPlot)]'); 
        timeIdx             = timeIdx(1):timeIdx(end);
        saveMat(:, iPlot)   = squeeze(mean(Tvalues.avg(:, timeIdx), 2));

    end
    
    % Save:
    if strcmp(job.regNames(iROI), 'Updatestd')
        figName         = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS16B_topo'));
        fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS16B_topo.csv'));
    end
    if strcmp(job.regNames(iROI), 'Updatedif')
        figName         = fullfile(dirs.root, sprintf('/Log/OutcomeLockedPaperPlots/FigS16C_topo'));
        fullFileName 	= fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/FigS16C_topo.csv'));
    end
        saveas(gcf, [figName '.png'])
    
    % Save source data file:
    saveMat         = horzcat(Tvalues.label, num2cell(saveMat));
    writetable(cell2table(saveMat), fullFileName, 'WriteVariableNames', 0);

    % Close:
    pause(2)
    close gcf
end
fprintf('Done :-)\n')

% END OF FILE.