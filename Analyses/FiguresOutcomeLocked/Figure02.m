% Figure02.m

% Plots for Figure 2A-H in manuscript.
% The same approach is used for Supplementary Figures S01, S06, S07. 
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2023.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresOutcomeLocked

clear all; close all; clc

% Set root directory:
cd(fileparts(which('figures_set_rootDir.m')));
dirs.root           = figures_set_rootDir();

% Set other directories:
dirs.log            = fullfile(dirs.root, 'Log');

% Add additional paths for scripts:
addpath(fullfile(dirs.root, 'Analyses/Behavior_Scripts/Matlab_Plots'));

% ----------------------------------------------------------------------- %
%% Alternative A: Load data:

% Execute for figures 2A, 2B, 2C.
% Skip for figures 2E, 2F, 2G.

% Retrieve recoded and pre-processed behavioral data:
out             = EEGfMRIPav_aggr_data(); % execute to get out
datType         = 'data'; % this is empirical data (not simulated data)

% Extract relevant data objects:
pGoCond         = out.pGoCond;
pCorrectCond    = out.pCorrectCond;
pStay           = out.stayGoOutVal;

% Data dimensions:
nSub            = size(pGoCond, 1);
nCond           = size(pGoCond, 2);
nRep            = size(pGoCond, 3);

fprintf('Finished loading behavioral data, already pre-processed\n');

% --> Move on to 'Selecting subjects' and then Figures 2A, ...

% ----------------------------------------------------------------------- %
%% Alternative B: Load simulations: 	

% Execute for figures 2E, 2F, 2G.
% Skip for figures 2A, 2B, 2C.

% Set model:
iMod        = 5; % 1-9 % model ID

% Model 1: Q-learner
% Model 2: Q-learner + Go bias.
% Model 3: Q-learner + Go bias + Pavlovian response bias.
% Model 4: Q-learner + Go-bias + Pavlovian learning bias.
% Model 5: Q-learner + Go-bias + Pavlovian response bias + Pavlovian learning bias.
% Model 6: action priming model.
% Model 7: M5 + single perseverance parameter.
% Model 8: M5 + dual (cue valence-dependent) perseverance parameter.
% Model 9: M5 + neutral outcome reinterpretation.

% Set fitting settings, simulation settings:
parType     = 'lap'; % lap hbi % fitting: LaPlace approximation or hierarchical Bayesian inference
simType     = 'osap'; % modSim osap % simulations: model simulations or one-step-ahead predictions

% Downstream settings:
datType     = sprintf('%s_%s_M%02d', simType, parType, iMod);

if strcmp(simType, 'osap')
    nIter   = 1;
else
    nIter   = 100;
end

% Directories to load from:
dirs.results    = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');
dirs.sims       = fullfile(dirs.results, 'Simulations');

% Input file name:
inputfile       = fullfile(dirs.sims, sprintf('%s_mod%02d_%s_iter%04d', simType, iMod, parType, nIter));
inputfile       = [inputfile '.mat'];

% Load simulations:
fprintf('Load %s for model %02d based on %s\n', simType, iMod, parType);
tmp             = load(inputfile);
simulations     = tmp.sim;

% Settings for assigning stimuli to conditions:
nSub            = 36;
nStim           = 16;
nCond           = 4;
nRep            = 40;
condMatrix      = [1 3 5 7; 2 4 6 8; 9 11 13 15; 10 12 14 16]; % assign stimuli to conditions

% ------------------------------------- %
% a) If simulations: Aggregate into pGoCond:

fprintf('Prepare pGoCond based on simulations\n');

% Average over iterations:
pGo         = squeeze(nanmean(simulations.pGo, 2));

% Sort into stimuli and then into conditions:
pGoStim     = nan(nSub, nStim, nRep);
pGoCond     = nan(nSub, nCond, nRep);

% Loop over subjects, extract and aggregate simulated data:
for iSub = 1:nSub
    
    % Extract subject data:
    subj = simulations.subj{iSub};
    
    % Sort into stimuli:
    for iStim = 1:nStim
        idx                     = subj.stimuli == iStim; % indices of trials where this stimuli appears
        pGoStim(iSub, iStim, :) = pGo(iSub, idx);
    end
    
    % Sort into conditions:
    for iCond = 1:nCond % iCond = 1;
        idx                     = condMatrix(:, iCond); % indices of stimuli in this condition
        pGoCond(iSub, iCond, :) = mean(pGoStim(iSub, idx, :), 2); % average over selected stimuli
    end
end

% ------------------------------------- %
% b) If simulations: Aggregate into pCorrectCond:

fprintf('Prepare pCorrectCond based on simulations\n');

% Sort into stimuli:
pCorrectStim = nan(nSub, nIter, nStim, nRep);

for iSub = 1:nSub % iSub = 1;
    
    % Extract subject data:
    subj = simulations.subj{iSub};
    
    % Sort into stimuli:
    for iIter = 1:nIter
        
        if strcmp(simType, 'osap')
            respIter    = squeeze(subj.actions); % single iteration
        elseif strcmp(simType, 'modSim')
            respIter    = squeeze(simulations.action(iSub, iIter, :)); % actions in this iteration
        else
            error('Invalid datType')
        end
        
        % Loop over stimuli, aggregate per stimulus:
        for iStim = 1:nStim
            idx         = subj.stimuli == iStim; % indices of trials where this stimuli appears
            respStim    = respIter(idx)'; % actions to these stimuli
            tmp         = squeeze(simulations.p(iSub, iIter, idx, :)); % extract nRep x nStim matrix
            linIdx      = sub2ind(size(tmp), 1:nRep, respStim); % linear indexing: action per rep
            pCorrectStim(iSub, iIter, iStim, :)   = tmp(linIdx); % apply linear indexing
        end
    end
end

% Average over iterations:
pCorrectStim = squeeze(nanmean(pCorrectStim, 2)); % average over iterations
    
% Sort into conditions:
pCorrectCond = nan(nSub, nCond, nRep);
for iSub = 1:nSub % iSub = 1;
    for iCond = 1:nCond % iCond = 1;
        idx                         = condMatrix(:, iCond); % indices of stimuli in this condition
        pCorrectCond(iSub, iCond, :) = mean(pCorrectStim(iSub, idx, :), 2); % average over selected stimuli
    end
end

% ------------------------------------- %
% c) If simulations: Prepare pStay:

fprintf('Prepare p(stay) based on simulations\n');

stayIter = nan(iSub, nIter, iCond); % initialize

for iSub = 1:nSub % iSub = 1;
    
    subj = simulations.subj{iSub}; % subject data here
    
    valVec      = ismember(subj.stimuli, [1 2 5 6 9 10 13 14]); % cue valence: win cue (1) or not (0)

    % ------------------------------------------------------------------- %
    % Count stimulus repetitions (same stimulus order for each iteration):
    
    % Initialize for counting stimulus repetition:
    nStim       = length(unique(subj.stimuli));
    nTrial      = length(subj.stimuli);
    stimCount   = zeros(nStim, 1);
    stimRep     = nan(nTrial, 1);
    
    for iTrial = 1:nTrial
        thisStim            = subj.stimuli(iTrial); % retrieve stimulus ID
        stimCount(thisStim) = stimCount(thisStim) + 1; % increment count for this stimulus ID
        stimRep(iTrial)     = stimCount(thisStim); % store stimulus count
    end
    
    % ------------------------------------------------------------------- % 
    % Loop over iterations:
    for iIter = 1:nIter % iIter = 1;

        if strcmp(simType,'osap') 
            actVec      = subj.actions;
            respVec     = ismember(subj.actions, [1 2]); % Go response (1) or not (0)
            outVec      = subj.outcome;
            
        elseif strcmp(simType,'modSim') 
            actVec      = squeeze(simulations.action(iSub, iIter, :));
            respVec     = squeeze(ismember(simulations.action(iSub, iIter, :), [1 2]));
            outVec      = squeeze(simulations.outcome(iSub, iIter, :));
            
        else
            error('Unknown simType')
        end
        

    % ------------------------------------------------------------------- % 
    % Compute p(Stay):
        
        stay = nan(nTrial, 1); %initialize
        for iTrial = 1:nTrial % iTrial = 1;
            
            thisStim    = subj.stimuli(iTrial); % retrieve stimulus identifier
            thisRep     = stimRep(iTrial); % retrieve cumulative number of repetitions of this stimulus
            nextTrial   = find(subj.stimuli == thisStim & stimRep == (thisRep + 1)); % next trial with same stimulus
            
            if thisRep < nRep % if there is still a "next trial" left
                thisAct         = actVec(iTrial); % action on next trial
                stay(iTrial)    = simulations.p(iSub, iIter, nextTrial, thisAct); % probably for iTrial action on next trial
            end
            
        end

        % --------------------------------------------------------------- % 
        % Recompute interpretation of outcome obtained:
        outVecRel               = 1 - outVec; % Avoid cues: 0-->1 becomes good, -1-->2 becomes bad; mind next line for completion
        outVecRel(valVec == 1)  = outVecRel(valVec ==1 ) + 1; % Win cues: one up (1-->0-->1 becomes 1, 0-->1-->2 becomes 2)
        outVecAll               = outVecRel + 2 * (valVec == 0); % 1-4: reward, no-reward, no-punishment, punishment 

        % --------------------------------------------------------------- % 
        % Loop over actions and outcomes, find trials, average:
        
        iCond = 0; % initialize condition count

        for iResp = [1 0] % response made
            for iOut = 1:4 % outcome obtained

                iCond = iCond + 1; % increment condition count

                % Find trials with this action and this outcome:
                idx                 = respVec == iResp & outVecAll == iOut;

                % Mean stay in this condition:
                stayIter(iSub, iIter, iCond)  = nanmean(stay(idx));

            end % end iResp
        end % end iOut
    end % end iIter 
end % end iSub

pStay = squeeze(nanmean(stayIter, 2)); % average across iterations
fprintf('Finished preparing data from simulations :-)\n');

% ----------------------------------------------------------------------- %
%% Select subjects:

% Run this section to generate 'validSubs' to be used in following sections

fprintf('Select subjects\n');

invalidSubs     = []; % keep all subjects
fprintf('Exclude subjects %s\n', num2str(invalidSubs));
validSubs       = setdiff(1:nSub, invalidSubs);
nSubValid       = length(validSubs);

% ======================================================================= %
% ======================================================================= %
% ======================================================================= %
% ======================================================================= %
%% Figures 2A & 2E: LEARNING CURVES:

fprintf('Plot learning curves per condition (pGoCond) based on %s data\n',datType);

% Select subjects, grand means:
nCond       = size(pGoCond, 2);
subMean     = squeeze(nanmean(pGoCond(validSubs, :, :), 2)); % average across conditions per subject
grandMean   = squeeze(nanmean(subMean, 1)); % average across subjects --> grand mean, single value

% Average per condition:
condMean    = nan(nCond, nRep);
condSE      = nan(nCond, nRep);

% Loop over conditions, compute mean per condition across subjects, compute
% standard error across subjects using Cousineau-Morey correction:
for iCond = 1:nCond
    condMean(iCond, :) = squeeze(nanmean(pGoCond(validSubs, iCond, :), 1));
    condSE(iCond, :)   = nCond / (nCond-1) * nanstd(squeeze(pGoCond(validSubs, iCond, :)) - ...
        subMean + repmat(grandMean, nSubValid, 1)) ./ sqrt(nSubValid); % Cousineau-Morey correction
end

% Fixed plotting settings:
% https://www.visualisingdata.com/2019/08/five-ways-to-design-for-red-green-colour-blindness/
colMat      = [0 113 116; 201 61 33; 0 113 116; 201 61 33] ./ 255; % Var 1: dark green, light green, orange, red
fontSize    = 38;
lineWidth   = 4;
xWidth      = 800;

% ------------------------------------------------------ %
% Start figure:

p = cell(nCond, 1);

figure('Position', [100 100 900 800], 'Color', 'white'); 

% ------------------------------------------------------ %
% Plot line with error shade:
for iCond = 1:nCond
    p{iCond} = boundedline(1:nRep, condMean(iCond, :), condSE(iCond, :), ...
        'cmap', colMat(iCond, :), 'alpha');
    set(p{iCond}, 'Linewidth', lineWidth) % adjust line width
    if iCond > nCond/2  % adjust line style for NoGo cues
        set(p{iCond}, 'Linestyle', '--');
    end
    
end

% Add plot features:
set(gca,'xlim', [0 nRep], 'ylim', [0 1], 'Linewidth', 3);
set(gca,'xtick', 0:10:nRep, 'ytick', 0:.2:1, 'FontSize', fontSize, 'Linewidth', 4);
xlabel('Trial', 'FontSize', fontSize, 'FontName', 'Arial');
ylabel('p(Go)', 'FontSize', fontSize, 'FontName', 'Arial');
box off

% Determine figure name:
if strcmp(datType, 'data')
    figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig2A');
else
    figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig2E');
end

% ------------------------------------------------------ %
% Save:
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data files:
if strcmp(datType, 'data')
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2A_G2W.csv'));
    csvwrite(fullFileName, squeeze(pGoCond(validSubs, 1, :)));
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2A_G2A.csv'));
    csvwrite(fullFileName, squeeze(pGoCond(validSubs, 2, :)));
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2A_NG2W.csv'));
    csvwrite(fullFileName, squeeze(pGoCond(validSubs, 3, :)));
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2A_NG2A.csv'));
    csvwrite(fullFileName, squeeze(pGoCond(validSubs, 4, :)));
else
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2E_G2W.csv'));
    csvwrite(fullFileName, squeeze(pGoCond(validSubs, 1, :)));
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2E_G2A.csv'));
    csvwrite(fullFileName, squeeze(pGoCond(validSubs, 2, :)));
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2E_NG2W.csv'));
    csvwrite(fullFileName, squeeze(pGoCond(validSubs, 3, :)));
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2E_NG2A.csv'));
    csvwrite(fullFileName, squeeze(pGoCond(validSubs, 4, :)));
end
% ======================================================================= %
% ======================================================================= %
% ======================================================================= %
%% Figures 2B & 2F: p(Go) BAR PLOTS

fprintf('Plot mean pGo and pCorrect per condition (pGoCond and pCorrectCond) based on %s data\n', datType);

% Select subjects, average across trials, grand means:
nCond               = size(pGoCond, 2);
pGoCondMean         = squeeze(nanmean(pGoCond(validSubs, :, :), 3)); % average across trials
pCorrectCondMean    = squeeze(nanmean(pCorrectCond(validSubs, :, :), 3)); % average across trials
subMean             = squeeze(nanmean(pGoCondMean, 2)); % average across conditions
grandMean           = squeeze(nanmean(subMean, 1)); % average across subjects

% Loop over conditions, compute mean per condition across subjects, compute
% standard error across subjects using Cousineau-Morey correction:
condMean            = nan(nCond, 1);
condMeanCorrect     = nan(nCond, 1);
condSE              = nan(nCond, 1);
for iCond = 1:nCond
    condMean(iCond)         = squeeze(nanmean(pGoCondMean(:, iCond), 1));
    condSE(iCond)           = nCond/(nCond-1) * nanstd(squeeze(pGoCondMean(:, iCond)) - ...
        subMean + repmat(grandMean, nSubValid, 1)) ./ sqrt(nSubValid);
    condMeanCorrect(iCond)  = squeeze(nanmean(pCorrectCondMean(:, iCond), 1));
end

% General plot settings:
colMat    = [0 113 116; 201 61 33; 0 113 116; 201 61 33] ./ 255; % dark green, light green, orange, red
colMatCor = [87 196 173; 240 174 102] ./ 255; % dark green, light green, orange, red
xLoc      = [1 2 3.5 4.5];
lineWidth = 4; capSize = 12; fontSize = 36;

% ------------------------------------------------------ %
% Start figure:

close all
figure('Position', [100 100 900 800], 'Color', 'white'); hold on

% ------------------------------------------------------ %
% a) Plot bars with errorbars:
for iCond = 1:nCond
    bar(xLoc(iCond), condMean(iCond), .75, 'FaceColor', colMat(iCond, :)); % bar plot
    errorbar(xLoc(iCond), condMean(iCond), condSE(iCond), ...
        'k', 'linestyle', 'none', 'linewidth', lineWidth, 'Capsize', capSize); % error bars
end

% ------------------------------------------------------ %
% b) Plot bars with correct Gos:
for iCond = 1:2
    bar(xLoc(iCond), condMeanCorrect(iCond), .75, 'FaceColor', colMatCor(iCond, :));
end

% ------------------------------------------------------ %
% c) Points per subject:
for iCond = 1:nCond
    s = scatter(repmat(xLoc(iCond), 1, nSubValid), pGoCondMean(:, iCond)', ...
        [], 'k', 'jitter', 'on', 'jitterAmount', 0.15); hold on
    set(s, 'MarkerEdgeColor', [0.4 0.4 0.4], 'linewidth', 3);
end

% ------------------------------------------------------ %
% Add plot features:
set(gca,'xlim', [.5 5.5], 'ylim', [0 1], 'xtick', 0:10:nRep,...
    'xtick', [1.5 4], 'xticklabel', {'Go','NoGo'}, 'ytick', 0:.2:1,...
    'FontSize', fontSize, 'FontName', 'Arial', 'FontWeight', 'normal', 'Linewidth', 4);
ylabel('p(Go)', 'FontSize', fontSize, 'FontName', 'Arial');
xlabel('Required Action', 'FontSize',fontSize,'FontName', 'Arial');
box off

% ------------------------------------------------------ %
% Save:

% Determine figure name:
if strcmp(datType, 'data')
    figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig2B');
else
    figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig2F');
end
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data files:
if strcmp(datType, 'data')
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2B_pGo.csv'));
    csvwrite(fullFileName, pGoCondMean);
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2B_pCorrect.csv'));
    csvwrite(fullFileName, pCorrectCondMean);
else
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2F_pGo.csv'));
    csvwrite(fullFileName, pGoCondMean);
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2F_pCorrect.csv'));
    csvwrite(fullFileName, pCorrectCondMean);
end

% ======================================================================= %
% ======================================================================= %
% ======================================================================= %
%% Figure 2C: Win-stay/ Lose-stay:

fprintf('Plot win-stay/lose-shift (pStay) based on %s data\n',datType);

% Grouping of conditions next to each other (see grouping below):
datatype = 'SalActOutVal'; % hierarchy: salience, action, outcome valence --> PAPER
% datatype = 'ActCueValOutVal'; % hierarchy: action, cue valence, outcome valence
% datatype = 'CueValOutValAct'; % hierarchy: cue valence, outcome valence, action
% datatype = 'SalOutValAct'; % hierarchy: salience, outcome valence, action

% General plot settings:
nCond     = size(pStay, 2);
condNames = {'Go Reward', 'Go No Reward', 'Go No Punishment', 'Go Punishment', 'NoGo Reward', 'NoGo No Reward', 'NoGo No Punishment', 'NoGo Punishment'};
colAll     = [0 113 116; 87 196 173; 240 174 102; 201 61 33] ./ 255; % Var 1: dark green, light green, orange, red

% Set x-axis locations:
xLoc      = nan(nCond, 1); % position of bars on x-axis
for iCond = 1:nCond
    if iCond <= nCond/2
        xLoc(iCond) = iCond; % first half
    else
        xLoc(iCond) = iCond + 0.5; % second half
    end
end

% ------------------------------------------------------ %
% Order all 8 conditions into the correct hierarchy:
if strcmp(datatype, 'SalActOutVal') % Option 1: Salience - Action - Outcome Valence   
    data        = pStay;
    condOrder   = [1 4 5 8 3 2 7 6]; % GoRew GoPun NoGoRew NoGoPun GoNoPun GoNoRew NoGoNoPun NoGoNoRew
    colMat     = colAll([1:nCond/2 1:nCond/2], :); % just as many colors as conditions
    % x-labels for action only:
    xTick       = [(xLoc(1)+xLoc(2))/2 (xLoc(3)+xLoc(4))/2 (xLoc(5)+xLoc(6))/2 (xLoc(7)+xLoc(8))/2];
    xTickLabel  = {'Go', 'NoGo', 'Go', 'NoGo'};
    xLabel      = 'Performed action';

elseif strcmp(datatype, 'ActCueValOutVal') % Option 2: Action - Cue Valence - Outcome Valence
    data        = pStay;
    condOrder   = [1 2 3 4 5 6 7 8]; % GoRew, GoNoRew, GoNoPun, GoPun, NoGoRew, NoGoNoRew, NoGoNoPun, NoGoPun
    colMat      = colAll([1:nCond/2 1:nCond/2], :); % just as many colors as conditions
    % x-labels for action only:
    xTick       = [(xLoc(2)+xLoc(3))/2 (xLoc(6)+xLoc(7))/2];
    xTickLabel  = {'Go', 'NoGo'};
    xLabel      = 'Performed action';
    
elseif strcmp(datatype,'CueValOutValAct') % Option 3: Cue Valence - Outcome Valence - Action
    data        = pStay;
    condOrder   = [1 5 2 6 3 7 4 8]; % GoRew, NoGoRew, GoNoRew, NoGoNoRew, GoNoPun, NoGoNoPun, GoPun, NoGoNoPun
    colMat      = colAll([1:nCond/2 1:nCond/2], :); % just as many colors as conditions
    % x-labels for valence only:
    xTick       = [(xLoc(2)+xLoc(3))/2 (xLoc(6)+xLoc(7))/2];
    xTickLabel  = {'Win ', 'Avoid'};
    xLabel      = 'Cue valence';
    
elseif strcmp(datatype,'SalOutValAct') % Option 4: Salience - Outcome Valence - Action
    data        = pStay;
    condOrder   = [1 5 4 8 2 6 3 7]; % GoRew, NoGoRew, GoPun, NoGoPun, GoNoRew, NoGoNoRew, GoNoPun, NoGoNoPun
    colMat      = colAll([1:nCond/2 1:nCond/2], :); % just as many colors as conditions
    % x-labels for salience only:
    xTick       = [(xLoc(2)+xLoc(3))/2 (xLoc(6)+xLoc(7))/2];
    xTickLabel  = {'Salient', 'Neutral'};
    xLabel      = 'Salience';
    
else
    error('Unknown data type')
end

% Loop over conditions, compute mean per condition across subjects, compute
% standard error across subjects using Cousineau-Morey correction:
nCond               = size(data, 2);
subMean             = squeeze(nanmean(data(validSubs, :), 2)); % average across conditions
grandMean           = squeeze(nanmean(subMean, 1)); % average across subjects
condMean            = nan(nCond, 1);
condSE              = nan(nCond, 1);
for iCond = 1:nCond
    condMean(iCond)         = squeeze(nanmean(data(validSubs, iCond), 1)); % average over subjects
    fprintf('Condition %s: M = %.02f\n', condNames{iCond}, condMean(iCond));
    condSE(iCond)           = nCond/(nCond-1)*nanstd(squeeze(data(validSubs,iCond)) - ...
        subMean + repmat(grandMean, nSubValid, 1)) ./ sqrt(nSubValid);
end

% Fixed plot settings:
xMax      = xLoc(nCond) + 0.5;
lineWidth = 4; capSize = 12; fontSize = 36;

% ------------------------------------------------------ %
% Start figure:

close all
p = cell(nCond, 1);
figure('Position', [100 100 900 800], 'Color', 'white'); hold on

% ------------------------------------------------------ %
% a) Plot bars with errorbars:
for iCond = 1:nCond
    condIdx = condOrder(iCond);
    p{iCond} = bar(xLoc(iCond), condMean(condIdx), .75, 'FaceColor', colMat(condIdx, :)); % bar plot
    errorbar(xLoc(iCond), condMean(condIdx), condSE(condIdx), ...
        'k', 'linestyle', 'none', 'linewidth', lineWidth, 'Capsize', capSize); % error bars
end

% ------------------------------------------------------ %
% b) Points:
for iCond = 1:nCond
    condIdx = condOrder(iCond);
    s = scatter(repmat(xLoc(iCond), 1, nSubValid), data(validSubs, condIdx)', [], ...
        'k', 'jitter' ,'on', 'jitterAmount', 0.15); hold on % was 0.05
    set(s, 'MarkerEdgeColor', [0.4 0.4 0.4], 'linewidth', 3); % was 1 
end
plot([0 xMax], ones(2, 1)*1/3, 'k--', 'linewidth', lineWidth); % line at chance

% ------------------------------------------------------ %
% Add plot features:
set(gca, 'xlim', [.5 xMax], 'ylim', [0 1],...
    'xtick', xTick, 'xticklabel', xTickLabel, 'ytick', 0:.2:1,...
    'FontSize', fontSize, 'FontName', 'Arial', 'FontWeight', 'normal', 'Linewidth', 4)
xlabel(xLabel, 'FontSize', fontSize, 'FontName', 'Arial');
ylabel('p(Stay)', 'FontSize', fontSize, 'FontName', 'Arial');

% ------------------------------------------------------ %
% Save:

% Determine figure name:
if strcmp(datType, 'data')
    figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig2C');
else
    figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig2G');
end
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data files:
if strcmp(datType, 'data')
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2C.csv'));
    csvwrite(fullFileName, data(validSubs, condOrder));
else
    fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2G.csv'));
    csvwrite(fullFileName, data(validSubs, condOrder));
end

% ======================================================================= %
% ======================================================================= %
% ======================================================================= %
%% Figure 2D: LOG MODEL EVIDENCE

% Connect data points of the same subject or not:
addLines    = false;

% Relevant settings:
modVec      = 1:5; 
modNames    = modVec;

% For x-axis labels of models and file name:
nMod        = length(modVec);
modLabels   = strcat('M', string(modNames));

% Directories:
dirs.results = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');
dirs.lap     = fullfile(dirs.results, 'LAP_Results');

% Construct names of input files:
fname_mod = struct([]);
for iMod = 1:nMod % keep absolute values for indexing
    modIdx = modVec(iMod); % retrieve model to be used
    fname_mod{iMod} = fullfile(dirs.lap,sprintf('lap_mod%02d.mat', modIdx));
end

% Load:
lME = nan(nSub, nMod); % initialize
for iMod = 1:nMod %  keep absolute values for indexing % nMod = 1;
    fname = load(fname_mod{iMod}); % load data
    cbm   = fname.cbm; % extract cbm object
    lME(:,iMod) = cbm.output.log_evidence;
end

% Loop over model, compute mean per model across subjects, compute
% standard error across subjects using Cousineau-Morey correction:
subMean             = squeeze(nanmean(lME(validSubs, :), 2)); % average across conditions
grandMean           = squeeze(nanmean(subMean, 1)); % average across subjects
condMean            = nan(nMod, 1);
condSE              = nan(nMod, 1);
for iMod = 1:nMod
    condMean(iMod)         = squeeze(nanmean(lME(:,iMod), 1)); % average over subjects
    fprintf('Model %d: mean log model evidence = %.02f\n', iMod, condMean(iMod));
    condSE(iMod)           = nMod / (nMod - 1) * nanstd(squeeze(lME(validSubs, iMod)) - ...
        subMean + repmat(grandMean, nSubValid, 1)) ./ sqrt(nSubValid);
end

% General plot settings:
colMat      = repmat([175 175 175] ./ 255, nMod, 1); % grey
xLoc        = 1:nMod;
xTickLabel  = modLabels;
lineWidth = 4; capSize = 12; fontSize = 36;

% ------------------------------------------------------ %
% Start plot:
figure('Position', [100 100 900 800], 'Color', 'white'); hold on

% ------------------------------------------------------ %
% a) Plot bars with errorbars:

for iCond = 1:nMod
    % Bars:
    bar(xLoc(iCond), condMean(iCond), .75, 'FaceColor', colMat(iCond, :)); % bar plot
    errorbar(xLoc(iCond), condMean(iCond), condSE(iCond), ...
        'k', 'linestyle', 'none', 'linewidth', lineWidth, 'Capsize', capSize); % error bars
    % Points:
    s = scatter(repmat(xLoc(iCond), 1, nSubValid), lME(validSubs, iCond)', ...
        [],'k', 'jitter', 'on', 'jitterAmount', 0.15); hold on % was 0.05
    set(s, 'MarkerEdgeColor', [0.4 0.4 0.4], 'linewidth', 3); % was 1 
end

% ------------------------------------------------------ %
% b) Lines:
if addLines
    for iSub = 1:nSubValid
        p = plot(xLoc, lME(iSub, :), 'k-', 'lineWidth', 1); hold on % was 0.05
        p.Color(4) = 0.20;
    end
end

% ------------------------------------------------------ %
% Y-axis limits:
leeway = 20;
yLim = [floor(min(lME(:))) - leeway, ceil(max(lME(:))) + leeway]; % yLim based on individual data points

% Further features:
set(gca,'xlim', [0 max(xLoc) + 0.5], 'ylim', yLim, ...
    'xtick',xLoc, 'xticklabel', xTickLabel, ...
    'FontSize', fontSize, 'FontName', 'Arial', 'FontWeight', 'normal', 'Linewidth', 4);
% Labels:
xlabel('Model','FontSize',fontSize,'FontName','Arial');
ylabel('log model evidence','FontSize',fontSize,'FontName','Arial');

% ------------------------------------------------------ %
% Save:
figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig2D');
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data files:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2D.csv'));
csvwrite(fullFileName, lME(validSubs, :));

% ======================================================================= %
% ======================================================================= %
% ======================================================================= %
%% Figure 2H: MODEL FREQUENCY AND EXCEEDANCE PROBABILITY:

% Select models:
modVec      = 1:5; 
modNames    = modVec;

% File to load:
modFile     = sprintf('hbi_mod%s', num2str(modVec, '_%02d')); % all models

% Directories:
dirs.results = fullfile(dirs.root, 'Log/Behavior/Modelling_CBM');
dirs.lap     = fullfile(dirs.results, 'LAP_Results');
dirs.hbi     = fullfile(dirs.results, 'HBI_Results');

% Add to directory to get final name:
fname_hbi   = fullfile(dirs.hbi, modFile); % use modName to get final name

% Load and extract:
f_hbi       = load([fname_hbi '.mat']);
cbm         = f_hbi.cbm;
nMod        = length(cbm.output.model_frequency);
modFreq     = cbm.output.model_frequency; % already normalized
PXP         = cbm.output.protected_exceedance_prob;

% ------------------------------------------------------ %
% Relative model frequency (normalized, divided by N) in red (blue);
% protected exceedance probability in red (right);
% small gap between bars;
% no SEs.

% Plot settings:
lineWidth   = 4; 
fontSize    = 36; % 25
xLoc        = 1:nMod;
xTickLabel  = strcat('M', string(modNames)); % actual names

% ------------------------------------------------------ %
% Start figure:

close all
p           = cell(nMod, 1);
figure('Position', [100 100 900 800], 'Color', 'white'); hold on

% Model frequency:
for iMod = 1:nMod
    p{iMod} = bar(xLoc(iMod) - 0.20, modFreq(iMod), .30, 'FaceColor', [0.1294 0.4000 0.6745]); % model frequency in blue
end

% Protected exceedance probability:
for iMod = 1:nMod
    p{iMod+nMod} = bar(xLoc(iMod) + 0.20, PXP(iMod), .30, 'FaceColor', [0.6980 0.0941 0.1686]); % PXP in red
end

% Add plot features:
set(gca, 'xlim', [0.5 nMod + 0.5] ,'ylim', [0 1], ...
    'xtick', xLoc, 'xticklabel', xTickLabel, ...
    'ytick', 0:.2:1, 'FontSize', fontSize, 'FontName', 'Arial', 'Linewidth', lineWidth);

% Labels:
xlabel('Model', 'FontSize', fontSize, 'FontName', 'Arial');

% Save:
figName = fullfile(dirs.root, 'Log/OutcomeLockedPaperPlots/Fig2H');
saveas(gcf, [figName '.png']);

% Close:
pause(3);
close gcf

% Save source data files:
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2H_modFreq.csv'));
csvwrite(fullFileName, modFreq);
fullFileName = fullfile(dirs.root, sprintf('Log/OutcomeLockedPaperPlots/Fig2H_PXP.csv'));
csvwrite(fullFileName, PXP);

% END OF FILE.