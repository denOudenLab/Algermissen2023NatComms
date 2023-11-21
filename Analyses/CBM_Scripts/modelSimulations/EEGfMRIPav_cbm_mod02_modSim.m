function [out] = EEGfMRIPav_cbm_mod02_modSim(parameters, subj)

% Standard Q-learning model with delta learning rule and Go bias.

% ----------------------------------------------------------------------- %
%% Retrieve parameters:

% 1) Feedback sensitivity:
nd_rho      = parameters(1); % normally-distributed rho
rho         = exp(nd_rho);

% 2) Learning rate:
nd_epsilon  = parameters(2);
epsilon     = 1 / (1 + exp(-nd_epsilon)); % epsilon (transformed to be between zero and one)

% 3) Go bias:
go_bias     = parameters(3);

% ----------------------------------------------------------------------- %
%% Unpack data:

% Extract task features:
stimuli     = subj.stimuli; % 1-16
reqactions  = subj.reqactions; % 1-3
feedback    = subj.feedback; % 0-1

% Data dimensions:
T       = length(stimuli); % number trials
C       = length(unique(stimuli)); % number stimuli
A       = length(unique(reqactions)); % number required actions

% To save outcome objects:
p       = nan(T, A);
Lik     = nan(T, 1);
PE      = nan(T, 1);
EV      = nan(T, A, C);

% Specific for model simulations (because probabilistic):
action  = nan(T, 1);
outcome = nan(T, 1);
stay    = nan(T, 1);

% For win-stay lose shift: count cumulative number presentations this
% stimulus:
cCount  = zeros(T, 1);
cRep    = nan(T, 1);

% Index whether valence of cue detected or not:
valenced= zeros(16, 1);

% Initialize Q-values:
q0      = repmat([1 1 -1 -1], 3, 4)/2;
q   	= q0*rho; % multiply with feedback sensitivity

% ----------------------------------------------------------------------- %
%% Loop over trials:

for t = 1:T    

    % Read info for the current trial:
    c   	= stimuli(t); % stimulus on this trial
    ra  	= reqactions(t); % required action on this trial
    f   	= feedback(t); % feedback validity on this trial
    v   	= q0(1, c); % valence of cue on this trial
    
    % Cumulative number presentations this stimulus:
    cCount(c)   = cCount(c) + 1; % increment count for this stimulus ID
    cRep(t)     = cCount(c); % store stimulus count

    % Retrieve Q-values of stimulus:
    w       = q(:, c) * valenced(c);
    
    % Add biases:
    w(1)    = w(1) + go_bias;
    w(2)    = w(2) + go_bias;
            
    % Softmax (turn Q-values into probabilities):
    pt 		= exp(w) / (exp(w(1)) + exp(w(2)) + exp(w(3)));
        
    % Choose actions:
    a       = randsample(3, 1, true, pt);
    
    % Determine outcome:
    if((a == ra && f == 1 && v > 0) || (a ~= ra && f == 0 && v > 0)); o = 1; end
    if((a == ra && f == 1 && v < 0) || (a ~= ra && f == 0 && v < 0)); o = 0; end
    if((a == ra && f == 0 && v > 0) || (a ~= ra && f == 1 && v > 0)); o = 0; end
    if((a == ra && f == 0 && v < 0) || (a ~= ra && f == 1 && v < 0)); o = -1; end

    % Update Q-values:
    delta    	= (rho * o) - q(a, c); % prediction error
    q(a, c)   	= q(a, c) + (epsilon * delta);    

    % Check if valence detected:
    if o ~= 0; valenced(c) = 1; end
    
    % Save objects for output:
    p(t, :)   	= pt;
    Lik(t)      = pt(a);
    action(t)   = a;
    outcome(t)  = o;
    PE(t)       = delta;
    EV(t, :, :) = q;
    
end

% ----------------------------------------------------------------------- %
%% Win stay-lose shift:

for t = 1:T % iTrial = 1;    
    c           = stimuli(t); % retrieve stimulus identifier
    repIdx      = cRep(t); % retrieve cumulative number of appearance of this stimulus
    nextT       = find(stimuli == c & cRep == (repIdx + 1)); % next trial with this stimulus
    if repIdx < max(cRep) % if there is still a "next trial" left
        stay(t) = double(action(t) == action(nextT)); % decide whether exact response (resp) repeated or not
    end
end

% ----------------------------------------------------------------------- %
%% Save as output object:

out.p       = p;
out.lik     = Lik;
out.PE      = PE;
out.EV      = EV;

out.action  = action;
out.outcome = outcome;
out.stay    = stay;

end % END OF FUNCTION.