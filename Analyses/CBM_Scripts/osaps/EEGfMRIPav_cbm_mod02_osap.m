function [out] = EEGfMRIPav_cbm_mod02_osap(parameters, subj)

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
stimuli = subj.stimuli; % 1-16
actions = subj.actions; % 1, 2, 3
outcome = subj.outcome; % 1, 0, -1

% Data dimensions:
T       = length(stimuli); % number trials
C       = length(unique(stimuli)); % number stimuli
A       = length(unique(reqactions)); % number required actions

% To save outcome objects:
p       = nan(T, A);
Lik     = nan(T, 1);
PE      = nan(T, 1);
EV      = nan(T, A, C);

% Index whether valence of cue detected or not:
valenced= zeros(16, 1);

% Initialize Q-values:
q0      = repmat([1 1 -1 -1], 3, 4)/2;
q   	= q0*rho; % multiply with feedback sensitivity

% ----------------------------------------------------------------------- %
%% Loop over trials:

for t = 1:T    

    % Read info for the current trial:
    c   = stimuli(t); % stimulus on this trial
    a 	= actions(t); % action on this trial
    o   = outcome(t); % outcome on this trial

    % Retrieve Q-values of stimulus:
    w       = q(:, c) * valenced(c);

    % Add biases:
    w(1)    = w(1) + go_bias;
    w(2)    = w(2) + go_bias;    

     % Softmax (turn Q-values into probabilities):
    pt      = exp(w) / (exp(w(1)) + exp(w(2)) + exp(w(3)));
       
    % Update Q-values:
    delta   = (rho*o) - q(a, c); % prediction error
    q(a, c) = q(a, c) + (epsilon*delta);    
        
    % Check if valence detected:
    if o ~= 0; valenced(c) = 1; end
    
    % Save objects for output:
    p(t, :)     = pt;
    Lik(t)      = pt(a);
    PE(t)       = delta;
    EV(t, :, :) = q;

end

% ----------------------------------------------------------------------- %
%% Save as output object:

out.p       = p;
out.lik     = Lik;
out.PE      = PE;
out.EV      = EV;

end % END OF FUNCTION.