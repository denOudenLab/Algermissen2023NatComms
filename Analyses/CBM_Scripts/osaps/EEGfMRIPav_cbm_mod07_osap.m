function [out] = EEGfMRIPav_cbm_mod07_osap(parameters, subj)

% Standard Q-learning model with delta learning rule and Go bias and
% Pavlovian response bias and Pavlovian learning bias and single
% perseveration parameter.

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

% 4) Pavlovian response bias:
pi_bias     = parameters(4);

% 5) Instrumental learning bias:
kappa_bias  = parameters(5);

% Transform:
biaseps     = nan(2, 1);
if epsilon < .5 % If default learning rate below 0.5
  biaseps(2)    = 1/ (1 + exp(-(nd_epsilon - kappa_bias))); % negative bias (Punishment after NoGo): subtract untransformed bias from untransformed epsilon, then transform
  biaseps(1)    = 2*epsilon - biaseps(2);                   % positive bias (Reward after Go): take difference transformed epsilon and transformed negative bias, add to transformed epsilon (= 2*transformed epsilon - transformed negative bias)
 else % If default learning rate above 0.5
  biaseps(1)    = 1 / (1 + exp(-(nd_epsilon + kappa_bias)));% positive bias (Reward after Go): add untransformed bias to untransformed epsilon, then transform
  biaseps(2)    = 2*epsilon - biaseps(1);                   % negative bias (Punishment after NoGo): take difference transformed epsilon and transformed positive bias, substract from transformed epsilon (= 2*transformed epsilon - transformed positive bias)
end

% 6) Choice perseveration parameter:
phi_bias    = parameters(6);

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

% Index of last action per cue:
lastAct = zeros(16, 1);

% Initialize Q-values:
q0      = repmat([1 1 -1 -1], 3, 4)/2;
q   	= q0*rho; % multiply with feedback sensitivity

% ----------------------------------------------------------------------- %
%% Loop over trials:

for t=1:T

    % Read info for the current trial:
    c    = stimuli(t); % stimulus on this trial
    a 	 = actions(t); % action on this trial
    o    = outcome(t); % outcome on this trial
    v    = q0(1,c); % valence of cue
    la   = lastAct(c); % last action to this cue

    % Retrieve Q-values of stimulus:
    w = q(:,c) * valenced(c);

    % Add biases:
    w(1) = w(1) + go_bias + valenced(c)*v*pi_bias;
    w(2) = w(2) + go_bias + valenced(c)*v*pi_bias;

    % Add choice perseveration:
    if la > 0 % if any last action available
         w(la)  = w(la) + phi_bias;
    end

    % Softmax (turn Q-values into probabilities):
    pt          = exp(w) / (exp(w(1)) + exp(w(2)) + exp(w(3)));

    % Update last action: 
    lastAct(c)  = a;
    
    % Determine learning rate:
    if ismember(a, [1, 2]) && o == 1
        eff_epsilon = biaseps(1);
    elseif a == 3 && o == -1
        eff_epsilon = biaseps(2);
    else
        eff_epsilon = epsilon;
    end
       
    % Update Q-values:
    delta   = (rho*o) - q(a, c); % prediction error
    q(a, c) = q(a, c) + (eff_epsilon*delta);    
        
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