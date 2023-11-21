function [loglik] = EEGfMRIPav_cbm_mod08(parameters,subj)

% Standard Q-learning model with delta learning rule and with Go bias,
% Pavlovian response bias, instrumental learning bias, and separate choice
% perseveration parameters for Win and Avoid cues.

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
  biaseps(2) = 1 / (1 + exp(-(nd_epsilon - kappa_bias)));   % negative bias (Punishment after NoGo): subtract untransformed bias from untransformed epsilon, then transform
  biaseps(1) = 2 * epsilon - biaseps(2);                    % positive bias (Reward after Go): take difference transformed epsilon and transformed negative bias, add to transformed epsilon (= 2*transformed epsilon - transformed negative bias)
 else           % If default learning rate above 0.5
  biaseps(1) = 1 / (1 + exp(-(nd_epsilon + kappa_bias)));   % positive bias (Reward after Go): add untransformed bias to untransformed epsilon, then transform
  biaseps(2) = 2 * epsilon - biaseps(1);                    % negative bias (Punishment after NoGo): take difference transformed epsilon and transformed positive bias, substract from transformed epsilon (= 2*transformed epsilon - transformed positive bias)
end

% 6) Choice perseveration parameters:
phi_Win     = parameters(6);
phi_Avoid   = parameters(7);

% ----------------------------------------------------------------------- %
%% Unpack data:

stimuli = subj.stimuli; % 1-16
actions = subj.actions; % 1, 2, 3
outcome = subj.outcome; % 1, 0, -1

% Number of trials:
T       = size(outcome, 1);

% To save probability of choice. Currently NaNs, will be filled below:
p       = nan(T, 1);

% Index whether valence of cue detected or not:
valenced= zeros(16, 1);

% Index of last action per cue:
lastAct = zeros(16, 1);

% Initialize Q-values:
q0      = repmat([1 1 -1 -1], 3, 4)/2;
q       = q0*rho; % multiply with feedback sensitivity

% ----------------------------------------------------------------------- %
%% Loop over trials:

for t = 1:T    
    
    % Read info for the current trial:
    c    = stimuli(t); % stimulus on this trial
    a    = actions(t); % action on this trial
    o    = outcome(t); % outcome on this trial
    v    = q0(1, c); % valence of cue
    la   = lastAct(c); % last action to this cue

    % Retrieve Q-values of stimulus:
    w       = q(:, c) * valenced(c);
    
    % Add biases:
    w(1)    = w(1) + go_bias + valenced(c)*v*pi_bias;
    w(2)    = w(2) + go_bias + valenced(c)*v*pi_bias;

    % Add choice perseveration:
    if la > 0 % if any last action available
         if v > 0 % for Win cues:
             w(la) = w(la) + phi_Win;
         else % for Avoid cues:
             w(la) = w(la) + phi_Avoid;
         end
    end
    
    % Softmax (turn Q-values into probabilities):
    pt          = exp(w) / (exp(w(1)) + exp(w(2)) + exp(w(3)));
    
    % Store probability of the chosen action:
    p(t)        = pt(a);

    % Update last action for this cue: 
    lastAct(c)  = a;
    
    % ------------------------------------------------------------------- %
    % Update Q-values:
    
    % Select learning rate:
    if ismember(a, [1, 2]) && o == 1
        eff_epsilon = biaseps(1);
    elseif a == 3 && o == -1
        eff_epsilon = biaseps(2);
    else
        eff_epsilon = epsilon;
    end
    
    % Update:
    delta    = (rho*o) - q(a, c); % prediction error
    q(a, c)  = q(a, c) + (eff_epsilon*delta);

    % Check if valence detected:
    if o ~= 0; valenced(c) = 1; end
    
end

% ----------------------------------------------------------------------- %
%% Compute log-likelihood (sum of log-probability of choice data given the parameters):

loglik  = sum(log(p + eps));

end % END OF FUNCTION.