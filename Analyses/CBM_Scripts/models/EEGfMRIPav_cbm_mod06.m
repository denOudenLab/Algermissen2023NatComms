function [loglik] = EEGfMRIPav_cbm_mod06(parameters,subj)

% Standard Q-learning model with delta learning rule and with Go bias,
% Pavlovian response bias, and separate learning rate for salient Gos
% (action priming model).

% ----------------------------------------------------------------------- %
%% Retrieve parameters:

% 1) Feedback sensitivity:
nd_rho      = parameters(1); % normally-distributed rho
rho         = exp(nd_rho);

% 2) Learning rate:
nd_epsilon  = parameters(2);
epsilon     = 1 / (1 + exp(-nd_epsilon)); % (transformed to be between zero and one)

% 3) Go bias:
go_bias     = parameters(3);

% 4) Pavlovian response bias:
pi_bias     = parameters(4);

% 5) Learning rate for Go actions with salient outcomes:
nd_epsilon_SalGo = parameters(5);
epsilon_SalGo = 1 / (1 + exp(-nd_epsilon_SalGo)); % (transformed to be between zero and one)

% ----------------------------------------------------------------------- %
%% Unpack data:

stimuli = subj.stimuli; % 1-16
actions = subj.actions; % 1 ,2, 3
outcome = subj.outcome; % 1, 0, -1

% Number of trials:
T       = size(outcome, 1);

% To save probability of choice. Currently NaNs, will be filled below:
p       = nan(T, 1);

% Index whether valence of cue detected or not:
valenced= zeros(16, 1);

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

    % Retrieve Q-values of stimulus:
    w       = q(:, c) * valenced(c);
    
    % Add biases:
    w(1)    = w(1) + go_bias + valenced(c)*v*pi_bias;
    w(2)    = w(2) + go_bias + valenced(c)*v*pi_bias;
    
    % Softmax (turn Q-values into probabilities):
    pt      = exp(w) / (exp(w(1)) + exp(w(2)) + exp(w(3)));
    
    % Store probability of the chosen action:
    p(t)    = pt(a);
        
    % Select learning rate:
    if ismember(a, [1, 2]) && ismember(o, [1, -1])
        eff_epsilon = epsilon_SalGo;
    else
        eff_epsilon = epsilon;
    end
    
    % Update:
    delta   = (rho*o) - q(a,c); % prediction error
    q(a, c) = q(a,c) + (eff_epsilon*delta);

    % Check if valence detected:
    if o ~= 0; valenced(c) = 1; end
    
end

% ----------------------------------------------------------------------- %
%% Compute log-likelihood (sum of log-probability of choice data given the parameters):

loglik  = sum(log(p + eps));

end % END OF FUNCTION.