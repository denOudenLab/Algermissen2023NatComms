function output = regress_baseline(input, x, startIdx, stopIdx, corMode)

% output = regress_baseline(input, x, startIdx, stopIdx, corMode)
%
% Remove baseline by performing regression over trials.
%
% INPUTSs:
% input         = matrix of numerics of n-dimensional format. 
%   First dimension is trials.
%   Arbitary number of intermediate dimensions (e.g. freq) collapsed into
%   separate bins.
%   Last dimension is time.
% x             = vector of numerics, single regressor in regression --
% indices of trials provided (default: 1 until number of trials provided).
% startIdx      = scalar integer, time index (last dimension of input
% matrix) from when onwards to do baseline extraction.
% stopIdx       =  scalar integer, time index (last dimension of input
% matrix) until when to do baseline extraction.
% corMode       = scalar string, type of baseline correction, either
% 'subtraction' or 'division'.
% 
% OUTPUT:
% output        = input matrix corrected for baseline
%
% RECOMMENDATIONS:
% Use 'subtraction' type baseline correction for time-domain data.
% Use 'division' type baseline correction for frequency-domain data,
% perform decibel-conversion (10*log10(output)) afterwards.

% ----------------------------------------------------------------------- %
%% Complete input settings:

if nargin < 2
    x           = 1:size(input, 1);
    fprintf('No trial indices for regressor provided---assume 1 until number of input trials\n');
end
if nargin < 3
    startIdx    = 1;
    fprintf('No start index for baseline extraction provided---assume 1\n');
end
if nargin < 4
    stopIdx     = 1;
    fprintf('No stop index for baseline extraction provided---assume 1\n');
end
if nargin < 5
    corMode     = 'subtraction';
    fprintf('No mode of baseline correction provided---perform subtraction\n');
end

% ----------------------------------------------------------------------- %
%% Check if input is cell:

if iscell(input)
    fprintf('Input is of type cell, reshape into matrix\n')
    nDim  = length(size(input{1})) + 1; % Determine which would be new extra dimension
    input = permute(cat(nDim, input{:}), [nDim 1:(nDim-1)]); % concatenate into final dimension, bring final dimension to start
end

% ----------------------------------------------------------------------- %
%% Assess matrix dimensions:

dimVec  = size(input);

nTrial  = dimVec(1); % first dimension: number trials
nTime   = dimVec(end); % last dmension: number time bins
nRest   = prod(dimVec) / nTrial / nTime; % number of any other bins in the middle

% Check x: must be column vector:
if size(x, 1) == 1
    x = x';
end    

% ----------------------------------------------------------------------- %
%% Reshape into 3D: trial x rest x time:

input_3D    = reshape(input, nTrial, nRest, nTime);

% ----------------------------------------------------------------------- %
%% Perform regression over trials:

fprintf('Perform regression-based baseline-correction over %d trials for %d separate bins\n', nTrial, nRest);

% Initialize empty template for baseline:
baseline    = nan(nTrial, nRest, nTime);

% ----------------------------------------------------------------------- %
%% Loop through nRest, extract baseline in specified window, fit regression over trials:

for iRest = 1:nRest % iRest = 1; % loop over all other dimensions (not trials or time):
       % Extract baseline from specified window, average:
       y = nanmean(input_3D(:, iRest, startIdx:stopIdx), 3);  % average in baselineTimings
       % Fit regression:
       p = polyfit(x, y, 1); % fit baseline as function of trial number
       % Predict baseline: 
       f = polyval(p, x); % predict values from regression for x
       % Fill into baseline template:
       baseline(:, iRest, :) = repmat(f, 1, 1, nTime); % save for all trials
end % end of iRest

% Subtract fitted baseline from actual data:
if strcmp(corMode,'subtraction')
    output_3D = input_3D - baseline; % element-wise subtraction
elseif strcmp(corMode,'division')
    output_3D = input_3D ./ baseline; % element-wise division
else
    error('Unknown correction mode');
end

% ----------------------------------------------------------------------- %
%% Reshape back into original format:

output = reshape(output_3D, size(input));

end % END OF FUNCTION.