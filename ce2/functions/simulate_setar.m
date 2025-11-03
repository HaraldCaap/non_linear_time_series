function [y, regimes] = simulate_setar(settings)
%SIMULATE_SETAR Simulates a SETAR(2; p_L, p_H) time series.
%   SETTINGS is a struct with fields:
%       n           - number of observations (excluding burn-in)
%       burn        - burn-in iterations discarded from the start
%       delay       - threshold delay d
%       threshold   - scalar threshold r
%       sigma       - innovation standard deviation
%       regimes     - cell array with two parameter vectors for each regime.
%                     Each parameter vector must be [c phi_1 ... phi_p].
%       init        - optional vector of initial conditions. If omitted the
%                     largest lag is filled with zeros.
%
%   Returns:
%       y        - simulated series of length SETTINGS.n
%       regimes  - indicator (1 or 2) of the active regime per observation
%
%   The routine supports different autoregressive orders in the two regimes.

requiredFields = {'n','burn','delay','threshold','sigma','regimes'};
for k = 1:numel(requiredFields)
    assert(isfield(settings, requiredFields{k}), ...
        'simulate_setar:MissingField', 'Missing field "%s".', requiredFields{k});
end

n          = settings.n;
burn       = settings.burn;
delay      = settings.delay;
r          = settings.threshold;
sigma      = settings.sigma;
regParams  = settings.regimes;

assert(iscell(regParams) && numel(regParams) == 2, ...
    'simulate_setar:BadRegimes', 'regimes must be a 1x2 cell array.');

pL = numel(regParams{1}) - 1;
pH = numel(regParams{2}) - 1;
p  = max([pL, pH, delay]);

nTotal = n + burn;
y = zeros(nTotal + p, 1);
regimes = ones(nTotal + p, 1);

if isfield(settings, 'init') && ~isempty(settings.init)
    init = settings.init(:);
    assert(numel(init) >= p, 'Initial vector too short.');
    y(1:p) = init(end-p+1:end);
end

for t = p + 1 : nTotal + p
    thresholdDriver = y(t - delay);
    if thresholdDriver <= r
        theta = regParams{1};
        arOrder = pL;
        regimes(t) = 1;
    else
        theta = regParams{2};
        arOrder = pH;
        regimes(t) = 2;
    end
    if arOrder == 0
        deterministic = theta(1);
    else
        lagged = y(t-1:-1:t-arOrder);
        deterministic = theta(1) + theta(2:end) * lagged;
    end
    y(t) = deterministic + sigma * randn();
end

% Drop burn-in and leading padding
startIdx = p + burn + 1;
y = y(startIdx:startIdx + n - 1);
regimes = regimes(startIdx:startIdx + n - 1);
end
