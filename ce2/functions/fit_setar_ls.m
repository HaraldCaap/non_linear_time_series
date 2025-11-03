function model = fit_setar_ls(y, pLow, pHigh, delay, thresholds)
%FIT_SETAR_LS Grid-search least squares estimation for SETAR(2; p_L, p_H).
%   MODEL = FIT_SETAR_LS(Y, pLow, pHigh, delay, thresholds) evaluates each
%   candidate threshold in the vector THRESHOLDS and returns the model with
%   the smallest sum of squared errors. Y must be a column vector.
%
%   The returned struct contains:
%       .threshold      - selected threshold
%       .coeffLow       - [c phi_1 ... phi_pLow]
%       .coeffHigh      - [c phi_1 ... phi_pHigh]
%       .residuals      - in-sample residuals
%       .regimeIndices  - regime indicator for each observation
%       .sse            - total sum of squared errors
%       .bic            - Bayesian information criterion
%       .nobs           - effective sample size used for estimation
%
%   Observations are dropped until the largest lag and delay are available.

assert(isvector(y), 'Input y must be a vector.');
y = y(:);

maxLag = max([pLow, pHigh, delay]);
T = numel(y);
startIdx = maxLag + 1;

Y = y(startIdx:end);
regDriver = y(startIdx - delay:end - delay);

XLow  = buildDesignMatrix(y, pLow, startIdx);
XHigh = buildDesignMatrix(y, pHigh, startIdx);

bestModel = struct('sse', inf);
thresholds = thresholds(:);
history = nan(numel(thresholds), 3);
idx = 0;

for r = thresholds'
    regimeLow = regDriver <= r;
    regimeHigh = ~regimeLow;

    idx = idx + 1;

    if sum(regimeLow) <= pLow || sum(regimeHigh) <= pHigh
        history(idx, :) = [r, NaN, NaN];
        continue; % Skip thresholds yielding too few observations
    end

    [betaLow, resLow] = olsSolve(XLow(regimeLow, :), Y(regimeLow));
    [betaHigh, resHigh] = olsSolve(XHigh(regimeHigh, :), Y(regimeHigh));

    sse = sum(resLow.^2) + sum(resHigh.^2);
    k = (pLow + 1) + (pHigh + 1) + 1; % parameters + threshold
    nEff = numel(Y);
    bic = nEff * log(sse / nEff) + k * log(nEff);
    history(idx, :) = [r, sse, bic];

    if sse < bestModel.sse
        bestModel.threshold = r;
        bestModel.coeffLow = betaLow(:)';
        bestModel.coeffHigh = betaHigh(:)';
        bestModel.sse = sse;
        bestModel.bic = bic;
        bestModel.nobs = nEff;
        residuals = zeros(size(Y));
        residuals(regimeLow) = resLow;
        residuals(regimeHigh) = resHigh;
        bestModel.residuals = residuals;
        regimes = ones(size(Y)) * 2;
        regimes(regimeLow) = 1;
        bestModel.regimeIndices = regimes;
    end
end

if ~isfield(bestModel, 'threshold')
    error('fit_setar_ls:NoValidThreshold', ...
          'No valid threshold produced enough observations per regime.');
end

model = bestModel;
model.history = history(1:idx, :);
end

function X = buildDesignMatrix(y, p, startIdx)
T = numel(y) - startIdx + 1;
if p == 0
    X = ones(T, 1);
    return;
end
X = ones(T, p + 1);
for j = 1:p
    X(:, j + 1) = y(startIdx - j:end - j);
end
end

function [beta, residuals] = olsSolve(X, y)
XtX = X' * X;
beta = XtX \ (X' * y);
residuals = y - X * beta;
end
