function yHat = setar_forecast(y, model, horizon)
%SETAR_FORECAST Produce iterative forecasts from an estimated SETAR model.
%   yHat = SETAR_FORECAST(y, model, horizon) returns an array of size
%   horizon-by-1 containing forecasts for steps 1:horizon ahead obtained by
%   iterating on the provided data vector y (latest observation last).
%
%   MODEL must be a struct returned by FIT_SETAR_LS and contains the fields
%   .threshold, .coeffLow, .coeffHigh. Forecasts use the last max(pLow,pHigh)
%   observations and the same delay as during estimation. The struct must
%   also contain .delay; if absent the function assumes delay = 1.

if ~isfield(model, 'delay')
    delay = 1;
else
    delay = model.delay;
end

pLow = numel(model.coeffLow) - 1;
pHigh = numel(model.coeffHigh) - 1;
p = max(pLow, pHigh);

buffer = y(:);
assert(numel(buffer) >= max(p, delay), ...
    'setar_forecast:InsufficientHistory', ...
    'Need at least %d past values.', max(p, delay));

buffer = buffer(:);
yHat = zeros(horizon, 1);

for h = 1:horizon
    driver = buffer(end - delay + 1);
    if driver <= model.threshold
        theta = model.coeffLow;
        order = pLow;
    else
        theta = model.coeffHigh;
        order = pHigh;
    end

    if order == 0
        yForecast = theta(1);
    else
        past = buffer(end:-1:end - order + 1);
        yForecast = theta(1) + theta(2:end) * past;
    end
    yHat(h) = yForecast;
    buffer = [buffer; yForecast]; %#ok<AGROW>
end
end
