%% Part 3 — Regime diagnostics and ergodic properties
% Using the estimated SETAR model we analyse the regime-specific behaviour
% and approximate long-run moments via simulation.

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

[folder, ~, ~] = fileparts(mfilename('fullpath'));

dataFile = fullfile(folder, 'data', 'part1_sim.mat');
modelFile = fullfile(folder, 'results', 'part2_model.mat');
assert(exist(dataFile, 'file') == 2, 'Run part1.m first');
assert(exist(modelFile, 'file') == 2, 'Run part2.m first');

Sdata = load(dataFile);
y = Sdata.y;
Smodel = load(modelFile);
model = Smodel.model;

%% Assign regimes according to the estimated threshold
T = numel(y);
pLow = numel(model.coeffLow) - 1;
pHigh = numel(model.coeffHigh) - 1;
delay = model.delay;
maxLag = max([pLow, pHigh, delay]);

regimeIdx = nan(T, 1);
regimeIdx(1:maxLag) = NaN; %#ok<NASGU>
for t = maxLag + 1:T
    if y(t - delay) <= model.threshold
        regimeIdx(t) = 1;
    else
        regimeIdx(t) = 2;
    end
end

validMask = ~isnan(regimeIdx);

%% Empirical regime statistics
stats.low.mean  = mean(y(regimeIdx == 1));
stats.low.var   = var(y(regimeIdx == 1));
stats.low.count = sum(regimeIdx == 1);

stats.high.mean  = mean(y(regimeIdx == 2));
stats.high.var   = var(y(regimeIdx == 2));
stats.high.count = sum(regimeIdx == 2);

%% Stability check (contractive slopes)
phiLow  = model.coeffLow(2:end);
phiHigh = model.coeffHigh(2:end);
stability.low  = all(abs(phiLow)  < 1);
stability.high = all(abs(phiHigh) < 1);

%% Monte Carlo approximation of ergodic moments using the estimated model
mc.nsim  = 5e4;
mc.burn  = 5e3;
mc.delay = delay;
mc.threshold = model.threshold;
mc.sigma = sqrt(model.sse / model.nobs);
mc.regimes = {model.coeffLow, model.coeffHigh};

[mcSeries, mcRegime] = simulate_setar(struct('n', mc.nsim, ...
    'burn', mc.burn, ...
    'delay', mc.delay, ...
    'threshold', mc.threshold, ...
    'sigma', mc.sigma, ...
    'regimes', mc.regimes));

longRun.mean = mean(mcSeries);
longRun.var  = var(mcSeries);
longRun.probLow  = mean(mcRegime == 1);
longRun.probHigh = mean(mcRegime == 2);

%% Save report
reportFile = fullfile(folder, 'results', 'part3_report.txt');
fid = fopen(reportFile, 'w');
fprintf(fid, 'Regime-specific sample moments (original data)\n');
fprintf(fid, 'Regime 1 mean: %.4f  variance: %.4f  count: %d\n', ...
    stats.low.mean, stats.low.var, stats.low.count);
fprintf(fid, 'Regime 2 mean: %.4f  variance: %.4f  count: %d\n', ...
    stats.high.mean, stats.high.var, stats.high.count);
fprintf(fid, '\nStability (|phi| < 1): Regime 1 = %d, Regime 2 = %d\n', ...
    stability.low, stability.high);
fprintf(fid, '\nMonte Carlo long-run estimates (estimated model)\n');
fprintf(fid, 'Mean: %.4f  Variance: %.4f\n', longRun.mean, longRun.var);
fprintf(fid, 'Pr(Regime 1): %.4f  Pr(Regime 2): %.4f\n', ...
    longRun.probLow, longRun.probHigh);
fclose(fid);

disp('--- Part 3 summary ---');
disp(fileread(reportFile));

%% Visualise empirical vs simulated distribution
figDir = fullfile(folder, 'figs');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

figure('Name', 'Part 3 distributions', 'Color', 'w');
hold on; grid on; box on;
histogram(y(validMask), 'Normalization', 'pdf', 'DisplayName', 'Observed');
histogram(mcSeries, 'Normalization', 'pdf', 'DisplayName', 'Simulated');
xlabel('$y$', 'Interpreter', 'latex'); ylabel('Density');
title('Empirical vs simulated stationary distribution', 'Interpreter', 'latex');
legend('Location', 'best');

print('-dpng', fullfile(figDir, 'part3_distribution.png'));
