%% Part 2 — Least squares estimation of the SETAR model
% We take the simulated data from Part 1 and perform a grid-search least
% squares fit over potential thresholds. The best model (in terms of SSE) is
% saved to disk for reuse in later parts.

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

%% Load simulated data
[folder, ~, ~] = fileparts(mfilename('fullpath'));
dataFile = fullfile(folder, 'data', 'part1_sim.mat');
assert(exist(dataFile, 'file') == 2, 'Run part1.m before part2.m');
S = load(dataFile);
y = S.y;

%% Estimation configuration
pLow = 1; pHigh = 1; delay = 1;

gridSize = 40;
threshGrid = quantile(y, linspace(0.10, 0.90, gridSize));
threshGrid = unique(threshGrid);  % remove duplicates caused by ties

%% Fit model
model = fit_setar_ls(y, pLow, pHigh, delay, threshGrid);
model.delay = delay;

%% Persist results
modelFile = fullfile(folder, 'results', 'part2_model.mat');
save(modelFile, 'model');

summaryFile = fullfile(folder, 'results', 'part2_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'SETAR(2;1,1) estimation summary\n');
fprintf(fid, 'Selected threshold: %.4f\n', model.threshold);
fprintf(fid, 'Low regime coefficients: [%.4f %.4f]\n', model.coeffLow);
fprintf(fid, 'High regime coefficients: [%.4f %.4f]\n', model.coeffHigh);
fprintf(fid, 'Sum of squared errors: %.4f\n', model.sse);
fprintf(fid, 'BIC: %.4f\n', model.bic);
fprintf(fid, 'Observations used: %d\n', model.nobs);
fclose(fid);

%% Display summary in console
disp('--- Part 2: Estimation complete ---');
disp(fileread(summaryFile));

%% Diagnostic plot of grid search history
figDir = fullfile(folder, 'figs');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

validMask = ~isnan(model.history(:,2));
figure('Name', 'Threshold grid search diagnostics', 'Color', 'w');
plot(model.history(validMask,1), model.history(validMask,2), 'o-', ...
    'LineWidth', 1.2); hold on;
yline(model.sse, 'r--', 'LineWidth', 1.2);
grid on;
xlabel('Threshold candidate r');
ylabel('Sum of squared errors');
title('SETAR threshold grid search');
legend({'Grid SSE', 'Selected SSE'}, 'Location', 'best');

print('-dpng', fullfile(figDir, 'part2_grid_search.png'));
