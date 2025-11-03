%% Part 5 — Residual bootstrap for threshold uncertainty
% We assess parameter uncertainty via a simple residual bootstrap. The
% procedure resamples centered residuals from the fitted SETAR model,
% generates synthetic series, re-estimates the model, and summarises the
% sampling distribution of key parameters.

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

[folder, ~, ~] = fileparts(mfilename('fullpath'));
dataFile = fullfile(folder, 'data', 'part1_sim.mat');
modelFile = fullfile(folder, 'results', 'part2_model.mat');
assert(exist(dataFile, 'file') == 2 && exist(modelFile, 'file') == 2, ...
    'Run part1.m and part2.m before part5.m');

Sdata = load(dataFile);
y = Sdata.y;
Smodel = load(modelFile);
model = Smodel.model;

pLow = numel(model.coeffLow) - 1;
pHigh = numel(model.coeffHigh) - 1;
delay = model.delay;
maxLag = max([pLow, pHigh, delay]);

residuals = model.residuals(:);
residuals = residuals - mean(residuals); % recentre
nEff = numel(residuals);

B = 200; % number of bootstrap replications
threshGrid = model.history(:,1);
threshGrid = unique(threshGrid(~isnan(threshGrid)));

bootThreshold = nan(B, 1);
bootCoeffLow = nan(B, numel(model.coeffLow));
bootCoeffHigh = nan(B, numel(model.coeffHigh));

fprintf('Running %d bootstrap replications...\n', B);

for b = 1:B
    yBoot = zeros(size(y));
    yBoot(1:maxLag) = y(1:maxLag); % seed with observed history
    for t = maxLag + 1:numel(y)
        driver = yBoot(t - delay);
        if driver <= model.threshold
            theta = model.coeffLow;
            order = pLow;
        else
            theta = model.coeffHigh;
            order = pHigh;
        end
        if order == 0
            det = theta(1);
        else
            past = yBoot(t-1:-1:t-order);
            det = theta(1) + theta(2:end) * past;
        end
        innov = residuals(randi(nEff));
        yBoot(t) = det + innov;
    end

    try
        bootModel = fit_setar_ls(yBoot, pLow, pHigh, delay, threshGrid);
        bootThreshold(b) = bootModel.threshold;
        bootCoeffLow(b, :) = bootModel.coeffLow;
        bootCoeffHigh(b, :) = bootModel.coeffHigh;
    catch ME %#ok<NASGU>
        % leave NaNs for failed replications
    end
end

fprintf('Bootstrap finished. Effective draws: %d\n', sum(~isnan(bootThreshold)));

%% Summaries
valid = ~isnan(bootThreshold);
ciLevel = [0.05 0.95];
threshCI = quantile(bootThreshold(valid), ciLevel);
coeffLowCI = quantile(bootCoeffLow(valid, :), ciLevel);
coeffHighCI = quantile(bootCoeffHigh(valid, :), ciLevel);

reportFile = fullfile(folder, 'results', 'part5_bootstrap_summary.txt');
fid = fopen(reportFile, 'w');
fprintf(fid, 'Bootstrap threshold mean: %.4f\n', mean(bootThreshold(valid)));
fprintf(fid, '95%% interval: [%.4f, %.4f]\n', threshCI(1), threshCI(2));

fprintf(fid, '\nRegime 1 coefficients (rows: lower/upper CI)\n');
fprintf(fid, '%.4f %.4f\n', coeffLowCI(1, :));
fprintf(fid, '%.4f %.4f\n', coeffLowCI(2, :));

fprintf(fid, '\nRegime 2 coefficients (rows: lower/upper CI)\n');
fprintf(fid, '%.4f %.4f\n', coeffHighCI(1, :));
fprintf(fid, '%.4f %.4f\n', coeffHighCI(2, :));
fclose(fid);

disp('--- Part 5 bootstrap summary ---');
disp(fileread(reportFile));

%% Visualise bootstrap distributions
figDir = fullfile(folder, 'figs');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

figure('Name', 'Part 5 bootstrap distributions', 'Color', 'w');
subplot(3,1,1);
histogram(bootThreshold(valid), 20, 'FaceColor', [0.2 0.4 0.8]);
grid on; xlabel('Threshold r'); ylabel('Frequency');
title('Bootstrap threshold distribution', 'Interpreter', 'latex');

subplot(3,1,2);
histogram(bootCoeffLow(valid,2), 20, 'FaceColor', [0.8 0.3 0.3]);
grid on; xlabel('$\phi_{1}^{(1)}$', 'Interpreter', 'latex'); ylabel('Frequency');
title('Regime 1 AR coefficient', 'Interpreter', 'latex');

subplot(3,1,3);
histogram(bootCoeffHigh(valid,2), 20, 'FaceColor', [0.3 0.7 0.3]);
grid on; xlabel('$\phi_{1}^{(2)}$', 'Interpreter', 'latex'); ylabel('Frequency');
title('Regime 2 AR coefficient', 'Interpreter', 'latex');

print('-dpng', fullfile(figDir, 'part5_bootstrap.png'));
