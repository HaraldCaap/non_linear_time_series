%% Part 1 — Simulating a SETAR(2;1,1) process
% This script reproduces the simulation step from Computer Exercise 1.
% We generate a SETAR(2;1,1) series with known parameters, visualize the
% regimes, and persist the data for later parts of the exercise.

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

%% Settings
disp('--- Part 1: Simulating SETAR(2;1,1) ---');

rng(20240513, 'twister');

params.n         = 1200;
params.burn      = 300;
params.delay     = 1;
params.threshold = 0.2;
params.sigma     = 0.5;
params.regimes   = { [0.4  0.6], ...   % Regime 1: y_t = 0.4 + 0.6 y_{t-1} + eps_t
                      [-0.3 0.9] };    % Regime 2: y_t = -0.3 + 0.9 y_{t-1} + eps_t

[y, regime] = simulate_setar(params);

%% Persist data for later parts
dataPath = fullfile(fileparts(mfilename('fullpath')), 'data', 'part1_sim.mat');
save(dataPath, 'y', 'regime', 'params');

fprintf('Saved simulated series to %s\n', dataPath);

%% Visualisation
figDir = fullfile(fileparts(mfilename('fullpath')), 'figs');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

t = (1:numel(y))';
figure('Name', 'Part 1: Simulated SETAR time series', 'Color', 'w');
subplot(2,1,1);
plot(t, y, 'k-', 'LineWidth', 1.0); grid on;
xlabel('t'); ylabel('$y_t$', 'Interpreter', 'latex');
title('Simulated SETAR(2;1,1) series', 'Interpreter', 'latex');

subplot(2,1,2);
stairs(t, regime, 'LineWidth', 1.0); grid on;
ylabel('Regime indicator'); xlabel('t');
yticks([1 2]); yticklabels({'Regime 1', 'Regime 2'});
title('Active regime (delay d = 1)', 'Interpreter', 'latex');

print('-dpng', fullfile(figDir, 'part1_simulation.png'));

disp('Part 1 completed.');
