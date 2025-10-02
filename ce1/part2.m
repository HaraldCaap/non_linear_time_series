close all
clear all

%% Simulation
r1 = 0;
rng(1)
N = 1000;
e = randn(N,1);
X = zeros(N,1);

a0 = [1.0, -1.0];
a1 = [0.6, 0.4];

% Initialize the first value of X based on r1
X(1) = sqrt(1 / (1 - a1(1)^2)) * e(1) * (r1 >= 0) + ...
       sqrt(1 / (1 - a1(2)^2)) * e(1) * (r1 < 0);

% Generate the time series
for n = 2:N
    X(n) = (a0(1) + a1(1) * X(n-1) + e(n)) * (X(n-1) <= r1) + ...
           (a0(2) + a1(2) * X(n-1) + e(n)) * (X(n-1) > r1);
end

figure(1)
plot(X)
title("SETAR(2;1;1)")

% Theoretical conditional expectation
x_theo = linspace(min(X), max(X), 100);
Mhat_theo = @(x) (a0(1) + a1(1) * x) .* (x <= r1) + ...
                 (a0(2) + a1(2) * x) .* (x > r1);

% Calculate an estimation of the conditional expectation
h = [0.1, 0.4, 0.8];
points = 100;
order = 1;
Mhat = zeros(length(h), points);
for i = 1:length(h)
    [x, Mhat(i,:)] = regsmooth1D([X(1:end-1) X(2:end)], points, order, h(i));
end

figure(2)
hold on
p1 = plot(x, Mhat, 'LineWidth', 2);
s1 = scatter(x_theo, Mhat_theo(x_theo), 'k+');
scatter(X(1:end-1), X(2:end), 'b', 'filled', 'MarkerFaceAlpha', 0.1);
legend("h=0.1", "h=0.4", "h=0.8", "Theoretical conditional mean", "Realizations")
xlabel("X_k")
ylabel("X_{k+1}")
title("Theoretical and estimated conditional mean")
grid on
hold off