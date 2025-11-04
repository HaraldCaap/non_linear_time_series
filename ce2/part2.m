%% SETAR(2;1;1)
N = 3000;
e = randn(N,1);
x = zeros(N,1);

% Parameters 
r1 = 0; % threshold
a0 = [-1.0, 1.0]; % coefficients a_0^J
a1 = [0.6, 0.4]; % coefficients a_1^J

% Initialize the process depending on the regime
if 0 <= r1
    x(1) = sqrt(1 / (1 - a1(1)^2)) * e(1);
else
    x(1) = sqrt(1 / (1 - a1(2)^2)) * e(1);
end

% Simulate the process
for n=2:N
    if x(n-1) <= r1
        x(n) = a0(1) + a1(1)*x(n-1) + e(n);
    else
        x(n) = a0(2) + a1(2)*x(n-1) + e(n);
    end
end

% Plot
figure;
plot(x)
yline(r1, 'k--')
xlabel("Time")
ylabel("x_t")
title("Simulation of SETAR(2;1;1)")

%%
r = 0;
d = 1;

Q = @(p, obs_low, obs_high) sum((x(d+1+obs_low-1:obs_high) ...
        - (p(1) + 0.6.*x(d+obs_low-1:obs_high-1)).*(x(1+obs_low-1:obs_high-d)<=r) ...
        - (p(2) + 0.4.*x(d+obs_low-1:obs_high-1)).*(x(1+obs_low-1:obs_high-d)>r)).^2);

t = linspace(-2, 2, 100);
s = linspace(-2, 2, 100);
[T, S] = meshgrid(t, s);

obs_lows = [1, 1, 1, 1001, 1001];
obs_highs = [3000, 300, 30, 1300, 1030];

figure;
tiled = tiledlayout(3, 2);
for i=1:length(obs_highs)
    Z = arrayfun(@(p1, p2) Q([p1, p2], obs_lows(i), obs_highs(i)), T, S);
    nexttile;
    imagesc(t, s, Z);
    axis xy;
    colorbar;
    hold on
    contour(T, S, Z, 40, 'k')  % 20 contour levels
    hold off
    xlabel('p1')
    ylabel('p2')
    title("Observations: " + obs_lows(i) + ":" + obs_highs(i));
end
tiled.TileSpacing = 'compact';
tiled.Padding = 'compact';