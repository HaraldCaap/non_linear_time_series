%% a

delta = 2^(-9);
theta = [0.7, 0.8, 3.0, -0.34];

sigmas = [0, 0.1, 0.2, 0.3, 0.4];

T = 100;
t = 0:delta:T;
N = length(t);

Y1 = zeros(length(sigmas), N);
Y2 = zeros(length(sigmas), N);
Y1(1) = -1.9;
Y2(1) = 1.2;

for i=1:length(sigmas)
    sigma = sigmas(i);
    for n=2:N
        Y1(i,n) = Y1(i,n-1) + theta(3) * (Y1(i,n-1) + Y2(i,n-1) - 1/3 * Y1(i,n-1)^3 + theta(4)) * delta + sigma*sqrt(delta)*randn(1);
        Y2(i,n) = Y2(i,n-1) - 1/theta(3) * (Y1(i,n-1) + theta(2)*Y2(i,n-1) - theta(1)) * delta;
    end
end

%% Quiver
[x1, x2] = meshgrid(linspace(min(Y1(:)), max(Y1(:)), 15), linspace(min(Y2(:)), max(Y2(:)), 15));
% Define derivatives
dx1 = theta(3) * (x1 + x2 - (1/3)*x1.^3 + theta(4));
dx2 = -(1/theta(3)) * (x1 + theta(2)*x2 - theta(1));

% Normalize vectors for better visualization
magnitude = sqrt(dx1.^2 + dx2.^2);
dx1n = dx1 ./ magnitude;
dx2n = dx2 ./ magnitude;

%%
figure;
t = tiledlayout(2, 3);
for i=1:length(sigmas)
    nexttile;
    hold on
    plot(Y1(i,:), Y2(i,:))
    % Possibly add quiver plot
    quiver(x1, x2, dx1n, dx2n, 0.6);
    xlabel("Y1")
    ylabel("Y2")
    axis tight;    % removes large margins
    title("\sigma = " + sigmas(i))
    hold off
end
t.TileSpacing = 'compact';
t.Padding = 'compact';

%% b

nbins = 100;

figure;
t = tiledlayout(2, 3);
for i=1:length(sigmas)
    nexttile;
    histogram2(Y1(i,:), Y2(i,:), nbins, 'Normalization','probability', 'FaceColor','flat')
    xlabel("Y1")
    ylabel("Y2")
    title("\sigma = " + sigmas(i))
end
t.TileSpacing = 'compact';
t.Padding = 'compact';


figure;
t = tiledlayout(2, 3);
for i=1:length(sigmas)
    nexttile;
    histogram2(Y1(i,:), Y2(i,:), nbins, 'DisplayStyle','tile','ShowEmptyBins','on')
    xlabel("Y1")
    ylabel("Y2")
    title("\sigma = " + sigmas(i))
end
t.TileSpacing = 'compact';
t.Padding = 'compact';