a11_grid = linspace(0, 1, 50);
a21_grid = linspace(0, 1, 50);

% Simulate process
rng(1)
N = 3000;
e_k = randn(N,1);
x_k = zeros(N,1);

% True parameters
a10 = -1; a11_true = 0.6;
a20 =  1; a21_true = 0.4;
r = 0;

x_k(1) = sqrt(1/(1 - a11_true^2)) * e_k(1);
for n = 2:N
    I1 = (x_k(n-1) <= r);
    I2 = (x_k(n-1) >  r);
    x_k(n) = (a10 + a11_true * x_k(n-1)) * I1 + (a20 + a21_true * x_k(n-1)) * I2 + e_k(n);
end

% Define index ranges for analysis
ranges = {1:3000, 1:300, 1:30, 1001:3000, 1001:1300, 1001:1030};

% Precompute global Zmin/Zmax
Z_all = cell(length(ranges),1);
Zmin = inf; Zmax = -inf;

for k = 1:length(ranges)
    idx = ranges{k};
    x_subset = x_k(idx);
    Z = zeros(length(a11_grid), length(a21_grid));
    for i = 1:length(a11_grid)
        for j = 1:length(a21_grid)
            theta = [a11_grid(i), a21_grid(j)];
            Z(i,j) = setar_sse_slope(theta, x_subset);
        end
    end
    Z = Z';
    Z_all{k} = Z;
    Zmin = min(Zmin, min(Z(:)));
    Zmax = max(Zmax, max(Z(:)));
end

levels = linspace(Zmin, Zmax, 15);

figure;
sgtitle('Q_N surfaces for different data ranges  —  \color{red}X\color{black} = True parameters,  ● = Estimated minimum', ...
        'FontWeight', 'bold', 'Interpreter', 'tex');

curvatures = zeros(length(ranges),1);

for k = 1:length(ranges)
    idx = ranges{k};
    Z = Z_all{k};

    subplot(2,3,k)
    imagesc(a11_grid, a21_grid, Z, [Zmin Zmax]);
    axis xy
    colorbar
    hold on
    contour(a11_grid, a21_grid, Z, levels, 'k', 'LineWidth', 0.8);

    % Plot true parameters
    plot(a11_true, a21_true, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    % Find and plot minimum
    [~, lin_idx] = min(Z(:));
    [row_idx, col_idx] = ind2sub(size(Z), lin_idx);
    phi11_min = a11_grid(col_idx);
    phi21_min = a21_grid(row_idx);
    plot(phi11_min, phi21_min, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'LineWidth', 1);

    % Numerical curvature measure (approx Hessian)
    [dx, dy] = gradient(Z, a11_grid(2)-a11_grid(1), a21_grid(2)-a21_grid(1));
    [dxx, dxy] = gradient(dx);
    [dyx, dyy] = gradient(dy);
    H = [dxx(row_idx,col_idx), dxy(row_idx,col_idx); dyx(row_idx,col_idx), dyy(row_idx,col_idx)];
    curvatures(k) = min(eig(H)); % smallest eigenvalue = flatness measure

    title(sprintf('Range %d:%d', idx(1), idx(end)));
    xlabel('\alpha_{1}');
    ylabel('\beta_{1}');
end

disp(table((1:length(ranges))', curvatures, 'VariableNames', {'RangeIndex','MinEigenvalue'}))

function sse = setar_sse_slope(p, x)
% SSE for SETAR(2;1;1) estimating slope parameters
phi11 = p(1);
phi21 = p(2);
phi10 = -1; 
phi20 = 1; 
r = 0;

N = length(x);
residuals = zeros(N-1,1);
for t = 2:N
    I1 = (x(t-1) <= r);
    I2 = (x(t-1) >  r);
    x_hat = (phi10 + phi11*x(t-1))*I1 + (phi20 + phi21*x(t-1))*I2;
    residuals(t-1) = x(t) - x_hat;
end
sse = 1/(N-1)*sum(residuals.^2);
end

%%
% ============================================================
% Plot the data and regime switches for 1:30 and 1001:1030
% ============================================================

ranges_regime = {1:300, 1001:1300};
titles = {'Range 1:30', 'Range 1001:1030'};

figure;
for k = 1:length(ranges_regime)
    idx = ranges_regime{k};
    x_subset = x_k(idx);
    x_prev = x_k(idx(1:end-1));

    % Regime indicators (for x_{t-1})
    I1 = (x_prev <= r);
    I2 = (x_prev >  r);

    % Time values and corresponding x_t (for t = 2,...)
    tvals = idx(2:end);
    x_vals = x_subset(2:end);

    subplot(2,1,k)
    hold on

    % Plot regime 1 points and connecting line
    plot(tvals, x_vals, 'o-', 'Color', [0.2 0.6 1]);

    % Threshold line
    yline(0, '--k', 'LineWidth', 1, 'DisplayName', 'Threshold r = 0');

    xlabel('Time index t');
    ylabel('x_t');
    title(['Simulated data and regime switches – ' titles{k}]);
    legend('show', 'Location', 'best');
    grid on
    box on
end