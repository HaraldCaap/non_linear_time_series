% part4
clear all; 
close all; 
% ===== LaTeX interpreters (safe defaults) =====
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% ===== Load & basic variables =====
T   = readtable('DataPart4.csv');   
Phi = T.Ph;
Ti  = T.Ti;
Te  = T.Te;
W   = T.W; % wind speed

dT  = Ti - Te;                      % ΔT = T^i - T^e

figure('Name','Sanity checks','Color','w');
subplot(1,3,1);
scatter(W, Phi, 10, 'filled'); grid on;
xlabel('$W$'); ylabel('$\Phi$'); title('Heat load vs wind $W$');

subplot(1,3,2);
scatter(W, dT, 10, 'filled'); grid on;
xlabel('$W$'); ylabel('$T^i - T^e$'); title('$\Delta T$ vs wind $W$');

subplot(1,3,3);
scatter(dT, Phi, 10, 'filled'); grid on;
xlabel('$T^i - T^e$'); ylabel('$\Phi$'); title('$\Phi$ vs $\Delta T$');

%% 1D approach 1A
mask1A   = isfinite(Phi) & isfinite(dT) & isfinite(W) & (abs(dT) > 0); % all data
Ua_obs1A = Phi(mask1A) ./ dT(mask1A);   % noisy observations of U_a(W)
W_use1A  = W(mask1A);

points1D_1A = 200;      % grid points for 1D smoother
order1D_1A  = 1;        % local linear
h1D_1A      = 0.20;     % 20% nearest-neighbour bandwidth

[xgrid_1A, Ua_hat_1A] = regsmooth1D([W_use1A Ua_obs1A], points1D_1A, order1D_1A, h1D_1A);

figure('Name','Ua(W) with dT color — 1A','Color','w'); 
scatter(W_use1A, Ua_obs1A, 12, dT(mask1A), 'filled', 'MarkerFaceAlpha', 0.35); hold on;
plot(xgrid_1A, Ua_hat_1A, 'k-', 'LineWidth', 2);
grid on; box on;
xlabel('Wind speed $W$', 'Interpreter','latex'); 
ylabel('$U_a(W)$', 'Interpreter','latex');
titleString = sprintf('All $| \\Delta T | > 0$) and h = %.1f', h1D_1A);
title(titleString, 'Interpreter', 'latex');cb = colorbar; 
cb.Label.String = 'Temperature difference $\Delta T$';
cb.Label.Interpreter = 'latex'; % Ensure LaTeX interpreter is usedcolormap turbo;


%% 1D approach different threshold 1B
dT_thr1B = 3;   % chosen |ΔT| threshold
mask1B   = isfinite(Phi) & isfinite(dT) & isfinite(W) & (abs(dT) > dT_thr1B);
Ua_obs1B = Phi(mask1B) ./ dT(mask1B);   
W_use1B  = W(mask1B);

points1D_1B = 200; 
order1D_1B  = 1; 
h1D_1B      = 0.20; 

[xgrid_1B, Ua_hat_1B] = regsmooth1D([W_use1B Ua_obs1B], points1D_1B, order1D_1B, h1D_1B);

figure('Name','Ua(W) with dT color — 1B','Color','w'); 
scatter(W_use1B, Ua_obs1B, 12, dT(mask1B), 'filled', 'MarkerFaceAlpha', 0.35); hold on;
plot(xgrid_1B, Ua_hat_1B, 'k-', 'LineWidth', 2);

grid on; box on;
xlabel('Wind speed $W$', 'Interpreter','latex'); 
ylabel('$U_a(W)$', 'Interpreter','latex');
title(sprintf('$|\\Delta T| > %d$', dT_thr1B), ...
      'Interpreter','latex');
cb = colorbar; 
cb.Label.String = 'Temperature difference $\Delta T$';
cb.Label.Interpreter = 'latex'; 
colormap turbo;


%% 1D approach weighted weights
mask = isfinite(Phi) & isfinite(dT) & isfinite(W) & (abs(dT) > 0.5);   % stricter |ΔT| floor
Ua_obs = Phi(mask) ./ dT(mask);
W_use  = W(mask);
dT_use = dT(mask);

% --- optional robust clipping of extreme Ua_obs ---
iqr_bounds = quantile(Ua_obs, [0.01 0.99]);               % gentle winsorization
Ua_obs = max(min(Ua_obs, iqr_bounds(2)), iqr_bounds(1));

% --- variance-aware weights: w ∝ (ΔT)^2 ---
w = dT_use.^2; 
w = w / median(w);                                       % scale for stability

points1D = 200;
order1D  = 1;
h1D      = 0.20;

% If regsmooth1D supports weights: regsmooth1D([x y w], ...)
% (If not, consider a weighted alternative or approximate by jittered replication.)
[xgrid_1D, Ua_hat_1D] = regsmooth1D([W_use Ua_obs w], points1D, order1D, h1D);

figure('Name','Ua(W) with dT color','Color','w'); 
scatter(W_use, Ua_obs, 12, dT_use, 'filled', 'MarkerFaceAlpha', 0.35); hold on;
plot(xgrid_1D, Ua_hat_1D, 'k-', 'LineWidth', 2);
grid on; box on;
xlabel('Wind speed $W$', 'Interpreter','latex'); 
ylabel('$U_a(W)$', 'Interpreter','latex');
title('Color-coded by $\Delta T$', 'Interpreter','latex');
cb = colorbar; 
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\Delta T$';
colormap turbo


%%
% ===== 2D conditional-parametric approach via regsmooth2D =====
% Model: Phi ≈ a(W) + b(W) * dT  (linear in dT, coefficients vary smoothly with W)
% We recover b(W) = U_a(W) by finite difference in the fitted surface along dT.
data   = [W, dT, Phi];              % [x, y, z] = [W, ΔT, Φ]

points2D = 150;                      % grid resolution (W × ΔT)
order2D  = 1;                        % local linear
h2D      = 0.20;                     % bandwidth (tune)
bound    = [min(W) max(W) min(dT) max(dT)];

% conditional_on_y = 1  -> weight only in x=W; linear in y=ΔT
[xGrid, dtGrid, phiFit, ~] = regsmooth2D(data, points2D, order2D, h2D, bound, 1);

% xGrid, yGrid are meshgrids (size points2D × points2D)
% zFit(row, col) corresponds to yGrid(row, col) and xGrid(row, col)

% ---- Extract slope wrt ΔT at each W (b(W) = U_a(W)) ----
% Since z ≈ a(W) + b(W)*ΔT (linear in ΔT), slope is constant in ΔT at fixed W.
% b(W)=E[∂ΔT∂Φ|​​W],
% Take two ΔT levels away from the edges to reduce boundary bias:
iy1 = max(2, floor(0.30 * points2D));
iy2 = min(points2D-1, ceil(0.70 * points2D));

% vectors for plotting
W_grid_vec  = xGrid(1, :);           % W grid along columns
dT1 = dtGrid(iy1, 1);
dT2 = dtGrid(iy2, 1);

% fitted Φ at these two ΔT levels (row slices)
Phi_y1 = phiFit(iy1, :);               % 1 × points2D (over W)
Phi_y2 = phiFit(iy2, :);               % 1 × points2D

Ua_hat_2D = (Phi_y2 - Phi_y1) ./ (dT2 - dT1);  % slope wrt ΔT = U_a(W)

% ===== Compare the estimates: 1A, 1B, and 2 =====
figure('Name','Comparison: 1A, 1B, and 2','Color','w');
hold on; grid on; box on;

% 2D conditional–parametric
plot(W_grid_vec, Ua_hat_2D, 'k-', 'LineWidth', 2, ...
     'DisplayName','2D cond-param $\hat{U}_a(W)$');

% 1D smooth with all |dT| > 0
plot(xgrid_1A, Ua_hat_1A, 'r--', 'LineWidth', 2, ...
     'DisplayName','1D smooth (1A)');

% 1D smooth with |dT| > dT_thr1B
plot(xgrid_1B, Ua_hat_1B, 'b-.', 'LineWidth', 2, ...
     'DisplayName', sprintf('1D smooth (1B, $|\\Delta T| > %d$)', dT_thr1B));

xlabel('Wind speed $W$','Interpreter','latex'); 
ylabel('$U_a(W)$','Interpreter','latex');
title('Estimated $U_a(W)$','Interpreter','latex');
legend('Location','best','Interpreter','latex');


figure('Name','Fitted surface: Φ(W,ΔT)','Color','w');
surf(xGrid, dtGrid, phiFit, 'EdgeColor','none'); colormap parula; view(135,30);
xlabel('$W$'); ylabel('$\Delta T$'); zlabel('$\hat{\Phi}$'); grid on;
title('Conditional-parametric fit: $\hat{\Phi}(W,\Delta T)$');


%% === Phi(W,0) = a(W) from the fitted surface (correct) ===
% Uses b(W) already computed as Ua_hat_2D and two interior ΔT slices (iy1, iy2).

% Intercept using each slice, then average (they should be near-identical if linear):
a_from_y1 = Phi_y1 - Ua_hat_2D .* dT1;   % a(W) via slice at ΔT = dT1
a_from_y2 = Phi_y2 - Ua_hat_2D .* dT2;   % a(W) via slice at ΔT = dT2
a_hat_2D  = 0.5 * (a_from_y1 + a_from_y2);

% Plot Φ(W,0) = a(W)
figure('Name','Phi(W,0) = a(W)','Color','w');
plot(W_grid_vec, a_hat_2D, 'LineWidth', 2);
grid on;
xlabel('Wind speed $W$','Interpreter','latex');
ylabel('$\hat{\Phi}(W,0)$','Interpreter','latex');
title('Intercept $a(W)$ from conditional-parametric fit','Interpreter','latex');
