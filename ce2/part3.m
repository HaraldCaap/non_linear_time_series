%% Doubly AR(1)-AR(1)
N = 500;

mu = 0.8;
a = 0.85;
delta = mu*(1-a);

sigma_eps = 1;
sigma_zeta = 0.1;

eps = sigma_eps*randn(N,1);
zeta = sigma_zeta*randn(N,1);
y = zeros(N,1);
phi = zeros(N,1);

phi(1) = mu + zeta(1);
y(1) = phi(1) + eps(1);

for t = 2:N
    phi(t) = a * phi(t-1) + delta + zeta(t);
    y(t) = phi(t)*y(t-1) + eps(t);
end

% Plot the Simulation
figName = ['mu' num2str(mu) '_a' num2str(a)];
figure('Name',figName);
%sgtitle('Simulation of Doubly AR(1)-AR(1) Model');

% === Subplot 1: Y_t with highlighted background ===
subplot(2,1,1)
hold on

% Determine regions where |phi| >= 1
idx = abs(phi) >= 1;

% Get start and end points of consecutive |phi| >= 1 segments
d_idx = diff([0; idx; 0]);
startPts = find(d_idx == 1);
endPts = find(d_idx == -1) - 1;

% Shade the regions
yl = [min(y) max(y)];
for i = 1:length(startPts)
    fill([startPts(i) endPts(i) endPts(i) startPts(i)], ...
         [yl(1) yl(1) yl(2) yl(2)], ...
         [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
end

% Plot the Y_t line on top
plot(y, 'b', 'LineWidth', 1.2)
ylabel('Y_t')
xlim([1 N])
%title('Y_t with Highlighted |phi_t| ≥ 1 Regions')
hold off

% === Subplot 2: Phi_t ===
subplot(2,1,2)
plot(phi, 'b', 'LineWidth', 1.2)
hold on
yline(1, 'k--', 'LineWidth', 1.2)
yline(-1, 'k--', 'LineWidth', 1.2)
ylabel('\Phi_t')
ylim([min(min(phi)*1.5, 0), max(max(phi)*1.2,0)])
hold off
