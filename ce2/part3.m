%% Doubly AR(1)-AR(1)
N = 500;

mu = 0.8;
a = 0.85;
delta = mu*(1-a);
sigma_eps = 4;
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

% Plot the generated time series
figure;
sgtitle('Generated Doubly AR(1)-AR(1) Time Series');
subplot(2,1,1)
plot(y);
ylabel('Y_t');
subplot(2,1,2)
plot(phi)
ylabel('Phi_t')
