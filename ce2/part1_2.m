clear all; close all; clc;

%% --- Simulate data 
N = 1000;
e_k = randn(N,1);
x_k = zeros(N,1);

% True parameters
r_true  = 0;
a10 = -1; a11 = 0.6;
a20 =  1; a21 = 0.4;

x_k(1) = sqrt(1/(1 - a11^2)) * e_k(1);
for n = 2:N
    I1 = (x_k(n-1) <= r_true);
    I2 = (x_k(n-1) >  r_true);
    x_k(n) = (a10 + a11 * x_k(n-1)) * I1 + (a20 + a21 * x_k(n-1)) * I2 + e_k(n);
end

figure;
plot(x_k); yline(r_true, 'k--', 'threshold');
title('Simulated SETAR(2;1;1) series');
xlabel('t'); ylabel('x_t');


%%
d_candidates = 1:5;                            % possible delay values
r_candidates = linspace(min(x_k), max(x_k), 50); % grid of thresholds
r_candidates = r_candidates(10:40);
best_loss = Inf;                               % initialize
best_params = struct();

% Store loss surface for visualization
loss_surface = zeros(length(d_candidates), length(r_candidates));

for i = 1:length(d_candidates)
    d = d_candidates(i);
    for j = 1:length(r_candidates)
        r = r_candidates(j);
        
        % Create indicator vectors
        I1 = double(x_k(1:end-d) <= r);
        I2 = double(x_k(1:end-d) >  r);
        
        % Construct big design matrix
        X = [I1, I1.*x_k(1:end-d), I2, I2.*x_k(1:end-d)];
        y = x_k((d+1):end);
        
        % OLS estimate
        theta_l_hat = (X' * X) \ (X' * y);
        
        % Residuals and variance
        res = y - X * theta_l_hat;
        sigma2_hat = (res' * res) / (length(y) - 4);
        
        % Covariance matrix and SEs
        cov_theta = sigma2_hat * inv(X' * X);
        se_theta = sqrt(diag(cov_theta));
        
        % Objective for given (r, d): minimize in theta_l = [alpha0, alpha1, beta0, beta1]
        Q = @(theta_l) setar_sse_l(theta_l, r, d, x_k);
        fval = Q(theta_l_hat);
        loss_surface(i,j) = fval;   % store for heatmap

        % Objective for given (r, d): minimize in theta_l = [alpha0, alpha1, beta0, beta1]
        %Q = @(theta_l) setar_sse_l(theta_l, r, d, x_k);
        %init_guess = [0, 0.5, 0, 0.5];
        
        %[theta_l_hat, fval] = fminsearch(Q, init_guess);
        %loss_surface(i,j) = fval;   % store for heatmap
        %se_theta = [0, 0, 0, 0];
        
        if fval < best_loss
            best_loss = fval;
            best_params = struct( ...
                'theta_l', theta_l_hat, ...
                'r', r, ...
                'd', d, ...
                'loss', fval, ...
                'se', se_theta ...
            );
        end
    end
end


%% --- Display best results ---
disp('Best parameter estimates:');
fprintf('alpha0 = %.3f (SE=%.3f), alpha1 = %.3f (SE=%.3f)\n', best_params.theta_l(1), best_params.se(1), best_params.theta_l(2), best_params.se(2));
fprintf('beta0  = %.3f (SE=%.3f), beta1  = %.3f (SE=%.3f)\n', best_params.theta_l(3), best_params.se(3), best_params.theta_l(4), best_params.se(4));
fprintf('r_hat  = %.3f, d_hat  = %d\n', best_params.r, best_params.d);


%% --- Plot loss surface ---
figure;
imagesc(r_candidates, d_candidates, loss_surface);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Threshold r');
ylabel('Delay d');
title('Loss surface Q(\theta_l) over (r,d)');
hold on;
plot(best_params.r, best_params.d, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
text(best_params.r, best_params.d, '  Optimum', 'Color', 'r', 'FontWeight', 'bold');
hold off;


%% --- Supporting function ---
function SSE = setar_sse_l(theta_l, r, d, x)
    % Extract parameters
    alpha0 = theta_l(1); alpha1 = theta_l(2);
    beta0  = theta_l(3); beta1  = theta_l(4);
    N = length(x);

    err = zeros(N - d, 1);
    for t = (d + 1):N
        if x(t - d) <= r
            x_hat = alpha0 + alpha1 * x(t - 1);
        else
            x_hat = beta0 + beta1 * x(t - 1);
        end
        err(t - d) = x(t) - x_hat;
    end
    
    SSE = sum(err.^2) / (N - d);   % mean squared error
end