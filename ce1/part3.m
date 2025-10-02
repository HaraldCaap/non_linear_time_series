close all
clear all

%% Simulation
r1 = 0;

N = 1000;
e = randn(N,1);
X = zeros(N,1);

a0 = [1.0, -1.0];
a1 = [0.6, 0.4];

if 0 <= r1
    X(1) = sqrt(1 / (1 - a1(1)^2)) * e(1);
else
    X(1) = sqrt(1 / (1 - a1(2)^2)) * e(1);
end

for n=2:N
    if X(n-1) <= r1
        X(n) = a0(1) + a1(1)*X(n-1) + e(n);
    else
        X(n) = a0(2) + a1(2)*X(n-1) + e(n);
    end
end

figure(1)
plot(X)
title("SETAR(2;1;1)")


%% 

% Hist regression prm.
bin = [-2 2];
n_bin = 20;
h = (bin(2) - bin(1))/n_bin;
bin_points = (bin(1)+h/2):h:(bin(2)-h/2);

% Hist regression
% the starting bin
cur_bin=[bin_points(1)-0.5*h bin_points(1)+0.5*h];
% init
lambda=zeros(n_bin,1);
gamma=zeros(n_bin,1);
f_hat=zeros(n_bin,1);

for i=1:n_bin
    index=(X(1:end-1)>cur_bin(1) & X(1:end-1)<=cur_bin(2));
    if (sum(index)>5)
        lambda(i) = sum( X(2:end).*index ) / sum(index);
        % f_hat is needed for confidence band calculation
        f_hat(i) = (n_bin*h)^(-1) * sum(index);
        gamma(i) = sum( (X(2:end) - lambda(i) ).^2 .* index ) / sum(index);
    else
        fprintf('Need more points')
        break
    end

    % move to next bin
    cur_bin=cur_bin+h;
end

% Make confidence bands for the cumulated function. Def. (3.10).
% 95% confidence band, c is found in table 3.1
c = 1.273;

Lambda = cumsum(lambda*h);
h_hat = zeros(n_bin,1);
for i=1:n_bin
    h_hat(i) = gamma(i)/f_hat(i);
end
H_hat = cumsum(h_hat*h);

H_hat_b = H_hat(n_bin);
Lambda_lower = Lambda - c .* n_bin.^(-0.5) .* H_hat_b.^(0.5) .* (1 + H_hat/H_hat_b);
Lambda_upper = Lambda + c .* n_bin.^(-0.5) .* H_hat_b.^(0.5) .* (1 + H_hat/H_hat_b);


x = linspace(bin(1), bin(2), 100);
lambda_theo = zeros(length(x),1);
for i=1:length(x)
    if x(i) <= r1
        lambda_theo(i) = a0(1) + a1(1)*x(i);
    else
        lambda_theo(i) = a0(2) + a1(2)*x(i);
    end
end

Lambda_theo = a0(1)*(min(x,r1) - bin(1)) + a1(1)/2 * (min(x,r1).^2-bin(1)^2) + ...
                + a0(2)*(max(x,r1) - r1) + a1(2)/2 * (max(x,r1).^2-r1^2);
Lambda_theo = Lambda_theo(:);
%Mhat_theo = @(x) (a0(1) + a1(1)*x) .* (x <= r1) + (a0(2) + a1(2)*x) .* (x > r1);
%Lambda_theo2 = cumtrapz(x, Mhat_theo(x));


figure(2)
hold on
plot(x, lambda_theo)
plot(bin_points, lambda)
%stairs([bin_points-0.5*h, bin_points(end)+0.5*h], [lambda; lambda(end)])
xline(0, 'k')
yline(0, 'k')
xlabel("X_k")
ylabel("X_{k+1}")
title("Conditional mean")
grid on
hold off

figure(3)
hold on
plot(x, Lambda_theo, 'k--')
plot(bin_points+0.5*h, Lambda, 'b')
plot(bin_points+0.5*h, Lambda_upper, 'r--')
plot(bin_points+0.5*h, Lambda_lower, 'r--')
title("Cumulative conditional mean")
xline(0, 'k')
yline(0, 'k')
grid on
hold off
