%% SETAR(2;1;1)
N = 10000;
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
% Iterate over multiple values of d
best_val = Inf;
for d = 1:10
    Q = @(p) sum((x(d+1:N) ...
        - (p(1) + p(2).*x(d:N-1)).*(x(1:N-d)<=p(5)) ...
        - (p(3) + p(4).*x(d:N-1)).*(x(1:N-d)>p(5))).^2);
    
    %Q = @(p) sum((x(d+1:N) ...
    %    - (-1 + p(1)*x(d:N-1)).*(x(1:N-d)<=r1) ...
    %    - (1 + p(2)*x(d:N-1)).*(x(1:N-d)>r1)).^2);

    options = optimset('MaxIter', 1000, 'MaxFunEvals', 10000, 'Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6);
    p0 = [-1, 0, 1, 0, 0];
    %p0 = [0,0];
    % Find minimum using fminsearch
    [p_opt, fval] = fminsearch(Q, p0, options);
    if fval < best_val
        best_val = fval;
        best_d = d;
        best_p = p_opt;
        best_r = r;
    end
end

%% 
% Iterate over multiple values of d and r
best_val = Inf;
for d = 1:10
    r_cands = linspace(min(x), max(x), 100);
    for r = r_cands
        Q = @(p) sum((x(d+1:N) ...
            - (p(1) + p(2)*x(d:N-1)).*(x(1:N-d)<=r) ...
            - (p(3) + p(4)*x(d:N-1)).*(x(1:N-d)>r)).^2);
        
        %Q = @(p) sum((x(d+1:N) ...
        %    - (-1 + p(1)*x(d:N-1)).*(x(1:N-d)<=r1) ...
        %    - (1 + p(2)*x(d:N-1)).*(x(1:N-d)>r1)).^2);
    
        options = optimset('MaxIter', 1000, 'MaxFunEvals', 10000, 'Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6);
        p0 = [-1, 0, 1, 0];
        %p0 = [0,0];
        % Find minimum using fminsearch
        [p_opt, fval] = fminsearch(Q, p0, options);
        if fval < best_val
            best_val = fval;
            best_d = d;
            best_p = p_opt;
            best_r = r;
        end
    end
end