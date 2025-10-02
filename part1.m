%% SETAR(2;1;1)
N = 100;
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
figure(1)
plot(x)
yline(r1, 'k--')
xlabel("Time")
ylabel("x_t")
title("Simulation of SETAR(2;1;1)")

%% IGAR(2;1) 
N = 100;
e = randn(N,1);
x = zeros(N,1);

% Regime probabilites
p1 = 0.8;
p2 = 1 - p1;

% Model parameters
a0 = [-1.0, 1.0];
a1 = [0.6, 0.4];

% Regimes
u = rand(N,1);
r = (u <= p2) + 1;

% Initialize the process depending on the regime
if r(1) == 1
    x(1) = sqrt(1 / (1 - a1(1)^2)) * e(1);
else
    x(1) = sqrt(1 / (1 - a1(2)^2)) * e(1);
end

% Simulate the process
for n=2:N
    if r(n) == 1
        x(n) = a0(1) + a1(1)*x(n-1) + e(n);
    else
        x(n) = a0(2) + a1(2)*x(n-1) + e(n);
    end
end

figure(2)
subplot(2,1,1)
plot(x)
ylabel("x_t")
title("IGAR(2;1)")

subplot(2,1,2)
stairs(r)
ylim([0.8, 2.2])
xlabel("Time")
ylabel("Regime")

%% MMAR(2;1)
N = 100;
e = randn(N,1);
x = zeros(N,1);

% Transition probabilites
P = [0.9, 0.5;
    0.1, 0.5];

% Predetermine regimes
u = rand(N,1);
r = zeros(N,1);
r(1) = 1;
for n=2:N
    state = r(n-1);
    other_state = max(1, (2-state)*2);
    if u(n) <= P(state, state)
        r(n) = state;
    else 
        r(n) = other_state;
    end
end

% Parameters
a0 = [-1.0, 1.0];
a1 = [0.6, 0.4];

x(1) = sqrt(1 / (1 - a1(1)^2)) * e(1);

% Simulate the process
for n=2:N
    if r(n) == 1
        x(n) = a0(1) + a1(1)*x(n-1) + e(n);
    else
        x(n) = a0(2) + a1(2)*x(n-1) + e(n);
    end
end

figure(3)
subplot(2,1,1)
plot(x)
title("MMAR(2;1)")
ylabel("x_t")

subplot(2,1,2)
stairs(r)
ylim([0.8, 2.2])
xlabel("Time")