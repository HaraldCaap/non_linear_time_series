%% a
clear all; close all;
%% Parameters
N = 500;
M = 20; % Number of simulations

a = 0.4;
sigma_v = 1;
sigma_e = 1;

%% Simulate
v = sigma_v * randn(N,M);
e = sigma_e * randn(N,M);

X = zeros(N,M);
Y = zeros(N,M);

for k=1:M
    X(1,k) = a*1/(1-a^2)*sigma_v+v(1,k);
    Y(1,k) = X(1,k) + e(1,k);

    for t=2:N
        X(t,k) = a*X(t-1,k) + v(t,k);
        Y(t,k) = X(t,k) + e(t,k);
    end
end

%%
%----------------------------
% Here comes the EKF algorithm
%----------------------------
% Initializing:
aInit=0.5;  % Initial value of "a"
aVarInit=1; % Initial variance on "a" try to var
sigma_v=1;  % sigma_v, try 10 and 1

Rv=[sigma_v 0; 0 0];  
Re=1;               % sigma_e
Ra=aVarInit;

X_hat = zeros(N,M);
A = zeros(N,M);
Avar = zeros(N,M);

for k=1:M 
    y = Y(:,k);

    zt=[0 aInit]';      % Initial value for z_t
    Pt=[Re 0;0 Ra];     % Initial variance of states
    Ht=[1 0];           % Differentiated observation function
    
    z = zeros(length(y),2);   % Allocates "z" to store values
    z(1,:) = zt;
    avar=zeros(length(y),1);   % Allocate "avar" to store the variance estimates of a
    avar(1) = Pt(2,2);
    
    % The actual estimation
    for i=1:length(y)-1
      ft=[zt(2)*zt(1); zt(2)];         % Transition function
      Ft=[zt(2) zt(1); 0 1];           % Differentiated transition function, Eq. 7.49
    
      Kt=Ft*Pt*Ht'* (Ht*Pt*Ht'+Re)^-1; % Kalman gain, Eq. 7.46
    
      zt=ft + Kt*(y(i)-zt(1));         % Predicted state estimate, Eq. 7.45
      Pt=Ft*Pt*Ft' + Rv - Kt*(Re+Ht*Pt*Ht')*Kt'; % Predicted state covariance estimate, Eq. 7.47 - 7.48
    
      % Write out
      z(i+1,:)=zt'; % Store the values of z_t in order to plot afterwards
      avar(i+1)=Pt(2,2);
    end
    X_hat(:,k) = z(:,1);  % Store the estimated state
    A(:,k) = z(:,2);      % Store the estimated parameter
    Avar(:,k) = avar;     % Store the variance estimates
end

%% Plots
figure; grid on; hold on;
sgtitle("a_{initial}="+aInit + ", \sigma_v^2="+ sigma_v + ", initial Var[a]="+aVarInit)
yline(a, 'k--', 'LineWidth', 2)
plot(A, 'b:')
plot(mean(A,2), 'b', 'LineWidth', 2)
ylabel("a_t")
hold off


%%
figure;
subplot(211)
plot(Y(:,1))
title("Y_t")
subplot(212)
plot(X(:,1))
title("X_t")

figure;
plot(x, 'LineWidth',1.2); hold on;
plot(y, 'LineWidth',0.8);
legend("True state x", "Observation y");
title('True state and observations'); grid on;

figure;
plot(x, 'LineWidth',1.2); hold on;
plot(z(:,1), 'LineWidth',0.8);
legend("True state x", "Predicted state");
title('Prediction vs truth'); grid on;

figure;
subplot(3,1,1);
plot(z(:,2));
hold on
yline(a, 'k--')
hold off
title('Parameter estimation')

subplot(3,1,2);
plot(x - z(:,1));
title('Prediction error (x - predicted x)'); grid on;

subplot(3,1,3);
plot(avar);
title('Predicted variance'); grid on;

figure;
plot(z(:,2));
hold on
plot(z(:,2)+2*avar, 'r--')
plot(z(:,2)-2*avar, 'r--')
yline(a, 'k--')
hold off
title('Parameter estimation with confidence interval')

%figure;
%histogram(x - z(:,1), 30);
%title('Histogram of prediction errors'); grid on;