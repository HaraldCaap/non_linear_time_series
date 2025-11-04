% EKF algorithm for use in Part 4 of computer exercise 2 in
% Advanced Time Series Analysis
%

% Simulate and store observations in "y"


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

zt=[0 aInit]';      % Initial value for z_t
Pt=[Re 0;0 Ra];     % Initial variance of states
Ht=[1 0];           % Differentiated observation function

z(length(y),2)=0;   % Allocates "z" to store values
avar=1:length(y);   % Allocate "avar" to store the variance estimates of a

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
