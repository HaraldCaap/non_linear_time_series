dataTable = readtable("DataPart5.csv");
X = dataTable.x;
X = X(:);

figure(1)
plot(X)
title('Time Series Data');
xlabel('Time');
ylabel('Value');

figure(2)
subplot(2,1,1);
autocorr(X);
title('ACF of the Data');

subplot(2,1,2);
parcorr(X);
title('PACF of the Data');
%% Remove mean
X = X - mean(X);

%% Fit ARMA model and check residuals
p = 2;
q = 0;

model = arima(p, 0, q);  % 0 for no differencing
fitModel = estimate(model, X);

residuals = infer(fitModel, X);

figure(3)
subplot(3,1,1);
plot(residuals);
title('Residuals');

subplot(3,1,2);
autocorr(residuals);
title('ACF of Residuals');

subplot(3,1,3);
parcorr(residuals);
title('PACF of Residuals');

%% Check LDF for non-linearity

% ldf
order = 1;
points = 100;
h = 0.6;
maxlag = 20;
[phi]=ldf(residuals, order, points, h, maxlag);

figure(4)
hold on
stem(0:maxlag, phi, 'k', 'filled')
yline(0, 'k')
grid on
hold off

%% Plot X_t against X_{t-2}

figure(5)
scatter(X(1:end-2), X(3:end));
xline(-4.12)
xline(-0.23)
xline(2.77)
ylabel("X_t")
xlabel("X_{t-2}")
grid on

figure(6)
scatter(residuals(1:end-2), residuals(3:end))
ylabel("e_t")
xlabel("e_{t-2}")
grid on

ldfone(residuals, 1, 100, h, 2);

%% Fit a SETAR

tt = [-4.1, -0.3, 2.8];
mdl1 = arima(2,0,0);
mdl2 = arima(2,0,0);
mdl3 = arima(2,0,0);
mdl4 = arima(2,0,0);

reg1 = X <= tt(1);
reg2 = tt(1) < X & X <= tt(2);
reg3 = tt(2) < X & X <= tt(3);
reg4 = tt(3) < X;

X1 = X(min(reg1+2, length(X)));
X2 = X(min(reg2+2, length(X)));
X3 = X(min(reg3+2, length(X)));
X4 = X(min(reg4+2, length(X)));

fitmdl1 = estimate(mdl1, X1);
fitmdl2 = estimate(mdl2, X2);
fitmdl3 = estimate(mdl3, X3);
fitmdl4 = estimate(mdl4, X4);

residuals_setar = [infer(fitmdl1, X1), 
    infer(fitmdl2, X2), 
    infer(fitmdl3, X3), 
    infer(fitmdl4, X4)];

figure(6)
subplot(3,1,1);
plot(residuals_setar);
title('Residuals');

subplot(3,1,2);
autocorr(residuals_setar);
title('ACF of Residuals');

subplot(3,1,3);
parcorr(residuals_setar);
title('PACF of Residuals');

