%% Part 5 — Next steps: ARMA, LDF on residuals (using ldf.m & ldfone.m)
% This script continues the analysis using the provided LDF utilities.
% It fits an ARMA model, computes residuals, applies LDF to the residuals,
% and visualizes a chosen lag's dependence with ldfone.m.

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
close all;
clear; 

%% 1) Load data (absolute path for robusbest_lag=2;

dataPath = 'DataPart5.csv';
T = readtable(dataPath);
vars = T.Properties.VariableNames;
icol = find(varfun(@isnumeric, T, 'OutputFormat','uniform'), 1, 'first');
if isempty(icol)
    error('No numeric column found in %s', dataPath);
end
y = T.(vars{icol});
y = y(:);
N = numel(y);

figure('Name','Raw series','Color','w');
plot(y,'b-'); grid on; xlabel('t'); ylabel('$y_t$');
title('Raw series');

%% 2) Stationarity quick check and demean
isStat = true;
if exist('adftest','file') == 2
    try
        isStat = adftest(y);
    catch
        isStat = true;
    end
end
if ~isStat
    y = diff(y);
    fprintf('Series differenced once to achieve stationarity.\n');
end
y = y - mean(y,'omitnan');

%% 3) ACF & PACF (visual guidance)

figure('Name','ACF/PACF of y','Color','w');
subplot(3,1,1); plot(y,'b-'); grid on; xlabel('t'); ylabel('$y_t$'); title("Y with zero mean")
subplot(3,1,2); autocorr(y, 'NumLags', 30); title('ACF of $y_t$');
subplot(3,1,3); parcorr(y, 'NumLags', 30); title('PACF of $y_t$');

%% 4) Fit ARMA(p,q) (user can tune p,q). Compute residuals
p = 2; q = 0; d = 0;  % Based on earlier selection; change if desired
Mdl = arima(p,d,q);
EstMdl = estimate(Mdl, y, 'Display','off');
[res,~,~] = infer(EstMdl, y);
res = res(:);
res = res(~isnan(res));
rng default

% Preprocess residuals (without touching teacher code): clip extremes and jitter slightly
res_use = res;
lo = prctile(res_use,1); hi = prctile(res_use,99);
res_use = min(max(res_use, lo), hi);
res_use = res_use + 1e-8*randn(size(res_use));

fprintf('\nEstimated ARMA(%d,%d):\n', p,q);
disp(EstMdl);

figure('Name',sprintf('ARMA(%d,%d) residuals',p,q),'Color','w');

subplot(3,2,1); 
plot(res,'b-'); grid on; title('Residuals','Interpreter','latex'); xlabel('t'); ylabel('$e_t$','Interpreter','latex');

subplot(3,2,2); 
histogram(res,'Normalization','pdf'); grid on; title('Residual histogram','Interpreter','latex');

subplot(3,2,3); 
autocorr(res,'NumLags',30); title('ACF of $e_t$','Interpreter','latex');

subplot(3,2,4); 
parcorr(res,'NumLags',30); title('PACF of $e_t$','Interpreter','latex');

subplot(3,2,5); 
autocorr(res.^2,'NumLags',30); title('ACF of $e_t^2$','Interpreter','latex');

subplot(3,2,6); 
parcorr(res.^2,'NumLags',30); title('PACF of $e_t^2$','Interpreter','latex');

%% 6) LDF on residuals (short + robust)
order  = 1; points = 101; h = 0.30; maxlag = min(20, numel(res_use)-1);
phi = [];
try, phi = ldf(res_use, order, points, h, maxlag); end
if isempty(phi)
    % simple histogram LDF fallback
    nBins = 20; maxlag = min(maxlag, max(10, floor(numel(res_use)/5)));
    phi = [1; arrayfun(@(k) hist_ldf(res_use,k,nBins), (1:maxlag)')];
    order = 0; h = NaN; points = nBins;
end

figure('Name','LDF on residuals','Color','w');
stem(0:maxlag, phi, 'filled'); grid on; yline(0,'k-');
xlabel('Lag k'); ylabel('$\phi(k)$');
title(sprintf('LDF on residuals (order=%d, h=%.2f)', order, h));

% Pick the lag with the largest absolute dependence (excluding k=0)
phiTail = phi(2:end);
phiTail(~isfinite(phiTail)) = 0;  % treat NaNs/Infs as 0 for selection
[~, idxMax] = max(abs(phiTail));
bestLag = idxMax;  % phi includes phi(0) at index 1
if isempty(bestLag) || bestLag < 1 || bestLag > maxlag || ~isfinite(phi(bestLag+1))
    warning('Could not determine best lag from LDF; defaulting to k=1');
    bestLag = 1;
end
fprintf('Max |LDF| at lag k=%d, value=%.3f\n', bestLag, phi(bestLag+1));

%% 7) Visualize dependence at a selectable lag
% Set lagToView to any k to inspect; leave [] to use bestLag
lagToView = [2];
if isempty(lagToView), lagToView = bestLag; end

ldfone(res, order, points, max(h,0.35), lagToView); % quiet plot

figure('Name','Dependence','Color','w');
xlag = res_use(1:end-lagToView); yt = res_use(1+lagToView:end);
edges = linspace(min(xlag), max(xlag), points+1);
centers = 0.5*(edges(1:end-1) + edges(2:end));
mHat = arrayfun(@(b) mean(yt(xlag>=edges(b) & xlag<edges(b+1))), 1:points)';
mHat = fillmissing(mHat,'nearest');
scatter(xlag, yt, 8, 'filled','MarkerFaceAlpha',0.08); hold on; grid on;
plot(centers, mHat, 'r-', 'LineWidth', 2); title(sprintf('k=%d: bin means', lagToView));
xlabel(sprintf('e_{t-%d}', lagToView)); ylabel('e_t');

%% 9) Nonlinear mean at k=2:SETAR
% Build lagged series
Y  = y(:);
Y1 = lagmatrix(Y,1);
Y2 = lagmatrix(Y,2);

% Common index (valid observations)
idx = isfinite(Y) & isfinite(Y1) & isfinite(Y2);
Yc  = Y(idx);
X1  = Y1(idx);
X2  = Y2(idx);

cand1 = -6:0.5:-2;                 % around -4
cand2 = -0.9:0.1:0.5;              % around -0.213
cand3 = 2.0:0.1:3.5;               % around 2.75
best = struct('BIC',Inf); minSeg = 25;

for t1 = cand1
  for t2 = cand2
    if t2 <= t1, continue; end
    for t3 = cand3
      if t3 <= t2, continue; end
      tau_try = [-Inf, t1, t2, t3, Inf];
      % fit OLS per region
      logL = 0; ok = true; kReg = 0; resAll = nan(size(Yc));
      for r = 1:4
        inR = (X2 > tau_try(r)) & (X2 <= tau_try(r+1));
        if nnz(inR) < minSeg, ok = false; break; end
        Xmat = [ones(nnz(inR),1) X1(inR) X2(inR)];
        Yr   = Yc(inR);
        br = Xmat\Yr; er = Yr - Xmat*br;
        sig2 = max(var(er,1), eps);
        logL = logL - 0.5*( nnz(inR)*(log(2*pi*sig2)+1) );
        kReg = kReg + 3 + 1;        % (c,phi1,phi2) + variance
        resAll(inR) = er;
      end
      if ~ok, continue; end
      [~,BIC] = aicbic(logL, kReg, numel(Yc));
      if BIC < best.BIC
        best = struct('BIC',BIC,'tau',tau_try,'logL',logL);
      end
    end
  end
end

% Use tuned breaks if found
if isfield(best,'tau')
  tau_vec = best.tau;
  fprintf('Tuned breaks: (%.3g, %.3g, %.3g)  BIC=%.2f\n', tau_vec(2),tau_vec(3),tau_vec(4), best.BIC);
else
  tau_vec = [-Inf, -4, -0.213, 2.75, Inf];
end
% Baseline AR(2) AIC/BIC via arima
mdl_base = arima(2,0,0);
try
    [est_base,~,logL_base] = estimate(mdl_base, Yc, 'Display','off');
catch
    est_base = mdl_base; logL_base = NaN;
end
% parameter count: const, var, 2 AR
k_base = 1 + 1 + 2; % constant + variance + ARs
[AIC_base,BIC_base] = aicbic(logL_base,k_base,numel(Yc));

% SETAR with 4 regions at k=2: [-inf,-4], [-4,-0.213], [-0.213,2.75], [2.75, inf]
%tau_vec = [-Inf, -4, -0.213, 2.75, Inf];
K = numel(tau_vec)-1;  % 4 regions
betas = zeros(K,3); ns = zeros(K,1); sig2 = zeros(K,1);
logL_setar = 0;
resSETAR = nan(size(Yc));
for r = 1:K
    inR = (X2 > tau_vec(r)) & (X2 <= tau_vec(r+1));
    Xmat = [ones(sum(inR),1) X1(inR) X2(inR)]; Yr = Yc(inR);
    if sum(inR) >= 5
        br = Xmat\Yr; er = Yr - Xmat*br;
        betas(r,:) = br(:)'; ns(r) = numel(Yr);
        sig2(r) = max(var(er,1), eps);
        logL_setar = logL_setar - 0.5*( ns(r)*(log(2*pi*sig2(r))+1) );
        resSETAR(inR) = er;
    else
        betas(r,:) = [NaN NaN NaN]; ns(r)=numel(Yr); sig2(r)=NaN;
    end
end
% params: (c,phi1,phi2) per regime + K variances
k_setar = 3*K + K;
[AIC_setar,BIC_setar] = aicbic(logL_setar,k_setar,numel(Yc));

fprintf('\nSETAR (4 regions at k=2) vs AR(2):\n');
fprintf('  AR(2):      AIC=%.2f  BIC=%.2f\n', AIC_base, BIC_base);
fprintf('  SETAR4:     AIC=%.2f  BIC=%.2f  (breaks = [%g  %g  %g])\n', AIC_setar, BIC_setar, tau_vec(2), tau_vec(3), tau_vec(4));

% Print regime coefficients
for r=1:K
    fprintf('  Reg %d (%.3g, %.3g]: c=%.4f, phi1=%.4f, phi2=%.4f  (n=%d)\n', r, tau_vec(r), tau_vec(r+1), betas(r,1), betas(r,2), betas(r,3), ns(r));
end

% Pretty print the piecewise SETAR(4) mean (k=2)
fprintf('\nSETAR(4) mean (k=2) as piecewise model:\n');
for r=1:K
    % interval text
    a = tau_vec(r); b = tau_vec(r+1);
    if isfinite(a), la = sprintf('%.3g',a); else, la = '-inf'; end
    if isfinite(b), rb = sprintf('%.3g',b); else, rb = '+inf'; end
    fprintf('  if x_{t-2} in (%s, %s]:  x_t = %.4f + %.4f x_{t-1} + %.4f x_{t-2}\n', ...
        la, rb, betas(r,1), betas(r,2), betas(r,3));
end

% Plot fitted lines per regime
figure('Name','SETAR(4) regimes (k=2)','Color','w');
scatter(X2, Yc, 10, 'filled','MarkerFaceAlpha',0.08); hold on; grid on;
x2grid = linspace(min(X2), max(X2), 400)';
x1med  = median(X1,'omitnan');
colors = lines(K);
for r=1:K
    mask = x2grid>tau_vec(r) & x2grid<=tau_vec(r+1);
    if any(mask) && all(isfinite(betas(r,:)))
        yline_r = betas(r,1) + betas(r,2)*x1med + betas(r,3)*x2grid(mask);
        plot(x2grid(mask), yline_r, '-', 'Color',colors(r,:), 'LineWidth', 2);
    end
end
for r=2:numel(tau_vec)-1
    xline(tau_vec(r),'k--');
end
xlabel('x_{t-2}'); ylabel('x_t'); title('SETAR(4) mean by regime');

% Residual diagnostics
figure('Name','SETAR(4) residual diagnostics','Color','w');
subplot(2,2,1); plot(resSETAR,'b-'); grid on; title('Residuals');
subplot(2,2,2); autocorr(resSETAR,'NumLags',30); title('ACF residuals');
subplot(2,2,3); autocorr(resSETAR.^2,'NumLags',30); title('ACF squared residuals');
subplot(2,2,4); histogram(resSETAR, 'Normalization','pdf'); grid on; title('Residual histogram');