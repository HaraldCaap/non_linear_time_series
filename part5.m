%% Part 5 — Next steps: ARMA, LDF on residuals (using ldf.m & ldfone.m)
% This script continues the analysis using the provided LDF utilities.
% It fits an ARMA model, computes residuals, applies LDF to the residuals,
% and visualizes a chosen lag's dependence with ldfone.m.

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
close all;
clear; 

%% Fig export 
figDir = fullfile(pwd, 'Figs');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

%% Load data
dataPath = '/Users/hcaap/Desktop/Skola/Olinjär Tidsserie/Computer Exercises/comp_ex_1_scrips_2018/matlab/DataPart5.csv';
T = readtable(dataPath);
y = T{:, varfun(@isnumeric, T, 'OutputFormat', 'uniform')}; % Extract first numeric column
y = y(:); % Ensure y is a column vector

fig = figure('Name', 'Raw series', 'Color', 'w');
plot(y, 'b-'); 
grid on; 
xlabel('t'); 
ylabel('$y_t$');
title('Raw series');
saveFigure(fig, figDir, 'P5_raw_series.png');

%% Demean
y = y - mean(y,'omitnan');

%% ACF & PACF (visual guidance)
fig = figure('Name','ACF/PACF of y','Color','w');
subplot(3,1,1); plot(y,'b-'); grid on; xlabel('t'); ylabel('$y_t$'); title("Y with zero mean")
subplot(3,1,2); autocorr(y, 'NumLags', 30); title('ACF of $y_t$');
subplot(3,1,3); parcorr(y, 'NumLags', 30); title('PACF of $y_t$');
saveFigure(fig, figDir, 'P5_acf_pacf_y.png');

%% Fit ARMA(p,q) 
p = 2; q = 0; d = 0;  
Mdl = arima(p,d,q);
EstMdl = estimate(Mdl, y, 'Display','off');

[res, V, ~] = infer(EstMdl, y);          % residuals + variance
z = res ./ sqrt(V);                       % standardized residuals

% --- Plots on standardized residuals ---
fig = figure('Name',sprintf('ARMA(%d,%d) standardized residuals',p,q),'Color','w');

subplot(3,2,1); 
plot(z,'b-'); grid on;
title('Standardized residuals','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$z_t$','Interpreter','latex');

subplot(3,2,2); 
histogram(z,'Normalization','pdf'); grid on; hold on
xg = linspace(min(z),max(z),400);
plot(xg, normpdf(xg,0,1), 'k-', 'LineWidth',1.5);
title('Histogram of $z_t$ with $N(0,1)$','Interpreter','latex');
xlabel('$z_t$','Interpreter','latex');

subplot(3,2,3); 
autocorr(z,'NumLags',30);
title('ACF of $z_t$','Interpreter','latex');

subplot(3,2,4); 
parcorr(z,'NumLags',30);
title('PACF of $z_t$','Interpreter','latex');

subplot(3,2,5); 
autocorr(z.^2,'NumLags',30);
title('ACF of $z_t^2$','Interpreter','latex');

subplot(3,2,6); 
parcorr(z.^2,'NumLags',30);
title('PACF of $z_t^2$','Interpreter','latex');

saveFigure(fig, figDir, sprintf('P5_arma_stdres_p%d_q%d.png', p, q));

%%
p=2; q=0; d=0;
[M,EstCov,logL,info]=estimate(arima(p,d,q),y,'Display','off');
c=M.Constant; ar=cell2mat(M.AR); if numel(ar)<p, ar(end+1:p)=0; end
sig2=M.Variance; N=numel(y); k=(~isempty(c))+p+q+1;   % +1 for sigma^2
[AIC,BIC]=aicbic(logL,k,N);

fprintf('\nEstimated AR(%d): y_t = %.4f + %.4f y_{t-1} + %.4f y_{t-2} + \\eta_t,  \\eta_t~N(0,%.4f)\n',...
    p,c,ar(1),ar(2),sig2);
fprintf('AIC = %.2f, BIC = %.2f, logL = %.2f\n',AIC,BIC,logL);

rts=roots([1 -ar(:).']); stbl=all(abs(rts)>1);
fprintf('Stability (|roots|>1): %d; roots: %s\n',stbl, mat2str(rts,4));

% LaTeX string (copy/paste)
latex_str = sprintf(['$$\\hat y_t = %.3f + %.3f\\,y_{t-1} + %.3f\\,y_{t-2} + \\eta_t,\\quad ',...
                     '\\eta_t\\sim\\mathcal{N}(0,%.3f)$$'], c, ar(1), ar(2), sig2);
disp(latex_str);



%% LDF on residuals
order  = 1; points = 100; h = 0.6; maxlag = 30;
phi = ldf(res, order, points, h, maxlag);

fig = figure('Name','LDF on residuals','Color','w');
stem(0:maxlag, phi, 'filled'); grid on; yline(0,'k-');
xlabel('Lag k'); ylabel('$\phi(k)$');
title(sprintf('LDF on residuals (order=%d, h=%.2f)', order, h));
saveFigure(fig, figDir, 'P5_ldf_residuals.png');

%% Visualize dependence at a selectable lag
lagToView = 2;

fig = figure('Name','Res Dependence','Color','w');
et = res(1:end-lagToView); 
et_2 = res(1+lagToView:end);
scatter(et, et_2, 8, 'filled', 'MarkerFaceAlpha', 1); 
hold on; 
grid on;

xlabel(sprintf('$e_{t-%d}$', lagToView), 'Interpreter', 'latex'); 
ylabel('$e_t$', 'Interpreter', 'latex');
title(sprintf("$e_{t-%d}$ versus $e_{t}$", lagToView), "Interpreter", "latex");
saveFigure(fig, figDir, sprintf('P5_e_dependence_k%d.png', lagToView));

fig= figure('Name','Data Dependence','Color','w');
yt = y(1:end-lagToView); 
yt_2 = res(1+lagToView:end);
scatter(yt, yt_2, 8, 'filled', 'MarkerFaceAlpha', 1); 
hold on; 
grid on;

xlabel(sprintf('$y_{t-%d}$', lagToView), 'Interpreter', 'latex'); 
ylabel('$y_t$', 'Interpreter', 'latex');
title(sprintf("$y_{t-%d}$ versus $y_{t}$", lagToView), "Interpreter", "latex");
if ~exist('EstMdl','var')||isempty(EstMdl), EstMdl=estimate(arima(2,0,0),Y,'Display','off'); end
c=EstMdl.Constant; phi=EstMdl.AR; if iscell(phi), phi=cell2mat(phi); end
phi=phi(:); if numel(phi)<2, phi(end+1:2)=0; end
yline_ar2 = c + phi(1)*x1m + phi(2)*x2;
plot(x2,yline_ar2,'k:','LineWidth',1.8,'DisplayName','AR(2) slice');
legend()
saveFigure(fig, figDir, sprintf('P5_y_dependence_k%d.png', lagToView));


%% Fit Setar model on data 
Y=y(:); Y1=lagmatrix(Y,1); Y2=lagmatrix(Y,2);
idx=isfinite(Y)&isfinite(Y1)&isfinite(Y2);
Yc=Y(idx); Y1=Y1(idx); Y2=Y2(idx);

tau_vec=[-Inf,-4.13,-0.231,2.77,Inf]; K=4; betas=zeros(K,3); logL_setar=0;
sigma2 = nan(K,1);
RSS    = nan(K,1);
n_r    = nan(K,1);

for r = 1:K
    inR = (Y2 > tau_vec(r)) & (Y2 <= tau_vec(r+1));
    X   = [ones(sum(inR),1), Y1(inR), Y2(inR)];
    Yr  = Yc(inR);

    betas(r,:) = (X\Yr)';

    er   = Yr - X*betas(r,:)';
    n_r(r)   = sum(inR);
    RSS(r)   = er' * er;
    sigma2(r)= RSS(r) / n_r(r);                      

    logL_setar = logL_setar - 0.5 * n_r(r) * (1 + log(2*pi*sigma2(r)));
end

%% Setar fit 
fig=figure('Name','SETAR(4) regimes (k=2)','Color','w');
hData = scatter(Y2,Yc,10,'filled','MarkerFaceAlpha',0.3,'DisplayName','Data'); hold on; grid on;
x2=linspace(min(Y1),max(Y2),400)'; x1m=median(Y1,'omitnan'); C=lines(K);

for r=1:K
    m=(x2>tau_vec(r))&(x2<=tau_vec(r+1));
    if any(m)&&all(isfinite(betas(r,:)))
        plot(x2(m),betas(r,1)+betas(r,2)*x1m+betas(r,3)*x2(m),...
            '-','Color',C(r,:),'LineWidth',2,'DisplayName',sprintf('Regime %d',r));
    end
end
for r=2:numel(tau_vec)-1, xline(tau_vec(r),'k--','HandleVisibility','off'); end
xlabel('$y_{t-2}$','Interpreter','latex'); ylabel('$y_t$','Interpreter','latex');
title('SETAR Model Regimes','Interpreter','latex');

legend('Location','best');
saveFigure(fig, figDir, 'P5_setar_model_dependence.png');
%%
N = numel(Y); 
yhat_1step = nan(N, 1); 
R = nan(N, 1);

for t = 3:N
    r = find((Y(t-2) > tau_vec(1:end-1)) & (Y(t-2) <= tau_vec(2:end)), 1);
    if ~isempty(r)
        yhat_1step(t) = betas(r, 1) + betas(r, 2) * Y(t-1) + betas(r, 3) * Y(t-2);
        R(t) = r;
    end
end

%% Residual of one-step
res = Y - yhat_1step;
res = res(:);
reg_idx = nan(size(res));
for t = 3:numel(Y)
    reg_idx(t) = find((Y(t-2)>tau_vec(1:end-1)) & (Y(t-2)<=tau_vec(2:end)), 1);
end

z = res;                                   % standardized residuals
for r = 1:K
    m = (reg_idx == r);
    z(m) = res(m) ./ sqrt(sigma2(r));      % divide by regime std
end

% === Diagnostics on standardized residuals ===
fig=figure('Name','SETAR standardized residual diagnostics','Color','w');
subplot(3,2,1); plot(z,'b-'); grid on; title('Standardized residuals','Interpreter','latex'); xlabel('t','Interpreter','latex'); ylabel('$z_t$','Interpreter','latex');
subplot(3,2,2); histogram(z,'Normalization','pdf'); grid on; hold on
xg = linspace(min(z),max(z),400); plot(xg, normpdf(xg,0,1), 'k-', 'LineWidth',1.5);
title('Histogram of $z_t$ with $N(0,1)$','Interpreter','latex');
subplot(3,2,3); autocorr(z,'NumLags',30);  title('ACF of $z_t$','Interpreter','latex');
subplot(3,2,4); parcorr(z,'NumLags',30);   title('PACF of $z_t$','Interpreter','latex');
subplot(3,2,5); autocorr(z.^2,'NumLags',30); title('ACF of $z_t^2$','Interpreter','latex');
subplot(3,2,6); parcorr(z.^2,'NumLags',30);  title('PACF of $z_t^2$','Interpreter','latex');
saveFigure(fig, figDir, 'P5_setar_model_resid.png');



%% Setar of different regions
N=numel(Y); yhat_1step=nan(N,1); R=nan(N,1);
for t=3:N
    r=find((Y(t-2)>tau_vec(1:end-1))&(Y(t-2)<=tau_vec(2:end)),1);
    if ~isempty(r)
        yhat_1step(t)=betas(r,1)+betas(r,2)*Y(t-1)+betas(r,3)*Y(t-2);
        R(t)=r;
    end
end

fig=figure('Name','SETAR 1-step predictions + regimes','Color','w');
subplot(2,1,1); hold on; grid on;
plot(Y,'k-','DisplayName','y_t');
plot(yhat_1step,'b-','LineWidth',1.2,'DisplayName','\hat{y}_{t|t-1}');
legend('Location','best'); xlabel('t'); ylabel('Level'); title('SETAR 1-step predictions');

subplot(2,1,2); hold on; grid on;
stairs(R,'LineWidth',1.2);
ylim([0.5 K+0.5]); yticks(1:K); yticklabels(compose('Regime %d',1:K));
xlabel('t'); ylabel('Regime'); title('Active regime (by $y_{t-2}$)','Interpreter','latex');


saveFigure(fig, figDir, 'P5_setar_1_step_pred.png');


%% Local helper
function saveFigure(figHandle, directory, fileName)
    if nargin < 1 || isempty(figHandle) || ~ishandle(figHandle)
        return;
    end
    if nargin < 2 || isempty(directory)
        directory = pwd;
    end
    if nargin < 3 || isempty(fileName)
        warning('saveFigure:MissingFileName', 'File name missing; skipping export.');
        return;
    end
    targetPath = fullfile(directory, fileName);
    try
        exportgraphics(figHandle, targetPath, 'Resolution', 300);
    catch
        try
            saveas(figHandle, targetPath);
        catch saveErr
            warning('saveFigure:ExportFailed', 'Could not export figure %s: %s', fileName, saveErr.message);
        end
    end
end
