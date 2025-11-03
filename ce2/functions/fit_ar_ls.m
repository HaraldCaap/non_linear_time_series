function model = fit_ar_ls(y, p)
%FIT_AR_LS Ordinary least squares estimation of an AR(p) model with intercept.
%   MODEL = FIT_AR_LS(y, p) returns a struct containing:
%       .coeff  - [c phi_1 ... phi_p]
%       .sigma2 - residual variance estimate
%       .resid  - in-sample residuals (length numel(y) - p)

assert(isvector(y), 'y must be a vector');
y = y(:);
T = numel(y);

if p < 0 || p ~= floor(p)
    error('p must be a non-negative integer');
end

if p == 0
    coeff = mean(y);
    resid = y - coeff;
    sigma2 = var(resid, 1);
    model = struct('coeff', coeff, 'sigma2', sigma2, 'resid', resid);
    return;
end

Y = y(p+1:end);
X = ones(T - p, p + 1);
for j = 1:p
    X(:, j + 1) = y(p + 1 - j : T - j);
end

coeff = X \ Y;
resid = Y - X * coeff;
sigma2 = (resid' * resid) / (T - p - (p + 1));

model = struct('coeff', coeff', 'sigma2', sigma2, 'resid', resid);
end
