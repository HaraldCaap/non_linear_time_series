# Computer Exercise 2 — SETAR analysis workflow

This folder implements all parts of the 2011 computer exercise on
Self-Exciting Threshold Autoregressive (SETAR) models. The scripts are
written in MATLAB/Octave and share a common simulation based on the
SETAR$(2;1,1)$ process used in Computer Exercise 1.

Each part consists of two components:

1. A concise theoretical recap describing the goal, the underlying
   mathematics, and important practical considerations.
2. A MATLAB script (`partX.m`) that carries out the requested computation
   and persists intermediate artefacts (figures, tables, `.mat` files) for
   later use.

## Part 1 — Simulating a SETAR$(2;1,1)$ process

**Theory.** A SETAR$(2;p_L,p_H)$ process with delay $d$ and threshold $r$ is
piecewise linear:
\[
  y_t = \begin{cases}
    c_1 + \sum_{j=1}^{p_L} \phi_{1,j} y_{t-j} + \varepsilon_t, & y_{t-d} \le r, \\
    c_2 + \sum_{j=1}^{p_H} \phi_{2,j} y_{t-j} + \varepsilon_t, & y_{t-d} > r,
  \end{cases}
\]
with innovations $\varepsilon_t \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0,\sigma^2)$. The
process is geometrically ergodic when each regime is contractive (eigenvalues
inside the unit circle) and the regimes communicate through the threshold. We
reuse the parameters from Exercise 1:
$\phi_{1,1}=0.6$, $\phi_{2,1}=0.9$, $c_1=0.4$, $c_2=-0.3$, $d=1$, $r=0.2$.

**Implementation.** `part1.m` draws a 1,200 point series (discarding 300 burn-in
observations), stores it in `data/part1_sim.mat`, and produces diagnostic plots
for the realised regimes.

## Part 2 — Least squares estimation with grid-search over thresholds

**Theory.** Threshold estimation under deterministic regime memberships is a
non-linear least squares problem. Conditional on a candidate threshold $r$, the
model reduces to two separate linear regressions. We minimise the global sum of
squared errors (SSE) over a grid $\mathcal{R}$ constructed from sample
quantiles of the threshold variable. Information criteria such as the BIC,
\(\mathrm{BIC} = n \log(\widehat\sigma^2) + k \log n\), provide secondary
model selection diagnostics.

**Implementation.** `part2.m` searches 40 quantile-based thresholds, estimates
$\widehat{\theta}_1$ and $\widehat{\theta}_2$ for the optimal $r$, saves the
model (`results/part2_model.mat`), and plots the SSE profile to verify the
uniqueness of the optimum.

## Part 3 — Regime diagnostics and ergodic properties

**Theory.** Once the model is estimated we can analyse regime-specific
behaviour by conditioning on the estimated memberships. The empirical regime
moments check whether the two subsystems capture distinct dynamics. Stability
requires $|\phi_{j}| < 1$ in each regime. Long-run properties (stationary mean,
regime probabilities) are approximated via Monte Carlo simulation using the
estimated model and innovation variance $\widehat{\sigma}^2 = \widehat{\text{SSE}}/n$.

**Implementation.** `part3.m` computes sample moments, validates the contraction
conditions, performs a long Monte Carlo run, and compares empirical vs
simulated distributions of $y_t$.

## Part 4 — Forecasting versus a linear AR benchmark

**Theory.** Forecasting in a SETAR model uses regime-dependent recursion. Given
observations $y_{t-1},\ldots$, the $h$-step forecast iterates the conditional
mean appropriate for the regime determined by $y_{t-d}$. We compare this
nonlinear predictor to a linear AR(1) model fitted by OLS. Forecast accuracy is
measured by root mean squared error (RMSE) both in-sample and on a hold-out
segment.

**Implementation.** `part4.m` generates recursive one-step forecasts for the
entire sample, computes RMSEs, and produces an out-of-sample comparison plot for
the last 200 observations.

## Part 5 — Residual bootstrap for threshold uncertainty

**Theory.** Analytical standard errors for the threshold are complicated, so a
residual bootstrap is an accessible alternative. We fix the estimated
parameters, resample centered residuals $\{\widehat\varepsilon_t\}$ with
replacement, simulate new series using the deterministic dynamics, and
re-estimate the model to approximate the sampling distribution of the threshold
and regime parameters. Percentile intervals summarize uncertainty.

**Implementation.** `part5.m` performs 200 bootstrap replications, stores
summary statistics in `results/part5_bootstrap_summary.txt`, and visualises the
bootstrap distributions for the threshold and leading autoregressive
coefficients.

## Running the workflow

Execute the scripts sequentially from MATLAB or Octave:

```matlab
cd('path/to/non_linear_time_series/ce2');
part1; % simulate and save data
part2; % estimate threshold model
part3; % diagnostics
part4; % forecasting comparison
part5; % bootstrap analysis
```

Figures are written to `figs/` and textual outputs to `results/`. All scripts
are deterministic conditional on the fixed RNG seeds and rely only on base
MATLAB/Octave functionality.
