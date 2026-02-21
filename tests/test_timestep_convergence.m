function r = test_timestep_convergence()
%TEST_TIMESTEP_CONVERGENCE  Temporal convergence check for explicit Euler.
%
%   Runs the Gray–Scott model to the same final time using three time step
%   sizes: dt, dt/2, dt/4 on the same spatial grid and same initial condition.
%   Computes successive differences in v:
%       d12 = ||v_dt   - v_dt/2||
%       d24 = ||v_dt/2 - v_dt/4||
%   Checks d24 < d12 and estimates observed order:
%       p ~ log(d12/d24)/log(2)
%   For explicit Euler, expected temporal order is ~1.

% Ensure output directory exists
thisDir = fileparts(mfilename("fullpath"));   % .../tests
rootDir = fullfile(thisDir, "..");            % project root
figDir  = fullfile(rootDir, "figures", "verification");
if ~exist(figDir, "dir")
    mkdir(figDir);
end

% --- Setup parameters (stable + short) ---
base = defaultParams();
base.Nx = 128;
base.Ny = 128;
base.T  = 50;
base.plotEvery = 0;     % disable plotting
base.savePngEvery = 0;  % disable PNG export

dtList = [0.2, 0.1, 0.05];  % dt, dt/2, dt/4 (must divide T)

% Use the same IC across all runs
rng(base.seed, "twister");
[U0, V0] = initialCondition(base);
u0 = U0(:);
v0 = V0(:);

% --- Run three simulations ---
[U1, V1] = runToFinal(base, dtList(1), u0, v0);
[U2, V2] = runToFinal(base, dtList(2), u0, v0);
[U4, V4] = runToFinal(base, dtList(3), u0, v0);

% --- Differences (use v-field) ---
dV12 = V1 - V2;
dV24 = V2 - V4;

hx = base.Lx / base.Nx;
hy = base.Ly / base.Ny;

[d12_L2, d12_Inf] = fieldNorms2D(dV12, hx, hy);
[d24_L2, d24_Inf] = fieldNorms2D(dV24, hx, hy);

% Observed time order indicator
pL2  = log(d12_L2 / d24_L2) / log(2);
pInf = log(d12_Inf / d24_Inf) / log(2);

fprintf("dt=%.4g: d12(L2)=%.3e, d24(L2)=%.3e, p(L2)=%.3f\n", dtList(1), d12_L2, d24_L2, pL2);
fprintf("dt=%.4g: d12(Inf)=%.3e, d24(Inf)=%.3e, p(Inf)=%.3f\n", dtList(1), d12_Inf, d24_Inf, pInf);

% --- Plot (log-log): successive differences at their associated dt ---
xVals = [dtList(1), dtList(2)];
yVals = [d12_L2, d24_L2];

fig = figure("Color", "w");
ax = axes(fig);

loglog(ax, xVals, yVals, "o-", "LineWidth", 1.6, "MarkerSize", 7);
grid(ax, "on");

xlabel(ax, "\Delta t", "Interpreter", "tex");
ylabel(ax, "Successive difference in v (L2)");
title(ax, sprintf("Time-step refinement (N=%d, T=%.1f): p=%.2f", base.Nx, base.T, pL2), ...
      "Interpreter", "none");

set(ax, "FontSize", 12, "LineWidth", 1.0, "Box", "on");
set(findall(fig, "Type", "text"), "FontSize", 12);

xlim(ax, [min(xVals)*0.9, max(xVals)*1.1]);
ylim(ax, [min(yVals)*0.9, max(yVals)*1.1]);

% --- Metadata ---
meta = struct();
meta.test   = "timestep_convergence";
meta.N      = [base.Nx base.Ny];
meta.T      = base.T;
meta.dtList = dtList;
meta.d12_L2 = d12_L2;  meta.d24_L2 = d24_L2;
meta.d12_Inf = d12_Inf; meta.d24_Inf = d24_Inf;
meta.pL2 = pL2; meta.pInf = pInf;

% --- Save figure with metadata ---
outFig = fullfile(figDir, "timestep_convergence.png");
saveFigWithMeta(fig, outFig, meta);

% Fail fast if saving didn't happen
assert(exist(outFig, "file") == 2, "Figure was not saved: %s", outFig);

% --- Result struct ---
r = makeResult("test_timestep_convergence", base);
r.dt = dtList(1);
r.T  = base.T;

r.metrics.d12_L2 = d12_L2;
r.metrics.d24_L2 = d24_L2;
r.metrics.pL2    = pL2;

r.metrics.d12_Inf = d12_Inf;
r.metrics.d24_Inf = d24_Inf;
r.metrics.pInf    = pInf;

% Pass conditions:
r.thresholds.requireDecrease = true;
r.thresholds.pBand = [0.7 1.3];

decrease = (d24_L2 < d12_L2) && (d24_Inf < d12_Inf);
orderOK  = (pL2 >= r.thresholds.pBand(1) && pL2 <= r.thresholds.pBand(2));

r.pass = decrease && orderOK;

r.figFiles = outFig;
r.notes = sprintf("d12(L2)=%.3e, d24(L2)=%.3e, pL2=%.3f; pInf=%.3f", d12_L2, d24_L2, pL2, pInf);
end

function [U, V] = runToFinal(base, dt, u0, v0)
%RUNTOFINAL  Run Gray–Scott from the same IC to time T with time step dt.
% Uses baseline sparse Laplacian (matrix mode) for verification.

p = base;
p.dt = dt;

Nt = p.T / p.dt;
if abs(Nt - round(Nt)) > 1e-12
    error("T must be divisible by dt. Got T=%.6g, dt=%.6g.", p.T, p.dt);
end
Nt = round(Nt);

g = buildGrid(p);
L = buildLaplacian2D(p, g);

S = [];   % matrix diffusion in verification tests

u = u0;
v = v0;

t = 0.0;
for n = 1:Nt
    [u, v, info] = eulerStep(u, v, L, S, p, t);

    if isfield(info,"hasNaNInf") && info.hasNaNInf
        error("NaN/Inf detected at step %d for dt=%.6g.", n, dt);
    end

    t = t + p.dt;
end

U = reshape(u, p.Ny, p.Nx);
V = reshape(v, p.Ny, p.Nx);
end
