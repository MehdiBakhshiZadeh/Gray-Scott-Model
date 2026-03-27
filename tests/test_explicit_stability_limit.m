function r = test_explicit_stability_limit()
%TEST_EXPLICIT_STABILITY_LIMIT  Stability-limit test for explicit Euler diffusion.
%
% This test verifies the conditional stability limit of the explicit Euler
% scheme for the diffusion part of the solver on a 2D periodic grid.
%
% The guideline asks, for explicit time integration methods, to investigate
% the critical time step and to compare behavior below and beyond that
% limit. This test does that in a controlled subcase of the implemented code.
%
% Strategy:
%   - Use a pure-diffusion subcase of the existing Gray-Scott solver:
%         F = 0, k = 0, v = 0  => reaction terms vanish identically
%   - Use a checkerboard perturbation in u, which corresponds to the
%     highest-frequency mode of the periodic 5-point Laplacian on an even grid.
%   - For this mode, the amplification factor of explicit Euler is known:
%         g = 1 - dt * Du * lambda_max
%     where
%         lambda_max = 4/hx^2 + 4/hy^2
%   - The critical time step is therefore
%         dt_crit = 2 / (Du * lambda_max)
%                 = 1 / (2*Du*(1/hx^2 + 1/hy^2))
%
% Expected behavior:
%   - dt <  dt_crit : perturbation decays
%   - dt =  dt_crit : perturbation magnitude is preserved
%   - dt >  dt_crit : perturbation grows (instability)
%
% Improvements over the basic version:
%   - Uses semilogy plots for clearer visualization of decay/growth.
%   - Checks agreement with theory over the full norm history, not only one step.

% --- Locate project folders robustly ---
thisDir = fileparts(mfilename("fullpath"));   % .../tests
rootDir = fileparts(thisDir);                 % project root
figDir  = fullfile(rootDir, "figures", "verification");

% Make this test runnable independently
if isempty(which("defaultParams"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(fullfile(rootDir, "tests")));
end

if ~exist(figDir, "dir")
    mkdir(figDir);
end

% --- Base parameters ---
base = defaultParams();
base.Nx = 128;
base.Ny = 128;
base.plotEvery = 0;
base.savePngEvery = 0;
base.diffusionMode = "matrix";

% Pure diffusion subcase of the same solver
base.F = 0.0;
base.k = 0.0;
base.Du = 2e-5;
base.Dv = 1e-5;   % irrelevant here because v is initialized to zero

% Keep periodic BC everywhere (default)
g = buildGrid(base);
L = buildLaplacian2D(base, g);

op = struct();
op.mode = "matrix";
op.L = L;
op.grid = g;

hx = g.hx;
hy = g.hy;

% Theoretical critical time step for explicit Euler + 2D diffusion
lambdaMax = 4/hx^2 + 4/hy^2;
dtCrit = 2 / (base.Du * lambdaMax);

% Test values: below, at, and above critical
factors = [0.90, 1.00, 1.10];
dtList  = factors * dtCrit;
labels  = ["0.9 dt_{crit}", "1.0 dt_{crit}", "1.1 dt_{crit}"];

% Fixed number of steps so the stability behavior is easy to compare
nSteps = 20;
steps  = 0:nSteps;

% Checkerboard perturbation = highest-frequency periodic mode
[Xi, Yi] = meshgrid(0:base.Nx-1, 0:base.Ny-1);
checker = (-1).^(Xi + Yi);

uMean = 1.0;
amp0  = 0.1;

U0 = uMean + amp0 * checker;
V0 = zeros(base.Ny, base.Nx);

u0 = U0(:);
v0 = V0(:);

% Containers
normHist   = zeros(numel(dtList), nSteps+1);
theoryHist = zeros(numel(dtList), nSteps+1);
gTheory    = zeros(size(dtList));
gObs       = zeros(size(dtList));
histErrAbs = zeros(size(dtList));
histErrRel = zeros(size(dtList));

for i = 1:numel(dtList)
    p = base;
    p.dt = dtList(i);
    p.T  = nSteps * p.dt;

    u = u0;
    v = v0;

    % Initial perturbation norm
    U = reshape(u, p.Ny, p.Nx);
    pert = U - uMean;
    normHist(i,1) = l2norm2D(pert, hx, hy);

    % Theoretical amplification
    gTheory(i) = abs(1 - p.dt * base.Du * lambdaMax);
    theoryHist(i,:) = normHist(i,1) * (gTheory(i) .^ steps);

    for n = 1:nSteps
        [u, v, info] = eulerStep(u, v, op, p, (n-1)*p.dt);
        if info.hasNaNInf
            error("NaN/Inf detected for dt = %.6g at step %d.", p.dt, n);
        end

        U = reshape(u, p.Ny, p.Nx);
        pert = U - uMean;
        normHist(i,n+1) = l2norm2D(pert, hx, hy);
    end

    % Observed one-step amplification in perturbation norm
    gObs(i) = normHist(i,2) / normHist(i,1);

    % Full-history agreement with theory
    histErrAbs(i) = max(abs(normHist(i,:) - theoryHist(i,:)));
    histErrRel(i) = histErrAbs(i) / normHist(i,1);

    fprintf(['case %d: dt = %.8g = %.2f dtCrit, ' ...
             'g_theory = %.6f, g_obs = %.6f, final/initial = %.6f, ' ...
             'histRelErr = %.3e\n'], ...
        i, p.dt, factors(i), gTheory(i), gObs(i), ...
        normHist(i,end)/normHist(i,1), histErrRel(i));
end

% --- Expected qualitative behavior ---
finalRatio = normHist(:,end) ./ normHist(:,1);

belowDecays = finalRatio(1) < 1.0;
atPreserved = abs(finalRatio(2) - 1.0) < 1e-10;
aboveGrows  = finalRatio(3) > 1.0;

% --- Theory agreement: one-step and full-history ---
tolTheory1Step = 1e-12;
tolTheoryHist  = 1e-10;

oneStepOK = all(abs(gObs - gTheory) < tolTheory1Step);
historyOK = all(histErrRel < tolTheoryHist);

% --- Plot ---
fig = figure("Color", "w");
ax = axes(fig);
hold(ax, "on");

normRatioHist = normHist ./ normHist(:,1);

semilogy(ax, steps, normRatioHist(1,:), "o-", "LineWidth", 1.6, "MarkerSize", 6);
semilogy(ax, steps, normRatioHist(2,:), "s-", "LineWidth", 1.6, "MarkerSize", 6);
semilogy(ax, steps, normRatioHist(3,:), "^-", "LineWidth", 1.6, "MarkerSize", 6);

grid(ax, "on");
xlabel(ax, "Euler step");
ylabel(ax, "normalized L2 norm");
title(ax, sprintf("Explicit Euler stability limit: dt_{crit}=%.6g", dtCrit), ...
    "Interpreter", "tex");
legend(ax, labels, "Location", "northwest");
set(ax, "FontSize", 12, "LineWidth", 1.0, "Box", "on");

% --- Metadata ---
meta = struct();
meta.test = "explicit_stability_limit";
meta.N = [base.Nx base.Ny];
meta.hx = hx;
meta.hy = hy;
meta.Du = base.Du;
meta.lambdaMax = lambdaMax;
meta.dtCrit = dtCrit;
meta.factors = factors;
meta.dtList = dtList;
meta.nSteps = nSteps;
meta.gTheory = gTheory;
meta.gObs = gObs;
meta.normHist = normHist;
meta.theoryHist = theoryHist;
meta.histErrAbs = histErrAbs;
meta.histErrRel = histErrRel;
meta.finalRatio = finalRatio;

outFig = fullfile(figDir, "explicit_stability_limit.png");
saveFigWithMeta(fig, outFig, meta);

assert(exist(outFig, "file") == 2, ...
    "Figure was not saved: %s", outFig);

% --- Result struct ---
r = makeResult("test_explicit_stability_limit", base);
r.metrics.N = base.Nx;
r.metrics.dtCrit = dtCrit;
r.metrics.dtList = dtList;
r.metrics.lambdaMax = lambdaMax;
r.metrics.gTheory = gTheory;
r.metrics.gObs = gObs;
r.metrics.finalRatio = finalRatio;
r.metrics.normHist = normHist;
r.metrics.theoryHist = theoryHist;
r.metrics.histErrAbs = histErrAbs;
r.metrics.histErrRel = histErrRel;

r.thresholds.behavior = "below decay, at preserve, above grow";
r.thresholds.oneStepTheoryTol = tolTheory1Step;
r.thresholds.historyTheoryTol = tolTheoryHist;

r.pass = belowDecays && atPreserved && aboveGrows && oneStepOK && historyOK;
r.figFiles = outFig;

r.notes = sprintf([ ...
    'N=%d; dtCrit=%.6g; finalRatios=[%.6f %.6f %.6f]; ' ...
    'gObs=[%.6f %.6f %.6f]; gTheory=[%.6f %.6f %.6f]; maxHistRelErr=%.3e'], ...
    base.Nx, dtCrit, ...
    finalRatio(1), finalRatio(2), finalRatio(3), ...
    gObs(1), gObs(2), gObs(3), ...
    gTheory(1), gTheory(2), gTheory(3), ...
    max(histErrRel));
end

function val = l2norm2D(A, hx, hy)
%L2NORM2D Discrete L2 norm on a uniform 2D grid.
val = sqrt(sum(A(:).^2) * hx * hy);
end