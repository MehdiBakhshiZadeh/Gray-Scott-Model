function r = test_grayscott_grid_refinement()
%TEST_GRAYSCOTT_GRID_REFINEMENT  Grid refinement study for Gray–Scott (explicit Euler).
%
%   This test performs a grid refinement study for the Gray–Scott model using
%   explicit Euler time integration. The time step is scaled as dt ~ h^2 to
%   limit time-discretization effects.
%
%   Because Gray–Scott is nonlinear and pattern-forming, strict second-order
%   convergence of the solution fields in L2 is not guaranteed (phase shifts
%   can dominate). Therefore, this test:
%     - checks that refinement reduces differences relative to the finest grid,
%     - reports observed convergence rates for reference,
%     - uses a simple, defensible pass criterion.

% --- Locate project folders robustly ---
thisDir = fileparts(mfilename("fullpath"));   % .../tests
rootDir = fullfile(thisDir, "..");
figDir  = fullfile(rootDir, "figures", "verification");

if ~exist(figDir, "dir")
    mkdir(figDir);
end

% --- Base parameters ---
base = defaultParams();
base.T = 1.0;
base.plotEvery = 0;     % disable plotting
base.savePngEvery = 0;  % disable PNG export

% Grid levels (finest grid used as reference)
Ns  = [64 128 256 512];
dt0 = 0.01;   % baseline dt for N=64 (coarsest grid)

U = cell(size(Ns));
V = cell(size(Ns));
grids = cell(size(Ns));

% ---- Run simulations on all grids ----
for i = 1:numel(Ns)
    p = base;
    p.Nx = Ns(i);
    p.Ny = Ns(i);

    % dt scaling anchored at N=64: dt = dt0 * (64/N)^2
    p.dt = dt0 * (64 / Ns(i))^2;

    % Adjust dt so that T is divisible by dt
    Nt = round(p.T / p.dt);
    p.dt = p.T / Nt;

    g = buildGrid(p);
    L = buildLaplacian2D(p, g);

    rng(p.seed, "twister");
    [U0, V0] = initialCondition(p);
    u = U0(:);
    v = V0(:);

    for n = 1:Nt
        [u, v, info] = eulerStep(u, v, L, p);
        if info.hasNaNInf
            error("NaN/Inf detected at N=%d, step=%d.", Ns(i), n);
        end
    end

    U{i} = reshape(u, p.Ny, p.Nx);
    V{i} = reshape(v, p.Ny, p.Nx);
    grids{i} = g;

    fprintf("Finished N=%d (dt=%.6g, Nt=%d)\n", Ns(i), p.dt, Nt);
end

% ---- Interpolate finest-grid reference onto coarser grids ----
refIdx = numel(Ns);
Uref = U{refIdx};
Vref = V{refIdx};
gref = grids{refIdx};
[Xref, Yref] = meshgrid(gref.x, gref.y);

Uref_on = cell(size(Ns));
Vref_on = cell(size(Ns));
Uref_on{refIdx} = Uref;
Vref_on{refIdx} = Vref;

interpMethod = "cubic";

for i = 1:refIdx-1
    gi = grids{i};
    [Xi, Yi] = meshgrid(gi.x, gi.y);

    Uref_on{i} = interp2(Xref, Yref, Uref, Xi, Yi, interpMethod);
    Vref_on{i} = interp2(Xref, Yref, Vref, Xi, Yi, interpMethod);

    fprintf("Interpolated reference onto N=%d (%s)\n", Ns(i), interpMethod);
end

% ---- Error norms on coarse grids ----
nCoarse = refIdx - 1;
hs = zeros(1, nCoarse);
Eu = zeros(1, nCoarse);
Ev = zeros(1, nCoarse);

for i = 1:nCoarse
    p = base;
    p.Nx = Ns(i);
    p.Ny = Ns(i);

    hx = p.Lx / p.Nx;
    hy = p.Ly / p.Ny;
    hs(i) = max(hx, hy);

    eU = U{i} - Uref_on{i};
    eV = V{i} - Vref_on{i};

    [Eu(i), ~] = fieldNorms2D(eU, hx, hy);
    [Ev(i), ~] = fieldNorms2D(eV, hx, hy);

    fprintf("N=%d: ||eU||_L2=%.3e, ||eV||_L2=%.3e\n", Ns(i), Eu(i), Ev(i));
end

% ---- Observed convergence rates ----
pU = convergenceRate(hs, Eu);
pV = convergenceRate(hs, Ev);

% ---- Plot ----
fig = figure("Color", "w");
loglog(hs, Eu, "o-", "LineWidth", 1.6, "MarkerSize", 7); hold on;
loglog(hs, Ev, "s-", "LineWidth", 1.6, "MarkerSize", 7);
grid on;
xlabel("h");
ylabel("L2 difference vs finest-grid reference");
title(sprintf("Gray–Scott grid refinement (Euler, T=%.1f): pU=%.2f, pV=%.2f", ...
    base.T, pU, pV), "Interpreter", "none");
legend("u", "v", "Location", "southwest");

% ---- Save figure with metadata ----
meta = struct();
meta.test = "grayscott_grid_refinement";
meta.Ns = Ns;
meta.hs = hs;
meta.Eu = Eu;
meta.Ev = Ev;
meta.interpMethod = interpMethod;
meta.T = base.T;
meta.dtScaling = "dt ~ h^2";

outFig = fullfile(figDir, "grayscott_grid_refinement.png");
saveFigWithMeta(fig, outFig, meta);

assert(exist(outFig, "file") == 2, ...
    "Figure was not saved: %s", outFig);

% ---- Result struct ----
r = makeResult("test_grayscott_grid_refinement", base);
r.metrics.orderU = pU;
r.metrics.orderV = pV;
r.metrics.hs = hs;
r.metrics.Eu = Eu;
r.metrics.Ev = Ev;

% Pass criteria (simple, defensible)
r.thresholds.minOrder = 0.8;

decreaseU = all(diff(Eu) < 0);
decreaseV = all(diff(Ev) < 0);
orderOK   = (pU >= r.thresholds.minOrder) && (pV >= r.thresholds.minOrder);

r.pass = decreaseU && decreaseV && orderOK;
r.figFiles = outFig;
r.notes = sprintf("Euler with dt~h^2, ref N=%d: pU=%.3f, pV=%.3f", Ns(refIdx), pU, pV);
end
