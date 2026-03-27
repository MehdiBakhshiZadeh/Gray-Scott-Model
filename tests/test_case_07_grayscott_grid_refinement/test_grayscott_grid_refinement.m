function r = test_grayscott_grid_refinement()
%TEST_GRAYSCOTT_GRID_REFINEMENT
% Grid refinement study for Gray-Scott with a smooth, grid-independent
% initial condition. Designed to measure observed convergence more cleanly
% than the original pattern-seeding setup.

% Locate this test folder and project root
thisFile = mfilename("fullpath");
thisDir  = fileparts(thisFile);   % .../tests/test_case_07_grayscott_grid_refinement
testsDir = fileparts(thisDir);    % .../tests
rootDir  = fileparts(testsDir);   % project root

% Make this test runnable independently
if isempty(which("defaultParams")) || isempty(which("eulerStep"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(testsDir));
end

% Local output folder for this test
outDir = fullfile(thisDir, "outputs");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

% --- Base parameters ---
base = defaultParams();
base.T = 0.2;           % shorter time reduces phase-separation effects
base.plotEvery = 0;
base.savePngEvery = 0;
base.icPerturb = 0.0;   % IMPORTANT: no random noise in this test

% Grid levels (finest grid used as reference)
Ns  = [64 128 256 512 1024 2048];
dt0 = 0.01;   % dt at N=64, then dt ~ h^2

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

    op = struct();
    op.mode = "matrix";
    op.L = L;
    op.grid = g;

    % Smooth, grid-independent IC in physical coordinates
    [X, Y] = meshgrid(g.x, g.y);
    [U0, V0] = smoothGrayScottIC(X, Y, p);

    u = U0(:);
    v = V0(:);

    for n = 1:Nt
        [u, v, info] = eulerStep(u, v, op, p);
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
title(sprintf("Gray-Scott smooth grid refinement (Euler, T=%.2f): pU=%.2f, pV=%.2f", ...
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
meta.icType = "smooth grid-independent Gaussian";

outFig = fullfile(outDir, "grayscott_grid_refinement.png");
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

r.thresholds.minOrder = 1.5;

decreaseU = all(diff(Eu) < 0);
decreaseV = all(diff(Ev) < 0);
orderOK   = (pU >= r.thresholds.minOrder) && (pV >= r.thresholds.minOrder);

r.pass = decreaseU && decreaseV && orderOK;
r.figFiles = outFig;
r.notes = sprintf("Smooth IC, Euler with dt~h^2, ref N=%d: pU=%.3f, pV=%.3f", ...
    Ns(refIdx), pU, pV);
end

function [U0, V0] = smoothGrayScottIC(X, Y, p)
% Smooth grid-independent initial condition in physical coordinates.

xc = 0.5 * p.Lx;
yc = 0.5 * p.Ly;
sigma = 0.08 * min(p.Lx, p.Ly);

R2 = (X - xc).^2 + (Y - yc).^2;
phi = exp(-R2 / (2 * sigma^2));

U0 = 1.0 - 0.50 * phi;
V0 = 0.25 * phi;
end