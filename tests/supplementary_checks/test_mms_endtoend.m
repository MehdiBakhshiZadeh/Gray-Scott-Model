function r = test_mms_endtoend()
%TEST_MMS_ENDTOEND  End-to-end MMS verification on a single grid.
%
%   This test verifies the full Gray--Scott solver using the Method of
%   Manufactured Solutions (MMS). Analytic source terms are applied so that
%   a known exact solution (u*, v*) is recovered by the numerical method.
%
%   Errors are evaluated at the final time T using weighted L2 and Linf norms.
%   One representative error field (v-error) is saved for reporting.

% Locate this test folder and project root
thisFile = mfilename("fullpath");
thisDir  = fileparts(thisFile);   % .../tests/supplementary_checks
testsDir = fileparts(thisDir);    % .../tests
rootDir  = fileparts(testsDir);   % project root

% Make this test runnable independently
if isempty(which("defaultParams")) || isempty(which("eulerStep"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(testsDir));
end

% Local output folder for this test group
outDir = fullfile(thisDir, "outputs");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

% Parameters (keep test short and stable)
p = defaultParams();
p.T = 1.0;           % short final time for MMS
p.dt = 0.01;         % small dt to reduce time-discretization error
p.plotEvery = 0;     % disable plotting
p.savePngEvery = 0;  % disable PNG export

% Build grid and operator
grid = buildGrid(p);
L = buildLaplacian2D(p, grid);

op = struct();
op.mode = "matrix";
op.L = L;
op.grid = grid;

% Build 2D coordinate arrays (Ny x Nx)
[X, Y] = meshgrid(grid.x, grid.y);

% Override grid coordinates for MMS source evaluation
op.grid.x = X;
op.grid.y = Y;

% Manufactured solution and source terms
m = mms_targets();
s = mms_sources(p);

% Source function used inside eulerStep (returns 2D arrays)
p.sourceFcn = @(x,y,t,pp) deal(s.su(x,y,t), s.sv(x,y,t));

% Initial condition = exact solution at t = 0
U0 = m.u(X, Y, 0);
V0 = m.v(X, Y, 0);
u = U0(:);
v = V0(:);

% Time loop
Nt = p.T / p.dt;
if abs(Nt - round(Nt)) > 1e-12
    error("T must be divisible by dt in this test.");
end
Nt = round(Nt);

t = 0.0;

for n = 1:Nt
    [u, v, info] = eulerStep(u, v, op, p, t);

    if isfield(info, "hasNaNInf") && info.hasNaNInf
        error("NaN/Inf detected during MMS run at step %d.", n);
    end

    t = t + p.dt;
end

% Exact solution at final time
Uex = m.u(X, Y, p.T);
Vex = m.v(X, Y, p.T);

uNum = reshape(u, p.Ny, p.Nx);
vNum = reshape(v, p.Ny, p.Nx);

errU = uNum - Uex;
errV = vNum - Vex;

% Grid spacing
hx = p.Lx / p.Nx;
hy = p.Ly / p.Ny;

% Error norms
[euL2,  euInf] = fieldNorms2D(errU, hx, hy);
[evL2,  evInf] = fieldNorms2D(errV, hx, hy);

% Save one representative error figure (v-field)
fig = figure("Color", "w");
imagesc(errV);
axis image;
colorbar;
title(sprintf("MMS error in v at T=%.2f (Nx=%d, dt=%.3g)", p.T, p.Nx, p.dt));

meta = struct();
meta.test = "mms_endtoend";
meta.grid = [p.Nx p.Ny];
meta.T    = p.T;
meta.dt   = p.dt;
meta.Du   = p.Du;
meta.Dv   = p.Dv;
meta.F    = p.F;
meta.k    = p.k;
meta.mms = struct();
meta.mms.tag = "mms_endtoend";
meta.mms.note = "Manufactured sources enabled (handles omitted in metadata)";

figPath = fullfile(outDir, "mms_endtoend_errV.png");
saveFigWithMeta(fig, figPath, meta);

% Standardized result struct
r = makeResult("test_mms_endtoend", p);
r.grid = [p.Nx, p.Ny];
r.dt   = p.dt;
r.T    = p.T;

r.metrics.errU_L2   = euL2;
r.metrics.errU_Linf = euInf;
r.metrics.errV_L2   = evL2;
r.metrics.errV_Linf = evInf;

r.thresholds.errU_L2_max = 5e-3;
r.thresholds.errV_L2_max = 5e-3;

r.pass = isfinite(euL2) && isfinite(euInf) && ...
         isfinite(evL2) && isfinite(evInf) && ...
         (euL2 < r.thresholds.errU_L2_max) && ...
         (evL2 < r.thresholds.errV_L2_max);

r.figFiles = figPath;
r.notes = sprintf("errU(L2)=%.3e, errV(L2)=%.3e", euL2, evL2);
end