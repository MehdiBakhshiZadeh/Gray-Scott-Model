function r = test_mms_operator_convergence()
%TEST_MMS_OPERATOR_CONVERGENCE  Spatial convergence of the discrete Laplacian.
%
%   This test verifies the spatial accuracy of the discrete Laplacian
%   operator using the Method of Manufactured Solutions (MMS).
%   The discrete operator L is applied to an exact smooth field u*, and the
%   result is compared against the analytic Laplacian ?²u* on a sequence
%   of refined grids. The observed convergence rate should be second order.

% Locate project folders relative to this file
thisDir = fileparts(mfilename("fullpath"));   % .../tests
rootDir = fullfile(thisDir, "..");            % project root
figDir  = fullfile(rootDir, "figures", "verification");

% Ensure output directory exists
if ~exist(figDir, "dir")
    mkdir(figDir);
end

base = defaultParams();
m = mms_targets();

Ns = [32 64 128 256];
hs = zeros(size(Ns));

errU_L2 = zeros(size(Ns));
errV_L2 = zeros(size(Ns));

t0 = 0.3;   % fixed evaluation time (arbitrary)

for i = 1:numel(Ns)
    p = base;
    p.Nx = Ns(i);
    p.Ny = Ns(i);

    hx = p.Lx / p.Nx;
    hy = p.Ly / p.Ny;
    hs(i) = max(hx, hy);

    g = buildGrid(p);
    L = buildLaplacian2D(p, g);

    [X, Y] = meshgrid(g.x, g.y);

    % Exact fields and analytic Laplacians
    Uex = m.u(X, Y, t0);
    Vex = m.v(X, Y, t0);

    lapU_ex = m.lapU(X, Y, t0);
    lapV_ex = m.lapV(X, Y, t0);

    % Discrete Laplacians from matrix operator
    lapU_num = reshape(L * Uex(:), p.Ny, p.Nx);
    lapV_num = reshape(L * Vex(:), p.Ny, p.Nx);

    eU = lapU_num - lapU_ex;
    eV = lapV_num - lapV_ex;

    [eU_L2, ~] = fieldNorms2D(eU, hx, hy);
    [eV_L2, ~] = fieldNorms2D(eV, hx, hy);

    errU_L2(i) = eU_L2;
    errV_L2(i) = eV_L2;

    fprintf("N=%d, err(lapU)_L2=%.3e, err(lapV)_L2=%.3e\n", ...
        Ns(i), eU_L2, eV_L2);
end

% Observed convergence rates
pU = convergenceRate(hs, errU_L2);
pV = convergenceRate(hs, errV_L2);

% Plot convergence
fig = figure;
loglog(hs, errU_L2, "o-"); hold on;
loglog(hs, errV_L2, "s-"); grid on;
xlabel("h");
ylabel("L2 error of discrete Laplacian");
title(sprintf("Discrete Laplacian convergence: pU=%.2f, pV=%.2f", pU, pV));
legend("u: ||L u* - \nabla^2 u*||_{L2}", ...
       "v: ||L v* - \nabla^2 v*||_{L2}", ...
       "Location", "southwest");

meta = struct();
meta.test   = "mms_operator_convergence";
meta.Ns     = Ns;
meta.hs     = hs;
meta.errU_L2 = errU_L2;
meta.errV_L2 = errV_L2;
meta.t0     = t0;
meta.orderU = pU;
meta.orderV = pV;

figPath = fullfile(figDir, "mms_operator_convergence.png");
saveFigWithMeta(fig, figPath, meta);

% Standardized result struct
r = makeResult("test_mms_operator_convergence", base);
r.metrics.orderU = pU;
r.metrics.orderV = pV;

r.thresholds.orderBand = [1.7 2.3];
r.pass = (pU >= 1.7 && pU <= 2.3) && (pV >= 1.7 && pV <= 2.3);

r.figFiles = figPath;
r.notes = sprintf("Observed orders: pU=%.3f, pV=%.3f", pU, pV);
end
