function r = test_mms_spatial_convergence()
%TEST_MMS_SPATIAL_CONVERGENCE  MMS grid refinement and observed spatial order.
%
%   This test performs a full MMS (end-to-end) simulation on a sequence of
%   refined grids and measures the error at the final time T.
%
%   To keep the time-discretization error small relative to spatial error,
%   the time step is chosen as:
%       dt = Cdt * h^2
%   and then adjusted so that T is exactly divisible by dt.
%
%   The observed convergence rate is computed from the L2 error of v.
%   For a second-order Laplacian and smooth manufactured solution, the
%   expected spatial order is approximately 2.

% Locate project folders relative to this file
thisDir = fileparts(mfilename("fullpath"));   % .../tests
rootDir = fullfile(thisDir, "..");            % project root
figDir  = fullfile(rootDir, "figures", "verification");

% Ensure output directory exists
if ~exist(figDir, "dir")
    mkdir(figDir);
end

base = defaultParams();
base.T = 1.0;
base.plotEvery = 0;     % disable plotting
base.savePngEvery = 0;  % disable PNG export

Ns  = [32 64 128 256];  % refinement levels
Cdt = 0.1;              % dt = Cdt*h^2 (keeps time error small)

errV_L2 = zeros(size(Ns));
hs      = zeros(size(Ns));

m = mms_targets();

for i = 1:numel(Ns)
    p = base;
    p.Nx = Ns(i);
    p.Ny = Ns(i);

    % Grid spacing (periodic grid on [0,L))
    hx = p.Lx / p.Nx;
    hy = p.Ly / p.Ny;
    h  = max(hx, hy);
    hs(i) = h;

    % Choose dt scaled with h^2 and ensure T divisible by dt
    p.dt = Cdt * h^2;
    Nt = round(p.T / p.dt);
    p.dt = p.T / Nt;   % adjust to divide T exactly

    g = buildGrid(p);
    L = buildLaplacian2D(p, g);

    % Build 2D coordinate arrays
    [X, Y] = meshgrid(g.x, g.y);

    % Forcing (manufactured source terms)
    s = mms_sources(p);
    p.grid = struct();
    p.grid.x = X;
    p.grid.y = Y;
    p.sourceFcn = @(x,y,t,pp) deal(s.su(x,y,t), s.sv(x,y,t));

    % Initial condition = exact at t=0
    U0 = m.u(X, Y, 0);
    V0 = m.v(X, Y, 0);
    u = U0(:);
    v = V0(:);

    % Time loop
    t = 0.0;
    for n = 1:Nt
        [u, v, info] = stepEuler(u, v, L, p, t);
        if info.hasNaNInf
            error("NaN/Inf detected at N=%d, step=%d.", Ns(i), n);
        end
        t = t + p.dt;
    end

    % Exact at final time
    Vex = m.v(X, Y, p.T);
    vNum = reshape(v, p.Ny, p.Nx);
    errV = vNum - Vex;

    [eL2, ~] = fieldNorms2D(errV, hx, hy);
    errV_L2(i) = eL2;

    fprintf("N=%d, dt=%.3g, errV_L2=%.3e\n", Ns(i), p.dt, eL2);
end

% Observed order (use L2 error)
p_est = convergenceRate(hs, errV_L2);

% Plot log-log
fig = figure;
loglog(hs, errV_L2, "o-"); grid on;
xlabel("h");
ylabel("||e_v||_{L2}");
title(sprintf("MMS spatial convergence (observed order %.2f)", p_est));

meta = struct();
meta.test = "mms_spatial_convergence";
meta.Ns = Ns;
meta.hs = hs;
meta.errV_L2 = errV_L2;
meta.T = base.T;
meta.Cdt = Cdt;
meta.observedOrder = p_est;

figPath = fullfile(figDir, "mms_spatial_convergence.png");
saveFigWithMeta(fig, figPath, meta);

% Standardized result
r = makeResult("test_mms_spatial_convergence", base);
r.T = base.T;
r.metrics.orderV_L2 = p_est;

% Also store the data (handy for runner output)
r.metrics.errV_L2_N32  = errV_L2(1);
r.metrics.errV_L2_N64  = errV_L2(2);
r.metrics.errV_L2_N128 = errV_L2(3);
r.metrics.errV_L2_N256 = errV_L2(4);

% Pass criterion: close to 2 (tolerant band)
r.thresholds.orderTarget = 2.0;
r.thresholds.orderBand   = [1.7 2.3];
r.pass = (p_est >= r.thresholds.orderBand(1)) && (p_est <= r.thresholds.orderBand(2));

r.figFiles = figPath;
r.notes = sprintf("Observed order (v, L2) = %.3f", p_est);
end
