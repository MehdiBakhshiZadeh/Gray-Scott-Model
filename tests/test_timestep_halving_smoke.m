function r = test_timestep_halving_smoke()
%TEST_TIMESTEP_HALVING_SMOKE  Compare solutions at T using dt, dt/2, dt/4.
% The differences should decrease as dt decreases (basic consistency check).

% Base parameters
p = defaultParams();

% Keep the test short and stable
p.T = 50;
p.dt = 0.2;
p.plotEvery = 0;     % disable plotting
p.savePngEvery = 0;  % disable PNG export

% Build spatial operators once (same grid/operator for all time step sizes)
grid = buildGrid(p);
L = buildLaplacian2D(p, grid);

% Initial condition (fixed seed)
rng(p.seed, "twister");
[U0, V0] = initialCondition(p);
u0 = U0(:);
v0 = V0(:);

% Run dt, dt/2, dt/4 (same spatial discretization, only dt changes)
p1 = p;

p2 = p;
p2.dt = p.dt / 2;

p4 = p;
p4.dt = p.dt / 4;

[u1, v1] = simulateToFinal(u0, v0, L, p1);
[u2, v2] = simulateToFinal(u0, v0, L, p2);
[u4, v4] = simulateToFinal(u0, v0, L, p4);

% Compare v-field (commonly used for Gray–Scott patterns)
d12 = v1 - v2;
d24 = v2 - v4;

rms12 = norm(d12, 2) / sqrt(numel(d12));
rms24 = norm(d24, 2) / sqrt(numel(d24));

inf12 = norm(d12, inf);
inf24 = norm(d24, inf);

% Standardized result struct
r = makeResult("test_timestep_halving_smoke", p);
r.grid = [p.Nx, p.Ny];
r.dt   = p.dt;
r.T    = p.T;

r.metrics.rms_dt_dt2   = rms12;
r.metrics.rms_dt2_dt4  = rms24;
r.metrics.inf_dt_dt2   = inf12;
r.metrics.inf_dt2_dt4  = inf24;

% Pass criteria:
% 1) all metrics are finite
% 2) smaller dt -> smaller difference (allow a small slack factor)
finiteOK = isfinite(rms12) && isfinite(rms24) && isfinite(inf12) && isfinite(inf24);

% Slack factor avoids false failures due to sensitivity of nonlinear dynamics.
slack = 1.05;
decreaseOK = (rms24 <= slack * rms12) && (inf24 <= slack * inf12);

r.pass = finiteOK && decreaseOK;
r.thresholds.slack = slack;

r.notes = sprintf("T=%.2f, dt=%.3f: rms12=%.3e, rms24=%.3e, inf12=%.3e, inf24=%.3e", ...
    p.T, p.dt, rms12, rms24, inf12, inf24);
end

function [u, v] = simulateToFinal(u, v, L, p)
%SIMULATETOFINAL  Minimal time loop without plotting/export.

Nt = p.T / p.dt;
if abs(Nt - round(Nt)) > 1e-12
    error("T must be divisible by dt in this test.");
end
Nt = round(Nt);

S = []; % matrix mode (no stencil operator here)

for n = 1:Nt
    t = (n-1) * p.dt;  % current time before the Euler update
    [u, v, info] = stepEuler(u, v, L, S, p, t);

    if isfield(info,"hasNaNInf") && info.hasNaNInf
        error("NaN/Inf detected at step %d.", n);
    end
end
end

