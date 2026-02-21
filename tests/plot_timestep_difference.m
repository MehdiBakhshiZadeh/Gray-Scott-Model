clear; clc; close all;

%PLOT_TIMESTEP_DIFFERENCE  Visualize the effect of halving the time step.
%
%   This script runs the Gray–Scott model up to time T using dt and dt/2
%   (same spatial grid), then visualizes the difference:
%       ?v = v_dt - v_dt/2
%
%   Output:
%     The figure and metadata are saved to:
%       figures/verification/timestep_halving_diffV.png
%
%   The script does not depend on the current working directory.

% Locate project folders relative to this script
thisDir = fileparts(mfilename("fullpath"));      % .../tests
rootDir = fullfile(thisDir, "..");              % project root

% Add required paths (absolute)
addpath(fullfile(rootDir, "src"));
addpath(fullfile(thisDir, "utils"));

% Ensure output directory exists
figDir = fullfile(rootDir, "figures", "verification");
if ~exist(figDir, "dir")
    mkdir(figDir);
end

% Parameters
p = defaultParams();
p.T = 50;
p.dt = 0.2;
p.plotEvery = 0;     % disable plotting
p.savePngEvery = 0;  % disable PNG export

% Spatial operators (fixed across runs)
grid = buildGrid(p);
L = buildLaplacian2D(p, grid);

% Initial condition (fixed seed)
rng(p.seed, "twister");
[U0, V0] = initialCondition(p);
u0 = U0(:);
v0 = V0(:);

% Run with dt
[uA, vA] = simulateToFinal(u0, v0, L, p);

% Run with dt/2
p2 = p;
p2.dt = p.dt / 2;
[uB, vB] = simulateToFinal(u0, v0, L, p2);

% Difference field (v)
dv = vA - vB;
DV = reshape(dv, p.Ny, p.Nx);

% Plot
fig = figure;
imagesc(DV);
axis image; colorbar;
title(sprintf("\\Delta v = v_{dt} - v_{dt/2} at T=%.0f (dt=%.3f)", p.T, p.dt));

% Save figure + metadata (consistent, reproducible)
meta = struct();
meta.script = "plot_timestep_difference";
meta.grid   = [p.Nx p.Ny];
meta.T      = p.T;
meta.dtA    = p.dt;
meta.dtB    = p2.dt;
meta.seed   = p.seed;

outPath = fullfile(figDir, "timestep_halving_diffV.png");
saveFigWithMeta(fig, outPath, meta);

% Print norms for reference
err2  = norm(dv, 2) / sqrt(numel(dv));
errInf = norm(dv, inf);
fprintf("Saved figures/verification/timestep_halving_diffV.png | RMS=%.3e, inf=%.3e\n", err2, errInf);

function [u, v] = simulateToFinal(u, v, L, p)
%SIMULATETOFINAL  Minimal time loop without plotting/export.

Nt = p.T / p.dt;
assert(abs(Nt - round(Nt)) < 1e-12, "T must be divisible by dt.");
Nt = round(Nt);

for n = 1:Nt
    [u, v, info] = eulerStep(u, v, L, p);
    if info.hasNaNInf
        error("NaN/Inf detected at step %d.", n);
    end
end
end
