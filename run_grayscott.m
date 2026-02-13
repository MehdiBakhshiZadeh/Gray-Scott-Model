clear; clc; close all;

% Path of this script (root of project)
rootDir = fileparts(mfilename("fullpath"));

% Ensure src is visible (absolute path, no cd side-effect)
addpath(fullfile(rootDir, "src"));

% Load default parameters
p = defaultParams();

% Fix random seed (single source of randomness for reproducibility)
rng(p.seed, "twister");

% Create a unique name for this run
timestamp = datestr(now, "yyyy-mm-dd_HH-MM-SS");
runName   = timestamp + "_" + p.caseName;

% Define output directory (absolute, safe)
outDir = fullfile(rootDir, "results", runName);
if ~exist(outDir, "dir")
    mkdir(outDir);
end

% Save metadata for reproducibility
meta.timestamp = timestamp;
meta.parameters = p;
meta.matlabVersion = version;
meta.platform = computer;
save(fullfile(outDir, "meta.mat"), "meta");

% Write a readable log file (safe close)
logPath = fullfile(outDir, "log.txt");
fid = fopen(logPath, "w");
assert(fid ~= -1, "Could not open log file: %s", logPath);
c = onCleanup(@() fclose(fid));

fprintf(fid, "Run name: %s\n", runName);
fprintf(fid, "Timestamp: %s\n", timestamp);
fprintf(fid, "MATLAB: %s\n", version);
fprintf(fid, "Platform: %s\n", computer);
fprintf(fid, "Grid: Nx=%d, Ny=%d\n", p.Nx, p.Ny);
fprintf(fid, "Time: dt=%g, T=%g\n", p.dt, p.T);
fprintf(fid, "Model: Du=%g, Dv=%g, F=%g, k=%g\n", p.Du, p.Dv, p.F, p.k);
fprintf(fid, "Output directory: %s\n", outDir);

% ---- Build grid and operators ----
grid = buildGrid(p);

mode = "matrix";
if isfield(p,"diffusionMode") && ~isempty(p.diffusionMode)
    mode = string(p.diffusionMode);
end

if mode == "stencil"
    S = buildStencil2D(p, grid);
    L = [];

elseif mode == "full"
    allowed = [32 64 128];
    assert(any(p.Nx == allowed) && any(p.Ny == allowed), ...
        "Full Laplacian only allowed for Nx,Ny in {32,64,128}. Got Nx=%d, Ny=%d.", p.Nx, p.Ny);

    L = full(buildLaplacian2D(p, grid));
    S = [];

else
    L = buildLaplacian2D(p, grid);
    S = [];
end



% Optional stability suggestion for explicit diffusion
Dmax = max(p.Du, p.Dv);
dt_stab = 1 / (2*Dmax*(1/grid.hx^2 + 1/grid.hy^2));
fprintf(fid, "Stability suggestion (explicit diffusion): dt <= %.6g\n", dt_stab);

if p.dt > dt_stab
    warning("dt=%.3g may be unstable for explicit diffusion (suggest dt <= %.3g).", p.dt, dt_stab);
end

% ---- Initial condition (2D -> vector) ----
[U0, V0] = initialCondition(p);
u = U0(:);
v = V0(:);

% ---- Time loop settings ----
Nt_exact = p.T / p.dt;
assert(abs(Nt_exact - round(Nt_exact)) < 1e-12, "T must be divisible by dt.");
Nt = round(Nt_exact);
t = 0.0;
fprintf(fid, "Steps: Nt=%d\n", Nt);
fprintf(fid, "plotEvery=%d, savePngEvery=%d\n", p.plotEvery, p.savePngEvery);

% ---- Optional live figure ----
doPlot = (isfield(p, "plotEvery") && p.plotEvery > 0);
if doPlot
    figure;
    plotState(u, v, p, t);
end

% ---- Time stepping loop ----
for n = 1:Nt
    [u, v, info] = stepEuler(u, v, L, S, p);
    t = n * p.dt;

    % Safety: stop if NaN/Inf appears
    if info.hasNaNInf
        error("NaN/Inf detected at step %d (t=%.3f). Stopping.", n, t);
    end

    % Live plot (avoid plotting every step)
    if doPlot && (mod(n, p.plotEvery) == 0 || n == 1)
        plotState(u, v, p, t);
    end

    % Export PNG frames occasionally (disable by setting savePngEvery=0)
    if p.savePngEvery > 0 && mod(n, p.savePngEvery) == 0
        exportFrame(u, v, p, t, outDir);
    end
end

% Always export final state
exportFrame(u, v, p, t, outDir);

% Optional: save final vectors (useful for restart / inspection)
save(fullfile(outDir, "final_state.mat"), "u", "v", "t", "grid", "p");

disp("Results saved in:");
disp(outDir);
