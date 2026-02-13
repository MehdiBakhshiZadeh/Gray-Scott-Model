function benchmark_solvers()
%BENCHMARK_SOLVERS  Benchmark protocol for solver comparisons.
%
% Compares three diffusion implementations:
%   - "sparse_matrix" : sparse Laplacian multiply (p.diffusionMode="matrix")
%   - "stencil"       : matrix-free stencil Laplacian (p.diffusionMode="stencil")
%   - "dense_matrix"  : dense/full Laplacian multiply (p.diffusionMode="full") [small grids only]
%
% Produces two timing modes:
%   1) compute_off : solver only (no visualization)
%   2) render_on   : solver + visualization (live display), update every plotEvery steps
%
% Timed region = ONLY the stepping loop for cfg.nsteps.
% Warm-up steps and setup (grid/operator/IC) are excluded from timing.
%
% Output:
%   results/benchmarks_<timestamp>/benchmark_results.mat
%   results/benchmarks_<timestamp>/benchmark_table.csv

root = fileparts(fileparts(mfilename('fullpath')));

% --- Make sure we use the CURRENT solver in src/, not the shadowed baseline copies ---
addpath(fullfile(root,'src'));
if exist(fullfile(root,'src','baseline'),'dir')
    % Remove baseline folder from path to avoid shadowing.
    % (If it wasn't on path, rmpath does nothing.)
    try
        rmpath(fullfile(root,'src','baseline'));
    catch
        % ignore
    end
end

% Grid policies
cfg.grids_all        = [32 32; 64 64; 128 128; 256 256; 512 512];
cfg.grids_dense_only = [32 32; 64 64; 128 128];

cfg.warmup  = 50;
cfg.nsteps  = 500;
cfg.repeats = 7;

cfg.seed = 1;       % fixed seed for reproducibility
cfg.plotEvery = 10; % fixed: update display every 10 steps (render_on mode)

% Solver definitions:
% label = what appears in the table
% mode  = what is passed into p.diffusionMode
solvers = struct( ...
    "label", {"sparse_matrix", "stencil", "dense_matrix"}, ...
    "mode",  {"matrix",       "stencil", "full"} ...
);

modes = ["compute_off","render_on"];

timestamp = datestr(now, "yyyy-mm-dd_HH-MM-SS");
outDir = fullfile(root, "results", "benchmarks_" + string(timestamp));
if ~exist(outDir, 'dir'); mkdir(outDir); end

% --- Canonical result row template (all fields, always present) ---
rowTemplate = struct( ...
    "solver","", ...
    "mode","", ...
    "Nx",NaN, ...
    "Ny",NaN, ...
    "warmup",NaN, ...
    "nsteps",NaN, ...
    "repeats",NaN, ...
    "dt_used",NaN, ...
    "plotEvery",NaN, ...
    "total_sec_median",NaN, ...
    "total_sec_min",NaN, ...
    "total_sec_iqr",NaN, ...
    "sec_per_step_median",NaN, ...
    "mem_bytes_est",NaN, ...
    "fps_median",NaN ...
);

rows = rowTemplate([]);   % empty struct WITH the right fields
row_i = 0;

% ===============================
% Main benchmark loops
% ===============================
for mi = 1:numel(modes)
    mode = modes(mi);

    for gi = 1:size(cfg.grids_all,1)
        Nx = cfg.grids_all(gi,1);
        Ny = cfg.grids_all(gi,2);

        % Compute dt_used once per grid (shared across solvers)
        dt_used = compute_dt_used(Nx, Ny, cfg);

        for si = 1:numel(solvers)
            solverLabel = string(solvers(si).label);
            solverMode  = string(solvers(si).mode);  % "matrix" | "stencil" | "full"

            % Enforce dense-only grids for "full"
            if solverMode == "full" && ~grid_is_allowed_dense(Nx, Ny, cfg)
                % Skip generating rows for dense outside allowed grids.
                continue;
            end

            % Repeat timings
            times = zeros(cfg.repeats,1);
            fps_vals = NaN(cfg.repeats,1);

            for r = 1:cfg.repeats
                switch mode
                    case "compute_off"
                        times(r) = run_one_compute_timed(Nx, Ny, cfg, solverMode, dt_used);
                    case "render_on"
                        [times(r), fps_vals(r)] = run_one_render_timed(Nx, Ny, cfg, solverMode, dt_used);
                    otherwise
                        error("Unknown mode: %s", mode);
                end
            end

            % Summary stats
            t_med = median(times);
            t_min = min(times);
            t_iqr = iqr(times);

            sec_per_step_med = t_med / cfg.nsteps;

            % Memory estimate depends on solverMode
            mem_bytes = estimate_memory_bytes(Nx, Ny, solverMode);

            % Record plotEvery for reporting
            pe = NaN;
            if mode == "render_on"
                pe = cfg.plotEvery;
            else
                pe = 0; % compute_off
            end

            row_i = row_i + 1;
            row = rowTemplate;

            row.solver = solverLabel;
            row.mode   = mode;
            row.Nx = Nx; row.Ny = Ny;
            row.warmup  = cfg.warmup;
            row.nsteps  = cfg.nsteps;
            row.repeats = cfg.repeats;
            row.dt_used = dt_used;
            row.plotEvery = pe;

            row.total_sec_median    = t_med;
            row.total_sec_min       = t_min;
            row.total_sec_iqr       = t_iqr;
            row.sec_per_step_median = sec_per_step_med;
            row.mem_bytes_est       = mem_bytes;

            if mode == "render_on"
                row.fps_median = median(fps_vals);
            else
                row.fps_median = NaN;
            end

            rows(row_i) = row;
        end
    end
end

% Save MAT
save(fullfile(outDir, "benchmark_results.mat"), "cfg", "rows");

% Save CSV
T = struct2table(rows);

% Stable column order for easier report tables
colOrder = { ...
    'solver','mode','Nx','Ny','warmup','nsteps','repeats','dt_used','plotEvery', ...
    'total_sec_median','total_sec_min','total_sec_iqr', ...
    'sec_per_step_median','mem_bytes_est','fps_median' ...
};
T = T(:, colOrder);

writetable(T, fullfile(outDir, "benchmark_table.csv"));

fprintf("Benchmark results saved in:\n%s\n", outDir);
end

% =========================================================
% Helpers
% =========================================================

function ok = grid_is_allowed_dense(Nx, Ny, cfg)
% Dense/full matrix Laplacian allowed only on cfg.grids_dense_only
ok = any(cfg.grids_dense_only(:,1) == Nx & cfg.grids_dense_only(:,2) == Ny);
end

function t_total = run_one_compute_timed(Nx, Ny, cfg, solverMode, dt_used)
% Compute mode: render OFF (solver only).
% Warm-up is executed but NOT timed. Setup is NOT timed.
% Timed region is only cfg.nsteps calls to stepEuler.

% Build params
p = defaultParams();
p.seed = cfg.seed;
p.Nx = Nx;
p.Ny = Ny;

if ~isfield(p,'bc') || isempty(p.bc)
    p.bc = "periodic";
end

p.plotEvery = 0;
p.savePngEvery = 0;

% Setup (NOT timed)
rng(p.seed, "twister");
grid = buildGrid(p);
p.grid = grid;

p.dt = dt_used;
p.diffusionMode = solverMode; % "matrix" | "stencil" | "full"
[L, S] = build_diffusion_operator(p, grid, solverMode);

[U0, V0] = initialCondition(p);
u = U0(:);
v = V0(:);

% Warm-up (NOT timed, NO plotting)
t = 0.0;
for i = 1:cfg.warmup
    [u, v, info] = stepEuler(u, v, L, S, p, t);
    if isfield(info,'hasNaNInf') && info.hasNaNInf
        error("Benchmark aborted: NaN/Inf during warmup (Nx=%d, Ny=%d, solverMode=%s).", Nx, Ny, solverMode);
    end
    t = t + p.dt;
end

% Timed steps (ONLY this is measured)
t0 = tic;
for i = 1:cfg.nsteps
    [u, v, info] = stepEuler(u, v, L, S, p, t);
    if isfield(info,'hasNaNInf') && info.hasNaNInf
        error("Benchmark aborted: NaN/Inf during timing (Nx=%d, Ny=%d, solverMode=%s).", Nx, Ny, solverMode);
    end
    t = t + p.dt;
end
t_total = toc(t0);
end

function [t_total, fps] = run_one_render_timed(Nx, Ny, cfg, solverMode, dt_used)
% Render mode: solver + visualization ON.
% Warm-up is NOT timed and does NOT plot.
% Timed region includes solver + rendering (plotEvery=cfg.plotEvery).

p = defaultParams();
p.seed = cfg.seed;
p.Nx = Nx;
p.Ny = Ny;

if ~isfield(p,'bc') || isempty(p.bc)
    p.bc = "periodic";
end

p.plotEvery = cfg.plotEvery;
p.savePngEvery = 0;

% Setup (NOT timed)
rng(p.seed, "twister");
grid = buildGrid(p);
p.grid = grid;

p.dt = dt_used;
p.diffusionMode = solverMode; % "matrix" | "stencil" | "full"
[L, S] = build_diffusion_operator(p, grid, solverMode);

[U0, V0] = initialCondition(p);
u = U0(:);
v = V0(:);

% Warm-up (NOT timed, NO plotting)
t = 0.0;
for i = 1:cfg.warmup
    [u, v, info] = stepEuler(u, v, L, S, p, t);
    if isfield(info,'hasNaNInf') && info.hasNaNInf
        error("Render benchmark aborted: NaN/Inf warmup (Nx=%d, Ny=%d, solverMode=%s).", Nx, Ny, solverMode);
    end
    t = t + p.dt;
end

% Timed region (solver + render)
plotEvery = max(1, round(cfg.plotEvery));

t0 = tic;
for i = 1:cfg.nsteps
    [u, v, info] = stepEuler(u, v, L, S, p, t);
    if isfield(info,'hasNaNInf') && info.hasNaNInf
        error("Render benchmark aborted: NaN/Inf timing (Nx=%d, Ny=%d, solverMode=%s).", Nx, Ny, solverMode);
    end

    if mod(i, plotEvery) == 0
        safePlotState(u, v, p, t);
        drawnow;
    end

    t = t + p.dt;
end
t_total = toc(t0);

fps = cfg.nsteps / t_total;
end

function [L, S] = build_diffusion_operator(p, grid, solverMode)
% Build diffusion operator according to solver choice.
% Setup cost is excluded from timing.

solverMode = string(solverMode);

if solverMode == "stencil"
    S = buildStencil2D(p, grid);
    L = [];

elseif solverMode == "full"
    % Dense Laplacian: build sparse then convert to full
    L = full(buildLaplacian2D(p, grid));
    S = [];

else
    % Default: sparse matrix Laplacian ("matrix")
    L = buildLaplacian2D(p, grid);
    S = [];
end
end

function dt_used = compute_dt_used(Nx, Ny, cfg)
% Compute grid-dependent stable dt used for explicit Euler diffusion.
% Policy: dt = min(p.dt, dt_stable) where dt_stable uses a safety factor.

p = defaultParams();
p.seed = cfg.seed;
p.Nx = Nx;
p.Ny = Ny;

if ~isfield(p,'bc') || isempty(p.bc)
    p.bc = "periodic";
end

grid = buildGrid(p);
dt_used = stable_dt_explicit_euler(p, grid);
end

function dt = stable_dt_explicit_euler(p, grid)
% Stability-limited dt for explicit Euler diffusion on a 2D periodic grid.
% Conservative factor.

hx = grid.hx;
hy = grid.hy;

Dmax = max([p.Du, p.Dv]);
dt_diff = 1 / (2 * Dmax * (1/hx^2 + 1/hy^2));

safety = 0.8;
dt_stable = safety * dt_diff;

% Policy: dt = min(user_dt, stable_dt)
dt = min(p.dt, dt_stable);
end

function safePlotState(u, v, p, t)
% Try common plotState signatures without crashing the benchmark.

U = reshape(u, p.Ny, p.Nx);
V = reshape(v, p.Ny, p.Nx);

try
    plotState(U, V, p, t);
    return;
catch
end

try
    plotState(U, V, p);
    return;
catch
end

try
    plotState(u, v, p, t);
    return;
catch
end

try
    plotState(u, v, p);
    return;
catch
end

% If none match, do nothing.
end

function mem_bytes = estimate_memory_bytes(Nx, Ny, solverMode)
% Rough estimate:
% - all solvers store u and v as double vectors
% - sparse matrix solver also stores sparse Laplacian (~5N nnz)
% - dense matrix solver stores full Laplacian (N^2 doubles)
% - stencil solver stores only a tiny struct S (negligible)

N = Nx*Ny;

% u and v as double
bytes_uv = 2 * N * 8;

solverMode = string(solverMode);

if solverMode == "matrix"
    nnz_est = 5*N; % approx 5-point stencil
    % Very rough sparse storage estimate:
    % values (8 bytes each) + row indices (8 bytes each) + col ptr (8*(N+1))
    bytes_sparse = nnz_est*(8 + 8) + (N+1)*8;
    mem_bytes = bytes_uv + bytes_sparse;

elseif solverMode == "full"
    % Full matrix: N-by-N doubles
    bytes_full = (N * N) * 8;
    mem_bytes = bytes_uv + bytes_full;

else
    % "stencil"
    mem_bytes = bytes_uv + 1024;
end
end
