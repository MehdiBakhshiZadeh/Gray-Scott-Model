function out = run_nongui_reference_dump()
%RUN_NONGUI_REFERENCE_DUMP  Create a non-GUI reference snapshot for GUI consistency.
%
% This function runs the Gray–Scott model in non-GUI mode using the same
% parameter pipeline as scripts/run_model_class.m and saves:
%   - initial (u0,v0)
%   - snapshot at step kSnap (us,vs)
%   - final at step nSteps (uf,vf)
%
% Output file:
%   <projectRoot>/results/consistency/nongui_reference_<timestamp>.mat

% ---- Locate project root and add src to path (no cd side effects) ----
thisFile = mfilename("fullpath");
rootDir  = fileparts(fileparts(fileparts(thisFile)));
addpath(genpath(fullfile(rootDir, "src")));

% ---- Load and finalize parameters exactly like non-GUI runner ----
p = defaultParams();
p = finalizeParams(p);

% ---- Consistency test settings (chosen to match your defaults) ----
% DefaultParams: dt=0.5, T=500 => Nt = 1000 steps exactly
settings.nSteps = round(p.T / p.dt);
settings.kSnap  = round(settings.nSteps / 2);

% Safety: ensure T divisible by dt (same check style as your script)
Nt_exact = p.T / p.dt;
assert(abs(Nt_exact - round(Nt_exact)) < 1e-12, "T must be divisible by dt.");

% Disable plotting / PNG export in this test run (numerics only)
p.plotEvery    = 0;
p.savePngEvery = 0;

% ---- Deterministic seed (matches your script) ----
seed = p.seed;
rng(seed, "twister");

% ---- Build model ----
grid  = buildGrid(p);
model = GrayScottModel(p, grid);

% ---- Capture initial state ----
snap.u0 = model.u;
snap.v0 = model.v;

% ---- Step loop ----
for n = 1:settings.nSteps
    info = model.step();

    if isfield(info, "hasNaNInf") && info.hasNaNInf
        error("NaN/Inf detected at step %d (t=%.3g).", model.n, model.t);
    end

    if n == settings.kSnap
        snap.us = model.u;
        snap.vs = model.v;
        snap.kSnap = n;
    end
end

% ---- Final state ----
snap.uf = model.u;
snap.vf = model.v;
snap.nSteps = settings.nSteps;

% ---- Save ----
outDir = fullfile(rootDir, "results", "consistency");
if ~exist(outDir, "dir"), mkdir(outDir); end

timestamp = datestr(now, "yyyy-mm-dd_HH-MM-SS");
refFile = fullfile(outDir, "nongui_reference_" + string(timestamp) + ".mat");

save(refFile, "p", "seed", "settings", "snap", "-v7.3");

out.refFile  = refFile;
out.settings = settings;

fprintf("Non-GUI reference dump saved:\n  %s\n", refFile);
end