function out = run_nongui_reference_dump()
%RUN_NONGUI_REFERENCE_DUMP  Create a non-GUI reference snapshot for GUI consistency.
%
% This function runs the Gray--Scott model in non-GUI mode using the same
% parameter pipeline as scripts/run_model_class.m and saves:
%   - initial (u0,v0)
%   - snapshot at step kSnap (us,vs)
%   - final at step nSteps (uf,vf)
%
% Output file:
%   tests/test_case_06_gui_nongui_consistency/outputs/nongui_reference_<timestamp>.mat

% Locate this test folder and project root
thisFile = mfilename("fullpath");
thisDir  = fileparts(thisFile);   % .../tests/test_case_06_gui_nongui_consistency
testsDir = fileparts(thisDir);    % .../tests
rootDir  = fileparts(testsDir);   % project root

% Make this function runnable independently
if isempty(which("defaultParams")) || isempty(which("GrayScottModel"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(testsDir));
end

% Local output folder for this test
outDir = fullfile(thisDir, "outputs");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

% ---- Load and finalize parameters exactly like non-GUI runner ----
p = defaultParams("pearson");
p = finalizeParams(p);

% ---- Consistency test settings (chosen to match your defaults) ----
% DefaultParams: dt=0.5, T=500 => Nt = 1000 steps exactly
settings.nSteps = round(p.T / p.dt);
settings.kSnap  = round(settings.nSteps / 2);

% Safety: ensure T divisible by dt
Nt_exact = p.T / p.dt;
assert(abs(Nt_exact - round(Nt_exact)) < 1e-12, "T must be divisible by dt.");

% Disable plotting / PNG export in this test run (numerics only)
p.plotEvery    = 0;
p.savePngEvery = 0;

% ---- Deterministic seed ----
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
timestamp = datestr(now, "yyyy-mm-dd_HH-MM-SS");
refFile = fullfile(outDir, "nongui_reference_" + string(timestamp) + ".mat");

save(refFile, "p", "seed", "settings", "snap", "-v7.3");

out.refFile  = refFile;
out.settings = settings;

fprintf("Non-GUI reference dump saved:\n  %s\n", refFile);
end