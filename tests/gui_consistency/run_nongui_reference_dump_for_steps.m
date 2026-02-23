function refFile = run_nongui_reference_dump_for_steps(rootDir, p, nStepsTarget)
%RUN_NONGUI_REFERENCE_DUMP_FOR_STEPS  Non-GUI reference dump for a fixed step count.
%
% This is a test helper (kept under tests/). It does not modify the main program.
%
% Inputs:
%   rootDir      : project root directory
%   p            : finalized parameter struct (already seeded outside)
%   nStepsTarget : number of steps to run
%
% Output:
%   refFile      : full path to saved MAT-file in results/consistency/

addpath(genpath(fullfile(rootDir, "src")));

% Deterministic seed (same as non-GUI script)
seed = p.seed;
rng(seed, "twister");

% Build model
grid  = buildGrid(p);
model = GrayScottModel(p, grid);

% Snapshot at midpoint
kSnap = round(nStepsTarget/2);

settings.nSteps = nStepsTarget;
settings.kSnap  = kSnap;

snap = struct();
snap.u0 = model.u;
snap.v0 = model.v;

for n = 1:nStepsTarget
    info = model.step();
    if isfield(info,"hasNaNInf") && info.hasNaNInf
        error("NaN/Inf detected during non-GUI reference at step %d.", model.n);
    end
    if n == kSnap
        snap.us = model.u;
        snap.vs = model.v;
        snap.kSnap = kSnap;
    end
end

snap.uf = model.u;
snap.vf = model.v;
snap.nSteps = nStepsTarget;

% Save
outDir = fullfile(rootDir, "results", "consistency");
if ~exist(outDir, "dir"), mkdir(outDir); end

timestamp = datestr(now, "yyyy-mm-dd_HH-MM-SS");
refFile = fullfile(outDir, "nongui_reference_" + string(timestamp) + ".mat");

save(refFile, "p", "seed", "settings", "snap", "-v7.3");
end