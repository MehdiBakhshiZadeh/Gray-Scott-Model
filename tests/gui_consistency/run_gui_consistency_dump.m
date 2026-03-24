function out = run_gui_consistency_dump()
%RUN_GUI_CONSISTENCY_DUMP
% Launch the GUI, use its CURRENT shipped setup, advance the GUI-owned model
% for the exact number of required steps, and save a GUI consistency dump
% MAT file for compare_gui_nongui.
%
% This helper:
% - uses the GUI startup defaults exactly as they are now
% - does NOT rely on private app.TEnd access
% - advances the model for exactly nSteps = p.T/p.dt
% - captures initial, midpoint snapshot, and final state
%
% Saved file format matches compare_gui_nongui:
%   p, seed, settings, snap
%
% Usage:
%   ref = run_nongui_reference_dump();
%   gui = run_gui_consistency_dump();
%   T   = compare_gui_nongui(ref.refFile, gui.guiFile);

% ---------------- Locate project root and add paths ----------------
thisFile = mfilename("fullpath");
rootDir  = fileparts(fileparts(fileparts(thisFile)));

addpath(genpath(fullfile(rootDir, "src")));
addpath(genpath(fullfile(rootDir, "tests")));
addpath(fullfile(rootDir, "scripts"));

% ---------------- Close any old app window ----------------
old = findall(0, "Type", "figure", "Tag", "GrayScottControllerUI");
if ~isempty(old)
    delete(old);
end

% ---------------- Launch GUI app ----------------
app = GrayScottController;
drawnow;

S0 = localAppStruct(app);
assert(isfield(S0, "model") && ~isempty(S0.model), ...
    "GUI helper: app.model is missing after startup.");

model = S0.model;
p     = model.p;
seed  = p.seed;

% GUI consistency logic assumes exact integer number of steps
Nt_exact = p.T / p.dt;
assert(abs(Nt_exact - round(Nt_exact)) < 1e-12, ...
    "GUI helper: p.T must be divisible by p.dt.");

settings.nSteps = round(Nt_exact);
settings.kSnap  = round(settings.nSteps / 2);

fprintf("GUI consistency helper started.\n");
fprintf("Current GUI setup:\n");
fprintf("  preset          : %s\n", string(p.preset));
fprintf("  diffusionMode   : %s\n", string(p.diffusionMode));
fprintf("  BC(left)        : %s\n", string(p.BC.left.type));
fprintf("  grid            : [%dx%d]\n", p.Nx, p.Ny);
fprintf("  dt              : %.6g\n", p.dt);
fprintf("  T               : %.6g\n", p.T);
fprintf("  nSteps target   : %d\n", settings.nSteps);
fprintf("  kSnap target    : %d\n", settings.kSnap);

% ---------------- Capture initial state from GUI-owned model ----------------
snap = struct();
snap.u0 = model.u;
snap.v0 = model.v;

% ---------------- Advance EXACTLY nSteps ----------------
for k = 1:settings.nSteps
    info = model.step();

    if isfield(info, "hasNaNInf") && info.hasNaNInf
        error("GUI helper: NaN/Inf encountered at step %d.", k);
    end

    if k == settings.kSnap
        snap.us   = model.u;
        snap.vs   = model.v;
        snap.kSnap = settings.kSnap;
    end
end

% ---------------- Capture final state ----------------
snap.uf = model.u;
snap.vf = model.v;
snap.n  = model.n;
snap.nSteps = settings.nSteps;

assert(snap.n == settings.nSteps, ...
    "GUI helper: final step mismatch (got n=%d, expected %d).", ...
    snap.n, settings.nSteps);

assert(isfield(snap, "us") && ~isempty(snap.us), ...
    "GUI helper: snapshot state was not captured.");
assert(isfield(snap, "vs") && ~isempty(snap.vs), ...
    "GUI helper: snapshot state was not captured.");

% ---------------- Save MAT file ----------------
outDir = fullfile(rootDir, "results", "consistency");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

timestamp = datestr(now, "yyyy-mm-dd_HH-MM-SS");
guiFile = fullfile(outDir, "gui_dump_" + string(timestamp) + ".mat");

save(guiFile, "p", "seed", "settings", "snap", "-v7.3");

fprintf("GUI consistency dump saved:\n  %s\n", guiFile);

% ---------------- Optional cleanup ----------------
try
    if isvalid(app.GrayScottContollerUIFigure)
        delete(app.GrayScottContollerUIFigure);
    end
catch
end

out.guiFile  = guiFile;
out.settings = settings;

end


% ========================= local helpers =========================

function S = localAppStruct(app)
warnState = warning("off", "MATLAB:structOnObject");
cleanupWarn = onCleanup(@() warning(warnState)); %#ok<NASGU>
S = struct(app);
end