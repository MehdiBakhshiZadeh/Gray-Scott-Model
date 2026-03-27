function T = compare_gui_nongui(refMatFile, guiMatFile)
%COMPARE_GUI_NONGUI  Compare GUI and non-GUI states and generate artifacts.
%
% Usage:
%   T = compare_gui_nongui(refMatFile, guiMatFile);
%   T = compare_gui_nongui();
%
% Inputs:
%   refMatFile : path to nongui_reference_*.mat
%   guiMatFile : path to gui_dump_*.mat
%
% If inputs are omitted, the function looks for them in:
%   tests/test_case_06_gui_nongui_consistency/inputs/
%
% Output:
%   T : table with error metrics and pass/fail vs tolerances.
%
% Files written:
%   - outputs/gui_consistency.png
%   - outputs/gui_consistency_table.csv

% -------- Locate this test folder and project root --------
thisFile = mfilename("fullpath");
thisDir  = fileparts(thisFile);   % .../tests/test_case_06_gui_nongui_consistency
testsDir = fileparts(thisDir);    % .../tests
rootDir  = fileparts(testsDir);   % project root

% -------- Make runnable independently --------
if isempty(which("saveFigWithMeta"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(testsDir));
end

% -------- Local input/output folders --------
inDir  = fullfile(thisDir, "inputs");
outDir = fullfile(thisDir, "outputs");

if ~exist(inDir, "dir")
    mkdir(inDir);
end
if ~exist(outDir, "dir")
    mkdir(outDir);
end

% -------- Resolve input files --------
if nargin < 1 || isempty(refMatFile)
    d = dir(fullfile(inDir, "nongui_reference_*.mat"));
    assert(~isempty(d), "No nongui_reference_*.mat file found in %s", inDir);
    [~, idx] = max([d.datenum]);
    refMatFile = fullfile(inDir, d(idx).name);
end

if nargin < 2 || isempty(guiMatFile)
    d = dir(fullfile(inDir, "gui_dump_*.mat"));
    assert(~isempty(d), "No gui_dump_*.mat file found in %s", inDir);
    [~, idx] = max([d.datenum]);
    guiMatFile = fullfile(inDir, d(idx).name);
end

% ---------------- Load ----------------
A = load(refMatFile);   % expects: p, seed, settings, snap
B = load(guiMatFile);   % expects: p, seed, settings, snap

% ---------------- Basic checks ----------------
assert(isfield(A,"snap") && isfield(B,"snap"), ...
    "Missing 'snap' struct in one of the MAT files.");
assert(isfield(A,"p") && isfield(B,"p"), ...
    "Missing 'p' parameter struct in one of the MAT files.");
assert(isfield(A,"seed") && isfield(B,"seed"), ...
    "Missing 'seed' in one of the MAT files.");

assert(A.seed == B.seed, "Seed mismatch (ref=%d, gui=%d).", A.seed, B.seed);

% Compare key parameters that must match for a strict test
mustMatch = ["Nx","Ny","dt","Du","Dv","F","k"];
for f = mustMatch
    assert(isfield(A.p,f) && isfield(B.p,f), "Missing parameter p.%s in one file.", f);
    if ~isequal(A.p.(f), B.p.(f))
        error("Parameter mismatch: p.%s (ref=%g, gui=%g).", f, A.p.(f), B.p.(f));
    end
end

% ---------------- Tolerances ----------------
tol.relL2  = 1e-12;
tol.maxAbs = 1e-12;
epsNorm = 1e-30;

metrics = @(xRef,xGui) struct( ...
    "relL2",  norm(xGui - xRef) / (norm(xRef) + epsNorm), ...
    "maxAbs", max(abs(xGui - xRef)) );

% ---------------- Snapshot comparison ----------------
haveSnapGUI = isfield(B.snap,"us") && ~isempty(B.snap.us) && ...
              isfield(B.snap,"vs") && ~isempty(B.snap.vs);

m_us = [];
m_vs = [];
if haveSnapGUI
    m_us = metrics(A.snap.us, B.snap.us);
    m_vs = metrics(A.snap.vs, B.snap.vs);
else
    warning("GUI snapshot (us/vs) is missing. Run past kSnap and export again.");
end

% ---------------- Final comparison ----------------
haveFinalGUI = isfield(B.snap,"uf") && ~isempty(B.snap.uf) && ...
               isfield(B.snap,"vf") && ~isempty(B.snap.vf);

haveStepCounter = isfield(B.snap,"n") && ~isempty(B.snap.n);

refNSteps = NaN;
if isfield(A.snap,"nSteps"), refNSteps = A.snap.nSteps; end

guiNSteps = NaN;
if isfield(B,"settings") && isfield(B.settings,"nSteps")
    guiNSteps = B.settings.nSteps;
elseif isfield(B.snap,"nSteps")
    guiNSteps = B.snap.nSteps;
end

sameFinal = haveFinalGUI && haveStepCounter && ...
            ((~isnan(refNSteps) && (B.snap.n == refNSteps)) || ...
             (~isnan(guiNSteps) && (B.snap.n == guiNSteps)));

m_uf = [];
m_vf = [];
if sameFinal
    m_uf = metrics(A.snap.uf, B.snap.uf);
    m_vf = metrics(A.snap.vf, B.snap.vf);
end

% ---------------- Results table ----------------
rows = {};

if haveSnapGUI
    rows(end+1,:) = {"snapshot","u", m_us.relL2, m_us.maxAbs, tol.relL2, tol.maxAbs, ...
        (m_us.relL2 <= tol.relL2) && (m_us.maxAbs <= tol.maxAbs)}; %#ok<AGROW>
    rows(end+1,:) = {"snapshot","v", m_vs.relL2, m_vs.maxAbs, tol.relL2, tol.maxAbs, ...
        (m_vs.relL2 <= tol.relL2) && (m_vs.maxAbs <= tol.maxAbs)}; %#ok<AGROW>
end

if sameFinal
    rows(end+1,:) = {"final","u", m_uf.relL2, m_uf.maxAbs, tol.relL2, tol.maxAbs, ...
        (m_uf.relL2 <= tol.relL2) && (m_uf.maxAbs <= tol.maxAbs)}; %#ok<AGROW>
    rows(end+1,:) = {"final","v", m_vf.relL2, m_vf.maxAbs, tol.relL2, tol.maxAbs, ...
        (m_vf.relL2 <= tol.relL2) && (m_vf.maxAbs <= tol.maxAbs)}; %#ok<AGROW>
else
    rows(end+1,:) = {"final","u", NaN, NaN, tol.relL2, tol.maxAbs, false}; %#ok<AGROW>
    rows(end+1,:) = {"final","v", NaN, NaN, tol.relL2, tol.maxAbs, false}; %#ok<AGROW>
end

T = cell2table(rows, "VariableNames", ...
    ["stage","field","relL2","maxAbs","tol_relL2","tol_maxAbs","pass"]);

disp(T);

% ---------------- Figure (one) ----------------
if haveSnapGUI
    Nx = A.p.Nx;
    Ny = A.p.Ny;
    dU = reshape(B.snap.us - A.snap.us, Ny, Nx);

    f = figure("Visible","off", "Color","w");
    ax = axes(f);

    imagesc(ax, dU);
    axis(ax, "image");
    colorbar(ax);

    title(ax, sprintf("GUI - nonGUI difference (snapshot u), relL2=%.3e, max=%.3e", ...
        m_us.relL2, m_us.maxAbs));
    xlabel(ax, "x-index");
    ylabel(ax, "y-index");

    set(ax, ...
        "Color", "w", ...
        "XColor", "k", ...
        "YColor", "k", ...
        "FontSize", 16, ...
        "LineWidth", 1.2, ...
        "Box", "on");

    set(get(ax, "Title"),  "Color", "k", "FontSize", 18);
    set(get(ax, "XLabel"), "Color", "k", "FontSize", 16);
    set(get(ax, "YLabel"), "Color", "k", "FontSize", 16);

    % Better contrast for a zero-difference field
    caxis(ax, [-1 1]);

    figPath = fullfile(outDir, "gui_consistency.png");

    meta = struct();
    meta.test = "gui_consistency";
    meta.stage = "snapshot_u_difference";
    meta.relL2 = m_us.relL2;
    meta.maxAbs = m_us.maxAbs;
    meta.refFile = refMatFile;
    meta.guiFile = guiMatFile;

    saveFigWithMeta(f, figPath, meta);
    close(f);
end

% ---------------- Table export ----------------
csvPath = fullfile(outDir, "gui_consistency_table.csv");
writetable(T, csvPath);

% ---------------- Summary ----------------
fprintf("Consistency table saved:\n  %s\n", csvPath);

if haveSnapGUI
    fprintf("Consistency figure saved:\n  %s\n", figPath);
else
    fprintf("Consistency figure not written (snapshot missing).\n");
end

end