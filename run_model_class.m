clc;
close all force;

% Project root (works even if you run this from a different folder).
thisFile = mfilename("fullpath");
if strlength(thisFile) == 0
    % Fallback when mfilename is empty (e.g., Run Section / unusual execution contexts)
    thisFile = matlab.desktop.editor.getActiveFilename();
end
assert(strlength(thisFile) > 0, "Cannot locate run_model_class.m on disk.");

rootDir = fileparts(thisFile);

% Ensure src is visible (absolute path, no cd side-effect)
addpath(genpath(fullfile(rootDir, "src")));

% Load default parameters + seed
p = defaultParams();           % "pearson"| "baseline"
p.diffusionMode = "stencil";   % "matrix" | "stencil" | "full"
p = finalizeParams(p);

% Create a unique name for this run
timestamp = datestr(now, "yyyy-mm-dd_HH-MM-SS");
runName   = timestamp + "_" + p.caseName;

% Output directory (absolute, created safely)
resultsDir = fullfile(rootDir, "results");
if ~exist(resultsDir, "dir")
    [ok,msg] = mkdir(resultsDir);
    assert(ok, "Failed to create results directory: %s", msg);
end

outDir = fullfile(resultsDir, runName);
if ~exist(outDir, "dir")
    [ok,msg] = mkdir(outDir);
    assert(ok, "Failed to create run directory: %s", msg);
end

% Save metadata for reproducibility
meta.timestamp     = timestamp;
meta.parameters    = p;
meta.matlabVersion = version;
meta.platform      = computer;
save(fullfile(outDir, "meta.mat"), "meta");

% Choose which field to plot (default: v)
field = "v";
if isfield(p, "plotField")
    field = lower(string(p.plotField));
end

% Write a readable log file
logPath = fullfile(outDir, "log.txt");
fid = fopen(logPath, "w");
assert(fid ~= -1, "Could not open log file: %s", logPath);
c = onCleanup(@() fclose(fid)); %#ok<NASGU> % script errors out, is closed automatically

fprintf(fid, "Run name: %s\n", runName);
fprintf(fid, "Timestamp: %s\n", timestamp);
fprintf(fid, "MATLAB: %s\n", version);
fprintf(fid, "Platform: %s\n", computer);
fprintf(fid, "Grid: Nx=%d, Ny=%d\n", p.Nx, p.Ny);
fprintf(fid, "Time: dt=%g, T=%g\n", p.dt, p.T);
fprintf(fid, "Model: Du=%g, Dv=%g, F=%g, k=%g\n", p.Du, p.Dv, p.F, p.k);
if isfield(p,"diffusionMode")
    fprintf(fid, "Diffusion mode: %s\n", string(p.diffusionMode));
end
fprintf(fid, "Output: plotEvery=%d (0 disables)\n", p.plotEvery);
fprintf(fid, "Output: savePngEvery=%d (0 disables)\n", p.savePngEvery);
fprintf(fid, "Plot field: %s\n", field);
fprintf(fid, "Output directory: %s\n", outDir);

% Create model object (builds Laplacian, and initial condition)
grid  = buildGrid(p);
model = GrayScottModel(p, grid);

% Time loop settings
Nt_exact = p.T / p.dt;
assert(abs(Nt_exact - round(Nt_exact)) < 1e-12, "T must be divisible by dt.");
Nt = round(Nt_exact);
fprintf(fid, "Steps: Nt=%d\n", Nt);

% Only create a figure when plotting is enabled
doAnyPlot = (p.plotEvery > 0);

hMainFig = [];
hMainAx  = [];
hMainIm  = [];
hMainCb  = [];

if doAnyPlot
    hMainFig = figure( ...
        "MenuBar", "none", ...
        "ToolBar", "none", ...
        "NumberTitle", "off", ...
        "Name", "Gray-Scott Simulation", ...
        "WindowStyle", "normal");

    hMainAx = axes("Parent", hMainFig);
end

% ---- Time stepping loop ----
for n = 1:Nt
    info = model.step();

    if info.hasNaNInf
        error("NaN/Inf detected at step %d (t=%.3f). Stopping.", model.n, model.t);
    end

    % ---- Plot cadence (controls ONLY visualization) ----
    doPlot = doAnyPlot && (mod(model.n, p.plotEvery) == 0 || model.n == 1);

    if doPlot
        [U, V] = model.getFields2D();

        switch field
            case "u"
                data2D = U;
                fieldName = "u";
            case "v"
                data2D = V;
                fieldName = "v";
            otherwise
                error("Unknown plotField='%s' (use 'u' or 'v')", field);
        end

        if isempty(hMainIm) || ~isgraphics(hMainIm)
            hMainIm = imagesc(hMainAx, data2D);
            axis(hMainAx, "image");
            hMainCb = colorbar(hMainAx); %#ok<NASGU>
        else
            set(hMainIm, "CData", data2D);
        end

        title(hMainAx, sprintf("%s-field | t = %.2f | F=%.3f k=%.3f Du=%.2g Dv=%.2g", ...
            fieldName, model.t, p.F, p.k, p.Du, p.Dv));

        drawnow;
    end

    % ---- Save cadence ----
    doSave = (p.savePngEvery > 0) && (mod(model.n, p.savePngEvery) == 0);
    if doSave
        exportFrame(model.u, model.v, p, model.t, outDir);

        % Force refresh on the visible simulation window after export
        if doAnyPlot && ~isempty(hMainFig) && isgraphics(hMainFig)
            drawnow;
        end
    end
end

% Always export final state once
exportFrame(model.u, model.v, p, model.t, outDir);

% Save final state for restart / inspection
finalState = struct();
finalState.u = model.u;
finalState.v = model.v;
finalState.t = model.t;
finalState.n = model.n;
finalState.p = p;

save(fullfile(outDir, "final_state.mat"), "finalState");

disp("Results saved in:");
disp(outDir);