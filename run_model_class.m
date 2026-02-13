clear; clc; close all;

% Root of project = folder where this script lives
rootDir = fileparts(mfilename("fullpath"));

% Ensure src is visible (absolute path, no cd side-effect)
addpath(fullfile(rootDir, "src"));

% Load default parameters + seed
% Use: p = defaultParams("pearson");  % Pearson preset
p = defaultParams();
p.diffusionMode = "matrix";
modelA = GrayScottModel(p);
while modelA.t < p.T
    info = modelA.step();
end
[UA, VA] = modelA.getFields2D();

p.diffusionMode = "stencil";
modelB = GrayScottModel(p);
while modelB.t < p.T
    info = modelB.step();
end
[UB, VB] = modelB.getFields2D();
rng(p.seed, "twister");

% Create a unique name for this run
timestamp = datestr(now, "yyyy-mm-dd_HH-MM-SS");
runName   = timestamp + "_" + p.caseName;

% Output directory (absolute)
outDir = fullfile(rootDir, "results", runName);
if ~exist(outDir, "dir")
    mkdir(outDir);
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
fprintf(fid, "Output: plotEvery=%d (0 disables)\n", p.plotEvery);
fprintf(fid, "Output: savePngEvery=%d (0 disables)\n", p.savePngEvery);
fprintf(fid, "Plot field: %s\n", field);
fprintf(fid, "Output directory: %s\n", outDir);

% Create model object (builds grid, Laplacian, and initial condition)
model = GrayScottModel(p);
model.reset();

% Time loop settings
Nt_exact = p.T / p.dt;
assert(abs(Nt_exact - round(Nt_exact)) < 1e-12, "T must be divisible by dt.");
Nt = round(Nt_exact);
fprintf(fid, "Steps: Nt=%d\n", Nt);

% Only create a figure when plotting is enabled
doAnyPlot = (p.plotEvery > 0);
if doAnyPlot
    figure;
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
                imagesc(U); fieldName = "u";
            case "v"
                imagesc(V); fieldName = "v";
            otherwise
                error("Unknown plotField='%s' (use 'u' or 'v')", field);
        end

        axis image; colorbar;
        title(sprintf("%s-field | t = %.2f | F=%.3f k=%.3f Du=%.2g Dv=%.2g", ...
            fieldName, model.t, p.F, p.k, p.Du, p.Dv));
        drawnow;
    end

    % ---- Save cadence (matches run_grayscott.m) ----
    doSave = (p.savePngEvery > 0) && (mod(model.n, p.savePngEvery) == 0);
    if doSave
        exportFrame(model.u, model.v, p, model.t, outDir);
    end
end

% Always export final state once
exportFrame(model.u, model.v, p, model.t, outDir);

% Optional: save final state for restart / inspection
save(fullfile(outDir, "final_state.mat"), "model", "p");

disp("Results saved in:");
disp(outDir);
