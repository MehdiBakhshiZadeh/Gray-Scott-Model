function r = test_pattern_sanity_pearson()
%TEST_PATTERN_SANITY_PEARSON  Pearson-style qualitative pattern sanity test.
%
% This test reproduces four canonical Pearson Gray-Scott cases
% (alpha, kappa, epsilon, theta) on a 256x256 periodic grid using
% dt = 1 and 200,000 time steps.
%
% The test is qualitative: it saves snapshots of the U and V fields
% at selected times and a final summary panel for visual inspection.

% Locate this test folder and project root
thisFile = mfilename("fullpath");
thisDir  = fileparts(thisFile);   % .../tests/test_case_09_pearson_pattern_sanity
testsDir = fileparts(thisDir);    % .../tests
rootDir  = fileparts(testsDir);   % project root

% Make this test runnable independently
if isempty(which("defaultParams")) || isempty(which("makeOperator"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(testsDir));
end

% Local output folder for this test
outDir = fullfile(thisDir, "outputs");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

% --- Pearson baseline setup ---
p0 = defaultParams("pearson");
p0.plotEvery = 0;
p0.savePngEvery = 0;

% Pearson-style long run
p0.dt = 1.0;
Nt_full = 200000;
p0.T = Nt_full * p0.dt;

% Snapshot times
tSnaps = [50000 100000 150000 200000];

% Pearson cases
cases = pearson_cases_alpha_kapa_epsilon_theta();

allFiles = strings(0,1);
finalU = cell(numel(cases),1);

% Global metadata
metaAll = struct();
metaAll.test = "pattern_sanity_pearson";
metaAll.grid = [p0.Nx p0.Ny];
metaAll.domain = [p0.Lx p0.Ly];
metaAll.Du = p0.Du;
metaAll.Dv = p0.Dv;
metaAll.BC = p0.BC;
metaAll.dt = p0.dt;
metaAll.Nt = Nt_full;
metaAll.T = p0.T;
metaAll.snapshotTimes = tSnaps;
metaAll.cases = cases;

fprintf("Running Pearson pattern sanity test\n");
fprintf("Grid: [%dx%d], dt=%.3g, Nt=%d, T=%.3g\n", p0.Nx, p0.Ny, p0.dt, Nt_full, p0.T);

% --- Run each Pearson case once ---
for icase = 1:numel(cases)
    p = p0;
    p.F = cases(icase).F;
    p.k = cases(icase).k;
    p.seed = 2000 + icase;

    fprintf("Case %d/%d: %s (F=%.6g, k=%.6g)\n", ...
        icase, numel(cases), cases(icase).tag, p.F, p.k);

    g  = buildGrid(p);
    op = makeOperator(p, g);

    [U0, V0] = initialCondition(p);
    u = U0(:);
    v = V0(:);

    snapSteps = round(tSnaps / p.dt);
    snapIdx = 1;

    for n = 1:Nt_full
        t = n * p.dt;
        [u, v, info] = eulerStep(u, v, op, p, t);

        if info.hasNaNInf
            error("NaN/Inf in Pearson case %s at step %d.", cases(icase).tag, n);
        end

        if snapIdx <= numel(snapSteps) && n == snapSteps(snapIdx)
            U = reshape(u, p.Ny, p.Nx);
            V = reshape(v, p.Ny, p.Nx);

            % --- Save U snapshot ---
            figU = figure("Visible","off", "Color","w");
            imagesc(U);
            axis image;
            colorbar;
            title(sprintf("Pearson %s: U, F=%.6g, k=%.6g, t=%d", ...
                cases(icase).tag, p.F, p.k, t));

            fnameU = fullfile(outDir, sprintf("%s_t%06d_U.png", cases(icase).tag, t));
            saveFigWithMeta(figU, fnameU, makeCaseMeta(p, cases(icase), t, "U"));
            close(figU);

            % --- Save V snapshot ---
            figV = figure("Visible","off", "Color","w");
            imagesc(V);
            axis image;
            colorbar;
            title(sprintf("Pearson %s: V, F=%.6g, k=%.6g, t=%d", ...
                cases(icase).tag, p.F, p.k, t));

            fnameV = fullfile(outDir, sprintf("%s_t%06d_V.png", cases(icase).tag, t));
            saveFigWithMeta(figV, fnameV, makeCaseMeta(p, cases(icase), t, "V"));
            close(figV);

            fprintf("  saved snapshots at t=%d\n", t);

            allFiles(end+1:end+2,1) = [string(fnameU); string(fnameV)];
            snapIdx = snapIdx + 1;
        end
    end

    % Store final U for summary
    finalU{icase} = reshape(u, p.Ny, p.Nx);
end

% --- Summary panel (final U only) ---
fig = figure("Visible","off", "Color","w");
nCases = numel(cases);
for icase = 1:nCases
    subplot(2,2,icase);
    imagesc(finalU{icase});
    axis image off;
    title(sprintf("%s (%s)", cases(icase).tag, cases(icase).label));
end
colormap(parula);
sgtitle("Pearson pattern sanity (U field, final time)");

summaryName = fullfile(outDir, "pearson_patterns_summary_U_final.png");
saveFigWithMeta(fig, summaryName, metaAll);
close(fig);

allFiles(end+1,1) = string(summaryName);

fprintf("Summary figure saved:\n  %s\n", summaryName);

% --- Result struct ---
r = makeResult("test_pattern_sanity_pearson", p0);
r.metrics.numCases = nCases;
r.metrics.numSnapshotsPerCase = numel(tSnaps);
r.metrics.totalImages = numel(allFiles);
r.metrics.Nt = Nt_full;
r.metrics.snapshotTimes = tSnaps;
r.pass = true;
r.figFiles = allFiles;
r.notes = sprintf("Pearson sanity run (full 200k steps) for cases: %s", ...
    strjoin([cases.tag], ", "));
end

function meta = makeCaseMeta(p, caseDef, t, fieldName)
%MAKECASEMETA  Metadata helper for saved Pearson figures.
meta = struct();
meta.test = "pattern_sanity_pearson";
meta.caseTag = caseDef.tag;
meta.caseLabel = caseDef.label;
meta.field = string(fieldName);
meta.time = t;
meta.grid = [p.Nx p.Ny];
meta.domain = [p.Lx p.Ly];
meta.Du = p.Du;
meta.Dv = p.Dv;
meta.F = p.F;
meta.k = p.k;
meta.dt = p.dt;
meta.seed = p.seed;
meta.BC = p.BC;
end