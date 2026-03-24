function r = test_matrix_stencil_endtoend_equivalence()
%TEST_MATRIX_STENCIL_ENDTOEND_EQUIVALENCE
% End-to-end comparison using GrayScottModel class:
% sparse-matrix diffusion vs stencil diffusion.
%
% Note:
% Diffusion operators use periodic connectivity. Physical boundary conditions
% are enforced after each time step and are not varied in this equivalence test.
%
% This version is report-friendly:
% - runnable independently
% - prints quantitative comparison metrics
% - stores metrics and setup in the returned result struct

% Locate project folders relative to this file
thisDir = fileparts(mfilename("fullpath"));   % .../tests
rootDir = fileparts(thisDir);                 % project root

% Make this test runnable independently
if isempty(which("defaultParams"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(fullfile(rootDir, "tests")));
end

testName = "test_matrix_stencil_endtoend_equivalence";
r = makeResult(testName, struct());
r.pass  = false;
r.notes = "";

try
    % --- Fixed run setup ---
    p = defaultParams();

    % Short run to keep it fast
    p.T  = 5;
    p.dt = 0.2;

    % Normalize parameters (solver mapping, BC defaults, etc.)
    p = finalizeParams(p);

    % Keep setup info for reporting
    Nx = p.Nx;
    Ny = p.Ny;
    T  = p.T;
    dt = p.dt;

    fprintf('Running end-to-end matrix vs stencil equivalence test\n');
    fprintf('Grid: [%dx%d], dt=%.3g, T=%.3g\n', Nx, Ny, dt, T);
    fprintf('Parameters: Du=%.3e, Dv=%.3e, F=%.3e, k=%.3e\n', ...
        p.Du, p.Dv, p.F, p.k);

    % If IC uses RNG, keep it reproducible
    rng(1, "twister");

    % Build grid once (same for both variants)
    grid = buildGrid(p);

    % --- Baseline (matrix) ---
    pA = p;
    pA.diffusionMode = "matrix";
    modelA = GrayScottModel(pA, grid);

    while modelA.t < pA.T
        info = modelA.step();
        if isfield(info,"hasNaNInf") && info.hasNaNInf
            error("Baseline run blew up (NaN/Inf).");
        end
    end
    [UA, VA] = modelA.getFields2D();

    % --- Stencil ---
    rng(1, "twister"); % reset RNG so IC matches exactly
    pB = p;
    pB.diffusionMode = "stencil";
    modelB = GrayScottModel(pB, grid);

    while modelB.t < pB.T
        info = modelB.step();
        if isfield(info,"hasNaNInf") && info.hasNaNInf
            error("Stencil run blew up (NaN/Inf).");
        end
    end
    [UB, VB] = modelB.getFields2D();

    % --- Compare final states ---
    uA = UA(:); vA = VA(:);
    uB = UB(:); vB = VB(:);

    du = uA - uB;
    dv = vA - vB;

    rel_u = norm(du) / max(1, norm(uA));
    rel_v = norm(dv) / max(1, norm(vA));
    max_u = max(abs(du));
    max_v = max(abs(dv));

    relTol = 1e-10;
    absTol = 1e-9;

    fprintf('rel_u = %.3e\n', rel_u);
    fprintf('rel_v = %.3e\n', rel_v);
    fprintf('max_u = %.3e\n', max_u);
    fprintf('max_v = %.3e\n', max_v);
    fprintf('Tolerances: relTol=%.1e, absTol=%.1e\n', relTol, absTol);

    assert(rel_u < relTol && rel_v < relTol, ...
        "rel mismatch: rel_u=%.3e, rel_v=%.3e (tol=%.1e)", rel_u, rel_v, relTol);
    assert(max_u < absTol && max_v < absTol, ...
        "max mismatch: max_u=%.3e, max_v=%.3e (tol=%.1e)", max_u, max_v, absTol);

    fprintf('PASS: matrix and stencil modes agree within tolerance.\n');

    r.pass = true;
    r.grid = [Nx, Ny];
    r.dt   = dt;
    r.T    = T;
    r.metrics.rel_u = rel_u;
    r.metrics.rel_v = rel_v;
    r.metrics.max_u = max_u;
    r.metrics.max_v = max_v;
    r.thresholds.relTol = relTol;
    r.thresholds.absTol = absTol;
    r.notes = sprintf([ ...
        'PASS; rel_u=%.3e; rel_v=%.3e; max_u=%.3e; max_v=%.3e'], ...
        rel_u, rel_v, max_u, max_v);

catch ME
    r.pass = false;
    r.notes = "Crashed: " + string(ME.message);
end
end