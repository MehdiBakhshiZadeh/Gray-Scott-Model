function r = test_boundary_conditions()
%TEST_BOUNDARY_CONDITIONS Quantitative verification of Dirichlet/Neumann BC.
%
% This test verifies that boundary conditions are enforced strongly after
% explicit Euler updates and remain stable over time across all solver modes:
%   - "matrix"  : sparse Laplacian
%   - "stencil" : matrix-free 5-point stencil
%   - "full"    : dense Laplacian
%
% The test performs:
%   1) Dirichlet constant enforcement check (multiple steps).
%   2) Neumann derivative consistency check (multiple steps).
%   3) Mixed BC robustness check (many steps with quantitative error tracking).
%
% The same BC specification is applied to both u and v.

% Locate project folders relative to this file
thisDir = fileparts(mfilename("fullpath"));   % .../tests
rootDir = fileparts(thisDir);                 % project root

% Make this test runnable independently
if isempty(which("defaultParams"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(fullfile(rootDir, "tests")));
end

modes = ["matrix","stencil","full"];

% Test-control parameters
nStepsShort = 20;      % for local Dirichlet / Neumann checks
nStepsMixed = 400;     % for long mixed-BC robustness check
tolBC       = 1e-12;   % strict boundary tolerance
probeCol    = 8;       % interior probe column for nontrivial evolution check

metrics = struct();

for m = 1:numel(modes)
    mode = modes(m);

    % --- Parameters and grid size ---
    p = defaultParams();
    p.diffusionMode = mode;
    p.Nx = 64;
    p.Ny = 64;
    p.dt = 0.25;
    p.T  = nStepsMixed * p.dt;
    p.seed = 1;   % reproducible reset()
    p = finalizeParams(p);

    % --- Build model ---
    grid  = buildGrid(p);
    model = GrayScottModel(p, grid);
    p = model.p;  % keep BCcache-enabled version

    % --- Smooth custom initial fields for the short tests ---
    [XX,YY] = meshgrid(grid.x, grid.y);
    Ubase = 0.80 + 0.10*sin(2*pi*XX).*cos(2*pi*YY);
    Vbase = 0.20 + 0.05*cos(2*pi*XX).*sin(2*pi*YY);

    % ============================================================
    % Test 1: Dirichlet constant enforcement (multiple steps)
    % ============================================================
    p = model.p;
    p.BC.left.type = "dirichlet";
    p.BC.left.dirichletMode = "constant";
    p.BC.left.dirichletValue.u = 0.123;
    p.BC.left.dirichletValue.v = 0.456;

    % Keep the other boundaries periodic
    p.BC.right.type  = "periodic";
    p.BC.top.type    = "periodic";
    p.BC.bottom.type = "periodic";

    model.p = p;
    model.reset();

    % Override default initial condition with a custom smooth field
    model.u = Ubase(:);
    model.v = Vbase(:);

    maxErrDirU = 0.0;
    maxErrDirV = 0.0;

    for k = 1:nStepsShort
        info = model.step();
        assert(~info.hasNaNInf, ...
            "NaN/Inf in Dirichlet BC test in mode %s at step %d.", mode, k);

        [U,V] = model.getFields2D();

        errU = max(abs(U(:,1) - 0.123));
        errV = max(abs(V(:,1) - 0.456));

        maxErrDirU = max(maxErrDirU, errU);
        maxErrDirV = max(maxErrDirV, errV);
    end

    assert(maxErrDirU < tolBC, ...
        "Dirichlet constant failed (u) in mode %s: maxErr=%g", mode, maxErrDirU);
    assert(maxErrDirV < tolBC, ...
        "Dirichlet constant failed (v) in mode %s: maxErr=%g", mode, maxErrDirV);

    % ============================================================
    % Test 2: Neumann derivative consistency (multiple steps)
    % ============================================================
    p = model.p;
    p.BC.left.type = "neumann";
    p.BC.left.neumannMode = "constant";
    p.BC.left.neumannValue.u = 1.0;
    p.BC.left.neumannValue.v = -2.0;

    % Keep the other boundaries periodic
    p.BC.right.type  = "periodic";
    p.BC.top.type    = "periodic";
    p.BC.bottom.type = "periodic";

    model.p = p;
    model.reset();

    % Use the same custom field
    model.u = Ubase(:);
    model.v = Vbase(:);

    hx = model.grid.hx;

    maxErrNeuU = 0.0;
    maxErrNeuV = 0.0;

    for k = 1:nStepsShort
        info = model.step();
        assert(~info.hasNaNInf, ...
            "NaN/Inf in Neumann BC test in mode %s at step %d.", mode, k);

        [U,V] = model.getFields2D();

        gu = (U(:,2) - U(:,1)) / hx;
        gv = (V(:,2) - V(:,1)) / hx;

        errGu = max(abs(gu - 1.0));
        errGv = max(abs(gv + 2.0));

        maxErrNeuU = max(maxErrNeuU, errGu);
        maxErrNeuV = max(maxErrNeuV, errGv);
    end

    assert(maxErrNeuU < tolBC, ...
        "Neumann failed (u) in mode %s: maxErr=%g", mode, maxErrNeuU);
    assert(maxErrNeuV < tolBC, ...
        "Neumann failed (v) in mode %s: maxErr=%g", mode, maxErrNeuV);

    % ============================================================
    % Test 3: Mixed BC robustness (long run)
    % ============================================================
    p = model.p;

    % Dirichlet left (constant)
    p.BC.left.type = "dirichlet";
    p.BC.left.dirichletMode = "constant";
    p.BC.left.dirichletValue.u = 0.2;
    p.BC.left.dirichletValue.v = 0.4;

    % Neumann right (homogeneous)
    p.BC.right.type = "neumann";
    p.BC.right.neumannMode = "constant";
    p.BC.right.neumannValue.u = 0.0;
    p.BC.right.neumannValue.v = 0.0;

    % Periodic top/bottom
    p.BC.top.type    = "periodic";
    p.BC.bottom.type = "periodic";

    model.p = p;
    model.reset();

    % Strong interior-boundary mismatch to make the test less trivial
    U0 = 0.95 * ones(p.Ny, p.Nx);
    V0 = 0.05 * ones(p.Ny, p.Nx);

    model.u = U0(:);
    model.v = V0(:);

    % Enforce BC immediately on the overwritten state
    [uTmp, vTmp] = applyBoundaryConditions(model.u, model.v, p, model.grid, p.BCcache);
    model.u = uTmp;
    model.v = vTmp;

    [Ustart, Vstart] = model.getFields2D();
    probeColUse = min(probeCol, p.Nx-1);
    Uprobe0 = Ustart(:, probeColUse);
    Vprobe0 = Vstart(:, probeColUse);

    maxErrMixDirU = 0.0;
    maxErrMixDirV = 0.0;
    maxErrMixNeuU = 0.0;
    maxErrMixNeuV = 0.0;
    maxProbeDeltaU = 0.0;
    maxProbeDeltaV = 0.0;

    for k = 1:nStepsMixed
        info = model.step();
        assert(~info.hasNaNInf, ...
            "NaN/Inf in mixed BC test in mode %s at step %d.", mode, k);

        [U,V] = model.getFields2D();

        % Dirichlet left must remain exact
        errDLu = max(abs(U(:,1) - 0.2));
        errDLv = max(abs(V(:,1) - 0.4));

        % Neumann right must remain zero-flux
        errNRu = max(abs((U(:,end) - U(:,end-1)) / hx));
        errNRv = max(abs((V(:,end) - V(:,end-1)) / hx));

        % Track maximum boundary errors over the full run
        maxErrMixDirU = max(maxErrMixDirU, errDLu);
        maxErrMixDirV = max(maxErrMixDirV, errDLv);
        maxErrMixNeuU = max(maxErrMixNeuU, errNRu);
        maxErrMixNeuV = max(maxErrMixNeuV, errNRv);

        % Track whether the interior evolves nontrivially
        maxProbeDeltaU = max(maxProbeDeltaU, max(abs(U(:,probeColUse) - Uprobe0)));
        maxProbeDeltaV = max(maxProbeDeltaV, max(abs(V(:,probeColUse) - Vprobe0)));
    end

    assert(maxErrMixDirU < tolBC, ...
        "Mixed BC Dirichlet drift (u) in mode %s: maxErr=%g", mode, maxErrMixDirU);
    assert(maxErrMixDirV < tolBC, ...
        "Mixed BC Dirichlet drift (v) in mode %s: maxErr=%g", mode, maxErrMixDirV);
    assert(maxErrMixNeuU < tolBC, ...
        "Mixed BC Neumann failed (u) in mode %s: maxErr=%g", mode, maxErrMixNeuU);
    assert(maxErrMixNeuV < tolBC, ...
        "Mixed BC Neumann failed (v) in mode %s: maxErr=%g", mode, maxErrMixNeuV);

    % The run should not be completely trivial in the interior
    assert(maxProbeDeltaU > 1e-8 || maxProbeDeltaV > 1e-8, ...
        "Mixed BC test in mode %s showed no measurable interior evolution at probe column %d.", ...
        mode, probeColUse);

    % Store metrics for report writing
    metrics.(mode).dirichlet.maxErrU = maxErrDirU;
    metrics.(mode).dirichlet.maxErrV = maxErrDirV;
    metrics.(mode).neumann.maxErrU   = maxErrNeuU;
    metrics.(mode).neumann.maxErrV   = maxErrNeuV;
    metrics.(mode).mixed.maxDirU     = maxErrMixDirU;
    metrics.(mode).mixed.maxDirV     = maxErrMixDirV;
    metrics.(mode).mixed.maxNeuU     = maxErrMixNeuU;
    metrics.(mode).mixed.maxNeuV     = maxErrMixNeuV;
    metrics.(mode).mixed.probeDeltaU = maxProbeDeltaU;
    metrics.(mode).mixed.probeDeltaV = maxProbeDeltaV;

    fprintf("[%s] Dirichlet (%d steps): maxErrU=%.3e, maxErrV=%.3e\n", ...
        mode, nStepsShort, maxErrDirU, maxErrDirV);
    fprintf("[%s] Neumann   (%d steps): maxErrU=%.3e, maxErrV=%.3e\n", ...
        mode, nStepsShort, maxErrNeuU, maxErrNeuV);
    fprintf("[%s] Mixed     (%d steps): maxDirU=%.3e, maxDirV=%.3e, maxNeuU=%.3e, maxNeuV=%.3e, probeDeltaU=%.3e, probeDeltaV=%.3e\n", ...
        mode, nStepsMixed, maxErrMixDirU, maxErrMixDirV, maxErrMixNeuU, maxErrMixNeuV, maxProbeDeltaU, maxProbeDeltaV);
end

r = makeResult("test_boundary_conditions", p);
r.grid = [p.Nx, p.Ny];
r.dt   = p.dt;
r.T    = p.T;
r.pass = true;
r.metrics = metrics;
r.thresholds.boundaryTol = tolBC;
r.thresholds.nStepsShort = nStepsShort;
r.thresholds.nStepsMixed = nStepsMixed;
r.notes = "Quantitative Dirichlet, Neumann, and mixed BC checks passed for matrix, stencil, and full modes.";

end