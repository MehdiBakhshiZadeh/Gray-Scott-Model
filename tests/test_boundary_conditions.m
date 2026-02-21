function test_boundary_conditions()
%TEST_BOUNDARY_CONDITIONS Quantitative verification of Dirichlet/Neumann BC.
%
% This test verifies that boundary conditions are enforced strongly after
% each explicit Euler update (Option A), and that enforcement is consistent
% across all solver modes:
%   - "matrix"  : sparse Laplacian
%   - "stencil" : matrix-free 5-point stencil
%   - "full"    : dense Laplacian (small grid only)
%
% The test performs:
%   1) Dirichlet constant enforcement check (1 step).
%   2) Neumann derivative consistency check (1 step).
%   3) Mixed BC robustness check (100 steps).
%
% The same BC specification is applied to both u and v.

modes = ["matrix","stencil","full"];
for m = 1:numel(modes)
    mode = modes(m);

    % --- Parameters and grid size (Option A from report) ---
    p = defaultParams();
    p.solver = modeToSolver(p, mode);

    p.Nx = 64;
    p.Ny = 64;

    % Keep runtime small but nontrivial
    p.dt = 0.25;
    p.T  = 100 * p.dt;

    p = finalizeParams(p);

    % Build model (creates BCcache from initial condition)
    grid  = buildGrid(p);
    model = GrayScottModel(p, grid);

    % Use model.p to ensure BCcache is present
    p = model.p;

    % ---------------- Test 1: Dirichlet constant (1 step) ----------------
    p.BC.left.type = "dirichlet";
    p.BC.left.dirichletMode = "constant";
    p.BC.left.dirichletValue.u = 0.123;
    p.BC.left.dirichletValue.v = 0.456;

    model.p = p;
    model.reset();
    model.step();

    [U,V] = model.getFields2D();
    errU = max(abs(U(:,1) - 0.123));
    errV = max(abs(V(:,1) - 0.456));

    assert(errU < 1e-12, "Dirichlet constant failed (u) in mode %s: err=%g", mode, errU);
    assert(errV < 1e-12, "Dirichlet constant failed (v) in mode %s: err=%g", mode, errV);

    % ---------------- Test 2: Neumann derivative (1 step) ----------------
    p = model.p; % reset p reference (keeps BCcache)
    p.BC.left.type = "neumann";
    p.BC.left.neumannMode = "constant";
    p.BC.left.neumannValue.u = 1.0;
    p.BC.left.neumannValue.v = -2.0;

    model.p = p;
    model.reset();
    model.step();

    [U,V] = model.getFields2D();
    hx = model.grid.hx;

    gu = (U(:,2) - U(:,1)) / hx;
    gv = (V(:,2) - V(:,1)) / hx;

    errGu = max(abs(gu - 1.0));
    errGv = max(abs(gv + 2.0));

    assert(errGu < 1e-12, "Neumann failed (u) in mode %s: err=%g", mode, errGu);
    assert(errGv < 1e-12, "Neumann failed (v) in mode %s: err=%g", mode, errGv);

    % ------------- Test 3: Mixed BC robustness (100 steps) --------------
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

    % Periodic top/bottom (default)
    p.BC.top.type    = "periodic";
    p.BC.bottom.type = "periodic";

    model.p = p;
    model.reset();

    for k = 1:100
        info = model.step();
        assert(~info.hasNaNInf, "NaN/Inf in mixed BC test in mode %s at step %d.", mode, k);

        [U,V] = model.getFields2D();

        % Dirichlet must remain exact
        errDLu = max(abs(U(:,1) - 0.2));
        errDLv = max(abs(V(:,1) - 0.4));
        assert(errDLu < 1e-12, "Mixed BC Dirichlet drift (u) in mode %s: err=%g", mode, errDLu);
        assert(errDLv < 1e-12, "Mixed BC Dirichlet drift (v) in mode %s: err=%g", mode, errDLv);

        % Neumann right: (U(:,end)-U(:,end-1))/hx = 0
        errNRu = max(abs((U(:,end) - U(:,end-1)) / hx));
        errNRv = max(abs((V(:,end) - V(:,end-1)) / hx));
        assert(errNRu < 1e-12, "Mixed BC Neumann failed (u) in mode %s: err=%g", mode, errNRu);
        assert(errNRv < 1e-12, "Mixed BC Neumann failed (v) in mode %s: err=%g", mode, errNRv);
    end
end

disp("test_boundary_conditions: PASSED");

end

% ===================== local helper =====================

function solver = modeToSolver(p, mode)
%MODETOSOLVER map diffusion mode name to p.solver legacy interface.
%
% Your code uses p.solver as a user-facing knob and finalizeParams maps it
% to p.diffusionMode. Keep this mapping centralized in the test.

mode = lower(string(mode));
switch mode
    case "matrix"
        solver = "sparse";
    case "stencil"
        solver = "stencil";
    case "full"
        solver = "dense";
    otherwise
        error("Unknown mode '%s'.", mode);
end
end

