function r = test_laplacian_stencil_equivalence()
%TEST_LAPLACIAN_STENCIL_EQUIVALENCE
% Compare sparse Laplacian vs stencil apply for periodic connectivity.
%
% Note:
% The diffusion operators are assembled with periodic connectivity. Physical
% boundary conditions (Dirichlet/Neumann) are enforced after each time step
% and are therefore not part of this operator equivalence test.
%
% This version is report-friendly:
% - runnable independently
% - prints quantitative worst-case errors
% - stores metrics in the returned result struct

% Locate project folders relative to this file
thisDir = fileparts(mfilename("fullpath"));   % .../tests
rootDir = fileparts(thisDir);                 % project root

% Make this test runnable independently
if isempty(which("defaultParams"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(fullfile(rootDir, "tests")));
end

testName = "test_laplacian_stencil_equivalence";
r = makeResult(testName, struct());   % standardized struct
r.pass  = false;                      % default fail until proven otherwise
r.notes = "";

try
    % Reproducible random trials
    rng(1, "twister");

    % Multiple sizes catch reshape/index mistakes
    sizes = [ ...
        32  32; ...
        64  48; ...
        128 128 ...
    ];

    nTrials = 3;
    relTol  = 1e-12;
    absTol  = 1e-10;

    % Per-case metrics
    caseMetrics = struct( ...
        "Nx", {}, "Ny", {}, "trial", {}, ...
        "relerr", {}, "maxabs", {} );

    % Worst-case trackers
    worstRelErr = -inf;
    worstAbsErr = -inf;
    worstRelNx  = NaN;
    worstRelNy  = NaN;
    worstRelTrial = NaN;
    worstAbsNx  = NaN;
    worstAbsNy  = NaN;
    worstAbsTrial = NaN;

    idx = 0;

    for s = 1:size(sizes,1)
        Nx = sizes(s,1);
        Ny = sizes(s,2);

        p = defaultParams();
        p.Nx = Nx;
        p.Ny = Ny;

        grid = buildGrid(p);
        L = buildLaplacian2D(p, grid);
        S = buildStencil2D(p, grid);

        N = Nx * Ny;

        for k = 1:nTrials
            u = randn(N,1);

            a = L * u;
            b = applyLaplacianStencil(u, S);

            diffVec = a - b;

            relerr = norm(diffVec) / max(1, norm(a));
            maxabs = max(abs(diffVec));

            idx = idx + 1;
            caseMetrics(idx).Nx = Nx;
            caseMetrics(idx).Ny = Ny;
            caseMetrics(idx).trial = k;
            caseMetrics(idx).relerr = relerr;
            caseMetrics(idx).maxabs = maxabs;

            if relerr > worstRelErr
                worstRelErr   = relerr;
                worstRelNx    = Nx;
                worstRelNy    = Ny;
                worstRelTrial = k;
            end

            if maxabs > worstAbsErr
                worstAbsErr   = maxabs;
                worstAbsNx    = Nx;
                worstAbsNy    = Ny;
                worstAbsTrial = k;
            end

            fprintf('[%dx%d, trial %d] relerr=%.3e, maxabs=%.3e\n', ...
                Nx, Ny, k, relerr, maxabs);

            assert(relerr < relTol, ...
                "relerr too large (Nx=%d,Ny=%d,trial=%d): %.3e (tol=%.1e)", ...
                Nx, Ny, k, relerr, relTol);

            assert(maxabs < absTol, ...
                "maxabs too large (Nx=%d,Ny=%d,trial=%d): %.3e (tol=%.1e)", ...
                Nx, Ny, k, maxabs, absTol);
        end
    end

    fprintf('\nWorst relative error: %.3e at [%dx%d], trial %d\n', ...
        worstRelErr, worstRelNx, worstRelNy, worstRelTrial);
    fprintf('Worst max absolute difference: %.3e at [%dx%d], trial %d\n', ...
        worstAbsErr, worstAbsNx, worstAbsNy, worstAbsTrial);
    fprintf('Tolerances: relTol=%.1e, absTol=%.1e\n', relTol, absTol);

    r.pass = true;
    r.grid = sizes;
    r.metrics.caseMetrics = caseMetrics;
    r.metrics.worstRelErr = worstRelErr;
    r.metrics.worstAbsErr = worstAbsErr;
    r.metrics.worstRelCase = [worstRelNx, worstRelNy, worstRelTrial];
    r.metrics.worstAbsCase = [worstAbsNx, worstAbsNy, worstAbsTrial];
    r.thresholds.relTol = relTol;
    r.thresholds.absTol = absTol;
    r.thresholds.nTrials = nTrials;
    r.notes = sprintf([ ...
        'PASS; worstRelErr=%.3e at [%dx%d], trial %d; ' ...
        'worstAbsErr=%.3e at [%dx%d], trial %d'], ...
        worstRelErr, worstRelNx, worstRelNy, worstRelTrial, ...
        worstAbsErr, worstAbsNx, worstAbsNy, worstAbsTrial);

catch ME
    r.pass  = false;
    r.notes = "Crashed: " + string(ME.message);
end
end