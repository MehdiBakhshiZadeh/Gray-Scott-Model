function r = test_laplacian_symmetry()
%TEST_LAPLACIAN_SYMMETRY  Check that the discrete Laplacian is symmetric.

% Locate this test folder and project root
thisFile = mfilename("fullpath");
thisDir  = fileparts(thisFile);   % .../tests/supplementary_checks
testsDir = fileparts(thisDir);    % .../tests
rootDir  = fileparts(testsDir);   % project root

% Make this test runnable independently
if isempty(which("defaultParams")) || isempty(which("buildLaplacian2D"))
    addpath(rootDir);
    addpath(genpath(fullfile(rootDir, "src")));
    addpath(genpath(testsDir));
end

% Local output folder for this test group (reserved for any future outputs)
outDir = fullfile(thisDir, "outputs");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

p = defaultParams();
grid = buildGrid(p);
L = buildLaplacian2D(p, grid);

A = L - L.';               % difference from transpose
infNorm = norm(A, inf);

% Standardized result struct
r = makeResult("test_laplacian_symmetry", p);
r.grid = [p.Nx, p.Ny];
if isfield(p, "dt"), r.dt = p.dt; end
r.T = p.T;

r.metrics.errLinf = infNorm;

% For periodic, centered finite differences, the Laplacian matrix should be
% symmetric; any asymmetry should be at floating-point roundoff level.
tol = 1e-12;
r.thresholds.errLinf = tol;

r.pass = (infNorm < tol);
r.notes = sprintf("||L - L^T||_inf = %.3e", infNorm);
end