function r = test_laplacian_symmetry()
%TEST_LAPLACIAN_SYMMETRY  Check that the discrete Laplacian is symmetric.

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
