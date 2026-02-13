function r = test_laplacian_constant()
%TEST_LAPLACIAN_CONSTANT  Laplacian of a constant field should be ~0.

p = defaultParams();
grid = buildGrid(p);
L = buildLaplacian2D(p, grid);

N = p.Nx * p.Ny;
u = ones(N, 1);
Lu = L * u;

infNorm = norm(Lu, inf);
twoNorm = norm(Lu, 2);

% Standardized result struct
r = makeResult("test_laplacian_constant", p);
r.grid = [p.Nx, p.Ny];
if isfield(p, "dt"), r.dt = p.dt; end
r.T = p.T;

r.metrics.errLinf  = infNorm;
r.metrics.errL2vec = twoNorm;   % vector 2-norm (not spatially weighted)

% For a periodic finite-difference Laplacian, the constant vector lies in the
% nullspace; any nonzero result should be at floating-point roundoff level.
tol = 1e-12;
r.thresholds.errLinf  = tol;
r.thresholds.errL2vec = tol;

r.pass = (infNorm < tol) && (twoNorm < tol);
r.notes = sprintf("||L*1||_inf=%.3e, ||L*1||_2=%.3e", infNorm, twoNorm);
end
