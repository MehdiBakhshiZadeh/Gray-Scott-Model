function L = buildLaplacian2D(p, grid)
%BUILDLAPLACIAN2D  Build the 2D Laplacian operator as a sparse matrix.
%
%   L = BUILDLAPLACIAN2D(p, grid) returns the sparse matrix that applies the
%   standard second-order (5-point) finite-difference Laplacian on a uniform,
%   periodic 2D grid.
%
%   Inputs:
%     p.Nx, p.Ny : number of grid points in x and y (integers >= 3)
%     p.bc       : boundary condition type (currently only "periodic")
%     grid.hx, grid.hy : grid spacings in x and y (positive)
%
%   Output:
%     L : sparse matrix of size (Nx*Ny) x (Nx*Ny) such that, if U is Ny-by-Nx
%         and u = U(:) (MATLAB column-major order), then L*u approximates
%         ?^2 U flattened into a vector.
%
%   Data layout convention:
%     U(j,i) corresponds to vector index k = j + (i-1)*Ny, with
%       j = 1..Ny (row index), i = 1..Nx (column index).
%     This convention determines the Kronecker-product ordering below.

% --- Input validation ---
requiredP = ["Nx","Ny","bc"];
for k = 1:numel(requiredP)
    assert(isfield(p, requiredP(k)), "buildLaplacian2D: missing parameter p.%s", requiredP(k));
end
assert(isfield(grid, "hx") && isfield(grid, "hy"), ...
    "buildLaplacian2D: grid must contain hx and hy.");

Nx = p.Nx;
Ny = p.Ny;

assert(mod(Nx,1) == 0 && mod(Ny,1) == 0, "buildLaplacian2D: Nx and Ny must be integers.");
assert(Nx >= 3 && Ny >= 3, "buildLaplacian2D: Nx and Ny must be >= 3.");
assert(grid.hx > 0 && grid.hy > 0, "buildLaplacian2D: hx and hy must be positive.");

bc = string(p.bc);
if bc ~= "periodic"
    error("buildLaplacian2D: only periodic boundary condition is implemented (p.bc = ""periodic"").");
end

% --- 1D Laplacians (unscaled stencils) ---
Lx1 = buildLaplacian1D_periodic(Nx); % acts along x (columns)
Ly1 = buildLaplacian1D_periodic(Ny); % acts along y (rows)

Ix = speye(Nx);
Iy = speye(Ny);

% --- Scale factors (include 1/h^2 here) ---
hx2 = grid.hx^2;
hy2 = grid.hy^2;

% --- 2D Laplacian (Kronecker sum consistent with u = U(:)) ---
L = (1/hy2) * kron(Ix, Ly1) + (1/hx2) * kron(Lx1, Iy);

end
