function S = buildStencil2D(p, grid)
%BUILDSTENCIL2D Build a matrix-free stencil representation of the 2D Laplacian.
%
% S stores only what is needed to apply the standard second-order 5-point
% Laplacian stencil on a uniform grid with the same data layout convention
% as buildLaplacian2D:
%   - U is Ny-by-Nx
%   - u is vectorized as u = U(:) in MATLAB column-major order
%
% Inputs:
%   p    : parameter struct (requires p.Nx, p.Ny)
%   grid : grid struct (requires grid.hx, grid.hy)
%
% Output:
%   S : struct with fields:
%       Nx, Ny     : grid sizes
%       inv_hx2    : 1/hx^2
%       inv_hy2    : 1/hy^2
%
% Note:
%   Boundary conditions (Dirichlet/Neumann) are enforced outside the stencil
%   after each time step. The stencil itself assumes periodic connectivity.

% --- Input validation ---
requiredP = ["Nx","Ny"];
for k = 1:numel(requiredP)
    assert(isfield(p, requiredP(k)), "buildStencil2D: missing parameter p.%s", requiredP(k));
end
assert(isfield(grid,"hx") && isfield(grid,"hy"), ...
    "buildStencil2D: grid must contain hx and hy.");

Nx = p.Nx;
Ny = p.Ny;

assert(Nx >= 3 && Ny >= 3, "buildStencil2D: Nx, Ny must be >= 3.");
assert(grid.hx > 0 && grid.hy > 0, "buildStencil2D: hx and hy must be positive.");

% --- Stencil data ---
S = struct();
S.Nx = Nx;
S.Ny = Ny;

S.inv_hx2 = 1/(grid.hx^2);
S.inv_hy2 = 1/(grid.hy^2);

end