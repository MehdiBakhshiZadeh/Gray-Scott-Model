function S = buildStencil2D(p, grid)
%BUILDSTENCIL2D  Build a matrix-free stencil representation of the 2D Laplacian.
%
% S stores only what is needed to apply the 5-point Laplacian stencil
% on a uniform grid with the same convention as buildLaplacian2D:
%   u is vectorized as U(:) with U size Ny-by-Nx.
%
% Inputs:
%   p    : parameter struct (requires p.Nx, p.Ny, p.bc)
%   grid : grid struct (requires grid.hx, grid.hy)
%
% Output:
%   S : struct with fields Nx, Ny, bc, inv_hx2, inv_hy2

% --- Input validation ---
requiredP = ["Nx","Ny","bc"];
for k = 1:numel(requiredP)
    assert(isfield(p, requiredP(k)), "buildStencil2D: missing parameter p.%s", requiredP(k));
end
assert(isfield(grid,"hx") && isfield(grid,"hy"), "buildStencil2D: grid must contain hx and hy.");

Nx = p.Nx; Ny = p.Ny;
assert(Nx >= 3 && Ny >= 3, "buildStencil2D: Nx, Ny must be >= 3.");

bc = string(p.bc);
if bc ~= "periodic"
    error("buildStencil2D: only periodic boundary condition is implemented (p.bc = ""periodic"").");
end

S.Nx = Nx;
S.Ny = Ny;
S.bc = bc;

S.inv_hx2 = 1/(grid.hx^2);
S.inv_hy2 = 1/(grid.hy^2);
end
