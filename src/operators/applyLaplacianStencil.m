function Lu = applyLaplacianStencil(u, S)
%APPLYLAPLACIANSTENCIL Matrix-free 2D Laplacian apply using a 5-point stencil.
%
%   Lu = applyLaplacianStencil(u, S) computes Lu = L*u without forming L,
%   where L matches the periodic-connectivity Laplacian assembled by
%   buildLaplacian2D and u = U(:) with U sized Ny-by-Nx (column-major).
%
% Inputs:
%   u : vectorized field (Nx*Ny x 1) or (1 x Nx*Ny)
%   S : stencil struct from buildStencil2D (fields: Nx, Ny, inv_hx2, inv_hy2)
%
% Output:
%   Lu : vectorized Laplacian (same size as u)
%
% Note:
%   Boundary conditions (Dirichlet/Neumann) are enforced outside the diffusion
%   operator after each time step. The stencil itself therefore uses periodic
%   connectivity unconditionally via circular shifts.

% --- Validate inputs ---
assert(isstruct(S), "applyLaplacianStencil: S must be a struct.");
required = ["Nx","Ny","inv_hx2","inv_hy2"];
for k = 1:numel(required)
    assert(isfield(S, required(k)), "applyLaplacianStencil: missing field S.%s", required(k));
end

Nx = S.Nx;
Ny = S.Ny;

u = u(:);
assert(numel(u) == Nx*Ny, "applyLaplacianStencil: u must have length Nx*Ny.");

% --- Reshape to 2D field (Ny-by-Nx) ---
U = reshape(u, Ny, Nx);

% --- Periodic neighbors (periodic connectivity) ---
U_e = circshift(U, [ 0, -1]);  % east  (x+)
U_w = circshift(U, [ 0,  1]);  % west  (x-)
U_n = circshift(U, [-1,  0]);  % north (y+)
U_s = circshift(U, [ 1,  0]);  % south (y-)

% --- 5-point Laplacian: (E - 2C + W)/hx^2 + (N - 2C + S)/hy^2 ---
LU = S.inv_hx2 * (U_e - 2*U + U_w) + S.inv_hy2 * (U_n - 2*U + U_s);

% --- Return vectorized ---
Lu = LU(:);
end