function Lu = applyLaplacianStencil(u, S)
%APPLYLAPLACIANSTENCIL  Matrix-free 2D Laplacian apply using a 5-point stencil.
%
%   Lu = applyLaplacianStencil(u, S) computes Lu = L*u without forming L,
%   where L matches buildLaplacian2D for periodic BC and u = U(:) with
%   U sized Ny-by-Nx (column-major).
%
% Inputs:
%   u : vectorized field (Nx*Ny x 1) or (1 x Nx*Ny)
%   S : stencil struct from buildStencil2D (fields: Nx, Ny, bc, inv_hx2, inv_hy2)
%
% Output:
%   Lu : vectorized Laplacian (same size as u)

% --- Validate inputs ---
assert(isstruct(S), "applyLaplacianStencil: S must be a struct.");
required = ["Nx","Ny","bc","inv_hx2","inv_hy2"];
for k = 1:numel(required)
    assert(isfield(S, required(k)), "applyLaplacianStencil: missing field S.%s", required(k));
end

Nx = S.Nx; Ny = S.Ny;

u = u(:);
assert(numel(u) == Nx*Ny, "applyLaplacianStencil: u must have length Nx*Ny.");

bc = string(S.bc);
if bc ~= "periodic"
    error("applyLaplacianStencil: only periodic BC is implemented.");
end

% --- Reshape to 2D field (Ny-by-Nx) ---
U = reshape(u, Ny, Nx);

% --- Periodic neighbors ---
U_e = circshift(U, [ 0, -1]);  % east  (j+1)
U_w = circshift(U, [ 0,  1]);  % west  (j-1)
U_n = circshift(U, [-1,  0]);  % north (i+1) in y-direction
U_s = circshift(U, [ 1,  0]);  % south (i-1)

% --- 5-point Laplacian: (E - 2C + W)/hx^2 + (N - 2C + S)/hy^2 ---
LU = S.inv_hx2 * (U_e - 2*U + U_w) + S.inv_hy2 * (U_n - 2*U + U_s);

% --- Return vectorized ---
Lu = LU(:);
end
