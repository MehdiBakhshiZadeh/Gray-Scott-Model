function grid = buildGrid(p)
%BUILDGRID  Construct a uniform 2D grid for the simulation domain.
%
%   grid = BUILDGRID(p) creates a uniformly-spaced grid on
%   [0, Lx) x [0, Ly) using Nx-by-Ny points.
%
%   Inputs:
%     p.Lx, p.Ly : domain lengths in x and y (positive)
%     p.Nx, p.Ny : number of grid points in x and y (positive integers)
%
%   Outputs (in struct grid):
%     grid.hx, grid.hy : grid spacings (periodic spacing: hx = Lx/Nx)
%     grid.x,  grid.y  : grid coordinates (row vectors)
%
%   Note:
%     For periodic grids, the endpoint is not duplicated, so hx = Lx/Nx
%     (not Lx/(Nx-1)).

% --- Basic parameter checks ---
required = ["Lx","Ly","Nx","Ny"];
for k = 1:numel(required)
    assert(isfield(p, required(k)), "buildGrid: missing parameter p.%s", required(k));
end

assert(p.Lx > 0 && p.Ly > 0, "buildGrid: Lx and Ly must be positive.");
assert(p.Nx > 0 && p.Ny > 0, "buildGrid: Nx and Ny must be positive.");
assert(mod(p.Nx,1) == 0 && mod(p.Ny,1) == 0, "buildGrid: Nx and Ny must be integers.");

% --- Uniform periodic grid on [0,L) ---
grid.hx = p.Lx / p.Nx;
grid.hy = p.Ly / p.Ny;

grid.x = (0:p.Nx-1) * grid.hx;  % row vector
grid.y = (0:p.Ny-1) * grid.hy;  % row vector

grid.Nx = p.Nx;
grid.Ny = p.Ny;
grid.Lx = p.Lx;
grid.Ly = p.Ly;

end
