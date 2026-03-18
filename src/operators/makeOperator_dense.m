function op = makeOperator_dense(p, grid)
%MAKEOPERATOR_DENSE  Dense (full) Laplacian operator (small grids only).

Nx = p.Nx; Ny = p.Ny;
N = Nx * Ny;
assert(N <= 128*128, ...
    "Dense Laplacian only allowed for small grids. Got Nx=%d, Ny=%d.", Nx, Ny);

op.mode = "full";
op.grid = grid;

Ls = buildLaplacian2D(p, grid); % sparse builder
op.L = full(Ls);                % dense
op.S = [];
end
