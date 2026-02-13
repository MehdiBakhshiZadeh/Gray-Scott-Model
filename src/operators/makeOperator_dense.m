function op = makeOperator_dense(p, grid)
%MAKEOPERATOR_DENSE  Dense (full) Laplacian operator (small grids only).

Nx = p.Nx; Ny = p.Ny;
allowed = [32 64 128];
assert(any(Nx == allowed) && any(Ny == allowed), ...
    "Dense Laplacian only allowed for Nx,Ny in {32,64,128}. Got Nx=%d, Ny=%d.", Nx, Ny);

op.mode = "full";
op.grid = grid;

Ls = buildLaplacian2D(p, grid); % sparse builder
op.L = full(Ls);                % dense
op.S = [];
end
