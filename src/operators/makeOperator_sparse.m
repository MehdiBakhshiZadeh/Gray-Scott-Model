function op = makeOperator_sparse(p, grid)
%MAKEOPERATOR_SPARSE  Sparse-matrix Laplacian operator.

op.mode = "matrix";
op.grid = grid;
op.L = buildLaplacian2D(p, grid);  % sparse
op.S = [];
end
