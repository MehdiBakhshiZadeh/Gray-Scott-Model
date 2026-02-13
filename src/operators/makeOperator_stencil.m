function op = makeOperator_stencil(p, grid)
%MAKEOPERATOR_STENCIL  Stencil-based operator.

op.mode = "stencil";
op.grid = grid;
op.L = [];
op.S = buildStencil2D(p, grid);
end
