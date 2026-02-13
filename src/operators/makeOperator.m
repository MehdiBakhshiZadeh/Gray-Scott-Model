function op = makeOperator(p, grid)
%MAKEOPERATOR  Dispatcher for diffusion operators (preferred: p.solver).
%
% Uses p.solver if present; falls back to p.diffusionMode for compatibility.

if isfield(p,"solver") && ~isempty(p.solver)
    s = lower(string(p.solver));
    switch s
        case "sparse"
            op = makeOperator_sparse(p, grid);
        case "dense"
            op = makeOperator_dense(p, grid);
        case "stencil"
            op = makeOperator_stencil(p, grid);
        otherwise
            error("Unknown p.solver='%s' (use sparse/dense/stencil)", s);
    end
    return;
end

% ---- Fallback (legacy) ----
mode = "matrix";
if isfield(p,"diffusionMode") && ~isempty(p.diffusionMode)
    mode = string(p.diffusionMode);
end

switch mode
    case "matrix"
        op = makeOperator_sparse(p, grid);
    case "full"
        op = makeOperator_dense(p, grid);
    case "stencil"
        op = makeOperator_stencil(p, grid);
    otherwise
        error("Unknown diffusionMode='%s'", mode);
end
end
