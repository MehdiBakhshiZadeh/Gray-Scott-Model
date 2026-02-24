function op = makeOperator(p, grid)
%MAKEOPERATOR  Dispatcher for diffusion operators (uses p.diffusionMode).

assert(isfield(p,"diffusionMode") && ~isempty(p.diffusionMode), ...
    "makeOperator: missing parameter p.diffusionMode.");
mode = lower(strtrim(string(p.diffusionMode)));

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
