function p = finalizeParams(p)
%FINALIZEPARAMS  Apply derived/default mappings after user edits.
%
% Ensures p.diffusionMode is consistent with p.solver.

if isfield(p, "solver") && ~isempty(p.solver)
    switch lower(string(p.solver))
        case "sparse"
            p.diffusionMode = "matrix";
        case "dense"
            p.diffusionMode = "full";
        case "stencil"
            p.diffusionMode = "stencil";
        otherwise
            error("Unknown p.solver='%s' (use sparse/dense/stencil)", string(p.solver));
    end
end
end
