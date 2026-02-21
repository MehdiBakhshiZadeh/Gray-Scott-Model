function p = finalizeParams(p)
%FINALIZEPARAMS Apply derived/default mappings after user edits.
%
% Ensures p.diffusionMode is consistent with p.solver and normalizes the
% boundary-condition specification p.BC (per-side periodic/dirichlet/neumann).

% ---- Solver ↔ diffusionMode mapping ----
if isfield(p, "solver") && ~isempty(p.solver)
    switch lower(string(p.solver))
        case "sparse"
            p.diffusionMode = "matrix";
        case "dense"
            p.diffusionMode = "full";
        case "stencil"
            p.diffusionMode = "stencil";
        otherwise
            error("finalizeParams: unknown p.solver='%s' (use sparse/dense/stencil)", string(p.solver));
    end
end

% ---- Boundary-condition normalization ----
p = normalizeBC(p);

end

% ======================= helper =======================

function p = normalizeBC(p)
%NORMALIZEBC Ensure p.BC exists and has all required subfields.

sides = ["left","right","bottom","top"];

% If p.BC is missing, create periodic defaults (non-breaking).
if ~isfield(p, "BC") || isempty(p.BC)
    p.BC = struct();
end

for i = 1:numel(sides)
    side = sides(i);

    if ~isfield(p.BC, side) || isempty(p.BC.(side))
        p.BC.(side) = struct();
    end

    bc = p.BC.(side);

    % type
    if ~isfield(bc, "type") || isempty(bc.type)
        bc.type = "periodic";
    end
    bc.type = lower(string(bc.type));

    if ~any(bc.type == ["periodic","dirichlet","neumann"])
        error("finalizeParams: invalid p.BC.%s.type='%s'.", side, bc.type);
    end

    % dirichletMode + values
    if ~isfield(bc, "dirichletMode") || isempty(bc.dirichletMode)
        bc.dirichletMode = "initial";
    end
    bc.dirichletMode = lower(string(bc.dirichletMode));

    if ~any(bc.dirichletMode == ["initial","constant"])
        error("finalizeParams: invalid p.BC.%s.dirichletMode='%s'.", side, bc.dirichletMode);
    end

    if ~isfield(bc, "dirichletValue") || isempty(bc.dirichletValue)
        bc.dirichletValue = struct();
    end
    if ~isfield(bc.dirichletValue, "u"), bc.dirichletValue.u = 0.0; end
    if ~isfield(bc.dirichletValue, "v"), bc.dirichletValue.v = 0.0; end

    % neumannMode + values
    if ~isfield(bc, "neumannMode") || isempty(bc.neumannMode)
        bc.neumannMode = "constant";
    end
    bc.neumannMode = lower(string(bc.neumannMode));

    if bc.neumannMode ~= "constant"
        error("finalizeParams: invalid p.BC.%s.neumannMode='%s'.", side, bc.neumannMode);
    end

    if ~isfield(bc, "neumannValue") || isempty(bc.neumannValue)
        bc.neumannValue = struct();
    end
    if ~isfield(bc.neumannValue, "u"), bc.neumannValue.u = 0.0; end
    if ~isfield(bc.neumannValue, "v"), bc.neumannValue.v = 0.0; end

    % Store back
    p.BC.(side) = bc;
end

end