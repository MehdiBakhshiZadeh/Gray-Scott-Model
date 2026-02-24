function p = finalizeParams(p)
%FINALIZEPARAMS Normalize parameters after user edits.
%
% Ensures p.diffusionMode is present and normalized, and normalizes the
% boundary-condition specification p.BC (per-side periodic/dirichlet/neumann).
%
% Canonical diffusion selector:
%   p.diffusionMode = "matrix" | "stencil" | "full"

% ---- Diffusion mode normalization ----
if ~isfield(p, "diffusionMode") || isempty(p.diffusionMode)
    error("finalizeParams: missing parameter p.diffusionMode.");
end
p.diffusionMode = lower(strtrim(string(p.diffusionMode)));

if ~any(p.diffusionMode == ["matrix","stencil","full"])
    error("finalizeParams: invalid p.diffusionMode='%s' (use matrix/stencil/full).", p.diffusionMode);
end

% ---- Boundary-condition normalization ----
p = normalizeBC(p);

end

% ======================= helper =======================

function p = normalizeBC(p)
%NORMALIZEBC Ensure p.BC exists and contains valid per-side settings.
% This function fills missing fields with safe defaults and validates
% supported options: periodic, dirichlet (initial/constant), neumann (constant).

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
    bc.type = lower(strtrim(string(bc.type)));

    if ~any(bc.type == ["periodic","dirichlet","neumann"])
        error("finalizeParams: invalid p.BC.%s.type='%s'.", side, bc.type);
    end

    % dirichletMode + values
    if ~isfield(bc, "dirichletMode") || isempty(bc.dirichletMode)
        bc.dirichletMode = "initial";
    end
    bc.dirichletMode = lower(strtrim(string(bc.dirichletMode)));

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
    bc.neumannMode = lower(strtrim(string(bc.neumannMode)));

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