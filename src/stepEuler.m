function [unew, vnew, info] = stepEuler(u, v, L, S, p, t)
%STEPEULER  One explicit Euler step for the Gray–Scott model.
%
% Inputs:
%   u, v : current state vectors (N x 1)
%   L    : sparse Laplacian operator (N x N)
%   p    : parameter struct containing Du, Dv, F, k, dt
%   t    : current time (optional, used only for manufactured solutions)
%
% Outputs:
%   unew, vnew : updated state vectors (N x 1)
%   info       : diagnostics (min/max values and NaN/Inf flag)

if nargin < 5
    t = 0;
end

assert(isequal(size(u), size(v)), "stepEuler: u and v must have the same size.");

% ---- Diffusion contributions ----
mode = "matrix";
if isfield(p,"diffusionMode") && ~isempty(p.diffusionMode)
    mode = string(p.diffusionMode);
end

if mode == "stencil"
    du_diff = p.Du * applyLaplacianStencil(u, S);
    dv_diff = p.Dv * applyLaplacianStencil(v, S);

elseif mode == "full"
    du_diff = p.Du * applyFullLaplacian(u, L);
    dv_diff = p.Dv * applyFullLaplacian(v, L);

else
    % "matrix" (sparse)
    du_diff = p.Du * (L * u);
    dv_diff = p.Dv * (L * v);
end



% ---- Reaction contributions (pointwise kinetics) ----
[du_reac, dv_reac] = reactionGrayScott(u, v, p);

% ---- Optional source terms (MMS verification only) ----
% These terms are used for the Method of Manufactured Solutions.
% They are inactive unless p.sourceFcn is provided.
su = zeros(size(u));
sv = su;

if isfield(p, 'sourceFcn') && ~isempty(p.sourceFcn)
    % Expect grid coordinates in p.grid.x and p.grid.y
    x = p.grid.x;
    y = p.grid.y;

    [su2D, sv2D] = p.sourceFcn(x, y, t, p);
    su = su2D(:);
    sv = sv2D(:);
end

% ---- Explicit Euler update ----
unew = u + p.dt * (du_diff + du_reac + su);
vnew = v + p.dt * (dv_diff + dv_reac + sv);

% ---- Diagnostics for stability and safety ----
info.u_min = min(unew);
info.u_max = max(unew);
info.v_min = min(vnew);
info.v_max = max(vnew);
info.hasNaNInf = any(~isfinite(unew)) || any(~isfinite(vnew));

if info.hasNaNInf
    warning("NaN or Inf detected in solution after Euler step.");
end

end
