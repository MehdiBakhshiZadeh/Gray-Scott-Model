function [unew, vnew, info] = eulerStep(u, v, op, p, t)
%EULERSTEP One explicit Euler step for the Gray–Scott model.
%
% Inputs:
%   u, v : current state vectors (N x 1)
%   op   : operator struct (contains op.mode, op.L/op.S, and op.grid)
%   p    : parameter struct (Du, Dv, F, k, dt, and optional p.BC / p.BCcache)
%   t    : current time (optional, used only for manufactured solutions)
%
% Outputs:
%   unew, vnew : updated state vectors (N x 1)
%   info       : diagnostics (min/max values and NaN/Inf flag)

if nargin < 5
    t = 0;
end

assert(isequal(size(u), size(v)), "eulerStep: u and v must have the same size.");

% ---- Diffusion contributions ----
[du_diff, dv_diff] = applyDiffusion(u, v, op, p);

% ---- Reaction contributions ----
[du_reac, dv_reac] = reactionGrayScott(u, v, p);

% ---- Optional source terms (MMS verification only) ----
su = zeros(size(u));
sv = su;

if isfield(p, "sourceFcn") && ~isempty(p.sourceFcn)
    x = op.grid.x;
    y = op.grid.y;
    [su2D, sv2D] = p.sourceFcn(x, y, t, p);
    su = su2D(:);
    sv = sv2D(:);
end

% ---- Explicit Euler update ----
unew = u + p.dt * (du_diff + du_reac + su);
vnew = v + p.dt * (dv_diff + dv_reac + sv);

% ---- Enforce boundary conditions (Option A: overwrite after update) ----
if isfield(p, "BC") && isfield(p, "BCcache")
    [unew, vnew] = applyBoundaryConditions(unew, vnew, p, op.grid, p.BCcache);
end

% ---- Diagnostics ----
info.u_min = min(unew);
info.u_max = max(unew);
info.v_min = min(vnew);
info.v_max = max(vnew);
info.hasNaNInf = any(~isfinite(unew)) || any(~isfinite(vnew));

if info.hasNaNInf
    warning("NaN or Inf detected in solution after Euler step.");
end
end