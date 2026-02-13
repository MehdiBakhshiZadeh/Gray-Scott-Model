function [fu, fv] = reactionGrayScott(u, v, p)
%REACTIONGRAYSCOTT Reaction terms of the Gray-Scott model.
%
% Inputs:
%   u, v : state vectors (size N x 1)
%   p    : parameter struct
%
% Outputs:
%   fu, fv : reaction contributions (size N x 1)

% Basic input checks to avoid silent dimension or parameter errors.
assert(isequal(size(u), size(v)), "reactionGrayScott: u and v must be same size.");
assert(isfield(p,"F") && isfield(p,"k"), "reactionGrayScott: p must contain F and k.");

% Reaction kinetics
uv2 = u .* (v.^2);

fu = -uv2 + p.F * (1-u);
fv = uv2 - (p.F + p.k) * v;
end