function [L2, Linf] = fieldNormsVec(a, hx, hy)
%FIELDNORMSVEC Norms for a vectorized 2D field.
%
%   NOTE (naming): This function is intended for vectorized 2D fields and
%   is effectively a "fieldNormsVector" utility. The name is kept as
%   FIELDNORMSVEC for backward compatibility with existing calls.
%
%   [L2, Linf] = fieldNormsVec(a) returns:
%     - L2   : unweighted Euclidean 2-norm of a(:)
%     - Linf : max-norm (infinity norm) of a(:)
%
%   [L2, Linf] = fieldNormsVec(a, hx, hy) returns a weighted L2 norm that
%   corresponds to a discrete approximation of the continuous L2 norm on a
%   uniform grid with spacings hx and hy:
%
%       L2 = sqrt( sum_{i} |a_i|^2 * hx * hy )
%
%   Inputs:
%     a  : vector (or array) of field values (vectorized field recommended)
%     hx : grid spacing in x (optional)
%     hy : grid spacing in y (optional)

Linf = max(abs(a));

if nargin < 2 || isempty(hx) || nargin < 3 || isempty(hy)
    L2 = norm(a, 2);
else
    L2 = sqrt(sum(a(:).^2) * hx * hy);
end
end
