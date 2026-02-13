function [L2, Linf] = fieldNorms2D(A, hx, hy)
%FIELDNORMS2D Weighted discrete norms for a 2D field A on a uniform grid.
%
%   L2   = sqrt( sum_{i,j} A(i,j)^2 * hx * hy )
%   Linf = max_{i,j} |A(i,j)|
%
%   If hx and hy are not provided (or are empty), the defaults hx=1 and hy=1
%   are used, which makes the L2 norm effectively unweighted:
%       L2 = sqrt( sum_{i,j} A(i,j)^2 )

if nargin < 2 || isempty(hx), hx = 1; end
if nargin < 3 || isempty(hy), hy = 1; end

Linf = max(abs(A(:)));

% Weighted L2 (area-weighted)
L2 = sqrt(sum(A(:).^2) * hx * hy);
end
