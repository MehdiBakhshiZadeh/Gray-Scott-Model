function y = applyFullLaplacian(x, Lfull)
%APPLYFULLLAPLACIAN Apply dense (full) Laplacian operator.
%
%   y = applyFullLaplacian(x, Lfull) computes y = Lfull * x.
%
% Inputs:
%   x     : state vector (N x 1)
%   Lfull : dense Laplacian matrix (N x N), double
%
% Output:
%   y     : (N x 1)

assert(isvector(x) && ~isempty(x), "applyFullLaplacian: x must be a nonempty vector.");
assert(ismatrix(Lfull) && ~issparse(Lfull), "applyFullLaplacian: Lfull must be a full (dense) matrix.");

x = x(:);
N = numel(x);
assert(size(Lfull,1) == N && size(Lfull,2) == N, ...
    "applyFullLaplacian: size mismatch. Lfull is %dx%d but x has length %d.", ...
    size(Lfull,1), size(Lfull,2), N);

y = Lfull * x;
end
