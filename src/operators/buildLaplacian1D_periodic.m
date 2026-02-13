function L1 = buildLaplacian1D_periodic(N)
%BUILDLAPLACIAN1D_PERIODIC  1D second-derivative matrix with periodic BC.
%
%   L1 = BUILDLAPLACIAN1D_PERIODIC(N) returns the N-by-N sparse matrix that
%   approximates the 1D operator d^2/dx^2 using the standard centered stencil:
%       u_xx(i) ? u(i-1) - 2u(i) + u(i+1)
%   with periodic wrap-around at the boundaries.
%
%   Important:
%     This function returns the *unscaled* stencil matrix.
%     To approximate the Laplacian on a grid with spacing hx, use:
%         (1/hx^2) * L1

% --- Input validation ---
assert(isscalar(N) && isfinite(N), "buildLaplacian1D_periodic: N must be a finite scalar.");
assert(mod(N,1) == 0 && N >= 3, "buildLaplacian1D_periodic: N must be an integer >= 3.");

e = ones(N, 1);

% Main tridiagonal part
L1 = spdiags([e, -2*e, e], [-1, 0, 1], N, N);

% Periodic wrap-around terms
L1(1, N) = 1;
L1(N, 1) = 1;
end
