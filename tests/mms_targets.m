function m = mms_targets()
%MMS_TARGETS  Exact manufactured solutions for MMS verification.
%
%   m = MMS_TARGETS() returns function handles for an exact, smooth,
%   periodic manufactured solution and its analytic derivatives.
%
%   The solution is defined on a periodic domain [0,1) x [0,1).
%   Inputs x, y are expected to be arrays (e.g., from meshgrid), and
%   t may be a scalar or array. All operations are elementwise.
%
%   Returned fields:
%     m.u(x,y,t)     : exact solution u
%     m.v(x,y,t)     : exact solution v
%     m.ut(x,y,t)    : time derivative of u
%     m.vt(x,y,t)    : time derivative of v
%     m.lapU(x,y,t)  : Laplacian of u
%     m.lapV(x,y,t)  : Laplacian of v

% --- Chosen wavenumbers (integers ensure periodicity) ---
kx = 2;
ky = 3;
omega = 1.0;   % time frequency

% --- Mean values (keep fields away from zero) ---
u0 = 1.0;
v0 = 0.5;

% --- Amplitudes ---
A = 0.10;
B = 0.10;

% Phase function (elementwise, periodic on [0,1)x[0,1))
phi = @(x,y,t) cos(2*pi*(kx.*x + ky.*y) + omega.*t);

% --- Exact solutions ---
m.u = @(x,y,t) u0 + A .* phi(x,y,t);
m.v = @(x,y,t) v0 + B .* phi(x,y,t);

% --- Time derivatives ---
m.ut = @(x,y,t) -A .* omega .* sin(2*pi*(kx.*x + ky.*y) + omega.*t);
m.vt = @(x,y,t) -B .* omega .* sin(2*pi*(kx.*x + ky.*y) + omega.*t);

% --- Laplacians ---
lapfactor = -(2*pi)^2 * (kx^2 + ky^2);
m.lapU = @(x,y,t) A .* lapfactor .* phi(x,y,t);
m.lapV = @(x,y,t) B .* lapfactor .* phi(x,y,t);

% --- Parameters for metadata / debugging ---
m.info = struct();
m.info.kx = kx;
m.info.ky = ky;
m.info.omega = omega;
m.info.A = A;
m.info.B = B;
m.info.u0 = u0;
m.info.v0 = v0;
end
