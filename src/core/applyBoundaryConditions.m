function [u, v] = applyBoundaryConditions(u, v, p, grid, BCcache)
%APPLYBOUNDARYCONDITIONS Enforce boundary conditions on the 2D state fields.
%
%   [u, v] = applyBoundaryConditions(u, v, p, BCcache) overwrites boundary
%   grid nodes according to p.BC (per-side periodic/dirichlet/neumann).
%
%   Boundary conventions (with U = reshape(u, Ny, Nx)):
%     left   : U(:,1)
%     right  : U(:,Nx)
%     bottom : U(1,:)
%     top    : U(Ny,:)
%
%   Dirichlet:
%     - "constant": boundary values set to user constants
%     - "initial" : boundary values frozen from the initial condition cache
%
%   Neumann (forward/backward one-sided difference, user value is du/dn):
%     left   : (U(:,2)  - U(:,1))  / hx = g_left
%     right  : (U(:,Nx) - U(:,Nx-1))/ hx = g_right
%     bottom : (U(2,:)  - U(1,:))  / hy = g_bottom
%     top    : (U(Ny,:) - U(Ny-1,:))/ hy = g_top
%
%   Periodic:
%     For Option A enforcement, periodic does not require overwriting
%     because the diffusion operator is already periodic. This function
%     leaves periodic sides unchanged.
%
% Inputs:
%   u, v     : state vectors (Nx*Ny x 1)
%   p        : parameter struct (requires Nx, Ny and grid spacing in p.grid)
%   BCcache  : struct of cached initial boundary values for "initial" mode
%
% Outputs:
%   u, v     : updated state vectors with BC enforced

% --- Basic validation ---
Nx = p.Nx; Ny = p.Ny;
assert(numel(u) == Nx*Ny && numel(v) == Nx*Ny, ...
    "applyBoundaryConditions: state size must be Nx*Ny.");
assert(isfield(p, "BC"), "applyBoundaryConditions: missing p.BC.");
assert(isfield(grid, "hx") && isfield(grid, "hy"), ...
    "applyBoundaryConditions: missing grid.hx / grid.hy.");

hx = grid.hx;
hy = grid.hy;

U = reshape(u, Ny, Nx);
V = reshape(v, Ny, Nx);

% Helper to fetch per-side struct safely
sides = ["left","right","bottom","top"];
for i = 1:numel(sides)
    side = sides(i);
    assert(isfield(p.BC, side), "applyBoundaryConditions: missing p.BC.%s", side);
end

% --- Dirichlet ---
applyDirichlet("left",   @(A,val) setLeft(A,val),   @(C) C.left);
applyDirichlet("right",  @(A,val) setRight(A,val),  @(C) C.right);
applyDirichlet("bottom", @(A,val) setBottom(A,val), @(C) C.bottom);
applyDirichlet("top",    @(A,val) setTop(A,val),    @(C) C.top);

% --- Neumann (du/dn = g) ---
applyNeumann("left",   @(A,g) setLeft(A,  A(:,2)   - hx*g));
applyNeumann("right",  @(A,g) setRight(A, A(:,Nx-1)+ hx*g));
applyNeumann("bottom", @(A,g) setBottom(A, A(2,:)  - hy*g));
applyNeumann("top",    @(A,g) setTop(A,    A(Ny-1,:)+ hy*g));

u = U(:);
v = V(:);

% ================== nested helpers ==================

    function applyDirichlet(sideName, setter, cacheGetter)
        bc = p.BC.(sideName);
        if string(bc.type) ~= "dirichlet"
            return;
        end

        mode = string(bc.dirichletMode);
        if mode == "constant"
            U = setter(U, bc.dirichletValue.u);
            V = setter(V, bc.dirichletValue.v);
        elseif mode == "initial"
            assert(nargin >= 4 && ~isempty(BCcache), ...
                "applyBoundaryConditions: BCcache required for dirichletMode='initial'.");
            C = cacheGetter(BCcache);
            U = setter(U, C.u);
            V = setter(V, C.v);
        else
            error("applyBoundaryConditions: unknown dirichletMode '%s' on %s.", mode, sideName);
        end
    end

    function applyNeumann(sideName, updater)
        bc = p.BC.(sideName);
        if string(bc.type) ~= "neumann"
            return;
        end

        % For now only "constant" is supported, but stored explicitly.
        mode = string(bc.neumannMode);
        if mode ~= "constant"
            error("applyBoundaryConditions: unknown neumannMode '%s' on %s.", mode, sideName);
        end

        U = updater(U, bc.neumannValue.u);
        V = updater(V, bc.neumannValue.v);
    end

    function A = setLeft(A, val),   A(:,1)  = val; end
    function A = setRight(A, val),  A(:,Nx) = val; end
    function A = setBottom(A, val), A(1,:)  = val; end
    function A = setTop(A, val),    A(Ny,:) = val; end

end