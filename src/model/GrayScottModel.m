classdef GrayScottModel < handle
    %GRAYSCOTTMODEL  Gray–Scott reaction–diffusion simulator (2D, periodic).
    %
    % Stores parameters, grid, diffusion operator, and state (u,v).
    % Provides:
    %   - step()        : advance one explicit Euler step
    %   - reset()       : restore initial condition and reset t,n
    %   - getFields2D() : return U,V as Ny-by-Nx arrays

    properties
        p       % parameter struct
        grid    % grid struct (hx, hy, x, y)

        op      % operator struct (mode, L or S, and op.grid)
        u       % state vector (Nx*Ny x 1)
        v       % state vector (Nx*Ny x 1)

        t = 0.0 % current time
        n = 0   % current step index
    end

    properties (Constant)
        eqLatex = "u_t = D_u \nabla^2 u - u v^2 + F(1-u),\quad " + ...
                  "v_t = D_v \nabla^2 v + u v^2 - (F+k)v";
    end

    methods
        function obj = GrayScottModel(p, grid)
            % Constructor: store parameters, grid, build operator, init state.
            if nargin < 1 || isempty(p)
                p = defaultParams();
            end
            GrayScottModel.validateParams(p);
            obj.p = p;

            if nargin < 2 || isempty(grid)
                error("GrayScottModel: grid must be provided. Use buildGrid(p).");
            end
            obj.grid = grid;

            % Build diffusion operator once
            obj.op = GrayScottModel.buildOperator(obj.p, obj.grid);

            % Initialize state
            obj.initState();
        end

        function reset(obj)
            % Reset initial condition and counters (keeps operator and params).
            obj.initState();
        end

        function info = step(obj)
            %STEP  Advance one explicit Euler step.
            %
            % Returns info struct with at least info.hasNaNInf.

            [obj.u, obj.v, info] = eulerStep(obj.u, obj.v, obj.op, obj.p, obj.t);

            obj.t = obj.t + obj.p.dt;
            obj.n = obj.n + 1;
        end

        function [U, V] = getFields2D(obj)
            Ny = obj.p.Ny;
            Nx = obj.p.Nx;

            assert(numel(obj.u) == Nx*Ny, "State vector size does not match Nx*Ny.");
            assert(numel(obj.v) == Nx*Ny, "State vector size does not match Nx*Ny.");

            U = reshape(obj.u, Ny, Nx);
            V = reshape(obj.v, Ny, Nx);
        end
    end

    methods (Access = private)
        function initState(obj)
            [U0, V0] = initialCondition(obj.p);
            obj.u = U0(:);
            obj.v = V0(:);

            % Cache initial boundary values for Dirichlet "initial" mode
            obj.p = GrayScottModel.buildBCcache(obj.p, U0, V0);

            obj.t = 0.0;
            obj.n = 0;
        end
    end

    methods (Static, Access = private)
        function validateParams(p)
            required = ["Nx","Ny","dt","T","Du","Dv","F","k","diffusionMode"];
            for i = 1:numel(required)
                assert(isfield(p, required(i)), "Missing parameter: p.%s", required(i));
            end
            assert(p.Nx > 0 && p.Ny > 0, "Nx and Ny must be positive.");
            assert(p.dt > 0 && p.T >= 0, "dt must be positive and T must be nonnegative.");
            assert(p.Du >= 0 && p.Dv >= 0, "Diffusion coefficients must be nonnegative.");
        end

        function p = buildBCcache(p, U0, V0)
            %BUILDBCCACHE Cache initial boundary values for Dirichlet "initial" mode.
            %
            % Stores boundary values from the initial condition so that boundaries can
            % remain fixed even as the interior evolves.
        
            if ~isfield(p, "BC")
                return;
            end
        
            Nx = p.Nx; Ny = p.Ny;
            assert(all(size(U0) == [Ny, Nx]) && all(size(V0) == [Ny, Nx]), ...
                "buildBCcache: U0/V0 must be Ny-by-Nx.");
        
            C = struct();
        
            C.left.u   = U0(:,1);   C.left.v   = V0(:,1);
            C.right.u  = U0(:,Nx);  C.right.v  = V0(:,Nx);
            C.bottom.u = U0(1,:);   C.bottom.v = V0(1,:);
            C.top.u    = U0(Ny,:);  C.top.v    = V0(Ny,:);
        
            p.BCcache = C;
        end

        function op = buildOperator(p, grid)
            %BUILDOPERATOR Create operator struct op used by eulerStep/applyDiffusion.
            %
            % Expected by applyDiffusion:
            %   op.mode = "stencil" | "full" | "matrix"
            %   op.L    = Laplacian matrix (sparse or full) or []
            %   op.S    = stencil struct or []
            % Also used by eulerStep MMS:
            %   op.grid.x, op.grid.y

            mode = string(p.diffusionMode);

            op = struct();
            op.mode = mode;
            op.grid = grid;

            % ---- IMPORTANT ----
            % You must have one of the operator builders below in your project.
            % If your function names differ, tell me and we'll adapt.
            %

            % Explicit builders (common pattern)
            switch mode
                case "matrix"
                    if exist('makeOperator_sparse','file') ~= 2
                        error("buildOperator: makeOperator_sparse.m not found on path.");
                    end
                    op = makeOperator_sparse(p, grid);
                case "full"
                    if exist('makeOperator_dense','file') ~= 2
                        error("buildOperator: makeOperator_dense.m not found on path.");
                    end
                    op = makeOperator_dense(p, grid);
                case "stencil"
                    if exist('makeOperator_stencil','file') ~= 2
                        error("buildOperator: makeOperator_stencil.m not found on path.");
                    end
                    op = makeOperator_stencil(p, grid);
                otherwise
                    error("buildOperator: unknown diffusionMode '%s'", mode);
            end

            % Ensure required fields exist
            if ~isfield(op,'mode'), op.mode = mode; end
            if ~isfield(op,'grid'), op.grid = grid; end
            if ~isfield(op,'L'), op.L = []; end
            if ~isfield(op,'S'), op.S = []; end
        end
    end
end
