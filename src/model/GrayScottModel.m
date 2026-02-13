classdef GrayScottModel < handle
    %GRAYSCOTTMODEL  Gray–Scott reaction–diffusion simulator (2D, periodic).
    %
    %   This class stores the model parameters, grid, Laplacian operator,
    %   and the current solution state (u,v). It provides:
    %     - step()   : advance one explicit Euler time step
    %     - reset()  : restore the initial condition and set t=0
    %     - getFields2D() : return u,v as Ny-by-Nx arrays for plotting
    %
    %   Typical use:
    %       p = defaultParams();
    %       model = GrayScottModel(p);
    %       for i = 1:1000
    %           info = model.step();
    %           if info.hasNaNInf, error("Solution blew up."); end
    %       end
    %       [U,V] = model.getFields2D();

    properties
        p       % parameter struct
        grid    % grid struct (hx, hy, x, y)

        u       % state vector (Nx*Ny x 1)
        v       % state vector (Nx*Ny x 1)

        t = 0.0 % current time
        n = 0   % current step index
    end

    properties (Constant)
        % LaTeX form of the Gray–Scott equations (for reports / documentation)
        eqLatex = "u_t = D_u \nabla^2 u - u v^2 + F(1-u),\quad " + ...
                  "v_t = D_v \nabla^2 v + u v^2 - (F+k)v";
    end

    methods
        function obj = GrayScottModel(p, grid)
            % Constructor: store parameters and build grid + operator.
            if nargin < 1 || isempty(p)
                p = defaultParams();
            end
            GrayScottModel.validateParams(p);
            obj.p = p;

            if nargin < 2 || isempty(grid)
                error("GrayScottModel: grid must be provided. Use buildGrid(p).");
            end
            obj.grid = grid;

            obj.initState();

        end


        function reset(obj)
            %RESET  Restore the initial condition and set t=0.

            obj.initState();
        end

        function [U, V] = getFields2D(obj)
            %GETFIELDS2D  Return U,V as Ny-by-Nx arrays.

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
            % Initialize (u,v) from the initial condition and reset time counters.
            [U0, V0] = initialCondition(obj.p);
            obj.u = U0(:);
            obj.v = V0(:);

            obj.t = 0.0;
            obj.n = 0;
        end
    end

    methods (Static, Access = private)
        function validateParams(p)
            % Minimal parameter validation to avoid confusing runtime errors.
            required = ["Nx","Ny","dt","T","Du","Dv","F","k"];
            for i = 1:numel(required)
                assert(isfield(p, required(i)), "Missing parameter: p.%s", required(i));
            end
            assert(p.Nx > 0 && p.Ny > 0, "Nx and Ny must be positive.");
            assert(p.dt > 0 && p.T >= 0, "dt must be positive and T must be nonnegative.");
            assert(p.Du >= 0 && p.Dv >= 0, "Diffusion coefficients must be nonnegative.");
        end
    end
end
