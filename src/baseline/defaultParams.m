function p = defaultParams(preset)
%DEFAULTPARAMS  Return default parameters for the Gray–Scott model.
%
%   p = DEFAULTPARAMS() returns the baseline parameter set.
%   p = DEFAULTPARAMS("pearson") returns a preset aligned with
%   Pearson's Gray–Scott pattern experiments.
%
%   Keeping all parameters in one place avoids hidden constants
%   in the code and makes simulation runs reproducible.

if nargin < 1
    preset = "baseline";
else
    preset = string(preset);
end

% Store preset name (normalized form)
p.preset = preset;

% --- Model parameters (Gray–Scott equations) ---
p.Du = 2e-5;    % diffusion coefficient of u
p.Dv = 1e-5;    % diffusion coefficient of v
p.F  = 0.04;    % feed rate
p.k  = 0.06;    % kill rate

% --- Domain / grid parameters ---
p.Lx = 1.0;     % domain length in x
p.Ly = 1.0;     % domain length in y
p.Nx = 128;     % number of grid points in x
p.Ny = 128;     % number of grid points in y

% --- Boundary condition ---
p.bc = "periodic";

% --- Time integration ---
% Smaller dt improves stability for explicit diffusion but increases runtime.
p.dt = 0.1;     % time step size
p.T  = 500.0;   % final simulation time

% --- Initial condition settings ---
p.seed       = 1;     % RNG seed for reproducible perturbations
p.icType     = "baseline"; % "baseline" | "pearson"
p.icPerturb  = 0.02;  % amplitude of random noise (if used)

% --- Output / visualization control ---
p.caseName      = "baseline"; % label used in results folder naming
p.plotField     = "v";        % field shown in live plots ("u" or "v")
p.plotEvery     = 10;         % update live plot every N steps (0 disables)
p.savePngEvery  = 50;        % export PNG every N steps (0 disables)

% --- Preset overrides ---
switch preset
    case "baseline"
        % Default parameters are used as defined above.

    case "pearson"
        % Pearson paper setup: larger domain and finer grid
        p.Lx = 2.5;  p.Ly = 2.5;
        p.Nx = 256;  p.Ny = 256;

        % Initial condition and visualization choices aligned with Pearson
        p.icType     = "pearson";
        p.icPerturb  = 0.01;   % ~1% noise as reported in the paper
        p.plotField  = "u";    % Pearson often visualizes the u-field

        % Keep dt conservative by default for stability
        % p.dt = 1.0;  % Enable only after stability testing

        p.caseName = "pearson";

    otherwise
        error("defaultParams: unknown preset '%s'. Use 'baseline' or 'pearson'.", preset);
end

end
