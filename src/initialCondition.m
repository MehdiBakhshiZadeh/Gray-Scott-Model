function [U, V] = initialCondition(p)
%INITIALCONDITION  Create the initial fields for the Gray–Scott model.
%
%   [U,V] = INITIALCONDITION(p) returns Ny-by-Nx arrays for the initial
%   concentrations U and V.
%
%   Supported initial-condition types:
%     - p.icType = "baseline":
%         A centered square patch is set to U=0.50 and V=0.25.
%     - p.icType = "pearson":
%         A centered 20-by-20 square patch is set to U=0.50 and V=0.25.
%
%   Noise:
%     If p.icPerturb > 0, Gaussian noise with amplitude p.icPerturb is added
%     to both U and V.
%
%   Notes:
%     This function does not set the random seed. If you require fully
%     reproducible noise, call rng(p.seed, ...) before calling this function.

% --- Input validation ---
required = ["Nx","Ny"];
for k = 1:numel(required)
    assert(isfield(p, required(k)), "initialCondition: missing parameter p.%s", required(k));
end
assert(mod(p.Nx,1) == 0 && mod(p.Ny,1) == 0, "initialCondition: Nx and Ny must be integers.");
assert(p.Nx > 0 && p.Ny > 0, "initialCondition: Nx and Ny must be positive.");

if ~isfield(p, "icPerturb") || isempty(p.icPerturb)
    error("initialCondition: missing parameter p.icPerturb.");
end
assert(p.icPerturb >= 0, "initialCondition: icPerturb must be nonnegative.");

Ny = p.Ny;
Nx = p.Nx;

% Baseline state
U = ones(Ny, Nx);
V = zeros(Ny, Nx);

% Center index (used for both presets)
cx = floor(Nx/2);
cy = floor(Ny/2);

icType = "baseline";
if isfield(p, "icType")
    icType = string(p.icType);
end

switch icType
    case "baseline"
        % Centered square patch (approx. 21x21 when r=10)
        r = 10;
        ix = max(1, cx-r) : min(Nx, cx+r);
        iy = max(1, cy-r) : min(Ny, cy+r);

        U(iy, ix) = 0.50;
        V(iy, ix) = 0.25;

    case "pearson"
        % Pearson-style 20x20 square patch at the center
        patch = 20;
        ix = (cx - patch/2) : (cx + patch/2 - 1);
        iy = (cy - patch/2) : (cy + patch/2 - 1);

        % Clamp indices if the grid is small
        ix = max(1, ix(1)) : min(Nx, ix(end));
        iy = max(1, iy(1)) : min(Ny, iy(end));

        U(iy, ix) = 0.50;
        V(iy, ix) = 0.25;

    otherwise
        error("initialCondition: unknown icType '%s'.", icType);
end

% Add small random perturbation
if p.icPerturb > 0
    U = U + p.icPerturb * randn(Ny, Nx);
    V = V + p.icPerturb * randn(Ny, Nx);
end

end
