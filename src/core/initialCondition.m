function [U, V] = initialCondition(p)
%INITIALCONDITION  Create the initial fields for the Gray–Scott model.
%
%   [U,V] = INITIALCONDITION(p) returns Ny-by-Nx arrays for the initial
%   concentrations U and V.
%
%   Supported initial-condition types:
%     - p.icType = "baseline":
%         A centered 10-by-10 square patch is set to U=0.50 and V=0.25.
%     - p.icType = "pearson":
%         A centered 20-by-20 square patch is set to U=0.50 and V=0.25.
%
%   Noise:
%     If p.icPerturb > 0, Gaussian noise with amplitude p.icPerturb is added
%     to both U and V.
%
%   Notes:
%     This function generates noise using randn. Reproducibility depends on
%     the current RNG state (the GrayScottModel constructor/reset sets rng(p.seed)
%     when p.seed is provided).

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

% Center index
cx = floor(Nx/2);
cy = floor(Ny/2);

icType = "baseline";
if isfield(p, "icType")
    icType = lower(strtrim(string(p.icType)));
end

switch icType
    case "baseline"
        patch = 10;   % 10x10 square

    case "pearson"
        patch = 20;   % 20x20 square

    otherwise
        error("initialCondition: unknown icType '%s'.", icType);
end

% Centered square patch
ix = (cx - patch/2) : (cx + patch/2 - 1);
iy = (cy - patch/2) : (cy + patch/2 - 1);

% Clamp indices if the grid is small
ix = max(1, ix(1)) : min(Nx, ix(end));
iy = max(1, iy(1)) : min(Ny, iy(end));

U(iy, ix) = 0.50;
V(iy, ix) = 0.25;

% Add small random perturbation
if p.icPerturb > 0
    U = U + p.icPerturb * randn(Ny, Nx);
    V = V + p.icPerturb * randn(Ny, Nx);
end

end
