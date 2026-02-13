function plotState(u, v, p, t)
%PLOTSTATE  Live visualization of the Gray–Scott fields (plotting only).
%
%   plotState(u, v, p, t) displays the selected field as an image.
%   This function is used only for visualization; it does not affect the
%   numerical solution.
%
% Inputs:
%   u, v : vectorized fields (N x 1), where N = p.Nx * p.Ny
%   p    : parameter struct (requires p.Nx, p.Ny; uses p.plotField if present)
%   t    : current time

% --- Input validation (clear error messages) ---
required = ["Nx","Ny"];
for k = 1:numel(required)
    assert(isfield(p, required(k)), "plotState: missing parameter p.%s", required(k));
end
assert(mod(p.Nx,1) == 0 && mod(p.Ny,1) == 0, "plotState: Nx and Ny must be integers.");
assert(p.Nx > 0 && p.Ny > 0, "plotState: Nx and Ny must be positive.");

N = p.Nx * p.Ny;
assert(numel(u) == N, "plotState: numel(u) must equal Nx*Ny.");
assert(numel(v) == N, "plotState: numel(v) must equal Nx*Ny.");

U = reshape(u, p.Ny, p.Nx);
V = reshape(v, p.Ny, p.Nx);

% Choose which field to plot (default: v)
field = "v";
if isfield(p, "plotField")
    field = string(p.plotField);
end

switch field
    case "u"
        data2D = U;
        fieldName = "u";
    case "v"
        data2D = V;
        fieldName = "v";
    otherwise
        error("plotState: unknown p.plotField='%s'. Use 'u' or 'v'.", field);
end

% --- Faster plotting: reuse graphics objects ---
persistent hIm hCb
if isempty(hIm) || ~isvalid(hIm)
    hIm = imagesc(data2D);
    axis image;
    hCb = colorbar;
else
    set(hIm, "CData", data2D);
end

% Avoid recreating colorbar each call
if isempty(hCb) || ~isvalid(hCb)
    hCb = colorbar;
end

title(sprintf("Gray-Scott: %s-field | t = %.2f | F=%.3f, k=%.3f", ...
    fieldName, t, p.F, p.k));

drawnow;
end
