function exportFrame(u, v, p, t, outDir)
%EXPORTFRAME  Save a snapshot of the Gray–Scott field as a PNG image.
%
%   exportFrame(u, v, p, t, outDir) writes one PNG file to outDir.
%   This function is intended for offline output only and does not affect
%   the numerical simulation state.
%
%   The exported field is selected using p.plotField ("u" or "v").
%   A new figure is not created for each call; an internal invisible
%   figure is reused for efficiency.
%
% Inputs:
%   u, v   : vectorized state fields (N x 1), where N = p.Nx * p.Ny
%   p      : parameter struct (requires p.Nx, p.Ny, optionally p.plotField)
%   t      : current simulation time
%   outDir : directory where the PNG file will be written

% --- Input validation ---
required = ["Nx","Ny"];
for k = 1:numel(required)
    assert(isfield(p, required(k)), "exportFrame: missing parameter p.%s", required(k));
end
assert(mod(p.Nx,1) == 0 && mod(p.Ny,1) == 0, "exportFrame: Nx and Ny must be integers.");
assert(p.Nx > 0 && p.Ny > 0, "exportFrame: Nx and Ny must be positive.");

N = p.Nx * p.Ny;
assert(numel(u) == N, "exportFrame: numel(u) must equal Nx*Ny.");
assert(numel(v) == N, "exportFrame: numel(v) must equal Nx*Ny.");

assert(isfolder(outDir), "exportFrame: output directory does not exist.");

U = reshape(u, p.Ny, p.Nx);
V = reshape(v, p.Ny, p.Nx);

% Choose which field to export (default: v)
field = "v";
if isfield(p, "plotField")
    field = lower(string(p.plotField));
end

switch field
    case "u"
        data2D = U;
        fieldName = "u";
    case "v"
        data2D = V;
        fieldName = "v";
    otherwise
        error("exportFrame: unknown p.plotField='%s'. Use 'u' or 'v'.", field);
end

% --- Reuse an invisible figure for efficiency ---
persistent hFig hAx hIm hCb

needCreate = isempty(hFig) || ~isgraphics(hFig) || ...
             isempty(hAx)  || ~isgraphics(hAx)  || ...
             isempty(hIm)  || ~isgraphics(hIm);

if needCreate
    hFig = figure("Visible","off", "Color","w");
    hAx  = axes(hFig);

    hIm = imagesc(hAx, data2D);
    axis(hAx, "image");
    hCb = colorbar(hAx); %#ok<NASGU>
else
    set(hIm, "CData", data2D);
end

title(hAx, sprintf("Gray-Scott: %s-field | t = %.2f | F=%.3f, k=%.3f", ...
    fieldName, t, p.F, p.k), "Interpreter", "none");

drawnow;

filename = sprintf("%s_t_%07.2f.png", fieldName, t);
filepath = fullfile(outDir, filename);

% More robust for PNG than saveas (especially with invisible figures)
print(hFig, filepath, "-dpng", "-r300");
end
