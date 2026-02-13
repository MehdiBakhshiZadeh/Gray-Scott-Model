function out = saveFigWithMeta(figHandle, filepath, meta)
%SAVEFIGWITHMETA Save a figure and a sidecar metadata .mat file.
% Compatible with MATLAB R2018b.
%
% out.imagePath : path to saved image
% out.metaPath  : path to saved metadata file

% -----------------------------
% Input handling
% -----------------------------
if nargin < 1 || isempty(figHandle)
    figHandle = gcf;
end
if nargin < 2
    error('filepath is required');
end
if nargin < 3 || isempty(meta)
    meta = struct();
end

% Force filepath to char (critical for R2018b)
if isstring(filepath)
    filepath = char(filepath);
end

% -----------------------------
% Ensure output directory exists
% -----------------------------
folder = fileparts(filepath);
if ~isempty(folder) && ~exist(folder, 'dir')
    mkdir(folder);
end

% -----------------------------
% Standardize figure appearance
% -----------------------------
set(figHandle, 'Color', 'w', 'Units', 'pixels');

% --- Dynamic size based on screen ---
screenSize = get(0, 'ScreenSize');   % [left bottom width height]
figWidth  = round(0.85 * screenSize(3));
figHeight = round(0.75 * screenSize(4));

set(figHandle, 'Position', [ ...
    round(0.075 * screenSize(3)), ...
    round(0.10  * screenSize(4)), ...
    figWidth, ...
    figHeight]);

% Force the window onto the visible screen
movegui(figHandle, 'onscreen');

axs = findall(figHandle, 'Type', 'axes');
for k = 1:numel(axs)
    ax = axs(k);

    % Axis ticks
    set(ax, ...
        'FontSize', 16, ...
        'LineWidth', 1.2, ...
        'Box', 'on');

    % Axis labels
    if ~isempty(get(ax, 'XLabel'))
        set(get(ax, 'XLabel'), 'FontSize', 18);
    end
    if ~isempty(get(ax, 'YLabel'))
        set(get(ax, 'YLabel'), 'FontSize', 18);
    end

    % Title
    if ~isempty(get(ax, 'Title'))
        set(get(ax, 'Title'), 'FontSize', 20);
    end
end

% Legends (if any)
lgd = findall(figHandle, 'Type', 'legend');
for k = 1:numel(lgd)
    set(lgd(k), 'FontSize', 16);
end

drawnow;  % finalize layout before export

% -----------------------------
% Save figure (PNG, high quality)
% -----------------------------
print(figHandle, filepath, '-dpng', '-r300');

% -----------------------------
% Save metadata sidecar
% -----------------------------
[p, name, ~] = fileparts(filepath);
metaPath = fullfile(p, [name '_meta.mat']);

if isstring(metaPath)
    metaPath = char(metaPath);
end

meta.savedAt = datestr(now, 31);  % ISO-like timestamp
meta.matlabVersion = version;
meta.platform = computer;

save(metaPath, 'meta');

% -----------------------------
% Output struct
% -----------------------------
out = struct();
out.imagePath = filepath;
out.metaPath  = metaPath;
end
