function out = saveFigWithMeta(figHandle, filepath, meta)
%SAVEFIGWITHMETA Save a figure and a sidecar metadata .mat file.
%
% out.imagePath : path to saved image
% out.metaPath  : path to saved metadata file
%
% This version:
% - forces report-friendly white styling
% - styles axes, labels, titles, legends, and colorbars
% - keeps legends optional
% - removes empty legends safely

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

% Force filepath to char
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
set(figHandle, ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'InvertHardcopy', 'off');

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

% -----------------------------
% Style normal axes
% -----------------------------
axs = findall(figHandle, 'Type', 'axes');

for k = 1:numel(axs)
    ax = axs(k);

    % Skip legend/colorbar axes in older MATLAB versions if they appear here
    axTag = '';
    try
        axTag = get(ax, 'Tag');
    catch
    end

    if strcmpi(axTag, 'legend') || strcmpi(axTag, 'Colorbar')
        continue;
    end

    set(ax, ...
        'Color', 'w', ...
        'XColor', 'k', ...
        'YColor', 'k', ...
        'GridColor', [0.75 0.75 0.75], ...
        'MinorGridColor', [0.88 0.88 0.88], ...
        'GridAlpha', 0.35, ...
        'MinorGridAlpha', 0.25, ...
        'FontSize', 16, ...
        'LineWidth', 1.2, ...
        'Box', 'on');

    % Axis labels
    xlab = get(ax, 'XLabel');
    if ~isempty(xlab) && isgraphics(xlab)
        set(xlab, 'FontSize', 18, 'Color', 'k');
    end

    ylab = get(ax, 'YLabel');
    if ~isempty(ylab) && isgraphics(ylab)
        set(ylab, 'FontSize', 18, 'Color', 'k');
    end

    % Title
    ttl = get(ax, 'Title');
    if ~isempty(ttl) && isgraphics(ttl)
        set(ttl, 'FontSize', 20, 'Color', 'k');
    end
end

% -----------------------------
% Style colorbars explicitly
% -----------------------------
cbs = findall(figHandle, 'Type', 'ColorBar');
for k = 1:numel(cbs)
    cb = cbs(k);
    set(cb, ...
        'Color', 'k', ...
        'FontSize', 16, ...
        'LineWidth', 1.0);
end

% -----------------------------
% Style legends if they exist
% Remove empty legends safely
% -----------------------------
lgd = findall(figHandle, 'Type', 'legend');
for k = 1:numel(lgd)
    L = lgd(k);

    keepLegend = true;

    % Check whether the legend actually has useful entries
    try
        strs = L.String;

        if ischar(strs)
            strs = cellstr(strs);
        elseif isstring(strs)
            strs = cellstr(strs);
        end

        if isempty(strs)
            keepLegend = false;
        else
            tmp = strings(numel(strs),1);
            for j = 1:numel(strs)
                tmp(j) = string(strs{j});
            end
            tmp = strtrim(tmp);
            if all(tmp == "")
                keepLegend = false;
            end
        end
    catch
        % If legend properties are inaccessible for some reason,
        % keep it and just style it.
        keepLegend = true;
    end

    if ~keepLegend
        delete(L);
        continue;
    end

    % Style only existing non-empty legends
    try
        set(L, ...
            'FontSize', 16, ...
            'TextColor', 'k', ...
            'Color', 'w', ...
            'EdgeColor', [0.2 0.2 0.2]);
    catch
        % Fallback for older MATLAB versions
        try
            set(L, ...
                'FontSize', 16, ...
                'Color', 'w', ...
                'EdgeColor', [0.2 0.2 0.2]);
        catch
        end
    end
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