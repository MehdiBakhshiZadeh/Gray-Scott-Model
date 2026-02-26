% run_gui.m
% Launches the Gray-Scott GUI safely from any working directory.

clc;

% Locate this file on disk
thisFile = mfilename("fullpath");
if strlength(thisFile) == 0
    thisFile = matlab.desktop.editor.getActiveFilename();
end
assert(strlength(thisFile) > 0, "Cannot locate run_gui.m on disk.");

% Project root = folder that contains this launcher file
rootDir = fileparts(thisFile);

% Add required folders temporarily
addpath(genpath(fullfile(rootDir, "src")));
addpath(fullfile(rootDir, "scripts"));

% Close any stale GUI instance from a previous run
old = findall(0, "Type", "figure", "Tag", "GrayScottControllerUI");
if ~isempty(old)
    delete(old);
end

% Launch GUI
GrayScottController;