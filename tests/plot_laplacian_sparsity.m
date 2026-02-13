clear; clc; close all;

%PLOT_LAPLACIAN_SPARSITY  Plot and save the sparsity pattern of the 2D Laplacian.
%
%   This script builds the 2D periodic Laplacian operator and visualizes its
%   sparsity pattern using spy(). The figure is saved to:
%       figures/verification/laplacian_sparsity.png
%
%   The script does not depend on the current working directory.

% Locate project folders relative to this script
thisDir = fileparts(mfilename("fullpath"));      % .../tests
rootDir = fullfile(thisDir, "..");              % project root

% Add src (absolute path)
addpath(fullfile(rootDir, "src"));

% Build Laplacian
p = defaultParams();
grid = buildGrid(p);
L = buildLaplacian2D(p, grid);

% Plot sparsity pattern
figure;
spy(L);
title("Sparsity pattern of 2D Laplacian (periodic BC)");

% Ensure output directories exist
figDir = fullfile(rootDir, "figures", "verification");
if ~exist(figDir, "dir")
    mkdir(figDir);
end

% Save figure for report
set(gcf, "PaperPositionMode", "auto");
print(gcf, fullfile(figDir, "laplacian_sparsity.png"), "-dpng", "-r200");

disp("Saved figures/verification/laplacian_sparsity.png");
