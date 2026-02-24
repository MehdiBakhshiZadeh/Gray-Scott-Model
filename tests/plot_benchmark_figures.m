function plot_benchmark_figures(csvPath, outDir)
%PLOT_BENCHMARK_FIGURES  Generate all benchmark plots
%
% R2018b compatible.

% --- Project root (one level above /tests) ---
thisFile = mfilename("fullpath");
testsDir = fileparts(thisFile);
rootDir  = fileparts(testsDir);

if nargin < 1 || isempty(csvPath)
    error("Provide path to benchmark_table.csv");
end
if nargin < 2 || isempty(outDir)
    outDir = "figures";
end

csvPath = char(csvPath);
outDir  = char(outDir);
csvPath = strrep(csvPath, "/", filesep);

if exist(csvPath,"file") == 2
    csvPathAbs = csvPath;
else
    csvPathAbs = fullfile(rootDir, csvPath);
end
assert(exist(csvPathAbs,"file")==2, "CSV not found: %s", csvPathAbs);

if exist(outDir,"dir") ~= 7
    outDir = fullfile(rootDir, outDir);
end
if exist(outDir,"dir") ~= 7
    mkdir(outDir);
end

T = readtable(csvPathAbs);

make_compute_loglog(T, outDir);
make_render_loglog(T, outDir);
make_memory_loglog(T, outDir);

fprintf("Saved figures to: %s\n", outDir);
end

function make_compute_loglog(T, outDir)

modeWanted    = "compute_off";
solversWanted = ["matrix","stencil","full"];

solverCol = string(T.diffusionMode);
modeCol   = string(T.mode);

Tc = T(modeCol == modeWanted,:);
assert(~isempty(Tc), "No compute_off data found.");

Nx = double(Tc.Nx);
Ny = double(Tc.Ny);
N  = Nx .* Ny;

Nu = sort(unique(N));
gridLabels = arrayfun(@(n) sprintf("%d^2", round(sqrt(n))), Nu, "UniformOutput", false);

fig = figure("Visible","off","Position",[100 100 900 520]);
ax  = axes("Parent",fig);
hold(ax,"on"); grid(ax,"on"); box(ax,"on");

for s = 1:numel(solversWanted)
    solver = solversWanted(s);
    Ts = Tc(solverCol(modeCol==modeWanted)==solver,:);
    if isempty(Ts), continue; end

    Ns = double(Ts.Nx).*double(Ts.Ny);
    y  = double(Ts.total_sec_median);
    e  = 0.5*double(Ts.total_sec_iqr);

    [Ns,ord] = sort(Ns);
    y = y(ord); e = e(ord);

    errorbar(ax,Ns,y,e,"o-","LineWidth",1.2,"DisplayName",char(solver));
end

set(ax,"XScale","log","YScale","log");
set(ax,"XTick",Nu,"XTickLabel",gridLabels);

xlabel(ax,"Grid size (N = N_x × N_y)","FontSize",14,"FontWeight","bold");
ylabel(ax,"Median runtime over 500 steps (s)","FontSize",14,"FontWeight","bold");
title(ax,"Compute-only benchmark (compute\_off)","FontSize",15,"FontWeight","bold");

set(ax,"FontSize",12,"FontWeight","bold","LineWidth",1.2);
lgd = legend(ax,"Location","best","Interpreter","none");
set(lgd,"FontSize",12,"FontWeight","bold");

set(fig,"PaperPositionMode","auto");
base = fullfile(outDir,"bench_compute_loglog");
print(fig,base,"-dpng","-r200");
print(fig,base,"-dpdf");
close(fig);
end

function make_render_loglog(T, outDir)

modeWanted    = "render_on";
solversWanted = ["matrix","stencil","full"];

solverCol = string(T.diffusionMode);
modeCol   = string(T.mode);

Tc = T(modeCol == modeWanted,:);
assert(~isempty(Tc), "No render_on data found.");

Nx = double(Tc.Nx);
Ny = double(Tc.Ny);
N  = Nx .* Ny;

Nu = sort(unique(N));
gridLabels = arrayfun(@(n) sprintf("%d^2", round(sqrt(n))), Nu, "UniformOutput", false);

fig = figure("Visible","off","Position",[100 100 900 520]);
ax  = axes("Parent",fig);
hold(ax,"on"); grid(ax,"on"); box(ax,"on");

for s = 1:numel(solversWanted)
    solver = solversWanted(s);
    Ts = Tc(solverCol(modeCol==modeWanted)==solver,:);
    if isempty(Ts), continue; end

    Ns = double(Ts.Nx).*double(Ts.Ny);
    y  = double(Ts.total_sec_median);
    e  = 0.5*double(Ts.total_sec_iqr);

    [Ns,ord] = sort(Ns);
    y = y(ord); e = e(ord);

    errorbar(ax,Ns,y,e,"o-","LineWidth",1.2,"DisplayName",char(solver));
end

set(ax,"XScale","log","YScale","log");
set(ax,"XTick",Nu,"XTickLabel",gridLabels);

xlabel(ax,"Grid size (N = N_x × N_y)","FontSize",14,"FontWeight","bold");
ylabel(ax,"Median runtime over 500 steps (s)","FontSize",14,"FontWeight","bold");
title(ax,"Live-display benchmark (render\_on)","FontSize",15,"FontWeight","bold");

set(ax,"FontSize",12,"FontWeight","bold","LineWidth",1.2);
lgd = legend(ax,"Location","best","Interpreter","none");
set(lgd,"FontSize",12,"FontWeight","bold");

set(fig,"PaperPositionMode","auto");
base = fullfile(outDir,"bench_render_loglog");
print(fig,base,"-dpng","-r200");
print(fig,base,"-dpdf");
close(fig);
end

function make_memory_loglog(T, outDir)

solversWanted = ["matrix","stencil","full"];
solverCol = string(T.diffusionMode);

Nx = double(T.Nx);
Ny = double(T.Ny);
N  = Nx .* Ny;

Nu = sort(unique(N));
gridLabels = arrayfun(@(n) sprintf("%d^2", round(sqrt(n))), Nu, "UniformOutput", false);

fig = figure("Visible","off","Position",[100 100 900 520]);
ax  = axes("Parent",fig);
hold(ax,"on"); grid(ax,"on"); box(ax,"on");

for s = 1:numel(solversWanted)
    solver = solversWanted(s);
    Ts = T(solverCol==solver,:);
    if isempty(Ts), continue; end

    Ns = double(Ts.Nx).*double(Ts.Ny);
    mem = double(Ts.mem_bytes_est)/(1024^2); % MB

    [Ns,ord] = sort(Ns);
    mem = mem(ord);

    plot(ax,Ns,mem,"o-","LineWidth",1.2,"DisplayName",char(solver));
end

set(ax,"XScale","log","YScale","log");
set(ax,"XTick",Nu,"XTickLabel",gridLabels);

xlabel(ax,"Grid size (N = N_x × N_y)","FontSize",14,"FontWeight","bold");
ylabel(ax,"Estimated memory (MB)","FontSize",14,"FontWeight","bold");
title(ax,"Memory scaling of diffusion realizations","FontSize",15,"FontWeight","bold");

set(ax,"FontSize",12,"FontWeight","bold","LineWidth",1.2);
lgd = legend(ax,"Location","northwest","Interpreter","none");
set(lgd,"FontSize",12,"FontWeight","bold");

set(fig,"PaperPositionMode","auto");
base = fullfile(outDir,"bench_memory_loglog");
print(fig,base,"-dpng","-r200");
print(fig,base,"-dpdf");
close(fig);
end
