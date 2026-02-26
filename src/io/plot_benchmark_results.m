function plot_benchmark_results(benchmarkDir)
%PLOT_BENCHMARK_RESULTS  Create performance plots from benchmark_table.csv
%
% Usage:
%   plot_benchmark_results("C:\...\results\benchmarks_YYYY-mm-dd_HH-MM-SS")
%
% Outputs (saved in benchmarkDir/figures):
%   runtime_vs_dof_compute_off.png
%   runtime_vs_dof_render_on.png
%   speedup_vs_dof_compute_off.png
%   speedup_vs_dof_render_on.png

assert(nargin==1 && (isstring(benchmarkDir) || ischar(benchmarkDir)), ...
    "Provide benchmarkDir as a string path to results/benchmarks_<timestamp>.");
benchmarkDir = string(benchmarkDir);

csvPath = fullfile(benchmarkDir, "benchmark_table.csv");
assert(exist(csvPath,"file")==2, "Missing file: %s", csvPath);

outFigDir = fullfile(benchmarkDir, "figures");
if ~exist(outFigDir,"dir"); mkdir(outFigDir); end

T = readtable(csvPath);

% Basic validation
reqCols = ["mode","diffusionMode","Nx","Ny","sec_per_step_median"];
for k = 1:numel(reqCols)
    assert(any(string(T.Properties.VariableNames)==reqCols(k)), ...
        "CSV missing required column: %s", reqCols(k));
end

% Compute DoF
T.DoF = double(T.Nx) .* double(T.Ny);

% Ensure consistent string columns
T.mode   = string(T.mode);
T.diffusionMode = lower(strtrim(T.diffusionMode));

% We will make plots per mode
modes = ["compute_off","render_on"];

for mi = 1:numel(modes)
    mode = modes(mi);

    Ti = T(T.mode == mode, :);
    assert(~isempty(Ti), "No rows found for mode=%s", mode);

    solversHere = unique(Ti.diffusionMode);
    needSolvers = ["matrix","stencil"];
    assert(all(ismember(needSolvers, solversHere)), ...
        "Mode=%s must contain both solvers: %s. Found: %s", ...
        mode, join(needSolvers,", "), join(solversHere,", "));

    % Split by solver
    Tm = Ti(Ti.diffusionMode == "matrix", :);
    Ts = Ti(Ti.diffusionMode == "stencil", :);

    assert(height(Tm)>0 && height(Ts)>0, ...
        "Need both solvers (matrix and stencil) for mode=%s.", mode);

    % Sort by DoF
    Tm = sortrows(Tm, "DoF");
    Ts = sortrows(Ts, "DoF");

    % Check matching DoF sets
    assert(isequal(Tm.DoF, Ts.DoF), ...
        "DoF mismatch between matrix and stencil rows for mode=%s.", mode);

    dof = Tm.DoF;

    % -------------------------
    % Plot 1: runtime vs DoF (sec/step median)
    % -------------------------
    fig1 = figure("Color","w","Visible","off");
    ax1 = axes(fig1);
    loglog(ax1, dof, Tm.sec_per_step_median, "o-", "LineWidth", 1.6, "MarkerSize", 7);
    hold(ax1,"on");
    loglog(ax1, dof, Ts.sec_per_step_median, "o-", "LineWidth", 1.6, "MarkerSize", 7);
    grid(ax1,"on");
    xlabel(ax1, "DoF = N_x N_y");
    ylabel(ax1, "Median time per step (s)");
    title(ax1, "Runtime vs DoF (" + mode + ")", "Interpreter","none");
    legend(ax1, {"matrix","stencil"}, "Location","northwest");
    set(ax1, "FontSize", 12, "LineWidth", 1.0, "Box", "on");

    out1 = fullfile(outFigDir, "runtime_vs_dof_" + mode + ".png");
    print(fig1, out1, "-dpng", "-r300");
    close(fig1);

    % -------------------------
    % Plot 2: speedup vs DoF
    % speedup = t_matrix / t_stencil (per step)
    % -------------------------
    speedup = Tm.sec_per_step_median ./ Ts.sec_per_step_median;

    fig2 = figure("Color","w","Visible","off");
    ax2 = axes(fig2);
    semilogx(ax2, dof, speedup, "o-", "LineWidth", 1.6, "MarkerSize", 7);
    grid(ax2,"on");
    xlabel(ax2, "DoF = N_x N_y");
    ylabel(ax2, "Speedup S = t_{matrix} / t_{stencil}");
    title(ax2, "Speedup vs DoF (" + mode + ")", "Interpreter","none");
    yline(ax2, 1.0, "--"); % reference line
    set(ax2, "FontSize", 12, "LineWidth", 1.0, "Box", "on");

    out2 = fullfile(outFigDir, "speedup_vs_dof_" + mode + ".png");
    print(fig2, out2, "-dpng", "-r300");
    close(fig2);
end

fprintf("Benchmark plots saved in:\n%s\n", outFigDir);
end
