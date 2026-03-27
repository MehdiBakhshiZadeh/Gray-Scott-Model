clear; clc; close all;

% RUN_ALL_TESTS  Execute the verification test suite.
%
%   This script runs the verification tests in the tests folder.
%   Each test is executed independently and returns a result structure
%   indicating pass/fail status and notes.
%
%   Outputs:
%     - results/tests/results_summary.mat
%
%   How to run:
%       run_all_tests
%
%   The script does not depend on the current working directory.

% Locate project folders from this file
thisFile    = mfilename("fullpath");
testsDir    = fileparts(thisFile);
projectRoot = fileparts(testsDir);

% Add required paths
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, "src")));
addpath(genpath(fullfile(projectRoot, "tests")));

% Create output folder safely
outDir = fullfile(projectRoot, "results", "tests");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

disp("=== Running verification test suite ===");

% Define test list
tests = { ...
    @test_laplacian_constant, ...
    @test_laplacian_symmetry, ...
    @test_boundary_conditions, ...
    @test_timestep_halving_smoke, ...
    @test_mms_endtoend, ...
    @test_mms_operator_convergence, ...
    @test_timestep_convergence, ...
    @test_explicit_stability_limit, ...
    @test_laplacian_stencil_equivalence, ...
    @test_matrix_stencil_endtoend_equivalence ...
};

results = cell(size(tests));
nPass = 0;

for i = 1:numel(tests)
    testName = func2str(tests{i});
    fprintf("\n--- Test %d/%d: %s ---\n", i, numel(tests), testName);

    try
        r = tests{i}();
    catch ME
        r = makeResult(testName, struct());
        r.pass  = false;
        whereText = "";
        if ~isempty(ME.stack)
            whereText = " [" + string(ME.stack(1).name) + ...
                ", line " + string(ME.stack(1).line) + "]";
        end
        r.notes = "Crashed: " + string(ME.message) + whereText;
    end

    results{i} = r;
    printResult(r);

    if isfield(r, "pass") && r.pass
        nPass = nPass + 1;
    end
end

fprintf("\nSummary: %d/%d tests PASS\n", nPass, numel(tests));

save(fullfile(outDir, "results_summary.mat"), "results");

disp("Verification results saved in results/tests");