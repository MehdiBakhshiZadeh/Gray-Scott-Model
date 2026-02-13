clear; clc; close all;

%RUN_ALL_TESTS  Execute the full verification test suite.
%
%   This script runs all verification tests defined in the 'tests' folder.
%   Each test is executed independently and returns a standardized result
%   structure indicating pass/fail status and notes.
%
%   Outputs:
%     - A timestamped results_summary.mat file is saved in:
%           results/verification_<timestamp>/
%
%   How to run:
%     From anywhere in MATLAB, simply execute:
%         run_all_tests
%
%   The script does not depend on the current working directory.

% Root of the tests folder
root = fileparts(mfilename("fullpath"));

% Add required paths (absolute, no pwd dependence)
addpath(fullfile(root, "..", "src"));
addpath(fullfile(root, "utils"));

disp("=== Running verification test suite ===");

% Define test list (order matters for readability only)
tests = { ...
    @test_laplacian_constant, ...
    @test_laplacian_symmetry, ...
    @test_timestep_halving_smoke, ...
    @test_mms_endtoend, ...
    @test_mms_operator_convergence, ...
    @test_timestep_convergence, ...
    @test_laplacian_stencil_equivalence, ...
    @test_variantB_endtoend_equivalence ...
};

results = cell(size(tests));
nPass = 0;

% Create output directory with timestamp
timestamp = datestr(now, "yyyy-mm-dd_HH-MM-SS");
outDir = fullfile(root, "..", "results", "verification_" + timestamp);
if ~exist(outDir, "dir")
    mkdir(outDir);
end

for i = 1:numel(tests)
    testName = func2str(tests{i});
    fprintf("\n--- Test %d/%d: %s ---\n", i, numel(tests), testName);

    try
        r = tests{i}();   % each test returns a standardized result struct
    catch ME
        r = makeResult(testName, struct());
        r.pass  = false;
        r.notes = "Crashed: " + string(ME.message);
    end

    results{i} = r;
    printResult(r);

    if isfield(r, "pass") && r.pass
        nPass = nPass + 1;
    end
end

fprintf("\nSummary: %d/%d tests PASS\n", nPass, numel(tests));

% Save results (timestamped, reproducible location)
save(fullfile(outDir, "results_summary.mat"), "results");

disp("Verification results saved in:");
disp(outDir);
