# Verification and Test Suite

This folder contains the verification tests, supplementary checks, benchmark workflow, and shared helper files used in the Gray--Scott PPP project.

Each major report-linked test or example has its own folder with its code, local inputs, local outputs, and a dedicated README file.

## Purpose

The material in this folder documents and reproduces the numerical checks and examples referenced in the report. The contents support three main purposes:

1. formal verification of the numerical implementation,
2. supplementary structural and consistency checks,
3. benchmark and application-oriented example workflows.

The main verification topics covered here are:

- manufactured-solution verification,
- operator-accuracy checks,
- temporal convergence checks,
- explicit-Euler stability-limit verification,
- boundary-condition verification,
- equivalence checks between diffusion implementations,
- GUI/non-GUI consistency checks,
- refinement behavior of the full Gray--Scott system,
- validation-oriented Pearson-type pattern sanity checks,
- solver benchmarking.

## Folder structure

The main folders are:

- `test_case_01_mms_operator/`  
  Test Case 1 from the report: MMS operator convergence of the discrete Laplacian.

- `test_case_02_timestep_convergence/`  
  Test Case 2 from the report: temporal convergence of the explicit Euler scheme.

- `test_case_03_boundary_conditions/`  
  Test Case 3 from the report: verification of boundary-condition enforcement.

- `test_case_04_matrix_stencil_operator_equivalence/`  
  Test Case 4 from the report: operator-level equivalence of matrix and stencil diffusion implementations.

- `test_case_05_matrix_stencil_endtoend/`  
  Test Case 5 from the report: end-to-end equivalence of matrix and stencil solver modes.

- `test_case_06_gui_nongui_consistency/`  
  Test Case 6 from the report: GUI and non-GUI consistency workflow.

- `test_case_07_grayscott_grid_refinement/`  
  Test Case 7 from the report: Gray--Scott grid-refinement study.

- `test_case_08_explicit_stability_limit/`  
  Test Case 8 from the report: explicit Euler stability-limit verification.

- `test_case_09_pearson_pattern_sanity/`  
  Test Case 9 from the report: Pearson-type pattern sanity check.

- `supplementary_checks/`  
  Additional smaller checks used in Section 8.3 of the report.

- `benchmark_solver_performance/`  
  Benchmark workflow used for the benchmark and solver-performance chapter of the report.

- `utils/`  
  Shared helper functions used by multiple tests.

- `run_all_tests.m`  
  Batch runner for the core fast verification suite.

## Recommended execution environment

Run the tests from the project root or from inside the relevant test folder. The code is written so that the individual tests can also be started independently.

A safe MATLAB setup is:

```matlab
projectRoot = '.../gray_scott_ppp';
cd(projectRoot)
addpath(projectRoot)
addpath(genpath(fullfile(projectRoot,'src')))
addpath(genpath(fullfile(projectRoot,'tests')))
```

## Main batch execution

The core fast verification suite can be executed with:

```matlab
run_all_tests
```

This runs the following tests:

1. `test_laplacian_constant`
2. `test_laplacian_symmetry`
3. `test_boundary_conditions`
4. `test_timestep_halving_smoke`
5. `test_mms_endtoend`
6. `test_mms_operator_convergence`
7. `test_timestep_convergence`
8. `test_explicit_stability_limit`
9. `test_laplacian_stencil_equivalence`
10.`test_matrix_stencil_endtoend_equivalence`

The batch runner saves:

- `results/tests/results_summary.mat`

## Tests and examples not included in `run_all_tests.m`

The following workflows are intentionally kept outside the main batch runner because they are longer-running or multi-step examples:

- `test_case_06_gui_nongui_consistency/`
- `test_case_07_grayscott_grid_refinement/`
- `test_case_09_pearson_pattern_sanity/`
- `benchmark_solver_performance/`

These should be executed separately using the instructions in their local README files.

## Report-linked map

The following list links each report item to the corresponding folder or script in this `tests` directory.

### Chapter 8: Testing for Verification

#### Section 8.2: Detailed Verification Test Cases

- Section 8.2.1, Test Case 1: MMS Operator Convergence of the Discrete Laplacian  
  Folder: `test_case_01_mms_operator/`  
  Main script: `test_mms_operator_convergence.m`

- Section 8.2.2, Test Case 2: Temporal Convergence of the Explicit Euler Scheme  
  Folder: `test_case_02_timestep_convergence/`  
  Main script: `test_timestep_convergence.m`

- Section 8.2.3, Test Case 3: Verification of Boundary Condition Enforcement  
  Folder: `test_case_03_boundary_conditions/`  
  Main script: `test_boundary_conditions.m`

- Section 8.2.4, Test Case 4: Matrix--Stencil Equivalence of the Diffusion Operator  
  Folder: `test_case_04_matrix_stencil_operator_equivalence/`  
  Main script: `test_laplacian_stencil_equivalence.m`

- Section 8.2.5, Test Case 5: End-to-End Equivalence of Matrix and Stencil Solver Modes  
  Folder: `test_case_05_matrix_stencil_endtoend/`  
  Main script: `test_matrix_stencil_endtoend_equivalence.m`

- Section 8.2.6, Test Case 6: GUI and Non-GUI Consistency  
  Folder: `test_case_06_gui_nongui_consistency/`  
  Main scripts:  
  - `run_nongui_reference_dump.m`  
  - `run_gui_consistency_dump.m`  
  - `compare_gui_nongui.m`

- Section 8.2.7, Test Case 7: Gray--Scott Grid Refinement Test  
  Folder: `test_case_07_grayscott_grid_refinement/`  
  Main script: `test_grayscott_grid_refinement.m`

- Section 8.2.8, Test Case 8: Explicit Euler Stability Limit  
  Folder: `test_case_08_explicit_stability_limit/`  
  Main script: `test_explicit_stability_limit.m`

#### Section 8.3: Supplementary Verification Checks

Folder: `supplementary_checks/`

Contained scripts:
- `test_laplacian_constant.m`
- `test_laplacian_symmetry.m`
- `test_mms_endtoend.m`
- `test_timestep_halving_smoke.m`

#### Section 8.4.1: Test Case 9: Pearson-Type Pattern Sanity Check

Folder: `test_case_09_pearson_pattern_sanity/`  
Main script: `test_pattern_sanity_pearson.m`

### Chapter 10: Benchmark and Solver Performance Study

Folder: `benchmark_solver_performance/`  
Main script: `benchmark_solvers.m`

## Inputs and outputs

In this organized structure, each major test or example uses its own local folders:

- `inputs/` for required input files, if any,
- `outputs/` for generated files such as figures, MAT files, tables, or CSV files.

Unless explicitly stated otherwise in a local README, most tests generate their required data internally and do not require external datasets.

Important notes:

- The GUI consistency workflow uses generated `.mat` dump files as its comparison inputs.
- The Pearson pattern sanity test and benchmark workflow may produce many output files and are therefore kept separate from the fast batch runner.
- Shared helper functions are stored in `utils/`.

## General conventions used in this folder

The current structure follows these rules:

- Each major report-linked test case has its own folder.
- Each major test folder is intended to contain:
  - code,
  - local `inputs/`,
  - local `outputs/`,
  - a local `README.md`.
- Core fast verification tests are included in `run_all_tests.m`.
- Longer workflows and examples are run separately.
- Tests are written so that they do not depend on the current working directory.
- Output files should be written into the local `outputs/` folder of the corresponding test or workflow.

## Summary

Use this folder as follows:

- run `run_all_tests.m` for the core verification suite,
- run the excluded long workflows separately from their own folders,
- check each local README for the exact report reference, run command, inputs, and outputs.

The numerical core of the project is verified primarily through the Chapter 8 test cases, while benchmark-related material is organized separately under the benchmark folder for Chapter 10.
