# Gray-Scott Reaction-Diffusion Project

This project implements and studies the two-dimensional Gray--Scott reaction-diffusion model in MATLAB. It includes a simulation engine, multiple diffusion backends, a graphical user interface built with App Designer, structured verification tests, and a benchmark workflow for solver comparison.

## Project structure

- `run_gui.m`  
  Launches the Gray--Scott graphical user interface.

- `run_model_class.m`  
  Runs one non-GUI simulation using the `GrayScottModel` class.

- `scripts/GrayScottController.mlapp`  
  App Designer GUI for interactive simulation, visualization, saving frames, and video export.

- `src/core/`  
  Core simulation utilities such as parameter handling, grid creation, initial conditions, Euler stepping, and boundary-condition application.

- `src/operators/`  
  Diffusion operators, including sparse-matrix, dense-matrix, stencil-based, and full-array realizations.

- `src/reaction/`  
  Gray--Scott reaction-term implementation.

- `src/model/GrayScottModel.m`  
  Main model class that stores the simulation state and advances the solution in time.

- `src/io/`  
  Plotting and frame export utilities.

- `tests/`  
  Structured verification, consistency, validation-oriented, and benchmark workflows. The `tests` folder now contains:
  - report-linked test-case folders,
  - supplementary verification checks,
  - benchmark workflow files,
  - shared helper utilities,
  - a top-level `README.md` describing the full test layout.

- `report/`  
  LaTeX report sources and report figures.

## Requirements

- MATLAB
- App Designer support for running the GUI (`.mlapp` file)

No additional external libraries are required. Depending on the local MATLAB installation, video export or figure handling may rely on built-in MATLAB components available in the user environment.

## How to run the GUI

From the project root, run:

```matlab
run_gui
```

This opens the Gray--Scott GUI.

## How to run one simulation without the GUI

From the project root, run:

```matlab
run_model_class
```

This creates a results folder under `results/` and stores logs, exported frames, and the final simulation state.

## Diffusion modes

The project supports multiple diffusion implementations through `p.diffusionMode`:

- `"matrix"`  
  sparse-matrix diffusion operator

- `"stencil"`  
  matrix-free stencil diffusion operator

- `"full"`  
  dense/full diffusion realization, mainly useful as a small-scale reference

These modes are used for comparison, benchmarking, and solver verification.

## Tests

The main automated verification runner is:

```matlab
tests/run_all_tests
```

or, if you are already inside the `tests` folder:

```matlab
run_all_tests
```

The core batch suite includes fast verification tests such as:

- constant-field Laplacian check
- Laplacian symmetry check
- manufactured-solution verification
- temporal convergence check
- explicit stability-limit check
- boundary-condition verification
- matrix/stencil operator equivalence
- end-to-end matrix/stencil equivalence

A detailed overview of the complete test structure is given in:

- `tests/README.md`

## Longer workflows outside the main batch runner

Some workflows are intentionally kept outside `run_all_tests.m` because they are longer-running, more specialized, or multi-step examples.

These include:

- GUI/non-GUI consistency workflow  
  `tests/test_case_06_gui_nongui_consistency/`

- Gray--Scott grid-refinement study  
  `tests/test_case_07_grayscott_grid_refinement/`

- Pearson-type pattern sanity check  
  `tests/test_case_09_pearson_pattern_sanity/`

- benchmark workflow  
  `tests/benchmark_solver_performance/`

Each of these folders contains its own local `README.md`.

## Benchmarking

Benchmarking is organized in:

- `tests/benchmark_solver_performance/`

The main benchmark script is:

- `benchmark_solvers.m`

This workflow compares the implemented solver modes with respect to runtime, rendering overhead, and memory usage. Generated benchmark outputs are written to the local `outputs/` folder of the benchmark workflow.

## Test and workflow outputs

Project outputs are written in different places depending on the workflow:

- `results/`  
  general simulation outputs from non-test program runs

- local `outputs/` folders inside the corresponding `tests/` subfolders  
  generated figures, MAT files, CSV files, and workflow-specific artifacts for verification and benchmark scripts

This folder-based organization is intended to make the submission self-contained and easier to inspect.

## Notes

- The project focuses on clarity, modularity, solver comparison, and reproducible verification.
- The GUI and the non-GUI solver share the same numerical core.
- The stencil solver is used as the recommended default solver because it provides a strong overall compromise between runtime, memory efficiency, and scalability.
- The `tests` folder is organized to match the report structure and submission requirements.

## Author

Mehdi Bakhshi Zadeh
