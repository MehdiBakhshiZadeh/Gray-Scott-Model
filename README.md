# Gray-Scott Reaction-Diffusion Project

This project implements and studies the 2D Gray-Scott reaction-diffusion model in MATLAB.  
It includes a simulation engine, multiple diffusion backends, a GUI built with App Designer, validation tests, and benchmarking scripts.

## Project structure

- `run_gui.m`  
  Launches the Gray-Scott graphical user interface.

- `run_model_class.m`  
  Runs one non-GUI simulation using the `GrayScottModel` class.

- `scripts/GrayScottController.mlapp`  
  App Designer GUI for interactive simulation, visualization, saving frames, and video export.

- `src/core/`  
  Core simulation utilities such as parameter handling, grid creation, initial conditions, Euler stepping, and boundary-condition application.

- `src/operators/`  
  Diffusion operators:
  - sparse matrix Laplacian
  - dense matrix Laplacian
  - stencil-based Laplacian
  - full-array Laplacian application

- `src/reaction/`  
  Gray-Scott reaction term.

- `src/model/GrayScottModel.m`  
  Main model class that stores the simulation state and advances the solution in time.

- `src/io/`  
  Plotting and frame export utilities.

- `tests/`  
  Verification, consistency, convergence, and benchmark scripts.

- `report/`  
  LaTeX report sources and report figures.

## Requirements

- MATLAB
- App Designer support for running the GUI (`.mlapp` file)

No additional toolboxes are required unless needed by your local MATLAB setup for video export or figure handling.

## How to run the GUI

From the project root, run:

```matlab
run_gui
```

This opens the Gray-Scott GUI.

## How to run one simulation without the GUI

From the project root, run:

```matlab
run_model_class
```

This creates a results folder under `results/` and stores logs, exported frames, and the final simulation state.

## Diffusion modes

The project supports multiple diffusion implementations through `p.diffusionMode`:

- `"matrix"`   : sparse/dense matrix-based Laplacian
- `"stencil"`  : stencil-based diffusion
- `"full"`     : direct full-array Laplacian application

These modes are used for comparison, benchmarking, and solver verification.

## Tests

The main automated test runner is:

```matlab
tests/run_all_tests
```

or, if you are already inside the `tests` folder:

```matlab
run_all_tests
```

The test suite includes checks for:

- Laplacian symmetry
- constant-field consistency
- stencil/matrix equivalence
- end-to-end matrix vs stencil agreement
- manufactured-solution verification
- timestep convergence and stability checks
- boundary-condition behavior
- GUI vs non-GUI consistency support scripts

## Long validation scripts

Some scripts in `tests/` are intended as long or qualitative validation runs rather than fast unit tests.

In particular:

- `test_pattern_sanity_pearson.m`

is a Pearson-style pattern validation script and may take a long time to run.

## Benchmarking

Benchmark and plotting scripts are available in `tests/`, including:

- `benchmark_solvers.m`
- `plot_benchmark_figures.m`
- `plot_laplacian_sparsity.m`
- `plot_timestep_difference.m`

These are used to compare runtime, memory behavior, and output consistency between different diffusion implementations.

## Output folders

Simulation outputs are typically written to:

- `results/`
- `figures/`
- `tests/results/`

Depending on the script, these may contain logs, images, benchmark tables, or exported videos.

## Notes

- The project focuses on clarity, modularity, and comparison of different diffusion implementations.
- The GUI and the non-GUI solver share the same simulation core.
- The stencil solver is included as a memory-efficient alternative to matrix-based diffusion.

## Author

Mehdi Bakhshi Zadeh