# Benchmark and Solver Performance Workflow

## Purpose

This folder contains the benchmark workflow used for the benchmark and solver-performance chapter of the report.

The benchmark compares the implemented diffusion realizations:
- sparse matrix (`matrix`),
- stencil (`stencil`),
- dense matrix (`full`, on restricted small grids only),

and measures performance in two modes:
- compute only,
- compute with live rendering.

## Report reference

- Chapter 10
- Benchmark and Solver Performance Study

This folder belongs to the benchmark chapter.

## Main file

- `benchmark_solvers.m`

## Required inputs

No external input files are required.

The benchmark generates all required data internally:
- grid list,
- solver modes,
- timing modes,
- warmup steps,
- repeat count,
- memory estimates and measurements.

## Output files

The benchmark writes its outputs to a timestamped subfolder inside:

- `outputs/benchmarks_<timestamp>/`

Typical generated files are:

- `outputs/benchmarks_<timestamp>/benchmark_results.mat`
- `outputs/benchmarks_<timestamp>/benchmark_table.csv`

## How to run

From the project root:

```matlab
benchmark_solvers
```

Or from this folder:

```matlab
benchmark_solvers
```

## Notes

- This workflow is not included in `run_all_tests.m`.
- It is intentionally separated from the verification batch because it belongs to the benchmark chapter and may take longer to run.
- The benchmark folder uses its own local `outputs/` directory to remain self-contained.
