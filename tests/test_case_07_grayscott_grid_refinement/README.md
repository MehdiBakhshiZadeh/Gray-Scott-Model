# Test Case 7: Gray--Scott Grid Refinement Test

## Purpose

This folder contains the files for Test Case 7 of the report, which investigates how the full nonlinear Gray--Scott solution changes under spatial grid refinement.

The test uses a smooth, grid-independent initial condition and compares coarser-grid solutions against an interpolated finest-grid reference solution.

## Report reference

- Chapter 8
- Section 8.2.7
- Test Case 7: Gray--Scott Grid Refinement Test
- Related report items:
  - Table 8.9
  - Figure 8.4

## Main file

- `test_grayscott_grid_refinement.m`

## Required inputs

No external input files are required.

The test generates all required data internally:
- smooth grid-independent initial condition,
- full grid sequence,
- finest-grid reference solution,
- interpolated reference data,
- refinement errors.

## Output files

The test writes its figure to:

- `outputs/grayscott_grid_refinement.png`

The function also returns a MATLAB result structure containing:
- observed refinement orders,
- grid-spacing values,
- error values for `u` and `v`,
- pass/fail result,
- notes.

## How to run

From the project root:

```matlab
r = test_grayscott_grid_refinement();
```

Or from this folder:

```matlab
r = test_grayscott_grid_refinement();
```

## Expected result

The expected behavior is:
- monotonic decrease of the coarse-grid error,
- observed refinement rates consistent with second-order spatial behavior,
- successful completion on all tested grids.

The implemented pass criterion requires:
- monotonic decrease of `Eu`,
- monotonic decrease of `Ev`,
- observed orders above the minimum threshold.

Current threshold:
- `minOrder = 1.5`

## Notes

- This test is intentionally kept outside `run_all_tests.m` because it is longer-running.
- The test uses `dt ~ h^2` to reduce the influence of time-discretization error.
