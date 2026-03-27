# Test Case 5: End-to-End Equivalence of Matrix and Stencil Solver Modes

## Purpose

This folder contains the files for Test Case 5 of the report, which verifies that the matrix-based and stencil-based diffusion modes produce the same final numerical solution in a full Gray--Scott simulation.

This is an end-to-end solver comparison, not only an isolated operator comparison.

## Report reference

- Chapter 8
- Section 8.2.5
- Test Case 5: End-to-End Equivalence of Matrix and Stencil Solver Modes
- Related report item:
  - Table 8.7

## Main file

- `test_matrix_stencil_endtoend_equivalence.m`

## Required inputs

No external input files are required.

The test generates all required data internally:
- default Gray--Scott parameters,
- same grid for both runs,
- same seed and initial condition,
- final-state difference measures.

## Output files

This test does not require a mandatory external output file.

The function returns a MATLAB result structure containing:
- relative differences in `u` and `v`,
- maximum absolute differences in `u` and `v`,
- tolerances,
- pass/fail result.

The local folder `outputs/` is provided for consistency with the project structure and can be used for optional saved summaries in future revisions.

## How to run

From the project root:

```matlab
r = test_matrix_stencil_endtoend_equivalence();
```

Or from this folder:

```matlab
r = test_matrix_stencil_endtoend_equivalence();
```

## Expected result

The matrix and stencil modes should produce numerically indistinguishable final states.

The implemented pass criterion requires:
- relative differences below `1e-10`,
- maximum absolute differences below `1e-9`.

## Notes

- The same seed is used for both runs so that the initial condition matches exactly.
- This test is included in `run_all_tests.m`.
