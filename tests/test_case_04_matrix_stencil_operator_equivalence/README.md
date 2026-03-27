# Test Case 4: Matrix--Stencil Equivalence of the Diffusion Operator

## Purpose

This folder contains the files for Test Case 4 of the report, which verifies that the sparse-matrix and stencil-based diffusion operators produce the same numerical result when applied to the same field on the same periodic grid.

The test compares the two operator realizations on several grid sizes and several random trials.

## Report reference

- Chapter 8
- Section 8.2.4
- Test Case 4: Matrix--Stencil Equivalence of the Diffusion Operator
- Related report item:
  - Table 8.6

## Main file

- `test_laplacian_stencil_equivalence.m`

## Required inputs

No external input files are required.

The test generates all required data internally:
- test grid sizes,
- random fields,
- sparse matrix operator,
- stencil operator,
- relative and absolute difference measures.

## Output files

This test does not require a mandatory external output file.

The function returns a MATLAB result structure containing:
- all case metrics,
- worst-case relative error,
- worst-case absolute error,
- tolerances,
- pass/fail result.

The local folder `outputs/` is provided for consistency with the project structure and can be used for optional saved summaries in future revisions.

## How to run

From the project root:

```matlab
r = test_laplacian_stencil_equivalence();
```

Or from this folder:

```matlab
r = test_laplacian_stencil_equivalence();
```

## Expected result

The sparse-matrix and stencil-based operators represent the same discrete Laplacian, so the differences should remain at floating-point roundoff level.

The implemented pass criterion requires:
- relative error below `1e-12`,
- maximum absolute difference below `1e-10`,
- successful completion for all tested grids and trials.

## Notes

- This is an operator-level equivalence test.
- Physical boundary conditions are not varied here; the comparison is performed for the periodic diffusion operator itself.
- This test is included in `run_all_tests.m`.
