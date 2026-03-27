# Test Case 1: MMS Operator Convergence of the Discrete Laplacian

## Purpose

This folder contains the files for Test Case 1 of the report, which verifies the spatial accuracy of the discrete Laplacian using the Method of Manufactured Solutions (MMS).

The test checks that the matrix-based finite-difference Laplacian converges with second-order accuracy on a sequence of refined periodic grids when applied to a smooth exact solution.

## Report reference

- Chapter 8
- Section 8.2.1
- Test Case 1: MMS Operator Convergence of the Discrete Laplacian
- Related report items:
  - Table 8.1
  - Figure 8.1

## Main file

- `test_mms_operator_convergence.m`

## Required inputs

No external input files are required.

The test generates all required data internally:
- manufactured exact fields,
- analytic Laplacians,
- grid sequence,
- error data.

## Output files

The test writes its figure to:

- `outputs/mms_operator_convergence.png`

The function also returns a MATLAB result structure containing:
- measured errors,
- observed convergence orders,
- pass/fail result,
- thresholds,
- notes.

## How to run

From the project root:

```matlab
r = test_mms_operator_convergence();
```

Or from this folder:

```matlab
r = test_timestep_convergence();
```

A safe MATLAB setup from the project root is:

```matlab
projectRoot = '.../gray_scott_ppp';
cd(projectRoot)
addpath(projectRoot)
addpath(genpath(fullfile(projectRoot,'src')))
addpath(genpath(fullfile(projectRoot,'tests')))
r = test_mms_operator_convergence();
```

## Expected result

The expected behavior is second-order spatial convergence of the discrete Laplacian.

This means:
- the discrete L2 error should decrease approximately by a factor of four when the grid spacing is halved,
- the measured order should be close to 2.

In the implemented pass criterion, the acceptable order band is:

- `1.7 <= pU <= 2.3`
- `1.7 <= pV <= 2.3`

## Notes

- The test uses periodic boundary conditions.
- The output figure is written to the local `outputs/` folder so that this test case remains self-contained.
- This test is included in `run_all_tests.m`.
