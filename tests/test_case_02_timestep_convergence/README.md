# Test Case 2: Temporal Convergence of the Explicit Euler Scheme

## Purpose

This folder contains the files for Test Case 2 of the report, which verifies the temporal behavior of the explicit Euler time-integration scheme used in the Gray--Scott solver.

The test runs the same problem to the same final time with three time-step sizes:
- `dt`
- `dt/2`
- `dt/4`

It then compares the final `v`-field using successive differences and estimates the observed temporal order.

## Report reference

- Chapter 8
- Section 8.2.2
- Test Case 2: Temporal Convergence of the Explicit Euler Scheme
- Related report items:
  - Table 8.2
  - Figure 8.2

## Main file

- `test_timestep_convergence.m`

## Required inputs

No external input files are required.

The test generates all required data internally:
- default Gray--Scott parameters,
- fixed initial condition using a reproducible seed,
- three time-step sizes,
- successive-difference error measures.

## Output files

The test writes its figure to:

- `outputs/timestep_convergence.png`

The function also returns a MATLAB result structure containing:
- L2 and Linf successive differences,
- observed temporal orders,
- pass/fail result,
- thresholds,
- notes.

## How to run

From the project root:

```matlab
r = test_timestep_convergence();
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
r = test_timestep_convergence();
```

## Expected result

The explicit Euler method is first-order accurate in time, so the expected behavior is:
- the finer-step difference is smaller than the coarser-step difference,
- the observed order is close to 1.

In the implemented pass criterion:
- both L2 and Linf successive differences must decrease,
- the L2-based observed order must lie in the band:

`0.7 <= pL2 <= 1.3`

## Notes

- The test uses the matrix diffusion mode.
- Plotting and PNG export from the solver are disabled during the run so that only the test figure is produced.
- The output figure is written to the local `outputs/` folder so that this test case remains self-contained.
- This test is included in `run_all_tests.m`.
