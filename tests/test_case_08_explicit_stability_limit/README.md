# Test Case 8: Explicit Euler Stability Limit

## Purpose

This folder contains the files for Test Case 8 of the report, which verifies the conditional stability limit of the explicit Euler method for the diffusion part of the solver.

The test uses a pure-diffusion subcase and compares behavior:
- below the critical time step,
- at the critical time step,
- above the critical time step.

## Report reference

- Chapter 8
- Section 8.2.8
- Test Case 8: Explicit Euler Stability Limit
- Related report items:
  - Table 8.10
  - Figure 8.5

## Main file

- `test_explicit_stability_limit.m`

## Required inputs

No external input files are required.

The test generates all required data internally:
- periodic grid,
- checkerboard perturbation,
- critical time-step value,
- three test cases,
- perturbation norm histories.

## Output files

The test writes its figure to:

- `outputs/explicit_stability_limit.png`

The function also returns a MATLAB result structure containing:
- theoretical amplification factors,
- observed amplification factors,
- norm histories,
- history errors relative to theory,
- pass/fail result.

## How to run

From the project root:

```matlab
r = test_explicit_stability_limit();
```

Or from this folder:

```matlab
r = test_explicit_stability_limit();
```

## Expected result

The expected behavior is:
- below the critical time step: perturbation decays,
- at the critical time step: perturbation is preserved,
- above the critical time step: perturbation grows.

The test also checks:
- one-step agreement with theory,
- full-history agreement with theory.

## Notes

- This test is included in `run_all_tests.m`.
- The figure is written into the local `outputs/` folder to keep the test self-contained.
