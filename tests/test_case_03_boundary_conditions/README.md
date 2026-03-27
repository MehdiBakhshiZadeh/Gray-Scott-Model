# Test Case 3: Verification of Boundary-Condition Enforcement

## Purpose

This folder contains the files for Test Case 3 of the report, which verifies that the implemented boundary-condition treatment is applied correctly and robustly in the Gray--Scott solver.

The test checks three aspects:
- Dirichlet constant enforcement over multiple time steps,
- Neumann derivative consistency over multiple time steps,
- mixed boundary-condition robustness in a longer run.

The test is carried out for all three implemented diffusion modes:
- `matrix`
- `stencil`
- `full`

## Report reference

- Chapter 8
- Section 8.2.3
- Test Case 3: Verification of Boundary Condition Enforcement
- Related report items:
  - Table 8.3
  - Table 8.4
  - Table 8.5

## Main file

- `test_boundary_conditions.m`

## Required inputs

No external input files are required.

The test generates all required data internally:
- smooth initial fields for the short checks,
- mixed-boundary long-run setup,
- Dirichlet and Neumann target values,
- boundary-condition error measures.

## Output files

This test does not require a mandatory external output file.

The function returns a MATLAB result structure containing:
- boundary-condition metrics for each solver mode,
- pass/fail result,
- thresholds,
- notes.

The local folder `outputs/` is provided for consistency with the project structure and can be used for optional saved summaries in future revisions.

## How to run

From the project root:

```matlab
r = test_boundary_conditions();
```

Or from this folder:

```matlab
r = test_boundary_conditions();
```

A safe MATLAB setup from the project root is:

```matlab
projectRoot = '.../gray_scott_ppp';
cd(projectRoot)
addpath(projectRoot)
addpath(genpath(fullfile(projectRoot,'src')))
addpath(genpath(fullfile(projectRoot,'tests')))
r = test_boundary_conditions();
```

## Expected result

The expected behavior is:
- Dirichlet boundaries remain fixed at the prescribed values,
- Neumann boundaries satisfy the prescribed derivative values,
- the mixed-boundary long run remains stable,
- no NaN or Inf values appear,
- the interior still evolves nontrivially.

The implemented pass criterion requires:
- all three solver modes to pass all checks,
- boundary errors to remain below the strict tolerance,
- measurable interior evolution in the mixed test.

## Notes

- This test is quantitative and not only visual.
- The test is designed to verify boundary handling across all implemented diffusion realizations.
- The local `outputs/` folder is present for structural consistency even though this test currently does not save a mandatory figure or table.
- This test is included in `run_all_tests.m`.
