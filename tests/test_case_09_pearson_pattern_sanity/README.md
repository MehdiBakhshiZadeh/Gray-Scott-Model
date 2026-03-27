# Test Case 9: Pearson-Type Pattern Sanity Check

## Purpose

This folder contains the files for Test Case 9 of the report, which provides a validation-oriented sanity check of long-time Gray--Scott pattern formation for standard Pearson-type parameter sets.

The test is qualitative. It reproduces several Pearson-style cases and saves representative snapshots for visual inspection.

## Report reference

- Chapter 8
- Section 8.4.1
- Test Case 9: Pearson-Type Pattern Sanity Check
- Related report items:
  - Figures 8.7, 8.8, 8.9, and 8.10

## Main file

- `test_pattern_sanity_pearson.m`

## Required inputs

No external input files are required.

The test generates all required data internally:
- Pearson preset baseline,
- four parameter cases,
- long-time solution history,
- snapshot times.

## Output files

The test writes multiple image files to the local `outputs/` folder, including:

- `alpha_t050000_U.png`
- `alpha_t050000_V.png`
- `...`
- `theta_t200000_U.png`
- `theta_t200000_V.png`
- `pearson_patterns_summary_U_final.png`

The exact set of files depends on the defined snapshot times and cases.

The function also returns a MATLAB result structure listing the generated files.

## How to run

From the project root:

```matlab
r = test_pattern_sanity_pearson();
```

Or from this folder:

```matlab
r = test_pattern_sanity_pearson();
```

## Expected result

This is a qualitative sanity test, so the expected result is:
- successful completion of all long runs,
- stable generated snapshots,
- clear visual differences between the Pearson-type cases,
- a meaningful summary panel.

## Notes

- This test is not included in `run_all_tests.m`.
- It is intentionally separated because it is a long qualitative workflow.
- The test is meant to support the report visually rather than provide a strict convergence proof.
