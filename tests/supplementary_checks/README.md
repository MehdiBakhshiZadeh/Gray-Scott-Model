# Supplementary Verification Checks

## Purpose

This folder contains the smaller supplementary checks referenced in Section 8.3 of the report.

These tests support the main verification strategy by checking additional structural and consistency properties of the implementation.

## Report reference

- Chapter 8
- Section 8.3: Supplementary Verification Checks

Contained report-linked items:
- Section 8.3.1: Constant-Field Laplacian Null Test
- Section 8.3.2: Symmetry Check of the Discrete Laplacian Matrix
- Section 8.3.3: End-to-End MMS Consistency Check
- Section 8.3.4: Time-Step Halving Smoke Test

## Contained files

- `test_laplacian_constant.m`
- `test_laplacian_symmetry.m`
- `test_mms_endtoend.m`
- `test_timestep_halving_smoke.m`

## Required inputs

No external input files are required.

All supplementary checks generate their needed data internally.

## Output files

Most supplementary checks do not require mandatory external files.

One exception is:

- `test_mms_endtoend.m`  
  writes:
  - `outputs/mms_endtoend_errV.png`

The other supplementary checks mainly return MATLAB result structures with metrics and pass/fail information.

## How to run

Each test can be run individually, for example:

```matlab
r = test_laplacian_constant();
r = test_laplacian_symmetry();
r = test_mms_endtoend();
r = test_timestep_halving_smoke();
```

These tests are also included in the batch suite:

```matlab
run_all_tests
```

## Notes

- This folder groups short supplementary checks that support the main Chapter 8 verification discussion.
- The local `outputs/` folder is used for any generated supplementary artifacts.
