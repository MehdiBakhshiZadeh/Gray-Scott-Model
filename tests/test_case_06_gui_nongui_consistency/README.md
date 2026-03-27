# Test Case 6: GUI and Non-GUI Consistency

## Purpose

This folder contains the files for Test Case 6 of the report, which verifies that a simulation run through the GUI produces the same numerical results as the corresponding non-GUI workflow.

This is a multi-step workflow:
1. create a non-GUI reference dump,
2. create a GUI-owned model dump,
3. compare the saved states quantitatively.

## Report reference

- Chapter 8
- Section 8.2.6
- Test Case 6: GUI and Non-GUI Consistency
- Related report items:
  - Table 8.8
  - Figure 8.3

## Main files

- `run_nongui_reference_dump.m`
- `run_gui_consistency_dump.m`
- `compare_gui_nongui.m`

## Required inputs

This workflow can run in two ways:

### Option 1: generate inputs directly
The first two scripts generate the required `.mat` input files automatically.

### Option 2: use local input files
If `compare_gui_nongui()` is called without arguments, it looks for:
- `inputs/nongui_reference_*.mat`
- `inputs/gui_dump_*.mat`

## Output files

Typical generated files are:

- `outputs/nongui_reference_<timestamp>.mat`
- `outputs/gui_dump_<timestamp>.mat`
- `outputs/gui_consistency_table.csv`
- `outputs/gui_consistency.png`

The comparison function also returns a MATLAB table of error metrics and pass/fail results.

## How to run

Recommended workflow:

```matlab
ref = run_nongui_reference_dump();
gui = run_gui_consistency_dump();
T   = compare_gui_nongui(ref.refFile, gui.guiFile);

% Alternative manual workflow:
% 1. Place nongui_reference_*.mat and gui_dump_*.mat in inputs/
% 2. Then run:
T = compare_gui_nongui();
```

Or, if input files are already available in `inputs/`:

```matlab
T = compare_gui_nongui();
```

## Expected result

The GUI and non-GUI workflows should agree exactly or to machine precision when the same setup is used.

The comparison checks:
- snapshot `u`,
- snapshot `v`,
- final `u`,
- final `v`.

The comparison uses strict tolerances:
- relative L2 difference below `1e-12`,
- maximum absolute difference below `1e-12`.

## Notes

- This workflow is not included in `run_all_tests.m`.
- It is intentionally kept separate because it is a multi-step consistency example rather than a short batch test.
- The folder is structured to be self-contained with `inputs/` and `outputs/`.
