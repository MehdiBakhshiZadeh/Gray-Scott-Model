# GUI / Non-GUI Consistency Test

This folder contains the scripts used to verify that the GUI execution path reproduces the same numerical results as the corresponding non-GUI workflow.

## Purpose

The goal of this test is to confirm that the graphical user interface does not alter the mathematical or numerical behavior of the Gray–Scott solver. The GUI should act only as a control and visualization layer.

The comparison checks:
- midpoint snapshot consistency
- final-state consistency
- agreement in both `u` and `v`

using strict tolerances on relative `L2` and maximum absolute differences.

## Files

### `run_nongui_reference_dump.m`
Runs the reference workflow without the GUI and saves a `.mat` file containing the parameters, seed, settings, and saved states required for comparison.

### `run_gui_consistency_dump.m`
Launches the GUI, uses the current shipped GUI setup, advances the GUI-owned model for the required number of steps, and saves a `.mat` file in the same comparison format as the non-GUI reference.

### `compare_gui_nongui.m`
Loads the non-GUI and GUI dump files, compares snapshot and final states, prints a result table, writes a CSV summary, and saves a difference figure.

## Recommended workflow

Run the following commands with the project paths available:

```matlab
ref = run_nongui_reference_dump();
gui = run_gui_consistency_dump();
T   = compare_gui_nongui(ref.refFile, gui.guiFile);