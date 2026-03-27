# Shared Test Utilities

## Purpose

This folder contains helper functions shared by multiple tests.

These files are not standalone report test cases. They support the implementation of the verification tests and workflows in the surrounding `tests` directory.

## Typical contents

Examples of shared helper functionality include:
- result-structure creation,
- printed result formatting,
- convergence-rate calculation,
- norm calculations,
- figure saving with metadata,
- Pearson case definitions,
- MMS helper functions.

## Notes

- Files in this folder are support utilities and should not normally be executed as standalone tests.
- They are used by multiple report-linked test cases and workflows.
