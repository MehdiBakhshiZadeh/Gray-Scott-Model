# Verification and Test Suite

This folder contains the verification, consistency, validation-oriented sanity, and benchmarking scripts used for the Gray–Scott PPP project.

## Purpose

The tests in this folder are intended to document and reproduce the main numerical checks reported in the verification chapter of the report. They cover manufactured-solution verification, operator accuracy, time-step convergence, explicit-Euler stability-limit behavior, boundary-condition enforcement, solver-mode equivalence, GUI/non-GUI consistency, supplementary structural checks, and a Pearson-type pattern sanity study for the nonlinear Gray–Scott system.

## Recommended execution environment

Run the tests from the project root or from inside the `tests` folder with the project paths available.

A safe MATLAB setup is:

```matlab
projectRoot = '.../gray_scott_ppp';
cd(projectRoot)
addpath(projectRoot)
addpath(genpath(fullfile(projectRoot,'src')))
addpath(genpath(fullfile(projectRoot,'tests')))