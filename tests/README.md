# Verification and Test Suite

This folder contains the verification, consistency, validation-oriented sanity, and benchmarking scripts used for the Gray–Scott PPP project.

## Purpose

The tests in this folder are intended to document and reproduce the main verification results reported in Chapter 7 of the report. They cover operator accuracy, time-step behavior, boundary-condition enforcement, solver-mode equivalence, GUI/non-GUI consistency, supplementary structural checks, and a long-time Pearson-type pattern sanity study.

## Recommended execution environment

Run the tests from the project root or from inside the `tests` folder with the project paths available.

A safe MATLAB setup is:

```matlab
projectRoot = '.../gray_scott_ppp';
cd(projectRoot)
addpath(projectRoot)
addpath(genpath(fullfile(projectRoot,'src')))
addpath(genpath(fullfile(projectRoot,'tests')))