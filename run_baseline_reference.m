function run_baseline_reference()
% Runs the frozen reference implementation from src/baseline
% (keeps it isolated from src/ to avoid accidental mixing).

root = fileparts(mfilename('fullpath'));

% Add baseline folder first, and remove src from path to avoid shadowing.
addpath(fullfile(root,'src','baseline'));
rmpath(fullfile(root,'src'));  % safe because this is only for reference runs

% Run baseline entry script
run(fullfile(root,'src','baseline','run_grayscott_baseline.m'));
end
