function r = makeResult(name, p)
%MAKERESULT Create a standardized test result struct.
%
%   r = makeResult(name, p) returns a struct with a consistent layout used by
%   all verification tests and by the test runner output formatter.
%
%   Fields:
%     r.name        : string test name
%     r.params      : parameter struct used for the test (may be empty)
%     r.grid        : [Nx Ny] grid size (NaN if not set)
%     r.dt          : time step used (NaN if not set)
%     r.T           : final time used (NaN if not set)
%     r.metrics     : struct of computed metrics (scalars preferred)
%     r.thresholds  : struct of thresholds/bands for metrics (optional)
%     r.pass        : logical pass/fail flag
%     r.figFiles    : string array of output figure paths
%     r.notes       : optional short human-readable note string

if nargin < 2
    p = struct();
end

r = struct();
r.name   = string(name);
r.params = p;               % store full params for reproducibility

r.grid = [NaN NaN];         % [Nx Ny]
r.dt   = NaN;
r.T    = NaN;

r.metrics    = struct();
r.thresholds = struct();

r.pass     = false;
r.figFiles = strings(0,1);
r.notes    = "";
end
