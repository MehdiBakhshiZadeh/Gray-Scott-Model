function printResult(r)
%PRINTRESULT Print one standardized test result in a compact, robust format.

% --- PASS / FAIL ---
status = "FAIL";
if isfield(r, "pass") && r.pass
    status = "PASS";
end

fprintf("[%s] %s\n", status, r.name);

% --- Header fields (ordered) ---
if isfield(r, "grid") && isnumeric(r.grid) && numel(r.grid) == 2 && all(isfinite(r.grid))
    fprintf("  grid: [%d %d]\n", r.grid(1), r.grid(2));
end

% Prefer r.T, fall back to r.tFinal
if isfield(r, "T") && isfinite(r.T)
    fprintf("  T: %g\n", r.T);
elseif isfield(r, "tFinal") && isfinite(r.tFinal)
    fprintf("  T: %g\n", r.tFinal);
end

if isfield(r, "dt") && isfinite(r.dt)
    fprintf("  dt: %g\n", r.dt);
end

% --- Metrics ---
if isfield(r, "metrics") && isstruct(r.metrics)
    fn = sort(fieldnames(r.metrics));
    for i = 1:numel(fn)
        k = fn{i};
        val = r.metrics.(k);

        fprintf("  %s: ", k);

        if isnumeric(val) && isscalar(val)
            fprintf("%g", val);
        elseif isnumeric(val)
            sz = size(val);
            fprintf("<numeric %dx%d>", sz(1), sz(2));
        elseif islogical(val)
            fprintf("%s", mat2str(val));
        elseif ischar(val) || isstring(val)
            fprintf("%s", string(val));
        else
            fprintf("<%s>", class(val));
        end

        % Threshold (if scalar)
        if isfield(r, "thresholds") && isfield(r.thresholds, k)
            thr = r.thresholds.(k);
            if isnumeric(thr) && isscalar(thr)
                fprintf("  (tol %g)", thr);
            end
        end

        fprintf("\n");
    end
end

% --- Figures (print only filename) ---
if isfield(r, "figFiles") && ~isempty(r.figFiles)
    for i = 1:numel(r.figFiles)
        [~, fname, ext] = fileparts(r.figFiles(i));
        fprintf("  fig: %s%s\n", fname, ext);
    end
end

% --- Notes ---
if isfield(r, "notes") && strlength(r.notes) > 0
    fprintf("  notes: %s\n", r.notes);
end
end
