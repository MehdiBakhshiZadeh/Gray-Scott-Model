function cases = pearson_cases_alpha_kapa_epsilon_theta()
%PEARSON_CASES_ALPHA_KAPPA_EPSILON_THETA Pearson Fig. 3 parameter points.
%
% This function is the single “source of truth” for the Pearson pattern
% sanity test. All labels are ASCII to avoid font/encoding issues that can
% display Greek letters as '?' on some systems.
%
% Output:
%   cases : struct array with fields
%       .tag   (string)  stable ASCII identifier (lowercase)
%       .label (string)  human-readable label (ASCII)
%       .F     (double)  feed rate
%       .k     (double)  kill rate

cases = repmat(struct("tag","","label","","F",NaN,"k",NaN), 4, 1);

cases(1).tag   = "alpha";
cases(1).label = "Alpha";
cases(1).F     = 0.0075294;
cases(1).k     = 0.04564;

cases(2).tag   = "kappa";
cases(2).label = "Kappa";
cases(2).F     = 0.034487;
cases(2).k     = 0.060378;

cases(3).tag   = "epsilon";
cases(3).label = "Epsilon";
cases(3).F     = 0.019707;
cases(3).k     = 0.055379;

cases(4).tag   = "theta";
cases(4).label = "Theta";
cases(4).F     = 0.040117;
cases(4).k     = 0.060515;

% Safety checks: stop early if any value is missing
for i = 1:numel(cases)
    if ~isfinite(cases(i).F) || ~isfinite(cases(i).k)
        error("Pearson case '%s' is missing F or k.", cases(i).tag);
    end
end
end
