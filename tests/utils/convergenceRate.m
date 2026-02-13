function [p_est, c0] = convergenceRate(h, E)
%CONVERGENCERATE Estimate convergence order from error vs grid spacing.
%
%   [p_est, c0] = convergenceRate(h, E) fits the model
%       E(h) ? C * h^p
%   using a least-squares fit in log-log space:
%       log(E) ? c0 + p * log(h)
%
%   Inputs:
%     h : vector of grid spacings (positive)
%     E : vector of errors (positive), same length as h
%
%   Outputs:
%     p_est : estimated order p
%     c0    : fitted intercept (log(C))

h = h(:);
E = E(:);

if numel(h) ~= numel(E)
    error("convergenceRate: h and E must have the same length.");
end

% Keep only valid points
mask = isfinite(h) & isfinite(E) & (h > 0) & (E > 0);
h = h(mask);
E = E(mask);

if numel(h) < 2
    error("convergenceRate: need at least two positive finite (h,E) pairs.");
end

% Sort by h (small -> large) for consistent behavior and easier debugging
[h, idx] = sort(h, "ascend");
E = E(idx);

% Fit log(E) = c0 + p*log(h)
c = polyfit(log(h), log(E), 1);
p_est = c(1);
c0    = c(2);
end
