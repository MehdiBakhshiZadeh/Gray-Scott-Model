function s = mms_sources(p)
%MMS_SOURCES  Manufactured source terms for MMS verification.
%
%   s = MMS_SOURCES(p) returns forcing terms su(x,y,t) and sv(x,y,t) such that
%   the manufactured solution (u*, v*) from MMS_TARGETS is an exact solution
%   of the forced Gray–Scott system.
%
%   These source terms are intended for verification via the Method of
%   Manufactured Solutions (MMS). They should be supplied to the solver
%   through p.sourceFcn.
%
%   Inputs:
%     p : parameter struct containing Du, Dv, F, k
%
%   Outputs (function handles):
%     s.su(x,y,t) : source term for u-equation
%     s.sv(x,y,t) : source term for v-equation
%     s.uExact    : exact manufactured solution u*
%     s.vExact    : exact manufactured solution v*
%     s.info      : metadata describing the manufactured solution

% --- Input validation ---
required = ["Du","Dv","F","k"];
for i = 1:numel(required)
    assert(isfield(p, required(i)), ...
        "mms_sources: missing required parameter p.%s", required(i));
end

Du = p.Du;
Dv = p.Dv;
F  = p.F;
k  = p.k;

% Manufactured target solution and derivatives
m = mms_targets();

% Reaction terms in the Gray–Scott model:
%   f(u,v) = -u*v^2 + F*(1-u)
%   g(u,v) =  u*v^2 - (F+k)*v
%
% Source terms are defined so that:
%   u_t = Du*lap(u) + f(u,v) + su
%   v_t = Dv*lap(v) + g(u,v) + sv

s.su = @(x,y,t) local_su(x,y,t);
s.sv = @(x,y,t) local_sv(x,y,t);

% Exact fields (for convenience in tests)
s.uExact = m.u;
s.vExact = m.v;
s.info   = m.info;

    function su = local_su(x,y,t)
        u  = m.u(x,y,t);
        v  = m.v(x,y,t);
        uv2 = u .* (v.^2);

        su = m.ut(x,y,t) ...
             - Du .* m.lapU(x,y,t) ...
             + uv2 ...
             - F .* (1 - u);
    end

    function sv = local_sv(x,y,t)
        u  = m.u(x,y,t);
        v  = m.v(x,y,t);
        uv2 = u .* (v.^2);

        sv = m.vt(x,y,t) ...
             - Dv .* m.lapV(x,y,t) ...
             - uv2 ...
             + (F + k) .* v;
    end
end
