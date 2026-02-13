function r = test_variantB_endtoend_equivalence()
%TEST_VARIANTB_ENDTOEND_EQUIVALENCE
% End-to-end comparison using GrayScottModel class:
% baseline sparse diffusion vs stencil diffusion.

testName = "test_variantB_endtoend_equivalence";
r = makeResult(testName, struct());
r.pass  = false;
r.notes = "";

try
    % --- Fixed run setup ---
    p = defaultParams();
    p.bc = "periodic";

    % Short run to keep it fast
    p.T  = 5;
    p.dt = 0.2;

    % If IC uses RNG, keep it reproducible
    rng(1);

    % --- Baseline (matrix) ---
    pA = p;
    pA.diffusionMode = "matrix";
    modelA = GrayScottModel(pA);
    while modelA.t < pA.T
        info = modelA.step();
        if isfield(info,"hasNaNInf") && info.hasNaNInf
            error("Baseline run blew up (NaN/Inf).");
        end
    end
    [UA, VA] = modelA.getFields2D();

    % --- Stencil ---
    rng(1); % reset RNG so IC matches exactly
    pB = p;
    pB.diffusionMode = "stencil";
    modelB = GrayScottModel(pB);
    while modelB.t < pB.T
        info = modelB.step();
        if isfield(info,"hasNaNInf") && info.hasNaNInf
            error("Stencil run blew up (NaN/Inf).");
        end
    end
    [UB, VB] = modelB.getFields2D();

    % --- Compare ---
    uA = UA(:); vA = VA(:);
    uB = UB(:); vB = VB(:);

    du = uA - uB;
    dv = vA - vB;

    rel_u = norm(du) / max(1, norm(uA));
    rel_v = norm(dv) / max(1, norm(vA));
    max_u = max(abs(du));
    max_v = max(abs(dv));

    relTol = 1e-10;
    absTol = 1e-9;

    assert(rel_u < relTol && rel_v < relTol, ...
        "rel mismatch: rel_u=%.3e, rel_v=%.3e (tol=%.1e)", rel_u, rel_v, relTol);
    assert(max_u < absTol && max_v < absTol, ...
        "max mismatch: max_u=%.3e, max_v=%.3e (tol=%.1e)", max_u, max_v, absTol);

    r.pass = true;
    r.metrics.rel_u = rel_u;
    r.metrics.rel_v = rel_v;
    r.metrics.max_u = max_u;
    r.metrics.max_v = max_v;
    r.thresholds.relTol = relTol;
    r.thresholds.absTol = absTol;
    r.notes = "PASS";

catch ME
    r.pass = false;
    r.notes = "Crashed: " + string(ME.message);
end
end
