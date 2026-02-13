function r = test_laplacian_stencil_equivalence()
%TEST_LAPLACIAN_STENCIL_EQUIVALENCE
% Compare sparse Laplacian vs stencil apply for periodic BC.

testName = "test_laplacian_stencil_equivalence";
r = makeResult(testName, struct());   % standardized struct
r.pass  = false;                      % default fail until proven otherwise
r.notes = "";

try
    % Multiple sizes catches reshape/index mistakes
    sizes = [ ...
        32  32; ...
        64  48; ...
        128 128 ...
    ];

    nTrials = 3;
    relTol  = 1e-12;
    absTol  = 1e-10;

    for s = 1:size(sizes,1)
        Nx = sizes(s,1);
        Ny = sizes(s,2);

        p = defaultParams();
        p.Nx = Nx;
        p.Ny = Ny;
        p.bc = "periodic";

        grid = buildGrid(p);
        L = buildLaplacian2D(p, grid);
        S = buildStencil2D(p, grid);

        N = Nx * Ny;

        for k = 1:nTrials
            u = randn(N,1);

            a = L*u;
            b = applyLaplacianStencil(u, S);

            diff = a - b;

            relerr = norm(diff) / max(1, norm(a));
            maxabs = max(abs(diff));

            assert(relerr < relTol, ...
                "relerr too large (Nx=%d,Ny=%d): %.3e (tol=%.1e)", ...
                Nx, Ny, relerr, relTol);

            assert(maxabs < absTol, ...
                "maxabs too large (Nx=%d,Ny=%d): %.3e (tol=%.1e)", ...
                Nx, Ny, maxabs, absTol);
        end
    end

    r.pass  = true;
    r.notes = "PASS";

catch ME
    r.pass  = false;
    r.notes = "Crashed: " + string(ME.message);
end
end
