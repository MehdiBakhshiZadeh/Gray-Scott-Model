function [du_diff, dv_diff] = applyDiffusion(u, v, op, p)
%APPLYDIFFUSION  Compute diffusion terms for Gray–Scott model.
%
% Inputs:
%   u, v : state vectors (N x 1)
%   op   : operator struct with fields:
%          op.mode  = "stencil" | "full" | "matrix"
%          op.L     = Laplacian matrix (sparse or full) or []
%          op.S     = stencil data structure or []
%   p    : parameters (expects p.Du, p.Dv)
%
% Outputs:
%   du_diff, dv_diff : diffusion contributions (N x 1)

mode = string(op.mode);

switch mode
    case "stencil"
        du_diff = p.Du * applyLaplacianStencil(u, op.S);
        dv_diff = p.Dv * applyLaplacianStencil(v, op.S);

    case "full"
        du_diff = p.Du * applyFullLaplacian(u, op.L);
        dv_diff = p.Dv * applyFullLaplacian(v, op.L);

    case "matrix"
        du_diff = p.Du * (op.L * u);
        dv_diff = p.Dv * (op.L * v);

    otherwise
        error("applyDiffusion: unknown op.mode='%s'", mode);
end
end
