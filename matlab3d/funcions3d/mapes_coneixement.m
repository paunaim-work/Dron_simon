function [Mdalt, Mcostat] = mapes_coneixement(cerc, alg, vis)
%MAPES_CONEIXEMENT  Projeccions (log10) del que sap el dron:
%   Mdalt (Ny x Nx, vista des de dalt) i Mcostat (Nz x Nx, vista x-z).
%   - infotaxis: probabilitat marginal de la posicio de la font
%     (sumant sobre l'eix que no es veu).
%   - gradient: maxim de la concentracio mitjana mesurada a cada cel.la.
switch alg
    case 'infotaxis'
        Mdalt   = log10(max(sum(cerc.P, 3), 1e-300));
        Mcostat = log10(max(reshape(sum(cerc.P, 1), size(cerc.P, 2), size(cerc.P, 3))', 1e-300));
    case 'gradient'
        M = cerc.mapaSuma ./ max(cerc.mapaN, 1);
        M(cerc.mapaN == 0) = 0;
        Mdalt   = log10(max(max(M, [], 3), vis.cmin));
        Mcostat = log10(max(reshape(max(M, [], 1), size(M, 2), size(M, 3))', vis.cmin));
end
end
