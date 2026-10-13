function ifx = infotaxis_actualitzar_3d(ifx, pos, k, dom)
%INFOTAXIS_ACTUALITZAR_3D  Actualitzacio bayesiana de P(r0) despres
%   d'observar k comptes a la cel.la pos (eq. 2 de l'article). La zona
%   dins del radi de deteccio passa a tenir probabilitat zero.

MU = mu_local_3d(ifx.K, pos, dom);
L  = versemblanca(k, MU, ifx.saturacio);

idx = cel_les_properes_3d(pos, dom);
Pnova = ifx.P .* L;
Pnova(idx) = 0;
s = sum(Pnova(:));

if s > 0 && isfinite(s)
    ifx.P = Pnova / s;
else
    P = ifx.P;
    P(idx) = 0;
    ifx.P = P / sum(P(:));
end
ifx.S = entropia(ifx.P);
end
