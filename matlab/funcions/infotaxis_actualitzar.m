function ifx = infotaxis_actualitzar(ifx, pos, k, dom)
%INFOTAXIS_ACTUALITZAR  Actualitzacio bayesiana de P(r0) (eq. 2 de
%   l'article) despres d'observar k comptes a la cel.la pos:
%
%       P_{t+1}(r0)  ~  P_t(r0) * P(k | r0)
%
%   A mes, la zona visitada (dins del radi de deteccio) te probabilitat
%   zero: la font no hi es, perque si no el dron l'hauria trobada.

MU = mu_local(ifx.K, pos, dom);
L  = versemblanca(k, MU, ifx.saturacio);

idx = cel_les_properes(pos, dom);
Pnova = ifx.P .* L;
Pnova(idx) = 0;
s = sum(Pnova(:));

if s > 0 && isfinite(s)
    ifx.P = Pnova / s;
else
    % Observacio incompatible amb el model (pot passar si el model es
    % dolent): nomes descartem la zona visitada.
    P = ifx.P;
    P(idx) = 0;
    ifx.P = P / sum(P(:));
end
ifx.S = entropia(ifx.P);
end
