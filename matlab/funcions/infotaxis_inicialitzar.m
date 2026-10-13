function ifx = infotaxis_inicialitzar(cfg, dom, posInicial)
%INFOTAXIS_INICIALITZAR  Estat inicial del cercador infotactic.
%
%   El dron no sap on es la font: la probabilitat a priori P(r0) es
%   uniforme sobre tot el domini (excepte la zona on comenca, perque
%   si la font hi fos ja l'hauria trobada).

[ifx.K, ifx.model] = model_taxa_deteccio(cfg, dom);
ifx.saturacio = cfg.infotaxis.saturacio;

% Intervals de recompte per a la planificacio: cada valor es el limit
% inferior d'un interval; l'ultim interval es obert (">= ultim valor").
% Si el sensor satura, no te sentit distingir per sobre de la saturacio.
b = unique(round(cfg.infotaxis.binsComptes(:)'));
b = b(b >= 0 & b < ifx.saturacio);
if isfinite(ifx.saturacio)
    b = unique([b ifx.saturacio]);
end
if isempty(b) || b(1) ~= 0
    b = [0 b];
end
ifx.bins = b;

P = ones(dom.Ny, dom.Nx);
P(cel_les_properes(posInicial, dom)) = 0;
ifx.P = P / sum(P(:));
ifx.S = entropia(ifx.P);
ifx.dS = zeros(size(dom.moviments, 1), 1);   % variacio esperada per a cada moviment
ifx.mode = 'infotaxis';
end
