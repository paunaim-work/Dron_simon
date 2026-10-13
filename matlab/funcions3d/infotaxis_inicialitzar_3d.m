function ifx = infotaxis_inicialitzar_3d(cfg, dom, posInicial)
%INFOTAXIS_INICIALITZAR_3D  Estat inicial del cercador infotactic 3D.
%   A priori, la font pot ser a qualsevol punt del volum (x, y i alcada)
%   amb la mateixa probabilitat.

[ifx.K, ifx.model] = model_taxa_deteccio_3d(cfg, dom);
ifx.saturacio = cfg.infotaxis.saturacio;

b = unique(round(cfg.infotaxis.binsComptes(:)'));
b = b(b >= 0 & b < ifx.saturacio);
if isfinite(ifx.saturacio)
    b = unique([b ifx.saturacio]);
end
if isempty(b) || b(1) ~= 0
    b = [0 b];
end
ifx.bins = b;

P = ones(dom.Ny, dom.Nx, dom.Nz);
P(cel_les_properes_3d(posInicial, dom)) = 0;
ifx.P = P / sum(P(:));
ifx.S = entropia(ifx.P);
ifx.dS = zeros(size(dom.moviments, 1), 1);
ifx.mode = 'infotaxis';
end
