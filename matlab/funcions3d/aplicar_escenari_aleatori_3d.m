function cfg = aplicar_escenari_aleatori_3d(cfg)
%APLICAR_ESCENARI_ALEATORI_3D  Escenari aleatori en 3D: sorteja el mateix
%   que la versio 2D (posicions horitzontals, vent, estabilitat, emissio;
%   vegeu aplicar_escenari_aleatori.m) i, a mes, l'alcada de la fuga i
%   l'alcada inicial del dron (cfg.atzar.alcadaFont, cfg.atzar.alcadaDron).
a = cfg.atzar;
if ~a.activat
    return;
end
cfg = aplicar_escenari_aleatori(cfg, false);
Lz = cfg.domini.Lz;
if a.alcadaFont
    r = min(max(a.rangAlcadaFont, 0), Lz);
    cfg.font.pos(3) = r(1) + diff(r) * rand;
end
if a.alcadaDron
    r = min(max(a.rangAlcadaDron, 0), Lz);
    cfg.dron.posInicial(3) = r(1) + diff(r) * rand;
end
if a.mostrar
    mostrar_escenari(cfg);
end
end
