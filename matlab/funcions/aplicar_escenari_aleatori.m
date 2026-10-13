function cfg = aplicar_escenari_aleatori(cfg, mostrar)
%APLICAR_ESCENARI_ALEATORI  Sorteja les condicions inicials de la prova
%   si cfg.atzar.activat es true: posicio de la fuga i del dron, direccio
%   i velocitat del vent, classe d'estabilitat i taxa d'emissio. Cada
%   parametre es pot sortejar o no per separat (cfg.atzar.font,
%   cfg.atzar.dron, cfg.atzar.ventDireccio, ...).
%
%   Les posicions se sortegen dins del domini (a una distancia
%   cfg.atzar.marge de les vores), amb el dron a una distancia minima de
%   la fuga i, si cfg.atzar.dronASotavent, a sotavent: dins d'un con de
%   +-cfg.atzar.angleMaxim graus al voltant de la direccio del vent, com
%   en un cas real en que la contaminacio s'ha detectat vent avall.
%
%   L'escenari te la seva propia llavor (cfg.atzar.llavor), independent de
%   la de la turbulencia i el soroll del sensor (cfg.llavor). Si es buida,
%   cada execucio dona un escenari diferent; la llavor usada es desa a
%   cfg.atzar.llavorUsada i s'escriu a la consola per poder-lo repetir.
%
%   El model intern del dron fa servir els valors sortejats (com si els
%   mesures amb un anemometre i conegues l'emissio). Nomes s'aleatoritzen
%   les dues primeres coordenades de les posicions (x, y).

if nargin < 2
    mostrar = cfg.atzar.mostrar;
end
a = cfg.atzar;
if ~a.activat
    return;
end

llavor = a.llavor;
if isempty(llavor)
    llavor = mod(floor(now * 1e7), 2^31 - 1);
end
cfg.atzar.llavorUsada = llavor;
fixar_llavor(llavor);

L = [cfg.domini.Lx cfg.domini.Ly];
m = a.marge;
font = cfg.font.pos;
dron = cfg.dron.posInicial;
dir  = cfg.vent.direccio;
algunaPosicio = a.ventDireccio || a.font || a.dron;

trobat = false;
for intent = 1:5000
    if a.ventDireccio
        dir = a.rangDireccio(1) + diff(a.rangDireccio) * rand;
    end
    if a.font
        font(1:2) = m + (L - 2 * m) .* rand(1, 2);
    end
    if a.dron
        dron(1:2) = m + (L - 2 * m) .* rand(1, 2);
    end
    v = dron(1:2) - font(1:2);
    dist = norm(v);
    ok = dist >= a.distanciaMinima;
    if ok && a.dronASotavent
        ok = dot(v, [cosd(dir) sind(dir)]) >= dist * cosd(a.angleMaxim);
    end
    if ok || ~algunaPosicio
        trobat = true;
        break;
    end
end
if ~trobat
    warning(['No s''ha trobat cap escenari amb el dron a mes de %g m i a sotavent ' ...
             'de la fuga: es fa servir l''ultim sortejat.'], a.distanciaMinima);
end
cfg.font.pos        = font;
cfg.dron.posInicial = dron;
cfg.vent.direccio   = mod(dir, 360);

if a.ventVelocitat
    U = a.rangVelocitat(1) + diff(a.rangVelocitat) * rand;
    cfg.gauss.u    = U;
    cfg.complexa.U = U;
end
if a.estabilitat
    cfg.gauss.estabilitat = a.classes(randi(numel(a.classes)));
end
cfg.atzar.factorEmissio = 1;
if a.emissio
    f = exp(log(a.rangEmissio(1)) + diff(log(a.rangEmissio)) * rand);
    cfg.gauss.Q    = cfg.gauss.Q * f;
    cfg.complexa.Q = cfg.complexa.Q * f;
    cfg.atzar.factorEmissio = f;
end

if mostrar
    mostrar_escenari(cfg);
end
end
