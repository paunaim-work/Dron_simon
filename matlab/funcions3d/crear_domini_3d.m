function dom = crear_domini_3d(cfg)
%CREAR_DOMINI_3D  Graella 3D de cerca. Les matrius tenen mida
%   Ny x Nx x Nz i s'indexen (iy, ix, iz); les cel.les estan centrades a
%   x = 0, dx, ..., Lx; y = 0, dx, ..., Ly; z = 0, dz, ..., Lz.

dx = cfg.domini.dx;
dz = cfg.domini.dz;
dom.dx = dx;
dom.dz = dz;
dom.x  = 0:dx:cfg.domini.Lx;
dom.y  = 0:dx:cfg.domini.Ly;
dom.z  = 0:dz:cfg.domini.Lz;
dom.Nx = numel(dom.x);
dom.Ny = numel(dom.y);
dom.Nz = numel(dom.z);

% Temps per pas (el dron avanca una cel.la, horitzontal o vertical)
dom.dt = dx / cfg.dron.velocitat;

% Temps de viatge, d'espera i finestra de mesura de cada pas (vegeu
% crear_domini.m i preparar_realista.m)
dom.tViatge = dom.dt;
dom.tEspera = 0;
dom.tMesura = dom.dt;
if isfield(cfg, 'realista') && cfg.realista.activat
    dom.dt      = cfg.realista.der.tPas;
    dom.tViatge = cfg.realista.der.tViatge;
    dom.tEspera = cfg.realista.der.tEspera;
    dom.tMesura = cfg.realista.der.tMesura;
end

dom.theta = cfg.vent.direccio * pi / 180;

% Graelles de visualitzacio: tall horitzontal fi i volum 3D
r = cfg.visual.resolucio;
dom.xv = 0:r:cfg.domini.Lx;
dom.yv = 0:r:cfg.domini.Ly;
[dom.Xv, dom.Yv] = meshgrid(dom.xv, dom.yv);
r3 = cfg.visual.resolucio3d;
dom.x3 = 0:r3(1):cfg.domini.Lx;
dom.y3 = 0:r3(1):cfg.domini.Ly;
dom.z3 = 0:r3(2):cfg.domini.Lz;
[dom.X3, dom.Y3, dom.Z3] = meshgrid(dom.x3, dom.y3, dom.z3);

% Moviments possibles [diy dix diz]: quiet, +x, -x, +y, -y, amunt, avall
dom.moviments = [0 0 0; 0 1 0; 0 -1 0; 1 0 0; -1 0 0; 0 0 1; 0 0 -1];

% Zona de deteccio (cilindre): la fuga es troba si el dron hi es a menys
% de radiDet en horitzontal i a menys de tolAlcada en vertical (per
% defecte 0: ha de ser exactament a la mateixa alcada que la fuga).
% Desplacaments [diy dix diz] de les cel.les d'aquesta zona:
if isempty(cfg.dron.radiDeteccio)
    dom.radiDet = 1.5 * dx;
else
    dom.radiDet = cfg.dron.radiDeteccio;
end
dom.tolAlcada = cfg.dron.toleranciaAlcada;
mh = floor(dom.radiDet / dx);
mv = floor(dom.tolAlcada / dz + 1e-9);
[DI, DJ, DK] = ndgrid(-mh:mh, -mh:mh, -mv:mv);
dins = hypot(DI * dx, DJ * dx) <= dom.radiDet + 1e-9 & abs(DK * dz) <= dom.tolAlcada + 1e-9;
dom.offDet = [DI(dins) DJ(dins) DK(dins)];
end
