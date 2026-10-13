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

% Radi de deteccio (esfera): desplacaments de les cel.les que el dron "veu"
if isempty(cfg.dron.radiDeteccio)
    dom.radiDet = 1.5 * dx;
else
    dom.radiDet = cfg.dron.radiDeteccio;
end
mh = floor(dom.radiDet / dx);
mv = floor(dom.radiDet / dz);
[DI, DJ, DK] = ndgrid(-mh:mh, -mh:mh, -mv:mv);
dins = sqrt((DI * dx).^2 + (DJ * dx).^2 + (DK * dz).^2) <= dom.radiDet + 1e-9;
dom.offDet = [DI(dins) DJ(dins) DK(dins)];
end
