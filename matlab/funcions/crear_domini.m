function dom = crear_domini(cfg)
%CREAR_DOMINI  Construeix la graella de cerca i el pas de temps.
%
%   Les cel.les estan centrades a x = 0, dx, 2dx, ..., Lx (idem per a y).
%   Les matrius sobre la graella tenen mida Ny x Nx i s'indexen (iy, ix).

dx = cfg.domini.dx;
dom.dx = dx;
dom.x  = 0:dx:cfg.domini.Lx;
dom.y  = 0:dx:cfg.domini.Ly;
dom.Nx = numel(dom.x);
dom.Ny = numel(dom.y);
[dom.X, dom.Y] = meshgrid(dom.x, dom.y);

% Temps que tarda el dron a creuar una cel.la (= durada d'un pas)
dom.dt = dx / cfg.dron.velocitat;

% Temps de viatge, d'espera i finestra de mesura de cada pas. En mode
% ideal el dron mesura mentre es mou i cada pas dura dx / velocitat; en
% mode realista pot aturar-se a mesurar (vegeu preparar_realista.m).
dom.tViatge = dom.dt;
dom.tEspera = 0;
dom.tMesura = dom.dt;
if isfield(cfg, 'realista') && cfg.realista.activat
    dom.dt      = cfg.realista.der.tPas;
    dom.tViatge = cfg.realista.der.tViatge;
    dom.tEspera = cfg.realista.der.tEspera;
    dom.tMesura = cfg.realista.der.tMesura;
end

% Direccio del vent mitja (radians) i vector unitari
dom.theta  = cfg.vent.direccio * pi / 180;
dom.uvent  = [cos(dom.theta) sin(dom.theta)];

% Graella fina per dibuixar la ploma
r = cfg.visual.resolucio;
dom.xv = 0:r:cfg.domini.Lx;
dom.yv = 0:r:cfg.domini.Ly;
[dom.Xv, dom.Yv] = meshgrid(dom.xv, dom.yv);

% Moviments possibles (diy, dix): quiet, +x, -x, +y, -y
dom.moviments = [0 0; 0 1; 0 -1; 1 0; -1 0];

% Radi de deteccio: si el dron es a menys d'aquesta distancia de la fuga,
% la troba (p. ex. la veu amb la camera). Desplacaments [diy dix] de les
% cel.les dins del radi:
if isempty(cfg.dron.radiDeteccio)
    dom.radiDet = 1.5 * dx;
else
    dom.radiDet = cfg.dron.radiDeteccio;
end
m = floor(dom.radiDet / dx);
[DI, DJ] = ndgrid(-m:m, -m:m);
dins = hypot(DI, DJ) * dx <= dom.radiDet + 1e-9;
dom.offDet = [DI(dins) DJ(dins)];
end
