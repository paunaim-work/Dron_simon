function pl = ploma_inicialitzar_3d(cfg, dom)
%PLOMA_INICIALITZAR_3D  Crea la ploma 3D (la "realitat" que mesura el dron).
%
%   'gaussiana' : ploma gaussiana estacionaria amb reflexio al terra.
%   'complexa'  : puffs gaussians 3D transportats per un vent amb perfil
%                 vertical, meandreig, rafegues i remolins; es dispersen
%                 per turbulencia (mes en horitzontal que en vertical),
%                 reboten al terra i decauen amb un temps de vida tau.

switch lower(cfg.tipusPloma)
    case 'gaussiana'
        pl = inicialitzar_gaussiana(cfg, dom);
    case 'complexa'
        pl = inicialitzar_complexa(cfg, dom);
    otherwise
        error('Tipus de ploma desconegut: %s (usa ''gaussiana'' o ''complexa'')', cfg.tipusPloma);
end
end

% ---------------------------------------------------------------------
function pl = inicialitzar_gaussiana(cfg, dom)
g = cfg.gauss;
pl.tipus = 'gaussiana';
pl.font  = cfg.font.pos;
pl.t     = 0;
pl.p.Q           = g.Q;
pl.p.u           = g.u;
pl.p.estabilitat = g.estabilitat;
pl.p.sigma0      = g.sigma0;
pl.p.theta       = dom.theta;
% Camp estacionari: es calcula un sol cop per dibuixar-lo en 3D
pl.camp3d = ploma_gaussiana_conc_3d(dom.X3 - pl.font(1), dom.Y3 - pl.font(2), ...
                                    dom.Z3, pl.font(3), pl.p);
end

% ---------------------------------------------------------------------
function pl = inicialitzar_complexa(cfg, dom)
c = cfg.complexa;
pl.tipus  = 'complexa';
pl.font   = cfg.font.pos;
pl.t      = 0;
pl.U      = c.U;
pl.zRef   = c.zRef;
pl.pPerfil = c.exponentPerfil;
pl.theta0 = dom.theta;

pl.nPerPas   = max(1, round(c.puffsPerSeg * dom.dt));
pl.massaPuff = c.Q * dom.dt / pl.nPerPas;
pl.sigma0    = c.sigma0;
pl.DpuffH    = c.DpuffH;
pl.DpuffV    = c.DpuffV;
pl.DturbH    = c.DturbH;
pl.DturbV    = c.DturbV;
pl.tau       = c.tau;

pl.sigmaTheta  = c.sigmaTheta * pi / 180;
pl.Ttheta      = c.Ttheta;
pl.sigmaRafega = c.sigmaRafega;
pl.Trafega     = c.Trafega;
pl.dTheta      = pl.sigmaTheta * randn;
pl.rafega      = 0;

M = c.nRemolins;
lam = c.lambdaMin + (c.lambdaMax - c.lambdaMin) * rand(M, 1);
ang = 2 * pi * rand(M, 1);
kmod = 2 * pi ./ lam;
pl.kx     = kmod .* cos(ang);
pl.ky     = kmod .* sin(ang);
pl.fase   = 2 * pi * rand(M, 1);
pl.sigmaA = c.Uremolins ./ kmod * sqrt(2 / max(M, 1));
pl.a      = pl.sigmaA .* randn(M, 1);
pl.Trem   = c.Tremolins;
pl.desp   = [0 0];

% Puffs (vectors columna): posicio, sigma^2 horitzontal i vertical, massa
pl.px   = zeros(0, 1);
pl.py   = zeros(0, 1);
pl.pz   = zeros(0, 1);
pl.ps2h = zeros(0, 1);
pl.ps2v = zeros(0, 1);
pl.pm   = zeros(0, 1);

marge = 60;
pl.lim = [dom.x(1) - marge, dom.x(end) + marge, dom.y(1) - marge, dom.y(end) + marge, ...
          dom.z(end) + marge];

if isempty(c.tEscalfament)
    tEsc = min(1.3 * hypot(cfg.domini.Lx, cfg.domini.Ly) / c.U, 5 * c.tau);
else
    tEsc = c.tEscalfament;
end
for i = 1:ceil(tEsc / dom.dt)
    pl = ploma_avancar_3d(pl, dom.dt);
end
pl.t = 0;
end
