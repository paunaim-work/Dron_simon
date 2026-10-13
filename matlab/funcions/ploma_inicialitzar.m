function pl = ploma_inicialitzar(cfg, dom)
%PLOMA_INICIALITZAR  Crea l'estat de la ploma (la "realitat" que el dron
%   mesura). Tipus disponibles:
%
%   'gaussiana' : ploma gaussiana estacionaria (camp mitja, suau). L'unica
%                 font d'aleatorietat es el soroll de recompte del sensor.
%
%   'complexa'  : model lagrangia de "puffs" (Farrell et al., 2002). La
%                 font emet petits nuvols gaussians que son transportats
%                 per un vent que canvia de direccio (meandreig), te
%                 rafegues i remolins de gran escala, es dispersen per
%                 turbulencia i decauen amb un temps de vida tau. El
%                 resultat es una ploma sinuosa, filamentosa i
%                 intermitent: lluny de la font les deteccions son
%                 esporadiques (limit "diluit" de l'article).

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
pl.p.H           = g.H;
pl.p.z           = g.z;
pl.p.estabilitat = g.estabilitat;
pl.p.theta       = dom.theta;
pl.p.sigma0      = g.sigma0;
% El camp es estacionari: es calcula un sol cop per dibuixar-lo
pl.campVisual = ploma_gaussiana_conc(dom.Xv - pl.font(1), dom.Yv - pl.font(2), pl.p);
end

% ---------------------------------------------------------------------
function pl = inicialitzar_complexa(cfg, dom)
c = cfg.complexa;
pl.tipus  = 'complexa';
pl.font   = cfg.font.pos;
pl.t      = 0;
pl.U      = c.U;
pl.theta0 = dom.theta;

% Emissio
pl.nPerPas   = max(1, round(c.puffsPerSeg * dom.dt));
pl.massaPuff = c.Q * dom.dt / pl.nPerPas;
pl.sigma0    = c.sigma0;
pl.Dpuff     = c.Dpuff;
pl.Dturb     = c.Dturb;
pl.tau       = c.tau;

% Vent variable: meandreig i rafegues (processos d'Ornstein-Uhlenbeck)
pl.sigmaTheta  = c.sigmaTheta * pi / 180;
pl.Ttheta      = c.Ttheta;
pl.sigmaRafega = c.sigmaRafega;
pl.Trafega     = c.Trafega;
pl.dTheta      = pl.sigmaTheta * randn;
pl.rafega      = 0;

% Remolins: funcio de corrent psi = sum a_m sin(k_m . (r - desp) + fase_m)
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
pl.desp   = [0 0];   % els remolins viatgen amb el vent mitja (turbulencia "congelada")

% Puffs (vectors columna)
pl.px  = zeros(0, 1);
pl.py  = zeros(0, 1);
pl.ps2 = zeros(0, 1);   % sigma^2 de cada puff
pl.pm  = zeros(0, 1);   % massa de cada puff

marge = 60;
pl.lim = [dom.x(1) - marge, dom.x(end) + marge, dom.y(1) - marge, dom.y(end) + marge];

% Escalfament: deixem que la ploma es desenvolupi abans que arribi el dron
if isempty(c.tEscalfament)
    tEsc = min(1.3 * hypot(cfg.domini.Lx, cfg.domini.Ly) / c.U, 5 * c.tau);
else
    tEsc = c.tEscalfament;
end
for i = 1:ceil(tEsc / dom.dt)
    pl = ploma_avancar(pl, dom.dt);
end
pl.t = 0;
end
