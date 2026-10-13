function [k, pl, sens] = pas_realista(pl, pIni, pFi, offsets, sens, cfg, dom, fAvancar, fConc)
%PAS_REALISTA  Simula un pas del dron amb les limitacions tecniques.
%
%   [k, pl, sens] = pas_realista(pl, pIni, pFi, offsets, sens, cfg, dom, fAvancar, fConc)
%
%   1) El dron viatja de pIni a pFi (posicions en m) durant dom.tViatge.
%   2) Si cal, s'atura dom.tEspera segons per mesurar.
%   Durant tot el pas la ploma evoluciona (fAvancar) i cada sensor (un per
%   fila d'offsets) llegeix la concentracio (fConc) al punt on es, amb:
%     - l'error del GPS (el punt real es una mica diferent),
%     - la dilucio i el soroll del corrent de les helices,
%     - el temps de resposta del PID (filtre de primer ordre),
%     - la saturacio del sensor.
%   Els comptes k s'acumulen durant la finestra de mesura (l'espera o, si
%   no s'espera, el viatge): k ~ Poisson(sensibilitat * integral(c) + fons * tMesura).
%   sens es l'estat intern dels sensors (lectura filtrada) entre passos.

R  = cfg.realista;
nS = size(offsets, 1);
D  = numel(pFi);
if isempty(sens)
    sens = zeros(nS, 1);              % el sensor surt d'aire net
end

eGPS = zeros(1, D);
if R.errorGPS.activat
    eGPS = R.errorGPS.valor * randn(1, D);
end

acc = zeros(nS, 1);
mesuraEnEspera = dom.tEspera > 0;

% --- 1. Viatge ---
nV = max(1, ceil(dom.tViatge / R.dtSub));
h  = dom.tViatge / nV;
for i = 1:nV
    pl = fAvancar(pl, h);
    p = pIni + (pFi - pIni) * (i / nV) + eGPS;
    [sens, lectura] = llegir(sens, pl, p, h, R, offsets, fConc);
    if ~mesuraEnEspera
        acc = acc + lectura * h;
    end
end

% --- 2. Espera aturat mesurant ---
if mesuraEnEspera
    nE = max(1, ceil(dom.tEspera / R.dtSub));
    h  = dom.tEspera / nE;
    for i = 1:nE
        pl = fAvancar(pl, h);
        [sens, lectura] = llegir(sens, pl, pFi + eGPS, h, R, offsets, fConc);
        acc = acc + lectura * h;
    end
end

mu = cfg.sensor.sensibilitat * acc + cfg.sensor.fons * dom.tMesura;
k  = poisson_aleatori(mu);

end

% ---------------------------------------------------------------------
function [s, lectura] = llegir(s, pl, p, h, R, offsets, fConc)
%LLEGIR  Una lectura de tots els sensors durant un interval h.
c = fConc(pl, p + offsets);
c = c(:);
if R.helices.activat
    c = c * (1 - R.helices.dilucio);
    if R.helices.soroll > 0
        sg = sqrt(log(1 + R.helices.soroll^2));     % soroll lognormal de mitjana 1
        c = c .* exp(sg * randn(size(c)) - sg^2 / 2);
    end
end
if R.respostaSensor.activat
    tau = R.respostaSensor.T90 / log(10);          % T90 = tau ln(10)
    s = s + (c - s) * (1 - exp(-h / tau));
else
    s = c;
end
lectura = s;
if R.saturacio.activat
    lectura = min(lectura, R.saturacio.valor);
end
end
