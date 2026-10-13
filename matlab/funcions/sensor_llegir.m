function [k, c] = sensor_llegir(pl, punts, cfg, dom)
%SENSOR_LLEGIR  Lectura del sensor (o sensors) del dron.
%
%   [k, c] = sensor_llegir(pl, punts, cfg, dom) retorna el nombre de
%   comptes k (una fila per punt de mesura, punts = [x y] en m) durant un
%   pas de temps dt, i la concentracio real c en aquests punts.
%
%   Model de sensor de recompte (com un comptador Geiger o un detector de
%   particules): les deteccions arriben a l'atzar amb una taxa mitjana
%   proporcional a la concentracio, mes un soroll de fons:
%
%       k ~ Poisson( mu ),   mu = dt * (sensibilitat * c + fons)

c  = ploma_concentracio(pl, punts(:, 1), punts(:, 2));
mu = dom.dt * (cfg.sensor.sensibilitat * c + cfg.sensor.fons);
k  = poisson_aleatori(mu);
end
