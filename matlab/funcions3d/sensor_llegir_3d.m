function [k, c] = sensor_llegir_3d(pl, punts, cfg, dom)
%SENSOR_LLEGIR_3D  Comptes del sensor (o sensors) als punts [x y z] (m)
%   durant un pas de temps: k ~ Poisson( dt*(sensibilitat*c + fons) ).
c  = ploma_concentracio_3d(pl, punts(:, 1), punts(:, 2), punts(:, 3));
mu = dom.dt * (cfg.sensor.sensibilitat * c + cfg.sensor.fons);
k  = poisson_aleatori(mu);
end
