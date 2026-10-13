function th = angle_vent_model(cfg, dom)
%ANGLE_VENT_MODEL  Direccio del vent (rad) que fa servir el model intern
%   del dron: la real mes l'error de l'anemometre en mode realista.
th = dom.theta;
if isfield(cfg, 'realista') && cfg.realista.activat
    th = th + cfg.realista.der.biaixVent * pi / 180;
end
end
