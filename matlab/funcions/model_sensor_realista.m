function c = model_sensor_realista(c, cfg)
%MODEL_SENSOR_REALISTA  Aplica al model intern del dron les limitacions
%   del sensor que coneix (esta calibrat): la dilucio mitjana de les
%   helices i la saturacio. Les fluctuacions, el retard del sensor i
%   l'error del GPS no hi son: el dron no les pot preveure.
if ~isfield(cfg, 'realista') || ~cfg.realista.activat
    return;
end
R = cfg.realista;
if R.helices.activat
    c = c * (1 - R.helices.dilucio);
end
if R.saturacio.activat
    c = min(c, R.saturacio.valor);
end
end
