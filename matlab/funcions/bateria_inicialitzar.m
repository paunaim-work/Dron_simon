function bat = bateria_inicialitzar(cfg, posBase)
%BATERIA_INICIALITZAR  Estat de la bateria (mode realista). Retorna [] si
%   la limitacio de bateria no es fa servir.
bat = [];
if ~isfield(cfg, 'realista') || ~cfg.realista.activat || ~cfg.realista.bateria.activat
    return;
end
B = cfg.realista.bateria;
bat.E0     = cfg.realista.der.energia0;      % J
bat.E      = bat.E0;
bat.Ph     = B.potencia;                     % W en estacionari
bat.Pm     = B.potencia * (1 + B.extraMoviment);   % W en moviment
bat.retorn = B.retornBase;
bat.v      = cfg.dron.velocitat;
bat.base   = posBase;                        % punt d'enlairament (m)
end
