function flot = flotacio_cfg(cfg)
%FLOTACIO_CFG  Parametres de flotabilitat del gas (o [] si no s'apliquen).
flot = [];
if isfield(cfg, 'realista') && cfg.realista.activat
    flot = cfg.realista.der.flot;
end
end
