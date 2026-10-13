function mostrar_escenari(cfg)
%MOSTRAR_ESCENARI  Escriu a la consola l'escenari sortejat.
if strcmpi(cfg.tipusPloma, 'gaussiana')
    U = cfg.gauss.u;
    extra = sprintf(', classe %s, Q = %.0f', cfg.gauss.estabilitat, cfg.gauss.Q);
else
    U = cfg.complexa.U;
    extra = sprintf(', Q = %.0f', cfg.complexa.Q);
end
fprintf('Escenari aleatori (llavor de l''escenari: %d)\n', cfg.atzar.llavorUsada);
fprintf('   fuga a [%s] m, dron a [%s] m (distancia %.0f m)\n', ...
        strtrim(sprintf('%.0f ', cfg.font.pos)), strtrim(sprintf('%.0f ', cfg.dron.posInicial)), ...
        norm(cfg.dron.posInicial - cfg.font.pos));
fprintf('   vent cap a %.0f graus a %.1f m/s%s\n', cfg.vent.direccio, U, extra);
end
