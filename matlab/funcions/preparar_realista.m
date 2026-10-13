function cfg = preparar_realista(cfg)
%PREPARAR_REALISTA  Calcula les magnituds que fa servir el mode realista
%   (si cfg.realista.activat) i les desa a cfg.realista.der:
%
%     tViatge, tEspera, tPas : temps de desplacament, d'espera i total d'un pas (s)
%     tMesura                : finestra de temps en que s'acumula la mesura (s)
%     flot                   : parametres de flotabilitat del gas (o buit)
%     biaixVent              : error de direccio del vent del model (graus)
%     energia0, autonomia    : energia util (J) i autonomia en estacionari (s)
%
%   S'ha de cridar despres de fixar la llavor (l'error de l'anemometre es
%   sorteja). Si R.mostrar, escriu un resum a la consola.

if ~isfield(cfg, 'realista') || ~cfg.realista.activat
    return;
end
R = cfg.realista;
d = struct();

% --- Temps de cada pas ---
if R.velocitat.activat
    cfg.dron.velocitat = R.velocitat.valor;
end
d.tViatge = cfg.domini.dx / cfg.dron.velocitat;
d.tEspera = 0;
if R.tempsMesura.activat
    d.tEspera = R.tempsMesura.valor;
end
d.tPas = d.tViatge + d.tEspera;
if d.tEspera > 0
    d.tMesura = d.tEspera;          % mesura amb el dron aturat
else
    d.tMesura = d.tViatge;          % mesura mentre es mou
end

% --- Gas: flotabilitat ---
d.flot = [];
d.gas  = [];
if R.gas.activat
    g = gas_propietats(R.gas.nom);
    gRed = 9.81 * (g.densitatAire - g.densitat) / g.densitatAire;   % gravetat reduida
    d.flot.w0  = sign(gRed) * sqrt(abs(gRed) * R.gas.midaFont);
    d.flot.tau = R.gas.tempsDilucio;
    d.gas = g;
end

% --- Error de l'anemometre ---
d.biaixVent = 0;
if R.errorVent.activat
    d.biaixVent = R.errorVent.valor * randn;
end

% --- Bateria ---
d.energia0  = R.bateria.fraccioUtil * R.bateria.energia * 3600;
d.autonomia = d.energia0 / R.bateria.potencia;

cfg.realista.der = d;

if R.mostrar
    fprintf('Mode REALISTA:\n');
    fprintf('   velocitat %.1f m/s | %.1f s per pas (%.1f s de viatge + %.1f s aturat mesurant)\n', ...
            cfg.dron.velocitat, d.tPas, d.tViatge, d.tEspera);
    if R.respostaSensor.activat
        fprintf('   sensor PID: T90 = %.1f s', R.respostaSensor.T90);
    else
        fprintf('   sensor PID: resposta instantania');
    end
    if R.saturacio.activat
        fprintf(', satura a %.0f ppm', R.saturacio.valor);
    end
    fprintf('\n');
    if R.helices.activat
        fprintf('   helices: dilucio %.0f %%, soroll %.0f %%\n', 100 * R.helices.dilucio, 100 * R.helices.soroll);
    end
    if R.gas.activat
        if d.flot.w0 > 0
            sentit = 'puja';
        else
            sentit = 'baixa';
        end
        fprintf('   gas: %s, %.2f kg/m^3 (aire %.2f): la ploma %s fins a %.0f m\n', d.gas.nom, ...
                d.gas.densitat, d.gas.densitatAire, sentit, abs(d.flot.w0 * d.flot.tau));
    end
    if R.bateria.activat
        fprintf('   bateria: %.0f Wh x %.0f %% -> autonomia %.0f min en estacionari', ...
                R.bateria.energia, 100 * R.bateria.fraccioUtil, d.autonomia / 60);
        if R.bateria.retornBase
            fprintf(' (amb reserva per tornar a la base)');
        end
        fprintf('\n');
    end
    if R.errorGPS.activat
        fprintf('   error GPS: %.1f m', R.errorGPS.valor);
    end
    if R.errorVent.activat
        fprintf(' | error de direccio del vent del model: %+.1f graus', d.biaixVent);
    end
    fprintf('\n');
end
end
