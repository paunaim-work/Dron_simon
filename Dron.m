classdef Dron
%DRON  El dron: els seus parametres, el seu sensor i els algorismes de cerca.
%
%   El dron NO sap on es la fuga: nomes sap el que llegeix el seu sensor.
%   Pot buscar-la amb dos algorismes:
%
%     'infotaxis'  Te un mapa de probabilitat de on pot ser la fuga i va
%                  alla on espera aprendre mes (Vergassola, Villermaux i
%                  Shraiman, Nature 2007).
%     'gradient'   Va cap on la concentracio augmenta. Si perd el senyal,
%                  busca la ploma en zig-zag de traves al vent.
%
%   Aquest fitxer nomes agrupa funcions. Es fan servir aixi:
%
%       par  = Dron.parametres();                                 % valors per defecte
%       dron = Dron.crear('infotaxis', '3d', escenari, par, ploma);
%       dron = Dron.pas(dron, ploma);           % mesura, decideix i es mou
%       Dron.veu_fuga(dron, posFuga)            % true si ja l'ha trobada
%
%   La variable "dron" es una struct amb tot el que cal saber del dron:
%   on es, el seu mapa, que esta fent...
%
%   El dron es mou per una graella de cel.les de mida par.midaCella: a cada
%   pas va a una cel.la veina o es queda quiet. Una cel.la s'identifica pels
%   seus indexs [iy ix iz] (fila, columna i alcada de les matrius).

methods (Static)

    function par = parametres()
    %PARAMETRES  Valors per defecte dels parametres del dron.

        % Vol
        par.velocitat     = 5;      % m/s
        par.midaCella     = 5;      % el que avanca a cada pas (m)
        par.autonomia     = 2500;   % temps maxim de vol (s)
        par.alcadaMinima  = 5;      % alcada minima de vol (m, nomes en 3D)
        par.alcadaMaxima  = 40;     % alcada maxima de vol (m, nomes en 3D)
        par.radiDeteccio  = 5;      % troba la fuga si hi es a menys d'aixo en x i en y (m)
        par.radiDeteccioZ = 2;      % ... i a menys d'aixo en alcada (m)

        % Sensor: comptador de deteccions (com un Geiger). En un pas de dt
        % segons compta  comptes ~ Poisson( dt * (sensibilitat * concentracio + fons) )
        par.sensibilitat = 1;       % comptes/s per unitat de concentracio
        par.fonsSensor   = 0.02;    % comptes/s de fons (falses deteccions)

        % Infotaxis
        par.intervalsComptes = [0 1 2 3 5 8 13 21 34 55 89];   % grups de comptes per planificar

        % Gradient
        par.separacioSensors  = 5;    % distancia dels sensors al centre del dron (m)
        par.suavitzat         = 0.6;  % suavitzat del gradient (0 = gens, 0.9 = molt)
        par.llindarDeteccio   = 2;    % comptes/s (sumant els sensors) per dir que detecta
        par.tempsMemoria      = 4;    % segons sense detectar abans de posar-se a buscar
        par.tramInicial       = 20;   % primer tram del zig-zag (m); cada tram es el doble
        par.canviAlcadaZigZag = 5;    % canvi d'alcada entre trams del zig-zag (m, nomes en 3D)
    end

    function dron = crear(algorisme, dimensio, escenari, par, ploma)
    %CREAR  Dron a punt de sortir de escenari.posInicialDron.
    %   algorisme: 'infotaxis' o 'gradient';  dimensio: '2d' o '3d'.
    %   La ploma nomes es fa servir per construir el model de l'infotaxis
    %   (el dron coneix el tipus de ploma i els seus parametres mitjans).
        dron.algorisme = algorisme;
        dron.par       = par;
        dron.dt        = par.midaCella / par.velocitat;   % durada d'un pas (s)
        dron.maxPassos = ceil(par.autonomia / dron.dt);    % passos que permet l'autonomia

        % Graella de cel.les per on es mou el dron
        g.x = 0:par.midaCella:escenari.midaX;
        g.y = 0:par.midaCella:escenari.midaY;
        if strcmp(dimensio, '3d')
            g.z = par.alcadaMinima:par.midaCella:par.alcadaMaxima;
        else
            g.z = escenari.posInicialDron(3);   % 2D: una sola alcada
        end
        g.mida = [numel(g.y) numel(g.x) numel(g.z)];   % [ny nx nz]
        % Moviments possibles [diy dix diz]: est, nord, oest, sud, amunt, avall, quiet
        g.moviments = [0 1 0; 1 0 0; 0 -1 0; -1 0 0; 0 0 1; 0 0 -1; 0 0 0];
        % Quantes cel.les veu a cada costat [iy ix iz] (si la fuga hi es, la troba)
        g.radiVisible = floor([par.radiDeteccio par.radiDeteccio par.radiDeteccioZ] / par.midaCella + 1e-9);
        dron.graella = g;

        % On es el dron: cel.la mes propera a la posicio inicial
        [~, iz]   = min(abs(g.z - escenari.posInicialDron(3)));
        dron.cel  = [round(escenari.posInicialDron(2) / par.midaCella) + 1, ...
                     round(escenari.posInicialDron(1) / par.midaCella) + 1, iz];
        dron.pos  = Dron.posicio(dron);
        dron.queFa = 'sortint';
        dron.concMesurada = 0;               % ultima concentracio mesurada

        % dron.mapa es el que sap el dron a cada cel.la de la graella
        if strcmp(algorisme, 'infotaxis')
            dron.desplSensors = [0 0 0];     % un sol sensor, al centre del dron
            dron.model = Dron.model_infotaxis(dron, ploma);
            dron.mapa  = ones(g.mida);       % probabilitat de que la fuga sigui a cada cel.la:
            dron.mapa(Dron.zona_visible(dron, dron.cel)) = 0;   % igual a tot arreu menys on ja veu
            dron.mapa  = dron.mapa / sum(dron.mapa(:));
        else
            % Sensors en creu: est, oest, nord, sud i, en 3D, a sobre i a sota
            despl = par.separacioSensors * [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
            dron.desplSensors = despl(1:4 + 2 * strcmp(dimensio, '3d'), :);
            dron.mapa = zeros(g.mida);       % concentracio maxima mesurada a cada cel.la
            dron.grad.filtrat           = [0 0 0];    % gradient suavitzat [dc/dx dc/dy dc/dz]
            dron.grad.passosSenseSenyal = Inf;        % passos seguits sense detectar
            dron.grad.centreZigZag      = dron.cel;   % on busca quan perd el senyal
            dron.grad.tram              = 0;          % tram del zig-zag (0 = no ha comencat)
            dron.grad.desti             = dron.cel;   % cel.la on acaba el tram actual
        end
    end

    function dron = pas(dron, ploma)
    %PAS  Un pas del dron: mesura, decideix on va i s'hi mou.
        p = dron.par;

        % 1. Mesura: comptes de cada sensor
        posSensors = dron.pos + dron.desplSensors;
        conc    = Ploma.concentracio(ploma, posSensors(:, 1), posSensors(:, 2), posSensors(:, 3));
        comptes = Dron.poisson_aleatori(dron.dt * (p.sensibilitat * conc + p.fonsSensor));
        % concentracio que el dron dedueix dels comptes (mitjana dels sensors)
        dron.concMesurada = mean(max(comptes / dron.dt - p.fonsSensor, 0) / p.sensibilitat);

        % 2. Decideix on va
        if strcmp(dron.algorisme, 'infotaxis')
            dron.mapa = Dron.infotaxis_actualitzar(dron, comptes);
            [moviment, dron.queFa] = Dron.infotaxis_decidir(dron);
        else
            c = dron.cel;
            dron.mapa(c(1), c(2), c(3)) = max(dron.mapa(c(1), c(2), c(3)), dron.concMesurada);
            [moviment, dron] = Dron.gradient_decidir(dron, comptes);
        end

        % 3. Es mou a una cel.la veina (o es queda quiet)
        dron.cel = dron.cel + moviment;
        dron.pos = Dron.posicio(dron);
    end

    function trobada = veu_fuga(dron, posFuga)
    %VEU_FUGA  true si la fuga es a tocar del dron: a menys de radiDeteccio
    %   en x i en y, i a menys de radiDeteccioZ en alcada.
        p = dron.par;
        trobada = all(abs(dron.pos - posFuga) <= [p.radiDeteccio p.radiDeteccio p.radiDeteccioZ]);
    end

end

methods (Static, Access = private)

    %% ------------------------------ EINES ------------------------------

    function pos = posicio(dron)
    %POSICIO  Posicio [x y z] (m) de la cel.la on es el dron.
        g = dron.graella;
        pos = [g.x(dron.cel(2)) g.y(dron.cel(1)) g.z(dron.cel(3))];
    end

    function visible = zona_visible(dron, cel)
    %ZONA_VISIBLE  Cel.les que el dron veu des de cel: si la fuga hi fos, la trobaria.
        g = dron.graella;
        r = g.radiVisible;
        visible = false(g.mida);
        visible(max(cel(1) - r(1), 1):min(cel(1) + r(1), g.mida(1)), ...
                max(cel(2) - r(2), 1):min(cel(2) + r(2), g.mida(2)), ...
                max(cel(3) - r(3), 1):min(cel(3) + r(3), g.mida(3))) = true;
    end

    function comptes = poisson_aleatori(mitjana)
    %POISSON_ALEATORI  Nombres aleatoris de Poisson amb la mitjana donada
    %   (sense necessitat de toolboxes).
        comptes = zeros(size(mitjana));
        for i = 1:numel(mitjana)
            if mitjana(i) < 30
                % Metode de Knuth: multiplica numeros a l'atzar fins a baixar de exp(-mitjana)
                producte = rand;
                while producte > exp(-mitjana(i))
                    comptes(i) = comptes(i) + 1;
                    producte = producte * rand;
                end
            else
                comptes(i) = max(0, round(mitjana(i) + sqrt(mitjana(i)) * randn));   % aproximacio normal
            end
        end
    end

    %% ---------------------------- INFOTAXIS ----------------------------

    function model = model_infotaxis(dron, ploma)
    %MODEL_INFOTAXIS  Comptes que el dron espera mesurar en un pas segons on
    %   sigui ell i on sigui la fuga. Nomes depen de la distancia horitzontal
    %   entre tots dos i de les dues alcades, aixi que es guarda en una taula:
    %       model(iyDron - iyFuga + ny, ixDron - ixFuga + nx, izDron, izFuga)
        g  = dron.graella;
        ny = g.mida(1);  nx = g.mida(2);  nz = g.mida(3);
        [dx, dy] = meshgrid((-(nx - 1):(nx - 1)) * dron.par.midaCella, ...
                            (-(ny - 1):(ny - 1)) * dron.par.midaCella);
        model = zeros(2 * ny - 1, 2 * nx - 1, nz, nz);
        for izDron = 1:nz
            for izFuga = 1:nz
                conc = Ploma.mitjana(ploma, dx, dy, g.z(izDron), g.z(izFuga));
                model(:, :, izDron, izFuga) = dron.dt * (dron.par.sensibilitat * conc + dron.par.fonsSensor);
            end
        end
    end

    function esperats = comptes_esperats(dron, cel)
    %COMPTES_ESPERATS  Comptes que el dron espera mesurar a la cel.la cel,
    %   per a cada possible posicio de la fuga (matriu de mida graella.mida).
        ny = dron.graella.mida(1);  nx = dron.graella.mida(2);
        esperats = reshape(dron.model(cel(1) + ny - 1:-1:cel(1), cel(2) + nx - 1:-1:cel(2), cel(3), :), ...
                           dron.graella.mida);
    end

    function S = entropia(prob)
    %ENTROPIA  Entropia de Shannon, S = -sum p ln p: mesura la incertesa d'un
    %   mapa de probabilitat (S = 0 vol dir que se sap del cert on es la fuga).
        prob = prob(prob > 0);
        S = -sum(prob .* log(prob));
    end

    function mapa = infotaxis_actualitzar(dron, comptes)
    %INFOTAXIS_ACTUALITZAR  Regla de Bayes despres de mesurar: cada possible
    %   posicio de la fuga es multiplica per la probabilitat (de Poisson)
    %   d'haver mesurat aquests comptes si la fuga fos alla. La zona que el
    %   dron veu queda descartada (si la fuga hi fos, ja l'hauria trobada).
        mapa     = dron.mapa;
        visible  = Dron.zona_visible(dron, dron.cel);
        esperats = Dron.comptes_esperats(dron, dron.cel);
        nou = mapa .* exp(-esperats + comptes * log(esperats) - gammaln(comptes + 1));
        nou(visible) = 0;
        if sum(nou(:)) > 0
            mapa = nou / sum(nou(:));
        else
            mapa(visible) = 0;               % mesura que el model no explica: nomes descartem la zona
            mapa = mapa / sum(mapa(:));
        end
    end

    function [moviment, queFa] = infotaxis_decidir(dron)
    %INFOTAXIS_DECIDIR  Tria el moviment que fa baixar mes l'entropia (la
    %   incertesa) esperada. Per a cada cel.la veina rj:
    %
    %       canvi(rj) = P(rj) * (-S)  +  (1 - P(rj)) * sum_k rho_k(rj) * dS_k
    %
    %   1r terme: amb probabilitat P(rj) hi troba la fuga (la incertesa passa a 0).
    %   2n terme: si no, mesurara k comptes amb probabilitat rho_k i la
    %   incertesa canviara dS_k. Els comptes s'agrupen en intervals.
        g    = dron.graella;
        prob = dron.mapa;
        S    = Dron.entropia(prob);
        canvi = Inf(7, 1);                   % canvi esperat d'entropia de cada moviment
        for m = 1:7
            celNova = dron.cel + g.moviments(m, :);
            if any(celNova < 1) || any(celNova > g.mida)
                continue;                    % fora de la zona (en 2D, amunt i avall)
            end
            visible    = Dron.zona_visible(dron, celNova);
            probTrobar = sum(prob(visible)); % probabilitat de trobar-hi la fuga
            probSiNo   = prob;               % mapa si NO la troba
            probSiNo(visible) = 0;
            if sum(probSiNo(:)) <= 0
                canvi(m) = -S;
                continue;
            end
            probSiNo = probSiNo / sum(probSiNo(:));
            esperats = Dron.comptes_esperats(dron, celNova);
            possible = probSiNo > 1e-15;     % nomes les posicions encara possibles (mes rapid)
            probSiNo = probSiNo(possible);
            esperats = esperats(possible);

            % Recorrem k = 0, 1, 2... sumant la llei de Poisson interval a interval
            intervals     = dron.par.intervalsComptes;
            canviMesura   = 0;
            probAcumulada = zeros(size(esperats));   % P(comptes < inici de l'interval)
            probK         = exp(-esperats);          % P(k comptes) per al k actual
            k = 0;
            for i = 1:numel(intervals)
                if i < numel(intervals)
                    probInterval = zeros(size(esperats));
                    while k < intervals(i + 1)
                        probInterval = probInterval + probK;
                        probK = probK .* esperats / (k + 1);
                        k = k + 1;
                    end
                    probAcumulada = probAcumulada + probInterval;
                else
                    probInterval = max(1 - probAcumulada, 0);   % ultim interval: la resta
                end
                conjunta   = probSiNo .* probInterval;
                probMesura = sum(conjunta(:));       % probabilitat de mesurar aquest interval
                if probMesura > 1e-300
                    canviMesura = canviMesura + probMesura * (Dron.entropia(conjunta / probMesura) - S);
                end
            end
            canvi(m) = probTrobar * (-S) + (1 - probTrobar) * canviMesura;
        end

        if min(canvi) > -1e-6
            % Ja no pot aprendre res mes: va directe al punt mes probable
            queFa = 'cap al punt mes probable';
            [~, idx] = max(prob(:));
            [iy, ix, iz] = ind2sub(g.mida, idx);
            moviment = Dron.pas_cap_a(dron.cel, [iy ix iz]);
        else
            queFa = 'infotaxis';
            [~, millor] = min(canvi + 1e-10 * abs(S) * rand(7, 1));   % desempat a l'atzar
            moviment = g.moviments(millor, :);
        end
    end

    %% ----------------------------- GRADIENT ----------------------------

    function [moviment, dron] = gradient_decidir(dron, comptes)
    %GRADIENT_DECIDIR  Estima el gradient amb els sensors (est, oest, nord,
    %   sud i, en 3D, a sobre i a sota) per diferencies finites:
    %
    %       dc/dx ~ (cEst - cOest) / 2d,   dc/dy ~ (cNord - cSud) / 2d,   dc/dz ~ (cSobre - cSota) / 2d
    %
    %   i el suavitza. Si te senyal recent, va cap on el gradient es mes gran,
    %   pero mai vent avall (sap d'on ve el vent i la fuga sempre es vent
    %   amunt). Si perd el senyal, busca la ploma en zig-zag de traves al vent
    %   al voltant d'on l'ha perduda (com fan les arnes): cada tram es el
    %   doble de llarg que l'anterior i, en 3D, a una altra alcada.
        p  = dron.par;
        g  = dron.graella;
        gr = dron.grad;
        conc = max(comptes(:)' / dron.dt - p.fonsSensor, 0) / p.sensibilitat;   % concentracio a cada sensor
        gradientNou = [0 0 0];
        gradientNou(1:numel(conc) / 2) = (conc(1:2:end) - conc(2:2:end)) / (2 * p.separacioSensors);
        passosMemoria = round(p.tempsMemoria / dron.dt);   % passos que recorda el senyal

        if sum(comptes) >= p.llindarDeteccio * dron.dt
            gr.passosSenseSenyal = 0;
            gr.filtrat = p.suavitzat * gr.filtrat + (1 - p.suavitzat) * gradientNou;
            gr.centreZigZag = dron.cel;      % si perd el senyal, el buscara al voltant d'aqui
            gr.tram = 0;
        else
            gr.passosSenseSenyal = gr.passosSenseSenyal + 1;
            if gr.passosSenseSenyal > passosMemoria
                gr.filtrat = [0 0 0];        % senyal perdut: oblida el gradient
            end
        end

        if gr.passosSenseSenyal <= passosMemoria
            dron.queFa = 'seguint el gradient';
            veines    = dron.cel + g.moviments(1:6, :);
            dinsZona  = all(veines >= 1, 2) & all(veines <= g.mida, 2);
            projeccio = g.moviments(1:6, [2 1 3]) * gr.filtrat' + 1e-12 * rand(6, 1);   % gradient en cada direccio
            projeccio(~dinsZona) = -Inf;
            projeccio(1) = -Inf;             % mai vent avall
            [~, millor] = max(projeccio);
            moviment = g.moviments(millor, :);
        else
            dron.queFa = 'buscant la ploma en zig-zag';
            if gr.tram == 0 || isequal(dron.cel, gr.desti)
                % Tram nou: cap a l'altre costat, el doble de lluny i a una altra alcada
                gr.tram  = gr.tram + 1;
                amplitud = max(1, round(p.tramInicial / p.midaCella)) * 2^(gr.tram - 1);    % en cel.les
                canviZ   = round(p.canviAlcadaZigZag / p.midaCella) * floor(gr.tram / 2);  % en cel.les
                costat   = (-1)^gr.tram;
                gr.desti = min(max(gr.centreZigZag + costat * [amplitud 0 canviZ], 1), g.mida);
            end
            moviment = Dron.pas_cap_a(dron.cel, gr.desti);
        end
        dron.grad = gr;
    end

    function moviment = pas_cap_a(cel, desti)
    %PAS_CAP_A  Un pas cap a la cel.la desti, per l'eix on falta mes.
        falta = desti - cel;
        [~, eix] = max(abs(falta));
        moviment = [0 0 0];
        moviment(eix) = sign(falta(eix));
    end

end
end
