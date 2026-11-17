classdef Dron < handle
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
%   Exemple:
%       par  = Dron.parametres();                         % valors per defecte
%       dron = Dron('infotaxis', '3d', escenari, par);    % constructor
%       dron.pas(camp);                      % mesura, decideix i es mou
%       dron.veuFuga(posFuga)                % true si ja l'ha trobada
%
%   "camp" es la ploma tal com surt del fitxer .nc (vegeu Ploma.llegir):
%   eixos x, y, z, t i la matriu c(ix, iy, iz, it). El dron nomes en llegeix
%   la concentracio als punts on te els sensors; no fa servir res mes de Ploma.m.
%
%   Com a Java: "< handle" fa que els objectes es passin per referencia, i
%   el primer argument dels metodes (obj) fa el paper de "this".
%
%   El dron es mou per una graella de cel.les de mida par.midaCella: a cada
%   pas va a una cel.la veina o es queda quiet. Una cel.la s'identifica pels
%   seus indexs [iy ix iz] (fila, columna i alcada de les matrius).

properties (SetAccess = private)     % es poden llegir des de fora pero no canviar
    algorisme            % 'infotaxis' o 'gradient'
    par                  % parametres (vegeu Dron.parametres)
    dt                   % durada d'un pas (s)
    maxPassos            % passos que permet l'autonomia
    graella              % cel.les per on es mou (struct amb x, y, z, mida...)
    cel                  % cel.la on es [iy ix iz]
    pos                  % posicio on es [x y z] (m)
    temps = 0            % segons que fa que vola
    mapa               % el que sap a cada cel.la: probabilitat de que hi sigui
                         % la fuga (infotaxis) o concentracio maxima mesurada (gradient)
    queFa = 'sortint'    % que esta fent, per mostrar-ho al dibuix
    concMesurada = 0     % ultima concentracio mesurada
end

properties (Access = private)        % nomes els fa servir la propia classe
    desplSensors         % on son els sensors respecte al centre del dron [dx dy dz]
    model                % infotaxis: comptes esperats segons on son el dron i la fuga
    % Nomes per a l'algorisme gradient:
    gradientFiltrat = [0 0 0]  % gradient suavitzat [dc/dx dc/dy dc/dz]
    passosSenseSenyal = Inf    % passos seguits sense detectar
    centreZigZag               % cel.la on busca quan perd el senyal
    tram = 0                   % tram del zig-zag (0 = no ha comencat)
    destiTram                  % cel.la on acaba el tram actual
end

methods

    function obj = Dron(algorisme, dimensio, escenari, par)
    %DRON  Constructor: dron a punt de sortir de escenari.posInicialDron.
        obj.algorisme = algorisme;
        obj.par       = par;
        obj.dt        = par.midaCella / par.velocitat;
        obj.maxPassos = ceil(par.autonomia / obj.dt);

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
        obj.graella = g;

        % Cel.la de sortida: la mes propera a la posicio inicial
        [~, iz] = min(abs(g.z - escenari.posInicialDron(3)));
        obj.cel = [round(escenari.posInicialDron(2) / par.midaCella) + 1, ...
                   round(escenari.posInicialDron(1) / par.midaCella) + 1, iz];
        obj.pos = obj.posicio(obj.cel);

        if strcmp(algorisme, 'infotaxis')
            obj.desplSensors = [0 0 0];      % un sol sensor, al centre del dron
            obj.model = obj.crearModel();
            mapaInicial = ones(g.mida);      % la fuga pot ser a qualsevol lloc...
            mapaInicial(obj.zonaVisible(obj.cel)) = 0;   % ... menys on ja veu
            obj.mapa = mapaInicial / sum(mapaInicial(:));
        else
            % Sensors en creu: est, oest, nord, sud i, en 3D, a sobre i a sota
            despl = par.separacioSensors * [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
            obj.desplSensors = despl(1:4 + 2 * strcmp(dimensio, '3d'), :);
            obj.mapa = zeros(g.mida);
            obj.centreZigZag = obj.cel;
            obj.destiTram = obj.cel;
        end
    end

    function pas(obj, camp)
    %PAS  Un pas del dron: mesura, decideix on va i s'hi mou.
        p = obj.par;

        % 1. Mesura: comptes de cada sensor
        conc    = Dron.llegirCamp(camp, obj.pos + obj.desplSensors, obj.temps);
        comptes = Dron.poissonAleatori(obj.dt * (p.sensibilitat * conc + p.fonsSensor));
        % concentracio que el dron dedueix dels comptes (mitjana dels sensors)
        obj.concMesurada = mean(max(comptes / obj.dt - p.fonsSensor, 0) / p.sensibilitat);

        % 2. Decideix on va
        if strcmp(obj.algorisme, 'infotaxis')
            obj.infotaxisActualitzar(comptes);
            moviment = obj.infotaxisDecidir();
        else
            c = obj.cel;
            obj.mapa(c(1), c(2), c(3)) = max(obj.mapa(c(1), c(2), c(3)), obj.concMesurada);
            moviment = obj.gradientDecidir(comptes);
        end

        % 3. Es mou a una cel.la veina (o es queda quiet)
        obj.cel = obj.cel + moviment;
        obj.pos = obj.posicio(obj.cel);
        obj.temps = obj.temps + obj.dt;
    end

    function trobada = veuFuga(obj, posFuga)
    %VEUFUGA  true si la fuga es a tocar del dron: a menys de radiDeteccio
    %   en x i en y, i a menys de radiDeteccioZ en alcada.
        p = obj.par;
        trobada = all(abs(obj.pos - posFuga) <= [p.radiDeteccio p.radiDeteccio p.radiDeteccioZ]);
    end

end

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
        % Model: com creu el dron que es la ploma de mitjana (vegeu concentracioModel).
        % Ha de tenir model.tipus ('gaussiana' o 'complexa') i els parametres
        % d'aquest tipus (els mateixos noms que a Ploma.parametres). Es posa des del main.
        par.model = [];

        % Gradient
        par.separacioSensors  = 5;    % distancia dels sensors al centre del dron (m)
        par.suavitzat         = 0.6;  % suavitzat del gradient (0 = gens, 0.9 = molt)
        par.llindarDeteccio   = 2;    % comptes/s (sumant els sensors) per dir que detecta
        par.tempsMemoria      = 4;    % segons sense detectar abans de posar-se a buscar
        par.tramInicial       = 20;   % primer tram del zig-zag (m); cada tram es el doble
        par.canviAlcadaZigZag = 5;    % canvi d'alcada entre trams del zig-zag (m, nomes en 3D)
    end

end

methods (Access = private)

    %% ------------------------------ EINES ------------------------------

    function pos = posicio(obj, cel)
    %POSICIO  Posicio [x y z] (m) de la cel.la cel.
        g = obj.graella;
        pos = [g.x(cel(2)) g.y(cel(1)) g.z(cel(3))];
    end

    function visible = zonaVisible(obj, cel)
    %ZONAVISIBLE  Cel.les que el dron veu des de cel: si la fuga hi fos, la trobaria.
        g = obj.graella;
        r = g.radiVisible;
        visible = false(g.mida);
        visible(max(cel(1) - r(1), 1):min(cel(1) + r(1), g.mida(1)), ...
                max(cel(2) - r(2), 1):min(cel(2) + r(2), g.mida(2)), ...
                max(cel(3) - r(3), 1):min(cel(3) + r(3), g.mida(3))) = true;
    end

    %% ---------------------------- INFOTAXIS ----------------------------

    function model = crearModel(obj)
    %CREARMODEL  Comptes que el dron espera mesurar en un pas segons on sigui
    %   ell i on sigui la fuga. Nomes depen de la distancia horitzontal entre
    %   tots dos i de les dues alcades, aixi que es guarda en una taula:
    %       model(iyDron - iyFuga + ny, ixDron - ixFuga + nx, izDron, izFuga)
        g  = obj.graella;
        ny = g.mida(1);  nx = g.mida(2);  nz = g.mida(3);
        [dx, dy] = meshgrid((-(nx - 1):(nx - 1)) * obj.par.midaCella, ...
                            (-(ny - 1):(ny - 1)) * obj.par.midaCella);
        model = zeros(2 * ny - 1, 2 * nx - 1, nz, nz);
        for izDron = 1:nz
            for izFuga = 1:nz
                conc = obj.concentracioModel(dx, dy, g.z(izDron), g.z(izFuga));
                model(:, :, izDron, izFuga) = obj.dt * (obj.par.sensibilitat * conc + obj.par.fonsSensor);
            end
        end
    end

    function c = concentracioModel(obj, dx, dy, zDron, zFuga)
    %CONCENTRACIOMODEL  Concentracio MITJANA que el dron creu que mesurara si
    %   es a (dx, dy) de la fuga, ell a l'alcada zDron i la fuga a zFuga.
    %   Es el "cap" del dron: no llegeix el .nc, nomes fa servir par.model.
    %
    %   - Gaussiana: la formula de la ploma gaussiana (vegeu Ploma.m).
    %   - Complexa: solucio mitjana de l'equacio d'adveccio-difusio
    %     (Vergassola et al., 2007). La turbulencia i el meandre del vent
    %     es resumeixen en dues difusivitats (Dh en horitzontal i Dv en
    %     vertical), i el terra es una fuga "imatge" sota terra:
    %
    %       c = Q / (4 pi sqrt(Dh Dv)) exp(dx U / 2Dh) [e^(-r1/abast)/r1 + e^(-r2/abast)/r2]
    %
    %     r1 i r2 son les distancies a la fuga i a la seva imatge (amb la
    %     z escalada per sqrt(Dh/Dv)), i abast es fins on arriba la ploma.
        m = obj.par.model;
        if strcmp(m.tipus, 'gaussiana')
            c  = zeros(size(dx));
            hi = dx > -m.midaFuga;           % vent amunt de la fuga no hi arriba res
            d  = max(dx(hi), 0);
            sy = sqrt((0.08 * d ./ sqrt(1 + 1e-4   * d)).^2 + m.midaFuga^2);
            sz = sqrt((0.06 * d ./ sqrt(1 + 1.5e-3 * d)).^2 + m.midaFuga^2);
            vertical = exp(-(zDron - zFuga)^2 ./ (2 * sz.^2)) + exp(-(zDron + zFuga)^2 ./ (2 * sz.^2));
            c(hi) = m.emissio ./ (2 * pi * m.vent * sy .* sz) .* exp(-dy(hi).^2 ./ (2 * sy.^2)) .* vertical;
        else
            Dh = m.turbulenciaH + m.creixementH + (m.vent * m.meandreGraus * pi / 180)^2 * m.meandreTemps;
            Dv = m.turbulenciaV + m.creixementV;
            abast = sqrt(Dh * m.tempsVida / (1 + m.vent^2 * m.tempsVida / (4 * Dh)));
            r1 = max(sqrt(dx.^2 + dy.^2 + Dh / Dv * (zDron - zFuga)^2), m.midaPuff);   % mai mes a prop
            r2 = max(sqrt(dx.^2 + dy.^2 + Dh / Dv * (zDron + zFuga)^2), m.midaPuff);   % que un puff
            c = m.emissio / (4 * pi * sqrt(Dh * Dv)) * ...
                (exp(dx * m.vent / (2 * Dh) - r1 / abast) ./ r1 + ...
                 exp(dx * m.vent / (2 * Dh) - r2 / abast) ./ r2);
        end
    end

    function esperats = comptesEsperats(obj, cel)
    %COMPTESESPERATS  Comptes que el dron espera mesurar a la cel.la cel,
    %   per a cada possible posicio de la fuga (matriu de mida graella.mida).
        ny = obj.graella.mida(1);  nx = obj.graella.mida(2);
        esperats = reshape(obj.model(cel(1) + ny - 1:-1:cel(1), cel(2) + nx - 1:-1:cel(2), cel(3), :), ...
                           obj.graella.mida);
    end

    function infotaxisActualitzar(obj, comptes)
    %INFOTAXISACTUALITZAR  Regla de Bayes despres de mesurar: cada possible
    %   posicio de la fuga es multiplica per la probabilitat (de Poisson)
    %   d'haver mesurat aquests comptes si la fuga fos alla. La zona que el
    %   dron veu queda descartada (si la fuga hi fos, ja l'hauria trobada).
        visible  = obj.zonaVisible(obj.cel);
        esperats = obj.comptesEsperats(obj.cel);
        nou = obj.mapa .* exp(-esperats + comptes * log(esperats) - gammaln(comptes + 1));
        nou(visible) = 0;
        if sum(nou(:)) == 0                  % mesura que el model no explica:
            nou = obj.mapa;                  % nomes descartem la zona visible
            nou(visible) = 0;
        end
        obj.mapa = nou / sum(nou(:));
    end

    function moviment = infotaxisDecidir(obj)
    %INFOTAXISDECIDIR  Tria el moviment que fa baixar mes l'entropia (la
    %   incertesa) esperada. Per a cada cel.la veina rj:
    %
    %       canvi(rj) = P(rj) * (-S)  +  (1 - P(rj)) * sum_k rho_k(rj) * dS_k
    %
    %   1r terme: amb probabilitat P(rj) hi troba la fuga (la incertesa passa a 0).
    %   2n terme: si no, mesurara k comptes amb probabilitat rho_k i la
    %   incertesa canviara dS_k. Els comptes s'agrupen en intervals.
        g    = obj.graella;
        prob = obj.mapa;
        S    = Dron.entropia(prob);
        canvi = Inf(7, 1);                   % canvi esperat d'entropia de cada moviment
        for m = 1:7
            celNova = obj.cel + g.moviments(m, :);
            if any(celNova < 1) || any(celNova > g.mida)
                continue;                    % fora de la zona (en 2D, amunt i avall)
            end
            visible    = obj.zonaVisible(celNova);
            probTrobar = sum(prob(visible)); % probabilitat de trobar-hi la fuga
            probSiNo   = prob;               % mapa si NO la troba
            probSiNo(visible) = 0;
            if sum(probSiNo(:)) <= 0
                canvi(m) = -S;
                continue;
            end
            probSiNo = probSiNo / sum(probSiNo(:));
            esperats = obj.comptesEsperats(celNova);
            possible = probSiNo > 1e-15;     % nomes les posicions encara possibles (mes rapid)
            probSiNo = probSiNo(possible);
            esperats = esperats(possible);

            % Recorrem k = 0, 1, 2... sumant la llei de Poisson interval a interval
            intervals     = obj.par.intervalsComptes;
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
            obj.queFa = 'cap al punt mes probable';
            [~, idx] = max(prob(:));
            [iy, ix, iz] = ind2sub(g.mida, idx);
            moviment = Dron.pasCapA(obj.cel, [iy ix iz]);
        else
            obj.queFa = 'infotaxis';
            [~, millor] = min(canvi + 1e-10 * abs(S) * rand(7, 1));   % desempat a l'atzar
            moviment = g.moviments(millor, :);
        end
    end

    %% ----------------------------- GRADIENT ----------------------------

    function moviment = gradientDecidir(obj, comptes)
    %GRADIENTDECIDIR  Estima el gradient amb els sensors (est, oest, nord,
    %   sud i, en 3D, a sobre i a sota) per diferencies finites:
    %
    %       dc/dx ~ (cEst - cOest) / 2d,   dc/dy ~ (cNord - cSud) / 2d,   dc/dz ~ (cSobre - cSota) / 2d
    %
    %   i el suavitza. Si te senyal recent, va cap on el gradient es mes gran,
    %   pero mai vent avall (sap d'on ve el vent i la fuga sempre es vent
    %   amunt). Si perd el senyal, busca la ploma en zig-zag de traves al vent
    %   al voltant d'on l'ha perduda (com fan les arnes): cada tram es el
    %   doble de llarg que l'anterior i, en 3D, a una altra alcada.
        p = obj.par;
        g = obj.graella;
        conc = max(comptes(:)' / obj.dt - p.fonsSensor, 0) / p.sensibilitat;   % concentracio a cada sensor
        gradientNou = [0 0 0];
        gradientNou(1:numel(conc) / 2) = (conc(1:2:end) - conc(2:2:end)) / (2 * p.separacioSensors);
        passosMemoria = round(p.tempsMemoria / obj.dt);   % passos que recorda el senyal

        if sum(comptes) >= p.llindarDeteccio * obj.dt
            obj.passosSenseSenyal = 0;
            obj.gradientFiltrat = p.suavitzat * obj.gradientFiltrat + (1 - p.suavitzat) * gradientNou;
            obj.centreZigZag = obj.cel;      % si perd el senyal, el buscara al voltant d'aqui
            obj.tram = 0;
        else
            obj.passosSenseSenyal = obj.passosSenseSenyal + 1;
            if obj.passosSenseSenyal > passosMemoria
                obj.gradientFiltrat = [0 0 0];   % senyal perdut: oblida el gradient
            end
        end

        if obj.passosSenseSenyal <= passosMemoria
            obj.queFa = 'seguint el gradient';
            veines    = obj.cel + g.moviments(1:6, :);
            dinsZona  = all(veines >= 1, 2) & all(veines <= g.mida, 2);
            projeccio = g.moviments(1:6, [2 1 3]) * obj.gradientFiltrat' + 1e-12 * rand(6, 1);   % gradient en cada direccio
            projeccio(~dinsZona) = -Inf;
            projeccio(1) = -Inf;             % mai vent avall
            [~, millor] = max(projeccio);
            moviment = g.moviments(millor, :);
        else
            obj.queFa = 'buscant la ploma en zig-zag';
            if obj.tram == 0 || isequal(obj.cel, obj.destiTram)
                % Tram nou: cap a l'altre costat, el doble de lluny i a una altra alcada
                obj.tram = obj.tram + 1;
                amplitud = max(1, round(p.tramInicial / p.midaCella)) * 2^(obj.tram - 1);    % en cel.les
                canviZ   = round(p.canviAlcadaZigZag / p.midaCella) * floor(obj.tram / 2);  % en cel.les
                costat   = (-1)^obj.tram;
                obj.destiTram = min(max(obj.centreZigZag + costat * [amplitud 0 canviZ], 1), g.mida);
            end
            moviment = Dron.pasCapA(obj.cel, obj.destiTram);
        end
    end

end

methods (Static, Access = private)

    function c = llegirCamp(camp, pos, t)
    %LLEGIRCAMP  El que llegeixen els sensors: concentracio del camp (la
    %   matriu del .nc) als punts pos (una fila [x y z] per sensor) a
    %   l'instant t. Agafa el punt de la graella i l'instant mes propers:
    %   - si t passa de l'ultim instant del fitxer, fa servir l'ultim (aixi
    %     tambe funciona amb fitxers d'un sol instant);
    %   - fora de la zona del fitxer, la concentracio es 0.
        it = Dron.mesProper(camp.t, t);
        c  = zeros(size(pos, 1), 1);
        for i = 1:size(pos, 1)
            x = pos(i, 1);  y = pos(i, 2);  z = pos(i, 3);
            if x >= camp.x(1) && x <= camp.x(end) && y >= camp.y(1) && y <= camp.y(end)
                c(i) = camp.c(Dron.mesProper(camp.x, x), Dron.mesProper(camp.y, y), ...
                              Dron.mesProper(camp.z, z), it);
            end
        end
    end

    function i = mesProper(v, valor)
    %MESPROPER  Index de l'element del vector v mes proper a valor.
        [~, i] = min(abs(v - valor));
    end

    function comptes = poissonAleatori(mitjana)
    %POISSONALEATORI  Nombres aleatoris de Poisson amb la mitjana donada
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

    function S = entropia(prob)
    %ENTROPIA  Entropia de Shannon, S = -sum p ln p: mesura la incertesa d'un
    %   mapa de probabilitat (S = 0 vol dir que se sap del cert on es la fuga).
        prob = prob(prob > 0);
        S = -sum(prob .* log(prob));
    end

    function moviment = pasCapA(cel, desti)
    %PASCAPA  Un pas cap a la cel.la desti, per l'eix on falta mes.
        falta = desti - cel;
        [~, eix] = max(abs(falta));
        moviment = [0 0 0];
        moviment(eix) = sign(falta(eix));
    end

end
end
