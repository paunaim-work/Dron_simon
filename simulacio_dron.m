function simulacio_dron()
%SIMULACIO_DRON  Dron autonom que busca l'origen d'una fuga (versio simplificada).
%
%   Simula una fuga de contaminant i un dron que en busca l'origen nomes
%   amb les lectures del seu sensor. Es pot triar:
%
%     - la dimensio:  '2d' (el dron vola sempre a l'alcada de la fuga i
%                          nomes busca en horitzontal)
%                     '3d' (el dron no sap l'alcada de la fuga i tambe es
%                          mou amunt i avall)
%     - la ploma:     'gaussiana' (estacionaria, model de Pasquill-Gifford)
%                     'complexa'  (turbulenta: "puffs" que el vent s'emporta)
%     - l'algorisme:  'infotaxis' (Vergassola, Villermaux i Shraiman, Nature 2007)
%                     'gradient'  (puja pel gradient de concentracio; si perd
%                                  el senyal, el busca en zig-zag de traves al vent)
%
%   Executa-la amb F5. Les opcions i els parametres son just a sota.

%% ============================= OPCIONS =============================
dimensio           = '3d';          % '2d' o '3d'
tipusPloma         = 'complexa';    % 'gaussiana' o 'complexa'
algorisme          = 'infotaxis';   % 'infotaxis' o 'gradient'
llavor             = 1;             % mateix numero = mateixa simulacio ([] = cada cop diferent)
mostrarAnimacio    = true;          % false = sense animacio (mes rapid)
dibuixaCadaNPassos = 2;             % redibuixa la figura cada N passos

%% ============================ PARAMETRES ============================
% Unitats: metres i segons. El vent bufa cap a +x i el terra es z = 0.
par.midaX = 400;                    % llargada de la zona de cerca (direccio del vent)
par.midaY = 300;                    % amplada de la zona de cerca (de traves al vent)
par.alcadaMaxima = 40;              % alcada maxima de vol (nomes en 3D; la minima es una cel.la)
par.midaCella = 5;                  % mida de cel.la = el que avanca el dron a cada pas
par.posFuga        = [60 150 10];   % on es la fuga [x y z] (el dron NO ho sap)
par.posInicialDron = [340 215 30];  % d'on surt el dron [x y z] (en 2D, z = alcada de la fuga)
par.velocitatDron  = 5;             % m/s
par.maxPassos      = 2500;          % limit de passos (autonomia)

% Sensor: comptador de deteccions (tipus Geiger).
% Comptes en un pas ~ Poisson( tempsPas * (sensibilitat * concentracio + fons) )
par.sensibilitat = 1;               % comptes/s per unitat de concentracio
par.fonsSensor   = 0.02;            % comptes/s de fons (falses deteccions)

% Ploma gaussiana
par.gaussEmissio  = 5000;           % taxa d'emissio
par.gaussVent     = 3;              % velocitat del vent (m/s)
par.gaussMidaFuga = 1;              % mida de la fuga (m)

% Ploma complexa
par.emissio       = 2000;           % taxa d'emissio
par.vent          = 2;              % velocitat mitjana del vent (m/s)
par.tempsVida     = 200;            % temps de vida del contaminant (s)
par.puffsPerSegon = 2;              % puffs emesos per segon
par.midaPuff      = 1.5;            % mida inicial d'un puff (m)
par.creixementH   = 0.4;            % rapidesa amb que s'eixampla un puff en horitzontal (m^2/s)
par.creixementV   = 0.2;            % ... i en vertical (m^2/s)
par.turbulenciaH  = 0.4;            % moviment aleatori dels puffs en horitzontal (m^2/s)
par.turbulenciaV  = 0.2;            % ... i en vertical (m^2/s)
par.meandreGraus  = 15;             % oscil.lacio de la direccio del vent (graus)
par.meandreTemps  = 40;             % temps tipic d'aquesta oscil.lacio (s)

% Infotaxis
par.intervalsComptes = [0 1 2 3 5 8 13 21 34 55 89];   % intervals de comptes per planificar

% Gradient
par.separacioSensors = 5;           % distancia dels sensors al centre del dron (m)
par.suavitzat        = 0.6;         % suavitzat del gradient (0 = cap, 0.9 = molt)
par.llindarDeteccio  = 2;           % comptes (suma dels sensors) per considerar que detecta
par.passosMemoria    = 4;           % passos sense detectar abans de posar-se a buscar
par.tramInicial      = 4;           % primer tram del zig-zag de cerca (cel.les); cada tram es el doble
%% ====================================================================

if ~isempty(llavor)
    rng(llavor);
end
usaInfotaxis  = strcmp(algorisme, 'infotaxis');
plomaComplexa = strcmp(tipusPloma, 'complexa');
tresD         = strcmp(dimensio, '3d');

% --- Graella: les matrius tenen mida ny x nx x nz i s'indexen (iy, ix, iz) ---
graella.x = 0:par.midaCella:par.midaX;
graella.y = 0:par.midaCella:par.midaY;
if tresD
    graella.z = par.midaCella:par.midaCella:par.alcadaMaxima;   % el dron no vola arran de terra
else
    graella.z = par.posFuga(3);                 % 2D: una sola alcada, la de la fuga
    par.posInicialDron(3) = par.posFuga(3);
end
graella.nx = numel(graella.x);
graella.ny = numel(graella.y);
graella.nz = numel(graella.z);
[graella.X, graella.Y] = meshgrid(graella.x, graella.y);
% Moviments possibles [diy dix diz]: est, nord, oest, sud, amunt, avall, quiet
graella.moviments = [0 1 0; 1 0 0; 0 -1 0; -1 0 0; 0 0 1; 0 0 -1; 0 0 0];
tempsPas = par.midaCella / par.velocitatDron;   % durada d'un pas (s)

celFuga = posicio_a_cella(par.posFuga, par, graella);          % cel.la [iy ix iz] de la fuga
celDron = posicio_a_cella(par.posInicialDron, par, graella);   % cel.la [iy ix iz] del dron

% --- Ploma: la complexa evoluciona 300 s abans que surti el dron ---
buit = zeros(0, 1);
ploma = struct('tipus', tipusPloma, 'x', buit, 'y', buit, 'z', buit, ...
               'sigma2H', buit, 'sigma2V', buit, 'massa', buit, 'angleVent', 0);
if plomaComplexa
    for i = 1:round(300 / tempsPas)
        ploma = avancar_ploma(ploma, par, tempsPas);
    end
end

% --- Estat inicial del cercador ---
if usaInfotaxis
    comptesModel = model_infotaxis(tipusPloma, par, graella, tempsPas);   % el que el dron espera mesurar
    probFuga = ones(graella.ny, graella.nx, graella.nz);                  % probabilitat de on es la fuga: uniforme
    probFuga(zona_visible(celDron, graella)) = 0;
    probFuga = probFuga / sum(probFuga(:));
else
    % Sensors en creu: est, oest, nord, sud i, en 3D, un a sobre i un a sota
    desplSensors = par.separacioSensors * [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
    desplSensors = desplSensors(1:4 + 2 * tresD, :);
    grad.filtrat           = [0 0 0];   % gradient suavitzat [dc/dx dc/dy dc/dz]
    grad.passosSenseSenyal = Inf;       % passos seguits sense detectar
    grad.centreZigZag      = celDron;   % cel.la al voltant de la qual busca quan no te senyal
    grad.tram              = 0;         % numero de tram del zig-zag (0 = encara no ha comencat)
    grad.desti             = celDron;   % cel.la on acaba el tram actual
    mapaMesures = zeros(graella.ny, graella.nx, graella.nz);   % concentracio maxima mesurada a cada cel.la
end

trajectoria  = nan(par.maxPassos + 1, 3);       % posicio del dron [x y z] a cada pas
trajectoria(1, :) = [graella.x(celDron(2)) graella.y(celDron(1)) graella.z(celDron(3))];
concMesurada = zeros(par.maxPassos, 1);         % concentracio mesurada a cada pas
trobada = false;

if mostrarAnimacio
    fig = crear_figura(par, graella, tresD, usaInfotaxis);
end

%% ============================== BUCLE ==============================
for pas = 1:par.maxPassos
    posDron = trajectoria(pas, :);              % posicio del dron [x y z] (m)

    % 1. El dron mesura i decideix on va
    if usaInfotaxis
        concReal = concentracio(ploma, par, posDron(1), posDron(2), posDron(3));
        comptes  = poisson_aleatori(tempsPas * (par.sensibilitat * concReal + par.fonsSensor));
        probFuga = infotaxis_actualitzar(probFuga, comptesModel, celDron, comptes, graella);
        [moviment, queFa] = infotaxis_decidir(probFuga, comptesModel, celDron, graella, par.intervalsComptes);
    else
        posSensors = posDron + desplSensors;
        concReal = concentracio(ploma, par, posSensors(:, 1), posSensors(:, 2), posSensors(:, 3));
        comptes  = poisson_aleatori(tempsPas * (par.sensibilitat * concReal + par.fonsSensor));
        [moviment, queFa, grad] = gradient_decidir(grad, comptes, celDron, graella, par, tempsPas);
    end
    % Concentracio que el dron dedueix dels comptes (mitjana dels sensors)
    concMesurada(pas) = mean(max(comptes / tempsPas - par.fonsSensor, 0) / par.sensibilitat);
    if ~usaInfotaxis
        mapaMesures(celDron(1), celDron(2), celDron(3)) = ...
            max(mapaMesures(celDron(1), celDron(2), celDron(3)), concMesurada(pas));
    end

    % 2. El dron es mou a una cel.la veina (o es queda quiet)
    celDron = celDron + moviment;
    posDron = [graella.x(celDron(2)) graella.y(celDron(1)) graella.z(celDron(3))];
    trajectoria(pas + 1, :) = posDron;
    % La fuga es troba si es a tocar (bloc de 3 x 3 cel.les) i a la mateixa alcada
    trobada = all(abs(celDron(1:2) - celFuga(1:2)) <= 1) && celDron(3) == celFuga(3);

    % 3. La ploma evoluciona (la gaussiana no canvia)
    if plomaComplexa
        ploma = avancar_ploma(ploma, par, tempsPas);
    end

    % 4. Dibuix
    if mostrarAnimacio && (trobada || mod(pas, dibuixaCadaNPassos) == 0)
        if ~ishandle(fig.finestra)
            break;
        end
        if usaInfotaxis
            coneixement = probFuga;
        else
            coneixement = mapaMesures;
        end
        titol = sprintf('%s | ploma %s | t = %.0f s | alcada %.0f m | %s', ...
                        upper(dimensio), tipusPloma, pas * tempsPas, posDron(3), queFa);
        actualitzar_figura(fig, par, graella, ploma, trajectoria(1:pas + 1, :), ...
                           coneixement, concMesurada(1:pas), tempsPas, titol);
    end
    if trobada
        break;
    end
end

%% ============================= RESULTAT =============================
recorregut = sum(sqrt(sum(diff(trajectoria(1:pas + 1, :)).^2, 2)));
if trobada
    resultat = sprintf('FUGA TROBADA en %.0f s (%.0f m recorreguts)', pas * tempsPas, recorregut);
else
    resultat = sprintf('Fuga NO trobada en %.0f s (distancia final %.0f m)', ...
                       pas * tempsPas, norm(trajectoria(pas + 1, :) - par.posFuga));
end
fprintf('%s, ploma %s, algorisme %s: %s\n', upper(dimensio), tipusPloma, algorisme, resultat);
if mostrarAnimacio && ishandle(fig.finestra)
    set(fig.titol, 'String', resultat);
end
end

function cel = posicio_a_cella(pos, par, graella)
%POSICIO_A_CELLA  Cel.la [iy ix iz] de la graella que correspon a la posicio [x y z].
[~, iz] = min(abs(graella.z - pos(3)));
cel = [round(pos(2) / par.midaCella) + 1, round(pos(1) / par.midaCella) + 1, iz];
end

%% ============================== PLOMA ==============================

function conc = concentracio(ploma, par, x, y, z)
%CONCENTRACIO  Concentracio real de la ploma als punts (x, y, z).
if strcmp(ploma.tipus, 'gaussiana')
    conc = ploma_gaussiana(x - par.posFuga(1), y - par.posFuga(2), z, par.posFuga(3), par);
else
    % Suma de tots els puffs. Cada puff es una campana de Gauss en 3D, mes
    % la seva imatge sota terra (el terra reflecteix el contaminant)
    pes        = ploma.massa ./ ((2 * pi)^1.5 * ploma.sigma2H .* sqrt(ploma.sigma2V));
    distHoritz = (x(:) - ploma.x').^2 + (y(:) - ploma.y').^2;             % (punts x puffs)
    vertical   = exp(-(z(:) - ploma.z').^2 ./ (2 * ploma.sigma2V')) + ...
                 exp(-(z(:) + ploma.z').^2 ./ (2 * ploma.sigma2V'));
    conc = reshape((exp(-distHoritz ./ (2 * ploma.sigma2H')) .* vertical) * pes, size(x));
end
end

function conc = ploma_gaussiana(aFavorVent, deTraves, z, alcadaFuga, par)
%PLOMA_GAUSSIANA  Ploma gaussiana estacionaria a una distancia aFavorVent
%   de la fuga, deTraves al vent i a l'alcada z, si la fuga es a alcadaFuga:
%
%                Q              n^2             (z-H)^2            (z+H)^2
%       c = ------------ exp(- ------) [ exp(- -------) + exp(- -------) ]
%           2 pi U sy sz       2 sy^2           2 sz^2             2 sz^2
%
%   sy i sz (amplada i gruix de la ploma) creixen amb la distancia segons
%   les formules de Briggs (1973) per a la classe d'estabilitat D. L'ultim
%   terme es la reflexio del contaminant al terra.
deTraves = deTraves + zeros(size(aFavorVent));
z        = z + zeros(size(aFavorVent));
conc     = zeros(size(aFavorVent));
arriba   = aFavorVent > -par.gaussMidaFuga;      % a sobrevent no hi arriba res
dist     = max(aFavorVent(arriba), 0);
amplada  = 0.08 * dist ./ sqrt(1 + 1e-4 * dist);
gruix    = 0.06 * dist ./ sqrt(1 + 1.5e-3 * dist);
amplada  = sqrt(amplada.^2 + par.gaussMidaFuga^2);   % la fuga te una mida inicial
gruix    = sqrt(gruix.^2 + par.gaussMidaFuga^2);
vertical = exp(-(z(arriba) - alcadaFuga).^2 ./ (2 * gruix.^2)) + ...
           exp(-(z(arriba) + alcadaFuga).^2 ./ (2 * gruix.^2));
conc(arriba) = par.gaussEmissio ./ (2 * pi * par.gaussVent * amplada .* gruix) ...
               .* exp(-deTraves(arriba).^2 ./ (2 * amplada.^2)) .* vertical;
end

function ploma = avancar_ploma(ploma, par, tempsPas)
%AVANCAR_PLOMA  Un pas de temps de la ploma complexa.
% La direccio del vent oscil.la a l'atzar (proces d'Ornstein-Uhlenbeck)
memoria = exp(-tempsPas / par.meandreTemps);
ploma.angleVent = ploma.angleVent * memoria + (par.meandreGraus * pi / 180) * sqrt(1 - memoria^2) * randn;
% Els puffs es mouen amb el vent i a l'atzar, s'eixamplen i es van esvaint
saltH = sqrt(2 * par.turbulenciaH * tempsPas);
saltV = sqrt(2 * par.turbulenciaV * tempsPas);
ploma.x = ploma.x + par.vent * cos(ploma.angleVent) * tempsPas + saltH * randn(size(ploma.x));
ploma.y = ploma.y + par.vent * sin(ploma.angleVent) * tempsPas + saltH * randn(size(ploma.y));
ploma.z = abs(ploma.z + saltV * randn(size(ploma.z)));     % rebot al terra
ploma.sigma2H = ploma.sigma2H + 2 * par.creixementH * tempsPas;
ploma.sigma2V = ploma.sigma2V + 2 * par.creixementV * tempsPas;
ploma.massa   = ploma.massa * exp(-tempsPas / par.tempsVida);
% La fuga emet puffs nous
nous = ones(round(par.puffsPerSegon * tempsPas), 1);
ploma.x = [ploma.x; par.posFuga(1) * nous];
ploma.y = [ploma.y; par.posFuga(2) * nous];
ploma.z = [ploma.z; par.posFuga(3) * nous];
ploma.sigma2H = [ploma.sigma2H; par.midaPuff^2 * nous];
ploma.sigma2V = [ploma.sigma2V; par.midaPuff^2 * nous];
ploma.massa   = [ploma.massa;   par.emissio * tempsPas / numel(nous) * nous];
% Eliminem els puffs que han sortit de la zona
fora = ploma.x > par.midaX + 60 | ploma.x < -60 | ploma.y > par.midaY + 60 | ploma.y < -60;
ploma.x(fora) = [];        ploma.y(fora) = [];        ploma.z(fora) = [];
ploma.sigma2H(fora) = [];  ploma.sigma2V(fora) = [];  ploma.massa(fora) = [];
end

function comptes = poisson_aleatori(mitjana)
%POISSON_ALEATORI  Nombres aleatoris de Poisson amb la mitjana donada (sense toolboxes).
comptes = zeros(size(mitjana));
for i = 1:numel(mitjana)
    if mitjana(i) < 30
        % Metode de Knuth: multiplica uniformes fins a baixar de exp(-mitjana)
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

%% ============================ INFOTAXIS ============================

function comptesModel = model_infotaxis(tipusPloma, par, graella, tempsPas)
%MODEL_INFOTAXIS  Model intern del dron: comptes que espera mesurar en un
%   pas si ell es a l'alcada z(izDron), la fuga a l'alcada z(izFuga) i el
%   desplacament horitzontal entre tots dos es (distX, distY). Es una taula
%   comptesModel(iyDron - iyFuga + ny, ixDron - ixFuga + nx, izDron, izFuga).
%
%   - Ploma gaussiana: el dron en coneix la formula (model exacte).
%   - Ploma complexa: solucio mitjana de l'equacio d'adveccio-difusio de
%     l'article, amb difusivitats efectives horitzontal (Dh) i vertical
%     (Dv) i una font imatge sota terra:
%
%       c = Q / (4 pi sqrt(Dh Dv)) exp(x U / 2Dh) [e^(-r1/abast)/r1 + e^(-r2/abast)/r2]
%
%     r1 i r2 son les distancies a la fuga i a la seva imatge (amb z
%     escalada per sqrt(Dh/Dv)). El dron NO coneix la turbulencia real.
[distX, distY] = meshgrid((-(graella.nx - 1):(graella.nx - 1)) * par.midaCella, ...
                          (-(graella.ny - 1):(graella.ny - 1)) * par.midaCella);
difusioH = par.turbulenciaH + par.creixementH + (par.vent * par.meandreGraus * pi / 180)^2 * par.meandreTemps;
difusioV = par.turbulenciaV + par.creixementV;
abast    = sqrt(difusioH * par.tempsVida / (1 + par.vent^2 * par.tempsVida / (4 * difusioH)));
comptesModel = zeros(2 * graella.ny - 1, 2 * graella.nx - 1, graella.nz, graella.nz);
for izDron = 1:graella.nz
    for izFuga = 1:graella.nz
        zDron = graella.z(izDron);
        zFuga = graella.z(izFuga);
        if strcmp(tipusPloma, 'gaussiana')
            conc = ploma_gaussiana(distX, distY, zDron, zFuga, par);
        else
            distFuga   = max(sqrt(distX.^2 + distY.^2 + difusioH / difusioV * (zDron - zFuga)^2), par.midaCella / 2);
            distImatge = max(sqrt(distX.^2 + distY.^2 + difusioH / difusioV * (zDron + zFuga)^2), par.midaCella / 2);
            conc = par.emissio / (4 * pi * sqrt(difusioH * difusioV)) * ...
                   (exp(distX * par.vent / (2 * difusioH) - distFuga / abast) ./ distFuga + ...
                    exp(distX * par.vent / (2 * difusioH) - distImatge / abast) ./ distImatge);
        end
        comptesModel(:, :, izDron, izFuga) = tempsPas * (par.sensibilitat * conc + par.fonsSensor);
    end
end
end

function esperats = comptes_esperats(comptesModel, cel, graella)
%COMPTES_ESPERATS  Comptes que el dron espera mesurar a la cel.la cel per a
%   cada possible posicio de la fuga a la graella (matriu ny x nx x nz).
esperats = reshape(comptesModel(cel(1) + graella.ny - 1:-1:cel(1), ...
                                cel(2) + graella.nx - 1:-1:cel(2), cel(3), :), ...
                   graella.ny, graella.nx, graella.nz);
end

function visible = zona_visible(cel, graella)
%ZONA_VISIBLE  Cel.les a tocar de cel (bloc de 3 x 3 a la mateixa alcada):
%   si la fuga hi es, el dron la veu.
visible = false(graella.ny, graella.nx, graella.nz);
visible(max(cel(1) - 1, 1):min(cel(1) + 1, graella.ny), ...
        max(cel(2) - 1, 1):min(cel(2) + 1, graella.nx), cel(3)) = true;
end

function incertesa = entropia(prob)
%ENTROPIA  Entropia de Shannon -sum p ln p: mesura la incertesa d'un mapa
%   de probabilitat (0 = se sap del cert on es la fuga).
prob = prob(prob > 0);
incertesa = -sum(prob .* log(prob));
end

function probFuga = infotaxis_actualitzar(probFuga, comptesModel, celDron, comptes, graella)
%INFOTAXIS_ACTUALITZAR  Regla de Bayes (eq. 2 de l'article) despres de
%   mesurar: cada possible posicio de la fuga es multiplica per la
%   probabilitat (de Poisson) d'haver mesurat aquests comptes si la fuga
%   fos alla. La zona que el dron veu queda descartada (la fuga no hi es).
esperats = comptes_esperats(comptesModel, celDron, graella);
probNova = probFuga .* exp(-esperats + comptes * log(esperats) - gammaln(comptes + 1));
probNova(zona_visible(celDron, graella)) = 0;
if sum(probNova(:)) > 0
    probFuga = probNova / sum(probNova(:));
else
    probFuga(zona_visible(celDron, graella)) = 0;    % mesura incompatible amb el model: nomes descartem la zona
    probFuga = probFuga / sum(probFuga(:));
end
end

function [moviment, queFa] = infotaxis_decidir(probFuga, comptesModel, celDron, graella, intervals)
%INFOTAXIS_DECIDIR  Tria el moviment que fa baixar mes l'entropia esperada
%   (eq. 1 de l'article):
%
%       canvi(r -> rj) = P(rj) [-S] + [1 - P(rj)] sum_k rho_k(rj) dS_k
%
%   Primer terme: probabilitat de trobar la fuga a rj (l'entropia passa a 0).
%   Segon terme: si no hi es, el dron mesurara k comptes amb probabilitat
%   rho_k i l'entropia canviara dS_k. Els comptes s'agrupen en intervals.
entropiaAra   = entropia(probFuga);
canviEntropia = Inf(7, 1);              % canvi esperat d'entropia per a cada moviment
for m = 1:7
    celNova = celDron + graella.moviments(m, :);
    if any(celNova < 1) || any(celNova > [graella.ny graella.nx graella.nz])
        continue;                       % fora de la zona (en 2D, amunt i avall)
    end
    visible    = zona_visible(celNova, graella);
    probTrobar = sum(probFuga(visible));    % probabilitat de trobar la fuga en arribar-hi
    probSiNo   = probFuga;
    probSiNo(visible) = 0;                  % mapa de probabilitat si NO la troba
    if sum(probSiNo(:)) <= 0
        canviEntropia(m) = -entropiaAra;
        continue;
    end
    probSiNo = probSiNo / sum(probSiNo(:));
    esperats = comptes_esperats(comptesModel, celNova, graella);
    possible = probSiNo > 1e-15;            % nomes compten les posicions encara possibles (mes rapid)
    probSiNo = probSiNo(possible);
    esperats = esperats(possible);

    % Recorrem k = 0, 1, 2... acumulant la llei de Poisson per intervals
    canviMesura   = 0;
    probAcumulada = zeros(size(esperats));  % P(comptes < limit actual)
    probK         = exp(-esperats);         % P(k comptes) per al k actual
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
            probInterval = max(1 - probAcumulada, 0);   % ultim interval: k >= intervals(end)
        end
        conjunta   = probSiNo .* probInterval;
        probMesura = sum(conjunta(:));      % probabilitat de mesurar comptes d'aquest interval
        if probMesura > 1e-300
            canviMesura = canviMesura + probMesura * (entropia(conjunta / probMesura) - entropiaAra);
        end
    end
    canviEntropia(m) = probTrobar * (-entropiaAra) + (1 - probTrobar) * canviMesura;
end

if min(canviEntropia) > -1e-6
    % Cap moviment aporta informacio (el dron ja "sap" on es la fuga):
    % vola directament cap al punt mes probable
    queFa = 'cap al punt mes probable';
    [~, mesProbable] = max(probFuga(:));
    [iy, ix, iz] = ind2sub([graella.ny graella.nx graella.nz], mesProbable);
    falta = [iy ix iz] - celDron;
    [~, eix] = max(abs(falta));
    moviment = [0 0 0];
    moviment(eix) = sign(falta(eix));
else
    queFa = 'infotaxis';
    [~, millor] = min(canviEntropia + 1e-10 * abs(entropiaAra) * rand(7, 1));    % desempat aleatori
    moviment = graella.moviments(millor, :);
end
end

%% ============================= GRADIENT =============================

function [moviment, queFa, grad] = gradient_decidir(grad, comptes, celDron, graella, par, tempsPas)
%GRADIENT_DECIDIR  Estima el gradient amb els comptes dels sensors (est,
%   oest, nord, sud i, en 3D, amunt i avall) per diferencies finites
%   centrades:
%
%       dc/dx ~ (cE - cO) / (2d),   dc/dy ~ (cN - cS) / (2d),   dc/dz ~ (cAmunt - cAvall) / (2d)
%
%   i el suavitza amb un filtre exponencial. Si hi ha senyal recent, el dron
%   va cap on el gradient es mes gran, pero mai vent avall (sap d'on ve el
%   vent, i la fuga sempre es vent amunt). Si perd el senyal, busca la
%   ploma en zig-zag de traves al vent al voltant d'on l'ha perduda (com
%   fan les arnes): cada tram es el doble de llarg que l'anterior i, en
%   3D, es fa a una alcada diferent.
conc = max(comptes(:)' / tempsPas - par.fonsSensor, 0) / par.sensibilitat;   % concentracio a cada sensor
gradientNou = [0 0 0];
gradientNou(1:numel(conc) / 2) = (conc(1:2:end) - conc(2:2:end)) / (2 * par.separacioSensors);
if sum(comptes) >= par.llindarDeteccio
    grad.passosSenseSenyal = 0;
    grad.filtrat = par.suavitzat * grad.filtrat + (1 - par.suavitzat) * gradientNou;
    grad.centreZigZag = celDron;        % si perd el senyal, el buscara al voltant d'aqui
    grad.tram = 0;
else
    grad.passosSenseSenyal = grad.passosSenseSenyal + 1;
    if grad.passosSenseSenyal > par.passosMemoria
        grad.filtrat = [0 0 0];         % senyal perdut: oblida el gradient
    end
end

limits = [graella.ny graella.nx graella.nz];
if grad.passosSenseSenyal <= par.passosMemoria
    queFa = 'seguint el gradient';
    veines    = celDron + graella.moviments(1:6, :);
    dinsZona  = all(veines >= 1, 2) & all(veines <= limits, 2);
    projeccio = graella.moviments(1:6, [2 1 3]) * grad.filtrat' + 1e-12 * rand(6, 1);   % gradient en cada direccio
    projeccio(~dinsZona) = -Inf;
    projeccio(1) = -Inf;                % mai vent avall: la fuga sempre es vent amunt
    [~, millor] = max(projeccio);
    moviment = graella.moviments(millor, :);
else
    queFa = 'buscant la ploma en zig-zag';
    if grad.tram == 0 || isequal(celDron, grad.desti)
        % Tram nou: cap a l'altre costat, el doble de lluny i a una altra alcada
        grad.tram = grad.tram + 1;
        amplitud  = par.tramInicial * 2^(grad.tram - 1);
        costat    = (-1)^grad.tram;
        desti     = grad.centreZigZag + [costat * amplitud, 0, costat * floor(grad.tram / 2)];
        grad.desti = min(max(desti, 1), limits);
    end
    falta = grad.desti - celDron;       % cel.les que falten per arribar al final del tram
    [~, eix] = max(abs(falta));
    moviment = [0 0 0];
    moviment(eix) = sign(falta(eix));
end
end

%% ============================== FIGURA ==============================

function fig = crear_figura(par, graella, tresD, usaInfotaxis)
%CREAR_FIGURA  Crea la finestra amb tots els panells i en retorna els elements
%   grafics, que despres s'actualitzen a actualitzar_figura:
%     - escena real (ploma, fuga, dron i trajectoria) i, en 3D, vista lateral
%     - el que sap el dron, vist des de dalt i, en 3D, de costat
%     - concentracio mesurada al llarg del temps
fig.finestra = figure('Color', 'w', 'Name', 'Dron: cerca de l''origen d''una fuga', ...
                      'Units', 'normalized', 'Position', [0.03 0.08 0.94 0.84]);
colormap(flipud(hot));
fig.tresD = tresD;
fig.usaInfotaxis = usaInfotaxis;
if tresD
    llocEscena = [0.05 0.44 0.42 0.48];  llocSaber = [0.56 0.56 0.40 0.36];  llocMesura = [0.56 0.07 0.40 0.14];
else
    llocEscena = [0.05 0.12 0.42 0.76];  llocSaber = [0.56 0.46 0.40 0.46];  llocMesura = [0.56 0.10 0.40 0.24];
end
if usaInfotaxis
    titolSaber = 'El que creu el dron: probabilitat de on es la fuga';
    escalaSaber = [0 1];
else
    titolSaber = 'El que sap el dron: concentracio que ha mesurat a cada lloc';
    escalaSaber = [-2 2.5];
end
estilDron = {'o', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k'};
estilFuga = {'p', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k'};

% Escena real. Al terra es pinta la ploma a l'alcada del dron.
subplot('Position', llocEscena);
fig.terra = surf(graella.x, graella.y, zeros(graella.ny, graella.nx), zeros(graella.ny, graella.nx), ...
                 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;  caxis([-2 2.5]);
plot3(par.posFuga(1) * [1 1], par.posFuga(2) * [1 1], [0 par.posFuga(3)], 'g-', 'LineWidth', 1.5);
plot3(par.posFuga(1), par.posFuga(2), par.posFuga(3), estilFuga{:});
fig.trajectoria = plot3(nan, nan, nan, 'b-', 'LineWidth', 1.5);
fig.palDron     = plot3([nan nan], [nan nan], [nan nan], 'k-');     % linia del dron fins a terra
fig.dron        = plot3(nan, nan, nan, estilDron{:});
axis([0 par.midaX 0 par.midaY 0 par.alcadaMaxima]);
xlabel('x (m)');  ylabel('y (m)');  zlabel('z (m)');
if tresD
    daspect([1 1 0.25]);  view(-30, 35);  grid on;
else
    daspect([1 1 1]);  view(2);  grid off;
end
fig.titol = title('');

% El que sap el dron, vist des de dalt
subplot('Position', llocSaber);
fig.saberDalt = imagesc(graella.x, graella.y, zeros(graella.ny, graella.nx));
set(gca, 'YDir', 'normal');
axis equal;  axis([0 par.midaX 0 par.midaY]);  caxis(escalaSaber);
hold on;
fig.dronSaberDalt = plot(nan, nan, estilDron{:});
xlabel('x (m)');  ylabel('y (m)');
title([titolSaber ' (des de dalt)']);

if tresD
    % Vista lateral de l'escena: tall vertical de la ploma que passa pel dron
    fig.zLateral = 0:par.alcadaMaxima;
    [fig.XLateral, fig.ZLateral] = meshgrid(graella.x, fig.zLateral);
    subplot('Position', [0.05 0.08 0.42 0.22]);
    fig.tallLateral = imagesc(graella.x, fig.zLateral, zeros(size(fig.XLateral)));
    set(gca, 'YDir', 'normal');
    axis([0 par.midaX 0 par.alcadaMaxima]);  caxis([-2 2.5]);
    hold on;
    plot(par.posFuga(1), par.posFuga(3), estilFuga{:});
    fig.trajLateral = plot(nan, nan, 'b-', 'LineWidth', 1.5);
    fig.dronLateral = plot(nan, nan, estilDron{:});
    xlabel('x (m)');  ylabel('z (m)');
    title('Vista lateral: ploma en el pla vertical que passa pel dron');

    % El que sap el dron, vist de costat
    subplot('Position', [0.56 0.31 0.40 0.13]);
    fig.saberCostat = imagesc(graella.x, graella.z, zeros(graella.nz, graella.nx));
    set(gca, 'YDir', 'normal');
    axis([0 par.midaX 0 par.alcadaMaxima + par.midaCella / 2]);  caxis(escalaSaber);
    hold on;
    fig.dronSaberCostat = plot(nan, nan, estilDron{:});
    ylabel('z (m)');
    title('El mateix, vist de costat');
end

% Concentracio mesurada al llarg del temps (escala logaritmica)
subplot('Position', llocMesura);
fig.mesura = plot(nan, nan, 'r-');
set(gca, 'YScale', 'log');
ylim([0.05 1000]);  grid on;
xlabel('Temps (s)');
title('Concentracio mesurada pel dron');
end

function actualitzar_figura(fig, par, graella, ploma, trajectoria, coneixement, concMesurada, tempsPas, titol)
%ACTUALITZAR_FIGURA  Redibuixa tots els panells amb l'estat actual.
%   coneixement es el que sap el dron (matriu ny x nx x nz): la
%   probabilitat de on es la fuga (infotaxis) o la concentracio que ha
%   mesurat a cada cel.la (gradient).
posDron = trajectoria(end, :);

% Escena real
set(fig.terra, 'CData', log10(concentracio(ploma, par, graella.X, graella.Y, posDron(3)) + 1e-3));
set(fig.trajectoria, 'XData', trajectoria(:, 1), 'YData', trajectoria(:, 2), 'ZData', trajectoria(:, 3));
set(fig.palDron, 'XData', posDron(1) * [1 1], 'YData', posDron(2) * [1 1], 'ZData', [0 posDron(3)]);
set(fig.dron, 'XData', posDron(1), 'YData', posDron(2), 'ZData', posDron(3));
set(fig.titol, 'String', {titol, 'Colors del terra: concentracio a l''alcada del dron'});

% El que sap el dron: des de dalt i de costat
if fig.usaInfotaxis
    saberDalt   = sum(coneixement, 3);                   % probabilitat sumada sobre totes les alcades
    saberCostat = reshape(sum(coneixement, 1), graella.nx, graella.nz)';   % ... i sobre totes les y
    saberDalt   = saberDalt / max(saberDalt(:));
    saberCostat = saberCostat / max(saberCostat(:));
else
    saberDalt   = log10(max(coneixement, [], 3) + 1e-3);
    saberCostat = log10(reshape(max(coneixement, [], 1), graella.nx, graella.nz)' + 1e-3);
end
set(fig.saberDalt, 'CData', saberDalt);
set(fig.dronSaberDalt, 'XData', posDron(1), 'YData', posDron(2));

if fig.tresD
    set(fig.tallLateral, 'CData', log10(concentracio(ploma, par, fig.XLateral, posDron(2), fig.ZLateral) + 1e-3));
    set(fig.trajLateral, 'XData', trajectoria(:, 1), 'YData', trajectoria(:, 3));
    set(fig.dronLateral, 'XData', posDron(1), 'YData', posDron(3));
    set(fig.saberCostat, 'CData', saberCostat);
    set(fig.dronSaberCostat, 'XData', posDron(1), 'YData', posDron(3));
end

% Concentracio mesurada (els zeros es dibuixen a la part de baix de l'escala)
set(fig.mesura, 'XData', (1:numel(concMesurada)) * tempsPas, 'YData', max(concMesurada, 0.05));
drawnow;
end
