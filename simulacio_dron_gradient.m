function simulacio_dron_gradient()
%SIMULACIO_DRON_GRADIENT  Versio preliminar del projecte "Dron autonom per
%   trobar l'origen d'una fuga".
%
%   Simula en 2D una ploma turbulenta de contaminant (model de "puffs") i
%   un dron que en busca l'origen seguint el gradient de concentracio (la
%   direccio en que la concentracio augmenta mes), mesurat amb 4 sensors
%   en creu.
%
%   Executa-la amb F5. Els parametres son just a sota.
%
%   En 20 proves amb llavors diferents el dron troba la fuga totes les
%   vegades: la meitat en menys de 5 minuts i la resta en fins a 21. Quan
%   triga, es perque passa de llarg la fuga, perd la ploma i l'ha de
%   retrobar explorant a l'atzar (limitacio coneguda del metode del gradient).

%% ============================ PARAMETRES ============================
llavor  = 5;                % mateix numero = mateixa simulacio ([] = cada cop diferent)
mostrar = true;             % false = sense animacio (mes rapid)

% Escenari (unitats: metres i segons). El vent bufa cap a +x.
P.L     = [400 300];        % mida de la zona [x y]
P.font  = [60 150];         % on es la fuga (el dron NO ho sap)
P.inici = [340 200];        % d'on surt el dron
P.U     = 2;                % velocitat mitjana del vent (m/s)
P.Q     = 200;              % intensitat de la fuga (massa emesa per segon)
P.tau   = 200;              % temps de vida del contaminant (s)

% Ploma: la fuga emet "puffs" (nuvolets gaussians) que el vent s'emporta
P.puffsPerSeg  = 2;         % puffs emesos per segon
P.sigma0       = 1.5;       % mida inicial d'un puff (m)
P.Dpuff        = 0.4;       % rapidesa amb que s'eixampla un puff (m^2/s)
P.Dturb        = 0.4;       % moviment aleatori dels puffs per la turbulencia (m^2/s)
P.sigmaMeandre = 15;        % oscil.lacio de la direccio del vent (graus)
P.Tmeandre     = 40;        % temps tipic d'aquesta oscil.lacio (s)

% Dron i sensors
vDron    = 5;               % velocitat del dron (m/s)
dt       = 1;               % pas de temps (s): a cada pas el dron avanca vDron*dt = 5 m
d        = 5;               % distancia de cada sensor al centre del dron (m)
soroll   = 0.05;            % soroll dels sensors (desviacio tipica)
llindar  = 0.1;             % concentracio a partir de la qual es considera que detecta
alfa     = 0.95;            % suavitzat del gradient: ~mitjana dels ultims 1/(1-alfa) = 20 passos
memoria  = 8;               % passos sense detectar res abans de tornar a explorar
pGir     = 0.15;            % en exploracio, probabilitat de girar a cada pas
radiFont = 7.5;             % a aquesta distancia el dron troba la fuga (m)
tMax     = 2500;            % temps maxim de vol (s)
cadaN    = 2;               % dibuixa cada N passos
%% ====================================================================

if ~isempty(llavor)
    rng(llavor);
end

% Direccions de moviment del dron: est, nord, oest, sud (en aquest ordre,
% aixi girar 90 graus es sumar o restar 1)
dirs = [1 0; 0 1; -1 0; 0 -1];

% --- Ploma inicial: la deixem evolucionar 300 s abans que surti el dron ---
pl = struct('x', zeros(0, 1), 'y', zeros(0, 1), 's2', zeros(0, 1), 'm', zeros(0, 1), 'th', 0);
for i = 1:300
    pl = avancar_ploma(pl, P, dt);
end

% --- Estat del dron ---
x  = P.inici;               % posicio
gf = [0 0];                 % gradient filtrat
senseSenyal = Inf;          % passos seguits sense detectar
k  = randi(4);              % direccio actual (fila de dirs)
trobada = false;
nMax = round(tMax / dt);
traj = nan(nMax + 1, 2);  traj(1, :) = x;
cMes = nan(nMax + 1, 1);
dist = nan(nMax + 1, 1);  dist(1) = norm(x - P.font);

% --- Figura ---
if mostrar
    [xv, yv] = meshgrid(0:5:P.L(1), 0:5:P.L(2));
    fig = figure('Color', 'w', 'Name', 'Dron: cerca per gradient', ...
                 'Units', 'normalized', 'Position', [0.05 0.1 0.9 0.75]);
    colormap(flipud(hot));
    subplot('Position', [0.05 0.1 0.5 0.8]);
    hImg = imagesc(0:5:P.L(1), 0:5:P.L(2), camp(pl, xv, yv));
    set(gca, 'YDir', 'normal');
    axis equal;
    axis([0 P.L(1) 0 P.L(2)]);
    caxis([-2 1.5]);
    cb = colorbar;
    ylabel(cb, 'log_{10} concentracio');
    hold on;
    plot(P.font(1), P.font(2), 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    hTraj = plot(x(1), x(2), 'b-', 'LineWidth', 1.5);
    hDron = plot(x(1), x(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');
    xlabel('x (m)'); ylabel('y (m)');
    hTit = title('');
    subplot('Position', [0.68 0.58 0.29 0.32]);
    hConc = plot(0, 0, 'r-');
    title('Concentracio mesurada (mitjana dels 4 sensors)'); grid on;
    subplot('Position', [0.68 0.1 0.29 0.32]);
    hDist = plot(0, dist(1), 'b-');
    title('Distancia a la fuga (m)'); xlabel('Temps (s)'); grid on;
end

%% ============================ BUCLE ============================
for n = 1:nMax
    t = n * dt;

    % 1. La ploma evoluciona
    pl = avancar_ploma(pl, P, dt);

    % 2. Mesura dels 4 sensors (est, oest, nord, sud) amb soroll
    sensors = x + [d 0; -d 0; 0 d; 0 -d];
    c = concentracio(sensors, pl);
    c = max(c + soroll * randn(4, 1), 0);
    cMes(n + 1) = mean(c);

    % 3. Gradient per diferencies finites centrades, suavitzat
    gNou = [c(1) - c(2), c(3) - c(4)] / (2 * d);
    if max(c) > llindar
        senseSenyal = 0;
        gf = alfa * gf + (1 - alfa) * gNou;
    else
        senseSenyal = senseSenyal + 1;
    end

    % 4. Decisio: seguir el gradient o explorar
    if senseSenyal <= memoria
        mode = 'seguint el gradient';
        [~, k] = max(dirs * gf');                   % direccio on el gradient es mes gran
    else
        mode = 'explorant';
        if rand < pGir
            k = mod(k - 1 + 2 * randi(2) - 3, 4) + 1;   % gira 90 graus a esquerra o dreta
        end
    end

    % 5. Moviment (si sortiria de la zona, fa mitja volta)
    xNou = x + vDron * dt * dirs(k, :);
    if any(xNou < 0) || any(xNou > P.L)
        k = mod(k + 1, 4) + 1;                      % direccio contraria
        xNou = x + vDron * dt * dirs(k, :);
    end
    x = xNou;
    traj(n + 1, :) = x;
    dist(n + 1) = norm(x - P.font);

    % 6. Ha arribat a la fuga?
    trobada = dist(n + 1) < radiFont;

    % 7. Dibuix
    if mostrar && (trobada || mod(n, cadaN) == 0)
        if ~ishandle(fig)
            break;
        end
        set(hImg, 'CData', camp(pl, xv, yv));
        set(hTraj, 'XData', traj(1:n + 1, 1), 'YData', traj(1:n + 1, 2));
        set(hDron, 'XData', x(1), 'YData', x(2));
        set(hConc, 'XData', (0:n) * dt, 'YData', cMes(1:n + 1));
        set(hDist, 'XData', (0:n) * dt, 'YData', dist(1:n + 1));
        set(hTit, 'String', sprintf('t = %.0f s | %s', t, mode));
        drawnow;
    end
    if trobada
        break;
    end
end

%% ============================ RESULTAT ============================
recorregut = sum(sqrt(sum(diff(traj(1:n + 1, :)).^2, 2)));
if trobada
    txt = sprintf('FONT TROBADA en %.0f s (%.0f m recorreguts)', n * dt, recorregut);
else
    txt = sprintf('Font NO trobada en %.0f s (distancia final %.0f m)', n * dt, dist(n + 1));
end
fprintf('%s\n', txt);
if mostrar && ishandle(fig)
    set(hTit, 'String', txt);
end
end

%% ======================== FUNCIONS AUXILIARS ========================

function c = concentracio(punts, pl)
%CONCENTRACIO  Concentracio als punts [x y] (una fila per punt): suma de
%   tots els puffs, cadascun amb forma de campana de Gauss.
if isempty(pl.x)
    c = zeros(size(punts, 1), 1);
    return;
end
dist2 = (punts(:, 1) - pl.x').^2 + (punts(:, 2) - pl.y').^2;   % (punts x puffs)
c = exp(-dist2 ./ (2 * pl.s2')) * (pl.m ./ (2 * pi * pl.s2));
end

function pl = avancar_ploma(pl, P, dt)
%AVANCAR_PLOMA  Un pas de temps de la ploma.
% La direccio del vent oscil.la a l'atzar (proces d'Ornstein-Uhlenbeck)
a = exp(-dt / P.Tmeandre);
pl.th = pl.th * a + (P.sigmaMeandre * pi / 180) * sqrt(1 - a^2) * randn;
u = P.U * cos(pl.th);
v = P.U * sin(pl.th);
% Els puffs es mouen amb el vent i a l'atzar, s'eixamplen i es van esvaint
sd = sqrt(2 * P.Dturb * dt);
pl.x  = pl.x + u * dt + sd * randn(size(pl.x));
pl.y  = pl.y + v * dt + sd * randn(size(pl.y));
pl.s2 = pl.s2 + 2 * P.Dpuff * dt;
pl.m  = pl.m * exp(-dt / P.tau);
% La fuga emet puffs nous
nNous = round(P.puffsPerSeg * dt);
pl.x  = [pl.x;  P.font(1) * ones(nNous, 1)];
pl.y  = [pl.y;  P.font(2) * ones(nNous, 1)];
pl.s2 = [pl.s2; P.sigma0^2 * ones(nNous, 1)];
pl.m  = [pl.m;  P.Q * dt / nNous * ones(nNous, 1)];
% Eliminem els puffs que han sortit de la zona
fora = pl.x > P.L(1) + 50 | pl.x < -50 | pl.y > P.L(2) + 50 | pl.y < -50;
pl.x(fora) = []; pl.y(fora) = []; pl.s2(fora) = []; pl.m(fora) = [];
end

function C = camp(pl, xv, yv)
%CAMP  Concentracio (log10) en una graella, nomes per dibuixar-la.
C = reshape(concentracio([xv(:) yv(:)], pl), size(xv));
C = log10(C + 1e-3);
end
