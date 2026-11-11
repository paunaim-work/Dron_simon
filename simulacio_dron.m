function simulacio_dron()
%SIMULACIO_DRON  Dron autonom que busca l'origen d'una fuga (versio simplificada).
%
%   Simula en 2D una fuga de contaminant i un dron que en busca l'origen
%   nomes amb les lectures del seu sensor. Es pot triar:
%
%     - la ploma:    'gaussiana' (estacionaria, model de Pasquill-Gifford)
%                    'complexa'  (turbulenta: "puffs" que el vent s'emporta)
%     - l'algorisme: 'infotaxis' (Vergassola, Villermaux i Shraiman, Nature 2007)
%                    'gradient'  (puja pel gradient de concentracio)
%
%   Executa-la amb F5. Les opcions i els parametres son just a sota.

%% ============================= OPCIONS =============================
tipusPloma = 'complexa';    % 'gaussiana' o 'complexa'
algorisme  = 'infotaxis';   % 'infotaxis' o 'gradient'
llavor     = 12;            % mateix numero = mateixa simulacio ([] = cada cop diferent)
mostrar    = true;          % false = sense animacio (mes rapid)
cadaN      = 2;             % dibuixa cada N passos

%% ============================ PARAMETRES ============================
% Unitats: metres i segons. El vent bufa cap a +x.
p.Lx = 400;  p.Ly = 300;    % mida de la zona de cerca
p.dx = 5;                   % mida de cel.la = el que avanca el dron a cada pas
p.font  = [60 150];         % on es la fuga (el dron NO ho sap)
p.inici = [340 215];        % d'on surt el dron
p.vDron = 5;                % velocitat del dron (m/s)
p.maxPassos = 2500;         % limit de passos (autonomia)

% Sensor: comptador de deteccions (tipus Geiger).
% Comptes en un pas ~ Poisson( dt * (sens * concentracio + fons) )
p.sens = 1;                 % comptes/s per unitat de concentracio
p.fons = 0.02;              % comptes/s de fons (falses deteccions)

% Ploma gaussiana
p.Qg = 5000;                % taxa d'emissio
p.Ug = 3;                   % velocitat del vent (m/s)
p.H  = 10;                  % alcada de la fuga (m)
p.z  = 10;                  % alcada de vol del dron (m)
p.sigma0g = 1;              % mida de la fuga (m)

% Ploma complexa
p.Q   = 200;                % taxa d'emissio
p.U   = 2;                  % velocitat mitjana del vent (m/s)
p.tau = 200;                % temps de vida del contaminant (s)
p.puffsPerSeg  = 2;         % puffs emesos per segon
p.sigma0       = 1.5;       % mida inicial d'un puff (m)
p.Dpuff        = 0.4;       % rapidesa amb que s'eixampla un puff (m^2/s)
p.Dturb        = 0.4;       % moviment aleatori dels puffs (m^2/s)
p.sigmaMeandre = 15;        % oscil.lacio de la direccio del vent (graus)
p.Tmeandre     = 40;        % temps tipic d'aquesta oscil.lacio (s)

% Infotaxis
p.bins = [0 1 2 3 5 8 13 21 34 55 89];   % intervals de comptes per planificar

% Gradient
p.d       = 5;              % distancia dels 4 sensors al centre del dron (m)
p.alfa    = 0.6;            % suavitzat del gradient (0 = cap, 0.9 = molt)
p.llindar = 2;              % comptes (suma dels 4 sensors) per considerar que detecta
p.memoria = 8;              % passos sense detectar abans de tornar a explorar
p.pGir    = 0.15;           % en exploracio, probabilitat de girar a cada pas
%% ====================================================================

if ~isempty(llavor)
    rng(llavor);
end
infotaxis = strcmp(algorisme, 'infotaxis');
complexa  = strcmp(tipusPloma, 'complexa');

% --- Graella: les matrius tenen mida Ny x Nx i s'indexen (iy, ix) ---
g.x  = 0:p.dx:p.Lx;
g.y  = 0:p.dx:p.Ly;
g.Nx = numel(g.x);
g.Ny = numel(g.y);
[g.X, g.Y] = meshgrid(g.x, g.y);
g.movs = [0 1; 1 0; 0 -1; -1 0; 0 0];    % [diy dix]: est, nord, oest, sud, quiet
dt = p.dx / p.vDron;                     % durada d'un pas (s)

celFont = round(p.font([2 1])  / p.dx) + 1;   % cel.la [iy ix] de la fuga
pos     = round(p.inici([2 1]) / p.dx) + 1;   % cel.la [iy ix] del dron

% --- Ploma: la complexa evoluciona 300 s abans que surti el dron ---
pl = struct('tipus', tipusPloma, 'x', zeros(0, 1), 'y', zeros(0, 1), ...
            's2', zeros(0, 1), 'm', zeros(0, 1), 'th', 0);
if complexa
    for i = 1:round(300 / dt)
        pl = avancar_ploma(pl, p, dt);
    end
end

% --- Estat inicial del cercador ---
if infotaxis
    K = model_infotaxis(tipusPloma, p, g, dt);   % comptes esperats segons el model del dron
    P = ones(g.Ny, g.Nx);                        % probabilitat de on es la fuga: uniforme
    P(zona(pos, g)) = 0;
    P = P / sum(P(:));
else
    G.gf    = [0 0];        % gradient filtrat [dc/dx dc/dy]
    G.sense = Inf;          % passos seguits sense detectar
    G.dir   = randi(4);     % direccio actual (fila de g.movs)
end

N = p.maxPassos;
traj    = nan(N + 1, 2);   traj(1, :) = [g.x(pos(2)) g.y(pos(1))];
comptes = zeros(N, 1);
trobada = false;

% --- Figura ---
if mostrar
    fig = figure('Color', 'w', 'Name', 'Dron: cerca de l''origen d''una fuga', ...
                 'Units', 'normalized', 'Position', [0.05 0.15 0.9 0.65]);
    colormap(flipud(hot));
    subplot('Position', [0.05 0.12 0.42 0.76]);
    hPloma = imagesc(g.x, g.y, log10(concentracio(pl, p, g.X, g.Y) + 1e-3));
    set(gca, 'YDir', 'normal');
    axis equal;  axis([0 p.Lx 0 p.Ly]);  caxis([-2 2.5]);
    hold on;
    plot(p.font(1), p.font(2), 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    hTraj = plot(traj(1, 1), traj(1, 2), 'b-', 'LineWidth', 1.5);
    hDron = plot(traj(1, 1), traj(1, 2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');
    xlabel('x (m)');  ylabel('y (m)');
    hTit = title('');
    subplot('Position', [0.55 0.12 0.42 0.76]);
    if infotaxis
        hP = imagesc(g.x, g.y, P / max(P(:)));
        set(gca, 'YDir', 'normal');
        axis equal;  axis([0 p.Lx 0 p.Ly]);  caxis([0 1]);
        hold on;
        hDron2 = plot(traj(1, 1), traj(1, 2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');
        xlabel('x (m)');
        title('El que creu el dron: probabilitat de on es la fuga');
    else
        hK = plot(0, 0, 'r-');
        xlabel('Temps (s)');  grid on;
        title('Comptes dels 4 sensors');
    end
end

%% ============================== BUCLE ==============================
for n = 1:N
    xy = traj(n, :);                         % posicio del dron (m)

    % 1. El dron mesura i decideix on va
    if infotaxis
        k = poisson(dt * (p.sens * concentracio(pl, p, xy(1), xy(2)) + p.fons));
        P = infotaxis_actualitzar(P, K, pos, k, g);
        [mov, estat] = infotaxis_decidir(P, K, pos, g, p.bins);
    else
        % 4 sensors en creu: est, oest, nord, sud
        k = poisson(dt * (p.sens * concentracio(pl, p, xy(1) + [p.d -p.d 0 0], xy(2) + [0 0 p.d -p.d]) + p.fons));
        [mov, estat, G] = gradient_decidir(G, k, pos, g, p, dt);
    end
    comptes(n) = sum(k);

    % 2. El dron es mou a una cel.la veina (o es queda quiet)
    pos = pos + mov;
    traj(n + 1, :) = [g.x(pos(2)) g.y(pos(1))];
    trobada = all(abs(pos - celFont) <= 1);  % la fuga es a tocar: trobada

    % 3. La ploma evoluciona (la gaussiana no canvia)
    if complexa
        pl = avancar_ploma(pl, p, dt);
    end

    % 4. Dibuix
    if mostrar && (trobada || mod(n, cadaN) == 0)
        if ~ishandle(fig)
            break;
        end
        if complexa
            set(hPloma, 'CData', log10(concentracio(pl, p, g.X, g.Y) + 1e-3));
        end
        set(hTraj, 'XData', traj(1:n + 1, 1), 'YData', traj(1:n + 1, 2));
        set(hDron, 'XData', traj(n + 1, 1), 'YData', traj(n + 1, 2));
        if infotaxis
            set(hP, 'CData', P / max(P(:)));
            set(hDron2, 'XData', traj(n + 1, 1), 'YData', traj(n + 1, 2));
        else
            set(hK, 'XData', (1:n) * dt, 'YData', comptes(1:n));
        end
        set(hTit, 'String', sprintf('Ploma %s | t = %.0f s | %s', tipusPloma, n * dt, estat));
        drawnow;
    end
    if trobada
        break;
    end
end

%% ============================= RESULTAT =============================
recorregut = sum(sqrt(sum(diff(traj(1:n + 1, :)).^2, 2)));
if trobada
    txt = sprintf('FONT TROBADA en %.0f s (%.0f m recorreguts)', n * dt, recorregut);
else
    txt = sprintf('Font NO trobada en %.0f s (distancia final %.0f m)', n * dt, norm(traj(n + 1, :) - p.font));
end
fprintf('Ploma %s, algorisme %s: %s\n', tipusPloma, algorisme, txt);
if mostrar && ishandle(fig)
    set(hTit, 'String', txt);
end
end

%% ============================== PLOMA ==============================

function c = concentracio(pl, p, x, y)
%CONCENTRACIO  Concentracio real de la ploma als punts (x, y).
if strcmp(pl.tipus, 'gaussiana')
    c = ploma_gaussiana(x - p.font(1), y - p.font(2), p);
else
    % Suma de tots els puffs, cadascun amb forma de campana de Gauss
    d2 = (x(:) - pl.x').^2 + (y(:) - pl.y').^2;             % (punts x puffs)
    c  = reshape(exp(-d2 ./ (2 * pl.s2')) * (pl.m ./ (2 * pi * pl.s2)), size(x));
end
end

function c = ploma_gaussiana(s, n, p)
%PLOMA_GAUSSIANA  Ploma gaussiana estacionaria a una distancia s a favor
%   del vent i n de traves respecte de la fuga:
%
%       c = Q / (2 pi U sy sz) * exp(-n^2 / (2 sy^2)) * [terme vertical]
%
%   sy i sz (amplada i gruix de la ploma) creixen amb s segons les formules
%   de Briggs (1973) per a la classe d'estabilitat D. El terme vertical te
%   en compte les alcades de la fuga i del dron i la reflexio al terra.
c = zeros(size(s));
m = s > -p.sigma0g;                          % a sobrevent no hi arriba res
sm = max(s(m), 0);
sy = 0.08 * sm ./ sqrt(1 + 1e-4 * sm);
sz = 0.06 * sm ./ sqrt(1 + 1.5e-3 * sm);
sy = sqrt(sy.^2 + p.sigma0g^2);              % la fuga te una mida inicial
sz = sqrt(sz.^2 + p.sigma0g^2);
vertical = exp(-(p.z - p.H)^2 ./ (2 * sz.^2)) + exp(-(p.z + p.H)^2 ./ (2 * sz.^2));
c(m) = p.Qg ./ (2 * pi * p.Ug * sy .* sz) .* exp(-n(m).^2 ./ (2 * sy.^2)) .* vertical;
end

function pl = avancar_ploma(pl, p, dt)
%AVANCAR_PLOMA  Un pas de temps de la ploma complexa.
% La direccio del vent oscil.la a l'atzar (proces d'Ornstein-Uhlenbeck)
a = exp(-dt / p.Tmeandre);
pl.th = pl.th * a + (p.sigmaMeandre * pi / 180) * sqrt(1 - a^2) * randn;
% Els puffs es mouen amb el vent i a l'atzar, s'eixamplen i es van esvaint
sd = sqrt(2 * p.Dturb * dt);
pl.x  = pl.x + p.U * cos(pl.th) * dt + sd * randn(size(pl.x));
pl.y  = pl.y + p.U * sin(pl.th) * dt + sd * randn(size(pl.y));
pl.s2 = pl.s2 + 2 * p.Dpuff * dt;
pl.m  = pl.m * exp(-dt / p.tau);
% La fuga emet puffs nous
nNous = round(p.puffsPerSeg * dt);
pl.x  = [pl.x;  p.font(1) * ones(nNous, 1)];
pl.y  = [pl.y;  p.font(2) * ones(nNous, 1)];
pl.s2 = [pl.s2; p.sigma0^2 * ones(nNous, 1)];
pl.m  = [pl.m;  p.Q * dt / nNous * ones(nNous, 1)];
% Eliminem els puffs que han sortit de la zona
fora = pl.x > p.Lx + 60 | pl.x < -60 | pl.y > p.Ly + 60 | pl.y < -60;
pl.x(fora) = [];  pl.y(fora) = [];  pl.s2(fora) = [];  pl.m(fora) = [];
end

function k = poisson(mu)
%POISSON  Nombres aleatoris de Poisson de mitjana mu (sense toolboxes).
k = zeros(size(mu));
for i = 1:numel(mu)
    if mu(i) < 30
        % Metode de Knuth: multiplica uniformes fins a baixar de exp(-mu)
        u = rand;
        while u > exp(-mu(i))
            k(i) = k(i) + 1;
            u = u * rand;
        end
    else
        k(i) = max(0, round(mu(i) + sqrt(mu(i)) * randn));   % aproximacio normal
    end
end
end

%% ============================ INFOTAXIS ============================

function K = model_infotaxis(tipusPloma, p, g, dt)
%MODEL_INFOTAXIS  Model intern del dron: comptes esperats en un pas quan el
%   dron es a un desplacament (DX, DY) de la fuga. K te mida
%   (2Ny-1) x (2Nx-1) i el desplacament zero es a l'element (Ny, Nx).
%
%   - Ploma gaussiana: el dron en coneix la formula (model exacte).
%   - Ploma complexa: solucio mitjana de l'equacio d'adveccio-difusio de
%     l'article, c = Q/(2 pi D) exp(s U/(2D)) K0(r/lambda), amb una
%     difusivitat efectiva D. El dron NO coneix la turbulencia real.
[DX, DY] = meshgrid((-(g.Nx - 1):(g.Nx - 1)) * p.dx, (-(g.Ny - 1):(g.Ny - 1)) * p.dx);
if strcmp(tipusPloma, 'gaussiana')
    c = ploma_gaussiana(DX, DY, p);
else
    D   = p.Dturb + p.Dpuff + (p.U * p.sigmaMeandre * pi / 180)^2 * p.Tmeandre;
    lam = sqrt(D * p.tau / (1 + p.U^2 * p.tau / (4 * D)));
    r   = max(hypot(DX, DY), p.dx / 2);
    % besselk(0, x, 1) = K0(x) exp(x): evita desbordaments numerics
    c   = p.Q / (2 * pi * D) * exp(DX * p.U / (2 * D) - r / lam) .* besselk(0, r / lam, 1);
end
K = dt * (p.sens * c + p.fons);
end

function MU = mu_local(K, q, g)
%MU_LOCAL  Comptes esperats a la cel.la q per a cada possible posicio de
%   la fuga a la graella (matriu Ny x Nx).
MU = K(q(1) + g.Ny - 1:-1:q(1), q(2) + g.Nx - 1:-1:q(2));
end

function m = zona(q, g)
%ZONA  Cel.les a tocar de q (bloc de 3 x 3): si la fuga hi es, el dron la veu.
m = false(g.Ny, g.Nx);
m(max(q(1) - 1, 1):min(q(1) + 1, g.Ny), max(q(2) - 1, 1):min(q(2) + 1, g.Nx)) = true;
end

function S = entropia(P)
%ENTROPIA  Entropia de Shannon S = -sum P ln P: mesura la incertesa.
P = P(P > 0);
S = -sum(P .* log(P));
end

function P = infotaxis_actualitzar(P, K, pos, k, g)
%INFOTAXIS_ACTUALITZAR  Regla de Bayes (eq. 2 de l'article) despres de
%   mesurar k comptes: P(r0) <- P(r0) * P(k | r0), amb P(k | r0) de Poisson.
%   La zona que el dron veu queda descartada (la fuga no hi es).
MU = mu_local(K, pos, g);
Pn = P .* exp(-MU + k * log(MU) - gammaln(k + 1));
Pn(zona(pos, g)) = 0;
if sum(Pn(:)) > 0
    P = Pn / sum(Pn(:));
else
    P(zona(pos, g)) = 0;             % mesura incompatible amb el model: nomes descartem la zona
    P = P / sum(P(:));
end
end

function [mov, estat] = infotaxis_decidir(P, K, pos, g, bins)
%INFOTAXIS_DECIDIR  Tria el moviment que fa baixar mes l'entropia esperada
%   (eq. 1 de l'article):
%
%       dS(r -> rj) = P(rj) [-S] + [1 - P(rj)] sum_k rho_k(rj) dS_k
%
%   Primer terme: probabilitat de trobar la fuga a rj (l'entropia passa a 0).
%   Segon terme: si no hi es, el dron rebra k comptes amb probabilitat rho_k
%   i l'entropia canviara dS_k. Els comptes s'agrupen en intervals (bins).
S0 = entropia(P);
dS = Inf(5, 1);
for j = 1:5
    q = pos + g.movs(j, :);
    if q(1) < 1 || q(2) < 1 || q(1) > g.Ny || q(2) > g.Nx
        continue;                    % fora de la zona
    end
    z  = zona(q, g);
    Pj = sum(P(z));                  % probabilitat de trobar la fuga en arribar a q
    Pp = P;
    Pp(z) = 0;                       % distribucio si NO la troba
    if sum(Pp(:)) <= 0
        dS(j) = -S0;
        continue;
    end
    Pp = Pp / sum(Pp(:));
    MU = mu_local(K, q, g);

    % Recorrem k = 0, 1, 2... acumulant la llei de Poisson per intervals
    canvi = 0;
    acc   = zeros(size(MU));         % P(k < limit actual)
    term  = exp(-MU);                % P(k) per al k actual
    k = 0;
    for ib = 1:numel(bins)
        if ib < numel(bins)
            L = zeros(size(MU));
            while k < bins(ib + 1)
                L    = L + term;
                term = term .* MU / (k + 1);
                k    = k + 1;
            end
            acc = acc + L;
        else
            L = max(1 - acc, 0);     % ultim interval: k >= bins(end)
        end
        w   = Pp .* L;
        rho = sum(w(:));             % probabilitat de mesurar comptes d'aquest interval
        if rho > 1e-300
            canvi = canvi + rho * (entropia(w / rho) - S0);
        end
    end
    dS(j) = Pj * (-S0) + (1 - Pj) * canvi;
end

if min(dS) > -1e-6
    % Cap moviment aporta informacio (el dron ja "sap" on es la fuga):
    % vola directament cap al punt mes probable
    estat = 'cap al punt mes probable';
    [~, im] = max(P(:));
    [iy, ix] = ind2sub(size(P), im);
    d = [iy ix] - pos;
    if abs(d(2)) >= abs(d(1))
        mov = [0 sign(d(2))];
    else
        mov = [sign(d(1)) 0];
    end
else
    estat = 'infotaxis';
    [~, j] = min(dS + 1e-10 * abs(S0) * rand(5, 1));    % desempat aleatori
    mov = g.movs(j, :);
end
end

%% ============================= GRADIENT =============================

function [mov, estat, G] = gradient_decidir(G, k, pos, g, p, dt)
%GRADIENT_DECIDIR  Estima el gradient amb els comptes k dels 4 sensors
%   (est, oest, nord, sud) per diferencies finites centrades:
%
%       dc/dx ~ (cE - cO) / (2d),    dc/dy ~ (cN - cS) / (2d)
%
%   i el suavitza amb un filtre exponencial. Si hi ha senyal recent, el
%   dron va cap on el gradient es mes gran; si no, explora a l'atzar.
c = max(k / dt - p.fons, 0) / p.sens;            % concentracio estimada a cada sensor
gNou = [c(1) - c(2), c(3) - c(4)] / (2 * p.d);
if sum(k) >= p.llindar
    G.sense = 0;
    G.gf = p.alfa * G.gf + (1 - p.alfa) * gNou;
else
    G.sense = G.sense + 1;
    if G.sense > p.memoria
        G.gf = [0 0];                            % senyal perdut: oblida el gradient
    end
end

q = pos + g.movs(1:4, :);                        % les 4 cel.les veines
valid = q(:, 1) >= 1 & q(:, 1) <= g.Ny & q(:, 2) >= 1 & q(:, 2) <= g.Nx;
if G.sense <= p.memoria
    estat = 'seguint el gradient';
    proj = g.movs(1:4, [2 1]) * G.gf' + 1e-12 * rand(4, 1);   % projeccio del gradient sobre cada direccio
    proj(~valid) = -Inf;
    [~, G.dir] = max(proj);
else
    estat = 'explorant';
    if rand < p.pGir
        G.dir = mod(G.dir - 1 + 2 * randi(2) - 3, 4) + 1;     % gira 90 graus a esquerra o dreta
    end
    if ~valid(G.dir)
        G.dir = mod(G.dir + 1, 4) + 1;                        % sortiria de la zona: mitja volta
    end
end
mov = g.movs(G.dir, :);
end
