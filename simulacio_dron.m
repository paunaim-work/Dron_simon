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
%                     'gradient'  (puja pel gradient de concentracio)
%
%   Executa-la amb F5. Les opcions i els parametres son just a sota.

%% ============================= OPCIONS =============================
dimensio   = '3d';          % '2d' o '3d'
tipusPloma = 'complexa';    % 'gaussiana' o 'complexa'
algorisme  = 'infotaxis';   % 'infotaxis' o 'gradient'
llavor     = 1;             % mateix numero = mateixa simulacio ([] = cada cop diferent)
mostrar    = true;          % false = sense animacio (mes rapid)
cadaN      = 2;             % dibuixa cada N passos

%% ============================ PARAMETRES ============================
% Unitats: metres i segons. El vent bufa cap a +x i el terra es z = 0.
p.Lx = 400;  p.Ly = 300;    % mida horitzontal de la zona de cerca
p.Lz = 40;                  % alcada maxima de vol (nomes en 3D; la minima es dx)
p.dx = 5;                   % mida de cel.la = el que avanca el dron a cada pas
p.font  = [60 150 10];      % on es la fuga [x y z] (el dron NO ho sap)
p.inici = [340 215 30];     % d'on surt el dron [x y z] (en 2D, z = alcada de la fuga)
p.vDron = 5;                % velocitat del dron (m/s)
p.maxPassos = 2500;         % limit de passos (autonomia)

% Sensor: comptador de deteccions (tipus Geiger).
% Comptes en un pas ~ Poisson( dt * (sens * concentracio + fons) )
p.sens = 1;                 % comptes/s per unitat de concentracio
p.fons = 0.02;              % comptes/s de fons (falses deteccions)

% Ploma gaussiana
p.Qg = 5000;                % taxa d'emissio
p.Ug = 3;                   % velocitat del vent (m/s)
p.sigma0g = 1;              % mida de la fuga (m)

% Ploma complexa
p.Q   = 2000;               % taxa d'emissio
p.U   = 2;                  % velocitat mitjana del vent (m/s)
p.tau = 200;                % temps de vida del contaminant (s)
p.puffsPerSeg  = 2;         % puffs emesos per segon
p.sigma0       = 1.5;       % mida inicial d'un puff (m)
p.DpuffH = 0.4;  p.DpuffV = 0.2;   % rapidesa amb que s'eixampla un puff, horitzontal i vertical (m^2/s)
p.DturbH = 0.4;  p.DturbV = 0.2;   % moviment aleatori dels puffs, horitzontal i vertical (m^2/s)
p.sigmaMeandre = 15;        % oscil.lacio de la direccio del vent (graus)
p.Tmeandre     = 40;        % temps tipic d'aquesta oscil.lacio (s)

% Infotaxis
p.bins = [0 1 2 3 5 8 13 21 34 55 89];   % intervals de comptes per planificar

% Gradient
p.d       = 5;              % distancia dels sensors al centre del dron (m)
p.alfa    = 0.6;            % suavitzat del gradient (0 = cap, 0.9 = molt)
p.llindar = 2;              % comptes (suma dels sensors) per considerar que detecta
p.memoria = 8;              % passos sense detectar abans de tornar a explorar
p.pGir    = 0.15;           % en exploracio, probabilitat de girar a cada pas
p.pVert   = 0.05;           % en exploracio 3D, probabilitat de canviar d'alcada a cada pas
%% ====================================================================

if ~isempty(llavor)
    rng(llavor);
end
infotaxis = strcmp(algorisme, 'infotaxis');
complexa  = strcmp(tipusPloma, 'complexa');
tresD     = strcmp(dimensio, '3d');

% --- Graella: les matrius tenen mida Ny x Nx x Nz i s'indexen (iy, ix, iz) ---
g.x = 0:p.dx:p.Lx;
g.y = 0:p.dx:p.Ly;
if tresD
    g.z = p.dx:p.dx:p.Lz;                % el dron no vola arran de terra
else
    g.z = p.font(3);                     % 2D: una sola alcada, la de la fuga
    p.inici(3) = p.font(3);
end
g.Nx = numel(g.x);  g.Ny = numel(g.y);  g.Nz = numel(g.z);
[g.X, g.Y] = meshgrid(g.x, g.y);
% Moviments [diy dix diz]: est, nord, oest, sud, amunt, avall, quiet
g.movs = [0 1 0; 1 0 0; 0 -1 0; -1 0 0; 0 0 1; 0 0 -1; 0 0 0];
dt = p.dx / p.vDron;                     % durada d'un pas (s)

celFont = cella(p.font, p, g);           % cel.la [iy ix iz] de la fuga
pos     = cella(p.inici, p, g);          % cel.la [iy ix iz] del dron

% --- Ploma: la complexa evoluciona 300 s abans que surti el dron ---
buit = zeros(0, 1);
pl = struct('tipus', tipusPloma, 'x', buit, 'y', buit, 'z', buit, ...
            's2h', buit, 's2v', buit, 'm', buit, 'th', 0);
if complexa
    for i = 1:round(300 / dt)
        pl = avancar_ploma(pl, p, dt);
    end
end

% --- Estat inicial del cercador ---
if infotaxis
    K = model_infotaxis(tipusPloma, p, g, dt);   % comptes esperats segons el model del dron
    P = ones(g.Ny, g.Nx, g.Nz);                  % probabilitat de on es la fuga: uniforme
    P(zona(pos, g)) = 0;
    P = P / sum(P(:));
else
    % Sensors en creu: est, oest, nord, sud i, en 3D, un a sobre i un a sota
    sensors = p.d * [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
    sensors = sensors(1:4 + 2 * tresD, :);
    G.gf    = [0 0 0];      % gradient filtrat [dc/dx dc/dy dc/dz]
    G.sense = Inf;          % passos seguits sense detectar
    G.dir   = randi(4);     % direccio horitzontal actual (fila de g.movs)
end

N = p.maxPassos;
traj    = nan(N + 1, 3);   traj(1, :) = [g.x(pos(2)) g.y(pos(1)) g.z(pos(3))];
comptes = zeros(N, 1);
trobada = false;

% --- Figura ---
if mostrar
    r = traj(1, :);
    fig = figure('Color', 'w', 'Name', 'Dron: cerca de l''origen d''una fuga', ...
                 'Units', 'normalized', 'Position', [0.05 0.15 0.9 0.65]);
    colormap(flipud(hot));
    % Esquerra: la realitat. Al terra es pinta la ploma a l'alcada del dron.
    subplot('Position', [0.05 0.12 0.42 0.76]);
    hPloma = surf(g.x, g.y, zeros(g.Ny, g.Nx), log10(concentracio(pl, p, g.X, g.Y, r(3)) + 1e-3), ...
                  'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;  caxis([-2 2.5]);
    plot3(p.font(1) * [1 1], p.font(2) * [1 1], [0 p.font(3)], 'g-', 'LineWidth', 1.5);
    plot3(p.font(1), p.font(2), p.font(3), 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    hTraj = plot3(r(1), r(2), r(3), 'b-', 'LineWidth', 1.5);
    hPal  = plot3(r(1) * [1 1], r(2) * [1 1], [0 r(3)], 'k-');
    hDron = plot3(r(1), r(2), r(3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');
    axis([0 p.Lx 0 p.Ly 0 p.Lz]);
    xlabel('x (m)');  ylabel('y (m)');  zlabel('z (m)');
    if tresD
        daspect([1 1 0.25]);  view(-30, 35);  grid on;
    else
        daspect([1 1 1]);  view(2);  grid off;
    end
    hTit = title('');
    % Dreta: el que sap el dron
    subplot('Position', [0.55 0.12 0.42 0.76]);
    if infotaxis
        hP = imagesc(g.x, g.y, sum(P, 3) / max(max(sum(P, 3))));
        set(gca, 'YDir', 'normal');
        axis equal;  axis([0 p.Lx 0 p.Ly]);  caxis([0 1]);
        hold on;
        hDron2 = plot(r(1), r(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');
        xlabel('x (m)');  ylabel('y (m)');
        title('El que creu el dron: probabilitat de on es la fuga (vista des de dalt)');
    else
        hK = plot(0, 0, 'r-');
        xlabel('Temps (s)');  grid on;
        title('Comptes dels sensors');
    end
end

%% ============================== BUCLE ==============================
for n = 1:N
    r = traj(n, :);                          % posicio del dron [x y z] (m)

    % 1. El dron mesura i decideix on va
    if infotaxis
        k = poisson(dt * (p.sens * concentracio(pl, p, r(1), r(2), r(3)) + p.fons));
        P = infotaxis_actualitzar(P, K, pos, k, g);
        [mov, estat] = infotaxis_decidir(P, K, pos, g, p.bins);
    else
        s = r + sensors;
        k = poisson(dt * (p.sens * concentracio(pl, p, s(:, 1), s(:, 2), s(:, 3)) + p.fons));
        [mov, estat, G] = gradient_decidir(G, k, pos, g, p, dt);
    end
    comptes(n) = sum(k);

    % 2. El dron es mou a una cel.la veina (o es queda quiet)
    pos = pos + mov;
    r = [g.x(pos(2)) g.y(pos(1)) g.z(pos(3))];
    traj(n + 1, :) = r;
    % La fuga es troba si es a tocar (bloc de 3 x 3 cel.les) i a la mateixa alcada
    trobada = all(abs(pos(1:2) - celFont(1:2)) <= 1) && pos(3) == celFont(3);

    % 3. La ploma evoluciona (la gaussiana no canvia)
    if complexa
        pl = avancar_ploma(pl, p, dt);
    end

    % 4. Dibuix
    if mostrar && (trobada || mod(n, cadaN) == 0)
        if ~ishandle(fig)
            break;
        end
        set(hPloma, 'CData', log10(concentracio(pl, p, g.X, g.Y, r(3)) + 1e-3));
        set(hTraj, 'XData', traj(1:n + 1, 1), 'YData', traj(1:n + 1, 2), 'ZData', traj(1:n + 1, 3));
        set(hPal,  'XData', r(1) * [1 1], 'YData', r(2) * [1 1], 'ZData', [0 r(3)]);
        set(hDron, 'XData', r(1), 'YData', r(2), 'ZData', r(3));
        if infotaxis
            set(hP, 'CData', sum(P, 3) / max(max(sum(P, 3))));
            set(hDron2, 'XData', r(1), 'YData', r(2));
        else
            set(hK, 'XData', (1:n) * dt, 'YData', comptes(1:n));
        end
        set(hTit, 'String', {sprintf('Ploma %s %s | t = %.0f s | alcada %.0f m | %s', ...
                                     tipusPloma, dimensio, n * dt, r(3), estat), ...
                             'Colors: concentracio a l''alcada del dron'});
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
fprintf('%s, ploma %s, algorisme %s: %s\n', upper(dimensio), tipusPloma, algorisme, txt);
if mostrar && ishandle(fig)
    set(hTit, 'String', txt);
end
end

function c = cella(r, p, g)
%CELLA  Cel.la [iy ix iz] de la graella que correspon a la posicio r = [x y z].
[~, iz] = min(abs(g.z - r(3)));
c = [round(r(2) / p.dx) + 1, round(r(1) / p.dx) + 1, iz];
end

%% ============================== PLOMA ==============================

function c = concentracio(pl, p, x, y, z)
%CONCENTRACIO  Concentracio real de la ploma als punts (x, y, z).
if strcmp(pl.tipus, 'gaussiana')
    c = ploma_gaussiana(x - p.font(1), y - p.font(2), z, p.font(3), p);
else
    % Suma de tots els puffs. Cada puff es una campana de Gauss en 3D, mes
    % la seva imatge sota terra (el terra reflecteix el contaminant)
    w  = pl.m ./ ((2 * pi)^1.5 * pl.s2h .* sqrt(pl.s2v));
    dh = (x(:) - pl.x').^2 + (y(:) - pl.y').^2;             % (punts x puffs)
    ez = exp(-(z(:) - pl.z').^2 ./ (2 * pl.s2v')) + exp(-(z(:) + pl.z').^2 ./ (2 * pl.s2v'));
    c  = reshape((exp(-dh ./ (2 * pl.s2h')) .* ez) * w, size(x));
end
end

function c = ploma_gaussiana(s, n, z, H, p)
%PLOMA_GAUSSIANA  Ploma gaussiana estacionaria a una distancia s a favor
%   del vent, n de traves i a l'alcada z, si la fuga es a l'alcada H:
%
%                Q              n^2             (z-H)^2            (z+H)^2
%       c = ------------ exp(- ------) [ exp(- -------) + exp(- -------) ]
%           2 pi U sy sz       2 sy^2           2 sz^2             2 sz^2
%
%   sy i sz (amplada i gruix de la ploma) creixen amb s segons les formules
%   de Briggs (1973) per a la classe d'estabilitat D. L'ultim terme es la
%   reflexio del contaminant al terra.
z = z + zeros(size(s));
c = zeros(size(s));
m = s > -p.sigma0g;                          % a sobrevent no hi arriba res
sm = max(s(m), 0);
sy = 0.08 * sm ./ sqrt(1 + 1e-4 * sm);
sz = 0.06 * sm ./ sqrt(1 + 1.5e-3 * sm);
sy = sqrt(sy.^2 + p.sigma0g^2);              % la fuga te una mida inicial
sz = sqrt(sz.^2 + p.sigma0g^2);
vertical = exp(-(z(m) - H).^2 ./ (2 * sz.^2)) + exp(-(z(m) + H).^2 ./ (2 * sz.^2));
c(m) = p.Qg ./ (2 * pi * p.Ug * sy .* sz) .* exp(-n(m).^2 ./ (2 * sy.^2)) .* vertical;
end

function pl = avancar_ploma(pl, p, dt)
%AVANCAR_PLOMA  Un pas de temps de la ploma complexa.
% La direccio del vent oscil.la a l'atzar (proces d'Ornstein-Uhlenbeck)
a = exp(-dt / p.Tmeandre);
pl.th = pl.th * a + (p.sigmaMeandre * pi / 180) * sqrt(1 - a^2) * randn;
% Els puffs es mouen amb el vent i a l'atzar, s'eixamplen i es van esvaint
sh = sqrt(2 * p.DturbH * dt);
sv = sqrt(2 * p.DturbV * dt);
pl.x   = pl.x + p.U * cos(pl.th) * dt + sh * randn(size(pl.x));
pl.y   = pl.y + p.U * sin(pl.th) * dt + sh * randn(size(pl.y));
pl.z   = abs(pl.z + sv * randn(size(pl.z)));     % rebot al terra
pl.s2h = pl.s2h + 2 * p.DpuffH * dt;
pl.s2v = pl.s2v + 2 * p.DpuffV * dt;
pl.m   = pl.m * exp(-dt / p.tau);
% La fuga emet puffs nous
nous = ones(round(p.puffsPerSeg * dt), 1);
pl.x   = [pl.x;   p.font(1) * nous];
pl.y   = [pl.y;   p.font(2) * nous];
pl.z   = [pl.z;   p.font(3) * nous];
pl.s2h = [pl.s2h; p.sigma0^2 * nous];
pl.s2v = [pl.s2v; p.sigma0^2 * nous];
pl.m   = [pl.m;   p.Q * dt / numel(nous) * nous];
% Eliminem els puffs que han sortit de la zona
fora = pl.x > p.Lx + 60 | pl.x < -60 | pl.y > p.Ly + 60 | pl.y < -60;
pl.x(fora) = [];    pl.y(fora) = [];    pl.z(fora) = [];
pl.s2h(fora) = [];  pl.s2v(fora) = [];  pl.m(fora) = [];
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
%   dron es a l'alcada z(iz), la fuga a l'alcada z(i0) i el desplacament
%   horitzontal entre tots dos es (DX, DY). Es una taula
%   K(iy - iy0 + Ny, ix - ix0 + Nx, iz, i0).
%
%   - Ploma gaussiana: el dron en coneix la formula (model exacte).
%   - Ploma complexa: solucio mitjana de l'equacio d'adveccio-difusio de
%     l'article, amb difusivitats efectives Dh (horitzontal) i Dv
%     (vertical) i una font imatge sota terra:
%
%       c = Q / (4 pi sqrt(Dh Dv)) exp(s U / 2Dh) [e^(-r1/lam)/r1 + e^(-r2/lam)/r2]
%
%     r1 i r2 son les distancies a la fuga i a la seva imatge (amb z
%     escalada per sqrt(Dh/Dv)). El dron NO coneix la turbulencia real.
[DX, DY] = meshgrid((-(g.Nx - 1):(g.Nx - 1)) * p.dx, (-(g.Ny - 1):(g.Ny - 1)) * p.dx);
Dh  = p.DturbH + p.DpuffH + (p.U * p.sigmaMeandre * pi / 180)^2 * p.Tmeandre;
Dv  = p.DturbV + p.DpuffV;
lam = sqrt(Dh * p.tau / (1 + p.U^2 * p.tau / (4 * Dh)));
K = zeros(2 * g.Ny - 1, 2 * g.Nx - 1, g.Nz, g.Nz);
for iz = 1:g.Nz
    for i0 = 1:g.Nz
        if strcmp(tipusPloma, 'gaussiana')
            c = ploma_gaussiana(DX, DY, g.z(iz), g.z(i0), p);
        else
            r1 = max(sqrt(DX.^2 + DY.^2 + Dh / Dv * (g.z(iz) - g.z(i0))^2), p.dx / 2);
            r2 = max(sqrt(DX.^2 + DY.^2 + Dh / Dv * (g.z(iz) + g.z(i0))^2), p.dx / 2);
            c  = p.Q / (4 * pi * sqrt(Dh * Dv)) * ...
                 (exp(DX * p.U / (2 * Dh) - r1 / lam) ./ r1 + exp(DX * p.U / (2 * Dh) - r2 / lam) ./ r2);
        end
        K(:, :, iz, i0) = dt * (p.sens * c + p.fons);
    end
end
end

function MU = mu_local(K, q, g)
%MU_LOCAL  Comptes esperats a la cel.la q per a cada possible posicio de
%   la fuga a la graella (matriu Ny x Nx x Nz).
MU = reshape(K(q(1) + g.Ny - 1:-1:q(1), q(2) + g.Nx - 1:-1:q(2), q(3), :), g.Ny, g.Nx, g.Nz);
end

function m = zona(q, g)
%ZONA  Cel.les a tocar de q (bloc de 3 x 3 a la mateixa alcada): si la
%   fuga hi es, el dron la veu.
m = false(g.Ny, g.Nx, g.Nz);
m(max(q(1) - 1, 1):min(q(1) + 1, g.Ny), max(q(2) - 1, 1):min(q(2) + 1, g.Nx), q(3)) = true;
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
dS = Inf(7, 1);
for j = 1:7
    q = pos + g.movs(j, :);
    if any(q < 1) || any(q > [g.Ny g.Nx g.Nz])
        continue;                    % fora de la zona (en 2D, amunt i avall)
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
    viu = Pp > 1e-15;                % nomes compten les posicions encara possibles (mes rapid)
    Pp = Pp(viu);
    MU = MU(viu);

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
    [iy, ix, iz] = ind2sub([g.Ny g.Nx g.Nz], im);
    d = [iy ix iz] - pos;
    [~, eix] = max(abs(d));
    mov = [0 0 0];
    mov(eix) = sign(d(eix));
else
    estat = 'infotaxis';
    [~, j] = min(dS + 1e-10 * abs(S0) * rand(7, 1));    % desempat aleatori
    mov = g.movs(j, :);
end
end

%% ============================= GRADIENT =============================

function [mov, estat, G] = gradient_decidir(G, k, pos, g, p, dt)
%GRADIENT_DECIDIR  Estima el gradient amb els comptes k dels sensors (est,
%   oest, nord, sud i, en 3D, amunt i avall) per diferencies finites
%   centrades:
%
%       dc/dx ~ (cE - cO) / (2d),   dc/dy ~ (cN - cS) / (2d),   dc/dz ~ (cAmunt - cAvall) / (2d)
%
%   i el suavitza amb un filtre exponencial. Si hi ha senyal recent, el
%   dron va cap on el gradient es mes gran; si no, explora a l'atzar.
c = max(k(:)' / dt - p.fons, 0) / p.sens;        % concentracio estimada a cada sensor
gNou = [0 0 0];
gNou(1:numel(c) / 2) = (c(1:2:end) - c(2:2:end)) / (2 * p.d);
if sum(k) >= p.llindar
    G.sense = 0;
    G.gf = p.alfa * G.gf + (1 - p.alfa) * gNou;
else
    G.sense = G.sense + 1;
    if G.sense > p.memoria
        G.gf = [0 0 0];                          % senyal perdut: oblida el gradient
    end
end

q = pos + g.movs(1:6, :);                        % les 6 cel.les veines
valid = all(q >= 1, 2) & q(:, 1) <= g.Ny & q(:, 2) <= g.Nx & q(:, 3) <= g.Nz;
if G.sense <= p.memoria
    estat = 'seguint el gradient';
    proj = g.movs(1:6, [2 1 3]) * G.gf' + 1e-12 * rand(6, 1);   % projeccio del gradient sobre cada direccio
    proj(~valid) = -Inf;
    [~, j] = max(proj);
    if j <= 4
        G.dir = j;
    end
else
    estat = 'explorant';
    if rand < p.pGir
        G.dir = mod(G.dir - 1 + 2 * randi(2) - 3, 4) + 1;     % gira 90 graus a esquerra o dreta
    end
    if ~valid(G.dir)
        G.dir = mod(G.dir + 1, 4) + 1;                        % sortiria de la zona: mitja volta
    end
    j = G.dir;
    vert = 4 + find(valid(5:6))';                             % moviments verticals possibles (cap en 2D)
    if ~isempty(vert) && rand < p.pVert
        j = vert(randi(numel(vert)));                         % de tant en tant canvia d'alcada
    end
end
mov = g.movs(j, :);
end
