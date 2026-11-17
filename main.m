function main()
%MAIN  Un dron autonom busca l'origen d'una fuga de contaminant.
%
%   Tria les opcions de sota i executa amb F5. El codi esta en 3 fitxers:
%
%     main.m   Aquest: opcions, escenari, simulacio i dibuix.
%     Ploma.m  Genera la ploma (4 versions) i la guarda en un fitxer .nc.
%     Dron.m   Classe Dron: parametres, sensor i algorismes de cerca.
%
%   Com funciona:  Ploma.generar -> fitxer .nc -> Ploma.llegir -> el dron vola

%% ============================== OPCIONS ==============================
opcions.dimensio     = '3d';         % '2d' o '3d'
opcions.tipusPloma   = 'complexa';   % 'gaussiana' o 'complexa'
opcions.algorisme    = 'infotaxis';  % 'infotaxis' o 'gradient'
opcions.generarPloma = true;         % false = fa servir el .nc que ja hi ha (no el torna a generar)
opcions.llavor       = 1;            % mateix numero = mateixa simulacio ([] = cada cop diferent)
opcions.animacio     = true;         % false = sense dibuix (molt mes rapid)
opcions.dibuixaCada  = 2;            % redibuixa cada tants segons de simulacio

%% ============================= ESCENARI ==============================
% Unitats: metres i segons. El vent bufa cap a +x i el terra es z = 0.
escenari.midaX          = 400;           % llargada de la zona (direccio del vent)
escenari.midaY          = 300;           % amplada de la zona (de traves al vent)
escenari.midaZ          = 50;            % alcada de la zona (nomes en 3D)
escenari.posFuga        = [60 150 10];   % on es la fuga [x y z] (el dron NO ho sap)
escenari.posInicialDron = [340 215 30];  % d'on surt el dron [x y z]

%% ============================ PARAMETRES =============================
% Els valors per defecte son a Ploma.m i Dron.m. Es poden canviar aqui, per exemple:
%   parDron.velocitat = 8;      parPloma.complexa.vent = 3;
parPloma = Ploma.parametres();
parDron  = Dron.parametres();
% El que el dron sap de la ploma (per a l'infotaxis): li diem el tipus i els seus parametres
parDron.model = parPloma.(opcions.tipusPloma);
parDron.model.tipus = opcions.tipusPloma;

%% ============================ SIMULACIO ==============================
if ~isempty(opcions.llavor)
    rng(opcions.llavor);
end
fitxer = sprintf('ploma_%s_%s.nc', opcions.tipusPloma, opcions.dimensio);
if opcions.generarPloma
    camp = Ploma.generar(opcions.tipusPloma, opcions.dimensio, escenari, parPloma);
    Ploma.guardar(camp, fitxer);
end
camp = Ploma.llegir(fitxer);                 % el dron treballa amb el que hi ha al fitxer

res = simular(opcions, escenari, camp, parDron);
fprintf('%s, ploma %s, algorisme %s: %s\n', upper(opcions.dimensio), ...
        opcions.tipusPloma, opcions.algorisme, res.text);
end


function res = simular(opcions, escenari, camp, parDron)
%SIMULAR  El dron busca la fuga dins del camp (la ploma llegida del .nc).
%   Retorna una struct amb: trobada (true/false), temps (s), recorregut (m),
%   distanciaFinal (m), trajectoria (posicio [x y z] a cada pas) i text.
if strcmp(opcions.dimensio, '2d')
    escenari.posInicialDron(3) = camp.z(1);  % en 2D vola a l'unica alcada del fitxer
end
dron = Dron(opcions.algorisme, opcions.dimensio, escenari, parDron);

trajectoria  = nan(dron.maxPassos + 1, 3);     % posicio del dron a cada pas
trajectoria(1, :) = dron.pos;
concMesurada = zeros(dron.maxPassos, 1);       % concentracio mesurada a cada pas
if opcions.animacio
    fig = crear_figura(camp, dron);
end

for pas = 1:dron.maxPassos
    dron.pas(camp);                              % el dron mesura, decideix i es mou
    trobada = dron.veuFuga(camp.posFuga);        % ha arribat a la fuga?

    trajectoria(pas + 1, :) = dron.pos;
    concMesurada(pas) = dron.concMesurada;

    if opcions.animacio && (trobada || mod(pas, max(1, round(opcions.dibuixaCada / dron.dt))) == 0)
        if ~ishandle(fig.finestra)
            break;                               % s'ha tancat la finestra
        end
        titol = sprintf('%s | ploma %s | t = %.0f s | alcada %.0f m | %s', upper(opcions.dimensio), ...
                        opcions.tipusPloma, dron.temps, dron.pos(3), dron.queFa);
        actualitzar_figura(fig, camp, dron, trajectoria(1:pas + 1, :), concMesurada(1:pas), titol);
    end
    if trobada
        break;
    end
end

% Resultat
res.trobada        = trobada;
res.temps          = dron.temps;
res.trajectoria    = trajectoria(1:pas + 1, :);
res.recorregut     = sum(sqrt(sum(diff(res.trajectoria).^2, 2)));
res.distanciaFinal = norm(dron.pos - camp.posFuga);
if trobada
    res.text = sprintf('FUGA TROBADA en %.0f s (%.0f m recorreguts)', res.temps, res.recorregut);
else
    res.text = sprintf('Fuga NO trobada en %.0f s (distancia final %.0f m)', res.temps, res.distanciaFinal);
end
if opcions.animacio && ishandle(fig.finestra)
    set(fig.titol, 'String', res.text);
end
end


%% =============================== DIBUIX ===============================

function fig = crear_figura(camp, dron)
%CREAR_FIGURA  Crea la finestra i en guarda els elements que es van actualitzant:
%   - l'escena real: ploma, fuga, dron i trajectoria (en 3D, tambe de costat)
%   - el que sap el dron, vist des de dalt
%   - la concentracio que mesura al llarg del temps
g     = dron.graella;
tresD = numel(camp.z) > 1;
f     = camp.posFuga;
fig.tresD = tresD;
fig.finestra = figure('Color', 'w', 'Name', 'Dron: cerca de l''origen d''una fuga', ...
                      'Units', 'normalized', 'Position', [0.03 0.08 0.94 0.84]);
colormap(flipud(hot));
estilDron = {'o', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k'};
estilFuga = {'p', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k'};

% Escena real. Al terra es pinta la ploma a l'alcada del dron (escala log10).
if tresD
    subplot('Position', [0.05 0.40 0.42 0.52]);
else
    subplot('Position', [0.05 0.12 0.42 0.76]);
end
buit = zeros(numel(camp.y), numel(camp.x));
fig.terra = surf(camp.x, camp.y, buit, buit, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;  caxis([-2 2.5]);
plot3(f(1) * [1 1], f(2) * [1 1], [0 f(3)], 'g-', 'LineWidth', 1.5);
plot3(f(1), f(2), f(3), estilFuga{:});
fig.trajectoria = plot3(nan, nan, nan, 'b-', 'LineWidth', 1.5);
fig.pal         = plot3(nan, nan, nan, 'k-');      % linia vertical del dron fins a terra
fig.dron        = plot3(nan, nan, nan, estilDron{:});
axis([camp.x(1) camp.x(end) camp.y(1) camp.y(end) 0 max(camp.z(end), f(3) + 1)]);
xlabel('x (m)');  ylabel('y (m)');  zlabel('z (m)');
if tresD
    daspect([1 1 0.25]);  view(-30, 35);  grid on;
else
    daspect([1 1 1]);  view(2);
end
fig.titol = title('');

% Vista lateral (nomes 3D): tall vertical de la ploma per la y del dron
if tresD
    subplot('Position', [0.05 0.08 0.42 0.22]);
    fig.lateral = imagesc(camp.x, camp.z, zeros(numel(camp.z), numel(camp.x)));
    set(gca, 'YDir', 'normal');  caxis([-2 2.5]);  hold on;
    plot(f(1), f(3), estilFuga{:});
    fig.trajLateral = plot(nan, nan, 'b-', 'LineWidth', 1.5);
    fig.dronLateral = plot(nan, nan, estilDron{:});
    axis([camp.x(1) camp.x(end) camp.z(1) camp.z(end)]);
    xlabel('x (m)');  ylabel('z (m)');
    title('Vista lateral: ploma en el pla vertical que passa pel dron');
end

% El que sap el dron, vist des de dalt
subplot('Position', [0.56 0.42 0.40 0.50]);
fig.mapa = imagesc(g.x, g.y, zeros(g.mida(1), g.mida(2)));
set(gca, 'YDir', 'normal');  axis equal;  axis([g.x(1) g.x(end) g.y(1) g.y(end)]);  hold on;
fig.dronMapa = plot(nan, nan, estilDron{:});
xlabel('x (m)');  ylabel('y (m)');
if strcmp(dron.algorisme, 'infotaxis')
    caxis([0 1]);
    title('El que creu el dron: probabilitat de on es la fuga');
else
    caxis([-2 2.5]);
    title('El que sap el dron: concentracio que ha mesurat a cada lloc');
end

% Concentracio mesurada al llarg del temps (escala logaritmica)
subplot('Position', [0.56 0.10 0.40 0.22]);
fig.mesura = plot(nan, nan, 'r-');
set(gca, 'YScale', 'log');  ylim([0.05 1000]);  grid on;
xlabel('Temps (s)');
title('Concentracio mesurada pel dron');
end


function actualitzar_figura(fig, camp, dron, trajectoria, concMesurada, titol)
%ACTUALITZAR_FIGURA  Redibuixa tots els panells amb l'estat actual.
p = dron.pos;
[~, it] = min(abs(camp.t - dron.temps));      % instant del fitxer mes proper
[~, iy] = min(abs(camp.y - p(2)));            % fila del fitxer mes propera al dron
[~, iz] = min(abs(camp.z - p(3)));            % alcada del fitxer mes propera al dron

% Escena real: tall horitzontal de la matriu a l'alcada del dron
set(fig.terra, 'CData', log10(double(camp.c(:, :, iz, it))' + 1e-3));
set(fig.trajectoria, 'XData', trajectoria(:, 1), 'YData', trajectoria(:, 2), 'ZData', trajectoria(:, 3));
set(fig.pal,  'XData', p(1) * [1 1], 'YData', p(2) * [1 1], 'ZData', [0 p(3)]);
set(fig.dron, 'XData', p(1), 'YData', p(2), 'ZData', p(3));
set(fig.titol, 'String', {titol, 'Colors del terra: concentracio a l''alcada del dron'});

% Vista lateral: tall vertical de la matriu per la y del dron
if fig.tresD
    set(fig.lateral, 'CData', log10(squeeze(double(camp.c(:, iy, :, it)))' + 1e-3));
    set(fig.trajLateral, 'XData', trajectoria(:, 1), 'YData', trajectoria(:, 3));
    set(fig.dronLateral, 'XData', p(1), 'YData', p(3));
end

% El que sap el dron (si hi ha diverses alcades, es resumeixen en una)
if strcmp(dron.algorisme, 'infotaxis')
    mapa = sum(dron.mapa, 3);                         % probabilitat sumant totes les alcades
    mapa = mapa / max(mapa(:));
else
    mapa = log10(max(dron.mapa, [], 3) + 1e-3);       % concentracio maxima a qualsevol alcada
end
set(fig.mapa, 'CData', mapa);
set(fig.dronMapa, 'XData', p(1), 'YData', p(2));

% Concentracio mesurada (els zeros es dibuixen a baix de tot de l'escala)
set(fig.mesura, 'XData', (1:numel(concMesurada)) * dron.dt, 'YData', max(concMesurada, 0.05));
drawnow;
end
