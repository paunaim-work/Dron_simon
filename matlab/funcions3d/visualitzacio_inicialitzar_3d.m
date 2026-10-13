function vis = visualitzacio_inicialitzar_3d(cfg, dom, pl, cerc, posFont, hist, n, alg)
%VISUALITZACIO_INICIALITZAR_3D  Crea la figura de l'animacio 3D.
%
%   Esquerra: escena 3D (ploma com a superficies d'igual concentracio,
%   font, dron, trajectoria i la seva ombra al terra, deteccions). Es pot
%   girar amb el ratoli.
%   Dreta, a dalt: tall horitzontal de la ploma a l'alcada del dron i el
%   que sap el dron vist des de dalt (x-y).
%   Dreta, al mig: el que sap el dron vist de costat (x-z).
%   Dreta, a baix: distancia a la font, alcada del dron i entropia/comptes.

vis.tancada = false;
vis.alg = alg;
vis.fig = figure('Name', 'Cerca 3D de l''origen d''una fuga amb un dron', ...
                 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', ...
                 'Position', [0.02 0.06 0.96 0.86]);

kc = cfg.sensor.sensibilitat * dom.dt;
vis.cmin = 0.02 / kc;
vis.clim = log10([vis.cmin 100 / kc]);
vis.nivells = cfg.visual.nivellsIso / kc;   % isosuperficies (en comptes esperats per pas)
vis.longFletxa = 0.08 * max(dom.x(end), dom.y(end));
Lx = dom.x(end); Ly = dom.y(end); Lz = dom.z(end);
p0 = hist.pos(1, :);

% ------------------------- Panell 1: escena 3D -------------------------
vis.ax1 = axes('Parent', vis.fig, 'Position', [0.06 0.17 0.41 0.72]);
hold(vis.ax1, 'on');
patch('Parent', vis.ax1, 'XData', [0 Lx Lx 0], 'YData', [0 0 Ly Ly], 'ZData', [0 0 0 0], ...
      'FaceColor', [0.93 0.93 0.89], 'EdgeColor', 'none');
vis.hIso(1) = patch('Parent', vis.ax1, 'Faces', [1 2 3], 'Vertices', nan(3, 3), ...
                    'FaceColor', [1 0.78 0.30], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
vis.hIso(2) = patch('Parent', vis.ax1, 'Faces', [1 2 3], 'Vertices', nan(3, 3), ...
                    'FaceColor', [0.85 0.15 0.05], 'EdgeColor', 'none', 'FaceAlpha', 0.55);
vis = actualitzar_isosuperficies(vis, ploma_camp_3d(pl, dom), dom);

vis.hOmbra = plot3(vis.ax1, p0(1), p0(2), 0, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1);
vis.hTraj  = plot3(vis.ax1, p0(1), p0(2), p0(3), '-', 'Color', [0.10 0.35 0.85], 'LineWidth', 1.8);
vis.hHits  = plot3(vis.ax1, NaN, NaN, NaN, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
vis.hIni   = plot3(vis.ax1, p0(1), p0(2), p0(3), '^', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
plot3(vis.ax1, posFont([1 1]), posFont([2 2]), [0 posFont(3)], '--', 'Color', [0.1 0.6 0.1]);
vis.hFont  = plot3(vis.ax1, posFont(1), posFont(2), posFont(3), 'p', 'MarkerSize', 18, ...
                   'MarkerFaceColor', [0.1 0.85 0.2], 'MarkerEdgeColor', 'k');
vis.hPal   = plot3(vis.ax1, p0([1 1]), p0([2 2]), [0 p0(3)], ':', 'Color', [0 0.55 0.7], 'LineWidth', 1.2);
vis.hDron  = plot3(vis.ax1, p0(1), p0(2), p0(3), 'o', 'MarkerSize', 11, ...
                   'MarkerFaceColor', [0 0.9 1], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

vis.baseVent = [0.09 * Lx, 0.88 * Ly, Lz];   % centre de la fletxa del vent
vis.hVent = plot3(vis.ax1, NaN, NaN, NaN, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2.5);
text(vis.baseVent(1), vis.baseVent(2), Lz + 0.12 * Lz, 'vent', 'Parent', vis.ax1, ...
     'FontWeight', 'bold', 'HorizontalAlignment', 'center');

switch alg
    case 'infotaxis'
        vis.hEst = plot3(vis.ax1, NaN, NaN, NaN, 'x', 'Color', [0.8 0 0.8], 'MarkerSize', 14, 'LineWidth', 2.5);
        hExtra = vis.hEst;
        txtExtra = 'Font mes probable';
    case 'gradient'
        vis.hGrad = plot3(vis.ax1, NaN, NaN, NaN, '-', 'Color', [0.8 0 0.8], 'LineWidth', 2.5);
        hExtra = vis.hGrad;
        txtExtra = 'Gradient estimat';
end

set(vis.ax1, 'XLim', [0 Lx], 'YLim', [0 Ly], 'ZLim', [0 Lz]);
daspect(vis.ax1, [1 1 1 / cfg.visual.exageracioVertical]);
view(vis.ax1, cfg.visual.vista);
grid(vis.ax1, 'on');
box(vis.ax1, 'on');
xlabel(vis.ax1, 'x (m)');
ylabel(vis.ax1, 'y (m)');
zlabel(vis.ax1, sprintf('alcada z (m)  [x%g]', cfg.visual.exageracioVertical));
try
    camlight(vis.ax1, 'headlight');
    set(vis.hIso, 'FaceLighting', 'gouraud');
catch
end
hl = legend(vis.ax1, [vis.hFont vis.hIni vis.hDron vis.hTraj vis.hHits hExtra vis.hIso(1) vis.hIso(2)], ...
            {'Font real', 'Inici', 'Dron', 'Trajectoria', 'Deteccions', txtExtra, ...
             sprintf('Ploma (%g comptes/pas)', cfg.visual.nivellsIso(1)), ...
             sprintf('Ploma (%g comptes/pas)', cfg.visual.nivellsIso(2))});
try
    set(hl, 'NumColumns', 4);
catch
end
set(hl, 'Position', [0.04 0.015 0.43 0.07]);

% -------------- Panell 2: tall horitzontal a l'alcada del dron --------------
vis.ax2 = axes('Parent', vis.fig, 'Position', [0.52 0.58 0.21 0.34]);
vis.zTall = p0(3);
vis.hTall = imagesc(dom.xv, dom.yv, log10(max(ploma_tall_3d(pl, dom, p0(3)), vis.cmin)), 'Parent', vis.ax2);
hold(vis.ax2, 'on');
colormap(vis.ax2, mapa_colors('ploma'));
set(vis.ax2, 'YDir', 'normal', 'CLim', vis.clim);
axis(vis.ax2, 'equal');
axis(vis.ax2, [0 Lx 0 Ly]);
vis.hTraj2 = plot(vis.ax2, p0(1), p0(2), '-', 'Color', [0.10 0.35 0.85], 'LineWidth', 1);
plot(vis.ax2, posFont(1), posFont(2), 'p', 'MarkerSize', 12, 'MarkerFaceColor', [0.1 0.85 0.2], 'MarkerEdgeColor', 'k');
vis.hDron2 = plot(vis.ax2, p0(1), p0(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.9 1], 'MarkerEdgeColor', 'k');

% ------------- Panells 3 i 4: el que sap el dron (dalt i costat) -------------
vis.ax3 = axes('Parent', vis.fig, 'Position', [0.77 0.58 0.21 0.34]);
vis.ax4 = axes('Parent', vis.fig, 'Position', [0.52 0.37 0.46 0.13]);
[Mdalt, Mcostat] = mapes_coneixement(cerc, alg, vis);
vis.hMapa3 = imagesc(dom.x, dom.y, Mdalt, 'Parent', vis.ax3);
vis.hMapa4 = imagesc(dom.x, dom.z, Mcostat, 'Parent', vis.ax4);
switch alg
    case 'infotaxis'
        colormap(vis.ax3, mapa_colors('probabilitat'));
        colormap(vis.ax4, mapa_colors('probabilitat'));
        title(vis.ax3, 'P(r_0) vista des de dalt (log_{10})');
        title(vis.ax4, 'P(r_0) vista de costat, x-z (log_{10})');
    case 'gradient'
        colormap(vis.ax3, mapa_colors('ploma'));
        colormap(vis.ax4, mapa_colors('ploma'));
        set([vis.ax3 vis.ax4], 'CLim', vis.clim);
        title(vis.ax3, 'Mesures del dron, des de dalt (log_{10})');
        title(vis.ax4, 'Mesures del dron, de costat x-z (log_{10})');
end
hold(vis.ax3, 'on');
hold(vis.ax4, 'on');
set([vis.ax3 vis.ax4], 'YDir', 'normal');
axis(vis.ax3, 'equal');
axis(vis.ax3, [0 Lx 0 Ly]);
set(vis.ax4, 'XLim', [0 Lx], 'YLim', [-dom.dz / 2, Lz + dom.dz / 2]);
ylabel(vis.ax4, 'z (m)');
vis.hTraj3 = plot(vis.ax3, p0(1), p0(2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
vis.hTraj4 = plot(vis.ax4, p0(1), p0(3), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
plot(vis.ax3, posFont(1), posFont(2), 'p', 'MarkerSize', 12, 'MarkerFaceColor', [0.1 0.85 0.2], 'MarkerEdgeColor', 'k');
plot(vis.ax4, posFont(1), posFont(3), 'p', 'MarkerSize', 12, 'MarkerFaceColor', [0.1 0.85 0.2], 'MarkerEdgeColor', 'k');
vis.hDron3 = plot(vis.ax3, p0(1), p0(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.9 1], 'MarkerEdgeColor', 'k');
vis.hDron4 = plot(vis.ax4, p0(1), p0(3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.9 1], 'MarkerEdgeColor', 'k');

% ------------------------- Panells 5 i 6: series -------------------------
vis.ax5 = axes('Parent', vis.fig, 'Position', [0.52 0.21 0.46 0.10]);
hold(vis.ax5, 'on');
vis.hDist = plot(vis.ax5, 0, hist.dist(1), '-', 'Color', [0.10 0.35 0.85], 'LineWidth', 1.5);
vis.hAlc  = plot(vis.ax5, 0, p0(3), '-', 'Color', [0.1 0.6 0.3], 'LineWidth', 1.5);
ylabel(vis.ax5, 'm');
legend(vis.ax5, {'Distancia a la font', 'Alcada del dron'}, 'Location', 'northeast');
grid(vis.ax5, 'on');

vis.ax6 = axes('Parent', vis.fig, 'Position', [0.52 0.06 0.46 0.10]);
switch alg
    case 'infotaxis'
        vis.hSerie = plot(vis.ax6, 0, hist.S(1), '-', 'Color', [0.55 0.1 0.6], 'LineWidth', 1.5);
        ylabel(vis.ax6, 'Entropia S');
        set(vis.ax6, 'YLim', [0 1.05 * max(hist.S(1), 1)]);
    case 'gradient'
        vis.hSerie = plot(vis.ax6, 0, hist.k(1), '-', 'Color', [0.85 0.3 0.1], 'LineWidth', 1);
        ylabel(vis.ax6, 'Comptes / pas');
end
xlabel(vis.ax6, 'Temps (s)');
grid(vis.ax6, 'on');

% ------------------------------ Video ------------------------------
vis.video = video_obrir(cfg, '3d');

vis = visualitzacio_actualitzar_3d(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, false, false);
end
