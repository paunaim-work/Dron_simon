function vis = visualitzacio_inicialitzar(cfg, dom, pl, cerc, posFont, hist, n, alg)
%VISUALITZACIO_INICIALITZAR  Crea la figura de l'animacio.
%
%   Panell gran (esquerra): la ploma real, la font, el dron, la seva
%   trajectoria i les deteccions.
%   Panell superior dret: el que "sap" el dron (probabilitat posterior de
%   la posicio de la font per a l'infotaxis, o mapa de mesures per al
%   gradient).
%   Panells inferiors drets: distancia a la font i entropia / comptes.

vis.tancada = false;
vis.alg = alg;
vis.fig = figure('Name', 'Cerca de l''origen d''una fuga amb un dron', ...
                 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', ...
                 'Position', [0.03 0.07 0.94 0.85]);

% Limits de color de la concentracio (escala logaritmica):
% de 0.02 a 100 comptes esperats per pas
kc = cfg.sensor.sensibilitat * dom.dt;
vis.cmin = 0.02 / kc;
vis.clim = log10([vis.cmin 100 / kc]);
vis.longFletxa = 0.08 * max(dom.x(end), dom.y(end));

% ------------------------- Panell 1: escena -------------------------
vis.ax1 = axes('Parent', vis.fig, 'Position', [0.05 0.17 0.50 0.72]);
C = ploma_camp(pl, dom);
vis.hImg = imagesc(dom.xv, dom.yv, log10(max(C, vis.cmin)), 'Parent', vis.ax1);
hold(vis.ax1, 'on');
set(vis.ax1, 'YDir', 'normal', 'CLim', vis.clim);
axis(vis.ax1, 'equal');
axis(vis.ax1, [dom.x(1) dom.x(end) dom.y(1) dom.y(end)]);
colormap(vis.ax1, mapa_colors('ploma'));
cb = colorbar(vis.ax1);
ylabel(cb, 'log_{10} concentracio');

p0 = hist.pos(1, :);
vis.hTraj = plot(vis.ax1, hist.pos(1:n, 1), hist.pos(1:n, 2), '-', 'Color', [0.10 0.35 0.85], 'LineWidth', 1.6);
vis.hHits = plot(vis.ax1, NaN, NaN, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
vis.hIni  = plot(vis.ax1, p0(1), p0(2), '^', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
vis.hFont = plot(vis.ax1, posFont(1), posFont(2), 'p', 'MarkerSize', 20, 'MarkerFaceColor', [0.1 0.85 0.2], 'MarkerEdgeColor', 'k');
vis.hDron = plot(vis.ax1, p0(1), p0(2), 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0 0.9 1], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Fletxa del vent (a dalt a l'esquerra)
vis.baseVent = [dom.x(1) + 0.06 * dom.x(end), dom.y(end) - 0.08 * dom.y(end)];
vis.hVent = plot(vis.ax1, NaN, NaN, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2.5);
text(vis.baseVent(1), vis.baseVent(2) + 0.05 * dom.y(end), 'vent', 'Parent', vis.ax1, 'FontWeight', 'bold');

switch alg
    case 'infotaxis'
        vis.hEst = plot(vis.ax1, NaN, NaN, 'x', 'Color', [0.8 0 0.8], 'MarkerSize', 14, 'LineWidth', 2.5);
        hLeg = [vis.hFont vis.hIni vis.hDron vis.hTraj vis.hHits vis.hEst];
        txtLeg = {'Font real', 'Inici', 'Dron', 'Trajectoria', 'Deteccions', 'Font mes probable'};
    case 'gradient'
        vis.hGrad = plot(vis.ax1, NaN, NaN, '-', 'Color', [0.8 0 0.8], 'LineWidth', 2.5);
        hLeg = [vis.hFont vis.hIni vis.hDron vis.hTraj vis.hHits vis.hGrad];
        txtLeg = {'Font real', 'Inici', 'Dron', 'Trajectoria', 'Deteccions', 'Gradient estimat'};
end
hl = legend(vis.ax1, hLeg, txtLeg, 'Orientation', 'horizontal');
set(hl, 'Position', [0.05 0.025 0.50 0.04]);
xlabel(vis.ax1, 'x (m)');
ylabel(vis.ax1, 'y (m)');

% ------------------- Panell 2: coneixement del dron -------------------
vis.ax2 = axes('Parent', vis.fig, 'Position', [0.62 0.50 0.34 0.42]);
switch alg
    case 'infotaxis'
        vis.hMapa = imagesc(dom.x, dom.y, log10(max(cerc.P, 1e-300)), 'Parent', vis.ax2);
        colormap(vis.ax2, mapa_colors('probabilitat'));
        title(vis.ax2, 'Probabilitat P(r_0) que la font sigui a cada punt (log_{10})');
    case 'gradient'
        vis.hMapa = imagesc(dom.x, dom.y, mapa_mesures_gradient(cerc, vis), 'Parent', vis.ax2);
        colormap(vis.ax2, mapa_colors('ploma'));
        set(vis.ax2, 'CLim', vis.clim);
        title(vis.ax2, 'Concentracio mesurada pel dron (log_{10})');
end
hold(vis.ax2, 'on');
set(vis.ax2, 'YDir', 'normal');
axis(vis.ax2, 'equal');
axis(vis.ax2, [dom.x(1) dom.x(end) dom.y(1) dom.y(end)]);
vis.hTraj2 = plot(vis.ax2, hist.pos(1:n, 1), hist.pos(1:n, 2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
plot(vis.ax2, posFont(1), posFont(2), 'p', 'MarkerSize', 14, 'MarkerFaceColor', [0.1 0.85 0.2], 'MarkerEdgeColor', 'k');
vis.hDron2 = plot(vis.ax2, p0(1), p0(2), 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0 0.9 1], 'MarkerEdgeColor', 'k');

% ------------------- Panells 3 i 4: series temporals -------------------
vis.ax3 = axes('Parent', vis.fig, 'Position', [0.62 0.27 0.34 0.14]);
vis.hDist = plot(vis.ax3, 0, hist.dist(1), '-', 'Color', [0.10 0.35 0.85], 'LineWidth', 1.5);
ylabel(vis.ax3, 'Distancia (m)');
set(vis.ax3, 'YLim', [0 1.05 * max(hist.dist(1), 1)]);
grid(vis.ax3, 'on');

vis.ax4 = axes('Parent', vis.fig, 'Position', [0.62 0.07 0.34 0.14]);
switch alg
    case 'infotaxis'
        vis.hSerie = plot(vis.ax4, 0, hist.S(1), '-', 'Color', [0.55 0.1 0.6], 'LineWidth', 1.5);
        ylabel(vis.ax4, 'Entropia S');
        set(vis.ax4, 'YLim', [0 1.05 * max(hist.S(1), 1)]);
    case 'gradient'
        vis.hSerie = plot(vis.ax4, 0, hist.k(1), '-', 'Color', [0.85 0.3 0.1], 'LineWidth', 1);
        ylabel(vis.ax4, 'Comptes / pas');
end
xlabel(vis.ax4, 'Temps (s)');
grid(vis.ax4, 'on');

% ------------------------------ Video ------------------------------
vis.video = [];
if cfg.visual.gravarVideo
    if exist('VideoWriter', 'file') || exist('VideoWriter', 'class')
        vis.video = VideoWriter(cfg.visual.fitxerVideo, 'MPEG-4');
        vis.video.FrameRate = 20;
        open(vis.video);
    else
        warning('VideoWriter no disponible: no es gravara el video.');
    end
end

vis = visualitzacio_actualitzar(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, false, false);
end
