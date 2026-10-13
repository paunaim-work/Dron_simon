function vis = visualitzacio_actualitzar(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, trobada, final)
%VISUALITZACIO_ACTUALITZAR  Redibuixa un fotograma de l'animacio.

if ~ishandle(vis.fig)
    vis.tancada = true;
    return;
end

t = (0:n - 1)' * dom.dt;
p = hist.pos(n, :);

% --- Ploma (nomes canvia la complexa) ---
if strcmp(pl.tipus, 'complexa')
    set(vis.hImg, 'CData', log10(max(ploma_camp(pl, dom), vis.cmin)));
    th = pl.theta0 + pl.dTheta;
else
    th = dom.theta;
end
v = vis.longFletxa * [cos(th) sin(th)];
set(vis.hVent, 'XData', fletxa(vis.baseVent - v / 2, v, 1), ...
               'YData', fletxa(vis.baseVent - v / 2, v, 2));

% --- Dron, trajectoria i deteccions ---
set(vis.hTraj, 'XData', hist.pos(1:n, 1), 'YData', hist.pos(1:n, 2));
set(vis.hTraj2, 'XData', hist.pos(1:n, 1), 'YData', hist.pos(1:n, 2));
d = hist.det(1:n);
set(vis.hHits, 'XData', hist.pos(d, 1), 'YData', hist.pos(d, 2));
set(vis.hDron, 'XData', p(1), 'YData', p(2));
set(vis.hDron2, 'XData', p(1), 'YData', p(2));

% --- Coneixement del dron ---
switch alg
    case 'infotaxis'
        LP = log10(max(cerc.P, 1e-300));
        mx = max(LP(:));
        set(vis.hMapa, 'CData', LP);
        set(vis.ax2, 'CLim', [mx - 4, mx + 1e-9]);
        [~, im] = max(cerc.P(:));
        [iy, ix] = ind2sub(size(cerc.P), im);
        set(vis.hEst, 'XData', dom.x(ix), 'YData', dom.y(iy));
        set(vis.hSerie, 'XData', t, 'YData', hist.S(1:n));
        estat = sprintf('Entropia S = %.2f   |   decisio: %s', hist.S(n), cerc.mode);
    case 'gradient'
        set(vis.hMapa, 'CData', mapa_mesures_gradient(cerc, vis));
        g = cerc.g;
        if any(g ~= 0)
            v = 1.5 * vis.longFletxa * g / norm(g);
            set(vis.hGrad, 'XData', fletxa(p, v, 1), 'YData', fletxa(p, v, 2));
        else
            set(vis.hGrad, 'XData', NaN, 'YData', NaN);
        end
        set(vis.hSerie, 'XData', t, 'YData', hist.k(1:n));
        estat = sprintf('Mode: %s', cerc.mode);
end
set(vis.hDist, 'XData', t, 'YData', hist.dist(1:n));
tmax = max(t(end), 10);
set(vis.ax3, 'XLim', [0 tmax], 'YLim', [0 1.05 * max([hist.dist(1:n); 1])]);
set(vis.ax4, 'XLim', [0 tmax]);

% --- Titol ---
linia1 = sprintf('Ploma %s  |  algorisme %s  |  t = %.0f s  (pas %d)', ...
                 cfg.tipusPloma, alg, t(end), n - 1);
if trobada
    linia2 = sprintf('FONT TROBADA en %.0f s (%d passos, %d deteccions)', t(end), n - 1, sum(d));
elseif final
    linia2 = sprintf('Font NO trobada (limit de passos). Distancia final: %.0f m', hist.dist(n));
else
    linia2 = estat;
end
title(vis.ax1, {linia1, linia2});

drawnow;
if cfg.visual.pausa > 0
    pause(cfg.visual.pausa);
end
vis.video = video_afegir(vis.video, vis.fig);
end

% ---------------------------------------------------------------------
function c = fletxa(base, v, coord)
%FLETXA  Coordenades d'una fletxa (tija + punta) com una sola linia.
punta = base + v;
L = norm(v);
if L == 0
    c = NaN;
    return;
end
u = v / L;
w = [-u(2) u(1)];
a = punta - 0.3 * L * u + 0.18 * L * w;
b = punta - 0.3 * L * u - 0.18 * L * w;
P = [base; punta; NaN NaN; a; punta; b];
c = P(:, coord);
end
