function vis = visualitzacio_actualitzar_3d(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, trobada, final)
%VISUALITZACIO_ACTUALITZAR_3D  Redibuixa un fotograma de l'animacio 3D.

if ~ishandle(vis.fig)
    vis.tancada = true;
    return;
end

t  = (0:n - 1)' * dom.dt;
p  = hist.pos(n, :);
P3 = hist.pos(1:n, :);

% --- Ploma 3D i vent (nomes canvien a la ploma complexa) ---
if strcmp(pl.tipus, 'complexa')
    vis = actualitzar_isosuperficies(vis, ploma_camp_3d(pl, dom), dom);
    th = pl.theta0 + pl.dTheta;
else
    th = dom.theta;
end
v = vis.longFletxa * [cos(th) sin(th) 0];
[fx, fy, fz] = fletxa3(vis.baseVent - v / 2, v);
set(vis.hVent, 'XData', fx, 'YData', fy, 'ZData', fz);

% --- Escena 3D ---
set(vis.hTraj, 'XData', P3(:, 1), 'YData', P3(:, 2), 'ZData', P3(:, 3));
set(vis.hOmbra, 'XData', P3(:, 1), 'YData', P3(:, 2), 'ZData', zeros(n, 1));
d = hist.det(1:n);
set(vis.hHits, 'XData', P3(d, 1), 'YData', P3(d, 2), 'ZData', P3(d, 3));
set(vis.hDron, 'XData', p(1), 'YData', p(2), 'ZData', p(3));
set(vis.hPal, 'XData', p([1 1]), 'YData', p([2 2]), 'ZData', [0 p(3)]);

% --- Tall horitzontal a l'alcada del dron ---
if strcmp(pl.tipus, 'complexa') || p(3) ~= vis.zTall
    set(vis.hTall, 'CData', log10(max(ploma_tall_3d(pl, dom, p(3)), vis.cmin)));
    vis.zTall = p(3);
end
title(vis.ax2, sprintf('Ploma a l''alcada del dron (z = %.0f m)', p(3)));
set(vis.hTraj2, 'XData', P3(:, 1), 'YData', P3(:, 2));
set(vis.hDron2, 'XData', p(1), 'YData', p(2));

% --- El que sap el dron ---
[Mdalt, Mcostat] = mapes_coneixement(cerc, alg, vis);
set(vis.hMapa3, 'CData', Mdalt);
set(vis.hMapa4, 'CData', Mcostat);
set(vis.hTraj3, 'XData', P3(:, 1), 'YData', P3(:, 2));
set(vis.hTraj4, 'XData', P3(:, 1), 'YData', P3(:, 3));
set(vis.hDron3, 'XData', p(1), 'YData', p(2));
set(vis.hDron4, 'XData', p(1), 'YData', p(3));

switch alg
    case 'infotaxis'
        mx = max(Mdalt(:));
        set(vis.ax3, 'CLim', [mx - 4, mx + 1e-9]);
        mx = max(Mcostat(:));
        set(vis.ax4, 'CLim', [mx - 4, mx + 1e-9]);
        [~, im] = max(cerc.P(:));
        [iy, ix, iz] = ind2sub(size(cerc.P), im);
        set(vis.hEst, 'XData', dom.x(ix), 'YData', dom.y(iy), 'ZData', dom.z(iz));
        set(vis.hSerie, 'XData', t, 'YData', hist.S(1:n));
        estat = sprintf('Entropia S = %.2f   |   decisio: %s', hist.S(n), cerc.mode);
    case 'gradient'
        g = cerc.g;
        if any(g ~= 0)
            [fx, fy, fz] = fletxa3(p, 1.5 * vis.longFletxa * g / norm(g));
            set(vis.hGrad, 'XData', fx, 'YData', fy, 'ZData', fz);
        else
            set(vis.hGrad, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
        end
        set(vis.hSerie, 'XData', t, 'YData', hist.k(1:n));
        estat = sprintf('Mode: %s', cerc.mode);
end

% --- Series temporals ---
set(vis.hDist, 'XData', t, 'YData', hist.dist(1:n));
set(vis.hAlc, 'XData', t, 'YData', P3(:, 3));
tmax = max(t(end), 10);
set(vis.ax5, 'XLim', [0 tmax], 'YLim', [0 1.05 * max([hist.dist(1:n); dom.z(end); 1])]);
set(vis.ax6, 'XLim', [0 tmax]);

% --- Titol ---
linia1 = sprintf('Ploma %s 3D  |  %s  |  t = %.0f s (pas %d)  |  dron a z = %.0f m', ...
                 cfg.tipusPloma, alg, t(end), n - 1, p(3));
if ~isnan(hist.bateria(n))
    linia1 = sprintf('%s  |  bateria %.0f %%', linia1, hist.bateria(n));
end
if trobada
    linia2 = sprintf('FONT TROBADA en %.0f s (%d passos, %d deteccions) | dron a z = %.0f m, fuga a z = %.0f m', ...
                     t(end), n - 1, sum(d), p(3), posFont(3));
elseif strcmp(vis.motiu, 'bateria')
    linia2 = sprintf('BATERIA BAIXA: el dron torna a la base. Distancia a la fuga: %.0f m', hist.dist(n));
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
function [X, Y, Z] = fletxa3(base, v)
%FLETXA3  Fletxa 3D (tija + punta en el pla que conte v i l'eix z, o el
%   pla horitzontal si v es vertical) com una sola linia amb NaN.
L = norm(v);
if L == 0
    X = NaN; Y = NaN; Z = NaN;
    return;
end
u = v / L;
w = cross(u, [0 0 1]);
if norm(w) < 1e-6
    w = [1 0 0];
end
w = w / norm(w);
punta = base + v;
a = punta - 0.3 * L * u + 0.18 * L * w;
b = punta - 0.3 * L * u - 0.18 * L * w;
P = [base; punta; NaN NaN NaN; a; punta; b];
X = P(:, 1); Y = P(:, 2); Z = P(:, 3);
end
