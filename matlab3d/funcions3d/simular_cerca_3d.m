function res = simular_cerca_3d(cfg)
%SIMULAR_CERCA_3D  Simula la cerca 3D de l'origen d'una fuga amb un dron.
%
%   res = simular_cerca_3d(cfg) amb cfg de configuracio_defecte_3d.m.
%   Retorna res.trobada, res.passos, res.temps (s), res.recorregut (m),
%   res.deteccions i res.hist (posicions [x y z], comptes, distancia,
%   entropia...).

% Condicions inicials aleatories (si cfg.atzar.activat)
cfg = aplicar_escenari_aleatori_3d(cfg);

% Llavor de la turbulencia i del soroll ([] -> se'n tria una i es desa)
if isempty(cfg.llavor)
    cfg.llavor = mod(floor(now * 1e7) + 7919, 2^31 - 1);
end
fixar_llavor(cfg.llavor);
dom = crear_domini_3d(cfg);
alg = lower(cfg.algorisme);

% La font se situa al centre d'una cel.la de la graella
celFont = posicio_a_cella_3d(cfg.font.pos, dom);
cfg.font.pos = cella_a_posicio_3d(celFont, dom);
posFont = cfg.font.pos;

pl  = ploma_inicialitzar_3d(cfg, dom);
pos = posicio_a_cella_3d(cfg.dron.posInicial, dom);

switch alg
    case 'infotaxis'
        cerc = infotaxis_inicialitzar_3d(cfg, dom, pos);
    case 'gradient'
        cerc = gradient_inicialitzar_3d(cfg, dom);
    otherwise
        error('Algorisme desconegut: %s (usa ''infotaxis'' o ''gradient'')', cfg.algorisme);
end

N = cfg.dron.maxPassos;
hist.pos  = nan(N + 1, 3);
hist.k    = zeros(N + 1, 1);
hist.det  = false(N + 1, 1);
hist.dist = nan(N + 1, 1);
hist.S    = nan(N + 1, 1);
hist.mode = zeros(N + 1, 1);

n = 1;
hist.pos(n, :) = cella_a_posicio_3d(pos, dom);
hist.dist(n)   = norm(hist.pos(n, :) - posFont);
trobada = hist.dist(n) <= dom.radiDet + 1e-9;
if ~trobada
    [hist, cerc] = mesurar_i_aprendre(alg, cerc, pl, pos, cfg, dom, hist, n);
end

vis = [];
if cfg.visual.mostrar
    vis = visualitzacio_inicialitzar_3d(cfg, dom, pl, cerc, posFont, hist, n, alg);
    % Si la simulacio s'interromp (error, Ctrl+C), el video es tanca igualment
    tancarVideo = onCleanup(@() video_tancar_segur(vis.video)); %#ok<NASGU>
end

for t = 1:N
    if trobada
        break;
    end
    pl = ploma_avancar_3d(pl, dom.dt);

    switch alg
        case 'infotaxis'
            [mov, cerc] = infotaxis_decidir_3d(cerc, pos, dom);
        case 'gradient'
            [mov, cerc] = gradient_decidir_3d(cerc, pos, dom);
    end
    pos = pos + mov;

    n = n + 1;
    hist.pos(n, :) = cella_a_posicio_3d(pos, dom);
    hist.dist(n)   = norm(hist.pos(n, :) - posFont);

    if hist.dist(n) <= dom.radiDet + 1e-9
        trobada = true;
        hist.S(n) = 0;
    else
        [hist, cerc] = mesurar_i_aprendre(alg, cerc, pl, pos, cfg, dom, hist, n);
    end

    if ~isempty(vis) && (trobada || mod(t, cfg.visual.cadaNPassos) == 0)
        vis = visualitzacio_actualitzar_3d(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, trobada, false);
        if vis.tancada
            fprintf('Finestra tancada: simulacio aturada.\n');
            break;
        end
    end
end

if ~isempty(vis)
    if ~vis.tancada
        vis = visualitzacio_actualitzar_3d(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, trobada, true);
    end
    visualitzacio_finalitzar(vis);
end

camps = fieldnames(hist);
for i = 1:numel(camps)
    hist.(camps{i}) = hist.(camps{i})(1:n, :);
end
res.trobada    = trobada;
res.passos     = n - 1;
res.temps      = (n - 1) * dom.dt;
res.recorregut = sum(sqrt(sum(diff(hist.pos, 1, 1).^2, 2)));
res.deteccions = sum(hist.det);
res.hist       = hist;
res.cfg        = cfg;
res.dom        = dom;
res.cercador   = cerc;
end

% ---------------------------------------------------------------------
function [hist, cerc] = mesurar_i_aprendre(alg, cerc, pl, pos, cfg, dom, hist, n)
p = cella_a_posicio_3d(pos, dom);
switch alg
    case 'infotaxis'
        k = sensor_llegir_3d(pl, p, cfg, dom);
        cerc = infotaxis_actualitzar_3d(cerc, pos, k, dom);
        hist.k(n)   = k;
        hist.det(n) = k > 0;
        hist.S(n)   = cerc.S;
    case 'gradient'
        k = sensor_llegir_3d(pl, p + cerc.offsets, cfg, dom);
        cerc = gradient_actualitzar_3d(cerc, pos, k, dom);
        hist.k(n)    = sum(k);
        hist.det(n)  = sum(k) >= cerc.llindar;
        hist.mode(n) = strcmp(cerc.mode, 'seguiment');
end
end
