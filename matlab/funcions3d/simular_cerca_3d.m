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
cfg = preparar_realista(cfg);           % limitacions tecniques (si s'activen)
realista = isfield(cfg, 'realista') && cfg.realista.activat;
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
hist.bateria = nan(N + 1, 1);  % bateria que queda (%) (mode realista)

n = 1;
hist.pos(n, :) = cella_a_posicio_3d(pos, dom);
hist.dist(n)   = norm(hist.pos(n, :) - posFont);
sens  = [];                     % estat del sensor (mode realista)
bat   = bateria_inicialitzar(cfg, hist.pos(1, :));
motiu = '';
trobada = font_trobada_3d(hist.pos(n, :), posFont, dom);
if ~trobada
    [hist, cerc, pl, sens] = mesurar_i_aprendre(alg, cerc, pl, pos, pos, sens, realista, cfg, dom, hist, n);
    [bat, ~, hist.bateria(n)] = bateria_actualitzar(bat, [0 0 0], hist.pos(n, :), dom);
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
    if ~realista
        pl = ploma_avancar_3d(pl, dom.dt);   % (en mode realista avanca dins de pas_realista)
    end

    switch alg
        case 'infotaxis'
            [mov, cerc] = infotaxis_decidir_3d(cerc, pos, dom);
        case 'gradient'
            [mov, cerc] = gradient_decidir_3d(cerc, pos, dom);
    end
    posAnt = pos;
    pos = pos + mov;

    n = n + 1;
    hist.pos(n, :) = cella_a_posicio_3d(pos, dom);
    hist.dist(n)   = norm(hist.pos(n, :) - posFont);

    if font_trobada_3d(hist.pos(n, :), posFont, dom)
        trobada = true;
        hist.S(n) = 0;
    else
        [hist, cerc, pl, sens] = mesurar_i_aprendre(alg, cerc, pl, posAnt, pos, sens, realista, cfg, dom, hist, n);
    end

    % Bateria: si no arriba per continuar i tornar a la base, s'acaba
    [bat, calTornar, hist.bateria(n)] = bateria_actualitzar(bat, mov, hist.pos(n, :), dom);
    if ~trobada && calTornar
        motiu = 'bateria';
    end

    if ~isempty(vis) && (trobada || ~isempty(motiu) || mod(t, cfg.visual.cadaNPassos) == 0)
        vis.motiu = motiu;
        vis = visualitzacio_actualitzar_3d(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, trobada, false);
        if vis.tancada
            fprintf('Finestra tancada: simulacio aturada.\n');
            motiu = 'finestra';
            break;
        end
    end
    if ~isempty(motiu)
        break;
    end
end
if trobada
    motiu = 'trobada';
elseif isempty(motiu)
    motiu = 'limit';
end

if ~isempty(vis)
    vis.motiu = motiu;
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
res.motiu      = motiu;          % 'trobada', 'bateria', 'limit' o 'finestra'
res.bateria    = hist.bateria(n);
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
function [hist, cerc, pl, sens] = mesurar_i_aprendre(alg, cerc, pl, posAnt, pos, sens, realista, cfg, dom, hist, n)
p = cella_a_posicio_3d(pos, dom);
if strcmp(alg, 'gradient')
    offsets = cerc.offsets;
else
    offsets = [0 0 0];
end
if realista
    fConc = @(pl, P) ploma_concentracio_3d(pl, P(:, 1), P(:, 2), P(:, 3));
    [k, pl, sens] = pas_realista(pl, cella_a_posicio_3d(posAnt, dom), p, offsets, sens, cfg, dom, ...
                                 @ploma_avancar_3d, fConc);
else
    k = sensor_llegir_3d(pl, p + offsets, cfg, dom);
end
switch alg
    case 'infotaxis'
        cerc = infotaxis_actualitzar_3d(cerc, pos, k, dom);
        hist.k(n)   = k;
        hist.det(n) = k > 0;
        hist.S(n)   = cerc.S;
    case 'gradient'
        cerc = gradient_actualitzar_3d(cerc, pos, k, dom);
        hist.k(n)    = sum(k);
        hist.det(n)  = sum(k) >= cerc.llindar;
        hist.mode(n) = strcmp(cerc.mode, 'seguiment');
end
end
