function res = simular_cerca(cfg)
%SIMULAR_CERCA  Simula la cerca de l'origen d'una fuga amb un dron.
%
%   res = simular_cerca(cfg) executa una cerca completa amb la
%   configuracio cfg (vegeu configuracio_defecte.m) i retorna:
%
%     res.trobada     true si el dron ha arribat a la font
%     res.passos      nombre de passos fets
%     res.temps       temps de cerca (s)
%     res.recorregut  distancia recorreguda (m)
%     res.hist        historial (posicions, comptes, distancia, entropia...)
%
%   A cada pas: (1) la ploma evoluciona, (2) el dron decideix on anar,
%   (3) es mou a una cel.la veina (o es queda quiet), (4) si es a menys
%   del radi de deteccio de la font la cerca acaba; si no, (5) mesura i
%   actualitza el que sap.
%
%   En mode realista (cfg.realista.activat) la ploma evoluciona durant el
%   pas mentre el dron es desplaca i s'atura a mesurar (pas_realista.m), i
%   la bateria pot obligar a acabar la missio (res.motiu = 'bateria').

% Condicions inicials aleatories (si cfg.atzar.activat)
cfg = aplicar_escenari_aleatori(cfg);

% Llavor de la turbulencia i del soroll: si es buida se'n tria una i es
% desa a res.cfg.llavor perque la prova es pugui repetir
if isempty(cfg.llavor)
    cfg.llavor = mod(floor(now * 1e7) + 7919, 2^31 - 1);
end
fixar_llavor(cfg.llavor);
cfg = preparar_realista(cfg);           % limitacions tecniques (si s'activen)
realista = isfield(cfg, 'realista') && cfg.realista.activat;
dom = crear_domini(cfg);
alg = lower(cfg.algorisme);

% La font se situa al centre d'una cel.la de la graella
celFont = posicio_a_cella(cfg.font.pos, dom);
cfg.font.pos = cella_a_posicio(celFont, dom);
posFont = cfg.font.pos;

pl  = ploma_inicialitzar(cfg, dom);
pos = posicio_a_cella(cfg.dron.posInicial, dom);

switch alg
    case 'infotaxis'
        cerc = infotaxis_inicialitzar(cfg, dom, pos);
    case 'gradient'
        cerc = gradient_inicialitzar(cfg, dom);
    otherwise
        error('Algorisme desconegut: %s (usa ''infotaxis'' o ''gradient'')', cfg.algorisme);
end

% --- Historial ---
N = cfg.dron.maxPassos;
hist.pos  = nan(N + 1, 2);   % posicio (m)
hist.k    = zeros(N + 1, 1); % comptes del sensor (suma de sensors per al gradient)
hist.det  = false(N + 1, 1); % hi ha hagut deteccio?
hist.dist = nan(N + 1, 1);   % distancia a la font (m)
hist.S    = nan(N + 1, 1);   % entropia (infotaxis)
hist.mode = zeros(N + 1, 1); % 1 = seguiment, 0 = exploracio (gradient)
hist.bateria = nan(N + 1, 1);  % bateria que queda (%) (mode realista)

n = 1;
hist.pos(n, :) = cella_a_posicio(pos, dom);
hist.dist(n)   = norm(hist.pos(n, :) - posFont);
sens  = [];                     % estat del sensor (mode realista)
bat   = bateria_inicialitzar(cfg, hist.pos(1, :));
motiu = '';
trobada = hist.dist(n) <= dom.radiDet + 1e-9;
if ~trobada
    [hist, cerc, pl, sens] = mesurar_i_aprendre(alg, cerc, pl, pos, pos, sens, realista, cfg, dom, hist, n);
    [bat, ~, hist.bateria(n)] = bateria_actualitzar(bat, [0 0], hist.pos(n, :), dom);
end

vis = [];
if cfg.visual.mostrar
    vis = visualitzacio_inicialitzar(cfg, dom, pl, cerc, posFont, hist, n, alg);
    % Si la simulacio s'interromp (error, Ctrl+C), el video es tanca igualment
    tancarVideo = onCleanup(@() video_tancar_segur(vis.video)); %#ok<NASGU>
end

% --- Bucle de cerca ---
for t = 1:N
    if trobada
        break;
    end
    if ~realista
        pl = ploma_avancar(pl, dom.dt);   % (en mode realista avanca dins de pas_realista)
    end

    switch alg
        case 'infotaxis'
            [mov, cerc] = infotaxis_decidir(cerc, pos, dom);
        case 'gradient'
            [mov, cerc] = gradient_decidir(cerc, pos, dom);
    end
    posAnt = pos;
    pos = pos + mov;

    n = n + 1;
    hist.pos(n, :) = cella_a_posicio(pos, dom);
    hist.dist(n)   = norm(hist.pos(n, :) - posFont);

    if hist.dist(n) <= dom.radiDet + 1e-9
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
        vis = visualitzacio_actualitzar(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, trobada, false);
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
        vis = visualitzacio_actualitzar(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, trobada, true);
    end
    visualitzacio_finalitzar(vis);
end

% --- Resultats ---
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
p = cella_a_posicio(pos, dom);
if strcmp(alg, 'gradient')
    offsets = cerc.offsets;
else
    offsets = [0 0];
end
if realista
    fConc = @(pl, P) ploma_concentracio(pl, P(:, 1), P(:, 2));
    [k, pl, sens] = pas_realista(pl, cella_a_posicio(posAnt, dom), p, offsets, sens, cfg, dom, ...
                                 @ploma_avancar, fConc);
else
    k = sensor_llegir(pl, p + offsets, cfg, dom);
end
switch alg
    case 'infotaxis'
        cerc = infotaxis_actualitzar(cerc, pos, k, dom);
        hist.k(n)   = k;
        hist.det(n) = k > 0;
        hist.S(n)   = cerc.S;
    case 'gradient'
        cerc = gradient_actualitzar(cerc, pos, k, dom);
        hist.k(n)    = sum(k);
        hist.det(n)  = sum(k) >= cerc.llindar;
        hist.mode(n) = strcmp(cerc.mode, 'seguiment');
end
end
