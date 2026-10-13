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

fixar_llavor(cfg.llavor);
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

n = 1;
hist.pos(n, :) = cella_a_posicio(pos, dom);
hist.dist(n)   = norm(hist.pos(n, :) - posFont);
trobada = hist.dist(n) <= dom.radiDet + 1e-9;
if ~trobada
    [hist, cerc] = mesurar_i_aprendre(alg, cerc, pl, pos, cfg, dom, hist, n);
end

vis = [];
if cfg.visual.mostrar
    vis = visualitzacio_inicialitzar(cfg, dom, pl, cerc, posFont, hist, n, alg);
end

% --- Bucle de cerca ---
for t = 1:N
    if trobada
        break;
    end
    pl = ploma_avancar(pl, dom.dt);

    switch alg
        case 'infotaxis'
            [mov, cerc] = infotaxis_decidir(cerc, pos, dom);
        case 'gradient'
            [mov, cerc] = gradient_decidir(cerc, pos, dom);
    end
    pos = pos + mov;

    n = n + 1;
    hist.pos(n, :) = cella_a_posicio(pos, dom);
    hist.dist(n)   = norm(hist.pos(n, :) - posFont);

    if hist.dist(n) <= dom.radiDet + 1e-9
        trobada = true;
        hist.S(n) = 0;
    else
        [hist, cerc] = mesurar_i_aprendre(alg, cerc, pl, pos, cfg, dom, hist, n);
    end

    if ~isempty(vis) && (trobada || mod(t, cfg.visual.cadaNPassos) == 0)
        vis = visualitzacio_actualitzar(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, trobada, false);
        if vis.tancada
            fprintf('Finestra tancada: simulacio aturada.\n');
            break;
        end
    end
end

if ~isempty(vis) && ~vis.tancada
    vis = visualitzacio_actualitzar(vis, cfg, dom, pl, cerc, posFont, hist, n, alg, trobada, true);
    visualitzacio_finalitzar(vis);
end

% --- Resultats ---
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
p = cella_a_posicio(pos, dom);
switch alg
    case 'infotaxis'
        k = sensor_llegir(pl, p, cfg, dom);
        cerc = infotaxis_actualitzar(cerc, pos, k, dom);
        hist.k(n)   = k;
        hist.det(n) = k > 0;
        hist.S(n)   = cerc.S;
    case 'gradient'
        k = sensor_llegir(pl, p + cerc.offsets, cfg, dom);
        cerc = gradient_actualitzar(cerc, pos, k, dom);
        hist.k(n)    = sum(k);
        hist.det(n)  = sum(k) >= cerc.llindar;
        hist.mode(n) = strcmp(cerc.mode, 'seguiment');
end
end
