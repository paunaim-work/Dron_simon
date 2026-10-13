%% SIMULACIO DE LA CERCA DE L'ORIGEN D'UNA FUGA AMB UN DRON AUTONOM
%
%  Simula una fuga (ploma de contaminant) i un dron que en busca l'origen
%  nomes a partir de les lectures del seu sensor. Es pot triar:
%
%    * El model de ploma:
%        'gaussiana' -> ploma gaussiana estacionaria (Pasquill-Gifford)
%        'complexa'  -> ploma turbulenta de "puffs" amb vent variable
%                       (meandreig, rafegues i remolins): sinuosa i
%                       intermitent
%
%    * L'algorisme de cerca:
%        'infotaxis' -> estrategia de l'article de Vergassola, Villermaux
%                       i Shraiman, Nature 445, 406-409 (2007): el dron
%                       es mou on espera guanyar mes informacio sobre la
%                       posicio de la font
%        'gradient'  -> quimiotaxi classica: el dron puja pel gradient de
%                       concentracio mesurat amb 4 sensors en creu
%
%  Executa aquest script (F5). Tots els parametres (domini, vent, sensor,
%  ploma, algorismes...) son a funcions/configuracio_defecte.m i es poden
%  canviar aqui mateix abans de cridar simular_cerca.

clear; close all; clc;
carpeta = fileparts(mfilename('fullpath'));
if isempty(carpeta)
    carpeta = pwd;
end
addpath(fullfile(carpeta, 'funcions'));

cfg = configuracio_defecte();

%% ===================== TRIA LES OPCIONS AQUI =====================
cfg.tipusPloma = 'gaussiana';   % 'gaussiana'  o  'complexa'
cfg.algorisme  = 'infotaxis';   % 'infotaxis'  o  'gradient'

% Escenari aleatori: true = cada execucio sorteja la posicio de la fuga,
% la del dron i el vent (direccio i velocitat); false = sempre el mateix
% escenari (el de configuracio_defecte.m). Que se sorteja i en quins
% marges es tria amb cfg.atzar.* (vegeu configuracio_defecte.m).
cfg.atzar.activat = true;
cfg.atzar.llavor  = [];         % [] = escenari nou cada cop; un numero = repeteix aquell escenari
cfg.llavor        = [];         % turbulencia i soroll del sensor: [] = diferents cada cop; un numero = repetible

% Posa-ho a true per triar les opcions per la finestra de comandes
preguntarOpcions = false;

% Visualitzacio
cfg.visual.mostrar      = true;
cfg.visual.cadaNPassos  = 1;      % puja'l (p. ex. 5) si l'animacio va lenta
cfg.visual.pausa        = 0.01;   % s entre fotogrames (alenteix l'animacio)
cfg.visual.gravarVideo  = false;  % true -> desa un video de l'animacio
cfg.visual.formatVideo  = 'mp4';  % 'mp4', 'avi' o 'gif'
cfg.visual.carpetaVideo = fullfile(carpeta, 'videos');   % on es desen els videos

% Exemples d'altres parametres que pots canviar:
% cfg.atzar.estabilitat = true;       % sorteja tambe la classe d'estabilitat
% cfg.atzar.emissio     = true;       % sorteja tambe la taxa d'emissio
% cfg.atzar.font        = false;      % no moguis la fuga (deixa-la a cfg.font.pos)
% cfg.font.pos          = [60 150];   % on es la fuga (m) si no se sorteja
% cfg.dron.posInicial   = [340 215];  % d'on surt el dron (m) si no se sorteja
% cfg.vent.direccio     = 0;          % cap on bufa el vent (graus) si no se sorteja
% cfg.gauss.estabilitat = 'D';        % classe d'estabilitat A..F
% cfg.sensor.sensibilitat = 0.3;      % sensor menys sensible -> senyal mes diluit
% cfg.infotaxis.saturacio = 1;        % sensor binari (deteccio si/no)
%% =================================================================

if preguntarOpcions
    op = input('Tipus de ploma: 1 = gaussiana, 2 = complexa  [1]: ');
    if isequal(op, 2)
        cfg.tipusPloma = 'complexa';
    else
        cfg.tipusPloma = 'gaussiana';
    end
    op = input('Algorisme: 1 = infotaxis (article), 2 = gradient  [1]: ');
    if isequal(op, 2)
        cfg.algorisme = 'gradient';
    else
        cfg.algorisme = 'infotaxis';
    end
    op = input('Escenari: 1 = aleatori, 2 = fix  [1]: ');
    cfg.atzar.activat = ~isequal(op, 2);
end

fprintf('Ploma: %s | Algorisme: %s\n', cfg.tipusPloma, cfg.algorisme);
res = simular_cerca(cfg);

%% Resum
fprintf('\n=============== RESULTAT ===============\n');
if res.trobada
    fprintf('Font TROBADA\n');
else
    fprintf('Font NO trobada (limit de %d passos)\n', cfg.dron.maxPassos);
end
fprintf('Temps de cerca:       %8.0f s (%.1f min)\n', res.temps, res.temps / 60);
fprintf('Passos:               %8d\n', res.passos);
fprintf('Distancia recorreguda:%8.0f m\n', res.recorregut);
fprintf('Distancia inicial:    %8.0f m\n', res.hist.dist(1));
fprintf('Deteccions:           %8d\n', res.deteccions);
fprintf('========================================\n');
if cfg.atzar.activat
    fprintf('Per repetir exactament aquesta prova:\n   cfg.atzar.llavor = %d;  cfg.llavor = %d;\n', ...
            res.cfg.atzar.llavorUsada, res.cfg.llavor);
else
    fprintf('Per repetir exactament aquesta prova:  cfg.llavor = %d;\n', res.cfg.llavor);
end
