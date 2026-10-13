%% SIMULACIO DE LA CERCA DE L'ORIGEN D'UNA FUGA AMB UN DRON AUTONOM
%
%  Simula una fuga (ploma de contaminant) i un dron que en busca l'origen
%  nomes a partir de les lectures del seu sensor. Es pot triar:
%
%    * La dimensio:
%        '2d' -> el dron vola a una alcada fixa (cerca en el pla x-y)
%        '3d' -> la fuga pot ser a qualsevol alcada i el dron tambe es mou
%                amunt i avall (cerca en el volum x-y-z)
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
%                       concentracio mesurat amb sensors en creu
%
%  Executa aquest script (F5). Tots els parametres son a
%  funcions/configuracio_defecte.m (2D) i funcions3d/configuracio_defecte_3d.m
%  (3D), i es poden canviar mes avall, a "ALTRES PARAMETRES".

clear; close all; clc;
carpeta = fileparts(mfilename('fullpath'));
if isempty(carpeta)
    carpeta = pwd;
end
addpath(fullfile(carpeta, 'funcions'));
addpath(fullfile(carpeta, 'funcions3d'));

%% ===================== TRIA LES OPCIONS AQUI =====================
dimensio   = '2d';          % '2d'  o  '3d'
tipusPloma = 'gaussiana';   % 'gaussiana'  o  'complexa'
algorisme  = 'infotaxis';   % 'infotaxis'  o  'gradient'

% Escenari aleatori: true = cada execucio sorteja la posicio de la fuga i
% del dron (en 3D tambe les alcades) i el vent; false = sempre el mateix
% escenari (el del fitxer de configuracio).
escenariAleatori = true;
llavorEscenari   = [];      % [] = escenari nou cada cop; un numero = repeteix aquell escenari
llavor           = [];      % turbulencia i soroll del sensor: [] = diferents cada cop; un numero = repetible

% Posa-ho a true per triar les opcions per la finestra de comandes
preguntarOpcions = false;

% Visualitzacio i video
mostrarAnimacio = true;
cadaNPassos     = 1;        % redibuixa cada N passos (puja'l si va lent; en 3D prova 3)
gravarVideo     = false;    % true -> desa un video a la carpeta matlab/videos
formatVideo     = 'mp4';    % 'mp4', 'avi' o 'gif'
%% =================================================================

if preguntarOpcions
    op = input('Dimensio: 1 = 2D, 2 = 3D  [1]: ');
    if isequal(op, 2), dimensio = '3d'; else, dimensio = '2d'; end
    op = input('Tipus de ploma: 1 = gaussiana, 2 = complexa  [1]: ');
    if isequal(op, 2), tipusPloma = 'complexa'; else, tipusPloma = 'gaussiana'; end
    op = input('Algorisme: 1 = infotaxis (article), 2 = gradient  [1]: ');
    if isequal(op, 2), algorisme = 'gradient'; else, algorisme = 'infotaxis'; end
    op = input('Escenari: 1 = aleatori, 2 = fix  [1]: ');
    escenariAleatori = ~isequal(op, 2);
end

switch lower(dimensio)
    case '2d'
        cfg = configuracio_defecte();
    case '3d'
        cfg = configuracio_defecte_3d();
    otherwise
        error('Dimensio desconeguda: %s (usa ''2d'' o ''3d'')', dimensio);
end
cfg.tipusPloma          = tipusPloma;
cfg.algorisme           = algorisme;
cfg.atzar.activat       = escenariAleatori;
cfg.atzar.llavor        = llavorEscenari;
cfg.llavor              = llavor;
cfg.visual.mostrar      = mostrarAnimacio;
cfg.visual.cadaNPassos  = cadaNPassos;
cfg.visual.gravarVideo  = gravarVideo;
cfg.visual.formatVideo  = formatVideo;
cfg.visual.carpetaVideo = fullfile(carpeta, 'videos');

%% ===================== ALTRES PARAMETRES (opcional) =====================
% Treu el % del davant per canviar-los. Les posicions son [x y] en 2D i
% [x y z] en 3D (nomes s'usen si l'escenari no es aleatori o si aquell
% parametre no se sorteja).
% cfg.atzar.estabilitat     = true;        % sorteja tambe la classe d'estabilitat
% cfg.atzar.emissio         = true;        % sorteja tambe la taxa d'emissio
% cfg.atzar.font            = false;       % no moguis la fuga (deixa-la a cfg.font.pos)
% cfg.font.pos              = [60 150];    % 2D: [x y];  3D: p. ex. [60 150 10]
% cfg.dron.posInicial       = [340 215];   % 2D: [x y];  3D: p. ex. [340 215 30]
% cfg.vent.direccio         = 0;           % cap on bufa el vent (graus)
% cfg.gauss.estabilitat     = 'D';         % classe d'estabilitat A..F
% cfg.sensor.sensibilitat   = 0.3;         % sensor menys sensible -> senyal mes diluit
% cfg.infotaxis.saturacio   = 1;           % sensor binari (deteccio si/no)
% cfg.visual.pausa          = 0.05;        % alenteix l'animacio (s entre fotogrames)
% cfg.dron.toleranciaAlcada = 5;           % (3D) diferencia d'alcada admesa per trobar la fuga
% cfg.visual.exageracioVertical = 4;       % (3D) exageracio de l'eix z al dibuix
%% =========================================================================

fprintf('Simulacio %s | Ploma: %s | Algorisme: %s\n', upper(dimensio), cfg.tipusPloma, cfg.algorisme);
if strcmpi(dimensio, '3d')
    res = simular_cerca_3d(cfg);
else
    res = simular_cerca(cfg);
end

%% Resum
fprintf('\n=============== RESULTAT ===============\n');
if res.trobada
    fprintf('Font TROBADA a [%s] m (dron a [%s] m)\n', strtrim(sprintf('%.0f ', res.cfg.font.pos)), ...
            strtrim(sprintf('%.0f ', res.hist.pos(end, :))));
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
    fprintf('Per repetir exactament aquesta prova:\n   llavorEscenari = %d;  llavor = %d;\n', ...
            res.cfg.atzar.llavorUsada, res.cfg.llavor);
else
    fprintf('Per repetir exactament aquesta prova:  llavor = %d;\n', res.cfg.llavor);
end
