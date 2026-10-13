%% SIMULACIO 3D DE LA CERCA DE L'ORIGEN D'UNA FUGA AMB UN DRON AUTONOM
%
%  Versio tridimensional de ../matlab/main_simulacio.m: la fuga pot ser
%  a qualsevol alcada, la ploma es 3D i el dron es mou en x, y i z.
%
%    * Ploma:      'gaussiana' (estacionaria, amb reflexio al terra)
%                  'complexa'  (puffs 3D amb vent variable i perfil vertical)
%    * Algorisme:  'infotaxis' (Vergassola et al., Nature 2007)
%                  'gradient'  (quimiotaxi amb 7 sensors)
%
%  Executa aquest script (F5). La vista 3D es pot girar amb el ratoli
%  (icona de rotacio de la barra d'eines de la figura).
%
%  Aquesta versio fa servir algunes funcions generals de la versio 2D
%  (../matlab/funcions): cal mantenir les carpetes matlab i matlab3d
%  juntes, una al costat de l'altra.

clear; close all; clc;
carpeta = fileparts(mfilename('fullpath'));
if isempty(carpeta)
    carpeta = pwd;
end
addpath(fullfile(carpeta, 'funcions3d'));
addpath(fullfile(carpeta, '..', 'matlab', 'funcions'));

cfg = configuracio_defecte_3d();

%% ===================== TRIA LES OPCIONS AQUI =====================
cfg.tipusPloma = 'gaussiana';   % 'gaussiana'  o  'complexa'
cfg.algorisme  = 'infotaxis';   % 'infotaxis'  o  'gradient'
cfg.llavor     = 1;             % canvia-la per veure una altra realitzacio ([] = aleatoria)

% Posa-ho a true per triar les opcions per la finestra de comandes
preguntarOpcions = false;

% Visualitzacio
cfg.visual.mostrar     = true;
cfg.visual.cadaNPassos = 1;      % puja'l (p. ex. 3) si l'animacio va lenta
cfg.visual.pausa       = 0;      % s entre fotogrames (alenteix l'animacio)
cfg.visual.gravarVideo = false;  % true -> desa cerca_dron_3d.mp4

% Exemples d'altres parametres que pots canviar:
% cfg.font.pos           = [60 150 10];  % on es la fuga [x y z] (m)
% cfg.dron.posInicial    = [340 215 30]; % d'on surt el dron [x y z] (m)
% cfg.domini.Lz          = 40;           % alcada maxima de vol (m)
% cfg.vent.direccio      = 0;            % cap on bufa el vent (graus)
% cfg.gauss.estabilitat  = 'D';          % classe d'estabilitat A..F
% cfg.visual.vista       = [-35 25];     % angle de la camera [azimut elevacio]
% cfg.visual.exageracioVertical = 4;     % exageracio de l'eix z al dibuix
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
end

fprintf('Ploma: %s 3D | Algorisme: %s\n', cfg.tipusPloma, cfg.algorisme);
res = simular_cerca_3d(cfg);

%% Resum
fprintf('\n=============== RESULTAT ===============\n');
if res.trobada
    fprintf('Font TROBADA a [%.0f %.0f %.0f] m\n', res.cfg.font.pos);
else
    fprintf('Font NO trobada (limit de %d passos)\n', cfg.dron.maxPassos);
end
fprintf('Temps de cerca:       %8.0f s (%.1f min)\n', res.temps, res.temps / 60);
fprintf('Passos:               %8d\n', res.passos);
fprintf('Distancia recorreguda:%8.0f m\n', res.recorregut);
fprintf('Distancia inicial:    %8.0f m\n', res.hist.dist(1));
fprintf('Deteccions:           %8d\n', res.deteccions);
fprintf('========================================\n');
