function cfg = configuracio_defecte_3d()
%CONFIGURACIO_DEFECTE_3D  Parametres per defecte de la simulacio 3D.
%
%   Igual que la versio 2D, pero ara la fuga pot estar a qualsevol
%   alcada, la ploma es tridimensional i el dron es mou tambe amunt i
%   avall. Unitats: metres (m) i segons (s).

% =====================================================================
%  OPCIONS PRINCIPALS
% =====================================================================
cfg.tipusPloma = 'gaussiana';   % 'gaussiana' | 'complexa'
cfg.algorisme  = 'infotaxis';   % 'infotaxis' (article) | 'gradient'
cfg.llavor     = 1;             % llavor aleatoria ([] = cada cop diferent)

% =====================================================================
%  DOMINI DE CERCA (graella 3D)
% =====================================================================
cfg.domini.Lx = 400;            % m
cfg.domini.Ly = 300;            % m
cfg.domini.Lz = 40;             % alcada maxima de vol (m); el terra es z = 0
cfg.domini.dx = 5;              % mida horitzontal de cel.la = pas del dron (m)
cfg.domini.dz = 5;              % mida vertical de cel.la (m)

% =====================================================================
%  FONT DE LA FUGA I DRON  ([x y z] en m)
% =====================================================================
cfg.font.pos          = [60 150 10];    % fuga a 10 m d'alcada (el dron NO ho sap)
cfg.dron.posInicial   = [340 215 30];   % el dron surt a 30 m d'alcada
cfg.dron.velocitat    = 5;              % m/s  -> dt = dx / velocitat
cfg.dron.maxPassos    = 2500;
cfg.dron.radiDeteccio = [];             % distancia horitzontal per trobar la fuga (m); [] = 1.5*dx
cfg.dron.toleranciaAlcada = 0;          % diferencia d'alcada admesa (m); 0 = ha d'arribar a la mateixa alcada

% =====================================================================
%  VENT MITJA
% =====================================================================
cfg.vent.direccio = 0;          % graus: direccio CAP ON bufa (0 = cap a +x)

% =====================================================================
%  ESCENARI ALEATORI (condicions inicials diferents a cada prova)
% =====================================================================
% Igual que a la versio 2D (vegeu aplicar_escenari_aleatori.m), mes les
% alcades de la fuga i del dron (aplicar_escenari_aleatori_3d.m).
cfg.atzar.activat         = false;
cfg.atzar.llavor          = [];        % llavor de l'escenari ([] = diferent cada cop)
cfg.atzar.mostrar         = true;
cfg.atzar.font            = true;      % posicio horitzontal de la fuga
cfg.atzar.alcadaFont      = true;      % alcada de la fuga
cfg.atzar.dron            = true;      % posicio horitzontal inicial del dron
cfg.atzar.alcadaDron      = true;      % alcada inicial del dron
cfg.atzar.ventDireccio    = true;
cfg.atzar.ventVelocitat   = true;
cfg.atzar.estabilitat     = false;
cfg.atzar.emissio         = false;
cfg.atzar.rangDireccio    = [0 360];   % graus
cfg.atzar.rangVelocitat   = [1.5 4.5]; % m/s
cfg.atzar.classes         = 'BCDE';
cfg.atzar.rangEmissio     = [0.5 2];
cfg.atzar.rangAlcadaFont  = [0 20];    % m
cfg.atzar.rangAlcadaDron  = [10 40];   % m
cfg.atzar.distanciaMinima = 150;       % distancia horitzontal minima dron-fuga (m)
cfg.atzar.marge           = 20;        % m
cfg.atzar.dronASotavent   = true;
cfg.atzar.angleMaxim      = 60;        % graus

% =====================================================================
%  SENSOR: k ~ Poisson( dt*(sensibilitat*c + fons) )
% =====================================================================
cfg.sensor.sensibilitat = 1;
cfg.sensor.fons         = 0.02;

% =====================================================================
%  PLOMA GAUSSIANA 3D (l'alcada de la font es font.pos(3))
% =====================================================================
cfg.gauss.Q           = 5000;   % taxa d'emissio (unitats/s)
cfg.gauss.u           = 3;      % velocitat del vent (m/s)
cfg.gauss.estabilitat = 'D';    % classe d'estabilitat A..F
cfg.gauss.sigma0      = 1;      % mida inicial de la ploma a la font (m)

% =====================================================================
%  PLOMA COMPLEXA 3D (puffs turbulents amb vent variable)
% =====================================================================
cfg.complexa.Q            = 2000;  % taxa d'emissio (unitats/s)
cfg.complexa.U            = 2;     % vent mitja a l'alcada de referencia (m/s)
cfg.complexa.zRef         = 10;    % alcada de referencia del vent (m)
cfg.complexa.exponentPerfil = 0.15; % perfil vertical U(z) = U (z/zRef)^p
cfg.complexa.puffsPerSeg  = 2;
cfg.complexa.sigma0       = 1.5;   % mida inicial d'un puff (m)
cfg.complexa.DpuffH       = 0.4;   % creixement horitzontal del puff (m^2/s)
cfg.complexa.DpuffV       = 0.2;   % creixement vertical del puff (m^2/s)
cfg.complexa.DturbH       = 0.4;   % difusio turbulenta horitzontal del centre (m^2/s)
cfg.complexa.DturbV       = 0.2;   % difusio turbulenta vertical del centre (m^2/s)
cfg.complexa.tau          = 200;   % temps de vida del contaminant (s)
cfg.complexa.sigmaTheta   = 15;    % meandreig de la direccio del vent (graus)
cfg.complexa.Ttheta       = 40;    % s
cfg.complexa.sigmaRafega  = 0.2;   % rafegues (fraccio de U)
cfg.complexa.Trafega      = 20;    % s
cfg.complexa.nRemolins    = 6;     % remolins horitzontals de gran escala (0 = cap)
cfg.complexa.Uremolins    = 0.6;   % m/s
cfg.complexa.lambdaMin    = 80;    % m
cfg.complexa.lambdaMax    = 250;   % m
cfg.complexa.Tremolins    = 60;    % s
cfg.complexa.tEscalfament = [];    % s; [] = auto

% =====================================================================
%  INFOTAXIS
% =====================================================================
cfg.infotaxis.saturacio   = Inf;   % 1 = sensor binari
cfg.infotaxis.binsComptes = [0 1 2 3 5 8 13 21 34 55 89];
cfg.infotaxis.Dmodel      = [];    % difusivitat horitzontal del model (m^2/s); [] = auto
cfg.infotaxis.DmodelV     = [];    % difusivitat vertical del model (m^2/s); [] = auto
cfg.infotaxis.factorModel = 1;

% =====================================================================
%  GRADIENT (7 sensors: centre, +-x, +-y, +-z)
% =====================================================================
cfg.gradient.separacio    = [];    % m; [] = dx
cfg.gradient.alfa         = 0.6;
cfg.gradient.memoria      = 8;
cfg.gradient.persistencia = 0.85;
cfg.gradient.probVertical = 0.15;  % en exploracio, prob. de triar un moviment vertical

% =====================================================================
%  VISUALITZACIO
% =====================================================================
cfg.visual.mostrar      = true;
cfg.visual.cadaNPassos  = 1;
cfg.visual.pausa        = 0.01;
cfg.visual.resolucio    = 2.5;       % resolucio del tall horitzontal (m)
cfg.visual.resolucio3d  = [5 2.5];   % resolucio [horitzontal vertical] de la ploma 3D (m)
cfg.visual.exageracioVertical = 4;   % l'eix z es veu 4 vegades mes gran
cfg.visual.nivellsIso  = [0.5 5];    % superficies de la ploma 3D (comptes esperats per pas)
cfg.visual.vista        = [-35 25];  % angle de la camera [azimut elevacio] (graus)
cfg.visual.gravarVideo  = false;     % desa un video de l'animacio
cfg.visual.formatVideo  = 'mp4';     % 'mp4' | 'avi' | 'gif'
cfg.visual.fpsVideo     = 15;
cfg.visual.carpetaVideo = '';        % '' = carpeta 'videos' dins la carpeta actual
cfg.visual.nomVideo     = '';        % '' = automatic (ploma, algorisme, data i hora)
end
