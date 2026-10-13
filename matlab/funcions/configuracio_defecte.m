function cfg = configuracio_defecte()
%CONFIGURACIO_DEFECTE  Parametres per defecte de la simulacio de cerca.
%
%   cfg = configuracio_defecte() retorna una estructura amb tots els
%   parametres del problema: domini, font, dron, vent, sensor, models de
%   ploma i algorismes de cerca. Modifica els camps que vulguis abans de
%   cridar simular_cerca(cfg).
%
%   Unitats: metres (m), segons (s). La concentracio esta en unitats
%   arbitraries (p. ex. g/m^3 o Bq/m^3); nomes importa a traves de la
%   sensibilitat del sensor.

% =====================================================================
%  OPCIONS PRINCIPALS
% =====================================================================
cfg.tipusPloma = 'gaussiana';   % 'gaussiana' | 'complexa'
cfg.algorisme  = 'infotaxis';   % 'infotaxis' (article) | 'gradient'
cfg.llavor     = 1;             % llavor aleatoria ([] = cada cop diferent)

% =====================================================================
%  DOMINI DE CERCA (graella regular)
% =====================================================================
cfg.domini.Lx = 400;            % amplada (m)
cfg.domini.Ly = 300;            % alcada del mapa (m)
cfg.domini.dx = 5;              % mida de cel.la = pas del dron (m)

% =====================================================================
%  FONT DE LA FUGA I DRON
% =====================================================================
cfg.font.pos          = [60 150];   % posicio real de la fuga (m) (el dron NO la coneix)
cfg.dron.posInicial   = [340 215];  % posicio inicial del dron (m)
cfg.dron.velocitat    = 5;          % m/s  -> temps per pas dt = dx / velocitat
cfg.dron.maxPassos    = 2500;       % limit de passos (autonomia de la bateria)
cfg.dron.radiDeteccio = [];         % la fuga es troba a menys d'aquesta distancia (m); [] = 1.5*dx

% =====================================================================
%  VENT MITJA (comu a les dues plomes)
% =====================================================================
cfg.vent.direccio = 0;          % graus: direccio CAP ON bufa (0 = cap a +x)

% =====================================================================
%  ESCENARI ALEATORI (condicions inicials diferents a cada prova)
% =====================================================================
% Si activat = true, abans de cada simulacio se sortegen els parametres
% marcats amb true (vegeu aplicar_escenari_aleatori.m). Si es false, es
% fan servir els valors fixos d'aquest fitxer.
cfg.atzar.activat         = false;
cfg.atzar.llavor          = [];        % llavor de l'escenari ([] = diferent cada cop)
cfg.atzar.mostrar         = true;      % escriu l'escenari sortejat a la consola
cfg.atzar.font            = true;      % posicio de la fuga
cfg.atzar.dron            = true;      % posicio inicial del dron
cfg.atzar.ventDireccio    = true;      % direccio del vent
cfg.atzar.ventVelocitat   = true;      % velocitat del vent
cfg.atzar.estabilitat     = false;     % classe d'estabilitat (ploma gaussiana)
cfg.atzar.emissio         = false;     % taxa d'emissio Q
cfg.atzar.rangDireccio    = [0 360];   % graus
cfg.atzar.rangVelocitat   = [1.5 4.5]; % m/s
cfg.atzar.classes         = 'BCDE';    % classes d'estabilitat possibles
cfg.atzar.rangEmissio     = [0.5 2];   % factor que multiplica Q
cfg.atzar.distanciaMinima = 150;       % distancia minima entre el dron i la fuga (m)
cfg.atzar.marge           = 20;        % distancia minima a les vores del domini (m)
cfg.atzar.dronASotavent   = true;      % el dron surt vent avall de la fuga...
cfg.atzar.angleMaxim      = 60;        % ...dins d'un con de +-angleMaxim graus

% =====================================================================
%  SENSOR (comptador de deteccions, tipus Geiger / detector de particules)
%  Nombre de comptes en un pas: k ~ Poisson( dt*(sensibilitat*c + fons) )
% =====================================================================
cfg.sensor.sensibilitat = 1;    % comptes/s per unitat de concentracio
cfg.sensor.fons         = 0.02; % comptes/s de fons (falses deteccions)

% =====================================================================
%  PLOMA GAUSSIANA (model estacionari classic de Pasquill-Gifford)
% =====================================================================
cfg.gauss.Q           = 5000;   % taxa d'emissio (unitats/s)
cfg.gauss.u           = 3;      % velocitat del vent (m/s)
cfg.gauss.H           = 10;     % alcada efectiva de la font (m)
cfg.gauss.z           = 10;     % alcada de vol del dron (m)
cfg.gauss.estabilitat = 'D';    % classe d'estabilitat atmosferica A..F
cfg.gauss.sigma0      = 1;      % mida inicial de la ploma a la font (m)

% =====================================================================
%  PLOMA COMPLEXA (model de "puffs" turbulents amb vent variable)
% =====================================================================
cfg.complexa.Q            = 200;   % taxa d'emissio (unitats/s)
cfg.complexa.U            = 2;     % velocitat mitjana del vent (m/s)
cfg.complexa.puffsPerSeg  = 2;     % puffs emesos per segon
cfg.complexa.sigma0       = 1.5;   % mida inicial d'un puff (m)
cfg.complexa.Dpuff        = 0.4;   % creixement del puff (difusio relativa, m^2/s)
cfg.complexa.Dturb        = 0.4;   % difusio turbulenta del centre del puff (m^2/s)
cfg.complexa.tau          = 200;   % temps de vida del contaminant (s)
% Meandreig: la direccio del vent oscil.la (proces d'Ornstein-Uhlenbeck)
cfg.complexa.sigmaTheta   = 15;    % desviacio de la direccio (graus)
cfg.complexa.Ttheta       = 40;    % temps de correlacio (s)
% Rafegues: la intensitat del vent oscil.la
cfg.complexa.sigmaRafega  = 0.2;   % fraccio de U
cfg.complexa.Trafega      = 20;    % s
% Remolins de gran escala (camp de velocitat de divergencia nul.la)
cfg.complexa.nRemolins    = 6;     % nombre de modes (0 = sense remolins)
cfg.complexa.Uremolins    = 0.6;   % velocitat tipica dels remolins (m/s)
cfg.complexa.lambdaMin    = 80;    % longitud d'ona minima (m)
cfg.complexa.lambdaMax    = 250;   % longitud d'ona maxima (m)
cfg.complexa.Tremolins    = 60;    % temps de correlacio (s)
cfg.complexa.tEscalfament = [];    % temps previ per desenvolupar la ploma (s); [] = auto

% =====================================================================
%  ALGORISME INFOTAXIS (Vergassola, Villermaux & Shraiman, Nature 2007)
% =====================================================================
% El dron necessita un MODEL intern de la taxa mitjana de deteccions
% R(r|r0). Per a la ploma gaussiana s'usa la mateixa formula gaussiana;
% per a la complexa, la solucio d'adveccio-difusio de l'article (eq. 5-7
% del material suplementari) amb una difusivitat efectiva.
cfg.infotaxis.saturacio   = Inf;   % comptes >= saturacio s'agrupen (1 = sensor binari si/no)
cfg.infotaxis.binsComptes = [0 1 2 3 5 8 13 21 34 55 89]; % intervals de comptes per planificar
cfg.infotaxis.Dmodel      = [];    % difusivitat efectiva del model (m^2/s); [] = auto
cfg.infotaxis.factorModel = 1;     % multiplicador de la taxa del model (error de calibratge)

% =====================================================================
%  ALGORISME DE GRADIENT (quimiotaxi)
% =====================================================================
% El dron porta 4 sensors en creu (a +-separacio en x i en y), estima el
% gradient de concentracio per diferencies finites i es mou cap on
% augmenta. Si no detecta res, fa un passeig aleatori correlacionat.
cfg.gradient.separacio    = [];    % distancia dels sensors al centre (m); [] = dx
cfg.gradient.alfa         = 0.6;   % filtre exponencial del gradient (0 = sense filtre)
cfg.gradient.memoria      = 8;     % passos sense deteccio abans de passar a exploracio
cfg.gradient.persistencia = 0.85;  % probabilitat de mantenir la direccio en exploracio

% =====================================================================
%  MODE REALISTA (limitacions tecniques del dron, fitxa tecnica)
% =====================================================================
% Tots els parametres i el seu significat: configuracio_realista.m
cfg.realista = configuracio_realista();

% =====================================================================
%  VISUALITZACIO
% =====================================================================
cfg.visual.mostrar      = true;    % animacio en directe
cfg.visual.cadaNPassos  = 1;       % redibuixa cada N passos (puja'l si va lent)
cfg.visual.pausa        = 0.01;    % pausa entre fotogrames (s)
cfg.visual.resolucio    = 2.5;     % resolucio del mapa de la ploma (m)
cfg.visual.gravarVideo  = false;   % desa un video de l'animacio
cfg.visual.formatVideo  = 'mp4';   % 'mp4' | 'avi' | 'gif'
cfg.visual.fpsVideo     = 15;      % fotogrames per segon del video
cfg.visual.carpetaVideo = '';      % carpeta on es desa ('' = carpeta 'videos' dins la carpeta actual)
cfg.visual.nomVideo     = '';      % nom del fitxer ('' = automatic, amb ploma, algorisme, data i hora)
end
