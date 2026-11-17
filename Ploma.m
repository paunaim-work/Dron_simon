classdef Ploma < handle
%PLOMA  La ploma de contaminant que surt de la fuga (4 versions).
%
%   Hi ha 4 versions, que surten de combinar dues opcions:
%
%     tipus     'gaussiana'  Estacionaria: no canvia amb el temps. Es la
%                            formula classica de Pasquill-Gifford.
%               'complexa'   Turbulenta: la fuga emet "puffs" (nuvolets de
%                            gas) que el vent s'emporta, s'eixamplen i es
%                            van esvaint.
%
%     dimensio  '3d'         La concentracio depen de x, y i z.
%               '2d'         Nomes es mira el pla horitzontal a l'alcada
%                            de la fuga (la z s'ignora).
%
%   Exemple:
%       par   = Ploma.parametres();                     % valors per defecte
%       ploma = Ploma('complexa', '3d', escenari, par); % constructor
%       ploma.avancar(dt);                              % passen dt segons
%       c = ploma.concentracio(x, y, z);                % concentracio a (x, y, z)
%
%   Com a Java: "< handle" fa que els objectes es passin per referencia, i
%   el primer argument dels metodes (obj) fa el paper de "this".
%
%   Unitats: metres i segons. El vent bufa cap a +x i el terra es z = 0.

properties (SetAccess = private)     % es poden llegir des de fora pero no canviar
    tipus                % 'gaussiana' o 'complexa'
    dimensio             % '2d' o '3d'
    posFuga              % [x y z] de la fuga (m)
    zona                 % [midaX midaY]: fora d'aqui s'esborren els puffs (m)
    par                  % parametres del tipus triat
    angleVent = 0        % direccio del vent (rad); 0 = cap a +x
    puffs                % puffs de la ploma complexa (struct: un valor per puff a cada camp)
end

methods

    function obj = Ploma(tipus, dimensio, escenari, par)
    %PLOMA  Constructor: ploma nova, just quan la fuga comenca a emetre.
    %   escenari ha de tenir posFuga, midaX i midaY.
        obj.tipus    = tipus;
        obj.dimensio = dimensio;
        obj.posFuga  = escenari.posFuga;
        obj.zona     = [escenari.midaX escenari.midaY];
        obj.par      = par.(tipus);          % nomes els parametres del tipus triat

        buit = zeros(0, 1);                  % encara no hi ha cap puff
        obj.puffs.x       = buit;            % posicio del centre (m)
        obj.puffs.y       = buit;
        obj.puffs.z       = buit;
        obj.puffs.sigma2H = buit;            % mida^2 en horitzontal (m^2)
        obj.puffs.sigma2V = buit;            % mida^2 en vertical (m^2)
        obj.puffs.massa   = buit;            % quantitat de contaminant
    end

    function avancar(obj, dt)
    %AVANCAR  Fa passar dt segons. La ploma gaussiana no canvia mai.
        if strcmp(obj.tipus, 'gaussiana')
            return;
        end
        p  = obj.par;
        pf = obj.puffs;

        % 1. La direccio del vent oscil.la a l'atzar: recorda una mica la
        %    direccio d'abans i hi suma soroll (proces d'Ornstein-Uhlenbeck)
        memoria = exp(-dt / p.meandreTemps);
        obj.angleVent = memoria * obj.angleVent + ...
                        (p.meandreGraus * pi / 180) * sqrt(1 - memoria^2) * randn;

        % 2. Els puffs es mouen amb el vent i a l'atzar (turbulencia),
        %    s'eixamplen i es van esvaint
        saltH = sqrt(2 * p.turbulenciaH * dt);
        saltV = sqrt(2 * p.turbulenciaV * dt);
        pf.x = pf.x + p.vent * cos(obj.angleVent) * dt + saltH * randn(size(pf.x));
        pf.y = pf.y + p.vent * sin(obj.angleVent) * dt + saltH * randn(size(pf.y));
        pf.z = abs(pf.z + saltV * randn(size(pf.z)));     % abs: reboten al terra
        pf.sigma2H = pf.sigma2H + 2 * p.creixementH * dt;
        pf.sigma2V = pf.sigma2V + 2 * p.creixementV * dt;
        pf.massa   = pf.massa * exp(-dt / p.tempsVida);

        % 3. La fuga emet puffs nous (com a minim un a cada pas)
        n = max(1, round(p.puffsPerSegon * dt));
        u = ones(n, 1);
        pf.x       = [pf.x;       obj.posFuga(1) * u];
        pf.y       = [pf.y;       obj.posFuga(2) * u];
        pf.z       = [pf.z;       obj.posFuga(3) * u];
        pf.sigma2H = [pf.sigma2H; p.midaPuff^2 * u];
        pf.sigma2V = [pf.sigma2V; p.midaPuff^2 * u];
        pf.massa   = [pf.massa;   p.emissio * dt / n * u];

        % 4. S'esborren els puffs que han sortit de la zona (amb 60 m de marge)
        fora = pf.x < -60 | pf.x > obj.zona(1) + 60 | pf.y < -60 | pf.y > obj.zona(2) + 60;
        pf.x(fora) = [];        pf.y(fora) = [];        pf.z(fora) = [];
        pf.sigma2H(fora) = [];  pf.sigma2V(fora) = [];  pf.massa(fora) = [];

        obj.puffs = pf;
    end

    function c = concentracio(obj, x, y, z)
    %CONCENTRACIO  Concentracio real de la ploma als punts (x, y, z).
    %   x, y i z poden ser numeros, vectors o matrius (de la mateixa mida).
        zero = zeros(size(x + y + z));       % perque x, y i z tinguin la mateixa mida
        x = x + zero;
        y = y + zero;
        z = z + zero;
        if strcmp(obj.dimensio, '2d')
            z(:) = obj.posFuga(3);           % 2D: sempre al pla de la fuga
        end

        if strcmp(obj.tipus, 'gaussiana')
            c = obj.formulaGaussiana(x - obj.posFuga(1), y - obj.posFuga(2), z, obj.posFuga(3));
        else
            c = obj.sumaPuffs(x, y, z);
        end
    end

    function c = mitjana(obj, dx, dy, zDron, zFuga)
    %MITJANA  Concentracio MITJANA quan el dron es a (dx, dy) de la fuga, ell
    %   a l'alcada zDron i la fuga a zFuga. L'infotaxis la fa servir com a
    %   model: nomes depen del tipus de ploma i dels seus parametres, no
    %   d'on es la fuga de veritat ni de la turbulencia d'aquell moment.
    %
    %   - Gaussiana: la mateixa formula (la ploma no canvia amb el temps).
    %   - Complexa: solucio mitjana de l'equacio d'adveccio-difusio
    %     (Vergassola et al., 2007). La turbulencia i el meandre del vent
    %     es resumeixen en dues difusivitats (Dh en horitzontal i Dv en
    %     vertical), i el terra es una fuga "imatge" sota terra:
    %
    %       c = Q / (4 pi sqrt(Dh Dv)) exp(dx U / 2Dh) [e^(-r1/abast)/r1 + e^(-r2/abast)/r2]
    %
    %     r1 i r2 son les distancies a la fuga i a la seva imatge (amb la
    %     z escalada per sqrt(Dh/Dv)), i abast es fins on arriba la ploma.
        if strcmp(obj.tipus, 'gaussiana')
            c = obj.formulaGaussiana(dx, dy, zDron, zFuga);
            return;
        end
        p  = obj.par;
        Dh = p.turbulenciaH + p.creixementH + (p.vent * p.meandreGraus * pi / 180)^2 * p.meandreTemps;
        Dv = p.turbulenciaV + p.creixementV;
        abast = sqrt(Dh * p.tempsVida / (1 + p.vent^2 * p.tempsVida / (4 * Dh)));
        r1 = max(sqrt(dx.^2 + dy.^2 + Dh / Dv * (zDron - zFuga)^2), p.midaPuff);   % mai mes a prop
        r2 = max(sqrt(dx.^2 + dy.^2 + Dh / Dv * (zDron + zFuga)^2), p.midaPuff);   % que un puff
        c = p.emissio / (4 * pi * sqrt(Dh * Dv)) * ...
            (exp(dx * p.vent / (2 * Dh) - r1 / abast) ./ r1 + ...
             exp(dx * p.vent / (2 * Dh) - r2 / abast) ./ r2);
    end

end

methods (Static)

    function par = parametres()
    %PARAMETRES  Valors per defecte dels parametres de cada tipus de ploma.

        % Ploma gaussiana
        par.gaussiana.emissio  = 5000;    % taxa d'emissio de la fuga
        par.gaussiana.vent     = 3;       % velocitat del vent (m/s)
        par.gaussiana.midaFuga = 1;       % mida de la fuga (m)

        % Ploma complexa
        par.complexa.emissio       = 2000;  % taxa d'emissio de la fuga
        par.complexa.vent          = 2;     % velocitat mitjana del vent (m/s)
        par.complexa.tempsVida     = 200;   % temps de vida del contaminant (s)
        par.complexa.puffsPerSegon = 2;     % puffs emesos per segon (com a minim 1 a cada pas)
        par.complexa.midaPuff      = 1.5;   % mida inicial d'un puff (m)
        par.complexa.creixementH   = 0.4;   % com de rapid s'eixampla un puff en horitzontal (m^2/s)
        par.complexa.creixementV   = 0.2;   % ... i en vertical (m^2/s)
        par.complexa.turbulenciaH  = 0.4;   % moviment aleatori dels puffs en horitzontal (m^2/s)
        par.complexa.turbulenciaV  = 0.2;   % ... i en vertical (m^2/s)
        par.complexa.meandreGraus  = 15;    % quant oscil.la la direccio del vent (graus)
        par.complexa.meandreTemps  = 40;    % temps tipic d'aquesta oscil.lacio (s)
    end

end

methods (Access = private)

    function c = formulaGaussiana(obj, dx, dy, z, zFuga)
    %FORMULAGAUSSIANA  Ploma gaussiana a dx metres vent avall de la fuga,
    %   dy de traves al vent i a l'alcada z, amb la fuga a l'alcada H = zFuga:
    %
    %                Q              dy^2             (z-H)^2            (z+H)^2
    %       c = ------------ exp(- ------) [ exp(- -------) + exp(- -------) ]
    %           2 pi U sy sz       2 sy^2           2 sz^2             2 sz^2
    %
    %   sy i sz (amplada i gruix de la ploma) creixen amb la distancia
    %   (formules de Briggs, 1973, atmosfera neutra). L'ultim terme es el
    %   contaminant que rebota al terra.
        p  = obj.par;
        dy = dy + zeros(size(dx));
        z  = z  + zeros(size(dx));
        c  = zeros(size(dx));
        hi = dx > -p.midaFuga;               % vent amunt de la fuga no hi arriba res
        d  = max(dx(hi), 0);
        sy = sqrt((0.08 * d ./ sqrt(1 + 1e-4   * d)).^2 + p.midaFuga^2);
        sz = sqrt((0.06 * d ./ sqrt(1 + 1.5e-3 * d)).^2 + p.midaFuga^2);
        vertical = exp(-(z(hi) - zFuga).^2 ./ (2 * sz.^2)) + exp(-(z(hi) + zFuga).^2 ./ (2 * sz.^2));
        c(hi) = p.emissio ./ (2 * pi * p.vent * sy .* sz) .* exp(-dy(hi).^2 ./ (2 * sy.^2)) .* vertical;
    end

    function c = sumaPuffs(obj, x, y, z)
    %SUMAPUFFS  Suma de tots els puffs. Cada puff es una campana de Gauss
    %   en 3D, mes la seva imatge sota terra (el terra fa rebotar el gas).
        pf       = obj.puffs;
        pes      = pf.massa ./ ((2 * pi)^1.5 * pf.sigma2H .* sqrt(pf.sigma2V));   % alcada de cada campana
        dist2    = (x(:) - pf.x').^2 + (y(:) - pf.y').^2;                        % (punts x puffs)
        vertical = exp(-(z(:) - pf.z').^2 ./ (2 * pf.sigma2V')) + ...
                   exp(-(z(:) + pf.z').^2 ./ (2 * pf.sigma2V'));
        c = reshape((exp(-dist2 ./ (2 * pf.sigma2H')) .* vertical) * pes, size(x));
    end

end
end
