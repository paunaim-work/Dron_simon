classdef Ploma
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
%   Aquest fitxer nomes agrupa funcions. Es fan servir aixi:
%
%       par   = Ploma.parametres();                         % valors per defecte
%       ploma = Ploma.crear('complexa', '3d', escenari, par);
%       ploma = Ploma.avancar(ploma, dt);                   % passen dt segons
%       c     = Ploma.concentracio(ploma, x, y, z);         % concentracio a (x, y, z)
%
%   La variable "ploma" es una struct amb tot el que cal saber de la ploma
%   en un instant (tipus, on es la fuga, puffs, direccio del vent...).
%
%   Unitats: metres i segons. El vent bufa cap a +x i el terra es z = 0.

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

    function ploma = crear(tipus, dimensio, escenari, par)
    %CREAR  Ploma nova, just quan la fuga comenca a emetre.
    %   tipus: 'gaussiana' o 'complexa';  dimensio: '2d' o '3d'.
    %   escenari ha de tenir posFuga ([x y z]), midaX i midaY.
        ploma.tipus     = tipus;
        ploma.dimensio  = dimensio;
        ploma.posFuga   = escenari.posFuga;
        ploma.zona      = [escenari.midaX escenari.midaY];   % fora d'aqui s'esborren els puffs
        ploma.par       = par.(tipus);       % nomes els parametres del tipus triat
        ploma.angleVent = 0;                 % direccio del vent (rad); 0 = cap a +x

        % Puffs (nomes la complexa): un valor per puff a cada vector
        buit = zeros(0, 1);
        ploma.puffs.x       = buit;          % posicio del centre (m)
        ploma.puffs.y       = buit;
        ploma.puffs.z       = buit;
        ploma.puffs.sigma2H = buit;          % mida^2 en horitzontal (m^2)
        ploma.puffs.sigma2V = buit;          % mida^2 en vertical (m^2)
        ploma.puffs.massa   = buit;          % quantitat de contaminant
    end

    function ploma = avancar(ploma, dt)
    %AVANCAR  Fa passar dt segons. La ploma gaussiana no canvia mai.
        if strcmp(ploma.tipus, 'gaussiana')
            return;
        end
        p  = ploma.par;
        pf = ploma.puffs;

        % 1. La direccio del vent oscil.la a l'atzar: recorda una mica la
        %    direccio d'abans i hi suma soroll (proces d'Ornstein-Uhlenbeck)
        memoria = exp(-dt / p.meandreTemps);
        ploma.angleVent = memoria * ploma.angleVent + ...
                          (p.meandreGraus * pi / 180) * sqrt(1 - memoria^2) * randn;

        % 2. Els puffs es mouen amb el vent i a l'atzar (turbulencia),
        %    s'eixamplen i es van esvaint
        saltH = sqrt(2 * p.turbulenciaH * dt);
        saltV = sqrt(2 * p.turbulenciaV * dt);
        pf.x = pf.x + p.vent * cos(ploma.angleVent) * dt + saltH * randn(size(pf.x));
        pf.y = pf.y + p.vent * sin(ploma.angleVent) * dt + saltH * randn(size(pf.y));
        pf.z = abs(pf.z + saltV * randn(size(pf.z)));     % abs: reboten al terra
        pf.sigma2H = pf.sigma2H + 2 * p.creixementH * dt;
        pf.sigma2V = pf.sigma2V + 2 * p.creixementV * dt;
        pf.massa   = pf.massa * exp(-dt / p.tempsVida);

        % 3. La fuga emet puffs nous (com a minim un a cada pas)
        n = max(1, round(p.puffsPerSegon * dt));
        u = ones(n, 1);
        pf.x       = [pf.x;       ploma.posFuga(1) * u];
        pf.y       = [pf.y;       ploma.posFuga(2) * u];
        pf.z       = [pf.z;       ploma.posFuga(3) * u];
        pf.sigma2H = [pf.sigma2H; p.midaPuff^2 * u];
        pf.sigma2V = [pf.sigma2V; p.midaPuff^2 * u];
        pf.massa   = [pf.massa;   p.emissio * dt / n * u];

        % 4. S'esborren els puffs que han sortit de la zona (amb 60 m de marge)
        fora = pf.x < -60 | pf.x > ploma.zona(1) + 60 | pf.y < -60 | pf.y > ploma.zona(2) + 60;
        pf.x(fora) = [];        pf.y(fora) = [];        pf.z(fora) = [];
        pf.sigma2H(fora) = [];  pf.sigma2V(fora) = [];  pf.massa(fora) = [];

        ploma.puffs = pf;
    end

    function c = concentracio(ploma, x, y, z)
    %CONCENTRACIO  Concentracio real de la ploma als punts (x, y, z).
    %   x, y i z poden ser numeros, vectors o matrius (de la mateixa mida).
        zero = zeros(size(x + y + z));       % perque x, y i z tinguin la mateixa mida
        x = x + zero;
        y = y + zero;
        z = z + zero;
        if strcmp(ploma.dimensio, '2d')
            z(:) = ploma.posFuga(3);         % 2D: sempre al pla de la fuga
        end

        if strcmp(ploma.tipus, 'gaussiana')
            c = Ploma.formula_gaussiana(x - ploma.posFuga(1), y - ploma.posFuga(2), ...
                                        z, ploma.posFuga(3), ploma.par);
        else
            c = Ploma.suma_puffs(ploma.puffs, x, y, z);
        end
    end

    function c = mitjana(ploma, dx, dy, zDron, zFuga)
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
        if strcmp(ploma.tipus, 'gaussiana')
            c = Ploma.formula_gaussiana(dx, dy, zDron, zFuga, ploma.par);
            return;
        end
        p  = ploma.par;
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

methods (Static, Access = private)

    function c = formula_gaussiana(dx, dy, z, zFuga, p)
    %FORMULA_GAUSSIANA  Ploma gaussiana a dx metres vent avall de la fuga,
    %   dy de traves al vent i a l'alcada z, amb la fuga a l'alcada H = zFuga:
    %
    %                Q              dy^2             (z-H)^2            (z+H)^2
    %       c = ------------ exp(- ------) [ exp(- -------) + exp(- -------) ]
    %           2 pi U sy sz       2 sy^2           2 sz^2             2 sz^2
    %
    %   sy i sz (amplada i gruix de la ploma) creixen amb la distancia
    %   (formules de Briggs, 1973, atmosfera neutra). L'ultim terme es el
    %   contaminant que rebota al terra.
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

    function c = suma_puffs(pf, x, y, z)
    %SUMA_PUFFS  Suma de tots els puffs. Cada puff es una campana de Gauss
    %   en 3D, mes la seva imatge sota terra (el terra fa rebotar el gas).
        pes      = pf.massa ./ ((2 * pi)^1.5 * pf.sigma2H .* sqrt(pf.sigma2V));   % alcada de cada campana
        dist2    = (x(:) - pf.x').^2 + (y(:) - pf.y').^2;                        % (punts x puffs)
        vertical = exp(-(z(:) - pf.z').^2 ./ (2 * pf.sigma2V')) + ...
                   exp(-(z(:) + pf.z').^2 ./ (2 * pf.sigma2V'));
        c = reshape((exp(-dist2 ./ (2 * pf.sigma2H')) .* vertical) * pes, size(x));
    end

end
end
