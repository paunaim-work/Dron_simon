classdef Ploma
%PLOMA  Genera la ploma de contaminant i la guarda en un fitxer .nc.
%
%   La ploma es una matriu de concentracions amb dimensio temporal:
%
%       camp.c(ix, iy, iz, it) = concentracio al punt (x(ix), y(iy), z(iz))
%                                a l'instant t(it)
%
%   "camp" es una struct amb els eixos x, y, z, t (vectors), la matriu c i
%   la posicio de la fuga posFuga. Es el mateix que hi ha dins del .nc.
%
%   Hi ha 4 versions, que surten de combinar dues opcions:
%
%     tipus     'gaussiana'  Estacionaria: no canvia amb el temps (un sol
%                            instant). Formula classica de Pasquill-Gifford.
%               'complexa'   Turbulenta: la fuga emet "puffs" (nuvolets de
%                            gas) que el vent s'emporta, s'eixamplen i es
%                            van esvaint.
%
%     dimensio  '3d'         Moltes alcades.
%               '2d'         Una sola alcada: la de la fuga.
%
%   Aquesta classe no crea objectes: nomes agrupa funcions (com una classe
%   de Java amb tot "static"). Exemple:
%
%       par  = Ploma.parametres();
%       camp = Ploma.generar('complexa', '3d', escenari, par);
%       Ploma.guardar(camp, 'ploma.nc');
%       camp = Ploma.llegir('ploma.nc');
%
%   Unitats: metres i segons. El vent bufa cap a +x i el terra es z = 0.

methods (Static)

    function par = parametres()
    %PARAMETRES  Valors per defecte dels parametres de la ploma.

        % Graella del fitxer (comuna a tots els tipus)
        par.resolucio  = 5;       % distancia entre punts de la graella (m)
        par.durada     = 1200;    % segons que es guarden (despres, es repeteix l'ultim instant)
        par.cadaSegons = 5;       % es guarda un instant cada tants segons
        par.tempsPrevi = 300;     % segons que fa que la fuga emet a l'instant t = 0
        par.pasTemps   = 1;       % pas de temps per moure els puffs (s)

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

    function camp = generar(tipus, dimensio, escenari, par)
    %GENERAR  Genera la ploma i la retorna com a matriu amb dimensio temporal.
    %   escenari ha de tenir midaX, midaY, midaZ i posFuga.
        f = escenari.posFuga;
        camp.x = 0:par.resolucio:escenari.midaX;
        camp.y = 0:par.resolucio:escenari.midaY;
        if strcmp(dimensio, '3d')
            camp.z = 0:par.resolucio:escenari.midaZ;
        else
            camp.z = f(3);                   % 2D: nomes l'alcada de la fuga
        end
        camp.posFuga = f;

        if strcmp(tipus, 'gaussiana')
            camp.t = 0;                      % no canvia amb el temps: un sol instant
            [X, Y, Z] = ndgrid(camp.x, camp.y, camp.z);
            camp.c = single(Ploma.formulaGaussiana(X - f(1), Y - f(2), Z, f(3), par.gaussiana));
        else
            camp.t = 0:par.cadaSegons:par.durada;
            camp.c = zeros(numel(camp.x), numel(camp.y), numel(camp.z), numel(camp.t), 'single');
            estat  = Ploma.puffsInicials();
            zona   = [escenari.midaX escenari.midaY];
            % Abans de t = 0 la fuga ja fa estona que emet
            for i = 1:round(par.tempsPrevi / par.pasTemps)
                estat = Ploma.avancarPuffs(estat, f, zona, par.complexa, par.pasTemps);
            end
            % Es guarda un instant cada par.cadaSegons
            for it = 1:numel(camp.t)
                if it > 1
                    for i = 1:round(par.cadaSegons / par.pasTemps)
                        estat = Ploma.avancarPuffs(estat, f, zona, par.complexa, par.pasTemps);
                    end
                end
                camp.c(:, :, :, it) = Ploma.sumaPuffs(estat, camp.x, camp.y, camp.z);
            end
        end
        camp.c(camp.c < 1e-4) = 0;           % el sensor no ho nota, i aixi el fitxer ocupa molt menys
    end

    function guardar(camp, fitxer)
    %GUARDAR  Guarda el camp en un fitxer .nc (NetCDF).
        if exist(fitxer, 'file')
            delete(fitxer);                  % nccreate no sobreescriu fitxers
        end
        n = size(camp.c);
        n(end + 1:4) = 1;                    % [nx ny nz nt], encara que nz o nt siguin 1
        nccreate(fitxer, 'x', 'Dimensions', {'x', n(1)});
        nccreate(fitxer, 'y', 'Dimensions', {'y', n(2)});
        nccreate(fitxer, 'z', 'Dimensions', {'z', n(3)});
        nccreate(fitxer, 't', 'Dimensions', {'t', n(4)});
        nccreate(fitxer, 'c', 'Dimensions', {'x', n(1), 'y', n(2), 'z', n(3), 't', n(4)}, ...
                 'Datatype', 'single', 'DeflateLevel', 5);   % comprimit
        ncwrite(fitxer, 'x', camp.x(:));         % (:) = com a columna
        ncwrite(fitxer, 'y', camp.y(:));
        ncwrite(fitxer, 'z', camp.z(:));
        ncwrite(fitxer, 't', camp.t(:));
        ncwrite(fitxer, 'c', camp.c);
        ncwriteatt(fitxer, 'x', 'units', 'm');
        ncwriteatt(fitxer, 'y', 'units', 'm');
        ncwriteatt(fitxer, 'z', 'units', 'm');
        ncwriteatt(fitxer, 't', 'units', 's');
        ncwriteatt(fitxer, '/', 'posFuga', camp.posFuga);
    end

    function camp = llegir(fitxer)
    %LLEGIR  Llegeix un camp d'un fitxer .nc (el que ha escrit guardar).
        camp.x = ncread(fitxer, 'x')';
        camp.y = ncread(fitxer, 'y')';
        camp.z = ncread(fitxer, 'z')';
        camp.t = ncread(fitxer, 't')';
        camp.c = ncread(fitxer, 'c');
        camp.posFuga = ncreadatt(fitxer, '/', 'posFuga');
        camp.posFuga = camp.posFuga(:)';
    end

end

methods (Static, Access = private)

    function c = formulaGaussiana(dx, dy, z, zFuga, p)
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
        c  = zeros(size(dx));
        hi = dx > -p.midaFuga;               % vent amunt de la fuga no hi arriba res
        d  = max(dx(hi), 0);
        sy = sqrt((0.08 * d ./ sqrt(1 + 1e-4   * d)).^2 + p.midaFuga^2);
        sz = sqrt((0.06 * d ./ sqrt(1 + 1.5e-3 * d)).^2 + p.midaFuga^2);
        vertical = exp(-(z(hi) - zFuga).^2 ./ (2 * sz.^2)) + exp(-(z(hi) + zFuga).^2 ./ (2 * sz.^2));
        c(hi) = p.emissio ./ (2 * pi * p.vent * sy .* sz) .* exp(-dy(hi).^2 ./ (2 * sy.^2)) .* vertical;
    end

    function estat = puffsInicials()
    %PUFFSINICIALS  Estat de la ploma complexa abans que la fuga emeti res:
    %   cap puff (un valor per puff a cada vector) i vent cap a +x.
        buit = zeros(0, 1);
        estat.x       = buit;                % posicio del centre de cada puff (m)
        estat.y       = buit;
        estat.z       = buit;
        estat.sigma2H = buit;                % mida^2 en horitzontal (m^2)
        estat.sigma2V = buit;                % mida^2 en vertical (m^2)
        estat.massa   = buit;                % quantitat de contaminant
        estat.angleVent = 0;                 % direccio del vent (rad); 0 = cap a +x
    end

    function e = avancarPuffs(e, posFuga, zona, p, dt)
    %AVANCARPUFFS  Fa passar dt segons a la ploma complexa.

        % 1. La direccio del vent oscil.la a l'atzar: recorda una mica la
        %    direccio d'abans i hi suma soroll (proces d'Ornstein-Uhlenbeck)
        memoria = exp(-dt / p.meandreTemps);
        e.angleVent = memoria * e.angleVent + (p.meandreGraus * pi / 180) * sqrt(1 - memoria^2) * randn;

        % 2. Els puffs es mouen amb el vent i a l'atzar (turbulencia),
        %    s'eixamplen i es van esvaint
        saltH = sqrt(2 * p.turbulenciaH * dt);
        saltV = sqrt(2 * p.turbulenciaV * dt);
        e.x = e.x + p.vent * cos(e.angleVent) * dt + saltH * randn(size(e.x));
        e.y = e.y + p.vent * sin(e.angleVent) * dt + saltH * randn(size(e.y));
        e.z = abs(e.z + saltV * randn(size(e.z)));       % abs: reboten al terra
        e.sigma2H = e.sigma2H + 2 * p.creixementH * dt;
        e.sigma2V = e.sigma2V + 2 * p.creixementV * dt;
        e.massa   = e.massa * exp(-dt / p.tempsVida);

        % 3. La fuga emet puffs nous (com a minim un a cada pas)
        n = max(1, round(p.puffsPerSegon * dt));
        u = ones(n, 1);
        e.x       = [e.x;       posFuga(1) * u];
        e.y       = [e.y;       posFuga(2) * u];
        e.z       = [e.z;       posFuga(3) * u];
        e.sigma2H = [e.sigma2H; p.midaPuff^2 * u];
        e.sigma2V = [e.sigma2V; p.midaPuff^2 * u];
        e.massa   = [e.massa;   p.emissio * dt / n * u];

        % 4. S'esborren els puffs que han sortit de la zona (amb 60 m de marge)
        fora = e.x < -60 | e.x > zona(1) + 60 | e.y < -60 | e.y > zona(2) + 60;
        e.x(fora) = [];        e.y(fora) = [];        e.z(fora) = [];
        e.sigma2H(fora) = [];  e.sigma2V(fora) = [];  e.massa(fora) = [];
    end

    function c = sumaPuffs(e, x, y, z)
    %SUMAPUFFS  Concentracio de tots els puffs a la graella x, y, z
    %   (resultat de mida [nx ny nz]). Cada puff es una campana de Gauss en
    %   3D, mes la seva imatge sota terra (el terra fa rebotar el gas).
    %   Truc: una campana 3D es el producte de tres campanes 1D (en x, en y
    %   i en z), aixi que es pot calcular amb productes de matrius.
        pes = e.massa ./ ((2 * pi)^1.5 * e.sigma2H .* sqrt(e.sigma2V));   % alcada de cada campana
        gx  = exp(-(x(:) - e.x').^2 ./ (2 * e.sigma2H'));                  % (nx x puffs)
        gy  = exp(-(y(:) - e.y').^2 ./ (2 * e.sigma2H'));                  % (ny x puffs)
        gz  = exp(-(z(:) - e.z').^2 ./ (2 * e.sigma2V')) + ...
              exp(-(z(:) + e.z').^2 ./ (2 * e.sigma2V'));                  % (nz x puffs)
        c = zeros(numel(x), numel(y), numel(z));
        for k = 1:numel(z)
            c(:, :, k) = gx * ((pes .* gz(k, :)') .* gy');
        end
    end

end
end
