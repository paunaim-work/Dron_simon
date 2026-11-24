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
%   I a mes hi ha la ploma REAL: una simulacio feta amb el model PALM (la
%   dispersio de PM10 al port de Tarragona), que ja ve en un fitxer .nc.
%   No es genera: importarPalm la llegeix i la posa en el mateix format.
%
%   Aquesta classe no crea objectes: nomes agrupa funcions (com una classe
%   de Java amb tot "static"). Exemple:
%
%       par  = Ploma.parametres();
%       camp = Ploma.generar('complexa', '3d', escenari, par);   % o be:
%       camp = Ploma.importarPalm('3d', par);                    % la real
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

        % Ploma real (simulacio PALM). Posicions en metres, en les coordenades del fitxer.
        par.real.fitxer       = 'tgn_port_dispersion_3d_N03_PM10_t0.nc';
        par.real.variable     = 'kc_PM10';           % concentracio de PM10 (kg/m^3)
        par.real.zona         = [237 538 0 500];     % tros que es llegeix: [xmin xmax ymin ymax]
        par.real.alcadaMaxima = 50;                  % ... fins a aquesta alcada
        par.real.posFuga      = [387.5 27.5 11.5];   % on es la fuga (el maxim de concentracio)
        par.real.direccioVent = 90;                  % cap a on bufa el vent (graus: 0 = +x, 90 = +y)
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

    function camp = importarPalm(dimensio, par)
    %IMPORTARPALM  Llegeix la ploma real (el .nc de PALM) i la deixa com les
    %   simulades (la struct camp), perque el dron la pugui fer servir igual:
    %
    %     1. Retalla la zona de cerca. El fitxer sencer son 750 x 750 x 202
    %        punts (1 m de resolucio): nomes se'n llegeix un tros.
    %     2. Fa la mitjana de cada punt amb els del seu voltant (un cub de
    %        5 m): el sensor no mesura en un punt, sino al llarg dels 5 m que
    %        el dron recorre en un pas.
    %     3. La gira perque el vent bufi cap a +x, que es el que suposa el
    %        dron. Nomes quarts de volta (rot90), aixi no cal interpolar.
    %     4. Canvia les unitats de la concentracio (kg/m^3) perque el sensor
    %        hi compti el mateix que a la ploma gaussiana (mateix flux Q/U).
    %
    %   Dins dels tancs i edificis el fitxer no te valor (NaN): hi posem c = 0.
        r = par.real;
        x = ncread(r.fitxer, 'x')';
        y = ncread(r.fitxer, 'y')';
        z = ncread(r.fitxer, 'zu_3d')';

        % 1. Retall: es llegeix nomes la zona i les alcades que volem
        ix = find(x >= r.zona(1) & x <= r.zona(2));
        iy = find(y >= r.zona(3) & y <= r.zona(4));
        iz = find(z > 0 & z <= r.alcadaMaxima);    % z = 0 es sota terra (i aixi el pas es 1 m)
        c  = ncread(r.fitxer, r.variable, [ix(1) iy(1) iz(1) 1], [numel(ix) numel(iy) numel(iz) 1]);
        c  = double(c);
        c(isnan(c) | c < 0) = 0;             % buits (tancs, edificis, sota terra) i soroll numeric

        % 2. Mitjana en un cub de 5 x 5 x 5 punts (convolucio amb un nucli constant)
        c = convn(c, ones(5, 5, 5) / 125, 'same');

        % 3. Gir de k quarts de volta (k > 0: sentit antihorari) fins que el vent va cap a +x.
        %    rot90 gira el pla (x, y) de totes les alcades alhora.
        k = -round(r.direccioVent / 90);
        c = rot90(c, k);
        % La fuga es gira igual: marquem la seva cel.la en una matriu i girem la matriu
        [~, ixFuga] = min(abs(x(ix) - r.posFuga(1)));
        [~, iyFuga] = min(abs(y(iy) - r.posFuga(2)));
        marca = false(numel(ix), numel(iy));
        marca(ixFuga, iyFuga) = true;
        [ixFuga, iyFuga] = find(rot90(marca, k));
        h = x(2) - x(1);                     % resolucio horitzontal (1 m)
        camp.x = (0:size(c, 1) - 1) * h;     % eixos nous, que comencen a 0
        camp.y = (0:size(c, 2) - 1) * h;
        camp.z = z(iz);
        camp.t = 0;                          % un sol instant: la ploma real esta "congelada"
        camp.posFuga = [camp.x(ixFuga) camp.y(iyFuga) r.posFuga(3)];

        % 4. Unitats. Si sumem la concentracio de tot un pla perpendicular al
        %    vent (x constant), I(x) = sum(c) * dy * dz, obtenim la massa que
        %    travessa el pla cada segon dividida per la velocitat del vent. En
        %    una ploma gaussiana val exactament Q/U, a qualsevol distancia.
        %    Multipliquem la ploma real perque la mediana de I(x) vent avall de
        %    la fuga sigui la Q/U de la gaussiana (la que el dron fa servir de model).
        dz = camp.z(2) - camp.z(1);          % pas vertical (1 m)
        I  = squeeze(sum(sum(c, 2), 3)) * h * dz;
        g  = par.gaussiana;
        c  = c * (g.emissio / g.vent) / median(I(camp.x > camp.posFuga(1)));

        if strcmp(dimensio, '2d')
            [~, nivell] = min(abs(camp.z - camp.posFuga(3)));   % 2D: nomes l'alcada de la fuga
            c = c(:, :, nivell);
            camp.z = camp.z(nivell);
        end
        camp.c = single(c);
        camp.c(camp.c < 1e-4) = 0;           % com a les plomes simulades
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
