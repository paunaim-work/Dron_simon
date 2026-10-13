function [K, info] = model_taxa_deteccio_3d(cfg, dom)
%MODEL_TAXA_DETECCIO_3D  Model intern del dron per a l'infotaxis 3D:
%   comptes esperats en un pas, mu(r | r0) = dt * R(r | r0).
%
%   En horitzontal mu nomes depen del desplacament (x-x0, y-y0), pero en
%   vertical depen de l'alcada del dron z i de la font z0 per separat (el
%   terra reflecteix el contaminant). Per aixo K es una taula 4D:
%
%       K(iy - iy0 + Ny, ix - ix0 + Nx, iz, iz0)
%
%   - Ploma gaussiana: formula gaussiana amb reflexio (model perfecte).
%   - Ploma complexa: solucio 3D de l'equacio d'adveccio-difusio de
%     l'article (material suplementari, eq. 6) amb difusivitat horitzontal
%     Dh i vertical Dv, i una font imatge sota terra (reflexio):
%
%        c = Q / (4 pi sqrt(Dh Dv)) exp(s V / 2Dh) [e^(-r1/lambda)/r1 + e^(-r2/lambda)/r2]
%
%     on r1 i r2 son les distancies (amb z escalada per sqrt(Dh/Dv)) a la
%     font i a la seva imatge, i lambda = sqrt(Dh tau / (1 + V^2 tau/(4 Dh))).

[DX, DY] = meshgrid((-(dom.Nx - 1):(dom.Nx - 1)) * dom.dx, ...
                    (-(dom.Ny - 1):(dom.Ny - 1)) * dom.dx);
K = zeros(2 * dom.Ny - 1, 2 * dom.Nx - 1, dom.Nz, dom.Nz);
info = struct();

theta = angle_vent_model(cfg, dom);     % en mode realista, amb l'error de l'anemometre
flot  = flotacio_cfg(cfg);              % densitat del gas (el dron sap quin gas busca)

switch lower(cfg.tipusPloma)
    case 'gaussiana'
        g = cfg.gauss;
        p.Q = g.Q; p.u = g.u; p.estabilitat = g.estabilitat;
        p.sigma0 = g.sigma0; p.theta = theta; p.flot = flot;
        for iz = 1:dom.Nz
            for iz0 = 1:dom.Nz
                c = ploma_gaussiana_conc_3d(DX, DY, dom.z(iz), dom.z(iz0), p);
                K(:, :, iz, iz0) = c;
            end
        end
        info.descripcio = 'gaussiana 3D (model exacte)';

    case 'complexa'
        cc = cfg.complexa;
        Dh = cfg.infotaxis.Dmodel;
        if isempty(Dh)
            Dh = cc.DturbH + cc.DpuffH + (cc.U * cc.sigmaTheta * pi / 180)^2 * cc.Ttheta;
            if cc.nRemolins > 0
                Dh = Dh + 0.5 * cc.Uremolins^2 * cc.Tremolins;
            end
        end
        Dv = cfg.infotaxis.DmodelV;
        if isempty(Dv)
            Dv = cc.DturbV + cc.DpuffV;
        end
        V   = cc.U;
        tau = cc.tau;
        lam = sqrt(Dh * tau / (1 + V^2 * tau / (4 * Dh)));
        a   = sqrt(Dh / Dv);                      % escalat vertical
        s   = DX * cos(theta) + DY * sin(theta);
        rh2 = DX.^2 + DY.^2;
        rmin = dom.dx / 2;
        % Densitat del gas: la ploma puja o baixa amb la distancia a sotavent
        dhS = desplacament_flotacio(max(s, 0) / V, flot);
        for iz = 1:dom.Nz
            for iz0 = 1:dom.Nz
                z0e = max(dom.z(iz0) + dhS, 0);
                r1 = max(sqrt(rh2 + (a * (dom.z(iz) - z0e)).^2), rmin);
                r2 = max(sqrt(rh2 + (a * (dom.z(iz) + z0e)).^2), rmin);
                c = cc.Q / (4 * pi * sqrt(Dh * Dv)) * ...
                    (exp(s * V / (2 * Dh) - r1 / lam) ./ r1 + exp(s * V / (2 * Dh) - r2 / lam) ./ r2);
                K(:, :, iz, iz0) = c;
            end
        end
        info.descripcio = sprintf('adveccio-difusio 3D (Dh = %.1f, Dv = %.2f m^2/s, lambda = %.1f m)', Dh, Dv, lam);
        info.Dh = Dh;
        info.Dv = Dv;
        info.lambda = lam;

    otherwise
        error('Tipus de ploma desconegut: %s', cfg.tipusPloma);
end

% Limitacions del sensor que el dron coneix (dilucio mitjana, saturacio)
K = model_sensor_realista(K, cfg);
% Comptes esperats durant la finestra de mesura
K = dom.tMesura * (cfg.infotaxis.factorModel * cfg.sensor.sensibilitat * K + cfg.sensor.fons);
end
