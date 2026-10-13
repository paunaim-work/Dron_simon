function [K, info] = model_taxa_deteccio(cfg, dom)
%MODEL_TAXA_DETECCIO  Model intern del dron per a l'infotaxis: nombre
%   esperat de comptes en un pas, mu(r | r0) = dt * R(r | r0), quan el dron
%   es a r i la font a r0.
%
%   Com que el vent mitja es uniforme, mu nomes depen del desplacament
%   r - r0. Es calcula un sol cop sobre tots els desplacaments possibles
%   de la graella: K te mida (2Ny-1) x (2Nx-1) i el desplacament zero es
%   a l'element (Ny, Nx). Vegeu mu_local.m.
%
%   - Ploma gaussiana: el dron coneix la forma gaussiana (model perfecte).
%   - Ploma complexa: el dron usa la solucio mitjana de l'equacio
%     d'adveccio-difusio de l'article (material suplementari, eq. 4-5):
%
%        c(r|r0) = Q/(2 pi D) exp(s V / (2D)) K0(|r - r0| / lambda)
%        lambda  = sqrt( D tau / (1 + V^2 tau / (4 D)) )
%
%     on s es la distancia a sotavent i K0 la funcio de Bessel modificada.
%     El model NO coneix el meandreig, les rafegues ni els remolins: nomes
%     una difusivitat efectiva D. Es la situacio realista de l'article
%     (robustesa davant d'errors de model).

[DX, DY] = meshgrid((-(dom.Nx - 1):(dom.Nx - 1)) * dom.dx, ...
                    (-(dom.Ny - 1):(dom.Ny - 1)) * dom.dx);
info = struct();

switch lower(cfg.tipusPloma)
    case 'gaussiana'
        g = cfg.gauss;
        p.Q = g.Q; p.u = g.u; p.H = g.H; p.z = g.z;
        p.estabilitat = g.estabilitat;
        p.theta = dom.theta;
        p.sigma0 = g.sigma0;
        c = ploma_gaussiana_conc(DX, DY, p);
        info.descripcio = 'gaussiana (model exacte)';

    case 'complexa'
        cc = cfg.complexa;
        D = cfg.infotaxis.Dmodel;
        if isempty(D)
            % Difusivitat efectiva: turbulencia + creixement dels puffs +
            % meandreig (U sigma_theta)^2 T_theta + remolins
            D = cc.Dturb + cc.Dpuff + (cc.U * cc.sigmaTheta * pi / 180)^2 * cc.Ttheta;
            if cc.nRemolins > 0
                D = D + 0.5 * cc.Uremolins^2 * cc.Tremolins;
            end
        end
        V   = cc.U;
        tau = cc.tau;
        lam = sqrt(D * tau / (1 + V^2 * tau / (4 * D)));
        s   = DX * cos(dom.theta) + DY * sin(dom.theta);
        d   = max(hypot(DX, DY), dom.dx / 2);
        % besselk(0, x, 1) = K0(x) exp(x): evita desbordaments numerics
        c = cc.Q / (2 * pi * D) * exp(s * V / (2 * D) - d / lam) .* besselk(0, d / lam, 1);
        info.descripcio = sprintf('adveccio-difusio (D = %.1f m^2/s, lambda = %.1f m)', D, lam);
        info.D = D;
        info.lambda = lam;

    otherwise
        error('Tipus de ploma desconegut: %s', cfg.tipusPloma);
end

K = dom.dt * (cfg.infotaxis.factorModel * cfg.sensor.sensibilitat * c + cfg.sensor.fons);
end
