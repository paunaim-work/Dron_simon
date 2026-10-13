function C = ploma_camp(pl, dom)
%PLOMA_CAMP  Concentracio sobre la graella fina de visualitzacio.

switch pl.tipus
    case 'gaussiana'
        C = pl.campVisual;
    case 'complexa'
        if isempty(pl.px)
            C = zeros(size(dom.Xv));
            return;
        end
        % Cada puff es separable en x i y: C = Ey * diag(w) * Ex
        w  = pl.pm ./ (2 * pi * pl.ps2);
        Ey = exp(-(dom.yv(:) - pl.py').^2 ./ (2 * pl.ps2'));   % Ny x Np
        Ex = exp(-(pl.px - dom.xv).^2 ./ (2 * pl.ps2));        % Np x Nx
        C  = (Ey .* w') * Ex;
end
end
