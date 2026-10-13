function C = ploma_camp_3d(pl, dom)
%PLOMA_CAMP_3D  Concentracio sobre la graella 3D de visualitzacio
%   (mida numel(y3) x numel(x3) x numel(z3)).

switch pl.tipus
    case 'gaussiana'
        C = pl.camp3d;
    case 'complexa'
        C = zeros(size(dom.X3));
        if isempty(pl.px)
            return;
        end
        % Cada puff es separable: C(y,x,z) = sum_p w_p Ey(y,p) Ex(p,x) Ez(z,p)
        w  = pl.pm ./ ((2 * pi)^1.5 * pl.ps2h .* sqrt(pl.ps2v));
        Ey = exp(-(dom.y3(:) - pl.py').^2 ./ (2 * pl.ps2h'));
        Ex = exp(-(pl.px - dom.x3).^2 ./ (2 * pl.ps2h));
        Ez = exp(-(dom.z3(:) - pl.pz').^2 ./ (2 * pl.ps2v')) ...
           + exp(-(dom.z3(:) + pl.pz').^2 ./ (2 * pl.ps2v'));
        for k = 1:numel(dom.z3)
            C(:, :, k) = (Ey .* (w' .* Ez(k, :))) * Ex;
        end
end
end
