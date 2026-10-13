function C = ploma_tall_3d(pl, dom, z)
%PLOMA_TALL_3D  Tall horitzontal de la ploma a l'alcada z (graella fina
%   dom.xv x dom.yv), el que es veuria volant a aquesta alcada.

switch pl.tipus
    case 'gaussiana'
        C = ploma_gaussiana_conc_3d(dom.Xv - pl.font(1), dom.Yv - pl.font(2), z, pl.font(3), pl.p);
    case 'complexa'
        if isempty(pl.px)
            C = zeros(size(dom.Xv));
            return;
        end
        w  = pl.pm ./ ((2 * pi)^1.5 * pl.ps2h .* sqrt(pl.ps2v));
        ez = exp(-(z - pl.pz).^2 ./ (2 * pl.ps2v)) + exp(-(z + pl.pz).^2 ./ (2 * pl.ps2v));
        Ey = exp(-(dom.yv(:) - pl.py').^2 ./ (2 * pl.ps2h'));
        Ex = exp(-(pl.px - dom.xv).^2 ./ (2 * pl.ps2h));
        C  = (Ey .* (w .* ez)') * Ex;
end
end
