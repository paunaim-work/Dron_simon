function c = ploma_concentracio_3d(pl, x, y, z)
%PLOMA_CONCENTRACIO_3D  Concentracio instantania als punts (x, y, z).

switch pl.tipus
    case 'gaussiana'
        c = ploma_gaussiana_conc_3d(x - pl.font(1), y - pl.font(2), z, pl.font(3), pl.p);
    case 'complexa'
        if isempty(pl.px)
            c = zeros(size(x));
            return;
        end
        % Puff gaussia 3D + la seva imatge sota terra (reflexio):
        % m / ((2pi)^(3/2) sh^2 sv) exp(-dh^2/2sh^2) [e^(-(z-zp)^2/2sv^2) + e^(-(z+zp)^2/2sv^2)]
        w  = pl.pm ./ ((2 * pi)^1.5 * pl.ps2h .* sqrt(pl.ps2v));
        dh = (x(:) - pl.px').^2 + (y(:) - pl.py').^2;
        ez = exp(-(z(:) - pl.pz').^2 ./ (2 * pl.ps2v')) + exp(-(z(:) + pl.pz').^2 ./ (2 * pl.ps2v'));
        c  = reshape((exp(-dh ./ (2 * pl.ps2h')) .* ez) * w, size(x));
end
end
