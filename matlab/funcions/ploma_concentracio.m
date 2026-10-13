function c = ploma_concentracio(pl, x, y)
%PLOMA_CONCENTRACIO  Concentracio instantania de la ploma als punts (x, y).

switch pl.tipus
    case 'gaussiana'
        c = ploma_gaussiana_conc(x - pl.font(1), y - pl.font(2), pl.p);
    case 'complexa'
        if isempty(pl.px)
            c = zeros(size(x));
            return;
        end
        % Suma de puffs gaussians 2D: m / (2 pi s^2) exp(-d^2 / (2 s^2))
        w  = pl.pm ./ (2 * pi * pl.ps2) .* atenuacio_flotacio(pl.page, pl.flot, pl.Hflot, ...
                                                          pl.sigma0^2 + 2 * pl.Dv * pl.page);
        d2 = (x(:) - pl.px').^2 + (y(:) - pl.py').^2;
        c  = reshape(exp(-d2 ./ (2 * pl.ps2')) * w, size(x));
end
end
