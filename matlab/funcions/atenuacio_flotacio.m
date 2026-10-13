function a = atenuacio_flotacio(t, flot, H, s2v)
%ATENUACIO_FLOTACIO  Fraccio de contaminant que queda a l'alcada de vol
%   (versio 2D, ploma complexa) quan la ploma puja o baixa per la densitat
%   del gas. t: temps des de la fuga (s); H: alcada de la fuga i del vol
%   en 2D (m) -un gas pesant no pot baixar per sota del terra-; s2v:
%   amplada vertical^2 de la ploma (m^2).
%
%       a = exp(-dh^2 / (2 s2v)),   dh = max(desplacament, -H)
if isempty(flot)
    a = ones(size(t));
    return;
end
dh = max(desplacament_flotacio(t, flot), -H);
a = exp(-dh.^2 ./ (2 * s2v));
end
