function dh = desplacament_flotacio(t, flot)
%DESPLACAMENT_FLOTACIO  Desplacament vertical (m) de la ploma per la
%   densitat del gas, un temps t (s) despres de sortir de la fuga.
%
%   Model simplificat: el gas surt amb una velocitat vertical
%   w0 = signe(g') sqrt(|g'| L), on g' = g (rho_aire - rho_gas) / rho_aire
%   es la "gravetat reduida" i L la mida de la fuga. A mesura que es barreja
%   amb l'aire la diferencia de densitat desapareix i la velocitat decau
%   amb un temps tau:  w(t) = w0 exp(-t/tau), i per tant
%
%       dh(t) = w0 tau (1 - exp(-t/tau))
%
%   Positiu = puja (gas lleuger), negatiu = baixa (gas pesant).
%   Si flot es buit (gas neutre o mode realista desactivat), dh = 0.
if isempty(flot)
    dh = zeros(size(t));
    return;
end
dh = flot.w0 * flot.tau * (1 - exp(-max(t, 0) / flot.tau));
end
