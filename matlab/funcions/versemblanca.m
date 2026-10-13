function L = versemblanca(k, MU, saturacio)
%VERSEMBLANCA  Probabilitat d'observar k comptes si la mitjana es MU,
%   P(k | r0), per a cada posicio candidata de la font.
%
%   Llei de Poisson: P(k) = MU^k exp(-MU) / k!
%   Si el sensor satura (k >= saturacio), nomes sabem que n'hi ha "molts":
%   P(k >= s) = 1 - sum_{j<s} P(j). Amb saturacio = 1 el sensor es
%   binari (deteccio si/no), com en molts treballs d'infotaxis.

if k < saturacio
    L = exp(-MU + k * log(MU) - gammaln(k + 1));
else
    acc  = zeros(size(MU));
    term = exp(-MU);
    for j = 0:saturacio - 1
        acc  = acc + term;
        term = term .* MU / (j + 1);
    end
    L = max(1 - acc, 0);
end
end
