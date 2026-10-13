function k = poisson_aleatori(mu)
%POISSON_ALEATORI  Mostres de Poisson amb mitjana mu (sense toolboxes).
%
%   Per a mu < 30 s'usa el metode de Knuth (producte d'uniformes); per a
%   mu gran, l'aproximacio normal arrodonida.

k = zeros(size(mu));

idx = find(mu < 30);
if ~isempty(idx)
    L  = exp(-mu(idx));
    p  = rand(size(idx));
    kk = zeros(size(idx));
    actiu = p > L;
    while any(actiu)
        kk(actiu) = kk(actiu) + 1;
        p(actiu)  = p(actiu) .* rand(nnz(actiu), 1);
        actiu = p > L;
    end
    k(idx) = kk;
end

gran = find(mu >= 30);
if ~isempty(gran)
    k(gran) = max(0, round(mu(gran) + sqrt(mu(gran)) .* randn(size(gran))));
end
end
