function MU = mu_local(K, q, dom)
%MU_LOCAL  Comptes esperats a la cel.la q = [iy ix] per a cada possible
%   posicio de la font r0 de la graella (matriu Ny x Nx).
%   MU(iy0, ix0) = K(iy - iy0 + Ny, ix - ix0 + Nx).
MU = K(q(1) + dom.Ny - 1:-1:q(1), q(2) + dom.Nx - 1:-1:q(2));
end
