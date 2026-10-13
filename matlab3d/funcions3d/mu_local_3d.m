function MU = mu_local_3d(K, q, dom)
%MU_LOCAL_3D  Comptes esperats a la cel.la q = [iy ix iz] per a cada
%   possible posicio de la font (matriu Ny x Nx x Nz).
MU = reshape(K(q(1) + dom.Ny - 1:-1:q(1), q(2) + dom.Nx - 1:-1:q(2), q(3), :), ...
             dom.Ny, dom.Nx, dom.Nz);
end
