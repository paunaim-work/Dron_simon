% Exemple: ploma gaussiana estacionària a z = H.
setup;
[X, Y] = meshgrid(linspace(1, 500, 300), linspace(-100, 100, 200));
C = ploma_gaussiana(X, Y, 10);
pcolor(X, Y, C); shading flat; colorbar;
xlabel('x (m)'); ylabel('y (m)'); title('Ploma gaussiana (z = H)');
