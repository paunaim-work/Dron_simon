function pl = ploma_avancar(pl, dt)
%PLOMA_AVANCAR  Fa evolucionar la ploma un pas de temps dt.
%   La ploma gaussiana es estacionaria; la complexa mou, fa creixer i
%   fa decaure els puffs, n'emet de nous i actualitza el vent.

pl.t = pl.t + dt;
if ~strcmp(pl.tipus, 'complexa')
    return;
end

% --- 1. Vent: processos d'Ornstein-Uhlenbeck (soroll amb memoria) ---
e = exp(-dt / pl.Ttheta);
pl.dTheta = pl.dTheta * e + pl.sigmaTheta * sqrt(1 - e^2) * randn;
e = exp(-dt / pl.Trafega);
pl.rafega = max(pl.rafega * e + pl.sigmaRafega * sqrt(1 - e^2) * randn, -0.9);
if ~isempty(pl.a)
    e = exp(-dt / pl.Trem);
    pl.a = pl.a * e + pl.sigmaA .* sqrt(1 - e^2) .* randn(size(pl.a));
end
pl.desp = pl.desp + pl.U * [cos(pl.theta0) sin(pl.theta0)] * dt;

% --- 2. Transport dels puffs existents (adveccio + difusio turbulenta) ---
if ~isempty(pl.px)
    [u, v] = ploma_vent(pl, pl.px, pl.py);
    sd = sqrt(2 * pl.Dturb * dt);
    pl.px  = pl.px + u * dt + sd * randn(size(pl.px));
    pl.py  = pl.py + v * dt + sd * randn(size(pl.py));
    pl.ps2 = pl.ps2 + 2 * pl.Dpuff * dt;      % el puff s'eixampla
    pl.pm  = pl.pm * exp(-dt / pl.tau);       % temps de vida finit
    pl.page = pl.page + dt;                   % edat del puff
end

% --- 3. Emissio de nous puffs (repartits al llarg del pas) ---
n = max(1, round(pl.puffsPerSeg * dt));
f = rand(n, 1);                                % fraccio del pas ja viatjada
[u0, v0] = ploma_vent(pl, pl.font(1), pl.font(2));
pl.px  = [pl.px;  pl.font(1) + u0 * f * dt];
pl.py  = [pl.py;  pl.font(2) + v0 * f * dt];
pl.ps2 = [pl.ps2; pl.sigma0^2 + 2 * pl.Dpuff * f * dt];
pl.pm  = [pl.pm;  pl.Q * dt / n * ones(n, 1)];
pl.page = [pl.page; f * dt];

% --- 4. Eliminem els puffs que surten del domini o s'han esvait ---
fora = pl.px < pl.lim(1) | pl.px > pl.lim(2) | pl.py < pl.lim(3) | pl.py > pl.lim(4) ...
     | pl.pm < 1e-4 * pl.massaPuff;
pl.px(fora) = [];
pl.py(fora) = [];
pl.ps2(fora) = [];
pl.pm(fora) = [];
pl.page(fora) = [];
end
