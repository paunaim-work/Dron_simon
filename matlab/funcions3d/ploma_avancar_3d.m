function pl = ploma_avancar_3d(pl, dt)
%PLOMA_AVANCAR_3D  Fa evolucionar la ploma 3D un pas de temps dt.

pl.t = pl.t + dt;
if ~strcmp(pl.tipus, 'complexa')
    return;
end

% --- 1. Vent: processos d'Ornstein-Uhlenbeck ---
e = exp(-dt / pl.Ttheta);
pl.dTheta = pl.dTheta * e + pl.sigmaTheta * sqrt(1 - e^2) * randn;
e = exp(-dt / pl.Trafega);
pl.rafega = max(pl.rafega * e + pl.sigmaRafega * sqrt(1 - e^2) * randn, -0.9);
if ~isempty(pl.a)
    e = exp(-dt / pl.Trem);
    pl.a = pl.a * e + pl.sigmaA .* sqrt(1 - e^2) .* randn(size(pl.a));
end
pl.desp = pl.desp + pl.U * [cos(pl.theta0) sin(pl.theta0)] * dt;

% --- 2. Transport dels puffs ---
if ~isempty(pl.px)
    [u, v] = ploma_vent_3d(pl, pl.px, pl.py, pl.pz);
    sh = sqrt(2 * pl.DturbH * dt);
    sv = sqrt(2 * pl.DturbV * dt);
    pl.px = pl.px + u * dt + sh * randn(size(pl.px));
    pl.py = pl.py + v * dt + sh * randn(size(pl.py));
    wv = 0;
    if ~isempty(pl.flot)                              % el gas puja o baixa (mode realista)
        wv = pl.flot.w0 * exp(-pl.page / pl.flot.tau);
    end
    pl.pz = abs(pl.pz + wv * dt + sv * randn(size(pl.pz)));   % rebot al terra (z = 0)
    pl.ps2h = pl.ps2h + 2 * pl.DpuffH * dt;
    pl.ps2v = pl.ps2v + 2 * pl.DpuffV * dt;
    pl.pm   = pl.pm * exp(-dt / pl.tau);
    pl.page = pl.page + dt;
end

% --- 3. Emissio de nous puffs ---
n = max(1, round(pl.puffsPerSeg * dt));
f = rand(n, 1);
[u0, v0] = ploma_vent_3d(pl, pl.font(1), pl.font(2), pl.font(3));
pl.px   = [pl.px;   pl.font(1) + u0 * f * dt];
pl.py   = [pl.py;   pl.font(2) + v0 * f * dt];
pl.pz   = [pl.pz;   pl.font(3) * ones(n, 1)];
pl.ps2h = [pl.ps2h; pl.sigma0^2 + 2 * pl.DpuffH * f * dt];
pl.ps2v = [pl.ps2v; pl.sigma0^2 + 2 * pl.DpuffV * f * dt];
pl.pm   = [pl.pm;   pl.Q * dt / n * ones(n, 1)];
pl.page = [pl.page; f * dt];

% --- 4. Eliminem els puffs que surten o s'han esvait ---
fora = pl.px < pl.lim(1) | pl.px > pl.lim(2) | pl.py < pl.lim(3) | pl.py > pl.lim(4) ...
     | pl.pz > pl.lim(5) | pl.pm < 1e-4 * pl.massaPuff;
pl.px(fora) = [];
pl.py(fora) = [];
pl.pz(fora) = [];
pl.ps2h(fora) = [];
pl.ps2v(fora) = [];
pl.pm(fora) = [];
pl.page(fora) = [];
end
