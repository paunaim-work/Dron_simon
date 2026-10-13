function R = configuracio_realista()
%CONFIGURACIO_REALISTA  Limitacions tecniques del dron (mode realista).
%
%   Valors per defecte presos de la fitxa tecnica del Dron Simon:
%   quadricopter de ~4 kg, bateria Li-ion 6S4P de 389 Wh (~390 W en vol
%   estacionari), sensor PID ION Science MiniPID 2 de 10,0 eV (0-100 ppm,
%   T90 < 8 s), tub de presa de mostra de 80 cm que surt fora del corrent
%   de les helices, GPS Holybro M10 i anemometre ultrasonic 3D.
%
%   Cada limitacio te un interruptor (activat = true/false) i un o mes
%   valors regulables. R.activat es l'interruptor general: si es false,
%   la simulacio es la ideal (com sempre).
%
%   En mode realista la concentracio de la simulacio s'interpreta en ppm.

R.activat = false;              % interruptor general del mode realista
R.mostrar = true;               % escriu el resum de les limitacions a la consola

% --- Velocitat de vol durant la cerca ---
% La fitxa recomana volar a poc a poc: amb T90 = 8 s, a 1 m/s la mesura
% queda "difuminada" en uns 8 m.
R.velocitat.activat = true;
R.velocitat.valor   = 2;        % m/s

% --- Temps aturat a cada punt per mesurar ---
% El dron s'atura en vol estacionari i acumula la mesura durant aquest
% temps (mes senyal i menys soroll, pero cada pas dura mes).
R.tempsMesura.activat = true;
R.tempsMesura.valor   = 8;      % s

% --- Temps de resposta del sensor PID ---
% La lectura segueix la concentracio real amb retard (filtre de primer
% ordre). T90 = temps que triga a arribar al 90 % d'un canvi brusc.
R.respostaSensor.activat = true;
R.respostaSensor.T90     = 8;   % s

% --- Corrent d'aire de les helices (downwash) ---
% Les helices empenyen aire cap avall i el barregen amb la mostra: la
% dilueixen (hi entra aire net) i hi afegeixen fluctuacions. Amb el tub de
% 80 cm que surt fora del corrent l'efecte es petit; sense tub seria molt
% mes gran (dilucio ~0.5).
R.helices.activat = true;
R.helices.dilucio = 0.10;       % fraccio d'aire net a la mostra (0 - 1)
R.helices.soroll  = 0.20;       % variabilitat relativa de la mesura (coeficient de variacio)

% --- Gas que es busca ---
% La densitat del gas fa que la ploma pugi (si es mes lleuger que l'aire)
% o baixi (si es mes pesant) fins que es dilueix. Opcions:
% 'heli', 'meta', 'benze', 'propa', 'co2' (vegeu gas_propietats.m).
R.gas.activat      = true;
R.gas.nom          = 'benze';
R.gas.midaFont     = 0.5;       % mida caracteristica de la fuga (m)
R.gas.tempsDilucio = 2;         % s: temps en que el gas es barreja amb l'aire i deixa de pujar o baixar

% --- Bateria i autonomia ---
% Autonomia = fraccioUtil * energia / potencia (fitxa: ~48 min en estacionari).
R.bateria.activat       = true;
R.bateria.energia       = 389;  % Wh
R.bateria.fraccioUtil   = 0.8;  % fraccio que es gasta (no es descarrega del tot)
R.bateria.potencia      = 390;  % W en vol estacionari (motors + electronica + sensors)
R.bateria.extraMoviment = 0.10; % consum addicional quan es desplaca (fraccio)
R.bateria.retornBase    = true; % guarda energia per tornar al punt d'enlairament

% --- Rang del sensor ---
R.saturacio.activat = true;
R.saturacio.valor   = 100;      % ppm: per sobre, la lectura satura

% --- Error de posicio del GPS ---
% La mesura es pren en un punt lleugerament diferent del que creu el dron.
R.errorGPS.activat = true;
R.errorGPS.valor   = 1.5;       % m (desviacio tipica)

% --- Error de l'anemometre en la direccio del vent ---
% El model intern del dron fa servir una direccio de vent desviada. Un
% anemometre ultrasonic sol tenir +-1-2 graus; al dron s'hi afegeix el
% corrent de les helices i el seu propi moviment.
R.errorVent.activat = true;
R.errorVent.valor   = 3;        % graus (desviacio tipica)

R.dtSub = 1;                    % s: pas de temps intern per simular la mesura
end
