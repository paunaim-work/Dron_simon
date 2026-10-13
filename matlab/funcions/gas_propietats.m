function g = gas_propietats(nom)
%GAS_PROPIETATS  Massa molar i densitat d'alguns gasos (a 20 C i 1 atm).
%
%   La densitat d'un gas ideal es rho = P M / (R T). Comparada amb la de
%   l'aire (1,20 kg/m^3) diu si la ploma tendeix a pujar o a baixar.
%
%     'heli'  : He,   4,0 g/mol -> 0,17 kg/m^3 (molt mes lleuger: puja)
%     'meta'  : CH4, 16,0 g/mol -> 0,67 kg/m^3 (gas natural: puja)
%     'benze' : C6H6, 78,1 g/mol -> 3,25 kg/m^3 (vapor, objectiu del projecte: baixa)
%     'propa' : C3H8, 44,1 g/mol -> 1,83 kg/m^3 (GLP: baixa)
%     'co2'   : CO2,  44,0 g/mol -> 1,83 kg/m^3 (baixa)

switch lower(nom)
    case 'heli'
        g.nom = 'heli (He)';             g.M = 4.003;
    case 'meta'
        g.nom = 'meta (CH4)';            g.M = 16.04;
    case 'benze'
        g.nom = 'benze (C6H6, vapor)';   g.M = 78.11;
    case 'propa'
        g.nom = 'propa (C3H8)';          g.M = 44.10;
    case 'co2'
        g.nom = 'dioxid de carboni (CO2)'; g.M = 44.01;
    otherwise
        error('Gas desconegut: %s (usa ''heli'', ''meta'', ''benze'', ''propa'' o ''co2'')', nom);
end
P = 101325;          % Pa
T = 293.15;          % K
Rg = 8.314;          % J/(mol K)
g.densitat     = P * g.M / 1000 / (Rg * T);       % kg/m^3
g.densitatAire = P * 28.96 / 1000 / (Rg * T);
end
