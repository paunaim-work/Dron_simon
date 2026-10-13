# Dron autònom per trobar l'origen d'una fuga

Simulació en MATLAB d'un dron que busca l'origen d'una fuga de contaminant (química o radioactiva) a partir només de les lectures del seu sensor. Es pot triar el **model de ploma** (gaussiana o complexa) i l'**algorisme de cerca** (infotaxis, de l'article de Vergassola *et al.* 2007, o gradient), i es veu en directe com el dron fa la cerca.

Projecte: *Disseny i programació d'un dron aeri autònom per a la mesura de contaminació i determinació del seu origen* (Ignasi Reichardt, 2026-2027).

| Ploma gaussiana + infotaxis | Ploma complexa + infotaxis |
|---|---|
| ![](figures/exemple_gaussiana_infotaxis.png) | ![](figures/exemple_complexa_infotaxis.png) |
| **Ploma gaussiana + gradient** | **Ploma complexa + gradient** |
| ![](figures/exemple_gaussiana_gradient.png) | ![](figures/exemple_complexa_gradient.png) |

## Com s'executa

1. Obre MATLAB (R2016b o posterior; no cal cap toolbox) a la carpeta `matlab/`.
2. Obre `main_simulacio.m` i tria les opcions al principi:
   ```matlab
   cfg.tipusPloma = 'gaussiana';   % 'gaussiana'  o  'complexa'
   cfg.algorisme  = 'infotaxis';   % 'infotaxis'  o  'gradient'
   cfg.llavor     = 1;             % canvia-la per veure una altra realització
   ```
3. Executa'l (F5). S'obre una finestra amb l'animació i, en acabar, surt un resum a la consola.

Si poses `preguntarOpcions = true`, el programa et demana la ploma i l'algorisme per la finestra de comandes. Per desar un vídeo `.mp4`, posa `cfg.visual.gravarVideo = true`.

Per comparar els algorismes estadísticament (mètode de **Monte Carlo**: repetir la simulació moltes vegades amb nombres aleatoris diferents i fer estadística dels resultats), executa `comparativa_montecarlo.m`. Fa 20 cerques per combinació sense animació i en mostra la taxa d'èxit i el temps de cerca.

També funciona amb GNU Octave (provat amb la versió 8.4). Com que Octave dibuixa més lentament, val la pena posar `cfg.visual.cadaNPassos = 5`.

## Què es veu a l'animació

- **Panell gran (esquerra):** la ploma real en escala logarítmica de colors, la font (estrella verda), el punt d'inici (triangle), el dron (cercle blau), la seva trajectòria i els punts on ha detectat contaminant (punts negres). La fletxa indica el vent (a la ploma complexa canvia de direcció). Amb infotaxis, la creu lila és on el dron creu que és més probable que hi hagi la font; amb gradient, la fletxa lila és el gradient que ha estimat.
- **Panell superior dret:** el que "sap" el dron. Amb infotaxis és la probabilitat P(r₀) que la font sigui a cada punt; es veu com es concentra a mesura que el dron rep informació. Amb gradient és el mapa de les concentracions que ha mesurat.
- **Panells inferiors drets:** la distància a la font i l'entropia (infotaxis) o el nombre de comptes del sensor (gradient).

## Estructura

```
matlab/
├── main_simulacio.m           ← script principal (tria ploma i algorisme)
├── comparativa_montecarlo.m   ← estadística de moltes cerques
└── funcions/
    ├── configuracio_defecte.m ← TOTS els paràmetres (amb comentaris)
    ├── simular_cerca.m        ← bucle de simulació
    ├── ploma_*.m, sigma_briggs.m   ← models de ploma
    ├── sensor_llegir.m, poisson_aleatori.m   ← model del sensor
    ├── infotaxis_*.m, model_taxa_deteccio.m, versemblanca.m, entropia.m
    ├── gradient_*.m
    └── visualitzacio_*.m, mapa_colors.m
```

## Models de ploma

### 1. Ploma gaussiana (`'gaussiana'`)

És el model clàssic de dispersió atmosfèrica: una font puntual que emet de manera contínua, amb un vent constant, dona un núvol estacionari (que no canvia amb el temps) amb forma de con. La concentració a l'alçada de vol del dron *z* és:

$$c(s,n)=\frac{Q}{2\pi u\,\sigma_y\sigma_z}\,e^{-\frac{n^2}{2\sigma_y^2}}\left[e^{-\frac{(z-H)^2}{2\sigma_z^2}}+e^{-\frac{(z+H)^2}{2\sigma_z^2}}\right]$$

on *s* és la distància a sotavent (en la direcció del vent), *n* la distància transversal, *Q* la taxa d'emissió, *u* la velocitat del vent i *H* l'alçada de la font. El segon terme del claudàtor representa la "reflexió" del contaminant al terra. Les amplades σ_y(s) i σ_z(s) creixen amb la distància segons la **classe d'estabilitat de Pasquill-Gifford** (de A, atmosfera molt turbulenta que dispersa molt, fins a F, atmosfera molt estable), amb les fórmules de Briggs (`sigma_briggs.m`). Per evitar que la concentració sigui infinita just a la font, aquesta té una mida inicial σ₀ = 1 m.

Aquest model representa el **camp mitjà**. L'única aleatorietat és el soroll de recompte del sensor, i per això és el cas "fàcil".

### 2. Ploma complexa (`'complexa'`)

És un **model lagrangià de *puffs*** (Farrell *et al.*, 2002). "Lagrangià" vol dir que se segueix el moviment de porcions de contaminant en comptes de calcular un camp sobre una graella. La font emet contínuament petits núvols gaussians (*puffs*) que:

- són **transportats pel vent**, que no és constant: la direcció *meandreja* (oscil·la) i la intensitat té *ràfegues*. Totes dues es modelen amb processos d'**Ornstein-Uhlenbeck**, un soroll aleatori "amb memòria" que varia de manera suau, amb un temps de correlació característic;
- són desviats per **remolins de gran escala**: un camp de velocitat sense divergència (que no crea ni destrueix aire), obtingut d'una funció de corrent ψ, que viatja amb el vent mitjà;
- es **dispersen per turbulència**: el centre fa un camí aleatori (difusió turbulenta) i el *puff* s'eixampla (difusió relativa);
- **decauen** amb un temps de vida τ, per desintegració radioactiva, deposició o reacció química.

El resultat és una ploma sinuosa i filamentosa. Lluny de la font, el contaminant arriba a pegats separats per zones buides: és el **límit diluït** de l'article, on les deteccions són esporàdiques i no hi ha un gradient útil.

## Model del sensor

El dron porta un sensor de **recompte** (com un comptador Geiger o un detector de partícules): durant cada pas de temps Δt compta un nombre aleatori de deteccions *k* que segueix una **llei de Poisson** (la llei dels esdeveniments aleatoris independents):

$$k\sim\text{Poisson}(\mu),\qquad \mu=\Delta t\,(\kappa\,c+b)$$

on κ és la sensibilitat del sensor i *b* el **soroll de fons** (deteccions falses o radioactivitat natural). Quan la concentració és baixa, el sensor sovint marca 0 encara que hi hagi contaminant: aquest soroll de recompte és el que fa difícil la cerca.

## Algorismes de cerca

El dron es mou per una graella de cel·les de mida dx (5 m per defecte). A cada pas pot anar a una de les 4 cel·les veïnes o quedar-se quiet, i triga Δt = dx / velocitat. La fuga es considera trobada quan el dron hi és a menys d'un **radi de detecció** (1,5·dx), a una distància on la podria veure amb una càmera.

### Infotaxis (`'infotaxis'`), de l'article

Idea: en comptes de seguir la concentració, el dron segueix la **informació**.

1. **Probabilitat posterior P(r₀).** El dron manté un mapa amb la probabilitat que la font sigui a cada punt r₀. Al principi és uniforme (no en sap res). La posterior és la probabilitat *actualitzada* amb les mesures. Cada cop que mesura *k* comptes, l'actualitza amb la **regla de Bayes** (eq. 2 de l'article):

   $$P_{t+1}(r_0)\propto P_t(r_0)\cdot\text{Poisson}\big(k;\,\mu(r\,|\,r_0)\big)$$

   És a dir, augmenta la probabilitat de les posicions de la font que expliquen bé el que ha mesurat. A més, la zona on és el dron passa a tenir probabilitat 0, perquè si la font hi fos ja l'hauria trobada.

2. **Entropia de Shannon** $S=-\sum P\ln P$. Mesura la incertesa: és màxima quan P és uniforme i val 0 quan se sap exactament on és la font. Segons l'article, el temps esperat de cerca creix com $e^{S}$, així que reduir l'entropia és acostar-se a trobar la font.

3. **Decisió (eq. 1 de l'article).** Per a cada moviment possible (4 veïns o quedar-se quiet), el dron calcula quant espera que baixi l'entropia:

   $$\overline{\Delta S}(r\to r_j)=P(r_j)\,[-S]+\big[1-P(r_j)\big]\sum_k\rho_k(r_j)\,\Delta S_k$$

   - El primer terme és l'**explotació**: la probabilitat de trobar la font just allà (l'entropia passaria a ser 0).
   - El segon terme és l'**exploració**: si no la troba, rebrà *k* comptes amb probabilitat ρ_k, i cada resultat canviaria l'entropia en ΔS_k.

   Tria el moviment amb el $\overline{\Delta S}$ més negatiu, és a dir, el que maximitza el guany d'informació esperat. Aquest equilibri entre explotar i explorar és el que produeix les trajectòries en ziga-zaga (*casting*) semblants a les de les arnes.

4. **Model intern** $\mu(r|r_0)$. El dron necessita saber quants comptes esperaria rebre a *r* si la font fos a *r₀* (`model_taxa_deteccio.m`):
   - Ploma gaussiana: fa servir la mateixa fórmula gaussiana (model perfecte).
   - Ploma complexa: fa servir la solució mitjana de l'**equació d'advecció-difusió** de l'article (material suplementari, eq. 5). "Advecció" és el transport pel vent i "difusió", l'escampament per turbulència:
     $$c(r|r_0)=\frac{Q}{2\pi D}\,e^{\frac{sV}{2D}}\,K_0\!\left(\frac{|r-r_0|}{\lambda}\right),\qquad \lambda=\sqrt{\frac{D\tau}{1+V^2\tau/4D}}$$
     on $K_0$ és la funció de Bessel modificada. Aquest model **no** coneix el meandreig, les ràfegues ni els remolins; només en té una difusivitat efectiva *D*. Serveix per comprovar la robustesa de l'algorisme davant d'errors de model, tal com fa l'article.

Detalls d'implementació:
- Per calcular el sumatori de l'eq. 1 sense un cost excessiu quan hi ha molts comptes, els recomptes s'agrupen en intervals (0, 1, 2, 3-4, 5-7, …) (`cfg.infotaxis.binsComptes`). L'actualització bayesiana sí que fa servir el valor exacte de *k*.
- Amb `cfg.infotaxis.saturacio = 1` el sensor és binari (detecta / no detecta), com en molts treballs posteriors a l'article.
- Si la posterior ja està tan concentrada que cap moviment aporta informació (entropia ≈ 0 lluny de la font), el dron vola directament cap al punt més probable. Sense aquest pas, quedaria aturat.

### Gradient (`'gradient'`), quimiotaxi clàssica

El dron porta 5 sensors: un al centre i 4 en creu (a +x, −x, +y, −y), com si n'hi hagués un a cada braç d'un quadricòpter. Estima el gradient (la direcció on la concentració augmenta més) per diferències finites:

$$\frac{\partial c}{\partial x}\approx\frac{c_E-c_O}{2d},\qquad\frac{\partial c}{\partial y}\approx\frac{c_N-c_S}{2d}$$

i el suavitza amb una mitjana mòbil exponencial per reduir el soroll. Es mou cap al veí en la direcció del gradient (mode *seguiment*). Si fa uns quants passos que no detecta res, fa un **passeig aleatori correlacionat** (manté la direcció amb una certa probabilitat i, si no, gira) fins que torna a trobar senyal (mode *exploració*).

Com explica el material suplementari de l'article (secció 1), mesurar un gradient amb soroll de recompte requereix molt de senyal, i el temps necessari creix exponencialment amb la distància a la font. Per això el gradient funciona bé a prop de la font d'una ploma suau, però falla en el límit diluït o amb plomes turbulentes.

## Resultats típics

Amb els paràmetres per defecte (domini de 400 × 300 m, dron a ~290 m de la fuga, 5 m/s), sobre 5 realitzacions:

| Ploma | Algorisme | Èxit | Temps de cerca |
|---|---|---|---|
| gaussiana | infotaxis | 5/5 | 120 – 260 s |
| gaussiana | gradient | 5/5 | 160 – 870 s |
| complexa | infotaxis | 5/5 | 220 – 350 s |
| complexa | gradient | 4/5 | 650 – 1900 s (o no la troba) |

L'infotaxis és més ràpid i més fiable, sobretot amb la ploma complexa, on el gradient sovint es perd. Per a estadístiques més sòlides, executa `comparativa_montecarlo.m`.

## Paràmetres principals (`configuracio_defecte.m`)

| Paràmetre | Defecte | Significat |
|---|---|---|
| `domini.Lx, Ly, dx` | 400, 300, 5 m | mida del domini i de la cel·la |
| `font.pos` | [60 150] m | posició real de la fuga (el dron no la coneix) |
| `dron.posInicial` | [340 215] m | punt de sortida |
| `dron.velocitat` | 5 m/s | defineix Δt = dx/v |
| `dron.maxPassos` | 2500 | límit d'autonomia |
| `vent.direccio` | 0° | cap on bufa el vent |
| `sensor.sensibilitat` | 1 | comptes/s per unitat de concentració (baixa-la per fer el senyal més diluït) |
| `sensor.fons` | 0.02 /s | deteccions falses |
| `gauss.estabilitat` | 'D' | classe A..F |
| `complexa.sigmaTheta` | 15° | amplitud del meandreig del vent |
| `complexa.nRemolins` | 6 | remolins de gran escala (0 = cap) |
| `infotaxis.saturacio` | Inf | 1 = sensor binari |
| `infotaxis.factorModel` | 1 | per simular un error de calibratge del model intern |
| `gradient.alfa` | 0.6 | filtre del gradient (0 = cap) |

## Idees per continuar el projecte

- Estudiar com depèn el temps de cerca de la distància inicial, de la sensibilitat del sensor o de la velocitat del dron (pla de treball de novembre i desembre).
- Comprovar la robustesa de l'infotaxis davant d'errors en el model intern (`factorModel`, `Dmodel`, direcció del vent mal estimada).
- Aprendre els paràmetres del medi a partir de les deteccions (secció 5 del material suplementari).
- Tenir en compte les correlacions temporals entre deteccions (secció 6, eq. 12-13), que són importants a la ploma complexa.
- Ampliar a 3D (alçada de vol variable) i a moviments en 8 direccions.

## Referències

- M. Vergassola, E. Villermaux, B. I. Shraiman, "'Infotaxis' as a strategy for searching without gradients", *Nature* **445**, 406–409 (2007), i el seu material suplementari.
- J. A. Farrell, J. Murlis, X. Long, W. Li, R. T. Cardé, "Filament-based atmospheric dispersion model to achieve short time-scale structure of odor plumes", *Environmental Fluid Mechanics* **2**, 143–169 (2002).
- G. A. Briggs, "Diffusion estimation for small emissions", ATDL Report 79 (1973).
- H. C. Berg, E. M. Purcell, "Physics of chemoreception", *Biophysical Journal* **20**, 193–219 (1977).
