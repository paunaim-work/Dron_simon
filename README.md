# Dron_simon

**Disseny i programació d'un dron aeri autònom per a la mesura de contaminació i determinació del seu origen**

Projecte 5 (curs 2026-2027) · Àmbit: enginyeria ambiental i teoria de la informació · Professor responsable: Ignasi Reichardt

## Descripció

L'objectiu és seleccionar els components per construir un dron aeri capaç de transportar un sistema de mesura de contaminació química i/o radioactiva fins a la zona afectada per una hipotètica fuita a l'aire, i dissenyar els algorismes que permetin que el dron trobi **l'origen de la fuita de manera autònoma**.

El dron mesura un camp escalar variable (la concentració del contaminant) i, a partir de les mesures, decideix iterativament cap on moure's per estimar la posició de la font (*source term estimation*, també anomenat *backtracking* de la ploma).

## Metodologia

1. Recerca bibliogràfica sobre drons per a la mesura de camps escalars variables i les seves característiques tècniques.
2. Mètodes estocàstics i generació de programari per simular **plomes** (*plumes*: la columna de contaminant dispersada pel vent) de complexitat creixent.
3. Mètodes de teoria de la informació (p. ex. *infotaxis*, maximització de l'entropia d'informació) o altres propostes per al *backtracking* de les plomes.

## Resultats acadèmics esperats

- Gestió d'un projecte de disseny i construcció d'un dispositiu de mesura autònom.
- Programari per produir dades sintètiques i analitzar-les.
- Comprensió de les plomes de contaminació (generació i mecanismes difusius).
- Ús de mètodes de Monte-Carlo aplicats a l'exploració aleatòria d'un camp escalar per determinar-ne l'origen.

## Pla de treball

| Mes | Tasques |
|---|---|
| Setembre | Estudi del context i treballs previs. Esquematització de les característiques del dron. |
| Octubre | Generació de plomes sintètiques. Disseny de l'algorisme d'estimació iterativa de l'origen. |
| Novembre | Aplicació de l'algorisme a plomes de complexitat creixent. Idoneïtat de les característiques del dron (velocitat, maniobrabilitat, autonomia). |
| Desembre | Limitacions de l'algorisme segons la complexitat de la ploma. Valoració econòmica. Conclusions. |

Detall a [`docs/pla_de_treball.md`](docs/pla_de_treball.md).

## Estructura del repositori

```
Dron_simon/
├── docs/                # Enunciat, pla de treball, notes i memòria
├── refs/                # Bibliografia (BibTeX i notes de lectura)
├── matlab/
│   ├── plomes/          # Simulació de plomes (MATLAB)
│   └── algorisme/       # Algorismes de cerca de l'origen (MATLAB)
├── python/
│   ├── dron_simon/      # Paquet Python (plomes, algorismes, utilitats)
│   └── tests/           # Tests
├── data/
│   ├── sintetic/        # Dades generades per simulació
│   └── reals/           # Dades mesurades (si n'hi ha)
└── figures/             # Figures generades
```

## Com començar

### Python (≥ 3.10)

```bash
cd python
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m pytest
python -m dron_simon.plomes     # exemple: ploma gaussiana estacionària
```

### MATLAB (R2021a o posterior)

```matlab
cd matlab
setup            % afegeix les carpetes al path
exemple_ploma    % exemple: ploma gaussiana estacionària
```

## Convencions

- Idioma: català (codi, comentaris, documentació).
- Unitats SI (m, s, kg, Bq); coordenades *x* en la direcció del vent, *y* transversal, *z* vertical.
- Les dades sintètiques grans no es versionen (vegeu `.gitignore`); es regeneren amb les llavors aleatòries fixades als scripts.
