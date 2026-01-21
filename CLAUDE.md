# fractall

Visualiseur de fractales portable en C utilisant SDL.

**Licence** : GPL-2.0 | **Auteur** : Arnaud VERHILLE (2001-2026) | **Version** : 1.0

## Compilation

```bash
./autogen.sh && ./configure && make && ./src/fractall
```

### Options de configuration

| Option | Description | Défaut |
|--------|-------------|--------|
| `--with-gmp` / `--without-gmp` | Haute précision (zooms profonds) | activé |
| `--disable-openmp` | Désactive le calcul parallèle | activé |
| `--disable-simd` | Désactive SIMD (SSE4.1/AVX) | activé |
| `--with-png` / `--without-png` | Export PNG avec métadonnées | activé |

**Dépendances** : SDL 1.2.0+, libm, GMP (opt), OpenMP (opt), libpng (opt)

## Utilisation

```bash
fractall [-x<largeur>] [-y<hauteur>] [-f] [-g<hauteur_gui>] [-nogui] [-h]
```

### Contrôles

| Touche | Effet |
|--------|-------|
| **F1-F12** | Von Koch, Dragon, Mandelbrot, Julia, Julia Sin, Newton, Phoenix, Buffalo, Burning Ship, Tricorn, Mandelbulb, Buddhabrot |
| **GUI** | 23 boutons (types 1-23) pour accéder à toutes les fractales |
| **Clic gauche** | Zoom / +1 itération (vectorielles) |
| **Clic gauche + glisser** | Zoom par sélection rectangulaire (escape-time) |
| **Clic droit** | Dézoom / -1 itération |
| **C** | Changer palette (9 palettes) |
| **R** | Répétitions gradient (2-40, pas de 2) |
| **S** | Screenshot PNG avec métadonnées |
| **Q/ESC** | Quitter |

## Architecture

```
src/
├── main.c              # Point d'entrée, événements, zoom par sélection
├── EscapeTime.[ch]     # Fractales escape-time (23 types)
├── colorization.[ch]   # Colorisation unifiée (9 palettes)
├── color_types.h       # Type color unifié
├── VonKoch.[ch]        # Fractales vectorielles
├── SDLGUI.[ch]         # Interface graphique
├── complexmath.[ch]    # Arithmétique complexe (double)
├── complexmath_gmp.[ch]# Arithmétique GMP
├── complexmath_simd.c  # SIMD (SSE4.1/AVX)
├── precision_detector.[ch] # Détection précision GMP
├── gmp_pool.[ch]       # Pool mémoire GMP
├── png_save.[ch]       # Export PNG avec métadonnées
└── SDL_gfx*            # Primitives graphiques
```

## Types de fractales

### Vectorielles (Types 1-2)

| Type | Nom | Iter max |
|------|-----|----------|
| 1 | Von Koch | 8 |
| 2 | Dragon | 20 |

### Escape-time (Types 3-23)

| Type | Nom | Fonction | Formule | Domaine | Defaults |
|------|-----|----------|---------|---------|----------|
| 3 | Mandelbrot | `Mendelbrot_Iteration()` | z₀=seed, z_{n+1}=z²+c | [-2.5,1.5]×[-1.5,1.5] | iter=9370, zoom=8 |
| 4 | Julia | `Julia_Iteration()` | z₀=zPixel, z_{n+1}=z²+seed | [-2,2]×[-1.5,1.5] | seed=(0.36,-0.08), iter=6250 |
| 5 | Julia Sin | `JuliaSin_Iteration()` | z₀=zPixel, z_{n+1}=seed×sin(z) | [-π,π]×[-2,2] | seed=(1,0.1) |
| 6 | Newton | `Newton_Iteration()` | z_{n+1}=((p-1)z^p+1)/(pz^{p-1}) | - | p=8 |
| 7 | Phoenix | `Phoenix_Iteration()` | z_{n+1}=z²+p1+p2×y, y_{n+1}=z | [-2,2]×[-1.5,1.5] | p1=0.567, p2=-0.5 |
| 8 | Buffalo | `Buffalo_Iteration()` | z_{n+1}=\|Re(z²)\|+i\|Im(z²)\|+c | [-2.5,1.5]×[-2,2] | iter=9370 |
| 9 | Barnsley J | `Barnsley1j_Iteration()` | z_{n+1}=(z±1)×seed selon Re(z) | [-4,4]×[-3,3] | seed=(1.1,0.6) |
| 10 | Barnsley M | `Barnsley1m_Iteration()` | z_{n+1}=(z±1)×c selon Re(z) | [-3,3]×[-2,2] | - |
| 11 | Magnet J | `Magnet1j_Iteration()` | z_{n+1}=(z²+(seed-1))²/(2z+(seed-2)) | [-2,2]² | seed=(1.63,-0.31) |
| 12 | Magnet M | `Magnet1m_Iteration()` | z_{n+1}=(z²+(c-1))²/(2z+(c-2)) | [-3,2]×[-2,2] | - |
| 13 | Burning Ship | `BurningShip_Iteration()` | z_{n+1}=(\|Re(z)\|+i\|Im(z)\|)²+c | [-2.5,1.5]×[-2,2] | - |
| 14 | Tricorn | `Tricorn_Iteration()` | z_{n+1}=conj(z)²+c | [-2.5,1.5]×[-1.5,1.5] | - |
| 15 | Mandelbulb | `Mandelbulb_Iteration()` | z_{n+1}=z⁸+c | [-1.5,1.5]² | - |
| 16 | Buddhabrot | `Buddhabrot_Draw()` | Densité trajectoires z²+c | [-2.5,1.5]×[-1.5,1.5] | iter=220 |
| 17 | Lyapunov | `Lyapunov_Draw()` | x_{n+1}=r×x(1-x), séq "BBBBBBAAAAAA" | [2.5,3.4]×[3.4,4] | iter=2000 |
| 18 | Perp. Burning Ship | `PerpendicularBurningShip_Iteration()` | z_{n+1}=(Re(z)-i\|Im(z)\|)²+c | [-2.5,1.5]×[-1.5,1.5] | iter=5000 |
| 19 | Celtic | `Celtic_Iteration()` | z_{n+1}=\|Re(z²)\|+i×Im(z²)+c | [-2,1]×[-1.5,1.5] | iter=5000 |
| 20 | Alpha Mandelbrot | `AlphaMandelbrot_Iteration()` | z_{n+1}=z²+(z²+c)²+c | [-2.5,1.5]×[-1.5,1.5] | iter=2000 |
| 21 | Pickover Stalks | `PickoverStalks_Iteration()` | z²+c avec trap min(\|Re\|,\|Im\|) | [-2,1]×[-1.5,1.5] | bailout=100, iter=1000 |
| 22 | Nova | `Nova_Iteration()` | z_{n+1}=z-a(z³-1)/(3z²)+c | [-3,3]×[-2,2] | bailout=20, iter=500 |
| 23 | Multibrot | `Multibrot_Iteration()` | z_{n+1}=z^{2.5}+c | [-2.5,1.5]×[-1.5,1.5] | iter=5000 |

**Notes** : Types 9-12, 17-23 accessibles uniquement via GUI. Bailout par défaut = 4.

## Palettes de couleurs (9)

| # | Nom | Gradient |
|---|-----|----------|
| 0 | SmoothFire | Noir → Rouge → Jaune → Blanc |
| 1 | SmoothOcean | Noir → Bleu → Cyan → Blanc |
| 2 | SmoothForest | Noir → Vert → Jaune/Vert → Blanc |
| 3 | SmoothViolet | Noir → Violet → Rose → Blanc |
| 4 | SmoothRainbow | Arc-en-ciel complet |
| 5 | SmoothSunset | Noir → Orange → Rouge → Violet → Bleu |
| 6 | SmoothPlasma | Bleu → Violet → Rose → Orange (défaut) |
| 7 | SmoothIce | Blanc → Cyan → Bleu → Noir |
| 8 | SmoothCosmic | Noir → Bleu nuit → Turquoise → Jaune → Orange → Rouge |

Interpolation continue avec alternance endroit/envers. `colorRepeat` : 40 (escape-time), 2 (Lyapunov).

## Structure principale

```c
typedef struct fractal_struct {
  int xpixel, ypixel;           // Dimensions pixels
  double xmin, xmax, ymin, ymax;// Plan complexe
  complex seed;                 // Paramètre Julia
  int iterationMax, bailout;    // Limites
  int zoomfactor, type;         // Zoom (2-8), type (1-23)
  int colorMode, colorRepeat;   // Palette (0-8), répétitions (2-40)
  int cmatrix_valid;            // Cache colorisation valide
  int *fmatrix;                 // Matrice itérations
  complex *zmatrix;             // Valeurs z finales
  color *cmatrix;               // Couleurs calculées
  fractal_cache cache;          // Cache zoom
#ifdef HAVE_GMP
  int use_gmp;                  // Mode GMP actif
  mp_bitcnt_t gmp_precision;    // Précision bits
  // ... structures GMP pré-allouées
#endif
} fractal;
```

## Optimisations

### Divergence Detection (DDp1/DDp2)
Rendu en 2 passes : grille 2×2 puis affinage. Évite calcul si 4 voisins identiques.

### OpenMP
Parallélisation multi-cœurs avec `#pragma omp parallel for`, ordonnancement `guided`. Incréments atomiques pour Buddhabrot.

### SIMD
- **SSE4.1** : 2 complexes/op (128 bits)
- **AVX** : 4 complexes/op (256 bits)
- Fonctions : `complex_mul_sse4()`, `complex_add_sse4()`, `complex_mag2_sse4()`, `complex_mul_avx()`

### GMP
Précision arbitraire automatique pour zooms profonds. Pool mémoire (`gmp_pool.[ch]`) et contextes pré-alloués.

## Export PNG

Métadonnées incluses : FractalType, XMin/XMax/YMin/YMax, Iterations, Bailout, ColorMode, ColorRepeat, JuliaSeed (types 4-5), dimensions, Software.

Fonctions : `SavePNGWithMetadata()`, `LoadPNGMetadata()`

## Barre d'état

`Type | Rep:N | Palette | Iter:N | Zoom | Centre (x,y) | Temps (ms) | [GMP N bits]`

## Flux d'exécution

```
main() → SDL_Init() → SDLGUI_Init() → Fractal_Init()
       → Boucle: SDL_WaitEvent() → EventCheck() → Clavier/Souris/Quit
```

## Fonctions clés

| Fonction | Fichier | Description |
|----------|---------|-------------|
| `Fractal_Draw()` | EscapeTime.c | Dispatch vers DDp1/DDp2 ou algorithmes spéciaux (16-17) |
| `Fractal_Init()` | EscapeTime.c | Initialisation avec defaults selon type |
| `Fractal_GetTypeName()` | EscapeTime.c | Nom du type (0-23) |
| `Buddhabrot_Draw()` | EscapeTime.c | Algorithme densité trajectoires |
| `Lyapunov_Draw()` | EscapeTime.c | Exposant de Lyapunov |
| `GetColorForIteration()` | colorization.c | Colorisation unifiée |
| `SDLGUI_Init()` | SDLGUI.c | Initialisation GUI 23 boutons |

## Notes techniques

- GUI : 51px haut, barre d'état : 20px bas
- Defines : `HAVE_OPENMP`, `HAVE_SSE4_1`, `HAVE_AVX`, `HAVE_AVX2`, `HAVE_PNG`, `HAVE_GMP`
- Zoom sélection : `main.c` (SDL_MOUSEBUTTONDOWN/MOTION/UP)

## TODO

- Historique zoom (undo/redo)
- Palettes personnalisables
- GUI : tooltips, indicateur sélection, panneau latéral
- Fractales 3D/Quaternions (architecture majeure requise)
