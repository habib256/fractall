# fractall

Visualiseur de fractales portable en C utilisant SDL.

**Licence** : GPL-2.0
**Auteur** : Arnaud VERHILLE (2001-2003)
**Version** : 0.5

## Compilation

```bash
./autogen.sh      # Si configure n'existe pas
./configure
make
./src/fractall
```

### Dépendances

- SDL 1.2.0+
- Bibliothèque mathématique standard (`-lm`)

## Utilisation

```bash
fractall [OPTIONS]
```

| Option | Description | Défaut |
|--------|-------------|--------|
| `-x<n>` | Largeur de la fenêtre | 800 |
| `-y<n>` | Hauteur de la fenêtre | 600 |
| `-f` | Mode plein écran | - |
| `-g<n>` | Hauteur du GUI | 51 |
| `-nogui` | Désactiver l'interface graphique | - |
| `-h`, `-help` | Afficher l'aide | - |

### Contrôles

| Touche/Action | Effet |
|---------------|-------|
| **F1** | Von Koch |
| **F2** | Dragon Fractal |
| **F3** | Mandelbrot |
| **F4** | Julia |
| **F5** | Julia Sin |
| **F6** | Newton |
| **F7** | Phoenix |
| **F8** | Sierpinski |
| **F9** | Burning Ship |
| **F10** | Tricorn |
| **F11** | Mandelbulb (2D, puissance 8) |
| **F12** | Buddhabrot |
| **Clic gauche** | Zoom / +1 itération (vectorielles) |
| **Clic droit** | Dézoom / -1 itération (vectorielles) |
| **C** | Changer de palette de couleurs |
| **S** | Capture d'écran (Screenshot.bmp) |
| **Q** / **ESC** | Quitter |

## Architecture

```
fractall/
├── src/                    # Code source
│   ├── main.c              # Point d'entrée, boucle événements
│   ├── EscapeTime.[ch]     # Fractales temps d'échappement
│   ├── VonKoch.[ch]        # Fractales vectorielles
│   ├── SDLGUI.[ch]         # Interface utilisateur
│   ├── complexmath.[ch]    # Arithmétique complexe
│   ├── exception.[ch]      # Gestion d'exceptions
│   └── SDL_gfx*            # Primitives graphiques SDL
├── docs/                   # Documentation
├── configure.ac            # Configuration autotools
└── Makefile.am             # Build automake
```

## Types de fractales

### Fractales vectorielles (récursives)

| Type | Nom | Itérations max |
|------|-----|----------------|
| 1 | Von Koch | 8 |
| 2 | Dragon Fractal | 20 |

Implémentées dans `VonKoch.c` avec tracé récursif.

### Fractales temps d'échappement

| Type | Nom | Formule |
|------|-----|---------|
| 3 | Mandelbrot | z(n+1) = z(n)² + c |
| 4 | Julia | z(n+1) = z(n)² + seed |
| 5 | Julia Sin | z(n+1) = c * sin(z(n)) |
| 6 | Newton | z(n+1) = ((p-1)*z^p + 1) / (p*z^(p-1)) |
| 7 | Phoenix | z(n+1) = z(n)² + p1 + p2*y(n) |
| 8 | Sierpinski | Transformation conditionnelle |
| 9 | Barnsley J1 | z(n+1) = (z±1)*c selon Re(z) |
| 10 | Barnsley M1 | Variante Mandelbrot de Barnsley |
| 11 | Magnet J1 | Formule magnétique (Julia) |
| 12 | Magnet M1 | Formule magnétique (Mandelbrot) |
| 13 | Burning Ship | z(n+1) = (\|Re(z)\| + i\|Im(z)\|)² + c |
| 14 | Tricorn | z(n+1) = conj(z)² + c |
| 15 | Mandelbulb | z(n+1) = z(n)⁸ + c |
| 16 | Buddhabrot | Densité des trajectoires d'échappement |

### Buddhabrot (algorithme de densité)

La Buddhabrot utilise un algorithme différent des autres fractales :
1. Échantillonnage aléatoire de points dans le plan complexe
2. Pour chaque point qui s'échappe (n'est PAS dans Mandelbrot), on trace sa trajectoire
3. On incrémente un compteur de densité pour chaque pixel traversé
4. La couleur finale est basée sur la densité (échelle logarithmique)

Cet algorithme révèle la "probabilité" qu'un point soit traversé par une trajectoire d'échappement, créant une image ressemblant à un Bouddha méditant.

## Palettes de couleurs

7 palettes disponibles (touche **C** pour cycler) :

| Mode | Description |
|------|-------------|
| Normal | Gradient arc-en-ciel (R, G, B cycliques) |
| Mono | Niveaux de gris |
| Fire | Noir → Rouge → Jaune → Blanc |
| Ocean | Noir → Bleu → Cyan → Blanc |
| Rainbow | Arc-en-ciel HSV avec smooth coloring |
| SmoothFire | Fire avec dégradés fluides (sans bandes) |
| SmoothOcean | Ocean avec dégradés fluides (sans bandes) |

Les palettes "Smooth" utilisent une interpolation continue basée sur |z| pour éliminer les bandes de couleur visibles.

## Barre d'état

Affiche en bas de l'écran :
- Type de fractale
- Palette active
- Facteur de zoom
- Coordonnées du centre
- Temps de rendu (ms)

## Structures de données principales

### `fractal` (EscapeTime.h)

```c
typedef struct {
  int xpixel, ypixel;       // Dimensions en pixels
  double xmin, xmax;        // Bornes du plan complexe
  double ymin, ymax;
  complex seed;             // Paramètre c pour Julia
  int iterationMax;         // Limite d'itérations
  int bailout;              // Seuil d'échappement
  int zoomfactor;           // Facteur de zoom
  int type;                 // Type de fractale (3-16)
  int colorMode;            // Palette (0-6: Normal, Mono, Fire, Ocean, Rainbow, SmoothFire, SmoothOcean)
  int *fmatrix;             // Matrice d'itérations
  complex *zmatrix;         // Valeurs z finales
  color *cmatrix;           // Couleurs calculées
} fractal;
```

### `gui` (SDLGUI.h)

Interface avec barre de menu (16 boutons) et barre d'état.

## Algorithmes d'optimisation

Le rendu utilise la **Divergence Detection** en deux passes :

1. **DDp1** : Calcul sur grille 2x2, interpolation des pixels intermédiaires
2. **DDp2** : Affinage des pixels manquants, évite le calcul si les 4 voisins ont la même valeur d'itération

Cette technique permet un aperçu rapide suivi d'un rendu complet.

## Flux d'exécution

```
main()
  ├── Parsing arguments
  ├── SDL_Init()
  ├── SDLGUI_Init()
  ├── Fractal_Init()
  └── Boucle principale
        ├── SDL_WaitEvent()
        └── EventCheck()
              ├── Clavier (F1-F12, S, Q)
              ├── Souris (zoom, GUI)
              └── Quit
```

## Fichiers sources

| Fichier | Description | Origine |
|---------|-------------|---------|
| `main.c` | Point d'entrée, gestion événements | Projet |
| `EscapeTime.c` | Calcul des fractales escape-time | Projet |
| `VonKoch.c` | Fractales vectorielles | Projet |
| `SDLGUI.c` | Interface graphique | Projet |
| `complexmath.c` | Arithmétique complexe | Akihiro Sato |
| `exception.c` | Gestion d'exceptions | Équipe LOGIN |
| `SDL_gfxPrimitives.c` | Primitives graphiques | SDL_gfx Team |
| `SDL_imageFilter.c` | Filtres d'image | SDL_gfx Team |
| `SDL_rotozoom.c` | Rotation/zoom surfaces | SDL_gfx Team |

## TODO

- Réutilisation de cmatrix lors du zoom (optimisation)
- Export PNG avec métadonnées
- Historique zoom (undo/redo)
- Palettes personnalisables (fichiers de configuration)

## Notes de développement

- Les types 1-2 sont des fractales vectorielles, les types 3-16 sont des fractales escape-time
- Les raccourcis F9-F12 mappent vers les types internes 13-16
- Le type 16 (Buddhabrot) utilise un algorithme de densité différent via `Buddhabrot_Draw()`
- La taille par défaut de la fenêtre est 800x600
- Le GUI occupe 51 pixels en haut, la barre d'état 20 pixels en bas
