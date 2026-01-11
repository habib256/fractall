# fractall

Visualiseur de fractales portable en C utilisant SDL.

**Licence** : GPL-2.0
**Auteur** : Arnaud VERHILLE (2001-2003)
**Version** : 0.5

## Compilation

```bash
./autogen.sh      # Si configure n'existe pas
./configure       # --with-gmp pour activer la haute précision
make
./src/fractall
```

### Dépendances

- SDL 1.2.0+
- Bibliothèque mathématique (`-lm`)
- GMP (optionnel) : précision arbitraire pour zooms profonds

## Utilisation

```bash
fractall [OPTIONS]
```

| Option | Description | Défaut |
|--------|-------------|--------|
| `-x<n>` | Largeur fenêtre | 800 |
| `-y<n>` | Hauteur fenêtre | 600 |
| `-f` | Plein écran | - |
| `-g<n>` | Hauteur GUI | 51 |
| `-nogui` | Sans interface | - |
| `-h` | Aide | - |

### Contrôles

| Touche | Effet |
|--------|-------|
| **F1-F8** | Von Koch, Dragon, Mandelbrot, Julia, Julia Sin, Newton, Phoenix, Sierpinski |
| **F9-F12** | Burning Ship, Tricorn, Mandelbulb, Buddhabrot |
| **Clic gauche** | Zoom / +1 itération (vectorielles) |
| **Clic droit** | Dézoom / -1 itération |
| **C** | Changer palette |
| **S** | Screenshot (Screenshot.bmp) |
| **Q/ESC** | Quitter |

## Architecture

```
src/
├── main.c           # Point d'entrée, événements
├── EscapeTime.[ch]  # Fractales escape-time
├── VonKoch.[ch]     # Fractales vectorielles
├── SDLGUI.[ch]      # Interface graphique
├── complexmath.[ch] # Arithmétique complexe
├── complexmath_gmp.[ch] # Arithmétique GMP (optionnel)
├── precision_detector.[ch] # Détection précision GMP
└── SDL_gfx*         # Primitives graphiques
```

## Types de fractales

### Vectorielles (récursives) - Types 1-2

| Type | Nom | Itérations max |
|------|-----|----------------|
| 1 | Von Koch | 8 |
| 2 | Dragon | 20 |

### Escape-time - Types 3-16

| Type | Nom | Formule |
|------|-----|---------|
| 3 | Mandelbrot | z(n+1) = z(n)² + c |
| 4 | Julia | z(n+1) = z(n)² + seed |
| 5 | Julia Sin | z(n+1) = c × sin(z(n)) |
| 6 | Newton | z(n+1) = ((p-1)×z^p + 1) / (p×z^(p-1)) |
| 7 | Phoenix | z(n+1) = z(n)² + p1 + p2×y(n) |
| 8 | Sierpinski | Transformation conditionnelle |
| 9 | Barnsley J | z(n+1) = (z±1)×c selon Re(z) |
| 10 | Barnsley M | Variante Mandelbrot de Barnsley |
| 11 | Magnet J | Formule magnétique (Julia) |
| 12 | Magnet M | Formule magnétique (Mandelbrot) |
| 13 | Burning Ship | z(n+1) = (\|Re(z)\| + i\|Im(z)\|)² + c |
| 14 | Tricorn | z(n+1) = conj(z)² + c |
| 15 | Mandelbulb | z(n+1) = z(n)⁸ + c |
| 16 | Buddhabrot | Densité des trajectoires d'échappement |

**Note** : F9-F12 mappent vers les types 13-16 (les types 9-12 sont accessibles via GUI).

### Buddhabrot

Algorithme de densité différent :
1. Échantillonnage aléatoire du plan complexe
2. Trace les trajectoires des points qui s'échappent
3. Incrémente un compteur de densité par pixel traversé
4. Couleur finale basée sur densité (échelle log)

## Palettes de couleurs

2 palettes disponibles (touche **C**) :

| Mode | Description |
|------|-------------|
| SmoothFire | Noir → Rouge → Jaune → Blanc (dégradés fluides, répétition 4×) |
| SmoothOcean | Noir → Bleu → Cyan → Blanc (dégradés fluides, répétition 4×) |

Les deux utilisent une interpolation continue basée sur |z| et alternent endroit/envers pour éviter les transitions brutales.

## Structure principale

### `fractal` (EscapeTime.h)

```c
typedef struct {
  int xpixel, ypixel;       // Dimensions en pixels
  double xmin, xmax;        // Bornes du plan complexe
  double ymin, ymax;
  complex seed;             // Paramètre c pour Julia
  int iterationMax;         // Limite d'itérations
  int bailout;              // Seuil d'échappement (généralement 4)
  int zoomfactor;           // Facteur de zoom (2-8 selon fractale)
  int type;                 // Type de fractale (1-16)
  int colorMode;            // 0=SmoothFire, 1=SmoothOcean
  int *fmatrix;             // Matrice d'itérations
  complex *zmatrix;         // Valeurs z finales
  color *cmatrix;           // Couleurs calculées
#ifdef HAVE_GMP
  int use_gmp;              // Utiliser GMP (1) ou double (0)
  mp_bitcnt_t gmp_precision; // Précision GMP en bits
  complex_gmp *zmatrix_gmp; // Matrice GMP optionnelle
  mpf_t xmin_gmp, xmax_gmp; // Coordonnées GMP
  mpf_t ymin_gmp, ymax_gmp;
#endif
} fractal;
```

## Optimisation du rendu

**Divergence Detection** en deux passes :

1. **DDp1** : Calcul grille 2×2, interpolation pixels intermédiaires
2. **DDp2** : Affinage des manquants, évite calcul si 4 voisins identiques

Permet un aperçu rapide suivi d'un rendu complet.

## Précision GMP

Avec `--with-gmp`, le programme détecte automatiquement quand la précision double (53 bits) devient insuffisante et bascule vers GMP. La précision augmente dynamiquement avec le niveau de zoom.

Fichiers clés :
- `precision_detector.[ch]` : Détection automatique du seuil
- `complexmath_gmp.[ch]` : Arithmétique complexe en précision arbitraire

## Barre d'état

Affiche : Type | Palette | Zoom | Centre (x,y) | Temps rendu (ms)

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
              ├── Clavier (F1-F12, C, S, Q)
              ├── Souris (zoom, GUI)
              └── Quit
```

## TODO

- Réutilisation de cmatrix lors du zoom
- Export PNG avec métadonnées
- Historique zoom (undo/redo)
- Palettes personnalisables

## Notes

- GUI : 51px en haut, barre d'état : 20px en bas
- `Buddhabrot_Draw()` pour le type 16
- `Fractal_GetTypeName()` retourne le nom selon le type
