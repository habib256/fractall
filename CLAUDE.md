# fractall

Visualiseur de fractales portable en C utilisant SDL.

**Licence** : GPL-2.0
**Auteur** : Arnaud VERHILLE (2001-2003)
**Version** : 1.0

## Compilation

```bash
./autogen.sh      # Si configure n'existe pas
./configure       # Options disponibles ci-dessous
make
./src/fractall
```

### Options de configuration

| Option | Description | Défaut |
|--------|-------------|--------|
| `--with-gmp` | Active la haute précision (zooms profonds) | activé |
| `--without-gmp` | Désactive GMP | - |
| `--disable-openmp` | Désactive le calcul parallèle OpenMP | activé |
| `--disable-simd` | Désactive les optimisations SIMD | activé |

### Dépendances

- SDL 1.2.0+
- Bibliothèque mathématique (`-lm`)
- GMP (optionnel) : précision arbitraire pour zooms profonds
- OpenMP (optionnel) : calcul parallèle multi-cœurs

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
| **F1-F7** | Von Koch, Dragon, Mandelbrot, Julia, Julia Sin, Newton, Phoenix |
| **F9-F12** | Burning Ship, Tricorn, Mandelbulb, Buddhabrot |
| **GUI** | 17 boutons pour sélectionner toutes les fractales (types 1-17) |
| **Clic gauche** | Zoom / +1 itération (vectorielles) |
| **Clic droit** | Dézoom / -1 itération |
| **C** | Changer palette |
| **R** | Changer nombre de répétitions du gradient (2, 4, 6, 8, 10, 12, 14, 16, 18, 20) |
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
├── complexmath_simd.c   # Arithmétique SIMD (SSE4.1/AVX)
├── precision_detector.[ch] # Détection précision GMP
└── SDL_gfx*         # Primitives graphiques
```

## Types de fractales

### Vectorielles (récursives) - Types 1-2

| Type | Nom | Itérations max |
|------|-----|----------------|
| 1 | Von Koch | 8 |
| 2 | Dragon | 20 |

### Escape-time - Types 3-17

| Type | Nom | Formule |
|------|-----|---------|
| 3 | Mandelbrot | z(n+1) = z(n)² + c |
| 4 | Julia | z(n+1) = z(n)² + seed |
| 5 | Julia Sin | z(n+1) = c × sin(z(n)) |
| 6 | Newton | z(n+1) = ((p-1)×z^p + 1) / (p×z^(p-1)) |
| 7 | Phoenix | z(n+1) = z(n)² + p1 + p2×y(n) |
| 8 | *(supprimé)* | - |
| 9 | Barnsley J | z(n+1) = (z±1)×c selon Re(z) |
| 10 | Barnsley M | Variante Mandelbrot de Barnsley |
| 11 | Magnet J | Formule magnétique (Julia) |
| 12 | Magnet M | Formule magnétique (Mandelbrot) |
| 13 | Burning Ship | z(n+1) = (\|Re(z)\| + i\|Im(z)\|)² + c |
| 14 | Tricorn | z(n+1) = conj(z)² + c |
| 15 | Mandelbulb | z(n+1) = z(n)⁸ + c |
| 16 | Buddhabrot | Densité des trajectoires d'échappement |
| 17 | Lyapunov Zircon City | Exposant de Lyapunov (séquence "BBBBBBAAAAAA", domaine [2.5, 3.4] × [3.4, 4.0]) |

**Note** : F9-F12 mappent vers les types 13-16 (les types 9-12 et 17 sont accessibles via GUI).

### Buddhabrot (Type 16)

Algorithme de densité différent :
1. Échantillonnage aléatoire du plan complexe
2. Trace les trajectoires des points qui s'échappent
3. Incrémente un compteur de densité par pixel traversé
4. Couleur finale basée sur densité (échelle log)

### Lyapunov (Type 17)

Algorithme basé sur l'exposant de Lyapunov de la suite logistique x_{n+1} = r_n × x_n × (1 - x_n) :

**Principe** :
1. Le paramètre r alterne entre a (axe X) et b (axe Y) selon la séquence "BBBBBBAAAAAA"
2. Phase de stabilisation (warmup) de 50 itérations
3. Calcul de l'exposant sur 2000 itérations par défaut
4. Coloration avec palettes (SmoothFire par défaut, 2 répétitions par défaut)
5. Domaine : [2.5, 3.4] × [3.4, 4.0] (image de référence classique "Zircon City")

**Note** : Les variantes Lyapunov 18-27 ont été supprimées. Seule Zircon City (type 17) est conservée.

## Palettes de couleurs

5 palettes disponibles (touche **C**) :

| Mode | Description |
|------|-------------|
| SmoothFire (0) | Noir → Rouge → Jaune → Blanc (dégradés fluides, répétition 4×) |
| SmoothOcean (1) | Noir → Bleu → Cyan → Blanc (dégradés fluides, répétition 4×) |
| SmoothForest (2) | Noir → Vert → Jaune → Blanc (dégradés fluides, répétition 4×) |
| SmoothViolet (3) | Noir → Violet → Rose → Blanc (dégradés fluides, répétition 4×) |
| SmoothRainbow (4) | Arc-en-ciel complet (Rouge → Orange → Jaune → Vert → Cyan → Bleu → Violet) |

Toutes utilisent une interpolation continue basée sur |z| et alternent endroit/envers pour éviter les transitions brutales.

**Note** : Buddhabrot a son propre algorithme de coloration et n'est pas affecté par le changement de palette. La fractale de Lyapunov (type 17) supporte maintenant toutes les palettes disponibles (touche C).

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
  int type;                 // Type de fractale (1-17)
  int colorMode;            // 0=SmoothFire, 1=SmoothOcean, 2=SmoothForest, 3=SmoothViolet, 4=SmoothRainbow
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

### Divergence Detection

Algorithme en deux passes :

1. **DDp1** : Calcul grille 2×2, interpolation pixels intermédiaires
2. **DDp2** : Affinage des manquants, évite calcul si 4 voisins identiques

Permet un aperçu rapide suivi d'un rendu complet.

### Parallélisation OpenMP

Avec OpenMP activé, le calcul des fractales est parallélisé sur tous les cœurs CPU :
- Boucles de rendu DDp1/DDp2 parallélisées avec `#pragma omp parallel for`
- Ordonnancement `guided` pour équilibrage dynamique de charge
- Variables thread-local pour éviter les race conditions
- Support spécial pour Buddhabrot (incréments atomiques sur la matrice de densité)
- Support spécial pour Lyapunov (parallélisation par lignes)

### Optimisations SIMD

Opérations vectorisées sur nombres complexes (`complexmath_simd.c`) :
- **SSE4.1** : 2 complexes traités simultanément (128 bits)
- **AVX** : 4 complexes traités simultanément (256 bits)
- Fonctions : `complex_mul_sse4()`, `complex_add_sse4()`, `complex_mag2_sse4()`, `complex_mul_avx()`

La détection des instructions SIMD est automatique à la compilation.

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
- 17 boutons dans le GUI (types 1-17)
- `Buddhabrot_Draw()` pour le type 16 (algorithme de densité)
- `Lyapunov_Draw()` pour le type 17 (algorithme d'exposant de Lyapunov - Zircon City)
- `Fractal_Draw()` détecte automatiquement les types 16-17 et appelle leurs fonctions spécialisées
- `Fractal_GetTypeName()` retourne le nom selon le type (0-17)
- OpenMP : `HAVE_OPENMP` défini dans `config.h` si disponible
- SIMD : `HAVE_SSE4_1`, `HAVE_AVX`, `HAVE_AVX2` définis selon le support CPU
