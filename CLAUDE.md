# fractall

Visualiseur de fractales portable en C utilisant SDL.

**Licence** : GPL-2.0
**Auteur** : Arnaud VERHILLE (2001-200)
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
| `-x<n>` | Largeur fenêtre | 1024 |
| `-y<n>` | Hauteur fenêtre | 768 |
| `-f` | Plein écran | - |
| `-g<n>` | Hauteur GUI | 51 |
| `-nogui` | Sans interface | - |
| `-h` | Aide | - |

### Contrôles

| Touche | Effet |
|--------|-------|
| **F1** | Von Koch |
| **F2** | Dragon |
| **F3** | Mandelbrot |
| **F4** | Julia |
| **F5** | Julia Sin |
| **F6** | Newton |
| **F7** | Phoenix |
| **F8** | Buffalo |
| **F9** | Burning Ship |
| **F10** | Tricorn |
| **F11** | Mandelbulb |
| **F12** | Buddhabrot |
| **GUI** | 23 boutons pour sélectionner toutes les fractales (types 1-23) |
| **Clic gauche** | Zoom / +1 itération (vectorielles) |
| **Clic droit** | Dézoom / -1 itération |
| **C** | Changer palette |
| **R** | Changer nombre de répétitions du gradient (2, 4, 6, 8, 10, 12, 14, 16, 18, 20) |
| **S** | Screenshot (Screenshot.bmp) |
| **Q/ESC** | Quitter |

**Note** : Les types 9-12 (Barnsley, Magnet), 17 (Lyapunov) et 18-23 (Perpendicular Burning Ship, Celtic, Alpha Mandelbrot, Pickover Stalks, Nova, Multibrot) sont accessibles uniquement via les boutons GUI.

## Architecture

```
src/
├── main.c           # Point d'entrée, événements
├── EscapeTime.[ch]  # Fractales escape-time (18 types)
├── colorization.[ch] # Système unifié de colorisation (8 palettes)
├── color_types.h    # Type color unifié
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

### Escape-time - Types 3-23

| Type | Nom | Formule |
|------|-----|---------|
| 3 | Mandelbrot | z(n+1) = z(n)² + c |
| 4 | Julia | z(n+1) = z(n)² + seed |
| 5 | Julia Sin | z(n+1) = c × sin(z(n)) |
| 6 | Newton | z(n+1) = ((p-1)×z^p + 1) / (p×z^(p-1)) |
| 7 | Phoenix | z(n+1) = z(n)² + p1 + p2×y(n) |
| 8 | Buffalo | z(n+1) = abs(Re(z²)) + i×abs(Im(z²)) + c |
| 9 | Barnsley J | z(n+1) = (z±1)×c selon Re(z) |
| 10 | Barnsley M | Variante Mandelbrot de Barnsley |
| 11 | Magnet J | Formule magnétique (Julia) |
| 12 | Magnet M | Formule magnétique (Mandelbrot) |
| 13 | Burning Ship | z(n+1) = (\|Re(z)\| + i\|Im(z)\|)² + c |
| 14 | Tricorn | z(n+1) = conj(z)² + c |
| 15 | Mandelbulb | z(n+1) = z(n)⁸ + c |
| 16 | Buddhabrot | Densité des trajectoires d'échappement |
| 17 | Lyapunov Zircon City | Exposant de Lyapunov (séquence "BBBBBBAAAAAA", domaine [2.5, 3.4] × [3.4, 4.0]) |
| 18 | Perpendicular Burning Ship | z(n+1) = (Re(z) - i×\|Im(z)\|)² + c |
| 19 | Celtic | z(n+1) = \|Re(z²)\| + i×Im(z²) + c |
| 20 | Alpha Mandelbrot | z(n+1) = z² + (z² + c)² + c |
| 21 | Pickover Stalks | z(n+1) = z² + c avec orbit trap min(\|Re(z)\|, \|Im(z)\|) |
| 22 | Nova | z(n+1) = z - a·(p(z)/p'(z)) + c (p(z) = z³ - 1) |
| 23 | Multibrot | z(n+1) = z^d + c (d réel non-entier, ex: 2.5) |

### Buffalo (Type 8)

Variation du Burning Ship avec `abs()` appliqué aux deux parties après le carré :
- Formule : `z(n+1) = abs(Re(z²)) + i×abs(Im(z²)) + c`
- Produit des formes symétriques ressemblant à un bison
- Accessible via **F8** ou bouton GUI

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
4. Coloration avec palettes (SmoothPlasma par défaut, 2 répétitions par défaut)
5. Domaine : [2.5, 3.4] × [3.4, 4.0] (image de référence classique "Zircon City")

**Note** : Les variantes Lyapunov 18-27 ont été supprimées. Seule Zircon City (type 17) est conservée.

## Palettes de couleurs

8 palettes disponibles (touche **C**) - système unifié de colorisation (`colorization.c`) :

| Mode | Description |
|------|-------------|
| SmoothFire (0) | Noir → Rouge → Jaune → Blanc (dégradés fluides, répétition configurable) |
| SmoothOcean (1) | Noir → Bleu → Cyan → Blanc (dégradés fluides, répétition configurable) |
| SmoothForest (2) | Noir → Vert foncé → Jaune/Vert clair → Blanc |
| SmoothViolet (3) | Noir → Violet foncé → Rose/Magenta → Blanc |
| SmoothRainbow (4) | Arc-en-ciel complet (Rouge → Orange → Jaune → Vert → Cyan → Bleu → Violet) |
| SmoothSunset (5) | Noir → Orange → Rouge → Violet → Bleu foncé |
| SmoothPlasma (6) | **NOUVEAU** - Bleu profond → Violet → Rose/Corail → Jaune/Orange |
| SmoothIce (7) | **NOUVEAU** - Blanc → Cyan clair → Bleu profond → Noir |

Toutes utilisent une interpolation continue basée sur des tables de gradient et alternent endroit/envers pour éviter les transitions brutales.

**Note** : Buddhabrot et Lyapunov supportent maintenant toutes les 8 palettes disponibles (touche C).

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
  int colorMode;            // 0=SmoothFire, 1=SmoothOcean, 2=SmoothForest, 3=SmoothViolet, 4=SmoothRainbow, 5=SmoothSunset, 6=SmoothPlasma (défaut), 7=SmoothIce
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

- Export PNG avec métadonnées
- Historique zoom (undo/redo)
- Palettes personnalisables
- Amélioration du GUI (voir section ci-dessous)
- Nouveaux types de fractales (voir PLAN_RECHERCHE_FRACTALES.md)

## Améliorations possibles : Nouveaux types de fractales

Un plan de recherche détaillé est disponible dans **PLAN_RECHERCHE_FRACTALES.md**.

### Catégories identifiées

#### Priorité HAUTE (faciles à implémenter) ✅ COMPLÉTÉ
1. **Perpendicular Burning Ship** (Type 18) ✅ - Variante du Burning Ship avec pliage perpendiculaire
   - Formule : `z(n+1) = (Re(z) - i×|Im(z)|)² + c`
   - Complexité : ⭐ Faible (similaire au type 13)

2. **Celtic Fractal** (Type 19) ✅ - Formes celtiques entrelacées
   - Formule : `z(n+1) = |Re(z²)| + i×Im(z²) + c`
   - Complexité : ⭐ Faible

3. **Alpha Mandelbrot** (Type 20) ✅ - Structures superposées
   - Formule : `z(n+1) = z² + (z² + c)² + c`
   - Complexité : ⭐ Faible

#### Priorité MOYENNE (intéressantes mais plus complexes) ✅ COMPLÉTÉ
4. **Pickover Stalks / Biomorphs** (Type 21) ✅ - Formes biologiques/organiques
   - Formule : `z(n+1) = z² + c` avec orbit trap `min(|Re(z)|, |Im(z)|)`
   - Complexité : ⭐⭐ Moyenne (système d'orbit trap)

5. **Nova Fractal** (Type 22) ✅ - Variante Newton avec spirales élégantes
   - Formule : `z(n+1) = z - a·(p(z)/p'(z)) + c` (p(z) = z³ - 1)
   - Complexité : ⭐⭐ Moyenne (dérivée polynomiale)

6. **Multibrot (puissances non-entières)** (Type 23) ✅ - Morphing entre formes
   - Formule : `z(n+1) = z^d + c` (d réel, ex: 2.5)
   - Complexité : ⭐⭐ Moyenne (gestion branch cuts)
   - Complexité : ⭐⭐ Moyenne (dérivée polynomiale)

6. **Multibrot (puissances non-entières)** - Morphing entre formes
   - Formule : `z(n+1) = z(n)^d + c` où d est réel (2.5, 3.7, etc.)
   - Complexité : ⭐⭐ Moyenne (gestion branch cuts)

#### Priorité BASSE (nécessitent modifications architecturales)
7. **Fractales 3D/Quaternions** - Mandelbox, Quaternion Julia
   - Complexité : ⭐⭐⭐ Élevée (rendu 3D requis)

### Ressources de recherche

- **FractalForums.com** : Base de données de formules
- **UltraFractal.com** : Bibliothèque de formules documentées
- **Paul Bourke** (paulbourke.net) : Articles et algorithmes
- **GitHub** : Implémentations open-source (kf2, Fractalshades)

Voir **PLAN_RECHERCHE_FRACTALES.md** pour le plan détaillé de recherche et d'implémentation.

## Notes

- GUI : 51px en haut, barre d'état : 20px en bas
- 23 boutons dans le GUI (types 1-23)
- `Buddhabrot_Draw()` pour le type 16 (algorithme de densité)
- `Lyapunov_Draw()` pour le type 17 (algorithme d'exposant de Lyapunov - Zircon City)
- `Fractal_Draw()` détecte automatiquement les types 16-17 et appelle leurs fonctions spécialisées
- `Fractal_GetTypeName()` retourne le nom selon le type (0-23)
- OpenMP : `HAVE_OPENMP` défini dans `config.h` si disponible
- SIMD : `HAVE_SSE4_1`, `HAVE_AVX`, `HAVE_AVX2` définis selon le support CPU

## Propositions d'amélioration du GUI

### État actuel

Le GUI actuel (`SDLGUI.c`) est minimaliste :
- Barre de 17 boutons carrés avec aperçu miniature de chaque fractale
- Barre de défilement horizontale (car tous les boutons ne rentrent pas)
- Barre d'état en bas affichant : Type | Palette | Zoom | Centre | Temps

### Améliorations proposées

#### 1. Tooltips sur les boutons
Afficher le nom de la fractale au survol du bouton (hover).
```c
// Nouvelle fonction dans SDLGUI.c
void SDLGUI_DrawTooltip(SDL_Surface* screen, int x, int y, const char* text);
```

#### 2. Indicateur de sélection
Mettre en évidence le bouton de la fractale actuellement sélectionnée (bordure colorée ou surbrillance).

#### 3. Panneau de contrôle latéral (optionnel)
Ajouter un panneau rétractable avec :
- Sélecteur de palette (liste déroulante ou boutons colorés)
- Slider pour colorRepeat (2-20)
- Affichage des coordonnées en temps réel
- Bouton Reset (retour aux coordonnées par défaut)

#### 4. Minimap de navigation
Petite fenêtre montrant la vue globale avec un rectangle indiquant la zone zoomée actuelle.

#### 5. Barre d'outils secondaire
Ajouter des boutons pour les actions fréquentes :
- [C] Palette suivante
- [R] Répétitions
- [S] Screenshot
- [←] Undo zoom
- [→] Redo zoom

#### 6. Raccourcis clavier affichés
Afficher les raccourcis sur les boutons ou dans les tooltips (ex: "F3" sur le bouton Mandelbrot).

#### 7. Mode plein écran amélioré
Cacher le GUI en mode plein écran avec possibilité de le faire apparaître au survol du haut de l'écran.

### Priorité suggérée

1. **Facile** : Indicateur de sélection, raccourcis affichés
2. **Moyen** : Tooltips, barre d'outils secondaire
3. **Complexe** : Panneau latéral, minimap, historique zoom
