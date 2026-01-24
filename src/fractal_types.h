/* fractal_types.h
   Structures de donnees pour les fractales
   released under GPL2
   Copyleft 2001-2026 VERHILLE Arnaud
*/

#ifndef FRACTAL_TYPES_H
#define FRACTAL_TYPES_H

#include "config.h"
#include "complexmath.h"
#include "color_types.h"
#ifdef HAVE_GMP
#include <gmp.h>
#include "complexmath_gmp.h"
#include "gmp_pool.h"
#endif

// Structure de cache pour reutilisation lors de zooms
typedef struct {
  double xmin_cached, xmax_cached, ymin_cached, ymax_cached;
  int* fmatrix_cached;
  color* cmatrix_cached;
  int cache_valid;
  int cache_xpixel, cache_ypixel;  // Dimensions du cache
} fractal_cache;

// Structure principale de fractale
typedef struct fractal_struct {
  int xpixel, ypixel;
  double xmin, xmax, ymin, ymax;
  complex seed;
  int iterationMax;
  int bailout;
  int zoomfactor;
  int type;
  int colorMode;   // 0-8 for different palettes (see colorization.h)
  int cmatrix_valid;   // 1 si cmatrix est valide pour le colorMode actuel
  int last_colorMode;  // colorMode lors du dernier calcul de cmatrix
  int colorRepeat;  // Nombre de repetitions du gradient de couleur (2-40, de 2 en 2)
  int last_colorRepeat;  // colorRepeat lors du dernier calcul de cmatrix
  double zoom_level;   // Niveau de zoom actuel pour detection de changement
  int *fmatrix;    // la matrice d'iteration
  complex *zmatrix;  // la matrice de la valeur de z a la derniere iteration
  color *cmatrix; // Une matrice de couleurs, soit la fractale finale.
  fractal_cache cache;  // Cache pour reutilisation lors de zooms
#ifdef HAVE_GMP
  int use_gmp;     // Flag pour utiliser GMP (1) ou double (0)
  mp_bitcnt_t gmp_precision; // Precision GMP en bits
  complex_gmp *zmatrix_gmp;  // Matrice GMP optionnelle
  mpf_t xmin_gmp, xmax_gmp, ymin_gmp, ymax_gmp;  // Coordonnees GMP pour precision
  gmp_iteration_context iteration_ctx;  // Contexte pre-alloue pour iterations
  gmp_mul_temps mul_temps;  // Temporaires de multiplication pre-alloues
#endif
} fractal;

// Resultat du calcul d'un point de fractale
typedef struct {
  int iteration;
  complex z;
} fractalresult;

#endif /* FRACTAL_TYPES_H */
