 /* ExcapeTime.h
  fractal formula header
  released under GPL2
  Copyleft 2001-2003 VERHILLE Arnaud
*/

#ifndef ESCAPETIME_H
#define ESCAPETIME_H

#include "config.h"
#include "SDL.h"
#include "complexmath.h"
#include "color_types.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#ifdef HAVE_GMP
#include <gmp.h>
#include "complexmath_gmp.h"
#include "gmp_pool.h"
#endif

// Structure de cache pour réutilisation lors de zooms
typedef struct {
  double xmin_cached, xmax_cached, ymin_cached, ymax_cached;
  int* fmatrix_cached;
  color* cmatrix_cached;
  int cache_valid;
  int cache_xpixel, cache_ypixel;  // Dimensions du cache
} fractal_cache;

// Create a new fractal type
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
  int colorRepeat;  // Nombre de répétitions du gradient de couleur (2-20, de 2 en 2)
  int last_colorRepeat;  // colorRepeat lors du dernier calcul de cmatrix
  double zoom_level;   // Niveau de zoom actuel pour détection de changement
  int *fmatrix;    // la matrice d'iteration
  complex *zmatrix;  // la matrice de la valeur de z a la derniere iteration
  color *cmatrix; // Une matrice de couleurs, soit la fractale finale.
  fractal_cache cache;  // Cache pour réutilisation lors de zooms
#ifdef HAVE_GMP
  int use_gmp;     // Flag pour utiliser GMP (1) ou double (0)
  mp_bitcnt_t gmp_precision; // Précision GMP en bits
  complex_gmp *zmatrix_gmp;  // Matrice GMP optionnelle
  mpf_t xmin_gmp, xmax_gmp, ymin_gmp, ymax_gmp;  // Coordonnées GMP pour précision
  gmp_iteration_context iteration_ctx;  // Contexte pré-alloué pour itérations
  gmp_mul_temps mul_temps;  // Temporaires de multiplication pré-alloués
#endif
} fractal;

// Fractal point Calcul result
typedef struct {
  int iteration;
  complex z;
} fractalresult;


// ******************
// Interface publique
// ******************

// Initialisation et destruction
 fractal Fractal_Init (int screenW, int screenH, int type);
 void Fractal_Destroy (fractal);

 // Calcul de la matrice de complexe z lors de leur derniere iteration
 void Fractal_CalculateMatrix (fractal*);
 void Fractal_CalculateMatrix_DDp1 (fractal*, SDL_Surface*, void*, int*, int, int, const char*);
 void Fractal_CalculateMatrix_DDp2 (fractal*, SDL_Surface*, void*, int*, int, int, const char*);

 // Calcul de la couleur
 void Fractal_CalculateColorMatrix (fractal*, SDL_Surface*, void*, int*, int, int); // Selecteur
 void FractalColorTest (fractal*);

// Color Formulae Utilities
color Fractal_ReadColorMatrix (fractal, int, int);
int Fractal_ReadColorMatrixRed (fractal, int, int);
int Fractal_ReadColorMatrixGreen (fractal, int, int);
int Fractal_ReadColorMatrixBlue (fractal, int, int);

Uint32 Fractal_Draw (SDL_Surface*, fractal, int, int, void* gui);
void Fractal_ChangeType (fractal* f, int type);

// Formulae Utilities
 fractalresult FormulaSelector (fractal f, complex zPixel);
 fractalresult Mendelbrot_Iteration (fractal, complex);
 fractalresult Julia_Iteration (fractal, complex);
 fractalresult JuliaSin_Iteration (fractal, complex);
 fractalresult Newton_Iteration (fractal, complex);
 fractalresult Phoenix_Iteration (fractal, complex);
fractalresult Barnsleyj1_Iteration (fractal, complex);
fractalresult Barnsleym1_Iteration (fractal, complex);
fractalresult BurningShip_Iteration (fractal, complex);
fractalresult Buffalo_Iteration (fractal, complex);
fractalresult Tricorn_Iteration (fractal, complex);
fractalresult Mandelbulb_Iteration (fractal, complex);
fractalresult PerpendicularBurningShip_Iteration (fractal, complex);
fractalresult Celtic_Iteration (fractal, complex);
fractalresult AlphaMandelbrot_Iteration (fractal, complex);
fractalresult PickoverStalks_Iteration (fractal, complex);
fractalresult Nova_Iteration (fractal, complex);
fractalresult Multibrot_Iteration (fractal, complex);

#ifdef HAVE_GMP
// GMP versions of iteration functions
fractalresult FormulaSelector_GMP (fractal f, complex_gmp zPixel);
fractalresult Mendelbrot_Iteration_GMP (fractal, complex_gmp);
fractalresult Julia_Iteration_GMP (fractal, complex_gmp);
fractalresult JuliaSin_Iteration_GMP (fractal, complex_gmp);
fractalresult Newton_Iteration_GMP (fractal, complex_gmp);
fractalresult Phoenix_Iteration_GMP (fractal, complex_gmp);
fractalresult Barnsleyj1_Iteration_GMP (fractal, complex_gmp);
fractalresult Barnsleym1_Iteration_GMP (fractal, complex_gmp);
fractalresult BurningShip_Iteration_GMP (fractal, complex_gmp);
fractalresult Buffalo_Iteration_GMP (fractal, complex_gmp);
fractalresult Tricorn_Iteration_GMP (fractal, complex_gmp);
fractalresult Mandelbulb_Iteration_GMP (fractal, complex_gmp);
fractalresult PerpendicularBurningShip_Iteration_GMP (fractal, complex_gmp);
fractalresult Celtic_Iteration_GMP (fractal, complex_gmp);
fractalresult AlphaMandelbrot_Iteration_GMP (fractal, complex_gmp);
fractalresult PickoverStalks_Iteration_GMP (fractal, complex_gmp);
fractalresult Nova_Iteration_GMP (fractal, complex_gmp);
fractalresult Multibrot_Iteration_GMP (fractal, complex_gmp);
#endif

 // Fractal Definition
 void Mendelbrot_def (fractal* f);
 void Julia_def (fractal* f);
 void JuliaSin_def (fractal* f);
 void Newton_def (fractal* f);
 void Phoenix_def (fractal* f);
 void Barnsley1j_def (fractal* f);
void Barnsley1m_def (fractal* f);
void Magnet1j_def (fractal* f);
void Magnet1m_def (fractal* f);
void BurningShip_def (fractal* f);
void Buffalo_def (fractal* f);
void Tricorn_def (fractal* f);
void Mandelbulb_def (fractal* f);
void PerpendicularBurningShip_def (fractal* f);
void Celtic_def (fractal* f);
void AlphaMandelbrot_def (fractal* f);
void PickoverStalks_def (fractal* f);
void Nova_def (fractal* f);
void Multibrot_def (fractal* f);
void Buddhabrot_def (fractal* f);
void Lyapunov_def (fractal* f);

// Buddhabrot special draw function (density algorithm)
// gui parameter can be NULL if no GUI progress display needed
Uint32 Buddhabrot_Draw (SDL_Surface*, fractal*, int, int, void* gui);

// Lyapunov special draw function (exponent-based coloring)
Uint32 Lyapunov_Draw (SDL_Surface*, fractal*, int, int, void* gui);

// Utilitaire pour obtenir le nom de la fractale selon son type
const char* Fractal_GetTypeName(int type);

#endif /* ESCAPETIME_H */

