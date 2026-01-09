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
#ifdef HAVE_GMP
#include <gmp.h>
#include "complexmath_gmp.h"
#endif

// A portable color type.
typedef struct {
  int r, g, b, a;
} color;

// Create a new fractal type
typedef struct {
  int xpixel, ypixel;
  double xmin, xmax, ymin, ymax;
  complex seed;
  int iterationMax;
  int bailout;
  int zoomfactor;
  int type;
  int colorMode;   // 0=Normal, 1=Monochrome, 2=Fire, 3=Ocean
  int *fmatrix;    // la matrice d'iteration
  complex *zmatrix;  // la matrice de la valeur de z a la derniere iteration
  color *cmatrix; // Une matrice de couleurs, soit la fractale finale.
#ifdef HAVE_GMP
  int use_gmp;     // Flag pour utiliser GMP (1) ou double (0)
  mp_bitcnt_t gmp_precision; // Précision GMP en bits
  complex_gmp *zmatrix_gmp;  // Matrice GMP optionnelle
  mpf_t xmin_gmp, xmax_gmp, ymin_gmp, ymax_gmp;  // Coordonnées GMP pour précision
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
 void FractalColorMonochrome (fractal*);
 void FractalColorNormal (fractal*);
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
 fractalresult Sierpinski_Iteration (fractal, complex);
fractalresult Barnsleyj1_Iteration (fractal, complex);
fractalresult Barnsleym1_Iteration (fractal, complex);
fractalresult BurningShip_Iteration (fractal, complex);
fractalresult Tricorn_Iteration (fractal, complex);
fractalresult Mandelbulb_Iteration (fractal, complex);

#ifdef HAVE_GMP
// GMP versions of iteration functions
fractalresult FormulaSelector_GMP (fractal f, complex_gmp zPixel);
fractalresult Mendelbrot_Iteration_GMP (fractal, complex_gmp);
fractalresult Julia_Iteration_GMP (fractal, complex_gmp);
fractalresult JuliaSin_Iteration_GMP (fractal, complex_gmp);
fractalresult Newton_Iteration_GMP (fractal, complex_gmp);
fractalresult Phoenix_Iteration_GMP (fractal, complex_gmp);
fractalresult Sierpinski_Iteration_GMP (fractal, complex_gmp);
fractalresult Barnsleyj1_Iteration_GMP (fractal, complex_gmp);
fractalresult Barnsleym1_Iteration_GMP (fractal, complex_gmp);
fractalresult BurningShip_Iteration_GMP (fractal, complex_gmp);
fractalresult Tricorn_Iteration_GMP (fractal, complex_gmp);
fractalresult Mandelbulb_Iteration_GMP (fractal, complex_gmp);
#endif

 // Fractal Definition
 void Mendelbrot_def (fractal* f);
 void Julia_def (fractal* f);
 void JuliaSin_def (fractal* f);
 void Newton_def (fractal* f);
 void Phoenix_def (fractal* f);
 void Sierpinski_def (fractal* f);
 void Barnsley1j_def (fractal* f);
void Barnsley1m_def (fractal* f);
void Magnet1j_def (fractal* f);
void Magnet1m_def (fractal* f);
void BurningShip_def (fractal* f);
void Tricorn_def (fractal* f);
void Mandelbulb_def (fractal* f);
void Buddhabrot_def (fractal* f);

// Buddhabrot special draw function (density algorithm)
// gui parameter can be NULL if no GUI progress display needed
Uint32 Buddhabrot_Draw (SDL_Surface*, fractal*, int, int, void* gui);

// Utilitaire pour obtenir le nom de la fractale selon son type
const char* Fractal_GetTypeName(int type);

#endif /* ESCAPETIME_H */

