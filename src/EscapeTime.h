/* EscapeTime.h
   Fractal formula header - API publique unifiee
   released under GPL2
   Copyleft 2001-2026 VERHILLE Arnaud
*/

#ifndef ESCAPETIME_H
#define ESCAPETIME_H

#include "config.h"
#include "SDL.h"

// Inclure les sous-headers modulaires
#include "fractal_types.h"
#include "fractal_iterations.h"
#include "fractal_definitions.h"
#include "fractal_special.h"

#ifdef HAVE_GMP
#include "fractal_iterations_gmp.h"
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

// ******************
// Interface publique
// ******************

// Initialisation et destruction
fractal Fractal_Init(int screenW, int screenH, int type);
void Fractal_Destroy(fractal f);

// Calcul de la matrice de complexe z lors de leur derniere iteration
void Fractal_CalculateMatrix(fractal* f);
void Fractal_CalculateMatrix_DDp1(fractal* f, SDL_Surface* canvas, void* gui, int* progress, int progressStart, int progressEnd, const char* fractalName, int decalageX, int decalageY);
void Fractal_CalculateMatrix_DDp2(fractal* f, SDL_Surface* canvas, void* gui, int* progress, int progressStart, int progressEnd, const char* fractalName, int decalageX, int decalageY);

// Calcul de la couleur
void Fractal_CalculateColorMatrix(fractal* f, SDL_Surface* canvas, void* gui, int* progress, int progressStart, int progressEnd);
void FractalColorTest(fractal* f);

// Color Formulae Utilities
color Fractal_ReadColorMatrix(fractal f, int x, int y);
int Fractal_ReadColorMatrixRed(fractal f, int x, int y);
int Fractal_ReadColorMatrixGreen(fractal f, int x, int y);
int Fractal_ReadColorMatrixBlue(fractal f, int x, int y);

// Rendu principal
Uint32 Fractal_Draw(SDL_Surface* canvas, fractal f, int decalageX, int decalageY, void* gui);
void Fractal_ChangeType(fractal* f, int type);

// Utilitaire pour obtenir le nom de la fractale selon son type
const char* Fractal_GetTypeName(int type);

#endif /* ESCAPETIME_H */
