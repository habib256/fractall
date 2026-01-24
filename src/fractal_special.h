/* fractal_special.h
   Declarations des algorithmes speciaux (non escape-time standard)
   released under GPL2
   Copyleft 2001-2026 VERHILLE Arnaud
*/

#ifndef FRACTAL_SPECIAL_H
#define FRACTAL_SPECIAL_H

#include "fractal_types.h"
#include "SDL.h"

// Buddhabrot - algorithme de densite par trajectoires (Type 16)
// gui parameter can be NULL if no GUI progress display needed
Uint32 Buddhabrot_Draw(SDL_Surface* canvas, fractal* f, int decalageX, int decalageY, void* gui);

// Lyapunov - coloration basee sur l'exposant de Lyapunov (Type 17)
Uint32 Lyapunov_Draw(SDL_Surface* canvas, fractal* f, int decalageX, int decalageY, void* gui);

// Nebulabrot - densite RGB avec 3 limites d'iterations differentes (Type 24)
// Red: 50 iter, Green: 500 iter, Blue: 5000 iter
Uint32 Nebulabrot_Draw(SDL_Surface* canvas, fractal* f, int decalageX, int decalageY, void* gui);

#endif /* FRACTAL_SPECIAL_H */
