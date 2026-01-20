/* png_save.h
 * Sauvegarde PNG avec métadonnées pour fractall
 * Released under GPL2
 */

#ifndef PNG_SAVE_H
#define PNG_SAVE_H

#include "EscapeTime.h"
#include "SDL.h"

// Sauvegarder une surface SDL en PNG avec métadonnées de fractale
// Retourne 0 en cas de succès, -1 en cas d'erreur
int SavePNGWithMetadata(SDL_Surface* surface, const char* filename, fractal* f, int typeFractale);

// Charger les métadonnées depuis un PNG (pour recalculer la fractale)
// Retourne 0 en cas de succès, -1 si les métadonnées sont absentes ou invalides
int LoadPNGMetadata(const char* filename, fractal* f, int* typeFractale);

#endif /* PNG_SAVE_H */
