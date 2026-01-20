/* SDLGUI.h
 * SDL General User Interface
 * VERHILLE Arnaud 2000
 * Released under GPL2
 */

#include "SDL.h"

// Include EscapeTime.h pour avoir accès à la structure fractal
// Utilisation d'une garde pour éviter les inclusions circulaires
#ifndef ESCAPETIME_H
#include "EscapeTime.h"
#endif

typedef struct {	// Definition d'un bouton
  char* name;
  char* pictureName;
} button;

typedef struct {  // Un menu utilisable ??
int x,y,w,h;
} menu;

typedef struct gui_struct {   // barre de boutton
  int x,y,w,h;      // Largeur et hauteur barre de boutton
  int stateH;	    // Bad hack to add state bar
  Uint32 bgcolor;    // La couleur de fond
  int buttonNumber;  // Le nombre de boutton dans ce menu
  int buttonSize;
  int xplace;        // xplace correspond a x=2 l'ecran
  int xplaceMax;     // Valeur Maximale de xplace
  int barh,barw;     // Largeur de la barre de deplacement
  int selectedType;  // Type de fractale actuellement sélectionné (1-23)
  int hoverButton;   // Bouton survolé (-1 si aucun)
} gui;

extern gui SDLGUI_Init (int, int, int, int, int, Uint32, int);
extern void SDLGUI_Destroy (gui*);
extern int SDLGUI_ButtonNumberRead (gui*, int);
extern void SDLGUI_Draw (SDL_Surface*, gui*);

void SDLGUI_MenuBar_Draw (SDL_Surface*, gui*);

extern void SDLGUI_Button_Draw (SDL_Surface*, gui*, int);
void SDLGUI_Draw3DBox (SDL_Surface *surface, int x, int y, int w, int h, Uint32 bgcolor, int type);
void SDLGUI_StateBar_Draw (SDL_Surface* screen, gui* g);
void SDLGUI_StateBar_Update (SDL_Surface* screen, gui* g, int type, int colorMode, double centerX, double centerY, double zoomFactor, Uint32 renderTime, fractal* f);
void SDLGUI_StateBar_Progress (SDL_Surface* screen, gui* g, int percent, const char* task);
