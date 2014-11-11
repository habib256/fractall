/* SDLGUI.c
 * SDL General User Interface
 * VERHILLE Arnaud 2000-2003
 * Released under GPL2
 */

#include <stdlib.h>
#include "SDLGUI.h"
#include "VonKoch.h"
#include "EscapeTime.h"


gui SDLGUI_Init(int x, int y, int w, int h,int stateH, Uint32 bgcolor, int buttonNumber) {
  gui g;

  g.x= x;
  g.y= y;
  g.w= w;
  g.h= h;

  g.stateH = stateH;

  g.bgcolor = bgcolor;

  g.buttonNumber = buttonNumber;
  g.buttonTab = (int *) malloc ((buttonNumber)* sizeof (button));

  g.barh = 8;  // Hauteur fixe
  g.buttonSize = g.h - (g.barh + 4) - 4;

  g.xplace = 0;  // For now, start by show the first button
  g.xplaceMax = (2 + (g.buttonSize + 2) * g.buttonNumber) - g.w;

  if (g.w >= ((g.buttonSize+2)*buttonNumber)) {
    g.barw = g.w - 4;
  } else {
    g.barw = g.w*(g.w-4)/((g.buttonSize+2)*buttonNumber);
  }

  return g;
}

void SDLGUI_Destroy(gui* g) {
  free (g->buttonTab);
}

int SDLGUI_ButtonNumberRead (gui* g, int x) {
  int buttonNb;
  x = x - 2;
  buttonNb = (int) (x / (g->buttonSize +2));
  return buttonNb;
}

button SDLGUI_Button_Init(char* name, char* pictureLocation) {
}

void SDLGUI_Draw(SDL_Surface* screen, gui* g) {
  int button, x;
  
  // Draw Container
  SDLGUI_Draw3DBox (screen, g->x, g->y, g->w, (g->buttonSize+4), g->bgcolor, 1);

  // Draw bar
  SDLGUI_MenuBar_Draw (screen, g);

  // Draw Widgets
  button = 1;
  x = 2;
  // Le premier bouton que l'on veut dessiner rentre t'il dans l'interface ?
  // D'abord on dessine les boutons qui rentrent dans l'interface
  while ((button <= g->buttonNumber) && (x <= (g->w-(g->buttonSize + 2)))) {
    SDLGUI_Button_Draw (screen, g, x);
    button++;
    x = x + (g->buttonSize + 2);
  }
  SDLGUI_StateBar_Draw (screen, g);

  SDL_UpdateRect (screen, 0, 0, screen->w, screen->h);

} 

void SDLGUI_StateBar_Draw (SDL_Surface* screen, gui* g) {
	SDLGUI_Draw3DBox (screen, 0, g->stateH, g->w, (g->h-g->stateH) , g->bgcolor, 1);

}

void SDLGUI_MenuBar_Draw (SDL_Surface* screen, gui* g) {
  SDLGUI_Draw3DBox (screen, g->x, (g->buttonSize+4), g->w, g->h - (g->buttonSize+4), g->bgcolor, 1);
  SDLGUI_Draw3DBox (screen, 2+(g->xplace*(g->w-(g->barw+4))/g->xplaceMax), g->buttonSize+6, g->barw, g->barh, g->bgcolor, 0);
}

void SDLGUI_Button_Draw (SDL_Surface* screen, gui* g, int xplace) {
  SDLGUI_Draw3DBox (screen, xplace, 2, g->buttonSize, g->buttonSize, g->bgcolor, 0);
  switch (SDLGUI_ButtonNumberRead (g, xplace)+1) {
  case 1:
  	VonKochDraw (screen, xplace+2 ,4,g->buttonSize+xplace-2,g->buttonSize-1,2 );
  break;
  case 2:
  	DragonFractDraw (screen, xplace+2 ,4,g->buttonSize+xplace-2,g->buttonSize-1,8 );
  break;
  default:
  {
  	fractal f;
	f = Fractal_Init (g->buttonSize-4, g->buttonSize-4,SDLGUI_ButtonNumberRead (g, xplace)+1);
	Fractal_Draw (screen, f,xplace+2, 4);
	Fractal_Destroy (f);
	}
	break;

  }
}


void SDLGUI_Draw3DBox (SDL_Surface *surface, int x, int y, int w, int h, Uint32 bgcolor, int type) {
	Uint32 color, colorBlack, colorWhite;
	int epaisseur = 1;
	SDL_Rect dstrect;

	colorBlack = SDL_MapRGB (surface->format, 0, 0, 0);
	colorWhite = SDL_MapRGB (surface->format, 0xFF, 0xFF, 0xFF);

	dstrect.x = x+epaisseur;
	dstrect.y = y+epaisseur;
	dstrect.w = w - 2*(epaisseur);
	dstrect.h = h - 2*(epaisseur);
	SDL_FillRect(surface, &dstrect, bgcolor);

	// franges Lumineuses
	if (type == 0) {
		color = colorWhite;
  } else { color = colorBlack; }

  dstrect.x = x;
  dstrect.y = y;
  dstrect.w = epaisseur;
  dstrect.h = h;
  SDL_FillRect(surface, &dstrect, color);

  dstrect.x = x;
  dstrect.y = y;
  dstrect.w = w;
  dstrect.h = epaisseur;
  SDL_FillRect(surface, &dstrect, color);

  // franges Sombres
  if (type == 0) {
	  color = colorBlack;
  } else { color = colorWhite; }

  dstrect.x = x+epaisseur;
  dstrect.y = y+(h-epaisseur);
  dstrect.w = w-epaisseur;
  dstrect.h = epaisseur;
  SDL_FillRect(surface, &dstrect, color);
  
  dstrect.x = x+(w-epaisseur);
  dstrect.y = y+epaisseur;
  dstrect.w = epaisseur;
  dstrect.h = h-epaisseur;
  SDL_FillRect(surface, &dstrect, color);
}



