/* SDLGUI.c
 * SDL General User Interface
 * VERHILLE Arnaud 2000-2003
 * Released under GPL2
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "SDLGUI.h"
#include "VonKoch.h"
#include "EscapeTime.h"
#include "SDL_gfxPrimitives.h"


gui SDLGUI_Init(int x, int y, int w, int h,int stateH, Uint32 bgcolor, int buttonNumber) {
  gui g;

  g.x= x;
  g.y= y;
  g.w= w;
  g.h= h;

  g.stateH = stateH;

  g.bgcolor = bgcolor;

  g.buttonNumber = buttonNumber;

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
  // Rien à libérer pour l'instant
}

int SDLGUI_ButtonNumberRead (gui* g, int x) {
  int buttonNb;
  x = x - 2;
  buttonNb = (int) (x / (g->buttonSize +2));
  return buttonNb;
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
  case 16:
  {
  	// Buddhabrot utilise son propre algorithme de rendu
  	fractal f;
	f = Fractal_Init (g->buttonSize-4, g->buttonSize-4, 16);
	Buddhabrot_Draw (screen, &f, xplace+2, 4, NULL);  // NULL = pas de progression pour le bouton
	Fractal_Destroy (f);
  }
  break;
  case 17:
  {
  	// Le type Lyapunov (17) utilise son propre algorithme de rendu
  	fractal f;
	f = Fractal_Init (g->buttonSize-4, g->buttonSize-4, SDLGUI_ButtonNumberRead (g, xplace)+1);
	Fractal_Draw (screen, f, xplace+2, 4, NULL);  // NULL = pas de progression pour le bouton
	Fractal_Destroy (f);
  }
  break;
  default:
  {
  	fractal f;
	f = Fractal_Init (g->buttonSize-4, g->buttonSize-4,SDLGUI_ButtonNumberRead (g, xplace)+1);
	Fractal_Draw (screen, f, xplace+2, 4, NULL);  // NULL = pas de progression pour le bouton
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

void SDLGUI_StateBar_Update (SDL_Surface* screen, gui* g, int type, int colorMode, double centerX, double centerY, double zoomFactor, Uint32 renderTime, fractal* f) {
	const char* typeNames[] = {"", "Von Koch", "Dragon", "Mandelbrot", "Julia", "Julia Sin",
		"Newton", "Phoenix", "", "Barnsley J", "Barnsley M",
		"Magnet J", "Magnet M", "Burning Ship", "Tricorn", "Mandelbulb", "Buddhabrot"};
	const char* paletteNames[] = {"SmoothFire", "SmoothOcean", "SmoothForest", "SmoothViolet", "SmoothRainbow"};
	char statusText[256];
	char precisionText[64];
	int precisionX;

	// Redessiner la barre d'état
	SDLGUI_StateBar_Draw(screen, g);

	// Formater le texte principal (gauche)
	if (type >= 1 && type <= 16) {
		// Afficher le zoom en notation adaptée selon la valeur
		if (zoomFactor >= 1e6) {
			snprintf(statusText, sizeof(statusText), "%s | %s | x%.2e | (%.6f, %.6f) | %dms",
				typeNames[type], paletteNames[colorMode], zoomFactor, centerX, centerY, renderTime);
		} else if (zoomFactor >= 1000) {
			snprintf(statusText, sizeof(statusText), "%s | %s | x%.0f | (%.6f, %.6f) | %dms",
				typeNames[type], paletteNames[colorMode], zoomFactor, centerX, centerY, renderTime);
		} else {
			snprintf(statusText, sizeof(statusText), "%s | %s | x%.1f | (%.3f, %.3f) | %dms",
				typeNames[type], paletteNames[colorMode], zoomFactor, centerX, centerY, renderTime);
		}
	} else {
		snprintf(statusText, sizeof(statusText), "Type %d | %s", type, paletteNames[colorMode]);
	}

	// Afficher le texte principal à gauche
	stringRGBA(screen, 5, g->stateH + 5, statusText, 0, 0, 0, 255);

	// Formater et afficher la précision en bas à droite
	if (f != NULL) {
#ifdef HAVE_GMP
		if (f->use_gmp) {
			snprintf(precisionText, sizeof(precisionText), "GMP %lu bits", (unsigned long)f->gmp_precision);
		} else {
			snprintf(precisionText, sizeof(precisionText), "double");
		}
#else
		snprintf(precisionText, sizeof(precisionText), "double");
#endif
		// Calculer la position X pour aligner à droite (approximatif, chaque caractère fait ~6 pixels)
		// Un peu moins à droite : on ajoute 20 pixels de marge
		precisionX = g->w - strlen(precisionText) * 6 - 25;
		if (precisionX < 5) precisionX = 5; // Éviter le débordement à gauche
		stringRGBA(screen, precisionX, g->stateH + 5, precisionText, 0, 0, 0, 255);
	}

	SDL_UpdateRect(screen, 0, g->stateH, g->w, screen->h - g->stateH);
}

void SDLGUI_StateBar_Progress (SDL_Surface* screen, gui* g, int percent, const char* task) {
	char progressText[128];
	int barWidth, barHeight, barX, barY;
	int fillWidth;
	SDL_Rect bgRect, fillRect;
	Uint32 bgColor, fillColor, borderColor;

	// Redessiner la barre d'état
	SDLGUI_StateBar_Draw(screen, g);

	// Couleurs
	bgColor = SDL_MapRGB(screen->format, 200, 200, 200);
	fillColor = SDL_MapRGB(screen->format, 80, 180, 80);
	borderColor = SDL_MapRGB(screen->format, 100, 100, 100);

	// Dimensions de la barre de progression
	barHeight = 12;
	barWidth = 200;
	barX = g->w - barWidth - 10;
	barY = g->stateH + 4;

	// Fond de la barre
	bgRect.x = barX;
	bgRect.y = barY;
	bgRect.w = barWidth;
	bgRect.h = barHeight;
	SDL_FillRect(screen, &bgRect, bgColor);

	// Bordure
	rectangleRGBA(screen, barX - 1, barY - 1, barX + barWidth, barY + barHeight, 100, 100, 100, 255);

	// Remplissage selon le pourcentage
	if (percent > 0) {
		fillWidth = (barWidth * percent) / 100;
		if (fillWidth > barWidth) fillWidth = barWidth;
		fillRect.x = barX;
		fillRect.y = barY;
		fillRect.w = fillWidth;
		fillRect.h = barHeight;
		SDL_FillRect(screen, &fillRect, fillColor);
	}

	// Texte : tâche et pourcentage
	snprintf(progressText, sizeof(progressText), "%s %d%%", task, percent);
	stringRGBA(screen, 5, g->stateH + 5, progressText, 0, 0, 0, 255);

	SDL_UpdateRect(screen, 0, g->stateH, g->w, screen->h - g->stateH);
}

