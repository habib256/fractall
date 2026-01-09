/*
*     FractAll, a portable fractal viewer
*                  Copyright (C) 2002-2003 by
*
*      Arnaud VERHILLE      (gist@wanadoo.fr)
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EscapeTime.h"
#include "SDLGUI.h"
#include "SDL_gfxPrimitives.h"
#include "VonKoch.h"
#ifdef HAVE_GMP
#include <gmp.h>
#include "precision_detector.h"
#endif



#define PACKAGE "fractall"
#ifdef VERSION
#undef VERSION
#endif
#define VERSION 0.5




typedef struct { // Conteneur graphique en 3 parties
	int y1, y2;
	int w,h;
} window;


// FUNCTIONS PROTOTYPES
//*********************

int EventFlush (SDL_Event*);
int EventCheck (SDL_Event*, SDL_Surface *, gui *, fractal *, int*, int*, window win, Uint32* renderTime);


// MAIN FUNCTION
//**************
int
main (int argc, char *argv[])
{
	
	// Initialisation des variables
	// ****************************
	
	SDL_Surface *screen;		// Surface principale de la fenetre
	
	int fullscreen = 0;		// Boolean for fullscreen or not
	int useGui = 1;			// Boolean for gui or not
	int screenW = 800, screenH = 600;	// Largeur et hauteur de la fenetre
	int menuH = 51;			// C'est la place que prend le gui
	int stateH = screenH - 20;		// Place Barre d etat
	int typeFractale = 1;		//Type de fractale par defaut
	int iteration = 1;   		// Nbr d iteration de la fractale par defaut
	
	SDL_Event event;		// Gestion des evenements
	
	fractal f;			// Allocate a new fractal.
	gui g;			// Allocate a new gui
	window win;		// Creation d un  nouveau conteneur
	Uint32 renderTime = 0;	// Temps de rendu en ms

	int i;			//tmp
	
	// Gestion de la ligne de commande
	// *******************************
	
	
	for (i = 1; (i < argc) && (*argv[i] == '-'); i++)
		switch (argv[i][1])
	{
      case 'h':
		  fprintf (stdout, "%s v %lf (GPL2) VERHILLE A.\n", PACKAGE, VERSION);
		  fprintf (stdout, "Usage : fractall [OPTION]...\n");
		  fprintf (stdout, "   -help Print this text\n");
		  fprintf (stdout,
			  "   -x Pixel width for main window (default : 640)\n");
		  fprintf (stdout,
			  "   -y Pixel height for main window (default : 480)\n");
		  fprintf (stdout, "   -g Gui height (default : 66 | min : 20)\n");
		  fprintf (stdout,
			  "   -f Set Fullscreen mode (Press ESC or Q to quit)\n");
		  fprintf (stdout,
			  "   -nogui No user interface (key selection still work)\n");
		  fprintf (stdout, "Example : fractall -f -x640 -y480\n");
		  // Shutdown All System
		  SDL_Quit ();
		  exit (0);
		  break;
      case 'x':
		  sscanf (argv[i], "-x%i", &screenW);
		  break;
      case 'y':
		  sscanf (argv[i], "-y%i", &screenH);
		  stateH = screenH - 20;
		  break;
      case 'f':
		  fullscreen = 1;
		  break;
      case 'g':
		  sscanf (argv[i], "-g%i", &menuH);
		  break;
      case 'n':
		  useGui = 0;
		  menuH = 0;
		  stateH = screenH;
      default:
		  printf ("Unknown parameter : -%c", argv[i][1]);
		  break;
	}
	
	// Initialisation de la Librairie SDL
	// **********************************
	if ((SDL_Init (SDL_INIT_VIDEO) == -1))
    {
		fprintf (stderr, "SDL Initialization error: %s.\n", SDL_GetError ());
		exit (-1);
    }
	atexit (SDL_Quit);
	/* Set the Window Manager Config */
	SDL_WM_SetCaption (PACKAGE, PACKAGE);
	/* Initialize the display on screenWxscreenH bpp-bits mode
	requesting a software surface         */
	if (fullscreen)
    {
		screen =
			SDL_SetVideoMode (screenW, screenH, 16,
			SDL_SWSURFACE | SDL_FULLSCREEN);
    }
	else
    {
		screen = SDL_SetVideoMode (screenW, screenH, 16, SDL_SWSURFACE);
    }
	if (screen == NULL)
    {
		fprintf (stderr, "Couldn't set %dx%dx16 Graphic Mode : %s\n", screenW,
			screenH, SDL_GetError ());
		exit (1);
    }
	printf ("SDL Initialized. \n");
	
	// On initialise le contenu de win, le menu et la fractale
	// **********************************************
	
	// Initialisation win
	win.y1 = menuH;
	win.y2 = stateH;
	win.w = screen->w;
	win.h = screen->h;
	
	// On initialise le menu si useGui = TRUE
	if (useGui)
		g = SDLGUI_Init (0, 0, win.w, win.y1, stateH, SDL_MapRGB (screen->format, 213, 214, 213), 16);	// 16 boutons (types 1-16)
	if (useGui)
		SDLGUI_Draw (screen, &g);
	
	// On initialise la fractale
	f = Fractal_Init (win.w, win.y2-win.y1,typeFractale);
	
	// Tracer les diverses fractales
	// *****************************
	
	switch (typeFractale)
    {
		// Trace de la mendelbrot
    case 1:
		VonKochDraw (screen,0,win.y1,win.w,win.y2, iteration);	// Trace de Von koch
		break;
		
    case 2:
		DragonFractDraw (screen,0,win.y1,win.w,win.y2, iteration);	// Trace du Dragon
		break;
    default:
    		Fractal_ChangeType (&f, typeFractale);
		Fractal_Draw (screen, f, 0, win.y1, (useGui ? (void*)&g : NULL));
		break;
    }
    
    SDL_UpdateRect (screen, 0, 0, screen->w, screen->h);

	//MAIN INFINITE BOUCLE
	for (;;)
    {
		EventFlush (&event);	// Event when Initializing
		while (SDL_WaitEvent (&event) >= 0)
		{
			EventCheck (&event, screen, &g, &f, &typeFractale, &iteration, win, &renderTime);
		}
    }
	return 1;			// Error, should never reach this line !
}

// Gestion des evenements
// **********************
// **********************

int
EventFlush (SDL_Event* event)
{
	while (SDL_PollEvent (event))
    {
    }
	//printf ("Event queue empty\n");
	return 0;
}

int EventCheck (SDL_Event* event, SDL_Surface* screen, gui* g, fractal* f,
				int* typeFractale, int* iteration, window win, Uint32* renderTime)
{
	double centerX, centerY;
	switch (event->type)
    {
		
		// Gestion des touches du clavier
		// ******************************
		
    case SDL_KEYDOWN:
		{
			if ((event->key.keysym.sym == SDLK_ESCAPE) || (event->key.keysym.sym ==
				SDLK_q))
			{
				// Shutdown All System
				printf ("Quit requested, quitting.\n");
				Fractal_Destroy (*f);
				if (g != NULL)
					SDLGUI_Destroy (g);			
				SDL_FreeSurface (screen);
				exit (0);
			}
			if (event->key.keysym.sym == SDLK_s)
			{			// Screenshot
				if (SDL_SaveBMP (screen, "Screenshot.bmp") == 1)
				{
					printf ("Cannot save the screenshot in BMP file\n");
				}
				else
				{
					printf ("Screenshot.bmp file saved successfully\n");
				}
			}
			if (event->key.keysym.sym == SDLK_c)
			{			// Changer palette de couleur
				const char* palettes[] = {"SmoothFire", "Rainbow", "SmoothOcean"};
				f->colorMode = (f->colorMode + 1) % 3;
				printf ("Palette: %s\n", palettes[f->colorMode]);
				if (*typeFractale >= 3) {
					if (*typeFractale == 16) {
						*renderTime = Buddhabrot_Draw (screen, f, 0, win.y1, g);
					} else {
						// Recalculer seulement la colorisation (sans progression car rapide)
						int dummyProgress = 0;
						Fractal_CalculateColorMatrix(f, screen, NULL, &dummyProgress, 0, 100);
						// Redessiner la fractale avec les nouvelles couleurs
						int i, j;
						Uint8 r, g_color, b;
						for (i=0; i<f->xpixel; i++) {
							for (j=0; j<f->ypixel; j++) {
								r=Fractal_ReadColorMatrixRed (*f,i,j);
								g_color=Fractal_ReadColorMatrixGreen (*f,i,j);
								b=Fractal_ReadColorMatrixBlue (*f,i,j);
								pixelRGBA(screen, (Sint16) (i), (Sint16) (j+win.y1),  r, g_color,  b, 255);
							}
						}
						SDL_UpdateRect (screen, 0, 0, screen->w, screen->h);
						*renderTime = 0; // Pas de nouveau calcul, juste changement de palette
					}
					centerX = (f->xmin + f->xmax) / 2;
					centerY = (f->ymin + f->ymax) / 2;
					if (g != NULL)
						SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
				}
			}

			if (event->key.keysym.sym == SDLK_F1)
			{
				*typeFractale =1;
				*iteration = 1;
				VonKochDraw (screen,0,win.y1,win.w,win.y2, *iteration);	// Trace de Von koch
			}

			if (event->key.keysym.sym == SDLK_F2)
			{
				*typeFractale =2;
				*iteration = 4;
				DragonFractDraw (screen,0,win.y1,win.w,win.y2, *iteration);	// Trace du Dragon
			}

			if (event->key.keysym.sym == SDLK_F3)
			{
				*typeFractale = 3;
				Fractal_ChangeType (f, 3);
				*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}
			if (event->key.keysym.sym == SDLK_F4)
			{
				*typeFractale = 4;
				Fractal_ChangeType (f, 4);
				*renderTime = Fractal_Draw (screen, *f,0, win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}
			if (event->key.keysym.sym == SDLK_F5)
			{
				*typeFractale = 5;
				Fractal_ChangeType (f, 5);
				*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}
			if (event->key.keysym.sym == SDLK_F6)
			{
				*typeFractale = 6;
				Fractal_ChangeType (f, 6);
				*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}
			if (event->key.keysym.sym == SDLK_F7)
			{
				*typeFractale = 7;
				Fractal_ChangeType (f, 7);
				*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}
			if (event->key.keysym.sym == SDLK_F8)
			{
				*typeFractale = 8;
				Fractal_ChangeType (f, 8);
				*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}
			if (event->key.keysym.sym == SDLK_F9)
			{
				*typeFractale = 13;
				Fractal_ChangeType (f, 13);
				*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}
			if (event->key.keysym.sym == SDLK_F10)
			{
				*typeFractale = 14;
				Fractal_ChangeType (f, 14);
				*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}
			if (event->key.keysym.sym == SDLK_F11)
			{
				*typeFractale = 15;
				Fractal_ChangeType (f, 15);
				*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}
			if (event->key.keysym.sym == SDLK_F12)
			{
				*typeFractale = 16;
				Fractal_ChangeType (f, 16);
				*renderTime = Buddhabrot_Draw (screen, f, 0, win.y1, g);
				centerX = (f->xmin + f->xmax) / 2; centerY = (f->ymin + f->ymax) / 2;
				if (g != NULL) SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
			}

		}
		break;
		
		// Gestion des evenements du window Manager
		// *****************************************
		
    case SDL_QUIT:
		{
			// Shutdown All System
			printf ("Quit requested, quitting.\n");
			Fractal_Destroy (*f);
			if (g != NULL)
				SDLGUI_Destroy (g);
			SDL_FreeSurface (screen);
			SDL_Quit ();
			exit (0);
		}
		break;
		
		// Gestion des evenements souris
		// ******************************
		
		
    case SDL_MOUSEBUTTONDOWN:	// On clique sur un bouton !!
		{
			if (g != NULL) {
				if (event->button.y <= g->h)	// C'est du boulot pour le GUI ?
				{
					int buttonNumber = SDLGUI_ButtonNumberRead (g, event->button.x);
					printf ("GUI Event n %d Detected ...\n", buttonNumber);
					*typeFractale = buttonNumber+1;
					switch (*typeFractale) {
					case 1:
					*typeFractale =1;
					*iteration = 1;
					VonKochDraw (screen,0,win.y1,win.w,win.y2, *iteration);
					break;
					case 2:
					*typeFractale =2;
					*iteration = 4;
					DragonFractDraw (screen,0,win.y1,win.w,win.y2, *iteration);
					break;
					case 16:
					{
					Fractal_ChangeType (f, 16);
					*renderTime = Buddhabrot_Draw (screen, f, 0, win.y1, g);
					centerX = (f->xmin + f->xmax) / 2;
					centerY = (f->ymin + f->ymax) / 2;
					SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
					}
					break;
					default:
					{
					Fractal_ChangeType (f, *typeFractale);
					*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
					centerX = (f->xmin + f->xmax) / 2;
					centerY = (f->ymin + f->ymax) / 2;
					SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
					}
					}
					return 0;
				}
			}
			if (event->button.y >= win.y2)	// On clique dans la barre d etat state !
			{
				printf ("Click sur la barre d etat ...\n");
				return 0;
			}
			
			
			//printf("typeFractale = %d \n",*typeFractale);
			switch (*typeFractale)	// Comportement pour quel type de fractale ?
			{
			case 1:
				{
					int i;
					i = *iteration;
					
					if (!(*iteration == 0) && (event->button.button == SDL_BUTTON_RIGHT))
						*iteration = *iteration -1;
					if (!(*iteration == 8) && (event->button.button == SDL_BUTTON_LEFT))
						*iteration = *iteration + 1;
					
					if (i != *iteration) {
						VonKochDraw (screen,0,win.y1,win.w,win.y2, *iteration);	// Trace de Von koch
						printf ("Draw VonKoch with %d iteration\n",*iteration);
					}
					
				}
				break;
				
			case 2:
				{
					int i;
					i = *iteration;
					
					if (!(*iteration == 0) && (event->button.button == SDL_BUTTON_RIGHT))
						*iteration = *iteration -1;
					if (!(*iteration == 20) && (event->button.button == SDL_BUTTON_LEFT))
						*iteration = *iteration + 1;
					
					if (i != *iteration) {
						Uint32 time; // Test de temps de calcul en ms
						
						time = SDL_GetTicks();
						DragonFractDraw (screen,0,win.y1,win.w,win.y2, *iteration);	// Trace du Dragon
						time = SDL_GetTicks() - time;
						printf ("Draw Dragon with %d iteration in %d ms\n",*iteration, time);
						
					}
					
				}
			break;
			default:
						{			// EscapeTime Fractal Zoom routine
					
					if (event->button.button == SDL_BUTTON_LEFT)
					{		// ZOOM
#ifdef HAVE_GMP
						if (f->use_gmp) {
							// Calculer le zoom directement en GMP pour préserver la précision
							mpf_t newCenterX_gmp, newCenterY_gmp, newxmax_gmp, newxmin_gmp, newymax_gmp, newymin_gmp;
							mpf_t range_x_gmp, range_y_gmp, half_range_x_gmp, half_range_y_gmp;
							mpf_t zoomfactor_gmp, two_gmp;
							mp_bitcnt_t prec = f->gmp_precision;
							
							mpf_init2(newCenterX_gmp, prec);
							mpf_init2(newCenterY_gmp, prec);
							mpf_init2(newxmax_gmp, prec);
							mpf_init2(newxmin_gmp, prec);
							mpf_init2(newymax_gmp, prec);
							mpf_init2(newymin_gmp, prec);
							mpf_init2(range_x_gmp, prec);
							mpf_init2(range_y_gmp, prec);
							mpf_init2(half_range_x_gmp, prec);
							mpf_init2(half_range_y_gmp, prec);
							mpf_init2(zoomfactor_gmp, prec);
							mpf_init2(two_gmp, prec);
							
							mpf_set_ui(zoomfactor_gmp, f->zoomfactor);
							mpf_set_ui(two_gmp, 2);
							
							// Calculer le centre du zoom à partir de la position du clic
							if (event->button.x >= 0 && event->button.x < f->xpixel &&
							    event->button.y >= win.y1 && event->button.y < win.y2) {
								mpf_t pixel_x_gmp, pixel_y_gmp, step_x_gmp, step_y_gmp;
								mpf_t range_x_temp, range_y_temp;
								mpf_init2(pixel_x_gmp, prec);
								mpf_init2(pixel_y_gmp, prec);
								mpf_init2(step_x_gmp, prec);
								mpf_init2(step_y_gmp, prec);
								mpf_init2(range_x_temp, prec);
								mpf_init2(range_y_temp, prec);
								
								mpf_set_ui(pixel_x_gmp, event->button.x);
								mpf_set_ui(pixel_y_gmp, event->button.y - win.y1);
								mpf_set_ui(step_x_gmp, f->xpixel);
								mpf_set_ui(step_y_gmp, f->ypixel);
								
								mpf_sub(range_x_temp, f->xmax_gmp, f->xmin_gmp);
								mpf_sub(range_y_temp, f->ymax_gmp, f->ymin_gmp);
								
								mpf_div(step_x_gmp, range_x_temp, step_x_gmp);
								mpf_mul(newCenterX_gmp, step_x_gmp, pixel_x_gmp);
								mpf_add(newCenterX_gmp, newCenterX_gmp, f->xmin_gmp);
								
								mpf_div(step_y_gmp, range_y_temp, step_y_gmp);
								mpf_mul(newCenterY_gmp, step_y_gmp, pixel_y_gmp);
								mpf_add(newCenterY_gmp, newCenterY_gmp, f->ymin_gmp);
								
								mpf_clear(pixel_x_gmp);
								mpf_clear(pixel_y_gmp);
								mpf_clear(step_x_gmp);
								mpf_clear(step_y_gmp);
								mpf_clear(range_x_temp);
								mpf_clear(range_y_temp);
							} else {
								// Si le clic est hors zone, zoomer sur le centre actuel
								mpf_add(newCenterX_gmp, f->xmin_gmp, f->xmax_gmp);
								mpf_div(newCenterX_gmp, newCenterX_gmp, two_gmp);
								mpf_add(newCenterY_gmp, f->ymin_gmp, f->ymax_gmp);
								mpf_div(newCenterY_gmp, newCenterY_gmp, two_gmp);
							}
							
							// Calculer la nouvelle plage : diviser par le facteur de zoom
							mpf_sub(range_x_gmp, f->xmax_gmp, f->xmin_gmp);
							mpf_sub(range_y_gmp, f->ymax_gmp, f->ymin_gmp);
							mpf_div(range_x_gmp, range_x_gmp, zoomfactor_gmp);
							mpf_div(range_y_gmp, range_y_gmp, zoomfactor_gmp);
							mpf_div(half_range_x_gmp, range_x_gmp, two_gmp);
							mpf_div(half_range_y_gmp, range_y_gmp, two_gmp);
							
							// Centrer la nouvelle vue sur le point cliqué
							mpf_add(newxmax_gmp, newCenterX_gmp, half_range_x_gmp);
							mpf_sub(newxmin_gmp, newCenterX_gmp, half_range_x_gmp);
							mpf_add(newymax_gmp, newCenterY_gmp, half_range_y_gmp);
							mpf_sub(newymin_gmp, newCenterY_gmp, half_range_y_gmp);
							
							// Mettre à jour les coordonnées GMP
							mpf_set(f->xmax_gmp, newxmax_gmp);
							mpf_set(f->xmin_gmp, newxmin_gmp);
							mpf_set(f->ymax_gmp, newymax_gmp);
							mpf_set(f->ymin_gmp, newymin_gmp);
							
							// Synchroniser vers double pour l'affichage
							f->xmax = mpf_get_d(newxmax_gmp);
							f->xmin = mpf_get_d(newxmin_gmp);
							f->ymax = mpf_get_d(newymax_gmp);
							f->ymin = mpf_get_d(newymin_gmp);
							
							printf ("Zoom x%d  to X= %f\tY= %f\n", f->zoomfactor,
								mpf_get_d(newCenterX_gmp), mpf_get_d(newCenterY_gmp));
							
							mpf_clear(newCenterX_gmp);
							mpf_clear(newCenterY_gmp);
							mpf_clear(newxmax_gmp);
							mpf_clear(newxmin_gmp);
							mpf_clear(newymax_gmp);
							mpf_clear(newymin_gmp);
							mpf_clear(range_x_gmp);
							mpf_clear(range_y_gmp);
							mpf_clear(half_range_x_gmp);
							mpf_clear(half_range_y_gmp);
							mpf_clear(zoomfactor_gmp);
							mpf_clear(two_gmp);
							
							precision_update_fractal(f);
							// Mettre à jour la précision des coordonnées GMP si nécessaire
							if (f->gmp_precision != prec) {
								mpf_set_prec(f->xmin_gmp, f->gmp_precision);
								mpf_set_prec(f->xmax_gmp, f->gmp_precision);
								mpf_set_prec(f->ymin_gmp, f->gmp_precision);
								mpf_set_prec(f->ymax_gmp, f->gmp_precision);
							}
							// Réallouer zmatrix_gmp si nécessaire
							if (f->use_gmp && f->zmatrix_gmp == NULL) {
								f->zmatrix_gmp = (complex_gmp *) malloc ((f->xpixel*f->ypixel)* sizeof (complex_gmp));
								if (f->zmatrix_gmp == NULL) {
									fprintf(stderr, "Erreur allocation mémoire GMP fractale\n");
								} else {
									for (int k = 0; k < f->xpixel * f->ypixel; k++) {
										complex_gmp_init(&f->zmatrix_gmp[k], f->gmp_precision);
									}
								}
							} else if (!f->use_gmp && f->zmatrix_gmp != NULL) {
								for (int k = 0; k < f->xpixel * f->ypixel; k++) {
									complex_gmp_clear(&f->zmatrix_gmp[k]);
								}
								free(f->zmatrix_gmp);
								f->zmatrix_gmp = NULL;
							}
						} else {
#endif
						double newCenterX, newCenterY, newxmax, newxmin, newymax,
							newymin;
						double range_x, range_y, half_range_x, half_range_y;
						
						// Calculer le centre du zoom à partir de la position du clic
						// Vérifier que le clic est dans la zone de la fractale
						if (event->button.x >= 0 && event->button.x < f->xpixel &&
						    event->button.y >= win.y1 && event->button.y < win.y2) {
							newCenterX =
								event->button.x * ((f->xmax - f->xmin) / f->xpixel) +
								f->xmin;
							newCenterY =
								(event->button.y - win.y1) * ((f->ymax - f->ymin) / f->ypixel) +
								f->ymin;
						} else {
							// Si le clic est hors zone, zoomer sur le centre actuel
							newCenterX = (f->xmin + f->xmax) / 2.0;
							newCenterY = (f->ymin + f->ymax) / 2.0;
						}
						
						printf ("Zoom x%d  to X= %f\tY= %f\n", f->zoomfactor,
							newCenterX, newCenterY);
						
						// Calculer la nouvelle plage : diviser par le facteur de zoom
						range_x = (f->xmax - f->xmin) / f->zoomfactor;
						range_y = (f->ymax - f->ymin) / f->zoomfactor;
						half_range_x = range_x / 2.0;
						half_range_y = range_y / 2.0;
						
						// Centrer la nouvelle vue sur le point cliqué
						newxmax = newCenterX + half_range_x;
						newxmin = newCenterX - half_range_x;
						newymax = newCenterY + half_range_y;
						newymin = newCenterY - half_range_y;
						
						f->xmax = newxmax;
						f->ymax = newymax;
						f->xmin = newxmin;
						f->ymin = newymin;
#ifdef HAVE_GMP
						precision_update_fractal(f);
						// Synchroniser vers GMP si nécessaire
						if (f->use_gmp) {
							mpf_set_d(f->xmin_gmp, f->xmin);
							mpf_set_d(f->xmax_gmp, f->xmax);
							mpf_set_d(f->ymin_gmp, f->ymin);
							mpf_set_d(f->ymax_gmp, f->ymax);
						}
						// Réallouer zmatrix_gmp si nécessaire
						if (f->use_gmp && f->zmatrix_gmp == NULL) {
							f->zmatrix_gmp = (complex_gmp *) malloc ((f->xpixel*f->ypixel)* sizeof (complex_gmp));
							if (f->zmatrix_gmp == NULL) {
								fprintf(stderr, "Erreur allocation mémoire GMP fractale\n");
							} else {
								for (int k = 0; k < f->xpixel * f->ypixel; k++) {
									complex_gmp_init(&f->zmatrix_gmp[k], f->gmp_precision);
								}
							}
						} else if (!f->use_gmp && f->zmatrix_gmp != NULL) {
							for (int k = 0; k < f->xpixel * f->ypixel; k++) {
								complex_gmp_clear(&f->zmatrix_gmp[k]);
							}
							free(f->zmatrix_gmp);
							f->zmatrix_gmp = NULL;
						}
						}
#endif
					}
					if (event->button.button == SDL_BUTTON_RIGHT)
					{		// UNZOOM
#ifdef HAVE_GMP
						if (f->use_gmp) {
							// Calculer le unzoom directement en GMP pour préserver la précision
							mpf_t newCenterX_gmp, newCenterY_gmp, newxmax_gmp, newxmin_gmp, newymax_gmp, newymin_gmp;
							mpf_t range_x_gmp, range_y_gmp, half_range_x_gmp, half_range_y_gmp;
							mpf_t zoomfactor_gmp, two_gmp;
							mp_bitcnt_t prec = f->gmp_precision;
							
							mpf_init2(newCenterX_gmp, prec);
							mpf_init2(newCenterY_gmp, prec);
							mpf_init2(newxmax_gmp, prec);
							mpf_init2(newxmin_gmp, prec);
							mpf_init2(newymax_gmp, prec);
							mpf_init2(newymin_gmp, prec);
							mpf_init2(range_x_gmp, prec);
							mpf_init2(range_y_gmp, prec);
							mpf_init2(half_range_x_gmp, prec);
							mpf_init2(half_range_y_gmp, prec);
							mpf_init2(zoomfactor_gmp, prec);
							mpf_init2(two_gmp, prec);
							
							mpf_set_ui(zoomfactor_gmp, f->zoomfactor);
							mpf_set_ui(two_gmp, 2);
							
							// Calculer le centre du zoom à partir de la position du clic
							if (event->button.x >= 0 && event->button.x < f->xpixel &&
							    event->button.y >= win.y1 && event->button.y < win.y2) {
								mpf_t pixel_x_gmp, pixel_y_gmp, step_x_gmp, step_y_gmp;
								mpf_t range_x_temp, range_y_temp;
								mpf_init2(pixel_x_gmp, prec);
								mpf_init2(pixel_y_gmp, prec);
								mpf_init2(step_x_gmp, prec);
								mpf_init2(step_y_gmp, prec);
								mpf_init2(range_x_temp, prec);
								mpf_init2(range_y_temp, prec);
								
								mpf_set_ui(pixel_x_gmp, event->button.x);
								mpf_set_ui(pixel_y_gmp, event->button.y - win.y1);
								mpf_set_ui(step_x_gmp, f->xpixel);
								mpf_set_ui(step_y_gmp, f->ypixel);
								
								mpf_sub(range_x_temp, f->xmax_gmp, f->xmin_gmp);
								mpf_sub(range_y_temp, f->ymax_gmp, f->ymin_gmp);
								
								mpf_div(step_x_gmp, range_x_temp, step_x_gmp);
								mpf_mul(newCenterX_gmp, step_x_gmp, pixel_x_gmp);
								mpf_add(newCenterX_gmp, newCenterX_gmp, f->xmin_gmp);
								
								mpf_div(step_y_gmp, range_y_temp, step_y_gmp);
								mpf_mul(newCenterY_gmp, step_y_gmp, pixel_y_gmp);
								mpf_add(newCenterY_gmp, newCenterY_gmp, f->ymin_gmp);
								
								mpf_clear(pixel_x_gmp);
								mpf_clear(pixel_y_gmp);
								mpf_clear(step_x_gmp);
								mpf_clear(step_y_gmp);
								mpf_clear(range_x_temp);
								mpf_clear(range_y_temp);
							} else {
								// Si le clic est hors zone, unzoomer depuis le centre actuel
								mpf_add(newCenterX_gmp, f->xmin_gmp, f->xmax_gmp);
								mpf_div(newCenterX_gmp, newCenterX_gmp, two_gmp);
								mpf_add(newCenterY_gmp, f->ymin_gmp, f->ymax_gmp);
								mpf_div(newCenterY_gmp, newCenterY_gmp, two_gmp);
							}
							
							// Calculer la nouvelle plage : multiplier par le facteur de zoom
							mpf_sub(range_x_gmp, f->xmax_gmp, f->xmin_gmp);
							mpf_sub(range_y_gmp, f->ymax_gmp, f->ymin_gmp);
							mpf_mul(range_x_gmp, range_x_gmp, zoomfactor_gmp);
							mpf_mul(range_y_gmp, range_y_gmp, zoomfactor_gmp);
							mpf_div(half_range_x_gmp, range_x_gmp, two_gmp);
							mpf_div(half_range_y_gmp, range_y_gmp, two_gmp);
							
							// Centrer la nouvelle vue sur le point cliqué
							mpf_add(newxmax_gmp, newCenterX_gmp, half_range_x_gmp);
							mpf_sub(newxmin_gmp, newCenterX_gmp, half_range_x_gmp);
							mpf_add(newymax_gmp, newCenterY_gmp, half_range_y_gmp);
							mpf_sub(newymin_gmp, newCenterY_gmp, half_range_y_gmp);
							
							// Mettre à jour les coordonnées GMP
							mpf_set(f->xmax_gmp, newxmax_gmp);
							mpf_set(f->xmin_gmp, newxmin_gmp);
							mpf_set(f->ymax_gmp, newymax_gmp);
							mpf_set(f->ymin_gmp, newymin_gmp);
							
							// Synchroniser vers double pour l'affichage
							f->xmax = mpf_get_d(newxmax_gmp);
							f->xmin = mpf_get_d(newxmin_gmp);
							f->ymax = mpf_get_d(newymax_gmp);
							f->ymin = mpf_get_d(newymin_gmp);
							
							printf ("UnZoom x%d  to X= %f\tY= %f\n", f->zoomfactor,
								mpf_get_d(newCenterX_gmp), mpf_get_d(newCenterY_gmp));
							
							mpf_clear(newCenterX_gmp);
							mpf_clear(newCenterY_gmp);
							mpf_clear(newxmax_gmp);
							mpf_clear(newxmin_gmp);
							mpf_clear(newymax_gmp);
							mpf_clear(newymin_gmp);
							mpf_clear(range_x_gmp);
							mpf_clear(range_y_gmp);
							mpf_clear(half_range_x_gmp);
							mpf_clear(half_range_y_gmp);
							mpf_clear(zoomfactor_gmp);
							mpf_clear(two_gmp);
							
							precision_update_fractal(f);
							// Mettre à jour la précision des coordonnées GMP si nécessaire
							if (f->gmp_precision != prec) {
								mpf_set_prec(f->xmin_gmp, f->gmp_precision);
								mpf_set_prec(f->xmax_gmp, f->gmp_precision);
								mpf_set_prec(f->ymin_gmp, f->gmp_precision);
								mpf_set_prec(f->ymax_gmp, f->gmp_precision);
							}
							// Réallouer zmatrix_gmp si nécessaire
							if (f->use_gmp && f->zmatrix_gmp == NULL) {
								f->zmatrix_gmp = (complex_gmp *) malloc ((f->xpixel*f->ypixel)* sizeof (complex_gmp));
								if (f->zmatrix_gmp == NULL) {
									fprintf(stderr, "Erreur allocation mémoire GMP fractale\n");
								} else {
									for (int k = 0; k < f->xpixel * f->ypixel; k++) {
										complex_gmp_init(&f->zmatrix_gmp[k], f->gmp_precision);
									}
								}
							} else if (!f->use_gmp && f->zmatrix_gmp != NULL) {
								for (int k = 0; k < f->xpixel * f->ypixel; k++) {
									complex_gmp_clear(&f->zmatrix_gmp[k]);
								}
								free(f->zmatrix_gmp);
								f->zmatrix_gmp = NULL;
							}
						} else {
#endif
						double newCenterX, newCenterY, newxmax, newxmin, newymax,
							newymin;
						double range_x, range_y, half_range_x, half_range_y;
						
						// Calculer le centre du zoom à partir de la position du clic
						// Vérifier que le clic est dans la zone de la fractale
						if (event->button.x >= 0 && event->button.x < f->xpixel &&
						    event->button.y >= win.y1 && event->button.y < win.y2) {
							newCenterX =
								event->button.x * ((f->xmax - f->xmin) / f->xpixel) +
								f->xmin;
							newCenterY =
								(event->button.y - win.y1) * ((f->ymax - f->ymin) / f->ypixel) +
								f->ymin;
						} else {
							// Si le clic est hors zone, unzoomer depuis le centre actuel
							newCenterX = (f->xmin + f->xmax) / 2.0;
							newCenterY = (f->ymin + f->ymax) / 2.0;
						}
						
						printf ("UnZoom x%d  to X= %f\tY= %f\n", f->zoomfactor,
							newCenterX, newCenterY);
						
						// Calculer la nouvelle plage : multiplier par le facteur de zoom
						range_x = (f->xmax - f->xmin) * f->zoomfactor;
						range_y = (f->ymax - f->ymin) * f->zoomfactor;
						half_range_x = range_x / 2.0;
						half_range_y = range_y / 2.0;
						
						// Centrer la nouvelle vue sur le point cliqué
						newxmax = newCenterX + half_range_x;
						newxmin = newCenterX - half_range_x;
						newymax = newCenterY + half_range_y;
						newymin = newCenterY - half_range_y;
						
						f->xmax = newxmax;
						f->ymax = newymax;
						f->xmin = newxmin;
						f->ymin = newymin;
#ifdef HAVE_GMP
						precision_update_fractal(f);
						// Synchroniser vers GMP si nécessaire
						if (f->use_gmp) {
							mpf_set_d(f->xmin_gmp, f->xmin);
							mpf_set_d(f->xmax_gmp, f->xmax);
							mpf_set_d(f->ymin_gmp, f->ymin);
							mpf_set_d(f->ymax_gmp, f->ymax);
						}
						// Réallouer zmatrix_gmp si nécessaire
						if (f->use_gmp && f->zmatrix_gmp == NULL) {
							f->zmatrix_gmp = (complex_gmp *) malloc ((f->xpixel*f->ypixel)* sizeof (complex_gmp));
							if (f->zmatrix_gmp == NULL) {
								fprintf(stderr, "Erreur allocation mémoire GMP fractale\n");
							} else {
								for (int k = 0; k < f->xpixel * f->ypixel; k++) {
									complex_gmp_init(&f->zmatrix_gmp[k], f->gmp_precision);
								}
							}
						} else if (!f->use_gmp && f->zmatrix_gmp != NULL) {
							for (int k = 0; k < f->xpixel * f->ypixel; k++) {
								complex_gmp_clear(&f->zmatrix_gmp[k]);
							}
							free(f->zmatrix_gmp);
							f->zmatrix_gmp = NULL;
						}
						}
#endif
					}
					// Buddhabrot utilise son propre algorithme de rendu
					if (*typeFractale == 16) {
						*renderTime = Buddhabrot_Draw (screen, f, 0, win.y1, g);
					} else {
						*renderTime = Fractal_Draw (screen, *f, 0,win.y1, g);
					}
					centerX = (f->xmin + f->xmax) / 2;
					centerY = (f->ymin + f->ymax) / 2;
					if (g != NULL)
						SDLGUI_StateBar_Update(screen, g, *typeFractale, f->colorMode, centerX, centerY, f->zoomfactor, *renderTime, f);
					SDL_UpdateRect (screen, 0, 0, screen->w, screen->h);

				}
				break;

			}
		}

		
			default:    // Default General du switch de traitement d'evenement
				{
					//printf ("Event non traite pour le moment\n");
					return 1;
				}
				break;
      }
	return 0;
    }
