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
#include "SDLGUI.h"
#include "SDL_gfxPrimitives.h"
#include "EscapeTime.h"
#include "VonKoch.h"



#define PACKAGE "fractall"
#define VERSION 0.5




typedef struct { // Conteneur graphique en 3 parties
	int y1, y2;
	int w,h;
} window;


// FUNCTIONS PROTOTYPES
//*********************

int EventFlush (SDL_Event*);
int EventCheck (SDL_Event*, SDL_Surface *, gui *, fractal *, int*, int*, window win);


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
	int screenW = 450, screenH = 400;	// Largeur et hauteur de la fenetre
	int menuH = 51;			// C'est la place que prend le gui
	int stateH = screenH - 20;		// Place Barre d etat
	int typeFractale = 1;		//Type de fractale par defaut
	int iteration = 1;   		// Nbr d iteration de la fractale par defaut
	
	SDL_Event event;		// Gestion des evenements
	
	fractal f;			// Allocate a new fractal.
	gui g;			// Allocate a new gui
	window win;		// Creation d un  nouveau conteneur
	
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
		g = SDLGUI_Init (0, 0, win.w, win.y1, stateH, SDL_MapRGB (screen->format, 213, 214, 213), 14);	// 8 bouttons par exemple
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
		Fractal_Draw (screen, f, 0,win.y1);
		break;
    }
    
    SDL_UpdateRect (screen, 0, 0, screen->w, screen->h);

	//MAIN INFINITE BOUCLE
	for (;;)
    {
		EventFlush (&event);	// Event when Initializing
		while (SDL_WaitEvent (&event) >= 0)
		{
			EventCheck (&event, screen, &g, &f, &typeFractale, &iteration,win);
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
				int* typeFractale, int* iteration, window win)
{	
	switch (event->type)
    {
		
		// Gestion des touches du clavier
		// ******************************
		
    case SDL_KEYDOWN:
		{
			if ((event->key.keysym.sym == SDLK_ESCAPE) | (event->key.keysym.sym ==
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
				Fractal_Draw (screen, *f, 0,win.y1);
			}
			if (event->key.keysym.sym == SDLK_F4)
			{
				*typeFractale = 4;
				Fractal_ChangeType (f, 4);
				Fractal_Draw (screen, *f,0, win.y1);
			}
			if (event->key.keysym.sym == SDLK_F5)
			{
				*typeFractale = 5;
				Fractal_ChangeType (f, 5);
				Fractal_Draw (screen, *f, 0,win.y1);
			}
			if (event->key.keysym.sym == SDLK_F6)
			{
				*typeFractale = 6;
				Fractal_ChangeType (f, 6);
				Fractal_Draw (screen, *f, 0,win.y1);
			}
			if (event->key.keysym.sym == SDLK_F7)
			{
				*typeFractale = 7;
				Fractal_ChangeType (f, 7);
				Fractal_Draw (screen, *f, 0,win.y1);
			}
			if (event->key.keysym.sym == SDLK_F8)
			{
				*typeFractale = 8;
				Fractal_ChangeType (f, 8);
				Fractal_Draw (screen, *f, 0,win.y1);
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
			if (g == NULL)
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
					default:
					{
					Fractal_ChangeType (f, *typeFractale);
					Fractal_Draw (screen, *f, 0,win.y1);
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
						double newCenterX, newCenterY, newxmax, newxmin, newymax,
							newymin;
						
						newCenterX =
							event->button.x * ((f->xmax - f->xmin) / f->xpixel) +
							f->xmin;
						newCenterY =
							(event->button.y - win.y1) * ((f->ymax - f->ymin) / f->ypixel) +
							f->ymin;
						printf ("Zoom x%d  to X= %f\tY= %f\n", f->zoomfactor,
							newCenterX, newCenterY);
						newxmax =
							(f->xmax - f->xmin) / (2 + f->zoomfactor) + newCenterX;
						newxmin =
							-(f->xmax - f->xmin) / (2 + f->zoomfactor) + newCenterX;
						newymax =
							(f->ymax - f->ymin) / (2 + f->zoomfactor) + newCenterY;
						newymin =
							-(f->ymax - f->ymin) / (2 + f->zoomfactor) + newCenterY;
						f->xmax = newxmax;
						f->ymax = newymax;
						f->xmin = newxmin;
						f->ymin = newymin;
					}
					if (event->button.button == SDL_BUTTON_RIGHT)
					{		// UNZOOM
						double newCenterX, newCenterY, newxmax, newxmin, newymax,
							newymin;
						
						newCenterX =
							event->button.x * ((f->xmax - f->xmin) / f->xpixel) +
							f->xmin;
						newCenterY =
							(event->button.y) * ((f->ymax - f->ymin) / f->ypixel) + f->ymin;
						printf ("UnZoom x%d  to X= %f\tY= %f\n", f->zoomfactor,
							newCenterX, newCenterY);
						newxmax =
							(f->xmax - f->xmin) * (f->zoomfactor / 2) + newCenterX;
						newxmin =
							-(f->xmax - f->xmin) * (f->zoomfactor / 2) + newCenterX;
						newymax =
							(f->ymax - f->ymin) * (f->zoomfactor / 2) + newCenterY;
						newymin =
							-(f->ymax - f->ymin) * (f->zoomfactor / 2) + newCenterY;
						f->xmax = newxmax;
						f->ymax = newymax;
						f->xmin = newxmin;
						f->ymin = newymin;
					}
					Fractal_Draw (screen, *f, 0,win.y1);
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
