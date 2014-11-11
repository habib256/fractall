/* VonKoch.c
   fractal formula src code
   released under GPL2
   Copyleft 2003 VERHILLE Arnaud
*/

#include <math.h>
#include "SDL_gfxPrimitives.h"
#include "VonKoch.h"


void DragonFractDraw (SDL_Surface* screen, int x1, int y1, int x2, int y2,int iteration) {
	struct point D, E;

	// Calcul de la position des points de depart
	D.x = ((0.338125) * (x2-x1))+x1;
	D.y = (0.208333 * (y2-y1))+y1;
	E.x = ((0.838125) * (x2-x1))+x1;
	//E.y = (0.708333 * (y2-y1));
	E.y = y2 -(0.31 * (y2-y1));

	boxRGBA(screen,(Sint16)x1,(Sint16)y1,(Sint16)x2, (Sint16)y2,0,0,0, 255);
	DragonFractRec (screen, D, E, iteration);
	SDL_UpdateRect (screen, x1, y1, x2-x1, y2-y1);
	}

void VonKochDraw (SDL_Surface* screen, int x1, int y1, int x2, int y2, int iteration) {
	struct point A, B, C;

	A.x = (0.5 * (x2-x1))+x1;
	A.y = (0.020833 * (y2-y1))+y1;
	B.x = (0.203125 * (x2-x1))+x1;
	B.y = ((0.708333) * (y2-y1))+y1;
	C.x = (0.796875 * (x2-x1))+x1;
	C.y = B.y;
	
	boxRGBA(screen,(Uint16)x1,(Uint16)y1,(Uint16)x2, (Uint16)y2,0,0,0, 255);
	VonKochRec (screen, A, B, iteration);
	VonKochRec (screen, C, A, iteration);
	VonKochRec (screen, B, C, iteration);
	SDL_UpdateRect (screen, x1, y1, x2-x1, y2-y1);
}

void DragonFractRec (SDL_Surface *screen, struct point A, struct point B, int n) {
	struct point C,D;
	
	if (n==0) {
		//Cas d'arret, on trace une ligne
	aalineRGBA(screen,(Sint16) A.x, (Sint16)A.y, (Sint16)B.x,(Sint16)B.y, 0xFF, 0xFF, 0xFF, 0xFF);
	} else {
		// On cree un nouveau point C qui est la rotation de 90°
		// de B par rapport au centre de [AB]
		// puis on applique recursivement la procedure sur les 2 nouveaux segments
		C.x = (A.x + B.x) / 2.0;
		C.y = (A.y + B.y) / 2.0;
		D = rotation (B, M_PI/2.0, C);
		DragonFractRec (screen, A, D, n-1);
		DragonFractRec (screen, B, D, n-1);
	}
}

// Trace du flocon de VonKoch
// A et B sont les deux extremites du segment
// et n le nombre d'iterations
void VonKochRec (SDL_Surface *screen,struct point A, struct point B, int n) {
	struct point C,D,E;
	
	if (n == 0) {
		//Cas d'arret, on trace une ligne
		aalineRGBA(screen,(Sint16) A.x, (Sint16)A.y,(Sint16)B.x,(Sint16)B.y, 0xFF, 0xFF, 0xFF, 0xFF);

	} else {
		C.x = A.x + ((B.x - A.x)/3.0);
		C.y = A.y + ((B.y - A.y)/3.0);
		// On obtient le point C = A + 1/3Vect(AB)
		D.x = C.x + ((B.x - A.x)/3.0);
		D.y = C.y + ((B.y - A.y)/3.0);
		// Obtention du point D = C + 1/3Vect(AB)
		E = rotation (B,(2.0*M_PI)/3.0, D);
		// E est la rotation de centre D, du point B et d'angle 120°
		VonKochRec(screen,A,C,n-1);
		VonKochRec(screen,C,E,n-1); // Appels recursifs sur les nouveaux segments crees
		VonKochRec(screen,E,D,n-1);
		VonKochRec(screen,D,B,n-1);
	}
}


// Fait la rotation du point p autour de "centre", d'angle theta (en radian)
struct point rotation (struct point p,double theta, struct point centre) {
	struct point Newp;
	struct point Temp;

	// Temp va contenir p, apres avoir centre la rotation sur (0,0)
	Temp.x = p.x - centre.x;
	Temp.y = p.y - centre.y;

	// On effectue la rotation autour de (0,0), apres produit avec la matrice rotation
	Newp.x = Temp.x * cos(theta) - Temp.y * sin(theta);
	Newp.y = Temp.x * sin(theta) + Temp.y * cos(theta);

	// Enfin, on retourne le point apres rotation, plus un decalage pour revenir dans le repere initial
	Newp.x = Newp.x + centre.x;
	Newp.y = Newp.y + centre.y;

	return Newp;
}
