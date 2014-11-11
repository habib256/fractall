/* Header file for constructing VonKoch like Vectoriel Fractals
   released under GPL2
   Copyleft 2003 VERHILLE Arnaud
*/

#include "SDL.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846 /* pi */
#endif  /* POUR microsoft VC6 */

struct point {
	double x;
	double y;
};

// Function Prototype
// *******************

// Les fonctions recursives
void DragonFractRec (SDL_Surface*,struct point, struct point, int);
void VonKochRec (SDL_Surface*,struct point, struct point, int);

// Mathematique geometrique
struct point rotation(struct point,double, struct point);

// Trace des fonctions de VonKoch
void DragonFractDraw (SDL_Surface*,int,int,int,int, int);
void VonKochDraw (SDL_Surface*, int,int,int,int,int );
