 /* ExcapeTime.h
   fractal formula header
   released under GPL2
   Copyleft 2001-2003 VERHILLE Arnaud
*/

#include "SDL.h"
#include "complexmath.h"

// A portable color type.
typedef struct {
  int r, g, b, a;
} color;

// Create a new fractal type
typedef struct {
  int xpixel, ypixel;
  double xmin, xmax, ymin, ymax;
  complex seed;
  int iterationMax;
  int bailout;
  int zoomfactor;
  int type;
  int colorMode;   // 0=Normal, 1=Monochrome, 2=Fire, 3=Ocean
  int *fmatrix;    // la matrice d'iteration
  complex *zmatrix;  // la matrice de la valeur de z a la derniere iteration
  color *cmatrix; // Une matrice de couleurs, soit la fractale finale.
} fractal;

// Fractal point Calcul result
typedef struct {
  int iteration;
  complex z;
} fractalresult;


// ******************
// Interface publique
// ******************

// Initialisation et destruction
 fractal Fractal_Init (int screenW, int screenH, int type);
 void Fractal_Destroy (fractal);

 // Calcul de la matrice de complexe z lors de leur derniere iteration
 void Fractal_CalculateMatrix (fractal*);
 void Fractal_CalculateMatrix_DDp1 (fractal*);
 void Fractal_CalculateMatrix_DDp2 (fractal*);

 // Calcul de la couleur
 void Fractal_CalculateColorMatrix (fractal*); // Selecteur
 void FractalColorMonochrome (fractal*);
 void FractalColorNormal (fractal*);
 void FractalColorTest (fractal*);

// Color Formulae Utilities
color Fractal_ReadColorMatrix (fractal, int, int);
int Fractal_ReadColorMatrixRed (fractal, int, int);
int Fractal_ReadColorMatrixGreen (fractal, int, int);
int Fractal_ReadColorMatrixBlue (fractal, int, int);

Uint32 Fractal_Draw (SDL_Surface*, fractal, int,int);
void Fractal_ChangeType (fractal* f, int type);

// Formulae Utilities
 fractalresult FormulaSelector (fractal f, complex zPixel);
 fractalresult Mendelbrot_Iteration (fractal, complex);
 fractalresult Julia_Iteration (fractal, complex);
 fractalresult JuliaSin_Iteration (fractal, complex);
 fractalresult Newton_Iteration (fractal, complex);
 fractalresult Phoenix_Iteration (fractal, complex);
 fractalresult Sierpinski_Iteration (fractal, complex);
fractalresult Barnsleyj1_Iteration (fractal, complex);
fractalresult Barnsleym1_Iteration (fractal, complex);
fractalresult BurningShip_Iteration (fractal, complex);
fractalresult Tricorn_Iteration (fractal, complex);
fractalresult Mandelbulb_Iteration (fractal, complex);

 // Fractal Definition
 void Mendelbrot_def (fractal* f);
 void Julia_def (fractal* f);
 void JuliaSin_def (fractal* f);
 void Newton_def (fractal* f);
 void Phoenix_def (fractal* f);
 void Sierpinski_def (fractal* f);
 void Barnsley1j_def (fractal* f);
void Barnsley1m_def (fractal* f);
void Magnet1j_def (fractal* f);
void Magnet1m_def (fractal* f);
void BurningShip_def (fractal* f);
void Tricorn_def (fractal* f);
void Mandelbulb_def (fractal* f);
void Buddhabrot_def (fractal* f);

// Buddhabrot special draw function (density algorithm)
// gui parameter can be NULL if no GUI progress display needed
Uint32 Buddhabrot_Draw (SDL_Surface*, fractal*, int, int, void* gui);


