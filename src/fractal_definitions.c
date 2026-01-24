/* fractal_definitions.c
   Fonctions de definition des parametres par defaut
   released under GPL2
   Copyleft 2001-2026 VERHILLE Arnaud
*/

#include "fractal_definitions.h"
#include "complexmath.h"
#include <math.h>

void Mendelbrot_def(fractal* f) {
	/* Valeurs de base pour la Mandelbrot */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 9370;
	f->zoomfactor = 8;
}

void Julia_def(fractal* f) {
	/* Valeurs de base pour la Julia */
	f->xmin = -2.0;
	f->xmax = +2.0;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = MakeComplex(0.36228, -0.0777);
	f->bailout = 4;
	f->iterationMax = 6250;
	f->zoomfactor = 4;
}

void JuliaSin_def(fractal* f) {
	/* Valeurs de base pour la JuliaSin */
	f->xmin = -M_PI;
	f->xmax = M_PI;
	f->ymin = -2;
	f->ymax = 2;
	f->seed = MakeComplex(1, 0.1);
	f->bailout = 4;
	f->iterationMax = 6250;
	f->zoomfactor = 4;
}

void Newton_def(fractal* f) {
	/* Valeurs de base pour la Newton */
	f->seed = MakeComplex(8, 0);
	f->xmin = -3;
	f->xmax = 3;
	f->ymin = -2;
	f->ymax = 2;
	f->bailout = 4;
	f->iterationMax = 1000;
	f->zoomfactor = 4;
}

void Phoenix_def(fractal* f) {
	/* Valeurs de base pour la Phoenix */
	f->xmin = -2.0;
	f->xmax = +2.0;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->bailout = 4;
	f->iterationMax = 9370;
	f->zoomfactor = 8;
}

void Barnsley1j_def(fractal* f) {
	/* Valeurs de base pour la Barnsley1j */
	f->xmin = -4;
	f->xmax = 4;
	f->ymin = -3;
	f->ymax = 3;
	f->seed = MakeComplex(1.1, 0.6);
	f->bailout = 4;
	f->iterationMax = 3120;
	f->zoomfactor = 4;
}

void Barnsley1m_def(fractal* f) {
	/* Valeurs de base pour la Barnsley1m */
	f->xmin = -3;
	f->xmax = 3;
	f->ymin = -2;
	f->ymax = 2;
	f->bailout = 4;
	f->iterationMax = 9370;
	f->zoomfactor = 8;
}

void Magnet1j_def(fractal* f) {
	/* Valeurs de base pour la Magnet1j */
	f->seed = MakeComplex(1.625458, -0.306159);
	f->xmin = -2;
	f->xmax = 2;
	f->ymin = -2;
	f->ymax = 2;
	f->bailout = 4;
	f->iterationMax = 9370;
	f->zoomfactor = 4;
}

void Magnet1m_def(fractal* f) {
	/* Valeurs de base pour la Magnet1m */
	f->xmin = -3;
	f->xmax = 2;
	f->ymin = -2;
	f->ymax = 2;
	f->bailout = 4;
	f->iterationMax = 9370;
	f->zoomfactor = 8;
}

void BurningShip_def(fractal* f) {
	/* Valeurs de base pour Burning Ship */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -2.0;
	f->ymax = 2.0;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 9370;
	f->zoomfactor = 8;
}

void Buffalo_def(fractal* f) {
	/* Valeurs de base pour Buffalo */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -2.0;
	f->ymax = 2.0;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 9370;
	f->zoomfactor = 8;
}

void Tricorn_def(fractal* f) {
	/* Valeurs de base pour Tricorn */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 9370;
	f->zoomfactor = 8;
}

void Mandelbulb_def(fractal* f) {
	/* Valeurs de base pour Mandelbulb */
	f->xmin = -1.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 9370;
	f->zoomfactor = 8;
}

void PerpendicularBurningShip_def(fractal* f) {
	/* Valeurs de base pour Perpendicular Burning Ship */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 5000;
	f->zoomfactor = 8;
}

void Celtic_def(fractal* f) {
	/* Valeurs de base pour Celtic Fractal */
	f->xmin = -2.0;
	f->xmax = 1.0;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 5000;
	f->zoomfactor = 8;
}

void AlphaMandelbrot_def(fractal* f) {
	/* Valeurs de base pour Alpha Mandelbrot */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 2000;
	f->zoomfactor = 8;
}

void PickoverStalks_def(fractal* f) {
	/* Valeurs de base pour Pickover Stalks */
	f->xmin = -2.0;
	f->xmax = 1.0;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 100.0;
	f->iterationMax = 1000;
	f->zoomfactor = 8;
}

void Nova_def(fractal* f) {
	/* Valeurs de base pour Nova Fractal */
	f->xmin = -3.0;
	f->xmax = 3.0;
	f->ymin = -2.0;
	f->ymax = 2.0;
	f->seed = ZeroSetofComplex();
	f->bailout = 20.0;
	f->iterationMax = 500;
	f->zoomfactor = 4;
}

void Multibrot_def(fractal* f) {
	/* Valeurs de base pour Multibrot (puissances non-entieres) */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 5000;
	f->zoomfactor = 8;
}

void Buddhabrot_def(fractal* f) {
	/* Valeurs de base pour Buddhabrot */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 220;
	f->zoomfactor = 4;
}

void Lyapunov_def(fractal* f) {
	/* Zircon City : Lyapunov fractal pour la sequence bbbbbbaaaaaa
	 * Domaine : (a, b) in [2.5, 3.4] x [3.4, 4.0] */
	f->xmin = 2.5;
	f->xmax = 3.4;
	f->ymin = 3.4;
	f->ymax = 4.0;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 2000;
	f->zoomfactor = 2;
}

void Nebulabrot_def(fractal* f) {
	/* Valeurs de base pour Nebulabrot */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 5000;
	f->zoomfactor = 4;
}
