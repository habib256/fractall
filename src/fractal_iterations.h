/* fractal_iterations.h
   Declarations des fonctions d'iteration (double precision)
   released under GPL2
   Copyleft 2001-2026 VERHILLE Arnaud
*/

#ifndef FRACTAL_ITERATIONS_H
#define FRACTAL_ITERATIONS_H

#include "fractal_types.h"

// Selecteur de formule (point d'entree principal)
fractalresult FormulaSelector(fractal f, complex zPixel);

// Fonctions d'iteration individuelles
fractalresult Mendelbrot_Iteration(fractal f, complex zPixel);
fractalresult Julia_Iteration(fractal f, complex zPixel);
fractalresult JuliaSin_Iteration(fractal f, complex zPixel);
fractalresult Newton_Iteration(fractal f, complex zPixel);
fractalresult Phoenix_Iteration(fractal f, complex zPixel);
fractalresult Buffalo_Iteration(fractal f, complex zPixel);
fractalresult Barnsleyj1_Iteration(fractal f, complex zPixel);
fractalresult Barnsleym1_Iteration(fractal f, complex zPixel);
fractalresult Magnet1j_Iteration(fractal f, complex zPixel);
fractalresult Magnet1m_Iteration(fractal f, complex zPixel);
fractalresult BurningShip_Iteration(fractal f, complex zPixel);
fractalresult Tricorn_Iteration(fractal f, complex zPixel);
fractalresult Mandelbulb_Iteration(fractal f, complex zPixel);
fractalresult PerpendicularBurningShip_Iteration(fractal f, complex zPixel);
fractalresult Celtic_Iteration(fractal f, complex zPixel);
fractalresult AlphaMandelbrot_Iteration(fractal f, complex zPixel);
fractalresult PickoverStalks_Iteration(fractal f, complex zPixel);
fractalresult Nova_Iteration(fractal f, complex zPixel);
fractalresult Multibrot_Iteration(fractal f, complex zPixel);

#endif /* FRACTAL_ITERATIONS_H */
