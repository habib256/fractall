/* fractal_iterations_gmp.h
   Declarations des fonctions d'iteration GMP (precision arbitraire)
   released under GPL2
   Copyleft 2001-2026 VERHILLE Arnaud
*/

#ifndef FRACTAL_ITERATIONS_GMP_H
#define FRACTAL_ITERATIONS_GMP_H

#include "config.h"

#ifdef HAVE_GMP

#include "fractal_types.h"
#include "complexmath_gmp.h"

// Selecteur de formule GMP
fractalresult FormulaSelector_GMP(fractal f, complex_gmp zPixel);

// Fonctions d'iteration GMP
fractalresult Mendelbrot_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Julia_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult JuliaSin_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Newton_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Phoenix_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Buffalo_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Barnsleyj1_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Barnsley1m_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Magnet1j_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Magnet1m_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult BurningShip_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Tricorn_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Mandelbulb_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult PerpendicularBurningShip_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Celtic_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult AlphaMandelbrot_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult PickoverStalks_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Nova_Iteration_GMP(fractal f, complex_gmp zPixel);
fractalresult Multibrot_Iteration_GMP(fractal f, complex_gmp zPixel);

#endif /* HAVE_GMP */

#endif /* FRACTAL_ITERATIONS_GMP_H */
