/* fractal_iterations.c
   Fonctions d'iteration double precision
   released under GPL2
   Copyleft 2001-2026 VERHILLE Arnaud
*/

#include "fractal_iterations.h"
#include "config.h"
#include <math.h>

#ifdef HAVE_GMP
#include "fractal_iterations_gmp.h"
#endif

// Calcul de la Mandelbrot
fractalresult Mendelbrot_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;

	z = f.seed;
	do { //z(n+1) = z(n)^2 + pixel
		complex zTemp;
		i++;
		zTemp = Mulz(z, z); // z(n)^2
		z = Addz(zTemp, zPixel); // Ajout de pixel
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de la Julia
fractalresult Julia_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;
	z = zPixel;
	do { //z(n+1) = z(n)^2 + seed
		complex zTemp;
		i++;
		zTemp = Mulz(z, z); // z(n)^2
		z = Addz(zTemp, f.seed); // ajout de seed
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de JuliaSin
fractalresult JuliaSin_Iteration(fractal f, complex zPixel) {
	/* JuliaSin
	 * z(0) = pixel;
	 * z(k+1) = c * sin(z(k))
	 */
	int i = 0;
	complex z;
	fractalresult result;
	z = zPixel;
	do {
		i++;
		z = Mulz(f.seed, sinz(z));
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Newton
fractalresult Newton_Iteration(fractal f, complex zPixel) {
	/*
	 * z(0) = pixel;
	 * z(n+1) = ((p-1)*z(n)^p + 1)/(p*z(n)^(p - 1)).
	 */
	int i = 0;
	complex z;
	fractalresult result;
	z = zPixel;
	do {
		int p = Rez(f.seed);  // Degree polynomial
		complex zQuot, zNum;
		i++;
		zNum = Addz((Mulz(MakeComplex(p-1, 0.0), Powz(z, MakeComplex(p, 0.0)))), MakeComplex(1.0, 0.0));
		zQuot = Mulz(Powz(z, MakeComplex(p-1, 0.0)), MakeComplex(p, 0.0));
		z = Divz(zNum, zQuot);
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Phoenix (degree 0 uniquement)
fractalresult Phoenix_Iteration(fractal f, complex zPixel) {
	/*
	 * z(0) = pixel, y(0) = 0;
	 * For degree = 0: z(n+1) = z(n)^2 + p1.x + p2.x*y(n), y(n+1) = z(n)
	 */
	int i = 0;
	complex z, y;
	complex zTemp, zpTemp;
	fractalresult result;
	int degree = 0;
	double p1 = 0.56667;
	double p2 = -0.5;

	z = zPixel;
	y = ZeroSetofComplex();
	do {
		i++;
		if (degree == 0) {
			zTemp = Mulz(z, z);
			zTemp = MakeComplex(Rez(zTemp) + p1, Imz(zTemp));
			zpTemp = ScalarTimesofComplex(p2, y);
			zTemp = Addz(zTemp, zpTemp);
			y = z;
			z = zTemp;
		}
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Barnsley1j
fractalresult Barnsleyj1_Iteration(fractal f, complex zPixel) {
	/* Barnsleyj1
	 * z(0) = pixel;
	 * if real(z) >= 0
	 * z(n+1) = (z-1)*c
	 * else
	 * z(n+1) = (z+1)*c */
	int i = 0;
	complex z, zTemp;
	fractalresult result;
	z = zPixel;
	do {
		i++;
		if (Rez(z) >= 0) {
			zTemp = MakeComplex((Rez(z) - 1), Imz(z));
			z = Mulz(zTemp, f.seed);
		} else {
			zTemp = MakeComplex((Rez(z) + 1), Imz(z));
			z = Mulz(zTemp, f.seed);
		}
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Barnsley1m
fractalresult Barnsleym1_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z, c, zTemp;
	fractalresult result;
	z = zPixel;
	c = zPixel;
	do {
		i++;
		if (Rez(z) >= 0) {
			zTemp = MakeComplex((Rez(z) - 1), Imz(z));
			z = Mulz(zTemp, c);
		} else {
			zTemp = MakeComplex((Rez(z) + 1), Imz(z));
			z = Mulz(zTemp, c);
		}
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Magnet1j
fractalresult Magnet1j_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z, N, Q;
	fractalresult result;
	z = zPixel;
	do {
		i++;
		N = Addz(Mulz(z, z), MakeComplex(Rez(f.seed) - 1, Imz(f.seed)));
		N = Mulz(N, N);
		Q = Addz(Mulz(MakeComplex(2, 0), z), MakeComplex(Rez(f.seed) - 2, Imz(f.seed)));
		z = Divz(N, Q);
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Magnet1m
fractalresult Magnet1m_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z, c, N, Q;
	fractalresult result;
	c = zPixel;
	z = ZeroSetofComplex();
	do {
		i++;
		N = Addz(Mulz(z, z), MakeComplex(Rez(c) - 1, Imz(c)));
		N = Mulz(N, N);
		Q = Addz(Mulz(MakeComplex(2, 0), z), MakeComplex(Rez(c) - 2, Imz(c)));
		z = Divz(N, Q);
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Burning Ship
fractalresult BurningShip_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;

	z = f.seed;
	do {
		complex zTemp;
		double re, im;
		i++;
		// z(n+1) = (|Re(z)| + i|Im(z)|)^2 + c
		re = fabs(Rez(z));
		im = fabs(Imz(z));
		zTemp = MakeComplex(re, im);
		zTemp = Mulz(zTemp, zTemp); // (|Re(z)| + i|Im(z)|)^2
		z = Addz(zTemp, zPixel); // Ajout de pixel
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Buffalo
fractalresult Buffalo_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;

	z = f.seed;
	do {
		complex zSq;
		double re_sq, im_sq;
		i++;
		// z^2 = (a+bi)^2 = (a^2-b^2) + 2abi
		zSq = Mulz(z, z);
		// Appliquer abs() aux deux parties
		re_sq = fabs(Rez(zSq));
		im_sq = fabs(Imz(zSq));
		z = MakeComplex(re_sq + Rez(zPixel), im_sq + Imz(zPixel));
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Tricorn
fractalresult Tricorn_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;

	z = f.seed;
	do {
		complex zTemp, zConj;
		i++;
		// z(n+1) = (conjugue(z))^2 + c
		zConj = MakeComplex(Rez(z), -Imz(z)); // Conjugue
		zTemp = Mulz(zConj, zConj); // (conjugue(z))^2
		z = Addz(zTemp, zPixel); // Ajout de pixel
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Mandelbulb (version 2D avec puissance 8)
fractalresult Mandelbulb_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;

	z = f.seed;
	do {
		complex zTemp;
		i++;
		// z(n+1) = z(n)^8 + c
		// Calcul de z^8 par multiplications successives
		zTemp = Mulz(z, z);       // z^2
		zTemp = Mulz(zTemp, zTemp); // z^4
		zTemp = Mulz(zTemp, zTemp); // z^8
		z = Addz(zTemp, zPixel); // z^8 + c
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Perpendicular Burning Ship
fractalresult PerpendicularBurningShip_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;

	z = f.seed;
	do {
		double x, y, y_abs;
		double x2, y2;
		i++;
		// z(n+1) = (Re(z) - i*|Im(z)|)^2 + c
		x = Rez(z);
		y = Imz(z);
		y_abs = fabs(y);

		// Calcul de (x - i*y_abs)^2 = x^2 - y_abs^2 - i*2*x*y_abs
		x2 = x * x;
		y2 = y_abs * y_abs;
		z = MakeComplex(x2 - y2 + Rez(zPixel), -2.0 * x * y_abs + Imz(zPixel));
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Celtic Fractal
fractalresult Celtic_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;

	z = f.seed;
	do {
		double x, y;
		double u, v;  // u = partie reelle de z^2, v = partie imaginaire de z^2
		i++;
		// Calcul de z^2 = (x + iy)^2 = (x^2 - y^2) + i(2xy)
		x = Rez(z);
		y = Imz(z);
		u = x * x - y * y;  // partie reelle de z^2
		v = 2.0 * x * y;    // partie imaginaire de z^2

		// Application de abs() a la partie reelle
		z = MakeComplex(fabs(u) + Rez(zPixel), v + Imz(zPixel));
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Alpha Mandelbrot (Nested/Composite)
fractalresult AlphaMandelbrot_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;

	z = f.seed;
	do {
		complex zSq, m, mSq;
		i++;
		// Calcul de z^2
		zSq = Mulz(z, z);
		// Calcul de m = z^2 + c
		m = Addz(zSq, zPixel);
		// Calcul de m^2
		mSq = Mulz(m, m);
		// z_{n+1} = z^2 + m^2 + c
		z = Addz(Addz(zSq, mSq), zPixel);
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Pickover Stalks / Biomorphs
fractalresult PickoverStalks_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;
	double trap_min = 1e10;  // Valeur initiale tres grande
	double trap_distance;
	const double trap_divisor = 0.03;  // Ajuste l'epaisseur des stalks

	z = f.seed;
	do {
		double re_abs, im_abs;
		i++;
		// Iteration Mandelbrot: z(n+1) = z^2 + c
		z = Addz(Mulz(z, z), zPixel);

		// Calcul de la distance au trap (cross: axes Re=0 et Im=0)
		re_abs = fabs(Rez(z));
		im_abs = fabs(Imz(z));
		trap_distance = (re_abs < im_abs) ? re_abs : im_abs;

		// Mise a jour de trap_min
		if (trap_distance < trap_min) {
			trap_min = trap_distance;
		}
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	// Stocker une valeur normalisee basee sur trap_min dans iteration
	if (trap_min > 1e-10) {
		double log_trap = -log(trap_min / trap_divisor);
		result.iteration = (int)(log_trap * 100.0);
		if (result.iteration < 0) result.iteration = 0;
		if (result.iteration >= f.iterationMax) result.iteration = f.iterationMax - 1;
	} else {
		// trap_min tres petit, point tres proche des axes
		result.iteration = f.iterationMax - 1;
	}

	result.z = z;
	return result;
}

// Calcul de Nova Fractal
fractalresult Nova_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z, z_prev;
	fractalresult result;
	complex a_relax = MakeComplex(1.0, 0.0);  // Parametre de relaxation
	int p = 3;  // Exposant polynomial (cubique par defaut)
	double conv_epsilon = 1e-7;  // Seuil de convergence
	double conv_epsilon_sq = conv_epsilon * conv_epsilon;

	// Pour Nova en mode Mandelbrot, utiliser z0 = 1
	z = MakeComplex(1.0, 0.0);
	z_prev = z;

	do {
		complex z_pow, z_pow_deriv, numerator, denominator, newton_step;
		i++;

		// Calcul de z^p
		z_pow = Powz(z, MakeComplex((double)p, 0.0));
		// Calcul de z^(p-1)
		z_pow_deriv = Powz(z, MakeComplex((double)(p-1), 0.0));

		// p(z) = z^p - 1
		numerator = Subz(z_pow, ESetofComplex());
		// p'(z) = p*z^(p-1)
		denominator = ScalarTimesofComplex((double)p, z_pow_deriv);

		// Eviter division par zero
		if (Magz(denominator) < 1e-10) {
			break;
		}

		// Newton step: p(z)/p'(z)
		newton_step = Divz(numerator, denominator);
		// Multiplier par a (relaxation)
		newton_step = Mulz(a_relax, newton_step);

		// z(n+1) = z - a*(p(z)/p'(z)) + c
		z_prev = z;
		z = Addz(Subz(z, newton_step), zPixel);

		// Detection de convergence
		double diff_sq = Magz2(Subz(z, z_prev));
		double z_sq = Magz2(z);
		double denom = (z_sq < 1.0) ? 1.0 : z_sq;

		if (diff_sq / denom < conv_epsilon_sq) {
			break;
		}

		// Echappement si |z| devient trop grand
		if (Magz(z) > f.bailout) {
			break;
		}
	} while (i < f.iterationMax);

	result.iteration = i;
	result.z = z;
	return result;
}

// Calcul de Multibrot (Puissances Non-entieres)
fractalresult Multibrot_Iteration(fractal f, complex zPixel) {
	int i = 0;
	complex z;
	fractalresult result;
	double d = 2.5;  // Exposant non-entier

	z = f.seed;
	do {
		complex z_pow;
		i++;

		// Calcul de z^d via exponentiation complexe
		z_pow = Powz(z, MakeComplex(d, 0.0));

		// Verifier si le resultat est valide (NaN ou Inf)
		if (isnan(Rez(z_pow)) || isnan(Imz(z_pow)) ||
		    isinf(Rez(z_pow)) || isinf(Imz(z_pow))) {
			break;
		}

		// z(n+1) = z^d + c
		z = Addz(z_pow, zPixel);
	} while ((i < f.iterationMax) && (Magz(z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}

// Selecteur de formule
fractalresult FormulaSelector(fractal f, complex zPixel) {
	fractalresult r;

#ifdef HAVE_GMP
	if (f.use_gmp) {
		complex_gmp zPixel_gmp = complex_to_gmp(zPixel, f.gmp_precision);
		r = FormulaSelector_GMP(f, zPixel_gmp);
		complex_gmp_clear(&zPixel_gmp);
		return r;
	}
#endif

	switch (f.type) {
	case 3:
		r = Mendelbrot_Iteration(f, zPixel);
		break;
	case 4:
		r = Julia_Iteration(f, zPixel);
		break;
	case 5:
		r = JuliaSin_Iteration(f, zPixel);
		break;
	case 6:
		r = Newton_Iteration(f, zPixel);
		break;
	case 7:
		r = Phoenix_Iteration(f, zPixel);
		break;
	case 8:
		r = Buffalo_Iteration(f, zPixel);
		break;
	case 9:
		r = Barnsleyj1_Iteration(f, zPixel);
		break;
	case 10:
		r = Barnsleym1_Iteration(f, zPixel);
		break;
	case 11:
		r = Magnet1j_Iteration(f, zPixel);
		break;
	case 12:
		r = Magnet1m_Iteration(f, zPixel);
		break;
	case 13:
		r = BurningShip_Iteration(f, zPixel);
		break;
	case 14:
		r = Tricorn_Iteration(f, zPixel);
		break;
	case 15:
		r = Mandelbulb_Iteration(f, zPixel);
		break;
	case 18:
		r = PerpendicularBurningShip_Iteration(f, zPixel);
		break;
	case 19:
		r = Celtic_Iteration(f, zPixel);
		break;
	case 20:
		r = AlphaMandelbrot_Iteration(f, zPixel);
		break;
	case 21:
		r = PickoverStalks_Iteration(f, zPixel);
		break;
	case 22:
		r = Nova_Iteration(f, zPixel);
		break;
	case 23:
		r = Multibrot_Iteration(f, zPixel);
		break;
	default:
		r = Mendelbrot_Iteration(f, zPixel);
		break;
	}

	return r;
}
