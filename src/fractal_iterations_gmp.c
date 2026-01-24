/* fractal_iterations_gmp.c
   Fonctions d'iteration GMP (precision arbitraire)
   released under GPL2
   Copyleft 2001-2026 VERHILLE Arnaud
*/

#include "config.h"

#ifdef HAVE_GMP

#include "fractal_iterations_gmp.h"
#include <math.h>

fractalresult Mendelbrot_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp;
	mpf_t mag2, bailout_mpf, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);

	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zTemp, prec);
	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);

	do {
		i++;
		complex_gmp_sq_to(&zTemp, z, &f.mul_temps);
		complex_gmp_add_to(&z, zTemp, zPixel);
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult Julia_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, seed_gmp;
	mpf_t mag2, bailout_mpf, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(zPixel, prec);
	complex_gmp_init(&zTemp, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);

	do {
		i++;
		complex_gmp_sq_to(&zTemp, z, &f.mul_temps);
		complex_gmp_add_to(&z, zTemp, seed_gmp);
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult JuliaSin_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, seed_gmp, sin_z;
	mpf_t mag, bailout_mpf;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(zPixel, prec);

	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);

	do {
		i++;
		sin_z = complex_gmp_sin(z, prec);
		complex_gmp old_z = z;
		z = complex_gmp_mul(seed_gmp, sin_z, prec);
		complex_gmp_clear(&old_z);
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&sin_z);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);

	return result;
}

fractalresult Newton_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zNum, zQuot, p_minus_one, p_val, one, p_minus_one_z;
	mpf_t mag, bailout_mpf, p_mpf, p_minus_one_mpf, zero_mpf;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	z = complex_gmp_copy(zPixel, prec);

	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(p_mpf, prec);
	mpf_init2(p_minus_one_mpf, prec);
	mpf_init2(zero_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_ui(zero_mpf, 0);

	double p_d = Rez(f.seed);
	mpf_set_d(p_mpf, p_d);
	mpf_set_d(p_minus_one_mpf, p_d - 1.0);

	one = complex_gmp_one(prec);
	p_val = complex_gmp_make(p_mpf, zero_mpf, prec);
	p_minus_one = complex_gmp_make(p_minus_one_mpf, zero_mpf, prec);

	do {
		complex_gmp z_pow_p, z_pow_p_minus_one;
		i++;

		z_pow_p = complex_gmp_pow_n(z, p_val, 0, prec);
		p_minus_one_z = complex_gmp_mul(p_minus_one, z_pow_p, prec);
		zNum = complex_gmp_add(p_minus_one_z, one, prec);

		z_pow_p_minus_one = complex_gmp_pow_n(z, p_minus_one, 0, prec);
		zQuot = complex_gmp_mul(p_val, z_pow_p_minus_one, prec);

		z = complex_gmp_div(zNum, zQuot, prec);

		complex_gmp_mag(mag, z);

		complex_gmp_clear(&z_pow_p);
		complex_gmp_clear(&z_pow_p_minus_one);
		complex_gmp_clear(&p_minus_one_z);
		complex_gmp_clear(&zNum);
		complex_gmp_clear(&zQuot);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&one);
	complex_gmp_clear(&p_val);
	complex_gmp_clear(&p_minus_one);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(p_mpf);
	mpf_clear(p_minus_one_mpf);
	mpf_clear(zero_mpf);

	return result;
}

fractalresult Phoenix_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, y, zTemp, zpTemp;
	mpf_t mag, bailout_mpf, p1, p2, z_real;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	z = complex_gmp_copy(zPixel, prec);
	y = complex_gmp_zero(prec);
	zTemp = complex_gmp_zero(prec);

	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(p1, prec);
	mpf_init2(p2, prec);
	mpf_init2(z_real, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_d(p1, 0.56667);
	mpf_set_d(p2, -0.5);

	do {
		i++;
		complex_gmp old_zTemp = zTemp;
		zTemp = complex_gmp_mul(z, z, prec);
		complex_gmp_clear(&old_zTemp);
		complex_gmp_get_real(z_real, zTemp);
		mpf_add(z_real, z_real, p1);
		complex_gmp_set_real(&zTemp, z_real);

		zpTemp = complex_gmp_scalar_mul(p2, y, prec);
		complex_gmp old_zTemp2 = zTemp;
		zTemp = complex_gmp_add(zTemp, zpTemp, prec);
		complex_gmp_clear(&old_zTemp2);

		complex_gmp old_y = y;
		complex_gmp old_z = z;
		y = z;
		z = zTemp;
		complex_gmp_clear(&old_y);
		complex_gmp_clear(&old_z);

		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zpTemp);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&y);
	complex_gmp_clear(&zTemp);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(p1);
	mpf_clear(p2);
	mpf_clear(z_real);

	return result;
}

fractalresult Barnsleyj1_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, seed_gmp, one;
	mpf_t mag, bailout_mpf, z_real, one_mpf;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(zPixel, prec);
	one = complex_gmp_one(prec);

	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(z_real, prec);
	mpf_init2(one_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_ui(one_mpf, 1);

	do {
		i++;
		complex_gmp_get_real(z_real, z);

		if (mpf_cmp(z_real, one_mpf) >= 0) {
			zTemp = complex_gmp_sub(z, one, prec);
		} else {
			zTemp = complex_gmp_add(z, one, prec);
		}
		complex_gmp old_z = z;
		z = complex_gmp_mul(zTemp, seed_gmp, prec);
		complex_gmp_clear(&old_z);
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zTemp);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	complex_gmp_clear(&one);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(z_real);
	mpf_clear(one_mpf);

	return result;
}

fractalresult Barnsley1m_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, c, zTemp, one;
	mpf_t mag, bailout_mpf, z_real, one_mpf;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	c = complex_gmp_copy(zPixel, prec);
	z = complex_gmp_zero(prec);
	one = complex_gmp_one(prec);

	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(z_real, prec);
	mpf_init2(one_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_ui(one_mpf, 1);

	do {
		i++;
		complex_gmp_get_real(z_real, z);

		if (mpf_cmp(z_real, one_mpf) >= 0) {
			zTemp = complex_gmp_sub(z, one, prec);
		} else {
			zTemp = complex_gmp_add(z, one, prec);
		}
		complex_gmp old_z = z;
		z = complex_gmp_mul(zTemp, c, prec);
		complex_gmp_clear(&old_z);
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zTemp);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&c);
	complex_gmp_clear(&one);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(z_real);
	mpf_clear(one_mpf);

	return result;
}

fractalresult Magnet1j_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, N, Q, seed_gmp, one, two, seed_minus_one, seed_minus_two;
	mpf_t mag, bailout_mpf;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(zPixel, prec);

	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);

	mpf_t two_mpf;
	mpf_init2(two_mpf, prec);
	mpf_set_ui(two_mpf, 2);
	one = complex_gmp_one(prec);
	two = complex_gmp_scalar_mul(two_mpf, one, prec);
	mpf_clear(two_mpf);

	seed_minus_one = complex_gmp_sub(seed_gmp, one, prec);
	seed_minus_two = complex_gmp_sub(seed_gmp, two, prec);

	do {
		i++;
		complex_gmp z_sq = complex_gmp_mul(z, z, prec);
		complex_gmp N_temp = complex_gmp_add(z_sq, seed_minus_one, prec);
		N = complex_gmp_mul(N_temp, N_temp, prec);
		complex_gmp_clear(&N_temp);

		complex_gmp two_z = complex_gmp_mul(two, z, prec);
		Q = complex_gmp_add(two_z, seed_minus_two, prec);

		complex_gmp old_z = z;
		z = complex_gmp_div(N, Q, prec);
		complex_gmp_clear(&old_z);

		complex_gmp_mag(mag, z);

		complex_gmp_clear(&z_sq);
		complex_gmp_clear(&two_z);
		complex_gmp_clear(&N);
		complex_gmp_clear(&Q);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	complex_gmp_clear(&one);
	complex_gmp_clear(&two);
	complex_gmp_clear(&seed_minus_one);
	complex_gmp_clear(&seed_minus_two);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);

	return result;
}

fractalresult Magnet1m_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, c, N, Q, one, two, c_minus_one, c_minus_two;
	mpf_t mag, bailout_mpf;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	c = complex_gmp_copy(zPixel, prec);
	z = complex_gmp_zero(prec);

	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);

	mpf_t two_mpf;
	mpf_init2(two_mpf, prec);
	mpf_set_ui(two_mpf, 2);
	one = complex_gmp_one(prec);
	two = complex_gmp_scalar_mul(two_mpf, one, prec);
	mpf_clear(two_mpf);

	c_minus_one = complex_gmp_sub(c, one, prec);
	c_minus_two = complex_gmp_sub(c, two, prec);

	do {
		i++;
		complex_gmp z_sq = complex_gmp_mul(z, z, prec);
		complex_gmp N_temp = complex_gmp_add(z_sq, c_minus_one, prec);
		N = complex_gmp_mul(N_temp, N_temp, prec);
		complex_gmp_clear(&N_temp);

		complex_gmp two_z = complex_gmp_mul(two, z, prec);
		Q = complex_gmp_add(two_z, c_minus_two, prec);

		complex_gmp old_z = z;
		z = complex_gmp_div(N, Q, prec);
		complex_gmp_clear(&old_z);

		complex_gmp_mag(mag, z);

		complex_gmp_clear(&z_sq);
		complex_gmp_clear(&two_z);
		complex_gmp_clear(&N);
		complex_gmp_clear(&Q);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&c);
	complex_gmp_clear(&one);
	complex_gmp_clear(&two);
	complex_gmp_clear(&c_minus_one);
	complex_gmp_clear(&c_minus_two);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);

	return result;
}

fractalresult BurningShip_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, zAbs;
	mpf_t mag2, bailout_mpf, abs_real, abs_imag, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zTemp, prec);
	complex_gmp_init(&zAbs, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(abs_real, prec);
	mpf_init2(abs_imag, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);

	do {
		i++;
		mpf_abs(abs_real, z.x);
		mpf_abs(abs_imag, z.y);

		mpf_set(zAbs.x, abs_real);
		mpf_set(zAbs.y, abs_imag);

		complex_gmp_sq_to(&zTemp, zAbs, &f.mul_temps);
		complex_gmp_add_to(&z, zTemp, zPixel);
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&zAbs);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(abs_real);
	mpf_clear(abs_imag);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult Buffalo_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zSq;
	mpf_t mag2, bailout_mpf, abs_real, abs_imag, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zSq, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(abs_real, prec);
	mpf_init2(abs_imag, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);

	do {
		i++;
		complex_gmp_sq_to(&zSq, z, &f.mul_temps);
		mpf_abs(abs_real, zSq.x);
		mpf_abs(abs_imag, zSq.y);
		mpf_add(z.x, abs_real, zPixel.x);
		mpf_add(z.y, abs_imag, zPixel.y);
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zSq);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(abs_real);
	mpf_clear(abs_imag);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult Tricorn_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, zConj;
	mpf_t mag2, bailout_mpf, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zTemp, prec);
	complex_gmp_init(&zConj, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);

	do {
		i++;
		mpf_set(zConj.x, z.x);
		mpf_neg(zConj.y, z.y);

		complex_gmp_sq_to(&zTemp, zConj, &f.mul_temps);
		complex_gmp_add_to(&z, zTemp, zPixel);
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&zConj);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult Mandelbulb_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, z2, z4, z8;
	mpf_t mag2, bailout_mpf, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&z2, prec);
	complex_gmp_init(&z4, prec);
	complex_gmp_init(&z8, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);

	do {
		i++;
		complex_gmp_sq_to(&z2, z, &f.mul_temps);
		complex_gmp_sq_to(&z4, z2, &f.mul_temps);
		complex_gmp_sq_to(&z8, z4, &f.mul_temps);
		complex_gmp_add_to(&z, z8, zPixel);
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&z2);
	complex_gmp_clear(&z4);
	complex_gmp_clear(&z8);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult PerpendicularBurningShip_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp;
	mpf_t mag2, bailout_mpf, temp1, temp2, x, y, y_abs, x2, y2, x_new, y_new;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zTemp, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_init2(x, prec);
	mpf_init2(y, prec);
	mpf_init2(y_abs, prec);
	mpf_init2(x2, prec);
	mpf_init2(y2, prec);
	mpf_init2(x_new, prec);
	mpf_init2(y_new, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);

	do {
		i++;
		mpf_set(x, z.x);
		mpf_set(y, z.y);
		mpf_abs(y_abs, y);

		mpf_mul(x2, x, x);
		mpf_mul(y2, y_abs, y_abs);
		mpf_sub(x_new, x2, y2);
		mpf_add(x_new, x_new, zPixel.x);

		mpf_mul(y_new, x, y_abs);
		mpf_mul_ui(y_new, y_new, 2);
		mpf_neg(y_new, y_new);
		mpf_add(y_new, y_new, zPixel.y);

		mpf_set(z.x, x_new);
		mpf_set(z.y, y_new);

		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);
	mpf_clear(x);
	mpf_clear(y);
	mpf_clear(y_abs);
	mpf_clear(x2);
	mpf_clear(y2);
	mpf_clear(x_new);
	mpf_clear(y_new);

	return result;
}

fractalresult Celtic_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zSq;
	mpf_t mag2, bailout_mpf, temp1, temp2, u, u_abs;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zSq, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_init2(u, prec);
	mpf_init2(u_abs, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);

	do {
		i++;
		complex_gmp_sq_to(&zSq, z, &f.mul_temps);

		mpf_abs(u_abs, zSq.x);
		mpf_add(z.x, u_abs, zPixel.x);
		mpf_add(z.y, zSq.y, zPixel.y);

		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zSq);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);
	mpf_clear(u);
	mpf_clear(u_abs);

	return result;
}

fractalresult AlphaMandelbrot_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zSq, m, mSq, zTemp;
	mpf_t mag2, bailout_mpf, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zSq, prec);
	complex_gmp_init(&m, prec);
	complex_gmp_init(&mSq, prec);
	complex_gmp_init(&zTemp, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);

	do {
		i++;
		complex_gmp_sq_to(&zSq, z, &f.mul_temps);
		complex_gmp_add_to(&m, zSq, zPixel);
		complex_gmp_sq_to(&mSq, m, &f.mul_temps);
		complex_gmp_add_to(&zTemp, zSq, mSq);
		complex_gmp_add_to(&z, zTemp, zPixel);
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zSq);
	complex_gmp_clear(&m);
	complex_gmp_clear(&mSq);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult PickoverStalks_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp;
	mpf_t mag2, bailout_mpf, temp1, temp2, re_abs, im_abs, trap_distance, trap_min, trap_divisor;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zTemp, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_init2(re_abs, prec);
	mpf_init2(im_abs, prec);
	mpf_init2(trap_distance, prec);
	mpf_init2(trap_min, prec);
	mpf_init2(trap_divisor, prec);

	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);
	mpf_set_d(trap_divisor, 0.03);
	mpf_set_d(trap_min, 1e10);

	do {
		i++;
		complex_gmp_sq_to(&zTemp, z, &f.mul_temps);
		complex_gmp_add_to(&z, zTemp, zPixel);

		mpf_abs(re_abs, z.x);
		mpf_abs(im_abs, z.y);
		if (mpf_cmp(re_abs, im_abs) < 0) {
			mpf_set(trap_distance, re_abs);
		} else {
			mpf_set(trap_distance, im_abs);
		}

		if (mpf_cmp(trap_distance, trap_min) < 0) {
			mpf_set(trap_min, trap_distance);
		}

		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	mpf_t threshold, ratio;
	mpf_init2(threshold, prec);
	mpf_init2(ratio, prec);
	mpf_set_d(threshold, 1e-10);

	if (mpf_cmp(trap_min, threshold) > 0) {
		mpf_div(ratio, trap_min, trap_divisor);
		double ratio_d = mpf_get_d(ratio);
		double log_val = -log(ratio_d) * 100.0;
		result.iteration = (int)log_val;
		if (result.iteration < 0) result.iteration = 0;
		if (result.iteration >= f.iterationMax) result.iteration = f.iterationMax - 1;
	} else {
		result.iteration = f.iterationMax - 1;
	}

	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);
	mpf_clear(re_abs);
	mpf_clear(im_abs);
	mpf_clear(trap_distance);
	mpf_clear(trap_min);
	mpf_clear(trap_divisor);
	mpf_clear(threshold);
	mpf_clear(ratio);

	return result;
}

fractalresult Nova_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, z_prev, z_pow, z_pow_deriv, numerator, denominator, newton_step, a_relax, one, p_val;
	mpf_t mag, bailout_mpf, conv_epsilon, conv_epsilon_sq, diff_sq, z_sq, denom, p_mpf, zero_mpf, a_real;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;

	int p = 3;
	double conv_epsilon_d = 1e-7;

	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(conv_epsilon, prec);
	mpf_init2(conv_epsilon_sq, prec);
	mpf_init2(diff_sq, prec);
	mpf_init2(z_sq, prec);
	mpf_init2(denom, prec);
	mpf_init2(p_mpf, prec);
	mpf_init2(zero_mpf, prec);
	mpf_init2(a_real, prec);

	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_d(conv_epsilon, conv_epsilon_d);
	mpf_mul(conv_epsilon_sq, conv_epsilon, conv_epsilon);
	mpf_set_ui(zero_mpf, 0);
	mpf_set_d(p_mpf, (double)p);
	mpf_set_d(a_real, 1.0);

	one = complex_gmp_one(prec);
	z = complex_gmp_copy(one, prec);
	z_prev = complex_gmp_copy(z, prec);
	a_relax = complex_gmp_make(a_real, zero_mpf, prec);
	p_val = complex_gmp_make(p_mpf, zero_mpf, prec);

	complex_gmp_init(&z_pow, prec);
	complex_gmp_init(&z_pow_deriv, prec);
	complex_gmp_init(&numerator, prec);
	complex_gmp_init(&denominator, prec);
	complex_gmp_init(&newton_step, prec);

	do {
		mpf_t p_minus_one_mpf;
		mpf_init2(p_minus_one_mpf, prec);
		mpf_sub_ui(p_minus_one_mpf, p_mpf, 1);
		complex_gmp p_minus_one_val = complex_gmp_make(p_minus_one_mpf, zero_mpf, prec);

		i++;

		complex_gmp_clear(&z_pow);
		complex_gmp_clear(&z_pow_deriv);
		complex_gmp_clear(&numerator);
		complex_gmp_clear(&denominator);
		complex_gmp_clear(&newton_step);

		z_pow = complex_gmp_pow_n(z, p_val, 0, prec);
		z_pow_deriv = complex_gmp_pow_n(z, p_minus_one_val, 0, prec);

		numerator = complex_gmp_sub(z_pow, one, prec);
		denominator = complex_gmp_mul(p_val, z_pow_deriv, prec);

		complex_gmp_mag(mag, denominator);
		if (mpf_cmp_d(mag, 1e-10) < 0) {
			complex_gmp_clear(&p_minus_one_val);
			mpf_clear(p_minus_one_mpf);
			break;
		}

		newton_step = complex_gmp_div(numerator, denominator, prec);
		complex_gmp old_newton = newton_step;
		newton_step = complex_gmp_mul(a_relax, newton_step, prec);
		complex_gmp_clear(&old_newton);

		complex_gmp_clear(&z_prev);
		z_prev = complex_gmp_copy(z, prec);
		complex_gmp old_z = z;
		z = complex_gmp_sub(z, newton_step, prec);
		complex_gmp_clear(&old_z);
		complex_gmp old_z2 = z;
		z = complex_gmp_add(z, zPixel, prec);
		complex_gmp_clear(&old_z2);

		complex_gmp diff = complex_gmp_sub(z, z_prev, prec);
		complex_gmp_mag2(diff_sq, diff);
		complex_gmp_mag2(z_sq, z);
		if (mpf_cmp_d(z_sq, 1.0) < 0) {
			mpf_set_ui(denom, 1);
		} else {
			mpf_set(denom, z_sq);
		}

		mpf_div(diff_sq, diff_sq, denom);
		if (mpf_cmp(diff_sq, conv_epsilon_sq) < 0) {
			complex_gmp_clear(&diff);
			complex_gmp_clear(&p_minus_one_val);
			mpf_clear(p_minus_one_mpf);
			break;
		}
		complex_gmp_clear(&diff);

		complex_gmp_mag(mag, z);
		if (mpf_cmp(mag, bailout_mpf) > 0) {
			complex_gmp_clear(&p_minus_one_val);
			mpf_clear(p_minus_one_mpf);
			break;
		}

		complex_gmp_clear(&p_minus_one_val);
		mpf_clear(p_minus_one_mpf);
	} while (i < f.iterationMax);

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&z_prev);
	complex_gmp_clear(&z_pow);
	complex_gmp_clear(&z_pow_deriv);
	complex_gmp_clear(&numerator);
	complex_gmp_clear(&denominator);
	complex_gmp_clear(&newton_step);
	complex_gmp_clear(&a_relax);
	complex_gmp_clear(&one);
	complex_gmp_clear(&p_val);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(conv_epsilon);
	mpf_clear(conv_epsilon_sq);
	mpf_clear(diff_sq);
	mpf_clear(z_sq);
	mpf_clear(denom);
	mpf_clear(p_mpf);
	mpf_clear(zero_mpf);
	mpf_clear(a_real);

	return result;
}

fractalresult Multibrot_Iteration_GMP(fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, z_pow, d_val;
	mpf_t mag2, bailout_mpf, temp1, temp2, d_mpf, zero_mpf;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	double d = 2.5;

	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&z_pow, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_init2(d_mpf, prec);
	mpf_init2(zero_mpf, prec);

	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);
	mpf_set_d(d_mpf, d);
	mpf_set_ui(zero_mpf, 0);
	d_val = complex_gmp_make(d_mpf, zero_mpf, prec);

	do {
		i++;

		z_pow = complex_gmp_pow(z, d_val, prec);

		if (complex_gmp_is_nan(z_pow) || complex_gmp_is_inf(z_pow)) {
			break;
		}

		complex_gmp old_z = z;
		z = complex_gmp_add(z_pow, zPixel, prec);
		complex_gmp_clear(&old_z);
		complex_gmp_clear(&z_pow);

		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&z_pow);
	complex_gmp_clear(&seed_gmp);
	complex_gmp_clear(&d_val);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);
	mpf_clear(d_mpf);
	mpf_clear(zero_mpf);

	return result;
}

fractalresult FormulaSelector_GMP(fractal f, complex_gmp zPixel) {
	fractalresult r;

	switch (f.type) {
	case 3:
		r = Mendelbrot_Iteration_GMP(f, zPixel);
		break;
	case 4:
		r = Julia_Iteration_GMP(f, zPixel);
		break;
	case 5:
		r = JuliaSin_Iteration_GMP(f, zPixel);
		break;
	case 6:
		r = Newton_Iteration_GMP(f, zPixel);
		break;
	case 7:
		r = Phoenix_Iteration_GMP(f, zPixel);
		break;
	case 8:
		r = Buffalo_Iteration_GMP(f, zPixel);
		break;
	case 9:
		r = Barnsleyj1_Iteration_GMP(f, zPixel);
		break;
	case 10:
		r = Barnsley1m_Iteration_GMP(f, zPixel);
		break;
	case 11:
		r = Magnet1j_Iteration_GMP(f, zPixel);
		break;
	case 12:
		r = Magnet1m_Iteration_GMP(f, zPixel);
		break;
	case 13:
		r = BurningShip_Iteration_GMP(f, zPixel);
		break;
	case 14:
		r = Tricorn_Iteration_GMP(f, zPixel);
		break;
	case 15:
		r = Mandelbulb_Iteration_GMP(f, zPixel);
		break;
	case 18:
		r = PerpendicularBurningShip_Iteration_GMP(f, zPixel);
		break;
	case 19:
		r = Celtic_Iteration_GMP(f, zPixel);
		break;
	case 20:
		r = AlphaMandelbrot_Iteration_GMP(f, zPixel);
		break;
	case 21:
		r = PickoverStalks_Iteration_GMP(f, zPixel);
		break;
	case 22:
		r = Nova_Iteration_GMP(f, zPixel);
		break;
	case 23:
		r = Multibrot_Iteration_GMP(f, zPixel);
		break;
	default:
		r = Mendelbrot_Iteration_GMP(f, zPixel);
		break;
	}

	return r;
}

#endif /* HAVE_GMP */
