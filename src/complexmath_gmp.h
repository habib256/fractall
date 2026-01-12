/*
 * Complexmath_gmp.h
 * Header file for Complex Number operations using GMP (GNU Multiple Precision)
 * Released under GPL2
 * Copyleft 2024
 */

#include <gmp.h>
#include "complexmath.h"

#ifndef COMPLEXGMPFLAG
#define COMPLEXGMPFLAG

// Structure pour nombre complexe avec GMP
typedef struct {
    mpf_t x;  // Partie réelle
    mpf_t y;  // Partie imaginaire
    int is_nan;  // Flag NaN (1 si NaN, 0 sinon)
} complex_gmp;

// Seuil pour détection d'infinité (2^GMP_INF_EXP)
#define GMP_INF_EXP 1000

// Initialisation et destruction
void complex_gmp_init(complex_gmp* z, mp_bitcnt_t prec);
void complex_gmp_clear(complex_gmp* z);
void complex_gmp_set_prec(complex_gmp* z, mp_bitcnt_t prec);

// Constructeurs
complex_gmp complex_gmp_make(mpf_t real, mpf_t imag, mp_bitcnt_t prec);
complex_gmp complex_gmp_zero(mp_bitcnt_t prec);
complex_gmp complex_gmp_one(mp_bitcnt_t prec);
complex_gmp complex_gmp_i(mp_bitcnt_t prec);
complex_gmp complex_gmp_inf(mp_bitcnt_t prec);
complex_gmp complex_gmp_nan(mp_bitcnt_t prec);
complex_gmp complex_gmp_copy(complex_gmp z, mp_bitcnt_t prec);

// Conversion
complex_gmp complex_to_gmp(complex z, mp_bitcnt_t prec);
complex gmp_to_complex(complex_gmp zg);

// Accès aux parties
void complex_gmp_get_real(mpf_t result, complex_gmp z);
void complex_gmp_get_imag(mpf_t result, complex_gmp z);
void complex_gmp_set_real(complex_gmp* z, mpf_t real);
void complex_gmp_set_imag(complex_gmp* z, mpf_t imag);

// Opérations arithmétiques
complex_gmp complex_gmp_add(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec);
complex_gmp complex_gmp_sub(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec);
complex_gmp complex_gmp_mul(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec);
complex_gmp complex_gmp_div(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec);
complex_gmp complex_gmp_neg(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_scalar_mul(mpf_t c, complex_gmp z, mp_bitcnt_t prec);

// Fonctions mathématiques
void complex_gmp_mag(mpf_t result, complex_gmp z);
void complex_gmp_mag2(mpf_t result, complex_gmp z);
void complex_gmp_arg(mpf_t result, complex_gmp z);

complex_gmp complex_gmp_sqrt(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_pow(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec);
complex_gmp complex_gmp_pow_n(complex_gmp z1, complex_gmp z2, int n, mp_bitcnt_t prec);
complex_gmp complex_gmp_sin(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_cos(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_tan(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_exp(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_sinh(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_cosh(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_tanh(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_log(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_arcsin(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_arccos(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_arctan(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_arcsinh(complex_gmp z, mp_bitcnt_t prec);
complex_gmp complex_gmp_arccosh(complex_gmp z, mp_bitcnt_t prec);

// Utilitaires
int complex_gmp_is_zero(complex_gmp z);
int complex_gmp_is_inf(complex_gmp z);
int complex_gmp_is_nan(complex_gmp z);
int complex_gmp_cmp_mag(mpf_t threshold, complex_gmp z);

// Structure pour temporaires de multiplication (évite allocations répétées)
typedef struct {
    mpf_t t1, t2, t3, t4;
    int initialized;
} gmp_mul_temps;

// Initialisation/destruction des temporaires
void gmp_mul_temps_init(gmp_mul_temps* temps, mp_bitcnt_t prec);
void gmp_mul_temps_clear(gmp_mul_temps* temps);

// Opérations in-place (pas d'allocation, écrit directement dans result)
void complex_gmp_add_to(complex_gmp* result, complex_gmp z1, complex_gmp z2);
void complex_gmp_sub_to(complex_gmp* result, complex_gmp z1, complex_gmp z2);
void complex_gmp_mul_to(complex_gmp* result, complex_gmp z1, complex_gmp z2, gmp_mul_temps* temps);
void complex_gmp_sq_to(complex_gmp* result, complex_gmp z, gmp_mul_temps* temps);
void complex_gmp_mag2_to(mpf_t result, complex_gmp z, mpf_t temp1, mpf_t temp2);

// Copie in-place (sans allocation)
void complex_gmp_copy_to(complex_gmp* dest, complex_gmp src);

#endif /* COMPLEXGMPFLAG */
