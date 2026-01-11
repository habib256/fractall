/*
 * Complexmath_gmp.c
 * Complex functions using GMP (GNU Multiple Precision)
 * Released under GPL2
 * Copyleft 2024
 */

#include <math.h>
#include <gmp.h>
#include "complexmath_gmp.h"
#include "exception.h"

// Initialisation et destruction
void complex_gmp_init(complex_gmp* z, mp_bitcnt_t prec) {
    mpf_init2(z->x, prec);
    mpf_init2(z->y, prec);
}

void complex_gmp_clear(complex_gmp* z) {
    mpf_clear(z->x);
    mpf_clear(z->y);
}

void complex_gmp_set_prec(complex_gmp* z, mp_bitcnt_t prec) {
    mpf_set_prec(z->x, prec);
    mpf_set_prec(z->y, prec);
}

// Constructeurs
complex_gmp complex_gmp_make(mpf_t real, mpf_t imag, mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_set(rz.x, real);
    mpf_set(rz.y, imag);
    return rz;
}

complex_gmp complex_gmp_zero(mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_set_ui(rz.x, 0);
    mpf_set_ui(rz.y, 0);
    return rz;
}

complex_gmp complex_gmp_one(mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_set_ui(rz.x, 1);
    mpf_set_ui(rz.y, 0);
    return rz;
}

complex_gmp complex_gmp_i(mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_set_ui(rz.x, 0);
    mpf_set_ui(rz.y, 1);
    return rz;
}

complex_gmp complex_gmp_inf(mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    // GMP n'a pas d'infini direct, on utilise un très grand nombre
    mpf_set_ui(rz.x, 1);
    mpf_mul_2exp(rz.x, rz.x, 1000); // 2^1000
    mpf_set_ui(rz.y, 0);
    return rz;
}

complex_gmp complex_gmp_nan(mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    // NaN simulé avec 0/0
    mpf_set_ui(rz.x, 0);
    mpf_set_ui(rz.y, 0);
    return rz;
}

complex_gmp complex_gmp_copy(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_set(rz.x, z.x);
    mpf_set(rz.y, z.y);
    return rz;
}

// Conversion
complex_gmp complex_to_gmp(complex z, mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_set_d(rz.x, z.x);
    mpf_set_d(rz.y, z.y);
    return rz;
}

complex gmp_to_complex(complex_gmp zg) {
    complex rz;
    rz.x = mpf_get_d(zg.x);
    rz.y = mpf_get_d(zg.y);
    return rz;
}

// Accès aux parties
void complex_gmp_get_real(mpf_t result, complex_gmp z) {
    mpf_set(result, z.x);
}

void complex_gmp_get_imag(mpf_t result, complex_gmp z) {
    mpf_set(result, z.y);
}

void complex_gmp_set_real(complex_gmp* z, mpf_t real) {
    mpf_set(z->x, real);
}

void complex_gmp_set_imag(complex_gmp* z, mpf_t imag) {
    mpf_set(z->y, imag);
}

// Opérations arithmétiques
complex_gmp complex_gmp_add(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_add(rz.x, z1.x, z2.x);
    mpf_add(rz.y, z1.y, z2.y);
    return rz;
}

complex_gmp complex_gmp_sub(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_sub(rz.x, z1.x, z2.x);
    mpf_sub(rz.y, z1.y, z2.y);
    return rz;
}

complex_gmp complex_gmp_mul(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t temp1, temp2, temp3, temp4;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(temp1, prec);
    mpf_init2(temp2, prec);
    mpf_init2(temp3, prec);
    mpf_init2(temp4, prec);
    
    // rz.x = z1.x * z2.x - z1.y * z2.y
    mpf_mul(temp1, z1.x, z2.x);
    mpf_mul(temp2, z1.y, z2.y);
    mpf_sub(rz.x, temp1, temp2);
    
    // rz.y = z1.x * z2.y + z1.y * z2.x
    mpf_mul(temp3, z1.x, z2.y);
    mpf_mul(temp4, z1.y, z2.x);
    mpf_add(rz.y, temp3, temp4);
    
    mpf_clear(temp1);
    mpf_clear(temp2);
    mpf_clear(temp3);
    mpf_clear(temp4);
    
    return rz;
}

complex_gmp complex_gmp_div(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t denom, temp1, temp2, temp3, temp4;
    int z2_is_zero;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(denom, prec);
    mpf_init2(temp1, prec);
    mpf_init2(temp2, prec);
    mpf_init2(temp3, prec);
    mpf_init2(temp4, prec);
    
    // Calcul du dénominateur: |z2|^2 = z2.x^2 + z2.y^2
    mpf_mul(temp1, z2.x, z2.x);
    mpf_mul(temp2, z2.y, z2.y);
    mpf_add(denom, temp1, temp2);
    
    z2_is_zero = (mpf_cmp_ui(denom, 0) == 0);
    
    if (z2_is_zero) {
        // Division par zéro
        complex_gmp_clear(&rz);
        rz = complex_gmp_nan(prec);
        mpf_clear(denom);
        mpf_clear(temp1);
        mpf_clear(temp2);
        mpf_clear(temp3);
        mpf_clear(temp4);
        THROW(EX_DIVZ);
        return rz;
    }
    
    // rz.x = (z1.x * z2.x + z1.y * z2.y) / denom
    mpf_mul(temp1, z1.x, z2.x);
    mpf_mul(temp2, z1.y, z2.y);
    mpf_add(temp3, temp1, temp2);
    mpf_div(rz.x, temp3, denom);
    
    // rz.y = (z1.y * z2.x - z1.x * z2.y) / denom
    mpf_mul(temp1, z1.y, z2.x);
    mpf_mul(temp2, z1.x, z2.y);
    mpf_sub(temp3, temp1, temp2);
    mpf_div(rz.y, temp3, denom);
    
    mpf_clear(denom);
    mpf_clear(temp1);
    mpf_clear(temp2);
    mpf_clear(temp3);
    mpf_clear(temp4);
    
    return rz;
}

complex_gmp complex_gmp_neg(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_neg(rz.x, z.x);
    mpf_neg(rz.y, z.y);
    return rz;
}

complex_gmp complex_gmp_scalar_mul(mpf_t c, complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    complex_gmp_init(&rz, prec);
    mpf_mul(rz.x, c, z.x);
    mpf_mul(rz.y, c, z.y);
    return rz;
}

// Fonctions mathématiques
void complex_gmp_mag(mpf_t result, complex_gmp z) {
    mpf_t temp1, temp2, sum;
    mp_bitcnt_t prec = mpf_get_prec(z.x);
    
    mpf_init2(temp1, prec);
    mpf_init2(temp2, prec);
    mpf_init2(sum, prec);
    
    mpf_mul(temp1, z.x, z.x);
    mpf_mul(temp2, z.y, z.y);
    mpf_add(sum, temp1, temp2);
    mpf_sqrt(result, sum);
    
    mpf_clear(temp1);
    mpf_clear(temp2);
    mpf_clear(sum);
}

void complex_gmp_mag2(mpf_t result, complex_gmp z) {
    mpf_t temp1, temp2;
    mp_bitcnt_t prec = mpf_get_prec(z.x);
    
    mpf_init2(temp1, prec);
    mpf_init2(temp2, prec);
    
    mpf_mul(temp1, z.x, z.x);
    mpf_mul(temp2, z.y, z.y);
    mpf_add(result, temp1, temp2);
    
    mpf_clear(temp1);
    mpf_clear(temp2);
}

void complex_gmp_arg(mpf_t result, complex_gmp z) {
    double x, y;
    x = mpf_get_d(z.x);
    y = mpf_get_d(z.y);
    mpf_set_d(result, atan2(y, x));
}

complex_gmp complex_gmp_sqrt(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t mag, arg, half_arg, cos_val, sin_val, sqrt_mag;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(mag, prec);
    mpf_init2(arg, prec);
    mpf_init2(half_arg, prec);
    mpf_init2(cos_val, prec);
    mpf_init2(sin_val, prec);
    mpf_init2(sqrt_mag, prec);
    
    complex_gmp_mag(mag, z);
    complex_gmp_arg(arg, z);
    
    // half_arg = arg / 2
    mpf_div_ui(half_arg, arg, 2);
    
    // sqrt_mag = sqrt(mag)
    mpf_sqrt(sqrt_mag, mag);
    
    // cos et sin de half_arg
    double half_arg_d = mpf_get_d(half_arg);
    mpf_set_d(cos_val, cos(half_arg_d));
    mpf_set_d(sin_val, sin(half_arg_d));
    
    // rz.x = sqrt_mag * cos(half_arg)
    mpf_mul(rz.x, sqrt_mag, cos_val);
    
    // rz.y = sqrt_mag * sin(half_arg)
    mpf_mul(rz.y, sqrt_mag, sin_val);
    
    mpf_clear(mag);
    mpf_clear(arg);
    mpf_clear(half_arg);
    mpf_clear(cos_val);
    mpf_clear(sin_val);
    mpf_clear(sqrt_mag);
    
    return rz;
}

complex_gmp complex_gmp_pow_n(complex_gmp z1, complex_gmp z2, int n, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t mag, arg, log_mag, exp_part, cos_part, sin_part, temp;
    double arg_d, exp_d, cos_d, sin_d;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(mag, prec);
    mpf_init2(arg, prec);
    mpf_init2(log_mag, prec);
    mpf_init2(exp_part, prec);
    mpf_init2(cos_part, prec);
    mpf_init2(sin_part, prec);
    mpf_init2(temp, prec);
    
    complex_gmp_mag(mag, z1);
    complex_gmp_arg(arg, z1);
    
    // Si mag == 0
    if (mpf_cmp_ui(mag, 0) == 0) {
        if (mpf_cmp_ui(z2.x, 0) == 0 && mpf_cmp_ui(z2.y, 0) == 0) {
            complex_gmp_clear(&rz);
            rz = complex_gmp_nan(prec);
            mpf_clear(mag);
            mpf_clear(arg);
            mpf_clear(log_mag);
            mpf_clear(exp_part);
            mpf_clear(cos_part);
            mpf_clear(sin_part);
            mpf_clear(temp);
            THROW(EX_POWZN);
            return rz;
        }
        if (mpf_cmp_ui(z2.y, 0) != 0) {
            complex_gmp_clear(&rz);
            rz = complex_gmp_nan(prec);
            mpf_clear(mag);
            mpf_clear(arg);
            mpf_clear(log_mag);
            mpf_clear(exp_part);
            mpf_clear(cos_part);
            mpf_clear(sin_part);
            mpf_clear(temp);
            THROW(EX_POWZN);
            return rz;
        }
        if (mpf_cmp_ui(z2.x, 0) != 0) {
            complex_gmp_clear(&rz);
            rz = complex_gmp_zero(prec);
            mpf_clear(mag);
            mpf_clear(arg);
            mpf_clear(log_mag);
            mpf_clear(exp_part);
            mpf_clear(cos_part);
            mpf_clear(sin_part);
            mpf_clear(temp);
            return rz;
        }
    }
    
    // Si z2 == 0
    if (mpf_cmp_ui(z2.x, 0) == 0 && mpf_cmp_ui(z2.y, 0) == 0) {
        complex_gmp_clear(&rz);
        rz = complex_gmp_one(prec);
        mpf_clear(mag);
        mpf_clear(arg);
        mpf_clear(log_mag);
        mpf_clear(exp_part);
        mpf_clear(cos_part);
        mpf_clear(sin_part);
        mpf_clear(temp);
        return rz;
    }
    
    // Calcul: z1^z2 = exp(z2 * log(z1))
    // log(z1) = log(mag) + i*arg
    double mag_d = mpf_get_d(mag);
    if (mag_d > 0) {
        mpf_set_d(log_mag, log(mag_d));
    } else {
        complex_gmp_clear(&rz);
        rz = complex_gmp_nan(prec);
        mpf_clear(mag);
        mpf_clear(arg);
        mpf_clear(log_mag);
        mpf_clear(exp_part);
        mpf_clear(cos_part);
        mpf_clear(sin_part);
        mpf_clear(temp);
        return rz;
    }
    
    // arg_total = arg + 2*n*PI
    arg_d = mpf_get_d(arg) + 2.0 * n * M_PI;
    
    // exp_part = mag^z2.x * exp(-z2.y * arg_total)
    double z2_x_d = mpf_get_d(z2.x);
    double z2_y_d = mpf_get_d(z2.y);
    exp_d = pow(mag_d, z2_x_d) * exp(-z2_y_d * arg_d);
    
    // th = z2.y * log(mag) + z2.x * arg_total
    double th = z2_y_d * log(mag_d) + z2_x_d * arg_d;
    
    cos_d = cos(th);
    sin_d = sin(th);
    
    mpf_set_d(rz.x, exp_d * cos_d);
    mpf_set_d(rz.y, exp_d * sin_d);
    
    mpf_clear(mag);
    mpf_clear(arg);
    mpf_clear(log_mag);
    mpf_clear(exp_part);
    mpf_clear(cos_part);
    mpf_clear(sin_part);
    mpf_clear(temp);
    
    return rz;
}

complex_gmp complex_gmp_pow(complex_gmp z1, complex_gmp z2, mp_bitcnt_t prec) {
    return complex_gmp_pow_n(z1, z2, 0, prec);
}

complex_gmp complex_gmp_sin(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t x, y, cosh_y, sinh_y, sin_x, cos_x;
    double x_d, y_d, cosh_y_d, sinh_y_d;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(x, prec);
    mpf_init2(y, prec);
    mpf_init2(cosh_y, prec);
    mpf_init2(sinh_y, prec);
    mpf_init2(sin_x, prec);
    mpf_init2(cos_x, prec);
    
    mpf_set(x, z.x);
    mpf_set(y, z.y);
    
    x_d = mpf_get_d(x);
    y_d = mpf_get_d(y);
    cosh_y_d = cosh(y_d);
    sinh_y_d = sinh(y_d);
    
    // rz.x = sin(x) * cosh(y)
    mpf_set_d(sin_x, sin(x_d));
    mpf_set_d(cosh_y, cosh_y_d);
    mpf_mul(rz.x, sin_x, cosh_y);
    
    // rz.y = cos(x) * sinh(y)
    mpf_set_d(cos_x, cos(x_d));
    mpf_set_d(sinh_y, sinh_y_d);
    mpf_mul(rz.y, cos_x, sinh_y);
    
    mpf_clear(x);
    mpf_clear(y);
    mpf_clear(cosh_y);
    mpf_clear(sinh_y);
    mpf_clear(sin_x);
    mpf_clear(cos_x);
    
    return rz;
}

complex_gmp complex_gmp_cos(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t x, y, cosh_y, sinh_y, cos_x, sin_x;
    double x_d, y_d, cosh_y_d, sinh_y_d;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(x, prec);
    mpf_init2(y, prec);
    mpf_init2(cosh_y, prec);
    mpf_init2(sinh_y, prec);
    mpf_init2(cos_x, prec);
    mpf_init2(sin_x, prec);
    
    mpf_set(x, z.x);
    mpf_set(y, z.y);
    
    x_d = mpf_get_d(x);
    y_d = mpf_get_d(y);
    cosh_y_d = cosh(y_d);
    sinh_y_d = sinh(y_d);
    
    // rz.x = cos(x) * cosh(y)
    mpf_set_d(cos_x, cos(x_d));
    mpf_set_d(cosh_y, cosh_y_d);
    mpf_mul(rz.x, cos_x, cosh_y);
    
    // rz.y = -sin(x) * sinh(y)
    mpf_set_d(sin_x, sin(x_d));
    mpf_set_d(sinh_y, sinh_y_d);
    mpf_mul(rz.y, sin_x, sinh_y);
    mpf_neg(rz.y, rz.y);
    
    mpf_clear(x);
    mpf_clear(y);
    mpf_clear(cosh_y);
    mpf_clear(sinh_y);
    mpf_clear(cos_x);
    mpf_clear(sin_x);
    
    return rz;
}

complex_gmp complex_gmp_tan(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp s, c, rz;
    s = complex_gmp_sin(z, prec);
    c = complex_gmp_cos(z, prec);
    rz = complex_gmp_div(s, c, prec);
    complex_gmp_clear(&s);
    complex_gmp_clear(&c);
    return rz;
}

complex_gmp complex_gmp_exp(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t exp_x, cos_y, sin_y;
    double x_d, y_d;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(exp_x, prec);
    mpf_init2(cos_y, prec);
    mpf_init2(sin_y, prec);
    
    x_d = mpf_get_d(z.x);
    y_d = mpf_get_d(z.y);
    
    mpf_set_d(exp_x, exp(x_d));
    mpf_set_d(cos_y, cos(y_d));
    mpf_set_d(sin_y, sin(y_d));
    
    // rz.x = exp(x) * cos(y)
    mpf_mul(rz.x, exp_x, cos_y);
    
    // rz.y = exp(x) * sin(y)
    mpf_mul(rz.y, exp_x, sin_y);
    
    mpf_clear(exp_x);
    mpf_clear(cos_y);
    mpf_clear(sin_y);
    
    return rz;
}

complex_gmp complex_gmp_sinh(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t x, y, sinh_x, cosh_x, cos_y, sin_y;
    double x_d, y_d;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(x, prec);
    mpf_init2(y, prec);
    mpf_init2(sinh_x, prec);
    mpf_init2(cosh_x, prec);
    mpf_init2(cos_y, prec);
    mpf_init2(sin_y, prec);
    
    x_d = mpf_get_d(z.x);
    y_d = mpf_get_d(z.y);
    
    mpf_set_d(sinh_x, sinh(x_d));
    mpf_set_d(cosh_x, cosh(x_d));
    mpf_set_d(cos_y, cos(y_d));
    mpf_set_d(sin_y, sin(y_d));
    
    // rz.x = sinh(x) * cos(y)
    mpf_mul(rz.x, sinh_x, cos_y);
    
    // rz.y = cosh(x) * sin(y)
    mpf_mul(rz.y, cosh_x, sin_y);
    
    mpf_clear(x);
    mpf_clear(y);
    mpf_clear(sinh_x);
    mpf_clear(cosh_x);
    mpf_clear(cos_y);
    mpf_clear(sin_y);
    
    return rz;
}

complex_gmp complex_gmp_cosh(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t x, y, sinh_x, cosh_x, cos_y, sin_y;
    double x_d, y_d;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(x, prec);
    mpf_init2(y, prec);
    mpf_init2(sinh_x, prec);
    mpf_init2(cosh_x, prec);
    mpf_init2(cos_y, prec);
    mpf_init2(sin_y, prec);
    
    x_d = mpf_get_d(z.x);
    y_d = mpf_get_d(z.y);
    
    mpf_set_d(sinh_x, sinh(x_d));
    mpf_set_d(cosh_x, cosh(x_d));
    mpf_set_d(cos_y, cos(y_d));
    mpf_set_d(sin_y, sin(y_d));
    
    // rz.x = cosh(x) * cos(y)
    mpf_mul(rz.x, cosh_x, cos_y);
    
    // rz.y = sinh(x) * sin(y)
    mpf_mul(rz.y, sinh_x, sin_y);
    
    mpf_clear(x);
    mpf_clear(y);
    mpf_clear(sinh_x);
    mpf_clear(cosh_x);
    mpf_clear(cos_y);
    mpf_clear(sin_y);
    
    return rz;
}

complex_gmp complex_gmp_tanh(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp sh, ch, rz;
    sh = complex_gmp_sinh(z, prec);
    ch = complex_gmp_cosh(z, prec);
    rz = complex_gmp_div(sh, ch, prec);
    complex_gmp_clear(&sh);
    complex_gmp_clear(&ch);
    return rz;
}

complex_gmp complex_gmp_log(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz;
    mpf_t mag, arg, log_mag;
    double mag_d, arg_d;
    
    complex_gmp_init(&rz, prec);
    mpf_init2(mag, prec);
    mpf_init2(arg, prec);
    mpf_init2(log_mag, prec);
    
    complex_gmp_mag(mag, z);
    
    if (mpf_cmp_ui(mag, 0) == 0) {
        complex_gmp_clear(&rz);
        rz = complex_gmp_nan(prec);
        mpf_clear(mag);
        mpf_clear(arg);
        mpf_clear(log_mag);
        THROW(EX_LOGZ);
        return rz;
    }
    
    mag_d = mpf_get_d(mag);
    complex_gmp_arg(arg, z);
    arg_d = mpf_get_d(arg);
    
    mpf_set_d(rz.x, log(mag_d));
    mpf_set_d(rz.y, arg_d);
    
    mpf_clear(mag);
    mpf_clear(arg);
    mpf_clear(log_mag);
    
    return rz;
}

complex_gmp complex_gmp_arcsin(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz, i, one, temp1, temp2, temp3;
    mpf_t i_val;
    
    i = complex_gmp_i(prec);
    one = complex_gmp_one(prec);
    
    // temp1 = i*z + sqrt(1 - z*z)
    temp1 = complex_gmp_mul(i, z, prec);
    temp2 = complex_gmp_mul(z, z, prec);
    temp3 = complex_gmp_sub(one, temp2, prec);
    temp2 = complex_gmp_sqrt(temp3, prec);
    temp3 = complex_gmp_add(temp1, temp2, prec);
    
    // rz = -i * log(temp3)
    rz = complex_gmp_log(temp3, prec);
    rz = complex_gmp_mul(i, rz, prec);
    rz = complex_gmp_neg(rz, prec);
    
    complex_gmp_clear(&i);
    complex_gmp_clear(&one);
    complex_gmp_clear(&temp1);
    complex_gmp_clear(&temp2);
    complex_gmp_clear(&temp3);
    
    return rz;
}

complex_gmp complex_gmp_arccos(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz, i, one, temp1, temp2;
    
    i = complex_gmp_i(prec);
    one = complex_gmp_one(prec);
    
    // temp1 = z + sqrt(z*z - 1)
    temp1 = complex_gmp_mul(z, z, prec);
    temp2 = complex_gmp_sub(temp1, one, prec);
    temp1 = complex_gmp_sqrt(temp2, prec);
    temp2 = complex_gmp_add(z, temp1, prec);
    
    // rz = -i * log(temp2)
    rz = complex_gmp_log(temp2, prec);
    rz = complex_gmp_mul(i, rz, prec);
    rz = complex_gmp_neg(rz, prec);
    
    complex_gmp_clear(&i);
    complex_gmp_clear(&one);
    complex_gmp_clear(&temp1);
    complex_gmp_clear(&temp2);
    
    return rz;
}

complex_gmp complex_gmp_arctan(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz, i, one, half_i, temp1, temp2, temp3;
    
    i = complex_gmp_i(prec);
    one = complex_gmp_one(prec);
    
    // half_i = i/2
    mpf_t half_val;
    mpf_init2(half_val, prec);
    mpf_set_ui(half_val, 1);
    mpf_div_ui(half_val, half_val, 2);
    half_i = complex_gmp_scalar_mul(half_val, i, prec);
    mpf_clear(half_val);
    
    // temp1 = (1 - i*z) / (1 + i*z)
    temp1 = complex_gmp_mul(i, z, prec);
    temp2 = complex_gmp_sub(one, temp1, prec);
    temp3 = complex_gmp_add(one, temp1, prec);
    temp1 = complex_gmp_div(temp2, temp3, prec);
    
    // rz = half_i * log(temp1)
    rz = complex_gmp_log(temp1, prec);
    rz = complex_gmp_mul(half_i, rz, prec);
    
    complex_gmp_clear(&i);
    complex_gmp_clear(&one);
    complex_gmp_clear(&half_i);
    complex_gmp_clear(&temp1);
    complex_gmp_clear(&temp2);
    complex_gmp_clear(&temp3);
    
    return rz;
}

complex_gmp complex_gmp_arcsinh(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz, one, temp1, temp2;
    
    one = complex_gmp_one(prec);
    
    // temp1 = z + sqrt(z*z + 1)
    temp1 = complex_gmp_mul(z, z, prec);
    temp2 = complex_gmp_add(temp1, one, prec);
    temp1 = complex_gmp_sqrt(temp2, prec);
    temp2 = complex_gmp_add(z, temp1, prec);
    
    // rz = log(temp2)
    rz = complex_gmp_log(temp2, prec);
    
    complex_gmp_clear(&one);
    complex_gmp_clear(&temp1);
    complex_gmp_clear(&temp2);
    
    return rz;
}

complex_gmp complex_gmp_arccosh(complex_gmp z, mp_bitcnt_t prec) {
    complex_gmp rz, one, temp1, temp2;
    
    one = complex_gmp_one(prec);
    
    // temp1 = z + sqrt(z*z - 1)
    temp1 = complex_gmp_mul(z, z, prec);
    temp2 = complex_gmp_sub(temp1, one, prec);
    temp1 = complex_gmp_sqrt(temp2, prec);
    temp2 = complex_gmp_add(z, temp1, prec);
    
    // rz = log(temp2)
    rz = complex_gmp_log(temp2, prec);
    
    complex_gmp_clear(&one);
    complex_gmp_clear(&temp1);
    complex_gmp_clear(&temp2);
    
    return rz;
}

// Utilitaires
int complex_gmp_is_zero(complex_gmp z) {
    return (mpf_cmp_ui(z.x, 0) == 0 && mpf_cmp_ui(z.y, 0) == 0);
}

int complex_gmp_is_inf(complex_gmp z) {
    // Vérification approximative (très grand nombre)
    mpf_t threshold;
    mpf_init2(threshold, mpf_get_prec(z.x));
    mpf_set_ui(threshold, 1);
    mpf_mul_2exp(threshold, threshold, 500); // 2^500
    
    int result = (mpf_cmp(z.x, threshold) > 0 || mpf_cmp(z.y, threshold) > 0);
    mpf_clear(threshold);
    return result;
}

int complex_gmp_is_nan(complex_gmp z) {
    // NaN simulé: vérifier si les deux parties sont zéro mais qu'on devrait avoir une valeur
    return 0; // Simplification
}

int complex_gmp_cmp_mag(mpf_t threshold, complex_gmp z) {
    mpf_t mag;
    mp_bitcnt_t prec = mpf_get_prec(z.x);
    int result;

    mpf_init2(mag, prec);
    complex_gmp_mag(mag, z);
    result = mpf_cmp(mag, threshold);
    mpf_clear(mag);

    return result;
}

// ============================================================================
// Fonctions in-place optimisées (Phase 1 - Optimisation GMP)
// ============================================================================

void gmp_mul_temps_init(gmp_mul_temps* temps, mp_bitcnt_t prec) {
    mpf_init2(temps->t1, prec);
    mpf_init2(temps->t2, prec);
    mpf_init2(temps->t3, prec);
    mpf_init2(temps->t4, prec);
    temps->initialized = 1;
}

void gmp_mul_temps_clear(gmp_mul_temps* temps) {
    if (temps->initialized) {
        mpf_clear(temps->t1);
        mpf_clear(temps->t2);
        mpf_clear(temps->t3);
        mpf_clear(temps->t4);
        temps->initialized = 0;
    }
}

void complex_gmp_add_to(complex_gmp* result, complex_gmp z1, complex_gmp z2) {
    mpf_add(result->x, z1.x, z2.x);
    mpf_add(result->y, z1.y, z2.y);
}

void complex_gmp_sub_to(complex_gmp* result, complex_gmp z1, complex_gmp z2) {
    mpf_sub(result->x, z1.x, z2.x);
    mpf_sub(result->y, z1.y, z2.y);
}

void complex_gmp_mul_to(complex_gmp* result, complex_gmp z1, complex_gmp z2, gmp_mul_temps* temps) {
    // result.x = z1.x * z2.x - z1.y * z2.y
    mpf_mul(temps->t1, z1.x, z2.x);
    mpf_mul(temps->t2, z1.y, z2.y);
    mpf_sub(result->x, temps->t1, temps->t2);

    // result.y = z1.x * z2.y + z1.y * z2.x
    mpf_mul(temps->t3, z1.x, z2.y);
    mpf_mul(temps->t4, z1.y, z2.x);
    mpf_add(result->y, temps->t3, temps->t4);
}

void complex_gmp_sq_to(complex_gmp* result, complex_gmp z, gmp_mul_temps* temps) {
    // Optimisé pour z² : result.x = z.x² - z.y², result.y = 2*z.x*z.y
    // Utilise seulement 3 multiplications au lieu de 4
    mpf_mul(temps->t1, z.x, z.x);      // z.x²
    mpf_mul(temps->t2, z.y, z.y);      // z.y²
    mpf_mul(temps->t3, z.x, z.y);      // z.x * z.y

    mpf_sub(result->x, temps->t1, temps->t2);  // z.x² - z.y²
    mpf_mul_ui(result->y, temps->t3, 2);       // 2 * z.x * z.y
}

void complex_gmp_mag2_to(mpf_t result, complex_gmp z, mpf_t temp1, mpf_t temp2) {
    // |z|² = z.x² + z.y² (sans sqrt)
    mpf_mul(temp1, z.x, z.x);
    mpf_mul(temp2, z.y, z.y);
    mpf_add(result, temp1, temp2);
}

void complex_gmp_copy_to(complex_gmp* dest, complex_gmp src) {
    mpf_set(dest->x, src.x);
    mpf_set(dest->y, src.y);
}
