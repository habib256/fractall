/*
 * precision_detector.c
 * Detection automatique de la précision nécessaire
 * Released under GPL2
 * Copyleft 2024
 */

#include <math.h>
#include "config.h"
#ifdef HAVE_GMP
#include <gmp.h>
#endif
#include "precision_detector.h"
#include "EscapeTime.h"

// Seuil de précision : quand la taille d'un pixel devient inférieure à 1e-14
// Cela correspond à environ 14 chiffres significatifs, limite pratique des double
#define PRECISION_THRESHOLD 1e-14

// Précision GMP minimale (64 bits)
#define GMP_PREC_MIN 64

// Précision GMP maximale (512 bits)
#define GMP_PREC_MAX 512

int precision_needs_gmp(fractal* f) {
    if (f == NULL) return 0;
    
    // Protection contre division par zéro
    if (f->xpixel <= 0 || f->ypixel <= 0) return 0;
    
    double pixel_size_x, pixel_size_y;
    
    // Calcul de la taille d'un pixel en unités complexes
    pixel_size_x = (f->xmax - f->xmin) / f->xpixel;
    pixel_size_y = (f->ymax - f->ymin) / f->ypixel;
    
    // Si la taille d'un pixel est inférieure au seuil, on a besoin de GMP
    if (pixel_size_x < PRECISION_THRESHOLD || pixel_size_y < PRECISION_THRESHOLD) {
        return 1;
    }
    
    return 0;
}

mp_bitcnt_t precision_calculate_gmp_bits(fractal* f) {
    if (f == NULL) return GMP_PREC_MIN;
    
    // Protection contre division par zéro
    if (f->xpixel <= 0 || f->ypixel <= 0) return GMP_PREC_MIN;
    
    double pixel_size_x, pixel_size_y;
    double min_pixel_size;
    double zoom_factor;
    mp_bitcnt_t prec;
    
    pixel_size_x = (f->xmax - f->xmin) / f->xpixel;
    pixel_size_y = (f->ymax - f->ymin) / f->ypixel;
    min_pixel_size = (pixel_size_x < pixel_size_y) ? pixel_size_x : pixel_size_y;
    
    // Calcul du facteur de zoom approximatif
    // Pour Mandelbrot standard: xmax-xmin ≈ 4.0
    // zoom_factor ≈ 4.0 / (xmax - xmin)
    double base_range = 4.0;
    // Protection contre division par zéro
    if (f->xmax - f->xmin <= 0.0) {
        return GMP_PREC_MIN;
    }
    zoom_factor = base_range / (f->xmax - f->xmin);
    
    // Précision basée sur le log2 du zoom
    // Formule: prec = 64 + log2(zoom_factor) * 8
    if (zoom_factor > 1.0) {
        prec = GMP_PREC_MIN + (mp_bitcnt_t)(log2(zoom_factor) * 8.0);
    } else {
        prec = GMP_PREC_MIN;
    }
    
    // Limiter entre min et max
    if (prec < GMP_PREC_MIN) prec = GMP_PREC_MIN;
    if (prec > GMP_PREC_MAX) prec = GMP_PREC_MAX;
    
    return prec;
}

void precision_update_fractal(fractal* f) {
    if (f == NULL) return;
    
    f->use_gmp = precision_needs_gmp(f);
    
    if (f->use_gmp) {
        f->gmp_precision = precision_calculate_gmp_bits(f);
    }
}

// Met à jour les structures GMP (mul_temps, iteration_ctx) après un changement de précision
void precision_update_gmp_structures(fractal* f) {
#ifdef HAVE_GMP
    if (f == NULL) return;
    
    if (f->use_gmp) {
        // Si mul_temps n'est pas initialisé, l'initialiser
        if (!f->mul_temps.initialized) {
            gmp_mul_temps_init(&f->mul_temps, f->gmp_precision);
        } else {
            // Si la précision a changé, mettre à jour mul_temps
            mpf_set_prec(f->mul_temps.t1, f->gmp_precision);
            mpf_set_prec(f->mul_temps.t2, f->gmp_precision);
            mpf_set_prec(f->mul_temps.t3, f->gmp_precision);
            mpf_set_prec(f->mul_temps.t4, f->gmp_precision);
        }
        
        // Si iteration_ctx n'est pas initialisé, l'initialiser
        if (!f->iteration_ctx.initialized) {
            gmp_iteration_context_init(&f->iteration_ctx, f->gmp_precision);
        } else {
            // Si la précision a changé, mettre à jour iteration_ctx
            gmp_iteration_context_set_precision(&f->iteration_ctx, f->gmp_precision);
        }
    } else {
        // Si use_gmp est désactivé, libérer les structures si elles sont initialisées
        if (f->mul_temps.initialized) {
            gmp_mul_temps_clear(&f->mul_temps);
        }
        if (f->iteration_ctx.initialized) {
            gmp_iteration_context_clear(&f->iteration_ctx);
        }
    }
#endif
}
