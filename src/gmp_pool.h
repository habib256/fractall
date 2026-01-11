/*
 * gmp_pool.h
 * Pool de mémoire GMP pour optimiser les allocations
 * Released under GPL2
 * Copyleft 2024
 */

#ifndef GMP_POOL_H
#define GMP_POOL_H

#ifdef HAVE_GMP

#include <gmp.h>
#include "complexmath_gmp.h"

/* Pool de mpf_t pré-alloués */
typedef struct {
    mpf_t* items;           /* Tableau de mpf_t pré-alloués */
    int* in_use;            /* Tableau de flags d'utilisation */
    int capacity;           /* Taille du pool */
    int count;              /* Nombre d'éléments actuellement utilisés */
    mp_bitcnt_t precision;  /* Précision des éléments */
} mpf_pool;

/* Pool de complex_gmp pré-alloués */
typedef struct {
    complex_gmp* items;     /* Tableau de complex_gmp pré-alloués */
    int* in_use;            /* Tableau de flags d'utilisation */
    int capacity;           /* Taille du pool */
    int count;              /* Nombre d'éléments actuellement utilisés */
    mp_bitcnt_t precision;  /* Précision des éléments */
} complex_gmp_pool;

/* Création et destruction de pools mpf_t */
mpf_pool* mpf_pool_create(int capacity, mp_bitcnt_t prec);
void mpf_pool_destroy(mpf_pool* pool);

/* Acquisition et libération d'éléments mpf_t */
mpf_t* mpf_pool_acquire(mpf_pool* pool);
void mpf_pool_release(mpf_pool* pool, mpf_t* item);

/* Réinitialise le pool (marque tout comme disponible) */
void mpf_pool_reset(mpf_pool* pool);

/* Création et destruction de pools complex_gmp */
complex_gmp_pool* complex_gmp_pool_create(int capacity, mp_bitcnt_t prec);
void complex_gmp_pool_destroy(complex_gmp_pool* pool);

/* Acquisition et libération d'éléments complex_gmp */
complex_gmp* complex_gmp_pool_acquire(complex_gmp_pool* pool);
void complex_gmp_pool_release(complex_gmp_pool* pool, complex_gmp* item);

/* Réinitialise le pool (marque tout comme disponible) */
void complex_gmp_pool_reset(complex_gmp_pool* pool);

/* Structure pour contexte d'itération pré-alloué */
typedef struct {
    complex_gmp z;           /* Variable z courante */
    complex_gmp zTemp;       /* Temporaire pour calculs */
    complex_gmp zPixel;      /* Copie du pixel */
    gmp_mul_temps mul_temps; /* Temporaires pour multiplication */
    mpf_t mag;               /* Magnitude */
    mpf_t mag_temps[2];      /* Temporaires pour magnitude² */
    mpf_t bailout;           /* Seuil d'échappement */
    int initialized;
} gmp_iteration_context;

/* Création et destruction du contexte d'itération */
void gmp_iteration_context_init(gmp_iteration_context* ctx, mp_bitcnt_t prec);
void gmp_iteration_context_clear(gmp_iteration_context* ctx);
void gmp_iteration_context_set_precision(gmp_iteration_context* ctx, mp_bitcnt_t prec);

#endif /* HAVE_GMP */

#endif /* GMP_POOL_H */
