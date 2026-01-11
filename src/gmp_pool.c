/*
 * gmp_pool.c
 * Pool de mémoire GMP pour optimiser les allocations
 * Released under GPL2
 * Copyleft 2024
 */

#include "config.h"

#ifdef HAVE_GMP

#include <stdlib.h>
#include <stdio.h>
#include "gmp_pool.h"

/* ============================================================================
 * Pool de mpf_t
 * ============================================================================ */

mpf_pool* mpf_pool_create(int capacity, mp_bitcnt_t prec) {
    int i;
    mpf_pool* pool = (mpf_pool*)malloc(sizeof(mpf_pool));
    if (!pool) {
        fprintf(stderr, "Erreur: impossible d'allouer mpf_pool\n");
        return NULL;
    }

    pool->items = (mpf_t*)malloc(capacity * sizeof(mpf_t));
    pool->in_use = (int*)calloc(capacity, sizeof(int));

    if (!pool->items || !pool->in_use) {
        fprintf(stderr, "Erreur: impossible d'allouer les éléments du pool\n");
        free(pool->items);
        free(pool->in_use);
        free(pool);
        return NULL;
    }

    pool->capacity = capacity;
    pool->count = 0;
    pool->precision = prec;

    /* Pré-allouer tous les éléments */
    for (i = 0; i < capacity; i++) {
        mpf_init2(pool->items[i], prec);
    }

    return pool;
}

void mpf_pool_destroy(mpf_pool* pool) {
    int i;
    if (!pool) return;

    for (i = 0; i < pool->capacity; i++) {
        mpf_clear(pool->items[i]);
    }

    free(pool->items);
    free(pool->in_use);
    free(pool);
}

mpf_t* mpf_pool_acquire(mpf_pool* pool) {
    int i;
    if (!pool) return NULL;

    /* Chercher un élément libre */
    for (i = 0; i < pool->capacity; i++) {
        if (!pool->in_use[i]) {
            pool->in_use[i] = 1;
            pool->count++;
            return &pool->items[i];
        }
    }

    /* Pool plein - allouer dynamiquement (fallback) */
    fprintf(stderr, "Warning: mpf_pool plein, allocation dynamique\n");
    mpf_t* new_item = (mpf_t*)malloc(sizeof(mpf_t));
    if (new_item) {
        mpf_init2(*new_item, pool->precision);
    }
    return new_item;
}

void mpf_pool_release(mpf_pool* pool, mpf_t* item) {
    int i;
    if (!pool || !item) return;

    /* Vérifier si l'item fait partie du pool */
    for (i = 0; i < pool->capacity; i++) {
        if (&pool->items[i] == item) {
            pool->in_use[i] = 0;
            pool->count--;
            return;
        }
    }

    /* Item alloué dynamiquement - le libérer */
    mpf_clear(*item);
    free(item);
}

void mpf_pool_reset(mpf_pool* pool) {
    int i;
    if (!pool) return;

    for (i = 0; i < pool->capacity; i++) {
        pool->in_use[i] = 0;
    }
    pool->count = 0;
}

/* ============================================================================
 * Pool de complex_gmp
 * ============================================================================ */

complex_gmp_pool* complex_gmp_pool_create(int capacity, mp_bitcnt_t prec) {
    int i;
    complex_gmp_pool* pool = (complex_gmp_pool*)malloc(sizeof(complex_gmp_pool));
    if (!pool) {
        fprintf(stderr, "Erreur: impossible d'allouer complex_gmp_pool\n");
        return NULL;
    }

    pool->items = (complex_gmp*)malloc(capacity * sizeof(complex_gmp));
    pool->in_use = (int*)calloc(capacity, sizeof(int));

    if (!pool->items || !pool->in_use) {
        fprintf(stderr, "Erreur: impossible d'allouer les éléments du pool\n");
        free(pool->items);
        free(pool->in_use);
        free(pool);
        return NULL;
    }

    pool->capacity = capacity;
    pool->count = 0;
    pool->precision = prec;

    /* Pré-allouer tous les éléments */
    for (i = 0; i < capacity; i++) {
        complex_gmp_init(&pool->items[i], prec);
    }

    return pool;
}

void complex_gmp_pool_destroy(complex_gmp_pool* pool) {
    int i;
    if (!pool) return;

    for (i = 0; i < pool->capacity; i++) {
        complex_gmp_clear(&pool->items[i]);
    }

    free(pool->items);
    free(pool->in_use);
    free(pool);
}

complex_gmp* complex_gmp_pool_acquire(complex_gmp_pool* pool) {
    int i;
    if (!pool) return NULL;

    /* Chercher un élément libre */
    for (i = 0; i < pool->capacity; i++) {
        if (!pool->in_use[i]) {
            pool->in_use[i] = 1;
            pool->count++;
            return &pool->items[i];
        }
    }

    /* Pool plein - allouer dynamiquement (fallback) */
    fprintf(stderr, "Warning: complex_gmp_pool plein, allocation dynamique\n");
    complex_gmp* new_item = (complex_gmp*)malloc(sizeof(complex_gmp));
    if (new_item) {
        complex_gmp_init(new_item, pool->precision);
    }
    return new_item;
}

void complex_gmp_pool_release(complex_gmp_pool* pool, complex_gmp* item) {
    int i;
    if (!pool || !item) return;

    /* Vérifier si l'item fait partie du pool */
    for (i = 0; i < pool->capacity; i++) {
        if (&pool->items[i] == item) {
            pool->in_use[i] = 0;
            pool->count--;
            return;
        }
    }

    /* Item alloué dynamiquement - le libérer */
    complex_gmp_clear(item);
    free(item);
}

void complex_gmp_pool_reset(complex_gmp_pool* pool) {
    int i;
    if (!pool) return;

    for (i = 0; i < pool->capacity; i++) {
        pool->in_use[i] = 0;
    }
    pool->count = 0;
}

/* ============================================================================
 * Contexte d'itération pré-alloué
 * ============================================================================ */

void gmp_iteration_context_init(gmp_iteration_context* ctx, mp_bitcnt_t prec) {
    if (!ctx) return;

    complex_gmp_init(&ctx->z, prec);
    complex_gmp_init(&ctx->zTemp, prec);
    complex_gmp_init(&ctx->zPixel, prec);
    gmp_mul_temps_init(&ctx->mul_temps, prec);
    mpf_init2(ctx->mag, prec);
    mpf_init2(ctx->mag_temps[0], prec);
    mpf_init2(ctx->mag_temps[1], prec);
    mpf_init2(ctx->bailout, prec);
    ctx->initialized = 1;
}

void gmp_iteration_context_clear(gmp_iteration_context* ctx) {
    if (!ctx || !ctx->initialized) return;

    complex_gmp_clear(&ctx->z);
    complex_gmp_clear(&ctx->zTemp);
    complex_gmp_clear(&ctx->zPixel);
    gmp_mul_temps_clear(&ctx->mul_temps);
    mpf_clear(ctx->mag);
    mpf_clear(ctx->mag_temps[0]);
    mpf_clear(ctx->mag_temps[1]);
    mpf_clear(ctx->bailout);
    ctx->initialized = 0;
}

void gmp_iteration_context_set_precision(gmp_iteration_context* ctx, mp_bitcnt_t prec) {
    if (!ctx || !ctx->initialized) return;

    complex_gmp_set_prec(&ctx->z, prec);
    complex_gmp_set_prec(&ctx->zTemp, prec);
    complex_gmp_set_prec(&ctx->zPixel, prec);
    mpf_set_prec(ctx->mul_temps.t1, prec);
    mpf_set_prec(ctx->mul_temps.t2, prec);
    mpf_set_prec(ctx->mul_temps.t3, prec);
    mpf_set_prec(ctx->mul_temps.t4, prec);
    mpf_set_prec(ctx->mag, prec);
    mpf_set_prec(ctx->mag_temps[0], prec);
    mpf_set_prec(ctx->mag_temps[1], prec);
    mpf_set_prec(ctx->bailout, prec);
}

#endif /* HAVE_GMP */
