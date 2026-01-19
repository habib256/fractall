/*
 * complexmath_simd.c
 * SIMD-optimized complex number operations
 * Released under GPL2
 * Copyleft 2025
 */

#include "complexmath.h"
#include "config.h"

#ifdef HAVE_SSE4_1
#include <smmintrin.h>

/* Multiplication de complexes vectorisée SSE4.1 (2 complexes à la fois) */
void complex_mul_sse4(complex* result, const complex* a, const complex* b, int count) {
    int i;
    for (i = 0; i < count - 1; i += 2) {
        // Charger 2 complexes de a et b
        __m128d a0 = _mm_loadu_pd((double*)&a[i]);
        __m128d a1 = _mm_loadu_pd((double*)&a[i+1]);
        __m128d b0 = _mm_loadu_pd((double*)&b[i]);
        __m128d b1 = _mm_loadu_pd((double*)&b[i+1]);
        
        // a0 = [a[i].x, a[i].y], b0 = [b[i].x, b[i].y]
        // Pour complex mul: (a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x)
        
        // Extraire les composantes
        __m128d a0_x = _mm_shuffle_pd(a0, a0, 0x0);  // [a[i].x, a[i].x]
        __m128d a0_y = _mm_shuffle_pd(a0, a0, 0x3);  // [a[i].y, a[i].y]
        __m128d b0_xy = _mm_shuffle_pd(b0, b0, 0x1); // [b[i].y, b[i].x]
        
        // Calcul: a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x
        __m128d mul1 = _mm_mul_pd(a0_x, b0);
        __m128d mul2 = _mm_mul_pd(a0_y, b0_xy);
        __m128d res0 = _mm_addsub_pd(mul1, mul2);
        
        // Même chose pour le deuxième complexe
        __m128d a1_x = _mm_shuffle_pd(a1, a1, 0x0);
        __m128d a1_y = _mm_shuffle_pd(a1, a1, 0x3);
        __m128d b1_xy = _mm_shuffle_pd(b1, b1, 0x1);
        
        __m128d mul3 = _mm_mul_pd(a1_x, b1);
        __m128d mul4 = _mm_mul_pd(a1_y, b1_xy);
        __m128d res1 = _mm_addsub_pd(mul3, mul4);
        
        // Stocker les résultats
        _mm_storeu_pd((double*)&result[i], res0);
        _mm_storeu_pd((double*)&result[i+1], res1);
    }
    
    // Traiter les éléments restants avec la version scalaire
    for (; i < count; i++) {
        result[i] = Mulz(a[i], b[i]);
    }
}

/* Addition de complexes vectorisée SSE4.1 (2 complexes à la fois) */
void complex_add_sse4(complex* result, const complex* a, const complex* b, int count) {
    int i;
    for (i = 0; i < count - 1; i += 2) {
        __m128d a0 = _mm_loadu_pd((double*)&a[i]);
        __m128d b0 = _mm_loadu_pd((double*)&b[i]);
        __m128d a1 = _mm_loadu_pd((double*)&a[i+1]);
        __m128d b1 = _mm_loadu_pd((double*)&b[i+1]);
        
        __m128d res0 = _mm_add_pd(a0, b0);
        __m128d res1 = _mm_add_pd(a1, b1);
        
        _mm_storeu_pd((double*)&result[i], res0);
        _mm_storeu_pd((double*)&result[i+1], res1);
    }
    
    for (; i < count; i++) {
        result[i] = Addz(a[i], b[i]);
    }
}

/* Magnitude au carré vectorisée SSE4.1 (2 complexes à la fois) */
void complex_mag2_sse4(double* result, const complex* a, int count) {
    int i;
    for (i = 0; i < count - 1; i += 2) {
        __m128d a0 = _mm_loadu_pd((double*)&a[i]);
        __m128d a1 = _mm_loadu_pd((double*)&a[i+1]);
        
        // |z|² = x² + y²
        __m128d sq0 = _mm_mul_pd(a0, a0);
        __m128d sq1 = _mm_mul_pd(a1, a1);
        
        // Addition horizontale: [x², y²] -> [x²+y², x²+y²]
        __m128d sum0 = _mm_hadd_pd(sq0, sq0);
        __m128d sum1 = _mm_hadd_pd(sq1, sq1);
        
        result[i] = _mm_cvtsd_f64(sum0);
        result[i+1] = _mm_cvtsd_f64(sum1);
    }
    
    for (; i < count; i++) {
        result[i] = Magz2(a[i]);
    }
}

#endif /* HAVE_SSE4_1 */

#ifdef HAVE_AVX
#include <immintrin.h>

/* Multiplication de complexes vectorisée AVX (4 complexes à la fois) */
void complex_mul_avx(complex* result, const complex* a, const complex* b, int count) {
    int i;
    for (i = 0; i < count - 3; i += 4) {
        // Charger 4 complexes (8 doubles)
        __m256d a_vec = _mm256_loadu_pd((double*)&a[i]);
        __m256d b_vec = _mm256_loadu_pd((double*)&b[i]);
        
        // Réorganiser pour calcul: [a0.x, a0.y, a1.x, a1.y] -> besoin de [a0.x, a0.x, a0.y, a0.y]
        // Pour simplifier, on traite par paires
        __m128d a0 = _mm256_extractf128_pd(a_vec, 0);
        __m128d a1 = _mm256_extractf128_pd(a_vec, 1);
        __m128d b0 = _mm256_extractf128_pd(b_vec, 0);
        __m128d b1 = _mm256_extractf128_pd(b_vec, 1);
        
        // Utiliser la version SSE pour chaque paire
        __m128d a0_x = _mm_shuffle_pd(a0, a0, 0x0);
        __m128d a0_y = _mm_shuffle_pd(a0, a0, 0x3);
        __m128d b0_xy = _mm_shuffle_pd(b0, b0, 0x1);
        
        __m128d mul1 = _mm_mul_pd(a0_x, b0);
        __m128d mul2 = _mm_mul_pd(a0_y, b0_xy);
        __m128d res0 = _mm_addsub_pd(mul1, mul2);
        
        __m128d a1_x = _mm_shuffle_pd(a1, a1, 0x0);
        __m128d a1_y = _mm_shuffle_pd(a1, a1, 0x3);
        __m128d b1_xy = _mm_shuffle_pd(b1, b1, 0x1);
        
        __m128d mul3 = _mm_mul_pd(a1_x, b1);
        __m128d mul4 = _mm_mul_pd(a1_y, b1_xy);
        __m128d res1 = _mm_addsub_pd(mul3, mul4);
        
        // Recombiner
        __m256d res_vec = _mm256_set_m128d(res1, res0);
        _mm256_storeu_pd((double*)&result[i], res_vec);
    }
    
    // Traiter les éléments restants
    for (; i < count; i++) {
        result[i] = Mulz(a[i], b[i]);
    }
}

#endif /* HAVE_AVX */
