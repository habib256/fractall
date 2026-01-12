/*
 * precision_detector.h
 * Detection automatique de la précision nécessaire pour les calculs de fractales
 * Released under GPL2
 * Copyleft 2024
 */

#ifndef PRECISION_DETECTOR_H
#define PRECISION_DETECTOR_H

#include "config.h"
#include "EscapeTime.h"

// Vérifie si GMP est nécessaire pour la précision actuelle
// Retourne 1 si GMP est nécessaire, 0 sinon
int precision_needs_gmp(fractal* f);

// Calcule la précision GMP nécessaire en bits
// Retourne la précision recommandée
mp_bitcnt_t precision_calculate_gmp_bits(fractal* f);

// Met à jour le flag use_gmp dans la structure fractal
void precision_update_fractal(fractal* f);

// Met à jour les structures GMP (mul_temps, iteration_ctx) après un changement de précision
// À appeler après precision_update_fractal
void precision_update_gmp_structures(fractal* f);

#endif /* PRECISION_DETECTOR_H */
