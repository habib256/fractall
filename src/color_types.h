/* color_types.h
 * Common color type definitions for fractal rendering
 * Released under GPL2
 * Copyleft 2001-2003 VERHILLE Arnaud
 */

#ifndef COLOR_TYPES_H
#define COLOR_TYPES_H

/* Unified color type used throughout the application.
 * Using int for components to allow potential overflow handling
 * during color calculations before clamping.
 */
typedef struct {
    int r, g, b, a;
} color;

#endif /* COLOR_TYPES_H */
