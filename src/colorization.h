/* colorization.h
 * Unified colorization system for fractal rendering
 * Released under GPL2
 * Copyleft 2001-2003 VERHILLE Arnaud
 */

#ifndef COLORIZATION_H
#define COLORIZATION_H

#include "SDL.h"
#include "color_types.h"

/* Forward declaration */
struct fractal_struct;
typedef struct fractal_struct fractal;
struct gui_struct;
typedef struct gui_struct gui;

/* Color stop for gradient definition */
typedef struct {
    float position;    /* [0.0, 1.0] */
    Uint8 r, g, b;
} color_stop;

#define MAX_COLOR_STOPS 16
#define NUM_PALETTES 9  /* 8 existing + 1 new SmoothCosmic */

/* Gradient table structure */
typedef struct {
    const char* name;
    int num_stops;
    color_stop stops[MAX_COLOR_STOPS];
} gradient_table;

/* Palette indices */
#define PALETTE_FIRE    0
#define PALETTE_OCEAN   1
#define PALETTE_FOREST  2
#define PALETTE_VIOLET  3
#define PALETTE_RAINBOW 4
#define PALETTE_SUNSET  5
#define PALETTE_PLASMA  6
#define PALETTE_ICE     7
#define PALETTE_COSMIC  8

/* Use unified color type from color_types.h */
typedef color colorization_color;

/* ******************
 * Public interface
 * ******************/

/* Initialize the colorization system */
void Colorization_Init(void);

/* Get a palette by index */
const gradient_table* Colorization_GetPalette(int index);

/* Get palette name by index */
const char* Colorization_GetPaletteName(int index);

/* Interpolate a color from a gradient at position t [0.0, 1.0] */
colorization_color Gradient_Interpolate(const gradient_table* g, double t);

/* Apply a gradient to a fractal's color matrix */
void Fractal_ApplyGradient(fractal* f, const gradient_table* g);

/* Calculate smooth iteration value for a pixel */
double Colorization_SmoothIteration(fractal* f, int i, int j);

/* HSV to RGB conversion */
colorization_color HSVtoRGB(double h, double s, double v);

/* Calculate color matrix for a fractal */
void Fractal_CalculateColorMatrix(fractal* f, SDL_Surface* canvas, void* gui, int* progress, int progressStart, int progressEnd);

/* Test color mode (debug) */
void FractalColorTest(fractal* f);

/* Read color matrix values */
colorization_color Fractal_ReadColorMatrix(fractal f, int i, int j);
int Fractal_ReadColorMatrixRed(fractal f, int i, int j);
int Fractal_ReadColorMatrixGreen(fractal f, int i, int j);
int Fractal_ReadColorMatrixBlue(fractal f, int i, int j);

#endif /* COLORIZATION_H */
