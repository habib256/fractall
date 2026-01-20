/* colorization.c
 * Unified colorization system for fractal rendering
 * Released under GPL2
 * Copyleft 2001-2003 VERHILLE Arnaud
 */

#include <math.h>
#include <string.h>
#include "colorization.h"
#include "EscapeTime.h"
#include "SDLGUI.h"
#include "config.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/* ******************
 * Palette definitions
 * ******************/

/* Fire palette: Black -> Red -> Yellow -> White */
static const gradient_table GRADIENT_FIRE = {
    "SmoothFire",
    4,
    {
        { 0.00f, 0,   0,   0   },  /* Black */
        { 0.33f, 255, 0,   0   },  /* Red */
        { 0.66f, 255, 255, 0   },  /* Yellow */
        { 1.00f, 255, 255, 255 }   /* White */
    }
};

/* Ocean palette: Black -> Blue -> Cyan -> White */
static const gradient_table GRADIENT_OCEAN = {
    "SmoothOcean",
    4,
    {
        { 0.00f, 0,   0,   0   },  /* Black */
        { 0.33f, 0,   0,   255 },  /* Blue */
        { 0.66f, 0,   255, 255 },  /* Cyan */
        { 1.00f, 255, 255, 255 }   /* White */
    }
};

/* Forest palette: Black -> Dark Green -> Yellow/Light Green -> White */
static const gradient_table GRADIENT_FOREST = {
    "SmoothForest",
    4,
    {
        { 0.00f, 0,   0,   0   },  /* Black */
        { 0.33f, 0,   180, 0   },  /* Dark Green */
        { 0.66f, 200, 255, 0   },  /* Yellow/Light Green */
        { 1.00f, 255, 255, 255 }   /* White */
    }
};

/* Violet palette: Black -> Dark Violet -> Pink/Magenta -> White */
static const gradient_table GRADIENT_VIOLET = {
    "SmoothViolet",
    4,
    {
        { 0.00f, 0,   0,   0   },  /* Black */
        { 0.33f, 128, 0,   200 },  /* Dark Violet */
        { 0.66f, 255, 100, 255 },  /* Pink/Magenta */
        { 1.00f, 255, 255, 255 }   /* White */
    }
};

/* Rainbow palette: Red -> Orange -> Yellow -> Green -> Cyan -> Blue -> Violet */
static const gradient_table GRADIENT_RAINBOW = {
    "SmoothRainbow",
    7,
    {
        { 0.000f, 255, 0,   0   },  /* Red */
        { 0.166f, 255, 165, 0   },  /* Orange */
        { 0.333f, 255, 255, 0   },  /* Yellow */
        { 0.500f, 0,   255, 0   },  /* Green */
        { 0.666f, 0,   255, 255 },  /* Cyan */
        { 0.833f, 0,   0,   255 },  /* Blue */
        { 1.000f, 180, 0,   255 }   /* Violet */
    }
};

/* Sunset palette: Black -> Orange -> Red -> Violet -> Dark Blue */
static const gradient_table GRADIENT_SUNSET = {
    "SmoothSunset",
    5,
    {
        { 0.00f, 0,   0,   0   },  /* Black */
        { 0.25f, 255, 140, 0   },  /* Orange */
        { 0.50f, 255, 0,   0   },  /* Red */
        { 0.75f, 255, 0,   200 },  /* Violet */
        { 1.00f, 55,  0,   255 }   /* Dark Blue */
    }
};

/* Plasma palette (NEW): Blue -> Violet -> Pink -> Orange */
static const gradient_table GRADIENT_PLASMA = {
    "SmoothPlasma",
    4,
    {
        { 0.00f, 13,  8,   135 },  /* Deep Blue */
        { 0.33f, 126, 3,   168 },  /* Violet */
        { 0.66f, 240, 87,  100 },  /* Pink/Coral */
        { 1.00f, 240, 230, 50  }   /* Yellow/Orange */
    }
};

/* Ice palette (NEW): White -> Cyan -> Deep Blue -> Black */
static const gradient_table GRADIENT_ICE = {
    "SmoothIce",
    4,
    {
        { 0.00f, 255, 255, 255 },  /* White */
        { 0.33f, 150, 230, 255 },  /* Light Cyan */
        { 0.66f, 30,  90,  200 },  /* Deep Blue */
        { 1.00f, 5,   10,  30  }   /* Near Black */
    }
};

/* Array of all palettes */
static const gradient_table* palettes[NUM_PALETTES] = {
    &GRADIENT_FIRE,
    &GRADIENT_OCEAN,
    &GRADIENT_FOREST,
    &GRADIENT_VIOLET,
    &GRADIENT_RAINBOW,
    &GRADIENT_SUNSET,
    &GRADIENT_PLASMA,
    &GRADIENT_ICE
};

/* ******************
 * Public functions
 * ******************/

void Colorization_Init(void) {
    /* Nothing to initialize for now - palettes are static */
}

const gradient_table* Colorization_GetPalette(int index) {
    if (index < 0 || index >= NUM_PALETTES) {
        return palettes[0];  /* Default to Fire */
    }
    return palettes[index];
}

const char* Colorization_GetPaletteName(int index) {
    if (index < 0 || index >= NUM_PALETTES) {
        return palettes[0]->name;
    }
    return palettes[index]->name;
}

/* Interpolate between two colors */
static colorization_color interpolate_rgb(const color_stop* a, const color_stop* b, double t) {
    colorization_color c;
    double factor;

    if (b->position == a->position) {
        factor = 0.0;
    } else {
        factor = (t - a->position) / (b->position - a->position);
    }

    c.r = (int)(a->r + factor * (b->r - a->r));
    c.g = (int)(a->g + factor * (b->g - a->g));
    c.b = (int)(a->b + factor * (b->b - a->b));
    c.a = 255;

    /* Clamp values */
    if (c.r < 0) c.r = 0; if (c.r > 255) c.r = 255;
    if (c.g < 0) c.g = 0; if (c.g > 255) c.g = 255;
    if (c.b < 0) c.b = 0; if (c.b > 255) c.b = 255;

    return c;
}

colorization_color Gradient_Interpolate(const gradient_table* g, double t) {
    int i;
    colorization_color c;
    const double EPSILON = 1e-9;

    /* Clamp t to [0, 1] */
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    /* Handle edge cases with small epsilon for floating point comparison */
    if (t <= g->stops[0].position + EPSILON) {
        c.r = g->stops[0].r;
        c.g = g->stops[0].g;
        c.b = g->stops[0].b;
        c.a = 255;
        return c;
    }
    if (t >= g->stops[g->num_stops - 1].position - EPSILON) {
        c.r = g->stops[g->num_stops - 1].r;
        c.g = g->stops[g->num_stops - 1].g;
        c.b = g->stops[g->num_stops - 1].b;
        c.a = 255;
        return c;
    }

    /* Find the two stops that bracket t */
    for (i = 0; i < g->num_stops - 1; i++) {
        /* Use < for upper bound to avoid matching multiple segments */
        if (t >= g->stops[i].position - EPSILON && t < g->stops[i + 1].position + EPSILON) {
            return interpolate_rgb(&g->stops[i], &g->stops[i + 1], t);
        }
    }

    /* Fallback: interpolate in the last segment (should rarely happen) */
    return interpolate_rgb(&g->stops[g->num_stops - 2], &g->stops[g->num_stops - 1], t);
}

/* Calculate smooth iteration value for a pixel */
double Colorization_SmoothIteration(fractal* f, int i, int j) {
    int iteration = *((f->fmatrix)+((f->xpixel*j)+i));
    complex z;
    double mag, log_zn, nu, smooth;

    /* For Lyapunov type (17), return normalized value directly */
    if (f->type == 17) {
        return (double)iteration / f->iterationMax;
    }

    z = *((f->zmatrix)+((f->xpixel*j)+i));
    mag = Magz(z);

    /* Need mag > e^2 (~7.4) for log(log(mag)/2/log(2)) to be valid (positive argument)
     * Also handle iteration >= iterationMax case */
    if (iteration >= f->iterationMax || mag < 7.4) {
        return (double)iteration / f->iterationMax;
    }

    /* Smooth coloring formula - at this point mag >= 7.4 so log_zn/log(2) > 1 */
    log_zn = log(mag) / 2.0;
    nu = log(log_zn / log(2.0)) / log(2.0);

    /* Additional safety check for numerical edge cases */
    if (!isfinite(nu)) {
        return (double)iteration / f->iterationMax;
    }

    smooth = iteration + 1 - nu;

    /* Normalize to [0, 1] */
    if (smooth < 0) smooth = 0;
    if (smooth > f->iterationMax) smooth = f->iterationMax;
    return smooth / f->iterationMax;
}

void Fractal_ApplyGradient(fractal* f, const gradient_table* g) {
    int xpixel = f->xpixel;
    int ypixel = f->ypixel;
    int iterMax = f->iterationMax;
    int j;

#ifdef HAVE_OPENMP
    #pragma omp parallel
    {
        #pragma omp for schedule(static) nowait
        for (j = 0; j < ypixel; j++) {
#else
    for (j = 0; j < ypixel; j++) {
#endif
        int i;
        for (i = 0; i < xpixel; i++) {
            int iteration = f->fmatrix[xpixel * j + i];
            double t, repeatCount, cycle, t_repeat;
            colorization_color c;
            color fc;

            /* Points in set: black */
            if (iteration >= iterMax) {
                fc.r = 0;
                fc.g = 0;
                fc.b = 0;
                fc.a = 255;
                f->cmatrix[xpixel * j + i] = fc;
                continue;
            }

            t = Colorization_SmoothIteration(f, i, j);

            /* Clamp t to avoid discontinuity at exactly 1.0 */
            if (t >= 1.0) t = 0.999999;

            /* Apply color repeat */
            repeatCount = (double)f->colorRepeat;
            cycle = floor(t * repeatCount);
            t_repeat = fmod(t * repeatCount, 1.0);

            /* Handle edge case where fmod returns 0 for exact integers */
            if (t_repeat < 0.000001 && t > 0.0) {
                t_repeat = 0.999999;
                cycle -= 1.0;
                if (cycle < 0) cycle = 0;
            }

            /* Alternate direction to avoid harsh transitions */
            if ((int)cycle % 2 == 1) {
                t_repeat = 1.0 - t_repeat;
            }

            c = Gradient_Interpolate(g, t_repeat);

            /* Convert colorization_color to color */
            fc.r = c.r;
            fc.g = c.g;
            fc.b = c.b;
            fc.a = c.a;
            f->cmatrix[xpixel * j + i] = fc;
        }
    }
#ifdef HAVE_OPENMP
    } /* End of parallel region */
#endif
}
