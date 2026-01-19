/*
 *     Complexmath.h
 *     Header file of Complex Number Trigonometric Functions
 *     based on Akihiro Sato work (27/Oct/1994)
 *
 *     Mostly rewritten by VERHILLE A. (2002-2003 GPL2)
 */

#include <math.h>
#include "exception.h"

#ifndef COMPLEXFLAG
#define COMPLEXFLAG /* COMPLEXFLAG */

typedef struct { // La structure d'un nombre complexe
   double x;
   double y;
}  complex;

#ifdef WIN32  /* POUR microsoft VC6 */
# define NAN 2745 // Attention A regler pour VC6 !!!
#endif  


#ifndef M_PI
#define M_PI	3.141592654 /* pi */
#endif  /* POUR microsoft VC6 */


#define NOT_A_NUMBER NAN

// Code d'erreur des exceptions complexes
#define EX_DIVZ 100
#define EX_POWZN 101
#define EX_TANZ 102
#define EX_TANHZ 103
#define EX_LOGZ 104
#define EX_ARCSINZ 105
#define EX_ARCCOSZ 106
#define EX_ARCTANZ 107
#define EX_ARCSINHZ 108
#define EX_ARCCOSHZ 109
#define EX_ARCTANHZ 110

#endif  /* COMPLEXFLAG */

// Set of Complex functions
extern complex MakeComplex(double, double); /* rz = real + i * imag */
extern complex ZeroSetofComplex(); /* rz = 0.0 + i * 0.0 */
extern complex ESetofComplex();   /* rz = 1.0 + i * 0.0 */
extern complex ISetofComplex();   /* rz = 0.0 + i * 1.0 */
extern complex InfSetofComplex(); /* rz = Infinity */
extern complex NaNSetofComplex(); /* rz = NaN + i * NaN */
extern complex SetofComplex(complex); /* rz = z */
extern complex ScalarTimesofComplex(double, complex); /* rz = c * z */

// Read and Neg Complex functions
extern double Rez(complex);
extern double Imz(complex);
extern complex Negz(complex);

// Operations and functions in Complex Space
extern complex Addz(complex, complex);
extern complex Subz(complex, complex);
extern complex Mulz(complex, complex);
extern complex  Divz(complex, complex);

extern double Magz(complex);
extern double Magz2(complex);
extern double Argz(complex);

extern double sign(double);

extern complex Sqrtz(complex);
extern complex Powz(complex, complex);
extern complex powzn(complex, complex, int);
extern complex sinz(complex);
extern complex cosz(complex);
extern complex tanz(complex);
extern complex expz(complex);
extern complex sinhz(complex);
extern complex coshz(complex);
extern complex tanhz(complex);
extern complex Logz(complex);
extern complex Arcsinz(complex);
extern complex Arccosz(complex);
extern complex Arctanz(complex);
extern complex Arcsinhz(complex);
extern complex Arccoshz(complex);
//extern complex Arctanhz(complex);

// SIMD-optimized batch operations
#ifdef HAVE_SSE4_1
void complex_mul_sse4(complex* result, const complex* a, const complex* b, int count);
void complex_add_sse4(complex* result, const complex* a, const complex* b, int count);
void complex_mag2_sse4(double* result, const complex* a, int count);
#endif

#ifdef HAVE_AVX
void complex_mul_avx(complex* result, const complex* a, const complex* b, int count);
#endif

/*  end  */
