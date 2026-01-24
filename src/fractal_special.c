/* fractal_special.c
   Algorithmes speciaux (non escape-time standard)
   - Buddhabrot (densite trajectoires)
   - Lyapunov (exposant de Lyapunov)
   - Nebulabrot (RGB densite)
   released under GPL2
   Copyleft 2001-2026 VERHILLE Arnaud
*/

#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fractal_special.h"
#include "fractal_types.h"
#include "colorization.h"
#include "complexmath.h"
#include "SDL_gfxPrimitives.h"
#include "SDLGUI.h"
#include "config.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#ifdef HAVE_GMP
#include "complexmath_gmp.h"
#endif

// Fonction utilitaire pour traiter les evenements SDL pendant les calculs
// Retourne 1 si l'utilisateur veut annuler (ESC, Q, ou fermeture fenetre), 0 sinon
static int check_events_and_cancel(void) {
    SDL_Event events[32];
    int i;
    int should_cancel = 0;

    SDL_PumpEvents();

    int num_all = SDL_PeepEvents(events, 32, SDL_GETEVENT, SDL_QUITMASK | SDL_KEYDOWNMASK);

    for (i = 0; i < num_all; i++) {
        if (events[i].type == SDL_QUIT) {
            should_cancel = 1;
            continue;
        }

        if (events[i].type == SDL_KEYDOWN) {
            if (events[i].key.keysym.sym == SDLK_ESCAPE ||
                events[i].key.keysym.sym == SDLK_q) {
                should_cancel = 1;
                continue;
            }
        }

        SDL_PushEvent(&events[i]);
    }

    return should_cancel;
}

// *************************************************************************
// *    Buddhabrot
// *************************************************************************

Uint32 Buddhabrot_Draw (SDL_Surface *canvas, fractal* f, int decalageX, int decalageY, void* guiPtr) {
	int i, j;
	Uint32 time;
	int numSamples;
	int maxDensity;
	int xpixel = f->xpixel;
	int ypixel = f->ypixel;
	int iterMax = f->iterationMax;
	int bailout = f->bailout;

	double xmin, ymin, xmax, ymax;
#ifdef HAVE_GMP
	if (f->use_gmp) {
		mp_bitcnt_t prec = f->gmp_precision > 0 ? f->gmp_precision : 64;
		mpf_init2(f->xmin_gmp, prec);
		mpf_init2(f->xmax_gmp, prec);
		mpf_init2(f->ymin_gmp, prec);
		mpf_init2(f->ymax_gmp, prec);

		mpf_set_d(f->xmin_gmp, f->xmin);
		mpf_set_d(f->xmax_gmp, f->xmax);
		mpf_set_d(f->ymin_gmp, f->ymin);
		mpf_set_d(f->ymax_gmp, f->ymax);

		xmin = mpf_get_d(f->xmin_gmp);
		ymin = mpf_get_d(f->ymin_gmp);
		xmax = mpf_get_d(f->xmax_gmp);
		ymax = mpf_get_d(f->ymax_gmp);
	} else {
		xmin = f->xmin;
		ymin = f->ymin;
		xmax = f->xmax;
		ymax = f->ymax;
	}
#else
	xmin = f->xmin;
	ymin = f->ymin;
	xmax = f->xmax;
	ymax = f->ymax;
#endif

	double xrange = xmax - xmin;
	double yrange = ymax - ymin;

	if (f == NULL || f->fmatrix == NULL || f->cmatrix == NULL) {
		fprintf(stderr, "Erreur: structure fractale invalide pour Buddhabrot\n");
		return 0;
	}

	if (xrange <= 0.0 || yrange <= 0.0 || xpixel <= 0 || ypixel <= 0 || iterMax <= 0) {
		fprintf(stderr, "Erreur: parametres invalides pour Buddhabrot\n");
		return 0;
	}

	time = SDL_GetTicks();

	printf("Calculating Buddhabrot (density algorithm)...\n");
#ifdef HAVE_OPENMP
	printf("Using OpenMP with %d threads\n", omp_get_max_threads());
#endif

	for (i = 0; i < xpixel * ypixel; i++) {
		f->fmatrix[i] = 0;
	}

	int pixels = xpixel * ypixel;
	if (pixels <= 640*480) {
		numSamples = pixels * 20;
	} else if (pixels <= 1024*768) {
		numSamples = pixels * 10;
	} else {
		numSamples = pixels * 5;
	}

#ifdef HAVE_GMP
	if (f->use_gmp) {
		numSamples = numSamples / 10;
		if (numSamples < 1000) numSamples = 1000;
	}
#endif

	if (numSamples < 1000) numSamples = 1000;
	if (numSamples > 50000000) numSamples = 50000000;

	int chunk_size = 1000;
	if (chunk_size > numSamples / 10) {
		chunk_size = numSamples / 10;
	}
	if (chunk_size < 500) chunk_size = 500;
	int num_chunks = (numSamples + chunk_size - 1) / chunk_size;

	if (guiPtr != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, 0, "Buddhabrot");
	}

#ifdef HAVE_OPENMP
	volatile int cancelled = 0;
	volatile int sample_counter = 0;
	int num_threads = omp_get_max_threads();
	printf("Buddhabrot: %d samples, %d threads\n", numSamples, num_threads);

	#pragma omp parallel shared(cancelled, sample_counter)
	{
		int thread_id = omp_get_thread_num();
		unsigned int seed = (unsigned int)(42 + thread_id * 12345);

		double *trajX = (double*) malloc(iterMax * sizeof(double));
		double *trajY = (double*) malloc(iterMax * sizeof(double));

		if (trajX != NULL && trajY != NULL) {
			while (!cancelled) {
				int my_sample;
				#pragma omp atomic capture
				my_sample = sample_counter++;

				if (my_sample >= numSamples) {
					break;
				}

				if (thread_id == 0 && my_sample % 500 == 0) {
					if (check_events_and_cancel()) {
						printf("Calcul Buddhabrot annule par l'utilisateur\n");
						cancelled = 1;
						break;
					}
					if (guiPtr != NULL) {
						int percent = (my_sample * 90) / numSamples;
						if (percent > 90) percent = 90;
						SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, "Buddhabrot");
					}
				}

				seed = seed * 1103515245 + 12345;
				double xg = ((double)(seed & 0x7FFFFFFF) / 2147483647.0) * xrange + xmin;
				seed = seed * 1103515245 + 12345;
				double yg = ((double)(seed & 0x7FFFFFFF) / 2147483647.0) * yrange + ymin;

				complex c = MakeComplex(xg, yg);
				complex z = ZeroSetofComplex();

				int escaped = 0;
				int iter;
				double mag2_z = 0.0;
				int early_exit_threshold = (iterMax < 50) ? iterMax / 2 : 50;

				for (iter = 0; iter < iterMax; iter++) {
					if (iter % 50 == 0 && cancelled) break;

					complex zTemp = Mulz(z, z);
					z = Addz(zTemp, c);

					double zx = Rez(z);
					double zy = Imz(z);
					if (isnan(zx) || isnan(zy) || isinf(zx) || isinf(zy)) {
						break;
					}

					mag2_z = zx * zx + zy * zy;

					if (iter == early_exit_threshold && mag2_z < 0.25) {
						break;
					}

					trajX[iter] = zx;
					trajY[iter] = zy;

					if (mag2_z > bailout * bailout) {
						escaped = 1;
						break;
					}
				}

				if (escaped && iter > 0 && !cancelled) {
					double scale_x = xpixel / xrange;
					double scale_y = ypixel / yrange;

					for (int traj_idx = 0; traj_idx < iter; traj_idx++) {
						if (isnan(trajX[traj_idx]) || isnan(trajY[traj_idx]) ||
						    isinf(trajX[traj_idx]) || isinf(trajY[traj_idx])) {
							continue;
						}

						double px_d = (trajX[traj_idx] - xmin) * scale_x;
						double py_d = (trajY[traj_idx] - ymin) * scale_y;

						if (px_d < 0.0 || px_d >= (double)xpixel ||
						    py_d < 0.0 || py_d >= (double)ypixel) {
							continue;
						}

						int px = (int)px_d;
						int py = (int)py_d;

						if (px >= 0 && px < xpixel && py >= 0 && py < ypixel) {
							int idx = py * xpixel + px;
							if (idx >= 0 && idx < xpixel * ypixel) {
								#pragma omp atomic
								f->fmatrix[idx]++;
							}
						}
					}
				}
			}
		}

		if (trajX != NULL) free(trajX);
		if (trajY != NULL) free(trajY);
	}

	if (cancelled) {
		printf("Buddhabrot: calcul interrompu a %d/%d echantillons\n", sample_counter, numSamples);
	}

	maxDensity = 1;
	for (i = 0; i < xpixel * ypixel; i++) {
		if (f->fmatrix[i] > maxDensity) {
			maxDensity = f->fmatrix[i];
		}
	}

	if (maxDensity > 0) {
		double logMaxDensity = log(1.0 + (double)maxDensity);
		if (logMaxDensity > 0.0 && !isnan(logMaxDensity) && !isinf(logMaxDensity)) {
			const gradient_table* gradient = Colorization_GetPalette(f->colorMode);
			if (gradient != NULL) {
				double inv_logMaxDensity = 1.0 / logMaxDensity;
				for (j = 0; j < ypixel; j++) {
					for (i = 0; i < xpixel; i++) {
						int idx = j * xpixel + i;
						if (idx >= 0 && idx < xpixel * ypixel) {
							double density = (double)f->fmatrix[idx];
							double normalized = log(1.0 + density) * inv_logMaxDensity;
							if (normalized < 0.0) normalized = 0.0;
							else if (normalized > 1.0) normalized = 1.0;
							if (isnan(normalized) || isinf(normalized)) normalized = 0.0;

							colorization_color cc = Gradient_Interpolate(gradient, normalized);
							f->cmatrix[idx].r = cc.r;
							f->cmatrix[idx].g = cc.g;
							f->cmatrix[idx].b = cc.b;
							f->cmatrix[idx].a = 255;
							pixelRGBA(canvas, (Sint16)(i + decalageX), (Sint16)(j + decalageY),
							          cc.r, cc.g, cc.b, 255);
						}
					}
				}
				SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);
			}
		}
	}
#else
	double *trajX = (double*) malloc(iterMax * sizeof(double));
	double *trajY = (double*) malloc(iterMax * sizeof(double));
	if (trajX == NULL || trajY == NULL) {
		fprintf(stderr, "Erreur allocation memoire trajectoire\n");
		return 0;
	}

	srand(42);

	for (int chunk = 0; chunk < num_chunks; chunk++) {
		if (check_events_and_cancel()) {
			printf("Calcul Buddhabrot annule par l'utilisateur\n");
			free(trajX);
			free(trajY);
			return SDL_GetTicks() - time;
		}

		int chunk_start = chunk * chunk_size;
		int chunk_end = (chunk_start + chunk_size < numSamples) ? chunk_start + chunk_size : numSamples;

		for (int sample = chunk_start; sample < chunk_end; sample++) {
			if (sample % 1000 == 0 && check_events_and_cancel()) {
				printf("Calcul Buddhabrot annule par l'utilisateur\n");
				free(trajX);
				free(trajY);
				return SDL_GetTicks() - time;
			}
			double xg = ((double)rand() / RAND_MAX) * xrange + xmin;
			double yg = ((double)rand() / RAND_MAX) * yrange + ymin;
			complex c = MakeComplex(xg, yg);
			complex z = ZeroSetofComplex();

			int escaped = 0;
			int iter;
			double mag2_z = 0.0;
			int early_exit_threshold = (iterMax < 50) ? iterMax / 2 : 50;

			for (iter = 0; iter < iterMax; iter++) {
				complex zTemp = Mulz(z, z);
				z = Addz(zTemp, c);

				double zx = Rez(z);
				double zy = Imz(z);
				if (isnan(zx) || isnan(zy) || isinf(zx) || isinf(zy)) {
					break;
				}

				mag2_z = zx * zx + zy * zy;

				if (iter == early_exit_threshold && mag2_z < 0.25) {
					break;
				}

				trajX[iter] = zx;
				trajY[iter] = zy;

				if (mag2_z > bailout * bailout) {
					escaped = 1;
					break;
				}
			}

			if (escaped && iter > 0) {
				double inv_xrange = 1.0 / xrange;
				double inv_yrange = 1.0 / yrange;

				for (int traj_idx = 0; traj_idx < iter; traj_idx++) {
					if (isnan(trajX[traj_idx]) || isnan(trajY[traj_idx]) ||
					    isinf(trajX[traj_idx]) || isinf(trajY[traj_idx])) {
						continue;
					}

					double px_d = (trajX[traj_idx] - xmin) * inv_xrange * xpixel;
					double py_d = (trajY[traj_idx] - ymin) * inv_yrange * ypixel;

					if (px_d < 0.0 || px_d >= (double)xpixel ||
					    py_d < 0.0 || py_d >= (double)ypixel) {
						continue;
					}

					int px = (int)px_d;
					int py = (int)py_d;

					if (px >= 0 && px < xpixel && py >= 0 && py < ypixel) {
						int idx = py * xpixel + px;
						if (idx >= 0 && idx < xpixel * ypixel) {
							f->fmatrix[idx]++;
						}
					}
				}
			}
		}

		if (chunk % 2 == 0 || chunk == num_chunks - 1) {
			maxDensity = 1;
			for (i = 0; i < xpixel * ypixel; i++) {
				if (f->fmatrix[i] > maxDensity) {
					maxDensity = f->fmatrix[i];
				}
			}

			if (maxDensity > 0) {
				double logMaxDensity = log(1.0 + (double)maxDensity);
				if (logMaxDensity > 0.0 && !isnan(logMaxDensity) && !isinf(logMaxDensity)) {
					const gradient_table* gradient = Colorization_GetPalette(f->colorMode);
					if (gradient != NULL) {
						double inv_logMaxDensity = 1.0 / logMaxDensity;

						for (j = 0; j < ypixel; j++) {
							for (i = 0; i < xpixel; i++) {
								int idx = j * xpixel + i;
								if (idx >= 0 && idx < xpixel * ypixel) {
									double density = (double)f->fmatrix[idx];
									double normalized = log(1.0 + density) * inv_logMaxDensity;

									if (normalized < 0.0) normalized = 0.0;
									else if (normalized > 1.0) normalized = 1.0;
									if (isnan(normalized) || isinf(normalized)) normalized = 0.0;

									colorization_color cc = Gradient_Interpolate(gradient, normalized);
									f->cmatrix[idx].r = cc.r;
									f->cmatrix[idx].g = cc.g;
									f->cmatrix[idx].b = cc.b;
									f->cmatrix[idx].a = 255;
									pixelRGBA(canvas, (Sint16)(i + decalageX), (Sint16)(j + decalageY),
									          cc.r, cc.g, cc.b, 255);
								}
							}
							if (j % 20 == 0 || j == ypixel - 1) {
								SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);
							}
						}
					}
				}
			}
		}

		if (guiPtr != NULL) {
			int percent = ((chunk + 1) * 90) / num_chunks;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, "Buddhabrot");
		}

		if (check_events_and_cancel()) {
			printf("Calcul Buddhabrot annule par l'utilisateur\n");
			free(trajX);
			free(trajY);
			return SDL_GetTicks() - time;
		}
	}

	free(trajX);
	free(trajY);
#endif

	if (guiPtr != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, 100, "Buddhabrot");
	}

	SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);

	time = SDL_GetTicks() - time;
	printf("Buddhabrot rendered in %d ms (%d samples)\n", time, numSamples);
	return time;
}

// *************************************************************************
// *    Lyapunov Fractal
// *************************************************************************

#define LYAP_MAX_SEQ_LEN 64
#define LYAP_BLOCK_SIZE  64
#define LYAP_MIN_DERIV   1e-10
#define LYAP_MIN_X       0.0001
#define LYAP_MAX_X       0.9999

static double Lyapunov_GetNormalizedValue(double lyap) {
	if (lyap < 0) {
		double t = -lyap;
		if (t > 2.0) t = 2.0;
		return (t / 2.0) * 0.85;
	} else {
		double t = lyap;
		if (t > 1.0) t = 1.0;
		return 0.85 + (t * 0.15);
	}
}

static inline double Lyapunov_ComputeExponent_Optimized(
    double a, double b,
    const int* seq_is_a,
    int seqLen,
    int warmup,
    int iterMax)
{
    double x = 0.5;
    double lyap = 0.0;
    int n, block_count;

    int warmup_full_cycles = warmup / seqLen;
    int warmup_remainder = warmup % seqLen;

    for (int cycle = 0; cycle < warmup_full_cycles; cycle++) {
        for (int s = 0; s < seqLen; s++) {
            double r = seq_is_a[s] ? a : b;
            x = r * x * (1.0 - x);
        }
        if (x < LYAP_MIN_X || x > LYAP_MAX_X) x = 0.5;
    }
    for (int s = 0; s < warmup_remainder; s++) {
        double r = seq_is_a[s] ? a : b;
        x = r * x * (1.0 - x);
    }
    if (x < LYAP_MIN_X || x > LYAP_MAX_X) x = 0.5;

    int total_blocks = iterMax / LYAP_BLOCK_SIZE;
    int remainder = iterMax % LYAP_BLOCK_SIZE;
    int seq_idx = 0;

    for (block_count = 0; block_count < total_blocks; block_count++) {
        double product = 1.0;
        int valid_count = 0;

        for (n = 0; n < LYAP_BLOCK_SIZE; n++) {
            double r = seq_is_a[seq_idx] ? a : b;
            x = r * x * (1.0 - x);

            double deriv = fabs(r * (1.0 - 2.0 * x));
            if (deriv > LYAP_MIN_DERIV) {
                product *= deriv;
                valid_count++;
            }

            if (x < LYAP_MIN_X || x > LYAP_MAX_X) x = 0.5;

            seq_idx++;
            if (seq_idx >= seqLen) seq_idx = 0;
        }

        if (valid_count > 0 && product > 0.0) {
            lyap += log(product);
        }
    }

    if (remainder > 0) {
        double product = 1.0;
        int valid_count = 0;

        for (n = 0; n < remainder; n++) {
            double r = seq_is_a[seq_idx] ? a : b;
            x = r * x * (1.0 - x);

            double deriv = fabs(r * (1.0 - 2.0 * x));
            if (deriv > LYAP_MIN_DERIV) {
                product *= deriv;
                valid_count++;
            }

            if (x < LYAP_MIN_X || x > LYAP_MAX_X) x = 0.5;

            seq_idx++;
            if (seq_idx >= seqLen) seq_idx = 0;
        }

        if (valid_count > 0 && product > 0.0) {
            lyap += log(product);
        }
    }

    return (iterMax > 0) ? (lyap / iterMax) : 0.0;
}

#ifdef HAVE_SSE4_1
#include <smmintrin.h>

static inline void Lyapunov_ComputeExponent_SSE(
    double a0, double b0,
    double a1, double b1,
    const int* seq_is_a,
    int seqLen,
    int warmup,
    int iterMax,
    double* lyap_out0,
    double* lyap_out1)
{
    __m128d x = _mm_set1_pd(0.5);
    __m128d lyap = _mm_setzero_pd();
    __m128d one = _mm_set1_pd(1.0);
    __m128d two = _mm_set1_pd(2.0);
    __m128d min_deriv = _mm_set1_pd(LYAP_MIN_DERIV);
    __m128d min_x = _mm_set1_pd(LYAP_MIN_X);
    __m128d max_x = _mm_set1_pd(LYAP_MAX_X);
    __m128d half = _mm_set1_pd(0.5);

    __m128d a_vec = _mm_set_pd(a1, a0);
    __m128d b_vec = _mm_set_pd(b1, b0);

    for (int n = 0; n < warmup; n++) {
        int s = n % seqLen;
        __m128d r = seq_is_a[s] ? a_vec : b_vec;
        __m128d one_minus_x = _mm_sub_pd(one, x);
        x = _mm_mul_pd(_mm_mul_pd(r, x), one_minus_x);
    }
    __m128d too_low = _mm_cmplt_pd(x, min_x);
    __m128d too_high = _mm_cmpgt_pd(x, max_x);
    __m128d out_of_range = _mm_or_pd(too_low, too_high);
    x = _mm_blendv_pd(x, half, out_of_range);

    int total_blocks = iterMax / LYAP_BLOCK_SIZE;
    int remainder = iterMax % LYAP_BLOCK_SIZE;
    int seq_idx = 0;

    for (int block = 0; block < total_blocks; block++) {
        __m128d product = one;

        for (int n = 0; n < LYAP_BLOCK_SIZE; n++) {
            __m128d r = seq_is_a[seq_idx] ? a_vec : b_vec;
            __m128d one_minus_x = _mm_sub_pd(one, x);
            x = _mm_mul_pd(_mm_mul_pd(r, x), one_minus_x);

            __m128d two_x = _mm_mul_pd(two, x);
            __m128d one_minus_2x = _mm_sub_pd(one, two_x);
            __m128d deriv = _mm_mul_pd(r, one_minus_2x);
            __m128d abs_mask = _mm_castsi128_pd(_mm_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
            deriv = _mm_and_pd(deriv, abs_mask);

            __m128d valid = _mm_cmpgt_pd(deriv, min_deriv);
            __m128d safe_deriv = _mm_blendv_pd(one, deriv, valid);
            product = _mm_mul_pd(product, safe_deriv);

            too_low = _mm_cmplt_pd(x, min_x);
            too_high = _mm_cmpgt_pd(x, max_x);
            out_of_range = _mm_or_pd(too_low, too_high);
            x = _mm_blendv_pd(x, half, out_of_range);

            seq_idx++;
            if (seq_idx >= seqLen) seq_idx = 0;
        }

        double prod0, prod1;
        _mm_storel_pd(&prod0, product);
        _mm_storeh_pd(&prod1, product);
        if (prod0 > 0.0) lyap = _mm_add_pd(lyap, _mm_set_pd(0.0, log(prod0)));
        if (prod1 > 0.0) lyap = _mm_add_pd(lyap, _mm_set_pd(log(prod1), 0.0));
    }

    if (remainder > 0) {
        __m128d product = one;

        for (int n = 0; n < remainder; n++) {
            __m128d r = seq_is_a[seq_idx] ? a_vec : b_vec;
            __m128d one_minus_x = _mm_sub_pd(one, x);
            x = _mm_mul_pd(_mm_mul_pd(r, x), one_minus_x);

            __m128d two_x = _mm_mul_pd(two, x);
            __m128d one_minus_2x = _mm_sub_pd(one, two_x);
            __m128d deriv = _mm_mul_pd(r, one_minus_2x);
            __m128d abs_mask = _mm_castsi128_pd(_mm_set1_epi64x(0x7FFFFFFFFFFFFFFFLL));
            deriv = _mm_and_pd(deriv, abs_mask);

            __m128d valid = _mm_cmpgt_pd(deriv, min_deriv);
            __m128d safe_deriv = _mm_blendv_pd(one, deriv, valid);
            product = _mm_mul_pd(product, safe_deriv);

            too_low = _mm_cmplt_pd(x, min_x);
            too_high = _mm_cmpgt_pd(x, max_x);
            out_of_range = _mm_or_pd(too_low, too_high);
            x = _mm_blendv_pd(x, half, out_of_range);

            seq_idx++;
            if (seq_idx >= seqLen) seq_idx = 0;
        }

        double prod0, prod1;
        _mm_storel_pd(&prod0, product);
        _mm_storeh_pd(&prod1, product);
        if (prod0 > 0.0) lyap = _mm_add_pd(lyap, _mm_set_pd(0.0, log(prod0)));
        if (prod1 > 0.0) lyap = _mm_add_pd(lyap, _mm_set_pd(log(prod1), 0.0));
    }

    double inv_iterMax = (iterMax > 0) ? (1.0 / iterMax) : 0.0;
    __m128d scale = _mm_set1_pd(inv_iterMax);
    lyap = _mm_mul_pd(lyap, scale);

    _mm_storel_pd(lyap_out0, lyap);
    _mm_storeh_pd(lyap_out1, lyap);
}
#endif /* HAVE_SSE4_1 */

static Uint32 Lyapunov_Draw_Sequence (SDL_Surface *canvas, fractal* f, int decalageX, int decalageY, void* guiPtr, const char* sequence, const char* fractalName) {
	int i, j;
	double xStep, yStep;
	Uint32 time;
	int seqLen = strlen(sequence);
	int warmup = 50;
	int iterMax = f->iterationMax;
	int xpixel = f->xpixel;
	int ypixel = f->ypixel;

	int seq_is_a[LYAP_MAX_SEQ_LEN];
	if (seqLen > LYAP_MAX_SEQ_LEN) seqLen = LYAP_MAX_SEQ_LEN;
	for (int s = 0; s < seqLen; s++) {
		seq_is_a[s] = (sequence[s] == 'A' || sequence[s] == 'a') ? 1 : 0;
	}

	time = SDL_GetTicks();
	printf("Calculating %s fractal (sequence: %s, optimized v2)...\n", fractalName, sequence);
	printf("Domain: x=[%.6f, %.6f] y=[%.6f, %.6f]\n", f->xmin, f->xmax, f->ymin, f->ymax);
#ifdef HAVE_OPENMP
	printf("Using OpenMP with %d threads\n", omp_get_max_threads());
#endif
#ifdef HAVE_SSE4_1
	printf("Using SSE4.1 SIMD optimization (2 pixels/op)\n");
#endif

	xStep = (f->xmax - f->xmin) / xpixel;
	yStep = (f->ymax - f->ymin) / ypixel;

	if (guiPtr != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, 0, fractalName);
	}

	int chunk_height = 32;
	int cancelled = 0;

#ifdef HAVE_OPENMP
	for (int chunk_start = 0; chunk_start < ypixel && !cancelled; chunk_start += chunk_height) {
		int chunk_end = (chunk_start + chunk_height < ypixel) ? chunk_start + chunk_height : ypixel;

		#pragma omp parallel for private(i) schedule(dynamic, 8)
		for (j = chunk_start; j < chunk_end; j++) {
			double b = f->ymin + j * yStep;

#ifdef HAVE_SSE4_1
			for (i = 0; i < xpixel - 1; i += 2) {
				double a0 = f->xmin + i * xStep;
				double a1 = f->xmin + (i + 1) * xStep;
				double lyap0, lyap1;

				Lyapunov_ComputeExponent_SSE(a0, b, a1, b, seq_is_a, seqLen,
				                              warmup, iterMax, &lyap0, &lyap1);

				double norm0 = Lyapunov_GetNormalizedValue(lyap0);
				f->fmatrix[j * xpixel + i] = (int)(norm0 * iterMax);
				f->zmatrix[j * xpixel + i].x = norm0 * 2.0;
				f->zmatrix[j * xpixel + i].y = 0.0;

				double norm1 = Lyapunov_GetNormalizedValue(lyap1);
				f->fmatrix[j * xpixel + i + 1] = (int)(norm1 * iterMax);
				f->zmatrix[j * xpixel + i + 1].x = norm1 * 2.0;
				f->zmatrix[j * xpixel + i + 1].y = 0.0;
			}
			if (xpixel & 1) {
				i = xpixel - 1;
				double a = f->xmin + i * xStep;
				double lyap = Lyapunov_ComputeExponent_Optimized(a, b, seq_is_a, seqLen, warmup, iterMax);
				double norm = Lyapunov_GetNormalizedValue(lyap);
				f->fmatrix[j * xpixel + i] = (int)(norm * iterMax);
				f->zmatrix[j * xpixel + i].x = norm * 2.0;
				f->zmatrix[j * xpixel + i].y = 0.0;
			}
#else
			for (i = 0; i < xpixel; i++) {
				double a = f->xmin + i * xStep;

				double lyap = Lyapunov_ComputeExponent_Optimized(a, b, seq_is_a, seqLen, warmup, iterMax);
				double normalizedValue = Lyapunov_GetNormalizedValue(lyap);

				f->fmatrix[j * xpixel + i] = (int)(normalizedValue * iterMax);
				f->zmatrix[j * xpixel + i].x = normalizedValue * 2.0;
				f->zmatrix[j * xpixel + i].y = 0.0;
			}
#endif
		}

		if (check_events_and_cancel()) {
			cancelled = 1;
			printf("Calcul Lyapunov annule par l'utilisateur\n");
		}
		if (guiPtr != NULL && !cancelled) {
			int percent = ((chunk_start + chunk_height) * 70) / ypixel;
			if (percent > 70) percent = 70;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
		}
	}
#else
	for (int chunk_start = 0; chunk_start < ypixel && !cancelled; chunk_start += chunk_height) {
		int chunk_end = (chunk_start + chunk_height < ypixel) ? chunk_start + chunk_height : ypixel;

		for (j = chunk_start; j < chunk_end; j++) {
			double b = f->ymin + j * yStep;

#ifdef HAVE_SSE4_1
			for (i = 0; i < xpixel - 1; i += 2) {
				double a0 = f->xmin + i * xStep;
				double a1 = f->xmin + (i + 1) * xStep;
				double lyap0, lyap1;

				Lyapunov_ComputeExponent_SSE(a0, b, a1, b, seq_is_a, seqLen,
				                              warmup, iterMax, &lyap0, &lyap1);

				double norm0 = Lyapunov_GetNormalizedValue(lyap0);
				f->fmatrix[j * xpixel + i] = (int)(norm0 * iterMax);
				f->zmatrix[j * xpixel + i].x = norm0 * 2.0;
				f->zmatrix[j * xpixel + i].y = 0.0;

				double norm1 = Lyapunov_GetNormalizedValue(lyap1);
				f->fmatrix[j * xpixel + i + 1] = (int)(norm1 * iterMax);
				f->zmatrix[j * xpixel + i + 1].x = norm1 * 2.0;
				f->zmatrix[j * xpixel + i + 1].y = 0.0;
			}
			if (xpixel & 1) {
				i = xpixel - 1;
				double a = f->xmin + i * xStep;
				double lyap = Lyapunov_ComputeExponent_Optimized(a, b, seq_is_a, seqLen, warmup, iterMax);
				double norm = Lyapunov_GetNormalizedValue(lyap);
				f->fmatrix[j * xpixel + i] = (int)(norm * iterMax);
				f->zmatrix[j * xpixel + i].x = norm * 2.0;
				f->zmatrix[j * xpixel + i].y = 0.0;
			}
#else
			for (i = 0; i < xpixel; i++) {
				double a = f->xmin + i * xStep;

				double lyap = Lyapunov_ComputeExponent_Optimized(a, b, seq_is_a, seqLen, warmup, iterMax);
				double normalizedValue = Lyapunov_GetNormalizedValue(lyap);

				f->fmatrix[j * xpixel + i] = (int)(normalizedValue * iterMax);
				f->zmatrix[j * xpixel + i].x = normalizedValue * 2.0;
				f->zmatrix[j * xpixel + i].y = 0.0;
			}
#endif
		}

		if (check_events_and_cancel()) {
			cancelled = 1;
			printf("Calcul Lyapunov annule par l'utilisateur\n");
		}
		if (guiPtr != NULL && !cancelled) {
			int percent = ((chunk_start + chunk_height) * 70) / ypixel;
			if (percent > 70) percent = 70;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
		}
	}
#endif

	if (cancelled) {
		time = SDL_GetTicks() - time;
		return time;
	}

	f->cmatrix_valid = 0;
	int progress = 0;
	Fractal_CalculateColorMatrix(f, canvas, guiPtr, &progress, 70, 90);

	for (j = 0; j < ypixel && !cancelled; j++) {
		for (i = 0; i < xpixel; i++) {
			color col = f->cmatrix[j * xpixel + i];
			pixelRGBA(canvas, (Sint16)(i + decalageX), (Sint16)(j + decalageY),
			          col.r, col.g, col.b, 255);
		}

		if (j % 50 == 0) {
			SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);
			if (guiPtr != NULL) {
				int percent = 90 + (j * 10) / ypixel;
				SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
			}
			if (check_events_and_cancel()) {
				cancelled = 1;
				printf("Rendu Lyapunov annule par l'utilisateur\n");
			}
		}
	}

	SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);

	time = SDL_GetTicks() - time;
	printf("%s fractal rendered in %d ms\n", fractalName, time);
	return time;
}

Uint32 Lyapunov_Draw (SDL_Surface *canvas, fractal* f, int decalageX, int decalageY, void* guiPtr) {
	return Lyapunov_Draw_Sequence(canvas, f, decalageX, decalageY, guiPtr, "BBBBBBAAAAAA", "Lyapunov Zircon City");
}

// *************************************************************************
// *    Nebulabrot - Methode RGB avec 3 limites d'iterations
// *************************************************************************

Uint32 Nebulabrot_Draw (SDL_Surface *canvas, fractal* f, int decalageX, int decalageY, void* guiPtr) {
	int i, j;
	Uint32 time;
	int xpixel = f->xpixel;
	int ypixel = f->ypixel;
	int bailout = f->bailout;

	const int ITER_R = 50;
	const int ITER_G = 500;
	const int ITER_B = 5000;
	const int ITER_MAX = ITER_B;

	double xmin = f->xmin;
	double ymin = f->ymin;
	double xmax = f->xmax;
	double ymax = f->ymax;
	double xrange = xmax - xmin;
	double yrange = ymax - ymin;

	if (f == NULL || f->fmatrix == NULL || f->cmatrix == NULL) {
		fprintf(stderr, "Erreur: structure fractale invalide pour Nebulabrot\n");
		return 0;
	}

	if (xrange <= 0.0 || yrange <= 0.0 || xpixel <= 0 || ypixel <= 0) {
		fprintf(stderr, "Erreur: parametres invalides pour Nebulabrot\n");
		return 0;
	}

	time = SDL_GetTicks();
	printf("Calculating Nebulabrot (RGB density algorithm)...\n");
	printf("  Red channel: %d iterations\n", ITER_R);
	printf("  Green channel: %d iterations\n", ITER_G);
	printf("  Blue channel: %d iterations\n", ITER_B);

	int matrixSize = xpixel * ypixel;
	int *densityR = (int*) calloc(matrixSize, sizeof(int));
	int *densityG = (int*) calloc(matrixSize, sizeof(int));
	int *densityB = (int*) calloc(matrixSize, sizeof(int));

	if (densityR == NULL || densityG == NULL || densityB == NULL) {
		fprintf(stderr, "Erreur allocation memoire pour Nebulabrot\n");
		if (densityR) free(densityR);
		if (densityG) free(densityG);
		if (densityB) free(densityB);
		return 0;
	}

	int pixels = xpixel * ypixel;
	int numSamples;
	if (pixels <= 640*480) {
		numSamples = pixels * 15;
	} else if (pixels <= 1024*768) {
		numSamples = pixels * 8;
	} else {
		numSamples = pixels * 4;
	}
	if (numSamples < 1000) numSamples = 1000;
	if (numSamples > 30000000) numSamples = 30000000;

	printf("Nebulabrot: %d samples\n", numSamples);

	if (guiPtr != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, 0, "Nebulabrot");
	}

#ifdef HAVE_OPENMP
	volatile int cancelled = 0;
	volatile int sample_counter = 0;
	int num_threads = omp_get_max_threads();
	printf("Using OpenMP with %d threads\n", num_threads);

	#pragma omp parallel shared(cancelled, sample_counter, densityR, densityG, densityB)
	{
		int thread_id = omp_get_thread_num();
		unsigned int seed = (unsigned int)(42 + thread_id * 12345);

		double *trajX = (double*) malloc(ITER_MAX * sizeof(double));
		double *trajY = (double*) malloc(ITER_MAX * sizeof(double));

		if (trajX != NULL && trajY != NULL) {
			while (!cancelled) {
				int my_sample;
				#pragma omp atomic capture
				my_sample = sample_counter++;

				if (my_sample >= numSamples) break;

				if (thread_id == 0 && my_sample % 500 == 0) {
					if (check_events_and_cancel()) {
						printf("Calcul Nebulabrot annule par l'utilisateur\n");
						cancelled = 1;
						break;
					}
					if (guiPtr != NULL) {
						int percent = (my_sample * 85) / numSamples;
						if (percent > 85) percent = 85;
						SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, "Nebulabrot");
					}
				}

				seed = seed * 1103515245 + 12345;
				double xg = ((double)(seed & 0x7FFFFFFF) / 2147483647.0) * xrange + xmin;
				seed = seed * 1103515245 + 12345;
				double yg = ((double)(seed & 0x7FFFFFFF) / 2147483647.0) * yrange + ymin;

				complex c = MakeComplex(xg, yg);
				complex z = ZeroSetofComplex();

				int escaped = 0;
				int iter;
				double mag2_z = 0.0;

				for (iter = 0; iter < ITER_MAX; iter++) {
					if (iter % 100 == 0 && cancelled) break;

					complex zTemp = Mulz(z, z);
					z = Addz(zTemp, c);

					double zx = Rez(z);
					double zy = Imz(z);
					if (isnan(zx) || isnan(zy) || isinf(zx) || isinf(zy)) break;

					mag2_z = zx * zx + zy * zy;
					trajX[iter] = zx;
					trajY[iter] = zy;

					if (mag2_z > bailout * bailout) {
						escaped = 1;
						break;
					}
				}

				if (escaped && iter > 0 && !cancelled) {
					double scale_x = xpixel / xrange;
					double scale_y = ypixel / yrange;

					int contributeR = (iter <= ITER_R);
					int contributeG = (iter <= ITER_G);
					int contributeB = (iter <= ITER_B);

					for (int traj_idx = 0; traj_idx < iter; traj_idx++) {
						if (isnan(trajX[traj_idx]) || isnan(trajY[traj_idx])) continue;

						double px_d = (trajX[traj_idx] - xmin) * scale_x;
						double py_d = (trajY[traj_idx] - ymin) * scale_y;

						if (px_d < 0.0 || px_d >= (double)xpixel ||
						    py_d < 0.0 || py_d >= (double)ypixel) continue;

						int px = (int)px_d;
						int py = (int)py_d;

						if (px >= 0 && px < xpixel && py >= 0 && py < ypixel) {
							int idx = py * xpixel + px;
							if (idx >= 0 && idx < matrixSize) {
								if (contributeR) {
									#pragma omp atomic
									densityR[idx]++;
								}
								if (contributeG) {
									#pragma omp atomic
									densityG[idx]++;
								}
								if (contributeB) {
									#pragma omp atomic
									densityB[idx]++;
								}
							}
						}
					}
				}
			}
		}

		if (trajX != NULL) free(trajX);
		if (trajY != NULL) free(trajY);
	}

	if (cancelled) {
		printf("Nebulabrot: calcul interrompu a %d/%d echantillons\n", sample_counter, numSamples);
	}
#else
	double *trajX = (double*) malloc(ITER_MAX * sizeof(double));
	double *trajY = (double*) malloc(ITER_MAX * sizeof(double));

	if (trajX == NULL || trajY == NULL) {
		fprintf(stderr, "Erreur allocation memoire trajectoire\n");
		free(densityR); free(densityG); free(densityB);
		return 0;
	}

	srand(42);

	for (int sample = 0; sample < numSamples; sample++) {
		if (sample % 1000 == 0) {
			if (check_events_and_cancel()) {
				printf("Calcul Nebulabrot annule par l'utilisateur\n");
				break;
			}
			if (guiPtr != NULL) {
				int percent = (sample * 85) / numSamples;
				SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, "Nebulabrot");
			}
		}

		double xg = ((double)rand() / RAND_MAX) * xrange + xmin;
		double yg = ((double)rand() / RAND_MAX) * yrange + ymin;

		complex c = MakeComplex(xg, yg);
		complex z = ZeroSetofComplex();

		int escaped = 0;
		int iter;
		double mag2_z = 0.0;

		for (iter = 0; iter < ITER_MAX; iter++) {
			complex zTemp = Mulz(z, z);
			z = Addz(zTemp, c);

			double zx = Rez(z);
			double zy = Imz(z);
			if (isnan(zx) || isnan(zy) || isinf(zx) || isinf(zy)) break;

			mag2_z = zx * zx + zy * zy;
			trajX[iter] = zx;
			trajY[iter] = zy;

			if (mag2_z > bailout * bailout) {
				escaped = 1;
				break;
			}
		}

		if (escaped && iter > 0) {
			double scale_x = xpixel / xrange;
			double scale_y = ypixel / yrange;

			int contributeR = (iter <= ITER_R);
			int contributeG = (iter <= ITER_G);
			int contributeB = (iter <= ITER_B);

			for (int traj_idx = 0; traj_idx < iter; traj_idx++) {
				if (isnan(trajX[traj_idx]) || isnan(trajY[traj_idx])) continue;

				double px_d = (trajX[traj_idx] - xmin) * scale_x;
				double py_d = (trajY[traj_idx] - ymin) * scale_y;

				if (px_d < 0.0 || px_d >= (double)xpixel ||
				    py_d < 0.0 || py_d >= (double)ypixel) continue;

				int px = (int)px_d;
				int py = (int)py_d;

				if (px >= 0 && px < xpixel && py >= 0 && py < ypixel) {
					int idx = py * xpixel + px;
					if (idx >= 0 && idx < matrixSize) {
						if (contributeR) densityR[idx]++;
						if (contributeG) densityG[idx]++;
						if (contributeB) densityB[idx]++;
					}
				}
			}
		}
	}

	free(trajX);
	free(trajY);
#endif

	if (guiPtr != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, 90, "Nebulabrot");
	}

	int maxR = 1, maxG = 1, maxB = 1;
	for (i = 0; i < matrixSize; i++) {
		if (densityR[i] > maxR) maxR = densityR[i];
		if (densityG[i] > maxG) maxG = densityG[i];
		if (densityB[i] > maxB) maxB = densityB[i];
	}

	printf("Max densities: R=%d, G=%d, B=%d\n", maxR, maxG, maxB);

	double logMaxR = log(1.0 + (double)maxR);
	double logMaxG = log(1.0 + (double)maxG);
	double logMaxB = log(1.0 + (double)maxB);

	for (j = 0; j < ypixel; j++) {
		for (i = 0; i < xpixel; i++) {
			int idx = j * xpixel + i;

			double normR = (logMaxR > 0) ? log(1.0 + (double)densityR[idx]) / logMaxR : 0;
			double normG = (logMaxG > 0) ? log(1.0 + (double)densityG[idx]) / logMaxG : 0;
			double normB = (logMaxB > 0) ? log(1.0 + (double)densityB[idx]) / logMaxB : 0;

			if (normR > 1.0) normR = 1.0;
			if (normG > 1.0) normG = 1.0;
			if (normB > 1.0) normB = 1.0;

			Uint8 r = (Uint8)(normR * 255);
			Uint8 g = (Uint8)(normG * 255);
			Uint8 b = (Uint8)(normB * 255);

			f->cmatrix[idx].r = r;
			f->cmatrix[idx].g = g;
			f->cmatrix[idx].b = b;
			f->cmatrix[idx].a = 255;

			pixelRGBA(canvas, (Sint16)(i + decalageX), (Sint16)(j + decalageY), r, g, b, 255);
		}

		if (j % 50 == 0 || j == ypixel - 1) {
			SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);
			if (guiPtr != NULL) {
				int percent = 90 + (j * 10) / ypixel;
				SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, "Nebulabrot");
			}
		}
	}

	SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);

	free(densityR);
	free(densityG);
	free(densityB);

	time = SDL_GetTicks() - time;
	printf("Nebulabrot rendered in %d ms (%d samples)\n", time, numSamples);
	return time;
}
