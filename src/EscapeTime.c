/* EscapeTime.c
fractal handling src code
released under GPL2
Copyleft 2001-2003 VERHILLE Arnaud
*/


#define _POSIX_C_SOURCE 200112L  // Pour posix_memalign
#include <stdio.h>
#include <stdlib.h>  // Pour malloc, free et posix_memalign
#include <string.h>  // Pour memcpy
#include <math.h>     // Pour fabs()
#include "SDL_gfxPrimitives.h"
#include "EscapeTime.h"
#include "SDLGUI.h"   // Pour la barre de progression
#include "colorization.h"  // Unified colorization system
#ifdef HAVE_GMP
#include "precision_detector.h"
#endif

// Fonction utilitaire pour traiter les événements SDL pendant les calculs
// Retourne 1 si l'utilisateur veut annuler (ESC, Q, ou fermeture fenêtre), 0 sinon
static int check_events_and_cancel(void) {
    SDL_Event event;
    SDL_PumpEvents();
    while (SDL_PollEvent(&event)) {
        if (event.type == SDL_QUIT) {
            return 1;
        }
        if (event.type == SDL_KEYDOWN) {
            if (event.key.keysym.sym == SDLK_ESCAPE || event.key.keysym.sym == SDLK_q) {
                return 1;
            }
        }
    }
    return 0;
}

// Fractal Functions
// *****************



fractal Fractal_Init (int screenW, int screenH, int type) {
	fractal f;
	f.xpixel = screenW;
	f.ypixel = screenH;
	f.type = type;
	f.colorMode = 6;  // SmoothPlasma par défaut
	f.cmatrix_valid = 0;
	f.last_colorMode = -1;
	f.last_colorRepeat = -1;
	// Répétition par défaut selon le type : 20 pour escape-time, 2 pour Lyapunov
	if (type == 17) {
		f.colorRepeat = 2;  // 2 répétitions pour Lyapunov
	} else {
		f.colorRepeat = 20;  // 20 répétitions pour escape-time
	}
	f.zoom_level = 1.0;  // Niveau de zoom initial
	// Initialiser le cache
	f.cache.cache_valid = 0;
	f.cache.fmatrix_cached = NULL;
	f.cache.cmatrix_cached = NULL;
	f.cache.cache_xpixel = 0;
	f.cache.cache_ypixel = 0;
#ifdef HAVE_GMP
	f.use_gmp = 0;
	f.gmp_precision = 64;
	f.zmatrix_gmp = NULL;
	f.iteration_ctx.initialized = 0;
	f.mul_temps.initialized = 0;
#endif
	Fractal_ChangeType (&f, type);
	
	// Allocation mémoire pour les matrices
	// Note: posix_memalign peut ne pas être disponible sur tous les systèmes
	// On utilise malloc standard pour la compatibilité
	f.fmatrix = (int *) malloc ((f.xpixel*f.ypixel)* sizeof (int));
	f.zmatrix = (complex *) malloc ((f.xpixel*f.ypixel)* sizeof (complex));
	f.cmatrix = (color *) malloc ((f.xpixel*f.ypixel)* sizeof (color));
	if (f.fmatrix == NULL || f.zmatrix == NULL || f.cmatrix == NULL) {
		fprintf(stderr, "Erreur allocation mémoire fractale\n");
		exit(1);
	}
	
	// Tentative d'alignement mémoire pour SIMD (optionnel, ne bloque pas si échoue)
#ifdef HAVE_AVX
	// Essayer d'allouer avec alignement, mais ne pas échouer si impossible
	// Note: posix_memalign nécessite _POSIX_C_SOURCE >= 200112L
	// Pour l'instant, on garde malloc standard pour éviter les problèmes de compatibilité
	// L'alignement peut être fait manuellement si nécessaire avec __attribute__((aligned(32)))
#endif
#ifdef HAVE_GMP
	precision_update_fractal(&f);
	precision_update_gmp_structures(&f);
	// Les coordonnées GMP sont déjà initialisées dans Fractal_ChangeType
	// Ici, on met juste à jour la précision si nécessaire
	if (f.use_gmp) {
		// Mettre à jour la précision si nécessaire
		mpf_set_prec(f.xmin_gmp, f.gmp_precision);
		mpf_set_prec(f.xmax_gmp, f.gmp_precision);
		mpf_set_prec(f.ymin_gmp, f.gmp_precision);
		mpf_set_prec(f.ymax_gmp, f.gmp_precision);
		// Re-synchroniser avec la bonne précision
		mpf_set_d(f.xmin_gmp, f.xmin);
		mpf_set_d(f.xmax_gmp, f.xmax);
		mpf_set_d(f.ymin_gmp, f.ymin);
		mpf_set_d(f.ymax_gmp, f.ymax);

		f.zmatrix_gmp = (complex_gmp *) malloc ((f.xpixel*f.ypixel)* sizeof (complex_gmp));
		if (f.zmatrix_gmp == NULL) {
			fprintf(stderr, "Erreur allocation mémoire GMP fractale\n");
			exit(1);
		}
		// Initialiser tous les éléments GMP
		for (int i = 0; i < f.xpixel * f.ypixel; i++) {
			complex_gmp_init(&f.zmatrix_gmp[i], f.gmp_precision);
		}
		// Initialiser le contexte d'itération et les temporaires de multiplication
		gmp_iteration_context_init(&f.iteration_ctx, f.gmp_precision);
		gmp_mul_temps_init(&f.mul_temps, f.gmp_precision);
	}
	// Les coordonnées GMP sont toujours initialisées maintenant (voir ci-dessus)
#endif
	return f;
}

void Fractal_Destroy (fractal f) {
	// Libération mémoire (posix_memalign utilise free)
	free (f.fmatrix);
	free (f.zmatrix);
	free (f.cmatrix);
	// Libérer le cache
	if (f.cache.fmatrix_cached != NULL) {
		free(f.cache.fmatrix_cached);
	}
	if (f.cache.cmatrix_cached != NULL) {
		free(f.cache.cmatrix_cached);
	}
#ifdef HAVE_GMP
	if (f.zmatrix_gmp != NULL) {
		for (int i = 0; i < f.xpixel * f.ypixel; i++) {
			complex_gmp_clear(&f.zmatrix_gmp[i]);
		}
		free (f.zmatrix_gmp);
	}
	// Libérer le contexte d'itération et les temporaires
	gmp_iteration_context_clear(&f.iteration_ctx);
	gmp_mul_temps_clear(&f.mul_temps);
	// Libérer les coordonnées GMP
	mpf_clear(f.xmin_gmp);
	mpf_clear(f.xmax_gmp);
	mpf_clear(f.ymin_gmp);
	mpf_clear(f.ymax_gmp);
#endif
}

// Fonctions de synchronisation des coordonnées
#ifdef HAVE_GMP
static void fractal_sync_coords_double_to_gmp(fractal* f) {
	if (f->use_gmp) {
		mpf_set_d(f->xmin_gmp, f->xmin);
		mpf_set_d(f->xmax_gmp, f->xmax);
		mpf_set_d(f->ymin_gmp, f->ymin);
		mpf_set_d(f->ymax_gmp, f->ymax);
	}
}

static void fractal_sync_coords_gmp_to_double(fractal* f) {
	if (f->use_gmp) {
		f->xmin = mpf_get_d(f->xmin_gmp);
		f->xmax = mpf_get_d(f->xmax_gmp);
		f->ymin = mpf_get_d(f->ymin_gmp);
		f->ymax = mpf_get_d(f->ymax_gmp);
	}
}
#endif

void Fractal_CalculateMatrix (fractal* f) {
	int i, j;
	complex zPixel;
	fractalresult result;

	printf ("Calculating Full EscapeTimeFractal Matrix  ...\n");

#ifdef HAVE_GMP
	if (f->use_gmp) {
		complex_gmp zPixel_gmp;
		mpf_t xg, yg;
		mp_bitcnt_t prec = f->gmp_precision;
		
		mpf_init2(xg, prec);
		mpf_init2(yg, prec);
		
		for (i=0; i<f->xpixel; i++) {
			for (j=0; j<f->ypixel; j++) {
				// Calcul des coordonnées en GMP directement depuis les coordonnées GMP
				mpf_t range_x, range_y, step_x, step_y, i_mpf, j_mpf;
				mpf_init2(range_x, prec);
				mpf_init2(range_y, prec);
				mpf_init2(step_x, prec);
				mpf_init2(step_y, prec);
				mpf_init2(i_mpf, prec);
				mpf_init2(j_mpf, prec);
				
				// Utiliser les coordonnées GMP directement
				mpf_sub(range_x, f->xmax_gmp, f->xmin_gmp);
				mpf_sub(range_y, f->ymax_gmp, f->ymin_gmp);
				mpf_set_ui(i_mpf, i);
				mpf_set_ui(j_mpf, j);
				mpf_set_ui(step_x, f->xpixel);
				mpf_set_ui(step_y, f->ypixel);
				
				// xg = (xmax-xmin)/xpixel * i + xmin
				mpf_div(step_x, range_x, step_x);
				mpf_mul(xg, step_x, i_mpf);
				mpf_add(xg, xg, f->xmin_gmp);
				
				// yg = (ymax-ymin)/ypixel * j + ymin
				mpf_div(step_y, range_y, step_y);
				mpf_mul(yg, step_y, j_mpf);
				mpf_add(yg, yg, f->ymin_gmp);
				
				zPixel_gmp = complex_gmp_make(xg, yg, prec);
				result = FormulaSelector_GMP(*f, zPixel_gmp);
				
				*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
				*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
				if (f->zmatrix_gmp != NULL) {
					complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+i], result.z, prec);
				}

				complex_gmp_clear(&zPixel_gmp);
				mpf_clear(range_x);
				mpf_clear(range_y);
				mpf_clear(step_x);
				mpf_clear(step_y);
				mpf_clear(i_mpf);
				mpf_clear(j_mpf);
			}
		}
		
		mpf_clear(xg);
		mpf_clear(yg);
		return;
	}
#endif

	for (i=0; i<f->xpixel; i++) {
		for (j=0; j<f->ypixel; j++) {
			double xg, yg; 
			xg = ((f->xmax-f->xmin)/f->xpixel)*i + f->xmin;
			yg = ((f->ymax-f->ymin)/f->ypixel)*j + f->ymin;
			zPixel = MakeComplex (xg, yg);
			
			result = FormulaSelector (*f, zPixel);
			*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
			*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
		}
	}
}

#ifdef HAVE_GMP
static void Fractal_CalculateMatrix_DDp1_GMP (fractal* f, SDL_Surface* canvas, void* guiPtr, int* progress, int progressStart, int progressEnd, const char* fractalName) {
	int i, j;
	mp_bitcnt_t prec = f->gmp_precision;
	int totalPixels, currentPixel = 0;
	int lastPercent = -1;
	int progressInterval;
	int cancelled = 0;

	printf ("Calculating EscapeTimeFractal Matrix with DDp1 (GMP) ...\n");
#ifdef HAVE_OPENMP
	int num_threads = omp_get_max_threads();
	printf("Using OpenMP with %d threads\n", num_threads);
#endif

	// Calculer le nombre total de pixels à traiter (grille 2x2)
	totalPixels = ((f->xpixel + 1) / 2) * ((f->ypixel + 1) / 2);
	progressInterval = totalPixels / 100;
	if (progressInterval < 1) progressInterval = 1;

	// Calculer les valeurs de step en double (thread-safe)
	// Lire les coordonnées GMP une seule fois avant la région parallèle
	double xmin_double = mpf_get_d(f->xmin_gmp);
	double xmax_double = mpf_get_d(f->xmax_gmp);
	double ymin_double = mpf_get_d(f->ymin_gmp);
	double ymax_double = mpf_get_d(f->ymax_gmp);
	double step_x_double = (xmax_double - xmin_double) / f->xpixel;
	double step_y_double = (ymax_double - ymin_double) / f->ypixel;

#ifdef HAVE_OPENMP
	int chunk_height = num_threads * 8;  // Au moins 4 itérations par thread (car j += 2)
	if (chunk_height < 32) chunk_height = 32;  // Minimum 32 lignes

	#pragma omp parallel
	{
		// Variables thread-local pour éviter allocations répétées et problèmes de thread-safety
		complex_gmp zPixel_gmp;
		mpf_t xg, yg, i_mpf, j_mpf, step_x_local, step_y_local, xmin_local, ymin_local;
		fractalresult result_local;
		fractal f_local;  // Copie locale de la structure fractal pour thread-safety
		int i, j;  // Variables de boucle thread-local

		// Initialiser f_local à zéro pour éviter que les mpf_t ne pointent vers de la mémoire invalide
		memset(&f_local, 0, sizeof(fractal));

		complex_gmp_init(&zPixel_gmp, prec);
		mpf_init2(xg, prec);
		mpf_init2(yg, prec);
		mpf_init2(i_mpf, prec);
		mpf_init2(j_mpf, prec);
		mpf_init2(step_x_local, prec);
		mpf_init2(step_y_local, prec);
		mpf_init2(xmin_local, prec);
		mpf_init2(ymin_local, prec);

		// Copier manuellement tous les champs de la structure fractal SAUF les structures GMP
		// (mul_temps, iteration_ctx, et les coordonnées mpf_t) car elles contiennent des mpf_t
		// qui ne peuvent pas être copiés directement et doivent être initialisées avec de la nouvelle mémoire.
		f_local.xpixel = f->xpixel;
		f_local.ypixel = f->ypixel;
		f_local.xmin = f->xmin;
		f_local.xmax = f->xmax;
		f_local.ymin = f->ymin;
		f_local.ymax = f->ymax;
		f_local.seed = f->seed;
		f_local.iterationMax = f->iterationMax;
		f_local.bailout = f->bailout;
		f_local.zoomfactor = f->zoomfactor;
		f_local.type = f->type;
		f_local.colorMode = f->colorMode;
		f_local.cmatrix_valid = f->cmatrix_valid;
		f_local.last_colorMode = f->last_colorMode;
		f_local.colorRepeat = f->colorRepeat;
		f_local.last_colorRepeat = f->last_colorRepeat;
		f_local.zoom_level = f->zoom_level;
		f_local.fmatrix = f->fmatrix;  // Partagé en écriture (mais chaque thread écrit à des indices différents)
		f_local.zmatrix = f->zmatrix;  // Partagé en écriture
		f_local.cmatrix = f->cmatrix;  // Partagé en écriture
		f_local.cache = f->cache;  // Structure simple, peut être copiée
		f_local.use_gmp = f->use_gmp;
		f_local.gmp_precision = f->gmp_precision;
		f_local.zmatrix_gmp = f->zmatrix_gmp;  // Partagé en écriture

		// Initialiser les coordonnées GMP thread-local à partir des valeurs double
		// (on évite de copier depuis f->xmin_gmp car cela nécessiterait de lire des mpf_t partagés)
		mpf_init2(f_local.xmin_gmp, prec);
		mpf_init2(f_local.xmax_gmp, prec);
		mpf_init2(f_local.ymin_gmp, prec);
		mpf_init2(f_local.ymax_gmp, prec);
		mpf_set_d(f_local.xmin_gmp, xmin_double);
		mpf_set_d(f_local.xmax_gmp, xmax_double);
		mpf_set_d(f_local.ymin_gmp, ymin_double);
		mpf_set_d(f_local.ymax_gmp, ymax_double);

		// Initialiser directement les structures GMP thread-local (sans copier)
		gmp_mul_temps_init(&f_local.mul_temps, prec);
		gmp_iteration_context_init(&f_local.iteration_ctx, prec);

		// Copier les valeurs partagées dans les variables thread-local (une seule fois par thread)
		mpf_set_d(step_x_local, step_x_double);
		mpf_set_d(step_y_local, step_y_double);
		mpf_set_d(xmin_local, xmin_double);
		mpf_set_d(ymin_local, ymin_double);

		// Traitement par chunks pour permettre la mise à jour de la progression et l'annulation
		for (int chunk_start = 0; chunk_start < f->ypixel && !cancelled; chunk_start += chunk_height) {
			int chunk_end = (chunk_start + chunk_height < f->ypixel) ?
			                chunk_start + chunk_height : f->ypixel;

			#pragma omp for schedule(guided) nowait
			for (j = chunk_start; j < chunk_end; j += 2) {
				for (i = 0; i < f->xpixel; i += 2) {
					// Calcul des coordonnées dans des variables temporaires puis copie dans zPixel_gmp
					mpf_set_ui(i_mpf, i);
					mpf_set_ui(j_mpf, j);

					mpf_mul(xg, step_x_local, i_mpf);
					mpf_add(xg, xg, xmin_local);

					mpf_mul(yg, step_y_local, j_mpf);
					mpf_add(yg, yg, ymin_local);

					mpf_set(zPixel_gmp.x, xg);
					mpf_set(zPixel_gmp.y, yg);
					zPixel_gmp.is_nan = 0;

					result_local = FormulaSelector_GMP(f_local, zPixel_gmp);

					*((f->fmatrix)+((f->xpixel*j)+i)) = result_local.iteration;
					*((f->zmatrix)+((f->xpixel*j)+i)) = result_local.z;
					if (f->zmatrix_gmp != NULL) {
						complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+i], result_local.z, prec);
					}

					// On complete pour un preview tout en evitant Segfault
					// Calculer les coordonnées correctes pour chaque pixel adjacent
					if ((i+1) < f->xpixel) {
						// Calculer les coordonnées pour (i+1, j)
						mpf_set_ui(i_mpf, i+1);
						mpf_mul(xg, step_x_local, i_mpf);
						mpf_add(xg, xg, xmin_local);
						mpf_set_ui(j_mpf, j);
						mpf_mul(yg, step_y_local, j_mpf);
						mpf_add(yg, yg, ymin_local);
						mpf_set(zPixel_gmp.x, xg);
						mpf_set(zPixel_gmp.y, yg);
						zPixel_gmp.is_nan = 0;

						result_local = FormulaSelector_GMP(f_local, zPixel_gmp);
						*((f->fmatrix)+((f->xpixel*j)+(i+1))) = result_local.iteration;
						*((f->zmatrix)+((f->xpixel*j)+(i+1))) = result_local.z;
						if (f->zmatrix_gmp != NULL) {
							complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+(i+1)], result_local.z, prec);
						}
					}
					if ((j+1) < f->ypixel) {
						// Calculer les coordonnées pour (i, j+1)
						mpf_set_ui(i_mpf, i);
						mpf_mul(xg, step_x_local, i_mpf);
						mpf_add(xg, xg, xmin_local);
						mpf_set_ui(j_mpf, j+1);
						mpf_mul(yg, step_y_local, j_mpf);
						mpf_add(yg, yg, ymin_local);
						mpf_set(zPixel_gmp.x, xg);
						mpf_set(zPixel_gmp.y, yg);
						zPixel_gmp.is_nan = 0;

						result_local = FormulaSelector_GMP(f_local, zPixel_gmp);
						*((f->fmatrix)+((f->xpixel*(j+1))+i)) = result_local.iteration;
						*((f->zmatrix)+((f->xpixel*(j+1))+i)) = result_local.z;
						if (f->zmatrix_gmp != NULL) {
							complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*(j+1))+i], result_local.z, prec);
						}
					}
					if ( ((i+1) < f->xpixel) && ((j+1) < f->ypixel)) {
						// Calculer les coordonnées pour (i+1, j+1)
						mpf_set_ui(i_mpf, i+1);
						mpf_mul(xg, step_x_local, i_mpf);
						mpf_add(xg, xg, xmin_local);
						mpf_set_ui(j_mpf, j+1);
						mpf_mul(yg, step_y_local, j_mpf);
						mpf_add(yg, yg, ymin_local);
						mpf_set(zPixel_gmp.x, xg);
						mpf_set(zPixel_gmp.y, yg);
						zPixel_gmp.is_nan = 0;

						result_local = FormulaSelector_GMP(f_local, zPixel_gmp);
						*((f->fmatrix)+((f->xpixel*(j+1))+(i+1))) = result_local.iteration;
						*((f->zmatrix)+((f->xpixel*(j+1))+(i+1))) = result_local.z;
						if (f->zmatrix_gmp != NULL) {
							complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*(j+1))+(i+1)], result_local.z, prec);
						}
					}
				}
			}

			#pragma omp barrier

			#pragma omp master
			{
				// Mise à jour de la progression
				if (guiPtr != NULL && progress != NULL) {
					int percent = progressStart +
						((chunk_start + chunk_height) * (progressEnd - progressStart)) / f->ypixel;
					if (percent > progressEnd) percent = progressEnd;
					SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
					*progress = percent;
				}

				// Vérifier annulation utilisateur
				if (check_events_and_cancel()) {
					printf("Calcul GMP DDp1 annulé par l'utilisateur\n");
					cancelled = 1;
				}
			}

			#pragma omp barrier  // Assure que tous les threads voient cancelled mis à jour
		}

		// Libération des variables thread-local
		complex_gmp_clear(&zPixel_gmp);
		mpf_clear(xg);
		mpf_clear(yg);
		mpf_clear(i_mpf);
		mpf_clear(j_mpf);
		mpf_clear(step_x_local);
		mpf_clear(step_y_local);
		mpf_clear(xmin_local);
		mpf_clear(ymin_local);
		gmp_mul_temps_clear(&f_local.mul_temps);
		gmp_iteration_context_clear(&f_local.iteration_ctx);
		// Libérer les coordonnées GMP thread-local
		mpf_clear(f_local.xmin_gmp);
		mpf_clear(f_local.xmax_gmp);
		mpf_clear(f_local.ymin_gmp);
		mpf_clear(f_local.ymax_gmp);
	} // Fin de la région parallèle
#else
	// Version séquentielle (sans OpenMP)
	complex_gmp zPixel_gmp;
	mpf_t xg, yg, i_mpf, j_mpf, step_x_local, step_y_local, xmin_local, ymin_local;
	fractalresult result;

	complex_gmp_init(&zPixel_gmp, prec);
	mpf_init2(xg, prec);
	mpf_init2(yg, prec);
	mpf_init2(i_mpf, prec);
	mpf_init2(j_mpf, prec);
	mpf_init2(step_x_local, prec);
	mpf_init2(step_y_local, prec);
	mpf_init2(xmin_local, prec);
	mpf_init2(ymin_local, prec);
	
	mpf_set_d(step_x_local, step_x_double);
	mpf_set_d(step_y_local, step_y_double);
	mpf_set_d(xmin_local, xmin_double);
	mpf_set_d(ymin_local, ymin_double);

	for (j=0; j<f->ypixel && !cancelled; j=j+2) {
		for (i=0; i<f->xpixel && !cancelled; i=i+2) {
			// Mise à jour de la progression et traitement des événements
			if (currentPixel % progressInterval == 0) {
				if (guiPtr != NULL) {
					int percent = progressStart + ((currentPixel * (progressEnd - progressStart)) / totalPixels);
					if (percent != lastPercent && percent <= progressEnd) {
						SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
						*progress = percent;
						lastPercent = percent;
					}
				}
				// Vérifier si l'utilisateur veut annuler
				if (check_events_and_cancel()) {
					printf("Calcul annulé par l'utilisateur\n");
					cancelled = 1;
					break;
				}
			}
			currentPixel++;

			// Calcul des coordonnées dans des variables temporaires puis copie dans zPixel_gmp
			mpf_set_ui(i_mpf, i);
			mpf_set_ui(j_mpf, j);

			mpf_mul(xg, step_x_local, i_mpf);
			mpf_add(xg, xg, xmin_local);

			mpf_mul(yg, step_y_local, j_mpf);
			mpf_add(yg, yg, ymin_local);

			mpf_set(zPixel_gmp.x, xg);
			mpf_set(zPixel_gmp.y, yg);
			zPixel_gmp.is_nan = 0;

			result = FormulaSelector_GMP(*f, zPixel_gmp);

			*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
			*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
			if (f->zmatrix_gmp != NULL) {
				complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+i], result.z, prec);
			}

			// On complete pour un preview tout en evitant Segfault
			// Calculer les coordonnées correctes pour chaque pixel adjacent
			if ((i+1) < f->xpixel) {
				// Calculer les coordonnées pour (i+1, j)
				mpf_set_ui(i_mpf, i+1);
				mpf_mul(xg, step_x_local, i_mpf);
				mpf_add(xg, xg, xmin_local);
				mpf_set_ui(j_mpf, j);
				mpf_mul(yg, step_y_local, j_mpf);
				mpf_add(yg, yg, ymin_local);
				mpf_set(zPixel_gmp.x, xg);
				mpf_set(zPixel_gmp.y, yg);
				zPixel_gmp.is_nan = 0;
				
				result = FormulaSelector_GMP(*f, zPixel_gmp);
				*((f->fmatrix)+((f->xpixel*j)+(i+1))) = result.iteration;
				*((f->zmatrix)+((f->xpixel*j)+(i+1))) = result.z;
				if (f->zmatrix_gmp != NULL) {
					complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+(i+1)], result.z, prec);
				}
			}
			if ((j+1) < f->ypixel) {
				// Calculer les coordonnées pour (i, j+1)
				mpf_set_ui(i_mpf, i);
				mpf_mul(xg, step_x_local, i_mpf);
				mpf_add(xg, xg, xmin_local);
				mpf_set_ui(j_mpf, j+1);
				mpf_mul(yg, step_y_local, j_mpf);
				mpf_add(yg, yg, ymin_local);
				mpf_set(zPixel_gmp.x, xg);
				mpf_set(zPixel_gmp.y, yg);
				zPixel_gmp.is_nan = 0;
				
				result = FormulaSelector_GMP(*f, zPixel_gmp);
				*((f->fmatrix)+((f->xpixel*(j+1))+i)) = result.iteration;
				*((f->zmatrix)+((f->xpixel*(j+1))+i)) = result.z;
				if (f->zmatrix_gmp != NULL) {
					complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*(j+1))+i], result.z, prec);
				}
			}
			if ( ((i+1) < f->xpixel) && ((j+1) < f->ypixel)) {
				// Calculer les coordonnées pour (i+1, j+1)
				mpf_set_ui(i_mpf, i+1);
				mpf_mul(xg, step_x_local, i_mpf);
				mpf_add(xg, xg, xmin_local);
				mpf_set_ui(j_mpf, j+1);
				mpf_mul(yg, step_y_local, j_mpf);
				mpf_add(yg, yg, ymin_local);
				mpf_set(zPixel_gmp.x, xg);
				mpf_set(zPixel_gmp.y, yg);
				zPixel_gmp.is_nan = 0;
				
				result = FormulaSelector_GMP(*f, zPixel_gmp);
				*((f->fmatrix)+((f->xpixel*(j+1))+(i+1))) = result.iteration;
				*((f->zmatrix)+((f->xpixel*(j+1))+(i+1))) = result.z;
				if (f->zmatrix_gmp != NULL) {
					complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*(j+1))+(i+1)], result.z, prec);
				}
			}
		}
	}

	// Libération des variables GMP après les boucles
	complex_gmp_clear(&zPixel_gmp);
	mpf_clear(xg);
	mpf_clear(yg);
	mpf_clear(i_mpf);
	mpf_clear(j_mpf);
	mpf_clear(step_x_local);
	mpf_clear(step_y_local);
	mpf_clear(xmin_local);
	mpf_clear(ymin_local);
#endif
}
#endif

void Fractal_CalculateMatrix_DDp1 (fractal* f, SDL_Surface* canvas, void* guiPtr, int* progress, int progressStart, int progressEnd, const char* fractalName) {
	// Invalider cmatrix car la matrice d'itérations va changer
	f->cmatrix_valid = 0;

#ifdef HAVE_GMP
	if (f->use_gmp) {
		Fractal_CalculateMatrix_DDp1_GMP(f, canvas, guiPtr, progress, progressStart, progressEnd, fractalName);
		return;
	}
#endif
	int j;
	int xpixel = f->xpixel;
	int ypixel = f->ypixel;
	double xmin = f->xmin;
	double ymin = f->ymin;
	double xstep = (f->xmax - f->xmin) / xpixel;
	double ystep = (f->ymax - f->ymin) / ypixel;

	printf ("Calculating EscapeTimeFractal Matrix with DDp1 ...\n");
#ifdef HAVE_OPENMP
	int num_threads = omp_get_max_threads();
	printf("Using OpenMP with %d threads\n", num_threads);
	// Granularité adaptative basée sur la taille de l'image
	int chunk_size = (ypixel / 2) / (num_threads * 4);
	if (chunk_size < 1) chunk_size = 1;
	if (chunk_size > 64) chunk_size = 64;
#else
	int num_threads = 1;
#endif

	// Mise à jour initiale de la progression
	if (guiPtr != NULL && progress != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressStart, fractalName);
		*progress = progressStart;
	}

#ifdef HAVE_OPENMP
	#pragma omp parallel
	{
		// Variables thread-local pour éviter allocations répétées
		complex zPixel_local;
		fractalresult result_local;
		int thread_id = omp_get_thread_num();
		
		#pragma omp for schedule(guided) nowait
		for (j = 0; j < ypixel; j += 2) {
#else
	for (j = 0; j < ypixel; j += 2) {
		complex zPixel_local;
		fractalresult result_local;
#endif
		int i;
		for (i = 0; i < xpixel; i += 2) {
			double xg = xstep * i + xmin;
			double yg = ystep * j + ymin;
			zPixel_local = MakeComplex(xg, yg);

			result_local = FormulaSelector(*f, zPixel_local);
			f->fmatrix[xpixel * j + i] = result_local.iteration;
			f->zmatrix[xpixel * j + i] = result_local.z;

			// On complete pour un preview tout en evitant Segfault
			if ((i + 1) < xpixel) {
				f->fmatrix[xpixel * j + (i + 1)] = result_local.iteration;
				f->zmatrix[xpixel * j + (i + 1)] = result_local.z;
			}
			if ((j + 1) < ypixel) {
				f->fmatrix[xpixel * (j + 1) + i] = result_local.iteration;
				f->zmatrix[xpixel * (j + 1) + i] = result_local.z;
			}
			if (((i + 1) < xpixel) && ((j + 1) < ypixel)) {
				f->fmatrix[xpixel * (j + 1) + (i + 1)] = result_local.iteration;
				f->zmatrix[xpixel * (j + 1) + (i + 1)] = result_local.z;
			}
		}
	}
#ifdef HAVE_OPENMP
	} // Fin de la région parallèle
#endif

	// Mise à jour finale de la progression
	if (guiPtr != NULL && progress != NULL) {
		*progress = progressEnd;
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressEnd, fractalName);
	}
}

#ifdef HAVE_GMP
static void Fractal_CalculateMatrix_DDp2_GMP (fractal* f, SDL_Surface* canvas, void* guiPtr, int* progress, int progressStart, int progressEnd, const char* fractalName) {
	int i, starti, j, compteur;
	int up, down, left, right;
	complex zup, zdown, zleft, zright;
	mp_bitcnt_t prec = f->gmp_precision;
	int progressMid = progressStart + (progressEnd - progressStart) / 2;

	printf ("Calculating EscapeTimeFractal Matrix with DDp2 (GMP) ...\n");
#ifdef HAVE_OPENMP
	int num_threads = omp_get_max_threads();
	printf("Using OpenMP with %d threads\n", num_threads);
#endif

	// Calculer les valeurs de step en double (thread-safe)
	// Lire les coordonnées GMP une seule fois avant la région parallèle
	double xmin_double = mpf_get_d(f->xmin_gmp);
	double xmax_double = mpf_get_d(f->xmax_gmp);
	double ymin_double = mpf_get_d(f->ymin_gmp);
	double ymax_double = mpf_get_d(f->ymax_gmp);
	double step_x_double = (xmax_double - xmin_double) / f->xpixel;
	double step_y_double = (ymax_double - ymin_double) / f->ypixel;

	// Mise à jour initiale de la progression
	if (guiPtr != NULL && progress != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressStart, fractalName);
		*progress = progressStart;
	}

	// Pass 1: Calcul initial des pixels impairs (parallélisable)
#ifdef HAVE_OPENMP
	int chunk_height_p1 = num_threads * 8;  // Au moins 4 itérations par thread (car j += 2)
	if (chunk_height_p1 < 32) chunk_height_p1 = 32;  // Minimum 32 lignes
	int cancelled_p1 = 0;

	#pragma omp parallel
	{
		// Variables thread-local
		complex_gmp zPixel_gmp;
		mpf_t xg, yg, i_mpf, j_mpf, step_x_local, step_y_local, xmin_local, ymin_local;
		fractalresult result_local;
		fractal f_local;  // Copie locale de la structure fractal pour thread-safety
		int i, j;  // Variables de boucle thread-local

		// Initialiser f_local à zéro pour éviter que les mpf_t ne pointent vers de la mémoire invalide
		memset(&f_local, 0, sizeof(fractal));

		complex_gmp_init(&zPixel_gmp, prec);
		mpf_init2(xg, prec);
		mpf_init2(yg, prec);
		mpf_init2(i_mpf, prec);
		mpf_init2(j_mpf, prec);
		mpf_init2(step_x_local, prec);
		mpf_init2(step_y_local, prec);
		mpf_init2(xmin_local, prec);
		mpf_init2(ymin_local, prec);

		// Copier manuellement tous les champs de la structure fractal SAUF les structures GMP
		f_local.xpixel = f->xpixel;
		f_local.ypixel = f->ypixel;
		f_local.xmin = f->xmin;
		f_local.xmax = f->xmax;
		f_local.ymin = f->ymin;
		f_local.ymax = f->ymax;
		f_local.seed = f->seed;
		f_local.iterationMax = f->iterationMax;
		f_local.bailout = f->bailout;
		f_local.zoomfactor = f->zoomfactor;
		f_local.type = f->type;
		f_local.colorMode = f->colorMode;
		f_local.cmatrix_valid = f->cmatrix_valid;
		f_local.last_colorMode = f->last_colorMode;
		f_local.colorRepeat = f->colorRepeat;
		f_local.last_colorRepeat = f->last_colorRepeat;
		f_local.zoom_level = f->zoom_level;
		f_local.fmatrix = f->fmatrix;
		f_local.zmatrix = f->zmatrix;
		f_local.cmatrix = f->cmatrix;
		f_local.cache = f->cache;
		f_local.use_gmp = f->use_gmp;
		f_local.gmp_precision = f->gmp_precision;
		f_local.zmatrix_gmp = f->zmatrix_gmp;

		// Initialiser les coordonnées GMP thread-local
		mpf_init2(f_local.xmin_gmp, prec);
		mpf_init2(f_local.xmax_gmp, prec);
		mpf_init2(f_local.ymin_gmp, prec);
		mpf_init2(f_local.ymax_gmp, prec);
		mpf_set_d(f_local.xmin_gmp, xmin_double);
		mpf_set_d(f_local.xmax_gmp, xmax_double);
		mpf_set_d(f_local.ymin_gmp, ymin_double);
		mpf_set_d(f_local.ymax_gmp, ymax_double);

		// Initialiser directement les structures GMP thread-local (sans copier)
		gmp_mul_temps_init(&f_local.mul_temps, prec);
		gmp_iteration_context_init(&f_local.iteration_ctx, prec);

		// Copier les valeurs partagées dans les variables thread-local (une seule fois par thread)
		mpf_set_d(step_x_local, step_x_double);
		mpf_set_d(step_y_local, step_y_double);
		mpf_set_d(xmin_local, xmin_double);
		mpf_set_d(ymin_local, ymin_double);

		// Traitement par chunks pour permettre la mise à jour de la progression et l'annulation
		for (int chunk_start = 1; chunk_start < f->ypixel && !cancelled_p1; chunk_start += chunk_height_p1) {
			int chunk_end = (chunk_start + chunk_height_p1 < f->ypixel) ?
			                chunk_start + chunk_height_p1 : f->ypixel;

			#pragma omp for schedule(guided) nowait
			for (j = chunk_start; j < chunk_end; j += 2) {
				for (i = 1; i < f->xpixel; i += 2) {
					mpf_set_ui(i_mpf, i);
					mpf_set_ui(j_mpf, j);

					mpf_mul(xg, step_x_local, i_mpf);
					mpf_add(xg, xg, xmin_local);

					mpf_mul(yg, step_y_local, j_mpf);
					mpf_add(yg, yg, ymin_local);

					mpf_set(zPixel_gmp.x, xg);
					mpf_set(zPixel_gmp.y, yg);
					zPixel_gmp.is_nan = 0;

					result_local = FormulaSelector_GMP(f_local, zPixel_gmp);

					*((f->fmatrix)+((f->xpixel*j)+i)) = result_local.iteration;
					*((f->zmatrix)+((f->xpixel*j)+i)) = result_local.z;
					if (f->zmatrix_gmp != NULL) {
						complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+i], result_local.z, prec);
					}
				}
			}

			#pragma omp barrier

			#pragma omp master
			{
				// Mise à jour de la progression (Pass 1 = première moitié)
				if (guiPtr != NULL && progress != NULL) {
					int progressMid = progressStart + (progressEnd - progressStart) / 2;
					int percent = progressStart +
						((chunk_start + chunk_height_p1 - 1) * (progressMid - progressStart)) / f->ypixel;
					if (percent > progressMid) percent = progressMid;
					SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
					*progress = percent;
				}

				// Vérifier annulation utilisateur
				if (check_events_and_cancel()) {
					printf("Calcul GMP DDp2 Pass 1 annulé par l'utilisateur\n");
					cancelled_p1 = 1;
				}
			}

			#pragma omp barrier  // Assure que tous les threads voient cancelled_p1 mis à jour
		}

		// Libération des variables thread-local
		complex_gmp_clear(&zPixel_gmp);
		mpf_clear(xg);
		mpf_clear(yg);
		mpf_clear(i_mpf);
		mpf_clear(j_mpf);
		mpf_clear(step_x_local);
		mpf_clear(step_y_local);
		mpf_clear(xmin_local);
		mpf_clear(ymin_local);
		gmp_mul_temps_clear(&f_local.mul_temps);
		gmp_iteration_context_clear(&f_local.iteration_ctx);
		// Libérer les coordonnées GMP thread-local
		mpf_clear(f_local.xmin_gmp);
		mpf_clear(f_local.xmax_gmp);
		mpf_clear(f_local.ymin_gmp);
		mpf_clear(f_local.ymax_gmp);
	} // Fin de la région parallèle pour Pass 1

	// Vérifier si annulé avant de continuer
	if (cancelled_p1) {
		return;
	}
#else
	// Version séquentielle Pass 1
	{
		complex_gmp zPixel_gmp;
		mpf_t xg, yg, i_mpf, j_mpf, step_x_local, step_y_local, xmin_local, ymin_local;
		fractalresult result;

		complex_gmp_init(&zPixel_gmp, prec);
		mpf_init2(xg, prec);
		mpf_init2(yg, prec);
		mpf_init2(i_mpf, prec);
		mpf_init2(j_mpf, prec);
		mpf_init2(step_x_local, prec);
		mpf_init2(step_y_local, prec);
		mpf_init2(xmin_local, prec);
		mpf_init2(ymin_local, prec);

		mpf_set_d(step_x_local, step_x_double);
		mpf_set_d(step_y_local, step_y_double);
		mpf_set_d(xmin_local, xmin_double);
		mpf_set_d(ymin_local, ymin_double);

		for (j=1; j<f->ypixel; j=j+2) {
			for (i=1; i<f->xpixel; i=i+2) {
				mpf_set_ui(i_mpf, i);
				mpf_set_ui(j_mpf, j);

				mpf_mul(xg, step_x_local, i_mpf);
				mpf_add(xg, xg, xmin_local);

				mpf_mul(yg, step_y_local, j_mpf);
				mpf_add(yg, yg, ymin_local);

				zPixel_gmp = complex_gmp_make(xg, yg, prec);
				result = FormulaSelector_GMP(*f, zPixel_gmp);

				*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
				*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
				if (f->zmatrix_gmp != NULL) {
					complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+i], result.z, prec);
				}

				complex_gmp_clear(&zPixel_gmp);
			}
		}

		mpf_clear(xg);
		mpf_clear(yg);
		mpf_clear(i_mpf);
		mpf_clear(j_mpf);
		mpf_clear(step_x_local);
		mpf_clear(step_y_local);
		mpf_clear(xmin_local);
		mpf_clear(ymin_local);
		complex_gmp_clear(&zPixel_gmp);
	}
#endif

	// Mise à jour de la progression après Pass 1
	if (guiPtr != NULL && progress != NULL) {
		*progress = progressMid;
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressMid, fractalName);
	}

	// Pass 2: Divergence Detection (DD) - parallélisable (chaque pixel écrit à son propre indice)
#ifdef HAVE_OPENMP
	int chunk_height_p2 = num_threads * 4;  // Au moins 4 itérations par thread
	if (chunk_height_p2 < 32) chunk_height_p2 = 32;  // Minimum 32 lignes
	int cancelled_p2 = 0;

	#pragma omp parallel
	{
		// Variables thread-local
		complex_gmp zPixel_gmp;
		mpf_t xg, yg, i_mpf, j_mpf, step_x_local, step_y_local, xmin_local, ymin_local;
		fractalresult result_local;
		complex zup_local, zdown_local, zleft_local, zright_local;
		fractal f_local;  // Copie locale de la structure fractal pour thread-safety
		int up, down, left, right, compteur, starti;
		int i, j;  // Variables de boucle thread-local

		// Initialiser f_local à zéro pour éviter que les mpf_t ne pointent vers de la mémoire invalide
		memset(&f_local, 0, sizeof(fractal));

		complex_gmp_init(&zPixel_gmp, prec);
		mpf_init2(xg, prec);
		mpf_init2(yg, prec);
		mpf_init2(i_mpf, prec);
		mpf_init2(j_mpf, prec);
		mpf_init2(step_x_local, prec);
		mpf_init2(step_y_local, prec);
		mpf_init2(xmin_local, prec);
		mpf_init2(ymin_local, prec);

		// Copier manuellement tous les champs de la structure fractal SAUF les structures GMP
		f_local.xpixel = f->xpixel;
		f_local.ypixel = f->ypixel;
		f_local.xmin = f->xmin;
		f_local.xmax = f->xmax;
		f_local.ymin = f->ymin;
		f_local.ymax = f->ymax;
		f_local.seed = f->seed;
		f_local.iterationMax = f->iterationMax;
		f_local.bailout = f->bailout;
		f_local.zoomfactor = f->zoomfactor;
		f_local.type = f->type;
		f_local.colorMode = f->colorMode;
		f_local.cmatrix_valid = f->cmatrix_valid;
		f_local.last_colorMode = f->last_colorMode;
		f_local.colorRepeat = f->colorRepeat;
		f_local.last_colorRepeat = f->last_colorRepeat;
		f_local.zoom_level = f->zoom_level;
		f_local.fmatrix = f->fmatrix;
		f_local.zmatrix = f->zmatrix;
		f_local.cmatrix = f->cmatrix;
		f_local.cache = f->cache;
		f_local.use_gmp = f->use_gmp;
		f_local.gmp_precision = f->gmp_precision;
		f_local.zmatrix_gmp = f->zmatrix_gmp;

		// Initialiser les coordonnées GMP thread-local à partir des valeurs double
		mpf_init2(f_local.xmin_gmp, prec);
		mpf_init2(f_local.xmax_gmp, prec);
		mpf_init2(f_local.ymin_gmp, prec);
		mpf_init2(f_local.ymax_gmp, prec);
		mpf_set_d(f_local.xmin_gmp, xmin_double);
		mpf_set_d(f_local.xmax_gmp, xmax_double);
		mpf_set_d(f_local.ymin_gmp, ymin_double);
		mpf_set_d(f_local.ymax_gmp, ymax_double);

		// Initialiser directement les structures GMP thread-local (sans copier)
		gmp_mul_temps_init(&f_local.mul_temps, prec);
		gmp_iteration_context_init(&f_local.iteration_ctx, prec);

		// Copier les valeurs partagées dans les variables thread-local (une seule fois par thread)
		mpf_set_d(step_x_local, step_x_double);
		mpf_set_d(step_y_local, step_y_double);
		mpf_set_d(xmin_local, xmin_double);
		mpf_set_d(ymin_local, ymin_double);

		// Traitement par chunks pour permettre la mise à jour de la progression et l'annulation
		for (int chunk_start = 0; chunk_start < f->ypixel && !cancelled_p2; chunk_start += chunk_height_p2) {
			int chunk_end = (chunk_start + chunk_height_p2 < f->ypixel) ?
			                chunk_start + chunk_height_p2 : f->ypixel;

			#pragma omp for schedule(guided) nowait
			for (j = chunk_start; j < chunk_end; j++) {
				if (fmod(j,2) == 0) {starti=1;} else {starti=0;}
				for (i = starti; i < f->xpixel; i += 2) {
					// Vérifier les bords des matrices pour éviter Seg Fault
					compteur=0;
					if (j != 0) {
						up = *((f->fmatrix)+((f->xpixel*(j-1))+i));
						zup_local = *((f->zmatrix)+((f->xpixel*(j-1))+i));
						compteur++;
					} else { up = 0; zup_local = ZeroSetofComplex();}
					if (j < (f->ypixel-1)) {
						down = *((f->fmatrix)+((f->xpixel*(j+1))+i));
						zdown_local = *((f->zmatrix)+((f->xpixel*(j+1))+i));
						compteur++;
					} else { down = 0; zdown_local = ZeroSetofComplex(); }
					if (i != 0) {
						left = *((f->fmatrix)+((f->xpixel*j)+(i-1)));
						zleft_local = *((f->zmatrix)+((f->xpixel*j)+(i-1)));
						compteur++;
					} else { left = 0; zleft_local = ZeroSetofComplex();}
					if (i < (f->xpixel-1)) {
						right = *((f->fmatrix)+((f->xpixel*j)+(i+1)));
						zright_local = *((f->zmatrix)+((f->xpixel*j)+(i+1)));
						compteur++;
					} else { right = 0; zright_local = ZeroSetofComplex();}

					if ((up == right) && (right == down) && (down == left) && (left == up)) {
						// Divergence identique : pas de calcul
						*((f->fmatrix)+((f->xpixel*j)+i)) = up;

						double rz, iz;
						rz = (Rez(zup_local)+Rez(zdown_local)+Rez(zleft_local)+Rez(zright_local))/compteur;
						iz = (Imz(zup_local)+Imz(zdown_local)+Imz(zleft_local)+Imz(zright_local))/compteur;
						*((f->zmatrix)+((f->xpixel*j)+i)) = MakeComplex(rz, iz);
						if (f->zmatrix_gmp != NULL) {
							complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+i], *((f->zmatrix)+((f->xpixel*j)+i)), prec);
						}

					} else {
						// Calcul nécessaire
						mpf_set_ui(i_mpf, i);
						mpf_set_ui(j_mpf, j);

						mpf_mul(xg, step_x_local, i_mpf);
						mpf_add(xg, xg, xmin_local);

						mpf_mul(yg, step_y_local, j_mpf);
						mpf_add(yg, yg, ymin_local);

						mpf_set(zPixel_gmp.x, xg);
						mpf_set(zPixel_gmp.y, yg);
						zPixel_gmp.is_nan = 0;

						result_local = FormulaSelector_GMP(f_local, zPixel_gmp);

						*((f->fmatrix)+((f->xpixel*j)+i)) = result_local.iteration;
						*((f->zmatrix)+((f->xpixel*j)+i)) = result_local.z;
						if (f->zmatrix_gmp != NULL) {
							complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+i], result_local.z, prec);
						}
					}
				}
			}

			#pragma omp barrier

			#pragma omp master
			{
				// Mise à jour de la progression (Pass 2 = deuxième moitié)
				if (guiPtr != NULL && progress != NULL) {
					int percent = progressMid +
						((chunk_start + chunk_height_p2) * (progressEnd - progressMid)) / f->ypixel;
					if (percent > progressEnd) percent = progressEnd;
					SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
					*progress = percent;
				}

				// Vérifier annulation utilisateur
				if (check_events_and_cancel()) {
					printf("Calcul GMP DDp2 Pass 2 annulé par l'utilisateur\n");
					cancelled_p2 = 1;
				}
			}

			#pragma omp barrier  // Assure que tous les threads voient cancelled_p2 mis à jour
		}

		// Libération des variables thread-local
		complex_gmp_clear(&zPixel_gmp);
		mpf_clear(xg);
		mpf_clear(yg);
		mpf_clear(i_mpf);
		mpf_clear(j_mpf);
		mpf_clear(step_x_local);
		mpf_clear(step_y_local);
		mpf_clear(xmin_local);
		mpf_clear(ymin_local);
		gmp_mul_temps_clear(&f_local.mul_temps);
		gmp_iteration_context_clear(&f_local.iteration_ctx);
		// Libérer les coordonnées GMP thread-local
		mpf_clear(f_local.xmin_gmp);
		mpf_clear(f_local.xmax_gmp);
		mpf_clear(f_local.ymin_gmp);
		mpf_clear(f_local.ymax_gmp);
	} // Fin de la région parallèle pour Pass 2
#else
	// Version séquentielle Pass 2
	{
		complex_gmp zPixel_gmp;
		mpf_t xg, yg, i_mpf, j_mpf, step_x_local, step_y_local, xmin_local, ymin_local;
		fractalresult result;

		complex_gmp_init(&zPixel_gmp, prec);
		mpf_init2(xg, prec);
		mpf_init2(yg, prec);
		mpf_init2(i_mpf, prec);
		mpf_init2(j_mpf, prec);
		mpf_init2(step_x_local, prec);
		mpf_init2(step_y_local, prec);
		mpf_init2(xmin_local, prec);
		mpf_init2(ymin_local, prec);

		mpf_set_d(step_x_local, step_x_double);
		mpf_set_d(step_y_local, step_y_double);
		mpf_set_d(xmin_local, xmin_double);
		mpf_set_d(ymin_local, ymin_double);

		for (j=0; j<f->ypixel; j++) {
			if (fmod(j,2) == 0) {starti=1;} else {starti=0;}
			for (i=starti; i<f->xpixel; i=i+2) {
				// Vérifier les bords des matrices pour éviter Seg Fault
				compteur=0;
				if (j != 0) {
					up = *((f->fmatrix)+((f->xpixel*(j-1))+i));
					zup = *((f->zmatrix)+((f->xpixel*(j-1))+i));
					compteur++;
				} else { up = 0; zup = ZeroSetofComplex();}
				if (j < (f->ypixel-1)) {
					down = *((f->fmatrix)+((f->xpixel*(j+1))+i));
					zdown = *((f->zmatrix)+((f->xpixel*(j+1))+i));
					compteur++;
				} else { down = 0; zdown = ZeroSetofComplex(); }
				if (i != 0) {
					left = *((f->fmatrix)+((f->xpixel*j)+(i-1)));
					zleft = *((f->zmatrix)+((f->xpixel*j)+(i-1)));
					compteur++;
				} else { left = 0; zleft = ZeroSetofComplex();}
				if (i < (f->xpixel-1)) {
					right = *((f->fmatrix)+((f->xpixel*j)+(i+1)));
					zright = *((f->zmatrix)+((f->xpixel*j)+(i+1)));
					compteur++;
				} else { right = 0; zright = ZeroSetofComplex();}

				if ((up == right) && (right == down) && (down == left) && (left == up)) {
					// Divergence identique : pas de calcul
					*((f->fmatrix)+((f->xpixel*j)+i)) = up;

					double rz, iz;
					rz = (Rez(zup)+Rez(zdown)+Rez(zleft)+Rez(zright))/compteur;
					iz = (Imz(zup)+Imz(zdown)+Imz(zleft)+Imz(zright))/compteur;
					*((f->zmatrix)+((f->xpixel*j)+i)) = MakeComplex(rz, iz);
					if (f->zmatrix_gmp != NULL) {
						complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+i], *((f->zmatrix)+((f->xpixel*j)+i)), prec);
					}

				} else {
					// Calcul nécessaire
					mpf_set_ui(i_mpf, i);
					mpf_set_ui(j_mpf, j);

					mpf_mul(xg, step_x_local, i_mpf);
					mpf_add(xg, xg, xmin_local);

					mpf_mul(yg, step_y_local, j_mpf);
					mpf_add(yg, yg, ymin_local);

					mpf_set(zPixel_gmp.x, xg);
					mpf_set(zPixel_gmp.y, yg);
					zPixel_gmp.is_nan = 0;
					result = FormulaSelector_GMP(*f, zPixel_gmp);

					*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
					*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
					if (f->zmatrix_gmp != NULL) {
						complex_to_gmp_inplace(&f->zmatrix_gmp[(f->xpixel*j)+i], result.z, prec);
					}
				}
			}
		}

		mpf_clear(xg);
		mpf_clear(yg);
		mpf_clear(i_mpf);
		mpf_clear(j_mpf);
		mpf_clear(step_x_local);
		mpf_clear(step_y_local);
		mpf_clear(xmin_local);
		mpf_clear(ymin_local);
		complex_gmp_clear(&zPixel_gmp);
	}
#endif

	// Mise à jour finale de la progression
	if (guiPtr != NULL && progress != NULL) {
		*progress = progressEnd;
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressEnd, fractalName);
	}
}
#endif

void Fractal_CalculateMatrix_DDp2 (fractal* f, SDL_Surface* canvas, void* guiPtr, int* progress, int progressStart, int progressEnd, const char* fractalName) {
	// Invalider cmatrix car DDp2 modifie aussi fmatrix
	f->cmatrix_valid = 0;

#ifdef HAVE_GMP
	if (f->use_gmp) {
		Fractal_CalculateMatrix_DDp2_GMP(f, canvas, guiPtr, progress, progressStart, progressEnd, fractalName);
		return;
	}
#endif
	int j;
	int xpixel = f->xpixel;
	int ypixel = f->ypixel;
	double xmin = f->xmin;
	double ymin = f->ymin;
	double xstep = (f->xmax - f->xmin) / xpixel;
	double ystep = (f->ymax - f->ymin) / ypixel;
	int progressMid = progressStart + (progressEnd - progressStart) / 2;

	printf ("Calculating EscapeTimeFractal Matrix with DDp2 ...\n");
#ifdef HAVE_OPENMP
	int num_threads = omp_get_max_threads();
	printf("Using OpenMP with %d threads\n", num_threads);
	// Granularité adaptative
	int chunk_size = (ypixel / 2) / (num_threads * 4);
	if (chunk_size < 1) chunk_size = 1;
	if (chunk_size > 64) chunk_size = 64;
#else
	int num_threads = 1;
#endif

	// Mise à jour initiale de la progression
	if (guiPtr != NULL && progress != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressStart, fractalName);
		*progress = progressStart;
	}

	// Pass 1: Calcul initial des pixels impairs (parallélisable)
#ifdef HAVE_OPENMP
	#pragma omp parallel
	{
		complex zPixel_local;
		fractalresult result_local;
		
		#pragma omp for schedule(guided) nowait
		for (j = 1; j < ypixel; j += 2) {
#else
	for (j = 1; j < ypixel; j += 2) {
		complex zPixel_local;
		fractalresult result_local;
#endif
		int i;
		for (i = 1; i < xpixel; i += 2) {
			double xg = xstep * i + xmin;
			double yg = ystep * j + ymin;
			zPixel_local = MakeComplex(xg, yg);

			result_local = FormulaSelector(*f, zPixel_local);
			f->fmatrix[xpixel * j + i] = result_local.iteration;
			f->zmatrix[xpixel * j + i] = result_local.z;
		}
	}
#ifdef HAVE_OPENMP
	} // Fin de la région parallèle pour Pass 1
#endif

	// Mise à jour de la progression après Pass 1
	if (guiPtr != NULL && progress != NULL) {
		*progress = progressMid;
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressMid, fractalName);
	}

	// Pass 2: Divergence Detection (parallélisable - chaque pixel écrit à son propre indice)
#ifdef HAVE_OPENMP
	#pragma omp parallel
	{
		complex zPixel_local;
		fractalresult result_local;
		complex zup_local, zdown_local, zleft_local, zright_local;
		
		#pragma omp for schedule(guided) nowait
		for (j = 0; j < ypixel; j++) {
#else
	for (j = 0; j < ypixel; j++) {
		complex zPixel_local;
		fractalresult result_local;
		complex zup_local, zdown_local, zleft_local, zright_local;
#endif
		int starti = (j % 2 == 0) ? 1 : 0;  // Lignes paires et impaires
		int i;
		for (i = starti; i < xpixel; i += 2) {
			int compteur = 0;
			int up, down, left, right;
			complex zup, zdown, zleft, zright;

			// Verifions les bords des matrices pour eviter Seg Fault
			if (j != 0) {
				up = f->fmatrix[xpixel * (j - 1) + i];
				zup = f->zmatrix[xpixel * (j - 1) + i];
				compteur++;
			} else {
				up = 0;
				zup = ZeroSetofComplex();
			}
			if (j < (ypixel - 1)) {
				down = f->fmatrix[xpixel * (j + 1) + i];
				zdown = f->zmatrix[xpixel * (j + 1) + i];
				compteur++;
			} else {
				down = 0;
				zdown = ZeroSetofComplex();
			}
			if (i != 0) {
				left = f->fmatrix[xpixel * j + (i - 1)];
				zleft = f->zmatrix[xpixel * j + (i - 1)];
				compteur++;
			} else {
				left = 0;
				zleft = ZeroSetofComplex();
			}
			if (i < (xpixel - 1)) {
				right = f->fmatrix[xpixel * j + (i + 1)];
				zright = f->zmatrix[xpixel * j + (i + 1)];
				compteur++;
			} else {
				right = 0;
				zright = ZeroSetofComplex();
			}

			/* Si les points adjacents ont meme valeur de divergence alors */
			if ((up == right) && (right == down) && (down == left) && (left == up)) {
				// OK, on ne calcule pas les Iterations
				f->fmatrix[xpixel * j + i] = up;

				// Et on calcule une moyenne pour la valeur de z
				double rz = (Rez(zup) + Rez(zdown) + Rez(zleft) + Rez(zright)) / compteur;
				double iz = (Imz(zup) + Imz(zdown) + Imz(zleft) + Imz(zright)) / compteur;
				f->zmatrix[xpixel * j + i] = MakeComplex(rz, iz);
			} else {
				// Sinon, je calcule
				double xg = xstep * i + xmin;
				double yg = ystep * j + ymin;
				zPixel_local = MakeComplex(xg, yg);

				result_local = FormulaSelector(*f, zPixel_local);
				f->fmatrix[xpixel * j + i] = result_local.iteration;
				f->zmatrix[xpixel * j + i] = result_local.z;
			}
		}
	}
#ifdef HAVE_OPENMP
	} // Fin de la région parallèle pour Pass 2
#endif

	// Mise à jour finale de la progression
	if (guiPtr != NULL && progress != NULL) {
		*progress = progressEnd;
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressEnd, fractalName);
	}
}



// Color Formulae Functions
// *************************

/* Un mode de couleur classique */

// Un mode de couleur bizarre pour tester que tout fonctionne avec les complexes
void FractalColorTest (fractal* f) {
	int i, j;
	int iteration, greyvalue;
	color c;

	for (j = 0; j < f->ypixel; j++) {
		for (i=0; i < f->xpixel; i++) {
			iteration = *((f->fmatrix)+((f->xpixel*j)+i));
			greyvalue = 255 - (iteration*255)/f->iterationMax;
			c.g = (int) fmod (greyvalue/(Rez (*((f->zmatrix)+((f->xpixel*j)+i)))), 255);
			c.r = (int) fmod (greyvalue, 255);
			c.b = (int) fmod (greyvalue/(Imz (*((f->zmatrix)+((f->xpixel*j)+i)))), 255);
			*((f->cmatrix)+((f->xpixel*j)+i)) = c;
		}
	}
}


/* Calcul smooth iteration pour coloring continu */
double Fractal_SmoothIteration(fractal* f, int i, int j) {
	int iteration = *((f->fmatrix)+((f->xpixel*j)+i));
	
	// Pour le type Lyapunov (17), retourner directement la valeur normalisée
	// qui est déjà stockée dans fmatrix comme "nombre d'itérations"
	if (f->type == 17) {
		return (double)iteration / f->iterationMax;
	}
	
	complex z = *((f->zmatrix)+((f->xpixel*j)+i));
	double mag = Magz(z);

	if (iteration >= f->iterationMax || mag < 1.0) {
		return (double)iteration / f->iterationMax;
	}

	// Formule smooth coloring
	double log_zn = log(mag) / 2.0;
	double nu = log(log_zn / log(2.0)) / log(2.0);
	double smooth = iteration + 1 - nu;

	return smooth / f->iterationMax;
}

/* Conversion HSV vers RGB */
color HSVtoRGB(double h, double s, double v) {
	color c;
	double r, g, b;
	int i;
	double f, p, q, t;

	if (s == 0) {
		r = g = b = v;
	} else {
		h = fmod(h, 360.0);
		if (h < 0) h += 360.0;
		h /= 60.0;
		i = (int)floor(h);
		f = h - i;
		p = v * (1 - s);
		q = v * (1 - s * f);
		t = v * (1 - s * (1 - f));

		switch (i) {
			case 0: r = v; g = t; b = p; break;
			case 1: r = q; g = v; b = p; break;
			case 2: r = p; g = v; b = t; break;
			case 3: r = p; g = q; b = v; break;
			case 4: r = t; g = p; b = v; break;
			default: r = v; g = p; b = q; break;
		}
	}

	c.r = (int)(r * 255);
	c.g = (int)(g * 255);
	c.b = (int)(b * 255);
	c.a = 255;
	return c;
}

/* Selection du mode de couleur - using unified colorization system */
//**********************************

void Fractal_CalculateColorMatrix (fractal* f, SDL_Surface* canvas, void* guiPtr, int* progress, int progressStart, int progressEnd) {
	const char* fractalName = Fractal_GetTypeName(f->type);
	const gradient_table* gradient;

	// Réutilisation de cmatrix si déjà valide pour le colorMode et colorRepeat actuel
	if (f->cmatrix_valid && f->last_colorMode == f->colorMode && f->last_colorRepeat == f->colorRepeat) {
		// Sauter le calcul, mettre à jour la progression directement
		if (guiPtr != NULL && progress != NULL) {
			*progress = progressEnd;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressEnd, fractalName);
		}
		return;
	}

	// Use the unified colorization system
	gradient = Colorization_GetPalette(f->colorMode);
	Fractal_ApplyGradient(f, gradient);

	// Marquer cmatrix comme valide
	f->cmatrix_valid = 1;
	f->last_colorMode = f->colorMode;
	f->last_colorRepeat = f->colorRepeat;

	// Mise à jour de la progression après colorisation
	if (guiPtr != NULL && progress != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressEnd, fractalName);
		*progress = progressEnd;
	}
}

color Fractal_ReadColorMatrix (fractal f, int i, int j) {
	color c;
	c = *((f.cmatrix)+((f.xpixel*j)+i));
	return c;
}
int Fractal_ReadColorMatrixRed (fractal f, int i, int j) {
	return (*((f.cmatrix)+((f.xpixel*j)+i))).r;
}
int Fractal_ReadColorMatrixBlue (fractal f, int i, int j) {
	return (*((f.cmatrix)+((f.xpixel*j)+i))).b;
}
int Fractal_ReadColorMatrixGreen (fractal f, int i, int j) {
	return (*((f.cmatrix)+((f.xpixel*j)+i))).g;
}



// Fonctions de gestion du cache
static void Fractal_SaveCache(fractal* f) {
	// Sauvegarder l'état actuel dans le cache
	if (f->cache.fmatrix_cached == NULL || 
	    f->cache.cache_xpixel != f->xpixel || 
	    f->cache.cache_ypixel != f->ypixel) {
		// Réallouer si nécessaire
		if (f->cache.fmatrix_cached != NULL) {
			free(f->cache.fmatrix_cached);
			free(f->cache.cmatrix_cached);
		}
		f->cache.fmatrix_cached = (int*) malloc(f->xpixel * f->ypixel * sizeof(int));
		f->cache.cmatrix_cached = (color*) malloc(f->xpixel * f->ypixel * sizeof(color));
		if (f->cache.fmatrix_cached == NULL || f->cache.cmatrix_cached == NULL) {
			fprintf(stderr, "Erreur allocation mémoire cache\n");
			return;
		}
		f->cache.cache_xpixel = f->xpixel;
		f->cache.cache_ypixel = f->ypixel;
	}
	
	// Copier les matrices
	memcpy(f->cache.fmatrix_cached, f->fmatrix, f->xpixel * f->ypixel * sizeof(int));
	memcpy(f->cache.cmatrix_cached, f->cmatrix, f->xpixel * f->ypixel * sizeof(color));
	
	// Sauvegarder les coordonnées
	f->cache.xmin_cached = f->xmin;
	f->cache.xmax_cached = f->xmax;
	f->cache.ymin_cached = f->ymin;
	f->cache.ymax_cached = f->ymax;
	f->cache.cache_valid = 1;
}

static int Fractal_CanReuseCache(fractal* f) {
	// Vérifier si on peut réutiliser le cache
	if (!f->cache.cache_valid || f->cache.fmatrix_cached == NULL) {
		return 0;
	}
	
	// Calculer le facteur de zoom
	double old_range_x = f->cache.xmax_cached - f->cache.xmin_cached;
	double old_range_y = f->cache.ymax_cached - f->cache.ymin_cached;
	double new_range_x = f->xmax - f->xmin;
	double new_range_y = f->ymax - f->ymin;
	
	double zoom_factor = old_range_x / new_range_x;
	
	// Réutiliser si zoom < 2× et même dimensions
	if (zoom_factor < 2.0 && zoom_factor > 0.5 && 
	    f->cache.cache_xpixel == f->xpixel && 
	    f->cache.cache_ypixel == f->ypixel) {
		return 1;
	}
	
	return 0;
}

Uint32 Fractal_Draw (SDL_Surface *canvas, fractal myfractal,int decalageX,int decalageY, void* guiPtr) {

	// Les fractales spéciales (Buddhabrot, Lyapunov) utilisent leur propre algorithme de rendu
	if (myfractal.type == 16) {
		return Buddhabrot_Draw(canvas, &myfractal, decalageX, decalageY, guiPtr);
	}
	if (myfractal.type == 17) {
		return Lyapunov_Draw(canvas, &myfractal, decalageX, decalageY, guiPtr);
	}

	int i, j;
	Uint8 r, g, b;
	Uint32 time; // Test de temps de calcul en ms
	int progress = 0;
	const char* fractalName = Fractal_GetTypeName(myfractal.type);

	time = SDL_GetTicks();

	// DDp1 : 0-25%
	Fractal_CalculateMatrix_DDp1 (&myfractal, canvas, guiPtr, &progress, 0, 25, fractalName);
	
	// Colorisation après DDp1 : 25-37.5%
	Fractal_CalculateColorMatrix (&myfractal, canvas, guiPtr, &progress, 25, 37);

	//Draw Demi-Fractal to the SDL_Buffer and colorize it
	for (i=0; i<myfractal.xpixel; i++) {
		for (j=0; j<myfractal.ypixel; j++) {
			r=Fractal_ReadColorMatrixRed (myfractal,i,j);
			g=Fractal_ReadColorMatrixGreen (myfractal,i,j);
			b=Fractal_ReadColorMatrixBlue (myfractal,i,j);
			pixelRGBA(canvas, (Sint16) (i+decalageX), (Sint16) (j+decalageY),  r, g,  b, 255);
		}
	}
	SDL_UpdateRect (canvas, 0, 0, canvas->w, canvas->h);
	
	// Mise à jour progression après affichage DDp1 : 37-37.5%
	if (guiPtr != NULL) {
		progress = 37;
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progress, fractalName);
	}

	// DDp2 : 37.5-87.5%
	Fractal_CalculateMatrix_DDp2 (&myfractal, canvas, guiPtr, &progress, 37, 87, fractalName);
	
	// Colorisation après DDp2 : 87.5-100%
	Fractal_CalculateColorMatrix (&myfractal, canvas, guiPtr, &progress, 87, 100);

	//Draw last Demi-Fractal to the SDL_Buffer and colorize it
	for (i=0; i<myfractal.xpixel; i++) {
		for (j=0; j<myfractal.ypixel; j++) {
			r=Fractal_ReadColorMatrixRed (myfractal,i,j);
			g=Fractal_ReadColorMatrixGreen (myfractal,i,j);
			b=Fractal_ReadColorMatrixBlue (myfractal,i,j);
			pixelRGBA(canvas, (Sint16) (i+decalageX), (Sint16) (j+decalageY),  r,g,  b, 255);
		}
	}
	SDL_UpdateRect (canvas, 0, 0, canvas->w, canvas->h);
	
	// Mise à jour progression finale : 100%
	if (guiPtr != NULL) {
		progress = 100;
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progress, fractalName);
	}
	
	// Sauvegarder le cache pour réutilisation lors de prochains zooms
	Fractal_SaveCache(&myfractal);
	
	time = SDL_GetTicks() - time;
	printf ("Time Elapsed : %d ms\n", time);
	return time;
}

void Fractal_ChangeType (fractal* f, int type) {
#ifdef HAVE_GMP
	// Libérer la mémoire GMP existante avant de changer de type
	if (f->zmatrix_gmp != NULL) {
		for (int i = 0; i < f->xpixel * f->ypixel; i++) {
			complex_gmp_clear(&f->zmatrix_gmp[i]);
		}
		free(f->zmatrix_gmp);
		f->zmatrix_gmp = NULL;
	}
	f->use_gmp = 0;
	f->gmp_precision = 64;
#endif

	f->type = type;
	printf ("On initialise une EscapeTimeFractal de type %d\n",type);

	// Définir la répétition par défaut selon le type
	if (type == 17) {
		f->colorRepeat = 2;  // 2 répétitions pour Lyapunov
	} else {
		f->colorRepeat = 20;  // 20 répétitions pour escape-time
	}

	switch (type)
	{
	case 3:
	Mendelbrot_def (f);
	break;
	case 4:
	Julia_def (f);
	break;
	case 5:
	JuliaSin_def (f);
	break;
	case 6:
	Newton_def (f);
	break;
	case 7:
	Phoenix_def (f);
	break;
	case 8:
	Buffalo_def (f);
	break;
	case 9:
	Barnsley1j_def (f);
	break;
	case 10:
	Barnsley1m_def (f);
	break;
	case 11:
	Magnet1j_def (f);
	break;
	case 12:
	Magnet1m_def (f);
	break;
	case 13:
	BurningShip_def (f);
	break;
	case 14:
	Tricorn_def (f);
	break;
	case 15:
	Mandelbulb_def (f);
	break;
	case 16:
	Buddhabrot_def (f);
	break;
	case 17:
	Lyapunov_def (f);
	break;
	case 18:
	PerpendicularBurningShip_def (f);
	break;
	case 19:
	Celtic_def (f);
	break;
	case 20:
	AlphaMandelbrot_def (f);
	break;
	case 21:
	PickoverStalks_def (f);
	break;
	case 22:
	Nova_def (f);
	break;
	case 23:
	Multibrot_def (f);
	break;
	default:
	Mendelbrot_def (f);
	}

#ifdef HAVE_GMP
	// Mettre à jour la précision après le changement de type
	precision_update_fractal(f);
	precision_update_gmp_structures(f);
	
	// Synchroniser les coordonnées double → GMP
	// Les coordonnées GMP peuvent ne pas être initialisées lors du premier appel
	// depuis Fractal_Init. On les initialise toujours ici pour éviter les segfaults.
	// mpf_init2 peut être appelé plusieurs fois sans problème
	if (f->use_gmp) {
		mpf_init2(f->xmin_gmp, f->gmp_precision);
		mpf_init2(f->xmax_gmp, f->gmp_precision);
		mpf_init2(f->ymin_gmp, f->gmp_precision);
		mpf_init2(f->ymax_gmp, f->gmp_precision);
		mpf_set_d(f->xmin_gmp, f->xmin);
		mpf_set_d(f->xmax_gmp, f->xmax);
		mpf_set_d(f->ymin_gmp, f->ymin);
		mpf_set_d(f->ymax_gmp, f->ymax);
	} else {
		// Initialiser avec précision minimale même si use_gmp est 0
		mpf_init2(f->xmin_gmp, 64);
		mpf_init2(f->xmax_gmp, 64);
		mpf_init2(f->ymin_gmp, 64);
		mpf_init2(f->ymax_gmp, 64);
		mpf_set_d(f->xmin_gmp, f->xmin);
		mpf_set_d(f->xmax_gmp, f->xmax);
		mpf_set_d(f->ymin_gmp, f->ymin);
		mpf_set_d(f->ymax_gmp, f->ymax);
	}
	
	// Réallouer zmatrix_gmp si nécessaire
	if (f->use_gmp) {
		if (f->zmatrix_gmp != NULL) {
			for (int i = 0; i < f->xpixel * f->ypixel; i++) {
				complex_gmp_clear(&f->zmatrix_gmp[i]);
			}
			free(f->zmatrix_gmp);
		}
		f->zmatrix_gmp = (complex_gmp *) malloc ((f->xpixel*f->ypixel)* sizeof (complex_gmp));
		if (f->zmatrix_gmp == NULL) {
			fprintf(stderr, "Erreur allocation mémoire GMP fractale\n");
			exit(1);
		}
		// Initialiser tous les éléments GMP
		for (int i = 0; i < f->xpixel * f->ypixel; i++) {
			complex_gmp_init(&f->zmatrix_gmp[i], f->gmp_precision);
		}
	}
#endif
}

// *************************************************************************
// *	Fonction utilitaire pour obtenir le nom de la fractale
// *************************************************************************

const char* Fractal_GetTypeName(int type) {
	const char* typeNames[] = {
		"",           // 0
		"Von Koch",   // 1
		"Dragon",     // 2
		"Mandelbrot", // 3
		"Julia",      // 4
		"Julia Sin",  // 5
		"Newton",     // 6
		"Phoenix",    // 7
		"Buffalo",    // 8
		"Barnsley J", // 9
		"Barnsley M", // 10
		"Magnet J",   // 11
		"Magnet M",   // 12
		"Burning Ship", // 13
		"Tricorn",    // 14
		"Mandelbulb", // 15
		"Buddhabrot", // 16
		"Lyapunov Zircon City",   // 17
		"Perpendicular Burning Ship", // 18
		"Celtic",     // 19
		"Alpha Mandelbrot",  // 20
		"Pickover Stalks",   // 21
		"Nova",              // 22
		"Multibrot"          // 23
	};

	if (type >= 0 && type <= 23) {
		return typeNames[type];
	}
	return "Unknown";
}

// *************************************************************************
// *	Voici les formules et les valeurs par default de quelques fractales
// *************************************************************************


// Calcul de la mendelbrot ...
fractalresult Mendelbrot_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;

	z = f.seed;
	do { //z(n+1) = z(n)^2 + pixel
		complex zTemp;
		i++;
		zTemp = Mulz (z,z); // z(n)^2
		z = Addz (zTemp, zPixel); //Ajout de pixel
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}
void Mendelbrot_def (fractal* f) {
	/* Valeurs de base pour la mandelbrot */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex ();
	f->bailout= 4;
	f->iterationMax= 9370;
	f->zoomfactor = 8;
}



// Calcul de la Julia
fractalresult Julia_Iteration ( fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;
	z = zPixel;
	do { //z(n+1) = z(n)^2 + seed
		complex zTemp;
		i++;
		zTemp = Mulz (z,z); // z(n)^2
		z = Addz (zTemp, f.seed); // ajout de seed
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}
void Julia_def (fractal* f) {
	/* Valeurs de base pour la julia  */
	f->xmin = -2.0;
	f->xmax = +2.0;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = MakeComplex (0.36228,- 0.0777);
	f->bailout= 4;
	f->iterationMax= 6250;
	f->zoomfactor = 4;
}




/* Calcul de JuliaSin */
fractalresult JuliaSin_Iteration ( fractal f, complex zPixel) {
	/* JuliaSin
	 * z(0) = pixel;
         * z(k+1) = c * sin(z(k))
	 */
	int i=0;
	complex z;
	fractalresult result;
	z = zPixel;
	do {
		i++;
		z = Mulz(f.seed, sinz (z));
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}
void JuliaSin_def (fractal* f) {
// Valeurs de base pour la JuliaSin
	f->xmin = -M_PI;
	f->xmax = M_PI;
	f->ymin = -2;
	f->ymax = 2;
	f->seed = MakeComplex (1,0.1);
	f->bailout= 4;
	f->iterationMax= 6250;
	f->zoomfactor = 4;
}




/* Calcul de la newton */
fractalresult Newton_Iteration ( fractal f, complex zPixel) {
	/*
	 * z(0) = pixel;
         * z(n+1) = ((p-1)*z(n)^p + 1)/(p*z(n)^(p - 1)).
	*/
	int i=0;
	complex z;
	fractalresult result;
	z = zPixel;
	do {
		int p = Rez(f.seed);  // Degree polynomial
		complex zQuot, zNum;
		i++;
		zNum = Addz ((Mulz(MakeComplex(p-1, 0.0),Powz(z,MakeComplex(p, 0.0)))),MakeComplex(1.0, 0.0));
		zQuot = Mulz (Powz (z, MakeComplex (p-1,0.0)),MakeComplex(p, 0.0));
		z = Divz (zNum, zQuot);
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}
void Newton_def (fractal* f) {
	/* Valeurs de base pour la newton */
	f->seed = MakeComplex (8,0);
	f->xmin = -3;
	f->xmax = 3;
	f->ymin = -2;
	f->ymax = 2;
	f->bailout= 4;
	f->iterationMax= 1000;
	f->zoomfactor = 4;
}




/* Calcul de la phoenix pour degree 0 uniquement !!! */
fractalresult Phoenix_Iteration ( fractal f, complex zPixel) {
	/*
	 * z(0) = pixel, y(0) = 0;
         * For degree = 0: z(n+1) = z(n)^2 + p1.x + p2.x*y(n), y(n+1) = z(n)
         * For degree >= 2:
         * z(n+1) = z(n)^degree + p1.x*z(n)^(degree-1) + p2.x*y(n), y(n+1) = z(n)
         * For degree <= -3:
         * z(n+1) = z(n)^|degree| + p1.x*z(n)^(|degree|-2) + p2.x*y(n), y(n+1) = z(n)
	*/
	int i=0;
	complex z,y;
	complex zTemp,zpTemp;
	fractalresult result;
	int degree = 0;
	double p1 = 0.56667;
	double p2 = -0.5;
	
	z = zPixel;
	y = ZeroSetofComplex();
	do {
		i++;
		if (degree == 0) {
			zTemp = Mulz(z,z);
			zTemp = MakeComplex(Rez(zTemp)+p1,Imz(zTemp));
			zpTemp = ScalarTimesofComplex(p2,y);
			zTemp = Addz (zTemp,zpTemp);
		y = z;
		z = zTemp;
		}
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}
void Phoenix_def (fractal* f) {
	/* Valeurs de base pour la phoenix 	*/
	f->xmin = -2.0;
	f->xmax = +2.0;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->bailout= 4;
	f->iterationMax= 9370;
	f->zoomfactor = 8;
}




/* Sierpinski supprimé - type 8 */

/* Calcul de Barnsley1j  */
fractalresult Barnsley1j_Iteration ( fractal f, complex zPixel) {
	/* Barnsleyj1
	 * z(0) = pixel;
         * if real(z) >= 0
         * z(n+1) = (z-1)*c
         * else
         * z(n+1) = (z+1)*c */
	int i=0;
	complex z, zTemp;
	fractalresult result;
	z = zPixel;
	do {
		i++;
		if (Rez(z) >= 0) {
		zTemp = MakeComplex((Rez (z) - 1),Imz(z));
		z = Mulz (zTemp, f.seed);
		} else {
		zTemp = MakeComplex((Rez (z) + 1),Imz(z));
		z = Mulz (zTemp, f.seed);
		}
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}
void Barnsley1j_def (fractal* f) {
// Valeurs de base pour la Barnsleyj1
	f->xmin = -4;
	f->xmax = 4;
	f->ymin = -3;
	f->ymax = 3;
	f->seed = MakeComplex (1.1,0.6);
	f->bailout= 4;
	f->iterationMax= 3120;
	f->zoomfactor = 4;
}


/* Calcul de Barnsley1m */
fractalresult Barnsley1m_Iteration ( fractal f, complex zPixel) {
	int i=0;
	complex z,c, zTemp;
	fractalresult result;
	z = zPixel;
	c= zPixel;
	do {
		i++;
		if (Rez(z) >= 0) {
		zTemp = MakeComplex((Rez (z) - 1),Imz(z));
		z = Mulz (zTemp, c);
		} else {
		zTemp = MakeComplex((Rez (z) + 1),Imz(z));
		z = Mulz (zTemp, c);
		}
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}
void Barnsley1m_def (fractal* f) {
// Valeurs de base pour la Barnsley1m
	f->xmin = -3;
	f->xmax = 3;
	f->ymin = -2;
	f->ymax = 2;
	f->bailout= 4;
	f->iterationMax= 9370;
	f->zoomfactor = 8;
}



/* Calcul de Magnet1j  */
fractalresult Magnet1j_Iteration ( fractal f, complex zPixel) {
	int i=0;
	complex z,N,Q;
	fractalresult result;
	z = zPixel;
	do {
		i++;
		N= Addz(Mulz(z,z),MakeComplex (Rez(f.seed)-1,Imz(f.seed)));
		N= Mulz(N,N);
		Q= Addz(Mulz(MakeComplex(2,0),z),MakeComplex (Rez(f.seed)-2,Imz(f.seed)));
		z = Divz(N,Q);
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}
void Magnet1j_def (fractal* f) {
// Valeurs de base pour la Magnet1j
	f->seed = MakeComplex ( 1.625458, -0.306159);
	f->xmin = -2;
	f->xmax = 2;
	f->ymin = -2;
	f->ymax = 2;
	f->bailout= 4;
	f->iterationMax= 9370;
	f->zoomfactor = 4;
}


/* Calcul de Magnet1m  */
fractalresult Magnet1m_Iteration ( fractal f, complex zPixel) {
	int i=0;
	complex z,c,N,Q;
	fractalresult result;
	c = zPixel;
	z= ZeroSetofComplex();
	do {
		i++;
		N= Addz(Mulz(z,z),MakeComplex (Rez(c)-1,Imz(c)));
		N= Mulz(N,N);
		Q= Addz(Mulz(MakeComplex(2,0),z),MakeComplex (Rez(c)-2,Imz(c)));
		z = Divz(N,Q);
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}
void Magnet1m_def (fractal* f) {
// Valeurs de base pour la Magnet1m
	f->xmin = -3;
	f->xmax = 2;
	f->ymin = -2;
	f->ymax = 2;
	f->bailout= 4;
	f->iterationMax= 9370;
	f->zoomfactor = 8;
}

/* Calcul de Burning Ship */
fractalresult BurningShip_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;
	
	z = f.seed;
	do {
		complex zTemp;
		double re, im;
		i++;
		// z(n+1) = (|Re(z)| + i|Im(z)|)² + c
		re = fabs(Rez(z));
		im = fabs(Imz(z));
		zTemp = MakeComplex(re, im);
		zTemp = Mulz(zTemp, zTemp); // (|Re(z)| + i|Im(z)|)²
		z = Addz(zTemp, zPixel); // Ajout de pixel
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	
	result.iteration = i;
	result.z = z;
	return result;
}
void BurningShip_def (fractal* f) {
	/* Valeurs de base pour Burning Ship */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -2.0;
	f->ymax = 2.0;
	f->seed = ZeroSetofComplex();
	f->bailout= 4;
	f->iterationMax= 9370;
	f->zoomfactor = 8;
}

/* Calcul de Buffalo
 * Formule: z(n+1) = abs(Re(z²)) + i*abs(Im(z²)) + c
 * C'est une variation du Burning Ship avec abs() sur les deux parties après le carré
 */
fractalresult Buffalo_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;

	z = f.seed;
	do {
		complex zSq;
		double re_sq, im_sq;
		i++;
		// z² = (a+bi)² = (a²-b²) + 2abi
		zSq = Mulz(z, z);
		// Appliquer abs() aux deux parties
		re_sq = fabs(Rez(zSq));
		im_sq = fabs(Imz(zSq));
		z = MakeComplex(re_sq + Rez(zPixel), im_sq + Imz(zPixel));
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));

	result.iteration = i;
	result.z = z;
	return result;
}
void Buffalo_def (fractal* f) {
	/* Valeurs de base pour Buffalo */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -2.0;
	f->ymax = 2.0;
	f->seed = ZeroSetofComplex();
	f->bailout= 4;
	f->iterationMax= 9370;
	f->zoomfactor = 8;
}

/* Calcul de Tricorn */
fractalresult Tricorn_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;
	
	z = f.seed;
	do {
		complex zTemp, zConj;
		i++;
		// z(n+1) = (conjugué(z))² + c
		zConj = MakeComplex(Rez(z), -Imz(z)); // Conjugué
		zTemp = Mulz(zConj, zConj); // (conjugué(z))²
		z = Addz(zTemp, zPixel); // Ajout de pixel
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	
	result.iteration = i;
	result.z = z;
	return result;
}
void Tricorn_def (fractal* f) {
	/* Valeurs de base pour Tricorn */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout= 4;
	f->iterationMax= 9370;
	f->zoomfactor = 8;
}

/* Calcul de Mandelbulb (version 2D avec puissance 8) */
fractalresult Mandelbulb_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;
	
	z = f.seed;
	do {
		complex zTemp;
		i++;
		// z(n+1) = z(n)^8 + c
		// Calcul de z^8 par multiplications successives
		zTemp = Mulz(z, z);      // z^2
		zTemp = Mulz(zTemp, zTemp); // z^4
		zTemp = Mulz(zTemp, zTemp); // z^8
		z = Addz(zTemp, zPixel); // z^8 + c
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	
	result.iteration = i;
	result.z = z;
	return result;
}
void Mandelbulb_def (fractal* f) {
	/* Valeurs de base pour Mandelbulb */
	f->xmin = -1.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout= 4;
	f->iterationMax= 9370;
	f->zoomfactor = 8;
}

/* Calcul de Perpendicular Burning Ship
 * Formule: z(n+1) = (Re(z) - i·|Im(z)|)² + c
 * En termes de parties réelle et imaginaire:
 * x_{n+1} = x_n² - |y_n|² + c_x
 * y_{n+1} = -2·x_n·|y_n| + c_y
 */
fractalresult PerpendicularBurningShip_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;
	
	z = f.seed;
	do {
		double x, y, y_abs;
		double x2, y2;
		i++;
		// z(n+1) = (Re(z) - i·|Im(z)|)² + c
		x = Rez(z);
		y = Imz(z);
		y_abs = fabs(y);
		
		// Calcul de (x - i·y_abs)² = x² - y_abs² - i·2·x·y_abs
		x2 = x * x;
		y2 = y_abs * y_abs;
		z = MakeComplex(x2 - y2 + Rez(zPixel), -2.0 * x * y_abs + Imz(zPixel));
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	
	result.iteration = i;
	result.z = z;
	return result;
}
void PerpendicularBurningShip_def (fractal* f) {
	/* Valeurs de base pour Perpendicular Burning Ship */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout= 4;
	f->iterationMax= 5000;
	f->zoomfactor = 8;
}

/* Calcul de Celtic Fractal
 * Formule: z(n+1) = |Re(z²)| + i·Im(z²) + c
 * En termes de parties réelle et imaginaire:
 * u = x² - y²        // partie réelle de z²
 * v = 2·x·y          // partie imaginaire de z²
 * x_{n+1} = |u| + c_x
 * y_{n+1} = v + c_y
 */
fractalresult Celtic_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;
	
	z = f.seed;
	do {
		double x, y;
		double u, v;  // u = partie réelle de z², v = partie imaginaire de z²
		i++;
		// Calcul de z² = (x + iy)² = (x² - y²) + i(2xy)
		x = Rez(z);
		y = Imz(z);
		u = x * x - y * y;  // partie réelle de z²
		v = 2.0 * x * y;    // partie imaginaire de z²
		
		// Application de abs() à la partie réelle
		z = MakeComplex(fabs(u) + Rez(zPixel), v + Imz(zPixel));
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	
	result.iteration = i;
	result.z = z;
	return result;
}
void Celtic_def (fractal* f) {
	/* Valeurs de base pour Celtic Fractal */
	f->xmin = -2.0;
	f->xmax = 1.0;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout= 4;
	f->iterationMax= 5000;
	f->zoomfactor = 8;
}

/* Calcul de Alpha Mandelbrot (Nested/Composite)
 * Formule: z(n+1) = z² + (z² + c)² + c
 * En termes de parties réelle et imaginaire:
 * m = z² + c
 * z_{n+1} = z² + m² + c
 */
fractalresult AlphaMandelbrot_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;
	
	z = f.seed;
	do {
		complex zSq, m, mSq;
		i++;
		// Calcul de z²
		zSq = Mulz(z, z);
		// Calcul de m = z² + c
		m = Addz(zSq, zPixel);
		// Calcul de m²
		mSq = Mulz(m, m);
		// z_{n+1} = z² + m² + c
		z = Addz(Addz(zSq, mSq), zPixel);
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	
	result.iteration = i;
	result.z = z;
	return result;
}
void AlphaMandelbrot_def (fractal* f) {
	/* Valeurs de base pour Alpha Mandelbrot */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout= 4;
	f->iterationMax= 2000;  // Divergence plus rapide, moins d'itérations nécessaires
	f->zoomfactor = 8;
}

/* Calcul de Pickover Stalks / Biomorphs
 * Formule de base: z(n+1) = z² + c (Mandelbrot)
 * Système d'orbit trap: distance_trap = min(|Re(z)|, |Im(z)|)
 * On stocke trap_min et l'utilisons pour la coloration
 * Pour la coloration, on utilise -log(trap_min) normalisé dans iteration
 */
fractalresult PickoverStalks_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;
	double trap_min = 1e10;  // Valeur initiale très grande
	double trap_distance;
	const double trap_divisor = 0.03;  // Ajuste l'épaisseur des stalks
	
	z = f.seed;
	do {
		double re_abs, im_abs;
		i++;
		// Itération Mandelbrot: z(n+1) = z² + c
		z = Addz(Mulz(z, z), zPixel);
		
		// Calcul de la distance au trap (cross: axes Re=0 et Im=0)
		re_abs = fabs(Rez(z));
		im_abs = fabs(Imz(z));
		trap_distance = (re_abs < im_abs) ? re_abs : im_abs;
		
		// Mise à jour de trap_min
		if (trap_distance < trap_min) {
			trap_min = trap_distance;
		}
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	
	// Stocker une valeur normalisée basée sur trap_min dans iteration
	// Utiliser -log(trap_min) pour créer l'effet de coloration souhaité
	// Normaliser pour que cela fonctionne avec le système de colorisation
	if (trap_min > 1e-10) {
		double log_trap = -log(trap_min / trap_divisor);
		// Normaliser pour que cela corresponde à un nombre d'itérations raisonnable
		// Utiliser une valeur entre 0 et iterationMax
		result.iteration = (int)(log_trap * 100.0);
		if (result.iteration < 0) result.iteration = 0;
		if (result.iteration >= f.iterationMax) result.iteration = f.iterationMax - 1;
	} else {
		// trap_min très petit, point très proche des axes
		result.iteration = f.iterationMax - 1;
	}
	
	result.z = z;
	return result;
}
void PickoverStalks_def (fractal* f) {
	/* Valeurs de base pour Pickover Stalks */
	f->xmin = -2.0;
	f->xmax = 1.0;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout= 100.0;  // Bailout plus élevé pour capturer plus de détails
	f->iterationMax= 1000;
	f->zoomfactor = 8;
}

/* Calcul de Nova Fractal
 * Formule: z(n+1) = z - a·(p(z)/p'(z)) + c
 * Pour p(z) = z^p - 1: p'(z) = p·z^(p-1)
 * Donc: z(n+1) = z - a·((z^p - 1)/(p·z^(p-1))) + c
 * Pour p=3 (cubique): z(n+1) = z - a·((z³ - 1)/(3·z²)) + c
 */
fractalresult Nova_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z, z_prev;
	fractalresult result;
	complex a_relax = MakeComplex(1.0, 0.0);  // Paramètre de relaxation (peut être ajusté)
	int p = 3;  // Exposant polynomial (cubique par défaut)
	double conv_epsilon = 1e-7;  // Seuil de convergence
	double conv_epsilon_sq = conv_epsilon * conv_epsilon;
	
	// Pour Nova en mode Mandelbrot, utiliser z0 = 1 (racine de z^3 - 1 = 0)
	// Cela évite la singularité à z=0 où p'(0) = 0 et correspond aux références standard
	z = MakeComplex(1.0, 0.0);
	z_prev = z;
	
	do {
		complex z_pow, z_pow_deriv, numerator, denominator, newton_step;
		i++;
		
		// Calcul de z^p
		z_pow = Powz(z, MakeComplex((double)p, 0.0));
		// Calcul de z^(p-1)
		z_pow_deriv = Powz(z, MakeComplex((double)(p-1), 0.0));
		
		// p(z) = z^p - 1
		numerator = Subz(z_pow, ESetofComplex());
		// p'(z) = p·z^(p-1)
		denominator = ScalarTimesofComplex((double)p, z_pow_deriv);
		
		// Éviter division par zéro
		if (Magz(denominator) < 1e-10) {
			break;
		}
		
		// Newton step: p(z)/p'(z)
		newton_step = Divz(numerator, denominator);
		// Multiplier par a (relaxation)
		newton_step = Mulz(a_relax, newton_step);
		
		// z(n+1) = z - a·(p(z)/p'(z)) + c
		z_prev = z;
		z = Addz(Subz(z, newton_step), zPixel);
		
		// Détection de convergence : |z_{n+1} - z_n| devient très petit
		double diff_sq = Magz2(Subz(z, z_prev));
		double z_sq = Magz2(z);
		double denom = (z_sq < 1.0) ? 1.0 : z_sq;
		
		if (diff_sq / denom < conv_epsilon_sq) {
			// Convergence détectée vers une racine
			break;
		}
		
		// Échappement si |z| devient trop grand
		if (Magz(z) > f.bailout) {
			break;
		}
	} while (i < f.iterationMax);
	
	result.iteration = i;
	result.z = z;
	return result;
}
void Nova_def (fractal* f) {
	/* Valeurs de base pour Nova Fractal */
	f->xmin = -3.0;
	f->xmax = 3.0;
	f->ymin = -2.0;
	f->ymax = 2.0;
	f->seed = ZeroSetofComplex();
	f->bailout= 20.0;  // Bailout réduit pour meilleure détection (échappement si |z| > 20)
	f->iterationMax= 500;  // Suffisant pour la convergence
	f->zoomfactor = 4;
}

/* Calcul de Multibrot (Puissances Non-entières)
 * Formule: z(n+1) = z^d + c
 * où d est un nombre réel (non-entier)
 * Utilise la branche principale: z^d = exp(d·log(z))
 * avec log(z) = ln(|z|) + i·arg(z), arg(z) ∈ (-π, π]
 */
fractalresult Multibrot_Iteration (fractal f, complex zPixel) {
	int i=0;
	complex z;
	fractalresult result;
	double d = 2.5;  // Exposant non-entier (exemple: morphing entre 2 et 3)
	
	z = f.seed;
	do {
		complex z_pow;
		i++;
		
		// Calcul de z^d via exponentiation complexe
		// z^d = exp(d·log(z))
		// Utilise la branche principale automatiquement via Logz
		// Powz(0, d) avec d > 0 retourne 0 correctement, donc pas besoin de vérification préalable
		z_pow = Powz(z, MakeComplex(d, 0.0));
		
		// Vérifier si le résultat est valide (NaN ou Inf)
		if (isnan(Rez(z_pow)) || isnan(Imz(z_pow)) || 
		    isinf(Rez(z_pow)) || isinf(Imz(z_pow))) {
			break;
		}
		
		// z(n+1) = z^d + c
		z = Addz(z_pow, zPixel);
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	
	result.iteration = i;
	result.z = z;
	return result;
}
void Multibrot_def (fractal* f) {
	/* Valeurs de base pour Multibrot (puissances non-entières) */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout= 4;
	f->iterationMax= 5000;
	f->zoomfactor = 8;
}

/* Buddhabrot - algorithme de densité par trajectoires */
void Buddhabrot_def (fractal* f) {
	/* Valeurs de base pour Buddhabrot */
	f->xmin = -2.5;
	f->xmax = 1.5;
	f->ymin = -1.5;
	f->ymax = 1.5;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 220;  // Plus d'itérations pour de meilleurs détails
	f->zoomfactor = 4;
}

/* Fonction de rendu spéciale pour Buddhabrot */
Uint32 Buddhabrot_Draw (SDL_Surface *canvas, fractal* f, int decalageX, int decalageY, void* guiPtr) {
	int i, j;
	Uint32 time;
	int numSamples;
	int maxDensity;
	int xpixel = f->xpixel;
	int ypixel = f->ypixel;
	int iterMax = f->iterationMax;
	int bailout = f->bailout;
	double xmin = f->xmin;
	double ymin = f->ymin;
	double xrange = f->xmax - f->xmin;
	double yrange = f->ymax - f->ymin;

	time = SDL_GetTicks();

	printf("Calculating Buddhabrot (density algorithm)...\n");
#ifdef HAVE_OPENMP
	printf("Using OpenMP with %d threads\n", omp_get_max_threads());
#endif

	// Réinitialiser la matrice de densité (on utilise fmatrix)
	for (i = 0; i < xpixel * ypixel; i++) {
		f->fmatrix[i] = 0;
	}

	// Nombre d'échantillons proportionnel à la surface
	numSamples = xpixel * ypixel * 50;
#ifdef HAVE_GMP
	// Réduire les échantillons avec GMP car le calcul est plus lent
	if (f->use_gmp) {
		numSamples = xpixel * ypixel * 5;
	}
#endif

	// Taille des chunks pour rendu incrémental (affichage progressif)
	int chunk_size = 1000;
	if (chunk_size > numSamples / 10) {
		chunk_size = numSamples / 10;
	}
	if (chunk_size < 100) chunk_size = 100;
	int num_chunks = (numSamples + chunk_size - 1) / chunk_size;

	// Mise à jour initiale de la progression
	if (guiPtr != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, 0, "Buddhabrot");
	}

	// Phase 1: Échantillonnage parallèle avec OpenMP et rendu incrémental
#ifdef HAVE_OPENMP
	// Rendu incrémental par chunks
	for (int chunk = 0; chunk < num_chunks; chunk++) {
		int chunk_start = chunk * chunk_size;
		int chunk_end = (chunk_start + chunk_size < numSamples) ? chunk_start + chunk_size : numSamples;
		
		#pragma omp parallel
		{
			// Variables thread-local
			unsigned int seed = 42 + omp_get_thread_num() * 1000 + chunk * 10000;
			double *trajX = (double*) malloc(iterMax * sizeof(double));
			double *trajY = (double*) malloc(iterMax * sizeof(double));

			if (trajX != NULL && trajY != NULL) {
				#pragma omp for schedule(dynamic, 100) nowait
				for (int sample = chunk_start; sample < chunk_end; sample++) {
					// Générer un point c aléatoire (thread-safe)
					double xg = ((double)rand_r(&seed) / RAND_MAX) * xrange + xmin;
					double yg = ((double)rand_r(&seed) / RAND_MAX) * yrange + ymin;
					complex c = MakeComplex(xg, yg);
					complex z = ZeroSetofComplex();

					// Itérer et stocker la trajectoire
					int escaped = 0;
					int iter;
					for (iter = 0; iter < iterMax; iter++) {
						complex zTemp = Mulz(z, z);
						z = Addz(zTemp, c);

						trajX[iter] = Rez(z);
						trajY[iter] = Imz(z);

						if (Magz(z) > bailout) {
							escaped = 1;
							break;
						}
					}

					// Si le point s'échappe, tracer sa trajectoire
					if (escaped) {
						for (int i = 0; i < iter; i++) {
							// Convertir coordonnées complexes en pixels
							int px = (int)((trajX[i] - xmin) / xrange * xpixel);
							int py = (int)((trajY[i] - ymin) / yrange * ypixel);

							// Vérifier les limites et incrémenter la densité (atomic pour éviter race condition)
							if (px >= 0 && px < xpixel && py >= 0 && py < ypixel) {
								#pragma omp atomic
								f->fmatrix[py * xpixel + px]++;
							}
						}
					}
				}
			}

			free(trajX);
			free(trajY);
		} // Fin région parallèle pour ce chunk
		
		// Après chaque chunk : normaliser, coloriser et afficher progressivement
		// Trouver la densité maximale actuelle
		maxDensity = 1;
		for (i = 0; i < xpixel * ypixel; i++) {
			if (f->fmatrix[i] > maxDensity) {
				maxDensity = f->fmatrix[i];
			}
		}
		
		// Convertir densité en couleurs et afficher progressivement
		if (maxDensity > 0) {
			double logMaxDensity = log(1 + maxDensity);
			const gradient_table* gradient = Colorization_GetPalette(f->colorMode);
			// Afficher toutes les 10 lignes pour performance
			for (j = 0; j < ypixel; j++) {
				for (i = 0; i < xpixel; i++) {
					double density = (double)f->fmatrix[j * xpixel + i];
					double normalized = log(1 + density) / logMaxDensity;

					colorization_color cc = Gradient_Interpolate(gradient, normalized);
					/* Update cmatrix for consistency */
					f->cmatrix[j * xpixel + i].r = cc.r;
					f->cmatrix[j * xpixel + i].g = cc.g;
					f->cmatrix[j * xpixel + i].b = cc.b;
					f->cmatrix[j * xpixel + i].a = 255;
					pixelRGBA(canvas, (Sint16)(i + decalageX), (Sint16)(j + decalageY),
					          cc.r, cc.g, cc.b, 255);
				}
				// Mise à jour partielle tous les 10 lignes
				if (j % 10 == 0 || j == ypixel - 1) {
					SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);
				}
			}
		}

		// Mise à jour de la progression
		if (guiPtr != NULL) {
			int percent = ((chunk + 1) * 90) / num_chunks; // 90% pour l'échantillonnage
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, "Buddhabrot");
		}
		
		// Vérifier annulation utilisateur
		if (check_events_and_cancel()) {
			printf("Calcul Buddhabrot annulé par l'utilisateur\n");
			break;
		}
	}
#else
	// Version séquentielle
	double *trajX = (double*) malloc(iterMax * sizeof(double));
	double *trajY = (double*) malloc(iterMax * sizeof(double));
	if (trajX == NULL || trajY == NULL) {
		fprintf(stderr, "Erreur allocation mémoire trajectoire\n");
		return 0;
	}

	srand(42);  // Seed fixe pour reproductibilité

	// Rendu incrémental par chunks (version séquentielle)
	for (int chunk = 0; chunk < num_chunks; chunk++) {
		int chunk_start = chunk * chunk_size;
		int chunk_end = (chunk_start + chunk_size < numSamples) ? chunk_start + chunk_size : numSamples;
		
		for (int sample = chunk_start; sample < chunk_end; sample++) {
			double xg = ((double)rand() / RAND_MAX) * xrange + xmin;
			double yg = ((double)rand() / RAND_MAX) * yrange + ymin;
			complex c = MakeComplex(xg, yg);
			complex z = ZeroSetofComplex();

			int escaped = 0;
			int iter;
			for (iter = 0; iter < iterMax; iter++) {
				complex zTemp = Mulz(z, z);
				z = Addz(zTemp, c);

				trajX[iter] = Rez(z);
				trajY[iter] = Imz(z);

				if (Magz(z) > bailout) {
					escaped = 1;
					break;
				}
			}

			if (escaped) {
				for (int i = 0; i < iter; i++) {
					int px = (int)((trajX[i] - xmin) / xrange * xpixel);
					int py = (int)((trajY[i] - ymin) / yrange * ypixel);

					if (px >= 0 && px < xpixel && py >= 0 && py < ypixel) {
						f->fmatrix[py * xpixel + px]++;
					}
				}
			}
		}
		
		// Après chaque chunk : normaliser, coloriser et afficher progressivement
		maxDensity = 1;
		for (i = 0; i < xpixel * ypixel; i++) {
			if (f->fmatrix[i] > maxDensity) {
				maxDensity = f->fmatrix[i];
			}
		}
		
		if (maxDensity > 0) {
			double logMaxDensity = log(1 + maxDensity);
			const gradient_table* gradient = Colorization_GetPalette(f->colorMode);
			for (j = 0; j < ypixel; j++) {
				for (i = 0; i < xpixel; i++) {
					double density = (double)f->fmatrix[j * xpixel + i];
					double normalized = log(1 + density) / logMaxDensity;

					colorization_color cc = Gradient_Interpolate(gradient, normalized);
					/* Update cmatrix for consistency */
					f->cmatrix[j * xpixel + i].r = cc.r;
					f->cmatrix[j * xpixel + i].g = cc.g;
					f->cmatrix[j * xpixel + i].b = cc.b;
					f->cmatrix[j * xpixel + i].a = 255;
					pixelRGBA(canvas, (Sint16)(i + decalageX), (Sint16)(j + decalageY),
					          cc.r, cc.g, cc.b, 255);
				}
				if (j % 10 == 0 || j == ypixel - 1) {
					SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);
				}
			}
		}

		// Mise à jour de la progression
		if (guiPtr != NULL) {
			int percent = ((chunk + 1) * 90) / num_chunks;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, "Buddhabrot");
		}
		
		// Vérifier annulation utilisateur
		if (check_events_and_cancel()) {
			printf("Calcul Buddhabrot annulé par l'utilisateur\n");
			free(trajX);
			free(trajY);
			return SDL_GetTicks() - time;
		}
	}

	free(trajX);
	free(trajY);
#endif

	// Phase finale : Affichage final complet
	if (guiPtr != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, 100, "Buddhabrot");
	}
	
	SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);

	time = SDL_GetTicks() - time;
	printf("Buddhabrot rendered in %d ms (%d samples)\n", time, numSamples);
	return time;
}

// *************************************************************************
// *	Lyapunov Fractal
// *************************************************************************

void Lyapunov_def (fractal* f) {
	/* Zircon City : Lyapunov fractal pour la séquence bbbbbbaaaaaa
	 * Domaine : (a, b) ∈ [2.5, 3.4] × [3.4, 4.0] */
	f->xmin = 2.5;
	f->xmax = 3.4;
	f->ymin = 3.4;
	f->ymax = 4.0;
	f->seed = ZeroSetofComplex();
	f->bailout = 4;
	f->iterationMax = 2000;  // Itérations pour calculer l'exposant
	f->zoomfactor = 2;
}

/* Fonction de normalisation de l'exposant de Lyapunov
 * Convertit l'exposant en une valeur normalisée (0.0-1.0) pour les palettes
 * - Zones stables (lyap < 0) : mapper vers 0.0-0.85
 * - Zones chaotiques (lyap >= 0) : mapper vers 0.85-1.0
 */
static double Lyapunov_GetNormalizedValue(double lyap) {
	if (lyap < 0) {
		// Exposant négatif = stable
		// Normaliser entre 0.0 et 0.85 (éviter 1.0 pour distinguer du chaos)
		double t = -lyap;
		if (t > 2.0) t = 2.0;  // Limiter à 2.0
		return (t / 2.0) * 0.85;  // 0.0 à 0.85
	} else {
		// Exposant positif = chaotique
		// Normaliser entre 0.85 et 1.0
		double t = lyap;
		if (t > 1.0) t = 1.0;  // Limiter à 1.0
		return 0.85 + (t * 0.15);  // 0.85 à 1.0
	}
}

/* Fonction générique de rendu pour Lyapunov avec séquence paramétrable
 * Calcule l'exposant de Lyapunov de la suite logistique x_{n+1} = r_n * x_n * (1 - x_n)
 * avec r alternant entre a (x) et b (y) selon la séquence fournie
 */
static Uint32 Lyapunov_Draw_Sequence (SDL_Surface *canvas, fractal* f, int decalageX, int decalageY, void* guiPtr, const char* sequence, const char* fractalName) {
	int i, j;
	double xStep, yStep;
	Uint32 time;
	int seqLen = strlen(sequence);
	int warmup = 50;  // Itérations de stabilisation
	int iterMax = f->iterationMax;
	int xpixel = f->xpixel;
	int ypixel = f->ypixel;

	time = SDL_GetTicks();
	printf("Calculating %s fractal (sequence: %s)...\n", fractalName, sequence);
	printf("Domain: x=[%.6f, %.6f] y=[%.6f, %.6f]\n", f->xmin, f->xmax, f->ymin, f->ymax);
#ifdef HAVE_OPENMP
	printf("Using OpenMP with %d threads\n", omp_get_max_threads());
#endif

	xStep = (f->xmax - f->xmin) / xpixel;
	yStep = (f->ymax - f->ymin) / ypixel;

	// Phase 1: Calcul parallèle des exposants et couleurs
#ifdef HAVE_OPENMP
	#pragma omp parallel for private(i) schedule(dynamic, 16)
#endif
	for (j = 0; j < ypixel; j++) {
		double b = f->ymin + j * yStep;  // Paramètre b sur l'axe Y

		for (i = 0; i < xpixel; i++) {
			double a = f->xmin + i * xStep;  // Paramètre a sur l'axe X
			double x = 0.5;  // Valeur initiale classique
			double lyap = 0.0;
			double r;
			int n;

			// Phase de stabilisation (warmup)
			for (n = 0; n < warmup; n++) {
				r = (sequence[n % seqLen] == 'A') ? a : b;
				x = r * x * (1.0 - x);
				if (x < 0.0001 || x > 0.9999) x = 0.5;  // Éviter les divergences
			}

			// Calcul de l'exposant de Lyapunov
			for (n = 0; n < iterMax; n++) {
				r = (sequence[n % seqLen] == 'A') ? a : b;
				x = r * x * (1.0 - x);

				// Éviter log(0) et les valeurs invalides
				double deriv = fabs(r * (1.0 - 2.0 * x));
				if (deriv > 0.0001) {
					lyap += log(deriv);
				}

				// Réinitialiser si x diverge
				if (x < 0.0001 || x > 0.9999) {
					x = 0.5;
				}
			}

			lyap /= iterMax;

			// Convertir l'exposant en valeur normalisée pour les palettes
			double normalizedValue = Lyapunov_GetNormalizedValue(lyap);
			
			// Stocker la valeur normalisée dans fmatrix comme "nombre d'itérations"
			// Multiplier par iterationMax pour que Fractal_SmoothIteration fonctionne correctement
			f->fmatrix[j * xpixel + i] = (int)(normalizedValue * iterMax);
			
			// Initialiser zmatrix avec une valeur fictive pour les types Lyapunov
			// On utilise une magnitude basée sur la valeur normalisée pour que Fractal_SmoothIteration
			// puisse fonctionner correctement
			complex z_fake;
			z_fake.x = normalizedValue * 2.0;  // Valeur entre 0 et 2.0
			z_fake.y = 0.0;
			f->zmatrix[j * xpixel + i] = z_fake;
		}
	}

	// Phase 2: Calcul de la colorisation avec les palettes
	// Réinitialiser cmatrix_valid pour forcer le recalcul
	f->cmatrix_valid = 0;
	int progress = 0;
	Fractal_CalculateColorMatrix(f, canvas, guiPtr, &progress, 0, 100);

	// Phase 3: Rendu SDL séquentiel (SDL n'est pas thread-safe)
	for (j = 0; j < ypixel; j++) {
		for (i = 0; i < xpixel; i++) {
			color col = f->cmatrix[j * xpixel + i];
			pixelRGBA(canvas, (Sint16)(i + decalageX), (Sint16)(j + decalageY),
			          col.r, col.g, col.b, 255);
		}

		// Mise à jour de l'affichage tous les 50 lignes
		if (j % 50 == 0) {
			SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);
			if (guiPtr != NULL) {
				int percent = (j * 100) / ypixel;
				SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
			}
		}
	}

	SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);

	time = SDL_GetTicks() - time;
	printf("%s fractal rendered in %d ms\n", fractalName, time);
	return time;
}

/* Fonction de rendu pour Lyapunov Zircon City (séquence "BBBBBBAAAAAA" - 6B puis 6A) */
Uint32 Lyapunov_Draw (SDL_Surface *canvas, fractal* f, int decalageX, int decalageY, void* guiPtr) {
	return Lyapunov_Draw_Sequence(canvas, f, decalageX, decalageY, guiPtr, "BBBBBBAAAAAA", "Lyapunov Zircon City");
}

/* Variantes Lyapunov 18-27 supprimées - seule Zircon City (type 17) est conservée */

#ifdef HAVE_GMP
// ******************************************************
// Versions GMP des fonctions d'itération
// ******************************************************

fractalresult Mendelbrot_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp;
	mpf_t mag2, bailout_mpf, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);

	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zTemp, prec);  // Pré-allouer zTemp
	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);  // Temporaires pour mag2
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);  // bailout² pour comparaison avec |z|²

	do {
		i++;
		// Utiliser fonctions in-place (pas d'allocation dans la boucle)
		complex_gmp_sq_to(&zTemp, z, &f.mul_temps);  // zTemp = z²
		complex_gmp_add_to(&z, zTemp, zPixel);       // z = zTemp + c
		complex_gmp_mag2_to(mag2, z, temp1, temp2);  // |z|² (sans sqrt)
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult Julia_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, seed_gmp;
	mpf_t mag2, bailout_mpf, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(zPixel, prec);
	complex_gmp_init(&zTemp, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);  // bailout²

	do {
		i++;
		complex_gmp_sq_to(&zTemp, z, &f.mul_temps);
		complex_gmp_add_to(&z, zTemp, seed_gmp);
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult JuliaSin_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, seed_gmp, sin_z;
	mpf_t mag, bailout_mpf;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(zPixel, prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	
	do {
		i++;
		sin_z = complex_gmp_sin(z, prec);
		complex_gmp old_z = z;
		z = complex_gmp_mul(seed_gmp, sin_z, prec);
		complex_gmp_clear(&old_z);
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&sin_z);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	
	return result;
}

fractalresult Newton_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zNum, zQuot, p_minus_one, p_val, one, p_minus_one_z;
	mpf_t mag, bailout_mpf, p_mpf, p_minus_one_mpf, zero_mpf;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	z = complex_gmp_copy(zPixel, prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(p_mpf, prec);
	mpf_init2(p_minus_one_mpf, prec);
	mpf_init2(zero_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_ui(zero_mpf, 0);
	
	double p_d = Rez(f.seed);
	mpf_set_d(p_mpf, p_d);
	mpf_set_d(p_minus_one_mpf, p_d - 1.0);
	
	one = complex_gmp_one(prec);
	p_val = complex_gmp_make(p_mpf, zero_mpf, prec);
	p_minus_one = complex_gmp_make(p_minus_one_mpf, zero_mpf, prec);
	
	do {
		complex_gmp z_pow_p, z_pow_p_minus_one;
		i++;
		
		// z^p
		z_pow_p = complex_gmp_pow_n(z, p_val, 0, prec);
		// (p-1)*z^p + 1
		p_minus_one_z = complex_gmp_mul(p_minus_one, z_pow_p, prec);
		zNum = complex_gmp_add(p_minus_one_z, one, prec);
		
		// z^(p-1)
		z_pow_p_minus_one = complex_gmp_pow_n(z, p_minus_one, 0, prec);
		// p*z^(p-1)
		zQuot = complex_gmp_mul(p_val, z_pow_p_minus_one, prec);
		
		// z = zNum / zQuot
		z = complex_gmp_div(zNum, zQuot, prec);
		
		complex_gmp_mag(mag, z);
		
		complex_gmp_clear(&z_pow_p);
		complex_gmp_clear(&z_pow_p_minus_one);
		complex_gmp_clear(&p_minus_one_z);
		complex_gmp_clear(&zNum);
		complex_gmp_clear(&zQuot);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&one);
	complex_gmp_clear(&p_val);
	complex_gmp_clear(&p_minus_one);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(p_mpf);
	mpf_clear(p_minus_one_mpf);
	mpf_clear(zero_mpf);
	
	return result;
}

fractalresult Phoenix_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, y, zTemp, zpTemp;
	mpf_t mag, bailout_mpf, p1, p2, z_real;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	z = complex_gmp_copy(zPixel, prec);
	y = complex_gmp_zero(prec);
	zTemp = complex_gmp_zero(prec);  // Initialiser zTemp
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(p1, prec);
	mpf_init2(p2, prec);
	mpf_init2(z_real, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_d(p1, 0.56667);
	mpf_set_d(p2, -0.5);
	
	do {
		i++;
		complex_gmp old_zTemp = zTemp;
		zTemp = complex_gmp_mul(z, z, prec);
		complex_gmp_clear(&old_zTemp);
		complex_gmp_get_real(z_real, zTemp);
		mpf_add(z_real, z_real, p1);
		complex_gmp_set_real(&zTemp, z_real);
		
		zpTemp = complex_gmp_scalar_mul(p2, y, prec);
		complex_gmp old_zTemp2 = zTemp;
		zTemp = complex_gmp_add(zTemp, zpTemp, prec);
		complex_gmp_clear(&old_zTemp2);
		
		// Sauvegarder les anciennes valeurs avant réassignation
		complex_gmp old_y = y;
		complex_gmp old_z = z;
		// Copier les valeurs (les structures mpf_t sont copiées, pas les pointeurs)
		y = z;  // y pointe maintenant vers la même mémoire que z
		z = zTemp;  // z pointe maintenant vers la mémoire de zTemp
		// Libérer les anciennes valeurs (y et z pointent maintenant vers d'autres mémoires)
		complex_gmp_clear(&old_y);
		complex_gmp_clear(&old_z);
		
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zpTemp);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&y);
	complex_gmp_clear(&zTemp);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(p1);
	mpf_clear(p2);
	mpf_clear(z_real);
	
	return result;
}

/* Sierpinski_Iteration_GMP supprimé - type 8 */

fractalresult Barnsley1j_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, seed_gmp, one;
	mpf_t mag, bailout_mpf, z_real, one_mpf;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(zPixel, prec);
	one = complex_gmp_one(prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(z_real, prec);
	mpf_init2(one_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_ui(one_mpf, 1);
	
	do {
		i++;
		complex_gmp_get_real(z_real, z);
		
		if (mpf_cmp(z_real, one_mpf) >= 0) {
			zTemp = complex_gmp_sub(z, one, prec);
		} else {
			zTemp = complex_gmp_add(z, one, prec);
		}
		complex_gmp old_z = z;
		z = complex_gmp_mul(zTemp, seed_gmp, prec);
		complex_gmp_clear(&old_z);
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zTemp);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	complex_gmp_clear(&one);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(z_real);
	mpf_clear(one_mpf);
	
	return result;
}

fractalresult Barnsley1m_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, c, zTemp, one;
	mpf_t mag, bailout_mpf, z_real, one_mpf;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	c = complex_gmp_copy(zPixel, prec);
	z = complex_gmp_zero(prec);
	one = complex_gmp_one(prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(z_real, prec);
	mpf_init2(one_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_ui(one_mpf, 1);
	
	do {
		i++;
		complex_gmp_get_real(z_real, z);
		
		if (mpf_cmp(z_real, one_mpf) >= 0) {
			zTemp = complex_gmp_sub(z, one, prec);
		} else {
			zTemp = complex_gmp_add(z, one, prec);
		}
		complex_gmp old_z = z;
		z = complex_gmp_mul(zTemp, c, prec);
		complex_gmp_clear(&old_z);
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zTemp);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&c);
	complex_gmp_clear(&one);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(z_real);
	mpf_clear(one_mpf);
	
	return result;
}

fractalresult BurningShip_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, zAbs;
	mpf_t mag2, bailout_mpf, abs_real, abs_imag, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zTemp, prec);
	complex_gmp_init(&zAbs, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(abs_real, prec);
	mpf_init2(abs_imag, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);  // bailout²

	do {
		i++;
		// |Re(z)| et |Im(z)| - accès direct sans copie
		mpf_abs(abs_real, z.x);
		mpf_abs(abs_imag, z.y);

		// zAbs = (|Re(z)|, |Im(z)|) - in-place
		mpf_set(zAbs.x, abs_real);
		mpf_set(zAbs.y, abs_imag);

		complex_gmp_sq_to(&zTemp, zAbs, &f.mul_temps);  // zTemp = zAbs²
		complex_gmp_add_to(&z, zTemp, zPixel);          // z = zTemp + c
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&zAbs);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(abs_real);
	mpf_clear(abs_imag);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

/* Buffalo GMP: z(n+1) = abs(Re(z²)) + i*abs(Im(z²)) + c */
fractalresult Buffalo_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zSq;
	mpf_t mag2, bailout_mpf, abs_real, abs_imag, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zSq, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(abs_real, prec);
	mpf_init2(abs_imag, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);  // bailout²

	do {
		i++;
		// zSq = z²
		complex_gmp_sq_to(&zSq, z, &f.mul_temps);
		// Appliquer abs() aux deux parties de z²
		mpf_abs(abs_real, zSq.x);
		mpf_abs(abs_imag, zSq.y);
		// z = (|Re(z²)|, |Im(z²)|) + c
		mpf_add(z.x, abs_real, zPixel.x);
		mpf_add(z.y, abs_imag, zPixel.y);
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zSq);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(abs_real);
	mpf_clear(abs_imag);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult Tricorn_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, zConj;
	mpf_t mag2, bailout_mpf, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&zTemp, prec);
	complex_gmp_init(&zConj, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);  // bailout²

	do {
		i++;
		// Conjugué: (x, -y) - in-place
		mpf_set(zConj.x, z.x);
		mpf_neg(zConj.y, z.y);

		complex_gmp_sq_to(&zTemp, zConj, &f.mul_temps);  // zTemp = conj(z)²
		complex_gmp_add_to(&z, zTemp, zPixel);           // z = zTemp + c
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&zTemp);
	complex_gmp_clear(&zConj);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult Mandelbulb_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, z2, z4, z8;
	mpf_t mag2, bailout_mpf, temp1, temp2;
	fractalresult result;

	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	complex_gmp_init(&z2, prec);
	complex_gmp_init(&z4, prec);
	complex_gmp_init(&z8, prec);

	mpf_init2(mag2, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(temp1, prec);
	mpf_init2(temp2, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_mul(bailout_mpf, bailout_mpf, bailout_mpf);  // bailout²

	do {
		i++;
		// z^8 par carrés successifs (in-place, pas d'allocation)
		complex_gmp_sq_to(&z2, z, &f.mul_temps);   // z²
		complex_gmp_sq_to(&z4, z2, &f.mul_temps);  // z⁴
		complex_gmp_sq_to(&z8, z4, &f.mul_temps);  // z⁸
		complex_gmp_add_to(&z, z8, zPixel);        // z = z⁸ + c
		complex_gmp_mag2_to(mag2, z, temp1, temp2);
	} while ((i < f.iterationMax) && (mpf_cmp(mag2, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&z2);
	complex_gmp_clear(&z4);
	complex_gmp_clear(&z8);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag2);
	mpf_clear(bailout_mpf);
	mpf_clear(temp1);
	mpf_clear(temp2);

	return result;
}

fractalresult Magnet1j_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, N, Q, seed_gmp, one, two, seed_minus_one, seed_minus_two;
	mpf_t mag, bailout_mpf;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(zPixel, prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	
	// Constantes
	mpf_t two_mpf;
	mpf_init2(two_mpf, prec);
	mpf_set_ui(two_mpf, 2);
	one = complex_gmp_one(prec);
	two = complex_gmp_scalar_mul(two_mpf, one, prec);
	mpf_clear(two_mpf);
	
	// seed - 1 et seed - 2
	seed_minus_one = complex_gmp_sub(seed_gmp, one, prec);
	seed_minus_two = complex_gmp_sub(seed_gmp, two, prec);
	
	do {
		i++;
		// N = (z^2 + seed - 1)^2
		complex_gmp z_sq = complex_gmp_mul(z, z, prec);
		complex_gmp N_temp = complex_gmp_add(z_sq, seed_minus_one, prec);
		N = complex_gmp_mul(N_temp, N_temp, prec);
		complex_gmp_clear(&N_temp);

		// Q = 2*z + seed - 2
		complex_gmp two_z = complex_gmp_mul(two, z, prec);
		Q = complex_gmp_add(two_z, seed_minus_two, prec);

		// z = N / Q (libérer l'ancien z d'abord)
		complex_gmp old_z = z;
		z = complex_gmp_div(N, Q, prec);
		complex_gmp_clear(&old_z);

		complex_gmp_mag(mag, z);

		complex_gmp_clear(&z_sq);
		complex_gmp_clear(&two_z);
		complex_gmp_clear(&N);
		complex_gmp_clear(&Q);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	complex_gmp_clear(&one);
	complex_gmp_clear(&two);
	complex_gmp_clear(&seed_minus_one);
	complex_gmp_clear(&seed_minus_two);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);

	return result;
}

fractalresult Magnet1m_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, c, N, Q, one, two, c_minus_one, c_minus_two;
	mpf_t mag, bailout_mpf;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	c = complex_gmp_copy(zPixel, prec);
	z = complex_gmp_zero(prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	
	// Constantes
	mpf_t two_mpf;
	mpf_init2(two_mpf, prec);
	mpf_set_ui(two_mpf, 2);
	one = complex_gmp_one(prec);
	two = complex_gmp_scalar_mul(two_mpf, one, prec);
	mpf_clear(two_mpf);
	
	// c - 1 et c - 2
	c_minus_one = complex_gmp_sub(c, one, prec);
	c_minus_two = complex_gmp_sub(c, two, prec);
	
	do {
		i++;
		// N = (z^2 + c - 1)^2
		complex_gmp z_sq = complex_gmp_mul(z, z, prec);
		complex_gmp N_temp = complex_gmp_add(z_sq, c_minus_one, prec);
		N = complex_gmp_mul(N_temp, N_temp, prec);
		complex_gmp_clear(&N_temp);

		// Q = 2*z + c - 2
		complex_gmp two_z = complex_gmp_mul(two, z, prec);
		Q = complex_gmp_add(two_z, c_minus_two, prec);

		// z = N / Q (libérer l'ancien z d'abord)
		complex_gmp old_z = z;
		z = complex_gmp_div(N, Q, prec);
		complex_gmp_clear(&old_z);

		complex_gmp_mag(mag, z);

		complex_gmp_clear(&z_sq);
		complex_gmp_clear(&two_z);
		complex_gmp_clear(&N);
		complex_gmp_clear(&Q);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));

	result.iteration = i;
	result.z = gmp_to_complex(z);

	complex_gmp_clear(&z);
	complex_gmp_clear(&c);
	complex_gmp_clear(&one);
	complex_gmp_clear(&two);
	complex_gmp_clear(&c_minus_one);
	complex_gmp_clear(&c_minus_two);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);

	return result;
}

fractalresult FormulaSelector_GMP (fractal f, complex_gmp zPixel) {
	fractalresult r;
	
	switch (f.type) {
	case 3:
		r = Mendelbrot_Iteration_GMP(f, zPixel);
		break;
	case 4:
		r = Julia_Iteration_GMP(f, zPixel);
		break;
	case 5:
		r = JuliaSin_Iteration_GMP(f, zPixel);
		break;
	case 6:
		r = Newton_Iteration_GMP(f, zPixel);
		break;
	case 7:
		r = Phoenix_Iteration_GMP(f, zPixel);
		break;
	case 8:
		r = Buffalo_Iteration_GMP(f, zPixel);
		break;
	case 9:
		r = Barnsley1j_Iteration_GMP(f, zPixel);
		break;
	case 10:
		r = Barnsley1m_Iteration_GMP(f, zPixel);
		break;
	case 11:
		r = Magnet1j_Iteration_GMP(f, zPixel);
		break;
	case 12:
		r = Magnet1m_Iteration_GMP(f, zPixel);
		break;
	case 13:
		r = BurningShip_Iteration_GMP(f, zPixel);
		break;
	case 14:
		r = Tricorn_Iteration_GMP(f, zPixel);
		break;
	case 15:
		r = Mandelbulb_Iteration_GMP(f, zPixel);
		break;
	default:
		r = Mendelbrot_Iteration_GMP(f, zPixel);
		break;
	}
	
	return r;
}

#endif /* HAVE_GMP */

// Selectionne la formule
fractalresult FormulaSelector (fractal f, complex zPixel) {
	fractalresult r;

#ifdef HAVE_GMP
	if (f.use_gmp) {
		complex_gmp zPixel_gmp = complex_to_gmp(zPixel, f.gmp_precision);
		r = FormulaSelector_GMP(f, zPixel_gmp);
		complex_gmp_clear(&zPixel_gmp);
		return r;
	}
#endif

	switch (f.type) {
	case 3:
	r = Mendelbrot_Iteration (f,zPixel) ;
	break;
	case 4:
        r= Julia_Iteration (f,zPixel);
	break;
	case 5:
	r=JuliaSin_Iteration (f,zPixel);
	break;
	case 6:
	r=Newton_Iteration (f,zPixel);
	break;
	case 7:
	r=Phoenix_Iteration (f,zPixel);
	break;
	case 8:
	r = Buffalo_Iteration(f, zPixel);
	break;
	case 9:
	r=Barnsley1j_Iteration (f,zPixel);
	break;
	case 10:
	r=Barnsley1m_Iteration (f,zPixel);
	break;
	case 11:
	r=Magnet1j_Iteration (f,zPixel);
	break;
	case 12:
	r=Magnet1m_Iteration (f,zPixel);
	break;
	case 13:
	r=BurningShip_Iteration (f,zPixel);
	break;
	case 14:
	r=Tricorn_Iteration (f,zPixel);
	break;
	case 15:
	r=Mandelbulb_Iteration (f,zPixel);
	break;
	case 18:
	r=PerpendicularBurningShip_Iteration (f,zPixel);
	break;
	case 19:
	r=Celtic_Iteration (f,zPixel);
	break;
	case 20:
	r=AlphaMandelbrot_Iteration (f,zPixel);
	break;
	case 21:
	r=PickoverStalks_Iteration (f,zPixel);
	break;
	case 22:
	r=Nova_Iteration (f,zPixel);
	break;
	case 23:
	r=Multibrot_Iteration (f,zPixel);
	break;

	}

		return r;
}
