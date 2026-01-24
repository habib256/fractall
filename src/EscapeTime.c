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
// Ne traite que les événements clavier et quit pour éviter d'interférer avec la souris
static int check_events_and_cancel(void) {
    SDL_Event events[32];
    int num_events, i;
    int should_cancel = 0;

    // Pomper les événements pour s'assurer qu'ils sont disponibles
    SDL_PumpEvents();

    // Récupérer TOUS les événements QUIT et clavier, consommer ceux qu'on veut
    // Utiliser SDL_GETEVENT pour consommer les événements immédiatement
    int num_all = SDL_PeepEvents(events, 32, SDL_GETEVENT, SDL_QUITMASK | SDL_KEYDOWNMASK);
    
    // Parcourir tous les événements
    for (i = 0; i < num_all; i++) {
        if (events[i].type == SDL_QUIT) {
            should_cancel = 1;
            // Ne pas remettre l'événement QUIT - il est consommé
            continue;
        }
        
        if (events[i].type == SDL_KEYDOWN) {
            if (events[i].key.keysym.sym == SDLK_ESCAPE ||
                events[i].key.keysym.sym == SDLK_q) {
                should_cancel = 1;
                // Ne pas remettre ESC/Q - ils sont consommés
                continue;
            }
        }
        
        // Remettre les autres événements pour traitement ultérieur
        SDL_PushEvent(&events[i]);
    }

    return should_cancel;
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
	// Répétition par défaut selon le type : 40 pour escape-time, 2 pour Lyapunov
	if (type == 17) {
		f.colorRepeat = 40;  // 40 répétitions par défaut (comme les autres fractales)
	} else {
		f.colorRepeat = 40;  // 40 répétitions pour escape-time
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
	int chunk_height = num_threads * 2;  // Chunks plus petits pour meilleure réactivité
	if (chunk_height < 16) chunk_height = 16;  // Minimum 16 lignes

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
	// Protection contre division par zéro
	if (xpixel <= 0 || ypixel <= 0) {
		fprintf(stderr, "Erreur: dimensions invalides pour DDp1 (xpixel=%d, ypixel=%d)\n", xpixel, ypixel);
		return;
	}
	double xstep = (f->xmax - f->xmin) / xpixel;
	double ystep = (f->ymax - f->ymin) / ypixel;

	printf ("Calculating EscapeTimeFractal Matrix with DDp1 ...\n");
#ifdef HAVE_OPENMP
	int num_threads = omp_get_max_threads();
	printf("Using OpenMP with %d threads\n", num_threads);
#endif

	// Mise à jour initiale de la progression
	if (guiPtr != NULL && progress != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressStart, fractalName);
		*progress = progressStart;
	}

	// Traitement par chunks pour permettre l'annulation et la mise à jour de la progression
	int chunk_height = 32;  // Traiter 32 lignes à la fois (16 itérations car j += 2)
	int cancelled = 0;

#ifdef HAVE_OPENMP
	for (int chunk_start = 0; chunk_start < ypixel && !cancelled; chunk_start += chunk_height) {
		int chunk_end = (chunk_start + chunk_height < ypixel) ? chunk_start + chunk_height : ypixel;

		#pragma omp parallel for schedule(guided)
		for (j = chunk_start; j < chunk_end; j += 2) {
			complex zPixel_local;
			fractalresult result_local;
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

		// Vérifier annulation et mettre à jour progression après chaque chunk
		if (check_events_and_cancel()) {
			cancelled = 1;
			printf("Calcul DDp1 annulé par l'utilisateur\n");
		}
		if (guiPtr != NULL && progress != NULL && !cancelled) {
			int percent = progressStart + ((chunk_start + chunk_height) * (progressEnd - progressStart)) / ypixel;
			if (percent > progressEnd) percent = progressEnd;
			*progress = percent;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
		}
	}
#else
	for (int chunk_start = 0; chunk_start < ypixel && !cancelled; chunk_start += chunk_height) {
		int chunk_end = (chunk_start + chunk_height < ypixel) ? chunk_start + chunk_height : ypixel;

		for (j = chunk_start; j < chunk_end; j += 2) {
			complex zPixel_local;
			fractalresult result_local;
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

		// Vérifier annulation et mettre à jour progression après chaque chunk
		if (check_events_and_cancel()) {
			cancelled = 1;
			printf("Calcul DDp1 annulé par l'utilisateur\n");
		}
		if (guiPtr != NULL && progress != NULL && !cancelled) {
			int percent = progressStart + ((chunk_start + chunk_height) * (progressEnd - progressStart)) / ypixel;
			if (percent > progressEnd) percent = progressEnd;
			*progress = percent;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
		}
	}
#endif

	// Mise à jour finale de la progression
	if (guiPtr != NULL && progress != NULL && !cancelled) {
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
	int chunk_height_p1 = num_threads * 2;  // Chunks plus petits pour meilleure réactivité
	if (chunk_height_p1 < 16) chunk_height_p1 = 16;  // Minimum 16 lignes
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
	int chunk_height_p2 = num_threads * 2;  // Chunks plus petits pour meilleure réactivité
	if (chunk_height_p2 < 16) chunk_height_p2 = 16;  // Minimum 16 lignes
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
						// Protection contre division par zéro (compteur devrait toujours être > 0 ici)
						if (compteur > 0) {
							rz = (Rez(zup_local)+Rez(zdown_local)+Rez(zleft_local)+Rez(zright_local))/compteur;
							iz = (Imz(zup_local)+Imz(zdown_local)+Imz(zleft_local)+Imz(zright_local))/compteur;
						} else {
							// Cas de sécurité : utiliser la moyenne des valeurs disponibles
							rz = (Rez(zup_local)+Rez(zdown_local)+Rez(zleft_local)+Rez(zright_local))/4.0;
							iz = (Imz(zup_local)+Imz(zdown_local)+Imz(zleft_local)+Imz(zright_local))/4.0;
						}
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
					// Protection contre division par zéro (compteur devrait toujours être > 0 ici)
					if (compteur > 0) {
						rz = (Rez(zup)+Rez(zdown)+Rez(zleft)+Rez(zright))/compteur;
						iz = (Imz(zup)+Imz(zdown)+Imz(zleft)+Imz(zright))/compteur;
					} else {
						// Cas de sécurité : utiliser la moyenne des valeurs disponibles
						rz = (Rez(zup)+Rez(zdown)+Rez(zleft)+Rez(zright))/4.0;
						iz = (Imz(zup)+Imz(zdown)+Imz(zleft)+Imz(zright))/4.0;
					}
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
	// Protection contre division par zéro
	if (xpixel <= 0 || ypixel <= 0) {
		fprintf(stderr, "Erreur: dimensions invalides pour DDp1 (xpixel=%d, ypixel=%d)\n", xpixel, ypixel);
		return;
	}
	double xstep = (f->xmax - f->xmin) / xpixel;
	double ystep = (f->ymax - f->ymin) / ypixel;
	int progressMid = progressStart + (progressEnd - progressStart) / 2;

	printf ("Calculating EscapeTimeFractal Matrix with DDp2 ...\n");
#ifdef HAVE_OPENMP
	int num_threads = omp_get_max_threads();
	printf("Using OpenMP with %d threads\n", num_threads);
#endif

	// Mise à jour initiale de la progression
	if (guiPtr != NULL && progress != NULL) {
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressStart, fractalName);
		*progress = progressStart;
	}

	// Traitement par chunks pour permettre l'annulation et la mise à jour de la progression
	int chunk_height = 32;  // Traiter 32 lignes à la fois
	int cancelled = 0;

	// Pass 1: Calcul initial des pixels impairs (parallélisable)
#ifdef HAVE_OPENMP
	for (int chunk_start = 1; chunk_start < ypixel && !cancelled; chunk_start += chunk_height) {
		int chunk_end = (chunk_start + chunk_height < ypixel) ? chunk_start + chunk_height : ypixel;

		#pragma omp parallel for schedule(guided)
		for (j = chunk_start; j < chunk_end; j += 2) {
			complex zPixel_local;
			fractalresult result_local;
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

		// Vérifier annulation et mettre à jour progression après chaque chunk
		if (check_events_and_cancel()) {
			cancelled = 1;
			printf("Calcul DDp2 Pass 1 annulé par l'utilisateur\n");
		}
		if (guiPtr != NULL && progress != NULL && !cancelled) {
			int percent = progressStart + ((chunk_start) * (progressMid - progressStart)) / ypixel;
			if (percent > progressMid) percent = progressMid;
			*progress = percent;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
		}
	}
#else
	for (int chunk_start = 1; chunk_start < ypixel && !cancelled; chunk_start += chunk_height) {
		int chunk_end = (chunk_start + chunk_height < ypixel) ? chunk_start + chunk_height : ypixel;

		for (j = chunk_start; j < chunk_end; j += 2) {
			complex zPixel_local;
			fractalresult result_local;
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

		// Vérifier annulation et mettre à jour progression après chaque chunk
		if (check_events_and_cancel()) {
			cancelled = 1;
			printf("Calcul DDp2 Pass 1 annulé par l'utilisateur\n");
		}
		if (guiPtr != NULL && progress != NULL && !cancelled) {
			int percent = progressStart + ((chunk_start) * (progressMid - progressStart)) / ypixel;
			if (percent > progressMid) percent = progressMid;
			*progress = percent;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
		}
	}
#endif

	// Mise à jour de la progression après Pass 1
	if (guiPtr != NULL && progress != NULL && !cancelled) {
		*progress = progressMid;
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressMid, fractalName);
	}

	// Pass 2: Divergence Detection (parallélisable - chaque pixel écrit à son propre indice)
#ifdef HAVE_OPENMP
	for (int chunk_start = 0; chunk_start < ypixel && !cancelled; chunk_start += chunk_height) {
		int chunk_end = (chunk_start + chunk_height < ypixel) ? chunk_start + chunk_height : ypixel;

		#pragma omp parallel for schedule(guided)
		for (j = chunk_start; j < chunk_end; j++) {
			complex zPixel_local;
			fractalresult result_local;
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
					double rz, iz;
					if (compteur > 0) {
						rz = (Rez(zup) + Rez(zdown) + Rez(zleft) + Rez(zright)) / compteur;
						iz = (Imz(zup) + Imz(zdown) + Imz(zleft) + Imz(zright)) / compteur;
					} else {
						rz = (Rez(zup) + Rez(zdown) + Rez(zleft) + Rez(zright)) / 4.0;
						iz = (Imz(zup) + Imz(zdown) + Imz(zleft) + Imz(zright)) / 4.0;
					}
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

		// Vérifier annulation et mettre à jour progression après chaque chunk
		if (check_events_and_cancel()) {
			cancelled = 1;
			printf("Calcul DDp2 Pass 2 annulé par l'utilisateur\n");
		}
		if (guiPtr != NULL && progress != NULL && !cancelled) {
			int percent = progressMid + ((chunk_start + chunk_height) * (progressEnd - progressMid)) / ypixel;
			if (percent > progressEnd) percent = progressEnd;
			*progress = percent;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
		}
	}
#else
	for (int chunk_start = 0; chunk_start < ypixel && !cancelled; chunk_start += chunk_height) {
		int chunk_end = (chunk_start + chunk_height < ypixel) ? chunk_start + chunk_height : ypixel;

		for (j = chunk_start; j < chunk_end; j++) {
			complex zPixel_local;
			fractalresult result_local;
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
					double rz, iz;
					if (compteur > 0) {
						rz = (Rez(zup) + Rez(zdown) + Rez(zleft) + Rez(zright)) / compteur;
						iz = (Imz(zup) + Imz(zdown) + Imz(zleft) + Imz(zright)) / compteur;
					} else {
						rz = (Rez(zup) + Rez(zdown) + Rez(zleft) + Rez(zright)) / 4.0;
						iz = (Imz(zup) + Imz(zdown) + Imz(zleft) + Imz(zright)) / 4.0;
					}
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

		// Vérifier annulation et mettre à jour progression après chaque chunk
		if (check_events_and_cancel()) {
			cancelled = 1;
			printf("Calcul DDp2 Pass 2 annulé par l'utilisateur\n");
		}
		if (guiPtr != NULL && progress != NULL && !cancelled) {
			int percent = progressMid + ((chunk_start + chunk_height) * (progressEnd - progressMid)) / ypixel;
			if (percent > progressEnd) percent = progressEnd;
			*progress = percent;
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
		}
	}
#endif

	// Mise à jour finale de la progression
	if (guiPtr != NULL && progress != NULL && !cancelled) {
		*progress = progressEnd;
		SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressEnd, fractalName);
	}
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
		// Vérifier chaque allocation individuellement pour éviter les fuites mémoire
		if (f->cache.fmatrix_cached == NULL || f->cache.cmatrix_cached == NULL) {
			fprintf(stderr, "Erreur allocation mémoire cache\n");
			// Libérer celle qui a réussi si l'autre a échoué
			if (f->cache.fmatrix_cached != NULL) {
				free(f->cache.fmatrix_cached);
				f->cache.fmatrix_cached = NULL;
			}
			if (f->cache.cmatrix_cached != NULL) {
				free(f->cache.cmatrix_cached);
				f->cache.cmatrix_cached = NULL;
			}
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
	
	// Protection contre division par zéro
	if (new_range_x <= 0.0 || old_range_x <= 0.0) {
		return 0;  // Cache non réutilisable
	}
	
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
		f->colorRepeat = 40;  // 40 répétitions par défaut (comme les autres fractales)
	} else {
		f->colorRepeat = 40;  // 40 répétitions pour escape-time
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
	case 24:
	Nebulabrot_def (f);
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
		"Multibrot",         // 23
		"Nebulabrot"         // 24
	};

	if (type >= 0 && type <= 24) {
		return typeNames[type];
	}
	return "Unknown";
}
