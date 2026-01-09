/* EscapeTime.c
fractal handling src code
released under GPL2
Copyleft 2001-2003 VERHILLE Arnaud
*/


#include <stdio.h>
#include <stdlib.h>  // Pour malloc et free
#include <math.h>     // Pour fabs()
#include "SDL_gfxPrimitives.h"
#include "EscapeTime.h"
#include "SDLGUI.h"   // Pour la barre de progression
#ifdef HAVE_GMP
#include "precision_detector.h"
#endif

// Fractal Functions
// *****************



fractal Fractal_Init (int screenW, int screenH, int type) {
	fractal f;
	f.xpixel = screenW;
	f.ypixel = screenH;
	f.type = type;
	f.colorMode = 0;
#ifdef HAVE_GMP
	f.use_gmp = 0;
	f.gmp_precision = 64;
	f.zmatrix_gmp = NULL;
#endif
	Fractal_ChangeType (&f, type);
	f.fmatrix = (int *) malloc ((f.xpixel*f.ypixel)* sizeof (int));
	f.zmatrix = (complex *) malloc ((f.xpixel*f.ypixel)* sizeof (complex));
	f.cmatrix = (color *) malloc ((f.xpixel*f.ypixel)* sizeof (color));
	if (f.fmatrix == NULL || f.zmatrix == NULL || f.cmatrix == NULL) {
		fprintf(stderr, "Erreur allocation mémoire fractale\n");
		exit(1);
	}
#ifdef HAVE_GMP
	precision_update_fractal(&f);
	if (f.use_gmp) {
		// Initialiser les coordonnées GMP
		mpf_init2(f.xmin_gmp, f.gmp_precision);
		mpf_init2(f.xmax_gmp, f.gmp_precision);
		mpf_init2(f.ymin_gmp, f.gmp_precision);
		mpf_init2(f.ymax_gmp, f.gmp_precision);
		// Synchroniser depuis les double
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
	} else {
		// Initialiser quand même pour éviter les problèmes (mais non utilisés)
		mpf_init2(f.xmin_gmp, 64);
		mpf_init2(f.xmax_gmp, 64);
		mpf_init2(f.ymin_gmp, 64);
		mpf_init2(f.ymax_gmp, 64);
	}
#endif
	return f;
}

void Fractal_Destroy (fractal f) {
	free (f.fmatrix);
	free (f.zmatrix);
	free (f.cmatrix);
#ifdef HAVE_GMP
	if (f.zmatrix_gmp != NULL) {
		for (int i = 0; i < f.xpixel * f.ypixel; i++) {
			complex_gmp_clear(&f.zmatrix_gmp[i]);
		}
		free (f.zmatrix_gmp);
	}
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
					complex_gmp_clear(&f->zmatrix_gmp[(f->xpixel*j)+i]);
					f->zmatrix_gmp[(f->xpixel*j)+i] = complex_to_gmp(result.z, prec);
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
	complex_gmp zPixel_gmp;
	mpf_t xg, yg;
	mp_bitcnt_t prec = f->gmp_precision;
	fractalresult result;
	int totalPixels, currentPixel = 0;
	int lastPercent = -1;
	int progressInterval;
	
	printf ("Calculating EscapeTimeFractal Matrix with DDp1 (GMP) ...\n");
	
	// Calculer le nombre total de pixels à traiter (grille 2x2)
	totalPixels = ((f->xpixel + 1) / 2) * ((f->ypixel + 1) / 2);
	progressInterval = totalPixels / 100;
	if (progressInterval < 1) progressInterval = 1;
	
	mpf_init2(xg, prec);
	mpf_init2(yg, prec);
	
	for (j=0; j<f->ypixel; j=j+2) {
		for (i=0; i<f->xpixel; i=i+2) {
			// Mise à jour de la progression
			if (guiPtr != NULL && (currentPixel % progressInterval == 0)) {
				int percent = progressStart + ((currentPixel * (progressEnd - progressStart)) / totalPixels);
				if (percent != lastPercent && percent <= progressEnd) {
					SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
					*progress = percent;
					lastPercent = percent;
				}
			}
			currentPixel++;
			// Calcul des coordonnées en GMP directement depuis les coordonnées GMP
			mpf_t range_x, range_y, step_x, step_y, i_mpf, j_mpf;
			mpf_init2(range_x, prec);
			mpf_init2(range_y, prec);
			mpf_init2(step_x, prec);
			mpf_init2(step_y, prec);
			mpf_init2(i_mpf, prec);
			mpf_init2(j_mpf, prec);
			
			mpf_sub(range_x, f->xmax_gmp, f->xmin_gmp);
			mpf_sub(range_y, f->ymax_gmp, f->ymin_gmp);
			mpf_set_ui(i_mpf, i);
			mpf_set_ui(j_mpf, j);
			mpf_set_ui(step_x, f->xpixel);
			mpf_set_ui(step_y, f->ypixel);
			
			mpf_div(step_x, range_x, step_x);
			mpf_mul(xg, step_x, i_mpf);
			mpf_add(xg, xg, f->xmin_gmp);
			
			mpf_div(step_y, range_y, step_y);
			mpf_mul(yg, step_y, j_mpf);
			mpf_add(yg, yg, f->ymin_gmp);
			
			zPixel_gmp = complex_gmp_make(xg, yg, prec);
			result = FormulaSelector_GMP(*f, zPixel_gmp);
			
			*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
			*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
			if (f->zmatrix_gmp != NULL) {
				complex_gmp_clear(&f->zmatrix_gmp[(f->xpixel*j)+i]);
				f->zmatrix_gmp[(f->xpixel*j)+i] = complex_to_gmp(result.z, prec);
			}
			
			// On complete pour un preview tout en evitant Segfault
			if ((i+1) < f->xpixel) {
				*((f->fmatrix)+((f->xpixel*j)+(i+1))) = result.iteration;
				*((f->zmatrix)+((f->xpixel*j)+(i+1))) = result.z;
				if (f->zmatrix_gmp != NULL) {
					complex_gmp_clear(&f->zmatrix_gmp[(f->xpixel*j)+(i+1)]);
					f->zmatrix_gmp[(f->xpixel*j)+(i+1)] = complex_to_gmp(result.z, prec);
				}
			}
			if ((j+1) < f->ypixel) {
				*((f->fmatrix)+((f->xpixel*(j+1))+i)) = result.iteration;
				*((f->zmatrix)+((f->xpixel*(j+1))+i)) = result.z;
				if (f->zmatrix_gmp != NULL) {
					complex_gmp_clear(&f->zmatrix_gmp[(f->xpixel*(j+1))+i]);
					f->zmatrix_gmp[(f->xpixel*(j+1))+i] = complex_to_gmp(result.z, prec);
				}
			}
			if ( ((i+1) < f->xpixel) && ((j+1) < f->ypixel)) {
				*((f->fmatrix)+((f->xpixel*(j+1))+(i+1))) = result.iteration;
				*((f->zmatrix)+((f->xpixel*(j+1))+(i+1))) = result.z;
				if (f->zmatrix_gmp != NULL) {
					complex_gmp_clear(&f->zmatrix_gmp[(f->xpixel*(j+1))+(i+1)]);
					f->zmatrix_gmp[(f->xpixel*(j+1))+(i+1)] = complex_to_gmp(result.z, prec);
				}
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
}
#endif

void Fractal_CalculateMatrix_DDp1 (fractal* f, SDL_Surface* canvas, void* guiPtr, int* progress, int progressStart, int progressEnd, const char* fractalName) {
#ifdef HAVE_GMP
	if (f->use_gmp) {
		Fractal_CalculateMatrix_DDp1_GMP(f, canvas, guiPtr, progress, progressStart, progressEnd, fractalName);
		return;
	}
#endif
	int i, j;
	double xg, yg;
	complex zPixel;
	fractalresult result;
	int totalPixels, currentPixel = 0;
	int lastPercent = -1;
	int progressInterval;
	
	printf ("Calculating EscapeTimeFractal Matrix with DDp1 ...\n");
	
	// Calculer le nombre total de pixels à traiter (grille 2x2)
	totalPixels = ((f->xpixel + 1) / 2) * ((f->ypixel + 1) / 2);
	progressInterval = totalPixels / 100;
	if (progressInterval < 1) progressInterval = 1;
	
	for (j=0; j<f->ypixel; j=j+2) {
		for (i=0; i<f->xpixel; i=i+2) {
			// Mise à jour de la progression
			if (guiPtr != NULL && (currentPixel % progressInterval == 0)) {
				int percent = progressStart + ((currentPixel * (progressEnd - progressStart)) / totalPixels);
				if (percent != lastPercent && percent <= progressEnd) {
					SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
					*progress = percent;
					lastPercent = percent;
				}
			}
			currentPixel++;
			
			xg = ((f->xmax-f->xmin)/f->xpixel)*i + f->xmin;
			yg = ((f->ymax-f->ymin)/f->ypixel)*j + f->ymin;
			zPixel = MakeComplex (xg, yg);
			
			result = FormulaSelector (*f,zPixel);
			*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
			*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
			
			// On complete pour un preview tout en evitant Segfault
			if ((i+1) < f->xpixel) {
				*((f->fmatrix)+((f->xpixel*j)+(i+1))) = result.iteration;
				*((f->zmatrix)+((f->xpixel*j)+(i+1))) = result.z;
			}
			if ((j+1) < f->ypixel) {
				*((f->fmatrix)+((f->xpixel*(j+1))+i)) = result.iteration;
				*((f->zmatrix)+((f->xpixel*(j+1))+i)) = result.z;
			}
			if ( ((i+1) < f->xpixel) && ((j+1) < f->ypixel)) {
				*((f->fmatrix)+((f->xpixel*(j+1))+(i+1))) = result.iteration;
				*((f->zmatrix)+((f->xpixel*(j+1))+(i+1))) = result.z;
			}
		}
	}
}

#ifdef HAVE_GMP
static void Fractal_CalculateMatrix_DDp2_GMP (fractal* f, SDL_Surface* canvas, void* guiPtr, int* progress, int progressStart, int progressEnd, const char* fractalName) {
	int i,starti, j, compteur;
	int up, down, left, right;
	complex zup, zdown, zleft, zright;
	complex_gmp zPixel_gmp;
	mpf_t xg, yg;
	mp_bitcnt_t prec = f->gmp_precision;
	fractalresult result;
	int totalPixels, currentPixel = 0;
	int lastPercent = -1;
	int progressInterval;
	int pass1Pixels, pass2Pixels;
	
	printf ("Calculating EscapeTimeFractal Matrix with DDp2 (GMP) ...\n");
	
	// Calculer le nombre total de pixels à traiter
	// Pass 1: pixels impairs (j=1, i=1, step 2)
	pass1Pixels = ((f->xpixel - 1) / 2) * ((f->ypixel - 1) / 2);
	// Pass 2: pixels restants (divergence detection)
	pass2Pixels = (f->xpixel * f->ypixel) / 2;
	totalPixels = pass1Pixels + pass2Pixels;
	progressInterval = totalPixels / 100;
	if (progressInterval < 1) progressInterval = 1;
	
	mpf_init2(xg, prec);
	mpf_init2(yg, prec);
	
	// Pass 1: Calcul initial
	for (j=1; j<f->ypixel; j=j+2) {
		for (i=1; i<f->xpixel; i=i+2) {
			// Mise à jour de la progression
			if (guiPtr != NULL && (currentPixel % progressInterval == 0)) {
				int percent = progressStart + ((currentPixel * (progressEnd - progressStart)) / totalPixels);
				if (percent != lastPercent && percent <= progressEnd) {
					SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
					*progress = percent;
					lastPercent = percent;
				}
			}
			currentPixel++;
			// Calcul des coordonnées en GMP directement depuis les coordonnées GMP
			mpf_t range_x, range_y, step_x, step_y, i_mpf, j_mpf;
			mpf_init2(range_x, prec);
			mpf_init2(range_y, prec);
			mpf_init2(step_x, prec);
			mpf_init2(step_y, prec);
			mpf_init2(i_mpf, prec);
			mpf_init2(j_mpf, prec);
			
			mpf_sub(range_x, f->xmax_gmp, f->xmin_gmp);
			mpf_sub(range_y, f->ymax_gmp, f->ymin_gmp);
			mpf_set_ui(i_mpf, i);
			mpf_set_ui(j_mpf, j);
			mpf_set_ui(step_x, f->xpixel);
			mpf_set_ui(step_y, f->ypixel);
			
			mpf_div(step_x, range_x, step_x);
			mpf_mul(xg, step_x, i_mpf);
			mpf_add(xg, xg, f->xmin_gmp);
			
			mpf_div(step_y, range_y, step_y);
			mpf_mul(yg, step_y, j_mpf);
			mpf_add(yg, yg, f->ymin_gmp);
			
			zPixel_gmp = complex_gmp_make(xg, yg, prec);
			result = FormulaSelector_GMP(*f, zPixel_gmp);
			
			*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
			*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
			if (f->zmatrix_gmp != NULL) {
				complex_gmp_clear(&f->zmatrix_gmp[(f->xpixel*j)+i]);
				f->zmatrix_gmp[(f->xpixel*j)+i] = complex_to_gmp(result.z, prec);
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

	// Pass 2: Dernier passage, cette fois on compare pour savoir si l'on peut se
	// permettre d'eviter le calcul, c'est la Divergence Detection (DD)
	
	for (j=0; j<f->ypixel; j++) {
		if (fmod (j,2) == 0) {starti=1;} else {starti=0;} // Lignes paires et impaires
		for (i=starti; i<f->xpixel; i=i+2) {
			// Mise à jour de la progression
			if (guiPtr != NULL && (currentPixel % progressInterval == 0)) {
				int percent = progressStart + ((currentPixel * (progressEnd - progressStart)) / totalPixels);
				if (percent != lastPercent && percent <= progressEnd) {
					SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
					*progress = percent;
					lastPercent = percent;
				}
			}
			currentPixel++;
			
			// Verifions les bords des matrices pour eviter Seg Fault
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
      
      /* Si les points adjacents ont meme valeur de divergence alors */
      if ((up == right) && (right == down) && (down == left) && (left == up)) {
		  
		  // OK, on ne calcule pas les Iterations
		  *((f->fmatrix)+((f->xpixel*j)+i)) = up;
		  
		  // Et on calcule une moyenne pour la valeur de z
		  double rz, iz;
		  rz = (Rez(zup)+Rez(zdown)+Rez(zleft)+Rez(zright))/compteur;
		  iz = (Imz(zup)+Imz(zdown)+Imz(zleft)+Imz(zright))/compteur;
		  *((f->zmatrix)+((f->xpixel*j)+i)) = MakeComplex (rz, iz);
		  if (f->zmatrix_gmp != NULL) {
			  complex_gmp_clear(&f->zmatrix_gmp[(f->xpixel*j)+i]);
			  f->zmatrix_gmp[(f->xpixel*j)+i] = complex_to_gmp(*((f->zmatrix)+((f->xpixel*j)+i)), prec);
		  }
		  
      } else { // Sinon, je calcule
		  // Calcul des coordonnées en GMP directement depuis les coordonnées GMP
		  mpf_t range_x, range_y, step_x, step_y, i_mpf, j_mpf;
		  mpf_init2(range_x, prec);
		  mpf_init2(range_y, prec);
		  mpf_init2(step_x, prec);
		  mpf_init2(step_y, prec);
		  mpf_init2(i_mpf, prec);
		  mpf_init2(j_mpf, prec);
		  
		  mpf_sub(range_x, f->xmax_gmp, f->xmin_gmp);
		  mpf_sub(range_y, f->ymax_gmp, f->ymin_gmp);
		  mpf_set_ui(i_mpf, i);
		  mpf_set_ui(j_mpf, j);
		  mpf_set_ui(step_x, f->xpixel);
		  mpf_set_ui(step_y, f->ypixel);
		  
		  mpf_div(step_x, range_x, step_x);
		  mpf_mul(xg, step_x, i_mpf);
		  mpf_add(xg, xg, f->xmin_gmp);
		  
		  mpf_div(step_y, range_y, step_y);
		  mpf_mul(yg, step_y, j_mpf);
		  mpf_add(yg, yg, f->ymin_gmp);
		  
		  zPixel_gmp = complex_gmp_make(xg, yg, prec);
		  result = FormulaSelector_GMP(*f, zPixel_gmp);
		  
		  *((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
		  *((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
		  if (f->zmatrix_gmp != NULL) {
			  complex_gmp_clear(&f->zmatrix_gmp[(f->xpixel*j)+i]);
			  f->zmatrix_gmp[(f->xpixel*j)+i] = complex_to_gmp(result.z, prec);
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
	}
	
	mpf_clear(xg);
	mpf_clear(yg);
}
#endif

void Fractal_CalculateMatrix_DDp2 (fractal* f, SDL_Surface* canvas, void* guiPtr, int* progress, int progressStart, int progressEnd, const char* fractalName) {
#ifdef HAVE_GMP
	if (f->use_gmp) {
		Fractal_CalculateMatrix_DDp2_GMP(f, canvas, guiPtr, progress, progressStart, progressEnd, fractalName);
		return;
	}
#endif
	int i,starti, j, compteur;
	int up, down, left, right;
	complex zup, zdown, zleft, zright;
	double xg, yg, rz, iz;
	complex zPixel;
	fractalresult result;
	int totalPixels, currentPixel = 0;
	int lastPercent = -1;
	int progressInterval;
	int pass1Pixels, pass2Pixels;
	
	printf ("Calculating EscapeTimeFractal Matrix with DDp2 ...\n");
	
	// Calculer le nombre total de pixels à traiter
	// Pass 1: pixels impairs (j=1, i=1, step 2)
	pass1Pixels = ((f->xpixel - 1) / 2) * ((f->ypixel - 1) / 2);
	// Pass 2: pixels restants (divergence detection)
	pass2Pixels = (f->xpixel * f->ypixel) / 2;
	totalPixels = pass1Pixels + pass2Pixels;
	progressInterval = totalPixels / 100;
	if (progressInterval < 1) progressInterval = 1;
	
	// Pass 1: Calcul initial
	for (j=1; j<f->ypixel; j=j+2) {
		for (i=1; i<f->xpixel; i=i+2) {
			// Mise à jour de la progression
			if (guiPtr != NULL && (currentPixel % progressInterval == 0)) {
				int percent = progressStart + ((currentPixel * (progressEnd - progressStart)) / totalPixels);
				if (percent != lastPercent && percent <= progressEnd) {
					SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
					*progress = percent;
					lastPercent = percent;
				}
			}
			currentPixel++;
			
			xg = ((f->xmax-f->xmin)/f->xpixel)*i + f->xmin;
			yg = ((f->ymax-f->ymin)/f->ypixel)*j + f->ymin;
			zPixel = MakeComplex (xg, yg);
			
			result = FormulaSelector (*f,zPixel);
			*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
			*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
		}
	}

	// Pass 2: Dernier passage, cette fois on compare pour savoir si l'on peut se
	// permettre d'eviter le calcul, c'est la Divergence Detection (DD)
	
	for (j=0; j<f->ypixel; j++) {
		if (fmod (j,2) == 0) {starti=1;} else {starti=0;} // Lignes paires et impaires
		for (i=starti; i<f->xpixel; i=i+2) {
			// Mise à jour de la progression
			if (guiPtr != NULL && (currentPixel % progressInterval == 0)) {
				int percent = progressStart + ((currentPixel * (progressEnd - progressStart)) / totalPixels);
				if (percent != lastPercent && percent <= progressEnd) {
					SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
					*progress = percent;
					lastPercent = percent;
				}
			}
			currentPixel++;
			
			// Verifions les bords des matrices pour eviter Seg Fault
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
      
      /* Si les points adjacents ont meme valeur de divergence alors */
      if ((up == right) && (right == down) && (down == left) && (left == up)) {
		  
		  // OK, on ne calcule pas les Iterations
		  *((f->fmatrix)+((f->xpixel*j)+i)) = up;
		  
		  // Et on calcule une moyenne pour la valeur de z
		  rz = (Rez(zup)+Rez(zdown)+Rez(zleft)+Rez(zright))/compteur;
		  iz = (Imz(zup)+Imz(zdown)+Imz(zleft)+Imz(zright))/compteur;
		  *((f->zmatrix)+((f->xpixel*j)+i)) = MakeComplex (rz, iz);
		  
      } else { // Sinon, je calcule
		  double xg, yg; 
		  xg = ((f->xmax-f->xmin)/f->xpixel)*i + f->xmin;
		  yg = ((f->ymax-f->ymin)/f->ypixel)*j + f->ymin;
		  zPixel = MakeComplex (xg, yg);
		  
		  result = FormulaSelector (*f,zPixel);
		  *((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
		  *((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
      }
		}
	}
}



// Color Formulae Functions
// *************************

/* Un mode de couleur classique */
void FractalColorNormal (fractal* f) {
	int i, j;
	int iteration, greyvalue;
	color c;

	for (j = 0; j < f->ypixel; j++) {
		for (i=0; i < f->xpixel; i++) {
			iteration = *((f->fmatrix)+((f->xpixel*j)+i));
			greyvalue = 255 - (iteration*255)/f->iterationMax;
			c.g = (int) fmod (greyvalue*3, 510);
			if (c.g > 255) { c.g = 255 - (c.g - 255); }
			c.r = (int) fmod (greyvalue, 510);
			if (c.r > 255) { c.r = 255 - (c.r - 255); }
			c.b = (int) fmod (greyvalue*2, 510);
			if (c.b > 255) { c.b = 255 - (c.b - 255); }
			*((f->cmatrix)+((f->xpixel*j)+i)) = c;
		}
	}
}

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

/* Un mode de couleur monochrome */
void FractalColorMonochrome (fractal* f) {
	int i, j;
	int iteration, greyvalue;
	color c;

	for (j = 0; j < f->ypixel; j++) {
		for (i=0; i < f->xpixel; i++) {
			iteration = *((f->fmatrix)+((f->xpixel*j)+i));
			greyvalue = 255 - (iteration*255)/f->iterationMax;
			c.g = greyvalue;
			c.r = greyvalue;
			c.b = greyvalue;
			*((f->cmatrix)+((f->xpixel*j)+i)) = c;
		}
	}
}

/* Palette Fire : noir -> rouge -> jaune -> blanc */
void FractalColorFire (fractal* f) {
	int i, j;
	int iteration;
	double t;
	color c;

	for (j = 0; j < f->ypixel; j++) {
		for (i=0; i < f->xpixel; i++) {
			iteration = *((f->fmatrix)+((f->xpixel*j)+i));
			t = (double)iteration / f->iterationMax;
			if (t < 0.33) {
				c.r = (int)(t * 3 * 255);
				c.g = 0;
				c.b = 0;
			} else if (t < 0.66) {
				c.r = 255;
				c.g = (int)((t - 0.33) * 3 * 255);
				c.b = 0;
			} else {
				c.r = 255;
				c.g = 255;
				c.b = (int)((t - 0.66) * 3 * 255);
			}
			*((f->cmatrix)+((f->xpixel*j)+i)) = c;
		}
	}
}

/* Palette Ocean : noir -> bleu -> cyan -> blanc */
void FractalColorOcean (fractal* f) {
	int i, j;
	int iteration;
	double t;
	color c;

	for (j = 0; j < f->ypixel; j++) {
		for (i=0; i < f->xpixel; i++) {
			iteration = *((f->fmatrix)+((f->xpixel*j)+i));
			t = (double)iteration / f->iterationMax;
			if (t < 0.33) {
				c.r = 0;
				c.g = 0;
				c.b = (int)(t * 3 * 255);
			} else if (t < 0.66) {
				c.r = 0;
				c.g = (int)((t - 0.33) * 3 * 255);
				c.b = 255;
			} else {
				c.r = (int)((t - 0.66) * 3 * 255);
				c.g = 255;
				c.b = 255;
			}
			*((f->cmatrix)+((f->xpixel*j)+i)) = c;
		}
	}
}

/* Calcul smooth iteration pour coloring continu */
double Fractal_SmoothIteration(fractal* f, int i, int j) {
	int iteration = *((f->fmatrix)+((f->xpixel*j)+i));
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

/* Palette Rainbow : arc-en-ciel HSV avec smooth coloring */
void FractalColorRainbow(fractal* f) {
	int i, j;
	double t;

	for (j = 0; j < f->ypixel; j++) {
		for (i = 0; i < f->xpixel; i++) {
			t = Fractal_SmoothIteration(f, i, j);
			// Cycle de couleurs (multiple de 360 pour répétition)
			*((f->cmatrix)+((f->xpixel*j)+i)) = HSVtoRGB(t * 360 * 3, 0.8, 1.0 - t * 0.3);
		}
	}
}

/* Palette Smooth Fire : version fluide de Fire */
void FractalColorSmoothFire(fractal* f) {
	int i, j;
	double t;
	color c;

	for (j = 0; j < f->ypixel; j++) {
		for (i = 0; i < f->xpixel; i++) {
			t = Fractal_SmoothIteration(f, i, j);
			if (t < 0.33) {
				c.r = (int)(t * 3 * 255);
				c.g = 0;
				c.b = 0;
			} else if (t < 0.66) {
				c.r = 255;
				c.g = (int)((t - 0.33) * 3 * 255);
				c.b = 0;
			} else {
				c.r = 255;
				c.g = 255;
				c.b = (int)((t - 0.66) * 3 * 255);
			}
			*((f->cmatrix)+((f->xpixel*j)+i)) = c;
		}
	}
}

/* Palette Smooth Ocean : version fluide de Ocean */
void FractalColorSmoothOcean(fractal* f) {
	int i, j;
	double t;
	color c;

	for (j = 0; j < f->ypixel; j++) {
		for (i = 0; i < f->xpixel; i++) {
			t = Fractal_SmoothIteration(f, i, j);
			if (t < 0.33) {
				c.r = 0;
				c.g = 0;
				c.b = (int)(t * 3 * 255);
			} else if (t < 0.66) {
				c.r = 0;
				c.g = (int)((t - 0.33) * 3 * 255);
				c.b = 255;
			} else {
				c.r = (int)((t - 0.66) * 3 * 255);
				c.g = 255;
				c.b = 255;
			}
			*((f->cmatrix)+((f->xpixel*j)+i)) = c;
		}
	}
}

/* Selection du mode de couleur */
//**********************************

void Fractal_CalculateColorMatrix (fractal* f, SDL_Surface* canvas, void* guiPtr, int* progress, int progressStart, int progressEnd) {
	int totalPixels = f->xpixel * f->ypixel;
	int currentPixel = 0;
	int lastPercent = -1;
	int progressInterval = totalPixels / 100;
	if (progressInterval < 1) progressInterval = 1;
	const char* fractalName = Fractal_GetTypeName(f->type);
	
	switch (f->colorMode) {
		case 0: FractalColorNormal(f); break;
		case 1: FractalColorMonochrome(f); break;
		case 2: FractalColorFire(f); break;
		case 3: FractalColorOcean(f); break;
		case 4: FractalColorRainbow(f); break;
		case 5: FractalColorSmoothFire(f); break;
		case 6: FractalColorSmoothOcean(f); break;
		default: FractalColorNormal(f);
	}
	
	// Mise à jour de la progression après colorisation
	// On simule la progression en parcourant les pixels
	if (guiPtr != NULL) {
		for (int j = 0; j < f->ypixel; j++) {
			for (int i = 0; i < f->xpixel; i++) {
				if (currentPixel % progressInterval == 0) {
					int percent = progressStart + ((currentPixel * (progressEnd - progressStart)) / totalPixels);
					if (percent != lastPercent && percent <= progressEnd) {
						SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, fractalName);
						*progress = percent;
						lastPercent = percent;
					}
				}
				currentPixel++;
			}
		}
		// Assurer que la progression atteint progressEnd
		if (*progress < progressEnd) {
			SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, progressEnd, fractalName);
			*progress = progressEnd;
		}
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



Uint32 Fractal_Draw (SDL_Surface *canvas, fractal myfractal,int decalageX,int decalageY, void* guiPtr) {

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
	Sierpinski_def (f);
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
	default:
	Mendelbrot_def (f);
	}

#ifdef HAVE_GMP
	// Mettre à jour la précision après le changement de type
	precision_update_fractal(f);
	
	// Synchroniser les coordonnées double → GMP
	// Les coordonnées GMP sont toujours initialisées dans Fractal_Init
	// On vérifie juste si la précision doit être mise à jour
	if (f->use_gmp) {
		// Les coordonnées GMP sont toujours initialisées dans Fractal_Init
		// On met juste à jour la précision si nécessaire
		mp_bitcnt_t current_prec = mpf_get_prec(f->xmin_gmp);
		if (current_prec != f->gmp_precision) {
			// Mettre à jour la précision si nécessaire
			mpf_set_prec(f->xmin_gmp, f->gmp_precision);
			mpf_set_prec(f->xmax_gmp, f->gmp_precision);
			mpf_set_prec(f->ymin_gmp, f->gmp_precision);
			mpf_set_prec(f->ymax_gmp, f->gmp_precision);
		}
		// Synchroniser les valeurs
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
		"Sierpinski", // 8
		"Barnsley J", // 9
		"Barnsley M", // 10
		"Magnet J",   // 11
		"Magnet M",   // 12
		"Burning Ship", // 13
		"Tricorn",    // 14
		"Mandelbulb", // 15
		"Buddhabrot"  // 16
	};
	
	if (type >= 0 && type <= 16) {
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
	f->iterationMax= 150;
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
	f->iterationMax= 100;
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
	f->iterationMax= 100;
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
	f->iterationMax= 30;
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
	f->iterationMax= 150;
	f->zoomfactor = 8;
}




/* Calcul de la Sierpinski */
fractalresult Sierpinski_Iteration (fractal f, complex zPixel) {
	/*
	 * z(n+1) = (2*x, 2*y-1) if y>0.5;
	 * else (2*x-1,2*y) if x> 0.5;
	 * else (2*x,2*y)
	 */
	int i=0;
	complex z;
	fractalresult result;
	z = zPixel;
	do {
		i++;
		if (Imz(z) > 0.5) {
			z = MakeComplex (2*Rez(z),(2*Imz(z))-1);
		} else {
			if (Rez(z) > 0.5) {
				z = MakeComplex ((2*Rez(z))-1,2*Imz(z));
			} else {
				z = MakeComplex (2*Rez(z),2*Imz(z));
			}
		}
	} while ((i < f.iterationMax) && ( Magz (z) < f.bailout));
	result.iteration = i;
	result.z = z;
	return result;
}
void Sierpinski_def (fractal* f) {
	/* Valeurs de base pour la Sierpinski */
	f->xmin = -1.333402669;
	f->xmax = 2.133402669;
	f->ymin = -0.9000520021;  // Corrigé : ymin doit être < ymax
	f->ymax = 1.700052002;     // Corrigé : ymax doit être > ymin
	f->bailout= 127;
	f->iterationMax= 149;
	f->zoomfactor = 4;
}

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
	f->iterationMax= 50;
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
	f->iterationMax= 150;
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
	f->iterationMax= 150;
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
	f->iterationMax= 150;
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
	f->iterationMax= 150;
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
	f->iterationMax= 150;
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
	f->iterationMax= 150;
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
	f->iterationMax = 100;  // Plus d'itérations pour de meilleurs détails
	f->zoomfactor = 4;
}

/* Fonction de rendu spéciale pour Buddhabrot */
Uint32 Buddhabrot_Draw (SDL_Surface *canvas, fractal* f, int decalageX, int decalageY, void* guiPtr) {
	int i, j, iter;
	int px, py;
	double xg, yg;
	complex c, z;
	Uint32 time;
	int numSamples;
	double maxDensity;
	color col;
	int lastPercent = -1;
	int progressInterval;

	// Tableaux pour stocker la trajectoire
	double *trajX, *trajY;

	time = SDL_GetTicks();

	printf("Calculating Buddhabrot (density algorithm)...\n");

	// Allouer mémoire pour trajectoire
	trajX = (double*) malloc(f->iterationMax * sizeof(double));
	trajY = (double*) malloc(f->iterationMax * sizeof(double));
	if (trajX == NULL || trajY == NULL) {
		fprintf(stderr, "Erreur allocation mémoire trajectoire\n");
		return 0;
	}

	// Réinitialiser la matrice de densité (on utilise fmatrix)
	for (i = 0; i < f->xpixel * f->ypixel; i++) {
		f->fmatrix[i] = 0;
	}

	// Nombre d'échantillons proportionnel à la surface
	numSamples = f->xpixel * f->ypixel * 50;

	// Intervalle de mise à jour de la progression (tous les 1%)
	progressInterval = numSamples / 100;
	if (progressInterval < 1) progressInterval = 1;

	// Echantillonnage aléatoire
	srand(42);  // Seed fixe pour reproductibilité

	for (int sample = 0; sample < numSamples; sample++) {
		// Mise à jour de la progression
		if (guiPtr != NULL && (sample % progressInterval == 0)) {
			int percent = (sample * 100) / numSamples;
			if (percent != lastPercent) {
				SDLGUI_StateBar_Progress(canvas, (gui*)guiPtr, percent, "Buddhabrot");
				lastPercent = percent;
			}
		}

		// Générer un point c aléatoire dans la zone
		xg = ((double)rand() / RAND_MAX) * (f->xmax - f->xmin) + f->xmin;
		yg = ((double)rand() / RAND_MAX) * (f->ymax - f->ymin) + f->ymin;
		c = MakeComplex(xg, yg);
		z = ZeroSetofComplex();

		// Itérer et stocker la trajectoire
		int escaped = 0;
		for (iter = 0; iter < f->iterationMax; iter++) {
			complex zTemp = Mulz(z, z);
			z = Addz(zTemp, c);

			trajX[iter] = Rez(z);
			trajY[iter] = Imz(z);

			if (Magz(z) > f->bailout) {
				escaped = 1;
				break;
			}
		}

		// Si le point s'échappe, tracer sa trajectoire
		if (escaped) {
			for (i = 0; i < iter; i++) {
				// Convertir coordonnées complexes en pixels
				px = (int)((trajX[i] - f->xmin) / (f->xmax - f->xmin) * f->xpixel);
				py = (int)((trajY[i] - f->ymin) / (f->ymax - f->ymin) * f->ypixel);

				// Vérifier les limites et incrémenter la densité
				if (px >= 0 && px < f->xpixel && py >= 0 && py < f->ypixel) {
					f->fmatrix[py * f->xpixel + px]++;
				}
			}
		}
	}

	// Trouver la densité maximale pour normalisation
	maxDensity = 1;
	for (i = 0; i < f->xpixel * f->ypixel; i++) {
		if (f->fmatrix[i] > maxDensity) {
			maxDensity = f->fmatrix[i];
		}
	}

	// Convertir densité en couleurs et afficher
	for (j = 0; j < f->ypixel; j++) {
		for (i = 0; i < f->xpixel; i++) {
			double density = (double)f->fmatrix[j * f->xpixel + i];
			// Utiliser une échelle logarithmique pour mieux voir les détails
			double normalized = log(1 + density) / log(1 + maxDensity);

			// Appliquer la palette selon colorMode
			switch (f->colorMode) {
				case 1: // Mono
					col.r = col.g = col.b = (int)(normalized * 255);
					break;
				case 2: // Fire
					if (normalized < 0.33) {
						col.r = (int)(normalized * 3 * 255);
						col.g = 0;
						col.b = 0;
					} else if (normalized < 0.66) {
						col.r = 255;
						col.g = (int)((normalized - 0.33) * 3 * 255);
						col.b = 0;
					} else {
						col.r = 255;
						col.g = 255;
						col.b = (int)((normalized - 0.66) * 3 * 255);
					}
					break;
				case 3: // Ocean
					if (normalized < 0.33) {
						col.r = 0;
						col.g = 0;
						col.b = (int)(normalized * 3 * 255);
					} else if (normalized < 0.66) {
						col.r = 0;
						col.g = (int)((normalized - 0.33) * 3 * 255);
						col.b = 255;
					} else {
						col.r = (int)((normalized - 0.66) * 3 * 255);
						col.g = 255;
						col.b = 255;
					}
					break;
				default: // Normal (violet/bleu style Buddhabrot classique)
					col.r = (int)(normalized * 180);
					col.g = (int)(normalized * 100);
					col.b = (int)(normalized * 255);
					break;
			}

			pixelRGBA(canvas, (Sint16)(i + decalageX), (Sint16)(j + decalageY),
			          col.r, col.g, col.b, 255);
		}
	}

	SDL_UpdateRect(canvas, 0, 0, canvas->w, canvas->h);

	free(trajX);
	free(trajY);

	time = SDL_GetTicks() - time;
	printf("Buddhabrot rendered in %d ms (%d samples)\n", time, numSamples);
	return time;
}

#ifdef HAVE_GMP
// ******************************************************
// Versions GMP des fonctions d'itération
// ******************************************************

fractalresult Mendelbrot_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp;
	mpf_t mag, bailout_mpf;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	
	z = complex_gmp_copy(seed_gmp, prec);
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	
	do {
		i++;
		zTemp = complex_gmp_mul(z, z, prec);
		z = complex_gmp_add(zTemp, zPixel, prec);
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zTemp);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	
	return result;
}

fractalresult Julia_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, seed_gmp;
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
		zTemp = complex_gmp_mul(z, z, prec);
		z = complex_gmp_add(zTemp, seed_gmp, prec);
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zTemp);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	
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
		z = complex_gmp_mul(seed_gmp, sin_z, prec);
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
		zTemp = complex_gmp_mul(z, z, prec);
		complex_gmp_get_real(z_real, zTemp);
		mpf_add(z_real, z_real, p1);
		complex_gmp_set_real(&zTemp, z_real);
		
		zpTemp = complex_gmp_scalar_mul(p2, y, prec);
		zTemp = complex_gmp_add(zTemp, zpTemp, prec);
		
		y = z;
		z = zTemp;
		
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

fractalresult Sierpinski_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, z_new;
	mpf_t mag, bailout_mpf, z_real, z_imag, two, one_mpf, half_mpf;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	z = complex_gmp_copy(zPixel, prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(z_real, prec);
	mpf_init2(z_imag, prec);
	mpf_init2(two, prec);
	mpf_init2(one_mpf, prec);
	mpf_init2(half_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	mpf_set_ui(two, 2);
	mpf_set_ui(one_mpf, 1);
	mpf_set_d(half_mpf, 0.5);
	
	do {
		i++;
		complex_gmp_get_real(z_real, z);
		complex_gmp_get_imag(z_imag, z);
		
		if (mpf_cmp(z_imag, half_mpf) > 0) {
			// z = (2*x, 2*y-1)
			mpf_t temp_real, temp_imag;
			mpf_init2(temp_real, prec);
			mpf_init2(temp_imag, prec);
			mpf_mul(temp_real, z_real, two);
			mpf_mul(temp_imag, z_imag, two);
			mpf_sub(temp_imag, temp_imag, one_mpf);
			z_new = complex_gmp_make(temp_real, temp_imag, prec);
			mpf_clear(temp_real);
			mpf_clear(temp_imag);
		} else if (mpf_cmp(z_real, half_mpf) > 0) {
			// z = (2*x-1, 2*y)
			mpf_t temp_real, temp_imag;
			mpf_init2(temp_real, prec);
			mpf_init2(temp_imag, prec);
			mpf_mul(temp_real, z_real, two);
			mpf_sub(temp_real, temp_real, one_mpf);
			mpf_mul(temp_imag, z_imag, two);
			z_new = complex_gmp_make(temp_real, temp_imag, prec);
			mpf_clear(temp_real);
			mpf_clear(temp_imag);
		} else {
			// z = (2*x, 2*y)
			mpf_t temp_real, temp_imag;
			mpf_init2(temp_real, prec);
			mpf_init2(temp_imag, prec);
			mpf_mul(temp_real, z_real, two);
			mpf_mul(temp_imag, z_imag, two);
			z_new = complex_gmp_make(temp_real, temp_imag, prec);
			mpf_clear(temp_real);
			mpf_clear(temp_imag);
		}
		
		complex_gmp_clear(&z);
		z = z_new;
		
		complex_gmp_mag(mag, z);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(z_real);
	mpf_clear(z_imag);
	mpf_clear(two);
	mpf_clear(one_mpf);
	mpf_clear(half_mpf);
	
	return result;
}

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
		z = complex_gmp_mul(zTemp, seed_gmp, prec);
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
		z = complex_gmp_mul(zTemp, c, prec);
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
	mpf_t mag, bailout_mpf, z_real, z_imag, abs_real, abs_imag;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(z_real, prec);
	mpf_init2(z_imag, prec);
	mpf_init2(abs_real, prec);
	mpf_init2(abs_imag, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	
	do {
		i++;
		complex_gmp_get_real(z_real, z);
		complex_gmp_get_imag(z_imag, z);
		
		// |Re(z)| et |Im(z)|
		mpf_abs(abs_real, z_real);
		mpf_abs(abs_imag, z_imag);
		
		zAbs = complex_gmp_make(abs_real, abs_imag, prec);
		zTemp = complex_gmp_mul(zAbs, zAbs, prec);
		z = complex_gmp_add(zTemp, zPixel, prec);
		
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zTemp);
		complex_gmp_clear(&zAbs);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(z_real);
	mpf_clear(z_imag);
	mpf_clear(abs_real);
	mpf_clear(abs_imag);
	
	return result;
}

fractalresult Tricorn_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp, zConj;
	mpf_t mag, bailout_mpf, z_imag;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_init2(z_imag, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	
	do {
		i++;
		// Conjugué: (x, -y)
		mpf_t z_real_conj;
		mpf_init2(z_real_conj, prec);
		complex_gmp_get_real(z_real_conj, z);
		complex_gmp_get_imag(z_imag, z);
		mpf_neg(z_imag, z_imag);
		zConj = complex_gmp_make(z_real_conj, z_imag, prec);
		mpf_clear(z_real_conj);
		
		zTemp = complex_gmp_mul(zConj, zConj, prec);
		z = complex_gmp_add(zTemp, zPixel, prec);
		
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zTemp);
		complex_gmp_clear(&zConj);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	mpf_clear(z_imag);
	
	return result;
}

fractalresult Mandelbulb_Iteration_GMP (fractal f, complex_gmp zPixel) {
	int i = 0;
	complex_gmp z, zTemp;
	mpf_t mag, bailout_mpf;
	fractalresult result;
	
	mp_bitcnt_t prec = f.gmp_precision;
	complex_gmp seed_gmp = complex_to_gmp(f.seed, prec);
	z = complex_gmp_copy(seed_gmp, prec);
	
	mpf_init2(mag, prec);
	mpf_init2(bailout_mpf, prec);
	mpf_set_ui(bailout_mpf, f.bailout);
	
	do {
		i++;
		// z^8 par multiplications successives
		zTemp = complex_gmp_mul(z, z, prec);      // z^2
		zTemp = complex_gmp_mul(zTemp, zTemp, prec); // z^4
		zTemp = complex_gmp_mul(zTemp, zTemp, prec); // z^8
		z = complex_gmp_add(zTemp, zPixel, prec);
		
		complex_gmp_mag(mag, z);
		complex_gmp_clear(&zTemp);
	} while ((i < f.iterationMax) && (mpf_cmp(mag, bailout_mpf) < 0));
	
	result.iteration = i;
	result.z = gmp_to_complex(z);
	
	complex_gmp_clear(&z);
	complex_gmp_clear(&seed_gmp);
	mpf_clear(mag);
	mpf_clear(bailout_mpf);
	
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
		N = complex_gmp_add(z_sq, seed_minus_one, prec);
		N = complex_gmp_mul(N, N, prec);
		
		// Q = 2*z + seed - 2
		complex_gmp two_z = complex_gmp_mul(two, z, prec);
		Q = complex_gmp_add(two_z, seed_minus_two, prec);
		
		// z = N / Q
		z = complex_gmp_div(N, Q, prec);
		
		complex_gmp_mag(mag, z);
		
		complex_gmp_clear(&z_sq);
		complex_gmp_clear(&two_z);
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
		N = complex_gmp_add(z_sq, c_minus_one, prec);
		N = complex_gmp_mul(N, N, prec);
		
		// Q = 2*z + c - 2
		complex_gmp two_z = complex_gmp_mul(two, z, prec);
		Q = complex_gmp_add(two_z, c_minus_two, prec);
		
		// z = N / Q
		z = complex_gmp_div(N, Q, prec);
		
		complex_gmp_mag(mag, z);
		
		complex_gmp_clear(&z_sq);
		complex_gmp_clear(&two_z);
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
		r = Sierpinski_Iteration_GMP(f, zPixel);
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
	r=Sierpinski_Iteration (f,zPixel);
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

	}

		return r;
}
