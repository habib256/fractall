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

// Fractal Functions
// *****************



fractal Fractal_Init (int screenW, int screenH, int type) {
	fractal f;
	f.xpixel = screenW;
	f.ypixel = screenH;
	f.type = type;
	f.colorMode = 0;
	Fractal_ChangeType (&f, type);
	f.fmatrix = (int *) malloc ((f.xpixel*f.ypixel)* sizeof (int));
	f.zmatrix = (complex *) malloc ((f.xpixel*f.ypixel)* sizeof (complex));
	f.cmatrix = (color *) malloc ((f.xpixel*f.ypixel)* sizeof (color));
	if (f.fmatrix == NULL || f.zmatrix == NULL || f.cmatrix == NULL) {
		fprintf(stderr, "Erreur allocation mémoire fractale\n");
		exit(1);
	}
	return f;
}

void Fractal_Destroy (fractal f) {
	free (f.fmatrix);
	free (f.zmatrix);
	free (f.cmatrix);
}

void Fractal_CalculateMatrix (fractal* f) {
	int i, j;
	complex zPixel;
	fractalresult result;

	printf ("Calculating Full EscapeTimeFractal Matrix  ...\n");

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

void Fractal_CalculateMatrix_DDp1 (fractal* f) {
	int i, j;
	double xg, yg;
	complex zPixel;
	fractalresult result;
	
	printf ("Calculating EscapeTimeFractal Matrix with DDp1 ...\n");
	
	for (j=0; j<f->ypixel; j=j+2) {
		for (i=0; i<f->xpixel; i=i+2) {
			
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

void Fractal_CalculateMatrix_DDp2 (fractal* f) {
	int i,starti, j, compteur;
	int up, down, left, right;
	complex zup, zdown, zleft, zright;
	double xg, yg, rz, iz;
	complex zPixel;
	fractalresult result;
	
	printf ("Calculating EscapeTimeFractal Matrix with DDp2 ...\n");
	
	for (j=1; j<f->ypixel; j=j+2) {
		for (i=1; i<f->xpixel; i=i+2) {
			
			xg = ((f->xmax-f->xmin)/f->xpixel)*i + f->xmin;
			yg = ((f->ymax-f->ymin)/f->ypixel)*j + f->ymin;
			zPixel = MakeComplex (xg, yg);
			
			result = FormulaSelector (*f,zPixel);
			*((f->fmatrix)+((f->xpixel*j)+i)) = result.iteration;
			*((f->zmatrix)+((f->xpixel*j)+i)) = result.z;
		}
	}

	// Dernier passage, cette fois on compare pour savoir si l'on peut se
	// permettre d'eviter le calcul, c'est la Divergence Detection (DD)
	
	for (j=0; j<f->ypixel; j++) {
		if (fmod (j,2) == 0) {starti=1;} else {starti=0;} // Lignes paires et impaires
		for (i=starti; i<f->xpixel; i=i+2) {
			
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

void Fractal_CalculateColorMatrix (fractal* f) {
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



Uint32 Fractal_Draw (SDL_Surface *canvas, fractal myfractal,int decalageX,int decalageY) {

	int i, j;
	Uint8 r, g, b;
	Uint32 time; // Test de temps de calcul en ms

	time = SDL_GetTicks();

	Fractal_CalculateMatrix_DDp1 (&myfractal);
	Fractal_CalculateColorMatrix (&myfractal);

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

	Fractal_CalculateMatrix_DDp2 (&myfractal);
	Fractal_CalculateColorMatrix (&myfractal);

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
	time = SDL_GetTicks() - time;
	printf ("Time Elapsed : %d ms\n", time);
	return time;
}

void Fractal_ChangeType (fractal* f, int type) {
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
	f->ymin = 1.700052002;
	f->ymax = -0.9000520021;
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

// Selectionne la formule
fractalresult FormulaSelector (fractal f, complex zPixel) {
	fractalresult r;

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
