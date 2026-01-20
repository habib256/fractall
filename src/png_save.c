/* png_save.c
 * Sauvegarde PNG avec métadonnées pour fractall
 * Released under GPL2
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "png_save.h"
#include "complexmath.h"
#include "EscapeTime.h"
#ifdef HAVE_GMP
#include "precision_detector.h"
#endif

#ifdef HAVE_PNG
#include <png.h>

// Sauvegarder une surface SDL en PNG avec métadonnées de fractale
int SavePNGWithMetadata(SDL_Surface* surface, const char* filename, fractal* f, int typeFractale) {
	FILE* fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Cannot open file %s for writing\n", filename);
		return -1;
	}
	
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr) {
		fclose(fp);
		return -1;
	}
	
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		png_destroy_write_struct(&png_ptr, NULL);
		fclose(fp);
		return -1;
	}
	
	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
		fclose(fp);
		return -1;
	}
	
	png_init_io(png_ptr, fp);
	
	// Configuration de l'image
	int width = surface->w;
	int height = surface->h;
	int color_type = PNG_COLOR_TYPE_RGB;
	int bit_depth = 8;
	
	png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, color_type,
	             PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	
	// Ajouter les métadonnées texte
	png_text text_ptr[20];
	int num_text = 0;
	char value_buffer[256];
	
	// Type de fractale
	snprintf(value_buffer, sizeof(value_buffer), "%d", typeFractale);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "FractalType";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	// Coordonnées
	snprintf(value_buffer, sizeof(value_buffer), "%.15e", f->xmin);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "XMin";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	snprintf(value_buffer, sizeof(value_buffer), "%.15e", f->xmax);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "XMax";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	snprintf(value_buffer, sizeof(value_buffer), "%.15e", f->ymin);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "YMin";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	snprintf(value_buffer, sizeof(value_buffer), "%.15e", f->ymax);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "YMax";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	// Paramètres de calcul
	snprintf(value_buffer, sizeof(value_buffer), "%d", f->iterationMax);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "Iterations";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	snprintf(value_buffer, sizeof(value_buffer), "%d", f->bailout);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "Bailout";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	// Paramètres de couleur
	snprintf(value_buffer, sizeof(value_buffer), "%d", f->colorMode);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "ColorMode";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	snprintf(value_buffer, sizeof(value_buffer), "%d", f->colorRepeat);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "ColorRepeat";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	// Seed Julia (si applicable)
	if (typeFractale == 4 || typeFractale == 5) {
		snprintf(value_buffer, sizeof(value_buffer), "%.15e", Rez(f->seed));
		text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
		text_ptr[num_text].key = "JuliaSeedX";
		text_ptr[num_text].text = strdup(value_buffer);
		text_ptr[num_text].text_length = strlen(value_buffer);
		num_text++;
		
		snprintf(value_buffer, sizeof(value_buffer), "%.15e", Imz(f->seed));
		text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
		text_ptr[num_text].key = "JuliaSeedY";
		text_ptr[num_text].text = strdup(value_buffer);
		text_ptr[num_text].text_length = strlen(value_buffer);
		num_text++;
	}
	
	// Dimensions de l'image
	snprintf(value_buffer, sizeof(value_buffer), "%d", f->xpixel);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "ImageWidth";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	snprintf(value_buffer, sizeof(value_buffer), "%d", f->ypixel);
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "ImageHeight";
	text_ptr[num_text].text = strdup(value_buffer);
	text_ptr[num_text].text_length = strlen(value_buffer);
	num_text++;
	
	// Version du programme
	text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[num_text].key = "Software";
	text_ptr[num_text].text = strdup("fractall 0.5");
	text_ptr[num_text].text_length = strlen(text_ptr[num_text].text);
	num_text++;
	
	png_set_text(png_ptr, info_ptr, text_ptr, num_text);
	
	// Libérer les chaînes dupliquées après png_set_text
	for (int i = 0; i < num_text; i++) {
		free(text_ptr[i].text);
	}
	
	// Écrire les informations de l'image
	png_write_info(png_ptr, info_ptr);
	
	// Convertir et écrire les données de l'image
	png_bytep* row_pointers = (png_bytep*)malloc(height * sizeof(png_bytep));
	if (!row_pointers) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
		fclose(fp);
		return -1;
	}
	
	// Allouer un buffer pour les données RGB
	int row_bytes = width * 3;
	png_bytep image_data = (png_bytep)malloc(height * row_bytes);
	if (!image_data) {
		free(row_pointers);
		png_destroy_write_struct(&png_ptr, &info_ptr);
		fclose(fp);
		return -1;
	}
	
	// Convertir depuis SDL vers RGB
	SDL_LockSurface(surface);
	for (int y = 0; y < height; y++) {
		row_pointers[y] = image_data + y * row_bytes;
		for (int x = 0; x < width; x++) {
			Uint32 pixel;
			// Accéder au pixel selon le format de la surface
			if (surface->format->BytesPerPixel == 4) {
				pixel = ((Uint32*)surface->pixels)[y * (surface->pitch / 4) + x];
			} else if (surface->format->BytesPerPixel == 2) {
				pixel = ((Uint16*)surface->pixels)[y * (surface->pitch / 2) + x];
			} else {
				pixel = ((Uint8*)surface->pixels)[y * surface->pitch + x];
			}
			Uint8 r, g, b;
			SDL_GetRGB(pixel, surface->format, &r, &g, &b);
			image_data[y * row_bytes + x * 3 + 0] = r;
			image_data[y * row_bytes + x * 3 + 1] = g;
			image_data[y * row_bytes + x * 3 + 2] = b;
		}
	}
	SDL_UnlockSurface(surface);
	
	png_write_image(png_ptr, row_pointers);
	png_write_end(png_ptr, NULL);
	
	// Nettoyage
	free(image_data);
	free(row_pointers);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	fclose(fp);
	
	return 0;
}

// Charger les métadonnées depuis un PNG
int LoadPNGMetadata(const char* filename, fractal* f, int* typeFractale) {
	FILE* fp = fopen(filename, "rb");
	if (!fp) {
		return -1;
	}
	
	// Vérifier la signature PNG
	png_byte header[8];
	fread(header, 1, 8, fp);
	if (png_sig_cmp(header, 0, 8)) {
		fclose(fp);
		return -1;
	}
	
	png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr) {
		fclose(fp);
		return -1;
	}
	
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		png_destroy_read_struct(&png_ptr, NULL, NULL);
		fclose(fp);
		return -1;
	}
	
	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		fclose(fp);
		return -1;
	}
	
	png_init_io(png_ptr, fp);
	png_set_sig_bytes(png_ptr, 8);
	png_read_info(png_ptr, info_ptr);
	
	// Lire les métadonnées texte
	png_textp text_ptr = NULL;
	int num_text = 0;
	png_get_text(png_ptr, info_ptr, &text_ptr, &num_text);
	
	if (num_text == 0 || !text_ptr) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		fclose(fp);
		return -1;
	}
	
	// Parser les métadonnées
	for (int i = 0; i < num_text; i++) {
		if (strcmp(text_ptr[i].key, "FractalType") == 0) {
			*typeFractale = atoi(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "XMin") == 0) {
			f->xmin = atof(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "XMax") == 0) {
			f->xmax = atof(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "YMin") == 0) {
			f->ymin = atof(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "YMax") == 0) {
			f->ymax = atof(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "Iterations") == 0) {
			f->iterationMax = atoi(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "Bailout") == 0) {
			f->bailout = atoi(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "ColorMode") == 0) {
			f->colorMode = atoi(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "ColorRepeat") == 0) {
			f->colorRepeat = atoi(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "JuliaSeedX") == 0) {
			double seedX = atof(text_ptr[i].text);
			double seedY = Imz(f->seed);
			f->seed = MakeComplex(seedX, seedY);
		} else if (strcmp(text_ptr[i].key, "JuliaSeedY") == 0) {
			double seedX = Rez(f->seed);
			double seedY = atof(text_ptr[i].text);
			f->seed = MakeComplex(seedX, seedY);
		} else if (strcmp(text_ptr[i].key, "ImageWidth") == 0) {
			f->xpixel = atoi(text_ptr[i].text);
		} else if (strcmp(text_ptr[i].key, "ImageHeight") == 0) {
			f->ypixel = atoi(text_ptr[i].text);
		}
	}
	
	// Mettre à jour les coordonnées GMP si nécessaire
#ifdef HAVE_GMP
	if (f->use_gmp) {
		mpf_set_d(f->xmin_gmp, f->xmin);
		mpf_set_d(f->xmax_gmp, f->xmax);
		mpf_set_d(f->ymin_gmp, f->ymin);
		mpf_set_d(f->ymax_gmp, f->ymax);
		precision_update_fractal(f);
		precision_update_gmp_structures(f);
	}
#endif
	
	png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
	fclose(fp);
	
	return 0;
}

#else /* HAVE_PNG not defined */

// Stub si PNG n'est pas disponible
int SavePNGWithMetadata(SDL_Surface* surface, const char* filename, fractal* f, int typeFractale) {
	fprintf(stderr, "PNG support not compiled in. Please install libpng and reconfigure.\n");
	return -1;
}

int LoadPNGMetadata(const char* filename, fractal* f, int* typeFractale) {
	fprintf(stderr, "PNG support not compiled in. Please install libpng and reconfigure.\n");
	return -1;
}

#endif /* HAVE_PNG */
