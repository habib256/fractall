## Process this file with automake to produce Makefile.in

## Created by Anjuta

AM_CPPFLAGS = \
	-DPACKAGE_LOCALE_DIR=\""$(localedir)"\" \
	-DPACKAGE_SRC_DIR=\""$(srcdir)"\" \
	-DPACKAGE_DATA_DIR=\""$(pkgdatadir)"\"

AM_CFLAGS =\
	 -g

bin_PROGRAMS = fractall

fractall_SOURCES = \
	main.c \
	complexmath.c \
	complexmath.h \
	color_types.h \
	colorization.c \
	colorization.h \
	EscapeTime.c \
	EscapeTime.h \
	exception.c \
	exception.h \
	SDL_gfxPrimitives.c \
	SDL_gfxPrimitives.h \
	SDL_gfxPrimitives_font.h \
	SDLGUI.c \
	SDLGUI.h \
	SDL_imageFilter.c \
	SDL_imageFilter.h \
	SDL_rotozoom.c \
	SDL_rotozoom.h \
	VonKoch.c \
	VonKoch.h \
	png_save.c \
	png_save.h

if HAVE_SSE4_1
fractall_SOURCES += complexmath_simd.c
AM_CFLAGS += -msse4.1
endif

if HAVE_AVX
AM_CFLAGS += -mavx
endif

if HAVE_AVX2
AM_CFLAGS += -mavx2
endif

if HAVE_GMP
fractall_SOURCES += \
	complexmath_gmp.c \
	complexmath_gmp.h \
	precision_detector.c \
	precision_detector.h \
	gmp_pool.c \
	gmp_pool.h
endif

fractall_LDFLAGS =  -lm

if HAVE_OPENMP
AM_CFLAGS += $(OPENMP_CFLAGS)
fractall_LDFLAGS += $(OPENMP_CFLAGS)
endif

if HAVE_PNG
fractall_LDFLAGS += $(PNG_LIBS)
endif

fractall_LDADD =  \
	$(IlmBase_LIBS) \
	$(ILMBASE_LIBS)
