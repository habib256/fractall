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
	VonKoch.h

fractall_LDFLAGS =  -lm

fractall_LDADD =  \
	$(IlmBase_LIBS) \
	$(ILMBASE_LIBS)
