dnl Process this file with autoconf to produce a configure script.
dnl Created by Anjuta application wizard.

AC_INIT(fractall, 0.1)

AC_CONFIG_HEADERS([config.h])

SDL_REQUIRED=1.2.0

AM_INIT_AUTOMAKE([1.11])

AM_SILENT_RULES([yes])

AC_PROG_CC
AC_PATH_XTRA





LT_INIT

dnl ***************************************************************************
dnl Check for SDL
dnl ***************************************************************************
AM_PATH_SDL($SDL_REQUIRED,
            :,
	    AC_MSG_ERROR([SDL version $SDL_REQUIRED not found]))
CXXFLAGS="$CXXFLAGS $SDL_CFLAGS"
CFLAGS="$CFLAGS $SDL_CFLAGS"
CPPFLAGS="$CPPFLAGS $SDL_CFLAGS"
LIBS="$LIBS $SDL_LIBS"

dnl ***************************************************************************
dnl Check for GMP (optional, for infinite zoom)
dnl ***************************************************************************
AC_ARG_WITH([gmp],
    [AS_HELP_STRING([--without-gmp], [Disable GMP support for infinite zoom])],
    [use_gmp=$withval],
    [use_gmp=yes])

if test "x$use_gmp" != "xno"; then
    AC_CHECK_HEADER([gmp.h],
        [AC_CHECK_LIB([gmp], [__gmpz_init],
            [GMP_LIBS="-lgmp"
             AC_DEFINE([HAVE_GMP], [1], [Define if GMP is available])
             have_gmp=yes],
            [have_gmp=no])],
        [have_gmp=no])
else
    have_gmp=no
fi

AM_CONDITIONAL([HAVE_GMP], [test "x$have_gmp" = "xyes"])

if test "x$have_gmp" = "xyes"; then
    LIBS="$LIBS $GMP_LIBS"
    AC_MSG_NOTICE([GMP support enabled for infinite zoom])
else
    AC_MSG_NOTICE([GMP support disabled - zoom will be limited by double precision])
fi

dnl ***************************************************************************
dnl Check for OpenMP (optional, for parallel computation)
dnl ***************************************************************************
AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--disable-openmp], [Disable OpenMP parallel support])],
    [enable_openmp=$enableval],
    [enable_openmp=yes])

if test "x$enable_openmp" != "xno"; then
    AC_OPENMP
    if test "x$OPENMP_CFLAGS" != "x"; then
        AC_DEFINE([HAVE_OPENMP], [1], [Define if OpenMP is available])
        have_openmp=yes
        AC_MSG_NOTICE([OpenMP support enabled])
    else
        have_openmp=no
        AC_MSG_NOTICE([OpenMP not found - parallel computation disabled])
    fi
else
    have_openmp=no
    AC_MSG_NOTICE([OpenMP support disabled by user])
fi

AM_CONDITIONAL([HAVE_OPENMP], [test "x$have_openmp" = "xyes"])

dnl ***************************************************************************
dnl Check for SIMD support (SSE4.1, AVX, AVX2)
dnl ***************************************************************************
AC_ARG_ENABLE([simd],
    [AS_HELP_STRING([--disable-simd], [Disable SIMD optimizations])],
    [enable_simd=$enableval],
    [enable_simd=yes])

have_sse4_1=no
have_avx=no
have_avx2=no

if test "x$enable_simd" != "xno"; then
    # Vérifier SSE4.1
    AC_MSG_CHECKING([for SSE4.1 support])
    AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([[
            #include <smmintrin.h>
        ]], [[
            __m128d a = _mm_setzero_pd();
            __m128d b = _mm_dp_pd(a, a, 0xFF);
        ]])],
        [have_sse4_1=yes
         AC_DEFINE([HAVE_SSE4_1], [1], [Define if SSE4.1 is available])
         AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])])
    
    # Vérifier AVX
    AC_MSG_CHECKING([for AVX support])
    AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([[
            #include <immintrin.h>
        ]], [[
            __m256d a = _mm256_setzero_pd();
            __m256d b = _mm256_add_pd(a, a);
        ]])],
        [have_avx=yes
         AC_DEFINE([HAVE_AVX], [1], [Define if AVX is available])
         AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])])
    
    # Vérifier AVX2
    if test "x$have_avx" = "xyes"; then
        AC_MSG_CHECKING([for AVX2 support])
        AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([[
                #include <immintrin.h>
            ]], [[
                __m256i a = _mm256_setzero_si256();
                __m256i b = _mm256_add_epi32(a, a);
            ]])],
            [have_avx2=yes
             AC_DEFINE([HAVE_AVX2], [1], [Define if AVX2 is available])
             AC_MSG_RESULT([yes])],
            [AC_MSG_RESULT([no])])
    fi
    
    if test "x$have_sse4_1" = "xyes" || test "x$have_avx" = "xyes"; then
        AC_MSG_NOTICE([SIMD support enabled])
    else
        AC_MSG_NOTICE([SIMD support disabled - no compatible instructions found])
    fi
else
    AC_MSG_NOTICE([SIMD support disabled by user])
fi

AM_CONDITIONAL([HAVE_SSE4_1], [test "x$have_sse4_1" = "xyes"])
AM_CONDITIONAL([HAVE_AVX], [test "x$have_avx" = "xyes"])
AM_CONDITIONAL([HAVE_AVX2], [test "x$have_avx2" = "xyes"])

dnl ***************************************************************************
dnl Check for libpng (optional, for PNG export with metadata)
dnl ***************************************************************************
AC_ARG_WITH([png],
    [AS_HELP_STRING([--without-png], [Disable PNG export support])],
    [use_png=$withval],
    [use_png=yes])

have_png=no
if test "x$use_png" != "xno"; then
    AC_CHECK_HEADER([png.h],
        [AC_CHECK_LIB([png], [png_create_write_struct],
            [PNG_LIBS="-lpng"
             AC_DEFINE([HAVE_PNG], [1], [Define if libpng is available])
             have_png=yes],
            [have_png=no])],
        [have_png=no])
fi

AM_CONDITIONAL([HAVE_PNG], [test "x$have_png" = "xyes"])

if test "x$have_png" = "xyes"; then
    LIBS="$LIBS $PNG_LIBS"
    AC_MSG_NOTICE([PNG export support enabled])
else
    AC_MSG_NOTICE([PNG export support disabled - BMP fallback will be used])
fi











AC_OUTPUT([
Makefile
src/Makefile

])
