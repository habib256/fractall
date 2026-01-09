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











AC_OUTPUT([
Makefile
src/Makefile

])
