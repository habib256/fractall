/* Gestionnaire d'exception inspire du magazine login HS5 */
/* exception.h */
/* VERHILLE Arnaud GPL2 */


#include <setjmp.h>

#ifndef EXCEPTIONFLAG
#define EXCEPTIONFLAG

typedef struct except_ctx {
  jmp_buf ctx;
  struct except_ctx *super;
} except;

extern except  *_exceptions;

#endif  /* EXCEPTIONFLAG */

// FUNCTION PROTOTYPES
// *******************
extern int push_exc();
extern int pop_exc(int);


/* Macros simples */

#define TRY { \
  int exc; \
  push_exc(); \
  exc=setjmp(_exceptions->ctx); \
  switch(exc) { \
  case 0:

#define CATCH(X) \
  break; \
  case (X): \
  exc=0; \

#define ENDTRY \
  break; \
  default : {}\
  } \
  pop_exc(exc); \
  }

#define THROW(X) longjmp(_exceptions->ctx, X)
