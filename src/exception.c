/* exception.c */

#include <stdio.h>
#include <stdlib.h>
#include "exception.h"

except *_exceptions = NULL;

int push_exc() {
  except *layer;
  layer = (except *) malloc (sizeof(except));
  layer->super = _exceptions;
  _exceptions=layer;
  return 0;
}

int pop_exc(int ex) {
  except *super;
  super = _exceptions->super;
  free(_exceptions);
  _exceptions = super;
  if (ex) {
	  if (_exceptions) {
      longjmp (_exceptions->ctx, ex);
	  return 0; 
	  } else {
      fprintf(stderr, "Exception non gérée : %d\n", ex);
      exit (-1);
	  return 0;
    }
	  return 0;
  }
	return 0;
}
